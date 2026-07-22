"""
Gene set enrichment analysis (GSEA) service.

Provides gseapy.prerank-based enrichment analysis as an alternative to the
Fisher's exact test path in enrichment_service.py. The dispatcher
(enrichment_service.run_enrichment) routes to this module based on the
user's method selection.
"""
import importlib
import math
import logging
from contextlib import contextmanager
from typing import Any, Dict, Iterator, List, Optional, Set

import numpy as np
import pandas as pd
import gseapy as gp

from config import Config
from services.enrichment_service import (
    EXCLUDED_NO_MAPPING,
    EXCLUDED_TOO_FEW_GENES,
    EXCLUDED_TOO_MANY_GENES,
    EXCLUDED_UNRESOLVED_MAPPING,
    KE_SUMMARY_ATTR,
    _build_ke_summary,
    normalise_unresolved_ke_pathways,
)

logger = logging.getLogger(__name__)

# --- NES normalisation regimes (issue #117) --------------------------------
#
# GSEA normalises an observed enrichment score by the mean of the *same-signed*
# tail of that gene set's permutation null. When the null lands entirely on the
# other side of zero that tail is empty, the normalisation is 0/0, and gseapy
# falls back to NES = 1.0 with p = 1.0 rather than declaring the cell
# undefined. Those are the values a thoroughly unenriched gene set produces, so
# a maximally enriched Key Event is reported as null — the one direction of
# error that cannot be caught by reading the table.
#
# The failure is continuous, not binary. A tail of one or two permutations
# still yields a finite NES, but its magnitude is a ratio against a sample of
# one or two and carries no information, even though the nominal p-value is
# sound. Testing NES == 1.0 catches only the first regime; the size of the
# same-signed tail catches both, and it is available from the null distribution
# gseapy already computes.
#
# Incidence rises with gene-set size: a larger set tightens the permutation
# null, and a tight null sitting off zero in a skewed ranking is exactly the
# one-signed case. So this gets *more* likely as KE curation improves, which is
# the opposite of what anyone would expect a curation improvement to do.

#: NES reported as a number only when at least this many permutations fell on
#: the same side of zero as the observed ES. Below 10 the mean being divided by
#: has a relative standard error above ~30%, so the magnitude is not a quantity
#: to rank or threshold on.
MIN_SAME_SIGNED_NULL = 10

#: Both tails healthy — NES is an estimate and may be read as one.
NES_OK = 'ok'
#: No permutation fell on the observed side of zero. The observed ES beat every
#: single permutation, so the empirical p is below 1/permutation_num; the NES
#: itself is not normalisable and is reported as NaN rather than as a number.
NES_BEYOND_RESOLUTION = 'beyond_permutation_resolution'
#: Fewer than MIN_SAME_SIGNED_NULL same-signed permutations. The p-value stands;
#: the NES magnitude is normalised against a handful of values and does not.
NES_UNSTABLE = 'unstable_normalisation'


def _build_ranking(df: pd.DataFrame) -> pd.Series:
    """Build a gene ranking Series for gseapy.prerank.

    Ranking metric: sign(log2FC) × −log10(pval).

    Edge cases:
    - Zero p-values are clamped to the smallest non-zero p-value in the
      dataset so no +inf values appear in the metric.
    - Ties are broken by ascending gene symbol (mergesort stable sort).

    Args:
        df: DataFrame with columns 'ID', 'log2FC', 'pval'.

    Returns:
        pd.Series indexed by upper-cased gene symbols, sorted descending by
        the ranking metric (highest-ranked genes first).
    """
    pvals = df['pval'].astype(float)
    nonzero_min = pvals[pvals > 0].min()
    pvals_clipped = pvals.clip(lower=nonzero_min)

    metric = np.sign(df['log2FC'].astype(float)) * -np.log10(pvals_clipped)

    s = pd.Series(metric.values, index=df['ID'].astype(str).str.upper())

    # Deterministic tie-breaking: sort by gene symbol ascending first
    # (mergesort preserves this order for equal metric values after the
    # second stable sort)
    s = s.sort_index(kind='mergesort')
    s = s.sort_values(ascending=False, kind='mergesort')

    return s


@contextmanager
def _capture_prerank_summaries(sink: List[Any]) -> Iterator[None]:
    """Collect the per-term GSEA summaries ``gseapy.prerank`` throws away.

    The permutation null distribution is the only place the empty-tail
    condition of #117 is visible, and it is computed inside the Rust extension.
    ``Prerank.run()`` binds its result to a local, hands the summaries to
    ``to_df`` and returns None, so nothing survives on the object gseapy hands
    back: ``res2d`` carries ES/NES/p/FDR but not ``esnull``.

    Rather than reimplement prerank, this wraps the extension entry point
    ``gseapy.gsea.prerank_rs`` for the duration of one call and keeps the
    summaries. Verified against gseapy 1.1.4 and 1.3.0, which both resolve that
    name as a module global at call time. If a future release renames it, or
    the caller selects the ``multilevel`` method (which routes through
    ``fgsea_rs`` instead), the sink stays empty and the analysis proceeds
    without the diagnostics rather than failing — the numbers are still
    gseapy's own.

    Args:
        sink: List the captured ``GSEASummary`` objects are appended to.

    Yields:
        None. The patch is reverted on exit, including on exception.
    """
    try:
        gsea_module = importlib.import_module('gseapy.gsea')
    except ImportError:  # pragma: no cover — gseapy is a hard dependency
        gsea_module = None
    original = getattr(gsea_module, 'prerank_rs', None)
    if original is None:
        logger.warning(
            "gseapy.gsea.prerank_rs not found — GSEA null diagnostics "
            "(issue #117) are unavailable for this run"
        )
        yield
        return

    def _wrapped(*args, **kwargs):
        result = original(*args, **kwargs)
        sink.extend(list(getattr(result, 'summaries', None) or []))
        return result

    setattr(gsea_module, 'prerank_rs', _wrapped)
    try:
        yield
    finally:
        setattr(gsea_module, 'prerank_rs', original)


def _null_tail_sizes(summaries: List[Any]) -> Dict[str, int]:
    """Count each term's same-signed permutation tail (issue #117).

    "Same-signed" is the side of zero the *observed* ES fell on, which is the
    side GSEA normalises against and the side the nominal p-value counts on.
    A count of zero means no permutation reached the observed side at all.

    Args:
        summaries: ``GSEASummary`` objects captured from ``prerank_rs``.

    Returns:
        Term → number of same-signed permutations. Terms whose summary cannot
        be read are omitted, which downstream reads as "not assessed" rather
        than as a count of zero.
    """
    sizes: Dict[str, int] = {}
    for summary in summaries:
        try:
            observed = float(summary.es)
            null = np.asarray(summary.esnull, dtype=float)
        except (AttributeError, TypeError, ValueError):
            continue
        null = null[np.isfinite(null)]
        if null.size == 0:
            continue
        same_signed = (null >= 0) if observed >= 0 else (null < 0)
        sizes[str(summary.term)] = int(same_signed.sum())
    return sizes


def apply_null_diagnostics(
    res: pd.DataFrame,
    same_signed_sizes: Dict[str, int],
    permutation_num: int,
) -> pd.DataFrame:
    """Classify each result row by its same-signed null tail and repair it (#117).

    Adds ``null_same_signed_n`` and ``nes_status``, and for the degenerate
    regime replaces the values gseapy invents. It is a separate function
    because the two regimes it distinguishes are the whole content of #117 and
    they must be testable without paying for a thousand permutations.

    An empty same-signed tail is not a missing result: the observed ES beat
    every permutation, so the empirical p is bounded above by
    ``1/permutation_num`` and the FDR limit is 0 — nothing in the pooled null
    is more extreme than a score no permutation reached. The NES genuinely has
    no value, being a ratio against an empty sample, so it is set to NaN and
    the status column carries the reason. NaN is the one thing that cannot be
    read as a weak result; the direction survives in the ES column beside it.

    A KE missing from ``same_signed_sizes`` keeps ``NaN`` there and status
    ``ok``: the diagnostics were not available, which is a different statement
    from a tail of zero, and it must not silently rewrite gseapy's numbers.

    Args:
        res: Result frame with at least 'KE', 'NES', 'p_value', 'FDR'.
        same_signed_sizes: KE → same-signed permutation count.
        permutation_num: Permutations run, i.e. the resolution of the p-value.

    Returns:
        pd.DataFrame: ``res``, mutated in place and returned for chaining.
    """
    res['null_same_signed_n'] = res['KE'].map(same_signed_sizes)
    tail = res['null_same_signed_n']
    beyond = tail == 0
    unstable = (tail > 0) & (tail < MIN_SAME_SIGNED_NULL)
    res['nes_status'] = NES_OK
    res.loc[unstable, 'nes_status'] = NES_UNSTABLE
    res.loc[beyond, 'nes_status'] = NES_BEYOND_RESOLUTION
    if beyond.any():
        res.loc[beyond, 'NES'] = np.nan
        res.loc[beyond, 'p_value'] = 1.0 / permutation_num
        res.loc[beyond, 'FDR'] = 0.0
        logger.warning(
            "GSEA: %d Key Event(s) had no same-signed permutation and were "
            "reported by gseapy as NES=1.0, p=1.0 — corrected to "
            "p < 1/%d with NES undefined (issue #117): %s. Incidence of this "
            "rises with gene-set size, so expect more of it as curation grows "
            "the Key Event gene sets.",
            int(beyond.sum()), permutation_num,
            sorted(res.loc[beyond, 'KE'].tolist()),
        )
    if unstable.any():
        logger.warning(
            "GSEA: %d Key Event(s) normalised against fewer than %d same-signed "
            "permutations — the p-value stands but the NES magnitude is not "
            "usable (issue #117): %s",
            int(unstable.sum()), MIN_SAME_SIGNED_NULL,
            sorted(res.loc[unstable, 'KE'].tolist()),
        )
    if not same_signed_sizes:
        logger.warning(
            "GSEA: permutation null diagnostics unavailable — degenerate "
            "one-signed nulls (issue #117) cannot be detected in this run"
        )
    return res


def run_gsea_analysis(
    df: pd.DataFrame,
    reference_sets: Dict[str, Set[str]],
    ke_list: Set[str],
    ke_title_map: Dict[str, str],
    *,
    min_size: int = 5,
    max_size: int = Config.GSEA_MAX_KE_GENES,
    permutation_num: int = 1000,
    seed: int = 42,
    unresolved_ke_pathways: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """Run GSEA (prerank) enrichment analysis for Key Events.

    Uses gseapy.prerank with locked parameters (D-01, SPEC § Constraints):
    - seed=42 and threads=1 guarantee deterministic results (D-15.1).
    - permutation_num=1000 is not user-adjustable this phase.
    - NES values are server-side clamped to ±3 before the result leaves
      Python (D-12; the frontend gradient also clamps as a second layer).

    Args:
        df: Processed gene expression DataFrame with columns 'ID', 'log2FC',
            'pval'. The 'significant' column (if present) is not used by GSEA.
        reference_sets: Mapping of KE_ID → set of gene symbols.
        ke_list: Set of KE IDs to include in this analysis.
        ke_title_map: Mapping of KE IDs to human-readable titles.
        min_size: Minimum number of overlapping genes for a KE to be included
            (keyword-only so locked constants cannot be overridden positionally).
        max_size: Maximum number of measured genes a KE may carry and still be
            tested (keyword-only). gseapy drops anything larger before it
            computes a statistic. The Broad convention is 500; 1000 is used
            here because three of the seven testable AOP 472 Key Events already
            exceed 500 measured genes, and dropping them would cost more than
            the dilution costs. A drop on this bound is accounted separately
            from a ``min_size`` drop (issue #120) — the two say opposite things
            about the Key Event.
        permutation_num: Number of permutations for the null distribution
            (keyword-only).
        seed: Random seed for reproducibility (keyword-only).
        unresolved_ke_pathways: Optional KE_ID -> pathway IDs that are mapped
            to that KE but could not be resolved to genes (issue #81). Same
            meaning and effect as in ``run_enrichment_analysis`` — ORA and
            GSEA report exclusions with one vocabulary.

    Returns:
        pd.DataFrame sorted by FDR ascending with columns:
        ['KE', 'Title', 'NES', 'nes_status', 'ES', 'p_value', 'FDR',
         'null_same_signed_n', 'lead_genes', 'total_KE_genes_in_dataset']

        ``nes_status`` is one of ``NES_OK``, ``NES_UNSTABLE`` or
        ``NES_BEYOND_RESOLUTION`` (issue #117). Under
        ``NES_BEYOND_RESOLUTION`` the NES is NaN — the quantity does not
        exist — while ``p_value`` carries the bound ``1/permutation_num``,
        ``FDR`` is 0.0 and ``ES`` still carries the direction. Under
        ``NES_UNSTABLE`` every number is gseapy's own; only the NES *magnitude*
        should not be ranked or thresholded on. ``null_same_signed_n`` is the
        evidence for the call and is NaN when the diagnostics could not be
        read.

    Raises:
        ValueError: When no reference gene sets are found for the selected AOP,
            or when gseapy returns an empty result set.
    """
    logger.info("Starting GSEA enrichment analysis")

    # Build ranking vector from gene expression data
    ranking = _build_ranking(df)

    # Filter and normalise reference gene sets — mirrors enrichment_service.py:43-52
    gene_sets = {
        ke: list(frozenset(g.strip().upper() for g in reference_sets[ke]))
        for ke in ke_list
        if ke in reference_sets
    }

    if not gene_sets:
        logger.error("No reference gene sets found for the selected AOP")
        raise ValueError("No reference gene sets found for the selected AOP")

    # Pre-compute overlap counts for total_KE_genes_in_dataset column
    # (gseapy res2d does not expose this directly)
    kes_to_overlap_count = {
        ke: len(set(genes) & set(ranking.index))
        for ke, genes in gene_sets.items()
    }

    # Run gseapy.prerank
    # threads=1 is the determinism guard per RESEARCH §Pitfall #2 — load-bearing for D-15.1
    summaries: List[Any] = []
    with _capture_prerank_summaries(summaries):
        pre_res = gp.prerank(
            rnk=ranking,
            gene_sets=gene_sets,
            outdir=None,
            min_size=min_size,
            max_size=max_size,
            permutation_num=permutation_num,
            seed=seed,
            threads=1,
            no_plot=True,
            verbose=False,
        )
    same_signed_sizes = _null_tail_sizes(summaries)

    res = pre_res.res2d.copy()

    # Rename gseapy wire-format columns (Title-cased) to our internal schema
    res = res.rename(columns={
        'Term': 'KE',
        'NOM p-val': 'p_value',
        'FDR q-val': 'FDR',
        'Lead_genes': 'lead_genes',
    })

    if res.empty:
        logger.error("No enrichment results generated by GSEA")
        raise ValueError("No enrichment results generated")

    # Normalise lead_genes to a Python list. gseapy emits a ';'-joined string;
    # downstream consumers (results table truncation, report export) expect a
    # list so the shape is consistent regardless of backend wire format.
    res['lead_genes'] = res['lead_genes'].apply(
        lambda s: [g for g in s.split(';') if g] if isinstance(s, str)
        else (list(s) if s is not None else [])
    )

    # Map KE IDs to human-readable titles (fall back to the KE ID itself)
    res['Title'] = res['KE'].map(ke_title_map).fillna(res['KE'])

    # Server-side NES clamp to ±3 per D-12 (frontend gradient also clamps
    # as a second layer so the wire payload stays interpretable)
    for col in ('ES', 'NES', 'p_value', 'FDR'):
        res[col] = pd.to_numeric(res[col], errors='coerce')
    res['NES'] = res['NES'].clip(lower=-3, upper=3)

    # Issue #117 — classify each cell by the size of its same-signed
    # permutation tail and repair the two regimes gseapy misreports.
    res = apply_null_diagnostics(res, same_signed_sizes, permutation_num)

    # Populate total_KE_genes_in_dataset from pre-computed overlap counts
    res['total_KE_genes_in_dataset'] = (
        res['KE'].map(kes_to_overlap_count).fillna(0).astype(int)
    )

    # Log KEs dropped by the size filters for parity with the Fisher per-KE
    # skip log (RESEARCH Open Question #2) — diff the input ke_list against
    # res['KE']. Issue #120: gseapy drops on *both* bounds and says which is
    # not recoverable from its output, so the two are separated here using the
    # overlap counts computed above. Reporting an over-large Key Event as
    # "fewer than min_size measured genes" is not a cosmetic slip: it sends the
    # reader to curate more genes onto the very set that is already too big.
    dropped = (set(ke_list) & set(reference_sets.keys())) - set(res['KE'])
    dropped_too_many = {
        ke for ke in dropped
        if gene_sets.get(ke) and kes_to_overlap_count.get(ke, 0) > max_size
    }
    dropped_too_few = dropped - dropped_too_many
    if dropped_too_few:
        logger.info(
            f"GSEA dropped {len(dropped_too_few)} KEs due to min_size={min_size} "
            f"filter: {sorted(dropped_too_few)}"
        )
    if dropped_too_many:
        logger.warning(
            "GSEA dropped %d KEs due to max_size=%d filter: %s. These are the "
            "best-covered Key Events in the analysis, not under-measured ones "
            "(issue #120)",
            len(dropped_too_many), max_size,
            sorted(
                f"{ke} ({kes_to_overlap_count.get(ke, 0)} measured genes)"
                for ke in dropped_too_many
            ),
        )

    # Issue #65 — same tested/excluded accounting the Fisher backend produces,
    # so the results page, network and reports read one shape regardless of
    # method.
    # Issue #81 — three distinct exclusions, not two: a KE with no curated
    # mapping, a KE whose mapping resolved to no genes, and a KE whose gene set
    # is known but under-measured in this dataset (zero measured genes included
    # — that is a coverage fact about the upload, not a curation gap).
    unresolved_map = normalise_unresolved_ke_pathways(unresolved_ke_pathways)
    unresolved_named = {}
    excluded_reasons = {}
    for ke in set(ke_list) - set(reference_sets.keys()):
        if unresolved_map.get(ke):
            excluded_reasons[ke] = EXCLUDED_UNRESOLVED_MAPPING
            unresolved_named[ke] = unresolved_map[ke]
        else:
            excluded_reasons[ke] = EXCLUDED_NO_MAPPING
    # Issue #120 — a fourth exclusion: a KE whose measured gene set exceeds the
    # GSEA size ceiling. It is dropped before any statistic is computed, and it
    # is the opposite fact about the KE from the min_size case.
    for ke in dropped:
        if not gene_sets.get(ke):
            excluded_reasons[ke] = EXCLUDED_UNRESOLVED_MAPPING
            if unresolved_map.get(ke):
                unresolved_named[ke] = unresolved_map[ke]
        elif ke in dropped_too_many:
            excluded_reasons[ke] = EXCLUDED_TOO_MANY_GENES
        else:
            excluded_reasons[ke] = EXCLUDED_TOO_FEW_GENES
    ke_summary = _build_ke_summary(
        ke_list, len(res), excluded_reasons, min_size,
        unresolved_pathways_by_ke=unresolved_named,
        max_ke_genes=max_size,
    )

    # Sort by FDR ascending; use KE ID as tie-breaker so row order is fully
    # deterministic regardless of the internal ordering gseapy returns.
    res = res.sort_values(['FDR', 'KE'], ascending=[True, True], kind='mergesort')

    out = res[['KE', 'Title', 'NES', 'nes_status', 'ES', 'p_value', 'FDR',
               'null_same_signed_n', 'lead_genes',
               'total_KE_genes_in_dataset']]
    # Attach on the returned object: pandas only propagates .attrs through
    # operations that call __finalize__, so set it last.
    out.attrs[KE_SUMMARY_ATTR] = ke_summary
    return out
