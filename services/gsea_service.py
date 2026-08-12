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
import threading
from contextlib import contextmanager
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple

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
# one or two and carries no information. Testing NES == 1.0 catches only the
# first regime; the size of the same-signed tail catches both, and it is
# available from the null distribution gseapy already computes.
#
# The *nominal p-value goes with it*. An earlier version of this module said
# the p-value survived a short tail. It does not: gseapy divides by the same
# same-signed tail it normalises against, so the p is quantised to multiples of
# 1/tail. Measured directly against gseapy 1.1.4 on a ranking of 200 genes, a
# term with a same-signed tail of 56 and one exceedance returned exactly
# 1/56 = 0.017857…, not 1/100 = 0.01, which is the pooled figure. A tail of two
# therefore admits only p ∈ {0, 0.5, 1}; a tail of one only p ∈ {0, 1}. Under
# NES_UNSTABLE both the magnitude and the p-value are coarse, and the coarseness
# is reported as ``p_value_resolution`` beside them.
#
# Incidence rises with gene-set size: a larger set tightens the permutation
# null, and a tight null sitting off zero in a skewed ranking is exactly the
# one-signed case. So this gets *more* likely as KE curation improves, which is
# the opposite of what anyone would expect a curation improvement to do.

#: NES reported as a number only when at least this many permutations fell on
#: the same side of zero as the observed ES. Below 10 the mean being divided by
#: has a relative standard error above ~30%, so the magnitude is not a quantity
#: to rank or threshold on — and the nominal p-value, computed on that same
#: tail, is quantised to steps of 1/10 or coarser.
MIN_SAME_SIGNED_NULL = 10

#: Both tails healthy — NES is an estimate and may be read as one.
NES_OK = 'ok'
#: No permutation fell on the observed side of zero. The observed ES beat every
#: single permutation, so the empirical p is below 1/permutation_num; the NES
#: itself is not normalisable and is reported as NaN rather than as a number.
NES_BEYOND_RESOLUTION = 'beyond_permutation_resolution'
#: Fewer than MIN_SAME_SIGNED_NULL same-signed permutations. Neither the NES
#: magnitude nor the nominal p-value is usable: both are computed against that
#: handful of values, so the p can only take the values k/tail. The achievable
#: resolution travels with the row in ``p_value_resolution``.
NES_UNSTABLE = 'unstable_normalisation'
#: The permutation null could not be read at all — the gseapy entry point this
#: module wraps was not called (renamed upstream, or the ``multilevel`` method
#: routed through ``fgsea_rs``). Distinct from NES_OK on purpose: "we looked and
#: the tails were healthy" and "we could not look" are different claims, and
#: collapsing the second into the first would let a gseapy upgrade reintroduce
#: the whole of #117 under a clean bill of health.
NES_UNDIAGNOSED = 'undiagnosed'

# --- Pooled-null FDR recomputation (issue #122) ----------------------------
#
# #117 repaired the degenerate row itself. It could not repair anybody else's
# q-value, because ``gsea_fdr`` had already run inside the Rust extension before
# this module saw a result — and that function pools across terms.
#
# GSEA's FDR is a ratio of two tail fractions (gseapy/algorithm.py:632-681).
# For an observed NES* >= 0:
#
#     FDR = [ #{NESnull >= NES*} / #{NESnull >= 0} ]
#           ------------------------------------------
#           [ #{NESobs  >= NES*} / #{NESobs  >= 0} ]
#
# Both pools are shared by every term. The substituted NES = 1.0 of a degenerate
# term does *not* enter the permutation pool in the numerator — a term whose
# null lies entirely on one side of zero contributes no same-signed null
# entries, which is correct. It enters the **observed** vector in the
# denominator, and that vector normalises everybody else's q.
#
# The bias is therefore one-sided and derivable. Let h = #{genuine NESobs >=
# NES*} and p = #{genuine NESobs >= 0}. Injecting one fabricated NES of 1.0:
#
#   * NES* > 1.0   — h unchanged, p becomes p+1, so the denominator fraction
#                    falls and the reported FDR is *inflated* (conservative).
#   * 0 <= NES* <= 1.0 — h becomes h+1 and p becomes p+1; since h <= p the
#                    fraction rises and the reported FDR is *deflated*
#                    (anti-conservative).
#
# So the strongest results are penalised and the weak-to-moderate ones
# flattered, systematically rather than randomly. Negatively-enriched terms are
# untouched: a fabricated +1.0 never enters the negative pool. On a grid of
# seven terms with one substitution the multiplier on the denominator is
# (p+1)/p ~ 1.15-1.3, i.e. a 15-30% shift in q for terms above NES 1.0 — not
# negligible against a 0.05 cutoff.
#
# The repair is to recompute the pooled FDR here, from the same permutation
# nulls #117 already captures, treating the degenerate terms honestly. See
# :func:`recompute_pooled_fdr` for why they are counted as +/-infinity rather
# than dropped.

#: A degenerate term's observed ES beat every permutation, so its normalised
#: score is unbounded on that side. Counting it at every same-signed threshold
#: is the faithful limit of the statistic and invents no magnitude.
DEGENERATE_INFINITY = 'infinity'
#: Drop degenerate terms from both pools. This is what issue #122 proposed; it
#: is available for the comparison the issue asks for, but it is not the default
#: — see :func:`recompute_pooled_fdr`.
DEGENERATE_EXCLUDE = 'exclude'

#: The q-value in this row was computed by this module from the captured
#: permutation nulls, with degenerate terms handled per :data:`DEGENERATE_INFINITY`.
FDR_RECOMPUTED = 'recomputed'
#: The q-value in this row is gseapy's own, untouched. Either the permutation
#: nulls could not be captured, or the recomputation declined this row.
FDR_GSEAPY = 'gseapy'

#: Serialises the process-global patch in :func:`_capture_prerank_summaries`.
#:
#: The capture rebinds ``gseapy.gsea.prerank_rs``, which is process-global, and
#: this application is concurrent: Flask runs threaded and every batch condition
#: is dispatched on its own ``threading.Thread``, so two ``run_gsea_analysis``
#: calls overlapping is routine rather than exotic. Without this lock the two
#: patches nest, every sink receives both runs' summaries, and the per-term
#: counts — keyed by KE ID, which is identical across conditions — are
#: last-write-wins. That does not merely lose the diagnostics: it writes one
#: condition's tail count onto another condition's row, and the repair in
#: :func:`apply_null_diagnostics` then rewrites that row's NES, p and FDR on the
#: strength of it. The lock covers the patch *and* the ``gp.prerank`` call it
#: exists for, so GSEA runs serialise against each other. gseapy is already
#: pinned to ``threads=1`` for determinism, so this costs wall-clock on
#: concurrent batches and nothing else.
_PRERANK_PATCH_LOCK = threading.Lock()


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
def _capture_prerank_summaries(sink: List[Any]) -> Iterator[Dict[str, bool]]:
    """Collect the per-term GSEA summaries ``gseapy.prerank`` throws away.

    The permutation null distribution is the only place the empty-tail
    condition of #117 is visible, and it is computed inside the Rust extension.
    ``Prerank.run()`` binds its result to a local, hands the summaries to
    ``to_df`` and returns None, so nothing survives on the object gseapy hands
    back: ``res2d`` carries ES/NES/p/FDR but not ``esnull``.

    Rather than reimplement prerank, this wraps the extension entry point
    ``gseapy.gsea.prerank_rs`` and keeps the summaries. Verified against gseapy
    1.1.4 and 1.3.0, which both resolve that name as a module global at call
    time. If a future release renames it, or the caller selects the
    ``multilevel`` method (which routes through ``fgsea_rs`` instead), the sink
    stays empty, ``captured`` stays False and every row is labelled
    :data:`NES_UNDIAGNOSED` — the analysis proceeds on gseapy's own numbers, but
    it does not claim to have checked them.

    The rebind is process-global, so the whole block runs under
    :data:`_PRERANK_PATCH_LOCK`; see that constant for what overlapping runs did
    without it. The restore is also identity-checked, because an unlocked nested
    patch used to leave the global bound to a dead wrapper appending to an
    abandoned list for the life of the process.

    Args:
        sink: List the captured ``GSEASummary`` objects are appended to.

    Yields:
        dict: ``{'captured': bool}`` — flipped to True once the wrapper has
        actually run, i.e. once the diagnostics are known to be available for
        this call. Read it after the block, not inside it. The patch is reverted
        on exit, including on exception.
    """
    state: Dict[str, bool] = {'captured': False}
    try:
        gsea_module = importlib.import_module('gseapy.gsea')
    except ImportError:  # pragma: no cover — gseapy is a hard dependency
        gsea_module = None

    # The lock is taken *before* the global is read. Reading it first and
    # locking afterwards is the same defect in slower motion: a run that arrives
    # while another holds the patch would capture the other run's wrapper as its
    # "original" and restore that on exit, leaving the global permanently bound
    # to a dead wrapper appending to an abandoned list.
    with _PRERANK_PATCH_LOCK:
        original = getattr(gsea_module, 'prerank_rs', None)
        if original is None:
            logger.error(
                "gseapy.gsea.prerank_rs not found — the GSEA permutation null "
                "cannot be inspected, so the degenerate one-signed null of "
                "issue #117 cannot be detected in this run. Every Key Event "
                "will be reported with nes_status='%s'. This usually means "
                "gseapy renamed the entry point; the wrapper needs updating.",
                NES_UNDIAGNOSED,
            )
            yield state
            return

        def _wrapped(*args, **kwargs):
            result = original(*args, **kwargs)
            sink.extend(list(getattr(result, 'summaries', None) or []))
            state['captured'] = True
            return result

        setattr(gsea_module, 'prerank_rs', _wrapped)
        try:
            yield state
        finally:
            if getattr(gsea_module, 'prerank_rs', None) is not _wrapped:
                logger.error(
                    "gseapy.gsea.prerank_rs was rebound by something other "
                    "than this capture while it was patched — restoring the "
                    "real function anyway. GSEA null diagnostics for this run "
                    "may be unreliable."
                )
            setattr(gsea_module, 'prerank_rs', original)


def _null_summaries(summaries: List[Any]) -> Dict[str, Tuple[float, np.ndarray]]:
    """Extract each term's observed ES and its permutation null (issues #117, #122).

    The captured ``GSEASummary`` objects are the only place the permutation
    distribution survives; :func:`_capture_prerank_summaries` explains why they
    have to be intercepted at all. #117 needed one integer per term and threw
    the arrays away, which is the sole reason the pooled-null recomputation of
    #122 could not be written against them.

    Args:
        summaries: ``GSEASummary`` objects captured from ``prerank_rs``.

    Returns:
        Term → ``(observed ES, finite permutation ES array)``. Terms whose
        summary cannot be read, or whose null contains no finite value, are
        omitted — downstream that reads as "not assessed", which is a different
        claim from a tail of zero.
    """
    out: Dict[str, Tuple[float, np.ndarray]] = {}
    for summary in summaries:
        try:
            observed = float(summary.es)
            null = np.asarray(summary.esnull, dtype=float)
        except (AttributeError, TypeError, ValueError):
            continue
        null = null[np.isfinite(null)]
        if null.size == 0:
            continue
        out[str(summary.term)] = (observed, null)
    return out


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
    for term, (observed, null) in _null_summaries(summaries).items():
        same_signed = (null >= 0) if observed >= 0 else (null < 0)
        sizes[term] = int(same_signed.sum())
    return sizes


def recompute_pooled_fdr(
    summaries: List[Any],
    degenerate: str = DEGENERATE_INFINITY,
) -> Dict[str, float]:
    """Recompute the pooled-null FDR without the fabricated NES (issue #122).

    A transcription of gseapy's own ``normalize`` (algorithm.py:567-609) and
    ``gsea_fdr`` (algorithm.py:632-681), not an independent derivation. That is
    deliberate: the point of the change is to remove one contaminating value
    from the pools, so every other step must stay bit-comparable to the
    implementation being replaced. The equivalence test in
    ``tests/test_gsea_service.py`` asserts exactly that — on a run with no
    degenerate term this function must reproduce gseapy's q-values to ~1e-12 —
    and it is the check that makes the change reviewable.

    **Why degenerate terms are counted as +/-infinity rather than dropped.**
    Issue #122 proposed excluding them. Excluding removes a term that genuinely
    belongs in the observed distribution, so it deflates the denominator for
    everybody else — the same class of error as the one being fixed, smaller and
    with the opposite sign. The observed ES beat every permutation, so the
    normalised score is unbounded on that side; counting it at every same-signed
    threshold is the faithful limit and requires no invented magnitude. Under
    that treatment the degenerate term's own q comes out as 0.0, which is the
    limit #117 already writes into the row by hand — the two agree rather than
    having to be reconciled. :data:`DEGENERATE_EXCLUDE` implements the other
    choice so the two can be compared on real data, which is what #122's second
    acceptance criterion asks for.

    Note that +/-infinity is used *inside the pooling only*. The reported NES
    for such a term stays NaN: the score is not normalisable, and printing an
    infinity in the results table would be a magnitude claim the data does not
    support.

    Args:
        summaries: ``GSEASummary`` objects captured from ``prerank_rs``.
        degenerate: :data:`DEGENERATE_INFINITY` or :data:`DEGENERATE_EXCLUDE`.

    Returns:
        Term → recomputed q-value, capped at 1.0 exactly as gseapy caps it.
        Empty when the nulls could not be read or are not rectangular, in which
        case the caller leaves gseapy's own numbers in place. Terms whose
        denominator would be zero are omitted for the same reason: gseapy
        substitutes 1e9 there, which is not a q-value and must not reach a
        results table.
    """
    if degenerate not in (DEGENERATE_INFINITY, DEGENERATE_EXCLUDE):
        raise ValueError(
            f"degenerate must be {DEGENERATE_INFINITY!r} or "
            f"{DEGENERATE_EXCLUDE!r}, got {degenerate!r}"
        )

    per_term = _null_summaries(summaries)
    if not per_term:
        return {}

    # Every term is permuted the same number of times, so a ragged set of rows
    # means the capture picked up two overlapping runs (the condition
    # _PRERANK_PATCH_LOCK exists to prevent) or a partially-read summary. Either
    # way the pools would be built from mismatched material, so decline rather
    # than pool them.
    widths = {null.size for _, null in per_term.values()}
    if len(widths) != 1:
        logger.error(
            "GSEA: permutation nulls are not rectangular (%s distinct widths: "
            "%s) — the pooled FDR of issue #122 cannot be recomputed and "
            "gseapy's own q-values stand. This usually means two runs' "
            "summaries were captured together.",
            len(widths), sorted(widths),
        )
        return {}

    terms = list(per_term)
    es = np.asarray([per_term[t][0] for t in terms], dtype=float)
    esnull = np.vstack([per_term[t][1] for t in terms])

    # gseapy/algorithm.py:596-607, transcribed. A fully-masked row — no
    # permutation on the observed side — leaves the mean masked; gseapy lets the
    # 0/0 through and the Rust path substitutes 1.0. Filling with NaN instead
    # keeps the degenerate rows identifiable so they can be replaced below,
    # rather than silently becoming a number.
    with np.errstate(invalid='ignore', divide='ignore'):
        pos_mean = np.ma.filled(
            np.ma.MaskedArray(esnull, mask=(esnull < 0)).mean(axis=1), np.nan
        ).astype(float)
        neg_mean = np.ma.filled(
            np.ma.MaskedArray(esnull, mask=(esnull >= 0)).mean(axis=1), np.nan
        ).astype(float)
        nes = np.where(es >= 0, es / pos_mean, -es / neg_mean)
        # The null rows stay well defined even for a degenerate term: its
        # permutations all sit on the other side of zero, so they normalise
        # against that side's mean, which exists. This is why the contamination
        # is confined to the observed vector.
        nesnull = np.where(
            esnull >= 0,
            esnull / pos_mean[:, np.newaxis],
            -esnull / neg_mean[:, np.newaxis],
        )

    degenerate_mask = ~np.isfinite(nes)
    if degenerate == DEGENERATE_INFINITY:
        nes = np.where(degenerate_mask, np.where(es >= 0, np.inf, -np.inf), nes)
        keep = np.ones(len(terms), dtype=bool)
    else:
        keep = ~degenerate_mask

    # A non-finite value in the *null* pool is not a regime this handles; drop
    # those entries rather than let np.sort park them at the end where
    # searchsorted would count them as the largest values in the distribution.
    nvals = np.sort(nesnull[np.isfinite(nesnull)].ravel())
    nnes = np.sort(nes[keep])
    if nvals.size == 0 or nnes.size == 0:
        return {}

    fdrs: Dict[str, float] = {}
    for i, term in enumerate(terms):
        if not keep[i]:
            continue
        value = nes[i]
        # gseapy/algorithm.py:649-677, transcribed. searchsorted over the two
        # sorted pools is gseapy's own formulation of "count at or beyond this
        # threshold, on this side of zero".
        if value >= 0:
            all_pos = nvals.size - int(np.searchsorted(nvals, 0, side='left'))
            all_higher = nvals.size - int(np.searchsorted(nvals, value, side='left'))
            nes_pos = nnes.size - int(np.searchsorted(nnes, 0, side='left'))
            nes_higher = nnes.size - int(np.searchsorted(nnes, value, side='left'))
        else:
            all_pos = int(np.searchsorted(nvals, 0, side='left'))
            all_higher = int(np.searchsorted(nvals, value, side='right'))
            nes_pos = int(np.searchsorted(nnes, 0, side='left'))
            nes_higher = int(np.searchsorted(nnes, value, side='right'))

        if all_pos == 0 or nes_pos == 0:
            # gseapy returns 1e9 here. That is a sentinel, not a q-value, and
            # this module's output goes straight into a results table and a
            # report, so the row is left with gseapy's number and labelled as
            # gseapy's rather than given a fabricated one.
            continue
        fdr = (all_higher / float(all_pos)) / (nes_higher / float(nes_pos))
        fdrs[term] = float(min(fdr, 1.0))

    return fdrs


def apply_null_diagnostics(
    res: pd.DataFrame,
    same_signed_sizes: Dict[str, int],
    permutation_num: int,
    diagnostics_captured: bool = True,
    recomputed_fdr: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    """Classify each result row by its same-signed null tail and repair it (#117).

    Adds ``null_same_signed_n``, ``nes_status`` and ``p_value_resolution``, and
    for the degenerate regime replaces the values gseapy invents. It is a
    separate function because the regimes it distinguishes are the whole content
    of #117 and they must be testable without paying for a thousand
    permutations.

    An empty same-signed tail is not a missing result: the observed ES beat
    every permutation, so the empirical p is bounded above by
    ``1/permutation_num`` and the FDR limit is 0 — nothing in the pooled null
    is more extreme than a score no permutation reached. The NES genuinely has
    no value, being a ratio against an empty sample, so it is set to NaN and
    the status column carries the reason. NaN is the one thing that cannot be
    read as a weak result; the direction survives in the ES column beside it.

    ``p_value_resolution`` is the finest p-value the row could have produced.
    gseapy divides the exceedance count by the *same-signed* tail, not by
    ``permutation_num``, so a row with a tail of two can only report
    p ∈ {0, 0.5, 1} — a p of 0.000 there means "zero of two", and a p of 0.600
    on a tail of five means "three of five". Reporting the step size beside the
    p-value is the only way a reader can tell those from a p estimated on a
    thousand permutations. Under :data:`NES_BEYOND_RESOLUTION` the column
    carries ``1/permutation_num``, the bound the substituted p-value is.

    A KE missing from ``same_signed_sizes`` is labelled
    :data:`NES_UNDIAGNOSED`, never ``ok``: the diagnostics were not available,
    which is a different statement from a healthy tail, and gseapy's numbers are
    left exactly as they came.

    Args:
        res: Result frame with at least 'KE', 'NES', 'p_value', 'FDR'.
        same_signed_sizes: KE → same-signed permutation count.
        permutation_num: Permutations run, i.e. the coarsest bound on the
            p-value when the same-signed tail is empty.
        diagnostics_captured: Whether the null-capture wrapper actually ran.
            False means the wrapper never fired — the failure is in this
            module's grip on gseapy, not in the data — and is logged as an
            error rather than a warning.
        recomputed_fdr: KE → q-value from :func:`recompute_pooled_fdr`
            (issue #122). Where a KE is present its q replaces gseapy's, which
            was pooled against a fabricated NES; where it is absent gseapy's
            number stands. Either way the row records which one it carries in
            ``fdr_source``, because "recomputed" and "gseapy's own" are
            different statistical claims and a reader thresholding on q is
            entitled to know which they are looking at.

    Returns:
        pd.DataFrame: ``res``, mutated in place and returned for chaining.
    """
    # Issue #122 — replace the pooled q before the regime classification below
    # reads it, so a degenerate row's hand-written limit is applied on top of
    # the recomputed value rather than being overwritten by it.
    recomputed = pd.Series(
        res['KE'].map(recomputed_fdr or {}), index=res.index, dtype=float
    )
    res['fdr_source'] = np.where(recomputed.notna(), FDR_RECOMPUTED, FDR_GSEAPY)
    res['FDR'] = recomputed.where(recomputed.notna(), res['FDR'])

    res['null_same_signed_n'] = res['KE'].map(same_signed_sizes)
    tail = res['null_same_signed_n']
    undiagnosed = tail.isna()
    beyond = tail == 0
    unstable = (tail > 0) & (tail < MIN_SAME_SIGNED_NULL)
    res['nes_status'] = NES_OK
    res.loc[unstable, 'nes_status'] = NES_UNSTABLE
    res.loc[beyond, 'nes_status'] = NES_BEYOND_RESOLUTION
    res.loc[undiagnosed, 'nes_status'] = NES_UNDIAGNOSED
    # The step size the nominal p-value was quantised to (see docstring).
    res['p_value_resolution'] = 1.0 / tail.where(tail > 0)
    res.loc[beyond, 'p_value_resolution'] = 1.0 / permutation_num
    if beyond.any():
        res.loc[beyond, 'NES'] = np.nan
        res.loc[beyond, 'p_value'] = 1.0 / permutation_num
        # The FDR limit is 0: nothing in the pooled null is more extreme than a
        # score no permutation reached. Written by hand only where the row was
        # not recomputed — under DEGENERATE_INFINITY the recomputation derives
        # the same 0.0 from the pools, and a disagreement between the two would
        # mean the infinity treatment had stopped being the limit it claims to
        # be, which is worth an error rather than a silent overwrite.
        untouched = beyond & (res['fdr_source'] == FDR_GSEAPY)
        res.loc[untouched, 'FDR'] = 0.0
        disagreed = res.loc[
            beyond & (res['fdr_source'] == FDR_RECOMPUTED) & (res['FDR'] != 0.0),
            'KE',
        ].tolist()
        if disagreed:
            logger.error(
                "GSEA: the recomputed pooled FDR (issue #122) is non-zero for "
                "%d Key Event(s) whose observed ES beat every permutation, "
                "where the limit is 0: %s. The infinity treatment in "
                "recompute_pooled_fdr is no longer producing the limit it "
                "documents.",
                len(disagreed), sorted(disagreed),
            )
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
            "permutations (issue #117): %s. Neither the NES magnitude nor the "
            "nominal p-value is usable there — gseapy computes the p on that "
            "same short tail, so it is quantised to steps of 1/tail (reported "
            "per row as p_value_resolution).",
            int(unstable.sum()), MIN_SAME_SIGNED_NULL,
            sorted(
                f"{ke} (tail {int(n)}, p in steps of {1.0 / int(n):.3g})"
                for ke, n in zip(
                    res.loc[unstable, 'KE'], res.loc[unstable, 'null_same_signed_n']
                )
            ),
        )
    if undiagnosed.any():
        logger.error(
            "GSEA: the permutation null could not be read for %d of %d Key "
            "Event(s)%s, so the degenerate one-signed null of issue #117 could "
            "not be ruled out for them. They are reported with "
            "nes_status='%s' — that is not a clean bill of health, and their "
            "NES and p-value are gseapy's unchecked output: %s",
            int(undiagnosed.sum()), len(res),
            "" if diagnostics_captured
            else " (the capture wrapper never ran — gseapy has probably "
                 "renamed its entry point, or method='multilevel' routed "
                 "around it)",
            NES_UNDIAGNOSED,
            sorted(res.loc[undiagnosed, 'KE'].tolist()),
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
        ['KE', 'Title', 'NES', 'nes_status', 'ES', 'p_value',
         'p_value_resolution', 'FDR', 'null_same_signed_n', 'lead_genes',
         'total_KE_genes_in_dataset']

        ``nes_status`` is one of ``NES_OK``, ``NES_UNSTABLE``,
        ``NES_BEYOND_RESOLUTION`` or ``NES_UNDIAGNOSED`` (issue #117). Under
        ``NES_BEYOND_RESOLUTION`` the NES is NaN — the quantity does not
        exist — while ``p_value`` carries the bound ``1/permutation_num``,
        ``FDR`` is 0.0 and ``ES`` still carries the direction. Under
        ``NES_UNSTABLE`` every number is gseapy's own, but neither the NES
        magnitude nor the p-value may be ranked or thresholded on: gseapy
        computes both against the same short same-signed tail, so the p-value
        is quantised to the step size in ``p_value_resolution``. Under
        ``NES_UNDIAGNOSED`` the null could not be read at all and nothing about
        the row has been checked. ``null_same_signed_n`` is the evidence for
        the call and is NaN when the diagnostics could not be read.

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
    with _capture_prerank_summaries(summaries) as capture_state:
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
    diagnostics_captured = bool(capture_state.get('captured'))
    same_signed_sizes = _null_tail_sizes(summaries)
    # Issue #122 — the pooled FDR gseapy returned was normalised against an
    # observed-NES vector containing a fabricated 1.0 for every degenerate term.
    # Recompute it from the same captured nulls; an empty result means the
    # nulls were unreadable and gseapy's own q-values stand, labelled as such.
    recomputed_fdr = recompute_pooled_fdr(summaries)

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
    res = apply_null_diagnostics(
        res, same_signed_sizes, permutation_num,
        diagnostics_captured=diagnostics_captured,
        recomputed_fdr=recomputed_fdr,
    )
    undiagnosed_kes = int((res['nes_status'] == NES_UNDIAGNOSED).sum())
    if recomputed_fdr:
        logger.info(
            "GSEA: recomputed the pooled-null FDR for %d of %d Key Event(s) "
            "(issue #122); the rest carry gseapy's own q-value.",
            int((res['fdr_source'] == FDR_RECOMPUTED).sum()), len(res),
        )
    else:
        logger.warning(
            "GSEA: the pooled-null FDR could not be recomputed (issue #122), "
            "so every q-value in this run is gseapy's own and was normalised "
            "against a pool that may contain a substituted NES."
        )

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
        # Issue #117 / R5: a run that could not inspect its own permutation null
        # has to say so where the reader is, not only in the server log.
        nes_undiagnosed_kes=undiagnosed_kes,
    )

    # Sort by FDR ascending; use KE ID as tie-breaker so row order is fully
    # deterministic regardless of the internal ordering gseapy returns.
    res = res.sort_values(['FDR', 'KE'], ascending=[True, True], kind='mergesort')

    out = res[['KE', 'Title', 'NES', 'nes_status', 'ES', 'p_value',
               'p_value_resolution', 'FDR', 'fdr_source', 'null_same_signed_n',
               'lead_genes', 'total_KE_genes_in_dataset']]
    # Attach on the returned object: pandas only propagates .attrs through
    # operations that call __finalize__, so set it last.
    out.attrs[KE_SUMMARY_ATTR] = ke_summary
    return out
