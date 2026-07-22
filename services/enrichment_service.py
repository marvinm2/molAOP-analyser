"""
Gene set enrichment analysis service.
"""
import math
import pandas as pd
import logging
from typing import Dict, Set, List, Any, Tuple, Optional
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from config import Config

logger = logging.getLogger(__name__)

# Issue #65 — key under which the tested/excluded KE accounting is attached to
# the result DataFrame's ``.attrs``. Using ``.attrs`` keeps the return type a
# plain DataFrame, so every existing call site (app.py:/analyze,
# services/batch_service.py, the GSEA backend) is unaffected; callers that want
# the accounting read it back with ``get_ke_summary()``.
KE_SUMMARY_ATTR = 'ke_summary'

# Issue #81 — the reasons a Key Event never reached a statistical test. These
# strings are the shared vocabulary of the ORA backend, the GSEA backend, the
# network's per-node ``excluded_reason`` and the accounting sentence, so all of
# them make the same claim about the same KE.
#
# ``no_mapping`` and ``unresolved_mapping`` used to be one counter, and that
# conflation pointed the reader at the wrong person: "no gene set mapped" says
# the curation is incomplete and sends them to the Builder to add a mapping,
# while a mapping that exists but resolves to no genes is a stale-reference-data
# bug in this tool. A KE whose gene set is known but barely measured in the
# uploaded dataset is a third thing again, and belongs with ``too_few_genes``
# (zero measured genes is "fewer than the minimum", not "not curated").
EXCLUDED_NO_MAPPING = 'no_mapping'
EXCLUDED_UNRESOLVED_MAPPING = 'unresolved_mapping'
EXCLUDED_TOO_FEW_GENES = 'too_few_genes'
EXCLUDED_ERROR = 'error'


def normalise_unresolved_ke_pathways(
    unresolved_ke_pathways: Optional[Dict[str, Any]],
) -> Dict[str, List[str]]:
    """Normalise a KE -> unresolvable-pathway-ID mapping (issue #81).

    Args:
        unresolved_ke_pathways: KE_ID -> iterable of pathway IDs that are
            curated and mapped to that KE but whose gene membership no
            configured source could resolve (see ``load_reference_sets``'s
            ``unresolved_out`` and the resolution entry's
            ``unresolved_pathways``). A plain string is accepted for a single
            pathway. None or empty yields an empty dict.

    Returns:
        KE_ID -> sorted list of unique, upper-cased pathway IDs, with entries
        that name no pathway dropped.
    """
    normalised: Dict[str, List[str]] = {}
    for ke, pathways in (unresolved_ke_pathways or {}).items():
        if not pathways:
            continue
        if isinstance(pathways, str):
            pathways = [pathways]
        ids = sorted({str(p).strip().upper() for p in pathways if str(p).strip()})
        if ids:
            normalised[str(ke)] = ids
    return normalised


def get_ke_summary(results: pd.DataFrame) -> Optional[Dict[str, Any]]:
    """Return the tested/excluded KE accounting attached to an enrichment result.

    Args:
        results: DataFrame returned by ``run_enrichment_analysis`` /
            ``run_gsea_analysis`` / ``run_enrichment``.

    Returns:
        The summary dict (see ``run_enrichment_analysis``) or None when the
        DataFrame carries no accounting (e.g. results restored from JSON, where
        pandas ``.attrs`` do not survive the round-trip).
    """
    try:
        return results.attrs.get(KE_SUMMARY_ATTR)
    except AttributeError:
        return None


def format_ke_summary(summary: Optional[Dict[str, Any]]) -> str:
    """Render a KE accounting summary as a one-line human-readable sentence.

    Example: ``"18 of 24 Key Events tested; 4 excluded (fewer than 5 measured
    genes), 2 excluded (no gene set mapped), 1 excluded (mapped, but no genes
    could be resolved: WP5477)"``. Shared by the results page, the single report
    and the batch report so all three word it identically.

    Issue #81: the mapped-but-unresolvable clause is deliberately separate from
    "no gene set mapped" and names the offending pathways. The two sentences
    ask different things of the reader — one is a curation gap to fill in the
    Builder, the other is stale reference data in this tool — and a summary
    that says the first when it means the second sends the reader to fix
    something that is not broken.

    Args:
        summary: Dict from ``get_ke_summary`` (or an equivalently-shaped dict
            rebuilt from a network payload). None/empty yields an empty string.
            Keys added after a summary was stored are read with defaults, so
            payloads written before #81 (which carry only the two older
            exclusion counters) still render.

    Returns:
        The formatted sentence, or '' when there is nothing to report.
    """
    if not summary:
        return ''
    parts = [
        f"{summary.get('tested', 0)} of {summary.get('total_kes', 0)} "
        f"Key Events tested"
    ]
    clauses = []
    too_few = summary.get('excluded_too_few_genes', 0)
    if too_few:
        clauses.append(
            f"{too_few} excluded (fewer than "
            f"{summary.get('min_ke_genes', Config.MIN_KE_GENES)} measured genes)"
        )
    no_mapping = summary.get('excluded_no_mapping', 0)
    if no_mapping:
        clauses.append(f"{no_mapping} excluded (no gene set mapped)")
    unresolved = summary.get('excluded_unresolved_mapping', 0)
    if unresolved:
        clauses.append(
            f"{unresolved} excluded (mapped, but no genes could be resolved"
            f"{_format_unresolved_pathways(summary.get('unresolved_pathways'))})"
        )
    errored = summary.get('excluded_error', 0)
    if errored:
        clauses.append(f"{errored} excluded (statistics could not be computed)")
    if clauses:
        parts.append('; ' + ', '.join(clauses))
    return ''.join(parts)


# How many unresolvable pathway IDs the accounting sentence names before it
# summarises the rest. Naming them is the point of the clause -- an unnamed
# count cannot be checked against the Builder -- but the sentence still has to
# stay a sentence.
MAX_NAMED_UNRESOLVED_PATHWAYS = 3


def _format_unresolved_pathways(pathways: Optional[List[str]]) -> str:
    """Render the ``": WP5477, WP5498 and 2 more"`` tail of the mapped-but-empty clause.

    Args:
        pathways: Pathway IDs that could not be resolved to genes, or None.

    Returns:
        The tail including its leading ``": "``, or '' when nothing is known.
    """
    if not pathways:
        return ''
    shown = list(pathways)[:MAX_NAMED_UNRESOLVED_PATHWAYS]
    tail = ', '.join(shown)
    remaining = len(pathways) - len(shown)
    if remaining > 0:
        tail += f" and {remaining} more"
    return f": {tail}"


# Below this fraction of the uploaded genes matching the reference universe,
# the chosen identifier column is almost certainly the wrong one (issue #69) --
# e.g. an Ensembl column tested against symbol-keyed gene sets. Real symbol
# datasets overlap far above this even when the AOP is a poor fit.
MIN_BACKGROUND_OVERLAP_FRACTION = 0.05


def assess_background_overlap(
    user_genes: Set[str],
    reference_sets: Dict[str, Set[str]],
) -> Dict[str, Any]:
    """Check that the uploaded identifiers actually match the gene sets.

    Issue #69: enrichment matches uppercase HGNC symbols. Point it at a column
    of Ensembl accessions and every test still runs -- the background looks
    plausible, the KE table is fully populated, and the overlap is zero. That
    produces a well-formed, entirely meaningless result, which is the failure
    mode most likely to reach a published figure unnoticed.

    Args:
        user_genes: Normalised (uppercase) identifiers from the upload.
        reference_sets: KE -> gene set mapping, before AOP filtering, so the
            check reflects the whole reference universe rather than one AOP.

    Returns:
        Dict with ``matched``, ``total``, ``fraction``, ``universe_size`` and
        ``is_suspect``. ``is_suspect`` is True when the overlap is below
        :data:`MIN_BACKGROUND_OVERLAP_FRACTION`, meaning the caller should warn
        rather than present the result as normal.
    """
    universe = set()
    for genes in reference_sets.values():
        universe.update(g.strip().upper() for g in genes)

    normalised = {str(g).strip().upper() for g in user_genes}
    matched = normalised & universe
    total = len(normalised)
    fraction = (len(matched) / total) if total else 0.0

    result = {
        'matched': len(matched),
        'total': total,
        'fraction': fraction,
        'universe_size': len(universe),
        'is_suspect': bool(total and universe and fraction < MIN_BACKGROUND_OVERLAP_FRACTION),
    }

    if result['is_suspect']:
        logger.warning(
            "Only %d of %d uploaded identifiers (%.2f%%) match the %d-gene reference "
            "universe -- the selected gene ID column is probably not gene symbols",
            len(matched), total, fraction * 100, len(universe),
        )
    return result


# Issue #70 — the three states a tested Key Event can be in once both tails of
# the Fisher test are run. 'ns' means "tested, neither tail significant", which
# is what the one-sided-only output used to conflate with depletion.
REPRESENTATION_ENRICHED = 'enriched'
REPRESENTATION_DEPLETED = 'depleted'
REPRESENTATION_NS = 'ns'

# Human-readable forms, shared by the results page, the reports and the
# documentation so all three word it identically.
REPRESENTATION_LABELS = {
    REPRESENTATION_ENRICHED: 'Enriched',
    REPRESENTATION_DEPLETED: 'Depleted',
    REPRESENTATION_NS: 'Not significant',
}


def _call_representation(
    fdr_enriched: float,
    fdr_depleted: float,
    fdr_cutoff: float,
) -> str:
    """Classify one tested KE as enriched, depleted or non-significant (#70).

    Args:
        fdr_enriched: BH-adjusted p-value of the over-representation tail.
        fdr_depleted: BH-adjusted p-value of the under-representation tail.
        fdr_cutoff: Cutoff both are compared against.

    Returns:
        One of the ``REPRESENTATION_*`` constants. Over-representation wins a
        tie: the two tails cannot both be significant except in degenerate
        tables, and over-representation is the headline statistic.
    """
    if fdr_enriched is not None and fdr_enriched < fdr_cutoff:
        return REPRESENTATION_ENRICHED
    if fdr_depleted is not None and fdr_depleted < fdr_cutoff:
        return REPRESENTATION_DEPLETED
    return REPRESENTATION_NS


def run_enrichment_analysis(
    df: pd.DataFrame,
    reference_sets: Dict[str, Set[str]],
    ke_list: Set[str],
    ke_title_map: Dict[str, str],
    gene_logfc_map: Optional[Dict[str, float]] = None,
    min_ke_genes: int = Config.MIN_KE_GENES,
    fdr_cutoff: float = Config.SIGNIFICANCE_FDR_CUTOFF,
    unresolved_ke_pathways: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """
    Run Fisher's exact test enrichment analysis for Key Events.

    Args:
        df: Processed gene expression dataframe with 'significant' column
        reference_sets: Dictionary mapping KE_ID to gene sets
        ke_list: Set of KE IDs to analyze
        ke_title_map: Mapping of KE IDs to titles
        gene_logfc_map: Optional mapping of gene to log2FC. When supplied,
            the result DataFrame gains a ``Direction`` column rendering
            ``"N↑ / M↓"`` -- counts of overlap genes with positive vs negative
            log2FC. AOP KEs are intrinsically directional (e.g. titles like
            "Down Regulation, HMGCS2" or "Inhibition, FoxA2") so curators need
            the observed direction alongside the direction-agnostic Fisher p
            to spot KEs whose overlap goes the wrong way. When ``None``, the
            column is omitted (default behaviour unchanged for callers that
            don't opt in).

        min_ke_genes: Minimum number of a KE's genes that must be measured in
            ``df`` for the KE to be tested (issue #65). Defaults to
            ``Config.MIN_KE_GENES``. KEs below it are excluded from the Fisher
            tests and from the BH denominator, and are recorded in the KE
            accounting summary instead of being dropped silently.

        unresolved_ke_pathways: Optional KE_ID -> pathway IDs that are mapped
            to that KE but whose gene membership no source could resolve
            (issue #81; see ``load_reference_sets``'s ``unresolved_out``).
            Such a KE is absent from ``reference_sets`` exactly like an
            uncurated one, so without this hint the two are indistinguishable
            here and both would be reported as "no gene set mapped" — telling
            the reader the curation is incomplete when in fact the mapping
            exists and this tool's reference data is stale.

        fdr_cutoff: BH-adjusted cutoff used to fill the ``Representation``
            column (issue #70). Both tails are always computed; this only
            decides where the pre-rendered ``enriched`` / ``depleted`` / ``ns``
            call is drawn. Consumers that let the user move the cutoff (the
            results-page selector, the network) recompute the call from
            ``FDR`` and ``FDR_depleted`` instead of reading this column.

    Returns:
        pd.DataFrame: Enrichment results sorted by FDR. Alongside the
        over-representation ``p_value``/``FDR``, every tested KE also carries
        ``p_value_depleted``/``FDR_depleted`` (the under-representation tail)
        and a ``Representation`` call of ``enriched``/``depleted``/``ns``
        (issue #70). The tested/excluded KE
        accounting is attached to ``.attrs[KE_SUMMARY_ATTR]`` (issue #65) and
        is read back with ``get_ke_summary()``::

            {
                'total_kes': 24,              # KEs in the AOP
                'tested': 18,                 # rows in the returned DataFrame
                'excluded_no_mapping': 2,     # no curated gene set at all
                'excluded_unresolved_mapping': 1,  # mapped, zero genes resolved
                'excluded_too_few_genes': 3,  # overlap below min_ke_genes
                'excluded_error': 0,          # contingency / Fisher failure
                'min_ke_genes': 5,
                'unresolved_pathways': ['WP5477'],
                'unresolved_pathways_by_ke': {'KE:1115': ['WP5477']},
                'excluded_reasons': {'KE:123': 'too_few_genes', ...},
            }
    """
    logger.info("Starting enrichment analysis")
    
    # Filter reference sets to only include KEs in the AOP
    # Pre-normalize gene symbols once to avoid repeated .strip().upper() in the loop
    filtered_reference_sets = {
        ke: frozenset(g.strip().upper() for g in genes)
        for ke, genes in reference_sets.items() if ke in ke_list
    }

    if not filtered_reference_sets:
        logger.error("No reference gene sets found for the selected AOP")
        raise ValueError("No reference gene sets found for the selected AOP")

    # Define gene universe and significance status
    all_genes = set(df['ID'])
    user_gene_status = dict(zip(df['ID'], df['significant']))
    
    # Use all measured genes as statistical background
    # This provides consistent comparisons across different AOPs
    sig_genes = df['significant'].sum()
    non_sig_genes = df.shape[0] - sig_genes
    
    logger.info(f"Statistical background: {len(all_genes)} total genes ({sig_genes} significant, {non_sig_genes} non-significant)")
    
    # Columnar result accumulation — avoids per-row dict allocation
    r_ke = []; r_title = []; r_total = []; r_sig_in = []; r_sig_not = []
    r_nsig_in = []; r_nsig_not = []; r_pct = []; r_odds = []; r_pval = []
    r_overlap = []; r_noverlap = []
    r_pval_dep = []  # issue #70 — the other tail (under-representation)
    r_direction = []  # populated only when gene_logfc_map is supplied
    emit_direction = gene_logfc_map is not None

    # Issue #65 — record why each KE never reached a Fisher test, so the
    # multiple-testing denominator is visible and "could not assess" is
    # distinguishable from "assessed and not enriched" downstream.
    # Issue #81 — and so that "not curated" is distinguishable from "curated,
    # but this tool could not resolve the mapping to genes".
    unresolved_map = normalise_unresolved_ke_pathways(unresolved_ke_pathways)
    unresolved_named: Dict[str, List[str]] = {}
    excluded_reasons: Dict[str, str] = {}
    for ke in ke_list:
        if ke in filtered_reference_sets:
            continue
        if unresolved_map.get(ke):
            excluded_reasons[ke] = EXCLUDED_UNRESOLVED_MAPPING
            unresolved_named[ke] = unresolved_map[ke]
        else:
            excluded_reasons[ke] = EXCLUDED_NO_MAPPING

    for ke, ref_genes in filtered_reference_sets.items():
        try:
            # A KE present in the reference sets but carrying no genes is
            # mapped — the mapping simply resolved to nothing (issue #81).
            if not ref_genes:
                logger.debug(f"KE {ke} is mapped but its gene set is empty")
                excluded_reasons[ke] = EXCLUDED_UNRESOLVED_MAPPING
                if unresolved_map.get(ke):
                    unresolved_named[ke] = unresolved_map[ke]
                continue

            # Find overlap between KE genes and user genes (ref_genes already normalized)
            ke_genes = ref_genes & all_genes

            # Skip KEs with too few genes for reliable statistics. Issue #81:
            # zero measured genes belongs here too — the KE has a gene set,
            # none of it was measured in this dataset, and calling that "no
            # gene set mapped" blamed the curator for the experiment's
            # coverage.
            if len(ke_genes) < max(min_ke_genes, 1):
                logger.debug(
                    f"Skipping KE {ke}: only {len(ke_genes)} of its "
                    f"{len(ref_genes)} genes measured "
                    f"(minimum {min_ke_genes} required)"
                )
                excluded_reasons[ke] = EXCLUDED_TOO_FEW_GENES
                continue

            sig_in_ke = {g for g in ke_genes if user_gene_status.get(g, False)}
            non_sig_in_ke = ke_genes - sig_in_ke

            # Create 2x2 contingency table
            a = len(sig_in_ke)                    # significant & in KE
            b = sig_genes - a                     # significant & not in KE
            c = len(non_sig_in_ke)                # non-significant & in KE
            d = non_sig_genes - c                 # non-significant & not in KE

            # Validate contingency table cells are non-negative
            if b < 0 or d < 0:
                logger.error(f"KE {ke}: invalid contingency table (a={a}, b={b}, c={c}, d={d}) — skipping")
                excluded_reasons[ke] = EXCLUDED_ERROR
                continue

            # Run Fisher's exact test (one-tailed, greater)
            odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

            # Issue #70 — the opposite tail. Over-representation stays the
            # headline statistic, but a KE whose gene set is conspicuously
            # spared while the rest of the transcriptome responds is a real,
            # interpretable result; tested on one tail only it is reported as
            # p ≈ 1 and reads exactly like "nothing here".
            _, pval_dep = fisher_exact([[a, b], [c, d]], alternative="less")

            # Add statistical warnings for edge cases
            if a == 0:
                logger.warning(f"KE {ke}: No significant genes found - result may not be meaningful")
            elif len(ke_genes) < 10:
                logger.debug(f"KE {ke}: Small gene set ({len(ke_genes)} genes) - interpret with caution")
            elif a / len(ke_genes) > 0.9:
                logger.debug(f"KE {ke}: Very high enrichment ratio ({a}/{len(ke_genes)}) - verify biological relevance")

            logger.debug(f"KE {ke}: {len(ke_genes)} total genes, {a} significant overlap (p={pval:.2e})")

            r_ke.append(ke)
            r_title.append(ke_title_map.get(ke, ke))
            r_total.append(len(ke_genes))
            r_sig_in.append(a)
            r_sig_not.append(b)
            r_nsig_in.append(c)
            r_nsig_not.append(d)
            r_pct.append(round((a / len(ke_genes)) * 100, 1) if ke_genes else 0)
            r_odds.append(round(odds, 4) if (not pd.isna(odds) and not math.isinf(odds)) else 'NA')
            r_pval.append(pval)
            r_pval_dep.append(pval_dep)
            r_overlap.append(", ".join(sorted(sig_in_ke)))
            r_noverlap.append(a)

            if emit_direction:
                # Observed-direction counts over the significant overlap genes.
                # log2FC == 0 (rare; only with threshold=0) counts as neither.
                n_up = 0
                n_down = 0
                for g in sig_in_ke:
                    fc = gene_logfc_map.get(g)
                    if fc is None:
                        continue
                    if fc > 0:
                        n_up += 1
                    elif fc < 0:
                        n_down += 1
                r_direction.append(f"{n_up}↑ / {n_down}↓")

        except Exception as e:
            logger.error(f"Error processing KE {ke}: {e}")
            excluded_reasons[ke] = EXCLUDED_ERROR
            continue

    if not r_ke:
        logger.error("No enrichment results generated")
        raise ValueError("No enrichment results generated")

    # Build DataFrame from columnar data and apply FDR correction.
    # Column order matters: results.html and shared_results.html render the
    # table by iterating dict keys, so `Direction` lands right after
    # `num_overlap` -- read alongside the overlap count it summarises.
    fdr_enriched = multipletests(r_pval, method="fdr_bh")[1]
    # Issue #70 — each tail is corrected within its own family of tests, so the
    # enrichment FDR keeps exactly the value it had before depletion testing
    # was added and the two columns stay independently interpretable.
    fdr_depleted = multipletests(r_pval_dep, method="fdr_bh")[1]

    columns = {
        'Title': r_title,
        'p_value': r_pval,
        'FDR': fdr_enriched,
        'Representation': [
            _call_representation(fdr_up, fdr_down, fdr_cutoff)
            for fdr_up, fdr_down in zip(fdr_enriched, fdr_depleted)
        ],
        'num_overlap': r_noverlap,
    }
    if emit_direction:
        columns['Direction'] = r_direction
    columns.update({
        'pct_sig_in_KE': r_pct,
        'total_KE_genes_in_dataset': r_total,
        'odds_ratio': r_odds,
        'overlap': r_overlap,
        'KE': r_ke,
        'sig_in_KE': r_sig_in,
        'sig_not_KE': r_sig_not,
        'non_sig_in_KE': r_nsig_in,
        'non_sig_not_KE': r_nsig_not,
        'p_value_depleted': r_pval_dep,
        'FDR_depleted': fdr_depleted,
    })
    df_results = pd.DataFrame(columns)
    df_results = df_results.sort_values("FDR")

    # Issue #65 — attach the KE accounting after the sort, since pandas only
    # propagates .attrs through operations that call __finalize__.
    df_results.attrs[KE_SUMMARY_ATTR] = _build_ke_summary(
        ke_list, len(r_ke), excluded_reasons, min_ke_genes,
        unresolved_pathways_by_ke=unresolved_named,
    )

    logger.info(
        f"Enrichment analysis completed: {len(df_results)} results generated "
        f"({format_ke_summary(df_results.attrs[KE_SUMMARY_ATTR])})"
    )

    return df_results


def _build_ke_summary(
    ke_list: Set[str],
    tested: int,
    excluded_reasons: Dict[str, str],
    min_ke_genes: int,
    unresolved_pathways_by_ke: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, Any]:
    """Assemble the tested/excluded KE accounting dict (issue #65).

    Args:
        ke_list: All KE IDs of the selected AOP (the accounting denominator).
        tested: Number of KEs that reached a statistical test — i.e. the number
            of rows in the result table, and the size of the BH denominator.
        excluded_reasons: KE_ID → one of the ``EXCLUDED_*`` reasons.
        min_ke_genes: Minimum measured-gene threshold that was applied.
        unresolved_pathways_by_ke: KE_ID → pathway IDs that are mapped to that
            KE but resolved to no genes (issue #81). Flattened into the
            summary's ``unresolved_pathways`` so the sentence can name them.

    Returns:
        Summary dict as documented on ``run_enrichment_analysis``.
    """
    reasons = list(excluded_reasons.values())
    unresolved_pathways = sorted({
        wp
        for wps in (unresolved_pathways_by_ke or {}).values()
        for wp in wps
    })
    return {
        'total_kes': len(set(ke_list)),
        'tested': tested,
        'excluded_no_mapping': sum(1 for r in reasons if r == EXCLUDED_NO_MAPPING),
        'excluded_unresolved_mapping': sum(
            1 for r in reasons if r == EXCLUDED_UNRESOLVED_MAPPING
        ),
        'excluded_too_few_genes': sum(1 for r in reasons if r == EXCLUDED_TOO_FEW_GENES),
        'excluded_error': sum(1 for r in reasons if r == EXCLUDED_ERROR),
        'min_ke_genes': min_ke_genes,
        # Issue #81 — named so the reader can check them against the Builder
        # rather than being told, wrongly, that the Key Event is uncurated.
        'unresolved_pathways': unresolved_pathways,
        # Issue #81 — kept per Key Event as well as flattened, so the network
        # can write the IDs onto the node that lost them and a summary rebuilt
        # from stored network JSON still names them.
        'unresolved_pathways_by_ke': dict(unresolved_pathways_by_ke or {}),
        'excluded_reasons': excluded_reasons,
    }

def build_ke_gene_mapping(
    reference_sets: Dict[str, Set[str]],
    ke_list: Set[str],
    gene_logfc_map: Dict[str, float],
    gene_significance_map: Dict[str, bool],
    gene_pvalue_raw_map: Optional[Dict[str, float]] = None,
    gene_pvalue_adj_map: Optional[Dict[str, float]] = None,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Build mapping of KE to gene data for network visualization.
    Only includes genes that exist in the user's dataset (have expression data).

    Args:
        reference_sets: Dictionary mapping KE_ID to gene sets
        ke_list: Set of KE IDs to include
        gene_logfc_map: Mapping of gene to log2FC values
        gene_significance_map: Mapping of gene to significance status
        gene_pvalue_raw_map: Optional mapping of gene to raw p-value (D-04 passthrough,
            Plan 11-03). When supplied, each per-gene dict gains a ``pvalue_raw`` key;
            genes absent from this map receive ``None`` for that key.
        gene_pvalue_adj_map: Optional mapping of gene to adjusted/FDR p-value
            (D-04 passthrough, Plan 11-03). Same fallback semantics as
            ``gene_pvalue_raw_map``.

    Returns:
        Dict mapping KE_ID to list of gene data dictionaries. Each entry always
        contains ``id``, ``log2FC``, and ``significant``. ``pvalue_raw`` and
        ``pvalue_adj`` keys appear only when the corresponding map is supplied;
        missing-from-map values become ``None``.
    """
    ke_gene_map = {}
    # Plan 11-03: when EITHER pvalue map is supplied, both `pvalue_raw` and
    # `pvalue_adj` keys are emitted on every per-gene entry (with None as the
    # fallback for genes absent from a given map). This keeps the per-gene
    # shape uniform for the downstream JS CSV writer (no per-row "key present?"
    # branching). When NEITHER map is supplied the existing 4-arg call shape
    # is preserved exactly -- only `id`, `log2FC`, and `significant` are emitted.
    emit_pvalue_keys = (gene_pvalue_raw_map is not None) or (gene_pvalue_adj_map is not None)
    raw_lookup = gene_pvalue_raw_map or {}
    adj_lookup = gene_pvalue_adj_map or {}

    for ke, genes in reference_sets.items():
        if ke in ke_list:
            gene_data = []
            for g in genes:
                # Only include genes that exist in the user's dataset
                if g in gene_logfc_map:
                    entry = {
                        "id": g,
                        "log2FC": gene_logfc_map[g],
                        "significant": bool(gene_significance_map.get(g, False)),
                    }
                    if emit_pvalue_keys:
                        entry["pvalue_raw"] = raw_lookup.get(g)
                        entry["pvalue_adj"] = adj_lookup.get(g)
                    gene_data.append(entry)
            ke_gene_map[ke] = gene_data

    total_genes = sum(len(genes) for genes in ke_gene_map.values())
    logger.debug(f"Built KE-gene mapping for {len(ke_gene_map)} KEs with {total_genes} genes (filtered to only include genes with expression data)")
    return ke_gene_map


def run_enrichment(
    method: str,
    df: pd.DataFrame,
    reference_sets: Dict[str, Set[str]],
    ke_list: Set[str],
    ke_title_map: Dict[str, str],
    **kwargs: Any,
) -> pd.DataFrame:
    """Dispatch enrichment to the Fisher or GSEA backend (D-02).

    Routes the call to ``run_enrichment_analysis`` (ORA / Fisher's exact) or
    ``run_gsea_analysis`` (gseapy.prerank) based on ``method``.  The GSEA
    import is lazy — inside the ``method == 'gsea'`` branch — so the Fisher-
    only hot path incurs zero gseapy import cost (mirrors the lazy-import
    pattern at services/data_service.py:200).

    Args:
        method: 'ora' for Fisher's exact (over-representation analysis) or
            'gsea' for rank-based GSEA via gseapy.prerank.
        df: Processed gene expression DataFrame passed through to the backend.
        reference_sets: KE_ID → gene set mapping passed through to the backend.
        ke_list: Set of KE IDs to analyse.
        ke_title_map: KE_ID → human-readable title mapping.
        **kwargs: Additional keyword arguments forwarded to the chosen backend.
            For ORA: ``gene_logfc_map`` drives the Direction column.
            For GSEA: ``gene_logfc_map`` is suppressed by the caller (D-14).
            ``unresolved_ke_pathways`` (issue #81) is accepted by both
            backends, so the exclusion accounting reads identically whichever
            method was used.

    Returns:
        pd.DataFrame from the chosen backend, sorted by FDR ascending.

    Raises:
        ValueError: If ``method`` is not 'ora' or 'gsea' (defence in depth —
            the route-level whitelist check at app.py:/analyze is the first
            layer; this raise is the second layer so direct service-layer
            callers cannot accidentally pass an unknown method).
    """
    if method == 'ora':
        return run_enrichment_analysis(df, reference_sets, ke_list, ke_title_map, **kwargs)
    if method == 'gsea':
        from services.gsea_service import run_gsea_analysis
        # gene_logfc_map is an ORA-only kwarg (drives the Direction column).
        # The /analyze route passes it as None under GSEA (D-14); drop it so
        # it is not forwarded to run_gsea_analysis, which does not accept it.
        kwargs.pop('gene_logfc_map', None)
        return run_gsea_analysis(df, reference_sets, ke_list, ke_title_map, **kwargs)
    raise ValueError(f"Unknown enrichment method: {method!r}. Expected 'ora' or 'gsea'.")