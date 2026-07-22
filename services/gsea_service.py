"""
Gene set enrichment analysis (GSEA) service.

Provides gseapy.prerank-based enrichment analysis as an alternative to the
Fisher's exact test path in enrichment_service.py. The dispatcher
(enrichment_service.run_enrichment) routes to this module based on the
user's method selection.
"""
import math
import logging
from typing import Any, Dict, Optional, Set

import numpy as np
import pandas as pd
import gseapy as gp

from services.enrichment_service import (
    EXCLUDED_NO_MAPPING,
    EXCLUDED_TOO_FEW_GENES,
    EXCLUDED_UNRESOLVED_MAPPING,
    KE_SUMMARY_ATTR,
    _build_ke_summary,
    normalise_unresolved_ke_pathways,
)

logger = logging.getLogger(__name__)


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


def run_gsea_analysis(
    df: pd.DataFrame,
    reference_sets: Dict[str, Set[str]],
    ke_list: Set[str],
    ke_title_map: Dict[str, str],
    *,
    min_size: int = 5,
    max_size: int = 1000,
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
        max_size: Maximum gene-set size passed to gseapy (keyword-only).
        permutation_num: Number of permutations for the null distribution
            (keyword-only).
        seed: Random seed for reproducibility (keyword-only).
        unresolved_ke_pathways: Optional KE_ID -> pathway IDs that are mapped
            to that KE but could not be resolved to genes (issue #81). Same
            meaning and effect as in ``run_enrichment_analysis`` — ORA and
            GSEA report exclusions with one vocabulary.

    Returns:
        pd.DataFrame sorted by FDR ascending with columns:
        ['KE', 'Title', 'NES', 'p_value', 'FDR', 'lead_genes',
         'total_KE_genes_in_dataset']

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
    res['NES'] = res['NES'].clip(lower=-3, upper=3)

    # Populate total_KE_genes_in_dataset from pre-computed overlap counts
    res['total_KE_genes_in_dataset'] = (
        res['KE'].map(kes_to_overlap_count).fillna(0).astype(int)
    )

    # Log KEs dropped by min_size filter for parity with Fisher per-KE skip log
    # (RESEARCH Open Question #2) — diff the input ke_list against res['KE']
    dropped = (set(ke_list) & set(reference_sets.keys())) - set(res['KE'])
    if dropped:
        logger.info(
            f"GSEA dropped {len(dropped)} KEs due to min_size={min_size} filter: "
            f"{sorted(dropped)}"
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
    for ke in dropped:
        if not gene_sets.get(ke):
            excluded_reasons[ke] = EXCLUDED_UNRESOLVED_MAPPING
            if unresolved_map.get(ke):
                unresolved_named[ke] = unresolved_map[ke]
        else:
            excluded_reasons[ke] = EXCLUDED_TOO_FEW_GENES
    ke_summary = _build_ke_summary(
        ke_list, len(res), excluded_reasons, min_size,
        unresolved_pathways_by_ke=unresolved_named,
    )

    # Sort by FDR ascending; use KE ID as tie-breaker so row order is fully
    # deterministic regardless of the internal ordering gseapy returns.
    res = res.sort_values(['FDR', 'KE'], ascending=[True, True], kind='mergesort')

    out = res[['KE', 'Title', 'NES', 'p_value', 'FDR', 'lead_genes',
               'total_KE_genes_in_dataset']]
    # Attach on the returned object: pandas only propagates .attrs through
    # operations that call __finalize__, so set it last.
    out.attrs[KE_SUMMARY_ATTR] = ke_summary
    return out
