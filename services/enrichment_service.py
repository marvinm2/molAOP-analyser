"""
Gene set enrichment analysis service.
"""
import pandas as pd
import logging
from typing import Dict, Set, List, Any, Tuple
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

def run_enrichment_analysis(
    df: pd.DataFrame, 
    reference_sets: Dict[str, Set[str]], 
    ke_list: Set[str],
    ke_title_map: Dict[str, str]
) -> pd.DataFrame:
    """
    Run Fisher's exact test enrichment analysis for Key Events.
    
    Args:
        df: Processed gene expression dataframe with 'significant' column
        reference_sets: Dictionary mapping KE_ID to gene sets
        ke_list: Set of KE IDs to analyze
        ke_title_map: Mapping of KE IDs to titles
    
    Returns:
        pd.DataFrame: Enrichment results sorted by FDR
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

    for ke, ref_genes in filtered_reference_sets.items():
        try:
            # Find overlap between KE genes and user genes (ref_genes already normalized)
            ke_genes = ref_genes & all_genes

            if not ke_genes:
                logger.debug(f"No overlap found for KE {ke}")
                continue

            # Skip KEs with too few genes for reliable statistics
            if len(ke_genes) < 5:
                logger.debug(f"Skipping KE {ke}: only {len(ke_genes)} genes (minimum 5 required)")
                continue

            sig_in_ke = {g for g in ke_genes if user_gene_status.get(g, False)}
            non_sig_in_ke = ke_genes - sig_in_ke

            # Create 2x2 contingency table
            a = len(sig_in_ke)                    # significant & in KE
            b = sig_genes - a                     # significant & not in KE
            c = len(non_sig_in_ke)                # non-significant & in KE
            d = non_sig_genes - c                 # non-significant & not in KE

            # Run Fisher's exact test (one-tailed, greater)
            odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

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
            r_odds.append(round(odds, 4) if not pd.isna(odds) else 'NA')
            r_pval.append(pval)
            r_overlap.append(", ".join(sorted(sig_in_ke)))
            r_noverlap.append(a)

        except Exception as e:
            logger.error(f"Error processing KE {ke}: {e}")
            continue

    if not r_ke:
        logger.error("No enrichment results generated")
        raise ValueError("No enrichment results generated")

    # Build DataFrame from columnar data and apply FDR correction
    df_results = pd.DataFrame({
        'Title': r_title,
        'p_value': r_pval,
        'FDR': multipletests(r_pval, method="fdr_bh")[1],
        'num_overlap': r_noverlap,
        'pct_sig_in_KE': r_pct,
        'total_KE_genes_in_dataset': r_total,
        'odds_ratio': r_odds,
        'overlap': r_overlap,
        'KE': r_ke,
        'sig_in_KE': r_sig_in,
        'sig_not_KE': r_sig_not,
        'non_sig_in_KE': r_nsig_in,
        'non_sig_not_KE': r_nsig_not,
    })
    df_results = df_results.sort_values("FDR")

    logger.info(f"Enrichment analysis completed: {len(df_results)} results generated")

    return df_results

def build_ke_gene_mapping(
    reference_sets: Dict[str, Set[str]], 
    ke_list: Set[str], 
    gene_logfc_map: Dict[str, float],
    gene_significance_map: Dict[str, bool]
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Build mapping of KE to gene data for network visualization.
    Only includes genes that exist in the user's dataset (have expression data).
    
    Args:
        reference_sets: Dictionary mapping KE_ID to gene sets
        ke_list: Set of KE IDs to include
        gene_logfc_map: Mapping of gene to log2FC values
        gene_significance_map: Mapping of gene to significance status
    
    Returns:
        Dict mapping KE_ID to list of gene data dictionaries
    """
    ke_gene_map = {}
    
    for ke, genes in reference_sets.items():
        if ke in ke_list:
            gene_data = []
            for g in genes:
                # Only include genes that exist in the user's dataset
                if g in gene_logfc_map:
                    gene_data.append({
                        "id": g,
                        "log2FC": gene_logfc_map[g],
                        "significant": bool(gene_significance_map.get(g, False))
                    })
            ke_gene_map[ke] = gene_data
    
    total_genes = sum(len(genes) for genes in ke_gene_map.values())
    logger.debug(f"Built KE-gene mapping for {len(ke_gene_map)} KEs with {total_genes} genes (filtered to only include genes with expression data)")
    return ke_gene_map