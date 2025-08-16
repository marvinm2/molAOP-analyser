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
    filtered_reference_sets = {ke: genes for ke, genes in reference_sets.items() if ke in ke_list}
    
    if not filtered_reference_sets:
        logger.error("No reference gene sets found for the selected AOP")
        raise ValueError("No reference gene sets found for the selected AOP")
    
    # Define gene universe and significance status
    all_genes = set(df['ID'])
    user_gene_status = dict(zip(df['ID'], df['significant']))
    
    sig_genes = df['significant'].sum()
    non_sig_genes = df.shape[0] - sig_genes
    
    logger.info(f"Gene universe: {len(all_genes)} total, {sig_genes} significant, {non_sig_genes} non-significant")
    
    results = []
    
    for ke, ref_genes in filtered_reference_sets.items():
        try:
            # Find overlap between KE genes and user genes
            ke_genes = {g.strip().upper() for g in ref_genes} & all_genes
            
            if not ke_genes:
                logger.debug(f"No overlap found for KE {ke}")
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
            
            logger.debug(f"KE {ke}: {len(ke_genes)} total genes, {a} significant overlap")
            
            results.append({
                'KE': ke,
                'Title': ke_title_map.get(ke, ke),
                'total_KE_genes_in_dataset': len(ke_genes),
                'sig_in_KE': a,
                'sig_not_KE': b,
                'non_sig_in_KE': c,
                'non_sig_not_KE': d,
                'pct_sig_in_KE': round((a / len(ke_genes)) * 100, 1) if ke_genes else 0,
                'odds_ratio': round(odds, 4) if not pd.isna(odds) else 'NA',
                'p_value': pval,
                'overlap': ", ".join(sorted(sig_in_ke)),
                'num_overlap': a
            })
            
        except Exception as e:
            logger.error(f"Error processing KE {ke}: {e}")
            continue
    
    if not results:
        logger.error("No enrichment results generated")
        raise ValueError("No enrichment results generated")
    
    # Convert to DataFrame and apply FDR correction
    df_results = pd.DataFrame(results)
    df_results["FDR"] = multipletests(df_results["p_value"], method="fdr_bh")[1]
    df_results = df_results.sort_values("FDR")
    
    # Reorder columns to prioritize most important information first
    column_order = [
        'Title',               # KE title (most important for user understanding)
        'p_value',            # Statistical significance 
        'FDR',                # Multiple testing corrected p-value
        'num_overlap',        # Number of overlapping genes
        'pct_sig_in_KE',      # Percentage enrichment
        'total_KE_genes_in_dataset',  # KE gene set size
        'odds_ratio',         # Effect size
        'overlap',            # Specific genes (detailed view)
        'KE',                 # KE ID (for reference)
        'sig_in_KE',          # Contingency table details
        'sig_not_KE',
        'non_sig_in_KE', 
        'non_sig_not_KE'
    ]
    
    df_results = df_results[column_order]
    
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
                gene_data.append({
                    "id": g,
                    "log2FC": gene_logfc_map.get(g, 0),
                    "significant": bool(gene_significance_map.get(g, False))
                })
            ke_gene_map[ke] = gene_data
    
    logger.debug(f"Built KE-gene mapping for {len(ke_gene_map)} KEs")
    return ke_gene_map