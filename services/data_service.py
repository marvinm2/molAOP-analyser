"""
Data processing service for gene expression data and reference sets.
"""
import pandas as pd
import logging
import re
from typing import Dict, Set, Tuple, List, Any
from scipy.stats import combine_pvalues
from config import Config
from exceptions import FileProcessingError, DataValidationError, AOPDataError

logger = logging.getLogger(__name__)

def load_and_validate_data(filepath: str, id_col: str, fc_col: str, pval_col: str) -> pd.DataFrame:
    """
    Load and validate gene expression data from file.
    
    Args:
        filepath: Path to the data file
        id_col: Column name for gene IDs
        fc_col: Column name for log2 fold change
        pval_col: Column name for p-values
    
    Returns:
        pd.DataFrame: Cleaned and validated dataframe
    
    Raises:
        ValueError: If required columns are missing or data is invalid
    """
    try:
        # Read file
        df_raw = pd.read_csv(filepath, sep=None, engine='python')
        logger.info(f"Loaded {len(df_raw)} rows from {filepath}")
        
        # Validate columns exist
        missing_cols = [col for col in [id_col, fc_col, pval_col] if col not in df_raw.columns]
        if missing_cols:
            raise DataValidationError(
                f"Missing required columns: {missing_cols}",
                field="columns",
                value=str(missing_cols)
            )
        
        # Expand if symbols are separated by '///'
        df_expanded = []
        for _, row in df_raw.iterrows():
            genes = str(row[id_col]).split('///')
            for gene in genes:
                df_expanded.append({
                    'ID': gene.strip().upper(),
                    'log2FC': row[fc_col],
                    'pval': row[pval_col]
                })
        
        df = pd.DataFrame(df_expanded).dropna()
        logger.info(f"Expanded to {len(df)} gene entries after processing")
        
        return df
        
    except DataValidationError:
        raise  # Re-raise validation errors as-is
    except pd.errors.EmptyDataError:
        raise FileProcessingError(f"File is empty: {filepath}", filename=filepath)
    except pd.errors.ParserError as e:
        raise FileProcessingError(f"Failed to parse file: {e}", filename=filepath)
    except Exception as e:
        logger.error(f"Failed to load data from {filepath}: {e}")
        raise FileProcessingError(f"Unexpected error loading file: {e}", filename=filepath)

def process_gene_expression(df: pd.DataFrame, logfc_threshold: float = 0.0) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Process gene expression data and combine duplicates.
    
    Args:
        df: Raw gene expression dataframe
        logfc_threshold: Threshold for significant fold change
    
    Returns:
        Tuple of (processed_dataframe, stats_dict)
    """
    # Combine duplicate genes
    grouped = df.groupby('ID')
    combined_rows = []
    
    for gene, group in grouped:
        if len(group) == 1:
            combined_rows.append(group.iloc[0])
        else:
            # Average log2FC and combine p-values using Fisher's method
            avg_fc = group['log2FC'].astype(float).mean()
            _, combined_p = combine_pvalues(group['pval'].astype(float), method='fisher')
            combined_rows.append({'ID': gene, 'log2FC': avg_fc, 'pval': combined_p})
    
    df_processed = pd.DataFrame(combined_rows)
    df_processed['ID'] = df_processed['ID'].astype(str).str.strip().str.upper()
    
    # Add significance flag
    df_processed['significant'] = (
        (df_processed['log2FC'].abs() >= logfc_threshold) & 
        (df_processed['pval'] <= Config.PVAL_CUTOFF)
    )
    
    # Calculate statistics
    n_sig = df_processed['significant'].sum()
    n_total = len(df_processed)
    n_non_sig = n_total - n_sig
    
    stats = {
        'total_genes': n_total,
        'significant_genes': n_sig,
        'non_significant_genes': n_non_sig,
        'logfc_threshold': logfc_threshold
    }
    
    logger.info(f"Processed {stats['total_genes']} unique genes: {stats['significant_genes']} significant, {stats['non_significant_genes']} non-significant")
    
    if n_sig == 0:
        logger.warning("No significant genes found. Check thresholds.")
    
    return df_processed, stats

def guess_id_type(gene_series: pd.Series) -> str:
    """
    Guess the type of gene identifiers.
    
    Args:
        gene_series: Series of gene identifiers
    
    Returns:
        str: Predicted ID type (HGNC, Ensembl, or Entrez)
    """
    hgnc_pattern = re.compile(r"^[A-Z0-9\-]+$")
    ensembl_pattern = re.compile(r"ENSG\d+")
    entrez_pattern = re.compile(r"^\d+$")
    
    hgnc, ensembl, entrez = 0, 0, 0
    
    for gene in gene_series.head(20):
        if ensembl_pattern.match(str(gene)):
            ensembl += 1
        elif entrez_pattern.match(str(gene)):
            entrez += 1
        elif hgnc_pattern.match(str(gene)):
            hgnc += 1
    
    scores = {"HGNC": hgnc, "Ensembl": ensembl, "Entrez": entrez}
    predicted_type = max(scores, key=scores.get)
    
    logger.info(f"Predicted gene ID type: {predicted_type} (scores: {scores})")
    return predicted_type

def load_aop_data(aop_id: str) -> Tuple[Set[str], pd.DataFrame, Dict[str, str], Dict[str, str]]:
    """
    Load AOP-specific data including KE list, edges, and metadata.
    
    Args:
        aop_id: AOP identifier
    
    Returns:
        Tuple of (ke_list, edges_df, ke_type_map, ke_title_map)
    """
    try:
        # Load AOP-KE mappings
        aop_ke_map = pd.read_csv('data/aop_ke_map.csv')
        ke_list = set(aop_ke_map[aop_ke_map['AOP_ID'] == aop_id]['KE_ID'].dropna())
        
        # Load AOP edges
        aop_ker_df = pd.read_csv('data/aop_ker_edges.csv')
        edges = aop_ker_df[aop_ker_df['AOP_ID'] == aop_id]
        
        # Load KE metadata
        ke_metadata_df = pd.read_csv("data/ke_metadata.csv")
        ke_type_map = dict(zip(ke_metadata_df["KE_ID"], ke_metadata_df["Type"]))
        ke_title_map = dict(zip(ke_metadata_df["KE_ID"], ke_metadata_df["Title"]))
        
        logger.info(f"Loaded AOP {aop_id} data: {len(ke_list)} KEs, {len(edges)} edges")
        
        return ke_list, edges, ke_type_map, ke_title_map
        
    except Exception as e:
        logger.error(f"Failed to load AOP data for {aop_id}: {e}")
        raise