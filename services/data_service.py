"""
Data processing service for gene expression data and reference sets.
"""
import pandas as pd
import logging
import re
from typing import Dict, Set, Tuple, List, Any, Optional
from scipy.stats import combine_pvalues
from config import Config
from exceptions import FileProcessingError, DataValidationError, AOPDataError

logger = logging.getLogger(__name__)

# Placeholders that pandas/R/Excel write for a missing gene symbol. Compared
# against the upper-cased, stripped identifier (issue #80).
_MISSING_ID_TOKENS = {'NAN', 'NA', 'NONE', 'NULL', '#N/A', '-', '.'}

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
        
        # Expand if symbols are separated by '///'.
        #
        # Issue #80: rows with no gene symbol used to survive this loop as the
        # literal identifier "NAN" (str(NaN).upper()), so every symbol-less row
        # in the upload collapsed into a single pseudo-gene that inflated the
        # background by one and, if it passed the significance threshold, took
        # part in every Fisher test while being unable to overlap any KE gene
        # set. Such rows are dropped up front and counted, so the discard is
        # reported rather than silently renamed.
        df_expanded = []
        dropped_rows = 0
        for _, row in df_raw.iterrows():
            raw_id = row[id_col]
            genes = [g.strip().upper() for g in str(raw_id).split('///')]
            genes = [g for g in genes if g and g not in _MISSING_ID_TOKENS]
            if not genes:
                dropped_rows += 1
                continue
            for gene in genes:
                df_expanded.append({
                    'ID': gene,
                    'log2FC': row[fc_col],
                    'pval': row[pval_col]
                })

        df = pd.DataFrame(df_expanded, columns=['ID', 'log2FC', 'pval']).dropna()

        # Carried through process_gene_expression() into its stats dict so the
        # count can be shown next to the background size.
        df.attrs['dropped_unidentified_rows'] = dropped_rows

        if dropped_rows:
            logger.info(
                f"Dropped {dropped_rows} row(s) without a usable gene identifier "
                f"in column '{id_col}'"
            )
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

def process_gene_expression(df: pd.DataFrame, logfc_threshold: float = 0.0, pval_threshold: Optional[float] = None) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Process gene expression data and combine duplicates.

    Args:
        df: Raw gene expression dataframe
        logfc_threshold: Threshold for significant fold change
        pval_threshold: Optional p-value threshold for the significance flag.
            When None, falls back to Config.PVAL_CUTOFF (the legacy default 0.05).
            Out-of-range or unparseable values should be coerced to None by the
            caller before invocation.

    Returns:
        Tuple of (processed_dataframe, stats_dict)
    """
    # Combine duplicate genes
    grouped = df.groupby('ID')
    combined_rows = []

    for gene, group in grouped:
        if len(group) == 1:
            # Convert Series to dict to maintain consistency
            row = group.iloc[0].to_dict()
            combined_rows.append(row)
        else:
            # Average log2FC and combine p-values using Fisher's method
            avg_fc = group['log2FC'].astype(float).mean()
            _, combined_p = combine_pvalues(group['pval'].astype(float), method='fisher')
            combined_rows.append({'ID': gene, 'log2FC': avg_fc, 'pval': combined_p})

    df_processed = pd.DataFrame(combined_rows)
    df_processed['ID'] = df_processed['ID'].astype(str).str.strip().str.upper()

    # Resolve effective p-value cutoff (user override or legacy default)
    effective_pval = pval_threshold if pval_threshold is not None else Config.PVAL_CUTOFF

    # Add significance flag
    df_processed['significant'] = (
        (df_processed['log2FC'].abs() >= logfc_threshold) &
        (df_processed['pval'] <= effective_pval)
    )

    # Calculate statistics
    n_sig = df_processed['significant'].sum()
    n_total = len(df_processed)
    n_non_sig = n_total - n_sig

    stats = {
        'total_genes': n_total,
        'significant_genes': n_sig,
        'non_significant_genes': n_non_sig,
        'logfc_threshold': logfc_threshold,
        'pval_threshold': effective_pval,
        # Issue #80: rows the loader discarded for having no gene identifier.
        # Reported alongside the background size so 'total_genes' can be quoted
        # as the number of measured genes and the discard stays visible.
        'dropped_unidentified_rows': int(df.attrs.get('dropped_unidentified_rows', 0)),
    }


    logger.info(
        f"Processed {stats['total_genes']} unique genes: {stats['significant_genes']} significant, "
        f"{stats['non_significant_genes']} non-significant "
        f"({stats['dropped_unidentified_rows']} row(s) dropped without a gene identifier)"
    )
    
    if n_sig == 0:
        logger.warning("No significant genes found. Check thresholds.")
    
    return df_processed, stats

def compute_logfc_percentiles(logfc: pd.Series, percentages: Tuple[int, ...] = (10, 20)) -> Dict[int, float]:
    """
    Compute |log2FC| thresholds for the "Top N%" quick-threshold buttons.

    Issue #64: this used to be done in the browser over the volcano payload,
    which app.py truncates to Config.MAX_GENES_DISPLAY. Any upload larger than
    that cap therefore had its percentile taken over an arbitrary subset. The
    calculation belongs server-side, where the full gene list is available.

    The maths deliberately mirrors the original client-side implementation
    (sort |log2FC| descending, index at N% of the row count) so results are
    unchanged for datasets that fit under the cap.

    Note this sets a fold-change threshold only -- the p-value criterion is
    still applied on top of it by process_gene_expression(), so the resulting
    DEG count is not N% of genes.

    Args:
        logfc: Series of log2 fold-change values (the full dataset).
        percentages: Percentages to compute thresholds for.

    Returns:
        Dict mapping each percentage to its rounded |log2FC| threshold. Empty
        if there is no usable data.
    """
    abs_logfc = pd.to_numeric(logfc, errors='coerce').abs().dropna()
    if abs_logfc.empty:
        return {}

    ordered = abs_logfc.sort_values(ascending=False).to_numpy()

    thresholds = {}
    for pct in percentages:
        index = int(len(ordered) * (pct / 100))
        # Guard the boundary: at 100%, or on very small datasets, the computed
        # index can land past the last element.
        index = min(index, len(ordered) - 1)
        thresholds[pct] = round(float(ordered[index]), 2)

    return thresholds

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

def _load_aop_data_csv(aop_id: str) -> Tuple[Set[str], pd.DataFrame, Dict[str, str], Dict[str, str]]:
    """Load AOP data from local CSV files.

    Args:
        aop_id: AOP identifier (e.g. "AOP:1")

    Returns:
        Tuple of (ke_list, edges_df, ke_type_map, ke_title_map)
    """
    aop_ke_map = pd.read_csv('data/aop_ke_map.csv')
    ke_list = set(aop_ke_map[aop_ke_map['AOP_ID'] == aop_id]['KE_ID'].dropna())

    aop_ker_df = pd.read_csv('data/aop_ker_edges.csv')
    edges = aop_ker_df[aop_ker_df['AOP_ID'] == aop_id]

    ke_metadata_df = pd.read_csv("data/ke_metadata.csv")
    ke_type_map = dict(zip(ke_metadata_df["KE_ID"], ke_metadata_df["Type"]))
    ke_title_map = dict(zip(ke_metadata_df["KE_ID"], ke_metadata_df["Title"]))

    return ke_list, edges, ke_type_map, ke_title_map


def load_aop_data(aop_id: str) -> Tuple[Set[str], pd.DataFrame, Dict[str, str], Dict[str, str]]:
    """
    Load AOP-specific data including KE list, edges, and metadata.

    Routes to local CSV files or AOP-Wiki SPARQL endpoint based on the
    ``source`` field in the CASE_STUDY_AOPS config entry for this AOP.

    Args:
        aop_id: AOP identifier (e.g. "AOP:105" or "NETWORK:kidney")

    Returns:
        Tuple of (ke_list, edges_df, ke_type_map, ke_title_map)
    """
    try:
        # Look up config entry for this AOP
        aop_config = None
        for entry in Config.CASE_STUDY_AOPS.values():
            if entry.get("id") == aop_id:
                aop_config = entry
                break

        source = (aop_config or {}).get("source", "sparql")

        if source == "sparql":
            from services.sparql_service import fetch_aop_ke_data_cached

            # Combined network: use the aop_ids list from config
            if aop_id.startswith("NETWORK:") and aop_config and aop_config.get("aop_ids"):
                all_ids = aop_config["aop_ids"]

                # Split IDs by source: some AOPs are CSV-only
                sparql_ids = []
                csv_ids = []
                for aid in all_ids:
                    aid_source = "sparql"
                    for entry in Config.CASE_STUDY_AOPS.values():
                        if entry.get("id") == aid:
                            aid_source = entry.get("source", "sparql")
                            break
                    if aid_source == "csv":
                        csv_ids.append(aid)
                    else:
                        sparql_ids.append(aid)

                # Fetch SPARQL AOPs
                ke_list, edges, ke_type_map, ke_title_map = fetch_aop_ke_data_cached(sparql_ids)

                # Merge in CSV-sourced AOPs
                for csv_aop_id in csv_ids:
                    csv_kes, csv_edges, csv_type_map, csv_title_map = _load_aop_data_csv(csv_aop_id)
                    ke_list |= csv_kes
                    edges = pd.concat([edges, csv_edges], ignore_index=True)
                    ke_type_map.update(csv_type_map)
                    ke_title_map.update(csv_title_map)
                    logger.info(f"Merged CSV AOP {csv_aop_id}: +{len(csv_kes)} KEs, +{len(csv_edges)} edges")
            else:
                sparql_ids = [aop_id]
                ke_list, edges, ke_type_map, ke_title_map = fetch_aop_ke_data_cached(sparql_ids)

            logger.info(f"Loaded AOP {aop_id} via SPARQL: {len(ke_list)} KEs, {len(edges)} edges")
            return ke_list, edges, ke_type_map, ke_title_map

        # Default: CSV path
        ke_list, edges, ke_type_map, ke_title_map = _load_aop_data_csv(aop_id)
        logger.info(f"Loaded AOP {aop_id} from CSV: {len(ke_list)} KEs, {len(edges)} edges")
        return ke_list, edges, ke_type_map, ke_title_map

    except Exception as e:
        logger.error(f"Failed to load AOP data for {aop_id}: {e}")
        raise