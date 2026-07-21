import logging

import pandas as pd

logger = logging.getLogger(__name__)

# Issue #60: minimum KE-mapping confidence. Builder mappings carry a
# confidence_level of High/Medium/Low; the analyser lets users restrict the
# enrichment to better-evidenced mappings. "all" is the default and keeps the
# pre-#60 behaviour (no filtering).
VALID_MIN_CONFIDENCE = ("all", "medium", "high")
DEFAULT_MIN_CONFIDENCE = "all"

# Ordinal rank per confidence level; a mapping qualifies when its rank is at
# least the rank of the requested threshold.
_CONFIDENCE_RANKS = {"low": 1, "medium": 2, "high": 3}

# Record/column keys that may carry the confidence value. The Builder API uses
# ``confidence_level``; the capitalised spellings are tolerated so CSV exports
# of the same data work unchanged.
CONFIDENCE_KEYS = ("confidence_level", "Confidence_Level", "Confidence", "confidence")


def confidence_rank(value):
    """Map a confidence label to its ordinal rank.

    Parameters
    ----------
    value : str or None
        Confidence label such as ``"High"``, ``"medium"`` or ``None``.

    Returns
    -------
    int or None
        ``3`` for High, ``2`` for Medium, ``1`` for Low; ``None`` when the value
        is missing, blank, or not a recognised level. ``None`` means "unknown"
        and callers treat it as a filter no-op (issue #60 graceful degradation).
    """
    if value is None:
        return None
    text = str(value).strip().lower()
    if not text or text in ('nan', 'none'):
        return None
    return _CONFIDENCE_RANKS.get(text)


def record_confidence_rank(record):
    """Return the confidence rank of a mapping record, or None when absent.

    Looks the value up under every known confidence key so both Builder API
    records (``confidence_level``) and CSV-derived records are handled.
    """
    for key in CONFIDENCE_KEYS:
        if key in record:
            rank = confidence_rank(record.get(key))
            if rank is not None:
                return rank
    return None


def filter_records_by_confidence(records, min_confidence=DEFAULT_MIN_CONFIDENCE):
    """Drop KE-mapping records below a minimum confidence threshold (issue #60).

    Records whose confidence is missing, blank, or unrecognised are **kept** —
    the filter is a no-op for them. That keeps the local ``data/KE-WP.csv``
    fallback (no confidence column) and any resource lacking the field fully
    usable, so an analysis is never silently emptied by the filter alone.

    Parameters
    ----------
    records : list[dict]
        Raw mapping records, e.g. from ``fetch_ke_wp_records()`` or
        ``load_ke_wp_records_csv()``.
    min_confidence : str
        One of ``"all"`` (default, no filtering), ``"medium"`` (Medium+High) or
        ``"high"``. Unknown values are treated as ``"all"``.

    Returns
    -------
    list[dict]
        The qualifying records as a new list; the input is not mutated.
    """
    threshold = confidence_rank(min_confidence)
    if min_confidence not in VALID_MIN_CONFIDENCE or threshold is None:
        # "all" (or an unrecognised threshold) -> no filtering.
        return list(records)

    kept = [r for r in records if (record_confidence_rank(r) or threshold) >= threshold]
    logger.info(
        "filter_records_by_confidence: kept %d of %d records at min_confidence=%s",
        len(kept), len(records), min_confidence,
    )
    return kept


def filter_ke_wp_df_by_confidence(ke_wp_df, min_confidence=DEFAULT_MIN_CONFIDENCE):
    """Drop below-threshold rows from a KE-WP DataFrame (issue #60).

    A no-op when ``min_confidence`` is ``"all"``/unknown, or when the frame
    carries no confidence column at all (the ``data/KE-WP.csv`` fallback shape).

    Parameters
    ----------
    ke_wp_df : pandas.DataFrame
        KE-WP mappings, optionally with one of :data:`CONFIDENCE_KEYS` as a column.
    min_confidence : str
        ``"all"``, ``"medium"`` or ``"high"``.

    Returns
    -------
    pandas.DataFrame
        The qualifying rows (the same frame when the filter is a no-op).
    """
    threshold = confidence_rank(min_confidence)
    if min_confidence not in VALID_MIN_CONFIDENCE or threshold is None:
        return ke_wp_df

    column = next((c for c in CONFIDENCE_KEYS if c in ke_wp_df.columns), None)
    if column is None:
        logger.info(
            "filter_ke_wp_df_by_confidence: no confidence column present; "
            "min_confidence=%s applied as a no-op", min_confidence
        )
        return ke_wp_df

    ranks = ke_wp_df[column].map(confidence_rank)
    # Unknown/missing confidence (NaN rank) is kept — graceful degradation.
    keep = ranks.isna() | (ranks >= threshold)
    # .copy() so downstream in-place column normalisation never writes to a slice.
    filtered = ke_wp_df[keep].copy()
    logger.info(
        "filter_ke_wp_df_by_confidence: kept %d of %d KE-WP rows at min_confidence=%s",
        len(filtered), len(ke_wp_df), min_confidence,
    )
    return filtered


def load_reference_sets(ke_wp_path='data/KE-WP.csv',
                        wp_gene_path='data/edges_wpid_to_gene.csv',
                        node_path='data/node_attributes.csv',
                        ke_wp_df=None,
                        min_confidence=DEFAULT_MIN_CONFIDENCE):
    """Build KE-to-gene reference sets from local CSV files or a pre-built DataFrame.

    Parameters
    ----------
    ke_wp_path : str
        Path to the KE-WP mapping CSV file. Ignored when ``ke_wp_df`` is provided.
    wp_gene_path : str
        Path to the WikiPathways-to-gene edge CSV file.
    node_path : str
        Path to the node attributes CSV file (gene ID to symbol mapping).
    ke_wp_df : pandas.DataFrame or None
        Optional pre-built DataFrame with columns ``KE_ID`` and ``WP_ID``.
        When provided, the CSV at ``ke_wp_path`` is not read. KE IDs must
        already be normalised to the ``"KE:55"`` colon format.
    min_confidence : str
        Issue #60 minimum mapping confidence — ``"all"`` (default, no
        filtering), ``"medium"`` or ``"high"``. Rows whose confidence is below
        the threshold are dropped before the WP-to-gene merge, so a KE keeps
        only the genes of its qualifying pathways and a KE whose mappings are
        all below threshold is absent from the result. The filter is a no-op
        when the mappings carry no confidence column (CSV fallback).

    Returns
    -------
    dict[str, set[str]]
        Mapping of KE ID -> set of uppercase gene symbols.
    """
    # Use the provided DataFrame or read from the CSV file
    if ke_wp_df is None:
        ke_wp_df = pd.read_csv(ke_wp_path)

    # Issue #60: apply the minimum-confidence threshold before any merging.
    ke_wp_df = filter_ke_wp_df_by_confidence(ke_wp_df, min_confidence)

    wp_gene_df = pd.read_csv(wp_gene_path)
    node_df = pd.read_csv(node_path)

    logger.debug("Sample gene_id values from wp_gene_df: %s", wp_gene_df['gene_id'].dropna().unique()[:10])
    logger.debug("Sample GeneID values from node_df: %s", node_df['GeneID'].dropna().unique()[:10])

    # Normalize formatting
    ke_wp_df['WP_ID'] = ke_wp_df['WP_ID'].astype(str).str.upper().str.strip()
    wp_gene_df['WPID'] = wp_gene_df['WPID'].astype(str).str.upper().str.strip()
    wp_gene_df['gene_id'] = wp_gene_df['gene_id'].astype(str).str.strip()
    node_df['GeneID'] = (
        node_df['GeneID']
        .dropna()                       # remove NaNs
        .astype(float)                 # ensure float first
        .astype(int)                   # then cast to int safely
        .astype(str)                   # convert to string
        .str.strip()                   # strip whitespace
    )


    node_df['GeneName'] = node_df['GeneName'].astype(str).str.strip()

    logger.debug("KE-WP mappings: %s", ke_wp_df.shape)
    logger.debug("Gene-pathway mappings: %s", wp_gene_df.shape)
    logger.debug("Node annotation mappings: %s", node_df.shape)

    # Merge gene → annotation
    wp_gene_annotated = wp_gene_df.merge(
        node_df[['GeneID', 'GeneName']],
        left_on='gene_id',
        right_on='GeneID',
        how='left'
    )

    logger.debug("After annotation merge: %s", wp_gene_annotated.shape)

    # Merge WP → KE
    merged = ke_wp_df.merge(
        wp_gene_annotated,
        left_on='WP_ID',
        right_on='WPID',
        how='inner'
    )

    pre_dropna_count = len(merged)
    logger.debug("Final merged KE-gene mapping before dropna: %s", merged.shape)

    # Drop entries without gene symbol
    merged = merged.dropna(subset=['GeneName'])

    logger.debug("Final merged KE-gene mapping after dropna: %s", merged.shape)

    # Report any KE-WP mappings lost in the inner join
    dropped_ke_wp = len(ke_wp_df) - merged['KE_ID'].nunique()
    if dropped_ke_wp > 0:
        logger.warning(
            "Inner join dropped %d of %d KE-WP mappings (no matching WikiPathways genes)",
            dropped_ke_wp, len(ke_wp_df)
        )

    # Group into dict: KE_ID → set of gene symbols
    reference_sets = (
        merged.groupby('KE_ID')['GeneName']
        .apply(lambda g: set(g.str.upper()))
        .to_dict()
    )

    return reference_sets

