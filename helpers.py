import logging

import pandas as pd

logger = logging.getLogger(__name__)


def load_reference_sets(ke_wp_path='data/KE-WP.csv',
                        wp_gene_path='data/edges_wpid_to_gene.csv',
                        node_path='data/node_attributes.csv',
                        ke_wp_df=None):
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

    Returns
    -------
    dict[str, set[str]]
        Mapping of KE ID -> set of uppercase gene symbols.
    """
    # Use the provided DataFrame or read from the CSV file
    if ke_wp_df is None:
        ke_wp_df = pd.read_csv(ke_wp_path)

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

