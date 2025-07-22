import pandas as pd

def load_reference_sets(ke_wp_path, wp_gene_path, node_path):
    # Load input files
    ke_wp_df = pd.read_csv(ke_wp_path)
    wp_gene_df = pd.read_csv(wp_gene_path)
    node_df = pd.read_csv(node_path)

    print("\n[DEBUG] Sample gene_id values from wp_gene_df:")
    print(wp_gene_df['gene_id'].dropna().unique()[:10])

    print("\n[DEBUG] Sample GeneID values from node_df:")
    print(node_df['GeneID'].dropna().unique()[:10])

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

    print("[DEBUG] KE-WP mappings:", ke_wp_df.shape)
    print("[DEBUG] Gene-pathway mappings:", wp_gene_df.shape)
    print("[DEBUG] Node annotation mappings:", node_df.shape)

    # Merge gene → annotation
    wp_gene_annotated = wp_gene_df.merge(
        node_df[['GeneID', 'GeneName']],
        left_on='gene_id',
        right_on='GeneID',
        how='left'
    )

    print("[DEBUG] After annotation merge:", wp_gene_annotated.shape)

    # Merge WP → KE
    merged = ke_wp_df.merge(
        wp_gene_annotated,
        left_on='WP_ID',
        right_on='WPID',
        how='inner'
    )

    print("[DEBUG] Final merged KE-gene mapping before dropna:", merged.shape)

    # Drop entries without gene symbol
    merged = merged.dropna(subset=['GeneName'])

    print("[DEBUG] Final merged KE-gene mapping after dropna:", merged.shape)

    # Group into dict: KE_ID → set of gene symbols
    reference_sets = (
        merged.groupby('KE_ID')['GeneName']
        .apply(lambda g: set(g.str.upper()))
        .to_dict()
    )

    return reference_sets

