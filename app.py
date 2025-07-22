from flask import Flask, render_template, request, redirect, url_for
import pandas as pd
import os
from scipy.stats import fisher_exact, combine_pvalues
from statsmodels.stats.multitest import multipletests
import json
import math

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)  # Create uploads folder if missing
from werkzeug.utils import secure_filename

ALLOWED_EXTENSIONS = {'csv', 'tsv', 'txt'}
#app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024  # 2 MB max
app.config['MAX_CONTENT_LENGTH'] = 10 * 1024 * 1024  # 10 MB


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

from flask import abort
from werkzeug.exceptions import RequestEntityTooLarge

@app.errorhandler(RequestEntityTooLarge)
def handle_file_too_large(e):
    return "File too large. Maximum upload size is 10 MB.", 413

CASE_STUDY_AOPS = {
    "DEMO": {"label": "---DEMO---", "enabled": False},
    "steatosis": {"id": "AOP:1", "label": "PXR activation leading to liver steatosis", "enabled": True},
    "VHP-CASES:": {"label": "---VHP CASES---", "enabled": False},
    "vhp-kidney": {"id": "AOP:2", "label": "DNA adduct formation leading to kidney failure", "enabled": False},
    "vhp-parkinson": {"id": "AOP:3", "label": "Calcium overload in dopaminergic neurons of the substantia nigra leading to parkinsonian motor deficits", "enabled": False},
    "vhp-thyroid": {"id": "AOP:4", "label": "Thyroid hormone-mediated neurodevelopmental toxicity", "enabled": False},
    "ORGAN-NETWORK": {"label": "---ORGAN NETWORKS---","enabled": False},
    "Liver-aop-network": {"id": "AOP:5", "label": "Liver AOP network", "enabled": False},
    "Brain-aop-network": {"id": "AOP:6", "label": "Brain AOP network", "enabled": False},
    "Kidney-aop-network": {"id": "AOP:7", "label": "Kidney AOP network", "enabled": False},
    "Lung-aop-network": {"id": "AOP:8", "label": "Lung AOP network", "enabled": False},
}


@app.route('/')
def index():
    return render_template('index.html', case_study_aops=CASE_STUDY_AOPS)

@app.route('/preview', methods=['POST'])
# def preview():
#     # Check if a file is uploaded or if we’re re-using one
#     file = request.files.get('gene_file')
#     filename = request.form.get("filename")
#     
#     if file and allowed_file(file.filename):
#         # New upload
#         filename = secure_filename(file.filename)
#         filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
#         file.save(filepath)
#     elif filename:
#         # Reusing existing file
#         filepath = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(filename))
#         if not os.path.exists(filepath):
#             return "Uploaded file not found.", 400
#     else:
#         return "No file uploaded or referenced.", 400
# 
#     # Load the file
#     try:
#         df = pd.read_csv(filepath, sep=None, engine='python')
#     except Exception as e:
#         return f"Failed to read file: {e}", 400
# 
#     preview_rows = df.head(5).to_dict(orient="records")
#     columns = df.columns.tolist()
#     # Assuming you already have: df = pd.read_csv(filepath, ...)
#     columns = df.columns.tolist()
#     preview_rows = df.head(5).to_dict(orient="records")
# 
#     # Prepare volcano data if all 3 columns are selected
#     id_col = request.form.get("id_column") or None
#     fc_col = request.form.get("fc_column") or None
#     pval_col = request.form.get("pval_column") or None
# 
#     logfc_threshold = request.form.get("logfc_threshold")
#     try:
#         logfc_threshold = float(logfc_threshold)
#     except (TypeError, ValueError):
#         logfc_threshold = None
#     import math
#     pval_cutoff = 0.05  # fixed for now
#     pval_y = -math.log10(pval_cutoff)
#     volcano_data = []
#     if id_col and fc_col and pval_col:
#         try:
#             df = df[[id_col, fc_col, pval_col]].dropna()
#             df.columns = ['ID', 'log2FC', 'pval']
#             df['ID'] = df['ID'].astype(str)
#             df['log2FC'] = pd.to_numeric(df['log2FC'], errors='coerce')
#             df['pval'] = pd.to_numeric(df['pval'], errors='coerce')
#             df = df.dropna()
#             df = df.head(10000)  # Limit to 10,000 points for performance
# 
#             volcano_data = df.to_dict(orient="records")
#         except Exception as e:
#             print("Volcano plot prep failed:", e)
#     print(f"Volcano data points: {len(volcano_data)}")
# 
#     return render_template(
#         "index.html",
#         preview=preview_rows,
#         columns=columns,
#         filename=filename,
#         volcano_data=volcano_data,
#         selected_columns={"id": id_col, "fc": fc_col, "pval": pval_col},
#         logfc_threshold=logfc_threshold,
#         pval_cutoff=pval_cutoff,
#         pval_y=pval_y,
#         case_study_aops=CASE_STUDY_AOPS
#     )


@app.route('/preview', methods=['POST'])
def preview():
    # Prefer uploaded file if available
    file = request.files.get('gene_file')
    demo_filename = request.form.get('demo_file')
    filename = request.form.get('filename')  # used when re-confirming column selection

    # Handle uploaded file
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

    # Handle demo file or re-used file
    elif demo_filename:
        filename = os.path.basename(demo_filename)  # Strips path
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)




    elif filename:
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(filename))

    else:
        return "No dataset provided.", 400

    # Read data
    if not os.path.exists(filepath):
        return f"File not found: {filepath}", 400
    try:
        df = pd.read_csv(filepath, sep=None, engine='python')
    except Exception as e:
        return f"Failed to read dataset: {e}", 400

    preview_rows = df.head(5).to_dict(orient="records")
    columns = df.columns.tolist()

    # Optional column selection (for volcano plot)
    id_col = request.form.get("id_column")
    fc_col = request.form.get("fc_column")
    pval_col = request.form.get("pval_column")

    logfc_threshold = request.form.get("logfc_threshold")
    try:
        logfc_threshold = float(logfc_threshold)
    except (TypeError, ValueError):
        logfc_threshold = None

    pval_cutoff = 0.05
    pval_y = -math.log10(pval_cutoff)
    volcano_data = []

    if id_col and fc_col and pval_col:
        try:
            df_v = df[[id_col, fc_col, pval_col]].dropna()
            df_v.columns = ['ID', 'log2FC', 'pval']
            df_v['ID'] = df_v['ID'].astype(str)
            df_v['log2FC'] = pd.to_numeric(df_v['log2FC'], errors='coerce')
            df_v['pval'] = pd.to_numeric(df_v['pval'], errors='coerce')
            df_v = df_v.dropna().head(10000)
            volcano_data = df_v.to_dict(orient="records")
        except Exception as e:
            print("Volcano plot generation failed:", e)

    return render_template(
        "index.html",
        preview=preview_rows,
        columns=columns,
        filename=filename,
        volcano_data=volcano_data,
        selected_columns={"id": id_col, "fc": fc_col, "pval": pval_col},
        logfc_threshold=logfc_threshold,
        pval_cutoff=pval_cutoff,
        pval_y=pval_y,
        case_study_aops=CASE_STUDY_AOPS
    )




import re

def guess_id_type(gene_series):
    hgnc_pattern = re.compile(r"^[A-Z0-9\-]+$")
    ensembl_pattern = re.compile(r"ENSG\d+")
    entrez_pattern = re.compile(r"^\d+$")

    hgnc, ensembl, entrez = 0, 0, 0
    for gene in gene_series.head(20):
        if ensembl_pattern.match(gene):
            ensembl += 1
        elif entrez_pattern.match(gene):
            entrez += 1
        elif hgnc_pattern.match(gene):
            hgnc += 1

    scores = {"HGNC": hgnc, "Ensembl": ensembl, "Entrez": entrez}
    return max(scores, key=scores.get)






@app.route('/analyze', methods=['POST'])
def analyze():
    # Get info from form
    filename = request.form.get('filename')
    id_col = request.form.get('id_column')
    fc_col = request.form.get('fc_column')
    pval_col = request.form.get('pval_column')
    # Default p-value cutoff
    pval_cutoff = 0.05

    aop_id = request.form.get('aop_selection')
    if not aop_id:
        return "Missing AOP selection", 400

    print(f"[INFO] Selected AOP: {aop_id}")

    # Load all KEs in this AOP
    aop_ke_map = pd.read_csv('data/aop_ke_map.csv')
    ke_list = set(aop_ke_map[aop_ke_map['AOP_ID'] == aop_id]['KE_ID'].dropna())

    # Load all edges for this AOP
    aop_ker_df = pd.read_csv('data/aop_ker_edges.csv')
    edges = aop_ker_df[aop_ker_df['AOP_ID'] == aop_id]

    
    # Restrict reference_sets to only those KEs
    filtered_reference_sets = {ke: genes for ke, genes in reference_sets.items() if ke in ke_list}

    ke_metadata_df = pd.read_csv("data/ke_metadata.csv")
    ke_type_map = dict(zip(ke_metadata_df["KE_ID"], ke_metadata_df["Type"]))
    ke_title_map = dict(zip(ke_metadata_df["KE_ID"], ke_metadata_df["Title"]))


    logfc_threshold = request.form.get("logfc_threshold")
    try:
        logfc_threshold = float(logfc_threshold)
    except (TypeError, ValueError):
        logfc_threshold = 0.0

    # Validate input
    if not (filename and id_col and fc_col and pval_col):
        return "Missing required inputs", 400

    filepath = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(filename))
    if not os.path.exists(filepath):
        return "Uploaded file not found", 400

    # Read file and drop missing values
    df_raw = pd.read_csv(filepath, sep=None, engine='python')

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


    # Combine duplicate genes
    grouped = df.groupby('ID')
    combined_rows = []
    for gene, group in grouped:
        if len(group) == 1:
            combined_rows.append(group.iloc[0])
        else:
            avg_fc = group['log2FC'].astype(float).mean()
            _, combined_p = combine_pvalues(group['pval'].astype(float), method='fisher')
            combined_rows.append({'ID': gene, 'log2FC': avg_fc, 'pval': combined_p})

    df = pd.DataFrame(combined_rows)
    

    df['ID'] = df['ID'].astype(str).str.strip().str.upper()

    # Print first few rows for inspection
    print("\n[INFO] Uploaded dataset preview:")
    print(df.head())

    # Check if THRSP is present
    print("\n[INFO] Is 'THRSP' in dataset?:", 'THRSP' in df['ID'].values)

    # Guess gene ID type for display
    id_type = guess_id_type(df['ID'])

    # Add significance flag
    df['significant'] = False
    if logfc_threshold is not None:
        df['significant'] = (df['log2FC'].abs() >= logfc_threshold) & (df['pval'] <= pval_cutoff)
    print("\n[INFO] Top significant genes:")
    print(df[df['significant']].head())
    n_sig = df['significant'].sum()
    if n_sig == 0:
        print("[WARNING] No significant genes found. Check thresholds.")
    # Build gene-level stats dict: gene → log2FC
    gene_logfc_map = df.set_index("ID")["log2FC"].to_dict()

    # Create a dictionary of gene → significant flag
    gene_significance_map = df.set_index("ID")["significant"].to_dict()

    # Build KE → list of gene data (with log2FC and significance)
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

    # Check if THRSP is significant
    if 'THRSP' in df['ID'].values:
        print("\n[INFO] THRSP entry:")
        print(df[df['ID'] == 'THRSP'])

    # Define gene universe and DEG status
    all_genes = set(df['ID'])
    user_gene_status = dict(zip(df['ID'], df['significant']))

    sig_genes = df['significant'].sum()
    non_sig_genes = df.shape[0] - sig_genes

    print(f"\n[INFO] Total genes: {len(all_genes)} | Significant: {sig_genes} | Non-significant: {non_sig_genes}")

    # Run enrichment
    results = []
    for ke, ref_genes in filtered_reference_sets.items():
        ke_genes = {g.strip().upper() for g in ref_genes} & all_genes

        sig_in_ke = {g for g in ke_genes if user_gene_status.get(g)}
        non_sig_in_ke = ke_genes - sig_in_ke

        a = len(sig_in_ke)                            # sig & in KE
        b = sig_genes - a                             # sig & not in KE
        c = len(non_sig_in_ke)                        # non-sig & in KE
        d = non_sig_genes - c                         # non-sig & not in KE

        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

        print(f"\n[INFO] {ke} - total KE genes in dataset: {len(ke_genes)}")
        print(f"       sig_in_KE: {a}, non_sig_in_KE: {c}")
        print(f"       sig_in_KE genes: {sig_in_ke}")

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
    if not results:
        print("[ERROR] No enrichment results generated — 'results' list is empty!")
            
    df_results = pd.DataFrame(results)
    df_results["FDR"] = multipletests(df_results["p_value"], method="fdr_bh")[1]
    df_results = df_results.sort_values("FDR")

    # Load AOP → KERs
    aop_ker_df = pd.read_csv('data/aop_ker_edges.csv')
    edges = aop_ker_df[aop_ker_df['AOP_ID'] == aop_id]

    # KE titles (optional)
    ke_titles = {ke: ke for ke in filtered_reference_sets.keys()}  # placeholder

    # Nodes = all KEs in results (you can also add others)
    cy_nodes = []
    for ke in ke_list:
        enrichment_row = df_results[df_results['KE'] == ke].squeeze() if ke in df_results['KE'].values else None
        label = ke_title_map.get(ke, ke)
        ke_type = ke_type_map.get(ke, "intermediate")

        logfc = enrichment_row['odds_ratio'] if enrichment_row is not None else 0
        sig = enrichment_row['FDR'] < 0.05 if enrichment_row is not None and 'FDR' in enrichment_row else False

        classes = []
        if sig:
            classes.append("significant")

        cy_nodes.append({
            "data": {
                "id": ke,
                "label": label,
                "logfc": logfc,
                "ke_type": ke_type
            },
            "classes": " ".join(classes)
        })



    cy_edges = []
    for _, row in edges.iterrows():
        cy_edges.append({
            "data": {
                "source": row['Source_KE'],
                "target": row['Target_KE'],
                "id": f"KER:{row['KER_ID']}"
            }
        })


    cy_network = {"nodes": cy_nodes, "edges": cy_edges}


    # Clean up
    # os.remove(filepath)

    return render_template(
        "results.html",
        table=df_results.to_dict(orient="records"),
        id_type=id_type,
        background_size=len(all_genes),
        threshold=logfc_threshold,
        network_json=json.dumps(cy_network),
        ke_gene_json=json.dumps(ke_gene_map),
        ke_type_map=json.dumps(ke_type_map),
        ke_title_map=json.dumps(ke_title_map),
    )



from helpers import load_reference_sets

# Load once at app startup
reference_sets = load_reference_sets(
    ke_wp_path='data/KE-WP.csv',
    wp_gene_path='data/edges_wpid_to_gene.csv',
    node_path='data/node_attributes.csv'
)
print(f"[INFO] Loaded {len(reference_sets)} KE sets from reference files.")
if reference_sets:
    first_ke = list(reference_sets.keys())[0]
    print(f"[INFO] Example KE set: {first_ke} with {len(reference_sets[first_ke])} genes.")

def build_cytoscape_network(aop_id, ker_file='data/aop_ker_edges.csv', metadata=None):
    ker_df = pd.read_csv(ker_file)
    edges = ker_df[ker_df['AOP_ID'] == aop_id]

    nodes = set(edges['Source_KE']) | set(edges['Target_KE'])

    cytoscape_nodes = []
    for ke in nodes:
        node_data = {
            "data": {"id": ke, "label": metadata.get(ke, ke) if metadata else ke}
        }
        cytoscape_nodes.append(node_data)

    cytoscape_edges = []
    for _, row in edges.iterrows():
        cytoscape_edges.append({
            "data": {
                "source": row['Source_KE'],
                "target": row['Target_KE'],
                "id": f"KER:{row['KER_ID']}"
            }
        })

    return {"nodes": cytoscape_nodes, "edges": cytoscape_edges}


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

