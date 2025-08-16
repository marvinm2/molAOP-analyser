from flask import Flask, render_template, request, abort, session, make_response, send_file
import pandas as pd
import os
import json
import math
import logging
from typing import Optional
from scipy.stats import fisher_exact, combine_pvalues
from statsmodels.stats.multitest import multipletests
from werkzeug.utils import secure_filename
from werkzeug.exceptions import RequestEntityTooLarge
from flask_wtf.csrf import CSRFProtect, generate_csrf

from config import Config, ExperimentMetadata
from validation import validate_form_data, validate_file_upload, log_validation_error
from helpers import load_reference_sets
from cache_manager import cache, cached_data_loader
from exceptions import AOPAnalysisError, format_error_response
from utils import cleanup_file, validate_file_path
from services.data_service import load_and_validate_data, process_gene_expression, guess_id_type, load_aop_data
from services.enrichment_service import run_enrichment_analysis, build_ke_gene_mapping
from services.network_service import build_cytoscape_network
from services.column_detector import column_detector
from services.gene_id_validator import gene_id_validator
from database import db_manager, init_database
from services.report_service import report_generator, ReportData, get_software_versions

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = Flask(__name__)
app.config.from_object(Config)
app.config['UPLOAD_FOLDER'] = Config.UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = Config.MAX_CONTENT_LENGTH  # Use higher limit for form data
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# Initialize CSRF protection
csrf = CSRFProtect(app)

@app.context_processor
def inject_csrf_token():
    """Make CSRF token available in all templates."""
    return dict(csrf_token=generate_csrf())

def allowed_file(filename):
    """Check if uploaded file has an allowed extension.
    
    Args:
        filename: Name of the uploaded file
        
    Returns:
        bool: True if file extension is allowed, False otherwise
    """
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in Config.ALLOWED_EXTENSIONS

@app.errorhandler(RequestEntityTooLarge)
def handle_file_too_large(e):
    """Handle file upload size limit exceeded error.
    
    Args:
        e: RequestEntityTooLarge exception
        
    Returns:
        tuple: Error message and HTTP status code
    """
    return f"File too large. Maximum upload size is {Config.MAX_FILE_SIZE // (1024*1024)} MB.", 413

@app.errorhandler(404)
def handle_not_found(e):
    """Handle 404 Not Found errors.
    
    Args:
        e: Not found exception
        
    Returns:
        tuple: Error message and HTTP status code
    """
    return "Page not found. Please check the URL and try again.", 404

@app.errorhandler(500)
def handle_internal_error(e):
    """Handle internal server errors with logging.
    
    Args:
        e: Internal server error exception
        
    Returns:
        tuple: Error message and HTTP status code
    """
    logger.error(f"Internal server error: {e}", exc_info=True)
    return "An internal server error occurred. Please try again later.", 500


def extract_metadata_from_request() -> Optional[ExperimentMetadata]:
    """Extract experiment metadata from request form data.
    
    Returns:
        ExperimentMetadata object if metadata fields are present, None otherwise
    """
    metadata_fields = ['dataset_id', 'stressor', 'dosing', 'owner', 'description']
    metadata_values = {}
    
    # Check if any metadata fields are present and non-empty
    has_metadata = False
    for field in metadata_fields:
        value = request.form.get(field, '').strip()
        if value:
            has_metadata = True
        metadata_values[field] = value
    
    if not has_metadata:
        return None
    
    try:
        return ExperimentMetadata(
            dataset_id=metadata_values.get('dataset_id', ''),
            stressor=metadata_values.get('stressor', ''),
            dosing=metadata_values.get('dosing', ''),
            owner=metadata_values.get('owner', ''),
            description=metadata_values.get('description', '')
        )
    except Exception as e:
        logger.warning(f"Failed to create metadata object: {e}")
        return None


@app.route('/')
def index():
    """Render the main application page with available AOP case studies.
    
    Returns:
        str: Rendered HTML template for the main page
    """
    return render_template('index.html', case_study_aops=Config.CASE_STUDY_AOPS)

@app.route('/preview', methods=['POST'])
def preview():
    """Handle file upload/selection and provide dataset preview with column auto-detection.
    
    Supports three input modes:
    1. New file upload via 'gene_file'
    2. Demo file selection via 'demo_file' 
    3. Re-processing existing file via 'filename'
    
    Performs automatic column detection to suggest gene ID, log2FC, and p-value columns.
    Generates volcano plot data when columns are selected.
    
    Returns:
        str: Rendered HTML template with dataset preview and column suggestions
    """
    # Extract and store experiment metadata if provided
    metadata = extract_metadata_from_request()
    if metadata:
        session['experiment_metadata'] = metadata.to_dict()
        logger.info(f"Stored experiment metadata for dataset: {metadata.dataset_id}")
    
    # Prefer uploaded file if available
    file = request.files.get('gene_file')
    demo_filename = request.form.get('demo_file')
    filename = request.form.get('filename')  # used when re-confirming column selection

    # Handle uploaded file
    if file and file.filename:
        # Validate uploaded file for security and format compliance
        is_valid, error_msg = validate_file_upload(file)
        if not is_valid:
            log_validation_error("file_upload", error_msg, {"filename": file.filename})
            return error_msg, 400
            
        # Secure the filename and save to upload directory
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

    # Handle demo file selection (pre-loaded datasets)
    elif demo_filename:
        filename = os.path.basename(demo_filename)  # Security: strip any path components
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        
        # Copy demo file from data directory to uploads directory
        source_path = os.path.join('data', filename)
        if not os.path.exists(source_path):
            return f"Demo file {filename} not found.", 404
        
        # Copy the file to uploads directory
        import shutil
        shutil.copy2(source_path, filepath)
        logger.info(f"Copied demo file {filename} from {source_path} to {filepath}")

    # Handle re-processing of existing uploaded file
    elif filename:
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(filename))

    else:
        return "No dataset provided.", 400

    # Load and parse dataset with automatic delimiter detection
    if not os.path.exists(filepath):
        return f"File not found: {filepath}", 400
    try:
        # Use pandas' automatic separator detection for CSV/TSV files
        df = pd.read_csv(filepath, sep=None, engine='python')
    except Exception as e:
        return f"Failed to read dataset: {e}", 400

    # Prepare dataset preview (first 5 rows) and extract column names
    preview_rows = df.head(5).to_dict(orient="records")
    columns = df.columns.tolist()
    
    # Auto-detect column types if not already specified by user
    column_suggestions = None
    id_col = request.form.get("id_column")    # User-selected gene ID column
    fc_col = request.form.get("fc_column")    # User-selected log2FC column
    pval_col = request.form.get("pval_column")  # User-selected p-value column
    
    # Run intelligent column detection if user hasn't selected all required columns
    if not (id_col and fc_col and pval_col):
        try:
            # Use ML-based column detection to identify gene IDs, log2FC, and p-values
            column_suggestions = column_detector.detect_columns(df)
            logger.info("Column auto-detection completed")
            
            # Auto-select high-confidence suggestions when user hasn't chosen
            if not id_col and column_suggestions.best_gene_id:
                id_col = column_suggestions.best_gene_id.column_name
                logger.info(f"Auto-selected gene ID column: {id_col}")
            
            if not fc_col and column_suggestions.best_log2fc:
                fc_col = column_suggestions.best_log2fc.column_name
                logger.info(f"Auto-selected log2FC column: {fc_col}")
                
            if not pval_col and column_suggestions.best_pvalue:
                pval_col = column_suggestions.best_pvalue.column_name
                logger.info(f"Auto-selected p-value column: {pval_col}")
                
        except Exception as e:
            logger.warning(f"Column auto-detection failed: {e}")
            column_suggestions = None

    # Parse log2FC threshold for significance filtering (optional)
    logfc_threshold = request.form.get("logfc_threshold")
    try:
        logfc_threshold = float(logfc_threshold)
    except (TypeError, ValueError):
        logfc_threshold = None

    # Prepare volcano plot configuration
    pval_cutoff = Config.PVAL_CUTOFF  # Statistical significance threshold (typically 0.05)
    pval_y = -math.log10(pval_cutoff)  # Y-axis position for significance line
    volcano_data = []  # Will store plot data points

    # Generate volcano plot data when all required columns are available
    if id_col and fc_col and pval_col:
        try:
            # Extract and standardize the three essential columns
            df_v = df[[id_col, fc_col, pval_col]].dropna()
            df_v.columns = ['ID', 'log2FC', 'pval']
            
            # Ensure proper data types for plotting
            df_v['ID'] = df_v['ID'].astype(str)  # Gene identifiers as strings
            df_v['log2FC'] = pd.to_numeric(df_v['log2FC'], errors='coerce')  # Numeric fold changes
            df_v['pval'] = pd.to_numeric(df_v['pval'], errors='coerce')  # Numeric p-values
            
            # Remove any remaining invalid data and limit for performance
            df_v = df_v.dropna().head(Config.MAX_GENES_DISPLAY)
            volcano_data = df_v.to_dict(orient="records")
        except Exception as e:
            logger.warning(f"Volcano plot generation failed: {e}")

    return render_template(
        "index.html",
        preview=preview_rows,
        columns=columns,
        filename=filename,
        volcano_data=volcano_data,
        selected_columns={"id": id_col, "fc": fc_col, "pval": pval_col},
        column_suggestions=column_suggestions,
        logfc_threshold=logfc_threshold,
        pval_cutoff=pval_cutoff,
        pval_y=pval_y,
        case_study_aops=Config.CASE_STUDY_AOPS
    )




import re

def guess_id_type(gene_series):
    """Determine the most likely gene identifier type in a pandas Series.
    
    Analyzes the first 20 gene identifiers to classify them as HGNC symbols,
    Ensembl IDs, or Entrez IDs based on pattern matching.
    
    Args:
        gene_series: pandas Series containing gene identifiers
        
    Returns:
        str: Most likely identifier type ('HGNC', 'Ensembl', or 'Entrez')
    """
    # Define regex patterns for common gene identifier formats
    hgnc_pattern = re.compile(r"^[A-Z0-9\-]+$")  # Uppercase letters, numbers, hyphens
    ensembl_pattern = re.compile(r"ENSG\d+")      # Ensembl gene format
    entrez_pattern = re.compile(r"^\d+$")         # Pure numeric IDs

    # Score each pattern against sample of identifiers
    hgnc, ensembl, entrez = 0, 0, 0
    for gene in gene_series.head(20):  # Sample first 20 for efficiency
        if ensembl_pattern.match(gene):
            ensembl += 1
        elif entrez_pattern.match(gene):
            entrez += 1
        elif hgnc_pattern.match(gene):
            hgnc += 1

    # Return the type with the highest match count
    scores = {"HGNC": hgnc, "Ensembl": ensembl, "Entrez": entrez}
    return max(scores, key=scores.get)






@app.route('/analyze', methods=['POST'])
def analyze():
    """Refactored analyze function using service modules."""
    try:
        # Update experiment metadata if provided
        metadata = extract_metadata_from_request()
        if metadata:
            session['experiment_metadata'] = metadata.to_dict()
            logger.info(f"Updated experiment metadata for analysis: {metadata.dataset_id}")
        
        # Get form parameters
        form_data = {
            'filename': request.form.get('filename'),
            'id_column': request.form.get('id_column'),
            'fc_column': request.form.get('fc_column'),
            'pval_column': request.form.get('pval_column'),
            'aop_selection': request.form.get('aop_selection'),
            'logfc_threshold': request.form.get('logfc_threshold', '0.0')
        }
        
        # Validate form data
        is_valid, errors = validate_form_data(form_data)
        if not is_valid:
            error_msg = "Validation errors: " + "; ".join(errors)
            log_validation_error("form_validation", error_msg, form_data)
            return error_msg, 400
        
        # Extract validated parameters
        filename = form_data['filename']
        id_col = form_data['id_column']
        fc_col = form_data['fc_column']
        pval_col = form_data['pval_column']
        aop_id = form_data['aop_selection']
            
        # Parse threshold
        try:
            logfc_threshold = float(request.form.get("logfc_threshold", 0.0))
        except (TypeError, ValueError):
            logfc_threshold = 0.0
            
        logger.info(f"Starting analysis for AOP: {aop_id}")
        
        # Validate file path
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(filename))
        if not validate_file_path(filepath):
            return "Invalid or missing file", 400
        
        # Load and process data using services
        df_raw = load_and_validate_data(filepath, id_col, fc_col, pval_col)
        df_processed, stats = process_gene_expression(df_raw, logfc_threshold)
        
        # Generate volcano plot data for PDF reports
        volcano_data = []
        try:
            # Use the processed data which already has standardized column names
            df_volcano = df_processed[['ID', 'log2FC', 'pval']].copy()
            
            # Ensure proper data types for plotting
            df_volcano['ID'] = df_volcano['ID'].astype(str)
            df_volcano['log2FC'] = pd.to_numeric(df_volcano['log2FC'], errors='coerce')
            df_volcano['pval'] = pd.to_numeric(df_volcano['pval'], errors='coerce')
            
            # Remove invalid data and limit for performance (first 2000 genes for volcano plot to reduce form size)
            df_volcano = df_volcano.dropna().head(2000)
            volcano_data = df_volcano.to_dict(orient="records")
            logger.info(f"Generated volcano plot data with {len(volcano_data)} genes")
            if volcano_data:
                logger.info(f"Volcano data sample: {volcano_data[0]}")
        except Exception as e:
            logger.warning(f"Volcano plot data generation failed: {e}")
            volcano_data = []
        
        # Load AOP data
        ke_list, edges, ke_type_map, ke_title_map = load_aop_data(aop_id)
        
        # Get cached reference sets and run enrichment analysis
        current_reference_sets = load_cached_reference_sets()
        enrichment_results = run_enrichment_analysis(
            df_processed, current_reference_sets, ke_list, ke_title_map
        )
        
        # Build network visualization data
        cy_network = build_cytoscape_network(
            ke_list, edges, enrichment_results, ke_title_map, ke_type_map
        )
        
        # Build gene mapping for interactive features
        gene_logfc_map = df_processed.set_index("ID")["log2FC"].to_dict()
        gene_significance_map = df_processed.set_index("ID")["significant"].to_dict()
        ke_gene_map = build_ke_gene_mapping(
            current_reference_sets, ke_list, gene_logfc_map, gene_significance_map
        )
        
        # Guess gene ID type for display
        id_type = guess_id_type(df_processed['ID'])
        
        # Clean up uploaded file
        cleanup_file(filepath)
        
        logger.info("Analysis completed successfully")
        
        # Get stored metadata for the results page
        stored_metadata = session.get('experiment_metadata', {})
        
        # Add analysis parameters to metadata for template
        stored_metadata['filename'] = filename
        stored_metadata['aop_id'] = aop_id
        # Find AOP label
        aop_label = aop_id
        for k, v in Config.CASE_STUDY_AOPS.items():
            if v.get('id') == aop_id:
                aop_label = v.get('label', aop_id)
                break
        stored_metadata['aop_label'] = aop_label
        stored_metadata['logfc_threshold'] = logfc_threshold
        stored_metadata['pval_cutoff'] = Config.PVAL_CUTOFF
        stored_metadata['id_column'] = id_col
        stored_metadata['fc_column'] = fc_col
        stored_metadata['pval_column'] = pval_col
        stored_metadata['significant_genes'] = len(df_processed[df_processed['significant'] == True])
        
        # Save experiment metadata and results to database
        try:
            if stored_metadata:
                analysis_params = {
                    'aop_id': aop_id,
                    'logfc_threshold': logfc_threshold,
                    'pval_cutoff': Config.PVAL_CUTOFF,
                    'id_column': id_col,
                    'fc_column': fc_col,
                    'pval_column': pval_col,
                }
                
                results_summary = {
                    'enrichment_table': enrichment_results.head(20).to_dict('records'),  # Save top 20 results
                    'gene_count': stats['total_genes'],
                    'significant_genes': len(df_processed[df_processed['significant'] == True]),
                }
                
                experiment_id = db_manager.save_experiment_metadata(
                    metadata=stored_metadata,
                    analysis_params=analysis_params,
                    results=results_summary
                )
                
                if experiment_id:
                    logger.info(f"Saved analysis results to database with ID: {experiment_id}")
                    stored_metadata['experiment_id'] = experiment_id
                    
        except Exception as e:
            # Don't fail the analysis if database save fails
            logger.warning(f"Failed to save experiment to database: {e}")
        
        # Convert pandas DataFrame to native Python types for JSON serialization
        enrichment_table = enrichment_results.to_dict(orient="records")
        
        # Convert any numpy/pandas types to native Python types
        for row in enrichment_table:
            for key, value in row.items():
                if pd.isna(value):
                    row[key] = None
                elif hasattr(value, 'item'):  # numpy scalar
                    row[key] = value.item()
                elif isinstance(value, (pd.Series, pd.DataFrame)):
                    row[key] = str(value)  # Convert complex pandas objects to strings
        
        # Pre-serialize JSON for template
        table_json_str = json.dumps(enrichment_table)
        logger.info(f"Template JSON: length={len(table_json_str)}, first 200 chars: {table_json_str[:200]}")
        
        return render_template(
            "results.html",
            table=enrichment_table,
            table_json=table_json_str,  # Pre-serialize for template
            volcano_data=volcano_data,
            volcano_json=json.dumps(volcano_data),  # Pre-serialize for template
            id_type=id_type,
            background_size=stats['total_genes'],
            threshold=logfc_threshold,
            network_json=json.dumps(cy_network),
            ke_gene_json=json.dumps(ke_gene_map),
            ke_type_map=json.dumps(ke_type_map),
            ke_title_map=json.dumps(ke_title_map),
            metadata=stored_metadata,
        )
        
    except AOPAnalysisError as e:
        logger.error(f"AOP analysis error: {e}")
        error_response = format_error_response(e, include_details=app.debug)
        return error_response["message"], 400
        
    except ValueError as e:
        logger.error(f"Validation error in analysis: {e}")
        return f"Analysis error: {str(e)}", 400
        
    except Exception as e:
        logger.error(f"Unexpected error in analysis: {e}", exc_info=True)
        return "An unexpected error occurred during analysis. Please try again or contact support.", 500


@app.route('/generate_report', methods=['POST'])
def generate_report():
    """Generate and download a comprehensive analysis report in PDF or HTML format."""
    try:
        # Get report format preference
        report_format = request.form.get('format', 'pdf').lower()
        
        # Get experiment metadata from session or form data
        metadata = session.get('experiment_metadata', {})
        logger.info(f"Session metadata for report generation: {list(metadata.keys()) if metadata else 'None'}")
        
        # If no session metadata, create from form data
        if not metadata:
            # Check if we have sufficient form data to generate a report
            form_filename = request.form.get('filename', '')
            form_aop_id = request.form.get('aop_id', '')
            logger.info(f"Form data check - filename: '{form_filename}', aop_id: '{form_aop_id}'")
            
            if not form_filename or not form_aop_id:
                logger.error(f"Missing required form data. Available form keys: {list(request.form.keys())}")
                return "No experiment data found. Please run an analysis first.", 400
            
            # Create minimal metadata from form data
            metadata = {
                'filename': form_filename,
                'dataset_id': request.form.get('dataset_id', 'Unknown Dataset'),
                'stressor': request.form.get('stressor', 'Unknown Stressor'),
                'dosing': request.form.get('dosing', 'Unknown Dosing'),
                'owner': request.form.get('owner', 'Unknown Owner'),
                'description': request.form.get('description', '')
            }
            logger.info(f"Created metadata from form data: {metadata}")
        
        # Extract report data from form/session
        try:
            enrichment_results_str = request.form.get('enrichment_results', '[]')
            logger.info(f"Enrichment results string (first 200 chars): {enrichment_results_str[:200]}")
            logger.info(f"Enrichment results string length: {len(enrichment_results_str)}")
            
            # Try to decode JSON, handling potential HTML entity issues
            if enrichment_results_str and enrichment_results_str != '[]':
                # Handle HTML entities that might have been escaped
                import html
                enrichment_results_str = html.unescape(enrichment_results_str)
                logger.info(f"After HTML unescape: {enrichment_results_str[:200]}")
                enrichment_results = json.loads(enrichment_results_str)
                logger.info(f"Successfully parsed {len(enrichment_results)} enrichment results")
                if enrichment_results:
                    logger.info(f"First result keys: {list(enrichment_results[0].keys())}")
            else:
                logger.warning("No enrichment results found in form data")
                enrichment_results = []
                
        except json.JSONDecodeError as e:
            logger.error(f"JSON decode error for enrichment_results: {e}")
            logger.error(f"Problematic JSON: {enrichment_results_str[:500]}")
            # Fallback: try to get enrichment results from session or recreate them
            enrichment_results = []
        
        # Get network PNG data and log it
        network_png_data = request.form.get('network_png', '')
        logger.info(f"Network PNG form data length: {len(network_png_data)}")
        if network_png_data:
            logger.info(f"Network PNG starts with: {network_png_data[:50]}...")
        else:
            logger.warning("No network PNG data received from form")
        
        report_data = ReportData(
            metadata=metadata,
            filename=request.form.get('filename', metadata.get('filename', 'unknown')),
            gene_count=int(request.form.get('gene_count', 0)),
            significant_genes=int(request.form.get('significant_genes', 0)),
            aop_id=request.form.get('aop_id', ''),
            aop_label=request.form.get('aop_label', ''),
            logfc_threshold=float(request.form.get('logfc_threshold', 0.0)),
            pval_cutoff=float(request.form.get('pval_cutoff', 0.05)),
            id_column=request.form.get('id_column', ''),
            fc_column=request.form.get('fc_column', ''),
            pval_column=request.form.get('pval_column', ''),
            id_type=request.form.get('id_type', ''),
            enrichment_results=enrichment_results,
            volcano_data=json.loads(request.form.get('volcano_data', '[]')) if request.form.get('volcano_data') else None,
            network_data=json.loads(request.form.get('network_data', '{}')) if request.form.get('network_data') else None,
            network_png=network_png_data,
            software_versions=get_software_versions()
        )
        
        # Generate report based on format
        if report_format == 'pdf':
            try:
                pdf_content = report_generator.generate_pdf_report(report_data)
                
                # Create response with PDF
                response = make_response(pdf_content)
                response.headers['Content-Type'] = 'application/pdf'
                response.headers['Content-Disposition'] = f'attachment; filename="AOP_Analysis_Report_{metadata.get("dataset_id", "report")}.pdf"'
                
                logger.info(f"PDF report generated for dataset: {metadata.get('dataset_id')}")
                return response
                
            except RuntimeError as e:
                logger.warning(f"PDF generation failed: {e}")
                return f"PDF generation not available: {e}", 400
                
        elif report_format == 'html':
            html_content = report_generator.generate_html_report(report_data)
            
            # Create response with HTML
            response = make_response(html_content)
            response.headers['Content-Type'] = 'text/html'
            response.headers['Content-Disposition'] = f'attachment; filename="AOP_Analysis_Report_{metadata.get("dataset_id", "report")}.html"'
            
            logger.info(f"HTML report generated for dataset: {metadata.get('dataset_id')}")
            return response
            
        else:
            return f"Unsupported report format: {report_format}", 400
            
    except ValueError as e:
        logger.error(f"Invalid report parameters: {e}")
        return f"Invalid report parameters: {e}", 400
        
    except Exception as e:
        logger.error(f"Report generation failed: {e}", exc_info=True)
        return "An error occurred while generating the report. Please try again.", 500


def load_cached_reference_sets():
    """Load AOP reference gene sets with intelligent caching.
    
    Loads the mapping between Key Events (KEs) and WikiPathways gene sets,
    with 1-hour TTL caching to improve performance on repeated requests.
    
    Returns:
        dict: Mapping of KE IDs to sets of associated genes
    """
    return cached_data_loader(
        cache_key="reference_sets",
        loader_func=lambda: load_reference_sets(
            ke_wp_path='data/KE-WP.csv',           # KE to WikiPathway mappings
            wp_gene_path='data/edges_wpid_to_gene.csv',  # WikiPathway to gene mappings
            node_path='data/node_attributes.csv'          # Gene annotation data
        ),
        ttl=3600  # Cache for 1 hour to balance freshness and performance
    )

# Load initial reference sets
reference_sets = load_cached_reference_sets()
logger.info(f"Loaded {len(reference_sets)} KE sets from reference files.")
if reference_sets:
    first_ke = list(reference_sets.keys())[0]
    logger.info(f"Example KE set: {first_ke} with {len(reference_sets[first_ke])} genes.")


if __name__ == '__main__':
    try:
        # Validate required data files exist
        Config.validate_data_files()
        logger.info("All required data files validated successfully")
        
        # Initialize database
        if init_database():
            logger.info("Database initialized successfully")
        else:
            logger.warning("Database initialization failed - continuing without persistence")
        
        # Start the application
        app.run(debug=True, host='0.0.0.0', port=5000)
        
    except FileNotFoundError as e:
        logger.error(f"Startup failed: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected startup error: {e}")
        raise

