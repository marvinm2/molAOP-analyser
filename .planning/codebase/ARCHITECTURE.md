# Architecture

**Analysis Date:** 2026-02-25

## Pattern Overview

**Overall:** Layered service-oriented architecture with Flask web framework

**Key Characteristics:**
- Modular service layer with single-responsibility functions
- Clear separation between request handling, data processing, and business logic
- Reference data caching with TTL expiration for performance
- Database persistence via SQLAlchemy ORM
- Stateless analysis pipeline with JSON serialization for inter-component communication

## Layers

**Presentation Layer:**
- Purpose: Render UI and handle HTTP responses
- Location: `templates/` (Jinja2 templates), `static/` (CSS, assets)
- Contains: HTML templates (`index.html`, `results.html`, `documentation.html`)
- Depends on: Flask context, session data, form data
- Used by: Flask routes in `app.py`

**Route/Controller Layer:**
- Purpose: Handle HTTP requests, validate input, orchestrate services
- Location: `app.py`
- Contains: Flask route handlers (`/`, `/preview`, `/analyze`, `/generate_report`, `/documentation`)
- Depends on: Validation layer, service layer, database layer, caching layer
- Used by: Web clients via HTTP

**Validation Layer:**
- Purpose: Sanitize and validate all user input before processing
- Location: `validation.py` (input validation), `exceptions.py` (error definitions)
- Contains: File validation, column name validation, threshold parsing, AOP selection validation
- Depends on: Config for allowed values
- Used by: Route handlers, service functions

**Service Layer:**
- Purpose: Encapsulate domain logic for specific tasks
- Location: `services/` directory
- Contains modular services:
  - `column_detector.py`: ML-based auto-detection of gene ID, log2FC, p-value columns with confidence scoring
  - `data_service.py`: Data loading, validation, normalization, gene ID type detection, AOP data loading
  - `enrichment_service.py`: Fisher's exact test enrichment analysis with FDR correction
  - `network_service.py`: Cytoscape.js network graph construction from enrichment results
  - `gene_id_validator.py`: Gene symbol validation and normalization against reference
  - `report_service.py`: PDF/HTML report generation with visualizations
- Depends on: Data helpers, reference sets, statistical libraries (scipy, statsmodels)
- Used by: Route handlers

**Data Processing Layer:**
- Purpose: Load, process, and prepare data for analysis
- Location: `helpers.py` (reference set loading), `data_service.py` (normalization)
- Contains:
  - Reference set construction from WikiPathways mappings
  - Gene expression normalization (uppercase, duplicate handling via Fisher's method)
  - Significance threshold filtering
- Depends on: CSV data files in `data/` directory
- Used by: Service layer, enrichment pipeline

**Data Access Layer:**
- Purpose: Persist experiment metadata and analysis results
- Location: `database.py`
- Contains: SQLAlchemy models (`ExperimentRecord`), database manager with CRUD operations
- Depends on: SQLite database (`molAOP_analyser.db`)
- Used by: Route handlers for experiment logging

**Caching Layer:**
- Purpose: Cache frequently accessed reference data to reduce load time
- Location: `cache_manager.py`
- Contains: In-memory cache with TTL expiration, cache key management
- Depends on: None (internal state)
- Used by: Route handlers loading reference sets (1-hour TTL)

**Configuration Layer:**
- Purpose: Centralize all configuration and metadata definitions
- Location: `config.py`
- Contains: File limits, allowed extensions, demo datasets, AOP case studies, experiment metadata dataclass
- Depends on: None (environment variables for secrets)
- Used by: All layers for configuration values

## Data Flow

**File Upload & Preview Flow:**

1. User uploads file or selects demo dataset via `/preview` POST
2. Route handler validates file (extension, size) via `validate_file_upload()` in `validation.py`
3. File copied to `uploads/` directory
4. `pd.read_csv()` loads file with automatic delimiter detection
5. `column_detector.detect_columns()` in `services/column_detector.py` analyzes column headers and data samples using regex patterns and statistical properties
6. User selects/confirms gene ID, log2FC, p-value columns
7. Volcano plot data generated if columns confirmed: `df[[id_col, fc_col, pval_col]]` → JSON → template
8. Preview rendered in `index.html` with column suggestions and sample data

**Analysis & Enrichment Flow:**

1. User submits analysis form via `/analyze` POST with form parameters
2. `validate_form_data()` validates all required fields
3. `load_and_validate_data()` from `services/data_service.py` loads file, extracts columns, expands '///' separated genes
4. `process_gene_expression()` normalizes gene symbols (uppercase), combines duplicates (mean log2FC, Fisher's p-value combination), marks significant genes based on thresholds
5. `load_aop_data()` loads KE list, edges, and metadata for selected AOP
6. `load_cached_reference_sets()` retrieves KE→gene mappings with 1-hour TTL from `cache_manager`
7. `run_enrichment_analysis()` from `services/enrichment_service.py` executes:
   - 2x2 contingency table construction for each KE
   - Fisher's exact test (one-tailed, greater alternative)
   - FDR correction via Benjamini-Hochberg method
   - Results sorted by FDR
8. `build_cytoscape_network()` from `services/network_service.py` constructs network JSON:
   - Nodes: KEs with FDR significance flag, odds ratio for coloring
   - Edges: KE relationships from AOP data, filtered to exist KEs
9. `build_ke_gene_mapping()` creates KE→genes mapping including expression data
10. Experiment metadata stored to database via `db_manager.save_experiment_metadata()`
11. Results rendered in `results.html` with interactive network, enrichment table, volcano plot

**Report Generation Flow:**

1. User clicks "Generate Report" on results page, uploads network PNG snapshot via `/api/upload_network_png`
2. PNG stored in `temp/` with UUID, returns temp_id to client
3. Form submitted to `/generate_report` with analysis parameters + temp_id
4. `report_generator.generate_pdf_report()` or `.generate_html_report()` constructs PDF/HTML with:
   - Metadata (dataset ID, stressor, dosing, owner, description)
   - Analysis parameters (AOP, thresholds, gene ID columns)
   - Enrichment results table (top KEs by FDR)
   - Volcano plot
   - Network visualization (embedded PNG)
5. Temp PNG file cleaned up
6. Report downloaded to client

**State Management:**

- Session storage via Flask `session` object: experiment metadata, timestamps
- Database storage: Complete experiment records with results summary (top 20 enrichment results)
- Cache storage: Reference sets (KE→genes) with 1-hour TTL
- File storage: Uploaded datasets in `uploads/`, temporary network PNGs in `temp/`

## Key Abstractions

**ExperimentMetadata (dataclass in `config.py`):**
- Purpose: Encapsulates experiment information for reports and database
- Examples: Dataset ID, stressor, dosing, owner, description, timestamps
- Pattern: Dataclass with `to_dict()` and `from_dict()` for serialization

**ColumnMatch and ColumnSuggestions (`services/column_detector.py`):**
- Purpose: Represent column detection results with confidence scores
- Pattern: Dataclass-based with ranked suggestions and metadata

**ColumnDetector (`services/column_detector.py`):**
- Purpose: Stateless column detection service with regex patterns and statistical analysis
- Pattern: Class methods analyze column headers, data types, and distributions independently

**ReportData (in `services/report_service.py`):**
- Purpose: Container for all data needed for report generation
- Pattern: Dataclass holding analysis parameters, results, network/volcano visualizations

**ExperimentRecord (SQLAlchemy model in `database.py`):**
- Purpose: Database representation of experiment metadata and results
- Pattern: ORM model with JSON serialization of nested results

**CacheManager (`cache_manager.py`):**
- Purpose: Generic in-memory cache with TTL and thread safety
- Pattern: Singleton instance with lock-protected access, pattern-based clearing

## Entry Points

**Web Application:**
- Location: `app.py` main block at line 795
- Triggers: `python app.py` or Flask development server
- Responsibilities: Starts Flask app on port 5010, initializes database, validates data files, loads initial reference sets

**Route: GET / (Index Page):**
- Location: `app.py:158-165`
- Triggers: Browser navigation to root URL
- Responsibilities: Render main upload/analysis interface with available AOPs

**Route: POST /preview (File Preview):**
- Location: `app.py:215-372`
- Triggers: File upload or demo file selection
- Responsibilities: Validate file, auto-detect columns, generate volcano plot data, return preview

**Route: POST /analyze (Run Analysis):**
- Location: `app.py:415-612`
- Triggers: Form submission with selected columns and AOP
- Responsibilities: Orchestrate complete analysis pipeline, persist results, render interactive results page

**Route: POST /generate_report (Report Generation):**
- Location: `app.py:615-765`
- Triggers: Report download button with format selection
- Responsibilities: Assemble report data, generate PDF/HTML, return as download

**Route: POST /api/upload_network_png (Network PNG Upload):**
- Location: `app.py:176-213`
- Triggers: Client-side network snapshot export
- Responsibilities: Decode base64 PNG, store in temp folder, return UUID for later retrieval

## Error Handling

**Strategy:** Layered error handling with user-friendly messages and technical logging

**Patterns:**

- **Custom Exception Hierarchy** (`exceptions.py`): Base `AOPAnalysisError` with specific subclasses
  - `DataValidationError`: Input validation failures
  - `FileProcessingError`: File reading/parsing failures
  - `EnrichmentAnalysisError`: Analysis computation failures
  - `AOPDataError`: Missing/invalid AOP reference data
  - `NetworkBuildError`: Network visualization failures

- **Error Formatting** (`exceptions.py:format_error_response()`): Returns user-friendly message with optional technical details

- **Route-Level Handling** (`app.py`): Try/except blocks catch specific exceptions, log errors, return HTTP responses with status codes

- **Validation-First** (`validation.py`): All user input validated before processing with detailed error messages

- **Graceful Degradation**: Database save failures don't block analysis (warning logged), network PNG upload failures handled independently

## Cross-Cutting Concerns

**Logging:** Configured in `app.py:31-36` with INFO level to `flask.log`
- Modules import logger via `logger = logging.getLogger(__name__)`
- Service modules log key steps: data loading, analysis progress, cache operations
- Error logging includes stack traces (exc_info=True)

**Validation:**
- Input validation in `validation.py` before route processing
- File upload validation: extension whitelist, size limits, file existence checks
- Column validation: length limits, dangerous character filtering
- Form data validation: required fields, threshold bounds, AOP selection against enabled list

**Authentication:** CSRF protection via `flask_wtf.csrf.CSRFProtect`
- Injected into all templates via `@app.context_processor` (line 49-52)
- Required for form submissions via `<input type="hidden" name="csrf_token">`

**File Security:**
- Filename sanitization via `werkzeug.utils.secure_filename()` before storage
- Path traversal prevention: reject `..` and absolute paths in demo file selection
- Upload folder restricted to `uploads/`, temp folder restricted to `temp/`
- Cleanup via `cleanup_file()` utility and `cleanup_old_temp_files()` for stale uploads

**Data Privacy:**
- Session storage restricted to experiment metadata (no raw gene data)
- Database stores only metadata and results summary (top 20 enrichment results), not full gene tables
- Uploaded files deleted after analysis completion

---

*Architecture analysis: 2026-02-25*
