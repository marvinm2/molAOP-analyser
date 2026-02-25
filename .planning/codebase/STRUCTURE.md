# Codebase Structure

**Analysis Date:** 2026-02-25

## Directory Layout

```
molAOP-analyser/
├── app.py                          # Flask application and route handlers
├── config.py                       # Configuration and metadata definitions
├── database.py                     # SQLAlchemy ORM models and database manager
├── helpers.py                      # Reference set loading and data merging
├── validation.py                   # Input validation utilities
├── exceptions.py                   # Custom exception classes
├── utils.py                        # General utility functions
├── cache_manager.py                # In-memory cache with TTL
├── services/                       # Service layer with domain logic
│   ├── __init__.py
│   ├── column_detector.py          # Auto-detection of gene ID/log2FC/p-value columns
│   ├── data_service.py             # Data loading, normalization, AOP data handling
│   ├── enrichment_service.py       # Fisher's exact test enrichment analysis
│   ├── network_service.py          # Cytoscape.js network construction
│   ├── gene_id_validator.py        # Gene symbol validation
│   └── report_service.py           # PDF/HTML report generation
├── data/                           # Reference datasets and demo files
│   ├── aop_ke_map.csv              # AOP ID → KE ID mappings
│   ├── aop_ker_edges.csv           # KE relationship edges (KER data)
│   ├── KE-WP.csv                   # KE ID → WikiPathways ID mappings
│   ├── edges_wpid_to_gene.csv      # WikiPathways ID → gene ID edges
│   ├── node_attributes.csv         # Gene ID → gene symbol mappings
│   ├── ke_metadata.csv             # KE ID → title/type metadata
│   ├── wikipathways_hsa_20240410.xgmml  # WikiPathways network graph
│   ├── GSE90122_SR12813.tsv        # Demo: PXR agonist 1 dataset
│   ├── GSE90122_TO90137.tsv        # Demo: PXR agonist 2 dataset
│   ├── Cisplatin_Kidney/           # Demo: Cisplatin toxicity (42 files)
│   │   ├── CSP_4hr_0.1uM.csv
│   │   ├── CSP_4hr_0.5uM.csv
│   │   └── ... (9 timepoints × 5 doses = 45 files total)
│   └── parse_xgmml_to_dataframe.py # XGMML parser utility
├── templates/                      # Jinja2 HTML templates
│   ├── index.html                  # Main upload and analysis interface
│   ├── results.html                # Interactive results with network visualization
│   └── documentation.html          # Documentation/help page
├── static/                         # CSS and static assets
│   ├── css/
│   │   └── style.css               # Main application styles with house colors
│   └── img/
│       └── logo.png
├── tests/                          # Pytest test suite
│   ├── __init__.py
│   ├── conftest.py                 # Pytest fixtures (flask_client, mocks)
│   ├── test_flask_routes.py        # Integration tests for Flask routes
│   ├── test_column_detector.py     # Unit tests for column detection
│   ├── test_database.py            # Unit tests for database operations
│   ├── test_report_service.py      # Unit tests for report generation
│   └── ...
├── uploads/                        # User-uploaded dataset files (transient)
├── temp/                           # Temporary network PNG snapshots (transient)
├── molAOP_analyser.db              # SQLite database (created at runtime)
├── pytest.ini                      # Pytest configuration
├── requirements.txt                # Python dependencies
├── Dockerfile                      # Docker image definition
├── docker-compose.yml              # Docker Compose orchestration
├── CLAUDE.md                       # Claude Code development guidelines
└── README.md                       # Project documentation
```

## Directory Purposes

**Root Level:**
- Purpose: Application entry point, configuration, and global utilities
- Contains: Flask app, database setup, validation, error handling
- Key files: `app.py` (main application), `config.py` (settings), `database.py` (persistence)

**services/**
- Purpose: Modular business logic encapsulation
- Contains: Stateless service classes and functions for specific domains
- Key files: Enrichment analysis, data processing, network visualization, report generation
- Pattern: Each service file focuses on one domain concern

**data/**
- Purpose: Reference datasets for enrichment analysis and demo data
- Contains: CSV mappings (KE→genes), WikiPathways data, demo gene expression datasets
- Key files: `aop_ke_map.csv` (core AOP reference), `edges_wpid_to_gene.csv` (KE→genes mapping)
- Committed: Yes (required for application function, 42 cisplatin demo files)

**templates/**
- Purpose: Jinja2 HTML templates for web interface
- Contains: Form interface (upload, analysis), results visualization, documentation
- Key files: `index.html` (80% of user interaction), `results.html` (Cytoscape network), `documentation.html`

**static/**
- Purpose: CSS and image assets
- Contains: Application styles (with official house colors), logo
- Key files: `style.css` (color definitions at top of file)

**tests/**
- Purpose: Test suite for quality assurance
- Contains: Unit tests, integration tests, Flask route tests
- Pattern: Parallel structure to main code (test_X.py matches X.py)
- Config: `conftest.py` provides Flask client, database fixtures

**uploads/**
- Purpose: Temporary storage of uploaded user files during analysis
- Committed: No
- Lifecycle: Files deleted after analysis completion

**temp/**
- Purpose: Temporary storage of network PNG snapshots from client
- Committed: No
- Lifecycle: Files deleted after report generation or expire after 1 hour

## Key File Locations

**Entry Points:**
- `app.py`: Main Flask application (starts server, initializes DB, loads reference sets)
- `app.py:795-816`: Application startup block with data validation and database initialization

**Configuration:**
- `config.py`: All configuration constants (file limits, demo datasets, AOPs, secret key)
- `config.py:32-90`: DEMO_DATASETS and CASE_STUDY_AOPS definitions

**Core Logic:**
- `app.py:415-612`: `/analyze` route orchestrating analysis pipeline
- `services/data_service.py:70-122`: Data processing with duplicate handling
- `services/enrichment_service.py:12-137`: Fisher's exact test implementation
- `services/network_service.py:10-92`: Cytoscape network JSON construction

**Data Loading:**
- `helpers.py:3-67`: Reference set construction from WikiPathways mappings
- `services/data_service.py:154-184`: AOP-specific data loading
- `cache_manager.py:65-84`: Cached data loading with TTL

**Testing:**
- `tests/conftest.py`: Pytest fixtures (Flask test client, mocks)
- `tests/test_flask_routes.py`: Route integration tests
- `tests/test_column_detector.py`: Column detection unit tests

## Naming Conventions

**Files:**
- Services: `{domain}_service.py` (e.g., `enrichment_service.py`)
- Tests: `test_{module}.py` (e.g., `test_data_service.py`)
- Data: Snake_case CSV names (e.g., `aop_ke_map.csv`)
- Templates: Kebab-case or descriptive lowercase (e.g., `index.html`, `results.html`)

**Directories:**
- Service modules: `services/` (lowercase)
- Data files: `data/` (lowercase)
- Templates: `templates/` (lowercase)
- Tests: `tests/` (lowercase)
- Uploaded files: `uploads/` (lowercase)
- Temporary files: `temp/` (lowercase)

**Classes:**
- Service classes: PascalCase with domain suffix (e.g., `ColumnDetector`, `DatabaseManager`)
- Dataclasses: PascalCase (e.g., `ExperimentMetadata`, `ColumnMatch`, `ReportData`)
- Exception classes: PascalCase ending with "Error" (e.g., `AOPAnalysisError`, `DataValidationError`)

**Functions:**
- Utilities: snake_case (e.g., `load_reference_sets()`, `validate_file_upload()`)
- Route handlers: snake_case (e.g., `index()`, `preview()`, `analyze()`)
- Service functions: snake_case (e.g., `run_enrichment_analysis()`, `build_cytoscape_network()`)

**Variables:**
- Constants: UPPER_SNAKE_CASE (e.g., `MAX_FILE_SIZE`, `PVAL_CUTOFF`)
- Local variables: snake_case (e.g., `enrichment_results`, `gene_logfc_map`)
- DataFrame columns: Standardized (e.g., `ID`, `log2FC`, `pval`, `significant`)

## Where to Add New Code

**New Feature (End-to-End):**
- Primary code: `services/{domain}_service.py` (create if new domain)
- Route handler: Add route to `app.py` with validation, service calls
- Tests: `tests/test_{domain}_service.py` for unit tests, add to `tests/test_flask_routes.py` for integration
- Template: Add HTML template in `templates/` if UI needed
- Config: Add settings to `config.py` if needed

**New Service/Module:**
- Implementation: `services/new_service.py`
- Exports: Stateless functions or single-instance classes
- Logging: Import logger and log key operations
- Testing: Corresponding `tests/test_new_service.py` with mocked dependencies
- Dependencies: Import from config, utils, exceptions; avoid importing from app.py (circular imports)

**New Reference Data:**
- CSV files: Place in `data/` directory with descriptive name
- Loading: Add loading logic to `helpers.py` or `services/data_service.py`
- Mapping: Update `load_reference_sets()` or create new loader function
- Validation: Add file existence check to `Config.validate_data_files()` if required

**Utilities:**
- Shared helpers: `utils.py` for general utility functions
- Validation rules: `validation.py` for input validation
- Error handling: Add exception class to `exceptions.py` if new error type

**Database Operations:**
- Models: Add SQLAlchemy model class to `database.py:Base` inheritance
- Manager methods: Add CRUD method to `DatabaseManager` class in `database.py`
- Usage: Call via `db_manager` singleton from route handlers

**Templates:**
- Location: `templates/{name}.html`
- Rendering: Call `render_template()` from route handler
- CSRF: Include `{{ csrf_token() }}` in forms
- Metadata injection: Use `@app.context_processor` in `app.py` if template-wide data needed

## Special Directories

**data/ (Reference Data):**
- Purpose: Contains immutable reference datasets for enrichment analysis
- Generated: No (source data for application)
- Committed: Yes (required for function)
- Lifecycle: Loaded once at startup, cached with 1-hour TTL
- Structure: CSV files (KE mappings), WikiPathways XGMML, demo datasets

**uploads/ (User Files):**
- Purpose: Temporary storage during analysis
- Generated: Yes (at runtime, user uploads)
- Committed: No (transient)
- Lifecycle: Created on file upload, deleted after analysis via `cleanup_file()`
- File types: CSV, TSV, TXT (validated via `Config.ALLOWED_EXTENSIONS`)

**temp/ (Temporary Data):**
- Purpose: Temporary storage of network PNG snapshots
- Generated: Yes (at runtime, client uploads)
- Committed: No (transient)
- Lifecycle: Created on network export (`/api/upload_network_png`), deleted after report generation or after 1 hour (`cleanup_old_temp_files()`)
- Format: Base64-encoded PNGs, filename format: `{uuid}.png`

**.planning/codebase/ (GSD Planning Documents):**
- Purpose: Architecture and patterns documentation for GSD command integration
- Generated: Yes (by GSD mapping command)
- Committed: Yes (guides future code generation)
- Contains: ARCHITECTURE.md, STRUCTURE.md, CONVENTIONS.md, TESTING.md, CONCERNS.md as applicable

---

*Structure analysis: 2026-02-25*
