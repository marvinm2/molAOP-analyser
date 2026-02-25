# Technology Stack

**Analysis Date:** 2026-02-25

## Languages

**Primary:**
- Python 3.11 - Main application language (specified in Dockerfile)

## Runtime

**Environment:**
- Python 3.11-slim (Docker base image)
- Flask WSGI web framework

**Package Manager:**
- pip - Python package manager
- Lockfile: Not detected (uses requirements.txt)

## Frameworks

**Core:**
- Flask - Web framework for HTTP routing, templating, and request handling
- Flask-WTF - CSRF protection for forms

**Data Processing:**
- pandas - Tabular data manipulation and analysis
- numpy - Numerical computing
- scipy - Statistical functions (Fisher's exact test, p-value combination)
- statsmodels - FDR correction (Benjamini-Hochberg method via `multipletests`)

**Visualization:**
- Plotly - Interactive charting for volcano plots and data visualization
- Cytoscape.js (frontend) - Interactive network visualization for AOP pathways

**Reporting & PDF Generation:**
- WeasyPrint - HTML-to-PDF conversion (primary report backend)
- ReportLab - PDF generation with tables, styles, and layout control
- svglib - SVG-to-PDF conversion support
- kaleido 0.2.1 - Plotly static image export (pinned version)

**Database:**
- SQLAlchemy - ORM for database models and queries
- SQLite - Database backend (file-based, `molAOP_analyser.db`)

**Testing:**
- pytest - Test runner and framework
- pytest-flask - Flask-specific fixtures and utilities
- pytest-mock - Mocking support for unit tests
- pytest-cov - Code coverage reporting

**Utilities:**
- Werkzeug - WSGI utilities (file upload handling, security filename validation)
- requests - HTTP client library (used in integration tests)

## Key Dependencies

**Critical:**
- scipy - Enables Fisher's exact test statistical enrichment analysis
- statsmodels - Required for FDR (Benjamini-Hochberg) p-value correction
- SQLAlchemy - Persistence layer for experiment metadata and analysis history
- pandas - Core data processing for gene expression datasets

**Infrastructure:**
- Flask - Web server and routing
- WeasyPrint/ReportLab - Report generation (PDF exports)
- Plotly - Interactive visualizations embedded in reports

## Configuration

**Environment:**
- `SECRET_KEY` - Flask session security (environment variable or 'dev-key-change-in-production' default)
- `FLASK_ENV` - Development/production mode flag (set via docker-compose)

**Build:**
- `Dockerfile` - Python 3.11-slim base with pip install
- `docker-compose.yml` - Single-service setup with volume mounts for development

**File Limits:**
- `MAX_FILE_SIZE`: 10 MB per uploaded file
- `MAX_CONTENT_LENGTH`: 50 MB for full form data (reports with embedded visualizations)
- `ALLOWED_EXTENSIONS`: csv, tsv, txt

## Platform Requirements

**Development:**
- Python 3.11 or Docker with Python 3.11-slim
- System dependencies: `build-essential` (for C compilation during pip install)
- Minimum 50 MB available for form data buffering

**Production:**
- Docker container (recommended)
- 5000 port exposure for Flask development server
- SQLite file storage for `molAOP_analyser.db`
- Temporary file storage for uploads and report generation (`uploads/`, `temp/` directories)

---

*Stack analysis: 2026-02-25*
