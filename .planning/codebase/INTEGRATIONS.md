# External Integrations

**Analysis Date:** 2026-02-25

## APIs & External Services

**WikiPathways Integration:**
- WikiPathways - Gene pathway annotation source
  - Data: `data/wikipathways_hsa_20240410.xgmml` (static XGMML network file, 18MB)
  - Purpose: Maps pathways (WP_ID) to gene symbols for enrichment gene sets
  - Client: Custom XGMML parser in `services/data_service.py`

**No External API Calls:**
- This application does NOT make live API calls to external services
- All reference data is pre-downloaded and bundled in the `data/` directory
- No webhook receivers or real-time integrations detected

## Data Storage

**Databases:**
- SQLite (local file-based)
  - Location: `molAOP_analyser.db`
  - Client: SQLAlchemy ORM
  - Purpose: Stores experiment metadata, analysis parameters, and enrichment results
  - Models: `ExperimentRecord` in `database.py`

**File Storage:**
- Local filesystem only
  - Uploads: `uploads/` directory (10MB max per file)
  - Temp files: `temp/` directory (network PNGs, report generation artifacts)
  - Cleanup: Automatic via `utils.cleanup_file()` after processing

**Caching:**
- In-memory cache (application runtime only)
  - Implementation: `CacheManager` in `cache_manager.py`
  - Purpose: Cache reference data (KE-gene mappings) and analysis results
  - TTL support: Configurable, default 3600 seconds (1 hour)
  - Thread-safe with lock-based synchronization

## Authentication & Identity

**Auth Provider:**
- None - This is an internal research tool
- CSRF Protection: Enabled via Flask-WTF
  - Token generation: `flask_wtf.csrf.generate_csrf()`
  - Applied to all forms in `app.py` context processor

**Session Management:**
- Flask sessions with `SECRET_KEY` configuration
- No user authentication system detected
- Sessions used for request-scoped state (e.g., analysis parameters, file paths)

## Monitoring & Observability

**Error Tracking:**
- None detected - No Sentry, Rollbar, or similar integration

**Logs:**
- Python logging module (stdlib)
  - Configuration: `basicConfig` in `app.py` (lines 31-35)
  - Format: `%(asctime)s - %(name)s - %(levelname)s - %(message)s`
  - Level: INFO
  - Loggers in: `app.py`, `database.py`, `cache_manager.py`, all `services/*.py`

**Debugging:**
- Debug print statements in `helpers.py` for data validation
- No external debugging service integration

## CI/CD & Deployment

**Hosting:**
- Docker container (recommended deployment)
  - Dockerfile: Python 3.11-slim base with pip dependencies
  - Exposed port: 5000 (Flask development server)

**CI Pipeline:**
- None detected - No GitHub Actions, GitLab CI, or Jenkins configuration

**Local Testing:**
- pytest for test execution
- Commands: `pytest`, `pytest -v`, `pytest --cov`

## Environment Configuration

**Required env vars:**
- `SECRET_KEY` - Flask session encryption (optional, defaults to 'dev-key-change-in-production')
- `FLASK_ENV` - Set to 'development' for docker-compose (optional)

**Secrets location:**
- Environment variables only (no .env file checked in)
- `.env` file present but NOT committed (listed in .gitignore)

## Webhooks & Callbacks

**Incoming:**
- None detected - No webhook endpoints

**Outgoing:**
- None detected - No outbound API calls or event notifications

## Data File Dependencies

**Reference Data (required for runtime):**
- `data/aop_ke_map.csv` - Maps AOP IDs to KE IDs
- `data/aop_ker_edges.csv` - KE relationship edges (for network topology)
- `data/KE-WP.csv` - KE → WikiPathways ID mappings
- `data/edges_wpid_to_gene.csv` - WikiPathways ID → gene ID mappings (715 KB)
- `data/node_attributes.csv` - Gene ID → gene symbol mappings (544 KB)
- `data/ke_metadata.csv` - KE titles and type classifications
- `data/wikipathways_hsa_20240410.xgmml` - XGMML network file (18 MB)

**Demo Datasets (bundled):**
- `data/GSE90122_TO90137.tsv` - PXR agonist study (7.5 MB)
- `data/GSE90122_SR12813.tsv` - PXR agonist study (7.4 MB)
- `data/Cisplatin_Kidney/` - 42 cisplatin toxicity datasets (4hr-72hr timepoints, 0.1-50μM doses)

All validated on startup via `Config.validate_data_files()` in `config.py` (lines 111-121).

---

*Integration audit: 2026-02-25*
