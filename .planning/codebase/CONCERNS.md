# Codebase Concerns

**Analysis Date:** 2026-02-25

## Tech Debt

**Debug Statements in Production Code:**
- Issue: Multiple debug print statements left in `helpers.py` that output to console instead of logger
- Files: `helpers.py` (lines 9-58)
- Impact: Pollutes stdout, no log level control, unprofessional output visible to users/monitoring systems
- Fix approach: Remove all `print()` statements and use `logger.debug()` instead. Already have logging configured in app.py

**Duplicate `guess_id_type()` Function:**
- Issue: Function defined in both `app.py` (lines 379-408) and `services/data_service.py` (lines 124-152)
- Files: `app.py`, `services/data_service.py`
- Impact: Code duplication creates maintenance burden; changes to one must be synchronized with the other
- Fix approach: Remove from `app.py` and always import from `services/data_service.py`

**Unused Legacy Network Function:**
- Issue: `build_cytoscape_network_legacy()` function exists but is not called anywhere
- Files: `services/network_service.py` (lines 94-149)
- Impact: Dead code increases maintenance burden and codebase size
- Fix approach: Remove deprecated function if not needed for backward compatibility, otherwise document why it's retained

**Hardcoded Paths in Data Loading:**
- Issue: Data file paths are hardcoded as relative paths throughout the codebase (`'data/...'`)
- Files: `app.py` (line 779-782), `services/data_service.py` (lines 166, 170, 174), `helpers.py` (implicit)
- Impact: Brittle deployment; fails if cwd differs, won't work in packaged/containerized environments
- Fix approach: Use `os.path.join(os.path.dirname(__file__), 'data', ...)` or configurable data directory

**Session-Based Metadata Storage:**
- Issue: Experiment metadata stored in Flask session without explicit cleanup
- Files: `app.py` (lines 232, 420-422, 517, 622-623)
- Impact: Sessions grow unbounded; session data persists across requests; memory leak in long-running deployments
- Fix approach: Use database for persistence instead of session; clear session after each analysis

**Inconsistent Error Handling in Services:**
- Issue: Services catch broad `Exception` without re-raising, logging at different levels
- Files: `services/enrichment_service.py` (line 103-105), `services/network_service.py` (line 147-149)
- Impact: Errors silently swallowed; inconsistent logging makes debugging difficult
- Fix approach: Standardize on raising specific custom exceptions, log errors consistently

## Known Bugs

**Unicode Character Handling in Demo Datasets:**
- Symptoms: Demo dataset names include Unicode characters (μM, °) that may cause issues in some browsers or file systems
- Files: `config.py` (lines 36-89)
- Trigger: Selecting demo datasets with non-ASCII characters; accessing via file path
- Workaround: Use numeric equivalents (uM instead of μM) or ensure UTF-8 support throughout stack

**JSON Decoding in Report Generation:**
- Symptoms: Report generation fails with JSON decode errors when enrichment results contain special characters
- Files: `app.py` (lines 649-672)
- Trigger: Gene names or KE titles with quotes, backslashes, or control characters; form submission and re-parsing
- Workaround: HTML unescape is applied (line 658) but may not handle all cases
- Fix approach: Pre-serialize JSON with `json.dumps()` before putting in form; deserialize once on receipt

**Metadata Upload Timestamp Not Set:**
- Symptoms: Database records may have missing upload_timestamp values
- Files: `database.py` (line 144-145)
- Trigger: When metadata.get('upload_timestamp') returns None or invalid format
- Workaround: Fallback to datetime.utcnow() is in place but uses UTC instead of configured timezone
- Fix approach: Always set timestamp in session when metadata is created (use ExperimentMetadata dataclass)

## Security Considerations

**Weak Default Secret Key:**
- Risk: Flask SECRET_KEY defaults to hardcoded string if not in environment variable
- Files: `config.py` (line 108)
- Current mitigation: Falls back to 'dev-key-change-in-production'
- Recommendations: Force environment variable in production; rotate key on deployment; add validation in config initialization

**File Path Traversal Prevention - Partial:**
- Risk: Demo file selection uses normpath() but could still access files outside data directory with symlinks
- Files: `app.py` (lines 257-273)
- Current mitigation: Checks for '..' and leading '/' in normalized path
- Recommendations: Use os.path.realpath() instead of normpath(); validate resolved path starts with data directory; use whitelist of allowed files

**SQL Injection in Database Search:**
- Risk: `search_experiments()` uses `.contains()` for LIKE queries without parameterization (but SQLAlchemy handles this)
- Files: `database.py` (lines 238-243)
- Current mitigation: SQLAlchemy ORM automatically parameterizes queries
- Recommendations: Document this assumption; never switch to raw SQL; add validation for search input length

**Base64 PNG Encoding Size Limit:**
- Risk: No size limit on base64-encoded PNG data submitted to `/api/upload_network_png`
- Files: `app.py` (lines 176-213)
- Current mitigation: General MAX_CONTENT_LENGTH limit (50MB) applies
- Recommendations: Add specific size check for PNG data; validate PNG format before saving; set reasonable limit on decoded size

**Temporary File Cleanup Race Condition:**
- Risk: Temporary PNG files cleaned before report generation completes if timing is tight
- Files: `app.py` (lines 185-186, 688-693)
- Current mitigation: Tries to clean up after reading; catches exceptions
- Recommendations: Use context managers or temp file library; guarantee cleanup only after reading; add temp file retention policy

## Performance Bottlenecks

**Reference Set Loading on Every Request:**
- Problem: `load_cached_reference_sets()` must reload large CSV files even with caching
- Files: `app.py` (lines 768-792), `helpers.py` (lines 3-67)
- Cause: Cache TTL is 1 hour; reference data is large; merge operations are expensive
- Improvement path: Increase TTL to 24 hours; pre-load on startup; use pickle or database for faster deserialization; consider lazy loading per KE instead of all-at-once

**Volcano Plot Data Generation - 10,000 Gene Limit:**
- Problem: MAX_GENES_DISPLAY = 10,000 but volcano plot limited to 2,000 (line 478) and preview to unlimited
- Files: `app.py` (lines 354, 478), `config.py` (line 19)
- Cause: Inconsistent limits; large volcano plots cause form submission delays
- Improvement path: Use pagination or sampling for display; pre-render to SVG; implement client-side filtering

**Full DataFrame in Memory:**
- Problem: Entire gene expression dataframe stays in memory throughout analysis
- Files: `services/data_service.py` (lines 45-56, 82-94)
- Cause: No streaming or chunked processing
- Improvement path: For very large datasets (>100K genes): stream processing, or only keep significant genes in memory

**DataFrame Grouping Without Index:**
- Problem: `process_gene_expression()` groups genes repeatedly without creating index first
- Files: `services/data_service.py` (line 82)
- Cause: Inefficient iteration over grouped data
- Improvement path: Set 'ID' as index before grouping; use vectorized operations for p-value combining

**Network Visualization - No Edge Case Optimization:**
- Problem: `build_cytoscape_network()` iterates through all edges twice (once per node lookup)
- Files: `services/network_service.py` (lines 34-45, 69-83)
- Cause: Linear lookup for enrichment results and node membership
- Improvement path: Create index/set before loop; pre-filter edges to valid KEs; use vectorized pandas operations

## Fragile Areas

**Column Detection and Selection:**
- Files: `services/column_detector.py`, `app.py` (lines 303-324)
- Why fragile: ML-based detection may fail on non-standard formats; auto-selection without user confirmation could select wrong columns
- Safe modification: Test edge cases (very few rows, mixed data types, unusual column names); add confidence thresholds; always allow manual override
- Test coverage: `tests/test_column_detector.py` exists but coverage gaps likely in edge cases

**Enrichment Analysis Edge Cases:**
- Files: `services/enrichment_service.py` (lines 52-105)
- Why fragile: Minimum 5 genes required (line 62) is arbitrary; contingency table can have zeros; FDR correction needs at least 2 results
- Safe modification: Add parameter for minimum gene threshold; handle edge cases (0/0 contingency tables, single result); log warnings
- Test coverage: Unit tests exist but edge case handling needs verification

**Report Generation - Multiple Optional Dependencies:**
- Files: `services/report_service.py` (lines 22-57)
- Why fragile: Three competing libraries (Plotly, WeasyPrint, ReportLab) with optional imports; fallbacks may not be fully tested
- Safe modification: Document which report format requires which library; test without optional libs; provide clear error messages when libraries missing
- Test coverage: `tests/test_report_service.py` exists but likely doesn't test all dependency combinations

**Database Session Management:**
- Files: `database.py` (lines 90-110, 116-120)
- Why fragile: Sessions created per operation; no connection pooling; no timeout handling; SQLAlchemy error handling catches broad exceptions
- Safe modification: Use context managers; implement proper session cleanup; add connection retry logic; handle timeout errors specifically
- Test coverage: `tests/test_database.py` exists but production load testing likely missing

**AOP Data Loading - No Validation:**
- Files: `services/data_service.py` (lines 154-184)
- Why fragile: No validation that loaded DataFrames have expected columns; no check that AOP ID actually exists
- Safe modification: Validate dataframe structure before returning; check ke_list not empty; validate edge consistency
- Test coverage: Integration tests likely assume well-formed data

## Scaling Limits

**In-Memory Cache - No Bounds:**
- Current capacity: Unbounded dictionary in `cache_manager.py`
- Limit: Grows until memory exhausted; no eviction policy for expired entries
- Scaling path: Implement LRU cache or size-based eviction; use Redis for distributed deployment; monitor cache growth

**Database - SQLite File-Based:**
- Current capacity: Single SQLite file on disk; no connection pooling
- Limit: Single writer at a time; locks during writes; not suitable for concurrent usage
- Scaling path: Migrate to PostgreSQL/MySQL for production; add connection pooling; implement read replicas

**Temporary File Accumulation:**
- Current capacity: Unbounded temp directory for network PNGs
- Limit: Disk space exhaustion; cleanup only every hour
- Scaling path: Implement daily/weekly cleanup; set disk quota; stream files directly to client instead of temp storage

**Form Data Size - 50MB Limit:**
- Current capacity: MAX_CONTENT_LENGTH = 50MB including enrichment results JSON
- Limit: Large analyses (1000+ KEs) may exceed limit; base64 PNG adds 33% overhead
- Scaling path: Stream report generation instead of form submission; use separate upload endpoint for network PNG; implement chunked encoding

## Dependencies at Risk

**WeasyPrint - System Dependencies:**
- Risk: Requires system libraries (libpango, libpangoft2, libfreetype6, libffi6) to generate PDFs
- Impact: Docker builds may fail; deployment requires OS-specific configuration
- Migration plan: Make PDF generation truly optional; provide warnings in logs if unavailable; implement pure-Python alternative using reportlab-only

**ReportLab - Complex SVG Support:**
- Risk: SVG support (`svglib`) is optional and fragile; may not render network visualizations correctly
- Impact: Network PNG fallback may be only working option; code path not well-tested
- Migration plan: Test all three report formats thoroughly; document which features work with which backends; consider single-backend approach

**pandas - Large File Performance:**
- Risk: `.read_csv()` with `engine='python'` is slower than 'c' engine; no streaming support
- Impact: Slow file loading for large datasets; entire file in memory
- Migration plan: Switch to 'c' engine; consider polars or dask for very large files; add progress indication

**scipy/statsmodels - Statistical Correctness:**
- Risk: FDR correction uses Benjamini-Hochberg; Fisher's exact test requires >= 1 observation in contingency table
- Impact: Edge cases (all significant genes, no significant genes) may not be handled correctly
- Migration plan: Add statistical validation; document assumptions; consider robust alternatives like permutation tests

## Missing Critical Features

**No User Authentication:**
- Problem: Application is fully open; no user tracking or access control
- Blocks: Cannot implement audit trails, user quotas, or personalization
- Workaround: Deploy behind authentication proxy; document for users that data is not private

**No Persistent Analysis History:**
- Problem: Database exists but results are not accessible through UI; only top 20 enrichment results saved
- Blocks: Users cannot browse past analyses; full enrichment results are lost after session ends
- Workaround: Save full enrichment results to database; implement results browser UI

**No Batch Processing:**
- Problem: Single dataset per request; no ability to analyze multiple files or AOPs in sequence
- Blocks: Cannot perform comparative analyses across conditions; users must upload manually multiple times
- Improvement path: Implement queue-based batch processing; add results comparison dashboard

**No Real-Time Analysis Feedback:**
- Problem: No progress indication during analysis; form submission waits for completion
- Blocks: Long analyses appear to hang; users cannot cancel
- Improvement path: Implement WebSocket/Server-Sent Events for progress updates; add job ID and status polling

## Test Coverage Gaps

**File Upload Security:**
- What's not tested: Malformed files, extremely large files, files with special characters in names
- Files: `validation.py`, `app.py` (preview route)
- Risk: Could accept invalid files or crash; filename injection possible
- Priority: High (security-relevant)

**Error Handler Routes:**
- What's not tested: 404, 500 handlers; behavior when required services fail
- Files: `app.py` (lines 77-100)
- Risk: Error pages may break the application further; cascading failures not prevented
- Priority: Medium (user experience)

**Database Failure Scenarios:**
- What's not tested: Database connection loss, transaction rollback, constraint violations
- Files: `database.py`
- Risk: Silent failures; data inconsistency; no retry logic
- Priority: High (data integrity)

**Report Generation with Missing Data:**
- What's not tested: Missing network_png, empty enrichment results, invalid volcano data
- Files: `services/report_service.py`, `app.py` (lines 708-726)
- Risk: Report generation fails silently or produces corrupt output
- Priority: Medium (feature completeness)

**Cache Expiration:**
- What's not tested: TTL expiration, concurrent access, cache hit/miss during updates
- Files: `cache_manager.py`
- Risk: Stale data served; race conditions in threaded environment
- Priority: Medium (data freshness)

**Gene ID Validation:**
- What's not tested: Invalid formats, non-existent gene IDs, mixed identifier types
- Files: `services/gene_id_validator.py`
- Risk: Silent failures to validate; incorrect enrichment results downstream
- Priority: High (analysis correctness)

**Column Type Guessing Edge Cases:**
- What's not tested: Columns with mostly valid but some invalid data; columns that match multiple patterns equally
- Files: `services/column_detector.py`
- Risk: Wrong column selected; analysis produces meaningless results
- Priority: High (user impact)

---

*Concerns audit: 2026-02-25*
