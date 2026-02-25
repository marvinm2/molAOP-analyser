# Coding Conventions

**Analysis Date:** 2026-02-25

## Naming Patterns

**Files:**
- Service classes: `snake_case` with `_service` suffix → `enrichment_service.py`, `column_detector.py`, `data_service.py`
- Modules/utilities: `snake_case` → `validation.py`, `helpers.py`, `exceptions.py`, `config.py`
- Tests: `test_` prefix followed by module name → `test_flask_routes.py`, `test_column_detector.py`, `test_database.py`

**Functions:**
- Snake case throughout → `load_reference_sets()`, `validate_file_upload()`, `run_enrichment_analysis()`
- Private functions: Single leading underscore → `_detect_gene_id_columns()`, `_score_column_name()`
- Public API functions: No leading underscore

**Variables:**
- Local variables and parameters: `snake_case` → `df_processed`, `ke_list`, `reference_sets`, `user_gene_status`
- Constants: `UPPER_CASE` → `GENE_ID_PATTERNS`, `MAX_FILE_SIZE`, `ALLOWED_EXTENSIONS`
- Loop variables: Single letters acceptable for simple iterations → `for ke, ref_genes in filtered_reference_sets.items()`

**Types/Classes:**
- Classes: `PascalCase` → `ColumnDetector`, `ColumnMatch`, `ColumnSuggestions`, `AOPAnalysisError`
- Dataclasses: `PascalCase` → `ExperimentMetadata`, `ReportData`, `ExperimentRecord`
- Exception classes: `PascalCase` with `Error` suffix → `DataValidationError`, `FileProcessingError`, `EnrichmentAnalysisError`

## Code Style

**Formatting:**
- No explicit formatter configured (no `.prettierrc` or equivalent)
- Follow Python standard conventions (PEP 8 compliant)
- Line length appears to follow standard conventions
- Indentation: 4 spaces (Python standard)

**Linting:**
- No `.flake8`, `.pylintrc`, or similar config found
- Code follows PEP 8 naming and structure conventions
- Focus on readability and maintainability rather than strict linting rules

## Import Organization

**Order:**
1. Standard library imports → `os`, `sys`, `logging`, `json`, `math`, `re`, `tempfile`, `datetime`, `uuid`, `base64`, `time`
2. Third-party imports → `pandas`, `numpy`, `scipy`, `flask`, `sqlalchemy`, `werkzeug`, `plotly`, etc.
3. Local/relative imports → From `config`, `validation`, `helpers`, `services`, `database`, `exceptions`

**Path Aliases:**
- Direct imports without aliases for most modules
- Relative imports within services: `from services.gene_id_validator import ...`
- Exception imports grouped together: `from exceptions import AOPAnalysisError, format_error_response`
- Conditional imports used for optional dependencies (e.g., report generation libraries)

## Error Handling

**Patterns:**
- Custom exception hierarchy rooted in `AOPAnalysisError` base class defined in `exceptions.py`
- Specific exception types for different contexts:
  - `DataValidationError` (with `field` and `value` attributes) → `services/data_service.py:38`
  - `FileProcessingError` (with `filename` and `line_number` attributes) → `services/data_service.py:63-68`
  - `EnrichmentAnalysisError` (with `ke_id` and `gene_count` attributes)
  - `AOPDataError` (with `aop_id` and `missing_data` attributes)
  - `NetworkBuildError` (with `node_count` and `edge_count` attributes)
  - `ConfigurationError` (with `config_key` attribute)

- Re-raise pattern: Catch and check exception type, then re-raise or wrap
  ```python
  except DataValidationError:
      raise  # Re-raise validation errors as-is
  except pd.errors.EmptyDataError:
      raise FileProcessingError(f"File is empty: {filepath}", filename=filepath)
  except Exception as e:
      logger.error(f"Failed to load data: {e}")
      raise FileProcessingError(f"Unexpected error: {e}", filename=filepath)
  ```
  See `services/data_service.py:60-68`

- User-friendly error messages: Centralized in `exceptions.py` USER_ERROR_MESSAGES dictionary with `format_error_response()` function
- Try-except blocks limited to specific operations; avoid bare `except` clauses

## Logging

**Framework:** Python's built-in `logging` module with `logging.getLogger(__name__)`

**Patterns:**
- Initialize logger at module level: `logger = logging.getLogger(__name__)` → Used in all major modules
- Log levels:
  - `logger.info()` for significant events → `"Starting enrichment analysis"`, `"Loaded {len(df_raw)} rows from {filepath}"`
  - `logger.debug()` for detailed flow information → `"No overlap found for KE {ke}"`
  - `logger.warning()` for edge cases → `"KE {ke}: No significant genes found - result may not be meaningful"`
  - `logger.error()` for failures with `exc_info=True` → `logger.error(f"Internal server error: {e}", exc_info=True)`

- Flask app logging configured in `app.py:32-36`:
  ```python
  logging.basicConfig(
      level=logging.INFO,
      format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
  )
  ```

## Comments

**When to Comment:**
- Function/method docstrings: Triple-quoted, always present → See `services/enrichment_service.py:18-29`
- Parameter explanations in docstrings: List each parameter with type and purpose
- Complex logic blocks: Inline comments for non-obvious operations → `services/enrichment_service.py:69-73` shows contingency table construction
- Edge cases: Document why certain conditions are handled → `services/enrichment_service.py:61-64`
- NOT: Don't comment obvious code (`x = 5  # Set x to 5` is redundant)

**JSDoc/TSDoc:**
- Not applicable (Python project, not JavaScript/TypeScript)
- Use Python docstring format: Google-style or NumPy-style (this codebase uses Google-style)
- Example from `services/enrichment_service.py:18-29`:
  ```python
  """
  Run Fisher's exact test enrichment analysis for Key Events.

  Args:
      df: Processed gene expression dataframe with 'significant' column
      reference_sets: Dictionary mapping KE_ID to gene sets
      ke_list: Set of KE IDs to analyze
      ke_title_map: Mapping of KE IDs to titles

  Returns:
      pd.DataFrame: Enrichment results sorted by FDR
  """
  ```

## Function Design

**Size:** Functions typically 20-100 lines; service methods are longer for multi-step operations
- `run_enrichment_analysis()` → ~120 lines (includes loop and statistical tests)
- `detect_columns()` → ~20 lines (orchestrator that calls private methods)
- `_detect_gene_id_columns()` → ~50-80 lines (detailed analysis)

**Parameters:**
- Positional parameters for required inputs
- Type hints included: `pd.DataFrame`, `Dict[str, Set[str]]`, `Optional[str]`
- Minimal default values (used sparingly for optional features)
- Example from `services/data_service.py:14`:
  ```python
  def load_and_validate_data(filepath: str, id_col: str, fc_col: str, pval_col: str) -> pd.DataFrame:
  ```

**Return Values:**
- Single return value (scalar or object) when possible
- Tuple returns for multiple related values: `(processed_dataframe, stats_dict)` → `services/data_service.py:70`
- DataFrame returns for data-heavy operations
- Dictionary returns for structured results
- Implicit `None` for functions with side effects (e.g., database writes)

## Module Design

**Exports:**
- Service classes instantiated as module-level singletons: `column_detector = ColumnDetector()`, `report_generator = ReportGenerator()`
- Functions exported directly: `run_enrichment_analysis`, `build_cytoscape_network`
- Exceptions exported from `exceptions.py` for use throughout application

**Barrel Files:**
- Not used; imports are explicit and specific
- `services/__init__.py` exists but is empty (services imported directly)
- Encourages direct imports: `from services.enrichment_service import run_enrichment_analysis` rather than `from services import run_enrichment_analysis`

## Data Structure Conventions

**DataFrames:**
- Column naming: Mixed case in data, normalized to specific names during processing
  - User data might have: `Gene_Symbol`, `log2FoldChange`, `padj`, `baseMean`
  - Normalized to: `ID` (gene identifier), `log2FC` (log fold change), `pval` (p-value), `significant` (boolean flag)
  - See `services/data_service.py:50-52`

- Consistent gene ID handling: Always uppercase and stripped
  ```python
  df_processed['ID'] = df_processed['ID'].astype(str).str.strip().str.upper()
  ```
  → `services/data_service.py:97`

**Configuration:**
- Config class constants for all application settings → `config.py:10-30`
- ExperimentMetadata dataclass for typed metadata structure → `config.py:70+`
- Explicit dictionary structures for reference data and results

---

*Convention analysis: 2026-02-25*
