# Testing Patterns

**Analysis Date:** 2026-02-25

## Test Framework

**Runner:**
- pytest [version in requirements.txt]
- Config: `pytest.ini`

**Assertion Library:**
- Standard pytest assertions (`assert` statements)

**Run Commands:**
```bash
pytest                          # Run all tests
pytest -v                       # Verbose output with test names
pytest --cov                    # Coverage report (configured for 80% minimum)
pytest tests/test_flask_routes.py  # Run specific test file
pytest -m unit                  # Run only unit tests
pytest -m integration           # Run only integration tests
pytest -m database              # Run only database tests
pytest -m web                   # Run only web/Flask tests
```

## Test File Organization

**Location:**
- Main test suite: `tests/` directory
- Test discovery: Files named `test_*.py` in `tests/` directory
- Standalone integration tests: Root level `test_integration.py` for end-to-end workflows

**Naming:**
- Test files: `test_{module_name}.py` → `test_flask_routes.py`, `test_column_detector.py`, `test_database.py`, `test_report_service.py`
- Test classes: `Test{ModuleName}` → `TestFlaskRoutes`, `TestColumnDetector`, `TestDatabaseOperations`
- Test functions: `test_{functionality}` → `test_index_route`, `test_gene_id_column_detection`, `test_save_experiment_metadata`

**Structure:**
```
tests/
├── conftest.py                 # Shared fixtures and configuration
├── test_flask_routes.py        # Integration/web tests
├── test_column_detector.py     # Unit tests
├── test_database.py            # Unit/database tests
├── test_report_service.py      # Unit tests
└── data/                       # Test data files (if needed)
```

## Test Structure

**Suite Organization:**
```python
@pytest.mark.unit
class TestColumnDetector:
    """Test column detection functionality."""

    def test_column_detector_initialization(self, column_detector):
        """Test ColumnDetector initialization."""
        assert column_detector.min_confidence == 0.3
        assert column_detector.high_confidence == 0.8
```
See `tests/test_column_detector.py:10-20`

**Patterns:**
- Class-based organization for related tests using `Test*` prefix
- Docstrings for classes and methods explaining what is tested
- One assertion per test preferred, multiple assertions acceptable for related checks
- Fixtures injected as parameters with `@pytest.fixture` decorator

**Setup/Teardown:**
- Fixture-based setup: `@pytest.fixture` decorator with optional `scope` parameter
  ```python
  @pytest.fixture(scope='session')
  def test_data_dir():
      """Fixture providing path to test data directory."""
      return os.path.join(os.path.dirname(__file__), 'data')
  ```
  See `tests/conftest.py:27-30`

- Fixture cleanup: `yield` pattern with cleanup code after yield
  ```python
  @pytest.fixture
  def temp_database():
      """Fixture providing a temporary test database."""
      with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
          db_path = f.name

      db_manager = DatabaseManager(db_url=f"sqlite:///{db_path}")
      db_manager.initialize()

      yield db_manager  # Setup complete, test runs

      # Cleanup
      try:
          os.unlink(db_path)
      except OSError:
          pass
  ```
  See `tests/conftest.py:86-102`

## Mocking

**Framework:** `unittest.mock` from Python standard library

**Patterns:**
- Import: `from unittest.mock import patch, MagicMock`
- Context manager usage: `with patch(...) as mock_obj:`
- Multiple patches: Use nested `with` statements or multiple patch decorators
  ```python
  with patch('app.load_and_validate_data') as mock_load, \
       patch('app.process_gene_expression') as mock_process, \
       patch('app.load_aop_data') as mock_aop_data:
      # Test code
  ```
  See `tests/test_flask_routes.py:55-66`

- MagicMock for complex return values:
  ```python
  mock_df = MagicMock()
  mock_df.head.return_value.to_dict.return_value = [
      {'Gene_Symbol': 'BRCA1', 'log2FoldChange': 2.5, 'padj': 0.001}
  ]
  mock_df.columns.tolist.return_value = ['Gene_Symbol', 'log2FoldChange', 'padj']
  ```
  See `tests/test_flask_routes.py:29-34`

- Patch function returns: `.return_value` attribute
  ```python
  mock_enrichment.return_value = MagicMock()
  mock_enrichment.return_value.to_dict.return_value = []
  ```
  See `tests/test_flask_routes.py:74-75`

**What to Mock:**
- External dependencies (file I/O, database, network calls)
- Flask-specific components when testing business logic in isolation
- Service calls when testing individual components
- Expensive operations (file processing, complex calculations) in unit tests

**What NOT to Mock:**
- Core data structures (DataFrames, dictionaries)
- Simple utility functions
- Application-specific dataclasses (ExperimentMetadata, ColumnMatch)
- Actual business logic being tested (enrichment analysis, column detection)

## Fixtures and Factories

**Test Data:**
Fixtures defined in `tests/conftest.py`:

```python
@pytest.fixture
def sample_gene_data():
    """Fixture providing sample gene expression data."""
    return pd.DataFrame({
        'Gene_Symbol': ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS'],
        'log2FoldChange': [2.5, -1.8, 3.2, -2.1, 1.9],
        'padj': [0.001, 0.005, 0.0001, 0.01, 0.03],
        'baseMean': [1000, 800, 1200, 600, 900]
    })
```
See `tests/conftest.py:33-41`

**Other common fixtures:**
- `sample_metadata()` → ExperimentMetadata instance → `tests/conftest.py:44-53`
- `sample_report_data()` → ReportData instance with complete test data → `tests/conftest.py:56-83`
- `temp_database()` → Temporary SQLite database for isolation → `tests/conftest.py:86-102`
- `flask_client()` → Flask test client with testing config → `tests/conftest.py:105-115`
- `authenticated_client()` → Flask client with session data → `tests/conftest.py:118+`
- `column_detector` → ColumnDetector instance → Instantiated in fixture
- `test_data_dir` → Path to test data directory → `tests/conftest.py:27-30`

**Location:**
- Central location: `tests/conftest.py`
- Session-scoped: `test_data_dir` (once per test session)
- Function-scoped: Most fixtures (fresh for each test)
- Available to all test files automatically via pytest discovery

## Coverage

**Requirements:**
- Minimum 80% coverage enforced via pytest configuration
- Configured in `pytest.ini`: `--cov-fail-under = 80`

**Coverage Exclusions:**
```
--cov-exclude=tests/*
--cov-exclude=test_*.py
```

**View Coverage:**
```bash
pytest --cov                    # Terminal report with missing lines
pytest --cov --cov-report=html  # HTML report in htmlcov/
pytest --cov --cov-report=term-missing  # Terminal with missing line numbers
```

## Test Types

**Unit Tests:**
- Scope: Individual functions/methods in isolation
- Marked with: `@pytest.mark.unit`
- Examples:
  - Column detection logic → `tests/test_column_detector.py:10-96`
  - Database model creation → `tests/test_database.py:16-32`
  - Report data structure → `tests/test_report_service.py:27-39`
- Use mocking to isolate from dependencies
- Fast execution (< 1 second per test typically)

**Integration Tests:**
- Scope: Multiple components working together
- Marked with: `@pytest.mark.integration`
- Examples:
  - Flask route + data loading + column detection → `tests/test_flask_routes.py:10`
  - Metadata capture + database storage + analysis
- Limited mocking; more realistic data flow
- Slower execution (1-5 seconds per test typically)

**E2E Tests:**
- Framework: None formally configured (tests/test_integration.py is standalone)
- Scope: Complete workflow from file upload through report generation
- Location: `test_integration.py` (not in pytest discovery, manual execution)
- Run: `python test_integration.py` or similar
- Tests actual Flask server running on `http://localhost:5000`
- Example flow:
  1. POST to `/preview` with demo file
  2. POST to `/analyze` with column selections
  3. POST to `/generate_report` with parameters
  4. Validate HTML/PDF output

## Common Patterns

**Async Testing:**
Not applicable - application is synchronous Flask/Pandas-based

**Error Testing:**
```python
def test_preview_route_missing_file(self, flask_client):
    """Test preview route with missing file."""
    response = flask_client.post('/preview', data={})

    assert response.status_code == 400
    assert b'No dataset provided' in response.data
```
See `tests/test_flask_routes.py:48-53`

- Assert HTTP status codes for error conditions
- Assert error messages in response data
- Use Flask test client for route error testing
- Validation errors return 400 Bad Request

**Database Testing:**
```python
def test_save_experiment_metadata(self, temp_database, sample_metadata):
    """Test saving experiment metadata to database."""
    db_manager = temp_database

    metadata_dict = sample_metadata.to_dict()
    metadata_dict['filename'] = 'test_file.csv'
    experiment_id = db_manager.save_experiment_metadata(metadata=metadata_dict)

    assert experiment_id is not None
    assert isinstance(experiment_id, int)
```
See `tests/test_database.py:48-58`

- Use `temp_database` fixture for isolation
- Convert dataclasses to dict for database operations
- Assert returned IDs and data integrity
- Each test gets fresh database via fixture cleanup

**Content Assertion:**
- Check for expected HTML elements: `assert b'Preview top 5 rows' in response.data`
- Check for data in response: `assert 'TEST001' in html_content`
- Use string/bytes conversion appropriately for Flask response data

## Test Markers

Pytest markers defined in `pytest.ini:22-27`:

```
markers =
    unit: Unit tests for individual components
    integration: Integration tests for workflows
    slow: Tests that take a long time to run
    database: Tests that require database access
    web: Tests that require web server
```

**Usage:**
```python
@pytest.mark.unit
class TestColumnDetector:
    ...

@pytest.mark.integration
@pytest.mark.web
class TestFlaskRoutes:
    ...

@pytest.mark.unit
@pytest.mark.database
class TestDatabaseOperations:
    ...
```

## Test Execution Strategy

**Fast feedback (unit only):**
```bash
pytest -m unit
```

**Complete pre-commit (unit + integration):**
```bash
pytest -m "unit or integration" --cov
```

**Subset by module:**
```bash
pytest tests/test_enrichment_service.py -v
```

**With markers:**
```bash
pytest -m "not slow" --cov  # Skip slow tests
pytest -m "database"         # Only database tests
```

---

*Testing analysis: 2026-02-25*
