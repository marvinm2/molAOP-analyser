"""
Test configuration and fixtures for molAOP Analyser test suite.

Provides reusable fixtures for database, Flask app, and test data.
"""

import os
import sys
import tempfile
import pytest
import pandas as pd
from datetime import datetime
from unittest.mock import Mock, patch

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import application components
from app import app as flask_app
from database import DatabaseManager, ExperimentRecord, Base
from config import Config, ExperimentMetadata
from services.report_service import ReportData
from services.column_detector import ColumnDetector
from services.gene_id_validator import GeneIDValidator


# ---------------------------------------------------------------------------
# Issue #104 — keeping the suite off the filesystem and off Chromium
# ---------------------------------------------------------------------------

# The real os.path.exists, captured before any test can patch it over.
_REAL_PATH_EXISTS = os.path.exists

# Extensions of the upload/demo files the route tests pretend are on disk.
_UPLOAD_SUFFIXES = ('.tsv', '.csv', '.txt', '.xlsx')


def uploaded_file_exists(path):
    """os.path.exists replacement that only lies about uploaded data files.

    Route tests need ``/analyze`` to believe the uploaded dataset is on disk,
    and the usual way to arrange that was ``patch('os.path.exists',
    return_value=True)``. That patch is process-wide: matplotlib's search for
    ``matplotlibrc``, kaleido's probe for its Chromium binary and every other
    library's existence check are all told "yes" for the duration of the test.
    ``test_full_workflow_integration`` fails outright when run alone for that
    reason (matplotlib opens the ``matplotlibrc`` it was told exists), and only
    passes in a full run because an earlier test imported matplotlib first —
    exactly the kind of order-dependence behind the intermittent hang in #104.

    This narrows the lie to what the tests actually need: a data file. Anything
    else falls through to the real check.

    Args:
        path: path being probed (str, bytes or os.PathLike).

    Returns:
        bool: True for data-file paths, otherwise the real answer.
    """
    try:
        text = os.fspath(path)
        if isinstance(text, bytes):
            text = text.decode('utf-8', 'replace')
    except TypeError:  # file descriptor or something exotic — defer
        return _REAL_PATH_EXISTS(path)
    return text.lower().endswith(_UPLOAD_SUFFIXES) or _REAL_PATH_EXISTS(path)


@pytest.fixture(autouse=True)
def stub_plotly_image_export(request):
    """Stop any test from spawning kaleido's headless Chromium (#104).

    ``Figure.to_image()`` launches a browser subprocess that can outlive the
    call; when it wedges, pytest sits in a blocking read with no deadline and
    the run stalls for as long as the harness allows. Nothing in the suite
    depends on the *content* of a rendered PNG, so the render is stubbed by
    default and the assertions are about whether a figure was produced at all.

    Opt out with ``@pytest.mark.real_image_render`` for a test that genuinely
    needs kaleido; it then runs against the production timeout in
    ``services.image_render``.
    """
    if 'real_image_render' in request.keywords:
        yield None
        return

    try:
        import plotly.graph_objects as go
    except ImportError:  # pragma: no cover — plotly is a hard dependency
        yield None
        return

    # Smallest valid PNG: a 1x1 transparent pixel. Callers only embed the bytes.
    png_stub = (
        b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01'
        b'\x08\x06\x00\x00\x00\x1f\x15\xc4\x89\x00\x00\x00\nIDATx\x9cc\x00\x01'
        b'\x00\x00\x05\x00\x01\r\n-\xb4\x00\x00\x00\x00IEND\xaeB`\x82'
    )
    with patch.object(go.Figure, 'to_image', autospec=True,
                      return_value=png_stub) as mock_to_image:
        yield mock_to_image


@pytest.fixture(scope='session')
def test_data_dir():
    """Fixture providing path to test data directory."""
    return os.path.join(os.path.dirname(__file__), 'data')


@pytest.fixture
def sample_gene_data():
    """Fixture providing sample gene expression data."""
    return pd.DataFrame({
        'Gene_Symbol': ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS'],
        'log2FoldChange': [2.5, -1.8, 3.2, -2.1, 1.9],
        'padj': [0.001, 0.005, 0.0001, 0.01, 0.03],
        'baseMean': [1000, 800, 1200, 600, 900]
    })


@pytest.fixture
def sample_metadata():
    """Fixture providing sample experiment metadata."""
    return ExperimentMetadata(
        dataset_id="TEST001",
        stressor="Test Chemical",
        dosing="10 µM for 24h",
        owner="Test User",
        description="Sample test experiment"
    )


@pytest.fixture
def sample_report_data(sample_metadata):
    """Fixture providing sample report data."""
    return ReportData(
        metadata=sample_metadata.to_dict(),
        filename="test_data.csv",
        gene_count=1000,
        significant_genes=150,
        aop_id="AOP:1",
        aop_label="Test AOP Pathway",
        logfc_threshold=1.0,
        pval_cutoff=0.05,
        id_column="Gene_Symbol",
        fc_column="log2FoldChange",
        pval_column="padj",
        id_type="HGNC",
        enrichment_results=[
            {
                'KE': 'KE:115',
                'Title': 'Test Key Event',
                'num_overlap': 25,
                'total_KE_genes_in_dataset': 100,
                'p_value': 0.001,
                'FDR': 0.01,
                'odds_ratio': 2.5,
                'pct_sig_in_KE': 25.0,
                'overlap': 'BRCA1, TP53, EGFR',
                'sig_in_KE': 25,
                'sig_not_KE': 125,
                'non_sig_in_KE': 75,
                'non_sig_not_KE': 775,
            }
        ]
    )


@pytest.fixture
def temp_database():
    """Fixture providing a temporary test database."""
    with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
        db_path = f.name
    
    # Create test database manager
    db_manager = DatabaseManager(db_url=f"sqlite:///{db_path}")
    db_manager.initialize()
    
    yield db_manager
    
    # Cleanup
    try:
        os.unlink(db_path)
    except OSError:
        pass


@pytest.fixture
def flask_client():
    """Fixture providing Flask test client."""
    # Configure app for testing
    flask_app.config['TESTING'] = True
    flask_app.config['WTF_CSRF_ENABLED'] = False  # Disable CSRF for tests
    flask_app.config['SECRET_KEY'] = 'test-secret-key'
    
    with flask_app.test_client() as client:
        with flask_app.app_context():
            yield client


@pytest.fixture
def authenticated_client(flask_client):
    """Fixture providing Flask client with session data."""
    with flask_client.session_transaction() as sess:
        sess['experiment_metadata'] = {
            'dataset_id': 'TEST001',
            'stressor': 'Test Chemical',
            'dosing': '10 µM',
            'owner': 'Test User'
        }
    return flask_client


@pytest.fixture
def mock_reference_sets():
    """Fixture providing mock reference gene sets."""
    return {
        'KE:115': {'BRCA1', 'TP53', 'EGFR'},
        'KE:116': {'MYC', 'KRAS', 'PIK3CA'},
        'KE:117': {'PTEN', 'AKT1', 'MTOR'}
    }


@pytest.fixture
def column_detector():
    """Fixture providing column detector instance."""
    return ColumnDetector()


@pytest.fixture
def gene_id_validator():
    """Fixture providing gene ID validator instance."""
    return GeneIDValidator()


@pytest.fixture(autouse=True)
def setup_test_environment(monkeypatch):
    """Automatically set up test environment for all tests."""
    # Mock file validation to avoid requiring actual data files
    monkeypatch.setattr('config.Config.validate_data_files', lambda: True)
    
    # Mock reference data loading
    mock_ref_sets = {
        'KE:115': {'BRCA1', 'TP53', 'EGFR'},
        'KE:116': {'MYC', 'KRAS'}
    }
    monkeypatch.setattr('helpers.load_reference_sets', lambda *args, **kwargs: mock_ref_sets)
    
    # Set test data directory
    test_dir = os.path.dirname(__file__)
    monkeypatch.setenv('TESTING', '1')


@pytest.fixture
def mock_enrichment_results():
    """Fixture providing mock enrichment analysis results."""
    return pd.DataFrame([
        {
            'KE_ID': 'KE:115',
            'KE_title': 'Nuclear receptor activation',
            'overlap_genes': 15,
            'ke_size': 50,
            'pvalue': 0.001,
            'pvalue_adjusted': 0.005,
            'enrichment_ratio': 3.0,
            'overlap_gene_names': 'BRCA1,TP53,EGFR'
        },
        {
            'KE_ID': 'KE:116', 
            'KE_title': 'Cell cycle disruption',
            'overlap_genes': 8,
            'ke_size': 40,
            'pvalue': 0.01,
            'pvalue_adjusted': 0.03,
            'enrichment_ratio': 2.0,
            'overlap_gene_names': 'MYC,KRAS'
        }
    ])


@pytest.fixture
def sample_volcano_data():
    """Fixture providing sample volcano plot data."""
    return [
        {'ID': 'BRCA1', 'log2FC': 2.5, 'pval': 0.001},
        {'ID': 'TP53', 'log2FC': -1.8, 'pval': 0.005},
        {'ID': 'EGFR', 'log2FC': 3.2, 'pval': 0.0001},
        {'ID': 'MYC', 'log2FC': -2.1, 'pval': 0.01},
        {'ID': 'KRAS', 'log2FC': 1.9, 'pval': 0.03}
    ]


# Test data creation helpers
@pytest.fixture
def create_test_csv():
    """Fixture factory for creating temporary CSV test files."""
    def _create_csv(data, filename=None):
        if filename is None:
            fd, filename = tempfile.mkstemp(suffix='.csv')
            os.close(fd)
        
        data.to_csv(filename, index=False)
        return filename
    
    return _create_csv


# Performance test fixtures
@pytest.fixture
def large_gene_dataset():
    """Fixture providing a large dataset for performance testing."""
    import numpy as np
    
    n_genes = 10000
    gene_names = [f"GENE_{i:05d}" for i in range(n_genes)]
    
    return pd.DataFrame({
        'Gene_Symbol': gene_names,
        'log2FoldChange': np.random.normal(0, 1.5, n_genes),
        'padj': np.random.uniform(0, 0.1, n_genes),
        'baseMean': np.random.lognormal(5, 2, n_genes)
    })


# Cleanup fixtures
@pytest.fixture(autouse=True) 
def cleanup_temp_files():
    """Automatically clean up temporary files after tests."""
    temp_files = []
    
    yield temp_files
    
    # Clean up any files that were created during tests
    for filepath in temp_files:
        try:
            if os.path.exists(filepath):
                os.unlink(filepath)
        except OSError:
            pass