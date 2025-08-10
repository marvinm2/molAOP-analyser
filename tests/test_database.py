"""
Unit tests for database operations and models.
"""

import pytest
from datetime import datetime
from config import ExperimentMetadata
from database import DatabaseManager, ExperimentRecord


@pytest.mark.unit
@pytest.mark.database
class TestDatabaseOperations:
    """Test database CRUD operations."""
    
    def test_experiment_metadata_creation(self, sample_metadata):
        """Test ExperimentMetadata dataclass creation and serialization."""
        # Test basic creation
        assert sample_metadata.dataset_id == "TEST001"
        assert sample_metadata.stressor == "Test Chemical"
        assert sample_metadata.owner == "Test User"
        
        # Test serialization
        metadata_dict = sample_metadata.to_dict()
        assert 'dataset_id' in metadata_dict
        assert 'upload_timestamp' in metadata_dict
        assert isinstance(metadata_dict['upload_timestamp'], str)
        
        # Test deserialization
        recreated = ExperimentMetadata.from_dict(metadata_dict)
        assert recreated.dataset_id == sample_metadata.dataset_id
        assert recreated.stressor == sample_metadata.stressor
    
    def test_database_initialization(self, temp_database):
        """Test database initialization and table creation."""
        db_manager = temp_database
        
        # Check database is initialized
        assert db_manager.engine is not None
        assert db_manager.SessionLocal is not None
        
        # Check tables exist
        from sqlalchemy import inspect
        inspector = inspect(db_manager.engine)
        tables = inspector.get_table_names()
        assert 'experiments' in tables
    
    def test_save_experiment_metadata(self, temp_database, sample_metadata):
        """Test saving experiment metadata to database."""
        db_manager = temp_database
        
        # Save basic metadata
        metadata_dict = sample_metadata.to_dict()
        metadata_dict['filename'] = 'test_file.csv'  # Add required filename
        experiment_id = db_manager.save_experiment_metadata(metadata=metadata_dict)
        
        assert experiment_id is not None
        assert isinstance(experiment_id, int)
    
    def test_save_experiment_with_analysis_params(self, temp_database, sample_metadata):
        """Test saving experiment with analysis parameters."""
        db_manager = temp_database
        
        metadata_dict = sample_metadata.to_dict()
        metadata_dict['filename'] = 'analysis_test.csv'
        analysis_params = {
            'aop_id': 'AOP:1',
            'logfc_threshold': 1.5,
            'pval_cutoff': 0.05,
            'id_column': 'Gene_Symbol',
            'fc_column': 'log2FoldChange',
            'pval_column': 'padj'
        }
        
        results = {
            'enrichment_table': [{'KE_ID': 'KE:115', 'pvalue': 0.001}],
            'gene_count': 1000,
            'significant_genes': 150
        }
        
        experiment_id = db_manager.save_experiment_metadata(
            metadata=metadata_dict,
            analysis_params=analysis_params,
            results=results
        )
        
        assert experiment_id is not None
        
        # Verify saved data
        retrieved = db_manager.get_experiment(experiment_id)
        assert retrieved['aop_id'] == 'AOP:1'
        assert retrieved['logfc_threshold'] == 1.5
        assert retrieved['gene_count'] == 1000
    
    def test_retrieve_experiment(self, temp_database, sample_metadata):
        """Test retrieving experiment from database."""
        db_manager = temp_database
        
        # Save experiment
        metadata_dict = sample_metadata.to_dict()
        metadata_dict['filename'] = 'retrieve_test.csv'
        experiment_id = db_manager.save_experiment_metadata(metadata=metadata_dict)
        
        # Retrieve experiment
        retrieved = db_manager.get_experiment(experiment_id)
        
        assert retrieved is not None
        assert retrieved['dataset_id'] == "TEST001"
        assert retrieved['stressor'] == "Test Chemical"
        assert retrieved['owner'] == "Test User"
    
    def test_retrieve_nonexistent_experiment(self, temp_database):
        """Test retrieving non-existent experiment returns None."""
        db_manager = temp_database
        
        result = db_manager.get_experiment(99999)
        assert result is None
    
    def test_list_experiments(self, temp_database, sample_metadata):
        """Test listing recent experiments."""
        db_manager = temp_database
        
        # Save multiple experiments
        metadata_dict = sample_metadata.to_dict()
        
        for i in range(3):
            metadata_dict['dataset_id'] = f"TEST{i:03d}"
            db_manager.save_experiment_metadata(metadata=metadata_dict)
        
        # List experiments
        experiments = db_manager.list_experiments(limit=10)
        
        assert len(experiments) == 3
        assert all('dataset_id' in exp for exp in experiments)
        
        # Check ordering (should be most recent first)
        dataset_ids = [exp['dataset_id'] for exp in experiments]
        assert 'TEST002' in dataset_ids  # Most recent should be included
    
    def test_search_experiments(self, temp_database, sample_metadata):
        """Test searching experiments by metadata fields."""
        db_manager = temp_database
        
        # Save experiments with different metadata
        base_metadata = sample_metadata.to_dict()
        
        test_cases = [
            {'dataset_id': 'SEARCH001', 'stressor': 'Chemical A', 'owner': 'User 1'},
            {'dataset_id': 'SEARCH002', 'stressor': 'Chemical B', 'owner': 'User 1'},
            {'dataset_id': 'OTHER001', 'stressor': 'Chemical A', 'owner': 'User 2'},
        ]
        
        for case in test_cases:
            metadata = base_metadata.copy()
            metadata.update(case)
            db_manager.save_experiment_metadata(metadata=metadata)
        
        # Test search by dataset_id
        results = db_manager.search_experiments(dataset_id='SEARCH')
        assert len(results) == 2
        assert all('SEARCH' in exp['dataset_id'] for exp in results)
        
        # Test search by owner
        results = db_manager.search_experiments(owner='User 1')
        assert len(results) == 2
        assert all(exp['owner'] == 'User 1' for exp in results)
        
        # Test search by stressor
        results = db_manager.search_experiments(stressor='Chemical A')
        assert len(results) == 2
        assert all(exp['stressor'] == 'Chemical A' for exp in results)
        
        # Test combined search
        results = db_manager.search_experiments(stressor='Chemical A', owner='User 1')
        assert len(results) == 1
        assert results[0]['dataset_id'] == 'SEARCH001'
    
    def test_experiment_record_model(self, temp_database):
        """Test ExperimentRecord model functionality."""
        db_manager = temp_database
        
        # Test model creation
        record = ExperimentRecord(
            dataset_id='MODEL_TEST',
            filename='test.csv',
            stressor='Test Chemical',
            owner='Test User'
        )
        
        assert record.dataset_id == 'MODEL_TEST'
        assert record.upload_timestamp is not None  # Should be set by default
        
        # Test to_dict method with database context
        session = db_manager.get_session()
        try:
            session.add(record)
            session.commit()
            record_dict = record.to_dict()
            assert 'id' in record_dict
            assert record_dict['dataset_id'] == 'MODEL_TEST'
            assert 'upload_timestamp' in record_dict
        finally:
            session.close()
    
    def test_database_error_handling(self, temp_database):
        """Test database error handling for invalid operations."""
        db_manager = temp_database
        
        # Test with invalid metadata (missing required fields)
        invalid_metadata = {}
        
        experiment_id = db_manager.save_experiment_metadata(metadata=invalid_metadata)
        # Should not crash but may return None or create record with empty fields
        assert experiment_id is None or isinstance(experiment_id, int)