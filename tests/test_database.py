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

        # Test to_dict method with database context
        session = db_manager.get_session()
        try:
            session.add(record)
            session.commit()
            assert record.upload_timestamp is not None  # Default set on INSERT
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


class TestMinConfidencePersistence:
    """Issue #60: the mapping-confidence threshold is stored for reproducibility."""

    def test_threshold_persisted_and_read_back(self, temp_database, sample_metadata):
        db_manager = temp_database
        metadata_dict = sample_metadata.to_dict()
        metadata_dict['filename'] = 'conf_test.csv'

        experiment_id = db_manager.save_experiment_metadata(
            metadata=metadata_dict,
            analysis_params={'aop_id': 'AOP:1', 'min_confidence': 'high'},
        )
        assert db_manager.get_experiment(experiment_id)['min_confidence'] == 'high'

    def test_missing_threshold_reads_back_as_all(self, temp_database, sample_metadata):
        """Pre-#60 rows (NULL column) are coerced to the historical behaviour."""
        db_manager = temp_database
        metadata_dict = sample_metadata.to_dict()
        metadata_dict['filename'] = 'legacy.csv'

        experiment_id = db_manager.save_experiment_metadata(
            metadata=metadata_dict,
            analysis_params={'aop_id': 'AOP:1'},
        )
        assert db_manager.get_experiment(experiment_id)['min_confidence'] == 'all'

    def test_migration_adds_the_column_to_legacy_tables(self, tmp_path):
        """The PRAGMA-then-ALTER migration is additive and idempotent."""
        from sqlalchemy import create_engine, text
        from database import _ensure_min_confidence_column

        db_path = tmp_path / 'legacy.db'
        engine = create_engine(f'sqlite:///{db_path}')
        with engine.connect() as conn:
            conn.execute(text('CREATE TABLE experiments (id INTEGER PRIMARY KEY)'))
            conn.execute(text('CREATE TABLE batches (id INTEGER PRIMARY KEY)'))
            conn.commit()

        _ensure_min_confidence_column(engine)
        _ensure_min_confidence_column(engine)  # second run must be a no-op

        for table in ('experiments', 'batches'):
            with engine.connect() as conn:
                cols = [row[1] for row in conn.execute(text(f'PRAGMA table_info({table})'))]
            assert cols.count('min_confidence') == 1

    def test_batch_record_stores_the_threshold(self, temp_database):
        from datetime import datetime, timedelta, timezone
        from database import BatchRecord

        session = temp_database.get_session()
        try:
            batch = BatchRecord(
                uuid='conf-batch-uuid',
                aop_id='AOP:1',
                min_confidence='medium',
                expires_at=datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(days=1),
            )
            session.add(batch)
            session.commit()
            stored = session.query(BatchRecord).filter_by(uuid='conf-batch-uuid').one()
            assert stored.min_confidence == 'medium'
        finally:
            session.close()

class TestDroppedRowsColumn:
    """Issue #103: the discarded symbol-less row count is persisted per condition."""

    def test_migration_adds_the_column_to_a_legacy_table(self, tmp_path):
        """PRAGMA-then-ALTER, additive and idempotent, as with #68/#69."""
        from sqlalchemy import create_engine, text
        from database import _ensure_dropped_rows_column

        db_path = tmp_path / 'legacy.db'
        engine = create_engine(f'sqlite:///{db_path}')
        with engine.connect() as conn:
            conn.execute(text('CREATE TABLE batch_conditions (id INTEGER PRIMARY KEY)'))
            conn.commit()

        _ensure_dropped_rows_column(engine)
        _ensure_dropped_rows_column(engine)  # second run must be a no-op

        with engine.connect() as conn:
            cols = [row[1] for row in conn.execute(text('PRAGMA table_info(batch_conditions)'))]
        assert cols.count('dropped_unidentified_rows') == 1

    def test_condition_defaults_to_null_not_zero(self, temp_database):
        """A condition run before the count existed must read as unrecorded.

        NULL and 0 are different claims: 0 says every row carried an identifier,
        NULL says nobody looked. The summary renders the former and not the
        latter, so the distinction has to survive storage.
        """
        from datetime import datetime, timedelta, timezone
        from database import BatchRecord, ConditionRecord

        session = temp_database.get_session()
        try:
            batch = BatchRecord(
                uuid='dropped-rows-batch',
                aop_id='AOP:1',
                expires_at=datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(days=1),
            )
            session.add(batch)
            session.commit()

            legacy = ConditionRecord(
                batch_id=batch.id, position=0, condition_label='legacy',
            )
            counted = ConditionRecord(
                batch_id=batch.id, position=1, condition_label='counted',
                dropped_unidentified_rows=875,
            )
            session.add_all([legacy, counted])
            session.commit()

            stored = {c.condition_label: c.dropped_unidentified_rows
                      for c in session.query(ConditionRecord).all()}
            assert stored['legacy'] is None
            assert stored['counted'] == 875
        finally:
            session.close()


class TestBackgroundRuleColumns:
    """Issue #132: the background rules are stored, and NULL means legacy.

    Every batch run before these columns existed used the loader's drop-any-NaN
    universe and an intersection across conditions. Reading NULL back as
    anything else would silently reinterpret results that are already published
    in shared links and exported reports.
    """

    def test_null_reads_back_as_the_pre_132_behaviour(self):
        from database import BatchRecord

        batch = BatchRecord()
        assert batch.effective_background_universe() == 'testable'
        assert batch.effective_background_harmonisation() == 'intersection'

    def test_stored_values_are_returned(self):
        from database import BatchRecord

        batch = BatchRecord(
            background_universe='measured', background_harmonisation='union'
        )
        assert batch.effective_background_universe() == 'measured'
        assert batch.effective_background_harmonisation() == 'union'

    def test_description_names_the_rules_not_just_the_size(self):
        from database import BatchRecord

        batch = BatchRecord(
            method='ora', background_universe='measured',
            background_harmonisation='union', harmonised_gene_count=18339,
        )
        text = batch.background_description()
        assert '18,339' in text
        assert 'any condition' in text

        legacy = BatchRecord(
            method='ora', background_universe='testable',
            background_harmonisation='intersection', harmonised_gene_count=4086,
        )
        legacy_text = legacy.background_description()
        assert '4,086' in legacy_text
        assert 'every condition' in legacy_text
        assert 'p-value' in legacy_text

    def test_gsea_has_no_background_to_describe(self):
        """GSEA scores a ranked list; quoting a background size would describe
        a quantity the analysis never used."""
        from database import BatchRecord

        batch = BatchRecord(method='gsea', harmonised_gene_count=4086)
        assert '4,086' not in batch.background_description()
        assert 'GSEA' in batch.background_description()

    def test_per_condition_describes_itself_without_a_size(self):
        from database import BatchRecord

        batch = BatchRecord(
            method='ora', background_harmonisation='per_condition',
            background_universe='measured', harmonised_gene_count=None,
        )
        assert 'per condition' in batch.background_description()

    def test_migration_is_idempotent_and_additive(self, tmp_path):
        """The columns are added to a database created without them."""
        from sqlalchemy import create_engine, text as sql_text
        from database import _ensure_background_rule_columns, Base

        db_path = tmp_path / 'legacy.db'
        engine = create_engine(f'sqlite:///{db_path}')
        Base.metadata.create_all(bind=engine)
        with engine.connect() as conn:
            conn.execute(sql_text('ALTER TABLE batches DROP COLUMN background_universe'))
            conn.execute(sql_text('ALTER TABLE batches DROP COLUMN background_harmonisation'))
            conn.commit()

        _ensure_background_rule_columns(engine)
        _ensure_background_rule_columns(engine)  # no-op on the second pass

        with engine.connect() as conn:
            cols = {row[1] for row in conn.execute(sql_text('PRAGMA table_info(batches)'))}
        assert {'background_universe', 'background_harmonisation'} <= cols
