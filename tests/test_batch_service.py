"""Tests for the batch analysis service."""
import os
import tempfile
from unittest.mock import patch

import pytest
import pandas as pd
from services.batch_service import (
    parse_cisplatin_filename,
    create_batch_upload_dir,
    cleanup_batch_upload_dir,
    validate_batch_columns,
    harmonise_backgrounds,
)


class TestParseCisplatinFilename:
    """Tests for cisplatin filename parsing."""

    def test_standard_filename(self):
        result = parse_cisplatin_filename('CSP_4hr_0.1uM.csv')
        assert result['timepoint'] == '4hr'
        assert result['dose'] == '0.1uM'
        assert result['condition_label'] == 'Cisplatin 4hr 0.1uM'

    def test_long_timepoint(self):
        result = parse_cisplatin_filename('CSP_72hr_50uM.csv')
        assert result['timepoint'] == '72hr'
        assert result['dose'] == '50uM'

    def test_non_matching_filename(self):
        result = parse_cisplatin_filename('random_file.csv')
        assert result['condition_label'] == 'random_file.csv'
        assert result['dose'] == ''
        assert result['timepoint'] == ''

    def test_partial_match_no_extension(self):
        result = parse_cisplatin_filename('CSP_4hr_0.1uM')
        assert result['dose'] == ''  # no .csv extension → no match

    def test_path_prefix_stripped(self):
        """Regex searches inside string, so path prefix is OK."""
        result = parse_cisplatin_filename('Cisplatin_Kidney/CSP_8hr_5uM.csv')
        assert result['timepoint'] == '8hr'
        assert result['dose'] == '5uM'


class TestBatchUploadDir:
    """Tests for batch upload directory management."""

    def test_create_and_cleanup(self):
        with tempfile.TemporaryDirectory() as tmp:
            batch_uuid, batch_dir = create_batch_upload_dir(tmp)
            assert os.path.isdir(batch_dir)
            assert batch_uuid in batch_dir

            cleanup_batch_upload_dir(tmp, batch_uuid)
            assert not os.path.isdir(batch_dir)

    def test_cleanup_nonexistent_is_noop(self):
        with tempfile.TemporaryDirectory() as tmp:
            cleanup_batch_upload_dir(tmp, 'nonexistent-uuid')  # should not raise


class TestValidateBatchColumns:
    """Tests for batch column validation."""

    def _write_csv(self, path, columns, rows=3):
        df = pd.DataFrame({c: range(rows) for c in columns})
        df.to_csv(path, index=False)
        return path

    def test_valid_batch(self, tmp_path):
        f1 = self._write_csv(str(tmp_path / 'a.csv'), ['Gene', 'logFC', 'pval'])
        f2 = self._write_csv(str(tmp_path / 'b.csv'), ['Gene', 'logFC', 'pval'])
        infos = [
            {'filepath': f1, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
            {'filepath': f2, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
        ]
        ok, err = validate_batch_columns(infos)
        assert ok
        assert err == ''

    def test_missing_column(self, tmp_path):
        f1 = self._write_csv(str(tmp_path / 'a.csv'), ['Gene', 'logFC'])  # missing pval
        infos = [{'filepath': f1, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'}]
        ok, err = validate_batch_columns(infos)
        assert not ok
        assert 'missing columns' in err.lower()

    def test_column_mismatch_across_files(self, tmp_path):
        f1 = self._write_csv(str(tmp_path / 'a.csv'), ['Gene', 'logFC', 'pval'])
        f2 = self._write_csv(str(tmp_path / 'b.csv'), ['Gene', 'logFC', 'pval'])
        infos = [
            {'filepath': f1, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
            {'filepath': f2, 'id_col': 'Gene', 'fc_col': 'FC', 'pval_col': 'pval'},  # mismatch
        ]
        ok, err = validate_batch_columns(infos)
        assert not ok
        assert 'mismatch' in err.lower()

    def test_empty_list(self):
        ok, err = validate_batch_columns([])
        assert not ok


class TestHarmoniseBackgrounds:
    """Tests for gene universe intersection across batch files."""

    def _write_gene_csv(self, path, genes, fc_vals=None, pvals=None):
        n = len(genes)
        df = pd.DataFrame({
            'Gene': genes,
            'logFC': fc_vals or [1.0] * n,
            'pval': pvals or [0.05] * n,
        })
        df.to_csv(path, index=False)
        return path

    def test_intersection_of_two_files(self, tmp_path):
        f1 = self._write_gene_csv(str(tmp_path / 'a.csv'), ['BRCA1', 'TP53', 'EGFR'])
        f2 = self._write_gene_csv(str(tmp_path / 'b.csv'), ['BRCA1', 'TP53', 'MYC'])
        infos = [
            {'filepath': f1, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
            {'filepath': f2, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
        ]
        harmonised, per_file = harmonise_backgrounds(infos)
        assert harmonised == {'BRCA1', 'TP53'}
        assert per_file[f1] == 3
        assert per_file[f2] == 3

    def test_no_overlap(self, tmp_path):
        f1 = self._write_gene_csv(str(tmp_path / 'a.csv'), ['BRCA1'])
        f2 = self._write_gene_csv(str(tmp_path / 'b.csv'), ['MYC'])
        infos = [
            {'filepath': f1, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
            {'filepath': f2, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
        ]
        harmonised, _ = harmonise_backgrounds(infos)
        assert harmonised == set()


class TestBatchMinConfidencePersistence:
    """Issue #60: _persist_and_launch_batch stores and applies the threshold."""

    @pytest.fixture
    def batch_app(self, temp_database, monkeypatch):
        """app module wired to a temp DB with the background run stubbed out."""
        import app as app_module

        monkeypatch.setattr(app_module, 'db_manager', temp_database)
        monkeypatch.setattr(app_module, 'run_batch', lambda *a, **k: None)
        return app_module

    @staticmethod
    def _launch(app_module, **overrides):
        kwargs = dict(
            batch_uuid='conf-parity-uuid',
            filenames=['a.csv'],
            condition_labels=['A'],
            doses=[''],
            timepoints=[''],
            id_col='Gene',
            fc_col='logFC',
            pval_col='pval',
            aop_id='AOP:1',
            logfc_threshold=0.0,
            pval_threshold=0.05,
            resources=['WikiPathways'],
            harmonised_genes={'BRCA1'},
            batch_name='conf batch',
            owner='',
            description='',
        )
        kwargs.update(overrides)
        return app_module._persist_and_launch_batch(**kwargs)

    def test_threshold_stored_on_the_batch_record(self, batch_app, temp_database):
        from database import BatchRecord

        with patch.object(batch_app, 'load_cached_reference_sets', return_value=({}, 'mock', [])):
            self._launch(batch_app, min_confidence='high')

        session = temp_database.get_session()
        try:
            batch = session.query(BatchRecord).filter_by(uuid='conf-parity-uuid').one()
            assert batch.min_confidence == 'high'
        finally:
            session.close()

    def test_threshold_applied_when_loading_the_reference_sets(self, batch_app):
        """Batch parity: the same shared loader call the single flow makes."""
        with patch.object(batch_app, 'load_cached_reference_sets', return_value=({}, 'mock', [])) as loader:
            self._launch(batch_app, min_confidence='medium')
        assert loader.call_args.kwargs['min_confidence'] == 'medium'

    def test_default_threshold_is_all(self, batch_app, temp_database):
        from database import BatchRecord

        with patch.object(batch_app, 'load_cached_reference_sets', return_value=({}, 'mock', [])) as loader:
            self._launch(batch_app)
        assert loader.call_args.kwargs['min_confidence'] == 'all'

        session = temp_database.get_session()
        try:
            batch = session.query(BatchRecord).filter_by(uuid='conf-parity-uuid').one()
            assert batch.min_confidence == 'all'
        finally:
            session.close()


class TestConditionDroppedRows:
    """Issue #103: each condition records the rows its file lost to a missing symbol."""

    def _write_file(self, path, symbol_less=2):
        """Write a small TSV where `symbol_less` rows carry no gene identifier."""
        rows = [{'Gene': f'G{i}', 'logFC': 2.0, 'pval': 0.001} for i in range(1, 6)]
        for _ in range(symbol_less):
            rows.append({'Gene': None, 'logFC': 1.0, 'pval': 0.01})
        pd.DataFrame(rows).to_csv(path, sep='\t', index=False)

    def _run(self, tmp_path, monkeypatch, temp_database, symbol_less):
        """Run one condition end-to-end and return its stored ConditionRecord."""
        import json
        from datetime import datetime, timedelta, timezone

        from config import Config
        from database import BatchRecord, ConditionRecord
        from services import batch_service

        upload_root = tmp_path / 'uploads'
        (upload_root / 'batch-103').mkdir(parents=True)
        monkeypatch.setattr(Config, 'UPLOAD_FOLDER', str(upload_root))
        self._write_file(upload_root / 'batch-103' / 'cond.tsv', symbol_less)

        session = temp_database.get_session()
        try:
            batch = BatchRecord(
                uuid='batch-103', status='running', aop_id='AOP:1',
                logfc_threshold=1.0, pval_cutoff=0.05,
                id_column='Gene', fc_column='logFC', pval_column='pval',
                expires_at=datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(days=1),
            )
            session.add(batch)
            session.flush()
            cond = ConditionRecord(
                batch_id=batch.id, position=0, filename='cond.tsv',
                condition_label='Cond', status='pending',
            )
            session.add(cond)
            session.commit()

            batch_service._run_condition(
                cond, batch,
                harmonised_genes={f'G{i}' for i in range(1, 6)},
                reference_sets={'KE:115': {f'G{i}' for i in range(1, 6)}},
                ke_list={'KE:115'},
                edges=pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID']),
                ke_type_map={'KE:115': 'MIE'},
                ke_title_map={'KE:115': 'Event 115'},
                db_session=session,
            )
            session.refresh(cond)
            return cond.status, cond.dropped_unidentified_rows
        finally:
            session.close()

    def test_discarded_rows_are_persisted(self, tmp_path, monkeypatch, temp_database):
        """The loader's count reaches the record, so the summary can show it."""
        status, dropped = self._run(tmp_path, monkeypatch, temp_database, symbol_less=2)

        assert status == 'complete'
        assert dropped == 2

    def test_clean_file_records_zero(self, tmp_path, monkeypatch, temp_database):
        """A recorded zero is a real claim — every row had an identifier."""
        status, dropped = self._run(tmp_path, monkeypatch, temp_database, symbol_less=0)

        assert status == 'complete'
        assert dropped == 0
