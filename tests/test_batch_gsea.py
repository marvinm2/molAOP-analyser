"""Tests for GSEA support in batch mode (issue #76).

Covers:
- the ``batches.method`` column, its idempotent migration and the NULL -> 'ora'
  coercion that keeps pre-#76 batches working
- method dispatch in ``_run_condition`` / ``run_batch``
- NES/FDR shaping in ``build_comparison_matrix`` and the matrix export helper
- the method whitelist on ``POST /batch/analyze``
- backward compatibility: a batch row written without a method still runs,
  compares and exports as ORA
"""
import json
import os
import sqlite3
import tempfile
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace

import pandas as pd
import pytest
from sqlalchemy import text

from database import BatchRecord, ConditionRecord, DatabaseManager, _ensure_method_column
from services.batch_service import run_batch
from services.comparison_service import (
    build_comparison_matrix,
    comparison_matrix_to_dataframe,
)


REFERENCE_SETS = {
    'KE:115': {f'G{i}' for i in range(1, 11)},    # G1..G10
    'KE:116': {f'G{i}' for i in range(11, 21)},   # G11..G20
}


def _write_dataset(path):
    """Write a small TSV: G1..G10 strongly up, G11..G20 flat."""
    rows = []
    for i in range(1, 21):
        sig = i <= 10
        rows.append({
            'gene': f'G{i}',
            'logFC': 2.0 if sig else 0.1,
            'pval': 0.001 if sig else 0.5,
        })
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)


def _seed_batch(db_manager, upload_root, method=None, n_conditions=2):
    """Create a batch (+ conditions and input files) with the given method.

    ``method=None`` writes the column as NULL, i.e. exactly what a batch
    created before issue #76 looks like on disk.
    """
    import uuid as _uuid
    batch_uuid = str(_uuid.uuid4())
    batch_dir = os.path.join(upload_root, batch_uuid)
    os.makedirs(batch_dir, exist_ok=True)

    harmonised = sorted({f'G{i}' for i in range(1, 21)})
    session = db_manager.get_session()
    try:
        batch = BatchRecord(
            uuid=batch_uuid, status='pending', aop_id='AOP:1', aop_label='Test AOP',
            method=method,
            logfc_threshold=1.0, pval_cutoff=0.05, selected_resources='WikiPathways',
            id_column='gene', fc_column='logFC', pval_column='pval',
            harmonised_background=json.dumps(harmonised),
            harmonised_gene_count=len(harmonised),
            batch_name='Test Batch',
            expires_at=datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(days=14),
        )
        session.add(batch)
        session.flush()
        for pos in range(n_conditions):
            fname = f'cond{pos}.tsv'
            _write_dataset(os.path.join(batch_dir, fname))
            session.add(ConditionRecord(
                batch_id=batch.id, position=pos, filename=fname,
                condition_label=f'Cond {pos}', status='pending',
            ))
        session.commit()
        batch_id = batch.id
    finally:
        session.close()
    return batch_id, batch_uuid


@pytest.fixture
def aop_stub(monkeypatch):
    """Stub load_aop_data so run_batch needs no real AOP topology."""
    def _fake_load_aop_data(aop_id):
        ke_list = {'KE:115', 'KE:116'}
        edges = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])
        ke_type_map = {'KE:115': 'MIE', 'KE:116': 'AO'}
        ke_title_map = {'KE:115': 'Event 115', 'KE:116': 'Event 116'}
        return ke_list, edges, ke_type_map, ke_title_map
    monkeypatch.setattr('services.data_service.load_aop_data', _fake_load_aop_data)


@pytest.fixture
def upload_root(tmp_path, monkeypatch):
    """Point Config.UPLOAD_FOLDER at a temp dir for the duration of a test."""
    root = tmp_path / 'uploads'
    root.mkdir()
    monkeypatch.setattr('config.Config.UPLOAD_FOLDER', str(root))
    return str(root)


class TestBatchMethodColumn:
    """Schema, migration and NULL coercion for batches.method."""

    def test_column_exists_on_fresh_database(self, temp_database):
        with temp_database.engine.connect() as conn:
            cols = {row[1] for row in conn.execute(text('PRAGMA table_info(batches)'))}
        assert 'method' in cols

    def test_effective_method_coerces_null_to_ora(self):
        assert BatchRecord(method=None).effective_method() == 'ora'
        assert BatchRecord(method='').effective_method() == 'ora'
        assert BatchRecord(method='gsea').effective_method() == 'gsea'
        assert BatchRecord(method='ora').effective_method() == 'ora'

    def test_migration_adds_column_to_legacy_database(self):
        """A batches table without the column gains it on startup, no data loss."""
        with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
            db_path = f.name
        try:
            # Hand-build a pre-#76 batches table and put a row in it.
            conn = sqlite3.connect(db_path)
            conn.execute(
                'CREATE TABLE batches ('
                'id INTEGER PRIMARY KEY, uuid TEXT, status TEXT, expires_at TIMESTAMP)'
            )
            conn.execute(
                "INSERT INTO batches (uuid, status, expires_at) "
                "VALUES ('legacy-uuid', 'complete', '2999-01-01 00:00:00')"
            )
            conn.commit()
            conn.close()

            manager = DatabaseManager(db_url=f'sqlite:///{db_path}')
            assert manager.initialize() is True

            with manager.engine.connect() as c:
                cols = {row[1] for row in c.execute(text('PRAGMA table_info(batches)'))}
                assert 'method' in cols
                row = c.execute(
                    text("SELECT method FROM batches WHERE uuid = 'legacy-uuid'")
                ).fetchone()
            # The pre-existing row survives; its method is the column default.
            assert row is not None
            assert (row[0] or 'ora') == 'ora'

            # Second call is a no-op rather than an error.
            _ensure_method_column(manager.engine)
        finally:
            try:
                os.unlink(db_path)
            except OSError:
                pass


class TestBatchExecutionDispatch:
    """run_batch honours the batch-level method."""

    def test_gsea_batch_produces_nes_columns(self, temp_database, aop_stub, upload_root):
        batch_id, _ = _seed_batch(temp_database, upload_root, method='gsea')

        run_batch(batch_id, temp_database.db_url, REFERENCE_SETS)

        session = temp_database.get_session()
        try:
            batch = session.query(BatchRecord).filter_by(id=batch_id).first()
            assert batch.status == 'complete'
            conditions = session.query(ConditionRecord).filter_by(batch_id=batch_id).all()
            assert conditions
            for cond in conditions:
                assert cond.status == 'complete'
                enrichment = json.loads(cond.enrichment_json)
                assert enrichment
                # GSEA rows carry NES and lead genes; ORA's odds_ratio is absent.
                for row in enrichment:
                    assert 'NES' in row
                    assert 'lead_genes' in row
                    assert 'odds_ratio' not in row
        finally:
            session.close()

    def test_legacy_batch_without_method_runs_as_ora(self, temp_database, aop_stub, upload_root):
        """Backward compatibility: NULL method must still produce an ORA table."""
        batch_id, _ = _seed_batch(temp_database, upload_root, method=None)

        run_batch(batch_id, temp_database.db_url, REFERENCE_SETS)

        session = temp_database.get_session()
        try:
            batch = session.query(BatchRecord).filter_by(id=batch_id).first()
            assert batch.status == 'complete'
            assert batch.method is None  # not silently backfilled
            conditions = session.query(ConditionRecord).filter_by(batch_id=batch_id).all()
            for cond in conditions:
                enrichment = json.loads(cond.enrichment_json)
                assert enrichment
                for row in enrichment:
                    assert 'odds_ratio' in row
                    assert 'NES' not in row
        finally:
            session.close()

    def test_explicit_ora_batch_matches_legacy_output(self, temp_database, aop_stub, upload_root):
        batch_id, _ = _seed_batch(temp_database, upload_root, method='ora')

        run_batch(batch_id, temp_database.db_url, REFERENCE_SETS)

        session = temp_database.get_session()
        try:
            cond = session.query(ConditionRecord).filter_by(batch_id=batch_id).first()
            enrichment = json.loads(cond.enrichment_json)
            assert any(r['KE'] == 'KE:115' for r in enrichment)
            assert all('odds_ratio' in r for r in enrichment)
        finally:
            session.close()


def _cond(label, entries, dose='', timepoint=''):
    """Build a stand-in ConditionRecord for comparison-matrix tests."""
    return SimpleNamespace(
        condition_label=label,
        enrichment_json=json.dumps(entries),
        ke_gene_json=None,
        gene_count=100,
        significant_genes=10,
        dose=dose,
        timepoint=timepoint,
    )


class TestComparisonMatrixGsea:
    """build_comparison_matrix NES handling."""

    def _gsea_conditions(self):
        c1 = _cond('Low', [
            {'KE': 'KE:1', 'Title': 'DNA damage', 'FDR': 0.01, 'NES': 1.8},
            {'KE': 'KE:2', 'Title': 'Mito dysfunction', 'FDR': 0.40, 'NES': 0.9},
        ])
        c2 = _cond('High', [
            {'KE': 'KE:1', 'Title': 'DNA damage', 'FDR': 0.001, 'NES': 2.4},
            {'KE': 'KE:2', 'Title': 'Mito dysfunction', 'FDR': 0.20, 'NES': 1.2},
        ])
        return [c1, c2]

    def test_gsea_matrix_exposes_nes(self):
        m = build_comparison_matrix(self._gsea_conditions(), method='gsea')

        assert m['method'] == 'gsea'
        assert len(m['nes_matrix']) == len(m['ke_labels'])
        row = m['nes_matrix'][m['ke_labels'].index('KE:1')]
        assert row == [1.8, 2.4]

    def test_nes_kept_for_non_significant_kes(self):
        """A sub-threshold but consistently signed NES is the point of GSEA."""
        m = build_comparison_matrix(self._gsea_conditions(), method='gsea')
        idx = m['ke_labels'].index('KE:2')
        # -log10(FDR) is blanked above the cutoff; NES is not.
        assert m['neg_log10_matrix'][idx] == [None, None]
        assert m['nes_matrix'][idx] == [0.9, 1.2]

    def test_ora_matrix_reports_ora_and_blank_nes(self):
        c1 = _cond('A', [{'KE': 'KE:1', 'Title': 'DNA damage', 'FDR': 0.01}])
        c2 = _cond('B', [{'KE': 'KE:1', 'Title': 'DNA damage', 'FDR': 0.30}])

        m = build_comparison_matrix([c1, c2])

        assert m['method'] == 'ora'
        assert m['nes_matrix'] == [[None, None]]
        # Existing ORA behaviour is untouched.
        assert m['fdr_matrix'] == [[0.01, 0.30]]

    def test_nes_matrix_row_order_matches_ke_labels(self):
        m = build_comparison_matrix(self._gsea_conditions(), method='gsea')
        for idx, ke in enumerate(m['ke_labels']):
            expected = {'KE:1': [1.8, 2.4], 'KE:2': [0.9, 1.2]}[ke]
            assert m['nes_matrix'][idx] == expected


class TestComparisonExportNes:
    """comparison_matrix_to_dataframe 'nes' selector."""

    def test_nes_export_columns(self):
        conds = TestComparisonMatrixGsea()._gsea_conditions()
        m = build_comparison_matrix(conds, method='gsea')

        df = comparison_matrix_to_dataframe(m, which='nes')

        assert list(df.columns) == ['Key Event ID', 'Key Event Title', 'Low', 'High']
        assert df.loc[df['Key Event ID'] == 'KE:1', 'High'].iloc[0] == 2.4

    def test_nes_export_blank_for_ora_matrix(self):
        c1 = _cond('A', [{'KE': 'KE:1', 'Title': 'T', 'FDR': 0.01}])
        m = build_comparison_matrix([c1])

        df = comparison_matrix_to_dataframe(m, which='nes')

        assert df['A'].isna().all()

    def test_unknown_selector_still_rejected(self):
        c1 = _cond('A', [{'KE': 'KE:1', 'Title': 'T', 'FDR': 0.01}])
        m = build_comparison_matrix([c1])
        with pytest.raises(ValueError):
            comparison_matrix_to_dataframe(m, which='bogus')

    def test_matrix_without_nes_key_exports_blank(self):
        """A matrix dict built before #76 has no nes_matrix key at all."""
        c1 = _cond('A', [{'KE': 'KE:1', 'Title': 'T', 'FDR': 0.01}])
        m = build_comparison_matrix([c1])
        m.pop('nes_matrix')

        df = comparison_matrix_to_dataframe(m, which='nes')

        assert df['A'].isna().all()


def _seed_complete_batch(db_manager, method, uuid):
    """Seed a completed batch whose conditions carry method-appropriate rows."""
    session = db_manager.get_session()
    try:
        batch = BatchRecord(
            uuid=uuid, status='complete', aop_id='AOP:1', aop_label='Test AOP',
            method=method, logfc_threshold=1.0, pval_cutoff=0.05,
            selected_resources='WikiPathways', id_column='gene',
            fc_column='logFC', pval_column='pval',
            harmonised_background=json.dumps(['G1', 'G2']), harmonised_gene_count=2,
            batch_name='Compare Test Batch',
            expires_at=datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(days=14),
        )
        session.add(batch)
        session.flush()
        if method == 'gsea':
            enrichment = [{'KE': 'KE:1', 'Title': 'Alpha', 'p_value': 0.0005,
                           'FDR': 0.001, 'NES': 2.1, 'lead_genes': ['TP53'],
                           'total_KE_genes_in_dataset': 20}]
        else:
            enrichment = [{'KE': 'KE:1', 'Title': 'Alpha', 'p_value': 0.0005,
                           'FDR': 0.001, 'num_overlap': 6, 'odds_ratio': 3.0,
                           'total_KE_genes_in_dataset': 20, 'Direction': '6up'}]
        network = {
            'nodes': [{'data': {'id': 'KE:1', 'label': 'Alpha', 'ke_type': 'MIE'}}],
            'edges': [],
        }
        for pos in range(2):
            session.add(ConditionRecord(
                batch_id=batch.id, position=pos, filename=f'c{pos}.tsv',
                condition_label=f'C{pos}', dose=f'{pos + 1}uM', timepoint='4hr',
                status='complete', gene_count=20, significant_genes=6,
                enrichment_json=json.dumps(enrichment),
                network_json=json.dumps(network),
                ke_gene_json=json.dumps({}),
            ))
        session.commit()
    finally:
        session.close()
    return uuid


class TestCompareViewMethod:
    """The compare page reflects the batch's method."""

    def test_gsea_batch_compare_page_shows_nes(self, flask_client, temp_database, monkeypatch):
        import app as app_module
        monkeypatch.setattr(app_module, 'db_manager', temp_database)
        uuid = _seed_complete_batch(
            temp_database, 'gsea', '22222222-2222-2222-2222-222222222222')

        resp = flask_client.get(f'/batch/{uuid}/compare')

        assert resp.status_code == 200
        body = resp.data.decode()
        assert 'GSEA (rank-based)' in body
        assert '"method": "gsea"' in body
        assert '"nes_matrix": [[2.1, 2.1]]' in body

    def test_legacy_batch_compare_page_stays_ora(self, flask_client, temp_database, monkeypatch):
        """A batch with a NULL method still renders the ORA comparison."""
        import app as app_module
        monkeypatch.setattr(app_module, 'db_manager', temp_database)
        uuid = _seed_complete_batch(
            temp_database, None, '33333333-3333-3333-3333-333333333333')

        resp = flask_client.get(f'/batch/{uuid}/compare')

        assert resp.status_code == 200
        body = resp.data.decode()
        assert "Fisher's exact (over-representation)" in body
        assert '"method": "ora"' in body


class TestBatchAnalyzeMethodWhitelist:
    """POST /batch/analyze validates the method field."""

    def test_invalid_method_rejected(self, flask_client):
        resp = flask_client.post('/batch/analyze', data={
            'batch_uuid': 'abc',
            'aop_selection': 'AOP:1',
            'id_col': 'gene', 'fc_col': 'logFC', 'pval_col': 'pval',
            'method': 'bootstrap',
        })
        assert resp.status_code == 400
        assert b'Invalid method' in resp.data

    def test_valid_methods_pass_the_whitelist(self, flask_client):
        """'ora'/'gsea' get past the method check (and fail later, on files)."""
        for method in ('ora', 'gsea', ''):
            resp = flask_client.post('/batch/analyze', data={
                'batch_uuid': 'abc',
                'aop_selection': 'AOP:1',
                'id_col': 'gene', 'fc_col': 'logFC', 'pval_col': 'pval',
                'method': method,
            })
            assert b'Invalid method' not in resp.data
