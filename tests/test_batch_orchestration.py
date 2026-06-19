"""End-to-end tests for the background batch orchestration (run_batch).

Covers (plan T-1/T-2):
- run_batch processes every condition, persists results, marks batch complete
- cancellation stops the run; a failing condition marks the batch failed
"""
import json
import os
from datetime import datetime, timedelta, timezone

import pandas as pd
import pytest

from database import BatchRecord, ConditionRecord
import services.batch_service as batch_service
from services.batch_service import run_batch


# Reference sets with >=5 genes per KE so Fisher's test is actually computed.
REFERENCE_SETS = {
    'KE:115': {f'G{i}' for i in range(1, 11)},    # G1..G10
    'KE:116': {f'G{i}' for i in range(11, 21)},   # G11..G20
}


def _write_dataset(path):
    """Write a small TSV: G1..G10 significant, G11..G20 not."""
    rows = []
    for i in range(1, 21):
        sig = i <= 10
        rows.append({
            'gene': f'G{i}',
            'logFC': 2.0 if sig else 0.1,
            'pval': 0.001 if sig else 0.5,
        })
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)


def _seed_batch(db_manager, upload_root, status='pending', n_conditions=2):
    """Create a batch + conditions in the DB and write their input files."""
    import uuid as _uuid
    batch_uuid = str(_uuid.uuid4())
    batch_dir = os.path.join(upload_root, batch_uuid)
    os.makedirs(batch_dir, exist_ok=True)

    harmonised = sorted({f'G{i}' for i in range(1, 21)})
    session = db_manager.get_session()
    try:
        batch = BatchRecord(
            uuid=batch_uuid, status=status, aop_id='AOP:1', aop_label='Test AOP',
            logfc_threshold=1.0, pval_cutoff=0.05, selected_resources='WikiPathways',
            id_column='gene', fc_column='logFC', pval_column='pval',
            harmonised_background=json.dumps(harmonised), harmonised_gene_count=len(harmonised),
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


class TestRunBatch:
    def test_full_run_marks_complete_and_persists(self, temp_database, aop_stub, upload_root):
        batch_id, _ = _seed_batch(temp_database, upload_root, status='pending')

        run_batch(batch_id, temp_database.db_url, REFERENCE_SETS)

        session = temp_database.get_session()
        try:
            batch = session.query(BatchRecord).filter_by(id=batch_id).first()
            assert batch.status == 'complete'
            assert batch.completed_at is not None
            conditions = session.query(ConditionRecord).filter_by(batch_id=batch_id).all()
            assert len(conditions) == 2
            for cond in conditions:
                assert cond.status == 'complete'
                assert cond.enrichment_json is not None
                assert cond.network_json is not None
                assert cond.ke_gene_json is not None
                assert cond.gene_count == 20
                assert cond.significant_genes == 10
                # KE:115 (G1..G10 all significant) must be enriched.
                enrichment = json.loads(cond.enrichment_json)
                assert any(r['KE'] == 'KE:115' for r in enrichment)
        finally:
            session.close()

    def test_cancellation_stops_run(self, temp_database, aop_stub, upload_root, monkeypatch):
        """Cooperative cancel arriving after condition 0 stops before condition 1."""
        batch_id, _ = _seed_batch(temp_database, upload_root, status='pending')

        real_run_condition = batch_service._run_condition

        def cancelling_run_condition(cond, batch, *args):
            db_session = args[-1]
            real_run_condition(cond, batch, *args)
            if cond.position == 0:
                # Simulate a cancel request landing mid-run.
                batch.status = 'cancelled'
                db_session.commit()
        monkeypatch.setattr(batch_service, '_run_condition', cancelling_run_condition)

        run_batch(batch_id, temp_database.db_url, REFERENCE_SETS)

        session = temp_database.get_session()
        try:
            batch = session.query(BatchRecord).filter_by(id=batch_id).first()
            assert batch.status == 'cancelled'
            conditions = sorted(
                session.query(ConditionRecord).filter_by(batch_id=batch_id).all(),
                key=lambda c: c.position,
            )
            assert conditions[0].status == 'complete'  # first ran
            assert conditions[1].status == 'pending'    # second skipped by cancel
        finally:
            session.close()

    def test_failing_condition_marks_failed(self, temp_database, aop_stub, upload_root, monkeypatch):
        batch_id, _ = _seed_batch(temp_database, upload_root, status='pending')

        def _boom(*args, **kwargs):
            raise RuntimeError('enrichment exploded')
        monkeypatch.setattr('services.enrichment_service.run_enrichment_analysis', _boom)

        run_batch(batch_id, temp_database.db_url, REFERENCE_SETS)

        session = temp_database.get_session()
        try:
            batch = session.query(BatchRecord).filter_by(id=batch_id).first()
            assert batch.status == 'failed'
            failed = session.query(ConditionRecord).filter_by(batch_id=batch_id, status='failed').all()
            assert len(failed) >= 1
            assert 'exploded' in (failed[0].error_message or '')
        finally:
            session.close()
