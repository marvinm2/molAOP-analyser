"""Tests for the one-click PXR batch demo and the configurable status redirect."""
import json
from datetime import datetime, timedelta, timezone

import pytest

import app as app_module
from database import BatchRecord, ConditionRecord


def _expiry():
    return datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(days=14)


@pytest.fixture
def demo_client(flask_client, temp_database, tmp_path, monkeypatch):
    """Flask client wired to a temp DB with run_batch / reference loading stubbed.

    The demo route copies the real shipped PXR TSVs and harmonises them (fast
    enough), but we stub the background thread and reference-set loading so the
    test stays deterministic and offline.
    """
    monkeypatch.setattr(app_module, 'db_manager', temp_database)
    monkeypatch.setattr(app_module, 'run_batch', lambda *a, **k: None)
    monkeypatch.setattr(app_module, 'load_cached_reference_sets', lambda *a, **k: ({}, 'mock'))
    monkeypatch.setattr('config.Config.UPLOAD_FOLDER', str(tmp_path))
    return flask_client


class TestBatchDemo:
    def test_creates_batch_and_redirects_to_compare(self, demo_client, temp_database):
        r = demo_client.post('/batch/demo')
        assert r.status_code == 302
        loc = r.headers['Location']
        assert '/batch/' in loc and loc.endswith('/progress?next=compare')

        uuid = loc.split('/batch/')[1].split('/progress')[0]
        session = temp_database.get_session()
        try:
            batch = session.query(BatchRecord).filter_by(uuid=uuid).first()
            assert batch is not None
            assert batch.aop_id == 'AOP:DEMO'
            assert batch.id_column == 'GENE_SYMBOL'
            assert batch.fc_column == 'logFC'
            assert batch.pval_column == 'adj.P.Val'
            assert batch.harmonised_gene_count > 0
            conditions = session.query(ConditionRecord).filter_by(batch_id=batch.id).order_by(ConditionRecord.position).all()
            assert [c.condition_label for c in conditions] == [
                'PXR agonist 1 (TO90137)', 'PXR agonist 2 (SR12813)'
            ]
        finally:
            session.close()

    def test_progress_page_renders(self, demo_client, temp_database):
        # Seed a running batch and fetch its progress page.
        session = temp_database.get_session()
        try:
            batch = BatchRecord(uuid='33333333-3333-3333-3333-333333333333', status='running',
                                aop_id='AOP:DEMO', batch_name='Demo', expires_at=_expiry())
            session.add(batch)
            session.commit()
        finally:
            session.close()
        r = demo_client.get('/batch/33333333-3333-3333-3333-333333333333/progress?next=compare')
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert 'batch-progress-container' in body
        assert 'next=compare' in body


class TestStatusRedirect:
    """The ?next param flips the completion redirect (default = summary)."""

    def _seed_complete(self, db_manager):
        session = db_manager.get_session()
        try:
            batch = BatchRecord(uuid='44444444-4444-4444-4444-444444444444', status='complete',
                                aop_id='AOP:DEMO', batch_name='Done', expires_at=_expiry())
            session.add(batch)
            session.commit()
        finally:
            session.close()
        return '44444444-4444-4444-4444-444444444444'

    def test_next_compare_redirects_to_compare(self, flask_client, temp_database, monkeypatch):
        monkeypatch.setattr(app_module, 'db_manager', temp_database)
        uuid = self._seed_complete(temp_database)
        r = flask_client.get(f'/batch/{uuid}/status?next=compare')
        assert r.headers.get('HX-Redirect') == f'/batch/{uuid}/compare'

    def test_default_redirects_to_summary(self, flask_client, temp_database, monkeypatch):
        monkeypatch.setattr(app_module, 'db_manager', temp_database)
        uuid = self._seed_complete(temp_database)
        r = flask_client.get(f'/batch/{uuid}/status')
        assert r.headers.get('HX-Redirect') == f'/batch/{uuid}/summary'
