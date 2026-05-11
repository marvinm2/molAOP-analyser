"""
Tests for the SharedResult model and shareable URL endpoints.

Covers:
- SharedResult.create() factory method
- cleanup_expired_shared_results() function
- POST /api/share endpoint
- GET /results/<uuid> route (valid, not-found, expired)
"""

import json
import uuid
import tempfile
import os
import pytest
from datetime import datetime, timedelta, timezone

from database import (
    DatabaseManager,
    SharedResult,
    Base,
    cleanup_expired_shared_results,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def temp_db():
    """Provide an isolated SQLite database with all tables created."""
    with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
        db_path = f.name

    manager = DatabaseManager(db_url=f"sqlite:///{db_path}")
    manager.initialize()

    yield manager

    # Teardown
    manager.get_session().close()
    try:
        os.unlink(db_path)
    except OSError:
        pass


@pytest.fixture
def db_session(temp_db):
    """Provide a fresh session bound to the isolated temp_db."""
    session = temp_db.get_session()
    yield session
    session.close()


@pytest.fixture
def sample_shared_result():
    """Return a SharedResult instance (not persisted) with valid data."""
    return SharedResult.create(
        enrichment_table=[{'KE': 'KE:1', 'p_value': 0.01, 'FDR': 0.05}],
        network={'nodes': [{'data': {'id': 'KE:1'}}], 'edges': []},
        ke_gene_map={'KE:1': [{'id': 'GENE1', 'log2FC': 1.5, 'significant': True}]},
        ke_type_map={'KE:1': 'MIE'},
        ke_title_map={'KE:1': 'Test Key Event'},
        aop_id='AOP:33',
        aop_label='Kidney Failure AOP',
        metadata={
            'dataset_id': 'CISPLATIN_24H',
            'gene_count': 5000,
            'significant_genes': 200,
        },
    )


# ---------------------------------------------------------------------------
# Unit tests: SharedResult model
# ---------------------------------------------------------------------------

class TestSharedResultModel:
    """Unit tests for the SharedResult SQLAlchemy model."""

    def test_shared_result_create(self, sample_shared_result):
        """SharedResult.create() returns instance with valid UUID and 30-day expiry."""
        record = sample_shared_result

        # UUID is 36 chars (xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx)
        assert len(record.uuid) == 36
        uuid.UUID(record.uuid)  # Raises ValueError if not valid UUID

        # Expires approximately 30 days from now
        delta = record.expires_at - datetime.now(timezone.utc).replace(tzinfo=None)
        assert 29 <= delta.days <= 30, f"Expected ~30 day expiry, got {delta.days} days"

        # Enrichment JSON is valid
        parsed = json.loads(record.enrichment_json)
        assert isinstance(parsed, list)
        assert parsed[0]['KE'] == 'KE:1'

        # Network JSON is valid
        net = json.loads(record.network_json)
        assert 'nodes' in net
        assert 'edges' in net

        # Metadata stored
        assert record.aop_id == 'AOP:33'
        assert record.aop_label == 'Kidney Failure AOP'
        assert record.dataset_label == 'CISPLATIN_24H'
        assert record.gene_count == 5000
        assert record.significant_genes == 200

    def test_shared_result_expiry_in_past(self):
        """SharedResult with past expires_at is treated as expired."""
        record = SharedResult.create(
            enrichment_table=[],
            network={'nodes': [], 'edges': []},
            ke_gene_map=None,
            ke_type_map=None,
            ke_title_map=None,
            aop_id='AOP:1',
            aop_label='Test AOP',
            metadata={},
        )
        # Override expires_at to the past
        record.expires_at = datetime.now(timezone.utc).replace(tzinfo=None) - timedelta(hours=1)

        assert record.expires_at < datetime.now(timezone.utc).replace(tzinfo=None), "Record should be expired"

    def test_cleanup_expired(self, db_session):
        """cleanup_expired_shared_results deletes expired rows and keeps valid ones."""
        # Create one expired and one valid record
        expired = SharedResult.create(
            enrichment_table=[],
            network={'nodes': [], 'edges': []},
            ke_gene_map=None, ke_type_map=None, ke_title_map=None,
            aop_id='AOP:1', aop_label='Expired AOP', metadata={},
        )
        expired.expires_at = datetime.now(timezone.utc).replace(tzinfo=None) - timedelta(days=1)

        valid = SharedResult.create(
            enrichment_table=[],
            network={'nodes': [], 'edges': []},
            ke_gene_map=None, ke_type_map=None, ke_title_map=None,
            aop_id='AOP:2', aop_label='Valid AOP', metadata={},
        )

        db_session.add(expired)
        db_session.add(valid)
        db_session.commit()

        expired_uuid = expired.uuid
        valid_uuid = valid.uuid

        # Run cleanup
        cleanup_expired_shared_results(db_session)

        # Expired row should be gone
        assert db_session.query(SharedResult).filter_by(uuid=expired_uuid).first() is None
        # Valid row should remain
        assert db_session.query(SharedResult).filter_by(uuid=valid_uuid).first() is not None


# ---------------------------------------------------------------------------
# Integration tests: Flask routes
# ---------------------------------------------------------------------------

class TestShareEndpoints:
    """Integration tests for POST /api/share and GET /results/<uuid> routes."""

    @pytest.fixture(autouse=True)
    def init_db(self, temp_db, monkeypatch):
        """Patch the app's db_manager with the isolated temp_db."""
        import app as flask_app_module
        monkeypatch.setattr(flask_app_module, 'db_manager', temp_db)

    def _share_payload(self):
        """Return a minimal valid payload for POST /api/share."""
        return {
            'enrichment_table': [{'KE': 'KE:1', 'p_value': 0.01}],
            'network': {'nodes': [], 'edges': []},
            'ke_gene_map': {},
            'ke_type_map': {'KE:1': 'MIE'},
            'ke_title_map': {'KE:1': 'Test KE'},
            'aop_id': 'AOP:33',
            'aop_label': 'Test AOP',
            'metadata': {
                'dataset_id': 'TEST001',
                'gene_count': 1000,
                'significant_genes': 50,
            },
        }

    def test_api_share_endpoint(self, flask_client):
        """POST /api/share with valid JSON returns 200 with a UUID."""
        payload = self._share_payload()
        resp = flask_client.post(
            '/api/share',
            data=json.dumps(payload),
            content_type='application/json',
        )
        assert resp.status_code == 200, f"Expected 200, got {resp.status_code}: {resp.data}"
        data = resp.get_json()
        assert 'uuid' in data
        assert len(data['uuid']) == 36
        uuid.UUID(data['uuid'])  # Validate UUID format

    def test_api_share_empty_body(self, flask_client):
        """POST /api/share with no body returns 400."""
        resp = flask_client.post(
            '/api/share',
            data='',
            content_type='application/json',
        )
        assert resp.status_code == 400

    def test_shared_results_route_valid(self, flask_client):
        """GET /results/<uuid> returns 200 for a valid, non-expired shared result."""
        # Create a share first
        payload = self._share_payload()
        create_resp = flask_client.post(
            '/api/share',
            data=json.dumps(payload),
            content_type='application/json',
        )
        assert create_resp.status_code == 200
        shared_uuid = create_resp.get_json()['uuid']

        # Retrieve it
        resp = flask_client.get(f'/results/{shared_uuid}')
        assert resp.status_code == 200
        html = resp.data.decode('utf-8')
        # Should contain enrichment table and network elements
        assert 'resultsTable' in html or 'KE:1' in html or 'network' in html.lower()
        # Should NOT contain volcano plot
        assert 'volcano' not in html.lower()
        # Should contain shared banner
        assert 'Shared Analysis' in html or 'shared' in html.lower()

    def test_shared_results_route_not_found(self, flask_client):
        """GET /results/<random-uuid> returns 404."""
        random_uuid = str(uuid.uuid4())
        resp = flask_client.get(f'/results/{random_uuid}')
        assert resp.status_code == 404
        html = resp.data.decode('utf-8')
        assert 'not found' in html.lower()

    def test_shared_results_route_expired(self, flask_client, temp_db):
        """GET /results/<uuid> for an expired record returns 410 with expiry message."""
        # Insert an expired record directly via the DB
        session = temp_db.get_session()
        try:
            record = SharedResult.create(
                enrichment_table=[],
                network={'nodes': [], 'edges': []},
                ke_gene_map=None, ke_type_map=None, ke_title_map=None,
                aop_id='AOP:99', aop_label='Expired AOP', metadata={},
            )
            record.expires_at = datetime.now(timezone.utc).replace(tzinfo=None) - timedelta(hours=1)
            session.add(record)
            session.commit()
            expired_uuid = record.uuid
        finally:
            session.close()

        resp = flask_client.get(f'/results/{expired_uuid}')
        assert resp.status_code == 410, f"Expected 410, got {resp.status_code}"
        html = resp.data.decode('utf-8')
        assert 'expired' in html.lower()
        assert 'Contact the owner' in html or 'expired' in html.lower()
