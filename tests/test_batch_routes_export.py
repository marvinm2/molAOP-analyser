"""Flask route tests for batch export, gene-tracking, and report endpoints.

Covers (plan T-6..T-10):
- /compare/export (csv + xlsx)
- /genes/export (long + summary), /genes/data
- /report (html + pdf)
- guards: non-complete batch -> 409, malformed UUID -> 400
"""
import json
from datetime import datetime, timedelta, timezone

import pytest

import app as app_module
from database import BatchRecord, ConditionRecord


def _expiry():
    return datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(days=14)


_NETWORK = {
    'nodes': [
        {'data': {'id': 'KE:1', 'label': 'Alpha', 'ke_type': 'MIE'}, 'classes': 'significant'},
        {'data': {'id': 'KE:2', 'label': 'Beta', 'ke_type': 'AO'}, 'classes': ''},
    ],
    'edges': [{'data': {'source': 'KE:1', 'target': 'KE:2', 'id': 'KER:1'}}],
}


def _seed_complete_batch(db_manager, status='complete'):
    session = db_manager.get_session()
    try:
        batch = BatchRecord(
            uuid='11111111-1111-1111-1111-111111111111', status=status,
            aop_id='AOP:1', aop_label='Test AOP', logfc_threshold=1.0, pval_cutoff=0.05,
            selected_resources='WikiPathways', id_column='gene', fc_column='logFC', pval_column='pval',
            harmonised_background=json.dumps(['G1', 'G2']), harmonised_gene_count=2,
            batch_name='Route Test Batch', expires_at=_expiry(),
        )
        session.add(batch)
        session.flush()
        enrichment = [{'KE': 'KE:1', 'Title': 'Alpha', 'p_value': 0.0005, 'FDR': 0.001,
                       'num_overlap': 6, 'odds_ratio': 3.0, 'total_KE_genes_in_dataset': 20,
                       'Direction': '6↑'}]
        ke_gene = {'KE:1': [{'id': 'TP53', 'log2FC': 2.0, 'significant': True},
                            {'id': 'EGFR', 'log2FC': 0.1, 'significant': False}]}
        for pos in range(2):
            session.add(ConditionRecord(
                batch_id=batch.id, position=pos, filename=f'c{pos}.tsv',
                condition_label=f'C{pos}', dose=f'{pos+1}uM', timepoint='4hr', status='complete',
                gene_count=20, significant_genes=6,
                enrichment_json=json.dumps(enrichment),
                network_json=json.dumps(_NETWORK),
                ke_gene_json=json.dumps(ke_gene),
            ))
        session.commit()
        uuid = batch.uuid
    finally:
        session.close()
    return uuid


@pytest.fixture
def seeded_client(flask_client, temp_database, monkeypatch):
    """Flask client wired to a temp DB seeded with one complete batch."""
    monkeypatch.setattr(app_module, 'db_manager', temp_database)
    uuid = _seed_complete_batch(temp_database)
    return flask_client, uuid


class TestCompareExport:
    def test_csv(self, seeded_client):
        client, uuid = seeded_client
        r = client.get(f'/batch/{uuid}/compare/export?fmt=csv&matrix=fdr')
        assert r.status_code == 200
        assert 'text/csv' in r.content_type
        assert 'attachment' in r.headers['Content-Disposition']
        body = r.get_data(as_text=True)
        assert 'Key Event ID' in body
        assert 'KE:1' in body

    def test_xlsx(self, seeded_client):
        client, uuid = seeded_client
        r = client.get(f'/batch/{uuid}/compare/export?fmt=xlsx&matrix=neglog10')
        assert r.status_code == 200
        assert 'spreadsheetml' in r.content_type
        assert r.get_data()[:2] == b'PK'  # xlsx is a zip

    def test_bad_format(self, seeded_client):
        client, uuid = seeded_client
        assert client.get(f'/batch/{uuid}/compare/export?fmt=txt').status_code == 400


class TestGenesExport:
    def test_long_csv(self, seeded_client):
        client, uuid = seeded_client
        r = client.get(f'/batch/{uuid}/genes/export?fmt=csv&view=long')
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert 'Gene_Symbol' in body
        assert 'TP53' in body

    def test_summary_csv(self, seeded_client):
        client, uuid = seeded_client
        r = client.get(f'/batch/{uuid}/genes/export?fmt=csv&view=summary')
        assert r.status_code == 200
        body = r.get_data(as_text=True)
        assert 'shared' in body
        assert 'log2FC: C0' in body

    def test_xlsx_two_sheets(self, seeded_client):
        client, uuid = seeded_client
        r = client.get(f'/batch/{uuid}/genes/export?fmt=xlsx')
        assert r.status_code == 200
        assert r.get_data()[:2] == b'PK'

    def test_genes_data_json(self, seeded_client):
        client, uuid = seeded_client
        r = client.get(f'/batch/{uuid}/genes/data?ke=KE:1')
        assert r.status_code == 200
        payload = r.get_json()
        assert payload['ke'] == 'KE:1'
        assert payload['condition_labels'] == ['C0', 'C1']
        assert any(g['Gene_Symbol'] == 'TP53' for g in payload['genes'])

    def test_genes_data_missing_ke_param(self, seeded_client):
        client, uuid = seeded_client
        assert client.get(f'/batch/{uuid}/genes/data').status_code == 400


class TestBatchReport:
    def test_html(self, seeded_client):
        client, uuid = seeded_client
        r = client.get(f'/batch/{uuid}/report?format=html')
        assert r.status_code == 200
        assert 'text/html' in r.content_type
        body = r.get_data(as_text=True)
        assert 'Route Test Batch' in body
        assert 'C0' in body and 'C1' in body

    def test_pdf(self, seeded_client):
        client, uuid = seeded_client
        r = client.get(f'/batch/{uuid}/report?format=pdf')
        assert r.status_code == 200
        assert r.content_type == 'application/pdf'
        assert r.get_data()[:4] == b'%PDF'

    def test_bad_format(self, seeded_client):
        client, uuid = seeded_client
        assert client.get(f'/batch/{uuid}/report?format=docx').status_code == 400


class TestGuards:
    def test_non_complete_returns_409(self, flask_client, temp_database, monkeypatch):
        monkeypatch.setattr(app_module, 'db_manager', temp_database)
        # Seed a running batch (overwrite the default complete one's status).
        session = temp_database.get_session()
        try:
            batch = BatchRecord(uuid='22222222-2222-2222-2222-222222222222', status='running',
                                aop_id='AOP:1', batch_name='Running', expires_at=_expiry())
            session.add(batch)
            session.commit()
        finally:
            session.close()
        uuid = '22222222-2222-2222-2222-222222222222'
        for path in ('/compare/export', '/genes/export', '/report'):
            assert flask_client.get(f'/batch/{uuid}{path}').status_code == 409
        assert flask_client.get(f'/batch/{uuid}/genes/data?ke=KE:1').status_code == 409

    def test_malformed_uuid_returns_400(self, flask_client):
        for path in ('/compare/export', '/genes/export', '/report', '/genes/data?ke=KE:1'):
            assert flask_client.get(f'/batch/not-a-uuid{path}').status_code == 400
