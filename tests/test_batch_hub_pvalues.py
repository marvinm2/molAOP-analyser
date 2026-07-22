"""End-to-end coverage for per-gene p-values in a batch run (issue #75).

The hub-gene panel on a batch condition page reads its "Adj. p-value" out of
the KE-gene payload stored on the ConditionRecord. The batch runner never built
the per-gene p-value maps, so every hub row on every condition page rendered an
em dash while the same file analysed on its own showed a p-value.

Seeding a ConditionRecord with p-values already in its ``ke_gene_json`` only
exercises the renderer, which was never broken. These tests run the real batch
pipeline over a real upload and then read the page, so the assertion is on the
producer.
"""
import json
import os
import re
import uuid as uuid_lib
from datetime import datetime, timedelta, timezone

import pandas as pd
import pytest

import app as app_module
from database import BatchRecord, ConditionRecord
from services.batch_service import run_batch


# Ten genes per KE so Fisher's test is actually computed (MIN_KE_GENES = 5).
REFERENCE_SETS = {
    'KE:115': {f'G{i}' for i in range(1, 11)},
    'KE:116': {f'G{i}' for i in range(11, 21)},
}


def _write_dataset(path, *, adj_column=True):
    """Write a limma-shaped TSV: G1..G10 significant, G11..G20 not.

    ``adj.P.Val`` carries a value distinct from ``P.Value`` so a test can tell
    which column reached the panel.
    """
    rows = []
    for i in range(1, 21):
        sig = i <= 10
        row = {
            'gene': f'G{i}',
            'logFC': 2.0 if sig else 0.1,
            'P.Value': 0.0001 if sig else 0.5,
        }
        if adj_column:
            row['adj.P.Val'] = 0.0004 if sig else 0.8
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)


@pytest.fixture
def aop_stub(monkeypatch):
    """Stub load_aop_data so the batch needs no real AOP topology."""
    def _fake_load_aop_data(aop_id):
        ke_list = {'KE:115', 'KE:116'}
        edges = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])
        return ke_list, edges, {'KE:115': 'MIE', 'KE:116': 'AO'}, \
            {'KE:115': 'Event 115', 'KE:116': 'Event 116'}
    monkeypatch.setattr('services.data_service.load_aop_data', _fake_load_aop_data)


@pytest.fixture
def upload_root(tmp_path, monkeypatch):
    root = tmp_path / 'uploads'
    root.mkdir()
    monkeypatch.setattr('config.Config.UPLOAD_FOLDER', str(root))
    return str(root)


@pytest.fixture
def ran_batch(temp_database, aop_stub, upload_root, monkeypatch):
    """Run a real one-condition batch and hand back the client and its UUID."""
    monkeypatch.setattr(app_module, 'db_manager', temp_database)

    def _run(*, adj_column=True, pval_column='adj.P.Val'):
        batch_uuid = str(uuid_lib.uuid4())
        batch_dir = os.path.join(upload_root, batch_uuid)
        os.makedirs(batch_dir, exist_ok=True)
        _write_dataset(os.path.join(batch_dir, 'cond0.tsv'), adj_column=adj_column)

        harmonised = sorted({f'G{i}' for i in range(1, 21)})
        session = temp_database.get_session()
        try:
            batch = BatchRecord(
                uuid=batch_uuid, status='pending', aop_id='AOP:1', aop_label='Test AOP',
                logfc_threshold=1.0, pval_cutoff=0.05,
                selected_resources='WikiPathways', min_confidence='all',
                id_column='gene', fc_column='logFC', pval_column=pval_column,
                harmonised_background=json.dumps(harmonised),
                harmonised_gene_count=len(harmonised),
                batch_name='Hub p-value batch',
                expires_at=datetime.now(timezone.utc).replace(tzinfo=None)
                + timedelta(days=14),
            )
            session.add(batch)
            session.flush()
            session.add(ConditionRecord(
                batch_id=batch.id, position=0, filename='cond0.tsv',
                condition_label='Cond 0', status='pending',
            ))
            session.commit()
            batch_id = batch.id
        finally:
            session.close()

        run_batch(batch_id, temp_database.db_url, REFERENCE_SETS)
        return batch_id, batch_uuid

    return _run


def _stored_ke_gene_map(db_manager, batch_id):
    session = db_manager.get_session()
    try:
        cond = session.query(ConditionRecord).filter_by(
            batch_id=batch_id, position=0
        ).first()
        assert cond.status == 'complete', cond.error_message
        return json.loads(cond.ke_gene_json)
    finally:
        session.close()


class TestStoredPvalues:
    """The batch runner must persist per-gene p-values, not just fold changes."""

    def test_condition_record_carries_adjusted_pvalues(self, ran_batch, temp_database):
        batch_id, _ = ran_batch()
        ke_gene_map = _stored_ke_gene_map(temp_database, batch_id)
        entries = {g['id']: g for g in ke_gene_map['KE:115']}
        assert entries['G1']['pvalue_adj'] == pytest.approx(0.0004)
        assert entries['G1']['pvalue_raw'] == pytest.approx(0.0001)

    def test_raw_only_upload_stores_the_raw_pvalue(self, ran_batch, temp_database):
        """A batch with no adjusted column must not pass a raw value off as FDR."""
        batch_id, _ = ran_batch(adj_column=False, pval_column='P.Value')
        ke_gene_map = _stored_ke_gene_map(temp_database, batch_id)
        entries = {g['id']: g for g in ke_gene_map['KE:115']}
        assert entries['G1']['pvalue_raw'] == pytest.approx(0.0001)
        assert entries['G1']['pvalue_adj'] is None


class TestHubPanelRendering:
    """The stored p-values must reach the hub table on the condition page."""

    def test_hub_row_shows_the_adjusted_pvalue(self, ran_batch, flask_client):
        _, batch_uuid = ran_batch()
        html = flask_client.get(
            f'/batch/{batch_uuid}/condition/0'
        ).get_data(as_text=True)
        row = re.search(r'<td>G1</td>.*?</tr>', html, re.DOTALL)
        assert row, 'G1 hub row not rendered'
        assert '4.00e-04' in row.group(0)
        assert '&mdash;' not in row.group(0)

    def test_hub_row_dashes_only_when_there_is_no_adjusted_column(
        self, ran_batch, flask_client
    ):
        _, batch_uuid = ran_batch(adj_column=False, pval_column='P.Value')
        html = flask_client.get(
            f'/batch/{batch_uuid}/condition/0'
        ).get_data(as_text=True)
        row = re.search(r'<td>G1</td>.*?</tr>', html, re.DOTALL)
        assert row and '&mdash;' in row.group(0)
