"""Row-order determinism of the driver-gene export (issue #82).

The long view used to inherit its order from dict iteration over each
condition's stored KE->gene map, so two runs over identical result data
produced CSVs that were equal only after sorting. These tests seed two
batches holding the same driver genes in different insertion orders and
assert the exports are byte-identical.
"""
import json
from datetime import datetime, timedelta, timezone

import pytest

import app as app_module
from database import BatchRecord, ConditionRecord


def _expiry():
    return datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(days=14)


_ENRICHMENT = [
    {'KE': 'KE:2', 'Title': 'Beta', 'p_value': 0.001, 'FDR': 0.01,
     'num_overlap': 2, 'odds_ratio': 2.0, 'total_KE_genes_in_dataset': 10},
    {'KE': 'KE:10', 'Title': 'Gamma', 'p_value': 0.002, 'FDR': 0.02,
     'num_overlap': 2, 'odds_ratio': 2.0, 'total_KE_genes_in_dataset': 10},
    {'KE': 'KE:9', 'Title': 'Delta', 'p_value': 0.003, 'FDR': 0.03,
     'num_overlap': 2, 'odds_ratio': 2.0, 'total_KE_genes_in_dataset': 10},
]


def _seed(db_manager, uuid, ke_gene_map):
    session = db_manager.get_session()
    try:
        batch = BatchRecord(
            uuid=uuid, status='complete', aop_id='AOP:1', aop_label='Test AOP',
            logfc_threshold=1.0, pval_cutoff=0.05, selected_resources='WikiPathways',
            id_column='gene', fc_column='logFC', pval_column='pval',
            harmonised_background=json.dumps(['TP53', 'EGFR', 'BRCA1']),
            harmonised_gene_count=3, batch_name='Order Test', expires_at=_expiry(),
        )
        session.add(batch)
        session.flush()
        for pos in range(2):
            session.add(ConditionRecord(
                batch_id=batch.id, position=pos, filename=f'c{pos}.tsv',
                condition_label=f'C{pos}', dose=f'{pos + 1}uM', timepoint='4hr',
                status='complete', gene_count=3, significant_genes=3,
                enrichment_json=json.dumps(_ENRICHMENT),
                ke_gene_json=json.dumps(ke_gene_map),
            ))
        session.commit()
    finally:
        session.close()


def _sig(gene, logfc):
    return {'id': gene, 'log2FC': logfc, 'significant': True}


# Same content, different insertion order for both KEs and genes.
_ORDER_A = {
    'KE:2': [_sig('TP53', 2.0), _sig('BRCA1', 1.5)],
    'KE:10': [_sig('EGFR', 1.1), _sig('TP53', 2.2)],
    'KE:9': [_sig('BRCA1', -1.7)],
}
_ORDER_B = {
    'KE:9': [_sig('BRCA1', -1.7)],
    'KE:10': [_sig('TP53', 2.2), _sig('EGFR', 1.1)],
    'KE:2': [_sig('BRCA1', 1.5), _sig('TP53', 2.0)],
}

_UUID_A = '22222222-2222-2222-2222-222222222222'
_UUID_B = '33333333-3333-3333-3333-333333333333'


@pytest.fixture
def two_batches(flask_client, temp_database, monkeypatch):
    monkeypatch.setattr(app_module, 'db_manager', temp_database)
    _seed(temp_database, _UUID_A, _ORDER_A)
    _seed(temp_database, _UUID_B, _ORDER_B)
    return flask_client


class TestGeneExportDeterminism:
    @pytest.mark.parametrize('view', ['long', 'summary'])
    def test_identical_result_data_gives_identical_bytes(self, two_batches, view):
        client = two_batches
        first = client.get(f'/batch/{_UUID_A}/genes/export?fmt=csv&view={view}')
        second = client.get(f'/batch/{_UUID_B}/genes/export?fmt=csv&view={view}')
        assert first.status_code == second.status_code == 200
        assert first.get_data() == second.get_data()

    def test_repeated_export_is_stable(self, two_batches):
        client = two_batches
        bodies = {
            client.get(f'/batch/{_UUID_A}/genes/export?fmt=csv&view=long').get_data()
            for _ in range(3)
        }
        assert len(bodies) == 1

    def test_ke_ids_sorted_numerically(self, two_batches):
        client = two_batches
        body = client.get(
            f'/batch/{_UUID_A}/genes/export?fmt=csv&view=long'
        ).get_data(as_text=True)
        rows = [line.split(',') for line in body.strip().splitlines()[1:]]
        ke_order = []
        for row in rows:
            if not ke_order or ke_order[-1] != row[0]:
                ke_order.append(row[0])
        assert ke_order == ['KE:2', 'KE:9', 'KE:10']

    def test_genes_sorted_within_ke(self, two_batches):
        client = two_batches
        body = client.get(
            f'/batch/{_UUID_A}/genes/export?fmt=csv&view=long'
        ).get_data(as_text=True)
        rows = [line.split(',') for line in body.strip().splitlines()[1:]]
        ke2 = [r[2] for r in rows if r[0] == 'KE:2']
        assert ke2 == sorted(ke2)
