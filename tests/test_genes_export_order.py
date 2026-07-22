"""Row-order determinism of the driver-gene tracking and export (issue #82).

``build_gene_tracking`` used to append its records in dict-iteration order over
each condition's stored ``ke_gene_json``, so the Genes tab payload, the batch
report and the CSV/Excel export all inherited an order that was an artefact of
how the JSON happened to be written. These tests pin the order at both ends:
directly on ``build_gene_tracking`` (the source) and through the export route
(the boundary), seeding batches that hold the same driver genes in different
insertion orders.
"""
import json
from datetime import datetime, timedelta, timezone

import pytest

import app as app_module
from database import BatchRecord, ConditionRecord
from services.comparison_service import build_gene_tracking


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


def _seed(db_manager, uuid, ke_gene_map, labels=('C0', 'C1')):
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
        for pos, label in enumerate(labels):
            # A list gives each condition its own map, so a gene can drive one
            # condition without driving the others.
            per_condition = (
                ke_gene_map[pos] if isinstance(ke_gene_map, (list, tuple))
                else ke_gene_map
            )
            session.add(ConditionRecord(
                batch_id=batch.id, position=pos, filename=f'c{pos}.tsv',
                condition_label=label, dose=label, timepoint='4hr',
                status='complete', gene_count=3, significant_genes=3,
                enrichment_json=json.dumps(_ENRICHMENT),
                ke_gene_json=json.dumps(per_condition),
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

# A dose series: label order and upload order disagree lexically.
_DOSE_MAP = {
    'KE:2': [_sig('TP53', 2.0), _sig('BRCA1', 1.5)],
    'KE:10': [_sig('EGFR', 1.1)],
}
_DOSE_LABELS = ('2uM', '10uM')

# Within KE:2, ZZZ drives both conditions and AAA only one. A plain
# (KE, gene) sort would list AAA first and silently drop the shared-first
# grouping the summary view exists to show.
_SHARED_MAPS = (
    {'KE:2': [_sig('ZZZ', 2.0), _sig('AAA', 1.5)]},
    {'KE:2': [_sig('ZZZ', 2.1)]},
)

_UUID_A = '22222222-2222-2222-2222-222222222222'
_UUID_B = '33333333-3333-3333-3333-333333333333'
_UUID_DOSE = '44444444-4444-4444-4444-444444444444'
_UUID_SHARED = '55555555-5555-5555-5555-555555555555'


@pytest.fixture
def two_batches(flask_client, temp_database, monkeypatch):
    monkeypatch.setattr(app_module, 'db_manager', temp_database)
    _seed(temp_database, _UUID_A, _ORDER_A)
    _seed(temp_database, _UUID_B, _ORDER_B)
    _seed(temp_database, _UUID_DOSE, _DOSE_MAP, labels=_DOSE_LABELS)
    _seed(temp_database, _UUID_SHARED, _SHARED_MAPS)
    return flask_client


def _rows(client, uuid, view='long'):
    body = client.get(
        f'/batch/{uuid}/genes/export?fmt=csv&view={view}'
    ).get_data(as_text=True)
    return [line.split(',') for line in body.strip().splitlines()[1:]]


class _StubCondition:
    """Minimal stand-in for a ConditionRecord for the unit-level tests."""

    def __init__(self, label, ke_gene_map):
        self.condition_label = label
        self.ke_gene_json = json.dumps(ke_gene_map)
        self.enrichment_json = json.dumps(_ENRICHMENT)


class TestGeneTrackingRecordOrder:
    """The source of the order, not just the export boundary."""

    def test_records_are_ordered_independently_of_json_key_order(self):
        first = build_gene_tracking([
            _StubCondition('C0', _ORDER_A), _StubCondition('C1', _ORDER_A),
        ])
        second = build_gene_tracking([
            _StubCondition('C0', _ORDER_B), _StubCondition('C1', _ORDER_B),
        ])
        assert first['records'] == second['records']

    def test_records_are_ke_numeric_then_gene_then_upload_position(self):
        tracking = build_gene_tracking([
            _StubCondition(label, _DOSE_MAP) for label in _DOSE_LABELS
        ])
        triples = [
            (r['KE_ID'], r['Gene_Symbol'], r['Condition'])
            for r in tracking['records']
        ]
        assert triples == [
            ('KE:2', 'BRCA1', '2uM'),
            ('KE:2', 'BRCA1', '10uM'),
            ('KE:2', 'TP53', '2uM'),
            ('KE:2', 'TP53', '10uM'),
            ('KE:10', 'EGFR', '2uM'),
            ('KE:10', 'EGFR', '10uM'),
        ]


class TestGeneExportDeterminism:
    @pytest.mark.parametrize('view', ['long', 'summary'])
    def test_identical_result_data_gives_identical_bytes(self, two_batches, view):
        client = two_batches
        first = client.get(f'/batch/{_UUID_A}/genes/export?fmt=csv&view={view}')
        second = client.get(f'/batch/{_UUID_B}/genes/export?fmt=csv&view={view}')
        assert first.status_code == second.status_code == 200
        assert first.get_data() == second.get_data()

    def test_summary_view_ke_order_is_numeric(self, two_batches):
        """The summary frame is sorted too, not merely stable by accident.

        The wide view was already insertion-order-independent, so byte equality
        between two batches proves nothing about it. Its KE order, however, was
        a lexical sort of the KE ID — 'KE:10' ahead of 'KE:2'.
        """
        rows = _rows(two_batches, _UUID_A, view='summary')
        ke_order = []
        for row in rows:
            if not ke_order or ke_order[-1] != row[0]:
                ke_order.append(row[0])
        assert ke_order == ['KE:2', 'KE:9', 'KE:10']

    def test_summary_view_keeps_shared_first_ordering(self, two_batches):
        """Sorting the export must not undo the summary's shared-first order.

        A gene driving several conditions is listed above one driving a single
        condition, which is the point of the summary view; re-sorting on Key
        Event and gene symbol alone would replace that with the alphabet.
        """
        rows = _rows(two_batches, _UUID_SHARED, view='summary')
        ke2 = [row[2] for row in rows if row[0] == 'KE:2']
        assert ke2 == ['ZZZ', 'AAA']

    def test_conditions_follow_upload_order_not_alphabetical(self, two_batches):
        """A dose series must not be reordered as text ('10uM' before '2uM')."""
        rows = _rows(two_batches, _UUID_DOSE)
        triples = [(r[0], r[2], r[3]) for r in rows]
        assert triples == [
            ('KE:2', 'BRCA1', '2uM'),
            ('KE:2', 'BRCA1', '10uM'),
            ('KE:2', 'TP53', '2uM'),
            ('KE:2', 'TP53', '10uM'),
            ('KE:10', 'EGFR', '2uM'),
            ('KE:10', 'EGFR', '10uM'),
        ]

    def test_ke_ids_sorted_numerically(self, two_batches):
        ke_order = []
        for row in _rows(two_batches, _UUID_A):
            if not ke_order or ke_order[-1] != row[0]:
                ke_order.append(row[0])
        assert ke_order == ['KE:2', 'KE:9', 'KE:10']

    def test_genes_sorted_within_ke(self, two_batches):
        rows = _rows(two_batches, _UUID_A)
        ke2 = [r[2] for r in rows if r[0] == 'KE:2']
        assert ke2 == sorted(ke2)
