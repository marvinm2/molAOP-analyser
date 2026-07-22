"""Route tests for /batch/<uuid>/condition/<n> (issues #74, #75).

The batch-condition page re-renders results.html from stored ConditionRecord
blobs. Everything the header states about the run therefore has to be carried
over from the BatchRecord explicitly — anything the route forgets is either
silently absent or, worse, replaced by a template default that names a real
resource. Issue #74 was the latter: a run over all three gene-set resources
reported "WikiPathways", which is a provenance claim a reader would carry into
a methods section.

The hub-gene assertions pin the other half of the contract (#75): whatever
per-gene adjusted p-values the stored ke_gene_json carries must reach the hub
table unaltered.
"""
import json
import re
from datetime import datetime, timedelta, timezone

import pytest

import app as app_module
from database import BatchRecord, ConditionRecord


def _expiry():
    return datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(days=14)


_NETWORK = {
    'nodes': [
        {'data': {'id': 'KE:1', 'label': 'Alpha', 'ke_type': 'MIE'}, 'classes': 'significant'},
        {'data': {'id': 'KE:2', 'label': 'Beta', 'ke_type': 'KE'}, 'classes': ''},
        {'data': {'id': 'KE:3', 'label': 'Gamma', 'ke_type': 'AO'}, 'classes': ''},
    ],
    'edges': [{'data': {'source': 'KE:1', 'target': 'KE:2', 'id': 'KER:1'}}],
}

_RESOLUTION = [
    {'resource': 'WikiPathways', 'status': 'loaded', 'source': 'api',
     'ke_count': 42, 'confidence_applied': True, 'error': None},
    {'resource': 'GO_BP', 'status': 'loaded', 'source': 'api',
     'ke_count': 7, 'confidence_applied': False, 'error': None},
    {'resource': 'Reactome', 'status': 'loaded', 'source': 'api',
     'ke_count': 5, 'confidence_applied': False, 'error': None},
]


def _gene(symbol, log2fc, significant, pvalue_adj=None):
    entry = {'id': symbol, 'log2FC': log2fc, 'significant': significant}
    if pvalue_adj is not None:
        entry['pvalue_raw'] = pvalue_adj / 2
        entry['pvalue_adj'] = pvalue_adj
    return entry


def _seed(db_manager, *, selected_resources='WikiPathways, GO_BP, Reactome',
          resolution=_RESOLUTION, with_pvalues=True):
    """Create one complete batch with a single complete condition."""
    session = db_manager.get_session()
    try:
        batch = BatchRecord(
            uuid='22222222-2222-2222-2222-222222222222', status='complete',
            aop_id='AOP:1', aop_label='Test AOP', logfc_threshold=1.0, pval_cutoff=0.05,
            selected_resources=selected_resources,
            min_confidence='all',
            resource_resolution=json.dumps(resolution) if resolution else None,
            id_column='gene', fc_column='logFC', pval_column='adj.P.Val',
            harmonised_background=json.dumps(['TP53', 'EGFR']), harmonised_gene_count=2,
            batch_name='Condition Route Batch', expires_at=_expiry(),
        )
        session.add(batch)
        session.flush()

        adj = 0.0004 if with_pvalues else None
        # TP53 sits in all three KEs, so it clears the hub threshold (3 KEs).
        ke_gene = {
            'KE:1': [_gene('TP53', 2.0, True, adj), _gene('EGFR', 0.1, False, 0.42 if with_pvalues else None)],
            'KE:2': [_gene('TP53', 2.0, True, adj)],
            'KE:3': [_gene('TP53', 2.0, True, adj)],
        }
        session.add(ConditionRecord(
            batch_id=batch.id, position=0, filename='c0.tsv',
            condition_label='C0', dose='1uM', timepoint='4hr', status='complete',
            gene_count=20, significant_genes=6,
            completed_at=datetime.now(timezone.utc).replace(tzinfo=None),
            enrichment_json=json.dumps([{
                'KE': 'KE:1', 'Title': 'Alpha', 'p_value': 0.0005, 'FDR': 0.001,
                'num_overlap': 6, 'odds_ratio': 3.0, 'total_KE_genes_in_dataset': 20,
                'Direction': '6↑', 'Representation': 'enriched',
            }]),
            network_json=json.dumps(_NETWORK),
            ke_gene_json=json.dumps(ke_gene),
            ke_type_map_json=json.dumps({'KE:1': 'MIE', 'KE:2': 'KE', 'KE:3': 'AO'}),
            ke_title_map_json=json.dumps({'KE:1': 'Alpha', 'KE:2': 'Beta', 'KE:3': 'Gamma'}),
        ))
        session.commit()
        return batch.uuid
    finally:
        session.close()


@pytest.fixture
def condition_client(flask_client, temp_database, monkeypatch):
    monkeypatch.setattr(app_module, 'db_manager', temp_database)

    def _make(**kwargs):
        uuid = _seed(temp_database, **kwargs)
        return flask_client, uuid

    return _make


def _resource_header_value(html):
    """Extract the rendered 'Gene Set Resources' metadata value."""
    match = re.search(
        r'Gene Set Resources</span>.*?metadata-summary__value">(.*?)</span>',
        html, re.DOTALL,
    )
    assert match, 'Gene Set Resources header not found on the page'
    return match.group(1).strip()


class TestResourceProvenance:
    """Issue #74: the header must state the batch's real resource selection."""

    def test_all_selected_resources_are_reported(self, condition_client):
        client, uuid = condition_client()
        r = client.get(f'/batch/{uuid}/condition/0')
        assert r.status_code == 200
        value = _resource_header_value(r.get_data(as_text=True))
        assert value == 'WikiPathways, GO_BP, Reactome'

    def test_resolution_provenance_line_is_rendered(self, condition_client):
        client, uuid = condition_client()
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        assert 'Gene set provenance' in html
        for resource in ('WikiPathways', 'GO_BP', 'Reactome'):
            assert resource in html

    def test_skipped_resource_raises_a_warning(self, condition_client):
        resolution = [
            dict(_RESOLUTION[0]),
            {'resource': 'Reactome', 'status': 'error', 'source': None,
             'ke_count': 0, 'confidence_applied': False, 'error': 'timeout'},
        ]
        client, uuid = condition_client(
            selected_resources='WikiPathways, Reactome', resolution=resolution
        )
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        assert 'unavailable, skipped' in html
        assert 'Fewer Key Events were testable' in html

    def test_missing_selection_is_not_reported_as_wikipathways(self, condition_client):
        """A batch predating #55 recorded no selection — say so, do not invent one."""
        client, uuid = condition_client(selected_resources=None, resolution=None)
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        value = _resource_header_value(html)
        assert value == 'Not recorded'
        assert 'WikiPathways' not in value


class TestHubGenePvalues:
    """Issue #75: stored per-gene adjusted p-values must reach the hub table."""

    def test_adjusted_pvalue_is_rendered_for_hub_genes(self, condition_client):
        client, uuid = condition_client(with_pvalues=True)
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        row = re.search(r'<td>TP53</td>.*?</tr>', html, re.DOTALL)
        assert row, 'TP53 hub row not rendered'
        assert '4.00e-04' in row.group(0)
        assert '&mdash;' not in row.group(0)

    def test_absent_pvalue_renders_as_a_dash(self, condition_client):
        """The em dash is only correct when the stored map genuinely has none."""
        client, uuid = condition_client(with_pvalues=False)
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        row = re.search(r'<td>TP53</td>.*?</tr>', html, re.DOTALL)
        assert row and '&mdash;' in row.group(0)
