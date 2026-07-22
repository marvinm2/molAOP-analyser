"""Route tests for /batch/<uuid>/condition/<n> and its report export (#74).

The batch-condition page re-renders results.html from stored ConditionRecord
blobs. Everything the header states about the run therefore has to be carried
over from the BatchRecord explicitly — anything the route forgets is either
silently absent or, worse, replaced by a template default that names a real
resource. Issue #74 was the latter: a run over all three gene-set resources
reported "WikiPathways", which is a provenance claim a reader would carry into
a methods section.

The same applies to the report exported from that page: a hidden field the form
does not post is a value the generator silently borrows from the session.

Per-gene p-values on the condition page (#75) are covered end-to-end in
tests/test_batch_hub_pvalues.py — asserting them against a hand-seeded
ke_gene_json would only exercise the renderer, which was never at fault.
"""
import json
import re
from datetime import datetime, timedelta, timezone
from unittest.mock import patch

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
          resolution=_RESOLUTION, min_confidence='all'):
    """Create one complete batch with a single complete condition."""
    session = db_manager.get_session()
    try:
        batch = BatchRecord(
            uuid='22222222-2222-2222-2222-222222222222', status='complete',
            aop_id='AOP:1', aop_label='Test AOP', logfc_threshold=1.0, pval_cutoff=0.05,
            selected_resources=selected_resources,
            min_confidence=min_confidence,
            resource_resolution=json.dumps(resolution) if resolution else None,
            id_column='gene', fc_column='logFC', pval_column='adj.P.Val',
            harmonised_background=json.dumps(['TP53', 'EGFR']), harmonised_gene_count=2,
            batch_name='Condition Route Batch', expires_at=_expiry(),
        )
        session.add(batch)
        session.flush()

        # TP53 sits in all three KEs, so it clears the hub threshold (3 KEs).
        ke_gene = {
            'KE:1': [_gene('TP53', 2.0, True, 0.0004), _gene('EGFR', 0.1, False, 0.42)],
            'KE:2': [_gene('TP53', 2.0, True, 0.0004)],
            'KE:3': [_gene('TP53', 2.0, True, 0.0004)],
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


def _metadata_value(html, label):
    """Extract a rendered value from the metadata-summary header by its label."""
    match = re.search(
        re.escape(label) + r'</span>.*?metadata-summary__value">(.*?)</span>',
        html, re.DOTALL,
    )
    assert match, f'{label!r} header not found on the page'
    return match.group(1).strip()


def _resource_header_value(html):
    """Extract the rendered 'Gene Set Resources' metadata value."""
    return _metadata_value(html, 'Gene Set Resources')


def _hidden_field(html, name):
    """Extract the value a report-form hidden input posts."""
    match = re.search(
        r'<input type="hidden" name="' + re.escape(name) + r'" value="(.*?)">', html
    )
    assert match, f'hidden field {name!r} not present in the report form'
    return match.group(1)


class TestResourceProvenance:
    """Issue #74: the header must state the batch's real resource selection."""

    def test_all_selected_resources_are_reported(self, condition_client):
        client, uuid = condition_client()
        r = client.get(f'/batch/{uuid}/condition/0')
        assert r.status_code == 200
        value = _resource_header_value(r.get_data(as_text=True))
        assert value == 'WikiPathways, GO_BP, Reactome'

    def test_resolution_provenance_line_is_rendered(self, condition_client):
        """The provenance line itself must name every resource and its KE count.

        Asserting each resource name appears *somewhere* on the page proves
        nothing — the header fixed above already lists all three. Scope the
        assertions to the provenance paragraph.
        """
        client, uuid = condition_client()
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        line = re.search(
            r'<strong>Gene set provenance:</strong>(.*?)</p>', html, re.DOTALL
        )
        assert line, 'Gene set provenance line not rendered'
        provenance = line.group(1)
        for resource, ke_count in (('WikiPathways', 42), ('GO_BP', 7), ('Reactome', 5)):
            assert resource in provenance
            assert f'({ke_count} KEs)' in provenance

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


class TestRunSettingsProvenance:
    """Issues #60/#74: the header must state the batch's own run settings."""

    def test_confidence_threshold_is_the_batch_threshold(self, condition_client):
        client, uuid = condition_client(min_confidence='high')
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        assert _metadata_value(html, 'Min. Mapping Confidence') == 'High only'

    def test_confidence_threshold_reaches_the_report_form(self, condition_client):
        client, uuid = condition_client(min_confidence='medium')
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        assert _hidden_field(html, 'min_confidence') == 'medium'

    def test_selected_resources_are_posted_to_the_report(self, condition_client):
        """Without this field the report falls back to the session (#74)."""
        client, uuid = condition_client()
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        assert _hidden_field(html, 'selected_resources') == 'WikiPathways, GO_BP, Reactome'

    def test_unrecorded_selection_is_posted_as_empty_not_wikipathways(
        self, condition_client
    ):
        client, uuid = condition_client(selected_resources=None, resolution=None)
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        assert _hidden_field(html, 'selected_resources') == ''

    def test_column_mapping_comes_from_the_batch(self, condition_client):
        client, uuid = condition_client()
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        assert _hidden_field(html, 'id_column') == 'gene'
        assert _hidden_field(html, 'fc_column') == 'logFC'
        assert _hidden_field(html, 'pval_column') == 'adj.P.Val'


class TestConditionMethod:
    """The page must render and export the method the batch actually ran."""

    def test_method_defaults_to_ora_without_the_column(self, condition_client):
        client, uuid = condition_client()
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        assert _hidden_field(html, 'method') == 'ora'
        assert 'Gene Expression Scale' in html

    def test_gsea_batch_renders_and_exports_gsea(self, condition_client, monkeypatch):
        """A BatchRecord carrying method='gsea' must not be reported as ORA.

        BatchRecord has no ``method`` column on this branch, so the attribute is
        planted on the model the route queries. That is the shape
        ``getattr(batch, 'method', None) or 'ora'`` has to cope with once the
        column lands on the other branch.
        """
        client, uuid = condition_client()

        # Make the batch look like a GSEA run without adding a column this
        # branch does not own.
        monkeypatch.setattr(BatchRecord, 'method', 'gsea', raising=False)
        html = client.get(f'/batch/{uuid}/condition/0').get_data(as_text=True)
        assert _hidden_field(html, 'method') == 'gsea'
        assert 'KE NES Scale' in html
        assert 'Gene Expression Scale' not in html


class TestReportResourceProvenance:
    """Issue #74: an exported report must not inherit another run's resources.

    ``generate_report`` fell back to ``session['experiment_metadata']`` when the
    form carried no resource list, and to the literal ``'WikiPathways'`` when
    the session was empty too. A PDF exported from a batch condition page
    therefore reported the resources of whatever single analysis the same
    browser session had last run.
    """

    _BASE_FORM = {
        'format': 'html',
        'filename': 'test.csv',
        'gene_count': '100',
        'significant_genes': '10',
        'aop_id': 'AOP:1',
        'aop_label': 'Test AOP',
        'logfc_threshold': '1.0',
        'pval_cutoff': '0.05',
        'id_column': 'gene',
        'fc_column': 'logFC',
        'pval_column': 'adj.P.Val',
        'id_type': 'HGNC',
        'enrichment_results': '[]',
    }

    def _capture(self, client, form):
        """POST the report form and return the ReportData the generator saw."""
        captured = {}

        def _fake_html_report(report_data):
            captured['data'] = report_data
            return '<html><body>Report</body></html>'

        with patch('app.report_generator.generate_html_report', _fake_html_report):
            response = client.post('/generate_report', data=form)
        assert response.status_code == 200, response.data[:200]
        return captured['data']

    def test_posted_resources_win_over_the_session(self, flask_client):
        with flask_client.session_transaction() as sess:
            sess['experiment_metadata'] = {
                'dataset_id': 'PREVIOUS RUN',
                'selected_resources': 'GO_BP',
            }
        form = dict(self._BASE_FORM, selected_resources='WikiPathways, Reactome')
        report = self._capture(flask_client, form)
        assert report.selected_resources == 'WikiPathways, Reactome'

    def test_empty_posted_resources_do_not_fall_back_to_the_session(self, flask_client):
        """A batch predating #55 posts an empty field. That is not consent to
        report the last single analysis's resources."""
        with flask_client.session_transaction() as sess:
            sess['experiment_metadata'] = {
                'dataset_id': 'PREVIOUS RUN',
                'selected_resources': 'GO_BP, Reactome',
            }
        form = dict(self._BASE_FORM, selected_resources='')
        report = self._capture(flask_client, form)
        assert report.selected_resources == 'Not recorded'

    def test_absent_field_and_empty_session_is_not_wikipathways(self, flask_client):
        report = self._capture(flask_client, dict(self._BASE_FORM))
        assert report.selected_resources == 'Not recorded'
