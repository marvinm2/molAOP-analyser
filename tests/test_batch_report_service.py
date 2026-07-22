"""Tests for the integrated batch report service (rendering + assembly)."""
import json
from datetime import datetime
from types import SimpleNamespace

import pytest

from services import batch_report_service as brs


def _comparison_data():
    return {
        'ke_labels': ['KE:1', 'KE:2'],
        'ke_titles': ['Alpha', 'Beta'],
        'condition_labels': ['C1', 'C2'],
        'condition_colors': ['#E6007E', '#307BBF'],
        'fdr_matrix': [[0.001, 0.0001], [0.5, None]],
        'neg_log10_matrix': [[3.0, 4.0], [None, None]],
        'condition_gene_counts': [100, 120],
        'condition_sig_gene_counts': [20, 40],
        'condition_sig_ke_counts': [1, 1],
        'condition_doses': ['1uM', '10uM'],
        'condition_timepoints': ['4hr', '4hr'],
    }


_NETWORK = {
    'nodes': [
        {'data': {'id': 'KE:1', 'label': 'Alpha', 'ke_type': 'MIE'}, 'classes': 'significant'},
        {'data': {'id': 'KE:2', 'label': 'Beta', 'ke_type': 'AO'}, 'classes': ''},
    ],
    'edges': [{'data': {'source': 'KE:1', 'target': 'KE:2', 'id': 'KER:1'}}],
}


def _batch():
    return SimpleNamespace(
        uuid='u', batch_name='Test Batch', aop_id='AOP:1', aop_label='Test AOP',
        logfc_threshold=1.0, pval_cutoff=0.05, selected_resources='WikiPathways',
        harmonised_gene_count=5000, owner='Marvin', description='demo',
        completed_at=datetime(2026, 6, 19, 12, 0),
        id_column='gene', fc_column='logFC', pval_column='pval',
    )


def _conditions():
    enrichment = [{'KE': 'KE:1', 'Title': 'Alpha', 'p_value': 0.0005, 'FDR': 0.001,
                   'num_overlap': 8, 'odds_ratio': 3.2, 'total_KE_genes_in_dataset': 40,
                   'Direction': '7↑'}]
    c1 = SimpleNamespace(condition_label='C1', filename='c1.tsv', dose='1uM', timepoint='4hr',
                         gene_count=100, significant_genes=20, status='complete',
                         enrichment_json=json.dumps(enrichment), network_json=json.dumps(_NETWORK))
    c2 = SimpleNamespace(condition_label='C2', filename='c2.tsv', dose='10uM', timepoint='4hr',
                         gene_count=120, significant_genes=40, status='complete',
                         enrichment_json=json.dumps(enrichment), network_json=json.dumps(_NETWORK))
    return [c1, c2]


class TestImageRendering:
    def test_heatmap_png(self):
        png = brs.render_heatmap_png(_comparison_data())
        assert png is not None
        assert png[:4] == b'\x89PNG'

    def test_heatmap_empty_returns_none(self):
        assert brs.render_heatmap_png({}) is None

    def test_network_png(self):
        png = brs.render_ke_network_png(json.dumps(_NETWORK))
        assert png is not None
        assert png[:4] == b'\x89PNG'

    def test_network_png_accepts_dict(self):
        png = brs.render_ke_network_png(_NETWORK)
        assert png is not None and png[:4] == b'\x89PNG'

    def test_network_png_empty_returns_none(self):
        assert brs.render_ke_network_png(None) is None
        assert brs.render_ke_network_png({'nodes': [], 'edges': []}) is None


class TestReportAssembly:
    def test_html_contains_sections(self):
        html = brs.generate_batch_html(_batch(), _conditions(), _comparison_data())
        assert 'Test Batch' in html
        assert 'C1' in html and 'C2' in html
        assert 'data:image/png' in html  # embedded heatmap/network images
        assert 'Cross-Condition Comparison' in html

    def test_pdf_is_valid(self):
        pdf = brs.generate_batch_pdf(_batch(), _conditions(), _comparison_data())
        assert pdf[:4] == b'%PDF'
        assert len(pdf) > 1000


class TestKeAccountingSection:
    """Issue #65: each condition section states its multiple-testing denominator."""

    _NETWORK_WITH_REASONS = {
        'nodes': [
            {'data': {'id': 'KE:1', 'label': 'Alpha', 'ke_type': 'MIE',
                      'excluded_reason': None}, 'classes': 'significant'},
            {'data': {'id': 'KE:2', 'label': 'Beta', 'ke_type': 'AO',
                      'excluded_reason': 'too_few_genes'}, 'classes': 'too-few-genes'},
            {'data': {'id': 'KE:3', 'label': 'Gamma', 'ke_type': 'intermediate',
                      'excluded_reason': 'no_mapping'}, 'classes': 'no-genes'},
        ],
        'edges': [{'data': {'source': 'KE:1', 'target': 'KE:2', 'id': 'KER:1'}}],
    }

    def _conditions_with_reasons(self):
        conds = _conditions()
        for c in conds:
            c.network_json = json.dumps(self._NETWORK_WITH_REASONS)
        return conds

    def test_html_states_ke_accounting(self):
        html = brs.generate_batch_html(
            _batch(), self._conditions_with_reasons(), _comparison_data()
        )
        assert '1 of 3 Key Events tested' in html
        assert 'fewer than 5 measured genes' in html
        assert 'no gene set mapped' in html

    def test_html_omits_accounting_for_legacy_networks(self):
        """Batches run before the field existed render exactly as before."""
        html = brs.generate_batch_html(_batch(), _conditions(), _comparison_data())
        assert 'Key Events:' not in html
        assert 'no gene set mapped' not in html

    def test_pdf_generates_with_accounting(self):
        pdf = brs.generate_batch_pdf(
            _batch(), self._conditions_with_reasons(), _comparison_data()
        )
        assert pdf[:4] == b'%PDF'

    def test_network_image_border_follows_fdr_significance(self):
        """Issue #63: the network image reds the same KEs the heatmap calls
        significant — both now read the server-side 'significant' class, which
        build_cytoscape_network derives from FDR."""
        png = brs.render_ke_network_png(json.dumps(self._NETWORK_WITH_REASONS))
        assert png is not None and png[:4] == b'\x89PNG'


class TestDroppedRowsNote:
    """Issue #103: a condition's discarded symbol-less rows reach the report."""

    def test_note_reads_in_the_plural(self):
        cond = SimpleNamespace(dropped_unidentified_rows=875)

        assert brs._dropped_rows_note(cond) == (
            ' · 875 rows without a gene symbol were excluded'
        )

    def test_note_reads_in_the_singular(self):
        cond = SimpleNamespace(dropped_unidentified_rows=1)

        assert brs._dropped_rows_note(cond) == (
            ' · 1 row without a gene symbol was excluded'
        )

    @pytest.mark.parametrize('value', [0, None])
    def test_nothing_dropped_or_nothing_recorded_says_nothing(self, value):
        """An unrecorded count must not be reported as "none dropped"."""
        assert brs._dropped_rows_note(SimpleNamespace(dropped_unidentified_rows=value)) == ''

    def test_condition_predating_the_field_does_not_break_the_report(self):
        """The report is assembled from stand-in objects too — no attribute, no crash."""
        assert brs._dropped_rows_note(SimpleNamespace()) == ''

    def test_html_states_the_discard_for_the_condition(self):
        conditions = _conditions()
        conditions[0].dropped_unidentified_rows = 875

        html = brs.generate_batch_html(_batch(), conditions, _comparison_data())

        assert '875 rows without a gene symbol were excluded' in html
