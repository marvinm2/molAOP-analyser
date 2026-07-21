"""Tests for network and comparison services."""
import json
import math
import numpy as np
import pytest
import pandas as pd
from types import SimpleNamespace

from config import Config
from services.network_service import build_cytoscape_network, ke_accounting_from_network
from services.comparison_service import build_comparison_matrix, CONDITION_PALETTE


class TestBuildCytoscapeNetwork:
    """Tests for Cytoscape network generation."""

    def _make_enrichment(self, ke_ids, fdrs):
        return pd.DataFrame({
            'KE': ke_ids,
            'FDR': fdrs,
            'odds_ratio': [2.0] * len(ke_ids),
            'p_value': fdrs,
            'num_overlap': [5] * len(ke_ids),
            'total_KE_genes_in_dataset': [20] * len(ke_ids),
            'pct_sig_in_KE': [25.0] * len(ke_ids),
        })

    def test_basic_network(self):
        ke_list = {'KE:1', 'KE:2'}
        edges = pd.DataFrame({
            'Source_KE': ['KE:1'], 'Target_KE': ['KE:2'],
            'KER_ID': ['KER:1'], 'AOP_ID': ['AOP:1'],
        })
        enrichment = self._make_enrichment(['KE:1', 'KE:2'], [0.01, 0.1])
        ke_title_map = {'KE:1': 'Event 1', 'KE:2': 'Event 2'}
        ke_type_map = {'KE:1': 'MIE', 'KE:2': 'AO'}

        result = build_cytoscape_network(ke_list, edges, enrichment, ke_title_map, ke_type_map)

        assert 'nodes' in result
        assert 'edges' in result
        assert len(result['nodes']) == 2
        assert len(result['edges']) == 1

    def test_empty_edges(self):
        ke_list = {'KE:1'}
        edges = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])
        enrichment = self._make_enrichment(['KE:1'], [0.01])

        result = build_cytoscape_network(ke_list, edges, enrichment, {'KE:1': 'E1'}, {'KE:1': 'MIE'})

        assert len(result['nodes']) == 1
        assert len(result['edges']) == 0

    def test_ke_without_enrichment(self):
        """KE in ke_list but not in enrichment results should still appear as a node."""
        ke_list = {'KE:1', 'KE:2'}
        edges = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])
        enrichment = self._make_enrichment(['KE:1'], [0.01])  # only KE:1

        result = build_cytoscape_network(ke_list, edges, enrichment,
                                         {'KE:1': 'E1', 'KE:2': 'E2'},
                                         {'KE:1': 'MIE', 'KE:2': 'AO'})

        assert len(result['nodes']) == 2
        # KE:2 should have non-significant default
        ke2_node = [n for n in result['nodes'] if n['data']['id'] == 'KE:2'][0]
        assert ke2_node['data'].get('significant') is False or not ke2_node.get('classes', '').count('significant')

    def test_node_data_includes_pvalue_and_fdr(self):
        """D-07: p_value and FDR are echoed onto each KE node's data payload
        so the existing Download Network (.cyjs) carries them automatically."""
        ke_list = {'KE:1', 'KE:2'}
        edges = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])
        enrichment = self._make_enrichment(['KE:1', 'KE:2'], [0.001, 0.5])
        result = build_cytoscape_network(
            ke_list, edges, enrichment,
            {'KE:1': 'E1', 'KE:2': 'E2'},
            {'KE:1': 'MIE', 'KE:2': 'AO'},
        )
        by_id = {n['data']['id']: n['data'] for n in result['nodes']}
        assert by_id['KE:1']['p_value'] == pytest.approx(0.001)
        assert by_id['KE:1']['fdr'] == pytest.approx(0.001)
        assert by_id['KE:2']['p_value'] == pytest.approx(0.5)
        assert by_id['KE:2']['fdr'] == pytest.approx(0.5)
        # Existing fields untouched (regression guard):
        assert by_id['KE:1']['logfc'] == 2.0  # odds_ratio surrogate, unchanged
        assert by_id['KE:1']['ke_type'] == 'MIE'
        assert by_id['KE:1']['label'] == 'E1'

    def test_node_data_pvalue_is_none_when_no_enrichment(self):
        """KE in ke_list but absent from enrichment_results gets None
        (keys present, values None — schema-stable for downstream consumers)."""
        ke_list = {'KE:1', 'KE:2'}
        enrichment = self._make_enrichment(['KE:1'], [0.01])  # KE:2 missing
        result = build_cytoscape_network(
            ke_list,
            pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID']),
            enrichment,
            {'KE:1': 'E1', 'KE:2': 'E2'},
            {'KE:1': 'MIE', 'KE:2': 'AO'},
        )
        ke2 = next(n for n in result['nodes'] if n['data']['id'] == 'KE:2')
        assert 'p_value' in ke2['data']
        assert 'fdr' in ke2['data']
        assert ke2['data']['p_value'] is None
        assert ke2['data']['fdr'] is None

    def test_node_data_nan_sanitised_to_none(self):
        """NaN values from enrichment must become None — json.dumps must succeed."""
        ke_list = {'KE:1'}
        enrichment = self._make_enrichment(['KE:1'], [float('nan')])
        result = build_cytoscape_network(
            ke_list,
            pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID']),
            enrichment,
            {'KE:1': 'E1'},
            {'KE:1': 'MIE'},
        )
        ke1 = next(n for n in result['nodes'] if n['data']['id'] == 'KE:1')
        assert ke1['data']['fdr'] is None
        # Critical: the entire dict must json-serialise (no NaN leaks)
        serialised = json.dumps(result)
        assert 'NaN' not in serialised

    def test_node_data_numpy_scalars_become_native(self):
        """numpy scalars (np.float64) must be unwrapped to native Python floats."""
        ke_list = {'KE:1'}
        enrichment = pd.DataFrame({
            'KE': ['KE:1'],
            'KE_Title': ['Event 1'],
            'p_value': np.array([np.float64(0.0042)]),
            'FDR': np.array([np.float64(0.0084)]),
            'odds_ratio': [2.0],
        })
        result = build_cytoscape_network(
            ke_list,
            pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID']),
            enrichment,
            {'KE:1': 'E1'},
            {'KE:1': 'MIE'},
        )
        ke1 = next(n for n in result['nodes'] if n['data']['id'] == 'KE:1')
        # Must be a native python float, not numpy scalar.
        assert isinstance(ke1['data']['p_value'], float)
        assert ke1['data']['p_value'] == pytest.approx(0.0042)
        assert isinstance(ke1['data']['fdr'], float)


class TestBuildComparisonMatrix:
    """Tests for cross-condition comparison matrix."""

    def _make_condition(self, label, enrichment_data, position=0):
        """Create a mock ConditionRecord-like object."""
        return SimpleNamespace(
            condition_label=label,
            enrichment_json=json.dumps(enrichment_data),
            position=position,
        )

    def test_basic_matrix(self):
        c1 = self._make_condition('Cond A', [
            {'KE': 'KE:1', 'FDR': 0.01, 'Title': 'Event 1'},
            {'KE': 'KE:2', 'FDR': 0.1, 'Title': 'Event 2'},
        ])
        c2 = self._make_condition('Cond B', [
            {'KE': 'KE:1', 'FDR': 0.05, 'Title': 'Event 1'},
            {'KE': 'KE:2', 'FDR': 0.001, 'Title': 'Event 2'},
        ])
        result = build_comparison_matrix([c1, c2])

        assert len(result['ke_labels']) == 2
        assert len(result['condition_labels']) == 2
        assert result['condition_labels'] == ['Cond A', 'Cond B']
        assert len(result['fdr_matrix']) == 2  # 2 KEs
        assert len(result['fdr_matrix'][0]) == 2  # 2 conditions per KE
        assert len(result['condition_colors']) == 2

    def test_empty_conditions(self):
        result = build_comparison_matrix([])
        assert result == {}

    def test_no_enrichment_data(self):
        c1 = self._make_condition('Empty', [])
        result = build_comparison_matrix([c1])
        assert result == {}

    def test_condition_colors_from_palette(self):
        conditions = [
            self._make_condition(f'C{i}', [{'KE': 'KE:1', 'FDR': 0.01, 'Title': 'E1'}])
            for i in range(5)
        ]
        result = build_comparison_matrix(conditions)
        assert result['condition_colors'] == CONDITION_PALETTE[:5]

    def test_nonsignificant_fdr_becomes_none_in_neg_log10(self):
        """FDR > 0.05 should produce None in the -log10 matrix."""
        c1 = self._make_condition('C1', [
            {'KE': 'KE:1', 'FDR': 0.5, 'Title': 'Non-sig'},
        ])
        result = build_comparison_matrix([c1])
        # The neg_log10 value for FDR=0.5 (non-significant) should be None
        assert result['neg_log10_matrix'][0][0] is None

    def test_malformed_json_skipped(self):
        """Condition with invalid JSON is gracefully skipped."""
        c1 = SimpleNamespace(condition_label='Bad', enrichment_json='not-json', position=0)
        c2 = self._make_condition('Good', [{'KE': 'KE:1', 'FDR': 0.01, 'Title': 'E1'}])
        result = build_comparison_matrix([c1, c2])
        assert len(result['condition_labels']) == 2
        assert len(result['ke_labels']) == 1


class TestNetworkSignificanceIsFdrDriven:
    """Issue #63: the red border is the BH-adjusted FDR, not the raw p-value."""

    def _enrichment(self, ke_ids, p_values, fdrs):
        return pd.DataFrame({
            'KE': ke_ids,
            'p_value': p_values,
            'FDR': fdrs,
            'odds_ratio': [2.0] * len(ke_ids),
            'num_overlap': [5] * len(ke_ids),
            'total_KE_genes_in_dataset': [20] * len(ke_ids),
        })

    def _build(self, enrichment, ke_list, **kwargs):
        return build_cytoscape_network(
            ke_list,
            pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID']),
            enrichment,
            {ke: ke for ke in ke_list},
            {ke: 'intermediate' for ke in ke_list},
            **kwargs,
        )

    def test_raw_p_below_but_fdr_above_cutoff_is_not_significant(self):
        """The regression this issue exists for: p=0.01, FDR=0.30 -> no red border."""
        enrichment = self._enrichment(['KE:1'], [0.01], [0.30])
        result = self._build(enrichment, {'KE:1'})
        node = result['nodes'][0]
        assert 'significant' not in node['classes']
        # The raw p stays visible in the payload — reported, never gating.
        assert node['data']['p_value'] == pytest.approx(0.01)
        assert node['data']['fdr'] == pytest.approx(0.30)

    def test_fdr_below_cutoff_is_significant(self):
        enrichment = self._enrichment(['KE:1'], [0.001], [0.02])
        node = self._build(enrichment, {'KE:1'})['nodes'][0]
        assert 'significant' in node['classes']

    def test_cutoff_is_a_parameter(self):
        """The cutoff is threadable, not a hardcoded literal."""
        enrichment = self._enrichment(['KE:1'], [0.01], [0.30])
        assert 'significant' in self._build(
            enrichment, {'KE:1'}, fdr_cutoff=0.5)['nodes'][0]['classes']
        assert 'significant' not in self._build(
            enrichment, {'KE:1'}, fdr_cutoff=0.01)['nodes'][0]['classes']

    def test_default_cutoff_matches_config(self):
        """Default is Config.SIGNIFICANCE_FDR_CUTOFF — the one app-wide cutoff."""
        assert Config.SIGNIFICANCE_FDR_CUTOFF == 0.05
        # Just above / just below the configured cutoff.
        enrichment = self._enrichment(['KE:1', 'KE:2'], [0.001, 0.001], [0.049, 0.051])
        by_id = {n['data']['id']: n for n in self._build(enrichment, {'KE:1', 'KE:2'})['nodes']}
        assert 'significant' in by_id['KE:1']['classes']
        assert 'significant' not in by_id['KE:2']['classes']

    def test_network_agrees_with_comparison_matrix(self):
        """Issue #63: the batch heatmap and the batch network in the SAME PDF
        must call the same Key Events significant."""
        rows = [
            {'KE': 'KE:1', 'Title': 'Alpha', 'p_value': 0.001, 'FDR': 0.01},   # sig both
            {'KE': 'KE:2', 'Title': 'Beta', 'p_value': 0.01, 'FDR': 0.30},     # raw p only
        ]
        enrichment = pd.DataFrame(rows)
        enrichment['odds_ratio'] = 2.0
        network = self._build(enrichment, {'KE:1', 'KE:2'})
        net_significant = {
            n['data']['id'] for n in network['nodes'] if 'significant' in n['classes']
        }

        cond = SimpleNamespace(condition_label='C1', enrichment_json=json.dumps(rows), position=0)
        matrix = build_comparison_matrix([cond])
        matrix_significant = {
            ke for ke, row in zip(matrix['ke_labels'], matrix['fdr_matrix'])
            if row[0] is not None and row[0] < 0.05
        }

        assert net_significant == matrix_significant == {'KE:1'}


class TestExcludedKeStyling:
    """Issue #65: untested KEs must not render like tested-but-unenriched ones."""

    def _build(self, ke_list, excluded_kes=None, reference_sets=None):
        enrichment = pd.DataFrame({
            'KE': ['KE:OK'], 'p_value': [0.001], 'FDR': [0.01], 'odds_ratio': [2.0],
        })
        return build_cytoscape_network(
            ke_list,
            pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID']),
            enrichment,
            {ke: ke for ke in ke_list},
            {ke: 'intermediate' for ke in ke_list},
            reference_sets=reference_sets,
            excluded_kes=excluded_kes,
        )

    def test_too_few_genes_gets_its_own_class(self):
        result = self._build(
            {'KE:OK', 'KE:SMALL'}, excluded_kes={'KE:SMALL': 'too_few_genes'}
        )
        by_id = {n['data']['id']: n for n in result['nodes']}
        assert 'too-few-genes' in by_id['KE:SMALL']['classes']
        assert by_id['KE:SMALL']['data']['excluded_reason'] == 'too_few_genes'
        # The tested KE must not pick up any exclusion styling.
        assert 'too-few-genes' not in by_id['KE:OK']['classes']
        assert by_id['KE:OK']['data']['excluded_reason'] is None

    def test_no_mapping_gets_the_existing_no_genes_class(self):
        result = self._build(
            {'KE:OK', 'KE:MISS'}, excluded_kes={'KE:MISS': 'no_mapping'}
        )
        by_id = {n['data']['id']: n for n in result['nodes']}
        assert 'no-genes' in by_id['KE:MISS']['classes']
        assert 'too-few-genes' not in by_id['KE:MISS']['classes']

    def test_no_genes_class_not_duplicated(self):
        """A KE with neither a gene set nor a mapping gets one class, not two."""
        result = self._build(
            {'KE:MISS'}, excluded_kes={'KE:MISS': 'no_mapping'}, reference_sets={}
        )
        assert result['nodes'][0]['classes'].split().count('no-genes') == 1

    def test_excluded_kes_optional(self):
        """Callers that pass nothing get the pre-issue payload shape (plus a
        None excluded_reason) — no crash, no styling."""
        result = self._build({'KE:OK'})
        assert result['nodes'][0]['classes'] == 'significant'
        assert result['nodes'][0]['data']['excluded_reason'] is None


class TestKeAccountingFromNetwork:
    """Issue #65: recover the accounting from a stored network payload."""

    def _network(self, reasons):
        return {
            'nodes': [
                {'data': {'id': ke, 'ke_type': 'intermediate', 'excluded_reason': r},
                 'classes': ''}
                for ke, r in reasons.items()
            ],
            'edges': [],
        }

    def test_counts_recovered_from_nodes(self):
        summary = ke_accounting_from_network(self._network({
            'KE:1': None, 'KE:2': None, 'KE:3': 'too_few_genes', 'KE:4': 'no_mapping',
        }))
        assert summary['total_kes'] == 4
        assert summary['tested'] == 2
        assert summary['excluded_too_few_genes'] == 1
        assert summary['excluded_no_mapping'] == 1

    def test_accepts_json_string(self):
        payload = json.dumps(self._network({'KE:1': None, 'KE:2': 'too_few_genes'}))
        assert ke_accounting_from_network(payload)['tested'] == 1

    def test_returns_none_for_legacy_networks(self):
        """Networks stored before the field existed yield None (line omitted)."""
        legacy = {'nodes': [{'data': {'id': 'KE:1', 'ke_type': 'MIE'}, 'classes': ''}],
                  'edges': []}
        assert ke_accounting_from_network(legacy) is None
        assert ke_accounting_from_network(None) is None
        assert ke_accounting_from_network('not-json') is None
