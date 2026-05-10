"""Tests for network and comparison services."""
import json
import math
import numpy as np
import pytest
import pandas as pd
from types import SimpleNamespace

from services.network_service import build_cytoscape_network
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
