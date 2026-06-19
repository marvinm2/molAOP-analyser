"""Tests for comparison export shaping and cross-condition gene tracking.

Covers (plan T-3/T-4/T-5):
- build_comparison_matrix new per-condition metadata + edge cases (H-1/H-2)
- comparison_matrix_to_dataframe wide shaping + CSV-injection guard
- build_gene_tracking / gene_tracking_to_dataframes shared vs specific drivers
"""
import json
from types import SimpleNamespace

import pandas as pd
import pytest

from services.comparison_service import (
    build_comparison_matrix,
    comparison_matrix_to_dataframe,
    build_gene_tracking,
    gene_tracking_to_dataframes,
    csv_guard,
)


def _condition(label, enrichment, ke_gene=None, position=0, **extra):
    """Build a ConditionRecord-like stub."""
    return SimpleNamespace(
        condition_label=label,
        enrichment_json=json.dumps(enrichment) if enrichment is not None else None,
        ke_gene_json=json.dumps(ke_gene) if ke_gene is not None else None,
        position=position,
        gene_count=extra.get('gene_count'),
        significant_genes=extra.get('significant_genes'),
        dose=extra.get('dose', ''),
        timepoint=extra.get('timepoint', ''),
    )


class TestComparisonMatrixMetadata:
    """T-3: new per-condition metadata keys + edge cases."""

    def test_metadata_keys_present_and_ordered(self):
        c1 = _condition('C1', [{'KE': 'KE:1', 'FDR': 0.001, 'Title': 'E1'}],
                        gene_count=100, significant_genes=20, dose='1uM', timepoint='4hr')
        c2 = _condition('C2', [{'KE': 'KE:1', 'FDR': 0.5, 'Title': 'E1'}],
                        gene_count=120, significant_genes=40, dose='10uM', timepoint='4hr', position=1)
        m = build_comparison_matrix([c1, c2])

        assert m['condition_gene_counts'] == [100, 120]
        assert m['condition_sig_gene_counts'] == [20, 40]
        assert m['condition_doses'] == ['1uM', '10uM']
        assert m['condition_timepoints'] == ['4hr', '4hr']
        # KE:1 is significant only in C1 (FDR 0.001 < 0.05), not C2 (0.5).
        assert m['condition_sig_ke_counts'] == [1, 0]

    def test_empty_condition_keeps_column_aligned(self):
        """H-2: a condition with no enrichment still occupies its column."""
        c1 = _condition('C1', [{'KE': 'KE:1', 'FDR': 0.001, 'Title': 'E1'}],
                        gene_count=100, significant_genes=20)
        c2 = _condition('C2', [], gene_count=0, significant_genes=0, position=1)
        m = build_comparison_matrix([c1, c2])

        assert m['condition_labels'] == ['C1', 'C2']
        assert len(m['condition_colors']) == 2
        # Every KE row must have exactly 2 columns (C2 back-filled with None).
        for row in m['fdr_matrix']:
            assert len(row) == 2
        # C2's cell is absent.
        assert m['fdr_matrix'][0][1] is None

    def test_duplicate_ke_does_not_crash(self):
        """H-1: duplicate (condition, KE) rows are tolerated (warning, no crash)."""
        c1 = _condition('C1', [
            {'KE': 'KE:1', 'FDR': 0.001, 'Title': 'E1'},
            {'KE': 'KE:1', 'FDR': 0.002, 'Title': 'E1 dup'},
        ])
        m = build_comparison_matrix([c1])
        assert m['ke_labels'] == ['KE:1']
        assert len(m['fdr_matrix']) == 1


class TestComparisonMatrixToDataFrame:
    """T-4: wide CSV/Excel shaping."""

    def _matrix(self):
        c1 = _condition('C1', [
            {'KE': 'KE:1', 'FDR': 0.001, 'Title': '=Evil'},
            {'KE': 'KE:2', 'FDR': 0.5, 'Title': 'Beta'},
        ])
        c2 = _condition('C2', [{'KE': 'KE:1', 'FDR': 0.0001, 'Title': 'Alpha'}], position=1)
        return build_comparison_matrix([c1, c2])

    def test_fdr_wide_shape(self):
        df = comparison_matrix_to_dataframe(self._matrix(), which='fdr')
        assert list(df.columns) == ['Key Event ID', 'Key Event Title', 'C1', 'C2']
        # KE:1 row present with both condition FDRs.
        ke1 = df[df['Key Event ID'] == 'KE:1'].iloc[0]
        assert ke1['C1'] == pytest.approx(0.001)
        assert ke1['C2'] == pytest.approx(0.0001)

    def test_neglog10_blanks_nonsignificant(self):
        df = comparison_matrix_to_dataframe(self._matrix(), which='neglog10')
        ke2 = df[df['Key Event ID'] == 'KE:2'].iloc[0]
        # FDR 0.5 is non-significant -> blank (NaN) in neglog10.
        assert pd.isna(ke2['C1'])

    def test_csv_injection_guard_on_title(self):
        df = comparison_matrix_to_dataframe(self._matrix(), which='fdr')
        titles = df['Key Event Title'].tolist()
        assert any(t.startswith("'=") for t in titles)

    def test_empty_matrix_returns_empty_df(self):
        assert comparison_matrix_to_dataframe({}, which='fdr').empty

    def test_invalid_selector_raises(self):
        with pytest.raises(ValueError):
            comparison_matrix_to_dataframe(self._matrix(), which='bogus')


class TestGeneTracking:
    """T-5: cross-condition driver-gene tracking."""

    def _conditions(self):
        c1 = _condition(
            'C1',
            [{'KE': 'KE:1', 'Title': 'Alpha', 'FDR': 0.001}],
            ke_gene={'KE:1': [
                {'id': 'TP53', 'log2FC': 2.1, 'significant': True},
                {'id': 'EGFR', 'log2FC': 0.1, 'significant': False},
                {'id': 'MYC', 'log2FC': 1.5, 'significant': True},
            ]},
        )
        c2 = _condition(
            'C2',
            [{'KE': 'KE:1', 'Title': 'Alpha', 'FDR': 0.0001}],
            ke_gene={'KE:1': [
                {'id': 'TP53', 'log2FC': 3.0, 'significant': True},
            ]},
            position=1,
        )
        return [c1, c2]

    def test_records_only_significant_drivers(self):
        t = build_gene_tracking(self._conditions())
        genes = {(r['Gene_Symbol'], r['Condition']) for r in t['records']}
        assert ('TP53', 'C1') in genes
        assert ('TP53', 'C2') in genes
        assert ('MYC', 'C1') in genes
        # EGFR is not significant -> never a driver.
        assert all(r['Gene_Symbol'] != 'EGFR' for r in t['records'])

    def test_shared_vs_specific_flags(self):
        t = build_gene_tracking(self._conditions())
        by_gene = {r['Gene_Symbol']: r for r in t['gene_ke_summary']}
        assert by_gene['TP53']['shared'] is True
        assert by_gene['TP53']['n_conditions_driver'] == 2
        assert by_gene['MYC']['shared'] is False
        assert by_gene['MYC']['n_conditions_driver'] == 1

    def test_log2fc_by_condition_alignment(self):
        t = build_gene_tracking(self._conditions())
        by_gene = {r['Gene_Symbol']: r for r in t['gene_ke_summary']}
        assert by_gene['TP53']['log2FC_by_condition'] == [2.1, 3.0]
        # MYC only in C1 -> second slot None.
        assert by_gene['MYC']['log2FC_by_condition'] == [1.5, None]

    def test_to_dataframes_columns(self):
        t = build_gene_tracking(self._conditions())
        frames = gene_tracking_to_dataframes(t)
        assert list(frames['long'].columns) == [
            'KE_ID', 'KE_Title', 'Gene_Symbol', 'Condition', 'log2FC', 'significant'
        ]
        assert 'log2FC: C1' in frames['summary'].columns
        assert 'log2FC: C2' in frames['summary'].columns


def test_csv_guard():
    assert csv_guard('=cmd') == "'=cmd"
    assert csv_guard('+1') == "'+1"
    assert csv_guard('normal') == 'normal'
    assert csv_guard(0.5) == 0.5
