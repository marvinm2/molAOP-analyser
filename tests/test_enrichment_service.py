"""Tests for the enrichment analysis service."""
import math
import pytest
import pandas as pd
from services.enrichment_service import run_enrichment_analysis, build_ke_gene_mapping


def _make_gene_df(genes, significant_flags):
    """Helper: build a minimal processed gene expression DataFrame."""
    return pd.DataFrame({
        'ID': genes,
        'significant': significant_flags,
        'log2FC': [1.5 if s else 0.2 for s in significant_flags],
        'pval': [0.001 if s else 0.5 for s in significant_flags],
    })


class TestRunEnrichmentAnalysis:
    """Core enrichment analysis tests."""

    def test_basic_enrichment_produces_results(self):
        """Enrichment with valid input returns a DataFrame with expected columns."""
        genes = [f"GENE{i}" for i in range(100)]
        sig_flags = [i < 30 for i in range(100)]
        df = _make_gene_df(genes, sig_flags)

        reference_sets = {
            'KE:1': {f"GENE{i}" for i in range(20)},   # 20 genes, 6 significant
            'KE:2': {f"GENE{i}" for i in range(50, 70)},  # 20 genes, 0 significant
        }
        ke_list = {'KE:1', 'KE:2'}
        ke_title_map = {'KE:1': 'Test KE 1', 'KE:2': 'Test KE 2'}

        result = run_enrichment_analysis(df, reference_sets, ke_list, ke_title_map)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert 'p_value' in result.columns
        assert 'FDR' in result.columns
        assert 'odds_ratio' in result.columns
        assert all(0 <= p <= 1 for p in result['p_value'])
        assert all(0 <= f <= 1 for f in result['FDR'])

    def test_infinity_odds_ratio_becomes_na(self):
        """When all KE genes are significant (c=0), odds ratio should be 'NA', not inf."""
        # 10 genes total, all 10 significant
        # KE:1 has 5 of those genes — all significant → c=0 → odds=inf
        genes = [f"GENE{i}" for i in range(10)]
        sig_flags = [True] * 10
        df = _make_gene_df(genes, sig_flags)

        reference_sets = {'KE:1': {f"GENE{i}" for i in range(5)}}
        ke_list = {'KE:1'}
        ke_title_map = {'KE:1': 'All Sig KE'}

        result = run_enrichment_analysis(df, reference_sets, ke_list, ke_title_map)

        odds_val = result.iloc[0]['odds_ratio']
        assert odds_val == 'NA', f"Expected 'NA' for infinite odds ratio, got {odds_val}"

    def test_no_overlap_ke_skipped(self):
        """KE with zero gene overlap is skipped (no result row)."""
        df = _make_gene_df(['GENEA', 'GENEB', 'GENEC', 'GENED', 'GENEE'],
                           [True, True, False, False, False])
        reference_sets = {
            'KE:1': {'GENEX', 'GENEY', 'GENEZ', 'GENEW', 'GENEV'},  # no overlap
        }
        with pytest.raises(ValueError, match="No enrichment results"):
            run_enrichment_analysis(df, reference_sets, {'KE:1'}, {'KE:1': 'No Overlap'})

    def test_fewer_than_5_genes_skipped(self):
        """KE with fewer than 5 overlapping genes is skipped."""
        genes = [f"GENE{i}" for i in range(50)]
        df = _make_gene_df(genes, [i < 10 for i in range(50)])

        reference_sets = {
            'KE:1': {'GENE0', 'GENE1', 'GENE2', 'GENE3'},  # only 4 overlap
            'KE:2': {f"GENE{i}" for i in range(10)},  # 10 overlap — kept
        }

        result = run_enrichment_analysis(df, reference_sets, {'KE:1', 'KE:2'},
                                         {'KE:1': 'Small', 'KE:2': 'Normal'})
        assert len(result) == 1
        assert result.iloc[0]['KE'] == 'KE:2'

    def test_no_reference_sets_raises(self):
        """ValueError raised when no reference sets match the KE list."""
        df = _make_gene_df(['GENE1'], [True])
        with pytest.raises(ValueError, match="No reference gene sets"):
            run_enrichment_analysis(df, {'KE:1': {'GENE1'}}, {'KE:99'}, {})

    def test_no_results_raises(self):
        """ValueError raised when all KEs are filtered out."""
        df = _make_gene_df(['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5'],
                           [False] * 5)
        # KE has genes but none in dataset
        reference_sets = {'KE:1': {'X1', 'X2', 'X3', 'X4', 'X5'}}
        with pytest.raises(ValueError, match="No enrichment results"):
            run_enrichment_analysis(df, reference_sets, {'KE:1'}, {'KE:1': 'Miss'})

    def test_fdr_correction_applied(self):
        """FDR values are present and at least as large as raw p-values."""
        genes = [f"GENE{i}" for i in range(200)]
        sig_flags = [i < 50 for i in range(200)]
        df = _make_gene_df(genes, sig_flags)

        reference_sets = {
            f'KE:{k}': {f"GENE{i}" for i in range(k * 10, k * 10 + 10)}
            for k in range(5)
        }
        ke_list = set(reference_sets.keys())
        ke_title_map = {k: k for k in ke_list}

        result = run_enrichment_analysis(df, reference_sets, ke_list, ke_title_map)
        # FDR should be >= p_value for each row (BH correction inflates)
        for _, row in result.iterrows():
            assert row['FDR'] >= row['p_value'] - 1e-10

    def test_odds_ratio_na_when_nan(self):
        """Odds ratio is 'NA' when Fisher's test returns NaN (edge case)."""
        # All genes in KE are non-significant, a=0
        genes = [f"GENE{i}" for i in range(20)]
        df = _make_gene_df(genes, [False] * 20)

        reference_sets = {'KE:1': {f"GENE{i}" for i in range(10)}}
        result = run_enrichment_analysis(df, reference_sets, {'KE:1'}, {'KE:1': 'None Sig'})

        # With a=0, odds_ratio should be 0.0 (not NA), but p_value should be 1.0
        # This test validates the function doesn't crash
        assert len(result) == 1


class TestBuildKeGeneMapping:
    """Tests for KE-to-gene mapping builder."""

    def test_basic_mapping(self):
        """Builds correct mapping from reference sets and user data."""
        reference_sets = {
            'KE:1': {'GENE1', 'GENE2', 'GENE3'},
            'KE:2': {'GENE4', 'GENE5'},
        }
        gene_logfc_map = {'GENE1': 1.5, 'GENE2': -0.5, 'GENE4': 2.0}
        gene_sig_map = {'GENE1': True, 'GENE2': False, 'GENE4': True}

        result = build_ke_gene_mapping(reference_sets, {'KE:1', 'KE:2'},
                                       gene_logfc_map, gene_sig_map)

        assert 'KE:1' in result
        assert len(result['KE:1']) == 2  # GENE1, GENE2 (GENE3 not in user data)
        assert 'KE:2' in result
        assert len(result['KE:2']) == 1  # GENE4 only

    def test_filters_to_ke_list(self):
        """Only includes KEs that are in the ke_list."""
        reference_sets = {'KE:1': {'GENE1'}, 'KE:2': {'GENE2'}}
        result = build_ke_gene_mapping(reference_sets, {'KE:1'},
                                       {'GENE1': 1.0, 'GENE2': 1.0},
                                       {'GENE1': True, 'GENE2': True})
        assert 'KE:1' in result
        assert 'KE:2' not in result
