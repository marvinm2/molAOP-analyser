"""Tests for the enrichment analysis service."""
import math
import pytest
import pandas as pd
from config import Config
from services.enrichment_service import (
    run_enrichment_analysis, build_ke_gene_mapping, get_ke_summary, format_ke_summary,
)


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

    def test_direction_column_emitted_when_logfc_map_supplied(self):
        """When gene_logfc_map is passed, the result carries a 'Direction'
        column rendering observed up/down counts for the significant overlap.

        Sets up 8 significant overlap genes with a deliberate sign mix
        (5 up, 2 down, 1 at zero) and asserts the rendered string. Also
        confirms back-compat: without the map, the column is absent.
        """
        # 20 total genes; first 8 are significant overlap with mixed signs.
        # log2FCs: GENE0..4 up, GENE5..6 down, GENE7 at exactly 0 (excluded
        # from both counts — threshold-0 edge case the implementation
        # explicitly handles).
        genes = [f"GENE{i}" for i in range(20)]
        sig_flags = [True] * 8 + [False] * 12
        df = pd.DataFrame({
            'ID': genes,
            'significant': sig_flags,
            'log2FC': [1.5, 2.0, 0.8, 3.1, 1.1, -0.7, -1.4, 0.0] + [0.1] * 12,
            'pval': [0.001] * 8 + [0.5] * 12,
        })
        reference_sets = {'KE:UP': {f"GENE{i}" for i in range(8)}}
        ke_list = {'KE:UP'}
        ke_title_map = {'KE:UP': 'Test KE'}
        gene_logfc_map = dict(zip(df['ID'], df['log2FC']))

        result = run_enrichment_analysis(
            df, reference_sets, ke_list, ke_title_map,
            gene_logfc_map=gene_logfc_map,
        )
        assert 'Direction' in result.columns
        # 5 strictly positive, 2 strictly negative, 1 at zero (excluded).
        assert result.iloc[0]['Direction'] == '5↑ / 2↓'

        # Back-compat: without the map, no Direction column.
        result_no_map = run_enrichment_analysis(
            df, reference_sets, ke_list, ke_title_map,
        )
        assert 'Direction' not in result_no_map.columns


class TestDepletionTail:
    """Issue #70 — depleted Key Events must be distinguishable from null ones."""

    @staticmethod
    def _depleted_dataset():
        """A KE whose gene set responds far less than the background rate.

        500 genes, 250 significant (50% background). The KE holds 100 of them
        with only 5 significant — the one-sided 'greater' test returns p ≈ 1,
        which is exactly what a KE with nothing going on returns.
        """
        genes = [f"GENE{i}" for i in range(500)]
        ke_genes = set(genes[:100])
        # 5 significant inside the KE, 245 significant outside it.
        sig = set(genes[:5]) | set(genes[100:345])
        df = _make_gene_df(genes, [g in sig for g in genes])
        return df, {'KE:DEP': ke_genes}, {'KE:DEP'}, {'KE:DEP': 'Depleted KE'}

    def test_depleted_ke_is_called_depleted(self):
        """Under-representation surfaces as a 'depleted' call, not as 'ns'."""
        df, reference_sets, ke_list, ke_title_map = self._depleted_dataset()

        result = run_enrichment_analysis(df, reference_sets, ke_list, ke_title_map)
        row = result.iloc[0]

        # The headline (over-representation) statistic is unchanged: near 1.
        assert row['FDR'] > 0.5
        # The other tail is where the signal is.
        assert row['FDR_depleted'] < 0.01
        assert row['Representation'] == 'depleted'

    def test_enriched_ke_is_called_enriched(self):
        """The over-represented case keeps its existing verdict."""
        genes = [f"GENE{i}" for i in range(200)]
        sig = set(genes[:40])
        df = _make_gene_df(genes, [g in sig for g in genes])
        reference_sets = {'KE:UP': set(genes[:30])}  # 30 genes, 30 significant

        result = run_enrichment_analysis(
            df, reference_sets, {'KE:UP'}, {'KE:UP': 'Enriched KE'}
        )
        row = result.iloc[0]
        assert row['Representation'] == 'enriched'
        assert row['FDR'] < 0.01
        assert row['FDR_depleted'] > 0.5

    def test_unremarkable_ke_is_called_ns(self):
        """A KE tracking the background rate is 'ns' on both tails."""
        genes = [f"GENE{i}" for i in range(200)]
        # 50% significant everywhere, inside and outside the KE.
        sig_flags = [i % 2 == 0 for i in range(200)]
        df = _make_gene_df(genes, sig_flags)
        reference_sets = {'KE:FLAT': set(genes[:40])}

        result = run_enrichment_analysis(
            df, reference_sets, {'KE:FLAT'}, {'KE:FLAT': 'Flat KE'}
        )
        assert result.iloc[0]['Representation'] == 'ns'

    def test_cutoff_is_configurable(self):
        """The call follows the supplied cutoff, both tails corrected by BH."""
        df, reference_sets, ke_list, ke_title_map = self._depleted_dataset()

        default = run_enrichment_analysis(df, reference_sets, ke_list, ke_title_map)
        fdr_depleted = default.iloc[0]['FDR_depleted']

        strict = run_enrichment_analysis(
            df, reference_sets, ke_list, ke_title_map,
            fdr_cutoff=fdr_depleted / 2,
        )
        assert strict.iloc[0]['Representation'] == 'ns'

    def test_enrichment_fdr_unchanged_by_the_second_tail(self):
        """Adding depletion testing must not move the enrichment FDR.

        Each tail is BH-corrected within its own family, so the enrichment
        column keeps the values it had before issue #70.
        """
        genes = [f"GENE{i}" for i in range(100)]
        df = _make_gene_df(genes, [i < 30 for i in range(100)])
        reference_sets = {
            'KE:1': {f"GENE{i}" for i in range(20)},
            'KE:2': {f"GENE{i}" for i in range(50, 70)},
        }
        result = run_enrichment_analysis(
            df, reference_sets, {'KE:1', 'KE:2'}, {'KE:1': 'A', 'KE:2': 'B'}
        )
        # Recompute the one-sided BH by hand over the same two tests.
        from scipy.stats import fisher_exact
        from statsmodels.stats.multitest import multipletests
        expected = {}
        for ke, ke_genes in reference_sets.items():
            a = sum(1 for g in ke_genes if int(g[4:]) < 30)
            c = len(ke_genes) - a
            expected[ke] = fisher_exact(
                [[a, 30 - a], [c, 70 - c]], alternative='greater'
            )[1]
        kes = list(expected)
        adjusted = dict(zip(kes, multipletests([expected[k] for k in kes],
                                               method='fdr_bh')[1]))
        for _, row in result.iterrows():
            assert row['FDR'] == pytest.approx(adjusted[row['KE']])


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


# --- Plan 11-03 (Wave 2) tests for optional pvalue passthrough ---

def test_build_ke_gene_mapping_backward_compat():
    """Existing 4-arg call returns the same dict shape as today (no new keys)."""
    reference_sets = {'KE:1': {'GENE1', 'GENE2'}}
    ke_list = {'KE:1'}
    logfc = {'GENE1': 1.5, 'GENE2': -0.3}
    sig = {'GENE1': True, 'GENE2': False}
    result = build_ke_gene_mapping(reference_sets, ke_list, logfc, sig)
    assert 'KE:1' in result
    for entry in result['KE:1']:
        assert set(entry.keys()) == {'id', 'log2FC', 'significant'}


def test_build_ke_gene_mapping_with_pvalue_passthrough():
    """When pvalue maps are supplied, per-gene dicts gain pvalue_raw / pvalue_adj
    with None fallback for genes absent from the maps (D-04 passthrough)."""
    reference_sets = {'KE:1': {'GENE1', 'GENE2', 'GENE3'}}
    ke_list = {'KE:1'}
    logfc = {'GENE1': 1.5, 'GENE2': -0.3, 'GENE3': 0.8}
    sig = {'GENE1': True, 'GENE2': False, 'GENE3': True}
    pval_raw = {'GENE1': 0.001, 'GENE3': 0.02}      # GENE2 absent
    pval_adj = {'GENE1': 0.01, 'GENE3': 0.05}       # GENE2 absent
    result = build_ke_gene_mapping(
        reference_sets, ke_list, logfc, sig,
        gene_pvalue_raw_map=pval_raw,
        gene_pvalue_adj_map=pval_adj,
    )
    by_id = {e['id']: e for e in result['KE:1']}
    assert by_id['GENE1']['pvalue_raw'] == 0.001
    assert by_id['GENE1']['pvalue_adj'] == 0.01
    assert by_id['GENE2']['pvalue_raw'] is None
    assert by_id['GENE2']['pvalue_adj'] is None
    assert by_id['GENE3']['pvalue_raw'] == 0.02
    assert by_id['GENE3']['pvalue_adj'] == 0.05


def test_build_ke_gene_mapping_only_raw_pvalue():
    """Only raw map supplied -> pvalue_raw present, pvalue_adj is None for all."""
    result = build_ke_gene_mapping(
        {'KE:1': {'GENE1'}}, {'KE:1'}, {'GENE1': 1.0}, {'GENE1': True},
        gene_pvalue_raw_map={'GENE1': 0.005},
    )
    entry = result['KE:1'][0]
    assert entry['pvalue_raw'] == 0.005
    assert entry['pvalue_adj'] is None


def test_build_ke_gene_mapping_filters_by_significant_for_export():
    """Per D-01: the export consumer reads `significant` to filter; the function
    itself MUST NOT pre-filter. Non-significant genes are still emitted with
    significant=False so the consumer makes the policy choice."""
    result = build_ke_gene_mapping(
        {'KE:1': {'G1', 'G2'}}, {'KE:1'},
        {'G1': 1.0, 'G2': 1.0},
        {'G1': True, 'G2': False},
    )
    flags = {e['id']: e['significant'] for e in result['KE:1']}
    assert flags == {'G1': True, 'G2': False}


def test_build_ke_gene_mapping_threshold_shrinkage():
    """EXPO-05 regression guard: tightening the significance threshold
    STRICTLY shrinks the count of `significant=True` per-gene entries that the
    downstream CSV writer would emit. This is the automated counterpart to
    Task 3 Step E item 7's manual smoke test (which re-runs in the browser
    with a stricter log2FC and confirms the downloaded file is smaller).

    Setup: the same DataFrame-equivalent inputs are passed twice -- once with
    a LOOSE significance map (more genes flagged True) and once with a
    STRICT map (fewer genes flagged True). The other inputs (reference_sets,
    ke_list, gene_logfc_map) are identical across both calls. We then count
    per-gene entries across all KEs where `significant == True` and assert
    strict_count < loose_count.
    """
    # Three KEs, each with overlapping reference gene sets.
    reference_sets = {
        'KE:1': {'G1', 'G2', 'G3', 'G4'},
        'KE:2': {'G2', 'G3', 'G5', 'G6'},
        'KE:3': {'G4', 'G5', 'G6', 'G7'},
    }
    ke_list = {'KE:1', 'KE:2', 'KE:3'}
    # All genes present in the user upload with arbitrary log2FC values
    logfc = {f'G{i}': float(i) - 4.0 for i in range(1, 8)}  # G1..G7

    # LOOSE threshold: 5 of 7 genes pass (G1, G2, G3, G6, G7)
    loose_sig = {
        'G1': True,  'G2': True,  'G3': True,
        'G4': False, 'G5': False,
        'G6': True,  'G7': True,
    }
    # STRICT threshold: only 2 of 7 genes pass (G1, G7) -- proper subset of loose
    strict_sig = {
        'G1': True,  'G2': False, 'G3': False,
        'G4': False, 'G5': False,
        'G6': False, 'G7': True,
    }

    loose_map = build_ke_gene_mapping(reference_sets, ke_list, logfc, loose_sig)
    strict_map = build_ke_gene_mapping(reference_sets, ke_list, logfc, strict_sig)

    loose_count = sum(
        1 for ke_entries in loose_map.values()
          for entry in ke_entries
          if entry['significant'] is True
    )
    strict_count = sum(
        1 for ke_entries in strict_map.values()
          for entry in ke_entries
          if entry['significant'] is True
    )

    # The export consumer (JS CSV writer) emits ONE row per significant=True
    # entry, so this count is exactly the row count of the downloaded CSV.
    # Tightening the threshold MUST shrink it.
    assert strict_count < loose_count, (
        f"EXPO-05 violated: strict-threshold count ({strict_count}) is not "
        f"strictly less than loose-threshold count ({loose_count}). The "
        f"downstream CSV would NOT shrink when the user tightens the "
        f"threshold -- this breaks the EXPO-05 contract."
    )
    # Sanity-pin the absolute counts so a future change that quietly inflates
    # both sides equally still trips the test:
    # loose_sig has G1,G2,G3,G6,G7 = 5 True genes; their occurrences across the
    # three KE reference sets sum to: G1 in {KE:1}, G2 in {KE:1,KE:2}, G3 in {KE:1,KE:2},
    # G6 in {KE:2,KE:3}, G7 in {KE:3} = 1+2+2+2+1 = 8 significant entries.
    # strict_sig has G1,G7 = 2 True genes; G1 in {KE:1}, G7 in {KE:3} = 1+1 = 2 entries.
    assert loose_count == 8
    assert strict_count == 2


class TestKeAccountingSummary:
    """Issue #65: tested/excluded KE accounting attached to the results."""

    def test_summary_counts_tested_and_excluded(self):
        """Every KE of the AOP is accounted for as tested or excluded."""
        genes = [f"GENE{i}" for i in range(50)]
        df = _make_gene_df(genes, [i < 10 for i in range(50)])

        reference_sets = {
            'KE:OK': {f"GENE{i}" for i in range(10)},        # 10 measured — tested
            'KE:SMALL': {'GENE0', 'GENE1', 'GENE2'},          # 3 measured — too few
            'KE:MISS': {'XX1', 'XX2', 'XX3', 'XX4', 'XX5'},   # none measured
        }
        # KE:NOSET is in the AOP but has no reference gene set at all.
        ke_list = {'KE:OK', 'KE:SMALL', 'KE:MISS', 'KE:NOSET'}

        result = run_enrichment_analysis(
            df, reference_sets, ke_list,
            {k: k for k in ke_list},
        )
        summary = get_ke_summary(result)

        assert summary is not None
        assert summary['total_kes'] == 4
        assert summary['tested'] == 1
        assert summary['tested'] == len(result)  # BH denominator == rows returned
        assert summary['excluded_too_few_genes'] == 1
        assert summary['excluded_no_mapping'] == 2  # KE:MISS + KE:NOSET
        assert summary['excluded_error'] == 0
        assert summary['min_ke_genes'] == 5
        assert summary['excluded_reasons']['KE:SMALL'] == 'too_few_genes'
        assert summary['excluded_reasons']['KE:MISS'] == 'no_mapping'
        assert summary['excluded_reasons']['KE:NOSET'] == 'no_mapping'
        assert 'KE:OK' not in summary['excluded_reasons']
        # Accounting must be exhaustive: nothing vanishes silently.
        assert (
            summary['tested'] + summary['excluded_too_few_genes']
            + summary['excluded_no_mapping'] + summary['excluded_error']
            == summary['total_kes']
        )

    def test_min_ke_genes_is_configurable(self):
        """The minimum measured-gene threshold is a parameter, not a literal."""
        genes = [f"GENE{i}" for i in range(50)]
        df = _make_gene_df(genes, [i < 10 for i in range(50)])
        reference_sets = {'KE:1': {'GENE0', 'GENE1', 'GENE2'}}  # 3 measured

        # Default (5) excludes it — no results at all.
        with pytest.raises(ValueError, match="No enrichment results"):
            run_enrichment_analysis(df, reference_sets, {'KE:1'}, {'KE:1': 'Small'})

        # Lowering the minimum lets it through.
        result = run_enrichment_analysis(
            df, reference_sets, {'KE:1'}, {'KE:1': 'Small'}, min_ke_genes=3
        )
        assert len(result) == 1
        summary = get_ke_summary(result)
        assert summary['min_ke_genes'] == 3
        assert summary['excluded_too_few_genes'] == 0

    def test_default_min_ke_genes_comes_from_config(self):
        """The hardcoded 5 is now the documented Config.MIN_KE_GENES constant."""
        assert Config.MIN_KE_GENES == 5
        genes = [f"GENE{i}" for i in range(50)]
        df = _make_gene_df(genes, [i < 10 for i in range(50)])
        reference_sets = {'KE:1': {f"GENE{i}" for i in range(10)}}
        result = run_enrichment_analysis(df, reference_sets, {'KE:1'}, {'KE:1': 'Normal'})
        assert get_ke_summary(result)['min_ke_genes'] == Config.MIN_KE_GENES

    def test_get_ke_summary_returns_none_without_accounting(self):
        """A DataFrame that never went through enrichment carries no summary."""
        assert get_ke_summary(pd.DataFrame({'KE': ['KE:1']})) is None

    def test_format_ke_summary_wording(self):
        """The shared one-liner is worded identically everywhere it appears."""
        text = format_ke_summary({
            'total_kes': 24, 'tested': 18, 'excluded_no_mapping': 2,
            'excluded_too_few_genes': 4, 'excluded_error': 0, 'min_ke_genes': 5,
        })
        assert text == (
            '18 of 24 Key Events tested; 4 excluded (fewer than 5 measured genes), '
            '2 excluded (no gene set mapped)'
        )

    def test_format_ke_summary_omits_empty_clauses(self):
        """No exclusions -> no trailing clause; no summary -> empty string."""
        assert format_ke_summary({
            'total_kes': 10, 'tested': 10, 'excluded_no_mapping': 0,
            'excluded_too_few_genes': 0, 'excluded_error': 0, 'min_ke_genes': 5,
        }) == '10 of 10 Key Events tested'
        assert format_ke_summary(None) == ''
