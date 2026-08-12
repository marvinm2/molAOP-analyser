"""Tests for the GSEA enrichment analysis service.

Covers D-15 determinism (D-15.1), ranking metric (D-15.2), and
zero-p clamp (D-15.3) as well as column shape, min_size filter, and
NES clamping.
"""
import importlib
import math
import logging
import pytest
import pandas as pd
import numpy as np
from services.enrichment_service import (
    EXCLUDED_TOO_FEW_GENES,
    EXCLUDED_TOO_MANY_GENES,
    KE_SUMMARY_ATTR,
    format_ke_summary,
)
from services.gsea_service import (
    DEGENERATE_EXCLUDE,
    DEGENERATE_INFINITY,
    FDR_GSEAPY,
    FDR_RECOMPUTED,
    MIN_SAME_SIGNED_NULL,
    NES_BEYOND_RESOLUTION,
    NES_OK,
    NES_UNDIAGNOSED,
    NES_UNSTABLE,
    _build_ranking,
    _capture_prerank_summaries,
    _null_summaries,
    _null_tail_sizes,
    apply_null_diagnostics,
    recompute_pooled_fdr,
    run_gsea_analysis,
)


def _make_ranked_df(genes, log2fcs, pvals):
    """Helper: build a minimal gene expression DataFrame for GSEA fixtures.

    Lets each test name specific (log2FC, pval) tuples that surface a known
    ranking-metric ordering (positive vs negative vs tied).
    """
    return pd.DataFrame({'ID': genes, 'log2FC': log2fcs, 'pval': pvals})


def _make_gsea_fixture(n_genes=50, n_sig=20, seed=0):
    """Build a synthetic fixture large enough for gseapy.prerank (min_size=5 overlap)."""
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i:03d}" for i in range(n_genes)]
    log2fcs = rng.normal(0, 1.5, n_genes).tolist()
    # Avoid zero pvalues in base fixture
    pvals = rng.uniform(1e-6, 1.0, n_genes).tolist()
    return pd.DataFrame({'ID': genes, 'log2FC': log2fcs, 'pval': pvals})


def _make_reference_sets(genes, ke_ids, overlap_size=10):
    """Build reference sets where each KE overlaps with the first `overlap_size` genes."""
    gene_list = list(genes)
    return {ke: set(gene_list[:overlap_size]) for ke in ke_ids}


class TestBuildRanking:
    """Unit tests for _build_ranking helper (D-15.2, D-15.3)."""

    def test_ranking_metric_signed_log_p(self):
        """D-15.2: ranking = sign(log2FC) * -log10(pval); ordering must match."""
        # Gene A: positive log2FC, very small pval → large positive metric
        # Gene B: negative log2FC, very small pval → large negative metric
        # Gene C: positive log2FC, large pval → small positive metric
        # Gene D: negative log2FC, large pval → small negative metric
        genes = ["GENEA", "GENEB", "GENEC", "GENED"]
        log2fcs = [2.0, -2.0, 0.5, -0.5]
        pvals = [1e-10, 1e-10, 0.5, 0.5]
        df = _make_ranked_df(genes, log2fcs, pvals)

        ranking = _build_ranking(df)

        assert isinstance(ranking, pd.Series)
        # Gene A should have the largest positive metric
        assert ranking["GENEA"] > 0
        assert ranking["GENEB"] < 0
        # Gene A (big positive) > Gene C (small positive)
        assert ranking["GENEA"] > ranking["GENEC"]
        # Gene B (big negative) < Gene D (small negative)
        assert ranking["GENEB"] < ranking["GENED"]
        # Top rank is Gene A
        assert ranking.index[0] == "GENEA"
        # Bottom rank is Gene B
        assert ranking.index[-1] == "GENEB"

    def test_ranking_zero_pvalue_clamp(self):
        """D-15.3: zero pvalue is clamped to smallest non-zero pval; no +inf in result."""
        genes = ["GENEUP", "GENEDOWN", "GENENORM"]
        log2fcs = [3.0, -1.0, 0.5]
        # GENEUP has pval=0.0 (should be clamped)
        pvals = [0.0, 1e-5, 0.1]
        df = _make_ranked_df(genes, log2fcs, pvals)

        ranking = _build_ranking(df)

        # No +inf values
        assert not any(math.isinf(v) for v in ranking.values)
        assert not any(math.isnan(v) for v in ranking.values)

        # The zero-pval gene's metric should equal sign(3.0) * -log10(1e-5)
        # because 1e-5 is the smallest nonzero pval
        expected_metric = 1.0 * -math.log10(1e-5)
        assert abs(ranking["GENEUP"] - expected_metric) < 1e-10

    def test_ties_broken_by_gene_symbol_ascending(self):
        """Genes with identical metrics resolve to ascending gene-symbol order."""
        # Three genes with identical (log2FC, pval) → identical metric values
        genes = ["GENEC", "GENEA", "GENEB"]
        log2fcs = [1.0, 1.0, 1.0]
        pvals = [0.01, 0.01, 0.01]
        df = _make_ranked_df(genes, log2fcs, pvals)

        ranking = _build_ranking(df)

        # After tie-breaking by ascending gene symbol (A < B < C)
        # but sort_values(ascending=False) means the first element is the highest metric
        # All metrics are equal, so order should be A, B, C (mergesort stable after sort_index)
        tied_genes = list(ranking.index)
        assert tied_genes == ["GENEA", "GENEB", "GENEC"]


class TestRunGseaAnalysis:
    """Integration tests for run_gsea_analysis (D-15.1, column shape, min_size, NES clamp)."""

    def _build_fixture(self):
        """50-gene fixture with two KEs of 10-gene overlap each."""
        df = _make_gsea_fixture(n_genes=50, seed=42)
        ke_ids = ["KE:001", "KE:002"]
        reference_sets = _make_reference_sets(df["ID"].tolist(), ke_ids, overlap_size=10)
        ke_list = set(ke_ids)
        ke_title_map = {"KE:001": "Key Event One", "KE:002": "Key Event Two"}
        return df, reference_sets, ke_list, ke_title_map

    def test_gsea_determinism(self):
        """D-15.1: two identical calls produce byte-identical DataFrames."""
        df, reference_sets, ke_list, ke_title_map = self._build_fixture()

        result1 = run_gsea_analysis(df, reference_sets, ke_list, ke_title_map)
        result2 = run_gsea_analysis(df, reference_sets, ke_list, ke_title_map)

        assert result1.reset_index(drop=True).equals(result2.reset_index(drop=True)), (
            "Two consecutive GSEA runs on identical input produced different results. "
            "threads=1 and seed=42 must be enforced."
        )

    def test_run_gsea_columns(self):
        """Returned DataFrame has exactly the expected columns in the right order."""
        df, reference_sets, ke_list, ke_title_map = self._build_fixture()

        result = run_gsea_analysis(df, reference_sets, ke_list, ke_title_map)

        expected_cols = ['KE', 'Title', 'NES', 'nes_status', 'ES', 'p_value',
                         'p_value_resolution', 'FDR', 'fdr_source',
                         'null_same_signed_n', 'lead_genes',
                         'total_KE_genes_in_dataset']
        assert list(result.columns) == expected_cols

    def test_min_size_filter(self, caplog):
        """KE with fewer than min_size=5 overlapping genes is excluded; a log message is emitted."""
        df = _make_gsea_fixture(n_genes=50, seed=42)
        gene_list = df["ID"].tolist()

        # KE:SMALL has only 3 overlapping genes — below min_size=5
        # KE:LARGE has 10 overlapping genes — kept
        reference_sets = {
            "KE:SMALL": set(gene_list[:3]),
            "KE:LARGE": set(gene_list[:10]),
        }
        ke_list = {"KE:SMALL", "KE:LARGE"}
        ke_title_map = {"KE:SMALL": "Too Small KE", "KE:LARGE": "Large KE"}

        with caplog.at_level(logging.INFO, logger="services.gsea_service"):
            result = run_gsea_analysis(df, reference_sets, ke_list, ke_title_map)

        # KE:SMALL should not appear in results
        assert "KE:SMALL" not in result["KE"].values
        # KE:LARGE should appear (enough overlap)
        assert "KE:LARGE" in result["KE"].values
        # A log message naming the dropped KE should have been emitted
        assert any("KE:SMALL" in record.message for record in caplog.records), (
            "Expected a log message naming the dropped KE:SMALL"
        )

    def test_nes_clamp(self):
        """NES values in returned DataFrame all satisfy -3 <= NES <= 3."""
        df, reference_sets, ke_list, ke_title_map = self._build_fixture()

        result = run_gsea_analysis(df, reference_sets, ke_list, ke_title_map)

        for nes in result["NES"]:
            assert -3.0 <= nes <= 3.0, f"NES value {nes} is outside ±3 clamp range"


def _one_signed_null_fixture(n_genes=200, set_size=100, heavy_tail=10):
    """A ranking whose permutation null lands entirely below zero (issue #117).

    The condition needs two things together: a gene set large relative to the
    ranking (half of it here), and a handful of genes carrying most of the
    ranking weight at the bottom. Any random set of that size then walks
    steeply down before it can recover, so every permuted enrichment score is
    negative — while the real set, sitting entirely at the top, scores +1.0.
    That is the maximally enriched case, and it is the one gseapy reports as
    NES = 1.0 with p = 1.0.

    Verified against origin/main (db31a90), which returns exactly
    ``NES 1.0, p_value 1.0`` for KE:A on this fixture.

    Returns:
        tuple: (df, reference_sets, ke_list, ke_title_map)
    """
    genes = [f"G{i:04d}" for i in range(n_genes)]
    log2fcs = [1.0] * set_size + [-1.0] * (n_genes - set_size)
    # The heavy tail is what tightens the null onto one side of zero.
    pvals = ([1e-4] * set_size
             + [0.5] * (n_genes - set_size - heavy_tail)
             + [1e-60] * heavy_tail)
    df = pd.DataFrame({'ID': genes, 'log2FC': log2fcs, 'pval': pvals})
    reference_sets = {'KE:A': set(genes[:set_size])}
    return df, reference_sets, {'KE:A'}, {'KE:A': 'Maximally enriched KE'}


class TestOneSignedNull:
    """Issue #117 — the degenerate and near-degenerate NES regimes."""

    def test_empty_same_signed_tail_is_not_reported_as_null(self, caplog):
        """The maximally enriched KE must not come back as NES 1.0 / p 1.0."""
        df, reference_sets, ke_list, ke_title_map = _one_signed_null_fixture()

        with caplog.at_level(logging.WARNING, logger="services.gsea_service"):
            result = run_gsea_analysis(
                df, reference_sets, ke_list, ke_title_map, permutation_num=100
            )

        row = result[result['KE'] == 'KE:A'].iloc[0]
        assert row['null_same_signed_n'] == 0, (
            "fixture no longer produces a one-signed null; the test is not "
            "exercising issue #117 any more"
        )
        assert row['nes_status'] == NES_BEYOND_RESOLUTION
        # The pair that reads as "tested, not enriched" must be gone.
        assert not (row['NES'] == 1.0 and row['p_value'] == 1.0)
        # NES is not normalisable and must not be a number that can be ranked.
        assert math.isnan(row['NES'])
        # p is the resolution bound, not 1.0; FDR is the limit, not mid-range.
        assert row['p_value'] == pytest.approx(1.0 / 100)
        assert row['FDR'] == 0.0
        # Direction survives, so the reader still knows which way it went.
        assert row['ES'] > 0
        assert any('issue #117' in r.getMessage() for r in caplog.records)

    def test_unstable_tail_keeps_the_pvalue_but_flags_the_magnitude(self):
        """A tail of 1-9 permutations is a real call with an unusable NES."""
        res = pd.DataFrame({
            'KE': ['KE:1', 'KE:2', 'KE:3'],
            'NES': [1.12, 2.4, 1.0],
            'p_value': [0.0, 0.001, 1.0],
            'FDR': [0.23, 0.001, 0.44],
        })

        out = apply_null_diagnostics(res, {'KE:1': 3, 'KE:2': 400, 'KE:3': 0}, 1000)

        unstable = out[out['KE'] == 'KE:1'].iloc[0]
        assert unstable['nes_status'] == NES_UNSTABLE
        # Nothing is rewritten: gseapy's p-value is sound in this regime.
        assert unstable['NES'] == 1.12
        assert unstable['p_value'] == 0.0
        assert unstable['FDR'] == 0.23
        assert out[out['KE'] == 'KE:2'].iloc[0]['nes_status'] == NES_OK
        assert out[out['KE'] == 'KE:3'].iloc[0]['nes_status'] == NES_BEYOND_RESOLUTION

    def test_value_based_detection_would_miss_the_unstable_regime(self):
        """A NES of exactly 1.0 is not the condition — the tail size is.

        Guards the choice of test: KE:1 carries a plausible NES of 1.12 and a
        p-value of 0.0, so no check on the reported values can find it, yet its
        magnitude rests on three permutations.
        """
        res = pd.DataFrame({
            'KE': ['KE:1'], 'NES': [1.12], 'p_value': [0.0], 'FDR': [0.23],
        })

        out = apply_null_diagnostics(res, {'KE:1': 3}, 1000)

        assert out.iloc[0]['nes_status'] == NES_UNSTABLE

    def test_missing_diagnostics_leave_the_numbers_alone(self):
        """No captured null means "not assessed", never "tail of zero"."""
        res = pd.DataFrame({
            'KE': ['KE:1'], 'NES': [1.0], 'p_value': [1.0], 'FDR': [0.44],
        })

        out = apply_null_diagnostics(res, {}, 1000)

        # gseapy's numbers are untouched...
        assert out.iloc[0]['NES'] == 1.0
        assert out.iloc[0]['p_value'] == 1.0
        assert math.isnan(out.iloc[0]['null_same_signed_n'])
        # ...but the row does NOT claim to have been checked. Labelling this
        # 'ok' is how a gseapy rename would reintroduce the whole of #117 under
        # a clean bill of health.
        assert out.iloc[0]['nes_status'] == NES_UNDIAGNOSED
        assert out.iloc[0]['nes_status'] != NES_OK

    def test_null_tail_sizes_counts_the_observed_side(self):
        """Same-signed means the side the observed ES fell on, not the majority."""
        class _Summary:
            def __init__(self, term, es, esnull):
                self.term, self.es, self.esnull = term, es, esnull

        sizes = _null_tail_sizes([
            _Summary('KE:POS', 0.33, [-0.5, -0.4, -0.3]),
            _Summary('KE:NEG', -0.33, [-0.5, 0.4, 0.3]),
            _Summary('KE:BROKEN', 'x', [0.1]),
        ])

        assert sizes['KE:POS'] == 0
        assert sizes['KE:NEG'] == 1
        assert 'KE:BROKEN' not in sizes

    def test_prerank_summaries_are_captured_and_the_patch_reverted(self):
        """The capture must not leave gseapy patched, whatever happens."""
        gsea_module = importlib.import_module('gseapy.gsea')
        original = gsea_module.prerank_rs

        sink = []
        with _capture_prerank_summaries(sink):
            assert gsea_module.prerank_rs is not original
        assert gsea_module.prerank_rs is original

        with pytest.raises(RuntimeError):
            with _capture_prerank_summaries(sink):
                raise RuntimeError('boom')
        assert gsea_module.prerank_rs is original

    def test_capture_populates_the_sink_on_a_real_run(self):
        """The wrapper actually sees the summaries gseapy discards."""
        df, reference_sets, ke_list, ke_title_map = _one_signed_null_fixture()

        result = run_gsea_analysis(
            df, reference_sets, ke_list, ke_title_map, permutation_num=100
        )

        assert result['null_same_signed_n'].notna().all(), (
            "null diagnostics were not captured — the gseapy entry point this "
            "wraps has probably been renamed"
        )


class _Summary:
    """Stand-in for gseapy's ``GSEASummary`` — term, observed ES, null ES."""

    def __init__(self, term, es, esnull):
        self.term, self.es, self.esnull = term, es, esnull


def _healthy_summaries(n_perm=40, seed=0):
    """Four terms whose permutation nulls straddle zero, so none is degenerate.

    Deterministic by construction: the null for each term is a fixed symmetric
    ladder rather than a draw, so the pooled distribution — and therefore every
    q-value computed from it — is identical on every run and on every platform.
    """
    step = np.linspace(-0.6, 0.6, n_perm)
    return [
        _Summary('KE:STRONG', 0.62, step * 1.0),
        _Summary('KE:MID', 0.24, step * 0.9),
        _Summary('KE:WEAK', 0.08, step * 1.1),
        _Summary('KE:DOWN', -0.31, step * 0.95),
    ]


class TestPooledNullFdr:
    """Issue #122 — the fabricated NES contaminating everybody else's q."""

    def test_matches_gseapy_exactly_when_no_term_is_degenerate(self):
        """The transcription must be a transcription, not a second opinion.

        This is the check that makes the whole change reviewable: with no
        degenerate term there is nothing to remove from the pools, so the
        recomputation has to land on gseapy's own numbers. Any drift here means
        the rewrite changed the statistic rather than repairing its input, and
        every q-value the tool reports would be quietly its own invention.
        """
        from gseapy.algorithm import gsea_fdr, normalize

        summaries = _healthy_summaries()
        es = np.array([s.es for s in summaries], dtype=float)
        esnull = np.vstack([np.asarray(s.esnull, dtype=float) for s in summaries])

        nes, nesnull = normalize(es, esnull)
        expected = gsea_fdr(nes, nesnull)

        got = recompute_pooled_fdr(summaries)

        assert set(got) == {s.term for s in summaries}
        for term, want in zip([s.term for s in summaries], expected):
            assert got[term] == pytest.approx(min(want, 1.0), abs=1e-12), (
                f"{term}: recomputed q diverged from gseapy's own on a run with "
                f"no degenerate term"
            )

    def test_substituted_nes_biases_other_terms_in_both_directions(self):
        """The defect, reproduced against gseapy's real gsea_fdr.

        A fifth term whose null lies entirely below zero while its observed ES
        is positive is the degenerate case: gseapy cannot normalise it and the
        Rust path substitutes NES = 1.0. That fabricated value lands in the
        *observed* vector, which is the denominator every other term's q is
        divided by — so it moves terms it has nothing to do with, and the
        direction depends on which side of 1.0 they sit.
        """
        from gseapy.algorithm import gsea_fdr, normalize

        summaries = _healthy_summaries() + [
            _Summary('KE:DEGENERATE', 0.71, np.linspace(-0.9, -0.2, 40)),
        ]
        es = np.array([s.es for s in summaries], dtype=float)
        esnull = np.vstack([np.asarray(s.esnull, dtype=float) for s in summaries])

        with np.errstate(invalid='ignore', divide='ignore'):
            nes, nesnull = normalize(es, esnull)
        # Model the Rust path exactly: the unnormalisable observed score becomes
        # 1.0, and that is what gseapy pools.
        nes = np.where(np.isfinite(nes), nes, 1.0)
        nes[-1] = 1.0
        contaminated = dict(zip([s.term for s in summaries], gsea_fdr(nes, nesnull)))

        repaired = recompute_pooled_fdr(summaries)

        # Something actually moved — otherwise the fixture is not exercising
        # the defect and the direction assertions below are vacuous.
        moved = [t for t in repaired if repaired[t] != pytest.approx(
            min(contaminated[t], 1.0), abs=1e-12)]
        assert moved, "the fabricated NES changed nobody's q — fixture is inert"

        for term in ('KE:STRONG', 'KE:MID', 'KE:WEAK'):
            observed_nes = nes[[s.term for s in summaries].index(term)]
            gseapy_q = min(contaminated[term], 1.0)
            if observed_nes > 1.0:
                assert gseapy_q >= repaired[term], (
                    f"{term} sits above NES 1.0, where the fabricated value "
                    f"inflates q; the repair must not raise it further"
                )
            else:
                assert gseapy_q <= repaired[term], (
                    f"{term} sits at or below NES 1.0, where the fabricated "
                    f"value deflates q; the repair must not lower it further"
                )

    def test_degenerate_term_lands_on_the_zero_limit(self):
        """The infinity treatment must reproduce #117's hand-written limit.

        ``apply_null_diagnostics`` writes FDR = 0.0 for a term whose observed ES
        beat every permutation. If the recomputation disagreed, the two repairs
        would be fighting each other and the row would carry whichever ran last.
        """
        summaries = _healthy_summaries() + [
            _Summary('KE:DEGENERATE', 0.71, np.linspace(-0.9, -0.2, 40)),
        ]

        got = recompute_pooled_fdr(summaries, degenerate=DEGENERATE_INFINITY)

        assert got['KE:DEGENERATE'] == 0.0

    def test_excluding_degenerate_terms_is_available_and_differs(self):
        """#122's own proposal, kept as a comparison rather than as the default.

        Excluding drops a term that genuinely belongs in the observed
        distribution, so it deflates the denominator for everyone else — the
        same class of error, smaller and opposite-signed. Shipping both is what
        lets the difference be measured on real data instead of argued about.
        """
        summaries = _healthy_summaries() + [
            _Summary('KE:DEGENERATE', 0.71, np.linspace(-0.9, -0.2, 40)),
        ]

        infinity = recompute_pooled_fdr(summaries, degenerate=DEGENERATE_INFINITY)
        excluded = recompute_pooled_fdr(summaries, degenerate=DEGENERATE_EXCLUDE)

        assert 'KE:DEGENERATE' in infinity
        assert 'KE:DEGENERATE' not in excluded
        assert any(
            infinity[t] != pytest.approx(excluded[t], abs=1e-12)
            for t in excluded
        ), "the two treatments produced identical pools — one of them is wrong"

    def test_ragged_nulls_decline_rather_than_pool_mismatched_material(self, caplog):
        """Unequal null widths mean two runs were captured together (#117's lock)."""
        summaries = [
            _Summary('KE:A', 0.5, np.linspace(-0.5, 0.5, 40)),
            _Summary('KE:B', 0.4, np.linspace(-0.5, 0.5, 25)),
        ]

        with caplog.at_level(logging.ERROR, logger="services.gsea_service"):
            assert recompute_pooled_fdr(summaries) == {}
        assert any('not rectangular' in r.getMessage() for r in caplog.records)

    def test_unreadable_summaries_leave_gseapy_in_charge(self):
        """No nulls means no recomputation — never a fabricated q of its own."""
        assert recompute_pooled_fdr([]) == {}
        assert recompute_pooled_fdr([_Summary('KE:X', 'not a number', [0.1])]) == {}

    def test_rejects_an_unknown_degenerate_treatment(self):
        """A typo in the mode must not silently select a default."""
        with pytest.raises(ValueError, match='degenerate must be'):
            recompute_pooled_fdr(_healthy_summaries(), degenerate='drop')

    def test_every_row_says_whose_q_value_it_carries(self):
        """`fdr_source` is the provenance #122's third criterion asks for."""
        res = pd.DataFrame({
            'KE': ['KE:1', 'KE:2'],
            'NES': [1.5, 0.9],
            'p_value': [0.01, 0.4],
            'FDR': [0.30, 0.44],
        })

        out = apply_null_diagnostics(
            res, {'KE:1': 400, 'KE:2': 400}, 1000,
            recomputed_fdr={'KE:1': 0.11},
        )

        assert out.iloc[0]['FDR'] == 0.11
        assert out.iloc[0]['fdr_source'] == FDR_RECOMPUTED
        # Untouched rows keep gseapy's number and say so, rather than being
        # presented as recomputed alongside rows that actually were.
        assert out.iloc[1]['FDR'] == 0.44
        assert out.iloc[1]['fdr_source'] == FDR_GSEAPY

    def test_no_recomputation_labels_the_whole_run_as_gseapy(self):
        """The pre-#122 path still works and does not claim to be the new one."""
        res = pd.DataFrame({
            'KE': ['KE:1'], 'NES': [1.0], 'p_value': [1.0], 'FDR': [0.44],
        })

        out = apply_null_diagnostics(res, {'KE:1': 400}, 1000)

        assert out.iloc[0]['FDR'] == 0.44
        assert out.iloc[0]['fdr_source'] == FDR_GSEAPY

    def test_null_summaries_keep_the_arrays_the_tail_count_throws_away(self):
        """#117 reduced each null to an integer; #122 needs the distribution."""
        got = _null_summaries([
            _Summary('KE:OK', 0.3, [-0.5, 0.4, np.inf]),
            _Summary('KE:BROKEN', 'x', [0.1]),
        ])

        assert 'KE:BROKEN' not in got
        es, null = got['KE:OK']
        assert es == 0.3
        # The non-finite entry is dropped, not carried into a pooled sort where
        # it would land at the top of the distribution.
        assert null.tolist() == [-0.5, 0.4]


class TestMaxSizeExclusion:
    """Issue #120 — a KE above the GSEA ceiling is dropped, and must say so."""

    def _fixture(self):
        genes = [f"G{i:04d}" for i in range(200)]
        df = pd.DataFrame({
            'ID': genes,
            'log2FC': [1.0] * 100 + [-1.0] * 100,
            'pval': [1e-4] * 100 + [0.5] * 100,
        })
        reference_sets = {'KE:BIG': set(genes[:60]), 'KE:OK': set(genes[:20])}
        return df, reference_sets, {'KE:BIG', 'KE:OK'}, {}

    def test_over_ceiling_ke_is_not_reported_as_under_measured(self, caplog):
        """60 measured genes must not be reported as "fewer than 5"."""
        df, reference_sets, ke_list, ke_title_map = self._fixture()

        with caplog.at_level(logging.INFO, logger="services.gsea_service"):
            result = run_gsea_analysis(
                df, reference_sets, ke_list, ke_title_map,
                permutation_num=50, max_size=30,
            )

        summary = result.attrs[KE_SUMMARY_ATTR]
        assert summary['excluded_reasons']['KE:BIG'] == EXCLUDED_TOO_MANY_GENES
        assert summary['excluded_too_many_genes'] == 1
        assert summary['excluded_too_few_genes'] == 0
        assert summary['max_ke_genes'] == 30
        # The log must name the real cause and the real count.
        messages = [r.getMessage() for r in caplog.records]
        assert any('max_size=30' in m and '60 measured genes' in m
                   for m in messages), messages
        assert not any('min_size' in m and 'KE:BIG' in m for m in messages)

    def test_accounting_sentence_names_the_ceiling(self):
        """The user-facing sentence must not say the KE was under-measured."""
        df, reference_sets, ke_list, ke_title_map = self._fixture()

        result = run_gsea_analysis(
            df, reference_sets, ke_list, ke_title_map,
            permutation_num=50, max_size=30,
        )

        sentence = format_ke_summary(result.attrs[KE_SUMMARY_ATTR])
        assert 'more than 30 measured genes' in sentence
        assert 'fewer than' not in sentence

    def test_min_size_drop_still_reads_as_too_few(self):
        """The existing min_size accounting is unchanged."""
        df, reference_sets, ke_list, ke_title_map = self._fixture()
        reference_sets = dict(reference_sets)
        reference_sets['KE:TINY'] = set(df['ID'][:3])
        ke_list = set(ke_list) | {'KE:TINY'}

        result = run_gsea_analysis(
            df, reference_sets, ke_list, ke_title_map, permutation_num=50,
        )

        summary = result.attrs[KE_SUMMARY_ATTR]
        assert summary['excluded_reasons']['KE:TINY'] == EXCLUDED_TOO_FEW_GENES
        assert summary['excluded_too_many_genes'] == 0
