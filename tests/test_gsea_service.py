"""Tests for the GSEA enrichment analysis service.

Covers D-15 determinism (D-15.1), ranking metric (D-15.2), and
zero-p clamp (D-15.3) as well as column shape, min_size filter, and
NES clamping.
"""
import math
import logging
import pytest
import pandas as pd
import numpy as np
from services.gsea_service import run_gsea_analysis, _build_ranking


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

        expected_cols = ['KE', 'Title', 'NES', 'p_value', 'FDR', 'lead_genes',
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
