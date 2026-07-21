"""Tests for the data processing service."""
import os
import tempfile
import pytest
import pandas as pd
from services.data_service import (
    load_and_validate_data,
    process_gene_expression,
    compute_logfc_percentiles,
)


class TestLoadAndValidateData:
    """Tests for file loading and validation."""

    def _write_tsv(self, path, data):
        df = pd.DataFrame(data)
        df.to_csv(path, sep='\t', index=False)
        return path

    def test_loads_valid_tsv(self, tmp_path):
        path = self._write_tsv(str(tmp_path / 'test.tsv'), {
            'Gene': ['BRCA1', 'TP53', 'EGFR'],
            'logFC': [1.5, -0.8, 2.1],
            'pval': [0.001, 0.05, 0.0001],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        assert len(df) == 3
        assert 'ID' in df.columns
        assert 'log2FC' in df.columns
        assert 'pval' in df.columns

    def test_missing_column_raises(self, tmp_path):
        path = self._write_tsv(str(tmp_path / 'test.tsv'), {
            'Gene': ['BRCA1'], 'logFC': [1.5],
        })
        with pytest.raises(Exception):
            load_and_validate_data(path, 'Gene', 'logFC', 'MISSING_COL')

    def test_expands_triple_slash_genes(self, tmp_path):
        """Genes separated by /// should be expanded into separate rows."""
        path = self._write_tsv(str(tmp_path / 'test.tsv'), {
            'Gene': ['BRCA1///TP53', 'EGFR'],
            'logFC': [1.5, 2.0],
            'pval': [0.001, 0.01],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        assert len(df) >= 3  # BRCA1, TP53, EGFR
        assert 'BRCA1' in df['ID'].values
        assert 'TP53' in df['ID'].values

    def test_gene_ids_uppercased(self, tmp_path):
        path = self._write_tsv(str(tmp_path / 'test.tsv'), {
            'Gene': ['brca1', 'tp53'],
            'logFC': [1.5, -0.8],
            'pval': [0.001, 0.05],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        assert all(g == g.upper() for g in df['ID'])


class TestProcessGeneExpression:
    """Tests for gene expression processing and significance assignment."""

    def _make_df(self, genes, fc_vals, pvals):
        return pd.DataFrame({'ID': genes, 'log2FC': fc_vals, 'pval': pvals})

    def test_basic_processing(self):
        df = self._make_df(['BRCA1', 'TP53', 'EGFR'],
                           [2.0, 0.1, -1.5], [0.001, 0.5, 0.01])
        result, stats = process_gene_expression(df, logfc_threshold=1.0)
        assert stats['total_genes'] == 3
        assert 'significant' in result.columns

    def test_threshold_zero_all_significant_if_low_pval(self):
        """With threshold=0, all genes with pval<=cutoff are significant."""
        df = self._make_df(['A', 'B', 'C'], [0.5, 0.1, 2.0], [0.01, 0.01, 0.01])
        result, stats = process_gene_expression(df, logfc_threshold=0.0)
        assert stats['significant_genes'] == 3

    def test_high_threshold_reduces_significant(self):
        df = self._make_df(['A', 'B', 'C'], [0.5, 1.5, 2.0], [0.01, 0.01, 0.01])
        result, stats = process_gene_expression(df, logfc_threshold=1.0)
        # Only B (1.5) and C (2.0) have |logFC| >= 1.0
        assert stats['significant_genes'] == 2

    def test_pval_threshold_honoured(self):
        """A custom p-value threshold controls which genes are flagged significant."""
        df = self._make_df(['A', 'B', 'C'], [2.0, 2.0, 2.0], [0.03, 0.07, 0.2])
        # Stricter than default 0.05: only A (0.03) qualifies.
        strict, strict_stats = process_gene_expression(df, logfc_threshold=0.0, pval_threshold=0.05)
        assert strict_stats['significant_genes'] == 1
        assert bool(strict[strict['ID'] == 'A']['significant'].iloc[0]) is True
        assert bool(strict[strict['ID'] == 'B']['significant'].iloc[0]) is False
        # Looser threshold lets B (0.07) through too.
        loose, loose_stats = process_gene_expression(df, logfc_threshold=0.0, pval_threshold=0.1)
        assert loose_stats['significant_genes'] == 2

    def test_pval_threshold_none_uses_default(self):
        """When pval_threshold is None, the legacy Config.PVAL_CUTOFF (0.05) applies."""
        df = self._make_df(['A', 'B'], [2.0, 2.0], [0.04, 0.06])
        result, stats = process_gene_expression(df, logfc_threshold=0.0, pval_threshold=None)
        assert stats['significant_genes'] == 1
        assert bool(result[result['ID'] == 'A']['significant'].iloc[0]) is True

    def test_duplicate_genes_combined(self):
        """Duplicate gene IDs should be combined via Fisher's method."""
        df = self._make_df(['BRCA1', 'BRCA1', 'TP53'],
                           [1.5, 2.0, 0.5], [0.01, 0.02, 0.05])
        result, stats = process_gene_expression(df, logfc_threshold=0.0)
        assert stats['total_genes'] == 2  # BRCA1 combined, TP53 separate
        brca1_row = result[result['ID'] == 'BRCA1']
        assert len(brca1_row) == 1
        # Combined logFC should be mean of 1.5 and 2.0
        assert abs(brca1_row.iloc[0]['log2FC'] - 1.75) < 0.01

    def test_empty_after_processing_warns(self):
        """All genes above p-value cutoff should yield zero significant genes."""
        df = self._make_df(['A', 'B'], [2.0, 2.0], [0.9, 0.9])
        result, stats = process_gene_expression(df, logfc_threshold=0.0)
        assert stats['significant_genes'] == 0


class TestComputeLogfcPercentiles:
    """Top N% quick-threshold calculation (issue #64)."""

    def test_known_percentiles(self):
        """Threshold is the |log2FC| at the Nth percentile, descending."""
        # 100 genes with |log2FC| of 1..100. Sorted descending, index 10 is the
        # 11th largest value (90.0) and index 20 the 21st (80.0).
        values = pd.Series(range(1, 101))
        result = compute_logfc_percentiles(values)
        assert result == {10: 90.0, 20: 80.0}

    def test_uses_absolute_values(self):
        """Down-regulated genes count toward the threshold too."""
        result = compute_logfc_percentiles(pd.Series([-5.0, -4.0, 3.0, 2.0, 1.0]))
        # |values| descending: 5, 4, 3, 2, 1 -> index 0 for both 10% and 20%.
        assert result[10] == 5.0

    def test_unaffected_by_display_cap(self):
        """
        The regression this fixes: percentiles must reflect the full dataset,
        not the truncated volcano payload.
        """
        from config import Config

        # A dataset larger than the display cap, with the largest fold changes
        # placed AFTER the cap so a truncated calculation cannot see them.
        n = Config.MAX_GENES_DISPLAY + 4000
        values = pd.Series(list(range(1, n + 1)))

        full = compute_logfc_percentiles(values)
        truncated = compute_logfc_percentiles(values.head(Config.MAX_GENES_DISPLAY))

        assert full != truncated
        assert full[10] == round(float(n - int(n * 0.10)), 2)

    def test_small_dataset_does_not_overrun(self):
        """Index must stay in range for tiny inputs."""
        assert compute_logfc_percentiles(pd.Series([1.0])) == {10: 1.0, 20: 1.0}

    def test_empty_input_returns_empty(self):
        """No usable data yields no thresholds rather than raising."""
        assert compute_logfc_percentiles(pd.Series([], dtype=float)) == {}
        assert compute_logfc_percentiles(pd.Series([None, None], dtype=object)) == {}
