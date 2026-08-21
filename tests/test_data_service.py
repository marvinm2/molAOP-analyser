"""Tests for the data processing service."""
import os
import tempfile
import pytest
import pandas as pd
from exceptions import DataValidationError
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


class TestMissingGeneSymbols:
    """Rows without a gene symbol must not become a pseudo-gene (issue #80)."""

    def _write_tsv(self, path, data):
        pd.DataFrame(data).to_csv(path, sep='\t', index=False)
        return path

    def test_missing_symbol_does_not_become_nan_gene(self, tmp_path):
        path = self._write_tsv(str(tmp_path / 'missing.tsv'), {
            'Gene': ['BRCA1', None, 'TP53', None],
            'logFC': [1.5, 3.0, -0.8, 2.2],
            'pval': [0.001, 0.001, 0.05, 0.001],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        assert 'NAN' not in set(df['ID'])
        assert set(df['ID']) == {'BRCA1', 'TP53'}

    def test_background_is_the_number_of_measured_genes(self, tmp_path):
        """The background must not be inflated by the discarded rows."""
        path = self._write_tsv(str(tmp_path / 'bg.tsv'), {
            'Gene': ['BRCA1', None, 'TP53', float('nan'), 'EGFR'],
            'logFC': [1.5, 3.0, -0.8, 2.2, 1.1],
            'pval': [0.001, 0.001, 0.05, 0.001, 0.001],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        _, stats = process_gene_expression(df, logfc_threshold=0.0)
        assert stats['total_genes'] == 3

    def test_dropped_row_count_is_reported(self, tmp_path):
        path = self._write_tsv(str(tmp_path / 'dropped.tsv'), {
            'Gene': ['BRCA1', None, None, 'TP53'],
            'logFC': [1.5, 3.0, 2.2, -0.8],
            'pval': [0.001, 0.001, 0.001, 0.05],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        assert df.attrs['dropped_unidentified_rows'] == 2
        _, stats = process_gene_expression(df, logfc_threshold=0.0)
        assert stats['dropped_unidentified_rows'] == 2

    def test_no_discards_reports_zero(self, tmp_path):
        path = self._write_tsv(str(tmp_path / 'clean.tsv'), {
            'Gene': ['BRCA1', 'TP53'],
            'logFC': [1.5, -0.8],
            'pval': [0.001, 0.05],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        _, stats = process_gene_expression(df, logfc_threshold=0.0)
        assert stats['dropped_unidentified_rows'] == 0

    def test_placeholder_strings_are_discarded(self, tmp_path):
        """Literal NA / nan / empty placeholders are not gene symbols."""
        path = self._write_tsv(str(tmp_path / 'placeholders.tsv'), {
            'Gene': ['BRCA1', 'NA', 'nan', '-', 'TP53'],
            'logFC': [1.5, 3.0, 2.0, 1.0, -0.8],
            'pval': [0.001, 0.001, 0.001, 0.001, 0.05],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        assert set(df['ID']) == {'BRCA1', 'TP53'}
        assert df.attrs['dropped_unidentified_rows'] == 3

    def test_partial_multi_symbol_row_keeps_real_symbols(self, tmp_path):
        """A '///' row with one usable symbol keeps it and is not counted as dropped."""
        path = self._write_tsv(str(tmp_path / 'partial.tsv'), {
            'Gene': ['BRCA1///NA', 'TP53'],
            'logFC': [1.5, -0.8],
            'pval': [0.001, 0.05],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        assert set(df['ID']) == {'BRCA1', 'TP53'}
        assert df.attrs['dropped_unidentified_rows'] == 0

    def test_no_usable_identifier_raises_named_validation_error(self, tmp_path):
        """A wrong ID column must be reported, not crash downstream.

        Dropping symbol-less rows can empty the frame outright (the realistic
        cause is a user pointing the wizard at the wrong column). An empty frame
        used to reach process_gene_expression(), whose groupby yields no rows,
        so the next read of ['ID'] raised a bare KeyError and the request fell
        through to the "unexpected error" catch-all.
        """
        path = self._write_tsv(str(tmp_path / 'wrongcol.tsv'), {
            'Gene': ['BRCA1', 'TP53'],
            'Empty': [None, None],
            'logFC': [1.5, -0.8],
            'pval': [0.001, 0.05],
        })
        with pytest.raises(DataValidationError) as excinfo:
            load_and_validate_data(path, 'Empty', 'logFC', 'pval')
        message = str(excinfo.value)
        assert 'Empty' in message
        assert 'identifier' in message.lower()

    def test_pandas_na_vocabulary_is_covered(self, tmp_path):
        """The placeholder set agrees with pandas' own missing-value words."""
        path = self._write_tsv(str(tmp_path / 'vocab.tsv'), {
            'Gene': ['BRCA1', 'N/A', 'null', 'None', '<NA>', '#N/A', '--', '?'],
            'logFC': [1.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            'pval': [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        assert set(df['ID']) == {'BRCA1'}
        assert df.attrs['dropped_unidentified_rows'] == 7

    def test_missing_symbol_row_cannot_be_significant(self, tmp_path):
        """A symbol-less row used to enter every Fisher test as 'NAN'."""
        path = self._write_tsv(str(tmp_path / 'sig.tsv'), {
            'Gene': ['BRCA1', None],
            'logFC': [1.5, 9.0],
            'pval': [0.001, 0.0001],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'pval')
        result, stats = process_gene_expression(df, logfc_threshold=1.0)
        assert stats['significant_genes'] == 1
        assert 'NAN' not in set(result['ID'])


class TestGeneUniverse:
    """Issue #132: which rows count as the gene universe.

    Read on `padj`, the loader's original drop-any-NaN rule removes every gene
    DESeq2's independent filtering withheld an adjusted p-value from. That
    selects against exactly the expression profile an induced stress response
    has — low at baseline, strongly up on treatment — so those genes cannot be
    counted in the Fisher background either.
    """

    def _write(self, tmp_path, rows):
        path = tmp_path / 'dge.csv'
        pd.DataFrame(rows).to_csv(path, index=False)
        return str(path)

    def test_testable_is_the_default_and_drops_missing_pvalues(self, tmp_path):
        path = self._write(tmp_path, {
            'Gene': ['BRCA1', 'HMOX1'],
            'logFC': [1.0, 4.47],
            'padj': [0.01, None],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'padj')
        assert set(df['ID']) == {'BRCA1'}

    def test_measured_keeps_a_gene_with_no_adjusted_pvalue(self, tmp_path):
        path = self._write(tmp_path, {
            'Gene': ['BRCA1', 'HMOX1'],
            'logFC': [1.0, 4.47],
            'padj': [0.01, None],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'padj', universe='measured')
        assert set(df['ID']) == {'BRCA1', 'HMOX1'}
        assert df.attrs['untestable_rows'] == 1

    def test_a_gene_with_no_pvalue_is_never_significant(self, tmp_path):
        """It belongs in the background, not among the hits — the asymmetry
        Fisher's exact test needs. NaN <= cutoff is False, so this needs no
        special handling, but it must stay true."""
        path = self._write(tmp_path, {
            'Gene': ['BRCA1', 'HMOX1'],
            'logFC': [1.0, 4.47],
            'padj': [0.01, None],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'padj', universe='measured')
        processed, stats = process_gene_expression(df, 0.0, pval_threshold=0.05)
        assert stats['total_genes'] == 2
        assert stats['significant_genes'] == 1
        assert stats['untestable_genes'] == 1
        assert not processed.loc[processed['ID'] == 'HMOX1', 'significant'].iloc[0]

    def test_a_missing_fold_change_drops_the_row_under_both_rules(self, tmp_path):
        """DESeq2 leaves log2FoldChange NA only for genes with no counts at all.
        An unmeasurable gene in the background deflates every enrichment."""
        path = self._write(tmp_path, {
            'Gene': ['BRCA1', 'ZEROCOUNT'],
            'logFC': [1.0, None],
            'padj': [0.01, None],
        })
        for universe in ('testable', 'measured'):
            df = load_and_validate_data(path, 'Gene', 'logFC', 'padj', universe=universe)
            assert set(df['ID']) == {'BRCA1'}, universe

    def test_duplicate_symbols_combine_over_finite_pvalues_only(self, tmp_path):
        """A symbol testable in one row must not be demoted to untestable
        because it shares a name with a row that is not."""
        path = self._write(tmp_path, {
            'Gene': ['MT1F', 'MT1F'],
            'logFC': [1.0, 2.0],
            'padj': [0.002, None],
        })
        df = load_and_validate_data(path, 'Gene', 'logFC', 'padj', universe='measured')
        processed, _ = process_gene_expression(df, 0.0, pval_threshold=0.05)
        assert len(processed) == 1
        assert processed['pval'].notna().all()
        assert processed['significant'].iloc[0]

    def test_unknown_universe_is_rejected(self, tmp_path):
        path = self._write(tmp_path, {'Gene': ['BRCA1'], 'logFC': [1.0], 'padj': [0.01]})
        with pytest.raises(DataValidationError, match='Unknown gene universe'):
            load_and_validate_data(path, 'Gene', 'logFC', 'padj', universe='everything')
