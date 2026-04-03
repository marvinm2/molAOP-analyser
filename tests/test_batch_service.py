"""Tests for the batch analysis service."""
import os
import tempfile
import pytest
import pandas as pd
from services.batch_service import (
    parse_cisplatin_filename,
    create_batch_upload_dir,
    cleanup_batch_upload_dir,
    validate_batch_columns,
    harmonise_backgrounds,
)


class TestParseCisplatinFilename:
    """Tests for cisplatin filename parsing."""

    def test_standard_filename(self):
        result = parse_cisplatin_filename('CSP_4hr_0.1uM.csv')
        assert result['timepoint'] == '4hr'
        assert result['dose'] == '0.1uM'
        assert result['condition_label'] == 'Cisplatin 4hr 0.1uM'

    def test_long_timepoint(self):
        result = parse_cisplatin_filename('CSP_72hr_50uM.csv')
        assert result['timepoint'] == '72hr'
        assert result['dose'] == '50uM'

    def test_non_matching_filename(self):
        result = parse_cisplatin_filename('random_file.csv')
        assert result['condition_label'] == 'random_file.csv'
        assert result['dose'] == ''
        assert result['timepoint'] == ''

    def test_partial_match_no_extension(self):
        result = parse_cisplatin_filename('CSP_4hr_0.1uM')
        assert result['dose'] == ''  # no .csv extension → no match

    def test_path_prefix_stripped(self):
        """Regex searches inside string, so path prefix is OK."""
        result = parse_cisplatin_filename('Cisplatin_Kidney/CSP_8hr_5uM.csv')
        assert result['timepoint'] == '8hr'
        assert result['dose'] == '5uM'


class TestBatchUploadDir:
    """Tests for batch upload directory management."""

    def test_create_and_cleanup(self):
        with tempfile.TemporaryDirectory() as tmp:
            batch_uuid, batch_dir = create_batch_upload_dir(tmp)
            assert os.path.isdir(batch_dir)
            assert batch_uuid in batch_dir

            cleanup_batch_upload_dir(tmp, batch_uuid)
            assert not os.path.isdir(batch_dir)

    def test_cleanup_nonexistent_is_noop(self):
        with tempfile.TemporaryDirectory() as tmp:
            cleanup_batch_upload_dir(tmp, 'nonexistent-uuid')  # should not raise


class TestValidateBatchColumns:
    """Tests for batch column validation."""

    def _write_csv(self, path, columns, rows=3):
        df = pd.DataFrame({c: range(rows) for c in columns})
        df.to_csv(path, index=False)
        return path

    def test_valid_batch(self, tmp_path):
        f1 = self._write_csv(str(tmp_path / 'a.csv'), ['Gene', 'logFC', 'pval'])
        f2 = self._write_csv(str(tmp_path / 'b.csv'), ['Gene', 'logFC', 'pval'])
        infos = [
            {'filepath': f1, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
            {'filepath': f2, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
        ]
        ok, err = validate_batch_columns(infos)
        assert ok
        assert err == ''

    def test_missing_column(self, tmp_path):
        f1 = self._write_csv(str(tmp_path / 'a.csv'), ['Gene', 'logFC'])  # missing pval
        infos = [{'filepath': f1, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'}]
        ok, err = validate_batch_columns(infos)
        assert not ok
        assert 'missing columns' in err.lower()

    def test_column_mismatch_across_files(self, tmp_path):
        f1 = self._write_csv(str(tmp_path / 'a.csv'), ['Gene', 'logFC', 'pval'])
        f2 = self._write_csv(str(tmp_path / 'b.csv'), ['Gene', 'logFC', 'pval'])
        infos = [
            {'filepath': f1, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
            {'filepath': f2, 'id_col': 'Gene', 'fc_col': 'FC', 'pval_col': 'pval'},  # mismatch
        ]
        ok, err = validate_batch_columns(infos)
        assert not ok
        assert 'mismatch' in err.lower()

    def test_empty_list(self):
        ok, err = validate_batch_columns([])
        assert not ok


class TestHarmoniseBackgrounds:
    """Tests for gene universe intersection across batch files."""

    def _write_gene_csv(self, path, genes, fc_vals=None, pvals=None):
        n = len(genes)
        df = pd.DataFrame({
            'Gene': genes,
            'logFC': fc_vals or [1.0] * n,
            'pval': pvals or [0.05] * n,
        })
        df.to_csv(path, index=False)
        return path

    def test_intersection_of_two_files(self, tmp_path):
        f1 = self._write_gene_csv(str(tmp_path / 'a.csv'), ['BRCA1', 'TP53', 'EGFR'])
        f2 = self._write_gene_csv(str(tmp_path / 'b.csv'), ['BRCA1', 'TP53', 'MYC'])
        infos = [
            {'filepath': f1, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
            {'filepath': f2, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
        ]
        harmonised, per_file = harmonise_backgrounds(infos)
        assert harmonised == {'BRCA1', 'TP53'}
        assert per_file[f1] == 3
        assert per_file[f2] == 3

    def test_no_overlap(self, tmp_path):
        f1 = self._write_gene_csv(str(tmp_path / 'a.csv'), ['BRCA1'])
        f2 = self._write_gene_csv(str(tmp_path / 'b.csv'), ['MYC'])
        infos = [
            {'filepath': f1, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
            {'filepath': f2, 'id_col': 'Gene', 'fc_col': 'logFC', 'pval_col': 'pval'},
        ]
        harmonised, _ = harmonise_backgrounds(infos)
        assert harmonised == set()
