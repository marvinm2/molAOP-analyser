"""
Regression tests for issue #69 — silent failure on the wrong ID / p-value column.

A DESeq2 export carries both an Ensembl 'GeneID' and an HGNC 'GeneSymbol', plus
raw 'pvalue' and adjusted 'padj'. The batch wizard used to auto-fill GeneID and
pvalue, and the run completed with a plausible-looking but empty result.
"""

import pandas as pd
import pytest

from services.gene_id_validator import gene_id_validator
from services.column_detector import column_detector
from services.enrichment_service import (
    assess_background_overlap,
    MIN_BACKGROUND_OVERLAP_FRACTION,
)


def _deseq2_frame():
    """A standard DESeq2 export layout."""
    return pd.DataFrame({
        'GeneID': ['ENSG00000000003', 'ENSG00000000005',
                   'ENSG00000000419', 'ENSG00000000457'],
        'GeneSymbol': ['TSPAN6', 'TNMD', 'DPM1', 'SCYL3'],
        'baseMean': [100.0, 5.0, 300.0, 50.0],
        'log2FoldChange': [1.2, -0.3, 2.4, 0.1],
        'lfcSE': [0.2, 0.4, 0.3, 0.5],
        'stat': [6.0, -0.75, 8.0, 0.2],
        'pvalue': [1e-9, 0.45, 1e-14, 0.83],
        'padj': [1e-7, 0.62, 1e-12, 0.91],
    })


@pytest.mark.unit
class TestIdentifierClassification:
    """The HGNC pattern is a catch-all and must not shadow specific types."""

    @pytest.mark.parametrize('ids,expected', [
        (['ENSG00000000003', 'ENSG00000000005', 'ENSG00000000419'], 'Ensembl'),
        (['ENSG00000000003.15', 'ENSG00000000005.6'], 'Ensembl'),
        (['TSPAN6', 'TNMD', 'DPM1'], 'HGNC'),
        (['7105', '64102', '8813'], 'NCBI'),
        (['NM_001005484', 'NM_015658'], 'RefSeq'),
        (['1007_s_at', '1053_at'], 'Probe'),
    ])
    def test_types_are_not_swallowed_by_hgnc(self, ids, expected):
        analysis = gene_id_validator.analyze_gene_column(pd.Series(ids))
        assert analysis.primary_type == expected
        assert analysis.confidence == pytest.approx(1.0)

    def test_only_symbols_are_enrichment_matchable(self):
        """Ensembl is a *supported* type but cannot be matched without translation."""
        assert 'HGNC' in gene_id_validator.ENRICHMENT_MATCHABLE_TYPES
        assert 'Ensembl' not in gene_id_validator.ENRICHMENT_MATCHABLE_TYPES


@pytest.mark.unit
class TestColumnDetection:

    def test_prefers_symbol_column_over_ensembl(self):
        """The reported bug: 'GeneID' won on column ordering."""
        result = column_detector.detect_columns(_deseq2_frame())
        assert result.best_gene_id.column_name == 'GeneSymbol'

    def test_ensembl_column_still_offered(self):
        """Ranked lower, not hidden — the user may have a reason to pick it."""
        result = column_detector.detect_columns(_deseq2_frame())
        names = [m.column_name for m in result.gene_id_suggestions]
        assert 'GeneID' in names
        by_name = {m.column_name: m.confidence for m in result.gene_id_suggestions}
        assert by_name['GeneSymbol'] > by_name['GeneID']

    def test_symbol_only_file_is_unaffected(self):
        """No regression when there is nothing to disambiguate."""
        df = _deseq2_frame().drop(columns=['GeneID'])
        assert column_detector.detect_columns(df).best_gene_id.column_name == 'GeneSymbol'

    def test_raw_and_adjusted_are_both_detected(self):
        result = column_detector.detect_columns(_deseq2_frame())
        assert result.best_pvalue.column_name == 'pvalue'
        assert result.best_pvalue_adj.column_name == 'padj'


@pytest.mark.unit
class TestBatchSuggestionPayload:
    """The batch wizard has one p-value field, so it must choose deliberately."""

    def _suggestions(self, df):
        from app import _batch_column_suggestions
        return _batch_column_suggestions(column_detector.detect_columns(df))

    def test_prefers_adjusted_pvalue(self):
        s = self._suggestions(_deseq2_frame())
        assert s['pval_col'] == 'padj'
        assert s['pval_is_adjusted'] is True
        assert s['pval_raw_col'] == 'pvalue'

    def test_falls_back_to_raw_when_no_adjusted_column(self):
        df = _deseq2_frame().drop(columns=['padj'])
        s = self._suggestions(df)
        assert s['pval_col'] == 'pvalue'
        assert s['pval_is_adjusted'] is False
        assert s['pval_adj_col'] is None

    def test_suggests_the_symbol_column(self):
        assert self._suggestions(_deseq2_frame())['id_col'] == 'GeneSymbol'


@pytest.mark.unit
class TestBackgroundOverlap:
    """Last line of defence: fail loudly when the IDs cannot match."""

    REFS = {
        'KE:1': {'TP53', 'BRCA1', 'DPM1', 'TSPAN6', 'SCYL3'},
        'KE:2': {'AKT1', 'MAPK8'},
    }

    def test_ensembl_upload_is_flagged(self):
        result = assess_background_overlap(
            {'ENSG00000000003', 'ENSG00000000005', 'ENSG00000000419'}, self.REFS)
        assert result['matched'] == 0
        assert result['is_suspect'] is True

    def test_symbol_upload_is_not_flagged(self):
        result = assess_background_overlap(
            {'TP53', 'BRCA1', 'AKT1', 'NOTAGENE1', 'NOTAGENE2'}, self.REFS)
        assert result['matched'] == 3
        assert result['is_suspect'] is False

    def test_case_and_whitespace_are_normalised(self):
        result = assess_background_overlap({' tp53 ', 'brca1'}, self.REFS)
        assert result['matched'] == 2
        assert result['is_suspect'] is False

    def test_threshold_boundary(self):
        """Just above the threshold must not be flagged."""
        universe = {f'GENE{i}' for i in range(100)}
        refs = {'KE:1': universe}
        # 6 of 100 uploaded genes match => 6%, above the 5% floor.
        genes = {f'GENE{i}' for i in range(6)} | {f'MISS{i}' for i in range(94)}
        result = assess_background_overlap(genes, refs)
        assert result['fraction'] > MIN_BACKGROUND_OVERLAP_FRACTION
        assert result['is_suspect'] is False

    def test_empty_inputs_do_not_flag(self):
        """An empty upload or empty reference set is a different problem."""
        assert assess_background_overlap(set(), self.REFS)['is_suspect'] is False
        assert assess_background_overlap({'TP53'}, {})['is_suspect'] is False


@pytest.mark.unit
class TestAnalyzeRouteBehaviour:
    """The two ways a wrong ID column reaches the user (issue #69)."""

    def _post(self, flask_client, ids):
        import io
        df = pd.DataFrame({
            'GENE': ids,
            'logFC': [2.5] * len(ids),
            'PVal': [0.0001] * len(ids),
        })
        buf = io.BytesIO(df.to_csv(index=False).encode())
        flask_client.post('/preview', data={
            'gene_file': (buf, 'x.csv'), 'id_col': 'GENE',
            'fc_col': 'logFC', 'pval_col': 'PVal', 'columns_confirmed': 'true',
        }, content_type='multipart/form-data')
        return flask_client.post('/analyze', data={
            'filename': 'x.csv', 'id_column': 'GENE', 'fc_column': 'logFC',
            'pval_column': 'PVal', 'aop_selection': 'AOP:472',
            'logfc_threshold': '1.0', 'method': 'ora',
            'resources': ['WikiPathways', 'GO_BP', 'Reactome'],
            'min_confidence': 'all',
        })

    def test_total_mismatch_gives_an_actionable_error(self, flask_client):
        """Was 'Analysis error: please check your input data and parameters'."""
        response = self._post(flask_client, [f'ENSG{i:011d}' for i in range(300)])
        assert response.status_code == 400
        body = response.get_data(as_text=True)
        assert 'matched the reference gene sets' in body
        assert 'HGNC gene symbols' in body

    def test_partial_mismatch_warns_on_the_results_page(self, flask_client):
        """Enough overlap to produce results, too little for them to mean anything."""
        from app import load_cached_reference_sets
        reference_sets, _ = load_cached_reference_sets(
            ['WikiPathways', 'GO_BP', 'Reactome'], min_confidence='all')
        aop_kes = ['KE:1115', 'KE:1194', 'KE:1392', 'KE:149', 'KE:177',
                   'KE:1825', 'KE:373', 'KE:1097', 'KE:759']
        symbols = sorted({g for k in aop_kes for g in reference_sets.get(k, set())})
        if len(symbols) < 20:
            pytest.skip('reference sets unavailable in this environment')

        ids = [f'ENSG{i:011d}' for i in range(580)] + symbols[:20]
        response = self._post(flask_client, ids)
        assert response.status_code == 200
        assert 'probably not meaningful' in response.get_data(as_text=True)
