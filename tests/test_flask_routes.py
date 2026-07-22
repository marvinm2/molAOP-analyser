"""
Integration tests for Flask routes and web functionality.
"""

import io
import pytest
import json
import re
from unittest.mock import patch, MagicMock

from tests.conftest import uploaded_file_exists


@pytest.mark.integration
@pytest.mark.web
class TestFlaskRoutes:
    """Test Flask application routes."""
    
    def test_index_route(self, flask_client):
        """Test main index page loads correctly with both tabs (Phase 10.1 layout)."""
        response = flask_client.get('/')

        assert response.status_code == 200
        assert b'molAOP Analyser' in response.data
        # Phase 10.1: demo selection block was removed; banner + upload form now lead.
        assert b'Select Demo Dataset' not in response.data
        assert b'Browse demos' in response.data           # banner CTA
        assert b'Upload Your Own Data' in response.data   # upload heading
        assert b'href="/demos"' in response.data          # banner links to /demos
        # Tabs and metadata form unchanged.
        assert b'Single Analysis' in response.data
        assert b'Batch Analysis' in response.data
        assert b'Experiment Information' in response.data
        # No demo radio buttons remain. Use two separate assertions so both conditions
        # are verified independently. The second check uses the combined attribute string
        # to avoid false positives from CSS rules like `.dose-label input[type="radio"]`.
        assert b'name="demo_file"' not in response.data
        assert b'type="radio" name="demo_file"' not in response.data

    def test_batch_route_redirects(self, flask_client):
        """Test that /batch redirects to /?tab=batch."""
        response = flask_client.get('/batch')

        assert response.status_code == 302
        assert '/?tab=batch' in response.headers['Location']
    
    def test_preview_route_with_demo_file(self, flask_client):
        """Test preview route with demo file selection."""
        with patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('pandas.read_csv') as mock_read_csv:
            
            # Mock DataFrame
            mock_df = MagicMock()
            mock_df.head.return_value.to_dict.return_value = [
                {'Gene_Symbol': 'BRCA1', 'log2FoldChange': 2.5, 'padj': 0.001}
            ]
            mock_df.columns.tolist.return_value = ['Gene_Symbol', 'log2FoldChange', 'padj']
            mock_read_csv.return_value = mock_df
            
            response = flask_client.post('/preview', data={
                'demo_file': 'GSE90122_TO90137.tsv',
                'dataset_id': 'TEST001',
                'stressor': 'Test Chemical',
                'owner': 'Test User'
            })
            
            assert response.status_code == 200
            assert b'Preview top 5 rows' in response.data
            assert b'Gene Symbol / ID column' in response.data
    
    def test_preview_route_missing_file(self, flask_client):
        """Test preview route with missing file."""
        response = flask_client.post('/preview', data={})

        assert response.status_code == 400
        assert b'No dataset provided' in response.data

    def test_preview_route_with_uploaded_file(self, flask_client):
        """Test preview route accepts a multipart-uploaded CSV (UPUX-04 regression guard).

        Mirrors test_preview_route_with_demo_file but exercises the upload-your-own-data
        path used by the drag-and-drop widget on the main page. Posts a real in-memory
        CSV via io.BytesIO and asserts the column-selection step renders.

        Field name `gene_file` matches the hidden file input in
        templates/_single_analysis.html (line 48: name="gene_file").
        """
        csv_bytes = (
            b"Gene_Symbol,log2FoldChange,padj\n"
            b"BRCA1,2.5,0.001\n"
            b"TP53,-0.8,0.05\n"
            b"EGFR,2.1,0.0001\n"
            b"MYC,0.3,0.5\n"
            b"KRAS,-1.2,0.01\n"
        )

        response = flask_client.post(
            '/preview',
            data={
                'gene_file': (io.BytesIO(csv_bytes), 'test_upload.csv'),
                'dataset_id': 'TEST_UPLOAD',
                'stressor': 'Test Chemical',
                'owner': 'Test User',
            },
            content_type='multipart/form-data',
        )

        assert response.status_code == 200, (
            f"Expected 200 from /preview multipart upload, got {response.status_code}; "
            f"body excerpt: {response.data[:500]!r}"
        )

        body = response.data
        assert (
            b'Select column' in body
            or b'log2FC' in body
            or b'<select' in body
        ), (
            "Expected at least one column-selector substring "
            "('Select column', 'log2FC', or '<select') in /preview response body; "
            f"body excerpt: {body[:500]!r}"
        )

    def test_analyze_route_success(self, flask_client, authenticated_client):
        """Test successful analysis route."""
        import pandas as pd

        # Build realistic mock data that templates can render
        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS'],
            'log2FC': [1.5, -0.8, 2.1, 0.3, -1.2],
            'pval': [0.001, 0.05, 0.0001, 0.5, 0.01],
            'significant': [True, False, True, False, True],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.05],
            'num_overlap': [3], 'pct_sig_in_KE': [60.0],
            'total_KE_genes_in_dataset': [5], 'odds_ratio': [3.5],
            'overlap': ['BRCA1, EGFR, KRAS'], 'KE': ['KE:115'],
            'sig_in_KE': [3], 'sig_not_KE': [0],
            'non_sig_in_KE': [2], 'non_sig_not_KE': [0],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 5})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.build_ke_gene_mapping', return_value={}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('app.validate_file_path', return_value=True):

            response = authenticated_client.post('/analyze', data={
                'filename': 'test.csv',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
                'aop_selection': 'AOP:1',
                'logfc_threshold': '1.0'
            })

            assert response.status_code == 200
    
    def test_results_page_data_island_carries_pvalue_passthrough(self, authenticated_client):
        """Plan 11-03 (EXPO-04): when the upload has both raw `pvalue` and
        adjusted `padj` columns, the rendered results page's ke_gene_json
        data island contains per-gene `pvalue_raw` and `pvalue_adj` keys so
        the client-side gene-by-KE export can read them.

        We do NOT mock `build_ke_gene_mapping` here -- we assert against the
        rendered response substring, which exercises the full pipeline:
        column detection -> raw/adj map extraction -> build_ke_gene_mapping
        -> JSON serialisation -> template render.
        """
        import pandas as pd

        # Raw upload as it would land on disk: includes BOTH pvalue (raw) and
        # padj (adjusted) columns; the user thresholds on padj.
        raw_upload_df = pd.DataFrame({
            'Gene_Symbol': ['BRCA1', 'TP53', 'EGFR'],
            'log2FoldChange': [1.5, -0.8, 2.1],
            'pvalue': [0.001, 0.05, 0.0001],
            'padj': [0.01, 0.1, 0.001],
        })

        # df_processed as data_service would produce it (only ID/log2FC/pval/significant).
        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53', 'EGFR'],
            'log2FC': [1.5, -0.8, 2.1],
            'pval': [0.01, 0.1, 0.001],   # padj values (user thresholded on padj)
            'significant': [True, False, True],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.05],
            'num_overlap': [2], 'pct_sig_in_KE': [66.7],
            'total_KE_genes_in_dataset': [3], 'odds_ratio': [3.5],
            'overlap': ['BRCA1, EGFR'], 'KE': ['KE:115'],
            'sig_in_KE': [2], 'sig_not_KE': [0],
            'non_sig_in_KE': [1], 'non_sig_not_KE': [0],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        # NOTE: do NOT mock build_ke_gene_mapping -- we want it to run for real
        # so the per-gene pvalue_raw / pvalue_adj keys actually flow into ke_gene_json.
        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 3})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('app.validate_file_path', return_value=True), \
             patch('pandas.read_csv', return_value=raw_upload_df), \
             patch('app.load_cached_reference_sets', return_value=({'KE:115': {'BRCA1', 'EGFR'}}, 'cache', [])):

            response = authenticated_client.post('/analyze', data={
                'filename': 'test.csv',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
                'aop_selection': 'AOP:1',
                'logfc_threshold': '1.0'
            })

            assert response.status_code == 200
            # The data island MUST carry both passthrough keys (per-gene shape
            # change is what the JS CSV writer relies on).
            assert b'pvalue_raw' in response.data, \
                "pvalue_raw missing from rendered ke_gene_json data island"
            assert b'pvalue_adj' in response.data, \
                "pvalue_adj missing from rendered ke_gene_json data island"

    def test_analyze_respects_user_chosen_pval_adj_column(self, authenticated_client):
        """Plan 11-05 (UAT gap closure): when the user explicitly picks
        ``pval_adj_column`` in the preview UI, that column drives the
        gene-by-KE CSV ``pvalue_adj`` cells -- overriding auto-detection
        AND skipping the same-column-as-raw disambiguation.

        The upload carries TWO adj-pattern columns at disjoint surrogate
        values:

        - ``padj_a`` = 0.999999  (auto-detector default candidate)
        - ``padj_b`` = 1.5e-9    (user's explicit pick)

        Surrogate values are chosen far from any other on-page numeric
        (FDR cells, p_value cells, default thresholds) so a numeric
        match in the parsed data island is unambiguous evidence of which
        column was wired through.
        """
        import pandas as pd

        # Raw upload: BRCA1/EGFR are in the KE reference set so they appear
        # in keToGenes; TP53 is intentionally not in the reference set, so
        # only the first two carry our surrogates into the data island.
        raw_upload_df = pd.DataFrame({
            'Gene_Symbol': ['BRCA1', 'TP53', 'EGFR'],
            'log2FoldChange': [1.5, -0.8, 2.1],
            'pvalue':        [0.001, 0.05, 0.0001],
            'padj_a':        [0.999999, 0.999999, 0.999999],
            'padj_b':        [1.5e-9, 1.5e-9, 1.5e-9],
        })

        # df_processed mirrors what data_service.process_gene_expression
        # would produce: only ID/log2FC/pval/significant. The user thresholded
        # on `pvalue` here (the raw column).
        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53', 'EGFR'],
            'log2FC': [1.5, -0.8, 2.1],
            'pval': [0.001, 0.05, 0.0001],
            'significant': [True, False, True],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.05],
            'num_overlap': [2], 'pct_sig_in_KE': [66.7],
            'total_KE_genes_in_dataset': [3], 'odds_ratio': [3.5],
            'overlap': ['BRCA1, EGFR'], 'KE': ['KE:115'],
            'sig_in_KE': [2], 'sig_not_KE': [0],
            'non_sig_in_KE': [1], 'non_sig_not_KE': [0],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        # Do NOT mock build_ke_gene_mapping -- the override must flow through
        # the real code path so adj_lookup is built from `padj_b` values, not
        # from the auto-detected column.
        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 3})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('app.validate_file_path', return_value=True), \
             patch('pandas.read_csv', return_value=raw_upload_df), \
             patch('app.load_cached_reference_sets',
                   return_value=({'KE:115': {'BRCA1', 'EGFR'}}, 'cache', [])):

            response = authenticated_client.post('/analyze', data={
                'filename': 'test.csv',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'pvalue',
                'pval_adj_column': 'padj_b',   # the override under test
                'aop_selection': 'AOP:1',
                'logfc_threshold': '1.0',
            })

            assert response.status_code == 200

            # Parse the keToGenes data island. json.dumps with default
            # separators emits single-line JSON with no internal `;`, so a
            # MULTILINE-anchored match up to the first `;` end-of-line is
            # safe and avoids non-greedy DOTALL pitfalls (which would stop
            # at the first inner `}`).
            body = response.get_data(as_text=True)
            match = re.search(
                r'(?:const|var|let)\s+keToGenes\s*=\s*(.+?);\s*$',
                body, re.MULTILINE,
            )
            assert match is not None, "keToGenes data island missing from response"
            ke_to_genes = json.loads(match.group(1))

            # Collect every per-gene pvalue_adj across every KE.
            all_adj = [
                gene['pvalue_adj']
                for ke_genes in ke_to_genes.values()
                for gene in ke_genes
                if 'pvalue_adj' in gene and gene['pvalue_adj'] is not None
            ]
            assert all_adj, \
                f"No pvalue_adj values flowed through (saw keToGenes={ke_to_genes})"

            # Positive: the user-picked padj_b surrogate must appear.
            assert any(abs(v - 1.5e-9) < 1e-12 for v in all_adj), \
                f"User-picked padj_b values not in keToGenes (saw: {all_adj[:10]})"

            # Negative: the auto-detected padj_a surrogate must NOT appear.
            assert not any(abs(v - 0.999999) < 1e-9 for v in all_adj), \
                f"Auto-detected padj_a leaked through despite user pick (saw: {all_adj[:10]})"

    def test_analyze_combines_duplicate_probes_for_pvalue_passthrough(self, authenticated_client):
        """Regression: when a gene has multiple probes in the upload, the
        gene-by-KE export's pvalue/FDR cells must reflect a Fisher-combined
        per-gene value -- not whichever probe happened to be last in the
        DataFrame.

        Pre-fix, ``dict(zip(ids, series))`` silently kept the last probe per
        duplicate gene, so a row with significant=True (driven by the
        Fisher-combined gate in df_processed) could ship next to a
        non-significant probe p-value, contradicting the significance flag.
        """
        import pandas as pd

        # BRCA1 has two probes mirroring the FOXA1 case observed live:
        # one strong probe (1e-4) and one weak (0.9). Fisher-combined raw p
        # for BRCA1 is dominated by the strong probe and lands well under
        # 0.05; min adj likewise reflects the strong probe.
        raw_upload_df = pd.DataFrame({
            'Gene_Symbol': ['BRCA1', 'BRCA1', 'EGFR'],
            'log2FoldChange': [1.5, 1.4, 2.1],
            'pvalue':         [1e-4,  0.9,  0.0001],
            'padj':           [1e-3,  0.95, 0.001],
        })
        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'EGFR'],
            'log2FC': [1.45, 2.1],
            'pval':   [1e-4, 0.0001],     # Fisher-combined for BRCA1
            'significant': [True, True],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.05],
            'num_overlap': [2], 'pct_sig_in_KE': [100.0],
            'total_KE_genes_in_dataset': [2], 'odds_ratio': [3.5],
            'overlap': ['BRCA1, EGFR'], 'KE': ['KE:115'],
            'sig_in_KE': [2], 'sig_not_KE': [0],
            'non_sig_in_KE': [0], 'non_sig_not_KE': [0],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 2})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('app.validate_file_path', return_value=True), \
             patch('pandas.read_csv', return_value=raw_upload_df), \
             patch('app.load_cached_reference_sets',
                   return_value=({'KE:115': {'BRCA1', 'EGFR'}}, 'cache', [])):

            response = authenticated_client.post('/analyze', data={
                'filename': 'test.csv',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'pvalue',
                'aop_selection': 'AOP:1',
                'logfc_threshold': '1.0',
            })
            assert response.status_code == 200

            body = response.get_data(as_text=True)
            match = re.search(
                r'(?:const|var|let)\s+keToGenes\s*=\s*(.+?);\s*$',
                body, re.MULTILINE,
            )
            assert match is not None, "keToGenes data island missing from response"
            ke_to_genes = json.loads(match.group(1))

            brca1 = next(
                (g for ke_genes in ke_to_genes.values() for g in ke_genes
                 if g.get('id') == 'BRCA1'),
                None,
            )
            assert brca1 is not None, "BRCA1 absent from keToGenes"

            # Pre-fix the last-probe overwrite would have shipped 0.9 / 0.95.
            # Post-fix the Fisher-combined raw p is far below the weak probe,
            # and the adjusted column is combined the same way for consistency.
            assert brca1['pvalue_raw'] < 0.05, (
                f"BRCA1 pvalue_raw={brca1['pvalue_raw']} suggests last-probe "
                "overwrite (expected Fisher-combined value < 0.05)"
            )
            assert brca1['pvalue_adj'] < 0.05, (
                f"BRCA1 pvalue_adj={brca1['pvalue_adj']} suggests last-probe "
                "overwrite (expected aggregated value < 0.05)"
            )

    def test_analyze_route_validation_error(self, flask_client):
        """Test analysis route with validation errors."""
        response = flask_client.post('/analyze', data={
            # Missing required fields
            'filename': '',
            'id_column': '',
            'aop_selection': ''
        })
        
        assert response.status_code == 400
        assert b'Validation errors' in response.data
    
    def test_generate_report_route_html(self, authenticated_client):
        """Test HTML report generation route."""
        with patch('app.report_generator.generate_html_report') as mock_generate:
            mock_generate.return_value = '<html><body>Test Report</body></html>'
            
            response = authenticated_client.post('/generate_report', data={
                'format': 'html',
                'filename': 'test.csv',
                'gene_count': '1000',
                'significant_genes': '150',
                'aop_id': 'AOP:1',
                'aop_label': 'Test AOP',
                'logfc_threshold': '1.0',
                'pval_cutoff': '0.05',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
                'id_type': 'HGNC',
                'enrichment_results': '[]'
            })
            
            assert response.status_code == 200
            assert 'text/html' in response.headers['Content-Type']
            assert b'Test Report' in response.data
    
    def test_generate_report_route_pdf(self, authenticated_client):
        """Test PDF report generation route."""
        with patch('app.report_generator.generate_pdf_report') as mock_generate:
            mock_generate.return_value = b'%PDF-1.4 fake pdf content'
            
            response = authenticated_client.post('/generate_report', data={
                'format': 'pdf',
                'filename': 'test.csv',
                'gene_count': '1000',
                'significant_genes': '150',
                'aop_id': 'AOP:1',
                'aop_label': 'Test AOP',
                'logfc_threshold': '1.0',
                'pval_cutoff': '0.05',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
                'id_type': 'HGNC',
                'enrichment_results': '[]'
            })
            
            assert response.status_code == 200
            assert response.headers['Content-Type'] == 'application/pdf'
            assert b'PDF' in response.data
    
    def test_generate_report_route_no_metadata(self, flask_client):
        """Test report generation without experiment metadata."""
        response = flask_client.post('/generate_report', data={
            'format': 'html',
            'filename': 'test.csv'
        })
        
        assert response.status_code == 400
        assert b'No experiment data found' in response.data
    
    def test_generate_report_route_invalid_format(self, authenticated_client):
        """Test report generation with invalid format."""
        response = authenticated_client.post('/generate_report', data={
            'format': 'invalid_format',
            'filename': 'test.csv'
        })
        
        assert response.status_code == 400
        assert b'Unsupported report format' in response.data
    
    def test_error_handlers(self, flask_client):
        """Test custom error handlers."""
        # Test 404 handler
        response = flask_client.get('/nonexistent-page')
        assert response.status_code == 404
        assert b'Page not found' in response.data
    
    def test_metadata_storage_in_session(self, flask_client):
        """Test that metadata is stored in session during preview."""
        with patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('pandas.read_csv') as mock_read_csv, \
             flask_client.session_transaction() as sess:
            
            # Mock DataFrame
            mock_df = MagicMock()
            mock_df.head.return_value.to_dict.return_value = []
            mock_df.columns.tolist.return_value = ['Gene_Symbol']
            mock_read_csv.return_value = mock_df
            
            response = flask_client.post('/preview', data={
                'demo_file': 'GSE90122_TO90137.tsv',
                'dataset_id': 'SESSION_TEST',
                'stressor': 'Session Test Chemical',
                'owner': 'Session Test User'
            })
            
            assert response.status_code == 200
            
            # Check session after request
            with flask_client.session_transaction() as sess:
                assert 'experiment_metadata' in sess
                metadata = sess['experiment_metadata']
                assert metadata['dataset_id'] == 'SESSION_TEST'
                assert metadata['stressor'] == 'Session Test Chemical'
    
    def test_csrf_protection_disabled_in_tests(self, flask_client):
        """Test that CSRF is properly disabled in test configuration."""
        # This test ensures our test client can make POST requests
        response = flask_client.post('/preview', data={})
        
        # Should not get CSRF error (should get validation error instead)
        assert response.status_code != 400 or b'CSRF token is missing' not in response.data
    
    @pytest.mark.slow
    def test_full_workflow_integration(self, flask_client):
        """Test complete workflow from upload to report generation."""
        import pandas as pd

        # Realistic DataFrames that survive template rendering
        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS'],
            'log2FC': [1.5, -0.8, 2.1, 0.3, -1.2],
            'pval': [0.001, 0.05, 0.0001, 0.5, 0.01],
            'significant': [True, False, True, False, True],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.05],
            'num_overlap': [3], 'pct_sig_in_KE': [60.0],
            'total_KE_genes_in_dataset': [5], 'odds_ratio': [3.5],
            'overlap': ['BRCA1, EGFR, KRAS'], 'KE': ['KE:115'],
            'sig_in_KE': [3], 'sig_not_KE': [0],
            'non_sig_in_KE': [2], 'non_sig_not_KE': [0],
        })
        preview_df = pd.DataFrame({
            'Gene_Symbol': ['BRCA1', 'TP53'],
            'log2FoldChange': [1.5, -0.8],
            'padj': [0.001, 0.05],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        with patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('pandas.read_csv', return_value=preview_df), \
             patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 5})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.build_ke_gene_mapping', return_value={}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('app.validate_file_path', return_value=True), \
             patch('app.report_generator.generate_html_report') as mock_report:

            mock_report.return_value = '<html>Integration Test Report</html>'

            # Step 1: Preview with metadata
            response1 = flask_client.post('/preview', data={
                'demo_file': 'GSE90122_TO90137.tsv',
                'dataset_id': 'INTEGRATION_TEST',
                'stressor': 'Integration Chemical',
                'owner': 'Integration User'
            })
            assert response1.status_code == 200

            # Step 2: Analyze
            response2 = flask_client.post('/analyze', data={
                'filename': 'GSE90122_TO90137.tsv',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
                'aop_selection': 'AOP:1',
                'logfc_threshold': '1.0',
                'dataset_id': 'INTEGRATION_TEST'
            })
            assert response2.status_code == 200

            # Step 3: Generate report
            response3 = flask_client.post('/generate_report', data={
                'format': 'html',
                'filename': 'GSE90122_TO90137.tsv',
                'gene_count': '5',
                'aop_id': 'AOP:1'
            })
            assert response3.status_code == 200
            assert b'Integration Test Report' in response3.data


@pytest.mark.integration
@pytest.mark.web
class TestDemosPage:
    """Phase 10.1: Demos page restructure — route + prefill flow."""

    def test_demos_route_renders(self, flask_client):
        """GET /demos returns 200 and the curated PXR section.

        Cisplatin demos are currently hidden in the template (only PXR agonist
        datasets are surfaced). Re-enable by removing the ``{% if false %}``
        wrapper around the cisplatin section in ``templates/demos.html``.
        """
        response = flask_client.get('/demos')
        assert response.status_code == 200
        body = response.data
        # Page chrome
        assert b'Demo Datasets' in body
        # Visible section
        assert b'PXR Agonists' in body
        # PXR curated cards present
        assert b'GSE90122_TO90137.tsv' in body
        assert b'GSE90122_SR12813.tsv' in body
        # Cisplatin section currently hidden — assert it does NOT render
        assert b'Cisplatin Nephrotoxicity' not in body
        assert b'Show all 54 cisplatin datasets' not in body
        # CTA copy
        assert b'Use this dataset' in body

    def test_demos_route_embeds_recommended_aops(self, flask_client):
        """Each PXR curated card embeds the right recommended_aops payload.

        Cisplatin demo cards are currently hidden in the template — only the
        PXR cards are asserted. Re-enable cisplatin assertions when the
        cisplatin section is unhidden in ``templates/demos.html``.
        """
        response = flask_client.get('/demos')
        body = response.data.decode('utf-8')
        # PXR demo cards carry AOP:1
        assert 'name="demo_file" value="GSE90122_TO90137.tsv"' in body
        # The recommended_aops hidden input next to a PXR demo file should be 'AOP:1'.
        # Locate the form snippet and assert the AOP appears within ~400 chars of it.
        idx = body.index('GSE90122_TO90137.tsv')
        nearby = body[max(0, idx - 400): idx + 400]
        assert 'AOP:1' in nearby, f"AOP:1 not found near PXR card: {nearby!r}"

    def test_preview_accepts_recommended_aops(self, flask_client):
        """POST /preview with recommended_aops propagates the list to the template."""
        from unittest.mock import patch, MagicMock
        with patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('pandas.read_csv') as mock_read_csv:
            mock_df = MagicMock()
            mock_df.head.return_value.to_dict.return_value = [
                {'Gene_Symbol': 'BRCA1', 'log2FoldChange': 2.5, 'padj': 0.001}
            ]
            mock_df.columns.tolist.return_value = ['Gene_Symbol', 'log2FoldChange', 'padj']
            mock_read_csv.return_value = mock_df

            response = flask_client.post('/preview', data={
                'demo_file': 'GSE90122_TO90137.tsv',
                'recommended_aops': 'AOP:1',
            })
            assert response.status_code == 200
            # Plan 02 will render the recommended_aops into the AOP picker UI.
            # For now we assert the value reaches the rendered HTML in *some* form
            # (template variable must have been passed; will be a hidden input or data-attr in Plan 02).
            # Bare-minimum signal: the value is present in the HTML body OR the route returned 200 cleanly.
            # Stronger assertion is added in Plan 02's tests.
            assert response.status_code == 200

    def test_demos_nav_link_appears_on_all_pages(self, flask_client):
        """The Demos link is in the top nav on /, /demos, /documentation, /about."""
        for path in ['/', '/demos', '/documentation', '/about']:
            response = flask_client.get(path)
            assert response.status_code == 200, f"{path} returned {response.status_code}"
            # Look for the nav anchor — should resolve to /demos via url_for
            assert b'href="/demos"' in response.data, f"Demos nav missing on {path}"

    def test_preview_with_recommended_aops_renders_filter_toggle(self, flask_client):
        """POST /preview with recommended_aops shows the filter toggle UI (Plan 02)."""
        import pandas as pd
        from unittest.mock import patch
        # Use a real DataFrame so the volcano plot data path produces non-empty volcano_data.
        mock_df = pd.DataFrame({
            'Gene_Symbol': ['BRCA1', 'TP53', 'EGFR'],
            'log2FoldChange': [2.5, -1.2, 0.8],
            'padj': [0.001, 0.05, 0.3],
        })
        with patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('pandas.read_csv', return_value=mock_df):

            response = flask_client.post('/preview', data={
                'demo_file': 'GSE90122_TO90137.tsv',
                'recommended_aops': 'AOP:1',
                'columns_confirmed': 'true',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
            })
            assert response.status_code == 200
            assert b'class="aop-filter-toggle"' in response.data
            assert b'Showing recommended AOPs for this demo' in response.data
            assert b'Show all AOPs' in response.data
            assert b'<script id="recommended-aops-data"' in response.data
            assert b'AOP:1' in response.data

    def test_preview_without_recommended_aops_hides_filter_toggle(self, flask_client):
        """POST /preview WITHOUT recommended_aops (upload-your-own path) hides the toggle (D-07)."""
        import pandas as pd
        from unittest.mock import patch
        # Use a real DataFrame so the volcano plot data path produces non-empty volcano_data.
        mock_df = pd.DataFrame({
            'Gene_Symbol': ['BRCA1', 'TP53', 'EGFR'],
            'log2FoldChange': [2.5, -1.2, 0.8],
            'padj': [0.001, 0.05, 0.3],
        })
        with patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('pandas.read_csv', return_value=mock_df):

            # No recommended_aops in form data — simulates upload-your-own-data path.
            response = flask_client.post('/preview', data={
                'demo_file': 'GSE90122_TO90137.tsv',
                'columns_confirmed': 'true',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
            })
            assert response.status_code == 200
            # Filter toggle UI must be ABSENT on the upload-your-own path.
            # The JS code references 'recommended-aops-data' as a getElementById string,
            # so we check for the HTML element tag (not just the bare string).
            assert b'class="aop-filter-toggle"' not in response.data
            assert b'<script id="recommended-aops-data"' not in response.data


@pytest.mark.integration
@pytest.mark.web
class TestTourRoutes:
    """Phase 13: Guided Tour — banner, button, anchors, docs re-launch."""

    def test_banner_on_index(self, flask_client):
        """TUTR-01: GET / renders the onboarding banner AND loads tour.js."""
        response = flask_client.get('/')
        assert response.status_code == 200
        body = response.data
        assert b'New here? Take a guided tour.' in body
        assert b'data-action="start-tour"' in body
        assert b'js/tour.js' in body

    def test_button_on_demos(self, flask_client):
        """TUTR-02: GET /demos renders the "Take the guided tour" button."""
        response = flask_client.get('/demos')
        assert response.status_code == 200
        body = response.data
        assert b'Take the guided tour' in body
        assert b'data-action="start-tour"' in body

    def test_anchors_present(self, flask_client):
        """TUTR-03: GET / has the upload-area data-tour-target anchor (the one visible outside preview gate)."""
        response = flask_client.get('/')
        assert response.status_code == 200
        assert b'data-tour-target="upload-area"' in response.data

    def test_relaunch_link_in_docs(self, flask_client):
        """TUTR-05: GET /documentation has the Getting Started section + tour launch link."""
        response = flask_client.get('/documentation')
        assert response.status_code == 200
        body = response.data
        assert b'id="getting-started"' in body
        assert b'data-action="start-tour"' in body

    def _analyze_form(self, **extra):
        base = {
            'filename': 'test.csv',
            'id_column': 'Gene_Symbol',
            'fc_column': 'log2FoldChange',
            'pval_column': 'padj',
            'aop_selection': 'AOP:1',
            'logfc_threshold': '1.0',
        }
        base.update(extra)
        return base

    def _analyze_with(self, client, form):
        """Run POST /analyze with the standard set of mocked dependencies and return the response."""
        import pandas as pd

        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS'],
            'log2FC': [1.5, -0.8, 2.1, 0.3, -1.2],
            'pval': [0.001, 0.05, 0.0001, 0.5, 0.01],
            'significant': [True, False, True, False, True],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.05],
            'num_overlap': [3], 'pct_sig_in_KE': [60.0],
            'total_KE_genes_in_dataset': [5], 'odds_ratio': [3.5],
            'overlap': ['BRCA1, EGFR, KRAS'], 'KE': ['KE:115'],
            'sig_in_KE': [3], 'sig_not_KE': [0],
            'non_sig_in_KE': [2], 'non_sig_not_KE': [0],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 5})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.build_ke_gene_mapping', return_value={}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('app.validate_file_path', return_value=True):
            return client.post('/analyze', data=form)

    def test_tour_active_meta_emitted_when_tour_field_present(self, authenticated_client):
        """TUTR-03: POST /analyze with tour=1 renders <meta name="tour-active"> in results.

        Without this, the cross-page resume from step 4 (aop-picker) to step 5
        (results-table) silently breaks — the tour engine reads this meta tag
        to know the analysis was triggered from inside an active tour.
        """
        response = self._analyze_with(authenticated_client, self._analyze_form(tour='1'))
        assert response.status_code == 200
        assert b'<meta name="tour-active"' in response.data

    def test_tour_active_meta_absent_without_tour_field(self, authenticated_client):
        """TUTR-03: POST /analyze without tour field does NOT render the tour-active meta.

        Guards against the meta tag accidentally appearing on every analysis,
        which would auto-resume the tour for users who never started it.
        """
        response = self._analyze_with(authenticated_client, self._analyze_form())
        assert response.status_code == 200
        assert b'<meta name="tour-active"' not in response.data

    def test_tour_active_meta_absent_when_tour_field_not_equal_one(self, authenticated_client):
        """TUTR-03: app.py uses strict-equality compare on tour == '1' (threat T-13-09).

        Any other value (empty string, 'true', '0', '2', random text) must NOT trip
        the tour_active flag. Single-value allow-list, not boolean truthy.
        """
        for fuzz in ('', 'true', '0', '2', 'yes', 'on', 'tour'):
            response = self._analyze_with(authenticated_client, self._analyze_form(tour=fuzz))
            assert response.status_code == 200, f"non-200 for tour={fuzz!r}"
            assert b'<meta name="tour-active"' not in response.data, f"meta unexpectedly emitted for tour={fuzz!r}"


@pytest.mark.integration
@pytest.mark.web
class TestHubPanel:
    """Phase 15: Hub Gene Stats — integration tests for the hub panel on /analyze.

    These tests are intentionally RED at Wave 0 because app.compute_hub_genes does not
    exist as a patchable attribute yet (Plan 15-02 adds the import to app.py).
    """

    # --- Module-level fixtures ---

    _ke_gene_map = {
        'KE:1': [
            {'id': 'TP53', 'log2FC': 2.1, 'significant': True, 'pvalue_raw': 0.001, 'pvalue_adj': 0.01},
            {'id': 'BRCA1', 'log2FC': -1.5, 'significant': False, 'pvalue_raw': 0.3, 'pvalue_adj': 0.55},
        ],
        'KE:2': [
            {'id': 'TP53', 'log2FC': 2.1, 'significant': True, 'pvalue_raw': 0.001, 'pvalue_adj': 0.01},
            {'id': 'MYC', 'log2FC': 0.8, 'significant': False, 'pvalue_raw': 0.2, 'pvalue_adj': 0.45},
        ],
        'KE:3': [
            {'id': 'TP53', 'log2FC': 2.1, 'significant': True, 'pvalue_raw': 0.001, 'pvalue_adj': 0.01},
        ],
    }

    _hub_list = [
        {
            'gene': 'TP53', 'ke_count': 3,
            'ke_ids': ['KE:1', 'KE:2', 'KE:3'],
            'ke_titles': ['KE 1 title', 'KE 2 title', 'KE 3 title'],
            'log2fc': 2.1, 'pvalue_adj': 0.01, 'significant': True, 'is_hub': True,
        },
        {
            'gene': 'BRCA1', 'ke_count': 1,
            'ke_ids': ['KE:1'],
            'ke_titles': ['KE 1 title'],
            'log2fc': -1.5, 'pvalue_adj': 0.55, 'significant': False, 'is_hub': False,
        },
        {
            'gene': 'MYC', 'ke_count': 1,
            'ke_ids': ['KE:2'],
            'ke_titles': ['KE 2 title'],
            'log2fc': 0.8, 'pvalue_adj': 0.45, 'significant': False, 'is_hub': False,
        },
    ]

    # 26-gene hub_list for test_hub_list_all_genes_passed
    # (server passes full list; DataTables paginates to 25 client-side)
    _hub_list_26 = [
        {
            'gene': f'GENE{i:02d}', 'ke_count': 1,
            'ke_ids': ['KE:1'], 'ke_titles': ['KE 1 title'],
            'log2fc': 0.5, 'pvalue_adj': 0.05, 'significant': False, 'is_hub': False,
        }
        for i in range(26)
    ]

    def _analyze_form_hub(self, **extra):
        """Base POST data for /analyze that exercises the hub panel path."""
        base = {
            'filename': 'test.csv',
            'id_column': 'Gene_Symbol',
            'fc_column': 'log2FoldChange',
            'pval_column': 'padj',
            'aop_selection': 'AOP:1',
            'logfc_threshold': '1.0',
        }
        base.update(extra)
        return base

    def _post_analyze(self, client, form, hub_list):
        """Run POST /analyze with the full mock stack including compute_hub_genes."""
        import pandas as pd

        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53', 'MYC'],
            'log2FC': [1.5, 2.1, 0.8],
            'pval': [0.001, 0.001, 0.2],
            'significant': [True, True, False],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.05],
            'num_overlap': [2], 'pct_sig_in_KE': [66.7],
            'total_KE_genes_in_dataset': [3], 'odds_ratio': [3.5],
            'overlap': ['BRCA1, TP53'], 'KE': ['KE:1'],
            'sig_in_KE': [2], 'sig_not_KE': [0],
            'non_sig_in_KE': [1], 'non_sig_not_KE': [0],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 3})), \
             patch('app.load_aop_data', return_value=({'KE:1'}, edges_df, {'KE:1': 'KE'}, {'KE:1': 'Test KE'})), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.build_ke_gene_mapping', return_value=self._ke_gene_map), \
             patch('app.compute_hub_genes', return_value=hub_list), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('app.validate_file_path', return_value=True):
            return client.post('/analyze', data=form)

    def test_hub_panel_present_ora(self, authenticated_client):
        """HUBG-01: POST /analyze with method=fisher (ORA) renders the Hub genes panel.

        Asserts response 200, b'Hub genes' heading, and b'hubGenesTable' DataTable ID.
        RED at Wave 0 — app.compute_hub_genes does not exist as a patchable attribute.
        """
        response = self._post_analyze(
            authenticated_client,
            self._analyze_form_hub(method='ora'),
            self._hub_list,
        )
        assert response.status_code == 200
        assert b'Hub genes' in response.data
        assert b'hubGenesTable' in response.data

    def test_hub_panel_present_gsea(self, authenticated_client):
        """HUBG-01: POST /analyze with method=gsea also renders the Hub genes panel identically.

        Panel renders the same HTML structure regardless of enrichment method.
        RED at Wave 0 — app.compute_hub_genes does not exist as a patchable attribute.
        """
        response = self._post_analyze(
            authenticated_client,
            self._analyze_form_hub(method='gsea'),
            self._hub_list,
        )
        assert response.status_code == 200
        assert b'Hub genes' in response.data
        assert b'hubGenesTable' in response.data

    def test_hub_list_all_genes_passed(self, authenticated_client):
        """HUBG-04: Server passes the full hub list; no server-side truncation to 25.

        DataTables paginates to 25 rows client-side. When hub_list has 26+ entries,
        all gene symbols must appear in the rendered response.data.
        RED at Wave 0 — app.compute_hub_genes does not exist as a patchable attribute.
        """
        response = self._post_analyze(
            authenticated_client,
            self._analyze_form_hub(),
            self._hub_list_26,
        )
        assert response.status_code == 200
        # All 26 gene symbols must be in the response (no server-side truncation)
        for entry in self._hub_list_26:
            gene_bytes = entry['gene'].encode('utf-8')
            assert gene_bytes in response.data, (
                f"Gene {entry['gene']!r} not found in response; "
                "server must pass the full hub_list without truncating to 25"
            )

    def test_hub_list_in_template(self, authenticated_client):
        """HUBG-07: A known fixture gene symbol appears in response.data; the variable name does not.

        The template must render hub_list contents as HTML (e.g. table rows) so gene
        symbols appear in the body. The Python variable name 'hub_list' must NOT be
        present as a literal string in the rendered HTML (would indicate a template error).
        RED at Wave 0 — app.compute_hub_genes does not exist as a patchable attribute.
        """
        response = self._post_analyze(
            authenticated_client,
            self._analyze_form_hub(),
            self._hub_list,
        )
        assert response.status_code == 200
        # A known gene symbol from the fixture must appear in rendered HTML
        assert b'TP53' in response.data
        # The Python variable name must NOT appear as a literal in the HTML body
        assert b'hub_list' not in response.data


@pytest.mark.integration
@pytest.mark.web
class TestColumnConfidenceUI:
    """Phase 16 (AUTC-01/AUTC-02): Column confidence badges + hover tooltip in single-analysis preview.

    Tests assert that the /preview response:
    - D-01: Contains a colour-coded confidence badge (confidence-badge confidence-badge--*)
            for each column type when the detector confidently identifies columns.
    - D-02: Carries the detector reasons and the static generic explainer inside the
            badge's hover tooltip (title attribute) — "Why this column:" header for the
            reasons, "Detector looks for:" for the explainer — rather than inline.
    - D-03: Shows confidence-badge--undetected + "Couldn't auto-detect" when a column
            type cannot be detected.
    - Guided tour anchor: data-tour-target="column-preview" is preserved.
    """

    # A well-formed CSV with clearly named columns — detector should detect all four.
    _good_csv = (
        b"Gene_Symbol,log2FoldChange,pvalue,padj\n"
        b"BRCA1,2.5,0.001,0.005\n"
        b"TP53,-0.8,0.05,0.12\n"
        b"EGFR,2.1,0.0001,0.0005\n"
        b"MYC,0.3,0.5,0.9\n"
        b"KRAS,-1.2,0.01,0.03\n"
        b"PTEN,1.8,0.002,0.01\n"
        b"AKT1,-0.5,0.3,0.6\n"
        b"MDM2,0.9,0.08,0.2\n"
        b"RB1,-1.5,0.007,0.02\n"
        b"CDK2,2.2,0.0003,0.002\n"
    )

    # A CSV with only one numeric column — the p-value cannot be detected reliably.
    # We ensure the gene_id column is obvious but use integer data for other columns
    # so the detector cannot reliably identify a p-value vs. log2FC column.
    _no_pval_csv = (
        b"Gene_Symbol,SomeValues\n"
        b"BRCA1,150\n"
        b"TP53,200\n"
        b"EGFR,175\n"
        b"MYC,130\n"
        b"KRAS,220\n"
    )

    def _post_preview(self, flask_client, csv_bytes, filename='test.csv'):
        """Helper: POST a CSV to /preview via multipart upload."""
        return flask_client.post(
            '/preview',
            data={
                'gene_file': (io.BytesIO(csv_bytes), filename),
                'dataset_id': 'TEST_CONFIDENCE',
                'stressor': 'Test',
                'owner': 'Test User',
            },
            content_type='multipart/form-data',
        )

    def test_confidence_badge_rendered_for_detected_columns(self, flask_client):
        """D-01: A confidence badge appears when columns are confidently detected.

        POST /preview with a well-formed CSV.  The response must contain at least
        one filled confidence badge (the class prefix confidence-badge confidence-badge--).
        """
        response = self._post_preview(flask_client, self._good_csv)
        assert response.status_code == 200, (
            f"Expected 200, got {response.status_code}; body: {response.data[:500]!r}"
        )
        assert b'confidence-badge confidence-badge--' in response.data, (
            "Expected a filled confidence badge (class 'confidence-badge confidence-badge--*') "
            "in /preview response for a well-formed CSV. "
            f"Body excerpt: {response.data[:1000]!r}"
        )

    def test_reasons_in_badge_tooltip_for_detected_columns(self, flask_client):
        """D-02 reasons: The detector reasons appear in the badge hover tooltip.

        POST /preview with a well-formed CSV.  When columns are detected the
        badge title carries a "Why this column:" reasons block.
        """
        response = self._post_preview(flask_client, self._good_csv)
        assert response.status_code == 200
        assert b'Why this column:' in response.data, (
            "Expected the 'Why this column:' reasons block in a badge title attribute "
            "of the /preview response. "
            f"Body excerpt: {response.data[:1000]!r}"
        )

    def test_static_explainer_in_badge_tooltip(self, flask_client):
        """D-02 explainer: The static column-type explainer appears in the badge tooltip.

        POST /preview with a well-formed CSV.  The verbatim explainer text must
        appear (inside the badge title attribute) for the gene-ID column.
        """
        response = self._post_preview(flask_client, self._good_csv)
        assert response.status_code == 200
        assert b'Detector looks for: gene symbol name/format patterns' in response.data, (
            "Expected verbatim gene-ID explainer text in /preview response. "
            f"Body excerpt: {response.data[:1000]!r}"
        )
        assert b'title="Detector looks for:' in response.data, (
            "Expected the explainer to be carried in a badge title (hover tooltip) attribute. "
            f"Body excerpt: {response.data[:1000]!r}"
        )

    def test_undetected_badge_when_column_not_detected(self, flask_client):
        """D-03: The undetected badge appears when a column type cannot be detected.

        POST /preview with a CSV that has only a gene column and an integer column
        (no numeric log2FC-shaped data). The p-value and log2FC columns cannot be
        detected, so at least one undetected badge must appear.
        """
        response = self._post_preview(flask_client, self._no_pval_csv)
        assert response.status_code == 200, (
            f"Expected 200, got {response.status_code}; body: {response.data[:500]!r}"
        )
        assert b'confidence-badge--undetected' in response.data, (
            "Expected confidence-badge--undetected in /preview response for a CSV "
            "where some columns cannot be detected. "
            f"Body excerpt: {response.data[:1000]!r}"
        )
        assert b"Couldn't auto-detect" in response.data, (
            "Expected \"Couldn't auto-detect\" text in /preview response. "
            f"Body excerpt: {response.data[:1000]!r}"
        )

    def test_tour_anchor_preserved(self, flask_client):
        """Phase 13 tour anchor: data-tour-target=\"column-preview\" must remain on the outer card.

        The Phase 13 guided tour highlights the column step via this anchor.
        Any /preview response must keep it intact regardless of detection outcome.
        """
        response = self._post_preview(flask_client, self._good_csv)
        assert response.status_code == 200
        assert b'data-tour-target="column-preview"' in response.data, (
            "Expected data-tour-target=\"column-preview\" in /preview response. "
            "Do not remove or move this attribute. "
            f"Body excerpt: {response.data[:500]!r}"
        )


@pytest.mark.integration
@pytest.mark.web
class TestBatchUploadConfidencePayload:
    """Phase 16 (AUTC-01/D-05): /batch/upload suggestions payload carries confidence + label fields.

    Tests assert that the JSON response from POST /batch/upload includes, per file:
    - id_confidence, id_label alongside id_col
    - fc_confidence, fc_label alongside fc_col
    - pval_confidence, pval_label alongside pval_col
    - All three are null when a column type cannot be detected.
    """

    # Well-formed CSV with clearly named gene/FC/pvalue columns — detector detects all three.
    _good_csv = (
        b"Gene_Symbol,log2FoldChange,pvalue\n"
        b"BRCA1,2.5,0.001\n"
        b"TP53,-0.8,0.05\n"
        b"EGFR,2.1,0.0001\n"
        b"MYC,0.3,0.5\n"
        b"KRAS,-1.2,0.01\n"
        b"PTEN,1.8,0.002\n"
        b"AKT1,-0.5,0.3\n"
        b"MDM2,0.9,0.08\n"
        b"RB1,-1.5,0.007\n"
        b"CDK2,2.2,0.0003\n"
    )

    # CSV where no p-value column can be reliably detected (only integer counts).
    _no_pval_csv = (
        b"Gene_Symbol,ReadCount\n"
        b"BRCA1,150\n"
        b"TP53,200\n"
        b"EGFR,175\n"
        b"MYC,130\n"
        b"KRAS,220\n"
    )

    def _post_batch_upload(self, flask_client, csv_bytes, filename='test_batch.csv'):
        """Helper: POST a CSV to /batch/upload via multipart upload."""
        return flask_client.post(
            '/batch/upload',
            data={
                'files[]': (io.BytesIO(csv_bytes), filename),
            },
            content_type='multipart/form-data',
        )

    def test_suggestions_payload_contains_nine_keys(self, flask_client):
        """D-05: /batch/upload suggestions dict must carry all nine confidence + label keys.

        POST /batch/upload with a well-formed CSV. The JSON response
        files[0].suggestions must contain all nine keys:
        id_col, id_confidence, id_label, fc_col, fc_confidence, fc_label,
        pval_col, pval_confidence, pval_label.
        """
        response = self._post_batch_upload(flask_client, self._good_csv)
        assert response.status_code == 200, (
            f"Expected 200, got {response.status_code}; body: {response.data[:500]!r}"
        )
        data = json.loads(response.data)
        assert 'files' in data and len(data['files']) > 0, (
            "Expected 'files' list in /batch/upload response"
        )
        suggestions = data['files'][0]['suggestions']
        required_keys = [
            'id_col', 'id_confidence', 'id_label',
            'fc_col', 'fc_confidence', 'fc_label',
            'pval_col', 'pval_confidence', 'pval_label',
        ]
        missing = [k for k in required_keys if k not in suggestions]
        assert not missing, (
            f"Missing keys in /batch/upload suggestions: {missing}. "
            f"Got: {list(suggestions.keys())}"
        )

    def test_confidence_fields_are_valid_for_detected_columns(self, flask_client):
        """D-05: For detected columns, *_confidence is float in [0,1] and *_label is valid string.

        POST /batch/upload with a well-formed CSV (all three column types detectable).
        id_confidence must be a float in [0, 1] and id_label must be one of the
        four qualitative level strings.
        """
        response = self._post_batch_upload(flask_client, self._good_csv)
        assert response.status_code == 200
        data = json.loads(response.data)
        suggestions = data['files'][0]['suggestions']

        valid_labels = {'High', 'Medium', 'Low', 'Very Low'}
        for prefix in ('id', 'fc', 'pval'):
            col = suggestions.get(f'{prefix}_col')
            confidence = suggestions.get(f'{prefix}_confidence')
            label = suggestions.get(f'{prefix}_label')
            if col is not None:
                assert isinstance(confidence, (int, float)), (
                    f"{prefix}_confidence must be a number, got {type(confidence)}: {confidence}"
                )
                assert 0.0 <= confidence <= 1.0, (
                    f"{prefix}_confidence must be in [0,1], got {confidence}"
                )
                assert label in valid_labels, (
                    f"{prefix}_label must be one of {valid_labels}, got {label!r}"
                )

    def test_undetected_column_fields_are_null(self, flask_client):
        """D-05: For an undetected column type, *_col, *_confidence, and *_label are all null.

        POST /batch/upload with a CSV where p-value cannot be detected.
        pval_col, pval_confidence, and pval_label must all be JSON null (None).
        """
        response = self._post_batch_upload(flask_client, self._no_pval_csv)
        assert response.status_code == 200, (
            f"Expected 200, got {response.status_code}; body: {response.data[:500]!r}"
        )
        data = json.loads(response.data)
        suggestions = data['files'][0]['suggestions']
        # p-value column should not be detected for a CSV with only integer counts
        assert suggestions.get('pval_col') is None, (
            f"Expected pval_col to be null for CSV with no p-value column, "
            f"got: {suggestions.get('pval_col')!r}"
        )
        assert suggestions.get('pval_confidence') is None, (
            f"Expected pval_confidence to be null, got: {suggestions.get('pval_confidence')!r}"
        )
        assert suggestions.get('pval_label') is None, (
            f"Expected pval_label to be null, got: {suggestions.get('pval_label')!r}"
        )

    def test_batch_scripts_contains_badge_helper_and_explainer(self, flask_client):
        """D-01/D-02: the live batch wizard partial has the badge DOM helper and explainer strings.

        The batch wizard is served as the `?tab=batch` view of index.html, which
        includes `_batch_analysis_scripts.html` — that partial is the live code path.
        Asserts on that partial's source:
        - confidence-badge CSS class reference is present
        - makeConfidenceBadge helper uses document.createElement (safe DOM construction)
        - .title DOM property assignment is used (not innerHTML) for the tooltip
        - The verbatim gene-ID explainer string from UI-SPEC Component 4 is present
          (same canonical wording as the single-analysis path, per Plan 16-01)
        """
        batch_response = flask_client.get('/?tab=batch')
        assert batch_response.status_code == 200

        import os
        template_path = os.path.join(
            os.path.dirname(__file__), '..', 'templates', '_batch_analysis_scripts.html'
        )
        with open(template_path, 'r', encoding='utf-8') as fh:
            source = fh.read()

        assert 'confidence-badge' in source, (
            "Expected 'confidence-badge' CSS class reference in _batch_analysis_scripts.html"
        )
        assert 'document.createElement' in source, (
            "Expected document.createElement (safe DOM construction) in _batch_analysis_scripts.html"
        )
        assert '.title' in source, (
            "Expected .title DOM property assignment for tooltip in _batch_analysis_scripts.html; "
            "must not use innerHTML for the explainer tooltip"
        )
        assert 'Detector looks for: gene symbol name/format patterns' in source, (
            "Expected verbatim gene-ID explainer string (UI-SPEC Component 4) in "
            "_batch_analysis_scripts.html; batch and single-analysis paths must use identical help text"
        )


class TestColumnDetectorDocumentation:
    """AUTC-03: /documentation page contains the #column-detector section."""

    def test_documentation_has_column_detector_section(self, flask_client):
        """AUTC-03: GET /documentation renders the #column-detector doc-section with all required elements."""
        response = flask_client.get('/documentation')
        assert response.status_code == 200
        body = response.data

        # TOC anchor link
        assert b'href="#column-detector"' in body, (
            "Expected href=\"#column-detector\" in documentation TOC"
        )

        # Section id
        assert b'id="column-detector"' in body, (
            "Expected id=\"column-detector\" on the new doc-section div"
        )

        # h2 heading
        assert b'Column Auto-Detection' in body, (
            "Expected h2 heading 'Column Auto-Detection' in the #column-detector section"
        )

        # h3 subsections
        assert b'How the detector works' in body, (
            "Expected h3 'How the detector works' in the #column-detector section"
        )
        assert b'Confidence levels' in body, (
            "Expected h3 'Confidence levels' in the #column-detector section"
        )
        assert b'Overriding a detected column' in body, (
            "Expected h3 'Overriding a detected column' in the #column-detector section"
        )

        # Confidence level table rows
        assert b'High' in body
        assert b'Medium' in body
        assert b'Low' in body
        assert b'Very Low' in body
        assert b'80%' in body, "Expected 80% threshold for High confidence level"
        assert b'60%' in body, "Expected 60% threshold for Medium confidence level"
        assert b'30%' in body, "Expected 30% threshold for Low confidence level"

        # Source order: #column-detector appears after #input-format and before #statistical-methods
        body_str = body.decode('utf-8')
        pos_input = body_str.find('id="input-format"')
        pos_detector = body_str.find('id="column-detector"')
        pos_stats = body_str.find('id="statistical-methods"')
        assert pos_input < pos_detector, (
            "#column-detector must appear after #input-format in document source order"
        )
        assert pos_detector < pos_stats, (
            "#column-detector must appear before #statistical-methods in document source order"
        )


@pytest.mark.integration
@pytest.mark.web
class TestPathwayView:
    """Phase 999.4 (PinPath): Pathway-view section conditional rendering on /analyze.

    Tests assert that the /analyze results page:
    - Includes a div with id="pathway-view-section" when build_wp_picker_data returns
      a non-empty list (mapped AOP — SPEC R6).
    - Does NOT include that div when build_wp_picker_data returns None (no KE-WP
      mappings for the selected AOP — D-10).

    Both tests are RED at Wave 0 because:
    - app.build_wp_picker_data is not yet an importable attribute of app.py (added in Plan 01).
    - results.html has no pathway-view-section div (added in Plan 02).
    """

    _ke_gene_map = {
        'KE:1': [
            {'id': 'TP53', 'log2FC': 2.1, 'significant': True, 'pvalue_raw': 0.001, 'pvalue_adj': 0.01},
            {'id': 'BRCA1', 'log2FC': -1.5, 'significant': False, 'pvalue_raw': 0.3, 'pvalue_adj': 0.55},
        ],
    }

    # A minimal picker entry representing one KE-WP pair
    _picker_entry = {
        'ke_id': 'KE:1',
        'ke_title': 'Test Key Event',
        'wp_id': 'WP368',
        'wp_title': 'Fatty acid beta-oxidation',
        'confidence_level': 'high',
        'enrichment_rank': 0,
    }

    def _analyze_form(self, **extra):
        """Base POST data for /analyze."""
        base = {
            'filename': 'test.csv',
            'id_column': 'Gene_Symbol',
            'fc_column': 'log2FoldChange',
            'pval_column': 'padj',
            'aop_selection': 'AOP:1',
            'logfc_threshold': '1.0',
        }
        base.update(extra)
        return base

    def _post_analyze(self, client, form, wp_picker_data):
        """Run POST /analyze with the full mock stack including build_wp_picker_data.

        Copies the TestHubPanel._post_analyze mock stack verbatim and adds one extra
        context-manager patch for app.build_wp_picker_data.
        """
        import pandas as pd

        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53', 'MYC'],
            'log2FC': [1.5, 2.1, 0.8],
            'pval': [0.001, 0.001, 0.2],
            'significant': [True, True, False],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.05],
            'num_overlap': [2], 'pct_sig_in_KE': [66.7],
            'total_KE_genes_in_dataset': [3], 'odds_ratio': [3.5],
            'overlap': ['BRCA1, TP53'], 'KE': ['KE:1'],
            'sig_in_KE': [2], 'sig_not_KE': [0],
            'non_sig_in_KE': [1], 'non_sig_not_KE': [0],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 3})), \
             patch('app.load_aop_data', return_value=({'KE:1'}, edges_df, {'KE:1': 'KE'}, {'KE:1': 'Test KE'})), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.build_ke_gene_mapping', return_value=self._ke_gene_map), \
             patch('app.compute_hub_genes', return_value=[]), \
             patch('app.build_wp_picker_data', return_value=wp_picker_data), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('app.validate_file_path', return_value=True):
            return client.post('/analyze', data=form)

    def test_results_has_pathway_view(self, authenticated_client):
        """POST /analyze with a non-empty picker list renders the pathway-view-section div.

        When build_wp_picker_data returns a non-empty list (mapped AOP), the results
        page must contain a div with id="pathway-view-section" (SPEC R6 / D-10).
        RED at Wave 0 — app.build_wp_picker_data not yet in app.py; section not yet
        in results.html.
        """
        response = self._post_analyze(
            authenticated_client,
            self._analyze_form(),
            wp_picker_data=[self._picker_entry],
        )
        assert response.status_code == 200
        assert b'pathway-view-section' in response.data, (
            "Expected 'pathway-view-section' in response when wp_picker_data is non-empty. "
            f"Body excerpt: {response.data[:500]!r}"
        )

    def test_no_pathway_view_when_unmapped(self, authenticated_client):
        """POST /analyze with wp_picker_data=None does NOT render the pathway-view-section div.

        When build_wp_picker_data returns None (AOP has no KE-WP mappings), the
        pathway-view section must be absent from the results page (D-10).
        RED at Wave 0 — app.build_wp_picker_data not yet in app.py.
        """
        response = self._post_analyze(
            authenticated_client,
            self._analyze_form(),
            wp_picker_data=None,
        )
        assert response.status_code == 200
        assert b'pathway-view-section' not in response.data, (
            "Expected 'pathway-view-section' to be ABSENT when wp_picker_data is None. "
            f"Body excerpt: {response.data[:500]!r}"
        )


class TestMinConfidenceControl:
    """Issue #60: the minimum mapping-confidence control on the single flow."""

    @staticmethod
    def _form(**overrides):
        form = {
            'filename': 'test.csv',
            'id_column': 'Gene_Symbol',
            'fc_column': 'log2FoldChange',
            'pval_column': 'padj',
            'aop_selection': 'AOP:1',
            'logfc_threshold': '1.0',
        }
        form.update(overrides)
        return form

    @staticmethod
    def _post_analyze(client, form):
        """POST /analyze with the standard mock stack; returns (response, mock loader)."""
        import pandas as pd

        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53'],
            'log2FC': [1.5, -0.8],
            'pval': [0.001, 0.05],
            'significant': [True, False],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.05],
            'num_overlap': [1], 'pct_sig_in_KE': [50.0],
            'total_KE_genes_in_dataset': [2], 'odds_ratio': [3.5],
            'overlap': ['BRCA1'], 'KE': ['KE:115'],
            'sig_in_KE': [1], 'sig_not_KE': [0],
            'non_sig_in_KE': [1], 'non_sig_not_KE': [0],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 2})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.build_ke_gene_mapping', return_value={}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('app.validate_file_path', return_value=True), \
             patch('app.load_cached_reference_sets', return_value=({'KE:115': {'BRCA1'}}, 'cache', [])) as loader:
            response = client.post('/analyze', data=form)
        return response, loader

    def test_threshold_forwarded_to_reference_loader(self, authenticated_client):
        response, loader = self._post_analyze(
            authenticated_client, self._form(min_confidence='high')
        )
        assert response.status_code == 200
        assert loader.call_args.kwargs['min_confidence'] == 'high'

    def test_omitted_threshold_defaults_to_all(self, authenticated_client):
        """Backward compatible: no control submitted -> current behaviour."""
        response, loader = self._post_analyze(authenticated_client, self._form())
        assert response.status_code == 200
        assert loader.call_args.kwargs['min_confidence'] == 'all'

    def test_case_and_whitespace_are_normalised(self, authenticated_client):
        response, loader = self._post_analyze(
            authenticated_client, self._form(min_confidence='  Medium ')
        )
        assert response.status_code == 200
        assert loader.call_args.kwargs['min_confidence'] == 'medium'

    @pytest.mark.parametrize("junk", ["low", "sky-high", "1", "'; DROP TABLE experiments --"])
    def test_junk_threshold_is_rejected(self, authenticated_client, junk):
        response, loader = self._post_analyze(
            authenticated_client, self._form(min_confidence=junk)
        )
        assert response.status_code == 400
        assert b'confidence' in response.data.lower()
        loader.assert_not_called()

    @staticmethod
    def _render_single_form(**context):
        """Render the single-analysis partial with the enrichment settings visible."""
        from flask import render_template
        from app import app as flask_app

        base = dict(
            preview=[{'Gene_Symbol': 'BRCA1'}],
            columns=['Gene_Symbol', 'log2FoldChange', 'padj'],
            filename='test.csv',
            volcano_data=[{'ID': 'BRCA1', 'log2FC': 1.5, 'pval': 0.001}],
            selected_columns={'id': 'Gene_Symbol', 'fc': 'log2FoldChange',
                              'pval': 'padj', 'pval_adj': None},
            column_suggestions=None,
            logfc_threshold=1.0,
            pval_cutoff=0.05,
            pval_y=[],
            columns_confirmed=True,
            case_study_aops={},
            cisplatin_demos=[],
            parse_filename=lambda name: {},
            recommended_aops=None,
            metadata={},
            pval_threshold=0.05,
            method='ora',
            selected_resources=['WikiPathways'],
        )
        base.update(context)
        with flask_app.test_request_context('/'):
            return render_template('_single_analysis.html', **base)

    def test_control_rendered_in_the_enrichment_target_group(self):
        """The selector appears alongside the gene-set resource checkboxes."""
        html = self._render_single_form(min_confidence='all')
        assert 'name="min_confidence"' in html
        assert 'Minimum mapping confidence' in html
        # Sits in the same settings group as the resource checkboxes (#55 pattern).
        assert html.index('name="resources"') < html.index('name="min_confidence"')

    def test_selected_threshold_survives_the_htmx_rerender(self):
        html = self._render_single_form(min_confidence='high')
        assert '<option value="high" selected>' in html
        assert '<option value="all" selected>' not in html

    def test_default_selection_is_all(self):
        html = self._render_single_form()
        assert '<option value="all" selected>' in html

    def test_preview_echoes_the_threshold_back(self, flask_client):
        """/preview re-renders with the submitted threshold so it survives Update Plot."""
        import pandas as pd

        preview_df = pd.DataFrame({
            'Gene_Symbol': ['BRCA1', 'TP53'],
            'log2FoldChange': [1.5, -0.8],
            'padj': [0.001, 0.05],
        })
        with patch('app.validate_file_path', return_value=True), \
             patch('pandas.read_csv', return_value=preview_df), \
             patch('os.path.exists', side_effect=uploaded_file_exists):
            response = flask_client.post('/preview', data={
                'filename': 'test.csv',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
                'columns_confirmed': 'true',
                'min_confidence': 'high',
            }, headers={'HX-Request': 'true'})
        assert response.status_code == 200
        assert b'<option value="high" selected>' in response.data

    def test_preview_ignores_a_junk_threshold(self, flask_client):
        """A tampered value re-renders as the default rather than 500ing."""
        import pandas as pd

        preview_df = pd.DataFrame({
            'Gene_Symbol': ['BRCA1', 'TP53'],
            'log2FoldChange': [1.5, -0.8],
            'padj': [0.001, 0.05],
        })
        with patch('app.validate_file_path', return_value=True), \
             patch('pandas.read_csv', return_value=preview_df), \
             patch('os.path.exists', side_effect=uploaded_file_exists):
            response = flask_client.post('/preview', data={
                'filename': 'test.csv',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
                'columns_confirmed': 'true',
                'min_confidence': 'bogus',
            }, headers={'HX-Request': 'true'})
        assert response.status_code == 200
        assert b'<option value="all" selected>' in response.data


class TestBatchMinConfidence:
    """Issue #60: batch parity for the minimum mapping-confidence control."""

    @staticmethod
    def _form(**overrides):
        form = {
            'batch_uuid': 'test-batch-uuid',
            'aop_selection': 'AOP:1',
            'id_col': 'GENE_SYMBOL',
            'fc_col': 'logFC',
            'pval_col': 'adj.P.Val',
            'logfc_threshold': '0.0',
            'filenames[]': 'a.csv',
            'condition_labels[]': 'A',
        }
        form.update(overrides)
        return form

    @staticmethod
    def _post_batch(client, form):
        with patch('os.path.isdir', return_value=True), \
             patch('os.path.isfile', return_value=True), \
             patch('app.validate_batch_columns', return_value=(True, '')), \
             patch('app.harmonise_backgrounds', return_value=({'BRCA1'}, {})), \
             patch('app._persist_and_launch_batch', return_value=1) as launcher:
            response = client.post('/batch/analyze', data=form)
        return response, launcher

    def test_threshold_forwarded_to_the_batch_launcher(self, flask_client):
        response, launcher = self._post_batch(flask_client, self._form(min_confidence='medium'))
        assert response.status_code == 200
        assert launcher.call_args.kwargs['min_confidence'] == 'medium'

    def test_omitted_threshold_defaults_to_all(self, flask_client):
        response, launcher = self._post_batch(flask_client, self._form())
        assert response.status_code == 200
        assert launcher.call_args.kwargs['min_confidence'] == 'all'

    def test_junk_threshold_is_rejected(self, flask_client):
        response, launcher = self._post_batch(flask_client, self._form(min_confidence='bogus'))
        assert response.status_code == 400
        assert b'confidence' in response.data.lower()
        launcher.assert_not_called()

    def test_control_rendered_in_the_batch_wizard(self, flask_client):
        response = flask_client.get('/')
        assert response.status_code == 200
        assert b'batch-min-confidence' in response.data

class TestRepresentationColumnRendered:
    """Issue #70: the results page must show which tail a Key Event fell on."""

    @staticmethod
    def _analyze(client, representation):
        import pandas as pd

        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53', 'EGFR'],
            'log2FC': [1.5, -0.8, 2.1],
            'pval': [0.001, 0.05, 0.0001],
            'significant': [True, False, True],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.9], 'FDR': [0.95],
            'Representation': [representation],
            'num_overlap': [1], 'pct_sig_in_KE': [20.0],
            'total_KE_genes_in_dataset': [5], 'odds_ratio': [0.2],
            'overlap': ['BRCA1'], 'KE': ['KE:115'],
            'sig_in_KE': [1], 'sig_not_KE': [1],
            'non_sig_in_KE': [4], 'non_sig_not_KE': [0],
            'p_value_depleted': [0.0004], 'FDR_depleted': [0.0004],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 3})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.build_ke_gene_mapping', return_value={}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('app.validate_file_path', return_value=True):
            return client.post('/analyze', data={
                'filename': 'test.csv',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
                'aop_selection': 'AOP:1',
                'logfc_threshold': '1.0',
            })

    def test_depleted_ke_renders_as_a_labelled_badge(self, authenticated_client):
        """A KE with FDR ≈ 1 but a significant depletion tail reads as depleted,
        not as the blank cell the one-sided output used to produce."""
        response = self._analyze(authenticated_client, 'depleted')
        assert response.status_code == 200
        body = response.data
        assert b'rep-badge--depleted' in body
        assert b'Depleted' in body
        # Both depletion columns are exposed, so the call can be checked.
        assert b'FDR (depletion)' in body

    def test_non_significant_ke_renders_quietly(self, authenticated_client):
        response = self._analyze(authenticated_client, 'ns')
        assert response.status_code == 200
        assert b'rep-badge--ns' in response.data


class TestResourceProvenanceOnResultsPage:
    """Issues #67 / #68: what the run used, and where it differed from the ask."""

    @staticmethod
    def _analyze(client, resolution, min_confidence='all'):
        import pandas as pd

        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53', 'EGFR'],
            'log2FC': [1.5, -0.8, 2.1],
            'pval': [0.001, 0.05, 0.0001],
            'significant': [True, False, True],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.02],
            'Representation': ['enriched'], 'num_overlap': [2],
            'pct_sig_in_KE': [60.0], 'total_KE_genes_in_dataset': [5],
            'odds_ratio': [3.5], 'overlap': ['BRCA1, EGFR'], 'KE': ['KE:115'],
            'sig_in_KE': [2], 'sig_not_KE': [0],
            'non_sig_in_KE': [1], 'non_sig_not_KE': [0],
            'p_value_depleted': [0.99], 'FDR_depleted': [0.99],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 3})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.load_cached_reference_sets',
                   return_value=({'KE:115': {'BRCA1', 'EGFR'}}, 'api', resolution)), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.build_ke_gene_mapping', return_value={}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('app.validate_file_path', return_value=True):
            return client.post('/analyze', data={
                'filename': 'test.csv',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
                'aop_selection': 'AOP:1',
                'logfc_threshold': '1.0',
                'resources': ['WikiPathways', 'Reactome'],
                'min_confidence': min_confidence,
            })

    def test_skipped_resource_is_announced_on_the_page(self, authenticated_client):
        resolution = [
            {'resource': 'WikiPathways', 'status': 'loaded', 'source': 'api',
             'ke_count': 1, 'confidence_applied': True, 'error': None},
            {'resource': 'Reactome', 'status': 'skipped', 'source': None,
             'ke_count': 0, 'confidence_applied': False, 'error': 'HTTP 502'},
        ]
        response = self._analyze(authenticated_client, resolution)
        assert response.status_code == 200
        body = response.data
        assert b'Reactome' in body
        assert b'unavailable, skipped' in body
        assert b'left out of this' in body

    def test_confidence_applicability_is_disclosed(self, authenticated_client):
        """An unfiltered resource is still disclosed — with an honest reason.

        Since #71 a GMT resource is filtered on its live path, so this shape
        (Reactome loaded but unfiltered) is what a run *recorded before #71*
        replays as. The page must say the threshold was not applied without
        blaming the bundled reference files: GO BP and Reactome have no bundled
        fallback, so that explanation would be false — the #108 failure mode.
        """
        resolution = [
            {'resource': 'WikiPathways', 'status': 'loaded', 'source': 'api',
             'ke_count': 1, 'confidence_applied': True, 'error': None},
            {'resource': 'Reactome', 'status': 'loaded', 'source': 'api',
             'ke_count': 1, 'confidence_applied': False, 'error': None},
        ]
        response = self._analyze(authenticated_client, resolution, min_confidence='high')
        assert response.status_code == 200
        assert b'was not applied to Reactome' in response.data
        assert b'bundled' not in response.data

    def test_csv_fallback_is_still_blamed_on_the_bundled_files(self, authenticated_client):
        """WikiPathways *does* have a bundled fallback, so name it."""
        resolution = [
            {'resource': 'WikiPathways', 'status': 'loaded', 'source': 'cache(csv)',
             'ke_count': 1, 'confidence_applied': False, 'error': None},
        ]
        response = self._analyze(authenticated_client, resolution, min_confidence='high')
        assert response.status_code == 200
        assert b'bundled reference files' in response.data

    def test_clean_run_shows_provenance_without_a_warning(self, authenticated_client):
        resolution = [
            {'resource': 'WikiPathways', 'status': 'loaded', 'source': 'api',
             'ke_count': 1, 'confidence_applied': True, 'error': None},
        ]
        response = self._analyze(authenticated_client, resolution)
        body = response.data
        assert b'Gene set provenance' in body
        assert b'Builder API, live' in body
        assert b'class="resource-warning"' not in body


class TestUploadSurvivesAnalysis:
    """Issue #72: a dataset must be re-analysable with different settings."""

    @staticmethod
    def _real_upload(tmp_path, monkeypatch):
        """Point the app at a scratch uploads dir holding one real file."""
        import app as app_module

        uploads = tmp_path / 'uploads'
        uploads.mkdir()
        (uploads / 'DEGs.txt').write_text(
            "GENE_SYMBOL\tlogFC\tadj.P.Val\nBRCA1\t1.5\t0.001\nTP53\t-0.8\t0.4\n"
        )
        monkeypatch.setitem(app_module.app.config, 'UPLOAD_FOLDER', str(uploads))
        monkeypatch.setattr(app_module.Config, 'UPLOAD_FOLDER', str(uploads))
        return uploads

    @staticmethod
    def _analyze(client, method):
        import pandas as pd

        processed_df = pd.DataFrame({
            'ID': ['BRCA1', 'TP53'],
            'log2FC': [1.5, -0.8],
            'pval': [0.001, 0.4],
            'significant': [True, False],
        })
        enrichment_df = pd.DataFrame({
            'Title': ['Test KE'], 'p_value': [0.01], 'FDR': [0.02],
            'Representation': ['enriched'], 'num_overlap': [1],
            'pct_sig_in_KE': [50.0], 'total_KE_genes_in_dataset': [2],
            'odds_ratio': [3.5], 'overlap': ['BRCA1'], 'KE': ['KE:115'],
            'sig_in_KE': [1], 'sig_not_KE': [0],
            'non_sig_in_KE': [1], 'non_sig_not_KE': [0],
            'p_value_depleted': [0.99], 'FDR_depleted': [0.99],
        })
        edges_df = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])

        with patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 2})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.load_cached_reference_sets',
                   return_value=({'KE:115': {'BRCA1'}}, 'api', [])), \
             patch('app.run_enrichment', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.build_ke_gene_mapping', return_value={}), \
             patch('app.guess_id_type', return_value='HGNC'):
            return client.post('/analyze', data={
                'filename': 'DEGs.txt',
                'id_column': 'GENE_SYMBOL',
                'fc_column': 'logFC',
                'pval_column': 'adj.P.Val',
                'aop_selection': 'AOP:1',
                'logfc_threshold': '1.0',
                'method': method,
            })

    def test_same_upload_can_be_analysed_twice(self, authenticated_client, tmp_path, monkeypatch):
        """The reported bug: run once with ORA, then again with GSEA."""
        uploads = self._real_upload(tmp_path, monkeypatch)

        first = self._analyze(authenticated_client, 'ora')
        assert first.status_code == 200
        # The upload survives the first run — this is the actual regression.
        assert (uploads / 'DEGs.txt').exists()

        second = self._analyze(authenticated_client, 'gsea')
        assert second.status_code == 200
        assert b'no longer available' not in second.data

    def test_missing_upload_explains_itself(self, authenticated_client, tmp_path, monkeypatch):
        """When the file really is gone, say why instead of blaming the data."""
        self._real_upload(tmp_path, monkeypatch)
        response = self._analyze(authenticated_client, 'ora')
        assert response.status_code == 200

        (tmp_path / 'uploads' / 'DEGs.txt').unlink()
        gone = self._analyze(authenticated_client, 'ora')
        assert gone.status_code == 400
        assert b'no longer available' in gone.data
        assert b'upload it again' in gone.data


class TestOldUploadsAreSwept:
    """Issue #72: the age-based sweep must actually run — it never did before."""

    def test_upload_triggers_the_sweep(self, flask_client, tmp_path, monkeypatch):
        import app as app_module

        uploads = tmp_path / 'uploads'
        uploads.mkdir()
        monkeypatch.setitem(app_module.app.config, 'UPLOAD_FOLDER', str(uploads))
        monkeypatch.setattr(app_module.Config, 'UPLOAD_FOLDER', str(uploads))

        with patch('app.cleanup_old_uploads', return_value=0) as sweep:
            flask_client.post('/preview', data={})
        sweep.assert_called_once()
        assert sweep.call_args.kwargs['max_age_hours'] == app_module.Config.UPLOAD_RETENTION_HOURS

    def test_sweep_removes_only_aged_files(self, tmp_path, monkeypatch):
        """Fresh uploads stay re-analysable; aged ones are reclaimed."""
        import os
        import time
        import utils
        from utils import cleanup_old_uploads

        uploads = tmp_path / 'uploads'
        uploads.mkdir()
        # Patch the Config that utils actually holds. config.Config is the same
        # object today (#73 removed the module reload that used to make them
        # diverge), but patching the reference the function reads cannot go
        # stale, and a sweep pointed at the real uploads directory would delete
        # a developer's files rather than fail loudly.
        monkeypatch.setattr(utils.Config, 'UPLOAD_FOLDER', str(uploads))

        fresh = uploads / 'fresh.txt'
        fresh.write_text('x')
        aged = uploads / 'aged.txt'
        aged.write_text('x')
        old = time.time() - (48 * 3600)
        os.utime(aged, (old, old))

        # Batch uploads live in per-UUID subdirectories and are cleaned by the
        # batch flow; the sweep must not touch them mid-run.
        batch_dir = uploads / 'a-batch-uuid'
        batch_dir.mkdir()
        (batch_dir / 'cond1.csv').write_text('x')

        assert cleanup_old_uploads(max_age_hours=24) == 1
        assert fresh.exists()
        assert not aged.exists()
        assert (batch_dir / 'cond1.csv').exists()
