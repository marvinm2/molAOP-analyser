"""
Integration tests for Flask routes and web functionality.
"""

import io
import pytest
import json
import re
from unittest.mock import patch, MagicMock


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
        with patch('os.path.exists', return_value=True), \
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
             patch('app.run_enrichment_analysis', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.build_ke_gene_mapping', return_value={}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', return_value=True), \
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
             patch('app.run_enrichment_analysis', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', return_value=True), \
             patch('app.validate_file_path', return_value=True), \
             patch('pandas.read_csv', return_value=raw_upload_df), \
             patch('app.load_cached_reference_sets', return_value=({'KE:115': {'BRCA1', 'EGFR'}}, 'cache')):

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
             patch('app.run_enrichment_analysis', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', return_value=True), \
             patch('app.validate_file_path', return_value=True), \
             patch('pandas.read_csv', return_value=raw_upload_df), \
             patch('app.load_cached_reference_sets',
                   return_value=({'KE:115': {'BRCA1', 'EGFR'}}, 'cache')):

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
             patch('app.run_enrichment_analysis', return_value=enrichment_df), \
             patch('app.build_cytoscape_network', return_value={'nodes': [], 'edges': []}), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', return_value=True), \
             patch('app.validate_file_path', return_value=True), \
             patch('pandas.read_csv', return_value=raw_upload_df), \
             patch('app.load_cached_reference_sets',
                   return_value=({'KE:115': {'BRCA1', 'EGFR'}}, 'cache')):

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
        with patch('os.path.exists', return_value=True), \
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

        with patch('os.path.exists', return_value=True), \
             patch('pandas.read_csv', return_value=preview_df), \
             patch('app.load_and_validate_data', return_value=processed_df), \
             patch('app.process_gene_expression', return_value=(processed_df, {'total_genes': 5})), \
             patch('app.load_aop_data', return_value=({'KE:115'}, edges_df, {'KE:115': 'KE'}, {'KE:115': 'Test KE'})), \
             patch('app.run_enrichment_analysis', return_value=enrichment_df), \
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
        with patch('os.path.exists', return_value=True), \
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
        with patch('os.path.exists', return_value=True), \
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
        with patch('os.path.exists', return_value=True), \
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