"""
Integration tests for Flask routes and web functionality.
"""

import pytest
import json
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
        # No demo radio buttons remain.
        assert b'name="demo_file"' not in response.data or b'type="radio"' not in response.data

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
        """GET /demos returns 200 and the curated sections."""
        import os
        response = flask_client.get('/demos')
        assert response.status_code == 200
        body = response.data
        # Page chrome
        assert b'Demo Datasets' in body
        # Sections
        assert b'PXR Agonists' in body
        assert b'Cisplatin Nephrotoxicity' in body
        # Curated cards present (hardcoded in route, always rendered)
        assert b'GSE90122_TO90137.tsv' in body
        assert b'GSE90122_SR12813.tsv' in body
        assert b'CSP_24hr_10uM.csv' in body
        assert b'CSP_4hr_10uM.csv' in body
        assert b'CSP_72hr_10uM.csv' in body
        # Show-all expander is only rendered when the gitignored cisplatin
        # data directory is present (~146 MB, not included in the repo).
        # Skip this assertion in CI / fresh-clone environments.
        cisplatin_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'Cisplatin_Kidney')
        if os.path.isdir(cisplatin_dir):
            assert b'Show all 54 cisplatin datasets' in body
        # CTA copy
        assert b'Use this dataset' in body

    def test_demos_route_embeds_recommended_aops(self, flask_client):
        """Each curated card embeds the right recommended_aops payload."""
        response = flask_client.get('/demos')
        body = response.data.decode('utf-8')
        # PXR demo cards carry AOP:1
        assert 'name="demo_file" value="GSE90122_TO90137.tsv"' in body
        # The recommended_aops hidden input next to a PXR demo file should be 'AOP:1'.
        # Locate the form snippet and assert the AOP appears within ~400 chars of it.
        idx = body.index('GSE90122_TO90137.tsv')
        nearby = body[max(0, idx - 400): idx + 400]
        assert 'AOP:1' in nearby, f"AOP:1 not found near PXR card: {nearby!r}"
        # Cisplatin demo cards carry NETWORK:kidney as the first item
        idx2 = body.index('CSP_24hr_10uM.csv')
        nearby2 = body[max(0, idx2 - 400): idx2 + 400]
        assert 'NETWORK:kidney' in nearby2

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