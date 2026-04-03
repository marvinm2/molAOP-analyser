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
        """Test main index page loads correctly with both tabs."""
        response = flask_client.get('/')

        assert response.status_code == 200
        assert b'molAOP Analyser' in response.data
        assert b'Select Demo Dataset' in response.data
        assert b'Experiment Information' in response.data  # Metadata form
        assert b'Single Analysis' in response.data  # Single tab
        assert b'Batch Analysis' in response.data  # Batch tab

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