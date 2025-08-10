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
        """Test main index page loads correctly."""
        response = flask_client.get('/')
        
        assert response.status_code == 200
        assert b'Molecular Adverse Outcome Pathway Analyser' in response.data
        assert b'Select Demo Dataset' in response.data
        assert b'Experiment Information' in response.data  # Metadata form
    
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
        with patch('app.load_and_validate_data') as mock_load, \
             patch('app.process_gene_expression') as mock_process, \
             patch('app.load_aop_data') as mock_aop_data, \
             patch('app.run_enrichment_analysis') as mock_enrichment, \
             patch('app.build_cytoscape_network') as mock_network, \
             patch('app.build_ke_gene_mapping') as mock_gene_mapping, \
             patch('app.guess_id_type') as mock_id_type, \
             patch('app.cleanup_file'), \
             patch('os.path.exists', return_value=True), \
             patch('app.validate_file_path', return_value=True):
            
            # Mock service responses
            mock_df = MagicMock()
            mock_df.__getitem__.return_value.__getitem__.return_value = 'HGNC'
            mock_load.return_value = mock_df
            mock_process.return_value = (mock_df, {'total_genes': 1000})
            mock_aop_data.return_value = (['KE:115'], [], {}, {})
            mock_enrichment.return_value = MagicMock()
            mock_enrichment.return_value.to_dict.return_value = []
            mock_network.return_value = {'nodes': [], 'edges': []}
            mock_gene_mapping.return_value = {}
            mock_id_type.return_value = 'HGNC'
            
            response = authenticated_client.post('/analyze', data={
                'filename': 'test.csv',
                'id_column': 'Gene_Symbol',
                'fc_column': 'log2FoldChange',
                'pval_column': 'padj',
                'aop_selection': 'AOP:1',
                'logfc_threshold': '1.0'
            })
            
            assert response.status_code == 200
            assert b'KE Enrichment Results' in response.data or b'results.html' in response.request.url
    
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
            assert response.headers['Content-Type'] == 'text/html; charset=utf-8'
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
        with patch('os.path.exists', return_value=True), \
             patch('pandas.read_csv') as mock_read_csv, \
             patch('app.load_and_validate_data') as mock_load, \
             patch('app.process_gene_expression') as mock_process, \
             patch('app.load_aop_data') as mock_aop_data, \
             patch('app.run_enrichment_analysis') as mock_enrichment, \
             patch('app.build_cytoscape_network') as mock_network, \
             patch('app.build_ke_gene_mapping') as mock_gene_mapping, \
             patch('app.guess_id_type') as mock_id_type, \
             patch('app.cleanup_file'), \
             patch('app.validate_file_path', return_value=True), \
             patch('app.report_generator.generate_html_report') as mock_report:
            
            # Setup mocks
            mock_df = MagicMock()
            mock_df.head.return_value.to_dict.return_value = [{'Gene_Symbol': 'BRCA1'}]
            mock_df.columns.tolist.return_value = ['Gene_Symbol', 'log2FoldChange', 'padj']
            mock_read_csv.return_value = mock_df
            
            mock_load.return_value = mock_df
            mock_process.return_value = (mock_df, {'total_genes': 1000})
            mock_aop_data.return_value = (['KE:115'], [], {}, {})
            mock_enrichment.return_value = MagicMock()
            mock_enrichment.return_value.to_dict.return_value = []
            mock_network.return_value = {'nodes': [], 'edges': []}
            mock_gene_mapping.return_value = {}
            mock_id_type.return_value = 'HGNC'
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
                'dataset_id': 'INTEGRATION_TEST'  # Carried from step 1
            })
            assert response2.status_code == 200
            
            # Step 3: Generate report
            response3 = flask_client.post('/generate_report', data={
                'format': 'html',
                'filename': 'GSE90122_TO90137.tsv',
                'gene_count': '1000',
                'aop_id': 'AOP:1'
            })
            assert response3.status_code == 200
            assert b'Integration Test Report' in response3.data