"""
Unit tests for report generation service.
"""

import pytest
from services.report_service import report_generator, ReportData, get_software_versions


@pytest.mark.unit
class TestReportService:
    """Test report generation functionality."""
    
    def test_software_versions(self):
        """Test software version collection."""
        versions = get_software_versions()
        
        assert isinstance(versions, dict)
        assert len(versions) > 0
        
        # Common packages should be included
        expected_packages = ['Flask', 'Pandas', 'NumPy', 'SciPy']
        for package in expected_packages:
            assert package in versions
            assert isinstance(versions[package], str)
            assert len(versions[package]) > 0
    
    def test_report_data_creation(self, sample_report_data):
        """Test ReportData dataclass creation."""
        report_data = sample_report_data
        
        # Check required fields
        assert report_data.metadata['dataset_id'] == 'TEST001'
        assert report_data.filename == 'test_data.csv'
        assert report_data.gene_count == 1000
        assert report_data.significant_genes == 150
        assert report_data.aop_id == 'AOP:1'
        
        # Check analysis timestamp is set
        assert report_data.analysis_timestamp is not None
    
    def test_html_report_generation(self, sample_report_data):
        """Test HTML report generation."""
        html_content = report_generator.generate_html_report(sample_report_data)
        
        # Check basic HTML structure
        assert html_content.startswith('<!DOCTYPE html>')
        assert '<html lang="en">' in html_content
        assert '</html>' in html_content
        
        # Check content includes key information
        assert 'TEST001' in html_content  # Dataset ID
        assert 'Test Chemical' in html_content  # Stressor
        assert 'Test AOP Pathway' in html_content  # AOP label
        assert 'Test Key Event' in html_content  # Enrichment results
        
        # Check styling is included
        assert '<style>' in html_content
        assert 'color: #29235C' in html_content  # Brand colors
        assert 'background-color: #307BBF' in html_content
    
    def test_report_sections(self, sample_report_data):
        """Test individual report sections are generated correctly."""
        html_content = report_generator.generate_html_report(sample_report_data)
        
        # Check all major sections are present
        assert 'Experiment Information' in html_content
        assert 'Input Data Summary' in html_content
        assert 'Analysis Settings' in html_content
        assert 'Key Event Enrichment Results' in html_content
        assert 'System Information' in html_content
        
        # Check specific data points
        assert '1000' in html_content  # Gene count
        assert '150' in html_content   # Significant genes
        assert 'HGNC' in html_content  # ID type
        assert '1.0' in html_content   # Threshold
    
    def test_enrichment_table_generation(self, sample_report_data):
        """Test enrichment results table generation."""
        html_content = report_generator.generate_html_report(sample_report_data)
        
        # Check table structure
        assert '<table class="enrichment-table">' in html_content
        assert '<thead>' in html_content
        assert '<tbody>' in html_content
        
        # Check column headers
        assert 'KE ID' in html_content
        assert 'Key Event Title' in html_content
        assert 'P-value' in html_content
        assert 'Adj. P-value' in html_content
        
        # Check data from sample
        assert 'KE:115' in html_content
        assert 'Test Key Event' in html_content
        assert '0.001' in html_content  # p-value
        assert '2.5' in html_content    # enrichment ratio
    
    def test_metadata_section_generation(self, sample_report_data):
        """Test metadata section includes all fields."""
        html_content = report_generator.generate_html_report(sample_report_data)
        
        # Check all metadata fields are displayed
        assert 'TEST001' in html_content        # Dataset ID
        assert 'Test Chemical' in html_content  # Stressor
        assert '10 µM for 24h' in html_content  # Dosing
        assert 'Test User' in html_content      # Owner
        assert 'Sample test experiment' in html_content  # Description
    
    def test_empty_enrichment_results(self):
        """Test report generation with empty enrichment results."""
        from datetime import datetime
        
        empty_report_data = ReportData(
            metadata={'dataset_id': 'EMPTY_TEST'},
            filename='empty_test.csv',
            gene_count=100,
            significant_genes=0,
            aop_id='AOP:1',
            aop_label='Empty Test AOP',
            logfc_threshold=1.0,
            pval_cutoff=0.05,
            id_column='Gene_Symbol',
            fc_column='log2FC',
            pval_column='pvalue',
            id_type='HGNC',
            enrichment_results=[],  # Empty results
            analysis_timestamp=datetime.now()
        )
        
        html_content = report_generator.generate_html_report(empty_report_data)
        
        # Should handle empty results gracefully
        assert 'No enrichment results available' in html_content
        assert 'EMPTY_TEST' in html_content  # Still show metadata
    
    def test_volcano_plot_without_plotly(self, sample_report_data, monkeypatch):
        """Test volcano plot generation when Plotly is not available."""
        # Mock Plotly as unavailable
        monkeypatch.setattr('services.report_service.PLOTLY_AVAILABLE', False)
        
        html_content = report_generator.generate_html_report(sample_report_data)
        
        # Should handle missing Plotly gracefully
        assert 'Volcano plot visualization not available' in html_content
        assert 'TEST001' in html_content  # Other sections should still work
    
    def test_pdf_report_without_weasyprint(self, sample_report_data, monkeypatch):
        """Test PDF generation when WeasyPrint is not available."""
        # Mock WeasyPrint as unavailable
        monkeypatch.setattr('services.report_service.WEASYPRINT_AVAILABLE', False)
        
        with pytest.raises(RuntimeError, match="WeasyPrint is required"):
            report_generator.generate_pdf_report(sample_report_data)
    
    @pytest.mark.slow
    def test_large_enrichment_results(self, sample_metadata):
        """Test report generation with large number of enrichment results."""
        from datetime import datetime
        
        # Create report data with many enrichment results
        large_results = []
        for i in range(100):
            large_results.append({
                'KE_ID': f'KE:{i:03d}',
                'KE_title': f'Key Event {i}',
                'overlap_genes': i % 10 + 1,
                'ke_size': (i % 5 + 1) * 20,
                'pvalue': 0.001 * (i + 1),
                'pvalue_adjusted': 0.01 * (i + 1),
                'enrichment_ratio': 2.0 + (i % 3)
            })
        
        large_report_data = ReportData(
            metadata=sample_metadata.to_dict(),
            filename='large_test.csv',
            gene_count=10000,
            significant_genes=1500,
            aop_id='AOP:1',
            aop_label='Large Test AOP',
            logfc_threshold=1.0,
            pval_cutoff=0.05,
            id_column='Gene_Symbol',
            fc_column='log2FC',
            pval_column='pvalue',
            id_type='HGNC',
            enrichment_results=large_results,
            analysis_timestamp=datetime.now()
        )
        
        html_content = report_generator.generate_html_report(large_report_data)
        
        # Should handle large results (only show top 20 in table)
        assert len(html_content) > 10000  # Should be substantial content
        assert 'KE:000' in html_content    # First result
        assert 'KE:019' in html_content    # 20th result (top 20)
        # KE:020 and beyond might not be in table but could be in summary
    
    def test_report_css_styling(self, sample_report_data):
        """Test that report includes proper CSS styling."""
        html_content = report_generator.generate_html_report(sample_report_data)
        
        # Check brand colors are used
        brand_colors = ['#29235C', '#307BBF', '#93D5F6', '#E6007E']
        for color in brand_colors:
            assert color in html_content
        
        # Check responsive design elements
        assert 'max-width' in html_content
        assert 'grid-template-columns' in html_content
        assert '@media print' in html_content  # Print styles
    
    def test_report_generator_error_handling(self):
        """Test error handling in report generation."""
        # Test with minimal/invalid data
        minimal_data = ReportData(
            metadata={},
            filename='',
            gene_count=0,
            significant_genes=0,
            aop_id='',
            aop_label='',
            logfc_threshold=0.0,
            pval_cutoff=0.05,
            id_column='',
            fc_column='',
            pval_column='',
            id_type='',
            enrichment_results=[]
        )
        
        # Should not crash with minimal data
        html_content = report_generator.generate_html_report(minimal_data)
        assert len(html_content) > 1000  # Should still generate basic structure
        assert '<html' in html_content