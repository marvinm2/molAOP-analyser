"""
Unit tests for report generation service.
"""

from dataclasses import replace

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
        
        # Common packages should be included. statsmodels and gseapy are the
        # packages performing the inferential work (BH correction and GSEA) and
        # are required for the report to be reproducible — see issue #66.
        expected_packages = [
            'Flask', 'Pandas', 'NumPy', 'SciPy', 'statsmodels', 'gseapy', 'Python',
        ]
        for package in expected_packages:
            assert package in versions
            assert isinstance(versions[package], str)
            assert len(versions[package]) > 0

    def test_software_versions_records_image_sha(self, monkeypatch):
        """The image SHA is reported when baked in, and omitted when absent."""
        monkeypatch.setenv('MOLAOP_IMAGE_SHA', 'abc1234')
        assert get_software_versions()['Image'] == 'abc1234'

        monkeypatch.delenv('MOLAOP_IMAGE_SHA', raising=False)
        assert 'Image' not in get_software_versions()


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
        assert html_content.strip().startswith('<!DOCTYPE html>')
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
        assert 'FDR' in html_content

        # Check data from sample
        assert 'KE:115' in html_content
        assert 'Test Key Event' in html_content
        assert '0.0010' in html_content  # p-value formatted
        assert '2.50' in html_content    # odds ratio formatted
    
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
        """Test PDF generation falls back to ReportLab when WeasyPrint unavailable."""
        monkeypatch.setattr('services.report_service.WEASYPRINT_AVAILABLE', False)

        # Should fall back to ReportLab, not raise
        pdf_bytes = report_generator.generate_pdf_report(sample_report_data)
        assert pdf_bytes is not None
        assert len(pdf_bytes) > 0
    
    @pytest.mark.slow
    def test_large_enrichment_results(self, sample_metadata):
        """Test report generation with large number of enrichment results."""
        from datetime import datetime
        
        # Create report data with many enrichment results
        large_results = []
        for i in range(100):
            large_results.append({
                'KE': f'KE:{i:03d}',
                'Title': f'Key Event {i}',
                'num_overlap': i % 10 + 1,
                'total_KE_genes_in_dataset': (i % 5 + 1) * 20,
                'p_value': 0.001 * (i + 1),
                'FDR': 0.01 * (i + 1),
                'odds_ratio': 2.0 + (i % 3),
                'pct_sig_in_KE': 25.0,
                'overlap': 'GENE1, GENE2',
                'sig_in_KE': i % 10 + 1,
                'sig_not_KE': 100,
                'non_sig_in_KE': 50,
                'non_sig_not_KE': 800,
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
        
        # Check primary brand colors used in report CSS
        assert '#29235C' in html_content  # Primary dark (text)
        assert '#307BBF' in html_content  # Primary blue (header)
        
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


class TestMinConfidenceInReport:
    """Issue #60: the chosen mapping-confidence threshold is reported."""

    def test_default_label_is_all_mappings(self, sample_report_data):
        assert sample_report_data.min_confidence == 'all'
        assert sample_report_data.min_confidence_label == 'All mappings'

    def test_labels_for_each_threshold(self, sample_report_data):
        sample_report_data.min_confidence = 'medium'
        assert sample_report_data.min_confidence_label == 'Medium and High only'
        sample_report_data.min_confidence = 'high'
        assert sample_report_data.min_confidence_label == 'High only'

    def test_threshold_rendered_in_the_html_report(self, sample_report_data):
        sample_report_data.min_confidence = 'high'
        html_content = report_generator.generate_html_report(sample_report_data)
        assert 'Minimum Mapping Confidence' in html_content
        assert 'High only' in html_content


@pytest.mark.unit
class TestKeAccountingInReport:
    """Issue #65: the report states the multiple-testing denominator."""

    _SUMMARY = {
        'total_kes': 24, 'tested': 18, 'excluded_no_mapping': 2,
        'excluded_too_few_genes': 4, 'excluded_error': 0, 'min_ke_genes': 5,
    }

    def _with_summary(self, sample_report_data, summary):
        return replace(sample_report_data, ke_summary=summary)

    def test_ke_summary_text_property(self, sample_report_data):
        rd = self._with_summary(sample_report_data, self._SUMMARY)
        assert rd.ke_summary_text == (
            '18 of 24 Key Events tested; 4 excluded (fewer than 5 measured genes), '
            '2 excluded (no gene set mapped)'
        )

    def test_html_report_includes_accounting(self, sample_report_data):
        html = report_generator.generate_html_report(
            self._with_summary(sample_report_data, self._SUMMARY)
        )
        assert 'Key Events tested' in html
        assert '18 of 24 Key Events tested' in html
        assert '4 excluded (fewer than 5 measured genes)' in html

    def test_html_report_omits_accounting_when_absent(self, sample_report_data):
        """Back-compat: ReportData without the summary renders as before."""
        assert sample_report_data.ke_summary is None
        html = report_generator.generate_html_report(sample_report_data)
        assert 'Key Events tested:</label>' not in html

    def test_pdf_report_generates_with_accounting(self, sample_report_data):
        """The PDF path must survive the extra parameter row."""
        pdf = report_generator._generate_reportlab_pdf(
            self._with_summary(sample_report_data, self._SUMMARY)
        )
        assert pdf.startswith(b'%PDF')

    def test_html_significance_call_uses_fdr(self, sample_report_data):
        """Issue #63: the enrichment table highlights rows on FDR, not raw p."""
        rd = replace(sample_report_data, enrichment_results=[
            dict(sample_report_data.enrichment_results[0], KE='KE:RAWONLY',
                 Title='Raw p only', p_value=0.01, FDR=0.30),
        ])
        html = report_generator.generate_html_report(rd)
        # The single data row must not be marked significant on its raw p.
        assert '<tr class="significant">' not in html
