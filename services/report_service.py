"""
Report generation service for the Molecular AOP Analyser.

Generates comprehensive PDF and HTML reports containing:
- Input data summary (gene count, volcano plot, filename)
- Analysis settings and parameters
- KE enrichment results (tables and visualizations)
- Software version information
"""

import os
import json
import base64
import logging
import math
from datetime import datetime
from typing import Dict, List, Any, Optional, Union
from dataclasses import dataclass
from io import BytesIO
import pandas as pd

try:
    import plotly.graph_objects as go
    import plotly.io as pio
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    
try:
    from weasyprint import HTML, CSS
    from weasyprint.html import get_html_document
    WEASYPRINT_AVAILABLE = True
except (ImportError, OSError) as e:
    WEASYPRINT_AVAILABLE = False
    HTML = None
    CSS = None

from config import Config, ExperimentMetadata

logger = logging.getLogger(__name__)

@dataclass
class ReportData:
    """Container for all data needed to generate a comprehensive report."""
    
    # Experiment metadata
    metadata: Dict[str, Any]
    
    # Input data summary
    filename: str
    gene_count: int
    significant_genes: int
    
    # Analysis parameters
    aop_id: str
    aop_label: str
    logfc_threshold: float
    pval_cutoff: float
    id_column: str
    fc_column: str
    pval_column: str
    id_type: str
    
    # Results
    enrichment_results: List[Dict[str, Any]]
    volcano_data: Optional[List[Dict[str, Any]]] = None
    network_data: Optional[Dict[str, Any]] = None
    
    # System information
    analysis_timestamp: datetime = None
    software_versions: Optional[Dict[str, str]] = None
    
    def __post_init__(self):
        """Set analysis timestamp if not provided."""
        if self.analysis_timestamp is None:
            self.analysis_timestamp = datetime.now()


class ReportGenerator:
    """Main report generation class supporting HTML and PDF formats."""
    
    def __init__(self):
        """Initialize the report generator."""
        self.template_dir = "templates"
        self.static_dir = "static"
        
    def generate_html_report(self, report_data: ReportData) -> str:
        """Generate an HTML report from the provided data.
        
        Args:
            report_data: Complete report data container
            
        Returns:
            HTML string containing the full report
        """
        try:
            return self._build_html_content(report_data)
        except Exception as e:
            logger.error(f"Failed to generate HTML report: {e}")
            raise
    
    def generate_pdf_report(self, report_data: ReportData) -> bytes:
        """Generate a PDF report from the provided data.
        
        Args:
            report_data: Complete report data container
            
        Returns:
            PDF content as bytes
            
        Raises:
            RuntimeError: If WeasyPrint is not available
        """
        if not WEASYPRINT_AVAILABLE:
            raise RuntimeError("WeasyPrint is required for PDF generation. Install with: pip install weasyprint")
        
        try:
            html_content = self._build_html_content(report_data)
            
            # Create PDF from HTML
            html_doc = HTML(string=html_content)
            pdf_bytes = html_doc.write_pdf()
            
            logger.info("PDF report generated successfully")
            return pdf_bytes
            
        except Exception as e:
            logger.error(f"Failed to generate PDF report: {e}")
            raise
    
    def _build_html_content(self, report_data: ReportData) -> str:
        """Build the complete HTML report content.
        
        Args:
            report_data: Report data container
            
        Returns:
            Complete HTML string
        """
        # Generate individual sections
        header_html = self._generate_header(report_data)
        metadata_html = self._generate_metadata_section(report_data)
        input_summary_html = self._generate_input_summary(report_data)
        analysis_params_html = self._generate_analysis_parameters(report_data)
        volcano_html = self._generate_volcano_section(report_data)
        enrichment_html = self._generate_enrichment_section(report_data)
        system_info_html = self._generate_system_info(report_data)
        
        # Combine all sections
        full_html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Molecular AOP Analysis Report</title>
            <style>
                {self._get_report_css()}
            </style>
        </head>
        <body>
            {header_html}
            <div class="container">
                {metadata_html}
                {input_summary_html}
                {analysis_params_html}
                {volcano_html}
                {enrichment_html}
                {system_info_html}
            </div>
        </body>
        </html>
        """
        
        return full_html
    
    def _generate_header(self, report_data: ReportData) -> str:
        """Generate report header section."""
        return f"""
        <header class="report-header">
            <h1>Molecular AOP Analysis Report</h1>
            <div class="report-meta">
                <p><strong>Dataset:</strong> {report_data.metadata.get('dataset_id', 'N/A')}</p>
                <p><strong>Generated:</strong> {report_data.analysis_timestamp.strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
        </header>
        """
    
    def _generate_metadata_section(self, report_data: ReportData) -> str:
        """Generate experiment metadata section."""
        metadata = report_data.metadata
        
        return f"""
        <section class="metadata-section">
            <h2>Experiment Information</h2>
            <div class="metadata-grid">
                <div class="metadata-item">
                    <label>Dataset ID:</label>
                    <span>{metadata.get('dataset_id', 'N/A')}</span>
                </div>
                <div class="metadata-item">
                    <label>Stressor/Treatment:</label>
                    <span>{metadata.get('stressor', 'N/A')}</span>
                </div>
                <div class="metadata-item">
                    <label>Dosing Information:</label>
                    <span>{metadata.get('dosing', 'N/A')}</span>
                </div>
                <div class="metadata-item">
                    <label>Owner/Contact:</label>
                    <span>{metadata.get('owner', 'N/A')}</span>
                </div>
                <div class="metadata-item">
                    <label>Upload Time:</label>
                    <span>{metadata.get('upload_timestamp', 'N/A')}</span>
                </div>
                {f'<div class="metadata-item full-width"><label>Description:</label><p>{metadata.get("description", "N/A")}</p></div>' if metadata.get('description') else ''}
            </div>
        </section>
        """
    
    def _generate_input_summary(self, report_data: ReportData) -> str:
        """Generate input data summary section."""
        return f"""
        <section class="input-summary">
            <h2>Input Data Summary</h2>
            <div class="summary-stats">
                <div class="stat-item">
                    <label>Filename:</label>
                    <span>{report_data.filename}</span>
                </div>
                <div class="stat-item">
                    <label>Total Genes:</label>
                    <span>{report_data.gene_count:,}</span>
                </div>
                <div class="stat-item">
                    <label>Significant Genes:</label>
                    <span>{report_data.significant_genes:,}</span>
                </div>
                <div class="stat-item">
                    <label>Gene ID Type:</label>
                    <span>{report_data.id_type}</span>
                </div>
            </div>
        </section>
        """
    
    def _generate_analysis_parameters(self, report_data: ReportData) -> str:
        """Generate analysis parameters section."""
        return f"""
        <section class="analysis-params">
            <h2>Analysis Settings</h2>
            <div class="params-grid">
                <div class="param-item">
                    <label>Selected AOP:</label>
                    <span>{report_data.aop_label} (ID: {report_data.aop_id})</span>
                </div>
                <div class="param-item">
                    <label>Log2FC Threshold:</label>
                    <span>{report_data.logfc_threshold}</span>
                </div>
                <div class="param-item">
                    <label>P-value Cutoff:</label>
                    <span>{report_data.pval_cutoff}</span>
                </div>
                <div class="param-item">
                    <label>Gene ID Column:</label>
                    <span>{report_data.id_column}</span>
                </div>
                <div class="param-item">
                    <label>Log2FC Column:</label>
                    <span>{report_data.fc_column}</span>
                </div>
                <div class="param-item">
                    <label>P-value Column:</label>
                    <span>{report_data.pval_column}</span>
                </div>
            </div>
        </section>
        """
    
    def _generate_volcano_section(self, report_data: ReportData) -> str:
        """Generate volcano plot section."""
        if not report_data.volcano_data or not PLOTLY_AVAILABLE:
            return """
            <section class="volcano-section">
                <h2>Volcano Plot</h2>
                <p class="note">Volcano plot visualization not available.</p>
            </section>
            """
        
        try:
            # Create volcano plot using Plotly
            volcano_html = self._create_volcano_plot(report_data.volcano_data, 
                                                   report_data.logfc_threshold,
                                                   report_data.pval_cutoff)
            
            return f"""
            <section class="volcano-section">
                <h2>Volcano Plot</h2>
                <div class="plot-description">
                    <p>Volcano plot showing log2 fold change vs. -log10(p-value) for all genes. 
                    Red points indicate significantly upregulated genes, blue points show significantly downregulated genes.</p>
                </div>
                <div class="plot-container">
                    {volcano_html}
                </div>
            </section>
            """
            
        except Exception as e:
            logger.warning(f"Failed to generate volcano plot: {e}")
            return """
            <section class="volcano-section">
                <h2>Volcano Plot</h2>
                <p class="note">Error generating volcano plot visualization.</p>
            </section>
            """
    
    def _generate_enrichment_section(self, report_data: ReportData) -> str:
        """Generate Key Event enrichment results section."""
        if not report_data.enrichment_results:
            return """
            <section class="enrichment-section">
                <h2>Key Event Enrichment Results</h2>
                <p class="note">No enrichment results available.</p>
            </section>
            """
        
        # Create table rows
        table_rows = ""
        for result in report_data.enrichment_results[:20]:  # Show top 20 results
            pval = result.get('pvalue', 0)
            adj_pval = result.get('pvalue_adjusted', pval)
            
            table_rows += f"""
            <tr class="{'significant' if adj_pval < 0.05 else ''}">
                <td>{result.get('KE_ID', 'N/A')}</td>
                <td>{result.get('KE_title', 'N/A')}</td>
                <td>{result.get('overlap_genes', 0)}</td>
                <td>{result.get('ke_size', 0)}</td>
                <td>{pval:.2e}</td>
                <td>{adj_pval:.2e}</td>
                <td>{result.get('enrichment_ratio', 0):.2f}</td>
            </tr>
            """
        
        return f"""
        <section class="enrichment-section">
            <h2>Key Event Enrichment Results</h2>
            <p class="section-description">
                Statistical enrichment analysis of significantly differentially expressed genes 
                against Key Event gene sets. Results are ranked by adjusted p-value.
            </p>
            
            <div class="table-container">
                <table class="enrichment-table">
                    <thead>
                        <tr>
                            <th>KE ID</th>
                            <th>Key Event Title</th>
                            <th>Overlap</th>
                            <th>KE Size</th>
                            <th>P-value</th>
                            <th>Adj. P-value</th>
                            <th>Enrichment Ratio</th>
                        </tr>
                    </thead>
                    <tbody>
                        {table_rows}
                    </tbody>
                </table>
            </div>
            
            <div class="results-summary">
                <p><strong>Total Key Events tested:</strong> {len(report_data.enrichment_results)}</p>
                <p><strong>Significantly enriched (adj. p < 0.05):</strong> 
                   {len([r for r in report_data.enrichment_results if r.get('pvalue_adjusted', 1) < 0.05])}</p>
            </div>
        </section>
        """
    
    def _generate_system_info(self, report_data: ReportData) -> str:
        """Generate system information and software versions section."""
        versions = report_data.software_versions or {}
        
        version_items = ""
        for package, version in versions.items():
            version_items += f"<li><strong>{package}:</strong> {version}</li>"
        
        return f"""
        <section class="system-info">
            <h2>System Information</h2>
            <div class="info-grid">
                <div class="info-item">
                    <label>Analysis Date:</label>
                    <span>{report_data.analysis_timestamp.strftime('%Y-%m-%d %H:%M:%S UTC')}</span>
                </div>
                <div class="info-item">
                    <label>Application:</label>
                    <span>Molecular AOP Analyser</span>
                </div>
            </div>
            
            {f'<div class="versions"><h3>Software Versions</h3><ul>{version_items}</ul></div>' if versions else ''}
        </section>
        """
    
    def _create_volcano_plot(self, volcano_data: List[Dict], threshold: float, pval_cutoff: float) -> str:
        """Create a volcano plot using Plotly and return as HTML string."""
        if not PLOTLY_AVAILABLE:
            return "<p>Plotly not available for plot generation.</p>"
        
        try:
            # Extract data for plotting
            x_values = [d['log2FC'] for d in volcano_data]
            y_values = [-math.log10(d['pval']) for d in volcano_data if d['pval'] > 0]
            gene_names = [d['ID'] for d in volcano_data]
            
            # Determine colors based on significance
            colors = []
            for d in volcano_data:
                fc = d['log2FC']
                pval = d['pval']
                
                if pval < pval_cutoff:
                    if fc >= threshold:
                        colors.append('red')  # Upregulated
                    elif fc <= -threshold:
                        colors.append('blue')  # Downregulated
                    else:
                        colors.append('green')  # Significant but low FC
                else:
                    colors.append('lightgray')  # Not significant
            
            # Create Plotly figure
            fig = go.Figure(data=go.Scatter(
                x=x_values,
                y=y_values,
                mode='markers',
                marker=dict(color=colors, size=4),
                text=gene_names,
                hovertemplate='Gene: %{text}<br>log2FC: %{x}<br>-log10(p): %{y}<extra></extra>'
            ))
            
            # Add threshold lines
            fig.add_hline(y=-math.log10(pval_cutoff), line_dash="dash", line_color="black")
            if threshold > 0:
                fig.add_vline(x=threshold, line_dash="dot", line_color="black")
                fig.add_vline(x=-threshold, line_dash="dot", line_color="black")
            
            fig.update_layout(
                title="Volcano Plot",
                xaxis_title="log2 Fold Change",
                yaxis_title="-log10(p-value)",
                width=700,
                height=500,
                showlegend=False
            )
            
            return fig.to_html(include_plotlyjs=False, div_id="volcano-plot")
            
        except Exception as e:
            logger.warning(f"Failed to create volcano plot: {e}")
            return f"<p>Error creating volcano plot: {e}</p>"
    
    def _get_report_css(self) -> str:
        """Get CSS styles for the report."""
        return """
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 0;
            color: #29235C;
            background-color: #f8f9fa;
        }
        
        .container {
            max-width: 1000px;
            margin: 0 auto;
            padding: 20px;
        }
        
        .report-header {
            background-color: #307BBF;
            color: white;
            padding: 30px;
            text-align: center;
            margin-bottom: 30px;
        }
        
        .report-header h1 {
            margin: 0 0 15px 0;
            font-size: 28px;
        }
        
        .report-meta p {
            margin: 5px 0;
            font-size: 14px;
        }
        
        section {
            background: white;
            margin-bottom: 30px;
            padding: 25px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        
        h2 {
            color: #307BBF;
            border-bottom: 2px solid #93D5F6;
            padding-bottom: 10px;
            margin-top: 0;
        }
        
        h3 {
            color: #29235C;
            margin-top: 20px;
        }
        
        .metadata-grid, .params-grid, .info-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }
        
        .metadata-item, .param-item, .info-item, .stat-item {
            display: flex;
            flex-direction: column;
        }
        
        .metadata-item.full-width {
            grid-column: 1 / -1;
        }
        
        .metadata-item label, .param-item label, .info-item label, .stat-item label {
            font-weight: bold;
            color: #29235C;
            margin-bottom: 5px;
        }
        
        .summary-stats {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }
        
        .table-container {
            overflow-x: auto;
            margin-top: 15px;
        }
        
        .enrichment-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }
        
        .enrichment-table th, .enrichment-table td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }
        
        .enrichment-table th {
            background-color: #307BBF;
            color: white;
            font-weight: bold;
        }
        
        .enrichment-table tr.significant {
            background-color: #fff3cd;
        }
        
        .enrichment-table tr:nth-child(even):not(.significant) {
            background-color: #f8f9fa;
        }
        
        .plot-container {
            margin-top: 15px;
            text-align: center;
        }
        
        .plot-description {
            font-style: italic;
            color: #666;
            margin-bottom: 15px;
        }
        
        .section-description {
            color: #666;
            font-style: italic;
            margin-bottom: 15px;
        }
        
        .results-summary {
            margin-top: 20px;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 5px;
            border-left: 4px solid #307BBF;
        }
        
        .note {
            color: #666;
            font-style: italic;
        }
        
        .versions ul {
            list-style-type: none;
            padding-left: 0;
        }
        
        .versions li {
            margin-bottom: 5px;
        }
        
        @media print {
            body { background-color: white; }
            section { box-shadow: none; border: 1px solid #ddd; }
            .report-header { background-color: #307BBF !important; }
        }
        """


# Global report generator instance
report_generator = ReportGenerator()


def get_software_versions() -> Dict[str, str]:
    """Get versions of key software packages."""
    versions = {}
    
    try:
        import flask
        versions['Flask'] = flask.__version__
    except (ImportError, AttributeError):
        pass
    
    try:
        import pandas
        versions['Pandas'] = pandas.__version__
    except (ImportError, AttributeError):
        pass
    
    try:
        import numpy
        versions['NumPy'] = numpy.__version__
    except (ImportError, AttributeError):
        pass
    
    try:
        import scipy
        versions['SciPy'] = scipy.__version__
    except (ImportError, AttributeError):
        pass
    
    if PLOTLY_AVAILABLE:
        try:
            import plotly
            versions['Plotly'] = plotly.__version__
        except AttributeError:
            pass
    
    return versions