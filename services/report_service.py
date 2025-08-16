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

try:
    from reportlab.lib.pagesizes import letter, A4
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak, Image
    from reportlab.platypus.flowables import HRFlowable
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.colors import Color, HexColor
    from reportlab.lib.units import inch
    from reportlab.lib import colors
    try:
        from reportlab.graphics import renderPDF
        from reportlab.graphics.shapes import Drawing
        from svglib.svglib import renderSVG
        SVG_SUPPORT = True
    except ImportError:
        SVG_SUPPORT = False
    REPORTLAB_AVAILABLE = True
except ImportError:
    REPORTLAB_AVAILABLE = False
    SVG_SUPPORT = False

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
    network_png: Optional[str] = None  # Base64 encoded PNG/JPEG data or SVG string
    
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
            RuntimeError: If no PDF generation method is available
        """
        # Try ReportLab first (pure Python, no system dependencies)
        if REPORTLAB_AVAILABLE:
            try:
                logger.info("Generating PDF using ReportLab")
                return self._generate_reportlab_pdf(report_data)
            except Exception as e:
                logger.warning(f"ReportLab PDF generation failed: {e}")
                # Fall through to WeasyPrint if available
        
        # Fallback to WeasyPrint if ReportLab fails
        if WEASYPRINT_AVAILABLE:
            try:
                logger.info("Falling back to WeasyPrint for PDF generation")
                html_content = self._build_html_content(report_data)
                html_doc = HTML(string=html_content)
                pdf_bytes = html_doc.write_pdf(optimize_images=True)
                
                if not pdf_bytes:
                    raise RuntimeError("WeasyPrint returned empty content")
                
                logger.info(f"WeasyPrint PDF generated successfully ({len(pdf_bytes)} bytes)")
                return pdf_bytes
                
            except Exception as e:
                logger.error(f"WeasyPrint PDF generation failed: {e}")
        
        # No PDF generation available
        error_msg = (
            "PDF generation not available. "
            "ReportLab failed and WeasyPrint requires system dependencies. "
            "Please try the HTML report format as an alternative."
        )
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
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
            # Map actual field names from enrichment service
            pval = result.get('p_value', 0)
            adj_pval = result.get('FDR', pval)
            ke_id = result.get('KE', 'N/A')
            title = result.get('Title', 'N/A')
            overlap_count = result.get('num_overlap', 0)
            ke_size = result.get('total_KE_genes_in_dataset', 0)
            odds_ratio = result.get('odds_ratio', 0)
            
            # Format odds ratio properly
            odds_ratio_str = f"{odds_ratio:.2f}" if isinstance(odds_ratio, (int, float)) and odds_ratio != 'NA' else str(odds_ratio)
            
            table_rows += f"""
            <tr class="{'significant' if adj_pval < 0.05 else ''}">
                <td>{ke_id}</td>
                <td>{title}</td>
                <td>{overlap_count}</td>
                <td>{ke_size}</td>
                <td>{pval:.2e if pval < 0.001 else pval:.4f}</td>
                <td>{adj_pval:.2e if adj_pval < 0.001 else adj_pval:.4f}</td>
                <td>{odds_ratio_str}</td>
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
                            <th># Overlap</th>
                            <th>KE Gene Count</th>
                            <th>P-value</th>
                            <th>FDR</th>
                            <th>Odds Ratio</th>
                        </tr>
                    </thead>
                    <tbody>
                        {table_rows}
                    </tbody>
                </table>
            </div>
            
            <div class="results-summary">
                <p><strong>Total Key Events tested:</strong> {len(report_data.enrichment_results)}</p>
                <p><strong>Significantly enriched (FDR < 0.05):</strong> 
                   {len([r for r in report_data.enrichment_results if r.get('FDR', 1) < 0.05])}</p>
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
    
    def _create_volcano_plot_image(self, volcano_data: List[Dict], threshold: float, pval_cutoff: float) -> bytes:
        """Create a volcano plot using Plotly and return as PNG image bytes for PDF inclusion."""
        if not PLOTLY_AVAILABLE:
            raise RuntimeError("Plotly not available for volcano plot generation")
        
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
                marker=dict(color=colors, size=3),
                text=gene_names,
                hovertemplate='Gene: %{text}<br>log2FC: %{x}<br>-log10(p): %{y}<extra></extra>'
            ))
            
            # Add threshold lines
            fig.add_hline(y=-math.log10(pval_cutoff), line_dash="dash", line_color="black", line_width=1)
            if threshold > 0:
                fig.add_vline(x=threshold, line_dash="dot", line_color="black", line_width=1)
                fig.add_vline(x=-threshold, line_dash="dot", line_color="black", line_width=1)
            
            fig.update_layout(
                title="Volcano Plot - Gene Expression Analysis",
                xaxis_title="log2 Fold Change",
                yaxis_title="-log10(p-value)",
                width=600,
                height=400,
                showlegend=False,
                plot_bgcolor='white',
                paper_bgcolor='white',
                font=dict(size=10)
            )
            
            fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
            fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
            
            # Export as PNG image bytes
            img_bytes = fig.to_image(format="png", width=600, height=400, scale=2)
            return img_bytes
            
        except Exception as e:
            logger.warning(f"Failed to create volcano plot image: {e}")
            raise
    
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
    
    def _generate_reportlab_pdf(self, report_data: ReportData) -> bytes:
        """Generate a PDF report using ReportLab (pure Python, no system dependencies).
        
        Args:
            report_data: Complete report data container
            
        Returns:
            PDF content as bytes
        """
        if not REPORTLAB_AVAILABLE:
            raise RuntimeError("ReportLab not available")
            
        buffer = BytesIO()
        
        # Create document with custom page size and margins
        doc = SimpleDocTemplate(
            buffer,
            pagesize=A4,
            rightMargin=72,
            leftMargin=72,
            topMargin=72,
            bottomMargin=18
        )
        
        # Get style sheet and define custom styles
        styles = getSampleStyleSheet()
        
        # Define house style colors
        primary_color = HexColor('#307BBF')
        dark_color = HexColor('#29235C')
        accent_color = HexColor('#E6007E')
        
        # Custom styles
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=24,
            spaceAfter=30,
            textColor=dark_color,
            alignment=1  # Center
        )
        
        heading_style = ParagraphStyle(
            'CustomHeading',
            parent=styles['Heading2'],
            fontSize=16,
            spaceAfter=12,
            textColor=primary_color,
            borderWidth=1,
            borderColor=primary_color,
            borderPadding=5
        )
        
        normal_style = ParagraphStyle(
            'CustomNormal',
            parent=styles['Normal'],
            fontSize=10,
            textColor=dark_color
        )
        
        # Build document content
        story = []
        
        # Header section
        story.append(Paragraph("Molecular AOP Analysis Report", title_style))
        story.append(Spacer(1, 12))
        
        # Report metadata
        metadata_text = f"""
        <b>Dataset:</b> {report_data.metadata.get('dataset_id', 'N/A')}<br/>
        <b>Generated:</b> {report_data.analysis_timestamp.strftime('%Y-%m-%d %H:%M:%S')}<br/>
        <b>AOP:</b> {report_data.aop_label}
        """
        story.append(Paragraph(metadata_text, normal_style))
        story.append(Spacer(1, 20))
        
        # Input Data Summary
        story.append(Paragraph("Input Data Summary", heading_style))
        
        summary_data = [
            ['Metric', 'Value'],
            ['Filename', report_data.filename],
            ['Total Genes', f"{report_data.gene_count:,}"],
            ['Significant Genes', f"{report_data.significant_genes:,}"],
            ['Gene ID Type', report_data.id_type],
            ['Log2FC Threshold', str(report_data.logfc_threshold)],
            ['P-value Cutoff', str(report_data.pval_cutoff)]
        ]
        
        summary_table = Table(summary_data, colWidths=[2.5*inch, 3.5*inch])
        summary_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), primary_color),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ]))
        
        story.append(summary_table)
        story.append(Spacer(1, 20))
        
        # Key Event Enrichment Results
        logger.info(f"ReportLab: enrichment_results present: {bool(report_data.enrichment_results)}")
        if report_data.enrichment_results:
            logger.info(f"ReportLab: Processing {len(report_data.enrichment_results)} enrichment results")
            story.append(Paragraph("Key Event Enrichment Results", heading_style))
            story.append(Paragraph("Statistical enrichment analysis results ranked by FDR", normal_style))
            story.append(Spacer(1, 12))
            
            # Prepare table data (show top 15 results for space)
            table_data = [['KE ID', 'Key Event Title', '# Overlap', 'P-value', 'FDR', 'Odds Ratio']]
            
            for result in report_data.enrichment_results[:15]:
                pval = result.get('p_value', 0)
                fdr = result.get('FDR', pval)
                odds_ratio = result.get('odds_ratio', 0)
                
                # Format values
                pval_str = f"{pval:.2e}" if pval < 0.001 else f"{pval:.4f}"
                fdr_str = f"{fdr:.2e}" if fdr < 0.001 else f"{fdr:.4f}"
                odds_str = f"{odds_ratio:.2f}" if isinstance(odds_ratio, (int, float)) and odds_ratio != 'NA' else str(odds_ratio)
                
                title = result.get('Title', 'N/A')
                if len(title) > 40:
                    title = title[:37] + "..."
                
                table_data.append([
                    result.get('KE', 'N/A'),
                    title,
                    str(result.get('num_overlap', 0)),
                    pval_str,
                    fdr_str,
                    odds_str
                ])
            
            # Create table with appropriate column widths
            col_widths = [0.8*inch, 2.2*inch, 0.8*inch, 0.8*inch, 0.8*inch, 0.8*inch]
            enrichment_table = Table(table_data, colWidths=col_widths)
            
            # Style the table
            table_style = [
                ('BACKGROUND', (0, 0), (-1, 0), primary_color),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 9),
                ('FONTSIZE', (0, 1), (-1, -1), 8),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 8),
                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ]
            
            # Highlight significant results
            for i, result in enumerate(report_data.enrichment_results[:15], 1):
                if result.get('FDR', 1) < 0.05:
                    table_style.append(('BACKGROUND', (0, i), (-1, i), colors.lightyellow))
            
            enrichment_table.setStyle(TableStyle(table_style))
            story.append(enrichment_table)
            
            # Summary statistics
            story.append(Spacer(1, 12))
            sig_count = len([r for r in report_data.enrichment_results if r.get('FDR', 1) < 0.05])
            summary_text = f"""
            <b>Total Key Events tested:</b> {len(report_data.enrichment_results)}<br/>
            <b>Significantly enriched (FDR &lt; 0.05):</b> {sig_count}<br/>
            <i>Note: Only top 15 results shown in table above</i>
            """
            story.append(Paragraph(summary_text, normal_style))
        
        # Volcano Plot Section
        logger.info(f"ReportLab: volcano_data present: {bool(report_data.volcano_data)}")
        if report_data.volcano_data:
            logger.info(f"ReportLab: volcano_data length: {len(report_data.volcano_data)}")
        if report_data.volcano_data and PLOTLY_AVAILABLE:
            try:
                story.append(Spacer(1, 30))
                story.append(Paragraph("Volcano Plot", heading_style))
                story.append(Paragraph("Gene expression analysis showing log2 fold change vs. statistical significance", normal_style))
                story.append(Spacer(1, 12))
                
                # Generate volcano plot image
                img_bytes = self._create_volcano_plot_image(
                    report_data.volcano_data,
                    report_data.logfc_threshold,
                    report_data.pval_cutoff
                )
                
                # Create image from bytes
                img_buffer = BytesIO(img_bytes)
                img = Image(img_buffer, width=5*inch, height=3.3*inch)
                story.append(img)
                
                # Add legend
                story.append(Spacer(1, 12))
                legend_text = """
                <b>Legend:</b><br/>
                • <font color="red">Red points</font>: Significantly upregulated genes (p &lt; {pval}, log2FC ≥ {fc})<br/>
                • <font color="blue">Blue points</font>: Significantly downregulated genes (p &lt; {pval}, log2FC ≤ -{fc})<br/>
                • <font color="green">Green points</font>: Significant but low fold change<br/>
                • <font color="gray">Gray points</font>: Non-significant genes<br/>
                • Horizontal dashed line: p-value cutoff ({pval})<br/>
                • Vertical dotted lines: log2FC thresholds (±{fc})
                """.format(pval=report_data.pval_cutoff, fc=report_data.logfc_threshold)
                story.append(Paragraph(legend_text, normal_style))
                
                logger.info("Added volcano plot to PDF report")
                
            except Exception as e:
                logger.warning(f"Failed to add volcano plot to PDF: {e}")
                story.append(Spacer(1, 30))
                story.append(Paragraph("Volcano Plot", heading_style))
                story.append(Paragraph("Error generating volcano plot visualization.", normal_style))
        
        elif report_data.volcano_data and not PLOTLY_AVAILABLE:
            story.append(Spacer(1, 30))
            story.append(Paragraph("Volcano Plot", heading_style))
            story.append(Paragraph("Volcano plot not available - Plotly not installed.", normal_style))
        
        # Network Visualization Section
        logger.info(f"ReportLab: network_png present: {bool(report_data.network_png)}")
        if report_data.network_png:
            logger.info(f"ReportLab: network_png length: {len(report_data.network_png)}")
            logger.info(f"ReportLab: network_png starts with: {report_data.network_png[:50]}...")
            
            try:
                story.append(Spacer(1, 30))
                story.append(Paragraph("AOP Network Visualization", heading_style))
                story.append(Paragraph("Adverse Outcome Pathway showing Key Events and their relationships", normal_style))
                story.append(Spacer(1, 12))
                
                # Handle network image data
                network_data = report_data.network_png
                
                if network_data.startswith('<svg') and SVG_SUPPORT:
                    # Handle SVG data directly
                    from io import StringIO
                    svg_buffer = StringIO(network_data)
                    drawing = renderSVG(svg_buffer)
                    if drawing:
                        drawing.width = 6*inch
                        drawing.height = 4*inch
                        story.append(drawing)
                    logger.info("Added SVG network visualization to PDF")
                else:
                    # Handle base64 encoded image data (PNG/JPEG)
                    import base64
                    # Remove data URL prefix if present
                    if network_data.startswith('data:image/'):
                        network_data = network_data.split(',')[1]
                    
                    # Create image from base64 data
                    img_bytes = base64.b64decode(network_data)
                    img_buffer = BytesIO(img_bytes)
                    
                    # Get original image dimensions to preserve aspect ratio
                    from PIL import Image as PILImage
                    pil_img = PILImage.open(BytesIO(img_bytes))
                    original_width, original_height = pil_img.size
                    
                    # Calculate dimensions preserving aspect ratio, max width 7 inches
                    max_width = 7*inch
                    aspect_ratio = original_height / original_width
                    img_width = min(max_width, 7*inch)
                    img_height = img_width * aspect_ratio
                    
                    img_buffer.seek(0)  # Reset buffer position
                    img = Image(img_buffer, width=img_width, height=img_height)
                    story.append(img)
                    logger.info(f"Added PNG network visualization to PDF ({original_width}x{original_height} -> {img_width:.0f}x{img_height:.0f})")
                
                # Add network legend
                story.append(Spacer(1, 12))
                legend_text = """
                <b>Network Legend:</b><br/>
                • <font color="#b3e6b3">Green circles</font>: Molecular Initiating Events (MIEs)<br/>
                • <font color="#ffd9b3">Orange circles</font>: Intermediate Key Events<br/>
                • <font color="#f4b3b3">Pink circles</font>: Adverse Outcomes (AOs)<br/>
                • <font color="red">Red border</font>: Significantly enriched Key Events (FDR &lt; 0.05)<br/>
                • Arrows indicate causal relationships between Key Events
                """
                story.append(Paragraph(legend_text, normal_style))
                
                # Add Network Statistics
                story.append(Spacer(1, 20))
                story.append(Paragraph("📊 Network Statistics", heading_style))
                
                # Extract statistics from network data if available
                if report_data.network_data and 'nodes' in report_data.network_data:
                    nodes = report_data.network_data['nodes']
                    edges = report_data.network_data['edges']
                    
                    # Count different node types
                    total_kes = len([n for n in nodes if not n.get('classes', '').startswith('gene')])
                    total_genes = len([n for n in nodes if n.get('classes', '').startswith('gene')])
                    sig_kes = len([n for n in nodes if 'significant' in n.get('classes', '') and not n.get('classes', '').startswith('gene')])
                    sig_genes = len([n for n in nodes if 'significant' in n.get('classes', '') and n.get('classes', '').startswith('gene')])
                    ke_edges = len([e for e in edges if not e.get('classes', '').startswith('gene-link')])
                    gene_edges = len([e for e in edges if e.get('classes', '').startswith('gene-link')])
                    
                    # Create statistics table
                    stats_data = [
                        ['Metric', 'Value'],
                        ['Total Key Events', str(total_kes)],
                        ['Significant Key Events', str(sig_kes)],
                        ['Total Gene Nodes', str(total_genes)],
                        ['Significant Gene Nodes', str(sig_genes)],
                        ['KE-KE Relationships', str(ke_edges)],
                        ['KE-Gene Connections', str(gene_edges)],
                        ['Total Network Edges', str(len(edges))]
                    ]
                    
                    stats_table = Table(stats_data, colWidths=[3*inch, 1.5*inch])
                    stats_table.setStyle(TableStyle([
                        # Header styling
                        ('BACKGROUND', (0, 0), (-1, 0), HexColor('#307BBF')),
                        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                        ('ALIGN', (1, 0), (1, -1), 'RIGHT'),  # Right-align values
                        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                        ('FONTSIZE', (0, 0), (-1, 0), 12),
                        
                        # Data styling
                        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
                        ('FONTSIZE', (0, 1), (-1, -1), 10),
                        ('GRID', (0, 0), (-1, -1), 1, colors.black),
                        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                        
                        # Alternating row colors
                        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, HexColor('#f9f9f9')])
                    ]))
                    
                    story.append(stats_table)
                    logger.info("Added network statistics to PDF report")
                
                logger.info("Added network visualization to PDF report")
                
            except Exception as e:
                logger.warning(f"Failed to add network visualization to PDF: {e}")
                story.append(Spacer(1, 30))
                story.append(Paragraph("AOP Network Visualization", heading_style))
                story.append(Paragraph("Error generating network visualization.", normal_style))
        
        # System Information
        story.append(Spacer(1, 30))
        story.append(Paragraph("System Information", heading_style))
        
        sys_info = f"""
        <b>Analysis Date:</b> {report_data.analysis_timestamp.strftime('%Y-%m-%d %H:%M:%S')}<br/>
        <b>Application:</b> Molecular AOP Analyser<br/>
        <b>Report Format:</b> PDF (ReportLab)
        """
        
        if report_data.software_versions:
            sys_info += "<br/><b>Software Versions:</b><br/>"
            for package, version in report_data.software_versions.items():
                sys_info += f"• {package}: {version}<br/>"
        
        story.append(Paragraph(sys_info, normal_style))
        
        # Build PDF
        doc.build(story)
        
        # Get PDF content
        buffer.seek(0)
        pdf_bytes = buffer.getvalue()
        buffer.close()
        
        logger.info(f"ReportLab PDF generated successfully ({len(pdf_bytes)} bytes)")
        return pdf_bytes


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
