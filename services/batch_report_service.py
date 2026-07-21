"""
Integrated batch report generation for the Molecular AOP Analyser.

Builds a single PDF or HTML document for a completed batch combining:
- batch metadata and software versions,
- a server-rendered cross-condition comparison heatmap,
- the wide comparison matrix table,
- one section per condition (metadata + enrichment table + KE network image).

Networks and the heatmap are rendered server-side (matplotlib/networkx and
plotly+kaleido respectively) so no browser is required. All image rendering is
wrapped in try/except with a graceful "image unavailable" fallback.
"""

import base64
import json
import logging
from datetime import datetime
from io import BytesIO
from typing import Any, Dict, List, Optional

from services.report_service import (
    report_generator,
    ReportData,
    get_software_versions,
    REPORTLAB_AVAILABLE,
)
from services.comparison_service import comparison_matrix_to_dataframe
from services.enrichment_service import format_ke_summary
from services.network_service import ke_accounting_from_network

logger = logging.getLogger(__name__)

# VHP4Safety-aligned KE node colours by type (matches the on-screen legend).
_KE_TYPE_COLORS = {
    'MIE': '#b3e6b3',           # green
    'molecular initiating event': '#b3e6b3',
    'intermediate': '#ffd9b3',  # orange
    'KE': '#ffd9b3',
    'AO': '#f4b3b3',            # pink
    'adverse outcome': '#f4b3b3',
}


def render_heatmap_png(comparison_data: dict) -> Optional[bytes]:
    """Render the cross-condition comparison heatmap as PNG bytes.

    Uses the -log10(FDR) matrix with the same yellow→red colourscale as the
    interactive compare page. Returns None on any failure (missing data,
    plotly/kaleido unavailable).

    Args:
        comparison_data: Dict from build_comparison_matrix.

    Returns:
        PNG bytes, or None if rendering is not possible.
    """
    if not comparison_data or not comparison_data.get('ke_labels'):
        return None
    try:
        import plotly.graph_objects as go

        z = comparison_data['neg_log10_matrix']
        x = comparison_data['condition_labels']
        y = comparison_data['ke_titles']

        fig = go.Figure(data=go.Heatmap(
            z=z, x=x, y=y,
            colorscale=[
                [0, '#ffffcc'], [0.25, '#feb24c'], [0.5, '#fd8d3c'],
                [0.75, '#e31a1c'], [1, '#800026'],
            ],
            colorbar=dict(title='-log10(FDR)'),
            hoverongaps=False,
        ))
        # Size scales with row/column count, bounded for sane PDFs.
        height = max(300, min(1200, 40 * len(y) + 120))
        width = max(500, min(1400, 120 * len(x) + 260))
        fig.update_layout(
            margin=dict(l=240, r=60, t=30, b=120),
            yaxis=dict(autorange='reversed'),
            xaxis=dict(tickangle=-30),
            paper_bgcolor='white', plot_bgcolor='#f0f0f0',
            height=height, width=width,
        )
        return fig.to_image(format='png', scale=2)
    except Exception as exc:
        logger.warning('render_heatmap_png failed: %s', exc)
        return None


def render_ke_network_png(network_json: Any) -> Optional[bytes]:
    """Render a KE-only AOP network skeleton as PNG bytes.

    Converts a Cytoscape network (KE nodes + KER edges) to a networkx DiGraph
    and draws it with matplotlib: node colour by KE type, red border for
    significant KEs. Gene nodes (if any) are excluded. Returns None on failure.

    Args:
        network_json: Parsed Cytoscape network dict, or a JSON string.

    Returns:
        PNG bytes, or None if rendering is not possible.
    """
    try:
        if isinstance(network_json, str):
            network_json = json.loads(network_json)
        if not network_json:
            return None

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import networkx as nx

        nodes = network_json.get('nodes', [])
        edges = network_json.get('edges', [])

        g = nx.DiGraph()
        node_colors = []
        edge_colors = []
        labels = {}

        for node in nodes:
            data = node.get('data', {})
            classes = node.get('classes', '') or ''
            # Skip gene-level nodes — KE skeleton only.
            if classes.startswith('gene') or data.get('ke_type') is None:
                continue
            ke_id = data.get('id')
            if not ke_id:
                continue
            g.add_node(ke_id)
            raw_label = data.get('label') or ke_id
            labels[ke_id] = raw_label if len(raw_label) <= 22 else raw_label[:19] + '…'
            node_colors.append(_KE_TYPE_COLORS.get(data.get('ke_type', 'intermediate'), '#ffd9b3'))
            edge_colors.append('#cc0000' if 'significant' in classes else '#999999')

        if g.number_of_nodes() == 0:
            return None

        for edge in edges:
            data = edge.get('data', {})
            src, tgt = data.get('source'), data.get('target')
            if src in g and tgt in g:
                g.add_edge(src, tgt)

        # Layout: prefer graphviz 'dot' (left-to-right DAG); fall back to spring.
        pos = None
        try:
            pos = nx.nx_agraph.graphviz_layout(g, prog='dot')
        except Exception:
            try:
                pos = nx.nx_pydot.graphviz_layout(g, prog='dot')
            except Exception:
                pos = nx.spring_layout(g, seed=42, k=1.5)

        n = g.number_of_nodes()
        fig_w = max(6, min(14, n * 1.1))
        fig_h = max(4, min(10, n * 0.7))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        nx.draw_networkx_edges(g, pos, ax=ax, edge_color='#aaaaaa',
                               arrows=True, arrowsize=12, width=1.2)
        nx.draw_networkx_nodes(g, pos, ax=ax, node_color=node_colors,
                               edgecolors=edge_colors, linewidths=2.0, node_size=1400)
        nx.draw_networkx_labels(g, pos, labels=labels, ax=ax, font_size=7)
        ax.axis('off')
        fig.tight_layout()

        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()
    except Exception as exc:
        logger.warning('render_ke_network_png failed: %s', exc)
        return None


def _png_to_data_uri(png_bytes: bytes) -> str:
    """Encode PNG bytes as a base64 data URI for HTML embedding."""
    b64 = base64.b64encode(png_bytes).decode('ascii')
    return f'data:image/png;base64,{b64}'


def _batch_meta_rows(batch, conditions, comparison_data) -> List[tuple]:
    """Build (label, value) rows describing the batch for the report header."""
    completed = batch.completed_at.strftime('%Y-%m-%d %H:%M') if batch.completed_at else 'N/A'
    return [
        ('AOP', f'{batch.aop_label or ""} ({batch.aop_id or "N/A"})'),
        ('Conditions', str(len(conditions))),
        ('Gene Set Resources', batch.selected_resources or 'WikiPathways'),
        ('Log2FC Threshold', str(batch.logfc_threshold)),
        ('P-value Cutoff', str(batch.pval_cutoff)),
        ('Harmonised Background', f'{batch.harmonised_gene_count or 0:,} genes'),
        ('Owner', batch.owner or 'N/A'),
        ('Completed', completed),
    ]


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

def generate_batch_html(batch, conditions, comparison_data: dict) -> str:
    """Generate the integrated batch report as a standalone HTML string.

    Args:
        batch: BatchRecord ORM object.
        conditions: List of completed ConditionRecord objects (position order).
        comparison_data: Dict from build_comparison_matrix.

    Returns:
        Complete HTML document string.
    """
    css = report_generator._get_report_css()
    versions = get_software_versions()
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Batch metadata table.
    meta_rows = ''.join(
        f'<div class="metadata-item"><label>{label}:</label><span>{value}</span></div>'
        for label, value in _batch_meta_rows(batch, conditions, comparison_data)
    )

    # Comparison heatmap.
    heatmap_png = render_heatmap_png(comparison_data)
    if heatmap_png:
        heatmap_html = (
            '<div class="plot-container">'
            f'<img src="{_png_to_data_uri(heatmap_png)}" style="max-width:100%;" '
            'alt="Comparison heatmap"/></div>'
        )
    else:
        heatmap_html = '<p class="note">Heatmap image unavailable.</p>'

    # Comparison matrix table (FDR, top 30 KEs).
    df = comparison_matrix_to_dataframe(comparison_data, which='fdr')
    if not df.empty:
        comp_table_html = df.head(30).to_html(index=False, na_rep='', border=0,
                                              classes='enrichment-table')
    else:
        comp_table_html = '<p class="note">No comparison data available.</p>'

    # Per-condition sections.
    condition_sections = []
    for cond in conditions:
        try:
            enrichment = json.loads(cond.enrichment_json) if cond.enrichment_json else []
        except (json.JSONDecodeError, TypeError):
            enrichment = []

        rd = _condition_report_data(batch, cond, enrichment)
        enrichment_html = report_generator._generate_enrichment_section(rd)
        # Issue #65: state the multiple-testing denominator for this condition.
        ke_accounting = format_ke_summary(ke_accounting_from_network(cond.network_json))

        net_png = render_ke_network_png(cond.network_json)
        if net_png:
            net_html = (
                '<div class="plot-container">'
                f'<img src="{_png_to_data_uri(net_png)}" style="max-width:100%;" '
                'alt="KE network"/></div>'
            )
        else:
            net_html = '<p class="note">Network image unavailable.</p>'

        meta_bits = []
        if cond.dose:
            meta_bits.append(f'Dose: {cond.dose}')
        if cond.timepoint:
            meta_bits.append(f'Time: {cond.timepoint}')
        meta_line = ' · '.join(meta_bits)

        condition_sections.append(f"""
        <section class="condition-section">
            <h2>Condition: {cond.condition_label}</h2>
            <p class="section-description">
                {meta_line}{' — ' if meta_line else ''}
                {cond.significant_genes or 0:,} significant of {cond.gene_count or 0:,} genes
            </p>
            {f'<p class="note">Key Events: {ke_accounting}</p>' if ke_accounting else ''}
            {enrichment_html}
            <h3>AOP Network</h3>
            {net_html}
        </section>
        """)

    version_items = ''.join(f'<li><strong>{p}:</strong> {v}</li>' for p, v in versions.items())

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Molecular AOP Batch Report</title>
    <style>{css}</style>
</head>
<body>
    <header class="report-header">
        <h1>Molecular AOP Batch Report</h1>
        <div class="report-meta">
            <p><strong>Batch:</strong> {batch.batch_name or 'Batch Analysis'}</p>
            <p><strong>Generated:</strong> {ts}</p>
        </div>
    </header>
    <div class="container">
        <section class="metadata-section">
            <h2>Batch Information</h2>
            <div class="metadata-grid">{meta_rows}</div>
            {f'<div class="metadata-item full-width"><label>Description:</label><p>{batch.description}</p></div>' if batch.description else ''}
        </section>
        <section class="comparison-section">
            <h2>Cross-Condition Comparison</h2>
            <p class="section-description">Significance (-log10 FDR) of each Key Event across conditions.</p>
            {heatmap_html}
            <div class="table-container">{comp_table_html}</div>
        </section>
        {''.join(condition_sections)}
        <section class="system-info">
            <h2>System Information</h2>
            <div class="info-grid">
                <div class="info-item"><label>Generated:</label><span>{ts}</span></div>
                <div class="info-item"><label>Application:</label><span>Molecular AOP Analyser</span></div>
            </div>
            {f'<div class="versions"><h3>Software Versions</h3><ul>{version_items}</ul></div>' if version_items else ''}
        </section>
    </div>
</body>
</html>"""


def _condition_report_data(batch, cond, enrichment: List[Dict]) -> ReportData:
    """Construct a minimal ReportData for one condition (for section reuse)."""
    return ReportData(
        metadata={'dataset_id': cond.condition_label},
        filename=cond.filename or cond.condition_label,
        gene_count=cond.gene_count or 0,
        significant_genes=cond.significant_genes or 0,
        aop_id=batch.aop_id or '',
        aop_label=batch.aop_label or '',
        logfc_threshold=batch.logfc_threshold,
        pval_cutoff=batch.pval_cutoff,
        id_column=batch.id_column or '',
        fc_column=batch.fc_column or '',
        pval_column=batch.pval_column or '',
        id_type='symbol',
        enrichment_results=enrichment,
        method='ora',
        selected_resources=batch.selected_resources or 'WikiPathways',
    )


# ---------------------------------------------------------------------------
# PDF report
# ---------------------------------------------------------------------------

def generate_batch_pdf(batch, conditions, comparison_data: dict) -> bytes:
    """Generate the integrated batch report as PDF bytes (ReportLab).

    Args:
        batch: BatchRecord ORM object.
        conditions: List of completed ConditionRecord objects (position order).
        comparison_data: Dict from build_comparison_matrix.

    Returns:
        PDF content as bytes.

    Raises:
        RuntimeError: If ReportLab is unavailable.
    """
    if not REPORTLAB_AVAILABLE:
        raise RuntimeError('ReportLab not available for PDF generation.')

    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import (
        SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak, Image,
    )
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.colors import HexColor, whitesmoke, black, lightyellow
    from reportlab.lib import colors as rl_colors
    from reportlab.lib.units import inch

    primary = HexColor('#307BBF')
    dark = HexColor('#29235C')

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle('BatchTitle', parent=styles['Heading1'], fontSize=22,
                                 spaceAfter=24, textColor=dark, alignment=1)
    heading_style = ParagraphStyle('BatchHeading', parent=styles['Heading2'], fontSize=15,
                                   spaceAfter=10, textColor=primary)
    normal_style = ParagraphStyle('BatchNormal', parent=styles['Normal'], fontSize=10, textColor=dark)

    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4, rightMargin=54, leftMargin=54,
                            topMargin=54, bottomMargin=36)
    story = []

    # Title + batch metadata.
    story.append(Paragraph('Molecular AOP Batch Report', title_style))
    story.append(Paragraph(f"<b>Batch:</b> {batch.batch_name or 'Batch Analysis'}", normal_style))
    story.append(Spacer(1, 12))

    meta_data = [['Field', 'Value']] + [list(r) for r in _batch_meta_rows(batch, conditions, comparison_data)]
    meta_table = Table(meta_data, colWidths=[2.2 * inch, 4.0 * inch])
    meta_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), primary),
        ('TEXTCOLOR', (0, 0), (-1, 0), whitesmoke),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 9),
        ('GRID', (0, 0), (-1, -1), 0.5, black),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [rl_colors.white, HexColor('#f4f6fb')]),
    ]))
    story.append(meta_table)
    story.append(Spacer(1, 20))

    # Comparison heatmap.
    story.append(Paragraph('Cross-Condition Comparison', heading_style))
    heatmap_png = render_heatmap_png(comparison_data)
    if heatmap_png:
        try:
            from PIL import Image as PILImage
            pil = PILImage.open(BytesIO(heatmap_png))
            ow, oh = pil.size
            max_w = 6.8 * inch
            img_w = min(max_w, 6.8 * inch)
            img_h = img_w * (oh / ow)
            story.append(Image(BytesIO(heatmap_png), width=img_w, height=img_h))
        except Exception as exc:
            logger.warning('batch pdf: heatmap embed failed: %s', exc)
            story.append(Paragraph('Heatmap image unavailable.', normal_style))
    else:
        story.append(Paragraph('Heatmap image unavailable.', normal_style))
    story.append(Spacer(1, 16))

    # Comparison matrix table (FDR, top 25).
    df = comparison_matrix_to_dataframe(comparison_data, which='fdr')
    if not df.empty:
        df = df.head(25)
        # Format floats compactly; blank for NaN.
        def _fmt(v):
            if v is None or (isinstance(v, float) and v != v):
                return ''
            if isinstance(v, float):
                return f'{v:.2e}' if v < 0.001 else f'{v:.4f}'
            return str(v)
        header = list(df.columns)
        comp_rows = [header] + [[_fmt(v) for v in row] for row in df.itertuples(index=False)]
        n_cols = len(header)
        comp_table = Table(comp_rows, repeatRows=1)
        comp_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), primary),
            ('TEXTCOLOR', (0, 0), (-1, 0), whitesmoke),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 6.5),
            ('GRID', (0, 0), (-1, -1), 0.4, black),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ]))
        story.append(comp_table)
    story.append(PageBreak())

    # Per-condition sections.
    for idx, cond in enumerate(conditions):
        try:
            enrichment = json.loads(cond.enrichment_json) if cond.enrichment_json else []
        except (json.JSONDecodeError, TypeError):
            enrichment = []

        meta_bits = []
        if cond.dose:
            meta_bits.append(f'Dose: {cond.dose}')
        if cond.timepoint:
            meta_bits.append(f'Time: {cond.timepoint}')
        meta_line = ' · '.join(meta_bits)

        story.append(Paragraph(f'Condition: {cond.condition_label}', heading_style))
        story.append(Paragraph(
            f"{meta_line}{' — ' if meta_line else ''}"
            f"{cond.significant_genes or 0:,} significant of {cond.gene_count or 0:,} genes",
            normal_style,
        ))
        # Issue #65: state the multiple-testing denominator for this condition.
        ke_accounting = format_ke_summary(ke_accounting_from_network(cond.network_json))
        if ke_accounting:
            story.append(Paragraph(f'Key Events: {ke_accounting}', normal_style))
        story.append(Spacer(1, 8))

        if enrichment:
            table = report_generator.build_enrichment_table_flowable(enrichment, is_gsea=False, max_rows=15)
            if table is not None:
                story.append(table)
        else:
            story.append(Paragraph('No enrichment results for this condition.', normal_style))
        story.append(Spacer(1, 12))

        net_png = render_ke_network_png(cond.network_json)
        if net_png:
            try:
                from PIL import Image as PILImage
                pil = PILImage.open(BytesIO(net_png))
                ow, oh = pil.size
                img_w = min(6.8 * inch, 6.8 * inch)
                img_h = img_w * (oh / ow)
                story.append(Paragraph('AOP Network', normal_style))
                story.append(Image(BytesIO(net_png), width=img_w, height=img_h))
            except Exception as exc:
                logger.warning('batch pdf: network embed failed for %s: %s', cond.condition_label, exc)
                story.append(Paragraph('Network image unavailable.', normal_style))
        else:
            story.append(Paragraph('Network image unavailable.', normal_style))

        # Page break between conditions (not after the last).
        if idx < len(conditions) - 1:
            story.append(PageBreak())

    # System info.
    story.append(Spacer(1, 20))
    story.append(Paragraph('System Information', heading_style))
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    sys_info = f'<b>Generated:</b> {ts}<br/><b>Application:</b> Molecular AOP Analyser<br/>'
    for p, v in get_software_versions().items():
        sys_info += f'• {p}: {v}<br/>'
    story.append(Paragraph(sys_info, normal_style))

    doc.build(story)
    buffer.seek(0)
    pdf_bytes = buffer.getvalue()
    buffer.close()
    logger.info('Batch PDF generated (%d bytes, %d conditions)', len(pdf_bytes), len(conditions))
    return pdf_bytes
