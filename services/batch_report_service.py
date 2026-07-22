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
from services.comparison_service import (
    comparison_matrix_to_dataframe,
    comparison_nes_display_dataframe,
)
from services.image_render import render_figure_png
from services.enrichment_service import format_ke_summary
from services.gsea_service import NES_BEYOND_RESOLUTION
from services.network_service import ke_accounting_from_network
from config import Config
from helpers import MIN_CONFIDENCE_LABELS

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


def batch_method(batch) -> str:
    """Return the enrichment method a batch was run with ('ora' or 'gsea').

    Issue #76. Reads ``BatchRecord.effective_method()`` where available so a
    batch created before the method column existed reads as ORA, and tolerates
    plain stand-in objects that only carry a ``method`` attribute.

    Args:
        batch: BatchRecord (or any object exposing ``effective_method`` or
            ``method``).

    Returns:
        str: 'gsea' when the batch was run with GSEA, otherwise 'ora'.
    """
    getter = getattr(batch, 'effective_method', None)
    if callable(getter):
        try:
            return getter() or 'ora'
        except Exception:  # pragma: no cover — provenance must never break a report
            pass
    return (getattr(batch, 'method', None) or 'ora')


def _mask_nes_by_significance(nes_matrix, fdr_matrix):
    """Blank NES cells that are not significant at the active cutoff (#107).

    On a diverging scale colour encodes direction, and for a cell nowhere near
    significance the sign of NES is not a stable quantity: the same conditions
    scored against a harmonised batch background and against their own
    per-file backgrounds returned opposite-signed NES for four of six KEs in
    the weakest condition, every one with FDR >= 0.48. Colouring those cells
    invites the reader to interpret a direction the data does not support, so
    they are left blank; the values remain in the matrix table and the exports.

    Args:
        nes_matrix: 2-D list of NES floats (rows=KEs, cols=conditions).
        fdr_matrix: the aligned FDR matrix. When absent, nothing can be
            assessed and the NES matrix is returned unchanged.

    Returns:
        list: a new matrix with non-significant cells set to None.
    """
    if not fdr_matrix:
        return nes_matrix
    cutoff = Config.SIGNIFICANCE_FDR_CUTOFF
    masked = []
    for row_idx, row in enumerate(nes_matrix):
        fdr_row = fdr_matrix[row_idx] if row_idx < len(fdr_matrix) else []
        masked.append([
            value if (col_idx < len(fdr_row)
                      and fdr_row[col_idx] is not None
                      and fdr_row[col_idx] < cutoff)
            else None
            for col_idx, value in enumerate(row)
        ])
    return masked


def _beyond_resolution_cells(comparison_data: dict):
    """Locate the cells whose NES could not be normalised (#117).

    Args:
        comparison_data: Dict from ``build_comparison_matrix``. Matrices stored
            before #117 carry no ``nes_status_matrix``, which yields no cells —
            those runs did not know, and marking nothing is the honest result.

    Returns:
        tuple: aligned ``(condition labels, KE titles)`` lists, one entry per
        beyond-resolution cell, ready for a Plotly scatter overlay.
    """
    status = (comparison_data or {}).get('nes_status_matrix') or []
    x_labels = (comparison_data or {}).get('condition_labels') or []
    y_labels = (comparison_data or {}).get('ke_titles') or []
    xs, ys = [], []
    for row_idx, row in enumerate(status):
        if row_idx >= len(y_labels) or not row:
            continue
        for col_idx, value in enumerate(row):
            if value != NES_BEYOND_RESOLUTION or col_idx >= len(x_labels):
                continue
            xs.append(x_labels[col_idx])
            ys.append(y_labels[row_idx])
    return xs, ys


def _dropped_rows_note(cond) -> str:
    """Describe the rows a condition's file lost to a missing gene symbol (#103).

    Args:
        cond: ConditionRecord. ``dropped_unidentified_rows`` is NULL for
            conditions run before the count was persisted.

    Returns:
        str: a leading-separator clause ready to append to the condition's
        gene-count line, or '' when nothing was dropped or nothing was
        recorded — an unrecorded count must not read as "none dropped".
    """
    dropped = getattr(cond, 'dropped_unidentified_rows', None)
    if not dropped:
        return ''
    noun = 'row' if dropped == 1 else 'rows'
    verb = 'was' if dropped == 1 else 'were'
    return f' · {dropped:,} {noun} without a gene symbol {verb} excluded'


def render_heatmap_png(comparison_data: dict) -> Optional[bytes]:
    """Render the cross-condition comparison heatmap as PNG bytes.

    Mirrors the interactive compare page: the -log10(FDR) matrix on a
    yellow→red scale for an ORA batch, and the signed NES on a diverging
    red/blue scale clamped to ±3 for a GSEA batch (issue #76), so the printed
    heatmap shows the same quantity the screen does. Returns None on any
    failure (missing data, plotly/kaleido unavailable).

    Args:
        comparison_data: Dict from build_comparison_matrix. Its ``method`` key
            selects the matrix; a matrix dict written before #76 has no such
            key and renders as ORA.

    Returns:
        PNG bytes, or None if rendering is not possible.
    """
    if not comparison_data or not comparison_data.get('ke_labels'):
        return None
    try:
        import plotly.graph_objects as go

        is_gsea = comparison_data.get('method') == 'gsea'
        x = comparison_data['condition_labels']
        y = comparison_data['ke_titles']

        if is_gsea:
            z = comparison_data.get('nes_matrix') or comparison_data['neg_log10_matrix']
            # Issue #107: colour direction only where direction is supported.
            # Applied here as well as on the compare page so the printed
            # heatmap and the screen cannot make different claims.
            if comparison_data.get('nes_matrix'):
                z = _mask_nes_by_significance(z, comparison_data.get('fdr_matrix'))
            heatmap_kwargs = dict(
                colorscale='RdBu',
                reversescale=True,  # red = positive NES (coordinated up-shift)
                zmin=-3, zmax=3, zmid=0,
                colorbar=dict(title='NES'),
            )
        else:
            z = comparison_data['neg_log10_matrix']
            heatmap_kwargs = dict(
                colorscale=[
                    [0, '#ffffcc'], [0.25, '#feb24c'], [0.5, '#fd8d3c'],
                    [0.75, '#e31a1c'], [1, '#800026'],
                ],
                colorbar=dict(title='-log10(FDR)'),
            )

        # Issue #107: the significance mask can blank every cell, and plotly
        # infers its axes from the data — an all-null z renders a featureless
        # box with no tick labels at all. A fully transparent companion trace
        # over the same categories anchors the axes without drawing anything,
        # so a blank heatmap still says which KEs and conditions it covers.
        axis_anchor = go.Heatmap(
            z=[[0] * len(x) for _ in y], x=x, y=y,
            opacity=0, showscale=False, hoverinfo='skip',
        )
        traces = [
            axis_anchor,
            go.Heatmap(z=z, x=x, y=y, hoverongaps=False, **heatmap_kwargs),
        ]
        # Issue #117: a Key Event that beat every permutation carries no NES, so
        # the mask above cannot colour it and it would print as a blank cell —
        # visually identical to a Key Event that was never tested, which is the
        # inverse of the result. Ring it instead of colouring it: there is no
        # magnitude to put on the scale, only a direction, and that is in the
        # legend entry rather than in a colour that would imply one.
        beyond_x, beyond_y = _beyond_resolution_cells(comparison_data)
        if is_gsea and beyond_x:
            traces.append(go.Scatter(
                x=beyond_x, y=beyond_y, mode='markers',
                marker=dict(symbol='circle-open', size=12,
                            line=dict(width=3, color='#111111'),
                            color='#111111'),
                name='beyond permutation resolution',
                showlegend=True, hoverinfo='skip',
            ))
        fig = go.Figure(data=traces)
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
        # Issue #104: bounded render — a wedged kaleido subprocess raises here
        # and the caller falls back to "image unavailable" instead of hanging.
        return render_figure_png(fig, format='png', scale=2)
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


def _resolution_text(batch) -> str:
    """Describe how the batch's gene-set resources resolved (issue #68).

    Imported lazily from ``app`` to avoid a circular import at module load.

    Args:
        batch: BatchRecord whose ``resource_resolution`` JSON is read.

    Returns:
        str: one-line provenance summary, or '' for batches written before the
        resolution was recorded.
    """
    try:
        from app import describe_resource_resolution, _parse_resource_resolution
        return describe_resource_resolution(
            _parse_resource_resolution(getattr(batch, 'resource_resolution', None))
        )
    except Exception as exc:  # pragma: no cover — provenance must never break a report
        logger.warning("Could not describe resource resolution: %s", exc)
        return ''


def _batch_resource_warnings(batch) -> List[str]:
    """Warnings about resources skipped or served from bundled files (#67, #68).

    Args:
        batch: BatchRecord whose ``resource_resolution`` JSON is read.

    Returns:
        list[str]: sentences to print alongside the batch metadata; empty when
        the run resolved exactly as requested (or predates the recording).
    """
    try:
        from app import resource_resolution_warnings, _parse_resource_resolution
        return resource_resolution_warnings(
            _parse_resource_resolution(getattr(batch, 'resource_resolution', None)),
            getattr(batch, 'min_confidence', None) or 'all',
        )
    except Exception as exc:  # pragma: no cover
        logger.warning("Could not build resource warnings: %s", exc)
        return []


def _batch_meta_rows(batch, conditions, comparison_data) -> List[tuple]:
    """Build (label, value) rows describing the batch for the report header."""
    completed = batch.completed_at.strftime('%Y-%m-%d %H:%M') if batch.completed_at else 'N/A'
    is_gsea = batch_method(batch) == 'gsea'
    return [
        ('AOP', f'{batch.aop_label or ""} ({batch.aop_id or "N/A"})'),
        ('Conditions', str(len(conditions))),
        # Issue #76: which statistic every number in this report came from.
        ('Statistical Method',
         'GSEA (rank-based)' if is_gsea else "Fisher's exact (over-representation)"),
        ('Gene Set Resources (requested)', batch.selected_resources or 'WikiPathways'),
        # Issue #68: what the resources actually resolved to. The requested list
        # cannot show a resource that was skipped or served from bundled files.
        ('Gene Set Provenance (used)', _resolution_text(batch) or 'not recorded'),
        # Issue #60: the mapping-confidence threshold the gene sets were built at.
        ('Min. Mapping Confidence',
         MIN_CONFIDENCE_LABELS.get(getattr(batch, 'min_confidence', None) or 'all', 'All mappings')),
        # Issue #76: under GSEA these thresholds do not enter the enrichment
        # statistic — they only define the per-condition significant-gene counts
        # and the driver genes on the comparison page. Say so rather than let
        # them read as the test's cutoffs.
        ('Log2FC Threshold',
         f'{batch.logfc_threshold} (gene counts only; GSEA ranks all genes)'
         if is_gsea else str(batch.logfc_threshold)),
        ('P-value Cutoff',
         f'{batch.pval_cutoff} (gene counts only; GSEA ranks all genes)'
         if is_gsea else str(batch.pval_cutoff)),
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

    # Issues #67 / #68 — where the run differed from the request.
    warning_html = ''.join(
        f'<p class="resource-warning"><strong>Note:</strong> {w}</p>'
        for w in _batch_resource_warnings(batch)
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

    # Comparison matrix table — NES for a GSEA batch, FDR otherwise (#76).
    is_gsea = batch_method(batch) == 'gsea'
    comparison_caption = (
        'Normalised enrichment score (NES, signed) of each Key Event across '
        'conditions. Red is a coordinated up-shift, blue a down-shift. The '
        'heatmap colours only Key Events reaching FDR &lt; '
        f'{Config.SIGNIFICANCE_FDR_CUTOFF} — above that cutoff the sign of NES '
        'is not stable, so colouring it would show noise as direction (#107). '
        # Issue #117 — the previous wording ("every NES value, significant or
        # not") was false: a Key Event that beat every permutation has no NES to
        # carry, and those are the strongest cells in the grid, not the weakest.
        'Not every cell has a NES: where no permutation reached the observed '
        'enrichment score there is nothing to normalise against, and the table '
        'marks those cells &#9650;/&#9660; beyond res. with the direction of '
        'the raw enrichment score (their FDR is its limit of zero). A trailing '
        '! marks a cell normalised against fewer than ten same-signed '
        'permutations — gseapy computes the nominal p-value on that same short '
        'tail, so neither its magnitude nor its p-value is fine-grained enough '
        'to rank on. A trailing ? marks a cell whose permutation null could not '
        'be inspected at all (#117).'
        if is_gsea else
        'Significance (-log10 FDR) of each Key Event across conditions.'
    )
    df = (
        comparison_nes_display_dataframe(comparison_data) if is_gsea
        else comparison_matrix_to_dataframe(comparison_data, which='fdr')
    )
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
                {_dropped_rows_note(cond)}
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
            {warning_html}
        </section>
        <section class="comparison-section">
            <h2>Cross-Condition Comparison</h2>
            <p class="section-description">{comparison_caption}</p>
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
        # Issue #103: NULL on conditions run before the count was persisted.
        # getattr, as elsewhere in this module, so a plain stand-in object with
        # only the fields a caller cares about still produces a report.
        dropped_rows=getattr(cond, 'dropped_unidentified_rows', None),
        aop_id=batch.aop_id or '',
        aop_label=batch.aop_label or '',
        logfc_threshold=batch.logfc_threshold,
        pval_cutoff=batch.pval_cutoff,
        id_column=batch.id_column or '',
        fc_column=batch.fc_column or '',
        pval_column=batch.pval_column or '',
        id_type='symbol',
        enrichment_results=enrichment,
        # Issue #76: the per-condition table is rendered from this, and reading
        # a GSEA row through the ORA columns prints 0 overlap / 0.00 odds ratio
        # and drops NES and the leading-edge genes.
        method=batch_method(batch),
        selected_resources=batch.selected_resources or 'WikiPathways',
        min_confidence=getattr(batch, 'min_confidence', None) or 'all',  # Issue #60
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

    # Issues #67 / #68 — a resource skipped or served from the bundled files
    # changes what every condition in this batch means. Say so in the report.
    for warning in _batch_resource_warnings(batch):
        story.append(Spacer(1, 6))
        story.append(Paragraph(f"<b>Note:</b> {warning}", normal_style))

    story.append(Spacer(1, 20))

    # Comparison heatmap.
    is_gsea = batch_method(batch) == 'gsea'
    story.append(Paragraph('Cross-Condition Comparison', heading_style))
    story.append(Paragraph(
        'Normalised enrichment score (NES, signed) per Key Event across '
        'conditions. The heatmap colours only Key Events reaching FDR &lt; '
        f'{Config.SIGNIFICANCE_FDR_CUTOFF}; above that cutoff the sign of NES is '
        'not stable, so colouring it would show noise as direction (#107). '
        # Issue #117 — see the HTML caption; the same false claim stood here.
        'Not every cell has a NES: where no permutation reached the observed '
        'enrichment score there is nothing to normalise against, and the table '
        'marks those cells "beyond res." with the direction of the raw '
        'enrichment score. A trailing ! marks a cell normalised against fewer '
        'than ten same-signed permutations, whose magnitude and nominal '
        'p-value are both computed on that short tail; a trailing ? marks a '
        'cell whose permutation null could not be inspected (#117).'
        if is_gsea else
        'Significance (-log10 FDR) per Key Event across conditions.',
        normal_style,
    ))
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

    # Comparison matrix table — NES for a GSEA batch, FDR otherwise (#76).
    # Under GSEA the cells are pre-rendered text so a Key Event with no NES is
    # marked rather than blank (#117); _fmt below passes strings through.
    df = (
        comparison_nes_display_dataframe(comparison_data, ascii_only=True)
        if is_gsea
        else comparison_matrix_to_dataframe(comparison_data, which='fdr')
    )
    if not df.empty:
        df = df.head(25)
        # Format floats compactly; blank for NaN. NES is a small signed number,
        # so it is printed with its sign at fixed precision rather than through
        # the FDR formatter, which would render -1.8 in scientific notation.
        def _fmt(v):
            if v is None or (isinstance(v, float) and v != v):
                return ''
            if isinstance(v, float):
                if is_gsea:
                    return f'{v:+.2f}'
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
            f"{cond.significant_genes or 0:,} significant of {cond.gene_count or 0:,} genes"
            f"{_dropped_rows_note(cond)}",
            normal_style,
        ))
        # Issue #65: state the multiple-testing denominator for this condition.
        ke_accounting = format_ke_summary(ke_accounting_from_network(cond.network_json))
        if ke_accounting:
            story.append(Paragraph(f'Key Events: {ke_accounting}', normal_style))
        story.append(Spacer(1, 8))

        if enrichment:
            # Issue #76: GSEA rows have no odds_ratio / num_overlap / Direction,
            # so the ORA columns printed zeros and dropped NES and lead genes.
            table = report_generator.build_enrichment_table_flowable(
                enrichment, is_gsea=is_gsea, max_rows=15
            )
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
