"""Shared template context for ``results.html``.

``results.html`` is rendered from more than one route: the live ``/analyze``
run, and the re-render of a stored batch condition. Both have to supply the
same set of template variables, and until now each route built that list by
hand. The lists drifted: the batch route was missing ``ke_summary`` and
``hub_list`` outright, which crashed the page (both reach ``| tojson``, and
Jinja's ``Undefined`` is not JSON-serialisable), and it silently omitted the
FDR selector, the KE accounting, the resource-provenance line and the
representation badges, because an undefined name is merely falsy everywhere
else.

``build_results_context()`` is the single definition of that contract. A route
supplies what it knows and the builder fills in the rest with the same
defaults, so a new template variable is added in one place rather than in
every caller.

Only variables the route has to vary belong here. ``representation_labels``
and ``resource_source_labels`` are Jinja globals (registered on
``app.jinja_env.globals``) and ``column_names`` is a ``{% set %}`` local in the
template itself; none of them are route state.
"""

import json
from typing import Any, Dict, List, Optional

from config import Config
from services.enrichment_service import format_ke_summary


def _as_json(value: Any) -> str:
    """Serialise to JSON, passing through a value that is already a JSON string.

    The live route pre-serialises some payloads while the batch route holds
    parsed structures; accepting both keeps the callers free of the
    distinction.
    """
    if isinstance(value, str):
        return value
    return json.dumps(value if value is not None else None)


def build_results_context(
    *,
    table: List[Dict[str, Any]],
    network: Any,
    metadata: Dict[str, Any],
    table_json: Optional[str] = None,
    volcano_data: Optional[List[Dict[str, Any]]] = None,
    id_type: str = "symbol",
    background_size: int = 0,
    threshold: float = 0.0,
    ke_gene_map: Any = None,
    ke_type_map: Any = None,
    ke_title_map: Any = None,
    method: str = "ora",
    pval_threshold: Optional[float] = None,
    tour_active: bool = False,
    hub_list: Optional[List[Dict[str, Any]]] = None,
    wp_picker_data: Any = None,
    ke_summary: Optional[Dict[str, Any]] = None,
    background_overlap: Any = None,
    dropped_rows: Optional[int] = None,
) -> Dict[str, Any]:
    """Assemble the full template context for ``results.html``.

    Args:
        table: Enrichment rows as native Python types.
        network: Cytoscape network — parsed dict or a JSON string.
        metadata: Run metadata dict consumed by the template header.
        table_json: Pre-serialised ``table``; derived from ``table`` when omitted.
        volcano_data: Volcano points; an empty list when the route has none
            (the batch re-render does not store them).
        id_type: Gene identifier type reported in the header.
        background_size: Measured-gene count used as the statistical background.
        threshold: log2 fold-change threshold applied to the run.
        ke_gene_map: KE -> gene records; parsed dict or JSON string.
        ke_type_map: KE -> node type; parsed dict or JSON string.
        ke_title_map: KE -> title; parsed dict or JSON string.
        method: ``'ora'`` or ``'gsea'`` — drives the method-specific columns.
        pval_threshold: Per-gene significance threshold, when the route set one.
        tour_active: Whether the guided tour should resume on load.
        hub_list: Ranked hub genes; an empty list renders no hub section.
        wp_picker_data: Pathway-view picker payload, or None to hide it.
        ke_summary: Issue #65 tested/excluded accounting. ``None`` omits the
            accounting line rather than failing.
        background_overlap: Issue #69 ID-column sanity-check result.
        dropped_rows: Issue #103 count of upload rows discarded for carrying no
            usable gene identifier. ``None`` means "not recorded" — a stored run
            that predates the counter — and renders nothing, which is distinct
            from a recorded zero.

    Returns:
        Keyword arguments for ``render_template('results.html', **context)``.
    """
    return {
        "table": table,
        "table_json": table_json if table_json is not None else json.dumps(table),
        "volcano_data": volcano_data or [],
        "volcano_json": json.dumps(volcano_data or []),
        "id_type": id_type,
        "background_size": background_size,
        "threshold": threshold,
        "network_json": _as_json(network),
        "ke_gene_json": _as_json(ke_gene_map if ke_gene_map is not None else {}),
        "ke_type_map": _as_json(ke_type_map if ke_type_map is not None else {}),
        "ke_title_map": _as_json(ke_title_map if ke_title_map is not None else {}),
        "metadata": metadata,
        "method": method or "ora",
        "pval_threshold": pval_threshold,
        "tour_active": tour_active,
        # Both of these reach `| tojson` in the template, so they must be
        # concrete values rather than an absent name.
        "hub_list": hub_list or [],
        "ke_summary": ke_summary,
        "ke_summary_text": format_ke_summary(ke_summary),
        "wp_picker_data": wp_picker_data,
        "fdr_cutoff": Config.SIGNIFICANCE_FDR_CUTOFF,
        "fdr_choices": Config.SIGNIFICANCE_FDR_CHOICES,
        "background_overlap": background_overlap,
        "dropped_rows": dropped_rows,
    }
