"""
Network visualization service for AOP pathways.
"""
import json
import math
import pandas as pd
import logging
from typing import Dict, Set, List, Any, Optional

from config import Config
from services.enrichment_service import (
    EXCLUDED_ERROR,
    EXCLUDED_NO_MAPPING,
    EXCLUDED_TOO_FEW_GENES,
    EXCLUDED_UNRESOLVED_MAPPING,
)

logger = logging.getLogger(__name__)


def _to_native(value):
    """Convert pandas/numpy values to JSON-safe Python natives.

    - NaN / inf -> None
    - numpy scalars (np.float64, np.int64, ...) -> native Python via .item()
    - None / native types -> passthrough

    Mirrors the sanitisation idiom at app.py:808-815 so values flowing into
    the cytoscape data payload survive json.dumps cleanly (and so the
    .cyjs Download Network export carries no `NaN` literals — see EXPO-06 /
    issue #50).
    """
    if value is None:
        return None
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    # Catches np.float64, np.int64, etc.
    if hasattr(value, 'item'):
        try:
            native = value.item()
            if isinstance(native, float) and (math.isnan(native) or math.isinf(native)):
                return None
            return native
        except (ValueError, AttributeError):
            return value
    return value


def build_cytoscape_network(
    ke_list: Set[str],
    edges: pd.DataFrame,
    enrichment_results: pd.DataFrame,
    ke_title_map: Dict[str, str],
    ke_type_map: Dict[str, str],
    reference_sets: Dict[str, Set[str]] = None,
    method: str = 'ora',
    fdr_cutoff: float = Config.SIGNIFICANCE_FDR_CUTOFF,
    excluded_kes: Optional[Dict[str, str]] = None,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Build Cytoscape.js network data structure.

    Args:
        ke_list: Set of KE IDs to include in network
        edges: DataFrame with network edges (Source_KE, Target_KE, KER_ID)
        enrichment_results: DataFrame with enrichment analysis results
        ke_title_map: Mapping of KE IDs to human-readable titles
        ke_type_map: Mapping of KE IDs to types (MIE, intermediate, AO)
        reference_sets: Optional KE_ID -> gene-set dict. When provided, KEs
            with no gene set get a `no-genes` class so the frontend can
            render them as smaller / muted nodes — making the curation gap
            visually obvious without hiding the KE from the AOP topology.
        method: 'ora' (default, Fisher's exact) or 'gsea'. Under GSEA, each
            KE node carries ``data.method='gsea'`` and ``data.nes=<clamped>``
            alongside ``data.logfc=0`` (neutral surrogate preserving the
            field for legacy .cyjs readers per D-10). Under ORA the payload
            is bytewise identical to pre-Phase-14 exports (D-10 back-compat).
        fdr_cutoff: BH-adjusted p-value below which a KE is called significant
            (issue #63). Defaults to ``Config.SIGNIFICANCE_FDR_CUTOFF`` — the
            same cutoff the enrichment table, the comparison matrix and the
            results-page FDR slider use. The raw Fisher p-value is echoed onto
            the node payload for transparency but never drives the call. The
            same cutoff decides the `depleted` class from ``FDR_depleted``
            (issue #70).
        excluded_kes: Optional KE_ID → exclusion reason map (issue #65), as
            produced by ``enrichment_service.get_ke_summary()['excluded_reasons']``.
            KEs excluded for ``'too_few_genes'`` get a `too-few-genes` class,
            KEs excluded for ``'no_mapping'`` get the `no-genes` class, and KEs
            excluded for ``'unresolved_mapping'`` get an `unresolved-mapping`
            class (issue #81), so "could not assess" never renders like
            "assessed but not enriched" and a mapped Key Event whose pathway
            could not be resolved never renders like an uncurated one.

    Returns:
        Dictionary with 'nodes' and 'edges' keys for Cytoscape.js
    """
    logger.info("Building Cytoscape network structure")

    # Build nodes
    cy_nodes = []
    for ke in ke_list:
        # Get enrichment data for this KE
        enrichment_row = enrichment_results[enrichment_results['KE'] == ke]

        if not enrichment_row.empty:
            enrichment_row = enrichment_row.iloc[0]
            odds_ratio = enrichment_row.get('odds_ratio', 0)
            # EXPO-06 / issue #50 — embed per-KE significance into the node
            # data payload so the existing client-side Download Network
            # (.cyjs) export carries it automatically. Note the DataFrame
            # column is lowercase 'p_value' (see enrichment_service.py:125).
            p_value = _to_native(enrichment_row.get('p_value'))
            fdr = _to_native(enrichment_row.get('FDR'))
            # Issue #63 — one notion of significance across the app: the
            # BH-adjusted FDR, not the raw p-value. _to_native() already
            # mapped NaN/inf to None, so a missing FDR is never significant.
            is_significant = fdr is not None and fdr < fdr_cutoff
            # Issue #70 — the depletion tail, so an under-represented KE is
            # styled as a result rather than dimmed like a null one. Absent on
            # GSEA results (NES already carries direction) and on results
            # produced before both tails were computed.
            fdr_depleted = _to_native(enrichment_row.get('FDR_depleted'))
            is_depleted = (
                not is_significant
                and fdr_depleted is not None
                and fdr_depleted < fdr_cutoff
            )
            # D-10/D-12: extract and belt-and-suspenders clamp NES under GSEA.
            # The gsea_service already clamped to ±3 but the network builder
            # must not trust upstream (two independent clamps per D-12 spec).
            nes = _to_native(enrichment_row.get('NES')) if method == 'gsea' else None
            if nes is not None:
                nes = max(-3.0, min(3.0, float(nes)))
        else:
            odds_ratio = 0
            is_significant = False
            is_depleted = False
            p_value = None
            fdr = None
            fdr_depleted = None
            nes = None

        # Get KE metadata
        label = ke_title_map.get(ke, ke)
        ke_type = ke_type_map.get(ke, "intermediate")
        has_gene_set = bool(reference_sets and reference_sets.get(ke))
        # Issue #65 — why this KE never reached a statistical test (if it didn't).
        excluded_reason = (excluded_kes or {}).get(ke)

        # Set CSS classes for styling
        classes = []
        if is_significant:
            classes.append("significant")
        if is_depleted:
            # Issue #70 — significantly under-represented. A distinct class, not
            # the significance border, so the reader can tell which way it went.
            classes.append("depleted")
        if excluded_reason == EXCLUDED_TOO_FEW_GENES:
            # Assessed-impossible, not assessed-and-null: fewer than the
            # minimum number of the KE's genes were measured (issue #65).
            classes.append("too-few-genes")
        if excluded_reason == EXCLUDED_UNRESOLVED_MAPPING:
            # Issue #81 — the KE *is* mapped; the mapping resolved to no genes.
            # Styled apart from the uncurated case because the two ask for
            # different repairs: one a Builder mapping, the other this tool's
            # reference data.
            classes.append("unresolved-mapping")
        elif excluded_reason == EXCLUDED_NO_MAPPING or (
            reference_sets is not None and not has_gene_set and not excluded_reason
        ):
            # No curated gene set at all. Note the reason, when present, is
            # authoritative: a mapped-but-unresolvable KE also has no gene set
            # here and must not be dressed up as a curation gap (issue #81).
            classes.append("no-genes")

        # D-10: build node payload with method-aware fields.
        # 'method' is always included so the frontend can branch on a single
        # attribute. Under GSEA, 'nes' is added and 'logfc' is set to 0 (a
        # neutral surrogate that renders mid-gradient grey for legacy .cyjs
        # readers that pre-date Phase 14 — preserves field shape per D-10).
        # Under ORA, 'logfc' carries the odds_ratio as before (no change).
        node_payload = {
            "id": ke,
            "label": label,
            "ke_type": ke_type,
            "has_gene_set": has_gene_set,
            "p_value": p_value,   # EXPO-06: raw Fisher p, reported not gating (#63)
            "fdr": fdr,           # EXPO-06: per-KE significance in .cyjs export
            # Issue #70: depletion tail, so a stored network (.cyjs, batch
            # report, shared link) can restyle for direction without the
            # enrichment table beside it. None on GSEA and on older results.
            "fdr_depleted": fdr_depleted,
            "method": method,     # D-10: frontend gradient selector
            # Issue #65: None when the KE was tested, otherwise the reason it
            # was not. Lets consumers of a stored network (batch report, .cyjs)
            # recover the tested/excluded accounting without a second payload.
            "excluded_reason": excluded_reason,
        }
        if method == 'gsea':
            node_payload["nes"] = nes
            node_payload["logfc"] = 0  # neutral surrogate per D-10 back-compat
        else:
            node_payload["logfc"] = odds_ratio  # ORA: odds ratio as color surrogate

        node_data = {
            "data": node_payload,
            "classes": " ".join(classes)
        }

        cy_nodes.append(node_data)
    
    # Build edges
    cy_edges = []
    for _, row in edges.iterrows():
        source_ke = row['Source_KE']
        target_ke = row['Target_KE']
        ker_id = row['KER_ID']
        
        # Only include edges where both nodes exist in our KE list
        if source_ke in ke_list and target_ke in ke_list:
            edge_data = {
                "data": {
                    "source": source_ke,
                    "target": target_ke,
                    "id": f"KER:{ker_id}"
                }
            }
            cy_edges.append(edge_data)
    
    network = {
        "nodes": cy_nodes,
        "edges": cy_edges
    }
    
    logger.info(f"Network built: {len(cy_nodes)} nodes, {len(cy_edges)} edges")
    
    return network


def ke_accounting_from_network(network_json: Any) -> Optional[Dict[str, Any]]:
    """Rebuild the tested/excluded KE accounting from a stored network (issue #65).

    Neither ConditionRecord (batch) nor SharedResult (public links) has a field
    for the enrichment summary dict, so the counts are recovered from the
    per-node ``excluded_reason`` written by ``build_cytoscape_network``. Returns
    None for networks predating that field, in which case the caller simply
    omits the accounting line.

    Args:
        network_json: Parsed Cytoscape network dict, or a JSON string.

    Returns:
        Summary dict shaped like ``enrichment_service.get_ke_summary()`` (minus
        the per-KE ``excluded_reasons`` map), or None when unavailable.
    """
    try:
        if isinstance(network_json, str):
            network_json = json.loads(network_json)
        if not network_json:
            return None
        ke_nodes = [
            n.get('data', {}) for n in network_json.get('nodes', [])
            if (n.get('data') or {}).get('ke_type') is not None
        ]
        if not ke_nodes or not any('excluded_reason' in d for d in ke_nodes):
            return None
        reasons = [d.get('excluded_reason') for d in ke_nodes]
        return {
            'total_kes': len(ke_nodes),
            'tested': sum(1 for r in reasons if not r),
            'excluded_no_mapping': sum(1 for r in reasons if r == EXCLUDED_NO_MAPPING),
            # Issue #81 — absent from networks stored before the distinction
            # existed, which is exactly right: those runs did not know it, and
            # a zero here keeps the clause out of the rendered sentence.
            'excluded_unresolved_mapping': sum(
                1 for r in reasons if r == EXCLUDED_UNRESOLVED_MAPPING
            ),
            'excluded_too_few_genes': sum(1 for r in reasons if r == EXCLUDED_TOO_FEW_GENES),
            'excluded_error': sum(1 for r in reasons if r == EXCLUDED_ERROR),
            'min_ke_genes': Config.MIN_KE_GENES,
        }
    except (json.JSONDecodeError, TypeError, AttributeError) as exc:
        logger.warning(f"ke_accounting_from_network failed: {exc}")
        return None
