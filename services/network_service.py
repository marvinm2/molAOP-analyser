"""
Network visualization service for AOP pathways.
"""
import math
import pandas as pd
import logging
from typing import Dict, Set, List, Any

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
            is_significant = enrichment_row.get('FDR', 1.0) < 0.05
            # EXPO-06 / issue #50 — embed per-KE significance into the node
            # data payload so the existing client-side Download Network
            # (.cyjs) export carries it automatically. Note the DataFrame
            # column is lowercase 'p_value' (see enrichment_service.py:125).
            p_value = _to_native(enrichment_row.get('p_value'))
            fdr = _to_native(enrichment_row.get('FDR'))
        else:
            odds_ratio = 0
            is_significant = False
            p_value = None
            fdr = None

        # Get KE metadata
        label = ke_title_map.get(ke, ke)
        ke_type = ke_type_map.get(ke, "intermediate")
        has_gene_set = bool(reference_sets and reference_sets.get(ke))

        # Set CSS classes for styling
        classes = []
        if is_significant:
            classes.append("significant")
        if reference_sets is not None and not has_gene_set:
            classes.append("no-genes")

        node_data = {
            "data": {
                "id": ke,
                "label": label,
                "logfc": odds_ratio,  # Using odds ratio as surrogate for color scaling
                "ke_type": ke_type,
                "has_gene_set": has_gene_set,
                "p_value": p_value,   # EXPO-06: per-KE significance in .cyjs export
                "fdr": fdr,           # EXPO-06: per-KE significance in .cyjs export
            },
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

def build_cytoscape_network_legacy(
    aop_id: str, 
    ker_file: str = 'data/aop_ker_edges.csv', 
    metadata: Dict[str, str] = None
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Legacy function for building Cytoscape network (kept for compatibility).
    
    Args:
        aop_id: AOP identifier
        ker_file: Path to KER edges file
        metadata: Optional metadata mapping for node labels
    
    Returns:
        Dictionary with 'nodes' and 'edges' keys for Cytoscape.js
    """
    logger.warning("Using legacy network building function")
    
    try:
        ker_df = pd.read_csv(ker_file)
        edges = ker_df[ker_df['AOP_ID'] == aop_id]
        
        # Extract unique nodes from edges
        nodes = set(edges['Source_KE']) | set(edges['Target_KE'])
        
        # Build nodes
        cytoscape_nodes = []
        for ke in nodes:
            node_data = {
                "data": {
                    "id": ke, 
                    "label": metadata.get(ke, ke) if metadata else ke
                }
            }
            cytoscape_nodes.append(node_data)
        
        # Build edges  
        cytoscape_edges = []
        for _, row in edges.iterrows():
            edge_data = {
                "data": {
                    "source": row['Source_KE'],
                    "target": row['Target_KE'], 
                    "id": f"KER:{row['KER_ID']}"
                }
            }
            cytoscape_edges.append(edge_data)
        
        return {
            "nodes": cytoscape_nodes, 
            "edges": cytoscape_edges
        }
        
    except Exception as e:
        logger.error(f"Failed to build legacy network for AOP {aop_id}: {e}")
        raise