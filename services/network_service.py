"""
Network visualization service for AOP pathways.
"""
import pandas as pd
import logging
from typing import Dict, Set, List, Any

logger = logging.getLogger(__name__)

def build_cytoscape_network(
    ke_list: Set[str],
    edges: pd.DataFrame, 
    enrichment_results: pd.DataFrame,
    ke_title_map: Dict[str, str],
    ke_type_map: Dict[str, str]
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Build Cytoscape.js network data structure.
    
    Args:
        ke_list: Set of KE IDs to include in network
        edges: DataFrame with network edges (Source_KE, Target_KE, KER_ID)
        enrichment_results: DataFrame with enrichment analysis results
        ke_title_map: Mapping of KE IDs to human-readable titles
        ke_type_map: Mapping of KE IDs to types (MIE, intermediate, AO)
    
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
        else:
            odds_ratio = 0
            is_significant = False
        
        # Get KE metadata
        label = ke_title_map.get(ke, ke)
        ke_type = ke_type_map.get(ke, "intermediate")
        
        # Set CSS classes for styling
        classes = []
        if is_significant:
            classes.append("significant")
        
        node_data = {
            "data": {
                "id": ke,
                "label": label,
                "logfc": odds_ratio,  # Using odds ratio as surrogate for color scaling
                "ke_type": ke_type
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