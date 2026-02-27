"""
SPARQL client for fetching AOP structure data from AOP-Wiki.

Queries the AOP-Wiki RDF endpoint for Key Events, Key Event Relationships,
and KE metadata (titles, types) for given AOP IDs.
"""
import logging
import re
from typing import Dict, Set, Tuple, List

import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON

from cache_manager import cache, CacheManager
from config import Config
from exceptions import AOPDataError

logger = logging.getLogger(__name__)


def _uri_to_ke_id(uri: str) -> str:
    """Convert AOP-Wiki URI to KE ID string.

    Example:
        https://identifiers.org/aop.events/709 -> KE:709
    """
    match = re.search(r'/aop\.events/(\d+)$', uri)
    if match:
        return f"KE:{match.group(1)}"
    raise ValueError(f"Cannot parse KE ID from URI: {uri}")


def _aop_id_to_uri(aop_id: str) -> str:
    """Convert AOP ID string to identifiers.org URI.

    Example:
        AOP:105 -> <https://identifiers.org/aop/105>
    """
    num = aop_id.replace("AOP:", "")
    return f"<https://identifiers.org/aop/{num}>"


def fetch_aop_ke_data(
    aop_ids: List[str],
    endpoint: str = None,
    timeout: int = None,
) -> Tuple[Set[str], pd.DataFrame, Dict[str, str], Dict[str, str]]:
    """Fetch KE data for one or more AOPs from the AOP-Wiki SPARQL endpoint.

    Args:
        aop_ids: List of AOP identifiers, e.g. ["AOP:105", "AOP:138"]
        endpoint: SPARQL endpoint URL (defaults to Config.SPARQL_ENDPOINT)
        timeout: Query timeout in seconds (defaults to Config.SPARQL_TIMEOUT)

    Returns:
        Tuple of (ke_set, edges_df, ke_type_map, ke_title_map)
        - ke_set: set of KE ID strings like {"KE:709", "KE:710"}
        - edges_df: DataFrame with columns [AOP_ID, Source_KE, Target_KE, KER_ID]
        - ke_type_map: dict mapping KE ID -> type (MIE/AO/intermediate)
        - ke_title_map: dict mapping KE ID -> title string

    Raises:
        AOPDataError: On network failure, timeout, or empty results
    """
    endpoint = endpoint or Config.SPARQL_ENDPOINT
    timeout = timeout or Config.SPARQL_TIMEOUT

    values_clause = " ".join(_aop_id_to_uri(aid) for aid in aop_ids)

    # Query 1: KEs with titles and types
    ke_query = f"""
    PREFIX aopo: <http://aopkb.org/aop_ontology#>
    PREFIX dc: <http://purl.org/dc/elements/1.1/>
    SELECT DISTINCT ?ke ?title
      (IF(BOUND(?is_mie), 'MIE', IF(BOUND(?is_ao), 'AO', 'intermediate')) AS ?type)
    WHERE {{
      VALUES ?aop {{ {values_clause} }}
      ?aop aopo:has_key_event ?ke .
      OPTIONAL {{ ?ke dc:title ?title }}
      OPTIONAL {{ ?aop aopo:has_molecular_initiating_event ?ke . BIND(true AS ?is_mie) }}
      OPTIONAL {{ ?aop aopo:has_adverse_outcome ?ke . BIND(true AS ?is_ao) }}
    }}
    """

    # Query 2: KER edges
    ker_query = f"""
    PREFIX aopo: <http://aopkb.org/aop_ontology#>
    SELECT DISTINCT ?upstream ?downstream WHERE {{
      VALUES ?aop {{ {values_clause} }}
      ?aop aopo:has_key_event_relationship ?ker .
      ?ker aopo:has_upstream_key_event ?upstream ;
           aopo:has_downstream_key_event ?downstream .
    }}
    """

    sparql = SPARQLWrapper(endpoint)
    sparql.setTimeout(timeout)
    sparql.setReturnFormat(JSON)

    # --- Run KE query ---
    try:
        sparql.setQuery(ke_query)
        ke_results = sparql.query().convert()
    except Exception as e:
        logger.error(f"SPARQL KE query failed for {aop_ids}: {e}")
        raise AOPDataError(
            f"Failed to fetch KE data from AOP-Wiki: {e}",
            aop_id=str(aop_ids),
            missing_data="ke_data",
        )

    bindings = ke_results.get("results", {}).get("bindings", [])
    if not bindings:
        raise AOPDataError(
            f"No Key Events found for {aop_ids} in AOP-Wiki",
            aop_id=str(aop_ids),
            missing_data="ke_list",
        )

    ke_set: Set[str] = set()
    ke_type_map: Dict[str, str] = {}
    ke_title_map: Dict[str, str] = {}

    # Priority for type conflicts: MIE > AO > intermediate
    type_priority = {"MIE": 3, "AO": 2, "intermediate": 1}

    for row in bindings:
        ke_uri = row["ke"]["value"]
        ke_id = _uri_to_ke_id(ke_uri)
        ke_set.add(ke_id)

        title = row.get("title", {}).get("value", ke_id)
        ke_title_map[ke_id] = title

        ke_type = row.get("type", {}).get("value", "intermediate")
        existing_priority = type_priority.get(ke_type_map.get(ke_id, ""), 0)
        new_priority = type_priority.get(ke_type, 1)
        if new_priority > existing_priority:
            ke_type_map[ke_id] = ke_type

    logger.info(f"SPARQL: fetched {len(ke_set)} KEs for {aop_ids}")

    # --- Run KER query ---
    try:
        sparql.setQuery(ker_query)
        ker_results = sparql.query().convert()
    except Exception as e:
        logger.error(f"SPARQL KER query failed for {aop_ids}: {e}")
        raise AOPDataError(
            f"Failed to fetch KER data from AOP-Wiki: {e}",
            aop_id=str(aop_ids),
            missing_data="ker_data",
        )

    edge_rows = []
    for row in ker_results.get("results", {}).get("bindings", []):
        upstream = _uri_to_ke_id(row["upstream"]["value"])
        downstream = _uri_to_ke_id(row["downstream"]["value"])
        edge_rows.append({
            "AOP_ID": ",".join(aop_ids),
            "upstream": upstream,
            "downstream": downstream,
        })

    edges_df = pd.DataFrame(edge_rows, columns=["AOP_ID", "upstream", "downstream"])
    # Rename to match network_service.build_cytoscape_network() expectations
    edges_df = edges_df.rename(columns={"upstream": "Source_KE", "downstream": "Target_KE"})
    # Add synthetic KER_ID column for network_service compatibility
    edges_df["KER_ID"] = range(len(edges_df))
    logger.info(f"SPARQL: fetched {len(edges_df)} edges for {aop_ids}")

    return ke_set, edges_df, ke_type_map, ke_title_map


def fetch_aop_ke_data_cached(
    aop_ids: List[str],
    endpoint: str = None,
    timeout: int = None,
) -> Tuple[Set[str], pd.DataFrame, Dict[str, str], Dict[str, str]]:
    """Cached wrapper around fetch_aop_ke_data.

    Uses the in-memory CacheManager with 1-hour TTL.
    Cache key is based on sorted AOP IDs.

    Args:
        aop_ids: List of AOP identifiers
        endpoint: SPARQL endpoint URL
        timeout: Query timeout in seconds

    Returns:
        Same as fetch_aop_ke_data
    """
    sorted_ids = sorted(aop_ids)
    cache_key = f"sparql_aop_data_{'_'.join(sorted_ids)}"

    return cache.get_or_set(
        cache_key,
        lambda: fetch_aop_ke_data(sorted_ids, endpoint, timeout),
        ttl=Config.CACHE_TTL,
    )
