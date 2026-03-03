"""
AOP Discovery Service.

Provides a three-tier fallback mechanism for fetching the list of all AOPs
with their KE counts and KE-gene mapping coverage:

    1. Disk cache (1-week TTL, shared across Gunicorn workers)
    2. Live fetch: SPARQL (all AOPs + KE membership) + Builder API (mapped KE IDs)
    3. Hardcoded fallback from Config.CASE_STUDY_AOPS

Each AOP entry has the shape::

    {
        "aop_id":         "AOP:619",
        "title":          "Some title",
        "ke_count":       12,
        "mapped_ke_count": 4,
    }

Mapped AOPs (mapped_ke_count > 0) are sorted first by mapped_ke_count descending;
unmapped AOPs follow sorted by aop_id.
"""

import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Set

from SPARQLWrapper import JSON, SPARQLWrapper

logger = logging.getLogger(__name__)

# Disk-cache key for the full AOP list.
AOP_LIST_CACHE_KEY = "aop_discovery_list_v1"

# 1-week TTL in seconds
AOP_LIST_CACHE_TTL = 7 * 24 * 3600


# ---------------------------------------------------------------------------
# URI helpers
# ---------------------------------------------------------------------------

def _uri_to_aop_id(uri: str) -> str:
    """Convert an identifiers.org AOP URI to a compact AOP ID string.

    Example::

        >>> _uri_to_aop_id("https://identifiers.org/aop/619")
        'AOP:619'

    Args:
        uri: Full identifiers.org AOP URI.

    Returns:
        Compact AOP ID like ``"AOP:619"``.

    Raises:
        ValueError: If the numeric AOP ID cannot be extracted from the URI.
    """
    match = re.search(r'/aop/(\d+)$', uri)
    if match:
        return f"AOP:{match.group(1)}"
    raise ValueError(f"Cannot parse AOP ID from URI: {uri}")


# ---------------------------------------------------------------------------
# SPARQL fetchers
# ---------------------------------------------------------------------------

def fetch_all_aops_from_sparql(
    endpoint: str = None,
    timeout: int = None,
) -> List[Dict]:
    """Fetch all AOPs with titles and KE counts from AOP-Wiki SPARQL.

    Runs a single aggregating SPARQL query that groups by AOP and counts
    distinct Key Events.  Returns only AOPs that have at least one KE.

    Args:
        endpoint: SPARQL endpoint URL.  Defaults to ``Config.SPARQL_ENDPOINT``.
        timeout: Query timeout in seconds.  Defaults to ``Config.SPARQL_TIMEOUT``.

    Returns:
        List of dicts with keys ``aop_id``, ``title``, ``ke_count``,
        ``mapped_ke_count`` (always 0 at this stage).

    Raises:
        Exception: On SPARQL network/timeout errors.
    """
    from config import Config  # local import to avoid circular deps

    endpoint = endpoint or Config.SPARQL_ENDPOINT
    timeout = timeout or Config.SPARQL_TIMEOUT

    query = """
PREFIX aopo: <http://aopkb.org/aop_ontology#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
SELECT ?aop ?title (COUNT(DISTINCT ?ke) AS ?ke_count)
WHERE {
  ?aop a aopo:AdverseOutcomePathway ;
       dc:title ?title ;
       aopo:has_key_event ?ke .
}
GROUP BY ?aop ?title
"""

    sparql = SPARQLWrapper(endpoint)
    sparql.setTimeout(timeout)
    sparql.setReturnFormat(JSON)
    sparql.setQuery(query)

    results = sparql.query().convert()
    bindings = results.get("results", {}).get("bindings", [])

    aops = []
    for row in bindings:
        try:
            aop_id = _uri_to_aop_id(row["aop"]["value"])
            title = row.get("title", {}).get("value", aop_id)
            ke_count = int(row.get("ke_count", {}).get("value", 0))
            aops.append({
                "aop_id": aop_id,
                "title": title,
                "ke_count": ke_count,
                "mapped_ke_count": 0,
            })
        except (ValueError, KeyError) as exc:
            logger.debug("Skipping malformed SPARQL row: %s — %s", row, exc)

    logger.info("SPARQL: fetched %d AOPs from AOP-Wiki", len(aops))
    return aops


def _build_ke_to_aop_map(
    endpoint: str = None,
    timeout: int = None,
) -> Dict[str, List[str]]:
    """Fetch all KE->AOP memberships from AOP-Wiki SPARQL in a single query.

    Args:
        endpoint: SPARQL endpoint URL.  Defaults to ``Config.SPARQL_ENDPOINT``.
        timeout: Query timeout in seconds.  Defaults to ``Config.SPARQL_TIMEOUT``.

    Returns:
        Dict mapping KE ID (e.g. ``"KE:709"``) to a list of AOP IDs
        (e.g. ``["AOP:38", "AOP:105"]``).

    Raises:
        Exception: On SPARQL network/timeout errors.
    """
    from config import Config  # local import to avoid circular deps

    endpoint = endpoint or Config.SPARQL_ENDPOINT
    timeout = timeout or Config.SPARQL_TIMEOUT

    query = """
PREFIX aopo: <http://aopkb.org/aop_ontology#>
SELECT DISTINCT ?ke ?aop
WHERE {
  ?aop a aopo:AdverseOutcomePathway ;
       aopo:has_key_event ?ke .
}
"""

    sparql = SPARQLWrapper(endpoint)
    sparql.setTimeout(timeout)
    sparql.setReturnFormat(JSON)
    sparql.setQuery(query)

    results = sparql.query().convert()
    bindings = results.get("results", {}).get("bindings", [])

    ke_to_aop: Dict[str, List[str]] = {}
    for row in bindings:
        try:
            # Parse KE URI: https://identifiers.org/aop.events/709 -> KE:709
            ke_uri = row["ke"]["value"]
            ke_match = re.search(r'/aop\.events/(\d+)$', ke_uri)
            if not ke_match:
                continue
            ke_id = f"KE:{ke_match.group(1)}"

            aop_id = _uri_to_aop_id(row["aop"]["value"])
            ke_to_aop.setdefault(ke_id, []).append(aop_id)
        except (ValueError, KeyError) as exc:
            logger.debug("Skipping malformed KE->AOP row: %s — %s", row, exc)

    logger.info(
        "SPARQL: built KE->AOP map with %d KEs covering %d unique AOPs",
        len(ke_to_aop),
        len({a for aops in ke_to_aop.values() for a in aops}),
    )
    return ke_to_aop


# ---------------------------------------------------------------------------
# Builder API fetcher
# ---------------------------------------------------------------------------

def fetch_mapped_ke_ids_from_builder(config) -> Set[str]:
    """Fetch all KE IDs that have Builder API gene-set mappings.

    Uses the existing ``api_service._make_api_session`` and
    ``api_service.fetch_all_ke_wp_mappings`` to retrieve all KE-WP
    mapping records, then extracts and normalises the KE IDs.

    If ``config.BUILDER_API_URL`` is empty, returns an empty set
    (graceful degradation — all AOPs will appear unmapped).

    Args:
        config: Application config object (must have ``BUILDER_API_URL``
                and ``BUILDER_API_TIMEOUT`` attributes).

    Returns:
        Set of normalised KE IDs like ``{'KE:55', 'KE:709'}``.
    """
    if not config.BUILDER_API_URL:
        logger.info(
            "BUILDER_API_URL is empty; skipping mapped KE fetch — "
            "all AOPs will show mapped_ke_count=0"
        )
        return set()

    try:
        from services.api_service import _make_api_session, fetch_all_ke_wp_mappings

        session = _make_api_session()
        records = fetch_all_ke_wp_mappings(
            session=session,
            base_url=config.BUILDER_API_URL,
            timeout=config.BUILDER_API_TIMEOUT,
        )

        # Normalise "KE 55" -> "KE:55"
        mapped_ids: Set[str] = set()
        for record in records:
            raw_ke_id = record.get("ke_id", "")
            if raw_ke_id:
                mapped_ids.add(raw_ke_id.replace(" ", ":"))

        logger.info(
            "Builder API: found %d mapped KE IDs across %d records",
            len(mapped_ids),
            len(records),
        )
        return mapped_ids

    except Exception as exc:
        logger.warning("Builder API fetch for mapped KEs failed: %s", exc)
        return set()


# ---------------------------------------------------------------------------
# AOP list builder
# ---------------------------------------------------------------------------

def build_aop_list(config) -> List[Dict]:
    """Build the complete AOP list with mapping coverage annotations.

    Fetches AOPs and KE membership from SPARQL and mapped KE IDs from
    the Builder API in parallel, then cross-references to compute
    ``mapped_ke_count`` for each AOP.

    Sorting:
        - Mapped AOPs (``mapped_ke_count > 0``) first, sorted by
          ``mapped_ke_count`` descending.
        - Unmapped AOPs second, sorted by numeric AOP ID ascending.

    Args:
        config: Application config object.

    Returns:
        List of AOP dicts with keys ``aop_id``, ``title``, ``ke_count``,
        ``mapped_ke_count``.

    Raises:
        Exception: If the SPARQL fetch fails (Builder API failure is
                   handled gracefully; SPARQL failure is propagated).
    """
    # Fetch SPARQL data and Builder API data in parallel
    aops: List[Dict] = []
    ke_to_aop: Dict[str, List[str]] = {}
    mapped_ke_ids: Set[str] = set()

    with ThreadPoolExecutor(max_workers=3) as executor:
        future_aops = executor.submit(fetch_all_aops_from_sparql)
        future_ke_map = executor.submit(_build_ke_to_aop_map)
        future_builder = executor.submit(fetch_mapped_ke_ids_from_builder, config)

        for future in as_completed([future_aops, future_ke_map, future_builder]):
            if future is future_aops:
                aops = future.result()  # may raise — propagated intentionally
            elif future is future_ke_map:
                try:
                    ke_to_aop = future.result()
                except Exception as exc:
                    logger.warning("KE->AOP map fetch failed: %s", exc)
            else:
                try:
                    mapped_ke_ids = future.result()
                except Exception as exc:
                    logger.warning("Builder API fetch failed: %s", exc)

    # Build AOP -> mapped KE count index from the KE->AOP membership map
    aop_to_mapped_count: Dict[str, int] = {}
    for ke_id, aop_ids in ke_to_aop.items():
        if ke_id in mapped_ke_ids:
            for aop_id in aop_ids:
                aop_to_mapped_count[aop_id] = aop_to_mapped_count.get(aop_id, 0) + 1

    # Annotate each AOP with its mapped KE count
    for aop in aops:
        aop["mapped_ke_count"] = aop_to_mapped_count.get(aop["aop_id"], 0)

    # Merge in configured AOPs from Config that are absent from SPARQL results
    from config import Config  # local import to avoid circular deps
    existing_ids = {a["aop_id"] for a in aops}
    for entry in Config.CASE_STUDY_AOPS.values():
        if not entry.get("enabled", False) or not entry.get("id"):
            continue
        aop_id = entry["id"]
        if aop_id not in existing_ids:
            aops.append({
                "aop_id": aop_id,
                "title": entry.get("label", aop_id),
                "ke_count": 0,
                "mapped_ke_count": 1,  # treat configured AOPs as mapped
            })
            logger.info("Merged configured AOP %s into discovery list", aop_id)

    # Sort: mapped first (by mapped_ke_count desc), then unmapped (by numeric AOP id)
    def _sort_key(aop: Dict):
        is_mapped = aop["mapped_ke_count"] > 0
        # Negative mapped count so higher counts sort first
        mapped_order = -aop["mapped_ke_count"]
        # Numeric AOP ID for stable secondary sort
        try:
            numeric_id = int(aop["aop_id"].replace("AOP:", ""))
        except ValueError:
            numeric_id = 999999
        return (not is_mapped, mapped_order, numeric_id)

    aops.sort(key=_sort_key)
    logger.info(
        "build_aop_list: %d total AOPs, %d with mapped KEs",
        len(aops),
        sum(1 for a in aops if a["mapped_ke_count"] > 0),
    )
    return aops


# ---------------------------------------------------------------------------
# Hardcoded fallback
# ---------------------------------------------------------------------------

def _hardcoded_aop_list() -> List[Dict]:
    """Return a minimal AOP list extracted from Config.CASE_STUDY_AOPS.

    Used as the last-resort fallback when both SPARQL and the cache are
    unavailable.  Only entries with ``enabled=True`` and an ``id`` field
    are included.  ``mapped_ke_count`` is set to 1 so these AOPs appear
    at the top of the typeahead (they are known-good configured entries).

    Returns:
        List of AOP dicts compatible with the standard shape.
    """
    from config import Config  # local import to avoid circular deps at module level

    result = []
    for entry in Config.CASE_STUDY_AOPS.values():
        if not entry.get("enabled", False):
            continue
        aop_id = entry.get("id")
        if not aop_id:
            continue
        result.append({
            "aop_id": aop_id,
            "title": entry.get("label", aop_id),
            "ke_count": 0,       # unknown without SPARQL
            "mapped_ke_count": 1,  # treat configured AOPs as having some mapping
        })

    logger.info("Hardcoded fallback: returning %d AOPs", len(result))
    return result


# ---------------------------------------------------------------------------
# Three-tier entry point
# ---------------------------------------------------------------------------

def get_aop_list(config) -> List[Dict]:
    """Return the AOP list using three-tier fallback.

    Tier 1 — Disk cache:
        If a valid cached list exists (TTL 1 week), return it immediately.

    Tier 2 — Live fetch:
        Call :func:`build_aop_list` to fetch from SPARQL + Builder API,
        cache the result, and return it.

    Tier 3 — Hardcoded config:
        If both cache and live fetch fail, return the list extracted from
        ``Config.CASE_STUDY_AOPS`` so the app always has something to show.

    Args:
        config: Application config object.

    Returns:
        Non-empty list of AOP dicts.
    """
    from cache_manager import get_reference_cache  # local import

    reference_cache = get_reference_cache(config.CACHE_DIR)

    # --- Tier 1: disk cache ---
    cached = reference_cache.get(AOP_LIST_CACHE_KEY)
    if cached is not None:
        logger.info(
            "AOP list loaded from disk cache (%d entries)", len(cached)
        )
        return cached

    # --- Tier 2: live fetch ---
    try:
        aops = build_aop_list(config)
        if aops:
            reference_cache.set(AOP_LIST_CACHE_KEY, aops, expire=AOP_LIST_CACHE_TTL)
            logger.info(
                "AOP list built and cached: %d entries (TTL %ds)",
                len(aops),
                AOP_LIST_CACHE_TTL,
            )
            return aops
    except Exception as exc:
        logger.error(
            "build_aop_list failed; falling back to hardcoded list: %s", exc
        )

    # --- Tier 3: hardcoded fallback ---
    fallback = _hardcoded_aop_list()
    logger.warning(
        "Using hardcoded AOP fallback (%d entries)", len(fallback)
    )
    return fallback
