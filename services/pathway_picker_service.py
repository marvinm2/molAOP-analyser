"""
WP picker data assembly service.

Provides the `build_wp_picker_data()` plain function that builds the WP picker
list for the results page Pathway view section.  It is called by the /analyze
route after `run_enrichment()` and `compute_hub_genes()` have completed, so
the enrichment FDR ranking is already available.

The function is intentionally free of any gene or expression data (SPEC R4 /
D-02 — no user log2FC / expression values belong in the picker entries).
"""
import logging
import re
from typing import Any, Dict, List, Optional, Set

import pandas as pd

logger = logging.getLogger(__name__)

# WikiPathways IDs are always "WP" followed by one or more digits. Builder API
# records are third-party data the analyser does not control, so the WP ID is
# validated server-side before it reaches the picker markup / iframe URL.
_WPID_RE = re.compile(r"^WP\d+$")


def build_wp_picker_data(
    ke_list: Set[str],
    ke_wp_records: List[Dict[str, Any]],
    enrichment_results: pd.DataFrame,
    ke_title_map: Dict[str, str],
    wp_title_map: Dict[str, str],
) -> List[Dict[str, Any]]:
    """Build the ordered WP picker list for the Pathway view section.

    Iterates over all KE-WP mapping records, keeps only those whose KE belongs
    to the selected AOP (``ke_list``), resolves titles and confidence level, and
    ranks entries by FDR-based enrichment rank so the most-enriched KE pathways
    appear first in the picker dropdown.

    KE ID normalisation (Pitfall 5): Builder API records carry ``ke_id`` in the
    space-separated format ``"KE 386"`` while the enrichment DataFrame and
    ``ke_list`` use the colon format ``"KE:386"``.  The function normalises by
    calling ``.replace(" ", ":")`` on every incoming KE ID before lookups.

    For CSV fallback records (only ``KE_ID`` / ``WP_ID`` columns present):
    - ``wp_title`` falls back to the WP ID itself (``wp_title_map.get(wp_id, wp_id)``).
    - ``confidence_level`` is ``""`` (empty string) when the source record has no
      ``confidence_level`` key, signalling the template to omit the badge entirely.

    Args:
        ke_list: Set of KE IDs (colon format, e.g. ``{"KE:386", "KE:123"}``) that
            belong to the selected AOP.  Records for other KEs are excluded.
        ke_wp_records: Raw KE-WP mapping records from either the Builder API
            (``fetch_ke_wp_records``) or the CSV fallback
            (``load_ke_wp_records_csv``).  API records carry ``ke_id``,
            ``pathway_id``, ``pathway_title``, and ``confidence_level``.  CSV
            records carry only ``KE_ID`` and ``WP_ID``.
        enrichment_results: DataFrame produced by ``run_enrichment()``, sorted
            ascending by FDR.  Must contain a ``KE`` column (colon-format IDs).
            Used to derive the 0-based enrichment rank for each KE.
        ke_title_map: Mapping from KE ID (colon format) to human-readable KE
            title.  When a KE is absent the raw ID is used as the title.
        wp_title_map: Mapping from WP ID (e.g. ``"WP368"``) to pathway title.
            When absent (CSV fallback) the WP ID itself is used as the title.

    Returns:
        A list of dicts, one per KE-WP record whose KE is in ``ke_list``,
        sorted ascending by ``enrichment_rank``.  Each dict contains exactly
        these six keys:

        - ``ke_id`` (str): KE ID in colon format (``"KE:386"``).
        - ``ke_title`` (str): Human-readable KE title from ``ke_title_map``.
        - ``wp_id`` (str): WikiPathways ID (``"WP368"``).
        - ``wp_title`` (str): Pathway title from ``wp_title_map``, or the WP ID.
        - ``confidence_level`` (str): Lowercased confidence level
          (``"high"``/``"medium"``/``"low"``), or ``""`` when absent.
        - ``enrichment_rank`` (int): 0-based index of the KE in
          ``enrichment_results``; ``len(enrichment_results)`` for KEs absent from
          enrichment (treated as least-enriched, placed last).

        Returns ``[]`` (an empty list, not ``None``) when no record matches
        ``ke_list`` — this is not an error, and lets
        ``{% if wp_picker_data %}`` be falsy in the template.
    """
    if not ke_wp_records:
        logger.debug("build_wp_picker_data: empty ke_wp_records — returning []")
        return []

    # Build enrichment rank lookup keyed on normalised KE id (colon format).
    # KEs absent from enrichment_results get a rank equal to len(enrichment_results)
    # so they sort last (least enriched).
    inf_rank = len(enrichment_results)
    enrichment_rank_map: Dict[str, int] = {
        row["KE"]: idx
        for idx, row in enumerate(enrichment_results.to_dict("records"))
    }

    entries: List[Dict[str, Any]] = []

    for record in ke_wp_records:
        # Resolve KE ID — API uses space format ("KE 386"), CSV uses colon ("KE:386")
        raw_ke_id = record.get("KE_ID") or record.get("ke_id", "")
        ke_id = raw_ke_id.replace(" ", ":")

        # Skip records whose KE is not in the selected AOP
        if ke_id not in ke_list:
            continue

        # Resolve WP ID — API uses "pathway_id", CSV uses "WP_ID"
        wp_id = record.get("WP_ID") or record.get("pathway_id", "")

        # Validate WP ID server-side (defence in depth — the client-side regex
        # in setPathway() only protects the iframe). A malformed pathway_id is
        # skipped rather than producing a silently-broken picker option.
        if not _WPID_RE.match(wp_id):
            logger.warning(
                "build_wp_picker_data: skipping malformed wp_id %r", wp_id
            )
            continue

        entry = {
            "ke_id": ke_id,
            "ke_title": ke_title_map.get(ke_id, ke_id),
            "wp_id": wp_id,
            "wp_title": wp_title_map.get(wp_id, wp_id),
            "confidence_level": (record.get("confidence_level") or "").lower(),
            "enrichment_rank": enrichment_rank_map.get(ke_id, inf_rank),
        }
        entries.append(entry)

    # Sort ascending by enrichment rank (most-enriched KE first)
    sorted_entries = sorted(entries, key=lambda x: x["enrichment_rank"])

    logger.debug(
        "build_wp_picker_data: %d entries built from %d records (%d KEs in AOP)",
        len(sorted_entries),
        len(ke_wp_records),
        len(ke_list),
    )

    return sorted_entries
