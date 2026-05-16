"""
Hub gene analysis service.

Provides gene-to-KE reverse mapping and hub ranking for the results page hub panel.
Consumed by the /analyze route immediately after build_ke_gene_mapping().

A hub gene is a gene symbol that appears in at least HUB_THRESHOLD distinct Key Events
within the analysed AOP. The output list is sorted descending by KE count (most
connected genes first) with gene symbol as a deterministic ascending tie-break.
"""
import logging
from typing import Dict, List, Any

logger = logging.getLogger(__name__)

# SPEC-mandated named constant (HUBG-06); value fixed at 3 — not user-adjustable.
HUB_THRESHOLD = 3


def compute_hub_genes(
    ke_gene_map: Dict[str, List[Dict[str, Any]]],
    ke_title_map: Dict[str, str],
) -> List[Dict[str, Any]]:
    """Invert ke_gene_map into a ranked hub-gene list.

    Iterates over every (KE_ID, gene_list) pair in ke_gene_map and builds a
    gene-keyed index. Each gene accumulates the count of KEs it belongs to
    together with the corresponding KE identifiers and titles. Expression
    metadata (log2fc, pvalue_adj, significant) is taken from the first
    occurrence of each gene symbol in the map, which is deterministic because
    build_ke_gene_mapping() always stores the same per-gene values regardless
    of which KE list is visited first.

    The result is sorted descending by ke_count, then ascending by gene symbol
    to guarantee a deterministic ordering when multiple genes share the same
    KE count. After sorting, ``is_hub`` is set to ``True`` for every entry
    whose ke_count meets or exceeds HUB_THRESHOLD.

    Args:
        ke_gene_map: Mapping from KE_ID to a list of gene dicts as produced by
            ``enrichment_service.build_ke_gene_mapping()``. Each gene dict must
            contain at least the ``id`` key (uppercase gene symbol). The keys
            ``log2FC``, ``significant``, and ``pvalue_adj`` are optional and
            default to ``0.0``, ``False``, and ``None`` respectively.
        ke_title_map: Mapping from KE_ID to a human-readable KE title. When a
            KE_ID is absent from this map the raw KE_ID is used as the title.

    Returns:
        A list of gene dicts, one per unique gene symbol, sorted as described
        above. Each dict contains exactly these eight keys:

        - ``gene`` (str): uppercase gene symbol.
        - ``ke_count`` (int): number of distinct KEs the gene belongs to.
        - ``ke_ids`` (list[str]): ordered list of KE identifiers.
        - ``ke_titles`` (list[str]): ordered list of KE titles (index-aligned
          with ``ke_ids``).
        - ``log2fc`` (float): log2 fold-change value from the input data.
        - ``pvalue_adj`` (float | None): adjusted p-value; may be ``None`` when
          the upload did not include an adjusted p-value column.
        - ``significant`` (bool): whether the gene passed the significance
          threshold during expression processing.
        - ``is_hub`` (bool): ``True`` when ke_count >= HUB_THRESHOLD.

        Returns ``[]`` when ke_gene_map is empty — this is not an error.
    """
    if not ke_gene_map:
        logger.debug("compute_hub_genes: empty ke_gene_map — returning empty list")
        return []

    gene_index: Dict[str, Dict] = {}

    for ke_id, genes in ke_gene_map.items():
        for g in genes:
            sym = g['id']
            if sym not in gene_index:
                gene_index[sym] = {
                    'gene': sym,
                    'ke_count': 0,
                    'ke_ids': [],
                    'ke_titles': [],
                    'log2fc': g.get('log2FC', 0.0),
                    'pvalue_adj': g.get('pvalue_adj'),         # None-safe: no coercion
                    'significant': bool(g.get('significant', False)),
                }
            gene_index[sym]['ke_count'] += 1
            gene_index[sym]['ke_ids'].append(ke_id)
            gene_index[sym]['ke_titles'].append(ke_title_map.get(ke_id, ke_id))

    hub_list = sorted(
        gene_index.values(),
        key=lambda x: (-x['ke_count'], x['gene']),  # desc ke_count, asc gene symbol
    )

    for entry in hub_list:
        entry['is_hub'] = entry['ke_count'] >= HUB_THRESHOLD

    logger.debug(
        "compute_hub_genes: %d unique genes, %d hubs (threshold=%d)",
        len(hub_list),
        sum(1 for e in hub_list if e['is_hub']),
        HUB_THRESHOLD,
    )

    return hub_list
