"""Comparison service for multi-condition batch analysis.

Provides data layer utilities for building comparison matrices from batch
enrichment results, suitable for heatmap, table, and network visualisations.
"""

import json
import logging
import numpy as np
import pandas as pd
from typing import Any

logger = logging.getLogger(__name__)

# VHP4Safety palette — fixed order by upload position.
# Sky Blue (#93D5F6) excluded (too similar to Light Blue #009FE3).
CONDITION_PALETTE = [
    '#E6007E',  # 1: Primary Magenta
    '#307BBF',  # 2: Primary Blue
    '#EB5B25',  # 3: Orange
    '#45A6B2',  # 4: Teal
    '#64358C',  # 5: Violet
    '#005A6C',  # 6: Dark Teal
    '#9A1C57',  # 7: Deep Magenta (fallback)
    '#29235C',  # 8: Primary Dark (fallback)
]


def build_comparison_matrix(conditions: list) -> dict[str, Any]:
    """Build a comparison matrix from a list of ConditionRecord objects.

    Pivots enrichment results for all conditions into a KE × condition matrix
    of FDR values and -log10(FDR) significance scores.  KE rows are sorted by
    mean significance (most significant first).  Columns are kept in upload-
    position order (not alphabetical) by explicit reindexing after pivot.

    Args:
        conditions: List of ConditionRecord ORM objects ordered by position.
                    Each must expose ``condition_label`` and ``enrichment_json``
                    (a JSON-serialised list of enrichment result dicts).

    Returns:
        A dict with keys:
            ke_labels        — list of KE ID strings in row order
            ke_titles        — list of human-readable KE titles
            condition_labels — list of condition label strings in upload order
            fdr_matrix       — 2-D list (rows=KEs, cols=conditions) of raw FDR
                               floats, None where data is absent
            neg_log10_matrix — 2-D list of -log10(FDR) floats, None where FDR
                               is absent or > 0.05 (non-significant)
            condition_colors — list of hex colour strings from CONDITION_PALETTE

        Returns an empty dict ``{}`` if no enrichment rows are found.
    """
    rows = []
    ke_title_map: dict[str, str] = {}

    for cond in conditions:
        if not cond.enrichment_json:
            continue
        try:
            enrichment_list = json.loads(cond.enrichment_json)
        except (json.JSONDecodeError, TypeError) as exc:
            logger.warning(
                'build_comparison_matrix: failed to parse enrichment_json '
                'for condition %s: %s', cond.condition_label, exc
            )
            continue

        for entry in enrichment_list:
            ke = entry.get('KE')
            fdr = entry.get('FDR')
            title = entry.get('Title', '')
            if ke is None or fdr is None:
                continue
            rows.append({
                'condition': cond.condition_label,
                'KE': ke,
                'Title': title,
                'FDR': float(fdr),
            })
            # Keep first-seen title per KE
            if ke not in ke_title_map:
                ke_title_map[ke] = title

    if not rows:
        logger.info('build_comparison_matrix: no enrichment rows found — returning empty dict')
        return {}

    df = pd.DataFrame(rows, columns=['condition', 'KE', 'Title', 'FDR'])

    # Pivot to KE × condition, keeping the first FDR value where there are
    # duplicate (condition, KE) pairs (should not happen in practice).
    pivot = df.pivot_table(index='KE', columns='condition', values='FDR', aggfunc='first')

    # Reindex columns to upload-position order to prevent pandas alphabetical sort.
    condition_labels = [c.condition_label for c in conditions]
    pivot = pivot.reindex(columns=condition_labels)

    # Compute -log10(FDR) matrix: None for absent or non-significant cells.
    def _neg_log10(fdr_val: Any) -> float | None:
        """Return -log10(fdr_val) or None for missing/non-significant values."""
        if pd.isna(fdr_val) or fdr_val > 0.05:
            return None
        return float(-np.log10(max(float(fdr_val), 1e-300)))

    neg_log10_pivot = pivot.map(_neg_log10)

    # Sort KE rows by mean -log10(FDR) significance descending (most significant first).
    # Fill NaN with 0 for sorting purposes only.
    mean_sig = neg_log10_pivot.fillna(0).mean(axis=1)
    sort_order = mean_sig.sort_values(ascending=False).index
    pivot = pivot.loc[sort_order]
    neg_log10_pivot = neg_log10_pivot.loc[sort_order]

    ke_labels = list(pivot.index)
    ke_titles = [ke_title_map.get(ke, ke) for ke in ke_labels]

    # Convert to nested Python lists, replacing NaN with None for JSON safety.
    def _to_list(frame: pd.DataFrame) -> list[list]:
        result = []
        for _, row in frame.iterrows():
            result.append([
                None if pd.isna(v) else v
                for v in row
            ])
        return result

    fdr_matrix = _to_list(pivot)
    neg_log10_matrix = _to_list(neg_log10_pivot)

    # Assign colours from palette in upload-position order.
    condition_colors = CONDITION_PALETTE[:len(conditions)]

    return {
        'ke_labels': ke_labels,
        'ke_titles': ke_titles,
        'condition_labels': condition_labels,
        'fdr_matrix': fdr_matrix,
        'neg_log10_matrix': neg_log10_matrix,
        'condition_colors': condition_colors,
    }
