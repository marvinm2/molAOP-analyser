"""Comparison service for multi-condition batch analysis.

Provides data layer utilities for building comparison matrices from batch
enrichment results, suitable for heatmap, table, and network visualisations.
"""

import json
import logging
import numpy as np
import pandas as pd
from typing import Any

from config import Config

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


def build_comparison_matrix(conditions: list, method: str = 'ora') -> dict[str, Any]:
    """Build a comparison matrix from a list of ConditionRecord objects.

    Pivots enrichment results for all conditions into a KE × condition matrix
    of FDR values and -log10(FDR) significance scores.  KE rows are sorted by
    mean significance (most significant first).  Columns are kept in upload-
    position order (not alphabetical) by explicit reindexing after pivot.

    For a GSEA batch (``method='gsea'``, issue #76) an additional ``nes_matrix``
    is produced from each condition's normalised enrichment score. NES is the
    quantity worth comparing across conditions — it is signed, so a coordinated
    upward shift in one condition and a downward shift in another are
    distinguishable, which a p-value-derived matrix cannot show. Unlike
    ``neg_log10_matrix`` it is *not* blanked at the significance cutoff: a
    consistent sub-threshold NES across a dose series is itself the signal, and
    this matrix is the record the table and the exports read.

    Issue #107 qualifies that for *colour* specifically. The sign of a
    sub-threshold NES is not a stable quantity — the same conditions run
    against a harmonised and a per-file background returned opposite signs for
    four of six KEs in the weakest condition, all with FDR >= 0.48 — so both
    heatmap renderers mask this matrix at the cutoff before colouring it. The
    masking lives with the rendering, not here: what is unsuitable for a
    diverging colour scale is still a legitimate number to read and to export.

    Args:
        conditions: List of ConditionRecord ORM objects ordered by position.
                    Each must expose ``condition_label`` and ``enrichment_json``
                    (a JSON-serialised list of enrichment result dicts).
        method: Enrichment method the batch was run with — ``'ora'`` (default,
                the value a pre-#76 batch reads back as) or ``'gsea'``.

    Returns:
        A dict with keys:
            ke_labels        — list of KE ID strings in row order
            ke_titles        — list of human-readable KE titles
            condition_labels — list of condition label strings in upload order
            method           — the method string this matrix was built for
            fdr_matrix       — 2-D list (rows=KEs, cols=conditions) of raw FDR
                               floats, None where data is absent
            neg_log10_matrix — 2-D list of -log10(FDR) floats, None where FDR
                               is absent or above Config.SIGNIFICANCE_FDR_CUTOFF
                               (non-significant)
            nes_matrix       — 2-D list of NES floats (GSEA batches only; all
                               None for ORA batches, which carry no NES)
            nes_status_matrix— 2-D list of ``nes_status`` strings (#117), None
                               where the condition carries no status. A None
                               NES with status
                               ``beyond_permutation_resolution`` is a Key Event
                               that beat every permutation, not an untested one
            es_matrix        — 2-D list of raw enrichment scores, the only
                               surviving direction for those cells
            condition_colors — list of hex colour strings from CONDITION_PALETTE

        Returns an empty dict ``{}`` if no enrichment rows are found.
    """
    is_gsea = method == 'gsea'
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
            nes = entry.get('NES')
            es = entry.get('ES')
            rows.append({
                'condition': cond.condition_label,
                'KE': ke,
                'Title': title,
                'FDR': float(fdr),
                'NES': float(nes) if nes is not None else None,
                # Issue #117 — a missing NES is not one thing. It is absent on
                # an untested Key Event and it is absent on a Key Event that
                # beat every permutation, and those are opposite results. The
                # status says which; the ES is the only surviving statement of
                # direction in the second case, and the cross-condition view is
                # the one place a reader looks for direction.
                'nes_status': entry.get('nes_status'),
                'ES': float(es) if es is not None else None,
            })
            # Keep first-seen title per KE
            if ke not in ke_title_map:
                ke_title_map[ke] = title

    if not rows:
        logger.info('build_comparison_matrix: no enrichment rows found — returning empty dict')
        return {}

    df = pd.DataFrame(
        rows,
        columns=['condition', 'KE', 'Title', 'FDR', 'NES', 'nes_status', 'ES'],
    )

    # H-1: a well-formed enrichment table has unique KEs per condition. If a
    # condition's blob contains duplicate (condition, KE) pairs, pivot_table's
    # aggfunc='first' silently keeps only the first FDR and drops the rest.
    # Warn so the data corruption is visible rather than silent.
    if df.duplicated(['condition', 'KE']).any():
        dup_count = int(df.duplicated(['condition', 'KE']).sum())
        logger.warning(
            'build_comparison_matrix: %d duplicate (condition, KE) row(s) found — '
            'keeping first FDR per pair; check upstream enrichment_json integrity.',
            dup_count,
        )

    # Pivot to KE × condition, keeping the first FDR value where there are
    # duplicate (condition, KE) pairs (should not happen in practice).
    pivot = df.pivot_table(index='KE', columns='condition', values='FDR', aggfunc='first')

    # Reindex columns to upload-position order to prevent pandas alphabetical sort.
    condition_labels = [c.condition_label for c in conditions]
    pivot = pivot.reindex(columns=condition_labels)

    # Compute -log10(FDR) matrix: None for absent or non-significant cells.
    def _neg_log10(fdr_val: Any) -> float | None:
        """Return -log10(fdr_val) or None for missing/non-significant values."""
        if pd.isna(fdr_val) or fdr_val > Config.SIGNIFICANCE_FDR_CUTOFF:
            return None
        return float(-np.log10(max(float(fdr_val), 1e-300)))

    neg_log10_pivot = pivot.map(_neg_log10)

    # NES pivot (#76). Built only for GSEA batches — an ORA enrichment table has
    # no NES column, so the pivot would be all-NaN and only add noise.
    if is_gsea:
        nes_pivot = df.pivot_table(
            index='KE', columns='condition', values='NES', aggfunc='first'
        ).reindex(columns=condition_labels)
        # Reindex rows too: a KE whose NES was NaN in every condition is dropped
        # by pivot_table, and the matrix must stay aligned with `pivot`.
        nes_pivot = nes_pivot.reindex(index=pivot.index)
        # Issue #117 — carried alongside, not folded in. A NES of None on a
        # beyond-resolution cell used to reach compare.html as an em-dash with
        # isSig false, i.e. the same rendering as a Key Event that was never
        # tested, which is the inverse of what happened. pivot_table drops
        # non-numeric columns, so the status goes through pivot() on the
        # de-duplicated frame instead.
        status_source = df.drop_duplicates(['condition', 'KE'])
        nes_status_pivot = status_source.pivot(
            index='KE', columns='condition', values='nes_status'
        ).reindex(columns=condition_labels).reindex(index=pivot.index)
        es_pivot = df.pivot_table(
            index='KE', columns='condition', values='ES', aggfunc='first'
        ).reindex(columns=condition_labels).reindex(index=pivot.index)
    else:
        nes_pivot = None
        nes_status_pivot = None
        es_pivot = None

    # Sort KE rows by mean -log10(FDR) significance descending (most significant first).
    # Fill NaN with 0 for sorting purposes only.
    mean_sig = neg_log10_pivot.fillna(0).mean(axis=1)
    sort_order = mean_sig.sort_values(ascending=False).index
    pivot = pivot.loc[sort_order]
    neg_log10_pivot = neg_log10_pivot.loc[sort_order]
    if nes_pivot is not None:
        nes_pivot = nes_pivot.loc[sort_order]
        nes_status_pivot = nes_status_pivot.loc[sort_order]
        es_pivot = es_pivot.loc[sort_order]

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
    nes_matrix = (
        _to_list(nes_pivot) if nes_pivot is not None
        else [[None] * len(condition_labels) for _ in ke_labels]
    )
    nes_status_matrix = (
        _to_list(nes_status_pivot) if nes_status_pivot is not None
        else [[None] * len(condition_labels) for _ in ke_labels]
    )
    es_matrix = (
        _to_list(es_pivot) if es_pivot is not None
        else [[None] * len(condition_labels) for _ in ke_labels]
    )

    # Assign colours from palette in upload-position order.
    condition_colors = CONDITION_PALETTE[:len(conditions)]

    # Per-condition metadata (aligned to condition_labels / upload-position order).
    # Read straight off the ConditionRecord objects — no extra queries.
    condition_gene_counts = [getattr(c, 'gene_count', None) for c in conditions]
    condition_sig_gene_counts = [getattr(c, 'significant_genes', None) for c in conditions]
    condition_doses = [getattr(c, 'dose', '') or '' for c in conditions]
    condition_timepoints = [getattr(c, 'timepoint', '') or '' for c in conditions]

    # Significant-KE count per condition: number of KEs with FDR below the
    # single app-wide cutoff (issue #63 — same rule as the network borders).
    # Derive from the pivot (columns already reindexed to condition_labels order).
    sig_ke_series = (pivot < Config.SIGNIFICANCE_FDR_CUTOFF).sum(axis=0).reindex(condition_labels).fillna(0)
    condition_sig_ke_counts = [int(v) for v in sig_ke_series]

    return {
        'ke_labels': ke_labels,
        'ke_titles': ke_titles,
        'condition_labels': condition_labels,
        'method': 'gsea' if is_gsea else 'ora',  # #76
        'fdr_matrix': fdr_matrix,
        'neg_log10_matrix': neg_log10_matrix,
        'nes_matrix': nes_matrix,  # #76 — all None for ORA batches
        # #117 — all None for ORA batches and for GSEA batches stored before the
        # diagnostics existed, which the frontend reads as "no claim made".
        'nes_status_matrix': nes_status_matrix,
        'es_matrix': es_matrix,
        'condition_colors': condition_colors,
        'condition_gene_counts': condition_gene_counts,
        'condition_sig_gene_counts': condition_sig_gene_counts,
        'condition_sig_ke_counts': condition_sig_ke_counts,
        'condition_doses': condition_doses,
        'condition_timepoints': condition_timepoints,
    }


# ---------------------------------------------------------------------------
# CSV / Excel export helpers (Feature A + Feature C)
# ---------------------------------------------------------------------------

# Leading characters that spreadsheet apps may interpret as a formula. Cells
# starting with any of these are prefixed with a single quote to neutralise
# CSV-injection (matches the client-side _csvEscape convention in results.html).
_CSV_INJECTION_CHARS = ('=', '+', '-', '@')


def csv_guard(value: Any) -> Any:
    """Neutralise CSV-injection for a single cell value.

    Strings beginning with a formula-trigger character are prefixed with a
    single quote. Non-string values are returned unchanged.

    Args:
        value: The cell value to sanitise.

    Returns:
        The sanitised value (string) or the original value unchanged.
    """
    if isinstance(value, str) and value and value[0] in _CSV_INJECTION_CHARS:
        return "'" + value
    return value


def comparison_matrix_to_dataframe(
    matrix_dict: dict, which: str = 'fdr'
) -> pd.DataFrame:
    """Shape a comparison matrix into a wide, download-ready DataFrame.

    Produces one row per Key Event and one column per condition (upload-position
    order), matching the on-screen comparison table. Two value matrices are
    selectable via ``which``:

      - ``'fdr'``      raw FDR floats (blank where the KE is absent in a condition)
      - ``'neglog10'`` -log10(FDR) (blank where non-significant or absent,
                       matching the heatmap's null semantics)
      - ``'nes'``      normalised enrichment score (#76; blank throughout for
                       an ORA batch, which produces no NES)
      - ``'nes_status'`` the per-cell NES regime (#117), so a blank in the NES
                       export can be told apart from an untested Key Event
      - ``'es'``       raw enrichment score, which carries direction even where
                       the NES could not be normalised (#117)

    The NES export is deliberately **not** masked at the significance cutoff the
    way the heatmap is (#107). A download is the raw record — the FDR matrix is
    downloadable alongside it, so a reader can apply any cutoff they like —
    whereas the heatmap has to choose a colour per cell and must not encode a
    direction the data cannot support.

    Key Event titles are passed through :func:`csv_guard` so a title beginning
    with a formula character cannot be interpreted as a formula on open.

    Args:
        matrix_dict: The dict returned by :func:`build_comparison_matrix`.
        which: ``'fdr'``, ``'neglog10'`` or ``'nes'``.

    Returns:
        A pandas DataFrame with columns ``Key Event ID``, ``Key Event Title``,
        then one column per condition label. Empty if ``matrix_dict`` is empty.
    """
    if not matrix_dict or not matrix_dict.get('ke_labels'):
        return pd.DataFrame()

    matrix_keys = {
        'fdr': 'fdr_matrix',
        'neglog10': 'neg_log10_matrix',
        'nes': 'nes_matrix',
        'nes_status': 'nes_status_matrix',
        'es': 'es_matrix',
    }
    if which not in matrix_keys:
        raise ValueError(
            f"Unknown matrix selector: {which!r} (expected one of "
            f"{', '.join(sorted(matrix_keys))})"
        )

    # nes_matrix is absent from matrices built before #76; treat it as all-blank.
    matrix = matrix_dict.get(matrix_keys[which])
    if matrix is None:
        matrix = [[None] * len(matrix_dict['condition_labels'])
                  for _ in matrix_dict['ke_labels']]
    ke_labels = matrix_dict['ke_labels']
    ke_titles = matrix_dict['ke_titles']
    condition_labels = matrix_dict['condition_labels']

    data = {
        'Key Event ID': [csv_guard(k) for k in ke_labels],
        'Key Event Title': [csv_guard(t) for t in ke_titles],
    }
    for col_idx, cond_label in enumerate(condition_labels):
        # None cells become NaN -> render as empty string in CSV / empty in Excel.
        data[cond_label] = [row[col_idx] for row in matrix]

    return pd.DataFrame(data)


#: Issue #117 — how a beyond-resolution cell reads in a report table. It has no
#: NES, so the cell carries the direction from the ES and says why there is no
#: number. Blank would be indistinguishable from an untested Key Event, which is
#: the opposite result.
BEYOND_RESOLUTION_UP = '▲ beyond res.'
BEYOND_RESOLUTION_DOWN = '▼ beyond res.'
BEYOND_RESOLUTION_UNKNOWN = '? beyond res.'
#: ASCII forms for ReportLab, whose Helvetica has no arrow glyphs and would
#: print a black box where the direction should be.
BEYOND_RESOLUTION_UP_ASCII = '+ beyond res.'
BEYOND_RESOLUTION_DOWN_ASCII = '- beyond res.'


def comparison_nes_display_dataframe(
    matrix_dict: dict, ascii_only: bool = False
) -> pd.DataFrame:
    """Shape the NES comparison matrix into report-ready **text** (#117).

    The numeric NES export (``which='nes'``) is the machine-readable record and
    stays numeric, so a beyond-resolution cell is blank in it — indistinguishable
    from a Key Event that was never tested. A report is read by a person, and
    both batch reports used to caption that table "every NES value, significant
    or not", which was false in exactly the direction that hides the strongest
    calls in the grid.

    This renders the same matrix as strings: a formatted NES where one exists, a
    direction-carrying marker where it does not, and a trailing marker on cells
    whose same-signed permutation tail was too short for either the magnitude or
    the nominal p-value to be fine-grained.

    Args:
        matrix_dict: The dict returned by :func:`build_comparison_matrix`.
        ascii_only: Use ``+``/``-`` instead of ``▲``/``▼`` for the direction
            marker. ReportLab's default Helvetica has no arrow glyphs and
            renders them as black boxes, so the PDF report passes True.

    Returns:
        A DataFrame with ``Key Event ID``, ``Key Event Title`` and one text
        column per condition. Empty if ``matrix_dict`` is empty.
    """
    if not matrix_dict or not matrix_dict.get('ke_labels'):
        return pd.DataFrame()

    ke_labels = matrix_dict['ke_labels']
    ke_titles = matrix_dict['ke_titles']
    condition_labels = matrix_dict['condition_labels']
    blank = [[None] * len(condition_labels) for _ in ke_labels]
    nes = matrix_dict.get('nes_matrix') or blank
    status = matrix_dict.get('nes_status_matrix') or blank
    es = matrix_dict.get('es_matrix') or blank

    def _cell(row: int, col: int) -> str:
        cell_status = status[row][col] if status[row] else None
        if cell_status == 'beyond_permutation_resolution':
            es_val = es[row][col] if es[row] else None
            if es_val is None:
                return BEYOND_RESOLUTION_UNKNOWN
            if es_val >= 0:
                return (BEYOND_RESOLUTION_UP_ASCII if ascii_only
                        else BEYOND_RESOLUTION_UP)
            return (BEYOND_RESOLUTION_DOWN_ASCII if ascii_only
                    else BEYOND_RESOLUTION_DOWN)
        value = nes[row][col] if nes[row] else None
        if value is None or (isinstance(value, float) and value != value):
            return ''
        text = f'{float(value):+.2f}'
        if cell_status == 'unstable_normalisation':
            return f'{text} !'
        if cell_status == 'undiagnosed':
            return f'{text} ?'
        return text

    data = {
        'Key Event ID': [csv_guard(k) for k in ke_labels],
        'Key Event Title': [csv_guard(t) for t in ke_titles],
    }
    for col_idx, cond_label in enumerate(condition_labels):
        data[cond_label] = [_cell(row, col_idx) for row in range(len(ke_labels))]
    return pd.DataFrame(data)


def build_gene_tracking(conditions: list) -> dict[str, Any]:
    """Track which genes drive each KE across batch conditions.

    A gene is a *driver* of a KE in a condition when it is present in that KE's
    reference set **and** significant in that condition (i.e. it appears in the
    enrichment overlap). The per-gene significance/log2FC is read from each
    condition's ``ke_gene_json`` (the richest source — carries log2FC and
    significance per gene per condition).

    To bound memory (H-3), only genes that are a driver in at least one
    condition are retained — the full KE membership is never materialised.

    Args:
        conditions: List of ConditionRecord objects ordered by position. Each
                    must expose ``condition_label``, ``ke_gene_json`` and
                    (optionally) ``enrichment_json`` for KE titles.

    Returns:
        Dict with keys:
            condition_labels — condition labels in upload-position order
            records          — tidy list, one dict per (KE, gene, condition)
                               driver event: KE_ID, KE_Title, Gene_Symbol,
                               Condition, log2FC, significant. Ordered by KE
                               (numerically, so KE:10 follows KE:9), then gene
                               symbol, then condition upload position (issue
                               #82) — the stored ``ke_gene_json`` is a JSON
                               object whose key order is an artefact of how it
                               was written, so an unordered list made both the
                               Genes tab payload and the batch report vary
                               between runs over identical data.
            gene_ke_summary  — derived list, one dict per (KE, gene) that is a
                               driver somewhere: KE_ID, KE_Title, Gene_Symbol,
                               n_conditions_driver, conditions_driver (list),
                               shared (bool, driver in >=2 conditions),
                               log2FC_by_condition (list aligned to condition order)
    """
    condition_labels = [c.condition_label for c in conditions]
    n_conditions = len(conditions)

    # KE title lookup, first-seen across conditions' enrichment tables.
    ke_title_map: dict[str, str] = {}
    for cond in conditions:
        if not getattr(cond, 'enrichment_json', None):
            continue
        try:
            for entry in json.loads(cond.enrichment_json):
                ke = entry.get('KE')
                if ke and ke not in ke_title_map:
                    ke_title_map[ke] = entry.get('Title', '')
        except (json.JSONDecodeError, TypeError):
            continue

    def _ke_sort_key(value: Any) -> tuple:
        """Order KE IDs by their trailing integer, so KE:10 follows KE:9.

        Kept local (rather than shared with the export route) so the ordering
        contract of this function is readable in one place; IDs with no numeric
        suffix sort before numbered ones within the same prefix.
        """
        text = str(value)
        cut = len(text)
        while cut > 0 and text[cut - 1].isdigit():
            cut -= 1
        if cut < len(text):
            return (text[:cut], int(text[cut:]))
        return (text, -1)

    # (sort key, record) pairs; the trailing counter makes the key total so the
    # sort never falls back to comparing dicts and never depends on input order.
    indexed_records: list[tuple[tuple, dict]] = []
    # summary keyed by (KE, gene) -> aggregation dict
    summary: dict[tuple, dict] = {}

    for col_idx, cond in enumerate(conditions):
        if not getattr(cond, 'ke_gene_json', None):
            continue
        try:
            ke_gene_map = json.loads(cond.ke_gene_json)
        except (json.JSONDecodeError, TypeError) as exc:
            logger.warning(
                'build_gene_tracking: failed to parse ke_gene_json for condition %s: %s',
                cond.condition_label, exc,
            )
            continue

        for ke, gene_entries in ke_gene_map.items():
            for entry in gene_entries:
                # Driver = significant in this condition (bounds the output).
                if not entry.get('significant'):
                    continue
                gene = entry.get('id')
                if not gene:
                    continue
                log2fc = entry.get('log2FC')
                title = ke_title_map.get(ke, ke)

                indexed_records.append((
                    (_ke_sort_key(ke), str(gene), col_idx, len(indexed_records)),
                    {
                        'KE_ID': ke,
                        'KE_Title': title,
                        'Gene_Symbol': gene,
                        'Condition': cond.condition_label,
                        'log2FC': log2fc,
                        'significant': True,
                    },
                ))

                key = (ke, gene)
                agg = summary.get(key)
                if agg is None:
                    agg = {
                        'KE_ID': ke,
                        'KE_Title': title,
                        'Gene_Symbol': gene,
                        'conditions_driver': [],
                        'log2FC_by_condition': [None] * n_conditions,
                    }
                    summary[key] = agg
                agg['conditions_driver'].append(cond.condition_label)
                agg['log2FC_by_condition'][col_idx] = log2fc

    gene_ke_summary = []
    for agg in summary.values():
        n_driver = len(agg['conditions_driver'])
        gene_ke_summary.append({
            'KE_ID': agg['KE_ID'],
            'KE_Title': agg['KE_Title'],
            'Gene_Symbol': agg['Gene_Symbol'],
            'n_conditions_driver': n_driver,
            'conditions_driver': agg['conditions_driver'],
            'shared': n_driver >= 2,
            'log2FC_by_condition': agg['log2FC_by_condition'],
        })

    # Stable, useful ordering: KE (numerically), then shared-first, then gene.
    gene_ke_summary.sort(
        key=lambda r: (
            _ke_sort_key(r['KE_ID']), -r['n_conditions_driver'], str(r['Gene_Symbol'])
        )
    )

    records = [record for _, record in sorted(indexed_records, key=lambda pair: pair[0])]

    return {
        'condition_labels': condition_labels,
        'records': records,
        'gene_ke_summary': gene_ke_summary,
    }


def gene_tracking_to_dataframes(tracking: dict) -> dict[str, pd.DataFrame]:
    """Shape gene-tracking output into download-ready DataFrames.

    Args:
        tracking: The dict returned by :func:`build_gene_tracking`.

    Returns:
        Dict with two DataFrames:
            'long'    — one row per (KE, gene, condition) driver event.
            'summary' — one row per (KE, gene) with a shared flag and one
                        log2FC column per condition (wide).
        Both have CSV-injection guards applied to free-text columns.
    """
    condition_labels = tracking.get('condition_labels', [])

    # Long / tidy table.
    long_rows = []
    for r in tracking.get('records', []):
        long_rows.append({
            'KE_ID': csv_guard(r['KE_ID']),
            'KE_Title': csv_guard(r['KE_Title']),
            'Gene_Symbol': csv_guard(r['Gene_Symbol']),
            'Condition': csv_guard(r['Condition']),
            'log2FC': r['log2FC'],
            'significant': r['significant'],
        })
    long_df = pd.DataFrame(
        long_rows,
        columns=['KE_ID', 'KE_Title', 'Gene_Symbol', 'Condition', 'log2FC', 'significant'],
    )

    # Summary / wide table: one log2FC column per condition.
    summary_rows = []
    for r in tracking.get('gene_ke_summary', []):
        row = {
            'KE_ID': csv_guard(r['KE_ID']),
            'KE_Title': csv_guard(r['KE_Title']),
            'Gene_Symbol': csv_guard(r['Gene_Symbol']),
            'n_conditions_driver': r['n_conditions_driver'],
            'shared': r['shared'],
        }
        for idx, cond_label in enumerate(condition_labels):
            row[f'log2FC: {cond_label}'] = r['log2FC_by_condition'][idx]
        summary_rows.append(row)

    summary_cols = ['KE_ID', 'KE_Title', 'Gene_Symbol', 'n_conditions_driver', 'shared'] + [
        f'log2FC: {c}' for c in condition_labels
    ]
    summary_df = pd.DataFrame(summary_rows, columns=summary_cols)

    return {'long': long_df, 'summary': summary_df}
