"""
Batch analysis service for multi-file enrichment runs.

Provides upload directory management, cisplatin filename parsing, demo file listing,
and core batch orchestration (background harmonisation, per-condition enrichment loop).
"""

import json
import logging
import math
import os
import re
import shutil
import uuid
from datetime import datetime, timezone
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

# Regex pattern for cisplatin kidney toxicity filenames.
# Matches patterns like: CSP_4hr_0.1uM.csv, CSP_72hr_50uM.csv
CISPLATIN_RE = re.compile(r'CSP_(\d+hr)_([\d.]+uM)\.csv$')


def parse_cisplatin_filename(filename: str) -> Dict[str, str]:
    """Extract dose, timepoint, and condition label from a cisplatin filename.

    Matches the naming convention used by the Cisplatin_Kidney demo dataset:
    CSP_{timepoint}_{dose}.csv (e.g. CSP_4hr_0.1uM.csv).

    Args:
        filename: The bare filename (not a full path) to parse.

    Returns:
        Dict with keys:
          - 'condition_label': Human-readable label (e.g. 'Cisplatin 4hr 0.1uM')
          - 'dose': Dose string (e.g. '0.1uM'), or '' if no match
          - 'timepoint': Timepoint string (e.g. '4hr'), or '' if no match
    """
    match = CISPLATIN_RE.search(filename)
    if match:
        timepoint = match.group(1)
        dose = match.group(2)
        return {
            'condition_label': f'Cisplatin {timepoint} {dose}',
            'dose': dose,
            'timepoint': timepoint,
        }
    # Non-matching filename — use filename as-is for the condition label
    return {
        'condition_label': filename,
        'dose': '',
        'timepoint': '',
    }


def create_batch_upload_dir(upload_folder: str) -> Tuple[str, str]:
    """Create a UUID-scoped upload directory for a new batch.

    Each batch gets its own subdirectory under uploads/ to prevent filename
    collisions across concurrent batches and to simplify cleanup.

    Args:
        upload_folder: Base upload directory path (typically Config.UPLOAD_FOLDER).

    Returns:
        Tuple of (batch_uuid, batch_dir_path) where batch_uuid is a new UUID
        string and batch_dir_path is the absolute/relative path to the created
        directory.
    """
    batch_uuid = str(uuid.uuid4())
    batch_dir = os.path.join(upload_folder, batch_uuid)
    os.makedirs(batch_dir, exist_ok=True)
    logger.info(f"Created batch upload directory: {batch_dir}")
    return batch_uuid, batch_dir


def cleanup_batch_upload_dir(upload_folder: str, batch_uuid: str) -> None:
    """Remove a batch's UUID-scoped upload directory and all its contents.

    Safe to call even if the directory does not exist (no-op in that case).

    Args:
        upload_folder: Base upload directory path (typically Config.UPLOAD_FOLDER).
        batch_uuid: UUID string identifying the batch whose directory should
                    be removed.
    """
    batch_dir = os.path.join(upload_folder, batch_uuid)
    if os.path.isdir(batch_dir):
        shutil.rmtree(batch_dir)
        logger.info(f"Cleaned up batch upload directory: {batch_dir}")
    else:
        logger.debug(f"Batch upload directory not found (already removed?): {batch_dir}")


def get_cisplatin_demo_files() -> List[Dict]:
    """Return available cisplatin demo files grouped by timepoint for the UI.

    Scans the data/Cisplatin_Kidney/ directory using Config.CISPLATIN_TIMEPOINTS
    and Config.CISPLATIN_DOSES to enumerate expected files. Only files that
    actually exist on disk are included.

    Returns:
        List of timepoint group dicts, each containing:
          - 'timepoint': Timepoint label string (e.g. '4hr')
          - 'files': List of dicts with 'path', 'label', 'dose' keys
            where 'path' is relative to data/ (e.g. 'Cisplatin_Kidney/CSP_4hr_0.1uM.csv')
    """
    # Import here to avoid circular imports at module load time
    from config import Config

    groups = []
    for timepoint in Config.CISPLATIN_TIMEPOINTS:
        files = []
        for dose in Config.CISPLATIN_DOSES:
            filename = f'CSP_{timepoint}_{dose}.csv'
            relative_path = f'Cisplatin_Kidney/{filename}'
            data_path = os.path.join('data', relative_path)
            if os.path.isfile(data_path):
                files.append({
                    'path': relative_path,
                    'label': filename,
                    'dose': dose,
                })
        if files:
            groups.append({
                'timepoint': timepoint,
                'files': files,
            })

    return groups


# ---------------------------------------------------------------------------
# Batch orchestration
# ---------------------------------------------------------------------------

def validate_batch_columns(file_infos: List[Dict]) -> Tuple[bool, str]:
    """Validate that all files in a batch share the same column mapping.

    Each entry in file_infos must have keys: ``filepath``, ``id_col``,
    ``fc_col``, ``pval_col``.  The function reads the first five rows of each
    file and checks that all named columns are present, and that the column
    names are consistent across all files.

    Args:
        file_infos: List of dicts with filepath and column name keys.

    Returns:
        ``(True, '')`` on success.
        ``(False, error_message)`` when validation fails.
    """
    if not file_infos:
        return False, 'No files provided for validation.'

    first_id = file_infos[0].get('id_col', '')
    first_fc = file_infos[0].get('fc_col', '')
    first_pval = file_infos[0].get('pval_col', '')

    for info in file_infos:
        filepath = info.get('filepath', '')
        id_col = info.get('id_col', '')
        fc_col = info.get('fc_col', '')
        pval_col = info.get('pval_col', '')

        # Column names must be consistent across all files
        if id_col != first_id or fc_col != first_fc or pval_col != first_pval:
            return (
                False,
                'Column structure mismatch across files — please ensure all '
                'files use the same column layout.',
            )

        try:
            df_head = pd.read_csv(filepath, sep=None, engine='python', nrows=5)
        except Exception as exc:
            return False, f'Could not read file {os.path.basename(filepath)}: {exc}'

        missing = [
            col for col in [id_col, fc_col, pval_col]
            if col not in df_head.columns
        ]
        if missing:
            return (
                False,
                f'File {os.path.basename(filepath)} is missing columns: '
                f'{", ".join(missing)}.',
            )

    return True, ''


def harmonise_backgrounds(
    file_infos: List[Dict],
) -> Tuple[Set[str], Dict[str, int]]:
    """Compute the intersection of gene universes across all batch files.

    Each entry in file_infos must have keys: ``filepath``, ``id_col``,
    ``fc_col``, ``pval_col``.

    The gene universe for each file is the full set of gene symbols (uppercased)
    in the 'ID' column returned by ``load_and_validate_data``.  The harmonised
    background is the intersection of all per-file universes, ensuring Fisher's
    exact test uses only genes measured in every condition.

    Args:
        file_infos: List of dicts describing each file and its column mapping.

    Returns:
        Tuple of:
          - ``harmonised_genes``: Set of gene symbols present in every file.
          - ``per_file_counts``: Dict mapping filepath to the raw gene count
            for that file (before intersection).
    """
    from services.data_service import load_and_validate_data  # avoid circular import at module level

    universes: List[Set[str]] = []
    per_file_counts: Dict[str, int] = {}

    for info in file_infos:
        filepath = info['filepath']
        id_col = info['id_col']
        fc_col = info['fc_col']
        pval_col = info['pval_col']

        df = load_and_validate_data(filepath, id_col, fc_col, pval_col)
        gene_universe = set(df['ID'].str.upper())
        universes.append(gene_universe)
        per_file_counts[filepath] = len(gene_universe)
        logger.info(
            f'harmonise_backgrounds: {os.path.basename(filepath)} — '
            f'{len(gene_universe)} unique genes'
        )

    if not universes:
        return set(), {}

    harmonised = universes[0]
    for universe in universes[1:]:
        harmonised = harmonised & universe

    logger.info(
        f'harmonise_backgrounds: intersection = {len(harmonised)} genes '
        f'across {len(file_infos)} files'
    )
    return harmonised, per_file_counts


def _run_condition(
    cond,
    batch,
    harmonised_genes: Set[str],
    reference_sets: Dict,
    ke_list,
    edges,
    ke_type_map: Dict[str, str],
    ke_title_map: Dict[str, str],
    db_session,
) -> None:
    """Run enrichment analysis for one ConditionRecord within a batch.

    Loads the file, applies the harmonised background filter, runs enrichment
    and network building, stores results on the ConditionRecord, and commits.

    Args:
        cond: ConditionRecord ORM object to process.
        batch: Parent BatchRecord (provides shared parameters).
        harmonised_genes: Intersection gene set used to filter each file.
        reference_sets: KE-to-gene reference dict from load_cached_reference_sets.
        ke_list: Set of KE IDs for the selected AOP.
        edges: DataFrame of KE-to-KE edges for the selected AOP.
        ke_type_map: Dict mapping KE IDs to type strings.
        ke_title_map: Dict mapping KE IDs to title strings.
        db_session: SQLAlchemy session to use for commits.
    """
    from services.data_service import load_and_validate_data, process_gene_expression
    from services.enrichment_service import run_enrichment_analysis, build_ke_gene_mapping
    from services.network_service import build_cytoscape_network

    cond.status = 'running'
    cond.started_at = datetime.now(timezone.utc).replace(tzinfo=None)
    db_session.commit()

    # Build filepath from the batch upload directory
    from config import Config
    filepath = os.path.join(Config.UPLOAD_FOLDER, batch.uuid, cond.filename)

    # Load raw data using the batch-level column mapping
    df_raw = load_and_validate_data(
        filepath,
        batch.id_column,
        batch.fc_column,
        batch.pval_column,
    )

    # Process gene expression (applies significance flags)
    df_processed, _ = process_gene_expression(df_raw, batch.logfc_threshold or 0.0)

    # Apply harmonisation: restrict to genes present in all batch files
    df_filtered = df_processed[df_processed['ID'].isin(harmonised_genes)].copy()

    # Recalculate gene counts after filtering
    total_genes = len(df_filtered)
    significant_genes = int(df_filtered['significant'].sum())

    # Run enrichment on the harmonised/filtered data
    enrichment_results = run_enrichment_analysis(
        df_filtered, reference_sets, ke_list, ke_title_map
    )

    # Build Cytoscape network
    cy_network = build_cytoscape_network(
        ke_list, edges, enrichment_results, ke_title_map, ke_type_map,
        reference_sets=reference_sets,
    )

    # Build KE gene membership maps for the drawer / results page
    gene_logfc_map = df_filtered.set_index('ID')['log2FC'].to_dict()
    gene_significance_map = df_filtered.set_index('ID')['significant'].to_dict()
    ke_gene_map = build_ke_gene_mapping(
        reference_sets, ke_list, gene_logfc_map, gene_significance_map
    )

    # Convert enrichment DataFrame to JSON-serialisable list
    enrichment_table = enrichment_results.to_dict(orient='records')
    for row in enrichment_table:
        for key, value in row.items():
            if isinstance(value, float) and (value != value or math.isinf(value)):  # NaN/inf
                row[key] = None
            elif hasattr(value, 'item'):  # numpy scalar
                row[key] = value.item()

    # Persist results on the ConditionRecord
    cond.enrichment_json = json.dumps(enrichment_table)
    cond.network_json = json.dumps(cy_network)
    cond.ke_gene_json = json.dumps(ke_gene_map)
    cond.ke_type_map_json = json.dumps(ke_type_map)
    cond.ke_title_map_json = json.dumps(ke_title_map)
    cond.gene_count = total_genes
    cond.significant_genes = significant_genes
    cond.status = 'complete'
    cond.completed_at = datetime.now(timezone.utc).replace(tzinfo=None)
    db_session.commit()

    logger.info(
        f'_run_condition: condition {cond.id} ({cond.condition_label}) complete — '
        f'{total_genes} genes, {significant_genes} significant'
    )


def run_batch(batch_id: int, db_url: str, reference_sets: Dict) -> None:
    """Orchestrate a full batch run: process each ConditionRecord sequentially.

    Creates its own SQLAlchemy session from ``db_url`` so it is safe to call
    from a background thread (SQLAlchemy sessions are not thread-safe).

    Sets batch.status to 'running', iterates over conditions in position order,
    calls ``_run_condition`` for each, then marks the batch 'complete'.
    On any unhandled exception in a condition, stops the batch and marks it
    'failed'.  Always cleans up the batch upload directory on exit.

    Args:
        batch_id: Primary key of the BatchRecord to process.
        db_url: SQLAlchemy database URL used to create a new session.
        reference_sets: KE-to-gene reference sets loaded by the app at startup.
    """
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker
    from sqlalchemy.pool import StaticPool
    from database import BatchRecord, ConditionRecord
    from services.data_service import load_aop_data
    from config import Config

    engine_kwargs = dict(echo=False)
    if 'sqlite' in db_url:
        engine_kwargs.update(poolclass=StaticPool, connect_args={"check_same_thread": False})
    else:
        engine_kwargs.update(pool_pre_ping=True, pool_size=5, max_overflow=10)
    engine = create_engine(db_url, **engine_kwargs)
    Session = sessionmaker(bind=engine)
    db_session = Session()

    try:
        batch = db_session.query(BatchRecord).filter_by(id=batch_id).first()
        if batch is None:
            logger.error(f'run_batch: BatchRecord {batch_id} not found.')
            return

        batch.status = 'running'
        db_session.commit()
        logger.info(f'run_batch: starting batch {batch.uuid} ({len(batch.conditions)} conditions)')

        # Load harmonised background from the stored JSON
        harmonised_genes: Set[str] = set(
            json.loads(batch.harmonised_background or '[]')
        )

        # Load AOP topology once for the whole batch
        ke_list, edges, ke_type_map, ke_title_map = load_aop_data(batch.aop_id)

        for cond in sorted(batch.conditions, key=lambda c: c.position):
            # Allow cooperative cancellation between conditions
            db_session.refresh(batch)
            if batch.status == 'cancelled':
                logger.info(f'run_batch: batch {batch.uuid} cancelled, stopping.')
                break

            try:
                _run_condition(
                    cond,
                    batch,
                    harmonised_genes,
                    reference_sets,
                    ke_list,
                    edges,
                    ke_type_map,
                    ke_title_map,
                    db_session,
                )
            except Exception as exc:
                logger.error(
                    f'run_batch: condition {cond.id} ({cond.condition_label}) failed: {exc}',
                    exc_info=True,
                )
                cond.status = 'failed'
                cond.error_message = str(exc)
                batch.status = 'failed'
                db_session.commit()
                return  # Stop-batch strategy on first failure

        # Mark complete if not already failed/cancelled
        if batch.status == 'running':
            batch.status = 'complete'
            batch.completed_at = datetime.now(timezone.utc).replace(tzinfo=None)
            db_session.commit()
            logger.info(f'run_batch: batch {batch.uuid} complete.')

    except Exception as exc:
        logger.error(f'run_batch: unexpected error for batch {batch_id}: {exc}', exc_info=True)
        try:
            db_session.rollback()
        except Exception:
            pass
    finally:
        db_session.close()
        # Clean up uploaded files regardless of outcome
        try:
            batch_obj = None
            with Session() as s:
                batch_obj = s.query(BatchRecord).filter_by(id=batch_id).first()
                batch_uuid = batch_obj.uuid if batch_obj else None
            if batch_uuid:
                cleanup_batch_upload_dir(Config.UPLOAD_FOLDER, batch_uuid)
        except Exception as cleanup_exc:
            logger.warning(f'run_batch: cleanup failed for batch {batch_id}: {cleanup_exc}')
