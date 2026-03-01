"""
Batch analysis service for multi-file enrichment runs.

Provides upload directory management, cisplatin filename parsing, and
demo file listing utilities used by batch analysis routes and orchestration.
"""

import os
import re
import shutil
import uuid
import logging
from typing import Dict, List, Optional, Tuple

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
