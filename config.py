"""
Configuration settings for the Molecular AOP Analyser application.
"""
import os
from dataclasses import dataclass
from typing import Optional
from datetime import datetime, timezone
import zoneinfo

class Config:
    # File upload settings
    MAX_FILE_SIZE = 10 * 1024 * 1024  # 10 MB
    MAX_CONTENT_LENGTH = 50 * 1024 * 1024  # 50 MB for form data (reports)
    ALLOWED_EXTENSIONS = {'csv', 'tsv', 'txt'}
    UPLOAD_FOLDER = 'uploads'
    TEMP_FOLDER = 'temp'  # Temporary files for network PNGs
    
    # Data processing settings
    MAX_GENES_DISPLAY = 10000
    PVAL_CUTOFF = 0.05
    
    # Required data files
    REQUIRED_DATA_FILES = [
        'data/aop_ke_map.csv',
        'data/aop_ker_edges.csv', 
        'data/KE-WP.csv',
        'data/edges_wpid_to_gene.csv',
        'data/node_attributes.csv',
        'data/ke_metadata.csv'
    ]
    
    # Demo datasets
    DEMO_DATASETS = {
        'GSE90122_TO90137.tsv': 'PXR agonist 1 – GSE90122_TO90137',
        'GSE90122_SR12813.tsv': 'PXR agonist 2 – GSE90122_SR12813',
        'Cisplatin_Kidney/CSP_4hr_0.1uM.csv': 'Cisplatin 4hr 0.1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_4hr_0.5uM.csv': 'Cisplatin 4hr 0.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_4hr_1uM.csv': 'Cisplatin 4hr 1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_4hr_2.5uM.csv': 'Cisplatin 4hr 2.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_4hr_5uM.csv': 'Cisplatin 4hr 5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_4hr_10uM.csv': 'Cisplatin 4hr 10μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_4hr_20uM.csv': 'Cisplatin 4hr 20μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_4hr_30uM.csv': 'Cisplatin 4hr 30μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_4hr_50uM.csv': 'Cisplatin 4hr 50μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_8hr_0.1uM.csv': 'Cisplatin 8hr 0.1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_8hr_0.5uM.csv': 'Cisplatin 8hr 0.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_8hr_1uM.csv': 'Cisplatin 8hr 1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_8hr_2.5uM.csv': 'Cisplatin 8hr 2.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_8hr_5uM.csv': 'Cisplatin 8hr 5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_8hr_10uM.csv': 'Cisplatin 8hr 10μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_8hr_20uM.csv': 'Cisplatin 8hr 20μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_8hr_30uM.csv': 'Cisplatin 8hr 30μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_8hr_50uM.csv': 'Cisplatin 8hr 50μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_16hr_0.1uM.csv': 'Cisplatin 16hr 0.1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_16hr_0.5uM.csv': 'Cisplatin 16hr 0.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_16hr_1uM.csv': 'Cisplatin 16hr 1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_16hr_2.5uM.csv': 'Cisplatin 16hr 2.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_16hr_5uM.csv': 'Cisplatin 16hr 5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_16hr_10uM.csv': 'Cisplatin 16hr 10μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_16hr_20uM.csv': 'Cisplatin 16hr 20μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_16hr_30uM.csv': 'Cisplatin 16hr 30μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_16hr_50uM.csv': 'Cisplatin 16hr 50μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_24hr_0.1uM.csv': 'Cisplatin 24hr 0.1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_24hr_0.5uM.csv': 'Cisplatin 24hr 0.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_24hr_1uM.csv': 'Cisplatin 24hr 1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_24hr_2.5uM.csv': 'Cisplatin 24hr 2.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_24hr_5uM.csv': 'Cisplatin 24hr 5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_24hr_10uM.csv': 'Cisplatin 24hr 10μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_24hr_20uM.csv': 'Cisplatin 24hr 20μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_24hr_30uM.csv': 'Cisplatin 24hr 30μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_24hr_50uM.csv': 'Cisplatin 24hr 50μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_48hr_0.1uM.csv': 'Cisplatin 48hr 0.1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_48hr_0.5uM.csv': 'Cisplatin 48hr 0.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_48hr_1uM.csv': 'Cisplatin 48hr 1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_48hr_2.5uM.csv': 'Cisplatin 48hr 2.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_48hr_5uM.csv': 'Cisplatin 48hr 5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_48hr_10uM.csv': 'Cisplatin 48hr 10μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_48hr_20uM.csv': 'Cisplatin 48hr 20μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_48hr_30uM.csv': 'Cisplatin 48hr 30μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_48hr_50uM.csv': 'Cisplatin 48hr 50μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_72hr_0.1uM.csv': 'Cisplatin 72hr 0.1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_72hr_0.5uM.csv': 'Cisplatin 72hr 0.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_72hr_1uM.csv': 'Cisplatin 72hr 1μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_72hr_2.5uM.csv': 'Cisplatin 72hr 2.5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_72hr_5uM.csv': 'Cisplatin 72hr 5μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_72hr_10uM.csv': 'Cisplatin 72hr 10μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_72hr_20uM.csv': 'Cisplatin 72hr 20μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_72hr_30uM.csv': 'Cisplatin 72hr 30μM – Kidney toxicity',
        'Cisplatin_Kidney/CSP_72hr_50uM.csv': 'Cisplatin 72hr 50μM – Kidney toxicity'
    }
    
    # AOP case studies
    CASE_STUDY_AOPS = {
        "DEMO": {"label": "---DEMO---", "enabled": False},
        "steatosis": {"id": "AOP:1", "label": "PXR activation leading to liver steatosis", "enabled": True},
        "VHP-CASES:": {"label": "---VHP CASES---", "enabled": False},
        "vhp-kidney": {"id": "AOP:2", "label": "DNA adduct formation leading to kidney failure", "enabled": True},
        "vhp-parkinson": {"id": "AOP:3", "label": "Calcium overload in dopaminergic neurons of the substantia nigra leading to parkinsonian motor deficits", "enabled": False},
        "vhp-thyroid": {"id": "AOP:4", "label": "Thyroid hormone-mediated neurodevelopmental toxicity", "enabled": False},
        "ORGAN-NETWORK": {"label": "---ORGAN NETWORKS---","enabled": False},
        "Liver-aop-network": {"id": "AOP:5", "label": "Liver AOP network", "enabled": False},
        "Brain-aop-network": {"id": "AOP:6", "label": "Brain AOP network", "enabled": False},
        "Kidney-aop-network": {"id": "AOP:7", "label": "Kidney AOP network", "enabled": False},
        "Lung-aop-network": {"id": "AOP:8", "label": "Lung AOP network", "enabled": False},
    }
    
    # Flask settings
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'dev-key-change-in-production'

    # Builder API settings
    BUILDER_API_URL = os.environ.get('BUILDER_API_URL', '')
    BUILDER_API_TIMEOUT = int(os.environ.get('BUILDER_API_TIMEOUT', '30'))
    CACHE_DIR = os.environ.get('CACHE_DIR', '/tmp/molaop_cache')
    CACHE_TTL = 3600  # 1 hour

    @classmethod
    def validate_data_files(cls):
        """Validate that all required data files exist."""
        missing_files = []
        for file_path in cls.REQUIRED_DATA_FILES:
            if not os.path.exists(file_path):
                missing_files.append(file_path)
        
        if missing_files:
            raise FileNotFoundError(f"Required data files missing: {', '.join(missing_files)}")
        
        return True


@dataclass
class ExperimentMetadata:
    """Dataclass for storing experiment metadata for reports."""
    dataset_id: str
    stressor: str
    dosing: str
    owner: str
    upload_timestamp: Optional[datetime] = None
    filename: Optional[str] = None
    description: Optional[str] = None
    
    def __post_init__(self):
        """Set upload timestamp if not provided."""
        if self.upload_timestamp is None:
            # Use Europe/Amsterdam timezone (CEST/CET)
            tz = zoneinfo.ZoneInfo("Europe/Amsterdam")
            self.upload_timestamp = datetime.now(tz)
    
    def to_dict(self) -> dict:
        """Convert metadata to dictionary for JSON serialization."""
        return {
            'dataset_id': self.dataset_id,
            'stressor': self.stressor,
            'dosing': self.dosing,
            'owner': self.owner,
            'upload_timestamp': self.upload_timestamp.isoformat() if self.upload_timestamp else None,
            'filename': self.filename,
            'description': self.description
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> 'ExperimentMetadata':
        """Create metadata from dictionary."""
        # Convert timestamp string back to datetime
        if data.get('upload_timestamp'):
            data['upload_timestamp'] = datetime.fromisoformat(data['upload_timestamp'])
        return cls(**data)