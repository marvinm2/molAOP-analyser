"""
Configuration settings for the Molecular AOP Analyser application.
"""
import os
from dataclasses import dataclass
from typing import Optional
from datetime import datetime

class Config:
    # File upload settings
    MAX_FILE_SIZE = 10 * 1024 * 1024  # 10 MB
    MAX_CONTENT_LENGTH = 50 * 1024 * 1024  # 50 MB for form data (reports)
    ALLOWED_EXTENSIONS = {'csv', 'tsv', 'txt'}
    UPLOAD_FOLDER = 'uploads'
    
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
        'GSE90122_SR12813.tsv': 'PXR agonist 2 – GSE90122_SR12813'
    }
    
    # AOP case studies
    CASE_STUDY_AOPS = {
        "DEMO": {"label": "---DEMO---", "enabled": False},
        "steatosis": {"id": "AOP:1", "label": "PXR activation leading to liver steatosis", "enabled": True},
        "VHP-CASES:": {"label": "---VHP CASES---", "enabled": False},
        "vhp-kidney": {"id": "AOP:2", "label": "DNA adduct formation leading to kidney failure", "enabled": False},
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
            self.upload_timestamp = datetime.now()
    
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