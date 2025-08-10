"""
Database models and setup for the Molecular AOP Analyser.

Provides SQLite-based persistence for experiment metadata and analysis results.
"""

from datetime import datetime
from typing import Optional, Dict, Any
import json
from sqlalchemy import create_engine, Column, Integer, String, DateTime, Text, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import SQLAlchemyError
import logging

logger = logging.getLogger(__name__)

Base = declarative_base()

class ExperimentRecord(Base):
    """Database model for storing experiment metadata and analysis results."""
    
    __tablename__ = 'experiments'
    
    # Primary key
    id = Column(Integer, primary_key=True, autoincrement=True)
    
    # Experiment metadata
    dataset_id = Column(String(255), nullable=False, index=True)
    filename = Column(String(255), nullable=True)  # Allow null for metadata-only records
    stressor = Column(String(500))
    dosing = Column(String(500))
    owner = Column(String(255))
    description = Column(Text)
    
    # Analysis parameters
    aop_id = Column(String(100))
    logfc_threshold = Column(Float)
    pval_cutoff = Column(Float)
    id_column = Column(String(100))
    fc_column = Column(String(100))
    pval_column = Column(String(100))
    
    # Analysis results summary (JSON)
    enrichment_results = Column(Text)  # JSON string of enrichment results
    gene_count = Column(Integer)
    significant_genes = Column(Integer)
    
    # Timestamps  
    upload_timestamp = Column(DateTime, default=datetime.utcnow, nullable=True)
    analysis_timestamp = Column(DateTime)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert record to dictionary for JSON serialization."""
        return {
            'id': self.id,
            'dataset_id': self.dataset_id,
            'filename': self.filename,
            'stressor': self.stressor,
            'dosing': self.dosing,
            'owner': self.owner,
            'description': self.description,
            'aop_id': self.aop_id,
            'logfc_threshold': self.logfc_threshold,
            'pval_cutoff': self.pval_cutoff,
            'id_column': self.id_column,
            'fc_column': self.fc_column,
            'pval_column': self.pval_column,
            'enrichment_results': json.loads(self.enrichment_results) if self.enrichment_results else None,
            'gene_count': self.gene_count,
            'significant_genes': self.significant_genes,
            'upload_timestamp': self.upload_timestamp.isoformat() if self.upload_timestamp else None,
            'analysis_timestamp': self.analysis_timestamp.isoformat() if self.analysis_timestamp else None,
        }


class DatabaseManager:
    """Manager class for database operations."""
    
    def __init__(self, db_url: str = "sqlite:///molAOP_analyser.db"):
        """Initialize database manager.
        
        Args:
            db_url: SQLAlchemy database URL (defaults to SQLite)
        """
        self.db_url = db_url
        self.engine = None
        self.SessionLocal = None
        
    def initialize(self):
        """Initialize database connection and create tables."""
        try:
            self.engine = create_engine(
                self.db_url,
                echo=False,  # Set to True for SQL debugging
                pool_pre_ping=True,  # Verify connections before use
            )
            
            # Create all tables
            Base.metadata.create_all(bind=self.engine)
            
            # Create session factory
            self.SessionLocal = sessionmaker(
                autocommit=False,
                autoflush=False,
                bind=self.engine
            )
            
            logger.info(f"Database initialized successfully: {self.db_url}")
            return True
            
        except SQLAlchemyError as e:
            logger.error(f"Failed to initialize database: {e}")
            return False
    
    def get_session(self):
        """Get a new database session."""
        if not self.SessionLocal:
            raise RuntimeError("Database not initialized. Call initialize() first.")
        return self.SessionLocal()
    
    def save_experiment_metadata(self, metadata: Dict[str, Any], 
                                analysis_params: Optional[Dict[str, Any]] = None,
                                results: Optional[Dict[str, Any]] = None) -> Optional[int]:
        """Save experiment metadata to database.
        
        Args:
            metadata: Experiment metadata dictionary
            analysis_params: Optional analysis parameters
            results: Optional analysis results
            
        Returns:
            Experiment record ID if successful, None otherwise
        """
        session = self.get_session()
        try:
            record = ExperimentRecord(
                dataset_id=metadata.get('dataset_id', ''),
                filename=metadata.get('filename', ''),
                stressor=metadata.get('stressor', ''),
                dosing=metadata.get('dosing', ''),
                owner=metadata.get('owner', ''),
                description=metadata.get('description', ''),
                upload_timestamp=datetime.fromisoformat(metadata['upload_timestamp']) 
                    if metadata.get('upload_timestamp') else datetime.utcnow()
            )
            
            # Add analysis parameters if provided
            if analysis_params:
                record.aop_id = analysis_params.get('aop_id')
                record.logfc_threshold = analysis_params.get('logfc_threshold')
                record.pval_cutoff = analysis_params.get('pval_cutoff')
                record.id_column = analysis_params.get('id_column')
                record.fc_column = analysis_params.get('fc_column')
                record.pval_column = analysis_params.get('pval_column')
                record.analysis_timestamp = datetime.utcnow()
            
            # Add results summary if provided
            if results:
                record.enrichment_results = json.dumps(results.get('enrichment_table', []))
                record.gene_count = results.get('gene_count')
                record.significant_genes = results.get('significant_genes')
            
            session.add(record)
            session.commit()
            
            experiment_id = record.id
            logger.info(f"Saved experiment record with ID: {experiment_id}")
            return experiment_id
            
        except SQLAlchemyError as e:
            session.rollback()
            logger.error(f"Failed to save experiment metadata: {e}")
            return None
        finally:
            session.close()
    
    def get_experiment(self, experiment_id: int) -> Optional[Dict[str, Any]]:
        """Retrieve experiment record by ID.
        
        Args:
            experiment_id: Experiment record ID
            
        Returns:
            Experiment record dictionary if found, None otherwise
        """
        session = self.get_session()
        try:
            record = session.query(ExperimentRecord).filter_by(id=experiment_id).first()
            return record.to_dict() if record else None
        except SQLAlchemyError as e:
            logger.error(f"Failed to retrieve experiment {experiment_id}: {e}")
            return None
        finally:
            session.close()
    
    def list_experiments(self, limit: int = 50, offset: int = 0) -> list:
        """List recent experiments.
        
        Args:
            limit: Maximum number of records to return
            offset: Number of records to skip
            
        Returns:
            List of experiment record dictionaries
        """
        session = self.get_session()
        try:
            records = session.query(ExperimentRecord)\
                           .order_by(ExperimentRecord.upload_timestamp.desc())\
                           .limit(limit).offset(offset).all()
            
            return [record.to_dict() for record in records]
            
        except SQLAlchemyError as e:
            logger.error(f"Failed to list experiments: {e}")
            return []
        finally:
            session.close()
    
    def search_experiments(self, dataset_id: Optional[str] = None,
                          owner: Optional[str] = None,
                          stressor: Optional[str] = None) -> list:
        """Search experiments by metadata fields.
        
        Args:
            dataset_id: Filter by dataset ID (partial match)
            owner: Filter by owner (partial match)  
            stressor: Filter by stressor (partial match)
            
        Returns:
            List of matching experiment record dictionaries
        """
        session = self.get_session()
        try:
            query = session.query(ExperimentRecord)
            
            if dataset_id:
                query = query.filter(ExperimentRecord.dataset_id.contains(dataset_id))
            if owner:
                query = query.filter(ExperimentRecord.owner.contains(owner))
            if stressor:
                query = query.filter(ExperimentRecord.stressor.contains(stressor))
            
            records = query.order_by(ExperimentRecord.upload_timestamp.desc()).all()
            return [record.to_dict() for record in records]
            
        except SQLAlchemyError as e:
            logger.error(f"Failed to search experiments: {e}")
            return []
        finally:
            session.close()


# Global database manager instance
db_manager = DatabaseManager()


def init_database() -> bool:
    """Initialize the database. Should be called during app startup.
    
    Returns:
        True if initialization successful, False otherwise
    """
    return db_manager.initialize()


def get_db_session():
    """Get a database session for manual operations."""
    return db_manager.get_session()