"""
Database models and setup for the Molecular AOP Analyser.

Provides SQLite-based persistence for experiment metadata and analysis results.
"""

from datetime import datetime, timedelta, timezone
from typing import Optional, Dict, Any
import json
import uuid as uuid_lib
from sqlalchemy import create_engine, Column, Integer, String, DateTime, Text, Float, ForeignKey, text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.exc import SQLAlchemyError
import logging
from config import Config

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
    method = Column(String(20))  # 'ora' or 'gsea'; NULL for pre-Phase-14 rows (coerced on read)
    logfc_threshold = Column(Float)
    pval_cutoff = Column(Float)
    id_column = Column(String(100))
    fc_column = Column(String(100))
    pval_column = Column(String(100))
    # Issue #55: comma-separated gene-set resources used (e.g. "WikiPathways, GO_BP").
    # NULL for pre-#55 rows (coerced to WikiPathways on read).
    selected_resources = Column(String(255))
    # Issue #60: minimum KE-mapping confidence used for enrichment
    # ('all'/'medium'/'high'). NULL for pre-#60 rows (coerced to 'all' on read).
    min_confidence = Column(String(20))
    # Issue #68: how each requested resource actually RESOLVED (JSON list of
    # per-resource dicts: source, status, ke_count, confidence_applied). The
    # column above records the request; this one records the outcome, which is
    # what makes a run reproducible. NULL for pre-#68 rows.
    resource_resolution = Column(Text)

    # Analysis results summary (JSON)
    enrichment_results = Column(Text)  # JSON string of enrichment results
    gene_count = Column(Integer)
    significant_genes = Column(Integer)
    
    # Timestamps  
    upload_timestamp = Column(DateTime, default=lambda: datetime.now(timezone.utc).replace(tzinfo=None), nullable=True)
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
            'method': self.method or 'ora',  # D-05: coerce NULL (pre-Phase-14 rows) to 'ora'
            'logfc_threshold': self.logfc_threshold,
            'pval_cutoff': self.pval_cutoff,
            'id_column': self.id_column,
            'fc_column': self.fc_column,
            'pval_column': self.pval_column,
            'selected_resources': self.selected_resources or 'WikiPathways',  # coerce NULL (pre-#55 rows)
            'min_confidence': self.min_confidence or 'all',  # coerce NULL (pre-#60 rows)
            'enrichment_results': json.loads(self.enrichment_results) if self.enrichment_results else None,
            'gene_count': self.gene_count,
            'significant_genes': self.significant_genes,
            'upload_timestamp': self.upload_timestamp.isoformat() if self.upload_timestamp else None,
            'analysis_timestamp': self.analysis_timestamp.isoformat() if self.analysis_timestamp else None,
        }


class SharedResult(Base):
    """Stores serialized analysis data for shareable UUID links.

    Created on-demand when user clicks 'Copy link' on results page.
    Expires after 30 days. Independent from ExperimentRecord.
    """

    __tablename__ = 'shared_results'

    uuid = Column(String(36), primary_key=True, default=lambda: str(uuid_lib.uuid4()))
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc).replace(tzinfo=None), nullable=False)
    expires_at = Column(DateTime, nullable=False)

    # Analysis identity
    aop_id = Column(String(100))
    aop_label = Column(String(500))

    # Serialized results (JSON Text columns)
    enrichment_json = Column(Text, nullable=False)
    network_json = Column(Text, nullable=False)
    ke_gene_json = Column(Text)
    ke_type_map_json = Column(Text)
    ke_title_map_json = Column(Text)

    # Display metadata
    dataset_label = Column(String(500))
    gene_count = Column(Integer)
    significant_genes = Column(Integer)
    method = Column(String(20))  # 'ora' or 'gsea'; NULL for pre-Phase-14 rows (coerced on read)

    @classmethod
    def create(cls, enrichment_table, network, ke_gene_map, ke_type_map,
               ke_title_map, aop_id, aop_label, metadata):
        """Factory method to create a SharedResult with 30-day expiry.

        Args:
            enrichment_table: List of enrichment result dicts
            network: Cytoscape network dict (nodes/edges)
            ke_gene_map: Dict mapping KE IDs to gene lists
            ke_type_map: Dict mapping KE IDs to type strings
            ke_title_map: Dict mapping KE IDs to title strings
            aop_id: AOP identifier string
            aop_label: Human-readable AOP label
            metadata: Dict with dataset_id, gene_count, significant_genes

        Returns:
            SharedResult instance with UUID and 30-day expiry
        """
        expiry = datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(days=30)
        return cls(
            uuid=str(uuid_lib.uuid4()),
            expires_at=expiry,
            aop_id=aop_id,
            aop_label=aop_label,
            enrichment_json=json.dumps(enrichment_table),
            network_json=json.dumps(network),
            ke_gene_json=json.dumps(ke_gene_map) if ke_gene_map else None,
            ke_type_map_json=json.dumps(ke_type_map) if ke_type_map else None,
            ke_title_map_json=json.dumps(ke_title_map) if ke_title_map else None,
            dataset_label=metadata.get('dataset_id', ''),
            gene_count=metadata.get('gene_count'),
            significant_genes=metadata.get('significant_genes'),
            method=metadata.get('method', 'ora'),  # D-06: persist method; legacy callers get 'ora' default
        )


def cleanup_expired_shared_results(session):
    """Delete shared results past their expiry date.

    Called lazily on each new share creation to avoid accumulating stale rows.

    Args:
        session: SQLAlchemy database session
    """
    session.query(SharedResult).filter(
        SharedResult.expires_at < datetime.now(timezone.utc).replace(tzinfo=None)
    ).delete(synchronize_session=False)
    session.commit()


class BatchRecord(Base):
    """Stores batch analysis state for multi-file enrichment runs.

    A batch groups multiple ConditionRecords (one per uploaded file) under a
    shared AOP selection, thresholds, and harmonised background gene set.
    Expires after BATCH_RETENTION_DAYS and is cleaned up lazily on new batch creation.
    """

    __tablename__ = 'batches'

    # Primary key
    id = Column(Integer, primary_key=True, autoincrement=True)

    # Identity
    uuid = Column(String(36), unique=True, nullable=False, default=lambda: str(uuid_lib.uuid4()))
    status = Column(String(20), default='pending')  # pending, running, complete, cancelled, failed

    # AOP selection
    aop_id = Column(String(100))
    aop_label = Column(String(500))

    # Shared analysis parameters
    # Issue #76: enrichment method applied to every condition in the batch
    # ('ora' or 'gsea'). NULL for rows written before #76 — coerced to 'ora'
    # on read so pre-existing batches keep opening, exporting and comparing
    # exactly as they did.
    #
    # The column deliberately lives here and not on ``batch_conditions``: the
    # method is a batch-level parameter, like logfc_threshold, pval_cutoff,
    # selected_resources and min_confidence. Cross-condition comparison is only
    # meaningful when every condition was scored the same way, so a per-
    # condition column could only ever encode an invalid state.
    method = Column(String(20))
    logfc_threshold = Column(Float)
    pval_cutoff = Column(Float)
    # Issue #55: comma-separated gene-set resources used (e.g. "WikiPathways, GO_BP").
    selected_resources = Column(String(255))
    # Issue #60: minimum KE-mapping confidence applied to every condition
    # ('all'/'medium'/'high'). NULL for pre-#60 rows (read as 'all').
    min_confidence = Column(String(20))
    # Issue #68: per-resource resolution shared by every condition in the batch
    # (JSON, see ExperimentRecord.resource_resolution). NULL for pre-#68 rows.
    resource_resolution = Column(Text)

    # Shared column mapping (same for all files in the batch)
    id_column = Column(String(100))
    fc_column = Column(String(100))
    pval_column = Column(String(100))

    # Harmonised background gene set (JSON array of sorted gene symbols)
    harmonised_background = Column(Text)
    harmonised_gene_count = Column(Integer)

    # Batch metadata
    batch_name = Column(String(500))
    owner = Column(String(255))
    description = Column(Text)

    # Timestamps
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc).replace(tzinfo=None), nullable=False)
    completed_at = Column(DateTime, nullable=True)
    expires_at = Column(DateTime, nullable=False)

    # Relationship to per-file condition records (ordered by upload position)
    conditions = relationship(
        'ConditionRecord',
        back_populates='batch',
        order_by='ConditionRecord.position',
    )

    def effective_method(self) -> str:
        """Return the enrichment method this batch was run with.

        Coerces NULL to ``'ora'`` (issue #76), the behaviour of every batch
        created before the method column existed. Callers should use this
        rather than reading ``.method`` directly so a legacy batch never
        reaches a GSEA-only code path.

        Returns:
            ``'ora'`` or ``'gsea'``.
        """
        return self.method or 'ora'


class ConditionRecord(Base):
    """Stores per-file analysis state within a batch run.

    Each uploaded file becomes one ConditionRecord attached to a BatchRecord.
    Tracks file metadata, execution status, enrichment results, and an optional
    deep-link UUID for the single-condition /results/<uuid> view.
    """

    __tablename__ = 'batch_conditions'

    # Primary key
    id = Column(Integer, primary_key=True, autoincrement=True)

    # Foreign key to parent batch
    batch_id = Column(Integer, ForeignKey('batches.id'), nullable=False, index=True)

    # Position within the batch (0-based upload order)
    position = Column(Integer, nullable=False)

    # File identity
    filename = Column(String(255))
    condition_label = Column(String(500), nullable=False)

    # Parsed condition metadata
    dose = Column(String(200))
    timepoint = Column(String(200))

    # Execution state
    status = Column(String(20), default='pending')  # pending, running, complete, failed
    error_message = Column(Text, nullable=True)

    # Gene counts from analysis
    gene_count = Column(Integer)
    significant_genes = Column(Integer)

    # Issue #69: fraction of this condition's identifiers that matched the
    # reference gene universe. NULL on rows written before the check existed.
    # A low value means the gene ID column was probably not gene symbols and
    # the condition's results, though well-formed, are empty of real overlap.
    id_match_fraction = Column(Float, nullable=True)

    # Serialised analysis outputs (JSON Text columns)
    enrichment_json = Column(Text)      # Full enrichment table as JSON
    network_json = Column(Text)         # Cytoscape network JSON
    ke_gene_json = Column(Text)         # KE-to-gene membership map JSON
    ke_type_map_json = Column(Text)     # KE ID to type string map JSON
    ke_title_map_json = Column(Text)    # KE ID to title string map JSON

    # Deep-link UUID for the shared /results/<uuid> page
    shared_result_uuid = Column(String(36), nullable=True)

    # Execution timestamps
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)

    # Relationship back to parent batch
    batch = relationship('BatchRecord', back_populates='conditions')


def cleanup_expired_batches(session):
    """Delete batches past their expires_at date.

    Called lazily on new batch creation to avoid accumulating stale rows.
    Child ConditionRecords are deleted explicitly before the parent BatchRecord
    to avoid relying on cascade behaviour.

    Args:
        session: SQLAlchemy database session
    """
    expired = session.query(BatchRecord).filter(
        BatchRecord.expires_at < datetime.now(timezone.utc).replace(tzinfo=None)
    ).all()
    for batch in expired:
        session.query(ConditionRecord).filter_by(
            batch_id=batch.id
        ).delete(synchronize_session=False)
        session.delete(batch)
    session.commit()


def _ensure_method_column(engine) -> None:
    """Idempotent PRAGMA-then-ALTER migration for the 'method' column.

    Adds ``method TEXT DEFAULT 'ora'`` to the ``experiments``,
    ``shared_results`` and ``batches`` tables when the column is absent.
    Safe to run on every startup — the PRAGMA check makes the ALTER a no-op
    on databases that already have the column (D-04, D-06, #76).

    ``batches`` joined the list with issue #76, which gave batch runs a
    method selector. SQLite's ``ADD COLUMN … DEFAULT 'ora'`` backfills every
    existing row with ``'ora'`` rather than leaving it NULL, so a batch created
    before the change reads back as ORA directly — which is what it was. NULL is
    still handled (a row written by another path, or an older migration that
    added the column without a default) and
    :meth:`BatchRecord.effective_method` coerces it to ``'ora'`` at every read
    site, so either way such a batch opens, runs, exports and compares as ORA.

    Security note: both table names and the column literal are module-internal
    constants — no user-supplied input is interpolated into the SQL statements.

    Args:
        engine: SQLAlchemy engine bound to the target database.
    """
    for table in ('experiments', 'shared_results', 'batches'):
        with engine.connect() as conn:
            result = conn.execute(text(f"PRAGMA table_info({table})"))
            existing_cols = {row[1] for row in result}
            if 'method' not in existing_cols:
                conn.execute(
                    text(f"ALTER TABLE {table} ADD COLUMN method TEXT DEFAULT 'ora'")
                )
                conn.commit()
                logger.info(f"Added 'method' column to '{table}' table (D-04 migration)")


def _ensure_selected_resources_column(engine) -> None:
    """Idempotent PRAGMA-then-ALTER migration for the 'selected_resources' column.

    Adds ``selected_resources TEXT`` to the ``experiments`` and ``batches``
    tables when absent (issue #55). Safe to run on every startup — the PRAGMA
    check makes the ALTER a no-op on databases that already have the column.

    Security note: the table names and column literal are module-internal
    constants — no user-supplied input is interpolated into the SQL statements.

    Args:
        engine: SQLAlchemy engine bound to the target database.
    """
    for table in ('experiments', 'batches'):
        with engine.connect() as conn:
            result = conn.execute(text(f"PRAGMA table_info({table})"))
            existing_cols = {row[1] for row in result}
            if 'selected_resources' not in existing_cols:
                conn.execute(
                    text(f"ALTER TABLE {table} ADD COLUMN selected_resources TEXT")
                )
                conn.commit()
                logger.info(
                    f"Added 'selected_resources' column to '{table}' table (#55 migration)"
                )


def _ensure_min_confidence_column(engine) -> None:
    """Idempotent PRAGMA-then-ALTER migration for the 'min_confidence' column.

    Adds ``min_confidence TEXT`` to the ``experiments`` and ``batches`` tables
    when absent (issue #60). Safe to run on every startup — the PRAGMA check
    makes the ALTER a no-op on databases that already have the column. Existing
    rows keep NULL, which is read back as ``'all'`` (the pre-#60 behaviour).

    Security note: the table names and column literal are module-internal
    constants — no user-supplied input is interpolated into the SQL statements.

    Args:
        engine: SQLAlchemy engine bound to the target database.
    """
    for table in ('experiments', 'batches'):
        with engine.connect() as conn:
            result = conn.execute(text(f"PRAGMA table_info({table})"))
            existing_cols = {row[1] for row in result}
            if 'min_confidence' not in existing_cols:
                conn.execute(
                    text(f"ALTER TABLE {table} ADD COLUMN min_confidence TEXT")
                )
                conn.commit()
                logger.info(
                    f"Added 'min_confidence' column to '{table}' table (#60 migration)"
                )


def _ensure_id_match_fraction_column(engine) -> None:
    """Idempotent PRAGMA-then-ALTER migration for 'id_match_fraction' (#69).

    Adds ``id_match_fraction REAL`` to ``batch_conditions`` when absent. Safe on
    every startup; existing rows keep NULL, which the UI reads as "not checked"
    rather than as a failed check.

    Security note: table and column names are module-internal constants — no
    user input is interpolated into the SQL.

    Args:
        engine: SQLAlchemy engine bound to the target database.
    """
    with engine.connect() as conn:
        result = conn.execute(text("PRAGMA table_info(batch_conditions)"))
        existing_cols = {row[1] for row in result}
        if 'id_match_fraction' not in existing_cols:
            conn.execute(
                text("ALTER TABLE batch_conditions ADD COLUMN id_match_fraction REAL")
            )
            conn.commit()
            logger.info(
                "Added 'id_match_fraction' column to 'batch_conditions' table (#69 migration)"
            )


def _ensure_resource_resolution_column(engine) -> None:
    """Idempotent PRAGMA-then-ALTER migration for 'resource_resolution' (#68).

    Adds ``resource_resolution TEXT`` to ``experiments`` and ``batches`` when
    absent. Existing rows keep NULL, which the UI reads as "not recorded" — the
    pre-#68 behaviour, where only the requested resource list was stored.

    Security note: table and column names are module-internal constants — no
    user input is interpolated into the SQL.

    Args:
        engine: SQLAlchemy engine bound to the target database.
    """
    for table in ('experiments', 'batches'):
        with engine.connect() as conn:
            result = conn.execute(text(f"PRAGMA table_info({table})"))
            existing_cols = {row[1] for row in result}
            if 'resource_resolution' not in existing_cols:
                conn.execute(
                    text(f"ALTER TABLE {table} ADD COLUMN resource_resolution TEXT")
                )
                conn.commit()
                logger.info(
                    f"Added 'resource_resolution' column to '{table}' table (#68 migration)"
                )


class DatabaseManager:
    """Manager class for database operations."""
    
    def __init__(self, db_url: Optional[str] = None):
        """Initialize database manager.

        Args:
            db_url: SQLAlchemy database URL. Defaults to Config.DATABASE_URL,
                which the deployed service points at a mounted volume so the
                database is not lost on redeploy.
        """
        self.db_url = db_url or Config.DATABASE_URL
        self.engine = None
        self.SessionLocal = None
        
    def initialize(self):
        """Initialize database connection and create tables."""
        try:
            engine_kwargs = dict(echo=False)
            if 'sqlite' in self.db_url:
                # SQLite does not support concurrent writers; use StaticPool
                # with check_same_thread=False for multi-threaded Flask/batch use.
                from sqlalchemy.pool import StaticPool
                engine_kwargs.update(
                    poolclass=StaticPool,
                    connect_args={"check_same_thread": False},
                )
            else:
                engine_kwargs.update(
                    pool_pre_ping=True,
                    pool_size=10,
                    max_overflow=20,
                )
            self.engine = create_engine(self.db_url, **engine_kwargs)
            
            # Create all tables
            Base.metadata.create_all(bind=self.engine)

            # Idempotent additive migration: add 'method' column to pre-existing
            # databases that lack it (D-04, D-06). No-op on fresh databases where
            # the column is already present from the model definition above.
            _ensure_method_column(self.engine)

            # Idempotent additive migration: add 'selected_resources' column to
            # pre-existing databases that lack it (#55). No-op on fresh databases.
            _ensure_selected_resources_column(self.engine)

            # Idempotent additive migration: add 'min_confidence' column to
            # pre-existing databases that lack it (#60). No-op on fresh databases.
            _ensure_min_confidence_column(self.engine)

            # Idempotent additive migration: add 'id_match_fraction' to
            # batch_conditions (#69). No-op on fresh databases.
            _ensure_id_match_fraction_column(self.engine)

            # Idempotent additive migration: add 'resource_resolution' to
            # experiments and batches (#68). No-op on fresh databases.
            _ensure_resource_resolution_column(self.engine)

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
                    if metadata.get('upload_timestamp') else datetime.now(timezone.utc).replace(tzinfo=None)
            )
            
            # Add analysis parameters if provided
            if analysis_params:
                record.aop_id = analysis_params.get('aop_id')
                record.method = analysis_params.get('method', 'ora')  # D-04: persist chosen method
                record.logfc_threshold = analysis_params.get('logfc_threshold')
                record.pval_cutoff = analysis_params.get('pval_cutoff')
                record.id_column = analysis_params.get('id_column')
                record.fc_column = analysis_params.get('fc_column')
                record.pval_column = analysis_params.get('pval_column')
                record.selected_resources = analysis_params.get('selected_resources')
                record.min_confidence = analysis_params.get('min_confidence')  # #60
                record.resource_resolution = analysis_params.get('resource_resolution')  # #68
                record.analysis_timestamp = datetime.now(timezone.utc).replace(tzinfo=None)
            
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