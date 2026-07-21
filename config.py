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

    # Enrichment statistics — issues #63 / #65
    #
    # SIGNIFICANCE_FDR_CUTOFF is the single definition of "a Key Event is
    # significant" used across the enrichment table, the red network borders,
    # the cross-condition comparison matrix and the reports. It gates the
    # Benjamini-Hochberg adjusted p-value (the `FDR` column) — never the raw
    # Fisher p-value, which over-calls enrichment across the 10-20 KEs of a
    # typical AOP.
    SIGNIFICANCE_FDR_CUTOFF = 0.05

    # The cutoffs offered by the results-page significance control. A discrete
    # set rather than a continuous slider: the conventional values are what
    # anyone actually reports, and a slider's step made everything between 0
    # and its first stop unreachable while leaving 0 (nothing significant) one
    # tick from the default. SIGNIFICANCE_FDR_CUTOFF must appear in this list.
    SIGNIFICANCE_FDR_CHOICES = (0.001, 0.01, 0.05, 0.10, 0.25)

    # MIN_KE_GENES is the minimum number of a KE's gene-set members that must
    # be measured in the uploaded dataset before that KE is tested. KEs below
    # it are excluded from the enrichment — and therefore from the BH
    # multiple-testing denominator — and are reported as "could not assess"
    # rather than "assessed and not enriched". Fisher's exact test on a handful
    # of genes has too little power for the result to be interpretable.
    MIN_KE_GENES = 5


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
    
    # SPARQL endpoint settings
    SPARQL_ENDPOINT = os.environ.get(
        'SPARQL_ENDPOINT',
        'https://aopwiki.rdf.bigcat-bioinformatics.org/sparql/'
    )
    SPARQL_TIMEOUT = int(os.environ.get('SPARQL_TIMEOUT', '30'))

    # Kidney AOP IDs for the combined kidney network.
    #
    # Curation criterion: this is a manually curated curator's working list —
    # there is no single reproducible rule. The set is a mix of AOPs whose
    # Adverse Outcome references kidney / renal failure and AOPs flagged as
    # relevant by VHP4Safety collaborators in the kidney case study. To
    # refresh, the maintainer reviews AOP-Wiki AOs touching renal endpoints
    # and re-confirms the list with the kidney case-study lead; there is no
    # script that reproduces this list automatically.
    #
    # Last reviewed: 2026-05-06
    # Maintainer: Marvin Martens (m.martens@maastrichtuniversity.nl)
    KIDNEY_AOP_IDS = [
        "AOP:33", "AOP:53", "AOP:105", "AOP:116", "AOP:128",
        "AOP:138", "AOP:177", "AOP:186", "AOP:256", "AOP:257",
        "AOP:258", "AOP:284", "AOP:384", "AOP:413", "AOP:437",
        "AOP:447", "AOP:472", "AOP:622",
    ]

    # Recommended AOPs per demo dataset for the /demos curated path (Phase 10.1, D-08/D-09).
    # Keys are demo file paths (relative to data/) or glob patterns matching DEMO_DATASETS keys.
    # Values are ordered AOP IDs (or NETWORK:* sentinels) presented as recommended on /preview
    # when the user lands via a /demos CTA. Soft restriction (D-06) — full list still reachable
    # via the "Show all AOPs" toggle owned by Plan 02.
    DEMO_AOP_RECOMMENDATIONS = {
        # PXR demos → locally curated AOP:DEMO (PXR activation → liver steatosis).
        # Namespaced away from AOP:N to avoid colliding with the AOP-Wiki id
        # space — AOP-Wiki's AOP:1 is a different pathway (uncharacterized
        # liver damage → HCC) and was overwriting our label in the typeahead.
        'GSE90122_TO90137.tsv': ['AOP:DEMO'],
        'GSE90122_SR12813.tsv': ['AOP:DEMO'],
        # Cisplatin demos → kidney network first + 6 individual kidney AOPs (D-08).
        # Order: NETWORK:kidney first (locked); remaining order is Claude's discretion
        # within the curator's most-cited kidney AOP list.
        'Cisplatin_Kidney/*': [
            'NETWORK:kidney',
            'AOP:472',
            'AOP:177',
            'AOP:447',
            'AOP:622',
            'AOP:258',
            'AOP:437',
        ],
    }

    @classmethod
    def get_recommended_aops(cls, demo_file: str) -> list:
        """Return recommended AOP IDs for a demo file path; [] if no match.

        Matches in two passes: exact key first, then fnmatch glob fallback.
        """
        import fnmatch
        if not demo_file:
            return []
        if demo_file in cls.DEMO_AOP_RECOMMENDATIONS:
            return list(cls.DEMO_AOP_RECOMMENDATIONS[demo_file])
        for pattern, aops in cls.DEMO_AOP_RECOMMENDATIONS.items():
            if '*' in pattern and fnmatch.fnmatch(demo_file, pattern):
                return list(aops)
        return []

    # AOP case studies
    CASE_STUDY_AOPS = {
        "DEMO": {"label": "---DEMO---", "enabled": False},
        "steatosis": {"id": "AOP:DEMO", "label": "PXR activation leading to liver steatosis", "enabled": True, "source": "csv"},
        "VHP-CASES:": {"label": "---VHP CASES---", "enabled": False},
        "vhp-kidney": {"id": "AOP:2", "label": "DNA adduct formation leading to kidney failure", "enabled": True, "source": "csv"},
        "vhp-parkinson": {"id": "AOP:3", "label": "Calcium overload in dopaminergic neurons of the substantia nigra leading to parkinsonian motor deficits", "enabled": False, "source": "csv"},
        "vhp-thyroid": {"id": "AOP:4", "label": "Thyroid hormone-mediated neurodevelopmental toxicity", "enabled": False, "source": "csv"},
        "KIDNEY-AOPS": {"label": "---KIDNEY AOPs---", "enabled": False},
        "kidney-33": {"id": "AOP:33", "label": "Kidney toxicity induced by activation of 5HT2C", "enabled": True, "source": "sparql"},
        "kidney-53": {"id": "AOP:53", "label": "ER agonism leading to reduced survival due to renal failure", "enabled": True, "source": "sparql"},
        "kidney-105": {"id": "AOP:105", "label": "Alpha2u-microglobulin cytotoxicity leading to renal tubular adenomas and carcinomas", "enabled": True, "source": "sparql"},
        "kidney-116": {"id": "AOP:116", "label": "Cytotoxicity leading to renal tubular adenomas and carcinomas", "enabled": True, "source": "sparql"},
        "kidney-128": {"id": "AOP:128", "label": "Kidney dysfunction by decreased thyroid hormone", "enabled": True, "source": "sparql"},
        "kidney-138": {"id": "AOP:138", "label": "OAT1 inhibition leading to renal failure and mortality", "enabled": True, "source": "sparql"},
        "kidney-177": {"id": "AOP:177", "label": "COX1 inhibition leading to renal failure and mortality", "enabled": True, "source": "sparql"},
        "kidney-186": {"id": "AOP:186", "label": "Unknown MIE leading to renal failure and mortality", "enabled": True, "source": "sparql"},
        "kidney-256": {"id": "AOP:256", "label": "Inhibition of mtDNA polymerase gamma leading to kidney toxicity", "enabled": True, "source": "sparql"},
        "kidney-257": {"id": "AOP:257", "label": "Receptor mediated endocytosis and lysosomal overload leading to kidney toxicity", "enabled": True, "source": "sparql"},
        "kidney-258": {"id": "AOP:258", "label": "Renal protein alkylation leading to kidney toxicity", "enabled": True, "source": "sparql"},
        "kidney-284": {"id": "AOP:284", "label": "Binding to SH-group proteins leading to chronic kidney disease", "enabled": True, "source": "sparql"},
        "kidney-384": {"id": "AOP:384", "label": "Hyperactivation of ACE/Ang-II/AT1R axis leading to chronic kidney disease", "enabled": True, "source": "sparql"},
        "kidney-413": {"id": "AOP:413", "label": "Oxidation of reduced glutathione leading to mortality via acute renal failure", "enabled": True, "source": "sparql"},
        "kidney-437": {"id": "AOP:437", "label": "Inhibition of mitochondrial ETC complexes leading to kidney toxicity", "enabled": True, "source": "sparql"},
        "kidney-447": {"id": "AOP:447", "label": "Inhibition of mitochondrial ETC leading to kidney failure", "enabled": True, "source": "sparql"},
        "kidney-622": {"id": "AOP:622", "label": "Calcineurin inhibitor induced nephrotoxicity leading to kidney failure", "enabled": True, "source": "sparql"},
        "ORGAN-NETWORK": {"label": "---ORGAN NETWORKS---", "enabled": False},
        "Kidney-aop-network": {"id": "NETWORK:kidney", "label": "Kidney AOP network (all 18 AOPs combined)", "enabled": True, "source": "sparql"},
        "Liver-aop-network": {"id": "AOP:5", "label": "Liver AOP network", "enabled": False},
        "Brain-aop-network": {"id": "AOP:6", "label": "Brain AOP network", "enabled": False},
        "Lung-aop-network": {"id": "AOP:8", "label": "Lung AOP network", "enabled": False},
    }

    # Wire up aop_ids for the kidney network entry
    CASE_STUDY_AOPS["Kidney-aop-network"]["aop_ids"] = KIDNEY_AOP_IDS  # noqa: E501

    # Flask settings.
    #
    # SECRET_KEY_FILE takes precedence over SECRET_KEY: the deployed service
    # points it at a Docker swarm secret (/run/secrets/...), so the key is
    # never exposed in the service's environment where `docker service
    # inspect` would print it. SECRET_KEY stays supported for local dev.
    #
    # Without either, a random key is generated per process — sessions then
    # do not survive a restart, and with more than one replica they would not
    # be shared between them.
    _secret = None
    _secret_file = os.environ.get('SECRET_KEY_FILE')
    if _secret_file:
        try:
            with open(_secret_file, 'r', encoding='utf-8') as fh:
                _secret = fh.read().strip()
        except OSError as exc:
            import warnings
            warnings.warn(f"SECRET_KEY_FILE {_secret_file!r} could not be read ({exc})")
        else:
            if not _secret:
                import warnings
                warnings.warn(f"SECRET_KEY_FILE {_secret_file!r} is empty")

    if not _secret:
        _secret = os.environ.get('SECRET_KEY')

    if not _secret:
        import warnings
        import secrets as _secrets
        warnings.warn("SECRET_KEY not set — using random key (sessions won't persist across restarts)")
        _secret = _secrets.token_hex(32)
    SECRET_KEY = _secret

    # Session cookie security
    SESSION_COOKIE_HTTPONLY = True
    SESSION_COOKIE_SAMESITE = 'Lax'
    PERMANENT_SESSION_LIFETIME = 86400  # 24 hours
    SESSION_COOKIE_SECURE = os.environ.get('SESSION_COOKIE_SECURE', '0') == '1'

    # Builder API settings
    BUILDER_API_URL = os.environ.get('BUILDER_API_URL', 'https://molaop-builder.vhp4safety.nl/')
    BUILDER_API_TIMEOUT = int(os.environ.get('BUILDER_API_TIMEOUT', '30'))
    # Reference-set disk cache. Issue #68: on /tmp the cache is container-local
    # and dies with every restart, so identical submissions minutes apart could
    # resolve to different gene-set sources. When the deployed service's mounted
    # volume is present the cache lives there instead, making cache behaviour
    # stable across restarts. CACHE_DIR overrides either way.
    CACHE_DIR = os.environ.get('CACHE_DIR') or (
        '/data/molaop_cache'
        if os.path.isdir('/data') and os.access('/data', os.W_OK)
        else '/tmp/molaop_cache'
    )
    CACHE_TTL = 3600  # 1 hour

    # Database location. Defaults to a file in the working directory, which in
    # a container is ephemeral — the deployed service overrides this to point
    # at a mounted volume so shared-result links and batch history survive a
    # redeploy. Keep the swarm service on `stop-first` when it does: SQLite on
    # a shared mount tolerates exactly one writer.
    DATABASE_URL = os.environ.get('DATABASE_URL', 'sqlite:///molAOP_analyser.db')

    # Batch analysis settings
    BATCH_MAX_FILES = 10
    BATCH_MAX_TOTAL_SIZE = 100 * 1024 * 1024  # 100 MB total batch limit
    BATCH_RETENTION_DAYS = 14  # Auto-clean batches older than 14 days
    BATCH_MIN_HARMONISED_GENES = 1000  # Minimum viable background for Fisher's test

    # Cisplatin batch demo files grouped by timepoint
    CISPLATIN_TIMEPOINTS = ['4hr', '8hr', '16hr', '24hr', '48hr', '72hr']
    CISPLATIN_DOSES = ['0.1uM', '0.5uM', '1uM', '2.5uM', '5uM', '10uM', '20uM', '30uM', '50uM']

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