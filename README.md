# Molecular AOP Analyser

[![CI/CD Pipeline](https://github.com/marvinm2/molAOP-analyser/actions/workflows/ci.yml/badge.svg)](https://github.com/marvinm2/molAOP-analyser/actions/workflows/ci.yml)
[![Code Quality](https://github.com/marvinm2/molAOP-analyser/actions/workflows/code-quality.yml/badge.svg)](https://github.com/marvinm2/molAOP-analyser/actions/workflows/code-quality.yml)
[![Publish Docker image](https://github.com/marvinm2/molAOP-analyser/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/marvinm2/molAOP-analyser/actions/workflows/docker-publish.yml)

This web application allows users to upload or select gene expression datasets and perform Key Event (KE) enrichment analysis in the context of Molecular Adverse Outcome Pathways (AOPs). The results are visualized in interactive tables and network diagrams with comprehensive reporting capabilities.

**Live instance:** https://molaop-analyser.vhp4safety.nl
**Container image:** `ghcr.io/marvinm2/molaop-analyser`

## Related

The analyser is the downstream half of a pair of tools.

| Tool | Role | Links |
|---|---|---|
| **[Molecular AOP Builder](https://github.com/marvinm2/molAOP-builder)** | Curator-facing app that builds and approves KE→WikiPathways, KE→GO and KE→Reactome mappings. **Produces** the mapping data. | [live](https://molaop-builder.vhp4safety.nl) · [API docs](https://molaop-builder.vhp4safety.nl/api/docs) |
| **Molecular AOP Analyser** (this repo) | User-facing app for KE enrichment of gene-expression data. **Consumes** the builder's approved mappings via its public REST API. | [live](https://molaop-analyser.vhp4safety.nl) |

The integration contract between them is documented in
[`docs/KE-MAPPING-API-REFERENCE.md`](docs/KE-MAPPING-API-REFERENCE.md). Changes to the
builder's API shape must be mirrored there and in `services/api_service.py`.

Both tools are part of the [VHP4Safety](https://platform.vhp4safety.nl) platform.

---

## Features

### Data Processing & Analysis

- Upload differential expression results or use demo datasets
- **Auto-detection** of gene ID, log2FC, and p-value columns with confidence indicators
- Interactive volcano plot visualization with customizable thresholds
- **Quick threshold options**: 0, 0.5, 1.0, 1.5, 2.0, Top 10%, Top 20% genes
- KE enrichment analysis using **Fisher's exact test or GSEA**, both with FDR correction
- Support for multiple gene ID formats and duplicate handling

### Batch Analysis

- Analyse multiple datasets in a single session using the batch wizard
- Upload up to 10 files at once or select from demo datasets
- Tag each file with condition labels, timepoints, and doses
- Shared analysis settings (AOP, threshold, p-value cutoff) applied across all files
- Per-file progress tracking with real-time status updates
- **Comparison view** for cross-condition analysis with heatmap and table views

### Visualization & Interactivity

- **Interactive AOP network** visualization using Cytoscape.js with:
  - Color-coded nodes by KE type (MIE, Intermediate, AO)
  - Gene expression coloring (-2 to +2 log2FC scale)
  - Significance borders for enriched KEs and genes
  - **Comprehensive legend** with node types and color scales
  - Add/toggle gene nodes, reset view, download options
- **Enhanced results table** with user-friendly column names
- **Hub genes panel** listing genes shared across multiple Key Events in the AOP
- **Pathway view** — embedded WikiPathways diagram for a Key Event's mapped pathway
- **Network statistics** dashboard showing node/edge counts
- Proper decimal formatting for p-values (scientific notation < 0.001)

### Report Generation

- **PDF and HTML reports** with complete analysis documentation
- Experiment metadata capture (dataset ID, stressor, dosing, owner, description)
- Network visualizations embedded in reports
- Exportable tables (CSV, Excel) and network data (JSON)

---

## Architecture Overview

The application is built on a **Flask + Jinja2 + HTMX** stack:

- **Backend**: Flask (Python) handles routing, data processing, and enrichment analysis
- **Frontend**: Jinja2 templates with HTMX for dynamic partial updates (e.g., batch progress polling)
- **Visualization**: Cytoscape.js for interactive AOP network graphs, Plotly.js for volcano plots
- **Database**: SQLite via SQLAlchemy for experiment persistence and analysis history
- **Styling**: Custom CSS using the VHP4Safety house style colour palette

### Service Architecture

The backend is organized into modular services under `services/`:

| Module | Responsibility |
|---|---|
| `enrichment_service.py` | Fisher's exact test, FDR correction, enrichment statistics; dispatches to GSEA |
| `gsea_service.py` | GSEA (gseapy prerank) as a threshold-free alternative to Fisher's exact |
| `network_service.py` | Cytoscape.js network graph generation |
| `hub_service.py` | Hub-gene computation — genes shared across multiple Key Events |
| `report_service.py` | PDF and HTML report generation |
| `data_service.py` | Data loading, normalization, and validation |
| `column_detector.py` | Auto-detection of gene ID, log2FC, and p-value columns |
| `gene_id_validator.py` | Gene symbol validation and normalization |
| `batch_service.py` | Multi-file batch analysis orchestration |
| `comparison_service.py` | Cross-condition comparison and heatmap generation |
| `aop_discovery_service.py` | AOP search and typeahead suggestions |
| `sparql_service.py` | SPARQL queries for AOP-Wiki data |
| `pathway_picker_service.py` | WikiPathways picker data for the results-page Pathway view |
| `api_service.py` | REST API client for the upstream molAOP Builder (KE→pathway mappings) |

### Data Flow

1. **Reference Set Loading**: `helpers.py` merges KE-WP mappings, WP-gene edges, and node attributes into KE-to-gene dictionaries
2. **Gene Expression Processing**: User data is normalized (uppercase symbols, duplicates combined via Fisher's method for p-values)
3. **Enrichment Analysis**: Fisher's exact test on 2x2 contingency tables (significant/non-significant x in-KE/not-in-KE)
4. **Network Visualization**: Results formatted as Cytoscape.js nodes and edges for interactive AOP pathway display

---

## Getting Started

### Prerequisites

- Python 3.11+
- pip (Python package manager)
- Optionally: Docker and Docker Compose

### Local Development Setup

```bash
# 1. Clone the repository
git clone https://github.com/marvinm2/molAOP-analyser.git
cd molAOP-analyser

# 2. Create and activate a virtual environment
python -m venv venv
source venv/bin/activate        # Linux/macOS
# venv\Scripts\activate         # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Run the application
python app.py
```

The app will be available at <http://localhost:5000>.

### Docker Compose (Recommended)

```bash
git clone https://github.com/marvinm2/molAOP-analyser.git
cd molAOP-analyser
docker-compose up --build
```

Open <http://localhost:5010> in your browser (`docker-compose.yml` publishes the container's port 5000 on host port 5010).

---

## Project Structure

```
molAOP-analyser/
├── app.py                          # Flask app — routes, request handling
├── config.py                       # Configuration (upload limits, demo datasets, AOP case studies)
├── database.py                     # SQLAlchemy models and database manager (SQLite)
├── helpers.py                      # Reference set loading (KE→gene mappings)
├── cache_manager.py                # Caching layer for reference data and results
├── validation.py                   # Input validation for datasets and parameters
├── exceptions.py                   # Custom exception classes
├── utils.py                        # Utility functions
│
├── services/                       # Modular service classes
│   ├── enrichment_service.py       # Fisher's exact test, FDR correction, GSEA dispatch
│   ├── gsea_service.py             # GSEA (gseapy prerank) enrichment
│   ├── network_service.py          # Cytoscape.js network generation
│   ├── hub_service.py              # Hub-gene computation
│   ├── report_service.py           # PDF/HTML report generation
│   ├── data_service.py             # Data processing and normalization
│   ├── column_detector.py          # Auto-detect gene ID, FC, p-value columns
│   ├── gene_id_validator.py        # Gene symbol validation
│   ├── batch_service.py            # Batch analysis orchestration
│   ├── comparison_service.py       # Cross-condition comparison
│   ├── aop_discovery_service.py    # AOP search and typeahead
│   ├── sparql_service.py           # SPARQL queries for AOP-Wiki
│   ├── pathway_picker_service.py   # WikiPathways picker for Pathway view
│   └── api_service.py              # REST API client for the molAOP Builder
│
├── templates/                      # Jinja2 templates
│   ├── base.html                   # Base layout (nav, footer, HTMX)
│   ├── index.html                  # Landing page (tabbed: single + batch)
│   ├── _single_analysis.html       # Single analysis form partial
│   ├── _single_analysis_scripts.html
│   ├── _batch_analysis.html        # Batch wizard partial
│   ├── _batch_analysis_scripts.html
│   ├── results.html                # Analysis results display
│   ├── shared_results.html         # Public shareable results page
│   ├── compare.html                # Batch comparison view
│   ├── batch_progress.html         # HTMX progress partial
│   ├── batch_summary.html          # Batch results summary
│   ├── documentation.html          # In-app documentation
│   ├── demos.html                  # Curated demo dataset landing page
│   └── about.html                  # About page
│
├── static/
│   ├── css/style.css               # Application styles (VHP4Safety palette)
│   ├── css/tour.css                # Guided tour styles
│   ├── js/tour.js                  # Guided tour logic
│   └── img/                        # Logos (molaop, vhp4safety) and favicons
│
├── data/                           # Reference data and demo datasets
│   ├── aop_ke_map.csv              # AOP → KE ID mappings
│   ├── aop_ker_edges.csv           # KE relationship edges within AOPs
│   ├── KE-WP.csv                   # KE → WikiPathways ID links
│   ├── edges_wpid_to_gene.csv      # WikiPathways → gene ID mappings
│   ├── node_attributes.csv         # Gene ID → symbol mappings
│   ├── ke_metadata.csv             # KE titles and type classifications
│   ├── wikipathways_hsa_20240410.xgmml  # WikiPathways network data
│   ├── GSE90122_SR12813.tsv        # PXR agonist demo dataset 1
│   ├── GSE90122_TO90137.tsv        # PXR agonist demo dataset 2
│   └── Cisplatin_Kidney/           # 54 cisplatin toxicity datasets
│       └── CSP_{time}_{dose}.csv   # (4–72 hr, 0.1–50 uM)
│
├── tests/                          # Pytest test suite
│   ├── conftest.py                 # Fixtures (Flask client, test data)
│   ├── test_flask_routes.py        # Integration tests for web routes
│   ├── test_column_detector.py     # Unit tests for column auto-detection
│   ├── test_enrichment_service.py  # Fisher's exact test and FDR correction
│   ├── test_gsea_service.py        # GSEA enrichment
│   ├── test_hub_service.py         # Hub-gene computation
│   ├── test_database.py            # Database model tests
│   ├── test_report_service.py      # Report generation tests
│   ├── test_aop_discovery.py       # AOP search tests
│   ├── test_shared_results.py      # Shared results feature tests
│   ├── test_batch_service.py       # Batch orchestration tests
│   ├── test_data_service.py        # Data loading and normalization
│   ├── test_network_comparison.py  # Cross-condition comparison matrix
│   └── test_pathway_picker_service.py  # WikiPathways picker data
│
├── uploads/                        # Temporary user uploads (gitignored)
├── molAOP_analyser.db              # SQLite database (auto-created)
├── requirements.txt                # Python dependencies
├── Dockerfile                      # Container build definition
├── docker-compose.yml              # Service orchestration
└── pytest.ini                      # Pytest configuration
```

---

## Demo Datasets

The app includes preloaded differential expression datasets:

### PXR Agonists

- `GSE90122_TO90137.tsv` — PXR agonist TO901317 gene expression data
- `GSE90122_SR12813.tsv` — PXR agonist SR12813 gene expression data

### Cisplatin Kidney Toxicity

54 datasets covering cisplatin exposure in kidney cells across:

- **Timepoints**: 4, 8, 16, 24, 48, 72 hours
- **Doses**: 0.1, 0.5, 1, 2.5, 5, 10, 20, 30, 50 uM

Naming pattern: `CSP_{timepoint}_{dose}.csv`

---

## Testing

```bash
# Run all tests
pytest

# Run a specific test file
pytest tests/test_flask_routes.py

# Verbose output
pytest -v

# With coverage report
pytest --cov
```

### Test Files

| File | Scope |
|---|---|
| `test_flask_routes.py` | Integration tests for all web routes and the full analysis workflow |
| `test_column_detector.py` | Unit tests for auto-detection of gene ID, FC, and p-value columns |
| `test_enrichment_service.py` | Fisher's exact test, FDR correction, significance thresholds |
| `test_gsea_service.py` | GSEA enrichment (gseapy prerank) |
| `test_hub_service.py` | Hub-gene computation and ranking |
| `test_database.py` | SQLAlchemy model and database manager tests |
| `test_report_service.py` | PDF and HTML report generation |
| `test_aop_discovery.py` | AOP typeahead search functionality |
| `test_shared_results.py` | Shared/public results link feature |
| `test_batch_service.py` | Batch file parsing and orchestration |
| `test_data_service.py` | Data loading and normalization |
| `test_network_comparison.py` | Cross-condition comparison matrix |
| `test_pathway_picker_service.py` | WikiPathways picker data assembly |

---

## Deployment

### Docker Compose

The recommended deployment method uses Docker Compose:

```bash
docker-compose up --build -d
```

This builds and starts the container, which runs the Flask app (`python app.py`). `docker-compose.yml` publishes the container's port 5000 on host port 5010.

### Data Directory

The `data/` directory contains all reference datasets and must be present at runtime. It is included in the Docker image during build.

### Database

The SQLite database (`molAOP_analyser.db`) is auto-created on first run. For persistent storage in Docker deployments, mount it as a volume:

```yaml
volumes:
  - ./molAOP_analyser.db:/app/molAOP_analyser.db
```

### Environment Variables

The application reads configuration from `config.py`. Key settings:

- `MAX_CONTENT_LENGTH`: Upload size limit (default 10 MB)
- `UPLOAD_FOLDER`: Directory for temporary file storage
- `SECRET_KEY`: Flask session secret (set via environment variable in production)
