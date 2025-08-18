# Molecular AOP KE Enrichment App

This web application allows users to upload or select gene expression datasets and perform Key Event (KE) enrichment analysis in the context of Molecular Adverse Outcome Pathways (AOPs). The results are visualized in interactive tables and network diagrams with comprehensive reporting capabilities.

---

## Features

### Data Processing & Analysis
- Upload differential expression results or use demo datasets  
- **Auto-detection** of gene ID, log2FC, and p-value columns with confidence indicators
- Interactive volcano plot visualization with customizable thresholds
- **Quick threshold options**: 0, 0.5, 1.0, 1.5, 2.0, Top 10%, Top 20% genes
- KE enrichment analysis using Fisher's exact test with FDR correction
- Support for multiple gene ID formats and duplicate handling

### Visualization & Interactivity
- **Interactive AOP network** visualization using Cytoscape.js with:
  - Color-coded nodes by KE type (MIE, Intermediate, AO)
  - Gene expression coloring (-2 to +2 log2FC scale)
  - Significance borders for enriched KEs and genes
  - **Comprehensive legend** with node types and color scales
  - Add/toggle gene nodes, reset view, download options
- **Enhanced results table** with user-friendly column names
- **Network statistics** dashboard showing node/edge counts
- Proper decimal formatting for p-values (scientific notation < 0.001)

### Report Generation
- **PDF and HTML reports** with complete analysis documentation
- Experiment metadata capture (dataset ID, stressor, dosing, owner, description)
- Network visualizations embedded in reports
- Exportable tables (CSV, Excel) and network data (JSON)

---

## Getting Started (Docker Compose)

### 1. Clone the repository

```bash
git clone https://github.com/marvinm2/molAOP-analyser.git
cd molAOP-analyser
```

### 2. Build and run the app
```bash
docker-compose up --build
```

### 3. Open in your browser
Go to: http://localhost:5000

---

## Project Structure
```bash
├── app.py                      # Flask app entry point
├── helpers.py                 # Utility functions
├── templates/                 # HTML templates
│   ├── index.html
│   └── results.html
├── static/
│   ├── css/style.css
│   └── img/logo.png
├── data/                      # Input mapping and demo datasets
│   └── *.csv, *.tsv, *.xgmml
├── uploads/                   # Temporary user uploads
├── requirements.txt           # Python dependencies
├── Dockerfile                 # Container build
└── docker-compose.yml         # Service orchestration
```
---

## Demo Datasets

The app includes preloaded differential expression datasets simulating chemical exposure:

`GSE90122_TO90137.tsv`: PXR agonist 1

`GSE90122_SR12813.tsv`: PXR agonist 2

---

## Development

To run locally without Docker:
```bash
pip install -r requirements.txt
python app.py
```
