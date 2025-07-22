# Molecular AOP KE Enrichment App

This web application allows users to upload or select gene expression datasets and perform Key Event (KE) enrichment analysis in the context of Molecular Adverse Outcome Pathways (AOPs). The results are visualized in interactive tables and network diagrams.

---

## Features

- Upload differential expression results or use demo datasets
- Column mapping for gene symbols, log2FC, and p-values
- Volcano plot visualization of DE genes
- KE enrichment analysis using Fisher's exact test
- AOP network visualization using Cytoscape.js
- Exportable results and interactive controls

---

## Getting Started (Docker Compose)

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/aop-ke-enrichment.git
cd aop-ke-enrichment
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
