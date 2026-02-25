# molAOP Analyser

## What This Is

A web-based tool for analysing transcriptomic data against molecular Adverse Outcome Pathways (AOPs). Researchers and regulators upload gene expression datasets, the tool performs Over-Representation Analysis (ORA) to calculate Key Event (KE) activity, and visualises results as interactive AOP networks. Part of the VHP4Safety platform alongside the molAOP Builder.

## Core Value

Researchers and regulators can quickly determine which Key Events in a molecular AOP are activated by their transcriptomic data, enabling evidence-based chemical safety assessment.

## Requirements

### Validated

- ✓ File upload for gene expression datasets (CSV/TSV/TXT, 10MB limit) — existing
- ✓ Demo dataset selection (PXR agonists, cisplatin kidney) — existing
- ✓ Auto-detection of gene ID, log2FC, and p-value columns with confidence scoring — existing
- ✓ Interactive threshold selection for log2FC and p-value — existing
- ✓ AOP case study selection (liver steatosis, kidney failure, organ networks) — existing
- ✓ Fisher's exact test enrichment analysis with FDR correction — existing
- ✓ Cytoscape.js interactive network visualisation — existing
- ✓ Enrichment results table with formatted p-values — existing
- ✓ Volcano plot visualisation — existing
- ✓ Experiment metadata capture and SQLite persistence — existing
- ✓ PDF/HTML report generation — existing
- ✓ CSRF protection and input validation — existing

### Active

- [ ] Live API integration with molAOP Builder for KE-gene set data
- [ ] UI/UX overhaul to match VHP4Safety/molAOP Builder aesthetic (Inter font, brand palette, component patterns)
- [ ] Improved network visualisation — better interactivity and performance for large networks
- [ ] Documentation page with usage guide
- [ ] Navigation component matching molAOP Builder style (sticky header, app branding)

### Out of Scope

- Report generation improvements — not a priority for this round
- OAuth/user authentication — not needed for analyser
- Real-time collaboration — single-user analysis tool
- GO term analysis — builder handles GO mappings separately

## Context

- **Existing codebase**: Fully functional Flask application with modular service architecture, SQLAlchemy persistence, Cytoscape.js visualisation
- **molAOP Builder**: Sister application providing the `/api/v1/mappings` REST API for KE-WP mapping data, KE gene sets via `/ke_genes/<ke_id>`, and AOP structure via `/get_aop_kes/<aop_id>`
- **Builder API features**: Paginated responses, CSV export, AOP filtering, CORS-enabled
- **Current data flow**: Static CSV files (aop_ke_map.csv, KE-WP.csv, edges_wpid_to_gene.csv, node_attributes.csv) provide reference data — these will be replaced by API calls
- **Demo datasets**: Local gene expression files (GSE90122_*.tsv, Cisplatin_Kidney/) remain local
- **Performance concern**: Large networks with dozens of KEs and their gene sets can cause Cytoscape.js slowdowns
- **Brand consistency**: Both tools share the VHP4Safety house style colour palette (#29235C primary dark, #E6007E magenta, #307BBF blue)

## Constraints

- **Tech stack**: Flask/Python — maintain consistency with existing codebase and builder
- **API dependency**: molAOP Builder must be running/accessible for live data
- **Browser support**: Modern browsers (Cytoscape.js requirement)
- **File size**: 10MB upload limit for gene expression datasets
- **Brand**: VHP4Safety house style colours and Inter typography must be used

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Live API calls over periodic sync | Real-time data ensures analyser always uses latest KE-gene mappings | — Pending |
| Keep demo datasets local | Demo data is static reference data, doesn't change with builder updates | — Pending |
| Match builder UI aesthetic | Visual consistency across VHP4Safety platform tools | — Pending |
| Deprioritise report generation | Focus on core analysis workflow and UI improvements first | — Pending |
| Keep SQLite persistence | Experiment history useful for reproducibility | — Pending |

---
*Last updated: 2026-02-25 after initialization*
