# molAOP Analyser

## What This Is

A web-based tool for analysing transcriptomic data against molecular Adverse Outcome Pathways (AOPs). Researchers and regulators upload gene expression datasets, the tool performs Over-Representation Analysis (ORA) to calculate Key Event (KE) activity, and visualises results as interactive AOP networks with progressive gene expansion. Part of the VHP4Safety platform alongside the molAOP Builder, with live API integration for KE-gene mappings.

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
- ✓ Live API integration with molAOP Builder for KE-gene set data — v1.0
- ✓ Graceful CSV fallback when Builder API unreachable — v1.0
- ✓ Cross-process disk caching of API responses — v1.0
- ✓ VHP4Safety brand palette and Inter typography on all pages — v1.0
- ✓ Sticky navigation header matching molAOP Builder style — v1.0
- ✓ Results page visual hierarchy with card layouts — v1.0
- ✓ Upload page restyled to match brand — v1.0
- ✓ Documentation page with usage guide — v1.0
- ✓ CSV export of enrichment results — v1.0
- ✓ Large network performance optimisation (per-KE progressive expansion) — v1.0
- ✓ PNG export for publication screenshots — v1.0
- ✓ Overlay zoom controls on network visualisation — v1.0

- ✓ AOP source-of-truth cleanup — hardcoded AOP:472 rows dropped, `kidney-472` manual entry removed, `KIDNEY_AOP_IDS` curation criterion documented — v3.0
- ✓ Demos page restructure — `/demos` route with curated PXR/cisplatin cards and recommended-AOP CTAs; main page leads with upload-your-data — v3.0
- ✓ Gene-by-KE CSV export with `KE_ID, KE_Title, Gene_Symbol, log2FC, pvalue, FDR` columns, respecting active thresholds — v3.0 (EXPO-04/05)
- ✓ Cytoscape network export carries per-KE p-value and FDR (closes issue #50) — v3.0 (EXPO-06)
- ✓ Column detector raw/adjusted p-value split with preview-UI picker — v3.0
- ✓ Real local file upload (multipart CSV/TSV/TXT) on upload page — v3.0 (UPUX-04)
- ✓ Tech-debt sweep: `datetime.utcnow()` deprecation removed (14 sites), duplicate `guess_id_type()` consolidated, debug prints migrated to `logger.debug()`, dead `build_cytoscape_network_legacy()` deleted — v3.0 (DEBT-01/02/03/04)

### Active

(None — pending next milestone definition)

### Shipped State after v3.0

- **GMT export** is intentionally **not** provided by the analyser. The upstream molAOP Builder is the single source of truth for KE→gene libraries and serves GMT exports directly (see `KE-MAPPING-API-REFERENCE.md`). The analyser-specific need ("which genes from THIS analysis drove each KE") is satisfied by the EXPO-04 gene-by-KE long-format CSV.

### Persistent Out of Scope

- **Species support (issue #7)** — reference data assumes human; mouse/rat support needs homology design.
- **Read-only AOP REST API (issue #16)** — AOP/KE data not unique to analyser; consumers should query AOP-Wiki SPARQL and the Builder API directly.
- **Analysis-submission API (issue #17)** — requires job-queue infrastructure; disproportionate for an interactive tool.
- **Secondary gene ID correction (issue #22)** — deferred pending mygene.info / HGNC resolver evaluation.
- **Guided tour / interactive tutorial (issue #18)** — documentation-style; defer.
- **Column auto-detection clarity (issue #21)** — small UX refinement; defer.

### Out of Scope

- Report generation improvements — not a priority for this round
- OAuth/user authentication — not needed for single-analyst tool; SQLite persistence covers reproducibility
- Real-time collaboration — single-user analysis tool; WebSocket complexity disproportionate
- GO term analysis — builder handles GO mappings separately
- AI/LLM analysis panel — hallucination risk unacceptable for regulatory tool
- GSEA (ranked list enrichment) — different statistical assumptions; ORA appropriate for binary lists
- Multi-omics input — no KE-gene mappings exist for non-transcriptomic identifiers
- Full pathway database browsing (KEGG/Reactome) — dilutes AOP-specific value proposition
- Command-line/API batch analysis — tool optimised for interactive regulatory review
- Offline mode — real-time API integration is core value

## Context

- **Shipped through v3.0** — ~10,439 LOC Python (146 tests passing)
- **Tech stack**: Flask, SQLAlchemy, Cytoscape.js 3.33.1, dagre, DataTables, Plotly, HTMX (batch progress)
- **Architecture**: Modular service layer (api_service, enrichment_service, network_service, report_service, data_service, column_detector, gene_id_validator, batch_service, comparison_service, aop_discovery_service, sparql_service)
- **API integration**: Builder API client with paginated fetching, exponential backoff, diskcache.FanoutCache for cross-process safety
- **Data flow**: API-first → disk cache → CSV fallback (three-tier with data source indicator); SPARQL queries to AOP-Wiki for dynamic AOP discovery
- **Network visualisation**: KE-only initial render with dagre LR layout, per-KE click-to-expand, shared gene dedup, overlay zoom, transparent PNG export
- **Comparison view**: Multi-condition heatmap, comparison table, multi-condition AOP network with pie-chart KE nodes, delta mode
- **Batch analysis**: 2-10 file sequential queue, SQLite-persisted state, harmonised backgrounds, HTMX progress polling
- **Demos**: `/demos` route with curated cards (PXR + cisplatin) and recommended-AOP CTAs; main `/` is upload-first
- **v3.0 cleanup**: All `datetime.utcnow()` migrated to timezone-aware; debug prints migrated to `logger.debug`; dead code removed; duplicate helpers consolidated

## Constraints

- **Tech stack**: Flask/Python — maintain consistency with existing codebase and builder
- **API dependency**: molAOP Builder must be running/accessible for live data (CSV fallback available)
- **Browser support**: Modern browsers (Cytoscape.js requirement)
- **File size**: 10MB upload limit for gene expression datasets
- **Brand**: VHP4Safety house style colours and Inter typography must be used

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Live API calls over periodic sync | Real-time data ensures analyser always uses latest KE-gene mappings | ✓ Good — works with diskcache layer |
| Keep demo datasets local | Demo data is static reference data, doesn't change with builder updates | ✓ Good |
| Match builder UI aesthetic | Visual consistency across VHP4Safety platform tools | ✓ Good — all pages use shared base.html |
| Deprioritise report generation | Focus on core analysis workflow and UI improvements first | ✓ Good — report generation still works, improvements deferred |
| Keep SQLite persistence | Experiment history useful for reproducibility | ✓ Good |
| diskcache.FanoutCache for cross-process caching | Avoids stale data across Gunicorn workers; file-backed survives restarts | ✓ Good — 8 shards, shared directory |
| Per-KE progressive gene expansion | Avoids browser freeze on large networks (AOP:2 ~696 genes) | ✓ Good — click-to-expand with shared gene dedup |
| Neutral grey background (#f5f6fa) | Sky blue (#93D5F6) too dominant with branded dark nav | ✓ Good |
| Three-section documentation page | Getting Started/FAQ removed to focus on essential content | ✓ Good |
| AOP:472 served from SPARQL only (no CSV fallback rows) | Single source of truth avoids drift; AOP:472 behaves like the other 17 dynamic kidney AOPs | ✓ Good (v3.0, Phase 10) |
| `KIDNEY_AOP_IDS` kept as curated static list with documented criterion | Dynamic derivation from AOP-Wiki would add a network dependency at startup; curation criterion is the meaningful contract | ✓ Good (v3.0, Phase 10) |
| GMT export served by Builder, not analyser | Builder is the upstream source of truth for KE→gene mappings; analyser would duplicate without value | ✓ Good (v3.0, Phase 11) |
| Gene-by-KE CSV (long format) over GMT for analyser exports | Long format carries per-gene log2FC/p-value/FDR — strictly more interrogable than intersection-style GMT | ✓ Good (v3.0, Phase 11) |
| Demos moved to dedicated `/demos` route with recommended-AOP CTAs | Cleans the main upload flow and surfaces curated AOP+demo pairings | ✓ Good (v3.0, Phase 10.1) |
| Multipart file upload over drag-drop-only widget | File-input is more accessible and works the same way as the demo flow | ✓ Good (v3.0, Phase 12) |

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd-transition`):
1. Requirements invalidated? → Move to Out of Scope with reason
2. Requirements validated? → Move to Validated with phase reference
3. New requirements emerged? → Add to Active
4. Decisions to log? → Add to Key Decisions
5. "What This Is" still accurate? → Update if drifted

**After each milestone** (via `/gsd-complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---
*Last updated: 2026-05-11 after v3.0 milestone shipped*
