# Requirements: molAOP Analyser

**Defined:** 2026-02-27
**Core Value:** Researchers and regulators can quickly determine which Key Events in a molecular AOP are activated by their transcriptomic data, enabling evidence-based chemical safety assessment.

## v2.0 Requirements

Requirements for v2.0 Batch Analysis & Dynamic AOPs milestone. Each maps to roadmap phases.

### AOP Discovery

- [x] **AOPD-01**: User can search AOPs with curated KE-gene mappings via free-text typeahead
- [x] **AOPD-02**: Dynamic AOP list fetched from Builder API (only AOPs with available KE-gene mappings)
- [x] **AOPD-03**: Curated (hardcoded) AOPs remain available alongside dynamic list
- [x] **AOPD-04**: Three-tier fallback for AOP data (cache → Builder API → hardcoded config)
- [x] **AOPD-05**: KE metadata (titles, types) enriched from SPARQL for selected AOP
- [x] **AOPD-06**: Builder API KE-gene mapping integration works for any AOP in the dynamic list

### Shareable URLs

- [x] **SHAR-01**: User can access saved analysis results via `/results/<uuid>` URL
- [x] **SHAR-02**: User can copy shareable link to clipboard from results page
- [x] **SHAR-03**: Shared results page displays enrichment table and network (without volcano plot)

### Network Interactivity

- [x] **NETW-01**: User can expand/collapse all gene nodes with one toggle button
- [x] **NETW-02**: Expand-all toggle disabled when total gene count exceeds 200
- [x] **NETW-03**: User can view KE detail panel on node click (title, type, FDR, gene list, AOP-Wiki link)
- [x] **NETW-04**: User can filter network by FDR significance threshold via slider

### Upload & UX

- [x] **UPUX-01**: Responsive layout at 768px and 1024px breakpoints
- [x] **UPUX-02**: Drag-and-drop file upload hint on upload page
- [x] **UPUX-03**: Network and results layout stacks vertically on tablet viewports

### Batch Analysis

- [x] **BATC-01**: User can upload multiple files (2-10) in a sequential queue
- [x] **BATC-02**: User can tag each file with condition label (dose, timepoint, treatment)
- [x] **BATC-03**: Per-file progress indicator during batch analysis
- [x] **BATC-04**: Batch state persisted in SQLite (not session cookies)
- [x] **BATC-05**: Background gene set harmonisation across conditions for valid comparison

### Comparison Visualisation

- [x] **COMP-01**: KE × condition heatmap showing -log10(FDR) significance
- [x] **COMP-02**: Comparison enrichment table (rows = KEs, columns = conditions)
- [x] **COMP-03**: Multi-condition AOP network overlay on shared KE skeleton
- [x] **COMP-04**: User can toggle between single-condition and comparison view
- [x] **COMP-05**: Delta mode showing change relative to user-selected reference condition
- [x] **COMP-06**: Condition colour legend using VHP4Safety palette

## v3.0 Requirements

Requirements for v3.0 Cleanup, Gene-Set Export & UX Gaps milestone.

### AOP Source-of-Truth Cleanup

- [x] **AOPD-07**: Hardcoded AOP:472 KE rows are removed from `data/aop_ke_map.csv` and `data/aop_ker_edges.csv`; AOP:472 data flows from SPARQL only
- [x] **AOPD-08**: `kidney-472` manual entry is removed from `Config.CASE_STUDY_AOPS`; AOP:472 surfaces only via dynamic discovery alongside the other kidney AOPs
- [x] **AOPD-09**: `KIDNEY_AOP_IDS` is either derived dynamically from AOP-Wiki (e.g. AOs related to renal failure) or explicitly justified in code comments as a curated static list with the curation criterion documented

### Gene-Set Export

- [ ] **EXPO-03**: User can export significant genes per KE in GMT format from the results page (one row per KE: `KE_ID<tab>description<tab>gene1<tab>gene2...`), suitable for GSEA/Enrichr import
- [ ] **EXPO-04**: User can export a gene-by-KE mapping table as CSV with columns `KE_ID, KE_Title, Gene_Symbol, log2FC, pvalue, FDR`
- [ ] **EXPO-05**: Gene-set exports respect the active significance thresholds (only genes meeting log2FC and p-value cutoffs are included)
- [ ] **EXPO-06**: Cytoscape network export (PNG metadata and any data-export JSON) includes per-KE p-value and FDR (issue #50)

### Upload UX

- [ ] **UPUX-04**: User can choose and upload a local CSV/TSV/TXT file from the upload page (currently demo-only despite the page advertising "upload your own data")

### Tech Debt

- [ ] **DEBT-01**: All `datetime.utcnow()` call sites use timezone-aware `datetime.now(datetime.UTC)` instead (database.py, app.py, tests); pytest run produces no `DeprecationWarning` from datetime
- [ ] **DEBT-02**: Duplicate `guess_id_type()` is removed from `app.py`; `services/data_service.py` is the single definition
- [ ] **DEBT-03**: Debug `print()` statements in `helpers.py` are replaced with `logger.debug()`
- [ ] **DEBT-04**: Unused `build_cytoscape_network_legacy()` is removed from `services/network_service.py`

## Future Requirements

Deferred to future release. Tracked but not in current roadmap.

### Deferred Interactivity

- **NETW-05**: Animated time-course playback across ordered conditions
- **NETW-06**: KE description tooltips from AOP-Wiki SPARQL on hover

### Deferred Export

- **EXPO-01**: Batch comparison CSV export (full heatmap matrix)
- **EXPO-02**: Network SVG export for publication use

### Deferred — Multi-Species & Data Quality

- **SPEC-01**: Global species setting (mouse/rat/human) so non-human transcriptomic datasets can be analysed (issue #7) — needs design (homology mapping vs species-specific reference sets)
- **GIDQ-01**: Secondary/withdrawn gene ID correction on upload via mygene.info or HGNC resolver (issue #22)
- **TUTR-01**: Interactive guided tour for new users (issue #18)
- **AUTC-01**: Improved clarity around column auto-detection (issue #21)

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Parallel batch analysis | No job queue infrastructure; sequential is reliable for 2-10 files |
| Per-condition separate Cytoscape.js instances | Memory and sync fragility; single instance with data properties is performant |
| Automatic condition detection from filenames | Fragile parsing; require user confirmation instead |
| Real-time WebSocket progress | Over-engineered for 5-50s batch runs; htmx polling sufficient |
| Persistent cross-session batch history | Requires user accounts; out of scope per constraints |
| Volcano plot on shared results page | Gene expression data not stored in DB; enrichment + network only |
| Read-only AOP REST API (issue #16) | AOP/KE data is not unique to the analyser — consumers should query AOP-Wiki SPARQL and the Builder API directly. A proxy here adds latency without adding value. |
| Analysis-submission REST API (issue #17) | Requires job-queue infrastructure (Celery/Redis); disproportionate for an interactive regulatory tool. |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| AOPD-01 | Phase 4 | Complete |
| AOPD-02 | Phase 4 | Complete |
| AOPD-03 | Phase 4 | Complete |
| AOPD-04 | Phase 4 | Complete |
| AOPD-05 | Phase 4 | Complete |
| AOPD-06 | Phase 4 | Complete |
| SHAR-01 | Phase 4 | Complete |
| SHAR-02 | Phase 4 | Complete |
| SHAR-03 | Phase 4 | Complete |
| NETW-01 | Phase 5 | Complete |
| NETW-02 | Phase 5 | Complete |
| NETW-03 | Phase 5 | Complete |
| NETW-04 | Phase 5 | Complete |
| UPUX-01 | Phase 5 | Complete |
| UPUX-02 | Phase 5 | Complete |
| UPUX-03 | Phase 5 | Complete |
| BATC-01 | Phase 6 | Complete |
| BATC-02 | Phase 6 | Complete |
| BATC-03 | Phase 6 | Complete |
| BATC-04 | Phase 6 | Complete |
| BATC-05 | Phase 6 | Complete |
| COMP-01 | Phase 7 | Complete |
| COMP-02 | Phase 7 | Complete |
| COMP-03 | Phase 7 | Complete |
| COMP-04 | Phase 7 | Complete |
| COMP-05 | Phase 7 | Complete |
| COMP-06 | Phase 7 | Complete |
| AOPD-07 | Phase 10 | Planned |
| AOPD-08 | Phase 10 | Planned |
| AOPD-09 | Phase 10 | Planned |
| EXPO-03 | Phase 11 | Planned |
| EXPO-04 | Phase 11 | Planned |
| EXPO-05 | Phase 11 | Planned |
| EXPO-06 | Phase 11 | Planned |
| UPUX-04 | Phase 12 | Planned |
| DEBT-01 | Phase 12 | Planned |
| DEBT-02 | Phase 12 | Planned |
| DEBT-03 | Phase 12 | Planned |
| DEBT-04 | Phase 12 | Planned |

**Coverage:**
- v2.0 requirements: 27 total — mapped 27, unmapped 0 ✓
- v3.0 requirements: 12 total — mapped 12, unmapped 0 ✓

---
*Requirements defined: 2026-02-27*
*Last updated: 2026-05-06 — v3.0 traceability added after roadmap creation*
