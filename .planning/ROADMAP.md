# Roadmap: molAOP Analyser

## Milestones

- ✅ **v1.0 MVP** — Phases 1-3 (shipped 2026-02-27)
- ✅ **v2.0 Batch Analysis & Dynamic AOPs** — Phases 4-9 (completed 2026-03-02)
- 🚧 **v3.0 Cleanup, Gene-Set Export & UX Gaps** — Phases 10-12 (in progress, started 2026-05-06)

## Phases

<details>
<summary>✅ v1.0 MVP (Phases 1-3) — SHIPPED 2026-02-27</summary>

- [x] Phase 1: API Integration (2/2 plans) — completed 2026-02-25
- [x] Phase 2: UI/UX and Documentation (3/3 plans) — completed 2026-02-25
- [x] Phase 3: Network Performance (2/2 plans) — completed 2026-02-25

</details>

### v2.0 Batch Analysis & Dynamic AOPs

**Milestone Goal:** Enable multi-dataset comparison and dynamic AOP discovery, transforming the tool from single-dataset analysis to cross-condition investigation.

- [x] **Phase 4: Shareable URLs + AOP Discovery** - UUID-based result sharing and SPARQL-driven dynamic AOP list replace hardcoded config (completed 2026-02-27)
- [x] **Phase 5: Network Interactivity + Responsive Layout** - Expand-all toggle, KE detail panel, FDR filter slider, and responsive CSS for tablet/laptop (completed 2026-02-28)
- [x] **Phase 6: Batch Analysis** - Sequential multi-file upload with condition tagging, SQLite state persistence, and harmonised background gene sets (completed 2026-03-01)
- [x] **Phase 7: Comparison Visualisation** - Side-by-side heatmap, comparison table, multi-condition AOP network overlay, and delta mode (completed 2026-03-02)
- [x] **Phase 8: Comparison UI Polish** - KE drawer unification, binary pie charts, delta mode banner (completed 2026-03-02)
- [x] **Phase 9: Deployment Hardening** - Batch nav link + DB init under WSGI (gap closure) (completed 2026-03-02)

### v3.0 Cleanup, Gene-Set Export & UX Gaps

**Milestone Goal:** Pay down post-v2.0 tech debt, drop dead AOP-472 fallback data, close two long-standing UX gaps, and ship gene-set export so users can take their KE-conditioned significant-gene sets into downstream pathway-analysis tools.

- [x] **Phase 10: AOP Source-of-Truth Cleanup** - Drop hardcoded AOP:472 rows and `kidney-472` config entry; document or dynamically derive `KIDNEY_AOP_IDS` (completed 2026-05-06)
- [ ] **Phase 10.1: Demos Page Restructure** *(INSERTED)* - Move demos to `/demos` with recommended-AOP options per demo type; main page leads with upload-your-data
- [x] **Phase 11: Gene-Set Export** - GMT and CSV export of significant genes per KE on the results page; per-KE p-value/FDR embedded in Cytoscape network exports — completed 2026-05-10 (EXPO-03 deferred-to-builder)
- [ ] **Phase 12: Upload Widget + Tech-Debt Sweep** - Real local file upload on upload page; datetime deprecation cleanup, dedup `guess_id_type()`, replace debug prints, drop dead legacy network builder

## Phase Details

### Phase 4: Shareable URLs + AOP Discovery
**Goal**: Users can share analysis results via a stable URL and select any AOP from a live dynamic list rather than a hardcoded menu
**Depends on**: Phase 3 (v1.0 complete)
**Requirements**: SHAR-01, SHAR-02, SHAR-03, AOPD-01, AOPD-02, AOPD-03, AOPD-04, AOPD-05, AOPD-06
**Success Criteria** (what must be TRUE):
  1. User can open `/results/<uuid>` in a browser and see the enrichment table and network from a previously run analysis
  2. User can click a "Copy link" button on the results page and paste a working URL into any browser
  3. User can type into a typeahead search box on the upload page and see matching AOPs from the live Builder API list
  4. The upload page AOP list loads and displays correctly even when the AOP-Wiki SPARQL endpoint is unreachable (hardcoded AOPs still appear)
  5. KE metadata (title, type) displayed in the results network is populated from SPARQL for any AOP selected from the dynamic list
**Plans**: TBD

Plans:
- [x] 04-01: UUID identity — migrate ExperimentRecord to UUID primary keys, add `/results/<uuid>` route and shared results template
- [x] 04-02: AOP discovery service — `aop_discovery_service.py` with three-tier fallback (cache → SPARQL → hardcoded), `GET /api/aops` endpoint, typeahead UI on upload page

### Phase 5: Network Interactivity + Responsive Layout
**Goal**: Users on tablet and laptop screens have a fully usable interface, and network interaction is richer with one-click gene expansion, detail panels, and significance filtering
**Depends on**: Phase 4
**Requirements**: NETW-01, NETW-02, NETW-03, NETW-04, UPUX-01, UPUX-02, UPUX-03
**Success Criteria** (what must be TRUE):
  1. User can click "Expand all" and all gene nodes appear at once; the button is disabled and labelled "Too many genes" when total genes exceed 200
  2. User clicks a KE node and a detail panel opens showing the KE title, type, FDR value, gene list, and a link to AOP-Wiki
  3. User drags an FDR threshold slider and non-significant KE nodes are immediately hidden from the network
  4. The upload page and results page display without horizontal scrolling on a 768px wide screen, with the network and results table stacked vertically
  5. The upload page shows a drag-and-drop target area with visual feedback when a file is dragged over it
**Plans**: TBD

Plans:
- [ ] 05-01: Network interactivity — expand-all toggle with gene count guard, KE detail drawer, FDR filter slider in `results.html` and network JS
- [ ] 05-02: Responsive layout + upload UX — CSS media queries at 768px and 1024px breakpoints, vertical stacking for tablet, drag-and-drop hint on `index.html`

### Phase 6: Batch Analysis
**Goal**: Users can upload 2-10 gene expression files, tag each with a condition label, and run enrichment sequentially, with all state persisted in SQLite and gene backgrounds harmonised for valid cross-condition comparison
**Depends on**: Phase 5
**Requirements**: BATC-01, BATC-02, BATC-03, BATC-04, BATC-05
**Success Criteria** (what must be TRUE):
  1. User can add up to 10 files to an upload queue on the batch page, assign a condition label to each, and start analysis with one button
  2. A per-file progress indicator updates as each file is processed, without the user needing to reload the page
  3. If the user closes and reopens the batch status page, the analysis state and results are still present (not lost on session expiry)
  4. After batch analysis completes, all condition results are comparable because gene backgrounds have been harmonised across the dataset files
**Plans**: TBD

Plans:
- [x] 06-01: Batch data model — `BatchRecord` and `ConditionRecord` SQLAlchemy models, session-UUID scoped upload directories, SQLite state persistence
- [x] 06-02: Batch orchestration — `batch_service.py` sequential enrichment loop, background gene set harmonisation, `batch.html` upload wizard
- [x] 06-03: Batch progress — htmx polling against `/batch/<batch_id>/status`, per-file progress indicator, summary page, condition results drill-through

### Phase 7: Comparison Visualisation
**Goal**: Users can see all batch conditions side-by-side on a shared AOP network and heatmap, toggle between absolute and delta modes, and understand relative significance across conditions at a glance
**Depends on**: Phase 6
**Requirements**: COMP-01, COMP-02, COMP-03, COMP-04, COMP-05, COMP-06
**Success Criteria** (what must be TRUE):
  1. User sees a Plotly heatmap with KEs as rows and conditions as columns, coloured by -log10(FDR) significance on the comparison page
  2. User sees a comparison enrichment table with one row per KE and one column per condition, showing FDR values
  3. User can click a toggle to switch between single-condition and comparison view on the AOP network, with KE nodes coloured by the condition palette
  4. User can select a reference condition and the heatmap and network switch to delta mode showing change relative to that baseline
  5. A condition colour legend using VHP4Safety palette colours is visible on the comparison page
**Plans**: 3 plans

Plans:
- [x] 07-01: Comparison service + route — `comparison_service.py` pandas pivot, Flask route, `compare.html` tabbed scaffold with control bar and VHP4Safety condition legend
- [x] 07-02: Heatmap + table — Plotly.js heatmap (KE x condition, -log10 FDR), sortable comparison table, absolute/delta mode for both views
- [x] 07-03: Network overlay — Cytoscape.js pie-chart KE nodes, single/comparison toggle, KE detail panel with bar chart, delta mode for network

### Phase 8: Comparison UI Polish

**Goal:** Polish the comparison visualisation UI: fix KE detail drawer sizing/overlap and missing AOP-Wiki link, reconsider pie chart slice sizing to use equal percentage per condition with FDR < 0.05 threshold, and improve delta mode network clarity
**Requirements**: POLISH-01, POLISH-02, POLISH-03
**Depends on:** Phase 7
**Plans:** 2/2 plans executed (COMPLETE)

Plans:
- [x] 08-01: Drawer CSS/HTML unification — migrate compare drawer to .ke-drawer class hierarchy, responsive width override, AOP-Wiki link styling
- [x] 08-02: Binary pie charts + delta mode — equal-sized binary pie slices (FDR < 0.05), negative-delta gray slices, delta mode banner

### Phase 9: Deployment Hardening

**Goal:** Close integration gaps from v2.0 audit — add batch analysis navigation entry point and fix database initialisation to work under any WSGI server (gunicorn, waitress, etc.)
**Requirements**: BATC-01 (discoverability), SHAR-01, BATC-04 (deployment robustness)
**Depends on:** Phase 8
**Gap Closure:** Closes NAV-BATCH and DB-INIT-GUNICORN from v2.0-MILESTONE-AUDIT.md
**Plans:** 1/1 plans complete

### Phase 10: AOP Source-of-Truth Cleanup
**Goal**: AOP:472 data flows from a single dynamic source (SPARQL + dynamic discovery) and `KIDNEY_AOP_IDS` is either derived dynamically or has its curation criterion documented in code
**Depends on**: Phase 9
**Requirements**: AOPD-07, AOPD-08, AOPD-09
**Success Criteria** (what must be TRUE):
  1. Selecting AOP:472 from the dropdown loads its KEs and KERs without any rows in `data/aop_ke_map.csv` or `data/aop_ker_edges.csv` referencing AOP:472
  2. `Config.CASE_STUDY_AOPS` contains no `kidney-472` entry, yet AOP:472 still appears in the dynamic AOP list alongside the other kidney AOPs when the SPARQL endpoint is reachable
  3. `KIDNEY_AOP_IDS` in `config.py` is either populated by a documented dynamic query (e.g. AOs related to renal failure) or is annotated with an explicit comment stating the curation criterion and how to refresh it
  4. Existing kidney case-study workflows (single + batch analysis on cisplatin demo) still produce the same enrichment results after the cleanup (no regression in KE coverage)
**Plans**: 2 plans

Plans:
- [x] 10-01-PLAN.md — Remove AOP:472 CSV rows + kidney-472 dict entry; annotate KIDNEY_AOP_IDS with curation provenance (autonomous) — completed 2026-05-06
- [x] 10-02-PLAN.md — Manual UAT regression gate on cisplatin demo with cache-clear discipline (checkpoint)

### Phase 10.1: Demos Page Restructure (INSERTED)

**Goal:** Demos move off the main page into a curated `/demos` page with recommended AOP options per demo type. The main page leads with "upload your data" and links to demos. The recommended-AOP mapping covers AOP:1 (PXR liver steatosis) for PXR demos and the curated kidney AOPs for cisplatin demos.
**Depends on:** Phase 10
**Requirements**: TBD (defined during discuss-phase)
**Success Criteria** (what must be TRUE):
  1. Visiting `/` shows an upload-your-data layout with no demo grid; a discoverable "Try a demo" link/CTA points to `/demos`
  2. Visiting `/demos` shows curated demo cards (PXR agonists + Cisplatin kidney exposure) and CTAs that prefill `/preview` with the dataset and a recommended AOP
  3. When the user lands on `/preview` via a demo CTA, the AOP picker shows that demo's recommended AOPs by default; a "Show all AOPs" toggle reveals the full list
  4. AOP:1 (PXR activation leading to liver steatosis) appears in the recommended set for the PXR demos
  5. The Kidney-aop-network and the curated single-AOP kidney entries appear in the recommended set for the Cisplatin demos
  6. No regression in the existing single-analysis flow when the user uploads their own data — the AOP picker is unrestricted in that path
**Open design questions** (resolved during discuss-phase, not pre-locked):
  - Cisplatin file count: keep all 54 dose×timepoint files, or curate a representative subset?
  - Liver combined network: leave `Liver-aop-network` stubbed in `config.py:146`, or curate a `LIVER_AOP_IDS` list as part of this phase?
  - Soft-vs-hard AOP restriction: "show all" toggle, or hide non-recommended AOPs entirely on the demo path?
  - Main-page CTA wording for demo discoverability
**Plans:** 3 plans

Plans:
- [x] 10.1-01-PLAN.md — Backend foundation: Config.DEMO_AOP_RECOMMENDATIONS, /demos route, demos.html, nav entry, /preview accepts recommended_aops (autonomous) — completed 2026-05-06
- [x] 10.1-02-PLAN.md — Frontend restructure: remove demo block, add banner, AOP picker filter toggle with sticky behavior, CSS modifiers (autonomous) — completed 2026-05-06
- [ ] 10.1-03-PLAN.md — Manual UAT walk-through (checkpoint:human-verify)

### Phase 11: Gene-Set Export
**Goal**: Users can export per-KE significant gene sets in formats consumable by downstream pathway-analysis tools, and Cytoscape network exports carry the per-KE significance metrics needed to interpret a static snapshot
**Depends on**: Phase 10
**Requirements**: EXPO-04, EXPO-05, EXPO-06 (EXPO-03 deferred-to-builder, see REQUIREMENTS.md)
**Success Criteria** (what must be TRUE):
  1. User clicks "Export gene-by-KE table (CSV)" on the results page and downloads a CSV with columns `KE_ID, KE_Title, Gene_Symbol, log2FC, pvalue, FDR`
  2. The export contains only genes that meet the active log2FC and p-value thresholds set on the analysis page (verifiable by re-running with stricter thresholds and confirming the file shrinks)
  3. The Cytoscape network data export (PNG metadata and any JSON/data export) includes per-KE p-value and FDR fields, so a downstream consumer can reconstruct significance without re-running the analysis (closes issue #50)
**Plans**: 4 plans
**UI hint**: yes

Plans:
- [x] 11-01-PLAN.md — Column detector raw vs adjusted p-value split (D-05) — services/column_detector.py + tests (autonomous, wave 1) — completed 2026-05-10
- [x] 11-02-PLAN.md — Embed per-KE p_value + fdr in cytoscape network export (D-07, EXPO-06, closes issue #50) — services/network_service.py + tests (autonomous, wave 1) — completed 2026-05-10
- [x] 11-03-PLAN.md — Gene-by-KE CSV export button + JS helper + Builder pointer (D-01/D-04/D-09/D-10/D-11/D-12, EXPO-04 + EXPO-05) — enrichment_service.py, app.py, results.html, shared_results.html + tests (autonomous, wave 2, depends on 11-01) — completed 2026-05-10
- [x] 11-04-PLAN.md — Requirement-doc cleanup: mark EXPO-03 deferred-to-builder, remove ROADMAP success criterion #1, add documentation pointer (autonomous, wave 1) — completed 2026-05-10

### Phase 12: Upload Widget + Tech-Debt Sweep
**Goal**: The upload page lets users actually upload their own CSV/TSV/TXT files, and the codebase is free of the long-standing tech-debt items (deprecated `datetime.utcnow()` calls, duplicate helpers, debug prints, dead legacy code)
**Depends on**: Phase 11
**Requirements**: UPUX-04, DEBT-01, DEBT-02, DEBT-03, DEBT-04
**Success Criteria** (what must be TRUE):
  1. User on the upload page can pick a local CSV/TSV/TXT file via a file-input or drag-and-drop, submit it, and proceed to the threshold/preview step exactly as the demo flow does today
  2. Running `pytest` produces zero `DeprecationWarning` entries originating from `datetime.utcnow()` (all call sites use `datetime.now(datetime.UTC)` or equivalent timezone-aware constructors)
  3. `grep -nR "def guess_id_type" .` returns exactly one definition (in `services/data_service.py`); `app.py` imports the helper rather than re-defining it
  4. `grep -nR "print(" helpers.py` returns no debug prints; equivalent `logger.debug()` calls exist where useful
  5. `build_cytoscape_network_legacy()` no longer exists in `services/network_service.py`, and the existing test suite still passes
**Plans**: TBD
**UI hint**: yes

## Progress

**Execution Order:**
Phases execute in numeric order: 4 → 5 → 6 → 7 → 8 → 9 → 10 → 10.1 → 11 → 12

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. API Integration | v1.0 | 2/2 | Complete | 2026-02-25 |
| 2. UI/UX and Documentation | v1.0 | 3/3 | Complete | 2026-02-25 |
| 3. Network Performance | v1.0 | 2/2 | Complete | 2026-02-25 |
| 4. Shareable URLs + AOP Discovery | v2.0 | 2/2 | Complete | 2026-02-27 |
| 5. Network Interactivity + Responsive Layout | v2.0 | 2/2 | Complete | 2026-02-28 |
| 6. Batch Analysis | v2.0 | 3/3 | Complete | 2026-03-01 |
| 7. Comparison Visualisation | v2.0 | 3/3 | Complete | 2026-03-02 |
| 8. Comparison UI Polish | v2.0 | 2/2 | Complete | 2026-03-02 |
| 9. Deployment Hardening | v2.0 | 1/1 | Complete | 2026-03-02 |
| 10. AOP Source-of-Truth Cleanup | v3.0 | 2/2 | Complete   | 2026-05-06 |
| 10.1. Demos Page Restructure (INSERTED) | v3.0 | 2/3 | In progress | - |
| 11. Gene-Set Export | v3.0 | 0/4 | Planned | - |
| 12. Upload Widget + Tech-Debt Sweep | v3.0 | 0/0 | Not started | - |
