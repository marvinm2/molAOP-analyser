# Roadmap: molAOP Analyser

## Milestones

- ✅ **v1.0 MVP** — Phases 1-3 (shipped 2026-02-27)
- 🚧 **v2.0 Batch Analysis & Dynamic AOPs** — Phases 4-7 (in progress)

## Phases

<details>
<summary>✅ v1.0 MVP (Phases 1-3) — SHIPPED 2026-02-27</summary>

- [x] Phase 1: API Integration (2/2 plans) — completed 2026-02-25
- [x] Phase 2: UI/UX and Documentation (3/3 plans) — completed 2026-02-25
- [x] Phase 3: Network Performance (2/2 plans) — completed 2026-02-25

</details>

### 🚧 v2.0 Batch Analysis & Dynamic AOPs (In Progress)

**Milestone Goal:** Enable multi-dataset comparison and dynamic AOP discovery, transforming the tool from single-dataset analysis to cross-condition investigation.

- [x] **Phase 4: Shareable URLs + AOP Discovery** - UUID-based result sharing and SPARQL-driven dynamic AOP list replace hardcoded config (completed 2026-02-27)
- [x] **Phase 5: Network Interactivity + Responsive Layout** - Expand-all toggle, KE detail panel, FDR filter slider, and responsive CSS for tablet/laptop (completed 2026-02-28)
- [x] **Phase 6: Batch Analysis** - Sequential multi-file upload with condition tagging, SQLite state persistence, and harmonised background gene sets (completed 2026-03-01)
- [ ] **Phase 7: Comparison Visualisation** - Side-by-side heatmap, comparison table, multi-condition AOP network overlay, and delta mode

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
**Plans**: TBD

Plans:
- [ ] 07-01: Comparison service + route — `comparison_service.py` pandas pivot, `GET /compare/<batch_id>` route, `compare.html` template scaffold
- [ ] 07-02: Comparison heatmap + table — Plotly heatmap (KE × condition, -log10 FDR), comparison enrichment table, condition colour legend
- [ ] 07-03: Comparison network overlay — multi-condition Cytoscape.js node styling on shared KE skeleton, single/comparison toggle, delta mode, locked dagre layout

## Progress

**Execution Order:**
Phases execute in numeric order: 4 → 5 → 6 → 7

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. API Integration | v1.0 | 2/2 | Complete | 2026-02-25 |
| 2. UI/UX and Documentation | v1.0 | 3/3 | Complete | 2026-02-25 |
| 3. Network Performance | v1.0 | 2/2 | Complete | 2026-02-25 |
| 4. Shareable URLs + AOP Discovery | v2.0 | 0/2 | Not started | - |
| 5. Network Interactivity + Responsive Layout | 2/2 | Complete   | 2026-02-28 | - |
| 6. Batch Analysis | 5/5 | Complete   | 2026-03-02 | 2026-03-01 |
| 7. Comparison Visualisation | v2.0 | 0/3 | Not started | - |
