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

- [ ] **NETW-01**: User can expand/collapse all gene nodes with one toggle button
- [ ] **NETW-02**: Expand-all toggle disabled when total gene count exceeds 200
- [ ] **NETW-03**: User can view KE detail panel on node click (title, type, FDR, gene list, AOP-Wiki link)
- [ ] **NETW-04**: User can filter network by FDR significance threshold via slider

### Upload & UX

- [ ] **UPUX-01**: Responsive layout at 768px and 1024px breakpoints
- [ ] **UPUX-02**: Drag-and-drop file upload hint on upload page
- [ ] **UPUX-03**: Network and results layout stacks vertically on tablet viewports

### Batch Analysis

- [ ] **BATC-01**: User can upload multiple files (2-10) in a sequential queue
- [ ] **BATC-02**: User can tag each file with condition label (dose, timepoint, treatment)
- [ ] **BATC-03**: Per-file progress indicator during batch analysis
- [ ] **BATC-04**: Batch state persisted in SQLite (not session cookies)
- [ ] **BATC-05**: Background gene set harmonisation across conditions for valid comparison

### Comparison Visualisation

- [ ] **COMP-01**: KE × condition heatmap showing -log10(FDR) significance
- [ ] **COMP-02**: Comparison enrichment table (rows = KEs, columns = conditions)
- [ ] **COMP-03**: Multi-condition AOP network overlay on shared KE skeleton
- [ ] **COMP-04**: User can toggle between single-condition and comparison view
- [ ] **COMP-05**: Delta mode showing change relative to user-selected reference condition
- [ ] **COMP-06**: Condition colour legend using VHP4Safety palette

## Future Requirements

Deferred to future release. Tracked but not in current roadmap.

### Deferred Interactivity

- **NETW-05**: Animated time-course playback across ordered conditions
- **NETW-06**: KE description tooltips from AOP-Wiki SPARQL on hover

### Deferred Export

- **EXPO-01**: Batch comparison CSV export (full heatmap matrix)
- **EXPO-02**: Network SVG export for publication use

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
| NETW-01 | Phase 5 | Pending |
| NETW-02 | Phase 5 | Pending |
| NETW-03 | Phase 5 | Pending |
| NETW-04 | Phase 5 | Pending |
| UPUX-01 | Phase 5 | Pending |
| UPUX-02 | Phase 5 | Pending |
| UPUX-03 | Phase 5 | Pending |
| BATC-01 | Phase 6 | Pending |
| BATC-02 | Phase 6 | Pending |
| BATC-03 | Phase 6 | Pending |
| BATC-04 | Phase 6 | Pending |
| BATC-05 | Phase 6 | Pending |
| COMP-01 | Phase 7 | Pending |
| COMP-02 | Phase 7 | Pending |
| COMP-03 | Phase 7 | Pending |
| COMP-04 | Phase 7 | Pending |
| COMP-05 | Phase 7 | Pending |
| COMP-06 | Phase 7 | Pending |

**Coverage:**
- v2.0 requirements: 27 total
- Mapped to phases: 27
- Unmapped: 0 ✓

---
*Requirements defined: 2026-02-27*
*Last updated: 2026-02-27 — traceability completed after roadmap creation*
