---
phase: 07-comparison-visualisation
plan: "01"
subsystem: comparison-page
tags: [comparison, batch, heatmap, cytoscape, flask-route, pandas, vhp4safety-palette]
dependency_graph:
  requires: [06-01, 06-02, 06-03]
  provides: [comparison-data-layer, compare-page-scaffold]
  affects: [07-02, 07-03]
tech_stack:
  added: [comparison_service.py]
  patterns: [pandas-pivot-table, neg-log10-transform, VHP4Safety-palette-assignment]
key_files:
  created:
    - services/comparison_service.py
    - templates/compare.html
  modified:
    - app.py
    - templates/batch_summary.html
    - static/css/style.css
decisions:
  - "CONDITION_PALETTE colours assigned in fixed order by upload position — alphabetical sort prevented by explicit pivot column reindex"
  - "Significance threshold for -log10 transform set at FDR <= 0.05; non-significant cells return None for JSON safety"
  - "KE rows sorted by mean -log10(FDR) significance descending via mean_sig Series"
  - "condition_colors sliced from CONDITION_PALETTE[:len(conditions)] — supports up to 8 conditions"
  - "compare.html tab switching dispatches CustomEvent('modechange') for Plans 02/03 to listen to"
metrics:
  duration: "~8 minutes"
  completed: "2026-03-02"
  tasks_completed: 2
  files_created: 2
  files_modified: 3
---

# Phase 07 Plan 01: Comparison Service and Page Scaffold Summary

**One-liner:** Comparison data service with pandas FDR pivot, -log10 transform, and tabbed compare.html page scaffold using VHP4Safety palette for condition colour coding.

## What Was Built

### Task 1: Comparison Service and Flask Route (commit: ba862a2)

Created `services/comparison_service.py` with:

- `CONDITION_PALETTE` — 8 VHP4Safety hex colours in fixed order by upload position (Sky Blue excluded as too similar to Light Blue)
- `build_comparison_matrix(conditions)` — pivots ConditionRecord enrichment JSON to KE × condition matrix, applies -log10(FDR) transform for significance visualisation, sorts KEs by mean significance descending, preserves upload-position column order via explicit reindex

Added to `app.py`:
- Import of `build_comparison_matrix` and `CONDITION_PALETTE` from the new service
- `GET /batch/<uuid>/compare` route (`batch_compare`) that queries the batch, builds the matrix, serialises to JSON, and renders `compare.html`

### Task 2: compare.html, batch_summary.html, style.css (commit: 385b6cf)

Created `templates/compare.html` extending `base.html`:
- Loads Plotly.js 2.30.0, Cytoscape 3.33.1 + dagre 0.8.5 + cytoscape-dagre 2.5.0 from CDN
- Breadcrumb navigation linking back to Batch Summary
- Page title card with subtitle showing AOP label and condition count
- Top control bar: Absolute/Delta mode toggle + reference condition dropdown (hidden in Absolute mode) + condition colour legend (horizontal pill row)
- Three-tab interface (Heatmap, Table, Network) with working JS tab switching
- Plotly.Plots.resize() called on heatmap tab activation to fix render-in-hidden-div issue
- Lazy cyInitialised flag for network tab (Plan 07-03 hook)
- CustomEvent('modechange') dispatched on mode/ref change for Plans 07-02/03 to listen to
- Empty state rendering when comparisonData has no ke_labels

Modified `templates/batch_summary.html`:
- "Compare All Conditions" button (class btn btn--primary) linking to `/batch/<uuid>/compare`
- Route guard in `batch_summary` view already prevents non-complete batches from showing this page

Added to `static/css/style.css` (before responsive section):
- `.compare-subtitle`, `.compare-controls`, `.compare-controls__mode`, `.compare-controls__ref`
- `.compare-legend`, `.compare-legend__item`, `.compare-legend__swatch`
- `.compare-tabs`, `.compare-tab`, `.compare-tab--active`
- `.compare-tab-panels`, `.compare-tab-panel`, `.batch-summary-actions`
- Responsive rules inside `@media (max-width: 768px)`: stacked controls, unindented legend, smaller tabs

## Deviations from Plan

None — plan executed exactly as written.

## Self-Check

### Files Created

- [x] `services/comparison_service.py` exists
- [x] `templates/compare.html` exists
- [x] `templates/batch_summary.html` modified (Compare All button added)
- [x] `app.py` modified (import + route added)
- [x] `static/css/style.css` modified (comparison CSS added)

### Commits

- [x] ba862a2 — feat(07-01): add comparison service and GET /batch/<uuid>/compare route
- [x] 385b6cf — feat(07-01): create compare.html scaffold with tabs, control bar, legend, and Compare All button

## Self-Check: PASSED
