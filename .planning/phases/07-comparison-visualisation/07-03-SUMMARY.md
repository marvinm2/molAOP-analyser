---
phase: 07-comparison-visualisation
plan: "03"
subsystem: comparison-page
tags: [cytoscape, pie-chart, comparison-network, delta-mode, ke-detail-drawer, dagre, vhp4safety-palette]
dependency_graph:
  requires:
    - phase: 07-01
      provides: compare.html scaffold, comparisonData JSON, modechange CustomEvent, tab switching, cyInitialised flag
    - phase: 07-02
      provides: computeDeltaMatrix(), hexToRgba(), heatmap and table rendering
  provides:
    - Cytoscape.js comparison network with pie-chart KE nodes (COMP-03)
    - Single-condition / comparison view toggle (COMP-04)
    - Delta mode network update without dagre re-layout (COMP-05)
    - KE detail drawer with per-condition horizontal bar chart
  affects: []
tech-stack:
  added: []
  patterns:
    - Cytoscape pie-N-background-color/-size/-opacity for multi-condition node overlay
    - Dagre LR layout computed once on first tab activation, positions saved in savedPositions{}
    - cy.batch() for efficient bulk node data update without layout re-run
    - Lazy Cytoscape initialisation inside tab-switch handler to avoid 0x0 container
    - buildComparisonElements() / buildComparisonStylesheet() functional builders passed to cy instance
    - Fixed-position drawer using CSS transition on right property (desktop) / bottom (mobile)
    - Plotly.newPlot inside KE drawer for per-condition horizontal bar chart
key-files:
  created: []
  modified:
    - app.py
    - templates/compare.html
key-decisions:
  - "Pie slice sizes normalised as equal shares: each condition gets 100/N percent of the full circle; slices sized by min(val,10)/10 * sliceMax"
  - "Delta mode clamps negative deltas to 0 for slice sizing (less significant than reference = no slice shown)"
  - "buildComparisonElements() skips gene-level nodes — only KE nodes (those with ke_type defined or present in ke_labels) included in comparison network"
  - "Single-condition view disables pie-size entirely (pie-size:0%), uses solid background-color per-node selectors from enrichment data"
  - "openKeDrawer() reads live currentMode and ref-condition at click time — avoids stale closure over init-time values"
  - "Drawer re-render triggered in updateNetwork() if drawer is open when mode changes"
  - "KE detail bar chart threshold line at x=1.301 (FDR=0.05) in absolute mode; x=0 (reference line) in delta mode"
requirements-completed: [COMP-03, COMP-04, COMP-05]
duration: "~3 minutes"
completed: "2026-03-02"
tasks_completed: 2
files_created: 0
files_modified: 2
---

# Phase 07 Plan 03: Comparison Network with Pie-Chart Nodes Summary

**Cytoscape.js comparison network with VHP4Safety-coloured pie-chart KE nodes, dagre LR layout locked after first render, single/comparison toggle, KE detail drawer with per-condition bar chart, and full delta mode integration across all three views.**

## Performance

- **Duration:** ~3 minutes
- **Started:** 2026-03-02T12:34:15Z
- **Completed:** 2026-03-02T12:37:51Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

### Task 1: Comparison network with pie-chart nodes and locked dagre layout (commit: 4bbc12c)

Modified `app.py`:
- `batch_compare` route now extracts the first complete condition's `network_json` and passes it as `network_json` template variable (provides the KE-level AOP skeleton for the comparison network)

Modified `templates/compare.html`:
- Added CSS for `.compare-network-controls` (view toggle bar) and `.compare-ke-drawer` (sliding detail panel, responsive: right-slide on desktop, bottom-slide on mobile)
- Added HTML: network controls bar with Comparison/Single Condition toggle buttons and condition dropdown; KE detail drawer div with close button, title, Plotly div, and AOP-Wiki link
- `buildComparisonElements(mode, refIndex)`: builds Cytoscape elements array with `pie_size_1`..`pie_size_N` data attrs; skips gene-level nodes; truncates labels to 25 chars + '...'
- `buildComparisonStylesheet(conditionCount)`: base node style with `pie-size: 100%`, grey background for non-significant portions; per-condition `pie-N-background-color/-size/-opacity`; KER edge style (thin grey arrows)
- `buildSingleConditionStylesheet(condIdx)`: disables pie (`pie-size: 0%`); adds per-node colour selectors for significant KEs in the chosen condition
- `initComparisonNetwork()`: guards with `cyInitialised` flag; creates Cytoscape instance on `#cy-compare` with dagre LR layout; saves node positions on `layoutstop`; wires tap handlers and view toggle buttons
- `updateNetwork(mode, refIndex)`: updates pie data in `cy.batch()` without re-running dagre; restores saved positions; refreshes open drawer if mode changes
- Tab-switch handler updated to call `initComparisonNetwork()` on first network tab activation
- `modechange` listener updated to call `updateNetwork()` when network is initialised

### Task 2: KE detail panel with per-condition bar chart and delta mode (commit: 4bbc12c)

Implemented as part of the same commit:
- `openKeDrawer(keId, mode, refIndex)`: renders per-condition horizontal bar chart with Plotly.newPlot; significance threshold line at x=1.301 (absolute mode) or x=0 (delta mode); AOP-Wiki link built from KE ID numeric part
- `closeKeDrawer()`: removes `.compare-ke-drawer--open` class
- Node tap handler reads live `currentMode` and `ref-condition.value` at click time to avoid stale closures
- Background tap closes drawer; close button closes drawer
- `updateNetwork()` calls `openKeDrawer()` to refresh chart if drawer is open during mode change

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Gene node exclusion using ke_type presence check**
- **Found during:** Task 1 implementation
- **Issue:** The plan described skipping gene-level nodes but didn't specify the exact condition. The network_json from `build_cytoscape_network()` mixes KE nodes and gene nodes. Gene nodes don't have `ke_type` in their data, but some KE nodes from the AOP skeleton might also lack it if they weren't in the enrichment results.
- **Fix:** Combined condition: `if (nodeData.ke_type === undefined && !keIndexMap.hasOwnProperty(keId)) return;` — skips nodes that have neither ke_type nor any enrichment data. KE nodes present in the enrichment results are always included regardless of ke_type.
- **Files modified:** templates/compare.html
- **Committed in:** 4bbc12c

**2. [Rule 1 - Bug] Stale mode closure in cy.on('tap') handler**
- **Found during:** Task 2 implementation
- **Issue:** Plan's pseudocode passed `mode` and `refIndex` from the init-time scope into the tap handler, meaning mode switches after init would still pass the initial 'absolute'/0 values to `openKeDrawer`.
- **Fix:** Tap handler reads live `currentMode` and `document.getElementById('ref-condition').value` at click time, then computes `liveRefIndex` from that.
- **Files modified:** templates/compare.html
- **Committed in:** 4bbc12c

---

**Total deviations:** 2 auto-fixed (Rule 1 — bug fixes; no scope creep)
**Impact on plan:** Both fixes essential for correctness; implementation otherwise follows plan exactly.

## COMP Requirements Completed

All 5 COMP requirements across Plans 01-03:
- **COMP-01**: Plotly heatmap with absolute/delta mode (Plan 02)
- **COMP-02**: Sortable comparison table with FDR values (Plan 02)
- **COMP-03**: Cytoscape comparison network with pie-chart KE nodes (this plan)
- **COMP-04**: Single-condition / comparison view toggle on network (this plan)
- **COMP-05**: Delta mode across all three views: heatmap, table, and network (this plan)

## Self-Check

### Files Modified

- [x] `app.py` modified — `batch_compare` route passes `network_json` to template
- [x] `templates/compare.html` modified — network rendering, drawer, pie-chart nodes

### Commits

- [x] 4bbc12c — feat(07-03): add Cytoscape comparison network with pie-chart KE nodes

## Self-Check: PASSED
