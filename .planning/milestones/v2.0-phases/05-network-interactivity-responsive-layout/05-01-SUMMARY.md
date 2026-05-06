---
phase: 05-network-interactivity-responsive-layout
plan: 01
subsystem: ui
tags: [cytoscape, javascript, css, flask, jinja2, drawer, fdr-filter, expand-all]

# Dependency graph
requires:
  - phase: 04-shareable-urls-aop-discovery
    provides: results.html template with keToGenes, network_json, ke_type_map, ke_title_map, table_json template variables
provides:
  - Expand All button with 200-gene guard that expands all KEs at once
  - KE detail drawer (right-side panel / mobile bottom sheet) on node click
  - FDR significance filter slider that grays out non-significant KE nodes
  - network-controls-bar above the Cytoscape canvas with all controls consolidated
affects: [06-batch-analysis, results-ui, network-visualization]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Cytoscape node class toggling for visual filtering (fdr-filtered) without layout disruption
    - BEM CSS naming for drawer sub-components (ke-drawer__*)
    - CSS transition on right/bottom properties for slide-in drawer vs bottom sheet responsive pattern
    - Global JS data lookups extracted from template variables once (fdrByKe, keTitleMap, keTypeMap)

key-files:
  created: []
  modified:
    - templates/results.html
    - static/css/style.css

key-decisions:
  - "FDR filter uses addClass/removeClass('fdr-filtered') with opacity:0.25 — nodes stay in layout, only become semi-transparent"
  - "Drawer uses CSS transition on 'right' property (desktop) and 'bottom' (mobile <768px) — no JS animation"
  - "expandAll() re-uses same batch/dedup pattern as the per-KE tap handler to avoid duplicates"
  - "tableData/keTitleMap/keTypeMap hoisted to global script scope so copyLinkBtn handler and drawer can share them without re-parsing template JSON"

patterns-established:
  - "Cytoscape filtering: always use addClass/removeClass, never node.hide() or node.remove() for non-destructive filtering"
  - "Drawer open state driven by CSS class (ke-drawer--open) added/removed via JS — no inline style manipulation"
  - "Network controls consolidated in .network-controls-bar above #cy, not below"

requirements-completed: [NETW-01, NETW-02, NETW-03, NETW-04]

# Metrics
duration: 5min
completed: 2026-02-28
---

# Phase 5 Plan 1: Network Interactivity — Expand-All, KE Detail Drawer, and FDR Filter Summary

**Cytoscape.js interactivity upgraded: expand-all with 200-gene guard, click-to-open KE detail drawer with gene list and AOP-Wiki link, and FDR slider that grays non-significant KE nodes without disrupting layout**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-28T14:46:35Z
- **Completed:** 2026-02-28T14:51:42Z
- **Tasks:** 1 of 1
- **Files modified:** 2

## Accomplishments
- Added `.network-controls-bar` above the Cytoscape container consolidating Expand All, Collapse All, FDR slider, Apply, Reset View, Export PNG, Download Network
- Implemented `expandAll()` function with 200-gene guard: disabled with "Too many genes (N)" label when total gene count exceeds limit
- Added KE detail drawer that slides in from right on KE node click, showing KE type badge, FDR value, title, AOP-Wiki link, and sorted gene list (significant first, then by |log2FC|)
- Drawer becomes a bottom sheet on screens <768px via CSS media query
- FDR Apply button toggles `fdr-filtered` class (opacity 0.25, gray background) on KE nodes above threshold — layout preserved
- Extended existing tap-node handler: already-expanded KEs open drawer immediately; new expansions open drawer after layoutstop

## Task Commits

Each task was committed atomically:

1. **Task 1: Add controls bar, expand-all toggle, KE detail drawer, FDR filter** - `942a276` (feat)

**Plan metadata:** _(to be added in final commit)_

## Files Created/Modified
- `/home/marvin/Documents/Services/molAOP-analyser/molAOP-analyser/templates/results.html` - Added controls bar HTML, drawer HTML, global JS data lookups, expandAll/openDrawer/closeDrawer/buildDrawerHTML/FDR handler functions, extended tap-node handler, added fdr-filtered Cytoscape style
- `/home/marvin/Documents/Services/molAOP-analyser/molAOP-analyser/static/css/style.css` - Added .network-controls-bar, .fdr-slider-group, .ke-drawer and all BEM sub-component styles, bottom-sheet media query

## Decisions Made
- FDR filter uses `addClass('fdr-filtered')` with `opacity: 0.25` — preserves node positions in DAG layout per plan anti-pattern guidance
- Drawer transitions animate `right` property (desktop) and `bottom` property (mobile) — no JS animation, no height/width animation
- Global JS variables `tableData`, `keTitleMap`, `keTypeMap` hoisted above all handlers so copyLinkBtn handler can reuse them without re-parsing Jinja JSON

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- `git stash pop` failed due to binary pycache file conflicts, which reverted my working tree changes. Re-applied all changes manually from scratch. No impact on final result.
- Pre-existing test failures in `test_flask_routes.py` (4 tests): `test_index_route`, `test_analyze_route_success`, `test_generate_report_route_html`, `test_full_workflow_integration` — all failing before my changes (confirmed by stash test). These are unrelated to results.html template changes (mock data returns empty list causing `table[0].keys()` UndefinedError in Jinja2).

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Network interactivity complete — Expand All, KE detail drawer, and FDR filter all implemented
- Phase 5 Plan 2 can proceed to tackle responsive layout if planned

---
*Phase: 05-network-interactivity-responsive-layout*
*Completed: 2026-02-28*
