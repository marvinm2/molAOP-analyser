---
phase: 03-network-performance
plan: 01
subsystem: ui
tags: [cytoscape, javascript, network-visualisation, dagre, performance]

# Dependency graph
requires:
  - phase: 02-ui-ux-and-documentation
    provides: results.html extended from base.html with metadata card and enrichment table
provides:
  - Cytoscape.js 3.33.1 with LR dagre layout and loading spinner on results.html
  - Per-KE click-to-expand gene nodes with shared gene deduplication
  - Collapse All button to restore KE-only dagre LR view
  - Spinner overlay for initial render and KE expansion
affects: [03-02-network-performance]

# Tech tracking
tech-stack:
  added: [cytoscape@3.33.1]
  patterns: [progressive-disclosure gene expansion, cy.batch() for bulk DOM updates, requestAnimationFrame before heavy layout work, promiseOn layoutstop for async spinner hide]

key-files:
  created: []
  modified: [templates/results.html]

key-decisions:
  - "Cytoscape.js upgraded from 3.26.0 to 3.33.1 for performance improvements"
  - "rankDir changed from TB to LR — left-to-right dagre better suits AOP pathway flow"
  - "Per-KE click-to-expand replaces bulk addGeneNodes() — avoids browser freeze on large networks (AOP:2 ~696 genes)"
  - "requestAnimationFrame used before cy.layout() call so spinner renders before heavy work begins"
  - "Shared gene nodes deduped by cy.getElementById check before cy.add — genes in multiple KEs show as one node with edges to all KEs"
  - "collapseAll() removes .gene nodes and restores dagre LR — no per-KE toggle, just full reset"

patterns-established:
  - "Spinner pattern: showSpinner() before heavy op, hideSpinner() on layoutstop or catch block"
  - "cy.batch() wraps all node/edge additions, layout() runs outside batch for correct positioning"
  - "Dedup by ID: cy.getElementById(id).length === 0 before cy.add() for both nodes and edges"

requirements-completed: [NET-01]

# Metrics
duration: 2min
completed: 2026-02-25
---

# Phase 3 Plan 01: Network Performance — Cytoscape Overhaul Summary

**Cytoscape.js 3.33.1 with LR dagre layout, spinner overlay, per-KE click-to-expand with shared gene dedup, and Collapse All button — replacing the old bulk addGeneNodes/toggleGenes pattern**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-25T18:45:06Z
- **Completed:** 2026-02-25T18:47:19Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Upgraded Cytoscape.js from 3.26.0 to 3.33.1
- Initial KE-only network renders in left-to-right dagre layout with a loading spinner that hides on layoutstop
- Clicking a KE node progressively expands its genes using cy.batch() with shared gene dedup; "Expanding KE..." spinner shown and hidden asynchronously via promiseOn layoutstop
- Large network delayed message ("Large network — this may take a moment") appears after 3 seconds of spinning
- "Collapse All" button removes all gene nodes and restores KE-only dagre LR layout
- report generation updated to expand all KEs programmatically before PNG capture

## Task Commits

Each task was committed atomically:

1. **Task 1: Add spinner overlay, upgrade Cytoscape.js, change layout to LR, implement per-KE expand and collapse all** - `34715f2` (feat)

## Files Created/Modified
- `templates/results.html` - Cytoscape.js CDN upgrade, spinner CSS/HTML, LR dagre layout, expandedKEs Set, tap handler with cy.batch(), collapseAll(), updated handleReportSubmit, removed addGeneNodes/toggleGenes/genesVisible

## Decisions Made
- Upgraded Cytoscape.js from 3.26.0 to 3.33.1 for latest performance and compatibility
- Changed dagre rankDir from TB to LR — horizontal flow more natural for AOP pathway diagrams where causality flows left to right
- Used requestAnimationFrame before layout call so the browser can paint the spinner before CPU-intensive layout work begins
- Per-KE expand chosen over bulk-add to avoid freezing the browser on large networks; user only pays the cost of expanding KEs they care about
- Shared gene nodes: cy.getElementById check prevents duplicate nodes; single gene node with edges to all KEs it belongs to shows cross-KE biological overlap

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- `git stash pop` conflicted on pycache binary files — resolved by checking out the pycache files from HEAD before re-applying stash. No impact on deliverable.
- 4 pre-existing test failures in test_flask_routes.py confirmed to exist before any changes were made (verified by stashing changes and running tests against original code).

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Per-KE expand and collapse all fully functional in templates/results.html
- Ready for Phase 03-02 (zoom controls overlay and PNG export improvements)
- No blockers

---
*Phase: 03-network-performance*
*Completed: 2026-02-25*
