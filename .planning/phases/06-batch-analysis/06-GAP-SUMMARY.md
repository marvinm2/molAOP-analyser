---
phase: 06-batch-analysis
plan: gap
subsystem: ui
tags: [htmx, batch, templates, javascript]

# Dependency graph
requires:
  - phase: 06-batch-analysis
    provides: batch analysis wizard, progress polling routes, htmx partials
provides:
  - htmx 2.0.4 loaded globally in base.html before block scripts
  - file card metadata showing row count and auto-detected column names after upload
affects: [07-comparison]

# Tech tracking
tech-stack:
  added: [htmx 2.0.4 (CDN with SRI hash)]
  patterns: [global CDN script before block scripts in base.html, progressive enhancement of file card metadata]

key-files:
  created: []
  modified:
    - templates/base.html
    - templates/batch.html

key-decisions:
  - "htmx loaded from unpkg CDN with SRI integrity hash for subresource integrity"
  - "htmx script placed before {% block scripts %} in base.html so all child templates inherit it"
  - "file-card__meta uses em-dash separator (U+2014) between parts, consistent with tagging step display"
  - "rowCount null-check uses != null (not !== null) to cover both null and undefined safely"

patterns-established:
  - "Global third-party libraries belong in base.html before {% block scripts %}, not in individual templates"
  - "File card meta built from parts array joined with em-dash — avoids conditional string concatenation"

requirements-completed: [BATC-01, BATC-02, BATC-03]

# Metrics
duration: 1min
completed: 2026-03-02
---

# Phase 06 Plan GAP: UAT Gap Closure Summary

**htmx 2.0.4 loaded globally in base.html and file cards augmented with row count and detected column names after upload**

## Performance

- **Duration:** 1 min
- **Started:** 2026-03-02T10:34:53Z
- **Completed:** 2026-03-02T10:35:43Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- htmx 2.0.4 script tag with SRI hash added to base.html, making htmx globally available on all pages
- Batch wizard progress modal, cancel polling, and HX-Redirect auto-navigation now function without ReferenceError
- File cards in the batch wizard now display row count and auto-detected column names (ID, FC, p-value) after the upload step completes

## Task Commits

Each task was committed atomically:

1. **Task 1: Load htmx CDN in base.html** - `643e812` (feat)
2. **Task 2: Show row count and detected columns in file cards** - `9e2c7a9` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified
- `templates/base.html` - Added htmx 2.0.4 CDN script tag with SRI integrity hash immediately before {% block scripts %}
- `templates/batch.html` - Updated renderFileList() to show rowCount and suggestions.id_col/fc_col/pval_col in file-card__meta

## Decisions Made
- htmx loaded from unpkg CDN with SRI integrity hash — same CDN approach used by the rest of the project for external libraries
- htmx script placed before `{% block scripts %}` in base.html so child templates get the global object without needing per-template script tags
- File card meta uses a `parts` array joined with em-dash (U+2014), matching the separator style already used in the tagging step display (line 889 of batch.html)
- `rowCount != null` used (not `!== null`) to safely handle both `null` and `undefined` from the upload response

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- All 5 UAT gaps (Tests 2, 6, 7, 8, 11) from Phase 06 verification are closed
- htmx global availability enables all htmx-dependent batch functionality (polling, cancel, redirect)
- File card metadata gives users immediate feedback on detected columns after upload
- Phase 07 (Comparison) can proceed; no outstanding blockers from Phase 06

## Self-Check

- templates/base.html contains htmx script tag: FOUND
- templates/batch.html contains entry.rowCount check: FOUND
- Commit 643e812 exists: FOUND
- Commit 9e2c7a9 exists: FOUND

## Self-Check: PASSED

---
*Phase: 06-batch-analysis*
*Completed: 2026-03-02*
