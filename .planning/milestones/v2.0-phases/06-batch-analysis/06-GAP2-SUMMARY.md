---
phase: 06-batch-analysis
plan: GAP2
subsystem: ui
tags: [flask, jinja2, javascript, pandas, batch-upload, file-preview]

# Dependency graph
requires:
  - phase: 06-batch-analysis
    provides: batch_upload() route with file_previews JSON response and file card rendering in batch.html
provides:
  - True file row count via byte-level line iteration (not capped at 200)
  - head_rows field in batch upload JSON response (first 5 data rows as list of dicts)
  - .file-card__head CSS and DOM rendering of 5-row preview table in file cards after upload
affects: [07-comparison]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Byte-level file line count: open(dest, 'rb') + sum(1 for _ in _f) - 1 for fast true row count without full CSV parse"
    - "head_rows pattern: df_head.head(5).to_dict(orient='records') returns list of column->value dicts safe for JSON"

key-files:
  created: []
  modified:
    - app.py
    - templates/batch.html

key-decisions:
  - "byte-level line count (open rb + sum) keeps nrows=200 parse for column detection performance, adds only one small file scan for count"
  - "head_rows returned as orient='records' — list of dicts lets JS render table with Object.keys/Object.values without column ordering issues"
  - "head table conditionally rendered only when entry.headRows is populated — no empty table before upload"

patterns-established:
  - "File card data flow: server adds head_rows to JSON response, JS stores on fileQueue entry, renderFileList() reads entry to build DOM"

requirements-completed: [BATC-01]

# Metrics
duration: 3min
completed: 2026-03-02
---

# Phase 06 GAP2: Batch File Card Metadata Fix Summary

**True row count via byte-level line scan and 5-row head preview table in batch upload file cards**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-02T11:10:57Z
- **Completed:** 2026-03-02T11:14:06Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Fixed row_count cap: was always <= 200 (len of nrows=200 DataFrame), now returns true file row count via byte-level line iteration (42,545 for GSE90122_SR12813.tsv)
- Added head_rows to batch upload JSON response: first 5 data rows as list of dicts, using df_head.head(5).to_dict(orient='records')
- File cards now render a scrollable .file-card__head table after upload step with correct column headers and 5 data rows for column verification

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix row count and add head_rows to batch_upload() response** - `5f595c8` (fix)
2. **Task 2: Render head preview table in batch.html file cards** - `18dbb68` (feat)

**Plan metadata:** (see final docs commit)

## Files Created/Modified
- `app.py` - Fixed row_count calculation (both uploaded-file and demo-file blocks); added head_rows to file_previews dicts; updated docstring
- `templates/batch.html` - Added .file-card__head CSS; stored entry.headRows from upload response; rendered 5-row preview table per file card in renderFileList()

## Decisions Made
- Kept nrows=200 in pd.read_csv — still correct for column detection performance; byte-level scan adds negligible overhead for true count
- Used orient='records' for head_rows — produces list of dicts, enabling JS Object.keys/Object.values table rendering without manual column index management
- Head table is conditionally rendered (entry.headRows && entry.headRows.length > 0) so pre-upload file cards are unchanged

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- GAP2 closes the last UAT issue from Test 4 (File Card Metadata After Upload)
- Phase 06 batch analysis fully complete; Phase 07 (comparison) can proceed
- No blockers introduced

## Self-Check: PASSED

- app.py: FOUND
- templates/batch.html: FOUND
- 06-GAP2-SUMMARY.md: FOUND
- Commit 5f595c8 (Task 1): FOUND
- Commit 18dbb68 (Task 2): FOUND

---
*Phase: 06-batch-analysis*
*Completed: 2026-03-02*
