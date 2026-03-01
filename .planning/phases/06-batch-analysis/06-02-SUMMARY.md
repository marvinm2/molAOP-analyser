---
phase: 06-batch-analysis
plan: 02
subsystem: ui, api, database
tags: [flask, batch, enrichment, harmonisation, threading, cytoscape, sqlite, jinja2]

requires:
  - phase: 06-01
    provides: BatchRecord, ConditionRecord, batch_service.py scaffolding, cisplatin demo utilities

provides:
  - validate_batch_columns() — column consistency check across all batch files
  - harmonise_backgrounds() — intersection gene universe across all batch files
  - run_batch() — background-thread-safe orchestration of sequential per-condition enrichment
  - _run_condition() — per-file load, harmonise-filter, enrich, network build, DB persist
  - id_column, fc_column, pval_column on BatchRecord — shared column mapping
  - GET /batch — batch upload wizard page
  - POST /batch/upload — multi-file receive, save to UUID dir, column detection, JSON preview
  - POST /batch/analyze — validate, harmonise, create DB records, launch background thread
  - GET /batch/<uuid>/status — polling endpoint for htmx progress modal
  - POST /batch/<uuid>/cancel — cooperative cancel flag on BatchRecord
  - templates/batch.html — three-step wizard (upload, tag conditions, settings)

affects:
  - 06-03 (batch results page — reads ConditionRecord.enrichment_json, network_json via /batch/<uuid>/status and /batch/<uuid>/results)

tech-stack:
  added:
    - threading.Thread (daemon background thread for batch analysis)
  patterns:
    - run_batch() accepts db_url string and creates its own SQLAlchemy session — thread-safe DB access
    - Background thread launched immediately; client polls /batch/<uuid>/status
    - harmonise_backgrounds() computes intersection of all per-file gene universes before any enrichment
    - _run_condition() filters df_processed to harmonised_genes set before passing to run_enrichment_analysis()

key-files:
  created:
    - templates/batch.html
  modified:
    - services/batch_service.py
    - database.py
    - app.py

key-decisions:
  - "run_batch() accepts db_url string instead of a session object — SQLAlchemy sessions are not thread-safe, so the background thread creates its own session via sessionmaker(bind=create_engine(db_url))"
  - "batch/upload and batch/analyze are CSRF-exempt — JS fetch() POSTs don't carry CSRF cookies in the same way as form submissions"
  - "Stop-batch strategy on first condition failure — sets batch.status='failed' and returns immediately rather than continuing with remaining conditions"
  - "Batch upload dir cleanup in run_batch() finally block — files removed regardless of success/failure after analysis"
  - "parse_cisplatin_filename passed as template global to batch.html via render_template kwarg so cisplatin demo checkboxes get pre-populated metadata"

patterns-established:
  - "Background thread pattern: launch daemon thread from route, return JSON immediately, client polls status endpoint"
  - "Harmonise-then-filter: compute intersection outside the enrichment loop, filter each df_processed inside _run_condition"

requirements-completed: [BATC-01, BATC-02, BATC-05]

duration: 5min
completed: 2026-03-01
---

# Phase 06 Plan 02: Batch Upload Wizard and Orchestration Summary

**Three-step batch upload wizard (templates/batch.html) with harmonise-then-enrich orchestration (batch_service.py) and five new Flask routes including background-threaded /batch/analyze**

## Performance

- **Duration:** 5 min
- **Started:** 2026-03-01T18:46:19Z
- **Completed:** 2026-03-01T18:52:05Z
- **Tasks:** 2
- **Files modified:** 3 (services/batch_service.py, database.py, app.py) + 1 created (templates/batch.html)

## Accomplishments
- Implemented validate_batch_columns(), harmonise_backgrounds(), run_batch(), and _run_condition() in batch_service.py — the complete batch orchestration pipeline
- Added id_column, fc_column, pval_column columns to BatchRecord for shared column mapping across all files in a batch
- Created a full three-step batch upload wizard in batch.html (file upload with drag-and-drop, condition tagging with auto-population for cisplatin demos, AOP selection with typeahead, threshold controls)
- Added five Flask routes covering the full batch lifecycle (upload, analyze, status polling, cancel)
- Batch analysis runs in a daemon thread so /batch/analyze returns immediately and the client can poll /batch/<uuid>/status for progress

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement batch orchestration functions in batch_service.py** - `8597b9a` (feat)
2. **Task 2: Create batch.html template and batch Flask routes** - `8c964c0` (feat)

**Plan metadata:** (created in next commit)

## Files Created/Modified
- `services/batch_service.py` — Added validate_batch_columns, harmonise_backgrounds, run_batch, _run_condition; added json/datetime/pandas imports
- `database.py` — Added id_column, fc_column, pval_column String(100) columns to BatchRecord
- `app.py` — Added threading import, batch service imports, five batch routes
- `templates/batch.html` — Three-step batch wizard with JS step navigation, drag-and-drop, demo checkboxes, AOP typeahead, progress modal, polling

## Decisions Made
- run_batch() accepts db_url string instead of a session object so it is thread-safe (SQLAlchemy sessions cannot be shared across threads)
- /batch/upload and /batch/analyze use @csrf.exempt since JS fetch() POSTs require explicit CSRF handling and these are API-style endpoints
- Stop-batch strategy: on first condition failure, batch.status='failed' and the loop exits (simpler than partial-success logic)
- Uploaded files are cleaned up in run_batch() finally block regardless of outcome
- parse_cisplatin_filename passed as a template kwarg so the Jinja2 template can call it to pre-populate demo checkbox metadata

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Added /batch/<uuid>/status and /batch/<uuid>/cancel routes**
- **Found during:** Task 2 (batch Flask routes)
- **Issue:** The plan specifies the batch runs in a background thread and returns immediately, but the polling endpoint (/batch/<uuid>/status) needed by the JS progress modal was not listed in the three routes specification — it would have left the wizard with no way to track progress
- **Fix:** Added GET /batch/<uuid>/status and POST /batch/<uuid>/cancel routes alongside the three specified routes
- **Files modified:** app.py
- **Verification:** Route list confirmed via Flask url_map inspection
- **Committed in:** 8c964c0 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 missing critical)
**Impact on plan:** The status and cancel routes are required for the progress modal JavaScript to function. The plan's success criteria include "Batch runs in background thread, updating ConditionRecord.status per file" — the polling route is what exposes this status to the client.

## Issues Encountered
- Pre-existing test failure in `tests/test_database.py::TestDatabaseOperations::test_experiment_record_model` (upload_timestamp is None before DB insert — SQLAlchemy column default only fires on INSERT, not Python object instantiation). This failure existed before this plan and is out of scope.

## Next Phase Readiness
- All batch orchestration infrastructure is complete for plan 06-03 (batch results page)
- ConditionRecord.enrichment_json, network_json, ke_gene_json, ke_type_map_json, ke_title_map_json are all populated by _run_condition() and readable via the status/results endpoints
- The polling pattern (JS setInterval + /batch/<uuid>/status) is wired in batch.html; 06-03 can extend with the /batch/<uuid>/results view

---
*Phase: 06-batch-analysis*
*Completed: 2026-03-01*
