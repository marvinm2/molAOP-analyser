---
phase: 09-deployment-hardening
plan: 01
subsystem: infra
tags: [flask, gunicorn, wsgi, navigation, database, sqlalchemy]

# Dependency graph
requires:
  - phase: 06-batch-analysis
    provides: batch_page route and BatchRecord/ConditionRecord DB tables
  - phase: 08-comparison-ui-polish
    provides: base.html global nav template
provides:
  - Batch Analysis link in global navigation accessible from every page
  - Module-level init_database() call — DB initialises on any WSGI import
affects: [deployment, all-templates, production-wsgi]

# Tech tracking
tech-stack:
  added: []
  patterns: [module-level WSGI initialisation pattern — DB setup at import time not __main__]

key-files:
  created: []
  modified:
    - templates/base.html
    - app.py

key-decisions:
  - "init_database() called at module level after csrf setup — runs under gunicorn/waitress/python app.py without double-initialisation"
  - "Batch Analysis link inserted between Home and Documentation using url_for('batch_page') — no hardcoded /batch path"
  - "Config.validate_data_files() remains in __main__ only per 09-RESEARCH pitfall 3 — not promoted to module level"

patterns-established:
  - "WSGI-safe initialisation: any setup needed before route handlers should run at module level, not inside __main__"

requirements-completed: [BATC-01, SHAR-01, BATC-04]

# Metrics
duration: 2min
completed: 2026-03-02
---

# Phase 9 Plan 1: Deployment Hardening — Nav and DB Init Summary

**Batch Analysis added to global nav via url_for('batch_page'), and init_database() promoted to module level so SQLite initialises under any WSGI server without the __main__ guard.**

## Performance

- **Duration:** 2 min
- **Started:** 2026-03-02T18:47:30Z
- **Completed:** 2026-03-02T18:48:57Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Global nav bar now contains a Batch Analysis link between Home and Documentation on every page
- Database initialisation is WSGI-safe — gunicorn/waitress importing the app module triggers init_database() automatically
- __main__ block simplified — no database calls, only data file validation and dev server start
- python app.py still works identically (no regression)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add Batch Analysis link to global navigation** - `f40c624` (feat)
2. **Task 2: Move init_database() to module-level startup** - `0433da7` (feat)

**Plan metadata:** (docs commit to follow)

## Files Created/Modified

- `templates/base.html` - Added Batch Analysis li element between Home and Documentation in site-nav__links
- `app.py` - Added module-level init_database() call after csrf setup; removed init_database() from __main__ block

## Decisions Made

- init_database() placed after `csrf = CSRFProtect(app)` on line 60 — ensures DB is ready before any route handler or context processor runs
- Double-initialisation avoided by removing the call from __main__ entirely
- Config.validate_data_files() kept in __main__ per plan guidance (data file paths are only guaranteed correct when running as the main process)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Both integration gaps from the v2.0 audit are closed
- App is deployment-safe: gunicorn `gunicorn app:app` will initialise the database on first worker import
- Global nav Batch Analysis link is live for all templates extending base.html

## Self-Check: PASSED

- templates/base.html: FOUND
- app.py: FOUND
- 09-01-SUMMARY.md: FOUND
- Commit f40c624: FOUND
- Commit 0433da7: FOUND

---
*Phase: 09-deployment-hardening*
*Completed: 2026-03-02*
