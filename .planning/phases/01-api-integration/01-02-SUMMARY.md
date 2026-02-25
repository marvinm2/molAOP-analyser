---
phase: 01-api-integration
plan: 02
subsystem: api
tags: [diskcache, fanout-cache, builder-api, csv-fallback, data-source-indicator]

# Dependency graph
requires:
  - phase: 01-01
    provides: Builder API client (fetch_reference_sets_from_api) and Config attributes (BUILDER_API_URL, CACHE_DIR, CACHE_TTL)
provides:
  - diskcache.FanoutCache-backed reference set cache shared across Gunicorn workers
  - API-first load_cached_reference_sets() with transparent CSV fallback
  - data_source tuple return passed through to results template
  - Data source indicator badge on results page (three states: live/cached/local)
  - BUILDER_API_URL placeholder in docker-compose.yml

affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    - API-first loading with silent CSV fallback on any exception
    - Function returning (data, source_label) tuple to propagate provenance to UI
    - FanoutCache at module level — each worker gets its own Python object pointing to same directory

key-files:
  created: []
  modified:
    - app.py
    - cache_manager.py
    - templates/results.html
    - docker-compose.yml

key-decisions:
  - "load_cached_reference_sets() returns (reference_sets, data_source) tuple so the analyze route can pass provenance to the template without a separate lookup"
  - "diskcache.FanoutCache created at module level in app.py using Config.CACHE_DIR; safe in Gunicorn pre-fork because all workers share the same directory"
  - "CSV fallback shown as 'Local reference files' not 'Builder API unavailable' — neutral tone, not alarming"
  - "Data source indicator placed on results page only; upload page (index.html) has no indicator per locked decision"

patterns-established:
  - "Source provenance pattern: functions that load from multiple backends return (data, source_label) so callers can expose origin to users"

requirements-completed: [API-02, API-03]

# Metrics
duration: 3min
completed: 2026-02-25
---

# Phase 01 Plan 02: API Integration Wiring Summary

**diskcache.FanoutCache-backed API-first reference set loading with transparent CSV fallback, and a three-state data source indicator on the results page**

## Performance

- **Duration:** 3 min
- **Started:** 2026-02-25T13:27:35Z
- **Completed:** 2026-02-25T13:30:40Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments

- Rewrote `load_cached_reference_sets()` to: check disk cache first, try Builder API second, fall back to CSV silently on any error; returns `(reference_sets, data_source)` tuple
- Added `get_reference_cache()` to `cache_manager.py` providing `diskcache.FanoutCache` with 8 shards for multi-process safety under Gunicorn
- Added data source indicator to results page showing live/cached/local state using VHP4Safety house style colors; index.html is untouched

## Task Commits

Each task was committed atomically:

1. **Task 1: Add diskcache reference cache and wire API-first loading into app.py** - `f968f8e` (feat)
2. **Task 2: Add data source indicator to results template** - `528fdd9` (feat)

## Files Created/Modified

- `app.py` - Added `get_reference_cache`/`fetch_reference_sets_from_api` imports; module-level `_reference_cache` and `REFERENCE_CACHE_KEY`; rewrote `load_cached_reference_sets()` to return `(dict, str)` tuple; updated analyze route to unpack tuple and set `stored_metadata['data_source']`
- `cache_manager.py` - Added `import diskcache` and new `get_reference_cache()` function; existing `CacheManager` and `cached_data_loader` unchanged
- `templates/results.html` - Added three-state data source indicator paragraph after "Background gene set size" line
- `docker-compose.yml` - Added `BUILDER_API_URL=` (empty default) to `ke-enrichment` service environment

## Decisions Made

- `load_cached_reference_sets()` returns a tuple so provenance flows to the template in a single call — no second lookup needed
- `diskcache.FanoutCache` created at module level once per worker process; all workers share the same `/tmp/molaop_cache` directory by design
- Empty `BUILDER_API_URL` raises `ValueError` (implemented in api_service.py in Plan 01), caught by `except Exception`, producing clean CSV fallback with log warning
- Neutral indicator language: "Local reference files" avoids alarming researchers when Builder API is not configured

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

Pre-existing test failures (12 tests in test_database, test_flask_routes, test_report_service) confirmed to predate this plan by stash-checking against the prior commit. No new regressions introduced.

## User Setup Required

None - `BUILDER_API_URL` defaults to empty, triggering CSV fallback automatically. Operators set the actual URL in docker-compose.yml or environment when the Builder service is available.

## Next Phase Readiness

- Full API integration layer is complete: API-first with disk cache and CSV fallback, data source visible to users
- Phase 1 (API Integration) is complete: Plan 01 built the client, Plan 02 wired it into the application
- Builder API bulk export endpoint verification remains an open concern (noted in STATE.md blockers)

---
*Phase: 01-api-integration*
*Completed: 2026-02-25*
