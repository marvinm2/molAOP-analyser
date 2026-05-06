---
phase: 04-shareable-urls-aop-discovery
plan: "02"
subsystem: api
tags: [sparql, typeahead, aop-discovery, diskcache, flask-endpoint]

requires:
  - phase: 04-01
    provides: UUID shareable results pattern and SharedResult model

provides:
  - GET /api/aops endpoint with query filtering and 50-item cap
  - AOP discovery service (services/aop_discovery_service.py) with three-tier fallback
  - Typeahead search UI replacing the hardcoded AOP dropdown
  - Dynamic AOP routing to SPARQL via load_aop_data() default source change

affects:
  - Phase 5 (batch analysis) — AOP selection now dynamic, not limited to CASE_STUDY_AOPS
  - Any phase using load_aop_data() — default source is now sparql, not csv

tech-stack:
  added: []
  patterns:
    - Three-tier cache pattern: disk cache (1-week TTL) -> live API fetch -> hardcoded config
    - Typeahead search via hidden input + text input pair, populated from REST endpoint
    - Parallel SPARQL+API fetches using ThreadPoolExecutor in service layer
    - Local imports inside function bodies to avoid circular dependency (data_service, aop_discovery_service)

key-files:
  created:
    - services/aop_discovery_service.py
    - tests/test_aop_discovery.py
  modified:
    - app.py
    - services/data_service.py
    - templates/index.html
    - static/css/style.css

key-decisions:
  - "get_aop_list() uses diskcache.FanoutCache (shared across Gunicorn workers) with 1-week TTL, keyed as aop_discovery_list_v1"
  - "load_aop_data() default source changed from csv to sparql so any dynamically selected AOP automatically routes through SPARQL"
  - "Typeahead hidden input retains the aop_id value; text input only shows the human-readable label — form submission uses hidden input"
  - "Empty /api/aops query returns only AOPs with mapped_ke_count > 0 (on-focus defaults); any query string searches all AOPs"

patterns-established:
  - "AOP discovery three-tier fallback: disk cache -> SPARQL+Builder API parallel fetch -> CASE_STUDY_AOPS hardcoded list"
  - "Typeahead search pattern: <input type=text> for display + <input type=hidden> for value, populated via fetch()"

requirements-completed:
  - AOPD-01
  - AOPD-02
  - AOPD-03
  - AOPD-04
  - AOPD-05
  - AOPD-06

duration: 6min
completed: 2026-02-27
---

# Phase 04 Plan 02: AOP Discovery Service + Typeahead UI Summary

**SPARQL+Builder API AOP discovery service with three-tier cache fallback, /api/aops REST endpoint, and typeahead search replacing the hardcoded AOP dropdown**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-27T19:05:22Z
- **Completed:** 2026-02-27T19:11Z
- **Tasks:** 5 of 5
- **Files modified:** 6

## Accomplishments

- Built `services/aop_discovery_service.py` with parallel SPARQL queries (all AOPs + KE membership) and Builder API fetch, three-tier fallback (disk cache -> live fetch -> hardcoded config), and 1-week TTL caching
- Added `GET /api/aops` Flask endpoint that returns mapped AOPs on empty query and filters by ID/title on typed query, capped at 50 results
- Replaced the hardcoded `<select id="aop_selection">` dropdown with a typeahead search box that queries `/api/aops` with debouncing and inline error/retry handling
- Fixed `load_aop_data()` default source from `"csv"` to `"sparql"` enabling any dynamically discovered AOP to be analyzed without requiring a config entry
- Added 21 pytest tests covering all service functions, endpoint behavior, and routing logic — all pass

## Task Commits

1. **Task 1: Create aop_discovery_service.py** - `0f1b79b` (feat)
2. **Task 2: Add GET /api/aops endpoint** - `4c225a4` (feat)
3. **Task 3: Fix load_aop_data() default source** - `7a479b6` (fix)
4. **Task 4: Replace dropdown with typeahead** - `49de7c5` (feat)
5. **Task 5: Write tests** - `5681304` (test)

## Files Created/Modified

- `services/aop_discovery_service.py` — AOP discovery service: SPARQL fetch, Builder API fetch, three-tier fallback, cache logic
- `tests/test_aop_discovery.py` — 21 tests for discovery service, /api/aops endpoint, and load_aop_data() routing
- `app.py` — Added `from services.aop_discovery_service import get_aop_list` and `GET /api/aops` route
- `services/data_service.py` — Single-line change: default source `"csv"` -> `"sparql"` in `load_aop_data()`
- `templates/index.html` — Replaced `<select>` with typeahead text+hidden input pair and JavaScript fetch logic
- `static/css/style.css` — Added `.typeahead-dropdown` and related styles

## Decisions Made

- Used `diskcache.FanoutCache` (already a project dependency) for the AOP list cache with 8-shard config for Gunicorn multi-worker safety
- `_hardcoded_aop_list()` sets `mapped_ke_count=1` for all configured AOPs so they appear in the on-focus defaults even when live APIs are unavailable
- Typeahead JavaScript is wrapped in an IIFE and guards against the case where elements are not present (the typeahead section only renders when `volcano_data` is set)
- `/api/aops` endpoint is not CSRF-exempt since it is a GET request with no state changes

## Deviations from Plan

None — plan executed exactly as written.

## Issues Encountered

- Test patching required targeting `services.sparql_service.fetch_aop_ke_data_cached` (not `services.data_service.*`) because `load_aop_data()` uses a local import inside the function body. Fixed by patching the correct module.
- `test_with_url_returns_normalised_ids` required patching `services.api_service` module-level attributes since `fetch_mapped_ke_ids_from_builder` does a local `from services.api_service import ...` inside the function.

## Self-Check: PASSED

Files verified:
- FOUND: services/aop_discovery_service.py
- FOUND: tests/test_aop_discovery.py

Commits verified:
- 0f1b79b: feat(04-02): add AOP discovery service with three-tier fallback
- 4c225a4: feat(04-02): add GET /api/aops endpoint for typeahead suggestions
- 7a479b6: fix(04-02): default load_aop_data() to SPARQL for dynamic AOPs
- 49de7c5: feat(04-02): replace AOP dropdown with typeahead search in index.html
- 5681304: test(04-02): add pytest suite for AOP discovery service and /api/aops endpoint

Tests: 21 new tests pass, 13 pre-existing failures unchanged.

## Next Phase Readiness

- Phase 4 is complete: UUID shareable links (04-01) and dynamic AOP discovery (04-02) are both done
- Phase 5 (batch analysis) can rely on AOP selection no longer being limited to CASE_STUDY_AOPS
- SPARQL routing is the default for all dynamic AOPs — no config changes needed for new AOPs

---
*Phase: 04-shareable-urls-aop-discovery*
*Completed: 2026-02-27*
