---
phase: 01-api-integration
plan: 01
subsystem: api
tags: [requests, diskcache, builder-api, ke-wp-mappings, pandas]

# Dependency graph
requires: []
provides:
  - Builder API client with paginated fetch, retry, and KE ID normalisation
  - Refactored load_reference_sets() accepting optional ke_wp_df DataFrame
  - Config settings BUILDER_API_URL, BUILDER_API_TIMEOUT, CACHE_DIR, CACHE_TTL

affects:
  - 01-02 (integration layer that wires API-first loading with CSV fallback)

# Tech tracking
tech-stack:
  added: [requests, diskcache]
  patterns:
    - Optional DataFrame injection into load_reference_sets() for dual CSV/API data sources
    - Retry-backed requests.Session via HTTPAdapter + urllib3 Retry
    - KE ID normalisation at API boundary (space to colon format)

key-files:
  created:
    - services/api_service.py
  modified:
    - helpers.py
    - config.py
    - requirements.txt

key-decisions:
  - "load_reference_sets() accepts ke_wp_df=None to avoid duplicating merge logic across CSV and API paths"
  - "KE ID normalisation (space to colon) happens in fetch_reference_sets_from_api() before DataFrame construction, not inside load_reference_sets()"
  - "Empty BUILDER_API_URL raises ValueError (disabled intent) not a network error, to allow caller to distinguish config vs runtime failures"
  - "Retry strategy: 3 retries, backoff_factor=1.0, status_forcelist=[429, 500, 502, 503, 504]"

patterns-established:
  - "API boundary normalisation: transform external format at fetch time, pass normalised data into shared pipeline"
  - "Dual-source loader: single function handles both file and in-memory DataFrame via optional param"

requirements-completed: [API-01]

# Metrics
duration: 6min
completed: 2026-02-25
---

# Phase 01 Plan 01: Builder API Client Summary

**Paginated Builder API client with retry/backoff and KE ID normalisation, plus dual-source load_reference_sets() accepting CSV or pre-built DataFrame**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-25T13:18:25Z
- **Completed:** 2026-02-25T13:24:00Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments

- Created `services/api_service.py` with `_make_api_session()`, `fetch_all_ke_wp_mappings()`, and `fetch_reference_sets_from_api()` providing full paginated API fetch with exponential backoff
- Refactored `helpers.py` `load_reference_sets()` to accept optional `ke_wp_df` DataFrame, enabling API-sourced data to flow through the same WP-to-gene merge pipeline as CSV data
- Added BUILDER_API_URL, BUILDER_API_TIMEOUT, CACHE_DIR, CACHE_TTL to Config class with environment variable support; added `requests` and `diskcache` to requirements.txt

## Task Commits

Each task was committed atomically:

1. **Task 1: Add diskcache dependency and Builder API config** - `12061bb` (chore)
2. **Task 2: Create API service and refactor helpers.py** - `84acb23` (feat)

## Files Created/Modified

- `services/api_service.py` - Builder API client: retry session, paginated KE-WP fetch, KE ID normalisation, reference set orchestration
- `helpers.py` - Added optional `ke_wp_df` parameter to `load_reference_sets()`; backward compatible with all existing CSV callers
- `config.py` - Added BUILDER_API_URL, BUILDER_API_TIMEOUT, CACHE_DIR, CACHE_TTL class attributes
- `requirements.txt` - Added `requests` and `diskcache` dependencies

## Decisions Made

- `load_reference_sets()` accepts `ke_wp_df=None` to keep all WP-to-gene merge logic in one place; both CSV and API paths converge here
- KE ID normalisation (`"KE 55"` to `"KE:55"`) occurs at the API fetch boundary before DataFrame construction, so the merge function always sees consistent format
- Empty `BUILDER_API_URL` raises `ValueError` (not a network error) so the integration layer can clearly distinguish "not configured" from "API down"

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required at this stage. `BUILDER_API_URL` defaults to empty, and Plan 02 will wire the fallback logic.

## Next Phase Readiness

- Builder API client is complete and importable
- `load_reference_sets()` is ready to receive API-sourced DataFrames
- Plan 02 can now implement the integration layer: try API, catch ValueError/network error, fall back to CSV

---
*Phase: 01-api-integration*
*Completed: 2026-02-25*
