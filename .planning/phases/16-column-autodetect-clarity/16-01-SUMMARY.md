---
phase: 16-column-autodetect-clarity
plan: "01"
subsystem: ui
tags: [flask, jinja2, css, confidence-badge, column-detector, autc-01, autc-02]

# Dependency graph
requires:
  - phase: 12-upload-widget-tech-debt-sweep
    provides: "column_detector singleton imported in app.py; _single_analysis.html column-hint blocks"
provides:
  - "confidence_label Jinja global exposing get_confidence_description to all templates"
  - ".confidence-badge CSS utility classes for High/Medium/Low/Very Low/Undetected states"
  - "Confidence badge + reasons list + static explainer in single-analysis preview column selectors"
  - "TestColumnConfidenceUI test class (5 tests) covering D-01/D-02/D-03 and tour anchor"
affects: [16-02-batch-wizard, phase-13-guided-tour]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "confidence_label Jinja global registered via app.jinja_env.globals for template access to detector methods"
    - "Split column-hint guard: data-stats + reasons inside best_* guard; explainer rendered unconditionally"
    - "TDD with real multipart CSV uploads to /preview — no mocking of detector, tests exercise full stack"

key-files:
  created: []
  modified:
    - "app.py — confidence_label Jinja global registration (AUTC-01)"
    - "static/css/style.css — .confidence-badge--* and .column-hint__reasons/.column-type-explainer CSS"
    - "templates/_single_analysis.html — confidence badges, reasons lists, static explainers in all 4 column groups"
    - "tests/test_flask_routes.py — TestColumnConfidenceUI (5 new tests)"

key-decisions:
  - "Register column_detector.get_confidence_description as a Jinja global (confidence_label) rather than pre-computing labels in the /preview route — keeps route unchanged and makes the helper available to any future template"
  - "Split .column-hint block: stats + reasons inside if-guard for detected columns; .column-type-explainer rendered unconditionally via an else branch (UI-SPEC Component 4 requirement)"
  - "TDD with real CSV multipart uploads (not mocked DataFrames) so tests exercise the actual column detector and produce genuine confidence/reasons data"

patterns-established:
  - "Jinja global pattern: app.jinja_env.globals['name'] = callable — use for detector/service methods that templates need without route changes"
  - "Confidence badge pattern: confidence-badge confidence-badge--{level} with level = confidence_label | lower | replace(' ', '-')"
  - "Explainer always-on pattern: unconditional .column-hint div with .column-type-explainer below the guarded stats+reasons block"

requirements-completed: [AUTC-01, AUTC-02]

# Metrics
duration: 12min
completed: 2026-05-16
---

# Phase 16 Plan 01: Column Autodetect Clarity (Single-Analysis) Summary

**Confidence badges with qualitative labels, percentages, detector reasons lists, and static column-type explainers surfaced in the single-analysis column-selector UI, satisfying AUTC-01 and AUTC-02 for the /preview path**

## Performance

- **Duration:** ~12 min
- **Started:** 2026-05-16T22:15:00Z
- **Completed:** 2026-05-16T22:27:00Z
- **Tasks:** 3 (Task 3 was TDD with RED + GREEN commits)
- **Files modified:** 4

## Accomplishments

- Registered `confidence_label` as a Jinja global so any template can call `confidence_label(float)` to get "High"/"Medium"/"Low"/"Very Low" without route changes
- Added 9 new CSS selectors (`.confidence-badge`, 5 modifier variants, `.column-hint__reasons`, `.column-hint__reasons li`, `.column-type-explainer`) using only VHP4Safety house-palette hex values
- Replaced silent "Auto-detected" spans with filled colour-coded confidence badges across all four column-selector groups in `_single_analysis.html`; added reasons bullet lists (max 5 + overflow label) and static per-column explainer text that renders unconditionally (even when undetected)
- Full test suite: 182 tests pass, 0 failures

## Task Commits

1. **Task 1: Register get_confidence_description as a Jinja global** — `df71abf` (feat)
2. **Task 2: Add confidence-badge and column-hint CSS rules** — `2040d85` (feat)
3. **Task 3: RED — failing TestColumnConfidenceUI tests** — `8f3125a` (test)
4. **Task 3: GREEN — render badges, reasons, explainers in template** — `84ed22d` (feat)

## Files Created/Modified

- `/home/marvin/Documents/Services/molAOP_services/molAOP-analyser/molAOP-analyser/.claude/worktrees/agent-acacb5d5a3f8f075a/app.py` — Added `confidence_label` Jinja global registration after `inject_csrf_token`
- `/home/marvin/Documents/Services/molAOP_services/molAOP-analyser/molAOP-analyser/.claude/worktrees/agent-acacb5d5a3f8f075a/static/css/style.css` — New CONFIDENCE BADGE block with 9 selectors; all hex from VHP4Safety palette; legacy `.confidence-*` classes preserved
- `/home/marvin/Documents/Services/molAOP_services/molAOP-analyser/molAOP-analyser/.claude/worktrees/agent-acacb5d5a3f8f075a/templates/_single_analysis.html` — Extended all four column-selector-group blocks with badge, reasons list, and static explainer; tour anchor preserved
- `/home/marvin/Documents/Services/molAOP_services/molAOP-analyser/molAOP-analyser/.claude/worktrees/agent-acacb5d5a3f8f075a/tests/test_flask_routes.py` — `TestColumnConfidenceUI` class with 5 tests (D-01/D-02/D-03 and tour anchor)

## Decisions Made

- `confidence_label` registered as Jinja global via `app.jinja_env.globals` rather than injecting labels into `template_vars` in the `/preview` route. This keeps the route untouched and makes the helper available to any future template (batch.html, documentation.html, etc.) without route changes.
- Column-hint structure split into two branches: `{% if best_* %}` for stats + reasons, `{% else %}` for undetected case. Both branches include the `.column-type-explainer` so it renders unconditionally per UI-SPEC Component 4.
- TDD tests use real multipart CSV POSTs (no mocked DataFrames) so the column detector actually runs and produces genuine `confidence`, `reasons`, and `data_analysis` values — this makes the tests meaningful behavioral coverage rather than mock-assertion tests.

## Deviations from Plan

None — plan executed exactly as written. All three tasks completed per specification. XSS guard confirmed: no `| safe` filter applied to `col`, `column_name`, or `reasons` values.

## Known Stubs

None.

## Threat Flags

None — no new network endpoints, auth paths, file access patterns, or schema changes introduced. Template changes render existing `ColumnMatch` dataclass fields through Jinja2 default auto-escaping (T-16-01 mitigation confirmed by acceptance criteria grep check).

## Issues Encountered

None — plan executed without blockers.

## Self-Check

- [x] `confidence_label` Jinja global registered and callable: verified with `python -c "import app; print(app.app.jinja_env.globals['confidence_label'](0.85), ...)"` — outputs `High Very Low`
- [x] CSS selectors present: grep count = 4 for key class names
- [x] All 5 TestColumnConfidenceUI tests pass
- [x] Full suite: 182 passed, 0 failed
- [x] Task commits exist: df71abf, 2040d85, 8f3125a, 84ed22d

## Self-Check: PASSED

## Next Phase Readiness

- `confidence_label` Jinja global is available for Plan 16-02 (batch wizard badges) without any additional app.py changes
- `.confidence-badge--*` CSS classes reusable in batch.html JS badge renderer (`classList.add('confidence-badge', 'confidence-badge--' + level)`)
- Single-analysis path (AUTC-01 + AUTC-02) complete; Plan 16-02 delivers the batch wizard half (D-04)

---
*Phase: 16-column-autodetect-clarity*
*Completed: 2026-05-16*
