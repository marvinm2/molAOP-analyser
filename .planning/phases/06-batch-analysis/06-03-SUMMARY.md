---
phase: 06-batch-analysis
plan: "03"
subsystem: batch-analysis
tags: [htmx, flask, templates, batch, progress, css]
dependency_graph:
  requires: ["06-02"]
  provides: ["batch-progress-ui", "batch-summary-page", "batch-condition-results", "htmx-polling"]
  affects: ["app.py", "templates/batch.html", "templates/results.html", "static/css/style.css"]
tech_stack:
  added: []
  patterns: ["htmx-polling-partial", "hx-redirect-on-complete", "breadcrumb-navigation"]
key_files:
  created:
    - templates/batch_progress.html
    - templates/batch_summary.html
  modified:
    - app.py
    - templates/batch.html
    - templates/results.html
    - static/css/style.css
decisions:
  - "batch_status returns htmx partial (not JSON) — htmx manages DOM swap; HX-Redirect triggers browser navigation on complete"
  - "batch.html switches from JS setInterval polling to htmx.process() on dynamically injected hx-get element"
  - "Cancel clears htmx attributes then re-processes container to deactivate htmx before sending cancel POST"
  - "results.html breadcrumb uses metadata.batch_uuid presence check — zero overhead for single-analysis results"
metrics:
  duration: "8 minutes"
  completed_date: "2026-03-01"
  tasks_completed: 2
  files_changed: 6
---

# Phase 06 Plan 03: Batch Progress, Summary, and Condition Results Summary

htmx-based batch progress polling, summary landing page, and individual condition results drill-through — completing the batch analysis user journey.

## What Was Built

### templates/batch_progress.html
htmx partial returned by `/batch/<uuid>/status`. Contains `hx-get` / `hx-trigger="every 1s"` / `hx-swap="outerHTML"` attributes only when `batch.status` is `running` or `pending` — polling stops automatically for terminal states. Renders per-file status rows using `.batch-progress__file--<status>` CSS modifiers.

### templates/batch_summary.html
Extends `base.html`. Displays batch metadata card, harmonisation note (blue callout explaining the background intersection), and a list of `.batch-condition-card` elements in upload order. Each complete condition has a "View Results" link to `/batch/<uuid>/condition/<position>`.

### app.py — four routes added/updated

| Route | Change |
|-------|--------|
| `GET /batch/<uuid>/status` | Refactored from JSON to htmx partial; returns `HX-Redirect` header on complete |
| `POST /batch/<uuid>/cancel` | No change — already correct |
| `GET /batch/<uuid>/summary` | New — renders `batch_summary.html`; redirects to `/batch` if not complete |
| `GET /batch/<uuid>/condition/<int:position>` | New — deserialises ConditionRecord JSON and renders `results.html` |

### templates/batch.html
`showProgressModal()` refactored: injects an htmx-polling div and calls `htmx.process()` to activate it, replacing the previous `setInterval` + JSON fetch approach. Cancel handler removes `hx-get` / `hx-trigger` attributes and calls `htmx.process()` before POSTing cancel.

### templates/results.html
Added breadcrumb block at top of `{% block content %}`: renders `<a class="breadcrumb-link">← Back to batch summary</a>` when `metadata.batch_uuid` is present. No effect on single-analysis results.

### static/css/style.css
Added batch CSS block before `@media` rules: `.batch-modal`, `.batch-modal__content`, `.batch-progress__file[--running|--complete|--failed]`, `.batch-meta`, `.batch-harmonisation-note`, `.batch-conditions`, `.batch-condition-card`, `.btn--small`, `.file-card`, `.tag-row`, `.breadcrumb-link`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed JS polling redirect URL (`/results` → `/summary`)**
- **Found during:** Task 2
- **Issue:** batch.html had `window.location.href = '/batch/' + uuid + '/results'` but the plan specifies the route is `/batch/<uuid>/summary`
- **Fix:** Replaced JS polling entirely with htmx approach; redirect now handled by `HX-Redirect` header from server
- **Files modified:** templates/batch.html, app.py
- **Commit:** f9ce908

**2. [Rule 2 - Missing functionality] Added `.breadcrumb-link` CSS class to style.css**
- **Found during:** Task 2
- **Issue:** Plan specified adding breadcrumb link HTML to results.html but the `.breadcrumb-link` CSS class was not included in the plan's CSS block
- **Fix:** Added `.breadcrumb-link` and `:hover` rules to the batch CSS block
- **Files modified:** static/css/style.css
- **Commit:** 446aa1f

## Self-Check: PASSED

Files verified present:
- templates/batch_progress.html — FOUND
- templates/batch_summary.html — FOUND
- static/css/style.css — contains .batch-modal, .batch-condition-card, .breadcrumb-link — FOUND

Commits verified:
- 446aa1f feat(06-03): add batch progress partial, summary template, and batch CSS — FOUND
- f9ce908 feat(06-03): add batch status, cancel, summary, and condition results routes — FOUND
