---
phase: 02-ui-ux-and-documentation
plan: 02
subsystem: ui
tags: [jinja2, flask, html, css, templates, documentation]

requires:
  - phase: 02-01
    provides: base.html template, style.css brand foundation, card component, site-nav

provides:
  - index.html extending base.html with card--narrow layout and migrated CSS classes
  - documentation.html extending base.html with three required sections only

affects:
  - 02-03 (results.html migration — same pattern)

tech-stack:
  added: []
  patterns:
    - "Jinja2 template inheritance: all pages extend base.html via {% extends %}"
    - "card--narrow wrapper for centred single-column page content"
    - "Page-specific CSS in {% block head_extra %} style tags"
    - "alert-band class for development-status notice"
    - "info-callout class for blue-left-border information boxes"
    - "dose-grid / dose-label classes for cisplatin selector grid"
    - "btn--option / btn--option-pct classes for threshold quick-option buttons"
    - "All scripts moved to {% block scripts %} at page bottom"

key-files:
  created: []
  modified:
    - templates/index.html
    - templates/documentation.html

key-decisions:
  - "alert-band replaces .subheader h2 — more semantic for a page-level notice outside the card"
  - "info-callout class defined in head_extra for the demo explanation box (reuses documentation.html info-box pattern with a distinct name to avoid collision)"
  - "Volcano plot script wrapped in {% if volcano_data %} guard inside {% block scripts %} — avoids JS reference errors when data absent"
  - "documentation.html scoped to exactly three sections: Input File Format, Statistical Methods, Interpreting Results — Getting Started, Understanding AOPs, Glossary, FAQ removed per plan"

patterns-established:
  - "Template inheritance: {% extends 'base.html' %} + {% block content %} + {% block scripts %}"
  - "Conditional sections (preview, volcano) each get their own card--narrow div"
  - "Page-specific utility classes defined in {% block head_extra %} rather than modifying global style.css"

requirements-completed: [UI-04, DOC-01]

duration: 3min
completed: 2026-02-25
---

# Phase 02 Plan 02: index.html and documentation.html Migration Summary

**index.html and documentation.html migrated to Jinja2 template inheritance via base.html, with card--narrow layout, migrated inline styles, and documentation scoped to three required sections**

## Performance

- **Duration:** 3 min
- **Started:** 2026-02-25T14:41:08Z
- **Completed:** 2026-02-25T14:44:48Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- index.html extends base.html — old DOCTYPE/header/footer/page-container boilerplate removed
- All content in index.html wrapped in card--narrow components (main card, preview card, volcano card)
- Inline styles migrated to CSS classes: info-callout, dose-grid, dose-label, btn--option, btn--option-pct, threshold-controls, metadata-grid
- All JS moved to {% block scripts %} block; volcano plot script guarded by {% if volcano_data %}
- documentation.html extends base.html — scoped to exactly three required sections (Input File Format, Statistical Methods, Interpreting Results)
- Getting Started, Understanding AOPs, Glossary, and FAQ sections removed per plan decision
- Quick navigation updated to list only the three remaining sections with correct anchor IDs
- Both routes (/ and /documentation) return HTTP 200

## Task Commits

Each task was committed atomically:

1. **Task 1: Refactor index.html to extend base.html with single centred card layout** - `1829bb2` (feat)
2. **Task 2: Refactor documentation.html to extend base.html with three required sections** - `bb4f0fe` (feat)

## Files Created/Modified
- `templates/index.html` - Migrated to Jinja2 inheritance, card layout, migrated styles, preserved all form/JS functionality
- `templates/documentation.html` - Migrated to Jinja2 inheritance, scoped to 3 sections, removed supplementary sections

## Decisions Made
- alert-band replaces the .subheader h2 element — more semantically appropriate for a page-level notice displayed outside the card
- info-callout class defined in head_extra for the demo explanation box; reuses the same visual pattern as documentation.html's info-box but with a distinct name to prevent CSS collision if both pages are ever on the same layout
- Volcano plot Plotly script wrapped in {% if volcano_data %} guard inside {% block scripts %} to prevent JS reference errors when the page is loaded without volcano data
- documentation.html scoped to exactly three sections per locked plan decisions; the removed sections (Getting Started, Understanding AOPs, Glossary, FAQ) were supplementary content not required for tool usage documentation

## Deviations from Plan

None — plan executed exactly as written.

## Issues Encountered
None.

## User Setup Required
None — no external service configuration required.

## Next Phase Readiness
- Two of four templates now fully migrated to base.html inheritance
- index.html and documentation.html complete — ready for results.html migration (02-03)
- Card--narrow pattern established for 02-03 to follow

---
*Phase: 02-ui-ux-and-documentation*
*Completed: 2026-02-25*
