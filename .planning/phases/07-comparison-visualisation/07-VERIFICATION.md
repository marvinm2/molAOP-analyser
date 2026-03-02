---
phase: 07-comparison-visualisation
verified: 2026-03-02T13:15:00Z
status: human_needed
score: 5/5 success criteria verified
re_verification:
  previous_status: gaps_found
  previous_score: 0/5 (all blocked by Jinja2 error)
  gaps_closed:
    - "Navigating to /batch/<uuid>/compare for a completed batch returns 200 with a rendered compare.html page — Jinja2 TemplateSyntaxError on line 342 fixed: @returns {Object} with keys trace, layout, and config (double braces removed)"
  gaps_remaining: []
  regressions: []
human_verification:
  - test: "Tab switching behaviour — Heatmap, Table, Network tabs"
    expected: "Clicking each tab shows its panel, hides the others, and updates the active styling. Switching to the heatmap tab triggers Plotly.Plots.resize(). Switching to the network tab for the first time calls initComparisonNetwork()."
    why_human: "JavaScript event wiring and Cytoscape lazy-init on tab switch cannot be verified statically."
  - test: "Absolute/delta mode toggle"
    expected: "Clicking Delta shows the reference condition dropdown; clicking Absolute hides it and dispatches CustomEvent('modechange') that updates the heatmap, table, and network."
    why_human: "CustomEvent dispatch and consumer coordination requires browser execution."
  - test: "Heatmap rendering — grey null cells and yellow-to-red significant gradient"
    expected: "Non-significant cells (FDR > 0.05) appear grey; significant cells show a yellow-to-red gradient scaled by -log10(FDR)."
    why_human: "Plotly.js rendering and colour scale cannot be verified without a browser."
  - test: "Cytoscape pie-chart KE nodes on Network tab"
    expected: "KE nodes show pie slices sized proportionally to -log10(FDR), coloured by the VHP4Safety condition palette. Non-significant portions show grey background. KER edges appear as thin grey arrows."
    why_human: "Cytoscape canvas rendering requires browser execution."
  - test: "Single-condition / comparison view toggle on network"
    expected: "Clicking 'Single Condition' disables pie display and colours nodes solid by condition significance. Selecting a different condition from the dropdown updates node colours."
    why_human: "Network view state toggling requires browser execution."
  - test: "KE detail drawer on node click"
    expected: "Tapping a KE node slides in the drawer showing KE title, a per-condition horizontal bar chart with significance threshold line, and an AOP-Wiki link."
    why_human: "Requires browser + Cytoscape tap event execution."
---

# Phase 07: Comparison Visualisation Verification Report

**Phase Goal:** Users can see all batch conditions side-by-side on a shared AOP network and heatmap, toggle between absolute and delta modes, and understand relative significance across conditions at a glance
**Verified:** 2026-03-02T13:15:00Z
**Status:** human_needed
**Re-verification:** Yes — after gap closure (Jinja2 fix, commit 775d0a9)

## Goal Achievement

### Observable Truths (from ROADMAP.md Success Criteria)

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | User sees a Plotly heatmap with KEs as rows and conditions as columns, coloured by -log10(FDR) significance | ? HUMAN NEEDED | Jinja2 fix confirmed (line 342 now `@returns {Object} with keys trace, layout, and config` — single braces). compare.html parses cleanly. Plotly.newPlot, updateHeatmap, computeDeltaMatrix all present and substantive. Requires browser to confirm render. |
| 2 | User sees a comparison enrichment table with one row per KE and one column per condition, showing FDR values | ? HUMAN NEEDED | renderComparisonTable, sortTableByColumn, toExponential all present. Page now renders. Requires browser to confirm. |
| 3 | User can click a toggle to switch between single-condition and comparison view on the AOP network | ? HUMAN NEEDED | btn-single-view, buildSingleConditionStylesheet, initComparisonNetwork all present. Page now renders. Requires browser to confirm. |
| 4 | User can select a reference condition and heatmap and network switch to delta mode | ? HUMAN NEEDED | computeDeltaMatrix, modechange CustomEvent dispatch and listener all present and wired. Page now renders. Requires browser to confirm. |
| 5 | A condition colour legend using VHP4Safety palette colours is visible on the comparison page | ? HUMAN NEEDED | CONDITION_PALETTE (8 colours), condition-legend div, renderLegend() function all verified present. Page now renders. Requires browser to confirm. |

**Score:** 5/5 truths have complete, wired, substantive implementation — all blocked by human verification only (no code defects remain).

**Previous gap resolved:** The sole blocker was a Jinja2 TemplateSyntaxError on line 342 of `compare.html`. The JSDoc `@returns` annotation previously read `@returns {{ trace: Object, layout: Object, config: Object }}` — Jinja2 interpreted the double curly braces as an invalid template expression. The fix (commit 775d0a9) rewrote this as `@returns {Object} with keys trace, layout, and config`, using single braces that Jinja2 ignores. All three key templates (`compare.html`, `batch_summary.html`, `base.html`) now parse without errors.

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `services/comparison_service.py` | Comparison data pivot, -log10 transform, CONDITION_PALETTE | VERIFIED | 143+ lines, fully substantive. build_comparison_matrix() pivots FDR matrix, applies neg-log10 transform, sorts KEs by significance, returns required dict. CONDITION_PALETTE has exactly 8 VHP4Safety colours. Import confirmed at app.py line 38. |
| `templates/compare.html` | Tabbed comparison page with control bar, legend, heatmap, table, network | VERIFIED | 1150+ lines. All required functions present (Plotly.newPlot, renderComparisonTable, initComparisonNetwork, buildComparisonStylesheet, openKeDrawer). Jinja2 parses cleanly — confirmed via `env.parse()`. No double-brace JSDoc patterns remain. |
| `templates/batch_summary.html` | Compare All button linking to /batch/<uuid>/compare | VERIFIED | Line 55: `<a href="/batch/{{ batch.uuid }}/compare" class="btn btn--primary">Compare All Conditions</a>` — parses cleanly. |
| `app.py` | GET /batch/<uuid>/compare route (batch_compare) | VERIFIED | Route at line 1431. Queries BatchRecord, handles 404 and non-complete redirect, builds comparison_data, serialises to JSON, extracts first_network (lines 1464-1474), renders compare.html with comparison_data_json and network_json. |
| `static/css/style.css` | Comparison page CSS classes | VERIFIED | All required classes present: .compare-subtitle, .compare-controls, .compare-legend, .compare-legend__item, .compare-legend__swatch, .compare-tabs, .compare-tab, .compare-tab--active, .batch-summary-actions, responsive rules at 768px. |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| templates/batch_summary.html | /batch/uuid/compare | href on Compare All button | WIRED | Line 55: `href="/batch/{{ batch.uuid }}/compare"` |
| app.py batch_compare route | services/comparison_service.py | build_comparison_matrix() call | WIRED | app.py line 38 imports; line 1461 calls build_comparison_matrix(conditions) |
| app.py batch_compare route | templates/compare.html | render_template with comparison_data_json | WIRED | Line 1482: render_template('compare.html', ..., comparison_data_json=comparison_data_json, network_json=...) |
| compare.html script | Plotly.js CDN | Plotly.newPlot('heatmap-div', ...) | WIRED | Line 436: Plotly.newPlot; CDN loaded in head_extra |
| compare.html modechange listener | Plotly.react and table re-render | CustomEvent 'modechange' | WIRED | Line 1122: addEventListener('modechange', ...) calls updateHeatmap() + updateTable() + updateNetwork() |
| compare.html network init | Cytoscape.js CDN | cytoscape({container, elements, style, layout}) | WIRED | Line 961 (approx): cyInstance = cytoscape({...}) in initComparisonNetwork() |
| compare.html modechange listener | network data update + cy.style().update() | modechange event -> updateNetwork() | WIRED | Line 927: cyInstance.style(...).update(); savedPositions restored on mode switch |
| compare.html network node tap | KE detail panel with bar chart | cy.on('tap', 'node') handler | WIRED | cyInstance.on('tap', 'node', ...) at line 976 calls openKeDrawer() |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| COMP-01 | 07-02 | KE x condition heatmap showing -log10(FDR) significance | VERIFIED (code) / HUMAN for render | Plotly.newPlot, computeDeltaMatrix, updateHeatmap present; compare.html parses; page reachable |
| COMP-02 | 07-02 | Comparison enrichment table (rows=KEs, columns=conditions) | VERIFIED (code) / HUMAN for render | renderComparisonTable, sortTableByColumn, toExponential all present |
| COMP-03 | 07-03 | Multi-condition AOP network overlay on shared KE skeleton | VERIFIED (code) / HUMAN for render | initComparisonNetwork, buildComparisonElements (via buildComparisonStylesheet), pie-1-background-color, dagre LR layout, savedPositions all present |
| COMP-04 | 07-03 | User can toggle between single-condition and comparison view | VERIFIED (code) / HUMAN for render | btn-single-view, buildSingleConditionStylesheet pattern, networkViewMode state present |
| COMP-05 | 07-03 | Delta mode showing change relative to user-selected reference condition | VERIFIED (code) / HUMAN for render | computeDeltaMatrix, modechange listener, RdBu colorscale, delta network update, drawer delta refresh all present |
| COMP-06 | 07-01 | Condition colour legend using VHP4Safety palette | VERIFIED (code) / HUMAN for render | CONDITION_PALETTE 8 colours (confirmed via import), .compare-legend CSS, renderLegend() function present |

All 6 requirements (COMP-01 through COMP-06) are accounted for. No orphaned requirements. REQUIREMENTS.md traceability table marks all as Complete.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| None | — | — | — | The previously blocking pattern (double-brace JSDoc on line 342) has been resolved. No remaining blockers or stubs detected. The `return {}` in comparison_service.py line 88 is a legitimate early-exit guard, not a stub. The `return null` patterns in compare.html are valid JavaScript null returns for missing data cases. |

### Human Verification Required

### 1. Plotly Heatmap Rendering

**Test:** Navigate to a completed batch's comparison page (`/batch/<uuid>/compare`). Verify the Heatmap tab shows a Plotly heatmap.
**Expected:** KEs appear as rows, conditions as columns. Non-significant cells (FDR > 0.05) appear grey. Significant cells show yellow-to-red gradient proportional to -log10(FDR). Hovering a cell shows KE title, condition name, and -log10(FDR) value.
**Why human:** Plotly.js canvas rendering and colour interpolation require browser execution.

### 2. Absolute/Delta Mode Toggle

**Test:** On the comparison page, click "Delta", select a reference condition. Verify heatmap, table, and network all update.
**Expected:** Heatmap switches to diverging blue-white-red RdBu scale with zmid=0. Table shows delta -log10(FDR) values with sign. Reference column is hidden when >= 3 conditions. Network pie slices update to delta magnitudes.
**Why human:** CustomEvent dispatch and cross-component coordination requires browser execution.

### 3. Sortable Comparison Table

**Test:** Switch to the Table tab. Verify FDR values display in scientific notation. Click a condition column header.
**Expected:** Cells with FDR < 0.05 show bold text and a 15%-opacity palette background tint. Clicking a header sorts rows by that condition ascending; clicking again reverses. Sort indicator (triangle) appears on the active column header.
**Why human:** DOM manipulation and sort state require browser execution.

### 4. Cytoscape Pie-Chart Network

**Test:** Click the Network tab on the comparison page.
**Expected:** Pie-chart KE nodes appear with condition-coloured slices sized by -log10(FDR). Non-significant portions show grey. KER edges appear as thin grey arrows. Dagre LR layout is computed.
**Why human:** Cytoscape.js canvas rendering and dagre layout require browser execution.

### 5. Single-Condition / Comparison View Toggle

**Test:** On the Network tab, click "Single Condition". Select different conditions from the dropdown.
**Expected:** Nodes switch from pie charts to solid colour (condition palette colour = significant, grey = non-significant). Switching back to "Comparison" restores pie charts without re-running dagre layout.
**Why human:** Cytoscape stylesheet switching and view state require browser execution.

### 6. KE Detail Drawer

**Test:** Click a KE node on the network. Verify the detail panel slides in.
**Expected:** Drawer shows KE title, horizontal bar chart with per-condition -log10(FDR) bars coloured by condition palette, dashed significance threshold line at x=1.301 (-log10(0.05)), and an AOP-Wiki link button for the KE. Clicking the close button or the network background closes the drawer.
**Why human:** Cytoscape tap events, Plotly bar chart rendering, and CSS slide animation require browser execution.

## Gaps Summary

No gaps remain. The single blocker from the initial verification — a Jinja2 `TemplateSyntaxError` caused by a JSDoc `@returns` annotation using double curly braces on line 342 of `compare.html` — has been fixed in commit 775d0a9. The fix changed `@returns {{ trace: Object, layout: Object, config: Object }}` to `@returns {Object} with keys trace, layout, and config`.

Jinja2 now parses all three templates (`compare.html`, `batch_summary.html`, `base.html`) without error. The `GET /batch/<uuid>/compare` route is reachable and will render the comparison page.

All COMP-01 through COMP-06 requirements have complete, substantive, and correctly-wired implementations. The remaining verification items (items 1-6 above) require browser execution and are standard human verification — none indicate code defects.

---

_Initial verified: 2026-03-02T12:43:22Z_
_Re-verified: 2026-03-02T13:15:00Z_
_Verifier: Claude (gsd-verifier)_
