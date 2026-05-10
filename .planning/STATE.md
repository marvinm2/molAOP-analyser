---
gsd_state_version: 1.0
milestone: v3.0
milestone_name: Cleanup, Gene-Set Export & UX Gaps
status: executing
stopped_at: Phase 11 UI-SPEC approved
last_updated: "2026-05-10T16:35:22.090Z"
last_activity: 2026-05-10 -- Phase 11 execution started
progress:
  total_phases: 10
  completed_phases: 2
  total_plans: 9
  completed_plans: 5
  percent: 56
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-05-06)

**Core value:** Researchers and regulators can quickly determine which Key Events in a molecular AOP are activated by their transcriptomic data
**Current focus:** Phase 11 — gene-set-export

## Current Position

Phase: 11 (gene-set-export) — EXECUTING
Plan: 1 of 4
Status: Executing Phase 11
Next recommended run: `/gsd-execute-phase 10.1`
Last activity: 2026-05-10 -- Phase 11 execution started

## Performance Metrics

**Velocity:**

- Total plans completed: 3 (v2.0)
- Average duration: 6 minutes
- Total execution time: 17 minutes

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 04 (Shareable URLs) | 2/2 | 12m | 6m |
| 05 (Network Interactivity) | 2/2 | 8m | 4m |
| 06 (Batch Analysis) | 4/4 | 16m | 4m |
| 07 (Comparison Visualisation) | 3/3 | 14m | 4.7m |
| 08 (Comparison UI Polish) | 1/1 | 1m | 1m |

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions logged in PROJECT.md Key Decisions table (9 decisions from v1.0, all good).

Recent decisions affecting current work:

- [v2.0 roadmap]: Phase 4 establishes UUID identity pattern reused by all subsequent phases — do not use integer experiment IDs
- [v2.0 roadmap]: Batch state must live in SQLite from the first commit — never prototype with session cookies (4KB overflow risk)
- [v2.0 roadmap]: Background gene set harmonisation ships with Phase 6 (batch), before comparison (Phase 7) is exposed to users
- [04-01]: UUID generated on-demand at /api/share click, not on every analysis run
- [04-01]: SharedResult table is independent from ExperimentRecord — no migration required
- [04-01]: Expired links return HTTP 410 (Gone) with friendly message, not 404
- [04-02]: get_aop_list() uses diskcache.FanoutCache with 1-week TTL (AOP_LIST_CACHE_KEY = aop_discovery_list_v1)
- [04-02]: load_aop_data() default source changed from csv to sparql — all dynamic AOPs use SPARQL automatically
- [04-02]: Typeahead hidden input carries the aop_id value; text input shows human-readable label only
- [05-01]: FDR filter uses addClass/removeClass('fdr-filtered') with opacity:0.25 — nodes stay in layout, only become semi-transparent
- [05-01]: Drawer uses CSS transition on 'right' property (desktop) and 'bottom' (mobile <768px) — no JS animation
- [05-01]: tableData/keTitleMap/keTypeMap hoisted to global script scope so copyLinkBtn handler and drawer share them without re-parsing template JSON
- [Phase 05-02]: Drop zone active state uses brand blue #307BBF solid border and #e8f4fb light blue tint
- [Phase 05-02]: .drop-zone__remove overrides global button styles with width:auto display:inline-block
- [Phase 05-02]: #cy uses height:500px !important on 768px to override inline style from results.html
- [Phase 05-02]: Upload form is separate from demo form with enctype=multipart/form-data POSTing to /preview
- [06-01]: No cascade='all, delete-orphan' on BatchRecord.conditions — explicit deletion in cleanup_expired_batches() to avoid accidental cascades during development
- [06-01]: BatchRecord and ConditionRecord tables auto-created by Base.metadata.create_all() on startup — no migration scripts needed
- [06-01]: parse_cisplatin_filename() returns filename-as-label for non-matching files, never raises exception
- [06-02]: run_batch() accepts db_url string instead of a session object — SQLAlchemy sessions are not thread-safe, background thread creates its own session
- [06-02]: /batch/upload and /batch/analyze use @csrf.exempt — JS fetch() POSTs require explicit CSRF handling; these are API-style endpoints
- [06-02]: Stop-batch strategy on first condition failure — simpler than partial-success tracking
- [06-02]: parse_cisplatin_filename passed as template kwarg to batch.html so Jinja2 can pre-populate demo checkbox metadata
- [06-03]: batch_status returns htmx partial (not JSON) — htmx manages DOM swap; HX-Redirect triggers browser navigation on complete
- [06-03]: batch.html switches from JS setInterval polling to htmx.process() on dynamically injected hx-get element
- [06-03]: Cancel clears htmx attributes then re-processes container to deactivate htmx before sending cancel POST
- [06-03]: results.html breadcrumb uses metadata.batch_uuid presence check — zero overhead for single-analysis results
- [06-GAP]: htmx loaded from unpkg CDN with SRI hash in base.html before {% block scripts %} — makes htmx global for all child templates
- [06-GAP]: file-card__meta uses parts array joined with em-dash (U+2014) — matches tagging step separator style; rowCount != null handles both null and undefined
- [06-GAP2]: byte-level line count (open rb + sum) gives true row count without full CSV parse; nrows=200 kept for column detection performance
- [06-GAP2]: head_rows uses orient='records' — list of dicts enables JS Object.keys/Object.values rendering without manual column index management
- [07-01]: CONDITION_PALETTE colours assigned in fixed order by upload position — alphabetical sort prevented by explicit pivot column reindex after pivot_table()
- [07-01]: FDR significance threshold at 0.05 for -log10 transform; non-significant cells return None for JSON-safe serialisation
- [07-01]: compare.html dispatches CustomEvent('modechange') so Plans 07-02/03 can listen without coupling
- [Phase 07]: Plotly.newPlot for initial heatmap render, Plotly.react for mode updates — avoids DOM re-creation
- [Phase 07]: data-col carries rendered column index not original condition index so sort works correctly when reference column is hidden in delta mode
- [07-03]: Pie slice sizes normalised as equal shares — each condition gets 100/N percent; slices sized by min(val,10)/10 * sliceMax
- [07-03]: Gene node exclusion in comparison network: skip nodes where ke_type is undefined AND keId not in ke_labels
- [07-03]: KE detail drawer reads live currentMode and ref-condition at tap time to avoid stale closures
- [08-01]: compare drawer inherits .ke-drawer class hierarchy — no duplicate CSS rules needed, only a single ID-scoped width override
- [08-01]: clamp(320px, 30vw, 500px) for responsive drawer width; calc(-1 * clamp(...)) for seamless hidden-state transition
- [08-01]: Mobile bottom-sheet behaviour inherited automatically from .ke-drawer @media query — no duplicate @media block required
- [08-02]: computeSliceData() takes valueMatrix as parameter — caller controls absolute vs delta matrix; helper stays pure
- [08-02]: Negative delta slices rendered at opacity 0.25 (muted gray) not hidden — preserves directional signal
- [08-02]: Delta banner update in modechange listener not setMode() — reference dropdown dispatches modechange without going through setMode()
- [08-02]: pie_opacity_N mapped via data(pie_opacity_N) in stylesheet so opacity survives Cytoscape style updates
- [Phase 09-deployment-hardening]: init_database() called at module level after csrf setup — runs under gunicorn/waitress/python app.py without double-initialisation
- [Phase 09-deployment-hardening]: Batch Analysis nav link uses url_for('batch_page') — no hardcoded /batch path, consistent with other nav items
- [10-01]: KIDNEY_AOP_IDS kept as static 18-element list with inline curation comment block (D-01/D-02) — no dynamic SPARQL derivation; comment contains 'manually curated', 'Last reviewed: 2026-05-06', 'Marvin Martens' verbatim
- [10-01]: AOP:472 CSV rows deleted in-place (9 from aop_ke_map, 8 from aop_ker_edges); strict scope per D-05 — AOP:1, AOP:2, AOP:3, AOP:100 rows untouched
- [10-01]: KER_ID gap (23–30) left intentionally — downstream code keys on (Source_KE, Target_KE) not contiguous KER_IDs, renumbering would be a noisy diff with no behavioural benefit
- [10-01]: kidney-472 entry removed from CASE_STUDY_AOPS but AOP:472 stays in KIDNEY_AOP_IDS — combined Kidney-aop-network view (line ~152 wiring) is unchanged
- [10.1-01]: DEMO_AOP_RECOMMENDATIONS uses fnmatch glob for Cisplatin_Kidney/* pattern — two-pass lookup (exact first, then glob fallback) via get_recommended_aops() classmethod
- [10.1-01]: demos.html uses Jinja2 for-loops; rendered output has 62 cards (2 PXR + 6 curated cisplatin + 54 in details), not 3 source occurrences
- [10.1-01]: /preview extended to read recommended_aops from form POST or query-string GET fallback; parsed list passed to template as recommended_aops variable
- [10.1-02]: Banner uses .card.card--banner modifier (additive) inheriting .card box-shadow — no new card shell needed
- [10.1-02]: aop-filter-toggle rendered inside {% if volcano_data %} — toggle only appears when AOP picker section is visible (correct)
- [10.1-02]: recommended-aops-data JSON data island placed outside {% if volcano_data %} but inside {% if recommended_aops %} — allows JS to read it even before volcano section renders
- [10.1-02]: Test mock updated from MagicMock to real pd.DataFrame — MagicMock did not survive df[[col1,col2,col3]].dropna() slicing needed for volcano_data generation
- [10.1-02]: Test assertion for upload-your-own path uses class="aop-filter-toggle" and <script id="recommended-aops-data" tag — bare strings appear in JS code as getElementById literals

### Roadmap Evolution

- Phase 8 added: Comparison UI Polish — fix KE detail drawer sizing/overlap, missing AOP-Wiki link, pie chart slice sizing, delta mode network clarity
- v3.0 milestone planned (2026-05-06): Phase 10 (AOP source-of-truth cleanup), Phase 11 (gene-set export GMT/CSV + Cytoscape p-value embedding), Phase 12 (upload widget + tech-debt sweep)
- Phase 10.1 inserted after Phase 10: Demos page restructure - move demos off main page (URGENT)

### Pending Todos

1 pending todo:

- **Separate demo setup from live analysis service** (ui) — dedicated demo page, clean main landing page

Resolved:

- ~~**Fix remaining test failures in database and report service**~~ — resolved by commit `0b7c297 fix: resolve all 9 remaining test failures` before Phase 10 began; pytest now reports 124 passed, 0 failed.

### Blockers/Concerns

- [Phase 4]: SPARQL cold-fetch latency for `fetch_all_aops()` unmeasured under production load — implement progressive enhancement (static list server-side, dynamic list on fetch success) as first step
- [Phase 6]: File lifetime management: UUID-scoped directories under uploads/<batch_uuid>/ resolved in 06-01 (create_batch_upload_dir / cleanup_batch_upload_dir)
- [Phase 7]: Cytoscape.js multi-condition node styling confirmed as pie-chart (pie-N-background-color/-size/-opacity) — resolved in 07-03
- [Phase 10]: AOPD-09 dynamic derivation of `KIDNEY_AOP_IDS` may not be feasible if AOP-Wiki lacks a clean "kidney AO" predicate — documenting the curation criterion is the acceptable fallback
- [Ongoing]: Builder API bulk export endpoint needs live verification; rate limit scope unclear for multi-worker production

## Session Continuity

Last session: 2026-05-10T12:33:21.869Z
Stopped at: Phase 11 UI-SPEC approved
Resume file: .planning/phases/11-gene-set-export/11-UI-SPEC.md
