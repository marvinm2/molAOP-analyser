---
gsd_state_version: 1.0
milestone: v2.0
milestone_name: Batch Analysis & Dynamic AOPs
status: unknown
last_updated: "2026-03-03T07:34:41.212Z"
progress:
  total_phases: 6
  completed_phases: 6
  total_plans: 15
  completed_plans: 15
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-27)

**Core value:** Researchers and regulators can quickly determine which Key Events in a molecular AOP are activated by their transcriptomic data
**Current focus:** Phase 8 — Comparison UI Polish

## Current Position

Phase: 9 of 9 (Deployment Hardening)
Plan: 1 of 1 in current phase (COMPLETE)
Status: Phase 9 complete — all plans done
Last activity: 2026-03-02 — 09-01 Batch Analysis nav link; init_database() moved to module level

Progress: [██████████] 100%

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

### Roadmap Evolution

- Phase 8 added: Comparison UI Polish — fix KE detail drawer sizing/overlap, missing AOP-Wiki link, pie chart slice sizing, delta mode network clarity

### Pending Todos

2 pending todos:
- **Fix remaining test failures in database and report service** (testing) — 9 pre-existing failures from shallow mocks
- **Separate demo setup from live analysis service** (ui) — dedicated demo page, clean main landing page

### Blockers/Concerns

- [Phase 4]: SPARQL cold-fetch latency for `fetch_all_aops()` unmeasured under production load — implement progressive enhancement (static list server-side, dynamic list on fetch success) as first step
- [Phase 6]: File lifetime management: UUID-scoped directories under uploads/<batch_uuid>/ resolved in 06-01 (create_batch_upload_dir / cleanup_batch_upload_dir)
- [Phase 7]: Cytoscape.js multi-condition node styling confirmed as pie-chart (pie-N-background-color/-size/-opacity) — resolved in 07-03
- [Ongoing]: Builder API bulk export endpoint needs live verification; rate limit scope unclear for multi-worker production

## Session Continuity

Last session: 2026-03-02
Stopped at: Completed 09-01-PLAN.md — Batch Analysis nav link and module-level init_database() for WSGI compatibility.
Resume file: None
