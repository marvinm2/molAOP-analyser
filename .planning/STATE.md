---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: Batch Analysis & Dynamic AOPs
status: active
last_updated: "2026-03-02T12:00:00.000Z"
progress:
  total_phases: 3
  completed_phases: 3
  total_plans: 10
  completed_plans: 10
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-27)

**Core value:** Researchers and regulators can quickly determine which Key Events in a molecular AOP are activated by their transcriptomic data
**Current focus:** Phase 7 — Comparison Visualisation

## Current Position

Phase: 7 of 7 (Comparison Visualisation)
Plan: 1 of 3 in current phase (COMPLETE)
Status: Active
Last activity: 2026-03-02 — 07-01 comparison data service, GET /batch/<uuid>/compare route, compare.html scaffold with tabs and condition colour legend

Progress: [█████████░] 78%

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
| 07 (Comparison Visualisation) | 1/3 | 8m | 8m |

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

### Pending Todos

None.

### Blockers/Concerns

- [Phase 4]: SPARQL cold-fetch latency for `fetch_all_aops()` unmeasured under production load — implement progressive enhancement (static list server-side, dynamic list on fetch success) as first step
- [Phase 6]: File lifetime management: UUID-scoped directories under uploads/<batch_uuid>/ resolved in 06-01 (create_batch_upload_dir / cleanup_batch_upload_dir)
- [Phase 7]: Cytoscape.js multi-condition node styling method (pie-chart vs border-colour) not yet determined — prototype required in Phase 7 planning before committing
- [Ongoing]: Builder API bulk export endpoint needs live verification; rate limit scope unclear for multi-worker production

## Session Continuity

Last session: 2026-03-02
Stopped at: Completed 07-01-PLAN.md — comparison data service, GET /batch/<uuid>/compare route, compare.html scaffold with tabs, control bar, and condition colour legend.
Resume file: None
