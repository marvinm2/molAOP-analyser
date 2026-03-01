---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: Batch Analysis & Dynamic AOPs
status: unknown
last_updated: "2026-03-01T18:52:05Z"
progress:
  total_phases: 3
  completed_phases: 2
  total_plans: 7
  completed_plans: 6
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-27)

**Core value:** Researchers and regulators can quickly determine which Key Events in a molecular AOP are activated by their transcriptomic data
**Current focus:** Phase 6 — Batch Analysis

## Current Position

Phase: 6 of 7 (Batch Analysis)
Plan: 2 of 3 in current phase (COMPLETE)
Status: Active
Last activity: 2026-03-01 — 06-02 Batch upload wizard, orchestration service, and Flask routes complete

Progress: [██████░░░░] 57%

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
| 06 (Batch Analysis) | 2/3 | 7m | 3.5m |

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

### Pending Todos

None.

### Blockers/Concerns

- [Phase 4]: SPARQL cold-fetch latency for `fetch_all_aops()` unmeasured under production load — implement progressive enhancement (static list server-side, dynamic list on fetch success) as first step
- [Phase 6]: File lifetime management: UUID-scoped directories under uploads/<batch_uuid>/ resolved in 06-01 (create_batch_upload_dir / cleanup_batch_upload_dir)
- [Phase 7]: Cytoscape.js multi-condition node styling method (pie-chart vs border-colour) not yet determined — prototype required in Phase 7 planning before committing
- [Ongoing]: Builder API bulk export endpoint needs live verification; rate limit scope unclear for multi-worker production

## Session Continuity

Last session: 2026-03-01
Stopped at: Completed 06-02-PLAN.md — Batch upload wizard (batch.html), orchestration service (batch_service.py), and Flask routes (/batch, /batch/upload, /batch/analyze) complete
Resume file: None
