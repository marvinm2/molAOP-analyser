---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: Batch Analysis & Dynamic AOPs
status: unknown
last_updated: "2026-02-28T14:59:11.322Z"
progress:
  total_phases: 2
  completed_phases: 2
  total_plans: 4
  completed_plans: 4
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-27)

**Core value:** Researchers and regulators can quickly determine which Key Events in a molecular AOP are activated by their transcriptomic data
**Current focus:** Phase 5 — Network Interactivity & Responsive Layout

## Current Position

Phase: 5 of 7 (Network Interactivity & Responsive Layout)
Plan: 2 of 2 in current phase (COMPLETE)
Status: Active
Last activity: 2026-02-28 — 05-02 Responsive layout and drag-and-drop upload complete

Progress: [████░░░░░░] 40%

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

### Pending Todos

None.

### Blockers/Concerns

- [Phase 4]: SPARQL cold-fetch latency for `fetch_all_aops()` unmeasured under production load — implement progressive enhancement (static list server-side, dynamic list on fetch success) as first step
- [Phase 6]: File lifetime management strategy for session-scoped upload directories not yet specified — resolve in Phase 6 planning before implementation
- [Phase 7]: Cytoscape.js multi-condition node styling method (pie-chart vs border-colour) not yet determined — prototype required in Phase 7 planning before committing
- [Ongoing]: Builder API bulk export endpoint needs live verification; rate limit scope unclear for multi-worker production

## Session Continuity

Last session: 2026-02-28
Stopped at: Completed 05-02-PLAN.md — Responsive layout (768px/1024px breakpoints) and drag-and-drop upload implemented
Resume file: None
