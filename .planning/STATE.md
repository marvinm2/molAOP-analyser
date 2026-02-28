---
gsd_state_version: 1.0
milestone: v2.0
milestone_name: Batch Analysis & Dynamic AOPs
status: active
last_updated: "2026-02-27"
progress:
  total_phases: 4
  completed_phases: 0
  total_plans: 10
  completed_plans: 0
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-27)

**Core value:** Researchers and regulators can quickly determine which Key Events in a molecular AOP are activated by their transcriptomic data
**Current focus:** Phase 5 — Network Interactivity & Responsive Layout

## Current Position

Phase: 5 of 7 (Network Interactivity & Responsive Layout)
Plan: 1 of 1 in current phase
Status: Active
Last activity: 2026-02-28 — 05-01 Network interactivity (expand-all, KE drawer, FDR filter) complete

Progress: [███░░░░░░░] 30%

## Performance Metrics

**Velocity:**
- Total plans completed: 3 (v2.0)
- Average duration: 6 minutes
- Total execution time: 17 minutes

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 04 (Shareable URLs) | 2/2 | 12m | 6m |
| 05 (Network Interactivity) | 1/1 | 5m | 5m |

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

### Pending Todos

None.

### Blockers/Concerns

- [Phase 4]: SPARQL cold-fetch latency for `fetch_all_aops()` unmeasured under production load — implement progressive enhancement (static list server-side, dynamic list on fetch success) as first step
- [Phase 6]: File lifetime management strategy for session-scoped upload directories not yet specified — resolve in Phase 6 planning before implementation
- [Phase 7]: Cytoscape.js multi-condition node styling method (pie-chart vs border-colour) not yet determined — prototype required in Phase 7 planning before committing
- [Ongoing]: Builder API bulk export endpoint needs live verification; rate limit scope unclear for multi-worker production

## Session Continuity

Last session: 2026-02-28
Stopped at: Completed 05-01-PLAN.md — Network interactivity (expand-all, KE drawer, FDR filter) implemented
Resume file: None
