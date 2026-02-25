# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-25)

**Core value:** Researchers and regulators can quickly determine which Key Events in a molecular AOP are activated by their transcriptomic data
**Current focus:** Phase 3 — Network Performance

## Current Position

Phase: 3 of 3 (Network Performance)
Plan: 2 of 2 in current phase — Plan 03-02 complete
Status: Phase 3 complete — all plans done
Last activity: 2026-02-25 — Completed 03-02 (overlay zoom controls, transparent PNG export)

Progress: [██████████] 100%

## Performance Metrics

**Velocity:**
- Total plans completed: 6
- Average duration: 2.8 min
- Total execution time: 0.28 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-api-integration | 2 | 9 min | 4.5 min |
| 02-ui-ux-and-documentation | 3 | 10 min | 3.3 min |
| 03-network-performance | 2 | 3 min | 1.5 min |

**Recent Trend:**
- Last 5 plans: 2 min, 3 min, 5 min, 2 min, 1 min
- Trend: consistent

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Live API calls over periodic sync: ensures analyser always uses latest KE-gene mappings
- Keep demo datasets local: demo data is static, does not change with Builder updates
- Critical: build CSV fallback before removing any static data files (rate limit + multi-worker cache risk)
- load_reference_sets() accepts ke_wp_df=None to avoid duplicating merge logic across CSV and API paths
- KE ID normalisation (space to colon) happens at API fetch boundary, not inside load_reference_sets()
- Empty BUILDER_API_URL raises ValueError (not network error) to distinguish config vs runtime failures
- [Phase 01-api-integration]: load_cached_reference_sets() returns (reference_sets, data_source) tuple so provenance flows to the template without a separate lookup
- [Phase 01-api-integration]: diskcache.FanoutCache at module level — all Gunicorn workers share same directory by design
- [Phase 01-api-integration]: Data source indicator uses neutral language: 'Local reference files' not error messaging for CSV fallback
- [Phase 02-ui-ux-and-documentation]: Background changed from sky blue (#93D5F6) to neutral light grey (#f5f6fa) — sky blue too dominant with branded dark nav
- [Phase 02-ui-ux-and-documentation]: Legacy header { display:none } kept rather than deleted to prevent bleed-through during incremental template migration
- [Phase 02-ui-ux-and-documentation]: .subheader reduced from red (#ff0000) to brand orange (#EB5B25) with reduced padding/font-size for subtle alert styling
- [Phase 02-ui-ux-and-documentation]: alert-band replaces .subheader h2 in index.html — more semantically appropriate for a page-level notice outside the card
- [Phase 02-ui-ux-and-documentation]: documentation.html scoped to three sections (Input File Format, Statistical Methods, Interpreting Results) — Getting Started, Understanding AOPs, Glossary, FAQ removed
- [Phase 02-ui-ux-and-documentation]: Volcano plot Plotly script wrapped in {% if volcano_data %} guard inside {% block scripts %} to prevent JS reference errors on initial page load
- [Phase 02-ui-ux-and-documentation]: CSV button uses btn--primary (blue) placed top-right of enrichment card header for immediate visibility without scrolling
- [Phase 02-ui-ux-and-documentation]: DataTables buttons bar preserved alongside branded CSV button for power user access to copy/excel/print
- [Phase 03-network-performance]: Cytoscape.js upgraded from 3.26.0 to 3.33.1; dagre rankDir changed from TB to LR for AOP pathway flow
- [Phase 03-network-performance]: Per-KE click-to-expand replaces bulk addGeneNodes() to avoid browser freeze on large networks (AOP:2 ~696 genes)
- [Phase 03-network-performance]: requestAnimationFrame before cy.layout() ensures spinner paints before CPU-intensive layout work begins
- [Phase 03-network-performance]: Shared gene nodes deduped by cy.getElementById check — single node with edges to all KEs shows cross-KE biological overlap
- [Phase 03-network-performance]: Zoom controls positioned absolute inside position:relative wrapper — no new wrapper needed; 1.25x zoom factor per step
- [Phase 03-network-performance]: downloadPNG() transparent by default (no bg param); report generation retains bg:'white' for PDF embedding

### Pending Todos

None.

### Blockers/Concerns

- Builder API bulk export endpoint (`/exports/gmt/ke-wp`) must be verified to exist on the live Builder instance before Phase 1 implementation begins; fallback is paginated `/api/v1/mappings` with aggressive caching
- Builder API rate limit scope (per-IP vs per-service) unclear; clarify before multi-worker production deployment
- Gunicorn worker count in Docker Compose needs confirmation to assess if in-memory cache failure (Pitfall 1) is already active in production

## Session Continuity

Last session: 2026-02-25
Stopped at: Completed 03-02-PLAN.md — Overlay zoom controls and transparent PNG export; Phase 3 complete
Resume file: None
