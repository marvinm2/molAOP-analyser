# Project Retrospective

*A living document updated after each milestone. Lessons feed forward into future planning.*

## Milestone: v1.0 — MVP

**Shipped:** 2026-02-27
**Phases:** 3 | **Plans:** 7 | **Commits:** 14

### What Was Built
- Builder API integration with paginated fetching, cross-process disk caching, and graceful CSV fallback
- VHP4Safety brand overhaul across all pages (sticky nav, Inter font, card layouts, brand palette)
- Documentation page (input format, statistical methods, result interpretation)
- Results page redesign with metadata summary card, styled enrichment table, CSV export
- Network performance optimisation (Cytoscape.js 3.33.1, dagre LR, per-KE progressive expansion)
- Overlay zoom controls and transparent PNG export

### What Worked
- Modular service architecture made API integration clean — new api_service.py slotted in without touching enrichment logic
- Shared base.html template eliminated boilerplate across 4 pages in a single plan
- Per-KE progressive gene expansion elegantly solved the large network freeze problem (AOP:2 ~696 genes)
- Three-tier data loading (API → cache → CSV) provides resilience without complexity

### What Was Inefficient
- Phase 1 plan checkboxes in ROADMAP.md not consistently updated (showed 1/2 when 2/2 were done)
- Pre-existing test failures carried through entire milestone without being fixed
- Some legacy CSS classes remain unused but weren't cleaned up

### Patterns Established
- diskcache.FanoutCache for cross-process caching in Gunicorn deployments
- `position: relative` wrapper pattern for overlay controls on Cytoscape containers
- Three-state data source indicator (api/cache/csv) with neutral fallback language
- `requestAnimationFrame` before CPU-intensive layout work to ensure spinner visibility

### Key Lessons
1. Per-KE progressive expansion is the right pattern for large AOP networks — bulk gene loading is never acceptable
2. Shared gene dedup via cy.getElementById prevents visual clutter and shows cross-KE biological overlap
3. Brand overhaul is most efficient when done with a shared base template first, then refactoring each page to extend it
4. CSV fallback should use neutral language ("Local reference files") not error-like messaging

### Cost Observations
- Model mix: primarily sonnet for execution, sonnet for planning/verification
- Execution velocity: 7 plans in ~22 minutes total (avg 3.1 min/plan)
- Notable: Phase 3 was fastest (1.5 min/plan avg) — template-only changes after Phase 2 established the structure

---

## Milestone: v3.0 — Cleanup, Gene-Set Export & UX Gaps

**Shipped:** 2026-05-11
**Phases:** 4 | **Plans:** 15 | **Commits:** 51 (6 days, 2026-05-06 → 2026-05-11)

### What Was Built
- AOP:472 fallback rows removed; AOP:472 now flows from SPARQL alongside the other dynamic kidney AOPs (Phase 10)
- `/demos` page restructure with recommended-AOP CTAs (PXR → AOP:1, cisplatin → kidney AOPs); main page leads with upload-your-data (Phase 10.1)
- Gene-by-KE CSV export (`KE_ID, KE_Title, Gene_Symbol, log2FC, pvalue, FDR`) respecting active significance thresholds (Phase 11, EXPO-04/05)
- Per-KE p-value/FDR embedded in Cytoscape network exports (Phase 11, EXPO-06, closes #50)
- Column detector raw/adjusted p-value split with preview-UI picker (Phase 11)
- Multipart CSV/TSV/TXT upload on the upload page with regression test (Phase 12, UPUX-04)
- `datetime.utcnow()` migration (14 sites) + duplicate `guess_id_type` consolidation + 8 debug `print()` → `logger.debug()` + 57-line dead `build_cytoscape_network_legacy()` deletion (Phase 12, DEBT-01/02/03/04)

### What Worked
- Tech debt was bundled into a single phase (12) with parallel wave execution, which finished in one day despite touching 4 unrelated areas
- Deferring EXPO-03 (GMT export) to the Builder was the right call — gene-by-KE CSV (EXPO-04) is strictly more interrogable than an intersection-style GMT
- Phase 10.1 was inserted between 10 and 11 mid-milestone and slotted in cleanly thanks to `--insert` phase tooling
- Wave-based plan ordering in Phase 12 (Wave 1 → 2 → 3) avoided REQUIREMENTS.md write conflicts on parallel autonomous executions
- The PXR demo grid being temporarily hidden in favour of cisplatin (then reversed) was easy to revert because demo cards live in a single `demos.html` block

### What Was Inefficient
- Phase 11-01 and 11-02 PLAN.md were closed via merge commits but never received written SUMMARY.md files — traceability table was the only completion record
- ROADMAP.md progress table got stale during the milestone (10.1 showed 2/3 In progress, 11 showed 0/4 Planned) — caught only at close
- AOPD-07/08/09 traceability rows said "Planned" though body checkboxes were ticked — divergence between two truth-sources of REQUIREMENTS.md
- Two 38-day-old todos (test failures, demo separation) sat in pending though v3.0 work had already resolved both — only triaged at close
- `.planning/debug/column-selection-reset.md` was status:resolved but the filename suggested otherwise; misled an initial progress scan

### Patterns Established
- **Deferred-to-builder** as an explicit requirement status alongside Planned / In Progress / Complete — captures upstream-of-truth decisions cleanly
- **Recommended-AOP mapping** in `Config.DEMO_AOP_RECOMMENDATIONS` — pairs demos with curated AOP sets via a single declarative dict
- **Wave-based plan execution** for autonomous phases with shared-file conflicts — explicit Wave 1/2/3 grouping in PLAN.md
- **Pre-close traceability sweep** as a milestone-close ritual — catch divergence between body checkboxes and traceability table before archiving

### Key Lessons
1. The molAOP Builder is the upstream source of truth for KE→gene libraries — duplicating that data in the analyser (as GMT) would create drift without value
2. Long-format gene-by-KE CSV beats intersection-style GMT when the analyser-specific need is "which genes from this analysis drove each KE"
3. Tech-debt bundles ship fast when the items are genuinely orthogonal (datetime, helpers dedup, prints, dead code) — wave scheduling handles the few shared-file collisions
4. Insert-phase tooling (`/gsd-phase --insert 10.1`) is the right way to handle "we need to do X before Y" — preserves numbering history
5. Run `audit-open` and a traceability sweep BEFORE drafting milestone close — both surface drift that's invisible from progress reports alone

### Cost Observations
- Model mix: opus for planning/verification, sonnet for execution (per `.planning/config.json`)
- Most efficient phase: Phase 12 — 5 plans across 4 tech-debt items + UPUX-04 closed in one working day via wave parallelism
- Most documentation rework: Phase 11 — EXPO-03 deferral cascaded through ROADMAP, REQUIREMENTS, MILESTONES, and a documentation pointer (Plan 11-04 dedicated to doc cleanup alone)

---

## Cross-Milestone Trends

### Process Evolution

| Milestone | Commits | Phases | Key Change |
|-----------|---------|--------|------------|
| v1.0 | 14 | 3 | Initial milestone — established modular service architecture and brand foundation |
| v3.0 | 51 | 4 | Tech-debt bundling with wave scheduling; explicit `Deferred-to-Builder` requirement state |

### Cumulative Quality

| Milestone | Tests Passing | Pre-existing Failures | New Failures |
|-----------|--------------|----------------------|--------------|
| v1.0 | 23/28 | 5 | 0 |
| v3.0 | 146/146 | 0 | 0 |

### Top Lessons (Verified Across Milestones)

1. Progressive rendering patterns (show KEs first, expand genes on demand) prevent browser freeze on large datasets
2. Three-tier data loading (API → cache → fallback) provides resilience without user-visible degradation
3. The Builder is the source of truth for KE→gene libraries — defer GMT-style exports upstream rather than duplicate them
4. Pre-close audit (`audit-open` + traceability sweep) catches doc drift that progress reports miss
