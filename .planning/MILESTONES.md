# Milestones

## v3.0 Cleanup, Gene-Set Export & UX Gaps (Shipped: 2026-05-11)

**Phases completed:** 4 phases, 15 plans
**Timeline:** 2026-05-06 → 2026-05-11 (6 days)
**Git range:** 51 commits, 48 files, +2400 / -426 lines
**Requirements:** 12/12 v3.0 requirements complete (EXPO-03 deferred-to-builder by design)

**Key accomplishments:**

- **AOP:472 hardcoded fallback removed** — kidney AOP:472 now flows from SPARQL alongside the other dynamic kidney AOPs; `KIDNEY_AOP_IDS` curation criterion documented in `config.py` (Phase 10, AOPD-07/08/09)
- **Demos page restructure** — new `/demos` route with curated PXR + cisplatin cards and recommended-AOP CTAs (AOP:1 for PXR, kidney AOPs for cisplatin); main page leads with "upload your data" (Phase 10.1)
- **Gene-by-KE CSV export** — new export with columns `KE_ID, KE_Title, Gene_Symbol, log2FC, pvalue, FDR` respecting active significance thresholds (Phase 11, EXPO-04/05)
- **Cytoscape network export carries significance** — per-KE p-value and FDR embedded in network exports; closes issue #50 (Phase 11, EXPO-06)
- **Column detector raw/adjusted p-value split** — separate detection patterns surfaced through the preview UI picker (Phase 11)
- **GMT export deferred-to-Builder** — analyser does not duplicate the upstream Builder's gene-set library; analyser-specific need satisfied by EXPO-04 gene-by-KE long-format CSV (Phase 11, EXPO-03)
- **Real local file upload** — multipart CSV/TSV/TXT upload on the main upload page with regression test (Phase 12, UPUX-04)
- **Tech debt swept** — 14 `datetime.utcnow()` call sites migrated to timezone-aware constructors; duplicate `guess_id_type` consolidated in `services/data_service.py`; 8 debug prints in `helpers.py` migrated to `logger.debug()`; 57-line dead `build_cytoscape_network_legacy()` deleted (Phase 12, DEBT-01/02/03/04)
- **Test suite stable** — all 146 tests pass; new multipart upload regression test added

**Known doc-hygiene debt:** Phase 11-01 and 11-02 PLAN.md files exist on disk but lack written SUMMARY.md (work is verified complete via traceability table and live in production: column detector split + cytoscape significance metadata).

---

## v2.0 Batch Analysis & Dynamic AOPs (Shipped: 2026-03-03)

**Phases completed:** 6 phases, 15 plans, 4 tasks

**Key accomplishments:**

- (none recorded)

---

## v1.0 MVP (Shipped: 2026-02-27)

**Phases completed:** 3 phases, 7 plans, 0 tasks

**Key accomplishments:**

- Live API integration with molAOP Builder for KE-gene mappings, with cross-process disk caching and graceful CSV fallback
- VHP4Safety brand overhaul — shared base template with sticky nav, Inter typography, and card-based layouts
- Documentation page with input format, statistical methods, and result interpretation
- Results page redesign with metadata summary card, styled enrichment table, and CSV export
- Cytoscape.js 3.33.1 with dagre LR layout, per-KE progressive gene expansion, and loading spinners
- Overlay zoom controls and transparent PNG export for publication use

---
