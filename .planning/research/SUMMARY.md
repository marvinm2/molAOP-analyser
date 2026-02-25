# Project Research Summary

**Project:** molAOP Analyser — Enhancement Round
**Domain:** Flask bioinformatics web application — AOP enrichment analysis with live API integration, UI modernisation, and network visualisation
**Researched:** 2026-02-25
**Confidence:** HIGH (existing codebase inspected directly; stack versions verified via PyPI/official docs; architecture patterns confirmed against official sources)

## Executive Summary

The molAOP Analyser is a working Flask/Python enrichment analysis tool requiring three targeted enhancements: live API integration with the molAOP Builder service, a VHP4Safety brand-compliant UI overhaul, and Cytoscape.js network performance improvements. Unlike a greenfield build, the core analysis pipeline (Fisher's exact test ORA, FDR correction, Cytoscape.js network, SQLite persistence) is already production-ready and must not be disrupted. The recommended approach keeps the existing service architecture entirely intact and introduces change at three narrow, independent seams: the reference data loading path (helpers.py + new api_client.py), the template layer (new base.html template inheritance shell), and the Cytoscape.js initialisation code (progressive rendering, performance flags). All three seams are architecturally isolated, enabling parallel development on separate branches.

The most important recommendation across all research: build the Builder API fallback to local CSV files before removing any static data files. The Builder API has a strict 100-requests/hour rate limit that is easily exhausted under multi-worker deployment, and the in-memory CacheManager used today is per-process, meaning each Gunicorn worker fetches its own copy independently. These two issues in combination — rate limit exhaustion and per-process cache isolation — represent the highest-probability production failure mode for this milestone. Both must be resolved in Phase 1 before any live API calls go to production.

The competitive landscape confirms that no existing tool combines AOP topology visualisation with transcriptomic expression overlay and pre-curated regulatory toxicology case studies — this combination is the tool's unique positioning within the VHP4Safety ecosystem. The UI overhaul and documentation page are critical for platform trust and adoption, not just aesthetics. The recommended UI stack (Tailwind CSS v4, htmx, Alpine.js) integrates naturally with Flask's server-rendered Jinja2 model without requiring a JavaScript build pipeline at runtime.

## Key Findings

### Recommended Stack

The new stack additions cover three distinct concerns, each with clear rationale. For API integration: `httpx` (async-ready HTTP client), `tenacity` (retry with exponential backoff and jitter), and `flask-caching` (shared cache with Redis backend for multi-worker safety). The existing `requests`-based patterns can be replaced with httpx at minimal migration cost. For UI: Tailwind CSS v4 maps VHP4Safety brand tokens in a single `@theme` block, htmx enables partial page updates without a JavaScript SPA framework, and Alpine.js handles residual interactive state. Both htmx and Alpine.js load from CDN — no Node build pipeline required at runtime. For network performance: `cytoscape-fcose` (force-directed layout, supersedes cose-bilkent), `cytoscape-dagre` (DAG-based hierarchical layout for AOP skeletons), and Cytoscape.js's experimental WebGL renderer (available in 3.33.1; delivers 5x FPS improvement but requires a canvas fallback).

**Core technologies:**
- `httpx==0.28.1`: HTTP client for molAOP Builder API calls — async-ready, nearly identical API to requests, type-annotated
- `tenacity==9.1.4`: Retry with exponential backoff for API calls — cleaner decorator API than urllib3 Retry, supports jitter and conditional retry
- `flask-caching==2.3.1`: Cache API responses in memory (SimpleCache for dev, Redis for production) — 5-minute TTL eliminates redundant Builder API calls
- `jinja2-fragments==1.11.0`: Render individual Jinja2 template blocks — required for htmx partial rendering pattern
- `Tailwind CSS v4`: Brand token mapping via CSS-first `@theme` directive — maps all VHP4Safety colours and Inter font as CSS variables
- `htmx 2.0.8`: Partial template updates without full page reload — pairs naturally with Jinja2 SSR, zero build step
- `Alpine.js 3.15.8`: Lightweight interactive UI state (dropdowns, modals, tooltips) — 15KB, declarative `x-data`, zero build step
- `cytoscape-fcose==2.2.0`: Force-directed layout for biological networks — official Cytoscape.js recommendation, spectral initialisation
- `cytoscape-dagre==2.5.0`: DAG-based hierarchical layout — preserves AOP pathway directionality for KE skeleton
- Cytoscape.js WebGL renderer (built-in, 3.31+): GPU-accelerated rendering — 1,200-node networks at 100+ FPS vs 20 FPS on canvas

**Version requirements:** tenacity 9.1.4 requires Python 3.10+; verify Docker base image version. All other additions are compatible with Python 3.8+.

### Expected Features

This is an enhancement round on a working tool. The four milestone deliverables are all confirmed P1 based on competitive analysis and user expectation research.

**Must have (this milestone):**
- Live Builder API integration — KE-gene mappings always current with AOP-Wiki state; no manual data refresh needed; primary differentiator for the VHP4Safety platform
- VHP4Safety brand UI overhaul — Inter typography, full secondary palette, sticky navigation; platform trust and adoption depend on visual continuity
- Network performance improvements — profile largest dataset (kidney AOP, cisplatin 72hr 50μM) first; apply progressive rendering, element batching, haystack edge style
- Documentation page — input format specification, statistical method explanation, result interpretation guide; table stakes for any new user arriving from the platform

**Should have (add after milestone P1 is stable):**
- Enrichment table CSV export — confirm if already present; trivial to add if not; first user complaint triggers this
- Network PNG/SVG export button — Cytoscape.js has built-in export; surface it with a clear UI button
- Shareable result URLs — SQLite persistence exists; expose a `/results/<experiment_id>` route; trigger: community request

**Defer (v2+):**
- Comparative dataset view — side-by-side AOP networks for different cisplatin timepoints; high user value, high complexity
- Structured KE metadata tooltips from AOP-Wiki API — inline KE descriptions on node hover; currently links to external AOP-Wiki
- AI/LLM analysis panel — deliberately avoided; AOP Net preprint (Aug 2025) shows the pattern but it remains unvalidated, high-maintenance, and high-hallucination-risk for regulatory toxicology

**Anti-features (do not build):**
- GSEA alongside ORA — doubles statistical implementation surface; ORA with Fisher's exact is statistically appropriate for binary significant/non-significant gene lists
- Multi-omics input — no clear KE-gene mappings exist for non-transcriptomic identifiers at this maturity level
- OAuth/user authentication — VHP4Safety platform authentication is out of scope; SQLite persistence covers reproducibility
- Full pathway database browsing (GO/KEGG/Reactome) — dilutes AOP focus; direct users to g:Profiler or Enrichr for general pathway analysis

### Architecture Approach

The architecture research confirmed that the correct strategy is surgical modification at three independent seams rather than any broad refactoring. The existing service layer (data_service, enrichment_service, network_service, column_detector, report_service, gene_id_validator) remains entirely unchanged. The critical architectural invariant is the `{KE_ID: set(gene_symbols)}` contract returned by `helpers.py` — this dict shape must not change because all downstream services depend on it. The new `services/api_client.py` (`BuilderAPIClient` class) is introduced as an independently testable HTTP client that `helpers.py` calls first, falling back to the existing CSV merge logic on any `BuilderAPIError`. Templates are restructured around a new `base.html` using Jinja2 inheritance, but all `render_template()` call signatures in `app.py` remain identical.

**Major components:**
1. `services/api_client.py` (NEW) — BuilderAPIClient with httpx/tenacity retry, timeout config, raises BuilderAPIError; independently testable with `pytest-httpx`
2. `helpers.py` (MODIFIED) — API-first reference set loading; optional `api_client` parameter; falls back to existing CSV merge logic on failure; contract unchanged
3. `templates/base.html` (NEW) — brand base layout with Inter font, VHP4Safety CSS custom properties, sticky navigation, Jinja2 block slots
4. `config.py` (EXTENDED) — adds `BUILDER_API_URL`, `BUILDER_API_TIMEOUT`, confidence level filter as explicit config constants
5. Cytoscape.js init in `results.html` (MODIFIED) — progressive gene layer (KE skeleton first, genes on tap), `cy.batch()` wrapping, haystack edges for KE-gene, WebGL renderer with canvas fallback

**Build order dependency:** Phase A (API client) and Phase B (UI overhaul) are independent and can be developed in parallel. Phase C (network performance) depends on Phase B because it modifies `results.html`, which needs the new base.html shell first.

### Critical Pitfalls

1. **In-memory cache does not survive multi-worker deployment** — The existing `CacheManager` is per-process. Under Gunicorn with multiple workers, each worker fetches its own API copy, exhausting the Builder's 100-requests/hour rate limit within minutes of startup. Fix before any live API call goes to production: replace `CacheManager` with Redis-backed `flask-caching`, or write the fetched payload to the shared SQLite database with a timestamp check.

2. **Rate limit exhaustion during full reference set fetch** — Building the full `KE → gene set` dict via paginated `/api/v1/mappings` may require 20-50+ requests per cold cache cycle. Use the bulk export endpoint (`/exports/gmt/ke-wp?min_confidence=High`) as the primary data source — one request returns the complete dataset. Reserve paginated endpoints for targeted AOP structure queries only.

3. **Startup failure when Builder API is unreachable** — `app.py` currently loads reference sets at module import time. Replacing local CSV reads with a network call makes the Builder a hard startup dependency. Implement stale-on-error cache semantics: serve any cached payload regardless of age on API failure, fall back to static CSVs if no cache exists, and move the initial fetch to a background thread so the application starts immediately.

4. **KE-gene data schema mismatch between CSV and API** — API returns `ke_id` as `"KE 55"` (string with space); existing CSV normalisation uses numeric `KE_ID`. A silent schema mismatch produces an empty `reference_sets` dict — enrichment analysis runs but finds no gene sets. Write an explicit schema validation step and a comparison test (CSV-derived vs API-derived reference sets for known AOPs) before deleting any CSV files.

5. **Dagre layout freeze on gene-expanded networks** — Running Dagre over the full expanded graph (KE nodes + hundreds of gene nodes) causes browser tab freezes. Never call a layout algorithm on the mixed KE+gene graph. Pre-position gene nodes programmatically around their parent KE anchor using `cy.position()`, wrapped in `cy.batch()`. Reserve Dagre layout for the KE-only skeleton only.

## Implications for Roadmap

Based on the combined research, the three milestone deliverables map naturally to three sequential phases. The dependency between Phase B and Phase C is the only hard ordering constraint. Phases A and B are parallelisable.

### Phase 1: API Integration Layer

**Rationale:** The API integration is the highest-value deliverable (always-current reference data) and has the most critical pre-conditions that must be in place before anything else. The shared cache problem (Pitfall 1), rate limit strategy (Pitfall 2), and startup resilience (Pitfall 3) must all be solved here. Building this first also ensures the CSV fallback remains intact and tested throughout the rest of the project.

**Delivers:** Live KE-gene reference data from molAOP Builder; graceful degradation to local CSVs when Builder is unavailable; rate-limit-safe cache strategy for multi-worker production.

**Addresses:** Live Builder API integration (FEATURES.md P1), stable current reference data (table stakes), Builder service availability fallback (dependency noted in FEATURES.md).

**Avoids:** In-memory cache multi-worker failure (Pitfall 1), rate limit exhaustion (Pitfall 2), startup failure on API unavailability (Pitfall 3), schema mismatch producing empty reference sets (Pitfall 4), API confidence level filter silently changing enrichment results (Pitfall 9).

**Stack additions:** `httpx`, `tenacity`, `flask-caching`, `flask-cors`, `pytest-httpx` (test support).

**Research flag:** STANDARD PATTERNS — API client with retry/fallback is well-documented. Builder API endpoints and response shape need verification against live Builder instance before writing transformation code; this is implementation validation, not research.

---

### Phase 2: UI/UX Overhaul and Documentation

**Rationale:** The UI overhaul is independent of Phase 1 and can be developed in parallel on a separate branch. However, the documentation page must use the same nav/header component established in this phase, so both belong together. CSS custom properties must be defined first (as the single source of truth for brand colours in both CSS and Cytoscape.js styles) before applying them to any template or network style.

**Delivers:** VHP4Safety brand-compliant UI across all pages; Inter typography; sticky navigation matching molAOP Builder style; documentation page covering input format, statistical method, and result interpretation.

**Addresses:** UI/UX overhaul (FEATURES.md P1), consistent platform branding (table stakes), documentation page (table stakes), responsive layout (table stakes).

**Avoids:** CSS brand overhaul breaking Cytoscape.js canvas sizing (Pitfall 7 — test `cy.width() > 0` after every layout change), brand colours inconsistent between CSS and Cytoscape styles (Pitfall 8 — use CSS custom properties as single source of truth), in-template JSON blobs growing too large (Pitfall 6 — refactor large payloads to API endpoints during this page restructure).

**Stack additions:** Tailwind CSS v4 (build-time, CLI), htmx 2.0.8 (CDN), Alpine.js 3.15.8 (CDN), `jinja2-fragments`, Inter font (self-hosted for GDPR compliance in EU deployments).

**Build order within phase:** CSS custom properties first → `base.html` → `index.html` → `results.html` → `documentation.html`. Visual regression check at each step against unchanged Flask routes.

**Research flag:** STANDARD PATTERNS — Tailwind v4, htmx, Alpine.js, and Jinja2 template inheritance are all well-documented. No deeper research needed.

---

### Phase 3: Network Visualisation Performance

**Rationale:** Network performance improvements require the new `results.html` shell from Phase 2, so this phase cannot proceed in parallel. Must be preceded by profiling actual large AOP datasets (kidney organ network, cisplatin 72hr 50μM) before implementing solutions — the bottleneck may be layout algorithm, DOM element count, or Cytoscape.js render loop. Avoid premature optimisation.

**Delivers:** Progressive gene layer rendering (KE skeleton on load, genes on click); browser-freeze-free gene node expansion; improved network aesthetics; performance verified against largest available test datasets.

**Addresses:** Network performance at AOP scale (FEATURES.md P1), network interactivity improvements (FEATURES.md launch criteria), network node/edge count within visible range (table stakes).

**Avoids:** Dagre layout freeze on gene-expanded networks (Pitfall 5 — pre-position genes programmatically, never run layout on full expanded graph), `cy.add()` without `cy.batch()` causing incremental redraws (Performance Traps table), function-based Cytoscape style values causing per-render recalculation (Architecture Anti-Pattern 4).

**Stack additions:** `cytoscape-fcose==2.2.0` (CDN), `cytoscape-dagre==2.5.0` (CDN), Cytoscape.js WebGL renderer (built-in, use with canvas fallback via try/catch).

**Research flag:** NEEDS VALIDATION — Profile actual large datasets before implementing solutions. WebGL renderer is experimental in 3.33.1; verify browser compatibility in target research environments before enabling by default.

---

### Phase 4: Post-Milestone Polish (v1.x)

**Rationale:** These features are confirmed desirable but not required for the milestone. Add after all P1 items are stable and verified in production.

**Delivers:** Enrichment table CSV download, network PNG/SVG export button, shareable result URLs via `/results/<experiment_id>` route.

**Addresses:** Export results as table (table stakes — confirm if already implemented), export network as image (table stakes), stable/shareable URLs (FEATURES.md P2).

**Research flag:** STANDARD PATTERNS — CSV export is trivial Flask; Cytoscape.js PNG export is built-in; SQLite experiment record already exists for URL routing. No research phase needed.

---

### Phase Ordering Rationale

- Phase 1 must not remove CSV fallback files until the API-derived `reference_sets` dict produces identical results to the CSV-derived version on at least two known AOP case studies (liver steatosis, kidney failure). Keep CSVs in the repository throughout the project.
- Phases 1 and 2 are independent; they can be developed on separate branches and merged in either order, but Phase 1's cache strategy must be confirmed before deployment.
- Phase 3 depends on Phase 2's `results.html` template restructure. Do not start Phase 3 until Phase 2 templates are merged.
- The three pitfalls with the highest recovery cost (Pitfall 1 in-memory cache, Pitfall 4 schema mismatch, Pitfall 3 startup failure) must be resolved in Phase 1 before any code reaches production.

### Research Flags

Phases needing validation during implementation:
- **Phase 1:** Verify actual Builder API response field names and pagination behavior against a live Builder instance before writing transformation code. The API reference doc (KE-MAPPING-API-REFERENCE.md) exists internally but field name conventions must be confirmed empirically (Pitfall 4 risk).
- **Phase 3:** Profile largest available test datasets (kidney organ network with cisplatin data) before writing any performance code. Identify the actual bottleneck — do not assume layout algorithm is the cause. WebGL renderer needs browser compatibility testing in target research environments.

Phases with standard patterns (no research phase needed):
- **Phase 2:** All UI technologies (Tailwind v4, htmx, Alpine.js, Jinja2 template inheritance) are well-documented with official sources. CSS custom properties pattern for Cytoscape style consistency is a single-file change.
- **Phase 4:** CSV export, Cytoscape.js PNG export, and Flask URL routing are standard patterns with no ambiguity.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | All versions verified via PyPI and official docs; rationale confirmed from multiple sources including official blogs |
| Features | MEDIUM-HIGH | Core milestones well-defined; competitor analysis from peer-reviewed papers where available; AOP Net (Aug 2025) is preprint only |
| Architecture | HIGH | Existing codebase inspected directly; patterns verified against official Flask and Cytoscape.js documentation; build order dependencies are deterministic |
| Pitfalls | MEDIUM-HIGH | Critical pitfalls verified across multiple consistent sources; UX pitfalls backed by lower-confidence WebSearch results but plausible and low-cost to prevent |

**Overall confidence:** HIGH

### Gaps to Address

- **Builder API bulk export endpoint availability:** Research identified `/exports/gmt/ke-wp` as the correct bulk data source (avoids rate limit exhaustion), but this endpoint must be verified to exist and return data in the expected format on the actual live Builder instance before Phase 1 implementation begins. If the endpoint is absent, the fallback strategy is paginated `/api/v1/mappings` with aggressive caching.
- **Builder API rate limit scope:** The 100-requests/hour limit is documented; it is unclear whether this is per-IP, per-service, or per-authenticated-client. Clarify with the Builder team before deploying multi-worker production configuration.
- **Gunicorn worker count in Docker Compose:** The current Docker Compose configuration may already run multiple workers. Confirm the worker count to assess whether the in-memory cache failure (Pitfall 1) is already active in the current production deployment.
- **WebGL renderer browser coverage:** The Cytoscape.js WebGL renderer is experimental. Research confirms the performance gains are significant, but browser compatibility in institutional research network environments (common in regulatory toxicology) needs empirical verification before enabling by default.
- **Inter font GDPR compliance:** Self-hosting Inter is recommended for EU deployments to avoid sending visitor IPs to Google Fonts CDN. Confirm the deployment region and institutional data handling policy before choosing CDN vs self-hosted.

## Sources

### Primary (HIGH confidence)
- PyPI — httpx 0.28.1, tenacity 9.1.4, Flask-Caching 2.3.1, jinja2-fragments 1.11.0, flask-cors 6.0.2 — versions verified
- Tailwind CSS official blog — v4.0 release January 22 2025, CSS-first config [https://tailwindcss.com/blog/tailwindcss-v4]
- Cytoscape.js official site — version 3.33.1, performance options [https://js.cytoscape.org/]
- Cytoscape.js blog — WebGL renderer preview, performance benchmarks [https://blog.js.cytoscape.org/2025/01/13/webgl-preview/]
- Cytoscape.js GitHub — performance.md official documentation
- Flask official docs — template inheritance, service layer patterns
- Alpine.js GitHub releases — v3.15.8 confirmed
- htmx.org — version 2.0.8 stable release notes
- Existing codebase — app.py, helpers.py, services/, config.py, cache_manager.py, database.py, exceptions.py inspected directly
- EnrichmentMap:RNASeq, Bioinformatics Advances 2025 — PMC peer-reviewed [https://pmc.ncbi.nlm.nih.gov/articles/PMC12373637/]
- g:Profiler 2023, Nucleic Acids Research — Oxford peer-reviewed
- WikiPathways 2024, Nucleic Acids Research — Oxford peer-reviewed

### Secondary (MEDIUM confidence)
- AOP-networkFinder, Bioinformatics Advances 2025 — Oxford peer-reviewed [https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf007/7972740]
- AOP-helpFinder 3.0, Bioinformatics Oxford 2025 — peer-reviewed [https://academic.oup.com/bioinformatics/article/41/7/btaf381/8177148]
- VHP4Safety platform tools — official platform [https://platform.vhp4safety.nl/tools]
- Flask multi-worker cache isolation — Medium / community sources (multiple consistent)
- htmx + Alpine.js production Flask usage — Hacker News + teclado.com case studies
- cytoscape-fcose npm package — npm verified; maintenance status noted as stable but not actively developed
- cytoscape-dagre GitHub v2.5.0 — version verified; stable

### Tertiary (LOW confidence)
- AOP Net preprint, bioRxiv August 2025 — not peer-reviewed [https://www.biorxiv.org/content/10.1101/2025.08.25.672239v1.full]; used only to inform anti-feature decision on LLM panel
- Network visualisation UX best practices — PMC survey [https://pmc.ncbi.nlm.nih.gov/articles/PMC5540468/]
- UX pitfalls (loading state, font readability) — WebSearch results; consistent with general UX practice but not domain-specific

---
*Research completed: 2026-02-25*
*Ready for roadmap: yes*
