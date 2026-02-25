# Roadmap: molAOP Analyser

## Overview

Three phases transform the working molAOP Analyser into a platform-integrated, brand-consistent tool: first making reference data live from the Builder API (with safe fallback), then overhauling the UI to match the VHP4Safety platform aesthetic and adding documentation, then optimising network visualisation performance for large AOP datasets. Phase 1 and Phase 2 are independent; Phase 3 depends on Phase 2 completing the results.html template restructure.

## Phases

- [x] **Phase 1: API Integration** - Replace static CSV reference data with live molAOP Builder API calls, with graceful CSV fallback and rate-limit-safe caching (completed 2026-02-25)
- [x] **Phase 2: UI/UX and Documentation** - Overhaul all pages to VHP4Safety brand standards, add sticky navigation, and publish a documentation page (completed 2026-02-25)
- [x] **Phase 3: Network Performance** - Optimise Cytoscape.js for large AOP networks and surface export controls (1/2 plans complete) (completed 2026-02-25)

## Phase Details

### Phase 1: API Integration
**Goal**: The analyser fetches current KE-gene mappings from the molAOP Builder API, degrades gracefully when the Builder is unreachable, and caches responses safely across all worker processes
**Depends on**: Nothing (first phase)
**Requirements**: API-01, API-02, API-03
**Success Criteria** (what must be TRUE):
  1. Running an analysis against a live Builder instance uses API-sourced gene sets, not local CSV files
  2. With the Builder API unreachable, analysis completes normally using local CSV files and no error is shown to the user
  3. Reloading the application under a multi-worker Gunicorn deployment does not trigger repeated Builder API calls within the cache TTL window
  4. The enrichment results produced from API-sourced data match the results produced from CSV-sourced data for at least two known AOP case studies (liver steatosis, kidney failure)
**Plans**: 2 plans
Plans:
- [ ] 01-01-PLAN.md -- API service, helpers refactor, config, and diskcache dependency
- [ ] 01-02-PLAN.md -- Wire API-first loading into app, disk cache integration, results page data source indicator

### Phase 2: UI/UX and Documentation
**Goal**: Every page uses the VHP4Safety brand palette and Inter typography, a sticky navigation header matches the molAOP Builder style, and a documentation page explains how to use the tool
**Depends on**: Phase 1
**Requirements**: UI-01, UI-02, UI-03, UI-04, DOC-01, EXP-01
**Success Criteria** (what must be TRUE):
  1. All pages render Inter font and apply VHP4Safety primary and secondary brand colours with no legacy styling visible
  2. A sticky navigation header is present on every page, displays app branding, and matches the molAOP Builder header pattern
  3. The results page presents enrichment data in clearly separated cards with consistent spacing and typography
  4. The documentation page is reachable from the navigation and covers input file format, statistical method, and result interpretation
  5. A user can download the enrichment results table as a CSV file from the results page
**Plans**: 3 plans
Plans:
- [x] 02-01-PLAN.md -- Base template (sticky nav, Inter font, footer), global CSS overhaul, About page and route
- [x] 02-02-PLAN.md -- Refactor index.html (single centred card) and documentation.html (3 required sections)
- [x] 02-03-PLAN.md -- Refactor results.html (metadata card, styled table card, CSV export button, analysis_date)

### Phase 3: Network Performance
**Goal**: Large AOP networks render without browser slowdown, and users can export the network as an image for publication use
**Depends on**: Phase 2
**Requirements**: NET-01, NET-02
**Success Criteria** (what must be TRUE):
  1. The kidney organ AOP network loaded with the cisplatin 72hr 50uM dataset renders to an interactive state in under three seconds on a modern laptop browser
  2. Expanding gene nodes on a KE in the largest available network does not freeze the browser tab
  3. A PNG export button is visible on the results page and downloads a full-resolution screenshot of the current network view
**Plans**: 2 plans
Plans:
- [x] 03-01-PLAN.md -- Spinner overlay, Cytoscape.js upgrade, dagre LR layout, per-KE click-to-expand, collapse all
- [ ] 03-02-PLAN.md -- Overlay zoom controls, transparent PNG export

## Progress

**Execution Order:**
Phases execute in numeric order: 1 → 2 → 3

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. API Integration | 1/2 | Complete    | 2026-02-25 |
| 2. UI/UX and Documentation | 3/3 | Complete   | 2026-02-25 |
| 3. Network Performance | 2/2 | Complete    | 2026-02-25 |
