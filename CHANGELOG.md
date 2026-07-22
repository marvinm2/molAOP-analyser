# Changelog

All notable changes to the Molecular AOP Analyser are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions correspond to development milestones; the deployed service tracks `main`
via `ghcr.io/marvinm2/molaop-analyser`.

## [Unreleased]

### Added

- **Batch analysis can now run GSEA, not only Fisher/ORA** (#76). Batch mode existed to
  compare conditions — compound, dose, timepoint — yet it was hardwired to the one method
  that is worst suited to that job: over-representation asks whether a thresholded gene list
  is unusually full of a Key Event's genes, and answers per condition in isolation. GSEA asks
  whether a Key Event's genes shift coherently in the ranking, which is exactly the quantity
  worth reading down a dose series. On a nine-condition platinum/AOP-472 dataset the two
  disagree informatively: KE 1194 (DNA damage) reaches significance in 2 of 9 conditions by
  ORA and 7 of 9 by GSEA, and KE 177 (mitochondrial dysfunction) is called *depleted* by ORA
  while GSEA shows coordinated positive enrichment — small per-gene effects with a consistent
  upward bias, which a threshold discards by construction. Getting that table before meant
  running every condition as a separate single analysis and assembling NES and FDR by hand,
  which also threw away the harmonised background that makes batch results comparable at all.

  The batch form gains the same method selector as the single-analysis form (selecting GSEA
  hides the thresholds, which it ignores), the method is stored on the batch rather than per
  condition so every condition is scored identically, and the comparison view switches to
  NES + FDR: the heatmap colours on the signed NES over a diverging scale centred at zero, so
  a coordinated up-shift and a down-shift are no longer both rendered as "significant", and
  the matrix keeps a sub-threshold NES rather than blanking it, because a consistent weak
  signal across a dose series is itself the observation. `batches.method` is added by the
  established idempotent PRAGMA-then-ALTER migration, and a NULL — every batch created before
  this change — reads back as `ora`, so existing batches still open, run, export and compare
  exactly as they did.

- **The AOP picker can now be populated from the molAOP Builder** (#3). The Builder gained
  `GET /api/v1/aops` (builder #207), which returns whole AOPs — title, KE count, and mapping
  coverage split across WikiPathways / GO / Reactome — rather than the bare KE IDs that were
  all `ke_aop_context` could offer. `services/aop_discovery_service.fetch_aops_from_builder`
  consumes it as a new third tier, so a SPARQL outage now leaves a picker holding **274
  curated AOPs with titles** instead of the five entries in `Config.CASE_STUDY_AOPS`. Cached
  for one hour rather than the usual week, since that list omits every unmapped AOP and
  should not outlive the outage that produced it.

  The live SPARQL cross-reference stays primary: the Builder's AOP membership comes from a
  precomputed snapshot, which had gone four months without a refresh and was undercounting
  (AOP 625: 15 KEs of 18) until builder #207 regenerated it.

- Repository tooling brought in line with the molAOP Builder: CI workflow (test matrix,
  flake8, coverage, startup smoke test), code-quality workflow (ruff/black/isort, bandit,
  pip-audit), Dependabot, CODEOWNERS, issue and pull-request templates, `LICENSE` (GPL-2.0),
  `pyproject.toml` and this changelog.
- Design-token layer in `static/css/main.css` mirroring the Builder's `:root` block —
  VHP4Safety palette, spacing, radius, shadow and z-index scales as CSS custom properties.

### Fixed

- **WikiPathways gene membership no longer comes only from a bundled snapshot** (#79).
  KE→pathway mappings were fetched live from the Builder, but pathway→gene membership was
  resolved against `data/edges_wpid_to_gene.csv` — a snapshot topping out at `WP5452`. Any
  pathway curated since resolved to **no genes at all**, silently, because the merge is an
  inner join. 8 of the 79 live-mapped pathways were affected; 8 Key Events lost every gene
  they had and were reported as *"no gene set mapped"*, which is the opposite of what
  happened — they are curated, and the tool could not resolve them. It also quietly deflated
  the WikiPathways coverage that a resource comparison reports, while the run advertised its
  provenance as `WikiPathways (Builder API, live)`: the *mappings* were live, the gene
  membership behind them was not.

  Membership now comes from the Builder's own KE-WP GMT export, per pathway, with the
  bundled CSV still covering anything the Builder does not answer for. Concretely: KE 1115
  (Reactive oxygen species) goes from excluded-entirely to 49 genes, KE 1392 from 33 to 73,
  KE 177 from 109 to 115, and the reference sets grow from 87 to 90 Key Events.

  Three pathways (`WP1234`, `WP3980`, `WP4010`) resolve in neither source — `WP1234` is a
  404 upstream. Those are now **named in a warning** on the results page and in the run's
  provenance rather than disappearing, so a Key Event that looks uncovered can be told apart
  from one that genuinely is.

### Changed

- Page-specific CSS moved out of the templates into `static/css/pages/`.
- `static/css/style.css` renamed to `static/css/main.css` for parity with the Builder.
- `docs/KE-MAPPING-API-REFERENCE.md` moved from the repository root.

### Fixed

- `pytest.ini` used a `[tool:pytest]` section header, which is only valid inside
  `setup.cfg`, so the entire test configuration — markers, testpaths, coverage settings,
  warning filters — was silently ignored. Coverage measurement and the marker registry now
  actually apply, and the coverage denominator excludes the test files themselves.
- `.planning/` is no longer tracked in git (it was already listed in `.gitignore`).

## [4.0] — 2026-05-18 — Onboarding, GSEA & Hub Genes

### Added

- **Guided tour** — six-step interactive walkthrough (upload → column detection → threshold
  UI → AOP picker → results table → Cytoscape network), reachable from the landing-page
  banner, `/demos` and `/documentation`. Skip / don't-show-again persists in `localStorage`
  (#18).
- **GSEA** as an alternative to Fisher's exact test — `gseapy.prerank` ranked enrichment with
  a method radio on the analysis form, NES/FDR results table, NES-scaled network colouring
  and method-aware reports (#52).
- **Hub genes panel** — sortable table ranking genes by Key Event participation within the
  selected AOP, with significant-only filter, panel→network highlight, ≥3-KE hub badges,
  CSV export and an opt-in Cytoscape gene-node layer (#53).
- **Column auto-detection clarity** — per-column confidence badges with percentages and
  detector-reason lists on both the single-analysis preview and the batch wizard, plus a
  `#column-detector` documentation section (#21).
- **WikiPathways viewer** on the single-analysis results page — lazy-loaded diagram with an
  enrichment-ordered picker and confidence badges. Context only; no expression overlay yet
  (#51).

### Fixed

- Shared GSEA links rendered as ORA because `/api/share` did not carry the analysis method.

## [3.0] — 2026-05-11 — Cleanup, Gene-Set Export & UX Gaps

### Added

- Gene-by-KE CSV export (`KE_ID, KE_Title, Gene_Symbol, log2FC, pvalue, FDR`) respecting the
  active significance thresholds.
- Real multipart file upload (CSV/TSV/TXT) on the main analysis page.
- `/demos` route with curated PXR and cisplatin cards and recommended-AOP calls to action.

### Changed

- AOP:472 now flows from SPARQL alongside the other kidney AOPs; the hardcoded fallback was
  removed.
- Column detector splits raw and adjusted p-value detection, surfaced in the preview picker.
- Cytoscape network exports carry per-KE p-value and FDR (#50).

### Removed

- Dead `build_cytoscape_network_legacy()`; duplicate `guess_id_type`; 14 naive
  `datetime.utcnow()` call sites replaced with timezone-aware constructors.

## [2.0] — 2026-03-03 — Batch Analysis & Dynamic AOPs

### Added

- Batch analysis wizard — up to 10 datasets per run, condition/timepoint/dose tagging, shared
  analysis settings, per-file progress via HTMX polling.
- Cross-condition comparison view with heatmap and matrix table.
- Dynamic AOP discovery — three-tier resolution (disk cache → AOP-Wiki SPARQL → bundled
  fallback) with typeahead and KE counts.

## [1.0] — 2026-02-27 — MVP

### Added

- Live integration with the molAOP Builder REST API for KE→gene mappings, with cross-process
  disk caching and graceful CSV fallback.
- VHP4Safety brand overhaul — shared base template, sticky navigation, Inter typography,
  card-based layouts.
- Documentation page covering input formats, statistical methods and result interpretation.
- Results page with metadata summary card, styled enrichment table and CSV export.
- Cytoscape.js network with dagre LR layout, per-KE progressive gene expansion, overlay zoom
  controls and transparent PNG export.
