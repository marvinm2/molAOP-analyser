# Changelog

All notable changes to the Molecular AOP Analyser are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions correspond to development milestones; the deployed service tracks `main`
via `ghcr.io/marvinm2/molaop-analyser`.

## [Unreleased]

### Added

- Repository tooling brought in line with the molAOP Builder: CI workflow (test matrix,
  flake8, coverage, startup smoke test), code-quality workflow (ruff/black/isort, bandit,
  pip-audit), Dependabot, CODEOWNERS, issue and pull-request templates, `LICENSE` (GPL-2.0),
  `pyproject.toml` and this changelog.
- Design-token layer in `static/css/main.css` mirroring the Builder's `:root` block —
  VHP4Safety palette, spacing, radius, shadow and z-index scales as CSS custom properties.

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
