# Changelog

All notable changes to the Molecular AOP Analyser are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions correspond to development milestones; the deployed service tracks `main`
via `ghcr.io/marvinm2/molaop-analyser`.

## [Unreleased]

### Added

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

- **Rows without a gene symbol no longer enter the background as a gene called `NAN`**
  (#80). The loader expanded multi-symbol rows with `str(row[id_col]).split('///')`, and
  `str()` on a missing value yields `'nan'` — uppercased into a literal identifier. Every
  symbol-less row in an upload therefore collapsed into one pseudo-gene. On a DESeq2 table
  of 13814 rows, 12939 of which carry a GeneSymbol and 12914 of which are distinct, the
  tool reported a background of 12915: the extra entry stood in for 875 discarded rows.
  Two things followed. The reported background was not the number of measured genes, so it
  could not honestly be quoted in a methods section; and if the symbol-less row passed the
  significance threshold, `NAN` counted as a significant gene that can never overlap a Key
  Event gene set, biasing every Fisher test very slightly toward the null.

  Rows whose identifier is null, empty or a missing-value placeholder are now dropped
  before expansion, so the pseudo-gene is gone and the reported background is the true
  count of measured genes. A `A///NA` row keeps its real symbol. The placeholder
  vocabulary is pandas' own `STR_NA_VALUES` plus the punctuation forms (`-`, `--`, `.`,
  `?`, …) that R, Excel and array-annotation pipelines write for "no symbol". The number
  of discarded rows is counted, but is **not yet shown anywhere in the interface** —
  reporting it next to the background size is still open under #80.

  A file in which *no* row yields a usable identifier — most often a wrong ID column —
  now fails with a validation error naming that column, instead of the generic
  "unexpected error" page it produced once the pseudo-gene stopped standing in for those
  rows.

- **The driver-gene export comes out in a deterministic order** (#82). `build_gene_tracking`
  appended its records in dict-iteration order over each condition's stored KE→gene map, so
  the row order of `genes_export.csv` was an artefact of how that JSON happened to be
  written; two runs over identical result data produced files that were equal only once
  sorted, while `compare_export.csv` was already byte-stable. That made the gene export
  unusable for diff-based reproducibility checks and generated spurious churn when exports
  were version-controlled. Records are now ordered at the source, and the export re-applies
  the same order: Key Event numerically (so `KE:10` follows `KE:9`), then gene symbol, then
  condition **in upload order** — sorting condition labels as text would put `10uM` before
  `2uM` and break a dose series. The summary view keeps its shared-first ordering (a gene
  driving three conditions stays above one driving a single condition).

  The export was the only affected artefact: the Genes tab payload was already sorted, and
  the batch report does not use these records.

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
