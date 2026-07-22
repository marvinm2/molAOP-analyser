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
  worth reading down a dose series, and it keeps the small per-gene effects with a consistent
  direction that a threshold discards by construction. Getting that table before meant
  running every condition as a separate single analysis and assembling NES and FDR by hand,
  which also threw away the harmonised background that makes batch results comparable at all.

  The batch form gains the same method selector as the single-analysis form, and the method
  is stored on the batch rather than per condition so every condition is scored identically.
  The log2FC and p-value thresholds stay on screen under GSEA: they do not enter the GSEA
  statistic, but the run still uses them to decide which genes count as significant, and
  those counts are what the comparison page's per-condition strip and the driver-gene tab
  are built from — so the form now says that instead of hiding controls whose numbers you
  are then shown.

  Everything downstream follows the batch's method rather than assuming Fisher:

  - **Comparison view** — heatmap and table switch to NES + FDR. The heatmap colours on the
    signed NES over a diverging scale centred at zero, so a coordinated up-shift and a
    down-shift are no longer both rendered as "significant", and the matrix keeps a
    sub-threshold NES rather than blanking it, because a consistent weak signal across a
    dose series is itself the observation. Clicking a condition header sorts the largest
    NES to the top (for FDR it still sorts the smallest to the top).
  - **Comparison network, delta mode** — the delta slices now show the *change in NES*
    against the reference condition on a red/blue divergent scale with the shift size in
    the shading, and the legend says so. Previously a downward shift was drawn the same way
    a loss of significance is drawn (muted grey), which read as weaker evidence when it
    meant the opposite direction.
  - **Matrix export** — NES joins FDR and -log10(FDR) in the export dropdown and in
    `GET /batch/<uuid>/compare/export?matrix=nes`, so what is on screen can be downloaded.
  - **Batch report (PDF and HTML)** — the per-condition tables print NES, the permutation
    FDR and the leading-edge genes, the comparison heatmap and matrix carry NES, and the
    header states the method. Before this, a GSEA batch's report was rendered through the
    Fisher columns and printed `0` overlap and `0.00` odds ratio for every Key Event while
    silently dropping NES and the leading-edge genes.

  `batches.method` is added by the established idempotent PRAGMA-then-ALTER migration.
  SQLite backfills existing rows with `ora`, and a NULL from any other path is read as `ora`
  too, so every batch created before this change still opens, runs, exports and compares
  exactly as it did.

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

- **Batch condition pages state the gene-set resources the batch actually used** (#74).
  `/batch/<uuid>/condition/<n>` re-renders the results page from stored condition data, but
  never carried the batch's resource selection into that context — and the template
  defaulted the missing value to `WikiPathways`. A batch run over all three resources
  therefore reported a single-resource provenance, confidently and wrongly, while the batch
  summary for the same run reported all three. Which resources contributed gene sets decides
  which Key Events were testable at all, so the header was misstating the analysis it sat on
  top of, in the one place a reader would copy from into a methods section. The condition
  page now carries both the requested selection (#55) and the per-resource resolution and
  warnings (#68), and the template no longer substitutes a real resource name for an absent
  one — an unrecorded selection reads as *"Not recorded"*.

  The same page was misstating two more of the run's settings, both of which travel into a
  methods section together with the resources. The mapping-confidence threshold read *"All
  mappings"* on every condition page, even one run at *"High only"* — and even directly above
  the warning banner that the stricter threshold had raised. And a report exported from a
  condition page reported no column mapping and a p-value cutoff of `0.05` whatever the batch
  had actually been run at. All of these now come from the batch record, as does the
  enrichment method: the page had hardcoded over-representation, so with batch GSEA (#76) it
  would have shown a GSEA run the ORA colour scale and exported it as ORA.

  **A report exported from a condition page reported the wrong resources too**, which is the
  form the claim is most likely to be quoted in. The report form posted no resource list at
  all, so the generator fell back to the browser session — meaning a PDF exported from a
  batch condition page named the resources of whatever single analysis had last been run in
  the same session, and named `WikiPathways` when there was no session to borrow from. The
  form now posts the run's own selection, always, and the report no longer invents a resource
  name: a run with nothing recorded reports *"Not recorded"*. The same invented default is
  gone from the stored-experiment JSON. It remains in the **batch report** (the multi-
  condition PDF), which is not covered by this change.

  The dataset identity leaked the same way: the report generator preferred the browser
  session over the values the page posted, so a PDF exported from a condition page carried
  the dataset name, stressor, dosing and owner of whatever single analysis had last been run
  in that session, under a page header naming the condition. The posted values now win.

- **Hub genes on a batch condition page show their p-value** (#75). The hub-gene panel reads
  its "Adj. p-value" from the per-gene data stored with each condition, and the batch runner
  never captured p-values at all — so every hub row of every condition of every batch showed
  an em dash, while the same file analysed on its own showed a number. Batch runs now capture
  the raw and adjusted p-value columns from each uploaded file the way a single analysis
  does, independently of the one column the batch thresholds on, and store them alongside the
  fold changes. A file that genuinely carries no adjusted p-value column still reads as a
  dash, which is the honest answer.

  This applies to batches run from now on. Conditions stored before this change hold no
  p-values and will keep showing a dash; re-run the batch to get them.

- **"No gene set mapped" no longer covers Key Events that are mapped** (#81). The
  per-condition accounting had two exclusion counters, and one of them was answering three
  different questions at once. A Key Event with no curated mapping, a Key Event whose
  mapping exists but could not be resolved to genes (#79), and a Key Event whose gene set is
  known but went unmeasured in the uploaded dataset were all reported as *"excluded (no gene
  set mapped)"*. They call for opposite responses: the first is a curation gap to fill in the
  molAOP Builder, the second is stale reference data in this tool and should be reported as a
  bug, the third is a fact about the experiment's coverage and about nobody's curation. The
  case that exposed it was KE 1115 on AOP:472, listed among the unmapped while carrying a
  live, high-confidence mapping to `WP5477` — a data-pipeline failure presented as missing
  curation. (#79 has since made `WP5477` resolvable, so that particular Key Event is now
  tested; the reporting bug it exposed outlived it.)

  Exclusions are now counted under three separate reasons — no curated mapping, mapped but
  unresolvable, and fewer than `MIN_KE_GENES` measured genes (zero measured genes included,
  since that is a coverage number, not a curation one) — and the mapped-but-unresolvable
  clause **names the pathway IDs** so the claim can be checked against the Builder. The
  wording still comes from the single authority `format_ke_summary()`, so the results page,
  the single report and the batch report say the same thing; the network styles the new
  reason apart from an uncurated Key Event, and both the results page and the shared-link
  view carry that style and legend entry. Summaries and networks stored before the split
  carry only the old reasons and render exactly as they did.

  Telling the two apart needs one piece of information that used to be thrown away: which
  Key Event each unresolvable pathway belongs to. The KE→pathway mappings are merged with
  pathway membership on an inner join, so a Key Event whose only pathway is unresolvable
  never reaches the reference sets at all and is indistinguishable, downstream, from one
  nobody has mapped. `load_reference_sets()` now builds that map *before* the merge and
  carries it on the returned gene sets, which is what lets both the single-analysis route
  and the batch runner hand it to the enrichment backend, and what puts the pathway IDs on
  the network nodes so a batch report or shared link rebuilt from stored network JSON still
  names them. A Key Event in that position now reads *"1 excluded (mapped, but no genes
  could be resolved: WP4010)"* rather than being counted among the uncurated.

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
