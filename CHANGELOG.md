# Changelog

All notable changes to the Molecular AOP Analyser are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions correspond to development milestones; the deployed service tracks `main`
via `ghcr.io/marvinm2/molaop-analyser`.

## [Unreleased]

## [6.0.0] - 2026-08-21

### Changed

- **The enrichment background is now a choice, and it is reported.** Two settings define it,
  kept separate because they are separate decisions and bundling them would hide that the
  second can undo the first (#132).

  **Gene universe** — `measured` (default) keeps every gene with a usable identifier and a
  fold change; `testable` keeps only genes with a value in the selected p-value column, which
  is what the loader always did. The old rule dropped every gene DESeq2's independent
  filtering withheld an adjusted p-value from, and that selects against exactly the profile
  an induced response has — low at baseline, strongly up on treatment. In a lead-exposure
  dataset `HMOX1` showed a 22-fold induction at *p* = 4.8e-5, was given no adjusted p-value,
  and so was absent from the background entirely: it could not be counted as a hit for any
  gene set, at any dose.

  **Across conditions** — `union` (default), `intersection` (the old behaviour), or
  `per_condition`. Intersecting reads as conservative and is not: each condition's universe
  has already been shaped by that condition's own filtering, so the intersection collapses the
  batch onto whichever condition was filtered hardest. On a four-dose batch it came out at
  4,086 genes, equal to the 50 µM contrast exactly — that set being a subset of all three
  others — against a union of 18,339.

  The consequence ran both ways. Inducible gene sets were hollowed out; constitutively
  expressed sets occupied a larger share of a smaller background, which inflated their
  expected overlap and deflated their enrichment.

  Both settings are stored on the batch and rendered wherever the background appears — the
  progress partial, the summary page and the batch report — so a gene count is never shown
  without the rules that produced it. **The tested gene list is unchanged**: a gene with no
  adjusted p-value belongs in the background, not among the hits.

  **Batches created before this keep their own rules** (`testable` + `intersection`) through
  NULL-coercing accessors, so existing shared links and exported reports reproduce exactly.

- **Harmonisation applies to ORA only.** It was previously applied before the backend
  dispatch, so a GSEA batch had its ranked list truncated to the genes shared with every
  other condition. GSEA has no 2×2 table and no background; that removed ranking information
  without matching anything. A code comment claiming it "is what makes the comparison
  legitimate" went with it (#132).

### Fixed

- **The gene-expression legend described a scale the code does not draw.** Two independent
  errors, both long-standing, both on the page results are read from. The legend was labelled
  −2 / +2 while `defaultLogFCColor()` clamps at **±1**, so every gene colour was read as half
  its actual intensity and a gene at +5 was indistinguishable from one at +1 with nothing
  saying the scale had run out. Separately, the legend's gradient bar was a five-stop RdYlBu
  ramp with a **yellow** midpoint, while the function interpolates blue → **white** → red; the
  endpoints agreed, which is why it went unnoticed for so long. Both are corrected on the
  results and shared-results pages, and `tests/test_colour_scale.py` now pins the labels and
  the CSS stops to the clamp. **No node colour changes** — only the description of them, so
  figures produced earlier remain valid.

  Surfaced by #51: the pathway overlay generates its ramp *from* the function, so it was
  correct on arrival and put two contradictory legends on one page.

- **The user guide said GSEA q-values come from `gseapy`.** Untrue since 5.0.0 (#122), which
  recomputes the pooled-null FDR here and ships an `fdr_source` column precisely so the two
  can be told apart. The statistical-methods section now describes what the tool does.

- **The "Software DOI" badge was not the concept DOI it claimed to be.** It pointed at
  `zenodo.org/badge/latestdoi/…`, which redirects to the newest *version* DOI
  (`10.5281/zenodo.21914318`), while the citation section told readers it was the concept DOI.
  Anyone following our own instructions cited a version-pinned identifier believing otherwise.
  The badge now resolves to the concept DOI `10.5281/zenodo.21914317`, and the section names
  both.

- **The hub-genes table could push the whole page sideways.** It had no horizontal-scroll
  wrapper while `results.css` applies `white-space: nowrap` to every cell, so a gene belonging
  to several Key Events made the "KE list" cell unbreakable. Affected desktop, not just narrow
  screens.

- **"Show gene nodes" could report a state the network did not have.** *Collapse All* removed
  every gene node without unticking the box, and *Expand All* added them without ticking it, so
  the checkbox and the graph disagreed and the next click on it appeared to do nothing. Both
  now keep it in sync, and toggling it scrolls the network into view — the control sits two
  cards below the thing it acts on. Unticking it also clears the expanded-Key-Event set, which
  it previously left stale, causing a later tap on a Key Event to open its drawer instead of
  re-adding its genes.

- **The GMT parsers now skip the Builder's provenance header explicitly.** The Builder prefixes
  every export with `#` lines whose own closing line states they are not gene sets. Both parsers
  survived them only because comment lines contain fewer than three tab-separated fields — so a
  single tab in a future header line would have been read as a gene set and silently polluted
  the reference data.

- Corrected several smaller claims: the documented network controls listed two buttons that do
  not exist and omitted the FDR cutoff selector, which is the single definition of significance;
  the batch tutorial pointed at a cisplatin demo panel that is hidden and whose data is not
  shipped; the Top 10% / 20% tooltips described the 10th and 20th percentiles when the code uses
  the 90th and 80th; the volcano legend hardcoded `p < 0.05` for a user-settable cutoff; the
  overlay's privacy note now says that a *binned* log2FC is what travels, matching the wording
  already on the results page; the About page said the tool performs Fisher's exact test and
  maps Key Events to WikiPathways alone, when GSEA has been first-class since 4.0 and GO BP and
  Reactome are both selectable; and the footer said 2025.

### Added

- **A "How to cite" section on the About page**, naming the three DOIs that are not
  interchangeable — this software, the Builder *dataset* the mappings come from, and the Builder
  application — plus the licence, the repository and a link to `/health`. None of this was
  reachable from the interface, for a tool whose README treats citation as a correctness matter.

- **A "Where the Gene Sets Come From" section in the user guide**, covering the three controls
  that change what a run reports and were previously undocumented: resource selection and what
  taking their union means for comparability, the minimum mapping confidence threshold and its
  one exception, and the provenance line that records how each resource was actually resolved —
  with the note that a methods section should quote the provenance, not the request.

### Changed

- The citation metadata now points at the Builder's real software DOI
  (`10.5281/zenodo.21914244`) and records this repository's own concept DOI
  (`10.5281/zenodo.21914317`) in `CITATION.cff`. Both landed after 5.0.0 was cut and were not
  previously recorded here.
- `docs/SOFTWARE-RELEASES.md`: the release-prep step now names **four** files, not three —
  `.zenodo.json`'s `version` was missing, and `tests/test_app_version.py` asserts it, so
  following the runbook literally failed CI. The DOI-backfill step now says which identifier
  belongs in which file, having previously been done in one place of three.
- `docs/KE-MAPPING-API-REFERENCE.md` re-verified against the live Builder (v2.9.0) and its
  header dated accordingly. Known upstream drift is listed rather than silently corrected: the
  `/dataset/*` and `/export/*` families return 503/500 live, `connection_type` is not
  WikiPathways-only, and the `go-mappings` CSV inserts `go_direction` at column 8, which breaks
  positional parsing tuned to the WikiPathways CSV.


- **Gene expression is overlaid on the embedded WikiPathways diagram (#51).** The pathway view
  showed a bare diagram; the genes of the selected Key Event that are significant in the uploaded
  dataset are now coloured by log2FC on the same blue-to-red ramp as the Cytoscape network, from
  the same `defaultLogFCColor()`, so the two views cannot disagree. An *Overlay expression data*
  checkbox (default on) returns the bare diagram, and the picker `<option>` carries `data-ke`
  because the colours come from the Key Event's gene set rather than from the pathway.

  This needed no rendering service. The WikiPathways Toolforge viewer already colours nodes from
  URL query parameters, where the parameter *name* is the colour and its value is a
  comma-separated node-selector list, and WikiPathways SVGs tag gene nodes with `HGNC_<SYMBOL>` —
  which matches this app's HGNC-only uppercase pipeline with no ID translation. The gene payload
  was already in the browser as `keToGenes`, so nothing extra is fetched and no new dependency,
  route or server-side rendering was introduced. The overlay had previously been deferred on the
  understanding that it required a PinPath R/Plumber microservice or `rpy2` in-container.

  **Only significant genes are coloured**, which makes an uncoloured node ambiguous — absent from
  the dataset, or measured and not significant. A counts line states both numbers and names the
  per-gene significance rule in force, via the new
  `results_context.describe_gene_significance_criterion()`. That rule differs by method: GSEA
  hides the threshold inputs, so `logfc_threshold` falls back to 0 and the flag becomes a bare
  p-value cut, qualifying many more genes than Fisher's exact test does on the same data.

  Colour bins must be **even** in number (8). With an odd count the middle bin's centre is exactly
  0.0, where `defaultLogFCColor()` returns its `#ccc` fallback rather than a ramp colour, so that
  bin coloured nothing and silently dropped significant genes with a near-zero fold change. The
  URL length cap measures the serialised query string rather than estimating a per-gene cost, and
  the reported dropped-gene count covers only genes the cap excluded.

  **This sends data to a third party.** With the overlay on, the gene symbols of the selected Key
  Event's significant genes travel to `pathway-viewer.toolforge.org` inside the diagram URL;
  log2FC values do not, only the colour bin. That reverses the pathway view's original design
  decision to keep expression data out of that URL, and the card discloses it next to the control
  that causes it.

### Fixed

- **The resource checkboxes no longer claim GO BP and Reactome have no confidence field.** Four
  labels reading `— no confidence field` sat directly above the paragraph stating that the
  threshold applies to all three resources. They were left behind when the threshold was extended
  to the GMT-backed resources (#71) and the prose making that claim was removed; the capability
  has been correct since 2026-07-22 and only the labels were wrong.

## [5.0.0] - 2026-08-13

### Changed

- **GSEA q-values are now computed here rather than taken from gseapy (#122).** This is a
  larger claim than the tool made before, so it is stated plainly: the FDR column is no
  longer a pass-through. A `fdr_source` column travels with every row, the results table
  and both reports render it, and where the recomputation could not run the row keeps
  gseapy's number and says so instead of being presented alongside recomputed ones.

  The reason is that a GSEA q-value is not a property of its own row. It is a ratio of tail
  fractions taken over pools shared by every Key Event in the run, so one Key Event's number
  can move another's. When a Key Event's observed enrichment score beats every permutation,
  gseapy cannot normalise it and substitutes NES = 1.0; #121 repaired that row's own values
  but could not touch anybody else's, because `gsea_fdr` had already run inside the Rust
  extension. The substituted value lands in the **observed** NES vector — the denominator
  every other Key Event's q is divided by.

  The bias is one-sided and derivable rather than random. With h = #{genuine NES ≥ NES\*} and
  p = #{genuine NES ≥ 0}, injecting one fabricated 1.0 leaves h alone and raises p when
  NES\* > 1.0 (q **inflated**, conservative), and raises both when 0 ≤ NES\* ≤ 1.0 (q
  **deflated**, anti-conservative). Negatively-enriched Key Events are untouched: a
  fabricated +1.0 never enters the negative pool.

  Measured on a six-Key-Event fixture containing one degenerate Key Event, a genuine
  non-degenerate neighbour at NES 1.03 moved from q = 0.364 to q = 0.273 — a 25% shift
  caused entirely by a different Key Event's substituted score, and in the direction the
  derivation predicts. On a run with no degenerate Key Event the recomputation reproduces
  gseapy's q-values exactly; that is asserted to ~1e-12 in
  `tests/test_gsea_service.py::TestPooledNullFdr` against gseapy's own `gsea_fdr`, and it
  is the check that makes the change reviewable.

  `recompute_pooled_fdr` is a transcription of gseapy's `normalize` and `gsea_fdr`, not an
  independent derivation — the point is to remove one contaminating value from the pools,
  so everything else must stay comparable to the implementation it replaces. Degenerate Key
  Events are counted at **±infinity** rather than excluded, which is what #122 proposed:
  excluding drops a Key Event that genuinely belongs in the observed distribution and so
  deflates the denominator for everyone else — the same class of error, smaller and
  opposite-signed. The observed score beat every permutation, so the normalised score is
  unbounded on that side and counting it at every same-signed threshold is the faithful
  limit, inventing no magnitude. Under that treatment the degenerate Key Event's own q comes
  out as 0.0, which is the limit #117 already writes into the row by hand — the two agree
  rather than having to be reconciled. `DEGENERATE_EXCLUDE` implements the other choice so
  the two can be compared on real data.

  The ±infinity is used inside the pooling only. The reported NES for such a Key Event stays
  NaN: the score is not normalisable, and printing an infinity in a results table would be a
  magnitude claim the data does not support.

### Fixed

- **Choosing an AOP and then touching a threshold silently reverted you to the demo's AOP.**
  On a demo dataset: press "Show all AOPs", pick something other than the recommended one,
  then click any log2FC quick button or "Update Plot" — and the picker snapped back to
  `AOP:DEMO` with the filter toggle back on "recommended". The threshold buttons auto-submit,
  so this happened without the user knowingly pressing anything, and the analysis then ran
  against an AOP they had not chosen.

  The threshold controls submit `#volcanoForm`; the AOP picker lives in `#enrichmentForm`.
  That split is why the resource checkboxes (#55) and the confidence select (#60) each needed
  a mirror — `aop_selection` and `aop_filter_mode` never got one, so they were not posted,
  `/preview` had nothing to echo, and the picker re-rendered with no value. The
  recommended-AOP auto-select then read the empty picker as "never chosen" and overwrote the
  pick, doing precisely what its own comment promised it would not:

  ```js
  // Only kicks in when the picker is empty (so we never clobber a manual selection …)
  if (recommendedAops.length === 1 && !hiddenInput.value && !searchInput.value) {
  ```

  The guard was correct in intent and could never hold, because after a swap the picker was
  always empty. Fixed by mirroring the selection, its label and the filter mode into
  `#volcanoForm`, and having `/preview` round-trip all three. Only the browser knows the
  display label — the server has the id — so the label rides along with it.

  The mirror is written both as hidden inputs and, authoritatively, into
  `evt.detail.parameters` on `htmx:configRequest`. Hidden inputs alone were not enough: they
  are only refreshed at the call sites that mutate the picker, so any other path that sets it
  (`tour.js`'s `preselectPxrAop`, for one) left them stale. This was caught by testing that
  path rather than reasoning about it.

- **A custom log2FC threshold was ignored unless you pressed "Update Plot" first.** The
  threshold inputs are in `#volcanoForm` and the Analyse button is in `#enrichmentForm`;
  `pval_threshold` had a hidden `#pval-mirror` kept in sync on every keystroke, `logfc_threshold`
  had no mirror at all. Typing `2.5` into Custom and pressing "Run Enrichment Analysis"
  therefore ran on the *previous* threshold while the box on screen read `2.5` — wrong numbers,
  no error, nothing on the results page to reveal it. Given the same mirror the p-value box has.

- **The demo datasets were unusable: "Confirm Column Selection" did nothing.** Picking a
  dataset on `/demos` rendered the preview correctly, and then the button below it was dead —
  no error, no swap, nothing. The demo path was the *only* way into the tool that broke this
  way; uploading your own file worked throughout, which is why it survived so long.

  `/preview` sweeps `uploads/` of anything older than `Config.UPLOAD_RETENTION_HOURS` before
  it resolves the requested file, and the demo copy was made with `shutil.copy2` — which
  preserves the source mtime. A demo file shipped in `data/` therefore arrived in `uploads/`
  carrying the mtime of whenever that file was last written, already far outside the retention
  window. The first request rendered fine (the copy existed for the length of that request);
  the second request — the Confirm button — swept it away *before* looking for it and answered
  `400 File not found`. HTMX does not swap on a 4xx, so the page simply sat there.

  The copy is now made with `shutil.copy` plus an explicit `os.utime`, so the retention clock
  starts at the copy rather than at the shipped file's age. The same trap hit the guided tour
  (which preloads a demo dataset) and any second `/preview` render such as "Update Plot".
  Batch demo copies use `copy2` too but are unaffected — batch directories expire from
  database timestamps, not mtime.

  Worth noting how the timing hid this in production: `COPY` carries the build context's
  mtimes into the image, so a freshly built image has fresh demo files and the demo flow works
  — until the image is a day old, at which point every demo click starts failing. The
  regression test ages the source file deliberately rather than trusting it, because a fresh
  git checkout stamps `data/` with the checkout time and would let the old behaviour pass on
  CI.

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

- **Nightly backup of the deployed database** (#114). `/data/molAOP_analyser.db` holds
  shared-result links and batch history — state that has survived redeploys by design since
  2026-07-21 — and nothing copied it anywhere. GlusterFS replica 2 is not a backup: it
  protects against a node dying, not against deletion, corruption or a bad schema migration,
  because both replicas take the same write. The service had taken three additive migrations
  with no recoverable prior state.

  `scripts/backup_db.py` uses SQLite's Online Backup API, which is safe while the app is
  writing, and runs as the Swarm cron job `cronjobs_molaop-analyser-backup` (daily 03:30 UTC,
  seven-day retention) off this repo's own image — the same pattern as the Builder's, so the
  backup logic stays versioned with the schema it belongs to. It reads the database file
  directly on GlusterFS rather than exec'ing into the running container, so it is
  node-independent.

  Three of its checks exist because the obvious ones pass on the failure they are meant to
  catch. SQLite treats a zero-length file as a valid empty database, so `integrity_check`
  returns `ok` on an aborted backup — a byte floor runs first. A structurally valid database
  with no rows is not a backup of *this* database — `experiments` and `batches` must be
  non-empty. And `Connection.backup()` retries a locked source *forever* rather than honouring
  the connection timeout, so the read lock is taken separately under a bounded wait: a job
  that hangs would, with `skip-running=true`, silently suppress every run after it. Any check
  that fails deletes the partial file, so a failed run never leaves something that sorts as
  the newest backup.

- Repository tooling brought in line with the molAOP Builder: CI workflow (test matrix,
  flake8, coverage, startup smoke test), code-quality workflow (ruff/black/isort, bandit,
  pip-audit), Dependabot, CODEOWNERS, issue and pull-request templates, `LICENSE` (GPL-2.0),
  `pyproject.toml` and this changelog.
- Design-token layer in `static/css/main.css` mirroring the Builder's `:root` block —
  VHP4Safety palette, spacing, radius, shadow and z-index scales as CSS custom properties.

### Changed

- **Python 3.11 → 3.14** (#88). Dependabot proposed this as a one-line Docker bump, which is
  the one way it cannot safely be done: `requirements.txt` and `requirements-dev.txt` are
  pip-compile locks resolved *inside* the image's base, and CI pins its interpreter to match,
  so moving the base alone would have left the locks and CI describing an interpreter the
  image no longer runs. Done here as the single change the docs prescribe — both locks
  regenerated in `python:3.14-bookworm`, the suite run there, and the Dockerfile plus all
  three `python-version` keys updated together.

  **The pins did not move.** Regenerating without `--upgrade` preserves existing pins where
  they still resolve, and every one did: the locks are byte-identical to their 3.11 versions
  apart from the `with Python 3.14` header. That is the desired outcome for an interpreter
  bump — one variable changes, so anything that breaks is attributable to Python and not to
  a simultaneous dependency upgrade. Notably `kaleido==0.2.1` and pycairo both install
  cleanly on 3.14.

  Verified before merge, in the 3.14 base: both locks compile, the shipped lock installs,
  **593 tests pass**, the image builds, and a container from it serves `/`, `/demos`,
  `/documentation` and `/about` at 200. The SQLite "unable to open database file" line in a
  mountless local run is pre-existing and identical on the 3.11 image — it disappears with
  the production `/data` mount.

### Removed

- **The "under development" banner on the landing page.** The orange band reading *"This tool
  is under development — results and interface are subject to change"* has been removed, along
  with its now-unused `.alert-band` rule in `static/css/pages/index.css`.

### Added

- **The mapping-confidence threshold now applies to GO BP and Reactome, not just
  WikiPathways** (#71). The threshold has been offered on both forms since #60, but it only
  ever filtered the Builder's curated WikiPathways mappings — a "High only" run over all
  three resources filtered one component and used the other two whole. #67 made the Analyser
  *disclose* that rather than fix it, because the fix was not available: the Builder's GMT
  `?min_confidence=` selected a **single tier**, so forwarding the user's choice would have
  dropped every high-confidence mapping when they asked for "Medium and High" — worse than
  not filtering at all.

  molAOP-builder#206 shipped on 2026-07-22 and the parameter now means *at-or-above*, with
  unrecorded confidences always included, so `services/api_service.fetch_gmt_reference_sets()`
  forwards it. Verified against the live Builder: `GO_BP` goes 12 KEs / 4,334 genes at `all`
  to 6 KEs / 1,791 at `high`, `Reactome` 8 / 755 to 4 / 144, and each level's Key Events are a
  subset of the one below — the nesting that distinguishes a threshold from the old tier
  selector.

  Two traps handled explicitly. Our vocabulary is `("all", "medium", "high")` and the
  Builder's whitelist is `{high, medium, low}`: `all` is **not** in it and returns HTTP 400,
  which would drop the resource for the whole run, so `gmt_min_confidence_param()` translates
  it to an omitted parameter and degrades anything unrecognised the same way. And the
  reference-set cache key already carried the threshold, which now matters materially rather
  than defensively — the three levels fetch genuinely different gene sets, so a shared key
  would serve a `high` run the sets cached for `all`.

  A run that falls back to the reference files bundled in the image is still unfiltered —
  those carry no confidence column — and both forms now say the threshold applies to all
  three resources with that one exception, instead of naming WikiPathways alone.

### Fixed

- **A maximally enriched Key Event is no longer reported as not enriched** (#117). GSEA
  normalises an observed enrichment score by the mean of the *same-signed* tail of that gene
  set's permutation null. When the null falls entirely on the other side of zero that tail is
  empty, the normalisation is undefined, and gseapy substitutes `NES = 1.0` with `p = 1.0` —
  the values a thoroughly unenriched gene set produces. The one case where every permutation
  failed to reach the observed score was therefore printed as the null result, and nothing in
  the table gave the reader a way to tell. On AOP 472 it hid "Increase, DNA damage" in
  carboplatin at 48 h, which an independent implementation puts among the strongest calls in
  the grid.

  The condition is now read off the null distribution rather than off the reported values: the
  size of the same-signed tail, not `NES == 1.0`. That matters because the failure is
  continuous. A tail of one or two permutations still yields a plausible-looking NES, but the
  magnitude is a ratio against a sample of two — and no check on the output values can find
  it. Results now carry a `nes_status` of `ok`, `unstable_normalisation`,
  `beyond_permutation_resolution` or `undiagnosed`, the same-signed count that produced the
  call, the achievable `p_value_resolution`, and the raw `ES` so direction survives even where
  the NES does not.

  **The p-value goes with the magnitude.** An earlier draft of this change said the nominal
  p-value survived a short same-signed tail and only the NES did not; that claim is retracted.
  gseapy divides the exceedance count by the *same-signed* tail, not by the permutation count,
  so a term with a tail of five can only report p in steps of 0.2 and a term with a tail of one
  only 0 or 1. Measured on gseapy 1.1.4, a term with a tail of 56 and one exceedance returns
  1/56, not 1/100 — and the AOP 472 grid contains a p of 0.600 that is three out of five and a
  p of 0.000 that is zero out of one. Under `unstable_normalisation` neither number may be
  ranked or thresholded on, the step size is reported per row as `p_value_resolution`, and the
  reports print the p-value with that qualifier attached.

  Where the tail is empty the NES is reported as undefined rather than as a number — it is
  genuinely not normalisable — while the p-value becomes the resolution bound `1/permutations`
  and the FDR its limit of zero, because no permutation was as extreme as the observation.
  The cell reads as maximally significant, which is what it is, instead of as null.

  **The distinction now reaches every view, not only the single-analysis table.** The batch
  comparison path — the cross-condition grid, its exports and both batch reports — received a
  recovered Key Event as a bare `NES = None` and rendered it as an em-dash: unshaded, sorted
  last, dropped from the heatmap and stripped even of the FDR every other cell carries, i.e.
  indistinguishable from a Key Event that was never tested and the exact inverse of the
  result. `nes_status` and `ES` now travel through `build_comparison_matrix` into the table,
  the heatmap (those cells are ringed rather than coloured — there is no magnitude to put on a
  diverging scale, only a direction) and both report tables. On the network the same cell
  reached `nesColor(null)`, whose `parseFloat(nes) || 0` folded it onto `#FFFFFF`, the colour
  reserved for a NES of exactly zero: a maximally enriched Key Event drawn white with a
  significance border, inherited by the `.cyjs` export and the network PNG in both reports. The
  node now carries `nes_status` and `es`, and the state has its own off-gradient treatment and
  legend entry. Two report captions that promised "every NES value, significant or not" were
  false for exactly these cells and have been corrected.

  **The diagnostics no longer fabricate statistics under concurrency.** The capture rebinds a
  process-global (`gseapy.gsea.prerank_rs`), and this app is threaded: Flask runs with
  `threaded=True` and every batch condition is dispatched on its own thread, so two overlapping
  runs are routine. Unserialised, the patches nested, both sinks received both runs' summaries,
  and — since the per-term counts are keyed by KE ID, which is identical across conditions —
  one condition's tail count landed on another condition's row, where the repair then rewrote
  that row's NES, p and FDR on the strength of it. The same nesting left the global permanently
  bound to a dead wrapper whenever the inner block outlived the outer one, which is the normal
  case. The patch and the `gp.prerank` call it exists for now run under one module-level lock,
  taken before the global is read, and the restore is identity-checked.

  **A failure to capture is no longer reported as a clean result.** If gseapy renames the entry
  point, or `method='multilevel'` routes through `fgsea_rs` instead, no null can be inspected.
  Those rows were labelled `ok` — so a future gseapy bump would have quietly reintroduced the
  original misreporting under a clean bill of health. They are now labelled `undiagnosed`,
  logged at error level naming the symbol, and counted into the run's own Key Event accounting
  sentence, which says in as many words that their NES and p-value are unchecked.

  Incidence rises with gene-set size, since a larger set tightens the null: this gets *more*
  likely as Key Event curation improves. Across the AOP 472 grid (5 gene-set pools × 9
  conditions, 315 cells) the six affected cells all carry 379–817 measured genes.

  Known residual, not addressed here: gseapy computes the FDR for the *unaffected* terms from
  a pooled null that still contains the degenerate terms' substituted `NES = 1.0`. Correcting
  that means recomputing the FDR outside gseapy, which is a larger change than this one and is
  tracked as #122.

- **A Key Event above the GSEA size ceiling is no longer reported as under-measured** (#120).
  `gseapy.prerank` discards a gene set larger than `max_size` before computing anything, and
  the exclusion accounting attributed every size-based drop to `min_size`. A Key Event with
  1200 measured genes was reported — on the results page, on the network and in both reports —
  as "excluded (fewer than 5 measured genes)", which points the reader at the opposite repair
  from the one needed. Over-ceiling drops now have their own exclusion reason and accounting
  clause, the log line names the real bound and the real gene count, and the ceiling itself is
  a named constant (`Config.GSEA_MAX_KE_GENES`) rather than a literal, so a run can report
  what it could not test. The value stays at 1000 rather than the Broad convention of 500:
  three of the seven testable AOP 472 Key Events already exceed 500 measured genes.

- **An unfiltered resource is no longer blamed on the bundled reference files** (#71). The
  disclosure added with #67 split its explanation on whether a resource appeared in the
  filterable list, which was the same thing as "has a bundled CSV fallback" only for as long
  as WikiPathways was the sole filterable resource. Once GO BP and Reactome became filterable,
  a run recorded before #71 — whose stored resolution marks them unfiltered — would have
  rendered *"served from the bundled reference files"* for resources that **have no bundled
  fallback and never did**. The reason is now attributed from the resource's own fallback path
  (`CSV_FALLBACK_RESOURCES`), so the bundled-files explanation is only ever given where it is
  true, and everything else says plainly that the threshold was not applied. This is the same
  class of defect as #108: a disclaimer that reads as diligence while stating something false.

- **The unresolved-pathway warning describes the AOP being analysed** (#108). The warning
  added with #79/#81 — *"N mapped pathway(s) could not be resolved to genes … their coverage
  is understated here"* — was computed over the whole reference universe rather than the
  selected AOP, so it fired on every run whether or not that run had lost anything. On
  AOP:472 it named `WP1234`, `WP3980` and `WP4010`, none of which is mapped to any of its Key
  Events; every one of its fourteen pathways resolved and no Key Event lost a single gene.
  This is the inverse of the failure #79 and #81 were about, and worse: rather than hiding a
  real gap it asserts one that does not exist, on the same page the numbers are read from. An
  author writing up that result would reasonably add a caveat that their coverage was
  understated, and the caveat would be false. A warning that fires unconditionally is also
  one readers learn to ignore — including on the runs where it is true.

  `scope_resolution_to_aop()` now narrows the resolution's `unresolved_pathways` and
  `unresolved_ke_pathways` to the run's own Key Events, and the block is suppressed entirely
  when nothing survives. Scoping happens where the run is created — the single-analysis route
  and `_persist_and_launch_batch()`, both of which already know the AOP — so the stored
  resolution is correct and the results page, the report, the batch summary, the condition
  pages and shared links all inherit it without repeating the lookup. An AOP whose topology
  cannot be loaded leaves the resolution unscoped rather than emptied: an over-broad warning
  is recoverable, a silently suppressed one is not. Runs stored before this change keep the
  unscoped accounting they were saved with.

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
