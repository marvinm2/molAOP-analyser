# Pitfalls Research

**Domain:** Flask enrichment analysis web app — API integration, UI/UX overhaul, Cytoscape.js network performance
**Researched:** 2026-02-25
**Confidence:** MEDIUM-HIGH (critical pitfalls verified across multiple sources; some UX findings LOW confidence where only WebSearch)

---

## Critical Pitfalls

### Pitfall 1: In-Memory Cache Does Not Survive Multi-Worker Deployment

**What goes wrong:**
The existing `CacheManager` in `cache_manager.py` is a plain Python dictionary protected by a threading lock. Once the API integration is live, the cached reference sets (KE→gene mappings) fetched from the Builder API are stored per-process. Under Gunicorn or uWSGI with multiple workers (the Docker Compose setup implies this production path), every worker fetches its own copy independently — so the Builder API receives `N_workers × N_cache_misses` requests per cache cycle, rapidly exhausting the 100-requests/hour rate limit.

**Why it happens:**
Developers test with Flask's single-process dev server (`python app.py`), where a shared in-memory dictionary works correctly. The leap to multi-worker production is transparent at the code level; nothing warns that `cache.set()` in worker 1 is invisible to worker 2. The existing code already calls `load_cached_reference_sets()` at module import time (`reference_sets = load_cached_reference_sets()`) — each worker will re-execute this import and independently trigger a cold fetch.

**How to avoid:**
Before activating live API calls, replace the in-memory `CacheManager` with a shared backend. Redis is the standard solution and is straightforward to add to the existing Docker Compose. Flask-Caching with `CACHE_TYPE = 'RedisCache'` is a drop-in replacement. Alternatively, since the reference data changes infrequently, write the fetched payload to a shared SQLite table or local JSON file with a timestamp check — the SQLite infrastructure already exists in `database.py`.

**Warning signs:**
- Builder API `429 Too Many Requests` errors appearing shortly after app startup or after container restarts
- Multiple identical log entries like `"Loaded N KE sets from reference files"` appearing within seconds of each other (one per worker)
- Cache hit rate is consistently zero in application logs despite repeated analysis runs

**Phase to address:**
API integration phase — implement before any live API call goes to production. Do not leave this until a "polish" phase.

---

### Pitfall 2: Rate Limit Exhaustion During Full Reference Set Fetch

**What goes wrong:**
The Builder API enforces 100 requests per hour per IP. The current `load_reference_sets()` fetches data from three static CSV files in a single pass. The equivalent API operation requires paginated calls to `/api/v1/mappings` (max 200 results per page), plus separate requests per KE or pathway to resolve gene memberships. For a dataset with hundreds of KE-WP mappings, building the full `KE → gene set` dictionary may require 20–50+ requests per cold fetch — consuming half or all of the hourly budget in one startup cycle.

**Why it happens:**
The API reference shows only per-endpoint rate limits; developers naturally model the API interaction as a single data fetch rather than accounting for pagination depth. The `/exports/gmt/ke-wp` bulk export endpoint exists precisely to address this but is easy to overlook when the natural integration target appears to be `/api/v1/mappings`.

**How to avoid:**
Use the bulk export endpoint (`/exports/gmt/ke-wp?min_confidence=High`) as the primary data source for the full KE→gene reference set. This returns the entire dataset in one HTTP request. Reserve the paginated `/api/v1/mappings` endpoint for targeted queries (e.g., fetching a specific AOP's KE list). Cache the bulk export result with a long TTL (6–24 hours) using a shared cache backend.

**Warning signs:**
- Analysis requests start failing with "No reference gene sets found" during busy periods
- Builder API logs show a spike of requests from the analyser IP shortly after startup
- `429` responses visible in Flask logs but appearing intermittent (some workers hit the limit, others have cached data)

**Phase to address:**
API integration phase — design the fetch strategy before writing any integration code.

---

### Pitfall 3: Analyser Startup Fails When Builder API Is Unreachable

**What goes wrong:**
`app.py` currently calls `reference_sets = load_cached_reference_sets()` at module import time (line 788). When the data source is live API calls instead of local CSV files, any network error during startup will prevent the Flask application from launching at all. This turns the Builder API into a hard dependency for the analyser's availability — if the Builder restarts, redeploys, or is temporarily unreachable, the analyser goes offline entirely.

**Why it happens:**
The startup-time data load pattern is appropriate for local CSV files (which are always present) but becomes a single point of failure when the source is a network resource. The pattern looks safe until the Builder has downtime.

**How to avoid:**
Implement a graceful degradation strategy:
1. On startup, attempt to load from shared cache (Redis or file-based) first. If a valid cached payload exists, use it regardless of age — log a warning that data may be stale.
2. If no cache exists, attempt the API fetch with a timeout (e.g., 5 seconds). On failure, fall back to the last-known-good CSV files shipped with the application.
3. Move the startup fetch to a background thread so the application becomes available immediately, completing the fetch asynchronously.

The existing `cache_manager.py` TTL logic already handles expiry correctly — extend it to support "stale-on-error" semantics.

**Warning signs:**
- Application container fails health checks immediately after a Builder redeploy
- Docker Compose shows the analyser exiting with a stack trace pointing to the reference set loader
- No clear error boundary in logs — the failure looks like an import error rather than a network timeout

**Phase to address:**
API integration phase — build fallback before removing the CSV files from the repository.

---

### Pitfall 4: KE–Gene Data Schema Mismatch Between CSV and API

**What goes wrong:**
The existing `helpers.py` `load_reference_sets()` merges three CSV files by joining on `WP_ID`, `WPID`, and `GeneID` fields with explicit normalisation (uppercase, strip, int→str casting). The Builder API returns `ke_id` (e.g., `"KE 55"`) and `pathway_id` (e.g., `"WP1234"`). These identifier formats differ from the existing CSV column conventions. If the API integration reuses the existing normalisation logic without verifying the actual field names and formats returned by the live API, the merge will silently produce an empty result set — the enrichment analysis will run but find no gene sets.

**Why it happens:**
The API reference document uses slightly different naming conventions to the CSV column headers (`ke_id` vs `KE_ID`, `pathway_id` vs `WP_ID`). When integrating quickly, developers often assume field names are consistent and skip validation of the actual API response shape.

**How to avoid:**
Write an explicit schema validation step as the first action after any API response. Assert the presence of expected keys before attempting any merge. Add an integration test that fetches one page from the live Builder API and validates the response structure. Do not delete the existing CSV files until the API-sourced data produces identical `reference_sets` dictionaries in a comparison test.

**Warning signs:**
- `reference_sets` dictionary has zero entries after switching to API data source
- Enrichment analysis raises `"No reference gene sets found"` immediately after API integration
- Merge produces a DataFrame with zero rows (the `dropna(subset=['GeneName'])` step silently removes everything)

**Phase to address:**
API integration phase — write the schema contract test before writing the transformation code.

---

## Moderate Pitfalls

### Pitfall 5: Dagre Layout Freeze on Gene-Expanded Networks

**What goes wrong:**
The current network renders KE nodes only (dozens of nodes at most), which Cytoscape.js handles without issue. The "Add Gene Nodes" button in the results template adds hundreds of gene nodes at runtime via `addGeneNodes()`. Running the Dagre layout algorithm over the full expanded graph (dozens of KEs + hundreds of gene nodes + thousands of edges) can freeze the browser tab for several seconds, sometimes indefinitely on lower-end hardware. Users interpret this as the application crashing.

**Why it happens:**
Dagre is a general-purpose DAG layout designed for directed hierarchies of moderate size. Its complexity scales non-linearly with node count. Gene nodes are typically positioned in large clusters around their parent KE — a layout algorithm is unnecessary for them; only relative positioning to the KE node matters.

**How to avoid:**
Do not run a full layout after gene node expansion. Instead:
1. Pre-position gene nodes programmatically around their parent KE node using the KE node's `cy.position()` as an anchor point, distributing genes in a radial or grid pattern without invoking any layout algorithm.
2. Use `cy.batch()` to add all gene nodes atomically, preventing incremental redraws.
3. Restrict Dagre layout to KE-only subgraphs; never run it on the full expanded network.

**Warning signs:**
- Browser tab becomes unresponsive for >3 seconds after clicking "Add Gene Nodes"
- Users reporting the network "crashes" when enabling gene visibility
- JavaScript `RangeError: Maximum call stack size exceeded` in browser console after expansion

**Phase to address:**
Network visualisation phase — implement programmatic gene positioning before enabling gene node expansion in the redesigned UI.

---

### Pitfall 6: In-Template JSON Blobs Break on Large Datasets

**What goes wrong:**
The current `analyze` route passes `network_json`, `ke_gene_json`, `volcano_json`, and `table_json` as pre-serialised JSON strings embedded directly into the Jinja2 template via `{{ ... | safe }}`. For large datasets (e.g., the cisplatin 72hr 50μM dataset with thousands of genes), these blobs can be several hundred kilobytes. The HTML page size grows proportionally, causing slow initial page loads and — in extreme cases — browser rendering delays. The report generation form compounds this by also embedding `enrichment_results` and `volcano_data` as form textarea content.

**Why it happens:**
Embedding JSON in HTML templates is the simplest pattern for passing server-computed data to client-side JavaScript. It works well at small scale. As datasets grow, no natural failure mode alerts the developer — pages just get slower.

**How to avoid:**
During the UI overhaul, refactor data delivery for the three largest payloads (`ke_gene_json`, `volcano_json`, `network_json`) to use dedicated API endpoints fetched by the client after page load. The Flask route returns a lightweight page shell; JavaScript fetches data from `/api/analysis/<id>/network` and `/api/analysis/<id>/genes` asynchronously. The experiment ID already exists in the database (`experiment_id`), making server-side data storage straightforward.

**Warning signs:**
- HTML source of results page exceeds 500KB
- Chrome DevTools shows the `results.html` document taking >2s to parse
- Form submission for report generation times out due to POST body size

**Phase to address:**
UI/UX overhaul phase — refactor template data delivery as part of the page restructure. Address before adding more data to the results view.

---

### Pitfall 7: CSS Brand Overhaul Breaks Cytoscape.js Canvas Sizing

**What goes wrong:**
Cytoscape.js determines its canvas size from the `#cy` container's dimensions at initialisation time. When a brand overhaul changes the page layout — particularly adding sticky headers, side panels, flex/grid wrappers, or altering `box-sizing` globally — the `#cy` div may have zero height or incorrect dimensions when Cytoscape initialises, rendering a blank or misshapen network canvas. This is a well-known CSS interaction with canvas-based libraries.

**Why it happens:**
The current CSS uses a fixed `height: 800px` on `#cy`. A new layout system that uses viewport-relative sizing (`vh`, `calc()`, flex `flex-grow`) may cause the container to have no intrinsic height until the browser completes layout, which may happen after Cytoscape's initialisation call. Adding a new sticky header also shifts the available viewport height.

**How to avoid:**
During the CSS migration, explicitly set an `min-height` on the `#cy` container that does not depend on parent flex sizing. After adding the new layout structure, test Cytoscape initialisation in isolation by calling `cy.invalidateSize()` and `cy.fit()` after a short timeout or after the page's `DOMContentLoaded` event is confirmed complete. Add a smoke test that measures `cy.width()` and `cy.height()` and asserts they are non-zero.

**Warning signs:**
- Cytoscape canvas is invisible (zero height) or shows network nodes clustered in the top-left corner
- `cy.width()` returns 0 in the browser console immediately after page load
- Network appears correctly only after manually resizing the browser window (which triggers a resize event)

**Phase to address:**
UI/UX overhaul phase — test Cytoscape canvas sizing explicitly as part of every layout change.

---

### Pitfall 8: Brand Colours Applied Inconsistently Between Cytoscape Styles and CSS

**What goes wrong:**
The VHP4Safety palette is applied in two entirely separate styling systems: CSS/Jinja2 templates and Cytoscape.js stylesheet objects defined in JavaScript. When the brand palette is updated (e.g., adding the full secondary palette, adjusting hex values), developers update one system but miss the other. The resulting network visualisation uses different shades of blue, magenta, or teal than the surrounding page, breaking the unified brand appearance that the overhaul is meant to establish.

**Why it happens:**
There is no single source of truth for brand colour values. CSS files define colours as hex literals. Cytoscape stylesheets define colours as hex literals in JavaScript. Both need updating independently, and there is no build-time mechanism to enforce consistency.

**How to avoid:**
Define brand colours once as CSS custom properties (`--primary-dark: #29235C`) in a `:root` block. For Cytoscape styles, read these values at initialisation time using `getComputedStyle(document.documentElement).getPropertyValue('--primary-dark')`. This makes CSS custom properties the single source of truth for both systems, eliminating the double-maintenance problem.

**Warning signs:**
- Network node colours visibly differ from the brand palette used in table headers, buttons, and the page header
- Pull requests changing brand colours in `style.css` do not touch Cytoscape stylesheet definitions
- Design review comments noting the network "looks different" from the rest of the tool

**Phase to address:**
UI/UX overhaul phase — implement CSS custom properties as the first step, before applying them to any UI element.

---

### Pitfall 9: API Confidence Level Filtering Changes Enrichment Results Silently

**What goes wrong:**
The Builder API supports `?confidence_level=High|Medium|Low` filtering. If the integration fetches only High-confidence KE-WP mappings (a reasonable default), the gene sets backing each KE will be smaller than the current CSV-derived sets (which include all mappings regardless of confidence). Enrichment p-values will change, and some KEs that previously appeared significant may no longer be — not because the biology changed, but because the reference set shrank. Researchers comparing old results to new results will see discrepancies with no explanation.

**Why it happens:**
Developers choose `confidence_level=High` as a sensible data quality filter without realising it changes the statistical background. The application shows no indication of which confidence level was used, so results appear to have changed arbitrarily.

**How to avoid:**
Make the confidence level filter an explicit, configurable application parameter stored in `Config`. Display the active confidence filter on the results page alongside background gene set size. If a default of `High` is chosen, document it in the UI and in the analysis report. Run a comparison test between CSV-derived and API-derived reference sets at each confidence level to quantify the difference before switching.

**Warning signs:**
- Enrichment results for known positive controls (e.g., the PXR agonist demo datasets against the liver steatosis AOP) change after API integration with no code changes to the analysis logic
- `total_KE_genes_in_dataset` values in the results table are systematically lower than before
- The number of KEs with sufficient gene overlap (≥5 genes) decreases noticeably

**Phase to address:**
API integration phase — decide and document the confidence level strategy before replacing CSV files.

---

## Technical Debt Patterns

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Keep in-memory `CacheManager` for API data | No Redis dependency in dev | Cache is per-process; rate limits exhausted in multi-worker production | Never — fix before live API |
| Inline Cytoscape stylesheet hex colours in JS | Simple, readable | Brand update requires touching multiple files; colours drift from CSS | Only in a single-file prototype |
| Embed full gene-expression JSON in HTML template | No additional endpoint needed | Page size grows with dataset; browser parse time increases | Acceptable for <200 genes; not for full cisplatin datasets |
| Run Dagre on full gene-expanded network | Single layout call | Browser freeze on any dataset with >100 genes; user-visible crash | Never — position genes programmatically |
| Delete CSV fallback files immediately after API integration | Cleaner repo | No recovery path when Builder API is down | Only after confirmed stable API with monitored uptime |

---

## Integration Gotchas

| Integration | Common Mistake | Correct Approach |
|-------------|----------------|------------------|
| Builder API `/api/v1/mappings` pagination | Fetch only first page (50 results), assume it is complete | Follow `pagination.next` until `null`; or use `/exports/gmt/ke-wp` for full dataset in one request |
| Builder API KE identifier format | Assume `KE_ID` from CSV matches `ke_id` from API | Validate: API returns `"KE 55"` (with space); existing code uses `KE_ID` column which may be numeric. Normalise explicitly. |
| Builder API rate limit | Treat 100 req/hr as generous headroom | 100 req/hr = ~1.67 req/min. A single multi-page fetch at startup + 3 workers = limit exhausted in minutes |
| Builder API `aop_id` filter | Assume it resolves KEs client-side | Builder resolves via AOP-Wiki SPARQL; network latency adds 200–500ms per AOP query |
| CSS `box-sizing` global reset | Apply `* { box-sizing: border-box }` in overhaul | Verify this does not change how Cytoscape computes its container dimensions; test with `cy.invalidateSize()` |
| VHP4Safety font (Inter) via CDN | Load from Google Fonts at render time | Self-host or use a preconnect hint; Google Fonts CDN adds 200–400ms on first load for researchers behind institutional firewalls |

---

## Performance Traps

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| Full Dagre layout on gene-expanded network | Browser freeze 3–30 seconds; tab unresponsive | Programmatic gene positioning around KE anchor points; never call layout on mixed KE+gene graph | Always — Dagre is O(n²) for dense graphs; manifests immediately with any cisplatin dataset |
| Cytoscape bezier edges with gradient colors | Network render drops to <5 FPS when scrolling | Use straight-line (`haystack`) curve style for KE-gene edges; reserve bezier for KE-KE edges only | >50 gene nodes visible simultaneously |
| Cold API fetch blocking Flask startup | Application takes 10–30s to start; health checks fail | Fetch reference data asynchronously after startup; serve cached data immediately | Every deployment cycle |
| JSON blobs >200KB in HTML template | `results.html` parse time >2s; report form submission fails | Move large payloads to dedicated API endpoints fetched post-load | Cisplatin datasets (full 72hr 50μM has ~15,000 genes) |
| Cytoscape `cy.add()` called in a loop without `cy.batch()` | Each `add()` triggers a redraw; adding 500 genes takes 5+ seconds | Wrap all element additions in `cy.batch(function(){ ... })` | >50 elements added in a single operation |

---

## Security Mistakes

| Mistake | Risk | Prevention |
|---------|------|------------|
| Trusting `ke_id` and `pathway_id` values from the Builder API without sanitisation before inserting into SQLite | Low (same-org API, no user content), but defence-in-depth principle applies | Validate that `ke_id` matches `KE \d+` pattern and `pathway_id` matches `WP\d+` before storing |
| Exposing Builder API base URL in client-side JavaScript | Enables direct requests from browser, bypassing the analyser's rate-limit management and caching | All Builder API calls must go through the Flask backend; never expose the Builder URL in templates |
| No timeout on outbound Builder API requests | A slow Builder response blocks the Flask worker thread indefinitely | Set `requests.get(url, timeout=10)` on all outbound calls; document timeout as a config constant |

---

## UX Pitfalls

| Pitfall | User Impact | Better Approach |
|---------|-------------|-----------------|
| No loading state shown while KE-gene network data is fetching from API | Researcher clicks "Analyse" and sees a blank results page for 5–15 seconds; assumes the tool crashed | Show a spinner with status text ("Fetching KE reference data…") immediately on form submit; update status as each pipeline stage completes |
| Brand overhaul removes or repositions the "Add Gene Nodes" button without announcement | Power users who discovered this feature can no longer find it; it is undocumented | Keep the gene expansion workflow prominent in the redesigned network panel; add a tooltip or inline help text |
| Significant KE threshold (FDR < 0.05) changed without version note | Regulatory users comparing runs before and after overhaul see different significant KEs; raises reproducibility concerns | Never change statistical defaults without a version changelog visible in the results page and analysis report |
| Network legend uses different colour definitions than the new brand palette | Researcher looks at the legend and the network node colours do not match; reduces trust | Generate legend swatches from the same CSS custom properties used by Cytoscape styles |
| Inter font renders differently at small sizes than the existing system font | Enrichment table values (p-values, gene lists) become harder to read at 13–14px in Inter | Use Inter for headings and labels; keep a monospace or system-ui fallback for dense data tables |

---

## "Looks Done But Isn't" Checklist

- [ ] **API integration complete:** Verify `reference_sets` dictionary produced from API contains the same KE IDs as the CSV-derived version on at least two AOP case studies (liver steatosis, kidney failure) — not just that the endpoint returns 200.
- [ ] **Cache working in production:** Confirm cache is shared across workers by deploying with 2+ Gunicorn workers and verifying only one API fetch occurs on startup — check Builder API request logs, not just analyser logs.
- [ ] **Rate limit headroom:** Simulate 10 consecutive analysis requests (different AOP selections) and verify no `429` responses are received from the Builder API.
- [ ] **Network canvas sized correctly after CSS overhaul:** Load results page in Chrome, Firefox, and Safari; verify `cy.width()` and `cy.height()` are non-zero without manual resize.
- [ ] **Gene node expansion does not freeze browser:** Add gene nodes on the kidney AOP with a cisplatin dataset; verify completion in <2 seconds on a mid-range laptop.
- [ ] **Brand colours consistent:** Export a screenshot of the results page; verify Cytoscape node colours match table header and button colours pixel-for-pixel using a colour picker.
- [ ] **Fallback to CSV when Builder is offline:** Stop the Builder API container; verify the analyser starts, loads reference data from cache or CSV fallback, and processes an analysis successfully.

---

## Recovery Strategies

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| In-memory cache exhausting rate limits in production | MEDIUM | Add Redis to Docker Compose; change `CacheManager` backend; redeploy. No data migration needed. |
| Dagre freeze shipped to production | MEDIUM | Hot-fix: disable "Add Gene Nodes" button client-side; deploy programmatic positioning fix; re-enable. |
| CSS overhaul breaks Cytoscape canvas | LOW | Add `cy.invalidateSize(); cy.fit();` on `window.addEventListener('resize', ...)` as immediate patch; fix root CSS cause in follow-up. |
| API schema mismatch produces empty reference sets | HIGH | Roll back to CSV data source immediately (keep CSV files in repo until confirmed equivalent); investigate field name mapping; re-test; re-deploy. |
| Confidence level filter silently changes enrichment results | MEDIUM | Re-fetch with `confidence_level=` matching the CSV scope (likely all levels); add config parameter; communicate change to users via results page notice. |

---

## Pitfall-to-Phase Mapping

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| In-memory cache breaks under multi-worker deployment | API integration | Deploy with 2+ workers; confirm single API fetch on startup via Builder logs |
| Rate limit exhaustion during full reference set fetch | API integration | Log all outbound API calls; assert <10 requests per analysis cycle |
| Startup fails when Builder API is unreachable | API integration | Stop Builder container; confirm analyser starts successfully within 10s |
| KE–gene data schema mismatch | API integration | Run comparison test: CSV-derived vs API-derived `reference_sets`; assert identical KE IDs |
| Dagre freeze on gene-expanded network | Network visualisation | Automated browser test: add gene nodes; assert completion in <2000ms |
| In-template JSON blobs too large | UI/UX overhaul | Measure `results.html` size for largest cisplatin dataset; assert <150KB HTML |
| CSS overhaul breaks Cytoscape canvas sizing | UI/UX overhaul | Assert `cy.width() > 0` after page load in Playwright smoke test |
| Brand colours inconsistent between CSS and Cytoscape | UI/UX overhaul | Visual regression test: screenshot comparison of network vs page chrome |
| API confidence level filter changes results silently | API integration | Run PXR agonist demo; compare significant KE list to CSV baseline; document any differences |

---

## Sources

- Cytoscape.js official performance documentation: [cytoscape/cytoscape.js performance.md](https://github.com/cytoscape/cytoscape.js/blob/master/documentation/md/performance.md)
- Cytoscape.js WebGL renderer preview (January 2025): [blog.js.cytoscape.org](https://blog.js.cytoscape.org/2025/01/13/webgl-preview/)
- Cytoscape.js discussion — performance at scale: [GitHub Discussion #3088](https://github.com/cytoscape/cytoscape.js/discussions/3088)
- Flask multi-worker cache isolation: [Medium — Sharing data across Gunicorn workers](https://medium.com/@jgleeee/sharing-data-across-workers-in-a-gunicorn-flask-application-2ad698591875)
- Flask concurrency and in-memory caching: [Medium — The Flask Concurrency Trap](https://medium.com/@engage.deepakagarwal/the-flask-concurrency-trap-why-the-quick-threading-fix-breaks-at-scale-8ea773c58a88)
- API rate limiting strategies: [blog.poespas.me — Flask + Redis caching](https://blog.poespas.me/posts/2025/03/02/handling-large-scale-api-caching-with-flask-and-redis/)
- Cytoscape.js dagre layout scaling: [cytoscape/cytoscape.js-dagre](https://github.com/cytoscape/cytoscape.js-dagre)
- Cytoscape.js layout strategy guide: [blog.js.cytoscape.org — Using layouts](https://blog.js.cytoscape.org/2020/05/11/layouts/)
- Data migration pitfalls from static to API sources: [Alation — Data Migration Challenges](https://www.alation.com/blog/5-pain-points-moving-data-to-cloud/)
- ETL pipeline design pitfalls: [Airbyte — ETL Pipeline Design Pitfalls](https://airbyte.com/data-engineering-resources/etl-pipeline-pitfalls-to-avoid)
- Builder API reference: `/home/marvin/Documents/Services/molAOP-analyser/molAOP-analyser/KE-MAPPING-API-REFERENCE.md` (internal)
- Existing codebase analysis: `app.py`, `cache_manager.py`, `helpers.py`, `services/network_service.py` (internal)

---
*Pitfalls research for: molAOP Analyser — API integration, UI/UX overhaul, Cytoscape.js network performance*
*Researched: 2026-02-25*
