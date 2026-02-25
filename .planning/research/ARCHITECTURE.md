# Architecture Research

**Domain:** Flask bioinformatics web application — enrichment analysis with live API integration, UI modernisation, and Cytoscape.js network visualisation
**Researched:** 2026-02-25
**Confidence:** HIGH (existing codebase inspected directly; patterns verified against official docs and community sources)

---

## Current Architecture (Baseline)

### System Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                        Browser (Client)                          │
│  index.html (upload/preview)   results.html (network + table)   │
│  Cytoscape.js canvas  │  Volcano plot (Plotly-style inline JS)  │
└─────────────────────┬───────────────────────────────────────────┘
                      │  HTTP (form POST / JSON)
┌─────────────────────▼───────────────────────────────────────────┐
│                    Flask (app.py)  — Route Layer                 │
│  /           /preview    /analyze    /generate_report            │
│  /documentation          /api/upload_network_png                 │
│  CSRF protection (flask-wtf)  │  Session management              │
└──────┬──────────────┬────────┴──────────────────────────────────┘
       │              │
┌──────▼──────┐  ┌────▼────────────────────────────────────────┐
│ validation  │  │                Services Layer                │
│ .py         │  │  data_service.py   enrichment_service.py     │
│ exceptions  │  │  network_service.py  column_detector.py       │
│ .py         │  │  report_service.py   gene_id_validator.py     │
└─────────────┘  └────────────────┬────────────────────────────┘
                                  │
┌─────────────────────────────────▼────────────────────────────┐
│                      Reference Data Layer                     │
│  helpers.py  load_reference_sets()  ←→  cache_manager.py     │
│  ┌─────────────────┐                                         │
│  │   data/ CSVs    │  (currently: static files on disk)      │
│  │  aop_ke_map.csv │                                         │
│  │  KE-WP.csv      │                                         │
│  │  edges_wpid_to_gene.csv                                   │
│  │  node_attributes.csv                                      │
│  │  aop_ker_edges.csv                                        │
│  │  ke_metadata.csv│                                         │
│  └─────────────────┘                                         │
└──────────────────────────────────────────────────────────────┘
┌──────────────────────────────────────────────────────────────┐
│                    Persistence Layer                          │
│  database.py (SQLAlchemy)  →  molAOP_analyser.db (SQLite)    │
│  ExperimentRecord model — metadata + enrichment summaries    │
└──────────────────────────────────────────────────────────────┘
```

### Current Component Boundaries

| Component | File(s) | Responsibility | Communicates With |
|-----------|---------|----------------|-------------------|
| Route layer | `app.py` | HTTP handling, session, CSRF, file I/O | All services, validation, cache |
| Data ingestion | `services/data_service.py` | Load/validate gene expression files, duplicate resolution via Fisher's method | Called by routes |
| Enrichment | `services/enrichment_service.py` | Fisher's exact test, FDR correction, contingency tables | Reference sets from cache |
| Network builder | `services/network_service.py` | Assembles Cytoscape.js node/edge JSON | Enrichment results, KE metadata |
| Column detector | `services/column_detector.py` | Heuristic/statistical column type inference | Called during /preview |
| Report generator | `services/report_service.py` | PDF and HTML report assembly | All analysis results |
| Gene ID validator | `services/gene_id_validator.py` | Gene symbol normalisation and validation | data_service |
| Reference loader | `helpers.py` | Merges CSV data into KE→gene dicts | data/ CSV files |
| Cache | `cache_manager.py` | In-memory TTL cache (thread-safe, Lock-guarded) | helpers.py, reference sets |
| Persistence | `database.py` | SQLAlchemy ORM, SQLite experiment history | Called from routes after analysis |
| Config | `config.py` | Upload limits, AOP case studies, demo datasets, ExperimentMetadata dataclass | All modules |
| Validation | `validation.py`, `exceptions.py` | Input guards, typed exception hierarchy | Routes |

---

## Target Architecture (After Three Milestones)

The three development areas — API integration, UI overhaul, and network performance — each affect a distinct layer. They are largely independent and can be built in parallel once the API client boundary is established.

### System Overview (Target)

```
┌────────────────────────────────────────────────────────────────────┐
│                        Browser (Client)                             │
│  base.html (Inter font, VHP4Safety palette, sticky nav)            │
│  index.html (modernised upload form)                               │
│  results.html (modernised results, network panel)                  │
│  Cytoscape.js — WebGL renderer (preview), haystack edges           │
│  Progressive gene layer: collapsed by default, expand on click     │
└──────────────────────┬─────────────────────────────────────────────┘
                       │
┌──────────────────────▼─────────────────────────────────────────────┐
│                    Flask Route Layer  (app.py — unchanged routes)   │
│  All existing routes preserved; templates swapped underneath        │
└──────┬───────────────────────────────────────────────────────────┬─┘
       │                                                           │
┌──────▼──────────┐         ┌──────────────────────────────────────▼─┐
│ validation.py   │         │           Services Layer                │
│ exceptions.py   │         │  (unchanged data_service, enrichment,   │
└─────────────────┘         │   network, column_detector, report)     │
                            └───────────────┬──────────────────────┬──┘
                                            │                      │
                            ┌───────────────▼──────────────┐       │
                            │    NEW: api_client.py         │       │
                            │  BuilderAPIClient class       │       │
                            │  - httpx or requests session  │       │
                            │  - retry with backoff         │       │
                            │  - timeout config             │       │
                            │  - raises BuilderAPIError     │       │
                            └───────────────┬───────────────┘       │
                                            │                       │
                            ┌───────────────▼──────────────┐  ┌────▼────────────────┐
                            │  MODIFIED: helpers.py         │  │  data/ CSVs         │
                            │  load_reference_sets()        │  │  (fallback only)    │
                            │  → tries API first            │  │  aop_ke_map.csv     │
                            │  → falls back to CSVs         │  │  KE-WP.csv etc.     │
                            │  → cache_manager unchanged    │  └─────────────────────┘
                            └───────────────────────────────┘
                            ┌──────────────────────────────────────────────────────┐
                            │  Persistence Layer (database.py — unchanged)          │
                            │  SQLite experiment history                            │
                            └──────────────────────────────────────────────────────┘
```

---

## Component Boundaries (Target State)

| Component | Responsibility | Communicates With | Change? |
|-----------|---------------|-------------------|---------|
| `app.py` routes | HTTP handling, session, CSRF, file I/O | Services, cache, validation | No |
| `services/data_service.py` | Gene expression loading, duplicate resolution | Routes | No |
| `services/enrichment_service.py` | Fisher's exact, FDR correction | Reference sets from cache | No |
| `services/network_service.py` | Cytoscape.js JSON assembly | Enrichment results | No |
| `services/column_detector.py` | Column type inference | Preview route | No |
| `services/report_service.py` | PDF/HTML generation | All results | No |
| `services/gene_id_validator.py` | Gene symbol normalisation | data_service | No |
| **NEW: `services/api_client.py`** | HTTP client for molAOP Builder API | External Builder API | New |
| **MODIFIED: `helpers.py`** | API-first reference set loading with CSV fallback | api_client.py, data/ | Modified |
| `cache_manager.py` | In-memory TTL cache | helpers.py | No |
| `database.py` | SQLite persistence | Routes | No |
| `config.py` | All configuration including Builder API URL | All modules | Extended |
| `validation.py`, `exceptions.py` | Input validation, typed exception hierarchy | Routes | Extended |
| **MODIFIED: `templates/`** | UI — base layout, index, results | Browser | Replaced |
| `static/css/style.css` | Visual styling | Browser | Replaced |

---

## Data Flow

### Current Reference Data Loading

```
app startup / first /analyze request
    ↓
load_cached_reference_sets()  [app.py]
    ↓
cache_manager.get_or_set("reference_sets", ..., ttl=3600)
    ↓
helpers.load_reference_sets(ke_wp_path, wp_gene_path, node_path)
    ↓
pd.read_csv() × 3 CSV files  [data/ on disk]
    ↓
pandas merge operations → {KE_ID: set(gene_symbols)} dict
    ↓
returned to analyze route → passed to enrichment_service
```

### Target Reference Data Loading (API-first with fallback)

```
app startup / first /analyze request
    ↓
load_cached_reference_sets()  [app.py — unchanged call]
    ↓
cache_manager.get_or_set("reference_sets", ..., ttl=3600)
    ↓
helpers.load_reference_sets()  [modified signature or env-driven]
    ↓
BuilderAPIClient.get_ke_gene_mappings()
    ├── HTTP GET /api/v1/mappings (paginated, with retries + backoff)
    ├── assembles {KE_ID: set(gene_symbols)} — same output shape as before
    └── on BuilderAPIError → falls back to pd.read_csv() CSV path (unchanged logic)
    ↓
{KE_ID: set(gene_symbols)} dict  [identical contract to current]
    ↓
cache_manager stores result → enrichment_service receives it unchanged
```

**Key invariant:** The dict contract `{KE_ID: set(gene_symbols)}` must not change. enrichment_service, network_service, and all callers remain untouched.

### Analysis Request Flow (unchanged)

```
POST /analyze
    ↓
validate_form_data()  →  load_and_validate_data()
    ↓
process_gene_expression()  →  load_aop_data()
    ↓
load_cached_reference_sets()  (may be cached)
    ↓
run_enrichment_analysis()  →  build_cytoscape_network()
    ↓
build_ke_gene_mapping()  →  db_manager.save_experiment_metadata()
    ↓
render_template("results.html", ...)  →  browser
```

### UI Modernisation Data Flow

```
base.html (new)
    ├── Inter font (CDN or local)
    ├── brand CSS variables (#29235C, #E6007E, #307BBF, etc.)
    ├── sticky navigation component
    └── block slots: head_extra, content, scripts_extra

index.html extends base.html
    └── upload form (same POST targets, same field names)

results.html extends base.html
    ├── enrichment table (same table_json, ke_gene_json data bindings)
    └── network panel (same network_json, modernised Cytoscape init)
```

Routes pass the same template variables; only the HTML shell changes.

---

## Architectural Patterns

### Pattern 1: API-First with Graceful CSV Fallback

**What:** `helpers.py` attempts live Builder API call first. On any `BuilderAPIError` or network failure, it silently falls back to reading the static CSVs already on disk. The caller (cache layer) sees no difference.

**When to use:** Any time an external dependency is optional but preferred. The CSVs represent a "last-known-good" snapshot of the same data.

**Trade-offs:** Slightly more complex helper; adds one HTTP dependency to startup. Fallback ensures analyser keeps working if Builder is down.

**Example:**
```python
# services/api_client.py
class BuilderAPIClient:
    def __init__(self, base_url: str, timeout: int = 10, max_retries: int = 3):
        self.base_url = base_url
        self.timeout = timeout
        self.max_retries = max_retries
        self.session = requests.Session()
        adapter = HTTPAdapter(max_retries=Retry(
            total=max_retries,
            backoff_factor=0.5,
            status_forcelist=[500, 502, 503, 504]
        ))
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

    def get_ke_gene_mappings(self) -> dict:
        """Fetches all KE→gene mappings from Builder API.
        Returns {KE_ID: set(gene_symbols)} — same shape as CSV path.
        Raises BuilderAPIError on any failure.
        """
        ...

# helpers.py (modified)
def load_reference_sets(ke_wp_path, wp_gene_path, node_path,
                        api_client=None):
    if api_client:
        try:
            return api_client.get_ke_gene_mappings()
        except BuilderAPIError as e:
            logger.warning(f"Builder API unavailable, falling back to CSVs: {e}")
    # existing CSV logic unchanged below
    ...
```

### Pattern 2: Jinja2 Template Inheritance for UI Swap

**What:** Introduce `base.html` with all brand variables, Inter font, sticky nav, and CSS block slots. All existing templates (`index.html`, `results.html`) become child templates that `extend base.html`. Route code and template variable bindings are untouched.

**When to use:** UI overhaul where Python route logic must not change but the entire visual shell needs replacing.

**Trade-offs:** One-time restructuring effort; after that, all pages share the brand base automatically.

**Build order implication:** Create `base.html` first. Then update `index.html`, then `results.html`. Each can be tested independently against the same Flask routes.

**Example:**
```jinja2
{# templates/base.html #}
<!DOCTYPE html>
<html lang="en">
<head>
  <link rel="stylesheet" href="{{ url_for('static', filename='css/style.css') }}">
  {% block head_extra %}{% endblock %}
</head>
<body>
  <nav class="vhp-nav"><!-- sticky brand nav --></nav>
  {% block content %}{% endblock %}
  {% block scripts_extra %}{% endblock %}
</body>
</html>

{# templates/results.html #}
{% extends "base.html" %}
{% block content %}
  <!-- enrichment table and network panel, same data bindings -->
{% endblock %}
```

### Pattern 3: Cytoscape.js Progressive Rendering

**What:** Separate the network into two render passes. Pass 1 renders only KE nodes and KER edges (the AOP skeleton) immediately on page load. Pass 2 adds gene nodes only when a KE is clicked, using the existing `ke_gene_json` data already in the page.

**When to use:** Networks where gene nodes dominate total node count, causing initial canvas congestion and slow layout.

**Trade-offs:** Adds JS event-handling logic; requires the gene layer to be toggleable. Already has the data available in `ke_gene_json` — no server round-trips needed.

**Example:**
```javascript
// Initialise with KE nodes only
const cy = cytoscape({
  container: document.getElementById('cy'),
  elements: networkData.nodes.filter(n => n.data.ke_type !== undefined)
             .concat(networkData.edges),
  style: cytoscapeStyles,
  layout: { name: 'dagre' },
  hideEdgesOnViewport: true,   // performance option
  textureOnViewport: true      // performance option
});

// Add genes on KE click
cy.on('tap', 'node[ke_type]', function(event) {
  const keId = event.target.id();
  const genes = keGeneData[keId] || [];
  cy.batch(() => {
    genes.forEach(g => {
      if (!cy.getElementById(g.id).length) {
        cy.add({ group: 'nodes', data: { id: g.id, label: g.id, log2FC: g.log2FC } });
        cy.add({ group: 'edges', data: { source: keId, target: g.id } });
      }
    });
  });
  cy.layout({ name: 'cola', animate: false }).run();
});
```

---

## Recommended Project Structure (Target)

```
molAOP-analyser/
├── app.py                     # Routes — unchanged
├── config.py                  # Extended: BUILDER_API_URL, BUILDER_API_TIMEOUT
├── helpers.py                 # Modified: API-first loading, CSV fallback
├── cache_manager.py           # Unchanged
├── database.py                # Unchanged
├── validation.py              # Minor: add BuilderAPIError to exception map
├── exceptions.py              # Add: BuilderAPIError class
├── utils.py                   # Unchanged
├── services/
│   ├── api_client.py          # NEW: BuilderAPIClient with retry adapter
│   ├── data_service.py        # Unchanged
│   ├── enrichment_service.py  # Unchanged
│   ├── network_service.py     # Unchanged
│   ├── column_detector.py     # Unchanged
│   ├── report_service.py      # Unchanged
│   └── gene_id_validator.py   # Unchanged
├── templates/
│   ├── base.html              # NEW: brand base layout (Inter, VHP4Safety palette)
│   ├── index.html             # REPLACED: extends base.html
│   ├── results.html           # REPLACED: extends base.html
│   └── documentation.html     # REPLACED: extends base.html
├── static/
│   └── css/
│       └── style.css          # REPLACED: CSS custom properties for brand palette
└── data/                      # Unchanged (CSV fallback)
```

### Structure Rationale

- **`services/api_client.py` as a new file:** Keeps HTTP client concerns out of `helpers.py`, which already has a clear responsibility (data assembly). The client is independently testable and mockable.
- **`helpers.py` modified, not replaced:** The CSV merge logic is stable and serves as the fallback. Adding an optional `api_client` parameter is a backwards-compatible extension.
- **Templates replaced, routes untouched:** Jinja2 inheritance means the Python side sees no change. All `render_template()` calls pass the same variable names — the new templates just present them differently.
- **`config.py` extended:** `BUILDER_API_URL` and `BUILDER_API_TIMEOUT` belong in config, not hardcoded in the client. Load from environment variable with a sensible default.

---

## Build Order (Dependencies Between Components)

```
Phase A: API Client Layer (no UI dependency)
  1. exceptions.py — add BuilderAPIError
  2. services/api_client.py — BuilderAPIClient with retry + timeout
  3. helpers.py — integrate api_client, add fallback logic
  4. config.py — add BUILDER_API_URL, BUILDER_API_TIMEOUT
  5. Unit tests for api_client (mock HTTP responses)
  6. Integration test: helpers.load_reference_sets() with Builder running vs. offline

Phase B: UI Modernisation (no API dependency, parallel with A)
  1. static/css/style.css — CSS custom properties, brand palette, Inter font
  2. templates/base.html — sticky nav, brand colours, block slots
  3. templates/index.html — extends base.html, same form fields
  4. templates/results.html — extends base.html, same data bindings
  5. templates/documentation.html — extends base.html
  6. Visual regression check: same routes, same data, new appearance

Phase C: Network Performance (depends on B for template context)
  1. Modify Cytoscape init in results.html: add hideEdgesOnViewport, textureOnViewport
  2. Implement progressive gene layer: KE skeleton first, genes on tap
  3. Set edge curve-style to haystack for KE-gene edges
  4. Test with large AOP (kidney network, organ network models)
  5. Evaluate WebGL renderer (provisional in Cytoscape.js 3.31+; test only)
```

**Key dependency:** Phase C requires B (needs the new results.html shell). Phases A and B are independent and can be developed simultaneously. Phases A and B have no dependencies on each other.

---

## Anti-Patterns

### Anti-Pattern 1: Replacing helpers.py signature entirely

**What people do:** Delete CSV logic from `helpers.py`, make it API-only, break the fallback path.

**Why it's wrong:** The Builder may be unreachable (different Docker container, network issue, Builder restart). Without CSVs, the analyser becomes completely non-functional during Builder downtime. The CSV files are a stable last-known-good dataset.

**Do this instead:** Keep the CSV merge code intact. Inject `api_client` as an optional parameter. Try API, fall back to CSV on any error.

### Anti-Pattern 2: Changing template variable names during UI overhaul

**What people do:** Rename `table_json` to `enrichment_json`, `network_json` to `cy_data` etc. to match a new naming convention in the new templates.

**Why it's wrong:** Any variable rename in templates requires a matching change in every `render_template()` call in `app.py`. This is unnecessary risk during a UI-only change and pollutes the diff.

**Do this instead:** Keep all `render_template()` keyword arguments identical. The new templates reference the same names — only the HTML structure changes.

### Anti-Pattern 3: Loading gene nodes eagerly in Cytoscape.js

**What people do:** Include all gene nodes in the initial `networkData` payload and let Cytoscape try to render hundreds of nodes at once on page load.

**Why it's wrong:** The `ke_gene_json` already contains all genes for every KE. If a network has 15 KEs each with 30 genes, that is 450 gene nodes on canvas before any interaction. Cytoscape.js canvas renderer becomes sluggish at this scale. The user wants to see the AOP structure first, not a gene cloud.

**Do this instead:** Use the progressive rendering pattern (Pattern 3 above). Genes appear only when the relevant KE is clicked. The `ke_gene_json` data is already available client-side — no extra server requests are needed.

### Anti-Pattern 4: Using function-based Cytoscape style values

**What people do:** Pass JavaScript functions as style property values to dynamically compute colors or sizes per-element.

**Why it's wrong:** Function values in Cytoscape.js stylesheet are the most expensive style type — they are called on every style recalculation. For a network with frequent interactions this causes visible lag.

**Do this instead:** Pre-encode the values in node data (e.g., `data.log2FC`, `data.significant`) and use `mapData()` or `data()` style accessors, which are computed once and cached.

---

## Integration Points

### External Services

| Service | Integration Pattern | Notes |
|---------|---------------------|-------|
| molAOP Builder `/api/v1/mappings` | HTTP GET via `BuilderAPIClient`, paginated, with retry adapter | Builder must have CORS enabled; pagination parameters need confirming against Builder API docs |
| molAOP Builder `/ke_genes/<ke_id>` | HTTP GET per KE (or bulk endpoint if available) | Per-KE endpoint may require N calls; bulk mappings endpoint preferred |
| molAOP Builder `/get_aop_kes/<aop_id>` | HTTP GET for AOP structure | Currently loaded from `aop_ke_map.csv` and `aop_ker_edges.csv`; these could also migrate to API |

### Internal Boundaries

| Boundary | Communication | Notes |
|----------|---------------|-------|
| `app.py` ↔ `helpers.py` | Direct function call, passes CSV paths | After change: also passes `api_client` instance or reads from config |
| `helpers.py` ↔ `api_client.py` | Direct method call: `BuilderAPIClient.get_ke_gene_mappings()` | api_client.py raises `BuilderAPIError`; helpers catches and falls back |
| `helpers.py` ↔ `cache_manager.py` | `cached_data_loader()` wraps `load_reference_sets()` | Cache contract unchanged; same TTL, same key |
| templates ↔ `app.py` | Jinja2 `render_template()` keyword args | Variable names must not change during UI overhaul |
| results.html ↔ Cytoscape.js | JSON blobs embedded in `<script>` tags (`network_json`, `ke_gene_json`) | Progressive loading reads from same data; no new endpoints needed |

---

## Scaling Considerations

| Scale | Architecture Adjustments |
|-------|--------------------------|
| Current (research tool, <100 concurrent users) | Monolith is correct. In-memory cache with Lock is fine. SQLite for experiment history is fine. |
| Moderate (100-1K users, platform adoption) | Replace `cache_manager.py` with Redis; replace SQLite with PostgreSQL. No architectural changes otherwise. |
| Large (1K+ users) | Worker queue (Celery) for long-running enrichment analysis; async HTTP client (httpx async) for Builder API calls; CDN for static assets. |

### Scaling Priorities (Current Scope)

1. **First bottleneck:** Builder API latency during reference set loading. **Fix:** 1-hour TTL cache already in place. Ensure cache is warm at startup (existing `load_cached_reference_sets()` call in `app.py`).
2. **Second bottleneck:** Cytoscape.js canvas rendering for large organ networks. **Fix:** Progressive gene layer + `hideEdgesOnViewport` + `haystack` edge style (Pattern 3).

---

## Sources

- Cytoscape.js performance documentation: [Performance Optimization | DeepWiki](https://deepwiki.com/cytoscape/cytoscape.js/8-performance-optimization) — HIGH confidence (official docs mirrored)
- Cytoscape.js official performance md: [cytoscape.js/documentation/md/performance.md](https://github.com/cytoscape/cytoscape.js/blob/master/documentation/md/performance.md) — HIGH confidence (official repo)
- Cytoscape.js WebGL renderer preview: [WebGL Renderer Preview blog post, Jan 2025](https://blog.js.cytoscape.org/2025/01/13/webgl-preview/) — HIGH confidence (official blog); status: provisional/not production-ready as of 2025-01
- Flask service layer pattern: [Architecture Patterns with Python — cosmicpython.com](https://www.cosmicpython.com/book/chapter_04_service_layer.html) — HIGH confidence (authoritative book)
- Flask template inheritance: [Flask Documentation — Templating](https://flask.palletsprojects.com/en/stable/templating/) — HIGH confidence (official docs)
- Circuit breaker / retry pattern: [PyBreaker PyPI](https://pypi.org/project/pybreaker/), [Resilient APIs — Medium](https://medium.com/@fahimad/resilient-apis-retry-logic-circuit-breakers-and-fallback-mechanisms-cfd37f523f43) — MEDIUM confidence (multiple consistent sources)
- Existing codebase: directly inspected `app.py`, `helpers.py`, `services/`, `config.py`, `cache_manager.py`, `database.py`, `exceptions.py` — HIGH confidence

---
*Architecture research for: molAOP Analyser — API integration, UI modernisation, network performance*
*Researched: 2026-02-25*
