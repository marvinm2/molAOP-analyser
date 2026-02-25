# Stack Research

**Domain:** Flask/Python web application — API integration, modern UI, biological network visualisation
**Researched:** 2026-02-25
**Confidence:** MEDIUM-HIGH (versions verified via PyPI/official docs; rationale drawn from official sources and verified community patterns)

---

## Scope

This document covers only the **three new concerns** for this milestone:

1. Live API integration with molAOP Builder (`/api/v1/mappings`, `/ke_genes/<ke_id>`, `/get_aop_kes/<aop_id>`)
2. UI/UX overhaul — Inter font, VHP4Safety brand palette, modern Flask component patterns
3. Cytoscape.js network visualisation performance for dozens of KE nodes + hundreds of gene nodes

The existing stack (Flask, scipy, pandas, Cytoscape.js canvas renderer, SQLite/SQLAlchemy) is **not re-evaluated here**.

---

## Recommended Stack

### API Integration Layer

| Library | Version | Purpose | Why Recommended |
|---------|---------|---------|-----------------|
| `httpx` | 0.28.1 | HTTP client for molAOP Builder API calls | Supports both sync and async, HTTP/2, fully type-annotated, requests-compatible API makes migration trivial. Flask runs synchronously, so sync client is used now; async path available if app moves to Quart/ASGI later. |
| `tenacity` | 9.1.4 | Retry with exponential backoff for API calls | Builder API may be unavailable during startup or restarts. Tenacity's decorator-based API keeps retry logic out of business logic. Supports both sync and async. Far more capable than urllib3's Retry for application-level retry. |
| `flask-caching` | 2.3.1 | Cache API responses in memory | Builder's `/api/v1/mappings` data changes infrequently (KE-gene sets are stable between builder edits). A 5-minute SimpleCache TTL eliminates redundant API calls per-analysis without adding Redis infrastructure. |
| `flask-cors` | 6.0.2 | CORS headers if analyser serves its own API endpoints | Currently only consuming Builder API, not serving one. But if the analyser exposes results endpoints for the platform, CORS is needed. Minimal cost to add now. |

**Rationale for httpx over `requests`:** The existing `requirements.txt` uses no HTTP client (all data is local CSV files). When adding one, httpx is the forward-looking choice: it has async support (not needed now, but zero migration cost later), HTTP/2, and first-class type annotations. The API is nearly identical to `requests` so existing team knowledge transfers directly. [Source: httpx docs, PyPI — HIGH confidence]

**Rationale for tenacity over urllib3 Retry:** urllib3's `Retry` with `HTTPAdapter` works for simple cases but requires more boilerplate and lacks jitter, conditional retry logic, and async support. Tenacity handles all of these with a clean decorator API. [Source: PyPI, verified with tenacity docs — HIGH confidence]

**API integration pattern (sync):**

```python
# services/builder_api_service.py
import httpx
from tenacity import retry, stop_after_attempt, wait_exponential
from flask_caching import Cache

cache = Cache()  # configured in app factory

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=10))
def fetch_ke_gene_sets(builder_base_url: str) -> dict:
    """Fetch KE→gene mappings from molAOP Builder API with retry."""
    with httpx.Client(timeout=10.0) as client:
        response = client.get(f"{builder_base_url}/api/v1/mappings")
        response.raise_for_status()
        return response.json()

@cache.cached(timeout=300, key_prefix="ke_gene_sets")
def get_ke_gene_sets_cached(builder_base_url: str) -> dict:
    return fetch_ke_gene_sets(builder_base_url)
```

---

### UI/UX Layer

| Technology | Version | Purpose | Why Recommended |
|-----------|---------|---------|-----------------|
| Tailwind CSS | 4.0 (January 2025) | Utility-first CSS framework for consistent brand styling | Zero-config automatic content detection works without touching Jinja2 template paths. CSS-first `@theme` directive maps VHP4Safety brand colours as CSS variables, accessible in templates with no JavaScript. Incremental builds 8x faster than v3. |
| Inter font | via Google Fonts CDN or self-hosted | Brand-mandated typography | Used by VHP4Safety platform. Free, open-source, no licensing concern. Google Fonts CDN adds one HTTP request; self-hosting avoids GDPR concerns for EU deployments. |
| htmx | 2.0.8 | Partial template updates without full page reload | Enables the analysis results panel and network to update without rewriting Flask to an SPA. Zero build step — include via CDN script tag. Pairs naturally with Jinja2 SSR. No JavaScript framework needed. |
| `jinja2-fragments` | 1.11.0 | Render individual Jinja2 template blocks | Required for htmx partial rendering pattern in Flask. `render_block()` from `jinja2_fragments.flask` is a drop-in replacement for `render_template`. Eliminates duplicated template logic for partial vs full-page responses. |
| Alpine.js | 3.15.8 | Lightweight JS for interactive UI state | Handles dropdown menus, modal toggles, tooltip visibility — things that htmx doesn't control. 15KB minified+gzipped. Declarative `x-data` attributes keep JS co-located with markup. No build step. |

**Why Tailwind v4 over Bootstrap:** Bootstrap's component classes impose visual opinions that conflict with a custom brand system. Tailwind's utility approach lets the VHP4Safety palette govern all visual decisions without overriding defaults. Tailwind v4's `@theme` maps brand tokens directly in CSS, making color palette maintenance a one-file change. [Source: tailwindcss.com/blog/tailwindcss-v4 — HIGH confidence]

**Why htmx + Alpine.js over React/Vue:** The existing codebase is Flask-rendered Jinja2. Introducing a JavaScript SPA framework would require: a separate build pipeline, API endpoints for all data, and maintaining two rendering systems. htmx extends HTML semantics rather than replacing them. Alpine.js handles the residual interactivity that pure HTML cannot. This stack is used commercially in production Flask applications (confirmed via Hacker News + teclado.com case studies — MEDIUM confidence).

**Why NOT Vue or React:** Overengineered for a primarily server-rendered analysis tool with moderate interactivity. Adds build complexity (webpack/vite), JS bundle weight, and split mental model between Jinja2 and component templates. Flask is not an API-first framework.

**VHP4Safety brand token mapping in Tailwind v4:**

```css
/* static/css/input.css */
@import "tailwindcss";

@theme {
  --color-primary-dark: #29235C;
  --color-primary-magenta: #E6007E;
  --color-primary-blue: #307BBF;
  --color-light-blue: #009FE3;
  --color-orange: #EB5B25;
  --color-sky-blue: #93D5F6;
  --color-deep-magenta: #9A1C57;
  --color-teal: #45A6B2;
  --color-purple: #B81178;
  --color-dark-teal: #005A6C;
  --color-violet: #64358C;
  --font-sans: "Inter", system-ui, sans-serif;
}
```

---

### Network Visualisation Performance

| Technology | Version | Purpose | Why Recommended |
|-----------|---------|---------|-----------------|
| Cytoscape.js | 3.33.1 | Core graph library (existing) | Latest stable; do not upgrade plugins without testing. |
| `cytoscape-fcose` | 2.2.0 | Force-directed layout for biological networks | Official recommendation from Cytoscape.js docs: "the first layout to try for force-directed layouts." Spectral initialisation + force refinement gives better results than random-init CoSE. Best for KE-gene compound networks. [Source: js.cytoscape.org — HIGH confidence] |
| `cytoscape-dagre` | 2.5.0 | DAG-based hierarchical layout | AOP networks have a directed acyclic structure (KEs flow upstream→downstream). Dagre produces cleaner hierarchical layouts for AOP pathways than force-directed. Use dagre for the AOP skeleton (KE→KE edges) and fcose for the gene satellite layer. |
| WebGL renderer (built-in, experimental) | Cytoscape.js 3.31+ | GPU-accelerated rendering | Available now in 3.33.1 via `renderer: { name: 'webgl' }`. Verified performance: 1,200-node/16,000-edge networks go from ~20 FPS to 100+ FPS. For molAOP networks (dozens of KEs + hundreds of genes = hundreds of nodes, thousands of edges), this is the correct solution. [Source: blog.js.cytoscape.org/2025/01/13/webgl-preview/ — HIGH confidence] |

**WebGL renderer limitations to know:**
- Only `straight`, `haystack`, and `bezier` edge types supported (unsupported types fall back to bezier — acceptable)
- No dashed lines, no edge gradients, no multiple edge labels
- Only triangle arrow shapes
- Nodes are fully supported (all styles render via sprite sheet)
- If any of these limitations are blocking, fall back to canvas renderer with `hideEdgesOnViewport: true` and `textureOnViewport: true` as a partial optimisation

**Recommended rendering configuration:**

```javascript
const cy = cytoscape({
  container: document.getElementById('cy'),
  renderer: { name: 'webgl' },          // GPU acceleration (3.31+)
  hideEdgesOnViewport: true,             // Fallback if WebGL unavailable
  textureOnViewport: true,
  wheelSensitivity: 0.3,
  pixelRatio: 'auto',
  // ... elements, style, layout
});

// Batch all element additions
cy.batch(() => {
  cy.add(elements);
});
```

**Layout strategy for molAOP networks:**
- Apply `dagre` for the AOP backbone (KE nodes + KE→KE edges) to preserve pathway directionality
- Apply `fcose` for the full compound graph including gene nodes (if compound graph needed)
- The existing KE→gene thin gray lines styling is compatible with WebGL renderer (solid color, straight/bezier edges)

---

## Supporting Tools (Development)

| Tool | Purpose | Notes |
|------|---------|-------|
| `pytest-httpx` | Mock httpx calls in tests | Required when mocking Builder API in test suite — drop-in for `responses` library but httpx-specific |
| Tailwind CSS CLI | Build utility CSS | `npx @tailwindcss/cli` — no Node project required. Add `tailwindcss` build step to Makefile/Docker. |

---

## Installation

```bash
# Python — API integration
pip install httpx==0.28.1 tenacity==9.1.4 Flask-Caching==2.3.1 flask-cors==6.0.2 jinja2-fragments==1.11.0

# Python — test support
pip install pytest-httpx

# Node — Tailwind CSS build only (not a runtime dependency)
npm install -D tailwindcss@4 @tailwindcss/cli

# JS via CDN (include in base template, no npm required)
# htmx 2.0.8: https://cdn.jsdelivr.net/npm/htmx.org@2.0.8/dist/htmx.min.js
# Alpine.js 3.15.8: https://cdn.jsdelivr.net/npm/alpinejs@3.15.8/dist/cdn.min.js
# Cytoscape.js-fcose 2.2.0: https://cdn.jsdelivr.net/npm/cytoscape-fcose@2.2.0/cytoscape-fcose.min.js
# Cytoscape.js-dagre 2.5.0: https://cdn.jsdelivr.net/npm/cytoscape-dagre@2.5.0/cytoscape-dagre.js
```

---

## Alternatives Considered

| Recommended | Alternative | When to Use Alternative |
|-------------|-------------|-------------------------|
| `httpx` | `requests` | Only if the team has strong objection to adding a new dependency and no async path is anticipated. Requests is unmaintained for HTTP/2 and lacks async. |
| `tenacity` | urllib3 `Retry` + `HTTPAdapter` | If retry logic is simple (same error type, constant backoff) and adding a dependency is undesirable. urllib3 Retry has no jitter and no async support. |
| `flask-caching` SimpleCache | Redis + flask-caching | When the app is multi-process (gunicorn workers > 1). SimpleCache is per-process; Redis is shared. For single-worker deployments, SimpleCache is sufficient. |
| Tailwind CSS v4 | Bootstrap 5 | If brand customisation is minimal and the team needs pre-built components faster. Bootstrap imposes visual opinions that conflict with custom brand tokens. |
| htmx + Alpine.js | React or Vue | If the app needs a rich interactive SPA with complex client-side state (e.g., live editing). Not applicable here. |
| WebGL renderer | Canvas + `hideEdgesOnViewport` | If network size stays reliably under ~200 nodes/2,000 edges. WebGL is experimental; canvas fallback is stable. Use both: try WebGL, fall back to canvas on unsupported browsers. |
| `cytoscape-fcose` | `cytoscape-cose-bilkent` | fcose supersedes cose-bilkent. Only use cose-bilkent if fcose has a regression for a specific graph type. |
| `cytoscape-dagre` | `cytoscape-elk` | ELK (Eclipse Layout Kernel) produces better layouts for very complex hierarchies but requires a WASM or server-side layout step. Dagre is pure JS and sufficient for AOP-scale networks. |

---

## What NOT to Use

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| `aiohttp` as HTTP client | async-only; Flask is sync. Using aiohttp in sync Flask requires running an event loop manually, which is error-prone. | `httpx` (sync+async) |
| `retrying` library | No longer maintained (last release 2016). Tenacity is its actively maintained successor. | `tenacity` |
| React/Vue/Svelte | Requires build pipeline, splits rendering between server (Jinja2) and client (component tree). No benefit for this tool's interaction model. | htmx + Alpine.js |
| Bootstrap 5 | Generic component styling conflicts with VHP4Safety brand tokens. Tailwind utilities require no override layer. | Tailwind CSS v4 |
| `cytoscape-cose-bilkent` | Superseded by fcose. fcose is faster and produces equivalent aesthetics. | `cytoscape-fcose` |
| Cytoscape.js WebGL renderer in production without canvas fallback | Experimental as of 3.33.1; may have browser-specific issues. | Use `try/catch` on renderer init and fall back to canvas |
| Google Fonts CDN (for EU deployments) | May violate GDPR by sending visitor IP to Google. | Self-host Inter from `fontsource` or downloaded via `npm install @fontsource/inter` |

---

## Stack Patterns by Variant

**If molAOP Builder is unavailable (offline/startup):**
- Serve stale cached API data from `flask-caching`
- Fall back to local CSV files as secondary data source (existing `load_reference_sets()` in helpers.py)
- Log warning but do not crash — graceful degradation

**If network > 500 nodes (very large AOP with many gene sets):**
- Enable WebGL renderer
- Reduce edge rendering: hide KE→gene edges by default, show on node hover/click only
- Use `cy.batch()` for all programmatic element additions
- Consider virtualised edge rendering: only render edges for visible viewport

**If browser does not support WebGL:**
- Detect via `cytoscape({ renderer: { name: 'webgl' } })` failure in try/catch
- Fall back to canvas renderer with `hideEdgesOnViewport: true` and `textureOnViewport: true`

**If Tailwind CLI build step is undesirable in Docker:**
- Use Tailwind CDN play script for development only: `<script src="https://cdn.tailwindcss.com"></script>` with inline config
- Note: CDN play script is 100KB+ and unsuitable for production — use CLI build in Dockerfile

---

## Version Compatibility

| Package | Compatible With | Notes |
|---------|-----------------|-------|
| `httpx==0.28.1` | Python 3.8+ | Requires Python 3.8+; project already on 3.x |
| `tenacity==9.1.4` | Python 3.10+ | Requires Python 3.10+; verify Docker base image |
| `Flask-Caching==2.3.1` | Flask 2.x, 3.x | Works with existing Flask version |
| `jinja2-fragments==1.11.0` | Jinja2 3.x, Flask 2.x/3.x | Drop-in; no Flask version constraint |
| `flask-cors==6.0.2` | Flask 2.x, 3.x | Major version bump from 3.x; check for breaking changes if previously used |
| `cytoscape-fcose==2.2.0` | Cytoscape.js 3.x | Tested with 3.33.x; no known issues |
| `cytoscape-dagre==2.5.0` | Cytoscape.js 3.x | Stable, last released 2022 — no active development but no regressions |
| Tailwind CSS v4 | Node 18+ | Build-time only dependency; no Node in production container |
| htmx 2.0.8 | Any modern browser | IE11 not supported — consistent with Cytoscape.js requirement |
| Alpine.js 3.15.8 | Any modern browser | IE11 not supported |

---

## Sources

- PyPI — `httpx 0.28.1` [https://pypi.org/project/httpx/] — version verified — HIGH confidence
- PyPI — `tenacity 9.1.4` [https://pypi.org/project/tenacity/] — version verified — HIGH confidence
- PyPI — `Flask-Caching 2.3.1` [https://pypi.org/project/Flask-Caching/] — version verified — HIGH confidence
- PyPI — `flask-cors 6.0.2` [https://pypi.org/project/flask-cors/] — version verified — HIGH confidence
- PyPI — `jinja2-fragments 1.11.0` [https://pypi.org/project/jinja2-fragments/] — version and Flask support verified — HIGH confidence
- Tailwind CSS official blog — v4.0 release January 22 2025, CSS-first config, auto content detection [https://tailwindcss.com/blog/tailwindcss-v4] — HIGH confidence
- Cytoscape.js official site — version 3.33.1, performance options [https://js.cytoscape.org/] — HIGH confidence
- Cytoscape.js blog — WebGL renderer preview, 3.31 release, performance benchmarks [https://blog.js.cytoscape.org/2025/01/13/webgl-preview/] — HIGH confidence
- GitHub — Alpine.js v3.15.8 releases [https://github.com/alpinejs/alpine/releases] — HIGH confidence
- htmx.org — version 2.0.8 stable [https://htmx.org/posts/2024-06-17-htmx-2-0-0-is-released/] — HIGH confidence
- npm — `cytoscape-fcose 2.2.0` [https://www.npmjs.com/package/cytoscape-fcose] — version verified, maintenance status noted — MEDIUM confidence (package not actively updated but stable)
- GitHub — `cytoscape-dagre v2.5.0` [https://github.com/cytoscape/cytoscape.js-dagre] — version verified — MEDIUM confidence

---

*Stack research for: molAOP Analyser — API integration, UI overhaul, network visualisation performance*
*Researched: 2026-02-25*
