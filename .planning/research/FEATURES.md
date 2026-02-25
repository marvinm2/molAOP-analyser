# Feature Research

**Domain:** Molecular AOP transcriptomics enrichment analysis web tool
**Researched:** 2026-02-25
**Confidence:** MEDIUM-HIGH (ecosystem surveyed via web search + official sources; some competitor analysis from published papers + live tools)

---

## Context: What Already Exists

The tool is in an enhancement round, not a greenfield build. Core analysis capabilities are already live:

**Already built (do not re-implement):**
- File upload with format detection (CSV/TSV/TXT)
- Auto-detection of gene ID, log2FC, p-value columns with confidence scoring
- Interactive log2FC and p-value threshold sliders
- Fisher's exact test ORA with Benjamini-Hochberg FDR
- Cytoscape.js interactive AOP network
- Enrichment results table with formatted p-values and FDR
- Volcano plot
- Experiment metadata capture and SQLite persistence
- PDF/HTML report generation
- CSRF protection and input validation
- Demo datasets (PXR agonist studies, 42 cisplatin datasets)

**Enhancement targets for this round:**
1. Live API integration with molAOP Builder (replace static CSV reference data)
2. UI/UX overhaul — VHP4Safety brand, Inter font, navigation component
3. Network visualisation improvements — performance and interactivity
4. Documentation page

---

## Feature Landscape

### Table Stakes (Users Expect These)

Features that, if missing, make the tool feel incomplete or broken. Users from comparable tools (Enrichr, g:Profiler, AOP-networkFinder, EnrichmentMap) will arrive with these expectations.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Stable, current reference data | KE-gene mappings go stale; researchers notice when pathways don't match AOP-Wiki | MEDIUM | **This is the API integration milestone** — live Builder API replaces static CSVs |
| Clear visual hierarchy in results | All major enrichment tools (Enrichr, g:Profiler) separate significant from non-significant hits visually | LOW | Color-coded significance in table + network already exists; review completeness |
| Network node/edge count within visible range | Tools with >300 nodes become unreadable; researchers leave | MEDIUM | Filter genes by expression data already implemented; verify performance at scale |
| Export results as table (CSV/TSV) | Every benchmarked tool (Enrichr, g:Profiler, AOP-networkFinder) provides table download | LOW | Confirm enrichment table is downloadable; trivial to add if not |
| Export network as image | AOP-networkFinder, EnrichmentMap, AOP-helpFinder all offer PNG/SVG export for publications | LOW | Cytoscape.js has built-in PNG export; verify it's exposed in UI |
| Stable URLs / shareable results | EnrichmentMap:RNASeq offers unique shareable URLs; researchers collaborate and cite analyses | MEDIUM | SQLite persistence exists; question is whether results are URL-addressable |
| Consistent branding with platform | VHP4Safety users navigate between tools; visual discontinuity signals "different product" | MEDIUM | **This is the UI overhaul milestone** — Inter font, brand palette, nav component |
| Basic usage documentation | Any new user arriving from VHP4Safety platform needs onboarding; Enrichr, Reactome all provide inline help | LOW | **This is the documentation page milestone** |
| Responsive layout on tablet | Researchers use lab iPads and conference laptops; broken layouts lose trust | LOW | Test existing layout; CSS fixes likely minor |

### Differentiators (Competitive Advantage)

Features that distinguish this tool from generic enrichment analysis platforms. These are where the tool earns its place in the VHP4Safety ecosystem.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| AOP-contextualised network (not just pathway lists) | Enrichr and g:Profiler return flat lists or gene-pathway networks; this tool shows KE causal chains as defined in AOP-Wiki — unique framing for regulatory toxicology | LOW (already built) | Preserve and highlight this distinction in documentation and UI |
| Live Builder API integration | Reference data is always current with AOP-Wiki state; no stale static files; first in class for molAOP tooling on VHP4Safety | MEDIUM | Replace helpers.py static CSV loading with API calls to Builder `/api/v1/mappings`, `/ke_genes/<ke_id>`, `/get_aop_kes/<aop_id>` |
| KE enrichment score displayed on AOP network topology | AOP-networkFinder shows topology without expression data; this tool overlays significance and expression direction on the pathway structure — mechanistically meaningful | MEDIUM (partially built) | Verify KE node coloring by enrichment p-value is clear; confirm gene expression color scale is readable |
| Multi-dataset cisplatin dose-response series | 42 datasets covering 4-72hr × 0.1-50μM enables temporal and dose comparisons; no other tool ships with comparable ready-to-use regulatory toxicology data | LOW (already built) | Document this capability prominently; it's a genuine differentiator for kidney AOP users |
| Case-study pre-configured AOPs | Liver steatosis, kidney failure pre-selected with one click; reduces time-to-insight for VHP4Safety case study workflows | LOW (already built) | Ensure AOP selector is prominent in UI overhaul |
| Experiment provenance and reproducibility | SQLite persistence with metadata (stressor, dosing, owner, dataset ID) supports regulatory reproducibility requirements better than tools that discard session state | MEDIUM (already built) | Surface experiment history in UI; consider export of experiment record |
| Integrated volcano plot | Most enrichment tools show enrichment table only; contextualising DE genes alongside KE significance in one view aids interpretation | LOW (already built) | Ensure volcano plot is not hidden; verify it loads reliably |
| Network performance at AOP scale | AOP-networkFinder and AOP-helpFinder use Cytoscape.js but struggle at large networks; targeted performance work for gene-rich AOPs is a concrete improvement | MEDIUM | **This is the network performance milestone** — see dependency notes |

### Anti-Features (Commonly Requested, Often Problematic)

| Anti-Feature | Why Requested | Why Problematic | Alternative |
|--------------|---------------|-----------------|-------------|
| Full AI/LLM analysis panel (like AOP Net) | AOP Net (August 2025 preprint) added LLM-powered pathway interpretation; impressive demo | Requires LLM API key management, cost, latency, hallucination risk; AOP Net is still preprint — not production-validated; adds significant maintenance surface | Add structured tooltips with curated KE descriptions from AOP-Wiki metadata instead; link to AOP-Wiki KE pages for authoritative context |
| OAuth/user authentication and accounts | Researchers want to save and revisit sessions; seems obvious | Adds security surface, infrastructure complexity, password management; VHP4Safety platform authentication is out of scope per PROJECT.md | SQLite persistence already solves reproducibility; consider session-based result URLs as lower-complexity alternative |
| Real-time collaboration (multi-user sessions) | Common request for lab sharing | WebSocket infrastructure, session locking, conflict resolution are disproportionate for a single-analyst tool; most analysis is solo or sequential | Shareable result URLs achieve 90% of the use case without architectural complexity |
| GSEA (ranked list enrichment) in addition to ORA | Researchers using tools like g:Profiler expect GSEA mode for expression-ranked datasets | GSEA requires different statistical assumptions, different reference gene sets, different result interpretation; mixing ORA and GSEA in one UI confuses users and doubles the statistical implementation surface | Keep ORA with Fisher's exact — it is statistically appropriate for binary significant/non-significant gene lists; document why ORA was chosen for AOP-KE enrichment |
| Multi-omics input (proteomics, metabolomics) | Omics integration is an active VHP4Safety research area | Requires separate gene-ID mapping layers, different normalization assumptions, different background set definitions; scope expands dramatically; no clear KE-gene mappings exist for non-transcriptomic identifiers at this maturity level | Stay transcriptomics-only; document input format clearly; link to other VHP4Safety multi-omics tools |
| Full pathway database browsing (Reactome/KEGG-style) | Users of g:Profiler and DAVID expect broad pathway library selection | The value of this tool is AOP-specific framing; adding GO/KEGG/Reactome databases dilutes AOP focus and creates a weaker general enrichment tool competing against Enrichr (100M+ queries) | Keep AOP/KE focus; this is the differentiator; direct users to g:Profiler or Enrichr for general pathway analysis |
| Command-line / API endpoint for batch analysis | Bioinformaticians want to script analyses | Adds REST API surface to maintain; security review needed; the tool is optimised for interactive regulatory review, not pipeline integration | Docker availability (already implemented) covers the power-user scripting need; document Docker usage |

---

## Feature Dependencies

```
[Live Builder API Integration]
    └──enables──> [Always-current KE-gene reference data]
    └──enables──> [Dynamic AOP list from Builder]
    └──requires──> [Builder service availability / fallback handling]

[Network Performance Improvements]
    └──enhances──> [Large AOP network rendering]
    └──requires──> [Understanding of current bottlenecks — profile first]
    └──enhances──> [Gene-node filtering already built — verify it works at scale]

[UI/UX Overhaul]
    └──requires──> [Brand palette already defined (CLAUDE.md)]
    └──requires──> [Navigation component pattern from molAOP Builder]
    └──enhances──> [Documentation page — consistent layout]
    └──blocks──> [Documentation page styling — share stylesheet]

[Documentation Page]
    └──blockedBy──> [UI/UX Overhaul — needs consistent nav/header component]
    └──requires──> [Content: workflow description, data format spec, statistical method explanation]

[Shareable Result URLs]
    └──requires──> [SQLite persistence (already built)]
    └──enhances──> [Experiment reproducibility differentiator]
```

### Dependency Notes

- **Builder API requires fallback**: If Builder is down, the analyser must gracefully degrade — either cached data or clear error messaging. Don't block analysis on API uptime.
- **UI overhaul before documentation page**: Documentation page should use the same nav/header component; build shared component first, apply to documentation second.
- **Network performance requires profiling before building**: The bottleneck could be layout algorithm, DOM element count, or Cytoscape.js render loop. Profile actual large AOP datasets before implementing solutions — avoid premature optimization.
- **API integration and UI overhaul are independent**: Can be developed in parallel on separate branches.

---

## MVP Definition for This Milestone

This is an enhancement round on a working tool. MVP = the three active requirements in PROJECT.md, plus documentation.

### Launch With (this milestone)

- [ ] **Live Builder API integration** — KE-gene sets fetched from Builder `/api/v1/mappings` and related endpoints, with graceful fallback to local static files when Builder is unavailable. This is the highest-value item: it ensures the tool stays current with AOP-Wiki without manual data refreshes.
- [ ] **UI/UX overhaul** — Inter typography, VHP4Safety brand palette applied consistently, sticky navigation matching Builder style, visual improvements to results page layout. Non-functional but critical for platform trust and adoption.
- [ ] **Network interactivity improvements** — At minimum: verify node filtering performance on largest available test dataset (kidney organ network, ~hundreds of nodes); apply any quick-win improvements (element batching, style caching, lazy gene node rendering). Full rewrite of network service is overkill unless profiling shows it's necessary.
- [ ] **Documentation page** — Single page covering: what the tool does, expected input format (gene symbol, log2FC, p-value columns), statistical method (ORA / Fisher's exact / FDR), how to interpret results, links to AOP-Wiki for KE details.

### Add After Validation (v1.x)

- [ ] **Shareable result URLs** — Requires verifying that SQLite experiment records contain enough data to reconstruct results; then expose a `/results/<experiment_id>` route. Trigger: user requests from VHP4Safety community.
- [ ] **Enrichment table CSV export** — If not already implemented; trivial addition once confirmed missing. Trigger: first user complaint about copying data.
- [ ] **Network PNG/SVG export button** — Surface existing Cytoscape.js export capability with a clear UI button if not already prominent.

### Future Consideration (v2+)

- [ ] **AOP-Wiki SPARQL live query** — Direct SPARQL endpoint queries to AOP-Wiki RDF for KE metadata and pathway structure, bypassing Builder API as intermediary. Dependency: Builder API covers this adequately for now.
- [ ] **Comparative analysis view** — Side-by-side AOP networks for two datasets (e.g., different cisplatin timepoints). High user value for dose-response studies but high implementation complexity.
- [ ] **Structured KE metadata tooltips** — Pull KE descriptions, evidence levels, and links from AOP-Wiki API and display on node hover. Currently nodes link to AOP-Wiki externally; inline context would be richer.

---

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Live Builder API integration | HIGH | MEDIUM | P1 |
| UI/UX overhaul (brand + nav) | HIGH | MEDIUM | P1 |
| Network performance at scale | HIGH | MEDIUM | P1 |
| Documentation page | MEDIUM | LOW | P1 |
| Enrichment table CSV export | MEDIUM | LOW | P2 |
| Network PNG/SVG export button | MEDIUM | LOW | P2 |
| Shareable result URLs | MEDIUM | MEDIUM | P2 |
| KE metadata tooltips from API | MEDIUM | MEDIUM | P3 |
| Comparative dataset view | HIGH | HIGH | P3 |
| AI/LLM analysis panel | LOW | HIGH | Defer indefinitely |

**Priority key:**
- P1: Must have for this milestone
- P2: Should have, add when P1 is stable
- P3: Nice to have, future consideration

---

## Competitor Feature Analysis

| Feature | AOP-networkFinder | AOP-helpFinder 3.0 | AOP Net (2025 preprint) | Enrichr | Our Approach |
|---------|-------------------|-------------------|------------------------|---------|--------------|
| AOP network visualization | Cytoscape.js export | Cytoscape.js embedded | Custom graph DB | None (flat list) | Cytoscape.js embedded, KE topology |
| Transcriptomic data overlay | No | No | Planned | Gene-pathway overlap only | Full expression coloring + enrichment significance |
| Live data from AOP-Wiki | SPARQL direct | SPARQL direct | Real-time sync | N/A | Via Builder API (abstraction layer) |
| Text mining / literature | No | Yes (PubMed) | Via LLM | No | Out of scope |
| Export formats | PNG, Cytoscape | JSON, Cytoscape | PNG, SVG, JSON, CSV, CX2 | CSV, API | PNG via Cytoscape.js, HTML/PDF reports |
| AI features | None | None | LLM Analysis Panel (7 modes) | None | Deliberately avoid — see anti-features |
| Reproducibility / sessions | Logs user actions | Per-user web session | Not specified | URL-based sharing | SQLite persistence + metadata |
| Pre-configured case studies | None | None | None | None | VHP4Safety AOPs (liver, kidney) |
| Demo datasets | None | None | None | None | 42 cisplatin + PXR datasets |

**Key observation**: No competitor combines AOP topology visualization with transcriptomic expression overlay AND pre-curated regulatory toxicology case studies. This combination is the unique positioning of molAOP Analyser.

---

## Sources

- [AOP-networkFinder paper, Bioinformatics Advances 2025](https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf007/7972740) — MEDIUM confidence (Oxford academic, peer-reviewed)
- [AOP Net preprint, bioRxiv August 2025](https://www.biorxiv.org/content/10.1101/2025.08.25.672239v1.full) — LOW confidence (preprint, not peer-reviewed)
- [AOP-helpFinder 3.0, Bioinformatics Oxford 2025](https://academic.oup.com/bioinformatics/article/41/7/btaf381/8177148) — MEDIUM confidence (peer-reviewed)
- [EnrichmentMap:RNASeq, Bioinformatics Advances 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12373637/) — HIGH confidence (peer-reviewed, PMC)
- [g:Profiler 2023 update, Nucleic Acids Research](https://academic.oup.com/nar/article/51/W1/W207/7152869) — HIGH confidence (Oxford NAR, official)
- [Enrichr platform, Maayan Lab](https://maayanlab.cloud/Enrichr/) — HIGH confidence (official tool)
- [WikiPathways 2024, Nucleic Acids Research](https://academic.oup.com/nar/article/52/D1/D679/7369835) — HIGH confidence (peer-reviewed)
- [VHP4Safety platform tools](https://platform.vhp4safety.nl/tools) — MEDIUM confidence (official platform, accessed 2026-02)
- [molAOP Analyser VHP4Safety cloud service](https://cloud.vhp4safety.nl/service/molaopanalyser.html) — HIGH confidence (official)
- Pathway enrichment visualization best practices from [PMC survey of network visualization tools](https://pmc.ncbi.nlm.nih.gov/articles/PMC5540468/) — MEDIUM confidence

---

*Feature research for: Molecular AOP Transcriptomics Enrichment Analysis Tool (enhancement round)*
*Researched: 2026-02-25*
