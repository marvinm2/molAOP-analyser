# Requirements: molAOP Analyser

**Defined:** 2026-02-25
**Core Value:** Researchers and regulators can quickly determine which Key Events in a molecular AOP are activated by their transcriptomic data

## v1 Requirements

### API Integration

- [x] **API-01**: App fetches KE-gene set mappings from Builder `/api/v1/mappings` endpoint instead of static CSVs
- [x] **API-02**: App falls back gracefully to local CSV files when Builder API is unreachable
- [x] **API-03**: API responses are cached to avoid repeated calls during a session

### UI/UX

- [x] **UI-01**: All pages use VHP4Safety brand palette (#29235C, #E6007E, #307BBF) and Inter typography
- [x] **UI-02**: Sticky navigation header matching molAOP Builder style with app branding
- [x] **UI-03**: Results page visual hierarchy improved (consistent card layouts, spacing, typography)
- [x] **UI-04**: Upload/analysis page restyled to match brand

### Network Visualisation

- [x] **NET-01**: Large AOP networks render without browser slowdown (profiled and optimised)
- [x] **NET-02**: PNG export button visible in results page for publication screenshots

### Export

- [x] **EXP-01**: User can download enrichment results table as CSV

### Documentation

- [x] **DOC-01**: Documentation page with usage guide, input format spec, and statistical method explanation

## v2 Requirements

### API Integration

- **API-04**: Dynamic AOP list fetched from Builder instead of hardcoded config

### UI/UX

- **UI-05**: Responsive layout for tablet/laptop screen compatibility

### Network Visualisation

- **NET-03**: Progressive gene node rendering (KE nodes first, expand genes on click)

### Results

- **RES-01**: Shareable result URLs via unique experiment IDs

## Out of Scope

| Feature | Reason |
|---------|--------|
| AI/LLM analysis panel | Hallucination risk unacceptable for regulatory tool; maintenance burden outweighs benefit |
| OAuth/user authentication | Not needed for single-analyst tool; SQLite persistence covers reproducibility |
| GSEA (ranked list enrichment) | Different statistical assumptions; ORA is appropriate for binary significant/non-significant lists |
| Multi-omics input (proteomics, metabolomics) | No KE-gene mappings exist for non-transcriptomic identifiers at this maturity level |
| Full pathway database browsing (KEGG/Reactome) | Dilutes AOP-specific value proposition; users have g:Profiler/Enrichr for general analysis |
| Command-line/API batch analysis | Tool optimised for interactive regulatory review; Docker covers scripting needs |
| Real-time collaboration | WebSocket complexity disproportionate; shareable URLs cover most use cases |
| Report generation improvements | Not a priority for this round per user decision |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| API-01 | Phase 1 | Complete |
| API-02 | Phase 1 | Complete |
| API-03 | Phase 1 | Complete |
| UI-01 | Phase 2 | Complete |
| UI-02 | Phase 2 | Complete |
| UI-03 | Phase 2 | Complete |
| UI-04 | Phase 2 | Complete |
| EXP-01 | Phase 2 | Complete |
| DOC-01 | Phase 2 | Complete |
| NET-01 | Phase 3 | Complete |
| NET-02 | Phase 3 | Complete |

**Coverage:**
- v1 requirements: 11 total
- Mapped to phases: 11
- Unmapped: 0

---
*Requirements defined: 2026-02-25*
*Last updated: 2026-02-25 after roadmap creation*
