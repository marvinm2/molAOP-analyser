# KE Mapping API Reference

> Context document for building a tool that fetches Key Event (KE) mappings from the Molecular AOP Builder application.

> *Last verified against the live Builder on **2026-08-18** (Builder **v2.9.0**). Earlier
> revisions tracked v2.8.0 / Phase 34 (Assessment Metadata Schema Parity); those changes were
> additive and left `services/api_service.py` unaffected, since it uses `.get()` patterns and
> ignores unknown keys.*
>
> **Known drift, not yet corrected here** (found 2026-08-18, none of it affects this
> application — all are integrator-facing and belong upstream):
> - The `/dataset/*` and `/export/*` endpoint families documented below return **503** and
>   **500** respectively on the live service. Do not wire against them without checking.
> - `connection_type` is described as WikiPathways-only; the live `go-mappings` records carry
>   it too, while Reactome does not.
> - The `go-mappings` JSON has ~19 fields against the 7 shown, and its CSV inserts
>   `go_direction` at column 8 — ahead of `suggestion_score` — so positional parsing tuned to
>   the WikiPathways CSV reads the wrong column.
> - GMT exports now begin with a `#` provenance block (resource, export revision,
>   source fingerprint, generation timestamp). `parse_gmt_reference_sets` and
>   `parse_gmt_pathway_gene_map` skip it explicitly as of 2026-08-18.

## What This System Is

The Molecular AOP Builder maps **Key Events** (KEs) from the AOP-Wiki to:
- **WikiPathways** (molecular pathway diagrams) -- stored as KE-WP mappings
- **Gene Ontology Biological Processes** (GO BP terms) -- stored as KE-GO mappings

All mappings are curator-approved with confidence levels (High/Medium/Low) and stable UUIDs.

---

## Public API v1

**Base URL:** `{host}/api/v1`
**Authentication:** None (fully public)
**Rate Limit:** 100 requests/hour per IP
**CORS:** Enabled (`Access-Control-Allow-Origin: *`)
**Methods:** GET only
**Response formats:** JSON (default) or CSV (`?format=csv`)

---

### GET /api/v1/mappings

List all approved KE-WikiPathways mappings. Supports filtering and pagination.

**Query Parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ke_id` | string | -- | Filter by KE ID (e.g. `KE 55`) |
| `pathway_id` | string | -- | Filter by WikiPathways ID (e.g. `WP1234`) |
| `confidence_level` | string | -- | `High`, `Medium`, or `Low` (case-insensitive) |
| `aop_id` | string | -- | Filter to KEs belonging to an AOP (numeric or `AOP 123` format) |
| `page` | integer | 1 | Page number |
| `per_page` | integer | 50 | Results per page (max 200) |
| `format` | string | json | Set to `csv` for CSV download |

**JSON Response:**

```json
{
  "data": [
    {
      "uuid": "a1b2c3d4-e5f6-...",
      "ke_id": "KE 55",
      "ke_name": "Increase, CYP2E1 expression",
      "pathway_id": "WP1234",
      "pathway_title": "Xenobiotic metabolism",
      "confidence_level": "High",
      "connection_type": "causative",
      "assessment": {
        "relationship": "causative",
        "basis": "known",
        "specificity": "specific",
        "coverage": "complete",
        "version": "v2"
      },
      "provenance": {
        "suggestion_score": 0.87,
        "approved_by": "curator_username",
        "approved_at": "2025-02-23T10:30:00"
      }
    }
  ],
  "pagination": {
    "page": 1,
    "per_page": 50,
    "total": 1234,
    "total_pages": 25,
    "next": "https://host/api/v1/mappings?page=2&per_page=50",
    "prev": null
  }
}
```

**Assessment block (Phase 34, v2.8.0+):**

| Field | Type | Allowed values |
|---|---|---|
| `assessment.relationship` | string \| null | `causative`, `responsive`, `bidirectional`, `unclear`, `null` |
| `assessment.basis` | string \| null | `known`, `likely`, `possible`, `uncertain`, `null` |
| `assessment.specificity` | string \| null | `specific`, `includes`, `loose`, `null` |
| `assessment.coverage` | string \| null | `complete`, `keysteps`, `minor`, `null` |
| `assessment.version` | string | `v1`, `v2` |

See "v1 vs v2 assessment row semantics" below for legacy-row behaviour.

**CSV columns:** `uuid, ke_id, ke_name, pathway_id, pathway_title, confidence_level, suggestion_score, approved_by, approved_at, proposed_by, connection_type, ke_aop_context, ke_bio_level, proposed_relationship, proposed_basis, proposed_specificity, proposed_coverage, assessment_version`

(Phase 34 added the trailing five columns at the END of the existing list for column-positional consumer back-compat.)

---

### GET /api/v1/mappings/{uuid}

Get a single KE-WP mapping by its stable UUID.

**Response:**

```json
{
  "data": {
    "uuid": "a1b2c3d4-e5f6-...",
    "ke_id": "KE 55",
    "ke_name": "Increase, CYP2E1 expression",
    "pathway_id": "WP1234",
    "pathway_title": "Xenobiotic metabolism",
    "confidence_level": "High",
    "connection_type": "causative",
    "assessment": {
      "relationship": "causative",
      "basis": "known",
      "specificity": "specific",
      "coverage": "complete",
      "version": "v2"
    },
    "provenance": {
      "suggestion_score": 0.87,
      "approved_by": "curator_username",
      "approved_at": "2025-02-23T10:30:00"
    }
  }
}
```

**Errors:** `404` if UUID not found.

---

### v1 vs v2 assessment row semantics

> **Pre-v1.6 mappings** (created before the Phase 34 schema migration) emit `assessment.version = "v1"` and all four answer fields (`relationship`, `basis`, `specificity`, `coverage`) as `null`. Consumers should treat these as "legacy rows with derived `confidence_level` only" — the four-question rubric was added in v1.6 and was not asked of pre-existing rows.
>
> **v1.6+ mappings** (created via the wired mapper UI on the builder) emit `assessment.version = "v2"` and at least one of the four answer fields non-null. The model-layer rule (`_classify_assessment_version`) flips a mapping to `v2` when ANY of the four answers is non-NULL — partial submissions during the Phase 34 → Phase 37 transition window are still classified `v2`.
>
> The four assessment columns are also exposed at the top level of the CSV export (`proposed_relationship`, `proposed_basis`, `proposed_specificity`, `proposed_coverage`, `assessment_version`) — same value set as the nested JSON keys.

**Sibling endpoint:** `GET /api/v1/reactome-mappings` (and `/{uuid}`) emits the identical `assessment` nested object and the same five trailing CSV columns — sibling parity for the rubric data is a Phase 34 invariant. (NOTE: `connection_type` is a WP-only top-level field — the `ke_reactome_mappings` table has no such column, only `ke_reactome_proposals.proposed_connection_type` exists. The four-question rubric in `assessment.*` is the cross-resource normalised representation.) The analyser does not currently consume the Reactome endpoint, but third-party consumers (e.g. `clusterProfiler` pipelines, FAIR exports) can rely on rubric-shape equivalence.

---

### GET /api/v1/go-mappings

List all approved KE-Gene Ontology (Biological Process) mappings.

**Query Parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ke_id` | string | -- | Filter by KE ID |
| `go_term_id` | string | -- | Filter by GO term ID (e.g. `GO:0006468`) |
| `confidence_level` | string | -- | `High`, `Medium`, or `Low` |
| `page` | integer | 1 | Page number |
| `per_page` | integer | 50 | Results per page (max 200) |
| `format` | string | json | Set to `csv` for CSV download |

**JSON Response:**

```json
{
  "data": [
    {
      "uuid": "b2c3d4e5-...",
      "ke_id": "KE 55",
      "ke_name": "Increase, CYP2E1 expression",
      "go_term_id": "GO:0006468",
      "go_term_name": "protein phosphorylation",
      "go_namespace": "biological_process",
      "confidence_level": "High",
      "provenance": {
        "suggestion_score": 0.92,
        "approved_by": "curator_username",
        "approved_at": "2025-02-23T10:30:00"
      }
    }
  ],
  "pagination": { ... }
}
```

**CSV columns:** `uuid, ke_id, ke_name, go_term_id, go_term_name, go_namespace, confidence_level, suggestion_score, approved_by, approved_at`

---

### GET /api/v1/go-mappings/{uuid}

Get a single KE-GO mapping by UUID. Same shape as list response `data[]` items.

---

### GET /api/v1/aops

List the AOPs the Builder knows about, each with its Key Event count and how many
of those KEs carry approved mappings, split per resource. Added 2026-07-21
(builder #207) so the analyser's AOP picker can be populated from the Builder
rather than only from AOP-Wiki SPARQL.

**Query params** (all optional):

| Param | Meaning |
|---|---|
| `mapped_only` | `true` returns only AOPs with at least one mapped KE |
| `q` | case-insensitive substring filter over `aop_id` and `aop_title` |
| `page` / `per_page` | pagination (default 50, max 200) |
| `format=csv` | CSV instead of JSON |

**Response** — sorted by `mapped_ke_count` descending:

```json
{
  "data": [
    {
      "aop_id": "AOP 625",
      "aop_title": "Increased 11\u03b2-Hydroxysteroid dehydrogenase type 1 ...",
      "ke_count": 18,
      "mapped_ke_count": 18,
      "wikipathways_ke_count": 18,
      "go_ke_count": 2,
      "reactome_ke_count": 2
    }
  ],
  "pagination": { "page": 1, "per_page": 50, "total": 274, "total_pages": 6, "next": "...", "prev": null }
}
```

> **`aop_id` uses a space, not a colon** — `"AOP 625"`, matching `ke_aop_context`
> on mapping records and `ke_id` (`"KE 1790"`). The analyser normalises to
> `"AOP:625"` in `services/aop_discovery_service.fetch_aops_from_builder`.

> **AOP membership is a snapshot, not live SPARQL.** It comes from the Builder's
> precomputed `data/ke_aop_membership.json`, the same source behind
> `ke_aop_context`. It went four months without a refresh before 2026-07-21, during
> which it undercounted KEs (AOP 625 reported 15 of 18) and omitted AOPs entirely
> (AOP 638). This is why the analyser uses the endpoint as a **fallback tier** and
> keeps the live SPARQL cross-reference as its primary source of coverage counts.
> If the numbers here disagree with AOP-Wiki, the snapshot is due a refresh
> (`scripts/precompute_ke_aop_membership.py`, then copy to the `/app/data` mount).

---

## Bulk Export Endpoints

These are on the main app (not under `/api/v1`), no auth required.

### GET /exports/gmt/ke-wp

Download KE-WP mappings as a GMT (Gene Matrix Transposed) file. The Builder
resolves each mapped pathway/term to its **member genes**, so the GMT contains
gene symbols (not pathway IDs) — ready for `fgsea` / `clusterProfiler`.

**Query params:** `?min_confidence=High|Medium|Low` (optional **threshold**)

`min_confidence` is a minimum, not a tier selector: `Medium` returns medium- and
high-confidence mappings, `High` returns high only, and `Low` is equivalent to
omitting the parameter. Mappings with no recorded confidence are always
included, so the filter can never silently empty an export.

> **Builder behaviour changed on 2026-07-22 (builder#206).** Before that date the
> parameter selected a *single* tier, so `min_confidence=Medium` returned medium
> only and dropped every high-confidence mapping — the opposite of what the name
> promises. This is why the Analyser's minimum-confidence control could not
> simply be forwarded to these endpoints (analyser#71): passing it through would
> have discarded the best-evidenced mappings.
>
> Two things to get right when wiring it up:
> - The Analyser's vocabulary is `("all", "medium", "high")`; the Builder's
>   whitelist is `{high, medium, low}`. `all` returns **HTTP 400** — translate it
>   to an omitted parameter rather than forwarding it verbatim.
> - `min_confidence` must be part of the reference-set cache key, or a run at one
>   threshold will serve gene sets cached at another.

**Response:** `text/plain` GMT file. One gene set per line, tab-separated:
```
<KE_id>_<KE_name>_<pathway_id>\t<pathway_title>\t<gene_1>\t<gene_2>\t...
```
e.g. `KE177_Increase_Mitochondrial_dysfunction_WP5241\tMitochondrial beta-oxidation\tECHS1\tACSL3\t...`.
A single KE typically spans several lines (one per mapped pathway/term); union
the gene columns across lines that share the same `KE<id>` token to obtain the
KE's full gene set.

### GET /exports/gmt/ke-go

Download KE-GO (GO Biological Process) mappings as GMT. Same params and
gene-level format as `ke-wp` (the descriptor column carries the `GO:nnnnnnn`
term id in place of the WP id).

### GET /exports/gmt/ke-reactome

Download KE-Reactome mappings as GMT. Same params and gene-level format as
`ke-wp` (descriptor column carries the `R-HSA-nnnnnn` pathway id).

**Analyser consumption:** the Molecular AOP Analyser builds its KE→gene
reference sets per resource. WikiPathways keeps its dedicated CSV-backed
pipeline; **GO BP and Reactome gene sets are sourced from `/exports/gmt/ke-go`
and `/exports/gmt/ke-reactome`** respectively (issue #55), parsed by
`services/api_service.py::fetch_gmt_reference_sets`. Users pick which
resource(s) feed enrichment in the single and batch analysis forms.

### GET /exports/rdf/ke-wp

Download KE-WP mappings as RDF/Turtle. **Response:** `text/turtle`

### GET /exports/rdf/ke-go

Download KE-GO mappings as RDF/Turtle.

### GET /export/{format_name}

Export full dataset. `format_name` can be: `json`, `jsonld`, `excel`, `parquet`, `csv`.

**Query params:**
- `metadata=true|false` -- include metadata (default true)
- `statistics=true|false` -- include stats sheet in Excel (default true)
- `provenance=true|false` -- include provenance info (default true)
- `compression=snappy|gzip` -- Parquet compression (default snappy)

### GET /export/formats

List available export formats with descriptions.

---

## Dataset Metadata Endpoints

No auth required.

### GET /dataset/metadata

```json
{
  "version": "2.0.0",
  "title": "KE-WikiPathways Curated Mapping Dataset",
  "total_mappings": 1234,
  "total_kes": 456,
  "total_pathways": 789,
  "confidence_distribution": { "high": 456, "medium": 345, "low": 433 },
  "last_updated": "2025-02-25T10:30:00Z"
}
```

### GET /dataset/versions

Version history of the dataset.

### GET /dataset/citation?format=apa|bibtex

Get a citation string for the dataset.

---

## Data Model

### KE-WP Mapping Fields

| Field | Type | Description |
|---|---|---|
| `uuid` | string | Stable unique identifier |
| `ke_id` | string | Key Event ID, e.g. `KE 55` |
| `ke_name` | string | Key Event title from AOP-Wiki |
| `pathway_id` | string | WikiPathways ID, e.g. `WP1234` |
| `pathway_title` | string | Pathway name |
| `confidence_level` | string | `High`, `Medium`, or `Low` (rubric-derived) |
| `connection_type` | string/null | Legacy v1 connection classification (`causative` / `responsive` / `other` / null). Superseded by `assessment.relationship` for v2 rows but kept populated for back-compat. |
| `assessment.relationship` | string/null | Phase 34 v1.6+: KER stance — `causative` / `responsive` / `bidirectional` / `unclear` / null |
| `assessment.basis` | string/null | Phase 34 v1.6+: evidence basis — `known` / `likely` / `possible` / `uncertain` / null |
| `assessment.specificity` | string/null | Phase 34 v1.6+: KE-to-pathway specificity — `specific` / `includes` / `loose` / null |
| `assessment.coverage` | string/null | Phase 34 v1.6+: pathway coverage of the KE — `complete` / `keysteps` / `minor` / null |
| `assessment.version` | string | `v1` (legacy, pre-v1.6 rows) or `v2` (rubric-aware rows) |
| `provenance.suggestion_score` | float/null | Hybrid scoring confidence (0.0-1.0) when suggested by the system |
| `provenance.approved_by` | string/null | Curator who approved the mapping |
| `provenance.approved_at` | string/null | ISO timestamp of approval |

### KE-GO Mapping Fields

| Field | Type | Description |
|---|---|---|
| `uuid` | string | Stable unique identifier |
| `ke_id` | string | Key Event ID |
| `ke_name` | string | Key Event title |
| `go_term_id` | string | GO term ID, e.g. `GO:0006468` |
| `go_term_name` | string | GO term name |
| `go_namespace` | string | Always `biological_process` |
| `confidence_level` | string | `High`, `Medium`, or `Low` |
| `provenance` | object | Same structure as KE-WP |

---

## Error Handling

All errors return JSON:

```json
{
  "error": "Human-readable error message"
}
```

| Status | Meaning |
|---|---|
| 400 | Bad request (invalid filter params, bad AOP ID) |
| 404 | UUID not found |
| 429 | Rate limit exceeded (check `Retry-After` header) |
| 500 | Server error |

---

## Integration Patterns

### Fetch all high-confidence KE-WP mappings

```bash
curl "https://host/api/v1/mappings?confidence_level=High&per_page=200"
```

To paginate through all results, follow `pagination.next` until it is `null`.

> **Analyser note (issue #60):** the `confidence_level` query parameter is
> exact-match, so it cannot express "Medium and above". The analyser therefore
> fetches every mapping once (it needs the full list for the pathway picker
> anyway) and applies the minimum-confidence threshold client-side in
> `helpers.filter_records_by_confidence()`. Comparison is case-insensitive —
> the live API returns lowercase values (`"high"`) while this document shows
> them title-cased. Records whose `confidence_level` is absent, `null`, blank,
> or outside `{high, medium, low}` are **kept** (the filter no-ops), so the
> local `data/KE-WP.csv` fallback and the confidence-free GMT exports keep
> working unchanged.

### Fetch mappings for a specific KE

```bash
curl "https://host/api/v1/mappings?ke_id=KE+55"
curl "https://host/api/v1/go-mappings?ke_id=KE+55"
```

### Fetch all mappings for an AOP

```bash
curl "https://host/api/v1/mappings?aop_id=23"
```

This resolves the AOP to its constituent KEs via AOP-Wiki SPARQL, then returns all mappings for those KEs.

### Bulk download as CSV

```bash
curl "https://host/api/v1/mappings?format=csv&per_page=200" -o mappings.csv
curl "https://host/api/v1/go-mappings?format=csv&per_page=200" -o go_mappings.csv
```

### Full dataset export

```bash
curl "https://host/export/json" -o dataset.json
curl "https://host/exports/rdf/ke-wp" -o ke_wp.ttl
```

---

## Key Identifiers

| Entity | Format | Example | Source |
|---|---|---|---|
| Key Event | `KE {number}` | `KE 55` | AOP-Wiki |
| WikiPathway | `WP{number}` | `WP1234` | WikiPathways |
| GO Term | `GO:{7 digits}` | `GO:0006468` | Gene Ontology |
| AOP | `{number}` or `AOP {number}` | `23` or `AOP 23` | AOP-Wiki |
| Mapping UUID | UUID v4 | `a1b2c3d4-e5f6-47a8-9b0c-d1e2f3a4b5c6` | Generated at creation |

---

## Rate Limits

The public API allows **100 requests per hour** per IP address. When exceeded, responses return `429` with a `Retry-After` header (seconds until reset).

For bulk data access, prefer the export endpoints (`/export/json`, `/exports/gmt/*`) which return complete datasets in a single request.
