# Bundled data: origin, licence and provenance

The repository's `LICENSE` (GPL-2.0-only) covers the **application source**. It does not
cover the files in this directory, which come from third parties and carry their own terms.
This file records, per file, where it came from and what licence it is under.

## Reference data

| File | Origin | Licence |
|---|---|---|
| `wikipathways_hsa_20240410.xgmml` | [WikiPathways](https://www.wikipathways.org/) human pathway collection, release **2024-04-10** (the date is in the filename) | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) — public domain dedication |
| `edges_wpid_to_gene.csv` | Derived from the XGMML above by `parse_xgmml_to_dataframe.py` | CC0 1.0, as a derivative of CC0 material |
| `node_attributes.csv` | Derived from the XGMML above | CC0 1.0, as above |
| `KE-WP.csv` | Snapshot of approved KE → WikiPathways mappings from the [molAOP Builder](https://molaop-builder.vhp4safety.nl) | CC0 1.0, except Key Event titles where present, which are reproduced from AOP-Wiki under CC BY-SA 4.0 |
| `aop_ke_map.csv`, `aop_ker_edges.csv`, `ke_metadata.csv` | AOP structure and Key Event metadata from [AOP-Wiki](https://aopwiki.org/) | Wiki content is **CC BY-SA 4.0** by default; attribute AOP-Wiki when re-using. Individual AOPs may be All Rights Reserved for a limited development period |

`KE-WP.csv` is a **fallback snapshot**, not the live source. The application fetches mappings
from the Builder's REST API and only falls back to this file when the API and the disk cache
are both unavailable; the provenance line on the results page records which was used. It is
therefore older than the live data by however long it has been since it was refreshed.

Two of these are pinned snapshots that drift from their upstreams and are known to:
`edges_wpid_to_gene.csv` tops out at `WP5452`, so any pathway curated since resolves to no
genes through this file alone — which is why the application prefers the Builder's GMT for
WikiPathways gene membership (see issue #79 and the note in `CLAUDE.md`).

## Demo expression datasets

`GSE90122_SR12813.tsv` and `GSE90122_TO90137.tsv` are differential-expression results derived
from **NCBI GEO series [GSE90122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90122)**:

> *Analysis of gene expression in human primary hepatocytes (HPPs) in response to the
> treatment of hPXR antagonist (SPA70) and agonists (rifampicin, SR12813, and TO90137)* —
> Homo sapiens, 24 samples, platform GPL14550, published 2017-10-02.

Each file corresponds to one PXR agonist named in its filename. The column set —
`logFC`, `t`, `B`, `P.Value`, `adj.P.Val` alongside the platform's probe annotation
(`SPOT_ID`, `GB_ACC`, `GENE_SYMBOL`, `GENE_NAME`, `UNIGENE_ID`, `SEQUENCE`) — is the standard
output of **GEO2R**, NCBI's limma-based web front end for GEO.

**The exact contrast and the limma settings behind these two files were not recorded at the
time and cannot now be reconstructed.** The presumption from the filenames is treated versus
vehicle control for the named agonist, but that is inference from the filename rather than a
recorded parameter, and the sample selection, any filtering, and the multiple-testing settings
are unknown. Anyone needing an exactly reproducible contrast should re-derive it from GSE90122
with GEO2R or limma directly and record the parameters; these two files are shipped as **demo
input for exercising the interface**, not as a reference analysis, and no result computed from
them should be reported as a finding about PXR agonism.

### Re-use terms

GEO does not apply a single blanket licence to submitted series. Data in GEO are generally
usable for research without restriction, and NCBI asks that the originating submitters and
the accession be credited; the underlying study should be cited when these data inform a
result. See [NCBI's data usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/).
The derived files here retain that expectation — cite GSE90122.

## If you are re-using this repository's data

- Cite **AOP-Wiki**, **WikiPathways** and **GSE90122** as appropriate to which files you
  touched.
- For the curated KE → pathway mappings, prefer the citable Zenodo deposit over this snapshot:
  [10.5281/zenodo.20184643](https://doi.org/10.5281/zenodo.20184643).
- Do not treat the GPL-2.0 in the repository root as covering any of it.
