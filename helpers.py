import logging
import os

import pandas as pd

from config import DATA_DIR

logger = logging.getLogger(__name__)

# Issue #60: minimum KE-mapping confidence. Builder mappings carry a
# confidence_level of High/Medium/Low; the analyser lets users restrict the
# enrichment to better-evidenced mappings. "all" is the default and keeps the
# pre-#60 behaviour (no filtering).
VALID_MIN_CONFIDENCE = ("all", "medium", "high")
DEFAULT_MIN_CONFIDENCE = "all"

# Human-readable labels for results pages and reports.
MIN_CONFIDENCE_LABELS = {
    "all": "All mappings",
    "medium": "Medium and High only",
    "high": "High only",
}

# Issue #123: how the per-Key-Event gene sets of the selected resources are
# combined. 'union' (a gene counts if ANY selected source has it) is the default
# and the only behaviour before #123. 'intersection' needs every selected source;
# 'at_least_<N>' needs N of them, of which union and intersection are the N=1
# and N=all endpoints. Stored in that canonical spelling; NULL reads as union.
SOURCE_COMBINATION_UNION = "union"
SOURCE_COMBINATION_INTERSECTION = "intersection"
SOURCE_COMBINATION_AT_LEAST = "at_least"
VALID_SOURCE_COMBINATION_MODES = (
    SOURCE_COMBINATION_UNION,
    SOURCE_COMBINATION_INTERSECTION,
    SOURCE_COMBINATION_AT_LEAST,
)
DEFAULT_SOURCE_COMBINATION = SOURCE_COMBINATION_UNION
_AT_LEAST_PREFIX = "at_least_"

# Ordinal rank per confidence level; a mapping qualifies when its rank is at
# least the rank of the requested threshold.
_CONFIDENCE_RANKS = {"low": 1, "medium": 2, "high": 3}

# Record/column keys that may carry the confidence value. The Builder API uses
# ``confidence_level``; the capitalised spellings are tolerated so CSV exports
# of the same data work unchanged.
CONFIDENCE_KEYS = ("confidence_level", "Confidence_Level", "Confidence", "confidence")


class ReferenceSets(dict):
    """KE_ID -> gene set mapping that remembers what it could not resolve.

    A plain ``dict`` in every respect, plus an attribute:
    ``unresolved_ke_pathways``, the KE_ID -> sorted pathway-ID map of curated
    mappings whose gene membership no configured source could resolve.

    Issue #81: that information is produced deep inside
    :func:`load_reference_sets`, but it is needed at the far end of the call
    chain — in the enrichment backends, which decide whether an untestable Key
    Event is *uncurated* or *unresolvable*. Between the two sit callers that
    only ever hand the gene sets on (``api_service.fetch_reference_sets_from_api``,
    ``batch_service.run_batch``), so carrying it on the mapping itself keeps it
    from being dropped by an intermediary that has no reason to know about it.
    Consumers should read it through :func:`unresolved_ke_pathways_for`, which
    tolerates a plain dict from a cache written before this existed.

    Issue #123 adds a second attribute, ``no_shared_genes_kes``: the Key Events
    that are mapped and have genes, but whose selected sources did not share
    enough of them to survive the requested source combination. Read it with
    :func:`no_shared_genes_kes_for`. Always empty under union.
    """

    def __init__(self, *args, unresolved_ke_pathways=None, no_shared_genes_kes=None,
                 gene_universe=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.unresolved_ke_pathways = {
            str(ke): sorted(pathways)
            for ke, pathways in (unresolved_ke_pathways or {}).items()
            if pathways
        }
        self.no_shared_genes_kes = {str(ke) for ke in (no_shared_genes_kes or ())}
        # Every gene any loaded source lists, before a narrowing combination
        # (#123). None means "the sets themselves are the universe", which is
        # true under union. Read it through reference_gene_universe_for.
        self.gene_universe = (
            frozenset(str(g).strip().upper() for g in gene_universe)
            if gene_universe is not None else None
        )


def unresolved_ke_pathways_for(reference_sets):
    """Read the KE -> unresolvable-pathway-ID map off a reference-set mapping.

    Parameters
    ----------
    reference_sets : dict
        Any KE_ID -> gene set mapping. A :class:`ReferenceSets` carries the
        map; a plain dict (an older cache entry, a hand-built test fixture)
        does not.

    Returns
    -------
    dict[str, list[str]]
        The map, or ``{}`` when the mapping does not carry one.
    """
    return dict(getattr(reference_sets, 'unresolved_ke_pathways', None) or {})


def no_shared_genes_kes_for(reference_sets):
    """Read the Key Events emptied by the source combination (issue #123).

    Parameters
    ----------
    reference_sets : dict
        Any KE_ID -> gene set mapping. Only a :class:`ReferenceSets` built by a
        non-union combination carries any.

    Returns
    -------
    set[str]
        The KE IDs, or an empty set when the mapping carries none.
    """
    return set(getattr(reference_sets, 'no_shared_genes_kes', None) or ())


def reference_gene_universe_for(reference_sets):
    """Every gene the loaded sources list, whatever the source combination.

    The #69 identifier check asks whether the upload holds gene symbols at all;
    that does not depend on how the sources were combined. Under intersection
    the combined sets can shrink to a handful of genes, so measuring against
    them reports a correct gene-symbol column as the wrong identifier type.

    Returns
    -------
    set[str]
        Uppercase gene symbols.
    """
    universe = getattr(reference_sets, 'gene_universe', None)
    if universe is not None:
        return set(universe)
    found = set()
    for genes in reference_sets.values():
        found.update(str(g).strip().upper() for g in genes)
    return found


def parse_source_combination(mode, min_n=None, n_sources=1):
    """Validate the source-combination form fields into one canonical spec (#123).

    Parameters
    ----------
    mode : str or None
        ``"union"``, ``"intersection"`` or ``"at_least"``; blank means union.
        Case and surrounding whitespace are ignored.
    min_n : str or int or None
        The N of "at least N sources". Read only when ``mode`` is ``at_least``.
    n_sources : int
        How many resources the user selected; N is validated against it.

    Returns
    -------
    str
        ``"union"``, ``"intersection"`` or ``"at_least_<N>"``. The endpoints
        collapse — N=1 is union and N=n_sources is intersection — and any mode
        over a single source is union, because with one source they are the same
        set. One spelling per gene set keeps stored runs comparable.

    Raises
    ------
    ValueError
        For an unknown mode, or an N that is not a whole number between 1 and
        ``n_sources``. The message is meant for the user.
    """
    mode = (mode or DEFAULT_SOURCE_COMBINATION).strip().lower()
    if mode not in VALID_SOURCE_COMBINATION_MODES:
        raise ValueError(
            "Source combination must be 'union', 'intersection' or 'at_least'."
        )
    n_sources = max(int(n_sources), 1)
    if mode == SOURCE_COMBINATION_AT_LEAST:
        try:
            required = int(str(min_n).strip())
        except (TypeError, ValueError):
            raise ValueError(
                "'At least N sources' needs N as a whole number."
            ) from None
        if not 1 <= required <= n_sources:
            raise ValueError(
                f"'At least N sources' needs N between 1 and the number of selected "
                f"sources — at most {n_sources} here. Got: {required}."
            )
    elif mode == SOURCE_COMBINATION_INTERSECTION:
        required = n_sources
    else:
        required = 1
    return _canonical_source_combination(required, n_sources)


def _canonical_source_combination(required, n_sources):
    """Spell "a gene needs ``required`` of ``n_sources`` sources" canonically."""
    if required <= 1 or n_sources <= 1:
        return SOURCE_COMBINATION_UNION
    if required >= n_sources:
        return SOURCE_COMBINATION_INTERSECTION
    return f"{_AT_LEAST_PREFIX}{required}"


def _at_least_n(spec):
    """Return N for an ``at_least_<N>`` spec with N >= 2, else None."""
    if not isinstance(spec, str) or not spec.startswith(_AT_LEAST_PREFIX):
        return None
    suffix = spec[len(_AT_LEAST_PREFIX):]
    if not suffix.isdigit() or int(suffix) < 2:
        return None
    return int(suffix)


def normalise_source_combination(value):
    """Read a stored source-combination spec, defaulting to union (#123).

    NULL — every run recorded before #123 — and anything unrecognised read as
    ``"union"``, which is what those runs did.
    """
    if value in (SOURCE_COMBINATION_UNION, SOURCE_COMBINATION_INTERSECTION):
        return value
    if _at_least_n(value) is not None:
        return value
    return DEFAULT_SOURCE_COMBINATION


def source_combination_label(spec):
    """Human-readable label for results pages, provenance lines and reports."""
    spec = normalise_source_combination(spec)
    if spec == SOURCE_COMBINATION_INTERSECTION:
        return "Intersection (genes in every selected source)"
    n = _at_least_n(spec)
    if n is not None:
        return f"At least {n} selected sources per gene"
    return "Union (genes in any selected source)"


def combine_source_sets(per_source_sets, spec=DEFAULT_SOURCE_COMBINATION):
    """Combine per-resource KE -> gene sets under a source combination (#123).

    Parameters
    ----------
    per_source_sets : list[dict]
        One KE_ID -> gene set mapping per resource that loaded. Not mutated.
    spec : str
        Canonical spec from :func:`parse_source_combination`.

    Returns
    -------
    tuple[dict, set]
        ``(combined, no_shared)``. Under union ``combined`` is built exactly as
        before #123 and ``no_shared`` is empty. Otherwise a gene is kept when at
        least the required number of sources list it for that Key Event — every
        source for intersection — counted over the sources that actually loaded,
        so a skipped resource cannot veto every gene. A Key Event that has genes
        in some source but keeps none is left out of ``combined`` and named in
        ``no_shared``, so the enrichment can report why it was not tested. A Key
        Event with no genes in any source is kept as an empty set, as union
        keeps it: that is a mapping that resolved to nothing (#81), not a
        disagreement between sources.
    """
    spec = normalise_source_combination(spec)
    n_sources = len(per_source_sets)
    if spec == SOURCE_COMBINATION_INTERSECTION:
        required = n_sources
    else:
        required = _at_least_n(spec) or 1

    combined = {}
    if required <= 1 or n_sources <= 1:
        # The pre-#123 merge, unchanged, so a union run is byte-identical.
        for sets in per_source_sets:
            for ke_id, genes in sets.items():
                combined.setdefault(ke_id, set()).update(genes)
        return combined, set()

    counts = {}
    for sets in per_source_sets:
        for ke_id, genes in sets.items():
            per_ke = counts.setdefault(ke_id, {})
            # A source listing a gene twice (case variants) still counts once.
            for gene in {str(g).strip().upper() for g in genes}:
                per_ke[gene] = per_ke.get(gene, 0) + 1

    no_shared = set()
    for ke_id, per_ke in counts.items():
        kept = {gene for gene, count in per_ke.items() if count >= required}
        if kept or not per_ke:
            combined[ke_id] = kept
        else:
            no_shared.add(ke_id)
    return combined, no_shared


def confidence_rank(value):
    """Map a confidence label to its ordinal rank.

    Parameters
    ----------
    value : str or None
        Confidence label such as ``"High"``, ``"medium"`` or ``None``.

    Returns
    -------
    int or None
        ``3`` for High, ``2`` for Medium, ``1`` for Low; ``None`` when the value
        is missing, blank, or not a recognised level. ``None`` means "unknown"
        and callers treat it as a filter no-op (issue #60 graceful degradation).
    """
    if value is None:
        return None
    text = str(value).strip().lower()
    if not text or text in ('nan', 'none'):
        return None
    return _CONFIDENCE_RANKS.get(text)


def record_confidence_rank(record):
    """Return the confidence rank of a mapping record, or None when absent.

    Looks the value up under every known confidence key so both Builder API
    records (``confidence_level``) and CSV-derived records are handled.
    """
    for key in CONFIDENCE_KEYS:
        if key in record:
            rank = confidence_rank(record.get(key))
            if rank is not None:
                return rank
    return None


def filter_records_by_confidence(records, min_confidence=DEFAULT_MIN_CONFIDENCE):
    """Drop KE-mapping records below a minimum confidence threshold (issue #60).

    Records whose confidence is missing, blank, or unrecognised are **kept** —
    the filter is a no-op for them. That keeps the local ``data/KE-WP.csv``
    fallback (no confidence column) and any resource lacking the field fully
    usable, so an analysis is never silently emptied by the filter alone.

    Parameters
    ----------
    records : list[dict]
        Raw mapping records, e.g. from ``fetch_ke_wp_records()`` or
        ``load_ke_wp_records_csv()``.
    min_confidence : str
        One of ``"all"`` (default, no filtering), ``"medium"`` (Medium+High) or
        ``"high"``. Unknown values are treated as ``"all"``.

    Returns
    -------
    list[dict]
        The qualifying records as a new list; the input is not mutated.
    """
    threshold = confidence_rank(min_confidence)
    if min_confidence not in VALID_MIN_CONFIDENCE or threshold is None:
        # "all" (or an unrecognised threshold) -> no filtering.
        return list(records)

    kept = [r for r in records if (record_confidence_rank(r) or threshold) >= threshold]
    logger.info(
        "filter_records_by_confidence: kept %d of %d records at min_confidence=%s",
        len(kept), len(records), min_confidence,
    )
    return kept


def filter_ke_wp_df_by_confidence(ke_wp_df, min_confidence=DEFAULT_MIN_CONFIDENCE):
    """Drop below-threshold rows from a KE-WP DataFrame (issue #60).

    A no-op when ``min_confidence`` is ``"all"``/unknown, or when the frame
    carries no confidence column at all (the ``data/KE-WP.csv`` fallback shape).

    Parameters
    ----------
    ke_wp_df : pandas.DataFrame
        KE-WP mappings, optionally with one of :data:`CONFIDENCE_KEYS` as a column.
    min_confidence : str
        ``"all"``, ``"medium"`` or ``"high"``.

    Returns
    -------
    pandas.DataFrame
        The qualifying rows (the same frame when the filter is a no-op).
    """
    threshold = confidence_rank(min_confidence)
    if min_confidence not in VALID_MIN_CONFIDENCE or threshold is None:
        return ke_wp_df

    column = next((c for c in CONFIDENCE_KEYS if c in ke_wp_df.columns), None)
    if column is None:
        logger.info(
            "filter_ke_wp_df_by_confidence: no confidence column present; "
            "min_confidence=%s applied as a no-op", min_confidence
        )
        return ke_wp_df

    ranks = ke_wp_df[column].map(confidence_rank)
    # Unknown/missing confidence (NaN rank) is kept — graceful degradation.
    keep = ranks.isna() | (ranks >= threshold)
    # .copy() so downstream in-place column normalisation never writes to a slice.
    filtered = ke_wp_df[keep].copy()
    logger.info(
        "filter_ke_wp_df_by_confidence: kept %d of %d KE-WP rows at min_confidence=%s",
        len(filtered), len(ke_wp_df), min_confidence,
    )
    return filtered


def load_reference_sets(ke_wp_path=os.path.join(DATA_DIR, 'KE-WP.csv'),
                        wp_gene_path=os.path.join(DATA_DIR, 'edges_wpid_to_gene.csv'),
                        node_path=os.path.join(DATA_DIR, 'node_attributes.csv'),
                        ke_wp_df=None,
                        min_confidence=DEFAULT_MIN_CONFIDENCE,
                        wp_gene_map=None,
                        unresolved_out=None,
                        unresolved_ke_out=None):
    """Build KE-to-gene reference sets from local CSV files or a pre-built DataFrame.

    Parameters
    ----------
    ke_wp_path : str
        Path to the KE-WP mapping CSV file. Ignored when ``ke_wp_df`` is provided.
    wp_gene_path : str
        Path to the WikiPathways-to-gene edge CSV file.
    node_path : str
        Path to the node attributes CSV file (gene ID to symbol mapping).
    ke_wp_df : pandas.DataFrame or None
        Optional pre-built DataFrame with columns ``KE_ID`` and ``WP_ID``.
        When provided, the CSV at ``ke_wp_path`` is not read. KE IDs must
        already be normalised to the ``"KE:55"`` colon format.
    min_confidence : str
        Issue #60 minimum mapping confidence — ``"all"`` (default, no
        filtering), ``"medium"`` or ``"high"``. Rows whose confidence is below
        the threshold are dropped before the WP-to-gene merge, so a KE keeps
        only the genes of its qualifying pathways and a KE whose mappings are
        all below threshold is absent from the result. The filter is a no-op
        when the mappings carry no confidence column (CSV fallback).
    wp_gene_map : dict[str, set[str]] or None
        Optional pathway-to-gene membership, normally the Builder's live copy
        from :func:`services.api_service.fetch_wp_pathway_gene_map`. Takes
        precedence over ``wp_gene_path`` per pathway; the CSV still covers any
        pathway the map lacks. Issue #79: the bundled CSV is a snapshot, so
        pathways curated after it was taken resolved to no genes at all.
    unresolved_out : list or None
        Optional sink. Pathway IDs that neither source can resolve are appended
        here so the caller can surface them. Without it such a pathway vanishes
        in the inner join, which is exactly the silence #79 reported: the Key
        Event is then reported as having no gene set mapped, when in truth it
        is mapped and the mapping could not be resolved.
    unresolved_ke_out : dict or None
        Optional sink for the same information attributed to Key Events:
        KE_ID -> sorted list of that KE's unresolvable pathway IDs, taken from
        the mappings *before* the inner merge drops them. Issue #81: a KE whose
        only pathway is unresolvable is absent from the returned sets entirely,
        so this is the only place the association still exists. The same map is
        always attached to the return value (see :class:`ReferenceSets`); the
        sink exists for callers that want it without depending on the type.

    Returns
    -------
    ReferenceSets
        Mapping of KE ID -> set of uppercase gene symbols, carrying the
        KE -> unresolvable-pathway-ID map on ``.unresolved_ke_pathways``.
    """
    # Use the provided DataFrame or read from the CSV file
    if ke_wp_df is None:
        ke_wp_df = pd.read_csv(ke_wp_path)

    # Issue #60: apply the minimum-confidence threshold before any merging.
    ke_wp_df = filter_ke_wp_df_by_confidence(ke_wp_df, min_confidence)

    wp_gene_df = pd.read_csv(wp_gene_path)
    node_df = pd.read_csv(node_path)

    logger.debug("Sample gene_id values from wp_gene_df: %s", wp_gene_df['gene_id'].dropna().unique()[:10])
    logger.debug("Sample GeneID values from node_df: %s", node_df['GeneID'].dropna().unique()[:10])

    # Normalize formatting
    ke_wp_df['WP_ID'] = ke_wp_df['WP_ID'].astype(str).str.upper().str.strip()
    wp_gene_df['WPID'] = wp_gene_df['WPID'].astype(str).str.upper().str.strip()
    wp_gene_df['gene_id'] = wp_gene_df['gene_id'].astype(str).str.strip()
    node_df['GeneID'] = (
        node_df['GeneID']
        .dropna()                       # remove NaNs
        .astype(float)                 # ensure float first
        .astype(int)                   # then cast to int safely
        .astype(str)                   # convert to string
        .str.strip()                   # strip whitespace
    )


    node_df['GeneName'] = node_df['GeneName'].astype(str).str.strip()

    logger.debug("KE-WP mappings: %s", ke_wp_df.shape)
    logger.debug("Gene-pathway mappings: %s", wp_gene_df.shape)
    logger.debug("Node annotation mappings: %s", node_df.shape)

    # Merge gene → annotation
    wp_gene_annotated = wp_gene_df.merge(
        node_df[['GeneID', 'GeneName']],
        left_on='gene_id',
        right_on='GeneID',
        how='left'
    )

    logger.debug("After annotation merge: %s", wp_gene_annotated.shape)

    # Issue #79: prefer the Builder's live pathway membership, per pathway, and
    # fall back to the bundled CSV for anything it does not carry. Rows are
    # appended rather than replacing the CSV wholesale so a partial Builder
    # response cannot lose pathways the snapshot does cover.
    if wp_gene_map:
        wanted = set(ke_wp_df['WP_ID'].unique())
        live_rows = [
            {'WPID': wp_id, 'GeneName': gene}
            for wp_id in wanted
            for gene in wp_gene_map.get(wp_id, ())
        ]
        if live_rows:
            live_df = pd.DataFrame(live_rows)
            covered = set(live_df['WPID'].unique())
            # Drop the snapshot's rows for pathways the Builder answered for, so
            # a pathway whose membership changed upstream does not end up as the
            # union of both versions.
            csv_only = wp_gene_annotated[~wp_gene_annotated['WPID'].isin(covered)]
            wp_gene_annotated = pd.concat(
                [csv_only[['WPID', 'GeneName']], live_df], ignore_index=True
            )
            logger.info(
                "WP gene membership: %d of %d mapped pathways resolved from the "
                "Builder, remainder from the bundled CSV",
                len(covered & wanted), len(wanted),
            )

    # Merge WP → KE
    merged = ke_wp_df.merge(
        wp_gene_annotated,
        left_on='WP_ID',
        right_on='WPID',
        how='inner'
    )

    pre_dropna_count = len(merged)
    logger.debug("Final merged KE-gene mapping before dropna: %s", merged.shape)

    # Drop entries without gene symbol
    merged = merged.dropna(subset=['GeneName'])

    logger.debug("Final merged KE-gene mapping after dropna: %s", merged.shape)

    # Issue #79: name the pathways that resolved to nothing, rather than
    # reporting a count of "lost mappings" that compared rows against unique
    # KEs and so could not be acted on. A pathway here is curated and mapped;
    # the analyser simply cannot say which genes it contains.
    resolvable = set(merged['WP_ID'].unique()) if len(merged) else set()
    unresolved = sorted(set(ke_wp_df['WP_ID'].unique()) - resolvable)
    unresolved_by_ke = {}
    if unresolved:
        logger.warning(
            "%d mapped pathway(s) could not be resolved to genes and are absent "
            "from the reference sets: %s",
            len(unresolved), ", ".join(unresolved),
        )
        if unresolved_out is not None:
            unresolved_out.extend(unresolved)

        # Issue #81: attribute the unresolvable pathways to their Key Events
        # here, on the pre-merge mappings. After the inner merge below the
        # association is gone — a KE whose only pathway is unresolvable does
        # not appear in the merged frame at all, so downstream it is
        # indistinguishable from a KE nobody has ever mapped. That is exactly
        # the misreport #81 is about, and this map is what fixes it.
        unresolved_set = set(unresolved)
        collected = {}
        for ke_id, wp_id in zip(ke_wp_df['KE_ID'], ke_wp_df['WP_ID']):
            if wp_id in unresolved_set:
                collected.setdefault(str(ke_id), set()).add(wp_id)
        unresolved_by_ke = {ke: sorted(wps) for ke, wps in collected.items()}
        logger.info(
            "%d Key Event(s) carry at least one unresolvable pathway mapping",
            len(unresolved_by_ke),
        )
        if unresolved_ke_out is not None:
            unresolved_ke_out.update(unresolved_by_ke)

    # Group into dict: KE_ID → set of gene symbols
    reference_sets = (
        merged.groupby('KE_ID')['GeneName']
        .apply(lambda g: set(g.str.upper()))
        .to_dict()
    )

    return ReferenceSets(reference_sets, unresolved_ke_pathways=unresolved_by_ke)

