"""
Builder API client service for fetching KE-WP mappings.

Provides paginated fetching with retry logic, KE ID normalisation,
and integration with the local reference set loading pipeline.
"""
import logging
import re

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from helpers import (
    DEFAULT_MIN_CONFIDENCE,
    VALID_MIN_CONFIDENCE,
    confidence_rank,
    filter_records_by_confidence,
    load_reference_sets,
    record_confidence_rank,
)

logger = logging.getLogger(__name__)

# Builder GMT export paths that already resolve pathways/terms to genes.
# WikiPathways keeps its own KE-WP pipeline (mappings API + confidence filter),
# but its GMT is still used as the *pathway-to-gene* source — see
# WP_GMT_PATH and fetch_wp_pathway_gene_map.
GMT_RESOURCE_PATHS = {
    "GO_BP": "exports/gmt/ke-go",
    "Reactome": "exports/gmt/ke-reactome",
}

# The KE-WP GMT. Not in GMT_RESOURCE_PATHS because WikiPathways reference sets
# are not built from it directly: the GMT carries no confidence level, so
# building KE sets from it would silently ignore min_confidence (#67). It is
# read only to learn which genes each *pathway* contains.
WP_GMT_PATH = "exports/gmt/ke-wp"

# Matches the leading "KE<number>" token of a GMT descriptor column, e.g.
# "KE177_Increase_Mitochondrial_dysfunction_WP5241" -> "177".
_KE_ID_RE = re.compile(r"^KE\s?(\d+)_")

# Matches the trailing pathway token of the same descriptor -> "WP5241".
_WP_ID_RE = re.compile(r"_(WP\d+)$", re.IGNORECASE)


def _make_api_session(retries=3, backoff_factor=1.0):
    """Create a requests Session with retry and backoff configured.

    Parameters
    ----------
    retries : int
        Number of retry attempts on transient failures (default 3).
    backoff_factor : float
        Exponential backoff multiplier; delays are 0s, 2s, 4s (default 1.0).

    Returns
    -------
    requests.Session
        Session with the retry adapter mounted for both http and https.
    """
    retry_strategy = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session = requests.Session()
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def fetch_all_ke_wp_mappings(session, base_url, timeout):
    """Fetch all KE-WP mapping records from the Builder API using pagination.

    Requests pages of up to 200 records from ``{base_url}/api/v1/mappings``
    until the ``pagination.next`` field is ``null``.

    Parameters
    ----------
    session : requests.Session
        An authenticated (or plain) session, typically from
        :func:`_make_api_session`.
    base_url : str
        Base URL of the Builder API (no trailing slash).
    timeout : int
        Request timeout in seconds.

    Returns
    -------
    list[dict]
        All records across all pages; each record contains at minimum
        ``ke_id`` and ``pathway_id`` fields.

    Raises
    ------
    requests.HTTPError
        On non-2xx responses after all retries are exhausted.
    """
    all_records = []
    url = f"{base_url.rstrip('/')}/api/v1/mappings"
    params = {"per_page": 200}
    page = 1

    while url is not None:
        logger.info("Fetching page %d from %s", page, url)
        response = session.get(url, params=params, timeout=timeout)
        response.raise_for_status()
        payload = response.json()

        records = payload.get("data", [])
        all_records.extend(records)
        logger.debug("Page %d returned %d records", page, len(records))

        # Advance to next page or stop
        pagination = payload.get("pagination", {})
        url = pagination.get("next")  # None when last page reached
        params = {}  # next URL already contains query params
        page += 1

    logger.info("Total records fetched from API: %d", len(all_records))
    return all_records


def fetch_ke_wp_records(config, min_confidence=DEFAULT_MIN_CONFIDENCE):
    """Fetch raw KE-WP mapping records from the Builder API for the WP picker.

    Returns the record list from ``fetch_all_ke_wp_mappings`` with each record's
    ``ke_id``, ``pathway_id``, ``pathway_title`` and ``confidence_level`` fields
    intact.  These records are consumed directly by ``build_wp_picker_data()``
    in ``services/pathway_picker_service.py``.  Issue #60 adds an optional
    minimum-confidence filter so the picker matches the mappings that were
    actually enriched.

    Parameters
    ----------
    config : Config
        Application config object providing ``BUILDER_API_URL`` and
        ``BUILDER_API_TIMEOUT`` attributes.
    min_confidence : str
        Issue #60 minimum mapping confidence (``"all"``/``"medium"``/``"high"``).
        Records below the threshold are dropped client-side; records without a
        confidence value are kept.

    Returns
    -------
    list[dict]
        Raw mapping records; each record contains at minimum ``ke_id``,
        ``pathway_id``, ``pathway_title``, and ``confidence_level``.

    Raises
    ------
    ValueError
        If ``config.BUILDER_API_URL`` is empty or falsy (API not configured).
    requests.HTTPError
        On non-2xx API responses after all retries are exhausted.
    """
    if not config.BUILDER_API_URL:
        raise ValueError("BUILDER_API_URL not configured; API integration is disabled")

    session = _make_api_session()
    records = fetch_all_ke_wp_mappings(
        session=session,
        base_url=config.BUILDER_API_URL,
        timeout=config.BUILDER_API_TIMEOUT,
    )

    logger.info("fetch_ke_wp_records: fetched %d raw KE-WP records from API", len(records))
    # Issue #60: apply the minimum-confidence threshold client-side. The Builder
    # confidence_level query param is exact-match only, which cannot express
    # "Medium and above", and the full record list is needed anyway.
    return filter_records_by_confidence(records, min_confidence)


def load_ke_wp_records_csv(csv_path="data/KE-WP.csv"):
    """Load KE-WP mapping records from the local CSV fallback file.

    Reads the two-column CSV at ``csv_path`` and returns the rows as a list of
    dicts each with keys ``KE_ID`` and ``WP_ID``.  This is the fallback record
    source used when the Builder API is unreachable, providing the same interface
    as ``fetch_ke_wp_records()`` but with less metadata (no ``pathway_title`` or
    ``confidence_level`` columns).

    Parameters
    ----------
    csv_path : str
        Path to the CSV file (default: ``"data/KE-WP.csv"``).

    Returns
    -------
    list[dict]
        List of dicts, each with keys ``KE_ID`` (e.g. ``"KE:245"``) and
        ``WP_ID`` (e.g. ``"WP2876"``).
    """
    df = pd.read_csv(csv_path)
    records = df.to_dict("records")
    logger.info("load_ke_wp_records_csv: loaded %d KE-WP records from %s", len(records), csv_path)
    return records


def fetch_reference_sets_from_api(config, min_confidence=DEFAULT_MIN_CONFIDENCE):
    """Orchestrate API fetching and build reference sets.

    Fetches all KE-WP mapping records from the Builder API, normalises
    KE IDs from ``"KE 55"`` (space) to ``"KE:55"`` (colon), constructs
    a DataFrame, and delegates to :func:`helpers.load_reference_sets` for
    the WP-to-gene merge.

    Parameters
    ----------
    config : Config
        Application config object providing ``BUILDER_API_URL``,
        ``BUILDER_API_TIMEOUT`` attributes.
    min_confidence : str
        Issue #60 minimum mapping confidence (``"all"``/``"medium"``/``"high"``).
        Mappings below the threshold are dropped before the gene sets are built,
        so a KE keeps only the genes of its qualifying pathways and a KE whose
        mappings are all below threshold disappears entirely.

    Returns
    -------
    tuple[dict[str, set[str]], list[str]]
        The KE ID -> gene symbol mapping (identical in shape to the CSV-based
        output of :func:`helpers.load_reference_sets`), and the sorted list of
        mapped pathway IDs that could not be resolved to genes by either the
        Builder or the bundled snapshot (issue #79). The second element is
        normally empty; when it is not, those Key Events are missing gene sets
        they are curated to have, and the caller must surface that.

    Raises
    ------
    ValueError
        If ``config.BUILDER_API_URL`` is empty or falsy, indicating that
        API integration is not configured (not a network failure).
    requests.HTTPError
        On non-2xx API responses after all retries are exhausted.
    """
    if not config.BUILDER_API_URL:
        raise ValueError("BUILDER_API_URL not configured; API integration is disabled")

    session = _make_api_session()
    records = fetch_all_ke_wp_mappings(
        session=session,
        base_url=config.BUILDER_API_URL,
        timeout=config.BUILDER_API_TIMEOUT,
    )

    # Issue #60: drop below-threshold mappings before the gene sets are built.
    qualifying = filter_records_by_confidence(records, min_confidence)

    # Normalise KE IDs: "KE 55" -> "KE:55"
    ke_wp_data = [
        {
            "KE_ID": record["ke_id"].replace(" ", ":"),
            "WP_ID": record["pathway_id"],
        }
        for record in qualifying
    ]
    ke_wp_df = pd.DataFrame(ke_wp_data, columns=["KE_ID", "WP_ID"])

    logger.info(
        "Built KE-WP DataFrame with %d rows from %d API records (min_confidence=%s)",
        len(ke_wp_df),
        len(records),
        min_confidence,
    )

    # Issue #79: resolve pathway membership from the Builder, not only from the
    # bundled snapshot. Without this a pathway curated after the snapshot was
    # taken contributes no genes, and the KE is reported as unmapped.
    wp_gene_map = fetch_wp_pathway_gene_map(config)
    unresolved = []
    reference_sets = load_reference_sets(
        ke_wp_df=ke_wp_df,
        wp_gene_map=wp_gene_map,
        unresolved_out=unresolved,
    )

    logger.info(
        "Reference sets loaded: %d KE gene sets produced from API data",
        len(reference_sets),
    )
    return reference_sets, sorted(set(unresolved))


def parse_gmt_pathway_gene_map(gmt_text):
    """Parse a KE-WP GMT export into a *pathway*-to-gene map.

    :func:`parse_gmt_reference_sets` collapses the same file to KE -> genes.
    This keeps the pathway dimension instead, because the WikiPathways pipeline
    needs to apply the confidence threshold to KE->pathway mappings first and
    only then look up each qualifying pathway's genes.

    Genes are unioned across every row naming the same pathway (a pathway
    mapped by several KEs appears several times, with identical membership).

    Parameters
    ----------
    gmt_text : str
        Raw GMT file contents.

    Returns
    -------
    dict[str, set[str]]
        Mapping of upper-case pathway ID (``"WP5477"``) -> gene symbols.
    """
    pathway_genes = {}
    for line in gmt_text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue  # descriptor + title but no genes
        match = _WP_ID_RE.search(fields[0])
        if not match:
            continue
        wp_id = match.group(1).upper()
        genes = {g.strip().upper() for g in fields[2:] if g.strip()}
        if not genes:
            continue
        pathway_genes.setdefault(wp_id, set()).update(genes)
    return pathway_genes


def fetch_wp_pathway_gene_map(config):
    """Fetch pathway-to-gene membership for WikiPathways from the Builder.

    The bundled ``data/edges_wpid_to_gene.csv`` is a snapshot, so any pathway
    curated after it was taken resolves to no genes at all — silently, because
    the merge is an inner join (issue #79). The Builder resolves pathway
    membership itself for its GMT export, so its copy is the fresher source.

    Parameters
    ----------
    config : Config
        Application config providing ``BUILDER_API_URL`` and
        ``BUILDER_API_TIMEOUT``.

    Returns
    -------
    dict[str, set[str]]
        Pathway ID -> gene symbols. Empty dict if the Builder is not
        configured or the export cannot be fetched — callers fall back to the
        bundled CSV rather than failing, since a stale gene set beats none.
    """
    if not config.BUILDER_API_URL:
        logger.info("BUILDER_API_URL not configured; skipping WP pathway-gene fetch")
        return {}

    url = f"{config.BUILDER_API_URL.rstrip('/')}/{WP_GMT_PATH}"
    try:
        session = _make_api_session()
        response = session.get(url, timeout=config.BUILDER_API_TIMEOUT)
        response.raise_for_status()
        pathway_genes = parse_gmt_pathway_gene_map(response.text)
        logger.info(
            "Loaded gene membership for %d WikiPathways pathways from %s",
            len(pathway_genes), url,
        )
        return pathway_genes
    except Exception as exc:
        logger.warning(
            "Could not fetch WikiPathways gene membership from the Builder (%s); "
            "falling back to the bundled CSV snapshot", exc
        )
        return {}

def parse_gmt_reference_sets(gmt_text):
    """Parse a Builder GMT export into KE-to-gene reference sets.

    The Builder GMT format is one gene set per line, tab-separated::

        KE<id>_<KE_name>_<pathway_id>\\t<pathway_title>\\t<gene>\\t<gene>...

    A single KE typically spans several rows (one per mapped pathway/term);
    genes are unioned across all rows that share the same KE.  KE IDs are
    normalised from the ``KE<id>`` token to ``"KE:<id>"`` so they match the
    convention used elsewhere (see :func:`fetch_reference_sets_from_api` and
    the AOP ``ke_list`` keys consumed by enrichment).

    Parameters
    ----------
    gmt_text : str
        Raw GMT file contents (``text/plain``).

    Returns
    -------
    dict[str, set[str]]
        Mapping of KE ID (``"KE:177"``) -> set of uppercase gene symbols.
    """
    reference_sets = {}
    for line in gmt_text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue  # descriptor + title but no genes -> nothing to add
        match = _KE_ID_RE.match(fields[0])
        if not match:
            continue
        ke_id = f"KE:{match.group(1)}"
        genes = {g.strip().upper() for g in fields[2:] if g.strip()}
        if not genes:
            continue
        reference_sets.setdefault(ke_id, set()).update(genes)
    return reference_sets


def gmt_min_confidence_param(min_confidence):
    """Translate our confidence vocabulary into the Builder's query param (#71).

    The two vocabularies do not match. Ours is
    :data:`helpers.VALID_MIN_CONFIDENCE` — ``("all", "medium", "high")``, where
    ``all`` is a sentinel meaning "do not filter". The Builder's whitelist is
    ``{high, medium, low}``, and ``all`` is **not in it**: forwarding it
    verbatim returns HTTP 400 and loses the resource for the whole run. So
    ``all`` becomes an omitted parameter, which under the post-builder#206
    semantics is exactly equivalent (``low`` is also equivalent, but omitting
    says what we mean).

    Args:
        min_confidence: a value from ``VALID_MIN_CONFIDENCE``. Anything
            unrecognised is treated as ``all`` — an unfiltered gene set is the
            safe failure here, since the alternative is a 400 that skips the
            resource entirely.

    Returns:
        dict: query params for ``requests``; empty when no filtering applies.
        The Builder matches case-insensitively (verified live 2026-07-22
        against `High`/`high`/`HIGH`), so our lowercase values go through
        unchanged.
    """
    if min_confidence in ("medium", "high"):
        return {"min_confidence": min_confidence}
    return {}


def fetch_gmt_reference_sets(config, resource, min_confidence=DEFAULT_MIN_CONFIDENCE):
    """Fetch KE-to-gene reference sets for a resource from the Builder GMT export.

    Used for resources whose genes are resolved Builder-side and exposed via the
    GMT export endpoints (GO Biological Process and Reactome).  WikiPathways is
    not handled here — it retains its dedicated CSV-backed pipeline.

    Parameters
    ----------
    config : Config
        Application config providing ``BUILDER_API_URL`` and
        ``BUILDER_API_TIMEOUT``.
    resource : str
        Resource key, one of :data:`GMT_RESOURCE_PATHS` (``"GO_BP"`` or
        ``"Reactome"``).
    min_confidence : str
        Issue #71 threshold, from :data:`helpers.VALID_MIN_CONFIDENCE`.
        Forwarded to the Builder as ``?min_confidence=`` via
        :func:`gmt_min_confidence_param`.

        This became possible only with **molAOP-builder#206** (2026-07-22).
        Before that the parameter selected a *single tier*, so passing
        ``medium`` returned medium-confidence mappings and dropped every
        high-confidence one — the opposite of what the name promises, and the
        reason the Analyser disclosed the limitation instead of forwarding the
        threshold. It now means at-or-above, and mappings with no recorded
        confidence are always included, so a threshold can never silently
        empty an export.

    Returns
    -------
    dict[str, set[str]]
        Mapping of KE ID -> set of gene symbols, identical in shape to the
        WikiPathways reference sets.

    Raises
    ------
    ValueError
        If ``config.BUILDER_API_URL`` is unset, or ``resource`` is unknown.
    requests.HTTPError
        On non-2xx responses after all retries are exhausted.
    """
    if not config.BUILDER_API_URL:
        raise ValueError("BUILDER_API_URL not configured; API integration is disabled")
    if resource not in GMT_RESOURCE_PATHS:
        raise ValueError(f"Unknown GMT resource: {resource!r}")

    url = f"{config.BUILDER_API_URL.rstrip('/')}/{GMT_RESOURCE_PATHS[resource]}"
    params = gmt_min_confidence_param(min_confidence)
    session = _make_api_session()
    response = session.get(url, params=params, timeout=config.BUILDER_API_TIMEOUT)
    response.raise_for_status()

    reference_sets = parse_gmt_reference_sets(response.text)
    logger.info(
        "Loaded %d KE gene sets for resource %s from GMT export %s (min_confidence=%s%s)",
        len(reference_sets),
        resource,
        url,
        min_confidence,
        "" if params else ", sent unfiltered",
    )
    return reference_sets
