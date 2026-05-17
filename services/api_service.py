"""
Builder API client service for fetching KE-WP mappings.

Provides paginated fetching with retry logic, KE ID normalisation,
and integration with the local reference set loading pipeline.
"""
import logging
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from helpers import load_reference_sets

logger = logging.getLogger(__name__)


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


def fetch_ke_wp_records(config):
    """Fetch raw KE-WP mapping records from the Builder API for the WP picker.

    Returns the full record list from ``fetch_all_ke_wp_mappings`` unchanged so
    that each record retains its ``ke_id``, ``pathway_id``, ``pathway_title``,
    and ``confidence_level`` fields.  These records are consumed directly by
    ``build_wp_picker_data()`` in ``services/pathway_picker_service.py``.

    Note: ``fetch_reference_sets_from_api()`` is intentionally left untouched —
    the reference-set loading pipeline is a separate concern from the picker-data
    pipeline (Option A from the research notes — additive, no refactor).

    Parameters
    ----------
    config : Config
        Application config object providing ``BUILDER_API_URL`` and
        ``BUILDER_API_TIMEOUT`` attributes.

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
    return records


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


def fetch_reference_sets_from_api(config):
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

    Returns
    -------
    dict[str, set[str]]
        Mapping of KE ID -> set of gene symbols, identical in shape to
        the CSV-based output of :func:`helpers.load_reference_sets`.

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

    # Normalise KE IDs: "KE 55" -> "KE:55"
    ke_wp_data = [
        {
            "KE_ID": record["ke_id"].replace(" ", ":"),
            "WP_ID": record["pathway_id"],
        }
        for record in records
    ]
    ke_wp_df = pd.DataFrame(ke_wp_data, columns=["KE_ID", "WP_ID"])

    logger.info(
        "Built KE-WP DataFrame with %d rows from %d API records",
        len(ke_wp_df),
        len(records),
    )

    reference_sets = load_reference_sets(ke_wp_df=ke_wp_df)

    logger.info(
        "Reference sets loaded: %d KE gene sets produced from API data",
        len(reference_sets),
    )
    return reference_sets
