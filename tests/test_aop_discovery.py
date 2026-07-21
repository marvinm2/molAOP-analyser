"""
Tests for the AOP discovery service and /api/aops endpoint.

Covers:
- _uri_to_aop_id URI conversion
- _hardcoded_aop_list structure and content
- get_aop_list three-tier fallback (both live fetch failure and cache hit)
- fetch_mapped_ke_ids_from_builder with empty BUILDER_API_URL
- /api/aops endpoint (no query, with query, error propagation)
- load_aop_data() default source routing for dynamic AOPs
"""

import os
import sys
import tempfile
from unittest.mock import MagicMock, patch

import pytest

# Ensure project root is on the path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_config(builder_url="", cache_dir=None):
    """Return a minimal mock config object for tests."""
    cfg = MagicMock()
    cfg.BUILDER_API_URL = builder_url
    cfg.BUILDER_API_TIMEOUT = 5
    cfg.CACHE_DIR = cache_dir or tempfile.mkdtemp()
    cfg.SPARQL_ENDPOINT = "https://aopwiki.rdf.bigcat-bioinformatics.org/sparql/"
    cfg.SPARQL_TIMEOUT = 10
    return cfg


# ---------------------------------------------------------------------------
# Unit: URI conversion
# ---------------------------------------------------------------------------

class TestUriToAopId:
    """Tests for _uri_to_aop_id."""

    def test_standard_uri(self):
        """Convert a typical identifiers.org AOP URI."""
        from services.aop_discovery_service import _uri_to_aop_id

        result = _uri_to_aop_id("https://identifiers.org/aop/619")
        assert result == "AOP:619"

    def test_low_id(self):
        """Convert a small numeric AOP ID."""
        from services.aop_discovery_service import _uri_to_aop_id

        assert _uri_to_aop_id("https://identifiers.org/aop/1") == "AOP:1"

    def test_invalid_uri_raises(self):
        """Raise ValueError for an unrecognised URI format."""
        from services.aop_discovery_service import _uri_to_aop_id

        with pytest.raises(ValueError):
            _uri_to_aop_id("https://example.com/not-an-aop/619")


# ---------------------------------------------------------------------------
# Unit: hardcoded fallback list
# ---------------------------------------------------------------------------

class TestHardcodedAopList:
    """Tests for _hardcoded_aop_list."""

    def test_returns_non_empty_list(self):
        """The hardcoded list must contain at least one entry."""
        from services.aop_discovery_service import _hardcoded_aop_list

        result = _hardcoded_aop_list()
        assert isinstance(result, list)
        assert len(result) > 0

    def test_entry_structure(self):
        """Every entry must have the four required keys."""
        from services.aop_discovery_service import _hardcoded_aop_list

        for entry in _hardcoded_aop_list():
            assert "aop_id" in entry
            assert "title" in entry
            assert "ke_count" in entry
            assert "mapped_ke_count" in entry

    def test_only_enabled_entries(self):
        """Only enabled CASE_STUDY_AOPS entries with an id are included."""
        from config import Config
        from services.aop_discovery_service import _hardcoded_aop_list

        expected_ids = {
            v["id"]
            for v in Config.CASE_STUDY_AOPS.values()
            if v.get("enabled") and v.get("id")
        }
        actual_ids = {e["aop_id"] for e in _hardcoded_aop_list()}
        assert actual_ids == expected_ids

    def test_mapped_ke_count_nonzero(self):
        """Hardcoded AOPs get mapped_ke_count=1 to appear in defaults."""
        from services.aop_discovery_service import _hardcoded_aop_list

        for entry in _hardcoded_aop_list():
            assert entry["mapped_ke_count"] >= 1


# ---------------------------------------------------------------------------
# Unit: fetch_mapped_ke_ids_from_builder
# ---------------------------------------------------------------------------

class TestFetchMappedKeIds:
    """Tests for fetch_mapped_ke_ids_from_builder."""

    def test_empty_url_returns_empty_set(self):
        """Empty BUILDER_API_URL must return an empty set without error."""
        from services.aop_discovery_service import fetch_mapped_ke_ids_from_builder

        config = _make_config(builder_url="")
        result = fetch_mapped_ke_ids_from_builder(config)
        assert isinstance(result, set)
        assert len(result) == 0

    def test_with_url_returns_normalised_ids(self):
        """When the Builder API returns records, KE IDs are normalised."""
        from services.aop_discovery_service import fetch_mapped_ke_ids_from_builder

        fake_records = [
            {"ke_id": "KE 55", "pathway_id": "WP1234"},
            {"ke_id": "KE 709", "pathway_id": "WP5678"},
        ]

        config = _make_config(builder_url="https://builder.example.com")

        # fetch_mapped_ke_ids_from_builder does a local import of api_service
        # inside the function body, so we patch at the api_service module level.
        import services.api_service as api_svc

        with patch.object(api_svc, "fetch_all_ke_wp_mappings", return_value=fake_records), \
             patch.object(api_svc, "_make_api_session", return_value=MagicMock()):
            result = fetch_mapped_ke_ids_from_builder(config)

        assert isinstance(result, set)
        assert "KE:55" in result
        assert "KE:709" in result


# ---------------------------------------------------------------------------
# Integration: three-tier fallback
# ---------------------------------------------------------------------------

class TestGetAopListFallback:
    """Tests for get_aop_list three-tier fallback logic."""

    def test_fallback_to_hardcoded_when_both_fail(self, tmp_path):
        """With no cache, no SPARQL and no Builder, the hardcoded list is returned.

        _make_config leaves BUILDER_API_URL empty, so the Builder tier returns
        nothing without touching the network and the chain reaches tier 4.
        """
        from services.aop_discovery_service import get_aop_list

        config = _make_config(cache_dir=str(tmp_path))

        # Make build_aop_list raise to simulate complete API failure
        with patch("services.aop_discovery_service.build_aop_list", side_effect=RuntimeError("SPARQL down")):
            result = get_aop_list(config)

        assert isinstance(result, list)
        assert len(result) > 0
        # All entries must have the standard keys
        for entry in result:
            assert "aop_id" in entry
            assert "mapped_ke_count" in entry

    def test_cache_hit_returns_immediately(self, tmp_path):
        """When the disk cache is primed, build_aop_list is not called."""
        import diskcache
        from services.aop_discovery_service import AOP_LIST_CACHE_KEY, get_aop_list

        config = _make_config(cache_dir=str(tmp_path))

        # Prime the disk cache manually
        fake_aops = [
            {"aop_id": "AOP:99", "title": "Cached AOP", "ke_count": 5, "mapped_ke_count": 2}
        ]
        dc = diskcache.FanoutCache(directory=str(tmp_path), shards=8)
        dc.set(AOP_LIST_CACHE_KEY, fake_aops)
        dc.close()

        build_called = []

        def _fake_build(cfg):
            build_called.append(True)
            return []

        with patch("services.aop_discovery_service.build_aop_list", side_effect=_fake_build):
            result = get_aop_list(config)

        assert result == fake_aops
        assert len(build_called) == 0, "build_aop_list should not be called on cache hit"

    def test_result_is_cached_after_build(self, tmp_path):
        """Successful build_aop_list result is stored in disk cache."""
        import diskcache
        from services.aop_discovery_service import AOP_LIST_CACHE_KEY, get_aop_list

        config = _make_config(cache_dir=str(tmp_path))

        fake_aops = [
            {"aop_id": "AOP:10", "title": "Live AOP", "ke_count": 3, "mapped_ke_count": 1}
        ]

        with patch("services.aop_discovery_service.build_aop_list", return_value=fake_aops):
            result = get_aop_list(config)

        assert result == fake_aops

        # Verify cache was written
        dc = diskcache.FanoutCache(directory=str(tmp_path), shards=8)
        cached = dc.get(AOP_LIST_CACHE_KEY)
        dc.close()
        assert cached == fake_aops


# ---------------------------------------------------------------------------
# Flask endpoint: /api/aops
# ---------------------------------------------------------------------------

class TestApiAopsEndpoint:
    """Tests for GET /api/aops."""

    @pytest.fixture
    def client(self):
        """Provide a CSRF-disabled Flask test client."""
        from app import app as flask_app

        flask_app.config["TESTING"] = True
        flask_app.config["WTF_CSRF_ENABLED"] = False
        flask_app.config["SECRET_KEY"] = "test-secret"

        with flask_app.test_client() as c:
            with flask_app.app_context():
                yield c

    def _mock_aop_list(self):
        """Return a small deterministic AOP list for tests."""
        return [
            {"aop_id": "AOP:33", "title": "Kidney toxicity induced by 5HT2C activation", "ke_count": 6, "mapped_ke_count": 4},
            {"aop_id": "AOP:105", "title": "Alpha2u kidney tubule adenoma", "ke_count": 8, "mapped_ke_count": 3},
            {"aop_id": "AOP:300", "title": "Liver steatosis via PPAR gamma", "ke_count": 5, "mapped_ke_count": 0},
        ]

    def test_endpoint_returns_json_array(self, client):
        """GET /api/aops with no query returns a JSON array."""
        with patch("app.get_aop_list", return_value=self._mock_aop_list()):
            resp = client.get("/api/aops")

        assert resp.status_code == 200
        data = resp.get_json()
        assert isinstance(data, list)

    def test_empty_query_returns_mapped_only(self, client):
        """GET /api/aops (no q) returns only AOPs with mapped_ke_count > 0."""
        with patch("app.get_aop_list", return_value=self._mock_aop_list()):
            resp = client.get("/api/aops")

        data = resp.get_json()
        assert all(a["mapped_ke_count"] > 0 for a in data)

    def test_query_filters_by_title(self, client):
        """GET /api/aops?q=kidney filters by title substring (case-insensitive)."""
        with patch("app.get_aop_list", return_value=self._mock_aop_list()):
            resp = client.get("/api/aops?q=kidney")

        data = resp.get_json()
        assert len(data) >= 1
        for entry in data:
            assert "kidney" in entry["title"].lower() or "kidney" in entry["aop_id"].lower()

    def test_query_filters_by_aop_id(self, client):
        """GET /api/aops?q=AOP:105 returns the matching entry."""
        with patch("app.get_aop_list", return_value=self._mock_aop_list()):
            resp = client.get("/api/aops?q=AOP:105")

        data = resp.get_json()
        assert any(a["aop_id"] == "AOP:105" for a in data)

    def test_results_capped_at_50(self, client):
        """Results are capped at 50 items regardless of list size."""
        big_list = [
            {"aop_id": f"AOP:{i}", "title": f"AOP number {i}", "ke_count": 5, "mapped_ke_count": 1}
            for i in range(100)
        ]
        with patch("app.get_aop_list", return_value=big_list):
            resp = client.get("/api/aops?q=AOP number")

        data = resp.get_json()
        assert len(data) <= 50

    def test_service_failure_returns_503(self, client):
        """When get_aop_list raises, endpoint returns HTTP 503 with error key."""
        with patch("app.get_aop_list", side_effect=RuntimeError("SPARQL exploded")):
            resp = client.get("/api/aops")

        assert resp.status_code == 503
        body = resp.get_json()
        assert "error" in body


# ---------------------------------------------------------------------------
# data_service: default source routing
# ---------------------------------------------------------------------------

class TestLoadAopDataDefaultSparql:
    """Verify that load_aop_data() routes dynamic AOPs to SPARQL."""

    def test_unknown_aop_uses_sparql(self):
        """An AOP not in CASE_STUDY_AOPS must route to fetch_aop_ke_data_cached."""
        import pandas as pd
        from services.data_service import load_aop_data

        fake_result = (
            {"KE:999"},
            pd.DataFrame(columns=["AOP_ID", "Source_KE", "Target_KE", "KER_ID"]),
            {"KE:999": "MIE"},
            {"KE:999": "Fake KE"},
        )

        # The import is done inside load_aop_data: "from services.sparql_service import ..."
        # Patch at the sparql_service module level so it's picked up on the local import.
        with patch("services.sparql_service.fetch_aop_ke_data_cached", return_value=fake_result) as mock_sparql:
            result = load_aop_data("AOP:99999")

        mock_sparql.assert_called_once_with(["AOP:99999"])
        assert result[0] == {"KE:999"}

    def test_csv_aop_still_uses_csv(self):
        """An AOP with source='csv' in CASE_STUDY_AOPS must NOT call SPARQL."""
        import pandas as pd
        from services.data_service import load_aop_data

        # AOP:DEMO is configured with source='csv'
        fake_ke_set = {"KE:115"}
        fake_edges = pd.DataFrame(columns=["AOP_ID", "Source_KE", "Target_KE"])
        fake_type = {"KE:115": "MIE"}
        fake_title = {"KE:115": "Some KE"}

        with patch("services.data_service._load_aop_data_csv",
                   return_value=(fake_ke_set, fake_edges, fake_type, fake_title)) as mock_csv, \
             patch("services.sparql_service.fetch_aop_ke_data_cached") as mock_sparql:

            result = load_aop_data("AOP:DEMO")

        mock_csv.assert_called_once_with("AOP:DEMO")
        mock_sparql.assert_not_called()

    def test_sparql_aop_uses_sparql(self):
        """An AOP explicitly configured with source='sparql' routes to SPARQL."""
        import pandas as pd
        from services.data_service import load_aop_data

        fake_result = (
            {"KE:200"},
            pd.DataFrame(columns=["AOP_ID", "Source_KE", "Target_KE", "KER_ID"]),
            {"KE:200": "AO"},
            {"KE:200": "SPARQL KE"},
        )

        # AOP:33 is configured as source='sparql'
        with patch("services.sparql_service.fetch_aop_ke_data_cached", return_value=fake_result) as mock_sparql:
            result = load_aop_data("AOP:33")

        mock_sparql.assert_called_once()
        assert result[0] == {"KE:200"}


# ---------------------------------------------------------------------------
# Builder-sourced AOP list (GET /api/v1/aops)
# ---------------------------------------------------------------------------

class TestFetchAopsFromBuilder:
    """Tests for fetch_aops_from_builder."""

    @staticmethod
    def _payload(rows, next_url=None):
        return {"data": rows, "pagination": {"next": next_url}}

    def test_empty_url_returns_empty_list(self):
        """No Builder configured must not raise, and must not fabricate AOPs."""
        from services.aop_discovery_service import fetch_aops_from_builder

        assert fetch_aops_from_builder(_make_config(builder_url="")) == []

    def test_normalises_ids_and_shape(self):
        """Builder's "AOP 625" becomes the analyser's "AOP:625"."""
        from services.aop_discovery_service import fetch_aops_from_builder

        rows = [{
            "aop_id": "AOP 625",
            "aop_title": "Increased 11B-HSD leading to metabolic syndrome",
            "ke_count": 18,
            "mapped_ke_count": 18,
            "wikipathways_ke_count": 18,
            "go_ke_count": 2,
            "reactome_ke_count": 2,
        }]
        session = MagicMock()
        session.get.return_value.json.return_value = self._payload(rows)

        import services.api_service as api_svc
        with patch.object(api_svc, "_make_api_session", return_value=session):
            result = fetch_aops_from_builder(_make_config(builder_url="https://b.example.com"))

        assert result == [{
            "aop_id": "AOP:625",
            "title": "Increased 11B-HSD leading to metabolic syndrome",
            "ke_count": 18,
            "mapped_ke_count": 18,
        }]

    def test_follows_pagination(self):
        """All pages are collected, not just the first."""
        from services.aop_discovery_service import fetch_aops_from_builder

        page1 = self._payload(
            [{"aop_id": "AOP 1", "aop_title": "One", "ke_count": 3, "mapped_ke_count": 1}],
            next_url="https://b.example.com/api/v1/aops?page=2",
        )
        page2 = self._payload(
            [{"aop_id": "AOP 2", "aop_title": "Two", "ke_count": 4, "mapped_ke_count": 2}]
        )
        session = MagicMock()
        session.get.return_value.json.side_effect = [page1, page2]

        import services.api_service as api_svc
        with patch.object(api_svc, "_make_api_session", return_value=session):
            result = fetch_aops_from_builder(_make_config(builder_url="https://b.example.com"))

        assert [a["aop_id"] for a in result] == ["AOP:1", "AOP:2"]

    def test_builder_failure_returns_empty_list(self):
        """A Builder outage degrades to empty, never an exception.

        get_aop_list relies on this: the Builder tier must be skippable so the
        hardcoded tier can still answer.
        """
        from services.aop_discovery_service import fetch_aops_from_builder

        session = MagicMock()
        session.get.side_effect = RuntimeError("connection refused")

        import services.api_service as api_svc
        with patch.object(api_svc, "_make_api_session", return_value=session):
            result = fetch_aops_from_builder(_make_config(builder_url="https://b.example.com"))

        assert result == []

    def test_rows_without_aop_id_are_skipped(self):
        """A malformed row must not produce an entry with an empty ID."""
        from services.aop_discovery_service import fetch_aops_from_builder

        rows = [
            {"aop_id": "", "aop_title": "No ID", "ke_count": 1, "mapped_ke_count": 1},
            {"aop_id": "AOP 7", "aop_title": "Fine", "ke_count": 2, "mapped_ke_count": 1},
        ]
        session = MagicMock()
        session.get.return_value.json.return_value = self._payload(rows)

        import services.api_service as api_svc
        with patch.object(api_svc, "_make_api_session", return_value=session):
            result = fetch_aops_from_builder(_make_config(builder_url="https://b.example.com"))

        assert [a["aop_id"] for a in result] == ["AOP:7"]


class TestBuilderFallbackTier:
    """get_aop_list falls through to the Builder when SPARQL is unavailable."""

    def test_builder_list_used_when_sparql_fails(self, tmp_path):
        """A SPARQL outage yields the Builder's AOPs, not five hardcoded ones."""
        from services.aop_discovery_service import get_aop_list

        builder_aops = [
            {"aop_id": "AOP:624", "title": "A", "ke_count": 18, "mapped_ke_count": 18},
            {"aop_id": "AOP:130", "title": "B", "ke_count": 11, "mapped_ke_count": 11},
        ]
        config = _make_config(cache_dir=str(tmp_path))

        with patch("services.aop_discovery_service.build_aop_list",
                   side_effect=RuntimeError("SPARQL down")), \
             patch("services.aop_discovery_service.fetch_aops_from_builder",
                   return_value=builder_aops):
            result = get_aop_list(config)

        assert result == builder_aops

    def test_builder_list_cached_with_short_ttl(self, tmp_path):
        """The degraded list is cached, but must not outlive the outage.

        It omits every unmapped AOP, so caching it for the full week would keep
        a partial list in place long after SPARQL recovered.
        """
        import diskcache
        from services.aop_discovery_service import (
            AOP_LIST_BUILDER_CACHE_TTL, AOP_LIST_CACHE_KEY, AOP_LIST_CACHE_TTL,
            get_aop_list,
        )

        assert AOP_LIST_BUILDER_CACHE_TTL < AOP_LIST_CACHE_TTL

        builder_aops = [{"aop_id": "AOP:1", "title": "A", "ke_count": 3, "mapped_ke_count": 2}]
        config = _make_config(cache_dir=str(tmp_path))

        with patch("services.aop_discovery_service.build_aop_list",
                   side_effect=RuntimeError("SPARQL down")), \
             patch("services.aop_discovery_service.fetch_aops_from_builder",
                   return_value=builder_aops):
            get_aop_list(config)

        # Must match get_reference_cache's FanoutCache sharding — a plain
        # Cache reads only shard 0 and would report a miss.
        dc = diskcache.FanoutCache(directory=str(tmp_path), shards=8)
        try:
            assert dc.get(AOP_LIST_CACHE_KEY) == builder_aops
        finally:
            dc.close()

    def test_hardcoded_tier_still_reachable_when_builder_empty(self, tmp_path):
        """Both SPARQL and Builder down still leaves a usable picker."""
        from services.aop_discovery_service import get_aop_list

        config = _make_config(cache_dir=str(tmp_path))

        with patch("services.aop_discovery_service.build_aop_list",
                   side_effect=RuntimeError("SPARQL down")), \
             patch("services.aop_discovery_service.fetch_aops_from_builder",
                   return_value=[]):
            result = get_aop_list(config)

        assert len(result) > 0
        assert all("aop_id" in a for a in result)
