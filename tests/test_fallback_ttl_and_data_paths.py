"""Fallback cache lifetime and cwd-independent data paths.

A bundled-CSV fallback cached for the full CACHE_TTL turned a builder outage of
seconds into an hour of analyses on the bundled files (2026-10-04: the analyser
restarted 3s before the builder). Fallback entries now expire after
Config.FALLBACK_CACHE_TTL, and the cache-age line still dates them correctly.

Reference data was read through paths relative to the working directory, which
only worked because the container happens to start in /app.
"""
import datetime as dt
import os
import subprocess
import sys
from unittest.mock import patch

import pytest

import app
from config import Config
from helpers import load_reference_sets

# conftest stubs both of these for every test; keep the real ones, captured at
# collection, because here the real file access is the subject.
_real_validate_data_files = Config.validate_data_files
_real_load_reference_sets = load_reference_sets

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class _RecordingCache:
    """A reference-cache double that remembers each write's expiry."""

    def __init__(self, store=None, expiry=None):
        self.store = dict(store or {})
        self.expiry = dict(expiry or {})
        self.writes = {}

    def get(self, key, expire_time=False):
        value = self.store.get(key)
        return (value, self.expiry.get(key)) if expire_time else value

    def set(self, key, value, expire=None):
        self.store[key] = value
        self.writes[key] = expire


@pytest.fixture(autouse=True)
def _clear_recorded_times():
    app._CACHE_FILL_TIMES.clear()
    yield
    app._CACHE_FILL_TIMES.clear()


class TestFallbackCacheTTL:

    def test_default_is_five_minutes_and_shorter_than_live(self):
        assert Config.FALLBACK_CACHE_TTL == 300
        assert Config.FALLBACK_CACHE_TTL < Config.CACHE_TTL

    def test_env_overrides_fallback_ttl(self):
        env = dict(os.environ, FALLBACK_CACHE_TTL="42")
        out = subprocess.run(
            [sys.executable, "-c",
             "from config import Config; print(Config.FALLBACK_CACHE_TTL)"],
            cwd=REPO_ROOT, env=env, capture_output=True, text=True, check=True,
        )
        assert out.stdout.strip() == "42"

    def test_csv_fallback_is_cached_for_the_fallback_ttl(self, monkeypatch):
        cache = _RecordingCache()
        monkeypatch.setattr(app, "_reference_cache", cache)
        with patch.object(app, "fetch_reference_sets_from_api",
                          side_effect=RuntimeError("connection refused")), \
             patch.object(app, "load_reference_sets",
                          return_value={"KE:1": {"A"}}):
            _, source, _, _ = app._load_wikipathways_reference_sets("all")

        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        assert source == "csv"
        assert cache.writes[key] == Config.FALLBACK_CACHE_TTL

    def test_live_api_load_keeps_the_full_ttl(self, monkeypatch):
        cache = _RecordingCache()
        monkeypatch.setattr(app, "_reference_cache", cache)
        with patch.object(app, "fetch_reference_sets_from_api",
                          return_value=({"KE:1": {"A"}}, [])):
            _, source, _, _ = app._load_wikipathways_reference_sets("all")

        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        assert source == "api"
        assert cache.writes[key] == Config.CACHE_TTL


class TestCacheAgeUsesTheEntrysTTL:
    """The fill time is expiry minus the TTL the entry was *written* with."""

    FILLED = dt.datetime(2026, 10, 4, 8, 30, tzinfo=dt.timezone.utc)

    def _load_cached(self, monkeypatch, source, ttl):
        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        cache = _RecordingCache(
            {key: ({"KE:1": {"A"}}, source, [], {})},
            {key: self.FILLED.timestamp() + ttl},
        )
        monkeypatch.setattr(app, "_reference_cache", cache)
        app._load_wikipathways_reference_sets("all")
        return app._cache_fill_time("WikiPathways", "all")

    def test_csv_entry_age_is_derived_from_the_fallback_ttl(self, monkeypatch):
        assert self._load_cached(
            monkeypatch, "csv", Config.FALLBACK_CACHE_TTL
        ) == "2026-10-04 08:30 UTC"

    def test_api_entry_age_is_derived_from_the_live_ttl(self, monkeypatch):
        assert self._load_cached(
            monkeypatch, "api", Config.CACHE_TTL
        ) == "2026-10-04 08:30 UTC"


class TestDataPathsIgnoreTheWorkingDirectory:

    @pytest.fixture(autouse=True)
    def _elsewhere(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)

    def test_required_data_files_validate(self):
        assert _real_validate_data_files() is True

    def test_aop_csv_loader(self):
        from services.data_service import _load_aop_data_csv

        ke_list, _, ke_type_map, ke_title_map = _load_aop_data_csv("AOP:DEMO")
        assert ke_list
        assert ke_type_map and ke_title_map

    def test_ke_wp_records_csv_default_path(self):
        from services.api_service import load_ke_wp_records_csv

        assert load_ke_wp_records_csv()

    def test_load_reference_sets_default_paths(self):
        assert _real_load_reference_sets()

    def test_wikipathways_csv_fallback(self, monkeypatch):
        monkeypatch.setattr(app, "_reference_cache", _RecordingCache())
        with patch.object(app, "fetch_reference_sets_from_api",
                          side_effect=RuntimeError("connection refused")):
            reference_sets, source, _, _ = app._load_wikipathways_reference_sets("all")
        assert source == "csv"
        assert reference_sets

    def test_preview_resolves_demo_file(self, flask_client, tmp_path, monkeypatch):
        uploads = tmp_path / "uploads"
        uploads.mkdir()
        monkeypatch.setitem(app.app.config, "UPLOAD_FOLDER", str(uploads))

        response = flask_client.post("/preview", data={
            "demo_file": "GSE90122_TO90137.tsv",
            "dataset_id": "TEST001",
            "stressor": "Test Chemical",
            "owner": "Test User",
        })

        assert response.status_code == 200, response.data[:200]
        assert (uploads / "GSE90122_TO90137.tsv").is_file()
