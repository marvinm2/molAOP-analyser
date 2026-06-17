"""Tests for multi-resource reference-set merging in app.load_cached_reference_sets (#55)."""
from unittest.mock import patch

import app


class TestLoadCachedReferenceSets:
    """Merge + graceful-degradation behaviour across gene-set resources."""

    def test_default_is_wikipathways_only(self):
        wp = {"KE:1": {"A", "B"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "csv")) as wp_loader, \
             patch.object(app, "_load_gmt_resource_reference_sets") as gmt_loader:
            sets, source = app.load_cached_reference_sets()
        wp_loader.assert_called_once()
        gmt_loader.assert_not_called()
        assert sets == {"KE:1": {"A", "B"}}
        # data_source reflects the WikiPathways component for picker compatibility.
        assert source == "csv"

    def test_unions_genes_across_resources(self):
        wp = {"KE:1": {"A", "B"}}
        go = {"KE:1": {"B", "C"}, "KE:2": {"D"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api")), \
             patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "api")):
            sets, source = app.load_cached_reference_sets(["WikiPathways", "GO_BP"])
        assert sets == {"KE:1": {"A", "B", "C"}, "KE:2": {"D"}}
        assert source == "api"  # WikiPathways source still drives data_source

    def test_does_not_mutate_source_sets(self):
        wp = {"KE:1": {"A"}}
        go = {"KE:1": {"B"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api")), \
             patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "api")):
            app.load_cached_reference_sets(["WikiPathways", "GO_BP"])
        # Cached per-resource sets are untouched by the merge.
        assert wp == {"KE:1": {"A"}}
        assert go == {"KE:1": {"B"}}

    def test_gmt_failure_degrades_to_working_resources(self):
        wp = {"KE:1": {"A"}}

        def _boom(resource):
            raise RuntimeError("builder GMT export unavailable")

        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api")), \
             patch.object(app, "_load_gmt_resource_reference_sets", side_effect=_boom):
            sets, source = app.load_cached_reference_sets(["WikiPathways", "Reactome"])
        # Reactome skipped; WikiPathways retained -> analysis still runs.
        assert sets == {"KE:1": {"A"}}
        assert source == "api"

    def test_unknown_resource_falls_back_to_default(self):
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "csv")) as wp_loader:
            sets, source = app.load_cached_reference_sets(["bogus"])
        wp_loader.assert_called_once()
        assert sets == {"KE:1": {"A"}}

    def test_gmt_only_selection_sets_gmt_source(self):
        go = {"KE:1": {"X"}}
        with patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "api")):
            sets, source = app.load_cached_reference_sets(["GO_BP"])
        assert sets == {"KE:1": {"X"}}
        # No WikiPathways component -> data_source is 'gmt' (picker falls back to CSV).
        assert source == "gmt"
