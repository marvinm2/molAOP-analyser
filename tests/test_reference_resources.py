"""Reference-set loading in app: multi-resource merging (#55) and the
minimum KE-mapping confidence threshold + cache keying (#60)."""
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

import app
from helpers import load_reference_sets
from services.api_service import parse_gmt_reference_sets


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

        def _boom(resource, min_confidence="all"):
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


class TestReferenceSetConfidenceFiltering:
    """Issue #60: how the threshold reshapes the KE gene sets themselves.

    Exercises helpers.load_reference_sets end to end against small on-disk CSVs
    so the assertions are about real gene sets, not mocks.
    """

    @pytest.fixture
    def gene_csvs(self, tmp_path):
        """WikiPathways gene edges + node annotations shared by these tests."""
        wp_gene = tmp_path / "edges.csv"
        wp_gene.write_text(
            "WPID,gene_id\n"
            "WP1,1\nWP1,2\n"   # high-confidence pathway for KE:1
            "WP2,3\n"          # low-confidence pathway for KE:1
            "WP3,4\n"          # medium-confidence pathway for KE:2
            "WP4,5\n"          # low-confidence pathway for KE:3 (its only one)
        )
        node = tmp_path / "nodes.csv"
        node.write_text(
            "GeneID,GeneName\n1,AAA\n2,BBB\n3,CCC\n4,DDD\n5,EEE\n"
        )
        return str(wp_gene), str(node)

    @staticmethod
    def _mappings(with_confidence=True):
        rows = [
            ("KE:1", "WP1", "high"),
            ("KE:1", "WP2", "low"),
            ("KE:2", "WP3", "medium"),
            ("KE:3", "WP4", "low"),
        ]
        data = {"KE_ID": [r[0] for r in rows], "WP_ID": [r[1] for r in rows]}
        if with_confidence:
            data["confidence_level"] = [r[2] for r in rows]
        return pd.DataFrame(data)

    def test_all_keeps_every_mapping(self, gene_csvs):
        wp_gene, node = gene_csvs
        sets = load_reference_sets(
            wp_gene_path=wp_gene, node_path=node,
            ke_wp_df=self._mappings(), min_confidence="all",
        )
        assert sets["KE:1"] == {"AAA", "BBB", "CCC"}
        assert sets["KE:2"] == {"DDD"}
        assert sets["KE:3"] == {"EEE"}

    def test_ke_keeps_high_mapping_and_loses_low_one(self, gene_csvs):
        """Per-KE multi-mapping filtering: KE:1 keeps WP1's genes, drops WP2's."""
        wp_gene, node = gene_csvs
        sets = load_reference_sets(
            wp_gene_path=wp_gene, node_path=node,
            ke_wp_df=self._mappings(), min_confidence="high",
        )
        assert sets["KE:1"] == {"AAA", "BBB"}
        assert "CCC" not in sets["KE:1"]

    def test_ke_with_only_below_threshold_mappings_drops_out(self, gene_csvs):
        """KE:3's single mapping is Low -> the KE is absent from enrichment."""
        wp_gene, node = gene_csvs
        sets = load_reference_sets(
            wp_gene_path=wp_gene, node_path=node,
            ke_wp_df=self._mappings(), min_confidence="medium",
        )
        assert "KE:3" not in sets
        # KE:2's Medium mapping survives at the Medium threshold.
        assert sets["KE:2"] == {"DDD"}
        assert "KE:2" not in load_reference_sets(
            wp_gene_path=wp_gene, node_path=node,
            ke_wp_df=self._mappings(), min_confidence="high",
        )

    def test_missing_confidence_column_is_a_no_op(self, gene_csvs):
        """CSV fallback shape (KE_ID/WP_ID only): the filter must not empty it."""
        wp_gene, node = gene_csvs
        sets = load_reference_sets(
            wp_gene_path=wp_gene, node_path=node,
            ke_wp_df=self._mappings(with_confidence=False), min_confidence="high",
        )
        assert sets["KE:1"] == {"AAA", "BBB", "CCC"}
        assert sets["KE:3"] == {"EEE"}

    def test_blank_confidence_values_are_kept(self, gene_csvs):
        """Rows with null/blank confidence survive the strictest threshold."""
        wp_gene, node = gene_csvs
        df = self._mappings()
        df["confidence_level"] = [None, "", "  ", "unknown"]
        sets = load_reference_sets(
            wp_gene_path=wp_gene, node_path=node,
            ke_wp_df=df, min_confidence="high",
        )
        assert sets["KE:1"] == {"AAA", "BBB", "CCC"}
        assert sets["KE:3"] == {"EEE"}

    def test_titlecase_confidence_ranks_like_lowercase(self, gene_csvs):
        """Live API sends 'high'; the API doc shows 'High'. Same gene sets either way."""
        wp_gene, node = gene_csvs
        lower = self._mappings()
        title = self._mappings()
        title["confidence_level"] = title["confidence_level"].str.capitalize()
        kwargs = dict(wp_gene_path=wp_gene, node_path=node, min_confidence="high")
        assert (
            load_reference_sets(ke_wp_df=lower, **kwargs)
            == load_reference_sets(ke_wp_df=title, **kwargs)
        )


class TestGmtResourcesUnaffectedByConfidence:
    """Issue #60: the confidence-free GMT resources stay fully intact."""

    GMT = (
        "KE1115_Increase_Reactive_oxygen_species_GO:0072593\tROS\tSOD1\tCAT\n"
        "KE1392_Oxidative_Stress_R-HSA-9818027\tOxidative stress\tNFE2L2\n"
    )

    @pytest.mark.parametrize("threshold", ["all", "medium", "high"])
    def test_gmt_sets_identical_at_every_threshold(self, threshold, monkeypatch):
        store = {}
        fake_cache = MagicMock()
        fake_cache.get.side_effect = lambda key: store.get(key)
        fake_cache.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake_cache)

        parsed = parse_gmt_reference_sets(self.GMT)
        with patch.object(app, "fetch_gmt_reference_sets", return_value=parsed):
            sets, _ = app.load_cached_reference_sets(["GO_BP"], min_confidence=threshold)
        assert sets == {"KE:1115": {"SOD1", "CAT"}, "KE:1392": {"NFE2L2"}}


class TestMinConfidenceThreading:
    """Issue #60: the confidence threshold reaches every resource loader."""

    def test_threshold_forwarded_to_resource_loaders(self):
        wp = {"KE:1": {"A"}}
        go = {"KE:1": {"B"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api")) as wp_loader, \
             patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "api")) as gmt_loader:
            app.load_cached_reference_sets(["WikiPathways", "GO_BP"], min_confidence="high")
        assert wp_loader.call_args.args[0] == "high"
        assert gmt_loader.call_args.args[1] == "high"

    def test_default_threshold_is_all(self):
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "csv")) as wp_loader:
            app.load_cached_reference_sets(["WikiPathways"])
        assert wp_loader.call_args.args[0] == "all"

    def test_junk_threshold_falls_back_to_all(self):
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "csv")) as wp_loader:
            app.load_cached_reference_sets(["WikiPathways"], min_confidence="'; DROP TABLE --")
        assert wp_loader.call_args.args[0] == "all"


class TestConfidenceScopedCacheKeys:
    """Issue #60: reference-set cache entries must never leak across thresholds."""

    def test_keys_differ_per_threshold(self):
        keys = {
            app._confidence_cache_key(app.REFERENCE_CACHE_KEY, level)
            for level in ("all", "medium", "high")
        }
        assert len(keys) == 3

    def test_gmt_keys_are_also_threshold_scoped(self):
        base = app._GMT_RESOURCE_CACHE_KEYS["GO_BP"]
        assert app._confidence_cache_key(base, "all") != app._confidence_cache_key(base, "high")

    def test_gmt_and_wikipathways_keys_never_collide(self):
        assert (
            app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "high")
            != app._confidence_cache_key(app._GMT_RESOURCE_CACHE_KEYS["Reactome"], "high")
        )

    def test_high_request_does_not_read_the_all_cache_entry(self, monkeypatch):
        """An 'all' entry in the disk cache must not satisfy a 'high' request."""
        store = {app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all"): ({"KE:1": {"A", "B"}}, "api")}
        fake_cache = MagicMock()
        fake_cache.get.side_effect = lambda key: store.get(key)
        fake_cache.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake_cache)

        with patch.object(app, "fetch_reference_sets_from_api", return_value={"KE:1": {"A"}}) as fetch:
            sets, source = app._load_wikipathways_reference_sets("high")

        # Cache miss for 'high' -> a fresh fetch at that threshold, not the 'all' set.
        fetch.assert_called_once()
        assert fetch.call_args.kwargs["min_confidence"] == "high"
        assert sets == {"KE:1": {"A"}}
        assert source == "api"
        # The 'all' entry is untouched and the 'high' result is cached separately.
        assert store[app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")][0] == {"KE:1": {"A", "B"}}
        assert store[app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "high")][0] == {"KE:1": {"A"}}

    def test_all_request_does_not_read_the_high_cache_entry(self, monkeypatch):
        """And the reverse: a 'high' entry must not satisfy an 'all' request."""
        store = {app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "high"): ({"KE:1": {"A"}}, "api")}
        fake_cache = MagicMock()
        fake_cache.get.side_effect = lambda key: store.get(key)
        fake_cache.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake_cache)

        with patch.object(app, "fetch_reference_sets_from_api", return_value={"KE:1": {"A", "B"}}):
            sets, _ = app._load_wikipathways_reference_sets("all")
        assert sets == {"KE:1": {"A", "B"}}

    def test_gmt_loader_caches_per_threshold(self, monkeypatch):
        store = {}
        fake_cache = MagicMock()
        fake_cache.get.side_effect = lambda key: store.get(key)
        fake_cache.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake_cache)

        # GMT exports carry no confidence field, so the sets are identical --
        # but they must still be stored under distinct, threshold-scoped keys.
        with patch.object(app, "fetch_gmt_reference_sets", return_value={"KE:1115": {"SOD1"}}) as fetch:
            app._load_gmt_resource_reference_sets("GO_BP", "all")
            app._load_gmt_resource_reference_sets("GO_BP", "high")
        assert fetch.call_count == 2
        assert len(store) == 2
