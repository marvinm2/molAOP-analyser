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
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "csv", [], {})) as wp_loader, \
             patch.object(app, "_load_gmt_resource_reference_sets") as gmt_loader:
            sets, source, _ = app.load_cached_reference_sets()
        wp_loader.assert_called_once()
        gmt_loader.assert_not_called()
        assert sets == {"KE:1": {"A", "B"}}
        # data_source reflects the WikiPathways component for picker compatibility.
        assert source == "csv"

    def test_unions_genes_across_resources(self):
        wp = {"KE:1": {"A", "B"}}
        go = {"KE:1": {"B", "C"}, "KE:2": {"D"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api", [], {})), \
             patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "api")):
            sets, source, _ = app.load_cached_reference_sets(["WikiPathways", "GO_BP"])
        assert sets == {"KE:1": {"A", "B", "C"}, "KE:2": {"D"}}
        assert source == "api"  # WikiPathways source still drives data_source

    def test_does_not_mutate_source_sets(self):
        wp = {"KE:1": {"A"}}
        go = {"KE:1": {"B"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api", [], {})), \
             patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "api")):
            app.load_cached_reference_sets(["WikiPathways", "GO_BP"])
        # Cached per-resource sets are untouched by the merge.
        assert wp == {"KE:1": {"A"}}
        assert go == {"KE:1": {"B"}}

    def test_gmt_failure_degrades_to_working_resources(self):
        wp = {"KE:1": {"A"}}

        def _boom(resource, min_confidence="all"):
            raise RuntimeError("builder GMT export unavailable")

        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api", [], {})), \
             patch.object(app, "_load_gmt_resource_reference_sets", side_effect=_boom):
            sets, source, _ = app.load_cached_reference_sets(["WikiPathways", "Reactome"])
        # Reactome skipped; WikiPathways retained -> analysis still runs.
        assert sets == {"KE:1": {"A"}}
        assert source == "api"

    def test_unknown_resource_falls_back_to_default(self):
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "csv", [], {})) as wp_loader:
            sets, source, _ = app.load_cached_reference_sets(["bogus"])
        wp_loader.assert_called_once()
        assert sets == {"KE:1": {"A"}}

    def test_gmt_only_selection_sets_gmt_source(self):
        go = {"KE:1": {"X"}}
        with patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "api")):
            sets, source, _ = app.load_cached_reference_sets(["GO_BP"])
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
        fake_cache.get.side_effect = lambda key, expire_time=False: (
            (store.get(key), None) if expire_time else store.get(key))
        fake_cache.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake_cache)

        parsed = parse_gmt_reference_sets(self.GMT)
        with patch.object(app, "fetch_gmt_reference_sets", return_value=parsed):
            sets, _, _ = app.load_cached_reference_sets(["GO_BP"], min_confidence=threshold)
        assert sets == {"KE:1115": {"SOD1", "CAT"}, "KE:1392": {"NFE2L2"}}


class TestMinConfidenceThreading:
    """Issue #60: the confidence threshold reaches every resource loader."""

    def test_threshold_forwarded_to_resource_loaders(self):
        wp = {"KE:1": {"A"}}
        go = {"KE:1": {"B"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api", [], {})) as wp_loader, \
             patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "api")) as gmt_loader:
            app.load_cached_reference_sets(["WikiPathways", "GO_BP"], min_confidence="high")
        assert wp_loader.call_args.args[0] == "high"
        assert gmt_loader.call_args.args[1] == "high"

    def test_default_threshold_is_all(self):
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "csv", [], {})) as wp_loader:
            app.load_cached_reference_sets(["WikiPathways"])
        assert wp_loader.call_args.args[0] == "all"

    def test_junk_threshold_falls_back_to_all(self):
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "csv", [], {})) as wp_loader:
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
        fake_cache.get.side_effect = lambda key, expire_time=False: (
            (store.get(key), None) if expire_time else store.get(key))
        fake_cache.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake_cache)

        with patch.object(app, "fetch_reference_sets_from_api", return_value=({"KE:1": {"A"}}, [])) as fetch:
            sets, source, _, _ = app._load_wikipathways_reference_sets("high")

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
        fake_cache.get.side_effect = lambda key, expire_time=False: (
            (store.get(key), None) if expire_time else store.get(key))
        fake_cache.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake_cache)

        with patch.object(app, "fetch_reference_sets_from_api", return_value=({"KE:1": {"A", "B"}}, [])):
            sets, _, _, _ = app._load_wikipathways_reference_sets("all")
        assert sets == {"KE:1": {"A", "B"}}

    def test_gmt_loader_caches_per_threshold(self, monkeypatch):
        store = {}
        fake_cache = MagicMock()
        fake_cache.get.side_effect = lambda key, expire_time=False: (
            (store.get(key), None) if expire_time else store.get(key))
        fake_cache.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake_cache)

        # GMT exports carry no confidence field, so the sets are identical --
        # but they must still be stored under distinct, threshold-scoped keys.
        with patch.object(app, "fetch_gmt_reference_sets", return_value={"KE:1115": {"SOD1"}}) as fetch:
            app._load_gmt_resource_reference_sets("GO_BP", "all")
            app._load_gmt_resource_reference_sets("GO_BP", "high")
        assert fetch.call_count == 2
        assert len(store) == 2


class TestCachePreservesOriginalSource:
    """Issue #68: a cache hop must not erase where the gene sets came from."""

    @staticmethod
    def _fake_cache(monkeypatch, store):
        fake_cache = MagicMock()
        fake_cache.get.side_effect = lambda key, expire_time=False: (
            (store.get(key), None) if expire_time else store.get(key))
        fake_cache.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake_cache)
        return fake_cache

    def test_cached_csv_fallback_still_reads_as_csv(self, monkeypatch):
        """The reproducibility hole: sets built during a Builder outage were
        reported as a plain cache hit one hop later, indistinguishable from a
        live API response."""
        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        self._fake_cache(monkeypatch, {key: ({"KE:1": {"A"}}, "csv")})

        sets, source, _, _ = app._load_wikipathways_reference_sets("all")
        assert sets == {"KE:1": {"A"}}
        assert source == "cache(csv)"

    def test_cached_api_response_reads_as_cached_api(self, monkeypatch):
        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        self._fake_cache(monkeypatch, {key: ({"KE:1": {"A"}}, "api")})

        _, source, _, _ = app._load_wikipathways_reference_sets("all")
        assert source == "cache(api)"

    def test_a_cached_csv_fallback_is_still_warned_about(self, monkeypatch):
        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        self._fake_cache(monkeypatch, {key: ({"KE:1": {"A"}}, "csv")})

        _, _, resolution = app.load_cached_reference_sets(["WikiPathways"])
        assert any("bundled" in w for w in app.resource_resolution_warnings(resolution))


class TestResourceResolutionRecording:
    """Issue #68: record what actually resolved, not just what was requested."""

    def test_loaded_resources_are_recorded_with_their_source(self):
        wp = {"KE:1": {"A"}}
        go = {"KE:1": {"B"}, "KE:2": {"C"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api", [], {})), \
             patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "cache(api)")):
            _, _, resolution = app.load_cached_reference_sets(["WikiPathways", "GO_BP"])

        by_resource = {e["resource"]: e for e in resolution}
        assert by_resource["WikiPathways"]["status"] == "loaded"
        assert by_resource["WikiPathways"]["source"] == "api"
        assert by_resource["WikiPathways"]["ke_count"] == 1
        assert by_resource["GO_BP"]["source"] == "cache(api)"
        assert by_resource["GO_BP"]["ke_count"] == 2

    def test_skipped_resource_is_recorded_and_warned_about(self):
        """The failure mode this issue exists for: a resource silently absent."""
        wp = {"KE:1": {"A"}}

        def _boom(resource, min_confidence="all"):
            raise RuntimeError("builder GMT export unavailable")

        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api", [], {})), \
             patch.object(app, "_load_gmt_resource_reference_sets", side_effect=_boom):
            _, _, resolution = app.load_cached_reference_sets(["WikiPathways", "Reactome"])

        skipped = [e for e in resolution if e["status"] == "skipped"]
        assert [e["resource"] for e in skipped] == ["Reactome"]
        assert "unavailable" in skipped[0]["error"]

        warnings = app.resource_resolution_warnings(resolution)
        assert any("Reactome" in w and "left out" in w for w in warnings)

    def test_csv_fallback_is_disclosed(self):
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "csv", [], {})):
            _, _, resolution = app.load_cached_reference_sets(["WikiPathways"])

        warnings = app.resource_resolution_warnings(resolution)
        assert any("bundled" in w for w in warnings)

    def test_clean_run_produces_no_warnings(self):
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api", [], {})):
            _, _, resolution = app.load_cached_reference_sets(["WikiPathways"])
        assert app.resource_resolution_warnings(resolution, "high") == []

    def test_description_reads_as_a_sentence(self):
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "cache(csv)", [], {})):
            _, _, resolution = app.load_cached_reference_sets(["WikiPathways"])
        assert app.describe_resource_resolution(resolution) == (
            "WikiPathways (bundled reference files, cached)"
        )

    def test_legacy_records_describe_as_empty(self):
        """Runs stored before the resolution existed omit the line, not guess."""
        assert app.describe_resource_resolution([]) == ''
        assert app.resource_resolution_warnings([]) == []
        assert app._parse_resource_resolution(None) == []
        assert app._parse_resource_resolution('not json') == []


class TestConfidenceApplicabilityIsRecorded:
    """Issue #67: the threshold only bites where mappings carry a confidence."""

    def test_gmt_resources_are_marked_unfiltered(self):
        wp = {"KE:1": {"A"}}
        go = {"KE:1": {"B"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api", [], {})), \
             patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "api")):
            _, _, resolution = app.load_cached_reference_sets(
                ["WikiPathways", "GO_BP"], min_confidence="high"
            )
        by_resource = {e["resource"]: e for e in resolution}
        assert by_resource["WikiPathways"]["confidence_applied"] is True
        assert by_resource["GO_BP"]["confidence_applied"] is False

        warnings = app.resource_resolution_warnings(resolution, "high")
        assert any("GO_BP" in w and "confidence field" in w for w in warnings)

    def test_csv_fallback_cannot_apply_the_threshold_either(self):
        """KE-WP.csv has no confidence column, so 'high' is a no-op there too."""
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "csv", [], {})):
            _, _, resolution = app.load_cached_reference_sets(
                ["WikiPathways"], min_confidence="high"
            )
        assert resolution[0]["confidence_applied"] is False

    def test_no_applicability_warning_at_the_default_threshold(self):
        """'All mappings' filters nothing, so there is nothing to caveat."""
        wp = {"KE:1": {"A"}}
        go = {"KE:1": {"B"}}
        with patch.object(app, "_load_wikipathways_reference_sets", return_value=(wp, "api", [], {})), \
             patch.object(app, "_load_gmt_resource_reference_sets", return_value=(go, "api")):
            _, _, resolution = app.load_cached_reference_sets(
                ["WikiPathways", "GO_BP"], min_confidence="all"
            )
        assert app.resource_resolution_warnings(resolution, "all") == []


# ---------------------------------------------------------------------------
# Issue #79: a mapped pathway the bundled CSV does not know about
# ---------------------------------------------------------------------------

class TestPathwayResolutionGap:
    """load_reference_sets must not silently drop unresolvable pathways.

    Reproduces the reported case: KE 1115 maps only to WP5477, which is curated
    in the Builder but absent from data/edges_wpid_to_gene.csv. Before the fix
    the inner join dropped it without trace, so the KE was reported as having no
    gene set mapped — the opposite of what had happened.
    """

    KE_WP = pd.DataFrame([
        {"KE_ID": "KE:1115", "WP_ID": "WP5477"},   # only in the Builder
        {"KE_ID": "KE:177", "WP_ID": "WP4324"},    # in the bundled CSV
    ])

    def _csv_files(self, tmp_path):
        """Write a minimal stand-in for the two bundled CSVs."""
        wp_gene = tmp_path / "edges.csv"
        wp_gene.write_text("WPID,gene_id,edge_id\nWP4324,111,e1\n")
        nodes = tmp_path / "nodes.csv"
        nodes.write_text("GeneID,GeneName\n111,TP53\n")
        return str(wp_gene), str(nodes)

    def test_pathway_absent_from_csv_is_reported_not_dropped(self, tmp_path):
        wp_gene_path, node_path = self._csv_files(tmp_path)
        unresolved = []

        sets = load_reference_sets(
            ke_wp_df=self.KE_WP.copy(),
            wp_gene_path=wp_gene_path,
            node_path=node_path,
            unresolved_out=unresolved,
        )

        assert unresolved == ["WP5477"]
        # The KE really does end up with no genes — the point is that we say so.
        assert "KE:1115" not in sets
        assert sets["KE:177"] == {"TP53"}

    def test_builder_membership_resolves_the_gap(self, tmp_path):
        """With the Builder's map, the same KE gets its genes and nothing is unresolved."""
        wp_gene_path, node_path = self._csv_files(tmp_path)
        unresolved = []

        sets = load_reference_sets(
            ke_wp_df=self.KE_WP.copy(),
            wp_gene_path=wp_gene_path,
            node_path=node_path,
            wp_gene_map={"WP5477": {"SOD1", "CAT"}},
            unresolved_out=unresolved,
        )

        assert unresolved == []
        assert sets["KE:1115"] == {"SOD1", "CAT"}
        # The CSV still covers pathways the Builder did not answer for.
        assert sets["KE:177"] == {"TP53"}

    def test_builder_membership_replaces_rather_than_unions_the_snapshot(self, tmp_path):
        """When both sources know a pathway, the Builder's membership wins.

        Unioning would resurrect genes removed upstream, which is worse than
        either source alone: the result would match no version of the pathway.
        """
        wp_gene_path, node_path = self._csv_files(tmp_path)

        sets = load_reference_sets(
            ke_wp_df=pd.DataFrame([{"KE_ID": "KE:177", "WP_ID": "WP4324"}]),
            wp_gene_path=wp_gene_path,
            node_path=node_path,
            wp_gene_map={"WP4324": {"MT-CO1"}},
        )

        assert sets["KE:177"] == {"MT-CO1"}, "stale TP53 should not survive"

    def test_unresolved_pathways_surface_as_a_warning(self):
        """The run must say so — that is the whole complaint in #79."""
        resolution = [{
            'resource': 'WikiPathways', 'status': 'loaded', 'source': 'api',
            'ke_count': 90, 'confidence_applied': True,
            'unresolved_pathways': ['WP1234', 'WP5477'],
            'unresolved_ke_pathways': {'KE:1115': ['WP5477']},
            'error': None,
        }]

        warnings = app.resource_resolution_warnings(resolution)
        joined = " ".join(warnings)
        assert "WP5477" in joined and "WP1234" in joined
        # Issue #81: the affected Key Events are named, and the sentence no
        # longer says they are reported as unmapped — they are not any more.
        assert "KE:1115" in joined
        assert "no gene set" not in joined
        assert "curated" in joined

    def test_no_warning_when_everything_resolves(self):
        resolution = [{
            'resource': 'WikiPathways', 'status': 'loaded', 'source': 'api',
            'ke_count': 90, 'confidence_applied': True,
            'unresolved_pathways': [], 'error': None,
        }]
        assert app.resource_resolution_warnings(resolution) == []


class TestUnresolvedWarningIsScopedToTheAOP:
    """Issue #108: the warning must describe the AOP being analysed.

    The resolution is built over the whole reference universe, so before this
    scoping every run claimed its coverage was understated — naming pathways
    belonging to Key Events in other AOPs entirely. Asserting a gap that does
    not exist is worse than the silent gap #79/#81 fixed: it travels into the
    caveats an author writes around the numbers on the same page.
    """

    def _resolution(self):
        return [{
            'resource': 'WikiPathways', 'status': 'loaded', 'source': 'api',
            'ke_count': 90, 'confidence_applied': True,
            'unresolved_pathways': ['WP1234', 'WP3980', 'WP4010'],
            'unresolved_ke_pathways': {
                'KE:840': ['WP1234'],
                'KE:344': ['WP3980'],
                'KE:1115': ['WP4010'],
            },
            'error': None,
        }]

    def test_pathways_from_other_aops_are_dropped(self):
        scoped = app.scope_resolution_to_aop(self._resolution(), {'KE:1115'})

        assert scoped[0]['unresolved_pathways'] == ['WP4010']
        assert scoped[0]['unresolved_ke_pathways'] == {'KE:1115': ['WP4010']}

    def test_warning_is_suppressed_when_this_aop_lost_nothing(self):
        """The AOP:472 case from the report: none of the three is its own."""
        scoped = app.scope_resolution_to_aop(
            self._resolution(), {'KE:1', 'KE:2', 'KE:3'}
        )

        assert scoped[0]['unresolved_pathways'] == []
        assert app.resource_resolution_warnings(scoped) == []

    def test_surviving_warning_names_only_this_aops_key_events(self):
        scoped = app.scope_resolution_to_aop(self._resolution(), {'KE:1115'})

        joined = " ".join(app.resource_resolution_warnings(scoped))
        assert "WP4010" in joined and "KE:1115" in joined
        for elsewhere in ("WP1234", "WP3980", "KE:840", "KE:344"):
            assert elsewhere not in joined

    def test_unknown_key_events_leave_the_resolution_alone(self):
        """Scoping we cannot perform must not silently empty the accounting.

        An over-broad warning is recoverable; one suppressed because the AOP
        topology failed to load hides a real gap with no trace.
        """
        original = self._resolution()

        assert app.scope_resolution_to_aop(original, None) == original
        assert app.scope_resolution_to_aop(original, set()) == original

    def test_the_source_resolution_is_not_mutated(self):
        """The resolution can come from a cache shared with the next caller."""
        original = self._resolution()

        app.scope_resolution_to_aop(original, {'KE:1115'})

        assert original[0]['unresolved_pathways'] == ['WP1234', 'WP3980', 'WP4010']
        assert set(original[0]['unresolved_ke_pathways']) == {
            'KE:840', 'KE:344', 'KE:1115'
        }

class TestCacheFormatUpgrade:
    """The #79 cache entry gained a third element; warm caches must survive.

    Production runs a persistent disk cache on the Gluster mount, so the first
    request after a deploy reads entries written by the previous build. The
    compatibility branch in _load_wikipathways_reference_sets is what stops
    those turning into a 500, and nothing exercised it.
    """

    def _cache(self, monkeypatch, store):
        fake = MagicMock()
        fake.get.side_effect = lambda key, expire_time=False: (
            (store.get(key), None) if expire_time else store.get(key))
        fake.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake)

    def test_pre_79_two_tuple_entry_still_loads(self, monkeypatch):
        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        self._cache(monkeypatch, {key: ({"KE:1": {"A"}}, "api")})

        sets, source, unresolved, unresolved_ke = app._load_wikipathways_reference_sets("all")

        assert sets == {"KE:1": {"A"}}
        assert source == "cache(api)"
        assert unresolved == [], "an old entry knows of no unresolved pathways"

    def test_new_three_tuple_entry_round_trips(self, monkeypatch):
        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        self._cache(monkeypatch, {key: ({"KE:1": {"A"}}, "api", ["WP5477"])})

        sets, source, unresolved, unresolved_ke = app._load_wikipathways_reference_sets("all")

        assert unresolved == ["WP5477"]
        assert source == "cache(api)"


class TestCacheAgeIsReported:
    """Issue #106: "cached" without an age cannot support a reproducibility claim."""

    TTL = 3600

    @staticmethod
    def _cache_with_expiry(monkeypatch, store, expiry_by_key):
        """A cache double that answers expire_time, as diskcache does."""
        fake = MagicMock()
        fake.get.side_effect = lambda key, expire_time=False: (
            (store.get(key), expiry_by_key.get(key)) if expire_time else store.get(key))
        fake.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
        monkeypatch.setattr(app, "_reference_cache", fake)

    @pytest.fixture(autouse=True)
    def _clear_recorded_times(self):
        """The fill-time record is module state; keep tests independent."""
        app._CACHE_FILL_TIMES.clear()
        yield
        app._CACHE_FILL_TIMES.clear()

    def test_cache_hit_records_the_fill_time(self, monkeypatch):
        """Fill time is the entry's expiry less the TTL it was written with."""
        import datetime as dt

        filled = dt.datetime(2026, 7, 22, 9, 14, tzinfo=dt.timezone.utc)
        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        self._cache_with_expiry(
            monkeypatch,
            {key: ({"KE:1": {"A"}}, "api")},
            {key: filled.timestamp() + app.Config.CACHE_TTL},
        )

        app._load_wikipathways_reference_sets("all")

        assert app._cache_fill_time("WikiPathways", "all") == "2026-07-22 09:14 UTC"

    def test_provenance_line_states_the_age(self, monkeypatch):
        """The line a methods section quotes must carry the vintage."""
        import datetime as dt

        filled = dt.datetime(2026, 7, 22, 9, 14, tzinfo=dt.timezone.utc)
        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        self._cache_with_expiry(
            monkeypatch,
            {key: ({"KE:1": {"A"}}, "api")},
            {key: filled.timestamp() + app.Config.CACHE_TTL},
        )

        _, _, resolution = app.load_cached_reference_sets(["WikiPathways"])

        assert resolution[0]["cached_at"] == "2026-07-22 09:14 UTC"
        assert app.describe_resource_resolution(resolution) == (
            "WikiPathways (Builder API, cached 2026-07-22 09:14 UTC)"
        )

    def test_a_live_load_claims_no_cache_age(self):
        """Nothing was cached, so there is nothing to date."""
        wp = {"KE:1": {"A"}}
        with patch.object(app, "_load_wikipathways_reference_sets",
                          return_value=(wp, "api", [], {})):
            _, _, resolution = app.load_cached_reference_sets(["WikiPathways"])

        assert resolution[0]["cached_at"] is None
        assert app.describe_resource_resolution(resolution) == (
            "WikiPathways (Builder API, live)"
        )

    def test_unknown_age_keeps_the_bare_wording(self, monkeypatch):
        """An entry written without an expiry dates to nothing — don't invent one."""
        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, "all")
        self._cache_with_expiry(monkeypatch, {key: ({"KE:1": {"A"}}, "csv")}, {})

        _, _, resolution = app.load_cached_reference_sets(["WikiPathways"])

        assert resolution[0]["cached_at"] is None
        assert app.describe_resource_resolution(resolution) == (
            "WikiPathways (bundled reference files, cached)"
        )

    def test_gmt_resources_report_their_own_age(self, monkeypatch):
        """GO BP and Reactome are cached separately, so each dates separately."""
        import datetime as dt

        filled = dt.datetime(2026, 7, 20, 6, 0, tzinfo=dt.timezone.utc)
        key = app._confidence_cache_key(app._GMT_RESOURCE_CACHE_KEYS["GO_BP"], "all")
        self._cache_with_expiry(
            monkeypatch,
            {key: ({"KE:1097": {"TP53"}}, "api")},
            {key: filled.timestamp() + app.Config.CACHE_TTL},
        )

        _, _, resolution = app.load_cached_reference_sets(["GO_BP"])

        assert resolution[0]["cached_at"] == "2026-07-20 06:00 UTC"

    def test_a_skipped_resource_carries_the_key_but_no_age(self):
        """Every entry has the field, so consumers need no special-casing."""
        with patch.object(app, "_load_gmt_resource_reference_sets",
                          side_effect=RuntimeError("HTTP 502")):
            _, _, resolution = app.load_cached_reference_sets(["Reactome"])

        assert resolution[0]["status"] == "skipped"
        assert resolution[0]["cached_at"] is None
