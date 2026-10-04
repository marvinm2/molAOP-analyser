"""Issue #123: how the selected gene-set sources combine per Key Event.

Until #123 the per-KE gene sets of every selected resource were always unioned.
The user can now ask for the intersection (a gene counts only if every selected
source has it) or for "at least N sources". Union stays the default, and a
union run must produce exactly what it produced before the option existed —
the golden tests below rebuild the pre-#123 merge inline and compare against it.
"""

from unittest.mock import patch

import pandas as pd
import pytest

import app
from helpers import (
    DEFAULT_SOURCE_COMBINATION,
    ReferenceSets,
    combine_source_sets,
    no_shared_genes_kes_for,
    normalise_source_combination,
    parse_source_combination,
    source_combination_label,
)
from services.enrichment_service import (
    EXCLUDED_NO_MAPPING,
    EXCLUDED_NO_SHARED_GENES,
    format_ke_summary,
    get_ke_summary,
    run_enrichment,
)
from services.network_service import build_cytoscape_network, ke_accounting_from_network

# Three overlapping sources, shaped like a real AOP: KE:1 overlaps partially,
# KE:2 is mapped by WikiPathways only, KE:3 is mapped everywhere but shares no
# gene, KE:4 is shared by exactly two sources, KE:5 resolved to nothing anywhere.
WP = {
    "KE:1": {f"G{i}" for i in range(0, 30)},
    "KE:2": {f"G{i}" for i in range(40, 50)},
    "KE:3": {"G60", "G61"},
    "KE:4": {f"G{i}" for i in range(70, 80)},
    "KE:5": set(),
}
GO = {
    "KE:1": {f"G{i}" for i in range(10, 40)},
    "KE:3": {"G62", "G63"},
    "KE:4": {f"G{i}" for i in range(70, 80)},
}
REACTOME = {
    "KE:1": {f"G{i}" for i in range(5, 25)},
    "KE:3": {"G64"},
    "KE:4": {"G99"},
}
SOURCES = {"WikiPathways": WP, "GO_BP": GO, "Reactome": REACTOME}


def _legacy_union(per_source):
    """The merge exactly as it stood on main before #123 (app.py)."""
    merged = {}
    for sets in per_source:
        for ke_id, genes in sets.items():
            merged.setdefault(ke_id, set()).update(genes)
    return merged


def _patched_loaders():
    """Patch both per-resource loaders to serve the three fixtures above."""

    def _gmt(resource, min_confidence="all"):
        return SOURCES[resource], "api"

    return (
        patch.object(
            app, "_load_wikipathways_reference_sets", return_value=(WP, "api", [], {})
        ),
        patch.object(app, "_load_gmt_resource_reference_sets", side_effect=_gmt),
    )


def _load(resources, **kwargs):
    wp_patch, gmt_patch = _patched_loaders()
    with wp_patch, gmt_patch:
        return app.load_cached_reference_sets(resources, **kwargs)


def _expression_df(n=120, n_sig=25):
    genes = [f"G{i}" for i in range(n)]
    return pd.DataFrame(
        {
            "ID": genes,
            "log2FC": [2.0 - i * 0.03 for i in range(n)],
            "pval": [1e-5 if i < n_sig else 0.4 for i in range(n)],
            "significant": [i < n_sig for i in range(n)],
        }
    )


ALL = ["WikiPathways", "GO_BP", "Reactome"]
KE_LIST = {"KE:1", "KE:2", "KE:3", "KE:4", "KE:5"}


class TestParseSourceCombination:
    """Form values -> one canonical spec, validated against the selection."""

    def test_default_is_union(self):
        assert DEFAULT_SOURCE_COMBINATION == "union"
        assert parse_source_combination(None, None, 3) == "union"
        assert parse_source_combination("", "", 3) == "union"

    def test_intersection(self):
        assert parse_source_combination("intersection", None, 3) == "intersection"

    def test_case_and_whitespace(self):
        assert parse_source_combination("  Intersection ", None, 2) == "intersection"

    def test_at_least_n(self):
        assert parse_source_combination("at_least", "2", 3) == "at_least_2"

    def test_at_least_endpoints_collapse_to_union_and_intersection(self):
        """N=1 *is* union and N=all *is* intersection — one spelling each."""
        assert parse_source_combination("at_least", "1", 3) == "union"
        assert parse_source_combination("at_least", "3", 3) == "intersection"

    def test_single_source_mode_is_irrelevant(self):
        """With one source every mode is the same set, so record it as union."""
        assert parse_source_combination("intersection", None, 1) == "union"

    def test_min_n_above_source_count_is_rejected(self):
        with pytest.raises(ValueError, match="at most 2"):
            parse_source_combination("at_least", "3", 2)

    def test_min_n_above_source_count_is_rejected_for_a_single_source(self):
        with pytest.raises(ValueError):
            parse_source_combination("at_least", "2", 1)

    @pytest.mark.parametrize("bad_n", ["0", "-1", "two", "", None, "1.5"])
    def test_bad_min_n_is_rejected(self, bad_n):
        with pytest.raises(ValueError):
            parse_source_combination("at_least", bad_n, 3)

    @pytest.mark.parametrize("junk", ["xor", "at_least_2", "'; DROP TABLE --"])
    def test_unknown_mode_is_rejected(self, junk):
        with pytest.raises(ValueError):
            parse_source_combination(junk, None, 3)

    def test_stored_values_normalise_with_null_as_union(self):
        """Rows written before #123 carry NULL and were unioned."""
        assert normalise_source_combination(None) == "union"
        assert normalise_source_combination("") == "union"
        assert normalise_source_combination("intersection") == "intersection"
        assert normalise_source_combination("at_least_2") == "at_least_2"
        assert normalise_source_combination("bogus") == "union"
        assert normalise_source_combination("at_least_1") == "union"

    def test_labels(self):
        assert "union" in source_combination_label("union").lower()
        assert "intersection" in source_combination_label("intersection").lower()
        assert "at least 2" in source_combination_label("at_least_2").lower()


class TestCombineSourceSets:
    """The combination itself, independent of the loaders."""

    def test_union_is_the_legacy_merge(self):
        combined, no_shared = combine_source_sets([WP, GO, REACTOME], "union")
        assert combined == _legacy_union([WP, GO, REACTOME])
        assert no_shared == set()

    def test_intersection_keeps_only_genes_in_every_source(self):
        combined, no_shared = combine_source_sets([WP, GO, REACTOME], "intersection")
        assert combined["KE:1"] == {f"G{i}" for i in range(10, 25)}
        # Mapped by one source only, or by all three without a shared gene.
        assert "KE:2" not in combined and "KE:3" not in combined
        assert "KE:4" not in combined
        assert no_shared == {"KE:2", "KE:3", "KE:4"}

    def test_at_least_two(self):
        combined, no_shared = combine_source_sets([WP, GO, REACTOME], "at_least_2")
        assert combined["KE:1"] == {f"G{i}" for i in range(5, 30)}
        assert combined["KE:4"] == {f"G{i}" for i in range(70, 80)}
        assert no_shared == {"KE:2", "KE:3"}

    def test_a_ke_with_no_genes_anywhere_is_not_a_sharing_failure(self):
        """KE:5 resolved to nothing; that stays the #81 'unresolved' case."""
        combined, no_shared = combine_source_sets([WP, GO, REACTOME], "intersection")
        assert combined.get("KE:5") == set()
        assert "KE:5" not in no_shared

    def test_inputs_are_not_mutated(self):
        before = {k: set(v) for k, v in WP.items()}
        combine_source_sets([WP, GO], "intersection")
        combine_source_sets([WP, GO], "union")
        assert WP == before

    def test_single_source_every_mode_is_identical(self):
        for spec in ("union", "intersection"):
            combined, no_shared = combine_source_sets([WP], spec)
            assert combined == _legacy_union([WP])
            assert no_shared == set()


class TestLoaderHonoursTheMode:
    """load_cached_reference_sets applies the mode and exposes the exclusions."""

    def test_union_golden_equivalence_with_main(self):
        """CRITICAL: union output is byte-for-byte what main produced."""
        default_sets, default_source, default_res = _load(ALL)
        union_sets, union_source, union_res = _load(ALL, source_combination="union")
        legacy = _legacy_union([WP, GO, REACTOME])
        assert dict(default_sets) == legacy
        assert dict(union_sets) == legacy
        assert union_source == default_source
        assert union_res == default_res
        assert no_shared_genes_kes_for(union_sets) == set()

    def test_intersection_sets_and_exclusions(self):
        sets, _, _ = _load(ALL, source_combination="intersection")
        assert sets["KE:1"] == {f"G{i}" for i in range(10, 25)}
        assert no_shared_genes_kes_for(sets) == {"KE:2", "KE:3", "KE:4"}

    def test_mode_switch_never_serves_another_modes_sets(self):
        """Per-resource caches are mode-independent; the combination is not
        cached, so alternating modes can never hand back the wrong sets."""
        a, _, _ = _load(ALL, source_combination="intersection")
        b, _, _ = _load(ALL, source_combination="union")
        c, _, _ = _load(ALL, source_combination="intersection")
        assert dict(b) == _legacy_union([WP, GO, REACTOME])
        assert dict(a) == dict(c) != dict(b)

    def test_invalid_spec_falls_back_to_union(self):
        sets, _, _ = _load(ALL, source_combination="bogus")
        assert dict(sets) == _legacy_union([WP, GO, REACTOME])

    def test_combination_counts_loaded_sources_only(self):
        """A skipped resource cannot veto every gene under intersection."""

        def _gmt(resource, min_confidence="all"):
            if resource == "Reactome":
                raise RuntimeError("down")
            return SOURCES[resource], "api"

        with (
            patch.object(
                app,
                "_load_wikipathways_reference_sets",
                return_value=(WP, "api", [], {}),
            ),
            patch.object(app, "_load_gmt_resource_reference_sets", side_effect=_gmt),
        ):
            sets, _, resolution = app.load_cached_reference_sets(
                ALL, source_combination="intersection"
            )
        assert sets["KE:1"] == {f"G{i}" for i in range(10, 30)}
        warnings = app.resource_resolution_warnings(
            resolution, source_combination="intersection"
        )
        assert any("Reactome" in w and "intersection" in w.lower() for w in warnings)


class TestEnrichmentGolden:
    """Union-mode enrichment results are identical to main, ORA and GSEA."""

    @pytest.mark.parametrize("method", ["ora", "gsea"])
    def test_union_results_identical(self, method):
        df = _expression_df()
        titles = {k: k for k in KE_LIST}
        legacy_sets = _legacy_union([WP, GO, REACTOME])
        union_sets, _, _ = _load(ALL, source_combination="union")
        kwargs = {"permutation_num": 50} if method == "gsea" else {}

        expected = run_enrichment(method, df, legacy_sets, KE_LIST, titles, **kwargs)
        actual = run_enrichment(
            method,
            df,
            union_sets,
            KE_LIST,
            titles,
            no_shared_genes_kes=no_shared_genes_kes_for(union_sets),
            **kwargs,
        )
        pd.testing.assert_frame_equal(actual, expected)
        assert get_ke_summary(actual) == get_ke_summary(expected)


class TestExclusionReporting:
    """A KE emptied by the combination is reported, not silently dropped."""

    @pytest.mark.parametrize("method", ["ora", "gsea"])
    def test_no_shared_genes_is_its_own_reason(self, method):
        df = _expression_df()
        sets, _, _ = _load(ALL, source_combination="intersection")
        kwargs = {"permutation_num": 50} if method == "gsea" else {}
        result = run_enrichment(
            method,
            df,
            sets,
            KE_LIST | {"KE:NOSET"},
            {},
            no_shared_genes_kes=no_shared_genes_kes_for(sets),
            **kwargs,
        )
        summary = get_ke_summary(result)
        for ke in ("KE:2", "KE:3", "KE:4"):
            assert summary["excluded_reasons"][ke] == EXCLUDED_NO_SHARED_GENES
        assert summary["excluded_reasons"]["KE:NOSET"] == EXCLUDED_NO_MAPPING
        assert summary["excluded_no_shared_genes"] == 3
        assert summary["excluded_no_mapping"] == 1
        sentence = format_ke_summary(summary)
        assert "3 excluded (no genes shared across the selected sources)" in sentence

    def test_accounting_survives_the_stored_network(self):
        """Batch report and shared links rebuild the counts from node payloads."""
        df = _expression_df()
        sets, _, _ = _load(ALL, source_combination="intersection")
        result = run_enrichment(
            "ora",
            df,
            sets,
            KE_LIST,
            {},
            no_shared_genes_kes=no_shared_genes_kes_for(sets),
        )
        summary = get_ke_summary(result)
        network = build_cytoscape_network(
            KE_LIST,
            pd.DataFrame(columns=["Source_KE", "Target_KE"]),
            result,
            {},
            {},
            reference_sets=sets,
            excluded_kes=summary["excluded_reasons"],
        )
        rebuilt = ke_accounting_from_network(network)
        assert rebuilt["excluded_no_shared_genes"] == 3
        node = next(n for n in network["nodes"] if n["data"]["id"] == "KE:3")
        assert "no-genes" in node["classes"]

    def test_reference_sets_carry_the_exclusions(self):
        sets = ReferenceSets({"KE:1": {"A"}}, no_shared_genes_kes={"KE:9"})
        assert no_shared_genes_kes_for(sets) == {"KE:9"}
        assert no_shared_genes_kes_for({"KE:1": {"A"}}) == set()


class TestProvenance:
    """The choice appears wherever the run's sources are described."""

    _RESOLUTION = [
        {"resource": r, "status": "loaded", "source": "api", "ke_count": 3} for r in ALL
    ]

    def test_provenance_line_names_the_combination(self):
        text = app.describe_resource_resolution(
            self._RESOLUTION, source_combination="intersection"
        )
        assert text.startswith("WikiPathways (Builder API")
        assert "intersection" in text.lower()

    def test_provenance_line_names_union_for_multi_source_runs(self):
        text = app.describe_resource_resolution(
            self._RESOLUTION, source_combination="union"
        )
        assert "union" in text.lower()

    def test_single_source_provenance_is_unchanged(self):
        one = self._RESOLUTION[:1]
        assert app.describe_resource_resolution(
            one, source_combination="union"
        ) == app.describe_resource_resolution(one)

    def test_provenance_without_a_combination_is_unchanged(self):
        assert "combined" not in app.describe_resource_resolution(self._RESOLUTION)

    def test_report_lists_the_combination(self, sample_report_data):
        from services.report_service import report_generator

        sample_report_data.selected_resources = "WikiPathways, GO_BP"
        sample_report_data.source_combination = "intersection"
        html = report_generator.generate_html_report(sample_report_data)
        assert "Source Combination" in html
        assert source_combination_label("intersection") in html


class TestPersistence:
    """Stored on experiments and batches; NULL reads back as union."""

    def test_experiment_round_trip(self, temp_database):
        exp_id = temp_database.save_experiment_metadata(
            metadata={"dataset_id": "d"},
            analysis_params={"aop_id": "AOP:1", "source_combination": "intersection"},
        )
        assert (
            temp_database.get_experiment(exp_id)["source_combination"] == "intersection"
        )

    def test_legacy_experiment_reads_as_union(self, temp_database):
        exp_id = temp_database.save_experiment_metadata(
            metadata={"dataset_id": "d"},
            analysis_params={"aop_id": "AOP:1"},
        )
        assert temp_database.get_experiment(exp_id)["source_combination"] == "union"

    def test_migration_is_idempotent(self, tmp_path):
        from sqlalchemy import create_engine, text

        from database import _ensure_source_combination_column

        engine = create_engine(f"sqlite:///{tmp_path / 'legacy.db'}")
        with engine.connect() as conn:
            conn.execute(text("CREATE TABLE experiments (id INTEGER PRIMARY KEY)"))
            conn.execute(text("CREATE TABLE batches (id INTEGER PRIMARY KEY)"))
            conn.commit()
        _ensure_source_combination_column(engine)
        _ensure_source_combination_column(engine)
        with engine.connect() as conn:
            for table in ("experiments", "batches"):
                cols = [r[1] for r in conn.execute(text(f"PRAGMA table_info({table})"))]
                assert cols.count("source_combination") == 1

    def test_batch_effective_value(self):
        from database import BatchRecord

        assert BatchRecord().effective_source_combination() == "union"
        assert (
            BatchRecord(source_combination="at_least_2").effective_source_combination()
            == "at_least_2"
        )


class TestSingleRoute:
    """/analyze parses, validates, forwards and records the choice."""

    @staticmethod
    def _form(**overrides):
        form = {
            "filename": "test.csv",
            "id_column": "Gene_Symbol",
            "fc_column": "log2FoldChange",
            "pval_column": "padj",
            "aop_selection": "AOP:1",
            "logfc_threshold": "1.0",
            "resources": ["WikiPathways", "GO_BP", "Reactome"],
        }
        form.update(overrides)
        return form

    @staticmethod
    def _post(client, form):
        processed_df = pd.DataFrame(
            {
                "ID": ["BRCA1", "TP53"],
                "log2FC": [1.5, -0.8],
                "pval": [0.001, 0.05],
                "significant": [True, False],
            }
        )
        enrichment_df = pd.DataFrame(
            {
                "Title": ["Test KE"],
                "p_value": [0.01],
                "FDR": [0.05],
                "num_overlap": [1],
                "pct_sig_in_KE": [50.0],
                "total_KE_genes_in_dataset": [2],
                "odds_ratio": [3.5],
                "overlap": ["BRCA1"],
                "KE": ["KE:115"],
                "sig_in_KE": [1],
                "sig_not_KE": [0],
                "non_sig_in_KE": [1],
                "non_sig_not_KE": [0],
            }
        )
        edges_df = pd.DataFrame(columns=["Source_KE", "Target_KE", "KER_ID", "AOP_ID"])
        resolution = [
            {"resource": r, "status": "loaded", "source": "api", "ke_count": 1}
            for r in ALL
        ]
        sets = ReferenceSets({"KE:115": {"BRCA1"}}, no_shared_genes_kes={"KE:9"})
        with (
            patch("app.load_and_validate_data", return_value=processed_df),
            patch(
                "app.process_gene_expression",
                return_value=(processed_df, {"total_genes": 2}),
            ),
            patch(
                "app.load_aop_data",
                return_value=(
                    {"KE:115"},
                    edges_df,
                    {"KE:115": "KE"},
                    {"KE:115": "Test KE"},
                ),
            ),
            patch("app.run_enrichment", return_value=enrichment_df) as enrich,
            patch(
                "app.build_cytoscape_network", return_value={"nodes": [], "edges": []}
            ),
            patch("app.build_ke_gene_mapping", return_value={}),
            patch("app.guess_id_type", return_value="HGNC"),
            patch("app.validate_file_path", return_value=True),
            patch(
                "app.load_cached_reference_sets", return_value=(sets, "api", resolution)
            ) as loader,
        ):
            response = client.post("/analyze", data=form)
        return response, loader, enrich

    def test_omitted_defaults_to_union(self, authenticated_client):
        response, loader, _ = self._post(authenticated_client, self._form())
        assert response.status_code == 200
        assert loader.call_args.kwargs["source_combination"] == "union"

    def test_intersection_forwarded_and_exclusions_reach_the_backend(
        self, authenticated_client
    ):
        response, loader, enrich = self._post(
            authenticated_client, self._form(source_combination="intersection")
        )
        assert response.status_code == 200
        assert loader.call_args.kwargs["source_combination"] == "intersection"
        assert enrich.call_args.kwargs["no_shared_genes_kes"] == {"KE:9"}
        html = response.data.decode()
        assert source_combination_label("intersection") in html

    def test_at_least_n_forwarded(self, authenticated_client):
        response, loader, _ = self._post(
            authenticated_client,
            self._form(source_combination="at_least", source_min_n="2"),
        )
        assert response.status_code == 200
        assert loader.call_args.kwargs["source_combination"] == "at_least_2"

    def test_min_n_above_source_count_is_a_400(self, authenticated_client):
        response, loader, _ = self._post(
            authenticated_client,
            self._form(source_combination="at_least", source_min_n="4"),
        )
        assert response.status_code == 400
        assert b"source" in response.data.lower()
        loader.assert_not_called()

    def test_junk_mode_is_a_400(self, authenticated_client):
        response, loader, _ = self._post(
            authenticated_client, self._form(source_combination="xor")
        )
        assert response.status_code == 400
        loader.assert_not_called()

    def test_single_source_records_union(self, authenticated_client):
        response, loader, _ = self._post(
            authenticated_client,
            self._form(resources=["WikiPathways"], source_combination="intersection"),
        )
        assert response.status_code == 200
        assert loader.call_args.kwargs["source_combination"] == "union"

    def test_choice_is_posted_to_the_report(self, authenticated_client):
        response, _, _ = self._post(
            authenticated_client, self._form(source_combination="intersection")
        )
        assert (
            'name="source_combination" value="intersection"' in response.data.decode()
        )


class TestBatchRoute:
    """Batch parity: same parsing, stored on the batch, applied to the load."""

    @staticmethod
    def _form(**overrides):
        form = {
            "batch_uuid": "test-batch-uuid",
            "aop_selection": "AOP:1",
            "id_col": "GENE_SYMBOL",
            "fc_col": "logFC",
            "pval_col": "adj.P.Val",
            "logfc_threshold": "0.0",
            "filenames[]": "a.csv",
            "condition_labels[]": "A",
            "resources": ["WikiPathways", "GO_BP"],
        }
        form.update(overrides)
        return form

    @staticmethod
    def _post(client, form):
        with (
            patch("os.path.isdir", return_value=True),
            patch("os.path.isfile", return_value=True),
            patch("app.validate_batch_columns", return_value=(True, "")),
            patch("app.harmonise_backgrounds", return_value=({"BRCA1"}, {})),
            patch("app._persist_and_launch_batch", return_value=1) as launcher,
        ):
            response = client.post("/batch/analyze", data=form)
        return response, launcher

    def test_forwarded(self, flask_client):
        response, launcher = self._post(
            flask_client, self._form(source_combination="intersection")
        )
        assert response.status_code == 200
        assert launcher.call_args.kwargs["source_combination"] == "intersection"

    def test_default_union(self, flask_client):
        response, launcher = self._post(flask_client, self._form())
        assert response.status_code == 200
        assert launcher.call_args.kwargs["source_combination"] == "union"

    def test_min_n_above_source_count_is_a_400(self, flask_client):
        response, launcher = self._post(
            flask_client, self._form(source_combination="at_least", source_min_n="3")
        )
        assert response.status_code == 400
        launcher.assert_not_called()

    def test_stored_and_applied(self, temp_database, monkeypatch):
        import app as app_module
        from database import BatchRecord

        monkeypatch.setattr(app_module, "db_manager", temp_database)
        monkeypatch.setattr(app_module, "run_batch", lambda *a, **k: None)
        with patch.object(
            app_module, "load_cached_reference_sets", return_value=({}, "mock", [])
        ) as loader:
            app_module._persist_and_launch_batch(
                batch_uuid="combo-uuid",
                filenames=["a.csv"],
                condition_labels=["A"],
                doses=[""],
                timepoints=[""],
                id_col="Gene",
                fc_col="logFC",
                pval_col="pval",
                aop_id="AOP:1",
                logfc_threshold=0.0,
                pval_threshold=0.05,
                resources=["WikiPathways", "GO_BP"],
                harmonised_genes={"BRCA1"},
                batch_name="b",
                owner="",
                description="",
                source_combination="intersection",
            )
        assert loader.call_args.kwargs["source_combination"] == "intersection"
        session = temp_database.get_session()
        try:
            batch = session.query(BatchRecord).filter_by(uuid="combo-uuid").one()
            assert batch.source_combination == "intersection"
        finally:
            session.close()

    def test_run_batch_reports_the_exclusion(
        self, temp_database, tmp_path, monkeypatch
    ):
        """The batch thread gets only the gene sets; the exclusion rides on them."""
        import json

        from database import ConditionRecord
        from services.batch_service import run_batch
        from tests.test_ke_exclusion_wiring import _seed_single_condition_batch

        upload_root = tmp_path / "uploads"
        upload_root.mkdir()
        monkeypatch.setattr("config.Config.UPLOAD_FOLDER", str(upload_root))
        monkeypatch.setattr(
            "services.data_service.load_aop_data",
            lambda aop_id: (
                {"KE:OK", "KE:SPLIT"},
                pd.DataFrame(columns=["Source_KE", "Target_KE", "KER_ID", "AOP_ID"]),
                {"KE:OK": "MIE", "KE:SPLIT": "AO"},
                {},
            ),
        )
        batch_id = _seed_single_condition_batch(temp_database, str(upload_root))
        sets = ReferenceSets(
            {"KE:OK": {f"G{i}" for i in range(1, 11)}},
            no_shared_genes_kes={"KE:SPLIT"},
        )
        run_batch(batch_id, temp_database.db_url, sets)

        session = temp_database.get_session()
        try:
            cond = session.query(ConditionRecord).filter_by(batch_id=batch_id).first()
            assert cond.status == "complete"
            node = next(
                n["data"]
                for n in json.loads(cond.network_json)["nodes"]
                if n["data"]["id"] == "KE:SPLIT"
            )
            assert node["excluded_reason"] == EXCLUDED_NO_SHARED_GENES
        finally:
            session.close()

    def test_control_rendered_in_the_batch_wizard(self, flask_client):
        response = flask_client.get("/")
        assert b"batch-source-combination" in response.data


class TestSingleFormControl:
    """Minimal control, beside the resource checkboxes, only for >1 source."""

    @staticmethod
    def _render(**context):
        from flask import render_template

        base = dict(
            preview=[{"Gene_Symbol": "BRCA1"}],
            columns=["Gene_Symbol", "log2FoldChange", "padj"],
            filename="test.csv",
            volcano_data=[{"ID": "BRCA1", "log2FC": 1.5, "pval": 0.001}],
            selected_columns={
                "id": "Gene_Symbol",
                "fc": "log2FoldChange",
                "pval": "padj",
                "pval_adj": None,
            },
            column_suggestions=None,
            logfc_threshold=1.0,
            pval_cutoff=0.05,
            pval_y=[],
            columns_confirmed=True,
            case_study_aops={},
            cisplatin_demos=[],
            parse_filename=lambda name: {},
            recommended_aops=None,
            metadata={},
            pval_threshold=0.05,
            method="ora",
            selected_resources=["WikiPathways"],
        )
        base.update(context)
        with app.app.test_request_context("/"):
            return render_template("_single_analysis.html", **base)

    def test_hidden_with_one_source(self):
        html = self._render()
        assert 'name="source_combination"' in html
        assert 'id="source-combination-group" hidden' in html

    def test_shown_with_several_sources(self):
        html = self._render(
            selected_resources=["WikiPathways", "GO_BP"],
            source_combination="intersection",
        )
        assert 'id="source-combination-group" hidden' not in html
        assert '<option value="intersection" selected>' in html
