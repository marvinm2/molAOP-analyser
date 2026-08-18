"""Tests for the Builder API client service (#55 GMT resources, #60 confidence filter)."""
import pytest
from unittest.mock import MagicMock, patch

from services.api_service import (
    parse_gmt_reference_sets,
    fetch_gmt_reference_sets,
    fetch_ke_wp_records,
    fetch_reference_sets_from_api,
    filter_records_by_confidence,
    confidence_rank,
    GMT_RESOURCE_PATHS,
)


# A small GMT fixture mirroring the live Builder exports: one gene set per line,
# tab-separated, with the descriptor column encoding "KE<id>_<name>_<pathway_id>".
SAMPLE_GMT = (
    "KE177_Increase_Mitochondrial_dysfunction_WP5241\tMito beta-oxidation\tECHS1\tACSL3\tACADL\n"
    "KE177_Increase_Mitochondrial_dysfunction_WP78\tTCA cycle\tDLD\tECHS1\n"  # ECHS1 repeats
    "KE55_Increase_Cell_injury_R-HSA-5357769\tCaspase activation\tCASP3\tCASP8\n"
    "\n"  # blank line -> skipped
    "KE999_No_genes_WP1\tEmpty set\n"  # descriptor + title only, no genes -> skipped
    "malformed_line_without_ke_token\tTitle\tGENEX\n"  # no KE token -> skipped
)


class TestParseGmtReferenceSets:
    """Tests for GMT text parsing into KE->gene reference sets."""

    def test_unions_genes_across_pathways_per_ke(self):
        result = parse_gmt_reference_sets(SAMPLE_GMT)
        # KE:177 unions WP5241 + WP78 genes, de-duplicating ECHS1.
        assert result["KE:177"] == {"ECHS1", "ACSL3", "ACADL", "DLD"}

    def test_separate_kes_kept_apart(self):
        result = parse_gmt_reference_sets(SAMPLE_GMT)
        assert result["KE:55"] == {"CASP3", "CASP8"}

    def test_ke_id_normalised_with_colon(self):
        result = parse_gmt_reference_sets(SAMPLE_GMT)
        assert all(k.startswith("KE:") for k in result)

    def test_blank_geneless_and_malformed_lines_skipped(self):
        result = parse_gmt_reference_sets(SAMPLE_GMT)
        # KE:999 had no genes; malformed line had no KE token.
        assert "KE:999" not in result
        assert all("GENEX" not in genes for genes in result.values())

    def test_genes_uppercased(self):
        result = parse_gmt_reference_sets("KE1_x_WP1\tTitle\tabc\tDef\n")
        assert result["KE:1"] == {"ABC", "DEF"}

    def test_empty_input_returns_empty_dict(self):
        assert parse_gmt_reference_sets("") == {}


# The provenance block the Builder actually emits, copied verbatim from
# GET /exports/gmt/ke-wp on 2026-08-18. Its own last line tells consumers that
# '#' lines are not gene sets.
BUILDER_HEADER = (
    "# molAOP Builder GMT export\n"
    "# resource: KE-WP\n"
    "# export-revision: 0e6d8234525b9882\n"
    "# source-fingerprint: 447:de059ae56c077366426bc3114b6dd55e|\n"
    "# confidence: all tiers\n"
    "# generated: 2026-08-18T06:40:56+00:00\n"
    "# Lines beginning with # are provenance, not gene sets.\n"
)


class TestGmtProvenanceHeader:
    """The Builder's '#' provenance header must never be read as a gene set.

    Before 2026-08-18 neither parser had an explicit '#' guard; they survived the
    header only because comment lines contain fewer than three tab-separated
    fields and were dropped by the geneless-line check. That is luck, not intent:
    one tab in a future header line — a tab-aligned key/value, say — would have
    been parsed as a descriptor plus genes and silently polluted the reference
    data, which is the failure mode that would be hardest to notice downstream.
    """

    def test_header_does_not_leak_into_reference_sets(self):
        with_header = BUILDER_HEADER + SAMPLE_GMT

        assert parse_gmt_reference_sets(with_header) == parse_gmt_reference_sets(SAMPLE_GMT)

    def test_header_alone_yields_nothing(self):
        assert parse_gmt_reference_sets(BUILDER_HEADER) == {}

    def test_tabbed_header_line_is_still_skipped(self):
        """The case the old field-count guard would have missed."""
        hostile = "# generated:\t2026-08-18\tby\tthe\tbuilder\n" + SAMPLE_GMT

        result = parse_gmt_reference_sets(hostile)

        assert result == parse_gmt_reference_sets(SAMPLE_GMT)

    def test_pathway_gene_map_skips_the_header_too(self):
        from services.api_service import parse_gmt_pathway_gene_map

        body = "KE177_Mito_WP5241\tFatty acid oxidation\tACSL3\tACADL\n"
        # A '#' line mentioning a WP id must not be mined for one.
        hostile = "# source-fingerprint:\tWP5241\t447\n" + BUILDER_HEADER + body

        assert parse_gmt_pathway_gene_map(hostile) == parse_gmt_pathway_gene_map(body)


class TestFetchGmtReferenceSets:
    """Tests for fetching + parsing a resource GMT export from the Builder."""

    def _config(self, url="https://builder.example.com"):
        cfg = MagicMock()
        cfg.BUILDER_API_URL = url
        cfg.BUILDER_API_TIMEOUT = 10
        return cfg

    def test_fetches_and_parses(self):
        mock_resp = MagicMock()
        mock_resp.text = SAMPLE_GMT
        mock_resp.raise_for_status = MagicMock()
        mock_session = MagicMock()
        mock_session.get.return_value = mock_resp

        import services.api_service as api_svc
        with patch.object(api_svc, "_make_api_session", return_value=mock_session):
            result = fetch_gmt_reference_sets(self._config(), "GO_BP")

        # Correct endpoint hit for the resource.
        called_url = mock_session.get.call_args[0][0]
        assert called_url.endswith(GMT_RESOURCE_PATHS["GO_BP"])
        assert result["KE:177"] == {"ECHS1", "ACSL3", "ACADL", "DLD"}

    def test_unknown_resource_raises(self):
        with pytest.raises(ValueError):
            fetch_gmt_reference_sets(self._config(), "NotAResource")

    def test_unset_builder_url_raises(self):
        with pytest.raises(ValueError):
            fetch_gmt_reference_sets(self._config(url=""), "Reactome")

    def _captured_params(self, min_confidence):
        """Run a fetch and return the query params actually sent."""
        mock_resp = MagicMock()
        mock_resp.text = SAMPLE_GMT
        mock_resp.raise_for_status = MagicMock()
        mock_session = MagicMock()
        mock_session.get.return_value = mock_resp

        import services.api_service as api_svc
        with patch.object(api_svc, "_make_api_session", return_value=mock_session):
            fetch_gmt_reference_sets(
                self._config(), "GO_BP", min_confidence=min_confidence
            )
        return mock_session.get.call_args.kwargs.get("params")

    def test_threshold_is_forwarded_to_the_builder(self):
        """Issue #71: the whole point — the GMT export does the filtering."""
        assert self._captured_params("high") == {"min_confidence": "high"}
        assert self._captured_params("medium") == {"min_confidence": "medium"}

    def test_all_is_sent_as_an_omitted_parameter_not_forwarded(self):
        """`all` is our sentinel and is NOT in the Builder's whitelist.

        Forwarding it verbatim returns HTTP 400, which would skip the resource
        for the whole run — a silently narrower analysis, not a visible error.
        """
        assert self._captured_params("all") == {}

    def test_unrecognised_threshold_degrades_to_unfiltered(self):
        """An unfiltered gene set beats a 400 that drops the resource."""
        assert self._captured_params("bogus") == {}
        assert self._captured_params(None) == {}


class TestGmtMinConfidenceParam:
    """Issue #71: the vocabulary mismatch that makes this translation necessary.

    Ours is ("all", "medium", "high"); the Builder's whitelist is
    {high, medium, low}. The overlap is partial in both directions, so the
    mapping is stated explicitly rather than assumed.
    """

    def test_our_vocabulary_maps_cleanly(self):
        from services.api_service import gmt_min_confidence_param

        assert gmt_min_confidence_param("high") == {"min_confidence": "high"}
        assert gmt_min_confidence_param("medium") == {"min_confidence": "medium"}
        assert gmt_min_confidence_param("all") == {}

    def test_never_emits_a_value_outside_the_builder_whitelist(self):
        """Anything we send must be one the Builder accepts, or the run loses
        the resource. Enumerated rather than spot-checked."""
        from helpers import VALID_MIN_CONFIDENCE
        from services.api_service import gmt_min_confidence_param

        builder_whitelist = {"high", "medium", "low"}
        for value in VALID_MIN_CONFIDENCE:
            params = gmt_min_confidence_param(value)
            if params:
                assert params["min_confidence"] in builder_whitelist


class TestFilterRecordsByConfidence:
    """Issue #60: minimum KE-mapping confidence filtering of raw records."""

    # Mirrors the live Builder payload: confidence_level values are lowercase.
    RECORDS = [
        {"ke_id": "KE 1", "pathway_id": "WP1", "confidence_level": "high"},
        {"ke_id": "KE 1", "pathway_id": "WP2", "confidence_level": "low"},
        {"ke_id": "KE 2", "pathway_id": "WP3", "confidence_level": "medium"},
        {"ke_id": "KE 3", "pathway_id": "WP4", "confidence_level": "low"},
    ]

    def test_all_is_a_no_op(self):
        assert filter_records_by_confidence(self.RECORDS, "all") == self.RECORDS

    def test_default_is_all(self):
        assert filter_records_by_confidence(self.RECORDS) == self.RECORDS

    def test_medium_keeps_medium_and_high(self):
        kept = filter_records_by_confidence(self.RECORDS, "medium")
        assert [r["pathway_id"] for r in kept] == ["WP1", "WP3"]

    def test_high_keeps_only_high(self):
        kept = filter_records_by_confidence(self.RECORDS, "high")
        assert [r["pathway_id"] for r in kept] == ["WP1"]

    def test_lowercase_and_titlecase_handled_identically(self):
        """The live API returns 'high'; the API doc shows 'High'. Both must rank alike."""
        for records in (
            [{"pathway_id": "WP1", "confidence_level": "high"}],
            [{"pathway_id": "WP1", "confidence_level": "High"}],
            [{"pathway_id": "WP1", "confidence_level": "  HIGH  "}],
        ):
            assert filter_records_by_confidence(records, "high") == records
            assert filter_records_by_confidence(records, "medium") == records

    @pytest.mark.parametrize("value", [None, "", "   ", "unknown", "Very High"])
    def test_missing_or_unrecognised_confidence_is_kept(self, value):
        """Graceful degradation: an unknown confidence never removes a mapping."""
        records = [{"pathway_id": "WP1", "confidence_level": value}]
        assert filter_records_by_confidence(records, "high") == records
        assert filter_records_by_confidence(records, "medium") == records

    def test_records_without_the_field_are_kept(self):
        """CSV fallback rows (KE_ID/WP_ID only) must survive any threshold."""
        records = [{"KE_ID": "KE:1", "WP_ID": "WP1"}]
        assert filter_records_by_confidence(records, "high") == records

    def test_unknown_threshold_falls_back_to_no_filtering(self):
        assert filter_records_by_confidence(self.RECORDS, "bogus") == self.RECORDS

    def test_input_is_not_mutated(self):
        original = [dict(r) for r in self.RECORDS]
        filter_records_by_confidence(self.RECORDS, "high")
        assert self.RECORDS == original

    def test_confidence_rank_ordering(self):
        assert confidence_rank("high") > confidence_rank("medium") > confidence_rank("low")
        assert confidence_rank(None) is None
        assert confidence_rank("bogus") is None


class TestFetchKeWpRecordsConfidence:
    """fetch_ke_wp_records applies the threshold to the fetched records (#60)."""

    RECORDS = [
        {"ke_id": "KE 1", "pathway_id": "WP1", "confidence_level": "high"},
        {"ke_id": "KE 2", "pathway_id": "WP2", "confidence_level": "low"},
    ]

    def _config(self):
        cfg = MagicMock()
        cfg.BUILDER_API_URL = "https://builder.example"
        cfg.BUILDER_API_TIMEOUT = 5
        return cfg

    def test_threshold_applied_to_fetched_records(self):
        with patch("services.api_service.fetch_all_ke_wp_mappings", return_value=self.RECORDS):
            kept = fetch_ke_wp_records(self._config(), min_confidence="high")
        assert [r["pathway_id"] for r in kept] == ["WP1"]

    def test_default_returns_every_record(self):
        with patch("services.api_service.fetch_all_ke_wp_mappings", return_value=self.RECORDS):
            kept = fetch_ke_wp_records(self._config())
        assert len(kept) == 2


class TestFetchReferenceSetsConfidence:
    """fetch_reference_sets_from_api pre-filters mappings before the gene merge (#60)."""

    RECORDS = [
        {"ke_id": "KE 1", "pathway_id": "WP1", "confidence_level": "high"},
        {"ke_id": "KE 1", "pathway_id": "WP2", "confidence_level": "low"},
        {"ke_id": "KE 3", "pathway_id": "WP4", "confidence_level": "low"},
    ]

    def _config(self):
        cfg = MagicMock()
        cfg.BUILDER_API_URL = "https://builder.example"
        cfg.BUILDER_API_TIMEOUT = 5
        return cfg

    def _captured_df(self, min_confidence):
        with patch("services.api_service.fetch_all_ke_wp_mappings", return_value=self.RECORDS), \
             patch("services.api_service.load_reference_sets", return_value={}) as loader:
            fetch_reference_sets_from_api(self._config(), min_confidence=min_confidence)
        return loader.call_args.kwargs["ke_wp_df"]

    def test_all_passes_every_mapping_through(self):
        df = self._captured_df("all")
        assert sorted(df["WP_ID"]) == ["WP1", "WP2", "WP4"]

    def test_high_keeps_only_high_confidence_mappings(self):
        df = self._captured_df("high")
        # KE:1 keeps WP1 only; KE:3 (all Low) disappears entirely.
        assert list(df["WP_ID"]) == ["WP1"]
        assert list(df["KE_ID"]) == ["KE:1"]


# ---------------------------------------------------------------------------
# Issue #79: pathway gene membership must not come only from the bundled CSV
# ---------------------------------------------------------------------------

class TestPathwayGeneMap:
    """parse_gmt_pathway_gene_map / fetch_wp_pathway_gene_map."""

    GMT = (
        "KE1115_Increase_Reactive_oxygen_species_WP5477\tOxidative stress\tSOD1\tCAT\n"
        "KE1392_Increase_Oxidative_Stress_WP5477\tOxidative stress\tSOD1\tGPX1\n"
        "KE177_Increase_Mitochondrial_dysfunction_WP5507\tMito\tNDUFA1\n"
        "KE999_No_genes_WP1\tEmpty\n"
        "not a descriptor\tignored\tXYZ\n"
    )

    def test_keys_on_pathway_not_ke(self):
        """The map is pathway-centric — that is what the KE-WP merge needs."""
        from services.api_service import parse_gmt_pathway_gene_map

        result = parse_gmt_pathway_gene_map(self.GMT)
        assert set(result) == {"WP5477", "WP5507"}

    def test_unions_genes_across_kes_sharing_a_pathway(self):
        """WP5477 is mapped by two KEs; membership is the union of both rows."""
        from services.api_service import parse_gmt_pathway_gene_map

        assert parse_gmt_pathway_gene_map(self.GMT)["WP5477"] == {"SOD1", "CAT", "GPX1"}

    def test_rows_without_genes_or_pathway_are_skipped(self):
        """A descriptor with no genes must not create an empty gene set."""
        from services.api_service import parse_gmt_pathway_gene_map

        result = parse_gmt_pathway_gene_map(self.GMT)
        assert "WP1" not in result
        assert all(genes for genes in result.values())

    def test_fetch_returns_empty_when_builder_unset(self):
        """No Builder configured falls back to the CSV rather than raising."""
        from services.api_service import fetch_wp_pathway_gene_map

        cfg = MagicMock()
        cfg.BUILDER_API_URL = ""
        assert fetch_wp_pathway_gene_map(cfg) == {}

    def test_fetch_swallows_builder_failure(self):
        """A Builder outage must degrade to the CSV, not fail the analysis."""
        from services.api_service import fetch_wp_pathway_gene_map

        cfg = MagicMock()
        cfg.BUILDER_API_URL = "https://b.example.com"
        cfg.BUILDER_API_TIMEOUT = 5
        session = MagicMock()
        session.get.side_effect = RuntimeError("connection refused")
        with patch("services.api_service._make_api_session", return_value=session):
            assert fetch_wp_pathway_gene_map(cfg) == {}
