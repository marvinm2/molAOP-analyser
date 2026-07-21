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
