"""Tests for the Builder API client service (issue #55 GMT resource loading)."""
import pytest
from unittest.mock import MagicMock, patch

from services.api_service import (
    parse_gmt_reference_sets,
    fetch_gmt_reference_sets,
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
