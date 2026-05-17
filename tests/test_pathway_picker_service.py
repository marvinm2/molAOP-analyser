"""Tests for the KE-WP picker data assembly service.

Unit tests for pathway_picker_service.build_wp_picker_data(), which assembles
the WP picker list from the intersection of the selected AOP's KEs, the available
KE-WP mappings, and the enrichment FDR ranking produced by run_enrichment().
These tests are intentionally RED at Wave 0 — services/pathway_picker_service.py
does not exist yet. They will turn GREEN once Plan 01 implements
build_wp_picker_data().
"""

import pandas as pd

from services.pathway_picker_service import build_wp_picker_data


def _make_picker_inputs(**overrides):
    """Helper: build the five-tuple (ke_list, ke_wp_records, enrichment_df, ke_title_map, wp_title_map).

    Returns default test data suitable for exercising build_wp_picker_data().
    Pass keyword overrides to replace individual components:
      ke_list, ke_wp_records, enrichment_df, ke_title_map, wp_title_map

    Default setup:
      - Two KEs in the AOP: KE:386 (rank 0, most enriched) and KE:123 (rank 1)
      - One API-style KE-WP record per KE, plus one record for an out-of-AOP KE
      - enrichment_df sorted ascending by FDR (as enrichment_service.py guarantees)
      - ke_title_map and wp_title_map populated for all IDs
    """
    ke_list = overrides.get('ke_list', {'KE:386', 'KE:123'})

    ke_wp_records = overrides.get('ke_wp_records', [
        {
            'ke_id': 'KE 386',           # API returns space-separated format
            'pathway_id': 'WP368',
            'pathway_title': 'Fatty acid beta-oxidation',
            'confidence_level': 'High',
        },
        {
            'ke_id': 'KE 123',
            'pathway_id': 'WP100',
            'pathway_title': 'Cell cycle',
            'confidence_level': 'Medium',
        },
        {
            'ke_id': 'KE 999',           # NOT in ke_list — should be filtered out
            'pathway_id': 'WP999',
            'pathway_title': 'Out-of-AOP pathway',
            'confidence_level': 'Low',
        },
    ])

    enrichment_df = overrides.get('enrichment_df', pd.DataFrame({
        'KE': ['KE:386', 'KE:123'],
        'FDR': [0.001, 0.05],
    }).sort_values('FDR').reset_index(drop=True))

    ke_title_map = overrides.get('ke_title_map', {
        'KE:386': 'Hepatic steatosis',
        'KE:123': 'Cell proliferation',
        'KE:999': 'Some other KE',
    })

    wp_title_map = overrides.get('wp_title_map', {
        'WP368': 'Fatty acid beta-oxidation',
        'WP100': 'Cell cycle',
    })

    return ke_list, ke_wp_records, enrichment_df, ke_title_map, wp_title_map


class TestBuildWpPickerData:
    """Unit tests for pathway_picker_service.build_wp_picker_data()."""

    def test_ordering(self):
        """Entries are sorted ascending by enrichment_rank (most-enriched KE first).

        KE:386 has FDR 0.001 (rank 0) and KE:123 has FDR 0.05 (rank 1).
        The WP entry for KE:386 must appear before the entry for KE:123.
        """
        ke_list, ke_wp_records, enrichment_df, ke_title_map, wp_title_map = _make_picker_inputs()

        result = build_wp_picker_data(ke_list, ke_wp_records, enrichment_df, ke_title_map, wp_title_map)

        assert isinstance(result, list)
        assert len(result) == 2
        # Most-enriched KE first
        assert result[0]['ke_id'] == 'KE:386'
        assert result[0]['enrichment_rank'] == 0
        assert result[1]['ke_id'] == 'KE:123'
        assert result[1]['enrichment_rank'] == 1
        # Verify key set on each entry
        expected_keys = {'ke_id', 'ke_title', 'wp_id', 'wp_title', 'confidence_level', 'enrichment_rank'}
        for entry in result:
            assert set(entry.keys()) == expected_keys, (
                f"Entry keys {set(entry.keys())} != expected {expected_keys}"
            )

    def test_ke_filtering(self):
        """KE-WP records whose ke_id is NOT in ke_list are excluded from the output.

        The default inputs include a record for KE:999 which is not in ke_list.
        It must not appear in the result.
        """
        ke_list, ke_wp_records, enrichment_df, ke_title_map, wp_title_map = _make_picker_inputs()

        result = build_wp_picker_data(ke_list, ke_wp_records, enrichment_df, ke_title_map, wp_title_map)

        result_ke_ids = [entry['ke_id'] for entry in result]
        assert 'KE:999' not in result_ke_ids, (
            "KE:999 is not in ke_list and must not appear in the picker output"
        )
        assert 'KE:386' in result_ke_ids
        assert 'KE:123' in result_ke_ids
        # Verify key set
        expected_keys = {'ke_id', 'ke_title', 'wp_id', 'wp_title', 'confidence_level', 'enrichment_rank'}
        for entry in result:
            assert set(entry.keys()) == expected_keys

    def test_confidence_normalisation(self):
        """A record with confidence_level 'High' yields an entry with confidence_level == 'high'.

        The API returns mixed-case values ('High'/'Medium'/'Low'); the picker must
        lowercase them for use as CSS class suffixes (.confidence-badge--high etc.).
        """
        ke_list = {'KE:386'}
        ke_wp_records = [
            {
                'ke_id': 'KE 386',
                'pathway_id': 'WP368',
                'pathway_title': 'Fatty acid beta-oxidation',
                'confidence_level': 'High',    # mixed-case from API
            },
        ]
        enrichment_df = pd.DataFrame({
            'KE': ['KE:386'],
            'FDR': [0.01],
        }).sort_values('FDR').reset_index(drop=True)
        ke_title_map = {'KE:386': 'Hepatic steatosis'}
        wp_title_map = {'WP368': 'Fatty acid beta-oxidation'}

        result = build_wp_picker_data(ke_list, ke_wp_records, enrichment_df, ke_title_map, wp_title_map)

        assert len(result) == 1
        assert result[0]['confidence_level'] == 'high', (
            f"Expected 'high' (lowercased), got {result[0]['confidence_level']!r}"
        )
        expected_keys = {'ke_id', 'ke_title', 'wp_id', 'wp_title', 'confidence_level', 'enrichment_rank'}
        assert set(result[0].keys()) == expected_keys

    def test_empty_mappings(self):
        """build_wp_picker_data returns [] (a list, not None) when no KE-WP record matches ke_list.

        If the selected AOP's KEs are not present in the ke_wp_records, the function
        must return an empty list so {% if wp_picker_data %} is falsy and the template
        section is hidden.
        """
        ke_list = {'KE:500', 'KE:501'}   # KEs in AOP
        ke_wp_records = [
            {
                'ke_id': 'KE 999',        # Not in ke_list
                'pathway_id': 'WP999',
                'pathway_title': 'Unrelated pathway',
                'confidence_level': 'Low',
            },
        ]
        enrichment_df = pd.DataFrame({
            'KE': ['KE:500', 'KE:501'],
            'FDR': [0.1, 0.2],
        }).sort_values('FDR').reset_index(drop=True)
        ke_title_map = {'KE:500': 'KE Five Hundred', 'KE:501': 'KE Five Hundred One'}
        wp_title_map = {}

        result = build_wp_picker_data(ke_list, ke_wp_records, enrichment_df, ke_title_map, wp_title_map)

        assert result == [], f"Expected [], got {result!r}"
        assert isinstance(result, list), "Must return a list, not None"

    def test_csv_fallback(self):
        """CSV-fallback records (KE_ID/WP_ID only) produce entries with WP ID as title and confidence_level == ''.

        The CSV fallback path uses data/KE-WP.csv which has only KE_ID and WP_ID columns.
        When pathway_title and confidence_level are absent:
        - wp_title should fall back to the WP ID itself (e.g. 'WP368')
        - confidence_level should be '' (empty string) — omit badge entirely in template
        """
        ke_list = {'KE:386'}
        # CSV-shape records — only KE_ID and WP_ID, no title or confidence
        ke_wp_records = [
            {
                'KE_ID': 'KE:386',        # CSV fallback uses KE_ID (uppercase key, colon format)
                'WP_ID': 'WP368',
            },
        ]
        enrichment_df = pd.DataFrame({
            'KE': ['KE:386'],
            'FDR': [0.02],
        }).sort_values('FDR').reset_index(drop=True)
        ke_title_map = {'KE:386': 'Hepatic steatosis'}
        wp_title_map = {}    # empty — no title lookup available in CSV path

        result = build_wp_picker_data(ke_list, ke_wp_records, enrichment_df, ke_title_map, wp_title_map)

        assert len(result) == 1
        entry = result[0]
        # WP title falls back to WP ID when wp_title_map is empty
        assert entry['wp_title'] == 'WP368', (
            f"Expected wp_title == 'WP368' (WP ID as fallback), got {entry['wp_title']!r}"
        )
        # confidence_level is empty string when source record has no confidence_level
        assert entry['confidence_level'] == '', (
            f"Expected confidence_level == '' for CSV fallback, got {entry['confidence_level']!r}"
        )
        expected_keys = {'ke_id', 'ke_title', 'wp_id', 'wp_title', 'confidence_level', 'enrichment_rank'}
        assert set(entry.keys()) == expected_keys
