"""Tests for the hub gene analysis service.

Unit tests for hub_service.compute_hub_genes(), which inverts the ke_gene_map
produced by enrichment_service.build_ke_gene_mapping() into a gene-ranked hub list.
These tests are intentionally RED at Wave 0 — services/hub_service.py does not exist yet.
"""

from services.hub_service import compute_hub_genes, HUB_THRESHOLD


def _make_ke_gene_map(ke_gene_pairs, significant_genes=None, pvalue_adj_map=None, log2fc_map=None):
    """Helper: build a minimal ke_gene_map dict for hub_service tests.

    Args:
        ke_gene_pairs: list of (ke_id, [gene_symbol, ...]) tuples.
        significant_genes: set of gene symbols that are significant (default: none are).
        pvalue_adj_map: dict mapping gene symbol -> pvalue_adj (default: all None).
        log2fc_map: dict mapping gene symbol -> log2FC (default: all 0.0).

    Returns:
        ke_gene_map shaped like build_ke_gene_mapping() output:
        {KE_ID: [{id, log2FC, significant, pvalue_raw, pvalue_adj}, ...]}
    """
    significant_genes = significant_genes or set()
    pvalue_adj_map = pvalue_adj_map or {}
    log2fc_map = log2fc_map or {}

    result = {}
    for ke_id, gene_symbols in ke_gene_pairs:
        result[ke_id] = []
        for sym in gene_symbols:
            result[ke_id].append({
                'id': sym,
                'log2FC': log2fc_map.get(sym, 0.0),
                'significant': sym in significant_genes,
                'pvalue_raw': None,
                'pvalue_adj': pvalue_adj_map.get(sym, None),
            })
    return result


class TestComputeHubGenes:
    """Unit tests for hub_service.compute_hub_genes()."""

    def test_ranking_order(self):
        """Genes are sorted descending by ke_count; ties sorted ascending by gene symbol."""
        # TP53 in 3 KEs, BRCA1 in 1 KE, MYC in 1 KE
        ke_gene_map = _make_ke_gene_map([
            ('KE:1', ['TP53', 'BRCA1']),
            ('KE:2', ['TP53', 'MYC']),
            ('KE:3', ['TP53']),
        ])
        ke_title_map = {'KE:1': 'KE 1 title', 'KE:2': 'KE 2 title', 'KE:3': 'KE 3 title'}

        result = compute_hub_genes(ke_gene_map, ke_title_map)

        assert isinstance(result, list)
        assert len(result) == 3
        # TP53 has the highest ke_count
        assert result[0]['gene'] == 'TP53'
        assert result[0]['ke_count'] == 3
        # Ties (BRCA1, MYC both ke_count=1) sorted ascending by gene symbol
        assert result[1]['gene'] == 'BRCA1'
        assert result[2]['gene'] == 'MYC'

    def test_output_has_all_columns(self):
        """Every result dict has exactly the eight required keys."""
        ke_gene_map = _make_ke_gene_map([
            ('KE:1', ['TP53', 'BRCA1']),
        ])
        ke_title_map = {'KE:1': 'KE 1 title'}

        result = compute_hub_genes(ke_gene_map, ke_title_map)

        expected_keys = {'gene', 'ke_count', 'ke_ids', 'ke_titles', 'log2fc', 'pvalue_adj', 'significant', 'is_hub'}
        for entry in result:
            assert set(entry.keys()) == expected_keys, (
                f"Entry keys {set(entry.keys())} != expected {expected_keys}"
            )

    def test_ke_membership_correct(self):
        """A gene in KE:1 and KE:2 has ke_ids == ['KE:1', 'KE:2'] and ke_titles matching the map."""
        ke_gene_map = _make_ke_gene_map([
            ('KE:1', ['TP53', 'EGFR']),
            ('KE:2', ['TP53']),
        ])
        ke_title_map = {'KE:1': 'KE One', 'KE:2': 'KE Two'}

        result = compute_hub_genes(ke_gene_map, ke_title_map)

        tp53 = next(r for r in result if r['gene'] == 'TP53')
        assert tp53['ke_count'] == 2
        assert set(tp53['ke_ids']) == {'KE:1', 'KE:2'}
        assert set(tp53['ke_titles']) == {'KE One', 'KE Two'}

        egfr = next(r for r in result if r['gene'] == 'EGFR')
        assert egfr['ke_count'] == 1
        assert egfr['ke_ids'] == ['KE:1']
        assert egfr['ke_titles'] == ['KE One']

    def test_expression_values_correct(self):
        """log2fc and pvalue_adj on a result entry equal the values supplied in ke_gene_map."""
        ke_gene_map = _make_ke_gene_map(
            [('KE:1', ['TP53'])],
            significant_genes={'TP53'},
            pvalue_adj_map={'TP53': 0.01},
            log2fc_map={'TP53': 2.1},
        )
        ke_title_map = {'KE:1': 'KE 1 title'}

        result = compute_hub_genes(ke_gene_map, ke_title_map)

        tp53 = next(r for r in result if r['gene'] == 'TP53')
        assert tp53['log2fc'] == 2.1
        assert tp53['pvalue_adj'] == 0.01

    def test_significant_flag(self):
        """The significant flag is propagated as-is from the input gene dict."""
        ke_gene_map = _make_ke_gene_map(
            [('KE:1', ['TP53', 'BRCA1'])],
            significant_genes={'TP53'},
        )
        ke_title_map = {'KE:1': 'KE 1 title'}

        result = compute_hub_genes(ke_gene_map, ke_title_map)

        by_gene = {r['gene']: r for r in result}
        assert by_gene['TP53']['significant'] is True
        assert by_gene['BRCA1']['significant'] is False

    def test_hub_threshold(self):
        """A gene with ke_count == HUB_THRESHOLD has is_hub == True."""
        assert HUB_THRESHOLD == 3, "SPEC mandates HUB_THRESHOLD = 3"

        ke_gene_map = _make_ke_gene_map([
            ('KE:1', ['TP53']),
            ('KE:2', ['TP53']),
            ('KE:3', ['TP53']),
        ])
        ke_title_map = {'KE:1': 'A', 'KE:2': 'B', 'KE:3': 'C'}

        result = compute_hub_genes(ke_gene_map, ke_title_map)

        tp53 = next(r for r in result if r['gene'] == 'TP53')
        assert tp53['ke_count'] == 3
        assert tp53['is_hub'] is True

    def test_below_hub_threshold(self):
        """A gene with ke_count == 2 (below HUB_THRESHOLD) has is_hub == False."""
        ke_gene_map = _make_ke_gene_map([
            ('KE:1', ['BRCA1']),
            ('KE:2', ['BRCA1']),
        ])
        ke_title_map = {'KE:1': 'A', 'KE:2': 'B'}

        result = compute_hub_genes(ke_gene_map, ke_title_map)

        brca1 = next(r for r in result if r['gene'] == 'BRCA1')
        assert brca1['ke_count'] == 2
        assert brca1['is_hub'] is False

    def test_method_independent_counts(self):
        """Calling compute_hub_genes twice with identical ke_gene_map yields identical ke_count values.

        Hub counts derive from ke_gene_map membership, never from enrichment method stats.
        """
        ke_gene_map = _make_ke_gene_map([
            ('KE:1', ['TP53', 'BRCA1']),
            ('KE:2', ['TP53', 'MYC']),
            ('KE:3', ['TP53']),
        ])
        ke_title_map = {'KE:1': 'A', 'KE:2': 'B', 'KE:3': 'C'}

        result1 = compute_hub_genes(ke_gene_map, ke_title_map)
        result2 = compute_hub_genes(ke_gene_map, ke_title_map)

        counts1 = {r['gene']: r['ke_count'] for r in result1}
        counts2 = {r['gene']: r['ke_count'] for r in result2}
        assert counts1 == counts2

    def test_none_pvalue_adj(self):
        """An input gene dict with pvalue_adj == None produces pvalue_adj == None in output; no exception raised."""
        ke_gene_map = _make_ke_gene_map(
            [('KE:1', ['TP53'])],
            pvalue_adj_map={},  # pvalue_adj will be None for TP53
        )
        ke_title_map = {'KE:1': 'KE 1 title'}

        result = compute_hub_genes(ke_gene_map, ke_title_map)

        tp53 = next(r for r in result if r['gene'] == 'TP53')
        assert tp53['pvalue_adj'] is None

    def test_empty_ke_gene_map(self):
        """compute_hub_genes({}, {}) returns [] and raises no exception."""
        result = compute_hub_genes({}, {})

        assert result == []
        assert isinstance(result, list)
