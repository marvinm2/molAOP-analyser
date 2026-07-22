"""Issue #81: a Key Event that is mapped is never reported as unmapped.

The exclusion accounting used to collapse three different findings into
"no gene set mapped":

* the Key Event has no curated mapping at all — a curation gap, fixed in the
  Builder;
* the Key Event is mapped but the pathway could not be resolved to genes — a
  stale-reference-data bug in the analyser (issue #79 plumbed the unresolvable
  pathway IDs through for exactly this reason);
* the Key Event's gene set is known but none of its genes were measured in the
  uploaded dataset — a coverage fact about the experiment.

Only the first is a curation gap, so these tests pin each of them to its own
counter and check that the mapped-but-unresolvable case names the pathways.
"""
import json

import pandas as pd
import pytest

from services.enrichment_service import (
    EXCLUDED_NO_MAPPING,
    EXCLUDED_TOO_FEW_GENES,
    EXCLUDED_UNRESOLVED_MAPPING,
    format_ke_summary,
    get_ke_summary,
    normalise_unresolved_ke_pathways,
    run_enrichment_analysis,
)
from services.network_service import build_cytoscape_network, ke_accounting_from_network


def _make_gene_df(n_genes=50, n_sig=10):
    """Build a minimal processed expression frame with `n_sig` significant genes."""
    genes = [f"GENE{i}" for i in range(n_genes)]
    return pd.DataFrame({
        'ID': genes,
        'significant': [i < n_sig for i in range(n_genes)],
        'log2FC': [1.0 if i < n_sig else 0.0 for i in range(n_genes)],
        'pval': [0.001 if i < n_sig else 0.5 for i in range(n_genes)],
    })


def _run(reference_sets, ke_list, **kwargs):
    """Run ORA over a fixed background and return the accounting summary."""
    df = _make_gene_df()
    result = run_enrichment_analysis(
        df, reference_sets, ke_list, {k: k for k in ke_list}, **kwargs
    )
    return get_ke_summary(result)


class TestExclusionReasonsAreDistinct:
    """The three exclusions are counted separately (issue #81)."""

    def test_mapped_but_unresolvable_is_not_a_curation_gap(self):
        """A KE whose pathway could not be resolved is never called unmapped."""
        reference_sets = {'KE:OK': {f"GENE{i}" for i in range(10)}}
        ke_list = {'KE:OK', 'KE:1115', 'KE:NOSET'}

        summary = _run(
            reference_sets, ke_list,
            unresolved_ke_pathways={'KE:1115': ['WP5477']},
        )

        assert summary['excluded_reasons']['KE:1115'] == EXCLUDED_UNRESOLVED_MAPPING
        assert summary['excluded_reasons']['KE:NOSET'] == EXCLUDED_NO_MAPPING
        assert summary['excluded_unresolved_mapping'] == 1
        assert summary['excluded_no_mapping'] == 1
        # The offending pathway is named, so the reader can check it against
        # the Builder rather than being sent to curate a mapping that exists.
        assert summary['unresolved_pathways'] == ['WP5477']

    def test_empty_reference_gene_set_counts_as_unresolvable(self):
        """A mapping that survives into the reference sets with no genes is mapped."""
        reference_sets = {
            'KE:OK': {f"GENE{i}" for i in range(10)},
            'KE:EMPTY': set(),
        }
        summary = _run(reference_sets, {'KE:OK', 'KE:EMPTY'})

        assert summary['excluded_reasons']['KE:EMPTY'] == EXCLUDED_UNRESOLVED_MAPPING
        assert summary['excluded_no_mapping'] == 0

    def test_unmeasured_gene_set_is_a_coverage_finding(self):
        """A known gene set with nothing measured is 'too few measured genes'."""
        reference_sets = {
            'KE:OK': {f"GENE{i}" for i in range(10)},
            'KE:MISS': {'ABSENT1', 'ABSENT2', 'ABSENT3'},
        }
        summary = _run(reference_sets, {'KE:OK', 'KE:MISS'})

        assert summary['excluded_reasons']['KE:MISS'] == EXCLUDED_TOO_FEW_GENES
        assert summary['excluded_no_mapping'] == 0
        assert summary['excluded_unresolved_mapping'] == 0

    def test_accounting_stays_exhaustive(self):
        """Every KE of the AOP is still accounted for exactly once."""
        reference_sets = {
            'KE:OK': {f"GENE{i}" for i in range(10)},
            'KE:SMALL': {'GENE0', 'GENE1'},
            'KE:MISS': {'ABSENT1'},
        }
        ke_list = {'KE:OK', 'KE:SMALL', 'KE:MISS', 'KE:NOSET', 'KE:1115'}
        summary = _run(
            reference_sets, ke_list,
            unresolved_ke_pathways={'KE:1115': ('WP5477',)},
        )

        assert (
            summary['tested']
            + summary['excluded_too_few_genes']
            + summary['excluded_no_mapping']
            + summary['excluded_unresolved_mapping']
            + summary['excluded_error']
            == summary['total_kes'] == 5
        )

    def test_unresolved_hint_for_a_tested_ke_is_ignored(self):
        """A KE that resolved to genes is tested, whatever the hint claims."""
        reference_sets = {'KE:OK': {f"GENE{i}" for i in range(10)}}
        summary = _run(
            reference_sets, {'KE:OK'},
            unresolved_ke_pathways={'KE:OK': ['WP1']},
        )

        assert summary['tested'] == 1
        assert summary['excluded_reasons'] == {}
        assert summary['unresolved_pathways'] == []


class TestNormaliseUnresolvedKePathways:
    """The hint accepts what the pipeline actually produces."""

    def test_none_and_empty_entries_are_dropped(self):
        assert normalise_unresolved_ke_pathways(None) == {}
        assert normalise_unresolved_ke_pathways({'KE:1': [], 'KE:2': None}) == {}

    def test_single_string_is_accepted(self):
        assert normalise_unresolved_ke_pathways({'KE:1': 'wp5477'}) == {'KE:1': ['WP5477']}

    def test_ids_are_deduplicated_and_sorted(self):
        assert normalise_unresolved_ke_pathways(
            {'KE:1': {'WP99', 'wp10', ' WP99 '}}
        ) == {'KE:1': ['WP10', 'WP99']}


class TestFormatKeSummaryWording:
    """One wording authority for the results page and both reports."""

    def test_unresolved_clause_names_the_pathways(self):
        text = format_ke_summary({
            'total_kes': 9, 'tested': 6, 'excluded_no_mapping': 1,
            'excluded_unresolved_mapping': 1, 'excluded_too_few_genes': 1,
            'excluded_error': 0, 'min_ke_genes': 5,
            'unresolved_pathways': ['WP5477'],
        })
        assert text == (
            '6 of 9 Key Events tested; 1 excluded (fewer than 5 measured genes), '
            '1 excluded (no gene set mapped), '
            '1 excluded (mapped, but no genes could be resolved: WP5477)'
        )

    def test_long_pathway_lists_are_summarised(self):
        text = format_ke_summary({
            'total_kes': 9, 'tested': 4, 'excluded_unresolved_mapping': 5,
            'min_ke_genes': 5,
            'unresolved_pathways': ['WP1', 'WP2', 'WP3', 'WP4', 'WP5'],
        })
        assert 'WP1, WP2, WP3 and 2 more' in text

    def test_unnamed_pathways_still_read_as_mapped(self):
        text = format_ke_summary({
            'total_kes': 9, 'tested': 8, 'excluded_unresolved_mapping': 1,
            'min_ke_genes': 5,
        })
        assert text.endswith('1 excluded (mapped, but no genes could be resolved)')

    def test_pre_issue_81_summary_is_unchanged(self):
        """Summaries stored before the split carry only the two old counters."""
        text = format_ke_summary({
            'total_kes': 24, 'tested': 18, 'excluded_no_mapping': 2,
            'excluded_too_few_genes': 4, 'excluded_error': 0, 'min_ke_genes': 5,
        })
        assert text == (
            '18 of 24 Key Events tested; 4 excluded (fewer than 5 measured genes), '
            '2 excluded (no gene set mapped)'
        )


def _network_for(excluded_kes, reference_sets, extra_kes=()):
    """Build a one-KE-per-reason network with a single tested KE."""
    ke_list = set(excluded_kes) | set(extra_kes) | {'KE:OK'}
    enrichment = pd.DataFrame({'KE': ['KE:OK'], 'p_value': [0.01], 'FDR': [0.01]})
    edges = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID'])
    return build_cytoscape_network(
        ke_list, edges, enrichment,
        {k: k for k in ke_list}, {k: 'intermediate' for k in ke_list},
        reference_sets=reference_sets,
        excluded_kes=excluded_kes,
    )


class TestNetworkStyling:
    """The network distinguishes the reasons it reports (issue #81)."""

    def test_unresolved_ke_is_not_styled_as_uncurated(self):
        network = _network_for(
            {'KE:1115': EXCLUDED_UNRESOLVED_MAPPING, 'KE:NOSET': EXCLUDED_NO_MAPPING},
            {'KE:OK': {'GENE0'}},
        )
        classes = {n['data']['id']: n['classes'] for n in network['nodes']}

        assert 'unresolved-mapping' in classes['KE:1115']
        assert 'no-genes' not in classes['KE:1115']
        assert 'no-genes' in classes['KE:NOSET']

    def test_unknown_reason_keeps_the_legacy_no_genes_styling(self):
        """Without an exclusion map, a KE with no gene set is still muted."""
        network = _network_for(
            {}, {'KE:OK': {'GENE0'}, 'KE:GAP': set()}, extra_kes={'KE:GAP'}
        )
        classes = {n['data']['id']: n['classes'] for n in network['nodes']}
        assert 'no-genes' in classes['KE:GAP']

    def test_accounting_recovered_from_a_stored_network(self):
        network = _network_for(
            {'KE:1115': EXCLUDED_UNRESOLVED_MAPPING, 'KE:NOSET': EXCLUDED_NO_MAPPING},
            {'KE:OK': {'GENE0'}},
        )
        summary = ke_accounting_from_network(json.dumps(network))

        assert summary['excluded_unresolved_mapping'] == 1
        assert summary['excluded_no_mapping'] == 1
        assert 'mapped, but no genes could be resolved' in format_ke_summary(summary)

    def test_pre_issue_81_network_reports_no_unresolved_kes(self):
        """A network stored before the split cannot claim a distinction it lacks."""
        stored = {
            'nodes': [
                {'data': {'id': 'KE:1', 'ke_type': 'MIE', 'excluded_reason': None}},
                {'data': {'id': 'KE:2', 'ke_type': 'AO', 'excluded_reason': 'no_mapping'}},
            ],
            'edges': [],
        }
        summary = ke_accounting_from_network(stored)

        assert summary['excluded_unresolved_mapping'] == 0
        assert summary['excluded_no_mapping'] == 1
        assert 'could not be resolved' not in format_ke_summary(summary)


class TestGseaParity:
    """ORA and GSEA report exclusions with one vocabulary."""

    def test_gsea_uses_the_same_reasons(self):
        gp = pytest.importorskip('gseapy')  # noqa: F841
        from services.gsea_service import run_gsea_analysis

        genes = [f"GENE{i:03d}" for i in range(60)]
        df = pd.DataFrame({
            'ID': genes,
            'log2FC': [(-1) ** i * (i / 10.0) for i in range(60)],
            'pval': [0.001 + i * 0.001 for i in range(60)],
        })
        reference_sets = {
            'KE:OK': set(genes[:20]),
            'KE:MISS': {'ABSENT1', 'ABSENT2'},
        }
        ke_list = {'KE:OK', 'KE:MISS', 'KE:NOSET', 'KE:1115'}

        result = run_gsea_analysis(
            df, reference_sets, ke_list, {k: k for k in ke_list},
            permutation_num=10,
            unresolved_ke_pathways={'KE:1115': ['WP5477']},
        )
        summary = get_ke_summary(result)

        assert summary['excluded_reasons']['KE:1115'] == EXCLUDED_UNRESOLVED_MAPPING
        assert summary['excluded_reasons']['KE:NOSET'] == EXCLUDED_NO_MAPPING
        assert summary['excluded_reasons']['KE:MISS'] == EXCLUDED_TOO_FEW_GENES
        assert summary['unresolved_pathways'] == ['WP5477']
