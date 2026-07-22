"""Issue #81, production wiring: does the running app actually build and forward
the KE -> unresolvable-pathway map?

The vocabulary for the three exclusion reasons is covered by
``test_ke_exclusion_reasons.py``. These tests cover the part that vocabulary is
useless without: the map has to be produced in ``helpers.load_reference_sets``
(the only place the association still exists, before the inner merge drops it),
survive every hand-off, and reach both enrichment call sites. Without that,
``unresolved_mapping`` is unreachable in production and KE:1115 on AOP:472 is
still reported as uncurated while carrying a live mapping to WP5477.
"""
import json
import os
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

import app
from tests.conftest import uploaded_file_exists
from helpers import (
    ReferenceSets,
    load_reference_sets,
    unresolved_ke_pathways_for,
)
from services.enrichment_service import (
    EXCLUDED_ERROR,
    EXCLUDED_NO_MAPPING,
    EXCLUDED_UNRESOLVED_MAPPING,
    format_ke_summary,
    get_ke_summary,
)
from services.network_service import (
    build_cytoscape_network,
    ke_accounting_from_network,
)

# The case in the issue: AOP:472's KE:1115 is mapped to WP5477, which no source
# in this repo can resolve to genes.
KE_UNRESOLVABLE = 'KE:1115'
WP_UNRESOLVABLE = 'WP5477'
KE_OK = 'KE:18'
WP_OK = 'WP3941'

RESOLVABLE_GENES = {f'G{i}' for i in range(1, 11)}


def _fake_wp_gene_csvs(tmp_path):
    """Minimal WP->gene and gene-annotation CSVs (neither one covers WP5477)."""
    wp_gene = tmp_path / 'edges.csv'
    node = tmp_path / 'nodes.csv'
    pd.DataFrame([
        {'WPID': WP_OK, 'gene_id': str(1000 + i), 'edge_id': i}
        for i in range(1, 11)
    ]).to_csv(wp_gene, index=False)
    pd.DataFrame([
        {'GeneID': 1000 + i, 'GeneName': f'G{i}'} for i in range(1, 11)
    ]).to_csv(node, index=False)
    return str(wp_gene), str(node)


class TestReferenceSetsCarryTheMap:
    """helpers.load_reference_sets is where the association still exists."""

    def test_unresolvable_pathway_is_attributed_to_its_key_event(self, tmp_path):
        wp_gene_path, node_path = _fake_wp_gene_csvs(tmp_path)
        ke_wp_df = pd.DataFrame([
            {'KE_ID': KE_OK, 'WP_ID': WP_OK},
            {'KE_ID': KE_UNRESOLVABLE, 'WP_ID': WP_UNRESOLVABLE},
        ])

        sets = load_reference_sets(
            ke_wp_df=ke_wp_df,
            wp_gene_path=wp_gene_path,
            node_path=node_path,
        )

        # The KE is absent from the sets entirely -- which is exactly why the
        # map is needed: nothing downstream could otherwise tell it from a KE
        # nobody has ever mapped.
        assert KE_UNRESOLVABLE not in sets
        assert sets[KE_OK] == RESOLVABLE_GENES
        assert unresolved_ke_pathways_for(sets) == {
            KE_UNRESOLVABLE: [WP_UNRESOLVABLE]
        }

    def test_sink_parameter_receives_the_same_map(self, tmp_path):
        wp_gene_path, node_path = _fake_wp_gene_csvs(tmp_path)
        sink = {}
        load_reference_sets(
            ke_wp_df=pd.DataFrame([
                {'KE_ID': KE_UNRESOLVABLE, 'WP_ID': WP_UNRESOLVABLE},
                {'KE_ID': KE_OK, 'WP_ID': WP_OK},
            ]),
            wp_gene_path=wp_gene_path,
            node_path=node_path,
            unresolved_ke_out=sink,
        )
        assert sink == {KE_UNRESOLVABLE: [WP_UNRESOLVABLE]}

    def test_a_key_event_keeps_its_resolvable_pathways(self, tmp_path):
        """A KE mapped to both a resolvable and an unresolvable pathway is
        still tested -- but the unresolvable one is recorded against it."""
        wp_gene_path, node_path = _fake_wp_gene_csvs(tmp_path)
        sets = load_reference_sets(
            ke_wp_df=pd.DataFrame([
                {'KE_ID': KE_OK, 'WP_ID': WP_OK},
                {'KE_ID': KE_OK, 'WP_ID': WP_UNRESOLVABLE},
            ]),
            wp_gene_path=wp_gene_path,
            node_path=node_path,
        )
        assert sets[KE_OK] == RESOLVABLE_GENES
        assert unresolved_ke_pathways_for(sets) == {KE_OK: [WP_UNRESOLVABLE]}

    def test_nothing_unresolvable_means_an_empty_map(self, tmp_path):
        wp_gene_path, node_path = _fake_wp_gene_csvs(tmp_path)
        sets = load_reference_sets(
            ke_wp_df=pd.DataFrame([{'KE_ID': KE_OK, 'WP_ID': WP_OK}]),
            wp_gene_path=wp_gene_path,
            node_path=node_path,
        )
        assert unresolved_ke_pathways_for(sets) == {}

    def test_a_plain_dict_is_tolerated(self):
        """Older cache entries and test fixtures are plain dicts."""
        assert unresolved_ke_pathways_for({KE_OK: {'A'}}) == {}


def _fake_reference_cache(monkeypatch):
    store = {}
    fake = MagicMock()
    fake.get.side_effect = lambda key: store.get(key)
    fake.set.side_effect = lambda key, value, expire=None: store.__setitem__(key, value)
    monkeypatch.setattr(app, '_reference_cache', fake)
    return store


class TestAppLoadersForwardTheMap:
    """app.py must carry the map from the Builder boundary to the caller."""

    def test_wikipathways_loader_returns_the_map(self, monkeypatch):
        _fake_reference_cache(monkeypatch)
        sets = ReferenceSets(
            {KE_OK: RESOLVABLE_GENES},
            unresolved_ke_pathways={KE_UNRESOLVABLE: [WP_UNRESOLVABLE]},
        )
        with patch.object(app, 'fetch_reference_sets_from_api',
                          return_value=(sets, [WP_UNRESOLVABLE])):
            _, source, unresolved, unresolved_ke = \
                app._load_wikipathways_reference_sets('all')

        assert source == 'api'
        assert unresolved == [WP_UNRESOLVABLE]
        assert unresolved_ke == {KE_UNRESOLVABLE: [WP_UNRESOLVABLE]}

    def test_the_map_survives_the_disk_cache(self, monkeypatch):
        store = _fake_reference_cache(monkeypatch)
        sets = ReferenceSets(
            {KE_OK: RESOLVABLE_GENES},
            unresolved_ke_pathways={KE_UNRESOLVABLE: [WP_UNRESOLVABLE]},
        )
        with patch.object(app, 'fetch_reference_sets_from_api',
                          return_value=(sets, [WP_UNRESOLVABLE])):
            app._load_wikipathways_reference_sets('all')

        assert len(store) == 1
        # Second call is served from the cache and must not lose the map.
        cached_sets, source, _, unresolved_ke = \
            app._load_wikipathways_reference_sets('all')
        assert source == 'cache(api)'
        assert unresolved_ke == {KE_UNRESOLVABLE: [WP_UNRESOLVABLE]}
        assert unresolved_ke_pathways_for(cached_sets) == {
            KE_UNRESOLVABLE: [WP_UNRESOLVABLE]
        }

    def test_pre_81_three_tuple_cache_entry_still_loads(self, monkeypatch):
        """A warm cache written before #81 must not blow up on deploy."""
        store = _fake_reference_cache(monkeypatch)
        key = app._confidence_cache_key(app.REFERENCE_CACHE_KEY, 'all')
        store[key] = ({KE_OK: {'A'}}, 'csv', [WP_UNRESOLVABLE])

        sets, source, unresolved, unresolved_ke = \
            app._load_wikipathways_reference_sets('all')
        assert sets == {KE_OK: {'A'}}
        assert source == 'cache(csv)'
        assert unresolved == [WP_UNRESOLVABLE]
        assert unresolved_ke == {}

    def test_merged_sets_carry_the_union_across_resources(self, monkeypatch):
        _fake_reference_cache(monkeypatch)
        wp = ReferenceSets({KE_OK: RESOLVABLE_GENES})
        with patch.object(
            app, '_load_wikipathways_reference_sets',
            return_value=(wp, 'api', [WP_UNRESOLVABLE],
                          {KE_UNRESOLVABLE: [WP_UNRESOLVABLE]}),
        ):
            merged, _, resolution = app.load_cached_reference_sets(['WikiPathways'])

        assert unresolved_ke_pathways_for(merged) == {
            KE_UNRESOLVABLE: [WP_UNRESOLVABLE]
        }
        assert resolution[0]['unresolved_ke_pathways'] == {
            KE_UNRESOLVABLE: [WP_UNRESOLVABLE]
        }


class TestAnalyseRouteClassifiesTheKeyEvent:
    """The end-to-end check #81 asks for, through the real /analyze handler."""

    @staticmethod
    def _dataset():
        rows = []
        for i in range(1, 21):
            sig = i <= 10
            rows.append({
                'ID': f'G{i}',
                'log2FC': 2.0 if sig else 0.1,
                'pval': 0.001 if sig else 0.5,
                'significant': sig,
            })
        return pd.DataFrame(rows)

    def test_ke1115_is_mapped_but_unresolvable_and_wp5477_is_named(
        self, authenticated_client, monkeypatch, tmp_path
    ):
        _fake_reference_cache(monkeypatch)
        wp_gene_path, node_path = _fake_wp_gene_csvs(tmp_path)

        # Only the network boundary is faked: the Builder says KE:1115 is
        # mapped to WP5477, and no source can resolve WP5477's membership.
        records = [
            {'ke_id': 'KE 18', 'pathway_id': WP_OK, 'confidence_level': 'High'},
            {'ke_id': 'KE 1115', 'pathway_id': WP_UNRESOLVABLE,
             'confidence_level': 'High'},
        ]
        real_load_reference_sets = app.load_reference_sets

        def _csv_backed(**kwargs):
            kwargs['wp_gene_path'] = wp_gene_path
            kwargs['node_path'] = node_path
            return real_load_reference_sets(**kwargs)

        df = self._dataset()
        edges = pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID'])
        ke_list = {KE_OK, KE_UNRESOLVABLE}

        with patch('services.api_service.fetch_all_ke_wp_mappings', return_value=records), \
             patch('services.api_service.fetch_wp_pathway_gene_map',
                   return_value={WP_OK: RESOLVABLE_GENES}), \
             patch('services.api_service.load_reference_sets', side_effect=_csv_backed), \
             patch('app.load_and_validate_data', return_value=df), \
             patch('app.process_gene_expression', return_value=(df, {'total_genes': 20})), \
             patch('app.load_aop_data', return_value=(
                 ke_list, edges,
                 {KE_OK: 'MIE', KE_UNRESOLVABLE: 'AO'},
                 {KE_OK: 'Event 18', KE_UNRESOLVABLE: 'Increase, Reactive oxygen species'},
             )), \
             patch('app.guess_id_type', return_value='HGNC'), \
             patch('app.cleanup_file'), \
             patch('os.path.exists', side_effect=uploaded_file_exists), \
             patch('app.validate_file_path', return_value=True):

            response = authenticated_client.post('/analyze', data={
                'filename': 'test.csv',
                'id_column': 'ID',
                'fc_column': 'log2FC',
                'pval_column': 'pval',
                'aop_selection': 'AOP:472',
                'logfc_threshold': '1.0',
            })

        assert response.status_code == 200
        body = response.get_data(as_text=True)
        # The whole point of #81: the sentence must say "mapped, but no genes
        # could be resolved", name WP5477, and not blame the curator.
        assert 'mapped, but no genes could be resolved' in body
        assert WP_UNRESOLVABLE in body
        assert 'excluded (no gene set mapped)' not in body


class TestBatchConditionsClassifyTheKeyEvent:
    """The batch thread is handed nothing but the gene sets (#81 finding 1)."""

    def test_run_batch_stores_the_reason_and_the_pathway_ids(
        self, temp_database, tmp_path, monkeypatch
    ):
        from database import ConditionRecord
        from services.batch_service import run_batch

        upload_root = tmp_path / 'uploads'
        upload_root.mkdir()
        monkeypatch.setattr('config.Config.UPLOAD_FOLDER', str(upload_root))

        def _fake_load_aop_data(aop_id):
            return (
                {KE_OK, KE_UNRESOLVABLE},
                pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID', 'AOP_ID']),
                {KE_OK: 'MIE', KE_UNRESOLVABLE: 'AO'},
                {KE_OK: 'Event 18', KE_UNRESOLVABLE: 'Increase, ROS'},
            )
        monkeypatch.setattr('services.data_service.load_aop_data', _fake_load_aop_data)

        batch_id = _seed_single_condition_batch(temp_database, str(upload_root))

        reference_sets = ReferenceSets(
            {KE_OK: RESOLVABLE_GENES},
            unresolved_ke_pathways={KE_UNRESOLVABLE: [WP_UNRESOLVABLE]},
        )
        run_batch(batch_id, temp_database.db_url, reference_sets)

        session = temp_database.get_session()
        try:
            cond = session.query(ConditionRecord).filter_by(batch_id=batch_id).first()
            assert cond.status == 'complete'
            network = json.loads(cond.network_json)
            node = next(
                n['data'] for n in network['nodes']
                if n['data']['id'] == KE_UNRESOLVABLE
            )
            assert node['excluded_reason'] == EXCLUDED_UNRESOLVED_MAPPING
            assert node['unresolved_pathways'] == [WP_UNRESOLVABLE]
        finally:
            session.close()


def _seed_single_condition_batch(db_manager, upload_root):
    import uuid as _uuid
    from datetime import datetime, timedelta, timezone
    from database import BatchRecord, ConditionRecord

    batch_uuid = str(_uuid.uuid4())
    batch_dir = os.path.join(upload_root, batch_uuid)
    os.makedirs(batch_dir, exist_ok=True)
    rows = [
        {'gene': f'G{i}', 'logFC': 2.0 if i <= 10 else 0.1,
         'pval': 0.001 if i <= 10 else 0.5}
        for i in range(1, 21)
    ]
    pd.DataFrame(rows).to_csv(
        os.path.join(batch_dir, 'cond0.tsv'), sep='\t', index=False
    )

    session = db_manager.get_session()
    try:
        batch = BatchRecord(
            uuid=batch_uuid, status='pending', aop_id='AOP:472',
            aop_label='Test AOP', logfc_threshold=1.0, pval_cutoff=0.05,
            selected_resources='WikiPathways', id_column='gene',
            fc_column='logFC', pval_column='pval',
            harmonised_background=json.dumps([f'G{i}' for i in range(1, 21)]),
            harmonised_gene_count=20, batch_name='Test Batch',
            expires_at=datetime.now(timezone.utc).replace(tzinfo=None)
            + timedelta(days=14),
        )
        session.add(batch)
        session.flush()
        session.add(ConditionRecord(
            batch_id=batch.id, position=0, filename='cond0.tsv',
            condition_label='Cond 0', status='pending',
        ))
        session.commit()
        return batch.id
    finally:
        session.close()


class TestStoredNetworksKeepThePathwayIds:
    """Finding 3: the batch report and shared links rebuild the sentence from
    the stored network, and the IDs are the actionable part."""

    @staticmethod
    def _network():
        results = pd.DataFrame({'KE': [KE_OK], 'p_value': [0.01], 'FDR': [0.02]})
        return build_cytoscape_network(
            {KE_OK, KE_UNRESOLVABLE},
            pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID']),
            results,
            {KE_OK: 'Event 18', KE_UNRESOLVABLE: 'Increase, ROS'},
            {KE_OK: 'MIE', KE_UNRESOLVABLE: 'AO'},
            reference_sets={KE_OK: RESOLVABLE_GENES},
            excluded_kes={KE_UNRESOLVABLE: EXCLUDED_UNRESOLVED_MAPPING},
            unresolved_ke_pathways={KE_UNRESOLVABLE: [WP_UNRESOLVABLE]},
        )

    def test_round_trip_names_the_pathways(self):
        summary = ke_accounting_from_network(self._network())
        assert summary['excluded_unresolved_mapping'] == 1
        assert summary['unresolved_pathways'] == [WP_UNRESOLVABLE]
        assert WP_UNRESOLVABLE in format_ke_summary(summary)

    def test_round_trip_survives_json_serialisation(self):
        summary = ke_accounting_from_network(json.dumps(self._network()))
        assert summary['unresolved_pathways'] == [WP_UNRESOLVABLE]

    def test_networks_without_the_ids_report_none(self):
        """Stored before #81: the clause renders without IDs, not wrongly."""
        network = self._network()
        for node in network['nodes']:
            node['data'].pop('unresolved_pathways', None)
        summary = ke_accounting_from_network(network)
        assert summary['unresolved_pathways'] == []
        assert 'mapped, but no genes could be resolved' in format_ke_summary(summary)

    def test_unrelated_nodes_gain_no_new_keys(self):
        node = next(
            n['data'] for n in self._network()['nodes'] if n['data']['id'] == KE_OK
        )
        assert 'unresolved_pathways' not in node


class TestErrorExclusionKeepsItsMutedStyling:
    """Finding 4: only 'unresolved_mapping' withholds the `no-genes` class."""

    @staticmethod
    def _classes(reason):
        network = build_cytoscape_network(
            {KE_UNRESOLVABLE},
            pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID']),
            pd.DataFrame({'KE': [], 'p_value': [], 'FDR': []}),
            {KE_UNRESOLVABLE: 'Increase, ROS'},
            {KE_UNRESOLVABLE: 'AO'},
            reference_sets={},
            excluded_kes={KE_UNRESOLVABLE: reason} if reason else None,
        )
        return network['nodes'][0]['classes'].split()

    def test_error_without_a_gene_set_is_still_muted(self):
        assert 'no-genes' in self._classes(EXCLUDED_ERROR)

    def test_no_mapping_is_muted(self):
        assert 'no-genes' in self._classes(EXCLUDED_NO_MAPPING)

    def test_unresolved_mapping_is_styled_apart(self):
        classes = self._classes(EXCLUDED_UNRESOLVED_MAPPING)
        assert 'unresolved-mapping' in classes
        assert 'no-genes' not in classes


class TestTheTwoResultViewsAgree:
    """Finding 2: shared_results.html renders the same networks results.html
    does, so a class one styles and the other does not is a silent regression."""

    @staticmethod
    def _read(name):
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        with open(os.path.join(here, 'templates', name), encoding='utf-8') as fh:
            return fh.read()

    @pytest.mark.parametrize('template', ['results.html', 'shared_results.html'])
    def test_every_exclusion_class_is_styled_and_in_the_legend(self, template):
        body = self._read(template)
        for selector in ('node.no-genes', 'node.too-few-genes',
                         'node.unresolved-mapping'):
            assert selector in body, f"{template} does not style {selector}"
        assert 'Mapped, but no genes could be resolved' in body

    def test_the_no_gene_set_legend_wording_is_identical(self):
        wording = ('No measurable gene set (not curated, or none of its genes '
                   'measured)')
        assert wording in self._read('results.html')
        assert wording in self._read('shared_results.html')
