"""Regression tests for the GSEA null-diagnostics mechanism (issue #117).

The classification itself is tested in ``test_gsea_service.py``. This file
covers the five ways the *mechanism around it* was shown to be wrong when the
first cut of #117 was reviewed:

R1  the capture patches a process-global and the app is concurrent, so two
    overlapping runs crossed each other's diagnostics and one leaked the patch;
R2  the batch / comparison path never learned about ``nes_status``, so the
    recovered cells stayed invisible in the cross-condition view;
R3  the network drew a recovered Key Event in the colour reserved for NES = 0;
R4  the docs and the UI claimed the nominal p-value survives a short
    same-signed tail — it does not, gseapy computes it on that same tail;
R5  a failure to capture was labelled ``ok``, so a gseapy rename would have
    reintroduced the whole bug under a clean bill of health.
"""
import importlib
import json
import math
import threading
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from config import Config
from services.comparison_service import (
    BEYOND_RESOLUTION_DOWN,
    BEYOND_RESOLUTION_UP,
    build_comparison_matrix,
    comparison_matrix_to_dataframe,
    comparison_nes_display_dataframe,
)
from services.enrichment_service import format_ke_summary
from services.gsea_service import (
    NES_BEYOND_RESOLUTION,
    NES_OK,
    NES_UNDIAGNOSED,
    NES_UNSTABLE,
    _capture_prerank_summaries,
    apply_null_diagnostics,
    run_gsea_analysis,
)
from services.network_service import build_cytoscape_network, ke_accounting_from_network
from services.report_service import format_nes, format_pvalue


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _one_signed_fixture(n_genes=200, set_size=100, heavy_tail=10, ke='KE:A'):
    """A ranking whose permutation null lands entirely below zero (#117).

    Same construction as ``test_gsea_service._one_signed_null_fixture``; kept
    local so the concurrency tests can vary the KE id and the set size.
    """
    genes = [f"G{i:04d}" for i in range(n_genes)]
    log2fcs = [1.0] * set_size + [-1.0] * (n_genes - set_size)
    pvals = ([1e-4] * set_size
             + [0.5] * (n_genes - set_size - heavy_tail)
             + [1e-60] * heavy_tail)
    df = pd.DataFrame({'ID': genes, 'log2FC': log2fcs, 'pval': pvals})
    return df, {ke: set(genes[:set_size])}, {ke}, {ke: 'Maximally enriched KE'}


def _healthy_fixture(ke='KE:A'):
    """A ranking whose permutation null straddles zero — the ordinary case."""
    rng = np.random.default_rng(7)
    genes = [f"H{i:04d}" for i in range(300)]
    df = pd.DataFrame({
        'ID': genes,
        'log2FC': rng.normal(0, 1.0, 300),
        'pval': rng.uniform(1e-3, 0.9, 300),
    })
    return df, {ke: set(genes[:40])}, {ke}, {ke: 'Ordinary KE'}


# ---------------------------------------------------------------------------
# R1 — the capture must not cross runs, and must not leak the patch
# ---------------------------------------------------------------------------

class TestCaptureUnderConcurrency:
    """R1: ``gseapy.gsea.prerank_rs`` is a process-global and this app is threaded."""

    def test_overlapping_captures_keep_their_own_summaries_and_restore(self):
        """Two overlapping captures must not see each other's summaries.

        Without a lock the two patches nest, so *both* sinks receive *both*
        runs' summaries. Since ``_null_tail_sizes`` is keyed by KE term — and
        the KE ids are identical across conditions — the second run's tail count
        then lands on the first run's row, where ``apply_null_diagnostics``
        rewrites that row's NES, p and FDR on the strength of it.

        The second defect is in the same block: the inner capture exits *after*
        the outer one here (the normal case, since conditions differ in size),
        so restoring "whatever was there on entry" leaves the global bound to a
        dead wrapper appending to an abandoned list for the life of the process.
        """
        gsea_module = importlib.import_module('gseapy.gsea')
        real = gsea_module.prerank_rs

        class _Result:
            def __init__(self, summaries):
                self.summaries = summaries

        def _fake_prerank_rs(tag, hold):
            time.sleep(hold)
            return _Result([tag])

        results = {}
        errors = []

        def _worker(tag, start_delay, hold):
            try:
                time.sleep(start_delay)
                sink = []
                # The sink itself is stored, not a copy of it: a nested wrapper
                # keeps appending to it long after its own block has returned,
                # and snapshotting here would miss exactly that.
                results[tag] = (sink, None)
                with _capture_prerank_summaries(sink) as state:
                    gsea_module.prerank_rs(tag, hold)
                results[tag] = (sink, state['captured'])
            except Exception as exc:  # pragma: no cover — surfaced by the assert
                errors.append(exc)

        gsea_module.prerank_rs = _fake_prerank_rs
        try:
            # 'first' finishes before 'second' — the out-of-order exit.
            threads = [
                threading.Thread(target=_worker, args=('first', 0.0, 0.15)),
                threading.Thread(target=_worker, args=('second', 0.05, 0.35)),
            ]
            for t in threads:
                t.start()
            for t in threads:
                t.join(timeout=30)
            assert not errors, errors

            # Read after both threads have finished, so a wrapper that outlives
            # its block and keeps appending is visible.
            assert results['first'] == (['first'], True), (
                "the first capture received another run's summaries — the "
                "process-global patch is not serialised"
            )
            assert results['second'] == (['second'], True), (
                "the second capture received another run's summaries — the "
                "process-global patch is not serialised"
            )
            assert gsea_module.prerank_rs is _fake_prerank_rs, (
                "the patch leaked: the module global is still bound to a "
                "wrapper appending to an abandoned sink"
            )
        finally:
            gsea_module.prerank_rs = real

    def test_two_gsea_runs_in_threads_get_their_own_tail_counts(self):
        """End to end: a degenerate and a healthy condition, run concurrently.

        Both use the same KE id, which is what makes the crossing invisible in
        production: every condition in a batch analyses the same Key Events.
        """
        degenerate = _one_signed_fixture()
        healthy = _healthy_fixture()

        serial = {
            'degenerate': run_gsea_analysis(*degenerate, permutation_num=100),
            'healthy': run_gsea_analysis(*healthy, permutation_num=100),
        }

        concurrent = {}
        errors = []

        def _worker(name, fixture):
            try:
                concurrent[name] = run_gsea_analysis(*fixture, permutation_num=100)
            except Exception as exc:  # pragma: no cover
                errors.append((name, exc))

        threads = [
            threading.Thread(target=_worker, args=('degenerate', degenerate)),
            threading.Thread(target=_worker, args=('healthy', healthy)),
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=120)
        assert not errors, errors

        for name in ('degenerate', 'healthy'):
            got = concurrent[name].iloc[0]
            want = serial[name].iloc[0]
            assert got['null_same_signed_n'] == want['null_same_signed_n'], (
                f"{name}: the same-signed tail count differs between a serial "
                "and a concurrent run — one condition's null was attributed to "
                "the other"
            )
            assert got['nes_status'] == want['nes_status']
            assert got['p_value'] == pytest.approx(want['p_value'])

        # And the two conditions really are in different regimes, so the test
        # would have caught a swap rather than comparing two identical rows.
        assert serial['degenerate'].iloc[0]['nes_status'] == NES_BEYOND_RESOLUTION
        assert serial['healthy'].iloc[0]['null_same_signed_n'] > 0

    def test_capture_is_reentrant_safe_against_a_foreign_patch(self):
        """A stranger rebinding the global mid-block still ends with the real one."""
        gsea_module = importlib.import_module('gseapy.gsea')
        real = gsea_module.prerank_rs

        with _capture_prerank_summaries([]):
            gsea_module.prerank_rs = lambda *a, **k: None
        assert gsea_module.prerank_rs is real


# ---------------------------------------------------------------------------
# R5 — a failed capture is not a clean result
# ---------------------------------------------------------------------------

class TestUndiagnosedIsNotOk:
    """R5: failure to capture must be distinguishable from a healthy null."""

    def test_a_renamed_entry_point_is_reported_loudly_not_silently(self, monkeypatch, caplog):
        """The capture itself, when the name it wraps is gone.

        Exercised on the context manager rather than through ``run_gsea_analysis``
        because gseapy's own ``Prerank.run`` resolves the same name: deleting it
        makes gseapy raise before the diagnostics matter. What this asserts is
        the contract the run depends on — no capture, ``captured`` False, and an
        error in the log naming the symbol.
        """
        gsea_module = importlib.import_module('gseapy.gsea')
        monkeypatch.delattr(gsea_module, 'prerank_rs', raising=True)

        sink = []
        with caplog.at_level('ERROR', logger='services.gsea_service'):
            with _capture_prerank_summaries(sink) as state:
                pass

        assert state['captured'] is False
        assert sink == []
        messages = ' '.join(r.getMessage() for r in caplog.records)
        assert 'prerank_rs' in messages
        assert NES_UNDIAGNOSED in messages

    def test_unreadable_summaries_label_every_row_undiagnosed(self, monkeypatch, caplog):
        """A capture that succeeds but yields nothing must not read as 'ok'.

        This is the shape a gseapy upgrade takes in practice: the wrapper still
        fires, but the summary objects no longer expose ``esnull``, so no tail
        count can be derived. Every row is then unchecked — and labelling those
        rows ``ok`` is exactly how #117 would come back under a clean bill of
        health.
        """
        import services.gsea_service as gsea_service
        monkeypatch.setattr(gsea_service, '_null_tail_sizes', lambda summaries: {})

        df, reference_sets, ke_list, ke_title_map = _healthy_fixture()
        with caplog.at_level('ERROR', logger='services.gsea_service'):
            result = run_gsea_analysis(
                df, reference_sets, ke_list, ke_title_map, permutation_num=50
            )

        assert set(result['nes_status']) == {NES_UNDIAGNOSED}
        assert NES_OK not in set(result['nes_status'])
        assert result['null_same_signed_n'].isna().all()
        assert result.attrs['ke_summary']['nes_undiagnosed_kes'] == len(result)
        messages = ' '.join(r.getMessage() for r in caplog.records)
        assert NES_UNDIAGNOSED in messages

    def test_undiagnosed_rows_are_surfaced_in_the_run_accounting(self):
        """The admission has to reach the reader, not only the server log."""
        res = pd.DataFrame({
            'KE': ['KE:1'], 'NES': [1.0], 'p_value': [1.0], 'FDR': [0.44],
        })
        apply_null_diagnostics(res, {}, 1000, diagnostics_captured=False)

        sentence = format_ke_summary({
            'tested': 1, 'total_kes': 1, 'nes_undiagnosed_kes': 1,
        })
        assert 'permutation null could not be inspected' in sentence
        assert 'unchecked' in sentence

    def test_accounting_sentence_is_silent_when_everything_was_diagnosed(self):
        """An ORA run must not grow a GSEA clause."""
        sentence = format_ke_summary({'tested': 5, 'total_kes': 5})
        assert 'permutation null' not in sentence


# ---------------------------------------------------------------------------
# R4 — the nominal p-value does NOT survive a short tail
# ---------------------------------------------------------------------------

class TestPValueIsCoarseOnAShortTail:
    """R4: the retracted claim, and the replacement for it."""

    def test_gseapy_computes_the_nominal_p_on_the_same_signed_tail(self):
        """The measurement the retraction rests on.

        If this ever fails, gseapy has changed its p-value denominator and the
        whole ``p_value_resolution`` column should be revisited — so it is
        asserted rather than left as a claim in a comment.
        """
        import gseapy as gp
        gsea_module = importlib.import_module('gseapy.gsea')

        genes = [f"G{i:04d}" for i in range(200)]
        df = pd.DataFrame({
            'ID': genes,
            'log2FC': [1.0] * 100 + [-1.0] * 100,
            'pval': [1e-4] * 100 + [0.5] * 90 + [1e-60] * 10,
        })
        metric = np.sign(df['log2FC']) * -np.log10(df['pval'].astype(float))
        rnk = pd.Series(metric.values, index=df['ID']).sort_values(ascending=False)

        sink = []
        with _capture_prerank_summaries(sink):
            res = gp.prerank(
                rnk=rnk, gene_sets={'K': list(genes[:10])}, outdir=None,
                min_size=5, max_size=1000, permutation_num=100, seed=1,
                threads=1, no_plot=True, verbose=False,
            )
        assert sink, "the null was not captured; the measurement cannot be made"

        summary = sink[0]
        es = float(summary.es)
        null = np.asarray(summary.esnull, dtype=float)
        same_signed = int((null >= 0).sum() if es >= 0 else (null < 0).sum())
        exceed = int((null >= es).sum() if es >= 0 else (null <= es).sum())
        reported = float(res.res2d.iloc[0]['NOM p-val'])

        assert same_signed > 0
        assert reported == pytest.approx(exceed / same_signed), (
            "gseapy's nominal p is no longer exceedances over the same-signed "
            "tail — p_value_resolution is computed on that assumption"
        )
        # The point of the retraction: the same-signed denominator is not the
        # permutation count, so the p-value is coarser than 1/permutation_num.
        assert same_signed < 100

    def test_unstable_rows_report_the_resolution_they_could_achieve(self):
        """A p of 0.6 on a tail of five is 3/5, and must say so."""
        res = pd.DataFrame({
            'KE': ['KE:TINY', 'KE:ONE', 'KE:FINE'],
            'NES': [1.1, 2.0, 1.5],
            'p_value': [0.6, 0.0, 0.012],
            'FDR': [0.6, 0.0, 0.02],
        })

        out = apply_null_diagnostics(res, {'KE:TINY': 5, 'KE:ONE': 1, 'KE:FINE': 500}, 1000)

        tiny = out[out['KE'] == 'KE:TINY'].iloc[0]
        assert tiny['nes_status'] == NES_UNSTABLE
        assert tiny['p_value_resolution'] == pytest.approx(0.2)
        one = out[out['KE'] == 'KE:ONE'].iloc[0]
        assert one['p_value_resolution'] == pytest.approx(1.0)
        fine = out[out['KE'] == 'KE:FINE'].iloc[0]
        assert fine['nes_status'] == NES_OK
        assert fine['p_value_resolution'] == pytest.approx(1 / 500)

    def test_beyond_resolution_reports_the_permutation_bound(self):
        res = pd.DataFrame({
            'KE': ['KE:MAX'], 'NES': [1.0], 'p_value': [1.0], 'FDR': [0.44],
        })

        out = apply_null_diagnostics(res, {'KE:MAX': 0}, 200)

        assert out.iloc[0]['p_value_resolution'] == pytest.approx(1 / 200)

    def test_report_qualifies_the_pvalue_not_only_the_nes(self):
        """The reports are what leaves the tool; the qualifier travels with it."""
        assert 'do not threshold' in format_pvalue(0.6, NES_UNSTABLE, 0.2)
        assert '0.2' in format_pvalue(0.6, NES_UNSTABLE, 0.2)
        assert format_pvalue(0.6, NES_OK) == '0.6000'
        assert 'bound' in format_pvalue(0.001, NES_BEYOND_RESOLUTION)
        assert 'not inspected' in format_pvalue(0.4, NES_UNDIAGNOSED)
        assert 'not inspected' in format_nes(1.2, NES_UNDIAGNOSED)

    def test_no_source_still_claims_the_pvalue_stands(self):
        """The retracted sentence must not survive anywhere it was written."""
        import pathlib
        root = pathlib.Path(__file__).resolve().parent.parent
        offenders = []
        for path in list(root.glob('services/*.py')) + list(root.glob('templates/*.html')) \
                + [root / 'CHANGELOG.md']:
            text = path.read_text(encoding='utf-8')
            for phrase in ('The p-value stands', 'the p-value stands',
                           'a sound p-value', 'the nominal p-value is\n    sound'):
                if phrase in text:
                    offenders.append((path.name, phrase))
        assert not offenders, offenders


# ---------------------------------------------------------------------------
# R2 — the cross-condition view
# ---------------------------------------------------------------------------

def _condition(label, enrichment, position=0):
    return SimpleNamespace(
        condition_label=label,
        enrichment_json=json.dumps(enrichment),
        ke_gene_json=None,
        position=position,
        gene_count=None,
        significant_genes=None,
        dose='',
        timepoint='',
    )


class TestComparisonCarriesTheStatus:
    """R2: the 5x9 grid in the manuscript is built from this path."""

    def _batch(self):
        beyond = {
            'KE': 'KE:1194', 'Title': 'Increase, DNA damage', 'FDR': 0.0,
            'NES': None, 'ES': 0.95, 'nes_status': NES_BEYOND_RESOLUTION,
        }
        ordinary = {
            'KE': 'KE:177', 'Title': 'Mitochondrial dysfunction', 'FDR': 0.01,
            'NES': -1.8, 'ES': -0.5, 'nes_status': NES_OK,
        }
        return [_condition('Carboplatin 48h', [beyond, ordinary])]

    def test_status_and_es_reach_the_matrix(self):
        m = build_comparison_matrix(self._batch(), method='gsea')

        row = m['ke_labels'].index('KE:1194')
        assert m['nes_matrix'][row][0] is None
        # Before this, a None NES here was indistinguishable from an untested KE.
        assert m['nes_status_matrix'][row][0] == NES_BEYOND_RESOLUTION
        assert m['es_matrix'][row][0] == pytest.approx(0.95)

        other = m['ke_labels'].index('KE:177')
        assert m['nes_status_matrix'][other][0] == NES_OK

    def test_ora_batches_carry_no_status(self):
        """An ORA batch has no NES and must not grow a phantom status."""
        conds = [_condition('C1', [{'KE': 'KE:1', 'Title': 'E', 'FDR': 0.01}])]
        m = build_comparison_matrix(conds, method='ora')
        assert m['nes_status_matrix'] == [[None]]
        assert m['es_matrix'] == [[None]]

    def test_report_table_marks_the_cell_instead_of_leaving_it_blank(self):
        """The captions claimed "every NES value" — this is what makes that true."""
        m = build_comparison_matrix(self._batch(), method='gsea')

        numeric = comparison_matrix_to_dataframe(m, which='nes')
        display = comparison_nes_display_dataframe(m)

        row = list(display['Key Event ID']).index('KE:1194')
        # The numeric export is still blank there — it is a numeric record.
        assert pd.isna(numeric['Carboplatin 48h'].iloc[row])
        # The report cell is not.
        assert display['Carboplatin 48h'].iloc[row] == BEYOND_RESOLUTION_UP
        other = list(display['Key Event ID']).index('KE:177')
        assert display['Carboplatin 48h'].iloc[other] == '-1.80'

    def test_display_marks_direction_from_the_es(self):
        cond = _condition('C1', [{
            'KE': 'KE:1', 'Title': 'Down', 'FDR': 0.0, 'NES': None,
            'ES': -0.9, 'nes_status': NES_BEYOND_RESOLUTION,
        }])
        m = build_comparison_matrix([cond], method='gsea')
        display = comparison_nes_display_dataframe(m)
        assert display['C1'].iloc[0] == BEYOND_RESOLUTION_DOWN
        ascii_display = comparison_nes_display_dataframe(m, ascii_only=True)
        assert ascii_display['C1'].iloc[0] == '- beyond res.'

    def test_unstable_and_undiagnosed_cells_are_flagged_in_the_report(self):
        cond = _condition('C1', [
            {'KE': 'KE:1', 'Title': 'A', 'FDR': 0.01, 'NES': 1.9, 'ES': 0.8,
             'nes_status': NES_UNSTABLE},
            {'KE': 'KE:2', 'Title': 'B', 'FDR': 0.02, 'NES': -1.2, 'ES': -0.4,
             'nes_status': NES_UNDIAGNOSED},
        ])
        m = build_comparison_matrix([cond], method='gsea')
        display = comparison_nes_display_dataframe(m)
        cells = dict(zip(display['Key Event ID'], display['C1']))
        assert cells['KE:1'] == '+1.90 !'
        assert cells['KE:2'] == '-1.20 ?'


# ---------------------------------------------------------------------------
# R3 — the network figure
# ---------------------------------------------------------------------------

class TestNetworkCarriesTheStatus:
    """R3: ``nesColor(null)`` renders white, the colour reserved for NES = 0."""

    def _network(self, nes_status=NES_BEYOND_RESOLUTION, nes=None, es=0.95):
        enrichment = pd.DataFrame([{
            'KE': 'KE:1194', 'Title': 'Increase, DNA damage', 'NES': nes,
            'ES': es, 'p_value': 0.001, 'FDR': 0.0, 'nes_status': nes_status,
        }])
        return build_cytoscape_network(
            ke_list={'KE:1194'},
            edges=pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID']),
            enrichment_results=enrichment,
            ke_title_map={'KE:1194': 'Increase, DNA damage'},
            ke_type_map={'KE:1194': 'intermediate'},
            reference_sets={'KE:1194': {'A', 'B'}},
            method='gsea',
        )

    def test_node_carries_status_and_direction(self):
        node = self._network()['nodes'][0]

        assert node['data']['nes'] is None
        # Without these two the frontend cannot tell a recovered Key Event from
        # a NES of exactly zero, and neither can the .cyjs export or the PNG.
        assert node['data']['nes_status'] == NES_BEYOND_RESOLUTION
        assert node['data']['es'] == pytest.approx(0.95)
        assert 'significant' in node['classes']

    def test_ordinary_gsea_node_is_unchanged_in_meaning(self):
        node = self._network(nes_status=NES_OK, nes=-1.8, es=-0.5)['nodes'][0]
        assert node['data']['nes'] == pytest.approx(-1.8)
        assert node['data']['nes_status'] == NES_OK

    def test_ora_nodes_carry_no_gsea_fields(self):
        enrichment = pd.DataFrame([{
            'KE': 'KE:1', 'Title': 'T', 'odds_ratio': 2.0,
            'p_value': 0.01, 'FDR': 0.02,
        }])
        network = build_cytoscape_network(
            ke_list={'KE:1'},
            edges=pd.DataFrame(columns=['Source_KE', 'Target_KE', 'KER_ID']),
            enrichment_results=enrichment,
            ke_title_map={'KE:1': 'T'},
            ke_type_map={'KE:1': 'intermediate'},
            reference_sets={'KE:1': {'A'}},
            method='ora',
        )
        assert 'nes_status' not in network['nodes'][0]['data']
        assert 'es' not in network['nodes'][0]['data']

    def test_undiagnosed_count_survives_into_a_stored_network(self):
        network = self._network(nes_status=NES_UNDIAGNOSED, nes=1.2, es=0.4)
        summary = ke_accounting_from_network(json.dumps(network))
        assert summary['nes_undiagnosed_kes'] == 1
        assert 'permutation null could not be inspected' in format_ke_summary(summary)
