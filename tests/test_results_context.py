"""Tests for the shared results.html template context.

The batch-condition route once hand-maintained its own argument list and drifted
behind the live /analyze route: it omitted `ke_summary` and `hub_list`, both of
which reach `| tojson`, so the page raised
"TypeError: Object of type Undefined is not JSON serializable" — while the
variables that are merely falsy when undefined (the FDR selector, the KE
accounting, the resource-provenance line) disappeared without any error at all.

These tests pin the contract: the builder must supply every name the template
reads, and the template must render from the builder's output alone.
"""

import re

import pytest

from services.results_context import build_results_context


TEMPLATE = "templates/results.html"


def _minimal_batch_context():
    """The context a stored batch condition can supply — no volcano, no picker."""
    return build_results_context(
        table=[{"KE": "KE:177", "KE_Title": "Increase, Mitochondrial dysfunction",
                "p_value": 0.99, "FDR": 0.99, "Representation": "depleted"}],
        network={"nodes": [], "edges": []},
        metadata={"dataset_id": "Cisplatin 24h", "aop_id": "AOP:472"},
        volcano_data=[],
        background_size=12546,
        threshold=0.0,
        ke_gene_map={},
        ke_type_map={},
        ke_title_map={},
        method="ora",
        pval_threshold=0.01,
        ke_summary=None,
        hub_list=None,
        wp_picker_data=None,
    )


def test_builder_supplies_every_name_the_template_reads():
    """Every top-level Jinja name in results.html must come from the builder.

    Loop variables, `{% set %}` locals and Jinja globals are excluded — they are
    not route state.
    """
    with open(TEMPLATE) as fh:
        tpl = fh.read()

    read_names = (
        set(re.findall(r"\{\{\s*([a-z_][a-z0-9_]*)", tpl))
        | set(re.findall(r"\{%\s*if\s+([a-z_][a-z0-9_]*)", tpl))
    )
    # Names bound inside the template rather than passed into it.
    template_local = set(re.findall(r"\{%\s*set\s+([a-z_][a-z0-9_]*)", tpl))
    template_local |= set(re.findall(r"\{%\s*for\s+([a-z_][a-z0-9_,\s]*?)\s+in\s", tpl))
    template_local = {n.strip() for chunk in template_local for n in chunk.split(",")}
    jinja_builtins = {
        "url_for", "range", "loop", "dict", "none", "request", "config",
        "session", "super", "get_flashed_messages", "csrf_token", "not",
        # registered on app.jinja_env.globals
        "representation_labels", "resource_source_labels", "confidence_label",
    }

    expected = read_names - template_local - jinja_builtins
    supplied = set(_minimal_batch_context())

    assert not (expected - supplied), (
        "results.html reads names the shared context does not supply: "
        f"{sorted(expected - supplied)}"
    )


@pytest.mark.parametrize("key", ["ke_summary", "hub_list"])
def test_tojson_names_are_always_concrete(key):
    """`ke_summary` and `hub_list` hit `| tojson`; Undefined there is a 500."""
    ctx = _minimal_batch_context()
    assert key in ctx
    # None is fine for ke_summary (the template guards on it); Undefined is not.
    assert ctx[key] is None or isinstance(ctx[key], (list, dict))


def test_defaults_fill_in_what_a_stored_condition_cannot_provide():
    ctx = _minimal_batch_context()
    assert ctx["volcano_json"] == "[]"
    assert ctx["hub_list"] == []
    assert ctx["ke_summary_text"] == ""
    assert ctx["fdr_cutoff"] is not None
    assert ctx["fdr_choices"]


def test_pre_serialised_payloads_pass_through_unchanged():
    """The live route pre-serialises `table_json`; it must not be double-encoded."""
    ctx = build_results_context(
        table=[{"KE": "KE:1"}],
        table_json='[{"KE": "KE:1"}]',
        network='{"nodes": []}',
        metadata={},
    )
    assert ctx["table_json"] == '[{"KE": "KE:1"}]'
    assert ctx["network_json"] == '{"nodes": []}'


def test_results_template_renders_from_builder_output_alone():
    """End-to-end guard: the batch-shaped context must render without raising."""
    from app import app

    with app.test_request_context():
        from flask import render_template
        html = render_template("results.html", **_minimal_batch_context())

    assert "Increase, Mitochondrial dysfunction" in html


# ---------------------------------------------------------------------------
# Issue #103 — the discarded symbol-less row count
# ---------------------------------------------------------------------------

def _render(**overrides):
    """Render results.html from a batch-shaped context with overrides applied."""
    from app import app

    ctx = _minimal_batch_context()
    ctx.update(overrides)
    with app.test_request_context():
        from flask import render_template
        return render_template("results.html", **ctx)


def test_dropped_rows_defaults_to_not_recorded():
    """A caller that says nothing gets None — 'not recorded', not zero."""
    assert _minimal_batch_context()["dropped_rows"] is None


def test_discarded_rows_are_reported_next_to_the_background():
    """The count the loader already had must be visible on the page (#103)."""
    html = _render(background_size=12914, dropped_rows=875)

    assert "12914 measured genes" in html
    assert "875 rows had no gene symbol and were excluded" in html


def test_single_discarded_row_reads_in_the_singular():
    html = _render(dropped_rows=1)

    assert "1 row had no gene symbol and was excluded" in html


@pytest.mark.parametrize("value", [0, None])
def test_nothing_is_claimed_when_nothing_was_dropped_or_recorded(value):
    """Zero needs no clause, and an unrecorded count must not read as zero."""
    html = _render(dropped_rows=value)

    assert "had no gene symbol" not in html
