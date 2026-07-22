"""Issue #117 follow-up — an untested Key Event must not be painted in the
GSEA gradient's off-gradient plum.

The GSEA colour rule `selector: 'node[method = "gsea"]'` calls
`nesColor(ele.data('nes'))` on every node carrying `method = "gsea"`. An
untested Key Event still carries that attribute (with `nes = null`), and
`nesColor(null)` returns `NES_UNDEFINED_COLOR` — the plum reserved for a Key
Event that beat every permutation. Because the gradient rule is ordered after
the `.no-genes` / `.too-few-genes` / `.unresolved-mapping` rules, it overrides
their muted fills, so an uncurated or untestable KE (including the MIE and AO,
unmapped by design) would render as the single most salient node in every GSEA
network figure and .cyjs export.

These tests pin the fix: a per-exclusion override scoped to `[method="gsea"]`,
ordered *after* the gradient rule, that restores the muted fill. They are
static assertions over the template's Cytoscape stylesheet because the defect
is one of rule ordering in the cascade, which no Python-level render exercises.
"""
import re
from pathlib import Path

import pytest

TEMPLATES = Path(__file__).resolve().parent.parent / "templates"


def _stylesheet_selectors(html: str):
    """Ordered list of every `selector: '...'` string in the file."""
    return re.findall(r"selector:\s*'([^']+)'", html)


@pytest.mark.parametrize("template", ["results.html", "shared_results.html"])
class TestUntestedGseaNodeKeepsItsMutedFill:
    def _selectors(self, template):
        html = (TEMPLATES / template).read_text()
        return html, _stylesheet_selectors(html)

    @pytest.mark.parametrize(
        "cls,fill",
        [("no-genes", "#cccccc"), ("unresolved-mapping", "#f7e6d5"),
         ("too-few-genes", "#e8e8e8")],
    )
    def test_each_exclusion_has_a_gsea_scoped_background_override(
        self, template, cls, fill
    ):
        html, _ = self._selectors(template)
        # The override rule exists, is scoped to GSEA, and re-asserts the fill.
        block = re.search(
            r"selector:\s*'node\." + re.escape(cls)
            + r"\[method = \"gsea\"\]'\s*,\s*style:\s*\{([^}]*)\}",
            html,
        )
        assert block, f"{template}: no [method=gsea] override for .{cls}"
        assert fill in block.group(1), (
            f"{template}: .{cls} GSEA override does not restore {fill}"
        )

    @pytest.mark.parametrize(
        "cls", ["no-genes", "unresolved-mapping", "too-few-genes"]
    )
    def test_the_override_is_ordered_after_the_gradient_rule(self, template, cls):
        _, selectors = self._selectors(template)
        gradient = selectors.index('node[method = "gsea"]')
        override = selectors.index(f'node.{cls}[method = "gsea"]')
        # Cytoscape resolves equal-specificity conflicts by document order, so
        # the muted fill only wins if it comes after the gradient rule.
        assert override > gradient, (
            f"{template}: .{cls} GSEA override at {override} must follow the "
            f"gradient rule at {gradient}, or the plum fill wins"
        )
