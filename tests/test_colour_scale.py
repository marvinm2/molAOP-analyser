"""The gene-expression legend must describe the ramp the code actually draws.

Until 2026-08-18 it did not, in two independent ways, and both survived every
review because nothing connected the legend to the function:

1. ``defaultLogFCColor()`` clamps at ±1, but the legend was labelled −2 / +2.
   A gene at log2FC 1.0 is already fully saturated red, so every reader
   mis-read every gene colour by a factor of two — and a +5 gene and a +1 gene
   are drawn identically with no indication that the scale has run out.
2. The legend's CSS gradient bar was a five-stop RdYlBu ramp with a **yellow**
   midpoint, while the function interpolates blue → **white** → red. The
   endpoints happened to agree, which is why nobody noticed: the error is only
   visible if you hold a mid-range node against the bar.

The expression overlay added in #51 generates its ramp *from* the function, so
it was correct on arrival — which put two mutually contradictory legends on the
same page and is what surfaced this.

These tests are deliberately textual. The colour logic lives in inline template
JavaScript (issue #77) with no JS test harness in this repo, so the cheap and
honest guard is to assert that the numbers and the stops in the markup still
match the constants in the script beside them. If the JS is ever extracted to
``static/js/``, point ``_script_text`` at the module instead of deleting this.
"""
import os
import re

import pytest

HERE = os.path.dirname(__file__)
ROOT = os.path.join(HERE, "..")

# Both pages embed their own copy of defaultLogFCColor() and their own legend.
# They are duplicates today (issue #77); test both so a fix to one that misses
# the other fails here rather than in someone's screenshot.
PAGES = ("templates/results.html", "templates/shared_results.html")
STYLESHEETS = ("static/css/pages/results.css", "static/css/pages/shared-results.css")

# `Math.max(-1, Math.min(1, fc))` — the clamp that defines where the ramp ends.
CLAMP_RE = re.compile(r"Math\.max\(\s*-(\d+(?:\.\d+)?)\s*,\s*Math\.min\(\s*(\d+(?:\.\d+)?)\s*,")

# The two endpoint colours interpolated by defaultLogFCColor(), as they appear
# in the function: interpolate(44, 255, t) etc. for the negative leg and
# interpolate(255, 215, t) etc. for the positive one.
NEGATIVE_LEG_RE = re.compile(r"interpolate\(44,\s*255,\s*t\)")
POSITIVE_LEG_RE = re.compile(r"interpolate\(255,\s*215,\s*t\)")


def _read(relpath):
    with open(os.path.join(ROOT, relpath), "r", encoding="utf-8") as fh:
        return fh.read()


@pytest.mark.parametrize("page", PAGES)
def test_clamp_is_symmetric_and_still_one(page):
    """The ramp saturates at ±1. If this changes, the labels below must too."""
    match = CLAMP_RE.search(_read(page))

    assert match, f"{page}: could not find defaultLogFCColor's clamp"
    low, high = match.group(1), match.group(2)
    assert low == high, f"{page}: asymmetric clamp -{low}/{high} — the legend assumes symmetry"
    assert high == "1", (
        f"{page}: the clamp is now ±{high}, but the legend is labelled ±1. "
        "Update both, and this assertion, together."
    )


@pytest.mark.parametrize("page", PAGES)
def test_legend_labels_match_the_clamp(page):
    """The bounds either side of the gradient bar must be the clamp, not double it."""
    html = _read(page)
    bar = html.index('class="network-legend__gradient-bar"')
    # The labels bracket the bar: one span before it, one after.
    before = html[max(0, bar - 400):bar]
    after = html[bar:bar + 400]

    assert "&minus;1" in before or ">-1<" in before, (
        f"{page}: lower legend bound is not −1. It read '-2' until 2026-08-18 "
        "while the code clamped at ±1."
    )
    assert ">+1<" in after, f"{page}: upper legend bound is not +1"
    assert ">-2<" not in before and ">+2<" not in after, (
        f"{page}: the ±2 labels are back; they do not match the ±1 clamp"
    )


@pytest.mark.parametrize("stylesheet", STYLESHEETS)
def test_gradient_bar_is_white_in_the_middle(stylesheet):
    """Zero is white in the function, so the legend bar must not show yellow."""
    css = _read(stylesheet)
    rule = css[css.index(".network-legend__gradient-bar"):]
    rule = rule[:rule.index("}")]

    assert "#ffffff" in rule.lower(), (
        f"{stylesheet}: the gradient bar has no white midpoint. "
        "defaultLogFCColor() interpolates through white at log2FC 0."
    )
    for banned in ("#ffffcc", "#fdae61", "#abd9e9"):
        assert banned not in rule.lower(), (
            f"{stylesheet}: {banned} is from the old RdYlBu bar, which no node "
            "is ever drawn in. The ramp is blue → white → red."
        )


@pytest.mark.parametrize("page", PAGES)
def test_gradient_bar_endpoints_match_the_function(page):
    """The bar's end colours must be the ones the function interpolates to.

    rgb(44,123,182) = #2c7bb6 at the negative end, rgb(215,25,28) = #d7191c at
    the positive one. Guards against the bar and the function drifting apart at
    the ends as they already did in the middle.
    """
    js = _read(page)
    assert NEGATIVE_LEG_RE.search(js), f"{page}: negative leg no longer starts at rgb(44,…)"
    assert POSITIVE_LEG_RE.search(js), f"{page}: positive leg no longer ends at rgb(215,…)"

    stylesheet = (
        "static/css/pages/results.css"
        if "shared" not in page
        else "static/css/pages/shared-results.css"
    )
    css = _read(stylesheet).lower()
    rule = css[css.index(".network-legend__gradient-bar"):]
    rule = rule[:rule.index("}")]
    assert "#2c7bb6" in rule, f"{stylesheet}: negative endpoint is not #2c7bb6"
    assert "#d7191c" in rule, f"{stylesheet}: positive endpoint is not #d7191c"
