"""The `hidden` attribute must actually hide, whatever main.css sets display to.

main.css gives select, number inputs and buttons `display: block`, which beats
the browser's own `[hidden] { display: none }`. In production that showed the
"At least N sources" box under union and intersection, and the drop zone's
Remove button before any file was chosen, although the templates hide both.
"""
import re
from pathlib import Path

MAIN_CSS = Path(__file__).resolve().parent.parent / "static" / "css" / "main.css"


def test_main_css_forces_hidden_to_display_none():
    css = MAIN_CSS.read_text(encoding="utf-8")
    rule = re.search(r"\[hidden\]\s*\{([^}]*)\}", css)
    assert rule, "main.css has no [hidden] rule"
    assert re.search(r"display\s*:\s*none\s*!important", rule.group(1))


def test_hidden_rule_precedes_the_display_block_rule():
    """Not required for !important, but keeps the intent next to its cause."""
    css = MAIN_CSS.read_text(encoding="utf-8")
    assert css.index("[hidden]") < css.index('select, input[type="number"], button {')
