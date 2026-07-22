"""Issue #104: server-side image rendering must have a deadline.

``Figure.to_image()`` drives a headless Chromium through kaleido. That
subprocess can wedge, and the calling thread then blocks with no deadline of
its own — the failure mode in #104, where a pytest run stalled at 15 minutes
against 40 seconds of CPU while a kaleido process tree stayed alive.

These tests cover the contract the callers rely on: a normal render passes its
arguments through untouched, and a render that never returns raises promptly
instead of hanging, so both report paths fall back to "image unavailable".
"""

import threading
import time

import pytest

from services.image_render import render_figure_png


class _FakeFigure:
    """Stand-in for a plotly Figure with a controllable ``to_image``."""

    def __init__(self, delay=0.0, result=b'PNG'):
        self.delay = delay
        self.result = result
        self.calls = []
        self.released = threading.Event()

    def to_image(self, **kwargs):
        self.calls.append(kwargs)
        if self.delay:
            # Emulate a wedged kaleido read: blocks until explicitly released.
            self.released.wait(self.delay)
        return self.result


class TestRenderFigurePng:
    """Bounded rendering behaviour."""

    def test_returns_bytes_and_forwards_kwargs(self):
        """A normal render returns the bytes and passes its kwargs through."""
        fig = _FakeFigure()

        out = render_figure_png(fig, width=600, height=400, scale=2)

        assert out == b'PNG'
        assert fig.calls == [{'format': 'png', 'width': 600, 'height': 400, 'scale': 2}]

    def test_defaults_to_png_format(self):
        """Callers that omit `format` still get a PNG (matches to_image usage)."""
        fig = _FakeFigure()

        render_figure_png(fig)

        assert fig.calls[0]['format'] == 'png'

    def test_wedged_render_raises_within_the_deadline(self):
        """A render that never returns raises TimeoutError, and raises quickly.

        This is the whole point: without the deadline the caller waits for the
        subprocess forever. The elapsed-time assertion is deliberately loose —
        it only has to prove the call did not wait for the 30s render.
        """
        fig = _FakeFigure(delay=30)

        started = time.monotonic()
        with pytest.raises(TimeoutError):
            render_figure_png(fig, timeout=0.2)
        elapsed = time.monotonic() - started

        assert elapsed < 5, f'timeout did not return promptly ({elapsed:.1f}s)'
        # Release the worker so it does not linger past the test.
        fig.released.set()

    def test_render_error_propagates_unchanged(self):
        """A genuine render failure is not masked as a timeout."""
        fig = _FakeFigure()
        fig.to_image = lambda **kwargs: (_ for _ in ()).throw(RuntimeError('no chrome'))

        with pytest.raises(RuntimeError, match='no chrome'):
            render_figure_png(fig)

    def test_timeout_defaults_to_config(self, monkeypatch):
        """With no explicit timeout the Config value bounds the render."""
        monkeypatch.setattr('services.image_render.Config.IMAGE_RENDER_TIMEOUT', 0.2)
        fig = _FakeFigure(delay=30)

        with pytest.raises(TimeoutError):
            render_figure_png(fig)

        fig.released.set()


class TestBatchHeatmapFallback:
    """The batch report degrades rather than hangs when a render times out."""

    def test_heatmap_returns_none_on_timeout(self, monkeypatch):
        """render_heatmap_png already swallows render failures — including this one."""
        from services import batch_report_service

        def _timeout(fig, **kwargs):
            raise TimeoutError('render did not complete')

        monkeypatch.setattr(batch_report_service, 'render_figure_png', _timeout)

        png = batch_report_service.render_heatmap_png({
            'method': 'ora',
            'ke_labels': ['KE:1'],
            'ke_titles': ['One'],
            'condition_labels': ['A'],
            'neg_log10_matrix': [[1.3]],
        })

        assert png is None
