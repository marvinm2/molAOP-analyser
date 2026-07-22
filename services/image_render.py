"""Bounded server-side image rendering for Plotly figures (issue #104).

``Figure.to_image()`` hands the figure to kaleido, which drives a headless
Chromium subprocess. That subprocess can wedge: it stays alive while the
calling thread sits in a blocking read with no deadline of its own, so a single
stuck render takes the request — or, in the test suite, the whole run — with it.
Issue #104 saw exactly that, a pytest run stalled at 15 minutes against 40
seconds of CPU with a live kaleido process tree.

Every ``to_image`` call therefore goes through :func:`render_figure_png`, which
gives the render a deadline. A wedged subprocess then produces a report with a
missing image and a warning in the log, which is what the callers already do for
every other render failure, instead of an unbounded hang.

The deadline cannot kill the orphaned Chromium — the worker thread stays blocked
until the process is reaped — but it does return control to the caller, which is
the difference between a degraded response and no response at all.
"""

import logging
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutureTimeoutError
from typing import Any, Optional

from config import Config

logger = logging.getLogger(__name__)


def render_figure_png(fig: Any, timeout: Optional[float] = None, **to_image_kwargs) -> bytes:
    """Render a Plotly figure to PNG bytes, giving kaleido a deadline.

    Args:
        fig: A ``plotly.graph_objects.Figure``.
        timeout: Seconds to wait before giving up. Defaults to
            ``Config.IMAGE_RENDER_TIMEOUT``.
        **to_image_kwargs: Passed straight through to ``fig.to_image()``
            (``format``, ``width``, ``height``, ``scale``, ...).

    Returns:
        bytes: the encoded image.

    Raises:
        TimeoutError: if the render did not finish within the deadline. Callers
            already treat a render failure as "image unavailable", so this
            surfaces through their existing handling.
        Exception: whatever ``to_image`` itself raised.
    """
    to_image_kwargs.setdefault('format', 'png')
    deadline = timeout if timeout is not None else Config.IMAGE_RENDER_TIMEOUT

    # A fresh executor per call: the point is to stop waiting on a wedged
    # render, and a shared pool would let one stuck worker starve the next
    # caller of its only thread. The thread is abandoned, not cancelled — see
    # the module docstring.
    executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix='kaleido')
    try:
        future = executor.submit(fig.to_image, **to_image_kwargs)
        try:
            return future.result(timeout=deadline)
        except FutureTimeoutError:
            logger.warning(
                'Image render exceeded %.0fs and was abandoned (#104); a kaleido '
                'subprocess may still be running', deadline
            )
            raise TimeoutError(
                f'Plotly image render did not complete within {deadline:.0f}s'
            )
    finally:
        # Never block on the abandoned worker: shutdown(wait=True) would
        # reintroduce the hang this function exists to prevent.
        executor.shutdown(wait=False)
