"""Deprecated: use biocompiler.engines.viennarna_fallback instead."""
import warnings

warnings.warn(
    "biocompiler.viennarna_fallback is deprecated — use biocompiler.engines.viennarna_fallback instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.engines.viennarna_fallback import *  # noqa: F401,F403
