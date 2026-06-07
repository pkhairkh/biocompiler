"""Deprecated: use biocompiler.optimizer.mutagenesis instead."""
import warnings

warnings.warn(
    "biocompiler.mutagenesis is deprecated — use biocompiler.optimizer.mutagenesis instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.optimizer.mutagenesis import *  # noqa: F401,F403
