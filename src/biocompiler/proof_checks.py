"""Deprecated: use biocompiler.provenance.proof_checks instead."""
import warnings

warnings.warn(
    "biocompiler.proof_checks is deprecated — use biocompiler.provenance.proof_checks instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.provenance.proof_checks import *  # noqa: F401,F403
