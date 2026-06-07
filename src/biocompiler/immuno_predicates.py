"""Deprecated: use biocompiler.immunogenicity.predicates instead."""
import warnings

warnings.warn(
    "biocompiler.immuno_predicates is deprecated — use biocompiler.immunogenicity.predicates instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.immunogenicity.predicates import *  # noqa: F401,F403
