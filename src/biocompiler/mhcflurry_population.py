"""Deprecated: use biocompiler.immunogenicity.mhcflurry_population instead."""
import warnings

warnings.warn(
    "biocompiler.mhcflurry_population is deprecated — use biocompiler.immunogenicity.mhcflurry_population instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.immunogenicity.mhcflurry_population import *  # noqa: F401,F403
