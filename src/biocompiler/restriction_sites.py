"""Deprecated: use biocompiler.sequence.restriction_sites instead."""
import warnings

warnings.warn(
    "biocompiler.restriction_sites is deprecated — use biocompiler.sequence.restriction_sites instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.restriction_sites import *  # noqa: F401,F403

__all__ = [
    "RESTRICTION_SITES",
    "MIN_SITE_LENGTH",
    "get_recognition_site",
    "expand_iupac_site",
    "get_eliminable_sites",
]
