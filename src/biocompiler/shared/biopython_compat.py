"""
Backward-compatibility shim: re-exports the biopython_compat module from
its new location in biocompiler.infrastructure.

The module was originally located at biocompiler.shared.biopython_compat
and was moved to biocompiler.infrastructure.biopython_compat during the
infrastructure consolidation. This shim preserves the old import path
for any external code or tests that still reference it.
"""

# Import the implementation module
from biocompiler.infrastructure import biopython_compat as _impl

# Re-export everything from the new location
from biocompiler.infrastructure.biopython_compat import *  # noqa: F401, F403

# Explicitly propagate __all__ so `from biocompiler.shared.biopython_compat import *`
# and `mod.__all__` both work correctly.
__all__ = _impl.__all__

# Also expose key names that tests import explicitly
from biocompiler.infrastructure.biopython_compat import (  # noqa: F401
    _check_biopython,
    to_seqrecord,
    from_seqrecord,
    optimize_to_seqrecord,
)
