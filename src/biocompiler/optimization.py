"""
BioCompiler Optimizer v10.0.0
==============================
Multi-step certified gene optimization pipeline with aggressive GT resolution.

.. deprecated::
    Direct imports from ``biocompiler.optimization`` are deprecated as of v12.
    Use ``biocompiler.optimizer`` instead.  This module is retained as a
    thin shim for backward compatibility and will be removed in v13.

    Migration guide::

        # Old (deprecated):
        from biocompiler.optimization import BioOptimizer, optimize_sequence

        # New:
        from biocompiler.optimizer import BioOptimizer, optimize_sequence

All public symbols are re-exported below for backward compatibility.
"""

import warnings

warnings.warn(
    "Direct imports from 'biocompiler.optimization' are deprecated. "
    "Use 'biocompiler.optimizer' instead. "
    "The 'biocompiler.optimization' shim will be removed in v13.",
    DeprecationWarning,
    stacklevel=2,
)

# Re-export everything from the optimizer subpackage
from .optimizer import *  # noqa: F401,F403

# Also re-export the __all__ list
from .optimizer import __all__  # noqa: F811
