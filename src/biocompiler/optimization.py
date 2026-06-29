"""
Backward-compatibility shim for ``biocompiler.optimization``.

The canonical module is :mod:`biocompiler.optimizer`.  This module re-exports
the most commonly-used names so that legacy ``from biocompiler.optimization
import ...`` statements continue to work.

New code should import from :mod:`biocompiler.optimizer` (or its submodules)
directly.
"""

# Re-export everything from the optimizer package
from biocompiler.optimizer import *  # noqa: F401,F403

# Explicitly re-export names that are frequently imported from this shim
from biocompiler.optimizer import (
    optimize_sequence,
    batch_optimize,
    BioOptimizer,
    OptimizationResult,
)

from biocompiler.optimizer.cai import (
    _BatchSwapScorer,
    _compute_cai_fast,
    _codon_to_index,
    _adaptiveness_to_array,
    _count_dinucs_fast,
    HAS_NUMBA,
)

from biocompiler.optimizer.utils import FullConstructResult

__all__ = [
    "optimize_sequence",
    "batch_optimize",
    "BioOptimizer",
    "OptimizationResult",
    "_BatchSwapScorer",
    "_compute_cai_fast",
    "_codon_to_index",
    "_adaptiveness_to_array",
    "_count_dinucs_fast",
    "HAS_NUMBA",
    "FullConstructResult",
]
