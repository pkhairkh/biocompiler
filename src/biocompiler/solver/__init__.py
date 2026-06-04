"""
BioCompiler CSP/SMT Solver — Constraint-based gene optimization.

Replaces the greedy optimizer with a proper constraint satisfaction approach:
- OR-Tools CP-SAT as primary solver (automaton constraints for forbidden substrings)
- Z3 SMT as fallback solver (rich theory support for continuous constraints)
- Greedy fallback when solvers are unavailable

Key advantage: Global optimality — all constraints are considered simultaneously,
not sequentially. Cross-codon constraints (restriction sites, splice sites, CpG)
are handled natively via automaton/table constraints.

Module structure:
    solver/
        __init__.py        — Public API (this file)
        types.py           — Shared data structures (SolverConfig, SolverResult, etc.)
        constraints.py     — Constraint builders (restriction sites, GC, etc.)
        automaton.py       — Automaton-based encoding for forbidden substrings
        dispatch.py        — Solver dispatch / fallback logic
        engine_ortools.py  — OR-Tools CP-SAT engine
        engine_z3.py       — Z3 SMT engine
        mus.py             — Minimal Unsatisfiable Subset extraction
        maxent_encoding.py — MaxEntScan splice scoring as CSP constraints
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, Optional, Union

if TYPE_CHECKING:
    from .engine_ortools import ORTOOLSEngine
    from .engine_z3 import Z3Engine
    from .engine_greedy import GreedyEngine

    EngineType = Union[ORTOOLSEngine, Z3Engine, GreedyEngine]

# ── Core types (always available, no heavy deps) ────────────────────────
from .types import (
    CodonVariable,
    ConstraintType,
    ConstraintViolation,
    CrossCodonSpliceConstraint,
    CSPModel,
    MUSReport,
    SolverBackend,
    SolverConfig,
    SolverResult,
    SolverStatus,
    SpliceConstraint,
)

# ── MaxEntScan encoding (always available, no heavy deps) ───────────────
from .maxent_encoding import (
    SpliceConstraintEncoder,
    precompute_splice_scores,
    build_splice_constraint_table,
    quantize_maxent_scores,
    encode_cross_codon_splice_context,
)

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# Solver availability checks
# ═══════════════════════════════════════════════════════════════════════

def is_solver_available(backend: Optional[SolverBackend] = None) -> bool:
    """Check whether a solver backend is available.

    Args:
        backend: Specific backend to check. If None, checks whether
            ANY solver (OR-Tools or Z3) is available.

    Returns:
        True if the requested backend (or any backend) is importable.
    """
    if backend is not None:
        return _check_backend(backend)
    return _check_backend(SolverBackend.ORTOOLS) or _check_backend(SolverBackend.Z3)


def _check_backend(backend: SolverBackend) -> bool:
    """Internal: check whether a single backend is importable."""
    if backend == SolverBackend.ORTOOLS:
        try:
            import ortools  # noqa: F401
            return True
        except ImportError:
            return False
    elif backend == SolverBackend.Z3:
        try:
            import z3  # noqa: F401
            return True
        except ImportError:
            return False
    elif backend == SolverBackend.GREEDY_FALLBACK:
        return True
    else:
        return False


# ═══════════════════════════════════════════════════════════════════════
# Main solver class
# ═══════════════════════════════════════════════════════════════════════

class CSPSolver:
    """Constraint Satisfaction Problem solver for gene optimization.

    Encodes codon selection as a CSP where each codon position is a
    finite-domain variable, and biological constraints (GC range, no
    restriction sites, no cryptic splice sites, minimize CpG, maximize
    CAI) are encoded as hard/soft constraints.

    Usage:
        >>> config = SolverConfig(backend=SolverBackend.ORTOOLS)
        >>> solver = CSPSolver(config)
        >>> result = solver.solve("MVLSPADKTNVKAAWGKVGA")
        >>> print(result.solved, result.cai, result.gc_content)

    Args:
        config: Solver configuration. If None, uses defaults.
    """

    def __init__(self, config: Optional[SolverConfig] = None) -> None:
        self.config = config or SolverConfig()
        self._engine: Optional[EngineType] = None

    def solve(
        self,
        protein: str,
        organism: str = "Homo_sapiens",
        **overrides: Any,
    ) -> SolverResult:
        """Solve the codon optimization CSP for the given protein.

        Args:
            protein: Amino acid sequence (single-letter codes).
            organism: Target organism for codon usage tables.
            **overrides: Optional config overrides for this solve call.

        Returns:
            SolverResult with the optimized sequence and diagnostics.
        """
        config = self._merge_overrides(**overrides)
        engine = self._get_engine(config)
        return engine.solve(protein, organism)

    def _merge_overrides(self, **overrides: Any) -> SolverConfig:
        """Create a copy of config with per-call overrides applied."""
        if not overrides:
            return self.config
        import dataclasses
        changes = {k: v for k, v in overrides.items() if v is not None}
        return dataclasses.replace(self.config, **changes)

    def _get_engine(self, config: SolverConfig) -> EngineType:
        """Instantiate the appropriate solver engine with fallback."""
        backend = config.backend

        if backend == SolverBackend.ORTOOLS:
            if _check_backend(SolverBackend.ORTOOLS):
                from .engine_ortools import ORTOOLSEngine
                return ORTOOLSEngine(config)
            logger.warning("OR-Tools not available, falling back to Z3")
            backend = SolverBackend.Z3

        if backend == SolverBackend.Z3:
            if _check_backend(SolverBackend.Z3):
                from .engine_z3 import Z3Engine
                return Z3Engine(config)
            logger.warning("Z3 not available, falling back to greedy")
            backend = SolverBackend.GREEDY_FALLBACK

        if backend == SolverBackend.GREEDY_FALLBACK:
            from .engine_greedy import GreedyEngine
            return GreedyEngine(config)

        raise RuntimeError(f"No solver engine available for backend {backend}")


# ═══════════════════════════════════════════════════════════════════════
# Convenience function
# ═══════════════════════════════════════════════════════════════════════

def solve(
    protein: str,
    organism: str = "Homo_sapiens",
    config: Optional[SolverConfig] = None,
    **overrides: Any,
) -> SolverResult:
    """Convenience function: solve a codon optimization CSP in one call."""
    solver = CSPSolver(config)
    return solver.solve(protein, organism, **overrides)


# ═══════════════════════════════════════════════════════════════════════
# Public API exports
# ═══════════════════════════════════════════════════════════════════════

__all__ = [
    # ── Main solver class ───────────────────────────────────
    "CSPSolver",

    # ── Result & config types ───────────────────────────────
    "SolverResult",
    "SolverConfig",
    "SolverBackend",
    "SolverStatus",
    "CSPModel",

    # ── Supporting types ────────────────────────────────────
    "CodonVariable",
    "ConstraintType",
    "ConstraintViolation",
    "CrossCodonSpliceConstraint",
    "MUSReport",
    "SpliceConstraint",

    # ── MaxEntScan encoding ─────────────────────────────────
    "SpliceConstraintEncoder",
    "precompute_splice_scores",
    "build_splice_constraint_table",
    "quantize_maxent_scores",
    "encode_cross_codon_splice_context",

    # ── Convenience functions ───────────────────────────────
    "solve",
    "is_solver_available",
]
