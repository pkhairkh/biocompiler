"""
BioCompiler CSP/SMT Solver — Constraint-based gene optimization.

Replaces the greedy optimizer with a proper constraint satisfaction approach:
- OR-Tools CP-SAT as primary solver (automaton constraints for forbidden substrings)
- Greedy fallback when OR-Tools is unavailable

Key advantage: Global optimality — all constraints are considered simultaneously,
not sequentially. Cross-codon constraints (restriction sites, splice sites, CpG)
are handled natively via automaton/table constraints.

Note: The deprecated Z3 SMT backend (``engine_z3.py``) was removed in the
second-pass cleanup. ``SolverBackend.Z3`` remains defined as an enum value for
backward compatibility, but selecting it now falls through to the greedy
fallback.

Module structure:
    solver/
        __init__.py            — Public API (this file)
        types.py               — Shared data structures (SolverConfig, SolverResult, etc.)
        constraints.py         — Constraint builders (restriction sites, GC, etc.)
        automaton.py           — Automaton-based encoding for forbidden substrings
        dispatch.py            — Solver dispatch / fallback logic
        engine_ortools.py      — OR-Tools CP-SAT engine
        engine_greedy.py       — Greedy fallback engine
        mus.py                 — Minimal Unsatisfiable Subset extraction
        maxent_encoding.py     — MaxEntScan splice scoring as CSP constraints
        enforcement.py         — Constraint enforcement mechanics
        scoring.py             — Constraint enforcement scoring & soft constraint scoring & Pareto analysis
        conflict_resolution.py — Hard constraint conflict detection & resolution
        conflict_provenance.py  — Constraint conflict resolution provenance tracking
        constraint_interaction.py — Constraint interaction map (CAI cost analysis)

All public symbols from sub-modules are re-exported here so that downstream
code can import from ``biocompiler.solver`` directly::

    from biocompiler.solver import CSPSolver, SolverConfig, solve_with_csp
    from biocompiler.solver import ConstraintEnforcer, SoftConstraintScorer

Direct imports from sub-modules (e.g. ``from biocompiler.solver.dispatch import
solve_with_csp``) still work but are discouraged — they may trigger
deprecation warnings in a future release.
"""

from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING, Any, Optional, Union

if TYPE_CHECKING:
    from .engine_ortools import ORTOOLSEngine
    from .engine_greedy import GreedyEngine

    # Z3 backend (engine_z3.py) was removed in second-pass cleanup.
    # SolverBackend.Z3 still exists as an enum value but selecting it falls
    # through to the greedy fallback at runtime.
    EngineType = Union[ORTOOLSEngine, GreedyEngine]

# ═══════════════════════════════════════════════════════════════════════
# Core types (always available, no heavy deps)
# ═══════════════════════════════════════════════════════════════════════
from .types import (
    CodonVariable,
    ConstraintPriority,
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintType,
    ConstraintViolation,
    CrossCodonSpliceConstraint,
    CSPModel,
    InfeasibilityReport,
    MUSReport,
    SatisfiedConstraint,
    SOLVER_VERDICT_MAP,
    SolverBackend,
    SolverBackendProtocol,
    SolverConfig,
    SolverResult,
    SolverStatus,
    SpliceConstraint,
    enforcement_to_verdicts,
)

# ═══════════════════════════════════════════════════════════════════════
# Constraint model (always available, no heavy deps)
# ═══════════════════════════════════════════════════════════════════════
from .constraints import (
    # Abstract base classes
    HardConstraint,
    SoftConstraint,
    # Hard constraints
    TranslationConstraint,
    NoRestrictionSiteConstraint,
    GCRangeConstraint,
    NoCrypticSpliceConstraint,
    NoCpGIslandConstraint,
    NoATTTAMotifConstraint,
    NoTRunConstraint,
    # Soft constraints / objectives
    MaximizeCAI,
    MinimizeCpG,
    MinimizeCodonPairBias,
    MinimizeMRNADG,
    # CSPModel (constraints variant — re-exported for convenience;
    # the types.CSPModel is also re-exported above)
    CSPModel as _CSPModelConstraints,  # noqa: F401 — re-exported via types already
    # Factory function
    build_csp_model,
    # Helper functions
    codon_gc_count,
    codon_contains_gt,
    codon_contains_ag,
    codon_contains_cpg,
    compute_gc_from_codons,
)

# ═══════════════════════════════════════════════════════════════════════
# MaxEntScan encoding (always available, no heavy deps)
# ═══════════════════════════════════════════════════════════════════════
from .maxent_encoding import (
    SpliceConstraintEncoder,
    precompute_splice_scores,
    build_splice_constraint_table,
    quantize_maxent_scores,
    encode_cross_codon_splice_context,
)

# ═══════════════════════════════════════════════════════════════════════
# Dispatch (lazy heavy deps — always importable)
# ═══════════════════════════════════════════════════════════════════════
from .dispatch import (
    solve_with_csp,
    solve_csp,
    get_csp_availability,
    csp_optimize,
    validate_csp_solution,
)

# ═══════════════════════════════════════════════════════════════════════
# Enforcement (always available, no heavy deps)
# ═══════════════════════════════════════════════════════════════════════
from .enforcement import (
    ConstraintEnforcer,
    EnforcementResult,
    ConflictResolution as EnforcementConflictResolution,
    enforcement_to_verdicts as enforcement_to_verdicts_enforcement,
)

# ═══════════════════════════════════════════════════════════════════════
# Scoring (always available, no heavy deps)
# ═══════════════════════════════════════════════════════════════════════
from .scoring import (
    ConstraintScorer,
    ScoringResult,
    SoftConstraintScorer,
    compute_pareto_frontier,
)

# ═══════════════════════════════════════════════════════════════════════
# Conflict resolution (always available, no heavy deps)
# ═══════════════════════════════════════════════════════════════════════
from .conflict_resolution import (
    ConstraintConflict,
    ConflictResolver,
    prioritize_constraints,
)

# ═══════════════════════════════════════════════════════════════════════
# Conflict provenance (always available, no heavy deps)
# ═══════════════════════════════════════════════════════════════════════
from .conflict_provenance import (
    ConflictProvenance,
    ConflictResolverWithProvenance,
)

# ═══════════════════════════════════════════════════════════════════════
# Constraint interaction map (always available, no heavy deps)
# ═══════════════════════════════════════════════════════════════════════
from .constraint_interaction import (
    ConstraintInteractionMap,
    InteractionInfo,
    print_interaction_report,
)

# ═══════════════════════════════════════════════════════════════════════
# CAI-Aware Constraint Resolver (added by findings)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .cai_aware_resolver import CAIAwareConstraintResolver
except ImportError:
    CAIAwareConstraintResolver = None  # type: ignore[assignment,misc]

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# Deprecation warnings for old import paths
# ═══════════════════════════════════════════════════════════════════════

def __getattr__(name: str) -> Any:
    """Lazy access with deprecation warnings for old import paths.

    Emits a DeprecationWarning when a user imports a symbol directly
    from a sub-module instead of from the ``biocompiler.solver`` package.
    """
    # Map of old sub-module paths to their canonical package-level names
    _deprecated_aliases: dict[str, str] = {
        # Old names that may have been used before the package re-exports
        # were consolidated.  These map old names → canonical names.
        "is_csp_available": "get_csp_availability",
    }

    if name in _deprecated_aliases:
        canonical = _deprecated_aliases[name]
        warnings.warn(
            f"biocompiler.solver.{name} is deprecated — use "
            f"biocompiler.solver.{canonical} instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        # Delegate to the canonical name
        return globals()[canonical]

    raise AttributeError(f"module 'biocompiler.solver' has no attribute {name!r}")


# ═══════════════════════════════════════════════════════════════════════
# Solver availability checks
# ═══════════════════════════════════════════════════════════════════════

def is_solver_available(backend: Optional[SolverBackend] = None) -> bool:
    """Check whether a solver backend is available.

    Args:
        backend: Specific backend to check. If None, checks whether
            any available backend (OR-Tools or greedy fallback) is usable.

    Returns:
        True if the requested backend (or any backend) is importable.
    """
    if backend is not None:
        return _check_backend(backend)
    return _check_backend(SolverBackend.ORTOOLS) or _check_backend(SolverBackend.GREEDY_FALLBACK)


def _check_backend(backend: SolverBackend) -> bool:
    """Internal: check whether a single backend is importable."""
    if backend == SolverBackend.ORTOOLS:
        try:
            import ortools  # noqa: F401
            return True
        except ImportError:
            return False
    elif backend == SolverBackend.Z3:
        # Z3 backend (engine_z3.py) was removed in second-pass cleanup.
        # The enum value is retained for backward compatibility but the
        # backend is no longer wired up.
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
        """Instantiate the appropriate solver engine with fallback.

        Backend selection order: OR-Tools → greedy fallback.
        The Z3 backend (``engine_z3.py``) was removed in the second-pass
        cleanup; an explicit ``SolverBackend.Z3`` request now falls through
        to the greedy fallback with a warning.
        """
        backend = config.backend

        if backend == SolverBackend.Z3:
            logger.warning(
                "Z3 backend was removed in second-pass cleanup "
                "(engine_z3.py deleted); falling back to greedy."
            )
            backend = SolverBackend.GREEDY_FALLBACK

        if backend == SolverBackend.ORTOOLS:
            if _check_backend(SolverBackend.ORTOOLS):
                from .engine_ortools import ORTOOLSEngine
                return ORTOOLSEngine(config)
            logger.warning("OR-Tools not available, falling back to greedy")
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

    # ── Result & config types (from .types) ─────────────────
    "SolverResult",
    "SolverConfig",
    "SolverBackend",
    "SolverStatus",
    "SolverBackendProtocol",
    "CSPModel",

    # ── Supporting types (from .types) ──────────────────────
    "CodonVariable",
    "ConstraintPriority",
    "ConstraintSpec",
    "ConstraintStrictness",
    "ConstraintType",
    "ConstraintViolation",
    "CrossCodonSpliceConstraint",
    "InfeasibilityReport",
    "MUSReport",
    "SatisfiedConstraint",
    "SOLVER_VERDICT_MAP",
    "SpliceConstraint",
    "enforcement_to_verdicts",

    # ── Constraint model (from .constraints) ────────────────
    "HardConstraint",
    "SoftConstraint",
    "TranslationConstraint",
    "NoRestrictionSiteConstraint",
    "GCRangeConstraint",
    "NoCrypticSpliceConstraint",
    "NoCpGIslandConstraint",
    "NoATTTAMotifConstraint",
    "NoTRunConstraint",
    "MaximizeCAI",
    "MinimizeCpG",
    "MinimizeMRNADG",
    "build_csp_model",
    "codon_gc_count",
    "codon_contains_gt",
    "codon_contains_ag",
    "codon_contains_cpg",
    "compute_gc_from_codons",

    # ── MaxEntScan encoding (from .maxent_encoding) ─────────
    "SpliceConstraintEncoder",
    "precompute_splice_scores",
    "build_splice_constraint_table",
    "quantize_maxent_scores",
    "encode_cross_codon_splice_context",

    # ── Dispatch (from .dispatch) ───────────────────────────
    "solve_with_csp",
    "get_csp_availability",
    "csp_optimize",
    "validate_csp_solution",

    # ── Enforcement (from .enforcement) ─────────────────────
    "ConstraintEnforcer",
    "EnforcementResult",
    "EnforcementConflictResolution",

    # ── Scoring (from .scoring) ─────────────────────────────
    "ConstraintScorer",
    "ScoringResult",
    "SoftConstraintScorer",
    "compute_pareto_frontier",

    # ── Conflict resolution (from .conflict_resolution) ─────
    "ConstraintConflict",
    "ConflictResolver",
    "prioritize_constraints",

    # ── Conflict provenance (from .conflict_provenance) ──────
    "ConflictProvenance",
    "ConflictResolverWithProvenance",

    # ── Constraint interaction map (from .constraint_interaction) ──
    "ConstraintInteractionMap",
    "InteractionInfo",
    "print_interaction_report",

    # ── CAI-Aware resolver (from .cai_aware_resolver) ────────
    "CAIAwareConstraintResolver",

    # ── Convenience functions ───────────────────────────────
    "solve",
    "is_solver_available",
]
