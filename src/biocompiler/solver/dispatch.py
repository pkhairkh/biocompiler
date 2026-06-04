"""
BioCompiler CSP Solver Dispatch
=================================

High-level API for the constraint-satisfaction solver pipeline.  This module
is the *only* entry point that downstream code (e.g. ``optimization.py``)
should import — it hides backend selection, model construction, fallback
logic, and result validation behind a small public surface.

Public API
----------
- ``solve_with_csp``   – main entry point for CSP optimisation
- ``is_csp_available``  – runtime backend-availability probe
- ``csp_optimize``      – convenience wrapper matching ``_greedy_optimize`` signature
- ``validate_csp_solution`` – post-solve constraint verification

Design principles
-----------------
1. **Always importable** — the module never raises ``ImportError`` at the
   top level.  Heavy backends (OR-Tools, Z3) are imported lazily and
   their absence is recorded at runtime.
2. **Graceful degradation** — if every backend is unavailable the caller
   receives a ``SolverResult`` with ``fallback_used=True`` so the greedy
   path can take over seamlessly.
3. **Validate then trust** — every solution is sanity-checked against all
   constraints before being returned.
"""

from __future__ import annotations

import importlib.util
import logging
import time
from typing import Any

from ..constants import AA_TO_CODONS
from ..organisms import SPECIES

from .types import (
    ConstraintStrictness,
    ConstraintViolation,
    CSPModel,
    SolverBackend,
    SolverConfig,
    SolverResult,
)
from .constraints import build_csp_model
from .mus import compute_mus, quick_feasibility_check

# Backend engines (optional — may not be installed)
_ORTOOLS_AVAILABLE: bool = False
_Z3_AVAILABLE: bool = False
_ORTOOLSEngine: Any = None
_Z3Engine: Any = None

try:
    from .engine_ortools import ORTOOLSEngine as _ORTOOLSEngine_cls
    _ORTOOLS_AVAILABLE = True
    _ORTOOLSEngine = _ORTOOLSEngine_cls
except ImportError:
    pass

try:
    from .engine_z3 import Z3Engine as _Z3Engine_cls
    _Z3_AVAILABLE = True
    _Z3Engine = _Z3Engine_cls
except ImportError:
    pass

logger = logging.getLogger(__name__)


def is_csp_available() -> dict[str, bool]:
    """Check which CSP backends are importable *without* actually importing them.

    Uses ``importlib.util.find_spec`` so that the heavy modules (ortools, z3)
    are never loaded into the running process just to test availability.

    Returns
    -------
    dict[str, bool]
        ``{"ortools": bool, "z3": bool, "any": bool}``
    """
    ortools_ok = importlib.util.find_spec("ortools") is not None
    z3_ok = importlib.util.find_spec("z3") is not None
    return {"ortools": ortools_ok, "z3": z3_ok, "any": ortools_ok or z3_ok}


def is_solver_available() -> bool:
    """Check whether any CSP solver backend is available.

    Returns
    -------
    bool
        True if at least one backend (OR-Tools or Z3) is importable.
    """
    return is_csp_available()["any"]


def _make_fallback_result(
    protein: str, organism: str, solve_time: float, reason: str,
) -> SolverResult:
    """Construct a SolverResult that signals the greedy path should be used."""
    return SolverResult(
        sequence="",
        solved=False,
        backend_used=SolverBackend.NONE,
        protein=protein,
        organism=organism,
        fallback_used=True,
        solve_time_seconds=solve_time,
        violations=[],
        metadata={"reason": reason},
    )


def solve_with_csp(
    protein: str,
    organism: str = "Homo_sapiens",
    config: SolverConfig | None = None,
    **kwargs: Any,
) -> SolverResult:
    """Main entry point for CSP-based codon optimization.

    Backend priority: OR-Tools → Z3 → fallback (greedy).

    Parameters
    ----------
    protein : str
        Amino-acid sequence (single-letter codes, e.g. ``"MVLSPADKTN"``).
    organism : str
        Target organism name.
    config : SolverConfig | None
        Full solver configuration.  If ``None``, a default is constructed.
    **kwargs
        Ignored when *config* is provided; otherwise reserved for future use.

    Returns
    -------
    SolverResult
        Optimized solution or a fallback indicator.

    Raises
    ------
    ValueError
        If *protein* is empty or contains invalid amino-acid codes.
    """
    start = time.monotonic()

    # Validate inputs
    if not protein or not protein.strip():
        raise ValueError("Protein sequence must be a non-empty string.")
    protein = protein.upper().strip()
    valid_aas = set(AA_TO_CODONS.keys())
    invalid = set(ch for ch in protein if ch not in valid_aas)
    if invalid:
        raise ValueError(f"Invalid amino-acid codes in protein: {invalid}.")
    if organism not in SPECIES:
        logger.warning("Organism %r not found in SPECIES; using default weights.", organism)

    # Build configuration
    if config is None:
        config = SolverConfig()

    # Build CSP model
    logger.info("Building CSP model for protein (%d aa), organism=%s", len(protein), organism)
    model = build_csp_model(protein, organism, config)

    # Quick feasibility check
    feasible, reason = quick_feasibility_check(model)
    if not feasible:
        logger.warning("CSP model infeasible: %s", reason)
        return _make_fallback_result(protein, organism, time.monotonic() - start, reason)

    # Try backends in priority order
    result: SolverResult | None = None

    # --- OR-Tools ---
    if _ORTOOLS_AVAILABLE and _ORTOOLSEngine is not None:
        logger.info("Attempting OR-Tools backend …")
        try:
            engine = _ORTOOLSEngine(config)
            result = engine.solve(model)
            if result is not None and result.solved:
                result.backend_used = SolverBackend.ORTOOLS
                logger.info("OR-Tools solved successfully.")
            else:
                logger.info("OR-Tools returned infeasible; trying next backend.")
                result = None
        except Exception:
            logger.warning("OR-Tools backend raised an exception", exc_info=True)
            result = None

    # --- Z3 ---
    if result is None and _Z3_AVAILABLE and _Z3Engine is not None:
        logger.info("Attempting Z3 backend …")
        try:
            engine = _Z3Engine(config)
            result = engine.solve(model)
            if result is not None and result.solved:
                result.backend_used = SolverBackend.Z3
                logger.info("Z3 solved successfully.")
            else:
                logger.info("Z3 returned infeasible.")
                result = None
        except Exception:
            logger.warning("Z3 backend raised an exception", exc_info=True)
            result = None

    # Handle total failure → fallback
    if result is None:
        logger.warning("All CSP backends failed; falling back to greedy.")
        return _make_fallback_result(
            protein, organism, time.monotonic() - start,
            "All CSP backends unavailable or infeasible",
        )

    # Post-solve validation
    if result.sequence:
        violations = validate_csp_solution(result.sequence, protein, config)
        result.violations = violations
        if violations:
            logger.warning(
                "CSP solution has %d constraint violation(s) — solver may be incorrect.",
                len(violations),
            )

    result.fallback_used = False
    result.solve_time_seconds = time.monotonic() - start
    return result


def csp_optimize(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: list[str] | None = None,
    cryptic_splice_threshold: float = 3.0,
    avoid_cpg: bool = True,
    avoid_attta: bool = True,
    max_homopolymer: int = 5,
    timeout_seconds: float = 120.0,
    **kwargs: Any,
) -> SolverResult:
    """Convenience wrapper matching the ``_greedy_optimize()`` signature.

    Translates the greedy-optimizer's parameter set into a ``SolverConfig``
    and delegates to :func:`solve_with_csp`.  On failure the returned
    ``SolverResult`` has ``fallback_used=True``.

    Parameters
    ----------
    protein : str
        Amino-acid sequence (single-letter codes).
    organism : str
        Target organism name.
    gc_lo, gc_hi : float
        Acceptable GC-content range ``[0, 1]``.
    restriction_sites : list[str] | None
        Restriction-enzyme recognition sequences to avoid.
    cryptic_splice_threshold : float
        Maximum allowed cryptic-splice score.
    avoid_cpg : bool
        Whether to avoid CpG dinucleotides.
    avoid_attta : bool
        Whether to avoid ATTTA instability motifs.
    max_homopolymer : int
        Maximum homopolymer run length allowed.
    timeout_seconds : float
        Solver wall-clock timeout.
    **kwargs
        Additional keyword arguments (reserved for future use).

    Returns
    -------
    SolverResult
        CSP solution, or a fallback indicator on failure.
    """
    try:
        config = SolverConfig(
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cryptic_splice_threshold=cryptic_splice_threshold,
            restriction_sites=restriction_sites or [],
            avoid_cpg=avoid_cpg,
            avoid_attta=avoid_attta,
            avoid_t_runs=max_homopolymer >= 6,
            timeout_seconds=timeout_seconds,
        )
    except Exception:
        logger.warning("Failed to build SolverConfig; returning fallback result.", exc_info=True)
        return _make_fallback_result(protein, organism, 0.0, "SolverConfig construction failed")

    try:
        return solve_with_csp(protein, organism, config=config)
    except Exception:
        logger.warning("solve_with_csp failed; returning fallback result.", exc_info=True)
        return _make_fallback_result(protein, organism, 0.0, "solve_with_csp raised an exception")


def validate_csp_solution(
    sequence: str,
    protein: str,
    config: SolverConfig,
) -> list[ConstraintViolation]:
    """Verify that a candidate sequence satisfies *all* constraints.

    This is a sanity check — if the solver is correct the returned list
    will be empty.  When violations are found they are reported so that
    the caller can decide whether to accept the solution, re-solve, or
    fall back.

    Iterates over every constraint in the CSP model and calls its
    ``check()`` method, then performs a translation-fidelity sanity check.

    Parameters
    ----------
    sequence : str
        Candidate DNA sequence (already back-translated).
    protein : str
        The original amino-acid sequence.
    config : SolverConfig
        Solver configuration (defines the active constraint set).

    Returns
    -------
    list[ConstraintViolation]
        Empty if fully valid; otherwise one entry per violated constraint.
    """
    violations: list[ConstraintViolation] = []

    if not sequence:
        violations.append(ConstraintViolation(
            constraint_name="non_empty_sequence",
            constraint_type=ConstraintStrictness.HARD,
            description="Sequence is empty — nothing to validate.",
            severity=1.0,
        ))
        return violations

    # Rebuild the model to iterate constraints
    model: CSPModel = build_csp_model(protein, config.organism, config)

    for constraint in model.constraints:
        try:
            satisfied = constraint.check(sequence)
        except Exception as exc:
            violations.append(ConstraintViolation(
                constraint_name=constraint.name,
                constraint_type=ConstraintStrictness.HARD,
                description=f"Constraint check raised {type(exc).__name__}: {exc}",
                severity=0.5,
            ))
            continue

        if not satisfied:
            violations.append(ConstraintViolation(
                constraint_name=constraint.name,
                constraint_type=ConstraintStrictness.HARD,
                description=f"Constraint '{constraint.name}' is not satisfied by the solution.",
                severity=0.8,
            ))

    # Extra sanity: translation fidelity
    from ..constants import CODON_TABLE

    if len(sequence) != len(protein) * 3:
        violations.append(ConstraintViolation(
            constraint_name="sequence_length",
            constraint_type=ConstraintStrictness.HARD,
            description=(
                f"Sequence length ({len(sequence)}) does not equal "
                f"protein length × 3 ({len(protein) * 3})."
            ),
            severity=1.0,
        ))
    else:
        codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
        for idx, (codon, expected_aa) in enumerate(zip(codons, protein)):
            actual_aa = CODON_TABLE.get(codon)
            if actual_aa != expected_aa:
                violations.append(ConstraintViolation(
                    constraint_name="translation_fidelity",
                    constraint_type=ConstraintStrictness.HARD,
                    description=f"Codon {idx} ({codon}) translates to '{actual_aa}', expected '{expected_aa}'.",
                    severity=1.0,
                ))

    return violations
