"""
BioCompiler CSP Solver Dispatch
=================================

High-level API for the constraint-satisfaction solver pipeline.  This module
is the *only* entry point that downstream code (e.g. ``optimization.py``)
should import — it hides backend selection, model construction, fallback
logic, and result validation behind a small public surface.

Public API
----------
- ``solve_with_csp``        – main entry point for CSP optimisation
- ``get_csp_availability``  – runtime backend-availability probe
- ``csp_optimize``          – convenience wrapper matching ``_greedy_optimize`` signature
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
import warnings
from typing import TYPE_CHECKING, Any

from ..constants import AA_TO_CODONS
from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism

from .types import (
    ConstraintPriority,
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintViolation,
    CSPModel,
    SolverBackend,
    SolverConfig,
    SolverResult,
)
from .constraints import build_csp_model, HardConstraint, SoftConstraint
from .mus import compute_mus, quick_feasibility_check, FeasibilityReport
from .scoring import ConstraintScorer
from .conflict_provenance import (
    ConflictProvenance,
    ConflictResolverWithProvenance,
)

if TYPE_CHECKING:
    # These imports exist ONLY for TYPE_CHECKING — they provide type
    # annotations (e.g. type[ORTOOLSEngine] below) without importing
    # the heavy backend modules at runtime.  The actual runtime
    # imports happen in the try-except blocks below.
    from .engine_ortools import ORTOOLSEngine
    from .engine_z3 import Z3Engine

# Backend engines (optional — may not be installed)
_ORTOOLS_AVAILABLE: bool = False
_Z3_AVAILABLE: bool = False
_ortools_engine: type[ORTOOLSEngine] | None = None  # type: ignore[valid-type]
_z3_engine: type[Z3Engine] | None = None  # type: ignore[valid-type]

logger = logging.getLogger(__name__)

try:
    from .engine_ortools import ORTOOLSEngine as _ortools_engine  # noqa: F811
    _ORTOOLS_AVAILABLE = True
except ImportError:
    logger.debug("OR-Tools engine not available")

try:
    from .engine_z3 import Z3Engine as _z3_engine  # noqa: F811
    _Z3_AVAILABLE = True
except ImportError:
    logger.debug("Z3 engine not available")


def get_csp_availability() -> dict[str, bool]:
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


def is_csp_available() -> dict[str, bool]:
    """Deprecated: use :func:`get_csp_availability` or :func:`is_solver_available` instead."""
    warnings.warn(
        "is_csp_available() is deprecated — use get_csp_availability() "
        "or is_solver_available() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return get_csp_availability()


def is_solver_available() -> bool:
    """Check whether any CSP solver backend is available.

    Returns
    -------
    bool
        True if at least one backend (OR-Tools or Z3) is importable.
    """
    return get_csp_availability()["any"]


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


def _try_backend(
    engine_cls: type,
    model: CSPModel,
    config: SolverConfig,
    backend_name: str,
    backend_enum: SolverBackend,
) -> SolverResult | None:
    """Try solving with a single backend engine.

    Encapsulates the try-except logic that was previously duplicated for
    OR-Tools and Z3.  Returns a solved ``SolverResult`` with
    ``backend_used`` set, or ``None`` if the backend failed or returned
    infeasible.

    Parameters
    ----------
    engine_cls : type
        The engine class to instantiate (``ORTOOLSEngine`` or ``Z3Engine``).
    model : CSPModel
        The CSP model to solve.
    config : SolverConfig
        Solver configuration passed to the engine constructor.
    backend_name : str
        Human-readable name for log messages (e.g. ``"OR-Tools"``).
    backend_enum : SolverBackend
        Enum value to set on ``result.backend_used`` on success.

    Returns
    -------
    SolverResult | None
        A solved result with ``backend_used`` set, or ``None``.
    """
    logger.info("Attempting %s backend …", backend_name)
    try:
        engine = engine_cls(config)
        result = engine.solve(model)
        if result is not None and result.solved:
            result.backend_used = backend_enum
            logger.info("%s solved successfully.", backend_name)
            return result
        elif result is not None and result.fallback_used:
            # Backend returned a fallback result (e.g. Z3 unavailable) —
            # treat as failure so dispatch can try the next backend.
            logger.info(
                "%s returned fallback (reason: %s); trying next backend.",
                backend_name,
                result.metadata.get("reason", "unknown"),
            )
            return None
        else:
            logger.info("%s returned infeasible; trying next backend.", backend_name)
            return None
    except Exception as e:
        logger.warning("%s backend raised %s: %s", backend_name, type(e).__name__, e)
        return None


def solve_with_csp(
    protein: str,
    organism: str = "Homo_sapiens",
    config: SolverConfig | None = None,
    track_provenance: bool = False,
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
    track_provenance : bool
        Whether to record conflict resolution provenance.  When ``True``,
        every conflict and its resolution is tracked and stored in
        ``result.metadata["conflict_provenance"]`` as a list of
        :class:`ConflictProvenance` instances.  Defaults to ``False``
        for backward compatibility.
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
    # Resolve organism name to canonical form using centralized resolution
    canonical = resolve_organism(organism)
    if canonical not in CODON_ADAPTIVENESS_TABLES:
        logger.warning("Organism %r not found in CODON_ADAPTIVENESS_TABLES; using default weights.", organism)

    # Build configuration
    if config is None:
        config = SolverConfig(organism=organism)
    else:
        # Ensure config.organism is set from the explicit parameter
        config = SolverConfig(
            organism=organism,
            backend=config.backend,
            timeout_seconds=config.timeout_seconds,
            max_codons=config.max_codons,
            gc_lo=config.gc_lo,
            gc_hi=config.gc_hi,
            cryptic_splice_threshold=config.cryptic_splice_threshold,
            donor_threshold=config.donor_threshold,
            acceptor_threshold=config.acceptor_threshold,
            n_quantize_bins=config.n_quantize_bins,
            restriction_sites=config.restriction_sites,
            avoid_cpg=config.avoid_cpg,
            avoid_attta=config.avoid_attta,
            avoid_t_runs=config.avoid_t_runs,
            cai_weight=config.cai_weight,
            cpg_weight=config.cpg_weight,
            mrna_dg_weight=config.mrna_dg_weight,
            optimize_codon_pair_bias=config.optimize_codon_pair_bias,
            codon_pair_bias_weight=config.codon_pair_bias_weight,
            cai_reference_set=config.cai_reference_set,
            verbose=config.verbose,
            auto_detect_organism_domain=config.auto_detect_organism_domain,
        )

    # Build CSP model
    logger.info("Building CSP model for protein (%d aa), organism=%s", len(protein), organism)
    model = build_csp_model(protein, organism, config)

    # Provenance tracker (no-op when track_provenance=False)
    provenance_resolver = ConflictResolverWithProvenance(
        track_provenance=track_provenance,
        organism=organism,
    )
    provenance_records: list[ConflictProvenance] = []

    # Quick feasibility check
    report = quick_feasibility_check(model)
    if not report.feasible:
        reason = "; ".join(report.impossible_constraints) if report.impossible_constraints else "Model infeasible"
        logger.warning("CSP model infeasible: %s", reason)

        # Record provenance for each impossible constraint as a
        # relaxation tradeoff
        if track_provenance:
            for impossible_name in report.impossible_constraints:
                provenance_records.append(
                    provenance_resolver.record_relaxation_provenance(
                        relaxed_constraint_name=impossible_name,
                        kept_constraint_name="<model_feasibility>",
                        positions_affected=[],
                        sequence="",
                        resolution_method="csp_backtrack",
                    )
                )

        result = _make_fallback_result(protein, organism, time.monotonic() - start, reason)
        if provenance_records:
            result.metadata["conflict_provenance"] = provenance_records
        return result

    # Detect and record conflicts before solving (provenance)
    if track_provenance:
        constraints = list(getattr(model, "constraints", []))
        _, conflict_provenance = provenance_resolver.resolve_conflicts(
            constraints, "",  # no sequence yet
        )
        provenance_records.extend(conflict_provenance)

    # Try backends in priority order
    result: SolverResult | None = None

    # If the user explicitly requested GREEDY_FALLBACK, skip OR-Tools and Z3
    skip_csp_backends = config.backend == SolverBackend.GREEDY_FALLBACK

    # Respect explicit backend selection.  When the user requests a
    # specific backend (Z3 or ORTOOLS), try that one first.  If it
    # fails, fall through to the other CSP backend before going to
    # greedy.  This ensures that ``config.backend=Z3`` actually uses
    # Z3 when it's available, rather than always trying OR-Tools first.
    if config.backend == SolverBackend.Z3:
        # User explicitly requested Z3 — try it first
        if not skip_csp_backends and _z3_engine is not None:
            result = _try_backend(_z3_engine, model, config, "Z3", SolverBackend.Z3)
        if result is None and not skip_csp_backends and _ortools_engine is not None:
            result = _try_backend(_ortools_engine, model, config, "OR-Tools", SolverBackend.ORTOOLS)
    else:
        # Default priority: OR-Tools → Z3
        if not skip_csp_backends and _ortools_engine is not None:
            result = _try_backend(_ortools_engine, model, config, "OR-Tools", SolverBackend.ORTOOLS)
        if not skip_csp_backends and result is None and _z3_engine is not None:
            result = _try_backend(_z3_engine, model, config, "Z3", SolverBackend.Z3)

    # If CSP backends failed or were skipped, try greedy fallback engine
    if result is None:
        logger.info(
            "Attempting greedy fallback engine for organism=%s protein_len=%d",
            organism, len(protein),
        )
        try:
            from .engine_greedy import GreedyEngine
            greedy_engine = GreedyEngine(config)
            result = greedy_engine.solve(model)
            if result is not None and result.solved:
                logger.info(
                    "Greedy fallback engine solved successfully for organism=%s",
                    organism,
                )
            else:
                logger.warning("Greedy fallback engine returned unsolved result")
        except Exception as e:
            logger.error(
                "Greedy fallback engine raised %s: %s",
                type(e).__name__, e,
            )

    # Handle total failure (even greedy failed)
    if result is None or not result.solved:
        logger.error(
            "All solver backends (including greedy) failed for organism=%s protein_len=%d. "
            "OR-Tools available: %s, Z3 available: %s",
            organism, len(protein),
            _ortools_engine is not None, _z3_engine is not None,
        )
        result = _make_fallback_result(
            protein, organism, time.monotonic() - start,
            "All solver backends (including greedy) unavailable or infeasible",
        )
        if provenance_records:
            result.metadata["conflict_provenance"] = provenance_records
        return result

    # Post-solve validation (returns violations + composite enforcement score)
    if result.sequence:
        violations, composite_score = validate_csp_solution(result.sequence, protein, config, organism)
        result.violations = violations
        result.metadata["composite_score"] = composite_score
        if violations:
            logger.warning(
                "CSP solution has %d constraint violation(s) — "
                "composite_score=%.4f — solver may be incorrect.",
                len(violations), composite_score,
            )
            # Record provenance for each violation — these represent
            # constraints that were implicitly relaxed by the solver
            if track_provenance:
                violation_provenance = provenance_resolver.record_violation_provenance(
                    violations, result.sequence,
                )
                provenance_records.extend(violation_provenance)

    # Preserve fallback_used=True from greedy engine; CSP backends set False
    if result.backend_used != SolverBackend.GREEDY_FALLBACK:
        result.fallback_used = False
    result.solve_time_seconds = time.monotonic() - start

    # Store provenance in result metadata (if any records were captured)
    if provenance_records:
        result.metadata["conflict_provenance"] = provenance_records

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
            organism=organism,
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


# Default priority mapping for well-known constraint types, used when
# the constraint object doesn't carry an explicit ConstraintPriority.
_DEFAULT_PRIORITY_MAP: dict[str, ConstraintPriority] = {
    "TranslationConstraint": ConstraintPriority.CRITICAL,
    "NoRestrictionSiteConstraint": ConstraintPriority.HIGH,
    "GCRangeConstraint": ConstraintPriority.MEDIUM,
    "NoCrypticSpliceConstraint": ConstraintPriority.MEDIUM,
    "NoCpGIslandConstraint": ConstraintPriority.MEDIUM,
    "NoATTTAMotifConstraint": ConstraintPriority.MEDIUM,
    "NoTRunConstraint": ConstraintPriority.LOW,
    "MaximizeCAI": ConstraintPriority.LOW,
    "MinimizeCpG": ConstraintPriority.LOW,
    "MinimizeMRNADG": ConstraintPriority.LOW,
}


def _check_constraint(
    constraint: HardConstraint | SoftConstraint | ConstraintSpec,
    sequence: str,
    strictness: ConstraintStrictness,
) -> ConstraintViolation | None:
    """Check a single constraint against a candidate sequence.

    Handles both :class:`HardConstraint` / :class:`SoftConstraint` instances
    and :class:`ConstraintSpec` objects.  All three types now support a
    ``.check()`` method — ``ConstraintSpec.check()`` dispatches to the
    appropriate verification logic based on ``constraint.ctype`` and
    ``constraint.params``.

    The returned :class:`ConstraintViolation` (if any) includes the
    constraint's enforcement priority and weight.

    Parameters
    ----------
    constraint : HardConstraint | SoftConstraint | ConstraintSpec
        The constraint to check.
    sequence : str
        Candidate DNA sequence.
    strictness : ConstraintStrictness
        Whether this is a HARD or SOFT constraint.

    Returns
    -------
    ConstraintViolation | None
        A violation if the constraint is not satisfied or raises an error;
        ``None`` if the constraint passes.
    """
    # Determine priority and weight for the violation
    if isinstance(constraint, ConstraintSpec):
        priority = constraint.priority
        weight = constraint.weight
    else:
        priority = _DEFAULT_PRIORITY_MAP.get(
            constraint.name, ConstraintPriority.MEDIUM
        )
        weight = 1.0

    try:
        satisfied = constraint.check(sequence)
    except Exception as exc:
        return ConstraintViolation(
            constraint_name=constraint.name,
            constraint_type=strictness,
            description=f"Constraint check raised {type(exc).__name__}: {exc}",
            severity=0.5 if strictness == ConstraintStrictness.HARD else 0.3,
            priority=priority,
            weight=weight,
        )

    if not satisfied:
        severity = 0.8 if strictness == ConstraintStrictness.HARD else 0.5
        label = "Constraint" if strictness == ConstraintStrictness.HARD else "Soft constraint"
        return ConstraintViolation(
            constraint_name=constraint.name,
            constraint_type=strictness,
            description=f"{label} '{constraint.name}' is not satisfied by the solution.",
            severity=severity,
            priority=priority,
            weight=weight,
        )

    return None


def validate_csp_solution(
    sequence: str,
    protein: str,
    config: SolverConfig,
    organism: str = "Homo_sapiens",
) -> tuple[list[ConstraintViolation], float]:
    """Verify that a candidate sequence satisfies *all* constraints.

    This is a sanity check — if the solver is correct the returned
    violation list will be empty.  When violations are found they are
    reported so that the caller can decide whether to accept the
    solution, re-solve, or fall back.

    In addition to the violation list, a **composite enforcement score**
    in [0.0, 1.0] is returned.  This score is computed by
    :class:`ConstraintScorer` and reflects how well the solution
    satisfies all constraints, weighted by each constraint's priority
    and weight:

    - 1.0 = all constraints satisfied
    - 0.0 = at least one CRITICAL constraint violated
    - Between = partial satisfaction with priority-weighted penalties

    Iterates over every constraint in the CSP model and calls its
    ``check()`` method, then performs a translation-fidelity sanity
    check.  Uses the :func:`_check_constraint` helper so that both
    ``HardConstraint``/``SoftConstraint`` instances (with ``.check()``)
    and ``ConstraintSpec`` objects are handled gracefully.

    Parameters
    ----------
    sequence : str
        Candidate DNA sequence (already back-translated).
    protein : str
        The original amino-acid sequence.
    config : SolverConfig
        Solver configuration (defines the active constraint set).
    organism : str
        Target organism name.

    Returns
    -------
    tuple[list[ConstraintViolation], float]
        A tuple of (violations, composite_score).
        - violations: Empty if fully valid; otherwise one entry per
          violated constraint.
        - composite_score: Enforcement score in [0.0, 1.0].
    """
    violations: list[ConstraintViolation] = []

    if not sequence:
        violations.append(ConstraintViolation(
            constraint_name="non_empty_sequence",
            constraint_type=ConstraintStrictness.HARD,
            description="Sequence is empty — nothing to validate.",
            severity=1.0,
            priority=ConstraintPriority.CRITICAL,
        ))
        return violations, 0.0

    # Rebuild the model to iterate constraints.  build_csp_model() returns
    # a constraints.CSPModel whose hard_constraints / soft_constraints are
    # actual HardConstraint / SoftConstraint instances with .check() methods.
    constraint_model = build_csp_model(protein, organism, config)

    # Check hard constraints
    seen_names: set[str] = set()
    for constraint in constraint_model.hard_constraints:
        violation = _check_constraint(constraint, sequence, ConstraintStrictness.HARD)
        if violation is not None:
            violations.append(violation)
            seen_names.add(violation.constraint_name)

    # Check soft constraints
    for constraint in constraint_model.soft_constraints:
        violation = _check_constraint(constraint, sequence, ConstraintStrictness.SOFT)
        if violation is not None:
            violations.append(violation)
            seen_names.add(violation.constraint_name)

    # Also iterate constraint_model.constraints if present (e.g. types.CSPModel
    # has a list[ConstraintSpec] attribute).  ConstraintSpec objects now have
    # a .check() method that dispatches based on ctype and params.
    # Skip duplicates — constraints.CSPModel.constraints derives from
    # hard_constraints and soft_constraints, so the same logical
    # constraint may appear in both lists.
    for constraint in getattr(constraint_model, "constraints", []):
        # Deduplicate by constraint name to avoid double-counting
        constraint_name = getattr(constraint, "name", None)
        if constraint_name and constraint_name in seen_names:
            continue

        strictness = (
            constraint.strictness
            if isinstance(constraint, ConstraintSpec)
            else ConstraintStrictness.HARD
        )
        violation = _check_constraint(constraint, sequence, strictness)
        if violation is not None:
            violations.append(violation)
            seen_names.add(violation.constraint_name)

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
            priority=ConstraintPriority.CRITICAL,
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
                    priority=ConstraintPriority.CRITICAL,
                ))

    # ── Compute composite enforcement score via ConstraintScorer ─────
    # Use the constraint specs from the model (which carry priority & weight)
    constraint_specs: list[ConstraintSpec] = list(
        getattr(constraint_model, "constraints", [])
    )
    scorer = ConstraintScorer()
    composite_score = scorer.score_solution(sequence, constraint_specs)

    # If we found violations not in the specs (e.g. translation fidelity),
    # factor them into the composite score
    has_critical_violation = any(
        v.priority == ConstraintPriority.CRITICAL for v in violations
    )
    if has_critical_violation:
        composite_score = 0.0

    # Sort violations by priority (CRITICAL first) then severity
    violations.sort(key=lambda v: (v.priority.rank, -v.severity))

    logger.info(
        "validate_csp_solution: %d violation(s), composite_score=%.4f",
        len(violations), composite_score,
    )

    return violations, composite_score
