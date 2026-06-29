"""OR-Tools CP-SAT Solver Engine for BioCompiler.

Primary solver backend using Google OR-Tools CP-SAT with automaton constraints
for cross-codon pattern avoidance. This is the default and most powerful backend
for constraint-based gene optimization.

Architecture:
    The engine encodes codon optimization as a Constraint Satisfaction Problem:
    - Decision variables: one IntVar per codon position (domain = synonymous codons)
    - Nucleotide extraction: element constraints link codon vars to base vars
    - Pattern avoidance: automaton (DFA) constraints on the nucleotide sequence
    - Cross-codon coordination: table constraints on adjacent codon pairs
    - Objective: maximize sum of log(CAI) weighted by codon adaptiveness

The CP-SAT solver handles:
    - Exact GC content range enforcement
    - Restriction site elimination via composite DFA
    - Cryptic splice site avoidance (GT/AG context + cross-codon)
    - ATTTA instability motif avoidance
    - T-run (6+ consecutive T) avoidance
    - CpG dinucleotide minimization
    - CAI maximization as objective

Fallback handling:
    - INFEASIBLE → MUS (Minimal Unsatisfiable Subset) report
    - MODEL_INVALID → greedy fallback
    - Timeout → best incumbent solution with fallback_used=True
"""

from __future__ import annotations

import logging
import math
import time
from dataclasses import dataclass
from typing import Any, Optional

from biocompiler.solver.types import (
    CodonVariable,
    InfeasibilityReport,
    MUSReport,
    SolverBackend,
    SolverConfig,
    SolverResult,
)
from .automaton import (
    build_composite_dfa,
    build_trun_dfa,
    dfa_to_ortools_format,
)
from biocompiler.shared.constants import (
    AA_TO_CODONS,
    BASE_REV,
    RESTRICTION_ENZYMES,
    INSTABILITY_MOTIF,
    reverse_complement,
)
from ..organisms import (
    CODON_ADAPTIVENESS_TABLES,
    SUPPORTED_ORGANISMS,
)

logger = logging.getLogger(__name__)

__all__ = [
    "ORTOOLSEngine",
    "CSPModel",
    "greedy_fallback",
]

# CAI scaling factor: multiply float log-CAI by this to get integer objective
_CAI_SCALE = 10000

# Maximum number of transitions in a DFA before we refuse to add the automaton
# constraint (CP-SAT has practical limits on automaton size)
_MAX_AUTOMATON_TRANSITIONS = 100000

# Default number of parallel search workers for CP-SAT solver
_DEFAULT_NUM_WORKERS = 8

# Threshold below which an objective bound is considered negative infinity
_NEGATIVE_INFINITY_THRESHOLD = -1e18

# Maximum protein length (in codons) for ATTTA automaton constraint
_ATTTA_AUTOMATON_MAX_CODONS = 50

# Maximum consecutive T bases allowed before T-run avoidance kicks in
_MAX_T_RUN_LENGTH = 5

# Minimum fraction of cross-codon pairs that must remain viable (1/N);
# if fewer than total_pairs // N pairs survive, skip the constraint
_MIN_VIABLE_PAIR_DENOMINATOR = 4

# Nucleotide alphabet size (A, C, G, T)
_NUCLEOTIDE_ALPHABET_SIZE = 4

# Default CAI weight for codons missing from adaptiveness tables
_DEFAULT_CAI_WEIGHT = 0.01

# Threshold protein length (in codons) for decomposed automaton encoding;
# sequences longer than this use chunked automaton constraints instead of
# a single large automaton for better solver performance
_DECOMPOSED_AUTOMATON_CODON_THRESHOLD = 500

# Threshold protein length (in codons) for reduced CAI scaling;
# long proteins can cause integer overflow in the CP-SAT objective,
# so we use a smaller scale factor
_LONG_PROTEIN_CODON_THRESHOLD = 1000

# Reduced CAI scaling factor for long proteins (avoids overflow)
_CAI_SCALE_LONG = 1000

# Threshold protein length (in codons) for very long proteins;
# proteins beyond this use an even smaller CAI scale to prevent
# integer overflow in the CP-SAT objective sum
_VERY_LONG_PROTEIN_CODON_THRESHOLD = 2000

# Very reduced CAI scaling factor for very long proteins (avoids overflow)
_CAI_SCALE_VERY_LONG = 100

# Chunk size (in nucleotides) for decomposed automaton encoding
_AUTOMATON_CHUNK_SIZE = 600


class ORTOOLSEngine:
    """OR-Tools CP-SAT solver engine for constraint-based gene optimization.

    This is the primary solver backend. It uses CP-SAT's automaton constraints
    to enforce forbidden nucleotide patterns (restriction sites, ATTTA motifs,
    T-runs) and table constraints for cross-codon coordination (splice sites,
    CpG avoidance).

    The automaton module builds DFAs that **accept** strings NOT containing
    forbidden patterns.  OR-Tools' ``AddAutomaton`` constraint requires the
    nucleotide sequence to end in one of the specified final (accepting)
    states.  Because our DFAs already designate safe states as accepting,
    we pass the DFA's accepting list directly.

    Example usage::

        from biocompiler.solver import SolverConfig, SolverBackend
        from biocompiler.solver.engine_ortools import ORTOOLSEngine

        config = SolverConfig(backend=SolverBackend.ORTOOLS)
        engine = ORTOOLSEngine(config)
        result = engine.solve(csp_model)
        if result.solved:
            print(f"Optimized: {result.sequence}, CAI={result.cai:.3f}")
    """

    def __init__(self, config: SolverConfig) -> None:
        """Initialize the OR-Tools engine with solver configuration.

        Args:
            config: Solver configuration controlling constraints, timeout,
                and objective weights.
        """
        self.config = config
        self._unavailable_reason: str = ""

    # ------------------------------------------------------------------
    # Static availability check
    # ------------------------------------------------------------------

    @staticmethod
    def is_available() -> bool:
        """Check whether OR-Tools is importable.

        Returns:
            True if ``ortools.sat.python.cp_model`` can be imported, False otherwise.
        """
        try:
            from ortools.sat.python import cp_model  # noqa: F401
            return True
        except ImportError:
            return False
        except Exception:
            # OR-Tools may be installed but have a broken native extension
            # (e.g. ABI mismatch, missing shared library).  Catch broadly
            # so the engine degrades gracefully instead of crashing.
            return False

    def _check_available(self) -> bool:
        """Check availability and store the reason if unavailable.

        Unlike the static :meth:`is_available`, this method stores the
        unavailability reason on ``self._unavailable_reason`` so that
        callers can diagnose why the engine is not usable.

        Returns:
            True if OR-Tools is usable, False otherwise.
        """
        try:
            from ortools.sat.python import cp_model  # noqa: F401
            self._unavailable_reason = ""
            return True
        except ImportError as exc:
            self._unavailable_reason = f"ImportError: {exc}"
            return False
        except Exception as exc:
            self._unavailable_reason = f"{type(exc).__name__}: {exc}"
            return False

    # ------------------------------------------------------------------
    # Main solve method
    # ------------------------------------------------------------------

    def solve(self, model: CSPModel) -> SolverResult:
        """Solve the CSP model using OR-Tools CP-SAT.

        Constructs a CP-SAT model with:
        1. IntVar per codon position (domain = synonymous codon indices)
        2. Nucleotide variable extraction via element constraints
        3. GC content range constraint
        4. Restriction site automaton constraint
        5. Splice site constraints (GT/AG avoidance tables)
        6. ATTTA motif automaton constraint
        7. T-run automaton constraint
        8. CpG avoidance constraints
        9. CAI maximization objective

        Args:
            model: The CSP model encoding the gene design problem.

        Returns:
            SolverResult with the optimized sequence and diagnostics.
            Check ``result.solved`` before using the sequence.
        """
        start_time = time.monotonic()

        # ── Availability check (stores reason for diagnostics) ────────
        if not self._check_available():
            logger.error("OR-Tools not available: %s", self._unavailable_reason)
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.ORTOOLS,
                fallback_used=True,
                solve_time_seconds=time.monotonic() - start_time,
                warnings=[
                    f"OR-Tools not available: {self._unavailable_reason}; "
                    "cannot solve with CP-SAT backend"
                ],
            )

        from ortools.sat.python import cp_model as ortools_cp

        # ── Safe attribute access with fallbacks ──────────────────────
        protein = getattr(model, "protein", None)
        if protein is None:
            # types.CSPModel uses protein_sequence instead of protein
            protein = getattr(model, "protein_sequence", "")
        organism = getattr(model, "organism", "") or getattr(self.config, "organism", "unknown")

        # Sync organism from model to config so that is_eukaryotic detection
        # works correctly even when only the model specifies the organism.
        if organism and organism != self.config.organism:
            try:
                object.__setattr__(self.config, "organism", organism)
            except (AttributeError, TypeError):
                # Frozen dataclass — cannot modify; log warning
                logger.warning(
                    "SolverConfig.organism='%s' differs from model organism='%s' "
                    "but config is frozen; eukaryotic detection may be incorrect",
                    self.config.organism, organism,
                )

        # Log model construction diagnostics
        codon_domains_attr = getattr(model, "codon_domains", {})
        n_codon_domains = len(codon_domains_attr) if codon_domains_attr else 0
        model_constraints = getattr(model, "constraints", [])
        logger.info(
            "OR-Tools engine diagnostics: protein_len=%d organism=%s "
            "codon_domains=%d model_constraints=%d "
            "config.gc_range=[%.2f,%.2f] config.avoid_cpg=%s "
            "config.avoid_t_runs=%s config.avoid_attta=%s "
            "config.cryptic_splice_threshold=%.1f",
            len(protein) if protein else 0, organism,
            n_codon_domains, len(model_constraints),
            self.config.gc_lo, self.config.gc_hi,
            self.config.avoid_cpg, self.config.avoid_t_runs,
            self.config.avoid_attta, self.config.cryptic_splice_threshold,
        )

        if not protein:
            logger.error("OR-Tools solve called with empty/missing protein attribute")
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.ORTOOLS,
                solve_time_seconds=time.monotonic() - start_time,
                warnings=["Model has no protein sequence (model.protein is empty/missing)"],
            )

        n_codons = len(protein)

        # ── Diagnostic logging ────────────────────────────────────────
        hard_constraints = getattr(model, "hard_constraints", [])
        logger.info(
            "OR-Tools solving: protein_len=%d organism=%s "
            "hard_constraints=%d config.avoid_cpg=%s config.avoid_t_runs=%s",
            n_codons, organism,
            len(hard_constraints) if hard_constraints else 0,
            self.config.avoid_cpg, self.config.avoid_t_runs,
        )

        if n_codons == 0:
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.ORTOOLS,
                solve_time_seconds=time.monotonic() - start_time,
                warnings=["Empty protein sequence"],
            )

        if n_codons > self.config.max_codons:
            logger.warning(
                "Protein length (%d) exceeds max_codons (%d); "
                "CP-SAT model too large, use greedy fallback",
                n_codons, self.config.max_codons,
            )
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.ORTOOLS,
                fallback_used=True,
                solve_time_seconds=time.monotonic() - start_time,
                warnings=[
                    f"Sequence length ({n_codons}) exceeds max_codons "
                    f"({self.config.max_codons}); use greedy fallback"
                ],
            )

        # CP-SAT practical limit warning: models with >2000 codons can
        # exhaust memory or take unreasonably long even with parallel
        # workers.  Log a warning but still attempt to solve.
        if n_codons > 2000:
            logger.warning(
                "Large protein (%d codons): CP-SAT solve may be slow "
                "or exceed memory; consider greedy fallback",
                n_codons,
            )

        # Build sorted codon lists per amino acid (CAI-descending order)
        adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})
        sorted_codons_per_aa: dict[str, list[str]] = {}
        for aa in set(protein):
            codons = AA_TO_CODONS.get(aa, [])
            sorted_codons_per_aa[aa] = sorted(
                codons, key=lambda c: adaptiveness.get(c, 0.0), reverse=True
            )

        try:
            # ── 1. Create CP-SAT model and decision variables ──────────
            cp_model = ortools_cp.CpModel()
            codon_vars: list[ortools_cp.IntVar] = []
            codon_domains: list[list[str]] = []  # parallel to codon_vars

            for i, aa in enumerate(protein):
                domain = sorted_codons_per_aa.get(aa, AA_TO_CODONS.get(aa, []))
                if not domain:
                    return SolverResult(
                        sequence="",
                        solved=False,
                        backend_used=SolverBackend.ORTOOLS,
                        solve_time_seconds=time.monotonic() - start_time,
                        warnings=[f"No codons for amino acid '{aa}' at position {i}"],
                    )
                var = cp_model.NewIntVar(0, len(domain) - 1, f"codon_{i}")
                codon_vars.append(var)
                codon_domains.append(domain)

            # ── 2. Build nucleotide variables via element constraints ──
            nucleotide_vars = self._build_nucleotide_variables(
                cp_model, codon_vars, codon_domains, protein
            )

            # ── 3. Add GC constraint ──────────────────────────────────
            gc_count_var = self._add_gc_constraint(
                cp_model, codon_vars, codon_domains, n_codons
            )
            if self.config.verbose:
                logger.info(
                    "VERBOSE: GC constraint active: gc_range=[%.2f, %.2f]",
                    self.config.gc_lo, self.config.gc_hi,
                )

            # ── 4. Add composite automaton constraint (restriction + ATTTA + T-run) ──
            # Combine all forbidden-nucleotide-pattern DFAs into ONE automaton
            # constraint for better solver performance and feasibility.
            num_constraints = self._add_composite_automaton_constraint(
                cp_model, nucleotide_vars, n_codons
            )

            # ── 5. Add splice site constraints ────────────────────────
            n_splice = self._add_splice_constraints(
                cp_model, codon_vars, codon_domains, protein
            )
            num_constraints += n_splice
            if self.config.verbose:
                logger.info(
                    "VERBOSE: Splice constraints: %d added "
                    "(threshold=%.1f, is_eukaryotic=%s)",
                    n_splice, self.config.cryptic_splice_threshold,
                    self.config.is_eukaryotic,
                )

            # ── 6. Add CpG avoidance constraints ──────────────────────
            n_cpg = self._add_cpg_constraints(
                cp_model, codon_vars, codon_domains, protein
            )
            num_constraints += n_cpg
            if self.config.verbose:
                logger.info(
                    "VERBOSE: CpG constraints: %d added (avoid_cpg=%s)",
                    n_cpg, self.config.avoid_cpg,
                )

            # ── 9. Set CAI maximization objective ─────────────────────
            # Use progressively reduced CAI scale for long/very-long proteins
            # to avoid integer overflow in the CP-SAT objective sum.
            if n_codons > _VERY_LONG_PROTEIN_CODON_THRESHOLD:
                cai_scale = _CAI_SCALE_VERY_LONG
            elif n_codons > _LONG_PROTEIN_CODON_THRESHOLD:
                cai_scale = _CAI_SCALE_LONG
            else:
                cai_scale = _CAI_SCALE
            if n_codons > _LONG_PROTEIN_CODON_THRESHOLD:
                logger.info(
                    "Long protein (%d codons): using reduced CAI scale=%d "
                    "to avoid integer overflow",
                    n_codons, cai_scale,
                )
            self._set_cai_objective(cp_model, codon_vars, codon_domains, organism, cai_scale)

        except Exception as build_exc:
            logger.error(
                "OR-Tools model construction failed: %s: %s",
                type(build_exc).__name__, build_exc, exc_info=True,
            )
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.ORTOOLS,
                fallback_used=True,
                solve_time_seconds=time.monotonic() - start_time,
                warnings=[
                    f"OR-Tools model construction failed "
                    f"({type(build_exc).__name__}: {build_exc}); "
                    "fall back to greedy optimizer"
                ],
            )

        # ── 10. Solve ─────────────────────────────────────────────────
        logger.info(
            "OR-Tools CP-SAT model built: %d variables, %d constraints, "
            "solving with timeout=%.1fs",
            n_codons, num_constraints, self.config.timeout_seconds,
        )

        solver = ortools_cp.CpSolver()
        solver.parameters.max_time_in_seconds = self.config.timeout_seconds
        solver.parameters.num_workers = _DEFAULT_NUM_WORKERS  # Use parallel search

        if self.config.verbose:
            solver.parameters.log_search_progress = True

        try:
            status = solver.Solve(cp_model)
        except Exception as solve_exc:
            logger.error(
                "OR-Tools Solve() raised %s: %s",
                type(solve_exc).__name__, solve_exc, exc_info=True,
            )
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.ORTOOLS,
                fallback_used=True,
                solve_time_seconds=time.monotonic() - start_time,
                num_constraints=num_constraints,
                num_variables=n_codons,
                warnings=[
                    f"OR-Tools Solve() raised {type(solve_exc).__name__}: {solve_exc}; "
                    "fall back to greedy optimizer"
                ],
            )

        solve_time = time.monotonic() - start_time

        # Map status to human-readable name for logging
        status_names = {
            ortools_cp.OPTIMAL: "OPTIMAL",
            ortools_cp.FEASIBLE: "FEASIBLE",
            ortools_cp.INFEASIBLE: "INFEASIBLE",
            ortools_cp.MODEL_INVALID: "MODEL_INVALID",
            ortools_cp.UNKNOWN: "UNKNOWN",
        }
        status_name = status_names.get(status, f"UNKNOWN({status})")
        logger.info(
            "OR-Tools solve completed: status=%s time=%.2fs "
            "variables=%d constraints=%d",
            status_name, solve_time, n_codons, num_constraints,
        )

        # ── 11. Process result ────────────────────────────────────────
        if status == ortools_cp.OPTIMAL or status == ortools_cp.FEASIBLE:
            # Extract solution
            sequence = self._extract_solution(
                solver, codon_vars, codon_domains, protein
            )

            # Validate: if extraction produced an empty sequence, treat as failure
            if not sequence:
                logger.error(
                    "OR-Tools: Solver reported FEASIBLE but solution "
                    "extraction produced empty sequence"
                )
                return SolverResult(
                    sequence="",
                    solved=False,
                    backend_used=SolverBackend.ORTOOLS,
                    fallback_used=True,
                    solve_time_seconds=solve_time,
                    num_constraints=num_constraints,
                    num_variables=n_codons,
                    warnings=[
                        "Solver reported feasible but solution extraction "
                        "yielded empty sequence — falling back"
                    ],
                )

            # Compute metrics
            gc_val = self._compute_gc(sequence)
            cai_val = self._compute_cai(sequence, protein, organism)
            obj_val = solver.ObjectiveValue()

            is_optimal = (status == ortools_cp.OPTIMAL)
            fallback = (status == ortools_cp.FEASIBLE and not is_optimal)

            result = SolverResult(
                sequence=sequence,
                solved=True,
                backend_used=SolverBackend.ORTOOLS,
                cai=cai_val,
                gc_content=gc_val,
                solve_time_seconds=solve_time,
                objective_value=obj_val / cai_scale,
                num_constraints=num_constraints,
                num_variables=n_codons,
                fallback_used=fallback,
                warnings=[] if is_optimal else [
                    "Solver found feasible solution but not proven optimal "
                    f"(timeout={self.config.timeout_seconds}s)"
                ],
            )
            return result

        elif status == ortools_cp.INFEASIBLE:
            # Compute MUS report
            mus = self._compute_mus(protein, organism, codon_domains)
            # Build InfeasibilityReport from MUS for user-facing diagnostics
            infeasibility = InfeasibilityReport.from_mus_report(
                mus, detection_method="cp_sat_infeasible"
            )
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.ORTOOLS,
                solve_time_seconds=solve_time,
                mus_report=mus,
                infeasibility_report=infeasibility,
                num_constraints=num_constraints,
                num_variables=n_codons,
                warnings=["Problem is INFEASIBLE — conflicting constraints detected"],
            )

        elif status == ortools_cp.MODEL_INVALID:
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.ORTOOLS,
                fallback_used=True,
                solve_time_seconds=solve_time,
                num_constraints=num_constraints,
                num_variables=n_codons,
                warnings=[
                    "CP-SAT model is invalid (likely constraint conflict); "
                    "fall back to greedy optimizer"
                ],
            )

        else:
            # UNKNOWN status — timeout or other
            # Try to extract any incumbent solution the solver may have found
            # before timing out. Both OPTIMAL and FEASIBLE are already handled
            # above; UNKNOWN means the solver could not prove feasibility or
            # infeasibility within the time limit, but may still have a
            # solution available.  Always attempt to return the best feasible
            # solution even on timeout.
            try:
                sequence = self._extract_solution(
                    solver, codon_vars, codon_domains, protein
                )
                if sequence:
                    gc_val = self._compute_gc(sequence)
                    cai_val = self._compute_cai(sequence, protein, organism)
                    # ObjectiveValue() may raise if no solution was found;
                    # fall back to NaN if unavailable.
                    try:
                        obj_val = solver.ObjectiveValue() / cai_scale
                    except Exception:
                        obj_val = float("nan")
                    return SolverResult(
                        sequence=sequence,
                        solved=True,
                        backend_used=SolverBackend.ORTOOLS,
                        cai=cai_val,
                        gc_content=gc_val,
                        solve_time_seconds=solve_time,
                        objective_value=obj_val,
                        num_constraints=num_constraints,
                        num_variables=n_codons,
                        fallback_used=True,
                        warnings=[
                            "Solver timed out; returning best incumbent solution"
                        ],
                    )
            except Exception:
                logger.warning("OR-Tools engine error", exc_info=True)

            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.ORTOOLS,
                fallback_used=True,
                solve_time_seconds=solve_time,
                num_constraints=num_constraints,
                num_variables=n_codons,
                warnings=["Solver returned UNKNOWN status (likely timeout)"],
            )

    # ==================================================================
    # Helper: Nucleotide variable construction
    # ==================================================================

    def _build_nucleotide_variables(
        self,
        cp_model: Any,
        codon_vars: list[Any],
        codon_domains: list[list[str]],
        protein: str,
    ) -> list[Any]:
        """Create nucleotide IntVars linked to codon variables via element constraints.

        For each codon position i and nucleotide offset j (0,1,2), we create
        a variable n_{i*3+j} with domain {0,1,2,3} (A=0, C=1, G=2, T=3) and
        link it to the codon variable using AddElement:
            codon_vars[i] = k  =>  n_{i*3+j} = encoding_of(domain[k][j])

        Args:
            cp_model: OR-Tools CpModel instance.
            codon_vars: List of codon IntVar variables.
            codon_domains: Parallel list of codon string lists (domains).
            protein: The amino acid sequence.

        Returns:
            List of nucleotide IntVar variables (length = 3 * len(protein)).
        """
        nucleotide_vars: list = []

        for i, (cvar, domain) in enumerate(zip(codon_vars, codon_domains)):
            # For each of the 3 nucleotide positions in this codon
            for j in range(3):
                nuc_var = cp_model.NewIntVar(0, 3, f"nuc_{i * 3 + j}")

                # Build the value table: for each codon index k,
                # what is the base encoding at position j?
                values = [BASE_REV[domain[k][j]] for k in range(len(domain))]

                # Element constraint: domain[cvar] selects which value
                cp_model.AddElement(cvar, values, nuc_var)
                nucleotide_vars.append(nuc_var)

        return nucleotide_vars

    # ==================================================================
    # Helper: GC content constraint
    # ==================================================================

    def _add_gc_constraint(
        self,
        cp_model: Any,
        codon_vars: list[Any],
        codon_domains: list[list[str]],
        n_codons: int,
    ) -> Any | None:
        """Add GC content range constraint to the model.

        For each codon variable, compute the GC count (0-3 G/C bases) as
        auxiliary variables. The total GC count must lie within
        [gc_lo * 3n, gc_hi * 3n].

        Uses element constraints to link codon index to GC count.

        Args:
            cp_model: OR-Tools CpModel instance.
            codon_vars: List of codon IntVar variables.
            codon_domains: Parallel list of codon string lists.
            n_codons: Number of codon positions.

        Returns:
            The total GC count IntVar, or None if no constraint added.
        """
        n_bases = n_codons * 3
        gc_lo = self.config.gc_lo
        gc_hi = self.config.gc_hi

        gc_count_vars: list = []
        for i, (cvar, domain) in enumerate(zip(codon_vars, codon_domains)):
            # GC count for each possible codon assignment
            gc_values = [
                sum(1 for ch in codon if ch in "GC")
                for codon in domain
            ]
            gc_var = cp_model.NewIntVar(0, 3, f"gc_{i}")
            cp_model.AddElement(cvar, gc_values, gc_var)
            gc_count_vars.append(gc_var)

        total_gc = cp_model.NewIntVar(0, n_bases, "total_gc")
        cp_model.Add(sum(gc_count_vars) == total_gc)

        # Range constraint
        lo = int(math.ceil(gc_lo * n_bases))
        hi = int(math.floor(gc_hi * n_bases))
        cp_model.Add(total_gc >= lo)
        cp_model.Add(total_gc <= hi)

        return total_gc

    # ==================================================================
    # Helper: Composite automaton constraint (restriction + ATTTA + T-run)
    # ==================================================================

    def _add_composite_automaton_constraint(
        self,
        cp_model: Any,
        nucleotide_vars: list[Any],
        n_codons: int,
    ) -> int:
        """Add a single composite automaton for all forbidden nucleotide patterns.

        Combines restriction sites, ATTTA motif, and T-run constraints into
        ONE automaton constraint for better solver performance. Multiple
        automaton constraints on the same variable sequence can cause
        infeasibility; a single combined constraint avoids this.

        The strategy is:
        1. Build a composite DFA for restriction sites (Aho-Corasick)
        2. Build a T-run DFA
        3. Use the composite DFA as the primary constraint and add
           T-run avoidance as a separate lightweight constraint

        Args:
            cp_model: CpModel instance.
            nucleotide_vars: List of nucleotide IntVar variables.
            n_codons: Number of codon positions.

        Returns:
            Number of constraints added.
        """
        num_constraints = 0

        # ── Collect all forbidden patterns ────────────────────────────
        all_patterns: list[str] = []

        # Restriction sites
        sites = self.config.restriction_sites
        if not sites and self.config.add_default_restriction_sites:
            # Default: avoid common 6+ bp restriction sites.
            # 4 bp sites (AluI, HaeIII, MboI, etc.) are excluded because
            # they appear too frequently to avoid in long sequences and
            # make the model infeasible.
            # Only add defaults when add_default_restriction_sites=True
            # (the default); when restriction_sites=[] is explicitly
            # passed alongside add_default_restriction_sites=False,
            # no restriction sites are added.
            sites = [
                seq for name, seq in RESTRICTION_ENZYMES.items()
                if all(ch in "ACGT" for ch in seq) and len(seq) >= 6
            ]

        for site in sites:
            site_upper = site.upper()
            if all(ch in "ACGT" for ch in site_upper):
                all_patterns.append(site_upper)
                rc = reverse_complement(site_upper)
                if rc != site_upper and all(ch in "ACGT" for ch in rc):
                    all_patterns.append(rc)

        # ATTTA motif and its reverse complement.
        # NOTE: ATTTA avoidance via automaton is disabled by default because
        # it interacts poorly with splice/cpG table constraints for long
        # proteins, often causing infeasibility. ATTTA avoidance is better
        # handled as a post-processing step by the greedy optimizer, which
        # uses context-aware motif detection. Set avoid_attta=True to enable
        # the automaton approach (works well for short proteins < 50 AA).
        if self.config.verbose:
            logger.info(
                "VERBOSE: ATTTA avoidance: avoid_attta=%s, "
                "n_codons=%d, max_codons_for_attta=%d, "
                "will_add_automaton=%s",
                self.config.avoid_attta, n_codons,
                _ATTTA_AUTOMATON_MAX_CODONS,
                self.config.avoid_attta and n_codons <= _ATTTA_AUTOMATON_MAX_CODONS,
            )
        if self.config.avoid_attta and n_codons <= _ATTTA_AUTOMATON_MAX_CODONS:
            all_patterns.append(INSTABILITY_MOTIF)
            rc_attta = reverse_complement(INSTABILITY_MOTIF)
            if rc_attta != INSTABILITY_MOTIF and all(ch in "ACGT" for ch in rc_attta):
                all_patterns.append(rc_attta)

        # Build composite DFA for all substring patterns
        if all_patterns:
            if self.config.verbose:
                logger.info(
                    "VERBOSE: Composite automaton: %d forbidden patterns: %s",
                    len(all_patterns), all_patterns[:10],
                )
            dfa_tuple = build_composite_dfa(all_patterns)
            transition_table, accepting_states = dfa_tuple

            if accepting_states and len(transition_table) * _NUCLEOTIDE_ALPHABET_SIZE <= _MAX_AUTOMATON_TRANSITIONS:
                init_state, trans_list, final_states = dfa_to_ortools_format(
                    transition_table, accepting_states
                )

                # For long sequences (>500aa), use decomposed (chunked) automaton
                # encoding to avoid the performance issues of a single large
                # automaton constraint on very long variable sequences.
                use_decomposed = n_codons > _DECOMPOSED_AUTOMATON_CODON_THRESHOLD
                if use_decomposed:
                    max_pat_len = max(len(p) for p in all_patterns)
                    # Overlap must be at least max_pattern_length to guarantee
                    # that any forbidden pattern spanning a chunk boundary is
                    # fully contained within at least one chunk. Using
                    # max_pat_len - 1 can miss patterns that start in the
                    # non-overlap tail of one chunk and extend exactly to the
                    # first character of the next chunk.
                    overlap = max_pat_len
                    logger.info(
                        "Using decomposed automaton encoding for long sequence "
                        "(%d codons > %d threshold): %d patterns, "
                        "chunk_size=%d, overlap=%d",
                        n_codons, _DECOMPOSED_AUTOMATON_CODON_THRESHOLD,
                        len(all_patterns), _AUTOMATON_CHUNK_SIZE, overlap,
                    )
                    num_constraints += self._add_automaton_constraint_decomposed(
                        cp_model, nucleotide_vars, init_state, final_states,
                        trans_list, overlap,
                    )
                else:
                    try:
                        cp_model.AddAutomaton(
                            nucleotide_vars, init_state, final_states, trans_list
                        )
                        num_constraints += 1
                    except Exception as e:
                        logger.warning("Failed to add composite automaton: %s", e)

        # T-run constraint: add as separate automaton
        if self.config.verbose:
            logger.info(
                "VERBOSE: T-run avoidance: avoid_t_runs=%s, "
                "max_t_run=%d, will_add=%s",
                self.config.avoid_t_runs, _MAX_T_RUN_LENGTH,
                self.config.avoid_t_runs,
            )
        # (This is a different type of constraint — it tracks consecutive
        # T count, not substring matching. It does not conflict with the
        # substring DFA.)
        if self.config.avoid_t_runs:
            if self.config.verbose:
                logger.info("VERBOSE: Adding T-run automaton constraint")
            trun_dfa = build_trun_dfa(max_t=_MAX_T_RUN_LENGTH)
            trun_trans, trun_accepting = trun_dfa

            if trun_accepting:
                init_state, trans_list, final_states = dfa_to_ortools_format(
                    trun_trans, trun_accepting
                )
                use_decomposed = n_codons > _DECOMPOSED_AUTOMATON_CODON_THRESHOLD
                if use_decomposed:
                    # T-run avoidance forbids max_t+1 consecutive T's,
                    # so the effective pattern length is max_t+1.
                    # Overlap must be at least that to catch boundary-spanning runs.
                    overlap = _MAX_T_RUN_LENGTH + 1
                    num_constraints += self._add_automaton_constraint_decomposed(
                        cp_model, nucleotide_vars, init_state, final_states,
                        trans_list, overlap,
                    )
                else:
                    try:
                        cp_model.AddAutomaton(
                            nucleotide_vars, init_state, final_states, trans_list
                        )
                        num_constraints += 1
                    except Exception as e:
                        logger.warning("Failed to add T-run automaton: %s", e)

        return num_constraints

    # ==================================================================
    # Helper: Decomposed (chunked) automaton constraint for long sequences
    # ==================================================================

    def _add_automaton_constraint_decomposed(
        self,
        cp_model: Any,
        nucleotide_vars: list[Any],
        init_state: int,
        final_states: list[int],
        trans_list: list[tuple[int, int, int]],
        overlap: int,
    ) -> int:
        """Add automaton constraints in decomposed (chunked) form for long sequences.

        Splits the nucleotide variable sequence into overlapping chunks and
        adds a separate automaton constraint for each chunk. The overlap
        ensures that forbidden patterns spanning chunk boundaries are detected
        by at least one chunk.

        This approach is essential for sequences > 500aa where a single large
        automaton constraint on the entire variable sequence causes poor solver
        performance due to the large propagation domain.

        Args:
            cp_model: CpModel instance.
            nucleotide_vars: List of nucleotide IntVar variables.
            init_state: Initial state of the DFA.
            final_states: List of accepting (final) states.
            trans_list: List of (from_state, symbol, to_state) transitions.
            overlap: Number of overlapping nucleotides between consecutive
                chunks. Should be at least (max_pattern_length - 1) to
                catch patterns that span chunk boundaries.

        Returns:
            Number of automaton constraints added.
        """
        n_nucleotides = len(nucleotide_vars)
        chunk_size = _AUTOMATON_CHUNK_SIZE
        stride = chunk_size - overlap
        num_constraints = 0

        start = 0
        while start < n_nucleotides:
            end = min(start + chunk_size, n_nucleotides)
            chunk_vars = nucleotide_vars[start:end]

            try:
                cp_model.AddAutomaton(
                    chunk_vars, init_state, final_states, trans_list
                )
                num_constraints += 1
            except Exception as e:
                logger.warning(
                    "Failed to add decomposed automaton chunk [%d:%d]: %s",
                    start, end, e,
                )

            if end >= n_nucleotides:
                break
            start += stride

        logger.info(
            "Decomposed automaton: %d chunks for %d nucleotides "
            "(chunk_size=%d, stride=%d, overlap=%d)",
            num_constraints, n_nucleotides, chunk_size, stride, overlap,
        )
        return num_constraints

    # ==================================================================
    # Helper: Splice site constraints
    # ==================================================================

    def _add_splice_constraints(
        self,
        cp_model: Any,
        codon_vars: list[Any],
        codon_domains: list[list[str]],
        protein: str,
    ) -> int:
        """Add cryptic splice site avoidance constraints.

        Within-codon: Forbid codon assignments that create GT (donor) or
        AG (acceptor) dinucleotides within a single codon, **when
        alternatives exist**. Amino acids like Valine that have only
        GT-containing codons are left unmodified (the cross-codon constraint
        may still help).

        Cross-codon: For adjacent codon pairs (i, i+1), add table constraints
        that forbid codon combinations where codon i ends with G and codon
        i+1 starts with T (creating cross-codon GT), or codon i ends with
        A and codon i+1 starts with G (creating cross-codon AG).

        Args:
            cp_model: CpModel instance.
            codon_vars: List of codon IntVar variables.
            codon_domains: Parallel list of codon string lists.
            protein: The amino acid sequence.

        Returns:
            Number of constraints added.
        """
        num_constraints = 0
        threshold = self.config.cryptic_splice_threshold

        if threshold <= 0:
            # Splice avoidance disabled
            return 0

        # Prokaryotic organisms lack spliceosomes; cryptic splice
        # sites are biologically irrelevant.
        if self.config.auto_detect_organism_domain and not self.config.is_eukaryotic:
            logger.info(
                "OR-Tools: Skipping splice constraints for prokaryotic "
                "organism '%s'",
                self.config.organism,
            )
            return 0

        # Within-codon: prefer codons without GT (donor).
        #
        # We focus on GT avoidance because:
        #   1. GT donor sites are stronger cryptic splice signals than AG acceptors
        #   2. AG avoidance is too restrictive — many common AAs (K, R, E, Q, G)
        #      have AG-containing codons and restricting all of them cascades
        #      into infeasibility with cross-codon constraints
        #   3. AG acceptor avoidance is better handled by the greedy optimizer
        #      with MaxEntScan context scoring
        #
        # Valine (V) has NO GT-free codons, so it is left unrestricted.
        # Cysteine (C) has only TGT and TGC; restricting TGT leaves only
        # TGC (ending in C), which creates cross-codon CpG issues. So we
        # also skip C unless we have 2+ GT-free alternatives.
        #
        # CAI-AWARE GT AVOIDANCE: For amino acids where SOME codons lack
        # GT, we only add the restriction if the best GT-free codon has
        # CAI weight >= 50% of the best GT-containing codon's weight.
        # This prevents excessive CAI degradation for amino acids like
        # Glycine (GGT w=1.0 vs GGA w=0.44 in yeast) where the CAI
        # penalty of avoiding GT would be catastrophic.
        cai_data_splice = CODON_ADAPTIVENESS_TABLES.get(
            self.config.organism, {}
        )
        _CAI_GT_AVOIDANCE_RATIO = 0.5

        for i, (cvar, domain) in enumerate(zip(codon_vars, codon_domains)):
            gt_free = [c for c in domain if "GT" not in c]
            gt_containing = [c for c in domain if "GT" in c]

            if not gt_free or not gt_containing:
                continue

            # Check CAI penalty of avoiding GT
            best_gt_free_cai = max(
                cai_data_splice.get(c, 0.01) for c in gt_free
            )
            best_gt_containing_cai = max(
                cai_data_splice.get(c, 0.01) for c in gt_containing
            )

            if best_gt_containing_cai > 0:
                cai_ratio = best_gt_free_cai / best_gt_containing_cai
            else:
                cai_ratio = 1.0

            # Only add GT restriction if the CAI penalty is acceptable
            if cai_ratio < _CAI_GT_AVOIDANCE_RATIO:
                continue

            # Only restrict if GT-free alternatives exist AND there are
            # at least 2 alternatives (to avoid cascading infeasibility
            # with cross-codon constraints)
            if len(gt_free) >= 2 and len(gt_free) < len(domain):
                allowed_indices = [domain.index(c) for c in gt_free]
                cp_model.AddAllowedAssignments(
                    [cvar],
                    [[idx] for idx in allowed_indices],
                )
                num_constraints += 1
            elif gt_free and len(gt_free) < len(domain):
                # Only 1 GT-free alternative — add constraint but check
                # if it would create cross-codon CpG issues (ends with C).
                # If so, skip to avoid infeasibility.
                sole_codon = gt_free[0]
                if sole_codon[-1] != "C":
                    # Safe: the sole codon does not end with C
                    allowed_indices = [domain.index(c) for c in gt_free]
                    cp_model.AddAllowedAssignments(
                        [cvar],
                        [[idx] for idx in allowed_indices],
                    )
                    num_constraints += 1
            # If no GT-free codons (e.g. Valine), leave domain unrestricted

        # Cross-codon: prefer codon pairs that avoid GT/AG at the boundary.
        #
        # A cross-codon GT arises when codon i ends with G and codon i+1
        # starts with T.  A cross-codon AG arises when codon i ends with
        # A and codon i+1 starts with G.  Both create cryptic splice sites.
        #
        # We only add the constraint when viable alternatives exist and at
        # least 25% of pairs remain. Overly restrictive pair constraints
        # cascade and make long sequences infeasible.
        for i in range(len(codon_vars) - 1):
            domain_i = codon_domains[i]
            domain_j = codon_domains[i + 1]
            total_pairs = len(domain_i) * len(domain_j)

            # Build allowed pairs: exclude (last=G, first=T) and (last=A, first=G)
            allowed_pairs = []
            for ki, codon_i in enumerate(domain_i):
                for kj, codon_j in enumerate(domain_j):
                    boundary = codon_i[-1] + codon_j[0]
                    if boundary != "GT" and boundary != "AG":
                        allowed_pairs.append([ki, kj])

            n_allowed = len(allowed_pairs)
            if 0 < n_allowed < total_pairs:
                # Only add if at least 25% of pairs remain viable
                min_required = max(1, total_pairs // _MIN_VIABLE_PAIR_DENOMINATOR)
                if n_allowed >= min_required:
                    cp_model.AddAllowedAssignments(
                        [codon_vars[i], codon_vars[i + 1]],
                        allowed_pairs,
                    )
                    num_constraints += 1

        return num_constraints

    # ==================================================================
    # Helper: CpG avoidance constraints
    # ==================================================================

    def _add_cpg_constraints(
        self,
        cp_model: Any,
        codon_vars: list[Any],
        codon_domains: list[list[str]],
        protein: str,
    ) -> int:
        """Add CpG dinucleotide avoidance constraints.

        Within-codon: Forbid codon assignments containing the CG dinucleotide
        if alternatives exist.

        Cross-codon: For adjacent codon pairs (i, i+1), forbid combinations
        where codon i ends with C and codon i+1 starts with G (creating a
        cross-codon CG dinucleotide).

        Args:
            cp_model: CpModel instance.
            codon_vars: List of codon IntVar variables.
            codon_domains: Parallel list of codon string lists.
            protein: The amino acid sequence.

        Returns:
            Number of constraints added.
        """
        if not self.config.avoid_cpg:
            return 0

        # CpG avoidance is primarily relevant for eukaryotic organisms
        # (DNA methylation silencing). Skip for prokaryotes.
        if self.config.auto_detect_organism_domain and not self.config.is_eukaryotic:
            logger.info(
                "OR-Tools: Skipping CpG constraints for prokaryotic "
                "organism '%s'",
                self.config.organism,
            )
            return 0

        num_constraints = 0

        # Within-codon: forbid codons containing CG
        for i, (cvar, domain) in enumerate(zip(codon_vars, codon_domains)):
            cg_free = [c for c in domain if "CG" not in c]
            if cg_free and len(cg_free) < len(domain):
                # Some codons contain CG — restrict to CG-free alternatives
                allowed_indices = [domain.index(c) for c in cg_free]
                cp_model.AddAllowedAssignments(
                    [cvar],
                    [[idx] for idx in allowed_indices],
                )
                num_constraints += 1

        # Cross-codon: forbid (C-ending codon, G-starting codon) pairs.
        # We only add the constraint when alternatives exist and at least
        # 25% of pairs remain viable, to avoid cascading infeasibility.
        for i in range(len(codon_vars) - 1):
            domain_i = codon_domains[i]
            domain_j = codon_domains[i + 1]
            total_pairs = len(domain_i) * len(domain_j)

            # Build allowed pairs: exclude (last=C, first=G)
            allowed_pairs = []
            for ki, codon_i in enumerate(domain_i):
                for kj, codon_j in enumerate(domain_j):
                    if not (codon_i[-1] == "C" and codon_j[0] == "G"):
                        allowed_pairs.append([ki, kj])

            n_allowed = len(allowed_pairs)
            if 0 < n_allowed < total_pairs:
                # Only add if at least 25% of pairs remain viable
                min_required = max(1, total_pairs // _MIN_VIABLE_PAIR_DENOMINATOR)
                if n_allowed >= min_required:
                    cp_model.AddAllowedAssignments(
                        [codon_vars[i], codon_vars[i + 1]],
                        allowed_pairs,
                    )
                    num_constraints += 1

        return num_constraints

    # ==================================================================
    # Helper: CAI objective
    # ==================================================================

    def _set_cai_objective(
        self,
        cp_model: Any,
        codon_vars: list[Any],
        codon_domains: list[list[str]],
        organism: str,
        cai_scale: int = _CAI_SCALE,
    ) -> None:
        """Set the CAI maximization objective on the CP-SAT model.

        For each codon position, computes the log(CAI) value for each possible
        codon assignment, scales to integers (×cai_scale), and creates
        auxiliary variables linked via element constraints. The objective is
        to maximize the sum of all CAI contribution terms.

        The CAI for a codon c encoding amino acid aa is defined as:
            w(c) = freq(c) / max_freq(aa)
        where freq is the per-thousand usage frequency. The objective
        maximizes Σ log(w(c_i)).

        The maximization objective is ALWAYS set, even when hard constraints
        dominate, to ensure the solver finds the FEASIBLE solution with
        MAXIMUM CAI rather than an arbitrary feasible solution.

        Args:
            cp_model: CpModel instance.
            codon_vars: List of codon IntVar variables.
            codon_domains: Parallel list of codon string lists.
            organism: Target organism name (e.g. "Homo_sapiens").
            cai_scale: Scaling factor for converting float log-CAI to integer
                objective. Use a smaller value (_CAI_SCALE_LONG=1000) for long
                proteins (>1000aa) to avoid integer overflow.
        """
        adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})
        cai_terms: list = []

        for i, (cvar, domain) in enumerate(zip(codon_vars, codon_domains)):
            # Compute scaled log-CAI for each codon in the domain
            cai_values = []
            for codon in domain:
                w = adaptiveness.get(codon, _DEFAULT_CAI_WEIGHT)  # Small default for missing data
                if w <= 0:
                    w = _DEFAULT_CAI_WEIGHT  # Avoid log(0)
                scaled = int(round(math.log(w) * cai_scale))
                cai_values.append(scaled)

            # Create auxiliary variable and link via element constraint
            min_val = min(cai_values)
            max_val = max(cai_values)

            if min_val == max_val:
                # All codons have the same CAI — no choice effect
                cai_var = cp_model.NewIntVar(min_val, max_val, f"cai_{i}")
                cp_model.Add(cai_var == min_val)
            else:
                cai_var = cp_model.NewIntVar(min_val, max_val, f"cai_{i}")
                cp_model.AddElement(cvar, cai_values, cai_var)

            cai_terms.append(cai_var)

        # Always set the CAI maximization objective to ensure the solver
        # seeks the FEASIBLE solution with MAXIMUM CAI, even when hard
        # constraints dominate. Without an objective, the solver may return
        # any arbitrary feasible solution with suboptimal CAI.
        if cai_terms:
            cp_model.Maximize(sum(cai_terms))
        else:
            # Edge case: no codon positions — set a trivial objective
            cp_model.Maximize(0)

    # ==================================================================
    # Helper: Solution extraction
    # ==================================================================

    def _extract_solution(
        self,
        solver: Any,
        codon_vars: list[Any],
        codon_domains: list[list[str]],
        protein: str,
    ) -> str:
        """Map solver variable values back to codon strings → DNA sequence.

        Args:
            solver: The CP-SAT solver with a feasible solution.
            codon_vars: List of codon IntVar variables.
            codon_domains: Parallel list of codon string lists.
            protein: The amino acid sequence.

        Returns:
            The optimized DNA sequence string.
        """
        codons = []
        for i, (cvar, domain) in enumerate(zip(codon_vars, codon_domains)):
            idx = solver.Value(cvar)
            # Clamp index to valid range (defense-in-depth)
            idx = max(0, min(idx, len(domain) - 1))
            codon = domain[idx]
            codons.append(codon)
        sequence = "".join(codons)

        logger.info(
            "OR-Tools solution extraction: %d codons extracted, "
            "sequence_len=%d expected_len=%d",
            len(codons), len(sequence), len(protein) * 3,
        )

        # Validate: non-empty sequence
        if not sequence:
            logger.error("OR-Tools: Solution extraction produced empty sequence")
            return ""

        # Validate: correct length
        expected_len = len(protein) * 3
        if len(sequence) != expected_len:
            logger.error(
                "OR-Tools: Sequence length (%d) != expected (%d)",
                len(sequence), expected_len,
            )

        # Validate: translation fidelity
        from biocompiler.shared.constants import CODON_TABLE
        for i, expected_aa in enumerate(protein):
            codon = sequence[i * 3 : i * 3 + 3]
            actual_aa = CODON_TABLE.get(codon)
            if actual_aa != expected_aa:
                logger.warning(
                    "OR-Tools: Codon %d (%s) translates to '%s', expected '%s'",
                    i, codon, actual_aa, expected_aa,
                )

        return sequence

    # ==================================================================
    # Helper: GC computation
    # ==================================================================

    @staticmethod
    def _compute_gc(sequence: str) -> float:
        """Compute GC content fraction of a DNA sequence.

        Args:
            sequence: DNA sequence string.

        Returns:
            GC fraction in [0.0, 1.0].
        """
        if not sequence:
            return 0.0
        gc = sum(1 for ch in sequence if ch in "GC")
        return round(gc / len(sequence), 4)

    # ==================================================================
    # Helper: CAI computation
    # ==================================================================

    @staticmethod
    def _compute_cai(sequence: str, protein: str, organism: str) -> float:
        """Compute the Codon Adaptation Index for a DNA sequence.

        CAI = geometric mean of relative adaptiveness values:
            CAI = (∏ w_i)^(1/N)
        where w_i is the relative adaptiveness of the codon at position i,
        and N is the number of codons.

        Args:
            sequence: DNA sequence string.
            protein: Amino acid sequence.
            organism: Target organism name.

        Returns:
            CAI value in [0.0, 1.0].
        """
        adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})
        if not sequence or not protein:
            return 0.0

        log_sum = 0.0
        n = 0
        for i, aa in enumerate(protein):
            codon = sequence[i * 3 : i * 3 + 3]
            w = adaptiveness.get(codon, _DEFAULT_CAI_WEIGHT)
            if w <= 0:
                w = _DEFAULT_CAI_WEIGHT
            log_sum += math.log(w)
            n += 1

        if n == 0:
            return 0.0

        cai = math.exp(log_sum / n)
        return round(min(max(cai, 0.0), 1.0), 4)

    # ==================================================================
    # Helper: MUS (Minimal Unsatisfiable Subset) computation
    # ==================================================================

    def _compute_mus(
        self,
        protein: str,
        organism: str,
        codon_domains: list[list[str]],
    ) -> MUSReport:
        """Compute a Minimal Unsatisfiable Subset for diagnosis.

        When the solver reports INFEASIBLE, this method attempts to identify
        the smallest set of conflicting constraints. Uses a heuristic
        approach that checks which constraint categories are likely
        conflicting based on static analysis of the protein sequence.

        Since full MUS computation (deletion-based) requires re-solving
        the model multiple times, we use a fast heuristic that checks:
        1. GC range feasibility
        2. Splice site vs. codon domain conflicts
        3. CpG vs. codon domain conflicts
        4. Cross-codon pair elimination

        Args:
            protein: Amino acid sequence.
            organism: Target organism.
            codon_domains: Codon domain lists.

        Returns:
            MUSReport with identified conflicts and suggested relaxations.
        """
        conflicting: list[str] = []
        suggestions: list[str] = []

        # Check GC range feasibility
        gc_lo = self.config.gc_lo
        gc_hi = self.config.gc_hi
        n_bases = len(protein) * 3

        # Compute min/max possible GC count
        min_gc = 0
        max_gc = 0
        for domain in codon_domains:
            gc_counts = [sum(1 for ch in c if ch in "GC") for c in domain]
            min_gc += min(gc_counts)
            max_gc += max(gc_counts)

        min_gc_frac = min_gc / n_bases if n_bases > 0 else 0
        max_gc_frac = max_gc / n_bases if n_bases > 0 else 1

        if gc_lo > max_gc_frac or gc_hi < min_gc_frac:
            conflicting.append("gc_content")
            suggestions.append(
                f"Widen GC range: achievable [{min_gc_frac:.2f}, {max_gc_frac:.2f}], "
                f"requested [{gc_lo:.2f}, {gc_hi:.2f}]"
            )

        # Check splice constraint feasibility
        if self.config.cryptic_splice_threshold > 0:
            splice_conflict = False
            for i, aa in enumerate(protein):
                domain = codon_domains[i]
                gt_free = [c for c in domain if "GT" not in c]
                ag_free = [c for c in domain if "AG" not in c]
                both_free = [c for c in domain if "GT" not in c and "AG" not in c]
                if not both_free and (not gt_free or not ag_free):
                    # This AA has unavoidable GT or AG
                    splice_conflict = True
                    break

            if splice_conflict:
                conflicting.append("cryptic_splice")
                suggestions.append(
                    "Some amino acids have no GT/AG-free codons; "
                    "consider reducing cryptic_splice_threshold or "
                    "accepting splice sites for Valine positions"
                )

            # Check cross-codon splice pairs
            for i in range(len(protein) - 1):
                domain_i = codon_domains[i]
                domain_j = codon_domains[i + 1]
                allowed = 0
                for ci in domain_i:
                    for cj in domain_j:
                        boundary = ci[-1] + cj[0]
                        if boundary != "GT" and boundary != "AG":
                            allowed += 1
                if allowed == 0:
                    conflicting.append(f"cross_splice_{i}_{i+1}")
                    suggestions.append(
                        f"Cross-codon splice conflict at positions {i}-{i+1} "
                        f"({protein[i]}-{protein[i+1]}): all codon pairs "
                        f"create GT or AG at boundary"
                    )

        # Check CpG constraint feasibility
        if self.config.avoid_cpg:
            for i in range(len(protein) - 1):
                domain_i = codon_domains[i]
                domain_j = codon_domains[i + 1]
                allowed = 0
                for ci in domain_i:
                    for cj in domain_j:
                        if not (ci[-1] == "C" and cj[0] == "G"):
                            allowed += 1
                if allowed == 0:
                    conflicting.append(f"cross_cpg_{i}_{i+1}")
                    suggestions.append(
                        f"Cross-codon CpG conflict at positions {i}-{i+1} "
                        f"({protein[i]}-{protein[i+1]}): consider "
                        f"disabling CpG avoidance"
                    )

        # If no specific conflicts found, report generic
        if not conflicting:
            conflicting = ["constraint_combination"]
            suggestions = [
                "Try relaxing one or more constraints: "
                "widen GC range, remove some restriction sites, "
                "or disable CpG/ATTTA avoidance"
            ]

        return MUSReport(
            conflicting_constraints=conflicting,
            explanation=(
                f"The constraint set is infeasible for protein of length "
                f"{len(protein)}. Conflicting constraints: "
                f"{', '.join(conflicting)}"
            ),
            suggested_relaxations=suggestions,
        )


# ======================================================================
# CSPModel convenience constructor
# ======================================================================

@dataclass
class CSPModel:
    """Constraint Satisfaction Problem model for gene optimization.

    Encodes the gene design problem as a CSP: given a protein sequence
    and a set of biological constraints, find a DNA sequence that
    encodes the protein while satisfying all hard constraints and
    maximizing the CAI objective.

    Attributes:
        protein: Amino acid sequence (string of single-letter codes).
        organism: Target organism for codon optimization.
        codon_variables: Pre-built codon decision variables (optional;
            the engine builds them if not provided).
    """

    protein: str
    organism: str = "Homo_sapiens"
    codon_variables: Optional[list[CodonVariable]] = None

    def __post_init__(self) -> None:
        """Validate the CSP model."""
        if not self.protein:
            raise ValueError("Protein sequence must be non-empty")
        if self.organism not in SUPPORTED_ORGANISMS:
            logger.warning(
                "Organism '%s' not in supported list %s; "
                "CAI values may be unavailable",
                self.organism, SUPPORTED_ORGANISMS,
            )


# ======================================================================
# Greedy fallback (always available)
# ======================================================================

def greedy_fallback(
    protein: str,
    organism: str,
    config: SolverConfig,
) -> SolverResult:
    """Simple greedy codon optimization fallback.

    Selects the highest-CAI codon for each amino acid position without
    global constraint solving. Used when the CP-SAT solver is unavailable
    or times out.

    This is NOT the main greedy optimizer (which is in optimization.py).
    This is a minimal fallback for the solver module.

    Args:
        protein: Amino acid sequence.
        organism: Target organism.
        config: Solver configuration.

    Returns:
        SolverResult with greedy-optimized sequence.
    """
    start_time = time.monotonic()
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})

    codons = []
    for aa in protein:
        domain = AA_TO_CODONS.get(aa, [])
        if not domain:
            codons.append("NNN")
            continue
        # Pick highest-CAI codon
        best = max(domain, key=lambda c: adaptiveness.get(c, 0.0))
        codons.append(best)

    sequence = "".join(codons)
    gc_val = ORTOOLSEngine._compute_gc(sequence)
    cai_val = ORTOOLSEngine._compute_cai(sequence, protein, organism)

    return SolverResult(
        sequence=sequence,
        solved=True,
        backend_used=SolverBackend.GREEDY_FALLBACK,
        cai=cai_val,
        gc_content=gc_val,
        solve_time_seconds=time.monotonic() - start_time,
        fallback_used=True,
        warnings=["Used greedy fallback (no constraint satisfaction)"],
    )
