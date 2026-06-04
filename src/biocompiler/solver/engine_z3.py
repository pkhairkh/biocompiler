"""
BioCompiler Z3 SMT Solver Engine
=================================

Fallback solver backend using Microsoft Z3 SMT solver.

Z3 excels at continuous/arithmetic constraints (MaxEntScan splice scores,
GC content optimization) but is slower than OR-Tools for automaton-style
pattern constraints (restriction site detection). Use this engine when
OR-Tools is unavailable or when the problem is dominated by arithmetic
constraints over the MaxEntScan scoring model.

Architecture
------------
1. **Variables**: One Z3 ``Int`` per codon position, domain [0, N) where
   N = number of synonymous codons for that amino acid.
2. **Nucleotide encoding**: Each codon variable is expanded into three
   Z3 integer expressions (one per base) using nested ``If`` chains.
3. **Constraints**:
   - GC content: sum of GC indicators ∈ [gc_lo·3n, gc_hi·3n]
   - Restriction sites: no window of |pattern| nucleotides equals the
     forbidden pattern (encoded as ``Not(And(...))`` per window)
   - Splice sites: If-Then-Else encoding of MaxEntScan thresholds
4. **Objective**: MAXIMIZE Σ log(CAI(codon_i)), scaled to integers
   for Z3 optimization (×10000 for 4-decimal precision).
5. **UNSAT core**: When the model is infeasible, extract an UNSAT core
   from Z3 and map it back to named constraints for MUS reporting.

Performance Notes
-----------------
- For proteins > 500 aa, Z3 may time out on restriction-site constraints.
  Consider relaxing the forbidden pattern list or increasing timeout.
- Splice site encoding is O(L) where L = sequence length, with a
  constant factor of ~9 If-Then-Else chains per GT position.

References
----------
- Yeo & Burge (2004) "Maximum entropy modeling of short sequence motifs
  with applications to RNA splicing", J Comput Biol 11(2-3):377-94.
- de Moura & Bjørner (2008) "Z3: An Efficient SMT Solver", TACAS.
"""

from __future__ import annotations

import math
import logging
import time
from typing import Optional

# ---------------------------------------------------------------------------
# Z3 import with graceful fallback
# ---------------------------------------------------------------------------
try:
    from z3 import (
        Optimize,
        Int,
        IntVal,
        If,
        And,
        Or,
        Not,
        Implies,
        Sum,
        ArithRef,
        BoolRef,
        sat,
        unsat,
        unknown,
    )
    _Z3_AVAILABLE = True
except ImportError:
    _Z3_AVAILABLE = False

from .types import (
    SolverConfig,
    SolverResult,
    SolverBackend,
    SolverBackendProtocol,
    CSPModel,
    ConstraintStrictness,
    MUSReport,
    ConstraintViolation,
)
from ..constants import AA_TO_CODONS, BASE_REV, IUPAC_EXPAND
from ..organisms import CODON_ADAPTIVENESS_TABLES

logger = logging.getLogger(__name__)

__all__ = ["Z3Engine"]

# Scaling factor for converting float CAI weights to Z3 integers
_CAI_SCALE: int = 10000  # 4 decimal places of precision
# Scale factor for MaxEntScan scores to convert to Z3-compatible integers
_SCORE_SCALE: int = 100  # 2 decimal places


def _precompute_donor_scores() -> list[list[float]]:
    """Pre-compute MaxEntScan donor PWM log-odds contributions.

    Returns a 9×4 matrix where entry [pwm_pos][base_idx] is the
    log-odds contribution of that base at that PWM position.
    Base index: A=0, C=1, G=2, T=3.
    """
    from ..maxentscan import DONOR_PWM_SCORE, BG_PROB
    _EPSILON = 0.001
    scores: list[list[float]] = []
    for pwm_pos in range(9):
        row: list[float] = []
        for base_idx in range(4):
            prob = max(DONOR_PWM_SCORE[pwm_pos][base_idx], _EPSILON)
            row.append(math.log2(prob / BG_PROB))
        scores.append(row)
    return scores


def _precompute_acceptor_scores() -> list[list[float]]:
    """Pre-compute MaxEntScan acceptor PWM log-odds contributions.

    Returns a 23×4 matrix where entry [pwm_pos][base_idx] is the
    log-odds contribution of that base at that PWM position.
    Base index: A=0, C=1, G=2, T=3.
    """
    from ..maxentscan import ACCEPTOR_PWM_SCORE, BG_PROB
    _EPSILON = 0.001
    scores: list[list[float]] = []
    for pwm_pos in range(23):
        row: list[float] = []
        for base_idx in range(4):
            prob = max(ACCEPTOR_PWM_SCORE[pwm_pos][base_idx], _EPSILON)
            row.append(math.log2(prob / BG_PROB))
        scores.append(row)
    return scores


# Pre-compute once at module load (cheap — just 9×4 + 23×4 floats)
_DONOR_SCORES = _precompute_donor_scores() if _Z3_AVAILABLE else []
_ACCEPTOR_SCORES = _precompute_acceptor_scores() if _Z3_AVAILABLE else []


class Z3Engine(SolverBackendProtocol):
    """Z3 SMT solver backend for BioCompiler codon optimization.

    This is the FALLBACK solver engine — use when OR-Tools is unavailable.
    Z3 is better for continuous constraints (MaxEntScan scores) but slower
    for automaton-style pattern constraints (restriction sites).

    Implements the :class:`SolverBackendProtocol` — any module that
    depends on a solver backend can use a Z3Engine instance via the
    protocol interface.

    Usage::

        config = SolverConfig(gc_lo=0.40, gc_hi=0.60)
        engine = Z3Engine(config, organism="Homo_sapiens")
        model = CSPModel(
            protein_sequence="MVLSPADKTNVKAAWGKVGA",
            codon_domains={0: ["ATG"], 1: ["GTT", "GTC", "GTA", "GTG"], ...},
            constraints=[],
        )
        result = engine.solve(model)
        if result.solved:
            print(result.sequence)

    Thread safety: NOT thread-safe. Create one engine per thread.
    """

    def __init__(
        self,
        config: SolverConfig,
        organism: str = "Homo_sapiens",
        seed: int = 0,
    ) -> None:
        """Initialize the Z3 solver engine.

        Args:
            config: Solver configuration parameters.
            organism: Target organism for codon usage tables.
                Defaults to ``"Homo_sapiens"``.
            seed: Deterministic random seed for the Z3 optimizer.
                Defaults to ``0``.  Use the same seed for reproducible
                results across runs.

        Raises:
            ImportError: If z3-solver is not installed (check
                :meth:`is_available` first to avoid this).
        """
        if not _Z3_AVAILABLE:
            raise ImportError(
                "z3-solver package is not installed. "
                "Install with: pip install z3-solver"
            )
        self.config = config
        self.organism = organism
        self._seed: int = seed
        self._optimizer: Optional[Optimize] = None

    # -------------------------------------------------------------------
    # Availability check
    # -------------------------------------------------------------------

    @staticmethod
    def is_available() -> bool:
        """Check whether the z3-solver package is importable.

        Returns:
            True if Z3 can be imported, False otherwise.
        """
        return _Z3_AVAILABLE

    # -------------------------------------------------------------------
    # Main solve method
    # -------------------------------------------------------------------

    def solve(self, model: CSPModel) -> SolverResult:
        """Solve the codon optimization CSP using Z3.

        Creates Z3 integer variables for each codon position, encodes
        all constraints (GC, restriction sites, splice sites), sets
        the CAI maximization objective, and invokes the Z3 optimizer.

        Args:
            model: The CSP model encoding the optimization problem.

        Returns:
            SolverResult with the optimized sequence on success, or
            an UNSAT/timeout/error report on failure.
        """
        start_time = time.perf_counter()
        config = model.config
        protein = model.protein_sequence.upper()

        if not protein:
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                solve_time_seconds=time.perf_counter() - start_time,
                warnings=["Empty protein sequence"],
            )

        # Get codon usage data for the organism
        cai_data = CODON_ADAPTIVENESS_TABLES.get(self.organism)
        if cai_data is None:
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                solve_time_seconds=time.perf_counter() - start_time,
                warnings=[f"Unknown organism: {self.organism}"],
            )

        n = len(protein)

        # Create optimizer
        optimizer = Optimize()
        optimizer.set("timeout", int(config.timeout_seconds * 1000))
        # Note: Z3's Optimize does not support 'random_seed' (only Solver does).
        # We use 'rlimit' as a resource limit proxy for reproducibility.
        try:
            optimizer.set("random_seed", self._seed)
        except Exception:
            # Not all Z3 builds support random_seed on Optimize; ignore.
            logger.debug("Z3 Optimize does not support random_seed, skipping")

        # --- Create codon variables and domains ---
        codon_vars: list[ArithRef] = []
        codon_lists: list[list[str]] = []

        for i, aa in enumerate(protein):
            # Use codon_domains if provided, else fall back to AA_TO_CODONS
            if model.codon_domains and i in model.codon_domains:
                codons_for_aa = model.codon_domains[i]
            else:
                codons_for_aa = AA_TO_CODONS.get(aa, [])

            if not codons_for_aa:
                return SolverResult(
                    sequence="",
                    solved=False,
                    backend_used=SolverBackend.Z3,
                    solve_time_seconds=time.perf_counter() - start_time,
                    warnings=[
                        f"No codons found for amino acid '{aa}' at position {i}"
                    ],
                )
            c_var = Int(f"codon_{i}")
            optimizer.add(c_var >= 0)
            optimizer.add(c_var < len(codons_for_aa))
            codon_vars.append(c_var)
            codon_lists.append(codons_for_aa)

        # --- Create nucleotide expressions ---
        nuc_exprs: list[ArithRef] = []
        for c_var, codons_for_aa in zip(codon_vars, codon_lists):
            nuc0, nuc1, nuc2 = self._encode_codon_as_nucleotides(
                c_var, codons_for_aa
            )
            nuc_exprs.extend([nuc0, nuc1, nuc2])

        total_nucs = len(nuc_exprs)  # = 3 * n

        # --- Add constraints ---
        constraint_names: list[str] = []

        # Track number of domain assertions for UNSAT core mapping
        # Each codon position adds 2 domain assertions (>=0, <N)
        n_domain_assertions = 2 * n

        # GC constraint
        constraint_names.extend(
            self._encode_gc_constraint(optimizer, nuc_exprs, protein, config)
        )

        # Restriction site constraints
        if config.restriction_sites:
            constraint_names.extend(
                self._encode_restriction_site_constraints(
                    optimizer, nuc_exprs, config.restriction_sites
                )
            )

        # Splice site constraints
        constraint_names.extend(
            self._encode_splice_constraints(optimizer, nuc_exprs, protein, config)
        )

        # --- Set optimization objective: MAXIMIZE sum of log(CAI) ---
        self._encode_cai_objective(optimizer, codon_vars, codon_lists)

        # --- Solve ---
        if config.verbose:
            logger.info(
                "Z3: Solving codon optimization for %d aa protein "
                "(%d nucleotide variables, %d constraints)",
                n, total_nucs, len(constraint_names),
            )

        result_status = optimizer.check()
        solve_time = time.perf_counter() - start_time

        if result_status == sat:
            z3_model = optimizer.model()
            return self._extract_solution(
                z3_model, codon_vars, codon_lists, protein,
                solve_time, len(constraint_names), n
            )
        elif result_status == unsat:
            mus = self._extract_unsat_core(optimizer, constraint_names, n_domain_assertions)
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                solve_time_seconds=solve_time,
                num_constraints=len(constraint_names),
                num_variables=n,
                mus_report=mus,
                warnings=["Model is unsatisfiable — conflicting constraints"],
            )
        elif result_status == unknown:
            reason = optimizer.reason_unknown()
            is_timeout = "timeout" in reason.lower() or "canceled" in reason.lower()
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                solve_time_seconds=solve_time,
                num_constraints=len(constraint_names),
                num_variables=n,
                warnings=[f"Solver returned unknown: {reason}"],
                fallback_used=is_timeout,
            )
        else:
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                solve_time_seconds=solve_time,
                warnings=[f"Unexpected Z3 result: {result_status}"],
            )

    # -------------------------------------------------------------------
    # Nucleotide encoding
    # -------------------------------------------------------------------

    def _encode_codon_as_nucleotides(
        self,
        z3_codon_var: ArithRef,
        codons_list: list[str],
    ) -> tuple[ArithRef, ArithRef, ArithRef]:
        """Encode a codon variable as three Z3 nucleotide expressions.

        For each of the three nucleotide positions within a codon, builds
        a nested If-Then-Else chain that maps the codon index to the
        corresponding base integer (A=0, C=1, G=2, T=3).

        Example for alanine (GCT, GCC, GCA, GCG), position 0::

            If(codon==0, 2, If(codon==1, 2, If(codon==2, 2, 2)))
            # All codons start with G → constant 2

        Args:
            z3_codon_var: Z3 Int variable for the codon choice.
            codons_list: List of synonymous codon strings for this AA.

        Returns:
            Tuple of three Z3 expressions (nuc0, nuc1, nuc2), each
            representing the base at that position as an integer.
        """
        n_codons = len(codons_list)

        def _build_ite(base_pos: int) -> ArithRef:
            """Build If-Then-Else chain for one nucleotide position."""
            if n_codons == 1:
                return IntVal(BASE_REV[codons_list[0][base_pos]])
            # Build from the last codon backwards
            expr: ArithRef = IntVal(BASE_REV[codons_list[-1][base_pos]])
            for k in range(n_codons - 2, -1, -1):
                expr = If(
                    z3_codon_var == k,
                    IntVal(BASE_REV[codons_list[k][base_pos]]),
                    expr,
                )
            return expr

        return _build_ite(0), _build_ite(1), _build_ite(2)

    # -------------------------------------------------------------------
    # GC constraint encoding
    # -------------------------------------------------------------------

    def _encode_gc_constraint(
        self,
        optimizer: Optimize,
        nuc_exprs: list[ArithRef],
        protein: str,
        config: SolverConfig,
    ) -> list[str]:
        """Encode GC content constraints.

        GC count = Σ If(nuc_i ∈ {C, G}, 1, 0)
        Constraint: gc_lo * 3n ≤ gc_count ≤ gc_hi * 3n

        Args:
            optimizer: Z3 Optimize instance to add constraints to.
            nuc_exprs: Z3 expressions for each nucleotide position.
            protein: Protein sequence.
            config: Solver configuration with GC bounds.

        Returns:
            List of constraint names added (for UNSAT core tracking).
        """
        total_nucs = 3 * len(protein)
        constraint_names: list[str] = []

        gc_indicators: list[ArithRef] = []
        for nuc in nuc_exprs:
            gc_indicators.append(If(Or(nuc == 1, nuc == 2), IntVal(1), IntVal(0)))

        gc_count = Sum(gc_indicators)

        gc_lo_bound = math.ceil(config.gc_lo * total_nucs)
        gc_hi_bound = math.floor(config.gc_hi * total_nucs)

        optimizer.add(gc_count >= gc_lo_bound)
        optimizer.add(gc_count <= gc_hi_bound)

        constraint_names.append(f"gc_lo_{config.gc_lo}")
        constraint_names.append(f"gc_hi_{config.gc_hi}")

        if config.verbose:
            logger.info(
                "Z3: GC constraint: %d <= gc_count <= %d (of %d nucs)",
                gc_lo_bound, gc_hi_bound, total_nucs,
            )

        return constraint_names

    # -------------------------------------------------------------------
    # Restriction site constraint encoding
    # -------------------------------------------------------------------

    def _encode_restriction_site_constraints(
        self,
        optimizer: Optimize,
        nuc_exprs: list[ArithRef],
        forbidden_patterns: list[str],
    ) -> list[str]:
        """Encode restriction site avoidance constraints.

        For each forbidden pattern of length k, adds constraints that
        no consecutive window of k nucleotides in the sequence equals
        the pattern. Each window is encoded as::

            Not(And(nuc_j == p[0], nuc_{j+1} == p[1], ...,
                    nuc_{j+k-1} == p[k-1]))

        Handles IUPAC ambiguity codes by expanding them to sets of bases.

        Args:
            optimizer: Z3 Optimize instance.
            nuc_exprs: Z3 expressions for each nucleotide position.
            forbidden_patterns: List of recognition sequences to avoid.

        Returns:
            List of constraint names added.
        """
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)

        for pattern in forbidden_patterns:
            pattern_upper = pattern.upper()
            pat_len = len(pattern_upper)
            if pat_len == 0 or pat_len > total_nucs:
                continue

            # Expand IUPAC codes in the pattern
            expanded_bases: list[list[int]] = []
            for base_char in pattern_upper:
                iupac_bases = IUPAC_EXPAND.get(base_char, base_char)
                base_ints = [BASE_REV[b] for b in iupac_bases if b in BASE_REV]
                if not base_ints:
                    break
                expanded_bases.append(base_ints)

            if len(expanded_bases) != pat_len:
                logger.warning(
                    "Z3: Skipping restriction pattern '%s' — "
                    "contains unresolvable IUPAC codes", pattern
                )
                continue

            # For each window position
            for j in range(total_nucs - pat_len + 1):
                base_eqs: list[BoolRef] = []
                for k in range(pat_len):
                    allowed = expanded_bases[k]
                    if len(allowed) == 1:
                        base_eqs.append(nuc_exprs[j + k] == allowed[0])
                    else:
                        base_eqs.append(
                            Or(*[nuc_exprs[j + k] == b for b in allowed])
                        )
                optimizer.add(Not(And(*base_eqs)))

            constraint_names.append(f"no_{pattern}")

            if self.config.verbose:
                logger.info(
                    "Z3: Added restriction site constraint for '%s' "
                    "(%d windows)", pattern, total_nucs - pat_len + 1
                )

        return constraint_names

    # -------------------------------------------------------------------
    # Splice site constraint encoding
    # -------------------------------------------------------------------

    def _encode_splice_constraints(
        self,
        optimizer: Optimize,
        nuc_exprs: list[ArithRef],
        protein: str,
        config: SolverConfig,
    ) -> list[str]:
        """Encode splice site avoidance constraints using MaxEntScan.

        For each potential GT (donor) position in the nucleotide sequence,
        encodes the MaxEntScan score as a Z3 expression and constrains
        it to be below the threshold. Similarly for AG (acceptor) sites.

        Donor encoding (9-mer, positions -3 to +6 relative to GT)::

            Implies(And(nuc_j == G, nuc_{j+1} == T),
                     donor_score_j < threshold)

        Acceptor encoding (23-mer, positions -20 to +3 relative to AG)::

            Implies(And(nuc_j == A, nuc_{j+1} == G),
                     acceptor_score_j < threshold)

        The score expressions use If-Then-Else chains over the
        pre-computed PWM log-odds contributions, scaled to integers
        for Z3 arithmetic.

        Args:
            optimizer: Z3 Optimize instance.
            nuc_exprs: Z3 expressions for each nucleotide position.
            protein: Protein sequence.
            config: Solver configuration with splice thresholds.

        Returns:
            List of constraint names added.
        """
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)
        donor_thresh = config.effective_donor_threshold
        acceptor_thresh = config.effective_acceptor_threshold

        # --- Donor site constraints (9-mer context around GT) ---
        for j in range(total_nucs - 1):
            ctx_start = j - 3
            ctx_end = j + 6  # exclusive
            if ctx_start < 0 or ctx_end > total_nucs:
                continue

            score_terms: list[ArithRef] = []
            for pwm_pos in range(9):
                nuc = nuc_exprs[ctx_start + pwm_pos]
                ss = [int(round(_DONOR_SCORES[pwm_pos][b] * _SCORE_SCALE))
                      for b in range(4)]
                score_terms.append(
                    If(nuc == 0, IntVal(ss[0]),
                       If(nuc == 1, IntVal(ss[1]),
                          If(nuc == 2, IntVal(ss[2]), IntVal(ss[3]))))
                )

            donor_score = Sum(score_terms)
            thresh_scaled = int(round(donor_thresh * _SCORE_SCALE))
            gt_detected = And(nuc_exprs[j] == 2, nuc_exprs[j + 1] == 3)  # G=2, T=3
            optimizer.add(Implies(gt_detected, donor_score < thresh_scaled))
            constraint_names.append(f"no_donor_splice_{j}")

        # --- Acceptor site constraints (23-mer context around AG) ---
        for j in range(total_nucs - 1):
            ctx_start = j - 20
            ctx_end = j + 4  # exclusive
            if ctx_start < 0 or ctx_end > total_nucs:
                continue

            score_terms: list[ArithRef] = []
            for pwm_pos in range(23):
                nuc = nuc_exprs[ctx_start + pwm_pos]
                ss = [int(round(_ACCEPTOR_SCORES[pwm_pos][b] * _SCORE_SCALE))
                      for b in range(4)]
                score_terms.append(
                    If(nuc == 0, IntVal(ss[0]),
                       If(nuc == 1, IntVal(ss[1]),
                          If(nuc == 2, IntVal(ss[2]), IntVal(ss[3]))))
                )

            acceptor_score = Sum(score_terms)
            thresh_scaled = int(round(acceptor_thresh * _SCORE_SCALE))
            ag_detected = And(nuc_exprs[j] == 0, nuc_exprs[j + 1] == 2)  # A=0, G=2
            optimizer.add(Implies(ag_detected, acceptor_score < thresh_scaled))
            constraint_names.append(f"no_acceptor_splice_{j}")

        if self.config.verbose:
            n_donors = sum(
                1 for j in range(total_nucs - 1)
                if j - 3 >= 0 and j + 6 <= total_nucs
            )
            n_acceptors = sum(
                1 for j in range(total_nucs - 1)
                if j - 20 >= 0 and j + 4 <= total_nucs
            )
            logger.info(
                "Z3: Added splice constraints (%d donor, %d acceptor positions)",
                n_donors, n_acceptors,
            )

        return constraint_names

    # -------------------------------------------------------------------
    # CAI objective encoding
    # -------------------------------------------------------------------

    def _encode_cai_objective(
        self,
        optimizer: Optimize,
        codon_vars: list[ArithRef],
        codon_lists: list[list[str]],
    ) -> None:
        """Encode the CAI maximization objective.

        CAI = (Π w_i)^(1/n), where w_i is the relative adaptiveness
        of the chosen codon. Maximizing CAI is equivalent to maximizing
        Σ log(w_i), which we scale to integers for Z3.

        For each codon position i, the CAI contribution is::

            If(codon_i == 0, log(cai_0) * SCALE,
               If(codon_i == 1, log(cai_1) * SCALE, ...))

        Args:
            optimizer: Z3 Optimize instance.
            codon_vars: Z3 Int variables for codon choices.
            codon_lists: Synonymous codon lists per position.
        """
        cai_data = CODON_ADAPTIVENESS_TABLES.get(self.organism, {})
        cai_terms: list[ArithRef] = []

        for c_var, codons_for_aa in zip(codon_vars, codon_lists):
            n_codons = len(codons_for_aa)
            log_cai_values: list[int] = []
            for codon in codons_for_aa:
                w = max(cai_data.get(codon, 0.01), 0.001)
                log_cai_values.append(int(round(math.log(w) * _CAI_SCALE)))

            if n_codons == 1:
                cai_terms.append(IntVal(log_cai_values[0]))
            else:
                expr: ArithRef = IntVal(log_cai_values[-1])
                for k in range(n_codons - 2, -1, -1):
                    expr = If(c_var == k, IntVal(log_cai_values[k]), expr)
                cai_terms.append(expr)

        optimizer.maximize(Sum(cai_terms))

        if self.config.verbose:
            logger.info(
                "Z3: Set CAI maximization objective for %d codon positions",
                len(codon_vars),
            )

    # -------------------------------------------------------------------
    # Solution extraction
    # -------------------------------------------------------------------

    def _extract_solution(
        self,
        z3_model: object,
        codon_vars: list[ArithRef],
        codon_lists: list[list[str]],
        protein: str,
        solve_time: float,
        num_constraints: int,
        num_variables: int,
    ) -> SolverResult:
        """Extract the codon assignment from a satisfiable Z3 model.

        Evaluates each codon variable in the Z3 model, maps the integer
        value back to a codon string, and computes the final sequence
        and metrics.

        Args:
            z3_model: The Z3 model from a successful solve.
            codon_vars: Z3 Int variables for codon choices.
            codon_lists: Synonymous codon lists per position.
            protein: Protein sequence.
            solve_time: Wall-clock solve time.
            num_constraints: Total constraints posted.
            num_variables: Total decision variables.

        Returns:
            SolverResult with the optimized sequence and metrics.
        """
        chosen_codons: list[str] = []
        for c_var, codons_for_aa in zip(codon_vars, codon_lists):
            val = z3_model.eval(c_var, model_completion=True)
            try:
                idx = val.as_long()
            except (AttributeError, ValueError) as exc:
                logger.warning(
                    "Z3: Could not extract integer value from model "
                    "evaluation for %s (got %r): %s — defaulting to 0",
                    c_var, val, exc,
                )
                idx = 0
            idx = max(0, min(idx, len(codons_for_aa) - 1))
            chosen_codons.append(codons_for_aa[idx])

        sequence = "".join(chosen_codons)

        # Compute GC content
        gc_count = sum(1 for base in sequence if base in ("G", "C"))
        gc_content = gc_count / len(sequence) if sequence else 0.0

        # Compute CAI
        cai_data = CODON_ADAPTIVENESS_TABLES.get(self.organism, {})
        log_cai_sum = sum(
            math.log(max(cai_data.get(codon, 0.01), 0.001))
            for codon in chosen_codons
        )
        n = len(chosen_codons)
        cai = math.exp(log_cai_sum / n) if n > 0 else 0.0

        return SolverResult(
            sequence=sequence,
            solved=True,
            backend_used=SolverBackend.Z3,
            cai=cai,
            gc_content=gc_content,
            solve_time_seconds=solve_time,
            num_constraints=num_constraints,
            num_variables=num_variables,
        )

    # -------------------------------------------------------------------
    # UNSAT core extraction
    # -------------------------------------------------------------------

    def _extract_unsat_core(
        self,
        optimizer: Optimize,
        constraint_names: list[str],
        n_domain_assertions: int = 0,
    ) -> MUSReport:
        """Extract UNSAT core from an unsatisfiable Z3 model.

        Z3's optimizer does not directly support UNSAT core extraction
        (only the Solver class does). We re-create the constraints using
        a Z3 Solver with tracked assertions to obtain the core.

        The core is then mapped back to human-readable constraint names
        for the MUS report.

        Args:
            optimizer: The Z3 Optimize instance (used to get assertions).
            constraint_names: List of constraint name strings that were
                added, in the order they were added.
            n_domain_assertions: Number of domain assertions added before
                the named constraints (used as an index offset).

        Returns:
            MUSReport with the minimal unsatisfiable subset.
        """
        violations: list[ConstraintViolation] = []
        conflicting_names: list[str] = []
        suggested_relaxations: list[str] = []

        try:
            from z3 import Solver, Bool as Z3Bool

            solver = Solver()
            assertions = optimizer.assertions()
            tracked: list[object] = []

            for i, assertion in enumerate(assertions):
                # Skip domain assertions — they are always satisfiable
                # on their own and don't correspond to named constraints
                if i < n_domain_assertions:
                    p = Z3Bool(f"domain_{i}")
                    solver.add(Implies(p, assertion))
                    tracked.append(p)
                    continue
                name_idx = i - n_domain_assertions
                name = constraint_names[name_idx] if name_idx < len(constraint_names) else f"constraint_{i}"
                p = Z3Bool(f"p_{name_idx}")
                solver.add(Implies(p, assertion))
                tracked.append(p)

            result = solver.check(*tracked)

            if result == unsat:
                core = solver.unsat_core()
                core_names: set[str] = set()
                for p in core:
                    try:
                        p_str = str(p)
                        if p_str.startswith("p_"):
                            name_idx = int(p_str[2:])
                            if name_idx < len(constraint_names):
                                name = constraint_names[name_idx]
                            else:
                                name = f"constraint_{name_idx}"
                            core_names.add(name)
                    except (ValueError, IndexError) as exc:
                        logger.debug(
                            "Z3: Could not parse UNSAT core element %r: %s",
                            p, exc,
                        )

                for name in sorted(core_names):
                    violation = self._classify_constraint(name)
                    violations.append(violation)
                    conflicting_names.append(name)

                # Suggest relaxations: splice constraints first, then
                # restriction sites, then GC (GC is usually hardest to relax)
                priority_order = {"splice_donor": 1, "splice_acceptor": 2,
                                  "restriction_site": 3, "gc": 4, "unknown": 5}
                sorted_violations = sorted(
                    violations,
                    key=lambda v: priority_order.get(
                        self._constraint_category(v.constraint_name), 5
                    ),
                    reverse=True,
                )
                suggested_relaxations = [v.constraint_name for v in sorted_violations]

        except Exception as e:
            logger.warning("Z3: UNSAT core extraction failed: %s", e)

        # Build explanation
        if violations:
            categories = set(
                self._constraint_category(v.constraint_name) for v in violations
            )
            explanation = (
                f"Unsatisfiable model — conflicting constraints in: "
                f"{', '.join(sorted(categories))}. "
                f"Try relaxing GC bounds, removing restriction enzymes, "
                f"or increasing splice score thresholds."
            )
        else:
            explanation = (
                "Unsatisfiable model — could not extract UNSAT core. "
                "Try relaxing constraints (wider GC range, fewer forbidden "
                "patterns, higher splice thresholds)."
            )

        return MUSReport(
            conflicting_constraints=conflicting_names,
            explanation=explanation,
            suggested_relaxations=suggested_relaxations,
        )

    # -------------------------------------------------------------------
    # Constraint classification helpers
    # -------------------------------------------------------------------

    @staticmethod
    def _constraint_category(name: str) -> str:
        """Extract the constraint category from a constraint name.

        Args:
            name: Constraint name string.

        Returns:
            Category string: 'gc', 'restriction_site', 'splice_donor',
            'splice_acceptor', or 'unknown'.
        """
        if name.startswith("gc_lo") or name.startswith("gc_hi"):
            return "gc"
        elif name.startswith("no_donor_splice"):
            return "splice_donor"
        elif name.startswith("no_acceptor_splice"):
            return "splice_acceptor"
        elif name.startswith("no_"):
            return "restriction_site"
        return "unknown"

    @staticmethod
    def _classify_constraint(name: str) -> ConstraintViolation:
        """Classify a constraint name into a ConstraintViolation.

        Parses the naming convention used by the encoding methods to
        determine the constraint type and description.

        Args:
            name: Constraint name string (e.g. 'gc_lo_0.4',
                'no_EcoRI', 'no_donor_splice_7').

        Returns:
            ConstraintViolation with parsed type and description.
        """
        category = Z3Engine._constraint_category(name)

        if category == "gc":
            bound = "lower" if name.startswith("gc_lo") else "upper"
            return ConstraintViolation(
                constraint_name=name,
                constraint_type=ConstraintStrictness.HARD,
                description=f"GC content {bound} bound constraint",
            )
        elif category == "splice_donor":
            pos_str = name.split("_")[-1]
            return ConstraintViolation(
                constraint_name=name,
                constraint_type=ConstraintStrictness.HARD,
                description=f"Donor splice site constraint at nucleotide position {pos_str}",
                positions=[int(pos_str)] if pos_str.isdigit() else [],
            )
        elif category == "splice_acceptor":
            pos_str = name.split("_")[-1]
            return ConstraintViolation(
                constraint_name=name,
                constraint_type=ConstraintStrictness.HARD,
                description=f"Acceptor splice site constraint at nucleotide position {pos_str}",
                positions=[int(pos_str)] if pos_str.isdigit() else [],
            )
        elif category == "restriction_site":
            pattern = name[3:]  # Remove 'no_' prefix
            return ConstraintViolation(
                constraint_name=name,
                constraint_type=ConstraintStrictness.HARD,
                description=f"Restriction site avoidance for pattern '{pattern}'",
            )
        else:
            return ConstraintViolation(
                constraint_name=name,
                constraint_type=ConstraintStrictness.HARD,
                description=f"Constraint: {name}",
            )
