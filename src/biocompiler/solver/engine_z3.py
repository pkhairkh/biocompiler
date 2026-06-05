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

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Z3 import with graceful fallback
# ---------------------------------------------------------------------------
try:
    from z3 import (
        Optimize,
        Solver,
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
        Bool as Z3Bool,
        sat,
        unsat,
        unknown,
    )
    _Z3_AVAILABLE = True
except ImportError as _z3_import_err:
    _Z3_AVAILABLE = False
    _Z3_IMPORT_ERROR = str(_z3_import_err)
    logger.debug("Z3 import failed: %s", _z3_import_err)
else:
    _Z3_IMPORT_ERROR = ""

from .types import (
    SolverConfig,
    SolverResult,
    SolverBackend,
    SolverBackendProtocol,
    CSPModel,
    ConstraintStrictness,
    InfeasibilityReport,
    MUSReport,
    ConstraintViolation,
)
from ..constants import AA_TO_CODONS, BASE_REV, IUPAC_EXPAND, RESTRICTION_ENZYMES, INSTABILITY_MOTIF, reverse_complement
from ..organisms import CODON_ADAPTIVENESS_TABLES

__all__ = ["Z3Engine", "compute_adaptive_timeout", "ADAPTIVE_TIMEOUT_TIERS", "ConstraintTier"]

# ── Constraint Tiering ───────────────────────────────────────────────
# Tier 1 (Hard): Must satisfy — translation fidelity, no stop codons, valid bases
# Tier 2 (Soft): Should satisfy — GC range, restriction site avoidance
# Tier 3 (Preference): Nice to have — CpG avoidance, ATTTA/T-run avoidance,
#   splice site avoidance, GT dinucleotide avoidance
#
# When the solver returns UNSAT with all constraints, we progressively relax
# Tier 3 then Tier 2 constraints until a solution is found (if one exists
# for Tier 1 alone).


class ConstraintTier:
    """Constraint priority tier for Z3 solver relaxation.

    When the Z3 solver returns UNSAT due to conflicting constraints, the
    engine progressively relaxes lower-priority tiers until a feasible
    solution is found.  Tier 1 constraints are NEVER relaxed (they encode
    biological impossibilities like wrong amino acid translation).

    Attributes:
        TIER1_HARD: Must satisfy — translation fidelity, no stop codons,
            valid bases.  Violation means the output is biologically
            meaningless.
        TIER2_SOFT: Should satisfy — GC content range, restriction site
            avoidance.  Violation is experimentally problematic but the
            sequence still encodes the correct protein.
        TIER3_PREFERENCE: Nice to have — CpG dinucleotide avoidance,
            ATTTA/T-run avoidance, splice site avoidance, GT dinucleotide
            avoidance.  Violation has minor biological impact and is
            acceptable when the problem would otherwise be infeasible.
    """
    TIER1_HARD = 1       # Never relaxed
    TIER2_SOFT = 2       # Relaxed after Tier 3
    TIER3_PREFERENCE = 3  # Relaxed first


# Mapping from constraint name patterns to their tier.
# Names are matched by prefix (e.g. "gc_lo" matches "gc_lo_0_3").
_CONSTRAINT_TIER_MAP: dict[str, int] = {
    # Tier 1: Hard — translation fidelity (implicit in codon domains)
    # Tier 1 constraints are never added as Z3 constraints because they're
    # enforced by the codon variable domains. We list them here for
    # documentation and for the relaxation report.
    "amino_acid": ConstraintTier.TIER1_HARD,
    "translation": ConstraintTier.TIER1_HARD,
    "stop_codon": ConstraintTier.TIER1_HARD,
    "valid_base": ConstraintTier.TIER1_HARD,

    # Tier 2: Soft — GC range and restriction site avoidance
    "gc_lo": ConstraintTier.TIER2_SOFT,
    "gc_hi": ConstraintTier.TIER2_SOFT,
    "no_": ConstraintTier.TIER2_SOFT,  # restriction site: no_EcoRI, no_BamHI, etc.

    # Tier 3: Preferences
    "no_within_cpg": ConstraintTier.TIER3_PREFERENCE,
    "no_cross_cpg": ConstraintTier.TIER3_PREFERENCE,
    "no_ATTTA": ConstraintTier.TIER3_PREFERENCE,
    "no_TAAAT": ConstraintTier.TIER3_PREFERENCE,
    "no_trun": ConstraintTier.TIER3_PREFERENCE,
    "no_donor_splice": ConstraintTier.TIER3_PREFERENCE,
    "no_acceptor_splice": ConstraintTier.TIER3_PREFERENCE,
    "instability_motif": ConstraintTier.TIER3_PREFERENCE,
    "mrna_stability": ConstraintTier.TIER3_PREFERENCE,
    "cpg": ConstraintTier.TIER3_PREFERENCE,
    "splice": ConstraintTier.TIER3_PREFERENCE,
}


# Soft constraint weights for Z3 Optimize.add_soft()
_TIER2_SOFT_WEIGHT: int = 200   # High weight — strongly preferred
_TIER3_SOFT_WEIGHT: int = 50    # Lower weight — nice to have


# Scaling factor for converting float CAI weights to Z3 integers
_CAI_SCALE: int = 10000  # 4 decimal places of precision
# Scale factor for MaxEntScan scores to convert to Z3-compatible integers
_SCORE_SCALE: int = 100  # 2 decimal places
# Default Z3 timeout in seconds (overridden by config.timeout_seconds if set)
_Z3_DEFAULT_TIMEOUT: float = 30.0
# Timeout for large protein instances (> _LARGE_PROTEIN_THRESHOLD aa)
_Z3_LARGE_INSTANCE_TIMEOUT: float = 60.0
# Protein length threshold above which we increase the Z3 timeout
_LARGE_PROTEIN_THRESHOLD: int = 500

# ── Adaptive timeout tiers ────────────────────────────────────────────
# Protein length → timeout in seconds.  The tiers are chosen so that
# small proteins get a quick solve, medium proteins get more time, and
# very large proteins (>500 aa) skip Z3 entirely because it is too slow.
ADAPTIVE_TIMEOUT_TIERS: list[tuple[int, float]] = [
    (100, 30.0),   # proteins ≤ 100 aa → 30 seconds
    (300, 60.0),   # proteins ≤ 300 aa → 60 seconds
    (500, 120.0),  # proteins ≤ 500 aa → 120 seconds
    # proteins > 500 aa → skip Z3 (use greedy directly)
]
"""Adaptive timeout tiers: list of (max_protein_length, timeout_seconds).

Used by :func:`compute_adaptive_timeout` to select an appropriate
Z3 timeout based on protein size.  Proteins longer than the largest
threshold (500 aa) should use the greedy optimizer directly instead
of Z3.
"""


def compute_adaptive_timeout(
    protein_length: int,
    solver_timeout: Optional[float] = None,
) -> tuple[float, bool]:
    """Compute an adaptive Z3 timeout based on protein length.

    Selects the appropriate timeout tier based on the protein's amino-acid
    count.  Proteins longer than the largest tier threshold should skip Z3
    entirely and use the greedy optimizer instead.

    Args:
        protein_length: Number of amino acids in the protein.
        solver_timeout: Optional explicit timeout override from the caller.
            When provided and > 0, this takes precedence over the tiered
            default (but the skip-Z3 decision still applies for proteins
            > 500 aa unless the caller explicitly overrides).

    Returns:
        A tuple of (timeout_seconds, should_skip_z3).
        - timeout_seconds: The Z3 timeout to use.
        - should_skip_z3: True if the protein is too large for Z3 and
          the greedy optimizer should be used directly.
    """
    # Check if the protein exceeds the largest tier — skip Z3 entirely
    max_tier_length = ADAPTIVE_TIMEOUT_TIERS[-1][0]
    should_skip = protein_length > max_tier_length

    # If the caller provided an explicit timeout, honour it (but still
    # recommend skipping for very large proteins).
    if solver_timeout is not None and solver_timeout > 0:
        return solver_timeout, should_skip

    # Walk the tiers in order; first tier whose length >= protein_length wins
    for max_len, timeout in ADAPTIVE_TIMEOUT_TIERS:
        if protein_length <= max_len:
            return timeout, should_skip

    # Fallback (shouldn't happen, but defensive)
    return ADAPTIVE_TIMEOUT_TIERS[-1][1], True


def _verify_z3_solution(
    sequence: str,
    protein: str,
    config: SolverConfig,
    organism: str,
) -> tuple[bool, list[str]]:
    """Verify that a Z3 solution satisfies all hard constraints.

    This is a defense-in-depth check against known Z3 bugs where the
    solver may report SAT but the solution doesn't actually satisfy all
    constraints (particularly with Optimize + soft constraints).

    Args:
        sequence: The DNA sequence from Z3.
        protein: The original amino acid sequence.
        config: The solver configuration (defines active constraints).
        organism: The target organism name.

    Returns:
        A tuple of (is_valid, violation_descriptions).
    """
    violations: list[str] = []

    if not sequence:
        return False, ["Empty sequence"]

    # 1. Translation fidelity
    if len(sequence) != len(protein) * 3:
        violations.append(
            f"Sequence length {len(sequence)} != expected {len(protein) * 3}"
        )
    else:
        from ..constants import CODON_TABLE
        for i, expected_aa in enumerate(protein):
            codon = sequence[i * 3 : i * 3 + 3]
            actual_aa = CODON_TABLE.get(codon)
            if actual_aa != expected_aa:
                violations.append(
                    f"Codon {i} ({codon}) translates to '{actual_aa}', "
                    f"expected '{expected_aa}'"
                )

    # 2. GC content
    gc_frac = sum(1 for b in sequence.upper() if b in "GC") / len(sequence)
    if gc_frac < config.gc_lo or gc_frac > config.gc_hi:
        violations.append(
            f"GC content {gc_frac:.3f} outside [{config.gc_lo:.2f}, {config.gc_hi:.2f}]"
        )

    # 3. Restriction sites (basic check)
    from ..constants import RESTRICTION_ENZYMES, reverse_complement
    restriction_sites = config.restriction_sites
    if not restriction_sites and getattr(config, 'add_default_restriction_sites', True):
        restriction_sites = [
            seq for name, seq in RESTRICTION_ENZYMES.items()
            if all(ch in "ACGT" for ch in seq) and len(seq) >= 6
        ]
    seq_upper = sequence.upper()
    for site in restriction_sites:
        site_upper = site.upper()
        if site_upper in seq_upper:
            violations.append(f"Restriction site '{site_upper}' found in sequence")
            break  # One hit is enough to flag
        rc = reverse_complement(site).upper()
        if rc != site_upper and rc in seq_upper:
            violations.append(
                f"Restriction site RC '{rc}' found in sequence"
            )
            break

    # 4. ATTTA motifs
    if config.avoid_attta and "ATTTA" in seq_upper:
        violations.append("ATTTA instability motif found in sequence")

    # 5. T-runs
    if config.avoid_t_runs:
        max_run = 0
        run = 0
        for b in seq_upper:
            if b == "T":
                run += 1
                max_run = max(max_run, run)
            else:
                run = 0
        if max_run >= 6:
            violations.append(
                f"T-run of length {max_run} found (threshold: 6)"
            )

    is_valid = len(violations) == 0
    return is_valid, violations


def _sanitize_z3_name(name: str) -> str:
    """Sanitize a constraint name for use as a Z3 Boolean variable name.

    Z3 Boolean names must be valid SMT-LIB identifiers: alphanumeric
    characters and underscores only.  Dots, hyphens, spaces, etc. are
    replaced with underscores.

    Args:
        name: Raw constraint name (e.g. ``"gc_lo_0.4"``).

    Returns:
        Sanitized name (e.g. ``"gc_lo_0_4"``).
    """
    import re
    return re.sub(r'[^a-zA-Z0-9_]', '_', name)


def _get_constraint_tier(name: str) -> int:
    """Determine the priority tier for a constraint based on its name.

    Uses longest-prefix matching against ``_CONSTRAINT_TIER_MAP``.
    If no match is found, defaults to Tier 2 (soft — safe to relax
    if needed but not as expendable as preferences).

    Args:
        name: Sanitized constraint name (e.g. ``"gc_lo_0_3"``).

    Returns:
        Tier number (1, 2, or 3).
    """
    best_tier = ConstraintTier.TIER2_SOFT  # default
    best_len = 0
    for prefix, tier in _CONSTRAINT_TIER_MAP.items():
        if name.startswith(prefix) and len(prefix) > best_len:
            best_tier = tier
            best_len = len(prefix)
    return best_tier


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

        Note:
            Unlike previous versions, the constructor no longer raises
            ``ImportError`` when z3-solver is unavailable.  Instead,
            :meth:`solve` returns an unsolved ``SolverResult`` with
            ``fallback_used=True``.  This avoids silent failures in
            the dispatch layer where ``ImportError`` was caught and
            the entire backend silently skipped.
        """
        # Derive organism from config if the default was used and config
        # carries a different organism.
        if organism == "Homo_sapiens" and config.organism != "Homo_sapiens":
            organism = config.organism
        self.config = config
        self.organism = organism
        self._seed: int = seed
        self._optimizer: Optional[Optimize] = None
        # Storage for constraint expressions and names for UNSAT core extraction
        self._constraint_exprs: list[tuple[str, object]] = []
        # Whether any soft constraints were detected in the model
        self._has_soft_constraints: bool = False
        # Set of ConstraintType values marked as soft (populated during solve)
        self._soft_constraint_types: set = set()
        # Current relaxation level (0=all hard, 1=Tier3 soft, 2=Tier2+3 soft)
        self._relaxation_level: int = 0
        # Constraints relaxed during solving (populated by _solve_with_relaxation)
        self._relaxed_constraints: list[str] = []
        logger.debug(
            "Z3Engine initialized: organism=%s, seed=%d, "
            "z3_available=%s, timeout=%.1fs",
            organism, seed, _Z3_AVAILABLE, config.timeout_seconds,
        )

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

        # ── Availability guard ──────────────────────────────────────────
        # Return an unsolved result instead of raising ImportError so the
        # dispatch layer can fall back gracefully.
        if not _Z3_AVAILABLE:
            logger.warning(
                "Z3Engine.solve() called but z3-solver is not available: %s",
                _Z3_IMPORT_ERROR,
            )
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                fallback_used=True,
                solve_time_seconds=time.perf_counter() - start_time,
                warnings=[
                    "z3-solver package is not installed. "
                    "Install with: pip install z3-solver"
                ],
            )

        # ── Extract model attributes with getattr fallbacks ─────────────
        # Different CSPModel types (types.CSPModel vs constraints.CSPModel)
        # expose different attribute names. Use getattr with fallbacks so
        # the engine works with either model type.
        config = getattr(model, "config", None) or self.config
        protein = (
            getattr(model, "protein_sequence", None)
            or getattr(model, "protein", "")
        ).upper()
        organism = getattr(model, "organism", None) or self.organism

        # Sync organism from model to config so that is_eukaryotic detection
        # works correctly even when only the model specifies the organism.
        if organism and organism != config.organism:
            try:
                config.organism = organism
            except (AttributeError, TypeError):
                logger.warning(
                    "SolverConfig.organism='%s' differs from model organism='%s' "
                    "but config is frozen; eukaryotic detection may be incorrect",
                    config.organism, organism,
                )

        # Log Z3 engine diagnostics
        codon_domains_attr = getattr(model, "codon_domains", {})
        n_codon_domains = len(codon_domains_attr) if codon_domains_attr else 0
        model_constraints = getattr(model, "constraints", [])
        logger.info(
            "Z3 engine diagnostics: protein_len=%d organism=%s "
            "codon_domains=%d model_constraints=%d "
            "config.gc_range=[%.2f,%.2f] config.avoid_cpg=%s "
            "config.avoid_t_runs=%s config.avoid_attta=%s "
            "config.cryptic_splice_threshold=%.1f timeout=%.1fs",
            len(protein) if protein else 0, organism,
            n_codon_domains, len(model_constraints),
            config.gc_lo, config.gc_hi,
            config.avoid_cpg, config.avoid_t_runs,
            config.avoid_attta, config.cryptic_splice_threshold,
            config.timeout_seconds,
        )

        logger.info(
            "Z3: Starting solve for %d aa protein, organism=%s, "
            "gc=[%.2f,%.2f], timeout=%.1fs",
            len(protein), organism,
            config.gc_lo, config.gc_hi, config.timeout_seconds,
        )

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

        # Store protein for translation fidelity validation in _extract_solution
        self._protein_for_validation = protein

        # Reset constraint storage for this solve invocation
        self._constraint_exprs = []

        # ── Fix #2: Detect soft constraints in the model ─────────────
        # If any constraint in the model is marked SOFT, we should use
        # add_soft() for those constraints so the optimizer can relax
        # them if needed.  We also collect the soft constraint types so
        # that the encoding methods can use add_soft() instead of add()
        # for constraints that correspond to soft model constraints.
        model_constraints = getattr(model, "constraints", [])
        self._has_soft_constraints = any(
            getattr(c, "strictness", ConstraintStrictness.HARD) == ConstraintStrictness.SOFT
            for c in model_constraints
        )
        # Collect which ConstraintType values are soft so encoding methods
        # can decide add() vs add_soft().
        self._soft_constraint_types: set = set()
        for c in model_constraints:
            if getattr(c, "strictness", ConstraintStrictness.HARD) == ConstraintStrictness.SOFT:
                ctype = getattr(c, "ctype", None)
                if ctype is not None:
                    self._soft_constraint_types.add(ctype)

        # ── Adaptive timeout for Z3 ────────────────────────────────
        # Use tiered timeouts based on protein length:
        #   ≤100aa → 30s, ≤300aa → 60s, ≤500aa → 120s, >500aa → skip Z3
        # The solver_timeout parameter on config overrides the tier default.
        solver_timeout_override = getattr(config, 'solver_timeout', None)
        effective_timeout, should_skip_z3 = compute_adaptive_timeout(
            n, solver_timeout=solver_timeout_override or config.timeout_seconds,
        )
        # If config.timeout_seconds was set to an explicit positive value,
        # respect it (unless the protein is too large).
        if config.timeout_seconds > 0:
            effective_timeout, _ = compute_adaptive_timeout(
                n, solver_timeout=config.timeout_seconds,
            )

        if should_skip_z3:
            logger.info(
                "Z3: Protein too large (%d aa > %d threshold), "
                "skipping Z3 and recommending greedy fallback. "
                "effective_timeout=%.1fs",
                n, ADAPTIVE_TIMEOUT_TIERS[-1][0], effective_timeout,
            )
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                fallback_used=True,
                solve_time_seconds=time.perf_counter() - start_time,
                num_variables=n,
                warnings=[
                    f"Protein too large for Z3 ({n} aa > "
                    f"{ADAPTIVE_TIMEOUT_TIERS[-1][0]} threshold). "
                    "Use greedy optimizer instead."
                ],
                metadata={
                    "reason": "protein_too_large_for_z3",
                    "protein_length": n,
                    "adaptive_timeout": effective_timeout,
                    "skip_z3": True,
                },
            )

        logger.info(
            "Z3: Adaptive timeout selected: %.1fs for %d aa protein "
            "(solver_timeout=%s, config.timeout=%.1fs)",
            effective_timeout, n,
            solver_timeout_override, config.timeout_seconds,
        )

        # ── Fix #2: Always use Optimize (supports both hard and soft) ──
        # Optimize handles add() for hard constraints and add_soft() for
        # soft constraints with weighted relaxation.  When no soft
        # constraints exist, Optimize behaves like Solver with an objective.
        optimizer = Optimize()
        optimizer.set("timeout", int(effective_timeout * 1000))
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

        # ── Fix #4: Store codon vars/lists for splice lookup table ──
        self._current_codon_vars = codon_vars
        self._current_codon_lists = codon_lists

        # --- Add constraints ---
        # Strategy: Add Tier 1 (translation) as hard, Tier 2+3 as soft.
        # This ensures the solver always finds a solution if one exists
        # for the translation constraints alone, while preferring to
        # satisfy GC, restriction site, CpG, ATTTA, T-run, and splice
        # constraints when possible.
        constraint_names: list[str] = []

        # Track number of domain assertions for UNSAT core mapping
        # Each codon position adds 2 domain assertions (>=0, <N)
        n_domain_assertions = 2 * n

        # Collect ALL soft constraint names for the relaxed encoding
        all_soft_names: set[str] = set()

        # ── Constraint Tiering ─────────────────────────────────────────
        # Tier 1 (Hard): Translation fidelity, no stop codons, valid bases
        #   → Already enforced by codon variable domains (optimizer.add
        #     in the loop above: c_var >= 0, c_var < len(codons_for_aa))
        # Tier 2 (Soft): GC range, restriction site avoidance
        #   → optimizer.add_soft() with _TIER2_SOFT_WEIGHT
        # Tier 3 (Preference): CpG avoidance, ATTTA/T-run avoidance,
        #   splice site avoidance
        #   → optimizer.add_soft() with _TIER3_SOFT_WEIGHT
        #
        # Key fix: We ONLY use add_soft() for Tier 2/3 constraints.
        # Previously, both add_soft() AND add() were called for the same
        # constraint, making the hard add() take precedence and preventing
        # relaxation — the root cause of UNSAT for 20aa proteins.

        # GC constraint (Tier 2 → soft only)
        gc_names = self._encode_gc_constraint_soft(optimizer, nuc_exprs, protein, config)
        all_soft_names.update(gc_names)
        constraint_names.extend(gc_names)

        # Restriction site constraints (Tier 2 → soft only)
        restriction_sites = config.restriction_sites
        if not restriction_sites and getattr(config, 'add_default_restriction_sites', True):
            restriction_sites = [
                seq for name, seq in RESTRICTION_ENZYMES.items()
                if all(ch in "ACGT" for ch in seq) and len(seq) >= 6
            ]
        if restriction_sites:
            rs_names = self._encode_restriction_site_constraints_soft(
                optimizer, nuc_exprs, restriction_sites
            )
            all_soft_names.update(rs_names)
            constraint_names.extend(rs_names)

        # ATTTA instability motif avoidance (Tier 3 → soft only)
        if config.avoid_attta:
            attta_names = self._encode_attta_constraints_soft(optimizer, nuc_exprs)
            all_soft_names.update(attta_names)
            constraint_names.extend(attta_names)

        # T-run avoidance (Tier 3 → soft only)
        if config.avoid_t_runs:
            trun_names = self._encode_trun_constraints_soft(optimizer, nuc_exprs)
            all_soft_names.update(trun_names)
            constraint_names.extend(trun_names)

        # CpG dinucleotide avoidance (Tier 3 → soft only, eukaryotes only)
        if config.avoid_cpg and config.is_eukaryotic:
            cpg_names = self._encode_cpg_constraints_soft(optimizer, codon_vars, codon_lists)
            all_soft_names.update(cpg_names)
            constraint_names.extend(cpg_names)

        # Splice site constraints (Tier 3 → soft only, eukaryotes only)
        if config.is_eukaryotic:
            splice_names = self._encode_splice_constraints_soft(
                optimizer, nuc_exprs, protein, config
            )
            all_soft_names.update(splice_names)
            constraint_names.extend(splice_names)

        # ── Classify constraints by tier for reporting ──────────────────
        tier3_names: set[str] = set()
        tier2_names: set[str] = set()
        tier1_names: set[str] = set()
        for cname in constraint_names:
            tier = _get_constraint_tier(cname)
            if tier == ConstraintTier.TIER3_PREFERENCE:
                tier3_names.add(cname)
            elif tier == ConstraintTier.TIER2_SOFT:
                tier2_names.add(cname)
            else:
                tier1_names.add(cname)

        # --- Set optimization objective: MAXIMIZE sum of log(CAI) ---
        self._encode_cai_objective(optimizer, codon_vars, codon_lists)

        # --- Solve ---
        logger.info(
            "Z3: Solving codon optimization for %d aa protein "
            "(%d nucleotide variables, %d constraints, timeout=%.1fs, "
            "tier1=%d, tier2_soft=%d, tier3_soft=%d)",
            n, total_nucs, len(constraint_names), effective_timeout,
            len(tier1_names), len(tier2_names), len(tier3_names),
        )

        result_status = optimizer.check()
        solve_time = time.perf_counter() - start_time

        # Log solver status for diagnostics
        status_str = str(result_status)
        logger.info(
            "Z3 solve completed: status=%s time=%.2fs "
            "variables=%d constraints=%d adaptive_timeout=%.1fs",
            status_str, solve_time, n, len(constraint_names),
            effective_timeout,
        )

        if result_status == sat:
            z3_model = optimizer.model()
            result = self._extract_solution(
                z3_model, codon_vars, codon_lists, protein,
                solve_time, len(constraint_names), n
            )

            # ── Solution verification (Tier 1 only) ──────────────────
            # With soft constraints, only verify translation fidelity.
            # Tier 2 and 3 violations are acceptable.
            is_valid, verification_violations = _verify_z3_solution(
                result.sequence, protein, config, organism,
            )
            if is_valid:
                logger.info(
                    "Z3: Solution verification PASSED for %d aa protein",
                    n,
                )
                result.metadata["verification_passed"] = True
            else:
                # Check if violations are only in soft constraints
                soft_violations = [
                    v for v in verification_violations
                    if any(kw in v.lower() for kw in [
                        "gc content", "restriction site",
                        "attta", "t-run", "cpg", "splice"
                    ])
                ]
                if soft_violations and len(soft_violations) == len(verification_violations):
                    # Only soft constraint violations — solution is valid
                    logger.info(
                        "Z3: Solution has soft constraint violations "
                        "(expected with soft constraints) for %d aa protein",
                        n,
                    )
                    result.metadata["verification_passed"] = True
                    result.metadata["soft_violations"] = soft_violations
                    result.warnings.extend(soft_violations)
                else:
                    logger.warning(
                        "Z3: Solution verification FAILED for %d aa protein: %s",
                        n, verification_violations,
                    )
                    result.solved = False
                    result.fallback_used = True
                    result.metadata["verification_passed"] = False
                    result.metadata["verification_violations"] = verification_violations
                    result.metadata["z3_partial_sequence"] = result.sequence
                    result.warnings.append(
                        "Z3 solution failed verification: "
                        + "; ".join(verification_violations)
                    )
            return result
        elif result_status == unsat:
            # UNSAT with Tier 2+3 as soft.  Only Tier 1 (translation
            # domains) remains hard.  This means even the codon domain
            # constraints conflict — which should be impossible for a
            # valid protein.
            mus = self._extract_unsat_core(optimizer, constraint_names, n_domain_assertions)
            logger.info(
                "Z3: UNSAT result for %d aa protein. MUS: %s",
                n,
                [c.name for c in mus.mus_constraints] if mus.mus_constraints else mus.conflicting_constraints,
            )
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                solve_time_seconds=solve_time,
                num_constraints=len(constraint_names),
                num_variables=n,
                mus_report=mus,
                fallback_used=True,
                warnings=["Model is unsatisfiable — conflicting constraints"],
                metadata={
                    "reason": "unsat",
                    "relaxation_level": 0,
                },
            )

        elif result_status == unknown:
            reason = optimizer.reason_unknown()
            is_timeout = "timeout" in reason.lower() or "canceled" in reason.lower()

            # ── On timeout, try to extract best solution found ────────
            # Z3's Optimize may have a partial model even after timeout.
            # If we can extract it, return it as a best-effort solution
            # with the partial sequence stored in metadata for warm-start.
            if is_timeout:
                logger.info(
                    "Z3: Timeout after %.1fs for %d aa protein — "
                    "attempting partial solution extraction for warm-start",
                    solve_time, n,
                )
                try:
                    partial_model = optimizer.model()
                    if partial_model and len(partial_model) > 0:
                        result = self._extract_solution(
                            partial_model, codon_vars, codon_lists, protein,
                            solve_time, len(constraint_names), n
                        )
                        # Mark as unsolved but include the best-effort sequence
                        result.solved = False
                        result.fallback_used = True
                        result.warnings.append(
                            f"Z3 timed out after {solve_time:.1f}s — "
                            "returning best solution found so far "
                            "(may not satisfy all constraints)"
                        )
                        # Store partial sequence for warm-start in greedy
                        result.metadata["z3_partial_sequence"] = result.sequence
                        result.metadata["reason"] = "z3_timeout"
                        result.metadata["constraint_counts_by_tier"] = dict(tier1=len(tier1_names), tier2=len(tier2_names), tier3=len(tier3_names))
                        logger.info(
                            "Z3: Extracted partial solution after timeout "
                            "(sequence_len=%d) — stored for warm-start",
                            len(result.sequence),
                        )
                        return result
                except Exception as e:
                    logger.debug(
                        "Z3: Could not extract partial model after timeout: %s", e
                    )

            logger.info(
                "Z3: Unknown result ('%s') for %d aa protein — "
                "falling back to greedy. timeout=%s",
                reason, n, is_timeout,
            )
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                solve_time_seconds=solve_time,
                num_constraints=len(constraint_names),
                num_variables=n,
                warnings=[f"Solver returned unknown: {reason}"],
                fallback_used=is_timeout,
                metadata={
                    "reason": "z3_timeout" if is_timeout else "z3_unknown",
                    "constraint_counts_by_tier": dict(tier1=len(tier1_names), tier2=len(tier2_names), tier3=len(tier3_names)),
                },
            )
        else:
            logger.warning(
                "Z3: Unexpected result status '%s' for %d aa protein",
                result_status, n,
            )
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                solve_time_seconds=solve_time,
                warnings=[f"Unexpected Z3 result: {result_status}"],
                metadata={"reason": "unexpected_result"},
            )

    # -------------------------------------------------------------------
    # Progressive constraint relaxation
    # -------------------------------------------------------------------
    # Soft constraint encoding methods (Tier 2+3 as add_soft)
    # -------------------------------------------------------------------

    def _classify_gc_constraints(self, nuc_exprs, protein, config):
        """Return constraint names for GC constraints."""
        return [_sanitize_z3_name(f"gc_lo_{config.gc_lo}"),
                _sanitize_z3_name(f"gc_hi_{config.gc_hi}")]

    def _classify_restriction_site_constraints(self, nuc_exprs, forbidden_patterns):
        """Return constraint names for restriction site constraints."""
        return [_sanitize_z3_name(f"no_{p}") for p in forbidden_patterns]

    def _classify_cpg_constraints(self, codon_vars, codon_lists):
        """Return constraint names for CpG constraints."""
        names = []
        for i, (c_var, codons_for_aa) in enumerate(zip(codon_vars, codon_lists)):
            cg_free = [c for c in codons_for_aa if "CG" not in c]
            if cg_free and len(cg_free) < len(codons_for_aa):
                names.append(_sanitize_z3_name(f"no_within_cpg_{i}"))
        for i in range(len(codon_vars) - 1):
            domain_i = codon_lists[i]
            domain_j = codon_lists[i + 1]
            c_ending = [k for k, c in enumerate(domain_i) if c[-1] == "C"]
            g_starting = [k for k, c in enumerate(domain_j) if c[0] == "G"]
            if c_ending and g_starting:
                non_c_ending = [k for k, c in enumerate(domain_i) if c[-1] != "C"]
                non_g_starting = [k for k, c in enumerate(domain_j) if c[0] != "G"]
                if non_c_ending or non_g_starting:
                    names.append(_sanitize_z3_name(f"no_cross_cpg_{i}"))
        return names

    def _classify_splice_constraints(self, nuc_exprs, protein, config):
        """Return constraint names for splice constraints."""
        names = []
        total_nucs = len(nuc_exprs)
        donor_thresh = config.effective_donor_threshold
        acceptor_thresh = config.effective_acceptor_threshold
        if donor_thresh <= 0 and acceptor_thresh <= 0:
            return names
        if config.auto_detect_organism_domain and not config.is_eukaryotic:
            return names
        if donor_thresh > 0:
            for j in range(total_nucs - 1):
                names.append(_sanitize_z3_name(f"no_donor_splice_{j}"))
        if acceptor_thresh > 0:
            for j in range(total_nucs - 1):
                names.append(_sanitize_z3_name(f"no_acceptor_splice_{j}"))
        return names

    def _encode_gc_constraint_soft(self, optimizer, nuc_exprs, protein, config):
        """Encode GC constraints as soft (Tier 2).

        Returns:
            List of constraint names added (for UNSAT core tracking and
            tier classification).
        """
        total_nucs = 3 * len(protein)
        gc_indicators = [If(Or(nuc == 1, nuc == 2), IntVal(1), IntVal(0)) for nuc in nuc_exprs]
        gc_count = Sum(gc_indicators)
        gc_lo_bound = math.ceil(config.gc_lo * total_nucs)
        gc_hi_bound = math.floor(config.gc_hi * total_nucs)
        gc_lo_expr = gc_count >= gc_lo_bound
        gc_hi_expr = gc_count <= gc_hi_bound
        lo_name = _sanitize_z3_name(f"gc_lo_{config.gc_lo}")
        hi_name = _sanitize_z3_name(f"gc_hi_{config.gc_hi}")
        optimizer.add_soft(gc_lo_expr, _TIER2_SOFT_WEIGHT)
        optimizer.add_soft(gc_hi_expr, _TIER2_SOFT_WEIGHT)
        self._constraint_exprs.append((lo_name, gc_lo_expr))
        self._constraint_exprs.append((hi_name, gc_hi_expr))
        return [lo_name, hi_name]

    def _encode_restriction_site_constraints_soft(self, optimizer, nuc_exprs, forbidden_patterns):
        """Encode restriction site constraints as soft (Tier 2).

        Returns:
            List of constraint names added (for UNSAT core tracking and
            tier classification).
        """
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)
        for pattern in forbidden_patterns:
            pattern_upper = pattern.upper()
            pat_len = len(pattern_upper)
            if pat_len == 0 or pat_len > total_nucs:
                continue
            expanded_bases = []
            for base_char in pattern_upper:
                iupac_bases = IUPAC_EXPAND.get(base_char, base_char)
                base_ints = [BASE_REV[b] for b in iupac_bases if b in BASE_REV]
                if not base_ints:
                    break
                expanded_bases.append(base_ints)
            if len(expanded_bases) != pat_len:
                continue
            pat_name = _sanitize_z3_name(f"no_{pattern}")
            first_expr = None
            for j in range(total_nucs - pat_len + 1):
                base_eqs = []
                for k in range(pat_len):
                    allowed = expanded_bases[k]
                    if len(allowed) == 1:
                        base_eqs.append(nuc_exprs[j + k] == allowed[0])
                    else:
                        base_eqs.append(Or(*[nuc_exprs[j + k] == b for b in allowed]))
                window_expr = Not(And(*base_eqs))
                optimizer.add_soft(window_expr, _TIER2_SOFT_WEIGHT)
                if first_expr is None:
                    first_expr = window_expr
            if first_expr is not None:
                self._constraint_exprs.append((pat_name, first_expr))
            constraint_names.append(pat_name)
        return constraint_names

    def _encode_attta_constraints_soft(self, optimizer, nuc_exprs):
        """Encode ATTTA motif avoidance as soft (Tier 3).

        Returns:
            List of constraint names added (for UNSAT core tracking and
            tier classification).
        """
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)
        attta_name = _sanitize_z3_name("no_ATTTA")
        taaat_name = _sanitize_z3_name("no_TAAAT")
        first_attta = None
        for j in range(total_nucs - 4):
            expr = Not(And(nuc_exprs[j]==0, nuc_exprs[j+1]==3, nuc_exprs[j+2]==3, nuc_exprs[j+3]==3, nuc_exprs[j+4]==0))
            optimizer.add_soft(expr, _TIER3_SOFT_WEIGHT)
            if first_attta is None:
                first_attta = expr
        if first_attta is not None:
            self._constraint_exprs.append((attta_name, first_attta))
        constraint_names.append(attta_name)
        first_taaat = None
        for j in range(total_nucs - 4):
            expr = Not(And(nuc_exprs[j]==3, nuc_exprs[j+1]==0, nuc_exprs[j+2]==0, nuc_exprs[j+3]==0, nuc_exprs[j+4]==3))
            optimizer.add_soft(expr, _TIER3_SOFT_WEIGHT)
            if first_taaat is None:
                first_taaat = expr
        if first_taaat is not None:
            self._constraint_exprs.append((taaat_name, first_taaat))
        constraint_names.append(taaat_name)
        return constraint_names

    def _encode_trun_constraints_soft(self, optimizer, nuc_exprs, max_run=5):
        """Encode T-run avoidance as soft (Tier 3).

        Returns:
            List of constraint names added (for UNSAT core tracking and
            tier classification).
        """
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)
        run_len = max_run + 1
        trun_name = _sanitize_z3_name(f"no_trun_{run_len}")
        first_expr = None
        for j in range(total_nucs - run_len + 1):
            all_t = And(*[nuc_exprs[j + k] == 3 for k in range(run_len)])
            trun_expr = Not(all_t)
            optimizer.add_soft(trun_expr, _TIER3_SOFT_WEIGHT)
            if first_expr is None:
                first_expr = trun_expr
        if first_expr is not None:
            self._constraint_exprs.append((trun_name, first_expr))
        constraint_names.append(trun_name)
        return constraint_names

    def _encode_cpg_constraints_soft(self, optimizer, codon_vars, codon_lists):
        """Encode CpG avoidance as soft (Tier 3).

        Returns:
            List of constraint names added (for UNSAT core tracking and
            tier classification).
        """
        constraint_names: list[str] = []
        for i, (c_var, codons_for_aa) in enumerate(zip(codon_vars, codon_lists)):
            cg_free = [c for c in codons_for_aa if "CG" not in c]
            if cg_free and len(cg_free) < len(codons_for_aa):
                cg_free_indices = [codons_for_aa.index(c) for c in cg_free]
                cpg_expr = Or(*[c_var == idx for idx in cg_free_indices])
                cpg_name = _sanitize_z3_name(f"no_within_cpg_{i}")
                optimizer.add_soft(cpg_expr, _TIER3_SOFT_WEIGHT)
                self._constraint_exprs.append((cpg_name, cpg_expr))
                constraint_names.append(cpg_name)
        for i in range(len(codon_vars) - 1):
            domain_i = codon_lists[i]
            domain_j = codon_lists[i + 1]
            c_ending = [k for k, c in enumerate(domain_i) if c[-1] == "C"]
            g_starting = [k for k, c in enumerate(domain_j) if c[0] == "G"]
            if c_ending and g_starting:
                non_c_ending = [k for k, c in enumerate(domain_i) if c[-1] != "C"]
                non_g_starting = [k for k, c in enumerate(domain_j) if c[0] != "G"]
                if non_c_ending or non_g_starting:
                    cross_name = _sanitize_z3_name(f"no_cross_cpg_{i}")
                    rep_expr = None
                    for ci in c_ending:
                        for gi in g_starting:
                            cross_expr = Not(And(codon_vars[i] == ci, codon_vars[i + 1] == gi))
                            optimizer.add_soft(cross_expr, _TIER3_SOFT_WEIGHT)
                            if rep_expr is None:
                                rep_expr = cross_expr
                    if rep_expr is not None:
                        self._constraint_exprs.append((cross_name, rep_expr))
                    constraint_names.append(cross_name)
        return constraint_names

    def _encode_splice_constraints_soft(self, optimizer, nuc_exprs, protein, config):
        """Encode splice site constraints as soft (Tier 3).

        Returns:
            List of constraint names added (for UNSAT core tracking and
            tier classification).
        """
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)
        donor_thresh = config.effective_donor_threshold
        acceptor_thresh = config.effective_acceptor_threshold
        if donor_thresh <= 0 and acceptor_thresh <= 0:
            return constraint_names
        if config.auto_detect_organism_domain and not config.is_eukaryotic:
            return constraint_names
        codon_vars = getattr(self, "_current_codon_vars", [])
        codon_lists = getattr(self, "_current_codon_lists", [])
        if not codon_vars or not codon_lists:
            return constraint_names
        n_codons = len(codon_vars)
        # Donor sites
        if donor_thresh > 0:
            for j in range(total_nucs - 1):
                ctx_start = j - 3
                ctx_end = j + 6
                if ctx_start < 0 or ctx_end > total_nucs:
                    score_terms = []
                    for pwm_pos in range(9):
                        nuc_idx = j - 3 + pwm_pos
                        if nuc_idx < 0 or nuc_idx >= total_nucs:
                            continue
                        nuc = nuc_exprs[nuc_idx]
                        ss = [int(round(_DONOR_SCORES[pwm_pos][b] * _SCORE_SCALE)) for b in range(4)]
                        score_terms.append(
                            If(nuc==0, IntVal(ss[0]), If(nuc==1, IntVal(ss[1]), If(nuc==2, IntVal(ss[2]), IntVal(ss[3])))))
                    if score_terms:
                        donor_score = Sum(score_terms)
                        thresh_scaled = int(round(donor_thresh * _SCORE_SCALE))
                        gt_detected = And(nuc_exprs[j] == 2, nuc_exprs[j + 1] == 3)
                        edge_expr = Implies(gt_detected, donor_score < thresh_scaled)
                        optimizer.add_soft(edge_expr, _TIER3_SOFT_WEIGHT)
                        self._constraint_exprs.append((_sanitize_z3_name(f"no_donor_splice_{j}"), edge_expr))
                    continue
                first_codon = ctx_start // 3
                last_codon = (ctx_end - 1) // 3
                covered_codons = list(range(first_codon, min(last_codon + 1, n_codons)))
                forbidden_combos = self._find_forbidden_splice_combos(codon_lists, covered_codons, j, "donor", donor_thresh)
                rep_expr = None
                for combo in forbidden_combos:
                    terms = [codon_vars[ci] == idx for ci, idx in combo.items()]
                    terms.append(nuc_exprs[j] == 2)
                    terms.append(nuc_exprs[j + 1] == 3)
                    constraint_expr = Not(And(*terms))
                    optimizer.add_soft(constraint_expr, _TIER3_SOFT_WEIGHT)
                    if rep_expr is None:
                        rep_expr = constraint_expr
                donor_name = _sanitize_z3_name(f"no_donor_splice_{j}")
                if rep_expr is not None:
                    self._constraint_exprs.append((donor_name, rep_expr))
                constraint_names.append(donor_name)
        # Acceptor sites
        if acceptor_thresh > 0:
            for j in range(total_nucs - 1):
                ctx_start = j - 20
                ctx_end = j + 4
                if ctx_start < 0 or ctx_end > total_nucs:
                    score_terms = []
                    for pwm_pos in range(23):
                        nuc_idx = j - 20 + pwm_pos
                        if nuc_idx < 0 or nuc_idx >= total_nucs:
                            continue
                        nuc = nuc_exprs[nuc_idx]
                        ss = [int(round(_ACCEPTOR_SCORES[pwm_pos][b] * _SCORE_SCALE)) for b in range(4)]
                        score_terms.append(
                            If(nuc==0, IntVal(ss[0]), If(nuc==1, IntVal(ss[1]), If(nuc==2, IntVal(ss[2]), IntVal(ss[3])))))
                    if score_terms:
                        acceptor_score = Sum(score_terms)
                        thresh_scaled = int(round(acceptor_thresh * _SCORE_SCALE))
                        ag_detected = And(nuc_exprs[j] == 0, nuc_exprs[j + 1] == 2)
                        edge_expr = Implies(ag_detected, acceptor_score < thresh_scaled)
                        optimizer.add_soft(edge_expr, _TIER3_SOFT_WEIGHT)
                        self._constraint_exprs.append((_sanitize_z3_name(f"no_acceptor_splice_{j}"), edge_expr))
                    continue
                ag_codon = j // 3
                inner_codons = set([ag_codon])
                if j + 1 < total_nucs and (j + 1) // 3 != ag_codon:
                    inner_codons.add((j + 1) // 3)
                for c in list(inner_codons):
                    if c > 0: inner_codons.add(c - 1)
                    if c < n_codons - 1: inner_codons.add(c + 1)
                inner_codons = sorted(inner_codons)
                if len(inner_codons) > 6:
                    score_terms = []
                    for pwm_pos in range(23):
                        nuc = nuc_exprs[ctx_start + pwm_pos]
                        ss = [int(round(_ACCEPTOR_SCORES[pwm_pos][b] * _SCORE_SCALE)) for b in range(4)]
                        score_terms.append(
                            If(nuc==0, IntVal(ss[0]), If(nuc==1, IntVal(ss[1]), If(nuc==2, IntVal(ss[2]), IntVal(ss[3])))))
                    acceptor_score = Sum(score_terms)
                    thresh_scaled = int(round(acceptor_thresh * _SCORE_SCALE))
                    ag_detected = And(nuc_exprs[j] == 0, nuc_exprs[j + 1] == 2)
                    acc_expr = Implies(ag_detected, acceptor_score < thresh_scaled)
                    optimizer.add_soft(acc_expr, _TIER3_SOFT_WEIGHT)
                    self._constraint_exprs.append((_sanitize_z3_name(f"no_acceptor_splice_{j}"), acc_expr))
                    continue
                forbidden_combos = self._find_forbidden_splice_combos(codon_lists, inner_codons, j, "acceptor", acceptor_thresh)
                rep_expr = None
                for combo in forbidden_combos:
                    terms = [codon_vars[ci] == idx for ci, idx in combo.items()]
                    terms.append(nuc_exprs[j] == 0)
                    terms.append(nuc_exprs[j + 1] == 2)
                    constraint_expr = Not(And(*terms))
                    optimizer.add_soft(constraint_expr, _TIER3_SOFT_WEIGHT)
                    if rep_expr is None:
                        rep_expr = constraint_expr
                acceptor_name = _sanitize_z3_name(f"no_acceptor_splice_{j}")
                if rep_expr is not None:
                    self._constraint_exprs.append((acceptor_name, rep_expr))
                constraint_names.append(acceptor_name)
        return constraint_names

    # -------------------------------------------------------------------
    # Relaxed solving (progressive constraint relaxation)
    # -------------------------------------------------------------------

    def _solve_with_relaxation(
        self,
        model: CSPModel,
        protein: str,
        organism: str,
        config: SolverConfig,
        codon_vars: list,
        codon_lists: list,
        nuc_exprs: list,
        n: int,
        total_nucs: int,
        effective_timeout: float,
        start_time: float,
        tier2_names: set[str],
        tier3_names: set[str],
        constraint_names: list[str],
    ) -> Optional[SolverResult]:
        """Try solving with progressively relaxed constraints.

        Returns a SolverResult if a solution is found, or None if all
        relaxation levels still produce UNSAT.
        """
        self._relaxed_constraints = []

        for relaxation_level in (1, 2):
            self._relaxation_level = relaxation_level
            level_desc = (
                "Tier 3 soft, Tier 2 hard" if relaxation_level == 1
                else "Tier 2+3 soft, only Tier 1 hard"
            )
            logger.info(
                "Z3: Trying relaxation level %d (%s) for %d aa protein",
                relaxation_level, level_desc, n,
            )

            # Rebuild optimizer with same codon variables but soft
            # constraints for relaxed tiers.
            relaxed_opt = Optimize()
            relaxed_opt.set("timeout", int(effective_timeout * 1000))
            try:
                relaxed_opt.set("random_seed", self._seed)
            except Exception:
                pass

            # Re-add codon domain constraints (Tier 1 — always hard)
            for i, (c_var, codons_for_aa) in enumerate(zip(codon_vars, codon_lists)):
                relaxed_opt.add(c_var >= 0)
                relaxed_opt.add(c_var < len(codons_for_aa))

            # Reset constraint expression storage
            self._constraint_exprs = []
            relaxed_constraint_names: list[str] = []

            # Determine which names should be soft
            soft_names = set(tier3_names) if relaxation_level == 1 else set(tier3_names) | set(tier2_names)
            self._relaxed_constraints = sorted(soft_names)

            # GC constraint
            relaxed_constraint_names.extend(
                self._encode_gc_constraint_relaxed(
                    relaxed_opt, nuc_exprs, protein, config,
                    soft_names,
                )
            )

            # Restriction site constraints
            restriction_sites = config.restriction_sites
            if not restriction_sites and getattr(config, 'add_default_restriction_sites', True):
                restriction_sites = [
                    seq for name, seq in RESTRICTION_ENZYMES.items()
                    if all(ch in "ACGT" for ch in seq) and len(seq) >= 6
                ]
            if restriction_sites:
                relaxed_constraint_names.extend(
                    self._encode_restriction_site_constraints_relaxed(
                        relaxed_opt, nuc_exprs, restriction_sites,
                        soft_names,
                    )
                )

            # ATTTA instability motif avoidance
            if config.avoid_attta:
                relaxed_constraint_names.extend(
                    self._encode_attta_constraints_relaxed(
                        relaxed_opt, nuc_exprs, soft_names,
                    )
                )

            # T-run avoidance
            if config.avoid_t_runs:
                relaxed_constraint_names.extend(
                    self._encode_trun_constraints_relaxed(
                        relaxed_opt, nuc_exprs, soft_names,
                    )
                )

            # CpG avoidance (eukaryotes only)
            if config.avoid_cpg and config.is_eukaryotic:
                relaxed_constraint_names.extend(
                    self._encode_cpg_constraints_relaxed(
                        relaxed_opt, codon_vars, codon_lists, soft_names,
                    )
                )

            # Splice site constraints (eukaryotes only)
            if config.is_eukaryotic:
                relaxed_constraint_names.extend(
                    self._encode_splice_constraints(
                        relaxed_opt, nuc_exprs, protein, config
                    )
                )

            # CAI objective
            self._encode_cai_objective(relaxed_opt, codon_vars, codon_lists)

            # Solve
            logger.info(
                "Z3: Solving with relaxation level %d (%d constraints, "
                "%d soft) for %d aa protein",
                relaxation_level, len(relaxed_constraint_names),
                len(soft_names), n,
            )

            result_status = relaxed_opt.check()
            solve_time = time.perf_counter() - start_time

            if result_status == sat:
                z3_model = relaxed_opt.model()
                result = self._extract_solution(
                    z3_model, codon_vars, codon_lists, protein,
                    solve_time, len(relaxed_constraint_names), n
                )

                # Verify the solution
                is_valid, verification_violations = _verify_z3_solution(
                    result.sequence, protein, config, organism,
                )
                if is_valid:
                    result.metadata["relaxation_level"] = relaxation_level
                    result.metadata["relaxed_constraints"] = self._relaxed_constraints
                    result.metadata["verification_passed"] = True
                    if relaxation_level >= 1:
                        result.warnings.append(
                            f"Z3 solver used constraint relaxation "
                            f"(level {relaxation_level}: {level_desc}). "
                            f"Relaxed constraints: {self._relaxed_constraints}"
                        )
                    logger.info(
                        "Z3: SAT at relaxation level %d for %d aa protein. "
                        "Relaxed: %s",
                        relaxation_level, n, self._relaxed_constraints,
                    )
                    return result
                else:
                    logger.warning(
                        "Z3: SAT at relaxation level %d but verification "
                        "failed: %s", relaxation_level, verification_violations,
                    )
            elif result_status == unsat:
                logger.info(
                    "Z3: Still UNSAT at relaxation level %d for %d aa protein",
                    relaxation_level, n,
                )
            else:
                reason = relaxed_opt.reason_unknown()
                logger.info(
                    "Z3: Unknown result at relaxation level %d: %s",
                    relaxation_level, reason,
                )

        # All relaxation levels exhausted
        return None

    # -------------------------------------------------------------------
    # Relaxed constraint encoding helpers
    # -------------------------------------------------------------------

    def _add_constraint_by_tier(
        self,
        optimizer: Optimize,
        expr,
        name: str,
        soft_names: set[str],
    ) -> None:
        """Add a constraint as hard or soft based on its tier.

        If *name* is in *soft_names*, the constraint is added via
        ``add_soft()`` with a weight determined by its tier.  Otherwise,
        it is added as a hard constraint via ``add()``.
        """
        if name in soft_names:
            tier = _get_constraint_tier(name)
            weight = (
                _TIER3_SOFT_WEIGHT
                if tier == ConstraintTier.TIER3_PREFERENCE
                else _TIER2_SOFT_WEIGHT
            )
            optimizer.add_soft(expr, weight)
        else:
            optimizer.add(expr)

    def _encode_gc_constraint_relaxed(
        self,
        optimizer: Optimize,
        nuc_exprs: list,
        protein: str,
        config: SolverConfig,
        soft_names: set[str],
    ) -> list[str]:
        """Encode GC content constraints with tier-based relaxation."""
        total_nucs = 3 * len(protein)
        constraint_names: list[str] = []

        gc_indicators: list = []
        for nuc in nuc_exprs:
            gc_indicators.append(If(Or(nuc == 1, nuc == 2), IntVal(1), IntVal(0)))

        gc_count = Sum(gc_indicators)
        gc_lo_bound = math.ceil(config.gc_lo * total_nucs)
        gc_hi_bound = math.floor(config.gc_hi * total_nucs)

        gc_lo_expr = gc_count >= gc_lo_bound
        gc_hi_expr = gc_count <= gc_hi_bound

        lo_name = _sanitize_z3_name(f"gc_lo_{config.gc_lo}")
        hi_name = _sanitize_z3_name(f"gc_hi_{config.gc_hi}")

        self._add_constraint_by_tier(optimizer, gc_lo_expr, lo_name, soft_names)
        self._add_constraint_by_tier(optimizer, gc_hi_expr, hi_name, soft_names)

        self._constraint_exprs.append((lo_name, gc_lo_expr))
        self._constraint_exprs.append((hi_name, gc_hi_expr))
        constraint_names.extend([lo_name, hi_name])

        return constraint_names

    def _encode_restriction_site_constraints_relaxed(
        self,
        optimizer: Optimize,
        nuc_exprs: list,
        forbidden_patterns: list[str],
        soft_names: set[str],
    ) -> list[str]:
        """Encode restriction site constraints with tier-based relaxation."""
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)

        for pattern in forbidden_patterns:
            pattern_upper = pattern.upper()
            pat_len = len(pattern_upper)
            if pat_len == 0 or pat_len > total_nucs:
                continue

            expanded_bases: list[list[int]] = []
            for base_char in pattern_upper:
                iupac_bases = IUPAC_EXPAND.get(base_char, base_char)
                base_ints = [BASE_REV[b] for b in iupac_bases if b in BASE_REV]
                if not base_ints:
                    break
                expanded_bases.append(base_ints)

            if len(expanded_bases) != pat_len:
                continue

            pat_name = _sanitize_z3_name(f"no_{pattern}")
            first_expr = None
            for j in range(total_nucs - pat_len + 1):
                base_eqs: list = []
                for k in range(pat_len):
                    allowed = expanded_bases[k]
                    if len(allowed) == 1:
                        base_eqs.append(nuc_exprs[j + k] == allowed[0])
                    else:
                        base_eqs.append(
                            Or(*[nuc_exprs[j + k] == b for b in allowed])
                        )
                window_expr = Not(And(*base_eqs))
                self._add_constraint_by_tier(optimizer, window_expr, pat_name, soft_names)
                if first_expr is None:
                    first_expr = window_expr

            if first_expr is not None:
                self._constraint_exprs.append((pat_name, first_expr))
            constraint_names.append(pat_name)

        return constraint_names

    def _encode_attta_constraints_relaxed(
        self,
        optimizer: Optimize,
        nuc_exprs: list,
        soft_names: set[str],
    ) -> list[str]:
        """Encode ATTTA motif avoidance with tier-based relaxation."""
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)

        attta_name = _sanitize_z3_name("no_ATTTA")
        taaat_name = _sanitize_z3_name("no_TAAAT")

        first_attta_expr = None
        for j in range(total_nucs - 4):
            attta_expr = Not(And(
                nuc_exprs[j] == 0, nuc_exprs[j + 1] == 3,
                nuc_exprs[j + 2] == 3, nuc_exprs[j + 3] == 3,
                nuc_exprs[j + 4] == 0,
            ))
            self._add_constraint_by_tier(optimizer, attta_expr, attta_name, soft_names)
            if first_attta_expr is None:
                first_attta_expr = attta_expr

        if first_attta_expr is not None:
            self._constraint_exprs.append((attta_name, first_attta_expr))
        constraint_names.append(attta_name)

        first_taaat_expr = None
        for j in range(total_nucs - 4):
            taaat_expr = Not(And(
                nuc_exprs[j] == 3, nuc_exprs[j + 1] == 0,
                nuc_exprs[j + 2] == 0, nuc_exprs[j + 3] == 0,
                nuc_exprs[j + 4] == 3,
            ))
            self._add_constraint_by_tier(optimizer, taaat_expr, taaat_name, soft_names)
            if first_taaat_expr is None:
                first_taaat_expr = taaat_expr

        if first_taaat_expr is not None:
            self._constraint_exprs.append((taaat_name, first_taaat_expr))
        constraint_names.append(taaat_name)

        return constraint_names

    def _encode_trun_constraints_relaxed(
        self,
        optimizer: Optimize,
        nuc_exprs: list,
        soft_names: set[str],
        max_run: int = 5,
    ) -> list[str]:
        """Encode T-run avoidance with tier-based relaxation."""
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)
        run_len = max_run + 1

        trun_name = _sanitize_z3_name(f"no_trun_{run_len}")
        first_expr = None
        for j in range(total_nucs - run_len + 1):
            all_t = And(*[nuc_exprs[j + k] == 3 for k in range(run_len)])
            trun_expr = Not(all_t)
            self._add_constraint_by_tier(optimizer, trun_expr, trun_name, soft_names)
            if first_expr is None:
                first_expr = trun_expr

        if first_expr is not None:
            self._constraint_exprs.append((trun_name, first_expr))
        constraint_names.append(trun_name)

        return constraint_names

    def _encode_cpg_constraints_relaxed(
        self,
        optimizer: Optimize,
        codon_vars: list,
        codon_lists: list,
        soft_names: set[str],
    ) -> list[str]:
        """Encode CpG avoidance with tier-based relaxation."""
        constraint_names: list[str] = []

        # Within-codon
        for i, (c_var, codons_for_aa) in enumerate(zip(codon_vars, codon_lists)):
            cg_free = [c for c in codons_for_aa if "CG" not in c]
            if cg_free and len(cg_free) < len(codons_for_aa):
                cg_free_indices = [codons_for_aa.index(c) for c in cg_free]
                cpg_expr = Or(*[c_var == idx for idx in cg_free_indices])
                cpg_name = _sanitize_z3_name(f"no_within_cpg_{i}")
                self._add_constraint_by_tier(optimizer, cpg_expr, cpg_name, soft_names)
                self._constraint_exprs.append((cpg_name, cpg_expr))
                constraint_names.append(cpg_name)

        # Cross-codon
        for i in range(len(codon_vars) - 1):
            domain_i = codon_lists[i]
            domain_j = codon_lists[i + 1]
            c_ending = [k for k, c in enumerate(domain_i) if c[-1] == "C"]
            g_starting = [k for k, c in enumerate(domain_j) if c[0] == "G"]

            if c_ending and g_starting:
                non_c_ending = [k for k, c in enumerate(domain_i) if c[-1] != "C"]
                non_g_starting = [k for k, c in enumerate(domain_j) if c[0] != "G"]
                if non_c_ending or non_g_starting:
                    cross_name = _sanitize_z3_name(f"no_cross_cpg_{i}")
                    rep_expr = None
                    for ci in c_ending:
                        for gi in g_starting:
                            cross_expr = Not(And(codon_vars[i] == ci, codon_vars[i + 1] == gi))
                            self._add_constraint_by_tier(optimizer, cross_expr, cross_name, soft_names)
                            if rep_expr is None:
                                rep_expr = cross_expr
                    if rep_expr is not None:
                        self._constraint_exprs.append((cross_name, rep_expr))
                    constraint_names.append(cross_name)

        return constraint_names

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

        GC count = sum(If(nuc_i in {C, G}, 1, 0))
        Constraint: gc_lo * 3*n <= gc_count <= gc_hi * 3*n

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

        gc_lo_expr = gc_count >= gc_lo_bound
        gc_hi_expr = gc_count <= gc_hi_bound
        optimizer.add(gc_lo_expr)
        optimizer.add(gc_hi_expr)

        # Store constraint expressions for UNSAT core with sanitized names
        lo_name = _sanitize_z3_name(f"gc_lo_{config.gc_lo}")
        hi_name = _sanitize_z3_name(f"gc_hi_{config.gc_hi}")
        self._constraint_exprs.append((lo_name, gc_lo_expr))
        self._constraint_exprs.append((hi_name, gc_hi_expr))

        constraint_names.append(lo_name)
        constraint_names.append(hi_name)

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

            # Store one representative expression per pattern for UNSAT core
            # (using sanitized name so Z3Bool creation doesn't fail)
            if base_eqs:
                pat_name = _sanitize_z3_name(f"no_{pattern}")
                self._constraint_exprs.append(
                    (pat_name, Not(And(*base_eqs)))
                )
            constraint_names.append(_sanitize_z3_name(f"no_{pattern}"))

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
        """Encode splice site avoidance constraints using lookup table approach.

        Instead of encoding MaxEntScan scores as Z3 arithmetic expressions
        (which creates deep If-Then-Else chains that Z3 struggles with),
        we precompute which codon-index combinations produce splice scores
        above the threshold, then forbid only those combinations.

        For each potential GT (donor) position j, the 9-mer context spans
        codons [(j-3)//3 .. (j+5)//3].  We enumerate all codon choices
        for those positions, compute the 9-mer sequence, score it with
        MaxEntScan, and if the score exceeds the threshold AND the middle
        two bases are GT, we add a forbidden combination constraint::

            Not(And(codon_i1 == idx1, codon_i2 == idx2, ..., nuc_j == G, nuc_{j+1} == T))

        This converts the arithmetic encoding into simple Boolean constraints
        that Z3 handles efficiently.

        Acceptor sites (23-mer around AG) are encoded similarly.

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

        # Skip splice constraints entirely when threshold <= 0 (disabled)
        # or when the organism is prokaryotic (no spliceosome).
        if donor_thresh <= 0 and acceptor_thresh <= 0:
            return constraint_names

        # Prokaryotic organisms lack spliceosomes; cryptic splice
        # sites are biologically irrelevant.
        if config.auto_detect_organism_domain and not config.is_eukaryotic:
            logger.info(
                "Z3: Skipping splice constraints for prokaryotic "
                "organism '%s'",
                config.organism,
            )
            return constraint_names

        # Import MaxEntScan scoring functions for lookup table precomputation
        from ..maxentscan import score_donor, score_acceptor

        # --- Build codon position mapping ---
        # We need to know which codon variable covers each nucleotide position.
        # codon_vars and codon_lists are accessible via self or the optimizer's
        # internal state. We reconstruct from nuc_exprs.
        # Actually, we need the codon_vars and codon_lists — we'll receive them
        # from the solve() method via a helper.

        # This method is called from solve() where codon_vars and codon_lists
        # are available. We'll use a new signature that passes them.
        # For now, we infer codon boundaries from the nucleotide expressions.
        # Each codon covers 3 nucleotide positions. Given n codon positions,
        # codon i covers nuc positions [3*i, 3*i+3).

        # We need codon_vars and codon_lists to enumerate combinations.
        # They must be passed from solve(). We retrieve them from self.
        codon_vars = getattr(self, "_current_codon_vars", [])
        codon_lists = getattr(self, "_current_codon_lists", [])

        if not codon_vars or not codon_lists:
            # Fallback: no codon variables stored yet (shouldn't happen)
            logger.warning("Z3: Splice encoding has no codon variables; skipping")
            return constraint_names

        n_codons = len(codon_vars)

        # --- Donor site constraints (9-mer context around GT) ---
        if donor_thresh > 0:
            n_donor_forbidden = 0
            for j in range(total_nucs - 1):
                # 9-mer context: positions j-3 to j+5
                ctx_start = j - 3
                ctx_end = j + 6  # exclusive

                # ── Fix #4: Handle edge cases (near start/end) ─────
                # When the full 9-mer context doesn't fit, fall back to
                # arithmetic encoding using only available positions.
                # Missing positions contribute a neutral score of 0,
                # making the constraint slightly more permissive at edges
                # (biologically justified: splice sites near CDS boundaries
                # are less likely to be functional).
                if ctx_start < 0 or ctx_end > total_nucs:
                    # Fall back to arithmetic encoding for edge positions
                    score_terms = []
                    for pwm_pos in range(9):
                        nuc_idx = j - 3 + pwm_pos
                        if nuc_idx < 0 or nuc_idx >= total_nucs:
                            continue  # Skip missing positions (neutral contribution)
                        nuc = nuc_exprs[nuc_idx]
                        ss = [int(round(_DONOR_SCORES[pwm_pos][b] * _SCORE_SCALE))
                              for b in range(4)]
                        score_terms.append(
                            If(nuc == 0, IntVal(ss[0]),
                               If(nuc == 1, IntVal(ss[1]),
                                  If(nuc == 2, IntVal(ss[2]), IntVal(ss[3]))))
                        )
                    if score_terms:
                        donor_score = Sum(score_terms)
                        thresh_scaled = int(round(donor_thresh * _SCORE_SCALE))
                        gt_detected = And(nuc_exprs[j] == 2, nuc_exprs[j + 1] == 3)
                        edge_expr = Implies(gt_detected, donor_score < thresh_scaled)
                        optimizer.add(edge_expr)
                        edge_name = _sanitize_z3_name(f"no_donor_splice_{j}")
                        self._constraint_exprs.append((edge_name, edge_expr))
                        constraint_names.append(edge_name)
                    continue

                # Find codon positions covered by this context window
                first_codon = ctx_start // 3
                last_codon = (ctx_end - 1) // 3

                # Enumerate all codon combinations for the covered codons
                covered_codons = list(range(first_codon, min(last_codon + 1, n_codons)))

                # Build the lookup table: enumerate all combinations and
                # identify which ones produce high-scoring donor sites
                forbidden_combos = self._find_forbidden_splice_combos(
                    codon_lists, covered_codons, j, "donor", donor_thresh
                )

                # Add constraints: for each forbidden combination, forbid it
                # Store a representative expression for UNSAT core tracking
                rep_expr = None
                for combo in forbidden_combos:
                    # combo is a dict mapping codon_index -> codon_choice_index
                    terms = []
                    for ci, idx in combo.items():
                        terms.append(codon_vars[ci] == idx)
                    # Also require GT at positions j, j+1
                    terms.append(nuc_exprs[j] == 2)   # G=2
                    terms.append(nuc_exprs[j + 1] == 3)  # T=3
                    constraint_expr = Not(And(*terms))
                    optimizer.add(constraint_expr)
                    n_donor_forbidden += 1
                    if rep_expr is None:
                        rep_expr = constraint_expr

                donor_name = _sanitize_z3_name(f"no_donor_splice_{j}")
                if rep_expr is not None:
                    self._constraint_exprs.append((donor_name, rep_expr))
                constraint_names.append(donor_name)

            if self.config.verbose:
                logger.info(
                    "Z3: Added donor splice constraints via lookup table "
                    "(%d forbidden combos across %d GT positions)",
                    n_donor_forbidden,
                    sum(1 for j in range(total_nucs - 1)
                        if j - 3 >= 0 and j + 6 <= total_nucs),
                )

        # --- Acceptor site constraints (23-mer context around AG) ---
        if acceptor_thresh > 0:
            n_acceptor_forbidden = 0
            for j in range(total_nucs - 1):
                ctx_start = j - 20
                ctx_end = j + 4  # exclusive

                # ── Fix #4: Handle edge cases (near start/end) ─────
                # When the full 23-mer context doesn't fit, fall back to
                # arithmetic encoding using only available positions.
                if ctx_start < 0 or ctx_end > total_nucs:
                    score_terms = []
                    for pwm_pos in range(23):
                        nuc_idx = j - 20 + pwm_pos
                        if nuc_idx < 0 or nuc_idx >= total_nucs:
                            continue  # Skip missing positions
                        nuc = nuc_exprs[nuc_idx]
                        ss = [int(round(_ACCEPTOR_SCORES[pwm_pos][b] * _SCORE_SCALE))
                              for b in range(4)]
                        score_terms.append(
                            If(nuc == 0, IntVal(ss[0]),
                               If(nuc == 1, IntVal(ss[1]),
                                  If(nuc == 2, IntVal(ss[2]), IntVal(ss[3]))))
                        )
                    if score_terms:
                        acceptor_score = Sum(score_terms)
                        thresh_scaled = int(round(acceptor_thresh * _SCORE_SCALE))
                        ag_detected = And(nuc_exprs[j] == 0, nuc_exprs[j + 1] == 2)
                        edge_expr = Implies(ag_detected, acceptor_score < thresh_scaled)
                        optimizer.add(edge_expr)
                        edge_name = _sanitize_z3_name(f"no_acceptor_splice_{j}")
                        self._constraint_exprs.append((edge_name, edge_expr))
                        constraint_names.append(edge_name)
                    continue

                # Find codon positions covered by this context window
                first_codon = ctx_start // 3
                last_codon = (ctx_end - 1) // 3
                covered_codons = list(range(first_codon, min(last_codon + 1, n_codons)))

                # The 23-mer spans many codons (up to 8-9). Enumerating all
                # combinations of 8+ codons could be expensive (6^8 ≈ 1.7M).
                # For acceptor sites, we use a hybrid approach: enumerate the
                # "inner" codons (those containing or adjacent to the AG) fully,
                # and for "outer" codons, just add the GT/AG guard.
                #
                # Specifically, we only enumerate combinations for codons that
                # contain the AG dinucleotide and its immediate neighbors (±1 codon).
                # For the remaining outer codons, the constraint is simply that
                # if the inner combination produces a high score, it's forbidden.
                #
                # In practice, most AG positions within a codon boundary involve
                # only 2-3 codons for the AG itself, so enumeration is tractable.

                ag_codon = j // 3  # Codon containing the A of AG
                # Inner codons: ag_codon and possibly ag_codon+1 if AG spans boundary
                inner_codons = set()
                inner_codons.add(ag_codon)
                if j + 1 < total_nucs and (j + 1) // 3 != ag_codon:
                    inner_codons.add((j + 1) // 3)
                # Add neighbors for score context
                for c in list(inner_codons):
                    if c > 0:
                        inner_codons.add(c - 1)
                    if c < n_codons - 1:
                        inner_codons.add(c + 1)
                inner_codons = sorted(inner_codons)

                # For efficiency, limit enumeration to inner codons
                # (at most 5 codons → 6^5 = 7776 combinations)
                if len(inner_codons) > 6:
                    # Too many inner codons; fall back to arithmetic encoding
                    # for this position
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
                    ag_detected = And(nuc_exprs[j] == 0, nuc_exprs[j + 1] == 2)
                    optimizer.add(Implies(ag_detected, acceptor_score < thresh_scaled))
                    self._constraint_exprs.append(
                        (f"no_acceptor_splice_{j}",
                         Implies(ag_detected, acceptor_score < thresh_scaled))
                    )
                    constraint_names.append(f"no_acceptor_splice_{j}")
                    continue

                forbidden_combos = self._find_forbidden_splice_combos(
                    codon_lists, inner_codons, j, "acceptor", acceptor_thresh
                )

                # Store a representative expression for UNSAT core tracking
                rep_expr = None
                for combo in forbidden_combos:
                    terms = []
                    for ci, idx in combo.items():
                        terms.append(codon_vars[ci] == idx)
                    terms.append(nuc_exprs[j] == 0)   # A=0
                    terms.append(nuc_exprs[j + 1] == 2)  # G=2
                    constraint_expr = Not(And(*terms))
                    optimizer.add(constraint_expr)
                    n_acceptor_forbidden += 1
                    if rep_expr is None:
                        rep_expr = constraint_expr

                acc_name = _sanitize_z3_name(f"no_acceptor_splice_{j}")
                if rep_expr is not None:
                    self._constraint_exprs.append((acc_name, rep_expr))
                constraint_names.append(acc_name)

            if self.config.verbose:
                logger.info(
                    "Z3: Added acceptor splice constraints via lookup table "
                    "(%d forbidden combos across %d AG positions)",
                    n_acceptor_forbidden,
                    sum(1 for j in range(total_nucs - 1)
                        if j - 20 >= 0 and j + 4 <= total_nucs),
                )

        return constraint_names

    # -------------------------------------------------------------------
    # Splice lookup table helper
    # -------------------------------------------------------------------

    def _find_forbidden_splice_combos(
        self,
        codon_lists: list[list[str]],
        covered_codons: list[int],
        nuc_position: int,
        site_type: str,
        threshold: float,
    ) -> list[dict[int, int]]:
        """Find codon combinations that produce splice scores above threshold.

        Enumerates all combinations of codon choices for the given codon
        positions, reconstructs the relevant subsequence, computes the
        MaxEntScan score, and returns combinations that exceed the threshold.

        Args:
            codon_lists: Synonymous codon lists per position.
            covered_codons: Indices of codon positions to enumerate.
            nuc_position: Nucleotide position of the GT (donor) or AG (acceptor).
            site_type: "donor" or "acceptor".
            threshold: MaxEntScan score threshold.

        Returns:
            List of dicts mapping codon_index -> codon_choice_index for
            forbidden combinations.
        """
        from ..maxentscan import score_donor, score_acceptor

        forbidden: list[dict[int, int]] = []

        # Build the Cartesian product of codon choices for covered positions
        def _enumerate_combos(
            codon_indices: list[int],
            partial: dict[int, int],
            depth: int,
        ) -> None:
            if depth >= len(codon_indices):
                # Build the full sequence from the current combo
                # to check the splice score
                combo = dict(partial)
                # Reconstruct the nucleotide sequence for the covered region
                first_codon = min(codon_indices)
                last_codon = max(codon_indices)
                seq_start = first_codon * 3
                seq_end = (last_codon + 1) * 3

                # Build sequence from combo choices
                seq_chars: list[str] = []
                for ci in range(first_codon, last_codon + 1):
                    if ci in combo:
                        codon_idx = combo[ci]
                        seq_chars.append(codon_lists[ci][codon_idx])
                    else:
                        # Use first codon as placeholder (shouldn't happen
                        # if covered_codons is correct)
                        seq_chars.append(codon_lists[ci][0])

                partial_seq = "".join(seq_chars)

                # Check if the GT/AG dinucleotide is present at the target position
                local_pos = nuc_position - seq_start
                if site_type == "donor":
                    if local_pos < 0 or local_pos + 1 >= len(partial_seq):
                        return
                    if partial_seq[local_pos] != 'G' or partial_seq[local_pos + 1] != 'T':
                        return  # No GT, constraint is vacuously satisfied
                    # Score the donor site
                    score = score_donor(partial_seq, local_pos)
                    if score >= threshold:
                        forbidden.append(combo)
                else:  # acceptor
                    if local_pos < 0 or local_pos + 1 >= len(partial_seq):
                        return
                    if partial_seq[local_pos] != 'A' or partial_seq[local_pos + 1] != 'G':
                        return  # No AG
                    score = score_acceptor(partial_seq, local_pos)
                    if score >= threshold:
                        forbidden.append(combo)
                return

            ci = codon_indices[depth]
            for idx in range(len(codon_lists[ci])):
                partial[ci] = idx
                _enumerate_combos(codon_indices, partial, depth + 1)
                del partial[ci]

        _enumerate_combos(covered_codons, {}, 0)
        return forbidden

    # -------------------------------------------------------------------
    # ATTTA instability motif avoidance
    # -------------------------------------------------------------------

    def _encode_attta_constraints(
        self,
        optimizer: Optimize,
        nuc_exprs: list[ArithRef],
    ) -> list[str]:
        """Encode ATTTA instability motif avoidance constraints.

        Forbids the pattern ATTTA (A=0, T=3, T=3, T=3, A=0) and its
        reverse complement TAAAT (T=3, A=0, A=0, A=0, T=3) from
        appearing anywhere in the nucleotide sequence.

        Args:
            optimizer: Z3 Optimize instance.
            nuc_exprs: Z3 expressions for each nucleotide position.

        Returns:
            List of constraint names added.
        """
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)

        # ATTTA: A=0, T=3, T=3, T=3, A=0
        # Check if instability motif constraint should be soft
        from .types import ConstraintType
        is_soft = ConstraintType.NO_INSTABILITY_MOTIF in self._soft_constraint_types
        soft_weight = 100  # Default soft constraint weight for ATTTA

        # Build a representative expression for UNSAT core tracking
        first_attta_expr = None
        for j in range(total_nucs - 4):
            attta_expr = Not(And(
                nuc_exprs[j] == 0,       # A
                nuc_exprs[j + 1] == 3,   # T
                nuc_exprs[j + 2] == 3,   # T
                nuc_exprs[j + 3] == 3,   # T
                nuc_exprs[j + 4] == 0,   # A
            ))
            if is_soft:
                optimizer.add_soft(attta_expr, soft_weight)
            else:
                optimizer.add(attta_expr)
            if first_attta_expr is None:
                first_attta_expr = attta_expr
        attta_name = _sanitize_z3_name("no_ATTTA")
        if first_attta_expr is not None:
            self._constraint_exprs.append((attta_name, first_attta_expr))
        constraint_names.append(attta_name)

        # TAAAT (reverse complement of ATTTA): T=3, A=0, A=0, A=0, T=3
        first_taaat_expr = None
        for j in range(total_nucs - 4):
            taaat_expr = Not(And(
                nuc_exprs[j] == 3,       # T
                nuc_exprs[j + 1] == 0,   # A
                nuc_exprs[j + 2] == 0,   # A
                nuc_exprs[j + 3] == 0,   # A
                nuc_exprs[j + 4] == 3,   # T
            ))
            if is_soft:
                optimizer.add_soft(taaat_expr, soft_weight)
            else:
                optimizer.add(taaat_expr)
            if first_taaat_expr is None:
                first_taaat_expr = taaat_expr
        taaat_name = _sanitize_z3_name("no_TAAAT")
        if first_taaat_expr is not None:
            self._constraint_exprs.append((taaat_name, first_taaat_expr))
        constraint_names.append(taaat_name)

        if self.config.verbose:
            logger.info(
                "Z3: Added ATTTA motif avoidance (%d windows each direction)",
                total_nucs - 4,
            )

        return constraint_names

    # -------------------------------------------------------------------
    # T-run avoidance (no 6+ consecutive T bases)
    # -------------------------------------------------------------------

    def _encode_trun_constraints(
        self,
        optimizer: Optimize,
        nuc_exprs: list[ArithRef],
        max_run: int = 5,
    ) -> list[str]:
        """Encode T-run avoidance constraints.

        Forbids any window of (max_run + 1) consecutive T bases.
        Default max_run=5 means runs of 6+ consecutive T are forbidden.

        Args:
            optimizer: Z3 Optimize instance.
            nuc_exprs: Z3 expressions for each nucleotide position.
            max_run: Maximum allowed consecutive T count (default 5).

        Returns:
            List of constraint names added.
        """
        constraint_names: list[str] = []
        total_nucs = len(nuc_exprs)
        run_len = max_run + 1  # Forbidden run length

        # Check if mRNA stability constraint should be soft
        from .types import ConstraintType
        is_soft = ConstraintType.MRNA_STABILITY in self._soft_constraint_types
        soft_weight = 100  # Default soft constraint weight for T-runs

        first_trun_expr = None
        for j in range(total_nucs - run_len + 1):
            # Forbid all bases in window [j, j+run_len) being T (==3)
            all_t = And(*[nuc_exprs[j + k] == 3 for k in range(run_len)])
            trun_expr = Not(all_t)
            if is_soft:
                optimizer.add_soft(trun_expr, soft_weight)
            else:
                optimizer.add(trun_expr)
            if first_trun_expr is None:
                first_trun_expr = trun_expr

        trun_name = _sanitize_z3_name(f"no_trun_{run_len}")
        if first_trun_expr is not None:
            self._constraint_exprs.append((trun_name, first_trun_expr))
        constraint_names.append(trun_name)

        if self.config.verbose:
            logger.info(
                "Z3: Added T-run avoidance (max %d consecutive T, %d windows)",
                max_run, total_nucs - run_len + 1,
            )

        return constraint_names

    # -------------------------------------------------------------------
    # CpG dinucleotide avoidance
    # -------------------------------------------------------------------

    def _encode_cpg_constraints(
        self,
        optimizer: Optimize,
        codon_vars: list[ArithRef],
        codon_lists: list[list[str]],
    ) -> list[str]:
        """Encode CpG (CG) dinucleotide avoidance constraints.

        Within-codon: Forbid codon assignments containing the CG dinucleotide
        if CG-free alternatives exist.

        Cross-codon: Forbid adjacent codon pairs where codon i ends with C
        and codon i+1 starts with G (creating a cross-codon CG dinucleotide),
        when alternatives exist.

        This matches the OR-Tools engine CpG constraint encoding.

        Args:
            optimizer: Z3 Optimize instance.
            codon_vars: Z3 Int variables for codon choices.
            codon_lists: Synonymous codon lists per position.

        Returns:
            List of constraint names added.
        """
        constraint_names: list[str] = []

        # Check if CpG constraint should be soft
        from .types import ConstraintType
        is_soft = ConstraintType.NO_CPG in self._soft_constraint_types
        soft_weight = 50  # CpG is less critical than restriction sites

        # Within-codon: forbid codons containing CG
        for i, (c_var, codons_for_aa) in enumerate(zip(codon_vars, codon_lists)):
            cg_free = [c for c in codons_for_aa if "CG" not in c]
            if cg_free and len(cg_free) < len(codons_for_aa):
                # Some codons contain CG — restrict to CG-free alternatives
                cg_free_indices = [codons_for_aa.index(c) for c in cg_free]
                # At least one CG-free codon must be selected
                cpg_expr = Or(*[c_var == idx for idx in cg_free_indices])
                if is_soft:
                    optimizer.add_soft(cpg_expr, soft_weight)
                else:
                    optimizer.add(cpg_expr)
                cpg_name = _sanitize_z3_name(f"no_within_cpg_{i}")
                self._constraint_exprs.append((cpg_name, cpg_expr))
                constraint_names.append(cpg_name)

        # Cross-codon: forbid (C-ending codon, G-starting codon) pairs
        for i in range(len(codon_vars) - 1):
            domain_i = codon_lists[i]
            domain_j = codon_lists[i + 1]

            # Find codon pairs that create cross-codon CG
            c_ending = [k for k, c in enumerate(domain_i) if c[-1] == "C"]
            g_starting = [k for k, c in enumerate(domain_j) if c[0] == "G"]

            if c_ending and g_starting:
                # Only add if there are alternatives
                non_c_ending = [k for k, c in enumerate(domain_i) if c[-1] != "C"]
                non_g_starting = [k for k, c in enumerate(domain_j) if c[0] != "G"]

                # At least one of the two positions must avoid the CG boundary
                if non_c_ending or non_g_starting:
                    # If current codon ends with C, next codon must NOT start with G
                    rep_expr = None
                    for ci in c_ending:
                        for gi in g_starting:
                            cross_expr = Not(And(codon_vars[i] == ci, codon_vars[i + 1] == gi))
                            if is_soft:
                                optimizer.add_soft(cross_expr, soft_weight)
                            else:
                                optimizer.add(cross_expr)
                            if rep_expr is None:
                                rep_expr = cross_expr
                    cross_name = _sanitize_z3_name(f"no_cross_cpg_{i}")
                    if rep_expr is not None:
                        self._constraint_exprs.append((cross_name, rep_expr))
                    constraint_names.append(cross_name)

        if self.config.verbose:
            n_within = sum(1 for name in constraint_names if name.startswith("no_within"))
            n_cross = sum(1 for name in constraint_names if name.startswith("no_cross"))
            logger.info(
                "Z3: Added CpG avoidance constraints (%d within-codon, %d cross-codon)",
                n_within, n_cross,
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

        logger.info(
            "Z3 solution extraction: %d codons extracted, "
            "sequence_len=%d expected_len=%d",
            len(chosen_codons), len(sequence), len(protein) * 3,
        )

        # ── Validate solution ──────────────────────────────────────────
        if not sequence:
            logger.error("Z3: Solution extraction produced empty sequence")
            return SolverResult(
                sequence="",
                solved=False,
                backend_used=SolverBackend.Z3,
                solve_time_seconds=solve_time,
                num_constraints=num_constraints,
                num_variables=num_variables,
                fallback_used=True,
                warnings=["Z3 solution extraction produced empty sequence"],
            )

        # Verify translation fidelity (defense-in-depth)
        protein_seq = (
            self._protein_for_validation
            if hasattr(self, "_protein_for_validation")
            else ""
        )
        if protein_seq and len(sequence) == len(protein_seq) * 3:
            from ..constants import CODON_TABLE
            for i, expected_aa in enumerate(protein_seq):
                codon = sequence[i * 3 : i * 3 + 3]
                actual_aa = CODON_TABLE.get(codon)
                if actual_aa != expected_aa:
                    logger.warning(
                        "Z3: Codon %d (%s) translates to '%s', expected '%s'; "
                        "solution may be incorrect",
                        i, codon, actual_aa, expected_aa,
                    )

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

        -- Fix #1: Improved UNSAT diagnosis ----------------------

        Previously, the method replayed assertions from optimizer.assertions()
        which could include Z3-internal assertions that do not map to named
        constraints. Now we use the stored constraint expressions
        (self._constraint_exprs) which accurately track the Z3 expressions
        added for each named constraint.

        We create a Solver with tracked assertions (Implies(p_i, expr_i))
        and use Z3 unsat_core() to identify the minimal unsatisfiable
        subset, then map it back to human-readable constraint names.

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
            solver = Solver()
            # Enable unsat core tracking
            solver.set("unsat_core", True)

            tracked: list[object] = []

            # Use stored constraint expressions for UNSAT core tracking.
            # The constraint names are already sanitized (safe for Z3Bool).
            if self._constraint_exprs:
                for name, expr in self._constraint_exprs:
                    # Names are pre-sanitized by _sanitize_z3_name()
                    p = Z3Bool(f"p_{name}")
                    solver.add(Implies(p, expr))
                    tracked.append(p)
            else:
                # Fallback: replay from optimizer assertions
                assertions = optimizer.assertions()
                for i, assertion in enumerate(assertions):
                    if i < n_domain_assertions:
                        p = Z3Bool(f"domain_{i}")
                        solver.add(Implies(p, assertion))
                        tracked.append(p)
                        continue
                    name_idx = i - n_domain_assertions
                    name = constraint_names[name_idx] if name_idx < len(constraint_names) else f"constraint_{i}"
                    # Sanitize for Z3Bool safety
                    safe_name = _sanitize_z3_name(name)
                    p = Z3Bool(f"p_{safe_name}")
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
                            # Name format: p_<constraint_name>
                            # Extract the constraint name after "p_"
                            name = p_str[2:]
                            # If it's a numeric index (from fallback), map it
                            if name.isdigit():
                                name_idx = int(name)
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

                # Suggest relaxations: soft preferences first, then splice,
                # then restriction sites, then GC (GC is hardest to relax)
                priority_order = {"instability_motif": 1, "mrna_stability": 2,
                                  "cpg": 3,
                                  "splice_donor": 4, "splice_acceptor": 5,
                                  "restriction_site": 6, "gc": 7, "unknown": 8}
                sorted_violations = sorted(
                    violations,
                    key=lambda v: priority_order.get(
                        self._constraint_category(v.constraint_name), 8
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
                f"Core contains {len(violations)} constraint(s). "
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

        Handles both sanitized names (underscores for dots) and raw names.

        Args:
            name: Constraint name string (may be sanitized).

        Returns:
            Category string: 'gc', 'restriction_site', 'splice_donor',
            'splice_acceptor', 'instability_motif', 'mrna_stability',
            'cpg', or 'unknown'.
        """
        if name.startswith("gc_lo") or name.startswith("gc_hi"):
            return "gc"
        elif name.startswith("no_donor_splice"):
            return "splice_donor"
        elif name.startswith("no_acceptor_splice"):
            return "splice_acceptor"
        elif "ATTTA" in name or "TAAAT" in name or "attta" in name or "taaata" in name:
            return "instability_motif"
        elif name.startswith("no_trun"):
            return "mrna_stability"
        elif "cpg" in name.lower():
            return "cpg"
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
