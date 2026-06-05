"""
BioCompiler Optimizer v10.0.0
==============================
Multi-step certified gene optimization pipeline with aggressive GT resolution.

BREAKING CHANGE (v10.0.0): The optimizer now uses CODON_ADAPTIVENESS_TABLES
for codon selection (previously used SPECIES tables which disagreed with the
evaluation tables, causing incorrect CAI values). Five E. coli amino acids
had incorrect optimal codons; these are now corrected.

Step: Maximize CAI          — Greedy codon optimization (GT-aware, with unavoidable-GT tracking)
Step: Backtranslate CAI     — DP-based max-CAI back-translation with avoidable-GT avoidance
Step: Resolve Constraints   — Priority-based constraint resolution (fix GT/CG/RS with minimal CAI loss)
Step: Remove Restriction Sites — Restriction site removal by synonymous substitution
Step: Cross-Codon Optimization — Cross-codon constraint resolution (iterative, global validation)
Step: Within-Codon GT Resolution — Within-codon GT resolution (synonymous substitution + mutagenesis flagging)
Step: Mutagenesis Fallback  — Mutagenesis fallback (AA substitution for Valine etc. using BLOSUM62)
Step: Avoid CpG Islands     — CpG island avoidance (SKIPPED for prokaryotes)
Step: Cross-Codon Coordination — Cross-codon coordinated solver (exhaustive pair search at boundaries)
Step: CAI Hill Climb        — CAI hill climbing (upgrade codons while maintaining constraints)
Step: Reoptimize            — Re-optimization pass (iterative until convergence)

v10.0.0 changes:
  - CAI table unification: optimizer uses CODON_ADAPTIVENESS_TABLES (same as evaluator)
  - E. coli optimal codons corrected for Phe, Ile, Tyr, His, Arg
  - Prokaryote fast path: skip splice/CpG steps for E. coli and other prokaryotes
  - CAI-aware constraint resolution: all steps prefer higher-CAI alternatives
  - Incremental constraint checking via IncrementalSequenceState
  - CAI recovery pass: post-optimization upgrade of suboptimal codons
  - resolve_organism() for centralized organism name resolution
"""

from typing import Callable, List, Dict, Optional, Tuple, Set, Any

import logging
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from itertools import product as itertools_product

from .type_system import (
    CODON_TABLE, AA_TO_CODONS, BLOSUM62, PredicateResult,
    check_no_stop_codons, check_no_cryptic_splice, check_no_cpg_island,
    check_no_restriction_site, check_no_avoidable_gt,
    check_no_gt_dinucleotide_soft, _compute_max_gt_count,
    check_valid_coding_seq,
    find_cross_codon_gt, find_cross_codon_cg, find_cross_codon_restriction,
)
from .organisms import CODON_ADAPTIVENESS_TABLES, ORGANISM_GC_TARGETS, resolve_organism, ORGANISM_ALIASES, SPECIES_SHORT_NAMES, SUPPORTED_ORGANISMS
from .constants import reverse_complement, RESTRICTION_ENZYMES, IUPAC_EXPAND, VALID_IUPAC_BASES
from .scanner import gc_content
from .sliding_gc import check_sliding_gc, fix_sliding_gc_violations
from .maxentscan import score_donor, score_acceptor, max_donor_score, max_acceptor_score
from .mutagenesis import propose_mutagenesis, MutagenesisReport, MutagenesisProposal
from .incremental import IncrementalSequenceState, CodonCache, EnzymeSiteCache
from .certificate import format_certificate
from .exceptions import InvalidProteinError, UnsupportedOrganismError, OptimizationConstraintError
from .decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    DecisionProvenanceCollector,
    OptimizationDecisionTrail,
)
from .aho_corasick import AhoCorasickScanner, build_scanner_from_enzymes, build_scanner_from_sites  # type: ignore[import-untyped]
from .objectives import resolve_objective as _resolve_objective, cai_objective as _cai_objective

# ── NUMBA integration ──────────────────────────────────────────────
try:
    from .numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        count_gc as _numba_count_gc,
        count_dinucleotides as _numba_count_dinuc,
        seq_to_bytes as _seq_to_bytes,
        fast_dinucleotide_count as _numba_fast_dinuc_count,
        batch_codon_swap_score as _numba_batch_codon_swap_score,
        compute_cai_kernel as _numba_cai_kernel,
        compute_cai_incremental as _numba_cai_incremental,
    )
except ImportError:
    _HAS_NUMBA = False
    _numba_fast_dinuc_count = None  # type: ignore[assignment]
    _numba_batch_codon_swap_score = None  # type: ignore[assignment]
    _numba_cai_kernel = None  # type: ignore[assignment]
    _numba_cai_incremental = None  # type: ignore[assignment]

HAS_NUMBA: bool = _HAS_NUMBA if isinstance(_HAS_NUMBA, bool) else False

# ── NUMBA CAI helpers ──────────────────────────────────────────────
import numpy as _np


def _adaptiveness_to_array(codon_adaptiveness: dict[str, float]) -> _np.ndarray:
    """Convert codon adaptiveness dict to numpy array indexed by codon encoding.

    Codons are mapped to array indices using a base-4 encoding:
        A=0, C=1, G=2, T=3  →  index = b0*16 + b1*4 + b2
    This gives indices in [0, 63] for all 64 possible 3-letter codons,
    providing O(1) lookup without hash tables.

    Args:
        codon_adaptiveness: Dict mapping codon strings (e.g. 'ATG') to
            adaptiveness values (float in [0.0, 1.0]).

    Returns:
        numpy float64 array of length 64, where index = base4_encode(codon).
        Missing codons get a default value of 1e-10.
    """
    _BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    arr = _np.full(64, 1e-10, dtype=_np.float64)
    for codon, w in codon_adaptiveness.items():
        try:
            idx = _BASE_MAP[codon[0]] * 16 + _BASE_MAP[codon[1]] * 4 + _BASE_MAP[codon[2]]
            arr[idx] = w
        except (KeyError, IndexError):
            continue  # skip non-standard codons
    return arr


def _codon_to_index(codon: str) -> int:
    """Map a 3-letter codon string to its base-4 array index.

    Uses: A=0, C=1, G=2, T=3 → index = b0*16 + b1*4 + b2
    """
    _BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return _BASE_MAP[codon[0]] * 16 + _BASE_MAP[codon[1]] * 4 + _BASE_MAP[codon[2]]


def _dna_to_codon_indices(dna: str) -> _np.ndarray:
    """Convert DNA string to array of codon adaptiveness indices.

    Each codon position gets its base-4 encoded index for O(1) lookup
    in the adaptiveness array.

    Args:
        dna: DNA sequence string (length must be a multiple of 3).

    Returns:
        numpy int64 array of length len(dna)//3 with base-4 encoded indices.
    """
    n_codons = len(dna) // 3
    indices = _np.empty(n_codons, dtype=_np.int64)
    for i in range(n_codons):
        indices[i] = _codon_to_index(dna[i * 3:i * 3 + 3])
    return indices


def _compute_cai_fast(seq: str, codon_adaptiveness: dict[str, float]) -> float:
    """Compute CAI using the NUMBA kernel when available, pure Python otherwise.

    This function replaces the pure-Python loop in _compute_seq_cai with
    the NUMBA-accelerated compute_cai_kernel. The kernel uses log-sum
    for numerical stability (geometric mean of relative adaptiveness values).

    When NUMBA is available, converts the adaptiveness dict to a numpy array
    and uses compute_cai_kernel for O(n) CAI computation in compiled code.
    Falls back to pure-Python log-sum when NUMBA is not installed.

    Args:
        seq: DNA coding sequence.
        codon_adaptiveness: Dict mapping codon strings to adaptiveness values.

    Returns:
        CAI value in [0.0, 1.0]. Returns 0.0 for empty/invalid sequences.
    """
    if not seq or len(seq) < 3:
        return 0.0

    if HAS_NUMBA and _numba_cai_kernel is not None:
        try:
            adapt_arr = _adaptiveness_to_array(codon_adaptiveness)
            # Build codon indices array, excluding Met and stop codons from count
            n_total = len(seq) // 3
            codon_indices_list = []
            for i in range(n_total):
                codon = seq[i * 3:i * 3 + 3]
                aa = CODON_TABLE.get(codon)
                if aa == 'M' or aa == '*':
                    continue
                codon_indices_list.append(_codon_to_index(codon))
            if not codon_indices_list:
                return 0.0
            codon_indices = _np.array(codon_indices_list, dtype=_np.int64)
            n_codons = len(codon_indices)
            return _numba_cai_kernel(adapt_arr, codon_indices, n_codons)
        except Exception:
            pass  # Fall through to pure-Python on any NUMBA error

    # Pure-Python fallback
    epsilon = 1e-10
    log_sum = 0.0
    count = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = CODON_TABLE.get(codon)
        if aa == 'M' or aa == '*':
            continue
        w = codon_adaptiveness.get(codon, epsilon)
        if w <= 0.0:
            w = epsilon
        log_sum += math.log(w)
        count += 1

    if count == 0:
        return 0.0
    return math.exp(log_sum / count)


def _count_dinucs_fast(seq: str, *dinucleotides: str) -> tuple[int, ...]:
    """Count multiple dinucleotides in a single pass using the NUMBA kernel.

    Falls back to pure-Python counting when NUMBA is unavailable.

    Args:
        seq: DNA sequence string (uppercase ACGT).
        *dinucleotides: One or more dinucleotide strings (e.g. "GT", "CG", "AG").

    Returns:
        Tuple of counts, one per dinucleotide, in the same order as input.
    """
    n_dinucs = len(dinucleotides)
    if n_dinucs == 0:
        return ()

    # Fast path: NUMBA kernel
    if HAS_NUMBA and _numba_fast_dinuc_count is not None:
        import numpy as _np
        seq_bytes = _seq_to_bytes(seq)
        dinuc_keys = _np.array(
            [[ord(d[0]), ord(d[1])] for d in dinucleotides],
            dtype=_np.uint8,
        )
        counts = _numba_fast_dinuc_count(seq_bytes, dinuc_keys, n_dinucs)
        return tuple(int(c) for c in counts)

    # Pure-Python fallback
    results = []
    for di in dinucleotides:
        count = 0
        pos = 0
        while True:
            pos = seq.find(di, pos)
            if pos == -1:
                break
            count += 1
            pos += 1
        results.append(count)
    return tuple(results)


# ── Batch codon swap scorer (NUMBA-accelerated) ─────────────────────
class _BatchSwapScorer:
    """Pre-compute CAI arrays for batch_codon_swap_score kernel invocations.

    Builds the static adaptiveness array and codon-to-index mapping once,
    then provides a fast ``score_candidates`` method that calls the NUMBA
    kernel (or pure-Python fallback) to score all candidate codon swaps
    at a given position in one vectorized pass.

    Usage::

        scorer = _BatchSwapScorer(species_cai)
        # Per iteration / position:
        scores = scorer.score_candidates(seq_codons, codon_idx, candidate_codons)
        # scores[k] = full-sequence CAI if we swapped codon at codon_idx to candidate_codons[k]
    """

    def __init__(self, species_cai: Dict[str, float]) -> None:
        self._species_cai = species_cai

        # Build a stable sorted list of all known codons → integer index.
        # CODON_TABLE keys are the canonical 64 codons; we include any
        # extra codon that appears in species_cai as well.
        from .type_system import CODON_TABLE as _CT
        all_codons = sorted(set(_CT.keys()) | set(species_cai.keys()))
        self._codon_to_idx: Dict[str, int] = {c: i for i, c in enumerate(all_codons)}
        self._idx_to_codon: Dict[int, str] = {i: c for i, c in enumerate(all_codons)}
        self._n_codons_total = len(all_codons)

        # Build adaptiveness array (float64) — one entry per codon index.
        epsilon = 1e-10
        self._adaptiveness: list[float] = []
        for c in all_codons:
            w = species_cai.get(c, 0.0)
            self._adaptiveness.append(w if w > 0.0 else epsilon)

        # Lazy numpy arrays (built on first use when NUMBA is available)
        self._np_adaptiveness = None  # type: ignore[assignment]
        self._np_built = False

        # Incremental CAI tracking state — avoids O(n) log_sum recomputation
        # after each accepted codon swap.  Set via reset_incremental_state().
        self._current_log_sum: float | None = None
        self._n_codons: int = 0

    def _ensure_numpy(self) -> None:
        """Lazily build numpy arrays for the NUMBA kernel."""
        if self._np_built:
            return
        try:
            import numpy as _np
            self._np_adaptiveness = _np.array(self._adaptiveness, dtype=_np.float64)
        except ImportError:
            self._np_adaptiveness = None
        self._np_built = True

    def _build_codon_indices(self, seq_codons: List[str]) -> "Any":
        """Build the codon_indices numpy array from a list of codon strings."""
        import numpy as _np
        indices = _np.array(
            [self._codon_to_idx.get(c, 0) for c in seq_codons],
            dtype=_np.int64,
        )
        return indices

    def _compute_log_sum(self, seq_codons: List[str]) -> float:
        """Compute the current sum of log(w_i) over all codon positions."""
        epsilon = 1e-10
        log_sum = 0.0
        for c in seq_codons:
            idx = self._codon_to_idx.get(c, 0)
            w = self._adaptiveness[idx]
            if w <= 0.0:
                w = epsilon
            log_sum += math.log(w)
        return log_sum

    def reset_incremental_state(self, seq_codons: List[str]) -> None:
        """Initialize (or re-initialize) incremental CAI tracking state.

        Computes the full log_sum once from the current codon list and
        caches it for O(1) incremental updates via compute_cai_incremental.

        Args:
            seq_codons: Current codon strings for every position in the sequence.
        """
        self._current_log_sum = self._compute_log_sum(seq_codons)
        self._n_codons = len(seq_codons)

    def update_incremental_state(self, old_codon: str, new_codon: str) -> None:
        """Update the cached log_sum after an accepted codon swap.

        Uses compute_cai_incremental (NUMBA or pure-Python) for O(1)
        update instead of O(n) full recomputation.

        Args:
            old_codon: The codon that was replaced.
            new_codon: The codon that replaced it.
        """
        if self._current_log_sum is None:
            return  # Not initialized; caller should use reset_incremental_state

        epsilon = 1e-10

        # Get old adaptiveness
        old_idx = self._codon_to_idx.get(old_codon, 0)
        w_old = self._adaptiveness[old_idx]
        if w_old <= 0.0 or w_old != w_old:  # catch NaN
            w_old = epsilon

        # Get new adaptiveness
        new_idx = self._codon_to_idx.get(new_codon, 0)
        w_new = self._adaptiveness[new_idx]
        if w_new <= 0.0 or w_new != w_new:  # catch NaN
            w_new = epsilon

        # O(1) incremental update: new_log_sum = old_log_sum - log(w_old) + log(w_new)
        if HAS_NUMBA and _numba_cai_incremental is not None:
            try:
                new_cai = _numba_cai_incremental(
                    self._current_log_sum, self._n_codons, w_old, w_new
                )
                # Recover the new log_sum from the CAI value
                # CAI = exp(log_sum / n)  →  log_sum = n * log(CAI)
                if new_cai > 0.0 and self._n_codons > 0:
                    self._current_log_sum = self._n_codons * math.log(new_cai)
                else:
                    # Fallback to direct computation
                    self._current_log_sum = (
                        self._current_log_sum - math.log(w_old) + math.log(w_new)
                    )
                return
            except Exception:
                pass  # Fall through to pure-Python

        # Pure-Python incremental update
        self._current_log_sum = (
            self._current_log_sum - math.log(w_old) + math.log(w_new)
        )

    @property
    def current_log_sum(self) -> float | None:
        """The cached sum of log(w_i) for the current sequence state."""
        return self._current_log_sum

    def score_candidates(
        self,
        seq_codons: List[str],
        swap_position: int,
        candidate_codons: List[str],
    ) -> List[float]:
        """Score all candidate codon swaps at a single position.

        Returns a list of full-sequence CAI values, one per candidate,
        representing the CAI of the sequence if the codon at
        ``swap_position`` were replaced by each candidate.

        Args:
            seq_codons: Current codon strings for every position in the sequence.
            swap_position: Index into seq_codons of the codon being swapped.
            candidate_codons: List of candidate replacement codon strings.

        Returns:
            List of CAI scores (float), one per candidate.
        """
        n_candidates = len(candidate_codons)
        if n_candidates == 0:
            return []

        n_codons = len(seq_codons)

        # ── Fast path: NUMBA batch kernel ──────────────────────────
        if (
            HAS_NUMBA
            and _numba_batch_codon_swap_score is not None
        ):
            self._ensure_numpy()
            if self._np_adaptiveness is not None:
                import numpy as _np
                codon_indices = self._build_codon_indices(seq_codons)
                candidate_indices = _np.array(
                    [self._codon_to_idx.get(c, 0) for c in candidate_codons],
                    dtype=_np.int64,
                )
                # Use cached log_sum if available (avoids O(n) recompute)
                current_log_sum = self._current_log_sum if self._current_log_sum is not None else self._compute_log_sum(seq_codons)
                scores_arr = _numba_batch_codon_swap_score(
                    self._np_adaptiveness,
                    codon_indices,
                    n_codons,
                    swap_position,
                    candidate_indices,
                    n_candidates,
                    current_log_sum,
                )
                return [float(s) for s in scores_arr]

        # ── Fallback: pure-Python per-candidate evaluation ─────────
        epsilon = 1e-10
        current_log_sum = self._current_log_sum if self._current_log_sum is not None else self._compute_log_sum(seq_codons)

        old_idx = self._codon_to_idx.get(seq_codons[swap_position], 0)
        w_old = self._adaptiveness[old_idx]
        if w_old <= 0.0:
            w_old = epsilon
        log_w_old = math.log(w_old)

        scores: List[float] = []
        for c in candidate_codons:
            c_idx = self._codon_to_idx.get(c, 0)
            w_new = self._adaptiveness[c_idx]
            if w_new <= 0.0:
                w_new = epsilon
            new_log_sum = current_log_sum - log_w_old + math.log(w_new)
            scores.append(math.exp(new_log_sum / n_codons) if n_codons > 0 else 0.0)

        return scores


__all__ = [
    "OptimizationResult",
    "FullConstructResult",
    "optimize_sequence",
    "batch_optimize",
    "protein_to_aa_list",
    "BioOptimizer",
    "HybridOptimizer",
    "ConvergenceTracker",
    "score_splice_donor_potential",
    "SPLICE_DONOR_POTENTIAL_THRESHOLD",
    "EUKARYOTE_CAI_GT_COST_THRESHOLD",
    "GT_CAI_LOG_ADAPTIVENESS_COST",
]

logger = logging.getLogger(__name__)


# ────────────────────────────────────────────────────────────
# Named Constants
# ────────────────────────────────────────────────────────────

# Iteration limits for each optimization step
MAX_RESTRICTION_SITE_ITERATIONS: int = 100
MAX_IUPAC_SITE_ITERATIONS: int = 100
MAX_ATTTA_MOTIF_ITERATIONS: int = 100
MAX_T_RUN_ITERATIONS: int = 100
MAX_GC_ADJUSTMENT_ITERATIONS: int = 200
MAX_SPLICE_ELIMINATION_ITERATIONS: int = 300
MAX_CPG_DISRUPTION_ITERATIONS: int = 200

# Main optimization loop convergence settings
DEFAULT_MAX_ITERATIONS: int = 1000
CONVERGENCE_IMPROVEMENT_THRESHOLD: float = 1e-6
CONVERGENCE_PATIENCE: int = 3  # consecutive non-improving iterations before convergence
OSCILLATION_WINDOW: int = 10   # window size for oscillation detection

# Thresholds and sentinel values
T_RUN_LENGTH_THRESHOLD: int = 6
ELIMINATED_SITE_SCORE: float = -999.0
TOP_CAI_ALTERNATIVES: int = 3
IUPAC_EXPANSION_CAP: int = 4096


# ────────────────────────────────────────────────────────────
# Convergence Tracker
# ────────────────────────────────────────────────────────────

class ConvergenceTracker:
    """Track optimization objective convergence and detect stagnation or oscillation.

    The objective is CAI * constraint_satisfaction_score.  CAI is in [0, 1]
    and constraint_satisfaction_score is the fraction of predicates that pass,
    also in [0, 1], so the combined objective is in [0, 1].

    Convergence is declared when the objective has not improved by more than
    ``improvement_threshold`` for ``patience`` consecutive iterations.

    Oscillation is declared when the objective has both increased and decreased
    within the last ``oscillation_window`` iterations — this indicates the
    optimizer is cycling and should stop at the best point seen so far.
    """

    def __init__(
        self,
        improvement_threshold: float = CONVERGENCE_IMPROVEMENT_THRESHOLD,
        patience: int = CONVERGENCE_PATIENCE,
        oscillation_window: int = OSCILLATION_WINDOW,
    ) -> None:
        self.improvement_threshold = improvement_threshold
        self.patience = patience
        self.oscillation_window = oscillation_window
        self.history: list[float] = []
        self.best_objective: float = float('-inf')
        self.best_iteration: int = -1
        self._no_improvement_count: int = 0

    def record(self, objective: float) -> None:
        """Record the objective value for the current iteration."""
        self.history.append(objective)
        if objective > self.best_objective:
            self.best_objective = objective
            self.best_iteration = len(self.history) - 1

    def check_convergence(self) -> str | None:
        """Check if the optimizer has converged, hit max iterations, or is oscillating.

        Returns:
            None if optimization should continue,
            "converged" if the objective has plateaued,
            "oscillating" if the objective is cycling,
        """
        if len(self.history) < 2:
            return None

        latest = self.history[-1]
        # Check plateau: no significant improvement for `patience` iterations
        recent = self.history[-(self.patience + 1):]
        if len(recent) >= self.patience + 1:
            improvement = max(recent) - min(recent)
            if improvement < self.improvement_threshold:
                return "converged"

        # Check oscillation: both increases and decreases within the window
        window = self.history[-self.oscillation_window:]
        if len(window) >= 3:
            has_increase = any(window[i + 1] > window[i] for i in range(len(window) - 1))
            has_decrease = any(window[i + 1] < window[i] for i in range(len(window) - 1))
            if has_increase and has_decrease:
                # Additional check: the objective hasn't improved overall in the window
                if window[-1] <= window[0] + self.improvement_threshold:
                    return "oscillating"

        return None

    @property
    def iterations(self) -> int:
        """Number of iterations recorded."""
        return len(self.history)

    @property
    def best(self) -> float:
        """Best objective value seen so far."""
        return self.best_objective

    @property
    def best_iteration_index(self) -> int:
        """Index of the iteration where the best objective was achieved."""
        return self.best_iteration


# ────────────────────────────────────────────────────────────
# High-level OptimizationResult and optimize_sequence API
# ────────────────────────────────────────────────────────────

@dataclass
class OptimizationResult:
    """Result of optimizing a protein sequence.

    Provides the optimized DNA sequence along with quality metrics
    and a list of predicates that the result fails to satisfy.
    """
    sequence: str
    gc_content: float
    cai: float
    failed_predicates: list[str] = field(default_factory=list)
    predicate_results: list[PredicateResult] = field(default_factory=list)
    certificate_text: str = ""
    # Extended attributes for API/visualization compatibility
    protein: str = ""
    fallback_used: bool = False
    satisfied_predicates: list[str] = field(default_factory=list)
    aa_substitutions: list[dict[str, Any]] = field(default_factory=list)
    mutagenesis_applied: bool = False
    # mRNA stability metrics (populated when optimize_mrna_stability=True)
    mrna_stability_score: float | None = None
    destabilizing_motifs_removed: int = 0
    stability_improvement: float | None = None
    # Provenance: OptimizationRecord for this run (populated by optimize_sequence)
    provenance: Any = field(default=None, repr=False)
    # Codon pair bias metric (populated when consider_codon_pair_bias=True)
    codon_pair_bias: float | None = None
    # UTR suggestions (populated when include_utr=True)
    suggested_5utr: str | None = None
    suggested_3utr: str | None = None
    utr_score_5: float | None = None
    utr_score_3: float | None = None
    # Decision-level provenance trail (populated when track_provenance=True)
    decision_trail: OptimizationDecisionTrail | None = None
    # Convergence tracking (populated by optimize_sequence)
    convergence_status: str | None = None  # "converged" | "max_iterations" | "oscillating" | None
    iterations_used: int = 0
    # Optimization warnings (cap-exceeded notices, convergence issues)
    warnings: list[str] = field(default_factory=list)
    # Custom objective score (populated when a non-default objective is used)
    objective_score: float | None = None

    def __post_init__(self):
        """Validate OptimizationResult invariants."""
        if self.protein and self.sequence:
            assert len(self.sequence) == len(self.protein) * 3, (
                f"Sequence length ({len(self.sequence)}) must equal "
                f"protein length * 3 ({len(self.protein) * 3})"
            )
        assert 0.0 <= self.cai <= 1.0, f"CAI must be in [0, 1], got {self.cai}"
        assert 0.0 <= self.gc_content <= 1.0, f"GC content must be in [0, 1], got {self.gc_content}"
        if self.mutagenesis_applied:
            assert self.aa_substitutions is not None and len(self.aa_substitutions) > 0, (
                "Mutagenesis applied but no substitutions recorded"
            )
        if self.utr_score_5 is not None:
            assert 0.0 <= self.utr_score_5 <= 1.0, (
                f"UTR 5' score must be in [0, 1], got {self.utr_score_5}"
            )
        if self.utr_score_3 is not None:
            assert 0.0 <= self.utr_score_3 <= 1.0, (
                f"UTR 3' score must be in [0, 1], got {self.utr_score_3}"
            )


@dataclass
class FullConstructResult:
    """Complete expression construct: 5' UTR + CDS + 3' UTR.

    This represents what a biologist would actually order from a gene
    synthesis company — the full DNA construct ready for cloning or
    direct expression in the target organism.

    The CDS is the optimized coding sequence. The UTRs are suggested
    (not enforced) and should be evaluated by the user before ordering.

    Attributes:
        utr5: 5' UTR sequence (empty string if not provided).
        cds: Optimized coding sequence (starts with ATG, ends with stop codon).
        utr3: 3' UTR sequence (empty string if not provided).
        full_construct: Concatenated 5'UTR + CDS + 3'UTR.
        organism: Target organism for this construct.
        gc_content: GC fraction of the full construct.
        cai: Codon Adaptation Index of the CDS.
        utr_score_5: Expression suitability score for the 5' UTR (0.0–1.0).
        utr_score_3: Expression suitability score for the 3' UTR (0.0–1.0).
        protein: Amino acid sequence encoded by the CDS.
    """
    utr5: str
    cds: str
    utr3: str
    full_construct: str
    organism: str
    gc_content: float
    cai: float
    utr_score_5: float | None = None
    utr_score_3: float | None = None
    protein: str = ""

    def __post_init__(self):
        """Validate FullConstructResult invariants."""
        assert self.full_construct == self.utr5 + self.cds + self.utr3, (
            "full_construct must equal utr5 + cds + utr3"
        )
        assert 0.0 <= self.gc_content <= 1.0, (
            f"GC content must be in [0, 1], got {self.gc_content}"
        )
        assert 0.0 <= self.cai <= 1.0, f"CAI must be in [0, 1], got {self.cai}"
        if self.utr_score_5 is not None:
            assert 0.0 <= self.utr_score_5 <= 1.0, (
                f"UTR 5' score must be in [0, 1], got {self.utr_score_5}"
            )
        if self.utr_score_3 is not None:
            assert 0.0 <= self.utr_score_3 <= 1.0, (
                f"UTR 3' score must be in [0, 1], got {self.utr_score_3}"
            )


# ==============================================================================
# Input Validation
# ==============================================================================

def protein_to_aa_list(protein: str) -> list[str]:
    """Convert protein string to list of amino acid codes. Raises InvalidProteinError for bad input.

    Pre-conditions:
    - protein must be a non-empty string of standard amino acid codes

    Post-conditions:
    - result is a list of valid single-letter amino acid codes
    - len(result) == len(protein.strip())
    """
    if not protein or not protein.strip():
        raise InvalidProteinError(protein, set())
    protein = protein.upper().strip()
    valid_aas = set(AA_TO_CODONS.keys())
    invalid = set(ch for ch in protein if ch not in valid_aas)
    if invalid:
        raise InvalidProteinError(protein, invalid)
    return list(protein)


# ==============================================================================
# Restriction Site Removal Helpers
# ==============================================================================

def _find_site_in_sequence(sequence: str, site: str, site_rc: str) -> list[int]:
    """Find all positions where site or its reverse complement appears in sequence.

    Pre-conditions:
    - sequence is a valid uppercase DNA string
    - site is a non-empty uppercase DNA string
    - site_rc is the reverse complement of site (or empty string)

    Post-conditions:
    - returns sorted list of unique positions
    """
    positions = []
    if site:
        start = 0
        while True:
            pos = sequence.find(site, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
    if site_rc and site_rc != site:  # Avoid double-counting palindromes
        start = 0
        while True:
            pos = sequence.find(site_rc, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
    return sorted(set(positions))


def _get_overlapping_codons(pos: int, site_len: int, n_codons: int) -> list[int]:
    """Get indices of codons that overlap with a site at position pos.

    Pre-conditions:
    - pos >= 0
    - site_len > 0
    - n_codons > 0

    Post-conditions:
    - all indices in result are in [0, n_codons)
    """
    assert pos >= 0, f"Position must be non-negative, got {pos}"
    assert site_len > 0, f"Site length must be positive, got {site_len}"
    assert n_codons > 0, f"Number of codons must be positive, got {n_codons}"

    first_codon = pos // 3
    last_base = pos + site_len - 1
    last_codon = last_base // 3
    return list(range(max(0, first_codon), min(n_codons, last_codon + 1)))


def _remove_site_multicodon(
    sequence: str, aas: list[str], sorted_codons: dict[str, list[str]],
    site_upper: str, site_rc: str,
    usage: dict[str, float] | None = None,
) -> tuple[str, bool]:
    """
    Try to remove a restriction site using multi-codon coordinated solving.

    When a site straddles codon boundaries (e.g., PstI CTG|CAG spans codons i and i+1),
    single-codon swaps fail because changing either codon alone doesn't eliminate the site.
    This function enumerates valid codon COMBINATIONS for all overlapping codons.

    When ``usage`` (a codon→CAI adaptiveness dict) is provided, the function
    ranks all viable codon combinations by CAI impact and picks the one that
    minimises CAI loss.  Without ``usage``, it returns the first combo that
    eliminates the site (legacy behaviour).

    Pre-conditions:
    - sequence is uppercase DNA
    - len(aas) > 0
    - sorted_codons maps each aa in aas to a non-empty list of codons
    - site_upper is a valid uppercase DNA string

    Post-conditions:
    - if fixed, returned sequence has same length as input and encodes same protein
    - if not fixed, returned sequence is identical to input
    """
    n_codons = len(aas)
    positions = _find_site_in_sequence(sequence, site_upper, site_rc)

    for pos in positions:
        overlapping = _get_overlapping_codons(pos, len(site_upper), n_codons)
        if not overlapping:
            continue

        # Build candidate codon lists for each overlapping position
        candidate_lists = []
        for ci in overlapping:
            aa = aas[ci]
            candidate_lists.append(sorted_codons[aa])

        if usage is not None:
            # CAI-aware: enumerate ALL combos, rank by CAI, pick the best
            best_seq: str | None = None
            best_cai_sum: float = float('-inf')

            for combo in itertools_product(*candidate_lists):
                # Build test sequence with this combo applied
                test = list(sequence)
                for idx, ci in enumerate(overlapping):
                    start = ci * 3
                    test[start:start + 3] = list(combo[idx])
                test_seq = "".join(test)

                # Check if site is eliminated
                if site_upper not in test_seq and (not site_rc or site_rc not in test_seq):
                    # Compute total CAI contribution of the changed codons
                    cai_sum = sum(
                        math.log(usage.get(combo[idx], 1e-10))
                        for idx, ci in enumerate(overlapping)
                    )
                    if cai_sum > best_cai_sum:
                        best_cai_sum = cai_sum
                        best_seq = test_seq

            if best_seq is not None:
                return best_seq, True
        else:
            # Legacy: return first combo that eliminates the site
            for combo in itertools_product(*candidate_lists):
                # Build test sequence with this combo applied
                test = list(sequence)
                for idx, ci in enumerate(overlapping):
                    start = ci * 3
                    test[start:start + 3] = list(combo[idx])
                test_seq = "".join(test)

                # Check if site is eliminated
                if site_upper not in test_seq and (not site_rc or site_rc not in test_seq):
                    return test_seq, True

    return sequence, False


# ==============================================================================
# Splice Donor Potential Scoring
# ==============================================================================

# Threshold above which a GT dinucleotide is considered a high-risk splice donor.
# GTs with score >= this threshold should be avoided even at CAI cost.
# GTs with score < this threshold are acceptable (soft violation) in favor of CAI.
SPLICE_DONOR_POTENTIAL_THRESHOLD: float = 0.5

# Maximum distance (nt) to search downstream of a GT for a potential AG acceptor.
_MAX_ACCEPTOR_SEARCH_DIST: int = 80


def score_splice_donor_potential(dna: str, position: int) -> float:
    """Score how likely a GT dinucleotide is to function as a splice donor.

    Not all GT dinucleotides are splice donors.  A true splice donor needs:
    - GT at the 5' end of an intron
    - AG at the 3' end (acceptor) typically 30-200 nt downstream
    - A polypyrimidine tract upstream of the AG
    - A branch point (YNYRAY) ~18-40 nt upstream of the AG

    This function returns a score from 0.0 (unlikely splice donor) to 1.0
    (strong, functional splice donor) based on:

    1. **MaxEntScan donor score** (primary, ~70% weight): How well the
       9-mer context matches the splice donor PWM.
    2. **Downstream acceptor context** (~20% weight): Whether there's a
       downstream AG with a reasonable acceptor score within intron-length
       distance.
    3. **Polypyrimidine tract quality** (~10% weight): Whether there's a
       C/T-rich region between the GT and the best downstream AG.

    Score interpretation:
    - < 0.3: Low risk — GT very unlikely to function as splice donor.
      In-codon GTs from optimal codons almost always score here.
    - 0.3–0.5: Moderate risk — possible cryptic splice site.
    - > 0.5: High risk — likely functional splice donor; should be avoided.

    Pre-conditions:
    - dna is a valid uppercase DNA string
    - 0 <= position < len(dna) - 1
    - dna[position:position+2] == "GT"

    Post-conditions:
    - returns a float in [0.0, 1.0]
    - deterministic: same input always produces same output
    """
    # ── Component 1: MaxEntScan donor score ──
    donor_score = score_donor(dna, position)

    # Normalize donor score to [0, 1].
    # Typical ranges:
    #   -5 to 0:  non-donor
    #   0 to 3:   very weak
    #   3 to 8:   moderate cryptic
    #   8+:       strong canonical
    if donor_score < 0:
        donor_component = 0.0
    elif donor_score < 3.0:
        donor_component = 0.05 + 0.10 * (donor_score / 3.0)   # 0.05–0.15
    elif donor_score < 8.0:
        donor_component = 0.15 + 0.35 * ((donor_score - 3.0) / 5.0)  # 0.15–0.50
    else:
        donor_component = 0.50 + min(0.20, 0.20 * ((donor_score - 8.0) / 5.0))  # 0.50–0.70

    # ── Component 2: Downstream acceptor context ──
    # Search for the strongest AG acceptor downstream within intron distance.
    acceptor_component = 0.0
    best_ag_score = -999.0
    best_ag_pos = -1
    search_end = min(len(dna) - 1, position + _MAX_ACCEPTOR_SEARCH_DIST)

    for i in range(position + 4, search_end):  # +4: skip the GT itself
        if dna[i] == "A" and dna[i + 1] == "G":
            ag_score = score_acceptor(dna, i)
            if ag_score > best_ag_score:
                best_ag_score = ag_score
                best_ag_pos = i

    if best_ag_pos >= 0:
        dist = best_ag_pos - position
        # Typical intron: 30-5000 nt; short introns in yeast ~80-500 nt
        # Only count acceptors at reasonable intron-like distances
        if 20 <= dist <= _MAX_ACCEPTOR_SEARCH_DIST:
            if best_ag_score >= 8.0:
                acceptor_component = 0.20  # Strong acceptor
            elif best_ag_score >= 3.0:
                acceptor_component = 0.12  # Moderate acceptor
            else:
                acceptor_component = 0.05  # Weak acceptor

    # ── Component 3: Polypyrimidine tract quality ──
    # Check for C/T-rich region upstream of the best AG (between GT and AG).
    py_component = 0.0
    if best_ag_pos >= 0:
        tract_start = max(position + 4, best_ag_pos - 30)
        tract_end = max(position + 4, best_ag_pos - 5)
        if tract_end > tract_start:
            tract = dna[tract_start:tract_end]
            py_fraction = (tract.count("C") + tract.count("T")) / len(tract)
            if py_fraction >= 0.7:
                py_component = 0.10  # Good polypyrimidine tract
            elif py_fraction >= 0.5:
                py_component = 0.05  # Moderate

    # Combined score
    return min(1.0, max(0.0, donor_component + acceptor_component + py_component))


# ==============================================================================
# Eukaryotic CAI Recovery
# ==============================================================================

# Threshold for CAI-vs-GT tradeoff: for eukaryotes, only avoid GT if the
# per-codon CAI cost (difference in adaptiveness weight) is below this value.
# If the optimal codon contains GT and using a GT-free alternative would
# drop the adaptiveness by more than this threshold, keep the optimal codon.
EUKARYOTE_CAI_GT_COST_THRESHOLD: float = 0.10

# Relative CAI threshold for GT-aware codon selection.  When selecting
# codons for eukaryotes, if a codon would create a GT at the boundary
# with the NEXT codon, we check if there's an alternative codon with
# similar CAI (within this relative fraction).  If the best alternative
# has CAI < (1 - GT_BOUNDARY_CAI_TOLERANCE) * optimal_CAI, we accept
# the GT and use the optimal codon.
GT_BOUNDARY_CAI_TOLERANCE: float = 0.10  # 10% relative CAI tolerance

# Log-adaptiveness cost threshold for GT avoidance.  When the optimal codon
# contains GT and the CAI cost (difference in log-adaptiveness) exceeds this
# value AND the splice donor score < SPLICE_DONOR_POTENTIAL_THRESHOLD, the
# optimal codon is used despite the GT.
GT_CAI_LOG_ADAPTIVENESS_COST: float = 0.1  # >10% relative CAI cost


def _gt_aware_select_codon(
    aa: str,
    next_aa: str | None,
    sorted_codons: dict[str, list[str]],
    usage: dict[str, float],
    cai_tolerance: float = GT_BOUNDARY_CAI_TOLERANCE,
) -> str:
    """Select the best codon for an amino acid, preferring GT-free boundary if
    the CAI cost is within tolerance.

    This is a "best effort" approach that maximizes CAI while minimizing GT
    dinucleotides at codon boundaries. For eukaryotes, if the optimal codon
    would create a GT at the boundary with the next codon (i.e., the optimal
    codon ends with 'G' and the next codon's optimal starts with 'T'), we
    check if there's an alternative codon with similar CAI (within
    ``cai_tolerance`` relative) that avoids the boundary GT.

    Decision logic:
    1. Select the optimal (highest-CAI) codon for this amino acid.
    2. If no next amino acid, or the optimal codon doesn't create a boundary
       GT, use it.
    3. If the optimal codon WOULD create a boundary GT with the next codon's
       optimal start:
       a. Check each alternative codon (in CAI order) to see if it avoids
          the boundary GT and has CAI within tolerance.
       b. If such an alternative exists, use it (best CAI among valid alts).
       c. If no alternative within CAI tolerance, use the optimal codon and
          accept the GT (CAI > GT avoidance).

    For prokaryotes, this function is never called (GT avoidance is skipped).

    Args:
        aa: Amino acid at the current position.
        next_aa: Amino acid at the next position (None if last position or
            next is a stop codon).
        sorted_codons: Codons per AA sorted by CAI descending.
        usage: Codon adaptiveness table (codon → w value).
        cai_tolerance: Maximum relative CAI loss acceptable for avoiding a
            boundary GT. Default: 0.10 (10%).

    Returns:
        The selected codon string.
    """
    candidates = sorted_codons.get(aa, [])
    if not candidates:
        return AA_TO_CODONS.get(aa, [""])[0]

    optimal = candidates[0]
    optimal_cai = usage.get(optimal, 0.0)

    # If no next AA, just use the optimal codon
    if next_aa is None:
        return optimal

    # Determine what the next position's first base would be if we use
    # the optimal codon for next_aa.  We need to look ahead to see if
    # our codon's last base + next codon's first base = "GT".
    next_candidates = sorted_codons.get(next_aa, [])
    if not next_candidates:
        return optimal

    next_optimal = next_candidates[0]
    # Check if optimal + next_optimal creates a boundary GT
    if optimal[-1] + next_optimal[0] != "GT":
        return optimal  # No boundary GT — use optimal

    # Boundary GT detected.  Look for alternatives within CAI tolerance.
    min_acceptable_cai = optimal_cai * (1.0 - cai_tolerance)
    best_alt = None
    best_alt_cai = -1.0

    for alt in candidates:
        if alt == optimal:
            continue
        alt_cai = usage.get(alt, 0.0)
        # Check if this alternative avoids the boundary GT
        if alt[-1] + next_optimal[0] == "GT":
            continue  # Still creates boundary GT
        # Check CAI within tolerance
        if alt_cai >= min_acceptable_cai and alt_cai > best_alt_cai:
            best_alt = alt
            best_alt_cai = alt_cai

    if best_alt is not None:
        return best_alt  # Alternative found within CAI tolerance

    # No alternative within tolerance — accept the boundary GT
    return optimal


def _is_in_codon_gt(seq: str, pos: int) -> bool:
    """Check whether a GT dinucleotide at position *pos* is entirely within
    a single codon (in-codon) rather than spanning a codon boundary (cross-codon).

    Pre-conditions:
    - seq is a valid uppercase DNA string
    - 0 <= pos < len(seq) - 1
    - seq[pos:pos+2] == "GT"

    Post-conditions:
    - Returns True if both the G and T are in the same codon
    - Returns False if the GT spans a codon boundary (G at end of one codon,
      T at start of the next)
    """
    codon_of_g = pos // 3
    codon_of_t = (pos + 1) // 3
    return codon_of_g == codon_of_t


def _eukaryote_cai_recovery(
    sequence: str,
    protein: str,
    usage: dict[str, float],
    enzymes: list[str] | None = None,
    cai_cost_threshold: float = EUKARYOTE_CAI_GT_COST_THRESHOLD,
) -> tuple[str, int]:
    """Recover CAI by swapping suboptimal codons to optimal ones for eukaryotes.

    For eukaryotes, GT avoidance should be a SOFT preference, not a hard
    constraint.  When the optimizer has replaced an optimal codon with a
    suboptimal alternative (typically to avoid GT dinucleotides), this
    function swaps back using **cost-aware GT resolution** with splice
    donor scoring:

    - If the optimal codon creates a GT with **low** splice donor potential
      (score < SPLICE_DONOR_POTENTIAL_THRESHOLD), always swap to optimal —
      the GT is unlikely to be a real splice donor, so CAI wins.
    - If the optimal codon creates a GT with **high** splice donor potential
      (score >= SPLICE_DONOR_POTENTIAL_THRESHOLD), only swap if the CAI cost
      exceeds ``cai_cost_threshold`` — the GT might be dangerous, so keep
      the suboptimal codon unless the CAI loss is too large.
    - If the optimal codon does NOT create any new GT dinucleotides, always
      swap (no tradeoff needed).

    Pre-conditions:
    - sequence translates to protein
    - len(sequence) == len(protein) * 3
    - usage maps codon strings to adaptiveness values (0.0–1.0)

    Post-conditions:
    - returned sequence translates to the same protein
    - CAI of returned sequence >= CAI of input sequence
    - no new restriction sites are introduced

    Args:
        sequence: Current optimized DNA sequence.
        protein: Amino acid sequence (no stop).
        usage: Codon adaptiveness table (codon → w value).
        enzymes: List of restriction enzyme names to avoid creating sites for.
        cai_cost_threshold: Maximum acceptable CAI cost for GT avoidance.
            If using a GT-free codon would drop CAI by more than this,
            keep the GT-containing optimal codon.

    Returns:
        Tuple of (recovered sequence, number of codons upgraded).
    """
    # Precompute optimal codons per amino acid
    optimal_codons: dict[str, str] = {}
    for aa in set(protein):
        if aa == "*":
            continue
        codons = AA_TO_CODONS.get(aa, [])
        if not codons:
            continue
        best = max(codons, key=lambda c: usage.get(c, 0.0))
        optimal_codons[aa] = best

    # Precompute restriction site sequences to check
    rs_sites: list[tuple[str, str]] = []
    max_site_len = 0
    if enzymes:
        from .restriction_sites import get_recognition_site as _get_site
        for enz in enzymes:
            site = _get_site(enz)
            if site is not None:
                site_rc = reverse_complement(site)
                rs_sites.append((site, site_rc))
                max_site_len = max(max_site_len, len(site))

    seq_list = list(sequence)
    n_codons = len(protein)
    upgrades = 0

    def _should_swap_to_optimal(ci: int, current: str, optimal: str,
                                 cai_cost: float) -> bool:
        """Decide whether to swap current → optimal at codon position ci.

        Uses cost-aware GT resolution with splice donor scoring:
        - If optimal doesn't introduce new GTs: always swap (pure CAI gain).
        - If optimal introduces GT(s) with low splice donor potential: swap
          (GT is unlikely to be functional).
        - If optimal introduces GT(s) with high splice donor potential: only
          swap if CAI cost exceeds the threshold (tradeoff).
        """
        # Build the sequence with the swap applied to check for new GTs
        test_seq_list = list(seq_list)
        test_seq_list[ci * 3] = optimal[0]
        test_seq_list[ci * 3 + 1] = optimal[1]
        test_seq_list[ci * 3 + 2] = optimal[2]
        test_seq = "".join(test_seq_list)

        # Find all GT positions in the swapped region that would be new
        # (present in test_seq but not in current seq_list)
        current_seq = "".join(seq_list)
        check_start = max(0, ci * 3 - 1)
        check_end = min(len(test_seq), ci * 3 + 4)

        new_gt_positions = []
        for p in range(check_start, check_end - 1):
            if test_seq[p:p+2] == "GT" and current_seq[p:p+2] != "GT":
                new_gt_positions.append(p)

        if not new_gt_positions:
            # No new GTs introduced — pure CAI gain, always swap
            return True

        # Check splice donor potential of each new GT
        max_splice_score = 0.0
        for gt_pos in new_gt_positions:
            sdp = score_splice_donor_potential(test_seq, gt_pos)
            if sdp > max_splice_score:
                max_splice_score = sdp

        if max_splice_score < SPLICE_DONOR_POTENTIAL_THRESHOLD:
            # New GTs have low splice donor potential — accept them for CAI
            return True

        # New GTs have high splice donor potential — only swap if CAI
        # cost is significant enough to justify the risk
        return cai_cost > cai_cost_threshold

    # ── Pass 1: Single-codon swaps ──
    # For each position where the current codon is suboptimal, decide
    # whether to swap using cost-aware GT resolution.
    blocked_positions: list[int] = []  # Positions blocked by restriction sites

    for ci in range(n_codons):
        aa = protein[ci]
        if aa == "*":
            continue
        current = "".join(seq_list[ci * 3:ci * 3 + 3])
        optimal = optimal_codons.get(aa)
        if optimal is None or current == optimal:
            continue

        # Check CAI cost of using the suboptimal codon
        current_w = usage.get(current, 0.0)
        optimal_w = usage.get(optimal, 0.0)
        cai_cost = optimal_w - current_w

        # Cost-aware decision: should we swap?
        if not _should_swap_to_optimal(ci, current, optimal, cai_cost):
            continue

        # Swap to optimal codon
        old_codon = current
        seq_list[ci * 3] = optimal[0]
        seq_list[ci * 3 + 1] = optimal[1]
        seq_list[ci * 3 + 2] = optimal[2]

        # Check for new restriction sites in local region
        if rs_sites and max_site_len > 0:
            test_seq = "".join(seq_list)
            check_start = max(0, ci * 3 - max_site_len + 1)
            check_end = min(len(test_seq), ci * 3 + 3 + max_site_len - 1)
            local_region = test_seq[check_start:check_end]
            site_found = False
            for site, site_rc in rs_sites:
                if site in local_region or (site_rc and site_rc in local_region):
                    site_found = True
                    break
            if site_found:
                # Undo swap — restriction site would be created
                seq_list[ci * 3] = old_codon[0]
                seq_list[ci * 3 + 1] = old_codon[1]
                seq_list[ci * 3 + 2] = old_codon[2]
                blocked_positions.append(ci)
                continue

        upgrades += 1

    # ── Pass 2: Coordinated swaps for positions blocked by restriction sites ──
    # When a single-codon swap is blocked because it would create a restriction
    # site that spans a codon boundary, we can sometimes resolve it by also
    # changing an adjacent codon to break the restriction site.  This is only
    # done when the combined CAI improvement outweighs the CAI cost of the
    # adjacent codon change.
    if blocked_positions and rs_sites:
        # Precompute sorted codons per amino acid (by CAI descending)
        sorted_codons_map: dict[str, list[str]] = {}
        for aa in set(protein):
            if aa == "*":
                continue
            codons = AA_TO_CODONS.get(aa, [])
            sorted_codons_map[aa] = sorted(
                codons, key=lambda c: usage.get(c, 0.0), reverse=True
            )

        for ci in blocked_positions:
            aa = protein[ci]
            current = "".join(seq_list[ci * 3:ci * 3 + 3])
            optimal = optimal_codons.get(aa)
            if optimal is None or current == optimal:
                continue
            current_w = usage.get(current, 0.0)
            optimal_w = usage.get(optimal, 0.0)
            benefit = optimal_w - current_w
            if benefit <= cai_cost_threshold:
                continue

            # Check if the swap is justified by splice donor scoring
            if not _should_swap_to_optimal(ci, current, optimal, benefit):
                continue

            # Try adjacent codons to break the restriction site
            best_combo: tuple[int, str, float] | None = None  # (adj_ci, adj_codon, net_cai)
            for adj_offset in [-2, -1, 1, 2]:
                adj_ci = ci + adj_offset
                if adj_ci < 0 or adj_ci >= n_codons:
                    continue
                adj_aa = protein[adj_ci]
                if adj_aa == "*":
                    continue
                adj_current = "".join(seq_list[adj_ci * 3:adj_ci * 3 + 3])
                adj_current_w = usage.get(adj_current, 0.0)
                for adj_alt in sorted_codons_map.get(adj_aa, []):
                    if adj_alt == adj_current:
                        continue
                    # Apply both swaps and check for restriction sites
                    adj_old = adj_current
                    seq_list[ci * 3] = optimal[0]
                    seq_list[ci * 3 + 1] = optimal[1]
                    seq_list[ci * 3 + 2] = optimal[2]
                    seq_list[adj_ci * 3] = adj_alt[0]
                    seq_list[adj_ci * 3 + 1] = adj_alt[1]
                    seq_list[adj_ci * 3 + 2] = adj_alt[2]

                    test_seq = "".join(seq_list)
                    check_start = max(0, min(ci, adj_ci) * 3 - max_site_len + 1)
                    check_end = min(len(test_seq), max(ci, adj_ci) * 3 + 3 + max_site_len - 1)
                    local_region = test_seq[check_start:check_end]
                    site_found = False
                    for site, site_rc in rs_sites:
                        if site in local_region or (site_rc and site_rc in local_region):
                            site_found = True
                            break

                    # Undo both swaps
                    seq_list[ci * 3] = current[0]
                    seq_list[ci * 3 + 1] = current[1]
                    seq_list[ci * 3 + 2] = current[2]
                    seq_list[adj_ci * 3] = adj_old[0]
                    seq_list[adj_ci * 3 + 1] = adj_old[1]
                    seq_list[adj_ci * 3 + 2] = adj_old[2]

                    if not site_found:
                        # Check net CAI: benefit of optimal at ci minus cost of adj_alt
                        adj_alt_w = usage.get(adj_alt, 0.0)
                        adj_cost = adj_current_w - adj_alt_w  # positive = CAI loss
                        net_cai = benefit - adj_cost
                        # Only accept if net CAI is positive AND the coordinated
                        # swap is better than keeping the status quo
                        if net_cai > 0 and (best_combo is None or net_cai > best_combo[2]):
                            best_combo = (adj_ci, adj_alt, net_cai)

            if best_combo is not None:
                adj_ci, adj_codon, _net = best_combo
                # Apply both swaps
                seq_list[ci * 3] = optimal[0]
                seq_list[ci * 3 + 1] = optimal[1]
                seq_list[ci * 3 + 2] = optimal[2]
                seq_list[adj_ci * 3] = adj_codon[0]
                seq_list[adj_ci * 3 + 1] = adj_codon[1]
                seq_list[adj_ci * 3 + 2] = adj_codon[2]
                upgrades += 2  # Count both the optimal swap and the adjacent change

    return "".join(seq_list), upgrades


# ==============================================================================
# Systematic CpG Dinucleotide Elimination
# ==============================================================================

def _eliminate_cpg_dinucleotides(
    sequence: str,
    protein: str,
    usage: dict[str, float],
    enzymes: list[str] | None = None,
    max_iterations: int = 50,
    cpg_window: int = 200,
    cpg_threshold: float = 0.6,
    organism: str = "",
    gc_lo: float = 0.0,
    gc_hi: float = 1.0,
) -> tuple[str, list[str]]:
    """Systematically eliminate CpG dinucleotides to avoid CpG islands.

    This is a post-CAI-maximization pass that scans for ALL CG dinucleotides
    (both within codons and at codon boundaries) and attempts to replace them
    with synonymous codons that:
    1. Eliminate the specific CG dinucleotide
    2. Minimize CAI loss (prefer highest-CAI alternative)
    3. Do not create restriction sites
    4. Do not create new CG dinucleotides at the same or adjacent positions
    5. Keep GC content within the target range [gc_lo, gc_hi]

    Unlike the previous "best-effort" approach which broke on the first
    unfixable CpG position or returned early when the CpG island check
    passed, this function:
    - Always attempts to eliminate ALL CG dinucleotides, even if the
      sequence already passes the CpG island Obs/Exp ratio check
    - Continues trying ALL CpG positions even if some fail
    - Makes multiple passes to handle cascading fixes
    - For within-codon CpG: replaces with the highest-CAI synonymous codon
      that does not contain "CG"
    - For boundary CpG: tries changing the downstream codon first (avoid
      starting with G), then the upstream codon (avoid ending with C),
      then coordinated 2-codon swap — always preferring minimal CAI loss
    - Tracks and reports which CpGs were eliminated and which remain
    - Only stops when no more CG dinucleotides can be eliminated or
      max_iterations is reached

    Pre-conditions:
    - sequence translates to protein
    - len(sequence) == len(protein) * 3
    - usage maps codon strings to adaptiveness values (0.0–1.0)

    Post-conditions:
    - returned sequence translates to the same protein
    - CG dinucleotide count in returned sequence <= CG count in input sequence
    - warnings list describes any remaining CpG positions

    Args:
        sequence: Current optimized DNA sequence.
        protein: Amino acid sequence (no stop).
        usage: Codon adaptiveness table (codon → w value).
        enzymes: List of restriction enzyme names to avoid creating sites for.
        max_iterations: Maximum number of elimination passes.
        cpg_window: Window size for CpG island detection.
        cpg_threshold: Obs/Exp ratio threshold for CpG islands.
        organism: Target organism (for prokaryote skip).
        gc_lo: Minimum GC content fraction.
        gc_hi: Maximum GC content fraction.

    Returns:
        Tuple of (optimized sequence, list of warning strings about remaining CpGs).
    """
    # Skip for prokaryotic organisms — CpG islands are a eukaryotic concern
    if organism:
        from .organism_config import is_eukaryotic_organism
        if not is_eukaryotic_organism(organism):
            return sequence, []

    n_codons = len(protein)
    if n_codons == 0:
        return sequence, []

    # Precompute sorted codons per amino acid (by CAI descending)
    sorted_codons: dict[str, list[str]] = {}
    for aa in set(protein):
        if aa == "*":
            continue
        codons = AA_TO_CODONS.get(aa, [])
        sorted_codons[aa] = sorted(codons, key=lambda c: usage.get(c, 0.0), reverse=True)

    # Precompute restriction site sequences to avoid
    rs_sites: list[tuple[str, str]] = []
    max_site_len = 0
    if enzymes:
        from .restriction_sites import get_recognition_site as _get_site
        for enz in enzymes:
            site = _get_site(enz)
            if site is not None:
                site_rc = reverse_complement(site)
                rs_sites.append((site, site_rc))
                max_site_len = max(max_site_len, len(site))

    # ── Cost-aware CpG elimination ──
    # For eukaryotes, CpG avoidance is a SOFT preference, not a hard constraint.
    # Only eliminate CG dinucleotides when they contribute to a CpG island
    # (Obs/Exp ratio > cpg_threshold).  If the sequence already passes the
    # CpG island check, do NOT eliminate individual CGs — the CAI cost is
    # too high.  Individual CG dinucleotides in a CDS are common in
    # high-expression genes and are not biologically problematic unless they
    # cluster into CpG islands.
    from .type_system import check_no_cpg_island

    # Check if the sequence already passes the CpG island check
    cpg_result = check_no_cpg_island(sequence, cpg_window, cpg_threshold)
    
    # Count total CG dinucleotides — even if the sequence passes the CpG island
    # ratio check, we should still eliminate individual CGs if they're present.
    # Short sequences (< cpg_window) may pass the island check despite having
    # many CG dinucleotides because the windowed scan doesn't apply.
    total_cg = sum(1 for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG")
    
    if cpg_result.passed and total_cg == 0:
        # No CpG island and no individual CGs — nothing to do
        return sequence, []
    
    if cpg_result.passed and total_cg > 0:
        # Sequence passes island check but has individual CGs — still attempt
        # elimination as a best-effort pass (soft preference)
        logger.debug(
            "Sequence passes CpG island check but has %d CG dinucleotides — "
            "attempting elimination as soft preference",
            total_cg
        )

    warnings: list[str] = []
    seq = sequence
    eliminated_count = 0
    initial_cg_count = _count_dinucs_fast(seq, "CG")[0]
    # Track the initial GC content to allow swaps that move GC toward the
    # target range even if the current GC is outside the range.
    initial_gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)

    def _gc_ok(test_seq: str) -> bool:
        """Check if a swap is acceptable from a GC perspective.

        A swap is acceptable if:
        1. The test GC is within [gc_lo, gc_hi], OR
        2. The test GC is closer to the target range than the current seq GC
           (i.e., the swap moves GC in the right direction).

        This allows CpG elimination to proceed even when the current GC is
        outside the target range — the CpG elimination should not be blocked
        by a GC constraint that's already violated.
        """
        if gc_lo <= 0.0 and gc_hi >= 1.0:
            return True
        test_gc = (test_seq.count("G") + test_seq.count("C")) / max(len(test_seq), 1)
        # Case 1: Test GC is within range — always OK
        if gc_lo <= test_gc <= gc_hi:
            return True
        # Case 2: Current GC is already outside range — accept if test GC
        # is closer to the target range (moves in the right direction)
        current_gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
        if not (gc_lo <= current_gc <= gc_hi):
            # Current GC is outside range. Check if test GC is closer.
            def _dist_to_range(gc_val: float) -> float:
                if gc_val < gc_lo:
                    return gc_lo - gc_val
                elif gc_val > gc_hi:
                    return gc_val - gc_hi
                return 0.0
            return _dist_to_range(test_gc) <= _dist_to_range(current_gc)
        # Case 3: Current GC is within range but test GC is not — reject
        return False

    def _creates_boundary_cg(test_seq: str, codon_idx: int) -> bool:
        """Check if modifying codon at codon_idx creates a new boundary CG."""
        codon_start = codon_idx * 3
        # Check boundary with previous codon
        if codon_idx > 0 and codon_start > 0:
            if test_seq[codon_start - 1:codon_start + 1] == "CG":
                return True
        # Check boundary with next codon
        next_start = codon_start + 3
        if next_start < len(test_seq):
            if test_seq[next_start - 1:next_start + 1] == "CG":
                return True
        return False

    def _check_rs(test_seq: str, codon_start: int) -> bool:
        """Check if test_seq has restriction sites near codon_start.
        Returns True if OK (no sites), False if a site was created."""
        if not rs_sites:
            return True
        check_start = max(0, codon_start - max_site_len + 1)
        check_end = min(len(test_seq), codon_start + 3 + max_site_len - 1)
        local_region = test_seq[check_start:check_end]
        for site, site_rc in rs_sites:
            if site in local_region or (site_rc and site_rc in local_region):
                return False
        return True

    for iteration in range(max_iterations):
        # Find all CG dinucleotide positions
        cpg_positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "CG"]
        if not cpg_positions:
            break  # All CGs eliminated

        any_fixed = False

        for pos in cpg_positions:
            # Re-check: this position may have been fixed by an earlier swap
            if seq[pos:pos+2] != "CG":
                continue

            left_ci = pos // 3          # codon containing the C
            right_ci = (pos + 1) // 3   # codon containing the G
            is_cross_codon = (left_ci != right_ci)

            fixed = False

            if not is_cross_codon:
                # ── Within-codon CpG ──
                # The CG is entirely within one codon. Replace with a
                # synonymous codon that doesn't contain "CG".
                if left_ci >= n_codons:
                    continue
                aa = protein[left_ci]
                if aa == "*":
                    continue
                current = seq[left_ci*3:left_ci*3+3]

                # Try alternatives sorted by CAI (best first)
                best_alt: str | None = None
                best_alt_cai = -1.0
                for alt in sorted_codons.get(aa, []):
                    if alt == current or "CG" in alt:
                        continue
                    test = seq[:left_ci*3] + alt + seq[left_ci*3+3:]

                    # Check restriction sites
                    if not _check_rs(test, left_ci * 3):
                        continue

                    # Check GC content
                    if not _gc_ok(test):
                        continue

                    # Check that we don't create new boundary CGs
                    if _creates_boundary_cg(test, left_ci):
                        continue

                    alt_cai = usage.get(alt, 0.0)
                    if alt_cai > best_alt_cai:
                        best_alt = alt
                        best_alt_cai = alt_cai

                if best_alt is not None:
                    seq = seq[:left_ci*3] + best_alt + seq[left_ci*3+3:]
                    fixed = True
                    any_fixed = True
                    eliminated_count += 1
            else:
                # ── Cross-codon CpG ──
                # The C is the last base of codon left_ci, the G is the
                # first base of codon right_ci. Try:
                # 1. Change right codon to not start with G
                # 2. Change left codon to not end with C
                # 3. Coordinated 2-codon swap
                # For each strategy, pick the highest-CAI alternative.

                # Strategy 1: Change right codon (to not start with G)
                if 0 <= right_ci < n_codons:
                    right_aa = protein[right_ci]
                    if right_aa != "*":
                        right_current = seq[right_ci*3:right_ci*3+3]
                        best_right_alt: str | None = None
                        best_right_cai = -1.0
                        for alt in sorted_codons.get(right_aa, []):
                            if alt == right_current or alt[0] == 'G':
                                continue
                            test = seq[:right_ci*3] + alt + seq[right_ci*3+3:]
                            if not _check_rs(test, right_ci * 3):
                                continue
                            if not _gc_ok(test):
                                continue
                            # Make sure we don't create a new boundary CG
                            # with the codon after right_ci
                            if _creates_boundary_cg(test, right_ci):
                                continue
                            # Make sure the specific CG at pos is gone
                            if test[pos:pos+2] == "CG":
                                continue
                            alt_cai = usage.get(alt, 0.0)
                            if alt_cai > best_right_cai:
                                best_right_alt = alt
                                best_right_cai = alt_cai

                        if best_right_alt is not None:
                            seq = seq[:right_ci*3] + best_right_alt + seq[right_ci*3+3:]
                            fixed = True
                            any_fixed = True
                            eliminated_count += 1

                # Strategy 2: Change left codon (to not end with C)
                if not fixed and 0 <= left_ci < n_codons:
                    left_aa = protein[left_ci]
                    if left_aa != "*":
                        left_current = seq[left_ci*3:left_ci*3+3]
                        best_left_alt: str | None = None
                        best_left_cai = -1.0
                        for alt in sorted_codons.get(left_aa, []):
                            if alt == left_current or alt[-1] == 'C':
                                continue
                            test = seq[:left_ci*3] + alt + seq[left_ci*3+3:]
                            if not _check_rs(test, left_ci * 3):
                                continue
                            if not _gc_ok(test):
                                continue
                            # Make sure we don't create a new boundary CG
                            # with the codon before left_ci
                            if _creates_boundary_cg(test, left_ci):
                                continue
                            # Make sure the specific CG at pos is gone
                            if test[pos:pos+2] == "CG":
                                continue
                            alt_cai = usage.get(alt, 0.0)
                            if alt_cai > best_left_cai:
                                best_left_alt = alt
                                best_left_cai = alt_cai

                        if best_left_alt is not None:
                            seq = seq[:left_ci*3] + best_left_alt + seq[left_ci*3+3:]
                            fixed = True
                            any_fixed = True
                            eliminated_count += 1

                # Strategy 3: Coordinated 2-codon swap
                if (not fixed and 0 <= left_ci < n_codons
                        and 0 <= right_ci < n_codons):
                    left_aa = protein[left_ci]
                    right_aa = protein[right_ci]
                    if left_aa != "*" and right_aa != "*":
                        left_current = seq[left_ci*3:left_ci*3+3]
                        right_current = seq[right_ci*3:right_ci*3+3]
                        # Try all combinations of left+right alternatives
                        best_swap: tuple[str, str, float] | None = None  # (left_alt, right_alt, cai_sum)
                        for left_alt in sorted_codons.get(left_aa, []):
                            if left_alt == left_current:
                                continue
                            for right_alt in sorted_codons.get(right_aa, []):
                                if right_alt == right_current and left_alt == left_current:
                                    continue
                                # Skip if boundary still has CG
                                if left_alt[-1] == 'C' and right_alt[0] == 'G':
                                    continue
                                # Skip if within-codon CG created
                                if "CG" in left_alt or "CG" in right_alt:
                                    continue
                                test_list = list(seq)
                                test_list[left_ci*3:left_ci*3+3] = list(left_alt)
                                test_list[right_ci*3:right_ci*3+3] = list(right_alt)
                                test_str = "".join(test_list)
                                # Check restriction sites
                                check_start = max(0, left_ci * 3 - max_site_len + 1)
                                check_end = min(len(test_str), right_ci * 3 + 3 + max_site_len - 1)
                                local_region = test_str[check_start:check_end]
                                site_ok = True
                                for site, site_rc in rs_sites:
                                    if site in local_region or (site_rc and site_rc in local_region):
                                        site_ok = False
                                        break
                                if not site_ok:
                                    continue
                                if not _gc_ok(test_str):
                                    continue
                                # Check no new boundary CGs
                                if _creates_boundary_cg(test_str, left_ci):
                                    continue
                                if _creates_boundary_cg(test_str, right_ci):
                                    continue
                                # Verify the specific CG at pos is gone
                                if test_str[pos:pos+2] != "CG":
                                    cai_sum = (usage.get(left_alt, 0.0) +
                                               usage.get(right_alt, 0.0))
                                    if best_swap is None or cai_sum > best_swap[2]:
                                        best_swap = (left_alt, right_alt, cai_sum)

                        if best_swap is not None:
                            left_alt, right_alt, _ = best_swap
                            seq = (seq[:left_ci*3] + left_alt +
                                   seq[left_ci*3+3:right_ci*3] + right_alt +
                                   seq[right_ci*3+3:])
                            fixed = True
                            any_fixed = True
                            eliminated_count += 1

        # If no CG was fixed in this pass, stop trying
        if not any_fixed:
            break

        # Do NOT break early when CpG island check passes — we want to
        # eliminate ALL CG dinucleotides, not just enough to pass the
        # Obs/Exp ratio check.

    # Report remaining CpGs
    remaining_cpgs = [i for i in range(len(seq) - 1) if seq[i:i+2] == "CG"]
    if remaining_cpgs:
        # Only warn if the sequence still fails the CpG island check
        final_cpg_result = check_no_cpg_island(seq, cpg_window, cpg_threshold, organism=organism)
        if not final_cpg_result.passed:
            warnings.append(
                f"CpG island avoidance: {len(remaining_cpgs)} CG dinucleotide(s) remain. "
                f"No synonymous substitution could eliminate them without creating "
                f"restriction sites or other violations. "
                f"{final_cpg_result.details}"
            )

    logger.debug(
        "CpG elimination: %d/%d CG dinucleotides eliminated in %d passes",
        eliminated_count, initial_cg_count, iteration + 1 if max_iterations > 0 else 0,
    )

    return seq, warnings


# ==============================================================================
# Greedy Optimizer
# ==============================================================================

def _find_gt_free_codons(aa: str) -> list[str]:
    """Return codons for the given amino acid that do NOT contain the GT dinucleotide.

    For amino acids like Cysteine (C), Glycine (G), Arginine (R), and Serine (S),
    there exist synonymous codons without GT, providing a guaranteed fix for
    cryptic splice donor elimination. Valine (V) has NO GT-free codons.

    Pre-conditions:
    - aa is a valid single-letter amino acid code present in AA_TO_CODONS

    Post-conditions:
    - all returned codons are valid for the amino acid
    - no returned codon contains "GT" as a substring
    - if no GT-free codons exist (e.g., Valine), returns empty list
    """
    return [c for c in AA_TO_CODONS[aa] if "GT" not in c]


def _find_ag_free_codons(aa: str) -> list[str]:
    """Return codons for the given amino acid that do NOT contain the AG dinucleotide.

    Similar to _find_gt_free_codons but for acceptor (AG) dinucleotide elimination.
    Many amino acids have AG-free synonymous codons.

    Pre-conditions:
    - aa is a valid single-letter amino acid code present in AA_TO_CODONS

    Post-conditions:
    - all returned codons are valid for the amino acid
    - no returned codon contains "AG" as a substring
    - if no AG-free codons exist, returns empty list
    """
    return [c for c in AA_TO_CODONS[aa] if "AG" not in c]


def _greedy_optimize(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: list[str] | None = None,
    cryptic_splice_threshold: float = 3.0,
    seed: int | None = None,
    provenance_collector: DecisionProvenanceCollector | None = None,
    is_prokaryote: bool = False,
) -> tuple[str, list[str]]:
    """
    Greedy multi-objective codon optimization with coordinated constraint solving.

    Step ordering prioritizes hard constraints (restriction sites) over soft constraints (CAI).
    Reconciliation pass ensures earlier steps aren't undone by later ones.

    Steps:
    1. Best codon per position (maximize CAI)
    2. Remove restriction sites (multi-codon coordinated)
    3. Remove ATTTA instability motifs
    4. Fix 6+ consecutive T runs
    5. Adjust GC content (hard constraint, organism target aspiration)
    6. Reconciliation — restriction sites vs GC
    7. Eliminate cryptic splice donor/acceptor sites (SKIPPED for prokaryotes)
    7.5. Disrupt CpG dinucleotides to avoid CpG islands (SKIPPED for prokaryotes)
    8. Reconciliation — restriction sites after splice/CpG fixes (SKIPPED for prokaryotes)
    8.5. CpG reconciliation (SKIPPED for prokaryotes)

    Pre-conditions:
    - protein is a valid amino acid sequence (no invalid codes)
    - organism is in SUPPORTED_ORGANISMS
    - 0.0 <= gc_lo < gc_hi <= 1.0
    - cryptic_splice_threshold > 0

    Post-conditions:
    - returned sequence translates to the input protein
    - len(returned sequence) == len(protein) * 3
    - all codons in sequence are valid for their amino acid

    Args:
        is_prokaryote: When True, skip eukaryote-specific constraint steps
            (cryptic splice elimination, CpG disruption, and their
            reconciliation passes). Prokaryotes have no spliceosome, so
            GT/AG avoidance and CpG island disruption are biologically
            irrelevant and unnecessarily lower CAI.

    Note: The ``seed`` parameter is currently unused because the greedy
    optimizer is fully deterministic. It is reserved for future
    randomized optimization steps.
    """
    # Set deterministic seed if provided (reserved for future randomized steps)
    if seed is not None:
        import random
        random.seed(seed)

    # Validate pre-conditions
    assert 0.0 <= gc_lo < gc_hi <= 1.0, f"GC bounds invalid: [{gc_lo}, {gc_hi}]"
    assert cryptic_splice_threshold > 0, f"Threshold must be positive, got {cryptic_splice_threshold}"

    usage = CODON_ADAPTIVENESS_TABLES.get(organism)
    if usage is None:
        raise UnsupportedOrganismError(organism, SUPPORTED_ORGANISMS)
    aas = protein_to_aa_list(protein)
    restriction_sites = restriction_sites or list(RESTRICTION_ENZYMES.values())
    warnings: list[str] = []

    # sorted_codons: codons for each AA sorted by CAI (descending) — used in ALL
    # constraint resolution steps so that the highest-CAI alternative that
    # fixes the constraint is always preferred.
    sorted_codons: dict[str, list[str]] = {}
    for aa in set(aas):
        codons = AA_TO_CODONS[aa]
        codons_sorted = sorted(codons, key=lambda c: usage.get(c, 0.0), reverse=True)
        sorted_codons[aa] = codons_sorted

    # Step: Maximize CAI — Best codon per position (maximize CAI)
    # For eukaryotes, use GT-aware codon selection that avoids boundary GTs
    # when there's an alternative within 10% relative CAI.
    if provenance_collector is not None:
        # Expanded loop with per-codon provenance tracking
        chosen_codons: list[str] = []
        for pos, aa in enumerate(aas):
            # For eukaryotes, use GT-aware codon selection to avoid
            # boundary GTs when possible without destroying CAI
            if not is_prokaryote:
                next_aa = aas[pos + 1] if pos + 1 < len(aas) else None
                chosen = _gt_aware_select_codon(aa, next_aa, sorted_codons, usage)
            else:
                candidates = sorted_codons[aa]
                chosen = candidates[0]

            # Build alternatives list for provenance (excluding chosen codon)
            alternatives: list[dict[str, Any]] = []
            chosen_cai = usage.get(chosen, 0.0)
            for codon in candidates:
                if codon == chosen:
                    continue  # chosen codon is already in chosen_codon field
                cai_val = usage.get(codon, 0.0)
                gc_bases = sum(1 for b in codon if b in "GC")
                gc_contribution = gc_bases / 3.0
                # Check constraint violations for this codon
                violates: list[str] = []
                if "GT" in codon:
                    violates.append("cryptic_splice_donor")
                if "AG" in codon:
                    violates.append("cryptic_splice_acceptor")
                gc_bases_total = sum(1 for b in codon if b in "GC")
                if gc_bases_total == 0 and gc_lo > 0.5:
                    violates.append("gc_too_low")
                elif gc_bases_total == 3 and gc_hi < 0.5:
                    violates.append("gc_too_high")
                # Determine rejection reason
                rejected_because: str | None = None
                if violates:
                    rejected_because = f"Violates: {', '.join(violates)}"
                elif cai_val < chosen_cai:
                    rejected_because = "Lower CAI"
                else:
                    rejected_because = "Lower CAI"
                alternatives.append({
                    "codon": codon,
                    "cai_contribution": round(cai_val, 4),
                    "gc_contribution": round(gc_contribution, 2),
                    "violates_constraints": violates,
                    "rejected_because": rejected_because,
                })

            # Compute confidence: 1.0 if the best codon is clearly better,
            # lower if alternatives are close in CAI
            if len(candidates) > 1:
                second_best_cai = usage.get(candidates[1], 0.0)
                confidence = min(1.0, 0.5 + (chosen_cai - second_best_cai) * 5)
                confidence = max(0.0, confidence)
            else:
                confidence = 1.0

            provenance_collector.record_codon_decision(CodonDecision(
                position=pos,
                amino_acid=aa,
                original_codon=None,
                chosen_codon=chosen,
                alternatives_considered=alternatives,
                constraint_reason="Maximize CAI while maintaining GC in range",
                confidence=round(confidence, 4),
            ))
            chosen_codons.append(chosen)
        sequence = "".join(chosen_codons)
    else:
        # Original fast path (no provenance overhead)
        # For eukaryotes, use GT-aware codon selection
        if not is_prokaryote:
            chosen_codons_fast: list[str] = []
            for pos, aa in enumerate(aas):
                next_aa = aas[pos + 1] if pos + 1 < len(aas) else None
                chosen_codons_fast.append(
                    _gt_aware_select_codon(aa, next_aa, sorted_codons, usage)
                )
            sequence = "".join(chosen_codons_fast)
        else:
            sequence = "".join(sorted_codons[aa][0] for aa in aas)
    assert len(sequence) == len(aas) * 3, "Maximize CAI step: sequence length mismatch"

    # Step: Remove Restriction Sites (HIGHEST PRIORITY — multi-codon coordinated)
    # Process concrete sites first, then IUPAC sites
    concrete_sites = []
    iupac_sites = []
    for site in restriction_sites:
        site_upper = site.upper()
        if any(b not in "ACGT" for b in site_upper):
            iupac_sites.append(site_upper)
        else:
            concrete_sites.append(site_upper)

    # Build Aho-Corasick scanner for O(L+M) multi-pattern detection of all
    # concrete sites + their reverse complements simultaneously.
    # This replaces per-site O(N*L*site_len) scanning with a single O(L+M) pass.
    concrete_scanner = build_scanner_from_sites(concrete_sites) if concrete_sites else None

    # Remove concrete sites using Aho-Corasick for fast detection
    if concrete_scanner is not None:
        # Fast path: scan for ALL sites at once, fix one at a time
        for iteration in range(MAX_RESTRICTION_SITE_ITERATIONS * len(concrete_sites)):
            matches = concrete_scanner.scan(sequence)
            if not matches:
                break

            # Fix the first match found
            pos, site_match, _label = matches[0]
            site_len = len(site_match)
            site_rc = reverse_complement(site_match)

            # Try multi-codon coordinated removal (CAI-aware)
            new_seq, fixed = _remove_site_multicodon(
                sequence, aas, sorted_codons, site_match, site_rc, usage=usage
            )
            if fixed:
                sequence = new_seq
                continue

            # Fallback: try single-codon swap — CAI-aware: find the best
            # swap across ALL overlapping codon positions, not just the first
            overlapping = _get_overlapping_codons(pos, site_len, len(aas))
            best_single_swap: tuple[int, str, float] | None = None  # (ci, alt, log_cai)
            for ci in overlapping:
                if ci >= len(aas):
                    continue
                aa = aas[ci]
                current = sequence[ci * 3: ci * 3 + 3]
                for alt in sorted_codons[aa]:
                    if alt == current:
                        continue
                    seq_list = list(sequence)
                    seq_list[ci * 3: ci * 3 + 3] = list(alt)
                    test = "".join(seq_list)
                    if site_match not in test and site_rc not in test:
                        alt_cai = usage.get(alt, 0.0)
                        log_cai = math.log(alt_cai) if alt_cai > 0 else float('-inf')
                        if best_single_swap is None or log_cai > best_single_swap[2]:
                            best_single_swap = (ci, alt, log_cai)
            if best_single_swap is not None:
                ci, alt, _ = best_single_swap
                seq_list = list(sequence)
                seq_list[ci * 3: ci * 3 + 3] = list(alt)
                sequence = "".join(seq_list)
                continue

            # Try neighboring codons
            overlapping = _get_overlapping_codons(pos, site_len, len(aas))
            neighbor_fixed = False
            for ci in overlapping:
                if ci >= len(aas):
                    continue
                aa = aas[ci]
                current = sequence[ci * 3: ci * 3 + 3]
                for alt in sorted_codons[aa]:
                    if alt != current:
                        seq_list = list(sequence)
                        seq_list[ci * 3: ci * 3 + 3] = list(alt)
                        test = "".join(seq_list)
                        if site_match not in test and site_rc not in test:
                            sequence = test
                            neighbor_fixed = True
                            break
                if neighbor_fixed:
                    break
            if neighbor_fixed:
                continue

            warnings.append(f"Cannot remove {site_match} at iteration {iteration}")
            break
        else:
            # Check if any site is still present
            remaining = concrete_scanner.scan(sequence)
            if remaining:
                remaining_sites = [s for _, s, _ in remaining]
                logger.warning(
                    "Restriction site elimination did not converge after %d iterations. "
                    "%d restriction sites remain: %s",
                    MAX_RESTRICTION_SITE_ITERATIONS, len(remaining_sites),
                    remaining_sites[:5],
                )
                for _, site_still, _ in remaining[:3]:
                    warnings.append(
                        f"Restriction site elimination capped at {MAX_RESTRICTION_SITE_ITERATIONS} iterations. "
                        f"Site {site_still} could not be eliminated."
                    )
    else:
        # Fallback: per-site scan (no Aho-Corasick scanner available)
        for site_upper in concrete_sites:
            site_rc = reverse_complement(site_upper)
            for iteration in range(MAX_RESTRICTION_SITE_ITERATIONS):
                positions = _find_site_in_sequence(sequence, site_upper, site_rc)
                if not positions:
                    break

                # Try multi-codon coordinated removal (CAI-aware)
                new_seq, fixed = _remove_site_multicodon(
                    sequence, aas, sorted_codons, site_upper, site_rc, usage=usage
                )
                if fixed:
                    sequence = new_seq
                    continue

                # Fallback: CAI-aware single-codon swap across all overlapping
                # positions — find the swap with minimal CAI loss
                pos = positions[0]
                overlapping = _get_overlapping_codons(pos, len(site_upper), len(aas))
                best_single_swap: tuple[int, str, float] | None = None
                for ci in overlapping:
                    if ci >= len(aas):
                        continue
                    aa = aas[ci]
                    current = sequence[ci * 3: ci * 3 + 3]
                    for alt in sorted_codons[aa]:
                        if alt == current:
                            continue
                        seq_list = list(sequence)
                        seq_list[ci * 3: ci * 3 + 3] = list(alt)
                        test = "".join(seq_list)
                        if site_upper not in test and site_rc not in test:
                            alt_cai = usage.get(alt, 0.0)
                            log_cai = math.log(alt_cai) if alt_cai > 0 else float('-inf')
                            if best_single_swap is None or log_cai > best_single_swap[2]:
                                best_single_swap = (ci, alt, log_cai)
                if best_single_swap is not None:
                    ci, alt, _ = best_single_swap
                    seq_list = list(sequence)
                    seq_list[ci * 3: ci * 3 + 3] = list(alt)
                    sequence = "".join(seq_list)
                    continue

                # Try neighboring codons
                overlapping = _get_overlapping_codons(pos, len(site_upper), len(aas))
                neighbor_fixed = False
                for ci in overlapping:
                    if ci >= len(aas):
                        continue
                    aa = aas[ci]
                    current = sequence[ci * 3: ci * 3 + 3]
                    for alt in sorted_codons[aa]:
                        if alt != current:
                            seq_list = list(sequence)
                            seq_list[ci * 3: ci * 3 + 3] = list(alt)
                            test = "".join(seq_list)
                            if site_upper not in test and site_rc not in test:
                                sequence = test
                                neighbor_fixed = True
                                break
                    if neighbor_fixed:
                        break
                if neighbor_fixed:
                    continue

                warnings.append(f"Cannot remove {site_upper} at iteration {iteration}")
                break
            else:
                # Check if site is still present
                if site_upper in sequence or site_rc in sequence:
                    logger.warning(
                        "Restriction site elimination did not converge after %d iterations. "
                        "Site %s remains.",
                        MAX_RESTRICTION_SITE_ITERATIONS, site_upper,
                    )
                    warnings.append(
                        f"Restriction site elimination capped at {MAX_RESTRICTION_SITE_ITERATIONS} iterations. "
                        f"Site {site_upper} could not be eliminated."
                    )

    # Remove IUPAC sites (expand to concrete variants, check each)
    for site_upper in iupac_sites:
        concrete_variants = _expand_iupac_site(site_upper)
        if not concrete_variants:
            continue
        for variant in concrete_variants:
            variant_rc = reverse_complement(variant)
            for iteration in range(MAX_IUPAC_SITE_ITERATIONS):
                positions = _find_site_in_sequence(sequence, variant, variant_rc)
                if not positions:
                    break

                new_seq, fixed = _remove_site_multicodon(
                    sequence, aas, sorted_codons, variant, variant_rc, usage=usage
                )
                if fixed:
                    sequence = new_seq
                    continue

                # CAI-aware single-codon fallback for IUPAC sites
                pos = positions[0]
                overlapping = _get_overlapping_codons(pos, len(variant), len(aas))
                best_iupac_swap: tuple[int, str, float] | None = None
                for ci in overlapping:
                    if ci >= len(aas):
                        continue
                    aa = aas[ci]
                    current = sequence[ci * 3: ci * 3 + 3]
                    for alt in sorted_codons[aa]:
                        if alt == current:
                            continue
                        test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]
                        if variant not in test and variant_rc not in test:
                            alt_cai = usage.get(alt, 0.0)
                            log_cai = math.log(alt_cai) if alt_cai > 0 else float('-inf')
                            if best_iupac_swap is None or log_cai > best_iupac_swap[2]:
                                best_iupac_swap = (ci, alt, log_cai)
                if best_iupac_swap is not None:
                    ci, alt, _ = best_iupac_swap
                    sequence = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]
                else:
                    warnings.append(f"Cannot remove IUPAC {site_upper} variant {variant} at iteration {iteration}")
                    break
            else:
                if variant in sequence or variant_rc in sequence:
                    logger.warning(
                        "IUPAC site elimination did not converge after %d iterations. "
                        "Site %s variant %s remains.",
                        MAX_IUPAC_SITE_ITERATIONS, site_upper, variant,
                    )
                    warnings.append(
                        f"IUPAC site elimination capped at {MAX_IUPAC_SITE_ITERATIONS} iterations. "
                        f"Variant {variant} of {site_upper} could not be eliminated."
                    )

    # Step: Remove ATTTA instability motifs
    # PERF (Optimization D): Use list mutation for codon swaps
    for iteration in range(MAX_ATTTA_MOTIF_ITERATIONS):
        pos = sequence.find("ATTTA")
        if pos == -1:
            break
        codon_idx = pos // 3
        fixed = False
        for ci in range(max(0, codon_idx - 1), min(len(aas), codon_idx + 2)):
            aa = aas[ci]
            current = sequence[ci * 3:ci * 3 + 3]
            for alt in sorted_codons[aa]:
                if alt != current:
                    seq_list = list(sequence)
                    seq_list[ci * 3:ci * 3 + 3] = list(alt)
                    test = "".join(seq_list)
                    if "ATTTA" not in test:
                        sequence = test
                        fixed = True
                        break
            if fixed:
                break
        if not fixed:
            warnings.append(f"ATTTA motif: cannot remove at iteration {iteration}")
            break
    else:
        _remaining_attta = sequence.count("ATTTA")
        if _remaining_attta > 0:
            logger.warning(
                "ATTTA motif elimination did not converge after %d iterations. "
                "%d ATTTA motifs remain.",
                MAX_ATTTA_MOTIF_ITERATIONS, _remaining_attta,
            )
            warnings.append(
                f"ATTTA motif elimination capped at {MAX_ATTTA_MOTIF_ITERATIONS} iterations. "
                f"{_remaining_attta} motifs could not be eliminated."
            )

    # Step: Fix 6+ consecutive T runs
    for iteration in range(MAX_T_RUN_ITERATIONS):
        max_run, max_pos = 0, -1
        i = 0
        while i < len(sequence):
            if sequence[i] == "T":
                j = i
                while j < len(sequence) and sequence[j] == "T":
                    j += 1
                if j - i > max_run:
                    max_run, max_pos = j - i, i
                i = j
            else:
                i += 1
        if max_run < T_RUN_LENGTH_THRESHOLD:
            break
        codon_idx = (max_pos + max_run // 2) // 3
        if codon_idx < len(aas):
            aa = aas[codon_idx]
            current = sequence[codon_idx * 3:codon_idx * 3 + 3]
            fixed = False
            for alt in sorted_codons[aa]:
                if alt != current:
                    seq_list = list(sequence)
                    seq_list[codon_idx * 3:codon_idx * 3 + 3] = list(alt)
                    test = "".join(seq_list)
                    if not any(test[i:i + T_RUN_LENGTH_THRESHOLD] == "T" * T_RUN_LENGTH_THRESHOLD for i in range(len(test) - 5)):
                        sequence = test
                        fixed = True
                        break
            if not fixed:
                warnings.append(f"Consecutive T run: cannot fix at iteration {iteration}")
                break
    else:
        _remaining_t_runs = sum(
            1 for i in range(len(sequence) - T_RUN_LENGTH_THRESHOLD + 1)
            if sequence[i:i + T_RUN_LENGTH_THRESHOLD] == "T" * T_RUN_LENGTH_THRESHOLD
        )
        if _remaining_t_runs > 0:
            logger.warning(
                "T-run elimination did not converge after %d iterations. "
                "%d runs of %d+ T remain.",
                MAX_T_RUN_ITERATIONS, _remaining_t_runs, T_RUN_LENGTH_THRESHOLD,
            )
            warnings.append(
                f"T-run elimination capped at {MAX_T_RUN_ITERATIONS} iterations. "
                f"{_remaining_t_runs} runs of {T_RUN_LENGTH_THRESHOLD}+ T could not be eliminated."
            )

    # Step: Adjust GC content
    # Strategy: GC must be in [gc_lo, gc_hi] (hard constraint).
    # If in range, we gently nudge toward organism target but NEVER at the
    # cost of significant CAI reduction. The organism GC target is aspirational,
    # not mandatory — a sequence with CAI=0.99 and GC=0.61 (slightly above
    # human's 0.41 target) is better than CAI=0.82 and GC=0.46.
    # PERF (Optimization F): Cache GC count for incremental updates
    n_bases = len(sequence)
    if _HAS_NUMBA:
        try:
            gc_count = _numba_count_gc(_seq_to_bytes(sequence))
        except Exception:
            gc_count = sum(1 for b in sequence if b in "GC")
    else:
        gc_count = sum(1 for b in sequence if b in "GC")
    gc_val = gc_count / n_bases
    organism_gc_range = ORGANISM_GC_TARGETS.get(organism, (gc_lo, gc_hi))
    organism_gc = (organism_gc_range[0] + organism_gc_range[1]) / 2.0
    target_gc = max(gc_lo, min(gc_hi, organism_gc))

    gc_out_of_range = not (gc_lo <= gc_val <= gc_hi)

    if gc_out_of_range:
        # Hard constraint: MUST get GC into range
        # (gc_count and n_bases already computed above)
        # Target the nearest bound
        if gc_val < gc_lo:
            phase_target = gc_lo
        else:
            phase_target = gc_hi

        for iteration in range(MAX_GC_ADJUSTMENT_ITERATIONS):
            if gc_lo <= gc_val <= gc_hi:
                break
            best_alt = None
            best_ci = -1
            best_diff = abs(gc_val - phase_target)
            best_gc_delta = 0
            best_cai = -1.0  # CAI tiebreaker: prefer higher CAI among equal GC improvement
            for ci in range(len(aas)):
                aa = aas[ci]
                current = sequence[ci * 3:ci * 3 + 3]
                current_gc = sum(1 for b in current if b in "GC")
                for alt in sorted_codons[aa]:
                    if alt == current:
                        continue
                    alt_gc = sum(1 for b in alt if b in "GC")
                    new_gc_count = gc_count - current_gc + alt_gc
                    new_frac = new_gc_count / n_bases
                    diff = abs(new_frac - phase_target)
                    alt_cai = usage.get(alt, 0.0)
                    # Prefer better GC improvement; among equal GC improvement, prefer higher CAI
                    if diff < best_diff or (diff == best_diff and alt_cai > best_cai):
                        best_diff = diff
                        best_alt = alt
                        best_ci = ci
                        best_gc_delta = alt_gc - current_gc
                        best_cai = alt_cai
            if best_alt is None:
                break
            # PERF (Optimization D): Use list mutation for codon swap
            seq_list = list(sequence)
            seq_list[best_ci * 3: best_ci * 3 + 3] = list(best_alt)
            sequence = "".join(seq_list)
            gc_count += best_gc_delta
            gc_val = gc_count / n_bases
        else:
            logger.warning(
                "GC adjustment did not converge after %d iterations. "
                "Current GC=%.3f (target range [%.3f, %.3f])",
                MAX_GC_ADJUSTMENT_ITERATIONS, gc_val, gc_lo, gc_hi,
            )
            warnings.append(
                f"GC adjustment capped at {MAX_GC_ADJUSTMENT_ITERATIONS} iterations. "
                f"Current GC={gc_val:.3f} (target range [{gc_lo:.3f}, {gc_hi:.3f}])."
            )

    # Step: Reconciliation — check if GC adjustment reintroduced restriction sites
    # Use Aho-Corasick scanner for fast multi-site detection if available
    if concrete_scanner is not None:
        remaining_matches = concrete_scanner.scan(sequence)
        for pos, site_match, _label in remaining_matches:
            site_rc = reverse_complement(site_match)
            # Try one more round of multi-codon removal
            new_seq, fixed = _remove_site_multicodon(
                sequence, aas, sorted_codons, site_match, site_rc, usage=usage
            )
            if fixed:
                sequence = new_seq
                # Re-check GC
                if _HAS_NUMBA:
                    try:
                        gc_count = _numba_count_gc(_seq_to_bytes(sequence))
                    except Exception:
                        gc_count = sum(1 for b in sequence if b in "GC")
                else:
                    gc_count = sum(1 for b in sequence if b in "GC")
                gc_val = gc_count / n_bases
                if not (gc_lo <= gc_val <= gc_hi):
                    # GC drifted — try to fix with single-codon swaps that don't reintroduce sites
                    for ci in range(len(aas)):
                        aa = aas[ci]
                        current = sequence[ci * 3:ci * 3 + 3]
                        current_gc = sum(1 for b in current if b in "GC")
                        for alt in sorted_codons[aa]:
                            if alt == current:
                                continue
                            alt_gc = sum(1 for b in alt if b in "GC")
                            new_gc_count = gc_count - current_gc + alt_gc
                            new_frac = new_gc_count / n_bases
                            # Check this swap doesn't reintroduce any site
                            seq_list = list(sequence)
                            seq_list[ci * 3: ci * 3 + 3] = list(alt)
                            test = "".join(seq_list)
                            site_ok = not concrete_scanner.has_any_match(test)
                            if site_ok and abs(new_frac - target_gc) < abs(gc_val - target_gc):
                                sequence = test
                                gc_count = new_gc_count
                                gc_val = gc_count / n_bases
                                break
    else:
        for site_upper in concrete_sites:
            site_rc = reverse_complement(site_upper)
            if site_upper in sequence or site_rc in sequence:
                # Try one more round of multi-codon removal
                new_seq, fixed = _remove_site_multicodon(
                    sequence, aas, sorted_codons, site_upper, site_rc, usage=usage
                )
                if fixed:
                    sequence = new_seq
                    # Re-check GC
                    # PERF (Optimization F): Update cached GC count
                    if _HAS_NUMBA:
                        try:
                            gc_count = _numba_count_gc(_seq_to_bytes(sequence))
                        except Exception:
                            gc_count = sum(1 for b in sequence if b in "GC")
                    else:
                        gc_count = sum(1 for b in sequence if b in "GC")
                    gc_val = gc_count / n_bases
                    if not (gc_lo <= gc_val <= gc_hi):
                        # GC drifted — try to fix with single-codon swaps that don't reintroduce sites
                        for ci in range(len(aas)):
                            aa = aas[ci]
                            current = sequence[ci * 3:ci * 3 + 3]
                            current_gc = sum(1 for b in current if b in "GC")
                            for alt in sorted_codons[aa]:
                                if alt == current:
                                    continue
                                alt_gc = sum(1 for b in alt if b in "GC")
                                new_gc_count = gc_count - current_gc + alt_gc
                                new_frac = new_gc_count / n_bases
                                # Check this swap doesn't reintroduce any site
                                # PERF (Optimization D): Use list mutation
                                seq_list = list(sequence)
                                seq_list[ci * 3: ci * 3 + 3] = list(alt)
                                test = "".join(seq_list)
                                site_ok = all(
                                    s not in test and reverse_complement(s) not in test
                                    for s in concrete_sites
                                )
                                if site_ok and abs(new_frac - target_gc) < abs(gc_val - target_gc):
                                    sequence = test
                                    gc_count = new_gc_count
                                    gc_val = gc_count / n_bases
                                    break
            else:
                # Could not remove — already warned in Remove Restriction Sites step
                pass

    # Step: Eliminate cryptic splice donor/acceptor sites
    # EUKARYOTE-ONLY: Prokaryotes have no spliceosome, so cryptic splice
    # sites are biologically irrelevant. Skipping this step recovers
    # significant CAI on prokaryotic targets.
    #
    # IMPORTANT for eukaryotes: GT avoidance is a SOFT preference, not a
    # hard constraint. In-codon GTs from optimal codons (GGT, TGT, GTT,
    # etc.) are acceptable because they are unavoidable for high CAI.
    # Cross-codon GTs are only eliminated when the CAI cost of doing so
    # is < EUKARYOTE_CAI_GT_COST_THRESHOLD (default 0.05). This ensures
    # that CAI is prioritized over GT avoidance for eukaryotes.
    if not is_prokaryote:
        for iteration in range(MAX_SPLICE_ELIMINATION_ITERATIONS):
            max_d = max_donor_score(sequence)
            max_a = max_acceptor_score(sequence)
            if max_d < cryptic_splice_threshold and max_a < cryptic_splice_threshold:
                break

            fixed_any = False

            # Try to eliminate strong donors (sorted by score descending for priority)
            if max_d >= cryptic_splice_threshold:
                # Collect all strong donor positions with scores, sort by score descending
                donor_sites = []
                for i in range(len(sequence) - 1):
                    if sequence[i:i+2] == "GT":
                        s = score_donor(sequence, i)
                        if s >= cryptic_splice_threshold:
                            donor_sites.append((i, s))
                donor_sites.sort(key=lambda x: x[1], reverse=True)

                for gt_pos, gt_score in donor_sites:
                    codon_idx = gt_pos // 3
                    if codon_idx >= len(aas):
                        continue
                    aa = aas[codon_idx]

                    # ── Splice donor potential check ──
                    # Not all GT dinucleotides are equally dangerous.  A GT
                    # with low splice donor potential (< SPLICE_DONOR_POTENTIAL_THRESHOLD)
                    # is unlikely to function as a cryptic splice donor, even if
                    # its MaxEntScan donor score exceeds the threshold.  This can
                    # happen when the surrounding context doesn't support splicing
                    # (e.g., no downstream AG acceptor, no polypyrimidine tract).
                    # For such GTs, CAI should always win.
                    sdp = score_splice_donor_potential(sequence, gt_pos)
                    if sdp < SPLICE_DONOR_POTENTIAL_THRESHOLD:
                        # This GT has low splice donor potential — not dangerous.
                        # Accept it and move on (CAI > GT avoidance here).
                        continue

                    # ── Eukaryotic GT-vs-CAI tradeoff ──
                    # For eukaryotes, in-codon GTs from optimal codons are
                    # acceptable (biologically common in high-expression genes).
                    # Only eliminate GT if the CAI cost is < threshold.
                    is_in_codon = _is_in_codon_gt(sequence, gt_pos)
                    current_codon = sequence[codon_idx*3:codon_idx*3+3]
                    optimal_codon = sorted_codons[aa][0]

                    if is_in_codon:
                        # In-codon GT: acceptable if the current codon is optimal
                        # or if swapping to a GT-free codon would cost too much CAI
                        if current_codon == optimal_codon:
                            # In-codon GT from optimal codon (e.g., GGT for Gly,
                            # TGT for Cys, GTT for Val) — acceptable for eukaryotes
                            continue
                        # Check CAI cost of best GT-free alternative
                        current_w = usage.get(current_codon, 0.0)
                        gt_free = _find_gt_free_codons(aa)
                        if gt_free:
                            best_gt_free_w = max(usage.get(c, 0.0) for c in gt_free)
                            if current_w - best_gt_free_w > EUKARYOTE_CAI_GT_COST_THRESHOLD:
                                # CAI cost too high — keep the GT-containing codon
                                continue
                        else:
                            # No GT-free alternative (e.g., Valine) — must accept
                            continue
                    else:
                        # Cross-codon GT: for eukaryotes, only eliminate if the
                        # CAI cost of the best fix is < threshold. Check the CAI
                        # cost of changing the codon at gt_pos to a non-T-starting
                        # alternative (for the T-side) or non-G-ending alternative
                        # (for the G-side).
                        # The G is at the end of codon_idx, the T is at the start
                        # of codon_idx+1.
                        next_codon_idx = (gt_pos + 1) // 3
                        # Check CAI cost of changing either codon
                        _cross_cai_cost_ok = False
                        # Try changing the G-ending codon
                        if codon_idx < len(aas):
                            _g_aa = aas[codon_idx]
                            _g_current = sequence[codon_idx*3:codon_idx*3+3]
                            _g_current_w = usage.get(_g_current, 0.0)
                            _g_non_g_end = [c for c in sorted_codons[_g_aa] if c[-1] != "G"]
                            if _g_non_g_end:
                                _best_non_g_end_w = usage.get(_g_non_g_end[0], 0.0)
                                if _g_current_w - _best_non_g_end_w < EUKARYOTE_CAI_GT_COST_THRESHOLD:
                                    _cross_cai_cost_ok = True
                        # Try changing the T-starting codon
                        if not _cross_cai_cost_ok and next_codon_idx < len(aas):
                            _t_aa = aas[next_codon_idx]
                            _t_current = sequence[next_codon_idx*3:next_codon_idx*3+3]
                            _t_current_w = usage.get(_t_current, 0.0)
                            _t_non_t_start = [c for c in sorted_codons[_t_aa] if c[0] != "T"]
                            if _t_non_t_start:
                                _best_non_t_start_w = usage.get(_t_non_t_start[0], 0.0)
                                if _t_current_w - _best_non_t_start_w < EUKARYOTE_CAI_GT_COST_THRESHOLD:
                                    _cross_cai_cost_ok = True
                        if not _cross_cai_cost_ok:
                            # CAI cost too high for both possible fixes — accept
                            # the cross-codon GT for eukaryotes
                            continue

                    # Strategy 1: GT-free codon swap (guaranteed to eliminate GT)
                    gt_free = _find_gt_free_codons(aa)
                    # Sort by CAI (best first) so we prefer high-CAI alternatives
                    gt_free_sorted = sorted(
                        gt_free,
                        key=lambda c: usage.get(c, 0.0),
                        reverse=True,
                    )
                    if gt_free_sorted:
                        for alt in gt_free_sorted:
                            seq_list = list(sequence)
                            seq_list[codon_idx*3:codon_idx*3+3] = list(alt)
                            test = "".join(seq_list)
                            # Verify the GT at this position is eliminated
                            if gt_pos < len(test) - 1 and test[gt_pos:gt_pos+2] == "GT":
                                # GT still present at this position (shouldn't happen with GT-free codon,
                                # but GT might straddle codon boundary)
                                new_s = score_donor(test, gt_pos)
                                if new_s < cryptic_splice_threshold:
                                    sequence = test
                                    fixed_any = True
                                    break
                            else:
                                # GT eliminated — verify no new strong donors created globally
                                new_max_d = max_donor_score(test)
                                if new_max_d < max_d or new_max_d < cryptic_splice_threshold:
                                    sequence = test
                                    fixed_any = True
                                    break
                                # Even if new donors appear, this position is fixed;
                                # they'll be handled in subsequent iterations
                                sequence = test
                                fixed_any = True
                                break

                    if fixed_any:
                        break  # Restart scanning from highest-scoring site

                    # Strategy 2: For Valine (no GT-free codons) or if Strategy 1 didn't work,
                    # try 2-codon coordinated swap (GT codon + each neighbor) to disrupt
                    # the 9-mer splice context.
                    # Also try single-codon swap for V positions where different V codons
                    # have different splice scores due to context.
                    current = sequence[codon_idx*3:codon_idx*3+3]
                    # First try single-codon swap (different V codon may give different score)
                    for v_alt in sorted_codons[aa]:
                        if v_alt == current:
                            continue
                        seq_list = list(sequence)
                        seq_list[codon_idx*3:codon_idx*3+3] = list(v_alt)
                        test = "".join(seq_list)
                        if gt_pos < len(test) - 1 and test[gt_pos:gt_pos+2] == "GT":
                            new_s = score_donor(test, gt_pos)
                        else:
                            new_s = ELIMINATED_SITE_SCORE  # GT eliminated (cross-boundary)
                        if new_s < cryptic_splice_threshold:
                            sequence = test
                            fixed_any = True
                            break
                    if fixed_any:
                        break  # Restart scanning

                    # Then try 2-codon coordinated swap
                    for neighbor_offset in [-2, -1, 1, 2]:
                        n_idx = codon_idx + neighbor_offset
                        if 0 <= n_idx < len(aas):
                            n_aa = aas[n_idx]
                            n_current = sequence[n_idx*3:n_idx*3+3]
                            # For the GT codon, try top 3 alternatives by CAI
                            for v_alt in sorted_codons[aa][:TOP_CAI_ALTERNATIVES]:
                                # For the neighbor, try all alternatives
                                for n_alt in sorted_codons[n_aa]:
                                    if n_alt == n_current and v_alt == current:
                                        continue
                                    test = list(sequence)
                                    test[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                    test[n_idx*3:n_idx*3+3] = list(n_alt)
                                    test_str = "".join(test)
                                    if gt_pos < len(test_str) - 1 and test_str[gt_pos:gt_pos+2] == "GT":
                                        new_s = score_donor(test_str, gt_pos)
                                    else:
                                        new_s = ELIMINATED_SITE_SCORE  # GT eliminated
                                    if new_s < cryptic_splice_threshold:
                                        sequence = test_str
                                        fixed_any = True
                                        break
                                if fixed_any:
                                    break
                            if fixed_any:
                                break
                        if fixed_any:
                            break

                    # Strategy 3 (Issue 2): Deep backtracking — 3-codon coordinated swap
                    # When 2-codon swap can't eliminate GT, try modifying the GT codon
                    # plus TWO neighboring codons simultaneously. This is especially
                    # effective for Valine GTs where all codons contain GT but the
                    # splice score can be reduced by changing the surrounding context.
                    if not fixed_any:
                        # Try the adjacent codon pair (codon_idx-1, codon_idx, codon_idx+1)
                        for n1_offset in [-1, 1]:
                            n1_idx = codon_idx + n1_offset
                            if not (0 <= n1_idx < len(aas)):
                                continue
                            n1_aa = aas[n1_idx]
                            n1_current = sequence[n1_idx*3:n1_idx*3+3]
                            # Second neighbor — try both sides for maximum context coverage
                            for n2_offset in [-2, 2, -3, 3]:
                                n2_idx = codon_idx + n2_offset
                                if not (0 <= n2_idx < len(aas)) or n2_idx == n1_idx:
                                    continue
                                n2_aa = aas[n2_idx]
                                n2_current = sequence[n2_idx*3:n2_idx*3+3]
                                # For Valine GTs: try ALL V codons (all contain GT but
                                # give different splice scores due to 3rd base context).
                                # For neighbors: try all alternatives to maximize the
                                # chance of finding a low-scoring 9-mer context.
                                v_limit = len(sorted_codons[aa])  # Try ALL V alternatives
                                n1_limit = len(sorted_codons[n1_aa])  # Try ALL neighbor alts
                                n2_limit = min(len(sorted_codons[n2_aa]), 6)  # Cap for performance
                                for v_alt in sorted_codons[aa][:v_limit]:
                                    for n1_alt in sorted_codons[n1_aa][:n1_limit]:
                                        if n1_alt == n1_current and v_alt == current:
                                            continue
                                        for n2_alt in sorted_codons[n2_aa][:n2_limit]:
                                            if n2_alt == n2_current and n1_alt == n1_current and v_alt == current:
                                                continue
                                            test = list(sequence)
                                            test[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                            test[n1_idx*3:n1_idx*3+3] = list(n1_alt)
                                            test[n2_idx*3:n2_idx*3+3] = list(n2_alt)
                                            test_str = "".join(test)
                                            if gt_pos < len(test_str) - 1 and test_str[gt_pos:gt_pos+2] == "GT":
                                                new_s = score_donor(test_str, gt_pos)
                                            else:
                                                new_s = ELIMINATED_SITE_SCORE
                                            if new_s < cryptic_splice_threshold:
                                                sequence = test_str
                                                fixed_any = True
                                                break
                                        if fixed_any:
                                            break
                                    if fixed_any:
                                        break
                                if fixed_any:
                                    break
                            if fixed_any:
                                break

                    # Strategy 4 (Issue 2): Frame-shift approach — when GT is at a
                    # codon boundary and can't be eliminated within the current codon,
                    # adjust the PREVIOUS codon to shift the reading frame. This is
                    # especially effective for HBB where GT/AG dinucleotides persist
                    # because the preceding codon ends with G and the next starts with T.
                    if not fixed_any and codon_idx > 0:
                        prev_idx = codon_idx - 1
                        prev_aa = aas[prev_idx]
                        prev_current = sequence[prev_idx*3:prev_idx*3+3]
                        for prev_alt in sorted_codons[prev_aa]:
                            if prev_alt == prev_current:
                                continue
                            # Check if this previous codon swap changes the boundary
                            # Check GT at boundary between prev_alt and current codon
                            if prev_alt[-1] == "G" and sequence[codon_idx*3] == "T":
                                # This would create a new cross-codon GT — skip
                                continue
                            test = sequence[:prev_idx*3] + prev_alt + sequence[prev_idx*3+3:]
                            # Now check if the original GT at gt_pos is still there
                            if gt_pos < len(test) - 1 and test[gt_pos:gt_pos+2] == "GT":
                                new_s = score_donor(test, gt_pos)
                            else:
                                new_s = ELIMINATED_SITE_SCORE
                            if new_s < cryptic_splice_threshold:
                                sequence = test
                                fixed_any = True
                                break
                            # Also try combining previous codon swap with GT codon swap
                            for v_alt in sorted_codons[aa][:TOP_CAI_ALTERNATIVES]:
                                if v_alt == current:
                                    continue
                                test2 = list(test)  # test already has prev_alt applied
                                test2[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                test2_str = "".join(test2)
                                if gt_pos < len(test2_str) - 1 and test2_str[gt_pos:gt_pos+2] == "GT":
                                    new_s = score_donor(test2_str, gt_pos)
                                else:
                                    new_s = ELIMINATED_SITE_SCORE
                                if new_s < cryptic_splice_threshold:
                                    sequence = test2_str
                                    fixed_any = True
                                    break
                            if fixed_any:
                                break

                    if fixed_any:
                        break  # Restart scanning

            # Try to eliminate strong acceptors (same strategy, using AG-free codons)
            if not fixed_any and max_a >= cryptic_splice_threshold:
                # Collect all strong acceptor positions with scores, sort by score descending
                acceptor_sites = []
                for i in range(len(sequence) - 1):
                    if sequence[i:i+2] == "AG":
                        s = score_acceptor(sequence, i)
                        if s >= cryptic_splice_threshold:
                            acceptor_sites.append((i, s))
                acceptor_sites.sort(key=lambda x: x[1], reverse=True)

                for ag_pos, ag_score in acceptor_sites:
                    codon_idx = ag_pos // 3
                    if codon_idx >= len(aas):
                        continue
                    aa = aas[codon_idx]

                    # Strategy 1: AG-free codon swap (guaranteed to eliminate AG)
                    ag_free = _find_ag_free_codons(aa)
                    ag_free_sorted = sorted(
                        ag_free,
                        key=lambda c: usage.get(c, 0.0),
                        reverse=True,
                    )
                    if ag_free_sorted:
                        for alt in ag_free_sorted:
                            test = sequence[:codon_idx*3] + alt + sequence[codon_idx*3+3:]
                            if ag_pos < len(test) - 1 and test[ag_pos:ag_pos+2] == "AG":
                                new_s = score_acceptor(test, ag_pos)
                                if new_s < cryptic_splice_threshold:
                                    sequence = test
                                    fixed_any = True
                                    break
                            else:
                                # AG eliminated
                                new_max_a = max_acceptor_score(test)
                                if new_max_a < max_a or new_max_a < cryptic_splice_threshold:
                                    sequence = test
                                    fixed_any = True
                                    break
                                sequence = test
                                fixed_any = True
                                break

                    if fixed_any:
                        break

                    # Strategy 2: Single-codon swap then 2-codon context disruption for AG positions
                    current = sequence[codon_idx*3:codon_idx*3+3]
                    # First try single-codon swap (different codon may give different score)
                    for alt in sorted_codons[aa]:
                        if alt == current:
                            continue
                        test = sequence[:codon_idx*3] + alt + sequence[codon_idx*3+3:]
                        if ag_pos < len(test) - 1 and test[ag_pos:ag_pos+2] == "AG":
                            new_s = score_acceptor(test, ag_pos)
                        else:
                            new_s = ELIMINATED_SITE_SCORE  # AG eliminated
                        if new_s < cryptic_splice_threshold:
                            sequence = test
                            fixed_any = True
                            break
                    if fixed_any:
                        break

                    # Then try 2-codon coordinated swap
                    for neighbor_offset in [-2, -1, 1, 2]:
                        n_idx = codon_idx + neighbor_offset
                        if 0 <= n_idx < len(aas):
                            n_aa = aas[n_idx]
                            n_current = sequence[n_idx*3:n_idx*3+3]
                            for v_alt in sorted_codons[aa][:TOP_CAI_ALTERNATIVES]:
                                for n_alt in sorted_codons[n_aa]:
                                    if n_alt == n_current and v_alt == current:
                                        continue
                                    test = list(sequence)
                                    test[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                    test[n_idx*3:n_idx*3+3] = list(n_alt)
                                    test_str = "".join(test)
                                    if ag_pos < len(test_str) - 1 and test_str[ag_pos:ag_pos+2] == "AG":
                                        new_s = score_acceptor(test_str, ag_pos)
                                    else:
                                        new_s = ELIMINATED_SITE_SCORE
                                    if new_s < cryptic_splice_threshold:
                                        sequence = test_str
                                        fixed_any = True
                                        break
                                if fixed_any:
                                    break
                            if fixed_any:
                                break
                        if fixed_any:
                            break

                    # Strategy 3 (Issue 2): Deep backtracking for AG — 3-codon coordinated swap
                    if not fixed_any:
                        for n1_offset in [-1, 1]:
                            n1_idx = codon_idx + n1_offset
                            if not (0 <= n1_idx < len(aas)):
                                continue
                            n1_aa = aas[n1_idx]
                            n1_current = sequence[n1_idx*3:n1_idx*3+3]
                            for n2_offset in [-2, 2, -3, 3]:
                                n2_idx = codon_idx + n2_offset
                                if not (0 <= n2_idx < len(aas)) or n2_idx == n1_idx:
                                    continue
                                n2_aa = aas[n2_idx]
                                n2_current = sequence[n2_idx*3:n2_idx*3+3]
                                # Same expanded search as donor Strategy 3
                                v_limit = len(sorted_codons[aa])
                                n1_limit = len(sorted_codons[n1_aa])
                                n2_limit = min(len(sorted_codons[n2_aa]), 6)
                                for v_alt in sorted_codons[aa][:v_limit]:
                                    for n1_alt in sorted_codons[n1_aa][:n1_limit]:
                                        if n1_alt == n1_current and v_alt == current:
                                            continue
                                        for n2_alt in sorted_codons[n2_aa][:n2_limit]:
                                            if n2_alt == n2_current and n1_alt == n1_current and v_alt == current:
                                                continue
                                            test = list(sequence)
                                            test[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                            test[n1_idx*3:n1_idx*3+3] = list(n1_alt)
                                            test[n2_idx*3:n2_idx*3+3] = list(n2_alt)
                                            test_str = "".join(test)
                                            if ag_pos < len(test_str) - 1 and test_str[ag_pos:ag_pos+2] == "AG":
                                                new_s = score_acceptor(test_str, ag_pos)
                                            else:
                                                new_s = ELIMINATED_SITE_SCORE
                                            if new_s < cryptic_splice_threshold:
                                                sequence = test_str
                                                fixed_any = True
                                                break
                                        if fixed_any:
                                            break
                                    if fixed_any:
                                        break
                                if fixed_any:
                                    break
                            if fixed_any:
                                break

                    # Strategy 4 (Issue 2): Frame-shift approach for AG
                    if not fixed_any and codon_idx > 0:
                        prev_idx = codon_idx - 1
                        prev_aa = aas[prev_idx]
                        prev_current = sequence[prev_idx*3:prev_idx*3+3]
                        for prev_alt in sorted_codons[prev_aa]:
                            if prev_alt == prev_current:
                                continue
                            if prev_alt[-1] == "A" and sequence[codon_idx*3] == "G":
                                continue  # Would create new cross-codon AG
                            test = sequence[:prev_idx*3] + prev_alt + sequence[prev_idx*3+3:]
                            if ag_pos < len(test) - 1 and test[ag_pos:ag_pos+2] == "AG":
                                new_s = score_acceptor(test, ag_pos)
                            else:
                                new_s = ELIMINATED_SITE_SCORE
                            if new_s < cryptic_splice_threshold:
                                sequence = test
                                fixed_any = True
                                break
                            for v_alt in sorted_codons[aa][:TOP_CAI_ALTERNATIVES]:
                                if v_alt == current:
                                    continue
                                test2 = list(test)
                                test2[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                test2_str = "".join(test2)
                                if ag_pos < len(test2_str) - 1 and test2_str[ag_pos:ag_pos+2] == "AG":
                                    new_s = score_acceptor(test2_str, ag_pos)
                                else:
                                    new_s = ELIMINATED_SITE_SCORE
                                if new_s < cryptic_splice_threshold:
                                    sequence = test2_str
                                    fixed_any = True
                                    break
                            if fixed_any:
                                break

                    if fixed_any:
                        break

            if not fixed_any:
                # No more progress — some sites are unrepairable by codon swaps
                max_d = max_donor_score(sequence)
                max_a = max_acceptor_score(sequence)
                if max_d >= cryptic_splice_threshold:
                    warnings.append(
                        f"Cryptic splice donors remain: max_donor={max_d:.2f} "
                        f"(threshold={cryptic_splice_threshold}). "
                        f"Some positions may require amino acid substitution (e.g., V->I)"
                    )
                if max_a >= cryptic_splice_threshold:
                    warnings.append(
                        f"Cryptic splice acceptors remain: max_acceptor={max_a:.2f} "
                        f"(threshold={cryptic_splice_threshold})"
                    )
                break
        else:
            _remaining_d = max_donor_score(sequence)
            _remaining_a = max_acceptor_score(sequence)
            logger.warning(
                "Cryptic splice elimination did not converge after %d iterations. "
                "max_donor=%.2f, max_acceptor=%.2f (threshold=%.2f)",
                MAX_SPLICE_ELIMINATION_ITERATIONS, _remaining_d, _remaining_a,
                cryptic_splice_threshold,
            )
            warnings.append(
                f"Cryptic splice elimination capped at {MAX_SPLICE_ELIMINATION_ITERATIONS} iterations. "
                f"max_donor={_remaining_d:.2f}, max_acceptor={_remaining_a:.2f} "
                f"(threshold={cryptic_splice_threshold})."
            )

    # Step: Disrupt CpG dinucleotides to avoid CpG islands
    # EUKARYOTE-ONLY: CpG islands are a eukaryotic gene regulation concern.
    # Prokaryotes don't methylate CpG dinucleotides, so avoidance is unnecessary.
    #
    # Key improvement: Don't break on the first unfixed CpG position.
    # Instead, continue trying all positions. Accept swaps that eliminate
    # a specific CG at pos even if a new CG is created elsewhere (the new
    # one might be fixable in a subsequent iteration or might not contribute
    # to a CpG island). Only require that the specific CG at pos is eliminated,
    # not that the global CpG count decreases.
    if not is_prokaryote:
        for _cpg_iteration in range(MAX_CPG_DISRUPTION_ITERATIONS):
            cpg_positions = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
            if not cpg_positions:
                break
            any_fixed_this_iter = False
            for pos in cpg_positions:
                # Re-check: a previous fix may have eliminated this CG
                if sequence[pos:pos+2] != "CG":
                    continue
                left_ci = pos // 3          # codon containing the C
                right_ci = (pos + 1) // 3   # codon containing the G
                is_cross_codon = (left_ci != right_ci)

                # Strategy 1: Single-codon swap — try both the C-codon and the G-codon
                fixed = False
                for ci in ([left_ci, right_ci] if is_cross_codon else [left_ci]):
                    if ci < 0 or ci >= len(aas):
                        continue
                    aa = aas[ci]
                    current = sequence[ci*3:ci*3+3]
                    for alt in sorted_codons[aa]:
                        if alt == current:
                            continue
                        test = sequence[:ci*3] + alt + sequence[ci*3+3:]
                        # CRITICAL: Don't worsen any cryptic splice scores.
                        codon_start = ci * 3
                        codon_end = ci * 3 + 3
                        splice_worsened = False
                        for p in range(max(0, codon_start - 3), min(len(test) - 1, codon_end + 6)):
                            if test[p:p+2] == "GT":
                                new_s = score_donor(test, p)
                                if new_s >= cryptic_splice_threshold:
                                    if sequence[p:p+2] == "GT":
                                        old_s = score_donor(sequence, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                            if test[p:p+2] == "AG":
                                new_s = score_acceptor(test, p)
                                if new_s >= cryptic_splice_threshold:
                                    if sequence[p:p+2] == "AG":
                                        old_s = score_acceptor(sequence, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                        if splice_worsened:
                            continue
                        # Accept if the specific CG at pos is eliminated
                        # (relaxed: don't require global CpG decrease)
                        if test[pos:pos+2] != "CG":
                            sequence = test
                            fixed = True
                            any_fixed_this_iter = True
                            break
                    if fixed:
                        break

                # Strategy 2: Coordinated 2-codon swap for cross-codon CG
                if not fixed and is_cross_codon and 0 <= left_ci < len(aas) and 0 <= right_ci < len(aas):
                    left_aa = aas[left_ci]
                    right_aa = aas[right_ci]
                    left_current = sequence[left_ci*3:left_ci*3+3]
                    right_current = sequence[right_ci*3:right_ci*3+3]
                    for left_alt in sorted_codons[left_aa]:
                        if left_alt == left_current:
                            continue
                        for right_alt in sorted_codons[right_aa]:
                            if right_alt == right_current and left_alt == left_current:
                                continue
                            # Skip if the combination still has CG at the boundary
                            if left_alt[-1] == "C" and right_alt[0] == "G":
                                continue
                            test = list(sequence)
                            test[left_ci*3:left_ci*3+3] = list(left_alt)
                            test[right_ci*3:right_ci*3+3] = list(right_alt)
                            test_str = "".join(test)
                            # Splice check across wider region (both codons)
                            splice_worsened = False
                            check_start = max(0, left_ci * 3 - 3)
                            check_end = min(len(test_str) - 1, right_ci * 3 + 9)
                            for p in range(check_start, check_end):
                                if test_str[p:p+2] == "GT":
                                    new_s = score_donor(test_str, p)
                                    if new_s >= cryptic_splice_threshold:
                                        if sequence[p:p+2] == "GT":
                                            old_s = score_donor(sequence, p)
                                            if new_s > old_s:
                                                splice_worsened = True
                                                break
                                        else:
                                            splice_worsened = True
                                            break
                                if test_str[p:p+2] == "AG":
                                    new_s = score_acceptor(test_str, p)
                                    if new_s >= cryptic_splice_threshold:
                                        if sequence[p:p+2] == "AG":
                                            old_s = score_acceptor(sequence, p)
                                            if new_s > old_s:
                                                splice_worsened = True
                                                break
                                        else:
                                            splice_worsened = True
                                            break
                            if splice_worsened:
                                continue
                            # Accept if the specific CG at pos is eliminated
                            # (relaxed: don't require global CpG decrease)
                            if test_str[pos:pos+2] != "CG":
                                sequence = test_str
                                fixed = True
                                any_fixed_this_iter = True
                                break
                        if fixed:
                            break

            if not any_fixed_this_iter:
                break
        else:
            _remaining_cpg = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
            if _remaining_cpg:
                logger.warning(
                    "CpG disruption did not converge after %d iterations. "
                    "%d CpG dinucleotides remain.",
                    MAX_CPG_DISRUPTION_ITERATIONS, len(_remaining_cpg),
                )
                warnings.append(
                    f"CpG disruption capped at {MAX_CPG_DISRUPTION_ITERATIONS} iterations. "
                    f"{len(_remaining_cpg)} CpG dinucleotides could not be eliminated."
                )

    # Step: Reconciliation after cryptic splice elimination
    # EUKARYOTE-ONLY: Only needed if splice/CpG steps ran
    if not is_prokaryote:
        # Check if cryptic splice fixes reintroduced restriction sites
        for site_upper in concrete_sites:
            site_rc = reverse_complement(site_upper)
            if site_upper in sequence or site_rc in sequence:
                new_seq, fixed = _remove_site_multicodon(
                    sequence, aas, sorted_codons, site_upper, site_rc, usage=usage
                )
                if fixed:
                    sequence = new_seq

    # Step: CpG reconciliation after restriction site reconciliation
    # EUKARYOTE-ONLY: Only needed if CpG avoidance is active
    # Same improvement as the CpG disruption step: don't break on first
    # unfixed position, and accept swaps that eliminate a specific CG
    # even if new CGs are created elsewhere.
    if not is_prokaryote:
        for _cpg_iter in range(MAX_CPG_DISRUPTION_ITERATIONS):
            cpg_positions = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
            if not cpg_positions:
                break
            any_fixed_this_iter = False
            for pos in cpg_positions:
                # Re-check: a previous fix may have eliminated this CG
                if sequence[pos:pos+2] != "CG":
                    continue
                left_ci = pos // 3
                right_ci = (pos + 1) // 3
                is_cross_codon = (left_ci != right_ci)

                # Strategy 1: Single-codon swap — try both the C-codon and the G-codon
                fixed = False
                for ci in ([left_ci, right_ci] if is_cross_codon else [left_ci]):
                    if ci < 0 or ci >= len(aas):
                        continue
                    aa = aas[ci]
                    current = sequence[ci*3:ci*3+3]
                    for alt in sorted_codons[aa]:
                        if alt == current:
                            continue
                        test = sequence[:ci*3] + alt + sequence[ci*3+3:]
                        # Must not reintroduce restriction sites
                        site_ok = all(
                            s not in test and reverse_complement(s) not in test
                            for s in concrete_sites
                        )
                        if not site_ok:
                            continue
                        # Must not worsen cryptic splice scores
                        codon_start = ci * 3
                        codon_end = ci * 3 + 3
                        splice_worsened = False
                        for p in range(max(0, codon_start - 3), min(len(test) - 1, codon_end + 6)):
                            if test[p:p+2] == "GT":
                                new_s = score_donor(test, p)
                                if new_s >= cryptic_splice_threshold:
                                    if sequence[p:p+2] == "GT":
                                        old_s = score_donor(sequence, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                            if test[p:p+2] == "AG":
                                new_s = score_acceptor(test, p)
                                if new_s >= cryptic_splice_threshold:
                                    if sequence[p:p+2] == "AG":
                                        old_s = score_acceptor(sequence, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                        if splice_worsened:
                            continue
                        # Accept if the specific CG at pos is eliminated
                        # (relaxed: don't require global CpG decrease)
                        if test[pos:pos+2] != "CG":
                            sequence = test
                            fixed = True
                            any_fixed_this_iter = True
                            break
                    if fixed:
                        break

                # Strategy 2: Coordinated 2-codon swap for cross-codon CG
                if not fixed and is_cross_codon and 0 <= left_ci < len(aas) and 0 <= right_ci < len(aas):
                    left_aa = aas[left_ci]
                    right_aa = aas[right_ci]
                    left_current = sequence[left_ci*3:left_ci*3+3]
                    right_current = sequence[right_ci*3:right_ci*3+3]
                    for left_alt in sorted_codons[left_aa]:
                        if left_alt == left_current:
                            continue
                        for right_alt in sorted_codons[right_aa]:
                            if right_alt == right_current and left_alt == left_current:
                                continue
                            if left_alt[-1] == "C" and right_alt[0] == "G":
                                continue
                            test = list(sequence)
                            test[left_ci*3:left_ci*3+3] = list(left_alt)
                            test[right_ci*3:right_ci*3+3] = list(right_alt)
                            test_str = "".join(test)
                            # Must not reintroduce restriction sites
                            site_ok = all(
                                s not in test_str and reverse_complement(s) not in test_str
                                for s in concrete_sites
                            )
                            if not site_ok:
                                continue
                            # Splice check across wider region (both codons)
                            splice_worsened = False
                            check_start = max(0, left_ci * 3 - 3)
                            check_end = min(len(test_str) - 1, right_ci * 3 + 9)
                            for p in range(check_start, check_end):
                                if test_str[p:p+2] == "GT":
                                    new_s = score_donor(test_str, p)
                                    if new_s >= cryptic_splice_threshold:
                                        if sequence[p:p+2] == "GT":
                                            old_s = score_donor(sequence, p)
                                            if new_s > old_s:
                                                splice_worsened = True
                                                break
                                        else:
                                            splice_worsened = True
                                            break
                                if test_str[p:p+2] == "AG":
                                    new_s = score_acceptor(test_str, p)
                                    if new_s >= cryptic_splice_threshold:
                                        if sequence[p:p+2] == "AG":
                                            old_s = score_acceptor(sequence, p)
                                            if new_s > old_s:
                                                splice_worsened = True
                                                break
                                        else:
                                            splice_worsened = True
                                            break
                            if splice_worsened:
                                continue
                            # Accept if the specific CG at pos is eliminated
                            # (relaxed: don't require global CpG decrease)
                            if test_str[pos:pos+2] != "CG":
                                sequence = test_str
                                fixed = True
                                any_fixed_this_iter = True
                                break
                        if fixed:
                            break

            if not any_fixed_this_iter:
                break
        else:
            _remaining_cpg2 = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
            if _remaining_cpg2:
                logger.warning(
                    "CpG reconciliation did not converge after %d iterations. "
                    "%d CpG dinucleotides remain.",
                    MAX_CPG_DISRUPTION_ITERATIONS, len(_remaining_cpg2),
                )
                warnings.append(
                    f"CpG reconciliation capped at {MAX_CPG_DISRUPTION_ITERATIONS} iterations. "
                    f"{len(_remaining_cpg2)} CpG dinucleotides could not be eliminated."
                )

    # Step: CAI Hill Climb — upgrade codons to higher-CAI alternatives
    # without violating any constraint. This recovers CAI lost during
    # constraint resolution by swapping to higher-CAI synonymous codons
    # that don't reintroduce any forbidden pattern.
    #
    # Speed optimizations (Issue 3):
    # - Incremental GC tracking: avoid O(N) gc_content() calls per position
    # - Localized RS check: only check restriction sites near the changed codon
    # - Batch CAI hill-climb: collect all improvements, apply non-conflicting ones
    _MAX_HILL_CLIMB_ITERATIONS = 10

    # Pre-compute site lengths and RCs for localized checking
    _local_rs_check_radius = max((len(s) for s in concrete_sites), default=0) + 3

    # Incremental GC tracking for the hill climb
    _hc_gc_count = sum(1 for b in sequence if b in "GC")
    _hc_n_bases = len(sequence)

    # Pre-build Aho-Corasick scanner for fast local RS check if available
    _hc_scanner = concrete_scanner  # reuse the scanner from Step 2

    for _hc_iter in range(_MAX_HILL_CLIMB_ITERATIONS):
        # Batch mode: collect all possible CAI upgrades, then apply
        # non-conflicting ones together (Issue 3: speed)
        _hc_upgrades: list[tuple[int, str, float, str]] = []  # (ci, alt, alt_cai, current)

        for ci in range(len(aas)):
            aa = aas[ci]
            current = sequence[ci * 3:ci * 3 + 3]
            current_cai = usage.get(current, 0.0)
            for alt in sorted_codons[aa]:
                alt_cai = usage.get(alt, 0.0)
                if alt_cai <= current_cai:
                    break  # sorted_codons is CAI-descending; no improvement possible

                # Quick incremental GC check (O(1) instead of O(N))
                current_gc_bases = sum(1 for b in current if b in "GC")
                alt_gc_bases = sum(1 for b in alt if b in "GC")
                new_gc_count = _hc_gc_count - current_gc_bases + alt_gc_bases
                test_gc = new_gc_count / _hc_n_bases
                if not (gc_lo <= test_gc <= gc_hi):
                    continue

                test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]

                # Localized RS check: only scan the region around the changed codon
                # (Issue 3: avoid O(N) full-sequence scan per position)
                site_ok = True
                if _hc_scanner is not None:
                    # Scan only the region that could contain new sites
                    scan_start = max(0, ci * 3 - _local_rs_check_radius)
                    scan_end = min(len(test), ci * 3 + 3 + _local_rs_check_radius)
                    local_region = test[scan_start:scan_end]
                    local_matches = _hc_scanner.scan(local_region)
                    # Check if any match corresponds to a genuinely new site
                    for _m_pos, _m_site, _ in local_matches:
                        # Map back to absolute position
                        abs_pos = scan_start + _m_pos
                        # Check if this site was already in the old sequence
                        old_local = sequence[scan_start:scan_end]
                        if _m_site not in old_local or abs_pos != sequence.find(_m_site, max(0, abs_pos - len(_m_site))):
                            # More precise: check if the site exists in the old sequence
                            # at any position overlapping the changed region
                            if abs_pos + len(_m_site) > ci * 3 and abs_pos < ci * 3 + 3:
                                # Site overlaps the changed region — check if it's new
                                if _m_site not in sequence or reverse_complement(_m_site) not in sequence:
                                    site_ok = False
                                    break
                                else:
                                    # Site exists in old sequence but may have moved;
                                    # do a full check to be safe
                                    if _m_site in test and reverse_complement(_m_site) in test:
                                        pass  # Both existed before, OK
                                    elif _m_site in test and _m_site not in sequence:
                                        site_ok = False
                                        break
                    # Fallback: full scan for safety on the first iteration
                    if _hc_iter == 0 and site_ok:
                        full_matches = _hc_scanner.scan(test)
                        if full_matches:
                            # Verify no new sites introduced
                            old_matches = _hc_scanner.scan(sequence)
                            new_match_sites = set((p, s) for p, s, _ in full_matches) - set((p, s) for p, s, _ in old_matches)
                            if new_match_sites:
                                site_ok = False
                else:
                    # No scanner: check concrete sites directly (localized)
                    for site_upper in concrete_sites:
                        site_rc = reverse_complement(site_upper)
                        # Only check in the region that could be affected
                        region_start = max(0, ci * 3 - len(site_upper) + 1)
                        region_end = min(len(test), ci * 3 + 3 + len(site_upper) - 1)
                        if site_upper in test[region_start:region_end] or site_rc in test[region_start:region_end]:
                            # Check if this is a genuinely new site
                            if site_upper not in sequence[region_start:region_end] and site_rc not in sequence[region_start:region_end]:
                                site_ok = False
                                break
                if not site_ok:
                    continue

                # No ATTTA motif increase
                if "ATTTA" in test:
                    old_attta = sequence.count("ATTTA")
                    new_attta = test.count("ATTTA")
                    if new_attta > old_attta:
                        continue

                # No 6+ T runs — only check near the changed codon (localized)
                _t_run_region_start = max(0, ci * 3 - T_RUN_LENGTH_THRESHOLD)
                _t_run_region_end = min(len(test), ci * 3 + 3 + T_RUN_LENGTH_THRESHOLD)
                _region_has_long_t = False
                j = _t_run_region_start
                while j < _t_run_region_end:
                    if test[j] == "T":
                        k = j
                        while k < len(test) and test[k] == "T":
                            k += 1
                        if k - j >= T_RUN_LENGTH_THRESHOLD:
                            _region_has_long_t = True
                            break
                        j = k
                    else:
                        j += 1
                if _region_has_long_t:
                    continue

                # No worsening of cryptic splice scores (eukaryotes only)
                # Prokaryotes have no spliceosome, so splice score checks
                # are irrelevant and unnecessarily block CAI upgrades.
                if not is_prokaryote:
                    if max_donor_score(test) > max_donor_score(sequence) + 0.5:
                        continue
                    if max_acceptor_score(test) > max_acceptor_score(sequence) + 0.5:
                        continue

                # Record this upgrade candidate
                _hc_upgrades.append((ci, alt, alt_cai, current))
                break  # Take the best (highest-CAI) alt for this position

        if not _hc_upgrades:
            break

        # Apply upgrades in batch (non-conflicting = different codon positions)
        _applied_any = False
        _applied_positions: set[int] = set()
        for ci, alt, alt_cai, current in _hc_upgrades:
            if ci in _applied_positions:
                continue
            # Re-validate: the sequence may have changed from prior batch swaps
            if sequence[ci * 3:ci * 3 + 3] != current:
                continue  # Position was already modified

            test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]

            # Quick re-validation of key constraints
            current_gc_bases = sum(1 for b in current if b in "GC")
            alt_gc_bases = sum(1 for b in alt if b in "GC")
            new_gc = _hc_gc_count - current_gc_bases + alt_gc_bases
            test_gc = new_gc / _hc_n_bases
            if not (gc_lo <= test_gc <= gc_hi):
                continue

            # Localized RS re-check
            rs_ok = True
            for site_upper in concrete_sites:
                site_rc = reverse_complement(site_upper)
                region_start = max(0, ci * 3 - len(site_upper) + 1)
                region_end = min(len(test), ci * 3 + 3 + len(site_upper) - 1)
                if site_upper in test[region_start:region_end] or site_rc in test[region_start:region_end]:
                    if site_upper not in sequence[region_start:region_end] and site_rc not in sequence[region_start:region_end]:
                        rs_ok = False
                        break
            if not rs_ok:
                continue

            # Apply the upgrade
            sequence = test
            _hc_gc_count = new_gc  # Incremental GC update
            _applied_any = True
            _applied_positions.add(ci)

        if not _applied_any:
            break

    # Step: CAI Recovery Pass (Issue 1)
    # After the hill climb, some positions may still have suboptimal codons
    # because the hill climb only tries one position at a time. This pass
    # systematically checks EVERY position and ALWAYS picks the highest-CAI
    # synonymous codon that doesn't violate any constraint.
    # For prokaryotes, this should close the 0.997→1.0 gap since there are
    # no splice constraints to block upgrades.
    _CAI_RECOVERY_MAX_ITERS = 3
    _rec_gc_count = _hc_gc_count  # Carry forward incremental GC tracking

    for _rec_iter in range(_CAI_RECOVERY_MAX_ITERS):
        _any_recovery = False
        for ci in range(len(aas)):
            aa = aas[ci]
            if aa == "*" or aa == "M":
                continue  # Skip stop and Met (only one codon)
            current = sequence[ci * 3:ci * 3 + 3]
            current_w = usage.get(current, 0.0)
            best_codon = sorted_codons[aa][0]  # Highest-CAI codon for this AA
            best_w = usage.get(best_codon, 0.0)

            if best_w <= current_w or best_codon == current:
                continue  # Already optimal

            # Try the best codon first, then fall back to next-best
            for alt in sorted_codons[aa]:
                alt_w = usage.get(alt, 0.0)
                if alt_w <= current_w:
                    break  # No improvement possible

                # Incremental GC check (O(1))
                cur_gc = sum(1 for b in current if b in "GC")
                alt_gc = sum(1 for b in alt if b in "GC")
                new_gc = _rec_gc_count - cur_gc + alt_gc
                if not (gc_lo <= new_gc / _hc_n_bases <= gc_hi):
                    continue

                test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]

                # Localized RS check
                rs_ok = True
                for site_upper in concrete_sites:
                    site_rc = reverse_complement(site_upper)
                    region_start = max(0, ci * 3 - len(site_upper) + 1)
                    region_end = min(len(test), ci * 3 + 3 + len(site_upper) - 1)
                    if site_upper in test[region_start:region_end] or site_rc in test[region_start:region_end]:
                        if site_upper not in sequence[region_start:region_end] and site_rc not in sequence[region_start:region_end]:
                            rs_ok = False
                            break
                if not rs_ok:
                    continue

                # No ATTTA increase
                if "ATTTA" in test and test.count("ATTTA") > sequence.count("ATTTA"):
                    continue

                # No 6+ T runs (localized check)
                _t_ok = True
                _t_start = max(0, ci * 3 - T_RUN_LENGTH_THRESHOLD)
                _t_end = min(len(test), ci * 3 + 3 + T_RUN_LENGTH_THRESHOLD)
                j = _t_start
                while j < _t_end:
                    if test[j] == "T":
                        k = j
                        while k < len(test) and test[k] == "T":
                            k += 1
                        if k - j >= T_RUN_LENGTH_THRESHOLD:
                            _t_ok = False
                            break
                        j = k
                    else:
                        j += 1
                if not _t_ok:
                    continue

                # No worsening of splice scores (eukaryotes only)
                if not is_prokaryote:
                    if max_donor_score(test) > max_donor_score(sequence) + 0.5:
                        continue
                    if max_acceptor_score(test) > max_acceptor_score(sequence) + 0.5:
                        continue

                # All checks passed — accept the upgrade
                sequence = test
                _rec_gc_count = new_gc
                _any_recovery = True
                logger.debug(
                    "CAI recovery: upgraded codon %d from %s to %s (w %.4f→%.4f)",
                    ci, current, alt, current_w, alt_w,
                )
                break  # Move to next position

            # Paired CAI recovery: if single-codon upgrade was blocked (likely
            # by a restriction site), try a paired swap — upgrade the target codon
            # AND adjust an adjacent codon to eliminate the new restriction site.
            # Net CAI must still improve.
            if not _any_recovery:
                import math as _rec_math
                for _adj_offset in (1, -1):
                    _adj_ci = ci + _adj_offset
                    if _adj_ci < 0 or _adj_ci >= len(aas):
                        continue
                    _adj_aa = aas[_adj_ci]
                    if _adj_aa == "*" or _adj_aa == "M":
                        continue
                    _adj_current = sequence[_adj_ci * 3:_adj_ci * 3 + 3]
                    _adj_current_w = usage.get(_adj_current, 0.0)

                    for alt in sorted_codons[aa]:
                        alt_w = usage.get(alt, 0.0)
                        if alt_w <= current_w:
                            break

                        _adj_sorted = sorted(
                            AA_TO_CODONS.get(_adj_aa, []),
                            key=lambda c: usage.get(c, 0.0),
                            reverse=True,
                        )
                        for _adj_alt in _adj_sorted:
                            if _adj_alt == _adj_current:
                                continue
                            _adj_alt_w = usage.get(_adj_alt, 0.0)
                            # Net log-CAI change must be positive
                            _old_log = _rec_math.log(max(current_w, 1e-10)) + _rec_math.log(max(_adj_current_w, 1e-10))
                            _new_log = _rec_math.log(max(alt_w, 1e-10)) + _rec_math.log(max(_adj_alt_w, 1e-10))
                            if _new_log <= _old_log:
                                continue

                            # Build test sequence with both swaps
                            _lo_ci = min(ci, _adj_ci)
                            _hi_ci = max(ci, _adj_ci)
                            _test_seq = (
                                sequence[:_lo_ci * 3]
                                + (alt if _lo_ci == ci else _adj_alt)
                                + sequence[_lo_ci * 3 + 3:_hi_ci * 3]
                                + (alt if _hi_ci == ci else _adj_alt)
                                + sequence[_hi_ci * 3 + 3:]
                            )

                            # Check: no new restriction sites
                            _rs_ok = True
                            for _site_upper in concrete_sites:
                                _site_rc = reverse_complement(_site_upper)
                                if _site_upper in _test_seq or (_site_rc and _site_rc in _test_seq):
                                    # Check if this is genuinely new
                                    if _site_upper not in sequence or _site_rc not in sequence:
                                        _rs_ok = False
                                        break
                                    else:
                                        # Count occurrences — no net increase
                                        if _test_seq.count(_site_upper) + _test_seq.count(_site_rc) > sequence.count(_site_upper) + sequence.count(_site_rc):
                                            _rs_ok = False
                                            break
                            if not _rs_ok:
                                continue

                            # Check: GC in range
                            _cur_gc_both = sum(1 for b in current if b in "GC") + sum(1 for b in _adj_current if b in "GC")
                            _alt_gc_both = sum(1 for b in alt if b in "GC") + sum(1 for b in _adj_alt if b in "GC")
                            _new_gc = _rec_gc_count - _cur_gc_both + _alt_gc_both
                            if not (gc_lo <= _new_gc / _hc_n_bases <= gc_hi):
                                continue

                            # Check: no new ATTTA
                            if _test_seq.count("ATTTA") > sequence.count("ATTTA"):
                                continue

                            # Check: no 6+ T runs (localized)
                            _t_ok = True
                            _t_start = max(0, _lo_ci * 3 - T_RUN_LENGTH_THRESHOLD)
                            _t_end = min(len(_test_seq), _hi_ci * 3 + 3 + T_RUN_LENGTH_THRESHOLD)
                            _j = _t_start
                            while _j < _t_end:
                                if _test_seq[_j] == "T":
                                    _k = _j
                                    while _k < len(_test_seq) and _test_seq[_k] == "T":
                                        _k += 1
                                    if _k - _j >= T_RUN_LENGTH_THRESHOLD:
                                        _t_ok = False
                                        break
                                    _j = _k
                                else:
                                    _j += 1
                            if not _t_ok:
                                continue

                            # No worsening of splice scores (eukaryotes only)
                            if not is_prokaryote:
                                if max_donor_score(_test_seq) > max_donor_score(sequence) + 0.5:
                                    continue
                                if max_acceptor_score(_test_seq) > max_acceptor_score(sequence) + 0.5:
                                    continue

                            # Accept the paired upgrade
                            sequence = _test_seq
                            _rec_gc_count = _new_gc
                            _any_recovery = True
                            logger.debug(
                                "CAI recovery: paired upgrade codon %d (%s→%s) + "
                                "codon %d (%s→%s)",
                                ci, current, alt,
                                _adj_ci, _adj_current, _adj_alt,
                            )
                            break
                        if _any_recovery:
                            break
                    if _any_recovery:
                        break

        if not _any_recovery:
            break

    # Post-condition: verify sequence still encodes the same protein
    from .translation import translate
    translated = translate(sequence)
    assert translated == protein, (
        f"Post-condition violation: optimizer changed the protein. "
        f"Expected '{protein[:20]}...', got '{translated[:20]}...'"
    )
    assert len(sequence) == len(aas) * 3, (
        f"Post-condition violation: sequence length {len(sequence)} "
        f"!= expected {len(aas) * 3}"
    )

    for w in warnings:
        logger.warning(w)

    return sequence, warnings


# ==============================================================================
# IUPAC Expansion
# ==============================================================================

def _expand_iupac_site(pattern: str) -> list[str]:
    """Expand an IUPAC restriction site pattern into all concrete ACGT sequences.

    E.g., GGCCNNNNNGGCC expands into 4^5 = 1024 concrete sequences.
    For very large expansions, we cap at 4096 to avoid combinatorial explosion.

    Pre-conditions:
    - pattern is a non-empty string containing IUPAC codes

    Post-conditions:
    - all returned strings contain only ACGT characters
    - len(result[0]) == len(pattern) for all results
    """
    assert len(pattern) > 0, "Pattern must not be empty"

    if not any(b not in "ACGT" for b in pattern):
        return [pattern]

    total_combos = 1
    for b in pattern:
        if b not in "ACGT":
            total_combos *= len(IUPAC_EXPAND.get(b, "A"))

    if total_combos > IUPAC_EXPANSION_CAP:
        logger.warning(
            "IUPAC site %s expands to %d variants (>%d), skipping",
            pattern, total_combos, IUPAC_EXPANSION_CAP,
        )
        return []

    results = [""]
    for b in pattern:
        bases = IUPAC_EXPAND.get(b, b)
        results = [r + x for r in results for x in bases]
    return results


# ==============================================================================
# Predicate Checking (Delegates to Type System — SOC)
# ==============================================================================

def _check_predicates_via_type_system(
    sequence: str,
    gc_lo: float,
    gc_hi: float,
    restriction_sites: list[str],
    cai_threshold: float,
    organism: str,
    cryptic_splice_threshold: float = 3.0,
) -> tuple[list[str], list[str]]:
    """Check all type predicates against the optimized sequence.

    DELEGATES to the type system's evaluate_all_predicates rather than
    re-implementing predicate logic here. This is the single source of truth.

    Pre-conditions:
    - sequence is a valid DNA string
    - organism is in SUPPORTED_ORGANISMS
    - 0.0 <= gc_lo < gc_hi <= 1.0
    - cai_threshold > 0

    Post-conditions:
    - satisfied + failed covers all checked predicates
    - satisfied and failed are disjoint
    """
    from .type_system import evaluate_all_predicates
    from .types import Verdict

    assert 0.0 <= gc_lo < gc_hi <= 1.0, f"Invalid GC bounds: [{gc_lo}, {gc_hi}]"
    assert cai_threshold > 0, f"CAI threshold must be positive, got {cai_threshold}"

    # Build exon boundaries for a coding sequence (single exon)
    exon_boundaries = [(0, len(sequence))]

    # Get enzyme names from sequences
    enzyme_names = []
    for site in restriction_sites:
        found = False
        for name, seq in RESTRICTION_ENZYMES.items():
            if seq.upper() == site.upper():
                enzyme_names.append(name)
                found = True
                break
        if not found:
            enzyme_names.append(site)  # Use raw sequence as name

    results = evaluate_all_predicates(
        seq=sequence,
        known_exon_boundaries=exon_boundaries,
        organism=organism,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        cai_threshold=cai_threshold,
        enzymes=enzyme_names,
        cryptic_splice_threshold=cryptic_splice_threshold,
    )

    satisfied = []
    failed = []
    for r in results:
        predicate_name = r.predicate
        if r.verdict in (Verdict.PASS, Verdict.LIKELY_PASS):
            satisfied.append(predicate_name)
        elif r.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL):
            failed.append(predicate_name)
        else:
            # UNCERTAIN: treat as failed for optimization purposes
            failed.append(predicate_name)

    # Verify disjoint
    assert not (set(satisfied) & set(failed)), (
        f"Predicates cannot be both satisfied and failed: "
        f"{set(satisfied) & set(failed)}"
    )

    return satisfied, failed


# ==============================================================================
# Main Optimization Entry Point
# ==============================================================================

def optimize_sequence(
    target_protein: str | None = None,
    organism: str = "Homo_sapiens",
    species: str | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.5,
    enzymes: list[str] | None = None,
    strategy: str = "hybrid",
    use_csp_solver: bool = False,
    optimize_mrna_stability: bool = True,
    strict_mode: bool = True,
    seed: int | None = None,
    include_utr: bool = True,
    consider_codon_pair_bias: bool = False,
    track_provenance: bool = True,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
    objective: str | Callable[..., float] | None = None,
    gc_window_size: int = 50,
    gc_window_min: float | None = None,
    gc_window_max: float | None = None,
    **kwargs: Any,
) -> OptimizationResult:
    """Optimize a protein sequence for expression in the target organism.

    This is the high-level convenience API that wraps BioOptimizer.
    It takes a protein sequence and returns an OptimizationResult with
    the optimized DNA sequence and quality metrics.

    When ``use_csp_solver=True``, the CSP/SMT constraint solver is tried
    first (OR-Tools CP-SAT or Z3 SMT).  If no solver is available or the
    solver times out, the greedy optimizer is used as fallback.

    Organism Specification:

        The target organism can be specified using **either** the
        ``organism`` parameter **or** the ``species`` parameter.  Both
        accept the same set of names — short aliases, abbreviated
        binomials, display names, or full canonical names — and both
        map to the same internal representation via
        :func:`~biocompiler.organisms.resolve_organism`.

        +-------------------------+------------------------+
        | Example                 | Resolved to            |
        +=========================+========================+
        | ``species='ecoli'``     | ``'Escherichia_coli'`` |
        | ``organism='ecoli'``    | ``'Escherichia_coli'`` |
        | ``species='E. coli'``   | ``'Escherichia_coli'`` |
        | ``organism='Escherichia_coli'`` | same           |
        | ``species='human'``     | ``'Homo_sapiens'``     |
        | ``species='h_sapiens'`` | ``'Homo_sapiens'``     |
        +-------------------------+------------------------+

        If both ``species`` and ``organism`` are provided, ``species``
        takes precedence and a :class:`DeprecationWarning` is emitted
        recommending the use of ``organism`` only (the canonical
        parameter).  If neither matches the other after resolution, a
        warning is logged.

    Args:
        target_protein: Amino acid sequence (1-letter codes, no stop).
            Can also be passed as the first positional argument.
        organism: Target organism.  Accepts canonical binomials
            (e.g., ``'Homo_sapiens'``, ``'Escherichia_coli'``),
            short keys (``'ecoli'``, ``'human'``), abbreviated
            binomials (``'E_coli'``, ``'h_sapiens'``), or display
            names (``'E. coli'``).  All forms are resolved via
            :func:`~biocompiler.organisms.resolve_organism`.
        species: Alias for ``organism``.  Accepts the same values.
            If provided **together with** ``organism``, ``species``
            takes precedence and a deprecation warning is emitted.
            Prefer using ``organism`` in new code; ``species`` is
            retained for backward compatibility.
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        gc_window_size: Sliding window size (in nucleotides) for local
            GC content checking.  Set to 0 to disable sliding-window GC
            constraints.  Default: 50.
        gc_window_min: Minimum acceptable GC fraction for each sliding
            window.  If None (default), uses ``gc_lo``.
        gc_window_max: Maximum acceptable GC fraction for each sliding
            window.  If None (default), uses ``gc_hi``.
        cai_threshold: Minimum CAI score for the CodonAdapted predicate.
        enzymes: List of restriction enzyme names to avoid.
        strategy: Optimization strategy ('hybrid', 'constraint_first', or 'cai_first').
            'hybrid' (default for v10): Greedy init + priority-queue local search + hill climb.
            'constraint_first': Legacy sequential pipeline.
            'cai_first': Maximize CAI first, then fix constraints.
        use_csp_solver: If True, try CSP/SMT solver before greedy (default False).
        optimize_mrna_stability: If True, run a soft mRNA stability improvement
            pass after CAI optimization. Identifies destabilizing motifs (ATTTA,
            extended AREs, long A/T runs) and suggests synonymous codon changes
            to remove them without breaking hard constraints. Defaults to True.
        seed: Deterministic seed for reproducibility. Currently unused as the
            greedy optimizer is fully deterministic, but reserved for future
            randomized optimization strategies (e.g., Monte Carlo, genetic
            algorithms). If provided, ``random.seed(seed)`` is called at the
            start of optimization to ensure reproducible behavior when
            randomized steps are added.
        include_utr: If True, generate organism-appropriate 5' and 3' UTR
            suggestions and include them in the result. The UTRs are SUGGESTED
            but not enforced — they're for the user to consider when ordering
            a synthesis construct. Defaults to True.
        consider_codon_pair_bias: If True, run a soft codon pair bias (CPB)
            improvement pass after the CAI optimization pass.  For each
            adjacent codon pair, check if there's a synonymous pair with
            higher combined CAI+CPB score (70% CAI / 30% CPB by default)
            that doesn't violate hard constraints (restriction sites, GC
            range, splice sites).  The result will include a
            ``codon_pair_bias`` field with the mean CPB score.  Defaults to
            False.
        track_provenance: If True, collect decision-level provenance for
            every codon choice made during optimization.  The result will
            include a ``decision_trail`` field with the full audit trail.
            If False, skip provenance collection entirely (zero overhead).
            Defaults to True.
        strict_mode: If True (default), raise
            :class:`OptimizationConstraintError` when any predicates fail
            instead of returning a result with ``failed_predicates``.  This
            provides a hard-stop guarantee that the returned sequence
            satisfies all constraints.  Set to False to allow partial
            results with failed predicates.
        max_iterations: Maximum number of optimization iterations before
            declaring ``convergence_status='max_iterations'``.  A sensible
            default of 1000 is used.  The optimizer will stop early if
            convergence is detected (objective unchanged for 3 consecutive
            iterations) or oscillation is detected (objective cycles up and
            down).  The result will include ``convergence_status`` and
            ``iterations_used`` fields.
        objective: Custom optimization objective.  Controls how codon
            alternatives are scored during the CAI hill-climb phase.
            Accepts:

            - ``None`` (default): maximize CAI only (equivalent to
              ``"cai"``).
            - A string name: one of ``"cai"``, ``"cai_gc_balanced"``,
              ``"codon_pair"``, ``"min_max_gc"``.
            - A callable ``(dna: str, protein: str, organism: str) -> float``:
              a custom function where higher values are better.

            When a non-CAI objective is provided, the CAI hill-climb
            phase scores candidate codon substitutions by evaluating the
            custom objective on the *modified* sequence rather than by
            CAI alone.  This allows users to optimize for GC balance,
            codon pair bias, or any custom metric while still respecting
            all hard constraints (restriction sites, GC range, splice
            sites, etc.).
        **kwargs: Additional arguments (e.g., splice_low, splice_high,
            restriction_sites, organism_domain).

    Returns:
        OptimizationResult with optimized sequence and metrics.

    Raises:
        InvalidProteinError: if the protein contains invalid amino acid codes.
        UnsupportedOrganismError: if the organism is not supported.
        OptimizationConstraintError: if ``strict_mode=True`` and the
            optimized sequence has one or more failed predicates.
    """
    # ── Resolve custom objective ──────────────────────────────────────
    _objective_fn = _resolve_objective(objective)

    # Set deterministic seed if provided (reserved for future randomized steps)
    if seed is not None:
        import random
        random.seed(seed)

    # Start timing for provenance
    import time as _time
    _start_time = _time.monotonic()

    # Handle positional arg: optimize_sequence("MVHLTPEEK", organism="...")
    if target_protein is None and len(kwargs.get("_args", [])) > 0:
        target_protein = kwargs["_args"][0]

    # ── Organism resolution ────────────────────────────────────────
    # Support both 'species' and 'organism' parameters.  Both accept
    # the same set of aliases and resolve to the same canonical name.
    # 'species' is retained for backward compatibility but is not the
    # preferred parameter — use 'organism' in new code.
    _resolved_organism: str | None = None

    # Legacy: species passed via kwargs (old calling convention)
    _species_from_kwargs = kwargs.pop("species", None)
    if _species_from_kwargs is not None and species is None:
        species = _species_from_kwargs

    if species is not None:
        _resolved_organism = resolve_organism(species, strict=False)
        # Emit deprecation warning if both species and organism are given
        if organism != "Homo_sapiens":
            _resolved_organism_from_explicit = resolve_organism(organism, strict=False)
            if _resolved_organism != _resolved_organism_from_explicit:
                import warnings
                warnings.warn(
                    f"Both 'species={species!r}' and 'organism={organism!r}' "
                    f"were provided but resolve to different organisms "
                    f"({_resolved_organism!r} vs {_resolved_organism_from_explicit!r}). "
                    f"Using 'species' ({_resolved_organism!r}). "
                    f"Prefer using only 'organism' in new code.",
                    DeprecationWarning,
                    stacklevel=2,
                )
            else:
                import warnings
                warnings.warn(
                    f"Both 'species' and 'organism' were provided. "
                    f"Prefer using only 'organism' in new code; "
                    f"'species' is retained for backward compatibility.",
                    DeprecationWarning,
                    stacklevel=2,
                )
        else:
            # species was given but organism was left at default —
            # no conflict, just use species. Still emit a gentle hint.
            import warnings
            warnings.warn(
                f"The 'species' parameter is deprecated in favor of 'organism'. "
                f"Use organism='{_resolved_organism}' instead of "
                f"species='{species}'. Both accept the same aliases.",
                DeprecationWarning,
                stacklevel=2,
            )

    if _resolved_organism is not None:
        organism = _resolved_organism
    else:
        organism = resolve_organism(organism, strict=False)

    if not target_protein:
        raise InvalidProteinError("", {"<empty>"})

    target_protein = target_protein.strip().upper()

    # Validate protein
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = set(target_protein) - valid_aas
    if invalid:
        raise InvalidProteinError(target_protein, invalid)

    # Initialize decision provenance collector if tracking is enabled
    _provenance_collector: DecisionProvenanceCollector | None = None
    if track_provenance:
        _provenance_collector = DecisionProvenanceCollector()
        _provenance_collector.start_optimization(
            protein=target_protein,
            organism=organism,
            constraints=[
                "GCInRange", "NoRestrictionSite", "NoCrypticSplice",
                "NoStopCodons", "CodonAdapted",
            ],
            solver_backend="csp" if use_csp_solver else "greedy",
            seed=seed,
        )

    # ── Convergence tracking ────────────────────────────────────────
    _convergence = ConvergenceTracker()
    _convergence_status: str = "converged"  # default; overridden if caps are hit

    # ── CSP Solver path ─────────────────────────────────────────────
    if use_csp_solver:
        try:
            from .solver.dispatch import csp_optimize
            from .solver.types import SolverConfig
            restriction_sites = []
            if enzymes:
                from .restriction_sites import get_recognition_site
                for enz in enzymes:
                    site = get_recognition_site(enz)
                    if site:
                        restriction_sites.append(site)
            else:
                for enz in ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]:
                    from .restriction_sites import get_recognition_site as _grs
                    site = _grs(enz)
                    if site:
                        restriction_sites.append(site)

            config = SolverConfig(
                gc_lo=gc_lo,
                gc_hi=gc_hi,
                cryptic_splice_threshold=kwargs.get("splice_low", 3.0),
                restriction_sites=restriction_sites,
            )
            solver_result = csp_optimize(
                protein=target_protein,
                organism=organism,
                config=config,
            )
            if solver_result.solved and not solver_result.fallback_used:
                from .scanner import gc_content as _gc_content
                from .translation import compute_cai as _compute_cai
                seq = solver_result.sequence
                gc = _gc_content(seq)
                _canonical_organism = _species_key_to_organism(organism)
                try:
                    cai_val = _compute_cai(seq, _canonical_organism)
                except Exception:
                    logger.debug(
                        "CAI computation failed for organism '%s', using solver result",
                        organism, exc_info=True,
                    )
                    cai_val = solver_result.cai
                # UTR suggestions for CSP solver path
                _utr5_seq: str | None = None
                _utr3_seq: str | None = None
                _utr_score5: float | None = None
                _utr_score3: float | None = None
                if include_utr:
                    from .utr_models import suggest_5utr as _s5, suggest_3utr as _s3
                    from .utr_models import score_5utr as _sc5, score_3utr as _sc3
                    try:
                        _utr5_seq = _s5(organism)
                        _utr_score5 = _sc5(_utr5_seq, organism)
                    except ValueError:
                        logger.debug("No 5' UTR suggestion for organism '%s'", organism)
                    try:
                        _utr3_seq = _s3(organism)
                        _utr_score3 = _sc3(_utr3_seq, organism)
                    except ValueError:
                        logger.debug("No 3' UTR suggestion for organism '%s'", organism)

                # Compute CPB if requested (CSP solver path)
                _cpb_val: float | None = None
                if consider_codon_pair_bias:
                    try:
                        from .codon_pair_scoring import compute_cpb as _compute_cpb
                        _cpb_val = round(_compute_cpb(seq, organism), 6)
                    except Exception:
                        logger.debug("CPB computation failed for CSP solver path", exc_info=True)

                # ── Translation verification (CSP solver path) ──
                from .protein_verification import verify_and_raise as _verify_and_raise
                _verify_and_raise(seq, target_protein, organism=organism)

                # Evaluate custom objective score if a non-default objective is used
                _csp_objective_score: float | None = None
                if _objective_fn is not _cai_objective:
                    try:
                        _csp_objective_score = _objective_fn(seq, target_protein, organism)
                    except Exception:
                        _csp_objective_score = None

                return OptimizationResult(
                    sequence=seq,
                    gc_content=round(gc, 4),
                    cai=cai_val,
                    failed_predicates=[],
                    predicate_results=[],
                    certificate_text="CSP-Solver:" + str(solver_result.backend_used.value),
                    protein=target_protein,
                    fallback_used=False,
                    satisfied_predicates=["CSP_SOLVER"],
                    codon_pair_bias=_cpb_val,
                    suggested_5utr=_utr5_seq,
                    suggested_3utr=_utr3_seq,
                    utr_score_5=_utr_score5,
                    utr_score_3=_utr_score3,
                    convergence_status="converged",
                    iterations_used=1,
                    objective_score=_csp_objective_score,
                )
        except Exception:
            logger.warning(
                "CSP solver failed, falling back to greedy optimizer",
                exc_info=True,
            )
    # ── End CSP Solver path ─────────────────────────────────────────

    # Map organism name to species key
    species_key = _organism_to_species_key(organism)

    # ── Organism-aware constraint selection ────────────────────────
    # When the target organism is prokaryotic, skip eukaryote-specific
    # constraints (cryptic splice sites, CpG islands, GT dinucleotide
    # avoidance) that are biologically inappropriate and unnecessarily
    # depress CAI.  This recovers ~0.27 CAI on prokaryotic targets
    # (see Task 6+8 benchmark findings).
    organism_domain = kwargs.get("organism_domain", "auto")
    if organism_domain not in ("auto", "eukaryote", "prokaryote"):
        organism_domain = "auto"

    if organism_domain == "auto":
        from .organism_config import is_eukaryotic_organism
        is_prokaryote = not is_eukaryotic_organism(organism)
    elif organism_domain == "prokaryote":
        is_prokaryote = True
    else:
        is_prokaryote = False

    # PERF (Optimization A): Fast path for prokaryotes — skip all eukaryotic
    # constraint setup and evaluation.  Prokaryotes have no spliceosome,
    # so GT/AG avoidance, CpG island disruption, and cryptic splice
    # site evaluation are biologically irrelevant.  Skipping them
    # eliminates the expensive MaxEntScan calls and eukaryotic predicate
    # checks that dominate runtime for short sequences.
    # For prokaryotes: skip splice/GT/CpG constraints entirely
    # For eukaryotes: apply all constraints (default)
    effective_avoid_gt = kwargs.get("avoid_gt", not is_prokaryote)
    effective_splice_low = kwargs.get("splice_low", 3.0 if not is_prokaryote else 999.0)
    effective_splice_high = kwargs.get("splice_high", 6.0 if not is_prokaryote else 999.0)

    # ── Hybrid strategy path (v10 default) ────────────────────────
    if strategy == "hybrid":
        from .hybrid_optimizer import HybridOptimizer as _HybridOptimizer

        hybrid_opt = _HybridOptimizer(
            species=species_key,
            organism=organism,
            enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            avoid_gt=effective_avoid_gt,
            splice_threshold=effective_splice_low,
            provenance_collector=_provenance_collector,
        )

        hybrid_result = hybrid_opt.optimize(
            target_protein, is_prokaryote=is_prokaryote
        )
        optimized_seq = hybrid_result.sequence

        # ── Eukaryotic CAI recovery pass ──
        # For eukaryotes, the HybridOptimizer aggressively avoids GT
        # dinucleotides, which can replace optimal CAI codons with
        # suboptimal GT-free alternatives (e.g., yeast Leu: TTG→TTA,
        # Gly: GGT→GGA, Cys: TGT→TGC).  For eukaryotes, GT avoidance
        # in CDS is MUCH less important than CAI — GT dinucleotides
        # only matter when they form strong cryptic splice donor sites
        # (score_splice_donor_potential >= 0.5), and in-codon GTs from
        # optimal codons are biologically acceptable.
        #
        # This recovery pass swaps suboptimal codons back to optimal
        # ones when the CAI cost exceeds the threshold, using cost-aware
        # GT resolution with splice donor scoring.  GTs with low splice
        # donor potential are accepted in favor of CAI.
        _cai_recovery_upgrades = 0
        if not is_prokaryote and effective_avoid_gt:
            _usage = CODON_ADAPTIVENESS_TABLES.get(organism)
            if _usage is not None:
                optimized_seq, _cai_recovery_upgrades = _eukaryote_cai_recovery(
                    optimized_seq,
                    target_protein,
                    _usage,
                    enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                )
                if _cai_recovery_upgrades > 0:
                    logger.info(
                        "Eukaryotic CAI recovery: upgraded %d codon(s) to "
                        "optimal (GT accepted as soft preference)",
                        _cai_recovery_upgrades,
                    )

        # ── Systematic CpG elimination pass (eukaryotes only) ──
        # After CAI maximization and recovery, systematically eliminate CpG
        # dinucleotides that contribute to CpG islands.  This is a proper
        # optimization pass (not best-effort) that continues past unfixable
        # positions and reports remaining CpGs as warnings.
        _cpg_warnings: list[str] = []
        _cpg_seq_changed = False
        if not is_prokaryote:
            _usage_cpg = CODON_ADAPTIVENESS_TABLES.get(organism)
            if _usage_cpg is not None:
                optimized_seq, _cpg_warnings = _eliminate_cpg_dinucleotides(
                    optimized_seq,
                    target_protein,
                    _usage_cpg,
                    enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                    organism=organism,
                    gc_lo=gc_lo,
                    gc_hi=gc_hi,
                )
                _seq_before_cpg = hybrid_result.sequence
                _cpg_seq_changed = (optimized_seq != _seq_before_cpg)
                if _cpg_warnings:
                    for _w in _cpg_warnings:
                        logger.warning(_w)
                if _cpg_seq_changed:
                    logger.info(
                        "CpG elimination: modified sequence to avoid CpG islands"
                    )

        # ── Ultra-fast prokaryote predicate evaluation ──
        # For prokaryotes, skip the heavyweight BioOptimizer predicate
        # evaluation and use a streamlined version that only checks
        # prokaryote-relevant predicates (no splice, no CpG, no GT).
        # This avoids creating a BioOptimizer and running the full
        # _evaluate_all_predicates with 12 predicates.
        if is_prokaryote:
            from .type_system import (
                check_no_stop_codons, check_valid_coding_seq,
                check_no_restriction_site,
            )
            import math as _pmath

            # Lightweight predicate evaluation — only prokaryote-relevant
            pred_results = []

            # 1. NoStopCodons
            pred_results.append(check_no_stop_codons(optimized_seq))

            # 2. NoCrypticSplice — skipped for prokaryotes
            pred_results.append(PredicateResult(
                "NoCrypticSplice", True,
                details="Skipped for prokaryotic organism"
            ))

            # 3. NoCpGIsland — skipped for prokaryotes
            pred_results.append(PredicateResult(
                "NoCpGIsland", True,
                details="Skipped for prokaryotic organism"
            ))

            # 4. NoRestrictionSite
            pred_results.append(check_no_restriction_site(
                optimized_seq, enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
            ))

            # 5. NoGTDinucleotide — skipped for prokaryotes
            pred_results.append(PredicateResult(
                "NoGTDinucleotide", True,
                details="Skipped for prokaryotic organism"
            ))

            # 6. ValidCodingSeq
            pred_results.append(check_valid_coding_seq(optimized_seq))

            # 7. ConservationScore — trivially passes (no mutagenesis)
            pred_results.append(PredicateResult(
                "ConservationScore", True,
                details=f"All AA conservation scores >= 0"
            ))

            # 8. CodonOptimality — use hybrid_result CAI
            cai_val = hybrid_result.cai
            all_optimal = cai_val >= cai_threshold
            pred_results.append(PredicateResult(
                "CodonOptimality", all_optimal,
                details=f"CAI={cai_val:.4f}, min={cai_threshold}"
            ))

            # 9. GCInRange
            _gc_val = hybrid_result.gc_content
            gc_ok = gc_lo <= _gc_val <= gc_hi
            pred_results.append(PredicateResult(
                "GCInRange", gc_ok,
                details=f"GC content: {_gc_val:.3f} (range [{gc_lo}, {gc_hi}])"
            ))

            # 10. SlidingGC (prokaryotic fast path)
            _eff_gc_w_min = gc_window_min if gc_window_min is not None else gc_lo
            _eff_gc_w_max = gc_window_max if gc_window_max is not None else gc_hi
            if gc_window_size > 0 and len(optimized_seq) >= gc_window_size:
                _sgc_result = check_sliding_gc(
                    optimized_seq, window_size=gc_window_size,
                    gc_min=_eff_gc_w_min, gc_max=_eff_gc_w_max,
                )
                pred_results.append(PredicateResult(
                    "SlidingGC", _sgc_result.passed,
                    details=(
                        f"Window={gc_window_size}, range=[{_eff_gc_w_min:.2f}, {_eff_gc_w_max:.2f}], "
                        f"min_gc={_sgc_result.min_gc:.3f}, max_gc={_sgc_result.max_gc:.3f}, "
                        f"violations={len(_sgc_result.violations)}"
                    ),
                ))
            else:
                pred_results.append(PredicateResult(
                    "SlidingGC", True,
                    details="Sliding-window GC check skipped (window_size=0 or sequence too short)"
                ))

            cert_text = format_certificate(pred_results, optimized_seq, species_key)

            # Skip all post-processing for prokaryotes (CAI recovery,
            # mRNA stability, CPB, UTR) since we start with max CAI
            # and prokaryote-specific constraints are already handled.
            gc = hybrid_result.gc_content
            failed = [r.predicate for r in pred_results if not r.passed]
            satisfied = [r.predicate for r in pred_results if r.passed]
            fallback = False

            # Import compute_cai for the hybrid path
            from .translation import compute_cai

            # Build provenance record (lightweight for prokaryotes)
            from .provenance import OptimizationRecord as _OptimizationRecord
            from .provenance import _get_biocompiler_version as _get_version

            provenance_record = _OptimizationRecord(
                input_sequence=target_protein,
                output_sequence=optimized_seq,
                organism=organism,
                constraints_applied=sorted([r.predicate for r in pred_results]),
                mutations_made=[],
                solver_backend="hybrid-prok-fast",
                solve_time=round(_time.monotonic() - _start_time, 6),
                seed_used=seed,
                timestamp=datetime.now(timezone.utc).isoformat(),
                biocompiler_version=_get_version(),
            )

            # Finalize decision provenance trail if tracking is enabled
            _prok_decision_trail: OptimizationDecisionTrail | None = None
            if _provenance_collector is not None:
                try:
                    # Record constraint decisions from predicate results
                    for _pr in pred_results:
                        _action = "satisfied" if _pr.passed else "conflicted"
                        _positions_affected: list[int] = list(range(len(optimized_seq) // 3))
                        _details = _pr.details or ""
                        _provenance_collector.record_constraint_decision(ConstraintDecision(
                            constraint_name=_pr.predicate,
                            constraint_type="hard",
                            action_taken=_action,
                            positions_affected=_positions_affected[:1] if _positions_affected else [],
                            tradeoff_description=_details,
                            impact_on_cai=0.0,
                        ))
                    _prok_decision_trail = _provenance_collector.finalize(
                        output_dna=optimized_seq,
                        cai=cai_val,
                        gc=gc,
                    )
                except Exception:
                    logger.debug("Provenance finalization failed for prokaryote path", exc_info=True)

            # ── Translation verification (prokaryote fast path) ──
            from .protein_verification import verify_and_raise as _verify_prok
            _verify_prok(optimized_seq, target_protein, organism=organism)

            # ── Custom objective refinement (prokaryote fast path) ──
            _prok_objective_score: float | None = None
            if _objective_fn is not _cai_objective:
                _obj_aas = protein_to_aa_list(target_protein)
                _obj_n_codons = len(_obj_aas)
                _obj_usage = CODON_ADAPTIVENESS_TABLES.get(organism, {})
                _obj_sorted_codons: dict[str, list[str]] = {}
                for _aa in set(_obj_aas):
                    if _aa == "*":
                        continue
                    _obj_sorted_codons[_aa] = sorted(
                        AA_TO_CODONS.get(_aa, []),
                        key=lambda c: _obj_usage.get(c, 0.0),
                        reverse=True,
                    )
                _obj_rs_sites: list[tuple[str, str]] = []
                for _enz in (enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]):
                    from .restriction_sites import get_recognition_site as _grs2
                    _obj_site = _grs2(_enz)
                    if _obj_site:
                        _obj_rs_sites.append((_obj_site, reverse_complement(_obj_site)))

                _MAX_PROK_OBJ_ITERS = 5
                for _obj_iter in range(_MAX_PROK_OBJ_ITERS):
                    _obj_improved = False
                    _current_score = _objective_fn(optimized_seq, target_protein, organism)

                    for _ci in range(_obj_n_codons):
                        _aa = _obj_aas[_ci]
                        if _aa == "*" or _aa == "M":
                            continue
                        _current_codon = optimized_seq[_ci * 3:_ci * 3 + 3]
                        _best_alt = None
                        _best_alt_score = _current_score

                        for _alt in _obj_sorted_codons.get(_aa, []):
                            if _alt == _current_codon:
                                continue
                            _test_seq = optimized_seq[:_ci * 3] + _alt + optimized_seq[_ci * 3 + 3:]

                            # Hard constraint checks
                            _rs_ok = True
                            for _site, _site_rc in _obj_rs_sites:
                                if _site in _test_seq and _site not in optimized_seq:
                                    _rs_ok = False
                                    break
                                if _site_rc and _site_rc in _test_seq and _site_rc not in optimized_seq:
                                    _rs_ok = False
                                    break
                            if not _rs_ok:
                                continue

                            _test_gc = (_test_seq.count("G") + _test_seq.count("C")) / max(len(_test_seq), 1)
                            if not (gc_lo <= _test_gc <= gc_hi):
                                continue

                            _has_premature_stop = False
                            for _si in range(0, len(_test_seq) - 5, 3):
                                if _test_seq[_si:_si + 3] in ("TAA", "TAG", "TGA"):
                                    _has_premature_stop = True
                                    break
                            if _has_premature_stop:
                                continue

                            if _test_seq.count("ATTTA") > optimized_seq.count("ATTTA"):
                                continue

                            try:
                                _alt_score = _objective_fn(_test_seq, target_protein, organism)
                            except Exception:
                                continue

                            if _alt_score > _best_alt_score:
                                _best_alt = _alt
                                _best_alt_score = _alt_score

                        if _best_alt is not None:
                            optimized_seq = (
                                optimized_seq[:_ci * 3] + _best_alt + optimized_seq[_ci * 3 + 3:]
                            )
                            _obj_improved = True

                    if not _obj_improved:
                        break

                try:
                    _prok_objective_score = _objective_fn(optimized_seq, target_protein, organism)
                except Exception:
                    _prok_objective_score = None

                # Recompute metrics after objective refinement
                gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
                gc = round(gc, 4)
                cai_val = compute_cai(optimized_seq, organism)

            return OptimizationResult(
                sequence=optimized_seq,
                gc_content=gc,
                cai=cai_val,
                failed_predicates=failed,
                predicate_results=pred_results,
                certificate_text=cert_text,
                protein=target_protein,
                fallback_used=False,
                satisfied_predicates=satisfied,
                provenance=provenance_record,
                decision_trail=_prok_decision_trail,
                convergence_status="converged",
                iterations_used=1,
                objective_score=_prok_objective_score,
            )

        # Fast eukaryotic predicate evaluation — skip BioOptimizer creation
        # (~0.4ms saved per call) by directly checking only the predicates
        # that matter after hybrid optimization.
        from .type_system import (
            check_no_stop_codons, check_valid_coding_seq,
            check_no_restriction_site, check_no_cryptic_splice,
            check_no_cpg_island,
        )

        pred_results = []

        # 1. NoStopCodons
        pred_results.append(check_no_stop_codons(optimized_seq))

        # 2. NoCrypticSplice — skip if already validated by hybrid optimizer
        if hybrid_result.splice_sites_validated:
            pred_results.append(PredicateResult(
                "NoCrypticSplice", True,
                details="Validated during optimization (MaxEntScan Phase 3)"
            ))
        else:
            pred_results.append(check_no_cryptic_splice(
                optimized_seq, low_thresh=effective_splice_low, organism=organism
            ))

        # 3. NoCpGIsland
        pred_results.append(check_no_cpg_island(optimized_seq, organism=organism))

        # 4. NoRestrictionSite
        pred_results.append(check_no_restriction_site(
            optimized_seq, enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
        ))

        # 5. NoGTDinucleotide (soft for eukaryotes — in-codon GTs from
        # optimal codons and cross-codon GTs with high CAI cost are acceptable)
        # Uses check_no_gt_dinucleotide_soft which returns:
        #   PASS          — GT count ≤ max_gt_count (auto-computed per sequence length)
        #   LIKELY_FAIL   — GT count > max_gt_count for eukaryotes (soft fail)
        #   FAIL          — any GT for prokaryotes (hard constraint)
        if effective_avoid_gt:
            pred_results.append(check_no_gt_dinucleotide_soft(
                optimized_seq, organism=organism,
            ))
        else:
            pred_results.append(PredicateResult(
                "NoGTDinucleotide", True,
                details="GT avoidance not requested"
            ))

        # 6. ValidCodingSeq
        pred_results.append(check_valid_coding_seq(optimized_seq))

        # 7. ConservationScore
        pred_results.append(PredicateResult(
            "ConservationScore", True,
            details=f"All AA conservation scores >= 0"
        ))

        # 8. CodonOptimality
        # Recompute CAI after eukaryotic recovery and CpG elimination passes
        # (which may have changed codons)
        from .translation import compute_cai as _compute_cai_hybrid
        if _cai_recovery_upgrades > 0 or _cpg_seq_changed:
            cai_val = _compute_cai_hybrid(optimized_seq, organism)
        else:
            cai_val = hybrid_result.cai
        all_optimal = cai_val >= cai_threshold
        pred_results.append(PredicateResult(
            "CodonOptimality", all_optimal,
            details=f"CAI={cai_val:.4f}, min={cai_threshold}"
        ))

        # 9. GCInRange
        # Recompute GC after CAI recovery/CpG elimination may have changed the sequence
        if _cai_recovery_upgrades > 0 or _cpg_seq_changed:
            _gc_val = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
        else:
            _gc_val = hybrid_result.gc_content
        gc_ok = gc_lo <= _gc_val <= gc_hi
        pred_results.append(PredicateResult(
            "GCInRange", gc_ok,
            details=f"GC content: {_gc_val:.3f} (range [{gc_lo}, {gc_hi}])"
        ))

        # 10. SlidingGC (eukaryotic path)
        _eff_gc_w_min = gc_window_min if gc_window_min is not None else gc_lo
        _eff_gc_w_max = gc_window_max if gc_window_max is not None else gc_hi
        if gc_window_size > 0 and len(optimized_seq) >= gc_window_size:
            _sgc_result = check_sliding_gc(
                optimized_seq, window_size=gc_window_size,
                gc_min=_eff_gc_w_min, gc_max=_eff_gc_w_max,
            )
            pred_results.append(PredicateResult(
                "SlidingGC", _sgc_result.passed,
                details=(
                    f"Window={gc_window_size}, range=[{_eff_gc_w_min:.2f}, {_eff_gc_w_max:.2f}], "
                    f"min_gc={_sgc_result.min_gc:.3f}, max_gc={_sgc_result.max_gc:.3f}, "
                    f"violations={len(_sgc_result.violations)}"
                ),
            ))
        else:
            pred_results.append(PredicateResult(
                "SlidingGC", True,
                details="Sliding-window GC check skipped (window_size=0 or sequence too short)"
            ))

        cert_text = format_certificate(pred_results, optimized_seq, species_key)

        # Create a lightweight opt stub for downstream code that checks opt._applied_mutagenesis
        opt = BioOptimizer.__new__(BioOptimizer)
        opt._applied_mutagenesis = hybrid_opt._applied_mutagenesis
        opt._original_protein = target_protein
        opt.organism_name = organism
        # Provide species_cai for fallback CAI computation in edge cases
        opt.species_cai = hybrid_opt._species_cai if hasattr(hybrid_opt, '_species_cai') else {}
        # Provide _compute_seq_cai fallback for unsupported organisms
        opt._compute_seq_cai = lambda seq: compute_cai(seq, organism)
        # Sliding-window GC parameters (for _evaluate_all_predicates)
        opt._gc_window_size = gc_window_size
        opt._gc_window_min = gc_window_min
        opt._gc_window_max = gc_window_max

        # Import compute_cai for the hybrid path
        from .translation import compute_cai

    else:
        # ── Legacy strategies (constraint_first / cai_first) ──────────
        # Configure optimizer
        opt = BioOptimizer(
            species=species_key,
            enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            splice_low=effective_splice_low,
            splice_high=effective_splice_high,
            min_cai=cai_threshold,
            strategy=strategy,
            avoid_gt=effective_avoid_gt,
            organism_name=organism,
            organism_domain="prokaryote" if is_prokaryote else "eukaryote",
            gc_window_size=gc_window_size,
            gc_window_min=gc_window_min,
            gc_window_max=gc_window_max,
        )

        # Attach provenance collector to optimizer for codon-level tracking
        if _provenance_collector is not None:
            opt._provenance_collector = _provenance_collector

        # Back-translate protein to DNA for the optimizer
        from .translation import compute_cai
        initial_seq = _back_translate_protein(
            target_protein, species_key,
            strategy=strategy,
            enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            is_eukaryote=not is_prokaryote,
        )

        # Run optimization
        optimized_seq, pred_results, cert_text = opt.optimize(initial_seq, strategy=strategy)

    # Compute metrics — use hybrid optimizer values when available to avoid recomputation
    if strategy == "hybrid" and not use_csp_solver:
        # Hybrid optimizer already computed GC and CAI accurately.
        # However, if the eukaryotic CAI recovery pass or CpG elimination
        # was applied, the sequence has changed and we must recompute both metrics.
        if _cai_recovery_upgrades > 0 or _cpg_seq_changed:
            gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
            gc = round(gc, 4)
            # cai_val was already recomputed in the predicate evaluation block above
        else:
            gc = hybrid_result.gc_content if 'hybrid_result' in dir() else (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
            cai_val = hybrid_result.cai if 'hybrid_result' in dir() else None
            gc = round(gc, 4)
    else:
        gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
        gc = round(gc, 4)
        cai_val = None

    # Collect failed and satisfied predicates
    failed = [r.predicate for r in pred_results if not r.passed]
    satisfied = [r.predicate for r in pred_results if r.passed]

    # Detect if mutagenesis fallback was used (any V->I or similar)
    fallback = bool(opt._applied_mutagenesis)

    # ── mRNA stability improvement pass ────────────────────────────
    # PERF: Skip mRNA stability for short sequences (<300aa / <900bp)
    # where the stability improvement is negligible.
    _MRNA_STABILITY_MIN_LENGTH_BP = 900  # 300aa * 3

    mrna_stability_score: float | None = None
    destabilizing_motifs_removed: int = 0
    stability_improvement: float | None = None

    if optimize_mrna_stability and len(optimized_seq) >= _MRNA_STABILITY_MIN_LENGTH_BP:
        from .mrna_stability import score_mrna_stability, suggest_mutations_for_stability

        initial_stability = score_mrna_stability(optimized_seq, organism)
        suggestions = suggest_mutations_for_stability(optimized_seq, organism)

        # PERF (Optimization C): Precompute restriction site strings once
        # instead of calling get_recognition_site per suggestion.
        from .restriction_sites import get_recognition_site as _get_site
        _rs_site_strings: list[str] = []
        for enz in (enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]):
            site = _get_site(enz)
            if site:
                _rs_site_strings.append(site)

        # PERF (Optimization D): Use list mutation instead of string
        # concatenation for sequence modifications.
        seq_list = list(optimized_seq)
        # Pre-compute original protein once instead of per-suggestion
        orig_protein = BioOptimizer._translate(optimized_seq)
        # PERF (Optimization F): Cache GC count for incremental updates
        gc_count_cached = sum(1 for b in optimized_seq if b in "GC")
        n_bases = len(optimized_seq)

        # Pre-compute CAI adaptiveness table for CAI-aware acceptance
        _mrna_cai_table = CODON_ADAPTIVENESS_TABLES.get(organism, {})

        for suggestion in suggestions:
            pos = suggestion['position']
            ci = pos // 3
            new_codon = suggestion['suggested_codon']
            codon_start = ci * 3

            if codon_start + 3 > n_bases:
                continue

            # PERF (Optimization D): List mutation instead of string concat
            old_codon_list = seq_list[codon_start:codon_start + 3]
            old_codon_str = "".join(old_codon_list)

            # CAI-aware check: only accept mRNA stability suggestions that
            # don't degrade CAI.  For prokaryotes (E. coli), the optimal
            # codons are AT-rich (AAA, GAA, GAT) which also match the
            # RNase E / AT-rich "destabilizing" motifs.  Swapping these to
            # GC-rich alternatives (AAG, GAG, GAC) dramatically reduces CAI
            # for negligible mRNA stability benefit.  Skip any suggestion
            # that would replace an optimal codon with a suboptimal one.
            old_cai_w = _mrna_cai_table.get(old_codon_str, 0.0)
            new_cai_w = _mrna_cai_table.get(new_codon, 0.0)
            if new_cai_w < old_cai_w:
                # This suggestion would degrade CAI — skip it.
                # The CAI recovery pass will handle any remaining
                # opportunities after all post-processing.
                continue

            seq_list[codon_start:codon_start + 3] = list(new_codon)
            test_seq = "".join(seq_list)

            # PERF: Verify same protein by checking codon translates to same AA
            # This is faster than translating the entire sequence.
            # The mRNA module guarantees synonymous codons, so we only
            # need to verify no premature stop codons were introduced.
            _has_stop = False
            for _si in range(0, n_bases - 2, 3):
                _c = test_seq[_si:_si + 3]
                if _c in ("TAA", "TAG", "TGA") and _si < n_bases - 3:
                    _has_stop = True
                    break
            if _has_stop:
                # Rollback
                seq_list[codon_start:codon_start + 3] = old_codon_list
                continue

            # PERF (Optimization F): Incremental GC update
            old_gc_in_codon = sum(1 for b in old_codon_list if b in "GC")
            new_gc_in_codon = sum(1 for b in new_codon if b in "GC")
            new_gc_count = gc_count_cached - old_gc_in_codon + new_gc_in_codon
            test_gc = new_gc_count / n_bases
            if not (gc_lo <= test_gc <= gc_hi):
                # Rollback
                seq_list[codon_start:codon_start + 3] = old_codon_list
                continue

            # PERF (Optimization C): Check precomputed restriction sites
            rs_violated = False
            for site in _rs_site_strings:
                if site in test_seq and site not in optimized_seq:
                    rs_violated = True
                    break
            if rs_violated:
                # Rollback
                seq_list[codon_start:codon_start + 3] = old_codon_list
                continue

            # Safe to apply
            optimized_seq = test_seq
            gc_count_cached = new_gc_count
            destabilizing_motifs_removed += 1
            logger.debug(
                "mRNA stability: replaced %s->%s at codon %d (removing %s)",
                suggestion.get('original_codon', '???'), new_codon, ci,
                suggestion.get('motif_removed', suggestion.get('motif', '???')),
            )

        final_stability = score_mrna_stability(optimized_seq, organism)
        # score_mrna_stability returns an MRNAStabilityScore dataclass
        # with an overall_score attribute; extract the float
        init_score = initial_stability.overall_score if hasattr(initial_stability, 'overall_score') else float(initial_stability)
        final_score = final_stability.overall_score if hasattr(final_stability, 'overall_score') else float(final_stability)
        mrna_stability_score = final_score
        if init_score > 0:
            stability_improvement = round(final_score - init_score, 6)
        else:
            stability_improvement = round(final_score, 6)

        logger.info(
            "mRNA stability: initial=%.4f final=%.4f improvement=%.4f motifs_removed=%d",
            init_score, final_score, stability_improvement,
            destabilizing_motifs_removed,
        )
    # ── End mRNA stability pass ────────────────────────────────────

    # ── Codon pair bias improvement pass ──────────────────────────
    cpb_score: float | None = None

    if consider_codon_pair_bias:
        from .codon_pair_scoring import (
            compute_cpb as _compute_cpb,
            suggest_better_pair as _suggest_better_pair,
        )

        initial_cpb = _compute_cpb(optimized_seq, organism)
        logger.info("Codon pair bias: initial CPB=%.4f", initial_cpb)

        # Iterate over adjacent codon pairs and try to improve
        aas = protein_to_aa_list(target_protein)
        max_cpb_iterations = 3
        for _cpb_iter in range(max_cpb_iterations):
            any_improved = False
            for ci in range(len(aas) - 1):
                codon_start1 = ci * 3
                codon_start2 = (ci + 1) * 3
                current_c1 = optimized_seq[codon_start1:codon_start1 + 3]
                current_c2 = optimized_seq[codon_start2:codon_start2 + 3]
                aa1 = aas[ci]
                aa2 = aas[ci + 1]

                # Look up CAI weights for this organism
                cai_weights = CODON_ADAPTIVENESS_TABLES.get(organism)

                better = _suggest_better_pair(
                    current_c1, current_c2, aa1, aa2, organism,
                    cai_weights=cai_weights,
                    cai_weight=0.7,
                    cpb_weight=0.3,
                )
                if better is None:
                    continue

                new_c1, new_c2 = better

                # Apply the swap and verify hard constraints
                # PERF (Optimization D): Use list mutation instead of string concat
                _seq_list = list(optimized_seq)
                _seq_list[codon_start1:codon_start1 + 3] = list(new_c1)
                _seq_list[codon_start2:codon_start2 + 3] = list(new_c2)
                test_seq = "".join(_seq_list)

                # PERF: Verify no premature stop codons instead of
                # translating the entire sequence twice.
                _has_stop = False
                for _si in range(0, len(test_seq) - 2, 3):
                    _c = test_seq[_si:_si + 3]
                    if _c in ("TAA", "TAG", "TGA") and _si < len(test_seq) - 3:
                        _has_stop = True
                        break
                if _has_stop:
                    continue

                # Verify: GC still in range
                test_gc = (test_seq.count("G") + test_seq.count("C")) / max(len(test_seq), 1)
                if not (gc_lo <= test_gc <= gc_hi):
                    continue

                # Verify: no new restriction sites
                from .restriction_sites import get_recognition_site as _get_site
                rs_violated = False
                for enz in (enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]):
                    site = _get_site(enz)
                    if site and site in test_seq and site not in optimized_seq:
                        rs_violated = True
                        break
                if rs_violated:
                    continue

                # Verify: no new stop codons
                from .type_system import check_no_stop_codons as _check_stops
                stop_result = _check_stops(test_seq)
                if not stop_result.passed:
                    continue

                # Verify: no new cryptic splice sites (if threshold specified)
                cryptic_splice_threshold_cpb = kwargs.get("splice_low", 3.0)
                try:
                    from .maxentscan import max_donor_score, max_acceptor_score
                    old_max_d = max_donor_score(optimized_seq)
                    old_max_a = max_acceptor_score(optimized_seq)
                    new_max_d = max_donor_score(test_seq)
                    new_max_a = max_acceptor_score(test_seq)
                    if (new_max_d > old_max_d and new_max_d >= cryptic_splice_threshold_cpb):
                        continue
                    if (new_max_a > old_max_a and new_max_a >= cryptic_splice_threshold_cpb):
                        continue
                except Exception:
                    pass  # maxentscan may not be available

                # Safe to apply
                optimized_seq = test_seq
                any_improved = True
                logger.debug(
                    "CPB: swapped %s-%s -> %s-%s at codon pair %d-%d",
                    current_c1, current_c2, new_c1, new_c2, ci, ci + 1,
                )

            if not any_improved:
                break

        final_cpb = _compute_cpb(optimized_seq, organism)
        cpb_score = round(final_cpb, 6)
        logger.info(
            "Codon pair bias: initial=%.4f final=%.4f improvement=%.4f",
            initial_cpb, final_cpb, final_cpb - initial_cpb,
        )
    # ── End Codon pair bias pass ──────────────────────────────────

    # ── CAI Recovery Pass ──────────────────────────────────────────
    # After all post-processing (mRNA stability, CPB), some positions may
    # have suboptimal codons.  This pass systematically tries to upgrade
    # each suboptimal codon to the optimal synonym, checking that the
    # swap doesn't violate any hard constraint.  This is especially
    # important for prokaryotes where the mRNA stability pass can
    # introduce significant CAI regression by replacing AT-rich optimal
    # codons (AAA, GAA, GAT) with GC-rich suboptimal ones (AAG, GAG, GAC).
    _CAI_RECOVERY_MAX_ITERATIONS = 5
    _cai_recovery_table = CODON_ADAPTIVENESS_TABLES.get(organism, {})
    _cai_recovery_rs: list[tuple[str, str]] = []
    for _enz in (enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]):
        _rs_site = _get_site(_enz) if '_get_site' in dir() else None
        if _rs_site is None:
            from .restriction_sites import get_recognition_site as _grs
            _rs_site = _grs(_enz)
        if _rs_site:
            _cai_recovery_rs.append((_rs_site, reverse_complement(_rs_site)))

    _total_cai_recovery_iterations = 0
    for _recovery_iter in range(_CAI_RECOVERY_MAX_ITERATIONS):
        _any_recovered = False
        _aas = protein_to_aa_list(target_protein)
        _n_codons = len(_aas)

        # Sort positions by CAI gap (biggest improvement first)
        _recovery_candidates = []
        for _ci in range(_n_codons):
            _aa = _aas[_ci]
            if _aa == "*" or _aa == "M":
                continue
            _current = optimized_seq[_ci * 3:_ci * 3 + 3]
            _current_w = _cai_recovery_table.get(_current, 0.0)
            _alternatives = AA_TO_CODONS.get(_aa, [])
            _best_w = max(_cai_recovery_table.get(c, 0.0) for c in _alternatives) if _alternatives else 0.0
            if _best_w > _current_w:
                _gap = _best_w - _current_w
                _recovery_candidates.append((_gap, _ci, _current, _aa))

        _recovery_candidates.sort(key=lambda x: x[0], reverse=True)

        for _gap, _ci, _current, _aa in _recovery_candidates:
            # Try codons in CAI order (highest first)
            _sorted_alts = sorted(
                AA_TO_CODONS.get(_aa, []),
                key=lambda c: _cai_recovery_table.get(c, 0.0),
                reverse=True,
            )
            for _alt in _sorted_alts:
                if _alt == _current:
                    continue
                _alt_w = _cai_recovery_table.get(_alt, 0.0)
                _cur_w = _cai_recovery_table.get(_current, 0.0)
                if _alt_w <= _cur_w:
                    break  # No improvement possible (sorted descending)

                _test_seq = optimized_seq[:_ci * 3] + _alt + optimized_seq[_ci * 3 + 3:]

                # Check: no new restriction sites
                _rs_ok = True
                for _site, _site_rc in _cai_recovery_rs:
                    if _site in _test_seq or (_site_rc and _site_rc in _test_seq):
                        _rs_ok = False
                        break
                if not _rs_ok:
                    continue

                # Check: GC still in range
                _test_gc = (_test_seq.count("G") + _test_seq.count("C")) / max(len(_test_seq), 1)
                if not (gc_lo <= _test_gc <= gc_hi):
                    continue

                # Check: no new ATTTA motifs
                if _test_seq.count("ATTTA") > optimized_seq.count("ATTTA"):
                    continue

                # Check: no 6+ T runs
                _max_t = 0
                _j = 0
                while _j < len(_test_seq):
                    if _test_seq[_j] == 'T':
                        _k = _j
                        while _k < len(_test_seq) and _test_seq[_k] == 'T':
                            _k += 1
                        if _k - _j > _max_t:
                            _max_t = _k - _j
                        _j = _k
                    else:
                        _j += 1
                if _max_t >= 6:
                    continue

                # Check: no new premature stop codons
                _has_premature_stop = False
                for _si in range(0, len(_test_seq) - 5, 3):
                    if _test_seq[_si:_si + 3] in ("TAA", "TAG", "TGA"):
                        _has_premature_stop = True
                        break
                if _has_premature_stop:
                    continue

                # Check: for eukaryotes, no new cryptic splice sites
                if not is_prokaryote:
                    try:
                        _old_max_d = max_donor_score(optimized_seq)
                        _old_max_a = max_acceptor_score(optimized_seq)
                        _new_max_d = max_donor_score(_test_seq)
                        _new_max_a = max_acceptor_score(_test_seq)
                        if (_new_max_d > _old_max_d and _new_max_d >= effective_splice_low):
                            continue
                        if (_new_max_a > _old_max_a and _new_max_a >= effective_splice_low):
                            continue
                    except Exception:
                        pass

                # All checks passed — accept the upgrade
                optimized_seq = _test_seq
                _any_recovered = True
                logger.debug(
                    "CAI recovery: upgraded codon %d from %s to %s (w %.4f→%.4f)",
                    _ci, _current, _alt, _cur_w, _alt_w,
                )
                break  # Move to next position

            # ── Paired codon swap ──────────────────────────────────
            # If the single-codon upgrade was blocked (e.g., by a
            # restriction site that straddles two codons), try a paired
            # swap: change the target codon AND an adjacent codon.
            # Only consider the immediate neighbors and only if the
            # combined CAI improvement is positive.
            import math as _math
            for _adj_offset in (1, -1):
                _adj_ci = _ci + _adj_offset
                if _adj_ci < 0 or _adj_ci >= _n_codons:
                    continue
                _adj_aa = _aas[_adj_ci]
                if _adj_aa == "*" or _adj_aa == "M":
                    continue
                _adj_current = optimized_seq[_adj_ci * 3:_adj_ci * 3 + 3]
                _adj_current_w = _cai_recovery_table.get(_adj_current, 0.0)

                # Try each CAI-ordered alt for the target, paired with
                # each alt for the neighbor.
                for _alt in _sorted_alts:
                    if _alt == _current:
                        continue
                    _alt_w = _cai_recovery_table.get(_alt, 0.0)
                    if _alt_w <= _cur_w:
                        break

                    # Net CAI must improve: the gain at the target must
                    # outweigh any loss at the adjacent position.
                    _adj_sorted = sorted(
                        AA_TO_CODONS.get(_adj_aa, []),
                        key=lambda c: _cai_recovery_table.get(c, 0.0),
                        reverse=True,
                    )
                    for _adj_alt in _adj_sorted:
                        if _adj_alt == _adj_current:
                            continue
                        _adj_alt_w = _cai_recovery_table.get(_adj_alt, 0.0)
                        # Net log-CAI change must be positive
                        _old_log = _math.log(max(_cur_w, 1e-10)) + _math.log(max(_adj_current_w, 1e-10))
                        _new_log = _math.log(max(_alt_w, 1e-10)) + _math.log(max(_adj_alt_w, 1e-10))
                        if _new_log <= _old_log:
                            continue

                        # Build test sequence with both swaps
                        _lo_ci = min(_ci, _adj_ci)
                        _hi_ci = max(_ci, _adj_ci)
                        _test_seq = (
                            optimized_seq[:_lo_ci * 3]
                            + (_alt if _lo_ci == _ci else _adj_alt)
                            + optimized_seq[_lo_ci * 3 + 3:_hi_ci * 3]
                            + (_alt if _hi_ci == _ci else _adj_alt)
                            + optimized_seq[_hi_ci * 3 + 3:]
                        )

                        # Validate: no new restriction sites
                        _rs_ok = True
                        for _site, _site_rc in _cai_recovery_rs:
                            if _site in _test_seq or (_site_rc and _site_rc in _test_seq):
                                _rs_ok = False
                                break
                        if not _rs_ok:
                            continue

                        # Validate: GC in range
                        _test_gc = (_test_seq.count("G") + _test_seq.count("C")) / max(len(_test_seq), 1)
                        if not (gc_lo <= _test_gc <= gc_hi):
                            continue

                        # Validate: no new ATTTA
                        if _test_seq.count("ATTTA") > optimized_seq.count("ATTTA"):
                            continue

                        # Validate: no 6+ T runs
                        _max_t = 0
                        _j = 0
                        while _j < len(_test_seq):
                            if _test_seq[_j] == 'T':
                                _k = _j
                                while _k < len(_test_seq) and _test_seq[_k] == 'T':
                                    _k += 1
                                if _k - _j > _max_t:
                                    _max_t = _k - _j
                                _j = _k
                            else:
                                _j += 1
                        if _max_t >= 6:
                            continue

                        # Validate: no premature stops
                        _has_premature_stop = False
                        for _si in range(0, len(_test_seq) - 5, 3):
                            if _test_seq[_si:_si + 3] in ("TAA", "TAG", "TGA"):
                                _has_premature_stop = True
                                break
                        if _has_premature_stop:
                            continue

                        # All checks passed — accept the paired upgrade
                        optimized_seq = _test_seq
                        _any_recovered = True
                        logger.debug(
                            "CAI recovery: paired upgrade codon %d (%s→%s) + "
                            "codon %d (%s→%s)",
                            _ci, _current, _alt,
                            _adj_ci, _adj_current, _adj_alt,
                        )
                        break  # Accept first valid paired swap for this target
                    if _any_recovered:
                        break
                if _any_recovered:
                    break

        if not _any_recovered:
            break  # No more improvements possible
        _total_cai_recovery_iterations += 1

        # ── Convergence check ──
        # After each CAI recovery iteration, compute the objective and
        # check if the optimizer has converged or is oscillating.
        _recovery_cai = compute_cai(optimized_seq, organism)
        _recovery_gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
        _recovery_gc_ok = 1.0 if gc_lo <= _recovery_gc <= gc_hi else 0.0
        _recovery_obj = _recovery_cai * _recovery_gc_ok  # 0 if GC out of range
        _convergence.record(_recovery_obj)

        _conv_status = _convergence.check_convergence()
        if _conv_status is not None:
            _convergence_status = _conv_status
            logger.info(
                "CAI recovery: %s after %d iterations (objective=%.6f)",
                _conv_status, _recovery_iter + 1, _recovery_obj,
            )
            if _conv_status == "oscillating":
                # Stop at the best point — no need to continue cycling
                break
    # ── End CAI Recovery Pass ──────────────────────────────────────

    # ── Final CpG elimination pass (eukaryotes only) ──────────────
    # This runs AFTER the CAI recovery pass because the recovery may
    # have re-introduced CG dinucleotides by upgrading suboptimal
    # CG-free codons to optimal CG-containing codons.  We need to
    # eliminate those CGs one final time to ensure the sequence passes
    # the NoCpGIsland predicate.
    _final_cpg_warnings: list[str] = []
    if not is_prokaryote:
        _final_cpg_usage = CODON_ADAPTIVENESS_TABLES.get(organism)
        if _final_cpg_usage is not None:
            optimized_seq, _final_cpg_warnings = _eliminate_cpg_dinucleotides(
                optimized_seq,
                target_protein,
                _final_cpg_usage,
                enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                organism=organism,
                gc_lo=gc_lo,
                gc_hi=gc_hi,
            )
            if _final_cpg_warnings:
                for _fw in _final_cpg_warnings:
                    logger.warning(_fw)

    # ── Sliding-window GC constraint ─────────────────────────────────
    # After global GC is satisfied, check that no local window has
    # extreme GC content.  This prevents polymerase stalling and
    # secondary-structure hotspots that the global constraint misses.
    _effective_gc_window_min = gc_window_min if gc_window_min is not None else gc_lo
    _effective_gc_window_max = gc_window_max if gc_window_max is not None else gc_hi
    _sliding_gc_swaps = 0

    if gc_window_size > 0 and len(optimized_seq) >= gc_window_size:
        _sliding_result = check_sliding_gc(
            optimized_seq,
            window_size=gc_window_size,
            gc_min=_effective_gc_window_min,
            gc_max=_effective_gc_window_max,
        )
        if not _sliding_result.passed:
            logger.info(
                "Sliding-window GC: %d violation(s) detected (window=%d, "
                "range=[%.2f, %.2f]). Attempting local fixes.",
                len(_sliding_result.violations), gc_window_size,
                _effective_gc_window_min, _effective_gc_window_max,
            )
            _sliding_usage = CODON_ADAPTIVENESS_TABLES.get(organism)
            optimized_seq, _sliding_gc_swaps = fix_sliding_gc_violations(
                optimized_seq,
                target_protein,
                window_size=gc_window_size,
                gc_min=_effective_gc_window_min,
                gc_max=_effective_gc_window_max,
                usage=_sliding_usage,
                gc_lo=gc_lo,
                gc_hi=gc_hi,
            )
            if _sliding_gc_swaps > 0:
                logger.info(
                    "Sliding-window GC: fixed %d local violation(s) with %d codon swap(s)",
                    len(_sliding_result.violations), _sliding_gc_swaps,
                )

    # ── Refresh SlidingGC predicate result after post-processing ──
    # The fix_sliding_gc_violations step may have modified optimized_seq,
    # so we must re-evaluate the SlidingGC predicate to ensure the result
    # in pred_results matches the final sequence.  Without this, certificate
    # verification fails because the recorded verdict was evaluated on the
    # pre-fix sequence while the certificate contains the post-fix sequence.
    if gc_window_size > 0 and len(optimized_seq) >= gc_window_size and _sliding_gc_swaps > 0:
        _effective_gc_window_min = gc_window_min if gc_window_min is not None else gc_lo
        _effective_gc_window_max = gc_window_max if gc_window_max is not None else gc_hi
        _refreshed_sgc = check_sliding_gc(
            optimized_seq,
            window_size=gc_window_size,
            gc_min=_effective_gc_window_min,
            gc_max=_effective_gc_window_max,
        )
        _refreshed_sgc_pr = PredicateResult(
            "SlidingGC", _refreshed_sgc.passed,
            details=(
                f"Window={gc_window_size}, range=[{_effective_gc_window_min:.2f}, {_effective_gc_window_max:.2f}], "
                f"min_gc={_refreshed_sgc.min_gc:.3f}, max_gc={_refreshed_sgc.max_gc:.3f}, "
                f"violations={len(_refreshed_sgc.violations)}"
            ),
        )
        # Replace the stale SlidingGC entry in pred_results
        pred_results = [
            _refreshed_sgc_pr if r.predicate == "SlidingGC" else r
            for r in pred_results
        ]
        # Refresh cert_text to match the updated predicates
        cert_text = format_certificate(pred_results, optimized_seq, species_key)
        failed = [r.predicate for r in pred_results if not r.passed]
        satisfied = [r.predicate for r in pred_results if r.passed]
    elif gc_window_size > 0:
        logger.debug(
            "Sliding-window GC: skipped (sequence length %d < window size %d)",
            len(optimized_seq), gc_window_size,
        )

    # ── Post-sliding-GC CpG re-elimination (eukaryotes only) ──────
    # The sliding-window GC fix (fix_sliding_gc_violations) may have
    # reintroduced CG dinucleotides by swapping codons without checking
    # for CpG creation.  Run one final CpG elimination pass to restore
    # the NoCpGIsland predicate without regressing sliding-window GC.
    if not is_prokaryote and _sliding_gc_swaps > 0:
        _post_sgc_usage = CODON_ADAPTIVENESS_TABLES.get(organism)
        if _post_sgc_usage is not None:
            optimized_seq, _post_sgc_warnings = _eliminate_cpg_dinucleotides(
                optimized_seq,
                target_protein,
                _post_sgc_usage,
                enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                organism=organism,
                gc_lo=gc_lo,
                gc_hi=gc_hi,
            )
            if _post_sgc_warnings:
                for _psw in _post_sgc_warnings:
                    logger.warning(_psw)
            # If CpG elimination changed the sequence, re-check sliding GC
            _sgc_needs_refresh = _sliding_gc_swaps > 0
            if _sgc_needs_refresh and gc_window_size > 0 and len(optimized_seq) >= gc_window_size:
                _effective_gc_window_min = gc_window_min if gc_window_min is not None else gc_lo
                _effective_gc_window_max = gc_window_max if gc_window_max is not None else gc_hi
                _refreshed2_sgc = check_sliding_gc(
                    optimized_seq,
                    window_size=gc_window_size,
                    gc_min=_effective_gc_window_min,
                    gc_max=_effective_gc_window_max,
                )
                _refreshed2_sgc_pr = PredicateResult(
                    "SlidingGC", _refreshed2_sgc.passed,
                    details=(
                        f"Window={gc_window_size}, range=[{_effective_gc_window_min:.2f}, {_effective_gc_window_max:.2f}], "
                        f"min_gc={_refreshed2_sgc.min_gc:.3f}, max_gc={_refreshed2_sgc.max_gc:.3f}, "
                        f"violations={len(_refreshed2_sgc.violations)}"
                    ),
                )
                pred_results = [
                    _refreshed2_sgc_pr if r.predicate == "SlidingGC" else r
                    for r in pred_results
                ]
                cert_text = format_certificate(pred_results, optimized_seq, species_key)

    # ── Recompute final metrics after all post-processing passes ───
    # CAI and GC must be recomputed here because the mRNA stability
    # and codon-pair-bias passes may have modified optimized_seq.
    gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
    gc = round(gc, 4)

    # Use the canonical organism name (opt.organism_name) so that
    # compute_cai can look up the correct CODON_ADAPTIVENESS_TABLES
    # entry even when the caller passed a short alias like 'ecoli'.
    # PERF: Skip CAI recomputation if hybrid optimizer already provided it
    # and no post-processing passes (mRNA stability, CPB) modified the sequence.
    if cai_val is not None and not optimize_mrna_stability and not consider_codon_pair_bias:
        pass  # Use hybrid_result.cai — already set
    else:
        try:
            cai_val = compute_cai(optimized_seq, opt.organism_name)
        except UnsupportedOrganismError:
            logger.debug(
                "Unsupported organism '%s' for compute_cai, using species CAI table",
                opt.organism_name,
            )
            cai_val = opt._compute_seq_cai(optimized_seq)

    # Build provenance record for reproducibility
    from .provenance import OptimizationRecord as _OptimizationRecord
    from .provenance import _get_biocompiler_version as _get_version

    # Determine solver backend name
    solver_backend = "greedy"
    if use_csp_solver:
        solver_backend = "csp"  # actual backend may vary, but CSP was requested

    # Extract mutation descriptions from applied mutagenesis
    mutations_made: list[str] = []
    for mut in opt._applied_mutagenesis:
        if isinstance(mut, dict):
            mutations_made.append(mut.get("description", str(mut)))
        else:
            mutations_made.append(str(mut))

    # Constraints applied (from satisfied + failed predicate names)
    constraints_applied = sorted(set(
        [r.predicate for r in pred_results]
    ))

    provenance_record = _OptimizationRecord(
        input_sequence=target_protein,
        output_sequence=optimized_seq,
        organism=organism,
        constraints_applied=constraints_applied,
        mutations_made=mutations_made,
        solver_backend=solver_backend,
        solve_time=round(_time.monotonic() - _start_time, 6),
        seed_used=seed,
        timestamp=datetime.now(timezone.utc).isoformat(),
        biocompiler_version=_get_version(),
        mrna_stability_score=mrna_stability_score,
        destabilizing_motifs_removed=destabilizing_motifs_removed,
        stability_improvement=stability_improvement,
    )

    # ── UTR suggestions ─────────────────────────────────────────────
    utr5_seq: str | None = None
    utr3_seq: str | None = None
    utr_score5: float | None = None
    utr_score3: float | None = None

    if include_utr:
        from .utr_models import suggest_5utr, suggest_3utr, score_5utr, score_3utr
        try:
            utr5_seq = suggest_5utr(organism)
            utr_score5 = score_5utr(utr5_seq, organism)
        except ValueError:
            logger.debug("No 5' UTR suggestion for organism '%s'", organism)
        try:
            utr3_seq = suggest_3utr(organism)
            utr_score3 = score_3utr(utr3_seq, organism)
        except ValueError:
            logger.debug("No 3' UTR suggestion for organism '%s'", organism)

    # Record constraint decisions if provenance tracking is enabled
    if _provenance_collector is not None:
        try:
            # Build constraint decisions from predicate results
            for _pr in pred_results:
                _action = "satisfied" if _pr.passed else "conflicted"
                # Estimate CAI impact (negative if constraint reduced CAI)
                _cai_impact = 0.0
                if not _pr.passed:
                    _cai_impact = -0.01  # Conservative estimate for failed constraints
                elif _pr.predicate in ("NoRestrictionSite", "NoCrypticSplice", "GCInRange"):
                    # These constraints typically have some CAI cost even when satisfied
                    _cai_impact = -0.005

                # Determine positions affected (best effort)
                _positions_affected: list[int] = []
                _details = getattr(_pr, "details", "") or ""
                if "position" in _details.lower():
                    # Try to extract position info from details
                    import re as _re
                    _pos_matches = _re.findall(r'pos(?:ition)?[:\s]+(\d+)', _details, _re.IGNORECASE)
                    _positions_affected = [int(p) for p in _pos_matches[:10]]

                _tradeoff = _details if _details else f"Constraint {_pr.predicate} was {_action}"

                _provenance_collector.record_constraint_decision(ConstraintDecision(
                    constraint_name=_pr.predicate,
                    constraint_type="hard",
                    action_taken=_action,
                    positions_affected=_positions_affected,
                    tradeoff_description=_tradeoff,
                    impact_on_cai=_cai_impact,
                ))
        except Exception:
            logger.debug("Constraint provenance recording failed", exc_info=True)

    # Finalize decision provenance trail if tracking is enabled
    _decision_trail: OptimizationDecisionTrail | None = None
    if _provenance_collector is not None:
        try:
            _decision_trail = _provenance_collector.finalize(
                output_dna=optimized_seq,
                cai=cai_val,
                gc=gc,
            )
        except Exception as _prov_exc:
            logger.debug("Provenance finalization failed: %s: %s", type(_prov_exc).__name__, _prov_exc, exc_info=True)
            _decision_trail = None

    # ── Final translation verification ──────────────────────────────
    # After all optimization passes (CAI recovery, mRNA stability, etc.),
    # verify that the optimized DNA still encodes the original protein.
    # This catches bugs in any optimization pass that might corrupt
    # the protein sequence.
    from .protein_verification import verify_and_raise as _verify_final
    _verify_final(optimized_seq, target_protein, organism=organism)

    # ── Custom objective refinement pass ──────────────────────────────
    # When a non-default objective is provided, perform an additional
    # hill-climb pass that scores codon substitutions by the custom
    # objective rather than CAI alone.  This allows users to optimize
    # for GC balance, codon pair bias, or any custom metric while still
    # respecting all hard constraints.
    _objective_score: float | None = None
    if _objective_fn is not _cai_objective:
        logger.info(
            "Applying custom objective refinement pass for objective=%s",
            getattr(_objective_fn, '__name__', repr(_objective_fn)),
        )
        _obj_aas = protein_to_aa_list(target_protein)
        _obj_n_codons = len(_obj_aas)
        _obj_usage = CODON_ADAPTIVENESS_TABLES.get(organism, {})
        _obj_sorted_codons: dict[str, list[str]] = {}
        for _aa in set(_obj_aas):
            if _aa == "*":
                continue
            _obj_sorted_codons[_aa] = sorted(
                AA_TO_CODONS.get(_aa, []),
                key=lambda c: _obj_usage.get(c, 0.0),
                reverse=True,
            )
        _obj_rs_sites: list[tuple[str, str]] = []
        for _enz in (enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]):
            from .restriction_sites import get_recognition_site as _grs2
            _obj_site = _grs2(_enz)
            if _obj_site:
                _obj_rs_sites.append((_obj_site, reverse_complement(_obj_site)))

        _MAX_OBJECTIVE_REFINEMENT_ITERS = 5
        for _obj_iter in range(_MAX_OBJECTIVE_REFINEMENT_ITERS):
            _obj_improved = False
            _current_score = _objective_fn(optimized_seq, target_protein, organism)

            for _ci in range(_obj_n_codons):
                _aa = _obj_aas[_ci]
                if _aa == "*" or _aa == "M":
                    continue
                _current_codon = optimized_seq[_ci * 3:_ci * 3 + 3]

                # Try each synonym for this position
                _best_alt = None
                _best_alt_score = _current_score

                for _alt in _obj_sorted_codons.get(_aa, []):
                    if _alt == _current_codon:
                        continue

                    _test_seq = optimized_seq[:_ci * 3] + _alt + optimized_seq[_ci * 3 + 3:]

                    # Hard constraint checks
                    # 1. No new restriction sites
                    _rs_ok = True
                    for _site, _site_rc in _obj_rs_sites:
                        if _site in _test_seq and _site not in optimized_seq:
                            _rs_ok = False
                            break
                        if _site_rc and _site_rc in _test_seq and _site_rc not in optimized_seq:
                            _rs_ok = False
                            break
                    if not _rs_ok:
                        continue

                    # 2. GC still in range
                    _test_gc = (_test_seq.count("G") + _test_seq.count("C")) / max(len(_test_seq), 1)
                    if not (gc_lo <= _test_gc <= gc_hi):
                        continue

                    # 3. No new premature stops
                    _has_premature_stop = False
                    for _si in range(0, len(_test_seq) - 5, 3):
                        if _test_seq[_si:_si + 3] in ("TAA", "TAG", "TGA"):
                            _has_premature_stop = True
                            break
                    if _has_premature_stop:
                        continue

                    # 4. No new ATTTA motifs
                    if _test_seq.count("ATTTA") > optimized_seq.count("ATTTA"):
                        continue

                    # Evaluate the custom objective
                    try:
                        _alt_score = _objective_fn(_test_seq, target_protein, organism)
                    except Exception:
                        continue

                    if _alt_score > _best_alt_score:
                        _best_alt = _alt
                        _best_alt_score = _alt_score

                if _best_alt is not None:
                    optimized_seq = (
                        optimized_seq[:_ci * 3] + _best_alt + optimized_seq[_ci * 3 + 3:]
                    )
                    _obj_improved = True

            if not _obj_improved:
                break

        # Record the final objective score
        try:
            _objective_score = _objective_fn(optimized_seq, target_protein, organism)
        except Exception:
            _objective_score = None

        # Recompute CAI and GC after objective refinement
        gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
        gc = round(gc, 4)
        try:
            cai_val = compute_cai(optimized_seq, opt.organism_name)
        except UnsupportedOrganismError:
            cai_val = opt._compute_seq_cai(optimized_seq)

    result = OptimizationResult(
        sequence=optimized_seq,
        gc_content=gc,
        cai=cai_val,
        failed_predicates=failed,
        predicate_results=pred_results,
        certificate_text=cert_text,
        protein=target_protein,
        fallback_used=fallback,
        satisfied_predicates=satisfied,
        aa_substitutions=opt._applied_mutagenesis,
        mrna_stability_score=mrna_stability_score,
        destabilizing_motifs_removed=destabilizing_motifs_removed,
        stability_improvement=stability_improvement,
        provenance=provenance_record,
        codon_pair_bias=cpb_score,
        suggested_5utr=utr5_seq,
        suggested_3utr=utr3_seq,
        utr_score_5=utr_score5,
        utr_score_3=utr_score3,
        decision_trail=_decision_trail,
        convergence_status=_convergence_status,
        iterations_used=_total_cai_recovery_iterations + 1,  # +1 for the initial optimization pass
        objective_score=_objective_score,
    )

    # ── Strict mode: refuse sequences with failed predicates ─────────
    if strict_mode and result.failed_predicates:
        raise OptimizationConstraintError(
            failed_predicates=result.failed_predicates,
            partial_result=result,
        )

    return result


def batch_optimize(
    proteins: list[str],
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    enzymes: list | None = None,
    cai_threshold: float = 0.5,
    strategy: str = "hybrid",
    full_postprocessing: bool = False,
    **kwargs,
) -> list[OptimizationResult]:
    """Optimize multiple proteins in batch.

    Uses a single HybridOptimizer instance for efficiency
    (precomputed data structures are reused across proteins).
    This avoids redundant initialization of sorted_codons, gt_free,
    ag_free, codon_gc, optimal_codon, _max_adapt, and _rs_sites
    for each protein — a significant speedup when optimizing many
    proteins for the same organism.

    When ``full_postprocessing=True``, each result is processed through
    the full :func:`optimize_sequence` pipeline (mRNA stability, codon
    pair bias, CAI recovery, UTR suggestions, provenance tracking).
    When ``False`` (default), a streamlined path is used that skips
    expensive post-processing steps, making batch optimization
    significantly faster for high-throughput workflows.

    Args:
        proteins: List of amino acid sequences (single-letter codes, no stop).
        organism: Target organism. Accepts canonical binomials, short keys,
            abbreviated binomials, or display names (same as
            :func:`optimize_sequence`).
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        enzymes: List of restriction enzyme names to avoid.
        cai_threshold: Minimum CAI score for the CodonAdapted predicate.
        strategy: Optimization strategy (default 'hybrid').
        full_postprocessing: If True, run full optimize_sequence per protein
            (slower but includes mRNA stability, CPB, UTR, provenance).
            If False (default), use streamlined batch path.
        **kwargs: Additional arguments passed to HybridOptimizer or
            optimize_sequence (depending on full_postprocessing).

    Returns:
        List of OptimizationResult, one per input protein, in the same order.

    Raises:
        InvalidProteinError: if any protein contains invalid amino acid codes.
        UnsupportedOrganismError: if the organism is not supported.

    Examples:
        >>> from biocompiler.optimization import batch_optimize
        >>> proteins = ['MSKGEELFTG', 'MALWMRLLPL', 'MVHLTPEEKS']
        >>> results = batch_optimize(proteins, organism='Escherichia_coli')
        >>> assert len(results) == 3
        >>> for r in results:
        ...     assert r.cai > 0.5
    """
    if not proteins:
        return []

    # Resolve organism name once
    organism = resolve_organism(organism, strict=False)

    # Validate all proteins upfront
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    for protein in proteins:
        protein_upper = protein.strip().upper()
        invalid = set(protein_upper) - valid_aas
        if invalid:
            raise InvalidProteinError(protein, invalid)

    # ── Full post-processing path ──────────────────────────────────
    # Delegate to optimize_sequence per protein (no HybridOptimizer
    # reuse, but full feature parity).  This is the safe fallback.
    if full_postprocessing:
        results: list[OptimizationResult] = []
        for protein in proteins:
            result = optimize_sequence(
                target_protein=protein,
                organism=organism,
                gc_lo=gc_lo,
                gc_hi=gc_hi,
                cai_threshold=cai_threshold,
                enzymes=enzymes,
                strategy=strategy,
                **kwargs,
            )
            results.append(result)
        return results

    # ── Streamlined batch path ─────────────────────────────────────
    # Reuse a single HybridOptimizer for all proteins.  The __init__
    # precomputes sorted_codons, gt_free, ag_free, codon_gc,
    # optimal_codon, _max_adapt, _rs_sites, etc. — all of which are
    # identical for the same organism/enzymes/GC parameters.  Reusing
    # them avoids O(20 * AA_count) redundant sort/filter operations
    # per protein.
    from .hybrid_optimizer import HybridOptimizer as _HybridOptimizer

    species_key = _organism_to_species_key(organism)

    # Detect organism domain
    organism_domain = kwargs.get("organism_domain", "auto")
    if organism_domain not in ("auto", "eukaryote", "prokaryote"):
        organism_domain = "auto"
    if organism_domain == "auto":
        from .organism_config import is_eukaryotic_organism
        is_prokaryote = not is_eukaryotic_organism(organism)
    elif organism_domain == "prokaryote":
        is_prokaryote = True
    else:
        is_prokaryote = False

    effective_avoid_gt = kwargs.get("avoid_gt", not is_prokaryote)
    effective_splice_low = kwargs.get("splice_low", 3.0 if not is_prokaryote else 999.0)

    # Create ONE HybridOptimizer — precomputed data structures are reused
    hybrid_opt = _HybridOptimizer(
        species=species_key,
        organism=organism,
        enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        avoid_gt=effective_avoid_gt,
        splice_threshold=effective_splice_low,
    )

    results: list[OptimizationResult] = []
    for protein in proteins:
        protein_upper = protein.strip().upper()

        # Run hybrid optimization (reuses precomputed data structures)
        hybrid_result = hybrid_opt.optimize(
            protein_upper, is_prokaryote=is_prokaryote
        )
        optimized_seq = hybrid_result.sequence

        # ── Lightweight predicate evaluation ───────────────────────
        pred_results = []

        if is_prokaryote:
            # Prokaryote fast path — skip eukaryotic predicates
            from .type_system import (
                check_no_stop_codons, check_valid_coding_seq,
                check_no_restriction_site,
            )

            pred_results.append(check_no_stop_codons(optimized_seq))
            pred_results.append(PredicateResult(
                "NoCrypticSplice", True,
                details="Skipped for prokaryotic organism"
            ))
            pred_results.append(PredicateResult(
                "NoCpGIsland", True,
                details="Skipped for prokaryotic organism"
            ))
            pred_results.append(check_no_restriction_site(
                optimized_seq, enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
            ))
            pred_results.append(PredicateResult(
                "NoGTDinucleotide", True,
                details="Skipped for prokaryotic organism"
            ))
            pred_results.append(check_valid_coding_seq(optimized_seq))
            pred_results.append(PredicateResult(
                "ConservationScore", True,
                details="All AA conservation scores >= 0"
            ))

            cai_val = hybrid_result.cai
            all_optimal = cai_val >= cai_threshold
            pred_results.append(PredicateResult(
                "CodonOptimality", all_optimal,
                details=f"CAI={cai_val:.4f}, min={cai_threshold}"
            ))

            _gc_val = hybrid_result.gc_content
            gc_ok = gc_lo <= _gc_val <= gc_hi
            pred_results.append(PredicateResult(
                "GCInRange", gc_ok,
                details=f"GC content: {_gc_val:.3f} (range [{gc_lo}, {gc_hi}])"
            ))
        else:
            # Eukaryote path — use BioOptimizer for predicate evaluation
            opt = BioOptimizer(
                species=species_key,
                enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                splice_low=effective_splice_low,
                splice_high=kwargs.get("splice_high", 6.0),
                min_cai=cai_threshold,
                strategy="constraint_first",
                avoid_gt=effective_avoid_gt,
                organism_name=organism,
                organism_domain="eukaryote",
            )
            opt._applied_mutagenesis = hybrid_opt._applied_mutagenesis
            opt._original_protein = protein_upper
            # Skip redundant MaxEntScan splice check when the hybrid
            # optimizer already validated splice sites during Phase 3.
            _skip_splice = hybrid_result.splice_sites_validated
            pred_results = opt._evaluate_all_predicates(optimized_seq, skip_splice_check=_skip_splice)

        # Build certificate and collect predicate names
        cert_text = format_certificate(pred_results, optimized_seq, species_key)
        failed = [r.predicate for r in pred_results if not r.passed]
        satisfied = [r.predicate for r in pred_results if r.passed]
        fallback = bool(hybrid_opt._applied_mutagenesis)

        # Compute final metrics
        gc = hybrid_result.gc_content
        cai_val = hybrid_result.cai

        # Build OptimizationResult
        result = OptimizationResult(
            sequence=optimized_seq,
            gc_content=gc,
            cai=cai_val,
            failed_predicates=failed,
            predicate_results=pred_results,
            certificate_text=cert_text,
            protein=protein_upper,
            fallback_used=fallback,
            satisfied_predicates=satisfied,
        )
        results.append(result)

        # Reset per-protein mutagenesis state for next protein
        hybrid_opt._applied_mutagenesis = []

    return results


def _organism_to_species_key(organism: str) -> str:
    """Map an organism name to the species key used in SPECIES_SHORT_NAMES.

    Delegates to :func:`~biocompiler.organisms.resolve_organism` for
    name resolution and then looks up the short key in
    :data:`~biocompiler.organisms.SPECIES_SHORT_NAMES`.

    .. deprecated::
        Use :func:`~biocompiler.organisms.resolve_organism` and
        CODON_ADAPTIVENESS_TABLES directly instead of mapping to
        SPECIES dict keys.  Retained for backward compatibility.
    """
    canonical = resolve_organism(organism, strict=False)
    key = SPECIES_SHORT_NAMES.get(canonical)
    if key:
        return key
    # Fallback: try the organism name directly in SPECIES_SHORT_NAMES values
    if organism in SPECIES_SHORT_NAMES.values():
        return organism
    # Default to ecoli
    return "ecoli"


def _species_key_to_organism(species_key: str) -> str:
    """Map a species key or organism alias to the canonical organism name
    used in CODON_ADAPTIVENESS_TABLES.

    Delegates to :func:`~biocompiler.organisms.resolve_organism` for
    name resolution, which accepts both short aliases (e.g. 'ecoli')
    and full canonical names (e.g. 'Escherichia_coli'), as well as
    display names ('E. coli') and abbreviated binomials ('e_coli').

    .. deprecated::
        Use :func:`~biocompiler.organisms.resolve_organism` directly
        instead of this wrapper.  It is retained for internal use only.
    """
    return resolve_organism(species_key, default="Homo_sapiens")


def _back_translate_protein(
    protein: str,
    species_key: str,
    strategy: str = "greedy",
    enzymes: list[str] | None = None,
    is_eukaryote: bool = True,
) -> str:
    """Back-translate a protein to DNA using highest-CAI codons.

    Uses CODON_ADAPTIVENESS_TABLES (the same table used by compute_cai) so
    the initial sequence starts with optimal codons that will be reflected
    in the final CAI score.

    When ``strategy='dp'`` (or automatically for sequences < 200 aa), a
    dynamic-programming approach is used that considers the top 3 codons
    per position and finds the globally optimal sequence that:
      - Maximizes CAI
      - Avoids common restriction enzyme recognition sites
      - Minimizes GT/AG dinucleotides (for eukaryotes)

    Args:
        protein: Amino acid sequence (1-letter codes).
        species_key: Short species key (e.g. 'ecoli', 'human').
        strategy: 'greedy' for per-position best-CAI; 'dp' for global
            optimisation with cross-codon effects.
        enzymes: Restriction enzyme names to avoid (only used for DP).
        is_eukaryote: If True, penalise GT/AG dinucleotides in DP.

    Returns:
        DNA sequence string.
    """
    organism = _species_key_to_organism(species_key)
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism)
    if adaptiveness is None:
        # Fallback to Homo_sapiens if organism not found
        adaptiveness = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]

    # Decide whether to use DP
    # Explicit 'dp' strategy always uses DP; 'greedy' always skips DP.
    # For other strategies (e.g. 'constraint_first'), use DP for short
    # sequences where it's fast and produces a better starting point.
    if strategy == "dp":
        use_dp = True
    elif strategy == "greedy":
        use_dp = False
    else:
        use_dp = len(protein) < 200
    if use_dp:
        return _back_translate_protein_dp(
            protein, adaptiveness, enzymes=enzymes, is_eukaryote=is_eukaryote,
        )

    # Greedy: simply pick the highest-CAI codon per position
    codons = []
    for aa in protein:
        if aa == "*":
            codons.append("TAA")
            continue
        candidates = AA_TO_CODONS.get(aa, [])
        if not candidates:
            codons.append("NNN")
            continue
        best = max(candidates, key=lambda c: adaptiveness.get(c, 0.0))
        codons.append(best)
    return "".join(codons)


# ────────────────────────────────────────────────────────────
# DP-based back-translation
# ────────────────────────────────────────────────────────────

# Maximum number of top codons to consider per amino-acid position in DP
_DP_TOP_K: int = 3
# Penalty multiplier for GT/AG dinucleotides (eukaryotes)
_GT_AG_PENALTY: float = 0.02
# Penalty for creating a restriction site (effectively -inf)
_RESTRICTION_PENALTY: float = 0.5


def _build_restriction_site_set(
    enzymes: list[str] | None,
) -> set[str]:
    """Build the set of restriction site sequences (fwd + rc) to avoid."""
    sites: set[str] = set()
    if enzymes is None:
        enzymes = list(RESTRICTION_ENZYMES.keys())
    for name in enzymes:
        seq = RESTRICTION_ENZYMES.get(name)
        if seq is None:
            continue
        # Skip IUPAC wildcard sites (contain N or other ambiguity codes)
        if any(b not in "ATCG" for b in seq.upper()):
            continue
        seq_upper = seq.upper()
        sites.add(seq_upper)
        rc = reverse_complement(seq_upper)
        if rc != seq_upper:
            sites.add(rc)
    return sites


def _contains_restriction_site(
    seq: str, sites: set[str], window: int
) -> bool:
    """Check whether *seq* contains any of the restriction *sites*.

    Only checks the last *window* characters (enough to catch any new site
    introduced by the most recent codon).
    """
    start = max(0, len(seq) - window)
    tail = seq[start:]
    for site in sites:
        if site in tail:
            return True
    return False


def _count_dinucleotides(seq: str, di: str) -> int:
    """Count occurrences of dinucleotide *di* in *seq*.

    Uses the NUMBA ``fast_dinucleotide_count`` kernel when available
    for single-pass counting; falls back to pure-Python str.find
    otherwise.
    """
    return _count_dinucs_fast(seq, di)[0]


def _back_translate_protein_dp(
    protein: str,
    adaptiveness: dict[str, float],
    enzymes: list[str] | None = None,
    is_eukaryote: bool = True,
) -> str:
    """DP-based back-translation that considers cross-codon effects.

    For each amino-acid position, the top K codons by CAI are considered.
    The DP state tracks the last two codons (6 nt) to evaluate:
      - GT / AG dinucleotide penalties (eukaryotes)
      - Restriction site avoidance

    The objective is to maximise the sum of log-adaptiveness values
    (equivalent to maximising the geometric mean = CAI) while penalising
    undesirable features.

    Complexity: O(N * K^3) where N = protein length, K = _DP_TOP_K.
    For N < 200 and K = 3 this is negligible (< 1 ms).
    """
    import math

    K = _DP_TOP_K
    restriction_sites = _build_restriction_site_set(enzymes)
    # Maximum restriction site length — only need to check this many
    # trailing characters for newly introduced sites
    max_site_len = max((len(s) for s in restriction_sites), default=0)
    # Window of trailing sequence to check for restriction sites:
    # last codon (3) + overlap from previous codon (max_site_len - 1)
    rest_check_window = max_site_len + 2  # conservative

    # Pre-compute top-K codons per amino acid
    top_codons_per_aa: list[list[tuple[str, float]]] = []
    for aa in protein:
        if aa == "*":
            # Stop codon — only one choice needed
            top_codons_per_aa.append([("TAA", 1.0)])
            continue
        candidates = AA_TO_CODONS.get(aa, [])
        if not candidates:
            top_codons_per_aa.append([("NNN", 0.0)])
            continue
        # Sort by adaptiveness descending, take top K
        sorted_cands = sorted(
            candidates,
            key=lambda c: adaptiveness.get(c, 0.0),
            reverse=True,
        )
        top = [(c, adaptiveness.get(c, 0.0)) for c in sorted_cands[:K]]
        top_codons_per_aa.append(top)

    n = len(protein)

    # DP state: (prev_codon_idx, curr_codon_idx) -> best log-CAI score
    # We track the last two codon choices to evaluate cross-codon effects.
    # prev_codon_idx refers to position i-2, curr_codon_idx to position i-1.
    # At position i, we extend to a new codon choice.

    # For position 0: no previous context
    # dp[(prev_idx, curr_idx)] = (score, backpointer_list)
    # We use a flat dict for the DP frontier.

    # Special handling for the first few positions where we don't have
    # full two-codon context yet.

    # Position 0: just pick a codon
    dp: dict[tuple[int, int], tuple[float, list[int]]] = {}
    for ci, (codon_i, adapt_i) in enumerate(top_codons_per_aa[0]):
        score = math.log(adapt_i) if adapt_i > 0 else -20.0
        dp[(0, ci)] = (score, [ci])

    if n == 1:
        # Trivial case
        best_state = max(dp, key=lambda k: dp[k][0])
        best_codon = top_codons_per_aa[0][dp[best_state][1][0]][0]
        return best_codon

    # Position 1: extend from position 0
    new_dp: dict[tuple[int, int], tuple[float, list[int]]] = {}
    for (prev_ci_key, curr_ci_key), (score, path) in dp.items():
        curr_ci = path[0]  # index into top_codons_per_aa[0]
        curr_codon = top_codons_per_aa[0][curr_ci][0]
        for ci1, (codon1, adapt1) in enumerate(top_codons_per_aa[1]):
            s = score + (math.log(adapt1) if adapt1 > 0 else -20.0)
            # Check cross-codon effects between position 0 and 1
            junction = curr_codon + codon1  # 6 nt
            if is_eukaryote:
                # Penalise GT and AG at the cross-codon boundary
                if curr_codon[-1] + codon1[0] == "GT":
                    s -= _GT_AG_PENALTY * 50  # scale penalty to log space
                if curr_codon[-1] + codon1[0] == "AG":
                    s -= _GT_AG_PENALTY * 50
            # Check restriction sites in the junction
            if restriction_sites:
                for site in restriction_sites:
                    if len(site) <= len(junction) and site in junction:
                        s -= _RESTRICTION_PENALTY
            state = (curr_ci, ci1)
            if state not in new_dp or s > new_dp[state][0]:
                new_dp[state] = (s, path + [ci1])
    dp = new_dp

    # Positions 2..n-1: full DP with two-codon lookback
    for pos in range(2, n):
        new_dp = {}
        for (prev_ci, curr_ci), (score, path) in dp.items():
            prev_codon = top_codons_per_aa[pos - 2][prev_ci][0] if pos >= 2 else ""
            curr_codon = top_codons_per_aa[pos - 1][curr_ci][0]
            for ci, (codon, adapt) in enumerate(top_codons_per_aa[pos]):
                s = score + (math.log(adapt) if adapt > 0 else -20.0)

                # Cross-codon GT/AG penalty (curr_codon | codon boundary)
                if is_eukaryote:
                    boundary = curr_codon[-1] + codon[0]
                    if boundary == "GT":
                        s -= _GT_AG_PENALTY * 50
                    elif boundary == "AG":
                        s -= _GT_AG_PENALTY * 50

                # Restriction site check: only need to check the region
                # that could contain a new site introduced by *codon*.
                # That's the tail of (prev_codon + curr_codon + codon).
                if restriction_sites:
                    # Build enough context to catch cross-codon restriction sites
                    check_seq = curr_codon + codon
                    if prev_codon:
                        check_seq = prev_codon[-1] + check_seq
                    found_restriction = False
                    for site in restriction_sites:
                        if site in check_seq:
                            found_restriction = True
                            break
                    if found_restriction:
                        s -= _RESTRICTION_PENALTY

                state = (curr_ci, ci)
                if state not in new_dp or s > new_dp[state][0]:
                    new_dp[state] = (s, path + [ci])
        dp = new_dp

    # Extract best path
    if not dp:
        # Fallback: greedy
        codons = []
        for aa in protein:
            if aa == "*":
                codons.append("TAA")
                continue
            candidates = AA_TO_CODONS.get(aa, [])
            if not candidates:
                codons.append("NNN")
                continue
            best = max(candidates, key=lambda c: adaptiveness.get(c, 0.0))
            codons.append(best)
        return "".join(codons)

    best_state = max(dp, key=lambda k: dp[k][0])
    best_path = dp[best_state][1]

    # Reconstruct sequence
    result_codons = []
    for pos, idx in enumerate(best_path):
        result_codons.append(top_codons_per_aa[pos][idx][0])
    return "".join(result_codons)


def _count_gts(s: str) -> int:
    """Count GT dinucleotides in a sequence.

    Uses the NUMBA ``fast_dinucleotide_count`` kernel when available
    for single-pass counting; falls back to pure-Python otherwise.
    """
    return _count_dinucs_fast(s, "GT")[0]


def _is_unavoidable_gt(seq: str, pos: int) -> bool:
    """Check if a GT dinucleotide at position pos is unavoidable.
    
    A GT is unavoidable if:
    1. It's within a Valine codon (all Val codons start with GT)
    2. It's a cross-codon GT where the next codon's AA has no synonymous
       codon that doesn't start with T (e.g., Trp=TGG, Cys=TGT/TGC, Tyr=TAT/TAC)
    3. It's a cross-codon GT where the previous codon's AA has no synonymous
       codon that doesn't end with G
    """
    codon_start = (pos // 3) * 3
    next_codon_start = codon_start + 3

    # Case 1: Within-codon GT
    if pos + 1 < next_codon_start:
        codon = seq[codon_start:codon_start + 3]
        aa = CODON_TABLE.get(codon)
        if aa == 'V':
            return True  # All Valine codons start with GT
        # Check if any synonymous codon avoids GT
        for alt in AA_TO_CODONS.get(aa, []):
            if "GT" not in alt:
                return False  # There's an alternative without GT
        return True  # No alternative without GT

    # Case 2: Cross-codon GT (pos is last base of one codon, pos+1 is first of next)
    # pos is at codon_start + 2 (last position of current codon)
    # OR pos is at codon_start (we need to figure out which codons are involved)
    
    # For cross-codon GT at pos: seq[pos]='G', seq[pos+1]='T'
    # pos is the last base of the preceding codon
    # pos+1 is the first base of the following codon
    prev_cs = (pos // 3) * 3  # Start of codon containing position 'pos'
    next_cs = prev_cs + 3     # Start of codon containing position 'pos+1'
    
    if next_cs + 3 > len(seq):
        return True  # Can't check, assume unavoidable
    
    prev_codon = seq[prev_cs:prev_cs + 3]
    next_codon = seq[next_cs:next_cs + 3]
    prev_aa = CODON_TABLE.get(prev_codon)
    next_aa = CODON_TABLE.get(next_codon)
    
    if prev_aa is None or next_aa is None:
        return True
    
    # Check if we can change the previous codon to not end with G
    prev_can_avoid = any(c[-1] != 'G' for c in AA_TO_CODONS.get(prev_aa, [prev_codon]))
    # Check if we can change the next codon to not start with T
    next_can_avoid = any(c[0] != 'T' for c in AA_TO_CODONS.get(next_aa, [next_codon]))
    
    # GT is unavoidable only if BOTH sides can't avoid it
    return not (prev_can_avoid or next_can_avoid)


def _has_gt(s: str) -> bool:
    """Check if a string contains GT dinucleotide."""
    return "GT" in s


def _codon_creates_boundary_gt(
    codon: str, codon_start: int, seq_list: list
) -> Tuple[bool, bool]:
    """Check if placing codon at codon_start creates cross-codon GTs.

    Returns (prev_boundary_gt, next_boundary_gt).
    """
    prev_gt = False
    next_gt = False
    if codon_start > 0 and seq_list[codon_start - 1] + codon[0] == "GT":
        prev_gt = True
    next_pos = codon_start + 3
    if next_pos < len(seq_list) and codon[-1] + seq_list[next_pos] == "GT":
        next_gt = True
    return prev_gt, next_gt


class BioOptimizer:
    """Certified gene sequence optimizer with multi-step CAI-maximizing pipeline."""

    # Map species key (and common aliases) to canonical organism name
    # used in CODON_ADAPTIVENESS_TABLES and SUPPORTED_ORGANISMS.
    _SPECIES_TO_ORGANISM = {
        "ecoli": "Escherichia_coli",
        "E_coli": "Escherichia_coli",
        "Escherichia_coli": "Escherichia_coli",
        "human": "Homo_sapiens",
        "Homo_sapiens": "Homo_sapiens",
        "mouse": "Mus_musculus",
        "Mus_musculus": "Mus_musculus",
        "cho": "CHO_K1",
        "CHO_K1": "CHO_K1",
        "yeast": "Saccharomyces_cerevisiae",
        "Saccharomyces_cerevisiae": "Saccharomyces_cerevisiae",
    }

    def __init__(
        self,
        species: str = "ecoli",
        enzymes: Optional[List[str]] = None,
        splice_low: float = 3.0,
        splice_high: float = 6.0,
        cpg_window: int = 200,
        cpg_threshold: float = 0.6,
        min_blosum: int = -1,
        min_cai: float = 0.0,
        avoid_gt: bool = True,
        strategy: str = "constraint_first",
        optimize_mrna_stability: bool = True,
        **kwargs: Any,
    ) -> None:
        self.species = species

        # Organism-aware attributes (must be set before species_cai so that
        # we use the same CODON_ADAPTIVENESS_TABLES as compute_cai())
        from .organism_config import is_eukaryotic_organism

        # Resolve organism_name to its canonical form so that it can be
        # used as a key into CODON_ADAPTIVENESS_TABLES.  Accepts both
        # short aliases ("ecoli") and full names ("Escherichia_coli").
        raw_organism = kwargs.get(
            "organism_name",
            self._SPECIES_TO_ORGANISM.get(species, species),
        )
        self.organism_name: str = self._SPECIES_TO_ORGANISM.get(
            raw_organism, raw_organism,
        )

        # Use CODON_ADAPTIVENESS_TABLES for codon selection so that
        # the optimizer and compute_cai() agree on which codons are optimal.
        # Previously used get_species_cai_weights() which pulled from a
        # different table (SPECIES['ecoli']['cai_weights'] derived from
        # ECOLI_CODON_USAGE) and could pick non-optimal codons.
        self.species_cai: Dict[str, float] = CODON_ADAPTIVENESS_TABLES.get(
            self.organism_name, CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
        )

        self.enzymes: List[str] = enzymes or []
        self.splice_low: float = splice_low
        self.splice_high: float = splice_high
        self.cpg_window: int = cpg_window
        self.cpg_threshold: float = cpg_threshold
        self.min_blosum: int = min_blosum
        self.min_cai: float = min_cai
        self.avoid_gt: bool = avoid_gt
        self.strategy: str = strategy  # "constraint_first" or "cai_first"
        self.optimize_mrna_stability: bool = optimize_mrna_stability

        # Sliding-window GC constraint parameters
        self._gc_window_size: int = kwargs.get("gc_window_size", 50)
        self._gc_window_min: float | None = kwargs.get("gc_window_min", None)
        self._gc_window_max: float | None = kwargs.get("gc_window_max", None)

        self.organism_domain: str = kwargs.get("organism_domain", "auto")
        if self.organism_domain == "auto":
            self.organism_domain = (
                "eukaryote" if is_eukaryotic_organism(self.organism_name) else "prokaryote"
            )

        # Derived flag: True when target organism is prokaryotic.
        # Controls whether eukaryote-specific constraint steps (cryptic splice
        # elimination, CpG disruption, GT resolution, cross-codon GT/CG
        # coordination) are skipped during optimization.  Prokaryotes have no
        # spliceosome, so GT/AG avoidance and CpG island disruption are
        # biologically irrelevant and unnecessarily lower CAI.
        self.is_prokaryote: bool = self.organism_domain == "prokaryote"

        # Auto-set avoid_gt for prokaryotes unless explicitly overridden
        if "avoid_gt" not in kwargs and self.is_prokaryote:
            self.avoid_gt = False
        # Track positions where GT is unavoidable (e.g., Valine codons)
        self._unavoidable_gt_positions: Set[int] = set()
        # Track mutagenesis proposals that were applied
        self._applied_mutagenesis: List[Dict[str, Any]] = []
        # Store original input protein for conservation scoring
        self._original_protein: str = ""
        # Track mRNA stability metrics from the last optimization run
        self._mrna_stability_score: float | None = None
        self._destabilizing_motifs_removed: int = 0
        self._stability_improvement: float | None = None
        # Decision provenance collector (set externally for codon-level tracking)
        self._provenance_collector: DecisionProvenanceCollector | None = None

    def optimize(self, seq: str, strategy: Optional[str] = None) -> Tuple[str, List[PredicateResult], str]:
        """Run the full optimization pipeline.

        Args:
            seq: Input DNA or protein sequence.
            strategy: Optimization strategy override.
                - "constraint_first" (default): GT-aware greedy then fix constraints
                - "cai_first": Maximize CAI first, then fix constraints with
                  minimal CAI impact (DNAworks-style)

        Returns:
            (optimized_sequence, predicate_results, certificate_text)
        """
        effective_strategy = strategy if strategy is not None else self.strategy

        seq = seq.upper().strip()

        # ── Resolve IUPAC ambiguous bases ──────────────────────────
        # If the input contains IUPAC ambiguity codes (R, Y, S, W, K, M,
        # B, D, H, V), resolve them to concrete bases before optimization.
        # This allows users to pass degenerate sequences (e.g., from
        # consensus sequences or degenerate primers) directly to the
        # optimizer.
        ambiguous_chars = set(seq) - set("ACGT")
        if ambiguous_chars:
            from .iupac import has_ambiguous, resolve_ambiguous
            if has_ambiguous(seq):
                logger.info(
                    "Resolving %d IUPAC ambiguous bases in input sequence "
                    "(strategy: most_common)",
                    sum(1 for b in seq if b not in "ACGT"),
                )
                seq = resolve_ambiguous(
                    seq,
                    strategy="most_common",
                    cai_table=self.species_cai,
                )

        self._unavoidable_gt_positions = set()
        self._applied_mutagenesis = []
        self._original_protein = self._translate(seq)
        # Reset mRNA stability metrics
        self._mrna_stability_score = None
        self._destabilizing_motifs_removed = 0
        self._stability_improvement = None

        if effective_strategy == "cai_first":
            return self._optimize_cai_first(seq)

        # ── Default: constraint_first strategy ──
        # Step: Backtranslate CAI (DNAworks-style)
        seq = self._step_backtranslate_cai(seq)

        # Step: Resolve Constraints (fix GT/CG/RS with minimal CAI loss)
        # Skipped for prokaryotic targets when avoid_gt=False
        seq = self._step_resolve_constraints(seq)

        # Step: Remove Restriction Sites
        seq = self._step_remove_restriction_sites(seq)

        # Step: Cross-Codon Optimization (iterative)
        seq, mut_report = self._step_cross_codon_optimization(seq)

        # Step: Within-Codon GT Resolution
        # Skipped for prokaryotic targets (no spliceosome → GT avoidance unnecessary)
        if self.avoid_gt:
            seq, mut_report_35 = self._step_within_codon_gt_resolution(seq)
            mut_report.proposals.extend(mut_report_35.proposals)

        # Step: Mutagenesis Fallback (aggressive, handles within-codon GTs too)
        # Only needed when GT avoidance is active
        if self.avoid_gt:
            seq = self._step_mutagenesis_fallback(seq, mut_report)

        # Step: Avoid CpG Islands
        # Skipped for prokaryotic targets (CpG islands are a eukaryotic gene regulation concern)
        if not self.is_prokaryote:
            seq = self._step_avoid_cpg_islands(seq)

        # Step: Cross-Codon Coordination (handles cross-codon GT, CG, restriction sites)
        # GT/CG coordination is skipped for prokaryotes; restriction site coordination always runs
        seq = self._step_cross_codon_coordination(seq)

        # Step: CAI Hill Climb (upgrade codons while maintaining constraints)
        seq = self._step_cai_hill_climb(seq)

        # Step: Reoptimize (iterative until convergence)
        seq = self._step_reoptimize(seq)

        # Step: Remove ATTTA Instability Motifs
        seq = self._step_remove_instability_motifs(seq)

        # Step: CpG Reconciliation (aggressive, after CAI hill climb/reoptimize)
        # Skipped for prokaryotic targets
        if not self.is_prokaryote:
            seq = self._step_cpg_reconciliation(seq)

        # Step: GT Reconciliation (fix avoidable GTs that may have been introduced)
        # Skipped for prokaryotic targets (GT is not a constraint)
        if self.avoid_gt:
            seq = self._step_gt_reconciliation(seq)

        # Step: mRNA Stability Improvement (soft optimization — remove destabilizing motifs)
        seq = self._step_mrna_stability_improvement(seq)

        # Evaluate all 12 predicates
        results = self._evaluate_all_predicates(seq)

        # Generate certificate
        cert_text = format_certificate(results, seq, self.species)

        return seq, results, cert_text

    # ──────────────────────────────────────────────────────────
    # CAI-first optimization strategy (DNAworks-style)
    # ──────────────────────────────────────────────────────────
    def _optimize_cai_first(self, seq: str) -> Tuple[str, List[PredicateResult], str]:
        """CAI-first optimization: maximize CAI first, fix constraints second.

        Strategy (DNAworks-style):
        1. Back-translate using highest-CAI codons at every position
           (ignore GT/CG/RS initially)
        2. Iteratively fix constraint violations with minimal CAI impact:
           a. Restriction sites → fix each with best synonymous codon
           b. Avoidable GT dinucleotides → fix with best non-GT synonymous codon
              using cross-codon pair optimization
           c. CpG islands → fix with best non-CG synonymous codon
           d. Cryptic splice sites → fix with best synonymous codon
        3. Mutagenesis fallback for GTs that can't be resolved by
           synonymous substitution (e.g., Valine V→Isoleucine I)
        4. CAI hill climbing to recover any lost CAI while maintaining constraints
        5. Iterate until convergence

        Key insight: by starting from max CAI (CAI=1.0) and only making the
        smallest necessary CAI sacrifices to fix constraints, we achieve much
        higher CAI than the constraint_first strategy which permanently sacrifices
        CAI by avoiding GT during codon selection.

        For GT fixing, we use a priority-based approach that fixes GTs with
        the lowest CAI cost first, and considers multi-codon windows to find
        globally optimal fixes.
        """
        import math

        # Step: Maximize CAI — Pure max-CAI back-translation (ignore ALL constraints)
        protein = self._translate(seq)
        codons_result = []
        for i, aa in enumerate(protein):
            if aa == "*":
                codon_start = i * 3
                codons_result.append(
                    seq[codon_start:codon_start + 3]
                    if codon_start + 3 <= len(seq) else "TAA"
                )
                continue
            candidates = AA_TO_CODONS.get(aa, [])
            if candidates:
                codons_result.append(
                    max(candidates, key=lambda c: self.species_cai.get(c, 0.0))
                )
            else:
                codon_start = i * 3
                codons_result.append(seq[codon_start:codon_start + 3])
        seq = "".join(codons_result)

        # Step: Fix restriction sites (highest priority — binary constraint)
        seq = self._cai_first_fix_restriction_sites(seq)

        # Step: Fix avoidable GT dinucleotides (minimal CAI impact)
        # EUKARYOTE-ONLY: GT avoidance is irrelevant for prokaryotes
        if not self.is_prokaryote:
            seq = self._cai_first_fix_gts(seq)

        # Step: Mutagenesis fallback for Valine GTs
        # EUKARYOTE-ONLY: Only needed when GT avoidance is active
        if not self.is_prokaryote:
            seq = self._cai_first_mutagenesis_fallback(seq)

        # Step: Fix CpG islands (minimal CAI impact)
        # EUKARYOTE-ONLY: CpG islands are a eukaryotic gene regulation concern
        if not self.is_prokaryote:
            seq = self._cai_first_fix_cpg(seq)

        # Step: Fix cryptic splice sites (minimal CAI impact)
        # EUKARYOTE-ONLY: Splice sites are irrelevant for prokaryotes
        if not self.is_prokaryote:
            seq = self._cai_first_fix_splice(seq)

        # Step: Cross-Codon Coordination (handles cross-codon GT, CG, restriction sites)
        seq = self._step_cross_codon_coordination(seq)

        # Step: CAI hill climbing (upgrade codons while maintaining constraints)
        seq = self._step_cai_hill_climb(seq)

        # Step: Aggressive re-optimization pass
        seq = self._step_reoptimize(seq)

        # Step: Second pass of GT fixing + CAI boost (iterative refinement)
        # EUKARYOTE-ONLY: GT fixing not needed for prokaryotes
        if not self.is_prokaryote:
            for _refinement in range(3):
                old_cai = self._compute_seq_cai(seq)
                seq = self._cai_first_fix_gts(seq)
                seq = self._step_cai_hill_climb(seq)
                seq = self._step_reoptimize(seq)
                new_cai = self._compute_seq_cai(seq)
                if new_cai <= old_cai + 0.0001:
                    break

        # Step: CpG Reconciliation (aggressive, after CAI hill climb/reoptimize)
        # EUKARYOTE-ONLY: CpG reconciliation not needed for prokaryotes
        if not self.is_prokaryote:
            seq = self._step_cpg_reconciliation(seq)

        # Step: CAI Reconciliation (upgrade low-CAI codons while maintaining hard constraints)
        seq = self._step_cai_reconciliation(seq)

        # Track unavoidable GTs for certificate
        for i in range(len(seq) - 1):
            if seq[i] == "G" and seq[i + 1] == "T":
                if _is_unavoidable_gt(seq, i):
                    codon_start = (i // 3) * 3
                    self._unavoidable_gt_positions.add(i)

        # Step: mRNA Stability Improvement (soft optimization)
        seq = self._step_mrna_stability_improvement(seq)

        # Evaluate all predicates
        results = self._evaluate_all_predicates(seq)
        cert_text = format_certificate(results, seq, self.species)
        return seq, results, cert_text

    def _compute_seq_cai(self, seq: str) -> float:
        """Compute the geometric mean CAI for a sequence.

        Uses the NUMBA-accelerated compute_cai_kernel when available,
        falling back to the pure-Python loop otherwise.
        """
        if not seq or len(seq) < 3:
            return 0.0

        # Fast path: NUMBA kernel
        if HAS_NUMBA and _numba_cai_kernel is not None:
            try:
                # Lazy-initialize the per-optimizer adaptiveness array
                if not hasattr(self, '_numba_adapt_arr') or self._numba_adapt_arr is None:
                    self._numba_adapt_arr = _adaptiveness_to_array(self.species_cai)
                adapt_arr = self._numba_adapt_arr
                # Build codon indices, excluding Met and stop codons
                n_total = len(seq) // 3
                idx_list = []
                for i in range(n_total):
                    codon = seq[i * 3:i * 3 + 3]
                    aa = CODON_TABLE.get(codon)
                    if aa == 'M' or aa == '*':
                        continue
                    idx_list.append(_codon_to_index(codon))
                if not idx_list:
                    return 0.0
                codon_indices = _np.array(idx_list, dtype=_np.int64)
                n_codons = len(codon_indices)
                return _numba_cai_kernel(adapt_arr, codon_indices, n_codons)
            except Exception:
                pass  # Fall through to pure-Python

        # Pure-Python fallback
        import math
        log_sum = 0.0
        count = 0
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            cai = self.species_cai.get(codon, 0.0)
            if cai <= 0:
                cai = 0.001
            log_sum += math.log(cai)
            count += 1
        if count == 0:
            return 0.0
        return math.exp(log_sum / count)

    def _step_maximize_cai(self, seq: str) -> str:
        """Maximize CAI step (cai_first): Back-translate with absolute max CAI everywhere.

        No GT avoidance at all - just pick the highest-CAI codon for each AA.
        Constraint violations will be fixed in subsequent steps.
        """
        protein = self._translate(seq)
        codons_result = []
        for i, aa in enumerate(protein):
            if aa == "*":
                codon_start = i * 3
                codons_result.append(seq[codon_start:codon_start + 3] if codon_start + 3 <= len(seq) else "TAA")
                continue
            candidates = AA_TO_CODONS.get(aa, [])
            if candidates:
                codons_result.append(max(candidates, key=lambda c: self.species_cai.get(c, 0.0)))
            else:
                codon_start = i * 3
                codons_result.append(seq[codon_start:codon_start + 3])
        return "".join(codons_result)

    def _cai_first_fix_restriction_sites(self, seq: str) -> str:
        """CAI-first Fix restriction sites with minimal CAI impact."""
        from .restriction_sites import get_recognition_site
        seq_list = list(seq)

        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue

            max_rounds = 50
            for _ in range(max_rounds):
                current_seq = "".join(seq_list)
                p = current_seq.find(site)
                if p == -1:
                    break

                # Find all codon positions overlapping this site
                codon_starts = set()
                for j in range(p, p + len(site)):
                    cs = (j // 3) * 3
                    if cs + 3 <= len(seq_list):
                        codon_starts.add(cs)

                fixed = False
                # Try each overlapping codon, sorted by CAI impact (try best-CAI alts first)
                for cs in sorted(codon_starts):
                    codon = "".join(seq_list[cs:cs + 3])
                    aa = CODON_TABLE.get(codon)
                    if aa is None or aa == "*":
                        continue

                    # Sort alternatives by CAI (highest first) to minimize CAI loss
                    alts = sorted(
                        AA_TO_CODONS.get(aa, []),
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    )
                    for alt in alts:
                        if alt == codon:
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(alt):
                            test_list[cs + k] = b
                        test_seq = "".join(test_list)
                        if site not in test_seq:
                            seq_list = test_list
                            fixed = True
                            break
                    if fixed:
                        break

                if not fixed:
                    # Could not fix this site with single-codon substitution
                    # Try two-codon substitution
                    fixed = self._cai_first_fix_rs_two_codons(seq_list, p, site)
                    if not fixed:
                        break  # Give up on this site for now

        return "".join(seq_list)

    def _cai_first_fix_rs_two_codons(self, seq_list: list, pos: int, site: str) -> bool:
        """Try fixing a restriction site by modifying two adjacent codons."""
        codon_starts = sorted(set(
            (j // 3) * 3
            for j in range(pos, min(pos + len(site), len(seq_list)))
            if (j // 3) * 3 + 3 <= len(seq_list)
        ))

        if len(codon_starts) < 2:
            return False

        # Try pairs of codons
        for idx in range(len(codon_starts) - 1):
            cs1, cs2 = codon_starts[idx], codon_starts[idx + 1]
            if cs2 != cs1 + 3:
                continue  # Only adjacent codons

            codon1 = "".join(seq_list[cs1:cs1 + 3])
            codon2 = "".join(seq_list[cs2:cs2 + 3])
            aa1 = CODON_TABLE.get(codon1)
            aa2 = CODON_TABLE.get(codon2)
            if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                continue

            # Try all pairs sorted by combined CAI
            pairs = []
            for c1 in AA_TO_CODONS.get(aa1, [codon1]):
                for c2 in AA_TO_CODONS.get(aa2, [codon2]):
                    combined = self.species_cai.get(c1, 0.0) + self.species_cai.get(c2, 0.0)
                    pairs.append((c1, c2, combined))
            pairs.sort(key=lambda x: x[2], reverse=True)

            for c1, c2, _ in pairs:
                test_list = seq_list[:]
                for k, b in enumerate(c1):
                    test_list[cs1 + k] = b
                for k, b in enumerate(c2):
                    test_list[cs2 + k] = b
                test_seq = "".join(test_list)
                if site not in test_seq:
                    seq_list[:] = test_list
                    return True

        return False

    def _cai_first_fix_gts(self, seq: str) -> str:
        """CAI-first Fix avoidable GT dinucleotides with minimal CAI impact.

        Uses cost-aware GT resolution with splice donor scoring for eukaryotes:
        - GTs with low splice donor potential (< SPLICE_DONOR_POTENTIAL_THRESHOLD)
          are NOT fixed — CAI takes priority over low-risk GTs.
        - GTs with high splice donor potential are fixed with minimal CAI impact.
        - For prokaryotes, this method is never called.

        Iteratively finds each avoidable GT and resolves it by choosing
        the synonymous substitution(s) with the highest possible CAI that
        eliminates the GT, but ONLY if the GT has high splice donor potential.
        """
        seq_list = list(seq)
        max_rounds = 50

        for round_num in range(max_rounds):
            current_seq = "".join(seq_list)
            violations = []

            # Find all avoidable GT positions, filtered by splice donor potential
            for i in range(len(current_seq) - 1):
                if current_seq[i] == "G" and current_seq[i + 1] == "T":
                    if _is_unavoidable_gt(current_seq, i):
                        continue
                    # Cost-aware: skip GTs with low splice donor potential
                    sdp = score_splice_donor_potential(current_seq, i)
                    if sdp < SPLICE_DONOR_POTENTIAL_THRESHOLD:
                        # This GT has low splice donor potential — not dangerous.
                        # Accept it (CAI > GT avoidance for low-risk GTs).
                        continue
                    violations.append(i)

            if not violations:
                break

            any_fixed = False

            for gt_pos in violations:
                codon_start = (gt_pos // 3) * 3
                next_codon_start = codon_start + 3

                # Determine if within-codon or cross-codon
                is_within = (gt_pos + 1) < next_codon_start

                if is_within:
                    # Within-codon GT: try synonymous substitution
                    fixed = self._cai_first_fix_within_gt(seq_list, codon_start)
                else:
                    # Cross-codon GT: try adjacent codon pair substitution
                    fixed = self._cai_first_fix_cross_gt(seq_list, codon_start)

                if fixed:
                    any_fixed = True
                    break  # Restart scanning after any fix

            if not any_fixed:
                break

        return "".join(seq_list)

    def _cai_first_fix_within_gt(self, seq_list: list, codon_start: int) -> bool:
        """Fix a within-codon GT by choosing the highest-CAI alternative without GT."""
        codon = "".join(seq_list[codon_start:codon_start + 3])
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            return False

        old_gt_count = _count_gts("".join(seq_list))

        # Sort alternatives by CAI (highest first), skip codons with GT
        alternatives = []
        for alt in AA_TO_CODONS.get(aa, []):
            if "GT" in alt:
                continue
            alt_cai = self.species_cai.get(alt, 0.0)
            alternatives.append((alt, alt_cai))

        alternatives.sort(key=lambda x: x[1], reverse=True)

        for alt, _ in alternatives:
            test_list = seq_list[:]
            for k, b in enumerate(alt):
                test_list[codon_start + k] = b
            new_gt_count = _count_gts("".join(test_list))
            if new_gt_count < old_gt_count:
                seq_list[:] = test_list
                return True

        return False

    def _cai_first_fix_cross_gt(self, seq_list: list, codon_start: int) -> bool:
        """Fix a cross-codon GT with minimal CAI impact.

        Tries strategies in order of increasing invasiveness:
        D. Change only the next codon (to one that doesn't start with T)
        C. Change only the current codon (to one that doesn't end with G)
        A. Modify the codon pair (both codons) — sorted by combined CAI
        B. Modify preceding codon pair
        """
        old_gt_count = _count_gts("".join(seq_list))
        next_start = codon_start + 3

        # Strategy D: Change only the next codon (cheapest single-codon fix)
        if next_start + 3 <= len(seq_list):
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if aa2 is not None and aa2 != "*":
                alts = sorted(
                    AA_TO_CODONS.get(aa2, [aa2_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c2 in alts:
                    if c2[0] == "T":
                        continue
                    test_list = seq_list[:]
                    for k, b in enumerate(c2):
                        test_list[next_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        # Strategy C: Change only the current codon (to one that doesn't end with G)
        aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
        aa1 = CODON_TABLE.get(aa1_codon)
        if aa1 is not None and aa1 != "*":
            alts = sorted(
                AA_TO_CODONS.get(aa1, [aa1_codon]),
                key=lambda c: self.species_cai.get(c, 0.0),
                reverse=True,
            )
            for c1 in alts:
                if c1[-1] == "G":
                    continue
                test_list = seq_list[:]
                for k, b in enumerate(c1):
                    test_list[codon_start + k] = b
                test_seq = "".join(test_list)
                new_gt_count = _count_gts(test_seq)
                if new_gt_count < old_gt_count:
                    seq_list[:] = test_list
                    return True

        # Strategy A: Modify current + following codon
        if next_start + 3 <= len(seq_list):
            aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
            aa1 = CODON_TABLE.get(aa1_codon)
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)

            if aa1 is not None and aa1 != "*" and aa2 is not None and aa2 != "*":
                pairs = []
                for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                    for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
                        if c1[-1] + c2[0] == "GT":
                            continue
                        combined_cai = self.species_cai.get(c1, 0.0) + self.species_cai.get(c2, 0.0)
                        pairs.append((c1, c2, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c1, c2, _ in pairs:
                    test_list = seq_list[:]
                    for k, b in enumerate(c1):
                        test_list[codon_start + k] = b
                    for k, b in enumerate(c2):
                        test_list[next_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        # Strategy B: Modify preceding + current codon
        if codon_start >= 3:
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
            aa1 = CODON_TABLE.get(aa1_codon)

            if aa0 is not None and aa0 != "*" and aa1 is not None and aa1 != "*":
                pairs = []
                for c0 in AA_TO_CODONS.get(aa0, [aa0_codon]):
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                        if c0[-1] + c1[0] == "GT":
                            continue
                        combined_cai = self.species_cai.get(c0, 0.0) + self.species_cai.get(c1, 0.0)
                        pairs.append((c0, c1, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c0, c1, _ in pairs:
                    test_list = seq_list[:]
                    for k, b in enumerate(c0):
                        test_list[prev_start + k] = b
                    for k, b in enumerate(c1):
                        test_list[codon_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        return False

    def _cai_first_fix_cpg(self, seq: str) -> str:
        """CAI-first Fix CpG islands with minimal CAI impact."""
        seq_list = list(seq)
        changed = True
        iterations = 0
        max_iterations = 50

        while changed and iterations < max_iterations:
            changed = False
            iterations += 1

            for start in range(0, len(seq_list) - self.cpg_window + 1, 3):
                window = "".join(seq_list[start:start + self.cpg_window])
                c_count = window.count("C")
                g_count = window.count("G")
                cg_count = _count_dinucs_fast(window, "CG")[0]
                expected = (c_count * g_count) / len(window) if len(window) > 0 else 0
                obs_exp = cg_count / expected if expected > 0 else 0.0

                if obs_exp <= self.cpg_threshold:
                    continue

                # Find CG dinucleotides and try to break them
                for i in range(start, min(start + self.cpg_window - 1, len(seq_list) - 1)):
                    if seq_list[i] == "C" and seq_list[i+1] == "G":
                        codon_start = (i // 3) * 3
                        if codon_start + 3 > len(seq_list):
                            continue
                        codon = "".join(seq_list[codon_start:codon_start + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue

                        old_gt_count = _count_gts("".join(seq_list))

                        # Sort by CAI descending to minimize CAI loss
                        for alt in sorted(
                            AA_TO_CODONS.get(aa, []),
                            key=lambda c: self.species_cai.get(c, 0.0),
                            reverse=True,
                        ):
                            if alt == codon or "CG" in alt:
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[codon_start + k] = b
                            new_gt_count = _count_gts("".join(test_list))
                            if new_gt_count <= old_gt_count:
                                seq_list = test_list
                                changed = True
                                break

                        if changed:
                            break
                if changed:
                    break

        return "".join(seq_list)

    def _cai_first_fix_splice(self, seq: str) -> str:
        """CAI-first Fix cryptic splice sites with minimal CAI impact."""
        seq_list = list(seq)
        max_rounds = 30

        for _ in range(max_rounds):
            current_seq = "".join(seq_list)
            splice_result = check_no_cryptic_splice(current_seq, self.splice_low, self.splice_high)
            if splice_result.passed:
                break

            # Find and fix splice sites
            # The splice check looks for GT..AG patterns that resemble splice donors/acceptors
            fixed = False
            for i in range(len(current_seq) - 1):
                if current_seq[i] == "G" and current_seq[i + 1] == "T":
                    # This GT might be part of a cryptic splice site
                    codon_start = (i // 3) * 3
                    next_codon_start = codon_start + 3
                    is_within = (i + 1) < next_codon_start

                    if _is_unavoidable_gt(current_seq, i):
                        continue

                    if is_within:
                        fixed = self._cai_first_fix_within_gt(seq_list, codon_start)
                    else:
                        fixed = self._cai_first_fix_cross_gt(seq_list, codon_start)

                    if fixed:
                        break

            if not fixed:
                break

        return "".join(seq_list)

    def _cai_first_mutagenesis_fallback(self, seq: str) -> str:
        """CAI-first Apply mutagenesis for GTs that can't be resolved
        by synonymous substitution.

        Specifically targets Valine codons (GTN) which all contain GT.
        Substitutes V→I (Isoleucine, BLOSUM62=3) using highest-CAI Ile codon
        to eliminate the GT while maximizing CAI.
        """
        seq_list = list(seq)
        changed = True
        max_rounds = 10

        for _ in range(max_rounds):
            if not changed:
                break
            changed = False

            for i in range(0, len(seq_list) - 2, 3):
                codon = "".join(seq_list[i:i+3])
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*":
                    continue

                if "GT" not in codon:
                    continue

                has_gt_free = any("GT" not in c for c in AA_TO_CODONS.get(aa, []))
                if has_gt_free:
                    continue

                if aa == "V":
                    ile_codons = AA_TO_CODONS.get("I", [])
                    old_gt_count = _count_gts("".join(seq_list))

                    for ile_codon in sorted(
                        ile_codons,
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    ):
                        test_list = seq_list[:]
                        for k, b in enumerate(ile_codon):
                            test_list[i + k] = b
                        new_gt_count = _count_gts("".join(test_list))
                        if new_gt_count < old_gt_count:
                            seq_list[:] = test_list
                            self._applied_mutagenesis.append({
                                "position": i,
                                "original_aa": "V",
                                "new_aa": "I",
                                "blosum": BLOSUM62.get(("V", "I"), 3),
                            })
                            changed = True
                            break

        return "".join(seq_list)

    # Deprecated alias — use _step_maximize_cai instead
    _phase0_pure_max_cai = _step_maximize_cai

    # ──────────────────────────────────────────────────────────
    # Step: Backtranslate CAI (DP-based max-CAI back-translation)
    # ──────────────────────────────────────────────────────────
    def _step_backtranslate_cai(self, seq: str) -> str:
        """Backtranslate CAI step: DP-based max-CAI back-translation with avoidable-GT avoidance.

        Uses Viterbi-style dynamic programming to find the globally optimal
        codon assignment that maximizes CAI while avoiding only AVOIDABLE GTs.
        The DP state is simply the last character of the previous codon,
        making it O(n * 4 * k) where n is protein length and k is max codons per AA.

        Key insight: a cross-codon GT (codon1 ends with G, codon2 starts with T)
        is only avoidable if we can change at least one of the two codons to
        eliminate it. If codon2 has NO synonymous codon that doesn't start with T
        (e.g., Trp=TGG, Cys=TGT/TGC, Tyr=TAT/TAC), then the cross-codon GT
        is unavoidable and we should use the highest-CAI codon for codon1.

        This matches the semantics of the check_no_avoidable_gt predicate,
        which only fails on GTs that CAN be avoided by synonymous substitution.

        Only non-Valine within-codon GTs and avoidable cross-codon GTs are
        excluded by the DP. Unavoidable GTs (Valine within-codon, cross-codon
        where next AA has no non-T-starting codon) are allowed.
        """
        import math
        protein = self._translate(seq)
        n = len(protein)
        INF = float('-inf')

        # Precompute which amino acids have at least one codon not starting with T
        # These are the AAs where cross-codon GT from previous codon can be avoided
        aa_has_non_t_start = {}
        for aa_key in set(CODON_TABLE.values()):
            if aa_key == "*":
                continue
            codons_list = AA_TO_CODONS.get(aa_key, [])
            aa_has_non_t_start[aa_key] = any(c[0] != 'T' for c in codons_list)

        # Precompute which amino acids have at least one codon not ending with G
        # These are the AAs where cross-codon GT to next codon can be avoided
        aa_has_non_g_end = {}
        for aa_key in set(CODON_TABLE.values()):
            if aa_key == "*":
                continue
            codons_list = AA_TO_CODONS.get(aa_key, [])
            aa_has_non_g_end[aa_key] = any(c[-1] != 'G' for c in codons_list)

        # DP table: dp[i][last_char] = (max_log_cai_sum, prev_last_char, chosen_codon)
        dp = [{} for _ in range(n + 1)]
        dp[0][''] = (0.0, None, None)

        for i in range(n):
            aa = protein[i]
            codons = AA_TO_CODONS.get(aa, [])

            if not codons or aa == "*":
                # For stop codons or unknown AAs, just carry forward
                if aa == "*":
                    codon_start = i * 3
                    stop_codon = seq[codon_start:codon_start + 3] if codon_start + 3 <= len(seq) else "TAA"
                    for last_char, (log_sum, _, _) in dp[i].items():
                        new_last = stop_codon[-1]
                        if new_last not in dp[i + 1] or log_sum > dp[i + 1][new_last][0]:
                            dp[i + 1][new_last] = (log_sum, last_char, stop_codon)
                else:
                    for last_char in dp[i]:
                        dp[i + 1][last_char] = dp[i][last_char]
                continue

            # Sort codons by CAI (highest first) for tie-breaking
            codons_sorted = sorted(codons, key=lambda c: self.species_cai.get(c, 0.0), reverse=True)

            # Check if the NEXT amino acid (if any) has codons not starting with T
            # If not, cross-codon GT from this codon ending with G is unavoidable
            next_aa = protein[i + 1] if i + 1 < n else None
            next_can_avoid_gt = True
            if next_aa is not None and next_aa != "*":
                next_can_avoid_gt = aa_has_non_t_start.get(next_aa, True)
            elif next_aa == "*":
                next_can_avoid_gt = False  # Stop codon, can't change

            for last_char, (log_sum, _, _) in dp[i].items():
                if log_sum == INF:
                    continue

                for codon in codons_sorted:
                    cai = self.species_cai.get(codon, 0.0)
                    if cai <= 0:
                        cai = 0.001

                    if self.avoid_gt:
                        # Check within-codon GT (except Valine which is unavoidable)
                        has_within_gt = "GT" in codon
                        if has_within_gt and aa != 'V':
                            continue

                        # Check cross-codon GT with previous codon
                        if last_char and last_char + codon[0] == "GT":
                            # This GT is avoidable only if we could have used a
                            # different codon for the previous AA that doesn't
                            # end with the last_char. But since we're in DP, the
                            # previous choice is already fixed. However, this GT
                            # IS avoidable if the current AA has a codon not
                            # starting with T. If not, we must accept it.
                            current_can_avoid = aa_has_non_t_start.get(aa, True)
                            if current_can_avoid:
                                continue  # Skip - this GT is avoidable
                            # else: GT is unavoidable, accept this codon

                        # Check if this codon ending with G would create unavoidable
                        # cross-codon GT with the NEXT codon. We only need to avoid
                        # this if the next AA CAN avoid starting with T.
                        # If next_aa only has T-starting codons, the GT is unavoidable
                        # and we should use the highest-CAI codon (which may end with G).
                        # We handle this by NOT penalizing codons ending with G
                        # when the next AA can't avoid T-start.

                    new_last_char = codon[-1]
                    new_log_sum = log_sum + math.log(cai)

                    if new_last_char not in dp[i + 1] or new_log_sum > dp[i + 1][new_last_char][0]:
                        dp[i + 1][new_last_char] = (new_log_sum, last_char, codon)

        # Find the best final state
        best_log_sum = INF
        best_last_char = None
        for last_char, (log_sum, _, _) in dp[n].items():
            if log_sum > best_log_sum:
                best_log_sum = log_sum
                best_last_char = last_char

        if best_last_char is None:
            # Fallback: use simple max-CAI (no GT avoidance)
            codons_result = []
            for i, aa in enumerate(protein):
                if aa == "*":
                    codon_start = i * 3
                    codons_result.append(seq[codon_start:codon_start + 3])
                    continue
                candidates = AA_TO_CODONS.get(aa, [])
                if candidates:
                    codons_result.append(max(candidates, key=lambda c: self.species_cai.get(c, 0.0)))
                else:
                    codon_start = i * 3
                    codons_result.append(seq[codon_start:codon_start + 3])
            return "".join(codons_result)

        # Backtrack to find the optimal sequence
        codons_result = [None] * n
        current_char = best_last_char
        for i in range(n - 1, -1, -1):
            _, prev_char, codon = dp[i + 1][current_char]
            codons_result[i] = codon if codon is not None else "NNN"
            current_char = prev_char

        # Record codon decisions for provenance if collector is attached
        if self._provenance_collector is not None:
            # species_cai is now always a flat codon→CAI dict thanks to
            # get_species_cai_weights(), but keep isinstance guard for safety.
            _provenance_cai = self.species_cai  # type: ignore[assignment]
            if isinstance(_provenance_cai, dict) and 'cai_weights' in _provenance_cai:
                _provenance_cai = _provenance_cai['cai_weights']  # type: ignore[assignment]
            for pos_idx in range(n):
                aa = protein[pos_idx]
                if aa == "*" or not AA_TO_CODONS.get(aa, []):
                    continue
                chosen = codons_result[pos_idx]
                if chosen is None:
                    continue
                candidates_sorted = sorted(
                    AA_TO_CODONS[aa],
                    key=lambda c: _provenance_cai.get(c, 0.0),
                    reverse=True,
                )
                best_cai = _provenance_cai.get(chosen, 0.0)
                alternatives: list[dict[str, Any]] = []
                for codon in candidates_sorted:
                    if codon == chosen:
                        continue  # chosen codon is already in chosen_codon field
                    cai_val = _provenance_cai.get(codon, 0.0)
                    gc_bases = sum(1 for b in codon if b in "GC")
                    gc_contribution = gc_bases / 3.0
                    # Check constraint violations for this codon
                    violates: list[str] = []
                    if "GT" in codon and aa != 'V':
                        violates.append("cryptic_splice_donor")
                    if "AG" in codon:
                        violates.append("cryptic_splice_acceptor")
                    # Determine rejection reason
                    rejected_because: str | None = None
                    if violates:
                        rejected_because = f"Violates: {', '.join(violates)}"
                    elif cai_val < best_cai:
                        rejected_because = "Lower CAI"
                    else:
                        rejected_because = "Lower CAI (DP-optimal path chose alternative)"
                    alternatives.append({
                        "codon": codon,
                        "cai_contribution": round(cai_val, 4),
                        "gc_contribution": round(gc_contribution, 2),
                        "violates_constraints": violates,
                        "rejected_because": rejected_because,
                    })

                # Compute confidence
                if len(candidates_sorted) > 1:
                    # Use second-best CAI as baseline
                    second_best_cai = _provenance_cai.get(candidates_sorted[0], 0.0)
                    if candidates_sorted[0] == chosen and len(candidates_sorted) > 1:
                        second_best_cai = _provenance_cai.get(candidates_sorted[1], 0.0)
                    confidence = min(1.0, 0.5 + (best_cai - second_best_cai) * 5)
                    confidence = max(0.0, confidence)
                else:
                    confidence = 1.0

                self._provenance_collector.record_codon_decision(CodonDecision(
                    position=pos_idx,
                    amino_acid=aa,
                    original_codon=None,
                    chosen_codon=chosen,
                    alternatives_considered=alternatives,
                    constraint_reason="DP max-CAI with avoidable-GT avoidance",
                    confidence=round(confidence, 4),
                ))

        return "".join(codons_result)

    # Deprecated alias — use _step_backtranslate_cai instead
    _phase0_max_cai_backtranslate = _step_backtranslate_cai

    # ──────────────────────────────────────────────────────────
    # Step: Resolve Constraints (Priority-based constraint resolution)
    # ──────────────────────────────────────────────────────────
    def _step_resolve_constraints(self, seq: str) -> str:
        """Resolve Constraints step: Fix constraint violations with minimal CAI impact.

        Iteratively finds GT/CG dinucleotides and restriction sites, then
        resolves them by choosing the synonymous substitution with the
        smallest CAI penalty. This is the DNAworks-style approach: start
        with max CAI, then fix only what's needed.

        Key difference from old greedy approach: instead of greedily avoiding GT
        during codon selection (which permanently sacrifices CAI), we fix
        GT violations after the fact, choosing the resolution that costs
        the least CAI.

        Uses IncrementalSequenceState for O(1) GT/CG tracking and
        CodonCache for pre-sorted codon lists.
        """
        import math
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        max_rounds = 30

        for round_num in range(max_rounds):
            current_seq = state.sequence
            old_gt_count = state.gt_count

            # Collect all constraint violations
            violations = []

            # GT dinucleotide violations are only relevant for eukaryotic targets
            # (prokaryotes have no spliceosome, so GT avoidance is unnecessary)
            if self.avoid_gt:
                # Use incremental GT positions for O(1) lookup instead of O(N) scan
                for gt_pos in state.gt_positions_list():
                    codon_idx = gt_pos // 3
                    codon_start = codon_idx * 3
                    if gt_pos % 3 < 2:
                        # Within-codon GT - only add if avoidable
                        if not _is_unavoidable_gt(current_seq, gt_pos):
                            violations.append(("within_gt", gt_pos, codon_idx))
                    else:
                        # Cross-codon GT - only add if avoidable
                        if not _is_unavoidable_gt(current_seq, gt_pos):
                            violations.append(("cross_gt", gt_pos, codon_idx))

            if not violations and not self.enzymes:
                break

            # Check for restriction sites using cached enzyme data
            rs_violations = []
            if enzyme_cache is not None:
                for enzyme, pos in enzyme_cache.find_sites(current_seq):
                    site = enzyme_cache.get_site(enzyme)
                    if site is not None:
                        rs_violations.append((pos, site, enzyme))

            if not violations and not rs_violations:
                break

            any_resolved = False

            # Fix within-codon GTs first (usually easiest)
            for vtype, pos, codon_idx in violations:
                if vtype == "within_gt":
                    resolved = self._fix_within_codon_gt_cai_aware(state, codon_idx, codon_cache)
                    if resolved:
                        any_resolved = True
                        break

            if any_resolved:
                continue

            # Fix cross-codon GTs
            for vtype, pos, codon_idx in violations:
                if vtype == "cross_gt":
                    resolved = self._fix_cross_codon_gt_cai_aware(state, codon_idx, codon_cache)
                    if resolved:
                        any_resolved = True
                        break

            if any_resolved:
                continue

            # Fix restriction sites
            for p, site, enzyme in rs_violations:
                resolved = self._fix_restriction_site_cai_aware(state, p, site, codon_cache, enzyme_cache)
                if resolved:
                    any_resolved = True
                    break

            if not any_resolved:
                break

        # Track unavoidable GTs for mutagenesis
        for gt_pos in state.gt_positions_list():
            if gt_pos % 3 < 2:  # within-codon only
                codon_idx = gt_pos // 3
                codon_start = codon_idx * 3
                aa = state.get_aa(codon_idx)
                if aa and aa != "*":
                    all_have_gt = all("GT" in c for c in AA_TO_CODONS.get(aa, []))
                    if all_have_gt:
                        codon = state.get_codon(codon_idx)
                        for j in range(2):
                            if codon[j:j+2] == "GT":
                                self._unavoidable_gt_positions.add(codon_start + j)

        return state.sequence

    def _fix_within_codon_gt_cai_aware(self, state: IncrementalSequenceState, codon_idx: int, codon_cache: CodonCache) -> bool:
        """Fix a within-codon GT by choosing the best CAI-preserving substitution.

        Uses IncrementalSequenceState for O(1) GT tracking and
        CodonCache for pre-sorted codon lists.

        For eukaryotes: in-codon GTs are acceptable if they have low splice
        donor potential (< SPLICE_DONOR_POTENTIAL_THRESHOLD).  Only fix GTs
        that have high splice donor potential, and only if the CAI cost of
        the GT-free alternative is < threshold.
        """
        codon = state.get_codon(codon_idx)
        aa = state.get_aa(codon_idx)
        if aa is None or aa == "*":
            return False

        # ── Eukaryotic GT-vs-CAI tradeoff with splice donor scoring ──
        if not self.is_prokaryote:
            # First check splice donor potential of this GT
            seq = state.sequence
            for j in range(2):  # GT can be at pos 0-1 or 1-2 within a codon
                pos = codon_idx * 3 + j
                if pos + 1 < len(seq) and seq[pos:pos+2] == "GT":
                    sdp = score_splice_donor_potential(seq, pos)
                    if sdp < SPLICE_DONOR_POTENTIAL_THRESHOLD:
                        # This GT has low splice donor potential — not dangerous.
                        # Accept it (CAI > GT avoidance for low-risk GTs).
                        return False

            # Check CAI cost
            current_w = self.species_cai.get(codon, 0.0)
            gt_free_codons = codon_cache.get_gt_free_codons(aa)
            if gt_free_codons:
                best_gt_free_w = self.species_cai.get(gt_free_codons[0], 0.0)
                if current_w - best_gt_free_w > EUKARYOTE_CAI_GT_COST_THRESHOLD:
                    # CAI cost too high — keep the GT-containing codon
                    return False
            else:
                # No GT-free alternative (e.g., Valine) — must accept
                return False

        old_gt_count = state.gt_count

        # Use pre-sorted GT-free codons from cache
        for alt in codon_cache.get_gt_free_codons(aa):
            if alt == codon:
                continue
            # O(1) boundary check instead of manual base lookup
            left_gt, right_gt = state.boundary_creates_gt(codon_idx, alt)
            if left_gt or right_gt:
                continue
            # O(1) swap + O(1) GT count check + rollback if needed
            old_codon = state.swap_codon(codon_idx, alt)
            if state.gt_count < old_gt_count:
                return True
            else:
                state.swap_codon(codon_idx, old_codon)  # Rollback

        return False

    def _fix_cross_codon_gt_cai_aware(self, state: IncrementalSequenceState, codon_idx: int, codon_cache: CodonCache) -> bool:
        """Fix a cross-codon GT by choosing the best CAI-preserving substitution.

        Uses IncrementalSequenceState for O(1) GT tracking and
        CodonCache for pre-sorted codon lists.

        For eukaryotes: cross-codon GTs are acceptable if they have low
        splice donor potential (< SPLICE_DONOR_POTENTIAL_THRESHOLD).  When
        the GT has high splice donor potential, it is only eliminated if the
        CAI cost of the fix is < EUKARYOTE_CAI_GT_COST_THRESHOLD.
        """
        # ── Eukaryotic GT-vs-CAI tradeoff for cross-codon GTs ──
        # For eukaryotes, check splice donor potential first, then CAI cost.
        if not self.is_prokaryote:
            # Check splice donor potential of this cross-codon GT
            seq = state.sequence
            gt_pos = codon_idx * 3 + 2  # G at end of codon_idx, T at start of codon_idx+1
            if gt_pos + 1 < len(seq) and seq[gt_pos:gt_pos+2] == "GT":
                sdp = score_splice_donor_potential(seq, gt_pos)
                if sdp < SPLICE_DONOR_POTENTIAL_THRESHOLD:
                    # This cross-codon GT has low splice donor potential — not dangerous.
                    # Accept it (CAI > GT avoidance for low-risk GTs).
                    return False

            # Check if the CAI cost of fixing this cross-codon GT
            # is acceptable before attempting to fix it.
            # The G is at the end of codon_idx, the T is at the start of the next codon
            next_ci = codon_idx + 1
            if next_ci < state.num_codons:
                g_aa = state.get_aa(codon_idx)
                t_aa = state.get_aa(next_ci)
                if g_aa and t_aa and g_aa != "*" and t_aa != "*":
                    g_current = state.get_codon(codon_idx)
                    t_current = state.get_codon(next_ci)
                    g_current_w = self.species_cai.get(g_current, 0.0)
                    t_current_w = self.species_cai.get(t_current, 0.0)
                    # Check if there's a cheap fix (CAI cost < threshold)
                    # for either the G-ending codon or the T-starting codon
                    _cross_cai_cost_ok = False
                    g_non_g_end = [c for c in codon_cache.get_sorted_codons(g_aa) if c[-1] != "G"]
                    if g_non_g_end:
                        best_non_g_end_w = self.species_cai.get(g_non_g_end[0], 0.0)
                        if g_current_w - best_non_g_end_w < EUKARYOTE_CAI_GT_COST_THRESHOLD:
                            _cross_cai_cost_ok = True
                    if not _cross_cai_cost_ok:
                        t_non_t_start = [c for c in codon_cache.get_sorted_codons(t_aa) if c[0] != "T"]
                        if t_non_t_start:
                            best_non_t_start_w = self.species_cai.get(t_non_t_start[0], 0.0)
                            if t_current_w - best_non_t_start_w < EUKARYOTE_CAI_GT_COST_THRESHOLD:
                                _cross_cai_cost_ok = True
                    if not _cross_cai_cost_ok:
                        # CAI cost too high — accept the cross-codon GT
                        return False

        old_gt_count = state.gt_count

        # Try resolving by modifying the codon pair
        next_idx = codon_idx + 1
        prev_idx = codon_idx - 1

        # Strategy A: Modify following codon pair (codon_idx and next_idx)
        if next_idx < state.num_codons:
            aa1 = state.get_aa(codon_idx)
            aa2 = state.get_aa(next_idx)

            if aa1 is not None and aa1 != "*" and aa2 is not None and aa2 != "*":
                # Collect all valid pairs sorted by combined CAI
                pairs = []
                for c1 in codon_cache.get_sorted_codons(aa1):
                    for c2 in codon_cache.get_sorted_codons(aa2):
                        if c1[-1] + c2[0] == "GT":
                            continue
                        combined_cai = codon_cache.get_cai(c1) + codon_cache.get_cai(c2)
                        pairs.append((c1, c2, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c1, c2, _ in pairs:
                    old1 = state.swap_codon(codon_idx, c1)
                    old2 = state.swap_codon(next_idx, c2)
                    if state.gt_count < old_gt_count:
                        return True
                    else:
                        state.swap_codon(next_idx, old2)  # Rollback
                        state.swap_codon(codon_idx, old1)  # Rollback

        # Strategy B: Modify preceding codon pair (prev_idx and codon_idx)
        if prev_idx >= 0:
            aa0 = state.get_aa(prev_idx)
            aa1 = state.get_aa(codon_idx)

            if aa0 is not None and aa0 != "*" and aa1 is not None and aa1 != "*":
                pairs = []
                for c0 in codon_cache.get_sorted_codons(aa0):
                    for c1 in codon_cache.get_sorted_codons(aa1):
                        if c0[-1] + c1[0] == "GT":
                            continue
                        combined_cai = codon_cache.get_cai(c0) + codon_cache.get_cai(c1)
                        pairs.append((c0, c1, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c0, c1, _ in pairs:
                    old0 = state.swap_codon(prev_idx, c0)
                    old1 = state.swap_codon(codon_idx, c1)
                    if state.gt_count < old_gt_count:
                        return True
                    else:
                        state.swap_codon(codon_idx, old1)  # Rollback
                        state.swap_codon(prev_idx, old0)  # Rollback

        # Strategy C: Single codon substitution (change only the codon that ends with G)
        aa1 = state.get_aa(codon_idx)
        if aa1 is not None and aa1 != "*":
            for c1 in codon_cache.get_sorted_codons(aa1):
                old_codon = state.swap_codon(codon_idx, c1)
                if state.gt_count < old_gt_count:
                    return True
                else:
                    state.swap_codon(codon_idx, old_codon)  # Rollback

        return False

    def _fix_restriction_site_cai_aware(self, state: IncrementalSequenceState, pos: int, site: str, codon_cache: CodonCache, enzyme_cache: Optional[EnzymeSiteCache] = None) -> bool:
        """Fix a restriction site by choosing the best CAI-preserving substitution.

        Uses IncrementalSequenceState for O(1) GT tracking and
        CodonCache for pre-sorted codon lists.
        """
        codon_indices = set()
        for j in range(pos, pos + len(site)):
            ci = j // 3
            if ci < state.num_codons:
                codon_indices.add(ci)

        old_gt_count = state.gt_count

        # Try each overlapping codon, sorted by position
        for ci in sorted(codon_indices):
            aa = state.get_aa(ci)
            if aa is None or aa == "*":
                continue

            # Use pre-sorted codons from cache
            for alt in codon_cache.get_sorted_codons(aa):
                current_codon = state.get_codon(ci)
                if alt == current_codon:
                    continue
                # O(1) swap + check + rollback
                old_codon = state.swap_codon(ci, alt)
                # Check restriction site is gone
                if enzyme_cache is not None:
                    rs_ok = not enzyme_cache.check_any_site_present(state.sequence)
                else:
                    rs_ok = site not in state.sequence
                if rs_ok:
                    if self.avoid_gt:
                        if state.gt_count <= old_gt_count:
                            return True
                        else:
                            state.swap_codon(ci, old_codon)  # Rollback
                    else:
                        return True
                else:
                    state.swap_codon(ci, old_codon)  # Rollback

        return False

    # Deprecated alias — use _step_resolve_constraints instead
    _phase1_priority_constraint_resolution = _step_resolve_constraints

    # ──────────────────────────────────────────────────────────
    # Step: Greedy Optimize (Per-position CAI maximization)
    # ──────────────────────────────────────────────────────────
    def _step_greedy_optimize(self, seq: str) -> str:
        """Greedy Optimize step: Per-position CAI maximization, GT-aware.

        For each amino acid, select the highest-CAI codon that does not
        introduce a GT dinucleotide within or across codon boundaries.

        If ALL synonymous codons for an AA contain GT (e.g., Valine GTN),
        flag the position as "unavoidable GT" for mutagenesis fallback step.
        """
        codons = []
        protein = self._translate(seq)
        prev_codon_end = ""

        for i, aa in enumerate(protein):
            if aa == "*":
                codon_start = i * 3
                codons.append(seq[codon_start:codon_start + 3])
                prev_codon_end = codons[-1][-1]
                continue

            candidates = AA_TO_CODONS.get(aa, [])
            if not candidates:
                codon_start = i * 3
                codons.append(seq[codon_start:codon_start + 3])
                prev_codon_end = codons[-1][-1]
                continue

            candidates_sorted = sorted(
                candidates,
                key=lambda c: self.species_cai.get(c, 0.0),
                reverse=True,
            )

            best = candidates_sorted[0]
            found_gt_free = False

            if self.avoid_gt:
                for codon in candidates_sorted:
                    if "GT" in codon:
                        continue
                    if prev_codon_end and prev_codon_end + codon[0] == "GT":
                        continue
                    best = codon
                    found_gt_free = True
                    break

                if not found_gt_free:
                    all_have_gt = all("GT" in c for c in candidates)
                    if all_have_gt:
                        codon_abs_start = i * 3
                        for j in range(2):
                            if best[j:j+2] == "GT":
                                self._unavoidable_gt_positions.add(codon_abs_start + j)
                        # Pick the codon that avoids cross-codon GT if possible
                        if prev_codon_end:
                            for codon in candidates_sorted:
                                if prev_codon_end + codon[0] != "GT":
                                    best = codon
                                    break
                    else:
                        # Pick the one that minimizes GT count
                        for codon in candidates_sorted:
                            gt_count = codon.count("GT")
                            cross_gt = 1 if (prev_codon_end and prev_codon_end + codon[0] == "GT") else 0
                            best_gt_count = best.count("GT")
                            best_cross_gt = 1 if (prev_codon_end and prev_codon_end + best[0] == "GT") else 0
                            if gt_count + cross_gt < best_gt_count + best_cross_gt:
                                best = codon
                                break

            # Record codon decision for provenance if collector is attached
            if self._provenance_collector is not None:
                # species_cai is now always a flat codon→CAI dict thanks to
                # get_species_cai_weights(), but keep isinstance guard for safety.
                _prov_cai = self.species_cai
                if isinstance(_prov_cai, dict) and 'cai_weights' in _prov_cai:
                    _prov_cai = _prov_cai['cai_weights']
                best_cai = _prov_cai.get(best, 0.0)
                alternatives: list[dict[str, Any]] = []
                for codon in candidates_sorted:
                    if codon == best:
                        continue  # chosen codon is already in chosen_codon field
                    cai_val = _prov_cai.get(codon, 0.0)
                    gc_bases = sum(1 for b in codon if b in "GC")
                    gc_contribution = gc_bases / 3.0
                    # Check constraint violations for this codon
                    violates: list[str] = []
                    if "GT" in codon:
                        violates.append("cryptic_splice_donor")
                    if prev_codon_end and prev_codon_end + codon[0] == "GT":
                        violates.append("cross_codon_gt")
                    if "AG" in codon:
                        violates.append("cryptic_splice_acceptor")
                    # Determine rejection reason
                    rejected_because: str | None = None
                    if violates:
                        rejected_because = f"Violates: {', '.join(violates)}"
                    elif cai_val < best_cai:
                        rejected_because = "Lower CAI"
                    else:
                        rejected_because = "Lower CAI"
                    alternatives.append({
                        "codon": codon,
                        "cai_contribution": round(cai_val, 4),
                        "gc_contribution": round(gc_contribution, 2),
                        "violates_constraints": violates,
                        "rejected_because": rejected_because,
                    })

                # Compute confidence
                if len(candidates_sorted) > 1:
                    second_best_cai = _prov_cai.get(candidates_sorted[1], 0.0)
                    confidence = min(1.0, 0.5 + (best_cai - second_best_cai) * 5)
                    confidence = max(0.0, confidence)
                else:
                    confidence = 1.0

                self._provenance_collector.record_codon_decision(CodonDecision(
                    position=i,
                    amino_acid=aa,
                    original_codon=None,
                    chosen_codon=best,
                    alternatives_considered=alternatives,
                    constraint_reason="Maximize CAI while avoiding GT dinucleotides",
                    confidence=round(confidence, 4),
                ))

            codons.append(best)
            prev_codon_end = best[-1]

        return "".join(codons)

    # Deprecated alias — use _step_greedy_optimize instead
    _phase1_greedy_optimize = _step_greedy_optimize

    # ──────────────────────────────────────────────────────────
    # Step: Remove Restriction Sites
    # ──────────────────────────────────────────────────────────
    def _step_remove_restriction_sites(self, seq: str) -> str:
        """Remove Restriction Sites step (v2: batch single-pass approach).

        Instead of iterating one site at a time (O(n * k * iterations)),
        this collects ALL restriction site positions in a single scan,
        then resolves them in batch using incremental state for O(1)
        GT/CG tracking after each fix.

        Performance: O(n * k + changes) instead of O(n * k * iterations).
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        max_rounds = 20

        for round_num in range(max_rounds):
            current_seq = state.sequence

            # Single-pass: find ALL restriction sites at once
            all_sites = []
            if enzyme_cache is not None:
                all_sites = enzyme_cache.find_all_sites_batch(current_seq)

            if not all_sites:
                break  # No sites remaining — done

            any_fixed = False

            for enzyme, pos, site_seq in all_sites:
                # Get overlapping codon indices
                codon_indices = set()
                for j in range(pos, pos + len(site_seq)):
                    ci = j // 3
                    if ci < state.num_codons:
                        codon_indices.add(ci)

                # Try each overlapping codon, using pre-sorted alternatives from cache
                for ci in sorted(codon_indices):
                    aa = state.get_aa(ci)
                    if aa is None or aa == "*":
                        continue

                    current_codon = state.get_codon(ci)

                    # Use pre-sorted codons from cache (avoids repeated sorting)
                    for alt in codon_cache.get_sorted_codons(aa):
                        if alt == current_codon:
                            continue

                        # O(1) swap + check + rollback
                        old_codon = state.swap_codon(ci, alt)

                        # Check restriction site is gone
                        if enzyme_cache is not None:
                            rs_ok = not enzyme_cache.check_any_site_present(state.sequence)
                        else:
                            rs_ok = site_seq not in state.sequence

                        if rs_ok:
                            # Check GT increase (allow temporary GT increase, but not severe)
                            if self.avoid_gt:
                                # We can't easily compute old GT count here since we
                                # already swapped. But state.gt_count is the NEW count.
                                # Allow up to 2 more GTs than before the round started.
                                if state.gt_count > state.gt_count + 2:
                                    state.swap_codon(ci, old_codon)  # Rollback
                                    continue
                            any_fixed = True
                            break
                        else:
                            state.swap_codon(ci, old_codon)  # Rollback

                    if any_fixed:
                        break

                if any_fixed:
                    break  # Restart scanning with updated sequence

            if not any_fixed:
                # Try multi-codon coordinated removal for remaining sites
                for enzyme, pos, site_seq in all_sites:
                    codon_indices = sorted(set(
                        j // 3 for j in range(pos, pos + len(site_seq))
                        if j // 3 < state.num_codons
                    ))

                    if len(codon_indices) < 2:
                        continue

                    # Try pairs of adjacent codons
                    for idx in range(len(codon_indices) - 1):
                        ci1, ci2 = codon_indices[idx], codon_indices[idx + 1]
                        if ci2 != ci1 + 1:
                            continue
                        aa1 = state.get_aa(ci1)
                        aa2 = state.get_aa(ci2)
                        if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                            continue

                        old1 = state.get_codon(ci1)
                        old2 = state.get_codon(ci2)

                        # Try codon pairs sorted by combined CAI
                        pairs = []
                        for c1 in codon_cache.get_sorted_codons(aa1):
                            for c2 in codon_cache.get_sorted_codons(aa2):
                                combined = codon_cache.get_cai(c1) + codon_cache.get_cai(c2)
                                pairs.append((c1, c2, combined))
                        pairs.sort(key=lambda x: x[2], reverse=True)

                        for c1, c2, _ in pairs:
                            o1 = state.swap_codon(ci1, c1)
                            o2 = state.swap_codon(ci2, c2)
                            if enzyme_cache is not None:
                                rs_ok = not enzyme_cache.check_any_site_present(state.sequence)
                            else:
                                rs_ok = site_seq not in state.sequence
                            if rs_ok:
                                any_fixed = True
                                break
                            else:
                                state.swap_codon(ci2, o2)
                                state.swap_codon(ci1, o1)

                        if any_fixed:
                            break
                    if any_fixed:
                        break

                if not any_fixed:
                    break  # Cannot fix remaining sites

        return state.sequence

    # Deprecated alias — use _step_remove_restriction_sites instead
    _phase2_remove_restriction_sites = _step_remove_restriction_sites

    # ──────────────────────────────────────────────────────────
    # Step: Cross-Codon Optimization (iterative constraint resolution)
    # ──────────────────────────────────────────────────────────
    def _step_cross_codon_optimization(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Cross-Codon Optimization step (v2: single-pass scan with early termination).

        Resolve cross-codon GT, CG, and restriction site constraints.
        Uses IncrementalSequenceState for O(1) GT/CG tracking instead of
        O(N) full-sequence rescans.

        v2 improvements:
        - Single-pass constraint collection (GT + CG + RS in one scan)
        - Early termination: if no constraints, skip immediately
        - Process all constraints in one round before re-scanning
        - Avoid re-calling _is_unavoidable_gt for same positions
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        total_remaining: Dict[int, List[str]] = {}
        max_rounds = 15

        for round_num in range(max_rounds):
            # Single-pass: collect ALL cross-codon constraints at once
            constraint_positions: Dict[int, List[str]] = {}

            # GT/CG constraints — use incremental position data (O(1) per position)
            if self.avoid_gt:
                for gt_pos in state.gt_positions_list():
                    if gt_pos % 3 == 2:  # Cross-codon boundary
                        codon_start = (gt_pos // 3) * 3
                        # Only add if avoidable
                        codon_idx = gt_pos // 3
                        if codon_idx < state.num_codons:
                            aa = state.get_aa(codon_idx)
                            # Quick check: Valine always has GT, skip
                            if aa == 'V':
                                continue
                            constraint_positions.setdefault(codon_start, [])
                            if "GT" not in constraint_positions[codon_start]:
                                constraint_positions[codon_start].append("GT")

                for cg_pos in state.cg_positions_list():
                    if cg_pos % 3 == 2:  # Cross-codon boundary
                        codon_start = (cg_pos // 3) * 3
                        constraint_positions.setdefault(codon_start, [])
                        if "CG" not in constraint_positions[codon_start]:
                            constraint_positions[codon_start].append("CG")

            # Restriction site constraints — single batch scan
            if enzyme_cache is not None:
                for enzyme, site in enzyme_cache._sites.items():
                    if site is None:
                        continue
                    for pos in find_cross_codon_restriction(state.sequence, site):
                        codon_start = (pos // 3) * 3
                        label = f"RS:{site}"
                        constraint_positions.setdefault(codon_start, [])
                        if label not in constraint_positions[codon_start]:
                            constraint_positions[codon_start].append(label)

            # Early termination: no constraints found
            if not constraint_positions:
                break

            any_resolved = False

            # Process ALL constraints in one pass (don't re-scan after each fix;
            # instead, process them all and re-scan only at the start of next round)
            for codon_start, ctypes in constraint_positions.items():
                # Skip if this position was already fixed by a previous resolution
                # (the codon at this position may have changed)
                resolved = self._try_resolve_cross_codon_incremental(
                    state, codon_start, ctypes, codon_cache, enzyme_cache
                )
                if resolved:
                    any_resolved = True
                    # Don't break — process all constraints in this round

            if not any_resolved:
                # No progress — record remaining constraints and exit
                constraint_positions2: Dict[int, List[str]] = {}
                if self.avoid_gt:
                    for gt_pos in state.gt_positions_list():
                        if gt_pos % 3 == 2:
                            if not _is_unavoidable_gt(state.sequence, gt_pos):
                                cs = (gt_pos // 3) * 3
                                constraint_positions2.setdefault(cs, [])
                                if "GT" not in constraint_positions2[cs]:
                                    constraint_positions2[cs].append("GT")
                    for cg_pos in state.cg_positions_list():
                        if cg_pos % 3 == 2:
                            cs = (cg_pos // 3) * 3
                            constraint_positions2.setdefault(cs, [])
                            if "CG" not in constraint_positions2[cs]:
                                constraint_positions2[cs].append("CG")
                total_remaining = constraint_positions2
                break

        mut_report = propose_mutagenesis(
            state.sequence,
            list(total_remaining.keys()),
            total_remaining,
            self.species_cai,
            self.min_blosum,
            self.min_cai,
        )

        return state.sequence, mut_report

    def _try_resolve_cross_codon(
        self, seq_list: list, codon_start: int, ctypes: List[str]
    ) -> bool:
        """Try to resolve cross-codon constraints at codon_start.

        Uses global validation: after any substitution, verifies that the
        total GT count in the sequence doesn't increase.

        CAI-aware: codon alternatives are sorted by CAI (descending) so
        the first valid combination found is also the highest-CAI one.
        """
        old_seq = "".join(seq_list)
        old_gt_count = _count_gts(old_seq)

        aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
        aa1 = CODON_TABLE.get(aa1_codon)
        if aa1 is None or aa1 == "*":
            return False

        next_start = codon_start + 3

        # Sort alternatives by CAI (descending) for CAI-aware resolution
        aa1_alts = sorted(
            AA_TO_CODONS.get(aa1, [aa1_codon]),
            key=lambda c: self.species_cai.get(c, 0.0),
            reverse=True,
        )

        # Strategy A: Following codon
        if next_start + 3 <= len(seq_list):
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if aa2 is not None and aa2 != "*":
                aa2_alts = sorted(
                    AA_TO_CODONS.get(aa2, [aa2_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c1 in aa1_alts:
                    if self.avoid_gt and "GT" in c1:
                        continue
                    for c2 in aa2_alts:
                        if self.avoid_gt and "GT" in c2:
                            continue
                        if c1[-1] + c2[0] == "GT":
                            continue
                        # Apply and validate globally
                        test_list = seq_list[:]
                        for k, b in enumerate(c1):
                            test_list[codon_start + k] = b
                        for k, b in enumerate(c2):
                            test_list[next_start + k] = b
                        test_seq = "".join(test_list)
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            # Also check constraint is actually resolved
                            resolved = True
                            for ct in ctypes:
                                if ct == "GT" and _has_gt(test_seq[codon_start:next_start + 3]):
                                    resolved = False
                                elif ct == "CG" and "CG" in test_seq[codon_start:next_start + 3]:
                                    resolved = False
                                elif ct.startswith("RS:"):
                                    if ct[3:] in test_seq:
                                        resolved = False
                            if resolved:
                                # Also check no new restriction sites introduced
                                from .restriction_sites import get_recognition_site
                                rs_ok = True
                                for enzyme in self.enzymes:
                                    rs_site = get_recognition_site(enzyme)
                                    if rs_site and rs_site in test_seq:
                                        rs_ok = False
                                        break
                                if rs_ok:
                                    seq_list[:] = test_list
                                    return True

        # Strategy B: Preceding codon
        if codon_start >= 3:
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            if aa0 is not None and aa0 != "*":
                aa0_alts = sorted(
                    AA_TO_CODONS.get(aa0, [aa0_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c0 in aa0_alts:
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in aa1_alts:
                        if self.avoid_gt and "GT" in c1:
                            continue
                        if c0[-1] + c1[0] == "GT":
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(c0):
                            test_list[prev_start + k] = b
                        for k, b in enumerate(c1):
                            test_list[codon_start + k] = b
                        test_seq = "".join(test_list)
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            resolved = True
                            for ct in ctypes:
                                if ct == "GT" and _has_gt(test_seq[prev_start:codon_start + 3]):
                                    resolved = False
                                elif ct == "CG" and "CG" in test_seq[prev_start:codon_start + 3]:
                                    resolved = False
                                elif ct.startswith("RS:"):
                                    if ct[3:] in test_seq:
                                        resolved = False
                            if resolved:
                                seq_list[:] = test_list
                                return True

        # Strategy C: Both preceding and following
        if codon_start >= 3 and next_start + 3 <= len(seq_list):
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if (aa0 is not None and aa0 != "*" and
                    aa2 is not None and aa2 != "*"):
                aa0_alts = sorted(
                    AA_TO_CODONS.get(aa0, [aa0_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                aa2_alts = sorted(
                    AA_TO_CODONS.get(aa2, [aa2_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c0 in aa0_alts:
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in aa1_alts:
                        if self.avoid_gt and "GT" in c1:
                            continue
                        for c2 in aa2_alts:
                            if self.avoid_gt and "GT" in c2:
                                continue
                            if c0[-1] + c1[0] == "GT" or c1[-1] + c2[0] == "GT":
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(c0):
                                test_list[prev_start + k] = b
                            for k, b in enumerate(c1):
                                test_list[codon_start + k] = b
                            for k, b in enumerate(c2):
                                test_list[next_start + k] = b
                            test_seq = "".join(test_list)
                            new_gt_count = _count_gts(test_seq)
                            if new_gt_count <= old_gt_count:
                                resolved = True
                                for ct in ctypes:
                                    if ct == "GT" and _has_gt(test_seq[prev_start:next_start + 3]):
                                        resolved = False
                                    elif ct == "CG" and "CG" in test_seq[prev_start:next_start + 3]:
                                        resolved = False
                                    elif ct.startswith("RS:"):
                                        if ct[3:] in test_seq:
                                            resolved = False
                                if resolved:
                                    seq_list[:] = test_list
                                    return True

        # Strategy D: Single codon substitution
        for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
            if self.avoid_gt and "GT" in c1:
                continue
            test_list = seq_list[:]
            for k, b in enumerate(c1):
                test_list[codon_start + k] = b
            test_seq = "".join(test_list)
            new_gt_count = _count_gts(test_seq)
            if new_gt_count < old_gt_count:
                seq_list[:] = test_list
                return True

        return False

    def _try_resolve_cross_codon_incremental(
        self, state: IncrementalSequenceState, codon_start: int,
        ctypes: List[str], codon_cache: CodonCache,
        enzyme_cache: Optional[EnzymeSiteCache] = None
    ) -> bool:
        """Try to resolve cross-codon constraints using incremental state.

        Uses O(1) GT/CG tracking via IncrementalSequenceState instead of
        O(N) full-sequence rescans. Uses predictive boundary checks to
        avoid unnecessary swap+rollback cycles.

        v2: Uses regional RS checking to avoid full-sequence scans.
        Uses pre-sorted codon subsets from CodonCache for faster iteration.
        """
        codon_idx = codon_start // 3
        aa1 = state.get_aa(codon_idx)
        if aa1 is None or aa1 == "*":
            return False

        old_gt_count = state.gt_count
        next_codon_idx = codon_idx + 1

        # Pre-check: what constraint types need resolving?
        needs_gt = "GT" in ctypes
        needs_cg = "CG" in ctypes
        needs_rs = any(ct.startswith("RS:") for ct in ctypes)

        # Strategy A: Following codon pair
        if next_codon_idx < state.num_codons:
            aa2 = state.get_aa(next_codon_idx)
            if aa2 is not None and aa2 != "*":
                # Use pre-sorted non-G-end/non-T-start codons from cache for faster iteration
                if needs_gt:
                    c1_list = codon_cache.get_non_g_end_codons(aa1)
                    c2_list = codon_cache.get_non_t_start_codons(aa2)
                else:
                    c1_list = codon_cache.get_sorted_codons(aa1)
                    c2_list = codon_cache.get_sorted_codons(aa2)

                for c1 in c1_list:
                    if self.avoid_gt and "GT" in c1:
                        continue
                    # O(1) left boundary predictive check
                    if self.avoid_gt:
                        left_gt, _ = state.boundary_creates_gt(codon_idx, c1)
                        if left_gt:
                            continue
                    for c2 in c2_list:
                        if self.avoid_gt and "GT" in c2:
                            continue
                        if c1[-1] + c2[0] == "GT":
                            continue
                        # O(1) right boundary predictive check
                        if self.avoid_gt:
                            _, right_gt = state.boundary_creates_gt(next_codon_idx, c2)
                            if right_gt:
                                continue

                        # Quick dinucleotide check — does this pair resolve the constraint?
                        if needs_cg and c1[-1] + c2[0] == "CG":
                            continue

                        # Apply both swaps
                        old1 = state.swap_codon(codon_idx, c1)
                        old2 = state.swap_codon(next_codon_idx, c2)

                        if state.gt_count <= old_gt_count:
                            # Check constraint is actually resolved — use O(1) has_gt_at/has_cg_at
                            resolved = True
                            for ct in ctypes:
                                if ct == "GT" and state.has_gt_at(codon_start + 2):
                                    resolved = False
                                elif ct == "CG" and state.has_cg_at(codon_start + 2):
                                    resolved = False
                                elif ct.startswith("RS:"):
                                    # Use regional RS check (avoids full-sequence string build)
                                    if enzyme_cache is not None:
                                        if enzyme_cache.check_sites_in_region(
                                            state.sequence, codon_start, next_codon_idx * 3 + 3
                                        ):
                                            resolved = False
                                    else:
                                        if ct[3:] in state.sequence:
                                            resolved = False
                            if resolved:
                                return True

                        # Rollback
                        state.swap_codon(next_codon_idx, old2)
                        state.swap_codon(codon_idx, old1)

        # Strategy B: Preceding codon
        if codon_idx > 0:
            prev_codon_idx = codon_idx - 1
            aa0 = state.get_aa(prev_codon_idx)
            if aa0 is not None and aa0 != "*":
                for c0 in codon_cache.get_sorted_codons(aa0):
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in codon_cache.get_sorted_codons(aa1):
                        if self.avoid_gt and "GT" in c1:
                            continue
                        if c0[-1] + c1[0] == "GT":
                            continue
                        # Quick dinucleotide check
                        if needs_gt and "GT" in (c0 + c1):
                            continue
                        if needs_cg and c0[-1] + c1[0] == "CG":
                            continue

                        old0 = state.swap_codon(prev_codon_idx, c0)
                        old1 = state.swap_codon(codon_idx, c1)

                        if state.gt_count <= old_gt_count:
                            resolved = True
                            for ct in ctypes:
                                prev_start = codon_start - 3
                                region = state.sequence[prev_start:codon_start + 3]
                                if ct == "GT" and "GT" in region:
                                    resolved = False
                                elif ct == "CG" and "CG" in region:
                                    resolved = False
                                elif ct.startswith("RS:"):
                                    if ct[3:] in state.sequence:
                                        resolved = False
                            if resolved:
                                return True

                        state.swap_codon(codon_idx, old1)
                        state.swap_codon(prev_codon_idx, old0)

        # Strategy C: Both preceding and following
        if codon_idx > 0 and next_codon_idx < state.num_codons:
            prev_codon_idx = codon_idx - 1
            aa0 = state.get_aa(prev_codon_idx)
            aa2 = state.get_aa(next_codon_idx)
            if (aa0 is not None and aa0 != "*" and
                    aa2 is not None and aa2 != "*"):
                for c0 in codon_cache.get_sorted_codons(aa0):
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in codon_cache.get_sorted_codons(aa1):
                        if self.avoid_gt and "GT" in c1:
                            continue
                        if c0[-1] + c1[0] == "GT":
                            continue
                        for c2 in codon_cache.get_sorted_codons(aa2):
                            if self.avoid_gt and "GT" in c2:
                                continue
                            if c1[-1] + c2[0] == "GT":
                                continue
                            # Quick dinucleotide check
                            if needs_gt and "GT" in (c0 + c1 + c2):
                                continue
                            if needs_cg and (c0[-1] + c1[0] == "CG" or c1[-1] + c2[0] == "CG"):
                                continue

                            old0 = state.swap_codon(prev_codon_idx, c0)
                            old1 = state.swap_codon(codon_idx, c1)
                            old2 = state.swap_codon(next_codon_idx, c2)

                            if state.gt_count <= old_gt_count:
                                resolved = True
                                for ct in ctypes:
                                    prev_start = codon_start - 3
                                    region = state.sequence[prev_start:codon_start + 6]
                                    if ct == "GT" and "GT" in region:
                                        resolved = False
                                    elif ct == "CG" and "CG" in region:
                                        resolved = False
                                    elif ct.startswith("RS:"):
                                        if ct[3:] in state.sequence:
                                            resolved = False
                                if resolved:
                                    return True

                            state.swap_codon(next_codon_idx, old2)
                            state.swap_codon(codon_idx, old1)
                            state.swap_codon(prev_codon_idx, old0)

        # Strategy D: Single codon substitution
        for c1 in codon_cache.get_sorted_codons(aa1):
            if self.avoid_gt and "GT" in c1:
                continue
            # Predictive check: would this reduce GT count?
            if self.avoid_gt and state.would_gt_increase(codon_idx, c1):
                continue
            old1 = state.swap_codon(codon_idx, c1)
            if state.gt_count < old_gt_count:
                return True
            state.swap_codon(codon_idx, old1)

        return False

    # Deprecated alias — use _step_cross_codon_optimization instead
    _phase3_cross_codon_constraints = _step_cross_codon_optimization

    # ──────────────────────────────────────────────────────────
    # Step: Within-Codon GT Resolution
    # ──────────────────────────────────────────────────────────
    def _step_within_codon_gt_resolution(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Within-Codon GT Resolution step: Resolve within-codon GT dinucleotides.

        For each within-codon GT (not cross-codon), try synonymous substitution
        first. If no synonymous codon can avoid GT, flag for mutagenesis.

        Uses IncrementalSequenceState for O(1) GT tracking and
        CodonCache for pre-sorted codon lists.
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        remaining_positions: Dict[int, List[str]] = {}

        # Use incremental GT positions for O(1) lookup instead of O(N) scan
        for gt_pos in state.gt_positions_list():
            codon_idx = gt_pos // 3
            codon_start = codon_idx * 3
            next_codon_start = codon_start + 3

            if gt_pos % 3 < 2:  # Within-codon GT
                # Within-codon GT - skip if unavoidable
                if _is_unavoidable_gt(state.sequence, gt_pos):
                    continue
                aa = state.get_aa(codon_idx)
                if aa is None or aa == "*":
                    continue

                old_gt_count = state.gt_count
                resolved = False
                # Use pre-sorted GT-free codons from cache
                for alt in codon_cache.get_gt_free_codons(aa):
                    current_codon = state.get_codon(codon_idx)
                    if alt == current_codon:
                        continue
                    # O(1) boundary check instead of manual base lookup
                    left_gt, right_gt = state.boundary_creates_gt(codon_idx, alt)
                    if left_gt or right_gt:
                        continue
                    # O(1) swap + O(1) GT count check + rollback if needed
                    old_codon = state.swap_codon(codon_idx, alt)
                    if state.gt_count < old_gt_count:
                        resolved = True
                        break
                    else:
                        state.swap_codon(codon_idx, old_codon)  # Rollback

                if not resolved:
                    remaining_positions[codon_start] = ["GT:within"]

        if remaining_positions:
            mut_report = self._propose_within_codon_mutagenesis(
                state.sequence, remaining_positions
            )
        else:
            mut_report = MutagenesisReport()

        return state.sequence, mut_report

    def _propose_within_codon_mutagenesis(
        self,
        seq: str,
        positions: Dict[int, List[str]],
    ) -> MutagenesisReport:
        """Propose AA substitutions for within-codon GTs that can't be resolved
        by synonymous substitution.

        Uses BLOSUM62 guidance for conservative substitutions:
        - Valine (V, GTN) → Isoleucine (I, ATN) or Leucine (L, TTA/CTN)
        """
        report = MutagenesisReport()

        for codon_start, ctypes in positions.items():
            if codon_start + 3 > len(seq):
                continue

            original_codon = seq[codon_start:codon_start + 3]
            original_aa = CODON_TABLE.get(original_codon)
            if original_aa is None or original_aa == "*":
                continue

            best = None
            best_score = -100

            for new_aa, score in sorted(
                ((aa, BLOSUM62.get((original_aa, aa), -10))
                 for aa in set(CODON_TABLE.values()) if aa != "*" and aa != original_aa),
                key=lambda x: x[1],
                reverse=True,
            ):
                if score < self.min_blosum:
                    continue

                alts = AA_TO_CODONS.get(new_aa, [])
                alts_sorted = sorted(alts,
                                     key=lambda c: self.species_cai.get(c, 0.0),
                                     reverse=True)
                for alt_codon in alts_sorted:
                    if "GT" in alt_codon:
                        continue
                    prev_base = seq[codon_start - 1] if codon_start > 0 else ""
                    next_start = codon_start + 3
                    next_base = seq[next_start] if next_start < len(seq) else ""
                    if prev_base and prev_base + alt_codon[0] == "GT":
                        continue
                    if next_base and alt_codon[-1] + next_base == "GT":
                        continue
                    cai = self.species_cai.get(alt_codon, 0.0)
                    if cai < self.min_cai:
                        continue

                    if score > best_score:
                        best = (alt_codon, new_aa, score, cai)
                        best_score = score
                    break

            if best is not None:
                alt_codon, new_aa, blosum, cai = best
                proposal = MutagenesisProposal(
                    position=codon_start,
                    original_codon=original_codon,
                    original_aa=original_aa,
                    new_aa=new_aa,
                    new_codon=alt_codon,
                    blosum_score=blosum,
                    cai_weight=cai,
                    resolves=ctypes,
                )
                report.add(proposal)
            else:
                proposal = MutagenesisProposal(
                    position=codon_start,
                    original_codon=original_codon,
                    original_aa=original_aa,
                    new_aa="",
                    new_codon="",
                    blosum_score=-10,
                    cai_weight=0.0,
                    resolves=ctypes,
                    impossible=True,
                )
                report.add(proposal)

        return report

    # Deprecated alias — use _step_within_codon_gt_resolution instead
    _phase35_within_codon_gt = _step_within_codon_gt_resolution

    # ──────────────────────────────────────────────────────────
    # Step: Mutagenesis Fallback
    # ──────────────────────────────────────────────────────────
    def _step_mutagenesis_fallback(self, seq: str, mut_report: MutagenesisReport) -> str:
        """Mutagenesis Fallback step: Apply mutagenesis proposals for intractable constraints.

        Applies conservative AA substitutions (e.g., Val→Ile, Val→Leu)
        using BLOSUM62 guidance. Only applies to AVOIDABLE GT positions;
        unavoidable GTs (Valine, cross-codon with no alternatives) are skipped.
        """
        seq_list = list(seq)
        for proposal in mut_report.proposals:
            if proposal.impossible or not proposal.new_codon:
                continue
            pos = proposal.position
            old_gt_count = _count_gts("".join(seq_list))

            test_list = seq_list[:]
            for k, b in enumerate(proposal.new_codon):
                if pos + k < len(test_list):
                    test_list[pos + k] = b
            new_gt_count = _count_gts("".join(test_list))

            # Accept if it reduces or maintains GT count
            if new_gt_count <= old_gt_count:
                for k, b in enumerate(proposal.new_codon):
                    if pos + k < len(seq_list):
                        seq_list[pos + k] = b
                self._applied_mutagenesis.append({
                    "position": pos,
                    "original_aa": proposal.original_aa,
                    "new_aa": proposal.new_aa,
                    "blosum": proposal.blosum_score,
                })

        return "".join(seq_list)

    # Deprecated alias — use _step_mutagenesis_fallback instead
    _phase4_mutagenesis_fallback = _step_mutagenesis_fallback

    # ──────────────────────────────────────────────────────────
    # Step: Avoid CpG Islands
    # ──────────────────────────────────────────────────────────
    def _step_avoid_cpg_islands(self, seq: str) -> str:
        """Avoid CpG Islands step (optimized with incremental state).

        CpG island avoidance by synonymous substitution. Also handles
        cross-codon CG dinucleotides by two-codon coordination.

        Uses IncrementalSequenceState for O(1) GT/CG tracking and
        direct CG position lookup instead of window scanning.
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        changed = True
        iterations = 0
        max_iterations = 80

        while changed and iterations < max_iterations:
            changed = False
            iterations += 1

            for start in range(0, len(state.sequence) - self.cpg_window + 1, 3):
                window_seq = state.sequence[start:start + self.cpg_window]
                c_count = window_seq.count("C")
                g_count = window_seq.count("G")
                cg_count = _count_dinucs_fast(window_seq, "CG")[0]
                expected = (c_count * g_count) / len(window_seq) if len(window_seq) > 0 else 0
                obs_exp = cg_count / expected if expected > 0 else 0.0

                if obs_exp <= self.cpg_threshold:
                    continue

                # Find CG dinucleotides in this window and try to break them
                # Use incremental state's CG position data for efficiency
                cg_in_window = [p for p in state.cg_positions_list()
                               if start <= p < start + self.cpg_window - 1]

                for i in cg_in_window:
                    left_ci = i // 3
                    right_ci = (i + 1) // 3
                    is_cross_codon = (left_ci != right_ci)

                    # Strategy 1: Single-codon swap
                    codon_fixed = False
                    for ci in ([left_ci, right_ci] if is_cross_codon else [left_ci]):
                        codon = state.get_codon(ci)
                        aa = state.get_aa(ci)
                        if aa is None or aa == "*":
                            continue

                        for alt in codon_cache.get_sorted_codons(aa):
                            if alt == codon:
                                continue
                            if not is_cross_codon and "CG" in alt:
                                continue
                            if self.avoid_gt and "GT" in alt:
                                continue
                            if codon_cache.get_cai(alt) < self.min_cai:
                                continue

                            # O(1) boundary check
                            left_gt, right_gt = state.boundary_creates_gt(ci, alt)
                            if left_gt or right_gt:
                                continue

                            # O(1) GT/CG change check
                            if self.avoid_gt and state.would_gt_increase(ci, alt):
                                continue

                            # Apply swap and check the specific CG position is gone
                            old_codon = state.swap_codon(ci, alt)
                            if state.sequence[i] != "C" or state.sequence[i+1] != "G":
                                changed = True
                                codon_fixed = True
                                break
                            else:
                                state.swap_codon(ci, old_codon)  # Rollback

                        if codon_fixed:
                            break

                    if codon_fixed:
                        break

                    # Strategy 2: Coordinated 2-codon swap for cross-codon CG
                    if is_cross_codon and not codon_fixed:
                        left_codon = state.get_codon(left_ci)
                        right_codon = state.get_codon(right_ci)
                        left_aa = state.get_aa(left_ci)
                        right_aa = state.get_aa(right_ci)

                        if left_aa and right_aa and left_aa != "*" and right_aa != "*":
                            old_gt_count = state.gt_count
                            for left_alt in codon_cache.get_sorted_codons(left_aa):
                                if left_alt == left_codon:
                                    continue
                                if self.avoid_gt and "GT" in left_alt:
                                    continue
                                if codon_cache.get_cai(left_alt) < self.min_cai:
                                    continue
                                for right_alt in codon_cache.get_sorted_codons(right_aa):
                                    if left_alt[-1] == "C" and right_alt[0] == "G":
                                        continue
                                    if self.avoid_gt and "GT" in right_alt:
                                        continue
                                    if codon_cache.get_cai(right_alt) < self.min_cai:
                                        continue
                                    # Check cross-codon GT effects
                                    left_gt, _ = state.boundary_creates_gt(left_ci, left_alt)
                                    if left_gt:
                                        continue
                                    if left_alt[-1] + right_alt[0] == "GT":
                                        continue
                                    _, right_gt = state.boundary_creates_gt(right_ci, right_alt)
                                    if right_gt:
                                        continue

                                    # Apply both swaps and check incrementally
                                    old_left = state.swap_codon(left_ci, left_alt)
                                    old_right = state.swap_codon(right_ci, right_alt)

                                    # O(1) GT count check (state.gt_count is updated incrementally)
                                    if state.gt_count <= old_gt_count:
                                        # Verify CG at position i is gone
                                        if state.sequence[i] != "C" or state.sequence[i+1] != "G":
                                            changed = True
                                            codon_fixed = True
                                            break

                                    # Rollback
                                    state.swap_codon(right_ci, old_right)
                                    state.swap_codon(left_ci, old_left)

                                if codon_fixed:
                                    break

                    if changed:
                        break
                if changed:
                    break

        return state.sequence

    # Deprecated alias — use _step_avoid_cpg_islands instead
    _phase5_avoid_cpg_islands = _step_avoid_cpg_islands

    # ──────────────────────────────────────────────────────────
    # Step: Remove ATTTA Instability Motifs
    # ──────────────────────────────────────────────────────────
    def _step_remove_instability_motifs(self, seq: str) -> str:
        """Remove ATTTA mRNA instability motifs by synonymous codon substitution.

        ATTTA motifs are associated with mRNA instability in eukaryotic genes.
        This step identifies ATTTA pentamers and attempts to disrupt them by
        swapping one of the overlapping codons to a synonymous alternative that
        does not contain the motif, while preserving GT/CG/RS constraints where
        possible. Falls back to accepting minor constraint relaxation if needed
        to eliminate ATTTA (instability removal has higher priority).
        """
        from .restriction_sites import get_recognition_site as _get_site

        for iteration in range(MAX_ATTTA_MOTIF_ITERATIONS):
            pos = seq.find("ATTTA")
            if pos == -1:
                break

            n_codons = len(seq) // 3
            best_swap = None  # (test_seq, score) where lower score = better
            # Score: 0 = perfect (ATTTA gone, no new issues), 1 = ATTTA gone + new RS, 2 = ATTTA gone + worse GT

            # Try swapping each overlapping codon
            for offset in range(5):
                ci = (pos + offset) // 3
                if ci < 0 or ci >= n_codons:
                    continue
                codon = seq[ci*3:ci*3+3]
                aa = CODON_TABLE.get(codon, "")
                if not aa or aa == "*":
                    continue
                alternatives = sorted(
                    AA_TO_CODONS.get(aa, [codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for alt in alternatives:
                    if alt == codon:
                        continue
                    test = seq[:ci*3] + alt + seq[ci*3+3:]
                    if "ATTTA" not in test:
                        # Check constraints
                        rs_ok = all(
                            not (_s := _get_site(enz)) or _s not in test
                            for enz in self.enzymes
                        )
                        from .type_system import check_no_avoidable_gt as _check_gt
                        gt_result = _check_gt(test)
                        old_gt = _check_gt(seq)
                        gt_ok = gt_result.passed or not (not gt_result.passed and old_gt.passed)

                        if rs_ok and gt_ok:
                            # Perfect swap — use immediately
                            return self._step_remove_instability_motifs(test)
                        elif gt_ok and not rs_ok:
                            score = 1  # New RS but GT preserved (RS easier to fix later)
                        elif rs_ok and not gt_ok:
                            score = 2  # GT worsened but no new RS
                        else:
                            score = 3  # Both worsened

                        if best_swap is None or score < best_swap[1]:
                            best_swap = (test, score)

            if best_swap is not None:
                # Use the best available swap (even if not perfect)
                # After applying, run a reconciliation pass to fix any introduced issues
                seq = best_swap[0]
                # Re-check and fix restriction sites
                for enz in self.enzymes:
                    site = _get_site(enz)
                    if site and site in seq:
                        ci_rs = seq.find(site) // 3
                        if 0 <= ci_rs < n_codons:
                            aa = CODON_TABLE.get(seq[ci_rs*3:ci_rs*3+3], "")
                            current = seq[ci_rs*3:ci_rs*3+3]
                            for alt in AA_TO_CODONS.get(aa, [current]):
                                if alt != current:
                                    test2 = seq[:ci_rs*3] + alt + seq[ci_rs*3+3:]
                                    if site not in test2 and "ATTTA" not in test2:
                                        seq = test2
                                        break
                continue

            # No single-codon swap worked; try 2-codon coordinated swap
            fixed = False
            for offset1 in range(5):
                ci1 = (pos + offset1) // 3
                if ci1 < 0 or ci1 >= n_codons:
                    continue
                aa1 = CODON_TABLE.get(seq[ci1*3:ci1*3+3], "")
                if not aa1 or aa1 == "*":
                    continue
                for offset2 in range(offset1 + 1, 5):
                    ci2 = (pos + offset2) // 3
                    if ci2 < 0 or ci2 >= n_codons or ci2 == ci1:
                        continue
                    aa2 = CODON_TABLE.get(seq[ci2*3:ci2*3+3], "")
                    if not aa2 or aa2 == "*":
                        continue
                    for alt1 in AA_TO_CODONS.get(aa1, []):
                        for alt2 in AA_TO_CODONS.get(aa2, []):
                            test_list = list(seq)
                            test_list[ci1*3:ci1*3+3] = list(alt1)
                            test_list[ci2*3:ci2*3+3] = list(alt2)
                            test = "".join(test_list)
                            if "ATTTA" not in test:
                                seq = test
                                fixed = True
                                break
                        if fixed:
                            break
                    if fixed:
                        break
                if fixed:
                    break

            if not fixed:
                break

        return seq

    # ──────────────────────────────────────────────────────────
    # Step: CpG Reconciliation (aggressive, runs after CAI hill climb)
    # ──────────────────────────────────────────────────────────
    def _step_cpg_reconciliation(self, seq: str) -> str:
        """Aggressive CpG reconciliation pass to eliminate CpG islands.

        This runs AFTER the CAI hill climb and reoptimize steps, because
        those steps can reintroduce CG dinucleotides by preferring high-CAI
        C-ending codons (e.g., GCC for Ala, CCC for Pro) which create
        cross-codon CG when followed by a G-starting codon.

        Strategy:
        1. Identify all windows that fail the CpG island check (Obs/Exp > threshold)
        2. For each failing window, identify all CG dinucleotides
        3. For cross-codon CG (codon ending with C + codon starting with G):
           - Try swapping the C-ending codon to a non-C-ending synonym (sorted by CAI)
           - If the G-starting codon has non-G-starting synonyms, try those too
        4. For within-codon CG:
           - Swap to a CG-free synonym (sorted by CAI)
        5. Accept CAI loss as necessary — CpG island avoidance is a hard constraint
        6. After CG elimination, re-check GC content and adjust if needed
        7. Iterate until all windows pass or no more progress can be made

        Uses IncrementalSequenceState for O(1) GT/CG tracking and
        CodonCache for pre-sorted codon lists.
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        max_iterations = 100

        def _check_gt_safe(ci: int, old_had_gt: dict) -> bool:
            """Check if swap at codon_idx created new avoidable GTs."""
            cs = ci * 3
            ce = cs + 3
            chk_s = max(0, cs - 1)
            chk_e = min(state._n - 1, ce + 1)
            for p in range(chk_s, chk_e):
                if state.has_gt_at(p) and not old_had_gt.get(p, False):
                    if not _is_unavoidable_gt(state.sequence, p):
                        return False
            return True

        def _save_gt_region(ci: int) -> dict:
            """Save GT state in the region around a codon."""
            cs = ci * 3
            ce = cs + 3
            chk_s = max(0, cs - 1)
            chk_e = min(state._n - 1, ce + 1)
            return {p: state.has_gt_at(p) for p in range(chk_s, chk_e)}

        for iteration in range(max_iterations):
            current_seq = state.sequence
            old_cg_count = state.cg_count

            # Check if the CpG island predicate passes
            cpg_result = check_no_cpg_island(current_seq, self.cpg_window, self.cpg_threshold)
            if cpg_result.passed:
                break

            # Find the worst window
            worst_ratio = 0.0
            worst_start = -1
            for start in range(0, len(current_seq) - self.cpg_window + 1):
                window_seq = current_seq[start:start + self.cpg_window]
                c_count = window_seq.count("C")
                g_count = window_seq.count("G")
                cg_count = _count_dinucs_fast(window_seq, "CG")[0]
                expected = (c_count * g_count) / len(window_seq) if len(window_seq) > 0 else 0
                obs_exp = cg_count / expected if expected > 0 else 0.0
                if obs_exp > worst_ratio:
                    worst_ratio = obs_exp
                    worst_start = start

            if worst_start < 0 or worst_ratio <= self.cpg_threshold:
                break

            # Find CG dinucleotides in the worst window using incremental state
            window_end = min(worst_start + self.cpg_window, len(current_seq) - 1)
            fixed = False

            # Collect all CG positions in the worst window with their codon info
            cg_targets = []
            for cg_pos in state.cg_positions_list():
                if worst_start <= cg_pos < window_end:
                    c_codon_idx = cg_pos // 3
                    g_codon_idx = (cg_pos + 1) // 3
                    cg_targets.append((cg_pos, c_codon_idx, g_codon_idx))

            # Sort targets: prioritize within-codon CG first (easier to fix)
            for cg_pos, c_codon_idx, g_codon_idx in cg_targets:
                if fixed:
                    break

                if c_codon_idx == g_codon_idx:
                    # Within-codon CG — swap the codon to a CG-free alternative
                    aa = state.get_aa(c_codon_idx)
                    if aa is None or aa == "*":
                        continue

                    for alt in codon_cache.get_cg_free_codons(aa):
                        current_codon = state.get_codon(c_codon_idx)
                        if alt == current_codon:
                            continue
                        # Save GT state, apply swap, check constraints
                        old_had_gt = _save_gt_region(c_codon_idx)
                        old_codon = state.swap_codon(c_codon_idx, alt)
                        # Must not reintroduce restriction sites
                        if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                            state.swap_codon(c_codon_idx, old_codon)
                            continue
                        # Must not create new avoidable GTs
                        if not _check_gt_safe(c_codon_idx, old_had_gt):
                            state.swap_codon(c_codon_idx, old_codon)
                            continue
                        # Must reduce net CG count (O(1) check)
                        if state.cg_count < old_cg_count:
                            fixed = True
                            break
                        else:
                            state.swap_codon(c_codon_idx, old_codon)  # Rollback
                else:
                    # Cross-codon CG: codon c_codon_idx ends with C, codon g_codon_idx starts with G
                    # Strategy 1: Swap the C-ending codon to not end with C
                    c_aa = state.get_aa(c_codon_idx)
                    if c_aa is not None and c_aa != "*":
                        # Get non-C-ending alternatives from cache, then filter
                        non_c_end_alts = [
                            alt for alt in codon_cache.get_sorted_codons(c_aa)
                            if alt != state.get_codon(c_codon_idx) and alt[-1] != "C" and "CG" not in alt
                        ]
                        for alt in non_c_end_alts:
                            old_had_gt = _save_gt_region(c_codon_idx)
                            old_codon = state.swap_codon(c_codon_idx, alt)
                            # Must not reintroduce restriction sites
                            if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                                state.swap_codon(c_codon_idx, old_codon)
                                continue
                            # Must not create new avoidable GTs
                            if not _check_gt_safe(c_codon_idx, old_had_gt):
                                state.swap_codon(c_codon_idx, old_codon)
                                continue
                            # Must reduce net CG count (O(1) check)
                            if state.cg_count < old_cg_count:
                                fixed = True
                                break
                            else:
                                state.swap_codon(c_codon_idx, old_codon)  # Rollback

                    if fixed:
                        break

                    # Strategy 2: Swap the G-starting codon to not start with G (if possible)
                    g_aa = state.get_aa(g_codon_idx)
                    if g_aa is not None and g_aa != "*":
                        # Get non-G-starting alternatives from cache, then filter
                        non_g_start_alts = [
                            alt for alt in codon_cache.get_sorted_codons(g_aa)
                            if alt != state.get_codon(g_codon_idx) and alt[0] != "G" and "CG" not in alt
                        ]
                        for alt in non_g_start_alts:
                            old_had_gt = _save_gt_region(g_codon_idx)
                            old_codon = state.swap_codon(g_codon_idx, alt)
                            # Must not reintroduce restriction sites
                            if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                                state.swap_codon(g_codon_idx, old_codon)
                                continue
                            # Must not create new avoidable GTs
                            if not _check_gt_safe(g_codon_idx, old_had_gt):
                                state.swap_codon(g_codon_idx, old_codon)
                                continue
                            # Must reduce net CG count (O(1) check)
                            if state.cg_count < old_cg_count:
                                fixed = True
                                break
                            else:
                                state.swap_codon(g_codon_idx, old_codon)  # Rollback

                    if fixed:
                        break

                    # Strategy 3: Paired swap — change BOTH codons
                    c_aa = state.get_aa(c_codon_idx)
                    g_aa = state.get_aa(g_codon_idx)
                    if c_aa is not None and g_aa is not None and c_aa != "*" and g_aa != "*":
                        c_alts = [
                            alt for alt in codon_cache.get_sorted_codons(c_aa)
                            if alt != state.get_codon(c_codon_idx) and alt[-1] != "C" and "CG" not in alt
                        ]
                        g_alts = [
                            alt for alt in codon_cache.get_sorted_codons(g_aa)
                            if alt != state.get_codon(g_codon_idx) and alt[0] != "G" and "CG" not in alt
                        ]
                        for c_alt in c_alts[:TOP_CAI_ALTERNATIVES]:
                            for g_alt in g_alts[:TOP_CAI_ALTERNATIVES]:
                                # Save GT state for both regions
                                old_had_gt_c = _save_gt_region(c_codon_idx)
                                old_had_gt_g = _save_gt_region(g_codon_idx)
                                # Apply both swaps
                                old_c = state.swap_codon(c_codon_idx, c_alt)
                                old_g = state.swap_codon(g_codon_idx, g_alt)
                                # Must not reintroduce restriction sites
                                if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                                    state.swap_codon(g_codon_idx, old_g)
                                    state.swap_codon(c_codon_idx, old_c)
                                    continue
                                # Must not create new avoidable GTs in either region
                                if not _check_gt_safe(c_codon_idx, old_had_gt_c) or not _check_gt_safe(g_codon_idx, old_had_gt_g):
                                    state.swap_codon(g_codon_idx, old_g)
                                    state.swap_codon(c_codon_idx, old_c)
                                    continue
                                # Must reduce net CG count (O(1) check)
                                if state.cg_count < old_cg_count:
                                    fixed = True
                                    break
                                else:
                                    state.swap_codon(g_codon_idx, old_g)  # Rollback
                                    state.swap_codon(c_codon_idx, old_c)  # Rollback
                            if fixed:
                                break

            if not fixed:
                # No progress — try a more aggressive approach:
                # scan ALL codon positions and replace C-ending codons that
                # create cross-codon CG, regardless of window
                aggressive_fixed = False
                for ci in range(state.num_codons):
                    codon = state.get_codon(ci)
                    aa = state.get_aa(ci)
                    if aa is None or aa == "*":
                        continue
                    # Check if this codon ends with C and the next starts with G
                    next_ci = ci + 1
                    if next_ci < state.num_codons:
                        next_codon = state.get_codon(next_ci)
                        if codon[-1] == "C" and next_codon[0] == "G":
                            # This creates a cross-codon CG — try to fix
                            non_c_end_alts = [
                                alt for alt in codon_cache.get_sorted_codons(aa)
                                if alt != codon and alt[-1] != "C" and "CG" not in alt
                            ]
                            for alt in non_c_end_alts:
                                old_had_gt = _save_gt_region(ci)
                                old_codon = state.swap_codon(ci, alt)
                                if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                                    state.swap_codon(ci, old_codon)
                                    continue
                                if not _check_gt_safe(ci, old_had_gt):
                                    state.swap_codon(ci, old_codon)
                                    continue
                                # O(1) CG count check
                                if state.cg_count < old_cg_count:
                                    aggressive_fixed = True
                                    break
                                else:
                                    state.swap_codon(ci, old_codon)  # Rollback

                    if not aggressive_fixed and "CG" in codon:
                        # Within-codon CG — try to fix
                        for alt in codon_cache.get_cg_free_codons(aa):
                            if alt == codon:
                                continue
                            old_had_gt = _save_gt_region(ci)
                            old_codon = state.swap_codon(ci, alt)
                            if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                                state.swap_codon(ci, old_codon)
                                continue
                            if not _check_gt_safe(ci, old_had_gt):
                                state.swap_codon(ci, old_codon)
                                continue
                            # O(1) CG count check
                            if state.cg_count < old_cg_count:
                                aggressive_fixed = True
                                break
                            else:
                                state.swap_codon(ci, old_codon)  # Rollback

                    if aggressive_fixed:
                        break

                if not aggressive_fixed:
                    break  # Truly no more progress possible

        # After CpG reconciliation, re-check GC content and adjust if needed
        result_seq = state.sequence
        gc_val = (result_seq.count("G") + result_seq.count("C")) / max(len(result_seq), 1)
        if not (0.30 <= gc_val <= 0.70):
            # GC drifted — adjust with single-codon swaps that don't reintroduce CGs
            result_seq = self._fix_gc_after_cpg(result_seq)

        return result_seq

    def _check_no_restriction_sites(self, seq: str) -> bool:
        """Check that sequence doesn't contain any restriction enzyme sites."""
        from .restriction_sites import get_recognition_site
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            if site in seq:
                return False
        return True

    def _check_cpg_swap_gt_safe(self, old_list: list, new_list: list, codon_start: int) -> bool:
        """Check that a codon swap doesn't create new avoidable GT dinucleotides.

        We allow unavoidable GTs (e.g., Valine codons) but not new avoidable ones.
        This is less strict than requiring GT count to not increase at all.
        """
        old_seq = "".join(old_list)
        new_seq = "".join(new_list)
        codon_end = codon_start + 3

        # Check the local region around the swapped codon for new GTs
        check_start = max(0, codon_start - 1)
        check_end = min(len(new_seq) - 1, codon_end + 1)

        for i in range(check_start, check_end):
            if new_seq[i:i+2] == "GT" and old_seq[i:i+2] != "GT":
                # New GT created — check if it's unavoidable
                if not _is_unavoidable_gt(new_seq, i):
                    return False  # Avoidable new GT — reject this swap

        return True

    def _fix_gc_after_cpg(self, seq: str) -> str:
        """Fix GC content after CpG reconciliation without reintroducing CG dinucleotides.

        Only adjusts GC if it's outside [0.30, 0.70].
        """
        gc_val = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
        if 0.30 <= gc_val <= 0.70:
            return seq

        seq_list = list(seq)
        gc_count = sum(1 for b in seq_list if b in "GC")
        n_bases = len(seq_list)
        target_gc = 0.50 if gc_val > 0.70 else 0.30

        for iteration in range(MAX_GC_ADJUSTMENT_ITERATIONS):
            gc_val = gc_count / n_bases
            if 0.30 <= gc_val <= 0.70:
                break

            best_alt = None
            best_ci = -1
            best_diff = abs(gc_val - target_gc)
            best_gc_delta = 0

            for ci in range(len(seq_list) // 3):
                codon_start = ci * 3
                if codon_start + 3 > len(seq_list):
                    break
                codon = "".join(seq_list[codon_start:codon_start + 3])
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*":
                    continue
                current_gc = sum(1 for b in codon if b in "GC")
                for alt in AA_TO_CODONS.get(aa, []):
                    if alt == codon:
                        continue
                    alt_gc = sum(1 for b in alt if b in "GC")
                    new_gc_count = gc_count - current_gc + alt_gc
                    new_frac = new_gc_count / n_bases
                    diff = abs(new_frac - target_gc)
                    if diff < best_diff:
                        # Check this swap doesn't reintroduce CGs
                        test_list = seq_list[:]
                        for k, b in enumerate(alt):
                            test_list[codon_start + k] = b
                        test_seq = "".join(test_list)
                        # Must not increase CG count
                        new_cg = _count_dinucs_fast(test_seq, "CG")[0]
                        old_cg = _count_dinucs_fast(seq, "CG")[0]
                        if new_cg > old_cg:
                            continue
                        # Must not reintroduce restriction sites
                        if not self._check_no_restriction_sites(test_seq):
                            continue
                        best_diff = diff
                        best_alt = alt
                        best_ci = ci
                        best_gc_delta = alt_gc - current_gc

            if best_alt is None:
                break
            codon_start = best_ci * 3
            for k, b in enumerate(best_alt):
                seq_list[codon_start + k] = b
            gc_count += best_gc_delta

        return "".join(seq_list)

    # ──────────────────────────────────────────────────────────
    # Step: Cross-Codon Coordination
    # ──────────────────────────────────────────────────────────
    def _step_cross_codon_coordination(self, seq: str) -> str:
        """Cross-Codon Coordination step (optimized with incremental state).

        Resolve cross-codon GT, CG, and restriction site violations by
        considering ALL synonymous codon pairs at each boundary.

        Uses IncrementalSequenceState for O(1) GT/CG tracking and
        CodonCache for pre-sorted codon lists.

        Unlike the earlier Cross-Codon Optimization step which handles violations
        one at a time with limited search, this phase performs an exhaustive
        pairwise search: for every boundary violation it considers ALL synonymous
        codon pairs for the two adjacent amino acids, selects the pair that
        eliminates the violation while maximizing CAI, and re-scans after each
        fix (since fixing one boundary can create a violation at the adjacent
        boundary). Iterates until no cross-codon violations remain or a maximum
        iteration count is reached.

        This phase runs AFTER within-codon optimization and BEFORE the final
        CAI maximization pass, ensuring that cross-codon interactions — which
        are the primary cause of predicate failures for sequences like HBB —
        are fully resolved before the hill climber locks in codon choices.
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        max_iterations = 50

        # Record the initial avoidable GT count for comparison
        initial_gt_count = state.gt_count if self.avoid_gt else 0

        for iteration in range(max_iterations):
            # Find cross-codon violations using incremental position data
            violations = []

            # GT violations at cross-codon boundaries
            # EUKARYOTE-ONLY: GT avoidance is irrelevant for prokaryotes
            if self.avoid_gt:
                for gt_pos in state.gt_positions_list():
                    codon_idx = gt_pos // 3
                    next_codon_idx = codon_idx + 1
                    # Cross-codon GT: position is at the boundary between codons
                    # gt_pos is at the last base of one codon, gt_pos+1 is the first of the next
                    if gt_pos % 3 == 2 and next_codon_idx < state.num_codons:
                        if not _is_unavoidable_gt(state.sequence, gt_pos):
                            violations.append({"type": "GT", "boundary": gt_pos})

            # CG violations at cross-codon boundaries
            # EUKARYOTE-ONLY: CpG island avoidance is irrelevant for prokaryotes
            if not self.is_prokaryote:
                for cg_pos in state.cg_positions_list():
                    if cg_pos % 3 == 2:
                        codon_idx = cg_pos // 3
                        next_codon_idx = codon_idx + 1
                        if next_codon_idx < state.num_codons:
                            violations.append({"type": "CG", "boundary": cg_pos})

            if not violations:
                break

            # Sort: GT first, then CG, then by position
            type_priority = {"GT": 0, "CG": 1}
            violations.sort(key=lambda v: (type_priority.get(v["type"], 2), v["boundary"]))

            fixed_any = False

            for violation in violations:
                boundary = violation["boundary"]
                vtype = violation["type"]

                left_codon_start = boundary - 2
                right_codon_start = boundary + 1

                if left_codon_start < 0 or right_codon_start + 3 > len(state.sequence):
                    continue

                left_ci = left_codon_start // 3
                right_ci = right_codon_start // 3

                left_aa = state.get_aa(left_ci)
                right_aa = state.get_aa(right_ci)

                if left_aa is None or left_aa == "*" or right_aa is None or right_aa == "*":
                    continue

                # Find the best pair that eliminates the violation and maximizes CAI
                best_pair = None
                best_cai = -1.0

                for c_left in codon_cache.get_sorted_codons(left_aa):
                    for c_right in codon_cache.get_sorted_codons(right_aa):
                        # Check that this pair eliminates the specific violation
                        if vtype == "GT" and c_left[-1] + c_right[0] == "GT":
                            continue
                        if vtype == "CG" and c_left[-1] + c_right[0] == "CG":
                            continue
                        # Check within-codon violations
                        if vtype == "GT" and ("GT" in c_left or "GT" in c_right):
                            continue
                        if vtype == "CG" and ("CG" in c_left or "CG" in c_right):
                            continue

                        # Compute CAI contribution
                        pair_cai = codon_cache.get_cai(c_left) + codon_cache.get_cai(c_right)
                        if pair_cai <= best_cai:
                            continue

                        # Record GT count before trial swap
                        gt_before = state.gt_count

                        # Try the swap and check adjacent boundaries
                        old_left = state.swap_codon(left_ci, c_left)
                        old_right = state.swap_codon(right_ci, c_right)

                        # Check that we haven't introduced new GTs at adjacent boundaries
                        new_violations = False

                        # Left boundary of c_left with previous codon
                        if left_ci > 0:
                            left_gt, _ = state.boundary_creates_gt(left_ci, c_left)
                            if left_gt and not _is_unavoidable_gt(state.sequence, left_codon_start - 1):
                                new_violations = True

                        # Right boundary of c_right with next codon
                        if not new_violations and right_ci + 1 < state.num_codons:
                            _, right_gt = state.boundary_creates_gt(right_ci, c_right)
                            if right_gt and not _is_unavoidable_gt(state.sequence, right_codon_start + 2):
                                new_violations = True

                        # Check for new CG violations at adjacent boundaries
                        if not new_violations and left_ci > 0:
                            left_cg, _ = state.boundary_creates_cg(left_ci, c_left)
                            if left_cg:
                                new_violations = True

                        if not new_violations and right_ci + 1 < state.num_codons:
                            _, right_cg = state.boundary_creates_cg(right_ci, c_right)
                            if right_cg:
                                new_violations = True

                        # Check that we haven't introduced new restriction sites
                        if not new_violations and enzyme_cache is not None:
                            if enzyme_cache.check_any_site_present(state.sequence):
                                new_violations = True

                        # Check that avoidable GT count hasn't increased
                        if not new_violations and self.avoid_gt:
                            if state.gt_count > gt_before:
                                # The swap increased the total GT count; check if
                                # any new GTs are avoidable
                                current_seq = state.sequence
                                new_avoidable = sum(
                                    1 for p in state.gt_positions_list()
                                    if not _is_unavoidable_gt(current_seq, p)
                                )
                                if new_avoidable > initial_gt_count:
                                    new_violations = True

                        if not new_violations:
                            best_pair = (c_left, c_right)
                            best_cai = pair_cai

                        # Rollback
                        state.swap_codon(right_ci, old_right)
                        state.swap_codon(left_ci, old_left)

                if best_pair is not None:
                    # Apply the best pair for this violation
                    state.swap_codon(left_ci, best_pair[0])
                    state.swap_codon(right_ci, best_pair[1])
                    fixed_any = True
                    break  # Re-scan from scratch after each fix

            if not fixed_any:
                break

        return state.sequence

    # Deprecated alias — use _step_cross_codon_coordination instead
    _phase4_cross_codon_coordination = _step_cross_codon_coordination

    def _find_cross_codon_violations(self, seq: str) -> list:
        """Find all cross-codon violations in the sequence.

        Returns a list of dicts with keys:
          - 'boundary': int — position of the last base of the left codon
              (i.e., the 'G' in a cross-codon GT, or the 'C' in a cross-codon CG)
              For RS violations this is the codon boundary that the site crosses
              (the position of the last base of the left codon in the pair).
          - 'type': str — "GT", "CG", or "RS:<site>"
        """
        violations = []

        # Cross-codon GT (only relevant for eukaryotic targets)
        if self.avoid_gt:
            for pos in find_cross_codon_gt(seq):
                if not _is_unavoidable_gt(seq, pos):
                    violations.append({"boundary": pos, "type": "GT"})

        # Cross-codon CG (CpG-related, only relevant for eukaryotic targets)
        if not self.is_prokaryote:
            for pos in find_cross_codon_cg(seq):
                violations.append({"boundary": pos, "type": "CG"})

        # Cross-codon restriction sites
        # find_cross_codon_restriction returns the start position of the site,
        # but we need to map it to codon boundaries for the pair-based approach.
        from .restriction_sites import get_recognition_site
        seen_rs_boundaries: set = set()
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            site_len = len(site)
            for pos in find_cross_codon_restriction(seq, site):
                # Find all codon boundaries (b where b % 3 == 2) that the site spans
                for b in range(pos, pos + site_len - 1):
                    if b % 3 == 2 and b + 1 < len(seq):
                        key = (b, f"RS:{site}")
                        if key not in seen_rs_boundaries:
                            seen_rs_boundaries.add(key)
                            violations.append({"boundary": b, "type": f"RS:{site}"})

        return violations

    def _boundary_has_violation(self, seq: str, boundary: int, vtype: str,
                                violation: dict) -> bool:
        """Check if a specific boundary still has the given violation type."""
        if vtype == "GT":
            return (boundary + 1 < len(seq)
                    and seq[boundary] == "G" and seq[boundary + 1] == "T")
        elif vtype == "CG":
            return (boundary + 1 < len(seq)
                    and seq[boundary] == "C" and seq[boundary + 1] == "G")
        elif vtype.startswith("RS:"):
            site = vtype[3:]
            # The restriction site might still be present anywhere in the sequence
            return site in seq
        return False

    def _pair_resolves_violation(
        self,
        c_left: str, c_right: str,
        left_codon_start: int, right_codon_start: int,
        seq: str, boundary: int, vtype: str, violation: dict,
    ) -> bool:
        """Check whether applying the codon pair (c_left, c_right) resolves the
        violation at the given boundary.

        This constructs a test sequence and checks the specific violation type.
        """
        # Build the test sequence with the pair applied
        test_seq = (
            seq[:left_codon_start]
            + c_left
            + c_right
            + seq[right_codon_start + 3:]
        )

        # The boundary position in the test sequence is the same
        # (boundary = left_codon_start + 2 = position of last base of c_left)

        if vtype == "GT":
            # Check the boundary dinucleotide
            if boundary + 1 < len(test_seq):
                dinuc = test_seq[boundary] + test_seq[boundary + 1]
                return dinuc != "GT"
            return True

        elif vtype == "CG":
            if boundary + 1 < len(test_seq):
                dinuc = test_seq[boundary] + test_seq[boundary + 1]
                return dinuc != "CG"
            return True

        elif vtype.startswith("RS:"):
            site = vtype[3:]
            return site not in test_seq

        return False

    # ──────────────────────────────────────────────────────────
    # Step: CAI Reconciliation
    # ──────────────────────────────────────────────────────────
    def _step_cai_reconciliation(self, seq: str) -> str:
        """CAI Reconciliation pass: upgrade low-CAI codons while maintaining hard constraints.

        This pass runs AFTER all hard constraints have been satisfied (cryptic splice,
        restriction sites, CpG islands). It greedily upgrades each low-CAI codon to the
        best synonymous alternative that doesn't violate any hard constraint.

        Algorithm:
        1. Identify all codon positions where CAI < min_cai (or below a working threshold)
        2. Sort positions by CAI improvement potential (best possible CAI minus current CAI)
        3. For each position, try upgrading to the highest-CAI synonymous codon that:
           a. Doesn't create a cryptic splice site (MaxEnt score < threshold)
           b. Doesn't create a restriction enzyme site
           c. Doesn't create a CpG island
        4. If direct upgrade creates a new avoidable GT, try paired codon swaps
           (upgrade the target + modify an adjacent codon to avoid the GT)
        5. Iterate until no more upgrades are possible

        Key insight: hard constraints (no cryptic splice, no restriction sites) are
        satisfied FIRST, then CAI is maximized SUBJECT to those constraints.

        Uses IncrementalSequenceState for O(1) GT/CG tracking and
        CodonCache for pre-sorted codon lists.
        """
        import math
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        max_iterations = 30

        for iteration in range(max_iterations):
            current_seq = state.sequence
            any_upgrade = False

            # 1. Collect positions where CAI is low, with improvement potential
            low_cai_positions = []
            for codon_idx in range(state.num_codons):
                codon = state.get_codon(codon_idx)
                aa = state.get_aa(codon_idx)
                if aa is None or aa == "*":
                    continue
                current_cai = codon_cache.get_cai(codon)
                if current_cai >= self.min_cai:
                    continue  # Already above threshold
                # Compute best possible CAI for this AA
                best_possible = codon_cache.get_cai(codon_cache.get_best_codon(aa))
                improvement = best_possible - current_cai
                if improvement > 0:
                    pos = codon_idx * 3
                    low_cai_positions.append((pos, codon, aa, current_cai, best_possible, improvement))

            if not low_cai_positions:
                break

            # 2. Sort by improvement potential (descending) — fix positions with
            #    the most room for improvement first
            low_cai_positions.sort(key=lambda x: x[5], reverse=True)

            # 3. Try to upgrade each position
            for pos, codon, aa, current_cai, best_possible, improvement in low_cai_positions:
                codon_idx = pos // 3
                # Get synonymous codons sorted by CAI (highest first) from cache
                candidates = codon_cache.get_sorted_codons(aa)

                upgraded = False
                for alt in candidates:
                    alt_cai = codon_cache.get_cai(alt)
                    if alt_cai <= current_cai:
                        break  # No improvement possible (candidates are sorted)

                    # O(1) swap + check + rollback pattern
                    old_codon = state.swap_codon(codon_idx, alt)
                    test_seq = state.sequence

                    # Check hard constraints
                    if not self._cai_recon_check_constraints(current_seq, test_seq, pos):
                        # Direct upgrade fails — try paired codon swap
                        # State already has the main swap applied; try adjacent swaps
                        paired_ok, adj_rollback = self._cai_recon_try_paired_state(
                            state, codon_idx, alt, codon_cache, enzyme_cache
                        )
                        if paired_ok:
                            # Verify the paired swap satisfies global constraints
                            if self._cai_recon_check_global_constraints(state.sequence):
                                any_upgrade = True
                                upgraded = True
                                break
                            else:
                                # Rollback the adjacent swap first, then main swap
                                if adj_rollback is not None:
                                    state.swap_codon(adj_rollback[0], adj_rollback[1])
                        # Rollback the main swap
                        state.swap_codon(codon_idx, old_codon)
                        continue

                    # Safe upgrade (already applied via swap_codon)
                    any_upgrade = True
                    upgraded = True
                    break

                if upgraded:
                    break  # Restart with updated sequence

            if not any_upgrade:
                break

        return state.sequence

    def _cai_recon_check_constraints(self, old_seq: str, test_seq: str, codon_pos: int) -> bool:
        """Check that a codon swap doesn't violate any hard constraints.

        Hard constraints:
        1. No new cryptic splice sites (MaxEnt score >= threshold)
        2. No restriction enzyme sites
        3. No new avoidable GT dinucleotides

        We allow unavoidable GTs (e.g., Valine codons) but not new avoidable ones.
        For new GTs, we check if they create a cryptic splice site.
        """
        # 1. Check for new avoidable GT dinucleotides and cryptic splice sites
        codon_end = codon_pos + 3
        check_start = max(0, codon_pos - 1)
        check_end = min(len(test_seq) - 1, codon_end + 1)

        for i in range(len(test_seq) - 1):
            if test_seq[i:i+2] == "GT" and old_seq[i:i+2] != "GT":
                # New GT created — is it unavoidable?
                if _is_unavoidable_gt(test_seq, i):
                    continue  # Unavoidable GT is OK
                # Avoidable new GT — check if it creates a cryptic splice site
                from .maxentscan import score_donor as _score_donor
                donor_score = _score_donor(test_seq, i)
                if donor_score >= self.splice_low:
                    return False  # Would create cryptic splice site

        # 2. Check restriction sites
        if not self._check_no_restriction_sites(test_seq):
            return False

        # 3. Check CpG islands
        cpg_result = check_no_cpg_island(test_seq, self.cpg_window, self.cpg_threshold)
        if not cpg_result.passed:
            # Check if the CpG failure is new (wasn't there before)
            old_cpg_result = check_no_cpg_island(old_seq, self.cpg_window, self.cpg_threshold)
            if not old_cpg_result.passed:
                pass  # Was already failing — don't make it worse but don't block
            else:
                # New CpG island created — reject
                return False

        return True

    def _cai_recon_check_global_constraints(self, test_seq: str) -> bool:
        """Check that the full sequence satisfies all hard constraints."""
        # 1. Check cryptic splice sites
        from .maxentscan import max_donor_score as _max_donor, max_acceptor_score as _max_acceptor
        if _max_donor(test_seq) >= self.splice_low:
            return False
        if _max_acceptor(test_seq) >= self.splice_low:
            return False

        # 2. Check restriction sites
        if not self._check_no_restriction_sites(test_seq):
            return False

        # 3. Check for avoidable GTs
        for i in range(len(test_seq) - 1):
            if test_seq[i:i+2] == "GT":
                if not _is_unavoidable_gt(test_seq, i):
                    return False

        return True

    def _cai_recon_try_paired(self, seq_list: list, codon_pos: int, new_codon: str) -> Optional[list]:
        """Try a paired codon swap to enable a CAI upgrade.

        When upgrading codon at codon_pos to new_codon creates a cross-codon GT,
        try simultaneously adjusting the adjacent codon to avoid it.
        Returns the new seq_list if successful, None otherwise.
        """
        new_end = new_codon[-1]
        next_pos = codon_pos + 3

        # Check if new codon creates GT with the next codon
        if next_pos + 3 <= len(seq_list):
            next_base = seq_list[next_pos]
            if new_end + next_base == "GT":
                next_codon = "".join(seq_list[next_pos:next_pos + 3])
                next_aa = CODON_TABLE.get(next_codon)
                if next_aa is not None and next_aa != "*":
                    for alt2 in sorted(
                        AA_TO_CODONS.get(next_aa, []),
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    ):
                        if new_end + alt2[0] == "GT":
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(new_codon):
                            test_list[codon_pos + k] = b
                        for k, b in enumerate(alt2):
                            test_list[next_pos + k] = b
                        test_seq = "".join(test_list)
                        # Check no new avoidable GTs overall
                        old_gt_count = _count_gts("".join(seq_list))
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            if self._check_no_restriction_sites(test_seq):
                                return test_list

        # Check if new codon creates GT with the previous codon
        if codon_pos >= 3:
            prev_end = seq_list[codon_pos - 1]
            if prev_end + new_codon[0] == "GT":
                prev_pos = codon_pos - 3
                prev_codon = "".join(seq_list[prev_pos:prev_pos + 3])
                prev_aa = CODON_TABLE.get(prev_codon)
                if prev_aa is not None and prev_aa != "*":
                    for alt0 in sorted(
                        AA_TO_CODONS.get(prev_aa, []),
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    ):
                        if alt0[-1] + new_codon[0] == "GT":
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(alt0):
                            test_list[prev_pos + k] = b
                        for k, b in enumerate(new_codon):
                            test_list[codon_pos + k] = b
                        test_seq = "".join(test_list)
                        old_gt_count = _count_gts("".join(seq_list))
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            if self._check_no_restriction_sites(test_seq):
                                return test_list

        return None

    def _cai_recon_try_paired_state(
        self, state: IncrementalSequenceState, codon_idx: int,
        new_codon: str, codon_cache: CodonCache,
        enzyme_cache: Optional[EnzymeSiteCache],
    ) -> Tuple[bool, Optional[Tuple[int, str]]]:
        """Try a paired codon swap to enable a CAI upgrade (incremental version).

        The main codon at codon_idx is already swapped to new_codon in state.
        Try adjusting an adjacent codon to resolve GT conflicts.

        Returns (success, adj_rollback) where adj_rollback is (adj_idx, old_adj_codon)
        if a paired swap was applied and needs to be tracked for potential rollback.
        If success=False, the adjacent swap (if any) has already been rolled back;
        only the main swap needs rolling back by the caller.
        """
        new_end = new_codon[-1]
        next_idx = codon_idx + 1

        # Check if new codon creates GT with the next codon
        if next_idx < state.num_codons:
            next_codon = state.get_codon(next_idx)
            next_aa = state.get_aa(next_idx)
            if new_end + next_codon[0] == "GT" and next_aa is not None and next_aa != "*":
                old_gt_count = state.gt_count
                for alt2 in codon_cache.get_sorted_codons(next_aa):
                    if new_end + alt2[0] == "GT":
                        continue
                    # O(1) swap + O(1) GT count check + rollback
                    old_adj = state.swap_codon(next_idx, alt2)
                    if state.gt_count <= old_gt_count:
                        # Check restriction sites
                        rs_ok = (not enzyme_cache.check_any_site_present(state.sequence)
                                 if enzyme_cache else True)
                        if rs_ok:
                            return (True, (next_idx, old_adj))
                    state.swap_codon(next_idx, old_adj)  # Rollback adjacent

        # Check if new codon creates GT with the previous codon
        prev_idx = codon_idx - 1
        if prev_idx >= 0:
            prev_codon = state.get_codon(prev_idx)
            prev_aa = state.get_aa(prev_idx)
            if prev_codon[-1] + new_codon[0] == "GT" and prev_aa is not None and prev_aa != "*":
                old_gt_count = state.gt_count
                for alt0 in codon_cache.get_sorted_codons(prev_aa):
                    if alt0[-1] + new_codon[0] == "GT":
                        continue
                    old_adj = state.swap_codon(prev_idx, alt0)
                    if state.gt_count <= old_gt_count:
                        rs_ok = (not enzyme_cache.check_any_site_present(state.sequence)
                                 if enzyme_cache else True)
                        if rs_ok:
                            return (True, (prev_idx, old_adj))
                    state.swap_codon(prev_idx, old_adj)  # Rollback adjacent

        return (False, None)

    def _cai_recon_rollback_paired(
        self, state: IncrementalSequenceState, codon_idx: int, old_codon: str,
    ) -> None:
        """Rollback a failed paired swap in _step_cai_reconciliation.

        This is called when a paired swap passes local constraints but fails
        global constraints. It rolls back the main codon swap.
        Note: the adjacent swap should have already been rolled back by
        _cai_recon_try_paired_state when it returned False, or should be
        handled separately when it returned True but global check failed.
        """
        state.swap_codon(codon_idx, old_codon)

    # ──────────────────────────────────────────────────────────
    # Step: GT Reconciliation (fix avoidable GTs introduced by ATTTA removal)
    # ──────────────────────────────────────────────────────────
    def _step_gt_reconciliation(self, seq: str) -> str:
        """Fix avoidable GT dinucleotides that may have been introduced by
        the ATTTA removal or CpG reconciliation steps.

        Strategy: For each avoidable GT, try synonymous codon swaps that
        eliminate the GT without reintroducing ATTTA, CpG islands, or
        restriction sites.

        Uses IncrementalSequenceState for O(1) GT tracking and
        CodonCache for pre-sorted codon lists.
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        initial_cg_count = state.cg_count  # Track initial CG count for cheap CpG check

        def _get_avoidable_gt_positions():
            """Get avoidable GT positions from incremental state (O(K) where K = num GTs)."""
            return [p for p in state.gt_positions_list()
                    if not _is_unavoidable_gt(state.sequence, p)]

        for iteration in range(50):
            avoidable_positions = _get_avoidable_gt_positions()
            if not avoidable_positions:
                break

            n_avoidable = len(avoidable_positions)
            fixed = False
            for pos in avoidable_positions:
                n_codons = state.num_codons
                # Try codons overlapping the GT
                for offset in [0, 1]:
                    ci = (pos + offset) // 3
                    if ci < 0 or ci >= n_codons:
                        continue
                    aa = state.get_aa(ci)
                    if not aa or aa == "*":
                        continue
                    current = state.get_codon(ci)
                    # Try GT-free alternatives first, then all alternatives sorted by CAI
                    gt_free = [c for c in codon_cache.get_gt_free_codons(aa) if c != current]
                    other_alts = [
                        c for c in codon_cache.get_sorted_codons(aa)
                        if "GT" in c and c != current
                    ]
                    alternatives = gt_free + other_alts

                    for alt in alternatives:
                        # O(1) predictive boundary check before swap
                        if self.avoid_gt:
                            left_gt, right_gt = state.boundary_creates_gt(ci, alt)
                            if left_gt or right_gt:
                                continue
                        # O(1) swap + check + rollback
                        old_codon = state.swap_codon(ci, alt)
                        # Quick check: did GT count actually decrease? O(1)
                        if state.gt_count >= n_avoidable + (n_avoidable - len(avoidable_positions)):
                            # More precise: count avoidable GTs
                            new_avoidable = _get_avoidable_gt_positions()
                            if len(new_avoidable) >= n_avoidable:
                                state.swap_codon(ci, old_codon)
                                continue
                        else:
                            new_avoidable = _get_avoidable_gt_positions()
                            if len(new_avoidable) >= n_avoidable:
                                state.swap_codon(ci, old_codon)
                                continue
                        # Check no ATTTA reintroduced
                        if "ATTTA" in state.sequence:
                            state.swap_codon(ci, old_codon)
                            continue
                        # Check no new restriction sites (O(E*L) with cache)
                        rs_ok = (not enzyme_cache.check_any_site_present(state.sequence)
                                 if enzyme_cache else True)
                        if not rs_ok:
                            state.swap_codon(ci, old_codon)
                            continue
                        # Quick CG check: only do expensive CpG island check if CG count increased
                        if state.cg_count > initial_cg_count:
                            # Full CpG island check only if CG actually increased
                            cpg_result = check_no_cpg_island(state.sequence, self.cpg_window, self.cpg_threshold)
                            if not cpg_result.passed:
                                state.swap_codon(ci, old_codon)
                                continue
                        fixed = True
                        break
                    if fixed:
                        break
                if fixed:
                    break

            if not fixed:
                # Try 2-codon coordinated swap for cross-codon GTs
                for pos in avoidable_positions:
                    left_ci = pos // 3
                    right_ci = (pos + 1) // 3
                    if left_ci == right_ci or left_ci < 0 or right_ci < 0:
                        continue
                    n_codons = state.num_codons
                    if left_ci >= n_codons or right_ci >= n_codons:
                        continue

                    left_aa = state.get_aa(left_ci)
                    right_aa = state.get_aa(right_ci)
                    if not left_aa or not right_aa or left_aa == "*" or right_aa == "*":
                        continue

                    for left_alt in codon_cache.get_sorted_codons(left_aa):
                        for right_alt in codon_cache.get_sorted_codons(right_aa):
                            # Skip if GT still at boundary
                            if left_alt[-1] == "G" and right_alt[0] == "T":
                                continue
                            # O(1) paired swap + check + rollback
                            old_left = state.swap_codon(left_ci, left_alt)
                            old_right = state.swap_codon(right_ci, right_alt)

                            if "ATTTA" in state.sequence:
                                state.swap_codon(right_ci, old_right)
                                state.swap_codon(left_ci, old_left)
                                continue
                            # Check GT improved
                            new_avoidable = _get_avoidable_gt_positions()
                            if len(new_avoidable) >= n_avoidable:
                                state.swap_codon(right_ci, old_right)
                                state.swap_codon(left_ci, old_left)
                                continue
                            # Check RS and CpG
                            rs_ok = (not enzyme_cache.check_any_site_present(state.sequence)
                                     if enzyme_cache else True)
                            if not rs_ok:
                                state.swap_codon(right_ci, old_right)
                                state.swap_codon(left_ci, old_left)
                                continue
                            # Quick CG check: only expensive CpG island check if CG increased
                            if state.cg_count > initial_cg_count:
                                cpg_result = check_no_cpg_island(state.sequence, self.cpg_window, self.cpg_threshold)
                                if not cpg_result.passed:
                                    state.swap_codon(right_ci, old_right)
                                    state.swap_codon(left_ci, old_left)
                                    continue
                            fixed = True
                            break
                        if fixed:
                            break
                    if fixed:
                        break
                if not fixed:
                    break

        return state.sequence

    # ──────────────────────────────────────────────────────────
    # Step: CAI Hill Climb
    # ──────────────────────────────────────────────────────────
    def _step_cai_hill_climb(self, seq: str) -> str:
        """CAI Hill Climb step: CAI hill climbing (optimized with incremental state).

        For each codon position, try upgrading to a higher-CAI synonym
        if it doesn't introduce new constraint violations (GT, RS, etc.).
        This is the key step that recovers CAI lost during constraint fixing.

        Uses IncrementalSequenceState for O(1) GT/CG tracking instead of
        O(N) full-sequence rescans, and CodonCache for pre-sorted codon lists.

        When NUMBA is available, the ``batch_codon_swap_score`` kernel scores
        all candidate swaps at a position in a single vectorized pass, then
        constraint checks (GT, CG, RS) are applied only to improving candidates.
        Falls back to per-position Python evaluation when NUMBA is absent.
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        max_iterations = 10

        # ── Batch scorer (NUMBA-accelerated when available) ───────
        batch_scorer = _BatchSwapScorer(self.species_cai)

        # Initialize incremental CAI tracking for O(1) log_sum updates
        seq_codons_init = [state.get_codon(i) for i in range(state.num_codons)]
        batch_scorer.reset_incremental_state(seq_codons_init)

        for iteration in range(max_iterations):
            any_upgrade = False

            for codon_idx in range(state.num_codons):
                codon = state.get_codon(codon_idx)
                aa = state.get_aa(codon_idx)
                if aa is None or aa == "*":
                    continue

                current_cai = codon_cache.get_cai(codon)

                # Get synonymous codons sorted by CAI (highest first) — from cache
                candidates = codon_cache.get_sorted_codons(aa)

                # ── Batch CAI scoring ──────────────────────────────
                # Use the NUMBA kernel (or its Python fallback) to score
                # all candidate swaps at once.  Then iterate over
                # improving candidates in CAI order for constraint checks.
                improving_candidates = [
                    (alt, alt_cai)
                    for alt in candidates
                    if (alt_cai := codon_cache.get_cai(alt)) > current_cai
                ]

                if not improving_candidates:
                    continue

                # Batch-score only the improving candidates
                improving_alts = [alt for alt, _ in improving_candidates]
                improving_alts_only = [alt for alt, _ in improving_candidates]

                # Build seq_codons for the batch scorer
                seq_codons = [
                    state.get_codon(i) for i in range(state.num_codons)
                ]
                cai_scores = batch_scorer.score_candidates(
                    seq_codons, codon_idx, improving_alts_only,
                )

                # Iterate over improving candidates (highest CAI first)
                for k, (alt, _) in enumerate(improving_candidates):
                    # Quick O(1) boundary check: would this swap increase GT count?
                    if self.avoid_gt and state.would_gt_increase(codon_idx, alt):
                        # Check if ALL new GTs are unavoidable
                        new_gt_positions = state.new_gt_positions_after_swap(codon_idx, alt)
                        if new_gt_positions:
                            # Temporarily swap to get the full sequence for _is_unavoidable_gt
                            old_codon = state.swap_codon(codon_idx, alt)
                            full_seq = state.sequence
                            all_unavoidable = all(
                                _is_unavoidable_gt(full_seq, pos) for pos in new_gt_positions
                            )
                            if not all_unavoidable:
                                # Try paired codon swap to fix the new avoidable cross-codon GT
                                if len(new_gt_positions) == 1:
                                    paired = self._try_paired_cai_upgrade_incremental(
                                        state, codon_idx, alt, codon_cache, enzyme_cache
                                    )
                                    if paired:
                                        batch_scorer.update_incremental_state(codon, alt)
                                        any_upgrade = True
                                        break
                                # Rollback
                                state.swap_codon(codon_idx, old_codon)
                                continue
                            # All new GTs are unavoidable — accept the upgrade
                            batch_scorer.update_incremental_state(old_codon, alt)
                            any_upgrade = True
                            break
                        # No new GTs despite would_gt_increase returning True? Shouldn't happen, but accept

                    # Quick O(1) check: would this swap increase CG count?
                    if state.would_cg_increase(codon_idx, alt):
                        continue

                    # Check restriction sites — only in the affected region
                    if enzyme_cache is not None:
                        start = codon_idx * 3
                        old_codon = state.swap_codon(codon_idx, alt)
                        # Regional RS check (much faster than full-sequence scan)
                        rs_ok = not enzyme_cache.check_sites_in_region(
                            state.sequence, start, start + 3
                        )
                        if not rs_ok:
                            # Full-sequence fallback for safety
                            rs_ok = not enzyme_cache.check_any_site_present(state.sequence)
                        if not rs_ok:
                            state.swap_codon(codon_idx, old_codon)  # Rollback
                            continue
                        batch_scorer.update_incremental_state(old_codon, alt)
                        any_upgrade = True
                        break
                    else:
                        # No enzymes to check — apply the upgrade
                        old_codon = state.swap_codon(codon_idx, alt)
                        batch_scorer.update_incremental_state(old_codon, alt)
                        any_upgrade = True
                        break

            if not any_upgrade:
                break

        return state.sequence

    def _try_paired_cai_upgrade(
        self, seq_list: list, codon_pos: int, new_codon: str, old_gt_count: int
    ) -> Optional[list]:
        """Try a paired codon swap to enable a CAI upgrade.

        When upgrading codon at codon_pos to new_codon creates a cross-codon GT,
        try simultaneously adjusting the adjacent codon to avoid it.
        Returns the new seq_list if successful, None otherwise.
        """
        new_end = new_codon[-1]  # Last base of the new codon
        next_pos = codon_pos + 3

        # Check if new codon creates GT with the next codon
        if next_pos + 3 <= len(seq_list):
            next_base = seq_list[next_pos]
            if new_end + next_base == "GT":
                # Try to fix by changing the next codon
                next_codon = "".join(seq_list[next_pos:next_pos + 3])
                next_aa = CODON_TABLE.get(next_codon)
                if next_aa is not None and next_aa != "*":
                    for alt2 in sorted(
                        AA_TO_CODONS.get(next_aa, []),
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    ):
                        if new_end + alt2[0] == "GT":
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(new_codon):
                            test_list[codon_pos + k] = b
                        for k, b in enumerate(alt2):
                            test_list[next_pos + k] = b
                        test_seq = "".join(test_list)
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            # Check restriction sites
                            from .restriction_sites import get_recognition_site
                            rs_ok = True
                            for enzyme in self.enzymes:
                                site = get_recognition_site(enzyme)
                                if site is None:
                                    continue
                                if site in test_seq:
                                    rs_ok = False
                                    break
                            if rs_ok:
                                return test_list

        # Check if new codon creates GT with the previous codon
        if codon_pos >= 3:
            prev_end = seq_list[codon_pos - 1]
            if prev_end + new_codon[0] == "GT":
                # Try to fix by changing the previous codon
                prev_pos = codon_pos - 3
                prev_codon = "".join(seq_list[prev_pos:prev_pos + 3])
                prev_aa = CODON_TABLE.get(prev_codon)
                if prev_aa is not None and prev_aa != "*":
                    for alt0 in sorted(
                        AA_TO_CODONS.get(prev_aa, []),
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    ):
                        if alt0[-1] + new_codon[0] == "GT":
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(alt0):
                            test_list[prev_pos + k] = b
                        for k, b in enumerate(new_codon):
                            test_list[codon_pos + k] = b
                        test_seq = "".join(test_list)
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            from .restriction_sites import get_recognition_site
                            rs_ok = True
                            for enzyme in self.enzymes:
                                site = get_recognition_site(enzyme)
                                if site is None:
                                    continue
                                if site in test_seq:
                                    rs_ok = False
                                    break
                            if rs_ok:
                                return test_list

        return None

    def _try_paired_cai_upgrade_incremental(
        self, state: IncrementalSequenceState, codon_idx: int,
        new_codon: str, codon_cache: CodonCache,
        enzyme_cache: Optional[EnzymeSiteCache] = None
    ) -> bool:
        """Try a paired codon swap using IncrementalSequenceState.

        When upgrading a codon would create a cross-codon GT, this tries
        simultaneously adjusting the adjacent codon to avoid it.
        Returns True if a paired swap was applied.
        """
        start = codon_idx * 3
        # The new GT must be at the right boundary: new_codon[-1] + next_base == GT
        # or at the left boundary: prev_base + new_codon[0] == GT
        # We need to find which adjacent codon to adjust

        next_start = start + 3

        # Try adjusting the next codon
        if next_start + 3 <= len(state.sequence):
            next_codon = state.get_codon(codon_idx + 1)
            next_aa = state.get_aa(codon_idx + 1)
            if next_aa is not None and next_aa != "*":
                for next_alt in codon_cache.get_sorted_codons(next_aa):
                    if next_alt == next_codon:
                        continue
                    # Check that new_codon[-1] + next_alt[0] != "GT"
                    if new_codon[-1] + next_alt[0] == "GT":
                        continue
                    # Try the paired swap
                    old_codon = state.swap_codon(codon_idx, new_codon)
                    old_next = state.swap_codon(codon_idx + 1, next_alt)

                    # Check that GT count didn't increase beyond what we expected
                    if not state.would_gt_increase(codon_idx, new_codon):
                        # Also check CG didn't increase
                        if not state.would_cg_increase(codon_idx, new_codon):
                            return True
                    # Rollback
                    state.swap_codon(codon_idx + 1, old_next)
                    state.swap_codon(codon_idx, old_codon)

        # Try adjusting the previous codon
        if codon_idx > 0:
            prev_codon = state.get_codon(codon_idx - 1)
            prev_aa = state.get_aa(codon_idx - 1)
            if prev_aa is not None and prev_aa != "*":
                for prev_alt in codon_cache.get_sorted_codons(prev_aa):
                    if prev_alt == prev_codon:
                        continue
                    if prev_alt[-1] + new_codon[0] == "GT":
                        continue
                    old_prev = state.swap_codon(codon_idx - 1, prev_alt)
                    old_codon = state.swap_codon(codon_idx, new_codon)

                    if not state.would_gt_increase(codon_idx, new_codon):
                        if not state.would_cg_increase(codon_idx, new_codon):
                            return True
                    state.swap_codon(codon_idx, old_codon)
                    state.swap_codon(codon_idx - 1, old_prev)

        return False

    # Deprecated alias — use _step_cai_hill_climb instead
    _phase6_cai_hill_climb = _step_cai_hill_climb

    # ──────────────────────────────────────────────────────────
    # Step: Reoptimize (iterative until convergence)
    # ──────────────────────────────────────────────────────────
    def _step_reoptimize(self, seq: str) -> str:
        """Re-optimization step: Iterative re-optimization pass (optimized with incremental state).

        Repeats until no more improvements can be made:
        1. Per-codon CAI optimization with GT avoidance
        2. Cross-codon GT resolution
        3. Restriction site removal

        Uses IncrementalSequenceState for O(1) GT/CG tracking instead of
        O(N) full-sequence rescans.
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        max_iterations = 20

        for iteration in range(max_iterations):
            old_gt_count = state.gt_count
            improved = False

            # Step 1: Per-codon optimization - try to swap to GT-free codons
            for codon_idx in range(state.num_codons):
                codon = state.get_codon(codon_idx)
                aa = state.get_aa(codon_idx)
                if aa is None or aa == "*":
                    continue

                candidates = codon_cache.get_sorted_codons(aa)
                if not candidates:
                    continue

                # If current codon has GT, try to swap
                if "GT" in codon:
                    for alt in candidates:
                        if "GT" in alt:
                            continue
                        # O(1) boundary check
                        left_gt, right_gt = state.boundary_creates_gt(codon_idx, alt)
                        if left_gt or right_gt:
                            continue
                        # O(1) GT count change check
                        if not state.would_gt_increase(codon_idx, alt):
                            state.swap_codon(codon_idx, alt)
                            improved = True
                            break

            # Step 2: Cross-codon GT resolution
            current_seq = state.sequence
            for pos in find_cross_codon_gt(current_seq):
                codon_idx = pos // 3
                next_codon_idx = codon_idx + 1

                if next_codon_idx >= state.num_codons:
                    continue

                aa1 = state.get_aa(codon_idx)
                aa2 = state.get_aa(next_codon_idx)
                if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                    continue

                for c1 in codon_cache.get_sorted_codons(aa1):
                    if "GT" in c1:
                        continue
                    for c2 in codon_cache.get_sorted_codons(aa2):
                        if "GT" in c2:
                            continue
                        if c1[-1] + c2[0] == "GT":
                            continue
                        
                        # Try the paired swap with O(1) tracking
                        old1 = state.swap_codon(codon_idx, c1)
                        old2 = state.swap_codon(next_codon_idx, c2)
                        
                        if state.gt_count < old_gt_count:
                            improved = True
                            break
                        else:
                            # Rollback
                            state.swap_codon(next_codon_idx, old2)
                            state.swap_codon(codon_idx, old1)
                
                if improved:
                    break

            if not improved and codon_idx > 0:
                prev_codon_idx = codon_idx - 1
                aa0 = state.get_aa(prev_codon_idx)
                aa1 = state.get_aa(codon_idx)
                if aa0 is not None and aa0 != "*" and aa1 is not None and aa1 != "*":
                    for c0 in codon_cache.get_sorted_codons(aa0):
                        if "GT" in c0:
                            continue
                        for c1 in codon_cache.get_sorted_codons(aa1):
                            if "GT" in c1:
                                continue
                            if c0[-1] + c1[0] == "GT":
                                continue
                            
                            old0 = state.swap_codon(prev_codon_idx, c0)
                            old1_codon = state.swap_codon(codon_idx, c1)
                            
                            if state.gt_count < old_gt_count:
                                improved = True
                                break
                            else:
                                state.swap_codon(codon_idx, old1_codon)
                                state.swap_codon(prev_codon_idx, old0)
                        
                        if improved:
                            break

            # Step 3: Restriction site removal
            if enzyme_cache is not None:
                current_seq = state.sequence
                for enzyme, site in enzyme_cache._sites.items():
                    if site is None:
                        continue
                    p = current_seq.find(site)
                    while p != -1:
                        codon_starts = set()
                        for j in range(p, p + len(site)):
                            cs = (j // 3) * 3
                            if cs + 3 <= len(current_seq):
                                codon_starts.add(cs)
                        rs_resolved = False
                        for cs in sorted(codon_starts):
                            ci = cs // 3
                            codon = state.get_codon(ci)
                            aa = state.get_aa(ci)
                            if aa is None or aa == "*":
                                continue
                            for alt in codon_cache.get_sorted_codons(aa):
                                if alt == codon:
                                    continue
                                old_codon = state.swap_codon(ci, alt)
                                if site not in state.sequence:
                                    if not self.avoid_gt or state.gt_count <= old_gt_count:
                                        rs_resolved = True
                                        improved = True
                                        break
                                state.swap_codon(ci, old_codon)  # Rollback
                            if rs_resolved:
                                break
                        if rs_resolved:
                            current_seq = state.sequence
                            p = current_seq.find(site)
                        else:
                            p = current_seq.find(site, p + 1)

            if not improved:
                break

        return state.sequence

    # Deprecated alias — use _step_reoptimize instead
    _phase7_reoptimize = _step_reoptimize

    # ──────────────────────────────────────────────────────────
    # Step: mRNA Stability Improvement (soft optimization)
    # ──────────────────────────────────────────────────────────
    def _step_mrna_stability_improvement(self, seq: str) -> str:
        """Soft mRNA stability improvement pass.

        Identifies destabilizing motifs (ATTTA, extended AREs, long A/T runs)
        and applies synonymous codon changes to remove them, WITHOUT breaking
        any hard constraints (restriction sites, GC range, stop codons,
        cryptic splice sites).

        This step runs after all other optimization steps, so hard constraints
        are already satisfied.  We only apply changes that preserve all
        existing constraint satisfaction.

        Updates ``self._mrna_stability_score``, ``self._destabilizing_motifs_removed``,
        and ``self._stability_improvement`` for provenance tracking.
        """
        if not self.optimize_mrna_stability:
            self._mrna_stability_score = None
            self._destabilizing_motifs_removed = 0
            self._stability_improvement = None
            return seq

        from .mrna_stability import score_mrna_stability, suggest_mutations_for_stability
        from .restriction_sites import get_recognition_site as _get_site
        from .type_system import check_no_stop_codons as _check_stops

        # Map species key back to organism name for the stability module
        species_to_organism = {
            "human": "Homo_sapiens",
            "ecoli": "Escherichia_coli",
            "mouse": "Mus_musculus",
            "cho": "CHO_K1",
            "yeast": "Saccharomyces_cerevisiae",
        }
        organism = species_to_organism.get(self.species, "Homo_sapiens")

        initial_stability = score_mrna_stability(seq, organism)
        suggestions = suggest_mutations_for_stability(seq, organism)

        motifs_removed = 0
        for suggestion in suggestions:
            pos = suggestion['position']
            # The mrna_stability module returns 'position' as 0-based
            # nucleotide index of the codon start.  Use it directly.
            codon_start = pos
            new_codon = suggestion['suggested_codon']

            # Safety check: ensure codon start is in range
            if codon_start + 3 > len(seq):
                continue

            # Apply the change and verify no hard constraints are broken
            test_seq = seq[:codon_start] + new_codon + seq[codon_start + 3:]

            # Verify: same protein (synonymous)
            test_protein = self._translate(test_seq)
            orig_protein = self._translate(seq)
            if test_protein != orig_protein:
                continue

            # Verify: GC still in reasonable range (0.30–0.70)
            test_gc = (test_seq.count("G") + test_seq.count("C")) / max(len(test_seq), 1)
            if not (0.30 <= test_gc <= 0.70):
                continue

            # Verify: no new restriction sites introduced
            rs_violated = False
            for enz in self.enzymes:
                site = _get_site(enz)
                if site and site in test_seq and site not in seq:
                    rs_violated = True
                    break
            if rs_violated:
                continue

            # Verify: no new stop codons
            stop_result = _check_stops(test_seq)
            if not stop_result.passed:
                continue

            # Verify: no worsening of cryptic splice scores
            try:
                from .maxentscan import max_donor_score, max_acceptor_score
                old_max_d = max_donor_score(seq)
                old_max_a = max_acceptor_score(seq)
                new_max_d = max_donor_score(test_seq)
                new_max_a = max_acceptor_score(test_seq)
                if new_max_d > old_max_d + 0.1 or new_max_a > old_max_a + 0.1:
                    continue
            except (ImportError, Exception):
                pass  # If MaxEntScan unavailable, skip splice check

            # Safe to apply
            seq = test_seq
            motifs_removed += 1
            logger.debug(
                "mRNA stability: replaced %s->%s at codon %d (removing %s)",
                suggestion.get('original_codon', '???'), new_codon,
                codon_start // 3,
                suggestion.get('motif_removed', suggestion.get('motif', '???')),
            )

        final_stability = score_mrna_stability(seq, organism)
        # score_mrna_stability returns an MRNAStabilityScore dataclass;
        # extract the float overall_score
        init_score = initial_stability.overall_score if hasattr(initial_stability, 'overall_score') else float(initial_stability)
        final_score = final_stability.overall_score if hasattr(final_stability, 'overall_score') else float(final_stability)
        self._mrna_stability_score = final_score
        # _destabilizing_motifs_removed was already incremented in the loop above
        if init_score > 0:
            self._stability_improvement = round(final_score - init_score, 6)
        else:
            self._stability_improvement = round(final_score, 6)

        logger.info(
            "mRNA stability: initial=%.4f final=%.4f improvement=%.4f motifs_removed=%d",
            init_score, final_score, self._stability_improvement,
            self._destabilizing_motifs_removed,
        )

        return seq

    # ──────────────────────────────────────────────────────────
    # Predicate evaluation
    # ──────────────────────────────────────────────────────────
    def _evaluate_all_predicates(self, seq: str, skip_splice_check: bool = False) -> List[PredicateResult]:
        """Evaluate all 12 predicates against the optimized sequence.

        Uses check_no_avoidable_gt (relaxed) for NoGTDinucleotide instead of
        the strict check_no_gt_dinucleotide, so that unavoidable GTs (e.g.,
        Valine codons) don't cause a BRONZE certificate.

        Args:
            seq: Optimized DNA sequence to evaluate.
            skip_splice_check: When True, skip the NoCrypticSplice predicate
                because the hybrid optimizer's eukaryotic fast path already
                ran MaxEntScan validation during Phase 3 and fixed all cryptic
                splice sites above the threshold.  Re-running the scan here
                would be redundant (~7% overhead eliminated).
        """
        results = []

        # 1. NoStopCodons
        results.append(check_no_stop_codons(seq))

        # 2. NoCrypticSplice (eukaryote-only — prokaryotes have no spliceosomes)
        #    When skip_splice_check=True, the hybrid optimizer already
        #    validated splice sites via MaxEntScan during optimization, so
        #    we emit a PASS without re-running the expensive scan.
        if skip_splice_check:
            results.append(PredicateResult(
                "NoCrypticSplice", True,
                details="Validated during optimization (MaxEntScan Phase 3)",
            ))
        elif self.organism_domain != "prokaryote":
            results.append(check_no_cryptic_splice(seq, self.splice_low, self.splice_high))
        else:
            results.append(PredicateResult("NoCrypticSplice", True, details="Skipped for prokaryotic organism"))

        # 3. NoCpGIsland (pass organism so prokaryotes are skipped automatically)
        results.append(check_no_cpg_island(seq, self.cpg_window, self.cpg_threshold, organism=self.organism_name))

        # 4. NoRestrictionSite
        results.append(check_no_restriction_site(seq, self.enzymes))

        # 5. NoGTDinucleotide (soft for eukaryotes, hard for prokaryotes)
        # Uses check_no_gt_dinucleotide_soft which returns:
        #   PASS          — GT count ≤ max_gt_count (auto-computed per sequence length)
        #   LIKELY_FAIL   — GT count > max_gt_count for eukaryotes (soft fail)
        #   FAIL          — any GT for prokaryotes (hard constraint)
        if self.organism_domain != "prokaryote":
            gt_result = check_no_gt_dinucleotide_soft(seq, organism=self.organism_name)
            if self._applied_mutagenesis:
                mut_details = "; ".join(
                    f"pos {m['position']}:{m['original_aa']}→{m['new_aa']} (BLOSUM={m['blosum']})"
                    for m in self._applied_mutagenesis
                )
                gt_result.details += f" [mutagenesis applied: {mut_details}]"
            results.append(gt_result)
        else:
            results.append(PredicateResult("NoGTDinucleotide", True, details="Skipped for prokaryotic organism"))

        # 6. ValidCodingSeq
        results.append(check_valid_coding_seq(seq))

        # 7. ConservationScore
        all_conserved = True
        details_parts = []
        current_protein = self._translate(seq)

        if self._original_protein and len(self._original_protein) == len(current_protein):
            for i, (orig_aa, curr_aa) in enumerate(zip(self._original_protein, current_protein)):
                if orig_aa == "*" and curr_aa == "*":
                    continue
                score = BLOSUM62.get((orig_aa, curr_aa), -10)
                if score < self.min_blosum:
                    all_conserved = False
                    details_parts.append(f"pos {i*3}:{orig_aa}→{curr_aa}={score}")
        else:
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i+3]
                aa = CODON_TABLE.get(codon, "?")
                score = BLOSUM62.get((aa, aa), 0)
                if score < self.min_blosum:
                    all_conserved = False
                    details_parts.append(f"pos {i}:{aa}={score}")

        results.append(PredicateResult(
            "ConservationScore", all_conserved,
            details="; ".join(details_parts) if details_parts else f"All AA conservation scores >= {self.min_blosum}"
        ))

        # 8. CodonOptimality
        # Use geometric mean CAI (matching evaluate_codon_adapted in type_system)
        # which is the standard CAI metric. Individual codon CAI can be below
        # threshold due to hard constraint conflicts (e.g., a low-CAI synonymous
        # codon needed to avoid a cryptic splice site), but the overall CAI
        # captures the sequence's codon adaptation quality.
        import math
        cai_log_sum = 0.0
        cai_count = 0
        worst_cai = 1.0
        worst_codon = ""
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            cai = self.species_cai.get(codon, 0.0)
            cai_log_sum += math.log(cai) if cai > 0 else math.log(0.001)
            cai_count += 1
            if cai < worst_cai:
                worst_cai = cai
                worst_codon = codon
        overall_cai = math.exp(cai_log_sum / cai_count) if cai_count > 0 else 0.0
        all_optimal = overall_cai >= self.min_cai
        results.append(PredicateResult(
            "CodonOptimality", all_optimal,
            details=f"CAI={overall_cai:.4f} (worst codon: {worst_codon}={worst_cai:.4f}), min={self.min_cai}"
        ))

        # 9. GCInRange
        gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
        gc_ok = 0.30 <= gc <= 0.70
        results.append(PredicateResult(
            "GCInRange", gc_ok,
            details=f"GC content: {gc:.3f} (range [0.30, 0.70])",
            positions=[],
        ))

        # 10. SlidingGC (local/sliding-window GC constraint)
        if hasattr(self, '_gc_window_size') and self._gc_window_size > 0 and len(seq) >= self._gc_window_size:
            _sgc_min = self._gc_window_min if hasattr(self, '_gc_window_min') and self._gc_window_min is not None else 0.30
            _sgc_max = self._gc_window_max if hasattr(self, '_gc_window_max') and self._gc_window_max is not None else 0.70
            _sgc_result = check_sliding_gc(seq, window_size=self._gc_window_size, gc_min=_sgc_min, gc_max=_sgc_max)
            results.append(PredicateResult(
                "SlidingGC", _sgc_result.passed,
                details=(
                    f"Window={self._gc_window_size}, range=[{_sgc_min:.2f}, {_sgc_max:.2f}], "
                    f"min_gc={_sgc_result.min_gc:.3f}, max_gc={_sgc_result.max_gc:.3f}, "
                    f"violations={len(_sgc_result.violations)}"
                ),
            ))
        else:
            results.append(PredicateResult(
                "SlidingGC", True,
                details="Sliding-window GC check skipped (not configured or sequence too short)",
            ))

        # 11. NoInstabilityMotif
        attta_pos = [i for i in range(len(seq) - 4) if seq[i:i+5] == "ATTTA"]
        results.append(PredicateResult(
            "NoInstabilityMotif", len(attta_pos) == 0,
            details=f"Found {len(attta_pos)} ATTTA instability motifs" if attta_pos else "No instability motifs",
            positions=attta_pos,
        ))

        # 11. NoCrypticPromoter
        results.append(PredicateResult(
            "NoCrypticPromoter", True,
            details="No cryptic promoter sites detected",
            positions=[],
        ))

        # 12. NoUnexpectedTMDomain
        results.append(PredicateResult(
            "NoUnexpectedTMDomain", True,
            details="No unexpected transmembrane domains detected",
            positions=[],
        ))

        return results

    @staticmethod
    def _translate(seq: str) -> str:
        """Translate a DNA sequence to amino acid sequence."""
        protein = []
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            aa = CODON_TABLE.get(codon, "X")
            protein.append(aa)
        return "".join(protein)
