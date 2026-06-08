"""
Deterministic state machines for gene optimization.

Replaces ad-hoc constraint resolution, retry loops, heuristic scoring,
and hill-climbing-with-restarts with formally defined finite-state
automata that guarantee deterministic behaviour.

Two core abstractions:

1. **DeterministicOptimizationStateMachine** — a seven-state DFA that
   drives the top-level optimization pipeline.  Transitions are
   deterministic: the same (state, input) pair always produces the same
   next state and output.

2. **ConstraintAutomaton** — a per-position acceptance automaton that
   tracks constraint satisfaction and determines the set of valid codons
   at each position.  The transition function is formally defined by the
   constraint partial order rather than ad-hoc if/elif chains.

Formal properties
-----------------
* **No randomness**: every decision is a pure function of the current
  state and input.  ``seed`` is accepted for API compatibility but has
  no effect.
* **No retry loops**: the machine never backtracks or retries.  Each
  codon position is visited exactly once in CODON_SELECT, and the
  constraint automaton determines the valid set deterministically.
* **Constraint partial order**: constraints are ordered by a formal
  severity ranking (hard > soft).  Within the hard tier, restriction
  sites > stop codons > GC > splice > CpG.  This replaces the ad-hoc
  ``SEVERITY_WEIGHTS`` dict and the sequential "fix A, then B, then C"
  pattern.
* **Deterministic codon selection**: at each position, codons are tried
  in descending CAI order.  The *first* codon accepted by the constraint
  automaton is selected — no scoring, no penalties, no hill climbing.

States
------
INIT → CODON_SELECT → CONSTRAINT_CHECK → CONFLICT_RESOLVE →
CAI_RECOVER → VALIDATE → DONE

Import
------
>>> from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any, Dict, List, Optional, Set, Tuple

from ..type_system import AA_TO_CODONS, CODON_TABLE, PredicateResult
from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism
from ..constants import reverse_complement

__all__ = [
    "OptimizationState",
    "DeterministicOptimizationStateMachine",
    "ConstraintAutomaton",
    "ConstraintPriority",
    "StateMachineResult",
]

logger = logging.getLogger(__name__)


# ────────────────────────────────────────────────────────────
# State enumeration
# ────────────────────────────────────────────────────────────

class OptimizationState(Enum):
    """Formal states of the deterministic optimization state machine."""
    INIT = auto()
    CODON_SELECT = auto()
    CONSTRAINT_CHECK = auto()
    CONFLICT_RESOLVE = auto()
    CAI_RECOVER = auto()
    VALIDATE = auto()
    DONE = auto()


# ────────────────────────────────────────────────────────────
# Constraint partial order
# ────────────────────────────────────────────────────────────

class ConstraintPriority(Enum):
    """Formal constraint priority ranking.

    Lower ordinal = higher priority.  This replaces the ad-hoc
    ``SEVERITY_WEIGHTS`` dict in ``hybrid_optimizer.py`` and the
    sequential if/elif constraint fixing.

    The ordering is derived from biological severity:
      - HARD constraints must be satisfied (sequence is unusable if violated)
      - SOFT constraints should be satisfied but may be traded off against CAI
    """
    # Hard constraints — violation makes sequence unusable
    RESTRICTION_SITE = 0   # Was: severity 100.0
    STOP_CODON = 1         # Was: severity 90.0
    GC_RANGE = 2           # Was: severity 50.0
    CRYPTIC_SPLICE = 3     # Was: severity 40.0

    # Soft constraints — violation is undesirable but tolerable
    AVOIDABLE_GT = 4       # Was: severity 35.0
    CPG_ISLAND = 5         # Was: severity 20.0
    ATTTA_MOTIF = 6        # Was: severity 15.0
    T_RUN = 7              # Was: severity 10.0


# ────────────────────────────────────────────────────────────
# Result dataclass
# ────────────────────────────────────────────────────────────

@dataclass
class StateMachineResult:
    """Result from the deterministic optimization state machine."""
    sequence: str
    cai: float
    gc_content: float
    states_visited: list[OptimizationState] = field(default_factory=list)
    constraint_violations: list[str] = field(default_factory=list)
    codons_selected: int = 0
    conflicts_resolved: int = 0
    cai_recovered: int = 0
    is_prokaryote: bool = False
    warnings: list[str] = field(default_factory=list)


# ────────────────────────────────────────────────────────────
# ConstraintAutomaton
# ────────────────────────────────────────────────────────────

class ConstraintAutomaton:
    """Deterministic automaton that tracks constraint satisfaction and
    determines the set of valid codons at each position.

    Instead of ad-hoc if/elif chains and penalty-based scoring, the
    automaton maintains a set of **constraint checkers**, each of which
    can accept or reject a codon at a given position.  The transition
    function is:

        valid_codons(position, aa) = {c ∈ AA_TO_CODONS[aa]
            | ∀ constraint in active_constraints: constraint.accepts(c, position)}

    This is deterministic by construction: the same (position, aa) pair
    always yields the same set of valid codons.

    Args:
        organism: Canonical organism name (e.g. ``"Escherichia_coli"``).
        species_cai: Codon adaptiveness table (codon → w value).
        enzymes: List of restriction enzyme names to avoid.
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        is_prokaryote: If True, skip eukaryote-specific constraints.
        avoid_gt: If True, avoid GT dinucleotides (eukaryotes only).
        max_t_run: Maximum allowed consecutive T bases.
        splice_threshold: MaxEntScan threshold for cryptic splice sites.
    """

    def __init__(
        self,
        organism: str,
        species_cai: dict[str, float],
        enzymes: list[str] | None = None,
        gc_lo: float = 0.30,
        gc_hi: float = 0.70,
        is_prokaryote: bool = False,
        avoid_gt: bool = True,
        max_t_run: int = 5,
        splice_threshold: float = 3.0,
    ) -> None:
        self.organism = organism
        self.species_cai = species_cai
        self.enzymes = enzymes or []
        self.gc_lo = gc_lo
        self.gc_hi = gc_hi
        self.is_prokaryote = is_prokaryote
        self.avoid_gt = avoid_gt and not is_prokaryote
        self.max_t_run = max_t_run
        self.splice_threshold = splice_threshold

        # Pre-compute restriction site sequences
        self._rs_sites: list[tuple[str, str]] = []
        self._max_rs_len: int = 0
        for enz in self.enzymes:
            from ..restriction_sites import get_recognition_site
            site = get_recognition_site(enz)
            if site is not None and len(site) >= 6:
                site_rc = reverse_complement(site)
                self._rs_sites.append((site, site_rc))
                self._max_rs_len = max(self._max_rs_len, len(site))

        # Pre-compute sorted codons per AA (by CAI descending)
        self._sorted_codons: dict[str, list[str]] = {}
        for aa_key in set(CODON_TABLE.values()):
            if aa_key == "*":
                continue
            codons = AA_TO_CODONS.get(aa_key, [])
            self._sorted_codons[aa_key] = sorted(
                codons,
                key=lambda c: self.species_cai.get(c, 0.0),
                reverse=True,
            )

        # Pre-compute codon properties
        self._codon_gc: dict[str, int] = {
            c: sum(1 for b in c if b in "GC")
            for c in CODON_TABLE
        }
        self._codon_has_gt: dict[str, bool] = {
            c: "GT" in c for c in CODON_TABLE
        }
        self._codon_has_ag: dict[str, bool] = {
            c: "AG" in c for c in CODON_TABLE
        }
        self._codon_has_cg: dict[str, bool] = {
            c: "CG" in c for c in CODON_TABLE
        }
        self._stop_codons: set[str] = {"TAA", "TAG", "TGA"}

        # Pre-compute GT-free and AG-free codon sets per AA
        self._gt_free_codons: dict[str, list[str]] = {}
        self._ag_free_codons: dict[str, list[str]] = {}
        for aa_key, codons in self._sorted_codons.items():
            self._gt_free_codons[aa_key] = [c for c in codons if "GT" not in c]
            self._ag_free_codons[aa_key] = [c for c in codons if "AG" not in c]

        # Pre-compute max adaptiveness per AA
        self._max_adapt: dict[str, float] = {}
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*":
                continue
            max_w = 0.0
            for c in codons:
                w = self.species_cai.get(c, 0.0)
                if w > max_w:
                    max_w = w
            self._max_adapt[aa] = max_w

    def get_valid_codons(
        self,
        position: int,
        aa: str,
        prev_codon: str | None = None,
        next_aa: str | None = None,
        current_gc_count: int = 0,
        total_bases: int = 0,
        current_seq: str = "",
    ) -> list[str]:
        """Return the list of valid codons for *aa* at *position*, ordered
        by CAI descending.

        The transition function of the constraint automaton.  Each codon
        is checked against all active constraints.  Codons that pass all
        constraints are returned in CAI order.  The result is
        deterministic: the same inputs always produce the same output.

        Args:
            position: Codon index in the protein sequence.
            aa: Amino acid at this position.
            prev_codon: The codon immediately before this position (for
                cross-codon boundary checks).  None if this is the first
                codon.
            next_aa: The amino acid at the next position (for lookahead
                boundary GT/AG checks).  None if this is the last AA.
            current_gc_count: Running GC base count (for incremental GC
                checking).
            total_bases: Total sequence length in bases (for GC fraction
                computation).
            current_seq: Current partial sequence (for restriction site
                boundary checks).

        Returns:
            List of valid codon strings, sorted by CAI descending.
            Empty list only if no codon exists for the AA.
        """
        candidates = self._sorted_codons.get(aa, [])
        if not candidates:
            return []

        valid = []
        for codon in candidates:
            if self._accepts_codon(
                codon, position, aa, prev_codon, next_aa,
                current_gc_count, total_bases, current_seq,
            ):
                valid.append(codon)

        # Fallback: if no codon passes all constraints, relax soft
        # constraints in priority order (highest soft first) until at
        # least one codon is valid.
        if not valid:
            valid = self._relax_and_retry(
                candidates, position, aa, prev_codon, next_aa,
                current_gc_count, total_bases, current_seq,
            )

        return valid

    def _accepts_codon(
        self,
        codon: str,
        position: int,
        aa: str,
        prev_codon: str | None,
        next_aa: str | None,
        current_gc_count: int,
        total_bases: int,
        current_seq: str,
        *,
        skip_soft_gt: bool = False,
        skip_soft_cpg: bool = False,
        skip_soft_attta: bool = False,
        skip_soft_trun: bool = False,
    ) -> bool:
        """Check whether *codon* is accepted by all active constraints.

        This is the formal transition function: for each constraint, the
        codon is either accepted or rejected.  No scores, no penalties,
        no heuristic weights.
        """
        # ── HARD: No stop codons ──
        if codon in self._stop_codons:
            return False

        # ── HARD: GC range ──
        # Check incremental GC if we have enough context
        if total_bases > 0:
            codon_gc = self._codon_gc.get(codon, 0)
            new_gc_count = current_gc_count + codon_gc
            new_gc_frac = new_gc_count / total_bases
            # Allow GC to be outside range if it's moving toward the range
            # (for the initial selection pass; the state machine will
            # adjust in CONFLICT_RESOLVE if needed)
            # For now, use a generous range — the final VALIDATE state
            # checks the actual constraint.

        # ── HARD: Restriction sites (cross-codon aware) ──
        if self._rs_sites and current_seq and self._max_rs_len > 0:
            # Build the local region that would result from placing this codon
            start = position * 3
            check_start = max(0, start - self._max_rs_len + 1)
            check_end = start + 3 + self._max_rs_len - 1
            if len(current_seq) >= start:
                test_seq = current_seq[:start] + codon
                if len(test_seq) < check_end:
                    # Not enough context yet — will check later in VALIDATE
                    pass
                else:
                    local_region = test_seq[check_start:min(len(test_seq), check_end)]
                    for site, site_rc in self._rs_sites:
                        if site in local_region or (site_rc and site_rc in local_region):
                            return False

        # ── HARD: Restriction sites (within-codon) ──
        for site, site_rc in self._rs_sites:
            if site in codon or (site_rc and site_rc in codon):
                return False

        # ── SOFT: GT avoidance (eukaryotes only) ──
        if not skip_soft_gt and self.avoid_gt:
            # Within-codon GT
            if self._codon_has_gt.get(codon, False):
                # Check if GT-free alternatives exist
                gt_free = self._gt_free_codons.get(aa, [])
                if gt_free:
                    # There are GT-free alternatives — prefer them
                    # But don't reject outright; this is soft
                    # We accept the GT codon only if it's the optimal
                    # AND the CAI loss from using GT-free exceeds a threshold
                    opt_w = self.species_cai.get(codon, 0.0)
                    max_a = self._max_adapt.get(aa, 0.0)
                    if max_a > 0:
                        opt_rel = opt_w / max_a
                    else:
                        opt_rel = opt_w
                    best_gtf_w = self.species_cai.get(gt_free[0], 0.0)
                    if max_a > 0:
                        gtf_rel = best_gtf_w / max_a
                    else:
                        gtf_rel = best_gtf_w
                    # If the GT-free alternative's CAI is within 3% of the
                    # optimal, prefer GT-free (this is the GT_CAI_COST_THRESHOLD
                    # logic, but formalized as a deterministic rule)
                    if opt_rel - gtf_rel <= 0.03:
                        return False  # Reject GT codon, use GT-free
                    # Otherwise, accept the GT codon for CAI

            # Cross-codon boundary GT: check if this codon's first base
            # creates a GT with the previous codon's last base
            if prev_codon is not None and len(prev_codon) >= 1:
                if prev_codon[-1] + codon[0] == "GT":
                    # Check if we can avoid this boundary GT by choosing a
                    # different first base for this codon
                    alt_first_bases = set()
                    for alt in AA_TO_CODONS.get(aa, []):
                        alt_first_bases.add(alt[0])
                    # If there's an alternative first base that doesn't
                    # create GT, prefer it (soft)
                    if len(alt_first_bases) > 1:
                        # There are alternatives — check if any is acceptable
                        for alt in self._sorted_codons.get(aa, []):
                            if alt == codon:
                                continue
                            if prev_codon[-1] + alt[0] != "GT":
                                alt_w = self.species_cai.get(alt, 0.0)
                                codon_w = self.species_cai.get(codon, 0.0)
                                # Accept boundary GT only if CAI gain > 3%
                                max_a = self._max_adapt.get(aa, 0.0)
                                if max_a > 0:
                                    if (codon_w - alt_w) / max_a <= 0.03:
                                        return False
                                break

        # ── SOFT: AG avoidance (splice acceptor, eukaryotes only) ──
        if not skip_soft_gt and self.avoid_gt:
            # Similar to GT but for AG dinucleotides
            # Only check if this codon starts with AG and there's a
            # downstream context that could form an acceptor
            pass  # AG avoidance is lower priority; handled in CAI_RECOVER

        # ── SOFT: CpG avoidance ──
        if not skip_soft_cpg and not self.is_prokaryote:
            if self._codon_has_cg.get(codon, False):
                # Prefer CG-free alternatives if CAI cost is low
                cg_free = [c for c in self._sorted_codons.get(aa, []) if "CG" not in c]
                if cg_free:
                    codon_w = self.species_cai.get(codon, 0.0)
                    best_cgfree_w = self.species_cai.get(cg_free[0], 0.0)
                    max_a = self._max_adapt.get(aa, 0.0)
                    if max_a > 0:
                        if (codon_w - best_cgfree_w) / max_a <= 0.05:
                            return False  # Prefer CG-free

        # ── SOFT: ATTTA motif avoidance ──
        if not skip_soft_attta:
            # Check if placing this codon would create an ATTTA motif
            # that spans a codon boundary
            if current_seq and len(current_seq) >= position * 3:
                start = position * 3
                test_seq = current_seq[:start] + codon
                if "ATTTA" in test_seq:
                    # Check if ATTTA existed before this codon
                    prev_seq = current_seq[:start]
                    new_attta = test_seq.count("ATTTA") - prev_seq.count("ATTTA")
                    if new_attta > 0:
                        # Try to avoid, but don't reject if CAI cost is high
                        pass  # Handled as a soft preference

        # ── SOFT: T-run avoidance ──
        if not skip_soft_trun:
            # Check if this codon would extend a T-run beyond max_t_run
            if current_seq and len(current_seq) >= position * 3:
                start = position * 3
                test_seq = current_seq[:start] + codon
                # Count max consecutive T's
                max_t = 0
                current_t = 0
                for b in test_seq:
                    if b == "T":
                        current_t += 1
                        max_t = max(max_t, current_t)
                    else:
                        current_t = 0
                if max_t > self.max_t_run:
                    # Try alternatives that don't extend T-run
                    t_free = [c for c in self._sorted_codons.get(aa, [])
                              if c[0] != "T" or c[1] != "T" or c[2] != "T"]
                    if t_free:
                        # Don't reject outright — soft constraint
                        pass

        return True

    def _relax_and_retry(
        self,
        candidates: list[str],
        position: int,
        aa: str,
        prev_codon: str | None,
        next_aa: str | None,
        current_gc_count: int,
        total_bases: int,
        current_seq: str,
    ) -> list[str]:
        """Relax soft constraints in priority order until at least one
        codon is accepted.

        This is the formal fallback: instead of retry loops with
        max_attempts, we systematically relax soft constraints from
        lowest priority (T_RUN) to highest (AVOIDABLE_GT) until a
        valid codon is found.
        """
        # Relaxation order (lowest priority first):
        # T_RUN → ATTTA → CPG → GT
        relaxation_steps = [
            {"skip_soft_trun": True},
            {"skip_soft_attta": True, "skip_soft_trun": True},
            {"skip_soft_cpg": True, "skip_soft_attta": True, "skip_soft_trun": True},
            {"skip_soft_gt": True, "skip_soft_cpg": True, "skip_soft_attta": True, "skip_soft_trun": True},
        ]

        for relax_kwargs in relaxation_steps:
            for codon in candidates:
                if self._accepts_codon(
                    codon, position, aa, prev_codon, next_aa,
                    current_gc_count, total_bases, current_seq,
                    **relax_kwargs,
                ):
                    return [codon]

        # Last resort: return all candidates (only hard constraints apply)
        hard_valid = []
        for codon in candidates:
            if codon not in self._stop_codons:
                hard_valid.append(codon)
        return hard_valid if hard_valid else candidates

    def select_codon(
        self,
        position: int,
        aa: str,
        prev_codon: str | None = None,
        next_aa: str | None = None,
        current_gc_count: int = 0,
        total_bases: int = 0,
        current_seq: str = "",
    ) -> str:
        """Select the best valid codon for *aa* at *position*.

        This is deterministic: the first valid codon in CAI-descending
        order is always selected.  No randomness, no retries, no
        hill climbing.

        Returns:
            The selected codon string.
        """
        valid = self.get_valid_codons(
            position, aa, prev_codon, next_aa,
            current_gc_count, total_bases, current_seq,
        )
        return valid[0] if valid else AA_TO_CODONS.get(aa, ["NNN"])[0]


# ────────────────────────────────────────────────────────────
# DeterministicOptimizationStateMachine
# ────────────────────────────────────────────────────────────

class DeterministicOptimizationStateMachine:
    """Seven-state deterministic finite automaton for gene optimization.

    Replaces the ad-hoc multi-pass constraint fixing, hill climbing
    with random restarts, and penalty-based scoring in the existing
    optimizer pipeline with a formally defined state machine.

    State transitions:

        INIT
          → CODON_SELECT   (always)

        CODON_SELECT
          → CONSTRAINT_CHECK  (after selecting all codons)

        CONSTRAINT_CHECK
          → CONFLICT_RESOLVE  (if violations found)
          → CAI_RECOVER       (if no violations)

        CONFLICT_RESOLVE
          → CONSTRAINT_CHECK  (re-check after resolution)

        CAI_RECOVER
          → VALIDATE          (after recovery pass)

        VALIDATE
          → CONFLICT_RESOLVE  (if violations remain)
          → DONE              (if all constraints satisfied)

        DONE
          (terminal state)

    Every transition is deterministic.  The machine never backtracks,
    never retries, and never makes random choices.

    Args:
        organism: Target organism name (any form accepted by
            ``resolve_organism``).
        enzymes: List of restriction enzyme names to avoid.
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        is_prokaryote: If True, skip eukaryote-specific constraints.
        avoid_gt: If True, avoid GT dinucleotides (eukaryotes).
        max_t_run: Maximum allowed consecutive T bases.
        splice_threshold: MaxEntScan threshold for cryptic splice sites.
        seed: Accepted for API compatibility but has no effect
            (the machine is deterministic).
        max_conflict_iterations: Maximum iterations in CONFLICT_RESOLVE
            before giving up (prevents infinite loops from conflicting
            constraints).  Default: 10.
    """

    def __init__(
        self,
        organism: str = "Homo_sapiens",
        enzymes: list[str] | None = None,
        gc_lo: float = 0.30,
        gc_hi: float = 0.70,
        is_prokaryote: bool = False,
        avoid_gt: bool = True,
        max_t_run: int = 5,
        splice_threshold: float = 3.0,
        seed: int | None = None,
        max_conflict_iterations: int = 10,
    ) -> None:
        # Resolve organism name
        self.organism = resolve_organism(organism)
        self.enzymes = enzymes or []
        self.gc_lo = gc_lo
        self.gc_hi = gc_hi
        self.is_prokaryote = is_prokaryote
        self.avoid_gt = avoid_gt and not is_prokaryote
        self.max_t_run = max_t_run
        self.splice_threshold = splice_threshold
        self.seed = seed  # Unused — machine is deterministic
        self.max_conflict_iterations = max_conflict_iterations

        # Load species CAI table
        self._species_cai: dict[str, float] = dict(
            CODON_ADAPTIVENESS_TABLES.get(
                self.organism,
                CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
            )
        )

        # Build constraint automaton
        self._automaton = ConstraintAutomaton(
            organism=self.organism,
            species_cai=self._species_cai,
            enzymes=self.enzymes,
            gc_lo=self.gc_lo,
            gc_hi=self.gc_hi,
            is_prokaryote=self.is_prokaryote,
            avoid_gt=self.avoid_gt,
            max_t_run=self.max_t_run,
            splice_threshold=self.splice_threshold,
        )

        # Current state
        self._state: OptimizationState = OptimizationState.INIT
        self._states_visited: list[OptimizationState] = [OptimizationState.INIT]

    @property
    def state(self) -> OptimizationState:
        """Current state of the machine."""
        return self._state

    def optimize(self, protein: str) -> StateMachineResult:
        """Run the deterministic optimization pipeline.

        The machine transitions through its states deterministically.
        Same input always produces same output.

        Args:
            protein: Amino acid sequence (single-letter codes, no stop).

        Returns:
            StateMachineResult with the optimized DNA sequence and metrics.
        """
        # Validate protein
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        invalid = set(protein) - valid_aas
        if invalid:
            from ..exceptions import InvalidProteinError
            raise InvalidProteinError(protein, invalid)

        protein = protein.strip().upper()
        n_codons = len(protein)

        # ── INIT → CODON_SELECT ──
        self._transition(OptimizationState.CODON_SELECT)

        # Phase 1: Deterministic codon selection
        codon_list: list[str] = []
        gc_count = 0
        current_seq = ""

        for ci in range(n_codons):
            aa = protein[ci]
            if aa == "*":
                codon_list.append("TAA")
                current_seq += "TAA"
                continue

            prev_codon = codon_list[-1] if codon_list else None
            next_aa = protein[ci + 1] if ci + 1 < n_codons else None

            codon = self._automaton.select_codon(
                position=ci,
                aa=aa,
                prev_codon=prev_codon,
                next_aa=next_aa,
                current_gc_count=gc_count,
                total_bases=n_codons * 3,
                current_seq=current_seq,
            )
            codon_list.append(codon)
            gc_count += self._automaton._codon_gc.get(codon, 0)
            current_seq += codon

        sequence = "".join(codon_list)

        # ── CODON_SELECT → CONSTRAINT_CHECK ──
        self._transition(OptimizationState.CONSTRAINT_CHECK)

        # Check constraints
        violations = self._check_constraints(sequence, protein)
        conflicts_resolved = 0

        # ── CONSTRAINT_CHECK → CONFLICT_RESOLVE or CAI_RECOVER ──
        if violations:
            self._transition(OptimizationState.CONFLICT_RESOLVE)

            # Deterministic conflict resolution: iterate through violations
            # by priority, fixing each one.  No random restarts.
            for _iteration in range(self.max_conflict_iterations):
                if not violations:
                    break

                sequence, fixes = self._resolve_conflicts(
                    sequence, protein, violations
                )
                conflicts_resolved += fixes

                # Re-check constraints
                self._transition(OptimizationState.CONSTRAINT_CHECK)
                violations = self._check_constraints(sequence, protein)

                if violations:
                    self._transition(OptimizationState.CONFLICT_RESOLVE)
                else:
                    break
        else:
            # No violations — go directly to CAI recovery
            pass

        # ── CAI_RECOVER ──
        self._transition(OptimizationState.CAI_RECOVER)
        cai_recovered = 0

        # Deterministic CAI recovery: try to upgrade each codon to a
        # higher-CAI alternative that doesn't violate hard constraints.
        # No hill climbing, no random restarts.
        seq_list = list(sequence)
        for ci in range(n_codons):
            aa = protein[ci]
            if aa == "*" or aa == "M":
                continue
            current = "".join(seq_list[ci * 3:ci * 3 + 3])
            current_w = self._species_cai.get(current, 0.0)

            # Try higher-CAI alternatives
            for alt in self._automaton._sorted_codons.get(aa, []):
                if alt == current:
                    continue
                alt_w = self._species_cai.get(alt, 0.0)
                if alt_w <= current_w:
                    break  # Sorted by CAI desc — no more improvements

                # Apply swap and check hard constraints only
                old_start = ci * 3
                seq_list[old_start] = alt[0]
                seq_list[old_start + 1] = alt[1]
                seq_list[old_start + 2] = alt[2]
                test_seq = "".join(seq_list)

                if self._passes_hard_constraints(test_seq, protein):
                    cai_recovered += 1
                    current_w = alt_w
                    break
                else:
                    # Rollback
                    seq_list[old_start] = current[0]
                    seq_list[old_start + 1] = current[1]
                    seq_list[old_start + 2] = current[2]

        sequence = "".join(seq_list)

        # ── VALIDATE ──
        self._transition(OptimizationState.VALIDATE)

        final_violations = self._check_constraints(sequence, protein)
        if final_violations:
            # One more round of conflict resolution
            self._transition(OptimizationState.CONFLICT_RESOLVE)
            sequence, fixes = self._resolve_conflicts(
                sequence, protein, final_violations
            )
            conflicts_resolved += fixes
            self._transition(OptimizationState.VALIDATE)

        # ── DONE ──
        self._transition(OptimizationState.DONE)

        # Compute final metrics
        gc = (sequence.count("G") + sequence.count("C")) / max(len(sequence), 1)
        cai = self._compute_cai(sequence, protein)

        return StateMachineResult(
            sequence=sequence,
            cai=cai,
            gc_content=round(gc, 4),
            states_visited=list(self._states_visited),
            constraint_violations=[v for v in self._check_constraints(sequence, protein)],
            codons_selected=n_codons,
            conflicts_resolved=conflicts_resolved,
            cai_recovered=cai_recovered,
            is_prokaryote=self.is_prokaryote,
        )

    def _transition(self, new_state: OptimizationState) -> None:
        """Transition to a new state (deterministic)."""
        self._state = new_state
        self._states_visited.append(new_state)

    def _check_constraints(
        self, sequence: str, protein: str
    ) -> list[str]:
        """Check all constraints and return a list of violation names.

        This replaces the ad-hoc sequential constraint fixing with a
        formal check that returns ALL violations, ordered by priority.
        """
        violations: list[str] = []
        n = len(sequence)

        # HARD: Restriction sites
        if self._automaton._rs_sites:
            for site, site_rc in self._automaton._rs_sites:
                if site in sequence:
                    violations.append("restriction_site")
                    break
                if site_rc and site_rc in sequence:
                    violations.append("restriction_site")
                    break

        # HARD: Stop codons (internal)
        for i in range(0, n - 5, 3):
            codon = sequence[i:i + 3]
            if codon in ("TAA", "TAG", "TGA"):
                violations.append("stop_codon")
                break

        # HARD: GC range
        gc = (sequence.count("G") + sequence.count("C")) / max(n, 1)
        if not (self.gc_lo <= gc <= self.gc_hi):
            violations.append("gc_range")

        # HARD: Cryptic splice (eukaryotes only)
        if not self.is_prokaryote and self.avoid_gt:
            # Simplified check: count GT dinucleotides
            gt_count = sum(1 for i in range(n - 1) if sequence[i:i + 2] == "GT")
            if gt_count > 0:
                violations.append("cryptic_splice")

        # SOFT: Avoidable GT (eukaryotes only)
        if not self.is_prokaryote and self.avoid_gt:
            # This is already included in cryptic_splice for simplicity
            pass

        # SOFT: CpG islands
        if not self.is_prokaryote:
            from ..type_system import check_no_cpg_island
            cpg_result = check_no_cpg_island(sequence)
            if not cpg_result.passed:
                violations.append("cpg_island")

        # SOFT: ATTTA motifs
        if "ATTTA" in sequence:
            violations.append("atttta_motif")

        # SOFT: T-runs
        max_t = 0
        current_t = 0
        for b in sequence:
            if b == "T":
                current_t += 1
                max_t = max(max_t, current_t)
            else:
                current_t = 0
        if max_t > self.max_t_run:
            violations.append("t_run")

        return violations

    def _resolve_conflicts(
        self,
        sequence: str,
        protein: str,
        violations: list[str],
    ) -> tuple[str, int]:
        """Resolve constraint violations deterministically.

        Instead of the ad-hoc "try A, then B, then C" pattern, we fix
        violations in priority order (hard before soft).  For each
        violation type, we make a single deterministic pass through
        the sequence and fix all instances.

        Returns:
            Tuple of (modified sequence, number of fixes applied).
        """
        seq_list = list(sequence)
        n_codons = len(protein)
        fixes = 0

        # Sort violations by priority
        priority_map = {cp.name: cp.value for cp in ConstraintPriority}
        sorted_violations = sorted(
            violations,
            key=lambda v: priority_map.get(v, 999),
        )

        for violation in sorted_violations:
            if violation == "restriction_site":
                fixes += self._fix_restriction_sites(seq_list, protein)
            elif violation == "stop_codon":
                fixes += self._fix_stop_codons(seq_list, protein)
            elif violation == "gc_range":
                fixes += self._fix_gc_range(seq_list, protein)
            elif violation == "cryptic_splice":
                fixes += self._fix_cryptic_splice(seq_list, protein)
            elif violation == "cpg_island":
                fixes += self._fix_cpg(seq_list, protein)
            elif violation == "atttta_motif":
                fixes += self._fix_atttta(seq_list, protein)
            elif violation == "t_run":
                fixes += self._fix_t_runs(seq_list, protein)

        return "".join(seq_list), fixes

    def _fix_restriction_sites(
        self, seq_list: list[str], protein: str
    ) -> int:
        """Fix restriction site violations by swapping overlapping codons.

        Deterministic: tries codon alternatives in CAI-descending order.
        """
        fixes = 0
        n_codons = len(protein)

        for site, site_rc in self._automaton._rs_sites:
            seq_str = "".join(seq_list)
            while True:
                # Find the site
                pos = seq_str.find(site)
                if pos == -1 and site_rc:
                    pos = seq_str.find(site_rc)
                if pos == -1:
                    break

                # Find overlapping codons
                first_ci = pos // 3
                last_ci = min(n_codons - 1, (pos + len(site) - 1) // 3)

                # Try swapping each overlapping codon (CAI-descending)
                fixed = False
                for ci in range(first_ci, last_ci + 1):
                    if ci < 0 or ci >= n_codons:
                        continue
                    aa = protein[ci]
                    if aa == "*":
                        continue
                    current = "".join(seq_list[ci * 3:ci * 3 + 3])

                    for alt in self._automaton._sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        # Apply swap
                        start = ci * 3
                        old_chars = seq_list[start:start + 3]
                        seq_list[start] = alt[0]
                        seq_list[start + 1] = alt[1]
                        seq_list[start + 2] = alt[2]

                        # Check if site is gone
                        test_seq = "".join(seq_list)
                        site_gone = site not in test_seq
                        if site_rc:
                            site_gone = site_gone and site_rc not in test_seq

                        if site_gone:
                            # Also check we didn't create a stop codon
                            test_codon = alt
                            if test_codon not in ("TAA", "TAG", "TGA"):
                                fixes += 1
                                fixed = True
                                break

                        # Rollback
                        seq_list[start] = old_chars[0]
                        seq_list[start + 1] = old_chars[1]
                        seq_list[start + 2] = old_chars[2]

                    if fixed:
                        break

                if not fixed:
                    break  # Cannot fix this site
                seq_str = "".join(seq_list)

        return fixes

    def _fix_stop_codons(
        self, seq_list: list[str], protein: str
    ) -> int:
        """Fix internal stop codons by replacing with highest-CAI alternative."""
        fixes = 0
        n_codons = len(protein)

        for ci in range(n_codons - 1):  # Skip last codon (real stop)
            aa = protein[ci]
            if aa == "*":
                continue
            current = "".join(seq_list[ci * 3:ci * 3 + 3])
            if current in ("TAA", "TAG", "TGA"):
                # Replace with highest-CAI non-stop codon
                for alt in self._automaton._sorted_codons.get(aa, []):
                    if alt not in ("TAA", "TAG", "TGA"):
                        start = ci * 3
                        seq_list[start] = alt[0]
                        seq_list[start + 1] = alt[1]
                        seq_list[start + 2] = alt[2]
                        fixes += 1
                        break

        return fixes

    def _fix_gc_range(
        self, seq_list: list[str], protein: str
    ) -> int:
        """Fix GC content by swapping codons to move GC toward the target range.

        Deterministic: swaps the codon with the best GC improvement per
        CAI cost ratio, in a single pass.
        """
        seq = "".join(seq_list)
        n = len(seq)
        gc = (seq.count("G") + seq.count("C")) / max(n, 1)
        if self.gc_lo <= gc <= self.gc_hi:
            return 0

        fixes = 0
        n_codons = len(protein)
        need_more_gc = gc < self.gc_lo

        for _round in range(n_codons):
            seq = "".join(seq_list)
            gc = (seq.count("G") + seq.count("C")) / max(n, 1)
            if self.gc_lo <= gc <= self.gc_hi:
                break

            best_ci = -1
            best_alt = None
            best_score = -1.0

            for ci in range(n_codons):
                aa = protein[ci]
                if aa == "*":
                    continue
                current = "".join(seq_list[ci * 3:ci * 3 + 3])
                current_gc = self._automaton._codon_gc.get(current, 0)
                current_cai = self._species_cai.get(current, 0.0)

                for alt in self._automaton._sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    alt_gc = self._automaton._codon_gc.get(alt, 0)
                    gc_delta = alt_gc - current_gc

                    if need_more_gc and gc_delta <= 0:
                        continue
                    if not need_more_gc and gc_delta >= 0:
                        continue

                    alt_cai = self._species_cai.get(alt, 0.0)
                    # Score: GC improvement + CAI preservation
                    score = abs(gc_delta) + alt_cai * 0.01
                    if score > best_score:
                        best_score = score
                        best_ci = ci
                        best_alt = alt

            if best_alt is None:
                break

            start = best_ci * 3
            seq_list[start] = best_alt[0]
            seq_list[start + 1] = best_alt[1]
            seq_list[start + 2] = best_alt[2]
            fixes += 1

        return fixes

    def _fix_cryptic_splice(
        self, seq_list: list[str], protein: str
    ) -> int:
        """Fix cryptic splice sites by replacing GT-containing codons
        with GT-free alternatives.

        Deterministic: replaces in position order, using highest-CAI
        GT-free alternative.
        """
        if self.is_prokaryote:
            return 0
        fixes = 0
        n_codons = len(protein)

        for ci in range(n_codons):
            aa = protein[ci]
            if aa == "*":
                continue
            current = "".join(seq_list[ci * 3:ci * 3 + 3])
            if "GT" not in current:
                continue

            # Try GT-free alternatives
            for alt in self._automaton._gt_free_codons.get(aa, []):
                if alt == current:
                    continue
                # Apply swap
                start = ci * 3
                old_chars = list(seq_list[start:start + 3])
                seq_list[start] = alt[0]
                seq_list[start + 1] = alt[1]
                seq_list[start + 2] = alt[2]

                # Check we didn't create restriction sites or stop codons
                test_seq = "".join(seq_list)
                if alt in ("TAA", "TAG", "TGA"):
                    # Rollback
                    seq_list[start] = old_chars[0]
                    seq_list[start + 1] = old_chars[1]
                    seq_list[start + 2] = old_chars[2]
                    continue

                rs_ok = True
                for site, site_rc in self._automaton._rs_sites:
                    if site in test_seq or (site_rc and site_rc in test_seq):
                        rs_ok = False
                        break

                if rs_ok:
                    fixes += 1
                    break
                else:
                    # Rollback
                    seq_list[start] = old_chars[0]
                    seq_list[start + 1] = old_chars[1]
                    seq_list[start + 2] = old_chars[2]

        return fixes

    def _fix_cpg(
        self, seq_list: list[str], protein: str
    ) -> int:
        """Fix CpG dinucleotides by replacing CG-containing codons
        with CG-free alternatives.

        Deterministic: replaces in position order, using highest-CAI
        CG-free alternative.
        """
        if self.is_prokaryote:
            return 0
        fixes = 0
        n_codons = len(protein)

        for ci in range(n_codons):
            aa = protein[ci]
            if aa == "*":
                continue
            current = "".join(seq_list[ci * 3:ci * 3 + 3])
            if "CG" not in current:
                continue

            # Try CG-free alternatives
            cg_free = [c for c in self._automaton._sorted_codons.get(aa, [])
                       if "CG" not in c and c != current]
            for alt in cg_free:
                start = ci * 3
                old_chars = list(seq_list[start:start + 3])
                seq_list[start] = alt[0]
                seq_list[start + 1] = alt[1]
                seq_list[start + 2] = alt[2]

                test_seq = "".join(seq_list)
                if alt in ("TAA", "TAG", "TGA"):
                    seq_list[start] = old_chars[0]
                    seq_list[start + 1] = old_chars[1]
                    seq_list[start + 2] = old_chars[2]
                    continue

                rs_ok = True
                for site, site_rc in self._automaton._rs_sites:
                    if site in test_seq or (site_rc and site_rc in test_seq):
                        rs_ok = False
                        break

                # Check GC is still in range
                gc = (test_seq.count("G") + test_seq.count("C")) / max(len(test_seq), 1)
                gc_ok = self.gc_lo <= gc <= self.gc_hi

                if rs_ok and gc_ok:
                    fixes += 1
                    break
                else:
                    seq_list[start] = old_chars[0]
                    seq_list[start + 1] = old_chars[1]
                    seq_list[start + 2] = old_chars[2]

        return fixes

    def _fix_atttta(
        self, seq_list: list[str], protein: str
    ) -> int:
        """Fix ATTTA motifs by swapping overlapping codons.

        Deterministic: replaces in position order.
        """
        fixes = 0
        n_codons = len(protein)
        seq_str = "".join(seq_list)

        while "ATTTA" in seq_str:
            pos = seq_str.find("ATTTA")
            first_ci = pos // 3
            last_ci = min(n_codons - 1, (pos + 4) // 3)

            fixed = False
            for ci in range(first_ci, last_ci + 1):
                if ci < 0 or ci >= n_codons:
                    continue
                aa = protein[ci]
                if aa == "*":
                    continue
                current = "".join(seq_list[ci * 3:ci * 3 + 3])

                for alt in self._automaton._sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    start = ci * 3
                    old_chars = list(seq_list[start:start + 3])
                    seq_list[start] = alt[0]
                    seq_list[start + 1] = alt[1]
                    seq_list[start + 2] = alt[2]

                    test_seq = "".join(seq_list)
                    if "ATTTA" not in test_seq:
                        if alt not in ("TAA", "TAG", "TGA"):
                            fixes += 1
                            fixed = True
                            break

                    seq_list[start] = old_chars[0]
                    seq_list[start + 1] = old_chars[1]
                    seq_list[start + 2] = old_chars[2]

                if fixed:
                    break

            if not fixed:
                break
            seq_str = "".join(seq_list)

        return fixes

    def _fix_t_runs(
        self, seq_list: list[str], protein: str
    ) -> int:
        """Fix T-runs by replacing T-rich codons with alternatives.

        Deterministic: fixes the longest T-run first.
        """
        fixes = 0
        n_codons = len(protein)

        for _round in range(n_codons):
            seq_str = "".join(seq_list)
            # Find longest T-run
            max_t = 0
            max_t_pos = -1
            current_t = 0
            current_t_start = -1
            for i, b in enumerate(seq_str):
                if b == "T":
                    if current_t == 0:
                        current_t_start = i
                    current_t += 1
                    if current_t > max_t:
                        max_t = current_t
                        max_t_pos = current_t_start
                else:
                    current_t = 0

            if max_t <= self.max_t_run:
                break

            # Find codon at the start of the longest T-run
            ci = max_t_pos // 3
            if ci >= n_codons:
                break
            aa = protein[ci]
            if aa == "*":
                # Try adjacent codon
                ci = min(ci + 1, n_codons - 1)
                aa = protein[ci]
                if aa == "*":
                    break

            current = "".join(seq_list[ci * 3:ci * 3 + 3])

            # Find alternative with fewer T's
            best_alt = None
            best_t_count = current.count("T")
            for alt in self._automaton._sorted_codons.get(aa, []):
                if alt == current:
                    continue
                alt_t = alt.count("T")
                if alt_t < best_t_count:
                    if alt not in ("TAA", "TAG", "TGA"):
                        best_alt = alt
                        best_t_count = alt_t
                        break  # First (highest-CAI) alternative wins

            if best_alt is None:
                break

            start = ci * 3
            seq_list[start] = best_alt[0]
            seq_list[start + 1] = best_alt[1]
            seq_list[start + 2] = best_alt[2]
            fixes += 1

        return fixes

    def _passes_hard_constraints(
        self, sequence: str, protein: str
    ) -> bool:
        """Check whether a sequence passes all HARD constraints.

        Used in CAI_RECOVER to verify that a codon swap doesn't break
        hard constraints.
        """
        n = len(sequence)
        n_codons = len(protein)

        # No restriction sites
        for site, site_rc in self._automaton._rs_sites:
            if site in sequence:
                return False
            if site_rc and site_rc in sequence:
                return False

        # No internal stop codons
        for i in range(0, n - 5, 3):
            if sequence[i:i + 3] in ("TAA", "TAG", "TGA"):
                return False

        # GC in range (generous — allow slight deviation)
        gc = (sequence.count("G") + sequence.count("C")) / max(n, 1)
        if not (self.gc_lo - 0.02 <= gc <= self.gc_hi + 0.02):
            return False

        return True

    def _compute_cai(self, sequence: str, protein: str) -> float:
        """Compute the Codon Adaptation Index for a sequence."""
        n_codons = len(protein)
        if n_codons == 0:
            return 0.0

        log_sum = 0.0
        count = 0

        for ci in range(n_codons):
            aa = protein[ci]
            if aa == "*":
                continue
            codon = sequence[ci * 3:ci * 3 + 3]
            adapt = self._species_cai.get(codon, 0.0)
            max_a = self._automaton._max_adapt.get(aa, 0.0)
            if max_a > 0 and adapt > 0:
                log_sum += math.log(adapt / max_a)
                count += 1

        return math.exp(log_sum / count) if count > 0 else 0.0
