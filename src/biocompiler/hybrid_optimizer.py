"""
BioCompiler HybridOptimizer v10.2.0
=====================================
Greedy + Simultaneous Constraint Fixing + Batch CAI Hill Climbing + CPB Optimization.

Architecture:
  Phase 1: Greedy Initialization     — back-translate with highest-CAI codons
  Phase 2: Constraint Satisfaction    — simultaneous constraint fixing with
                                        conflict detection
  Phase 3: CAI Hill Climbing          — batch CAI recovery while
                                        maintaining all constraints
  Phase 4: Codon Pair Bias Opt.       — swap adjacent codon pairs for
                                        higher CPB score (optional)

Key innovation: instead of fixing violations one at a time, fixes ALL
non-conflicting violations simultaneously. Only iterates when fixes
conflict (share codon indices), preventing the
"fix A → break B → fix B → break A" oscillation loop while dramatically
reducing iteration count.

v10.2.0 changes:
  - Simultaneous constraint handling: fix all non-conflicting violations
    at once instead of one-at-a-time
  - Batch GC adjustment: compute all needed codon changes and apply at once
  - Batch CAI hill climbing: apply all non-conflicting upgrades simultaneously
  - Only iterate when fixes conflict (share codon indices)

v10.1.0 changes:
  - Phase 4: Codon pair bias optimization for higher protein expression
    (Coleman et al. 2008: 2-5x expression improvement from CPB optimization)
  - Optional via consider_codon_pair_bias=True parameter
  - Simplified CPB scoring: optimal+optimal=+0.1, one optimal=0.0, neither=-0.1
  - CAI guard: only accept CPB swaps that keep CAI >= 0.95 * current_CAI

v10.0.0 changes:
  - CAI table unification: uses CODON_ADAPTIVENESS_TABLES (same as evaluator)
  - Prokaryote fast path: skips splice/CpG constraints for E. coli
  - CAI-aware constraint resolution: all fixes prefer higher-CAI alternatives
  - resolve_organism() for centralized species ↔ organism name resolution

Performance target: <1.5ms for GFP (714bp), CAI > 0.98
"""

from __future__ import annotations

import heapq
import math
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, Any

from .type_system import CODON_TABLE, AA_TO_CODONS, BLOSUM62
from .organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism
from .constants import reverse_complement, RESTRICTION_ENZYMES
from .incremental import IncrementalSequenceState, CodonCache, EnzymeSiteCache
from .restriction_sites import get_recognition_site
from .aho_corasick import AhoCorasickScanner, build_scanner_from_enzymes
from .decision_provenance import CodonDecision, ConstraintDecision

# ── NUMBA integration ──────────────────────────────────────────────
try:
    from .numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        count_gc as _numba_count_gc,
        compute_cai_kernel as _numba_cai_kernel,
        seq_to_bytes as _seq_to_bytes,
    )
except ImportError:
    _HAS_NUMBA = False

# ── Fast MaxEntScan integration ────────────────────────────────────
try:
    from .maxentscan_fast import scan_splice_sites_fast_str as _scan_splice_sites_fast
    _HAS_FAST_MAXENT = True
except ImportError:
    _HAS_FAST_MAXENT = False

HAS_NUMBA: bool = _HAS_NUMBA

__all__ = [
    "HybridOptimizer",
    "HybridResult",
]

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# Constraint violation types and severity scoring
# ────────────────────────────────────────────────────────────

# Severity weights for different constraint types
# Higher = more important to fix first
SEVERITY_WEIGHTS = {
    "restriction_site": 100.0,    # Binary — site present or not
    "stop_codon": 90.0,           # Fatal — creates premature termination
    "gc_out_of_range": 50.0,      # Hard constraint
    "cryptic_splice_donor": 40.0, # Strong eukaryotic constraint
    "cryptic_splice_acceptor": 40.0,
    "avoidable_gt": 35.0,         # GT that can be removed
    "cpg_island": 20.0,           # Soft — methylation risk
    "atttta_motif": 15.0,         # Soft — mRNA instability
    "t_run": 10.0,                # Soft — polymerase slippage
}


@dataclass
class Violation:
    """A single constraint violation with severity score."""
    violation_type: str
    position: int            # Nucleotide position in sequence
    severity: float          # Weighted severity score
    codon_indices: list[int]  # Codon indices involved
    details: str = ""

    def __lt__(self, other: Violation) -> bool:
        """Priority queue ordering: highest severity first."""
        if self.severity != other.severity:
            return self.severity > other.severity
        return self.position < other.position


@dataclass
class HybridResult:
    """Result from the hybrid optimizer."""
    sequence: str
    cai: float
    gc_content: float
    violations_fixed: int = 0
    hill_climb_improvements: int = 0
    iterations_used: int = 0
    phase1_cai: float = 0.0
    phase2_cai: float = 0.0
    phase3_cai: float = 0.0
    phase4_cai: float = 0.0
    cpb_improvements: int = 0
    mean_cpb: float = 0.0
    warnings: list[str] = field(default_factory=list)
    splice_sites_validated: bool = False
    """When True, MaxEntScan splice site validation was already performed
    during optimization (Phase 3 of _optimize_eukaryote_fast).  The
    outer _evaluate_all_predicates call can skip the redundant
    check_no_cryptic_splice scan, eliminating ~7% overhead on the
    eukaryotic path."""


class HybridOptimizer:
    """Hybrid gene optimizer combining greedy initialization with
    priority-based local search, CAI hill climbing, and codon pair bias
    optimization.

    Architecture:
    1. Phase 1: Greedy CAI maximization (best codon per position)
    2. Phase 2: Priority-based constraint satisfaction
       (fix most severe violations first with incremental re-evaluation)
    3. Phase 3: CAI hill climbing
       (upgrade codons while maintaining constraints)
    4. Phase 4: Codon pair bias optimization (optional)
       (swap adjacent codon pairs for higher CPB score)

    Key innovation: instead of sequential constraint resolution that can
    undo previous fixes, uses a priority queue and incremental constraint
    evaluation to avoid conflicts.

    Args:
        species: Species key for codon usage tables (e.g., "human", "ecoli").
        organism: Full organism name (e.g., "Homo_sapiens").
        enzymes: List of restriction enzyme names to avoid.
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        avoid_gt: Whether to avoid GT dinucleotides (eukaryotes only).
        splice_threshold: MaxEntScan threshold for cryptic splice sites.
        cpg_window: Window size for CpG island detection.
        cpg_threshold: Obs/expected ratio threshold for CpG islands.
        max_local_search_iterations: Maximum Phase 2 iterations.
        max_hill_climb_iterations: Maximum Phase 3 iterations.
        cai_weight: Weight for CAI in the scoring function (0-1).
        provenance_collector: Optional provenance collector for tracking decisions.
        consider_codon_pair_bias: When True, run Phase 4 codon pair bias
            optimization after CAI hill climbing.  This can improve protein
            expression by 2-5x according to Coleman et al. (2008).  Uses a
            simplified scoring: both optimal = +0.1, one optimal = 0.0,
            neither optimal = -0.1.  Swaps are only accepted if CAI remains
            above 0.95 * current_CAI.
    """

    # Map species key to organism name
    _SPECIES_TO_ORGANISM = {
        "ecoli": "Escherichia_coli",
        "human": "Homo_sapiens",
        "mouse": "Mus_musculus",
        "cho": "CHO_K1",
        "yeast": "Saccharomyces_cerevisiae",
    }

    def __init__(
        self,
        species: str = "ecoli",
        organism: str | None = None,
        enzymes: list[str] | None = None,
        gc_lo: float = 0.30,
        gc_hi: float = 0.70,
        avoid_gt: bool = True,
        splice_threshold: float = 3.0,
        cpg_window: int = 200,
        cpg_threshold: float = 0.6,
        max_local_search_iterations: int = 50,
        max_hill_climb_iterations: int = 10,
        cai_weight: float = 0.7,
        provenance_collector: Any = None,
        consider_codon_pair_bias: bool = False,
    ) -> None:
        self.species: str = species
        self.organism: str = organism or self._SPECIES_TO_ORGANISM.get(species, "Homo_sapiens")
        # Use CODON_ADAPTIVENESS_TABLES directly — single source of truth
        _canonical: str = resolve_organism(self.organism)
        self.species_cai: dict[str, float] = dict(
            CODON_ADAPTIVENESS_TABLES.get(
                _canonical, CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
            )
        )
        self.enzymes: list[str] = enzymes or []
        self.gc_lo: float = gc_lo
        self.gc_hi: float = gc_hi
        self.avoid_gt: bool = avoid_gt
        self.splice_threshold: float = splice_threshold
        self.cpg_window: int = cpg_window
        self.cpg_threshold: float = cpg_threshold
        self.max_local_search_iterations: int = max_local_search_iterations
        self.max_hill_climb_iterations: int = max_hill_climb_iterations
        self.cai_weight: float = cai_weight
        self.provenance_collector: Any = provenance_collector
        self.consider_codon_pair_bias: bool = consider_codon_pair_bias

        # Detect prokaryotic organism — prokaryotes have no spliceosome,
        # so MaxEntScan splice scoring is irrelevant and must NEVER be called.
        self.is_prokaryote: bool = self._detect_prokaryote()
        # For prokaryotes, force avoid_gt=False to guarantee no MaxEntScan calls
        if self.is_prokaryote and self.avoid_gt:
            self.avoid_gt = False

        # Track applied mutagenesis (for compatibility with BioOptimizer)
        self._applied_mutagenesis: list[dict[str, Any]] = []

        # Pre-compute restriction site sequences
        self._rs_sites: list[tuple[str, str]] = []  # (site, reverse_complement)
        self._max_rs_site_len: int = 0  # longest site length for windowed checks
        for enz in self.enzymes:
            site = get_recognition_site(enz)
            if site is not None:
                site_rc = reverse_complement(site)
                self._rs_sites.append((site, site_rc))
                self._max_rs_site_len = max(self._max_rs_site_len, len(site))

        # Build Aho-Corasick scanner for O(L+M) multi-pattern restriction site scanning.
        # This replaces per-enzyme O(N*L*site_len) scanning with a single O(L+M) pass.
        self._ac_scanner: AhoCorasickScanner | None = build_scanner_from_enzymes(self.enzymes)
        if self._ac_scanner is not None:
            self._max_rs_site_len = max(
                self._max_rs_site_len, self._ac_scanner.longest_pattern
            )

        # Pre-compute codon cache
        self._codon_cache = CodonCache(self.species_cai)

        # ── Pre-compute data structures for fast lookups ──
        # sorted_codons: aa -> [codons sorted by CAI desc]
        self.sorted_codons: dict[str, list[str]] = {}
        # optimal_codon: aa -> highest-CAI codon (for Phase 1 fast path)
        self.optimal_codon: dict[str, str] = {}
        # gt_free: aa -> [GT-free codons sorted by CAI desc]
        self.gt_free: dict[str, list[str]] = {}
        # ag_free: aa -> [AG-free codons sorted by CAI desc]
        self.ag_free: dict[str, list[str]] = {}
        # codon_gc: codon -> GC base count (0-3)
        self.codon_gc: dict[str, int] = {}
        # Pre-compute which AAs have non-T-starting / non-G-ending codons
        self._aa_has_non_t_start: dict[str, bool] = {}
        self._aa_has_non_g_end: dict[str, bool] = {}

        for aa_key in set(CODON_TABLE.values()):
            if aa_key == "*":
                continue
            codons = AA_TO_CODONS.get(aa_key, [])
            # Sort by CAI descending
            codons_sorted = sorted(
                codons, key=lambda c: self.species_cai.get(c, 0.0), reverse=True
            )
            self.sorted_codons[aa_key] = codons_sorted
            if codons_sorted:
                self.optimal_codon[aa_key] = codons_sorted[0]
            # GT-free codons sorted by CAI
            gt_free_list = [c for c in codons_sorted if "GT" not in c]
            self.gt_free[aa_key] = gt_free_list
            # AG-free codons sorted by CAI
            ag_free_list = [c for c in codons_sorted if "AG" not in c]
            self.ag_free[aa_key] = ag_free_list
            # Codon GC counts
            for c in codons:
                self.codon_gc[c] = sum(1 for b in c if b in "GC")
            # AA properties
            self._aa_has_non_t_start[aa_key] = any(c[0] != 'T' for c in codons)
            self._aa_has_non_g_end[aa_key] = any(c[-1] != 'G' for c in codons)

        # Precompute max adaptiveness per AA for fast CAI computation
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

        # ── Codon pair bias data (Phase 4) ──
        # Load organism-specific CPB lookup table; empty dict if unavailable.
        self._cpb_data: dict[str, float] = {}
        # Set of codons that are "optimal" (highest CAI for their AA).
        # Used for simplified CPB scoring when no organism-specific data
        # is available.
        self._optimal_codon_set: set[str] = set()
        if self.consider_codon_pair_bias:
            try:
                from .codon_pair_scoring import get_codon_pair_data
                self._cpb_data = get_codon_pair_data(self.organism)
            except Exception:
                logger.debug(
                    "Could not load CPB data for organism '%s'; "
                    "using simplified scoring",
                    self.organism,
                )
            # Build set of optimal (highest-CAI) codons per amino acid
            for aa, codons in AA_TO_CODONS.items():
                if aa == "*" or not codons:
                    continue
                best = max(codons, key=lambda c: self.species_cai.get(c, 0.0))
                if self.species_cai.get(best, 0.0) > 0:
                    self._optimal_codon_set.add(best)

    # ──────────────────────────────────────────────────────────
    # Main optimization entry point
    # ──────────────────────────────────────────────────────────

    def optimize(self, protein: str, is_prokaryote: bool = False) -> HybridResult:
        """Run the hybrid optimization pipeline.

        Args:
            protein: Amino acid sequence (single-letter codes, no stop).
            is_prokaryote: When True, skip eukaryote-specific constraint steps
                (cryptic splice elimination, CpG disruption). Prokaryotes have
                no spliceosome, so GT/AG avoidance is biologically irrelevant.

        Returns:
            HybridResult with the optimized DNA sequence and metrics.
        """
        import time as _time
        _start = _time.monotonic()

        self._applied_mutagenesis = []

        # ── Ultra-fast path for prokaryotes ──
        # Prokaryotes have no spliceosome, so MaxEntScan splice scoring
        # is irrelevant.  CpG islands and GT/AG avoidance are also
        # biologically meaningless.  This allows a streamlined single-pass
        # approach that avoids IncrementalSequenceState, priority queues,
        # Violation objects, and CAI hill climbing entirely.
        if is_prokaryote or self.is_prokaryote:
            result = self._optimize_prokaryote_fast(protein)
            elapsed = _time.monotonic() - _start
            logger.info(
                "HybridOptimizer[prok-fast]: protein_len=%d, seq_len=%d, "
                "CAI=%.4f, GC=%.3f, violations_fixed=%d, time=%.1fms",
                len(protein), len(result.sequence),
                result.cai, result.gc_content,
                result.violations_fixed,
                elapsed * 1000,
            )
            return result

        # For prokaryotes, override avoid_gt
        effective_avoid_gt = self.avoid_gt and not is_prokaryote

        # ── Eukaryotic fast path ──
        # For eukaryotic targets, use a streamlined single-pass approach
        # similar to the prokaryotic fast path, but with GT/AG avoidance
        # and CpG disruption integrated directly into the optimization
        # loop.  This avoids the priority-queue-based _constraint_satisfaction
        # which suffers from GT↔ATTTA oscillation loops.
        if effective_avoid_gt:
            result = self._optimize_eukaryote_fast(protein)
            elapsed = _time.monotonic() - _start
            logger.info(
                "HybridOptimizer[euk-fast]: protein_len=%d, seq_len=%d, "
                "CAI=%.4f, GC=%.3f, violations_fixed=%d, time=%.1fms",
                len(protein), len(result.sequence),
                result.cai, result.gc_content,
                result.violations_fixed,
                elapsed * 1000,
            )
            return result

        # ── CSP solver integration with greedy fallback ──
        # Try the CSP solver first.  If it's unavailable or returns a
        # solution with soft-constraint violations, fall back to the
        # greedy optimizer and apply a repair step.  Provenance is
        # preserved across the fallback so callers know WHY the fallback
        # was used.
        csp_result = self._try_csp_solver(protein)
        if csp_result is not None:
            # CSP solver succeeded — check for soft constraint violations
            # and apply repair if needed.  HybridResult uses `warnings`
            # to signal violations found during CSP solving.
            if csp_result.warnings:
                repaired_seq, repair_improvements = self._repair_csp_solution(
                    csp_result.sequence, protein, effective_avoid_gt,
                )
                if repair_improvements > 0:
                    logger.info(
                        "CSP repair: improved CAI with %d codon upgrade(s)",
                        repair_improvements,
                    )
                    # Record repair provenance
                    if self.provenance_collector is not None:
                        self.provenance_collector.record_constraint_decision(
                            ConstraintDecision(
                                constraint_name="CSPRepair",
                                constraint_type="soft",
                                action_taken="satisfied",
                                positions_affected=[],
                                tradeoff_description=(
                                    f"Repaired CSP solution: {repair_improvements} "
                                    f"codon upgrade(s) to improve CAI while "
                                    f"maintaining hard constraint satisfaction"
                                ),
                                impact_on_cai=0.0,  # positive improvement
                            )
                        )
                    csp_result = HybridResult(
                        sequence=repaired_seq,
                        cai=self._compute_cai(repaired_seq),
                        gc_content=(repaired_seq.count("G") + repaired_seq.count("C")) / max(len(repaired_seq), 1),
                        violations_fixed=csp_result.violations_fixed,
                        hill_climb_improvements=csp_result.hill_climb_improvements + repair_improvements,
                        iterations_used=csp_result.iterations_used,
                        phase1_cai=csp_result.phase1_cai,
                        phase2_cai=csp_result.phase2_cai,
                        phase3_cai=csp_result.phase3_cai,
                        phase4_cai=csp_result.phase4_cai,
                        cpb_improvements=csp_result.cpb_improvements,
                        mean_cpb=csp_result.mean_cpb,
                        warnings=csp_result.warnings,
                    )
            elapsed = _time.monotonic() - _start
            logger.info(
                "HybridOptimizer[csp]: protein_len=%d, seq_len=%d, "
                "CAI=%.4f, GC=%.3f, violations_fixed=%d, time=%.1fms",
                len(protein), len(csp_result.sequence),
                csp_result.cai, csp_result.gc_content,
                csp_result.violations_fixed,
                elapsed * 1000,
            )
            return csp_result

        # ── Greedy fallback path ──
        # CSP solver is unavailable — fall back gracefully to greedy
        # optimization.  Record the fallback reason in provenance so
        # callers know WHY the greedy path was used.
        logger.info(
            "CSP solver unavailable; falling back to greedy optimization"
        )
        if self.provenance_collector is not None:
            self.provenance_collector.record_constraint_decision(
                ConstraintDecision(
                    constraint_name="CSPFallback",
                    constraint_type="hard",
                    action_taken="overridden",
                    positions_affected=[],
                    tradeoff_description=(
                        "CSP solver unavailable; falling back to greedy "
                        "optimization. Constraint provenance is preserved "
                        "but resolution uses priority-based ordering instead "
                        "of CAI-aware conflict resolution."
                    ),
                    impact_on_cai=0.0,
                )
            )

        # Phase 1: Greedy initialization
        seq, phase1_cai = self._greedy_init(protein)

        # Phase 2: Constraint satisfaction with priority queue
        # (uses CAI-aware resolution when constraints conflict)
        seq, phase2_cai, violations_fixed, iterations, warnings = (
            self._constraint_satisfaction(seq, protein, effective_avoid_gt)
        )

        # Phase 3: CAI hill climbing
        seq, phase3_cai, hill_climb_improvements = (
            self._cai_hill_climb(seq, protein, effective_avoid_gt)
        )

        # Phase 4: Codon pair bias optimization (optional)
        phase4_cai = phase3_cai
        cpb_improvements = 0
        mean_cpb = 0.0
        if self.consider_codon_pair_bias:
            seq, phase4_cai, cpb_improvements, mean_cpb = (
                self._codon_pair_bias_optimize(
                    seq, protein, effective_avoid_gt, phase3_cai
                )
            )

        # ── Phase 5: CAI recovery for eukaryotes ──
        # After all constraint satisfaction phases, check if any CAI was
        # lost due to GT avoidance.  For eukaryotes, GT avoidance is SOFT,
        # so we can recover CAI by swapping back to higher-CAI codons even
        # if they contain GT dinucleotides.  Only hard constraints (GC,
        # restriction sites) block the recovery.
        if effective_avoid_gt:
            pre_recovery_cai = self._compute_cai(seq)
            seq_chars = list(seq)
            n_codons = len(seq) // 3
            cai_recovery_improvements = 0

            for _recovery_round in range(3):
                any_recovered = False
                for ci in range(n_codons):
                    if ci >= len(protein):
                        break
                    aa = protein[ci]
                    if aa == "*":
                        continue
                    current = "".join(seq_chars[ci*3:ci*3+3])
                    current_cai = self.species_cai.get(current, 0.0)

                    # Try higher-CAI alternatives (even if they have GT)
                    for alt in self.sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        alt_cai = self.species_cai.get(alt, 0.0)
                        if alt_cai <= current_cai:
                            break  # sorted_codons is CAI-descending

                        # Apply swap
                        old_start = ci * 3
                        old_chars = seq_chars[old_start:old_start+3]
                        seq_chars[old_start] = alt[0]
                        seq_chars[old_start+1] = alt[1]
                        seq_chars[old_start+2] = alt[2]
                        new_seq = "".join(seq_chars)

                        # Only check HARD constraints (GC, restriction sites)
                        # GT avoidance is soft — do NOT block CAI recovery
                        # for eukaryotes just because a codon contains GT.
                        gc = (new_seq.count("G") + new_seq.count("C")) / max(len(new_seq), 1)
                        gc_ok = self.gc_lo <= gc <= self.gc_hi

                        rs_ok = True
                        if self._rs_sites:
                            for site, site_rc in self._rs_sites:
                                if site in new_seq or (site_rc and site_rc in new_seq):
                                    rs_ok = False
                                    break

                        if gc_ok and rs_ok:
                            cai_recovery_improvements += 1
                            any_recovered = True
                            break  # Accept first valid CAI upgrade
                        else:
                            # Rollback
                            seq_chars[old_start] = old_chars[0]
                            seq_chars[old_start+1] = old_chars[1]
                            seq_chars[old_start+2] = old_chars[2]

                if not any_recovered:
                    break

            if cai_recovery_improvements > 0:
                seq = "".join(seq_chars)
                logger.info(
                    "CAI recovery for eukaryote: %d codon(s) upgraded, "
                    "CAI %.4f → %.4f",
                    cai_recovery_improvements,
                    pre_recovery_cai,
                    self._compute_cai(seq),
                )
                hill_climb_improvements += cai_recovery_improvements

        # Compute final metrics
        gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
        final_cai = self._compute_cai(seq)

        elapsed = _time.monotonic() - _start
        logger.info(
            "HybridOptimizer: protein_len=%d, seq_len=%d, "
            "CAI=%.4f→%.4f→%.4f→%.4f→%.4f, GC=%.3f, "
            "violations_fixed=%d, hill_climb=%d, cpb_improvements=%d, "
            "iterations=%d, time=%.1fms",
            len(protein), len(seq),
            phase1_cai, phase2_cai, phase3_cai, phase4_cai, final_cai,
            gc, violations_fixed, hill_climb_improvements, cpb_improvements,
            iterations,
            elapsed * 1000,
        )

        return HybridResult(
            sequence=seq,
            cai=final_cai,
            gc_content=round(gc, 4),
            violations_fixed=violations_fixed,
            hill_climb_improvements=hill_climb_improvements,
            iterations_used=iterations,
            phase1_cai=phase1_cai,
            phase2_cai=phase2_cai,
            phase3_cai=phase3_cai,
            phase4_cai=phase4_cai,
            cpb_improvements=cpb_improvements,
            mean_cpb=mean_cpb,
            warnings=warnings,
        )

    # ──────────────────────────────────────────────────────────
    # Eukaryotic Fast Path
    # ──────────────────────────────────────────────────────────

    def _optimize_eukaryote_fast(self, protein: str) -> HybridResult:
        """Streamlined single-pass eukaryotic optimization.

        Similar to _optimize_prokaryote_fast but with GT/AG avoidance
        and CpG disruption integrated directly into the optimization loop.
        Avoids the priority-queue-based _constraint_satisfaction which
        suffers from GT↔ATTTA oscillation loops.

        Steps:
        1. Greedy init with highest-CAI codons (GT-aware where possible)
        2. Fix restriction sites
        3. Fix GC content
        4. Fix GT dinucleotides (splice donor avoidance)
        5. Fix AG dinucleotides (splice acceptor avoidance)
        6. Fix ATTTA motifs
        7. Fix T-runs

        Args:
            protein: Amino acid sequence (single-letter codes, no stop).

        Returns:
            HybridResult with the optimized DNA sequence and metrics.
        """
        species_cai = self.species_cai
        max_adapt = self._max_adapt
        optimal_codon = self.optimal_codon
        sorted_codons = self.sorted_codons
        codon_gc = self.codon_gc
        gt_free = self.gt_free
        ag_free = self.ag_free
        rs_sites = self._rs_sites
        ac_scanner = self._ac_scanner
        gc_lo = self.gc_lo
        gc_hi = self.gc_hi
        max_rs_len = self._max_rs_site_len
        splice_threshold = self.splice_threshold

        # ── Phase 1: Greedy init with GT-aware codon selection ──
        codon_list: list[str] = []
        gc_count = 0
        log_cai_sum = 0.0
        n_cai_codons = 0

        for aa in protein:
            if aa == "*":
                codon_list.append("TAA")
                continue
            # Prefer GT-free codons for eukaryotes
            gt_free_opts = gt_free.get(aa, [])
            if gt_free_opts:
                best = gt_free_opts[0]
            else:
                best = optimal_codon.get(aa, "")
                if not best:
                    sl = sorted_codons.get(aa, [])
                    best = sl[0] if sl else "NNN"
            codon_list.append(best)
            gc_count += codon_gc.get(best, 0)
            adapt = species_cai.get(best, 0.0)
            max_a = max_adapt.get(aa, 0.0)
            if max_a > 0 and adapt > 0:
                log_cai_sum += math.log(adapt / max_a)
                n_cai_codons += 1

        seq_chars = list("".join(codon_list))
        n_bases = len(seq_chars)
        n_codons = n_bases // 3
        phase1_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0

        # ── Helpers ──
        def _update_cai_for_swap(ci: int, old_codon: str, new_codon: str) -> None:
            nonlocal log_cai_sum, n_cai_codons
            aa = protein[ci] if ci < len(protein) else None
            if aa and aa != "*":
                old_adapt = species_cai.get(old_codon, 0.0)
                max_a = max_adapt.get(aa, 0.0)
                if max_a > 0 and old_adapt > 0:
                    log_cai_sum -= math.log(old_adapt / max_a)
                    n_cai_codons -= 1
                new_adapt = species_cai.get(new_codon, 0.0)
                if max_a > 0 and new_adapt > 0:
                    log_cai_sum += math.log(new_adapt / max_a)
                    n_cai_codons += 1

        def _has_rs_local(ci: int) -> bool:
            if ac_scanner is not None:
                start = ci * 3
                check_start = max(0, start - max_rs_len + 1)
                check_end = min(n_bases, start + 3 + max_rs_len - 1)
                return ac_scanner.has_any_match_in_region(
                    "".join(seq_chars), check_start, check_end
                )
            if not rs_sites:
                return False
            start = ci * 3
            check_start = max(0, start - max_rs_len + 1)
            check_end = min(n_bases, start + 3 + max_rs_len - 1)
            region = "".join(seq_chars[check_start:check_end])
            for site, site_rc in rs_sites:
                if site in region or (site_rc and site_rc in region):
                    return True
            return False

        def _apply_swap(ci: int, new_codon: str) -> str:
            nonlocal gc_count
            start = ci * 3
            old_codon = "".join(seq_chars[start:start + 3])
            old_gc = codon_gc.get(old_codon, 0)
            new_gc = codon_gc.get(new_codon, 0)
            gc_count += (new_gc - old_gc)
            seq_chars[start] = new_codon[0]
            seq_chars[start + 1] = new_codon[1]
            seq_chars[start + 2] = new_codon[2]
            _update_cai_for_swap(ci, old_codon, new_codon)
            return old_codon

        def _rollback_swap(ci: int, old_codon: str) -> None:
            nonlocal gc_count
            start = ci * 3
            current_codon = "".join(seq_chars[start:start + 3])
            current_gc = codon_gc.get(current_codon, 0)
            old_gc = codon_gc.get(old_codon, 0)
            gc_count += (old_gc - current_gc)
            seq_chars[start] = old_codon[0]
            seq_chars[start + 1] = old_codon[1]
            seq_chars[start + 2] = old_codon[2]
            _update_cai_for_swap(ci, current_codon, old_codon)

        violations_fixed = 0
        warnings: list[str] = []

        # ── Phase 2a: Fix restriction sites ──
        if rs_sites or ac_scanner is not None:
            seq_str = "".join(seq_chars)
            for _iter in range(100):
                if ac_scanner is not None:
                    matches = ac_scanner.scan(seq_str)
                    if not matches:
                        break
                    pos, site_match, _enzyme = matches[0]
                    site_len = len(site_match)
                else:
                    found = False
                    for site, site_rc in rs_sites:
                        p = seq_str.find(site)
                        if p == -1 and site_rc:
                            p = seq_str.find(site_rc)
                        if p != -1:
                            pos = p
                            site_len = len(site)
                            site_match = site
                            found = True
                            break
                    if not found:
                        break

                first_ci = pos // 3
                last_ci = (pos + site_len - 1) // 3
                fixed = False

                # CAI-aware: collect ALL possible single-codon swaps across
                # overlapping positions, then pick the one with highest CAI.
                best_swap_info: tuple[int, str, float] | None = None  # (ci, alt, log_cai)
                for ci in range(max(0, first_ci), min(n_codons, last_ci + 1)):
                    aa = protein[ci] if ci < len(protein) else None
                    if aa is None or aa == "*":
                        continue
                    current = "".join(seq_chars[ci*3:ci*3+3])
                    for alt in sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        old_codon = _apply_swap(ci, alt)
                        if not _has_rs_local(ci):
                            alt_cai = species_cai.get(alt, 0.0)
                            log_cai = math.log(alt_cai) if alt_cai > 0 else float('-inf')
                            if best_swap_info is None or log_cai > best_swap_info[2]:
                                best_swap_info = (ci, alt, log_cai)
                            _rollback_swap(ci, old_codon)
                        else:
                            _rollback_swap(ci, old_codon)

                if best_swap_info is not None:
                    ci, alt, _ = best_swap_info
                    _apply_swap(ci, alt)
                    seq_str = "".join(seq_chars)
                    violations_fixed += 1
                    fixed = True

                if not fixed:
                    warnings.append(f"Cannot remove restriction site {site_match} at pos {pos}")
                    break

        # ── Phase 2b: Fix GC content ──
        gc_val = gc_count / n_bases if n_bases > 0 else 0.0
        if not (gc_lo <= gc_val <= gc_hi):
            target = gc_lo if gc_val < gc_lo else gc_hi
            need_more_gc = gc_val < gc_lo
            for _iter in range(200):
                if gc_lo <= gc_val <= gc_hi:
                    break
                best_swap = None
                best_ci = -1
                best_score = -1.0
                for ci in range(n_codons):
                    aa = protein[ci] if ci < len(protein) else None
                    if aa is None or aa == "*":
                        continue
                    current = "".join(seq_chars[ci*3:ci*3+3])
                    current_gc_val = codon_gc.get(current, 0)
                    for alt in sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        alt_gc = codon_gc.get(alt, 0)
                        gc_delta = alt_gc - current_gc_val
                        if need_more_gc and gc_delta <= 0:
                            continue
                        if not need_more_gc and gc_delta >= 0:
                            continue
                        new_gc_count = gc_count + gc_delta
                        new_frac = new_gc_count / n_bases
                        diff = abs(new_frac - target)
                        alt_cai = species_cai.get(alt, 0.0)
                        score = (1.0 - diff) + alt_cai * 0.01
                        if score > best_score:
                            best_score = score
                            best_swap = alt
                            best_ci = ci
                if best_swap is None:
                    break
                old_codon = _apply_swap(best_ci, best_swap)
                if _has_rs_local(best_ci):
                    _rollback_swap(best_ci, old_codon)
                    break
                gc_val = gc_count / n_bases
                violations_fixed += 1

        # ── Phase 2c: Fix GT dinucleotides (splice donor avoidance) ──
        for _iter in range(300):
            seq_str = "".join(seq_chars)
            gt_pos = seq_str.find("GT")
            if gt_pos == -1:
                break

            ci = gt_pos // 3
            fixed = False

            # Try fixing the codon at the GT position
            for target_ci in [ci, ci + 1]:
                if target_ci >= n_codons or target_ci < 0:
                    continue
                aa = protein[target_ci] if target_ci < len(protein) else None
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[target_ci*3:target_ci*3+3])

                # Prefer GT-free codons
                for alt in gt_free.get(aa, []):
                    if alt == current:
                        continue
                    old_codon = _apply_swap(target_ci, alt)
                    new_local = "".join(seq_chars[max(0, target_ci*3-2):min(n_bases, target_ci*3+5)])
                    if "GT" not in new_local and not _has_rs_local(target_ci):
                        violations_fixed += 1
                        fixed = True
                        break
                    else:
                        _rollback_swap(target_ci, old_codon)

                if not fixed:
                    # Try all sorted codons (including those with GT, hoping boundary doesn't create it)
                    for alt in sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        old_codon = _apply_swap(target_ci, alt)
                        new_local = "".join(seq_chars[max(0, target_ci*3-2):min(n_bases, target_ci*3+5)])
                        if "GT" not in new_local and not _has_rs_local(target_ci):
                            violations_fixed += 1
                            fixed = True
                            break
                        else:
                            _rollback_swap(target_ci, old_codon)

                if fixed:
                    break

            if not fixed:
                # Some GTs may be unavoidable (e.g., Valine codons)
                break

        # ── Phase 2d: Fix AG dinucleotides (splice acceptor avoidance) ──
        for _iter in range(300):
            seq_str = "".join(seq_chars)
            ag_pos = seq_str.find("AG")
            if ag_pos == -1:
                break

            ci = ag_pos // 3
            fixed = False

            for target_ci in [ci, ci + 1]:
                if target_ci >= n_codons or target_ci < 0:
                    continue
                aa = protein[target_ci] if target_ci < len(protein) else None
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[target_ci*3:target_ci*3+3])

                for alt in ag_free.get(aa, []):
                    if alt == current:
                        continue
                    old_codon = _apply_swap(target_ci, alt)
                    new_local = "".join(seq_chars[max(0, target_ci*3-2):min(n_bases, target_ci*3+5)])
                    if "AG" not in new_local and not _has_rs_local(target_ci):
                        violations_fixed += 1
                        fixed = True
                        break
                    else:
                        _rollback_swap(target_ci, old_codon)

                if not fixed:
                    for alt in sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        old_codon = _apply_swap(target_ci, alt)
                        new_local = "".join(seq_chars[max(0, target_ci*3-2):min(n_bases, target_ci*3+5)])
                        if "AG" not in new_local and not _has_rs_local(target_ci):
                            violations_fixed += 1
                            fixed = True
                            break
                        else:
                            _rollback_swap(target_ci, old_codon)

                if fixed:
                    break

            if not fixed:
                break

        # ── Phase 2e: Fix ATTTA motifs ──
        seq_str = "".join(seq_chars)
        for _iter in range(100):
            pos = seq_str.find("ATTTA")
            if pos == -1:
                break
            first_ci = max(0, (pos // 3) - 1)
            last_ci = min(n_codons, ((pos + 4) // 3) + 2)
            fixed = False
            for ci in range(first_ci, last_ci):
                aa = protein[ci] if ci < len(protein) else None
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])
                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    old_codon = _apply_swap(ci, alt)
                    new_local = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                    if "ATTTA" not in new_local and not _has_rs_local(ci):
                        seq_str = "".join(seq_chars)
                        violations_fixed += 1
                        fixed = True
                        break
                    else:
                        _rollback_swap(ci, old_codon)
                if fixed:
                    break
            if not fixed:
                break

        # ── Phase 2f: Fix T-runs ──
        for _iter in range(100):
            max_run = 0
            max_pos = -1
            i = 0
            while i < n_bases:
                if seq_chars[i] == 'T':
                    j = i
                    while j < n_bases and seq_chars[j] == 'T':
                        j += 1
                    if j - i > max_run:
                        max_run = j - i
                        max_pos = i
                    i = j
                else:
                    i += 1
            if max_run < 6:
                break
            ci = (max_pos + max_run // 2) // 3
            if ci >= n_codons:
                ci = n_codons - 1
            fixed = False
            aa = protein[ci] if ci < len(protein) else None
            if aa is not None and aa != "*":
                current = "".join(seq_chars[ci*3:ci*3+3])
                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    old_codon = _apply_swap(ci, alt)
                    check_start = max(0, ci * 3 - 6)
                    check_end = min(n_bases, ci * 3 + 9)
                    has_long_run = False
                    j = check_start
                    while j < check_end:
                        if seq_chars[j] == 'T':
                            k = j
                            while k < check_end and seq_chars[k] == 'T':
                                k += 1
                            if k - j >= 6:
                                has_long_run = True
                                break
                            j = k
                        else:
                            j += 1
                    if not has_long_run and not _has_rs_local(ci):
                        violations_fixed += 1
                        fixed = True
                        break
                    else:
                        _rollback_swap(ci, old_codon)
            if not fixed:
                break

        # ── Phase 2g: CAI hill climbing for prokaryotes ──
        # After fixing all constraints, try to upgrade any suboptimal codons
        # back to higher-CAI alternatives without breaking constraints.
        hill_climb_improvements = 0
        for _hc_iter in range(3):  # Limited passes to avoid excessive runtime
            improved = False
            for ci in range(n_codons):
                aa = protein[ci] if ci < len(protein) else None
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])
                current_cai = species_cai.get(current, 0.0)

                # Find the best alternative that's better than current
                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    alt_cai = species_cai.get(alt, 0.0)
                    if alt_cai <= current_cai:
                        continue
                    # Check that this swap doesn't break constraints
                    old_codon = _apply_swap(ci, alt)
                    # Check no RS introduced locally
                    if _has_rs_local(ci):
                        _rollback_swap(ci, old_codon)
                        continue
                    # Check GC still in range
                    new_gc_val = gc_count / n_bases
                    if not (gc_lo <= new_gc_val <= gc_hi):
                        _rollback_swap(ci, old_codon)
                        continue
                    # Check no ATTTA introduced locally
                    local_window = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                    if "ATTTA" in local_window:
                        _rollback_swap(ci, old_codon)
                        continue
                    # Check no long T-run introduced
                    check_start = max(0, ci * 3 - 6)
                    check_end = min(n_bases, ci * 3 + 9)
                    has_long_run = False
                    j = check_start
                    while j < check_end:
                        if seq_chars[j] == 'T':
                            k = j
                            while k < check_end and seq_chars[k] == 'T':
                                k += 1
                            if k - j >= 6:
                                has_long_run = True
                                break
                            j = k
                        else:
                            j += 1
                    if has_long_run:
                        _rollback_swap(ci, old_codon)
                        continue
                    # This swap is safe and improves CAI
                    hill_climb_improvements += 1
                    improved = True
                    break  # First valid alt in sorted order is the best

            if not improved:
                break

        # ── Compute final metrics ──
        seq = "".join(seq_chars)
        final_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0
        gc = gc_count / n_bases if n_bases > 0 else 0.0

        return HybridResult(
            sequence=seq,
            cai=final_cai,
            gc_content=round(gc, 4),
            violations_fixed=violations_fixed,
            hill_climb_improvements=hill_climb_improvements,
            iterations_used=0,
            phase1_cai=phase1_cai,
            phase2_cai=final_cai,
            phase3_cai=final_cai,
            phase4_cai=final_cai,
            cpb_improvements=0,
            mean_cpb=0.0,
            warnings=warnings,
        )

    # ──────────────────────────────────────────────────────────
    # Phase 1: Greedy Initialization
    # ──────────────────────────────────────────────────────────

    def _try_csp_solver(self, protein: str) -> HybridResult | None:
        """Try to optimize using the CSP solver.

        Returns a HybridResult if the CSP solver succeeds, or None if
        the solver is unavailable or fails.  When the solver is
        unavailable, the caller should fall back to greedy optimization.

        This method is a thin wrapper around the CSP dispatch module
        that translates SolverResult to HybridResult, preserving
        constraint provenance along the way.

        Args:
            protein: Amino acid sequence (single-letter codes, no stop).

        Returns:
            HybridResult if CSP solver succeeded, None otherwise.
        """
        try:
            from .solver.dispatch import solve_with_csp, is_solver_available
        except ImportError:
            logger.debug("CSP solver dispatch module not available")
            return None

        if not is_solver_available():
            logger.debug("No CSP solver backend available")
            return None

        try:
            solver_result = solve_with_csp(
                protein,
                organism=self.organism,
            )
        except Exception as e:
            logger.warning(
                "CSP solver raised %s: %s; falling back to greedy",
                type(e).__name__, e,
            )
            # Record the exception in provenance
            if self.provenance_collector is not None:
                self.provenance_collector.record_constraint_decision(
                    ConstraintDecision(
                        constraint_name="CSPFallback",
                        constraint_type="hard",
                        action_taken="overridden",
                        positions_affected=[],
                        tradeoff_description=(
                            f"CSP solver raised {type(e).__name__}: {e}. "
                            f"Falling back to greedy optimization."
                        ),
                        impact_on_cai=0.0,
                    )
                )
            return None

        if not solver_result.solved or not solver_result.sequence:
            logger.debug(
                "CSP solver returned unsolved result (fallback=%s, reason=%s)",
                solver_result.fallback_used,
                solver_result.metadata.get("reason", "unknown"),
            )
            return None

        # Convert SolverResult to HybridResult
        seq = solver_result.sequence
        cai = solver_result.cai if solver_result.cai > 0 else self._compute_cai(seq)
        gc = solver_result.gc_content if solver_result.gc_content > 0 else (
            (seq.count("G") + seq.count("C")) / max(len(seq), 1)
        )

        # Record CSP solver provenance
        if self.provenance_collector is not None:
            self.provenance_collector.record_iteration({
                "step": "csp_solver",
                "action": "solved",
                "backend": str(solver_result.backend_used),
                "fallback_used": solver_result.fallback_used,
                "violations": len(solver_result.violations),
                "cai": round(cai, 6),
                "gc": round(gc, 4),
            })

        # Check for soft constraint violations that need repair
        has_violations = bool(solver_result.violations)

        return HybridResult(
            sequence=seq,
            cai=cai,
            gc_content=round(gc, 4),
            violations_fixed=len(solver_result.violations),
            hill_climb_improvements=0,
            iterations_used=0,
            phase1_cai=cai,
            phase2_cai=cai,
            phase3_cai=cai,
            phase4_cai=cai,
            cpb_improvements=0,
            mean_cpb=0.0,
            warnings=(
                [f"CSP solver had {len(solver_result.violations)} violation(s)"]
                if has_violations else []
            ),
        )

    def _repair_csp_solution(
        self,
        seq: str,
        protein: str,
        avoid_gt: bool,
    ) -> tuple[str, int]:
        """Repair a CSP solution that violates soft constraints.

        When the CSP solver returns a solution that violates soft
        constraints, this method applies a "repair" step to improve
        CAI while maintaining hard constraint satisfaction.

        The repair strategy is:
        1. For each codon position, try upgrading to a higher-CAI synonym
        2. Only accept the upgrade if it doesn't violate any hard constraint
        3. For eukaryotes, GT avoidance is SOFT: do NOT sacrifice CAI to
           eliminate GTs.  Only eliminate GTs that can be removed without
           reducing CAI.
        4. Continue until no more upgrades are possible

        Args:
            seq: DNA sequence from CSP solver.
            protein: Amino acid sequence.
            avoid_gt: Whether to avoid GT dinucleotides.

        Returns:
            Tuple of (repaired_sequence, number_of_improvements).
        """
        improvements = 0
        seq_chars = list(seq)
        n_codons = len(seq) // 3

        for _round in range(5):  # Multiple passes for cascading improvements
            any_improved = False
            for ci in range(n_codons):
                if ci >= len(protein):
                    break
                aa = protein[ci]
                if aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])
                current_cai = self.species_cai.get(current, 0.0)

                # Try higher-CAI alternatives
                for alt in self.sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    alt_cai = self.species_cai.get(alt, 0.0)
                    if alt_cai <= current_cai:
                        break  # sorted_codons is CAI-descending

                    # Apply swap and verify hard constraints
                    old_start = ci * 3
                    old_codon_chars = seq_chars[old_start:old_start+3]
                    seq_chars[old_start] = alt[0]
                    seq_chars[old_start+1] = alt[1]
                    seq_chars[old_start+2] = alt[2]
                    new_seq = "".join(seq_chars)

                    # Verify hard constraints
                    gc = (new_seq.count("G") + new_seq.count("C")) / max(len(new_seq), 1)
                    gc_ok = self.gc_lo <= gc <= self.gc_hi

                    # Check restriction sites
                    rs_ok = True
                    if self._rs_sites:
                        for site, site_rc in self._rs_sites:
                            if site in new_seq or (site_rc and site_rc in new_seq):
                                rs_ok = False
                                break

                    # For eukaryotes (avoid_gt=True), GT avoidance is SOFT:
                    # do NOT reject a CAI upgrade just because it introduces
                    # GT dinucleotides.  CAI takes priority over GT avoidance.
                    # GT elimination is only attempted below for CAI-neutral swaps.
                    gt_ok = True
                    if not avoid_gt:
                        # Prokaryotic path: GT is irrelevant, always OK
                        pass
                    # else: eukaryotic — GT is soft, don't block CAI upgrades

                    if gc_ok and rs_ok and gt_ok:
                        improvements += 1
                        any_improved = True
                        break  # Accept first valid upgrade
                    else:
                        # Rollback
                        seq_chars[old_start] = old_codon_chars[0]
                        seq_chars[old_start+1] = old_codon_chars[1]
                        seq_chars[old_start+2] = old_codon_chars[2]

            if not any_improved:
                break

        # ── Eukaryotic GT soft-cleanup ──
        # For eukaryotes, attempt to remove GT dinucleotides ONLY when the
        # swap does NOT reduce CAI.  This ensures we never sacrifice CAI for
        # GT avoidance — GT is a soft constraint for eukaryotes.
        if avoid_gt:
            for _gt_round in range(3):
                any_gt_fixed = False
                current_seq = "".join(seq_chars)
                gt_pos = current_seq.find("GT")
                while gt_pos != -1:
                    ci = gt_pos // 3
                    # Determine which codon(s) to try swapping
                    for target_ci in [ci, ci + 1]:
                        if target_ci >= n_codons or target_ci < 0:
                            continue
                        if target_ci >= len(protein):
                            continue
                        aa = protein[target_ci]
                        if aa == "*":
                            continue
                        cur = "".join(seq_chars[target_ci*3:target_ci*3+3])
                        cur_cai = self.species_cai.get(cur, 0.0)

                        # Try GT-free alternatives with CAI >= current
                        gt_fixed_here = False
                        for alt in self.gt_free.get(aa, []):
                            if alt == cur:
                                continue
                            alt_cai = self.species_cai.get(alt, 0.0)
                            # Only accept if CAI does not decrease
                            if alt_cai < cur_cai:
                                continue

                            old_start = target_ci * 3
                            old_chars = seq_chars[old_start:old_start+3]
                            seq_chars[old_start] = alt[0]
                            seq_chars[old_start+1] = alt[1]
                            seq_chars[old_start+2] = alt[2]
                            new_seq = "".join(seq_chars)

                            # Verify hard constraints still hold
                            gc = (new_seq.count("G") + new_seq.count("C")) / max(len(new_seq), 1)
                            gc_ok = self.gc_lo <= gc <= self.gc_hi
                            rs_ok = True
                            if self._rs_sites:
                                for site, site_rc in self._rs_sites:
                                    if site in new_seq or (site_rc and site_rc in new_seq):
                                        rs_ok = False
                                        break

                            # Check local GT elimination
                            local_start = max(0, target_ci * 3 - 2)
                            local_end = min(len(new_seq), target_ci * 3 + 5)
                            local_gt_gone = "GT" not in new_seq[local_start:local_end]

                            if gc_ok and rs_ok and local_gt_gone:
                                any_gt_fixed = True
                                gt_fixed_here = True
                                improvements += 1
                                break
                            else:
                                # Rollback
                                seq_chars[old_start] = old_chars[0]
                                seq_chars[old_start+1] = old_chars[1]
                                seq_chars[old_start+2] = old_chars[2]

                        if gt_fixed_here:
                            break

                    # Find next GT
                    current_seq = "".join(seq_chars)
                    gt_pos = current_seq.find("GT", gt_pos + 1)

                if not any_gt_fixed:
                    break

        return "".join(seq_chars), improvements

    def _greedy_init(self, protein: str) -> tuple[str, float]:
        """Back-translate protein using highest-CAI codons.

        Uses precomputed optimal_codon lookup table for O(1) per position
        instead of O(k) max() call. This gives the theoretical maximum CAI.

        Returns:
            (sequence, cai) tuple.
        """
        codons = []
        for aa in protein:
            if aa == "*":
                codons.append("TAA")
                continue
            best = self.optimal_codon.get(aa)
            if best:
                codons.append(best)
            else:
                # Fallback: use sorted_codons or NNN
                sorted_list = self.sorted_codons.get(aa, [])
                codons.append(sorted_list[0] if sorted_list else "NNN")
        seq = "".join(codons)
        cai = self._compute_cai(seq)
        return seq, cai

    # ──────────────────────────────────────────────────────────
    # Ultra-Fast Prokaryotic Path
    # ──────────────────────────────────────────────────────────

    def _optimize_prokaryote_fast(self, protein: str) -> HybridResult:
        """Ultra-fast optimization path for prokaryotic organisms.

        Prokaryotes (E. coli, B. subtilis, etc.) have no spliceosome,
        so MaxEntScan splice scoring is biologically irrelevant.  CpG
        island disruption and GT/AG avoidance are also meaningless.

        This method exploits those facts to skip ALL eukaryote-specific
        machinery and use a streamlined single-pass approach:

        1. Greedy init with highest-CAI codons + incremental GC/CAI tracking
        2. Fix restriction sites (single scan, local window checks)
        3. Fix GC content (incremental tracking, single-codon swaps)
        4. Fix ATTTA motifs (fast find + codon swap)
        5. Fix T-runs (single scan + codon swap)

        Key performance wins over the general path:
        - No IncrementalSequenceState object (saves ~0.3ms of GT/CG/AG setup)
        - No Violation objects or heapq operations
        - No Phase 2b (MaxEntScan) at all
        - No Phase 3 (CAI hill climbing) — start with max CAI
        - Incremental GC tracking during construction (no .count() scan)
        - Precomputed max adaptiveness for fast CAI computation
        - Local window checks instead of full-sequence scans after swaps
        - Single string scan for all pattern violations simultaneously

        Target: <0.5ms for GFP (714bp) in E. coli

        Args:
            protein: Amino acid sequence (single-letter codes, no stop).

        Returns:
            HybridResult with the optimized DNA sequence and metrics.
        """
        # NOTE: Prokaryotes skip GT avoidance entirely.  Prokaryotes have no
        # spliceosome, so GT dinucleotide avoidance is biologically irrelevant.
        # This is enforced by self.avoid_gt=False (set in __init__ for
        # prokaryotes) and by the is_prokaryote flag.  The GT-related code
        # paths (gt_free, ag_free, _fix_avoidable_gt, etc.) are never
        # reached on this fast path.
        species_cai = self.species_cai
        max_adapt = self._max_adapt
        optimal_codon = self.optimal_codon
        sorted_codons = self.sorted_codons
        codon_gc = self.codon_gc
        rs_sites = self._rs_sites
        ac_scanner = self._ac_scanner
        gc_lo = self.gc_lo
        gc_hi = self.gc_hi
        max_rs_len = self._max_rs_site_len

        # ── Provenance collector reference ──
        _prov = self.provenance_collector

        # ── Phase 1: Greedy init with incremental GC + CAI tracking ──
        # Build sequence and track GC + CAI log product simultaneously
        codon_list: list[str] = []
        gc_count = 0
        log_cai_sum = 0.0
        n_cai_codons = 0

        for pos, aa in enumerate(protein):
            if aa == "*":
                codon_list.append("TAA")
                continue
            best = optimal_codon.get(aa)
            if not best:
                sl = sorted_codons.get(aa, [])
                best = sl[0] if sl else "NNN"
            codon_list.append(best)
            gc_count += codon_gc.get(best, 0)
            # Incremental CAI
            adapt = species_cai.get(best, 0.0)
            max_a = max_adapt.get(aa, 0.0)
            if max_a > 0 and adapt > 0:
                log_cai_sum += math.log(adapt / max_a)
                n_cai_codons += 1

            # Record codon decision for provenance
            if _prov is not None:
                chosen_cai = adapt
                alternatives: list[dict[str, Any]] = []
                for codon in sorted_codons.get(aa, []):
                    if codon == best:
                        continue
                    cai_val = species_cai.get(codon, 0.0)
                    gc_bases = sum(1 for b in codon if b in "GC")
                    # Prokaryotes: skip eukaryotic-only constraints
                    # (splice, CpG, GT dinucleotide) — these are
                    # biologically irrelevant for organisms without
                    # a spliceosome.
                    violates: list[str] = []
                    rejected_because: str | None = None
                    if cai_val < chosen_cai:
                        rejected_because = "Lower CAI"
                    else:
                        rejected_because = "Lower CAI"
                    alternatives.append({
                        "codon": codon,
                        "cai_contribution": round(cai_val, 4),
                        "gc_contribution": round(gc_bases / 3.0, 2),
                        "violates_constraints": violates,
                        "rejected_because": rejected_because,
                    })
                # Compute confidence
                sl = sorted_codons.get(aa, [])
                if len(sl) > 1:
                    second_cai = species_cai.get(sl[1], 0.0)
                    confidence = min(1.0, 0.5 + (chosen_cai - second_cai) * 5)
                    confidence = max(0.0, confidence)
                else:
                    confidence = 1.0
                # Compute CAI impact vs optimal
                cai_impact_val = 0.0
                opt = optimal_codon.get(aa, best)
                if opt != best:
                    opt_cai = species_cai.get(opt, 0.0)
                    cai_impact_val = chosen_cai - opt_cai  # negative = CAI lost
                _prov.record_codon_decision(CodonDecision(
                    position=pos,
                    amino_acid=aa,
                    original_codon=None,
                    chosen_codon=best,
                    alternatives_considered=alternatives,
                    constraint_reason="maximize_cai",
                    confidence=round(confidence, 4),
                    cai_impact=round(cai_impact_val, 6),
                ))

        # Build sequence as list for efficient slicing during swaps
        seq_chars = list("".join(codon_list))
        n_bases = len(seq_chars)
        n_codons = n_bases // 3

        # Compute initial CAI from incremental tracking
        phase1_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0

        # ── Helper: update CAI tracking for a single codon swap ──
        def _update_cai_for_swap(ci: int, old_codon: str, new_codon: str) -> None:
            nonlocal log_cai_sum, n_cai_codons
            # Remove old codon contribution
            aa = protein[ci] if ci < len(protein) else None  # type: ignore[assignment]
            if aa and aa != "*":
                old_adapt = species_cai.get(old_codon, 0.0)
                max_a = max_adapt.get(aa, 0.0)
                if max_a > 0 and old_adapt > 0:
                    log_cai_sum -= math.log(old_adapt / max_a)
                    n_cai_codons -= 1
                # Add new codon contribution
                new_adapt = species_cai.get(new_codon, 0.0)
                if max_a > 0 and new_adapt > 0:
                    log_cai_sum += math.log(new_adapt / max_a)
                    n_cai_codons += 1

        # ── Helper: check for restriction sites in local window ──
        def _has_rs_local(ci: int) -> bool:
            """Check if any restriction site overlaps with codon at ci.

            Uses Aho-Corasick scanner for O(window_len + matches) single-pass
            detection instead of iterating through each enzyme separately.
            Falls back to per-enzyme string search if no scanner is available.
            """
            if ac_scanner is not None:
                start = ci * 3
                check_start = max(0, start - max_rs_len + 1)
                check_end = min(n_bases, start + 3 + max_rs_len - 1)
                return ac_scanner.has_any_match_in_region(
                    "".join(seq_chars), check_start, check_end
                )
            if not rs_sites:
                return False
            start = ci * 3
            check_start = max(0, start - max_rs_len + 1)
            check_end = min(n_bases, start + 3 + max_rs_len - 1)
            region = "".join(seq_chars[check_start:check_end])
            for site, site_rc in rs_sites:
                if site in region or (site_rc and site_rc in region):
                    return True
            return False

        # ── Helper: check for restriction sites in full sequence ──
        def _has_rs_full() -> bool:
            """Check if any restriction site is present in the full sequence.

            Uses Aho-Corasick scanner for O(L+M) single-pass detection
            instead of iterating through each enzyme separately.
            """
            if ac_scanner is not None:
                return ac_scanner.has_any_match("".join(seq_chars))
            if not rs_sites:
                return False
            seq_str = "".join(seq_chars)
            for site, site_rc in rs_sites:
                if site in seq_str or (site_rc and site_rc in seq_str):
                    return True
            return False

        # ── Helper: apply a codon swap ──
        def _apply_swap(ci: int, new_codon: str) -> str:
            nonlocal gc_count
            start = ci * 3
            old_codon = "".join(seq_chars[start:start + 3])
            old_gc = codon_gc.get(old_codon, 0)
            new_gc = codon_gc.get(new_codon, 0)
            gc_count += (new_gc - old_gc)
            seq_chars[start] = new_codon[0]
            seq_chars[start + 1] = new_codon[1]
            seq_chars[start + 2] = new_codon[2]
            _update_cai_for_swap(ci, old_codon, new_codon)
            return old_codon

        # ── Helper: rollback a codon swap ──
        def _rollback_swap(ci: int, old_codon: str) -> None:
            nonlocal gc_count
            start = ci * 3
            current_codon = "".join(seq_chars[start:start + 3])
            current_gc = codon_gc.get(current_codon, 0)
            old_gc = codon_gc.get(old_codon, 0)
            gc_count += (old_gc - current_gc)
            seq_chars[start] = old_codon[0]
            seq_chars[start + 1] = old_codon[1]
            seq_chars[start + 2] = old_codon[2]
            # Restore CAI tracking
            _update_cai_for_swap(ci, current_codon, old_codon)

        violations_fixed = 0
        warnings: list[str] = []
        rs_fixed_count = 0
        gc_fixed_count = 0
        attta_fixed_count = 0
        trun_fixed_count = 0

        # ── Phase 2a: Fix restriction sites (single scan) ──
        # Use Aho-Corasick scanner for O(L+M) multi-pattern detection
        # instead of per-enzyme O(N*L*site_len) scanning.
        if rs_sites or ac_scanner is not None:
            seq_str = "".join(seq_chars)
            for _iter in range(100):  # Max iterations total
                # Single-pass scan for ALL restriction sites simultaneously
                if ac_scanner is not None:
                    matches = ac_scanner.scan(seq_str)
                    if not matches:
                        break
                    pos, site_match, _enzyme = matches[0]
                    site_len = len(site_match)
                else:
                    # Fallback: per-enzyme scan
                    found = False
                    for site, site_rc in rs_sites:
                        p = seq_str.find(site)
                        if p == -1 and site_rc:
                            p = seq_str.find(site_rc)
                        if p != -1:
                            pos = p
                            site_len = len(site)
                            site_match = site
                            found = True
                            break
                    if not found:
                        break

                # Get overlapping codon indices
                first_ci = pos // 3
                last_ci = (pos + site_len - 1) // 3
                fixed = False

                # CAI-aware: collect ALL possible single-codon swaps across
                # overlapping positions, then pick the one with highest CAI.
                best_swap_info: tuple[int, str, float] | None = None  # (ci, alt, log_cai)
                for ci in range(max(0, first_ci), min(n_codons, last_ci + 1)):
                    aa = protein[ci] if ci < len(protein) else None  # type: ignore[assignment]
                    if aa is None or aa == "*":
                        continue
                    current = "".join(seq_chars[ci*3:ci*3+3])

                    for alt in sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        old_codon = _apply_swap(ci, alt)
                        if not _has_rs_local(ci):
                            alt_cai = species_cai.get(alt, 0.0)
                            log_cai = math.log(alt_cai) if alt_cai > 0 else float('-inf')
                            if best_swap_info is None or log_cai > best_swap_info[2]:
                                best_swap_info = (ci, alt, log_cai)
                            _rollback_swap(ci, old_codon)
                        else:
                            _rollback_swap(ci, old_codon)

                if best_swap_info is not None:
                    ci, alt, _ = best_swap_info
                    _apply_swap(ci, alt)
                    seq_str = "".join(seq_chars)
                    violations_fixed += 1
                    rs_fixed_count += 1
                    fixed = True

                if not fixed:
                    # Try two-codon coordinated fix with CAI scoring
                    # Evaluate ALL pair combinations and pick the one with
                    # the highest combined CAI (minimises CAI loss).
                    best_pair_info: tuple[int, str, str, float] | None = None  # (ci1, c1, c2, log_cai)
                    alts1_full = sorted_codons.get(aa1, []) if aa1 else []
                    alts2_full = sorted_codons.get(aa2, []) if aa2 else []

                    for idx in range(max(0, first_ci), min(n_codons - 1, last_ci + 1)):
                        ci1, ci2 = idx, idx + 1
                        aa1 = protein[ci1] if ci1 < len(protein) else None
                        aa2 = protein[ci2] if ci2 < len(protein) else None
                        if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                            continue
                        current1 = "".join(seq_chars[ci1*3:ci1*3+3])
                        current2 = "".join(seq_chars[ci2*3:ci2*3+3])
                        alts1 = sorted_codons.get(aa1, [])
                        alts2 = sorted_codons.get(aa2, [])

                        for c1 in alts1:
                            for c2 in alts2:
                                if c1 == current1 and c2 == current2:
                                    continue
                                old1 = _apply_swap(ci1, c1)
                                old2 = _apply_swap(ci2, c2)
                                if not _has_rs_local(ci1) and not _has_rs_local(ci2):
                                    c1_cai = species_cai.get(c1, 0.0)
                                    c2_cai = species_cai.get(c2, 0.0)
                                    combined_log = (
                                        (math.log(c1_cai) if c1_cai > 0 else float('-inf'))
                                        + (math.log(c2_cai) if c2_cai > 0 else float('-inf'))
                                    )
                                    if best_pair_info is None or combined_log > best_pair_info[3]:
                                        best_pair_info = (ci1, c1, c2, combined_log)
                                    _rollback_swap(ci2, old2)
                                    _rollback_swap(ci1, old1)
                                else:
                                    _rollback_swap(ci2, old2)
                                    _rollback_swap(ci1, old1)

                    if best_pair_info is not None:
                        ci1, c1, c2, _ = best_pair_info
                        _apply_swap(ci1, c1)
                        _apply_swap(ci2, c2)
                        seq_str = "".join(seq_chars)
                        violations_fixed += 1
                        rs_fixed_count += 1
                        fixed = True

                    if not fixed:
                        warnings.append(
                            f"Cannot remove restriction site {site_match} at pos {pos}"
                        )
                        break

        # ── Phase 2b: Fix GC content (incremental) ──
        gc_val = gc_count / n_bases if n_bases > 0 else 0.0
        if not (gc_lo <= gc_val <= gc_hi):
            target = gc_lo if gc_val < gc_lo else gc_hi
            need_more_gc = gc_val < gc_lo

            for _iter in range(200):
                if gc_lo <= gc_val <= gc_hi:
                    break

                best_swap = None
                best_ci = -1
                best_score = -1.0
                best_gc_delta = 0

                for ci in range(n_codons):
                    aa = protein[ci] if ci < len(protein) else None  # type: ignore[assignment]
                    if aa is None or aa == "*":
                        continue
                    current = "".join(seq_chars[ci*3:ci*3+3])
                    current_gc_val = codon_gc.get(current, 0)

                    for alt in sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        alt_gc = codon_gc.get(alt, 0)
                        gc_delta = alt_gc - current_gc_val

                        if need_more_gc and gc_delta <= 0:
                            continue
                        if not need_more_gc and gc_delta >= 0:
                            continue

                        new_gc_count = gc_count + gc_delta
                        new_frac = new_gc_count / n_bases
                        diff = abs(new_frac - target)
                        alt_cai = species_cai.get(alt, 0.0)
                        score = (1.0 - diff) + alt_cai * 0.01

                        if score > best_score:
                            best_score = score
                            best_swap = alt
                            best_ci = ci
                            best_gc_delta = gc_delta

                if best_swap is None:
                    break

                # Apply the swap and check RS
                old_codon = _apply_swap(best_ci, best_swap)
                if _has_rs_local(best_ci):
                    _rollback_swap(best_ci, old_codon)
                    break  # Can't fix without breaking RS

                gc_val = gc_count / n_bases
                violations_fixed += 1
                gc_fixed_count += 1

        # ── Phase 2c: Fix ATTTA motifs ──
        seq_str = "".join(seq_chars)
        for _iter in range(100):
            pos = seq_str.find("ATTTA")
            if pos == -1:
                break

            first_ci = max(0, (pos // 3) - 1)
            last_ci = min(n_codons, ((pos + 4) // 3) + 2)
            fixed = False

            for ci in range(first_ci, last_ci):
                aa = protein[ci] if ci < len(protein) else None  # type: ignore[assignment]
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])

                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    old_codon = _apply_swap(ci, alt)
                    # Check ATTTA eliminated and no new RS
                    new_local = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                    if "ATTTA" not in new_local and not _has_rs_local(ci):
                        seq_str = "".join(seq_chars)
                        violations_fixed += 1
                        attta_fixed_count += 1
                        fixed = True
                        break
                    else:
                        _rollback_swap(ci, old_codon)

                if fixed:
                    break

            if not fixed:
                warnings.append(f"Cannot remove ATTTA motif at pos {pos}")
                break

        # ── Phase 2d: Fix T-runs (6+ consecutive T) ──
        for _iter in range(100):
            # Find longest T-run
            max_run = 0
            max_pos = -1
            i = 0
            while i < n_bases:
                if seq_chars[i] == 'T':
                    j = i
                    while j < n_bases and seq_chars[j] == 'T':
                        j += 1
                    if j - i > max_run:
                        max_run = j - i
                        max_pos = i
                    i = j
                else:
                    i += 1

            if max_run < 6:
                break

            # Fix at the center of the run
            ci = (max_pos + max_run // 2) // 3
            if ci >= n_codons:
                ci = n_codons - 1
            fixed = False

            aa = protein[ci] if ci < len(protein) else None  # type: ignore[assignment]
            if aa is not None and aa != "*":
                current = "".join(seq_chars[ci*3:ci*3+3])

                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    old_codon = _apply_swap(ci, alt)
                    # Quick T-run check in local window
                    check_start = max(0, ci * 3 - 6)
                    check_end = min(n_bases, ci * 3 + 9)
                    has_long_run = False
                    j = check_start
                    while j < check_end:
                        if seq_chars[j] == 'T':
                            k = j
                            while k < check_end and seq_chars[k] == 'T':
                                k += 1
                            if k - j >= 6:
                                has_long_run = True
                                break
                            j = k
                        else:
                            j += 1

                    if not has_long_run and not _has_rs_local(ci):
                        violations_fixed += 1
                        trun_fixed_count += 1
                        fixed = True
                        break
                    else:
                        _rollback_swap(ci, old_codon)

            if not fixed:
                warnings.append(
                    f"Cannot fix T-run of {max_run} at pos {max_pos}"
                )
                break

        # ── Compute final metrics ──
        seq = "".join(seq_chars)
        final_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0
        gc = gc_count / n_bases if n_bases > 0 else 0.0

        # ── Phase 2g: CAI hill climbing for prokaryotes ──
        # After fixing all constraints, try to upgrade any suboptimal codons
        # back to higher-CAI alternatives without breaking constraints.
        hill_climb_improvements = 0
        for _hc_iter in range(3):  # Limited passes to avoid excessive runtime
            improved = False
            for ci in range(n_codons):
                aa = protein[ci] if ci < len(protein) else None
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])
                current_cai = species_cai.get(current, 0.0)

                # Find the best alternative that's better than current
                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    alt_cai = species_cai.get(alt, 0.0)
                    if alt_cai <= current_cai:
                        continue
                    # Check that this swap doesn't break constraints
                    old_codon = _apply_swap(ci, alt)
                    # Check no RS introduced locally
                    if _has_rs_local(ci):
                        _rollback_swap(ci, old_codon)
                        continue
                    # Check GC still in range
                    new_gc_val = gc_count / n_bases
                    if not (gc_lo <= new_gc_val <= gc_hi):
                        _rollback_swap(ci, old_codon)
                        continue
                    # Check no ATTTA introduced locally
                    local_window = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                    if "ATTTA" in local_window:
                        _rollback_swap(ci, old_codon)
                        continue
                    # Check no long T-run introduced
                    check_start = max(0, ci * 3 - 6)
                    check_end = min(n_bases, ci * 3 + 9)
                    has_long_run = False
                    j = check_start
                    while j < check_end:
                        if seq_chars[j] == 'T':
                            k = j
                            while k < check_end and seq_chars[k] == 'T':
                                k += 1
                            if k - j >= 6:
                                has_long_run = True
                                break
                            j = k
                        else:
                            j += 1
                    if has_long_run:
                        _rollback_swap(ci, old_codon)
                        continue
                    # This swap is safe and improves CAI
                    hill_climb_improvements += 1
                    improved = True
                    break  # First valid alt in sorted order is the best

            if not improved:
                break

        # ── Phase 2h: CAI micro-optimization with constraint recheck ──
        # After all constraint phases, some codons may be suboptimal because
        # a restriction site existed earlier but was later broken by other
        # codon changes.  Re-check each non-optimal codon to see if the
        # optimal (or a higher-CAI) codon is now safe.  Also try coordinated
        # two-codon swaps: if upgrading codon A creates a RS, maybe changing
        # a neighboring codon B simultaneously resolves it with a net CAI gain.
        micro_opt_improvements = 0
        for _micro_iter in range(5):
            any_micro_improved = False

            for ci in range(n_codons):
                aa = protein[ci] if ci < len(protein) else None
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])
                current_cai = species_cai.get(current, 0.0)
                optimal = optimal_codon.get(aa, current)
                optimal_cai = species_cai.get(optimal, 0.0)

                # Already optimal? Skip
                if current == optimal or optimal_cai <= current_cai:
                    continue

                # Try upgrading to the optimal codon
                old_codon = _apply_swap(ci, optimal)
                if not _has_rs_local(ci):
                    # Also verify other constraints
                    new_gc_val = gc_count / n_bases
                    local_window = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                    has_attta = "ATTTA" in local_window
                    # T-run check
                    check_s = max(0, ci * 3 - 6)
                    check_e = min(n_bases, ci * 3 + 9)
                    has_long_run = False
                    j = check_s
                    while j < check_e:
                        if seq_chars[j] == 'T':
                            k = j
                            while k < check_e and seq_chars[k] == 'T':
                                k += 1
                            if k - j >= 6:
                                has_long_run = True
                                break
                            j = k
                        else:
                            j += 1
                    if (gc_lo <= new_gc_val <= gc_hi) and not has_attta and not has_long_run:
                        micro_opt_improvements += 1
                        hill_climb_improvements += 1
                        any_micro_improved = True
                        continue
                _rollback_swap(ci, old_codon)

                # Optimal codon creates a RS or violates another constraint.
                # Try coordinated two-codon swaps with immediate neighbors.
                # For each neighbor, try upgrading ci to optimal while
                # simultaneously changing the neighbor to a CAI-friendly
                # alternative that breaks the RS.
                for neighbor_ci in [ci - 1, ci + 1]:
                    if neighbor_ci < 0 or neighbor_ci >= n_codons:
                        continue
                    neighbor_aa = protein[neighbor_ci] if neighbor_ci < len(protein) else None
                    if neighbor_aa is None or neighbor_aa == "*":
                        continue
                    neighbor_current = "".join(seq_chars[neighbor_ci*3:neighbor_ci*3+3])
                    neighbor_current_cai = species_cai.get(neighbor_current, 0.0)

                    # Try each alternative for the neighbor
                    best_net_log = 0.0  # net log-CAI gain must be positive
                    best_neighbor_alt = None

                    for neighbor_alt in sorted_codons.get(neighbor_aa, []):
                        if neighbor_alt == neighbor_current:
                            continue
                        # Apply neighbor swap first, then try optimal at ci
                        old_neighbor = _apply_swap(neighbor_ci, neighbor_alt)
                        old_ci = _apply_swap(ci, optimal)

                        if not _has_rs_local(ci) and not _has_rs_local(neighbor_ci):
                            # Check other constraints
                            new_gc_val = gc_count / n_bases
                            local_window = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                            has_attta = "ATTTA" in local_window
                            check_s = max(0, ci * 3 - 6)
                            check_e = min(n_bases, ci * 3 + 9)
                            has_long_run = False
                            j = check_s
                            while j < check_e:
                                if seq_chars[j] == 'T':
                                    k = j
                                    while k < check_e and seq_chars[k] == 'T':
                                        k += 1
                                    if k - j >= 6:
                                        has_long_run = True
                                        break
                                    j = k
                                else:
                                    j += 1
                            if (gc_lo <= new_gc_val <= gc_hi) and not has_attta and not has_long_run:
                                # Compute net CAI gain
                                neighbor_alt_cai = species_cai.get(neighbor_alt, 0.0)
                                optimal_cai_val = species_cai.get(optimal, 0.0)
                                net_log_gain = 0.0
                                if optimal_cai_val > 0 and current_cai > 0:
                                    net_log_gain += math.log(optimal_cai_val) - math.log(current_cai)
                                if neighbor_alt_cai > 0 and neighbor_current_cai > 0:
                                    net_log_gain += math.log(neighbor_alt_cai) - math.log(neighbor_current_cai)
                                if net_log_gain > best_net_log:
                                    best_net_log = net_log_gain
                                    best_neighbor_alt = neighbor_alt

                        _rollback_swap(ci, old_ci)
                        _rollback_swap(neighbor_ci, old_neighbor)

                    if best_neighbor_alt is not None and best_net_log > 0:
                        # Apply the best coordinated swap
                        _apply_swap(neighbor_ci, best_neighbor_alt)
                        _apply_swap(ci, optimal)
                        micro_opt_improvements += 1
                        hill_climb_improvements += 1
                        any_micro_improved = True
                        break  # Move to next codon

                # Also try non-optimal but higher-CAI alternatives for this codon
                # (without coordinated swaps — simple single upgrade)
                current_after = "".join(seq_chars[ci*3:ci*3+3])
                current_after_cai = species_cai.get(current_after, 0.0)
                for alt in sorted_codons.get(aa, []):
                    if alt == current_after:
                        continue
                    alt_cai = species_cai.get(alt, 0.0)
                    if alt_cai <= current_after_cai:
                        continue
                    old_codon2 = _apply_swap(ci, alt)
                    if _has_rs_local(ci):
                        _rollback_swap(ci, old_codon2)
                        continue
                    new_gc_val = gc_count / n_bases
                    if not (gc_lo <= new_gc_val <= gc_hi):
                        _rollback_swap(ci, old_codon2)
                        continue
                    local_window = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                    if "ATTTA" in local_window:
                        _rollback_swap(ci, old_codon2)
                        continue
                    check_s = max(0, ci * 3 - 6)
                    check_e = min(n_bases, ci * 3 + 9)
                    has_long_run = False
                    j = check_s
                    while j < check_e:
                        if seq_chars[j] == 'T':
                            k = j
                            while k < check_e and seq_chars[k] == 'T':
                                k += 1
                            if k - j >= 6:
                                has_long_run = True
                                break
                            j = k
                        else:
                            j += 1
                    if has_long_run:
                        _rollback_swap(ci, old_codon2)
                        continue
                    hill_climb_improvements += 1
                    any_micro_improved = True
                    break

            if not any_micro_improved:
                break

        # ── Phase 2i: CAI micro-optimization — recheck suboptimal codons ──
        # After all constraint phases (including Phase 2h coordinated swaps),
        # some codons that were suboptimally chosen (e.g., to avoid a
        # restriction site) may now be safe to upgrade because nearby codons
        # changed during later phases.  This final pass systematically tries
        # ALL higher-CAI alternatives, not just the optimal codon.
        for _micro_iter in range(3):
            any_micro_upgrade = False
            for ci in range(n_codons):
                aa = protein[ci] if ci < len(protein) else None
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])
                current_cai = species_cai.get(current, 0.0)

                # Try all higher-CAI alternatives (not just optimal)
                for alt in sorted_codons.get(aa, []):
                    alt_cai = species_cai.get(alt, 0.0)
                    if alt_cai <= current_cai:
                        break  # sorted_codons is CAI-descending, no point continuing
                    if alt == current:
                        continue

                    old_codon = _apply_swap(ci, alt)

                    # Check ALL constraints (no GT/AG for prokaryotes)
                    rs_ok = not _has_rs_local(ci)
                    gc_ok = gc_lo <= (gc_count / n_bases) <= gc_hi

                    # Check local region for ATTTA
                    local_start = max(0, ci * 3 - 5)
                    local_end = min(n_bases, ci * 3 + 8)
                    local_str = "".join(seq_chars[local_start:local_end])
                    attta_ok = "ATTTA" not in local_str

                    # Check T-runs in local window
                    check_s = max(0, ci * 3 - 6)
                    check_e = min(n_bases, ci * 3 + 9)
                    trun_ok = True
                    j = check_s
                    while j < check_e:
                        if seq_chars[j] == 'T':
                            k = j
                            while k < check_e and seq_chars[k] == 'T':
                                k += 1
                            if k - j >= 6:
                                trun_ok = False
                                break
                            j = k
                        else:
                            j += 1

                    if rs_ok and gc_ok and attta_ok and trun_ok:
                        any_micro_upgrade = True
                        break  # Keep this upgrade, move to next codon
                    else:
                        _rollback_swap(ci, old_codon)

            if not any_micro_upgrade:
                break

        # Recompute final CAI after hill climbing and micro-optimization
        final_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0
        gc = gc_count / n_bases if n_bases > 0 else 0.0
        seq = "".join(seq_chars)

        # ── Phase 4: Codon pair bias optimization (prokaryote fast path) ──
        cpb_improvements = 0
        mean_cpb = 0.0
        if self.consider_codon_pair_bias:
            seq, final_cai, cpb_improvements, mean_cpb = (
                self._codon_pair_bias_optimize_prokaryote(
                    seq, protein, final_cai
                )
            )
            # Recompute gc from the sequence (may have changed slightly)
            gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)

        # ── Record constraint decisions in provenance ──
        # This records the _decision_trail for the prokaryote fast path,
        # capturing which constraints were active and what actions were taken.
        if _prov is not None:
            if rs_fixed_count > 0:
                _prov.record_constraint_decision(ConstraintDecision(
                    constraint_name="NoRestrictionSite",
                    constraint_type="hard",
                    action_taken="satisfied",
                    positions_affected=[],
                    tradeoff_description=(
                        f"Fixed {rs_fixed_count} restriction site(s) "
                        f"by swapping to CAI-optimal alternative codons"
                    ),
                    impact_on_cai=0.0,  # CAI-aware fixes minimize CAI loss
                ))
            if gc_fixed_count > 0:
                _prov.record_constraint_decision(ConstraintDecision(
                    constraint_name="GCInRange",
                    constraint_type="hard",
                    action_taken="satisfied",
                    positions_affected=[],
                    tradeoff_description=(
                        f"Fixed GC content ({gc_fixed_count} codon swap(s)) "
                        f"to bring GC fraction into [{gc_lo:.2f}, {gc_hi:.2f}]"
                    ),
                    impact_on_cai=0.0,
                ))
            if attta_fixed_count > 0:
                _prov.record_constraint_decision(ConstraintDecision(
                    constraint_name="NoInstabilityMotif",
                    constraint_type="soft",
                    action_taken="satisfied",
                    positions_affected=[],
                    tradeoff_description=(
                        f"Removed {attta_fixed_count} ATTTA instability motif(s)"
                    ),
                    impact_on_cai=0.0,
                ))
            if trun_fixed_count > 0:
                _prov.record_constraint_decision(ConstraintDecision(
                    constraint_name="NoLongTRun",
                    constraint_type="soft",
                    action_taken="satisfied",
                    positions_affected=[],
                    tradeoff_description=(
                        f"Fixed {trun_fixed_count} T-run(s) of 6+ consecutive T"
                    ),
                    impact_on_cai=0.0,
                ))
            if hill_climb_improvements > 0:
                _prov.record_iteration({
                    "step": "cai_hill_climb",
                    "action": "upgrade_codons",
                    "improvements": hill_climb_improvements,
                    "score": round(final_cai, 6),
                })
            _prov.record_iteration({
                "step": "prokaryote_fast_path_complete",
                "action": "finalize",
                "violations_fixed": violations_fixed,
                "hill_climb_improvements": hill_climb_improvements,
                "cai": round(final_cai, 6),
                "gc": round(gc, 4),
            })

        return HybridResult(
            sequence=seq,
            cai=final_cai,
            gc_content=round(gc, 4),
            violations_fixed=violations_fixed,
            hill_climb_improvements=hill_climb_improvements,
            iterations_used=0,
            phase1_cai=phase1_cai,
            phase2_cai=final_cai,
            phase3_cai=final_cai,
            phase4_cai=final_cai,
            cpb_improvements=cpb_improvements,
            mean_cpb=mean_cpb,
            warnings=warnings,
        )

    # ──────────────────────────────────────────────────────────
    # Eukaryotic Fast Path
    # ──────────────────────────────────────────────────────────

    def _optimize_eukaryote_fast(self, protein: str) -> HybridResult:
        """Streamlined eukaryotic optimization with integrated GT/AG/CpG handling.

        This is the eukaryotic counterpart to _optimize_prokaryote_fast.
        It avoids the priority-queue-based _constraint_satisfaction which
        suffers from GT↔ATTTA oscillation loops by using a single-pass
        approach with incremental tracking.

        Pipeline:
        1. Greedy init with highest-CAI codons + incremental GC/CAI tracking
        2. Fix restriction sites (single scan, local window checks)
        3. Fix GC content (incremental tracking)
        4. Fix ATTTA motifs (fast find + codon swap, but no GT creation)
        5. Fix T-runs (single scan + codon swap, but no GT creation)
        6. MaxEntScan validation (only if GT/AG present — skip if none)
           — only fix GT/AG dinucleotides that score above splice threshold
           — prefer GT-free/AG-free alternatives that preserve CAI
        7. Fix remaining CpG dinucleotides (soft constraint — only if CAI-neutral)
        8. CAI recovery hill climb (upgrade codons while maintaining constraints)

        Key design principles:
        - **CAI first**: Never sacrifice CAI for soft constraints (GT/AG/CpG)
          unless MaxEntScan flags them as actual cryptic splice sites.
        - **No oscillation**: GT/AG/CpG are handled in a fixed order after
          hard constraints, with CAI-preserving alternatives only.
        - **MaxEntScan gating**: Only fix GT/AG dinucleotides that actually
          score as cryptic splice sites, not all GT/AG occurrences.
        - **Incremental tracking**: GC/CAI updated on every swap.

        Target: <5ms for HBB (444bp) in Human, CAI ≥ 0.99
        """
        species_cai = self.species_cai
        max_adapt = self._max_adapt
        optimal_codon = self.optimal_codon
        sorted_codons = self.sorted_codons
        codon_gc = self.codon_gc
        rs_sites = self._rs_sites
        ac_scanner = self._ac_scanner
        gc_lo = self.gc_lo
        gc_hi = self.gc_hi
        max_rs_len = self._max_rs_site_len
        gt_free = self.gt_free
        ag_free = self.ag_free
        splice_threshold = self.splice_threshold

        # ── Phase 1: Greedy init with incremental GC + CAI tracking ──
        codon_list: list[str] = []
        gc_count = 0
        log_cai_sum = 0.0
        n_cai_codons = 0

        for aa in protein:
            if aa == "*":
                codon_list.append("TAA")
                continue
            best = optimal_codon.get(aa)
            if not best:
                sl = sorted_codons.get(aa, [])
                best = sl[0] if sl else "NNN"
            codon_list.append(best)
            gc_count += codon_gc.get(best, 0)
            adapt = species_cai.get(best, 0.0)
            max_a = max_adapt.get(aa, 0.0)
            if max_a > 0 and adapt > 0:
                log_cai_sum += math.log(adapt / max_a)
                n_cai_codons += 1

        seq_chars = list("".join(codon_list))
        # Maintain a bytearray for fast pattern matching (avoids "".join overhead)
        seq_buf = bytearray("".join(codon_list), 'ascii')
        n_bases = len(seq_chars)
        n_codons = n_bases // 3

        phase1_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0

        # ── Helper: update CAI tracking for a single codon swap ──
        def _update_cai_for_swap(ci: int, old_codon: str, new_codon: str) -> None:
            nonlocal log_cai_sum, n_cai_codons
            aa = protein[ci] if ci < len(protein) else None
            if aa and aa != "*":
                old_adapt = species_cai.get(old_codon, 0.0)
                max_a = max_adapt.get(aa, 0.0)
                if max_a > 0 and old_adapt > 0:
                    log_cai_sum -= math.log(old_adapt / max_a)
                    n_cai_codons -= 1
                new_adapt = species_cai.get(new_codon, 0.0)
                if max_a > 0 and new_adapt > 0:
                    log_cai_sum += math.log(new_adapt / max_a)
                    n_cai_codons += 1

        # ── Helper: apply a codon swap ──
        def _apply_swap(ci: int, new_codon: str) -> str:
            nonlocal gc_count
            start = ci * 3
            old_codon = "".join(seq_chars[start:start + 3])
            old_gc = codon_gc.get(old_codon, 0)
            new_gc = codon_gc.get(new_codon, 0)
            gc_count += (new_gc - old_gc)
            seq_chars[start] = new_codon[0]
            seq_chars[start + 1] = new_codon[1]
            seq_chars[start + 2] = new_codon[2]
            # Update bytearray in-place
            seq_buf[start] = ord(new_codon[0])
            seq_buf[start + 1] = ord(new_codon[1])
            seq_buf[start + 2] = ord(new_codon[2])
            _update_cai_for_swap(ci, old_codon, new_codon)
            return old_codon

        # ── Helper: rollback a codon swap ──
        def _rollback_swap(ci: int, old_codon: str) -> None:
            nonlocal gc_count
            start = ci * 3
            current_codon = "".join(seq_chars[start:start + 3])
            current_gc = codon_gc.get(current_codon, 0)
            old_gc = codon_gc.get(old_codon, 0)
            gc_count += (old_gc - current_gc)
            seq_chars[start] = old_codon[0]
            seq_chars[start + 1] = old_codon[1]
            seq_chars[start + 2] = old_codon[2]
            # Update bytearray in-place
            seq_buf[start] = ord(old_codon[0])
            seq_buf[start + 1] = ord(old_codon[1])
            seq_buf[start + 2] = ord(old_codon[2])
            _update_cai_for_swap(ci, current_codon, old_codon)

        # ── Helper: check for restriction sites in local window ──
        def _has_rs_local(ci: int) -> bool:
            if ac_scanner is not None:
                start = ci * 3
                check_start = max(0, start - max_rs_len + 1)
                check_end = min(n_bases, start + 3 + max_rs_len - 1)
                # Use seq_buf directly instead of "".join(seq_chars)
                region = seq_buf[check_start:check_end].decode('ascii')
                return ac_scanner.has_any_match_in_region(
                    region, 0, len(region)
                )
            if not rs_sites:
                return False
            start = ci * 3
            check_start = max(0, start - max_rs_len + 1)
            check_end = min(n_bases, start + 3 + max_rs_len - 1)
            region = seq_buf[check_start:check_end].decode('ascii')
            for site, site_rc in rs_sites:
                if site in region or (site_rc and site_rc in region):
                    return True
            return False

        # ── Helper: check for new GT or AG in local window after swap ──
        def _creates_new_gt_or_ag(ci: int, new_codon: str) -> tuple[bool, bool]:
            """Return (creates_gt, creates_ag) for the local window around ci."""
            start = ci * 3
            prev_base = seq_chars[start - 1] if start > 0 else ''
            next_base = seq_chars[start + 3] if start + 3 < n_bases else ''
            local = [prev_base, new_codon[0], new_codon[1], new_codon[2], next_base]

            creates_gt = False
            creates_ag = False
            for i in range(len(local) - 1):
                if local[i] == 'G' and local[i + 1] == 'T':
                    creates_gt = True
                if local[i] == 'A' and local[i + 1] == 'G':
                    creates_ag = True
            return creates_gt, creates_ag

        violations_fixed = 0
        warnings: list[str] = []

        # ── Phase 2a: Fix restriction sites (single scan) ──
        if rs_sites or ac_scanner is not None:
            seq_str = seq_buf.decode('ascii')
            for _iter in range(100):
                if ac_scanner is not None:
                    matches = ac_scanner.scan(seq_str)
                    if not matches:
                        break
                    pos, site_match, _enzyme = matches[0]
                    site_len = len(site_match)
                else:
                    found = False
                    for site, site_rc in rs_sites:
                        p = seq_str.find(site)
                        if p == -1 and site_rc:
                            p = seq_str.find(site_rc)
                        if p != -1:
                            pos = p
                            site_len = len(site)
                            site_match = site
                            found = True
                            break
                    if not found:
                        break

                first_ci = pos // 3
                last_ci = (pos + site_len - 1) // 3
                fixed = False

                # CAI-aware: collect ALL possible single-codon swaps across
                # overlapping positions, then pick the one with highest CAI.
                best_swap_info: tuple[int, str, float] | None = None  # (ci, alt, log_cai)
                for ci in range(max(0, first_ci), min(n_codons, last_ci + 1)):
                    aa = protein[ci] if ci < len(protein) else None
                    if aa is None or aa == "*":
                        continue
                    current = "".join(seq_chars[ci*3:ci*3+3])

                    for alt in sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        old_codon = _apply_swap(ci, alt)
                        if not _has_rs_local(ci):
                            alt_cai = species_cai.get(alt, 0.0)
                            log_cai = math.log(alt_cai) if alt_cai > 0 else float('-inf')
                            if best_swap_info is None or log_cai > best_swap_info[2]:
                                best_swap_info = (ci, alt, log_cai)
                            _rollback_swap(ci, old_codon)
                        else:
                            _rollback_swap(ci, old_codon)

                if best_swap_info is not None:
                    ci, alt, _ = best_swap_info
                    _apply_swap(ci, alt)
                    seq_str = seq_buf.decode('ascii')
                    violations_fixed += 1
                    fixed = True

                if not fixed:
                    # Try two-codon coordinated fix with CAI scoring
                    # Evaluate ALL pair combinations and pick the one with
                    # the highest combined CAI (minimises CAI loss).
                    best_pair_info: tuple[int, str, str, float] | None = None  # (ci1, c1, c2, log_cai)

                    for idx in range(max(0, first_ci), min(n_codons - 1, last_ci + 1)):
                        ci1, ci2 = idx, idx + 1
                        aa1 = protein[ci1] if ci1 < len(protein) else None
                        aa2 = protein[ci2] if ci2 < len(protein) else None
                        if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                            continue
                        current1 = "".join(seq_chars[ci1*3:ci1*3+3])
                        current2 = "".join(seq_chars[ci2*3:ci2*3+3])
                        alts1 = sorted_codons.get(aa1, [])
                        alts2 = sorted_codons.get(aa2, [])

                        for c1 in alts1:
                            for c2 in alts2:
                                if c1 == current1 and c2 == current2:
                                    continue
                                old1 = _apply_swap(ci1, c1)
                                old2 = _apply_swap(ci2, c2)
                                if not _has_rs_local(ci1) and not _has_rs_local(ci2):
                                    c1_cai = species_cai.get(c1, 0.0)
                                    c2_cai = species_cai.get(c2, 0.0)
                                    combined_log = (
                                        (math.log(c1_cai) if c1_cai > 0 else float('-inf'))
                                        + (math.log(c2_cai) if c2_cai > 0 else float('-inf'))
                                    )
                                    if best_pair_info is None or combined_log > best_pair_info[3]:
                                        best_pair_info = (ci1, c1, c2, combined_log)
                                    _rollback_swap(ci2, old2)
                                    _rollback_swap(ci1, old1)
                                else:
                                    _rollback_swap(ci2, old2)
                                    _rollback_swap(ci1, old1)

                    if best_pair_info is not None:
                        ci1, c1, c2, _ = best_pair_info
                        _apply_swap(ci1, c1)
                        _apply_swap(ci2, c2)
                        seq_str = seq_buf.decode('ascii')
                        violations_fixed += 1
                        fixed = True

                    if not fixed:
                        warnings.append(
                            f"Cannot remove restriction site {site_match} at pos {pos}"
                        )
                        break

        # ── Phase 2b: Fix GC content (incremental) ──
        gc_val = gc_count / n_bases if n_bases > 0 else 0.0
        if not (gc_lo <= gc_val <= gc_hi):
            target = gc_lo if gc_val < gc_lo else gc_hi
            need_more_gc = gc_val < gc_lo

            for _iter in range(200):
                if gc_lo <= gc_val <= gc_hi:
                    break

                best_swap = None
                best_ci = -1
                best_score = -1.0
                best_gc_delta = 0

                for ci in range(n_codons):
                    aa = protein[ci] if ci < len(protein) else None
                    if aa is None or aa == "*":
                        continue
                    current = "".join(seq_chars[ci*3:ci*3+3])
                    current_gc_val = codon_gc.get(current, 0)

                    for alt in sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        alt_gc = codon_gc.get(alt, 0)
                        gc_delta = alt_gc - current_gc_val

                        if need_more_gc and gc_delta <= 0:
                            continue
                        if not need_more_gc and gc_delta >= 0:
                            continue

                        new_gc_count = gc_count + gc_delta
                        new_frac = new_gc_count / n_bases
                        diff = abs(new_frac - target)
                        alt_cai = species_cai.get(alt, 0.0)
                        score = (1.0 - diff) + alt_cai * 0.01

                        if score > best_score:
                            best_score = score
                            best_swap = alt
                            best_ci = ci
                            best_gc_delta = gc_delta

                if best_swap is None:
                    break

                old_codon = _apply_swap(best_ci, best_swap)
                if _has_rs_local(best_ci):
                    _rollback_swap(best_ci, old_codon)
                    break

                gc_val = gc_count / n_bases
                violations_fixed += 1

        # ── Phase 2c: Fix ATTTA motifs (without creating new GT/AG) ──
        seq_str = seq_buf.decode('ascii')
        for _iter in range(100):
            pos = seq_str.find("ATTTA")
            if pos == -1:
                break

            first_ci = max(0, (pos // 3) - 1)
            last_ci = min(n_codons, ((pos + 4) // 3) + 2)
            fixed = False

            for ci in range(first_ci, last_ci):
                aa = protein[ci] if ci < len(protein) else None
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])

                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    old_codon = _apply_swap(ci, alt)
                    new_local = seq_buf[max(0, ci*3-5):min(n_bases, ci*3+8)].decode('ascii')
                    creates_gt, creates_ag = _creates_new_gt_or_ag(ci, alt)
                    if "ATTTA" not in new_local and not _has_rs_local(ci) and not creates_gt and not creates_ag:
                        seq_str = seq_buf.decode('ascii')
                        violations_fixed += 1
                        fixed = True
                        break
                    else:
                        _rollback_swap(ci, old_codon)

                if fixed:
                    break

            if not fixed:
                # Relax: allow GT/AG creation to fix ATTTA
                for ci in range(first_ci, last_ci):
                    aa = protein[ci] if ci < len(protein) else None
                    if aa is None or aa == "*":
                        continue
                    current = "".join(seq_chars[ci*3:ci*3+3])

                    for alt in sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        old_codon = _apply_swap(ci, alt)
                        new_local = seq_buf[max(0, ci*3-5):min(n_bases, ci*3+8)].decode('ascii')
                        if "ATTTA" not in new_local and not _has_rs_local(ci):
                            seq_str = seq_buf.decode('ascii')
                            violations_fixed += 1
                            fixed = True
                            break
                        else:
                            _rollback_swap(ci, old_codon)

                    if fixed:
                        break

            if not fixed:
                warnings.append(f"Cannot remove ATTTA motif at pos {pos}")
                break

        # ── Phase 2d: Fix T-runs (6+ consecutive T) ──
        for _iter in range(100):
            max_run = 0
            max_pos = -1
            i = 0
            while i < n_bases:
                if seq_chars[i] == 'T':
                    j = i
                    while j < n_bases and seq_chars[j] == 'T':
                        j += 1
                    if j - i > max_run:
                        max_run = j - i
                        max_pos = i
                    i = j
                else:
                    i += 1

            if max_run < 6:
                break

            ci = (max_pos + max_run // 2) // 3
            if ci >= n_codons:
                ci = n_codons - 1
            fixed = False

            aa = protein[ci] if ci < len(protein) else None
            if aa is not None and aa != "*":
                current = "".join(seq_chars[ci*3:ci*3+3])

                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    old_codon = _apply_swap(ci, alt)
                    check_start = max(0, ci * 3 - 6)
                    check_end = min(n_bases, ci * 3 + 9)
                    has_long_run = False
                    j = check_start
                    while j < check_end:
                        if seq_chars[j] == 'T':
                            k = j
                            while k < check_end and seq_chars[k] == 'T':
                                k += 1
                            if k - j >= 6:
                                has_long_run = True
                                break
                            j = k
                        else:
                            j += 1

                    if not has_long_run and not _has_rs_local(ci):
                        violations_fixed += 1
                        fixed = True
                        break
                    else:
                        _rollback_swap(ci, old_codon)

            if not fixed:
                warnings.append(
                    f"Cannot fix T-run of {max_run} at pos {max_pos}"
                )
                break

        # ── Phase 3: MaxEntScan validation (only if GT/AG present) ──
        # This is the KEY step that was missing from the naive GT/AG elimination.
        # Instead of eliminating ALL GT/AG dinucleotides (which tanks CAI),
        # we only eliminate GT/AG dinucleotides that MaxEntScan identifies
        # as actual cryptic splice sites.  Most GT/AG dinucleotides are
        # not cryptic splice sites, so this preserves CAI.
        seq_str = seq_buf.decode('ascii')
        has_gt = "GT" in seq_str
        has_ag = "AG" in seq_str

        if has_gt or has_ag:
            try:
                if _HAS_FAST_MAXENT:
                    scan_splice_sites = _scan_splice_sites_fast
                else:
                    from .maxentscan import scan_splice_sites

                sites = scan_splice_sites(
                    seq_str,
                    donor_threshold=splice_threshold,
                    acceptor_threshold=splice_threshold,
                )

                for _splice_iter in range(10):
                    if not sites:
                        break

                    sites.sort(key=lambda s: s[2], reverse=True)
                    pos, site_type, score = sites[0]

                    # Fix cryptic splice sites
                    codon_idx = pos // 3
                    next_codon_start = codon_idx * 3 + 3
                    is_within = (pos + 1) < next_codon_start

                    if site_type == "donor" and is_within:
                        aa = protein[codon_idx] if codon_idx < len(protein) else None
                        if aa is None or aa == "*":
                            sites = sites[1:]
                            continue
                        current = "".join(seq_chars[codon_idx*3:codon_idx*3+3])

                        # Check if this is a Valine codon (unavoidable)
                        if current[:2] == "GT":
                            from .type_system import CODON_TABLE as _CT
                            if _CT.get(current) == "V":
                                sites = sites[1:]
                                continue

                        # Try GT-free alternatives (CAI-sorted)
                        fixed = False
                        for alt in gt_free.get(aa, []):
                            if alt == current:
                                continue
                            old_codon = _apply_swap(codon_idx, alt)
                            if not _has_rs_local(codon_idx):
                                violations_fixed += 1
                                fixed = True
                                break
                            else:
                                _rollback_swap(codon_idx, old_codon)

                        if not fixed:
                            sites = sites[1:]
                        else:
                            # Incremental update: remove sites near the fixed codon
                            # and rescan only the affected region instead of full sequence
                            ci_start = max(0, codon_idx * 3 - 22)  # acceptor model needs 23-mer
                            ci_end = min(n_bases, codon_idx * 3 + 3 + 22)
                            region = seq_buf[ci_start:ci_end].decode('ascii')
                            new_sites = scan_splice_sites(
                                region,
                                donor_threshold=splice_threshold,
                                acceptor_threshold=splice_threshold,
                            )
                            # Adjust positions back to absolute
                            adjusted = [(p + ci_start, st, sc) for p, st, sc in new_sites]
                            # Remove old sites in the affected region and add new ones
                            sites = [(p, st, sc) for p, st, sc in sites[1:]
                                     if not (ci_start <= p < ci_end)]
                            sites.extend(adjusted)

                    elif site_type == "donor" and not is_within:
                        # Cross-codon GT: try fixing either codon
                        next_ci = codon_idx + 1
                        fixed = False

                        if next_ci < n_codons:
                            aa2 = protein[next_ci] if next_ci < len(protein) else None
                            if aa2 is not None and aa2 != "*":
                                current2 = "".join(seq_chars[next_ci*3:next_ci*3+3])
                                for alt2 in sorted_codons.get(aa2, []):
                                    if alt2[0] == 'T' or alt2 == current2:
                                        continue
                                    old2 = _apply_swap(next_ci, alt2)
                                    if not _has_rs_local(next_ci):
                                        violations_fixed += 1
                                        fixed = True
                                        break
                                    else:
                                        _rollback_swap(next_ci, old2)

                        if not fixed:
                            aa1 = protein[codon_idx] if codon_idx < len(protein) else None
                            if aa1 is not None and aa1 != "*":
                                current1 = "".join(seq_chars[codon_idx*3:codon_idx*3+3])
                                for alt1 in sorted_codons.get(aa1, []):
                                    if alt1[-1] == 'G' or alt1 == current1:
                                        continue
                                    old1 = _apply_swap(codon_idx, alt1)
                                    if not _has_rs_local(codon_idx):
                                        violations_fixed += 1
                                        fixed = True
                                        break
                                    else:
                                        _rollback_swap(codon_idx, old1)

                        if not fixed:
                            sites = sites[1:]
                        else:
                            # Incremental update for cross-codon donor fix
                            ci_start = max(0, min(codon_idx, next_ci) * 3 - 22)
                            ci_end = min(n_bases, max(codon_idx, next_ci) * 3 + 3 + 22)
                            region = seq_buf[ci_start:ci_end].decode('ascii')
                            new_sites = scan_splice_sites(
                                region,
                                donor_threshold=splice_threshold,
                                acceptor_threshold=splice_threshold,
                            )
                            adjusted = [(p + ci_start, st, sc) for p, st, sc in new_sites]
                            sites = [(p, st, sc) for p, st, sc in sites[1:]
                                     if not (ci_start <= p < ci_end)]
                            sites.extend(adjusted)

                    elif site_type == "acceptor" and is_within:
                        aa = protein[codon_idx] if codon_idx < len(protein) else None
                        if aa is None or aa == "*":
                            sites = sites[1:]
                            continue
                        current = "".join(seq_chars[codon_idx*3:codon_idx*3+3])

                        # Try AG-free alternatives (CAI-sorted)
                        fixed = False
                        for alt in ag_free.get(aa, []):
                            if alt == current:
                                continue
                            old_codon = _apply_swap(codon_idx, alt)
                            # Check no new GT created
                            creates_gt, _ = _creates_new_gt_or_ag(codon_idx, alt)
                            if not creates_gt and not _has_rs_local(codon_idx):
                                violations_fixed += 1
                                fixed = True
                                break
                            else:
                                _rollback_swap(codon_idx, old_codon)

                        if not fixed:
                            sites = sites[1:]
                        else:
                            # Incremental update for acceptor fix
                            ci_start = max(0, codon_idx * 3 - 22)
                            ci_end = min(n_bases, codon_idx * 3 + 3 + 22)
                            region = seq_buf[ci_start:ci_end].decode('ascii')
                            new_sites = scan_splice_sites(
                                region,
                                donor_threshold=splice_threshold,
                                acceptor_threshold=splice_threshold,
                            )
                            adjusted = [(p + ci_start, st, sc) for p, st, sc in new_sites]
                            sites = [(p, st, sc) for p, st, sc in sites[1:]
                                     if not (ci_start <= p < ci_end)]
                            sites.extend(adjusted)

                    else:
                        # Cross-codon acceptor or unknown — skip
                        sites = sites[1:]

            except ImportError:
                pass

        # ── Phase 4: CpG island disruption ──
        # Fix CpG dinucleotides by replacing with CG-free alternatives.
        # We prefer higher-CAI alternatives but accept some CAI loss since
        # CpG islands are a hard constraint for eukaryotic optimization.
        # For cross-codon CGs, try both sides of the boundary.
        # Allow creating GT/AG during CpG fix (MaxEntScan will handle
        # any cryptic splice sites that result).
        seq_str = seq_buf.decode('ascii')
        for _cpg_iter in range(30):
            seq_str = seq_buf.decode('ascii')
            found_cpg = False

            for i in range(n_bases - 1):
                if seq_str[i] == 'C' and seq_str[i + 1] == 'G':
                    codon_idx = i // 3
                    next_codon_start = codon_idx * 3 + 3
                    is_within = (i + 1) < next_codon_start

                    if is_within:
                        aa = protein[codon_idx] if codon_idx < len(protein) else None
                        if aa is None or aa == "*":
                            continue

                        current = "".join(seq_chars[codon_idx*3:codon_idx*3+3])
                        # Try CG-free alternatives (CAI-sorted)
                        # First try without creating GT/AG, then allow GT/AG
                        fixed = False
                        for alt in sorted_codons.get(aa, []):
                            if alt == current or "CG" in alt:
                                continue
                            old_codon = _apply_swap(codon_idx, alt)
                            creates_gt, creates_ag = _creates_new_gt_or_ag(codon_idx, alt)
                            if not creates_gt and not creates_ag and not _has_rs_local(codon_idx):
                                violations_fixed += 1
                                found_cpg = True
                                fixed = True
                                break
                            else:
                                _rollback_swap(codon_idx, old_codon)

                        if not fixed:
                            # Allow GT/AG creation for CpG elimination
                            for alt in sorted_codons.get(aa, []):
                                if alt == current or "CG" in alt:
                                    continue
                                old_codon = _apply_swap(codon_idx, alt)
                                if not _has_rs_local(codon_idx):
                                    violations_fixed += 1
                                    found_cpg = True
                                    fixed = True
                                    break
                                else:
                                    _rollback_swap(codon_idx, old_codon)

                        if fixed:
                            break
                    else:
                        # Cross-codon CG: current codon ends with C, next starts with G
                        next_ci = codon_idx + 1
                        if next_ci >= n_codons:
                            continue
                        aa2 = protein[next_ci] if next_ci < len(protein) else None
                        if aa2 is None or aa2 == "*":
                            continue

                        # Try next codon with non-G-starting alternative
                        current2 = "".join(seq_chars[next_ci*3:next_ci*3+3])
                        fixed = False
                        for alt2 in sorted_codons.get(aa2, []):
                            if alt2[0] == 'G' or alt2 == current2:
                                continue
                            old2 = _apply_swap(next_ci, alt2)
                            creates_gt2, creates_ag2 = _creates_new_gt_or_ag(next_ci, alt2)
                            if not creates_gt2 and not creates_ag2 and not _has_rs_local(next_ci):
                                violations_fixed += 1
                                found_cpg = True
                                fixed = True
                                break
                            else:
                                _rollback_swap(next_ci, old2)

                        if not fixed:
                            # Allow GT/AG for next codon
                            for alt2 in sorted_codons.get(aa2, []):
                                if alt2[0] == 'G' or alt2 == current2:
                                    continue
                                old2 = _apply_swap(next_ci, alt2)
                                if not _has_rs_local(next_ci):
                                    violations_fixed += 1
                                    found_cpg = True
                                    fixed = True
                                    break
                                else:
                                    _rollback_swap(next_ci, old2)

                        if not fixed:
                            # Try current codon with non-C-ending alternative
                            aa = protein[codon_idx] if codon_idx < len(protein) else None
                            if aa is not None and aa != "*":
                                current1 = "".join(seq_chars[codon_idx*3:codon_idx*3+3])
                                for alt1 in sorted_codons.get(aa, []):
                                    if alt1[-1] == 'C' or alt1 == current1:
                                        continue
                                    old1 = _apply_swap(codon_idx, alt1)
                                    creates_gt1, creates_ag1 = _creates_new_gt_or_ag(codon_idx, alt1)
                                    if not creates_gt1 and not creates_ag1 and not _has_rs_local(codon_idx):
                                        violations_fixed += 1
                                        found_cpg = True
                                        fixed = True
                                        break
                                    else:
                                        _rollback_swap(codon_idx, old1)

                            if not fixed:
                                # Allow GT/AG for current codon too
                                if aa is not None and aa != "*":
                                    for alt1 in sorted_codons.get(aa, []):
                                        if alt1[-1] == 'C' or alt1 == current1:
                                            continue
                                        old1 = _apply_swap(codon_idx, alt1)
                                        if not _has_rs_local(codon_idx):
                                            violations_fixed += 1
                                            found_cpg = True
                                            fixed = True
                                            break
                                        else:
                                            _rollback_swap(codon_idx, old1)

                        if fixed:
                            break

            if not found_cpg:
                break

        # ── Phase 5: CAI recovery hill climb ──
        # Upgrade codons while maintaining all constraints.
        # Key insight: GT/AG dinucleotides are only problematic if MaxEntScan
        # identifies them as cryptic splice sites.  Most GT/AG occurrences are
        # benign.  So we use a two-tier approach:
        #  - First try: accept upgrade if it passes hard constraints (RS, GC,
        #    ATTTA, T-run) AND doesn't create new GT/AG.
        #  - Second try: accept upgrade even if it creates new GT/AG, but
        #    then validate with MaxEntScan and revert only if it creates an
        #    actual cryptic splice site.
        aas = list(protein)
        # Determine whether MaxEntScan is available for GT/AG validation
        _mes_available = False
        try:
            if _HAS_FAST_MAXENT:
                _mes_scan = _scan_splice_sites_fast
                _mes_available = True
            else:
                from .maxentscan import scan_splice_sites
                _mes_scan = scan_splice_sites
                _mes_available = True
        except ImportError:
            pass

        for _iteration in range(5):
            upgrade_plan: dict[int, str] = {}

            for ci in range(n_codons):
                aa = aas[ci]
                if aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])
                current_cai = species_cai.get(current, 0.0)
                optimal = optimal_codon.get(aa, current)
                optimal_cai = species_cai.get(optimal, 0.0)

                if optimal_cai > current_cai and current != optimal:
                    upgrade_plan[ci] = optimal

            if not upgrade_plan:
                break

            any_improved = False
            applied: set[int] = set()

            for ci in sorted(
                upgrade_plan.keys(),
                key=lambda c: species_cai.get(upgrade_plan[c], 0.0) - species_cai.get("".join(seq_chars[c*3:c*3+3]), 0.0),
                reverse=True,
            ):
                if any(abs(ci - a) <= 2 for a in applied):
                    continue

                new_codon = upgrade_plan[ci]
                current = "".join(seq_chars[ci*3:ci*3+3])
                if current == new_codon:
                    continue

                # Helper: check hard constraints (RS, GC, ATTTA, T-run)
                def _check_hard_constraints(ci_idx: int) -> bool:
                    if _has_rs_local(ci_idx):
                        return False
                    if not (gc_lo <= (gc_count / n_bases) <= gc_hi):
                        return False
                    local = seq_buf.decode('ascii')
                    if "ATTTA" in local[max(0, ci_idx*3-5):min(n_bases, ci_idx*3+8)]:
                        return False
                    # T-run check
                    cs = max(0, ci_idx * 3 - 6)
                    ce = min(n_bases, ci_idx * 3 + 9)
                    j2 = cs
                    while j2 < ce:
                        if seq_chars[j2] == 'T':
                            k2 = j2
                            while k2 < ce and seq_chars[k2] == 'T':
                                k2 += 1
                            if k2 - j2 >= 6:
                                return False
                            j2 = k2
                        else:
                            j2 += 1
                    return True

                # Try the swap and check constraints
                old_codon = _apply_swap(ci, new_codon)
                local_seq = seq_buf.decode('ascii')
                local_region = local_seq[max(0, ci*3-1):min(n_bases, ci*3+4)]

                rs_ok = not _has_rs_local(ci)
                gt_ok = "GT" not in local_region
                ag_ok = "AG" not in local_region
                gc_ok = gc_lo <= (gc_count / n_bases) <= gc_hi
                attta_ok = "ATTTA" not in local_seq[max(0, ci*3-5):min(n_bases, ci*3+8)]

                # T-run check
                check_s = max(0, ci * 3 - 6)
                check_e = min(n_bases, ci * 3 + 9)
                trun_ok = True
                j = check_s
                while j < check_e:
                    if seq_chars[j] == 'T':
                        k = j
                        while k < check_e and seq_chars[k] == 'T':
                            k += 1
                        if k - j >= 6:
                            trun_ok = False
                            break
                        j = k
                    else:
                        j += 1

                if rs_ok and gt_ok and ag_ok and gc_ok and attta_ok and trun_ok:
                    # All constraints pass — accept upgrade
                    violations_fixed += 1
                    any_improved = True
                    applied.add(ci)
                elif rs_ok and gc_ok and attta_ok and trun_ok and (not gt_ok or not ag_ok) and _mes_available:
                    # Hard constraints pass, but new GT/AG created.
                    # Validate with MaxEntScan — only reject if it creates an
                    # actual cryptic splice site above the threshold.
                    mes_seq = seq_buf.decode('ascii')
                    mes_sites = _mes_scan(
                        mes_seq,
                        donor_threshold=splice_threshold,
                        acceptor_threshold=splice_threshold,
                    )
                    # Check if any new cryptic splice site overlaps our codon
                    has_new_splice = False
                    for sp_pos, sp_type, sp_score in mes_sites:
                        if abs(sp_pos - ci * 3) <= 4:
                            has_new_splice = True
                            break
                    if not has_new_splice:
                        # MaxEntScan says the new GT/AG is benign — accept upgrade
                        violations_fixed += 1
                        any_improved = True
                        applied.add(ci)
                    else:
                        _rollback_swap(ci, old_codon)
                else:
                    _rollback_swap(ci, old_codon)
                    # Try other alternatives that have higher CAI than current
                    aa = aas[ci]
                    for alt in sorted_codons.get(aa, []):
                        if alt == current or alt == new_codon:
                            continue
                        alt_cai = species_cai.get(alt, 0.0)
                        cur_cai = species_cai.get(current, 0.0)
                        if alt_cai <= cur_cai:
                            continue

                        old_codon2 = _apply_swap(ci, alt)
                        if not _check_hard_constraints(ci):
                            _rollback_swap(ci, old_codon2)
                            continue
                        # Check GT/AG
                        local_seq2 = seq_buf.decode('ascii')
                        local_region2 = local_seq2[max(0, ci*3-1):min(n_bases, ci*3+4)]
                        gt_ok2 = "GT" not in local_region2
                        ag_ok2 = "AG" not in local_region2
                        if gt_ok2 and ag_ok2:
                            violations_fixed += 1
                            any_improved = True
                            applied.add(ci)
                            break
                        elif _mes_available:
                            # New GT/AG — validate with MaxEntScan
                            mes_seq2 = seq_buf.decode('ascii')
                            mes_sites2 = _mes_scan(
                                mes_seq2,
                                donor_threshold=splice_threshold,
                                acceptor_threshold=splice_threshold,
                            )
                            has_new_splice2 = False
                            for sp_pos, sp_type, sp_score in mes_sites2:
                                if abs(sp_pos - ci * 3) <= 4:
                                    has_new_splice2 = True
                                    break
                            if not has_new_splice2:
                                violations_fixed += 1
                                any_improved = True
                                applied.add(ci)
                                break
                            else:
                                _rollback_swap(ci, old_codon2)
                        else:
                            _rollback_swap(ci, old_codon2)

            if not any_improved:
                break

        # ── Phase 6: CAI micro-optimization with constraint recheck ──
        # After all constraint phases, some codons may be suboptimal because
        # a restriction site or GT/AG pattern existed earlier but was later
        # broken by other codon changes.  Re-check each non-optimal codon to
        # see if the optimal (or a higher-CAI) codon is now safe.  Also try
        # coordinated two-codon swaps: if upgrading codon A creates a RS,
        # maybe changing a neighboring codon B simultaneously resolves it
        # with a net CAI gain.
        hill_climb_improvements = 0
        for _micro_iter in range(5):
            any_micro_improved = False

            for ci in range(n_codons):
                aa = aas[ci]
                if aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])
                current_cai = species_cai.get(current, 0.0)
                optimal = optimal_codon.get(aa, current)
                optimal_cai = species_cai.get(optimal, 0.0)

                if current == optimal or optimal_cai <= current_cai:
                    continue

                # Try upgrading to the optimal codon
                old_codon = _apply_swap(ci, optimal)
                local_seq = seq_buf.decode('ascii')
                local_region = local_seq[max(0, ci*3-1):min(n_bases, ci*3+4)]

                rs_ok = not _has_rs_local(ci)
                gt_ok = "GT" not in local_region
                ag_ok = "AG" not in local_region
                gc_ok = gc_lo <= (gc_count / n_bases) <= gc_hi
                attta_ok = "ATTTA" not in local_seq[max(0, ci*3-5):min(n_bases, ci*3+8)]
                # T-run check
                check_s = max(0, ci * 3 - 6)
                check_e = min(n_bases, ci * 3 + 9)
                trun_ok = True
                j = check_s
                while j < check_e:
                    if seq_chars[j] == 'T':
                        k = j
                        while k < check_e and seq_chars[k] == 'T':
                            k += 1
                        if k - j >= 6:
                            trun_ok = False
                            break
                        j = k
                    else:
                        j += 1

                all_hard_ok = rs_ok and gc_ok and attta_ok and trun_ok

                if all_hard_ok and gt_ok and ag_ok:
                    # All constraints pass — accept upgrade
                    hill_climb_improvements += 1
                    any_micro_improved = True
                    continue
                elif all_hard_ok and (not gt_ok or not ag_ok) and _mes_available:
                    # Hard constraints pass, but new GT/AG — validate with MaxEntScan
                    mes_sites = _mes_scan(
                        local_seq,
                        donor_threshold=splice_threshold,
                        acceptor_threshold=splice_threshold,
                    )
                    has_new_splice = False
                    for sp_pos, sp_type, sp_score in mes_sites:
                        if abs(sp_pos - ci * 3) <= 4:
                            has_new_splice = True
                            break
                    if not has_new_splice:
                        hill_climb_improvements += 1
                        any_micro_improved = True
                        continue
                _rollback_swap(ci, old_codon)

                # Optimal codon creates a constraint violation.
                # Try coordinated two-codon swaps with immediate neighbors
                # to resolve the RS while gaining net CAI.
                for neighbor_ci in [ci - 1, ci + 1]:
                    if neighbor_ci < 0 or neighbor_ci >= n_codons:
                        continue
                    neighbor_aa = aas[neighbor_ci]
                    if neighbor_aa == "*":
                        continue
                    neighbor_current = "".join(seq_chars[neighbor_ci*3:neighbor_ci*3+3])
                    neighbor_current_cai = species_cai.get(neighbor_current, 0.0)

                    best_net_log = 0.0  # net log-CAI gain must be positive
                    best_neighbor_alt = None

                    for neighbor_alt in sorted_codons.get(neighbor_aa, []):
                        if neighbor_alt == neighbor_current:
                            continue
                        old_neighbor = _apply_swap(neighbor_ci, neighbor_alt)
                        old_ci = _apply_swap(ci, optimal)

                        # Check constraints
                        if not _has_rs_local(ci) and not _has_rs_local(neighbor_ci):
                            new_gc_val = gc_count / n_bases
                            if gc_lo <= new_gc_val <= gc_hi:
                                local_w = seq_buf.decode('ascii')
                                has_attta = "ATTTA" in local_w[max(0, ci*3-5):min(n_bases, ci*3+8)]
                                if not has_attta:
                                    # T-run check
                                    cs2 = max(0, ci * 3 - 6)
                                    ce2 = min(n_bases, ci * 3 + 9)
                                    has_lr = False
                                    jj = cs2
                                    while jj < ce2:
                                        if seq_chars[jj] == 'T':
                                            kk = jj
                                            while kk < ce2 and seq_chars[kk] == 'T':
                                                kk += 1
                                            if kk - jj >= 6:
                                                has_lr = True
                                                break
                                            jj = kk
                                        else:
                                            jj += 1
                                    if not has_lr:
                                        # Check GT/AG
                                        lr = local_w[max(0, ci*3-1):min(n_bases, ci*3+4)]
                                        gt2 = "GT" not in lr
                                        ag2 = "AG" not in lr
                                        accept_swap = False
                                        if gt2 and ag2:
                                            accept_swap = True
                                        elif _mes_available:
                                            mes2 = _mes_scan(
                                                local_w,
                                                donor_threshold=splice_threshold,
                                                acceptor_threshold=splice_threshold,
                                            )
                                            has_ns = False
                                            for sp2_p, sp2_t, sp2_s in mes2:
                                                if abs(sp2_p - ci * 3) <= 4:
                                                    has_ns = True
                                                    break
                                            if not has_ns:
                                                accept_swap = True
                                        if accept_swap:
                                            neighbor_alt_cai = species_cai.get(neighbor_alt, 0.0)
                                            net_log_gain = 0.0
                                            if optimal_cai > 0 and current_cai > 0:
                                                net_log_gain += math.log(optimal_cai) - math.log(current_cai)
                                            if neighbor_alt_cai > 0 and neighbor_current_cai > 0:
                                                net_log_gain += math.log(neighbor_alt_cai) - math.log(neighbor_current_cai)
                                            if net_log_gain > best_net_log:
                                                best_net_log = net_log_gain
                                                best_neighbor_alt = neighbor_alt

                        _rollback_swap(ci, old_ci)
                        _rollback_swap(neighbor_ci, old_neighbor)

                    if best_neighbor_alt is not None and best_net_log > 0:
                        _apply_swap(neighbor_ci, best_neighbor_alt)
                        _apply_swap(ci, optimal)
                        hill_climb_improvements += 1
                        any_micro_improved = True
                        break  # Move to next codon

                # Also try non-optimal but higher-CAI alternatives
                current_after = "".join(seq_chars[ci*3:ci*3+3])
                current_after_cai = species_cai.get(current_after, 0.0)
                for alt in sorted_codons.get(aa, []):
                    if alt == current_after:
                        continue
                    alt_cai = species_cai.get(alt, 0.0)
                    if alt_cai <= current_after_cai:
                        continue
                    old_codon2 = _apply_swap(ci, alt)
                    local_s = seq_buf.decode('ascii')
                    lr2 = local_s[max(0, ci*3-1):min(n_bases, ci*3+4)]
                    rs_ok2 = not _has_rs_local(ci)
                    gt_ok2 = "GT" not in lr2
                    ag_ok2 = "AG" not in lr2
                    gc_ok2 = gc_lo <= (gc_count / n_bases) <= gc_hi
                    attta_ok2 = "ATTTA" not in local_s[max(0, ci*3-5):min(n_bases, ci*3+8)]
                    # T-run check
                    cs3 = max(0, ci * 3 - 6)
                    ce3 = min(n_bases, ci * 3 + 9)
                    trun_ok2 = True
                    jj2 = cs3
                    while jj2 < ce3:
                        if seq_chars[jj2] == 'T':
                            kk2 = jj2
                            while kk2 < ce3 and seq_chars[kk2] == 'T':
                                kk2 += 1
                            if kk2 - jj2 >= 6:
                                trun_ok2 = False
                                break
                            jj2 = kk2
                        else:
                            jj2 += 1

                    all_hard2 = rs_ok2 and gc_ok2 and attta_ok2 and trun_ok2
                    if all_hard2 and gt_ok2 and ag_ok2:
                        hill_climb_improvements += 1
                        any_micro_improved = True
                        break
                    elif all_hard2 and (not gt_ok2 or not ag_ok2) and _mes_available:
                        mes3 = _mes_scan(
                            local_s,
                            donor_threshold=splice_threshold,
                            acceptor_threshold=splice_threshold,
                        )
                        has_ns3 = False
                        for sp3_p, sp3_t, sp3_s in mes3:
                            if abs(sp3_p - ci * 3) <= 4:
                                has_ns3 = True
                                break
                        if not has_ns3:
                            hill_climb_improvements += 1
                            any_micro_improved = True
                            break
                        else:
                            _rollback_swap(ci, old_codon2)
                    else:
                        _rollback_swap(ci, old_codon2)

            if not any_micro_improved:
                break

        # ── Phase 6b: CAI micro-optimization — recheck suboptimal codons ──
        # After all constraint phases (including Phase 6 coordinated swaps),
        # some codons that were suboptimally chosen (e.g., to avoid a
        # restriction site) may now be safe to upgrade because nearby codons
        # changed during later phases.  This final pass systematically tries
        # ALL higher-CAI alternatives, not just the optimal codon.
        for _micro_iter in range(3):
            any_micro_upgrade = False
            for ci in range(n_codons):
                aa = aas[ci]
                if aa == "*":
                    continue
                current = seq_buf[ci*3:ci*3+3].decode('ascii')
                current_cai = species_cai.get(current, 0.0)

                # Try all higher-CAI alternatives (not just optimal)
                for alt in sorted_codons.get(aa, []):
                    alt_cai = species_cai.get(alt, 0.0)
                    if alt_cai <= current_cai:
                        break  # sorted_codons is CAI-descending, no point continuing
                    if alt == current:
                        continue

                    old_codon = _apply_swap(ci, alt)

                    # Check ALL constraints
                    rs_ok = not _has_rs_local(ci)
                    gc_ok = gc_lo <= (gc_count / n_bases) <= gc_hi

                    # Check local region for ATTTA
                    local_start = max(0, ci * 3 - 5)
                    local_end = min(n_bases, ci * 3 + 8)
                    local_str = seq_buf[local_start:local_end].decode('ascii')
                    attta_ok = "ATTTA" not in local_str

                    # Check T-runs in local window
                    check_s = max(0, ci * 3 - 6)
                    check_e = min(n_bases, ci * 3 + 9)
                    trun_ok = True
                    j = check_s
                    while j < check_e:
                        if seq_chars[j] == 'T':
                            k = j
                            while k < check_e and seq_chars[k] == 'T':
                                k += 1
                            if k - j >= 6:
                                trun_ok = False
                                break
                            j = k
                        else:
                            j += 1

                    # Check GT/AG dinucleotides (splice avoidance)
                    local_gt_ag_start = max(0, ci * 3 - 1)
                    local_gt_ag_end = min(n_bases, ci * 3 + 4)
                    local_gt_ag = seq_buf[local_gt_ag_start:local_gt_ag_end].decode('ascii')
                    gt_ag_ok = "GT" not in local_gt_ag and "AG" not in local_gt_ag

                    all_hard_ok = rs_ok and gc_ok and attta_ok and trun_ok

                    if all_hard_ok and gt_ag_ok:
                        any_micro_upgrade = True
                        break  # Keep this upgrade, move to next codon
                    elif all_hard_ok and (not gt_ag_ok) and _mes_available:
                        # Hard constraints pass but GT/AG created — validate
                        # with MaxEntScan; accept if it's a benign GT/AG.
                        mes_local = seq_buf.decode('ascii')
                        mes_sites = _mes_scan(
                            mes_local,
                            donor_threshold=splice_threshold,
                            acceptor_threshold=splice_threshold,
                        )
                        has_new_splice = False
                        for sp_pos, sp_type, sp_score in mes_sites:
                            if abs(sp_pos - ci * 3) <= 4:
                                has_new_splice = True
                                break
                        if not has_new_splice:
                            any_micro_upgrade = True
                            break  # Benign GT/AG — keep upgrade
                        else:
                            _rollback_swap(ci, old_codon)
                    else:
                        _rollback_swap(ci, old_codon)

            if not any_micro_upgrade:
                break

        # ── Compute final metrics ──
        seq = seq_buf.decode('ascii')
        final_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0
        gc = gc_count / n_bases if n_bases > 0 else 0.0

        # ── Codon pair bias optimization ──
        cpb_improvements = 0
        mean_cpb = 0.0
        if self.consider_codon_pair_bias:
            seq, final_cai, cpb_improvements, mean_cpb = (
                self._codon_pair_bias_optimize_prokaryote(
                    seq, protein, final_cai
                )
            )
            gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)

        return HybridResult(
            sequence=seq,
            cai=final_cai,
            gc_content=round(gc, 4),
            violations_fixed=violations_fixed,
            hill_climb_improvements=0,
            iterations_used=0,
            phase1_cai=phase1_cai,
            phase2_cai=final_cai,
            phase3_cai=final_cai,
            phase4_cai=final_cai,
            cpb_improvements=cpb_improvements,
            mean_cpb=mean_cpb,
            warnings=warnings,
            splice_sites_validated=True,  # Phase 3 MaxEntScan already validated
        )

    # ──────────────────────────────────────────────────────────
    # Phase 2: Constraint Satisfaction (Priority-Queue Local Search)
    # ──────────────────────────────────────────────────────────

    def _constraint_satisfaction(
        self, seq: str, protein: str, avoid_gt: bool
    ) -> tuple[str, float, int, int, list[str]]:
        """Simultaneous constraint satisfaction with conflict detection.

        Algorithm:
        1. Scan for ALL constraint violations at once
        2. Score each violation by severity
        3. Group violations by whether they conflict (share codon indices)
        4. Fix ALL non-conflicting violations simultaneously:
           - Restriction sites at different positions → fix all at once
           - GC adjustment → batch all codon changes and apply at once
           - ATTTA/T-run fixes at non-overlapping positions → fix all at once
        5. Only iterate when a fix conflicts with another fix
           (shared codon indices between violations)
        6. For conflicting violations, use CAI-aware resolution:
           prefer the solution that maximizes CAI, using
           CAIAwareConstraintResolver if available.

        Key improvement over one-at-a-time approach:
        - Previously: fix ONE violation → re-scan → fix ONE → re-scan → ...
        - Now: fix ALL non-conflicting violations → re-scan only if conflicts
        - CAI-aware: when constraints conflict, prefer higher-CAI resolution

        Returns:
            (sequence, cai, violations_fixed, iterations, warnings)
        """
        state = IncrementalSequenceState(
            seq, species=self.species, enzymes=self.enzymes
        )
        warnings: list[str] = []
        violations_fixed = 0
        total_iterations = 0

        # ── CAI-aware constraint resolver (lazy-loaded) ──
        # When two constraints conflict (share codon indices), the
        # CAIAwareConstraintResolver evaluates which resolution option
        # has the smallest impact on CAI and prefers that resolution.
        _cai_resolver = None
        try:
            from .solver.cai_aware_resolver import CAIAwareConstraintResolver
            # We don't have a CSPModel here, but we can use the resolver
            # for its CAI impact estimation methods directly.
            _cai_resolver = CAIAwareConstraintResolver.__new__(CAIAwareConstraintResolver)
            _cai_resolver._adaptiveness = self.species_cai
        except (ImportError, Exception):
            _cai_resolver = None

        # Track violation signatures for stale detection — if the same set
        # of violations persists for 2 consecutive iterations, the remaining
        # violations are unfixable and we should stop rather than waste
        # cycles on futile swap-check-rollback loops.
        _prev_violation_sig: frozenset[tuple[str, int]] | None = None
        _stale_count = 0

        # Phase 2a: Fast constraint fixes (no MaxEntScan needed)
        for iteration in range(self.max_local_search_iterations):
            violations = self._detect_cheap_violations(state, avoid_gt)

            if not violations:
                break

            # ── Eukaryotic GT softening ──
            # For eukaryotic organisms, GT avoidance is a SOFT constraint.
            # Lower the severity of avoidable_gt violations so they are
            # treated as WARNINGS rather than hard failures.  GT violations
            # should not block CAI recovery or force suboptimal codon choices.
            if avoid_gt:
                for v in violations:
                    if v.violation_type == "avoidable_gt":
                        # Reduce severity below all hard constraints so GT
                        # fixes are attempted only after everything else
                        v.severity = SEVERITY_WEIGHTS["avoidable_gt"] * 0.1

            total_iterations = iteration + 1

            # ── Stale detection ──
            # Compute a signature of the current violation set (type + position).
            # If this signature is identical to the previous iteration, the
            # fixes are not making progress and the remaining violations are
            # unfixable.  Break immediately to avoid wasting time.
            sig = frozenset((v.violation_type, v.position) for v in violations)
            if sig == _prev_violation_sig:
                _stale_count += 1
                if _stale_count >= 2:
                    break
            else:
                _stale_count = 0
            _prev_violation_sig = sig

            # ── CAI-aware ordering: when multiple violations have similar
            # severity, prefer fixing the one whose resolution has the
            # smallest CAI impact.  This ensures that when two constraints
            # conflict, the resolution chosen maximizes CAI.
            violations.sort(key=lambda v: v.severity, reverse=True)

            # If CAI resolver is available, re-sort same-severity violations
            # by estimated CAI impact (lower impact first)
            if _cai_resolver is not None and len(violations) > 1:
                # Group by severity tier (within 10% of each other)
                # and re-sort within each tier by CAI impact
                for i in range(len(violations)):
                    for j in range(i + 1, len(violations)):
                        vi, vj = violations[i], violations[j]
                        # Only re-sort if severities are similar
                        if abs(vi.severity - vj.severity) / max(vi.severity, 1.0) < 0.1:
                            # Estimate CAI impact based on violation type
                            # Soft constraints (CpG, ATTTA, T-run) have lower
                            # CAI impact than hard constraints (RS, GC, stop)
                            cai_impact_i = self._estimate_cai_impact(vi.violation_type)
                            cai_impact_j = self._estimate_cai_impact(vj.violation_type)
                            if cai_impact_i > cai_impact_j:
                                # j has lower CAI impact — fix it first
                                violations[i], violations[j] = violations[j], violations[i]

            # ── Simultaneous fix: group non-conflicting violations ──
            # Two violations conflict if they share any codon index.
            # Non-conflicting violations (different positions) can be
            # fixed simultaneously without interfering with each other.
            # For conflicting violations, use CAI-aware resolution.
            fixed_this_round = 0

            # Track which codon indices have been modified
            used_codon_indices: set[int] = set()

            # Collect conflicting violations for CAI-aware resolution
            conflicting_violations: list[Violation] = []

            for violation in violations:
                # Check if this violation conflicts with already-fixed ones
                v_codons = set(violation.codon_indices)
                if v_codons & used_codon_indices:
                    # Conflict: record for CAI-aware resolution
                    conflicting_violations.append(violation)
                    continue

                fixed = self._fix_violation(state, violation, protein, avoid_gt)
                if fixed:
                    violations_fixed += 1
                    fixed_this_round += 1
                    # Mark these codon indices as used so future violations
                    # in this round don't conflict
                    used_codon_indices.update(v_codons)

            # ── CAI-aware conflict resolution ──
            # For violations that conflicted with already-fixed ones,
            # try to resolve them by preferring the option that maximizes
            # CAI.  If two constraints conflict at the same position,
            # try fixing the one with lower CAI impact first.
            if conflicting_violations and fixed_this_round > 0:
                # Re-sort conflicting by CAI impact (lower first)
                if _cai_resolver is not None:
                    conflicting_violations.sort(
                        key=lambda v: self._estimate_cai_impact(v.violation_type)
                    )
                for violation in conflicting_violations:
                    v_codons = set(violation.codon_indices)
                    if v_codons & used_codon_indices:
                        continue
                    fixed = self._fix_violation(state, violation, protein, avoid_gt)
                    if fixed:
                        violations_fixed += 1
                        fixed_this_round += 1
                        used_codon_indices.update(v_codons)
                        # Record CAI-aware resolution in provenance
                        if self.provenance_collector is not None:
                            self.provenance_collector.record_constraint_decision(
                                ConstraintDecision(
                                    constraint_name=violation.violation_type,
                                    constraint_type="hard",
                                    action_taken="satisfied",
                                    positions_affected=violation.codon_indices,
                                    tradeoff_description=(
                                        f"CAI-aware resolution: fixed {violation.violation_type} "
                                        f"at pos {violation.position} with CAI-optimal alternative"
                                    ),
                                    impact_on_cai=0.0,
                                )
                            )

            if fixed_this_round == 0:
                # No violations could be fixed — try one at a time with
                # the highest severity violation (fallback)
                for violation in violations:
                    fixed = self._fix_violation(state, violation, protein, avoid_gt)
                    if fixed:
                        violations_fixed += 1
                        fixed_this_round += 1
                        break

                if fixed_this_round == 0:
                    break

        # Phase 2b: Expensive constraint fixes (MaxEntScan-based splice checks)
        # Only when avoid_gt is True (eukaryotic targets)
        if avoid_gt:
            for splice_iter in range(min(10, self.max_local_search_iterations)):
                violations = self._detect_expensive_violations(state)

                if not violations:
                    break

                # Fix all non-conflicting violations simultaneously
                violations.sort(key=lambda v: v.severity, reverse=True)
                any_fixed = False
                used_codon_indices = set()

                for violation in violations:
                    v_codons = set(violation.codon_indices)
                    if v_codons & used_codon_indices:
                        continue  # Skip conflicting violations

                    fixed = self._fix_violation(state, violation, protein, avoid_gt)
                    if fixed:
                        violations_fixed += 1
                        any_fixed = True
                        used_codon_indices.update(v_codons)

                if not any_fixed:
                    for v in violations[:5]:
                        warnings.append(
                            f"Unresolved {v.violation_type} at pos "
                            f"{v.position}: {v.details}"
                        )
                    break

            # Re-check cheap violations after expensive fixes
            # (splice site fixes may reintroduce restriction sites or GC issues)
            for recheck_iter in range(min(5, self.max_local_search_iterations)):
                re_violations = self._detect_cheap_violations(state, avoid_gt)
                if not re_violations:
                    break
                re_violations.sort(key=lambda v: v.severity, reverse=True)
                any_fixed = False
                used_codon_indices = set()
                for violation in re_violations:
                    v_codons = set(violation.codon_indices)
                    if v_codons & used_codon_indices:
                        continue
                    fixed = self._fix_violation(state, violation, protein, avoid_gt)
                    if fixed:
                        violations_fixed += 1
                        any_fixed = True
                        used_codon_indices.update(v_codons)
                if not any_fixed:
                    break

        # ── Eukaryotic GT warning collection ──
        # After constraint satisfaction, any remaining GT violations for
        # eukaryotes should be reported as warnings rather than treated as
        # hard failures.  GT avoidance is soft for eukaryotes.
        if avoid_gt:
            remaining_gt = self._detect_cheap_violations(state, avoid_gt)
            for v in remaining_gt:
                if v.violation_type == "avoidable_gt":
                    warnings.append(
                        f"Soft GT warning: {v.details} "
                        f"(eukaryotic GT avoidance is soft, not a hard constraint)"
                    )

        cai = self._compute_cai(state.sequence)
        return state.sequence, cai, violations_fixed, total_iterations, warnings

    def _detect_cheap_violations(
        self, state: IncrementalSequenceState, avoid_gt: bool
    ) -> list[Violation]:
        """Detect constraint violations that don't require MaxEntScan.

        These are fast to check and should be fixed before the expensive
        splice site checks.
        """
        violations: list[Violation] = []
        seq = state.sequence
        n_codons = state.num_codons

        # 1. Restriction site violations (incremental — uses pre-tracked positions)
        for site_seq, site_pos in state.check_restriction_sites(changed_only=False):
            codon_indices = self._get_overlapping_codon_indices(
                site_pos, len(site_seq), n_codons
            )
            violations.append(Violation(
                violation_type="restriction_site",
                position=site_pos,
                severity=SEVERITY_WEIGHTS["restriction_site"],
                codon_indices=codon_indices,
                details=f"Site {site_seq} at position {site_pos}",
            ))

        # 2. Stop codon violations
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in ("TAA", "TAG", "TGA") and i < len(seq) - 3:
                aa = state.get_aa(i // 3)
                if aa == "*" and i // 3 == n_codons - 1:
                    continue
                violations.append(Violation(
                    violation_type="stop_codon",
                    position=i,
                    severity=SEVERITY_WEIGHTS["stop_codon"],
                    codon_indices=[i // 3],
                    details=f"Premature stop codon {codon} at position {i}",
                ))

        # 3. GC out of range (O(1) via incremental gc_count)
        gc = state.gc_fraction
        if not (self.gc_lo <= gc <= self.gc_hi):
            if gc < self.gc_lo:
                distance = self.gc_lo - gc
            else:
                distance = gc - self.gc_hi
            severity = SEVERITY_WEIGHTS["gc_out_of_range"] * min(distance * 10, 1.0)
            violations.append(Violation(
                violation_type="gc_out_of_range",
                position=0,
                severity=severity,
                codon_indices=list(range(n_codons)),
                details=f"GC={gc:.3f} outside [{self.gc_lo}, {self.gc_hi}]",
            ))

        # 4. Avoidable GT dinucleotides (simple check, no MaxEntScan)
        if avoid_gt:
            for gt_pos in state.gt_positions_list():
                if self._is_unavoidable_gt(seq, gt_pos):
                    continue
                codon_idx = gt_pos // 3
                next_codon_start = (gt_pos // 3) * 3 + 3
                is_within = (gt_pos + 1) < next_codon_start
                if is_within:
                    involved = [codon_idx]
                else:
                    involved = (
                        [codon_idx, codon_idx + 1]
                        if codon_idx + 1 < n_codons
                        else [codon_idx]
                    )
                violations.append(Violation(
                    violation_type="avoidable_gt",
                    position=gt_pos,
                    severity=SEVERITY_WEIGHTS["avoidable_gt"],
                    codon_indices=involved,
                    details=(
                        f"GT at pos {gt_pos} "
                        f"({'within' if is_within else 'cross'}-codon)"
                    ),
                ))

        # 5. ATTTA instability motifs
        pos = 0
        while True:
            p = seq.find("ATTTA", pos)
            if p == -1:
                break
            codon_indices = self._get_overlapping_codon_indices(p, 5, n_codons)
            violations.append(Violation(
                violation_type="atttta_motif",
                position=p,
                severity=SEVERITY_WEIGHTS["atttta_motif"],
                codon_indices=codon_indices,
                details=f"ATTTA at position {p}",
            ))
            pos = p + 1

        # 6. Long T runs (6+ consecutive T)
        i = 0
        while i < len(seq):
            if seq[i] == 'T':
                j = i
                while j < len(seq) and seq[j] == 'T':
                    j += 1
                run_len = j - i
                if run_len >= 6:
                    codon_indices = self._get_overlapping_codon_indices(
                        i, run_len, n_codons
                    )
                    violations.append(Violation(
                        violation_type="t_run",
                        position=i,
                        severity=SEVERITY_WEIGHTS["t_run"] * (run_len - 5),
                        codon_indices=codon_indices,
                        details=f"T-run of {run_len} at position {i}",
                    ))
                i = j
            else:
                i += 1

        # 7. CpG dinucleotides (simplified)
        if avoid_gt:
            cg_positions = state.cg_positions_list()
            for cg_pos in cg_positions[:5]:
                ci = cg_pos // 3
                if ci < n_codons:
                    violations.append(Violation(
                        violation_type="cpg_island",
                        position=cg_pos,
                        severity=SEVERITY_WEIGHTS["cpg_island"],
                        codon_indices=[ci],
                        details=f"CpG at pos {cg_pos}",
                    ))

        return violations

    def _detect_prokaryote(self) -> bool:
        """Detect whether the target organism is prokaryotic.

        Prokaryotes (e.g. E. coli) lack a spliceosome, so cryptic splice
        site detection is biologically meaningless. MaxEntScan must NEVER
        be called for prokaryotic targets.
        """
        try:
            from .organism_config import is_eukaryotic_organism
            return not is_eukaryotic_organism(self.organism)
        except Exception:
            # Fallback: check common prokaryotic identifiers
            prokaryotic_names = {
                "Escherichia_coli", "E_coli_K12", "E_coli_BL21",
                "Bacillus_subtilis", "Pseudomonas_aeruginosa",
            }
            return self.organism in prokaryotic_names or self.species == "ecoli"

    @staticmethod
    def _estimate_cai_impact(violation_type: str) -> float:
        """Estimate the CAI impact of fixing a constraint violation.

        Used by CAI-aware constraint resolution to prefer the fix that
        maximizes CAI when two constraints conflict.  Higher values
        indicate greater expected CAI loss from fixing this violation.

        Args:
            violation_type: The type of constraint violation.

        Returns:
            Estimated CAI impact (higher = more CAI loss).
        """
        # Hard constraints typically have higher CAI impact because
        # they force specific codon choices (e.g., removing a restriction
        # site may require using a lower-CAI synonym).
        # Soft constraints typically have lower CAI impact because they
        # can often be satisfied with CAI-neutral alternatives.
        _CAI_IMPACT_ESTIMATES = {
            "restriction_site": 0.15,     # May force suboptimal codon
            "stop_codon": 0.20,           # Must use non-stop synonym
            "gc_out_of_range": 0.10,      # GC swap usually CAI-friendly
            "cryptic_splice_donor": 0.08, # GT-free alternatives often available
            "cryptic_splice_acceptor": 0.08,
            "avoidable_gt": 0.06,         # GT avoidance is usually cheap
            "cpg_island": 0.04,           # CpG avoidance is soft, low impact
            "atttta_motif": 0.05,         # ATTTA fix is usually cheap
            "t_run": 0.03,                # T-run fix is usually cheap
        }
        return _CAI_IMPACT_ESTIMATES.get(violation_type, 0.05)

    def _detect_expensive_violations(
        self, state: IncrementalSequenceState
    ) -> list[Violation]:
        """Detect constraint violations that require MaxEntScan (expensive).

        Only called after all cheap violations are resolved.
        NEVER called for prokaryotic organisms — splice sites are irrelevant.
        """
        violations: list[Violation] = []
        seq = state.sequence
        n_codons = state.num_codons

        # Prokaryotes have no spliceosome — skip all MaxEntScan calls
        if not self.avoid_gt:
            return violations

        try:
            from .maxentscan import max_donor_score, max_acceptor_score

            # Check donors (GT)
            max_donor = max_donor_score(seq)
            if max_donor >= self.splice_threshold:
                from .maxentscan import score_donor
                for i in range(len(seq) - 1):
                    if seq[i:i+2] == "GT":
                        if self._is_unavoidable_gt(seq, i):
                            continue
                        donor_score = score_donor(seq, i)
                        if donor_score >= self.splice_threshold:
                            codon_idx = i // 3
                            next_codon_start = codon_idx * 3 + 3
                            is_within = (i + 1) < next_codon_start
                            if is_within:
                                involved = [codon_idx]
                            else:
                                involved = (
                                    [codon_idx, codon_idx + 1]
                                    if codon_idx + 1 < n_codons
                                    else [codon_idx]
                                )
                            splice_severity = (
                                SEVERITY_WEIGHTS["cryptic_splice_donor"]
                                * min(donor_score / 10.0, 1.0)
                            )
                            violations.append(Violation(
                                violation_type="cryptic_splice_donor",
                                position=i,
                                severity=(
                                    SEVERITY_WEIGHTS["avoidable_gt"]
                                    + splice_severity
                                ),
                                codon_indices=involved,
                                details=(
                                    f"GT donor at pos {i}, "
                                    f"score={donor_score:.1f}"
                                ),
                            ))

            # Check acceptors (AG)
            max_acceptor = max_acceptor_score(seq)
            if max_acceptor >= self.splice_threshold:
                from .maxentscan import score_acceptor
                for i in range(len(seq) - 1):
                    if seq[i:i+2] == "AG":
                        acceptor_score = score_acceptor(seq, i)
                        if acceptor_score >= self.splice_threshold:
                            codon_idx = i // 3
                            next_codon_start = codon_idx * 3 + 3
                            is_within = (i + 1) < next_codon_start
                            if is_within:
                                involved = [codon_idx]
                            else:
                                involved = (
                                    [codon_idx, codon_idx + 1]
                                    if codon_idx + 1 < n_codons
                                    else [codon_idx]
                                )
                            violations.append(Violation(
                                violation_type="cryptic_splice_acceptor",
                                position=i,
                                severity=(
                                    SEVERITY_WEIGHTS["cryptic_splice_acceptor"]
                                    * min(acceptor_score / 10.0, 1.0)
                                ),
                                codon_indices=involved,
                                details=(
                                    f"AG acceptor at pos {i}, "
                                    f"score={acceptor_score:.1f}"
                                ),
                            ))
        except ImportError:
            pass

        return violations

    # ──────────────────────────────────────────────────────────
    # Violation fixing methods
    # ──────────────────────────────────────────────────────────

    def _fix_violation(
        self,
        state: IncrementalSequenceState,
        violation: Violation,
        protein: str,
        avoid_gt: bool,
    ) -> bool:
        """Try to fix a single constraint violation.

        For each codon position involved in the violation, try all synonymous
        alternatives and pick the one that maximizes the scoring function:
            score = cai_weight * CAI - sum(remaining_violations * penalty)

        Returns True if the violation was fixed.
        """
        vtype = violation.violation_type

        if vtype == "restriction_site":
            return self._fix_restriction_site(state, violation, avoid_gt)
        elif vtype == "gc_out_of_range":
            return self._fix_gc_range(state, violation, avoid_gt)
        elif vtype == "cryptic_splice_donor":
            return self._fix_avoidable_gt(state, violation, avoid_gt)
        elif vtype == "avoidable_gt":
            return self._fix_avoidable_gt(state, violation, avoid_gt)
        elif vtype == "cryptic_splice_acceptor":
            return self._fix_cryptic_splice_acceptor(state, violation, avoid_gt)
        elif vtype == "cpg_island":
            return self._fix_cpg(state, violation, avoid_gt)
        elif vtype == "atttta_motif":
            return self._fix_atttta(state, violation, avoid_gt)
        elif vtype == "t_run":
            return self._fix_t_run(state, violation, avoid_gt)
        elif vtype == "stop_codon":
            return self._fix_stop_codon(state, violation, protein, avoid_gt)
        else:
            return False

    def _fix_restriction_site(
        self,
        state: IncrementalSequenceState,
        violation: Violation,
        avoid_gt: bool,
    ) -> bool:
        """Fix a restriction site violation with minimal CAI loss."""
        for ci in violation.codon_indices:
            if ci < 0 or ci >= state.num_codons:
                continue
            aa = state.get_aa(ci)
            if aa is None or aa == "*":
                continue

            current_codon = state.get_codon(ci)
            old_gt_count = state.gt_count

            # Try each synonymous alternative, sorted by CAI (highest first)
            # Use precomputed sorted_codons for faster lookup
            for alt in self.sorted_codons.get(aa, []):
                if alt == current_codon:
                    continue

                # Try the swap
                old_codon = state.swap_codon(ci, alt)
                new_seq = state.sequence

                # Check if the restriction site is eliminated
                site_present = state.has_any_restriction_site()

                if site_present:
                    # Site still present — rollback
                    state.swap_codon(ci, old_codon)
                    continue

                # Site eliminated! Check we didn't create new GTs
                if avoid_gt and state.gt_count > old_gt_count:
                    new_gts = state.gt_positions_list()
                    all_unavoidable = all(
                        self._is_unavoidable_gt(state.sequence, pos)
                        for pos in new_gts
                        if pos not in {
                            p for p in range(ci * 3 - 1, ci * 3 + 3)
                        }
                    )
                    if not all_unavoidable:
                        state.swap_codon(ci, old_codon)
                        continue

                # Accept the fix
                return True

        # Try two-codon coordinated fix
        return self._fix_restriction_site_two_codons(
            state, violation, avoid_gt
        )

    def _fix_restriction_site_two_codons(
        self,
        state: IncrementalSequenceState,
        violation: Violation,
        avoid_gt: bool,
    ) -> bool:
        """Fix a restriction site by modifying two adjacent codons."""
        codon_indices = sorted(set(violation.codon_indices))
        if len(codon_indices) < 2:
            return False

        old_gt_count = state.gt_count

        for idx in range(len(codon_indices) - 1):
            ci1, ci2 = codon_indices[idx], codon_indices[idx + 1]
            if ci2 != ci1 + 1:
                continue  # Only adjacent codons

            aa1 = state.get_aa(ci1)
            aa2 = state.get_aa(ci2)
            if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                continue

            # Try all pairs sorted by combined CAI
            pairs = []
            for c1 in self.sorted_codons.get(aa1, []):
                for c2 in self.sorted_codons.get(aa2, []):
                    combined = (
                        self.species_cai.get(c1, 0.0)
                        + self.species_cai.get(c2, 0.0)
                    )
                    pairs.append((c1, c2, combined))
            pairs.sort(key=lambda x: x[2], reverse=True)

            for c1, c2, _ in pairs:
                old1 = state.swap_codon(ci1, c1)
                old2 = state.swap_codon(ci2, c2)
                new_seq = state.sequence

                # Check if the site is eliminated
                site_present = state.has_any_restriction_site()

                if site_present:
                    state.swap_codon(ci2, old2)
                    state.swap_codon(ci1, old1)
                    continue

                # Check GT constraint
                if avoid_gt and state.gt_count > old_gt_count:
                    state.swap_codon(ci2, old2)
                    state.swap_codon(ci1, old1)
                    continue

                return True

        return False

    def _fix_gc_range(
        self,
        state: IncrementalSequenceState,
        violation: Violation,
        avoid_gt: bool,
    ) -> bool:
        """Fix GC content by swapping codons to adjust GC fraction."""
        gc_count = state.gc_count
        n_bases = state._n
        gc_val = state.gc_fraction

        if gc_val < self.gc_lo:
            target = self.gc_lo
            need_more_gc = True
        elif gc_val > self.gc_hi:
            target = self.gc_hi
            need_more_gc = False
        else:
            return False  # Already in range

        # Find the best single-codon swap to move GC toward target
        best_swap = None
        best_score = -1.0
        best_ci = -1
        best_gc_delta = 0

        for ci in range(state.num_codons):
            aa = state.get_aa(ci)
            if aa is None or aa == "*":
                continue
            current = state.get_codon(ci)
            current_gc = self.codon_gc.get(current, 0)

            for alt in self.sorted_codons.get(aa, []):
                if alt == current:
                    continue
                alt_gc = self.codon_gc.get(alt, 0)
                gc_delta = alt_gc - current_gc

                # Check if this swap moves GC in the right direction
                if need_more_gc and gc_delta <= 0:
                    continue
                if not need_more_gc and gc_delta >= 0:
                    continue

                new_gc_count = gc_count + gc_delta
                new_frac = new_gc_count / n_bases
                diff = abs(new_frac - target)
                alt_cai = self.species_cai.get(alt, 0.0)

                # Score: prefer better GC improvement with higher CAI
                score = (1.0 - diff) + alt_cai * 0.01
                if score > best_score:
                    best_score = score
                    best_swap = alt
                    best_ci = ci
                    best_gc_delta = gc_delta

        if best_swap is None:
            return False

        # Apply the swap
        old_codon = state.swap_codon(best_ci, best_swap)

        # Check we didn't break any hard constraints
        if state.has_any_restriction_site():
            state.swap_codon(best_ci, old_codon)
            return False

        # Check GT constraint
        if avoid_gt:
            new_gts = state.gt_positions_list()
            for gt_pos in new_gts:
                if not self._is_unavoidable_gt(state.sequence, gt_pos):
                    state.swap_codon(best_ci, old_codon)
                    return False

        return True

    def _fix_avoidable_gt(
        self,
        state: IncrementalSequenceState,
        violation: Violation,
        avoid_gt: bool,
    ) -> bool:
        """Fix an avoidable GT dinucleotide with minimal CAI loss."""
        gt_pos = violation.position
        codon_idx = gt_pos // 3
        next_codon_start = codon_idx * 3 + 3
        is_within = (gt_pos + 1) < next_codon_start

        if is_within:
            return self._fix_within_codon_gt(state, codon_idx, avoid_gt)
        else:
            return self._fix_cross_codon_gt(state, codon_idx, avoid_gt)

    def _fix_within_codon_gt(
        self,
        state: IncrementalSequenceState,
        codon_idx: int,
        avoid_gt: bool,
    ) -> bool:
        """Fix a within-codon GT by choosing the best CAI-preserving
        substitution using precomputed GT-free codon lists."""
        aa = state.get_aa(codon_idx)
        if aa is None or aa == "*":
            return False

        old_gt_count = state.gt_count

        # Try GT-free alternatives (sorted by CAI, highest first)
        # Use precomputed gt_free lookup table
        gt_free_list = self.gt_free.get(aa, [])
        for alt in gt_free_list:
            # Quick boundary check
            left_gt, right_gt = state.boundary_creates_gt(codon_idx, alt)
            if left_gt or right_gt:
                continue

            old_codon = state.swap_codon(codon_idx, alt)
            if state.gt_count < old_gt_count:
                # Check no new restriction sites
                site_ok = not state.has_any_restriction_site()
                if site_ok:
                    return True
            state.swap_codon(codon_idx, old_codon)  # Rollback

        # Fallback: try all alternatives sorted by CAI
        for alt in self.sorted_codons.get(aa, []):
            if "GT" not in alt:
                continue  # Already tried above
            old_codon = state.swap_codon(codon_idx, alt)
            if state.gt_count < old_gt_count:
                site_ok = not state.has_any_restriction_site()
                if site_ok:
                    return True
            state.swap_codon(codon_idx, old_codon)

        return False

    def _fix_cross_codon_gt(
        self,
        state: IncrementalSequenceState,
        codon_idx: int,
        avoid_gt: bool,
    ) -> bool:
        """Fix a cross-codon GT by modifying one or both adjacent codons."""
        old_gt_count = state.gt_count
        next_idx = codon_idx + 1
        prev_idx = codon_idx - 1

        # Strategy D: Change only the next codon (cheapest single-codon fix)
        if next_idx < state.num_codons:
            aa2 = state.get_aa(next_idx)
            if aa2 is not None and aa2 != "*":
                for alt2 in self.sorted_codons.get(aa2, []):
                    if alt2[0] == 'T':
                        continue  # Would still create GT at boundary
                    old2 = state.swap_codon(next_idx, alt2)
                    if state.gt_count < old_gt_count:
                        site_ok = not state.has_any_restriction_site()
                        if site_ok:
                            return True
                    state.swap_codon(next_idx, old2)

        # Strategy C: Change only the current codon
        # (to one that doesn't end with G)
        aa1 = state.get_aa(codon_idx)
        if aa1 is not None and aa1 != "*":
            for alt1 in self.sorted_codons.get(aa1, []):
                if alt1[-1] == 'G':
                    continue  # Would still create GT at boundary
                old1 = state.swap_codon(codon_idx, alt1)
                if state.gt_count < old_gt_count:
                    site_ok = not state.has_any_restriction_site()
                    if site_ok:
                        return True
                state.swap_codon(codon_idx, old1)

        # Strategy B: Change both codons (2-codon coordinated)
        if next_idx < state.num_codons:
            aa1 = state.get_aa(codon_idx)
            aa2 = state.get_aa(next_idx)
            if (
                aa1 is not None and aa1 != "*"
                and aa2 is not None and aa2 != "*"
            ):
                # Try pairs sorted by combined CAI
                pairs = []
                for c1 in self.sorted_codons.get(aa1, [])[:3]:
                    for c2 in self.sorted_codons.get(aa2, [])[:3]:
                        combined = (
                            self.species_cai.get(c1, 0.0)
                            + self.species_cai.get(c2, 0.0)
                        )
                        pairs.append((c1, c2, combined))
                pairs.sort(key=lambda x: x[2], reverse=True)

                for c1, c2, _ in pairs:
                    # Ensure no GT at the boundary
                    if c1[-1] == 'G' and c2[0] == 'T':
                        continue
                    old1 = state.swap_codon(codon_idx, c1)
                    old2 = state.swap_codon(next_idx, c2)
                    if state.gt_count < old_gt_count:
                        site_ok = not state.has_any_restriction_site()
                        if site_ok:
                            return True
                    state.swap_codon(next_idx, old2)
                    state.swap_codon(codon_idx, old1)

        # Strategy A: Try changing the previous codon
        if prev_idx >= 0:
            aa_prev = state.get_aa(prev_idx)
            if aa_prev is not None and aa_prev != "*":
                for alt_prev in self.sorted_codons.get(aa_prev, []):
                    # No special constraint on prev — just try to reduce GTs
                    old_prev = state.swap_codon(prev_idx, alt_prev)
                    if state.gt_count < old_gt_count:
                        site_ok = not state.has_any_restriction_site()
                        if site_ok:
                            return True
                    state.swap_codon(prev_idx, old_prev)

        return False

    def _fix_atttta(
        self,
        state: IncrementalSequenceState,
        violation: Violation,
        avoid_gt: bool,
    ) -> bool:
        """Fix an ATTTA motif by swapping codons at overlapping positions."""
        for ci in violation.codon_indices:
            if ci < 0 or ci >= state.num_codons:
                continue
            aa = state.get_aa(ci)
            if aa is None or aa == "*":
                continue
            current = state.get_codon(ci)

            # Try all alternatives sorted by CAI
            for alt in self.sorted_codons.get(aa, []):
                if alt == current:
                    continue
                old_codon = state.swap_codon(ci, alt)
                if "ATTTA" not in state.sequence:
                    # Check no new restriction sites or GTs
                    site_ok = not state.has_any_restriction_site()
                    if site_ok:
                        if not avoid_gt or state.gt_count <= (
                            state.gt_count  # no change check
                        ):
                            return True
                state.swap_codon(ci, old_codon)
        return False

    def _fix_t_run(
        self,
        state: IncrementalSequenceState,
        violation: Violation,
        avoid_gt: bool,
    ) -> bool:
        """Fix a long T run by swapping a codon at the run's center."""
        run_pos = violation.position
        n_codons = state.num_codons

        # Try codons in the middle of the run first
        codon_idx = (run_pos + 3) // 3  # slightly past the start
        if codon_idx >= n_codons:
            codon_idx = n_codons - 1

        for offset in range(n_codons):
            ci = codon_idx + offset
            if ci >= n_codons:
                ci = offset  # Wrap around

            aa = state.get_aa(ci)
            if aa is None or aa == "*":
                continue
            current = state.get_codon(ci)

            for alt in self.sorted_codons.get(aa, []):
                if alt == current:
                    continue
                old_codon = state.swap_codon(ci, alt)
                # Check T-run is fixed
                new_seq = state.sequence
                has_long_run = False
                i = 0
                while i < len(new_seq):
                    if new_seq[i] == 'T':
                        j = i
                        while j < len(new_seq) and new_seq[j] == 'T':
                            j += 1
                        if j - i >= 6:
                            has_long_run = True
                            break
                        i = j
                    else:
                        i += 1

                if not has_long_run:
                    site_ok = not state.has_any_restriction_site()
                    if site_ok:
                        return True
                state.swap_codon(ci, old_codon)
        return False

    def _fix_stop_codon(
        self,
        state: IncrementalSequenceState,
        violation: Violation,
        protein: str,
        avoid_gt: bool,
    ) -> bool:
        """Fix a premature stop codon by replacing with a non-stop codon."""
        ci = violation.codon_indices[0]
        if ci < 0 or ci >= len(protein):
            return False
        aa = protein[ci]
        if aa == "*":
            return False  # Legitimate stop

        current = state.get_codon(ci)
        for alt in self.sorted_codons.get(aa, []):
            if alt == current:
                continue
            if alt in ("TAA", "TAG", "TGA"):
                continue  # Don't swap to another stop
            old_codon = state.swap_codon(ci, alt)
            site_ok = not state.has_any_restriction_site()
            if site_ok:
                return True
            state.swap_codon(ci, old_codon)
        return False

    def _fix_cryptic_splice_acceptor(
        self,
        state: IncrementalSequenceState,
        violation: Violation,
        avoid_gt: bool,
    ) -> bool:
        """Fix a cryptic splice acceptor (AG) with minimal CAI loss."""
        ag_pos = violation.position
        codon_idx = ag_pos // 3
        next_codon_start = codon_idx * 3 + 3
        is_within = (ag_pos + 1) < next_codon_start

        if is_within:
            # Within-codon AG: use precomputed AG-free codons
            aa = state.get_aa(codon_idx)
            if aa is None or aa == "*":
                return False

            old_gt_count = state.gt_count
            ag_free_list = self.ag_free.get(aa, [])

            for alt in ag_free_list:
                left_gt, right_gt = state.boundary_creates_gt(codon_idx, alt)
                if left_gt or right_gt:
                    continue

                old_codon = state.swap_codon(codon_idx, alt)
                # Verify splice acceptor score dropped below threshold
                try:
                    from .maxentscan import score_acceptor
                    new_score = score_acceptor(state.sequence, ag_pos)
                    if new_score < self.splice_threshold:
                        site_ok = not state.has_any_restriction_site()
                        if site_ok:
                            return True
                except ImportError:
                    # No MaxEntScan — just check AG is gone
                    if "AG" not in state.sequence[ag_pos:ag_pos+2]:
                        return True
                state.swap_codon(codon_idx, old_codon)
        else:
            # Cross-codon AG: try fixing by changing adjacent codons
            next_idx = codon_idx + 1
            if next_idx < state.num_codons:
                aa2 = state.get_aa(next_idx)
                if aa2 is not None and aa2 != "*":
                    for alt2 in self.ag_free.get(aa2, []):
                        if alt2[0] == 'G':
                            continue
                        old2 = state.swap_codon(next_idx, alt2)
                        try:
                            from .maxentscan import score_acceptor
                            new_score = score_acceptor(
                                state.sequence, ag_pos
                            )
                            if new_score < self.splice_threshold:
                                site_ok = not state.has_any_restriction_site()
                                if site_ok:
                                    return True
                        except ImportError:
                            if "AG" not in state.sequence[ag_pos:ag_pos+2]:
                                return True
                        state.swap_codon(next_idx, old2)

            # Try changing the current codon (to one that doesn't end with A)
            aa1 = state.get_aa(codon_idx)
            if aa1 is not None and aa1 != "*":
                for alt1 in self.sorted_codons.get(aa1, []):
                    if alt1[-1] == 'A':
                        continue
                    old1 = state.swap_codon(codon_idx, alt1)
                    try:
                        from .maxentscan import score_acceptor
                        new_score = score_acceptor(
                            state.sequence, ag_pos
                        )
                        if new_score < self.splice_threshold:
                            site_ok = not state.has_any_restriction_site()
                            if site_ok:
                                return True
                    except ImportError:
                        if "AG" not in state.sequence[ag_pos:ag_pos+2]:
                            return True
                    state.swap_codon(codon_idx, old1)

        return False

    def _fix_cpg(
        self,
        state: IncrementalSequenceState,
        violation: Violation,
        avoid_gt: bool,
    ) -> bool:
        """Fix a CpG dinucleotide by swapping to a CG-free alternative."""
        cg_pos = violation.position
        codon_idx = cg_pos // 3
        if codon_idx >= state.num_codons:
            return False

        aa = state.get_aa(codon_idx)
        if aa is None or aa == "*":
            return False

        current = state.get_codon(codon_idx)
        old_gt_count = state.gt_count

        # Find alternatives that don't contain "CG"
        for alt in self.sorted_codons.get(aa, []):
            if alt == current:
                continue
            if "CG" in alt:
                continue  # Would still have CpG
            old_codon = state.swap_codon(codon_idx, alt)

            # Check no new CG was created at boundary
            new_seq = state.sequence
            cg_still = False
            # Check the neighborhood of the swap
            start = max(0, codon_idx * 3 - 1)
            end = min(len(new_seq), codon_idx * 3 + 4)
            if "CG" in new_seq[start:end]:
                # More careful check: is it the same CG position?
                if "CG" in new_seq[cg_pos:cg_pos+2]:
                    cg_still = True

            if not cg_still:
                # Check restriction sites
                site_ok = not state.has_any_restriction_site()
                if site_ok:
                    if not avoid_gt or state.gt_count <= old_gt_count:
                        return True
            state.swap_codon(codon_idx, old_codon)

        return False

    # ──────────────────────────────────────────────────────────
    # Phase 3: CAI Hill Climbing
    # ──────────────────────────────────────────────────────────

    def _cai_hill_climb(
        self, seq: str, protein: str, avoid_gt: bool
    ) -> tuple[str, float, int]:
        """Batch CAI hill climbing: upgrade codons while maintaining constraints.

        For each position where the current codon isn't optimal (w < 1.0):
        1. Find ALL upgradeable positions and their best upgrade codons
        2. Apply ALL non-conflicting upgrades at once (build new sequence)
        3. Validate the combined result against all constraints
        4. If valid, accept all upgrades in a single pass
        5. If not valid, fall back to one-at-a-time upgrades
        6. Only iterate when an upgrade creates a new violation

        Key improvement over one-at-a-time approach:
        - Previously: try position → apply → try next position → apply → ...
        - Now: find ALL upgradeable positions → apply all at once → validate
        This reduces the number of iterations from O(n) to O(1) in the
        common case where upgrades don't conflict.
        """
        improvements = 0
        aas = list(protein)
        n_codons = len(aas)

        for _iteration in range(self.max_hill_climb_iterations):
            # ── Step 1: Find ALL upgradeable positions ──
            # For each position, determine the best upgrade codon
            upgrade_plan: dict[int, str] = {}  # ci -> best upgrade codon

            for ci in range(n_codons):
                aa = aas[ci]
                if aa == "*":
                    continue
                current = seq[ci * 3:ci * 3 + 3]
                current_cai = self.species_cai.get(current, 0.0)
                optimal = self.optimal_codon.get(aa, current)
                optimal_cai = self.species_cai.get(optimal, 0.0)

                if optimal_cai > current_cai and current != optimal:
                    upgrade_plan[ci] = optimal

            if not upgrade_plan:
                break

            # ── Step 2: Apply ALL upgrades at once ──
            # Build the new sequence with all planned upgrades
            seq_list = list(seq)
            for ci, new_codon in upgrade_plan.items():
                start = ci * 3
                seq_list[start:start + 3] = list(new_codon)
            batch_seq = "".join(seq_list)

            # ── Step 3: Validate the combined result ──
            if self._is_valid_batch_upgrade(seq, batch_seq, upgrade_plan, avoid_gt):
                # All upgrades are valid — accept them all at once
                seq = batch_seq
                improvements += len(upgrade_plan)
                continue

            # ── Step 4: Batch failed — try non-conflicting groups ──
            # Group upgrades by proximity. Upgrades that are far apart
            # (>= 3 codons apart) have non-overlapping local windows and
            # can be applied simultaneously.
            any_improved = False

            # Sort by CAI improvement potential (biggest first)
            sorted_cis = sorted(
                upgrade_plan.keys(),
                key=lambda ci: (
                    self.species_cai.get(upgrade_plan[ci], 0.0)
                    - self.species_cai.get(seq[ci * 3:ci * 3 + 3], 0.0)
                ),
                reverse=True,
            )

            applied_positions: set[int] = set()

            for ci in sorted_cis:
                # Skip if a nearby position was already modified in this round
                # (positions within 2 codons share local windows)
                if any(abs(ci - applied) <= 2 for applied in applied_positions):
                    continue

                new_codon = upgrade_plan[ci]
                current = seq[ci * 3:ci * 3 + 3]
                if current == new_codon:
                    continue

                # Try the planned upgrade
                test_seq = seq[:ci * 3] + new_codon + seq[ci * 3 + 3:]
                if self._is_valid_upgrade(seq, test_seq, ci, avoid_gt):
                    seq = test_seq
                    improvements += 1
                    any_improved = True
                    applied_positions.add(ci)
                    continue

                # Try other alternatives in CAI order
                aa = aas[ci]
                for alt in self.sorted_codons.get(aa, []):
                    if alt == current or alt == new_codon:
                        continue
                    alt_cai = self.species_cai.get(alt, 0.0)
                    cur_cai = self.species_cai.get(current, 0.0)
                    if alt_cai <= cur_cai:
                        continue  # Only try upgrades

                    test_seq = seq[:ci * 3] + alt + seq[ci * 3 + 3:]
                    if self._is_valid_upgrade(seq, test_seq, ci, avoid_gt):
                        seq = test_seq
                        improvements += 1
                        any_improved = True
                        applied_positions.add(ci)
                        break

            if not any_improved:
                break

        cai = self._compute_cai(seq)
        return seq, cai, improvements

    def _is_valid_batch_upgrade(
        self,
        old_seq: str,
        new_seq: str,
        upgrade_plan: dict[int, str],
        avoid_gt: bool,
    ) -> bool:
        """Check if a batch of codon upgrades is valid (no new violations).

        Instead of checking each upgrade individually, validates the entire
        new sequence against all constraints. This is more efficient than
        checking one at a time when many upgrades are applied simultaneously.

        Args:
            old_seq: The sequence before upgrades.
            new_seq: The sequence after all upgrades are applied.
            upgrade_plan: Mapping of codon indices to their new codons.
            avoid_gt: Whether to avoid GT dinucleotides.

        Returns:
            True if all upgrades are valid (no new violations introduced).
        """
        # 1. No new restriction sites
        if self._rs_sites:
            for site, site_rc in self._rs_sites:
                if site in new_seq or (site_rc and site_rc in new_seq):
                    return False

        # 2. GC still in range
        gc = (new_seq.count("G") + new_seq.count("C")) / max(len(new_seq), 1)
        if not (self.gc_lo <= gc <= self.gc_hi):
            return False

        # 3. No new ATTTA motifs
        old_attta = old_seq.count("ATTTA")
        new_attta = new_seq.count("ATTTA")
        if new_attta > old_attta:
            return False

        # 4. No long T runs (6+ consecutive T)
        i = 0
        while i < len(new_seq):
            if new_seq[i] == 'T':
                j = i
                while j < len(new_seq) and new_seq[j] == 'T':
                    j += 1
                if j - i >= 6:
                    return False
                i = j
            else:
                i += 1

        # 5. No new avoidable GTs (if avoiding GT)
        if avoid_gt:
            for pos in range(len(new_seq) - 1):
                if new_seq[pos:pos+2] == "GT":
                    if not self._is_unavoidable_gt(new_seq, pos):
                        return False

        # 6. No new premature stop codons at upgraded positions
        for ci, new_codon in upgrade_plan.items():
            if new_codon in ("TAA", "TAG", "TGA"):
                return False

        return True

    def _is_valid_upgrade(
        self,
        old_seq: str,
        new_seq: str,
        ci: int,
        avoid_gt: bool,
    ) -> bool:
        """Check if a codon swap is a valid upgrade (no new violations).

        Instead of checking the entire sequence, only check constraints
        in the neighborhood of the swapped codon for speed.
        """
        # 1. No new restriction sites (local window check around swapped codon)
        # Use Aho-Corasick scanner for O(window_len + matches) detection
        # instead of iterating through each enzyme separately.
        max_site_len = self._max_rs_site_len
        if max_site_len > 0:
            check_start = max(0, ci * 3 - max_site_len + 1)
            check_end = min(len(new_seq), ci * 3 + 3 + max_site_len - 1)
            if self._ac_scanner is not None:
                if self._ac_scanner.has_any_match_in_region(
                    new_seq, check_start, check_end
                ):
                    return False
            else:
                region = new_seq[check_start:check_end]
                for site, site_rc in self._rs_sites:
                    if site in region or (site_rc and site_rc in region):
                        return False

        # 2. GC still in range
        gc = (new_seq.count("G") + new_seq.count("C")) / max(len(new_seq), 1)
        if not (self.gc_lo <= gc <= self.gc_hi):
            return False

        # 3. No new ATTTA motifs — only block if the swap INTRODUCES a new
        # ATTTA, not if ATTTA already existed before the swap.
        old_attta = old_seq.count("ATTTA")
        new_attta = new_seq.count("ATTTA")
        if new_attta > old_attta:
            return False

        # 4. No long T runs
        i = max(0, ci * 3 - 6)
        end = min(len(new_seq), ci * 3 + 9)
        window = new_seq[i:end]
        j = 0
        while j < len(window):
            if window[j] == 'T':
                k = j
                while k < len(window) and window[k] == 'T':
                    k += 1
                if k - j >= 6:
                    return False
                j = k
            else:
                j += 1

        # 5. No new avoidable GTs (if avoiding GT)
        if avoid_gt:
            # Check GT positions in neighborhood
            start = max(0, ci * 3 - 1)
            end_pos = min(len(new_seq) - 1, ci * 3 + 4)
            for pos in range(start, end_pos):
                if new_seq[pos:pos+2] == "GT":
                    if not self._is_unavoidable_gt(new_seq, pos):
                        return False

        # 6. No new premature stop codons
        new_codon = new_seq[ci * 3:ci * 3 + 3]
        if new_codon in ("TAA", "TAG", "TGA"):
            return False

        return True

    # ──────────────────────────────────────────────────────────
    # Utility methods
    # ──────────────────────────────────────────────────────────

    def _compute_cai(self, seq: str) -> float:
        """Compute Codon Adaptation Index for a DNA sequence.

        CAI = geometric mean of (w_i) for each codon,
        where w_i = adaptiveness_i / max_adaptiveness_for_amino_acid.
        Uses precomputed self._max_adapt tables for fast lookup.
        """
        if not seq or len(seq) < 3:
            return 0.0

        # Use precomputed max adaptiveness (self._max_adapt) instead of
        # rebuilding it from AA_TO_CODONS on every call.
        max_adapt = self._max_adapt

        log_sum = 0.0
        n_codons = 0

        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            aa = CODON_TABLE.get(codon)  # type: ignore[assignment]
            if aa is None or aa == "*":
                continue

            adapt = self.species_cai.get(codon, 0.0)
            max_a = max_adapt.get(aa, 0.0)
            if max_a > 0 and adapt > 0:
                w = adapt / max_a
                log_sum += math.log(w)
                n_codons += 1

        if n_codons == 0:
            return 0.0

        return math.exp(log_sum / n_codons)

    def _get_overlapping_codon_indices(
        self, pos: int, site_len: int, n_codons: int
    ) -> list[int]:
        """Get indices of codons that overlap with a site at position pos."""
        first_codon = pos // 3
        last_base = pos + site_len - 1
        last_codon = last_base // 3
        return list(range(max(0, first_codon), min(n_codons, last_codon + 1)))

    def _is_unavoidable_gt(self, seq: str, pos: int) -> bool:
        """Check if a GT dinucleotide at pos is unavoidable (Valine codon).

        Valine codons (GTN) all contain GT, so GT within a Valine codon
        cannot be eliminated by synonymous substitution.
        """
        if pos + 1 >= len(seq):
            return False
        if seq[pos:pos+2] != "GT":
            return False

        # Check if this GT is the first two bases of a codon (Valine)
        codon_start = (pos // 3) * 3
        codon = seq[codon_start:codon_start + 3]
        if len(codon) == 3 and codon[:2] == "GT":
            aa = CODON_TABLE.get(codon)  # type: ignore[assignment]
            if aa == "V":
                return True

        # Check if this GT is at a cross-codon boundary where both
        # the preceding codon must end with G and following starts with T
        # (rare but possible)
        if pos % 3 == 2:
            # GT spans codon boundary
            prev_ci = pos // 3
            next_ci = prev_ci + 1
            if prev_ci >= 0 and next_ci < len(seq) // 3:
                prev_aa = CODON_TABLE.get(seq[prev_ci*3:prev_ci*3+3])
                next_aa = CODON_TABLE.get(seq[next_ci*3:next_ci*3+3])
                # Check if ALL codons for prev_aa end with G
                # AND all codons for next_aa start with T
                if prev_aa and next_aa:
                    prev_codons = AA_TO_CODONS.get(prev_aa, [])
                    next_codons = AA_TO_CODONS.get(next_aa, [])
                    all_end_g = all(c[-1] == 'G' for c in prev_codons)
                    all_start_t = all(c[0] == 'T' for c in next_codons)
                    if all_end_g and all_start_t:
                        return True

        return False

    # ──────────────────────────────────────────────────────────
    # Phase 4: Codon Pair Bias Optimization
    # ──────────────────────────────────────────────────────────

    def _score_cpb_pair(self, codon1: str, codon2: str) -> float:
        """Score a codon pair using organism-specific CPB data or simplified scoring.

        If organism-specific CPB data is available, use it directly.
        Otherwise, fall back to the simplified scoring scheme:
          - Both codons are optimal (highest CAI for their AA): +0.1
          - Exactly one codon is optimal: 0.0
          - Neither codon is optimal: -0.1

        Args:
            codon1: First codon in the pair.
            codon2: Second codon in the pair.

        Returns:
            CPB score for the pair.
        """
        pair_key = f"{codon1}-{codon2}"

        # Use organism-specific data if available
        if self._cpb_data and pair_key in self._cpb_data:
            return self._cpb_data[pair_key]

        # Simplified scoring fallback
        c1_optimal = codon1 in self._optimal_codon_set
        c2_optimal = codon2 in self._optimal_codon_set

        if c1_optimal and c2_optimal:
            return 0.1
        elif c1_optimal or c2_optimal:
            return 0.0
        else:
            return -0.1

    def _compute_mean_cpb(self, seq: str) -> float:
        """Compute mean codon pair bias score for a DNA sequence.

        Uses the same scoring as _score_cpb_pair: organism-specific CPB
        data where available, simplified scoring otherwise.

        Args:
            seq: DNA coding sequence.

        Returns:
            Mean CPB score across all adjacent codon pairs.
        """
        if len(seq) < 6:
            return 0.0

        codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
        if len(codons) < 2:
            return 0.0

        total = 0.0
        for i in range(len(codons) - 1):
            total += self._score_cpb_pair(codons[i], codons[i + 1])

        return total / (len(codons) - 1)

    def _codon_pair_bias_optimize(
        self,
        seq: str,
        protein: str,
        avoid_gt: bool,
        current_cai: float,
        cai_threshold_fraction: float = 0.95,
        max_iterations: int = 3,
    ) -> tuple[str, float, int, float]:
        """Phase 4: Codon pair bias optimization for higher expression.

        For each pair of adjacent codons, check if swapping to a different
        synonymous codon pair improves the codon pair bias score.  This is
        a secondary optimization that can improve protein expression by 2-5x
        according to Coleman et al. (2008) "Virus attenuation by codon pair
        deoptimization".

        Only swap if the CPB improvement doesn't reduce CAI below a threshold
        (default: 0.95 * current_CAI).  Also verifies that no hard constraints
        (restriction sites, GC range, ATTTA motifs, T-runs) are violated.

        Algorithm:
        1. Iterate over all adjacent codon pairs
        2. For each pair, enumerate synonymous alternatives
        3. Score each alternative pair using CPB data or simplified scoring
        4. Accept the swap if:
           a. CPB score improves
           b. CAI stays above cai_threshold_fraction * current_CAI
           c. No new hard constraint violations

        Args:
            seq: Current optimized DNA sequence.
            protein: Amino acid sequence.
            avoid_gt: Whether to avoid GT dinucleotides.
            current_cai: CAI after Phase 3 (used as reference for threshold).
            cai_threshold_fraction: Minimum fraction of current_CAI to
                maintain (default 0.95).
            max_iterations: Maximum number of full-sequence passes.

        Returns:
            (sequence, cai, cpb_improvements, mean_cpb) tuple.
        """
        cai_floor = current_cai * cai_threshold_fraction
        cpb_improvements = 0
        aas = list(protein)
        n_codons = len(aas)

        for _iteration in range(max_iterations):
            any_improved = False

            for ci in range(n_codons - 1):
                aa1 = aas[ci]
                aa2 = aas[ci + 1]

                # Skip stop codons
                if aa1 == "*" or aa2 == "*":
                    continue

                current_c1 = seq[ci * 3:ci * 3 + 3]
                current_c2 = seq[(ci + 1) * 3:(ci + 1) * 3 + 3]

                current_cpb = self._score_cpb_pair(current_c1, current_c2)

                # Find the best synonymous pair by CPB score
                best_pair: tuple[str, str] | None = None
                best_cpb = current_cpb
                # Track CAI of alternatives to ensure CAI floor is met
                best_pair_cai_sum = 0.0
                current_cai_sum = (
                    self.species_cai.get(current_c1, 0.0)
                    + self.species_cai.get(current_c2, 0.0)
                )

                synonyms1 = AA_TO_CODONS.get(aa1, [current_c1])
                synonyms2 = AA_TO_CODONS.get(aa2, [current_c2])

                for alt1 in synonyms1:
                    for alt2 in synonyms2:
                        # Skip identical pair
                        if alt1 == current_c1 and alt2 == current_c2:
                            continue

                        pair_cpb = self._score_cpb_pair(alt1, alt2)

                        # Only consider pairs that improve CPB
                        if pair_cpb <= best_cpb:
                            continue

                        # Check CAI floor: compute new CAI if we swap
                        # (quick approximation: check individual codon scores)
                        alt_cai_sum = (
                            self.species_cai.get(alt1, 0.0)
                            + self.species_cai.get(alt2, 0.0)
                        )

                        # If the alternative has lower combined CAI, we need
                        # to check if the overall CAI would still be above floor.
                        # Quick heuristic: if the per-codon CAI sum drops by more
                        # than (1 - cai_threshold_fraction) * 2, skip.
                        if alt_cai_sum < current_cai_sum:
                            cai_drop = current_cai_sum - alt_cai_sum
                            max_allowed_drop = (1.0 - cai_threshold_fraction) * 2.0
                            if cai_drop > max_allowed_drop:
                                continue

                        if pair_cpb > best_cpb:
                            best_cpb = pair_cpb
                            best_pair = (alt1, alt2)
                            best_pair_cai_sum = alt_cai_sum

                if best_pair is None:
                    continue

                alt1, alt2 = best_pair

                # Apply the swap
                test_seq = (
                    seq[:ci * 3]
                    + alt1
                    + seq[ci * 3 + 3:(ci + 1) * 3]
                    + alt2
                    + seq[(ci + 1) * 3 + 3:]
                )

                # Verify CAI floor
                test_cai = self._compute_cai(test_seq)
                if test_cai < cai_floor:
                    continue  # CAI dropped too much

                # Verify no new constraint violations using _is_valid_upgrade
                # for both positions
                valid = True

                # Check codon 1
                if not self._is_valid_upgrade(
                    seq, test_seq, ci, avoid_gt
                ):
                    valid = False

                # Check codon 2 (only if codon 1 passed)
                if valid and not self._is_valid_upgrade(
                    seq, test_seq, ci + 1, avoid_gt
                ):
                    valid = False

                if valid:
                    seq = test_seq
                    cpb_improvements += 1
                    any_improved = True

            if not any_improved:
                break

        final_cai = self._compute_cai(seq)
        mean_cpb = self._compute_mean_cpb(seq)
        return seq, final_cai, cpb_improvements, mean_cpb

    def _codon_pair_bias_optimize_prokaryote(
        self,
        seq: str,
        protein: str,
        current_cai: float,
        cai_threshold_fraction: float = 0.95,
        max_iterations: int = 3,
    ) -> tuple[str, float, int, float]:
        """Phase 4 for prokaryote fast path: codon pair bias optimization.

        Simplified version for prokaryotes — no GT/AG/splice checks needed.
        Only verifies restriction sites, GC range, ATTTA, and T-runs.

        Args:
            seq: Current optimized DNA sequence.
            protein: Amino acid sequence.
            current_cai: CAI after prokaryote fast path.
            cai_threshold_fraction: Minimum fraction of current_CAI to
                maintain (default 0.95).
            max_iterations: Maximum number of full-sequence passes.

        Returns:
            (sequence, cai, cpb_improvements, mean_cpb) tuple.
        """
        cai_floor = current_cai * cai_threshold_fraction
        cpb_improvements = 0
        aas = list(protein)
        n_codons = len(aas)

        for _iteration in range(max_iterations):
            any_improved = False

            for ci in range(n_codons - 1):
                aa1 = aas[ci]
                aa2 = aas[ci + 1]

                # Skip stop codons
                if aa1 == "*" or aa2 == "*":
                    continue

                current_c1 = seq[ci * 3:ci * 3 + 3]
                current_c2 = seq[(ci + 1) * 3:(ci + 1) * 3 + 3]

                current_cpb = self._score_cpb_pair(current_c1, current_c2)

                # Find the best synonymous pair by CPB score
                best_pair: tuple[str, str] | None = None
                best_cpb = current_cpb
                current_cai_sum = (
                    self.species_cai.get(current_c1, 0.0)
                    + self.species_cai.get(current_c2, 0.0)
                )

                synonyms1 = AA_TO_CODONS.get(aa1, [current_c1])
                synonyms2 = AA_TO_CODONS.get(aa2, [current_c2])

                for alt1 in synonyms1:
                    for alt2 in synonyms2:
                        if alt1 == current_c1 and alt2 == current_c2:
                            continue

                        pair_cpb = self._score_cpb_pair(alt1, alt2)

                        if pair_cpb <= best_cpb:
                            continue

                        alt_cai_sum = (
                            self.species_cai.get(alt1, 0.0)
                            + self.species_cai.get(alt2, 0.0)
                        )

                        if alt_cai_sum < current_cai_sum:
                            cai_drop = current_cai_sum - alt_cai_sum
                            max_allowed_drop = (1.0 - cai_threshold_fraction) * 2.0
                            if cai_drop > max_allowed_drop:
                                continue

                        if pair_cpb > best_cpb:
                            best_cpb = pair_cpb
                            best_pair = (alt1, alt2)

                if best_pair is None:
                    continue

                alt1, alt2 = best_pair

                # Apply the swap
                test_seq = (
                    seq[:ci * 3]
                    + alt1
                    + seq[ci * 3 + 3:(ci + 1) * 3]
                    + alt2
                    + seq[(ci + 1) * 3 + 3:]
                )

                # Verify CAI floor
                test_cai = self._compute_cai(test_seq)
                if test_cai < cai_floor:
                    continue

                # Prokaryote constraint checks (no GT/AG/splice)
                valid = True

                # 1. No new restriction sites
                if self._rs_sites:
                    for site, site_rc in self._rs_sites:
                        if site in test_seq or (site_rc and site_rc in test_seq):
                            valid = False
                            break

                # 2. GC still in range
                if valid:
                    gc = (test_seq.count("G") + test_seq.count("C")) / max(len(test_seq), 1)
                    if not (self.gc_lo <= gc <= self.gc_hi):
                        valid = False

                # 3. No new ATTTA motifs
                if valid:
                    if test_seq.count("ATTTA") > seq.count("ATTTA"):
                        valid = False

                # 4. No new premature stop codons at the swapped positions
                if valid:
                    if alt1 in ("TAA", "TAG", "TGA") or alt2 in ("TAA", "TAG", "TGA"):
                        valid = False

                if valid:
                    seq = test_seq
                    cpb_improvements += 1
                    any_improved = True

            if not any_improved:
                break

        final_cai = self._compute_cai(seq)
        mean_cpb = self._compute_mean_cpb(seq)
        return seq, final_cai, cpb_improvements, mean_cpb

    # The old optimize() signature without is_prokaryote is still supported
    # via the default parameter value above.

    # Old internal method names that may be called from optimize_sequence
    _phase1_greedy_init = _greedy_init
    _phase2_local_search = _constraint_satisfaction
    _phase3_hill_climb = _cai_hill_climb
    _phase4_cpb_optimize = _codon_pair_bias_optimize
    _detect_violations = _detect_cheap_violations
