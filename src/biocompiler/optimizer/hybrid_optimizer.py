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

import math
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, Any

from ..type_system import CODON_TABLE, AA_TO_CODONS, BLOSUM62
from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism
from ..constants import reverse_complement, RESTRICTION_ENZYMES
from ..incremental import IncrementalSequenceState, CodonCache, EnzymeSiteCache
from ..restriction_sites import get_recognition_site
from ..aho_corasick import AhoCorasickScanner, build_scanner_from_enzymes
from ..decision_provenance import CodonDecision, ConstraintDecision

# Import shared types from hybrid_types (avoids circular imports)
from .hybrid_types import (
    HybridResult, Violation, SEVERITY_WEIGHTS,
    GT_CAI_COST_THRESHOLD as _GT_CAI_COST_THRESHOLD,
)

# Import extracted submodules for delegation
from .hybrid_prokaryote import (
    _optimize_prokaryote_fast as _prok_optimize,
    _codon_pair_bias_optimize_prokaryote as _prok_cpb_optimize,
)
from .hybrid_eukaryote import _optimize_eukaryote_fast as _euk_optimize
from .hybrid_constraints import (
    _constraint_satisfaction as _cs_constraint_satisfaction,
    _detect_cheap_violations as _cs_detect_cheap,
    _detect_expensive_violations as _cs_detect_expensive,
    _estimate_cai_impact as _cs_estimate_cai_impact,
    _fix_violation as _cs_fix_violation,
    _fix_restriction_site as _cs_fix_rs,
    _fix_restriction_site_two_codons as _cs_fix_rs_two,
    _fix_gc_range as _cs_fix_gc,
    _fix_avoidable_gt as _cs_fix_gt,
    _fix_within_codon_gt as _cs_fix_within_gt,
    _fix_cross_codon_gt as _cs_fix_cross_gt,
    _fix_atttta as _cs_fix_atttta,
    _fix_t_run as _cs_fix_trun,
    _fix_stop_codon as _cs_fix_stop,
    _fix_cryptic_splice_acceptor as _cs_fix_splice_acc,
    _fix_cpg as _cs_fix_cpg,
)
from .hybrid_hillclimb import (
    _cai_hill_climb as _hc_cai_hill_climb,
    _is_valid_batch_upgrade as _hc_is_valid_batch,
    _is_valid_upgrade as _hc_is_valid_upgrade,
)

# ── NUMBA integration ──────────────────────────────────────────────
try:
    from ..numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        count_gc as _numba_count_gc,
        compute_cai_kernel as _numba_cai_kernel,
        seq_to_bytes as _seq_to_bytes,
    )
except ImportError:
    _HAS_NUMBA = False

# ── Fast MaxEntScan integration ────────────────────────────────────
try:
    from ..maxentscan_fast import scan_splice_sites_fast_str as _scan_splice_sites_fast
    _HAS_FAST_MAXENT = True
except ImportError:
    _HAS_FAST_MAXENT = False

HAS_NUMBA: bool = _HAS_NUMBA

__all__ = [
    "HybridOptimizer",
    "HybridResult",
]

logger = logging.getLogger(__name__)

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
                from ..codon_pair_scoring import get_codon_pair_data
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
            from ..solver.dispatch import solve_with_csp, is_solver_available
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
    # Delegated methods — extracted to optimizer submodules
    # ──────────────────────────────────────────────────────────

    def _optimize_prokaryote_fast(self, protein: str) -> HybridResult:
        """Ultra-fast optimization path for prokaryotic organisms.
        Delegated to optimizer.hybrid_prokaryote module.
        """
        return _prok_optimize(self, protein)

    def _optimize_eukaryote_fast(self, protein: str) -> HybridResult:
        """Streamlined eukaryotic optimization with integrated GT/AG/CpG handling.
        Delegated to optimizer.hybrid_eukaryote module.
        """
        return _euk_optimize(self, protein)

    def _constraint_satisfaction(
        self, seq: str, protein: str, avoid_gt: bool
    ) -> tuple[str, float, int, int, list[str]]:
        """Simultaneous constraint satisfaction with conflict detection.
        Delegated to optimizer.hybrid_constraints module.
        """
        return _cs_constraint_satisfaction(self, seq, protein, avoid_gt)

    def _detect_cheap_violations(
        self, state: IncrementalSequenceState, avoid_gt: bool
    ) -> list[Violation]:
        """Detect cheap constraint violations.
        Delegated to optimizer.hybrid_constraints module.
        """
        return _cs_detect_cheap(self, state, avoid_gt)

    def _detect_expensive_violations(
        self, state: IncrementalSequenceState
    ) -> list[Violation]:
        """Detect expensive constraint violations (MaxEntScan-based).
        Delegated to optimizer.hybrid_constraints module.
        """
        return _cs_detect_expensive(self, state)

    @staticmethod
    def _estimate_cai_impact(violation_type: str) -> float:
        """Estimate the CAI impact of fixing a constraint violation.
        Delegated to optimizer.hybrid_constraints module.
        """
        return _cs_estimate_cai_impact(violation_type)

    def _fix_violation(
        self, state: IncrementalSequenceState, violation: Violation,
        protein: str, avoid_gt: bool,
    ) -> bool:
        """Try to fix a single constraint violation.
        Delegated to optimizer.hybrid_constraints module.
        """
        return _cs_fix_violation(self, state, violation, protein, avoid_gt)

    def _fix_restriction_site(
        self, state: IncrementalSequenceState, violation: Violation, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_rs(self, state, violation, avoid_gt)

    def _fix_restriction_site_two_codons(
        self, state: IncrementalSequenceState, violation: Violation, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_rs_two(self, state, violation, avoid_gt)

    def _fix_gc_range(
        self, state: IncrementalSequenceState, violation: Violation, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_gc(self, state, violation, avoid_gt)

    def _fix_avoidable_gt(
        self, state: IncrementalSequenceState, violation: Violation, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_gt(self, state, violation, avoid_gt)

    def _fix_within_codon_gt(
        self, state: IncrementalSequenceState, codon_idx: int, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_within_gt(self, state, codon_idx, avoid_gt)

    def _fix_cross_codon_gt(
        self, state: IncrementalSequenceState, codon_idx: int, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_cross_gt(self, state, codon_idx, avoid_gt)

    def _fix_atttta(
        self, state: IncrementalSequenceState, violation: Violation, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_atttta(self, state, violation, avoid_gt)

    def _fix_t_run(
        self, state: IncrementalSequenceState, violation: Violation, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_trun(self, state, violation, avoid_gt)

    def _fix_stop_codon(
        self, state: IncrementalSequenceState, violation: Violation,
        protein: str, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_stop(self, state, violation, protein, avoid_gt)

    def _fix_cryptic_splice_acceptor(
        self, state: IncrementalSequenceState, violation: Violation, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_splice_acc(self, state, violation, avoid_gt)

    def _fix_cpg(
        self, state: IncrementalSequenceState, violation: Violation, avoid_gt: bool,
    ) -> bool:
        return _cs_fix_cpg(self, state, violation, avoid_gt)

    def _cai_hill_climb(
        self, seq: str, protein: str, avoid_gt: bool
    ) -> tuple[str, float, int]:
        """Batch CAI hill climbing.
        Delegated to optimizer.hybrid_hillclimb module.
        """
        return _hc_cai_hill_climb(self, seq, protein, avoid_gt)

    def _is_valid_batch_upgrade(
        self, old_seq: str, new_seq: str,
        upgrade_plan: dict[int, str], avoid_gt: bool,
    ) -> bool:
        return _hc_is_valid_batch(self, old_seq, new_seq, upgrade_plan, avoid_gt)

    def _is_valid_upgrade(
        self, old_seq: str, new_seq: str, ci: int, avoid_gt: bool,
    ) -> bool:
        return _hc_is_valid_upgrade(self, old_seq, new_seq, ci, avoid_gt)

    def _codon_pair_bias_optimize_prokaryote(
        self, seq: str, protein: str, current_cai: float,
        cai_threshold_fraction: float = 0.95, max_iterations: int = 3,
    ) -> tuple[str, float, int, float]:
        """Phase 4 for prokaryote fast path: codon pair bias optimization.
        Delegated to optimizer.hybrid_prokaryote module.
        """
        return _prok_cpb_optimize(
            self, seq, protein, current_cai, cai_threshold_fraction, max_iterations
        )
    def _detect_prokaryote(self) -> bool:
        """Detect whether the target organism is prokaryotic.

        Prokaryotes (e.g. E. coli) lack a spliceosome, so cryptic splice
        site detection is biologically meaningless. MaxEntScan must NEVER
        be called for prokaryotic targets.
        """
        try:
            from ..organisms.config import is_eukaryotic_organism
            return not is_eukaryotic_organism(self.organism)
        except Exception:
            # Fallback: check common prokaryotic identifiers
            prokaryotic_names = {
                "Escherichia_coli", "E_coli_K12", "E_coli_BL21",
                "Bacillus_subtilis", "Pseudomonas_aeruginosa",
            }
            return self.organism in prokaryotic_names or self.species == "ecoli"
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
        """Check if a GT dinucleotide at pos is unavoidable or CAI-critical.

        A GT is considered unavoidable if:
        1. It's within a Valine codon (all Val codons start with GT)
        2. It's at a cross-codon boundary where both codons force GT

        A GT is considered CAI-critical if:
        3. It's within a codon where the optimal (highest-CAI) codon
           contains GT, and the best GT-free alternative would cost
           > _GT_CAI_COST_THRESHOLD (3%) in relative adaptiveness.

        CAI-critical GTs are accepted because the CAI cost of avoiding
        them is too high. For eukaryotes, in-codon GTs from optimal
        codons are biologically acceptable, so only trivial CAI
        sacrifices are warranted for GT avoidance.
        """
        if pos + 1 >= len(seq):
            return False
        if seq[pos:pos+2] != "GT":
            return False

        # Check if this GT is within a codon (positions 0-1 or 1-2)
        codon_start = (pos // 3) * 3
        next_codon_start = codon_start + 3
        is_within = (pos + 1) < next_codon_start

        if is_within:
            codon = seq[codon_start:codon_start + 3]
            if len(codon) == 3 and "GT" in codon:
                aa = CODON_TABLE.get(codon)  # type: ignore[assignment]

                # Valine: ALL codons start with GT — truly unavoidable
                if aa == "V":
                    return True

                # CAI-critical check: is this a within-codon GT where the
                # optimal codon contains GT and CAI loss from avoiding it
                # would be > _GT_CAI_COST_THRESHOLD (3%) in relative
                # adaptiveness?
                gt_free_list = self.gt_free.get(aa, [])
                if gt_free_list and aa in self.sorted_codons:
                    optimal = self.sorted_codons[aa][0]
                    if "GT" in optimal:
                        # Optimal codon contains GT — check CAI cost
                        opt_w = self.species_cai.get(optimal, 0.0)
                        best_gtf_w = self.species_cai.get(gt_free_list[0], 0.0)
                        max_a = self._max_adapt.get(aa, 0.0)
                        if max_a > 0:
                            opt_rel = opt_w / max_a
                            best_gtf_rel = best_gtf_w / max_a
                        else:
                            opt_rel = opt_w
                            best_gtf_rel = best_gtf_w
                        if opt_rel - best_gtf_rel > _GT_CAI_COST_THRESHOLD:
                            return True  # GT is CAI-critical

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

