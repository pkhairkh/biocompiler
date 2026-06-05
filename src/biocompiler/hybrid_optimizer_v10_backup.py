"""
BioCompiler HybridOptimizer v11.0.0
=====================================
Greedy + Priority-Based Local Search + CAI Hill Climbing hybrid solver.

Architecture:
  Phase 1: Greedy Initialization     — back-translate with highest-CAI codons
  Phase 2: Constraint Satisfaction    — priority-queue local search with
                                        incremental evaluation
  Phase 3: CAI Hill Climbing          — aggressive CAI recovery while
                                        maintaining all constraints

Key innovation: instead of sequential constraint resolution that can undo
previous fixes, uses a priority queue and incremental constraint evaluation
to avoid conflicts. After each fix, only affected positions are re-evaluated,
preventing the "fix A → break B → fix B → break A" oscillation loop.

Performance target: <1.5ms for GFP (714bp), CAI > 0.98
"""

from __future__ import annotations

import heapq
import math
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, Any

from .type_system import CODON_TABLE, AA_TO_CODONS, BLOSUM62
from .organisms import CODON_ADAPTIVENESS_TABLES, get_species_cai_weights
from .constants import reverse_complement, RESTRICTION_ENZYMES
from .incremental import IncrementalSequenceState, CodonCache, EnzymeSiteCache
from .restriction_sites import get_recognition_site

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
    codon_indices: list      # Codon indices involved
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
    warnings: list = field(default_factory=list)


class HybridOptimizer:
    """Hybrid gene optimizer combining greedy initialization with
    priority-based local search and CAI hill climbing.

    Architecture:
    1. Phase 1: Greedy CAI maximization (best codon per position)
    2. Phase 2: Priority-based constraint satisfaction
       (fix most severe violations first with incremental re-evaluation)
    3. Phase 3: CAI hill climbing
       (upgrade codons while maintaining constraints)

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
    ) -> None:
        self.species = species
        self.organism = organism or self._SPECIES_TO_ORGANISM.get(species, "Homo_sapiens")
        self.species_cai: dict[str, float] = get_species_cai_weights(species)
        self.enzymes = enzymes or []
        self.gc_lo = gc_lo
        self.gc_hi = gc_hi
        self.avoid_gt = avoid_gt
        self.splice_threshold = splice_threshold
        self.cpg_window = cpg_window
        self.cpg_threshold = cpg_threshold
        self.max_local_search_iterations = max_local_search_iterations
        self.max_hill_climb_iterations = max_hill_climb_iterations
        self.cai_weight = cai_weight
        self.provenance_collector = provenance_collector

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
        for enz in self.enzymes:
            site = get_recognition_site(enz)
            if site is not None:
                site_rc = reverse_complement(site)
                self._rs_sites.append((site, site_rc))

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
        # For prokaryotes, override avoid_gt
        effective_avoid_gt = self.avoid_gt and not is_prokaryote

        # Phase 1: Greedy initialization
        seq, phase1_cai = self._greedy_init(protein)

        # Phase 2: Constraint satisfaction with priority queue
        seq, phase2_cai, violations_fixed, iterations, warnings = (
            self._constraint_satisfaction(seq, protein, effective_avoid_gt)
        )

        # Phase 3: CAI hill climbing
        seq, phase3_cai, hill_climb_improvements = (
            self._cai_hill_climb(seq, protein, effective_avoid_gt)
        )

        # Compute final metrics
        gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
        final_cai = self._compute_cai(seq)

        elapsed = _time.monotonic() - _start
        logger.info(
            "HybridOptimizer: protein_len=%d, seq_len=%d, "
            "CAI=%.4f→%.4f→%.4f→%.4f, GC=%.3f, "
            "violations_fixed=%d, hill_climb=%d, iterations=%d, "
            "time=%.1fms",
            len(protein), len(seq),
            phase1_cai, phase2_cai, phase3_cai, final_cai,
            gc, violations_fixed, hill_climb_improvements, iterations,
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
            warnings=warnings,
        )

    # ──────────────────────────────────────────────────────────
    # Phase 1: Greedy Initialization
    # ──────────────────────────────────────────────────────────

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
    # Phase 2: Constraint Satisfaction (Priority-Queue Local Search)
    # ──────────────────────────────────────────────────────────

    def _constraint_satisfaction(
        self, seq: str, protein: str, avoid_gt: bool
    ) -> tuple[str, float, int, int, list[str]]:
        """Priority-based constraint satisfaction with incremental evaluation.

        Algorithm:
        1. Scan for ALL constraint violations at once
        2. Score each violation by severity
           (restriction site > stop codon > ATTTA > T-run > GC > splice > CpG)
        3. Use a PRIORITY QUEUE (heapq) to process violations in severity order
        4. For each violation, find the CAI-optimal fix:
           - Enumerate all codon alternatives at overlapping positions
           - For each alternative, compute:
             score = CAI_weight * cai_delta - constraint_penalty
           - Pick the alternative with the best score
        5. After each fix, only re-evaluate constraints at affected positions
           (incremental update)
        6. Stop when no violations remain or no improvement possible

        Optimizations:
        - Phase 2a: Fix cheap violations first (restriction sites, GC, ATTTA,
          T-runs) using fast checks without MaxEntScan
        - Phase 2b: Fix expensive violations (cryptic splice) using MaxEntScan
          only after all cheap violations are resolved

        Returns:
            (sequence, cai, violations_fixed, iterations, warnings)
        """
        state = IncrementalSequenceState(seq)
        warnings: list[str] = []
        violations_fixed = 0

        # Phase 2a: Fast constraint fixes (no MaxEntScan needed)
        for iteration in range(self.max_local_search_iterations):
            violations = self._detect_cheap_violations(state, avoid_gt)

            if not violations:
                break

            # Build priority queue — highest severity first
            heap: list[Violation] = []
            for v in violations:
                heapq.heappush(heap, v)

            any_fixed = False
            while heap:
                violation = heapq.heappop(heap)
                fixed = self._fix_violation(state, violation, protein, avoid_gt)
                if fixed:
                    violations_fixed += 1
                    any_fixed = True
                    break

            if not any_fixed:
                break

        # Phase 2b: Expensive constraint fixes (MaxEntScan-based splice checks)
        # Only when avoid_gt is True (eukaryotic targets)
        if avoid_gt:
            max_splice_iterations = min(10, self.max_local_search_iterations)
            splice_iter = 0
            for splice_iter in range(max_splice_iterations):
                violations = self._detect_expensive_violations(state)

                if not violations:
                    break

                # Fix all violations found in this scan (batch mode)
                any_fixed = False
                violations.sort(key=lambda v: v.severity, reverse=True)
                for violation in violations:
                    fixed = self._fix_violation(state, violation, protein, avoid_gt)
                    if fixed:
                        violations_fixed += 1
                        any_fixed = True

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
                heap = []
                for v in re_violations:
                    heapq.heappush(heap, v)
                any_fixed = False
                while heap:
                    violation = heapq.heappop(heap)
                    fixed = self._fix_violation(state, violation, protein, avoid_gt)
                    if fixed:
                        violations_fixed += 1
                        any_fixed = True
                        break
                if not any_fixed:
                    break

        cai = self._compute_cai(state.sequence)
        total_iter = iteration + 1 if violations_fixed > 0 else 0
        return state.sequence, cai, violations_fixed, total_iter, warnings

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

        # 1. Restriction site violations
        for site, site_rc in self._rs_sites:
            pos = 0
            while True:
                p = seq.find(site, pos)
                if p == -1:
                    break
                codon_indices = self._get_overlapping_codon_indices(
                    p, len(site), n_codons
                )
                violations.append(Violation(
                    violation_type="restriction_site",
                    position=p,
                    severity=SEVERITY_WEIGHTS["restriction_site"],
                    codon_indices=codon_indices,
                    details=f"Site {site} at position {p}",
                ))
                pos = p + 1

            if site_rc and site_rc != site:
                pos = 0
                while True:
                    p = seq.find(site_rc, pos)
                    if p == -1:
                        break
                    codon_indices = self._get_overlapping_codon_indices(
                        p, len(site_rc), n_codons
                    )
                    violations.append(Violation(
                        violation_type="restriction_site",
                        position=p,
                        severity=SEVERITY_WEIGHTS["restriction_site"],
                        codon_indices=codon_indices,
                        details=f"Site RC {site_rc} at position {p}",
                    ))
                    pos = p + 1

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

        # 3. GC out of range
        gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
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

        # Skip if GT avoidance is disabled (no splice checking needed)
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
                site_present = False
                for site, site_rc in self._rs_sites:
                    if site in new_seq or (site_rc and site_rc in new_seq):
                        site_present = True
                        break

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
                site_present = False
                for site, site_rc in self._rs_sites:
                    if site in new_seq or (site_rc and site_rc in new_seq):
                        site_present = True
                        break

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
        seq = state.sequence
        gc_count = sum(1 for b in seq if b in "GC")
        n_bases = len(seq)
        gc_val = gc_count / n_bases

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
        new_seq = state.sequence
        for site, site_rc in self._rs_sites:
            if site in new_seq or (site_rc and site_rc in new_seq):
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
                site_ok = True
                for site, site_rc in self._rs_sites:
                    if site in state.sequence or (
                        site_rc and site_rc in state.sequence
                    ):
                        site_ok = False
                        break
                if site_ok:
                    return True
            state.swap_codon(codon_idx, old_codon)  # Rollback

        # Fallback: try all alternatives sorted by CAI
        for alt in self.sorted_codons.get(aa, []):
            if "GT" not in alt:
                continue  # Already tried above
            old_codon = state.swap_codon(codon_idx, alt)
            if state.gt_count < old_gt_count:
                site_ok = True
                for site, site_rc in self._rs_sites:
                    if site in state.sequence or (
                        site_rc and site_rc in state.sequence
                    ):
                        site_ok = False
                        break
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
                        site_ok = True
                        for site, site_rc in self._rs_sites:
                            if site in state.sequence or (
                                site_rc and site_rc in state.sequence
                            ):
                                site_ok = False
                                break
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
                    site_ok = True
                    for site, site_rc in self._rs_sites:
                        if site in state.sequence or (
                            site_rc and site_rc in state.sequence
                        ):
                            site_ok = False
                            break
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
                        site_ok = True
                        for site, site_rc in self._rs_sites:
                            if site in state.sequence or (
                                site_rc and site_rc in state.sequence
                            ):
                                site_ok = False
                                break
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
                        site_ok = True
                        for site, site_rc in self._rs_sites:
                            if site in state.sequence or (
                                site_rc and site_rc in state.sequence
                            ):
                                site_ok = False
                                break
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
                    site_ok = True
                    for site, site_rc in self._rs_sites:
                        if site in state.sequence or (
                            site_rc and site_rc in state.sequence
                        ):
                            site_ok = False
                            break
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
                    site_ok = True
                    for site, site_rc in self._rs_sites:
                        if site in new_seq or (
                            site_rc and site_rc in new_seq
                        ):
                            site_ok = False
                            break
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
            site_ok = True
            for site, site_rc in self._rs_sites:
                if site in state.sequence or (
                    site_rc and site_rc in state.sequence
                ):
                    site_ok = False
                    break
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
                        site_ok = True
                        for site, site_rc in self._rs_sites:
                            if site in state.sequence or (
                                site_rc and site_rc in state.sequence
                            ):
                                site_ok = False
                                break
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
                                site_ok = True
                                for site, site_rc in self._rs_sites:
                                    if site in state.sequence or (
                                        site_rc and site_rc in state.sequence
                                    ):
                                        site_ok = False
                                        break
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
                            site_ok = True
                            for site, site_rc in self._rs_sites:
                                if site in state.sequence or (
                                    site_rc and site_rc in state.sequence
                                ):
                                    site_ok = False
                                    break
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
                site_ok = True
                for site, site_rc in self._rs_sites:
                    if site in new_seq or (
                        site_rc and site_rc in new_seq
                    ):
                        site_ok = False
                        break
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
        """CAI hill climbing: upgrade codons while maintaining constraints.

        For each position where the current codon isn't optimal (w < 1.0):
        1. Try swapping to the optimal codon
        2. If no new violations introduced, accept
        3. If violations introduced, try next-best codon
        4. Repeat until convergence

        Uses a greedy approach: try positions in order of potential CAI
        improvement (biggest gap between current and optimal first).
        """
        improvements = 0
        aas = list(protein)
        n_codons = len(aas)

        for _iteration in range(self.max_hill_climb_iterations):
            any_improved = False

            # Build a priority list: positions sorted by CAI improvement potential
            improvement_candidates = []
            for ci in range(n_codons):
                aa = aas[ci]
                if aa == "*":
                    continue
                current = seq[ci * 3:ci * 3 + 3]
                current_cai = self.species_cai.get(current, 0.0)
                optimal = self.optimal_codon.get(aa, current)
                optimal_cai = self.species_cai.get(optimal, 0.0)

                if optimal_cai > current_cai:
                    improvement_candidates.append(
                        (optimal_cai - current_cai, ci, optimal)
                    )

            # Sort by improvement potential (biggest first)
            improvement_candidates.sort(key=lambda x: x[0], reverse=True)

            for _, ci, optimal in improvement_candidates:
                aa = aas[ci]
                current = seq[ci * 3:ci * 3 + 3]
                if current == optimal:
                    continue  # Already optimal

                # Try the optimal codon first
                test_seq = seq[:ci * 3] + optimal + seq[ci * 3 + 3:]
                if self._is_valid_upgrade(
                    seq, test_seq, ci, avoid_gt
                ):
                    seq = test_seq
                    improvements += 1
                    any_improved = True
                    continue

                # Try other alternatives in CAI order
                for alt in self.sorted_codons.get(aa, []):
                    if alt == current or alt == optimal:
                        continue
                    alt_cai = self.species_cai.get(alt, 0.0)
                    cur_cai = self.species_cai.get(current, 0.0)
                    if alt_cai <= cur_cai:
                        continue  # Only try upgrades

                    test_seq = seq[:ci * 3] + alt + seq[ci * 3 + 3:]
                    if self._is_valid_upgrade(
                        seq, test_seq, ci, avoid_gt
                    ):
                        seq = test_seq
                        improvements += 1
                        any_improved = True
                        break

            if not any_improved:
                break

        cai = self._compute_cai(seq)
        return seq, cai, improvements

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
        # 1. No new restriction sites
        for site, site_rc in self._rs_sites:
            if site in new_seq or (site_rc and site_rc in new_seq):
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
        """
        if not seq or len(seq) < 3:
            return 0.0

        # Precompute max adaptiveness per amino acid
        max_adapt: dict[str, float] = {}
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*":
                continue
            max_w = 0.0
            for c in codons:
                w = self.species_cai.get(c, 0.0)
                if w > max_w:
                    max_w = w
            max_adapt[aa] = max_w

        log_sum = 0.0
        n_codons = 0

        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            aa = CODON_TABLE.get(codon)
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
            aa = CODON_TABLE.get(codon)
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
    # Backward compatibility: old method names as aliases
    # ──────────────────────────────────────────────────────────

    # The old optimize() signature without is_prokaryote is still supported
    # via the default parameter value above.

    # Old internal method names that may be called from optimize_sequence
    _phase1_greedy_init = _greedy_init
    _phase2_local_search = _constraint_satisfaction
    _phase3_hill_climb = _cai_hill_climb
    _detect_violations = _detect_cheap_violations
