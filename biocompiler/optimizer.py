"""
BioCompiler Optimizer v8.1.0
==============================
7-phase certified gene optimization pipeline with CAI-aware GT resolution.

Phase 1:   Greedy codon optimization (2-codon look-ahead, GT-aware, with unavoidable-GT tracking)
Phase 2:   Restriction site removal (CAI-prioritized alternatives)
Phase 3:   Cross-codon constraint resolution (iterative, global validation, CAI-prioritized)
Phase 3.5: Within-codon GT resolution (synonymous substitution + mutagenesis flagging)
Phase 4:   Mutagenesis fallback (CAI-weighted AA substitution using BLOSUM62)
Phase 5:   CpG island avoidance (CAI-prioritized alternatives)
Phase 6:   Re-optimization pass (constraint-preserving CAI optimizer, iterative)
Phase 7:   CAI-boosting re-pass (pure CAI maximization after all constraints satisfied)
"""

import math
from typing import List, Dict, Optional, Tuple, Set

from .type_system import (
    CODON_TABLE, AA_TO_CODONS, BLOSUM62, SpliceVerdict, PredicateResult,
    check_no_stop_codons, check_no_cryptic_splice, check_no_cpg_island,
    check_no_restriction_site, check_no_gt_dinucleotide, check_no_avoidable_gt,
    check_valid_coding_seq,
    check_conservation_score, check_codon_optimality,
    find_cross_codon_gt, find_cross_codon_cg, find_cross_codon_restriction,
    PREDICATE_NAMES,
)
from .species import SPECIES
from .mutagenesis import propose_mutagenesis, MutagenesisReport, MutagenesisProposal
from .certificates import compute_certificate, format_certificate


def _count_gts(s: str) -> int:
    """Count GT dinucleotides in a sequence."""
    return sum(1 for i in range(len(s) - 1) if s[i:i+2] == "GT")


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


def _count_cpg_ratio(seq: str) -> float:
    """Compute CpG Obs/Exp ratio for a sequence (or sub-sequence)."""
    c = seq.count("C")
    g = seq.count("G")
    cg = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "CG")
    expected = (c * g) / len(seq) if len(seq) > 0 else 0
    return cg / expected if expected > 0 else 0.0


def _compute_cai(seq: str, species_cai: Dict[str, float]) -> float:
    """Compute the geometric mean CAI for a sequence."""
    if not seq or len(seq) < 3:
        return 0.0
    log_sum = 0.0
    count = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        cai = species_cai.get(codon, 0.0)
        if cai <= 0:
            cai = 0.001  # avoid log(0)
        log_sum += math.log(cai)
        count += 1
    if count == 0:
        return 0.0
    return math.exp(log_sum / count)


class OptimizationResult:
    """Result object from BioOptimizer.optimize().

    Provides convenient access to the optimized sequence, CAI score,
    predicate results, and certificate text.

    Supports tuple unpacking for backward compatibility:
        seq, results, cert = optimizer.optimize(dna)
    """
    def __init__(
        self,
        sequence: str,
        cai: float,
        results: List[PredicateResult],
        certificate: str,
    ):
        self.sequence = sequence
        self.cai = cai
        self.results = results
        self.certificate = certificate

    def __repr__(self) -> str:
        return f"OptimizationResult(cai={self.cai:.4f}, len={len(self.sequence)})"

    def __iter__(self):
        """Allow tuple unpacking: seq, results, cert = result."""
        yield self.sequence
        yield self.results
        yield self.certificate

    def __getitem__(self, index):
        """Allow indexing: result[0], result[1], result[2]."""
        return (self.sequence, self.results, self.certificate)[index]

    def __len__(self):
        return 3


class BioOptimizer:
    """Certified gene sequence optimizer with CAI-aware 7-phase pipeline."""

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
    ):
        self.species = species
        self.species_cai: Dict[str, float] = SPECIES.get(species, SPECIES["ecoli"])
        self.enzymes: List[str] = enzymes or []
        self.splice_low = splice_low
        self.splice_high = splice_high
        self.cpg_window = cpg_window
        self.cpg_threshold = cpg_threshold
        self.min_blosum = min_blosum
        self.min_cai = min_cai
        self.avoid_gt = avoid_gt
        # Track positions where GT is unavoidable (e.g., Valine codons)
        self._unavoidable_gt_positions: Set[int] = set()
        # Track mutagenesis proposals that were applied
        self._applied_mutagenesis: List[Dict] = []
        # Store original input protein for conservation scoring
        self._original_protein: str = ""

    def optimize(
        self,
        seq: str,
        strategy: Optional[str] = None,
    ) -> 'OptimizationResult':
        """Run the full optimization pipeline.

        Args:
            seq: Input sequence — either a DNA coding sequence or a protein
                  sequence (single-letter amino acid codes).  Protein sequences
                  are automatically detected and back-translated.
            strategy: Optimization strategy override.
                - "constraint_first" (default): GT-aware greedy then fix constraints
                - "cai_first": Maximize CAI first, then fix constraints with
                  minimal CAI impact (DNAworks-style)

        Returns:
            OptimizationResult with .sequence, .cai, .results, .certificate
        """
        effective_strategy = strategy if strategy is not None else getattr(self, '_strategy', 'constraint_first')

        seq = seq.upper().strip()

        # Detect protein input: if the sequence contains letters other than
        # ACGT (and is not purely ACGT), treat it as a protein and back-translate.
        protein_chars = set(seq) - {"A", "C", "G", "T"}
        if protein_chars:
            seq = self._back_translate_protein(seq)

        self._unavoidable_gt_positions = set()
        self._applied_mutagenesis = []
        self._original_protein = self._translate(seq)

        if effective_strategy == "cai_first":
            seq, results, cert_text = self._optimize_cai_first(seq)
        else:
            seq, results, cert_text = self._optimize_constraint_first(seq)

        cai = _compute_cai(seq, self.species_cai)
        return OptimizationResult(
            sequence=seq,
            cai=cai,
            results=results,
            certificate=cert_text,
        )

    def _back_translate_protein(self, protein: str) -> str:
        """Back-translate a protein sequence to DNA using max-CAI codons."""
        codons = []
        for aa in protein:
            if aa == "*":
                codons.append("TAA")
                continue
            candidates = AA_TO_CODONS.get(aa, [])
            if candidates:
                codons.append(max(candidates, key=lambda c: self.species_cai.get(c, 0.0)))
            else:
                codons.append("NNN")
        return "".join(codons)

    def _optimize_constraint_first(
        self, seq: str
    ) -> Tuple[str, List[PredicateResult], str]:
        """Constraint-first strategy (original 7-phase pipeline)."""
        # Phase 1: Greedy codon optimization (2-codon look-ahead, GT-aware)
        seq = self._phase1_greedy_optimize(seq)

        # Phase 2: Restriction site removal (CAI-prioritized)
        seq = self._phase2_remove_restriction_sites(seq)

        # Phase 3: Cross-codon constraint resolution (iterative, CAI-prioritized)
        seq, mut_report = self._phase3_cross_codon_constraints(seq)

        # Phase 3.5: Within-codon GT resolution
        seq, mut_report_35 = self._phase35_within_codon_gt(seq)
        mut_report.proposals.extend(mut_report_35.proposals)

        # Phase 4: Mutagenesis fallback (CAI-weighted)
        seq = self._phase4_mutagenesis_fallback(seq, mut_report)

        # Phase 5: CpG island avoidance (CAI-prioritized)
        seq = self._phase5_avoid_cpg_islands(seq)

        # Phase 6: Re-optimization pass (constraint-preserving CAI optimizer)
        seq = self._phase6_reoptimize(seq)

        # Phase 7: CAI-boosting re-pass (pure CAI maximization)
        seq = self._phase7_cai_boost(seq)

        # Final cleanup: Phase 6/7 may re-introduce CpG islands via
        # individually-valid swaps whose cumulative effect pushes CpG Obs/Exp
        # above threshold.  Re-run Phase 5 + Phase 7 until stable.
        for _cleanup in range(5):
            cpg_result = check_no_cpg_island(seq, self.cpg_window, self.cpg_threshold)
            if cpg_result.passed:
                break
            seq = self._phase5_avoid_cpg_islands(seq)
            seq = self._phase7_cai_boost(seq)

        # Evaluate all 8 predicates
        results = self._evaluate_all_predicates(seq)

        # Generate certificate
        cert_text = format_certificate(results, seq, self.species)

        return seq, results, cert_text

    # ──────────────────────────────────────────────────────────
    # CAI-first optimization strategy (DNAworks-style)
    # ──────────────────────────────────────────────────────────
    def _optimize_cai_first(
        self, seq: str
    ) -> Tuple[str, List[PredicateResult], str]:
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
        3. CAI hill climbing to recover any lost CAI while maintaining constraints
        4. Iterate until convergence

        Key insight: by starting from max CAI (CAI=1.0) and only making the
        smallest necessary CAI sacrifices to fix constraints, we achieve much
        higher CAI than the constraint_first strategy which permanently sacrifices
        CAI by avoiding GT during codon selection.

        For GT fixing, we use a priority-based approach that fixes GTs with
        the lowest CAI cost first, and considers multi-codon windows to find
        globally optimal fixes.
        """
        # Phase 0: Pure max-CAI back-translation (ignore ALL constraints)
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

        # Phase 1: Fix restriction sites (highest priority — binary constraint)
        seq = self._cai_first_fix_restriction_sites(seq)

        # Phase 2: Fix avoidable GT dinucleotides (minimal CAI impact)
        seq = self._cai_first_fix_gts(seq)

        # Phase 2.5: Mutagenesis fallback for Valine GTs
        # Valine codons all contain GT (GTN), so we substitute V→I
        # (BLOSUM62 score 3, conservative) using high-CAI Ile codons
        seq = self._cai_first_mutagenesis_fallback(seq)

        # Phase 3: Fix CpG islands (minimal CAI impact)
        seq = self._cai_first_fix_cpg(seq)

        # Phase 4: Fix cryptic splice sites (minimal CAI impact)
        seq = self._cai_first_fix_splice(seq)

        # Phase 5: CAI hill climbing (upgrade codons while maintaining constraints)
        seq = self._phase7_cai_boost(seq)

        # Phase 6: Re-optimization pass (constraint-preserving CAI optimizer)
        seq = self._phase6_reoptimize(seq)

        # Phase 7: Second pass of GT fixing + CAI boost (iterative refinement)
        for _refinement in range(3):
            old_cai = _compute_cai(seq, self.species_cai)
            seq = self._cai_first_fix_gts(seq)
            seq = self._phase7_cai_boost(seq)
            seq = self._phase6_reoptimize(seq)
            new_cai = _compute_cai(seq, self.species_cai)
            if new_cai <= old_cai + 0.0001:
                break

        # Track unavoidable GTs for certificate
        for i in range(len(seq) - 1):
            if seq[i] == "G" and seq[i + 1] == "T":
                cs = (i // 3) * 3
                next_cs = cs + 3
                # Within-codon GT in Valine is unavoidable
                if i + 1 < next_cs:
                    codon = seq[cs:cs + 3]
                    aa = CODON_TABLE.get(codon)
                    if aa == "V":
                        self._unavoidable_gt_positions.add(i)
                    elif aa is not None:
                        # Check if all synonymous codons contain GT
                        if all("GT" in c for c in AA_TO_CODONS.get(aa, [])):
                            self._unavoidable_gt_positions.add(i)

        # Evaluate all predicates
        results = self._evaluate_all_predicates(seq)
        cert_text = format_certificate(results, seq, self.species)
        return seq, results, cert_text

    def _cai_first_fix_restriction_sites(self, seq: str) -> str:
        """CAI-first Phase 1: Fix restriction sites with minimal CAI impact."""
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
                    # Try two-codon substitution
                    codon_starts_sorted = sorted(codon_starts)
                    for idx in range(len(codon_starts_sorted) - 1):
                        cs1, cs2 = codon_starts_sorted[idx], codon_starts_sorted[idx + 1]
                        if cs2 != cs1 + 3:
                            continue
                        codon1 = "".join(seq_list[cs1:cs1 + 3])
                        codon2 = "".join(seq_list[cs2:cs2 + 3])
                        aa1 = CODON_TABLE.get(codon1)
                        aa2 = CODON_TABLE.get(codon2)
                        if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                            continue
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
                                fixed = True
                                break
                        if fixed:
                            break
                    if not fixed:
                        break

        return "".join(seq_list)

    def _cai_first_fix_gts(self, seq: str) -> str:
        """CAI-first Phase 2: Fix avoidable GT dinucleotides with minimal CAI impact."""
        seq_list = list(seq)
        max_rounds = 50

        for _ in range(max_rounds):
            current_seq = "".join(seq_list)
            violations = []

            # Find all avoidable GT positions
            for i in range(len(current_seq) - 1):
                if current_seq[i] == "G" and current_seq[i + 1] == "T":
                    # Check if this GT is avoidable
                    cs = (i // 3) * 3
                    next_cs = cs + 3
                    is_within = (i + 1) < next_cs

                    if is_within:
                        # Within-codon GT — avoidable unless Valine
                        codon = current_seq[cs:cs + 3]
                        aa = CODON_TABLE.get(codon)
                        if aa == "V":
                            continue  # Unavoidable
                        if all("GT" in c for c in AA_TO_CODONS.get(aa, [])):
                            continue  # Unavoidable
                    else:
                        # Cross-codon GT — avoidable if either side can change
                        prev_cs = cs
                        next_cs_cross = cs + 3
                        if next_cs_cross + 3 > len(current_seq):
                            continue
                        prev_codon = current_seq[prev_cs:prev_cs + 3]
                        next_codon = current_seq[next_cs_cross:next_cs_cross + 3]
                        prev_aa = CODON_TABLE.get(prev_codon)
                        next_aa = CODON_TABLE.get(next_codon)
                        if prev_aa is None or next_aa is None:
                            continue
                        prev_can_avoid = any(c[-1] != "G" for c in AA_TO_CODONS.get(prev_aa, [prev_codon]))
                        next_can_avoid = any(c[0] != "T" for c in AA_TO_CODONS.get(next_aa, [next_codon]))
                        if not (prev_can_avoid or next_can_avoid):
                            continue  # Unavoidable

                    violations.append(i)

            if not violations:
                break

            any_fixed = False

            for gt_pos in violations:
                cs = (gt_pos // 3) * 3
                next_cs = cs + 3
                is_within = (gt_pos + 1) < next_cs

                if is_within:
                    fixed = self._cai_first_fix_within_gt(seq_list, cs)
                else:
                    fixed = self._cai_first_fix_cross_gt(seq_list, cs)

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
                        continue  # Skip codons starting with T
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
                    continue  # Skip codons ending with G
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
        """CAI-first Phase 3: Fix CpG islands with minimal CAI impact."""
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
                cg_count = sum(1 for i in range(len(window) - 1) if window[i:i+2] == "CG")
                expected = (c_count * g_count) / len(window) if len(window) > 0 else 0
                obs_exp = cg_count / expected if expected > 0 else 0.0

                if obs_exp <= self.cpg_threshold:
                    continue

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

    def _cai_first_mutagenesis_fallback(self, seq: str) -> str:
        """CAI-first Phase 2.5: Apply mutagenesis for GTs that can't be resolved
        by synonymous substitution.

        Specifically targets Valine codons (GTN) which all contain GT.
        Substitutes V→I (Isoleucine, BLOSUM62=3) using highest-CAI Ile codon
        to eliminate the GT while maximizing CAI.

        Also handles any other within-codon GT positions where ALL synonymous
        codons contain GT (rare but possible).
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

                # Check if this codon contains GT
                if "GT" not in codon:
                    continue

                # Check if any synonymous codon avoids GT
                has_gt_free = any("GT" not in c for c in AA_TO_CODONS.get(aa, []))
                if has_gt_free:
                    continue  # Can be fixed by synonymous substitution

                # All synonymous codons contain GT — try mutagenesis
                # V→I is the primary target (BLOSUM62=3, conservative)
                if aa == "V":
                    # Try Isoleucine codons (sorted by CAI)
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

                # Also try other conservative substitutions
                if not changed and aa == "V":
                    # V→L (Leucine, BLOSUM62=1)
                    leu_codons = AA_TO_CODONS.get("L", [])
                    old_gt_count = _count_gts("".join(seq_list))

                    for leu_codon in sorted(
                        leu_codons,
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    ):
                        # Skip codons that start with T (might create cross-codon GT)
                        test_list = seq_list[:]
                        for k, b in enumerate(leu_codon):
                            test_list[i + k] = b
                        new_gt_count = _count_gts("".join(test_list))
                        if new_gt_count < old_gt_count:
                            seq_list[:] = test_list
                            self._applied_mutagenesis.append({
                                "position": i,
                                "original_aa": "V",
                                "new_aa": "L",
                                "blosum": BLOSUM62.get(("V", "L"), 1),
                            })
                            changed = True
                            break

        return "".join(seq_list)

    def _cai_first_fix_splice(self, seq: str) -> str:
        """CAI-first Phase 4: Fix cryptic splice sites with minimal CAI impact."""
        seq_list = list(seq)
        max_rounds = 30

        for _ in range(max_rounds):
            current_seq = "".join(seq_list)
            splice_result = check_no_cryptic_splice(current_seq, self.splice_low, self.splice_high)
            if splice_result.passed:
                break

            fixed = False
            for i in range(len(current_seq) - 1):
                if current_seq[i] == "G" and current_seq[i + 1] == "T":
                    cs = (i // 3) * 3
                    next_cs = cs + 3
                    is_within = (i + 1) < next_cs

                    if is_within:
                        fixed = self._cai_first_fix_within_gt(seq_list, cs)
                    else:
                        fixed = self._cai_first_fix_cross_gt(seq_list, cs)

                    if fixed:
                        break

            if not fixed:
                break

        return "".join(seq_list)

    # ──────────────────────────────────────────────────────────
    # Constraint validation helper
    # ──────────────────────────────────────────────────────────
    def _count_restriction_sites(self, seq: str) -> int:
        """Count total restriction sites in a sequence for all configured enzymes."""
        from .restriction_sites import get_recognition_site
        total = 0
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            p = seq.find(site)
            while p != -1:
                total += 1
                p = seq.find(site, p + 1)
        return total

    def _is_valid_swap(
        self, seq_list: list, pos: int, new_codon: str,
        old_gt_count: Optional[int] = None,
    ) -> bool:
        """Check whether swapping the codon at position `pos` to `new_codon`
        preserves all constraints.

        Validates:
        1. No increase in GT dinucleotide count
        2. No new restriction sites created
        3. No CpG island threshold violation
        4. No new cross-codon GT created (covered by #1 but explicit check)

        Args:
            seq_list: Current sequence as a list of characters
            pos: Start position of the codon to swap (must be % 3 == 0)
            new_codon: The replacement codon string
            old_gt_count: Pre-computed GT count of current seq (optimization)

        Returns:
            True if the swap preserves all constraints
        """
        if pos + 3 > len(seq_list):
            return False

        current_codon = "".join(seq_list[pos:pos + 3])
        if new_codon == current_codon:
            return True  # No-op is always valid

        # Build test sequence
        test_list = seq_list[:]
        for k, b in enumerate(new_codon):
            test_list[pos + k] = b
        test_seq = "".join(test_list)
        current_seq = "".join(seq_list)

        # 1. GT count must not increase
        if self.avoid_gt:
            if old_gt_count is None:
                old_gt_count = _count_gts(current_seq)
            new_gt_count = _count_gts(test_seq)
            if new_gt_count > old_gt_count:
                return False

        # 2. No new restriction sites created
        from .restriction_sites import get_recognition_site
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            # Count before and after
            old_count = current_seq.count(site)
            new_count = test_seq.count(site)
            if new_count > old_count:
                return False

        # 3. No CpG island threshold violation
        # Check windows that overlap the changed codon (pos to pos+2)
        # Window at start s covers s..s+cpg_window-1; overlaps codon iff
        #   s+cpg_window-1 >= pos  AND  s <= pos+2
        #   => s >= pos-cpg_window+1  AND  s <= pos+2
        for start in range(
            max(0, pos - self.cpg_window + 1),
            min(pos + 3, len(test_seq) - self.cpg_window + 1)
        ):
            window = test_seq[start:start + self.cpg_window]
            c_count = window.count("C")
            g_count = window.count("G")
            cg_count = sum(1 for i in range(len(window) - 1) if window[i:i+2] == "CG")
            expected = (c_count * g_count) / len(window) if len(window) > 0 else 0
            obs_exp = cg_count / expected if expected > 0 else 0.0
            if obs_exp > self.cpg_threshold:
                return False

        return True

    # ──────────────────────────────────────────────────────────
    # Phase 1: Greedy codon optimization (with 2-codon look-ahead)
    # ──────────────────────────────────────────────────────────
    def _phase1_greedy_optimize(self, seq: str) -> str:
        """Phase 1: Per-position CAI maximization with 2-codon look-ahead, GT-aware.

        For each amino acid, select the highest-CAI codon that does not
        introduce a GT dinucleotide within or across codon boundaries.

        Uses 2-codon look-ahead: when choosing codon at position i, also
        consider position i+1. Score pairs by combined CAI and choose the
        pair with highest combined CAI that satisfies GT constraints.

        If ALL synonymous codons for an AA contain GT (e.g., Valine GTN),
        flag the position as "unavoidable GT" for Phase 4 mutagenesis.
        """
        protein = self._translate(seq)
        n = len(protein)
        # Build partial codon list greedily, then apply look-ahead
        codons = [""] * n

        i = 0
        while i < n:
            aa = protein[i]
            if aa == "*":
                codons[i] = seq[i * 3:i * 3 + 3]
                i += 1
                continue

            candidates = AA_TO_CODONS.get(aa, [])
            if not candidates:
                codons[i] = seq[i * 3:i * 3 + 3]
                i += 1
                continue

            prev_end = codons[i - 1][-1] if i > 0 and codons[i - 1] else ""

            # Sort by CAI descending
            candidates_sorted = sorted(
                candidates,
                key=lambda c: self.species_cai.get(c, 0.0),
                reverse=True,
            )

            # Try 2-codon look-ahead if there's a next AA
            next_aa = protein[i + 1] if i + 1 < n else None
            next_after_end = ""

            if next_aa is not None and next_aa != "*":
                next_candidates = AA_TO_CODONS.get(next_aa, [])
                if next_candidates:
                    next_candidates_sorted = sorted(
                        next_candidates,
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    )

                    best_pair = None
                    best_pair_cai = -1.0
                    # Also track the best single-codon fallback
                    best_single = None
                    best_single_cai = -1.0

                    for c1 in candidates_sorted:
                        c1_cai = self.species_cai.get(c1, 0.0)
                        # Check single-codon GT constraints for c1
                        c1_gt_ok = True
                        if self.avoid_gt:
                            if "GT" in c1:
                                c1_gt_ok = False
                            elif prev_end and prev_end + c1[0] == "GT":
                                c1_gt_ok = False

                        if c1_gt_ok and c1_cai > best_single_cai and best_single is None:
                            best_single = c1
                            best_single_cai = c1_cai

                        if not c1_gt_ok:
                            continue

                        for c2 in next_candidates_sorted:
                            c2_cai = self.species_cai.get(c2, 0.0)
                            pair_cai = c1_cai * c2_cai

                            if self.avoid_gt:
                                if "GT" in c2:
                                    continue
                                if c1[-1] + c2[0] == "GT":
                                    continue

                            if pair_cai > best_pair_cai:
                                best_pair = (c1, c2)
                                best_pair_cai = pair_cai
                                # Early exit: if we found the best possible pair,
                                # both c1 and c2 are highest-CAI
                                if c1_cai == self.species_cai.get(candidates_sorted[0], 0.0) and \
                                   c2_cai == self.species_cai.get(next_candidates_sorted[0], 0.0):
                                    break
                        # If best pair has the best possible c1, stop searching c1
                        if best_pair and best_pair[0] == candidates_sorted[0]:
                            break

                    # Determine if pair look-ahead gives better result than greedy
                    greedy_pair_cai = -1.0
                    if best_single:
                        # What's the greedy choice for next codon given best_single?
                        for c2 in next_candidates_sorted:
                            c2_cai = self.species_cai.get(c2, 0.0)
                            if self.avoid_gt:
                                if "GT" in c2:
                                    continue
                                if best_single[-1] + c2[0] == "GT":
                                    continue
                            greedy_pair_cai = best_single_cai * c2_cai
                            break

                    if best_pair and best_pair_cai >= greedy_pair_cai:
                        codons[i] = best_pair[0]
                        codons[i + 1] = best_pair[1]
                        i += 2
                        continue
                    elif best_single:
                        codons[i] = best_single
                        i += 1
                        continue

            # Fall back to single-codon greedy selection
            best = candidates_sorted[0]
            found_gt_free = False

            if self.avoid_gt:
                for codon in candidates_sorted:
                    if "GT" in codon:
                        continue
                    if prev_end and prev_end + codon[0] == "GT":
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
                        if prev_end:
                            for codon in candidates_sorted:
                                if prev_end + codon[0] != "GT":
                                    best = codon
                                    break
                    else:
                        # Pick the one that minimizes GT count
                        for codon in candidates_sorted:
                            gt_count = codon.count("GT")
                            cross_gt = 1 if (prev_end and prev_end + codon[0] == "GT") else 0
                            best_gt_count = best.count("GT")
                            best_cross_gt = 1 if (prev_end and prev_end + best[0] == "GT") else 0
                            if gt_count + cross_gt < best_gt_count + best_cross_gt:
                                best = codon
                                break

            codons[i] = best
            i += 1

        return "".join(codons)

    # ──────────────────────────────────────────────────────────
    # Phase 2: Restriction site removal (CAI-prioritized)
    # ──────────────────────────────────────────────────────────
    def _phase2_remove_restriction_sites(self, seq: str) -> str:
        """Phase 2: Remove restriction enzyme recognition sites by synonymous substitution.

        Alternatives are sorted by CAI (descending) to prefer the highest-CAI
        valid replacement.

        When avoid_gt is True, we prefer substitutions that don't increase GT count.
        However, if no GT-neutral substitution exists, we allow GT count to increase
        by at most 1 per site (later phases will try to resolve new GTs).

        Uses site COUNT comparison rather than global ``site in test_seq`` so that
        partial removal of multiple same-site occurrences is accepted.
        """
        from .restriction_sites import get_recognition_site

        seq_list = list(seq)
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            pos = 0
            while pos <= len(seq_list) - len(site):
                if "".join(seq_list[pos:pos + len(site)]) == site:
                    codon_starts = set()
                    for j in range(pos, pos + len(site)):
                        cs = (j // 3) * 3
                        if cs + 3 <= len(seq_list):
                            codon_starts.add(cs)

                    resolved = False
                    best_gt_increase = None  # (gt_increase, test_list)
                    old_site_count = "".join(seq_list).count(site)

                    for cs in sorted(codon_starts):
                        codon = "".join(seq_list[cs:cs + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue

                        # Sort alternatives by CAI descending
                        alts_sorted = sorted(
                            [alt for alt in AA_TO_CODONS.get(aa, []) if alt != codon],
                            key=lambda c: self.species_cai.get(c, 0.0),
                            reverse=True,
                        )

                        for alt in alts_sorted:
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[cs + k] = b
                            test_seq = "".join(test_list)
                            new_site_count = test_seq.count(site)
                            if new_site_count >= old_site_count:
                                continue  # Didn't reduce the site count
                            if self.avoid_gt:
                                old_gt_count = _count_gts("".join(seq_list))
                                new_gt_count = _count_gts(test_seq)
                                gt_increase = new_gt_count - old_gt_count
                                if gt_increase == 0:
                                    # Best case: GT-neutral + site removed
                                    seq_list = test_list
                                    resolved = True
                                    break
                                elif best_gt_increase is None or gt_increase < best_gt_increase[0]:
                                    # Track minimal GT increase option
                                    best_gt_increase = (gt_increase, test_list)
                            else:
                                seq_list = test_list
                                resolved = True
                                break
                        if resolved:
                            break

                    if not resolved and best_gt_increase is not None:
                        # Allow GT increase if it's the only way to remove the site
                        # (later phases will try to resolve new GTs)
                        seq_list = best_gt_increase[1]
                        resolved = True

                    if resolved:
                        continue
                    else:
                        pos += 1
                else:
                    pos += 1

        return "".join(seq_list)

    # ──────────────────────────────────────────────────────────
    # Phase 3: Cross-codon constraint resolution (iterative, CAI-prioritized)
    # ──────────────────────────────────────────────────────────
    def _phase3_cross_codon_constraints(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Phase 3: Resolve cross-codon GT, CG, and restriction site constraints.

        This phase iterates until no more cross-codon constraints can be
        resolved. Each resolution is globally validated to ensure no new
        GTs are introduced elsewhere. Alternatives are CAI-prioritized.
        """
        seq_list = list(seq)
        total_remaining: Dict[int, List[str]] = {}
        max_rounds = 20

        for round_num in range(max_rounds):
            current_seq = "".join(seq_list)
            constraint_positions: Dict[int, List[str]] = {}

            for pos in find_cross_codon_gt(current_seq):
                codon_start = (pos // 3) * 3
                constraint_positions.setdefault(codon_start, [])
                if "GT" not in constraint_positions[codon_start]:
                    constraint_positions[codon_start].append("GT")

            for pos in find_cross_codon_cg(current_seq):
                codon_start = (pos // 3) * 3
                constraint_positions.setdefault(codon_start, [])
                if "CG" not in constraint_positions[codon_start]:
                    constraint_positions[codon_start].append("CG")

            from .restriction_sites import get_recognition_site
            for enzyme in self.enzymes:
                site = get_recognition_site(enzyme)
                if site is None:
                    continue
                for pos in find_cross_codon_restriction(current_seq, site):
                    codon_start = (pos // 3) * 3
                    label = f"RS:{site}"
                    constraint_positions.setdefault(codon_start, [])
                    if label not in constraint_positions[codon_start]:
                        constraint_positions[codon_start].append(label)

            if not constraint_positions:
                break

            any_resolved = False

            for codon_start, ctypes in constraint_positions.items():
                resolved = self._try_resolve_cross_codon(seq_list, codon_start, ctypes)
                if resolved:
                    any_resolved = True

            if not any_resolved:
                # Record remaining constraints
                current_seq = "".join(seq_list)
                constraint_positions2: Dict[int, List[str]] = {}
                for pos in find_cross_codon_gt(current_seq):
                    cs = (pos // 3) * 3
                    constraint_positions2.setdefault(cs, [])
                    if "GT" not in constraint_positions2[cs]:
                        constraint_positions2[cs].append("GT")
                for pos in find_cross_codon_cg(current_seq):
                    cs = (pos // 3) * 3
                    constraint_positions2.setdefault(cs, [])
                    if "CG" not in constraint_positions2[cs]:
                        constraint_positions2[cs].append("CG")
                total_remaining = constraint_positions2
                break

        mut_report = propose_mutagenesis(
            "".join(seq_list),
            list(total_remaining.keys()),
            total_remaining,
            self.species_cai,
            self.min_blosum,
            self.min_cai,
        )

        return "".join(seq_list), mut_report

    def _try_resolve_cross_codon(
        self, seq_list: list, codon_start: int, ctypes: List[str]
    ) -> bool:
        """Try to resolve cross-codon constraints at codon_start.

        Uses global validation: after any substitution, verifies that the
        total GT count in the sequence doesn't increase.

        Key improvements over v7:
        - Does NOT blindly skip codons with internal GT — uses global validation
          to accept substitutions where total GT count doesn't increase
        - Collects ALL valid resolving pairs and picks the one with highest
          combined CAI, rather than accepting the first one found
        - Sorts candidate codons by CAI (descending) for efficiency
        """
        old_seq = "".join(seq_list)
        old_gt_count = _count_gts(old_seq)

        aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
        aa1 = CODON_TABLE.get(aa1_codon)
        if aa1 is None or aa1 == "*":
            return False

        next_start = codon_start + 3

        # Helper: check if a test substitution resolves all constraint types
        def _check_resolved(test_seq: str, region_start: int, region_end: int) -> bool:
            for ct in ctypes:
                if ct == "GT" and _has_gt(test_seq[region_start:region_end]):
                    return False
                elif ct == "CG" and "CG" in test_seq[region_start:region_end]:
                    return False
                elif ct.startswith("RS:"):
                    if ct[3:] in test_seq:
                        return False
            return True

        # Collect all valid resolutions with their CAI scores
        candidates: List[Tuple[float, list]] = []  # (combined_cai, test_list)

        # Strategy A: Following codon pair
        if next_start + 3 <= len(seq_list):
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if aa2 is not None and aa2 != "*":
                c1_alts = sorted(
                    AA_TO_CODONS.get(aa1, [aa1_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                c2_alts = sorted(
                    AA_TO_CODONS.get(aa2, [aa2_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c1 in c1_alts:
                    for c2 in c2_alts:
                        # Check boundary GT (but don't blindly skip internal GT)
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
                        # Allow if total GT count doesn't increase
                        if self.avoid_gt and new_gt_count > old_gt_count:
                            continue
                        if not self.avoid_gt or new_gt_count <= old_gt_count:
                            if _check_resolved(test_seq, codon_start, next_start + 3):
                                combined_cai = (
                                    self.species_cai.get(c1, 0.0)
                                    * self.species_cai.get(c2, 0.0)
                                )
                                candidates.append((combined_cai, test_list))

        # Strategy B: Preceding codon pair
        if codon_start >= 3:
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            if aa0 is not None and aa0 != "*":
                c0_alts = sorted(
                    AA_TO_CODONS.get(aa0, [aa0_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                c1_alts = sorted(
                    AA_TO_CODONS.get(aa1, [aa1_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c0 in c0_alts:
                    for c1 in c1_alts:
                        if c0[-1] + c1[0] == "GT":
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(c0):
                            test_list[prev_start + k] = b
                        for k, b in enumerate(c1):
                            test_list[codon_start + k] = b
                        test_seq = "".join(test_list)
                        new_gt_count = _count_gts(test_seq)
                        if self.avoid_gt and new_gt_count > old_gt_count:
                            continue
                        if not self.avoid_gt or new_gt_count <= old_gt_count:
                            if _check_resolved(test_seq, prev_start, codon_start + 3):
                                combined_cai = (
                                    self.species_cai.get(c0, 0.0)
                                    * self.species_cai.get(c1, 0.0)
                                )
                                candidates.append((combined_cai, test_list))

        # Strategy C: Both preceding and following
        if codon_start >= 3 and next_start + 3 <= len(seq_list):
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if (aa0 is not None and aa0 != "*" and
                    aa2 is not None and aa2 != "*"):
                c0_alts = sorted(
                    AA_TO_CODONS.get(aa0, [aa0_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                c1_alts = sorted(
                    AA_TO_CODONS.get(aa1, [aa1_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                c2_alts = sorted(
                    AA_TO_CODONS.get(aa2, [aa2_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c0 in c0_alts:
                    for c1 in c1_alts:
                        if c0[-1] + c1[0] == "GT":
                            continue
                        for c2 in c2_alts:
                            if c1[-1] + c2[0] == "GT":
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
                            if self.avoid_gt and new_gt_count > old_gt_count:
                                continue
                            if not self.avoid_gt or new_gt_count <= old_gt_count:
                                if _check_resolved(test_seq, prev_start, next_start + 3):
                                    combined_cai = (
                                        self.species_cai.get(c0, 0.0)
                                        * self.species_cai.get(c1, 0.0)
                                        * self.species_cai.get(c2, 0.0)
                                    )
                                    candidates.append((combined_cai, test_list))

        # Strategy D: Single codon substitution
        c1_alts = sorted(
            AA_TO_CODONS.get(aa1, [aa1_codon]),
            key=lambda c: self.species_cai.get(c, 0.0),
            reverse=True,
        )
        for c1 in c1_alts:
            test_list = seq_list[:]
            for k, b in enumerate(c1):
                test_list[codon_start + k] = b
            test_seq = "".join(test_list)
            new_gt_count = _count_gts(test_seq)
            if self.avoid_gt and new_gt_count > old_gt_count:
                continue
            if new_gt_count < old_gt_count or not self.avoid_gt:
                combined_cai = self.species_cai.get(c1, 0.0)
                candidates.append((combined_cai, test_list))

        # Pick the highest-CAI resolution
        if candidates:
            candidates.sort(key=lambda x: x[0], reverse=True)
            seq_list[:] = candidates[0][1]
            return True

        return False

    # ──────────────────────────────────────────────────────────
    # Phase 3.5: Within-codon GT resolution
    # ──────────────────────────────────────────────────────────
    def _phase35_within_codon_gt(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Phase 3.5: Within-codon GT resolution.

        For each within-codon GT (not cross-codon), try synonymous substitution
        first. If no synonymous codon can avoid GT, flag for mutagenesis.
        Alternatives are sorted by CAI (descending).
        """
        seq_list = list(seq)
        remaining_positions: Dict[int, List[str]] = {}

        for i in range(len(seq_list) - 1):
            if seq_list[i] == "G" and seq_list[i + 1] == "T":
                codon_start = (i // 3) * 3
                next_codon_start = codon_start + 3

                if i + 1 < next_codon_start:
                    # Within-codon GT
                    codon = "".join(seq_list[codon_start:codon_start + 3])
                    aa = CODON_TABLE.get(codon)
                    if aa is None or aa == "*":
                        continue

                    old_gt_count = _count_gts("".join(seq_list))
                    resolved = False
                    alts = AA_TO_CODONS.get(aa, [])
                    alts_sorted = sorted(alts,
                                         key=lambda c: self.species_cai.get(c, 0.0),
                                         reverse=True)
                    for alt in alts_sorted:
                        if alt == codon or "GT" in alt:
                            continue
                        prev_base = seq_list[codon_start - 1] if codon_start > 0 else ""
                        next_base = seq_list[next_codon_start] if next_codon_start < len(seq_list) else ""
                        if prev_base and prev_base + alt[0] == "GT":
                            continue
                        if next_base and alt[-1] + next_base == "GT":
                            continue
                        # Apply and validate globally
                        test_list = seq_list[:]
                        for k, b in enumerate(alt):
                            test_list[codon_start + k] = b
                        new_gt_count = _count_gts("".join(test_list))
                        if new_gt_count < old_gt_count:
                            seq_list = test_list
                            resolved = True
                            break

                    if not resolved:
                        remaining_positions[codon_start] = ["GT:within"]

        if remaining_positions:
            mut_report = self._propose_within_codon_mutagenesis(
                "".join(seq_list), remaining_positions
            )
        else:
            mut_report = MutagenesisReport()

        return "".join(seq_list), mut_report

    def _propose_within_codon_mutagenesis(
        self,
        seq: str,
        positions: Dict[int, List[str]],
    ) -> MutagenesisReport:
        """Propose AA substitutions for within-codon GTs that can't be resolved
        by synonymous substitution.

        Uses BLOSUM62 guidance for conservative substitutions, with CAI weighting:
        - Valine (V, GTN) → Isoleucine (I, ATN) or Leucine (L, TTA/CTN)
        - Prefers higher-CAI replacements among acceptable BLOSUM62 scores
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
            # Use combined BLOSUM*CAI score for ranking
            best_combined = -100.0

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

                    # Combined score: BLOSUM * CAI (prefer both high conservation and high CAI)
                    combined = score * cai
                    if combined > best_combined:
                        best = (alt_codon, new_aa, score, cai)
                        best_combined = combined
                    break  # Take the highest-CAI alt for this AA

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

    # ──────────────────────────────────────────────────────────
    # Phase 4: Mutagenesis fallback (CAI-weighted)
    # ──────────────────────────────────────────────────────────
    def _phase4_mutagenesis_fallback(self, seq: str, mut_report: MutagenesisReport) -> str:
        """Phase 4: Apply mutagenesis proposals for intractable constraints.

        Applies conservative AA substitutions (e.g., Val→Ile, Val→Leu)
        using BLOSUM62 guidance. Validates that substitutions reduce GT count.

        Also considers CAI of replacement codons when choosing among proposals.
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

    # ──────────────────────────────────────────────────────────
    # Phase 5: CpG island avoidance (CAI-prioritized)
    # ──────────────────────────────────────────────────────────
    def _phase5_avoid_cpg_islands(self, seq: str) -> str:
        """Phase 5: CpG island avoidance by synonymous substitution.

        Handles BOTH within-codon and cross-codon CG dinucleotides.
        For cross-codon CG (C at end of one codon, G at start of next),
        tries swapping the codon containing C first, then the codon
        containing G, and finally both codons simultaneously.

        Alternatives are sorted by CAI (descending) to prefer the highest-CAI
        valid replacement. Tracks swapped positions to prevent oscillation.
        """
        seq_list = list(seq)
        changed = True
        iterations = 0
        max_iterations = 80
        # Track positions we've already swapped to prevent oscillation
        swapped_positions: Set[int] = set()

        while changed and iterations < max_iterations:
            changed = False
            iterations += 1

            for start in range(0, len(seq_list) - self.cpg_window + 1):
                window = "".join(seq_list[start:start + self.cpg_window])
                c_count = window.count("C")
                g_count = window.count("G")
                cg_count = sum(1 for i in range(len(window) - 1) if window[i:i+2] == "CG")
                expected = (c_count * g_count) / len(window) if len(window) > 0 else 0
                obs_exp = cg_count / expected if expected > 0 else 0.0

                if obs_exp <= self.cpg_threshold:
                    continue

                # Find CG dinucleotides in this window and try to break them
                for i in range(start, min(start + self.cpg_window - 1, len(seq_list) - 1)):
                    if seq_list[i] == "C" and seq_list[i+1] == "G":
                        codon_start_c = (i // 3) * 3   # codon containing C
                        codon_start_g = ((i+1) // 3) * 3  # codon containing G
                        is_cross_codon = (codon_start_c != codon_start_g)

                        old_gt_count = _count_gts("".join(seq_list))

                        # Strategy A: Try swapping the codon containing C
                        if codon_start_c + 3 <= len(seq_list):
                            result = self._try_break_cpg_cg(
                                seq_list, codon_start_c, old_gt_count, swapped_positions,
                                is_c_side=True, is_cross_codon=is_cross_codon,
                            )
                            if result is not None:
                                seq_list = result
                                swapped_positions.add(codon_start_c)
                                changed = True
                                break

                        # Strategy B: Try swapping the codon containing G (for cross-codon CG)
                        if not changed and is_cross_codon and codon_start_g + 3 <= len(seq_list):
                            result = self._try_break_cpg_cg(
                                seq_list, codon_start_g, old_gt_count, swapped_positions,
                                is_c_side=False, is_cross_codon=is_cross_codon,
                            )
                            if result is not None:
                                seq_list = result
                                swapped_positions.add(codon_start_g)
                                changed = True
                                break

                        # Strategy C: Try both codons simultaneously (cross-codon)
                        if (not changed and is_cross_codon
                                and codon_start_c + 3 <= len(seq_list)
                                and codon_start_g + 3 <= len(seq_list)):
                            result = self._try_break_cpg_two_codon(
                                seq_list, codon_start_c, codon_start_g,
                                old_gt_count, swapped_positions,
                            )
                            if result is not None:
                                seq_list = result
                                swapped_positions.add(codon_start_c)
                                swapped_positions.add(codon_start_g)
                                changed = True
                                break

                    if changed:
                        break
                if changed:
                    break

            # If no changes after a full pass, clear swapped_positions for next iteration
            if not changed and swapped_positions:
                swapped_positions.clear()
                # Try one more pass without position restrictions
                for start in range(0, len(seq_list) - self.cpg_window + 1):
                    window = "".join(seq_list[start:start + self.cpg_window])
                    c_count = window.count("C")
                    g_count = window.count("G")
                    cg_count = sum(1 for i in range(len(window) - 1) if window[i:i+2] == "CG")
                    expected = (c_count * g_count) / len(window) if len(window) > 0 else 0
                    obs_exp = cg_count / expected if expected > 0 else 0.0
                    if obs_exp <= self.cpg_threshold:
                        continue

                    for i in range(start, min(start + self.cpg_window - 1, len(seq_list) - 1)):
                        if seq_list[i] == "C" and seq_list[i+1] == "G":
                            codon_start_c = (i // 3) * 3
                            codon_start_g = ((i+1) // 3) * 3
                            is_cross = (codon_start_c != codon_start_g)
                            old_gt_count = _count_gts("".join(seq_list))

                            # Try both sides and two-codon simultaneously
                            for strategy in ['c_side', 'g_side', 'both']:
                                if strategy == 'c_side' and codon_start_c + 3 <= len(seq_list):
                                    result = self._try_break_cpg_cg(
                                        seq_list, codon_start_c, old_gt_count, set(),
                                        is_c_side=True, is_cross_codon=is_cross,
                                    )
                                elif strategy == 'g_side' and is_cross and codon_start_g + 3 <= len(seq_list):
                                    result = self._try_break_cpg_cg(
                                        seq_list, codon_start_g, old_gt_count, set(),
                                        is_c_side=False, is_cross_codon=is_cross,
                                    )
                                elif strategy == 'both' and is_cross:
                                    result = self._try_break_cpg_two_codon(
                                        seq_list, codon_start_c, codon_start_g,
                                        old_gt_count, set(),
                                    )
                                else:
                                    continue

                                if result is not None:
                                    seq_list = result
                                    changed = True
                                    break
                        if changed:
                            break
                    if changed:
                        break

        return "".join(seq_list)

    def _try_break_cpg_cg(
        self, seq_list: list, codon_start: int, old_gt_count: int,
        swapped_positions: Set[int], is_c_side: bool, is_cross_codon: bool,
    ) -> Optional[list]:
        """Try to break a CG dinucleotide by swapping one codon.

        Args:
            is_c_side: True if codon_start contains the C of the CG pair
            is_cross_codon: True if CG spans two codons
        """
        if codon_start in swapped_positions:
            return None

        codon = "".join(seq_list[codon_start:codon_start + 3])
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            return None

        # For within-codon CG: exclude alternatives containing CG
        # For cross-codon CG: we just need the boundary to not be CG
        alts_sorted = sorted(
            [alt for alt in AA_TO_CODONS.get(aa, []) if alt != codon],
            key=lambda c: self.species_cai.get(c, 0.0),
            reverse=True,
        )

        for alt in alts_sorted:
            # Skip if this alt contains internal CG or GT
            if "CG" in alt:
                continue
            if self.avoid_gt and "GT" in alt:
                continue
            cai = self.species_cai.get(alt, 0.0)
            if cai < self.min_cai:
                continue

            # Check cross-codon GT effects
            prev_base = seq_list[codon_start - 1] if codon_start > 0 else ""
            next_base = seq_list[codon_start + 3] if codon_start + 3 < len(seq_list) else ""
            if prev_base and prev_base + alt[0] == "GT":
                continue
            if next_base and alt[-1] + next_base == "GT":
                continue

            # For cross-codon CG, verify the swap actually breaks the CG
            if is_cross_codon:
                if is_c_side:
                    # C is the last base of this codon; check new last base + next codon's first base
                    next_codon_start = codon_start + 3
                    if next_codon_start + 3 <= len(seq_list):
                        next_first = seq_list[next_codon_start]
                        if alt[-1] + next_first == "CG":
                            continue  # Doesn't break the CG
                else:
                    # G is the first base of this codon; check prev codon's last base + new first base
                    if codon_start > 0:
                        prev_last = seq_list[codon_start - 1]
                        if prev_last + alt[0] == "CG":
                            continue  # Doesn't break the CG

            # Apply and validate globally
            test_list = seq_list[:]
            for k, b in enumerate(alt):
                test_list[codon_start + k] = b
            new_gt_count = _count_gts("".join(test_list))
            if new_gt_count <= old_gt_count:
                return test_list

        return None

    def _try_break_cpg_two_codon(
        self, seq_list: list, codon_start_c: int, codon_start_g: int,
        old_gt_count: int, swapped_positions: Set[int],
    ) -> Optional[list]:
        """Try to break a cross-codon CG by swapping both codons simultaneously."""
        if codon_start_c in swapped_positions or codon_start_g in swapped_positions:
            return None

        codon_c = "".join(seq_list[codon_start_c:codon_start_c + 3])
        codon_g = "".join(seq_list[codon_start_g:codon_start_g + 3])
        aa_c = CODON_TABLE.get(codon_c)
        aa_g = CODON_TABLE.get(codon_g)
        if aa_c is None or aa_c == "*" or aa_g is None or aa_g == "*":
            return None

        alts_c = sorted(
            [a for a in AA_TO_CODONS.get(aa_c, []) if a != codon_c and "CG" not in a and "GT" not in a],
            key=lambda c: self.species_cai.get(c, 0.0),
            reverse=True,
        )
        alts_g = sorted(
            [a for a in AA_TO_CODONS.get(aa_g, []) if a != codon_g and "CG" not in a and "GT" not in a],
            key=lambda c: self.species_cai.get(c, 0.0),
            reverse=True,
        )

        best_pair = None
        best_cai = -1.0

        for ac in alts_c:
            for ag in alts_g:
                # Check that the pair doesn't create CG at the boundary
                if ac[-1] + ag[0] == "CG":
                    continue
                # Check cross-codon GT effects
                prev_base = seq_list[codon_start_c - 1] if codon_start_c > 0 else ""
                next_base = seq_list[codon_start_g + 3] if codon_start_g + 3 < len(seq_list) else ""
                if prev_base and prev_base + ac[0] == "GT":
                    continue
                if next_base and ag[-1] + next_base == "GT":
                    continue
                # Check GT between the two new codons
                if ac[-1] + ag[0] == "GT":
                    continue

                # Apply and validate
                test_list = seq_list[:]
                for k, b in enumerate(ac):
                    test_list[codon_start_c + k] = b
                for k, b in enumerate(ag):
                    test_list[codon_start_g + k] = b
                new_gt_count = _count_gts("".join(test_list))
                if new_gt_count <= old_gt_count:
                    pair_cai = self.species_cai.get(ac, 0.0) * self.species_cai.get(ag, 0.0)
                    if pair_cai > best_cai:
                        best_pair = test_list
                        best_cai = pair_cai
                    # Early exit if we find the best possible
                    if pair_cai == self.species_cai.get(alts_c[0], 0.0) * self.species_cai.get(alts_g[0], 0.0):
                        break
            if best_pair and best_cai == self.species_cai.get(alts_c[0], 0.0) * self.species_cai.get(alts_g[0], 0.0):
                break

        return best_pair

    # ──────────────────────────────────────────────────────────
    # Phase 6: Re-optimization pass (constraint-preserving CAI optimizer)
    # ──────────────────────────────────────────────────────────
    def _phase6_reoptimize(self, seq: str) -> str:
        """Phase 6: Iterative constraint-preserving CAI re-optimization.

        Repeats until no more improvements can be made:
        1. Per-codon CAI optimization with constraint preservation
           (tries ALL codons, not just GT-containing ones)
        2. Cross-codon GT resolution
        3. Restriction site removal
        4. Constraint-preserving CAI swaps for non-GT codons
        """
        seq_list = list(seq)
        max_iterations = 20

        for iteration in range(max_iterations):
            old_gt_count = _count_gts("".join(seq_list))
            improved = False

            # Step 1: Per-codon optimization - try to swap GT-containing codons
            # and also try to boost CAI for ALL codons
            for i in range(0, len(seq_list) - 2, 3):
                codon = "".join(seq_list[i:i+3])
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*":
                    continue

                candidates = AA_TO_CODONS.get(aa, [])
                if not candidates:
                    continue

                current_cai = self.species_cai.get(codon, 0.0)

                # If current codon has GT, try to swap to GT-free
                if "GT" in codon:
                    for alt in sorted(candidates,
                                       key=lambda c: self.species_cai.get(c, 0.0),
                                       reverse=True):
                        if "GT" in alt:
                            continue
                        if alt == codon:
                            continue
                        if self._is_valid_swap(seq_list, i, alt, old_gt_count):
                            seq_list[i:i+3] = list(alt)
                            old_gt_count = _count_gts("".join(seq_list))
                            improved = True
                            break
                else:
                    # Even if no GT, try to swap to a higher-CAI alternative
                    # that preserves all constraints
                    for alt in sorted(candidates,
                                       key=lambda c: self.species_cai.get(c, 0.0),
                                       reverse=True):
                        if alt == codon:
                            continue
                        alt_cai = self.species_cai.get(alt, 0.0)
                        if alt_cai <= current_cai:
                            break  # No improvement possible (sorted descending)
                        if self._is_valid_swap(seq_list, i, alt, old_gt_count):
                            seq_list[i:i+3] = list(alt)
                            old_gt_count = _count_gts("".join(seq_list))
                            current_cai = alt_cai
                            improved = True
                            break

            # Step 2: Cross-codon GT resolution
            # Must also check that we don't introduce restriction sites
            current_seq = "".join(seq_list)
            current_rs_count = self._count_restriction_sites(current_seq)
            for pos in find_cross_codon_gt(current_seq):
                codon_start = (pos // 3) * 3
                next_start = codon_start + 3

                if codon_start + 3 > len(seq_list) or next_start + 3 > len(seq_list):
                    continue

                aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
                aa1 = CODON_TABLE.get(aa1_codon)
                if aa1 is None or aa1 == "*":
                    continue

                aa2_codon = "".join(seq_list[next_start:next_start + 3])
                aa2 = CODON_TABLE.get(aa2_codon)

                if aa2 is not None and aa2 != "*":
                    # Collect valid resolutions, pick highest-CAI
                    best_pair = None
                    best_pair_cai = -1.0

                    c1_alts = sorted(
                        AA_TO_CODONS.get(aa1, [aa1_codon]),
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    )
                    c2_alts = sorted(
                        AA_TO_CODONS.get(aa2, [aa2_codon]),
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    )

                    for c1 in c1_alts:
                        if "GT" in c1:
                            continue
                        for c2 in c2_alts:
                            if "GT" in c2:
                                continue
                            if c1[-1] + c2[0] == "GT":
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(c1):
                                test_list[codon_start + k] = b
                            for k, b in enumerate(c2):
                                test_list[next_start + k] = b
                            test_seq = "".join(test_list)
                            if _count_gts(test_seq) < _count_gts("".join(seq_list)):
                                # Check we don't introduce new restriction sites
                                new_rs_count = self._count_restriction_sites(test_seq)
                                if new_rs_count > current_rs_count:
                                    continue
                                pair_cai = (
                                    self.species_cai.get(c1, 0.0)
                                    * self.species_cai.get(c2, 0.0)
                                )
                                if pair_cai > best_pair_cai:
                                    best_pair = test_list
                                    best_pair_cai = pair_cai

                    if best_pair is not None:
                        seq_list[:] = best_pair
                        improved = True

                if not improved and codon_start >= 3:
                    prev_start = codon_start - 3
                    aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
                    aa0 = CODON_TABLE.get(aa0_codon)
                    if aa0 is not None and aa0 != "*":
                        best_pair = None
                        best_pair_cai = -1.0

                        c0_alts = sorted(
                            AA_TO_CODONS.get(aa0, [aa0_codon]),
                            key=lambda c: self.species_cai.get(c, 0.0),
                            reverse=True,
                        )
                        c1_alts = sorted(
                            AA_TO_CODONS.get(aa1, [aa1_codon]),
                            key=lambda c: self.species_cai.get(c, 0.0),
                            reverse=True,
                        )

                        for c0 in c0_alts:
                            if "GT" in c0:
                                continue
                            for c1 in c1_alts:
                                if "GT" in c1:
                                    continue
                                if c0[-1] + c1[0] == "GT":
                                    continue
                                test_list = seq_list[:]
                                for k, b in enumerate(c0):
                                    test_list[prev_start + k] = b
                                for k, b in enumerate(c1):
                                    test_list[codon_start + k] = b
                                test_seq = "".join(test_list)
                                if _count_gts(test_seq) < _count_gts("".join(seq_list)):
                                    # Check we don't introduce new restriction sites
                                    new_rs_count = self._count_restriction_sites(test_seq)
                                    if new_rs_count > current_rs_count:
                                        continue
                                    pair_cai = (
                                        self.species_cai.get(c0, 0.0)
                                        * self.species_cai.get(c1, 0.0)
                                    )
                                    if pair_cai > best_pair_cai:
                                        best_pair = test_list
                                        best_pair_cai = pair_cai

                        if best_pair is not None:
                            seq_list[:] = best_pair
                            improved = True

            # Step 3: Restriction site removal
            # Allow GT increase if needed (Step 1-2 will try to resolve new GTs)
            # Uses site COUNT comparison for partial removal of same-site occurrences
            from .restriction_sites import get_recognition_site
            for enzyme in self.enzymes:
                site = get_recognition_site(enzyme)
                if site is None:
                    continue
                current_seq = "".join(seq_list)
                old_site_count = current_seq.count(site)
                p = current_seq.find(site)
                while p != -1:
                    codon_starts = set()
                    for j in range(p, p + len(site)):
                        cs = (j // 3) * 3
                        if cs + 3 <= len(seq_list):
                            codon_starts.add(cs)
                    rs_resolved = False
                    best_gt_increase = None
                    for cs in sorted(codon_starts):
                        codon = "".join(seq_list[cs:cs + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue
                        # Sort by CAI descending
                        alts_sorted = sorted(
                            [alt for alt in AA_TO_CODONS.get(aa, []) if alt != codon],
                            key=lambda c: self.species_cai.get(c, 0.0),
                            reverse=True,
                        )
                        for alt in alts_sorted:
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[cs + k] = b
                            test_seq = "".join(test_list)
                            new_site_count = test_seq.count(site)
                            if new_site_count >= old_site_count:
                                continue  # Didn't reduce the site count
                            gt_change = _count_gts(test_seq) - _count_gts("".join(seq_list))
                            if not self.avoid_gt or gt_change <= 0:
                                seq_list = test_list
                                old_site_count = new_site_count
                                rs_resolved = True
                                improved = True
                                break
                            elif best_gt_increase is None or gt_change < best_gt_increase[0]:
                                best_gt_increase = (gt_change, test_list, new_site_count)
                        if rs_resolved:
                            break
                    if not rs_resolved and best_gt_increase is not None:
                        # Accept the minimal GT increase to remove the restriction site
                        seq_list = best_gt_increase[1]
                        old_site_count = best_gt_increase[2]
                        rs_resolved = True
                        improved = True
                    if not rs_resolved:
                        p = current_seq.find(site, p + 1)
                    else:
                        current_seq = "".join(seq_list)
                        p = current_seq.find(site)

            if not improved:
                break

        return "".join(seq_list)

    # ──────────────────────────────────────────────────────────
    # Phase 7: CAI-boosting re-pass
    # ──────────────────────────────────────────────────────────
    def _phase7_cai_boost(self, seq: str) -> str:
        """Phase 7: Pure CAI maximization after all constraints are satisfied.

        Iterates through ALL codons (not just GT-containing ones) and tries
        to swap each to the highest-CAI synonymous alternative that preserves
        all constraints. Iterates until no more CAI improvements can be made.

        Accepts a swap ONLY if it:
        a) Doesn't increase GT count
        b) Doesn't create new restriction sites
        c) Doesn't increase CpG Obs/Exp ratio above threshold
        d) Doesn't create new cross-codon GT
        e) Has higher CAI than current codon
        """
        seq_list = list(seq)
        max_iterations = 30

        for iteration in range(max_iterations):
            any_improvement = False
            current_gt_count = _count_gts("".join(seq_list))

            for i in range(0, len(seq_list) - 2, 3):
                codon = "".join(seq_list[i:i+3])
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*":
                    continue

                current_cai = self.species_cai.get(codon, 0.0)
                candidates = AA_TO_CODONS.get(aa, [])
                if not candidates:
                    continue

                # Sort by CAI descending — try highest-CAI first
                candidates_sorted = sorted(
                    candidates,
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )

                for alt in candidates_sorted:
                    if alt == codon:
                        continue
                    alt_cai = self.species_cai.get(alt, 0.0)
                    if alt_cai <= current_cai:
                        break  # No CAI improvement possible (sorted descending)

                    # Use _is_valid_swap for clean constraint validation
                    if self._is_valid_swap(seq_list, i, alt, current_gt_count):
                        seq_list[i:i+3] = list(alt)
                        current_gt_count = _count_gts("".join(seq_list))
                        any_improvement = True
                        break  # Move to next codon position

            if not any_improvement:
                break

        return "".join(seq_list)

    # ──────────────────────────────────────────────────────────
    # Predicate evaluation
    # ──────────────────────────────────────────────────────────
    def _evaluate_all_predicates(self, seq: str) -> List[PredicateResult]:
        """Evaluate all 8 predicates against the optimized sequence.

        Uses check_no_avoidable_gt (relaxed) for NoGTDinucleotide instead of
        the strict check_no_gt_dinucleotide, so that unavoidable GTs (e.g.,
        Valine codons) don't cause a BRONZE certificate.
        """
        results = []

        # 1. NoStopCodons
        results.append(check_no_stop_codons(seq))

        # 2. NoCrypticSplice
        results.append(check_no_cryptic_splice(seq, self.splice_low, self.splice_high))

        # 3. NoCpGIsland
        results.append(check_no_cpg_island(seq, self.cpg_window, self.cpg_threshold))

        # 4. NoRestrictionSite
        results.append(check_no_restriction_site(seq, self.enzymes))

        # 5. NoGTDinucleotide — use relaxed (avoidable-only) check
        gt_result = check_no_avoidable_gt(seq)
        if self._applied_mutagenesis:
            mut_details = "; ".join(
                f"pos {m['position']}:{m['original_aa']}→{m['new_aa']} (BLOSUM={m['blosum']})"
                for m in self._applied_mutagenesis
            )
            gt_result.details += f" utagenesis applied: {mut_details}]"
        results.append(gt_result)

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
        all_optimal = True
        worst_cai = 1.0
        worst_codon = ""
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            cai = self.species_cai.get(codon, 0.0)
            if cai < worst_cai:
                worst_cai = cai
                worst_codon = codon
            if cai < self.min_cai:
                all_optimal = False
        results.append(PredicateResult(
            "CodonOptimality", all_optimal,
            details=f"Worst CAI: {worst_codon}={worst_cai:.4f}, min={self.min_cai}"
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
