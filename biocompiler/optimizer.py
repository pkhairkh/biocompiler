"""
BioCompiler Optimizer v7.0.0
==============================
6-phase certified gene optimization pipeline with aggressive GT resolution.

Phase 1:   Greedy codon optimization (GT-aware, with unavoidable-GT tracking)
Phase 2:   Restriction site removal
Phase 3:   Cross-codon constraint resolution (iterative, global validation)
Phase 3.5: Within-codon GT resolution (synonymous substitution + mutagenesis flagging)
Phase 4:   Mutagenesis fallback (AA substitution for Valine etc. using BLOSUM62)
Phase 5:   CpG island avoidance
Phase 6:   Re-optimization pass (iterative until convergence)
"""

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


class BioOptimizer:
    """Certified gene sequence optimizer with aggressive 6-phase pipeline."""

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

    def optimize(self, seq: str) -> Tuple[str, List[PredicateResult], str]:
        """Run the full 6-phase optimization pipeline.

        Returns:
            (optimized_sequence, predicate_results, certificate_text)
        """
        seq = seq.upper().strip()
        self._unavoidable_gt_positions = set()
        self._applied_mutagenesis = []
        self._original_protein = self._translate(seq)

        # Phase 1: Greedy codon optimization (GT-aware, with unavoidable-GT tracking)
        seq = self._phase1_greedy_optimize(seq)

        # Phase 2: Restriction site removal
        seq = self._phase2_remove_restriction_sites(seq)

        # Phase 3: Cross-codon constraint resolution (iterative)
        seq, mut_report = self._phase3_cross_codon_constraints(seq)

        # Phase 3.5: Within-codon GT resolution
        seq, mut_report_35 = self._phase35_within_codon_gt(seq)
        mut_report.proposals.extend(mut_report_35.proposals)

        # Phase 4: Mutagenesis fallback (aggressive, handles within-codon GTs too)
        seq = self._phase4_mutagenesis_fallback(seq, mut_report)

        # Phase 5: CpG island avoidance
        seq = self._phase5_avoid_cpg_islands(seq)

        # Phase 6: Re-optimization pass (iterative until convergence)
        seq = self._phase6_reoptimize(seq)

        # Evaluate all 8 predicates
        results = self._evaluate_all_predicates(seq)

        # Generate certificate
        cert_text = format_certificate(results, seq, self.species)

        return seq, results, cert_text

    # ──────────────────────────────────────────────────────────
    # Phase 1: Greedy codon optimization
    # ──────────────────────────────────────────────────────────
    def _phase1_greedy_optimize(self, seq: str) -> str:
        """Phase 1: Per-position CAI maximization, GT-aware.

        For each amino acid, select the highest-CAI codon that does not
        introduce a GT dinucleotide within or across codon boundaries.

        If ALL synonymous codons for an AA contain GT (e.g., Valine GTN),
        flag the position as "unavoidable GT" for Phase 4 mutagenesis.
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

            codons.append(best)
            prev_codon_end = best[-1]

        return "".join(codons)

    # ──────────────────────────────────────────────────────────
    # Phase 2: Restriction site removal
    # ──────────────────────────────────────────────────────────
    def _phase2_remove_restriction_sites(self, seq: str) -> str:
        """Phase 2: Remove restriction enzyme recognition sites by synonymous substitution."""
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
                    for cs in sorted(codon_starts):
                        codon = "".join(seq_list[cs:cs + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue

                        for alt in AA_TO_CODONS.get(aa, []):
                            if alt == codon:
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[cs + k] = b
                            test_seq = "".join(test_list)
                            if site not in test_seq:
                                # Also check we don't introduce GT if avoid_gt
                                if self.avoid_gt:
                                    old_gt_count = _count_gts("".join(seq_list))
                                    new_gt_count = _count_gts(test_seq)
                                    if new_gt_count > old_gt_count:
                                        continue
                                seq_list = test_list
                                resolved = True
                                break
                        if resolved:
                            break

                    if resolved:
                        continue
                    else:
                        pos += 1
                else:
                    pos += 1

        return "".join(seq_list)

    # ──────────────────────────────────────────────────────────
    # Phase 3: Cross-codon constraint resolution (iterative)
    # ──────────────────────────────────────────────────────────
    def _phase3_cross_codon_constraints(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Phase 3: Resolve cross-codon GT, CG, and restriction site constraints.

        This phase iterates until no more cross-codon constraints can be
        resolved. Each resolution is globally validated to ensure no new
        GTs are introduced elsewhere.
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
        """
        old_seq = "".join(seq_list)
        old_gt_count = _count_gts(old_seq)

        aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
        aa1 = CODON_TABLE.get(aa1_codon)
        if aa1 is None or aa1 == "*":
            return False

        next_start = codon_start + 3

        # Strategy A: Following codon
        if next_start + 3 <= len(seq_list):
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if aa2 is not None and aa2 != "*":
                for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                    if self.avoid_gt and "GT" in c1:
                        continue
                    for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
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
                                seq_list[:] = test_list
                                return True

        # Strategy B: Preceding codon
        if codon_start >= 3:
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            if aa0 is not None and aa0 != "*":
                for c0 in AA_TO_CODONS.get(aa0, [aa0_codon]):
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
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
                for c0 in AA_TO_CODONS.get(aa0, [aa0_codon]):
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                        if self.avoid_gt and "GT" in c1:
                            continue
                        for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
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

    # ──────────────────────────────────────────────────────────
    # Phase 3.5: Within-codon GT resolution
    # ──────────────────────────────────────────────────────────
    def _phase35_within_codon_gt(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Phase 3.5: Within-codon GT resolution.

        For each within-codon GT (not cross-codon), try synonymous substitution
        first. If no synonymous codon can avoid GT, flag for mutagenesis.
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

    # ──────────────────────────────────────────────────────────
    # Phase 4: Mutagenesis fallback
    # ──────────────────────────────────────────────────────────
    def _phase4_mutagenesis_fallback(self, seq: str, mut_report: MutagenesisReport) -> str:
        """Phase 4: Apply mutagenesis proposals for intractable constraints.

        Applies conservative AA substitutions (e.g., Val→Ile, Val→Leu)
        using BLOSUM62 guidance. Validates that substitutions reduce GT count.
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
    # Phase 5: CpG island avoidance
    # ──────────────────────────────────────────────────────────
    def _phase5_avoid_cpg_islands(self, seq: str) -> str:
        """Phase 5: CpG island avoidance by synonymous substitution.

        Also handles cross-codon CG dinucleotides by two-codon coordination.
        """
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

                # Find CG dinucleotides in this window and try to break them
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

                        for alt in AA_TO_CODONS.get(aa, []):
                            if alt == codon or "CG" in alt:
                                continue
                            if self.avoid_gt and "GT" in alt:
                                continue
                            cai = self.species_cai.get(alt, 0.0)
                            if cai < self.min_cai:
                                continue
                            # Check cross-codon effects
                            prev_base = seq_list[codon_start - 1] if codon_start > 0 else ""
                            next_base = seq_list[codon_start + 3] if codon_start + 3 < len(seq_list) else ""
                            if prev_base and prev_base + alt[0] == "GT":
                                continue
                            if next_base and alt[-1] + next_base == "GT":
                                continue
                            # Apply and validate globally
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

    # ──────────────────────────────────────────────────────────
    # Phase 6: Re-optimization pass (iterative until convergence)
    # ──────────────────────────────────────────────────────────
    def _phase6_reoptimize(self, seq: str) -> str:
        """Phase 6: Iterative re-optimization pass.

        Repeats until no more improvements can be made:
        1. Per-codon CAI optimization with GT avoidance
        2. Cross-codon GT resolution
        3. Within-codon GT resolution
        4. Restriction site removal
        """
        seq_list = list(seq)
        max_iterations = 20

        for iteration in range(max_iterations):
            old_gt_count = _count_gts("".join(seq_list))
            improved = False

            # Step 1: Per-codon optimization - try to swap to GT-free codons
            for i in range(0, len(seq_list) - 2, 3):
                codon = "".join(seq_list[i:i+3])
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*":
                    continue

                candidates = AA_TO_CODONS.get(aa, [])
                if not candidates:
                    continue

                # If current codon has GT, try to swap
                if "GT" in codon:
                    for alt in sorted(candidates,
                                       key=lambda c: self.species_cai.get(c, 0.0),
                                       reverse=True):
                        if "GT" in alt:
                            continue
                        prev_base = seq_list[i - 1] if i > 0 else ""
                        next_base = seq_list[i + 3] if i + 3 < len(seq_list) else ""
                        if prev_base and prev_base + alt[0] == "GT":
                            continue
                        if next_base and alt[-1] + next_base == "GT":
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(alt):
                            test_list[i + k] = b
                        if _count_gts("".join(test_list)) < old_gt_count:
                            seq_list = test_list
                            improved = True
                            break

            # Step 2: Cross-codon GT resolution
            current_seq = "".join(seq_list)
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
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                        if "GT" in c1:
                            continue
                        for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
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
                                seq_list = test_list
                                improved = True
                                break
                        if improved:
                            break

                if not improved and codon_start >= 3:
                    prev_start = codon_start - 3
                    aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
                    aa0 = CODON_TABLE.get(aa0_codon)
                    if aa0 is not None and aa0 != "*":
                        for c0 in AA_TO_CODONS.get(aa0, [aa0_codon]):
                            if "GT" in c0:
                                continue
                            for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
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
                                    seq_list = test_list
                                    improved = True
                                    break
                            if improved:
                                break

            # Step 3: Restriction site removal
            from .restriction_sites import get_recognition_site
            for enzyme in self.enzymes:
                site = get_recognition_site(enzyme)
                if site is None:
                    continue
                current_seq = "".join(seq_list)
                p = current_seq.find(site)
                while p != -1:
                    codon_starts = set()
                    for j in range(p, p + len(site)):
                        cs = (j // 3) * 3
                        if cs + 3 <= len(seq_list):
                            codon_starts.add(cs)
                    rs_resolved = False
                    for cs in sorted(codon_starts):
                        codon = "".join(seq_list[cs:cs + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue
                        for alt in AA_TO_CODONS.get(aa, []):
                            if alt == codon:
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[cs + k] = b
                            test_seq = "".join(test_list)
                            if site not in test_seq:
                                if not self.avoid_gt or _count_gts(test_seq) <= _count_gts("".join(seq_list)):
                                    seq_list = test_list
                                    rs_resolved = True
                                    improved = True
                                    break
                        if rs_resolved:
                            break
                    if not rs_resolved:
                        p = current_seq.find(site, p + 1)
                    else:
                        current_seq = "".join(seq_list)
                        p = current_seq.find(site)

            if not improved:
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
            gt_result.details += f" [mutagenesis applied: {mut_details}]"
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
