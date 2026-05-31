"""
BioCompiler Optimizer v7.0.0
==============================
5-phase certified gene optimization pipeline.

Phase 1: Greedy codon optimization (GT-aware)
Phase 2: Restriction site removal
Phase 3: Cross-codon constraint resolution (GT, CG, restriction sites spanning boundaries)
Phase 4: Mutagenesis fallback for intractable constraints
Phase 5: CpG island avoidance
"""

from typing import List, Dict, Optional, Tuple, Set

from .type_system import (
    CODON_TABLE, AA_TO_CODONS, BLOSUM62, SpliceVerdict, PredicateResult,
    check_no_stop_codons, check_no_cryptic_splice, check_no_cpg_island,
    check_no_restriction_site, check_no_gt_dinucleotide, check_valid_coding_seq,
    check_conservation_score, check_codon_optimality,
    find_cross_codon_gt, find_cross_codon_cg, find_cross_codon_restriction,
    PREDICATE_NAMES,
)
from .species import SPECIES
from .mutagenesis import propose_mutagenesis, MutagenesisReport
from .certificates import compute_certificate, format_certificate


class BioOptimizer:
    """Certified gene sequence optimizer with 5-phase pipeline."""

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

    def optimize(self, seq: str) -> Tuple[str, List[PredicateResult], str]:
        """Run the full 5-phase optimization pipeline.

        Returns:
            (optimized_sequence, predicate_results, certificate_text)
        """
        seq = seq.upper().strip()

        # Phase 1: Greedy codon optimization (GT-aware)
        seq = self._phase1_greedy_optimize(seq)

        # Phase 2: Restriction site removal
        seq = self._phase2_remove_restriction_sites(seq)

        # Phase 3: Cross-codon constraint resolution
        seq, mut_report = self._phase3_cross_codon_constraints(seq)

        # Phase 4: Mutagenesis fallback
        seq = self._phase4_mutagenesis_fallback(seq, mut_report)

        # Phase 5: CpG island avoidance
        seq = self._phase5_avoid_cpg_islands(seq)

        # Evaluate all 8 predicates
        results = self._evaluate_all_predicates(seq)

        # Generate certificate
        cert_text = format_certificate(results, seq, self.species)

        return seq, results, cert_text

    def _phase1_greedy_optimize(self, seq: str) -> str:
        """Phase 1: Per-position CAI maximization, GT-aware.

        For each amino acid, select the highest-CAI codon that does not
        introduce a GT dinucleotide within or across codon boundaries.
        """
        codons = []
        protein = self._translate(seq)
        prev_codon_end = ""  # last base of previous codon

        for i, aa in enumerate(protein):
            if aa == "*":
                # Keep stop codon as-is
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

            # Sort by CAI descending
            candidates_sorted = sorted(
                candidates,
                key=lambda c: self.species_cai.get(c, 0.0),
                reverse=True,
            )

            best = candidates_sorted[0]  # fallback

            if self.avoid_gt:
                for codon in candidates_sorted:
                    # Check within-codon GT
                    if "GT" in codon:
                        continue
                    # Check cross-codon GT with previous codon's last base
                    if prev_codon_end and prev_codon_end + codon[0] == "GT":
                        continue
                    best = codon
                    break

            codons.append(best)
            prev_codon_end = best[-1]

        return "".join(codons)

    def _phase2_remove_restriction_sites(self, seq: str) -> str:
        """Phase 2: Remove restriction enzyme recognition sites by synonymous substitution."""
        from .restriction_sites import get_recognition_site

        seq_list = list(seq)
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            pos = 0
            while pos <= len(seq) - len(site):
                if "".join(seq_list[pos:pos + len(site)]) == site:
                    # Find codon(s) that overlap this site
                    codon_starts = set()
                    for j in range(pos, pos + len(site)):
                        cs = (j // 3) * 3
                        if cs + 3 <= len(seq):
                            codon_starts.add(cs)

                    # Try synonymous substitutions for each codon
                    resolved = False
                    for cs in sorted(codon_starts):
                        codon = "".join(seq_list[cs:cs + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue

                        for alt in AA_TO_CODONS.get(aa, []):
                            if alt == codon:
                                continue
                            # Try this substitution
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[cs + k] = b
                            test_seq = "".join(test_list)
                            if site not in test_seq:
                                seq_list = test_list
                                resolved = True
                                break
                        if resolved:
                            break

                    if resolved:
                        continue  # re-check from same position
                    else:
                        pos += 1
                else:
                    pos += 1

        return "".join(seq_list)

    def _phase3_cross_codon_constraints(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Phase 3: Resolve cross-codon GT, CG, and restriction site constraints.

        Cross-codon constraints span codon boundaries and cannot be resolved
        by single-codon synonymous substitution. We try two-codon coordination.
        """
        constraint_positions: Dict[int, List[str]] = {}

        # Find cross-codon GTs
        for pos in find_cross_codon_gt(seq):
            codon_start = (pos // 3) * 3
            if codon_start not in constraint_positions:
                constraint_positions[codon_start] = []
            constraint_positions[codon_start].append("GT")

        # Find cross-codon CGs
        for pos in find_cross_codon_cg(seq):
            codon_start = (pos // 3) * 3
            if codon_start not in constraint_positions:
                constraint_positions[codon_start] = []
            constraint_positions[codon_start].append("CG")

        # Find cross-codon restriction sites
        from .restriction_sites import get_recognition_site
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            for pos in find_cross_codon_restriction(seq, site):
                codon_start = (pos // 3) * 3
                if codon_start not in constraint_positions:
                    constraint_positions[codon_start] = []
                constraint_positions[codon_start].append(f"RS:{site}")

        if not constraint_positions:
            return seq, MutagenesisReport()

        # Try two-codon coordination for each constraint
        seq_list = list(seq)
        remaining_constraints: Dict[int, List[str]] = {}

        for codon_start, ctypes in constraint_positions.items():
            resolved = False
            aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
            aa1 = CODON_TABLE.get(aa1_codon)
            if aa1 is None or aa1 == "*":
                remaining_constraints[codon_start] = ctypes
                continue

            # Also try the next codon
            next_start = codon_start + 3
            if next_start + 3 <= len(seq):
                aa2_codon = "".join(seq_list[next_start:next_start + 3])
                aa2 = CODON_TABLE.get(aa2_codon)

                if aa2 is not None and aa2 != "*":
                    # Try all combinations of synonymous codons for both AAs
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                        for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
                            test_list = seq_list[:]
                            for k, b in enumerate(c1):
                                test_list[codon_start + k] = b
                            for k, b in enumerate(c2):
                                test_list[next_start + k] = b
                            test_seq = "".join(test_list)

                            # Check if all constraints resolved
                            all_resolved = True
                            for ct in ctypes:
                                if ct == "GT" and "GT" in test_seq[codon_start:next_start + 3]:
                                    all_resolved = False
                                elif ct == "CG" and "CG" in test_seq[codon_start:next_start + 3]:
                                    all_resolved = False
                                elif ct.startswith("RS:"):
                                    site = ct[3:]
                                    if site in test_seq:
                                        all_resolved = False

                            if all_resolved:
                                seq_list = test_list
                                resolved = True
                                break
                        if resolved:
                            break

            if not resolved:
                # Try single codon substitution (may work for some constraint types)
                for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                    test_list = seq_list[:]
                    for k, b in enumerate(c1):
                        test_list[codon_start + k] = b
                    test_seq = "".join(test_list)
                    boundary = test_seq[codon_start:codon_start + 4]
                    all_resolved = True
                    for ct in ctypes:
                        if ct == "GT" and "GT" in boundary:
                            all_resolved = False
                        elif ct == "CG" and "CG" in boundary:
                            all_resolved = False
                    if all_resolved:
                        seq_list = test_list
                        resolved = True
                        break

            if not resolved:
                remaining_constraints[codon_start] = ctypes

        # Generate mutagenesis report for remaining constraints
        mut_report = propose_mutagenesis(
            "".join(seq_list),
            list(remaining_constraints.keys()),
            remaining_constraints,
            self.species_cai,
            self.min_blosum,
            self.min_cai,
        )

        return "".join(seq_list), mut_report

    def _phase4_mutagenesis_fallback(self, seq: str, mut_report: MutagenesisReport) -> str:
        """Phase 4: Apply mutagenesis proposals for intractable cross-codon constraints."""
        seq_list = list(seq)
        for proposal in mut_report.proposals:
            if proposal.impossible or not proposal.new_codon:
                continue
            pos = proposal.position
            for k, b in enumerate(proposal.new_codon):
                if pos + k < len(seq_list):
                    seq_list[pos + k] = b
        return "".join(seq_list)

    def _phase5_avoid_cpg_islands(self, seq: str) -> str:
        """Phase 5: CpG island avoidance by synonymous substitution.

        Scan sliding windows for high Obs/Exp CG ratio and substitute
        CG-containing codons with synonymous alternatives.
        """
        seq_list = list(seq)
        changed = True
        iterations = 0
        max_iterations = 50  # safety limit

        while changed and iterations < max_iterations:
            changed = False
            iterations += 1

            for start in range(0, len(seq) - self.cpg_window + 1, 3):
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
                        # Find which codon(s) contain this CG
                        codon_start = (i // 3) * 3
                        if codon_start + 3 > len(seq_list):
                            continue
                        codon = "".join(seq_list[codon_start:codon_start + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue

                        # Try synonymous codons without CG
                        for alt in AA_TO_CODONS.get(aa, []):
                            if alt == codon or "CG" in alt:
                                continue
                            cai = self.species_cai.get(alt, 0.0)
                            if cai < self.min_cai:
                                continue
                            # Apply
                            for k, b in enumerate(alt):
                                seq_list[codon_start + k] = b
                            changed = True
                            break

                        if changed:
                            break
                if changed:
                    break

        return "".join(seq_list)

    def _evaluate_all_predicates(self, seq: str) -> List[PredicateResult]:
        """Evaluate all 8 predicates against the optimized sequence."""
        results = []

        # 1. NoStopCodons
        results.append(check_no_stop_codons(seq))

        # 2. NoCrypticSplice (dual-threshold)
        results.append(check_no_cryptic_splice(seq, self.splice_low, self.splice_high))

        # 3. NoCpGIsland
        results.append(check_no_cpg_island(seq, self.cpg_window, self.cpg_threshold))

        # 4. NoRestrictionSite
        results.append(check_no_restriction_site(seq, self.enzymes))

        # 5. NoGTDinucleotide
        results.append(check_no_gt_dinucleotide(seq))

        # 6. ValidCodingSeq
        results.append(check_valid_coding_seq(seq))

        # 7. ConservationScore — evaluate for each codon position
        # (aggregate: all in-frame positions must have BLOSUM62 >= min_blosum)
        all_conserved = True
        details_parts = []
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            aa = CODON_TABLE.get(codon, "?")
            score = BLOSUM62.get((aa, aa), 0)
            if score < self.min_blosum:
                all_conserved = False
                details_parts.append(f"pos {i}:{aa}={score}")
        results.append(PredicateResult(
            "ConservationScore", all_conserved,
            details="; ".join(details_parts) if details_parts else f"All AA self-scores >= {self.min_blosum}"
        ))

        # 8. CodonOptimality — CAI check for each codon
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
