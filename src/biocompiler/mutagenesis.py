"""
BioCompiler Mutagenesis Engine v7.0.0
=======================================
Proposes conservative amino acid substitutions using BLOSUM62 scoring
to resolve cross-codon constraints (GT, CG, restriction sites).

Two failure modes:
  - "chose_poorly": A conservative substitution exists but was not found by the algorithm
  - "impossible": No substitution can resolve the constraint without violating conservation
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Set, Tuple

from .type_system import (
    CODON_TABLE, AA_TO_CODONS, BLOSUM62,
    check_conservation_score, check_codon_optimality
)


@dataclass
class MutagenesisProposal:
    """A proposed amino acid substitution to resolve a constraint."""
    position: int                # 0-based codon position
    original_codon: str
    original_aa: str
    new_aa: str
    new_codon: str
    blosum_score: int
    cai_weight: float
    resolves: List[str]          # which constraints this resolves
    chose_poorly: bool = False   # True if a better sub existed but was missed
    impossible: bool = False     # True if no sub can resolve


@dataclass
class MutagenesisReport:
    """Collection of mutagenesis proposals for a sequence."""
    proposals: List[MutagenesisProposal] = field(default_factory=list)

    def add(self, proposal: MutagenesisProposal) -> None:
        self.proposals.append(proposal)

    @property
    def has_chose_poorly(self) -> bool:
        return any(p.chose_poorly for p in self.proposals)

    @property
    def has_impossible(self) -> bool:
        return any(p.impossible for p in self.proposals)


def propose_mutagenesis(
    seq: str,
    constraint_positions: List[int],
    constraint_types: Dict[int, List[str]],
    species_cai: Dict[str, float],
    min_blosum: int = -1,
    min_cai: float = 0.0,
) -> MutagenesisReport:
    """Propose amino acid substitutions to resolve cross-codon constraints.

    For each constrained position, find the best conservative substitution
    (by BLOSUM62 score, then by CAI weight) that avoids the constraint.

    Args:
        seq: DNA sequence
        constraint_positions: Codon start positions with constraints
        constraint_types: Map of position -> list of constraint names
        species_cai: CAI weights for the target species
        min_blosum: Minimum BLOSUM62 score for acceptable substitution
        min_cai: Minimum CAI weight for the new codon

    Returns:
        MutagenesisReport with proposals
    """
    report = MutagenesisReport()

    for pos in constraint_positions:
        if pos % 3 != 0:
            # Find the codon start
            codon_start = (pos // 3) * 3
        else:
            codon_start = pos

        if codon_start + 3 > len(seq):
            continue

        original_codon = seq[codon_start:codon_start + 3]
        original_aa = CODON_TABLE.get(original_codon)

        if original_aa is None or original_aa == "*":
            continue

        # Try all synonymous codons first (no AA change)
        synonymous_candidates = []
        for alt_codon in AA_TO_CODONS.get(original_aa, []):
            if alt_codon == original_codon:
                continue
            cai = species_cai.get(alt_codon, 0.0)
            synonymous_candidates.append((alt_codon, original_aa, BLOSUM62[(original_aa, original_aa)], cai))

        # Sort synonymous by CAI (highest first)
        synonymous_candidates.sort(key=lambda x: x[3], reverse=True)

        # Check if any synonymous codon resolves the constraint
        resolved = False
        for alt_codon, aa, blosum, cai in synonymous_candidates:
            if _resolves_constraints(alt_codon, codon_start, seq, constraint_types.get(codon_start, [])):
                proposal = MutagenesisProposal(
                    position=codon_start,
                    original_codon=original_codon,
                    original_aa=original_aa,
                    new_aa=aa,
                    new_codon=alt_codon,
                    blosum_score=blosum,
                    cai_weight=cai,
                    resolves=constraint_types.get(codon_start, []),
                )
                report.add(proposal)
                resolved = True
                break

        if resolved:
            continue

        # No synonymous codon works — try conservative AA substitutions
        substitution_candidates = []
        for new_aa in BLOSUM62:
            if new_aa == original_aa:
                continue
            blosum = BLOSUM62.get((original_aa, new_aa), -10)
            if blosum < min_blosum:
                continue
            for alt_codon in AA_TO_CODONS.get(new_aa, []):
                cai = species_cai.get(alt_codon, 0.0)
                if cai < min_cai:
                    continue
                substitution_candidates.append((alt_codon, new_aa, blosum, cai))

        # Sort by BLOSUM62 (highest first), then CAI
        substitution_candidates.sort(key=lambda x: (x[2], x[3]), reverse=True)

        if not substitution_candidates:
            # Impossible: no substitution can resolve
            proposal = MutagenesisProposal(
                position=codon_start,
                original_codon=original_codon,
                original_aa=original_aa,
                new_aa="",
                new_codon="",
                blosum_score=-10,
                cai_weight=0.0,
                resolves=constraint_types.get(codon_start, []),
                impossible=True,
            )
            report.add(proposal)
            continue

        # Find the best candidate that resolves the constraint
        best = None
        for alt_codon, new_aa, blosum, cai in substitution_candidates:
            if _resolves_constraints(alt_codon, codon_start, seq, constraint_types.get(codon_start, [])):
                best = (alt_codon, new_aa, blosum, cai)
                break

        if best is None:
            # Check for "chose_poorly": is there a substitution that resolves
            # but we ranked it lower than the best that doesn't?
            for alt_codon, new_aa, blosum, cai in substitution_candidates:
                if _resolves_constraints(alt_codon, codon_start, seq, constraint_types.get(codon_start, [])):
                    proposal = MutagenesisProposal(
                        position=codon_start,
                        original_codon=original_codon,
                        original_aa=original_aa,
                        new_aa=new_aa,
                        new_codon=alt_codon,
                        blosum_score=blosum,
                        cai_weight=cai,
                        resolves=constraint_types.get(codon_start, []),
                        chose_poorly=False,
                    )
                    report.add(proposal)
                    best = (alt_codon, new_aa, blosum, cai)
                    break

        if best is None:
            proposal = MutagenesisProposal(
                position=codon_start,
                original_codon=original_codon,
                original_aa=original_aa,
                new_aa="",
                new_codon="",
                blosum_score=-10,
                cai_weight=0.0,
                resolves=constraint_types.get(codon_start, []),
                impossible=True,
            )
            report.add(proposal)
        else:
            alt_codon, new_aa, blosum, cai = best
            # Check if a better BLOSUM score was available but didn't resolve
            chose_poorly = False
            if substitution_candidates:
                top_blosum = substitution_candidates[0][2]
                if top_blosum > blosum:
                    chose_poorly = True

            proposal = MutagenesisProposal(
                position=codon_start,
                original_codon=original_codon,
                original_aa=original_aa,
                new_aa=new_aa,
                new_codon=alt_codon,
                blosum_score=blosum,
                cai_weight=cai,
                resolves=constraint_types.get(codon_start, []),
                chose_poorly=chose_poorly,
            )
            report.add(proposal)

    return report


def _resolves_constraints(
    new_codon: str, codon_start: int, seq: str, constraint_names: List[str]
) -> bool:
    """Check if substituting new_codon at codon_start resolves all listed constraints."""
    # Build the modified sequence
    seq_list = list(seq)
    for i, base in enumerate(new_codon):
        seq_list[codon_start + i] = base
    modified = "".join(seq_list)

    for cname in constraint_names:
        if cname == "GT":
            # Check if the GT that was at this position is gone
            region_start = max(0, codon_start - 1)
            region_end = min(len(modified), codon_start + 4)
            if "GT" in modified[region_start:region_end]:
                return False
        elif cname == "CG":
            region_start = max(0, codon_start - 1)
            region_end = min(len(modified), codon_start + 4)
            if "CG" in modified[region_start:region_end]:
                return False
        elif cname.startswith("RS:"):
            site = cname[3:]
            if site in modified:
                return False

    return True


# ────────────────────────────────────────────────────────────
# GT-mandatory amino acids and cryptic donor analysis
# ────────────────────────────────────────────────────────────

# Amino acids where ALL synonymous codons contain GT
# Valine (V): GTT, GTC, GTA, GTG — all start with GT
GT_MANDATORY_AAS: Set[str] = {"V"}


def find_unrepairable_cryptic_donors(
    seq: str,
    protein: str,
    organism: str = "Homo_sapiens",
    threshold: float = 3.0,
) -> List[Tuple[int, int, str, float, bool, bool]]:
    """Find cryptic splice donors that cannot be repaired by synonymous substitution.

    A cryptic splice donor is "unrepairable" if all synonymous codons for the
    amino acid at that position also contain GT, meaning no synonymous
    substitution can eliminate the GT dinucleotide. This is only possible
    for GT-mandatory amino acids like Valine (all codons: GTN).

    Args:
        seq: Optimized DNA sequence.
        protein: Target protein sequence (1-letter codes).
        organism: Target organism for CAI lookup.
        threshold: MaxEntScan score threshold for considering a GT a cryptic donor.

    Returns:
        List of tuples (position, codon_index, amino_acid, score, fixable, gt_mandatory)
        where:
        - position: 0-based nucleotide position of the GT
        - codon_index: 0-based codon index in the protein
        - amino_acid: single-letter amino acid code
        - score: MaxEntScan donor score at this position
        - fixable: True if a synonymous codon can eliminate the GT
        - gt_mandatory: True if the AA has no GT-free synonymous codons
    """
    from .maxentscan import score_donor
    from .species import SPECIES

    seq = seq.upper()
    results = []

    for i in range(len(seq) - 1):
        if seq[i:i+2] != "GT":
            continue

        # Score this GT as a donor
        donor_score = score_donor(seq, i)
        if donor_score <= -50.0:
            donor_score = 0.0

        if donor_score < threshold:
            continue

        # Determine codon context
        codon_start = (i // 3) * 3
        next_codon_start = codon_start + 3
        codon_idx = codon_start // 3

        if codon_idx >= len(protein):
            continue

        aa = protein[codon_idx]

        # Determine if within-codon or cross-codon
        is_within = (i + 1) < next_codon_start

        gt_mandatory = aa in GT_MANDATORY_AAS

        if is_within:
            # Within-codon GT: check if any synonymous codon avoids GT
            has_gt_free = any("GT" not in c for c in AA_TO_CODONS.get(aa, []))
            fixable = has_gt_free
        else:
            # Cross-codon GT: check if adjacent codons can be changed
            # to eliminate the boundary GT
            fixable = True  # Cross-codon GTs are generally fixable
            if codon_start + 3 <= len(seq):
                prev_cs = codon_start
                next_cs = next_codon_start
                if next_cs + 3 <= len(seq):
                    prev_codon = seq[prev_cs:prev_cs + 3]
                    next_codon = seq[next_cs:next_cs + 3]
                    prev_aa = CODON_TABLE.get(prev_codon)
                    next_aa = CODON_TABLE.get(next_codon)
                    if prev_aa and next_aa:
                        prev_can_avoid = any(c[-1] != 'G' for c in AA_TO_CODONS.get(prev_aa, [prev_codon]))
                        next_can_avoid = any(c[0] != 'T' for c in AA_TO_CODONS.get(next_aa, [next_codon]))
                        fixable = prev_can_avoid or next_can_avoid

        results.append((i, codon_idx, aa, donor_score, fixable, gt_mandatory))

    return results
