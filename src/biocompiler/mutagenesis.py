"""
BioCompiler Type-Directed Protein Mutagenesis Engine

When the type system proves that NO codon assignment can satisfy all
predicates for a given protein, this engine proposes conservative
amino acid substitutions that WOULD make satisfaction possible.

This is the key innovation: the type predicate doesn't just VERIFY —
it DIRECTS design across the central dogma boundary (DNA->RNA->Protein).
The derivation is not documentation; it's a repair directive.

Architecture:
    Optimizer -> Type System -> [FAIL] -> Mutagenesis Engine
                                              |
                                    Modified Protein -> Optimizer -> Type System -> [PASS]

Key insight: Valine (V) is the ONLY amino acid whose codons ALL contain
the GT dinucleotide (GTT, GTC, GTA, GTG). When a V position creates
an unrepairable cryptic splice donor, the ONLY solution is V->I (or
V->L, V->A). No codon swap can ever fix it.

BLOSUM62 scores for key substitutions:
  V->I: +3 (very conservative, similar hydrophobicity & size)
  V->L: +1 (conservative, both hydrophobic)
  V->A:  0 (moderate, Ala smaller)
  L->I: +2 (conservative, similar properties)
  L->V: +1 (conservative)
  L->M: +2 (conservative)

GT-mandatory vs optimizer-weakness distinction:
- GT-mandatory: an amino acid position where ALL synonymous codons contain
  the GT dinucleotide (e.g., Valine). No codon swap can remove GT, so
  amino acid substitution is the ONLY option.
- AG-mandatory: similarly, an amino acid position where ALL synonymous
  codons contain the AG dinucleotide.
- Optimizer-weakness: a position where the amino acid HAS codons without
  the problematic dinucleotide, but the optimizer chose a suboptimal codon.
  This is the optimizer's problem, not a mutagenesis target.

Key functions:
- is_gt_mandatory(aa): Returns True if all codons for the amino acid
  contain the GT dinucleotide (i.e., the amino acid is GT-mandatory).
- is_ag_mandatory(aa): Returns True if all codons for the amino acid
  contain the AG dinucleotide (i.e., the amino acid is AG-mandatory).
- diagnose_optimizer_weakness(sequence, protein, threshold): Identifies
  positions where the predicate fails NOT because of amino acid identity
  (GT/AG-mandatory) but because the optimizer chose a weak codon — these
  positions should be re-optimized, not substituted.

Substitution policy:
- Only propose AA substitutions for truly GT-mandatory or AG-mandatory
  positions. Optimizer-weakness positions should be handled by re-running
  the optimizer with better codon selection, not by mutagenesis.
- This avoids unnecessary protein modifications when a simple codon swap
  would suffice.

Separation of Concerns:
- Impossibility detection (find_unrepairable_*) queries the type system's
  scoring functions (MaxEntScan) but does NOT re-implement type checking.
- Substitution proposal is purely about amino acid properties (BLOSUM62).
- The optimizer handles the re-optimization loop; mutagenesis only proposes
  substitutions and applies them to the protein string.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional, Callable

from .constants import AA_TO_CODONS, CODON_TABLE
from .maxentscan import score_donor, score_acceptor, max_donor_score, max_acceptor_score
from .translation import compute_cai
from .scanner import gc_content
from .types import Verdict, TypeCheckResult
from .exceptions import InvalidProteinError

logger = logging.getLogger(__name__)


# ==============================================================================
# BLOSUM62 Substitution Matrix
# Source: Henikoff & Henikoff (1992) PNAS 89(22):10915-9
# ==============================================================================

BLOSUM62: dict[str, dict[str, int]] = {
    "A": {"A":  4, "R": -1, "N": -2, "D": -2, "C":  0, "Q": -1, "E": -1, "G":  0, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S":  1, "T":  0, "W": -3, "Y": -2, "V":  0},
    "R": {"A": -1, "R":  5, "N":  0, "D": -2, "C": -3, "Q":  1, "E":  0, "G": -2, "H":  0, "I": -3, "L": -2, "K":  2, "M": -1, "F": -3, "P": -2, "S": -1, "T": -1, "W": -3, "Y": -2, "V": -3},
    "N": {"A": -2, "R":  0, "N":  6, "D":  1, "C": -3, "Q":  0, "E":  0, "G":  0, "H":  1, "I": -3, "L": -3, "K":  0, "M": -2, "F": -3, "P": -2, "S":  1, "T":  0, "W": -4, "Y": -2, "V": -3},
    "D": {"A": -2, "R": -2, "N":  1, "D":  6, "C": -3, "Q":  0, "E":  2, "G": -1, "H": -1, "I": -3, "L": -4, "K": -1, "M": -3, "F": -3, "P": -1, "S":  0, "T": -1, "W": -4, "Y": -3, "V": -3},
    "C": {"A":  0, "R": -3, "N": -3, "D": -3, "C":  9, "Q": -3, "E": -4, "G": -3, "H": -3, "I": -1, "L": -1, "K": -3, "M": -1, "F": -2, "P": -3, "S": -1, "T": -1, "W": -2, "Y": -2, "V": -1},
    "Q": {"A": -1, "R":  1, "N":  0, "D":  0, "C": -3, "Q":  5, "E":  2, "G": -2, "H":  0, "I": -3, "L": -2, "K":  1, "M":  0, "F": -3, "P": -1, "S":  0, "T": -1, "W": -2, "Y": -1, "V": -2},
    "E": {"A": -1, "R":  0, "N":  0, "D":  2, "C": -4, "Q":  2, "E":  5, "G": -2, "H":  0, "I": -3, "L": -3, "K":  1, "M": -2, "F": -3, "P": -1, "S":  0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "G": {"A":  0, "R": -2, "N":  0, "D": -1, "C": -3, "Q": -2, "E": -2, "G":  6, "H": -2, "I": -4, "L": -4, "K": -2, "M": -3, "F": -3, "P": -2, "S":  0, "T": -2, "W": -2, "Y": -3, "V": -3},
    "H": {"A": -2, "R":  0, "N":  1, "D": -1, "C": -3, "Q":  0, "E":  0, "G": -2, "H":  8, "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2, "Y":  2, "V": -3},
    "I": {"A": -1, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -3, "E": -3, "G": -4, "H": -3, "I":  4, "L":  2, "K": -3, "M":  1, "F":  0, "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V":  3},
    "L": {"A": -1, "R": -2, "N": -3, "D": -4, "C": -1, "Q": -2, "E": -3, "G": -4, "H": -3, "I":  2, "L":  4, "K": -2, "M":  2, "F":  0, "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V":  1},
    "K": {"A": -1, "R":  2, "N":  0, "D": -1, "C": -3, "Q":  1, "E":  1, "G": -2, "H": -1, "I": -3, "L": -2, "K":  5, "M": -1, "F": -3, "P": -1, "S":  0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "M": {"A": -1, "R": -1, "N": -2, "D": -3, "C": -1, "Q":  0, "E": -2, "G": -3, "H": -2, "I":  1, "L":  2, "K": -1, "M":  5, "F":  0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V":  1},
    "F": {"A": -2, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -3, "E": -3, "G": -3, "H": -1, "I":  0, "L":  0, "K": -3, "M":  0, "F":  6, "P": -4, "S": -2, "T": -2, "W":  1, "Y":  3, "V": -1},
    "P": {"A": -1, "R": -2, "N": -2, "D": -1, "C": -3, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -4, "P":  7, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2},
    "S": {"A":  1, "R": -1, "N":  1, "D":  0, "C": -1, "Q":  0, "E":  0, "G":  0, "H": -1, "I": -2, "L": -2, "K":  0, "M": -1, "F": -2, "P": -1, "S":  4, "T":  1, "W": -3, "Y": -2, "V": -2},
    "T": {"A":  0, "R": -1, "N":  0, "D": -1, "C": -1, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S":  1, "T":  5, "W": -2, "Y": -2, "V":  0},
    "W": {"A": -3, "R": -3, "N": -4, "D": -4, "C": -2, "Q": -2, "E": -3, "G": -2, "H": -2, "I": -3, "L": -2, "K": -3, "M": -1, "F":  1, "P": -4, "S": -3, "T": -2, "W": 11, "Y":  2, "V": -3},
    "Y": {"A": -2, "R": -2, "N": -2, "D": -3, "C": -2, "Q": -1, "E": -2, "G": -3, "H":  2, "I": -1, "L": -1, "K": -2, "M": -1, "F":  3, "P": -3, "S": -2, "T": -2, "W":  2, "Y":  7, "V": -1},
    "V": {"A":  0, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -2, "E": -2, "G": -3, "H": -3, "I":  3, "L":  1, "K": -2, "M":  1, "F": -1, "P": -2, "S": -2, "T":  0, "W": -3, "Y": -1, "V":  4},
}


# ==============================================================================
# Amino Acid Properties for Substitution Guidance
# ==============================================================================

# Amino acids grouped by physicochemical properties
HYDROPHOBIC = {"A", "V", "I", "L", "M", "F", "W", "P"}
POLAR_UNCHARGED = {"S", "T", "N", "Q", "C", "Y", "G"}
POSITIVELY_CHARGED = {"K", "R", "H"}
NEGATIVELY_CHARGED = {"D", "E"}

# Standard amino acid set for validation
STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")

# Amino acids whose ALL codons contain the GT splice donor dinucleotide
# Valine: GTT, GTC, GTA, GTG — ALL contain GT
GT_MANDATORY_AAS: set[str] = set()

# Pre-compute: which amino acids have GT in ALL their codons?
for _aa, _codons in AA_TO_CODONS.items():
    if all("GT" in codon for codon in _codons):
        GT_MANDATORY_AAS.add(_aa)

# Amino acids whose ALL codons contain the AG splice acceptor dinucleotide
# Check: which amino acids have AG in ALL their codons?
AG_MANDATORY_AAS: set[str] = set()
for _aa, _codons in AA_TO_CODONS.items():
    if all("AG" in codon for codon in _codons):
        AG_MANDATORY_AAS.add(_aa)

# Minimum BLOSUM62 score for a substitution to be considered "conservative"
CONSERVATIVE_THRESHOLD = 0  # BLOSUM62 >= 0

# Minimum BLOSUM62 score for "very conservative" (structure-preserving)
VERY_CONSERVATIVE_THRESHOLD = 2


# ==============================================================================
# Data Structures
# ==============================================================================

@dataclass(frozen=True)
class AASubstitution:
    """A proposed amino acid substitution with rationale.

    Invariants:
    - original_aa and substitute_aa are single standard amino acid codes
    - original_aa != substitute_aa
    - position >= 0
    - blosum62_score == BLOSUM62[original_aa][substitute_aa]
    """
    position: int          # 0-indexed position in protein
    original_aa: str       # original amino acid
    substitute_aa: str     # proposed substitution
    blosum62_score: int    # BLOSUM62 substitution score
    reason: str            # why this substitution is proposed
    predicate_addressed: str  # which predicate this fixes
    cai_impact: float = 0.0   # predicted CAI change (negative = loss)
    confidence: str = "high"   # high, medium, low

    def __post_init__(self):
        """Validate AASubstitution invariants."""
        assert self.position >= 0, f"Position must be non-negative, got {self.position}"
        assert self.original_aa in STANDARD_AAS, (
            f"Original AA must be standard, got '{self.original_aa}'"
        )
        assert self.substitute_aa in STANDARD_AAS, (
            f"Substitute AA must be standard, got '{self.substitute_aa}'"
        )
        assert self.original_aa != self.substitute_aa, (
            f"Self-substitution not allowed: {self.original_aa}->{self.substitute_aa}"
        )
        assert self.confidence in ("high", "medium", "low"), (
            f"Confidence must be high/medium/low, got '{self.confidence}'"
        )

    @property
    def is_conservative(self) -> bool:
        return self.blosum62_score >= CONSERVATIVE_THRESHOLD

    @property
    def is_very_conservative(self) -> bool:
        return self.blosum62_score >= VERY_CONSERVATIVE_THRESHOLD


@dataclass
class MutagenesisResult:
    """Result of type-directed mutagenesis.

    Invariants:
    - len(original_protein) == len(modified_protein)
    - n_substitutions <= len(original_protein)
    - 0.0 <= protein_identity_pct <= 100.0
    """
    original_protein: str
    modified_protein: str
    substitutions: list[AASubstitution]
    iterations: int
    all_predicates_pass: bool
    predicate_improvement: dict[str, str]  # predicate -> "PASS" or "STILL_FAIL"
    cai_before: float = 0.0
    cai_after: float = 0.0
    gc_before: float = 0.0
    gc_after: float = 0.0

    def __post_init__(self):
        """Validate MutagenesisResult invariants."""
        assert len(self.original_protein) == len(self.modified_protein), (
            f"Protein length changed: {len(self.original_protein)} -> "
            f"{len(self.modified_protein)}"
        )
        assert self.iterations >= 0, f"Iterations must be non-negative, got {self.iterations}"

    @property
    def n_substitutions(self) -> int:
        return len(self.substitutions)

    @property
    def protein_identity_pct(self) -> float:
        """Percent identity between original and modified protein."""
        if not self.original_protein:
            return 0.0
        matches = sum(1 for a, b in zip(self.original_protein, self.modified_protein) if a == b)
        return 100.0 * matches / len(self.original_protein)


# ==============================================================================
# Core Analysis Functions
# ==============================================================================

def find_unrepairable_cryptic_donors(
    sequence: str,
    protein: str,
    organism: str,
    threshold: float = 3.0,
) -> list[tuple[int, int, str, float, bool, bool]]:
    """
    Find cryptic splice donor sites that CANNOT be eliminated by codon swaps.

    Returns 6-tuples: (seq_pos, codon_idx, aa, score, fixable, gt_mandatory)
    - gt_mandatory: True if ALL codons for this AA contain GT (truly unrepairable
      at codon level, e.g. Valine). False if GT-free codons exist but optimizer
      didn't use them (optimizer weakness, not a mutagenesis case).

    Pre-conditions:
    - sequence is uppercase DNA, length == len(protein) * 3
    - protein contains only standard amino acid codes
    - threshold > 0

    Post-conditions:
    - all positions in result are valid indices into sequence
    - all scores are >= threshold (only strong donors are reported)
    - gt_mandatory is True iff the amino acid is in GT_MANDATORY_AAS
    """
    assert threshold > 0, f"Threshold must be positive, got {threshold}"
    assert len(sequence) == len(protein) * 3, (
        f"Sequence length ({len(sequence)}) must equal protein length * 3 "
        f"({len(protein) * 3})"
    )

    unrepairable = []
    for i in range(len(sequence) - 1):
        if sequence[i:i+2] == "GT":
            s = score_donor(sequence, i)
            if s >= threshold:
                codon_idx = i // 3
                if codon_idx >= len(protein):
                    continue
                aa = protein[codon_idx]
                codons = AA_TO_CODONS.get(aa, [])

                # Check if ANY alternative codon can bring score below threshold
                fixable = False
                current_codon = sequence[codon_idx*3:codon_idx*3+3]
                for alt_codon in codons:
                    if alt_codon == current_codon:
                        continue
                    test = sequence[:codon_idx*3] + alt_codon + sequence[codon_idx*3+3:]
                    new_score = score_donor(test, i)
                    if new_score < threshold:
                        fixable = True
                        break

                gt_mandatory = is_gt_mandatory(aa)
                unrepairable.append((i, codon_idx, aa, s, fixable, gt_mandatory))

    return unrepairable


def find_unrepairable_cryptic_acceptors(
    sequence: str,
    protein: str,
    organism: str,
    threshold: float = 3.0,
) -> list[tuple[int, int, str, float, bool, bool]]:
    """
    Find cryptic splice acceptor sites that CANNOT be eliminated by codon swaps.

    Returns 6-tuples: (seq_pos, codon_idx, aa, score, fixable, ag_mandatory)
    - ag_mandatory: True if ALL codons for this AA contain AG (truly unrepairable
      at codon level). False if AG-free codons exist but optimizer didn't use
      them (optimizer weakness, not a mutagenesis case).

    Pre-conditions:
    - sequence is uppercase DNA, length == len(protein) * 3
    - protein contains only standard amino acid codes
    - threshold > 0

    Post-conditions:
    - all positions in result are valid indices into sequence
    - all scores are >= threshold (only strong acceptors are reported)
    - ag_mandatory is True iff the amino acid is in AG_MANDATORY_AAS
    """
    assert threshold > 0, f"Threshold must be positive, got {threshold}"
    assert len(sequence) == len(protein) * 3, (
        f"Sequence length ({len(sequence)}) must equal protein length * 3 "
        f"({len(protein) * 3})"
    )

    unrepairable = []
    for i in range(len(sequence) - 1):
        if sequence[i:i+2] == "AG":
            s = score_acceptor(sequence, i)
            if s >= threshold:
                codon_idx = i // 3
                if codon_idx >= len(protein):
                    continue
                aa = protein[codon_idx]
                codons = AA_TO_CODONS.get(aa, [])

                current_codon = sequence[codon_idx*3:codon_idx*3+3]
                fixable = False
                for alt_codon in codons:
                    if alt_codon == current_codon:
                        continue
                    test = sequence[:codon_idx*3] + alt_codon + sequence[codon_idx*3+3:]
                    new_score = score_acceptor(test, i)
                    if new_score < threshold:
                        fixable = True
                        break

                ag_mandatory = is_ag_mandatory(aa)
                unrepairable.append((i, codon_idx, aa, s, fixable, ag_mandatory))

    return unrepairable


def propose_substitutions(
    aa: str,
    predicate_name: str,
    reason: str,
) -> list[tuple[str, int, str]]:
    """
    Propose conservative amino acid substitutions for a given amino acid.

    Pre-conditions:
    - aa is a standard single-letter amino acid code
    - predicate_name is a non-empty string

    Post-conditions:
    - no self-substitutions (aa not in result)
    - all BLOSUM62 scores >= -1
    - results sorted by BLOSUM62 score descending
    """
    assert aa in STANDARD_AAS, f"Invalid amino acid: '{aa}'"
    assert predicate_name, "Predicate name must not be empty"

    candidates = []
    scores = BLOSUM62.get(aa, {})

    for sub_aa, score in scores.items():
        if sub_aa == aa:
            continue
        if score < -1:
            continue  # Skip very unfavorable substitutions

        # Check if substitution actually helps with the predicate
        sub_codons = AA_TO_CODONS.get(sub_aa, [])

        if "CrypticSplice" in predicate_name or "SpliceCorrect" in predicate_name:
            # For splice predicates, the substitution must avoid the problematic
            # dinucleotide. Check if sub_aa has GT-free or AG-free codons.
            if "donor" in reason.lower() or "GT" in reason:
                has_gt_free = any("GT" not in c for c in sub_codons)
                if not has_gt_free and aa in GT_MANDATORY_AAS:
                    confidence = "low"
                elif has_gt_free:
                    confidence = "high"
                else:
                    confidence = "medium"
            elif "acceptor" in reason.lower() or "AG" in reason:
                has_ag_free = any("AG" not in c for c in sub_codons)
                confidence = "high" if has_ag_free else "low"
            else:
                confidence = "medium"
        else:
            # For non-splice predicates, rank by BLOSUM62 only
            if score >= VERY_CONSERVATIVE_THRESHOLD:
                confidence = "high"
            elif score >= CONSERVATIVE_THRESHOLD:
                confidence = "medium"
            else:
                confidence = "low"

        candidates.append((sub_aa, score, confidence))

    # Sort by BLOSUM62 score descending
    candidates.sort(key=lambda x: -x[1])
    return candidates


def apply_substitution(
    protein: str,
    position: int,
    new_aa: str,
) -> str:
    """Apply a single amino acid substitution to a protein sequence.

    Pre-conditions:
    - protein is a valid amino acid sequence
    - new_aa is a standard amino acid code
    - 0 <= position < len(protein) (or position is out of bounds, returns unchanged)

    Post-conditions:
    - len(result) == len(protein)
    - result[position] == new_aa (if position was in bounds)
    - all other positions unchanged
    """
    assert new_aa in STANDARD_AAS, f"Invalid amino acid: '{new_aa}'"
    p = list(protein)
    if position < 0 or position >= len(p):
        logger.warning(
            "Substitution position %d out of bounds [0, %d), protein unchanged",
            position, len(p),
        )
        return protein
    p[position] = new_aa
    result = "".join(p)
    assert len(result) == len(protein), "Substitution changed protein length"
    return result


# ==============================================================================
# Convenience Queries for GT/AG Mandation
# ==============================================================================

def is_gt_mandatory(aa: str) -> bool:
    """Check if ALL codons for an amino acid contain the GT dinucleotide.

    This is the core biological insight that drives mutagenesis: Valine (V)
    is the only standard amino acid where every synonymous codon contains GT,
    making it impossible to eliminate cryptic splice donor sites by codon
    swapping alone.

    Pre-conditions:
    - aa is a single standard amino acid code (uppercase)

    Post-conditions:
    - returns True iff every codon for aa contains "GT"
    - returns False for non-standard or stop codons
    """
    if aa not in AA_TO_CODONS:
        return False
    return all("GT" in codon for codon in AA_TO_CODONS[aa])


def is_ag_mandatory(aa: str) -> bool:
    """Check if ALL codons for an amino acid contain the AG dinucleotide.

    Analogous to is_gt_mandatory but for the AG splice acceptor dinucleotide.
    If an amino acid's codons all contain AG, cryptic splice acceptor sites
    at that position cannot be eliminated by codon swapping.

    Pre-conditions:
    - aa is a single standard amino acid code (uppercase)

    Post-conditions:
    - returns True iff every codon for aa contains "AG"
    - returns False for non-standard or stop codons
    """
    if aa not in AA_TO_CODONS:
        return False
    return all("AG" in codon for codon in AA_TO_CODONS[aa])


def diagnose_optimizer_weakness(
    sequence: str,
    protein: str,
    organism: str = "",
    threshold: float = 3.0,
) -> list[dict]:
    """Find positions where optimizer chose GT-containing codons when GT-free alternatives exist.

    These are NOT mutagenesis candidates — they're optimizer bugs.
    Returns a list of dicts with position, current_codon, gt_free_alternatives, and donor_score.

    Pre-conditions:
    - sequence is uppercase DNA, length == len(protein) * 3
    - protein contains only standard amino acid codes
    - threshold > 0

    Post-conditions:
    - each result dict has keys: seq_pos, codon_idx, aa, current_codon,
      gt_free_alternatives, donor_score
    - gt_free_alternatives is a non-empty list of codons without GT
    - donor_score >= threshold for all results
    - position is a 0-indexed codon position
    """
    assert threshold > 0, f"Threshold must be positive, got {threshold}"
    assert len(sequence) == len(protein) * 3, (
        f"Sequence length ({len(sequence)}) must equal protein length * 3 "
        f"({len(protein) * 3})"
    )

    weak_positions = []
    for i in range(len(sequence) - 1):
        if sequence[i:i+2] == "GT":
            s = score_donor(sequence, i)
            if s >= threshold:
                codon_idx = i // 3
                if codon_idx >= len(protein):
                    continue
                aa = protein[codon_idx]
                current_codon = sequence[codon_idx*3:codon_idx*3+3]
                gt_free = [c for c in AA_TO_CODONS.get(aa, []) if "GT" not in c]

                if gt_free and aa not in GT_MANDATORY_AAS:
                    weak_positions.append({
                        "seq_pos": i,
                        "codon_idx": codon_idx,
                        "aa": aa,
                        "current_codon": current_codon,
                        "gt_free_alternatives": gt_free,
                        "donor_score": s,
                    })

    return weak_positions


def force_gt_free_reoptimization(
    sequence: str,
    protein: str,
    organism: str,
    threshold: float = 3.0,
) -> str:
    """Aggressively force GT-free codons at all non-Valine positions with strong cryptic donors.

    This function is called when NoCrypticSplice fails and there are non-Valine
    positions with strong cryptic donors where GT-free synonymous codons exist.
    Unlike the optimizer's Phase 7 which uses a cautious acceptance check, this
    function UNCONDITIONALLY swaps to the highest-CAI GT-free codon at every
    non-GT-mandatory position with a donor score >= threshold.

    It does NOT propose amino acid substitutions — it only does synonymous
    codon swaps, so the protein is unchanged.

    Pre-conditions:
    - sequence is uppercase DNA, length == len(protein) * 3
    - protein contains only standard amino acid codes
    - organism is a supported organism name
    - threshold > 0

    Post-conditions:
    - returned sequence has the same length as input
    - returned sequence encodes the same protein
    - all non-GT-mandatory positions with donor score >= threshold have been
      swapped to GT-free codons (where available)
    - GT-mandatory positions (Valine) are left unchanged
    """
    assert threshold > 0, f"Threshold must be positive, got {threshold}"
    assert len(sequence) == len(protein) * 3, (
        f"Sequence length ({len(sequence)}) must equal protein length * 3 "
        f"({len(protein) * 3})"
    )

    from .organisms import CODON_ADAPTIVENESS_TABLES
    usage = CODON_ADAPTIVENESS_TABLES.get(organism, {})
    aas = list(protein)

    for iteration in range(5):  # Max 5 rounds for cross-boundary GT propagation
        swapped = False
        # Collect all strong donor positions
        donor_sites = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == "GT":
                s = score_donor(sequence, i)
                if s >= threshold:
                    donor_sites.append((i, s))

        if not donor_sites:
            break

        # Sort by score descending to fix worst sites first
        donor_sites.sort(key=lambda x: x[1], reverse=True)

        for gt_pos, gt_score in donor_sites:
            codon_idx = gt_pos // 3
            if codon_idx >= len(aas):
                continue
            aa = aas[codon_idx]

            # Skip GT-mandatory amino acids (Valine) — can't fix with codon swaps
            if is_gt_mandatory(aa):
                continue

            # Get GT-free alternatives, sorted by CAI
            gt_free = [c for c in AA_TO_CODONS.get(aa, []) if "GT" not in c]
            if not gt_free:
                continue  # Shouldn't happen for non-GT-mandatory, but defensive

            gt_free_sorted = sorted(
                gt_free,
                key=lambda c: usage.get(c, 0.0),
                reverse=True,
            )

            # Unconditionally swap to the best GT-free codon
            current_codon = sequence[codon_idx*3:codon_idx*3+3]
            if current_codon == gt_free_sorted[0]:
                # Already using best GT-free codon — try next best
                if len(gt_free_sorted) > 1:
                    new_codon = gt_free_sorted[1]
                else:
                    continue  # Only one GT-free codon and it's already selected
            else:
                new_codon = gt_free_sorted[0]

            sequence = sequence[:codon_idx*3] + new_codon + sequence[codon_idx*3+3:]
            swapped = True

        if not swapped:
            break

    return sequence


# ==============================================================================
# Main Mutagenesis Loop
# ==============================================================================

def type_directed_mutagenesis(
    protein: str,
    organism: str,
    failed_predicates: list[str],
    optimize_fn: Callable,
    max_iterations: int = 10,
    max_substitutions: int = 30,
    min_blosum62: int = -1,
    cryptic_splice_threshold: float = 3.0,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    **optimize_kwargs,
) -> MutagenesisResult:
    """
    Type-directed protein mutagenesis: when codon optimization fails, propose
    conservative amino acid substitutions to make constraint satisfaction possible.

    This is the core loop:
    1. Optimize the current protein
    2. Check which predicates still fail
    3. For each failing predicate, analyze WHY it fails at the codon level
    4. If failure is due to amino acid identity (e.g., Valine's GT), propose
       conservative substitutions
    5. Apply the best substitution, re-optimize, iterate

    Pre-conditions:
    - protein is a valid amino acid sequence (non-empty, standard codes only)
    - organism is a supported organism name
    - failed_predicates is a non-empty list of predicate names
    - optimize_fn is a callable that accepts (target_protein, organism, ...)
      and returns an OptimizationResult
    - max_iterations > 0
    - max_substitutions > 0
    - 0.0 <= gc_lo < gc_hi <= 1.0

    Post-conditions:
    - len(result.modified_protein) == len(protein)
    - all substitutions are conservative (BLOSUM62 >= min_blosum62)
    - result.iterations <= max_iterations
    """
    # Validate pre-conditions
    assert protein, "Protein must not be empty"
    assert all(aa in STANDARD_AAS for aa in protein), (
        f"Protein contains non-standard amino acids: "
        f"{set(aa for aa in protein if aa not in STANDARD_AAS)}"
    )
    assert failed_predicates, "Must have at least one failed predicate"
    assert max_iterations > 0, f"Max iterations must be positive, got {max_iterations}"
    assert max_substitutions > 0, f"Max substitutions must be positive, got {max_substitutions}"
    assert 0.0 <= gc_lo < gc_hi <= 1.0, f"Invalid GC bounds: [{gc_lo}, {gc_hi}]"
    assert cryptic_splice_threshold > 0, f"Threshold must be positive: {cryptic_splice_threshold}"

    current_protein = protein
    all_substitutions: list[AASubstitution] = []
    predicate_improvement: dict[str, str] = {}

    # Initial optimization
    initial_result = optimize_fn(
        target_protein=current_protein,
        organism=organism,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        cryptic_splice_threshold=cryptic_splice_threshold,
        **optimize_kwargs,
    )
    cai_before = initial_result.cai
    gc_before = initial_result.gc_content

    iterations_used = 0
    for iteration in range(max_iterations):
        iterations_used = iteration + 1

        # Re-optimize current protein
        opt_result = optimize_fn(
            target_protein=current_protein,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cryptic_splice_threshold=cryptic_splice_threshold,
            **optimize_kwargs,
        )
        sequence = opt_result.sequence

        # Check if all predicates pass now
        if not opt_result.failed_predicates:
            break

        if len(all_substitutions) >= max_substitutions:
            logger.warning(
                "Mutagenesis: reached max substitutions (%d), stopping",
                max_substitutions,
            )
            break

        # Analyze failures and find unrepairable positions
        proposed_subs: list[AASubstitution] = []

        for pred in opt_result.failed_predicates:
            # Handle cryptic splice donor sites
            if "CrypticSplice" in pred or "NoCrypticSplice" in pred:
                unrepairable_donors = find_unrepairable_cryptic_donors(
                    sequence, current_protein, organism,
                    threshold=cryptic_splice_threshold,
                )

                for seq_pos, codon_idx, aa, score, fixable, gt_mandatory in unrepairable_donors:
                    if fixable:
                        # The optimizer should handle this — skip
                        continue
                    if not gt_mandatory:
                        # This position HAS GT-free codons but optimizer didn't use them
                        # This is an optimizer bug, not a mutagenesis case — log but skip
                        logger.warning(
                            "Mutagenesis: position %d (AA=%s) has GT-free codons but optimizer "
                            "didn't use them. This is an optimizer issue, not mutagenesis.",
                            codon_idx, aa
                        )
                        continue

                    # Only propose substitutions for GT-mandatory positions (Valine)
                    reason = (
                        f"Cryptic donor at seq pos {seq_pos} (MaxEntScan={score:.2f}): "
                        f"{aa} at codon {codon_idx} has GT in ALL codons (GT-mandatory)"
                    )

                    candidates = propose_substitutions(aa, "NoCrypticSplice", reason)

                    for sub_aa, blosum, confidence in candidates:
                        if blosum < min_blosum62:
                            continue
                        # Check if substitution would eliminate GT at this position
                        sub_codons = AA_TO_CODONS.get(sub_aa, [])
                        has_gt_free = any("GT" not in c for c in sub_codons)

                        proposed_subs.append(AASubstitution(
                            position=codon_idx,
                            original_aa=aa,
                            substitute_aa=sub_aa,
                            blosum62_score=blosum,
                            reason=reason,
                            predicate_addressed="NoCrypticSplice",
                            confidence=confidence if has_gt_free else "low",
                        ))
                        break  # Take the best (highest BLOSUM62) substitution

            # Handle cryptic splice acceptor sites
            if "acceptor" in pred.lower():
                unrepairable_acceptors = find_unrepairable_cryptic_acceptors(
                    sequence, current_protein, organism,
                    threshold=cryptic_splice_threshold,
                )

                for seq_pos, codon_idx, aa, score, fixable, ag_mandatory in unrepairable_acceptors:
                    if fixable:
                        continue
                    if not ag_mandatory:
                        # This position HAS AG-free codons but optimizer didn't use them
                        # This is an optimizer bug, not a mutagenesis case — log but skip
                        logger.warning(
                            "Mutagenesis: position %d (AA=%s) has AG-free codons but optimizer "
                            "didn't use them. This is an optimizer issue, not mutagenesis.",
                            codon_idx, aa
                        )
                        continue

                    # Only propose substitutions for AG-mandatory positions
                    reason = (
                        f"Cryptic acceptor at seq pos {seq_pos} (MaxEntScan={score:.2f}): "
                        f"{aa} at codon {codon_idx} has AG in ALL codons (AG-mandatory)"
                    )

                    candidates = propose_substitutions(aa, "NoCrypticSplice", reason)
                    for sub_aa, blosum, confidence in candidates:
                        if blosum < min_blosum62:
                            continue
                        proposed_subs.append(AASubstitution(
                            position=codon_idx,
                            original_aa=aa,
                            substitute_aa=sub_aa,
                            blosum62_score=blosum,
                            reason=reason,
                            predicate_addressed="NoCrypticSplice",
                            confidence=confidence,
                        ))
                        break

        if not proposed_subs:
            # No more substitutions can be proposed
            logger.info(
                "Mutagenesis: no more substitutions possible at iteration %d",
                iteration,
            )
            break

        # Deduplicate by position, keeping highest BLOSUM62 for each
        best_by_position: dict[int, AASubstitution] = {}
        for sub in proposed_subs:
            if sub.position not in best_by_position or sub.blosum62_score > best_by_position[sub.position].blosum62_score:
                best_by_position[sub.position] = sub

        # Apply substitutions (most conservative first)
        sorted_subs = sorted(best_by_position.values(), key=lambda s: -s.blosum62_score)
        applied_this_round = 0

        for sub in sorted_subs:
            if len(all_substitutions) >= max_substitutions:
                break
            # Only apply high-confidence substitutions first
            if sub.confidence in ("high", "medium"):
                current_protein = apply_substitution(
                    current_protein, sub.position, sub.substitute_aa
                )
                all_substitutions.append(sub)
                applied_this_round += 1
                predicate_improvement[sub.predicate_addressed] = "ADDRESSED"

        if applied_this_round == 0:
            # No progress — stop
            logger.info("Mutagenesis: no progress at iteration %d", iteration)
            break

        logger.info(
            "Mutagenesis iteration %d: applied %d substitutions (%d total)",
            iteration, applied_this_round, len(all_substitutions),
        )

    # Final optimization
    final_result = optimize_fn(
        target_protein=current_protein,
        organism=organism,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        cryptic_splice_threshold=cryptic_splice_threshold,
        **optimize_kwargs,
    )

    # Record final predicate status
    for pred in final_result.failed_predicates:
        if pred not in predicate_improvement:
            predicate_improvement[pred] = "STILL_FAIL"
    for pred in final_result.satisfied_predicates:
        if pred not in predicate_improvement:
            predicate_improvement[pred] = "PASS"

    # Post-condition: protein length must be preserved
    assert len(current_protein) == len(protein), (
        f"Post-condition violation: protein length changed from "
        f"{len(protein)} to {len(current_protein)}"
    )

    return MutagenesisResult(
        original_protein=protein,
        modified_protein=current_protein,
        substitutions=all_substitutions,
        iterations=iterations_used,
        all_predicates_pass=len(final_result.failed_predicates) == 0,
        predicate_improvement=predicate_improvement,
        cai_before=cai_before,
        cai_after=final_result.cai,
        gc_before=gc_before,
        gc_after=final_result.gc_content,
    )
