"""
BioCompiler Deimmunization Engine v7.2.0
==========================================
Reduces protein immunogenicity while preserving function by disrupting
T-cell epitopes through conservative amino acid substitutions.

Algorithm:
  1. Compute current immunogenicity score
  2. If below target_score -> done
  3. Find strongest T-cell epitopes
  4. For each epitope, find positions where mutation reduces MHC binding most
  5. Rank mutations by: binding_reduction * (1 - |blosum62|/4) * (1 - max(0, ddg)/5)
  6. Apply best mutation (if it passes filters: blosum62 >= blosum62_min, ddg < max_ddg)
  7. Repeat until target_score reached or max_mutations applied

References:
  - BLOSUM62: Henikoff & Henikoff (1992) PNAS 89:10915
  - Kyte-Doolittle hydropathy: Kyte & Doolittle (1982) J Mol Biol 157:105
  - MHC binding prediction: NetMHCpan approach (pseudo-sequence + scoring)
  - CamSol: Sormanni et al. (2015) J Mol Biol 427:478
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

from .constants import BLOSUM62, HYDROPATHY

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# Standard amino acids
# ────────────────────────────────────────────────────────────

_STANDARD_AAS = list("ARNDCQEGHILKMFPSTWYV")

# MHC class I peptide length (canonical 9-mer)
_MHC_PEPTIDE_LENGTH = 9

# Default MHC alleles for common organisms
_DEFAULT_MHC_ALLELES: dict[str, list[str]] = {
    "Homo_sapiens": [
        "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01",
        "HLA-B*07:02", "HLA-B*08:01", "HLA-B*27:05",
        "HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01",
    ],
    "Mus_musculus": [
        "H2-Kb", "H2-Db", "H2-IAb",
    ],
    "Escherichia_coli": [],
    "CHO_K1": [
        "MHC-C*01:01", "MHC-C*02:01",
    ],
    "Saccharomyces_cerevisiae": [],
}

# Pseudo-position-specific scoring matrix for MHC class I binding
# Simplified: positions 1, 2, and C-terminal (9) are anchor positions
# Higher positive values = stronger binding preference
_MHC_ANCHOR_POSITIONS = {0, 1, 8}  # 0-based positions in 9-mer
_MHC_ANCHOR_WEIGHT = 2.0
_MHC_NONANCHOR_WEIGHT = 1.0


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class DeimmunizationResult:
    """Result of a deimmunization optimization run."""

    original_protein: str
    optimized_protein: str
    mutations_applied: list[dict]  # [{position, wildtype, mutant, epitope_removed, ddg, blosum62}]
    original_immunogenicity: float  # immunogenicity score before
    optimized_immunogenicity: float  # immunogenicity score after
    original_t_cell_epitopes: int
    optimized_t_cell_epitopes: int
    stability_preserved: bool  # True if sum of all ddG < threshold
    iterations: int
    success: bool  # True if target_score was reached
    method: str  # algorithm name


@dataclass
class EpitopeMutation:
    """A mutation that disrupts a T-cell epitope."""

    position: int  # 0-based position in protein
    wildtype: str  # original amino acid
    mutant: str  # substitution amino acid
    epitope_disrupted: str  # peptide that was an epitope
    binding_reduction: float  # how much binding score decreased
    blosum62: int  # conservation score
    ddg_estimate: float  # estimated stability impact
    solubility_impact: float  # estimated solubility change


# ────────────────────────────────────────────────────────────
# Internal helpers
# ────────────────────────────────────────────────────────────

def _estimate_mhc_binding_score(peptide: str, allele: str = "") -> float:
    """Estimate MHC binding score for a 9-mer peptide.

    Uses a simplified position-specific scoring approach based on
    BLOSUM62 conservation at anchor positions and hydropathy preferences.
    Higher score = stronger predicted binding = more immunogenic.

    This is a heuristic; production use should call out to NetMHCpan.

    Args:
        peptide: 9-mer amino acid string.
        allele: MHC allele name (used for allele-specific adjustments).

    Returns:
        Estimated binding score (arbitrary units, higher = stronger binding).
    """
    if len(peptide) != _MHC_PEPTIDE_LENGTH:
        # For non-9-mers, use a sliding window approach
        if len(peptide) < _MHC_PEPTIDE_LENGTH:
            return 0.0
        # Score all possible 9-mers within the longer peptide
        max_score = 0.0
        for i in range(len(peptide) - _MHC_PEPTIDE_LENGTH + 1):
            score = _estimate_mhc_binding_score(
                peptide[i:i + _MHC_PEPTIDE_LENGTH], allele
            )
            max_score = max(max_score, score)
        return max_score

    score = 0.0
    for pos, aa in enumerate(peptide):
        if aa not in BLOSUM62:
            continue
        # Self-conservation at this position
        self_score = BLOSUM62[aa][aa]

        if pos in _MHC_ANCHOR_POSITIONS:
            # Anchor positions: hydrophobic and large AAs preferred
            hydro = HYDROPATHY.get(aa, 0.0)
            # Strong binding: hydrophobic at anchors
            anchor_bonus = max(0.0, hydro) * 0.5
            # AAs with high self-conservation at anchors are preferred
            score += (self_score / 4.0 + anchor_bonus) * _MHC_ANCHOR_WEIGHT
        else:
            # Non-anchor: more permissive, still use conservation
            score += (self_score / 6.0) * _MHC_NONANCHOR_WEIGHT

    # Allele-specific adjustments
    if "A*02" in allele:
        # HLA-A*02:01 prefers Leu/Met/Val at position 2, Val/Leu at C-term
        if peptide[1] in "LMV":
            score += 1.5
        if peptide[-1] in "VLIA":
            score += 1.0
    elif "A*01" in allele:
        # HLA-A*01 prefers Tyr/Phe at position 2
        if peptide[1] in "YFW":
            score += 1.2
    elif "B*07" in allele:
        # HLA-B*07 prefers Pro at position 2
        if peptide[1] == "P":
            score += 1.5
    elif "B*27" in allele:
        # HLA-B*27 prefers Arg at position 2
        if peptide[1] == "R":
            score += 2.0
    elif "DRB1" in allele:
        # MHC class II: different anchor pattern
        # Hydrophobic at position 1, small at position 6
        if len(peptide) >= 6 and peptide[0] in "AVLIMFYW":
            score += 1.0
        if len(peptide) >= 6 and peptide[5] in "GASPN":
            score += 0.8

    return max(0.0, score)


def _compute_binding_score_for_region(
    protein: str, start: int, end: int, mhc_alleles: list[str] | None = None
) -> float:
    """Compute the maximum MHC binding score for a protein region.

    Scores all possible 9-mers overlapping the [start, end) region
    against all specified alleles.

    Args:
        protein: Full protein sequence.
        start: Start position (0-based) of region.
        end: End position (0-based, exclusive) of region.
        mhc_alleles: MHC alleles to test.

    Returns:
        Maximum binding score across all 9-mers and alleles.
    """
    if mhc_alleles is None:
        mhc_alleles = _DEFAULT_MHC_ALLELES.get("Homo_sapiens", [])

    max_score = 0.0
    # Scan all 9-mers that overlap the region
    scan_start = max(0, start - _MHC_PEPTIDE_LENGTH + 1)
    scan_end = min(len(protein) - _MHC_PEPTIDE_LENGTH + 1, end)

    for i in range(scan_start, scan_end):
        peptide = protein[i:i + _MHC_PEPTIDE_LENGTH]
        if len(peptide) < _MHC_PEPTIDE_LENGTH:
            continue
        # Check if this 9-mer overlaps the target region
        pep_end = i + _MHC_PEPTIDE_LENGTH
        if pep_end <= start or i >= end:
            continue
        for allele in mhc_alleles:
            score = _estimate_mhc_binding_score(peptide, allele)
            max_score = max(max_score, score)

    return max_score


def _estimate_ddg(wildtype: str, mutant: str) -> float:
    """Estimate ΔΔG from BLOSUM62 score and hydropathy change.

    Uses the empirical correlation:
      - Lower BLOSUM62 score -> higher ΔΔG (less conservative = more destabilizing)
      - Large hydropathy change -> higher ΔΔG

    Returns:
        Estimated ΔΔG in kcal/mol (positive = destabilizing).
    """
    blosum = BLOSUM62.get(wildtype, {}).get(mutant, -4)

    # BLOSUM62-based estimate: map [-4, 4] -> [0, ~5] kcal/mol
    # Higher blosum = more conservative = less destabilizing
    blosum_component = max(0.0, (-blosum + 4) * 0.3)

    # Hydropathy change: large hydro changes destabilize
    hydro_wt = HYDROPATHY.get(wildtype, 0.0)
    hydro_mt = HYDROPATHY.get(mutant, 0.0)
    hydro_change = abs(hydro_wt - hydro_mt)
    hydro_component = hydro_change * 0.2

    # Proline substitutions in non-Pro positions are especially destabilizing
    proline_penalty = 0.0
    if mutant == "P" and wildtype != "P":
        proline_penalty = 1.0
    elif wildtype == "P" and mutant != "P":
        proline_penalty = 0.5

    ddg = blosum_component + hydro_component + proline_penalty
    return round(ddg, 3)


def _estimate_solubility_impact(wildtype: str, mutant: str) -> float:
    """Estimate solubility impact of a substitution using CamSol-like heuristics.

    Uses intrinsic solubility propensity based on hydropathy and charge.
    Negative values = decreased solubility, positive = increased.

    Returns:
        Estimated solubility change (arbitrary units).
    """
    # More negative hydropathy = more hydrophilic = better solubility
    hydro_wt = HYDROPATHY.get(wildtype, 0.0)
    hydro_mt = HYDROPATHY.get(mutant, 0.0)

    # Change in hydrophilicity: more negative hydro = more soluble
    solubility_change = (hydro_wt - hydro_mt) * 0.15

    # Charged residues promote solubility
    charged_positive = {"K", "R"}
    charged_negative = {"D", "E"}
    charged = charged_positive | charged_negative

    if mutant in charged and wildtype not in charged:
        solubility_change += 0.3
    elif wildtype in charged and mutant not in charged:
        solubility_change -= 0.3

    # Disulfide-breaking mutations (C -> anything) decrease stability
    if wildtype == "C" and mutant != "C":
        solubility_change -= 0.2

    return round(solubility_change, 3)


def _get_mhc_alleles(organism: str) -> list[str]:
    """Get default MHC alleles for an organism."""
    return _DEFAULT_MHC_ALLELES.get(organism, [])


def _organism_to_species(organism: str) -> str:
    """Convert organism name to species name used by the immunogenicity module.

    The immunogenicity module uses 'human'/'mouse' style names,
    while deimmunize() uses 'Homo_sapiens' style names.
    """
    mapping = {
        "Homo_sapiens": "human",
        "Mus_musculus": "mouse",
        "Escherichia_coli": "ecoli",
        "CHO_K1": "cho",
        "Saccharomyces_cerevisiae": "yeast",
    }
    return mapping.get(organism, organism.lower().split("_")[0] if "_" in organism else organism.lower())


# ────────────────────────────────────────────────────────────
# Public API
# ────────────────────────────────────────────────────────────

def compute_mutation_impact(
    protein: str, position: int, mutant_aa: str
) -> dict:
    """Compute the impact of a single amino acid mutation.

    Assesses the effect on MHC binding, stability, and solubility.

    Args:
        protein: Protein sequence (1-letter codes).
        position: 0-based position to mutate.
        mutant_aa: Substitution amino acid.

    Returns:
        Dictionary with keys:
          - binding_impact: list of dicts describing affected epitopes
          - stability_impact: estimated ΔΔG
          - solubility_impact: estimated solubility change
          - blosum62: BLOSUM62 score for the substitution
    """
    if not protein or position < 0 or position >= len(protein):
        return {
            "binding_impact": [],
            "stability_impact": 0.0,
            "solubility_impact": 0.0,
            "blosum62": 0,
        }

    wildtype = protein[position]
    if wildtype not in BLOSUM62 or mutant_aa not in BLOSUM62:
        return {
            "binding_impact": [],
            "stability_impact": 0.0,
            "solubility_impact": 0.0,
            "blosum62": 0,
        }

    # Compute BLOSUM62 score
    blosum = BLOSUM62[wildtype][mutant_aa]

    # Estimate ΔΔG
    ddg = _estimate_ddg(wildtype, mutant_aa)

    # Estimate solubility impact
    sol_impact = _estimate_solubility_impact(wildtype, mutant_aa)

    # Compute binding impact: find all 9-mer epitopes overlapping this position
    binding_impact = []
    for i in range(max(0, position - _MHC_PEPTIDE_LENGTH + 1),
                   min(len(protein) - _MHC_PEPTIDE_LENGTH + 1, position + 1)):
        original_peptide = protein[i:i + _MHC_PEPTIDE_LENGTH]
        if len(original_peptide) < _MHC_PEPTIDE_LENGTH:
            continue

        # Only consider 9-mers that include the mutation position
        if not (i <= position < i + _MHC_PEPTIDE_LENGTH):
            continue

        # Create mutated peptide
        mut_pos_in_peptide = position - i
        mutated_peptide = (
            original_peptide[:mut_pos_in_peptide]
            + mutant_aa
            + original_peptide[mut_pos_in_peptide + 1:]
        )

        # Score before and after
        orig_score = _estimate_mhc_binding_score(original_peptide)
        mut_score = _estimate_mhc_binding_score(mutated_peptide)

        if orig_score > 0:
            binding_impact.append({
                "epitope": original_peptide,
                "start": i,
                "end": i + _MHC_PEPTIDE_LENGTH,
                "original_binding": round(orig_score, 3),
                "mutated_binding": round(mut_score, 3),
                "binding_reduction": round(orig_score - mut_score, 3),
            })

    return {
        "binding_impact": binding_impact,
        "stability_impact": ddg,
        "solubility_impact": sol_impact,
        "blosum62": blosum,
    }


def find_epitope_disrupting_mutations(
    protein: str,
    epitope_start: int,
    epitope_end: int,
    mhc_alleles: list[str] | None = None,
    blosum62_min: int = 0,
) -> list[EpitopeMutation]:
    """Find mutations that disrupt a specific T-cell epitope region.

    Tries all 19 possible substitutions at each position within the
    epitope region, computes binding score change, and filters by
    BLOSUM62 conservation.

    Args:
        protein: Full protein sequence.
        epitope_start: Start position (0-based) of epitope.
        epitope_end: End position (0-based, exclusive) of epitope.
        mhc_alleles: MHC alleles to test binding against.
        blosum62_min: Minimum BLOSUM62 score for acceptable substitutions.

    Returns:
        List of EpitopeMutation objects, sorted by binding_reduction
        (most disruption first).
    """
    if mhc_alleles is None:
        mhc_alleles = _DEFAULT_MHC_ALLELES.get("Homo_sapiens", [])

    # Compute original binding score for the region
    original_binding = _compute_binding_score_for_region(
        protein, epitope_start, epitope_end, mhc_alleles
    )

    mutations: list[EpitopeMutation] = []

    for pos in range(max(0, epitope_start), min(len(protein), epitope_end)):
        wildtype = protein[pos]
        if wildtype not in BLOSUM62:
            continue

        for mutant in _STANDARD_AAS:
            if mutant == wildtype:
                continue
            if mutant not in BLOSUM62:
                continue

            blosum = BLOSUM62[wildtype][mutant]
            if blosum < blosum62_min:
                continue

            # Build mutated protein
            mutated_protein = protein[:pos] + mutant + protein[pos + 1:]

            # Compute new binding score for the region
            new_binding = _compute_binding_score_for_region(
                mutated_protein, epitope_start, epitope_end, mhc_alleles
            )

            binding_reduction = original_binding - new_binding

            # Only include mutations that reduce binding
            if binding_reduction <= 0:
                continue

            # Estimate ΔΔG
            ddg = _estimate_ddg(wildtype, mutant)

            # Estimate solubility impact
            sol_impact = _estimate_solubility_impact(wildtype, mutant)

            # Identify the epitope peptide (the 9-mer with highest binding)
            epitope_peptide = protein[
                max(0, pos - _MHC_PEPTIDE_LENGTH + 1):
                pos + _MHC_PEPTIDE_LENGTH
            ]
            # Use the specific 9-mer that includes this position
            best_pep_start = max(0, pos - _MHC_PEPTIDE_LENGTH + 1)
            # Find the 9-mer with highest binding score
            best_score = 0.0
            best_peptide = ""
            for pep_start in range(
                max(0, pos - _MHC_PEPTIDE_LENGTH + 1),
                min(len(protein) - _MHC_PEPTIDE_LENGTH + 1, pos + 1)
            ):
                if not (pep_start <= pos < pep_start + _MHC_PEPTIDE_LENGTH):
                    continue
                pep = protein[pep_start:pep_start + _MHC_PEPTIDE_LENGTH]
                if len(pep) < _MHC_PEPTIDE_LENGTH:
                    continue
                sc = max(
                    _estimate_mhc_binding_score(pep, allele)
                    for allele in (mhc_alleles or [""])
                )
                if sc > best_score:
                    best_score = sc
                    best_peptide = pep

            mutations.append(EpitopeMutation(
                position=pos,
                wildtype=wildtype,
                mutant=mutant,
                epitope_disrupted=best_peptide,
                binding_reduction=round(binding_reduction, 3),
                blosum62=blosum,
                ddg_estimate=ddg,
                solubility_impact=sol_impact,
            ))

    # Sort by binding_reduction descending (most disruption first)
    mutations.sort(key=lambda m: m.binding_reduction, reverse=True)
    return mutations


def rank_deimmunization_mutations(
    protein: str,
    mhc_alleles: list[str] | None = None,
    blosum62_min: int = 0,
) -> list[EpitopeMutation]:
    """Find all possible deimmunization mutations across all epitopes.

    Identifies T-cell epitopes in the protein, then for each epitope
    finds mutations that reduce MHC binding. Results are ranked by a
    combined score considering binding reduction, conservation, and stability.

    Args:
        protein: Protein sequence.
        mhc_alleles: MHC alleles to test. Defaults to Homo_sapiens.
        blosum62_min: Minimum BLOSUM62 score for acceptable substitutions.

    Returns:
        Ranked list of EpitopeMutation objects (best first).
    """
    # Lazy import to avoid circular dependency
    try:
        from biocompiler.immunogenicity import compute_immunogenicity
        # Use immunogenicity module's T-cell epitope data
        species = _organism_to_species("Homo_sapiens")
        result = compute_immunogenicity(protein, mhc_alleles=mhc_alleles, species=species)
        if result.num_t_cell_epitopes == 0:
            return []
    except ImportError:
        pass
    except Exception:
        pass

    # Always use internal epitope scanning for detailed position data
    epitopes = _scan_t_cell_epitopes(protein, mhc_alleles)

    if not epitopes:
        return []

    all_mutations: list[EpitopeMutation] = []

    for epitope in epitopes:
        start = epitope.get("start", epitope.get("position", 0))
        end = epitope.get("end", start + _MHC_PEPTIDE_LENGTH)
        score = epitope.get("score", epitope.get("binding_score", 0.0))

        # Only consider epitopes with meaningful binding
        if score <= 0:
            continue

        epi_mutations = find_epitope_disrupting_mutations(
            protein, start, end, mhc_alleles, blosum62_min
        )
        all_mutations.extend(epi_mutations)

    # Rank by combined score:
    # binding_reduction * (1 - |blosum62|/4) * (1 - max(0, ddg)/5)
    def _combined_score(m: EpitopeMutation) -> float:
        conservation_factor = 1.0 - abs(m.blosum62) / 4.0
        stability_factor = 1.0 - max(0.0, m.ddg_estimate) / 5.0
        # Ensure factors are non-negative
        conservation_factor = max(0.0, conservation_factor)
        stability_factor = max(0.0, stability_factor)
        return m.binding_reduction * conservation_factor * stability_factor

    all_mutations.sort(key=_combined_score, reverse=True)

    # Deduplicate: keep the best mutation per position
    seen_positions: dict[int, EpitopeMutation] = {}
    for mut in all_mutations:
        key = (mut.position, mut.mutant)
        if key not in seen_positions:
            seen_positions[key] = mut

    deduped = sorted(seen_positions.values(), key=_combined_score, reverse=True)
    return deduped


def _scan_t_cell_epitopes(
    protein: str,
    mhc_alleles: list[str] | None = None,
) -> list[dict]:
    """Scan protein for T-cell epitopes using internal heuristic scoring.

    This is a fallback when the immunogenicity module is not available.
    Identifies 9-mer peptides with high predicted MHC binding.

    Args:
        protein: Protein sequence.
        mhc_alleles: MHC alleles to test.

    Returns:
        List of dicts with keys: start, end, peptide, score.
    """
    if mhc_alleles is None:
        mhc_alleles = _DEFAULT_MHC_ALLELES.get("Homo_sapiens", [])

    epitopes = []

    for i in range(len(protein) - _MHC_PEPTIDE_LENGTH + 1):
        peptide = protein[i:i + _MHC_PEPTIDE_LENGTH]

        # Score against each allele
        max_score = 0.0
        best_allele = ""
        for allele in mhc_alleles:
            score = _estimate_mhc_binding_score(peptide, allele)
            if score > max_score:
                max_score = score
                best_allele = allele

        # Threshold: scores above ~3.0 are potential epitopes
        if max_score > 3.0:
            epitopes.append({
                "start": i,
                "end": i + _MHC_PEPTIDE_LENGTH,
                "peptide": peptide,
                "score": round(max_score, 3),
                "allele": best_allele,
            })

    # Sort by score descending
    epitopes.sort(key=lambda e: e["score"], reverse=True)
    return epitopes


def _compute_immunogenicity_score(
    protein: str,
    organism: str = "Homo_sapiens",
    mhc_alleles: list[str] | None = None,
) -> float:
    """Compute overall immunogenicity score for a protein.

    The score is the average of the top epitope binding scores,
    normalized to [0, 1] range. Higher = more immunogenic.

    Args:
        protein: Protein sequence.
        organism: Target organism.
        mhc_alleles: MHC alleles to test.

    Returns:
        Immunogenicity score in [0, 1].
    """
    # Try the immunogenicity module first
    try:
        from biocompiler.immunogenicity import compute_immunogenicity
        species = _organism_to_species(organism)
        result = compute_immunogenicity(protein, mhc_alleles=mhc_alleles, species=species)
        return result.immunogenicity_score
    except ImportError:
        pass
    except Exception:
        pass  # Fall back to internal heuristic

    if mhc_alleles is None:
        mhc_alleles = _get_mhc_alleles(organism)

    epitopes = _scan_t_cell_epitopes(protein, mhc_alleles)

    if not epitopes:
        return 0.0

    # Take top 5 epitope scores (or fewer if < 5)
    top_scores = [e["score"] for e in epitopes[:5]]

    # Normalize: typical strong binder scores are ~5-10
    # Map to [0, 1] using sigmoid-like transformation
    avg_score = sum(top_scores) / len(top_scores)
    normalized = 1.0 / (1.0 + math.exp(-0.5 * (avg_score - 5.0)))

    return round(normalized, 4)


def _count_t_cell_epitopes(
    protein: str,
    organism: str = "Homo_sapiens",
    mhc_alleles: list[str] | None = None,
) -> int:
    """Count the number of predicted T-cell epitopes in a protein.

    Args:
        protein: Protein sequence.
        organism: Target organism.
        mhc_alleles: MHC alleles to test.

    Returns:
        Number of predicted epitopes.
    """
    try:
        from biocompiler.immunogenicity import compute_immunogenicity
        species = _organism_to_species(organism)
        result = compute_immunogenicity(protein, mhc_alleles=mhc_alleles, species=species)
        return result.num_t_cell_epitopes
    except ImportError:
        pass
    except Exception:
        pass  # Fall back to internal heuristic

    if mhc_alleles is None:
        mhc_alleles = _get_mhc_alleles(organism)

    epitopes = _scan_t_cell_epitopes(protein, mhc_alleles)
    return len(epitopes)


def _compute_solubility_score(protein: str) -> float:
    """Compute intrinsic solubility score using CamSol heuristic.

    Falls back to the camsol module if available.

    Args:
        protein: Protein sequence.

    Returns:
        Solubility score (arbitrary units, higher = more soluble).
    """
    try:
        from biocompiler.camsol import compute_intrinsic_solubility
        result = compute_intrinsic_solubility(protein)
        # Handle SolubilityResult object (camsol module returns objects)
        if hasattr(result, 'solubility_score'):
            return result.solubility_score
        return float(result)
    except ImportError:
        pass
    except Exception:
        pass  # Fall back to internal heuristic

    # Simplified CamSol-like intrinsic solubility
    # Based on average hydrophilicity and charge distribution
    if not protein:
        return 0.0

    score = 0.0
    for aa in protein:
        hydro = HYDROPATHY.get(aa, 0.0)
        # More negative hydropathy = more hydrophilic = better solubility
        score -= hydro * 0.1
        # Charged residues help solubility
        if aa in "DEKR":
            score += 0.2
        # Aggregation-prone residues
        if aa in "VFILMYW":
            score -= 0.15

    # Normalize by length
    score /= len(protein)
    return round(score, 4)


def deimmunize(
    protein: str,
    organism: str = "Homo_sapiens",
    target_score: float = 0.3,
    max_mutations: int = 10,
    blosum62_min: int = 0,
    max_ddg: float = 2.0,
    preserve_positions: list[int] | None = None,
) -> DeimmunizationResult:
    """Iteratively reduce protein immunogenicity by disrupting T-cell epitopes.

    Algorithm:
      1. Compute current immunogenicity score
      2. If below target_score -> done
      3. Find strongest T-cell epitopes
      4. For each epitope, find positions where mutation reduces MHC binding most
      5. Rank mutations by: binding_reduction * (1 - |blosum62|/4) * (1 - max(0, ddg)/5)
      6. Apply best mutation (if it passes filters: blosum62 >= blosum62_min, ddg < max_ddg)
      7. Repeat until target_score reached or max_mutations applied

    Args:
        protein: Protein sequence (1-letter amino acid codes).
        organism: Target organism for MHC allele selection.
        target_score: Target immunogenicity score (lower = less immunogenic).
        max_mutations: Maximum number of mutations to apply.
        blosum62_min: Minimum BLOSUM62 conservation score for allowed substitutions.
        max_ddg: Maximum allowed ΔΔG per mutation (kcal/mol).
        preserve_positions: List of 0-based positions that must not be mutated
            (e.g., active site residues).

    Returns:
        DeimmunizationResult with optimized protein and detailed metrics.
    """
    preserve_set = set(preserve_positions) if preserve_positions else set()
    mhc_alleles = _get_mhc_alleles(organism)

    # Compute initial metrics
    original_immunogenicity = _compute_immunogenicity_score(
        protein, organism, mhc_alleles
    )
    original_epitope_count = _count_t_cell_epitopes(protein, organism, mhc_alleles)

    current_protein = protein
    mutations_applied: list[dict] = []
    total_ddg = 0.0
    iteration = 0

    while iteration < max_mutations:
        # Check if target is reached
        current_score = _compute_immunogenicity_score(
            current_protein, organism, mhc_alleles
        )
        if current_score <= target_score:
            break

        # Find T-cell epitopes
        epitopes = _scan_t_cell_epitopes(current_protein, mhc_alleles)
        if not epitopes:
            break

        # Collect all possible mutations for the strongest epitopes
        # Process epitopes in order of strength (strongest first)
        all_candidates: list[EpitopeMutation] = []

        for epitope in epitopes:
            start = epitope["start"]
            end = epitope["end"]

            epi_mutations = find_epitope_disrupting_mutations(
                current_protein, start, end, mhc_alleles, blosum62_min
            )
            all_candidates.extend(epi_mutations)

        if not all_candidates:
            logger.info("No disruptive mutations found for remaining epitopes")
            break

        # Rank by combined score
        def _combined_score(m: EpitopeMutation) -> float:
            conservation_factor = 1.0 - abs(m.blosum62) / 4.0
            stability_factor = 1.0 - max(0.0, m.ddg_estimate) / 5.0
            conservation_factor = max(0.0, conservation_factor)
            stability_factor = max(0.0, stability_factor)
            return m.binding_reduction * conservation_factor * stability_factor

        all_candidates.sort(key=_combined_score, reverse=True)

        # Find the best mutation that passes all filters
        best_mutation = None
        for candidate in all_candidates:
            # Skip preserved positions
            if candidate.position in preserve_set:
                continue

            # Check BLOSUM62 filter
            if candidate.blosum62 < blosum62_min:
                continue

            # Check ΔΔG filter
            if candidate.ddg_estimate >= max_ddg:
                continue

            # Check cumulative ΔΔG
            if total_ddg + candidate.ddg_estimate >= max_ddg * 2:
                continue

            # Check if this position was already mutated
            already_mutated = any(
                m["position"] == candidate.position for m in mutations_applied
            )
            if already_mutated:
                continue

            best_mutation = candidate
            break

        if best_mutation is None:
            logger.info("No mutation passes all filters")
            break

        # Apply the mutation
        current_protein = (
            current_protein[:best_mutation.position]
            + best_mutation.mutant
            + current_protein[best_mutation.position + 1:]
        )

        total_ddg += best_mutation.ddg_estimate

        mutations_applied.append({
            "position": best_mutation.position,
            "wildtype": best_mutation.wildtype,
            "mutant": best_mutation.mutant,
            "epitope_removed": best_mutation.epitope_disrupted,
            "ddg": best_mutation.ddg_estimate,
            "blosum62": best_mutation.blosum62,
            "binding_reduction": best_mutation.binding_reduction,
            "solubility_impact": best_mutation.solubility_impact,
        })

        iteration += 1
        logger.debug(
            "Iteration %d: %s%d%s (ddG=%.2f, BLOSUM62=%d, binding_reduction=%.3f)",
            iteration,
            best_mutation.wildtype,
            best_mutation.position + 1,  # 1-based for display
            best_mutation.mutant,
            best_mutation.ddg_estimate,
            best_mutation.blosum62,
            best_mutation.binding_reduction,
        )

    # Compute final metrics
    optimized_immunogenicity = _compute_immunogenicity_score(
        current_protein, organism, mhc_alleles
    )
    optimized_epitope_count = _count_t_cell_epitopes(
        current_protein, organism, mhc_alleles
    )

    # Check stability: all ΔΔG values summed < max_ddg * len(mutations)
    stability_threshold = max_ddg * max(len(mutations_applied), 1)
    stability_preserved = total_ddg < stability_threshold

    success = optimized_immunogenicity <= target_score

    return DeimmunizationResult(
        original_protein=protein,
        optimized_protein=current_protein,
        mutations_applied=mutations_applied,
        original_immunogenicity=original_immunogenicity,
        optimized_immunogenicity=optimized_immunogenicity,
        original_t_cell_epitopes=original_epitope_count,
        optimized_t_cell_epitopes=optimized_epitope_count,
        stability_preserved=stability_preserved,
        iterations=iteration,
        success=success,
        method="iterative_epitope_disruption",
    )


def validate_deimmunized_protein(
    protein: str,
    original_protein: str,
    organism: str = "Homo_sapiens",
) -> dict:
    """Validate that a deimmunized protein meets quality criteria.

    Checks:
      a. Reduced immunogenicity
      b. Preserved stability (cumulative ΔΔG < threshold)
      c. Preserved solubility (CamSol score not significantly decreased)
      d. All mutations are conservative (BLOSUM62 check)

    Args:
        protein: Deimmunized protein sequence.
        original_protein: Original (pre-deimmunization) protein sequence.
        organism: Target organism.

    Returns:
        Validation report dict with keys:
          - immunogenicity_reduced: bool
          - original_immunogenicity: float
          - optimized_immunogenicity: float
          - stability_preserved: bool
          - total_ddg: float
          - solubility_preserved: bool
          - original_solubility: float
          - optimized_solubility: float
          - solubility_change: float
          - all_mutations_conservative: bool
          - mutations: list of dicts with validation per mutation
          - overall_valid: bool
    """
    mhc_alleles = _get_mhc_alleles(organism)

    # Check immunogenicity
    original_immuno = _compute_immunogenicity_score(
        original_protein, organism, mhc_alleles
    )
    optimized_immuno = _compute_immunogenicity_score(
        protein, organism, mhc_alleles
    )
    immunogenicity_reduced = optimized_immuno < original_immuno

    # Check solubility
    original_solubility = _compute_solubility_score(original_protein)
    optimized_solubility = _compute_solubility_score(protein)
    solubility_change = optimized_solubility - original_solubility
    # Allow up to 1.0 decrease in solubility score
    solubility_preserved = solubility_change > -1.0

    # Analyze mutations
    mutations = []
    total_ddg = 0.0
    all_conservative = True
    min_blosum_threshold = -2  # Mutations below this are non-conservative

    for i, (wt_aa, mt_aa) in enumerate(zip(original_protein, protein)):
        if wt_aa == mt_aa:
            continue

        blosum = BLOSUM62.get(wt_aa, {}).get(mt_aa, -4)
        ddg = _estimate_ddg(wt_aa, mt_aa)
        total_ddg += ddg

        is_conservative = blosum >= min_blosum_threshold
        if not is_conservative:
            all_conservative = False

        mutations.append({
            "position": i,
            "wildtype": wt_aa,
            "mutant": mt_aa,
            "blosum62": blosum,
            "ddg_estimate": ddg,
            "conservative": is_conservative,
        })

    # Stability: cumulative ΔΔG should be < 5.0 kcal/mol
    stability_preserved = total_ddg < 5.0

    # Overall validation
    overall_valid = (
        immunogenicity_reduced
        and stability_preserved
        and solubility_preserved
        and all_conservative
    )

    return {
        "immunogenicity_reduced": immunogenicity_reduced,
        "original_immunogenicity": original_immuno,
        "optimized_immunogenicity": optimized_immuno,
        "stability_preserved": stability_preserved,
        "total_ddg": round(total_ddg, 3),
        "solubility_preserved": solubility_preserved,
        "original_solubility": original_solubility,
        "optimized_solubility": optimized_solubility,
        "solubility_change": round(solubility_change, 4),
        "all_mutations_conservative": all_conservative,
        "mutations": mutations,
        "overall_valid": overall_valid,
    }
