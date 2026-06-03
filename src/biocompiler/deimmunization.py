"""
BioCompiler Deimmunization Engine v7.5.0
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
  - MHC binding prediction: PSSM-based scoring (immunogenicity module)
  - CamSol: Sormanni et al. (2015) J Mol Biol 427:478
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

from .constants import BLOSUM62, HYDROPATHY, STANDARD_AAS
from .engine_base import (
    BaseEngineResult,
    EngineTimer,
    MutationResult,
    validate_protein_sequence,
)
from .exceptions import ImmunogenicityError
from .immunogenicity import (
    predict_mhc_i_binding,
    predict_mhc_ii_binding,
    predict_t_cell_epitopes,
    score_peptide_pssm,
)

logger = logging.getLogger(__name__)


__all__ = [
    "DeimmunizationResult",
    "EpitopeMutation",
    "compute_mutation_impact",
    "deimmunize",
    "find_epitope_disrupting_mutations",
    "rank_deimmunization_mutations",
    "validate_deimmunized_protein",
]


# ────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────

# MHC class I peptide length (canonical 9-mer)
_MHC_I_PEPTIDE_LENGTH = 9

# Default MHC alleles for common organisms.
# Note: alleles not present in the immunogenicity module's PSSMs will be
# silently skipped by score_peptide_pssm and predict_t_cell_epitopes.
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


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class DeimmunizationResult(BaseEngineResult):
    """Result of a deimmunization optimization run.

    Inherits from BaseEngineResult for unified API compatibility.
    Domain-specific fields are preserved for backward compatibility.

    Unified field mapping:
      sequence → optimized_protein
      primary_score → optimized_immunogenicity
      classification → 'deimmunized' | 'partially_deimmunized' | 'failed'
      engine_name → 'deimmunization'
      primary_score_label → 'immunogenicity_score'
    """

    # Override base class required fields with defaults for keyword-arg compat
    sequence: str = ""
    primary_score: float = 0.0
    classification: str = ""
    success: bool = False
    error: str | None = None
    execution_time_s: float = 0.0
    engine_name: str = "deimmunization"
    primary_score_label: str = "immunogenicity_score"

    # Domain-specific fields (backward-compatible)
    original_protein: str = ""
    optimized_protein: str = ""
    mutations_applied: list[dict] = field(default_factory=list)  # [{position, wildtype, mutant, epitope_removed, ddg, blosum62}]
    original_immunogenicity: float = 0.0  # immunogenicity score before
    optimized_immunogenicity: float = 0.0  # immunogenicity score after
    original_t_cell_epitopes: int = 0
    optimized_t_cell_epitopes: int = 0
    stability_preserved: bool = False  # True if sum of all ddG < threshold
    iterations: int = 0
    method: str = "iterative_epitope_disruption"  # algorithm name

    def __post_init__(self):
        # Sync unified base fields from domain-specific fields
        if not self.sequence and self.optimized_protein:
            object.__setattr__(self, 'sequence', self.optimized_protein)
        if self.primary_score == 0.0 and self.optimized_immunogenicity != 0.0:
            object.__setattr__(self, 'primary_score', self.optimized_immunogenicity)
        elif self.optimized_immunogenicity == 0.0 and self.primary_score != 0.0:
            object.__setattr__(self, 'optimized_immunogenicity', self.primary_score)
        if not self.classification:
            if self.success:
                label = "deimmunized"
            elif self.optimized_immunogenicity < self.original_immunogenicity:
                label = "partially_deimmunized"
            else:
                label = "failed"
            object.__setattr__(self, 'classification', label)

    @property
    def immunogenicity_score(self) -> float:
        """Unified API alias for optimized_immunogenicity."""
        return self.optimized_immunogenicity

    @property
    def mutations(self) -> list[dict]:
        """Unified API alias for mutations_applied."""
        return self.mutations_applied


@dataclass
class EpitopeMutation(MutationResult):
    """Deimmunization-specific mutation result, compatible with MutationResult.

    Subclass of MutationResult with deimmunization-appropriate defaults.
    Use to_mutation_result() to convert to a plain MutationResult.
    """
    score_type: str = "immunogenicity"
    engine: str = "deimmunization"
    recommendation: str = "deimmunizing"

    def to_mutation_result(self) -> MutationResult:
        """Convert to a plain MutationResult."""
        return MutationResult(
            position=self.position,
            original=self.original,
            mutant=self.mutant,
            delta_score=self.delta_score,
            score_type=self.score_type,
            engine=self.engine,
            recommendation=self.recommendation,
            description=self.description,
            details=self.details,
        )


# ────────────────────────────────────────────────────────────
# Internal helpers
# ────────────────────────────────────────────────────────────

def _compute_binding_score_for_region(
    protein: str, start: int, end: int, mhc_alleles: list[str] | None = None
) -> float:
    """Compute the maximum MHC binding score for a protein region.

    Scores all possible 9-mers overlapping the [start, end) region
    against all specified alleles using score_peptide_pssm from the
    immunogenicity module.

    Args:
        protein: Full protein sequence.
        start: Start position (0-based) of region.
        end: End position (0-based, exclusive) of region.
        mhc_alleles: MHC alleles to test.

    Returns:
        Maximum binding score across all 9-mers and alleles (0.0–1.0).
    """
    if mhc_alleles is None:
        mhc_alleles = _DEFAULT_MHC_ALLELES.get("Homo_sapiens", [])

    max_score = 0.0
    # Scan all 9-mers that overlap the region
    scan_start = max(0, start - _MHC_I_PEPTIDE_LENGTH + 1)
    scan_end = min(len(protein) - _MHC_I_PEPTIDE_LENGTH + 1, end)

    for i in range(scan_start, scan_end):
        peptide = protein[i:i + _MHC_I_PEPTIDE_LENGTH]
        if len(peptide) < _MHC_I_PEPTIDE_LENGTH:
            continue
        # Check if this 9-mer overlaps the target region
        pep_end = i + _MHC_I_PEPTIDE_LENGTH
        if pep_end <= start or i >= end:
            continue
        for allele in mhc_alleles:
            score = score_peptide_pssm(peptide, allele)
            max_score = max(max_score, score)

    return max_score


def _estimate_ddg(
    wildtype: str,
    mutant: str,
    protein: str | None = None,
    position: int | None = None,
) -> float:
    """Estimate ΔΔG from BLOSUM62 score and hydropathy change.

    Delegates to ``foldx.empirical_stability()`` when the full protein
    context is available, comparing stability before and after the
    substitution.  Falls back to a BLOSUM62-based heuristic otherwise.

    Args:
        wildtype: Wild-type amino acid.
        mutant: Mutant amino acid.
        protein: Full protein sequence (optional, for FoldX delegation).
        position: 0-based position of the mutation (optional, for FoldX delegation).

    Returns:
        Estimated ΔΔG in kcal/mol (positive = destabilizing).
    """
    # Delegate to FoldX when full protein context is available
    if protein is not None and position is not None:
        try:
            from . import foldx as _foldx
            mutated_protein = protein[:position] + mutant + protein[position + 1:]
            orig_result = _foldx.empirical_stability(protein)
            mut_result = _foldx.empirical_stability(mutated_protein)
            return round(mut_result.stability_kcal - orig_result.stability_kcal, 3)
        except ImportError:
            logger.debug("FoldX module unavailable, falling back to BLOSUM62 heuristic")
        except Exception:
            logger.debug("FoldX empirical_stability failed, falling back to BLOSUM62 heuristic")

    # Local BLOSUM62 heuristic fallback
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


def _estimate_solubility_impact(
    wildtype: str,
    mutant: str,
    protein: str | None = None,
    position: int | None = None,
) -> float:
    """Estimate solubility impact of a substitution using CamSol.

    Delegates to ``camsol.compute_intrinsic_solubility()`` when the
    full protein context is available, computing the score before and
    after the substitution.  Falls back to a hydropathy/charge heuristic.

    Negative values = decreased solubility, positive = increased.

    Args:
        wildtype: Wild-type amino acid.
        mutant: Mutant amino acid.
        protein: Full protein sequence (optional, for CamSol delegation).
        position: 0-based position of the mutation (optional, for CamSol delegation).

    Returns:
        Estimated solubility change (arbitrary units).
    """
    # Delegate to CamSol when full protein context is available
    if protein is not None and position is not None:
        try:
            from .camsol import compute_intrinsic_solubility as _camsol_score
            mutated_protein = protein[:position] + mutant + protein[position + 1:]
            orig_score = _camsol_score(protein).intrinsic_score
            mut_score = _camsol_score(mutated_protein).intrinsic_score
            return round(mut_score - orig_score, 3)
        except ImportError:
            logger.debug("CamSol module unavailable, falling back to heuristic")
        except Exception:
            logger.debug("CamSol compute_intrinsic_solubility failed, falling back to heuristic")

    # Local heuristic fallback
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


def _filter_binder_epitopes(epitopes: list[dict]) -> list[dict]:
    """Filter epitope predictions to only include strong and moderate binders."""
    return [
        e for e in epitopes
        if e.get("binding_class", "") in ("strong_binder", "moderate_binder")
    ]


# ────────────────────────────────────────────────────────────
# Public API
# ────────────────────────────────────────────────────────────

def compute_mutation_impact(
    protein: str, position: int, mutant_aa: str,
    mhc_alleles: list[str] | None = None,
) -> dict:
    """Compute the impact of a single amino acid mutation.

    Assesses the effect on MHC binding, stability, and solubility.
    Uses score_peptide_pssm from the immunogenicity module for
    MHC binding prediction.

    Args:
        protein: Protein sequence (1-letter codes).
        position: 0-based position to mutate.
        mutant_aa: Substitution amino acid.
        mhc_alleles: MHC alleles to test binding against.
            Defaults to Homo_sapiens alleles.

    Returns:
        Dictionary with keys:
          - binding_impact: list of dicts describing affected epitopes
          - stability_impact: estimated ΔΔG
          - solubility_impact: estimated solubility change
          - blosum62: BLOSUM62 score for the substitution
    """
    if mhc_alleles is None:
        mhc_alleles = _DEFAULT_MHC_ALLELES.get("Homo_sapiens", [])

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
    ddg = _estimate_ddg(wildtype, mutant_aa, protein, position)

    # Estimate solubility impact
    sol_impact = _estimate_solubility_impact(wildtype, mutant_aa, protein, position)

    # Compute binding impact: find all 9-mer epitopes overlapping this position
    binding_impact = []
    for i in range(max(0, position - _MHC_I_PEPTIDE_LENGTH + 1),
                   min(len(protein) - _MHC_I_PEPTIDE_LENGTH + 1, position + 1)):
        original_peptide = protein[i:i + _MHC_I_PEPTIDE_LENGTH]
        if len(original_peptide) < _MHC_I_PEPTIDE_LENGTH:
            continue

        # Only consider 9-mers that include the mutation position
        if not (i <= position < i + _MHC_I_PEPTIDE_LENGTH):
            continue

        # Create mutated peptide
        mut_pos_in_peptide = position - i
        mutated_peptide = (
            original_peptide[:mut_pos_in_peptide]
            + mutant_aa
            + original_peptide[mut_pos_in_peptide + 1:]
        )

        # Score before and after across all alleles
        for allele in mhc_alleles:
            orig_score = score_peptide_pssm(original_peptide, allele)
            mut_score = score_peptide_pssm(mutated_peptide, allele)

            if orig_score > 0:
                binding_impact.append({
                    "epitope": original_peptide,
                    "start": i,
                    "end": i + _MHC_I_PEPTIDE_LENGTH,
                    "allele": allele,
                    "original_binding": round(orig_score, 4),
                    "mutated_binding": round(mut_score, 4),
                    "binding_reduction": round(orig_score - mut_score, 4),
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
) -> list[MutationResult]:
    """Find mutations that disrupt a specific T-cell epitope region.

    Tries all 19 possible substitutions at each position within the
    epitope region, computes binding score change using score_peptide_pssm
    from the immunogenicity module, and filters by BLOSUM62 conservation.

    Args:
        protein: Full protein sequence.
        epitope_start: Start position (0-based) of epitope.
        epitope_end: End position (0-based, exclusive) of epitope.
        mhc_alleles: MHC alleles to test binding against.
        blosum62_min: Minimum BLOSUM62 score for acceptable substitutions.

    Returns:
        List of MutationResult objects, sorted by score
        (binding_reduction, most disruption first).
    """
    if mhc_alleles is None:
        mhc_alleles = _DEFAULT_MHC_ALLELES.get("Homo_sapiens", [])

    # Compute original binding score for the region
    original_binding = _compute_binding_score_for_region(
        protein, epitope_start, epitope_end, mhc_alleles
    )

    mutations: list[MutationResult] = []

    for pos in range(max(0, epitope_start), min(len(protein), epitope_end)):
        wildtype = protein[pos]
        if wildtype not in BLOSUM62:
            continue

        for mutant in STANDARD_AAS:
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
            ddg = _estimate_ddg(wildtype, mutant, protein, pos)

            # Estimate solubility impact
            sol_impact = _estimate_solubility_impact(wildtype, mutant, protein, pos)

            # Find the 9-mer with highest binding score that includes this position
            best_score = 0.0
            best_peptide = ""
            for pep_start in range(
                max(0, pos - _MHC_I_PEPTIDE_LENGTH + 1),
                min(len(protein) - _MHC_I_PEPTIDE_LENGTH + 1, pos + 1)
            ):
                if not (pep_start <= pos < pep_start + _MHC_I_PEPTIDE_LENGTH):
                    continue
                pep = protein[pep_start:pep_start + _MHC_I_PEPTIDE_LENGTH]
                if len(pep) < _MHC_I_PEPTIDE_LENGTH:
                    continue
                for allele in mhc_alleles:
                    sc = score_peptide_pssm(pep, allele)
                    if sc > best_score:
                        best_score = sc
                        best_peptide = pep

            mutations.append(EpitopeMutation(
                position=pos,
                original=wildtype,
                mutant=mutant,
                delta_score=round(binding_reduction, 4),
                description=f"Disrupts epitope {best_peptide} at position {pos}",
                details={
                    "epitope_disrupted": best_peptide,
                    "binding_reduction": round(binding_reduction, 4),
                    "blosum62": blosum,
                    "ddg_estimate": ddg,
                    "solubility_impact": sol_impact,
                },
            ))

    # Sort by score descending (most disruption first)
    mutations.sort(key=lambda m: m.score, reverse=True)
    return mutations


def rank_deimmunization_mutations(
    protein: str,
    mhc_alleles: list[str] | None = None,
    blosum62_min: int = 0,
) -> list[MutationResult]:
    """Find all possible deimmunization mutations across all epitopes.

    Identifies T-cell epitopes in the protein using predict_t_cell_epitopes
    from the immunogenicity module, then for each epitope finds mutations
    that reduce MHC binding. Results are ranked by a combined score
    considering binding reduction, conservation, and stability.

    Args:
        protein: Protein sequence.
        mhc_alleles: MHC alleles to test. Defaults to Homo_sapiens.
        blosum62_min: Minimum BLOSUM62 score for acceptable substitutions.

    Returns:
        Ranked list of MutationResult objects (best first).
    """
    # Use immunogenicity module's T-cell epitope data for screening
    epitopes = predict_t_cell_epitopes(protein, mhc_alleles)
    binder_epitopes = _filter_binder_epitopes(epitopes)

    if not binder_epitopes:
        return []

    all_mutations: list[MutationResult] = []

    for epitope in binder_epitopes:
        start = epitope.get("start", epitope.get("position", 0))
        end = epitope.get("end", start + _MHC_I_PEPTIDE_LENGTH)
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
    def _combined_score(m: MutationResult) -> float:
        blosum62 = m.details.get("blosum62", 0)
        ddg = m.details.get("ddg_estimate", 0.0)
        conservation_factor = 1.0 - abs(blosum62) / 4.0
        stability_factor = 1.0 - max(0.0, ddg) / 5.0
        # Ensure factors are non-negative
        conservation_factor = max(0.0, conservation_factor)
        stability_factor = max(0.0, stability_factor)
        return m.score * conservation_factor * stability_factor

    all_mutations.sort(key=_combined_score, reverse=True)

    # Deduplicate: keep the best mutation per position
    seen_positions: dict[tuple[int, str], MutationResult] = {}
    for mut in all_mutations:
        key = (mut.position, mut.mutant)
        if key not in seen_positions:
            seen_positions[key] = mut

    deduped = sorted(seen_positions.values(), key=_combined_score, reverse=True)
    return deduped


def _compute_immunogenicity_score(
    protein: str,
    organism: str = "Homo_sapiens",
    mhc_alleles: list[str] | None = None,
) -> float:
    """Compute overall immunogenicity score for a protein.

    The score is based on the average of the top epitope binding scores.
    Higher = more immunogenic.

    Args:
        protein: Protein sequence.
        organism: Target organism.
        mhc_alleles: MHC alleles to test.

    Returns:
        Immunogenicity score in [0, 1].
    """
    # Try the immunogenicity module's compute_immunogenicity first
    try:
        from .immunogenicity import compute_immunogenicity
        species = _organism_to_species(organism)
        result = compute_immunogenicity(protein, mhc_alleles=mhc_alleles, species=species)
        return result.immunogenicity_score
    except ImportError:
        logger.warning("immunogenicity.compute_immunogenicity unavailable, using local fallback")
    except Exception:
        logger.warning("immunogenicity.compute_immunogenicity failed, using local fallback")

    # Fallback: use predict_t_cell_epitopes and compute score from epitope data
    if mhc_alleles is None:
        mhc_alleles = _get_mhc_alleles(organism)

    epitopes = predict_t_cell_epitopes(protein, mhc_alleles)
    binder_epitopes = _filter_binder_epitopes(epitopes)

    if not binder_epitopes:
        return 0.0

    # Take top 5 epitope scores (or fewer if < 5)
    top_scores = [e["score"] for e in binder_epitopes[:5]]

    # Scores from predict_t_cell_epitopes are already normalized to [0, 1]
    avg_score = sum(top_scores) / len(top_scores)

    return round(min(1.0, avg_score), 4)


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
        Number of predicted epitopes (strong and moderate binders).
    """
    try:
        from .immunogenicity import compute_immunogenicity
        species = _organism_to_species(organism)
        result = compute_immunogenicity(protein, mhc_alleles=mhc_alleles, species=species)
        return result.num_t_cell_epitopes
    except ImportError:
        logger.warning("immunogenicity.compute_immunogenicity unavailable, using local fallback")
    except Exception:
        logger.warning("immunogenicity.compute_immunogenicity failed, using local fallback")

    if mhc_alleles is None:
        mhc_alleles = _get_mhc_alleles(organism)

    epitopes = predict_t_cell_epitopes(protein, mhc_alleles)
    binder_epitopes = _filter_binder_epitopes(epitopes)
    return len(binder_epitopes)


def _compute_solubility_score(protein: str) -> float:
    """Compute intrinsic solubility score using CamSol heuristic.

    Falls back to the camsol module if available.

    Args:
        protein: Protein sequence.

    Returns:
        Solubility score (arbitrary units, higher = more soluble).
    """
    try:
        from .camsol import compute_intrinsic_solubility
        result = compute_intrinsic_solubility(protein)
        # Handle SolubilityResult object (camsol module returns objects)
        if hasattr(result, 'overall_score'):
            return result.overall_score
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
    mhc_alleles: dict | None = None,
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
        mhc_alleles: Optional dict mapping allele names to their details
            (overrides organism-based defaults).

    Returns:
        DeimmunizationResult with optimized protein and detailed metrics.

    Raises:
        ImmunogenicityError: If the protein sequence is invalid.
    """
    with EngineTimer() as timer:
        # Validate protein sequence
        try:
            protein = validate_protein_sequence(protein, "Deimmunization")
        except ValueError as exc:
            raise ImmunogenicityError(str(exc)) from exc

        preserve_set = set(preserve_positions) if preserve_positions else set()
        if mhc_alleles is not None:
            allele_list = list(mhc_alleles.keys()) if isinstance(mhc_alleles, dict) else list(mhc_alleles)
        else:
            allele_list = _get_mhc_alleles(organism)

        # Compute initial metrics
        original_immunogenicity = _compute_immunogenicity_score(
            protein, organism, allele_list
        )
        original_epitope_count = _count_t_cell_epitopes(protein, organism, allele_list)

        current_protein = protein
        mutations_applied: list[dict] = []
        total_ddg = 0.0
        iteration = 0

        while iteration < max_mutations:
            # Check if target is reached
            current_score = _compute_immunogenicity_score(
                current_protein, organism, allele_list
            )
            if current_score <= target_score:
                break

            # Find T-cell epitopes using predict_t_cell_epitopes
            all_epitopes = predict_t_cell_epitopes(current_protein, allele_list)
            epitopes = _filter_binder_epitopes(all_epitopes)
            if not epitopes:
                break

            # Collect all possible mutations for the strongest epitopes
            # Process epitopes in order of strength (strongest first)
            all_candidates: list[MutationResult] = []

            for epitope in epitopes:
                start = epitope["start"]
                end = epitope["end"]

                epi_mutations = find_epitope_disrupting_mutations(
                    current_protein, start, end, allele_list, blosum62_min
                )
                all_candidates.extend(epi_mutations)

            if not all_candidates:
                logger.info("No disruptive mutations found for remaining epitopes")
                break

            # Rank by combined score
            def _combined_score(m: MutationResult) -> float:
                blosum62 = m.details.get("blosum62", 0)
                ddg = m.details.get("ddg_estimate", 0.0)
                conservation_factor = 1.0 - abs(blosum62) / 4.0
                stability_factor = 1.0 - max(0.0, ddg) / 5.0
                conservation_factor = max(0.0, conservation_factor)
                stability_factor = max(0.0, stability_factor)
                return m.score * conservation_factor * stability_factor

            all_candidates.sort(key=_combined_score, reverse=True)

            # Find the best mutation that passes all filters
            best_mutation = None
            for candidate in all_candidates:
                ddg_est = candidate.details.get("ddg_estimate", 0.0)
                blosum62_val = candidate.details.get("blosum62", 0)

                # Skip preserved positions
                if candidate.position in preserve_set:
                    continue

                # Check BLOSUM62 filter
                if blosum62_val < blosum62_min:
                    continue

                # Check ΔΔG filter
                if ddg_est >= max_ddg:
                    continue

                # Check cumulative ΔΔG
                if total_ddg + ddg_est >= max_ddg * 2:
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

            ddg_est = best_mutation.details.get("ddg_estimate", 0.0)
            total_ddg += ddg_est

            mutations_applied.append({
                "position": best_mutation.position,
                "wildtype": best_mutation.original,
                "mutant": best_mutation.mutant,
                "epitope_removed": best_mutation.details.get("epitope_disrupted", ""),
                "ddg": ddg_est,
                "blosum62": best_mutation.details.get("blosum62", 0),
                "binding_reduction": best_mutation.score,
                "solubility_impact": best_mutation.details.get("solubility_impact", 0.0),
            })

            iteration += 1
            logger.debug(
                "Iteration %d: %s%d%s (ddG=%.2f, BLOSUM62=%d, binding_reduction=%.4f)",
                iteration,
                best_mutation.original,
                best_mutation.position + 1,  # 1-based for display
                best_mutation.mutant,
                ddg_est,
                best_mutation.details.get("blosum62", 0),
                best_mutation.score,
            )

        # Compute final metrics
        optimized_immunogenicity = _compute_immunogenicity_score(
            current_protein, organism, allele_list
        )
        optimized_epitope_count = _count_t_cell_epitopes(
            current_protein, organism, allele_list
        )

        # Check stability: all ΔΔG values summed < max_ddg * len(mutations)
        stability_threshold = max_ddg * max(len(mutations_applied), 1)
        stability_preserved = total_ddg < stability_threshold

        success = optimized_immunogenicity <= target_score

    result = DeimmunizationResult(
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
        execution_time_s=round(timer.elapsed, 4),
    )

    logger.info(
        "Deimmunization complete: %d mutations, immunogenicity %.4f -> %.4f, "
        "epitopes %d -> %d, %.2fs",
        iteration,
        original_immunogenicity,
        optimized_immunogenicity,
        original_epitope_count,
        optimized_epitope_count,
        round(timer.elapsed, 4),
    )

    return result


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

    Raises:
        ImmunogenicityError: If either protein sequence is invalid.
    """
    # Validate protein sequences
    try:
        protein = validate_protein_sequence(protein, "Deimmunization")
        original_protein = validate_protein_sequence(original_protein, "Deimmunization")
    except ValueError as exc:
        raise ImmunogenicityError(str(exc)) from exc

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

    if not overall_valid:
        reasons = []
        if not immunogenicity_reduced:
            reasons.append("immunogenicity not reduced")
        if not stability_preserved:
            reasons.append("stability not preserved")
        if not solubility_preserved:
            reasons.append("solubility not preserved")
        if not all_conservative:
            reasons.append("non-conservative mutations present")
        logger.warning(
            "Deimmunized protein validation failed: %s", "; ".join(reasons)
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
