"""Data classes and TypedDicts for immunogenicity prediction results.

Split out of ``core.py`` (W8-a refactor).
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List, Optional, TypedDict

from ..engines.base import BaseEngineResult, MutationResult
from biocompiler.shared.types import Verdict

# ``MHCPredictionResult.population_coverage`` reads the POPULATION_COVERAGE
# table defined in ``_constants``; import it (and siblings) for safety.
from ._constants import *  # noqa: F401,F403
from ._constants import POPULATION_COVERAGE  # noqa: F401

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# MHC binding: data classes
# ═══════════════════════════════════════════════════════════════════════════


@dataclass
class TCellEpitopeDict(TypedDict, total=False):
    """Structured type for T-cell epitope prediction results."""

    start: int
    end: int
    peptide: str
    score: float
    allele: str
    binding_class: str


class BCellEpitopeDict(TypedDict, total=False):
    """Structured type for B-cell epitope prediction results (backward-compat)."""

    start: int
    end: int
    peptide: str
    score: float
    method: str
    antigenic: bool


@dataclass
class MHCBindingResult:
    """Result of a single peptide-MHC binding prediction."""

    allele: str
    peptide: str
    start_position: int
    end_position: int
    binding_score: float  # 0 to 1, 1 = strong binder
    ic50_nm: float | None  # estimated IC50 in nM, if computable
    binding_class: str  # "strong_binder" | "moderate_binder" | "weak_binder" | "non_binder"
    anchor_residues: dict[int, str]  # position -> AA at anchor positions
    anchor_scores: dict[int, float]  # position -> binding contribution
    method: str = "pssm_fallback"  # prediction method: "pssm_fallback" | "mhcflurry" | "mhcflurry_presentation" | "precomputed_lookup" | "netmhcpan"
    rank: Optional[float] = None  # percentile rank from prediction tool, if available
    confidence: float = 0.5  # confidence in prediction: 0.95 (netmhcpan), 0.85 (mhcflurry), 0.7 (precomputed), 0.5 (pssm_fallback)

    # Removed: score (alias for binding_score), position (alias
    # for start_position) backward-compat properties.


@dataclass
class MHCPredictionResult:
    """Aggregated MHC binding prediction for a protein."""

    protein: str
    mhc_i_results: list[MHCBindingResult]
    mhc_ii_results: list[MHCBindingResult]
    strong_binders: int
    weak_binders: int
    non_binders: int
    binding_profile: dict[str, float]  # allele -> max binding score
    # EngineResult protocol fields
    method: str = "immunogenicity_pssm"
    success: bool = True
    error: str | None = None
    execution_time_s: float = 0.0

    @property
    def predictions(self) -> list[MHCBindingResult]:
        """All predictions combined (MHC-I + MHC-II)."""
        return self.mhc_i_results + self.mhc_ii_results

    # Removed: class_i_results / class_ii_results aliases —
    # use mhc_i_results / mhc_ii_results directly.

    @property
    def binders(self) -> list[MHCBindingResult]:
        """All results classified as strong or moderate binders."""
        return [
            r for r in self.predictions
            if r.binding_class in ("strong_binder", "moderate_binder")
        ]

    @property
    def binding_rate(self) -> float:
        """Fraction of peptides that are binders (strong or moderate)."""
        total = len(self.predictions)
        if total == 0:
            return 0.0
        return len(self.binders) / total

    @property
    def population_coverage(self) -> float:
        """Estimated population coverage based on binding profile.

        Uses the fraction of alleles that have at least one strong
        or moderate binder (IC50 < 500 nM, i.e. binding_score > 0.5),
        weighted by population frequency.
        """
        if not self.binding_profile:
            return 0.0
        total_coverage = 0.0
        for allele, max_score in self.binding_profile.items():
            if max_score > 0.5:
                cov = POPULATION_COVERAGE.get(allele, {})
                freq = cov.get("Caucasian", 0.0)
                total_coverage += freq / 100.0
        return min(1.0, total_coverage)


# ═══════════════════════════════════════════════════════════════════════════
# MHC binding: utility functions
# ═══════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════


@dataclass
class PeptideResult:
    """Result of scoring a single peptide against an MHC allele."""
    position: int
    peptide: str
    ic50_nm: float
    binding_class: str


@dataclass
class ImmunogenicityPrediction:
    """Result of predicting immunogenicity for a single peptide-allele pair."""
    allele: str
    peptide: str
    ic50_nm: float
    binding_class: str
    method: str
    confidence: str


@dataclass
class EpitopeRegion:
    """A predicted B-cell epitope region."""

    start: int          # 0-indexed start position
    end: int            # exclusive end position
    peptide: str
    score: float        # 0 to 1
    method: str         # "kolaskar_tongaonkar", "bepipred", "parker_hydrophilicity",
                        # "chou_fasman", "eea", "consensus", "conformational"
    is_linear: bool
    properties: dict = field(default_factory=dict)


@dataclass
class EpitopePredictionResult:
    """Combined B-cell epitope prediction result."""

    protein: str
    linear_epitopes: list[EpitopeRegion]
    conformational_epitopes: list[EpitopeRegion]   # empty if no structure provided
    per_residue_score: list[float]                  # combined epitope propensity per residue
    epitope_coverage: float                         # fraction of residues in predicted epitopes
    methods_used: list[str]


# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: internal helpers
# ═══════════════════════════════════════════════════════════════════════════


@dataclass
class ImmunogenicityResult(BaseEngineResult):
    """Result of immunogenicity scoring for a protein sequence.

    Inherits from :class:`BaseEngineResult` for unified API access.
    Domain-specific field aliases are preserved as properties for
    backward compatibility.

    Unified fields (from BaseEngineResult):
      - sequence: protein sequence (alias: protein)
      - primary_score: overall immunogenicity score (alias: immunogenicity_score, overall_score)
      - classification: risk category (alias: risk_class, immunogenicity_class)
      - mutations: deimmunization suggestions (alias: deimmunization_mutations, deimmunization_candidates)
      - engine_name: "immunogenicity"
      - primary_score_label: "immunogenicity"

    Engine-specific fields:
      - t_cell_score: T-cell epitope contribution
      - b_cell_score: B-cell epitope contribution
      - t_cell_epitopes: predicted T-cell epitopes
      - b_cell_epitopes: predicted B-cell epitopes
    """

    # Engine-specific fields
    t_cell_score: float = 0.0  # T-cell epitope contribution
    b_cell_score: float = 0.0  # B-cell epitope contribution
    t_cell_epitopes: List[TCellEpitopeDict] = field(default_factory=list)  # predicted T-cell epitopes
    b_cell_epitopes: List[BCellEpitopeDict] = field(default_factory=list)  # predicted B-cell epitopes

    # Override BaseEngineResult fields with defaults for convenience
    success: bool = True
    error: Optional[str] = None
    execution_time_s: float = 0.0
    engine_name: str = "immunogenicity"
    primary_score_label: str = "immunogenicity"
    mutations: List[MutationResult] = field(default_factory=list)

    # ── Honesty fields (TIGHTEN-4) ─────────────────────────────────
    # The default PSSMs in this module are hand-crafted approximations,
    # NOT real binding data.  Every ImmunogenicityResult therefore carries
    # a `verdict` indicating whether the scores are trustworthy enough
    # to issue a PASS/FAIL claim, plus a human-readable `message` and
    # a machine-readable `reason` / `data_source` pair.
    #
    # Default state: UNCERTAIN / fabricated_scores / guessed_pssm.
    # Only `compute_immunogenicity(use_real_data=True)` with a real
    # NetMHCpan / MHCflurry backend can produce PASS / FAIL here.
    verdict: Verdict = Verdict.UNCERTAIN
    reason: str = "fabricated_scores"
    message: str = (
        "Immunogenicity scores use approximate PSSMs, NOT real binding data. "
        "Verdict is UNCERTAIN. Install NetMHCpan or MHCflurry for verified "
        "predictions."
    )
    data_source: str = "guessed_pssm"
    scores: dict = field(default_factory=dict)

    # ---- Backward-compatible property aliases ----

    @property
    def protein(self) -> str:
        """Alias for sequence (backward compat)."""
        return self.sequence

    @protein.setter
    def protein(self, value: str) -> None:
        self.sequence = value

    @property
    def overall_score(self) -> float:
        """Alias for primary_score (backward compat)."""
        return self.primary_score

    @overall_score.setter
    def overall_score(self, value: float) -> None:
        self.primary_score = value

    @property
    def immunogenicity_score(self) -> float:
        """Alias for primary_score."""
        return self.primary_score

    @immunogenicity_score.setter
    def immunogenicity_score(self, value: float) -> None:
        self.primary_score = value

    @property
    def immunogenicity_class(self) -> str:
        """Alias for classification (backward compat)."""
        return self.classification

    @immunogenicity_class.setter
    def immunogenicity_class(self, value: str) -> None:
        self.classification = value

    @property
    def risk_class(self) -> str:
        """Alias for classification."""
        return self.classification

    @risk_class.setter
    def risk_class(self, value: str) -> None:
        self.classification = value

    @property
    def deimmunization_candidates(self) -> List[MutationResult]:
        """Alias for mutations (backward compat)."""
        return self.mutations

    @deimmunization_candidates.setter
    def deimmunization_candidates(self, value: List[MutationResult]) -> None:
        self.mutations = value

    @property
    def deimmunization_mutations(self) -> List[MutationResult]:
        """Alias for mutations."""
        return self.mutations

    @deimmunization_mutations.setter
    def deimmunization_mutations(self, value: List[MutationResult]) -> None:
        self.mutations = value

    @property
    def method(self) -> str:
        """Alias for engine_name (backward compat)."""
        return self.engine_name

    @method.setter
    def method(self, value: str) -> None:
        self.engine_name = value

    @property
    def confidence_level(self) -> str:
        """Accuracy confidence level for the immunogenicity prediction.

        Returns one of:
          - ``"high"`` -- NetMHCpan mode (AUC-ROC 0.85-0.95)
          - ``"medium"`` -- PSSM mode with strong anchor matches
          - ``"low"`` -- PSSM mode or B-cell epitope only
        """
        method = self.engine_name
        if "netmhcpan" in method.lower():
            return "high"
        # PSSM-based — check if we have any strong binders to assess
        if self.t_cell_epitopes:
            # If we have epitope predictions, at least PSSM found something
            return "medium"
        return "low"



__all__ = [
    "MHCBindingResult",
    "MHCPredictionResult",
    "TCellEpitopeDict",
    "BCellEpitopeDict",
    "PeptideResult",
    "ImmunogenicityPrediction",
    "EpitopeRegion",
    "EpitopePredictionResult",
    "ImmunogenicityResult",
]
