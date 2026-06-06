"""Type stubs for biocompiler.immunogenicity — MHC binding, B-cell epitope, and immunogenicity scoring."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, List, Optional, TypedDict


# ────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────

IMMUNOGENICITY_PSSM_AUC_ROC_LOW: float
IMMUNOGENICITY_PSSM_AUC_ROC_HIGH: float
IMMUNOGENICITY_NETMHCPAN_AUC_ROC_LOW: float
IMMUNOGENICITY_NETMHCPAN_AUC_ROC_HIGH: float
IMMUNOGENICITY_BCELL_AUC_ROC: float

IC50_STRONG_BINDER_THRESHOLD: float
IC50_MODERATE_BINDER_THRESHOLD: float
IC50_WEAK_BINDER_THRESHOLD: float
IC50_LOG_INTERCEPT: float
IC50_LOG_SLOPE: float

PSSM_UNKNOWN_AA_SCORE: float
PSSM_CONTRAST_POWER: float

HYDROPHOBICITY_OFFSET: float
HYDROPHOBICITY_RANGE: float

T_CELL_WEIGHT: float
B_CELL_WEIGHT: float

IMMUNOGENICITY_LOW_THRESHOLD: float
IMMUNOGENICITY_HIGH_THRESHOLD: float

MAX_DEIMMUNIZATION_CANDIDATES: int
MHC_II_CORE_LENGTH: int

CONF_EPITOPE_NEIGHBOR_CUTOPT_ANGSTROM: float
CONF_EPITOPE_MAX_NEIGHBORS: int

EPITOPE_DENSITY_CLUSTER_DISTANCE: float

DEFAULT_MHC_I_ALLELES: list[str]
DEFAULT_MHC_II_ALLELES: list[str]
POPULATION_COVERAGE: dict[str, dict[str, float]]
PRECOMPUTED_BINDERS: dict[tuple[str, str], float]

MHC_I_PSSM: dict[str, list[dict[str, float]]]
MHC_II_PSSM: dict[str, list[dict[str, float]]]

ANTIGENICITY_SCALE: dict[str, float]
PARKER_SCALE: dict[str, float]
CHOU_FASMAN_TURN: dict[str, float]
EMINI_SCALE: dict[str, float]
ALL_SCALES: dict[str, dict[str, float]]


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

class TCellEpitopeDict(TypedDict):
    allele: str
    position: int
    peptide: str
    score: float
    ic50: float
    binding_class: str


class BCellEpitopeDict(TypedDict):
    method: str
    position: int
    peptide: str
    score: float


@dataclass
class MHCBindingResult:
    allele: str
    peptide: str
    score: float
    ic50: float
    binding_class: str
    method: str


@dataclass
class MHCPredictionResult:
    mhc_i_results: List[MHCBindingResult]
    mhc_ii_results: List[MHCBindingResult]


@dataclass
class EpitopeRegion:
    start: int
    end: int
    epitope_type: str
    score: float
    details: str


@dataclass
class EpitopePredictionResult:
    t_cell_epitopes: List[TCellEpitopeDict]
    b_cell_epitopes: List[BCellEpitopeDict]


@dataclass
class ImmunogenicityResult:
    score: float
    risk_level: str
    t_cell_epitopes: List[TCellEpitopeDict]
    b_cell_epitopes: List[BCellEpitopeDict]
    recommendations: List[str]


# ────────────────────────────────────────────────────────────
# Offline prediction API
# ────────────────────────────────────────────────────────────

@dataclass
class PeptideResult:
    peptide: str
    score: float
    ic50: float
    binding_class: str
    allele: str


@dataclass
class ImmunogenicityPrediction:
    overall_score: float
    risk_level: str
    peptide_results: List[PeptideResult]
    recommendations: List[str]


# ────────────────────────────────────────────────────────────
# Core functions
# ────────────────────────────────────────────────────────────

def clear_cache() -> None: ...
def score_peptide_pssm(peptide: str, allele: str, pssm: list[dict[str, float]] | None = ...) -> float: ...
def binding_score_to_ic50(score: float) -> float: ...
def classify_binding(ic50: float) -> str: ...


# ────────────────────────────────────────────────────────────
# MHC binding prediction
# ────────────────────────────────────────────────────────────

def predict_mhc_i_binding(protein: str, alleles: list[str] = ..., peptide_length: int = ..., use_netmhcpan: bool = ...) -> List[MHCBindingResult]: ...
def predict_mhc_ii_binding(protein: str, alleles: list[str] = ..., use_netmhcpan: bool = ...) -> List[MHCBindingResult]: ...
def predict_all(protein: str, mhc_i_alleles: list[str] = ..., mhc_ii_alleles: list[str] = ..., use_netmhcpan: bool = ...) -> MHCPredictionResult: ...
def predict_t_cell_epitopes(protein: str, alleles: list[str] = ..., **kwargs: Any) -> List[TCellEpitopeDict]: ...


# ────────────────────────────────────────────────────────────
# B-cell epitope prediction
# ────────────────────────────────────────────────────────────

def predict_kolaskar_tongaonkar(protein: str, window: int = ...) -> List[BCellEpitopeDict]: ...
def predict_parker_hydrophilicity(protein: str, window: int = ...) -> List[BCellEpitopeDict]: ...
def predict_chou_fasman_beta_turn(protein: str, window: int = ...) -> List[BCellEpitopeDict]: ...
def predict_eea(protein: str, window: int = ...) -> List[BCellEpitopeDict]: ...
def predict_bepipred_like(protein: str, threshold: float = ...) -> List[BCellEpitopeDict]: ...
def predict_conformational_epitopes(pdb_content: str, cutoff: float = ...) -> List[BCellEpitopeDict]: ...
def predict_epitopes(protein: str, methods: list[str] = ..., **kwargs: Any) -> EpitopePredictionResult: ...
def compute_surface_accessibility_approx(protein: str) -> list[float]: ...
def predict_b_cell_epitopes(protein: str, **kwargs: Any) -> List[BCellEpitopeDict]: ...


# ────────────────────────────────────────────────────────────
# Combined immunogenicity scoring
# ────────────────────────────────────────────────────────────

def compute_immunogenicity(protein: str, alleles: list[str] = ..., use_netmhcpan: bool = ..., **kwargs: Any) -> ImmunogenicityResult: ...
def find_deimmunization_mutations(protein: str, alleles: list[str] = ..., max_mutations: int = ..., use_netmhcpan: bool = ..., **kwargs: Any) -> List[dict[str, Any]]: ...
def suggest_mutations(protein: str, alleles: list[str] = ..., **kwargs: Any) -> List[dict[str, Any]]: ...
def compute_epitope_density(protein: str, alleles: list[str] = ...) -> float: ...
def compute_immunogenicity_batch(proteins: list[str], alleles: list[str] = ..., **kwargs: Any) -> List[ImmunogenicityResult]: ...


# ────────────────────────────────────────────────────────────
# Offline prediction API
# ────────────────────────────────────────────────────────────

def predict_immunogenicity(protein: str, alleles: list[str] = ..., **kwargs: Any) -> ImmunogenicityPrediction: ...
def scan_peptides(protein: str, alleles: list[str] = ..., peptide_length: int = ...) -> List[PeptideResult]: ...
