"""Immunogenicity scoring, MHC binding prediction, and B-cell epitope prediction.

This module consolidates three formerly separate modules:

* **MHC binding** (formerly ``mhc_binding.py``): Predicts peptide-MHC
  binding affinity using position-specific scoring matrices (PSSMs) derived
  from known binding motifs in the Immune Epitope Database (IEDB).
  When ``use_netmhcpan=True``, the NetMHCpan web API is tried first
  for more accurate predictions; the PSSM heuristic serves as fallback.
* **B-cell epitope prediction** (formerly ``epitope.py``): Linear and
  conformational B-cell epitope prediction using multiple classical scales
  and methods (Kolaskar-Tongaonkar, Parker hydrophilicity, Chou-Fasman
  beta-turn, Emini surface accessibility, BepiPred-like composite, and
  conformational epitope prediction from PDB structure).
* **Immunogenicity scoring** (original ``immunogenicity.py``): Combined
  T-cell / B-cell immunogenicity scoring and deimmunization mutation
  suggestions.

Accuracy and Confidence
----------------------
**PSSM-based MHC binding prediction** (default, offline mode):
  - Expected AUC-ROC: 0.60–0.75 for MHC-I binding classification
  - This is significantly below state-of-the-art methods
  - PSSMs capture anchor position preferences but miss subtle
    peptide-MHC interaction features
  - IC50 estimates are rough approximations (log-linear mapping)
  - Binding classification thresholds (50/500/5000 nM) are standard
    but PSSM-derived IC50 values have high uncertainty

**NetMHCpan-based prediction** (when ``use_netmhcpan=True``):
  - Expected AUC-ROC: 0.85–0.95 for MHC-I binding (NetMHCpan 4.1)
  - This is the gold standard for computational MHC binding prediction
  - Requires API connectivity to the NetMHCpan web service
  - Falls back to PSSM if API is unavailable

**B-cell epitope prediction:**
  - Classical scale-based methods (Kolaskar-Tongaonkar, Parker, etc.)
    have typical AUC-ROC of 0.55–0.65
  - Performance varies significantly by epitope type and protein
  - Conformational epitope prediction (when PDB available) is more
    reliable than linear epitope prediction

**Deimmunization mutation suggestions:**
  - Confidence depends on the underlying binding prediction method
  - PSSM-based suggestions: **LOW** confidence
  - NetMHCpan-based suggestions: **MEDIUM-HIGH** confidence
  - Always verify experimentally before clinical use

**Upgrade path:**
  - Replace PSSMs with a neural network-based method for offline use
  - Add MHCflurry as an alternative offline predictor
  - Integrate BepiPred-2.0 for B-cell epitope prediction

  **Confidence levels:**
    - NetMHCpan mode: **HIGH** for MHC-I, **MEDIUM** for MHC-II
    - PSSM mode, strong anchor matches: **MEDIUM**
    - PSSM mode, weak anchor matches: **LOW**
    - B-cell epitope (linear): **LOW**
    - B-cell epitope (conformational with PDB): **MEDIUM**

All predictions are sequence-based heuristics and do not replace
experimental validation. When NetMHCpan integration is enabled
(``use_netmhcpan=True``), predictions use the NetMHCpan 4.1 API,
which provides significantly more accurate binding affinity estimates
than the PSSM-based approach.

References
----------
- Kolaskar & Tongaonkar, FEBS Lett 1990; 276:172-174
- Parker et al., Biochemistry 1986; 25:5424-5432
- Chou & Fasman, Biochemistry 1974; 13:222-245
- Emini et al., J Virol 1985; 55:836-839
- Reynisson et al., Nucleic Acids Res 2020; 48:W449 (NetMHCpan 4.1)
- O'Donnell et al., Bioinformatics 2018; 34:2696 (MHCflurry)
- Jespersen et al., Nucleic Acids Res 2017; 45:W39 (BepiPred-2.0)
"""
from __future__ import annotations

import concurrent.futures
import hashlib
import logging
import math
from dataclasses import dataclass, field

from typing import List, Optional, TypedDict

from .constants import BLOSUM62, DEFAULT_MHC_PEPTIDE_LENGTH, HYDROPATHY, STANDARD_AAS

from .engine_base import (
    BaseEngineResult,
    BatchResult,
    EngineTimer,
    MutationResult,
    validate_protein_sequence,
)
from .exceptions import ImmunogenicityError

__all__ = [
    "ImmunogenicityError",
    "MHCBindingResult",
    "MHCPredictionResult",
    "ImmunogenicityResult",
    "EpitopeRegion",
    "EpitopePredictionResult",
    "TCellEpitopeDict",
    "BCellEpitopeDict",
    "DEFAULT_MHC_I_ALLELES",
    "DEFAULT_MHC_II_ALLELES",
    "POPULATION_COVERAGE",
    "ANTIGENICITY_SCALE",
    "PARKER_SCALE",
    "CHOU_FASMAN_TURN",
    "EMINI_SCALE",
    "ALL_SCALES",
    "MHC_I_PSSM",
    "MHC_II_PSSM",
    "clear_cache",
    "score_peptide_pssm",
    "binding_score_to_ic50",
    "classify_binding",
    "predict_mhc_i_binding",
    "predict_mhc_ii_binding",
    "predict_all",
    "predict_t_cell_epitopes",
    "predict_kolaskar_tongaonkar",
    "predict_parker_hydrophilicity",
    "predict_chou_fasman_beta_turn",
    "predict_eea",
    "predict_bepipred_like",
    "predict_conformational_epitopes",
    "predict_epitopes",
    "compute_surface_accessibility_approx",
    "predict_b_cell_epitopes",
    "compute_immunogenicity",
    "find_deimmunization_mutations",
    "compute_immunogenicity_batch",
    "IMMUNOGENICITY_PSSM_AUC_ROC_LOW",
    "IMMUNOGENICITY_PSSM_AUC_ROC_HIGH",
    "IMMUNOGENICITY_NETMHCPAN_AUC_ROC_LOW",
    "IMMUNOGENICITY_NETMHCPAN_AUC_ROC_HIGH",
    "IMMUNOGENICITY_BCELL_AUC_ROC",
    # Binding classification thresholds (nM)
    "IC50_STRONG_BINDER_THRESHOLD",
    "IC50_MODERATE_BINDER_THRESHOLD",
    "IC50_WEAK_BINDER_THRESHOLD",
    # IC50 mapping constants
    "IC50_LOG_INTERCEPT",
    "IC50_LOG_SLOPE",
    # PSSM scoring constants
    "PSSM_UNKNOWN_AA_SCORE",
    "PSSM_CONTRAST_POWER",
    # Hydrophobicity normalization
    "HYDROPHOBICITY_OFFSET",
    "HYDROPHOBICITY_RANGE",
    # Immunogenicity scoring weights
    "T_CELL_WEIGHT",
    "B_CELL_WEIGHT",
    # Immunogenicity classification thresholds
    "IMMUNOGENICITY_LOW_THRESHOLD",
    "IMMUNOGENICITY_HIGH_THRESHOLD",
    # Deimmunization limits
    "MAX_DEIMMUNIZATION_CANDIDATES",
    # MHC-II core peptide length
    "MHC_II_CORE_LENGTH",
    # Conformational epitope constants
    "CONF_EPITOPE_NEIGHBOR_CUTOPT_ANGSTROM",
    "CONF_EPITOPE_MAX_NEIGHBORS",
]

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# Accuracy constants
# ═══════════════════════════════════════════════════════════════════════════

#: Expected AUC-ROC lower bound for PSSM-based MHC binding prediction
#: PSSMs capture anchor preferences but miss subtle interaction features
IMMUNOGENICITY_PSSM_AUC_ROC_LOW: float = 0.60

#: Expected AUC-ROC upper bound for PSSM-based MHC binding prediction
IMMUNOGENICITY_PSSM_AUC_ROC_HIGH: float = 0.75

#: Expected AUC-ROC lower bound for NetMHCpan-based MHC binding prediction
#: NetMHCpan 4.1 (Reynisson et al., Nucleic Acids Res 2020)
IMMUNOGENICITY_NETMHCPAN_AUC_ROC_LOW: float = 0.85

#: Expected AUC-ROC upper bound for NetMHCpan-based MHC binding prediction
IMMUNOGENICITY_NETMHCPAN_AUC_ROC_HIGH: float = 0.95

#: Expected AUC-ROC for B-cell epitope prediction (classical scales)
#: Linear epitope methods typically perform poorly
IMMUNOGENICITY_BCELL_AUC_ROC: float = 0.60

# ═══════════════════════════════════════════════════════════════════════════
# MHC binding classification & scoring constants
# ═══════════════════════════════════════════════════════════════════════════

#: IC50 threshold (nM) for strong binder classification
IC50_STRONG_BINDER_THRESHOLD: float = 50.0

#: IC50 threshold (nM) for moderate binder classification
IC50_MODERATE_BINDER_THRESHOLD: float = 500.0

#: IC50 threshold (nM) for weak binder classification
IC50_WEAK_BINDER_THRESHOLD: float = 5000.0

#: Log-linear IC50 mapping intercept (IC50 = 10^(intercept - slope * score))
IC50_LOG_INTERCEPT: float = 3.949

#: Log-linear IC50 mapping slope
IC50_LOG_SLOPE: float = 2.5

#: Default PSSM score for amino acids not present in a position row
PSSM_UNKNOWN_AA_SCORE: float = 0.3

#: Power exponent applied to normalised PSSM score to increase contrast
PSSM_CONTRAST_POWER: float = 2.0

# ═══════════════════════════════════════════════════════════════════════════
# Hydrophobicity & immunogenicity scoring constants
# ═══════════════════════════════════════════════════════════════════════════

#: Offset used to normalise hydrophobicity scores to [0, 1]
HYDROPHOBICITY_OFFSET: float = 4.5

#: Range divisor used to normalise hydrophobicity scores to [0, 1]
HYDROPHOBICITY_RANGE: float = 9.0

#: Weight of T-cell epitope contribution in overall immunogenicity score
T_CELL_WEIGHT: float = 0.6

#: Weight of B-cell epitope contribution in overall immunogenicity score
B_CELL_WEIGHT: float = 0.4

#: Overall immunogenicity score below which risk class is "low"
IMMUNOGENICITY_LOW_THRESHOLD: float = 0.3

#: Overall immunogenicity score at or above which risk class is "high"
IMMUNOGENICITY_HIGH_THRESHOLD: float = 0.6

# ═══════════════════════════════════════════════════════════════════════════
# Deimmunization & MHC-II constants
# ═══════════════════════════════════════════════════════════════════════════

#: Maximum number of deimmunization mutation candidates to return
MAX_DEIMMUNIZATION_CANDIDATES: int = 200

#: Core peptide length for MHC-II binding prediction (9-mer core)
MHC_II_CORE_LENGTH: int = 9

# ═══════════════════════════════════════════════════════════════════════════
# Conformational epitope prediction constants
# ═══════════════════════════════════════════════════════════════════════════

#: C-alpha neighbor distance cutoff (Angstrom) for surface residue identification
CONF_EPITOPE_NEIGHBOR_CUTOPT_ANGSTROM: float = 12.0

#: Maximum C-alpha neighbors before a residue is classified as buried
CONF_EPITOPE_MAX_NEIGHBORS: int = 15

# ═══════════════════════════════════════════════════════════════════════════
# Amino-acid constants — derived from shared constants.py
# ═══════════════════════════════════════════════════════════════════════════

_STANDARD_AA_SET: set[str] = set(STANDARD_AAS)


# ═══════════════════════════════════════════════════════════════════════════
# MHC allele defaults and population coverage
# ═══════════════════════════════════════════════════════════════════════════

DEFAULT_MHC_I_ALLELES: list[str] = [
    "HLA-A*02:01",
    "HLA-A*01:01",
    "HLA-A*03:01",
    "HLA-A*24:02",
    "HLA-B*07:02",
    "HLA-B*08:01",
]

DEFAULT_MHC_II_ALLELES: list[str] = [
    "HLA-DRB1*01:01",
    "HLA-DRB1*04:01",
    "HLA-DRB1*07:01",
]

POPULATION_COVERAGE: dict[str, dict[str, float]] = {
    "HLA-A*02:01": {
        "Caucasian": 28.0,
        "Asian": 10.0,
        "African": 5.0,
        "Hispanic": 18.0,
    },
    "HLA-A*01:01": {
        "Caucasian": 16.0,
        "Asian": 2.0,
        "African": 3.0,
        "Hispanic": 8.0,
    },
    "HLA-A*03:01": {
        "Caucasian": 14.0,
        "Asian": 3.0,
        "African": 4.0,
        "Hispanic": 7.0,
    },
    "HLA-A*24:02": {
        "Caucasian": 10.0,
        "Asian": 15.0,
        "African": 3.0,
        "Hispanic": 12.0,
    },
    "HLA-B*07:02": {
        "Caucasian": 12.0,
        "Asian": 4.0,
        "African": 6.0,
        "Hispanic": 6.0,
    },
    "HLA-B*08:01": {
        "Caucasian": 10.0,
        "Asian": 1.0,
        "African": 3.0,
        "Hispanic": 4.0,
    },
    "HLA-DRB1*01:01": {
        "Caucasian": 10.0,
        "Asian": 5.0,
        "African": 4.0,
        "Hispanic": 6.0,
    },
    "HLA-DRB1*04:01": {
        "Caucasian": 15.0,
        "Asian": 8.0,
        "African": 3.0,
        "Hispanic": 12.0,
    },
    "HLA-DRB1*07:01": {
        "Caucasian": 17.0,
        "Asian": 6.0,
        "African": 8.0,
        "Hispanic": 10.0,
    },
}

# ═══════════════════════════════════════════════════════════════════════════
# MHC binding: PSSM construction
# ═══════════════════════════════════════════════════════════════════════════


def _make_pssm_row(
    preferred: dict[str, float] | None = None,
    disfavored: dict[str, float] | None = None,
    default: float = 1.0,
) -> dict[str, float]:
    """Build a single PSSM position row.

    Parameters
    ----------
    preferred : dict mapping AA -> score for preferred residues
    disfavored : dict mapping AA -> score for disfavored residues
    default : score for residues not mentioned in *preferred* or *disfavored*
    """
    row: dict[str, float] = {aa: default for aa in STANDARD_AAS}
    if preferred:
        for aa, score in preferred.items():
            if aa in row:
                row[aa] = score
    if disfavored:
        for aa, score in disfavored.items():
            if aa in row:
                row[aa] = score
    return row


def _build_mhc_i_pssms() -> dict[str, list[dict[str, float]]]:
    """Construct PSSMs for common MHC-I alleles."""

    pssms: dict[str, list[dict[str, float]]] = {}

    # --- HLA-A*02:01 ---
    pssms["HLA-A*02:01"] = [
        _make_pssm_row(
            preferred={"L": 1.2, "M": 1.2, "I": 1.2, "V": 1.2, "A": 1.1, "F": 1.1},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.5, "R": 0.5},
        ),
        _make_pssm_row(
            preferred={"L": 2.0, "M": 2.0, "I": 1.8, "V": 1.8},
            disfavored={"D": 0.3, "E": 0.3, "K": 0.3, "R": 0.3, "P": 0.4},
            default=0.8,
        ),
        _make_pssm_row(
            preferred={"L": 1.1, "V": 1.1, "A": 1.1},
            disfavored={"P": 0.6},
        ),
        _make_pssm_row(
            preferred={"K": 1.1, "R": 1.1},
            disfavored={"P": 0.6},
        ),
        _make_pssm_row(
            preferred={"A": 1.1, "V": 1.1, "I": 1.1},
            disfavored={"P": 0.6},
        ),
        _make_pssm_row(
            preferred={"V": 1.1, "I": 1.1, "L": 1.1},
            disfavored={"P": 0.6},
        ),
        _make_pssm_row(
            preferred={"L": 1.1, "I": 1.1, "V": 1.1},
            disfavored={"P": 0.6},
        ),
        _make_pssm_row(
            preferred={"A": 1.1, "V": 1.1, "L": 1.1},
            disfavored={"P": 0.6},
        ),
        _make_pssm_row(
            preferred={"V": 1.5, "L": 1.5, "I": 1.3, "A": 1.2},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.4, "R": 0.4, "P": 0.4},
            default=0.8,
        ),
    ]

    # --- HLA-A*01:01 ---
    pssms["HLA-A*01:01"] = [
        _make_pssm_row(
            preferred={"A": 1.1, "S": 1.1},
            disfavored={"W": 0.5, "R": 0.5},
        ),
        _make_pssm_row(
            preferred={"T": 1.8, "S": 1.6, "D": 1.5, "E": 1.5},
            disfavored={"L": 0.4, "I": 0.4, "V": 0.5, "F": 0.4},
            default=0.8,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"Y": 1.8, "F": 1.6},
            disfavored={"K": 0.4, "R": 0.4, "D": 0.5, "E": 0.5},
            default=0.8,
        ),
    ]

    # --- HLA-A*03:01 ---
    pssms["HLA-A*03:01"] = [
        _make_pssm_row(
            preferred={"A": 1.1, "S": 1.1},
        ),
        _make_pssm_row(
            preferred={"V": 1.8, "I": 1.8, "L": 1.6, "M": 1.6},
            disfavored={"D": 0.4, "E": 0.4, "N": 0.5, "Q": 0.5},
            default=0.8,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 2.0, "R": 1.8, "H": 1.4},
            disfavored={"D": 0.3, "E": 0.3, "S": 0.5, "T": 0.5},
            default=0.7,
        ),
    ]

    # --- HLA-A*24:02 ---
    pssms["HLA-A*24:02"] = [
        _make_pssm_row(
            preferred={"Y": 1.2, "F": 1.1},
        ),
        _make_pssm_row(
            preferred={"Y": 2.0, "F": 2.0, "W": 1.8},
            disfavored={"D": 0.3, "E": 0.3, "K": 0.3, "R": 0.3, "P": 0.4},
            default=0.7,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"F": 1.5, "L": 1.5, "I": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.4, "R": 0.4},
            default=0.8,
        ),
    ]

    # --- HLA-B*07:02 ---
    pssms["HLA-B*07:02"] = [
        _make_pssm_row(
            preferred={"A": 1.1, "P": 1.1},
        ),
        _make_pssm_row(
            preferred={"P": 2.0, "A": 1.8},
            disfavored={"D": 0.3, "E": 0.3, "K": 0.3, "R": 0.3, "W": 0.4},
            default=0.7,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.5, "I": 1.5, "V": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.4, "R": 0.4},
            default=0.8,
        ),
    ]

    # --- HLA-B*08:01 ---
    pssms["HLA-B*08:01"] = [
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 1.8, "R": 1.8},
            disfavored={"D": 0.3, "E": 0.3, "P": 0.4, "G": 0.5},
            default=0.8,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.5, "I": 1.3, "V": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.5},
            default=0.8,
        ),
    ]

    return pssms


def _build_mhc_ii_pssms() -> dict[str, list[dict[str, float]]]:
    """Construct PSSMs for common MHC-II alleles (core 9-mer)."""

    pssms: dict[str, list[dict[str, float]]] = {}

    # --- HLA-DRB1*01:01 ---
    pssms["HLA-DRB1*01:01"] = [
        _make_pssm_row(
            preferred={"F": 1.8, "Y": 1.7, "W": 1.6, "L": 1.5, "I": 1.4, "V": 1.4, "M": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.5, "R": 0.5},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.6, "S": 1.4, "T": 1.4, "N": 1.3, "G": 1.2},
            disfavored={"W": 0.5, "F": 0.6, "Y": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.5, "I": 1.4, "V": 1.4, "M": 1.3, "F": 1.3},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.5},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 1.3, "R": 1.3, "N": 1.2, "Q": 1.2, "E": 1.1, "D": 1.1},
        ),
    ]

    # --- HLA-DRB1*04:01 ---
    pssms["HLA-DRB1*04:01"] = [
        _make_pssm_row(
            preferred={"F": 1.8, "Y": 1.7, "W": 1.6, "L": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.5, "R": 0.5},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"D": 1.8, "E": 1.6},
            disfavored={"K": 0.4, "R": 0.4, "W": 0.5},
            default=0.8,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.6, "S": 1.4, "G": 1.3, "N": 1.2},
            disfavored={"W": 0.5, "F": 0.6, "Y": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.4, "I": 1.3, "V": 1.3, "F": 1.3},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.5},
            default=0.9,
        ),
    ]

    # --- HLA-DRB1*07:01 ---
    pssms["HLA-DRB1*07:01"] = [
        _make_pssm_row(
            preferred={"F": 1.6, "Y": 1.5, "L": 1.4, "I": 1.3, "V": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.5, "R": 0.5},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.4, "S": 1.3, "T": 1.3, "N": 1.2},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.4, "I": 1.3, "V": 1.3, "F": 1.3},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.5},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 1.3, "R": 1.3, "N": 1.2, "Q": 1.2},
        ),
    ]

    return pssms


# ═══════════════════════════════════════════════════════════════════════════
# Lazy PSSM initialization
# ═══════════════════════════════════════════════════════════════════════════

MHC_I_PSSM: dict[str, list[dict[str, float]]] = {}
MHC_II_PSSM: dict[str, list[dict[str, float]]] = {}
_pssm_built: bool = False


def _ensure_pssms_built() -> None:
    """Build PSSM dicts on first access (lazy initialization)."""
    global MHC_I_PSSM, MHC_II_PSSM, _pssm_built
    if not _pssm_built:
        MHC_I_PSSM = _build_mhc_i_pssms()
        MHC_II_PSSM = _build_mhc_ii_pssms()
        _pssm_built = True


def _get_mhc_i_pssms() -> dict[str, list[dict[str, float]]]:
    """Return MHC-I PSSMs, building them lazily on first call."""
    _ensure_pssms_built()
    return MHC_I_PSSM


def _get_mhc_ii_pssms() -> dict[str, list[dict[str, float]]]:
    """Return MHC-II PSSMs, building them lazily on first call."""
    _ensure_pssms_built()
    return MHC_II_PSSM


# Eagerly build PSSMs at module load so MHC_I_PSSM / MHC_II_PSSM are
# populated when tests access them directly.
_ensure_pssms_built()


# Removed in v7.5.0: MHC_I_PREFERENCES, MHC_II_PREFERENCES (backward-compat
# aliases derived from PSSMs — use MHC_I_PSSM / MHC_II_PSSM instead).
# Removed in v7.5.0: _DEFAULT_MHC_I_ALLELES, _DEFAULT_MHC_II_ALLELES,
# _DEFAULT_MHC_ALLELES (private backward-compat aliases).

# ═══════════════════════════════════════════════════════════════════════════
# Prediction cache
# ═══════════════════════════════════════════════════════════════════════════

_prediction_cache: dict[tuple[str, str, int], list[MHCBindingResult]] = {}


def clear_cache() -> None:
    """Clear the MHC binding prediction cache.

    PSSMs are built once and kept across cache clears.
    """
    global _prediction_cache
    _prediction_cache.clear()
    logger.info("Immunogenicity prediction cache cleared")


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

    # Removed in v7.5.0: score (alias for binding_score), position (alias
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

    # Removed in v7.5.0: class_i_results / class_ii_results aliases —
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


def score_peptide_pssm(
    peptide: str,
    pssm: list[dict[str, float]] | str,
) -> float:
    """Compute binding score from a PSSM.

    Uses the geometric mean of position-specific scores and normalises
    to the [0, 1] range.

    Parameters
    ----------
    peptide : str
        Amino-acid sequence of the peptide (must equal PSSM length).
    pssm : list[dict[str, float]] or str
        Position-specific scoring matrix, one dict per position.
        Alternatively, an allele name string (e.g. ``"HLA-A*02:01"``)
        which will be looked up in :data:`MHC_I_PSSM` and
        :data:`MHC_II_PSSM`.

    Returns
    -------
    float
        Normalised binding score in [0, 1].
    """
    if isinstance(pssm, str):
        allele = pssm
        lookup = _get_mhc_i_pssms().get(allele) or _get_mhc_ii_pssms().get(allele)
        if lookup is None:
            logger.debug("No PSSM for allele %s — returning 0.0", allele)
            return 0.0
        pssm = lookup

    if len(peptide) != len(pssm):
        logger.warning(
            "Peptide length %d does not match PSSM length %d — returning 0.0",
            len(peptide),
            len(pssm),
        )
        return 0.0

    scores: list[float] = []
    for i, aa in enumerate(peptide):
        aa_upper = aa.upper()
        if aa_upper not in pssm[i]:
            scores.append(PSSM_UNKNOWN_AA_SCORE)
        else:
            scores.append(pssm[i][aa_upper])

    log_sum = sum(math.log(max(s, 1e-10)) for s in scores)
    geo_mean = math.exp(log_sum / len(scores))

    max_scores: list[float] = []
    min_scores: list[float] = []
    for pos_dict in pssm:
        vals = list(pos_dict.values())
        max_scores.append(max(vals))
        min_scores.append(min(vals))
    max_log_sum = sum(math.log(max(s, 1e-10)) for s in max_scores)
    max_geo_mean = math.exp(max_log_sum / len(max_scores))
    min_log_sum = sum(math.log(max(s, 1e-10)) for s in min_scores)
    min_geo_mean = math.exp(min_log_sum / len(min_scores))

    if max_geo_mean <= min_geo_mean:
        return 0.0

    raw = (geo_mean - min_geo_mean) / (max_geo_mean - min_geo_mean)
    raw = max(0.0, min(1.0, raw))

    normalised = raw ** PSSM_CONTRAST_POWER

    return max(0.0, min(1.0, normalised))


def binding_score_to_ic50(score: float) -> float:
    """Map a binding score to an estimated IC50 (nM) using a log-linear mapping.

    Parameters
    ----------
    score : float
        Normalised binding score in [0, 1].

    Returns
    -------
    float
        Estimated IC50 in nM.

    Notes
    -----
    Effective formula: IC50 = 10 ** (IC50_LOG_INTERCEPT - IC50_LOG_SLOPE * score), calibrated so:
      - score ~0.9 -> ~50 nM (strong)
      - score ~0.5 -> ~500 nM (moderate)
      - score ~0.1 -> ~5000 nM (weak)
    """
    clamped = max(0.0, min(1.0, score))
    return 10.0 ** (IC50_LOG_INTERCEPT - IC50_LOG_SLOPE * clamped)


def classify_binding(ic50: float) -> str:
    """Classify a peptide by its IC50 value.

    Parameters
    ----------
    ic50 : float
        IC50 in nM.

    Returns
    -------
    str
        One of ``"strong_binder"``, ``"moderate_binder"``,
        ``"weak_binder"``, ``"non_binder"``.
    """
    if ic50 < IC50_STRONG_BINDER_THRESHOLD:
        return "strong_binder"
    elif ic50 <= IC50_MODERATE_BINDER_THRESHOLD:
        return "moderate_binder"
    elif ic50 <= IC50_WEAK_BINDER_THRESHOLD:
        return "weak_binder"
    else:
        return "non_binder"


def _identify_anchor_positions(
    peptide: str,
    pssm: list[dict[str, float]],
    threshold: float = 2.5,
) -> tuple[dict[int, str], dict[int, float]]:
    """Identify anchor residues in a peptide relative to a PSSM.

    An anchor position is one where the position has high selectivity
    (the ratio of max/min score in the PSSM row exceeds *threshold*).

    Returns
    -------
    anchor_residues : dict[int, str]
        Position index -> amino acid at that anchor position.
    anchor_scores : dict[int, float]
        Position index -> the PSSM score at that position.
    """
    anchor_residues: dict[int, str] = {}
    anchor_scores: dict[int, float] = {}

    for i, aa in enumerate(peptide):
        aa_upper = aa.upper()
        row = pssm[i]
        row_values = list(row.values())
        selectivity = max(row_values) / max(min(row_values), 1e-10)
        if selectivity >= threshold:
            score = row.get(aa_upper, 0.5)
            anchor_residues[i] = aa_upper
            anchor_scores[i] = score

    return anchor_residues, anchor_scores


# ═══════════════════════════════════════════════════════════════════════════
# MHC binding: prediction functions
# ═══════════════════════════════════════════════════════════════════════════


def predict_mhc_i_binding(
    protein: str,
    alleles: list[str] | None = None,
    peptide_length: int = DEFAULT_MHC_PEPTIDE_LENGTH,
    use_netmhcpan: bool = False,
    use_mhcflurry: bool = False,
) -> list[MHCBindingResult]:
    """Predict MHC class I binding for overlapping peptides.

    Parameters
    ----------
    protein : str
        Full protein amino-acid sequence.
    alleles : list[str] or None
        MHC-I alleles to evaluate. Defaults to :data:`DEFAULT_MHC_I_ALLELES`.
    peptide_length : int
        Length of peptides to extract (default 9).
    use_netmhcpan : bool
        If True, try the NetMHCpan web API first for more accurate
        predictions.  Falls back to MHCflurry or PSSM if unavailable.
        Default False.
    use_mhcflurry : bool
        If True, try MHCflurry as an offline neural-network predictor
        (AUC 0.80-0.85).  Used as intermediate tier between NetMHCpan
        (0.85-0.95) and PSSM (0.60-0.75).  Falls back to PSSM if
        MHCflurry is not installed.  Default False.

    Returns
    -------
    list[MHCBindingResult]
        Binding predictions for every peptide x allele combination.

    Prediction hierarchy
    --------------------
    NetMHCpan (online, AUC 0.85-0.95)
      → MHCflurry (offline NN, AUC 0.80-0.85)
        → PSSM (offline heuristic, AUC 0.60-0.75)
    """
    if alleles is None:
        alleles = DEFAULT_MHC_I_ALLELES

    if not protein or peptide_length < 1:
        return []

    # Try NetMHCpan if requested
    if use_netmhcpan:
        try:
            from .netmhcpan import NetMHCpanClient
            client = NetMHCpanClient()
            results = client.batch_predict(
                protein, alleles, epitope_lengths=[peptide_length],
            )
            # Convert netmhcpan MHCBindingResult objects to the
            # immunogenicity module's MHCBindingResult format
            converted = []
            for r in results:
                converted.append(MHCBindingResult(
                    allele=r.allele,
                    peptide=r.peptide,
                    start_position=r.start_position,
                    end_position=r.end_position,
                    binding_score=r.binding_score,
                    ic50_nm=r.ic50_nm,
                    binding_class=r.binding_class,
                    anchor_residues=r.anchor_residues,
                    anchor_scores=r.anchor_scores,
                ))
            logger.info(
                "MHC-I prediction via NetMHCpan: %d results for %d alleles, "
                "protein length %d",
                len(converted), len(alleles), len(protein),
            )
            return converted
        except Exception as exc:
            logger.warning(
                "NetMHCpan API failed, falling back to MHCflurry/PSSM: %s", exc,
            )

    # Try MHCflurry if requested (offline NN predictor, AUC 0.80-0.85)
    if use_mhcflurry:
        try:
            from .mhcflurry_adapter import MHCflurryClient, is_mhcflurry_available
            if is_mhcflurry_available():
                client = MHCflurryClient()
                results = client.batch_predict(
                    protein, alleles, epitope_lengths=[peptide_length],
                )
                converted = []
                for r in results:
                    converted.append(MHCBindingResult(
                        allele=r.allele,
                        peptide=r.peptide,
                        start_position=r.start_position,
                        end_position=r.end_position,
                        binding_score=r.binding_score,
                        ic50_nm=r.ic50_nm,
                        binding_class=r.binding_class,
                        anchor_residues=r.anchor_residues,
                        anchor_scores=r.anchor_scores,
                    ))
                logger.info(
                    "MHC-I prediction via MHCflurry: %d results for %d alleles, "
                    "protein length %d",
                    len(converted), len(alleles), len(protein),
                )
                return converted
        except Exception as exc:
            logger.warning(
                "MHCflurry prediction failed, falling back to PSSM: %s", exc,
            )

    # PSSM-based prediction (original implementation, also serves as fallback)
    # Check cache
    cache_key = (hashlib.sha256(protein.encode()).hexdigest(), ",".join(alleles), peptide_length)
    if cache_key in _prediction_cache:
        return _prediction_cache[cache_key]

    results: list[MHCBindingResult] = []
    mhc_i_pssms = _get_mhc_i_pssms()

    for allele in alleles:
        pssm = mhc_i_pssms.get(allele)
        if pssm is None:
            logger.debug("No PSSM for allele %s — skipping", allele)
            continue

        if len(pssm) != peptide_length:
            logger.debug(
                "PSSM length %d does not match peptide_length %d for %s — skipping",
                len(pssm),
                peptide_length,
                allele,
            )
            continue

        for start in range(len(protein) - peptide_length + 1):
            peptide = protein[start : start + peptide_length]

            if any(c.upper() not in _STANDARD_AA_SET for c in peptide):
                continue

            score = score_peptide_pssm(peptide, pssm)
            ic50 = binding_score_to_ic50(score)
            binding_class = classify_binding(ic50)
            anchor_residues, anchor_scores = _identify_anchor_positions(peptide, pssm)

            results.append(
                MHCBindingResult(
                    allele=allele,
                    peptide=peptide,
                    start_position=start,
                    end_position=start + peptide_length - 1,
                    binding_score=round(score, 6),
                    ic50_nm=round(ic50, 2),
                    binding_class=binding_class,
                    anchor_residues=anchor_residues,
                    anchor_scores={k: round(v, 4) for k, v in anchor_scores.items()},
                )
            )

    _prediction_cache[cache_key] = results

    logger.info(
        "MHC-I prediction: %d results for %d alleles, protein length %d",
        len(results),
        len(alleles),
        len(protein),
    )
    return results


def predict_mhc_ii_binding(
    protein: str,
    alleles: list[str] | None = None,
    peptide_length: int = 15,
    use_netmhcpan: bool = False,
) -> list[MHCBindingResult]:
    """Predict MHC class II binding for overlapping 15-mer peptides.

    MHC-II binding is evaluated by scanning all possible 9-mer core
    registers within each 15-mer peptide and keeping the best-scoring
    core.

    Parameters
    ----------
    protein : str
        Full protein amino-acid sequence.
    alleles : list[str] or None
        MHC-II alleles to evaluate. Defaults to :data:`DEFAULT_MHC_II_ALLELES`.
    peptide_length : int
        Length of peptides to extract (default 15).
    use_netmhcpan : bool
        If True, try the NetMHCpan web API first for more accurate
        predictions.  Falls back to PSSM if the API is unavailable or
        fails.  Default False (use PSSM heuristic only).

    Returns
    -------
    list[MHCBindingResult]
        Binding predictions for every peptide x allele combination.
    """
    if alleles is None:
        alleles = DEFAULT_MHC_II_ALLELES

    if not protein or peptide_length < 9:
        return []

    # Try NetMHCpan if requested
    if use_netmhcpan:
        try:
            from .netmhcpan import NetMHCpanClient
            client = NetMHCpanClient()
            results = client.batch_predict(
                protein, alleles, epitope_lengths=[peptide_length],
            )
            converted = []
            for r in results:
                converted.append(MHCBindingResult(
                    allele=r.allele,
                    peptide=r.peptide,
                    start_position=r.start_position,
                    end_position=r.end_position,
                    binding_score=r.binding_score,
                    ic50_nm=r.ic50_nm,
                    binding_class=r.binding_class,
                    anchor_residues=r.anchor_residues,
                    anchor_scores=r.anchor_scores,
                ))
            logger.info(
                "MHC-II prediction via NetMHCpan: %d results for %d alleles, "
                "protein length %d",
                len(converted), len(alleles), len(protein),
            )
            return converted
        except Exception as exc:
            logger.warning(
                "NetMHCpan API failed, falling back to PSSM: %s", exc,
            )

    # PSSM-based prediction (original implementation, also serves as fallback)
    # Check cache
    cache_key = (hashlib.sha256(protein.encode()).hexdigest(), ",".join(alleles), -peptide_length)
    if cache_key in _prediction_cache:
        return _prediction_cache[cache_key]

    results: list[MHCBindingResult] = []
    mhc_ii_pssms = _get_mhc_ii_pssms()

    for allele in alleles:
        pssm = mhc_ii_pssms.get(allele)
        if pssm is None:
            logger.debug("No PSSM for allele %s — skipping", allele)
            continue

        if len(pssm) != MHC_II_CORE_LENGTH:
            logger.debug(
                "PSSM length %d != core length %d for %s — skipping",
                len(pssm),
                MHC_II_CORE_LENGTH,
                allele,
            )
            continue

        for start in range(len(protein) - peptide_length + 1):
            peptide = protein[start : start + peptide_length]

            if any(c.upper() not in _STANDARD_AA_SET for c in peptide):
                continue

            best_score = 0.0
            best_core = peptide[:MHC_II_CORE_LENGTH]
            best_core_offset = 0

            for core_start in range(peptide_length - MHC_II_CORE_LENGTH + 1):
                core = peptide[core_start : core_start + MHC_II_CORE_LENGTH]
                score = score_peptide_pssm(core, pssm)
                if score > best_score:
                    best_score = score
                    best_core = core
                    best_core_offset = core_start

            ic50 = binding_score_to_ic50(best_score)
            binding_class = classify_binding(ic50)

            anchor_residues, anchor_scores = _identify_anchor_positions(
                best_core, pssm
            )

            adjusted_anchors: dict[int, str] = {
                k + best_core_offset: v for k, v in anchor_residues.items()
            }
            adjusted_scores: dict[int, float] = {
                k + best_core_offset: round(v, 4)
                for k, v in anchor_scores.items()
            }

            results.append(
                MHCBindingResult(
                    allele=allele,
                    peptide=peptide,
                    start_position=start,
                    end_position=start + peptide_length - 1,
                    binding_score=round(best_score, 6),
                    ic50_nm=round(ic50, 2),
                    binding_class=binding_class,
                    anchor_residues=adjusted_anchors,
                    anchor_scores=adjusted_scores,
                )
            )

    _prediction_cache[cache_key] = results

    logger.info(
        "MHC-II prediction: %d results for %d alleles, protein length %d",
        len(results),
        len(alleles),
        len(protein),
    )
    return results


def predict_all(
    protein: str,
    mhc_i_alleles: list[str] | None = None,
    mhc_ii_alleles: list[str] | None = None,
    use_netmhcpan: bool = False,
) -> MHCPredictionResult:
    """Run both MHC-I and MHC-II predictions and aggregate results.

    Parameters
    ----------
    protein : str
        Full protein amino-acid sequence.
    mhc_i_alleles : list[str] or None
        MHC-I alleles (defaults to :data:`DEFAULT_MHC_I_ALLELES`).
    mhc_ii_alleles : list[str] or None
        MHC-II alleles (defaults to :data:`DEFAULT_MHC_II_ALLELES`).
    use_netmhcpan : bool
        If True, try the NetMHCpan web API first for more accurate
        predictions.  Falls back to PSSM if the API is unavailable or
        fails.  Default False (use PSSM heuristic only).

    Returns
    -------
    MHCPredictionResult
        Aggregated binding prediction.
    """
    with EngineTimer() as timer:
        mhc_i_results = predict_mhc_i_binding(
            protein, alleles=mhc_i_alleles, use_netmhcpan=use_netmhcpan,
        )
        mhc_ii_results = predict_mhc_ii_binding(
            protein, alleles=mhc_ii_alleles, use_netmhcpan=use_netmhcpan,
        )

        all_results = mhc_i_results + mhc_ii_results

        strong_binders = sum(
            1 for r in all_results if r.binding_class == "strong_binder"
        )
        moderate_binders = sum(
            1 for r in all_results if r.binding_class == "moderate_binder"
        )
        weak_binders = sum(
            1 for r in all_results if r.binding_class == "weak_binder"
        )
        non_binders = sum(
            1 for r in all_results if r.binding_class == "non_binder"
        )

        binding_profile: dict[str, float] = {}
        for r in all_results:
            if r.allele not in binding_profile or r.binding_score > binding_profile[r.allele]:
                binding_profile[r.allele] = round(r.binding_score, 6)

        result = MHCPredictionResult(
            protein=protein,
            mhc_i_results=mhc_i_results,
            mhc_ii_results=mhc_ii_results,
            strong_binders=strong_binders + moderate_binders,
            weak_binders=weak_binders,
            non_binders=non_binders,
            binding_profile=binding_profile,
        )

    result.execution_time_s = round(timer.elapsed, 4)

    logger.info(
        "predict_all: %d MHC-I, %d MHC-II results; "
        "strong+moderate=%d, weak=%d, non=%d (%.2fs)",
        len(mhc_i_results),
        len(mhc_ii_results),
        result.strong_binders,
        result.weak_binders,
        result.non_binders,
        result.execution_time_s,
    )
    return result


# ═══════════════════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════════════════


def _validate_protein(protein: str) -> str:
    """Validate and normalise a protein sequence.

    Uses the shared ``validate_protein_sequence`` from engine_base, then
    raises ``ImmunogenicityError`` on failure instead of ``ValueError``.
    """
    try:
        return validate_protein_sequence(protein, "Immunogenicity")
    except ValueError as exc:
        raise ImmunogenicityError(str(exc)) from exc


def _peptide_hydrophobicity_score(peptide: str) -> float:
    """Score the hydrophobicity of a peptide core (0-1 range).

    Hydrophobic cores favour MHC binding.
    """
    if len(peptide) < 3:
        return 0.0
    core = peptide[1:-1]
    if not core:
        return 0.0
    avg_hydro = sum(HYDROPATHY.get(aa, 0.0) for aa in core) / len(core)
    normalised = (avg_hydro + HYDROPHOBICITY_OFFSET) / HYDROPHOBICITY_RANGE
    return max(0.0, min(1.0, normalised))


def _peptide_charge_score(peptide: str) -> float:
    """Score the charge balance of a peptide (0-1 range).

    A mix of charged and neutral residues favours MHC binding.
    """
    charged = sum(1 for aa in peptide if aa in "DEKRH")
    neutral = len(peptide) - charged
    if len(peptide) == 0:
        return 0.0
    ratio = min(charged, neutral) / max(charged, neutral, 1)
    return min(1.0, ratio * 1.5)


def _score_peptide_for_allele(peptide: str, allele: str) -> float:
    """Score a peptide against an MHC allele using PSSM.

    For MHC-I, the peptide must match the PSSM length.
    For MHC-II, scans all 9-mer cores within the peptide.
    """
    mhc_i = _get_mhc_i_pssms()
    mhc_ii = _get_mhc_ii_pssms()
    if allele in mhc_i:
        return score_peptide_pssm(peptide, allele)
    elif allele in mhc_ii:
        best_score = 0.0
        for offset in range(len(peptide) - 8):
            core = peptide[offset : offset + 9]
            s = score_peptide_pssm(core, allele)
            best_score = max(best_score, s)
        return best_score
    return 0.0


# ═══════════════════════════════════════════════════════════════════════════
# T-cell epitope prediction
# ═══════════════════════════════════════════════════════════════════════════


def predict_t_cell_epitopes(
    protein: str,
    mhc_alleles: list[str] | None = None,
    peptide_length: int = DEFAULT_MHC_PEPTIDE_LENGTH,
) -> list[TCellEpitopeDict]:
    """Predict T-cell epitopes in a protein sequence.

    Uses PSSM-based scoring for MHC-I (9-mers) and MHC-II
    (15-mers with 9-mer core scanning) alleles.

    Parameters
    ----------
    protein : str
        Amino-acid sequence (one-letter codes).
    mhc_alleles : list[str] | None
        MHC alleles to evaluate.  Defaults to all alleles in
        :data:`DEFAULT_MHC_I_ALLELES` + :data:`DEFAULT_MHC_II_ALLELES`.
    peptide_length : int
        Length of the sliding window for MHC-I peptides (default 9).
        MHC-II peptides are always scanned as 15-mers internally.

    Returns
    -------
    list[TCellEpitopeDict]
        Each dict contains: start, end, peptide, score, allele,
        binding_class.
    """
    protein = _validate_protein(protein)

    mhc_i_pssms = _get_mhc_i_pssms()
    mhc_ii_pssms = _get_mhc_ii_pssms()

    if mhc_alleles is not None:
        mhc_i_alleles = [a for a in mhc_alleles if a in mhc_i_pssms]
        mhc_ii_alleles = [a for a in mhc_alleles if a in mhc_ii_pssms]
        unrecognised = set(mhc_alleles) - set(mhc_i_alleles) - set(mhc_ii_alleles)
        for allele in unrecognised:
            logger.warning("Unrecognised MHC allele: %s — skipping", allele)
    else:
        mhc_i_alleles = DEFAULT_MHC_I_ALLELES
        mhc_ii_alleles = DEFAULT_MHC_II_ALLELES

    epitopes: list[TCellEpitopeDict] = []

    # MHC-I predictions
    if mhc_i_alleles:
        mhc_i_results = predict_mhc_i_binding(protein, mhc_i_alleles, peptide_length)
        for r in mhc_i_results:
            epitopes.append(TCellEpitopeDict(
                start=r.start_position,
                end=r.end_position + 1,  # exclusive
                peptide=r.peptide,
                score=round(r.binding_score, 4),
                allele=r.allele,
                binding_class=r.binding_class,
            ))

    # MHC-II predictions
    if mhc_ii_alleles:
        mhc_ii_results = predict_mhc_ii_binding(protein, mhc_ii_alleles)
        for r in mhc_ii_results:
            epitopes.append(TCellEpitopeDict(
                start=r.start_position,
                end=r.end_position + 1,  # exclusive
                peptide=r.peptide,
                score=round(r.binding_score, 4),
                allele=r.allele,
                binding_class=r.binding_class,
            ))

    epitopes.sort(key=lambda e: e["score"], reverse=True)
    return epitopes


# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: amino acid scales
# ═══════════════════════════════════════════════════════════════════════════

ANTIGENICITY_SCALE: dict[str, float] = {
    # Kolaskar-Tongaonkar antigenicity scale
    "A": 1.064, "R": 1.008, "N": 0.873, "D": 1.026,
    "C": 1.412, "E": 0.895, "Q": 1.091, "G": 0.842,
    "H": 1.105, "I": 1.142, "L": 1.170, "K": 0.933,
    "M": 1.207, "F": 1.279, "P": 0.658, "S": 0.772,
    "T": 0.789, "W": 1.190, "Y": 1.161, "V": 1.132,
}

PARKER_SCALE: dict[str, float] = {
    # Parker hydrophilicity scale (Parker et al. 1986)
    "A": -1.01, "R": 1.40, "N": 0.81, "D": 1.21,
    "C": -1.20, "E": 1.64, "Q": 0.96, "G": -0.16,
    "H": 0.56, "I": -1.42, "L": -1.42, "K": 1.73,
    "M": -1.27, "F": -1.42, "P": 0.26, "S": 0.52,
    "T": -0.19, "W": -1.07, "Y": -0.31, "V": -1.07,
}

CHOU_FASMAN_TURN: dict[str, float] = {
    # Chou-Fasman beta-turn propensity (Chou & Fasman 1974)
    "A": 0.060, "R": 0.095, "N": 0.147, "D": 0.161,
    "C": 0.108, "E": 0.056, "Q": 0.098, "G": 0.102,
    "H": 0.140, "I": 0.043, "L": 0.053, "K": 0.101,
    "M": 0.068, "F": 0.059, "P": 0.301, "S": 0.120,
    "T": 0.086, "W": 0.077, "Y": 0.114, "V": 0.050,
}

EMINI_SCALE: dict[str, float] = {
    # Emini surface probability scale (Emini et al. 1985)
    "A": 0.510, "R": 1.008, "N": 0.849, "D": 0.628,
    "C": 0.358, "E": 0.977, "Q": 0.993, "G": 0.471,
    "H": 0.873, "I": 0.296, "L": 0.332, "K": 1.027,
    "M": 0.411, "F": 0.328, "P": 0.709, "S": 0.643,
    "T": 0.549, "W": 0.307, "Y": 0.361, "V": 0.265,
}

_FLEXIBILITY_SCALE: dict[str, float] = {
    # Flexibility scale used by the BepiPred-like composite method
    "A": 0.360, "R": 0.530, "N": 0.460, "D": 0.510,
    "C": 0.350, "E": 0.500, "Q": 0.490, "G": 0.540,
    "H": 0.320, "I": 0.460, "L": 0.370, "K": 0.470,
    "M": 0.300, "F": 0.310, "P": 0.510, "S": 0.510,
    "T": 0.440, "W": 0.310, "Y": 0.420, "V": 0.390,
}

ALL_SCALES: dict[str, dict[str, float]] = {
    "kolaskar_tongaonkar": ANTIGENICITY_SCALE,
    "parker_hydrophilicity": PARKER_SCALE,
    "chou_fasman": CHOU_FASMAN_TURN,
    "eea": EMINI_SCALE,
    "bepipred_flexibility": _FLEXIBILITY_SCALE,
}

# 3-letter to 1-letter amino acid mapping (for PDB parsing)
_THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Removed in v7.5.0: ANTIGENICITY_PROPENSITY (alias for ANTIGENICITY_SCALE).

# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: data classes
# ═══════════════════════════════════════════════════════════════════════════


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


def _sliding_window_average(values: list[float], window: int) -> list[float]:
    """Compute sliding window average of a list of values."""
    n = len(values)
    if n == 0 or window <= 0:
        return []
    half = window // 2
    result: list[float] = []
    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        result.append(sum(values[start:end]) / (end - start))
    return result


def _normalize_01(scores: list[float]) -> list[float]:
    """Min-max normalize a list of scores to [0, 1]."""
    if not scores:
        return []
    min_s = min(scores)
    max_s = max(scores)
    if max_s - min_s > 1e-10:
        return [(s - min_s) / (max_s - min_s) for s in scores]
    return [0.5] * len(scores)


def _find_regions(
    norm_scores: list[float],
    protein: str,
    method: str,
    threshold: float = 0.5,
    min_length: int = 6,
    is_linear: bool = True,
    extra_props: dict | None = None,
    raw_scores: list[float] | None = None,
) -> list[EpitopeRegion]:
    """Find contiguous regions where normalized scores exceed threshold.

    Args:
        norm_scores: Per-residue scores normalized to [0, 1].
        protein: Amino acid sequence.
        method: Method name for the EpitopeRegion.
        threshold: Normalized score threshold.
        min_length: Minimum region length in residues.
        is_linear: Whether the epitope is linear.
        extra_props: Additional properties to include.
        raw_scores: Raw (un-normalized) scores for property recording.

    Returns:
        List of EpitopeRegion objects.
    """
    if not norm_scores:
        return []

    regions: list[EpitopeRegion] = []
    in_region = False
    region_start = 0

    for i, s in enumerate(norm_scores):
        if s >= threshold and not in_region:
            in_region = True
            region_start = i
        elif s < threshold and in_region:
            in_region = False
            if i - region_start >= min_length:
                _add_region(
                    regions, norm_scores, raw_scores, protein,
                    region_start, i, method, is_linear, extra_props,
                )

    if in_region:
        i = len(norm_scores)
        if i - region_start >= min_length:
            _add_region(
                regions, norm_scores, raw_scores, protein,
                region_start, i, method, is_linear, extra_props,
            )

    return regions


def _add_region(
    regions: list[EpitopeRegion],
    norm_scores: list[float],
    raw_scores: list[float] | None,
    protein: str,
    start: int,
    end: int,
    method: str,
    is_linear: bool,
    extra_props: dict | None,
) -> None:
    """Append a single EpitopeRegion to the list."""
    length = end - start
    avg_score = sum(norm_scores[start:end]) / length
    props: dict = {}
    if raw_scores is not None:
        props["raw_score_avg"] = round(sum(raw_scores[start:end]) / length, 4)
    if extra_props:
        props.update(extra_props)
    regions.append(EpitopeRegion(
        start=start,
        end=end,
        peptide=protein[start:end],
        score=round(avg_score, 4),
        method=method,
        is_linear=is_linear,
        properties=props,
    ))


# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: prediction methods
# ═══════════════════════════════════════════════════════════════════════════


def predict_kolaskar_tongaonkar(
    protein: str,
    window: int = 7,
    threshold: float = 1.0,
) -> list[EpitopeRegion]:
    """Predict B-cell epitopes using the Kolaskar-Tongaonkar antigenicity scale.

    Kolaskar & Tongaonkar (1990) developed an antigenicity scale based on
    the frequency of amino acids in known antigenic determinants. Regions
    with higher antigenicity scores are more likely to be epitopes.

    Args:
        protein: Amino acid sequence (1-letter codes).
        window: Sliding window size for averaging.
        threshold: Antigenicity threshold (default 1.0).

    Returns:
        List of EpitopeRegion predictions.
    """
    if not protein:
        return []

    protein = protein.upper()

    raw_scores = [ANTIGENICITY_SCALE.get(aa, 1.0) for aa in protein]
    smoothed = _sliding_window_average(raw_scores, window)
    norm_scores = _normalize_01(smoothed)

    min_s = min(smoothed) if smoothed else 0.0
    max_s = max(smoothed) if smoothed else 1.0
    if max_s - min_s > 1e-10:
        norm_threshold = (threshold - min_s) / (max_s - min_s)
        norm_threshold = max(0.0, min(1.0, norm_threshold))
    else:
        norm_threshold = 0.5

    return _find_regions(
        norm_scores, protein, "kolaskar_tongaonkar",
        threshold=norm_threshold,
        min_length=6,
        is_linear=True,
        extra_props={"window": window, "threshold": threshold},
        raw_scores=smoothed,
    )


def predict_parker_hydrophilicity(
    protein: str,
    window: int = 7,
) -> list[EpitopeRegion]:
    """Predict B-cell epitopes using the Parker hydrophilicity scale.

    Parker et al. (1986) showed that hydrophilic regions of proteins
    tend to be on the surface and are more likely to be B-cell epitopes.

    Args:
        protein: Amino acid sequence (1-letter codes).
        window: Sliding window size for averaging.

    Returns:
        List of EpitopeRegion predictions for hydrophilic regions.
    """
    if not protein:
        return []

    protein = protein.upper()

    raw_scores = [PARKER_SCALE.get(aa, 0.0) for aa in protein]
    smoothed = _sliding_window_average(raw_scores, window)
    norm_scores = _normalize_01(smoothed)

    return _find_regions(
        norm_scores, protein, "parker_hydrophilicity",
        threshold=0.5,
        min_length=6,
        is_linear=True,
        extra_props={"window": window, "threshold": "above_mean"},
        raw_scores=smoothed,
    )


def predict_chou_fasman_beta_turn(
    protein: str,
    window: int = 7,
) -> list[EpitopeRegion]:
    """Predict B-cell epitopes using Chou-Fasman beta-turn propensity.

    Chou & Fasman (1974) showed that beta-turns are often surface-exposed
    and correspond to B-cell epitope regions.

    Args:
        protein: Amino acid sequence (1-letter codes).
        window: Sliding window size for averaging.

    Returns:
        List of EpitopeRegion predictions for high turn-propensity regions.
    """
    if not protein:
        return []

    protein = protein.upper()

    raw_scores = [CHOU_FASMAN_TURN.get(aa, 0.0) for aa in protein]
    smoothed = _sliding_window_average(raw_scores, window)
    norm_scores = _normalize_01(smoothed)

    return _find_regions(
        norm_scores, protein, "chou_fasman",
        threshold=0.5,
        min_length=6,
        is_linear=True,
        extra_props={"window": window, "threshold": "above_mean"},
        raw_scores=smoothed,
    )


def predict_eea(
    protein: str,
) -> list[EpitopeRegion]:
    """Predict B-cell epitopes using Emini surface accessibility (EEA).

    Emini et al. (1985) developed a surface probability method based on
    the statistical analysis of surface accessibility in protein structures.

    Args:
        protein: Amino acid sequence (1-letter codes).

    Returns:
        List of EpitopeRegion predictions for surface-accessible regions.
    """
    if not protein:
        return []

    protein = protein.upper()
    n = len(protein)
    window = 6
    half = window // 2

    raw_scores: list[float] = []
    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        log_sum = 0.0
        count = 0
        for j in range(start, end):
            prob = EMINI_SCALE.get(protein[j], 0.5)
            log_sum += math.log(max(prob, 1e-10))
            count += 1
        raw_scores.append(math.exp(log_sum / count) if count > 0 else 0.5)

    norm_scores = _normalize_01(raw_scores)

    min_s = min(raw_scores) if raw_scores else 0.0
    max_s = max(raw_scores) if raw_scores else 1.0
    if max_s >= 1.0 and max_s - min_s > 1e-10:
        norm_threshold = (1.0 - min_s) / (max_s - min_s)
        norm_threshold = max(0.0, min(1.0, norm_threshold))
        raw_threshold = 1.0
    else:
        norm_threshold = 0.5
        raw_threshold = round(
            (max_s + min_s) / 2.0, 4
        ) if raw_scores else 1.0

    return _find_regions(
        norm_scores, protein, "eea",
        threshold=norm_threshold,
        min_length=6,
        is_linear=True,
        extra_props={"window": window, "threshold": raw_threshold},
        raw_scores=raw_scores,
    )


def predict_bepipred_like(
    protein: str,
    window: int = 9,
) -> list[EpitopeRegion]:
    """Predict B-cell epitopes using a simplified BepiPred-like composite method.

    Combines three properties into a composite score:
        Score = 0.4 * hydrophilicity + 0.3 * flexibility + 0.3 * surface_accessibility

    Args:
        protein: Amino acid sequence (1-letter codes).
        window: Sliding window size for smoothing.

    Returns:
        List of EpitopeRegion predictions.
    """
    if not protein:
        return []

    protein = protein.upper()

    hydro_raw = [PARKER_SCALE.get(aa, 0.0) for aa in protein]
    flex_raw = [_FLEXIBILITY_SCALE.get(aa, 0.4) for aa in protein]
    surf_raw = [EMINI_SCALE.get(aa, 0.5) for aa in protein]

    hydro_smooth = _sliding_window_average(hydro_raw, window)
    flex_smooth = _sliding_window_average(flex_raw, window)
    surf_smooth = _sliding_window_average(surf_raw, window)

    hydro_norm = _normalize_01(hydro_smooth)
    flex_norm = _normalize_01(flex_smooth)
    surf_norm = _normalize_01(surf_smooth)

    composite: list[float] = [
        0.4 * h + 0.3 * f + 0.3 * s
        for h, f, s in zip(hydro_norm, flex_norm, surf_norm)
    ]

    regions: list[EpitopeRegion] = []
    in_region = False
    region_start = 0
    min_length = 6

    for i, s in enumerate(composite):
        if s >= 0.5 and not in_region:
            in_region = True
            region_start = i
        elif s < 0.5 and in_region:
            in_region = False
            if i - region_start >= min_length:
                length = i - region_start
                avg = sum(composite[region_start:i]) / length
                h_avg = sum(hydro_norm[region_start:i]) / length
                f_avg = sum(flex_norm[region_start:i]) / length
                s_avg = sum(surf_norm[region_start:i]) / length
                regions.append(EpitopeRegion(
                    start=region_start,
                    end=i,
                    peptide=protein[region_start:i],
                    score=round(avg, 4),
                    method="bepipred",
                    is_linear=True,
                    properties={
                        "hydrophilicity_avg": round(h_avg, 4),
                        "flexibility_avg": round(f_avg, 4),
                        "surface_avg": round(s_avg, 4),
                        "window": window,
                        "weights": {
                            "hydrophilicity": 0.4,
                            "flexibility": 0.3,
                            "surface": 0.3,
                        },
                    },
                ))

    if in_region:
        i = len(composite)
        if i - region_start >= min_length:
            length = i - region_start
            avg = sum(composite[region_start:i]) / length
            h_avg = sum(hydro_norm[region_start:i]) / length
            f_avg = sum(flex_norm[region_start:i]) / length
            s_avg = sum(surf_norm[region_start:i]) / length
            regions.append(EpitopeRegion(
                start=region_start,
                end=i,
                peptide=protein[region_start:i],
                score=round(avg, 4),
                method="bepipred",
                is_linear=True,
                properties={
                    "hydrophilicity_avg": round(h_avg, 4),
                    "flexibility_avg": round(f_avg, 4),
                    "surface_avg": round(s_avg, 4),
                    "window": window,
                    "weights": {
                        "hydrophilicity": 0.4,
                        "flexibility": 0.3,
                        "surface": 0.3,
                    },
                },
            ))

    return regions


def predict_conformational_epitopes(
    pdb_string: str,
    distance_cutoff: float = 6.0,
) -> list[EpitopeRegion]:
    """Predict conformational B-cell epitopes from a PDB structure.

    Identifies surface patches on the protein structure and scores them
    by hydrophilicity, charge, and flexibility. Surface residues are
    identified by having fewer than 15 C-alpha neighbors within 12 A.
    Adjacent surface residues (within distance_cutoff in sequence) are
    clustered into patches.

    Args:
        pdb_string: PDB file content as a string.
        distance_cutoff: Maximum sequence gap to cluster adjacent surface
                        residues into patches (default 6.0).

    Returns:
        List of EpitopeRegion predictions with is_linear=False.
    """
    if not pdb_string:
        return []

    ca_atoms: list[tuple[int, float, float, float, str]] = []

    for line in pdb_string.splitlines():
        line = line.rstrip()
        if not line.startswith("ATOM"):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue

        try:
            resnum = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except (ValueError, IndexError) as exc:
            logger.debug(
                "Skipping malformed PDB ATOM line (resnum/coords): %s", exc,
            )
            continue

        resname = line[17:20].strip() if len(line) >= 20 else ""
        aa = _THREE_TO_ONE.get(resname, "X")
        ca_atoms.append((resnum, x, y, z, aa))

    if len(ca_atoms) < 3:
        logger.warning(
            "Too few C-alpha atoms (%d) in PDB for conformational epitope prediction",
            len(ca_atoms),
        )
        return []

    surface_indices: set[int] = set()
    for i, (_, xi, yi, zi, _) in enumerate(ca_atoms):
        neighbors = 0
        for j, (_, xj, yj, zj, _) in enumerate(ca_atoms):
            if i == j:
                continue
            dist_sq = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if dist_sq <= CONF_EPITOPE_NEIGHBOR_CUTOPT_ANGSTROM ** 2:
                neighbors += 1
                if neighbors >= CONF_EPITOPE_MAX_NEIGHBORS:
                    break
        if neighbors < CONF_EPITOPE_MAX_NEIGHBORS:
            surface_indices.add(i)

    if not surface_indices:
        return []

    sorted_surface = sorted(surface_indices)
    patches: list[list[int]] = []
    current_patch = [sorted_surface[0]]

    for k in range(1, len(sorted_surface)):
        if sorted_surface[k] - sorted_surface[k - 1] <= distance_cutoff:
            current_patch.append(sorted_surface[k])
        else:
            if len(current_patch) >= 3:
                patches.append(current_patch)
            current_patch = [sorted_surface[k]]

    if len(current_patch) >= 3:
        patches.append(current_patch)

    charged_aas = {"R", "K", "D", "E"}
    flexible_aas = {"G", "P"}

    epitopes: list[EpitopeRegion] = []

    for patch in patches:
        aas: list[str] = []
        hydro_sum = 0.0
        charge_count = 0
        flex_count = 0

        for idx in patch:
            _, _, _, _, aa = ca_atoms[idx]
            aas.append(aa)
            hydro_sum += PARKER_SCALE.get(aa, 0.0)
            if aa in charged_aas:
                charge_count += 1
            if aa in flexible_aas:
                flex_count += 1

        n_res = len(patch)
        hydro_avg = hydro_sum / n_res
        charge_frac = charge_count / n_res
        flex_frac = flex_count / n_res

        hydro_norm = max(0.0, min(1.0, (hydro_avg + 1.5) / 3.3))

        score = 0.4 * hydro_norm + 0.3 * charge_frac + 0.3 * flex_frac

        first_resnum = ca_atoms[patch[0]][0]
        last_resnum = ca_atoms[patch[-1]][0]
        start = first_resnum - 1
        end = last_resnum

        epitopes.append(EpitopeRegion(
            start=start,
            end=end,
            peptide="".join(aas),
            score=round(score, 4),
            method="conformational",
            is_linear=False,
            properties={
                "hydrophilicity_avg": round(hydro_avg, 4),
                "charge_fraction": round(charge_frac, 4),
                "flexibility_fraction": round(flex_frac, 4),
                "surface_residue_count": n_res,
                "pdb_residue_range": (first_resnum, last_resnum),
            },
        ))

    return epitopes


# Method dispatch
_METHOD_MAP: dict[str, object] = {
    "kolaskar_tongaonkar": predict_kolaskar_tongaonkar,
    "parker_hydrophilicity": predict_parker_hydrophilicity,
    "chou_fasman_beta_turn": predict_chou_fasman_beta_turn,
    "eea": predict_eea,
    "bepipred": predict_bepipred_like,
}


# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: combined prediction
# ═══════════════════════════════════════════════════════════════════════════


def predict_epitopes(
    protein: str,
    pdb_string: str | None = None,
    methods: list[str] | None = None,
) -> EpitopePredictionResult:
    """Run multiple B-cell epitope prediction methods and combine results.

    Default methods: kolaskar_tongaonkar, parker_hydrophilicity,
    chou_fasman_beta_turn. If a PDB structure is provided,
    conformational epitope prediction is also run.

    Per-residue scores are computed as the average across all methods
    that predict each residue as an epitope. Consensus epitope regions
    are those predicted by >= 2 methods.

    Args:
        protein: Amino acid sequence (1-letter codes).
        pdb_string: Optional PDB file content for conformational prediction.
        methods: List of method names to use. Defaults to
                 ["kolaskar_tongaonkar", "parker_hydrophilicity",
                  "chou_fasman_beta_turn"].

    Returns:
        EpitopePredictionResult with combined predictions.
    """
    if not protein:
        return EpitopePredictionResult(
            protein="",
            linear_epitopes=[],
            conformational_epitopes=[],
            per_residue_score=[],
            epitope_coverage=0.0,
            methods_used=[],
        )

    protein = protein.upper()
    n = len(protein)

    if methods is None:
        methods = [
            "kolaskar_tongaonkar",
            "parker_hydrophilicity",
            "chou_fasman_beta_turn",
        ]

    valid_methods = [m for m in methods if m in _METHOD_MAP]
    if not valid_methods:
        logger.warning("No valid methods specified for epitope prediction")
        return EpitopePredictionResult(
            protein=protein,
            linear_epitopes=[],
            conformational_epitopes=[],
            per_residue_score=[0.0] * n,
            epitope_coverage=0.0,
            methods_used=[],
        )

    all_linear_epitopes: list[EpitopeRegion] = []
    residue_method_count: list[int] = [0] * n
    residue_score_sum: list[float] = [0.0] * n

    for method_name in valid_methods:
        func = _METHOD_MAP[method_name]
        try:
            epitopes = func(protein)  # type: ignore[operator]
            all_linear_epitopes.extend(epitopes)

            for ep in epitopes:
                for pos in range(ep.start, ep.end):
                    if 0 <= pos < n:
                        residue_method_count[pos] += 1
                        residue_score_sum[pos] += ep.score
        except Exception as e:
            logger.warning("Method %s failed: %s", method_name, e)

    per_residue_score: list[float] = []
    for i in range(n):
        if residue_method_count[i] > 0:
            per_residue_score.append(
                residue_score_sum[i] / residue_method_count[i]
            )
        else:
            per_residue_score.append(0.0)

    consensus_epitopes: list[EpitopeRegion] = []
    in_region = False
    region_start = 0
    min_consensus = 2
    min_length = 6

    for i in range(n):
        if residue_method_count[i] >= min_consensus and not in_region:
            in_region = True
            region_start = i
        elif residue_method_count[i] < min_consensus and in_region:
            in_region = False
            if i - region_start >= min_length:
                length = i - region_start
                avg_score = sum(per_residue_score[region_start:i]) / length
                method_count_avg = sum(
                    residue_method_count[region_start:i]
                ) / length
                consensus_epitopes.append(EpitopeRegion(
                    start=region_start,
                    end=i,
                    peptide=protein[region_start:i],
                    score=round(avg_score, 4),
                    method="consensus",
                    is_linear=True,
                    properties={
                        "method_count_avg": round(method_count_avg, 4),
                        "contributing_methods": min_consensus,
                    },
                ))

    if in_region:
        i = n
        if i - region_start >= min_length:
            length = i - region_start
            avg_score = sum(per_residue_score[region_start:i]) / length
            method_count_avg = sum(
                residue_method_count[region_start:i]
            ) / length
            consensus_epitopes.append(EpitopeRegion(
                start=region_start,
                end=i,
                peptide=protein[region_start:i],
                score=round(avg_score, 4),
                method="consensus",
                is_linear=True,
                properties={
                    "method_count_avg": round(method_count_avg, 4),
                    "contributing_methods": min_consensus,
                },
            ))

    final_linear = list(all_linear_epitopes) + consensus_epitopes

    conformational_epitopes: list[EpitopeRegion] = []
    if pdb_string:
        try:
            conformational_epitopes = predict_conformational_epitopes(pdb_string)
        except Exception as e:
            logger.warning("Conformational epitope prediction failed: %s", e)

    epitope_residues: set[int] = set()
    for ep in all_linear_epitopes:
        for pos in range(ep.start, ep.end):
            epitope_residues.add(pos)
    for ep in consensus_epitopes:
        for pos in range(ep.start, ep.end):
            epitope_residues.add(pos)
    for ep in conformational_epitopes:
        for pos in range(ep.start, ep.end):
            epitope_residues.add(pos)

    coverage = len(epitope_residues) / n if n > 0 else 0.0

    methods_used = list(valid_methods)
    if pdb_string:
        methods_used.append("conformational")
    if consensus_epitopes:
        methods_used.append("consensus")

    return EpitopePredictionResult(
        protein=protein,
        linear_epitopes=final_linear,
        conformational_epitopes=conformational_epitopes,
        per_residue_score=per_residue_score,
        epitope_coverage=round(coverage, 4),
        methods_used=methods_used,
    )


# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: backward compatibility
# ═══════════════════════════════════════════════════════════════════════════


def compute_surface_accessibility_approx(protein: str) -> list[float]:
    """Approximate relative surface accessibility per residue.

    Based on amino-acid type and local flexibility. The result is a
    per-residue value in the range [0, 1].

    .. deprecated::
        Prefer :func:`predict_eea` for Emini surface accessibility or
        :func:`predict_epitopes` for combined B-cell predictions.

    Parameters
    ----------
    protein : str
        Amino-acid sequence.

    Returns
    -------
    list[float]
        Per-residue surface accessibility estimate.
    """
    protein = _validate_protein(protein)
    n = len(protein)
    if n == 0:
        return []

    _surface_base: dict[str, float] = {
        "A": 0.45, "C": 0.30, "D": 0.75, "E": 0.78, "F": 0.35,
        "G": 0.55, "H": 0.65, "I": 0.30, "K": 0.80, "L": 0.30,
        "M": 0.40, "N": 0.70, "P": 0.70, "Q": 0.72, "R": 0.78,
        "S": 0.65, "T": 0.60, "V": 0.30, "W": 0.40, "Y": 0.55,
    }

    _flexibility: dict[str, float] = {
        "A": 0.35, "C": 0.25, "D": 0.45, "E": 0.45, "F": 0.20,
        "G": 0.60, "H": 0.35, "I": 0.20, "K": 0.40, "L": 0.20,
        "M": 0.30, "N": 0.45, "P": 0.55, "Q": 0.40, "R": 0.40,
        "S": 0.45, "T": 0.40, "V": 0.20, "W": 0.20, "Y": 0.30,
    }

    accessibility: list[float] = []
    for i, aa in enumerate(protein):
        base = _surface_base.get(aa, 0.40)
        flex = _flexibility.get(aa, 0.30)

        win_start = max(0, i - 2)
        win_end = min(n, i + 3)
        local_flex = sum(
            _flexibility.get(protein[j], 0.30) for j in range(win_start, win_end)
        ) / (win_end - win_start)

        terminal_boost = 0.0
        if i < 3 or i >= n - 3:
            terminal_boost = 0.15 * (1.0 - min(i, n - 1 - i) / 3.0)

        combined = 0.60 * base + 0.30 * local_flex + 0.10 * flex + terminal_boost
        accessibility.append(max(0.0, min(1.0, combined)))

    return accessibility


def predict_b_cell_epitopes(
    protein: str,
    method: str = "kolaskar_tongaonkar",
) -> list[BCellEpitopeDict]:
    """Predict B-cell epitopes.

    .. deprecated::
        Prefer :func:`predict_kolaskar_tongaonkar` or
        :func:`predict_epitopes` for richer results.

    Parameters
    ----------
    protein : str
        Amino-acid sequence.
    method : str
        Prediction method (currently only "kolaskar_tongaonkar").

    Returns
    -------
    list[BCellEpitopeDict]
        Each dict: start, end, peptide, score, antigenic.
    """
    protein = _validate_protein(protein)

    if method != "kolaskar_tongaonkar":
        raise ImmunogenicityError(
            f"Unsupported B-cell epitope method: {method!r}. "
            "Only 'kolaskar_tongaonkar' is supported."
        )

    regions = predict_kolaskar_tongaonkar(protein)

    # Convert EpitopeRegion objects to dicts for backward compatibility
    result: list[BCellEpitopeDict] = []
    for r in regions:
        avg_prop = r.score
        result.append(BCellEpitopeDict(
            start=r.start,
            end=r.end,
            peptide=r.peptide,
            score=round(avg_prop, 4),
            antigenic=avg_prop >= 0.5,
        ))

    if not result:
        logger.info(
            "predict_b_cell_epitopes: no antigenic regions found above threshold "
            "for protein of length %d — returning empty list",
            len(protein),
        )

    return result


# ═══════════════════════════════════════════════════════════════════════════
# Combined immunogenicity scoring
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


def compute_immunogenicity(
    protein: str,
    mhc_alleles: list[str] | None = None,
    organism: str = "Homo_sapiens",
) -> ImmunogenicityResult:
    """Compute combined immunogenicity score for a protein.

    Runs both T-cell and B-cell epitope prediction and combines
    the results into a single score.

    Scoring formula::

        overall = 0.6 * t_cell_score + 0.4 * b_cell_score

    where:
    - t_cell_score = max epitope score (capped at 1.0)
    - b_cell_score = epitope coverage (fraction of residues in predicted epitopes)

    Classification:
    - low:      overall < 0.3
    - moderate: 0.3 <= overall < 0.6
    - high:     overall >= 0.6

    Parameters
    ----------
    protein : str
        Amino-acid sequence.
    mhc_alleles : list[str] | None
        MHC alleles for T-cell epitope prediction.

    Returns
    -------
    ImmunogenicityResult
    """
    try:
        protein = _validate_protein(protein)
    except ImmunogenicityError as exc:
        return ImmunogenicityResult(
            sequence=protein if protein else "",
            primary_score=0.0,
            classification="low",
            success=False,
            error=str(exc),
        )

    with EngineTimer() as timer:
        # T-cell prediction
        t_epitopes = predict_t_cell_epitopes(protein, mhc_alleles)
        if t_epitopes:
            t_cell_score = min(1.0, max(e["score"] for e in t_epitopes))
        else:
            t_cell_score = 0.0

        # B-cell prediction using epitope.py's predict_epitopes
        b_result = predict_epitopes(protein)
        b_epitopes_converted: list[BCellEpitopeDict] = [
            BCellEpitopeDict(
                start=ep.start,
                end=ep.end,
                peptide=ep.peptide,
                score=ep.score,
                method=ep.method,
            )
            for ep in b_result.linear_epitopes
        ]
        b_cell_score = b_result.epitope_coverage

        # Combined score
        overall_score = T_CELL_WEIGHT * t_cell_score + B_CELL_WEIGHT * b_cell_score
        overall_score = max(0.0, min(1.0, overall_score))

        # Classification
        if overall_score < IMMUNOGENICITY_LOW_THRESHOLD:
            immuno_class = "low"
        elif overall_score < IMMUNOGENICITY_HIGH_THRESHOLD:
            immuno_class = "moderate"
        else:
            immuno_class = "high"

        # Deimmunization candidates
        deimm_candidates = find_deimmunization_mutations(protein)

        return ImmunogenicityResult(
            sequence=protein,
            primary_score=round(overall_score, 4),
            classification=immuno_class,
            t_cell_score=round(t_cell_score, 4),
            b_cell_score=round(b_cell_score, 4),
            t_cell_epitopes=t_epitopes,
            b_cell_epitopes=b_epitopes_converted,
            mutations=deimm_candidates,
            execution_time_s=round(timer.elapsed, 4),
        )


# ═══════════════════════════════════════════════════════════════════════════
# Deimmunization mutation finding
# ═══════════════════════════════════════════════════════════════════════════


def find_deimmunization_mutations(
    protein: str,
    epitope_threshold: float = 0.7,
    blosum62_min: int = 0,
    organism: str = "Homo_sapiens",
) -> list[MutationResult]:
    """Find mutations that may reduce immunogenicity.

    For each T-cell epitope scoring above *epitope_threshold*,
    considers every position within the epitope and evaluates
    all 19 possible substitutions.  A substitution that reduces
    the epitope binding score and satisfies the BLOSUM62
    conservation threshold is returned as a candidate.

    Parameters
    ----------
    protein : str
        Amino-acid sequence.
    epitope_threshold : float
        Only consider epitopes with score > this threshold.
    blosum62_min : int
        Minimum BLOSUM62 substitution score (conservative mutations).

    Returns
    -------
    list[MutationResult]
        Mutation suggestions that reduce immunogenicity, sorted by
        largest binding score reduction.
    """
    protein = _validate_protein(protein)

    # Get T-cell epitopes using PSSM-based scoring
    t_epitopes = predict_t_cell_epitopes(protein)

    # Filter to strong epitopes
    strong_epitopes = [
        e for e in t_epitopes if e["score"] > epitope_threshold
    ]

    if not strong_epitopes:
        return []

    # Deduplicate: track which (position, allele) combos we've already scored
    seen: set[tuple[int, str]] = set()
    candidates: list[MutationResult] = []

    mhc_i_pssms = _get_mhc_i_pssms()
    mhc_ii_pssms = _get_mhc_ii_pssms()

    for epi in strong_epitopes:
        allele = epi["allele"]
        original_score = epi["score"]
        start = epi["start"]
        end = epi["end"]  # exclusive
        peptide = epi["peptide"]

        for pos in range(start, end):
            if pos >= len(protein):
                continue
            key = (pos, allele)
            if key in seen:
                continue
            seen.add(key)

            wildtype = protein[pos]

            for mutant in sorted(_STANDARD_AA_SET):
                if mutant == wildtype:
                    continue

                # Check BLOSUM62 conservation
                blosum_score = BLOSUM62.get(wildtype, {}).get(mutant, -10)
                if blosum_score < blosum62_min:
                    continue

                # Build mutated protein and re-score the epitope region
                mutated_protein = protein[:pos] + mutant + protein[pos + 1 :]

                # Re-score using PSSM-based scoring
                if allele in mhc_i_pssms:
                    if end <= len(mutated_protein):
                        new_peptide = mutated_protein[start:end]
                        new_score = score_peptide_pssm(new_peptide, allele)
                    else:
                        new_score = original_score
                elif allele in mhc_ii_pssms:
                    # MHC-II: 15-mer with 9-mer core scanning
                    mhc_ii_window = 15
                    pep_end = start + mhc_ii_window
                    if pep_end <= len(mutated_protein):
                        new_peptide = mutated_protein[start:pep_end]
                        new_score = _score_peptide_for_allele(new_peptide, allele)
                    else:
                        new_score = original_score
                else:
                    new_score = original_score

                score_change = new_score - original_score

                # Only keep substitutions that reduce binding
                if score_change < 0:
                    candidates.append(
                        MutationResult(
                            position=pos,
                            original=wildtype,
                            mutant=mutant,
                            delta_score=round(-score_change, 4),  # positive = improvement
                            score_type="immunogenicity",
                            engine="immunogenicity",
                            recommendation="deimmunizing",
                            description=(
                                f"{wildtype}{pos+1}{mutant}: reduces {allele} "
                                f"binding by {abs(score_change):.4f}"
                            ),
                            details={
                                "epitope": peptide,
                                "binding_score_change": round(score_change, 4),
                                "blosum62": blosum_score,
                                "protein_preserved": blosum_score >= 0,
                                "allele": allele,
                            },
                        )
                    )

    # Sort by largest improvement (highest delta_score), then by BLOSUM62
    candidates.sort(key=lambda c: (-c.delta_score, -c.details.get("blosum62", 0)))

    # Limit to top candidates
    return candidates[:MAX_DEIMMUNIZATION_CANDIDATES]


# ═══════════════════════════════════════════════════════════════════════════
# Batch API
# ═══════════════════════════════════════════════════════════════════════════


def compute_immunogenicity_batch(
    sequences: list[str],
    max_workers: int | None = None,
    **kwargs,
) -> BatchResult[ImmunogenicityResult]:
    """Compute immunogenicity scores for multiple sequences in parallel.

    Uses ``concurrent.futures.ThreadPoolExecutor`` for parallelism.

    Parameters
    ----------
    sequences : list[str]
        List of protein amino-acid sequences.
    **kwargs
        Additional keyword arguments passed to :func:`compute_immunogenicity`
        (e.g. ``mhc_alleles``).

    Returns
    -------
    BatchResult[ImmunogenicityResult]
        Batch result containing one result per input sequence, in the same order.
    """
    logger.info("compute_immunogenicity_batch: processing %d sequences", len(sequences))

    results: list[ImmunogenicityResult] = []
    errors: list[str] = []

    with EngineTimer() as batch_timer:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_idx = {
                executor.submit(compute_immunogenicity, seq, **kwargs): i
                for i, seq in enumerate(sequences)
            }
            result_map: dict[int, ImmunogenicityResult] = {}
            for future in concurrent.futures.as_completed(future_to_idx):
                idx = future_to_idx[future]
                try:
                    result_map[idx] = future.result()
                except Exception as exc:
                    logger.error(
                        "compute_immunogenicity_batch: sequence %d failed: %s",
                        idx, exc,
                    )
                    result_map[idx] = ImmunogenicityResult(
                        sequence=sequences[idx] if sequences[idx] else "",
                        primary_score=0.0,
                        classification="low",
                        success=False,
                        error=str(exc),
                    )
                    errors.append(f"sequence {idx}: {exc}")

            for i in range(len(sequences)):
                results.append(result_map[i])

    logger.info(
        "compute_immunogenicity_batch: completed %d/%d successfully",
        sum(1 for r in results if r.success),
        len(results),
    )

    return BatchResult[ImmunogenicityResult](
        results=results,
        errors=errors,
        total_time_s=round(batch_timer.elapsed, 4),
    )
