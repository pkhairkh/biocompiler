"""Module-level constants for the immunogenicity subpackage.

Scoring thresholds, MHC allele defaults, population-coverage table,
precomputed binder database, and the standard-amino-acid set.

Split out of ``core.py`` (W8-a refactor).
"""
from __future__ import annotations

import concurrent.futures
import hashlib
import logging
import math
import warnings
from dataclasses import dataclass, field
from typing import List, Optional, TypedDict

from biocompiler.shared.constants import (
    BLOSUM62,
    DEFAULT_MHC_PEPTIDE_LENGTH,
    HYDROPATHY,
    STANDARD_AAS,
)
from ..engines.base import (
    BaseEngineResult,
    BatchResult,
    EngineTimer,
    MutationResult,
    validate_protein_sequence,
)
from biocompiler.shared.exceptions import ImmunogenicityError
from biocompiler.shared.types import Verdict

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
#: DEPRECATED: This constant will be removed in a future version.
#: Use MHCflurry or NetMHCpan for IC50 predictions instead.
IC50_LOG_INTERCEPT: float = 3.949

#: Log-linear IC50 mapping slope
#: DEPRECATED: This constant will be removed in a future version.
#: Use MHCflurry or NetMHCpan for IC50 predictions instead.
IC50_LOG_SLOPE: float = 2.5

#: Default PSSM score for amino acids not present in a position row
PSSM_UNKNOWN_AA_SCORE: float = 0.3

#: Power exponent applied to normalised PSSM score to increase contrast
PSSM_CONTRAST_POWER: float = 4.0

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
# Precomputed binding database
# ═══════════════════════════════════════════════════════════════════════════

#: Precomputed MHC-peptide binding affinities from IEDB.
PRECOMPUTED_BINDERS: dict[tuple[str, str], float] = {
    ("HLA-A*02:01", "GILGFVFTL"): 5.0,
    ("HLA-A*02:01", "LLFGYPVYV"): 12.0,
    ("HLA-A*02:01", "YLDVGIVLTL"): 35.0,
    ("HLA-A*02:01", "ELAGIGILTV"): 20.0,
    ("HLA-A*02:01", "IMDQVPFSV"): 45.0,
    ("HLA-A*02:01", "FLPSDFFPSV"): 50.0,
    ("HLA-A*01:01", "EADPTGHSY"): 30.0,
    ("HLA-A*01:01", "YLDVGIVLTL"): 55.0,
    ("HLA-A*03:01", "GILGFVFTL"): 200.0,
    ("HLA-A*03:01", "RLRAEAQVK"): 25.0,
    ("HLA-A*03:01", "KVLEYVIKV"): 40.0,
    ("HLA-B*07:02", "RPPIFIRRL"): 15.0,
    ("HLA-B*07:02", "IPFVSLLKP"): 60.0,
    ("HLA-B*08:01", "RAKFKQLL"): 80.0,
    ("HLA-B*08:01", "FLRGRAYGL"): 25.0,
    ("HLA-DRB1*01:01", "PKYVKQNTLKLAT"): 50.0,
    ("HLA-DRB1*04:01", "YVKQNTLKLAT"): 80.0,
    ("HLA-DRB1*07:01", "PKYVKQNTLKLAT"): 120.0,
}

__all__ = [
    "IMMUNOGENICITY_PSSM_AUC_ROC_LOW",
    "IMMUNOGENICITY_PSSM_AUC_ROC_HIGH",
    "IMMUNOGENICITY_NETMHCPAN_AUC_ROC_LOW",
    "IMMUNOGENICITY_NETMHCPAN_AUC_ROC_HIGH",
    "IMMUNOGENICITY_BCELL_AUC_ROC",
    "IC50_STRONG_BINDER_THRESHOLD",
    "IC50_MODERATE_BINDER_THRESHOLD",
    "IC50_WEAK_BINDER_THRESHOLD",
    "IC50_LOG_INTERCEPT",
    "IC50_LOG_SLOPE",
    "PSSM_UNKNOWN_AA_SCORE",
    "PSSM_CONTRAST_POWER",
    "HYDROPHOBICITY_OFFSET",
    "HYDROPHOBICITY_RANGE",
    "T_CELL_WEIGHT",
    "B_CELL_WEIGHT",
    "IMMUNOGENICITY_LOW_THRESHOLD",
    "IMMUNOGENICITY_HIGH_THRESHOLD",
    "MAX_DEIMMUNIZATION_CANDIDATES",
    "MHC_II_CORE_LENGTH",
    "CONF_EPITOPE_NEIGHBOR_CUTOPT_ANGSTROM",
    "CONF_EPITOPE_MAX_NEIGHBORS",
    "DEFAULT_MHC_I_ALLELES",
    "DEFAULT_MHC_II_ALLELES",
    "POPULATION_COVERAGE",
    "PRECOMPUTED_BINDERS",
]
