"""MHC-I and MHC-II binding affinity prediction.

Predicts peptide-MHC binding affinity using position-specific scoring
matrices (PSSMs) derived from known binding motifs in the Immune Epitope
Database (IEDB). Provides both MHC class I (9-mer peptides) and MHC
class II (15-mer peptides with 9-mer core scanning) predictions.

This module is standalone and does not import from other biocompiler modules.
"""
from __future__ import annotations

import math
import logging
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Amino acid alphabet
# ---------------------------------------------------------------------------
AMINO_ACIDS: list[str] = list("ACDEFGHIKLMNPQRSTVWY")
_AA_INDEX: dict[str, int] = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# ---------------------------------------------------------------------------
# Default allele lists
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Population coverage of common alleles (% by ethnicity)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Helper: build a PSSM row (dict[str, float]) with defaults
# ---------------------------------------------------------------------------


def _make_pssm_row(
    preferred: dict[str, float] | None = None,
    disfavored: dict[str, float] | None = None,
    default: float = 1.0,
) -> dict[str, float]:
    """Build a single PSSM position row.

    Parameters
    ----------
    preferred : dict mapping AA → score for preferred residues
    disfavored : dict mapping AA → score for disfavored residues
    default : score for residues not mentioned in *preferred* or *disfavored*
    """
    row: dict[str, float] = {aa: default for aa in AMINO_ACIDS}
    if preferred:
        for aa, score in preferred.items():
            if aa in row:
                row[aa] = score
    if disfavored:
        for aa, score in disfavored.items():
            if aa in row:
                row[aa] = score
    return row


# ---------------------------------------------------------------------------
# MHC-I PSSMs  (9 positions x 20 AAs per allele)
# Based on known binding motifs from IEDB
# ---------------------------------------------------------------------------


def _build_mhc_i_pssms() -> dict[str, list[dict[str, float]]]:
    """Construct PSSMs for common MHC-I alleles."""

    pssms: dict[str, list[dict[str, float]]] = {}

    # --- HLA-A*02:01 ---
    # Anchor positions: P2 (L, M, I, V preferred), P9 (V, L, I preferred)
    # P1: hydrophobic preferred
    pssms["HLA-A*02:01"] = [
        # P1 - hydrophobic preferred
        _make_pssm_row(
            preferred={"L": 1.2, "M": 1.2, "I": 1.2, "V": 1.2, "A": 1.1, "F": 1.1},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.5, "R": 0.5},
            default=1.0,
        ),
        # P2 - primary anchor, strong preference for L, M, I, V
        _make_pssm_row(
            preferred={"L": 2.0, "M": 2.0, "I": 1.8, "V": 1.8},
            disfavored={"D": 0.3, "E": 0.3, "K": 0.3, "R": 0.3, "P": 0.4},
            default=0.8,
        ),
        # P3
        _make_pssm_row(
            preferred={"L": 1.1, "V": 1.1, "A": 1.1},
            disfavored={"P": 0.6},
            default=1.0,
        ),
        # P4
        _make_pssm_row(
            preferred={"K": 1.1, "R": 1.1},
            disfavored={"P": 0.6},
            default=1.0,
        ),
        # P5
        _make_pssm_row(
            preferred={"A": 1.1, "V": 1.1, "I": 1.1},
            disfavored={"P": 0.6},
            default=1.0,
        ),
        # P6
        _make_pssm_row(
            preferred={"V": 1.1, "I": 1.1, "L": 1.1},
            disfavored={"P": 0.6},
            default=1.0,
        ),
        # P7
        _make_pssm_row(
            preferred={"L": 1.1, "I": 1.1, "V": 1.1},
            disfavored={"P": 0.6},
            default=1.0,
        ),
        # P8
        _make_pssm_row(
            preferred={"A": 1.1, "V": 1.1, "L": 1.1},
            disfavored={"P": 0.6},
            default=1.0,
        ),
        # P9 - secondary anchor, V, L, I preferred
        _make_pssm_row(
            preferred={"V": 1.5, "L": 1.5, "I": 1.3, "A": 1.2},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.4, "R": 0.4, "P": 0.4},
            default=0.8,
        ),
    ]

    # --- HLA-A*01:01 ---
    # Anchor: P2 (T, S, D, E), P9 (Y, F)
    pssms["HLA-A*01:01"] = [
        # P1
        _make_pssm_row(
            preferred={"A": 1.1, "S": 1.1},
            disfavored={"W": 0.5, "R": 0.5},
            default=1.0,
        ),
        # P2 - anchor
        _make_pssm_row(
            preferred={"T": 1.8, "S": 1.6, "D": 1.5, "E": 1.5},
            disfavored={"L": 0.4, "I": 0.4, "V": 0.5, "F": 0.4},
            default=0.8,
        ),
        # P3
        _make_pssm_row(default=1.0),
        # P4
        _make_pssm_row(default=1.0),
        # P5
        _make_pssm_row(default=1.0),
        # P6
        _make_pssm_row(default=1.0),
        # P7
        _make_pssm_row(default=1.0),
        # P8
        _make_pssm_row(default=1.0),
        # P9 - anchor, Y, F preferred
        _make_pssm_row(
            preferred={"Y": 1.8, "F": 1.6},
            disfavored={"K": 0.4, "R": 0.4, "D": 0.5, "E": 0.5},
            default=0.8,
        ),
    ]

    # --- HLA-A*03:01 ---
    # Anchor: P2 (V, I, L, M), P9 (K, R, H) basic C-terminus
    pssms["HLA-A*03:01"] = [
        # P1
        _make_pssm_row(
            preferred={"A": 1.1, "S": 1.1},
            default=1.0,
        ),
        # P2 - anchor
        _make_pssm_row(
            preferred={"V": 1.8, "I": 1.8, "L": 1.6, "M": 1.6},
            disfavored={"D": 0.4, "E": 0.4, "N": 0.5, "Q": 0.5},
            default=0.8,
        ),
        # P3
        _make_pssm_row(default=1.0),
        # P4
        _make_pssm_row(default=1.0),
        # P5
        _make_pssm_row(default=1.0),
        # P6
        _make_pssm_row(default=1.0),
        # P7
        _make_pssm_row(default=1.0),
        # P8
        _make_pssm_row(default=1.0),
        # P9 - anchor, basic residues
        _make_pssm_row(
            preferred={"K": 2.0, "R": 1.8, "H": 1.4},
            disfavored={"D": 0.3, "E": 0.3, "S": 0.5, "T": 0.5},
            default=0.7,
        ),
    ]

    # --- HLA-A*24:02 ---
    # Anchor: P2 (Y, F, W), P9 (F, L, I)
    pssms["HLA-A*24:02"] = [
        # P1
        _make_pssm_row(
            preferred={"Y": 1.2, "F": 1.1},
            default=1.0,
        ),
        # P2 - primary anchor, aromatic
        _make_pssm_row(
            preferred={"Y": 2.0, "F": 2.0, "W": 1.8},
            disfavored={"D": 0.3, "E": 0.3, "K": 0.3, "R": 0.3, "P": 0.4},
            default=0.7,
        ),
        # P3
        _make_pssm_row(default=1.0),
        # P4
        _make_pssm_row(default=1.0),
        # P5
        _make_pssm_row(default=1.0),
        # P6
        _make_pssm_row(default=1.0),
        # P7
        _make_pssm_row(default=1.0),
        # P8
        _make_pssm_row(default=1.0),
        # P9 - secondary anchor, hydrophobic
        _make_pssm_row(
            preferred={"F": 1.5, "L": 1.5, "I": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.4, "R": 0.4},
            default=0.8,
        ),
    ]

    # --- HLA-B*07:02 ---
    # Anchor: P2 (P, A), P9 (L, I, V)
    pssms["HLA-B*07:02"] = [
        # P1
        _make_pssm_row(
            preferred={"A": 1.1, "P": 1.1},
            default=1.0,
        ),
        # P2 - primary anchor, Pro/Ala
        _make_pssm_row(
            preferred={"P": 2.0, "A": 1.8},
            disfavored={"D": 0.3, "E": 0.3, "K": 0.3, "R": 0.3, "W": 0.4},
            default=0.7,
        ),
        # P3
        _make_pssm_row(default=1.0),
        # P4
        _make_pssm_row(default=1.0),
        # P5
        _make_pssm_row(default=1.0),
        # P6
        _make_pssm_row(default=1.0),
        # P7
        _make_pssm_row(default=1.0),
        # P8
        _make_pssm_row(default=1.0),
        # P9 - secondary anchor, hydrophobic
        _make_pssm_row(
            preferred={"L": 1.5, "I": 1.5, "V": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.4, "R": 0.4},
            default=0.8,
        ),
    ]

    # --- HLA-B*08:01 ---
    # Anchor: P2 (K, R), P9 (L, I, V)
    pssms["HLA-B*08:01"] = [
        # P1
        _make_pssm_row(default=1.0),
        # P2 - primary anchor, basic
        _make_pssm_row(
            preferred={"K": 1.8, "R": 1.8},
            disfavored={"D": 0.3, "E": 0.3, "P": 0.4, "G": 0.5},
            default=0.8,
        ),
        # P3
        _make_pssm_row(default=1.0),
        # P4
        _make_pssm_row(default=1.0),
        # P5
        _make_pssm_row(default=1.0),
        # P6
        _make_pssm_row(default=1.0),
        # P7
        _make_pssm_row(default=1.0),
        # P8
        _make_pssm_row(default=1.0),
        # P9 - secondary anchor, hydrophobic
        _make_pssm_row(
            preferred={"L": 1.5, "I": 1.3, "V": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.5},
            default=0.8,
        ),
    ]

    return pssms


MHC_I_PSSM: dict[str, list[dict[str, float]]] = _build_mhc_i_pssms()

# ---------------------------------------------------------------------------
# MHC-II PSSMs  (9 positions x 20 AAs per allele, for core binding region)
# ---------------------------------------------------------------------------


def _build_mhc_ii_pssms() -> dict[str, list[dict[str, float]]]:
    """Construct PSSMs for common MHC-II alleles (core 9-mer)."""

    pssms: dict[str, list[dict[str, float]]] = {}

    # --- HLA-DRB1*01:01 ---
    # P1: hydrophobic (F, Y, W, L, I, V, M)
    # P4: small (A, S, T, N)
    # P6: hydrophobic (L, I, V, M, F)
    # P9: polar/charged
    pssms["HLA-DRB1*01:01"] = [
        # P1 - hydrophobic pocket
        _make_pssm_row(
            preferred={"F": 1.8, "Y": 1.7, "W": 1.6, "L": 1.5, "I": 1.4, "V": 1.4, "M": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.5, "R": 0.5},
            default=0.9,
        ),
        # P2
        _make_pssm_row(default=1.0),
        # P3
        _make_pssm_row(default=1.0),
        # P4 - small residues
        _make_pssm_row(
            preferred={"A": 1.6, "S": 1.4, "T": 1.4, "N": 1.3, "G": 1.2},
            disfavored={"W": 0.5, "F": 0.6, "Y": 0.6},
            default=0.9,
        ),
        # P5
        _make_pssm_row(default=1.0),
        # P6 - hydrophobic
        _make_pssm_row(
            preferred={"L": 1.5, "I": 1.4, "V": 1.4, "M": 1.3, "F": 1.3},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.5},
            default=0.9,
        ),
        # P7
        _make_pssm_row(default=1.0),
        # P8
        _make_pssm_row(default=1.0),
        # P9 - polar/charged
        _make_pssm_row(
            preferred={"K": 1.3, "R": 1.3, "N": 1.2, "Q": 1.2, "E": 1.1, "D": 1.1},
            default=1.0,
        ),
    ]

    # --- HLA-DRB1*04:01 ---
    # P1: aromatic (F, Y, W)
    # P4: acidic (D, E)
    # P6: small (A, S, G, N)
    # P9: hydrophobic
    pssms["HLA-DRB1*04:01"] = [
        # P1 - aromatic pocket
        _make_pssm_row(
            preferred={"F": 1.8, "Y": 1.7, "W": 1.6, "L": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.5, "R": 0.5},
            default=0.9,
        ),
        # P2
        _make_pssm_row(default=1.0),
        # P3
        _make_pssm_row(default=1.0),
        # P4 - acidic
        _make_pssm_row(
            preferred={"D": 1.8, "E": 1.6},
            disfavored={"K": 0.4, "R": 0.4, "W": 0.5},
            default=0.8,
        ),
        # P5
        _make_pssm_row(default=1.0),
        # P6 - small
        _make_pssm_row(
            preferred={"A": 1.6, "S": 1.4, "G": 1.3, "N": 1.2},
            disfavored={"W": 0.5, "F": 0.6, "Y": 0.6},
            default=0.9,
        ),
        # P7
        _make_pssm_row(default=1.0),
        # P8
        _make_pssm_row(default=1.0),
        # P9 - hydrophobic
        _make_pssm_row(
            preferred={"L": 1.4, "I": 1.3, "V": 1.3, "F": 1.3},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.5},
            default=0.9,
        ),
    ]

    # --- HLA-DRB1*07:01 ---
    # P1: hydrophobic/aromatic
    # P4: small/neutral
    # P6: hydrophobic
    # P9: polar
    pssms["HLA-DRB1*07:01"] = [
        # P1 - hydrophobic
        _make_pssm_row(
            preferred={"F": 1.6, "Y": 1.5, "L": 1.4, "I": 1.3, "V": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.5, "R": 0.5},
            default=0.9,
        ),
        # P2
        _make_pssm_row(default=1.0),
        # P3
        _make_pssm_row(default=1.0),
        # P4 - small/neutral
        _make_pssm_row(
            preferred={"A": 1.4, "S": 1.3, "T": 1.3, "N": 1.2},
            default=0.9,
        ),
        # P5
        _make_pssm_row(default=1.0),
        # P6 - hydrophobic
        _make_pssm_row(
            preferred={"L": 1.4, "I": 1.3, "V": 1.3, "F": 1.3},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.5},
            default=0.9,
        ),
        # P7
        _make_pssm_row(default=1.0),
        # P8
        _make_pssm_row(default=1.0),
        # P9 - polar
        _make_pssm_row(
            preferred={"K": 1.3, "R": 1.3, "N": 1.2, "Q": 1.2},
            default=1.0,
        ),
    ]

    return pssms


MHC_II_PSSM: dict[str, list[dict[str, float]]] = _build_mhc_ii_pssms()

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


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

    @property
    def score(self) -> float:
        """Alias for ``binding_score`` for backward compatibility."""
        return self.binding_score

    @property
    def position(self) -> int:
        """Alias for ``start_position`` for backward compatibility."""
        return self.start_position


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

    @property
    def predictions(self) -> list[MHCBindingResult]:
        """All predictions combined (MHC-I + MHC-II)."""
        return self.mhc_i_results + self.mhc_ii_results

    @property
    def class_i_results(self) -> list[MHCBindingResult]:
        """Alias for ``mhc_i_results`` for backward compatibility."""
        return self.mhc_i_results

    @property
    def class_ii_results(self) -> list[MHCBindingResult]:
        """Alias for ``mhc_ii_results`` for backward compatibility."""
        return self.mhc_ii_results

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
        # Sum coverage of alleles that have binders
        total_coverage = 0.0
        for allele, max_score in self.binding_profile.items():
            # IC50 < 500 nM corresponds to score > 0.5
            if max_score > 0.5:
                cov = POPULATION_COVERAGE.get(allele, {})
                # Use Caucasian frequency as default reference
                freq = cov.get("Caucasian", 0.0)
                total_coverage += freq / 100.0
        # Simplified: does not account for HLA linkage disequilibrium
        return min(1.0, total_coverage)

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------


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
    # Allow allele name string as shorthand
    if isinstance(pssm, str):
        allele = pssm
        lookup = MHC_I_PSSM.get(allele) or MHC_II_PSSM.get(allele)
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

    # Collect position-specific scores
    scores: list[float] = []
    for i, aa in enumerate(peptide):
        aa_upper = aa.upper()
        if aa_upper not in pssm[i]:
            # Unknown amino acid -> small penalty
            scores.append(0.3)
        else:
            scores.append(pssm[i][aa_upper])

    # Geometric mean
    log_sum = sum(math.log(max(s, 1e-10)) for s in scores)
    geo_mean = math.exp(log_sum / len(scores))

    # Normalise: theoretical max = geometric mean when every position uses
    # the highest possible score for that position.
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

    # Normalise to [0, 1] using the range [min, max]
    # A score of 0 means worst possible peptide, 1 means perfect
    raw = (geo_mean - min_geo_mean) / (max_geo_mean - min_geo_mean)
    raw = max(0.0, min(1.0, raw))

    # Power transform increases contrast between binders and non-binders
    # power=2 provides moderate contrast: random→~0.3, good→0.7, perfect→1.0
    CONTRAST_POWER = 2.0
    normalised = raw ** CONTRAST_POWER

    return max(0.0, min(1.0, normalised))


def binding_score_to_ic50(score: float) -> float:
    """Map a binding score to an estimated IC50 (nM) using a sigmoid.

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
    Uses a log-linear mapping derived from the sigmoid formulation
    IC50 = 50000 / (1 + exp(5 * (score - 0.5))) with parameters
    calibrated so that:
      - score ~0.9 → ~50 nM (strong)
      - score ~0.5 → ~500 nM (moderate)
      - score ~0.1 → ~5000 nM (weak)

    The effective formula is:
        IC50 = 10 ** (3.949 - 2.5 * score)
    """
    clamped = max(0.0, min(1.0, score))
    return 10.0 ** (3.949 - 2.5 * clamped)


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
    if ic50 < 50:
        return "strong_binder"
    elif ic50 <= 500:
        return "moderate_binder"
    elif ic50 <= 5000:
        return "weak_binder"
    else:
        return "non_binder"


# ---------------------------------------------------------------------------
# Anchor-position identification
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# MHC-I binding prediction
# ---------------------------------------------------------------------------


def predict_mhc_i_binding(
    protein: str,
    alleles: list[str] | None = None,
    peptide_length: int = 9,
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

    Returns
    -------
    list[MHCBindingResult]
        Binding predictions for every peptide x allele combination.
    """
    if alleles is None:
        alleles = DEFAULT_MHC_I_ALLELES

    if not protein or peptide_length < 1:
        return []

    results: list[MHCBindingResult] = []

    for allele in alleles:
        pssm = MHC_I_PSSM.get(allele)
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

            # Skip peptides with non-standard characters
            if any(c.upper() not in _AA_INDEX for c in peptide):
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

    logger.info(
        "MHC-I prediction: %d results for %d alleles, protein length %d",
        len(results),
        len(alleles),
        len(protein),
    )
    return results


# ---------------------------------------------------------------------------
# MHC-II binding prediction
# ---------------------------------------------------------------------------


def predict_mhc_ii_binding(
    protein: str,
    alleles: list[str] | None = None,
    peptide_length: int = 15,
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

    Returns
    -------
    list[MHCBindingResult]
        Binding predictions for every peptide x allele combination.
    """
    if alleles is None:
        alleles = DEFAULT_MHC_II_ALLELES

    if not protein or peptide_length < 9:
        return []

    core_length = 9  # MHC-II core binding region

    results: list[MHCBindingResult] = []

    for allele in alleles:
        pssm = MHC_II_PSSM.get(allele)
        if pssm is None:
            logger.debug("No PSSM for allele %s — skipping", allele)
            continue

        if len(pssm) != core_length:
            logger.debug(
                "PSSM length %d != core length %d for %s — skipping",
                len(pssm),
                core_length,
                allele,
            )
            continue

        for start in range(len(protein) - peptide_length + 1):
            peptide = protein[start : start + peptide_length]

            # Skip peptides with non-standard characters
            if any(c.upper() not in _AA_INDEX for c in peptide):
                continue

            # Scan all 9-mer cores within the 15-mer
            best_score = 0.0
            best_core = peptide[:core_length]
            best_core_offset = 0

            for core_start in range(peptide_length - core_length + 1):
                core = peptide[core_start : core_start + core_length]
                score = score_peptide_pssm(core, pssm)
                if score > best_score:
                    best_score = score
                    best_core = core
                    best_core_offset = core_start

            ic50 = binding_score_to_ic50(best_score)
            binding_class = classify_binding(ic50)

            # Identify anchor residues relative to the best core
            anchor_residues, anchor_scores = _identify_anchor_positions(
                best_core, pssm
            )

            # Adjust anchor positions to be relative to the full peptide
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

    logger.info(
        "MHC-II prediction: %d results for %d alleles, protein length %d",
        len(results),
        len(alleles),
        len(protein),
    )
    return results


# ---------------------------------------------------------------------------
# Combined prediction
# ---------------------------------------------------------------------------


def predict_all(
    protein: str,
    mhc_i_alleles: list[str] | None = None,
    mhc_ii_alleles: list[str] | None = None,
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

    Returns
    -------
    MHCPredictionResult
        Aggregated binding prediction.
    """
    mhc_i_results = predict_mhc_i_binding(protein, alleles=mhc_i_alleles)
    mhc_ii_results = predict_mhc_ii_binding(protein, alleles=mhc_ii_alleles)

    all_results = mhc_i_results + mhc_ii_results

    # Count binder classes
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

    # Build binding profile: allele -> max binding score
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

    logger.info(
        "predict_all: %d MHC-I, %d MHC-II results; "
        "strong+moderate=%d, weak=%d, non=%d",
        len(mhc_i_results),
        len(mhc_ii_results),
        result.strong_binders,
        result.weak_binders,
        result.non_binders,
    )
    return result
