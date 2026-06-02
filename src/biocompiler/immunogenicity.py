"""Immunogenicity scoring module for BioCompiler.

Provides T-cell and B-cell epitope prediction, combined immunogenicity
scoring, and deimmunization mutation suggestions for therapeutic protein
engineering.

All predictions are sequence-based heuristics and do not replace
experimental validation or structure-based tools such as NetMHCpan.
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Amino-acid constants
# ---------------------------------------------------------------------------

_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# ---------------------------------------------------------------------------
# Antigenicity propensity (Kolaskar-Tongaonkar, 1990)
# ---------------------------------------------------------------------------

ANTIGENICITY_PROPENSITY: dict[str, float] = {
    "C": 1.06, "W": 1.19, "M": 1.34, "H": 1.24, "Y": 0.88,
    "F": 1.31, "Q": 0.93, "L": 1.34, "I": 1.31, "P": 0.49,
    "V": 1.14, "D": 1.01, "T": 0.77, "A": 0.87, "N": 0.82,
    "G": 0.48, "S": 0.64, "E": 1.01, "K": 1.01, "R": 0.95,
}

# ---------------------------------------------------------------------------
# Kyte-Doolittle hydrophobicity (simplified lookup for scoring)
# ---------------------------------------------------------------------------

_HYDROPHOBICITY: dict[str, float] = {
    "A":  1.8, "C":  2.5, "D": -3.5, "E": -3.5, "F":  2.8,
    "G": -0.4, "H": -3.2, "I":  4.5, "K": -3.9, "L":  3.8,
    "M":  1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V":  4.2, "W": -0.9, "Y": -1.3,
}

# ---------------------------------------------------------------------------
# MHC-I position-specific preferences
# ---------------------------------------------------------------------------

MHC_I_PREFERENCES: dict[str, dict[int, dict[str, float]]] = {
    "HLA-A*02:01": {
        2: {"L": 1.0, "M": 0.9, "I": 0.85, "V": 0.8, "A": 0.3, "T": 0.2},
        9: {"V": 1.0, "L": 0.9, "I": 0.85, "A": 0.7, "T": 0.3, "S": 0.2},
    },
    "HLA-A*24:02": {
        2: {"Y": 1.0, "F": 0.9, "W": 0.8, "L": 0.3, "I": 0.2},
        9: {"F": 1.0, "L": 0.85, "I": 0.8, "Y": 0.5, "W": 0.4},
    },
    "HLA-B*07:02": {
        2: {"P": 1.0, "A": 0.8, "S": 0.3, "T": 0.2},
        9: {"L": 1.0, "I": 0.85, "V": 0.8, "A": 0.4, "M": 0.3},
    },
}

# ---------------------------------------------------------------------------
# MHC-II position-specific preferences
# ---------------------------------------------------------------------------

MHC_II_PREFERENCES: dict[str, dict[int, dict[str, float]]] = {
    "HLA-DRB1*01:01": {
        1: {"W": 1.0, "F": 0.95, "Y": 0.9, "L": 0.85, "I": 0.8, "V": 0.5, "M": 0.6},
        4: {"A": 1.0, "S": 0.9, "T": 0.85, "G": 0.5, "N": 0.4},
        6: {"W": 1.0, "F": 0.9, "Y": 0.85, "L": 0.8, "I": 0.75, "V": 0.5, "M": 0.6},
        9: {"W": 1.0, "F": 0.95, "Y": 0.9, "L": 0.85, "I": 0.8, "V": 0.5, "M": 0.6},
    },
    "HLA-DRB1*04:01": {
        1: {"Y": 1.0, "F": 0.95, "W": 0.9, "L": 0.5, "I": 0.4},
        4: {"D": 1.0, "E": 0.9, "N": 0.5, "Q": 0.4},
        6: {"S": 1.0, "T": 0.9, "A": 0.5, "G": 0.4},
        9: {"W": 1.0, "F": 0.95, "Y": 0.9, "L": 0.85, "I": 0.8, "V": 0.5, "M": 0.6},
    },
}

# ---------------------------------------------------------------------------
# BLOSUM62 substitution matrix (20x20)
# ---------------------------------------------------------------------------

BLOSUM62: dict[str, dict[str, int]] = {
    "A": {"A":  4, "R": -1, "N": -2, "D": -2, "C":  0, "Q": -1, "E": -1, "G":  0, "H": -2, "I": -1, "L": -1, "K": -1, "M":  1, "F": -2, "P": -1, "S":  1, "T":  0, "W": -3, "Y": -2, "V":  0},
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
    "M": {"A":  1, "R": -1, "N": -2, "D": -3, "C": -1, "Q":  0, "E": -2, "G": -3, "H": -2, "I":  1, "L":  2, "K": -1, "M":  5, "F":  0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V":  1},
    "F": {"A": -2, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -3, "E": -3, "G": -3, "H": -1, "I":  0, "L":  0, "K": -3, "M":  0, "F":  6, "P": -4, "S": -2, "T": -2, "W":  1, "Y":  3, "V": -1},
    "P": {"A": -1, "R": -2, "N": -2, "D": -1, "C": -3, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -4, "P":  7, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2},
    "S": {"A":  1, "R": -1, "N":  1, "D":  0, "C": -1, "Q":  0, "E":  0, "G":  0, "H": -1, "I": -2, "L": -2, "K":  0, "M": -1, "F": -2, "P": -1, "S":  4, "T":  1, "W": -3, "Y": -2, "V": -2},
    "T": {"A":  0, "R": -1, "N":  0, "D": -1, "C": -1, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S":  1, "T":  5, "W": -2, "Y": -2, "V":  0},
    "W": {"A": -3, "R": -3, "N": -4, "D": -4, "C": -2, "Q": -2, "E": -3, "G": -2, "H": -2, "I": -3, "L": -2, "K": -3, "M": -1, "F":  1, "P": -4, "S": -3, "T": -2, "W": 11, "Y":  2, "V": -3},
    "Y": {"A": -2, "R": -2, "N": -2, "D": -3, "C": -2, "Q": -1, "E": -2, "G": -3, "H":  2, "I": -1, "L": -1, "K": -2, "M": -1, "F":  3, "P": -3, "S": -2, "T": -2, "W":  2, "Y":  7, "V": -1},
    "V": {"A":  0, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -2, "E": -2, "G": -3, "H": -3, "I":  3, "L":  1, "K": -2, "M":  1, "F": -1, "P": -2, "S": -2, "T":  0, "W": -3, "Y": -1, "V":  4},
}

# ---------------------------------------------------------------------------
# Default MHC alleles
# ---------------------------------------------------------------------------

_DEFAULT_MHC_I_ALLELES = [
    "HLA-A*02:01",
    "HLA-A*24:02",
    "HLA-B*07:02",
]

_DEFAULT_MHC_II_ALLELES = [
    "HLA-DRB1*01:01",
    "HLA-DRB1*04:01",
]

_DEFAULT_MHC_ALLELES = _DEFAULT_MHC_I_ALLELES + _DEFAULT_MHC_II_ALLELES

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class ImmunogenicityResult:
    """Result of immunogenicity scoring for a protein sequence."""

    protein: str
    overall_score: float  # 0 (not immunogenic) to 1 (highly immunogenic)
    immunogenicity_class: str  # "low", "moderate", "high"
    t_cell_score: float  # T-cell epitope contribution
    b_cell_score: float  # B-cell epitope contribution
    t_cell_epitopes: list[dict]  # predicted T-cell epitopes
    b_cell_epitopes: list[dict]  # predicted B-cell epitopes
    deimmunization_candidates: list[dict]  # suggested mutations
    method: str = "sequence_based"


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _validate_protein(protein: str) -> str:
    """Validate and normalise a protein sequence."""
    protein = protein.strip().upper()
    if not protein:
        raise ValueError("Protein sequence must not be empty")
    invalid = set(protein) - _STANDARD_AA
    if invalid:
        raise ValueError(
            f"Protein contains non-standard amino acids: {sorted(invalid)}"
        )
    return protein


def _is_mhc_i_allele(allele: str) -> bool:
    """Return True if *allele* is a recognised MHC class I allele."""
    return allele in MHC_I_PREFERENCES


def _is_mhc_ii_allele(allele: str) -> bool:
    """Return True if *allele* is a recognised MHC class II allele."""
    return allele in MHC_II_PREFERENCES


def _peptide_hydrophobicity_score(peptide: str) -> float:
    """Score the hydrophobicity of a peptide core (0-1 range).

    Hydrophobic cores favour MHC binding.
    """
    if len(peptide) < 3:
        return 0.0
    # Core excludes anchor positions
    core = peptide[1:-1]
    if not core:
        return 0.0
    avg_hydro = sum(_HYDROPHOBICITY.get(aa, 0.0) for aa in core) / len(core)
    # Normalise: typical range is roughly -4.5 to +4.5
    normalised = (avg_hydro + 4.5) / 9.0
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
    # Perfect balance when charged ~ neutral
    return min(1.0, ratio * 1.5)


# ---------------------------------------------------------------------------
# T-cell epitope prediction
# ---------------------------------------------------------------------------


def _score_mhc_i_peptide(
    peptide: str,
    allele: str,
) -> float:
    """Score a 9-mer peptide against an MHC-I allele.

    Combines anchor position preferences, hydrophobicity, and charge
    balance into a single 0-1 binding score.
    """
    prefs = MHC_I_PREFERENCES.get(allele)
    if prefs is None:
        return 0.0
    if len(peptide) != 9:
        return 0.0

    # Anchor contribution (weight 0.55)
    anchor_score = 0.0
    anchor_count = 0
    for pos, aa_prefs in prefs.items():
        if 1 <= pos <= len(peptide):
            aa = peptide[pos - 1]
            anchor_score += aa_prefs.get(aa, 0.0)
            anchor_count += 1
    if anchor_count > 0:
        anchor_score /= anchor_count

    # Hydrophobicity contribution (weight 0.25)
    hydro_score = _peptide_hydrophobicity_score(peptide)

    # Charge contribution (weight 0.20)
    charge_score = _peptide_charge_score(peptide)

    raw = 0.55 * anchor_score + 0.25 * hydro_score + 0.20 * charge_score
    return max(0.0, min(1.0, raw))


def _score_mhc_ii_peptide(
    peptide: str,
    allele: str,
) -> float:
    """Score a peptide (9-mer core) against an MHC-II allele.

    MHC-II binding cores are typically 9 residues but the flanking
    residues are variable.  We score the 9-residue core starting at
    each possible offset and return the best score.
    """
    prefs = MHC_II_PREFERENCES.get(allele)
    if prefs is None:
        return 0.0

    best = 0.0
    for offset in range(len(peptide) - 8):
        core = peptide[offset : offset + 9]
        anchor_score = 0.0
        anchor_count = 0
        for pos, aa_prefs in prefs.items():
            if 1 <= pos <= len(core):
                aa = core[pos - 1]
                anchor_score += aa_prefs.get(aa, 0.0)
                anchor_count += 1
        if anchor_count > 0:
            anchor_score /= anchor_count
        hydro_score = _peptide_hydrophobicity_score(core)
        charge_score = _peptide_charge_score(core)
        raw = 0.55 * anchor_score + 0.25 * hydro_score + 0.20 * charge_score
        raw = max(0.0, min(1.0, raw))
        if raw > best:
            best = raw
    return best


def _binding_class(score: float) -> str:
    """Classify a binding score."""
    if score > 0.8:
        return "strong_binder"
    if score >= 0.5:
        return "weak_binder"
    return "non_binder"


def predict_t_cell_epitopes(
    protein: str,
    mhc_alleles: list[str] | None = None,
    peptide_length: int = 9,
) -> list[dict]:
    """Predict T-cell epitopes in a protein sequence.

    Uses position-specific scoring for MHC-I (9-mers) and MHC-II
    (15-mers with 9-mer core scanning) alleles.

    Parameters
    ----------
    protein : str
        Amino-acid sequence (one-letter codes).
    mhc_alleles : list[str] | None
        MHC alleles to evaluate.  Defaults to HLA-A*02:01,
        HLA-A*24:02, HLA-B*07:02, HLA-DRB1*01:01.
    peptide_length : int
        Length of the sliding window for MHC-I peptides (default 9).
        MHC-II peptides are always scanned as 15-mers internally.

    Returns
    -------
    list[dict]
        Each dict contains: start, end, peptide, score, allele,
        binding_class.
    """
    protein = _validate_protein(protein)
    alleles = mhc_alleles if mhc_alleles is not None else _DEFAULT_MHC_ALLELES

    epitopes: list[dict] = []

    for allele in alleles:
        if _is_mhc_i_allele(allele):
            # MHC-I: sliding window of *peptide_length* (default 9)
            for i in range(len(protein) - peptide_length + 1):
                pep = protein[i : i + peptide_length]
                score = _score_mhc_i_peptide(pep, allele)
                epitopes.append(
                    {
                        "start": i,
                        "end": i + peptide_length,
                        "peptide": pep,
                        "score": round(score, 4),
                        "allele": allele,
                        "binding_class": _binding_class(score),
                    }
                )
        elif _is_mhc_ii_allele(allele):
            # MHC-II: 15-mer sliding window with 9-mer core scoring
            mhc_ii_window = 15
            for i in range(len(protein) - mhc_ii_window + 1):
                pep = protein[i : i + mhc_ii_window]
                score = _score_mhc_ii_peptide(pep, allele)
                epitopes.append(
                    {
                        "start": i,
                        "end": i + mhc_ii_window,
                        "peptide": pep,
                        "score": round(score, 4),
                        "allele": allele,
                        "binding_class": _binding_class(score),
                    }
                )
        else:
            logger.warning("Unrecognised MHC allele: %s — skipping", allele)

    # Sort by score descending
    epitopes.sort(key=lambda e: e["score"], reverse=True)
    return epitopes


# ---------------------------------------------------------------------------
# B-cell epitope prediction
# ---------------------------------------------------------------------------


def compute_surface_accessibility_approx(protein: str) -> list[float]:
    """Approximate relative surface accessibility per residue.

    Based on amino-acid type and local flexibility.  Glycine and
    proline confer additional flexibility which increases the
    likelihood of surface exposure.  The result is a per-residue
    value in the range [0, 1].

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

    # Base surface propensity per AA type
    _surface_base: dict[str, float] = {
        "A": 0.45, "C": 0.30, "D": 0.75, "E": 0.78, "F": 0.35,
        "G": 0.55, "H": 0.65, "I": 0.30, "K": 0.80, "L": 0.30,
        "M": 0.40, "N": 0.70, "P": 0.70, "Q": 0.72, "R": 0.78,
        "S": 0.65, "T": 0.60, "V": 0.30, "W": 0.40, "Y": 0.55,
    }

    # Flexibility contribution per AA (G, P = more flexible)
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

        # Local flexibility: average of window ±2
        win_start = max(0, i - 2)
        win_end = min(n, i + 3)
        local_flex = sum(
            _flexibility.get(protein[j], 0.30) for j in range(win_start, win_end)
        ) / (win_end - win_start)

        # Terminal residues are more exposed
        terminal_boost = 0.0
        if i < 3 or i >= n - 3:
            terminal_boost = 0.15 * (1.0 - min(i, n - 1 - i) / 3.0)

        combined = 0.60 * base + 0.30 * local_flex + 0.10 * flex + terminal_boost
        accessibility.append(max(0.0, min(1.0, combined)))

    return accessibility


def predict_b_cell_epitopes(
    protein: str,
    method: str = "kolaskar_tongaonkar",
) -> list[dict]:
    """Predict B-cell epitopes using the Kolaskar-Tongaonkar method.

    Computes per-residue antigenicity propensity using a sliding
    window (7 residues) and thresholds at 1.0 to identify antigenic
    regions.

    Parameters
    ----------
    protein : str
        Amino-acid sequence.
    method : str
        Prediction method (currently only "kolaskar_tongaonkar").

    Returns
    -------
    list[dict]
        Each dict: start, end, peptide, score, antigenic.
    """
    protein = _validate_protein(protein)

    if method != "kolaskar_tongaonkar":
        raise ValueError(
            f"Unsupported B-cell epitope method: {method!r}. "
            "Only 'kolaskar_tongaonkar' is supported."
        )

    window_size = 7
    threshold = 1.0
    n = len(protein)

    if n < window_size:
        # Short protein: score the whole thing
        avg_prop = sum(ANTIGENICITY_PROPENSITY.get(aa, 0.5) for aa in protein) / n
        return [
            {
                "start": 0,
                "end": n,
                "peptide": protein,
                "score": round(avg_prop, 4),
                "antigenic": avg_prop >= threshold,
            }
        ]

    epitopes: list[dict] = []
    for i in range(n - window_size + 1):
        window = protein[i : i + window_size]
        avg_prop = sum(ANTIGENICITY_PROPENSITY.get(aa, 0.5) for aa in window) / window_size
        epitopes.append(
            {
                "start": i,
                "end": i + window_size,
                "peptide": window,
                "score": round(avg_prop, 4),
                "antigenic": avg_prop >= threshold,
            }
        )

    # Merge overlapping antigenic windows into contiguous epitopes
    merged: list[dict] = []
    antigenic_windows = [e for e in epitopes if e["antigenic"]]

    if not antigenic_windows:
        # Return the top-scoring windows even if below threshold
        epitopes.sort(key=lambda e: e["score"], reverse=True)
        return epitopes[:10]

    # Group overlapping windows
    current_start = antigenic_windows[0]["start"]
    current_end = antigenic_windows[0]["end"]
    current_scores = [antigenic_windows[0]["score"]]

    for epi in antigenic_windows[1:]:
        if epi["start"] < current_end:
            # Overlapping — extend
            current_end = max(current_end, epi["end"])
            current_scores.append(epi["score"])
        else:
            # Emit merged epitope
            avg_score = sum(current_scores) / len(current_scores)
            merged.append(
                {
                    "start": current_start,
                    "end": current_end,
                    "peptide": protein[current_start:current_end],
                    "score": round(avg_score, 4),
                    "antigenic": True,
                }
            )
            current_start = epi["start"]
            current_end = epi["end"]
            current_scores = [epi["score"]]

    # Last group
    avg_score = sum(current_scores) / len(current_scores)
    merged.append(
        {
            "start": current_start,
            "end": current_end,
            "peptide": protein[current_start:current_end],
            "score": round(avg_score, 4),
            "antigenic": True,
        }
    )

    return merged


# ---------------------------------------------------------------------------
# Combined immunogenicity scoring
# ---------------------------------------------------------------------------


def compute_immunogenicity(
    protein: str,
    mhc_alleles: list[str] | None = None,
) -> ImmunogenicityResult:
    """Compute combined immunogenicity score for a protein.

    Runs both T-cell and B-cell epitope prediction and combines
    the results into a single score.

    Scoring formula::

        overall = 0.6 * t_cell_score + 0.4 * b_cell_score

    where:
    - t_cell_score = max epitope score (capped at 1.0)
    - b_cell_score = fraction of surface residues that are antigenic

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
    protein = _validate_protein(protein)

    # T-cell prediction
    t_epitopes = predict_t_cell_epitopes(protein, mhc_alleles)
    if t_epitopes:
        t_cell_score = min(1.0, max(e["score"] for e in t_epitopes))
    else:
        t_cell_score = 0.0

    # B-cell prediction
    b_epitopes = predict_b_cell_epitopes(protein)
    surface = compute_surface_accessibility_approx(protein)

    if surface and b_epitopes:
        # Identify antigenic residues (covered by antigenic B-cell epitopes)
        antigenic_residues = set()
        for epi in b_epitopes:
            if epi.get("antigenic", False):
                for pos in range(epi["start"], epi["end"]):
                    antigenic_residues.add(pos)

        # Surface residues = those with accessibility > 0.5
        surface_residue_count = sum(1 for s in surface if s > 0.5)
        if surface_residue_count > 0:
            surface_antigenic_count = sum(
                1 for i in antigenic_residues if i < len(surface) and surface[i] > 0.5
            )
            b_cell_score = surface_antigenic_count / surface_residue_count
        else:
            b_cell_score = 0.0
    else:
        b_cell_score = 0.0

    # Combined score
    overall_score = 0.6 * t_cell_score + 0.4 * b_cell_score
    overall_score = max(0.0, min(1.0, overall_score))

    # Classification
    if overall_score < 0.3:
        immuno_class = "low"
    elif overall_score < 0.6:
        immuno_class = "moderate"
    else:
        immuno_class = "high"

    # Deimmunization candidates
    deimm_candidates = find_deimmunization_mutations(protein)

    return ImmunogenicityResult(
        protein=protein,
        overall_score=round(overall_score, 4),
        immunogenicity_class=immuno_class,
        t_cell_score=round(t_cell_score, 4),
        b_cell_score=round(b_cell_score, 4),
        t_cell_epitopes=t_epitopes,
        b_cell_epitopes=b_epitopes,
        deimmunization_candidates=deimm_candidates,
        method="sequence_based",
    )


# ---------------------------------------------------------------------------
# Deimmunization mutation finding
# ---------------------------------------------------------------------------


def find_deimmunization_mutations(
    protein: str,
    epitope_threshold: float = 0.7,
    blosum62_min: int = 0,
) -> list[dict]:
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
    list[dict]
        Each dict: position, wildtype, mutant, epitope,
        binding_score_change, blosum62, protein_preserved.
    """
    protein = _validate_protein(protein)

    # Get T-cell epitopes for all default alleles
    t_epitopes = predict_t_cell_epitopes(protein)

    # Filter to strong epitopes
    strong_epitopes = [
        e for e in t_epitopes if e["score"] > epitope_threshold
    ]

    if not strong_epitopes:
        return []

    # Deduplicate: track which (position, allele) combos we've already scored
    seen: set[tuple[int, str]] = set()
    candidates: list[dict] = []

    for epi in strong_epitopes:
        allele = epi["allele"]
        original_score = epi["score"]
        start = epi["start"]
        end = epi["end"]
        peptide = epi["peptide"]

        for pos in range(start, end):
            if pos >= len(protein):
                continue
            key = (pos, allele)
            if key in seen:
                continue
            seen.add(key)

            wildtype = protein[pos]

            for mutant in sorted(_STANDARD_AA):
                if mutant == wildtype:
                    continue

                # Check BLOSUM62 conservation
                blosum_score = BLOSUM62.get(wildtype, {}).get(mutant, -10)
                if blosum_score < blosum62_min:
                    continue

                # Build mutated protein and re-score the epitope region
                mutated_protein = protein[:pos] + mutant + protein[pos + 1 :]

                # Determine the peptide to re-score based on allele type
                if _is_mhc_i_allele(allele):
                    # MHC-I: 9-mer starting at same position
                    if end <= len(mutated_protein):
                        new_peptide = mutated_protein[start:end]
                        new_score = _score_mhc_i_peptide(new_peptide, allele)
                    else:
                        new_score = original_score
                elif _is_mhc_ii_allele(allele):
                    # MHC-II: 15-mer starting at same position
                    mhc_ii_window = 15
                    pep_end = start + mhc_ii_window
                    if pep_end <= len(mutated_protein):
                        new_peptide = mutated_protein[start:pep_end]
                        new_score = _score_mhc_ii_peptide(new_peptide, allele)
                    else:
                        new_score = original_score
                else:
                    new_score = original_score

                score_change = new_score - original_score

                # Only keep substitutions that reduce binding
                if score_change < 0:
                    candidates.append(
                        {
                            "position": pos,
                            "wildtype": wildtype,
                            "mutant": mutant,
                            "epitope": peptide,
                            "binding_score_change": round(score_change, 4),
                            "blosum62": blosum_score,
                            "protein_preserved": blosum_score >= 0,
                        }
                    )

    # Sort by largest binding score reduction, then by BLOSUM62 (most conservative first)
    candidates.sort(key=lambda c: (c["binding_score_change"], -c["blosum62"]))

    # Limit to top candidates to avoid overwhelming output
    return candidates[:200]
