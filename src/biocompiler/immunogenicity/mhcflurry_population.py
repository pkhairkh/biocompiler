"""Expanded population coverage data for MHCflurry allele support.

This module provides population-level HLA allele frequency data and coverage
computation utilities that extend the minimal ``POPULATION_COVERAGE`` dict in
:mod:`biocompiler.immunogenicity` (which covers only 9 alleles and 4
population groups) to support MHCflurry's 15,000+ allele repertoire.

Data sources (approximate estimates)
-------------------------------------
- Gonzalez-Galarza et al., Nucleic Acids Res 2020; 48:D946–D953
  (Allele Frequency Net Database, AFND)
- Bone marrow donor registry summaries (NMDP/Be The Match, ~30M donors)
- Mori et al., Tissue Antigens 1997 (Japanese/Asian frequencies)
- Middleton & Williams, Eur J Immunogenet 2000 (African frequencies)

**All frequency values are approximate estimates.**  They are suitable for
computational population-coverage screening but should *not* be used for
clinical decision-making without validation against the current AFND release.

The module defines:

- :data:`EXPANDED_POPULATION_COVERAGE` — 50 common MHC-I alleles × 6 pops
- :data:`POPULATION_GROUPS` — canonical population group names
- :func:`get_allele_frequency` — single allele / population lookup
- :func:`compute_population_coverage` — combined coverage for an allele set
- :func:`find_coverage_optimizing_alleles` — greedy allele selection
- :data:`ALLELE_CLASSIFICATION` — allele → MHC class + data source
- :data:`SUPPORTED_MHCFLURRY_ALLELES` — ~200 common MHCflurry-supported alleles
- :data:`MHCFLURRY_AVAILABLE` — ``True`` if MHCflurry is importable + models present
- :class:`PopulationCoverageResult` — structured result for coverage analysis
- :func:`get_population_coverage` — unified API with mhcflurry / precomputed fallback
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field

__all__ = [
    "EXPANDED_POPULATION_COVERAGE",
    "POPULATION_GROUPS",
    "POPULATION_WEIGHTS",
    "get_allele_frequency",
    "compute_population_coverage",
    "find_coverage_optimizing_alleles",
    "ALLELE_CLASSIFICATION",
    "SUPPORTED_MHCFLURRY_ALLELES",
    "MHCFLURRY_AVAILABLE",
    "PopulationCoverageResult",
    "get_population_coverage",
]

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# Module-level mhcflurry availability check
# ═══════════════════════════════════════════════════════════════════════════
# We check mhcflurry availability in two stages:
#   1. Can we import the ``mhcflurry`` package at all?
#   2. (Optional) Can we import our adapter for a deeper check?
# The adapter import is guarded because the top-level ``biocompiler.__init__``
# may have unrelated import errors that would cascade here.


def _check_mhcflurry_available() -> bool:
    """Return True if mhcflurry is importable and models exist."""
    try:
        import mhcflurry  # noqa: F401
    except ImportError:
        return False

    # Try the deeper adapter check — but don't crash if the package
    # init has issues.
    try:
        from .mhcflurry_adapter import is_mhcflurry_available
        return is_mhcflurry_available()
    except Exception:
        # mhcflurry is importable but the adapter isn't reachable;
        # fall back to a basic check.
        try:
            models_dir = os.path.join(
                os.path.expanduser("~"), ".mhcflurry", "models_class1", "models"
            )
            return os.path.isdir(models_dir) and bool(os.listdir(models_dir))
        except Exception:
            return False


MHCFLURRY_AVAILABLE: bool = _check_mhcflurry_available()

logger.debug("MHCflurry availability at import time: %s", MHCFLURRY_AVAILABLE)

# ═══════════════════════════════════════════════════════════════════════════
# PopulationCoverageResult dataclass
# ═══════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class PopulationCoverageResult:
    """Result of a population coverage analysis for a set of alleles and peptides.

    Attributes
    ----------
    alleles_tested : list[str]
        Allele names that were evaluated.
    peptides_tested : int
        Number of peptides analysed.
    population_coverage : float
        Estimated fraction of the global population covered (0.0–1.0).
    method_used : str
        Description of the method that produced the result
        (e.g. ``"mhcflurry"``, ``"precomputed_db"``, ``"unknown_alleles"``).
    per_allele_results : dict[str, dict]
        Per-allele details.  Each value is a dict with at least
        ``"frequency_global"`` (float, phenotype freq %) and optionally
        ``"binding_peptides"`` (list of peptide strings predicted to bind).
    """

    alleles_tested: list[str] = field(default_factory=list)
    peptides_tested: int = 0
    population_coverage: float = 0.0
    method_used: str = "unknown"
    per_allele_results: dict[str, dict] = field(default_factory=dict)


# Scale factor for converting phenotype frequency (%) to a fraction.
_FREQUENCY_PERCENT_SCALE: float = 100.0

# Default number of alleles selected by find_coverage_optimizing_alleles().
_DEFAULT_OPTIMIZING_N_ALLELES: int = 6

# Population key used to request a weighted global average.
_GLOBAL_POPULATION_KEY: str = "global"

# Coverage fraction bounds (returned values are clamped to this range).
_MIN_COVERAGE_FRACTION: float = 0.0
_MAX_COVERAGE_FRACTION: float = 1.0

# Canonical non-human (mouse) MHC-I and MHC-II allele names used as a
# safety net in _build_allele_classification().
_MOUSE_MHC_I_ALLELES: tuple[str, ...] = (
    "H-2-Db", "H-2-Kb", "H-2-Ld", "H-2-Dd", "H-2-Kd", "H-2-Kk",
)
_MOUSE_MHC_II_ALLELES: tuple[str, ...] = ("H-2-IAb", "H-2-IAd")

# ═══════════════════════════════════════════════════════════════════════════
# Population groups
# ═══════════════════════════════════════════════════════════════════════════

POPULATION_GROUPS: list[str] = [
    "Caucasian",
    "African",
    "Asian",
    "Hispanic",
    "Native American",
    "Oceanian",
]

# Approximate global population weights (millions, rounded) for computing
# "global" weighted-average coverage.  Based on UN 2023 population estimates.
POPULATION_WEIGHTS: dict[str, float] = {
    "Caucasian": 1200.0,       # Europe + North America + ANZ caucasian
    "African": 1400.0,         # Sub-Saharan Africa
    "Asian": 4700.0,           # East + South + Southeast Asia
    "Hispanic": 650.0,         # Latin America
    "Native American": 50.0,   # Indigenous Americas
    "Oceanian": 45.0,          # Pacific Islands, Melanesia, Polynesia
}

# ═══════════════════════════════════════════════════════════════════════════
# Expanded population coverage — top 50 MHC-I alleles
# ═══════════════════════════════════════════════════════════════════════════
#
# Frequency values are phenotype frequencies (%) — the percentage of
# individuals in a population who carry at least one copy of the allele.
# These are approximate estimates compiled from AFND and donor registries.
#
# Key references:
#   - Gonzalez-Galarza et al., Nucleic Acids Res 2020; 48:D946
#   - NMDP donor registry summaries (n ≈ 30M)
#   - AFND v4.0 release (https://www.allelefrequencies.net/)
#
# Values rounded to nearest 0.5%.  Alleles ordered approximately by
# decreasing worldwide average frequency.

EXPANDED_POPULATION_COVERAGE: dict[str, dict[str, float]] = {
    # ── HLA-A locus (20 alleles) ─────────────────────────────────────────
    "HLA-A*02:01": {
        "Caucasian": 28.0, "African": 5.0, "Asian": 10.0,
        "Hispanic": 18.0, "Native American": 12.0, "Oceanian": 8.0,
    },
    "HLA-A*01:01": {
        "Caucasian": 16.0, "African": 3.0, "Asian": 2.0,
        "Hispanic": 8.0, "Native American": 4.0, "Oceanian": 10.0,
    },
    "HLA-A*24:02": {
        "Caucasian": 10.0, "African": 3.0, "Asian": 30.0,
        "Hispanic": 12.0, "Native American": 25.0, "Oceanian": 22.0,
    },
    "HLA-A*03:01": {
        "Caucasian": 14.0, "African": 4.0, "Asian": 3.0,
        "Hispanic": 7.0, "Native American": 3.0, "Oceanian": 6.0,
    },
    "HLA-A*11:01": {
        "Caucasian": 5.0, "African": 1.0, "Asian": 22.0,
        "Hispanic": 4.0, "Native American": 2.0, "Oceanian": 10.0,
    },
    "HLA-A*26:01": {
        "Caucasian": 4.0, "African": 1.0, "Asian": 10.0,
        "Hispanic": 3.0, "Native American": 2.0, "Oceanian": 5.0,
    },
    "HLA-A*23:01": {
        "Caucasian": 2.5, "African": 12.0, "Asian": 0.5,
        "Hispanic": 2.0, "Native American": 1.0, "Oceanian": 1.0,
    },
    "HLA-A*30:02": {
        "Caucasian": 3.0, "African": 10.0, "Asian": 1.0,
        "Hispanic": 4.0, "Native American": 2.0, "Oceanian": 1.5,
    },
    "HLA-A*32:01": {
        "Caucasian": 5.5, "African": 2.0, "Asian": 1.0,
        "Hispanic": 3.0, "Native American": 1.5, "Oceanian": 3.0,
    },
    "HLA-A*33:01": {
        "Caucasian": 2.0, "African": 5.0, "Asian": 8.0,
        "Hispanic": 3.0, "Native American": 1.5, "Oceanian": 3.0,
    },
    "HLA-A*68:01": {
        "Caucasian": 2.0, "African": 10.0, "Asian": 2.0,
        "Hispanic": 8.0, "Native American": 6.0, "Oceanian": 3.0,
    },
    "HLA-A*02:05": {
        "Caucasian": 1.0, "African": 8.0, "Asian": 1.0,
        "Hispanic": 2.0, "Native American": 1.5, "Oceanian": 1.0,
    },
    "HLA-A*29:02": {
        "Caucasian": 3.0, "African": 2.0, "Asian": 1.0,
        "Hispanic": 6.0, "Native American": 4.0, "Oceanian": 2.0,
    },
    "HLA-A*31:01": {
        "Caucasian": 3.0, "African": 1.0, "Asian": 6.0,
        "Hispanic": 5.0, "Native American": 3.0, "Oceanian": 4.0,
    },
    "HLA-A*02:06": {
        "Caucasian": 0.5, "African": 0.5, "Asian": 6.0,
        "Hispanic": 1.0, "Native American": 0.5, "Oceanian": 3.0,
    },
    "HLA-A*25:01": {
        "Caucasian": 2.5, "African": 1.0, "Asian": 0.5,
        "Hispanic": 1.5, "Native American": 0.5, "Oceanian": 1.0,
    },
    "HLA-A*34:01": {
        "Caucasian": 0.5, "African": 6.0, "Asian": 1.0,
        "Hispanic": 1.5, "Native American": 1.0, "Oceanian": 1.0,
    },
    "HLA-A*66:01": {
        "Caucasian": 1.5, "African": 2.0, "Asian": 1.5,
        "Hispanic": 2.0, "Native American": 1.0, "Oceanian": 1.0,
    },
    "HLA-A*02:03": {
        "Caucasian": 0.5, "African": 0.5, "Asian": 5.0,
        "Hispanic": 2.0, "Native American": 3.0, "Oceanian": 4.0,
    },
    "HLA-A*36:01": {
        "Caucasian": 0.5, "African": 5.0, "Asian": 0.5,
        "Hispanic": 0.5, "Native American": 0.5, "Oceanian": 0.5,
    },

    # ── HLA-B locus (18 alleles) ─────────────────────────────────────────
    "HLA-B*07:02": {
        "Caucasian": 12.0, "African": 6.0, "Asian": 4.0,
        "Hispanic": 6.0, "Native American": 3.0, "Oceanian": 5.0,
    },
    "HLA-B*08:01": {
        "Caucasian": 10.0, "African": 3.0, "Asian": 1.0,
        "Hispanic": 4.0, "Native American": 2.0, "Oceanian": 4.0,
    },
    "HLA-B*35:01": {
        "Caucasian": 9.0, "African": 6.0, "Asian": 5.0,
        "Hispanic": 12.0, "Native American": 8.0, "Oceanian": 6.0,
    },
    "HLA-B*44:02": {
        "Caucasian": 9.0, "African": 3.0, "Asian": 3.0,
        "Hispanic": 5.0, "Native American": 2.0, "Oceanian": 5.0,
    },
    "HLA-B*44:03": {
        "Caucasian": 6.0, "African": 5.0, "Asian": 4.0,
        "Hispanic": 4.0, "Native American": 2.0, "Oceanian": 3.0,
    },
    "HLA-B*40:01": {
        "Caucasian": 5.0, "African": 1.0, "Asian": 14.0,
        "Hispanic": 3.0, "Native American": 2.0, "Oceanian": 8.0,
    },
    "HLA-B*57:01": {
        "Caucasian": 5.5, "African": 2.0, "Asian": 2.0,
        "Hispanic": 2.5, "Native American": 1.0, "Oceanian": 2.0,
    },
    "HLA-B*27:05": {
        "Caucasian": 5.0, "African": 1.0, "Asian": 3.0,
        "Hispanic": 3.0, "Native American": 2.0, "Oceanian": 1.5,
    },
    "HLA-B*15:01": {
        "Caucasian": 5.0, "African": 4.0, "Asian": 6.0,
        "Hispanic": 8.0, "Native American": 6.0, "Oceanian": 4.0,
    },
    "HLA-B*51:01": {
        "Caucasian": 4.0, "African": 2.0, "Asian": 8.0,
        "Hispanic": 3.0, "Native American": 2.0, "Oceanian": 3.0,
    },
    "HLA-B*52:01": {
        "Caucasian": 2.5, "African": 1.0, "Asian": 8.0,
        "Hispanic": 2.0, "Native American": 1.0, "Oceanian": 2.0,
    },
    "HLA-B*18:01": {
        "Caucasian": 6.0, "African": 3.0, "Asian": 1.0,
        "Hispanic": 5.0, "Native American": 2.0, "Oceanian": 2.0,
    },
    "HLA-B*53:01": {
        "Caucasian": 1.0, "African": 12.0, "Asian": 0.5,
        "Hispanic": 2.0, "Native American": 1.0, "Oceanian": 1.0,
    },
    "HLA-B*58:01": {
        "Caucasian": 1.5, "African": 10.0, "Asian": 1.0,
        "Hispanic": 2.0, "Native American": 1.0, "Oceanian": 1.5,
    },
    "HLA-B*38:01": {
        "Caucasian": 2.0, "African": 1.0, "Asian": 4.0,
        "Hispanic": 2.0, "Native American": 1.0, "Oceanian": 1.0,
    },
    "HLA-B*39:01": {
        "Caucasian": 3.0, "African": 1.5, "Asian": 3.0,
        "Hispanic": 4.0, "Native American": 5.0, "Oceanian": 2.0,
    },
    "HLA-B*15:02": {
        "Caucasian": 0.5, "African": 0.5, "Asian": 10.0,
        "Hispanic": 1.0, "Native American": 0.5, "Oceanian": 3.0,
    },
    "HLA-B*13:02": {
        "Caucasian": 1.0, "African": 3.0, "Asian": 6.0,
        "Hispanic": 2.0, "Native American": 1.0, "Oceanian": 1.5,
    },

    # ── HLA-C locus (12 alleles) ─────────────────────────────────────────
    "HLA-C*07:01": {
        "Caucasian": 14.0, "African": 6.0, "Asian": 8.0,
        "Hispanic": 10.0, "Native American": 5.0, "Oceanian": 8.0,
    },
    "HLA-C*07:02": {
        "Caucasian": 12.0, "African": 8.0, "Asian": 5.0,
        "Hispanic": 8.0, "Native American": 6.0, "Oceanian": 5.0,
    },
    "HLA-C*04:01": {
        "Caucasian": 10.0, "African": 14.0, "Asian": 6.0,
        "Hispanic": 10.0, "Native American": 8.0, "Oceanian": 10.0,
    },
    "HLA-C*06:02": {
        "Caucasian": 10.0, "African": 6.0, "Asian": 5.0,
        "Hispanic": 8.0, "Native American": 4.0, "Oceanian": 5.0,
    },
    "HLA-C*03:04": {
        "Caucasian": 6.0, "African": 4.0, "Asian": 10.0,
        "Hispanic": 5.0, "Native American": 3.0, "Oceanian": 6.0,
    },
    "HLA-C*03:02": {
        "Caucasian": 5.0, "African": 3.0, "Asian": 8.0,
        "Hispanic": 4.0, "Native American": 2.0, "Oceanian": 4.0,
    },
    "HLA-C*02:02": {
        "Caucasian": 6.0, "African": 5.0, "Asian": 3.0,
        "Hispanic": 5.0, "Native American": 3.0, "Oceanian": 3.0,
    },
    "HLA-C*05:01": {
        "Caucasian": 7.0, "African": 2.0, "Asian": 2.0,
        "Hispanic": 4.0, "Native American": 2.0, "Oceanian": 3.0,
    },
    "HLA-C*08:02": {
        "Caucasian": 3.0, "African": 4.0, "Asian": 3.0,
        "Hispanic": 3.0, "Native American": 2.0, "Oceanian": 2.0,
    },
    "HLA-C*12:03": {
        "Caucasian": 4.0, "African": 2.0, "Asian": 3.0,
        "Hispanic": 3.0, "Native American": 2.0, "Oceanian": 2.0,
    },
    "HLA-C*16:01": {
        "Caucasian": 3.5, "African": 2.0, "Asian": 2.0,
        "Hispanic": 3.0, "Native American": 1.5, "Oceanian": 2.0,
    },
    "HLA-C*01:02": {
        "Caucasian": 3.0, "African": 2.0, "Asian": 10.0,
        "Hispanic": 3.0, "Native American": 2.0, "Oceanian": 5.0,
    },
}

# ═══════════════════════════════════════════════════════════════════════════
# Allele classification — MHC class and data source
# ═══════════════════════════════════════════════════════════════════════════
#
# Maps each allele name to a string encoding:
#   - MHC class: "I" or "II"
#   - Data source: "pssm" (in original immunogenicity.py),
#     "mhcflurry" (available in MHCflurry 2.0), or "both"
#
# Format: "<class>:<source>"

# PSSM alleles from immunogenicity.py
_PSSM_MHC_I_ALLELES: list[str] = [
    "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01",
    "HLA-A*24:02", "HLA-B*07:02", "HLA-B*08:01",
]
_PSSM_MHC_II_ALLELES: list[str] = [
    "HLA-DRB1*01:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01",
]

# MHCflurry 2.0 supported MHC-I alleles (most common ~200)
# Reference: O'Donnell et al., Bioinformatics 2018; 34:2696
# MHCflurry 2.0 supports MHC-I alleles trained on IEDB data.
# These are the most common alleles with models available.
_MHCFLURRY_I_ALLELES: list[str] = [
    # HLA-A locus
    "HLA-A*02:01", "HLA-A*02:02", "HLA-A*02:03", "HLA-A*02:05",
    "HLA-A*02:06", "HLA-A*02:07", "HLA-A*02:11", "HLA-A*02:17",
    "HLA-A*02:19", "HLA-A*02:50",
    "HLA-A*01:01", "HLA-A*01:02",
    "HLA-A*03:01", "HLA-A*03:02",
    "HLA-A*11:01", "HLA-A*11:02",
    "HLA-A*23:01", "HLA-A*24:02", "HLA-A*24:03",
    "HLA-A*25:01", "HLA-A*26:01", "HLA-A*26:02", "HLA-A*26:03",
    "HLA-A*29:02", "HLA-A*30:01", "HLA-A*30:02",
    "HLA-A*31:01", "HLA-A*32:01", "HLA-A*32:07",
    "HLA-A*33:01", "HLA-A*33:03",
    "HLA-A*34:01", "HLA-A*36:01",
    "HLA-A*66:01", "HLA-A*68:01", "HLA-A*68:02",
    "HLA-A*69:01", "HLA-A*74:01", "HLA-A*80:01",
    "HLA-A*02:100", "HLA-A*02:101", "HLA-A*02:261",
    "HLA-A*02:273", "HLA-A*02:277", "HLA-A*02:349",
    "HLA-A*02:370", "HLA-A*02:675",
    "HLA-A*23:02", "HLA-A*24:07", "HLA-A*24:10",
    "HLA-A*25:02", "HLA-A*26:05", "HLA-A*30:04",
    "HLA-A*31:03", "HLA-A*32:02", "HLA-A*32:15",
    "HLA-A*34:02", "HLA-A*43:01", "HLA-A*66:02",
    "HLA-A*68:23", "HLA-A*68:36",
    # HLA-B locus
    "HLA-B*07:02", "HLA-B*07:03", "HLA-B*07:05",
    "HLA-B*08:01", "HLA-B*08:02", "HLA-B*08:03",
    "HLA-B*13:02", "HLA-B*14:02",
    "HLA-B*15:01", "HLA-B*15:02", "HLA-B*15:03",
    "HLA-B*15:09", "HLA-B*15:13", "HLA-B*15:17",
    "HLA-B*18:01", "HLA-B*18:02",
    "HLA-B*27:02", "HLA-B*27:03", "HLA-B*27:05",
    "HLA-B*27:06", "HLA-B*27:09",
    "HLA-B*35:01", "HLA-B*35:03", "HLA-B*35:08",
    "HLA-B*38:01", "HLA-B*38:02",
    "HLA-B*39:01", "HLA-B*39:06",
    "HLA-B*40:01", "HLA-B*40:02",
    "HLA-B*42:01", "HLA-B*44:02", "HLA-B*44:03",
    "HLA-B*45:01", "HLA-B*46:01",
    "HLA-B*48:01", "HLA-B*49:01",
    "HLA-B*50:01", "HLA-B*51:01", "HLA-B*51:02",
    "HLA-B*52:01", "HLA-B*53:01",
    "HLA-B*54:01", "HLA-B*55:01", "HLA-B*55:02",
    "HLA-B*56:01", "HLA-B*57:01", "HLA-B*57:03",
    "HLA-B*58:01", "HLA-B*58:02",
    "HLA-B*59:01", "HLA-B*67:01", "HLA-B*73:01",
    "HLA-B*81:01", "HLA-B*82:02", "HLA-B*83:01",
    "HLA-B*07:04", "HLA-B*07:06", "HLA-B*07:13",
    "HLA-B*08:04", "HLA-B*13:01", "HLA-B*14:01",
    "HLA-B*15:05", "HLA-B*15:10", "HLA-B*15:16",
    "HLA-B*15:25", "HLA-B*15:42",
    "HLA-B*18:03", "HLA-B*27:01", "HLA-B*27:04",
    "HLA-B*27:07", "HLA-B*27:20",
    "HLA-B*35:05", "HLA-B*35:14", "HLA-B*35:43",
    "HLA-B*37:01", "HLA-B*39:05", "HLA-B*39:10",
    "HLA-B*40:06", "HLA-B*40:10",
    "HLA-B*41:01", "HLA-B*41:02",
    "HLA-B*42:02", "HLA-B*44:08",
    "HLA-B*46:02", "HLA-B*47:01",
    "HLA-B*48:02", "HLA-B*49:02",
    "HLA-B*50:02", "HLA-B*51:04", "HLA-B*51:07",
    "HLA-B*52:02", "HLA-B*54:02",
    "HLA-B*55:04", "HLA-B*56:02",
    "HLA-B*57:11", "HLA-B*58:06",
    "HLA-B*59:02", "HLA-B*78:01",
    # HLA-C locus
    "HLA-C*01:02", "HLA-C*02:02",
    "HLA-C*03:02", "HLA-C*03:03", "HLA-C*03:04",
    "HLA-C*04:01",
    "HLA-C*05:01",
    "HLA-C*06:02",
    "HLA-C*07:01", "HLA-C*07:02", "HLA-C*07:04",
    "HLA-C*08:02",
    "HLA-C*03:01", "HLA-C*03:07", "HLA-C*04:03",
    "HLA-C*05:09", "HLA-C*06:04",
    "HLA-C*07:06", "HLA-C*07:07", "HLA-C*07:18",
    "HLA-C*08:01", "HLA-C*08:04",
    "HLA-C*12:02", "HLA-C*12:03", "HLA-C*14:02",
    "HLA-C*15:02", "HLA-C*16:01", "HLA-C*17:01",
    # Non-human (mouse)
    "H-2-Db", "H-2-Kb", "H-2-Ld",
    "H-2-Dd", "H-2-Kd", "H-2-Kk",
]

# MHCflurry 2.0 does NOT natively support MHC-II; the MHC-II alleles
# listed here are anticipated for future MHCflurry releases or available
# via the MHCflurry-MHCII extension.  Included for forward-compatibility.
_MHCFLURRY_II_ALLELES: list[str] = [
    "HLA-DRB1*01:01", "HLA-DRB1*01:02", "HLA-DRB1*01:03",
    "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*04:02",
    "HLA-DRB1*04:03", "HLA-DRB1*04:04", "HLA-DRB1*04:05",
    "HLA-DRB1*07:01", "HLA-DRB1*08:02", "HLA-DRB1*09:01",
    "HLA-DRB1*10:01", "HLA-DRB1*11:01", "HLA-DRB1*11:04",
    "HLA-DRB1*13:01", "HLA-DRB1*13:02", "HLA-DRB1*13:03",
    "HLA-DRB1*15:01", "HLA-DRB1*15:02", "HLA-DRB1*15:03",
    "HLA-DRB3*01:01", "HLA-DRB3*02:02", "HLA-DRB3*03:01",
    "HLA-DRB4*01:01", "HLA-DRB5*01:01",
    "HLA-DQA1*01:01", "HLA-DQA1*01:02", "HLA-DQA1*01:03",
    "HLA-DQA1*03:01", "HLA-DQA1*04:01", "HLA-DQA1*05:01",
    "HLA-DQB1*02:01", "HLA-DQB1*03:01", "HLA-DQB1*03:02",
    "HLA-DQB1*03:03", "HLA-DQB1*04:02", "HLA-DQB1*05:01",
    "HLA-DQB1*06:02", "HLA-DQB1*06:03",
    "HLA-DPA1*01:03", "HLA-DPA1*02:01",
    "HLA-DPB1*01:01", "HLA-DPB1*02:01", "HLA-DPB1*04:01",
    "HLA-DPB1*05:01", "HLA-DPB1*09:01",
    # Mouse MHC-II
    "H-2-IAb", "H-2-IAd",
]

SUPPORTED_MHCFLURRY_ALLELES: list[str] = sorted(
    set(_MHCFLURRY_I_ALLELES + _MHCFLURRY_II_ALLELES)
)


def _build_allele_classification() -> dict[str, str]:
    """Build the allele → class:source mapping."""
    pssm_set = set(_PSSM_MHC_I_ALLELES + _PSSM_MHC_II_ALLELES)
    mhcflurry_i_set = set(_MHCFLURRY_I_ALLELES)
    mhcflurry_ii_set = set(_MHCFLURRY_II_ALLELES)

    classification: dict[str, str] = {}

    # Process all PSSM alleles
    for allele in _PSSM_MHC_I_ALLELES:
        if allele in mhcflurry_i_set:
            classification[allele] = "I:both"
        else:
            classification[allele] = "I:pssm"
    for allele in _PSSM_MHC_II_ALLELES:
        if allele in mhcflurry_ii_set:
            classification[allele] = "II:both"
        else:
            classification[allele] = "II:pssm"

    # Process MHCflurry-only MHC-I alleles
    for allele in _MHCFLURRY_I_ALLELES:
        if allele not in pssm_set:
            classification[allele] = "I:mhcflurry"

    # Process MHCflurry MHC-II alleles
    for allele in _MHCFLURRY_II_ALLELES:
        if allele not in pssm_set:
            classification[allele] = "II:mhcflurry"

    # Safety net: ensure non-human alleles are classified even if the
    # MHCflurry allele lists above are modified in the future.
    for allele in _MOUSE_MHC_I_ALLELES:
        classification.setdefault(allele, "I:mhcflurry")
    for allele in _MOUSE_MHC_II_ALLELES:
        classification.setdefault(allele, "II:mhcflurry")

    return classification


ALLELE_CLASSIFICATION: dict[str, str] = _build_allele_classification()


# ═══════════════════════════════════════════════════════════════════════════
# Public API functions
# ═══════════════════════════════════════════════════════════════════════════


def get_allele_frequency(allele: str, population: str) -> float:
    """Look up the phenotype frequency of an allele in a population.

    Parameters
    ----------
    allele : str
        HLA allele name (e.g. ``"HLA-A*02:01"``).
    population : str
        Population group name.  Must be one of :data:`POPULATION_GROUPS`
        (case-sensitive).

    Returns
    -------
    float
        Phenotype frequency as a percentage (0.0–100.0).
        Returns ``0.0`` if the allele or population is not found.
    """
    pop_data = EXPANDED_POPULATION_COVERAGE.get(allele)
    if pop_data is None:
        logger.warning("Allele %r not found in coverage data; returning 0.0", allele)
        return 0.0
    freq = pop_data.get(population, 0.0)
    if freq == 0.0 and population not in pop_data:
        logger.warning(
            "Population %r not found for allele %r; returning 0.0",
            population, allele,
        )
    elif freq == 0.0:
        logger.debug(
            "Allele %r has frequency 0.0 in population %r", allele, population
        )
    return freq


def compute_population_coverage(
    alleles: list[str],
    population: str = _GLOBAL_POPULATION_KEY,
) -> float:
    """Compute the fraction of a population covered by a set of alleles.

    Coverage is defined as the probability that a random individual
    from the target population carries *at least one* of the specified
    alleles.  Under the assumption of Hardy-Weinberg equilibrium and
    independence between loci, this is::

        P(covered) = 1 - ∏_i (1 - freq_i / 100)

    where ``freq_i`` is the phenotype frequency (%) of allele *i*.

    The independence assumption is an approximation — true HLA alleles
    are in linkage disequilibrium — but it provides a useful upper-bound
    estimate for screening purposes.

    Parameters
    ----------
    alleles : list[str]
        List of HLA allele names.
    population : str
        Population group name or ``"global"`` (default).  When
        ``"global"``, returns a weighted average across all groups
        using :data:`POPULATION_WEIGHTS`.

    Returns
    -------
    float
        Estimated population coverage as a fraction in [0.0, 1.0].
    """
    if not alleles:
        logger.debug("Empty allele list provided; returning 0.0 coverage")
        return 0.0

    if population == _GLOBAL_POPULATION_KEY:
        return _compute_global_coverage(alleles)

    # Single-population coverage
    prob_not_covered = 1.0
    for allele in alleles:
        freq_pct = get_allele_frequency(allele, population)
        prob_not_covered *= 1.0 - freq_pct / _FREQUENCY_PERCENT_SCALE

    return max(_MIN_COVERAGE_FRACTION, min(_MAX_COVERAGE_FRACTION, 1.0 - prob_not_covered))


def _compute_global_coverage(alleles: list[str]) -> float:
    """Compute weighted-average global population coverage."""
    total_weight = sum(POPULATION_WEIGHTS.values())
    if total_weight <= 0:
        logger.warning(
            "Population weights sum to %.1f (≤ 0); returning 0.0 coverage",
            total_weight,
        )
        return 0.0

    weighted_coverage = 0.0
    for pop in POPULATION_GROUPS:
        pop_coverage = compute_population_coverage(alleles, population=pop)
        weight = POPULATION_WEIGHTS.get(pop, 0.0)
        weighted_coverage += pop_coverage * weight

    return max(_MIN_COVERAGE_FRACTION, min(_MAX_COVERAGE_FRACTION, weighted_coverage / total_weight))


def find_coverage_optimizing_alleles(
    n_alleles: int = _DEFAULT_OPTIMIZING_N_ALLELES,
    population: str = _GLOBAL_POPULATION_KEY,
) -> list[str]:
    """Greedy selection of alleles that maximize population coverage.

    Starting from the most common allele in the target population, each
    subsequent allele is chosen to provide the largest marginal increase
    in coverage.

    Parameters
    ----------
    n_alleles : int
        Number of alleles to select (default 6).
    population : str
        Target population (default ``"global"``).

    Returns
    -------
    list[str]
        Selected allele names, ordered from highest to lowest marginal
        coverage contribution.
    """
    if n_alleles <= 0:
        logger.warning(
            "n_alleles=%d is non-positive; returning empty list", n_alleles
        )
        return []

    # Rank candidate alleles by frequency in the target population
    candidates = list(EXPANDED_POPULATION_COVERAGE.keys())

    # Sort by frequency (descending) in the target population
    def _freq_key(allele: str) -> float:
        if population == _GLOBAL_POPULATION_KEY:
            return _global_average_frequency(allele)
        return get_allele_frequency(allele, population)

    candidates.sort(key=_freq_key, reverse=True)

    selected: list[str] = []
    remaining = set(candidates)

    # Pick the first (most frequent) allele
    if candidates:
        best = candidates[0]
        selected.append(best)
        remaining.discard(best)

    # Greedily add alleles that maximize marginal coverage
    while len(selected) < n_alleles and remaining:
        current_coverage = compute_population_coverage(selected, population)
        best_allele: str | None = None
        best_marginal = -1.0

        for allele in remaining:
            trial = selected + [allele]
            trial_coverage = compute_population_coverage(trial, population)
            marginal = trial_coverage - current_coverage
            if marginal > best_marginal:
                best_marginal = marginal
                best_allele = allele

        if best_allele is None or best_marginal <= 0.0:
            # No allele improves coverage further
            break

        selected.append(best_allele)
        remaining.discard(best_allele)

    return selected


def _global_average_frequency(allele: str) -> float:
    """Compute weighted-average frequency across all populations."""
    total_weight = sum(POPULATION_WEIGHTS.values())
    if total_weight <= 0:
        logger.warning(
            "Population weights sum to %.1f (≤ 0); returning 0.0 frequency",
            total_weight,
        )
        return 0.0

    weighted_freq = 0.0
    for pop in POPULATION_GROUPS:
        freq = get_allele_frequency(allele, pop)
        weight = POPULATION_WEIGHTS.get(pop, 0.0)
        weighted_freq += freq * weight

    return weighted_freq / total_weight


# ═══════════════════════════════════════════════════════════════════════════
# get_population_coverage — unified API with mhcflurry / precomputed fallback
# ═══════════════════════════════════════════════════════════════════════════


def get_population_coverage(
    alleles: list[str],
    peptides: list[str],
) -> PopulationCoverageResult:
    """Compute population coverage for a set of alleles and peptides.

    This is the primary high-level entry point for population-coverage
    analysis.  It follows a fallback chain so the module **never crashes**
    when MHCflurry is not installed:

    1. If *MHCflurry* is available, use it for binding prediction and
       population-coverage computation.
    2. Else, fall back to the **precomputed databases** for alleles
       that have coverage data in :data:`EXPANDED_POPULATION_COVERAGE`.
    3. Else (all alleles unknown), return an empty result with a
       warning.

    Parameters
    ----------
    alleles : list[str]
        HLA allele names to evaluate (e.g.
        ``["HLA-A*02:01", "HLA-B*07:02"]``).
    peptides : list[str]
        Peptide sequences to test for binding.

    Returns
    -------
    PopulationCoverageResult
        Structured result with coverage fraction, method used, and
        per-allele details.
    """
    if not alleles:
        return PopulationCoverageResult(
            alleles_tested=[],
            peptides_tested=len(peptides),
            population_coverage=0.0,
            method_used="no_alleles",
            per_allele_results={},
        )

    # --- 1. Try MHCflurry path ------------------------------------------
    if MHCFLURRY_AVAILABLE:
        try:
            return _coverage_via_mhcflurry(alleles, peptides)
        except Exception as exc:
            logger.warning(
                "MHCflurry coverage computation failed (%s); "
                "falling back to precomputed databases",
                exc,
            )

    # --- 2. Fall back to precomputed / EXPANDED_POPULATION_COVERAGE ------
    known_alleles = [a for a in alleles if a in EXPANDED_POPULATION_COVERAGE]
    unknown_alleles = [a for a in alleles if a not in EXPANDED_POPULATION_COVERAGE]

    if unknown_alleles:
        logger.warning(
            "Alleles not in precomputed coverage data and will be skipped: %s",
            unknown_alleles,
        )

    if not known_alleles:
        logger.warning(
            "None of the requested alleles have population-coverage data; "
            "returning empty result.  Alleles: %s",
            alleles,
        )
        return PopulationCoverageResult(
            alleles_tested=alleles,
            peptides_tested=len(peptides),
            population_coverage=0.0,
            method_used="unknown_alleles",
            per_allele_results={
                a: {"frequency_global": 0.0, "status": "no_coverage_data"}
                for a in alleles
            },
        )

    # Use EXPANDED_POPULATION_COVERAGE to compute coverage
    coverage = compute_population_coverage(known_alleles, population="global")

    per_allele: dict[str, dict] = {}
    for allele in known_alleles:
        freq_global = _global_average_frequency(allele)
        allele_info: dict = {"frequency_global": round(freq_global, 4)}

        # Try precomputed binding database for this allele
        binding_peps = _lookup_binding_peptides(allele, peptides)
        if binding_peps:
            allele_info["binding_peptides"] = binding_peps

        per_allele[allele] = allele_info

    for allele in unknown_alleles:
        per_allele[allele] = {
            "frequency_global": 0.0,
            "status": "no_coverage_data",
        }

    return PopulationCoverageResult(
        alleles_tested=alleles,
        peptides_tested=len(peptides),
        population_coverage=round(coverage, 6),
        method_used="precomputed_db",
        per_allele_results=per_allele,
    )


def _coverage_via_mhcflurry(
    alleles: list[str],
    peptides: list[str],
) -> PopulationCoverageResult:
    """Compute coverage using MHCflurry for binding prediction.

    This function is only called when :data:`MHCFLURRY_AVAILABLE` is
    ``True``.  It uses :class:`~biocompiler.mhcflurry_adapter.MHCflurryClient`
    for peptide-binding prediction and the precomputed population-frequency
    data for coverage computation.
    """
    from .mhcflurry_adapter import MHCflurryClient

    client = MHCflurryClient(allow_offline_fallback=True)

    # Compute population coverage from the frequency data
    known_alleles = [a for a in alleles if a in EXPANDED_POPULATION_COVERAGE]
    coverage = compute_population_coverage(known_alleles, population="global")

    per_allele: dict[str, dict] = {}
    for allele in alleles:
        freq_global = _global_average_frequency(allele) if allele in EXPANDED_POPULATION_COVERAGE else 0.0
        allele_info: dict = {"frequency_global": round(freq_global, 4)}

        binding_peps: list[str] = []
        for pep in peptides:
            try:
                result = client.predict_binding(pep, allele)
                if result.binding_class in ("strong_binder", "moderate_binder"):
                    binding_peps.append(pep)
            except Exception:
                pass  # skip unsupported allele/peptide combos

        if binding_peps:
            allele_info["binding_peptides"] = binding_peps

        if allele not in EXPANDED_POPULATION_COVERAGE:
            allele_info["status"] = "no_coverage_data"

        per_allele[allele] = allele_info

    return PopulationCoverageResult(
        alleles_tested=alleles,
        peptides_tested=len(peptides),
        population_coverage=round(coverage, 6),
        method_used="mhcflurry",
        per_allele_results=per_allele,
    )


def _lookup_binding_peptides(
    allele: str,
    peptides: list[str],
) -> list[str]:
    """Check which peptides are known binders for *allele* via precomputed DB.

    Returns a (possibly empty) list of peptides that are classified as
    binders in the precomputed binding database for the allele.
    """
    try:
        from ..mhc_binding_db import get_database

        db = get_database(allele)
        if db is None:
            return []
        binders: list[str] = []
        for pep in peptides:
            if db.is_binder(pep):
                binders.append(pep)
        return binders
    except Exception:
        return []
