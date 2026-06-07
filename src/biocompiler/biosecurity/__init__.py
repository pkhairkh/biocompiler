"""
Biosecurity Screening Module for BioCompiler

Detects hazardous biological sequences before optimization to prevent
accidental or intentional design of harmful constructs.

Screening categories:
  - Select agent toxins (ricin, abrin, botulinum, shiga, diphtheria,
    tetanus, cholera, anthrax EF/LF)
  - Viral surface proteins (influenza HA/NA, SARS-CoV-2 spike, HIV env,
    Ebola GP)
  - Antibiotic resistance markers (blaTEM, nptII, aac(6'), cat, tetA/M/O)
  - Oncogenes and growth factors (MYC, RAS, EGFR, VEGF, BRAF, etc.)

Approach:
  - Short peptide motifs (8-12 aa) for toxin and viral protein detection
  - Nucleotide patterns (15-21 nt) for resistance marker detection
  - Risk-level classification: none, low, medium, high, critical

References:
  - CDC Select Agent Program (42 CFR Part 73)
  - Australia Group Common Control List
  - WHO Laboratory Biosafety Manual, 4th ed. (2020)
  - Nucleotide signatures from CARD (Comprehensive Antibiotic Resistance
    Database, https://card.mcmaster.ca)
"""

from __future__ import annotations

# ── Types ─────────────────────────────────────────────────────────────────────
from .types import (
    BiosecurityMode,
    BiosecurityReport,
    BiosecurityScreeningResult,
    HazardMatch,
    MatchType,
    RiskLevel,
    StrandType,
)

# ── Data tables ───────────────────────────────────────────────────────────────
from .hazard_signatures import (
    HAZARD_SIGNATURE_COUNT,
    _DNA_SIGNATURES,
    _HAZARD_SIGNATURES,
    _PROTEIN_SIGNATURES,
)
from .pathogen_signatures import (
    _MOTIF_TO_PATHOGEN,
    _PATHOGEN_SIGNATURES,
)

# ── K-mer / risk helpers ─────────────────────────────────────────────────────
from .kmer_similarity import (
    _KMER_SIZE,
    _RISK_ORDER,
    _SIMILARITY_THRESHOLD,
    _compute_kmer_similarity,
    _extract_kmers,
    _max_risk,
)

# ── Fuzzy matching ───────────────────────────────────────────────────────────
from .fuzzy_matching import (
    _COMPLEMENT,
    _fuzzy_match_edit_distance,
    _fuzzy_match_hamming,
    _hamming_distance,
    _levenshtein_distance,
    reverse_complement,
)

# ── Screening API ─────────────────────────────────────────────────────────────
from .screening import (
    _build_recommendations,
    check_biosecurity_before_optimize,
    get_biosecurity_mode,
    screen_hazardous_sequence,
    sig_risk_for_match,
)

# ── Re-export BiosecurityError from exceptions ───────────────────────────────
from ..exceptions import BiosecurityError

__all__ = [
    "BiosecurityReport",
    "BiosecurityScreeningResult",
    "HazardMatch",
    "BiosecurityError",
    "BiosecurityMode",
    "MatchType",
    "StrandType",
    "RiskLevel",
    "screen_hazardous_sequence",
    "check_biosecurity_before_optimize",
    "get_biosecurity_mode",
    "sig_risk_for_match",
    "reverse_complement",
    "_hamming_distance",
    "_levenshtein_distance",
    "_fuzzy_match_hamming",
    "_fuzzy_match_edit_distance",
    # Expose the database size for testing/validation
    "HAZARD_SIGNATURE_COUNT",
    # Legacy pathogen signature exports (used by integration tests)
    "_PATHOGEN_SIGNATURES",
    "_KMER_SIZE",
    "_SIMILARITY_THRESHOLD",
    "_extract_kmers",
    "_compute_kmer_similarity",
]
