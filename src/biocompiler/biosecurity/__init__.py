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

Screening pipeline (two-pass):
  Pass 1 — Motif matching (always runs):
    - Short peptide motifs (8-12 aa) for toxin and viral protein detection
    - Nucleotide patterns (15-21 nt) for resistance marker detection
    - Fuzzy matching (Hamming distance, edit distance) for near-match detection
    - Risk-level classification: none, low, medium, high, critical
  Pass 2 — BLAST homology search (optional, when BLAST+ is available):
    - If NCBI BLAST+ is installed and a local database is configured
      (``BIOCOMPILER_BLAST_DB_PATH``), a homology search runs as a second
      pass after motif matching to catch novel variants missed by exact
      motif matching.
    - If BLAST+ is not available, this step is skipped gracefully with no
      effect on the motif-based results.

References:
  - CDC Select Agent Program (42 CFR Part 73)
  - Australia Group Common Control List
  - WHO Laboratory Biosafety Manual, 4th ed. (2020)
  - Nucleotide signatures from CARD (Comprehensive Antibiotic Resistance
    Database, https://card.mcmaster.ca)
"""

from __future__ import annotations

# ── Types ─────────────────────────────────────────────────────────────────────
from biocompiler.biosecurity.types import (
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

# ── BLAST integration ────────────────────────────────────────────────────────
from .blast_integration import check_biosecurity_blast

# ── Re-export BiosecurityError from exceptions ───────────────────────────────
from biocompiler.shared.exceptions import BiosecurityError

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
    "check_biosecurity_blast",
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
