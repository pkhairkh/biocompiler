"""
BioCompiler Translation Engine — Deterministic FST

FIXES from toy model:
- Multi-organism CAI support
- Handles partial codons at end of sequence
- Detailed translation metadata
- Consistent species/organism parameter naming via resolve_organism

v1.0.0 changes (BREAKING):
- CAI table unification: compute_cai() now uses CODON_ADAPTIVENESS_TABLES
  exclusively (previously the optimizer used SPECIES tables which disagreed,
  causing incorrect CAI values). CAI values will differ from legacy versions — they
  are now correct.
- Both 'species' and 'organism' parameters accepted and resolved via
  resolve_organism(). 'species' emits a DeprecationWarning when used.
"""

from __future__ import annotations

import logging
import math
import warnings
from typing import TypedDict

from biocompiler.shared.constants import CODON_TABLE, START_CODON, STOP_CODONS, reverse_complement
from biocompiler.sequence.scanner import validate_dna_sequence
from ..organisms import (
    CODON_ADAPTIVENESS_TABLES, SUPPORTED_ORGANISMS,
    resolve_organism,
)
from biocompiler.shared.exceptions import UnsupportedOrganismError


class PartialCodonError(ValueError):
    """Raised when a DNA sequence length is not a multiple of 3.

    This is a specific, catchable exception for partial codons at the
    end of a sequence, raised by :func:`translate_with_confidence`.
    It subclasses :class:`ValueError` so that existing ``except ValueError``
    handlers continue to work, but allows callers to distinguish partial-codon
    errors from other value errors.

    Attributes:
        sequence_length: Length of the input sequence.
        remainder: Number of trailing bases that do not form a complete codon.
    """

    def __init__(self, sequence_length: int, remainder: int):
        self.sequence_length = sequence_length
        self.remainder = remainder
        super().__init__(
            f"Sequence length {sequence_length} is not a multiple of 3; "
            f"{remainder} trailing base(s) do not form a complete codon. "
            f"Trim or pad the sequence before calling translate_with_confidence()."
        )

# ── NUMBA integration ──────────────────────────────────────────────
try:
    from biocompiler.optimizer.numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        compute_cai_kernel as _numba_cai_kernel,
        seq_to_bytes as _seq_to_bytes,
    )
except ImportError:
    _HAS_NUMBA = False

HAS_NUMBA: bool = _HAS_NUMBA

__all__ = [
    "translate",
    "translate_with_confidence",
    "compute_cai",
    "find_orfs",
    "ORFResult",
    "DEFAULT_MIN_ORF_LENGTH_AA",
    "BACTERIAL_START_CODONS",
    "STANDARD_START_CODON",
    "PartialCodonError",
]

logger = logging.getLogger(__name__)

# --- Named constants (replacing magic numbers) ---

# Number of nucleotides per codon
_CODON_LENGTH: int = 3

# Default minimum ORF length in amino acids for find_orfs()
DEFAULT_MIN_ORF_LENGTH_AA: int = 30

# Decimal places for rounding CAI values
_CAI_ROUND_PRECISION: int = 4

# Epsilon floor for zero-adaptiveness codons in CAI computation
_ZERO_ADAPTIVENESS_EPSILON: float = 1e-10

# ── Start codon constants ──────────────────────────────────────────
# Standard (canonical) start codon used in all organisms
STANDARD_START_CODON: str = "ATG"

# Bacterial start codons per NCBI translation table 1:
#   ATG — canonical start (Met), used by ~83% of E. coli genes
#   GTG — alternative start (codes Val internally, Met at initiation),
#         used by ~14% of E. coli genes
#   TTG — alternative start (codes Leu internally, Met at initiation),
#         used by ~3% of E. coli genes
#
# At the ribosomal initiation site, the initiator fMet-tRNA recognises
# all three codons and inserts formyl-methionine (fMet).  Internally
# (i.e., not at position 1), GTG codes for Val and TTG codes for Leu
# per the standard genetic code.
#
# Usage:
#   - translate() and translate_with_confidence() use the standard
#     genetic code only — they do NOT perform start-codon recognition.
#   - find_orfs(start_codons=BACTERIAL_START_CODONS) recognises all
#     three as starts and translates the first codon as 'M'.
#
# Source: NCBI Genetic Codes Table 1; Blattner et al. (1997) Science 277:1453;
#   RCSB PDB summary of initiation codon usage in E. coli K-12.
BACTERIAL_START_CODONS: frozenset[str] = frozenset({"ATG", "GTG", "TTG"})


class ORFResult(TypedDict):
    """Structured type for ORF results returned by find_orfs."""

    start: int
    end: int
    frame: int
    strand: str
    protein: str
    length: int


def translate(sequence: str, to_stop: bool = True) -> str:
    """
    Translate a DNA sequence to a protein sequence via the standard genetic code.

    This is a deterministic FST: codon → amino acid mapping.
    Translates from frame 0 until a stop codon or end of sequence.

    Args:
        sequence: DNA coding sequence
        to_stop: if True, stop at first stop codon; if False, include stops as '*'

    Returns:
        Protein sequence as single-letter amino acid codes.
    """
    sequence = validate_dna_sequence(sequence)
    if not sequence:
        return ""

    # Warn about partial codons
    if len(sequence) % _CODON_LENGTH != 0:
        logger.warning(
            "Sequence length %d is not a multiple of %d; last %d base(s) will be ignored",
            len(sequence), _CODON_LENGTH, len(sequence) % _CODON_LENGTH,
        )

    protein: list[str] = []
    for i in range(0, len(sequence) - (_CODON_LENGTH - 1), _CODON_LENGTH):
        codon = sequence[i:i + _CODON_LENGTH]
        aa = CODON_TABLE.get(codon)
        if aa is None:
            logger.warning("Unknown codon '%s' at position %d — mapping to 'X'", codon, i)
            aa = "X"
        if aa == "*" and to_stop:
            break
        protein.append(aa)
    return "".join(protein)


def translate_with_confidence(sequence: str, to_stop: bool = True) -> str:
    """Translate a DNA sequence, raising an error for unknown codons.

    This is the strict variant of :func:`translate`.  Unlike
    :func:`translate`, which silently maps unknown codons to ``'X'``
    and logs a warning for partial codons, this function **raises**
    errors in both cases:

    - **Partial codons** at the end of the sequence raise
      :class:`PartialCodonError` (a specific subclass of
      :class:`ValueError`).  This prevents silent data loss from
      truncated sequences.
    - **Unknown codons** raise :class:`ValueError` with a descriptive
      message.

    Use this function when silent fallback is undesirable,
    e.g., in optimisation pipelines where every codon must be valid
    and the sequence must be complete.

    Bacterial start codon handling
        This function translates codons strictly according to the
        standard genetic code (NCBI table 1).  At the **initiation**
        position, bacteria may use GTG or TTG as alternative start
        codons, both translated as **M** (methionine) by the initiator
        fMet-tRNA.  However, this low-level translation function does
        **not** perform ORF detection or start-codon recognition —
        it maps each codon to its standard internal translation.
        For ORF-aware translation with bacterial start codon support,
        use :func:`find_orfs` with
        ``start_codons=BACTERIAL_START_CODONS``.

    Args:
        sequence: DNA coding sequence.  Length must be a multiple of 3.
        to_stop: If True, stop at first stop codon; if False, include
            stops as ``'*'``.

    Returns:
        Protein sequence as single-letter amino acid codes.

    Raises:
        PartialCodonError: If the sequence length is not a multiple of 3
            (i.e., trailing bases do not form a complete codon).
        ValueError: If the sequence contains a codon that is not in the
            standard genetic code (NCBI table 1).
    """
    sequence = validate_dna_sequence(sequence)
    if not sequence:
        return ""

    # Raise specific error for partial codons
    remainder = len(sequence) % _CODON_LENGTH
    if remainder != 0:
        raise PartialCodonError(len(sequence), remainder)

    protein: list[str] = []
    for i in range(0, len(sequence) - (_CODON_LENGTH - 1), _CODON_LENGTH):
        codon = sequence[i:i + _CODON_LENGTH]
        aa = CODON_TABLE.get(codon)
        if aa is None:
            raise ValueError(
                f"Unknown codon '{codon}' at position {i} — "
                f"not in the standard genetic code (NCBI translation table 1). "
                f"Use translate() for a lenient variant that maps unknown codons to 'X'."
            )
        if aa == "*" and to_stop:
            break
        protein.append(aa)
    return "".join(protein)


def compute_cai(
    sequence: str,
    organism: str = "Homo_sapiens",
    species: str | None = None,
) -> float:
    """
    Compute Codon Adaptation Index (CAI) for a coding sequence.

    CAI = geometric mean of relative adaptiveness values of codons used.
    This is a DETERMINISTIC computation: same sequence → same CAI.

    Organism Specification:

        The target organism can be specified using **either** the
        ``organism`` parameter **or** the ``species`` parameter.  Both
        accept the same set of names — short aliases, abbreviated
        binomials, display names, or full canonical names — and both
        map to the same internal representation via
        :func:`~biocompiler.organisms.resolve_organism`.

        If both ``species`` and ``organism`` are provided, ``species``
        takes precedence and a :class:`DeprecationWarning` is emitted.

    Args:
        sequence: DNA coding sequence (starts with ATG)
        organism: Organism name.  Accepts canonical binomials
            (e.g., ``'Homo_sapiens'``, ``'Escherichia_coli'``),
            short keys (``'ecoli'``, ``'human'``), abbreviated
            binomials (``'E_coli'``, ``'h_sapiens'``), or display
            names (``'E. coli'``).  All forms are resolved via
            :func:`~biocompiler.organisms.resolve_organism`.
        species: Alias for ``organism``.  Accepts the same values.
            If provided **together with** ``organism``, ``species``
            takes precedence and a deprecation warning is emitted.
            Prefer using ``organism`` in new code; ``species`` is
            retained for backward compatibility.

    Returns:
        CAI value in [0.0, 1.0]. Returns 0.0 for empty/invalid sequences.

    Raises:
        UnsupportedOrganismError: if the organism is not supported.
    """
    # ── Organism resolution ────────────────────────────────────────
    if species is not None:
        resolved = resolve_organism(species, strict=False)
        if organism != "Homo_sapiens":
            resolved_explicit = resolve_organism(organism, strict=False)
            if resolved != resolved_explicit:
                warnings.warn(
                    f"Both 'species={species!r}' and 'organism={organism!r}' "
                    f"were provided but resolve to different organisms "
                    f"({resolved!r} vs {resolved_explicit!r}). "
                    f"Using 'species' ({resolved!r}). "
                    f"Prefer using only 'organism' in new code.",
                    DeprecationWarning,
                    stacklevel=2,
                )
            else:
                warnings.warn(
                    f"Both 'species' and 'organism' were provided. "
                    f"Prefer using only 'organism' in new code; "
                    f"'species' is retained for backward compatibility.",
                    DeprecationWarning,
                    stacklevel=2,
                )
        else:
            warnings.warn(
                f"The 'species' parameter is deprecated in favor of 'organism'. "
                f"Use organism='{resolved}' instead of "
                f"species='{species}'. Both accept the same aliases.",
                DeprecationWarning,
                stacklevel=2,
            )
        organism = resolved
    else:
        organism = resolve_organism(organism, strict=False)

    sequence = validate_dna_sequence(sequence)
    if not sequence:
        return 0.0
    if organism not in SUPPORTED_ORGANISMS:
        raise UnsupportedOrganismError(organism, SUPPORTED_ORGANISMS)

    adaptiveness = CODON_ADAPTIVENESS_TABLES[organism]

    # ── NUMBA-accelerated path ───────────────────────────────────────
    if _HAS_NUMBA:
        try:
            import numpy as np
            # Build flat adaptiveness array and index mapping for NUMBA kernel
            # First pass: collect codon adaptiveness values and build index array
            codon_list: list[str] = []
            adapt_values: list[float] = []
            codon_to_idx: dict[str, int] = {}
            idx_counter = 0

            # Build unique codon → index mapping
            for i in range(0, len(sequence) - (_CODON_LENGTH - 1), _CODON_LENGTH):
                codon = sequence[i:i + _CODON_LENGTH]
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*" or aa == "M":
                    continue
                if codon not in codon_to_idx:
                    codon_to_idx[codon] = idx_counter
                    w = adaptiveness.get(codon, 0.0)
                    if w <= 0:
                        w = _ZERO_ADAPTIVENESS_EPSILON
                    adapt_values.append(w)
                    idx_counter += 1
                codon_list.append(codon)

            if not codon_list:
                return 0.0

            adapt_array = np.array(adapt_values, dtype=np.float64)
            indices = np.array([codon_to_idx[c] for c in codon_list], dtype=np.int64)
            n_codons = len(codon_list)

            cai = _numba_cai_kernel(adapt_array, indices, n_codons)
            return round(cai, _CAI_ROUND_PRECISION)
        except Exception:
            # Fallback to pure-Python if NUMBA fails at runtime
            pass

    # ── Pure-Python path (original implementation) ───────────────────
    ratios: list[float] = []

    for i in range(0, len(sequence) - (_CODON_LENGTH - 1), _CODON_LENGTH):
        codon = sequence[i:i + _CODON_LENGTH]
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*" or aa == "M":
            continue
        w = adaptiveness.get(codon, 0.0)
        if w <= 0:
            w = _ZERO_ADAPTIVENESS_EPSILON  # Floor for zero-adaptiveness codons
        ratios.append(w)

    if not ratios:
        return 0.0

    # Geometric mean via log-based computation
    log_sum: float = sum(math.log(r) for r in ratios)
    cai = math.exp(log_sum / len(ratios))
    return round(cai, _CAI_ROUND_PRECISION)


def find_orfs(
    sequence: str,
    min_length_aa: int = DEFAULT_MIN_ORF_LENGTH_AA,
    start_codons: set[str] | frozenset[str] | None = None,
) -> list[ORFResult]:
    """Find all Open Reading Frames in all 6 frames (3 forward + 3 reverse complement).

    This is a production feature — the toy model did not do ORF finding at all.

    Bacterial start codons (ATG, GTG, TTG) are recognised when the
    *start_codons* argument is set to :data:`BACTERIAL_START_CODONS`.
    By default only the canonical ATG start codon is used.

    When an alternative start codon (GTG or TTG) initiates an ORF, the
    first amino acid is translated as **M** (methionine), consistent
    with the biological behaviour of bacterial initiation where the
    initiator fMet-tRNA recognises all three codons.

    Args:
        sequence: DNA sequence
        min_length_aa: minimum ORF length in amino acids
        start_codons: Set of codons to treat as start codons.  Defaults
            to ``{STANDARD_START_CODON}`` (i.e. ``{"ATG"}``).  Pass
            :data:`BACTERIAL_START_CODONS` to include GTG and TTG as
            bacterial alternative starts.

    Returns:
        List of ORF dicts with keys: start, end, frame, strand, protein, length
    """
    sequence = validate_dna_sequence(sequence)
    orfs: list[ORFResult] = []

    # Resolve start codon set
    if start_codons is None:
        start_codons = {STANDARD_START_CODON}
    else:
        start_codons = set(start_codons)

    def _find_forward_orfs(seq: str, strand: str) -> list[ORFResult]:
        """Scan a single-strand sequence for ORFs in all 3 forward frames."""
        found: list[ORFResult] = []
        for frame in range(_CODON_LENGTH):
            i: int = frame
            while i < len(seq) - (_CODON_LENGTH - 1):
                codon = seq[i:i + _CODON_LENGTH]
                if codon in start_codons:
                    # Found start codon — translate until stop.
                    # Bacterial alternative starts (GTG, TTG) encode
                    # Met at the initiation position, not their internal
                    # amino acid identity (Val / Leu).
                    protein: list[str] = []
                    orf_start: int = i
                    j: int = i
                    while j < len(seq) - (_CODON_LENGTH - 1):
                        c = seq[j:j + _CODON_LENGTH]
                        if j == orf_start and c in BACTERIAL_START_CODONS:
                            # Initiator Met — even for GTG/TTG
                            aa = "M"
                        else:
                            aa = CODON_TABLE.get(c)
                            if aa is None:
                                logger.warning(
                                    "Unknown codon '%s' at position %d in ORF scan — mapping to 'X'",
                                    c, j,
                                )
                                aa = "X"
                        if aa == "*":
                            break
                        protein.append(aa)
                        j += _CODON_LENGTH
                    if len(protein) >= min_length_aa:
                        found.append({
                            "start": orf_start,
                            "end": j + _CODON_LENGTH,
                            "frame": frame,
                            "strand": strand,
                            "protein": "".join(protein),
                            "length": len(protein),
                        })
                    i = j + _CODON_LENGTH
                else:
                    i += _CODON_LENGTH
        return found

    # Forward strand
    orfs.extend(_find_forward_orfs(sequence, "+"))
    # Reverse complement strand
    rc_seq = reverse_complement(sequence)
    rc_orfs = _find_forward_orfs(rc_seq, "-")
    # Convert RC positions back to original coordinates
    seq_len = len(sequence)
    for orf in rc_orfs:
        orig_end = seq_len - orf["start"]
        orig_start = seq_len - orf["end"]
        orf["start"] = orig_start
        orf["end"] = orig_end
    orfs.extend(rc_orfs)

    orfs.sort(key=lambda o: o["start"])
    return orfs
