"""
BioCompiler Protein Verification Module
========================================

Double-checks that optimized DNA sequences actually encode the input protein.
This is a critical safety net: after optimization, the translated protein
must match the original input exactly — any mismatch indicates a bug in the
optimizer, the codon table, or the constraint resolution pipeline.

Usage::

    from biocompiler.protein_verification import verify_translation

    result = verify_translation(dna, expected_protein)
    if not result.is_valid:
        for mm in result.mismatches:
            print(f"Position {mm.position}: expected {mm.expected}, got {mm.actual} (codon {mm.codon_used})")

The module also integrates with the optimizer via :func:`verify_and_raise`,
which raises :class:`TranslationVerificationError` on failure.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from ..constants import CODON_TABLE, STOP_CODONS

__all__ = [
    "PositionMismatch",
    "VerificationResult",
    "verify_translation",
    "verify_and_raise",
]

logger = logging.getLogger(__name__)


@dataclass
class PositionMismatch:
    """A single amino acid position where the translated protein differs from expected.

    Attributes:
        position: 0-based amino acid index where the mismatch occurs.
        expected: The expected single-letter amino acid code.
        actual: The actual single-letter amino acid code (or ``'*'`` for stop).
        codon_used: The 3-letter DNA codon that produced the mismatch.
    """

    position: int
    expected: str
    actual: str
    codon_used: str

    def __str__(self) -> str:
        return (
            f"Position {self.position}: expected '{self.expected}', "
            f"got '{self.actual}' (codon: {self.codon_used})"
        )


@dataclass
class VerificationResult:
    """Result of verifying that a DNA sequence encodes the expected protein.

    Attributes:
        is_valid: True if all checks pass (no mismatches, no premature stops,
            proper length, and the translated protein matches expected).
        matches_expected: True if the translated protein exactly matches the
            expected protein (ignoring terminal stop codons).
        mismatches: List of position-level mismatches between translated and
            expected protein.
        has_premature_stop: True if a stop codon was found before the end of
            the expected protein.
        has_stop_codon: True if the last codon of the DNA is a stop codon.
        length_correct: True if the DNA length is consistent with encoding the
            expected protein (with or without a terminal stop codon).
        translated_protein: The full translated protein sequence (including any
            premature stop as ``'*'``; excluding the terminal stop codon if present).
    """

    is_valid: bool
    matches_expected: bool
    mismatches: list[PositionMismatch]
    has_premature_stop: bool
    has_stop_codon: bool
    length_correct: bool
    translated_protein: str


def verify_translation(
    dna: str,
    expected_protein: str,
    organism: str = "",
) -> VerificationResult:
    """Double-check that the optimized DNA sequence encodes exactly the input protein.

    This function translates the DNA using the standard genetic code and performs
    the following checks:

    1. **Length check**: DNA length should be ``3 * len(expected_protein)`` (no stop)
       or ``3 * len(expected_protein) + 3`` (with terminal stop codon).
    2. **Position-by-position comparison**: Every amino acid is compared.
    3. **Premature stop codon detection**: Any stop codon before the end of the
       expected protein region is flagged.
    4. **Terminal stop codon check**: Whether the sequence ends with a stop codon.

    Args:
        dna: The optimized DNA sequence (uppercase A/C/G/T).
        expected_protein: The expected amino acid sequence (single-letter codes,
            no stop codon).
        organism: Optional organism name (reserved for future organism-specific
            genetic code support; currently unused).

    Returns:
        A :class:`VerificationResult` with detailed mismatch information.
    """
    dna = dna.upper().strip()
    expected_protein = expected_protein.upper().strip()

    # ── Determine expected lengths ─────────────────────────────────────
    protein_len = len(expected_protein)
    dna_len_no_stop = protein_len * 3
    dna_len_with_stop = dna_len_no_stop + 3
    dna_len = len(dna)

    # Length check: accept either with or without terminal stop codon
    length_correct = dna_len in (dna_len_no_stop, dna_len_with_stop)
    # If length is neither, we still try to translate what we can

    # ── Translate the DNA codon by codon ──────────────────────────────
    translated_chars: list[str] = []
    codons_used: list[str] = []

    # Determine how many codons to translate
    num_codons = dna_len // 3
    for i in range(num_codons):
        codon = dna[i * 3 : i * 3 + 3]
        if len(codon) < 3:
            break  # Partial codon at end — shouldn't happen with valid input
        aa = CODON_TABLE.get(codon)
        if aa is None:
            # Unknown codon — map to 'X' and flag as mismatch
            logger.warning("Unknown codon '%s' at position %d during verification", codon, i * 3)
            aa = "X"
        translated_chars.append(aa)
        codons_used.append(codon)

    # The full translated string (may include '*' for stop codons)
    full_translated = "".join(translated_chars)

    # ── Determine if terminal stop codon is present ───────────────────
    has_stop_codon = False
    if codons_used:
        last_codon = codons_used[-1]
        last_aa = translated_chars[-1]
        if last_aa == "*" and dna_len == dna_len_with_stop:
            has_stop_codon = True

    # ── Extract the protein portion (excluding terminal stop if present) ──
    if has_stop_codon:
        # Remove the terminal stop from the protein comparison
        protein_translation = full_translated[:-1]
        protein_codons = codons_used[:-1]
    else:
        protein_translation = full_translated
        protein_codons = codons_used

    # ── Check for premature stops ─────────────────────────────────────
    has_premature_stop = False
    for i, aa in enumerate(protein_translation):
        if aa == "*":
            has_premature_stop = True
            break

    # ── Position-by-position comparison ───────────────────────────────
    mismatches: list[PositionMismatch] = []
    matches_expected = True

    # Compare up to the length of whichever is shorter
    compare_len = min(len(protein_translation), protein_len)

    for i in range(compare_len):
        actual_aa = protein_translation[i]
        expected_aa = expected_protein[i]
        codon = protein_codons[i] if i < len(protein_codons) else "???"

        if actual_aa != expected_aa:
            matches_expected = False
            mismatches.append(
                PositionMismatch(
                    position=i,
                    expected=expected_aa,
                    actual=actual_aa,
                    codon_used=codon,
                )
            )

    # If lengths differ, flag remaining positions as mismatches
    if len(protein_translation) < protein_len:
        # Translated protein is shorter than expected
        matches_expected = False
        for i in range(len(protein_translation), protein_len):
            mismatches.append(
                PositionMismatch(
                    position=i,
                    expected=expected_protein[i],
                    actual="<missing>",
                    codon_used="---",
                )
            )
    elif len(protein_translation) > protein_len:
        # Translated protein is longer than expected (extra AAs after expected end)
        matches_expected = False
        for i in range(protein_len, len(protein_translation)):
            codon = protein_codons[i] if i < len(protein_codons) else "???"
            actual = protein_translation[i]
            mismatches.append(
                PositionMismatch(
                    position=i,
                    expected="<end>",
                    actual=actual if actual != "*" else "*",
                    codon_used=codon,
                )
            )

    # ── Overall validity ──────────────────────────────────────────────
    is_valid = (
        matches_expected
        and not has_premature_stop
        and length_correct
    )

    return VerificationResult(
        is_valid=is_valid,
        matches_expected=matches_expected,
        mismatches=mismatches,
        has_premature_stop=has_premature_stop,
        has_stop_codon=has_stop_codon,
        length_correct=length_correct,
        translated_protein=protein_translation,
    )


def verify_and_raise(
    dna: str,
    expected_protein: str,
    organism: str = "",
) -> VerificationResult:
    """Verify translation and raise :class:`TranslationVerificationError` on failure.

    This is the integration hook for the optimizer's return path. After
    optimization completes, call this function with the optimized DNA and the
    original protein. If verification fails, a :class:`TranslationVerificationError`
    is raised with detailed position-level information.

    Args:
        dna: The optimized DNA sequence.
        expected_protein: The original input protein sequence.
        organism: Optional organism name.

    Returns:
        The :class:`VerificationResult` if verification passes.

    Raises:
        TranslationVerificationError: If verification fails (mismatches,
            premature stops, or incorrect length).
    """
    from ..exceptions import TranslationVerificationError

    result = verify_translation(dna, expected_protein, organism)

    if not result.is_valid:
        raise TranslationVerificationError(
            mismatches=result.mismatches,
            has_premature_stop=result.has_premature_stop,
            has_stop_codon=result.has_stop_codon,
            length_correct=result.length_correct,
            translated_protein=result.translated_protein,
            expected_protein=expected_protein,
            dna_sequence=dna,
        )

    return result
