"""
BioCompiler Exception Hierarchy

Production-grade error handling — no raw asserts.
Every failure mode has a specific, catchable exception.

Extended with:
- FileFormatError for FASTA/GenBank parsing failures
- SplicingError for NDFST computation failures
- Better InvalidSequenceError with position context
- Typed unsat_core in OptimizationError
- EngineError base class for unified engine error handling
- ESMFoldError, FoldXError refactored to inherit from EngineError
- CamSolError for CamSol solubility engine errors
- ImmunogenicityError for immunogenicity engine errors
- BiosecurityError for biosecurity screening hazards that block optimization
"""

from typing import Any


class BioCompilerError(Exception):
    """Base exception for all BioCompiler errors."""
    pass


class InvalidSequenceError(BioCompilerError):
    """Raised when a DNA/RNA sequence contains invalid characters."""
    def __init__(self, sequence: str, invalid_chars: set[str]):
        self.sequence = sequence
        self.invalid_chars = invalid_chars
        # Find positions of invalid characters for better diagnostics
        positions = [i for i, c in enumerate(sequence) if c in invalid_chars]
        pos_str = ", ".join(str(p) for p in positions[:10])
        if len(positions) > 10:
            pos_str += f", ... ({len(positions)} total)"
        super().__init__(
            f"Invalid DNA bases found in sequence: {invalid_chars}. "
            f"Only A, C, G, T, N are permitted. "
            f"Positions: [{pos_str}]"
        )


class CertificateGenerationError(BioCompilerError):
    """Raised when certificate generation is attempted with failing predicates."""
    def __init__(self, failures: list[Any]):
        self.failures = failures
        details = "; ".join(
            f"{r.predicate}={r.verdict.value}" for r in failures
        )
        super().__init__(
            f"Cannot generate certificate: {len(failures)} predicate(s) failed. "
            f"{details}"
        )


class CertificateVerificationError(BioCompilerError):
    """Raised when certificate verification fails."""
    def __init__(self, reasons: list[str]):
        self.reasons = reasons
        super().__init__(
            f"Certificate verification failed: {'; '.join(reasons)}"
        )


class UnknownPredicateError(BioCompilerError):
    """Raised when an unregistered predicate is encountered."""
    def __init__(self, predicate_name: str):
        self.predicate_name = predicate_name
        super().__init__(
            f"Unknown predicate: '{predicate_name}'. "
            f"Register it with the predicate registry before use."
        )


class OptimizationError(BioCompilerError):
    """Raised when sequence optimization fails to find a valid solution."""
    def __init__(self, reason: str, unsat_core: list[str] | None = None):
        self.unsat_core = unsat_core
        super().__init__(f"Optimization failed: {reason}")


class UnsupportedOrganismError(BioCompilerError):
    """Raised when an organism is not supported for codon optimization."""
    def __init__(self, organism: str, available: list[str]):
        self.organism = organism
        self.available = available
        super().__init__(
            f"Unsupported organism: '{organism}'. "
            f"Available: {available}"
        )


class InvalidProteinError(BioCompilerError):
    """Raised when a protein sequence contains invalid amino acid codes."""
    def __init__(self, protein: str, invalid_chars: set[str]):
        self.protein = protein
        self.invalid_chars = invalid_chars
        super().__init__(
            f"Invalid amino acid(s) in protein: {invalid_chars}. "
            f"Only standard single-letter codes are permitted."
        )


class FileFormatError(BioCompilerError):
    """Raised when a file cannot be parsed (FASTA, GenBank, etc.)."""
    def __init__(self, path: str, format_name: str, reason: str):
        self.path = path
        self.format_name = format_name
        self.reason = reason
        super().__init__(
            f"Cannot parse {format_name} file '{path}': {reason}"
        )


class SplicingError(BioCompilerError):
    """Raised when NDFST computation fails."""
    def __init__(self, reason: str):
        self.reason = reason
        super().__init__(f"Splicing computation error: {reason}")


class MutagenesisError(BioCompilerError):
    """Raised when type-directed mutagenesis fails or encounters invalid state."""
    def __init__(self, reason: str, substitutions_applied: int = 0):
        self.reason = reason
        self.substitutions_applied = substitutions_applied
        super().__init__(
            f"Mutagenesis error: {reason} "
            f"({substitutions_applied} substitutions applied before failure)"
        )


# ==============================================================================
# Engine-specific errors
# ==============================================================================

class EngineError(BioCompilerError):
    """Base exception for all engine-specific errors.

    All computational engine errors (ESMFold, FoldX, CamSol,
    Immunogenicity, etc.) inherit from this class, allowing callers
    to catch any engine failure with a single ``except EngineError``
    handler while still being able to catch specific engine errors
    individually.
    """

    def __init__(self, reason: str, engine: str = "unknown"):
        self.reason = reason
        self.engine = engine
        # Subclasses may set _message before calling super().__init__;
        # only set the default format if the subclass hasn't already.
        if not hasattr(self, "_message"):
            self._message = f"[{engine}] {reason}"
        super().__init__(self._message)

    def __str__(self) -> str:
        return self._message


class ESMFoldError(EngineError):
    """Raised when ESMFold structure prediction fails.

    This covers API communication errors, invalid protein sequences,
    parsing failures of returned PDB data, and timeout conditions
    during structure prediction.
    """

    def __init__(self, reason: str, protein: str | None = None):
        self.reason = reason
        self.protein = protein
        self._message = f"ESMFold prediction failed: {reason}"
        if protein:
            self._message += f" (protein length={len(protein)})"
        super().__init__(self._message, engine="ESMFold")

    def __str__(self) -> str:
        return self._message


class FoldXError(EngineError):
    """Raised when FoldX stability analysis fails.

    This covers command execution errors, missing FoldX installation,
    PDB processing failures, and energy computation errors.
    """

    def __init__(self, reason: str, command: str | None = None):
        self.reason = reason
        self.command = command
        self._message = f"FoldX error: {reason}"
        if command:
            self._message += f" (command: {command})"
        super().__init__(self._message, engine="FoldX")

    def __str__(self) -> str:
        return self._message


class CamSolError(EngineError):
    """Raised when CamSol solubility prediction fails.

    This covers errors during solubility score computation, invalid
    protein input for the CamSol engine, and failures in the
    solubility-guided optimization loop.
    """

    def __init__(self, reason: str, protein: str | None = None):
        self.reason = reason
        self.protein = protein
        self._message = f"CamSol error: {reason}"
        if protein:
            self._message += f" (protein length={len(protein)})"
        super().__init__(self._message, engine="CamSol")

    def __str__(self) -> str:
        return self._message


class ImmunogenicityError(EngineError):
    """Raised when immunogenicity prediction or optimization fails.

    This covers errors during epitope scoring, T-cell binding
    prediction failures, and errors in the immunogenicity-guided
    de-immunization optimization loop.
    """

    def __init__(self, reason: str, protein: str | None = None):
        self.reason = reason
        self.protein = protein
        self._message = f"Immunogenicity error: {reason}"
        if protein:
            self._message += f" (protein length={len(protein)})"
        super().__init__(self._message, engine="Immunogenicity")

    def __str__(self) -> str:
        return self._message


class OptimizationConstraintError(BioCompilerError):
    """Raised when optimization cannot satisfy all constraints in strict mode.

    In strict mode (the default), the optimizer refuses to return sequences
    that have failed predicates.  This exception provides the caller with
    the list of failed predicates, the partial result (for inspection), and
    a suggestion to relax the mode.

    Attributes:
        failed_predicates: List of predicate names that could not be satisfied.
        partial_result: The :class:`OptimizationResult` that was produced but
            not returned, so callers can still inspect it if desired.
    """

    def __init__(self, failed_predicates: list[str], partial_result: Any = None):
        self.failed_predicates = failed_predicates
        self.partial_result = partial_result
        super().__init__(
            f"Optimization failed {len(failed_predicates)} predicate(s) in strict mode: "
            f"{', '.join(failed_predicates)}. "
            f"Set strict_mode=False to allow partial results."
        )


# ==============================================================================
# Safety-specific errors
# ==============================================================================

class BiosecurityError(BioCompilerError):
    """Raised when a biosecurity screen flags a hazardous sequence.

    This covers detection of known toxins, virulence factors,
    or other hazardous biological sequences during optimization
    or screening. Carries structured information about the
    risk level, flagged categories, and individual matches.

    Supports multiple calling conventions:

    1. **Report form** (recommended): pass a :class:`BiosecurityReport`
       object as the first argument.  The report's attributes are
       automatically extracted for the error message.
    2. **Legacy form**: pass ``reason`` (str) plus optional
       ``risk_level``, ``flagged_categories``, and ``matches``.
    3. **Keyword form**: pass ``protein``, ``flagged_pathogens``,
       ``risk_levels``, and ``match_details`` as keyword arguments.

    Attributes:
        report: The :class:`BiosecurityReport` that triggered the error
            (only set when the report-form constructor is used).
        reason: Human-readable description of the hazard.
        risk_level: One of "none", "low", "medium", "high", "critical".
        flagged_categories: Categories that were flagged.
        matches: Individual hazard matches.
        protein: The protein sequence that was screened.
        flagged_pathogens: Names of pathogens that were flagged.
        risk_levels: Risk levels for each flagged pathogen.
        match_details: Details about each match.
    """

    def __init__(
        self,
        reason_or_report: "str | BiosecurityReport | None" = None,  # noqa: F821
        risk_level: str | None = None,
        flagged_categories: list[str] | None = None,
        matches: list[Any] | None = None,
        *,
        protein: str | None = None,
        flagged_pathogens: list[str] | None = None,
        risk_levels: list[str] | None = None,
        match_details: list[str] | None = None,
    ):
        # Detect report-form: if the first argument is a BiosecurityReport
        # (has is_hazardous, risk_level attributes), extract info from it.
        report = None
        if reason_or_report is not None and hasattr(reason_or_report, "is_hazardous") and hasattr(reason_or_report, "risk_level"):
            report = reason_or_report
            self.reason = f"risk_level={report.risk_level}"
            self.risk_level = report.risk_level
            self.flagged_categories = report.flagged_categories
            self.matches = report.matches
        else:
            self.reason = reason_or_report or "hazardous sequence detected"
            self.risk_level = risk_level
            self.flagged_categories = flagged_categories or []
            self.matches = matches or []

        self.report = report

        # Keyword-form attributes
        self.protein = protein or ""
        self.flagged_pathogens = flagged_pathogens or []
        self.risk_levels = risk_levels or []
        self.match_details = match_details or []

        # If risk_level not already set from report/legacy form, derive from risk_levels
        if self.risk_level is None and self.risk_levels:
            priority = {"critical": 4, "high": 3, "medium": 2, "low": 1, "none": 0}
            self.risk_level = max(self.risk_levels, key=lambda r: priority.get(r.lower(), 0))

        # Derive flagged_categories from report if not set and report is available
        if not self.flagged_categories and self.report is not None:
            if hasattr(self.report, "flagged_categories"):
                self.flagged_categories = list(self.report.flagged_categories)

        # Build message
        detail = f" (risk_level={self.risk_level})" if self.risk_level else ""
        if self.flagged_categories:
            detail += f" categories={self.flagged_categories}"
        if self.flagged_pathogens:
            detail += f" pathogens={self.flagged_pathogens}"
        if self.risk_levels:
            detail += f" risk_levels={self.risk_levels}"

        msg = f"BIOSECURITY ALERT: Optimization BLOCKED — {self.reason}{detail}"
        if self.flagged_pathogens:
            msg += "\nConsult your institution's biosafety officer before proceeding."
        super().__init__(msg)


class TranslationVerificationError(BioCompilerError):
    """Raised when optimized DNA does not encode the expected protein.

    This is a critical safety check: after optimization, the translated
    protein must match the original input exactly.  Any mismatch indicates
    a bug in the optimizer, the codon table, or the constraint resolution
    pipeline.

    Supports two calling conventions:

    1. **Rich form** (from :func:`verify_and_raise`): pass structured
       position-level mismatch data, premature stop flag, length info, etc.
    2. **Simple form** (backward compat): pass a reason string and optional
       mismatch list.

    Attributes:
        mismatches: List of position-level mismatches between the translated
            and expected protein.  Each element is either a
            :class:`~biocompiler.protein_verification.PositionMismatch`
            instance or a dict with keys ``position``, ``expected``,
            ``actual``, ``codon_used``.
        has_premature_stop: True if a stop codon was found before the end
            of the expected protein.
        has_stop_codon: True if the DNA ends with a stop codon.
        length_correct: True if the DNA length is consistent with the
            expected protein length.
        translated_protein: The protein obtained by translating the DNA.
        expected_protein: The protein that was expected.
        dna_sequence: The DNA sequence that failed verification.
    """

    def __init__(
        self,
        mismatches: list[Any] | None = None,
        has_premature_stop: bool = False,
        has_stop_codon: bool = False,
        length_correct: bool = True,
        translated_protein: str = "",
        expected_protein: str = "",
        dna_sequence: str = "",
        *,
        # Backward compat: allow simple (reason, mismatches, translated_protein) form
        reason: str | None = None,
    ):
        # Support simple form: TranslationVerificationError("message", mismatches=..., translated_protein=...)
        if reason is not None:
            # Simple form called with positional reason
            self.mismatches = mismatches or []
            self.has_premature_stop = has_premature_stop
            self.has_stop_codon = has_stop_codon
            self.length_correct = length_correct
            self.translated_protein = translated_protein
            self.expected_protein = expected_protein
            self.dna_sequence = dna_sequence
            n = len(self.mismatches)
            detail = f" ({n} mismatch(es))" if n else ""
            super().__init__(f"Translation verification failed: {reason}{detail}")
            return

        # Rich form: from verify_and_raise
        self.mismatches = mismatches or []
        self.has_premature_stop = has_premature_stop
        self.has_stop_codon = has_stop_codon
        self.length_correct = length_correct
        self.translated_protein = translated_protein
        self.expected_protein = expected_protein
        self.dna_sequence = dna_sequence

        # Build a human-readable message
        parts: list[str] = []
        if self.mismatches:
            n = len(self.mismatches)
            parts.append(f"{n} amino acid mismatch(es)")
            for mm in self.mismatches[:5]:
                if isinstance(mm, dict):
                    parts.append(
                        f"  pos {mm['position']}: expected '{mm['expected']}', "
                        f"got '{mm['actual']}' (codon {mm['codon_used']})"
                    )
                else:
                    parts.append(f"  {mm}")
            if n > 5:
                parts.append(f"  ... and {n - 5} more")
        if has_premature_stop:
            parts.append("premature stop codon detected")
        if not length_correct:
            expected_len = len(expected_protein) * 3
            parts.append(
                f"DNA length {len(dna_sequence)} != expected "
                f"{expected_len} or {expected_len + 3}"
            )

        message = (
            f"Translation verification failed: {'; '.join(parts)}. "
            f"Translated: '{translated_protein[:50]}{'...' if len(translated_protein) > 50 else ''}', "
            f"Expected: '{expected_protein[:50]}{'...' if len(expected_protein) > 50 else ''}'"
        )
        super().__init__(message)



