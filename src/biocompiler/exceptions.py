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
