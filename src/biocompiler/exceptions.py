"""
BioCompiler Exception Hierarchy

Production-grade error handling — no raw asserts.
Every failure mode has a specific, catchable exception.
"""


class BioCompilerError(Exception):
    """Base exception for all BioCompiler errors."""
    pass


class InvalidSequenceError(BioCompilerError):
    """Raised when a DNA/RNA sequence contains invalid characters."""
    def __init__(self, sequence: str, invalid_chars: set):
        self.sequence = sequence
        self.invalid_chars = invalid_chars
        super().__init__(
            f"Invalid DNA bases found in sequence: {invalid_chars}. "
            f"Only A, C, G, T, N are permitted."
        )


class CertificateGenerationError(BioCompilerError):
    """Raised when certificate generation is attempted with failing predicates."""
    def __init__(self, failures: list):
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
    def __init__(self, reasons: list):
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
    def __init__(self, reason: str, unsat_core: list | None = None):
        self.unsat_core = unsat_core
        super().__init__(f"Optimization failed: {reason}")


class UnsupportedOrganismError(BioCompilerError):
    """Raised when an organism is not supported for codon optimization."""
    def __init__(self, organism: str, available: list):
        self.organism = organism
        self.available = available
        super().__init__(
            f"Unsupported organism: '{organism}'. "
            f"Available: {available}"
        )


class InvalidProteinError(BioCompilerError):
    """Raised when a protein sequence contains invalid amino acid codes."""
    def __init__(self, protein: str, invalid_chars: set):
        self.protein = protein
        self.invalid_chars = invalid_chars
        super().__init__(
            f"Invalid amino acid(s) in protein: {invalid_chars}. "
            f"Only standard single-letter codes are permitted."
        )
