"""
BioCompiler — Machine-Verified Gene Design

A compiler framework for human protein synthesis using intermediate
representations. Pipeline:

  Scanner → NDFST Splicing → Translation → Type Check → Certificate → Verify

All computation is DETERMINISTIC: same input always produces identical output.
"""

__version__ = "2.2.0"

from .types import (
    Verdict,
    Token,
    PositionRange,
    SpliceIsoform,
    TypeCheckResult,
    Certificate,
    three_valued_and,
    three_valued_or,
    combined_verdict,
)
from .exceptions import (
    BioCompilerError,
    InvalidSequenceError,
    CertificateGenerationError,
    CertificateVerificationError,
    UnknownPredicateError,
    OptimizationError,
    UnsupportedOrganismError,
    InvalidProteinError,
    FileFormatError,
    SplicingError,
)
from .scanner import (
    validate_dna_sequence,
    gc_content,
    scan_sequence,
)
from .translation import (
    translate,
    compute_cai,
    find_orfs,
)
from .splicing import compute_splice_isoforms
from .type_system import (
    evaluate_no_cryptic_splice,
    evaluate_splice_correct,
    evaluate_gc_in_range,
    evaluate_codon_adapted,
    evaluate_no_restriction_site,
    evaluate_in_frame,
    evaluate_no_instability_motif,
    evaluate_no_cpg_island,
    evaluate_all_predicates,
    registry as predicate_registry,
)
from .certificate import (
    generate_certificate,
    verify_certificate,
)
from .maxentscan import (
    score_donor,
    score_acceptor,
    scan_splice_sites,
    max_donor_score,
    max_acceptor_score,
)
from .optimization import optimize_sequence
from .grammar_loader import load_grammar, grammar_to_predicate_params

__all__ = [
    # Version
    "__version__",
    # Types
    "Verdict", "Token", "PositionRange", "SpliceIsoform",
    "TypeCheckResult", "Certificate",
    "three_valued_and", "three_valued_or", "combined_verdict",
    # Exceptions
    "BioCompilerError", "InvalidSequenceError",
    "CertificateGenerationError", "CertificateVerificationError",
    "UnknownPredicateError", "OptimizationError",
    "UnsupportedOrganismError", "InvalidProteinError",
    "FileFormatError", "SplicingError",
    # Scanner
    "validate_dna_sequence", "gc_content", "scan_sequence",
    # Translation
    "translate", "compute_cai", "find_orfs",
    # Splicing
    "compute_splice_isoforms",
    # Type System
    "evaluate_no_cryptic_splice", "evaluate_splice_correct",
    "evaluate_gc_in_range", "evaluate_codon_adapted",
    "evaluate_no_restriction_site", "evaluate_in_frame",
    "evaluate_no_instability_motif", "evaluate_no_cpg_island",
    "evaluate_all_predicates", "predicate_registry",
    # Certificate
    "generate_certificate", "verify_certificate",
    # MaxEntScan
    "score_donor", "score_acceptor", "scan_splice_sites",
    "max_donor_score", "max_acceptor_score",
    # Optimization
    "optimize_sequence",
    # Grammar
    "load_grammar", "grammar_to_predicate_params",
]
