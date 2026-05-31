"""
BioCompiler — Machine-Verified Gene Design

A compiler framework for human protein synthesis using intermediate
representations. Pipeline:

  Scanner → NDFST Splicing → Translation → Type Check → Certificate → Verify

All computation is DETERMINISTIC: same input always produces identical output.
"""

__version__ = "7.2.0"

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

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
    MutagenesisError,
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
    analyze_codon_at_position,
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
from .optimization import optimize_sequence, OptimizationResult
from .grammar_loader import load_grammar, grammar_to_predicate_params
from .export import export_fasta, export_genbank, export_genbank_with_certificate, export_multi_fasta
from .report import generate_report
from .benchmark import run_benchmarks, BenchmarkReport
from .organism_db import OrganismDatabase, get_database
from .tissue_data import get_tissue_weights, list_available_tissues, add_custom_tissue
from .dna_chisel_compat import compare_optimizers, run_comparative_benchmark
from .dataset_validation import run_dataset_validation, DatasetValidationReport
from .import_seq import import_fasta, import_genbank, import_sequence
from .biopython_compat import to_seqrecord, from_seqrecord, optimize_to_seqrecord
from .jupyter import display_sequence, display_optimization_result, display_type_check, plot_gc_content, plot_codon_usage
from .mutagenesis import (
    MutagenesisResult,
    AASubstitution,
    BLOSUM62,
    GT_MANDATORY_AAS,
    AG_MANDATORY_AAS,
    is_gt_mandatory,
    is_ag_mandatory,
    diagnose_optimizer_weakness,
    force_gt_free_reoptimization,
    type_directed_mutagenesis,
    find_unrepairable_cryptic_donors,
    find_unrepairable_cryptic_acceptors,
    propose_substitutions,
    apply_substitution,
)

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
    "FileFormatError", "SplicingError", "MutagenesisError",
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
    "evaluate_all_predicates", "analyze_codon_at_position", "predicate_registry",
    # Certificate
    "generate_certificate", "verify_certificate",
    # MaxEntScan
    "score_donor", "score_acceptor", "scan_splice_sites",
    "max_donor_score", "max_acceptor_score",
    # Optimization
    "optimize_sequence", "OptimizationResult",
    # Grammar
    "load_grammar", "grammar_to_predicate_params",
    # Export
    "export_fasta", "export_genbank", "export_genbank_with_certificate", "export_multi_fasta",
    # Report
    "generate_report",
    # Benchmark
    "run_benchmarks", "BenchmarkReport",
    # Database
    "OrganismDatabase", "get_database",
    # Tissue Data (GTEx)
    "get_tissue_weights", "list_available_tissues", "add_custom_tissue",
    # DNA Chisel Compatibility
    "compare_optimizers", "run_comparative_benchmark",
    # Dataset Validation
    "run_dataset_validation", "DatasetValidationReport",
    # Import
    "import_fasta", "import_genbank", "import_sequence",
    # BioPython Interop
    "to_seqrecord", "from_seqrecord", "optimize_to_seqrecord",
    # Jupyter Integration
    "display_sequence", "display_optimization_result", "display_type_check",
    "plot_gc_content", "plot_codon_usage",
    # Mutagenesis Engine
    "MutagenesisResult", "AASubstitution", "BLOSUM62",
    "GT_MANDATORY_AAS", "AG_MANDATORY_AAS",
    "is_gt_mandatory", "is_ag_mandatory", "diagnose_optimizer_weakness",
    "force_gt_free_reoptimization",
    "type_directed_mutagenesis", "find_unrepairable_cryptic_donors",
    "find_unrepairable_cryptic_acceptors", "propose_substitutions",
    "apply_substitution",
]
