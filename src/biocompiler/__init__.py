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

# Type system imports
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

try:
    from .certificate import generate_certificate, verify_certificate
except ImportError:
    pass

try:
    from .maxentscan import (
        score_donor, score_acceptor, scan_splice_sites,
        max_donor_score, max_acceptor_score,
    )
except ImportError:
    pass

from .optimization import BioOptimizer, optimize_sequence, OptimizationResult

try:
    from .grammar_loader import load_grammar, grammar_to_predicate_params, load_builtin_grammar, list_builtin_grammars
except ImportError:
    pass

try:
    from .export import export_fasta, export_genbank, export_genbank_with_certificate, export_multi_fasta
except ImportError:
    pass

try:
    from .report import generate_report
except ImportError:
    pass

try:
    from .benchmark import run_benchmarks, BenchmarkReport
except ImportError:
    try:
        from .benchmark import run_benchmark as run_benchmarks
        BenchmarkReport = None
    except ImportError:
        run_benchmarks = None
        BenchmarkReport = None

try:
    from .organism_db import OrganismDatabase, get_database
except ImportError:
    pass

try:
    from .tissue_data import get_tissue_weights, list_available_tissues, add_custom_tissue
except ImportError:
    pass

try:
    from .dna_chisel_compat import compare_optimizers, run_comparative_benchmark
except ImportError:
    compare_optimizers = None
    run_comparative_benchmark = None

try:
    from .tool_comparison import (
        run_head_to_head_benchmark, compare_tools, HeadToHeadReport,
        format_head_to_head_text, format_head_to_head_json,
        is_dna_chisel_available,
    )
except ImportError:
    pass

try:
    from .dataset_validation import run_dataset_validation, DatasetValidationReport
except ImportError:
    pass

try:
    from .import_seq import import_fasta, import_genbank, import_sequence
except ImportError:
    pass

try:
    from .biopython_compat import to_seqrecord, from_seqrecord, optimize_to_seqrecord
except ImportError:
    pass

try:
    from .jupyter import display_sequence, display_optimization_result, display_type_check, plot_gc_content, plot_codon_usage
except ImportError:
    pass

try:
    from .mutagenesis import (
        MutagenesisResult, AASubstitution, BLOSUM62,
        GT_MANDATORY_AAS, AG_MANDATORY_AAS,
        is_gt_mandatory, is_ag_mandatory, diagnose_optimizer_weakness,
        force_gt_free_reoptimization,
        type_directed_mutagenesis, find_unrepairable_cryptic_donors,
        find_unrepairable_cryptic_acceptors, propose_substitutions,
        apply_substitution,
    )
except ImportError:
    pass

__all__ = [
    "__version__",
    "Verdict", "Token", "PositionRange", "SpliceIsoform",
    "TypeCheckResult", "Certificate",
    "three_valued_and", "three_valued_or", "combined_verdict",
    "BioCompilerError", "InvalidSequenceError",
    "CertificateGenerationError", "CertificateVerificationError",
    "UnknownPredicateError", "OptimizationError",
    "UnsupportedOrganismError", "InvalidProteinError",
    "FileFormatError", "SplicingError", "MutagenesisError",
    "validate_dna_sequence", "gc_content", "scan_sequence",
    "translate", "compute_cai", "find_orfs",
    "compute_splice_isoforms",
    "generate_certificate", "verify_certificate",
    "score_donor", "score_acceptor", "scan_splice_sites",
    "max_donor_score", "max_acceptor_score",
    "BioOptimizer",
    "load_grammar", "grammar_to_predicate_params",
    "load_builtin_grammar", "list_builtin_grammars",
    "export_fasta", "export_genbank", "export_genbank_with_certificate", "export_multi_fasta",
    "generate_report",
    "run_benchmarks", "BenchmarkReport",
    "OrganismDatabase", "get_database",
    "get_tissue_weights", "list_available_tissues", "add_custom_tissue",
    "compare_optimizers", "run_comparative_benchmark",
    "run_head_to_head_benchmark", "compare_tools", "HeadToHeadReport",
    "format_head_to_head_text", "format_head_to_head_json",
    "is_dna_chisel_available",
    "run_dataset_validation", "DatasetValidationReport",
    "import_fasta", "import_genbank", "import_sequence",
    "to_seqrecord", "from_seqrecord", "optimize_to_seqrecord",
    "display_sequence", "display_optimization_result", "display_type_check",
    "plot_gc_content", "plot_codon_usage",
    "MutagenesisResult", "AASubstitution", "BLOSUM62",
    "GT_MANDATORY_AAS", "AG_MANDATORY_AAS",
    "is_gt_mandatory", "is_ag_mandatory", "diagnose_optimizer_weakness",
    "force_gt_free_reoptimization",
    "type_directed_mutagenesis", "find_unrepairable_cryptic_donors",
    "find_unrepairable_cryptic_acceptors", "propose_substitutions",
    "apply_substitution",
]
