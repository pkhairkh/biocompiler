"""
BioCompiler Shared Kernel
=========================

Domain-agnostic primitives shared by ALL subpackages.  This package
must NOT import from any other biocompiler subpackage (sequence/,
optimizer/, solver/, etc.) to avoid circular dependencies.

Contents:
  - types.py          : Core data structures (Verdict, Token, Certificate, …)
  - exceptions.py     : Exception hierarchy (BioCompilerError, …)
  - constants.py      : Biological constants (codon table, restriction enzymes, …)
  - five_valued_logic : 5-valued logic bridge to Lean4 3-valued model
  - thread_safety     : ThreadSafeDict, ThreadSafeDefaultDict, ThreadSafeLazy
"""

# Re-export key symbols for convenient access:
#   from biocompiler.shared import Verdict, BioCompilerError, CODON_TABLE, …

from biocompiler.shared.types import (  # noqa: F401
    SLOTMode,
    Verdict,
    five_valued_and,
    five_valued_or,
    three_valued_and,
    three_valued_or,
    combined_verdict,
    PositionRange,
    Token,
    SpliceIsoform,
    TypeCheckResult,
    Certificate,
)

from biocompiler.shared.exceptions import (  # noqa: F401
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
    EngineError,
    ESMFoldError,
    FoldXError,
    CamSolError,
    ImmunogenicityError,
    OptimizationConstraintError,
    BiosecurityError,
    TranslationVerificationError,
)

from biocompiler.shared.constants import (  # noqa: F401
    CODON_TABLE,
    STOP_CODONS,
    START_CODON,
    AA_TO_CODONS,
    DONOR_CONSENSUS,
    ACCEPTOR_CONSENSUS,
    KOZAK_CONSENSUS,
    INSTABILITY_MOTIF,
    MIN_INTRON_LENGTHS,
    MIN_INTRON_LENGTH,
    POLYPYRIMIDINE_WINDOW,
    POLYPYRIMIDINE_THRESHOLD,
    RESTRICTION_ENZYMES,
    BASE_MAP,
    BASE_REV,
    COMPLEMENT,
    IUPAC_EXPAND,
    VALID_IUPAC_BASES,
    reverse_complement,
    STANDARD_AAS,
    STANDARD_AAS_BLOSUM_ORDER,
    BLOSUM62,
    HYDROPATHY,
    HYDROPHOBIC_AAS,
    DEFAULT_ENGINE_TIMEOUT,
    DEFAULT_BATCH_SIZE,
    DEFAULT_SOLUBILITY_WINDOW,
    DEFAULT_SOLUBILITY_SMOOTHING,
    DEFAULT_MHC_PEPTIDE_LENGTH,
)

from biocompiler.shared.five_valued_logic import (  # noqa: F401
    FiveValuedVerdict,
    to_three_valued,
    confidence_score,
    combine_verdicts,
    verify_five_valued_soundness,
    refinement_is_sound,
    five_valued_not,
)

from .thread_safety import (  # noqa: F401
    ThreadSafeDict,
    ThreadSafeDefaultDict,
    ThreadSafeLazy,
)
