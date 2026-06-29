"""
BioCompiler Expression Subpackage
=================================

Modules for predicting and scoring gene expression:
translation, codon adaptation, mRNA stability, UTR models, etc.

This package may import from ``shared`` and ``sequence`` modules
but NOT from ``optimizer``.
"""

from biocompiler.expression.translation import (
    translate,
    translate_with_confidence,
    compute_cai,
    find_orfs,
    ORFResult,
    DEFAULT_MIN_ORF_LENGTH_AA,
    BACTERIAL_START_CODONS,
    STANDARD_START_CODON,
    PartialCodonError,
    HAS_NUMBA,
)

from biocompiler.optimizer.cai import (
    _compute_cai_fast,
    _count_dinucs_fast,
    _BatchSwapScorer,
    HAS_NUMBA as CAI_HAS_NUMBA,
    USE_NUMBA as CAI_USE_NUMBA,
)

from biocompiler.expression.tai import (
    compute_tai,
    calculate_tai,
    compute_codon_weights,
    TRNA_GENE_COPIES,
    WOBBLE_RULES,
    WOBBLE_EFFICIENCY,
    SUPPORTED_ORGANISMS_TAI,
)

from biocompiler.expression.mrna_stability import (
    STABILITY_MOTIFS,
    MRNAStabilityScore,
    score_mrna_stability,
    compute_mrna_half_life_score,
    predict_mrna_stability,
    suggest_mutations_for_stability,
)

from biocompiler.expression.codon_pair_scoring import (
    compute_cpb,
    compute_cpb_score,
    estimate_cpb_from_codon_freq,
    get_codon_pair_data,
    score_codon_pair,
    suggest_better_pair,
)

from biocompiler.expression.utr_models import (
    UTRConfig,
    ORGANISM_UTR_CONFIGS,
    AVAILABLE_ORGANISMS,
    score_5utr,
    score_3utr,
    suggest_5utr,
    suggest_3utr,
)

from biocompiler.expression.expression_predictor import (
    ExpressionPrediction,
    predict_expression,
    ExpressionPredictor,
)

from biocompiler.expression.sharp_li_tables import (
    get_sharp_li_table,
    set_sharp_li_table,
    _SharpLiState,
    # Re-exported constants from organisms.sharp_li_reference
    _CODON_TABLE,
    _AA_TO_CODONS,
    _STOP_CODONS,
    ECOLI_SHARP_LI_REFERENCE_GENES,
    ECOLI_SHARP_LI_CODON_USAGE,
    ECOLI_SHARP_LI_CAI_WEIGHTS,
    YEAST_SHARP_LI_REFERENCE_GENES,
    YEAST_SHARP_LI_CODON_USAGE,
    YEAST_SHARP_LI_CAI_WEIGHTS,
    SHARP_LI_PUBLISHED_CAI,
    _REFERENCE_WEIGHTS,
    SHARP_LI_REFERENCE_GENES,
    SHARP_LI_CODON_USAGE,
    SHARP_LI_CAI_WEIGHTS,
    compute_cai_with_reference,
    get_sharp_li_cai_weights,
)

__all__ = [
    # translation
    "translate",
    "translate_with_confidence",
    "compute_cai",
    "find_orfs",
    "ORFResult",
    "DEFAULT_MIN_ORF_LENGTH_AA",
    "BACTERIAL_START_CODONS",
    "STANDARD_START_CODON",
    "PartialCodonError",
    "HAS_NUMBA",
    # cai
    "_compute_cai_fast",
    "_count_dinucs_fast",
    "_BatchSwapScorer",
    "CAI_HAS_NUMBA",
    "CAI_USE_NUMBA",
    # tai
    "compute_tai",
    "calculate_tai",
    "compute_codon_weights",
    "TRNA_GENE_COPIES",
    "WOBBLE_RULES",
    "WOBBLE_EFFICIENCY",
    "SUPPORTED_ORGANISMS_TAI",
    # mrna_stability
    "STABILITY_MOTIFS",
    "MRNAStabilityScore",
    "score_mrna_stability",
    "compute_mrna_half_life_score",
    "predict_mrna_stability",
    "suggest_mutations_for_stability",
    # codon_pair_scoring
    "compute_cpb",
    "compute_cpb_score",
    "estimate_cpb_from_codon_freq",
    "get_codon_pair_data",
    "score_codon_pair",
    "suggest_better_pair",
    # utr_models
    "UTRConfig",
    "ORGANISM_UTR_CONFIGS",
    "AVAILABLE_ORGANISMS",
    "score_5utr",
    "score_3utr",
    "suggest_5utr",
    "suggest_3utr",
    # expression_predictor
    "ExpressionPrediction",
    "predict_expression",
    "ExpressionPredictor",
    # sharp_li_tables
    "get_sharp_li_table",
    "set_sharp_li_table",
    "_SharpLiState",
    "compute_cai_with_reference",
    "get_sharp_li_cai_weights",
]
