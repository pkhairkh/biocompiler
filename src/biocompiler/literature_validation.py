"""Deprecated: use biocompiler.validation.literature_validation instead."""
import warnings

warnings.warn(
    "biocompiler.literature_validation is deprecated — use biocompiler.validation.literature_validation instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.validation.literature_validation import *  # noqa: F401,F403

__all__ = [
    "LiteratureCase",
    "ValidationResult",
    "DomainReport",
    "CAIValidationResult",
    "SCID_CASES",
    "THALASSEMIA_CASES",
    "AGGREGATION_CASES",
    "IMMUNOGENICITY_CASES",
    "ALL_LITERATURE_CASES",
    "EXTENDED_PUBLISHED_CAI",
    "evaluate_case",
    "run_literature_validation",
    "format_literature_validation_report",
    "validate_cai_against_published",
    "compare_reference_sets",
    "correlation_analysis",
    "cross_validate_with_dnachisel",
    "multi_source_cai_comparison",
    "ALPHA_SYNUCLEIN_FULL",
    "ALPHA_SYNUCLEIN_NAC",
    "AMYLOID_BETA_42",
    "EPO_MATURE",
    "FACTOR_VIII_A2",
    "HBB_EXON1_PLUS_IVS1_WT",
    "HBB_IVS1_110_CONTEXT",
    "HBB_IVS1_1_MUTANT",
    "HBB_IVS1_5_MUTANT",
    "HGH_MATURE",
    "HSA_DOMAIN",
    "HUNTINGTIN_EXON1",
    "IL2RG_CDNA_FRAGMENT",
    "INTERFERON_ALPHA",
    "MLV_LTR_PROMOTER",
    "RAG1_CDNA_FRAGMENT",
    "UBIQUITIN",
]
