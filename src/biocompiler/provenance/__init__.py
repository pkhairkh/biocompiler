"""
BioCompiler Provenance Subpackage
==================================

Decision audit trail, certificate generation, cryptographic integrity,
and verification for gene design optimization.

This package should minimize dependencies beyond ``shared`` modules.
"""

from .tracker import (
    DecisionRecord,
    ProvenanceTracker,
    OptimizationProvenance,
    OptimizationRecord,
    ProvenanceStore,
    generate_provenance_report,
    # Decision category constants
    DECISION_CATEGORY_CAI,
    DECISION_CATEGORY_GT_AVOIDANCE,
    DECISION_CATEGORY_GC_CONTENT,
    DECISION_CATEGORY_RESTRICTION_SITE,
    DECISION_CATEGORY_SPLICE_PREVENTION,
    DECISION_CATEGORY_MUTATION,
    DECISION_CATEGORY_CONSTRAINT_RELAXATION,
    DECISION_CATEGORY_OTHER,
    ALL_DECISION_CATEGORIES,
)

from biocompiler.provenance.certificate import (
    generate_certificate,
    verify_certificate,
    compute_certificate,
    format_certificate,
    VERSION as CERT_VERSION,
    _CERTIFICATE_VERSION,
    _REQUIRED_INPUT_PARAM_KEYS,
    _compute_certificate_hash,
    _compute_gc_content,
    _CURRENT_HASH_VERSION,
    _HASH_ALGORITHM,
    _V2_HASH_PARAM_KEYS,
)

from biocompiler.provenance.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    OptimizationDecisionTrail,
    DecisionProvenanceCollector,
    ProvenanceStore as DecisionProvenanceStore,
)

from .crypto import (
    sign_record,
    verify_record,
    ProvenanceIntegrityError,
    get_provenance_secret,
)

from .reporting import (
    ProvenanceQuery,
    ProvenanceReport,
    explain_position,
)

from biocompiler.provenance.slot_verification import (
    SLOT_PREDICATES,
    is_slot_predicate,
    VerificationEvidence,
    verify_no_cryptic_splice,
    verify_no_cryptic_promoter,
    verify_no_unexpected_tm_domain,
    verify_mrna_secondary_structure,
    verify_co_translational_folding,
    verify_conservation_score,
    verify_codon_optimality,
    verify_structure_predicate,
    verify_stability_predicate,
    verify_solubility_predicate,
    verify_immunogenicity_predicate,
    verify_soundness,
    SoundnessReport,
    SoundnessResult,
)

from biocompiler.provenance.report import generate_report

from biocompiler.provenance.proof_checks import (
    assert_conservative_safe,
    assert_verified_evidence,
    assert_and_pass_iff,
    assert_no_slot_in_pass_list,
    assert_valine_gt_invariant,
    assert_synonymous_preserves_translation,
    assert_verdict_refines,
)

from biocompiler.provenance.decision_provenance import _get_biocompiler_version

from biocompiler.provenance.runtime_evidence import (
    run_all_evidence_checks,
    EvidenceCheckResult,
)

__all__ = [
    # tracker (was provenance.py)
    "DecisionRecord",
    "ProvenanceTracker",
    "OptimizationProvenance",
    "OptimizationRecord",
    "ProvenanceStore",
    "generate_provenance_report",
    "DECISION_CATEGORY_CAI",
    "DECISION_CATEGORY_GT_AVOIDANCE",
    "DECISION_CATEGORY_GC_CONTENT",
    "DECISION_CATEGORY_RESTRICTION_SITE",
    "DECISION_CATEGORY_SPLICE_PREVENTION",
    "DECISION_CATEGORY_MUTATION",
    "DECISION_CATEGORY_CONSTRAINT_RELAXATION",
    "DECISION_CATEGORY_OTHER",
    "ALL_DECISION_CATEGORIES",
    # certificate
    "generate_certificate",
    "verify_certificate",
    "compute_certificate",
    "format_certificate",
    # decision_provenance
    "CodonDecision",
    "ConstraintDecision",
    "OptimizationDecisionTrail",
    "DecisionProvenanceCollector",
    # crypto
    "sign_record",
    "verify_record",
    "ProvenanceIntegrityError",
    "get_provenance_secret",
    # reporting
    "ProvenanceQuery",
    "ProvenanceReport",
    "explain_position",
    # slot_verification
    "SLOT_PREDICATES",
    "is_slot_predicate",
    "VerificationEvidence",
    "verify_no_cryptic_splice",
    "verify_no_cryptic_promoter",
    "verify_no_unexpected_tm_domain",
    "verify_mrna_secondary_structure",
    "verify_co_translational_folding",
    "verify_conservation_score",
    "verify_codon_optimality",
    "verify_structure_predicate",
    "verify_stability_predicate",
    "verify_solubility_predicate",
    "verify_immunogenicity_predicate",
    "verify_soundness",
    "SoundnessReport",
    "SoundnessResult",
    # report
    "generate_report",
    # proof_checks
    "assert_conservative_safe",
    "assert_verified_evidence",
    "assert_and_pass_iff",
    "assert_no_slot_in_pass_list",
    "assert_valine_gt_invariant",
    "assert_synonymous_preserves_translation",
    "assert_verdict_refines",
    # runtime_evidence (NEW: W2-A3)
    "run_all_evidence_checks",
    "EvidenceCheckResult",
]
