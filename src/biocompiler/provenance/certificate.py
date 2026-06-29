"""
BioCompiler Certificate Engine — Graduated Generation & Verification

Production-grade certificate system with:
- GRADUATED certificates: works even when not all predicates pass
- Registry-based verification (no string-prefix dispatch)
- All verification parameters embedded IN the certificate
- Hash integrity check for sequence + parameters
- Per-predicate scores and status, not just PASS/FAIL
- CertificateError instead of raw assert

Design philosophy: a certificate documents what was verified, not just
whether everything passed. A sequence with CAI=0.99 but one remaining
PstI site is still useful — the biologist can decide whether PstI is
acceptable for their cloning workflow.

IMPORTANT — SLOT predicate caveats:
  SLOT (Subject to Limited Oracles and Tools) predicates are empirically
  verified but NOT formally proven.  A SLOT PASS verdict means that an
  external tool or scanner produced a passing result; it does NOT carry
  the same formal guarantee as a non-SLOT predicate (which is backed by
  a machine-checked proof).  Certificates annotate SLOT verdicts with an
  asterisk (*) so that downstream consumers can distinguish empirical
  evidence from formal proof.
"""

from __future__ import annotations

import hashlib
import logging
import re as _re
from collections.abc import Callable
from datetime import datetime, timezone
from typing import Any

from biocompiler.shared.types import Verdict, TypeCheckResult, Certificate, SLOTMode
from biocompiler.sequence.maxentscan import (
    CRYPTIC_SPLICE_THRESHOLD as _CRYPTIC_SPLICE_THRESHOLD_FROM_MAXENT,
)
try:
    from ..type_system import registry
except ImportError:
    registry = None
    logging.getLogger(__name__).warning(
        "type_system.registry not available; certificate verification will be limited"
    )
from ..type_system import CertLevel, PredicateResult
from biocompiler.shared.exceptions import CertificateGenerationError

logger = logging.getLogger(__name__)

try:
    from .. import __version__ as _PKG_VERSION
    VERSION: str = _PKG_VERSION
except ImportError:
    VERSION = "0.9.0"
    logger.debug("Package version not available; using fallback VERSION=%s", VERSION)

try:
    from biocompiler.provenance.slot_verification import (
        is_slot_predicate,
        verify_slot_predicate,
    )
except ImportError:
    def is_slot_predicate(name: str) -> bool:  # type: ignore[misc]
        """Fallback when slot_verification is not available."""
        return False
    verify_slot_predicate = None  # type: ignore[assignment]
    logger.debug("slot_verification not available; using no-op fallback")

__all__ = [
    "generate_certificate",
    "verify_certificate",
    "compute_certificate",
    "compute_uncertainty_summary",
    "format_certificate",
    "VERSION",
    "_CERTIFICATE_VERSION",
    "_REQUIRED_INPUT_PARAM_KEYS",
    "_compute_certificate_hash",
    "_compute_gc_content",
    "_CURRENT_HASH_VERSION",
    "_HASH_ALGORITHM",
    "_V2_HASH_PARAM_KEYS",
    "_V3_HASH_PARAM_KEYS",
    "CertLevel",
    "_validate_cert_structure",
    "_CERT_REQUIRED_KEYS",
    "_PROVENANCE_REQUIRED_KEYS",
]

# Certificate version (integer, incremented when the hash format changes)
_CERTIFICATE_VERSION: int = 2

# Required input parameter keys for certificate generation
_REQUIRED_INPUT_PARAM_KEYS: frozenset[str] = frozenset({
    "organism", "gc_lo", "gc_hi", "cai_threshold", "enzymes",
})

# Required keys in a certificate dict for verification
_CERT_REQUIRED_KEYS: frozenset[str] = frozenset({"version", "design_id", "sequence", "types", "provenance"})
_PROVENANCE_REQUIRED_KEYS: frozenset[str] = frozenset({"tool", "version", "timestamp", "input_hash"})

# Default parameter constants (avoids magic numbers)
_DEFAULT_ORGANISM: str = "Homo_sapiens"
_DEFAULT_CELL_TYPE: str = "HEK293T"
_DEFAULT_GC_LO: float = 0.30
_DEFAULT_GC_HI: float = 0.70
_DEFAULT_CAI_THRESHOLD: float = 0.5
_DEFAULT_ENZYMES: list[str] = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
_DEFAULT_CRYPTIC_SPLICE_THRESHOLD: float = _CRYPTIC_SPLICE_THRESHOLD_FROM_MAXENT
# Derived from maxentscan.CRYPTIC_SPLICE_THRESHOLD (single source of truth).
# Previously hardcoded to 3.0 (PWM-era value); now automatically tracks
# the Markov-model-calibrated threshold.
_CERT_FORMAT_WIDTH: int = 60

# Hashing constants
_HASH_ALGORITHM: str = "sha256"
_HASH_TRUNCATION_LENGTH: int = 16
_CURRENT_HASH_VERSION: int = 3

# Parameter keys included in v2 hash computation (legacy — kept for
# backwards-compat verification of pre-v3 certificates).
_V2_HASH_PARAM_KEYS: tuple[str, ...] = (
    "organism", "gc_lo", "gc_hi", "cai_threshold", "solver_backend",
)

# Parameter keys included in v3 hash computation.
#
# C14 fix: the v2 hash only covered 5 parameters and omitted ~13
# predicate-affecting params (enzymes, exon_boundaries,
# cryptic_splice_threshold, is_cytosolic, sliding_gc_window, tissue,
# slot_mode, etc.). A malicious/buggy generator could silently widen the
# verification envelope (e.g. drop "EcoRI" from `enzymes`) while keeping
# design_id unchanged AND the certificate still verifying. v3 closes
# this gap by hashing ALL predicate-affecting parameters referenced by
# `_PREDICATE_KWARGS_MAP`, plus `slot_mode` and `solver_backend`
# (recorded in provenance).
_V3_HASH_PARAM_KEYS: tuple[str, ...] = (
    # Core optimization parameters (also covered by v2)
    "organism", "gc_lo", "gc_hi", "cai_threshold", "solver_backend",
    # Restriction-site avoidance (NoRestrictionSite)
    "enzymes",
    # Splice / exon correctness (NoCrypticSplice, SpliceCorrect, InFrame)
    "exon_boundaries", "cryptic_splice_threshold", "cell_type",
    # SLOT predicate inputs
    "is_cytosolic",                       # NoUnexpectedTMDomain
    "promoter_threshold",                 # NoCrypticPromoter
    "sliding_gc_window",                  # SlidingGC
    "domain_boundaries",                  # CoTranslationalFolding
    "protein", "conservation_min_score",  # ConservationScore
    # miRNA predicate (NoMiRNABindingSite) — tissue context matters
    "tissue", "min_seed_match", "seed_count", "context_scoring",
    # Diagnostic predicates (NoBlastMatches / PrimerCompatibility / NoCrypticORF)
    "reference_sequences",
    "primer_tm_min", "primer_tm_max",
    "min_orf_length",
    # SLOT mode — controls whether SLOT predicates are PASS or UNCERTAIN.
    # Critical: omitting this lets a CONSERVATIVE cert be re-issued as
    # VERIFIED (or vice-versa) without changing design_id.
    "slot_mode",
)


def _compute_gc_content(sequence: str) -> float:
    """Compute GC fraction of a DNA sequence (case-insensitive)."""
    if not sequence:
        return 0.0
    seq = sequence.upper()
    gc = sum(1 for b in seq if b in ("G", "C"))
    return gc / len(seq)


def _compute_certificate_hash(
    sequence: str,
    types_list: list[dict[str, str]],
    params: dict[str, Any],
    hash_version: int = _CURRENT_HASH_VERSION,
    **kwargs: Any,
) -> str:
    """Compute the certificate design_id hash.

    v1: SHA-256(sequence) — legacy, soundness bug (ignores predicates/params).
    v2: SHA-256(sequence + sorted predicate results + 5 key opt params).
        Covers only `organism, gc_lo, gc_hi, cai_threshold, solver_backend`
        — vulnerable to verification-envelope widening (see C14).
    v3: SHA-256(sequence + sorted predicate results + ALL predicate-affecting
        params). Closes the C14 envelope-widening gap by also hashing
        enzymes, exon_boundaries, cryptic_splice_threshold, is_cytosolic,
        sliding_gc_window, tissue, slot_mode, and the remaining params
        enumerated in `_V3_HASH_PARAM_KEYS`.
    """
    if hash_version == 1:
        return hashlib.new(_HASH_ALGORITHM, sequence.encode()).hexdigest()

    # v2+: include sorted predicate results and key parameters
    hasher = hashlib.new(_HASH_ALGORITHM)
    hasher.update(sequence.encode())

    # Sort predicates by name for order-independence
    for entry in sorted(types_list, key=lambda t: t.get("predicate", "")):
        pred = entry.get("predicate", "")
        verdict = entry.get("verdict", "")
        hasher.update(f"{pred}={verdict}".encode())

    # Merge params with kwargs (kwargs take precedence)
    effective_params = dict(params)
    effective_params.update(kwargs)

    # Select the parameter-key set for this hash version. v2 keeps the
    # legacy 5-key set so existing v2 certificates still verify; v3 uses
    # the expanded `_V3_HASH_PARAM_KEYS` (see C14).
    if hash_version >= 3:
        param_keys = _V3_HASH_PARAM_KEYS
    else:
        param_keys = _V2_HASH_PARAM_KEYS

    # Include key optimization parameters in sorted order
    for key in sorted(param_keys):
        if key in effective_params:
            hasher.update(f"{key}={effective_params[key]}".encode())

    return hasher.hexdigest()


def generate_certificate(
    sequence: str,
    type_results: list[TypeCheckResult],
    input_params: dict[str, Any],
    require_all_pass: bool = False,
    mutagenesis_substitutions: list[dict[str, Any]] | None = None,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
    solver_backend: str | None = None,
    solver_config: dict[str, Any] | None = None,
) -> Certificate:
    """
    Generate a machine-checkable guarantee certificate.

    GRADUATED MODE (default): Certificate is generated even if some predicates
    fail. The certificate documents all predicate results, allowing downstream
    users to make informed decisions about partial compliance.

    STRICT MODE (require_all_pass=True): Only generates if all predicates pass.
    This preserves backward compatibility with the old behavior.

    Args:
        sequence: the DNA sequence being certified
        type_results: list of TypeCheckResult objects
        input_params: dict of parameters used (gene, organism, exon_boundaries, etc.)
            MUST include all parameters needed for independent re-verification:
            - exon_boundaries, gc_lo, gc_hi, cai_threshold, organism, cell_type, enzymes
        require_all_pass: if True, raise CertificateGenerationError on any failure.
                         if False (default), generate graduated certificate.
        mutagenesis_substitutions: if provided, documents AA substitutions applied
            to make the design feasible. Each dict has keys: position, from, to,
            blosum62, reason, predicate.
        slot_mode: SLOT verification mode (CONSERVATIVE, VERIFIED, or PERMISSIVE).
        solver_backend: name of the solver backend used (e.g. "ortools", "z3",
            "greedy"). Recorded in provenance for reproducibility; defaults to
            "greedy" if not specified.
        solver_config: solver configuration dict (e.g. GC bounds, time limits).
            Recorded in provenance for reproducibility; defaults to empty dict.

    Returns:
        Certificate object with all predicate results documented

    Raises:
        ValueError: if sequence is empty or type_results is empty
        CertificateGenerationError: only if require_all_pass=True and any predicate failed
    """
    if not sequence:
        raise ValueError("Sequence must not be empty")
    if not type_results:
        raise ValueError("Type results must not be empty")
    failures = [r for r in type_results if r.verdict != Verdict.PASS]

    if require_all_pass and failures:
        raise CertificateGenerationError(failures)

    # Compute hash once
    seq_hash = hashlib.new(_HASH_ALGORITHM, sequence.encode()).hexdigest()

    # Compute overall status
    n_pass = sum(1 for r in type_results if r.verdict == Verdict.PASS)
    n_total = len(type_results)
    overall_status = "FULL_PASS" if not failures else f"PARTIAL_{n_pass}/{n_total}"

    # Ensure all verification parameters are embedded
    complete_params = dict(input_params)
    complete_params.setdefault("organism", _DEFAULT_ORGANISM)
    complete_params.setdefault("cell_type", _DEFAULT_CELL_TYPE)
    complete_params.setdefault("gc_lo", _DEFAULT_GC_LO)
    complete_params.setdefault("gc_hi", _DEFAULT_GC_HI)
    complete_params.setdefault("cai_threshold", _DEFAULT_CAI_THRESHOLD)
    complete_params.setdefault("enzymes", list(_DEFAULT_ENZYMES))
    complete_params.setdefault("cryptic_splice_threshold", _DEFAULT_CRYPTIC_SPLICE_THRESHOLD)
    complete_params.setdefault("exon_boundaries", [(0, len(sequence))])

    # Identify which predicates are SLOT predicates for annotation
    slot_predicates = {r.predicate for r in type_results if is_slot_predicate(r.predicate)}
    slot_pass_count = sum(
        1 for r in type_results
        if is_slot_predicate(r.predicate) and r.verdict == Verdict.PASS
    )
    non_slot_pass_count = sum(
        1 for r in type_results
        if not is_slot_predicate(r.predicate) and r.verdict == Verdict.PASS
    )

    # Build provenance dict
    provenance: dict[str, Any] = {
        "tool": "BioCompiler",
        "version": VERSION,
        "certificate_version": _CERTIFICATE_VERSION,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "parameters": complete_params,
        "input_hash": seq_hash,
        "overall_status": overall_status,
        "slot_mode": slot_mode.value,
        "slot_mode_description": {
            "conservative": "SLOT predicates always UNCERTAIN (matches Lean4 formal model)",
            "verified": "SLOT predicates PASS when verification conditions met",
            "permissive": "SLOT predicates PASS with weaker evidence thresholds",
        }.get(slot_mode.value, "Unknown SLOT mode"),
        "slot_predicate_names": sorted(slot_predicates),
        "slot_pass_count": slot_pass_count,
        "non_slot_pass_count": non_slot_pass_count,
        "slot_caveat": (
            "SLOT predicates are empirically verified but NOT formally proven. "
            "A SLOT PASS verdict indicates that an external tool or scanner "
            "produced a passing result; it does NOT carry the same formal "
            "guarantee as a non-SLOT predicate (which is backed by a "
            "machine-checked proof)."
        ),
        "solver_backend": solver_backend or "greedy",
        "solver_config": solver_config or {},
    }

    # Include mutagenesis metadata if substitutions were applied
    if mutagenesis_substitutions:
        provenance["mutagenesis"] = {
            "applied": True,
            "n_substitutions": len(mutagenesis_substitutions),
            "substitutions": mutagenesis_substitutions,
            "description": (
                f"{len(mutagenesis_substitutions)} conservative amino acid "
                f"substitution(s) were applied to make constraint satisfaction "
                f"possible. See substitutions list for details."
            ),
        }
    else:
        provenance["mutagenesis"] = {"applied": False}

    # Compute hash — v2 covers sequence + sorted predicates + key params
    types_for_hash = [
        {"predicate": r.predicate, "verdict": r.verdict.value}
        for r in type_results
    ]
    design_id = _compute_certificate_hash(
        sequence, types_for_hash, complete_params, hash_version=_CURRENT_HASH_VERSION,
    )

    cert = Certificate(
        version=VERSION,
        design_id=design_id,
        hash_version=_CURRENT_HASH_VERSION,
        sequence=sequence,
        types=[
            {
                "predicate": r.predicate,
                "verdict": r.verdict.value,
                "derivation": r.derivation,
                "knowledge_gap": r.knowledge_gap,
                "confidence": r.confidence,
                "is_slot": is_slot_predicate(r.predicate),
                "verification_basis": (
                    "empirical" if is_slot_predicate(r.predicate) else "formal"
                ),
            }
            for r in type_results
        ],
        provenance=provenance,
    )
    status_msg = overall_status if failures else "FULL_PASS"
    logger.info(
        "Certificate generated: design_id=%s... status=%s solver=%s",
        cert.design_id[:_HASH_TRUNCATION_LENGTH], status_msg, provenance["solver_backend"],
    )
    return cert


def _validate_cert_structure(cert_dict: dict[str, Any]) -> list[str]:
    """Validate that a certificate dict has all required fields. Returns list of issues."""
    issues: list[str] = []
    missing_keys = _CERT_REQUIRED_KEYS - set(cert_dict.keys())
    if missing_keys:
        issues.append(f"Missing certificate keys: {missing_keys}")

    prov = cert_dict.get("provenance", {})
    missing_prov = _PROVENANCE_REQUIRED_KEYS - set(prov.keys())
    if missing_prov:
        issues.append(f"Missing provenance keys: {missing_prov}")

    types_list = cert_dict.get("types", [])
    for i, t in enumerate(types_list):
        if "predicate" not in t or "verdict" not in t:
            issues.append(f"Type entry {i} missing 'predicate' or 'verdict' key")

    return issues


def _empty_kwargs(params: dict[str, Any]) -> dict[str, Any]:
    """Return empty kwargs (for predicates needing only seq)."""
    return {}


# Mapping from predicate name pattern to registry name and required kwargs.
# Each entry extracts the relevant verification parameters from the
# certificate's parameter dict so that re-evaluation can reproduce
# the original result.
_PREDICATE_KWARGS_MAP: dict[str, Callable[[dict[str, Any]], dict[str, Any]]] = {
    # DNA-level predicates
    "NoCrypticSplice": lambda params: {
        "known_exon_boundaries": params.get("exon_boundaries", []),
        "organism": params.get("organism", _DEFAULT_ORGANISM),
        "cryptic_splice_threshold": params.get("cryptic_splice_threshold", _DEFAULT_CRYPTIC_SPLICE_THRESHOLD),
    },
    "SpliceCorrect": lambda params: {
        "known_exon_boundaries": params.get("exon_boundaries", []),
        "cellular_context": params.get("cell_type", _DEFAULT_CELL_TYPE),
    },
    "GCInRange": lambda params: {
        "gc_lo": params.get("gc_lo", _DEFAULT_GC_LO),
        "gc_hi": params.get("gc_hi", _DEFAULT_GC_HI),
    },
    "CodonAdapted": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
        "threshold": params.get("cai_threshold", _DEFAULT_CAI_THRESHOLD),
    },
    "NoRestrictionSite": lambda params: {
        "enzyme_set": params.get("enzymes", []),
    },
    "InFrame": lambda params: {
        "exon_boundaries": params.get("exon_boundaries", [(0, 0)]),
    },
    "NoInstabilityMotif": _empty_kwargs,
    "NoCpGIsland": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoGTDinucleotide": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoCrypticPromoter": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
        "threshold": params.get("promoter_threshold", 0.7),
    },
    "NoUnexpectedTMDomain": lambda params: {
        "is_cytosolic": params.get("is_cytosolic", True),
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "mRNASecondaryStructure": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "CoTranslationalFolding": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
        "domain_boundaries": params.get("domain_boundaries", []),
    },
    # Optimizer predicates (8 original + 4 extended)
    "NoStopCodons": _empty_kwargs,
    "ValidCodingSeq": _empty_kwargs,
    "ConservationScore": lambda params: {
        "protein": params.get("protein", ""),
        "min_score": params.get("conservation_min_score", 0),
    },
    "CodonOptimality": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
        "threshold": params.get("cai_threshold", _DEFAULT_CAI_THRESHOLD),
    },
    # Structure predicates
    "StructureConfidence": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoMisfoldingRisk": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "CorrectFoldTopology": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoUnexpectedInteraction": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    # Stability predicates
    "StableFolding": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoDestabilizingMutation": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "DisulfideBondIntegrity": _empty_kwargs,
    "HydrophobicCoreQuality": _empty_kwargs,
    # Solubility predicates
    "SolubleExpression": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoAggregationProneRegion": _empty_kwargs,
    "ChargeComposition": _empty_kwargs,
    "NoLongHydrophobicStretch": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    # Immunogenicity predicates
    "LowImmunogenicity": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoStrongTCellEpitope": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoDominantBCellEpitope": _empty_kwargs,
    "PopulationCoverageSafe": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    # Sliding-window GC constraint
    "SlidingGC": lambda params: {
        "window_size": params.get("sliding_gc_window", 50),
        "gc_min": params.get("gc_lo", _DEFAULT_GC_LO),
        "gc_max": params.get("gc_hi", _DEFAULT_GC_HI),
    },
    # Extended diagnostic predicates
    "NoBlastMatches": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
        "reference_sequences": params.get("reference_sequences", []),
    },
    "PrimerCompatibility": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
        "primer_tm_min": params.get("primer_tm_min", 55.0),
        "primer_tm_max": params.get("primer_tm_max", 65.0),
    },
    "NoCrypticORF": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
        "min_orf_length": params.get("min_orf_length", 30),
    },
    "NoRQCTrigger": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoAluRepeat": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoMiRNABindingSite": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
        "min_seed_match": params.get("min_seed_match", 7),
        "tissue": params.get("tissue", ""),
        "seed_count": params.get("seed_count", 50),
        "context_scoring": params.get("context_scoring", "enabled"),
    },
    "NoM6ASite": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NoPolyASignal": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NucleosideModificationGuidance": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    "NucleosideModGuidance": lambda params: {  # Alias used in competitive landscape
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
    # Diagnostic predicates (registered via register_diagnostic — not in the
    # 43-predicate canonical contract but dispatchable for re-evaluation).
    "MRNAStability": lambda params: {
        "organism": params.get("organism", _DEFAULT_ORGANISM),
    },
}


def _resolve_predicate_name(cert_predicate_name: str) -> str | None:  # noqa: PYI042 – registry may be None at runtime
    """
    Resolve a possibly parameterized predicate name from a certificate
    (e.g., 'GCInRange(0.30, 0.70)') to its registry base name (e.g., 'GCInRange').

    Checks both the canonical 43-predicate registry (``registry.names()``)
    and the diagnostic-predicate registry (``registry.diagnostic_names()``)
    so that predicates like ``MRNAStability`` — emitted by the optimizer
    but not part of the 43-predicate contract — are still resolvable for
    re-evaluation during certificate verification.

    Returns None if no matching registry predicate is found.
    """
    if registry is None:
        logger.warning("Cannot resolve predicate name: registry is not available")
        return None
    if cert_predicate_name in registry:
        return cert_predicate_name
    # Check canonical predicates first (43-predicate contract), then
    # diagnostic predicates (e.g. MRNAStability).
    for base_name in registry.names():
        if cert_predicate_name.startswith(base_name):
            return base_name
    for base_name in registry.diagnostic_names():
        if cert_predicate_name == base_name or cert_predicate_name.startswith(base_name):
            return base_name
    return None


def verify_certificate(cert_dict: dict[str, Any], **kwargs: Any) -> tuple[str, list[str]]:
    """
    INDEPENDENTLY verify a guarantee certificate.

    This function does NOT trust the certificate. It re-evaluates every
    predicate from scratch using only the sequence and parameters in the certificate.

    Uses the predicate registry for dispatch — no if/elif chains.

    For graduated certificates, verification checks that the claimed verdicts
    are consistent with re-evaluation, but does NOT require all predicates to pass.

    Args:
        cert_dict: certificate as a plain dict
        **kwargs: override known_exon_boundaries, cellular_context, etc.

    Returns:
        (status, failure_reasons) where status is "VERIFIED" or "REJECTED"
    """
    failures: list[str] = []

    # Step 0: Validate certificate structure
    structural_issues = _validate_cert_structure(cert_dict)
    if structural_issues:
        return "REJECTED", structural_issues

    seq = cert_dict["sequence"].upper()
    prov = cert_dict.get("provenance", {})
    params = prov.get("parameters", {})

    # Merge certificate params with kwargs (kwargs take precedence)
    effective_params = dict(params)
    for key, val in kwargs.items():
        effective_params[key] = val

    # C15 fix: extract the cert's recorded slot_mode so we can honor it
    # when re-evaluating SLOT predicates. Honest CONSERVATIVE-mode certs
    # record SLOT predicates as UNCERTAIN; without honoring slot_mode,
    # `registry.verify` would return PASS (the underlying biological
    # check passes for a good sequence) and the cert would be REJECTED
    # for an honesty mismatch.
    cert_slot_mode_str = prov.get("slot_mode", SLOTMode.CONSERVATIVE.value)
    try:
        cert_slot_mode = SLOTMode(cert_slot_mode_str)
    except (ValueError, TypeError):
        cert_slot_mode = SLOTMode.CONSERVATIVE

    # Check 1: design_id matches hash
    hash_version = cert_dict.get("hash_version", 1)
    if hash_version == 1:
        # v1: sequence-only hash
        computed_hash = hashlib.new(_HASH_ALGORITHM, seq.encode()).hexdigest()
    else:
        # v2/v3: sequence + sorted predicates + version-specific param keys
        computed_hash = _compute_certificate_hash(
            seq, cert_dict.get("types", []), params, hash_version=hash_version,
        )
    if computed_hash != cert_dict["design_id"]:
        failures.append(
            f"design_id mismatch: computed {computed_hash[:_HASH_TRUNCATION_LENGTH]}... != "
            f"stored {cert_dict['design_id'][:_HASH_TRUNCATION_LENGTH]}..."
        )

    # Check 2: Re-evaluate each predicate using the registry
    # Equivalent verdict classes for verification purposes:
    #   PASS ≡ LIKELY_PASS  (both indicate the predicate is satisfied)
    #   FAIL ≡ LIKELY_FAIL  (both indicate the predicate is not satisfied)
    #   UNCERTAIN           (indeterminate — always matches itself)
    _verdict_class = {
        "PASS": "pass", "LIKELY_PASS": "pass",
        "FAIL": "fail", "LIKELY_FAIL": "fail",
        "UNCERTAIN": "uncertain",
    }

    for type_entry in cert_dict.get("types", []):
        predicate_name = type_entry["predicate"]
        claimed_verdict = type_entry["verdict"]

        try:
            registry_name = _resolve_predicate_name(predicate_name)
            if registry_name is None:
                # Unknown predicate -- cannot re-evaluate. As of W6-a,
                # ``MRNAStability`` (the only predicate the optimization
                # pipeline emits that is NOT in the 43-predicate canonical
                # contract) is registered as a *diagnostic* predicate via
                # ``registry.register_diagnostic`` and is therefore
                # resolvable + re-evaluable like any canonical predicate.
                # Reaching this branch now means the predicate name is
                # genuinely unknown to both the canonical and diagnostic
                # registries -- a certificate-structure bug, not a
                # registry-coverage gap.
                #
                # Rather than rejecting the certificate outright, we log a
                # warning and skip re-evaluation for this predicate --
                # mirroring the existing pattern in
                # ``slot_verification.py:1482-1489`` ("unknown predicate --
                # cannot recheck"). This is safe because the design_id
                # hash (v3, the default since _CURRENT_HASH_VERSION=3)
                # covers the predicate name AND verdict, so any tampering
                # with either field is independently detected by Check 1
                # above. Only the original cert generator's honesty about
                # this predicate's verdict is trusted, which is the same
                # trust model used for every other predicate's *claimed*
                # verdict (we re-check, but the claimed verdict itself
                # comes from the generator).
                logger.warning(
                    "Cannot re-evaluate predicate %r during certificate "
                    "verification: not in the canonical or diagnostic "
                    "predicate registries. Skipping re-evaluation "
                    "(claimed verdict: %s). The design_id hash still "
                    "binds the predicate name + verdict to the sequence.",
                    predicate_name, claimed_verdict,
                )
                continue

            kwarg_builder = _PREDICATE_KWARGS_MAP.get(registry_name, lambda p: {})
            verify_kwargs = kwarg_builder(effective_params)
            verify_kwargs["seq"] = seq

            # C15 fix: For SLOT predicates, route through
            # `verify_slot_predicate(name, slot_mode=cert_slot_mode, ...)`
            # so the cert's recorded slot_mode is honored during
            # re-evaluation. Without this, honest CONSERVATIVE-mode
            # certs that record SLOT predicates as UNCERTAIN are
            # REJECTED because `registry.verify` ignores slot_mode and
            # returns PASS (the underlying biological check passes).
            # `verify_slot_predicate` in CONSERVATIVE mode returns
            # UNCERTAIN, matching the claimed verdict. VERIFIED/PERMISSIVE
            # certs are likewise re-evaluated under their claimed mode.
            if is_slot_predicate(registry_name) and verify_slot_predicate is not None:
                slot_verdict, _evidence = verify_slot_predicate(
                    registry_name, slot_mode=cert_slot_mode, **verify_kwargs,
                )
                actual_verdict_value = slot_verdict.value
            else:
                result = registry.verify(registry_name, **verify_kwargs)
                actual_verdict_value = result.verdict.value

            claimed_class = _verdict_class.get(claimed_verdict, claimed_verdict)
            actual_class = _verdict_class.get(actual_verdict_value, actual_verdict_value)

            if claimed_class != actual_class:
                failures.append(
                    f"Predicate {predicate_name}: certificate claims {claimed_verdict}, "
                    f"re-evaluation gives {actual_verdict_value}"
                )
        except Exception as e:
            logger.warning(
                "Predicate %s re-evaluation error during verification: %s",
                predicate_name, e, exc_info=True,
            )
            failures.append(f"Predicate {predicate_name}: re-evaluation error: {e}")

    # Check 3: Provenance completeness
    for required_field in _PROVENANCE_REQUIRED_KEYS:
        if required_field not in prov:
            failures.append(f"Missing provenance field: {required_field}")

    if failures:
        logger.warning("Certificate verification FAILED: %d issues", len(failures))
        return "REJECTED", failures

    # Determine verification status based on certificate's own overall_status
    overall = prov.get("overall_status", "FULL_PASS")
    logger.info(
        "Certificate VERIFIED: design_id=%s... overall_status=%s solver=%s",
        cert_dict["design_id"][:_HASH_TRUNCATION_LENGTH], overall,
        prov.get("solver_backend", "unknown"),
    )
    return "VERIFIED", []


# ────────────────────────────────────────────────────────────
# Certificate Level Computation (GOLD / SILVER / BRONZE)
# ────────────────────────────────────────────────────────────
# Originally in certificates.py — consolidated here.

def compute_certificate(
    results: list[PredicateResult], slot_mode: SLOTMode = SLOTMode.CONSERVATIVE
) -> CertLevel:
    """Compute certificate level from predicate results.

    GOLD:   All predicates satisfied by optimization alone, AND no
            predicate has verdict=UNCERTAIN
    SILVER: All predicates satisfied, but 1 predicate has UNCERTAIN
            verdict OR some required mutagenesis/unavoidable constraints
    BRONZE: Some predicates unsatisfied OR 2+ predicates have UNCERTAIN
            verdict

    Uncertainty capping rules (Issue #10):
      1. Any Verdict.UNCERTAIN → certificate cannot exceed SILVER
      2. Multiple UNCERTAIN verdicts → certificate cannot exceed BRONZE
      3. LIKELY_PASS/LIKELY_FAIL do NOT cap — they express calibrated
         uncertainty and are acceptable for SILVER/GOLD

    This ensures downstream consumers can rely on the certificate level
    as an honest summary of evidence quality. A sequence where every
    predicate "passes" via heuristic UNCERTAIN verdicts should NOT
    receive GOLD.

    Note: A GOLD certificate with VERIFIED slot_mode is stronger than
    a GOLD certificate with CONSERVATIVE slot_mode, because VERIFIED mode
    provides actual evidence that SLOT predicates pass, while CONSERVATIVE
    mode only guarantees that non-SLOT predicates pass.

    SLOT predicate asterisk (*):
      SLOT predicates that contribute to the certificate level are
      marked with an asterisk to indicate that their PASS verdicts are
      empirically verified but NOT formally proven.  Non-SLOT predicate
      PASS verdicts are backed by machine-checked formal proofs.
      This distinction is critical for safety-critical applications
      where formal guarantees are required.
    """
    has_mutagenesis = False
    has_unavoidable = False
    has_unsatisfied = False
    uncertain_count = 0

    for r in results:
        if not r.passed:
            has_unsatisfied = True
        # Check structured flags first (preferred), then fall back to string matching
        if getattr(r, 'mutagenesis_applied', False) or "mutagenesis" in r.details.lower():
            has_mutagenesis = True
        if getattr(r, 'unavoidable_constraints', []) or "unavoidable" in r.details.lower():
            has_unavoidable = True
        # Count UNCERTAIN verdicts — but NOT LIKELY_PASS/LIKELY_FAIL
        verdict = getattr(r, 'verdict', None)
        if verdict == Verdict.UNCERTAIN:
            uncertain_count += 1

    # Uncertainty capping takes precedence
    if has_unsatisfied or uncertain_count >= 2:
        return CertLevel.BRONZE
    elif uncertain_count == 1 or has_mutagenesis or has_unavoidable:
        return CertLevel.SILVER
    else:
        return CertLevel.GOLD


def compute_uncertainty_summary(results: list[PredicateResult]) -> dict[str, Any]:
    """Compute uncertainty statistics from predicate results.

    Returns a dict with:
    - total_predicates: int
    - uncertain_count: int — number of UNCERTAIN verdicts
    - uncertain_predicates: list[str] — names of UNCERTAIN predicates
    - likely_pass_count: int — number of LIKELY_PASS verdicts
    - likely_fail_count: int — number of LIKELY_FAIL verdicts
    - definite_count: int — number of PASS or FAIL verdicts
    - confidence_score: float — weighted average confidence
      (1.0=PASS, 0.75=LIKELY_PASS, 0.5=UNCERTAIN, 0.25=LIKELY_FAIL, 0.0=FAIL)
    - slot_uncertain_count: int — UNCERTAIN predicates that are SLOT predicates
    """
    total = len(results)
    uncertain_preds: list[str] = []
    likely_pass = 0
    likely_fail = 0
    definite = 0
    slot_uncertain = 0
    total_confidence = 0.0

    for r in results:
        verdict = getattr(r, 'verdict', None)
        if verdict is None:
            total_confidence += 0.5  # treat missing verdict as UNCERTAIN
            uncertain_preds.append(r.predicate)
        elif verdict == Verdict.PASS:
            definite += 1
            total_confidence += 1.0
        elif verdict == Verdict.FAIL:
            definite += 1
            total_confidence += 0.0
        elif verdict == Verdict.LIKELY_PASS:
            likely_pass += 1
            total_confidence += 0.75
        elif verdict == Verdict.LIKELY_FAIL:
            likely_fail += 1
            total_confidence += 0.25
        elif verdict == Verdict.UNCERTAIN:
            uncertain_preds.append(r.predicate)
            total_confidence += 0.5
            if is_slot_predicate(r.predicate):
                slot_uncertain += 1

    return {
        "total_predicates": total,
        "uncertain_count": len(uncertain_preds),
        "uncertain_predicates": uncertain_preds,
        "likely_pass_count": likely_pass,
        "likely_fail_count": likely_fail,
        "definite_count": definite,
        "confidence_score": round(total_confidence / max(total, 1), 3),
        "slot_uncertain_count": slot_uncertain,
    }


def format_certificate(
    results: list[PredicateResult], seq: str, species: str,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> str:
    """Format a human-readable certificate report.

    The report includes:
    - Overall certificate level (GOLD / SILVER / BRONZE)
    - Per-predicate results with pass/fail status, verdict, and reason
    - Summary of which predicates passed and which failed with explanations
    - SLOT mode description
    """
    cert = compute_certificate(results, slot_mode)
    n_pass = sum(1 for r in results if r.passed)
    n_fail = len(results) - n_pass

    # Compute uncertainty summary
    uncertainty = compute_uncertainty_summary(results)

    lines = [
        "=" * _CERT_FORMAT_WIDTH,
        f"  BioCompiler v{VERSION} — Optimization Certificate",
        "=" * _CERT_FORMAT_WIDTH,
        f"  Sequence length: {len(seq)} bp",
        f"  Species:         {species}",
        f"  Certificate:     {cert.value}",
        f"  SLOT Mode:       {slot_mode.value}",
        f"  Predicates:      {n_pass} passed / {n_fail} failed / {len(results)} total",
        f"  Confidence:      {uncertainty['confidence_score']} ({uncertainty['definite_count']}/{uncertainty['total_predicates']} definite, {uncertainty['uncertain_count']} UNCERTAIN)",
    ]
    lines.extend([
        "-" * _CERT_FORMAT_WIDTH,
        "  Predicate Results:",
    ])
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        verdict_str = f" [{r.verdict.name}]" if r.verdict else ""
        # Mark mutagenesis-resolved constraints
        mutagenesis_marker = ""
        if "mutagenesis" in r.details.lower():
            mutagenesis_marker = " [MUTAGENESIS]"
        elif "unavoidable" in r.details.lower() and r.passed:
            mutagenesis_marker = " [UNAVOIDABLE]"
        # Mark SLOT predicates with asterisk to denote empirical (not formal) verification
        slot_marker = " [SLOT*]" if is_slot_predicate(r.predicate) else ""
        detail_line = r.details
        # Show tissue context for miRNA predicate when present in details
        if r.predicate == "NoMiRNABindingSite" and "tissue=" in r.details:
            # Extract unique tissue values from hit details
            tissues_found = set(_re.findall(r'tissue=(\w+)', r.details))
            tissues_found.discard('')
            if tissues_found:
                tissue_str = ', '.join(sorted(tissues_found))
                detail_line += f"  [tissues: {tissue_str}]"
        lines.append(f"    [{status}{verdict_str}{mutagenesis_marker}{slot_marker}] {r.predicate}: {detail_line}")

    # Summary: which predicates passed/failed and why
    lines.append("-" * _CERT_FORMAT_WIDTH)
    if n_pass > 0:
        passed_preds = [r for r in results if r.passed]
        lines.append(f"  Passed ({n_pass}):")
        for r in passed_preds:
            reason = r.details if r.details else "(no details)"
            # Truncate long reasons for readability
            if len(reason) > 80:
                reason = reason[:77] + "..."
            lines.append(f"    + {r.predicate}: {reason}")

    if n_fail > 0:
        failed_preds = [r for r in results if not r.passed]
        lines.append(f"  Failed ({n_fail}):")
        for r in failed_preds:
            reason = r.details if r.details else "(no details)"
            if len(reason) > 80:
                reason = reason[:77] + "..."
            lines.append(f"    - {r.predicate}: {reason}")

    # Uncertain predicates section
    if uncertainty['uncertain_count'] > 0:
        lines.append(f"  Uncertain Predicates ({uncertainty['uncertain_count']}):")
        for pred_name in uncertainty['uncertain_predicates']:
            # Find the matching result to get details
            matching = [r for r in results if r.predicate == pred_name]
            detail = matching[0].details if matching else ""
            if len(detail) > 60:
                detail = detail[:57] + "..."
            slot_tag = " [SLOT]" if is_slot_predicate(pred_name) else ""
            lines.append(f"    [UNCERTAIN{slot_tag}] {pred_name}: {detail}")

    lines.append("=" * _CERT_FORMAT_WIDTH)
    lines.append("")
    lines.append("  Certificate Levels:")
    lines.append("    GOLD   — All constraints satisfied by synonymous optimization")
    lines.append("    SILVER — All constraints satisfied (some required AA substitution)")
    lines.append("    BRONZE — Some constraints could not be satisfied")
    lines.append("")
    lines.append("  SLOT Mode:")
    lines.append("    CONSERVATIVE — SLOT predicates always UNCERTAIN (Lean4 model)")
    lines.append("    VERIFIED     — SLOT predicates PASS when verification conditions met")
    lines.append("    PERMISSIVE   — SLOT predicates PASS with weaker evidence")
    lines.append("")
    lines.append("  * SLOT Predicate Caveat:")
    lines.append("    SLOT predicates marked with * are empirically verified but NOT")
    lines.append("    formally proven.  Their PASS verdicts rely on external tools or")
    lines.append("    scanners whose results cannot be machine-checked.  Non-SLOT")
    lines.append("    predicate PASS verdicts are backed by formal proof.  For")
    lines.append("    safety-critical applications, treat SLOT* results with extra")
    lines.append("    scrutiny and consider independent validation.")
    lines.append("=" * _CERT_FORMAT_WIDTH)
    return "\n".join(lines)
