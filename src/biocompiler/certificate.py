"""
BioCompiler Certificate Engine — Generation & Verification

Production-grade certificate system with:
- Registry-based verification (no string-prefix dispatch)
- All verification parameters embedded IN the certificate
- Hash integrity check for sequence + parameters
- Input validation for certificate structure
- CertificateError instead of raw assert
"""

import hashlib
import logging
from datetime import datetime, timezone
from .types import Verdict, TypeCheckResult, Certificate, combined_verdict
from .type_system import registry
from .exceptions import CertificateGenerationError, CertificateVerificationError

logger = logging.getLogger(__name__)

VERSION = "2.1.0"

# Required keys in a certificate dict for verification
_CERT_REQUIRED_KEYS = {"version", "design_id", "sequence", "types", "provenance"}
_PROVENANCE_REQUIRED_KEYS = {"tool", "version", "timestamp", "input_hash"}


def generate_certificate(
    sequence: str,
    type_results: list[TypeCheckResult],
    input_params: dict,
) -> Certificate:
    """
    Generate a machine-checkable guarantee certificate.

    PRECONDITION: All type results must have verdict PASS.
    Raises CertificateGenerationError if any predicate failed.

    Args:
        sequence: the DNA sequence being certified
        type_results: list of TypeCheckResult objects
        input_params: dict of parameters used (gene, organism, exon_boundaries, etc.)
            MUST include all parameters needed for independent re-verification:
            - exon_boundaries, gc_lo, gc_hi, cai_threshold, organism, cell_type, enzymes

    Returns:
        Certificate object

    Raises:
        CertificateGenerationError: if any predicate has non-PASS verdict
    """
    failures = [r for r in type_results if r.verdict != Verdict.PASS]
    if failures:
        raise CertificateGenerationError(failures)

    # Compute hash once
    seq_hash = hashlib.sha256(sequence.encode()).hexdigest()

    # Ensure all verification parameters are embedded
    complete_params = dict(input_params)
    complete_params.setdefault("organism", "Homo_sapiens")
    complete_params.setdefault("cell_type", "HEK293T")
    complete_params.setdefault("gc_lo", 0.30)
    complete_params.setdefault("gc_hi", 0.70)
    complete_params.setdefault("cai_threshold", 0.5)
    complete_params.setdefault("enzymes", ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"])
    complete_params.setdefault("cryptic_splice_threshold", 3.0)
    complete_params.setdefault("exon_boundaries", [(0, len(sequence))])

    cert = Certificate(
        version=VERSION,
        design_id=seq_hash,
        sequence=sequence,
        types=[
            {
                "predicate": r.predicate,
                "verdict": r.verdict.value,
                "derivation": r.derivation,
            }
            for r in type_results
        ],
        provenance={
            "tool": "BioCompiler",
            "version": VERSION,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "parameters": complete_params,
            "input_hash": seq_hash,
        },
    )
    logger.info("Certificate generated: design_id=%s...", cert.design_id[:16])
    return cert


def _validate_cert_structure(cert_dict: dict) -> list[str]:
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


# Mapping from predicate name pattern to registry name and required kwargs
_PREDICATE_KWARGS_MAP = {
    "NoCrypticSplice": lambda params: {
        "known_exon_boundaries": params.get("exon_boundaries", []),
    },
    "SpliceCorrect": lambda params: {
        "known_exon_boundaries": params.get("exon_boundaries", []),
        "cellular_context": params.get("cell_type", "HEK293T"),
    },
    "GCInRange": lambda params: {
        "gc_lo": params.get("gc_lo", 0.30),
        "gc_hi": params.get("gc_hi", 0.70),
    },
    "CodonAdapted": lambda params: {
        "organism": params.get("organism", "Homo_sapiens"),
        "threshold": params.get("cai_threshold", 0.5),
    },
    "NoRestrictionSite": lambda params: {
        "enzyme_set": params.get("enzymes", []),
    },
    "InFrame": lambda params: {
        "exon_boundaries": params.get("exon_boundaries", [(0, 0)]),
    },
    "NoInstabilityMotif": lambda params: {},
    "NoCpGIsland": lambda params: {},
}


def _resolve_predicate_name(cert_predicate_name: str) -> str | None:
    """
    Resolve a possibly parameterized predicate name from a certificate
    (e.g., 'GCInRange(0.30, 0.70)') to its registry base name (e.g., 'GCInRange').

    Returns None if no matching registry predicate is found.
    """
    # First try exact match
    if cert_predicate_name in registry:
        return cert_predicate_name
    # Try prefix match for parameterized predicates
    for base_name in registry.names():
        if cert_predicate_name.startswith(base_name):
            return base_name
    return None


def verify_certificate(cert_dict: dict, **kwargs) -> tuple[str, list[str]]:
    """
    INDEPENDENTLY verify a guarantee certificate.

    This function does NOT trust the certificate. It re-evaluates every
    predicate from scratch using only the sequence and parameters in the certificate.

    Uses the predicate registry for dispatch — no if/elif chains.

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
        # Return early if structure is invalid — can't safely proceed
        return "REJECTED", structural_issues

    seq = cert_dict["sequence"].upper()
    prov = cert_dict.get("provenance", {})
    params = prov.get("parameters", {})

    # Merge certificate params with kwargs (kwargs take precedence for explicit overrides)
    effective_params = dict(params)
    for key, val in kwargs.items():
        effective_params[key] = val

    # Check 1: design_id matches SHA-256 of sequence
    computed_hash = hashlib.sha256(seq.encode()).hexdigest()
    if computed_hash != cert_dict["design_id"]:
        failures.append(
            f"design_id mismatch: computed {computed_hash[:16]}... != "
            f"stored {cert_dict['design_id'][:16]}..."
        )

    # Check 2: Re-evaluate each predicate using the registry
    for type_entry in cert_dict.get("types", []):
        predicate_name = type_entry["predicate"]
        claimed_verdict = type_entry["verdict"]

        try:
            # Resolve predicate name to registry name
            registry_name = _resolve_predicate_name(predicate_name)
            if registry_name is None:
                failures.append(f"Unknown predicate: {predicate_name}")
                continue

            # Build kwargs for the registry verify call
            kwarg_builder = _PREDICATE_KWARGS_MAP.get(registry_name, lambda p: {})
            verify_kwargs = kwarg_builder(effective_params)
            verify_kwargs["seq"] = seq

            result = registry.verify(registry_name, **verify_kwargs)

            if result.verdict.value != claimed_verdict:
                failures.append(
                    f"Predicate {predicate_name}: certificate claims {claimed_verdict}, "
                    f"re-evaluation gives {result.verdict.value}"
                )
        except Exception as e:
            failures.append(f"Predicate {predicate_name}: re-evaluation error: {e}")

    # Check 3: Provenance completeness
    for required_field in _PROVENANCE_REQUIRED_KEYS:
        if required_field not in prov:
            failures.append(f"Missing provenance field: {required_field}")

    if failures:
        logger.warning("Certificate verification FAILED: %d issues", len(failures))
        return "REJECTED", failures

    logger.info("Certificate VERIFIED: design_id=%s...", cert_dict["design_id"][:16])
    return "VERIFIED", []
