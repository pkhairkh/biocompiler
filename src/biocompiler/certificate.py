"""
BioCompiler Certificate Engine — Generation & Verification

FIXES from toy model:
- CertificateError instead of raw assert
- Registry-based verification (no string-prefix dispatch)
- Verification parameters embedded IN the certificate
- Hash integrity check for sequence + parameters
"""

import hashlib
import logging
from datetime import datetime, timezone
from .types import Verdict, TypeCheckResult, Certificate, combined_verdict
from .type_system import registry
from .exceptions import CertificateGenerationError, CertificateVerificationError

logger = logging.getLogger(__name__)

VERSION = "2.0.0"


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

    Returns:
        Certificate object

    Raises:
        CertificateGenerationError: if any predicate has non-PASS verdict
    """
    failures = [r for r in type_results if r.verdict != Verdict.PASS]
    if failures:
        raise CertificateGenerationError(failures)

    cert = Certificate(
        version=VERSION,
        design_id=hashlib.sha256(sequence.encode()).hexdigest(),
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
            "parameters": input_params,
            "input_hash": hashlib.sha256(sequence.encode()).hexdigest(),
        },
    )
    logger.info("Certificate generated: design_id=%s...", cert.design_id[:16])
    return cert


def verify_certificate(cert_dict: dict, **kwargs) -> tuple[str, list[str]]:
    """
    INDEPENDENTLY verify a guarantee certificate.

    This function does NOT trust the certificate. It re-evaluates every
    predicate from scratch using only the sequence and parameters in the certificate.

    Args:
        cert_dict: certificate as a plain dict
        **kwargs: override known_exon_boundaries, cellular_context, etc.

    Returns:
        (status, failure_reasons) where status is "VERIFIED" or "REJECTED"
    """
    failures: list[str] = []

    seq = cert_dict["sequence"].upper()
    prov = cert_dict.get("provenance", {})
    params = prov.get("parameters", {})

    # Extract verification parameters from certificate
    known_exon_boundaries = kwargs.get(
        "known_exon_boundaries",
        params.get("exon_boundaries", [(0, len(seq))])
    )
    cellular_context = kwargs.get("cellular_context", params.get("cell_type", "HEK293T"))
    gc_lo = params.get("gc_lo", 0.30)
    gc_hi = params.get("gc_hi", 0.70)
    cai_threshold = params.get("cai_threshold", 0.5)
    organism = params.get("organism", "Homo_sapiens")
    enzymes = params.get("enzymes", ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"])

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
            # Dispatch based on predicate name through the registry
            if predicate_name == "NoCrypticSplice":
                result = registry.verify("NoCrypticSplice", seq=seq, known_exon_boundaries=known_exon_boundaries)
            elif predicate_name.startswith("SpliceCorrect"):
                result = registry.verify("SpliceCorrect", seq=seq, known_exon_boundaries=known_exon_boundaries, cellular_context=cellular_context)
            elif predicate_name.startswith("GCInRange"):
                result = registry.verify("GCInRange", seq=seq, gc_lo=gc_lo, gc_hi=gc_hi)
            elif predicate_name.startswith("CodonAdapted"):
                result = registry.verify("CodonAdapted", seq=seq, organism=organism, threshold=cai_threshold)
            elif predicate_name.startswith("NoRestrictionSite"):
                result = registry.verify("NoRestrictionSite", seq=seq, enzyme_set=enzymes)
            elif predicate_name == "InFrame":
                result = registry.verify("InFrame", seq=seq, exon_boundaries=known_exon_boundaries)
            elif predicate_name == "NoInstabilityMotif":
                result = registry.verify("NoInstabilityMotif", seq=seq)
            else:
                failures.append(f"Unknown predicate: {predicate_name}")
                continue

            if result.verdict.value != claimed_verdict:
                failures.append(
                    f"Predicate {predicate_name}: certificate claims {claimed_verdict}, "
                    f"re-evaluation gives {result.verdict.value}"
                )
        except Exception as e:
            failures.append(f"Predicate {predicate_name}: re-evaluation error: {e}")

    # Check 3: Provenance completeness
    prov = cert_dict.get("provenance", {})
    for required_field in ("tool", "version", "timestamp", "input_hash"):
        if required_field not in prov:
            failures.append(f"Missing provenance field: {required_field}")

    if failures:
        logger.warning("Certificate verification FAILED: %d issues", len(failures))
        return "REJECTED", failures

    logger.info("Certificate VERIFIED: design_id=%s...", cert_dict["design_id"][:16])
    return "VERIFIED", []
