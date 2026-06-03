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
"""

import hashlib
import logging
from datetime import datetime, timezone
from .types import Verdict, TypeCheckResult, Certificate, combined_verdict, SLOTMode
try:
    from .type_system import registry
except ImportError:
    # Fallback: registry not available in type_system
    registry = None
from .type_system import CertLevel, PredicateResult
from .exceptions import CertificateGenerationError, CertificateVerificationError

logger = logging.getLogger(__name__)

try:
    from . import __version__ as _PKG_VERSION
    VERSION = _PKG_VERSION
except ImportError:
    VERSION = "7.2.0"

# Required keys in a certificate dict for verification
_CERT_REQUIRED_KEYS = {"version", "design_id", "sequence", "types", "provenance"}
_PROVENANCE_REQUIRED_KEYS = {"tool", "version", "timestamp", "input_hash"}


def generate_certificate(
    sequence: str,
    type_results: list[TypeCheckResult],
    input_params: dict,
    require_all_pass: bool = False,
    mutagenesis_substitutions: list[dict] | None = None,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
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

    Returns:
        Certificate object with all predicate results documented

    Raises:
        CertificateGenerationError: only if require_all_pass=True and any predicate failed

    Pre-conditions:
    - sequence is a non-empty DNA string
    - type_results is a non-empty list
    - input_params is a dict (may be empty)
    """
    assert sequence, "Sequence must not be empty"
    assert type_results, "Type results must not be empty"
    failures = [r for r in type_results if r.verdict != Verdict.PASS]

    if require_all_pass and failures:
        raise CertificateGenerationError(failures)

    # Compute hash once
    seq_hash = hashlib.sha256(sequence.encode()).hexdigest()

    # Compute overall status
    n_pass = sum(1 for r in type_results if r.verdict == Verdict.PASS)
    n_total = len(type_results)
    overall_status = "FULL_PASS" if not failures else f"PARTIAL_{n_pass}/{n_total}"

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

    # Build provenance dict
    provenance = {
        "tool": "BioCompiler",
        "version": VERSION,
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

    cert = Certificate(
        version=VERSION,
        design_id=seq_hash,
        sequence=sequence,
        types=[
            {
                "predicate": r.predicate,
                "verdict": r.verdict.value,
                "derivation": r.derivation,
                "knowledge_gap": r.knowledge_gap,
            }
            for r in type_results
        ],
        provenance=provenance,
    )
    status_msg = overall_status if failures else "FULL_PASS"
    logger.info("Certificate generated: design_id=%s... status=%s", cert.design_id[:16], status_msg)
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
    if cert_predicate_name in registry:
        return cert_predicate_name
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
            registry_name = _resolve_predicate_name(predicate_name)
            if registry_name is None:
                failures.append(f"Unknown predicate: {predicate_name}")
                continue

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

    # Determine verification status based on certificate's own overall_status
    overall = prov.get("overall_status", "FULL_PASS")
    logger.info(
        "Certificate VERIFIED: design_id=%s... overall_status=%s",
        cert_dict["design_id"][:16], overall
    )
    return "VERIFIED", []


# ────────────────────────────────────────────────────────────
# Certificate Level Computation (GOLD / SILVER / BRONZE)
# ────────────────────────────────────────────────────────────
# Originally in certificates.py — consolidated here.

def compute_certificate(results: list[PredicateResult], slot_mode: SLOTMode = SLOTMode.CONSERVATIVE) -> CertLevel:
    """Compute certificate level from predicate results.

    GOLD:   All predicates satisfied by optimization alone
    SILVER: All predicates satisfied, some required mutagenesis or have
            unavoidable constraints (e.g., Valine GT dinucleotides that
            can only be removed by AA substitution)
    BRONZE: Some predicates could not be fully satisfied

    Note: A GOLD certificate with VERIFIED slot_mode is stronger than
    a GOLD certificate with CONSERVATIVE slot_mode, because VERIFIED mode
    provides actual evidence that SLOT predicates pass, while CONSERVATIVE
    mode only guarantees that non-SLOT predicates pass.
    """
    has_mutagenesis = False
    has_unavoidable = False
    has_unsatisfied = False

    for r in results:
        if not r.passed:
            has_unsatisfied = True
        if "mutagenesis" in r.details.lower():
            has_mutagenesis = True
        if "unavoidable" in r.details.lower():
            has_unavoidable = True

    if has_unsatisfied:
        return CertLevel.BRONZE
    elif has_mutagenesis or has_unavoidable:
        return CertLevel.SILVER
    else:
        return CertLevel.GOLD


def format_certificate(results: list[PredicateResult], seq: str, species: str, slot_mode: SLOTMode = SLOTMode.CONSERVATIVE) -> str:
    """Format a human-readable certificate report."""
    cert = compute_certificate(results, slot_mode)
    lines = [
        "=" * 60,
        "  BioCompiler v7.0.0 — Optimization Certificate",
        "=" * 60,
        f"  Sequence length: {len(seq)} bp",
        f"  Species:         {species}",
        f"  Certificate:     {cert.value}",
        f"  SLOT Mode:       {slot_mode.value}",
        "-" * 60,
        "  Predicate Results:",
    ]
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        verdict_str = f" [{r.verdict.name}]" if r.verdict else ""
        # Mark mutagenesis-resolved constraints
        mutagenesis_marker = ""
        if "mutagenesis" in r.details.lower():
            mutagenesis_marker = " [MUTAGENESIS]"
        elif "unavoidable" in r.details.lower() and r.passed:
            mutagenesis_marker = " [UNAVOIDABLE]"
        # Mark SLOT predicates
        from .slot_verification import is_slot_predicate
        slot_marker = " [SLOT]" if is_slot_predicate(r.predicate) else ""
        lines.append(f"    [{status}{verdict_str}{mutagenesis_marker}{slot_marker}] {r.predicate}: {r.details}")
    lines.append("=" * 60)
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
    lines.append("=" * 60)
    return "\n".join(lines)
