"""Type stubs for biocompiler.certificate — Certificate generation and verification."""

from __future__ import annotations

from typing import Any

from .types import Certificate, TypeCheckResult, SLOTMode
from .type_system import CertLevel, PredicateResult


# ────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────

VERSION: str


# ────────────────────────────────────────────────────────────
# Core functions
# ────────────────────────────────────────────────────────────

def generate_certificate(
    sequence: str,
    type_results: list[TypeCheckResult],
    input_params: dict[str, Any],
    require_all_pass: bool = ...,
    mutagenesis_substitutions: list[dict[str, Any]] | None = ...,
    slot_mode: SLOTMode = ...,
    solver_backend: str | None = ...,
    solver_config: dict[str, Any] | None = ...,
) -> Certificate:
    """Generate a machine-checkable guarantee certificate.

    Args:
        sequence: the DNA sequence being certified
        type_results: list of TypeCheckResult objects
        input_params: dict of parameters used for verification
        require_all_pass: if True, raise on any failure
        mutagenesis_substitutions: AA substitutions applied
        slot_mode: SLOT verification mode
        solver_backend: name of the solver backend used
        solver_config: solver configuration dict

    Returns:
        Certificate object with all predicate results documented

    Raises:
        ValueError: if sequence is empty or type_results is empty
        CertificateGenerationError: if require_all_pass=True and any predicate failed
    """
    ...


def verify_certificate(cert_dict: dict[str, Any], **kwargs: Any) -> tuple[str, list[str]]:
    """Independently verify a guarantee certificate.

    Args:
        cert_dict: certificate as a plain dict
        **kwargs: override verification parameters

    Returns:
        (status, failure_reasons) where status is "VERIFIED" or "REJECTED"
    """
    ...


def compute_certificate(
    results: list[PredicateResult],
    slot_mode: SLOTMode = ...,
) -> CertLevel:
    """Compute certificate level (GOLD/SILVER/BRONZE) from predicate results."""
    ...


def format_certificate(
    results: list[PredicateResult],
    seq: str,
    species: str,
    slot_mode: SLOTMode = ...,
) -> str:
    """Format a human-readable certificate report."""
    ...
