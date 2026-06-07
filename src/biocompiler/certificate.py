"""Deprecated: use biocompiler.provenance.certificate instead."""
import warnings

warnings.warn(
    "biocompiler.certificate is deprecated — use biocompiler.provenance.certificate instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.provenance.certificate import *  # noqa: F401,F403

__all__ = [
    "generate_certificate",
    "verify_certificate",
    "compute_certificate",
    "format_certificate",
    "VERSION",
    "_CERTIFICATE_VERSION",
    "_REQUIRED_INPUT_PARAM_KEYS",
    "_compute_certificate_hash",
    "_compute_gc_content",
    "_CURRENT_HASH_VERSION",
    "_HASH_ALGORITHM",
    "_V2_HASH_PARAM_KEYS",
    "CertLevel",
    "_validate_cert_structure",
    "_CERT_REQUIRED_KEYS",
    "_PROVENANCE_REQUIRED_KEYS",
]
