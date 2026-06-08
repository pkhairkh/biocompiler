"""Deprecated: use biocompiler.provenance.crypto instead."""
import warnings

warnings.warn(
    "biocompiler.provenance_crypto is deprecated — use biocompiler.provenance.crypto instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.provenance.crypto import *  # noqa: F401,F403

__all__ = [
    "sign_record",
    "verify_record",
    "ProvenanceIntegrityError",
    "get_provenance_secret",
    "_canonical_json",
    "is_mandatory_provenance",
]
