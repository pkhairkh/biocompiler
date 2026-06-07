"""Deprecated: use biocompiler.immunogenicity.netmhcpan instead."""
import warnings

warnings.warn(
    "biocompiler.netmhcpan is deprecated — use biocompiler.immunogenicity.netmhcpan instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.immunogenicity.netmhcpan import *  # noqa: F401,F403

__all__ = [
    "_rank_to_binding_score",
    "_predict_binding_local",
    "_extract_error_message",
    "_extract_result_url",
    "_get_default_cache",
    "_get_default_client",
    "_is_mhc_ii_allele",
    "_looks_like_netmhcpan_output",
    "_parse_data_line",
]
