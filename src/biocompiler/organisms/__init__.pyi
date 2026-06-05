"""Type stubs for biocompiler.organisms — public API surface."""

from __future__ import annotations

from typing import TypedDict


# ────────────────────────────────────────────────────────────
# TypedDict for SPECIES entries
# ────────────────────────────────────────────────────────────

class SpeciesEntry(TypedDict):
    """Structure for per-species codon optimization data."""
    cai_weights: dict[str, float]
    codon_usage_validation: bool


# ────────────────────────────────────────────────────────────
# Registry tables
# ────────────────────────────────────────────────────────────

CODON_USAGE_TABLES: dict[str, dict[str, tuple[str, float, float, int]]]
CODON_ADAPTIVENESS_TABLES: dict[str, dict[str, float]]
PREFERRED_CODON_TABLES: dict[str, dict[str, str]]
SUPPORTED_ORGANISMS: list[str]
ORGANISM_GC_TARGETS: dict[str, tuple[float, float]]

# Sharp-Li CAI reference set
SHARP_LI_ADAPTIVENESS_TABLES: dict[str, dict[str, float]]
SHARP_LI_CODON_USAGE: dict[str, dict[str, float]]
SHARP_LI_CAI_WEIGHTS: dict[str, dict[str, float]]
SHARP_LI_REFERENCE_GENES: dict[str, dict[str, dict[str, float]]]
SHARP_LI_PUBLISHED_CAI: dict[str, dict[str, float]]

# Organism name aliases
HUMAN: str
E_COLI: str
MOUSE: str
CHO: str
YEAST: str

# Organism name resolution maps
ORGANISM_ALIASES: dict[str, str]
SPECIES_SHORT_NAMES: dict[str, str]

# Legacy backward-compatible API (DEPRECATED)
SPECIES: dict[str, SpeciesEntry]
ECOLI_CAI: dict[str, float]
HUMAN_CAI: dict[str, float]
MOUSE_CAI: dict[str, float]
CHO_CAI: dict[str, float]
YEAST_CAI: dict[str, float]


# ────────────────────────────────────────────────────────────
# Public functions
# ────────────────────────────────────────────────────────────

def resolve_organism(
    name: str,
    *,
    default: str | None = ...,
    strict: bool = ...,
) -> str:
    """Resolve any organism name, alias, or short key to the canonical name.

    Args:
        name: Any organism identifier.
        default: Value to return when name cannot be resolved.
        strict: When True, raise an error instead of falling back.

    Returns:
        The canonical organism name (e.g. 'Escherichia_coli').

    Raises:
        UnsupportedOrganismError: If strict is True and name is not a known alias.
    """
    ...


def get_species_cai_weights(species_key: str) -> dict[str, float]:
    """Return CAI weights for a species. DEPRECATED: use CODON_ADAPTIVENESS_TABLES directly."""
    ...


def compute_cai_weights(usage: dict[str, float]) -> dict[str, float]:
    """Compute CAI weights from codon usage. Most frequent codon per AA = 1.0."""
    ...


def get_sharp_li_adaptiveness_tables() -> dict[str, dict[str, float]]:
    """Return the Sharp-Li adaptiveness tables, building them lazily."""
    ...


def validate_cai_tables() -> list[str]:
    """Verify all CAI tables are internally consistent.

    Returns:
        List of error description strings. Empty list means all tables are valid.
    """
    ...
