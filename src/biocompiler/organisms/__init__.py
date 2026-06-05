"""
BioCompiler Organism Data

Multi-organism codon usage tables loaded from structured data.
Extensible — add a new organism by adding a new module.

Supported organisms:
- Homo_sapiens (Human)
- Escherichia_coli (E. coli)
- Mus_musculus (Mouse)
- CHO_K1 (Chinese Hamster Ovary)
- Saccharomyces_cerevisiae (Yeast)

Migrated from species.py (v9.2.0):
- compute_cai_weights, ECOLI_CODON_USAGE, ECOLI_CAI
- HUMAN_CODON_USAGE_SIMPLE, HUMAN_CAI, SPECIES

v10.0.0 changes (BREAKING):
- CAI table unification: SPECIES['ecoli']['cai_weights'] now uses
  E_COLI_CODON_ADAPTIVENESS (from CODON_ADAPTIVENESS_TABLES) instead of
  the legacy ECOLI_CODON_USAGE. This ensures the optimizer and evaluator
  agree on optimal codons. CAI values will differ from v9.x (now correct).
- E. coli optimal codons corrected for 5 amino acids (Phe, Ile, Tyr, His, Arg)
  where per_thousand data previously identified the wrong codon as optimal.
- ORGANISM_ALIASES / SPECIES_SHORT_NAMES: centralized name resolution maps.
- resolve_organism(): single source of truth for organism name resolution.
- ORGANISM_GC_TARGETS now uses per-organism (gc_lo, gc_hi) ranges
- SPECIES entries include codon_usage_validation flag
- All re-exports are explicitly typed

v10.1.0 changes (Task 27 — organisms cleanup):
- CODON_ADAPTIVENESS_TABLES is now the single source of truth for CAI weights.
- SPECIES dict is DEPRECATED; derived FROM CODON_ADAPTIVENESS_TABLES.
- get_species_cai_weights() is DEPRECATED; use CODON_ADAPTIVENESS_TABLES
  with resolve_organism() instead.
- ECOLI_CAI, HUMAN_CAI, MOUSE_CAI, CHO_CAI, YEAST_CAI are DEPRECATED;
  use CODON_ADAPTIVENESS_TABLES directly.
- validate_cai_tables() verifies SPECIES consistency with CODON_ADAPTIVENESS_TABLES.
"""

from __future__ import annotations

import warnings
from typing import Any, TypedDict

from .human import HUMAN_CODON_USAGE as HUMAN_CODON_USAGE
from .human import HUMAN_CODON_ADAPTIVENESS as HUMAN_CODON_ADAPTIVENESS
from .human import HUMAN_PREFERRED_CODONS as HUMAN_PREFERRED_CODONS
from .human import HUMAN_CODON_USAGE_SIMPLE as HUMAN_CODON_USAGE_SIMPLE
from .human import HUMAN_CODON_PAIR_BIAS as HUMAN_CODON_PAIR_BIAS
from .human import HUMAN_EXPRESSION_OPTIMIZATION_PARAMS as HUMAN_EXPRESSION_OPTIMIZATION_PARAMS
from .human import HUMAN_UTR_MODELS as HUMAN_UTR_MODELS
from .e_coli import E_COLI_CODON_USAGE as E_COLI_CODON_USAGE
from .e_coli import E_COLI_CODON_ADAPTIVENESS as E_COLI_CODON_ADAPTIVENESS
from .e_coli import E_COLI_PREFERRED_CODONS as E_COLI_PREFERRED_CODONS
from .e_coli import E_COLI_CODON_PAIR_BIAS as E_COLI_CODON_PAIR_BIAS
from .e_coli import E_COLI_EXPRESSION_OPTIMIZATION_PARAMS as E_COLI_EXPRESSION_OPTIMIZATION_PARAMS
from .e_coli import ECOLI_CODON_USAGE as ECOLI_CODON_USAGE
from .e_coli import compute_codon_pair_bias as compute_codon_pair_bias
from .mouse import MOUSE_CODON_USAGE as MOUSE_CODON_USAGE
from .mouse import MOUSE_CODON_ADAPTIVENESS as MOUSE_CODON_ADAPTIVENESS
from .mouse import MOUSE_PREFERRED_CODONS as MOUSE_PREFERRED_CODONS
from .mouse import MOUSE_CODON_PAIR_BIAS as MOUSE_CODON_PAIR_BIAS
from .mouse import MOUSE_EXPRESSION_OPTIMIZATION_PARAMS as MOUSE_EXPRESSION_OPTIMIZATION_PARAMS
from .mouse import MOUSE_UTR_MODELS as MOUSE_UTR_MODELS
from .cho import CHO_CODON_USAGE as CHO_CODON_USAGE
from .cho import CHO_CODON_ADAPTIVENESS as CHO_CODON_ADAPTIVENESS
from .cho import CHO_PREFERRED_CODONS as CHO_PREFERRED_CODONS
from .yeast import YEAST_CODON_USAGE as YEAST_CODON_USAGE
from .yeast import YEAST_CODON_ADAPTIVENESS as YEAST_CODON_ADAPTIVENESS
from .yeast import YEAST_PREFERRED_CODONS as YEAST_PREFERRED_CODONS
from ..organism_db import OrganismDatabase as OrganismDatabase
from ..organism_db import get_database as get_database

__all__: list[str] = [
    # TypedDict for SPECIES entries
    "SpeciesEntry",
    # Registry tables
    "CODON_USAGE_TABLES",
    "CODON_ADAPTIVENESS_TABLES",
    "PREFERRED_CODON_TABLES",
    "SUPPORTED_ORGANISMS",
    "ORGANISM_GC_TARGETS",
    # Sharp-Li CAI reference set
    "SHARP_LI_ADAPTIVENESS_TABLES",
    "SHARP_LI_CODON_USAGE",
    "SHARP_LI_CAI_WEIGHTS",
    "SHARP_LI_REFERENCE_GENES",
    "SHARP_LI_PUBLISHED_CAI",
    "get_sharp_li_cai_weights",
    "compute_cai_with_reference",
    # Organism name aliases
    "HUMAN",
    "E_COLI",
    "MOUSE",
    "CHO",
    "YEAST",
    # Organism name resolution (species ↔ organism)
    "resolve_organism",
    "ORGANISM_ALIASES",
    "SPECIES_SHORT_NAMES",
    # Legacy backward-compatible API (DEPRECATED)
    "SPECIES",
    "get_species_cai_weights",
    "ECOLI_CAI",
    "HUMAN_CAI",
    "MOUSE_CAI",
    "CHO_CAI",
    "YEAST_CAI",
    "compute_cai_weights",
    # Validation
    "validate_cai_tables",
    # Per-organism re-exports
    "HUMAN_CODON_USAGE",
    "HUMAN_CODON_ADAPTIVENESS",
    "HUMAN_PREFERRED_CODONS",
    "HUMAN_CODON_USAGE_SIMPLE",
    "HUMAN_CODON_PAIR_BIAS",
    "HUMAN_EXPRESSION_OPTIMIZATION_PARAMS",
    "HUMAN_UTR_MODELS",
    "E_COLI_CODON_USAGE",
    "E_COLI_CODON_ADAPTIVENESS",
    "E_COLI_PREFERRED_CODONS",
    "E_COLI_CODON_PAIR_BIAS",
    "E_COLI_EXPRESSION_OPTIMIZATION_PARAMS",
    "ECOLI_CODON_USAGE",
    "compute_codon_pair_bias",
    "MOUSE_CODON_USAGE",
    "MOUSE_CODON_ADAPTIVENESS",
    "MOUSE_PREFERRED_CODONS",
    "MOUSE_CODON_PAIR_BIAS",
    "MOUSE_EXPRESSION_OPTIMIZATION_PARAMS",
    "MOUSE_UTR_MODELS",
    "CHO_CODON_USAGE",
    "CHO_CODON_ADAPTIVENESS",
    "CHO_PREFERRED_CODONS",
    "YEAST_CODON_USAGE",
    "YEAST_CODON_ADAPTIVENESS",
    "YEAST_PREFERRED_CODONS",
    # External re-exports
    "OrganismDatabase",
    "get_database",
]

# ────────────────────────────────────────────────────────────
# Explicit type annotations for re-exported names
# These serve as documentation; the actual types come from the
# imports above.  mypy may flag these as redefinitions, so we
# suppress that check.
# ────────────────────────────────────────────────────────────
HUMAN_CODON_USAGE: dict[str, tuple[str, float, float, int]]  # type: ignore[no-redef]
HUMAN_CODON_ADAPTIVENESS: dict[str, float]  # type: ignore[no-redef]
HUMAN_PREFERRED_CODONS: dict[str, str]  # type: ignore[no-redef]
HUMAN_CODON_USAGE_SIMPLE: dict[str, float]  # type: ignore[no-redef]
HUMAN_CODON_PAIR_BIAS: dict[str, float]  # type: ignore[no-redef]
HUMAN_EXPRESSION_OPTIMIZATION_PARAMS: dict[str, object]  # type: ignore[no-redef]
HUMAN_UTR_MODELS: dict[str, object]  # type: ignore[no-redef]
E_COLI_CODON_USAGE: dict[str, tuple[str, float, float, int]]  # type: ignore[no-redef]
E_COLI_CODON_ADAPTIVENESS: dict[str, float]  # type: ignore[no-redef]
E_COLI_PREFERRED_CODONS: dict[str, str]  # type: ignore[no-redef]
E_COLI_CODON_PAIR_BIAS: dict[str, float]  # type: ignore[no-redef]
E_COLI_EXPRESSION_OPTIMIZATION_PARAMS: dict[str, object]  # type: ignore[no-redef]
ECOLI_CODON_USAGE: dict[str, float]  # type: ignore[no-redef]
MOUSE_CODON_USAGE: dict[str, tuple[str, float, float, int]]  # type: ignore[no-redef]
MOUSE_CODON_ADAPTIVENESS: dict[str, float]  # type: ignore[no-redef]
MOUSE_PREFERRED_CODONS: dict[str, str]  # type: ignore[no-redef]
MOUSE_CODON_PAIR_BIAS: dict[str, float]  # type: ignore[no-redef]
MOUSE_EXPRESSION_OPTIMIZATION_PARAMS: dict[str, object]  # type: ignore[no-redef]
MOUSE_UTR_MODELS: dict[str, object]  # type: ignore[no-redef]
CHO_CODON_USAGE: dict[str, tuple[str, float, float, int]]  # type: ignore[no-redef]
CHO_CODON_ADAPTIVENESS: dict[str, float]  # type: ignore[no-redef]
CHO_PREFERRED_CODONS: dict[str, str]  # type: ignore[no-redef]
YEAST_CODON_USAGE: dict[str, tuple[str, float, float, int]]  # type: ignore[no-redef]
YEAST_CODON_ADAPTIVENESS: dict[str, float]  # type: ignore[no-redef]
YEAST_PREFERRED_CODONS: dict[str, str]  # type: ignore[no-redef]

# ════════════════════════════════════════════════════════════
# SINGLE SOURCE OF TRUTH: CODON_ADAPTIVENESS_TABLES
#
# All CAI weights for codon optimization and evaluation are
# derived from CODON_ADAPTIVENESS_TABLES.  The legacy SPECIES
# dict and per-organism _CAI aliases are DEPRECATED and kept
# only for backward compatibility — they are derived FROM
# CODON_ADAPTIVENESS_TABLES, never independently maintained.
# ════════════════════════════════════════════════════════════

# Registry: organism name -> codon data
CODON_USAGE_TABLES: dict[str, dict[str, tuple[str, float, float, int]]] = {
    "Homo_sapiens": HUMAN_CODON_USAGE,
    "human": HUMAN_CODON_USAGE,  # alias for resolve_organism convenience
    "Escherichia_coli": E_COLI_CODON_USAGE,
    "e_coli": E_COLI_CODON_USAGE,  # alias for resolve_organism convenience
    "Mus_musculus": MOUSE_CODON_USAGE,
    "mouse": MOUSE_CODON_USAGE,  # alias for resolve_organism convenience
    "CHO_K1": CHO_CODON_USAGE,
    "cho": CHO_CODON_USAGE,  # alias for resolve_organism convenience
    "Saccharomyces_cerevisiae": YEAST_CODON_USAGE,
    "yeast": YEAST_CODON_USAGE,  # alias for resolve_organism convenience
}

CODON_ADAPTIVENESS_TABLES: dict[str, dict[str, float]] = {
    "Homo_sapiens": HUMAN_CODON_ADAPTIVENESS,
    "human": HUMAN_CODON_ADAPTIVENESS,  # alias for resolve_organism convenience
    "Escherichia_coli": E_COLI_CODON_ADAPTIVENESS,
    "e_coli": E_COLI_CODON_ADAPTIVENESS,  # alias for resolve_organism convenience
    "Mus_musculus": MOUSE_CODON_ADAPTIVENESS,
    "mouse": MOUSE_CODON_ADAPTIVENESS,  # alias for resolve_organism convenience
    "CHO_K1": CHO_CODON_ADAPTIVENESS,
    "cho": CHO_CODON_ADAPTIVENESS,  # alias for resolve_organism convenience
    "Saccharomyces_cerevisiae": YEAST_CODON_ADAPTIVENESS,
    "yeast": YEAST_CODON_ADAPTIVENESS,  # alias for resolve_organism convenience
}

PREFERRED_CODON_TABLES: dict[str, dict[str, str]] = {
    "Homo_sapiens": HUMAN_PREFERRED_CODONS,
    "human": HUMAN_PREFERRED_CODONS,  # alias for resolve_organism convenience
    "Escherichia_coli": E_COLI_PREFERRED_CODONS,
    "e_coli": E_COLI_PREFERRED_CODONS,  # alias for resolve_organism convenience
    "Mus_musculus": MOUSE_PREFERRED_CODONS,
    "mouse": MOUSE_PREFERRED_CODONS,  # alias for resolve_organism convenience
    "CHO_K1": CHO_PREFERRED_CODONS,
    "cho": CHO_PREFERRED_CODONS,  # alias for resolve_organism convenience
    "Saccharomyces_cerevisiae": YEAST_PREFERRED_CODONS,
    "yeast": YEAST_PREFERRED_CODONS,  # alias for resolve_organism convenience
}

SUPPORTED_ORGANISMS: list[str] = list(CODON_USAGE_TABLES.keys())

# Organism-specific GC content target ranges (gc_lo, gc_hi)
# Source: genome-wide coding-sequence composition studies
# These are aspirational targets, not hard constraints — the optimizer
# will only nudge toward them when GC is already in range.
# Previously hardcoded to 0.30/0.70 for all organisms; now uses
# organism-specific ranges derived from observed coding GC distributions.
ORGANISM_GC_TARGETS: dict[str, tuple[float, float]] = {
    "Escherichia_coli": (0.45, 0.55),            # E. coli K-12 coding GC ~50.8%
    "Homo_sapiens": (0.40, 0.60),                # Human coding GC ~52.3%
    "Saccharomyces_cerevisiae": (0.35, 0.45),    # Yeast coding GC ~38.3%
    "Mus_musculus": (0.40, 0.55),                # Mouse coding GC ~42.0%
    "CHO_K1": (0.40, 0.60),                      # Chinese Hamster coding GC ~44%
}

# ────────────────────────────────────────────────────────────
# Sharp-Li (1987) CAI reference adaptiveness tables
#
# These tables use the codon frequency data from the original
# Sharp & Li (1987) highly-expressed gene reference sets,
# sourced via the cai_validated module.  They differ from the
# Kazusa-derived CODON_ADAPTIVENESS_TABLES because the
# underlying reference gene composition is different:
#
#   - Kazusa: broader collection of high-expression genes from
#     the Kazusa Codon Usage Database
#   - Sharp-Li: the specific 24 highly-expressed E. coli genes
#     (or equivalent subsets for other organisms) used in the
#     original CAI paper
#
# The differences are most pronounced for E. coli, where the
# Sharp-Li reference set yields CAI values closer to those
# reported in the 1987 paper.
# ────────────────────────────────────────────────────────────

SHARP_LI_ADAPTIVENESS_TABLES: dict[str, dict[str, float]] = {}
# Lazily populated on first access to avoid circular import at module level.
# The actual computation is deferred until SHARP_LI_ADAPTIVENESS_TABLES
# is accessed through get_sharp_li_adaptiveness_tables().


def get_sharp_li_adaptiveness_tables() -> dict[str, dict[str, float]]:
    """Return the Sharp-Li adaptiveness tables, building them lazily.

    Computes the tables on first call and caches the result in the
    module-level SHARP_LI_ADAPTIVENESS_TABLES dict.

    Returns:
        Dict mapping organism name to {codon: adaptiveness} dicts.
    """
    global SHARP_LI_ADAPTIVENESS_TABLES
    if SHARP_LI_ADAPTIVENESS_TABLES:
        return SHARP_LI_ADAPTIVENESS_TABLES

    from ..benchmarking.cai_validated import _REFERENCE_SETS

    tables: dict[str, dict[str, float]] = {}
    for organism, ref in _REFERENCE_SETS.items():
        adaptiveness: dict[str, float] = {}
        for _aa, codon_freqs in ref.items():
            if not codon_freqs:
                continue
            max_freq = max(codon_freqs.values())
            for codon, freq in codon_freqs.items():
                w = freq / max_freq if max_freq > 0 else 0.0
                # Floor at 0.01 as per Sharp & Li recommendation
                adaptiveness[codon] = max(w, 0.01)
        tables[organism] = adaptiveness

    # Add ALL aliases from ORGANISM_ALIASES to ensure consistency
    # with CODON_ADAPTIVENESS_TABLES.  This replaces the old
    # hard-coded _ALIASES dict so that every name recognised by
    # resolve_organism() also works as a key here.
    for alias, canonical in ORGANISM_ALIASES.items():
        if canonical in tables and alias not in tables:
            tables[alias] = tables[canonical]

    SHARP_LI_ADAPTIVENESS_TABLES = tables
    return tables


# Aliases for backward compat
HUMAN: str = "Homo_sapiens"
E_COLI: str = "Escherichia_coli"
MOUSE: str = "Mus_musculus"
CHO: str = "CHO_K1"
YEAST: str = "Saccharomyces_cerevisiae"

# ────────────────────────────────────────────────────────────
# Organism name resolution: species ↔ organism aliases
# ────────────────────────────────────────────────────────────

# Mapping of all known organism name aliases to the canonical name
# used as a key in CODON_ADAPTIVENESS_TABLES, CODON_USAGE_TABLES, etc.
#
# This centralises the alias resolution so that every module in
# BioCompiler uses the same mapping, instead of each file
# maintaining its own partial copy.
#
# Alias categories:
#   1. Short informal names:  'ecoli', 'human', 'mouse', 'cho', 'yeast'
#   2. Abbreviated binomials:  'E_coli', 'e_coli', 'H_sapiens', 'h_sapiens',
#                               'M_musculus', 'm_musculus', 'S_cerevisiae',
#                               's_cerevisiae'
#   3. Display names:          'E. coli', 'H. sapiens', 'M. musculus',
#                               'S. cerevisiae'
#   4. Canonical binomials:    'Escherichia_coli', 'Homo_sapiens', etc.
#   5. Legacy synonyms:        'CHO' (maps to CHO_K1)
ORGANISM_ALIASES: dict[str, str] = {
    # ── Escherichia coli ──
    "ecoli": "Escherichia_coli",
    "E_coli": "Escherichia_coli",
    "e_coli": "Escherichia_coli",
    "E. coli": "Escherichia_coli",
    "E. Coli": "Escherichia_coli",
    "e. coli": "Escherichia_coli",
    "Escherichia_coli": "Escherichia_coli",
    "Escherichia coli": "Escherichia_coli",
    # ── Homo sapiens ──
    "human": "Homo_sapiens",
    "H_sapiens": "Homo_sapiens",
    "h_sapiens": "Homo_sapiens",
    "H. sapiens": "Homo_sapiens",
    "H. Sapiens": "Homo_sapiens",
    "h. sapiens": "Homo_sapiens",
    "Homo_sapiens": "Homo_sapiens",
    "Homo sapiens": "Homo_sapiens",
    # ── Mus musculus ──
    "mouse": "Mus_musculus",
    "M_musculus": "Mus_musculus",
    "m_musculus": "Mus_musculus",
    "M. musculus": "Mus_musculus",
    "M. Musculus": "Mus_musculus",
    "m. musculus": "Mus_musculus",
    "Mus_musculus": "Mus_musculus",
    "Mus musculus": "Mus_musculus",
    # ── CHO (Chinese Hamster Ovary) ──
    "cho": "CHO_K1",
    "CHO": "CHO_K1",
    "CHO_K1": "CHO_K1",
    "Cricetulus_griseus": "CHO_K1",
    # ── Saccharomyces cerevisiae ──
    "yeast": "Saccharomyces_cerevisiae",
    "S_cerevisiae": "Saccharomyces_cerevisiae",
    "s_cerevisiae": "Saccharomyces_cerevisiae",
    "S. cerevisiae": "Saccharomyces_cerevisiae",
    "S. Cerevisiae": "Saccharomyces_cerevisiae",
    "s. cerevisiae": "Saccharomyces_cerevisiae",
    "Saccharomyces_cerevisiae": "Saccharomyces_cerevisiae",
    "Saccharomyces cerevisiae": "Saccharomyces_cerevisiae",
    "saccharomyces_cerevisiae": "Saccharomyces_cerevisiae",
    "saccharomyces cerevisiae": "Saccharomyces_cerevisiae",
}

# ────────────────────────────────────────────────────────────
# Post-construction: ensure ALL aliases from ORGANISM_ALIASES
# are available as keys in the three registry tables.
#
# The dict literals above include the 5 primary short-name aliases
# (human, e_coli, mouse, cho, yeast) for discoverability, but
# ORGANISM_ALIASES defines many more (ecoli, E_coli, H_sapiens,
# E. coli, etc.).  This loop adds any remaining aliases so that
# every name in ORGANISM_ALIASES works as a direct lookup key.
# ────────────────────────────────────────────────────────────
for _alias, _canonical in ORGANISM_ALIASES.items():
    if _alias not in CODON_ADAPTIVENESS_TABLES and _canonical in CODON_ADAPTIVENESS_TABLES:
        CODON_ADAPTIVENESS_TABLES[_alias] = CODON_ADAPTIVENESS_TABLES[_canonical]
    if _alias not in CODON_USAGE_TABLES and _canonical in CODON_USAGE_TABLES:
        CODON_USAGE_TABLES[_alias] = CODON_USAGE_TABLES[_canonical]
    if _alias not in PREFERRED_CODON_TABLES and _canonical in PREFERRED_CODON_TABLES:
        PREFERRED_CODON_TABLES[_alias] = PREFERRED_CODON_TABLES[_canonical]

# Reverse mapping: canonical organism name → primary short species key.
# Used to map from organism names back to the SPECIES dict keys.
SPECIES_SHORT_NAMES: dict[str, str] = {
    "Escherichia_coli": "ecoli",
    "Homo_sapiens": "human",
    "Mus_musculus": "mouse",
    "CHO_K1": "cho",
    "Saccharomyces_cerevisiae": "yeast",
}


def resolve_organism(
    name: str,
    *,
    default: str | None = None,
    strict: bool = False,
) -> str:
    """Resolve any organism name, alias, or short key to the canonical name.

    This is the single source of truth for organism name resolution
    across BioCompiler.  All public APIs (``optimize_sequence``,
    ``compute_cai``, REST endpoints, etc.) should call this function
    instead of maintaining their own partial alias maps.

    Accepted input forms:

    +-------------------------+------------------------+
    | Example input           | Resolved output        |
    +=========================+========================+
    | ``'ecoli'``             | ``'Escherichia_coli'`` |
    | ``'E. coli'``           | ``'Escherichia_coli'`` |
    | ``'e_coli'``            | ``'Escherichia_coli'`` |
    | ``'Escherichia_coli'``  | ``'Escherichia_coli'`` |
    | ``'human'``             | ``'Homo_sapiens'``     |
    | ``'H. sapiens'``        | ``'Homo_sapiens'``     |
    | ``'h_sapiens'``         | ``'Homo_sapiens'``     |
    | ``'Homo_sapiens'``      | ``'Homo_sapiens'``     |
    +-------------------------+------------------------+

    Args:
        name: Any organism identifier — short key, abbreviated
            binomial, display name with periods, or full canonical
            binomial with underscores.
        default: Value to return when *name* cannot be resolved.
            If ``None`` (the default) and *strict* is ``False``,
            the function returns ``name`` unchanged as a fallback.
            If *strict* is ``True`` and no *default* is given, an
            :class:`~.exceptions.UnsupportedOrganismError` is raised.
        strict: When ``True``, raise an error instead of falling
            back to the unresolved *name*.

    Returns:
        The canonical organism name (e.g. ``'Escherichia_coli'``).

    Raises:
        UnsupportedOrganismError: If *strict* is ``True`` and
            *name* is not a known alias.
    """
    if not name or not name.strip():
        if default is not None:
            return default
        if strict:
            from ..exceptions import UnsupportedOrganismError
            raise UnsupportedOrganismError(
                name or "", list(ORGANISM_ALIASES.keys())
            )
        return "Homo_sapiens"  # sensible default

    canonical = ORGANISM_ALIASES.get(name)
    if canonical is not None:
        return canonical

    # If already a canonical name (present in CODON_ADAPTIVENESS_TABLES
    # but not in the alias dict for some reason), return as-is.
    if name in CODON_ADAPTIVENESS_TABLES:
        return name

    # Unresolved name
    if default is not None:
        return default
    if strict:
        from ..exceptions import UnsupportedOrganismError
        raise UnsupportedOrganismError(name, list(ORGANISM_ALIASES.keys()))

    # Non-strict: return as-is (backward compat — caller's responsibility)
    return name

# ────────────────────────────────────────────────────────────
# Migrated from species.py — backward-compatible API
# ────────────────────────────────────────────────────────────

# Pseudocount for missing codons in CAI weight computation.
# Prevents zero-frequency codons from getting a weight of 0.0,
# which would make the geometric mean CAI zero for any sequence
# containing that codon.
_DEFAULT_MISSING_CODON_FREQ: float = 0.1


def compute_cai_weights(usage: dict[str, float]) -> dict[str, float]:
    """Compute CAI weights from codon usage. Most frequent codon per AA = 1.0.

    Args:
        usage: Dict mapping codon strings to per-thousand frequency values.

    Returns:
        Dict mapping codon strings to CAI weight (0.0–1.0).
    """
    from ..type_system import AA_TO_CODONS
    weights: dict[str, float] = {}
    for aa, codons in AA_TO_CODONS.items():
        if aa == "*":
            continue
        freqs = [usage.get(c, _DEFAULT_MISSING_CODON_FREQ) for c in codons]
        max_freq = max(freqs) if freqs else 1.0
        for codon, freq in zip(codons, freqs):
            weights[codon] = freq / max_freq if max_freq > 0 else 0.0
    return weights


# ────────────────────────────────────────────────────────────
# DEPRECATED: Per-organism CAI aliases
#
# These are kept for backward compatibility only.
# Use CODON_ADAPTIVENESS_TABLES directly instead.
#
# Previously ECOLI_CAI used the legacy ECOLI_CODON_USAGE (different source)
# and HUMAN_CAI used HUMAN_CODON_USAGE_SIMPLE, which could produce
# different optimal codon rankings than CODON_ADAPTIVENESS_TABLES.
# Now they all point to CODON_ADAPTIVENESS_TABLES entries.
# ────────────────────────────────────────────────────────────
ECOLI_CAI: dict[str, float] = E_COLI_CODON_ADAPTIVENESS  # deprecated

HUMAN_CAI: dict[str, float] = HUMAN_CODON_ADAPTIVENESS  # deprecated

MOUSE_CAI: dict[str, float] = MOUSE_CODON_ADAPTIVENESS  # deprecated

CHO_CAI: dict[str, float] = CHO_CODON_ADAPTIVENESS  # deprecated

YEAST_CAI: dict[str, float] = YEAST_CODON_ADAPTIVENESS  # deprecated


# TypedDict for SPECIES entries — provides both CAI weights and
# a codon_usage_validation flag indicating whether the organism's
# codon usage data has been validated against reference databases.
class SpeciesEntry(TypedDict):
    """Structure for per-species codon optimization data.

    Attributes:
        cai_weights: Codon Adaptation Index weights (0.0–1.0), mapping
            codon strings to their relative adaptiveness values.
        codon_usage_validation: Whether the codon usage table for this
            organism has been validated against reference databases
            (e.g., Kazusa Codon Usage Database, CoCoPUTs).
    """
    cai_weights: dict[str, float]
    codon_usage_validation: bool


# ────────────────────────────────────────────────────────────
# DEPRECATED: SPECIES dict
#
# SPECIES is a legacy structure used by the optimizer for codon selection.
# Now that we've unified on CODON_ADAPTIVENESS_TABLES, SPECIES is kept
# for backward compatibility only.  It is derived FROM
# CODON_ADAPTIVENESS_TABLES, not independently maintained.
#
# New code should use CODON_ADAPTIVENESS_TABLES directly:
#   from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
#   weights = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
#
# To resolve organism names, use resolve_organism():
#   from biocompiler.organisms import resolve_organism, CODON_ADAPTIVENESS_TABLES
#   canonical = resolve_organism("ecoli")   # → "Escherichia_coli"
#   weights = CODON_ADAPTIVENESS_TABLES[canonical]
# ────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────
# DEPRECATED: _DeprecatedSpeciesDict
#
# A dict subclass that wraps the legacy SPECIES data but emits
# DeprecationWarning on __getitem__ access.  New code should
# use CODON_ADAPTIVENESS_TABLES directly with resolve_organism().
#
# The nesting trap: SPECIES["ecoli"] returns a SpeciesEntry
# (dict with "cai_weights" and "codon_usage_validation" keys),
# NOT a flat codon→weight dict.  Any code that does
# SPECIES["ecoli"]["ATG"] will get a KeyError, and
# SPECIES["ecoli"].get("ATG", 0.0) will return 0.0.
# Use get_species_cai_weights() or CODON_ADAPTIVENESS_TABLES instead.
# ────────────────────────────────────────────────────────────

class _DeprecatedSpeciesDict(dict):  # type: ignore[misc]
    """Dict subclass that warns on item access to encourage migration.

    The SPECIES dict has a nesting issue: SPECIES[key] returns a
    SpeciesEntry (with 'cai_weights' and 'codon_usage_validation'
    keys), not a flat codon→weight dict.  This wrapper emits a
    DeprecationWarning on every ``__getitem__`` call to steer
    callers toward CODON_ADAPTIVENESS_TABLES or
    get_species_cai_weights().
    """

    _warned_keys: set[str]  # Track keys we've already warned about

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._warned_keys: set[str] = set()

    def __getitem__(self, key: str) -> SpeciesEntry:  # type: ignore[override]
        if key not in self._warned_keys:
            warnings.warn(
                f"SPECIES[{key!r}] is deprecated — use "
                f"CODON_ADAPTIVENESS_TABLES with resolve_organism() "
                f"instead.  Note: SPECIES[{key!r}] returns a nested "
                f"SpeciesEntry, not a flat codon→weight dict.",
                DeprecationWarning,
                stacklevel=2,
            )
            self._warned_keys.add(key)
        return super().__getitem__(key)


_SPECIES_RAW: dict[str, SpeciesEntry] = {
    short_key: {
        "cai_weights": dict(CODON_ADAPTIVENESS_TABLES[canonical_name]),
        "codon_usage_validation": True,
    }
    for canonical_name, short_key in SPECIES_SHORT_NAMES.items()
}

SPECIES: dict[str, SpeciesEntry] = _DeprecatedSpeciesDict(_SPECIES_RAW)


def get_species_cai_weights(species_key: str) -> dict[str, float]:
    """Return CAI weights for a species.  DEPRECATED: use CODON_ADAPTIVENESS_TABLES directly.

    This function is kept for backward compatibility only.  New code
    should use CODON_ADAPTIVENESS_TABLES with resolve_organism()::

        from biocompiler.organisms import resolve_organism, CODON_ADAPTIVENESS_TABLES
        organism = resolve_organism("ecoli")
        weights = CODON_ADAPTIVENESS_TABLES[organism]

    Args:
        species_key: Short species name used as a SPECIES dict key
            (e.g. ``"ecoli"``, ``"human"``).  Also accepts canonical
            organism names (e.g. ``"Escherichia_coli"``) and other
            aliases recognised by resolve_organism().

    Returns:
        Dict mapping codon strings to their CAI weight values.
    """
    warnings.warn(
        "get_species_cai_weights is deprecated — use CODON_ADAPTIVENESS_TABLES "
        "directly with resolve_organism()",
        DeprecationWarning,
        stacklevel=2,
    )
    organism = resolve_organism(species_key)
    if organism in CODON_ADAPTIVENESS_TABLES:
        return dict(CODON_ADAPTIVENESS_TABLES[organism])
    # Fallback to ecoli for unknown organisms
    return dict(CODON_ADAPTIVENESS_TABLES["Escherichia_coli"])

# ────────────────────────────────────────────────────────────
# Sharp & Li (1987) CAI reference set (E. coli)
# ────────────────────────────────────────────────────────────
# Import from the dedicated sharp_li_reference module which provides the
# original 24-gene E. coli reference set from Sharp & Li (1987).
# This produces CAI values directly comparable to the published paper
# (e.g. lacZ≈0.27, trpA≈0.84, recA≈0.76).
#
# The import is placed after compute_cai_weights is defined above
# to avoid circular import issues (sharp_li_reference imports
# compute_cai_weights from this module).

from .sharp_li_reference import SHARP_LI_REFERENCE_GENES as SHARP_LI_REFERENCE_GENES
from .sharp_li_reference import SHARP_LI_CODON_USAGE as SHARP_LI_CODON_USAGE
from .sharp_li_reference import SHARP_LI_CAI_WEIGHTS as SHARP_LI_CAI_WEIGHTS
from .sharp_li_reference import SHARP_LI_PUBLISHED_CAI as SHARP_LI_PUBLISHED_CAI
from .sharp_li_reference import get_sharp_li_cai_weights as get_sharp_li_cai_weights
from .sharp_li_reference import compute_cai_with_reference as compute_cai_with_reference


# ────────────────────────────────────────────────────────────
# CAI Table Validation
#
# Validates that all CODON_ADAPTIVENESS_TABLES entries are
# internally consistent: for each amino acid with multiple
# synonymous codons, exactly one codon should have w = 1.0
# (the most frequent), and all others should have w < 1.0.
# Also validates that the SPECIES dict is consistent with
# CODON_ADAPTIVENESS_TABLES (since SPECIES is now derived
# from it, this should always pass, but the check is kept
# as a safety net).
# ────────────────────────────────────────────────────────────

def validate_cai_tables() -> list[str]:
    """Verify all CAI tables are internally consistent.

    For each organism in CODON_ADAPTIVENESS_TABLES, checks that:
    1. For every amino acid with multiple synonymous codons,
       exactly one codon has adaptiveness w = 1.0.
    2. All other synonymous codons have w < 1.0.
    3. No codon has w > 1.0 or w < 0.
    4. The CODON_USAGE_TABLES fractions are consistent with
       the per-thousand values (same ranking within each AA).
    5. The SPECIES dict (derived from CODON_ADAPTIVENESS_TABLES)
       is consistent with its source tables.

    Returns:
        List of error description strings.  An empty list means
        all tables are valid.
    """
    from ..type_system import AA_TO_CODONS

    errors: list[str] = []

    # 1. Validate adaptiveness tables
    for org_name, adaptiveness in CODON_ADAPTIVENESS_TABLES.items():
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*":
                continue  # skip stop codons
            if len(codons) == 1:
                continue  # Met, Trp — only one codon

            # Check that exactly one codon has w=1.0
            optimal = [c for c in codons if adaptiveness.get(c, 0) == 1.0]
            if len(optimal) != 1:
                errors.append(
                    f"{org_name} {aa}: expected 1 optimal codon with w=1.0, "
                    f"got {len(optimal)}: {optimal}"
                )

            # Check no codon exceeds 1.0 or is negative
            for c in codons:
                w = adaptiveness.get(c, 0)
                if w > 1.0:
                    errors.append(
                        f"{org_name} {aa} {c}: w={w} > 1.0"
                    )
                if w < 0:
                    errors.append(
                        f"{org_name} {aa} {c}: w={w} < 0"
                    )

    # 2. Validate CODON_USAGE_TABLES internal consistency
    for org_name, usage in CODON_USAGE_TABLES.items():
        # Group codons by amino acid
        aa_codons: dict[str, list[tuple[str, float, float]]] = {}
        for codon, (aa, frac, per_thousand, _count) in usage.items():
            if aa == "*":
                continue
            aa_codons.setdefault(aa, []).append((codon, frac, per_thousand))

        for aa, codon_data in aa_codons.items():
            if len(codon_data) <= 1:
                continue

            # Check that the codon with the highest fraction also has
            # the highest per-thousand value
            by_frac = sorted(codon_data, key=lambda x: x[1], reverse=True)
            by_pt = sorted(codon_data, key=lambda x: x[2], reverse=True)

            if by_frac[0][0] != by_pt[0][0]:
                errors.append(
                    f"{org_name} {aa}: fraction says {by_frac[0][0]} is optimal "
                    f"(frac={by_frac[0][1]}) but per_thousand says "
                    f"{by_pt[0][0]} is optimal (pt={by_pt[0][2]})"
                )

    # 3. Validate SPECIES dict consistency with CODON_ADAPTIVENESS_TABLES
    #    (Since SPECIES is now derived FROM CODON_ADAPTIVENESS_TABLES,
    #     this should always pass, but we keep it as a safety net.)
    #    SPECIES_SHORT_NAMES maps canonical_name → short_key
    #    Use _SPECIES_RAW to avoid triggering deprecation warnings.
    for canonical_name, short_key in SPECIES_SHORT_NAMES.items():
        if short_key not in _SPECIES_RAW:
            errors.append(f"SPECIES missing key '{short_key}'")
            continue
        if canonical_name not in CODON_ADAPTIVENESS_TABLES:
            errors.append(f"CODON_ADAPTIVENESS_TABLES missing '{canonical_name}'")
            continue

        species_weights = _SPECIES_RAW[short_key]["cai_weights"]
        adapt_weights = CODON_ADAPTIVENESS_TABLES[canonical_name]

        # Check that the optimal codon (w=1.0) for each AA matches
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*" or len(codons) <= 1:
                continue
            species_optimal = [c for c in codons if species_weights.get(c, 0) == 1.0]
            adapt_optimal = [c for c in codons if adapt_weights.get(c, 0) == 1.0]
            if species_optimal != adapt_optimal:
                errors.append(
                    f"SPECIES['{short_key}'] {aa}: optimal={species_optimal} "
                    f"but CODON_ADAPTIVENESS_TABLES['{canonical_name}'] "
                    f"optimal={adapt_optimal}"
                )

    return errors
