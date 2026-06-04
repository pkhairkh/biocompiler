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

v9.2.0 changes:
- ORGANISM_GC_TARGETS now uses per-organism (gc_lo, gc_hi) ranges
- SPECIES entries include codon_usage_validation flag
- All re-exports are explicitly typed
"""

from __future__ import annotations

from typing import TypedDict

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
    # Legacy backward-compatible API
    "SPECIES",
    "ECOLI_CAI",
    "HUMAN_CAI",
    "MOUSE_CAI",
    "CHO_CAI",
    "YEAST_CAI",
    "compute_cai_weights",
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
# ────────────────────────────────────────────────────────────
HUMAN_CODON_USAGE: dict[str, tuple[str, float, float, int]]
HUMAN_CODON_ADAPTIVENESS: dict[str, float]
HUMAN_PREFERRED_CODONS: dict[str, str]
HUMAN_CODON_USAGE_SIMPLE: dict[str, float]
HUMAN_CODON_PAIR_BIAS: dict[str, float]
HUMAN_EXPRESSION_OPTIMIZATION_PARAMS: dict[str, object]
HUMAN_UTR_MODELS: dict[str, object]
E_COLI_CODON_USAGE: dict[str, tuple[str, float, float, int]]
E_COLI_CODON_ADAPTIVENESS: dict[str, float]
E_COLI_PREFERRED_CODONS: dict[str, str]
E_COLI_CODON_PAIR_BIAS: dict[str, float]
E_COLI_EXPRESSION_OPTIMIZATION_PARAMS: dict[str, object]
ECOLI_CODON_USAGE: dict[str, float]
MOUSE_CODON_USAGE: dict[str, tuple[str, float, float, int]]
MOUSE_CODON_ADAPTIVENESS: dict[str, float]
MOUSE_PREFERRED_CODONS: dict[str, str]
MOUSE_CODON_PAIR_BIAS: dict[str, float]
MOUSE_EXPRESSION_OPTIMIZATION_PARAMS: dict[str, object]
MOUSE_UTR_MODELS: dict[str, object]
CHO_CODON_USAGE: dict[str, tuple[str, float, float, int]]
CHO_CODON_ADAPTIVENESS: dict[str, float]
CHO_PREFERRED_CODONS: dict[str, str]
YEAST_CODON_USAGE: dict[str, tuple[str, float, float, int]]
YEAST_CODON_ADAPTIVENESS: dict[str, float]
YEAST_PREFERRED_CODONS: dict[str, str]

# Registry: organism name -> codon data
CODON_USAGE_TABLES: dict[str, dict[str, tuple[str, float, float, int]]] = {
    "Homo_sapiens": HUMAN_CODON_USAGE,
    "Escherichia_coli": E_COLI_CODON_USAGE,
    "Mus_musculus": MOUSE_CODON_USAGE,
    "CHO_K1": CHO_CODON_USAGE,
    "Saccharomyces_cerevisiae": YEAST_CODON_USAGE,
}

CODON_ADAPTIVENESS_TABLES: dict[str, dict[str, float]] = {
    "Homo_sapiens": HUMAN_CODON_ADAPTIVENESS,
    "Escherichia_coli": E_COLI_CODON_ADAPTIVENESS,
    "Mus_musculus": MOUSE_CODON_ADAPTIVENESS,
    "CHO_K1": CHO_CODON_ADAPTIVENESS,
    "Saccharomyces_cerevisiae": YEAST_CODON_ADAPTIVENESS,
}

PREFERRED_CODON_TABLES: dict[str, dict[str, str]] = {
    "Homo_sapiens": HUMAN_PREFERRED_CODONS,
    "Escherichia_coli": E_COLI_PREFERRED_CODONS,
    "Mus_musculus": MOUSE_PREFERRED_CODONS,
    "CHO_K1": CHO_PREFERRED_CODONS,
    "Saccharomyces_cerevisiae": YEAST_PREFERRED_CODONS,
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
    SHARP_LI_ADAPTIVENESS_TABLES = tables
    return tables


# Aliases for backward compat
HUMAN: str = "Homo_sapiens"
E_COLI: str = "Escherichia_coli"
MOUSE: str = "Mus_musculus"
CHO: str = "CHO_K1"
YEAST: str = "Saccharomyces_cerevisiae"

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


ECOLI_CAI: dict[str, float] = compute_cai_weights(ECOLI_CODON_USAGE)

HUMAN_CAI: dict[str, float] = compute_cai_weights(HUMAN_CODON_USAGE_SIMPLE)

MOUSE_CAI: dict[str, float] = MOUSE_CODON_ADAPTIVENESS

CHO_CAI: dict[str, float] = CHO_CODON_ADAPTIVENESS

YEAST_CAI: dict[str, float] = YEAST_CODON_ADAPTIVENESS


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


SPECIES: dict[str, SpeciesEntry] = {
    "ecoli": {"cai_weights": ECOLI_CAI, "codon_usage_validation": True},
    "human": {"cai_weights": HUMAN_CAI, "codon_usage_validation": True},
    "mouse": {"cai_weights": MOUSE_CAI, "codon_usage_validation": True},
    "cho": {"cai_weights": CHO_CAI, "codon_usage_validation": True},
    "yeast": {"cai_weights": YEAST_CAI, "codon_usage_validation": True},
}

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
