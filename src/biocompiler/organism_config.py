"""
BioCompiler Organism Configuration

Organism-specific intelligence for gene optimization.  Each entry
captures the parameters that guide codon choice, GC targeting,
homopolymer avoidance, and mRNA stability modelling so that the
optimizer can make biologically informed decisions.

Configurations are keyed by a short, stable identifier (e.g.
``"E_coli_K12"``) that is separate from the display / taxonomy name
stored inside each :class:`OrganismConfig` instance.

Public API
----------
- :data:`ORGANISM_CONFIGS` — built-in registry
- :func:`get_organism_config` — lookup with safe fallback
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

__all__: list[str] = [
    "OrganismConfig",
    "ORGANISM_CONFIGS",
    "get_organism_config",
    "is_eukaryotic_organism",
    "auto_detect_organism_domain",
]

logger = logging.getLogger(__name__)

# ─── Data model ────────────────────────────────────────────────────


@dataclass(frozen=True)
class OrganismConfig:
    """Organism-specific parameters for gene optimization.

    Attributes:
        name: Human-readable organism name (e.g. ``"Escherichia coli K-12"``).
        gc_target_lo: Lower bound of the target GC-content window (fraction).
        gc_target_hi: Upper bound of the target GC-content window (fraction).
        codon_usage_validated: Whether the codon-usage table has been
            validated against an external reference (Kazusa, CoCoPUTs, …).
        rbs_calculator_available: Whether a Ribosome Binding Site
            calculator (e.g. Salis Lab RBS Calculator) is available for
            this organism.
        preferred_codons: Mapping from single-letter amino-acid code to
            the preferred codon for that residue (e.g. ``{"L": "CTG"}``).
        avoided_motifs: DNA sequence motifs that must be eliminated from
            the designed construct (e.g. restriction sites, cryptic
            splice signals, instability motifs).
        max_homopolymer_run: Maximum allowed length of a single-base
            repeat (e.g. ``4`` rejects ``"AAAA"``).
        mrna_degradation_model: Fidelity of the mRNA degradation model
            used during optimisation.  One of ``"none"``, ``"simple"``,
            ``"detailed"``.
    """

    name: str
    gc_target_lo: float
    gc_target_hi: float
    codon_usage_validated: bool
    rbs_calculator_available: bool
    preferred_codons: dict[str, str] = field(default_factory=dict)
    avoided_motifs: list[str] = field(default_factory=list)
    max_homopolymer_run: int = 6
    mrna_degradation_model: str = "none"
    domain: str = "eukaryote"

    # ── Derived helpers ─────────────────────────────────────────────

    @property
    def is_eukaryote(self) -> bool:
        """Return ``True`` when the organism belongs to Domain Eukarya."""
        return self.domain == "eukaryote"


# ─── Preferred-codon tables (AA → codon) ───────────────────────────
# Sourced from the organisms sub-package at import time so that
# preferred_codons stay in sync with the canonical codon-usage data.

def _build_preferred_codons(module_name: str) -> dict[str, str]:
    """Lazily import *module_name* from ``biocompiler.organisms`` and
    return its ``*_PREFERRED_CODONS`` mapping.

    Falls back to an empty dict if the import fails, which keeps the
    config usable even when optional data packages are missing.
    """
    try:
        import importlib
        mod = importlib.import_module(f".organisms.{module_name}", __package__)
        # Each organism module exports <PREFIX>_PREFERRED_CODONS
        for attr in dir(mod):
            if attr.endswith("_PREFERRED_CODONS"):
                return getattr(mod, attr)
    except ImportError:
        logger.debug("Could not import preferred codons from organisms.%s", module_name)
    return {}


# ─── Motif lists ───────────────────────────────────────────────────

_E_COLI_AVOIDED_MOTIFS: list[str] = [
    # Common E. coli restriction sites
    "GAATTC",   # EcoRI
    "GGATCC",   # BamHI
    "AAGCTT",   # HindIII
    "TCTAGA",   # XbaI
    # mRNA instability motifs
    "ATTTA",
]

_HUMAN_AVOIDED_MOTIFS: list[str] = [
    # Cryptic splice donor / acceptor signals
    "GTAGGT",
    "NCAGGT",  # requires fuzzy matching; listed as-is for reference
    # mRNA instability motifs
    "ATTTA",
    "TTATTTAT",  # AU-rich element variant
]

_MAMMALIAN_AVOIDED_MOTIFS: list[str] = [
    # Common cryptic splice signals
    "GTAGGT",
    # mRNA instability
    "ATTTA",
]

_YEAST_AVOIDED_MOTIFS: list[str] = [
    # S. cerevisiae does not have strong restriction-site concerns
    # for standard cloning, but AU-rich instability elements matter.
    "ATTTA",
]


# ─── Built-in registry ────────────────────────────────────────────

ORGANISM_CONFIGS: dict[str, OrganismConfig] = {
    # ── E. coli K-12 ──────────────────────────────────────────────
    "E_coli_K12": OrganismConfig(
        name="Escherichia coli K-12",
        gc_target_lo=0.45,
        gc_target_hi=0.55,
        codon_usage_validated=True,
        rbs_calculator_available=False,
        domain="prokaryote",
        preferred_codons=_build_preferred_codons("e_coli"),
        avoided_motifs=_E_COLI_AVOIDED_MOTIFS,
        max_homopolymer_run=5,
        mrna_degradation_model="simple",
    ),

    # ── E. coli BL21(DE3) ─────────────────────────────────────────
    "E_coli_BL21": OrganismConfig(
        name="Escherichia coli BL21(DE3)",
        gc_target_lo=0.45,
        gc_target_hi=0.55,
        codon_usage_validated=True,
        rbs_calculator_available=True,
        domain="prokaryote",
        preferred_codons=_build_preferred_codons("e_coli"),
        avoided_motifs=_E_COLI_AVOIDED_MOTIFS,
        max_homopolymer_run=5,
        mrna_degradation_model="simple",
    ),

    # ── Homo sapiens ──────────────────────────────────────────────
    "Homo_sapiens": OrganismConfig(
        name="Homo sapiens",
        gc_target_lo=0.40,
        gc_target_hi=0.60,
        codon_usage_validated=True,
        rbs_calculator_available=False,
        domain="eukaryote",
        preferred_codons=_build_preferred_codons("human"),
        avoided_motifs=_HUMAN_AVOIDED_MOTIFS,
        max_homopolymer_run=5,
        mrna_degradation_model="detailed",
    ),

    # ── Saccharomyces cerevisiae ──────────────────────────────────
    "Saccharomyces_cerevisiae": OrganismConfig(
        name="Saccharomyces cerevisiae",
        gc_target_lo=0.35,
        gc_target_hi=0.45,
        codon_usage_validated=True,
        rbs_calculator_available=False,
        domain="eukaryote",
        preferred_codons=_build_preferred_codons("yeast"),
        avoided_motifs=_YEAST_AVOIDED_MOTIFS,
        max_homopolymer_run=5,
        mrna_degradation_model="simple",
    ),

    # ── Mus musculus ──────────────────────────────────────────────
    "Mus_musculus": OrganismConfig(
        name="Mus musculus",
        gc_target_lo=0.40,
        gc_target_hi=0.55,
        codon_usage_validated=False,
        rbs_calculator_available=False,
        domain="eukaryote",
        preferred_codons=_build_preferred_codons("mouse"),
        avoided_motifs=_MAMMALIAN_AVOIDED_MOTIFS,
        max_homopolymer_run=5,
        mrna_degradation_model="simple",
    ),

    # ── CHO-K1 ────────────────────────────────────────────────────
    "CHO_K1": OrganismConfig(
        name="CHO-K1 (Cricetulus griseus)",
        gc_target_lo=0.40,
        gc_target_hi=0.60,
        codon_usage_validated=False,
        rbs_calculator_available=False,
        domain="eukaryote",
        preferred_codons=_build_preferred_codons("cho"),
        avoided_motifs=_MAMMALIAN_AVOIDED_MOTIFS,
        max_homopolymer_run=5,
        mrna_degradation_model="simple",
    ),
}


# ─── Lookup helper ─────────────────────────────────────────────────

# Fallback config returned when the caller requests an unknown organism.
# Uses deliberately permissive defaults so that the optimizer can still
# run (albeit without organism-specific intelligence).
_FALLBACK_CONFIG = OrganismConfig(
    name="Unknown organism (fallback)",
    gc_target_lo=0.30,
    gc_target_hi=0.70,
    codon_usage_validated=False,
    rbs_calculator_available=False,
    preferred_codons={},
    avoided_motifs=[],
    max_homopolymer_run=6,
    mrna_degradation_model="none",
)


def get_organism_config(organism: str) -> OrganismConfig:
    """Look up an :class:`OrganismConfig` by key.

    Args:
        organism: Configuration key (e.g. ``"E_coli_K12"``).
            Also accepts legacy aliases such as ``"Escherichia_coli"``
            or ``"ecoli"`` for backward compatibility.

    Returns:
        The matching :class:`OrganismConfig`, or the fallback config
        with a warning if *organism* is not found.

    Examples::

        cfg = get_organism_config("E_coli_K12")
        print(cfg.gc_target_lo)   # 0.45

        # Legacy alias works too
        cfg = get_organism_config("Escherichia_coli")
    """
    # Direct hit
    if organism in ORGANISM_CONFIGS:
        return ORGANISM_CONFIGS[organism]

    # ── Legacy alias resolution ────────────────────────────────
    _ALIASES: dict[str, str] = {
        "Escherichia_coli": "E_coli_K12",
        "ecoli": "E_coli_K12",
        "E_coli": "E_coli_K12",
        "human": "Homo_sapiens",
        "mouse": "Mus_musculus",
        "cho": "CHO_K1",
        "yeast": "Saccharomyces_cerevisiae",
    }
    canonical = _ALIASES.get(organism)
    if canonical and canonical in ORGANISM_CONFIGS:
        logger.info(
            "Resolved organism alias %r -> %r", organism, canonical,
        )
        return ORGANISM_CONFIGS[canonical]

    # ── Fallback ───────────────────────────────────────────────
    available = list(ORGANISM_CONFIGS.keys())
    logger.warning(
        "Unknown organism %r; using fallback config. "
        "Available: %s",
        organism,
        available,
    )
    return _FALLBACK_CONFIG


def is_eukaryotic_organism(name: str) -> bool:
    """Return ``True`` if the organism identified by *name* is eukaryotic.

    Uses :func:`get_organism_config` internally so that legacy aliases
    are resolved and the safe fallback (eukaryote) is returned for
    unknown organisms.

    Args:
        name: Organism key (e.g. ``"E_coli_K12"``, ``"human"``) or
            any legacy alias accepted by :func:`get_organism_config`.

    Returns:
        ``True`` when the resolved config has ``domain == "eukaryote"``.

    Examples::

        >>> is_eukaryotic_organism("E_coli_K12")
        False
        >>> is_eukaryotic_organism("Homo_sapiens")
        True
        >>> is_eukaryotic_organism("unknown_organism")
        True   # fallback defaults to eukaryote
    """
    return get_organism_config(name).is_eukaryote


def auto_detect_organism_domain(name: str) -> str:
    """Return the domain of life for the given organism name.

    Uses :func:`get_organism_config` internally so that legacy aliases
    are resolved.  Falls back to ``"eukaryote"`` for unknown organisms.

    Args:
        name: Organism key or legacy alias.

    Returns:
        One of ``"eukaryote"``, ``"prokaryote"``, or ``"archaea"``.

    Examples::

        >>> auto_detect_organism_domain("E_coli_K12")
        'prokaryote'
        >>> auto_detect_organism_domain("Homo_sapiens")
        'eukaryote'
    """
    return get_organism_config(name).domain
