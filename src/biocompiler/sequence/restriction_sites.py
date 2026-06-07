"""
BioCompiler Restriction Site Database v9.2.1
===============================================
Common restriction enzyme recognition sequences.

This module provides a curated dictionary of restriction enzyme names mapped
to their recognition sequences (5'→3'), along with a case-insensitive lookup
function. Both non-degenerate (pure ACGT) and degenerate (IUPAC ambiguity)
recognition sequences are included.

Key features:
- ``MIN_SITE_LENGTH``: only sites >= 6 bp are considered for restriction site
  elimination (4-cutters are too frequent to eliminate practically).
- ``expand_iupac_site()``: expands a degenerate recognition sequence into all
  concrete ACGT sequences, enabling exact-match scanning (e.g. Aho-Corasick).
- ``get_eliminable_sites()``: returns only sites >= ``MIN_SITE_LENGTH`` that
  are candidates for elimination during gene optimization.
"""

from __future__ import annotations

import logging

__all__: list[str] = [
    "RESTRICTION_SITES",
    "MIN_SITE_LENGTH",
    "get_recognition_site",
    "expand_iupac_site",
    "get_eliminable_sites",
]

logger = logging.getLogger(__name__)

# Minimum recognition-site length for elimination during gene optimization.
# Sites shorter than 6 bp (4-cutters like AluI, HaeIII) occur too frequently
# in typical coding sequences to be practical elimination targets.
MIN_SITE_LENGTH: int = 6

RESTRICTION_SITES: dict[str, str] = {
    # ── 6-cutters (most common cloning enzymes) ──────────────────────────
    "EcoRI":   "GAATTC",
    "BamHI":   "GGATCC",
    "HindIII": "AAGCTT",
    "XhoI":    "CTCGAG",
    "XbaI":    "TCTAGA",
    "SalI":    "GTCGAC",
    "PstI":    "CTGCAG",
    "SphI":    "GCATGC",
    "KpnI":    "GGTACC",
    "SacI":    "GAGCTC",
    "NcoI":    "CCATGG",
    "NdeI":    "CATATG",
    "BglII":   "AGATCT",
    "ClaI":    "ATCGAT",
    "EcoRV":   "GATATC",
    "SmaI":    "CCCGGG",
    "SpeI":    "ACTAGT",
    "NheI":    "GCTAGC",
    "ApaI":    "GGGCCC",
    "MluI":    "ACGCGT",
    "MfeI":    "CAATTG",
    "AgeI":    "ACCGGT",
    "BsiWI":   "CGTACG",
    "BsrGI":   "TGTACA",
    "AflII":   "CTTAAG",
    "ScaI":    "AGTACT",
    "NsiI":    "ATGCAT",
    "DraI":    "TTTAAA",
    "HpaI":    "GTTAAC",
    "BclI":    "TGATCA",
    "SspI":    "AATATT",
    "StuI":    "AGGCCT",
    "NarI":    "GGCGCC",
    "AvrII":   "CCTAGG",
    "AatII":   "GACGTC",
    "BstBI":   "TTCGAA",
    "PmlI":    "CACGTG",
    "NgoMIV":  "GCCGGC",
    "BspEI":   "TCCGGA",
    # ── 8-cutters (rare-cutter enzymes) ──────────────────────────────────
    "NotI":    "GCGGCCGC",
    "SbfI":    "CCTGCAGG",
    "AscI":    "GGCGCGCC",
    "PmeI":    "GTTTAAAC",
    "FseI":    "GGCCGGCC",
    "PacI":    "TTAATTAA",
    # ── 4-cutters (frequent cutters) ─────────────────────────────────────
    # Included for reference/completeness, but excluded from elimination
    # by get_eliminable_sites() since they are too short to eliminate.
    "AluI":    "AGCT",
    "HaeIII":  "GGCC",
    "MspI":    "CCGG",
    "TaqI":    "TCGA",
    "Sau3AI":  "GATC",
    # ── Degenerate (IUPAC) recognition sites ──────────────────────────────
    # These contain ambiguity codes and require expand_iupac_site() or
    # IUPAC-aware matching for correct detection.
    "SfiI":    "GGCCNNNNNGGCC",   # 13 bp with 5 N spacers
    "HincII":  "GTYRAC",          # 6 bp, Y=C/T, R=A/G
}

# Lowercase-keyed lookup for case-insensitive enzyme name resolution.
_LOWER_LOOKUP: dict[str, str] = {
    name.lower(): seq for name, seq in RESTRICTION_SITES.items()
}


def get_recognition_site(enzyme: str) -> str | None:
    """Get recognition sequence for a restriction enzyme (case-insensitive).

    Args:
        enzyme: Name of the restriction enzyme (e.g. ``"EcoRI"``, ``"ecori"``).

    Returns:
        Recognition sequence string, or ``None`` if the enzyme is unknown.
    """
    result: str | None = _LOWER_LOOKUP.get(enzyme.lower())
    if result is None:
        logger.debug("Unknown restriction enzyme requested: %r", enzyme)
    return result


# ── IUPAC expansion ────────────────────────────────────────────────────────────

# Canonical IUPAC ambiguity code → set of concrete bases
# Covers all 15 degenerate codes plus the 4 deterministic bases.
_IUPAC_MAP: dict[str, str] = {
    "A": "A",  "C": "C",  "G": "G",  "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGT",
}


def expand_iupac_site(site: str) -> list[str]:
    """Expand an IUPAC recognition sequence into all concrete ACGT sequences.

    For example, ``"GTYRAC"`` expands to::

        ["GTCGAC", "GTCAGC", "GTTGAC", "GTTAGC"]

    which are the four concrete sequences that HincII recognises.

    This is essential for feeding degenerate recognition sites into the
    Aho-Corasick automaton, which only handles exact string matching.

    Args:
        site: Recognition sequence that may contain IUPAC ambiguity codes
            (uppercase or lowercase).

    Returns:
        Sorted list of all concrete ACGT sequences matching the pattern.
        Returns ``[site.upper()]`` as a singleton if *site* is already
        pure ACGT.
    """
    site_upper = site.upper()
    # Fast path: pure ACGT
    if all(b in "ACGT" for b in site_upper):
        return [site_upper]

    expansions: list[str] = [""]
    for base in site_upper:
        options = _IUPAC_MAP.get(base)
        if options is None:
            raise ValueError(
                f"Unknown IUPAC code '{base}' in site '{site}'. "
                f"Valid codes: {sorted(_IUPAC_MAP)}"
            )
        expansions = [prefix + opt for prefix in expansions for opt in options]
    expansions.sort()
    return expansions


def get_eliminable_sites(min_length: int = MIN_SITE_LENGTH) -> dict[str, str]:
    """Return restriction sites eligible for elimination during optimization.

    Filters to sites whose recognition sequence is >= *min_length* base pairs.
    Short sites (4-cutters) are excluded because they occur too frequently
    in coding sequences to be practical elimination targets.

    IUPAC sites (e.g. ``"GTYRAC"``) are included as long as their total
    length meets the threshold.  Use :func:`expand_iupac_site` to expand
    them into concrete sequences for exact-match scanning.

    Args:
        min_length: Minimum recognition-site length (default 6).

    Returns:
        Dictionary of enzyme_name -> recognition_sequence for eliminable sites.
    """
    return {
        name: seq
        for name, seq in RESTRICTION_SITES.items()
        if len(seq) >= min_length
    }
