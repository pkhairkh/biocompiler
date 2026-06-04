"""
BioCompiler Restriction Site Database v9.2.0
===============================================
Common restriction enzyme recognition sequences.

This module provides a curated dictionary of restriction enzyme names mapped
to their recognition sequences (5'→3'), along with a case-insensitive lookup
function. Only non-degenerate (pure ACGT) recognition sequences are included;
enzymes requiring IUPAC ambiguity codes (e.g. SfiI with GGCCNNNNNGGCC) must
be handled via the IUPAC expansion utilities in ``constants.py``.
"""

from __future__ import annotations

import logging

__all__: list[str] = [
    "RESTRICTION_SITES",
    "get_recognition_site",
]

logger = logging.getLogger(__name__)

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
    "AluI":    "AGCT",
    "HaeIII":  "GGCC",
    "MspI":    "CCGG",
    "TaqI":    "TCGA",
    "Sau3AI":  "GATC",
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
