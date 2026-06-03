"""
BioCompiler Restriction Site Database v7.0.0
===============================================
Common restriction enzyme recognition sequences.
"""

from typing import Dict, Optional

RESTRICTION_SITES: Dict[str, str] = {
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
    "NotI":    "GCGGCCGC",
    "BglII":   "AGATCT",
    "ClaI":    "ATCGAT",
    "EcoRV":   "GATATC",
    "SmaI":    "CCCGGG",
    "SpeI":    "ACTAGT",
    "NheI":    "GCTAGC",
    "ApaI":    "GGGCCC",
    "AluI":    "AGCT",
    "HaeIII":  "GGCC",
    "MspI":    "CCGG",
    "TaqI":    "TCGA",
    "Sau3AI":  "GATC",
    "SbfI":    "CCTGCAGG",
    "AscI":    "GGCGCGCC",
    "PmeI":    "GTTTAAAC",
    "FseI":    "GGCCGGCC",
    "PacI":    "TTAATTAA",
}


def get_recognition_site(enzyme: str) -> Optional[str]:
    """Get recognition sequence for a restriction enzyme (case-insensitive)."""
    return RESTRICTION_SITES.get(enzyme) or RESTRICTION_SITES.get(enzyme.lower())
