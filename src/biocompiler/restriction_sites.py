"""
BioCompiler Restriction Site Database v7.0.0
===============================================
Common restriction enzyme recognition sequences.
"""

RESTRICTION_SITES: dict[str, str] = {
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

# Lowercase-keyed lookup for case-insensitive enzyme name resolution.
_LOWER_LOOKUP: dict[str, str] = {name.lower(): seq for name, seq in RESTRICTION_SITES.items()}


def get_recognition_site(enzyme: str) -> str | None:
    """Get recognition sequence for a restriction enzyme (case-insensitive).

    Args:
        enzyme: Name of the restriction enzyme (e.g. ``"EcoRI"``, ``"ecori"``).

    Returns:
        Recognition sequence string, or ``None`` if the enzyme is unknown.
    """
    return _LOWER_LOOKUP.get(enzyme.lower())
