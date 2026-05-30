"""
BioCompiler Organism Data

Multi-organism codon usage tables loaded from structured data.
Extensible — add a new organism by adding a new module.
"""

from .human import HUMAN_CODON_USAGE, HUMAN_CODON_ADAPTIVENESS, HUMAN_PREFERRED_CODONS
from .e_coli import E_COLI_CODON_USAGE, E_COLI_CODON_ADAPTIVENESS, E_COLI_PREFERRED_CODONS

# Registry: organism name -> codon data
CODON_USAGE_TABLES: dict[str, dict[str, tuple[str, float, float, int]]] = {
    "Homo_sapiens": HUMAN_CODON_USAGE,
    "Escherichia_coli": E_COLI_CODON_USAGE,
}

CODON_ADAPTIVENESS_TABLES: dict[str, dict[str, float]] = {
    "Homo_sapiens": HUMAN_CODON_ADAPTIVENESS,
    "Escherichia_coli": E_COLI_CODON_ADAPTIVENESS,
}

PREFERRED_CODON_TABLES: dict[str, dict[str, str]] = {
    "Homo_sapiens": HUMAN_PREFERRED_CODONS,
    "Escherichia_coli": E_COLI_PREFERRED_CODONS,
}

SUPPORTED_ORGANISMS = list(CODON_USAGE_TABLES.keys())

# Alias for backward compat
HUMAN = "Homo_sapiens"
E_COLI = "Escherichia_coli"
