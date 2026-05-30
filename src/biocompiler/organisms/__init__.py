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
"""

from .human import HUMAN_CODON_USAGE, HUMAN_CODON_ADAPTIVENESS, HUMAN_PREFERRED_CODONS
from .e_coli import E_COLI_CODON_USAGE, E_COLI_CODON_ADAPTIVENESS, E_COLI_PREFERRED_CODONS
from .mouse import MOUSE_CODON_USAGE, MOUSE_CODON_ADAPTIVENESS, MOUSE_PREFERRED_CODONS
from .cho import CHO_CODON_USAGE, CHO_CODON_ADAPTIVENESS, CHO_PREFERRED_CODONS
from .yeast import YEAST_CODON_USAGE, YEAST_CODON_ADAPTIVENESS, YEAST_PREFERRED_CODONS

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

SUPPORTED_ORGANISMS = list(CODON_USAGE_TABLES.keys())

# Organism-specific natural GC content targets
# Source: Genome GC content from NCBI Genome Assembly data
ORGANISM_GC_TARGETS: dict[str, float] = {
    "Homo_sapiens": 0.41,
    "Escherichia_coli": 0.52,
    "Mus_musculus": 0.42,
    "CHO_K1": 0.43,
    "Saccharomyces_cerevisiae": 0.38,
}

# Aliases for backward compat
HUMAN = "Homo_sapiens"
E_COLI = "Escherichia_coli"
MOUSE = "Mus_musculus"
CHO = "CHO_K1"
YEAST = "Saccharomyces_cerevisiae"
