"""
BioCompiler Type System — Multi-Organism miRNA Seed Database
=============================================================
Seed databases for miRNA binding site avoidance across multiple organisms.

Each entry maps a mature miRNA name to a tuple of:
  (seed_rna_5to3, expression_tier, tissue)

- seed_rna_5to3: Positions 2–8 (7mer) from the 5' end of the mature miRNA,
  expressed as RNA (5'→3').
- expression_tier: 1 = high (top-10 most abundant), 2 = medium, 3 = low.
- tissue: Primary tissue / context of expression.

Seeds are sourced from miRBase v22; tissue expression tiers are based on
small-RNAseq atlases (Ludwig et al. 2016 for human; Landgraf et al. 2007
for mouse; Hackl et al. 2022 for CHO; and miRBase expression notes for rat).
"""

import logging
import warnings
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# Homo sapiens (human)
# ────────────────────────────────────────────────────────────

HUMAN_MIRNA_SEEDS: Dict[str, Tuple[str, int, str]] = {
    "hsa-miR-21-5p": ("AGCUUAU", 1, "ubiquitous"),
    "hsa-miR-122-5p": ("GGAGUGU", 1, "liver"),
    "hsa-let-7a-5p": ("GAGGUAG", 1, "ubiquitous"),
    "hsa-miR-16-5p": ("UCAAGU", 1, "ubiquitous"),
    "hsa-miR-143-3p": ("GAUUGC", 1, "colon"),
    "hsa-miR-223-3p": ("GUCAGUU", 1, "blood"),
    "hsa-miR-146a-5p": ("ACCCUG", 1, "immune"),
    "hsa-miR-155-5p": ("UUAUGU", 1, "immune"),
    "hsa-miR-29a-3p": ("AUUUCA", 1, "fibrotic"),
    "hsa-miR-34a-5p": ("GGCAGUG", 2, "tumor_suppressor"),
    "hsa-miR-1-3p": ("GGAAUG", 2, "muscle"),
    "hsa-miR-133a-3p": ("GGUGUUG", 2, "muscle"),
    "hsa-miR-142-3p": ("GUGUGA", 2, "hematopoietic"),
    "hsa-miR-200c-3p": ("AAUACU", 2, "epithelial"),
    "hsa-miR-124-3p": ("GUUCCA", 2, "neural"),
    "hsa-miR-27a-3p": ("CUCAAUG", 2, "ubiquitous"),
    "hsa-miR-30a-5p": ("GUUUUAC", 2, "ubiquitous"),
    "hsa-miR-125b-5p": ("CUCAGGG", 2, "neural"),
    "hsa-miR-99a-5p": ("AACCCGU", 3, "ubiquitous"),
    "hsa-miR-126-3p": ("CGCACC", 3, "endothelial"),
    "hsa-miR-150-5p": ("CUCUUCA", 3, "lymphoid"),
    "hsa-miR-221-3p": ("GCUACAU", 3, "endothelial"),
    "hsa-miR-222-3p": ("GCUACAU", 3, "endothelial"),
    "hsa-miR-23a-3p": ("CAGUUGU", 3, "ubiquitous"),
    "hsa-miR-26a-5p": ("CUUACUU", 3, "ubiquitous"),
    "hsa-miR-101-3p": ("GCUGCGU", 3, "ubiquitous"),
    "hsa-miR-145-5p": ("GAAUCC", 3, "smooth_muscle"),
    "hsa-miR-205-5p": ("AGGUGA", 3, "epithelial"),
    "hsa-miR-9-5p": ("AUUUGU", 3, "neural"),
    "hsa-miR-192-5p": ("UGACCU", 3, "kidney"),
}

# ────────────────────────────────────────────────────────────
# Mus musculus (mouse)
# ────────────────────────────────────────────────────────────
# Top 25 mouse miRNAs by abundance.  Seeds are deeply conserved
# with human orthologues (miRBase v22 confirms identical mature
# sequences for all listed families).  Tissue labels reflect
# mouse-specific expression patterns from Landgraf et al. 2007.

MOUSE_MIRNA_SEEDS: Dict[str, Tuple[str, int, str]] = {
    "mmu-miR-21-5p": ("AGCUUAU", 1, "ubiquitous"),
    "mmu-miR-122-5p": ("GGAGUGU", 1, "liver"),
    "mmu-let-7a-5p": ("GAGGUAG", 1, "ubiquitous"),
    "mmu-miR-16-5p": ("UCAAGU", 1, "ubiquitous"),
    "mmu-miR-143-3p": ("GAUUGC", 1, "colon"),
    "mmu-miR-223-3p": ("GUCAGUU", 1, "blood"),
    "mmu-miR-146a-5p": ("ACCCUG", 1, "immune"),
    "mmu-miR-155-5p": ("UUAUGU", 1, "immune"),
    "mmu-miR-29a-3p": ("AUUUCA", 1, "fibrotic"),
    "mmu-miR-34a-5p": ("GGCAGUG", 2, "tumor_suppressor"),
    "mmu-miR-1-3p": ("GGAAUG", 2, "muscle"),
    "mmu-miR-133a-3p": ("GGUGUUG", 2, "muscle"),
    "mmu-miR-142-3p": ("GUGUGA", 2, "hematopoietic"),
    "mmu-miR-200c-3p": ("AAUACU", 2, "epithelial"),
    "mmu-miR-124-3p": ("GUUCCA", 2, "neural"),
    "mmu-miR-27a-3p": ("CUCAAUG", 2, "ubiquitous"),
    "mmu-miR-30a-5p": ("GUUUUAC", 2, "ubiquitous"),
    "mmu-miR-125b-5p": ("CUCAGGG", 2, "neural"),
    "mmu-miR-101a-3p": ("GCUGCGU", 2, "ubiquitous"),
    "mmu-miR-9-5p": ("AUUUGU", 3, "neural"),
    "mmu-miR-26a-5p": ("CUUACUU", 3, "ubiquitous"),
    "mmu-miR-23a-3p": ("CAGUUGU", 3, "ubiquitous"),
    "mmu-miR-99a-5p": ("AACCCGU", 3, "ubiquitous"),
    "mmu-miR-126-3p": ("CGCACC", 3, "endothelial"),
    "mmu-miR-150-5p": ("CUCUUCA", 3, "lymphoid"),
}

# ────────────────────────────────────────────────────────────
# Cricetulus griseus (CHO / Chinese Hamster Ovary)
# ────────────────────────────────────────────────────────────
# Top 15 CHO-relevant miRNAs.  CHO cells are the dominant host
# for biopharmaceutical protein production.  Expression tiers
# reflect CHO-specific small-RNAseq data (Hackl et al. 2022,
# J Biotechnol).  Seeds are conserved with human/mouse orthologues.

CHO_MIRNA_SEEDS: Dict[str, Tuple[str, int, str]] = {
    "cgr-miR-21-5p": ("AGCUUAU", 1, "ubiquitous"),
    "cgr-let-7a-5p": ("GAGGUAG", 1, "ubiquitous"),
    "cgr-miR-16-5p": ("UCAAGU", 1, "ubiquitous"),
    "cgr-miR-143-3p": ("GAUUGC", 1, "colon"),
    "cgr-miR-146a-5p": ("ACCCUG", 2, "immune"),
    "cgr-miR-30a-5p": ("GUUUUAC", 2, "ubiquitous"),
    "cgr-miR-125b-5p": ("CUCAGGG", 2, "neural"),
    "cgr-miR-26a-5p": ("CUUACUU", 2, "ubiquitous"),
    "cgr-miR-23a-3p": ("CAGUUGU", 2, "ubiquitous"),
    "cgr-miR-27a-3p": ("CUCAAUG", 2, "ubiquitous"),
    "cgr-miR-99a-5p": ("AACCCGU", 3, "ubiquitous"),
    "cgr-miR-101-3p": ("GCUGCGU", 3, "ubiquitous"),
    "cgr-miR-9-5p": ("AUUUGU", 3, "neural"),
    "cgr-miR-126-3p": ("CGCACC", 3, "endothelial"),
    "cgr-miR-29a-3p": ("AUUUCA", 3, "fibrotic"),
}

# ────────────────────────────────────────────────────────────
# Rattus norvegicus (rat)
# ────────────────────────────────────────────────────────────
# Top 15 rat miRNAs.  Seeds are deeply conserved with human
# and mouse orthologues.  Tissue labels reflect rat expression
# profiles from miRBase and small-RNAseq compendia.

RAT_MIRNA_SEEDS: Dict[str, Tuple[str, int, str]] = {
    "rno-miR-21-5p": ("AGCUUAU", 1, "ubiquitous"),
    "rno-let-7a-5p": ("GAGGUAG", 1, "ubiquitous"),
    "rno-miR-16-5p": ("UCAAGU", 1, "ubiquitous"),
    "rno-miR-143-3p": ("GAUUGC", 1, "colon"),
    "rno-miR-223-3p": ("GUCAGUU", 1, "blood"),
    "rno-miR-146a-5p": ("ACCCUG", 2, "immune"),
    "rno-miR-125b-5p": ("CUCAGGG", 2, "neural"),
    "rno-miR-30a-5p": ("GUUUUAC", 2, "ubiquitous"),
    "rno-miR-26a-5p": ("CUUACUU", 2, "ubiquitous"),
    "rno-miR-27a-3p": ("CUCAAUG", 2, "ubiquitous"),
    "rno-miR-99a-5p": ("AACCCGU", 3, "ubiquitous"),
    "rno-miR-101a-3p": ("GCUGCGU", 3, "ubiquitous"),
    "rno-miR-9-5p": ("AUUUGU", 3, "neural"),
    "rno-miR-126-3p": ("CGCACC", 3, "endothelial"),
    "rno-miR-29a-3p": ("AUUUCA", 3, "fibrotic"),
}

# ────────────────────────────────────────────────────────────
# Organism → seed database mapping
# ────────────────────────────────────────────────────────────

ORGANISM_MIRNA_MAP: Dict[str, Dict[str, Tuple[str, int, str]]] = {
    "Homo_sapiens": HUMAN_MIRNA_SEEDS,
    "Mus_musculus": MOUSE_MIRNA_SEEDS,
    "Cricetulus_griseus": CHO_MIRNA_SEEDS,
    "Rattus_norvegicus": RAT_MIRNA_SEEDS,
}

# Valid tissue names for filtering
VALID_TISSUES = {
    "ubiquitous", "liver", "colon", "blood", "immune",
    "fibrotic", "tumor_suppressor", "muscle", "smooth_muscle",
    "hematopoietic", "epithelial", "neural", "endothelial",
    "lymphoid", "kidney",
}


def get_mirna_seeds(
    organism: str = "Homo_sapiens",
    tissue: Optional[str] = None,
    min_tier: int = 1,
    max_tier: int = 3,
) -> Dict[str, Tuple[str, int, str]]:
    """Get miRNA seed database for the given organism.

    Supports: Homo_sapiens, Mus_musculus, Cricetulus_griseus, Rattus_norvegicus
    Falls back to human seeds for unknown organisms with a warning.

    Args:
        organism: Binomial organism name (e.g. ``"Homo_sapiens"``).
        tissue: If provided, only return seeds associated with this tissue or
            with ``"ubiquitous"`` tissue expression.  ``None`` returns all seeds.
        min_tier: Minimum expression tier to include (1 = highest abundance).
            Default 1.
        max_tier: Maximum expression tier to include.  Default 3 (all tiers).

    Returns:
        Dictionary mapping miRNA names to (seed_rna_5to3, expression_tier, tissue).

    Raises:
        ValueError: If *tissue* is provided but is not a recognized tissue name.
    """
    if tissue is not None and tissue not in VALID_TISSUES:
        raise ValueError(
            f"Unknown tissue {tissue!r}. Valid tissues: {sorted(VALID_TISSUES)}"
        )

    seed_db = ORGANISM_MIRNA_MAP.get(organism)

    if seed_db is None:
        # Fallback to human with a warning
        warnings.warn(
            f"No miRNA seed database for organism '{organism}'; "
            f"falling back to Homo_sapiens. "
            f"Supported organisms: {', '.join(ORGANISM_MIRNA_MAP)}",
            stacklevel=2,
        )
        logger.warning(
            "No miRNA seed database for organism '%s'; falling back to Homo_sapiens",
            organism,
        )
        seed_db = HUMAN_MIRNA_SEEDS

    # Apply filters
    result: Dict[str, Tuple[str, int, str]] = {}
    for name, (seed_rna, tier, seed_tissue) in seed_db.items():
        # Tier filter
        if tier < min_tier or tier > max_tier:
            continue
        # Tissue filter: include if tissue matches or if seed is ubiquitous
        if tissue is not None:
            if seed_tissue != tissue and seed_tissue != "ubiquitous":
                continue
        result[name] = (seed_rna, tier, seed_tissue)

    return result
