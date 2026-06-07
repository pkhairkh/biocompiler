"""
tRNA Adaptation Index (tAI)
============================

Calculate the tRNA Adaptation Index for DNA coding sequences.

tAI is a measure of codon usage bias that accounts for tRNA gene copy
numbers (a proxy for tRNA abundance) and wobble base pairing rules.
It was introduced by dos Reis et al. (2004) and is considered a more
biologically meaningful metric than CAI because:

1. **tRNA gene copy numbers**: tAI uses the actual number of tRNA genes
   in the genome (from GtRNAdb) as a proxy for tRNA abundance, rather
   than codon frequencies which can be influenced by mutational bias.
2. **Wobble base pairing**: tAI accounts for the fact that not all
   anticodon-codon pairings are equal — some tRNAs can read multiple
   codons through wobble base pairing, with different efficiencies.
3. **Expression prediction**: tAI better predicts expression levels
   for heterologous genes because it captures the translational
   capacity of the host organism's tRNA pool.

The tAI for a gene is defined as the geometric mean of the relative
adaptiveness values of its codons:

    tAI = (∏_{i=1}^{L} w_i)^{1/L}

where w_i is the relative adaptiveness of codon i, defined as:

    w_i = ∑_{j} T_{aa,j} × s(codon, anticodon_j)

where T_{aa,j} is the gene copy number of the j-th tRNA for amino acid
aa, and s(codon, anticodon_j) is the wobble efficiency factor for the
pairing between codon and anticodon j.

Data sources
------------
tRNA gene copy numbers are from the Genomic tRNA Database (GtRNAdb):
    Chan, P.P. & Lowe, T.M. (2016) Nucleic Acids Research 44:D184-D189

References
----------
dos Reis, M., Savva, R. & Wernisch, L. (2004). Solving the riddle of
codon usage preferences: a test for translational selection.
*Nucleic Acids Research*, 32(17), 5036-5044.
doi:10.1093/nar/gkh834
"""

from __future__ import annotations

import math
import warnings
from typing import Optional

from ..constants import CODON_TABLE
from ..organisms import resolve_organism
from ..organisms.tai_data import (
    TRNA_GENE_COPIES,
    WOBBLE_EFFICIENCY,
    WOBBLE_RULES,
    SUPPORTED_ORGANISMS_TAI,
    compute_tai_weights,
)

__all__ = [
    "compute_tai",
    "calculate_tai",
    "TRNA_GENE_COPIES",
    "WOBBLE_RULES",
    "WOBBLE_EFFICIENCY",
    "SUPPORTED_ORGANISMS_TAI",
    "compute_codon_weights",
    "compute_tai_and_cai",
    "optimize_for_tai",
]


# Canonical organism name mapping for tRNA data
_TAI_ORGANISM_ALIASES: dict[str, str] = {
    # E. coli
    "Escherichia_coli": "e_coli",
    "e_coli": "e_coli",
    "ecoli": "e_coli",
    "E. coli": "e_coli",
    "E_coli": "e_coli",
    # Human
    "Homo_sapiens": "human",
    "human": "human",
    "H. sapiens": "human",
    "h_sapiens": "human",
    "H_sapiens": "human",
    # Yeast
    "Saccharomyces_cerevisiae": "yeast",
    "yeast": "yeast",
    "S. cerevisiae": "yeast",
    "s_cerevisiae": "yeast",
    "S_cerevisiae": "yeast",
    # Mouse
    "Mus_musculus": "mouse",
    "mouse": "mouse",
    "M. musculus": "mouse",
    "m_musculus": "mouse",
    "M_musculus": "mouse",
    # CHO-K1
    "CHO_K1": "cho",
    "cho": "cho",
    "CHO": "cho",
    "Cricetulus_griseus": "cho",
}


# Epsilon floor for zero-adaptiveness codons
_TAI_EPSILON: float = 1e-10


def _resolve_tai_organism(organism: str) -> str:
    """Resolve an organism name to a tAI data key.

    Parameters
    ----------
    organism : str
        Any organism name or alias.

    Returns
    -------
    str
        Key in TRNA_GENE_COPIES.

    Raises
    ------
    ValueError
        If no tRNA data is available for the organism.
    """
    # Try direct alias resolution
    key = _TAI_ORGANISM_ALIASES.get(organism)
    if key is not None and key in TRNA_GENE_COPIES:
        return key

    # Try resolve_organism then alias
    resolved = resolve_organism(organism, strict=False)
    key = _TAI_ORGANISM_ALIASES.get(resolved)
    if key is not None and key in TRNA_GENE_COPIES:
        return key

    # Try lowercase
    key = _TAI_ORGANISM_ALIASES.get(organism.lower())
    if key is not None and key in TRNA_GENE_COPIES:
        return key

    # Direct lookup
    if organism in TRNA_GENE_COPIES:
        return organism

    raise ValueError(
        f"No tRNA gene copy data available for organism '{organism}'. "
        f"Available organisms: {list(TRNA_GENE_COPIES.keys())}"
    )


def _dna_to_rna(dna: str) -> str:
    """Convert DNA sequence to RNA (T → U)."""
    return dna.upper().replace("T", "U")


def _compute_codon_weight(
    rna_codon: str,
    trna_copies: dict[str, int],
) -> float:
    """Compute the raw weight for a single codon.

    The weight is the sum over all anticodons that can read this codon,
    of (tRNA_gene_copies × wobble_efficiency):

        W(codon) = Σ_j T(j) × s(codon, anticodon_j)

    Parameters
    ----------
    rna_codon : str
        Codon in RNA format (U instead of T).
    trna_copies : dict
        tRNA gene copy numbers {anticodon: count}.

    Returns
    -------
    float
        Raw weight for this codon.
    """
    wobble_list = WOBBLE_RULES.get(rna_codon, [])
    total_weight = 0.0

    for anticodon, efficiency in wobble_list:
        copies = trna_copies.get(anticodon, 0)
        total_weight += copies * efficiency

    return total_weight


def compute_codon_weights(
    organism: str,
) -> dict[str, float]:
    """Compute relative adaptiveness weights for all codons for an organism.

    For each codon, the relative adaptiveness is:

        w(codon) = W(codon) / W_max(amino_acid)

    where W(codon) is the raw weight and W_max is the maximum raw weight
    among all synonymous codons for the same amino acid.

    Parameters
    ----------
    organism : str
        Organism name.  Accepts canonical binomials, short keys,
        or any alias recognised by ``_resolve_tai_organism()``.

    Returns
    -------
    dict[str, float]
        Mapping of RNA codon → relative adaptiveness value in [0, 1].
    """
    tai_key = _resolve_tai_organism(organism)
    return compute_tai_weights(tai_key)


def compute_tai(
    sequence: str,
    organism: str = "Escherichia_coli",
    species: str | None = None,
    *,
    skip_stop: bool = True,
    skip_met: bool = True,
) -> float:
    """Compute the tRNA Adaptation Index (tAI) for a DNA coding sequence.

    tAI is the geometric mean of relative adaptiveness values for all
    codons in the sequence (excluding Met and stop codons by default).

    The tAI differs from CAI in that it uses tRNA gene copy numbers
    (a proxy for tRNA abundance) and wobble pairing rules, rather than
    codon frequency data from highly expressed genes.

    Organism Specification:

        The target organism can be specified using **either** the
        ``organism`` parameter **or** the ``species`` parameter.  Both
        accept the same set of names — short aliases, abbreviated
        binomials, display names, or full canonical names — and both
        map to the same internal representation via
        :func:`~biocompiler.organisms.resolve_organism`.

        If both ``species`` and ``organism`` are provided, ``species``
        takes precedence and a :class:`DeprecationWarning` is emitted.

    Args:
        sequence: DNA coding sequence (length must be a multiple of 3).
        organism: Organism name.  Accepts canonical binomials
            (e.g., ``'Homo_sapiens'``, ``'Escherichia_coli'``),
            short keys (``'ecoli'``, ``'human'``), abbreviated
            binomials (``'E_coli'``, ``'h_sapiens'``), or display
            names (``'E. coli'``).  All forms are resolved via
            :func:`~biocompiler.organisms.resolve_organism`.
        species: Alias for ``organism``.  Accepts the same values.
            If provided **together with** ``organism``, ``species``
            takes precedence and a deprecation warning is emitted.
            Prefer using ``organism`` in new code; ``species`` is
            retained for backward compatibility.
        skip_stop: Whether to exclude stop codons from the calculation
            (default True).
        skip_met: Whether to exclude the Met (ATG) codon from the
            calculation (default True, following the CAI convention).

    Returns:
        tAI value in [0, 1]. Returns 0.0 for empty or invalid sequences.

    Raises:
        ValueError: If the DNA length is not a multiple of 3.
        ValueError: If no tRNA data is available for the organism.

    References
    ----------
    dos Reis, M., Savva, R. & Wernisch, L. (2004). Solving the riddle of
    codon usage preferences: a test for translational selection.
    *Nucleic Acids Research*, 32(17), 5036-5044.
    """
    # ── Organism resolution (same pattern as compute_cai) ────────────
    if species is not None:
        resolved = resolve_organism(species, strict=False)
        if organism != "Escherichia_coli":
            resolved_explicit = resolve_organism(organism, strict=False)
            if resolved != resolved_explicit:
                warnings.warn(
                    f"Both 'species={species!r}' and 'organism={organism!r}' "
                    f"were provided but resolve to different organisms "
                    f"({resolved!r} vs {resolved_explicit!r}). "
                    f"Using 'species' ({resolved!r}). "
                    f"Prefer using only 'organism' in new code.",
                    DeprecationWarning,
                    stacklevel=2,
                )
            else:
                warnings.warn(
                    f"Both 'species' and 'organism' were provided. "
                    f"Prefer using only 'organism' in new code; "
                    f"'species' is retained for backward compatibility.",
                    DeprecationWarning,
                    stacklevel=2,
                )
        else:
            warnings.warn(
                f"The 'species' parameter is deprecated in favor of 'organism'. "
                f"Use organism='{resolved}' instead of "
                f"species='{species}'. Both accept the same aliases.",
                DeprecationWarning,
                stacklevel=2,
            )
        organism = resolved
    else:
        organism = resolve_organism(organism, strict=False)

    return calculate_tai(
        sequence,
        organism,
        skip_stop=skip_stop,
        skip_met=skip_met,
    )


def calculate_tai(
    dna: str,
    organism: str,
    *,
    skip_stop: bool = True,
    skip_met: bool = True,
) -> float:
    """Calculate the tRNA Adaptation Index (tAI) for a DNA sequence.

    tAI is the geometric mean of relative adaptiveness values for all
    codons in the sequence (excluding Met and stop codons by default).

    The tAI differs from CAI in that it uses tRNA gene copy numbers
    (a proxy for tRNA abundance) and wobble pairing rules, rather than
    codon frequency data from highly expressed genes.

    Parameters
    ----------
    dna : str
        DNA coding sequence (length must be a multiple of 3).
    organism : str
        Target organism name. Must have tRNA data in TRNA_GENE_COPIES.
        Currently supports: ``"e_coli"``, ``"human"``, ``"yeast"``,
        ``"mouse"``, ``"cho"`` (and their canonical/alias names).
    skip_stop : bool
        Whether to exclude stop codons from the calculation (default True).
    skip_met : bool
        Whether to exclude the Met (ATG) codon from the calculation
        (default True, following the CAI convention).

    Returns
    -------
    float
        tAI value in [0, 1].  Returns 0.0 for empty or invalid sequences.

    Raises
    ------
    ValueError
        If the DNA length is not a multiple of 3.
    ValueError
        If no tRNA data is available for the organism.

    References
    ----------
    dos Reis, M., Savva, R. & Wernisch, L. (2004). Solving the riddle of
    codon usage preferences: a test for translational selection.
    *Nucleic Acids Research*, 32(17), 5036-5044.
    """
    if not dna or len(dna) < 3:
        return 0.0

    dna = dna.upper().strip()
    if len(dna) % 3 != 0:
        raise ValueError(
            f"DNA sequence length ({len(dna)}) is not a multiple of 3"
        )

    # Resolve organism and get tRNA data
    tai_key = _resolve_tai_organism(organism)

    # Compute codon weights
    weights = compute_tai_weights(tai_key)

    # Compute tAI as geometric mean
    log_sum: float = 0.0
    count: int = 0

    for i in range(0, len(dna), 3):
        dna_codon = dna[i:i + 3]
        rna_codon = _dna_to_rna(dna_codon)

        aa = CODON_TABLE.get(dna_codon)
        if aa is None:
            continue
        if skip_stop and aa == "*":
            continue
        if skip_met and aa == "M":
            continue

        # Get relative adaptiveness
        w = weights.get(rna_codon, 0.0)
        if w <= 0:
            w = _TAI_EPSILON
        log_sum += math.log(w)
        count += 1

    if count == 0:
        return 0.0

    tai = math.exp(log_sum / count)
    return round(tai, 4)


def compute_tai_and_cai(dna: str, organism: str = "Escherichia_coli",
                        **kwargs) -> dict:
    """Compute both tAI and CAI for a DNA sequence.

    Convenience function that returns both metrics in a single call,
    along with a correlation indicator.

    Args:
        dna: DNA sequence string.
        organism: Target organism name or alias.

    Returns:
        Dictionary with keys:
            - 'tai': tRNA Adaptation Index value (0.0-1.0)
            - 'cai': Codon Adaptation Index value (0.0-1.0)
            - 'correlation': +1.0 if both metrics agree on direction,
              -1.0 if they disagree.
    """
    from ..translation import compute_cai
    tai_val = compute_tai(dna, organism=organism, **kwargs)
    cai_val = compute_cai(dna, organism=organism)
    # Simple correlation: both high or both low → +1, otherwise -1
    mid_tai = 0.5
    mid_cai = 0.5
    correlation = 1.0 if (tai_val >= mid_tai) == (cai_val >= mid_cai) else -1.0
    return {"tai": tai_val, "cai": cai_val, "correlation": correlation}


def optimize_for_tai(protein: str, organism: str = "Escherichia_coli") -> str:
    """Optimize a protein sequence for maximum tAI.

    Selects the codon with the highest tAI weight for each amino acid.

    Args:
        protein: Amino acid sequence string.
        organism: Target organism name or alias.

    Returns:
        Optimized DNA sequence string.

    Raises:
        ValueError: If protein contains invalid amino acid codes.
    """
    from ..constants import CODON_TABLE
    weights = compute_codon_weights(organism)
    # Build reverse map: amino acid → list of (codon_dna, weight)
    aa_to_codons: dict[str, list[tuple[str, float]]] = {}
    for codon_dna, aas in CODON_TABLE.items():
        for aa in aas:
            if aa == "*" or aa == "M":
                continue
            codon_rna = codon_dna.replace("T", "U")
            w = weights.get(codon_rna, 0.0)
            aa_to_codons.setdefault(aa, []).append((codon_dna, w))

    result_codons = []
    for aa in protein:
        if aa == "M":
            result_codons.append("ATG")
            continue
        if aa == "*":
            result_codons.append("TAA")
            continue
        candidates = aa_to_codons.get(aa)
        if candidates is None:
            raise ValueError(f"Invalid amino acid: {aa}")
        best_codon = max(candidates, key=lambda x: x[1])[0]
        result_codons.append(best_codon)

    return "".join(result_codons)
