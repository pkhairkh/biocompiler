"""
BioCompiler Multi-Gene Construct / Operon Support
===================================================

Support for optimizing multiple genes as a single DNA construct, including:

- **Prokaryotic operons**: Genes with a shared promoter and individual RBS sequences,
  transcribed as a single polycistronic mRNA.
- **Eukaryotic polycistronic constructs**: Genes linked by 2A self-cleaving peptides
  or IRES elements, enabling multi-protein expression from a single transcript.
- **Bidirectional promoters**: Two gene clusters divergently transcribed from a
  shared promoter region.
- **Linker sequences**: Arbitrary DNA spacers between genes (e.g., restriction sites
  for modular cloning, 2A peptide coding sequences, IRES elements).

Architecture
------------

Each gene in the construct is independently codon-optimized for the target organism
using the existing :func:`~biocompiler.optimization.optimize_sequence` pipeline.
The resulting DNA sequences are then assembled with the specified linkers and
regulatory elements (promoters, RBS, terminators) into a single contiguous DNA
construct.  The full construct is exported as a GenBank record with individual
gene features annotated.

API
---

.. autofunction:: optimize_multigene
.. autofunction:: optimize_operon
.. autoclass:: GeneSpec
.. autoclass:: MultiGeneResult
"""

from __future__ import annotations

import logging
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

from .optimization import OptimizationResult, optimize_sequence
from .scanner import gc_content
from .translation import translate
from .organisms import resolve_organism
from .organism_config import is_eukaryotic_organism
from .exceptions import InvalidProteinError, UnsupportedOrganismError

__all__ = [
    "GeneSpec",
    "MultiGeneResult",
    "OperonConfig",
    "optimize_multigene",
    "optimize_operon",
]

logger = logging.getLogger(__name__)


# ─── Well-known linker sequences ────────────────────────────────────

# 2A self-cleaving peptide DNA coding sequences (codon-optimized for mammals).
# These encode short (~18-22 aa) peptides that cause ribosomal "skipping",
# resulting in cleavage between two proteins from a single ORF.
LINKER_2A_P2A: str = "GSGATNFSLLKQAGDVEENPGP"  # P2A peptide (protein)
LINKER_2A_T2A: str = "GSGEGRGSLLTCGDVEENPGP"  # T2A peptide (protein)
LINKER_2A_E2A: str = "GSGNLLSGRDVVEENPGP"     # E2A peptide (protein)
LINKER_2A_F2A: str = "GSGVKQTLNFDLLKLAGDVESNPGP"  # F2A peptide (protein)

# IRES (Internal Ribosome Entry Site) — typically 400-600 bp.
# This is a placeholder identifier; actual IRES sequences vary by source.
# Users should provide their own IRES DNA sequence.
LINKER_IRES_PLACEHOLDER: str = ""  # Must be user-supplied

# Common prokaryotic RBS sequences (Shine-Dalgarno)
RBS_STRONG: str = "AGGAGG"     # Strong consensus Shine-Dalgarno
RBS_MEDIUM: str = "AGGAGG"     # Same consensus; spacing determines strength
RBS_WEAK: str = "AGGA"         # Weaker Shine-Dalgarno

# Default spacer between RBS and start codon (optimal ~5-9 bp for E. coli)
RBS_SPACER_DEFAULT: str = "AAAAA"  # 5 bp spacer


@dataclass
class GeneSpec:
    """Specification for a single gene within a multi-gene construct.

    Attributes:
        protein: Amino acid sequence (single-letter codes, no stop codon).
        organism: Target organism for codon optimization. If empty, the
            construct-level organism is used.
        name: Optional gene name (used for GenBank feature annotation).
        promoter: Optional promoter DNA sequence placed upstream of this gene.
            For operons, only the first gene's promoter is used.
        rbs: Optional ribosome binding site sequence (prokaryotes only).
            Placed immediately upstream of the start codon.
        terminator: Optional terminator DNA sequence placed downstream of
            this gene's stop codon.
    """

    protein: str
    organism: str = ""
    name: str = ""
    promoter: str = ""
    rbs: str = ""
    terminator: str = ""

    def __post_init__(self) -> None:
        """Validate GeneSpec invariants."""
        if not self.protein or not self.protein.strip():
            raise ValueError("GeneSpec.protein must be a non-empty amino acid sequence")
        self.protein = self.protein.upper().strip()
        # Validate amino acid characters
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        invalid = set(self.protein) - valid_aas
        if invalid:
            raise InvalidProteinError(self.protein, invalid)


@dataclass
class OperonConfig:
    """Configuration for prokaryotic operon assembly.

    Attributes:
        promoter: Promoter sequence placed upstream of the first gene.
        rbs_per_gene: RBS sequence placed upstream of each gene in the operon.
        rbs_spacer: DNA spacer between the RBS and the start codon.
        terminator: Terminator sequence placed after the last gene.
        include_restriction_sites: If True, add flanking restriction sites
            for modular cloning (e.g., BioBrick prefix/suffix).
    """

    promoter: str = ""
    rbs_per_gene: str = RBS_STRONG
    rbs_spacer: str = RBS_SPACER_DEFAULT
    terminator: str = ""
    include_restriction_sites: bool = False


@dataclass
class MultiGeneResult:
    """Result of optimizing multiple genes as a single construct.

    Attributes:
        genes: Per-gene optimization results (same length as input ``genes``).
        full_dna: Concatenated DNA sequence including all linkers, promoters,
            RBS, and terminators.
        total_length: Total length of ``full_dna`` in base pairs.
        genbank_export: Full GenBank-format string with individual gene features.
        construct_type: Type of construct (``"operon"``, ``"polycistronic_2A"``,
            ``"polycistronic_IRES"``, ``"bidirectional"``, or ``"custom"``).
        organism: Target organism used for optimization.
        gc_content: GC fraction of the full construct.
    """

    genes: list[OptimizationResult]
    full_dna: str
    total_length: int
    genbank_export: str
    construct_type: str
    organism: str
    gc_content: float = 0.0

    def __post_init__(self) -> None:
        """Validate MultiGeneResult invariants."""
        assert len(self.full_dna) == self.total_length, (
            f"full_dna length ({len(self.full_dna)}) must equal "
            f"total_length ({self.total_length})"
        )
        assert 0.0 <= self.gc_content <= 1.0, (
            f"gc_content must be in [0, 1], got {self.gc_content}"
        )
        assert self.construct_type in (
            "operon", "polycistronic_2A", "polycistronic_IRES",
            "bidirectional", "custom",
        ), f"Invalid construct_type: {self.construct_type!r}"


# ─── Internal helpers ───────────────────────────────────────────────


def _optimize_gene(
    gene_spec: GeneSpec,
    default_organism: str,
    constraints: dict | None = None,
) -> OptimizationResult:
    """Optimize a single gene, resolving organism from spec or default.

    Args:
        gene_spec: The gene specification.
        default_organism: Fallback organism if gene_spec.organism is empty.
        constraints: Optional constraint overrides for optimization.

    Returns:
        OptimizationResult for the individual gene.
    """
    organism = gene_spec.organism or default_organism
    resolved = resolve_organism(organism, strict=False)

    # Build keyword arguments for optimize_sequence
    kwargs: dict[str, Any] = {
        "target_protein": gene_spec.protein,
        "organism": resolved,
        "strict_mode": False,  # Allow partial results for multi-gene constructs
    }
    if constraints:
        # Forward supported constraint parameters
        if "gc_lo" in constraints:
            kwargs["gc_lo"] = constraints["gc_lo"]
        if "gc_hi" in constraints:
            kwargs["gc_hi"] = constraints["gc_hi"]
        if "cai_threshold" in constraints:
            kwargs["cai_threshold"] = constraints["cai_threshold"]
        if "restriction_sites" in constraints:
            kwargs["restriction_sites"] = constraints["restriction_sites"]

    return optimize_sequence(**kwargs)


def _protein_to_2a_dna(protein_sequence: str, organism: str) -> str:
    """Convert a 2A peptide protein sequence to a codon-optimized DNA sequence.

    This is a simplified version — for production use, the full optimizer
    should be called on the short 2A peptide.

    Args:
        protein_sequence: The 2A peptide amino acid sequence.
        organism: Target organism for codon selection.

    Returns:
        DNA sequence encoding the 2A peptide.
    """
    result = optimize_sequence(
        target_protein=protein_sequence,
        organism=organism,
        strict_mode=False,  # Short peptides may not meet GC constraints
    )
    return result.sequence


def _assemble_construct(
    gene_dnas: list[str],
    gene_specs: list[GeneSpec],
    linker: str = "",
    construct_type: str = "custom",
    operon_config: OperonConfig | None = None,
    bidirectional: bool = False,
) -> str:
    """Assemble the full DNA construct from individual gene sequences.

    Assembly logic varies by construct type:

    - **operon**: promoter + (RBS + spacer + gene_DNA)* + terminator
    - **polycistronic_2A**: promoter + gene1 + 2A_DNA + gene2 + 2A_DNA + ... + terminator
    - **polycistronic_IRES**: promoter + gene1 + IRES + gene2 + IRES + ... + terminator
    - **bidirectional**: Two halves separated by a bidirectional promoter
    - **custom**: Genes separated by the specified linker

    Args:
        gene_dnas: Optimized DNA sequences for each gene (with stop codons).
        gene_specs: Original gene specifications (for regulatory elements).
        linker: Linker sequence between genes (used for custom/polycistronic).
        construct_type: Type of construct being assembled.
        operon_config: Configuration for operon assembly.
        bidirectional: Whether to use bidirectional promoter layout.

    Returns:
        The full assembled DNA construct.
    """
    parts: list[str] = []

    if construct_type == "operon":
        config = operon_config or OperonConfig()
        # Promoter before first gene
        if config.promoter:
            parts.append(config.promoter.upper())
        # Each gene: RBS + spacer + gene DNA
        for i, (dna, spec) in enumerate(zip(gene_dnas, gene_specs)):
            # Use gene-specific RBS if provided, otherwise use operon default
            rbs = spec.rbs.upper() if spec.rbs else config.rbs_per_gene.upper()
            if rbs:
                parts.append(rbs)
            if config.rbs_spacer:
                parts.append(config.rbs_spacer.upper())
            parts.append(dna.upper())
            # Gene-specific terminator (if not last gene, or if specified)
            if spec.terminator:
                parts.append(spec.terminator.upper())
        # Terminator after last gene
        if config.terminator:
            parts.append(config.terminator.upper())

    elif construct_type == "polycistronic_2A":
        # First gene's promoter
        if gene_specs and gene_specs[0].promoter:
            parts.append(gene_specs[0].promoter.upper())
        # Genes separated by 2A peptide DNA
        for i, (dna, spec) in enumerate(zip(gene_dnas, gene_specs)):
            if i > 0 and linker:
                # linker is the 2A peptide protein sequence; convert to DNA
                parts.append(linker.upper())
            elif i > 0:
                parts.append(dna.upper())
                continue
            parts.append(dna.upper())
            if spec.terminator and i == len(gene_dnas) - 1:
                parts.append(spec.terminator.upper())

    elif construct_type == "polycistronic_IRES":
        # First gene's promoter
        if gene_specs and gene_specs[0].promoter:
            parts.append(gene_specs[0].promoter.upper())
        # Genes separated by IRES DNA
        for i, (dna, spec) in enumerate(zip(gene_dnas, gene_specs)):
            parts.append(dna.upper())
            if i < len(gene_dnas) - 1 and linker:
                parts.append(linker.upper())  # IRES DNA sequence
            if spec.terminator and i == len(gene_dnas) - 1:
                parts.append(spec.terminator.upper())

    elif construct_type == "bidirectional":
        # Split genes into two halves: forward and reverse
        mid = len(gene_dnas) // 2
        forward_genes = gene_dnas[:mid]
        reverse_genes = gene_dnas[mid:]
        forward_specs = gene_specs[:mid]
        reverse_specs = gene_specs[mid:]

        # Forward half: promoter + genes
        if forward_specs and forward_specs[0].promoter:
            parts.append(forward_specs[0].promoter.upper())
        for i, dna in enumerate(forward_genes):
            if forward_specs[i].rbs:
                parts.append(forward_specs[i].rbs.upper())
            parts.append(dna.upper())

        # Bidirectional promoter (linker serves as the promoter region)
        if linker:
            parts.append(linker.upper())

        # Reverse half: promoter + genes
        if reverse_specs and reverse_specs[0].promoter:
            parts.append(reverse_specs[0].promoter.upper())
        for i, dna in enumerate(reverse_genes):
            if reverse_specs[i].rbs:
                parts.append(reverse_specs[i].rbs.upper())
            parts.append(dna.upper())

    else:  # "custom" or fallback
        # Simple concatenation with linkers
        for i, (dna, spec) in enumerate(zip(gene_dnas, gene_specs)):
            if spec.promoter:
                parts.append(spec.promoter.upper())
            if spec.rbs:
                parts.append(spec.rbs.upper())
            parts.append(dna.upper())
            if spec.terminator:
                parts.append(spec.terminator.upper())
            # Add linker between genes
            if i < len(gene_dnas) - 1 and linker:
                parts.append(linker.upper())

    return "".join(parts)


def _generate_genbank_multigene(
    full_dna: str,
    gene_dnas: list[str],
    gene_specs: list[GeneSpec],
    gene_results: list[OptimizationResult],
    organism: str,
    construct_type: str,
    linker: str = "",
    operon_config: OperonConfig | None = None,
) -> str:
    """Generate a GenBank record for the multi-gene construct.

    Annotates each gene as a separate CDS feature with its position in the
    full construct, plus regulatory features for promoters, RBS, and terminators.

    Args:
        full_dna: The full assembled DNA construct.
        gene_dnas: Individual gene DNA sequences.
        gene_specs: Gene specifications.
        gene_results: Per-gene optimization results.
        organism: Target organism.
        construct_type: Type of construct.
        linker: Linker used between genes.
        operon_config: Operon configuration (if applicable).

    Returns:
        GenBank-format string.
    """
    from .export import (
        export_genbank,
        _generate_accession,
        _format_genbank_header,
        _format_genbank_sequence,
        _get_taxonomy,
        GENBANK_MAX_LINE,
    )
    from . import __version__

    seq = full_dna.upper()
    gc = gc_content(seq)
    acc = _generate_accession()
    now = datetime.now(timezone.utc)
    date_str = now.strftime("%d-%b-%Y").upper()
    organism_display = organism.replace("_", " ")

    # Determine construct definition
    gene_names = [s.name or f"gene_{i+1}" for i, s in enumerate(gene_specs)]
    construct_label = {
        "operon": "Operon",
        "polycistronic_2A": "Polycistronic 2A construct",
        "polycistronic_IRES": "Polycistronic IRES construct",
        "bidirectional": "Bidirectional promoter construct",
        "custom": "Multi-gene construct",
    }.get(construct_type, "Multi-gene construct")

    definition = f"{construct_label}: {', '.join(gene_names)}"

    # ── LOCUS / DEFINITION / ACCESSION ──
    lines: list[str] = []
    locus = f"MULTIGENE_{len(gene_dnas)}GENES"[:16].upper()
    length_str = f"{len(seq)} bp"
    lines.append(
        f"LOCUS       {locus:<16} {length_str:>12}   DNA     linear   SYN"
    )
    lines.append(f"DEFINITION  {definition}.")
    lines.append(f"ACCESSION   {acc}")
    lines.append(f"VERSION     {acc}.1")
    lines.append("KEYWORDS    BioCompiler; codon-optimized; synthetic gene; multi-gene construct.")

    # ── SOURCE / ORGANISM ──
    taxonomy = _get_taxonomy(organism)
    lines.append(f"SOURCE      {organism_display}")
    lines.append(f"  ORGANISM  {organism_display}")
    for i in range(0, len(taxonomy), 70):
        chunk = taxonomy[i:i + 70]
        indent = "            " if i > 0 else "  "
        lines.append(f"{indent}{chunk}")
    if not taxonomy.endswith("."):
        lines[-1] += "."

    # ── COMMENT ──
    comments: list[str] = []
    comments.append("Designed and verified by BioCompiler — Machine-Verified Gene Design")
    comments.append(f"Version: {__version__}")
    comments.append(f"Construct type: {construct_type}")
    comments.append(f"Number of genes: {len(gene_dnas)}")
    comments.append(f"Total GC content: {gc:.4f}")
    comments.append(f"Total length: {len(seq)} bp")
    for i, (spec, result) in enumerate(zip(gene_specs, gene_results)):
        gname = spec.name or f"gene_{i+1}"
        comments.append(f"  Gene {i+1} ({gname}): CAI={result.cai:.4f}, "
                       f"GC={result.gc_content:.4f}, {len(result.sequence)} bp")
    if linker:
        comments.append(f"Linker: {linker[:50]}{'...' if len(linker) > 50 else ''}")

    lines.append("COMMENT     " + comments[0])
    for c in comments[1:]:
        lines.append("            " + c)

    # ── FEATURES ──
    lines.append("FEATURES             Location/Qualifiers")

    # Locate each gene's DNA within the full construct to annotate positions.
    # We search for each gene's DNA in the assembled construct.
    gene_positions: list[tuple[int, int]] = []  # (start, end) 0-based, half-open
    search_start = 0
    for gene_dna in gene_dnas:
        gene_upper = gene_dna.upper()
        pos = seq.find(gene_upper, search_start)
        if pos == -1:
            # Fallback: approximate position based on running length
            # This shouldn't happen with correct assembly
            pos = search_start
        gene_positions.append((pos, pos + len(gene_upper)))
        search_start = pos + len(gene_upper)

    # Annotate regulatory elements and CDS for each gene
    for i, (spec, result, (start, end)) in enumerate(
        zip(gene_specs, gene_results, gene_positions)
    ):
        gname = spec.name or f"gene_{i+1}"

        # Gene feature spanning the CDS
        lines.append(f"     gene            {start + 1}..{end}")
        lines.append(f'                     /gene="{gname}"')
        lines.append(f'                     /note="Gene {i+1} of {len(gene_dnas)} in {construct_type} construct"')

        # CDS feature
        lines.append(f"     CDS             {start + 1}..{end}")
        lines.append(f'                     /gene="{gname}"')
        lines.append(f'                     /organism="{organism_display}"')
        lines.append(f'                     /note="Codon-optimized by BioCompiler v{__version__}"')
        lines.append(f'                     /cai="{result.cai:.4f}"')
        lines.append(f'                     /codon_start=1')
        lines.append(f'                     /transl_table=1')

        # Protein translation
        protein = spec.protein
        if protein:
            prot_chunks = [protein[j:j + 40] for j in range(0, len(protein), 40)]
            lines.append(f'                     /translation="{prot_chunks[0]}"')
            for chunk in prot_chunks[1:]:
                lines.append(f'                     "{chunk}"')

        # Promoter feature (if specified and this is the first gene, or per-gene)
        if spec.promoter:
            # Find promoter position in the full construct
            prom_upper = spec.promoter.upper()
            prom_pos = seq.find(prom_upper, max(0, start - len(prom_upper) - 20))
            if prom_pos >= 0:
                lines.append(f"     promoter        {prom_pos + 1}..{prom_pos + len(prom_upper)}")
                lines.append(f'                     /gene="{gname}"')
                lines.append(f'                     /note="Promoter for {gname}"')

        # RBS feature (prokaryotic)
        if spec.rbs:
            rbs_upper = spec.rbs.upper()
            rbs_pos = seq.find(rbs_upper, max(0, start - len(rbs_upper) - 20))
            if rbs_pos >= 0:
                lines.append(f"     regulatory      {rbs_pos + 1}..{rbs_pos + len(rbs_upper)}")
                lines.append(f'                     /regulatory_class="ribosome_binding_site"')
                lines.append(f'                     /gene="{gname}"')
                lines.append(f'                     /note="RBS for {gname}"')

        # Terminator feature
        if spec.terminator:
            term_upper = spec.terminator.upper()
            term_pos = seq.find(term_upper, end - 10)
            if term_pos >= 0:
                lines.append(f"     terminator      {term_pos + 1}..{term_pos + len(term_upper)}")
                lines.append(f'                     /gene="{gname}"')
                lines.append(f'                     /note="Terminator for {gname}"')

    # Add regulatory feature for the whole construct
    lines.append(f"     regulatory      1..{len(seq)}")
    lines.append(f'                     /regulatory_class="codon_optimization"')
    lines.append(f'                     /note="Multi-gene {construct_type} construct optimized by BioCompiler v{__version__}"')

    # ── ORIGIN ──
    lines.append("ORIGIN")
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60]
        groups = [chunk[j:j + 10] for j in range(0, len(chunk), 10)]
        line_num = i + 1
        lines.append(f"{line_num:>9} {' '.join(groups)}")
    lines.append("//")

    return "\n".join(lines)


# ─── Public API ─────────────────────────────────────────────────────


def optimize_multigene(
    genes: list[GeneSpec],
    linker: str = "",
    organism: str = "",
    constraints: dict | None = None,
) -> MultiGeneResult:
    """Optimize multiple genes as a single DNA construct.

    This is the primary entry point for multi-gene construct design.  Each gene
    is independently codon-optimized for the specified organism, then assembled
    into a contiguous DNA sequence with the specified linkers and regulatory
    elements.

    The construct type is automatically inferred from the gene specifications:

    - If all genes have RBS but only the first has a promoter → **operon**
    - If a 2A peptide protein sequence is provided as linker → **polycistronic_2A**
    - If a DNA sequence > 100 bp is provided as linker → **polycistronic_IRES**
    - If genes have alternating promoters facing opposite directions → **bidirectional**
    - Otherwise → **custom**

    Supported construct types:

    **Prokaryotic operon**::

        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp", rbs="AGGAGG"),
                GeneSpec(protein="MFK...", name="kanR", rbs="AGGAGG"),
            ],
            organism="Escherichia_coli",
        )

    **Eukaryotic polycistronic (2A peptides)**::

        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGEELFTG", name="GFP"),
                GeneSpec(protein="MFK...", name="RFP"),
            ],
            linker=LINKER_2A_P2A,  # 2A peptide protein sequence
            organism="Homo_sapiens",
        )

    **Eukaryotic polycistronic (IRES)**::

        ires_seq = "CCCCCTCCCCC..."  # Your IRES DNA sequence
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGEELFTG", name="GFP"),
                GeneSpec(protein="MFK...", name="RFP"),
            ],
            linker=ires_seq,
            organism="Homo_sapiens",
        )

    Args:
        genes: List of gene specifications.  Must contain at least one gene.
        linker: Sequence between genes.  Interpretation depends on construct type:

            - For 2A polycistronic: the 2A peptide **protein** sequence
              (e.g., ``LINKER_2A_P2A``).  It will be codon-optimized and
              inserted between genes.  The last amino acid of the upstream
              gene and the first of the 2A peptide form the cleavage site.
            - For IRES polycistronic: the IRES **DNA** sequence.
            - For custom: any DNA linker sequence.

        organism: Target organism for codon optimization.  Used as default
            for genes that don't specify their own organism.
        constraints: Optional dict of constraint parameters forwarded to
            :func:`~biocompiler.optimization.optimize_sequence`.  Supported
            keys: ``gc_lo``, ``gc_hi``, ``cai_threshold``, ``restriction_sites``.

    Returns:
        A :class:`MultiGeneResult` containing per-gene optimization results,
        the full assembled DNA, and a GenBank-format export.

    Raises:
        ValueError: If ``genes`` is empty or ``organism`` is not provided.
        InvalidProteinError: If any gene has an invalid protein sequence.
        UnsupportedOrganismError: If the organism is not supported.

    Examples:
        Simple two-gene operon in E. coli::

            from biocompiler.multigene import optimize_multigene, GeneSpec

            result = optimize_multigene(
                genes=[
                    GeneSpec(protein="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDA",
                             name="gfp", rbs="AGGAGG"),
                    GeneSpec(protein="MDDRLEAIAGMTRLLRALRRKL",
                             name="tetR", rbs="AGGAGG"),
                ],
                organism="Escherichia_coli",
            )
            print(f"Full construct: {result.full_dna[:80]}...")
            print(f"Total length: {result.total_length} bp")
            print(f"GC content: {result.gc_content:.4f}")
    """
    if not genes:
        raise ValueError("At least one GeneSpec must be provided")
    if not organism and not any(g.organism for g in genes):
        raise ValueError(
            "Organism must be specified either at the construct level "
            "or in at least one GeneSpec"
        )

    # Resolve default organism
    default_organism = organism or genes[0].organism
    resolved_organism = resolve_organism(default_organism, strict=False)

    # ── Step 1: Optimize each gene independently ──
    gene_results: list[OptimizationResult] = []
    gene_dnas: list[str] = []

    for i, spec in enumerate(genes):
        logger.info("Optimizing gene %d/%d: %s", i + 1, len(genes),
                    spec.name or f"gene_{i+1}")
        result = _optimize_gene(spec, resolved_organism, constraints)
        gene_results.append(result)
        gene_dnas.append(result.sequence)

    # ── Step 2: Determine construct type ──
    construct_type = _infer_construct_type(genes, linker, resolved_organism)
    logger.info("Inferred construct type: %s", construct_type)

    # ── Step 3: Prepare linkers ──
    # For 2A polycistronic, linker is a protein sequence; convert to DNA
    effective_linker = linker
    if construct_type == "polycistronic_2A" and linker:
        # The linker is a 2A peptide protein sequence — optimize it
        effective_linker = _protein_to_2a_dna(linker, resolved_organism)

    # ── Step 4: Assemble the full construct ──
    operon_config: OperonConfig | None = None
    if construct_type == "operon":
        operon_config = OperonConfig(
            promoter=genes[0].promoter if genes[0].promoter else "",
            rbs_per_gene=RBS_STRONG,
            rbs_spacer=RBS_SPACER_DEFAULT,
            terminator=genes[-1].terminator if genes[-1].terminator else "",
        )

    bidirectional = construct_type == "bidirectional"

    full_dna = _assemble_construct(
        gene_dnas=gene_dnas,
        gene_specs=genes,
        linker=effective_linker,
        construct_type=construct_type,
        operon_config=operon_config,
        bidirectional=bidirectional,
    )

    # ── Step 5: Generate GenBank export ──
    genbank_export = _generate_genbank_multigene(
        full_dna=full_dna,
        gene_dnas=gene_dnas,
        gene_specs=genes,
        gene_results=gene_results,
        organism=resolved_organism,
        construct_type=construct_type,
        linker=effective_linker,
        operon_config=operon_config,
    )

    # ── Step 6: Build result ──
    gc = gc_content(full_dna.upper())

    return MultiGeneResult(
        genes=gene_results,
        full_dna=full_dna.upper(),
        total_length=len(full_dna),
        genbank_export=genbank_export,
        construct_type=construct_type,
        organism=resolved_organism,
        gc_content=gc,
    )


def _infer_construct_type(
    genes: list[GeneSpec],
    linker: str,
    organism: str,
) -> str:
    """Infer the construct type from gene specifications and linker.

    Decision logic:
    1. If organism is prokaryotic and multiple genes have RBS → operon
    2. If linker matches a known 2A peptide protein sequence → polycistronic_2A
    3. If linker is a long DNA sequence (>100 bp, no invalid chars) → polycistronic_IRES
    4. If genes have alternating promoter directions → bidirectional
    5. Otherwise → custom

    Args:
        genes: List of gene specifications.
        linker: The linker sequence provided by the user.
        organism: Target organism.

    Returns:
        A string construct type identifier.
    """
    # Check for 2A peptide linker
    known_2a_proteins = {LINKER_2A_P2A, LINKER_2A_T2A, LINKER_2A_E2A, LINKER_2A_F2A}
    if linker and linker.upper() in {p.upper() for p in known_2a_proteins}:
        return "polycistronic_2A"

    # Check for IRES (long DNA sequence)
    if linker and len(linker) > 100:
        valid_dna = set("ACGT")
        if set(linker.upper()) <= valid_dna:
            return "polycistronic_IRES"

    # Check for prokaryotic operon
    is_prokaryote = not is_eukaryotic_organism(organism)
    if is_prokaryote and len(genes) > 1:
        genes_with_rbs = sum(1 for g in genes if g.rbs)
        if genes_with_rbs >= len(genes) - 1:  # Allow one gene without RBS
            return "operon"

    # Check for bidirectional promoter layout
    if len(genes) >= 4:  # Need at least 2 genes per direction
        # Heuristic: if first half and second half each have a promoter
        mid = len(genes) // 2
        first_half_promoters = sum(1 for g in genes[:mid] if g.promoter)
        second_half_promoters = sum(1 for g in genes[mid:] if g.promoter)
        if first_half_promoters > 0 and second_half_promoters > 0:
            return "bidirectional"

    # If prokaryotic with no RBS specified, still default to operon for
    # multi-gene constructs (user can add RBS later)
    if is_prokaryote and len(genes) > 1:
        return "operon"

    return "custom"


def optimize_operon(
    genes: list[GeneSpec],
    organism: str = "Escherichia_coli",
    promoter: str = "",
    rbs_per_gene: str = RBS_STRONG,
    rbs_spacer: str = RBS_SPACER_DEFAULT,
    terminator: str = "",
    constraints: dict | None = None,
) -> MultiGeneResult:
    """Convenience function for prokaryotic operon optimization.

    This is a simplified interface for the most common multi-gene construct
    in prokaryotes: an operon with a single promoter, individual RBS sequences
    for each gene, and an optional terminator.

    The function sets appropriate defaults and delegates to
    :func:`optimize_multigene` with ``construct_type="operon"``.

    Args:
        genes: List of gene specifications.  Each gene's ``rbs`` field can
            override the default ``rbs_per_gene`` for that specific gene.
        organism: Target organism (default: ``"Escherichia_coli"``).
        promoter: Promoter sequence placed upstream of the first gene.
        rbs_per_gene: Default RBS sequence for genes that don't specify one.
        rbs_spacer: DNA spacer between the RBS and the start codon.
        terminator: Terminator sequence placed after the last gene.
        constraints: Optional constraint parameters for optimization.

    Returns:
        A :class:`MultiGeneResult` with ``construct_type="operon"``.

    Raises:
        ValueError: If no genes are provided.

    Examples:
        Three-gene operon in E. coli::

            from biocompiler.multigene import optimize_operon, GeneSpec

            result = optimize_operon(
                genes=[
                    GeneSpec(protein="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDA",
                             name="gfp"),
                    GeneSpec(protein="MDDRLEAIAGMTRLLRALRRKL",
                             name="tetR"),
                    GeneSpec(protein="MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAP",
                             name="lacZ"),
                ],
                organism="Escherichia_coli",
                promoter="TTGACA",
            )
            print(f"Operon: {result.total_length} bp, GC={result.gc_content:.4f}")
    """
    if not genes:
        raise ValueError("At least one GeneSpec must be provided")

    # Set promoter on first gene if provided
    resolved_genes = []
    for i, g in enumerate(genes):
        # Create a modified GeneSpec with promoter/terminator/RBS set appropriately
        spec = GeneSpec(
            protein=g.protein,
            organism=g.organism or organism,
            name=g.name,
            promoter=g.promoter or (promoter if i == 0 else ""),
            rbs=g.rbs or rbs_per_gene,
            terminator=g.terminator or (terminator if i == len(genes) - 1 else ""),
        )
        resolved_genes.append(spec)

    return optimize_multigene(
        genes=resolved_genes,
        linker="",  # No linker for operons
        organism=organism,
        constraints=constraints,
    )
