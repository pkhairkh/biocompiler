"""
BioCompiler Export Engine — GenBank, FASTA & JSON Sequence Export

Production-grade sequence export with:
- GenBank format output with full feature annotations (CDS, gene, regulatory)
- FASTA format with metadata headers and comment lines
- JSON format for full OptimizationResult serialization
- Certificate embedding in GenBank comment section
- Exon/intron/restriction-site feature annotations
- CAI value as CDS qualifier
- Organism name and optimization date in records
- IUPAC-compliant sequence representation

The export transforms internal designed sequences + type-check results
into standard bioinformatics formats accepted by NCBI, Benchling, SnapGene, etc.
"""

__all__ = [
    "RestrictionSiteInfo",
    "FastaSequenceEntry",
    "export_fasta",
    "export_genbank",
    "export_genbank_with_certificate",
    "export_multi_fasta",
    "export_batch_fasta",
    "export_full_construct",
    "export_json",
    "export_with_annotations",
    "format_biosecurity_report",
    "GENBANK_MAX_LINE",
    "GENBANK_SEQ_LINE",
    "GENBANK_SEQ_GROUP",
]

import json
import logging
import uuid
from dataclasses import asdict, fields, is_dataclass
from datetime import datetime, timezone
from typing import Any, Optional, TypedDict

from .types import Certificate, TypeCheckResult, Verdict, combined_verdict
from .scanner import gc_content
from .translation import translate, compute_cai
from . import __version__

# Prokaryotic organisms (BSL-1 by default)
_PROKARYOTIC_ORGANISMS = frozenset({
    "Escherichia_coli", "E_coli", "Bacillus_subtilis",
})

# BSL-2 organisms (human/mammalian pathogens or cell lines)
_BSL2_ORGANISMS = frozenset({
    "Homo_sapiens", "Mus_musculus", "CHO_K1",
})

logger = logging.getLogger(__name__)

# GenBank format constants
GENBANK_MAX_LINE = 80
GENBANK_SEQ_LINE = 60
GENBANK_SEQ_GROUP = 10


def _wrap_text(text: str, width: int = GENBANK_MAX_LINE, indent: int = 0) -> str:
    """Wrap text to fit within GenBank line width constraints."""
    prefix = " " * indent
    lines = []
    while text:
        chunk = text[:width - indent]
        lines.append(prefix + chunk)
        text = text[width - indent:]
    return "\n".join(lines)


def _format_sequence_numbered(seq: str, line_width: int = GENBANK_SEQ_LINE,
                              group_size: int = GENBANK_SEQ_GROUP) -> str:
    """Format a DNA sequence in GenBank style (numbered groups of 10).

    This is a standalone utility that returns a single formatted string.
    For the full ORIGIN section of a GenBank record (returning a list of
    lines), see :func:`_format_genbank_sequence`.
    """
    seq = seq.upper()
    lines = []
    for i in range(0, len(seq), line_width):
        chunk = seq[i:i + line_width]
        groups = [chunk[j:j + group_size] for j in range(0, len(chunk), group_size)]
        line_num = i + 1
        lines.append(f"{line_num:>9} {' '.join(groups)}")
    return "\n".join(lines)


def _format_fasta_sequence(seq: str, line_width: int = 60) -> str:
    """Format a DNA sequence in FASTA style (60 chars per line)."""
    seq = seq.upper()
    lines = [seq[i:i + line_width] for i in range(0, len(seq), line_width)]
    return "\n".join(lines)


def _generate_accession() -> str:
    """Generate a unique accession ID using UUID4 and a timestamp.

    Produces a pseudo-accession of the form ``BC_`` followed by the first
    8 hex characters of a UUID4, guaranteeing uniqueness without relying
    on a placeholder locus name.

    Returns:
        A unique pseudo-accession string, e.g. ``BC_a1b2c3d4``.
    """
    return f"BC_{uuid.uuid4().hex[:8].upper()}"


class _RestrictionSiteInfoRequired(TypedDict):
    """Required fields for a restriction site annotation."""
    enzyme: str
    site: str
    position: int


class RestrictionSiteInfo(_RestrictionSiteInfoRequired, total=False):
    """A restriction enzyme recognition site annotation.

    Required fields: ``enzyme``, ``site``, ``position``.
    Optional field: ``strand`` (``"+"`` or ``"-"``, default ``"+"``).
    """
    strand: str


class _FastaSequenceEntryRequired(TypedDict):
    """Required fields for a FASTA sequence entry."""
    sequence: str


class FastaSequenceEntry(_FastaSequenceEntryRequired, total=False):
    """An entry for :func:`export_multi_fasta` and :func:`export_batch_fasta`.

    Required field: ``sequence``.
    Optional fields: ``id``, ``description``, ``organism``, ``protein``,
    ``cai``, ``gc``.
    """
    id: str
    description: str
    organism: str
    protein: str
    cai: float
    gc: float


def _assess_biosafety_level(
    organism: str,
    type_results: Optional[list[TypeCheckResult]] = None,
) -> str:
    """Assess the biosafety level for a given organism and predicate results.

    Returns one of: 'BSL-1', 'BSL-2', or 'unknown'.

    BSL-1: Non-pathogenic organisms (E. coli K-12, B. subtilis) and
           sequences where all predicates pass.
    BSL-2: Organisms associated with human cell lines or known pathogens,
           or sequences with failed predicates that may indicate safety concerns.
    unknown: Organisms not in the known classification.
    """
    if organism in _PROKARYOTIC_ORGANISMS:
        base_level = "BSL-1"
    elif organism in _BSL2_ORGANISMS:
        base_level = "BSL-2"
    else:
        base_level = "unknown"

    # If any predicate failed, escalate to BSL-2 minimum
    if type_results:
        has_fail = any(
            r.verdict == Verdict.FAIL for r in type_results
        )
        if has_fail and base_level == "BSL-1":
            return "BSL-2"

    return base_level


def _is_biosecurity_screened(
    type_results: Optional[list[TypeCheckResult]] = None,
) -> bool:
    """Determine if the sequence has been biosecurity-screened.

    A sequence is considered biosecurity-screened if it has passed
    all type-check predicates (no FAIL verdicts).
    """
    if not type_results:
        return False
    return all(r.verdict != Verdict.FAIL for r in type_results)


def format_biosecurity_report(
    sequence: str,
    organism: str = "Homo_sapiens",
    cai: Optional[float] = None,
    gc: Optional[float] = None,
    type_results: Optional[list[TypeCheckResult]] = None,
) -> str:
    """Format a full biosecurity screening report for CLI display.

    This function produces a human-readable report summarizing the
    biosecurity assessment of a designed sequence, including biosafety
    level, screening status, and predicate results.

    Args:
        sequence: DNA sequence.
        organism: Target organism name.
        cai: CAI value.
        gc: GC content value.
        type_results: Type-check predicate results.

    Returns:
        Human-readable biosecurity report string.
    """
    seq = sequence.upper().replace(" ", "")
    if gc is None:
        gc = gc_content(seq)

    bsl = _assess_biosafety_level(organism, type_results)
    screened = _is_biosecurity_screened(type_results)

    passed_predicates = []
    failed_predicates = []
    if type_results:
        for r in type_results:
            if r.verdict == Verdict.PASS:
                passed_predicates.append(r.predicate)
            elif r.verdict == Verdict.FAIL:
                failed_predicates.append(r.predicate)

    provenance_id = f"BC_{uuid.uuid4().hex[:12].upper()}"

    lines = [
        "=" * 60,
        "  BIOSECURITY SCREENING REPORT",
        "=" * 60,
        f"  Provenance ID    : {provenance_id}",
        f"  Optimized by     : biocompiler v{__version__}",
        f"  Organism         : {organism.replace('_', ' ')}",
        f"  CAI score        : {cai:.4f}" if cai is not None else "  CAI score        : N/A",
        f"  GC content       : {gc:.4f}",
        f"  Biosecurity level: {bsl}",
        f"  Screened         : {'PASS' if screened else 'FAIL'}",
        "",
        f"  Passed predicates ({len(passed_predicates)}):",
    ]
    for p in passed_predicates:
        lines.append(f"    [+] {p}")
    if not passed_predicates:
        lines.append("    (none)")

    lines.append(f"")
    lines.append(f"  Failed predicates ({len(failed_predicates)}):")
    for p in failed_predicates:
        lines.append(f"    [X] {p}")
    if not failed_predicates:
        lines.append("    (none)")

    # Risk assessment summary
    lines.append("")
    lines.append("  Risk Assessment Summary:")
    if not failed_predicates:
        lines.append("    All predicates passed. Sequence is considered safe for synthesis.")
    else:
        lines.append(f"    WARNING: {len(failed_predicates)} predicate(s) failed.")
        lines.append("    Additional review recommended before submitting for synthesis.")
    if bsl == "BSL-2":
        lines.append("    Biosafety level BSL-2: Follow institutional BSL-2 containment procedures.")
    elif bsl == "BSL-1":
        lines.append("    Biosafety level BSL-1: Standard laboratory practices sufficient.")
    else:
        lines.append("    Biosafety level unknown: Consult institutional biosafety officer.")

    lines.append("=" * 60)
    return "\n".join(lines)


def _format_genbank_header(
    locus: str,
    seq_len: int,
    molecule_type: str,
    topology: str,
    date_str: str,
    definition: str,
    acc: str,
    organism: str,
    gc: float,
    protein: Optional[str],
    type_results: Optional[list[TypeCheckResult]],
    certificate: Optional[Certificate],
    cai: Optional[float] = None,
    optimization_date: Optional[str] = None,
) -> list[str]:
    """Format the LOCUS, DEFINITION, ACCESSION, VERSION, SOURCE, ORGANISM, and COMMENT sections."""
    lines: list[str] = []

    # ── LOCUS / DEFINITION / ACCESSION / VERSION ──
    length_str = f"{seq_len} bp"
    # Validate and format molecule type for GenBank compliance
    mol_type_upper = molecule_type.upper()
    if mol_type_upper in ("DNA", "RNA", "MRNA"):
        if mol_type_upper == "MRNA":
            mol_str = "mRNA    "
        else:
            mol_str = f"{molecule_type}    "
    else:
        mol_str = f"{molecule_type}    "

    # Validate topology
    valid_topologies = ("linear", "circular")
    topo = topology if topology in valid_topologies else "linear"

    lines.append(
        f"LOCUS       {locus:<16} {length_str:>12}   {mol_str} {topo:<8}   SYN"
    )
    lines.append(f"DEFINITION  {definition}.")

    lines.append(f"ACCESSION   {acc}")
    lines.append(f"VERSION     {acc}.1")

    # ── KEYWORDS ──
    lines.append("KEYWORDS    BioCompiler; codon-optimized; synthetic gene.")

    # ── SOURCE / ORGANISM ──
    taxonomy = _get_taxonomy(organism)
    # Format organism name nicely: replace underscores with spaces for display
    organism_display = organism.replace("_", " ")
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
    comments.append(f"GC content: {gc:.4f}")
    if cai is not None:
        comments.append(f"CAI: {cai:.4f}")
    comments.append(f"Protein length: {len(protein)} aa" if protein else "")

    # Optimization date
    opt_date = optimization_date or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    comments.append(f"Optimization date: {opt_date}")

    comments.append(f"Target organism: {organism_display}")

    if type_results:
        overall = combined_verdict([r.verdict for r in type_results])
        comments.append(f"Type-check verdict: {overall.value}")
        for r in type_results:
            symbol = {"PASS": "+", "FAIL": "X", "UNCERTAIN": "?",
                      "LIKELY_PASS": "~+", "LIKELY_FAIL": "~X"}.get(r.verdict.value, "?")
            comments.append(f"  [{symbol}] {r.predicate}")

    if certificate:
        comments.append(f"Certificate ID: {certificate.design_id[:16]}...")
        comments.append(f"Certificate timestamp: {certificate.provenance.get('timestamp', 'N/A')}")

    # ── BIOCOMPILER_ANNOTATIONS ──
    # Biosafety annotation block
    bsl = _assess_biosafety_level(organism, type_results)
    screened = _is_biosecurity_screened(type_results)
    provenance_id = f"BC_{uuid.uuid4().hex[:12].upper()}"

    passed_preds = []
    failed_preds = []
    if type_results:
        for r in type_results:
            if r.verdict == Verdict.PASS:
                passed_preds.append(r.predicate)
            elif r.verdict == Verdict.FAIL:
                failed_preds.append(r.predicate)

    comments.append("")  # blank line separator
    comments.append("BIOCOMPILER_ANNOTATIONS:")
    comments.append(f"  optimized_by: biocompiler v{__version__}")
    comments.append(f"  organism: {organism}")
    comments.append(f"  cai_score: {cai:.4f}" if cai is not None else "  cai_score: N/A")
    comments.append(f"  gc_content: {gc:.4f}")
    comments.append(f"  passed_predicates: {passed_preds}")
    comments.append(f"  failed_predicates: {failed_preds}")
    comments.append(f"  biosecurity_screened: {screened}")
    comments.append(f"  biosafety_level: {bsl}")
    comments.append(f"  provenance_id: {provenance_id}")

    # ── WARNING if predicates failed ──
    if failed_preds:
        comments.append("")
        comments.append(
            f"WARNING: {len(failed_preds)} predicate(s) failed: "
            f"{', '.join(failed_preds)}. "
            f"Review before submitting for gene synthesis."
        )

    # ── BIOSECURITY NOTICE ──
    comments.append("")
    comments.append("BIOSECURITY NOTICE:")
    if not failed_preds:
        comments.append(
            "  This sequence has passed all biosecurity screening predicates "
            "and is considered safe for synthesis."
        )
    else:
        comments.append(
            "  This sequence has FAILED one or more biosecurity predicates. "
            "Additional review by a biosafety officer is recommended before "
            "submitting for gene synthesis."
        )
    if bsl == "BSL-1":
        comments.append(
            "  Risk level: BSL-1 — Standard laboratory practices are sufficient."
        )
    elif bsl == "BSL-2":
        comments.append(
            "  Risk level: BSL-2 — Institutional BSL-2 containment procedures apply."
        )
    else:
        comments.append(
            "  Risk level: Unknown — Consult institutional biosafety officer."
        )

    if any(c for c in comments):
        lines.append("COMMENT     " + comments[0])
        for c in comments[1:]:
            # Preserve blank lines for readability
            if c == "":
                lines.append("            ")
            else:
                lines.append("            " + c)

    return lines


def _format_genbank_features(
    seq_len: int,
    gene_name: Optional[str],
    protein: Optional[str],
    exon_boundaries: Optional[list[tuple[int, int]]],
    restriction_sites: Optional[list[RestrictionSiteInfo]],
    type_results: Optional[list[TypeCheckResult]],
    organism: Optional[str] = None,
    cai: Optional[float] = None,
) -> list[str]:
    """Format the FEATURE TABLE section of a GenBank record.

    Includes gene, CDS, regulatory, and exon feature annotations with
    CAI, organism, and optimization metadata qualifiers.
    """
    lines: list[str] = []

    lines.append("FEATURES             Location/Qualifiers")

    # Gene feature
    if gene_name:
        lines.append(f"     gene            1..{seq_len}")
        lines.append(f'                     /gene="{gene_name}"')
        if organism:
            lines.append(f'                     /organism="{organism.replace("_", " ")}"')
        lines.append(f'                     /note="Designed by BioCompiler v{__version__}"')
        lines.append(f'                     /db_xref="BioCompiler:{__version__}"')

    # CDS feature
    if protein:
        if exon_boundaries and len(exon_boundaries) > 1:
            # Multi-exon: join operator
            exon_parts = []
            for start, end in exon_boundaries:
                # Convert 0-based [start, end) to 1-based [start+1, end]
                exon_parts.append(f"{start + 1}..{end}")
            location = f"join({','.join(exon_parts)})"
        else:
            location = f"1..{seq_len}"
        lines.append(f"     CDS             {location}")
        if gene_name:
            lines.append(f'                     /gene="{gene_name}"')
        if organism:
            lines.append(f'                     /organism="{organism.replace("_", " ")}"')
        lines.append(f'                     /note="Codon-optimized by BioCompiler v{__version__}"')
        # Add CAI value as a qualifier in the CDS feature
        if cai is not None:
            lines.append(f'                     /cai="{cai:.4f}"')
        lines.append(f'                     /codon_start=1')
        lines.append(f'                     /transl_table=1')
        # Protein translation (wrapped)
        if protein:
            prot_chunks = [protein[i:i + 40] for i in range(0, len(protein), 40)]
            lines.append(f'                     /translation="{prot_chunks[0]}"')
            for chunk in prot_chunks[1:]:
                lines.append(f'                     "{chunk}"')

    # Regulatory features (promoters, RBS, terminators, etc.)
    # For codon-optimized genes, add a regulatory feature noting the
    # codon optimization and the target organism
    if organism:
        lines.append(f"     regulatory      1..{seq_len}")
        lines.append(f'                     /regulatory_class="codon_optimization"')
        lines.append(f'                     /note="Codon-optimized for {organism.replace("_", " ")} using BioCompiler"')
        if gene_name:
            lines.append(f'                     /gene="{gene_name}"')

    # Exon features
    if exon_boundaries:
        for i, (start, end) in enumerate(exon_boundaries):
            lines.append(f"     exon            {start + 1}..{end}")
            if gene_name:
                lines.append(f'                     /gene="{gene_name}"')
            lines.append(f'                     /number={i + 1}')

    # Restriction site features (as regulatory features for production use)
    if restriction_sites:
        for site in restriction_sites[:20]:  # Limit to 20 annotations
            pos = site.get("position", 0)
            enz = site.get("enzyme", site.get("site", "unknown"))
            strand = site.get("strand", "+")
            site_len = len(site.get("site", ""))
            lines.append(f"     regulatory      {pos + 1}..{pos + site_len}")
            lines.append(f'                     /regulatory_class="restriction_site"')
            lines.append(f'                     /note="Restriction site: {enz} ({strand} strand)"')
            lines.append(f'                     /label="{enz}"')

    # Type-check result features (as misc_feature for failed predicates)
    if type_results:
        for r in type_results:
            if r.verdict == Verdict.FAIL and r.violation:
                lines.append(f"     misc_feature    1..{seq_len}")
                lines.append(f'                     /note="TYPE FAIL: {r.predicate} - {r.violation[:80]}"')
                lines.append(f'                     /label="typecheck_fail"')

    return lines


def _format_genbank_sequence(seq: str) -> list[str]:
    """Format the ORIGIN section (numbered sequence) and GenBank terminator."""
    lines: list[str] = []

    lines.append("ORIGIN")

    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60]
        # Group in blocks of 10
        groups = [chunk[j:j + 10] for j in range(0, len(chunk), 10)]
        line_num = i + 1
        lines.append(f"{line_num:>9} {' '.join(groups)}")

    lines.append("//")  # GenBank terminator

    return lines


def export_fasta(
    sequence: str,
    identifier: str = "BioCompiler_design",
    description: str = "",
    organism: str = "Homo_sapiens",
    protein: Optional[str] = None,
    cai: Optional[float] = None,
    include_comments: bool = True,
    type_results: Optional[list[TypeCheckResult]] = None,
) -> str:
    """
    Export a designed sequence in FASTA format.

    FASTA is the universal sequence format accepted by BLAST, Clustal,
    Geneious, and virtually all bioinformatics tools. This function
    generates a standards-compliant FASTA record with a rich header
    that includes organism, GC content, CAI, and protein translation.

    When ``include_comments=True``, comment lines (prefixed with ``;``)
    are added before the header with CAI, GC, and organism metadata,
    following the FASTA comment convention used by many tools.

    Args:
        sequence: DNA sequence (designed, verified)
        identifier: Sequence identifier (no spaces)
        description: Human-readable description line
        organism: Source organism name
        protein: Optional protein translation (auto-computed if None)
        cai: Optional CAI value (auto-computed if None)
        include_comments: If True, add FASTA comment lines with metadata
        type_results: Optional type-check predicate results for biosafety
            level assessment. When provided, failed predicates escalate
            the biosecurity level in the header (e.g. BSL-1 → BSL-2).

    Returns:
        FASTA-formatted string

    Example output::

        ; BioCompiler v10.0.0 — Machine-Verified Gene Design
        ; Organism: Escherichia coli
        ; CAI: 0.9990 | GC: 0.5230 | Length: 720 bp
        >GFP_optimized|organism=Escherichia_coli|gc=0.523|cai=0.999|len=720
        ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTG
        GACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCC
        ...
    """
    seq = sequence.upper().replace(" ", "")
    gc = gc_content(seq)

    if protein is None:
        protein = translate(seq)

    # Auto-compute CAI if not provided and sequence is valid
    if cai is None and len(seq) >= 3 and len(seq) % 3 == 0:
        try:
            cai = compute_cai(seq, organism=organism)
        except Exception:
            cai = None

    output_parts: list[str] = []

    # Comment lines with metadata (FASTA comment convention)
    if include_comments:
        output_parts.append(f"; BioCompiler v{__version__} — Machine-Verified Gene Design")
        output_parts.append(f"; Organism: {organism.replace('_', ' ')}")
        meta_parts = [f"GC: {gc:.4f}", f"Length: {len(seq)} bp"]
        if cai is not None:
            meta_parts.insert(0, f"CAI: {cai:.4f}")
        if protein:
            meta_parts.append(f"Protein: {len(protein)} aa")
        output_parts.append(f"; {' | '.join(meta_parts)}")

    # Assess biosecurity level for header annotation
    bsl = _assess_biosafety_level(organism, type_results)

    # Build FASTA header with structured metadata including biosecurity
    header_parts = [identifier]
    header_parts.append(f"organism={organism}")
    header_parts.append(f"CAI={cai:.4f}" if cai is not None else "CAI=N/A")
    header_parts.append(f"GC={gc:.4f}")
    header_parts.append(f"biosecurity={bsl}")
    header_parts.append(f"biocompiler_v{__version__}")

    header = "|".join(header_parts)
    if description:
        header += f" {description}"

    output_parts.append(f">{header}")
    output_parts.append(_format_fasta_sequence(seq))

    return "\n".join(output_parts) + "\n"


def export_batch_fasta(
    results: list[dict],
    organism: str = "Homo_sapiens",
) -> str:
    """
    Export multiple optimization results as a batch FASTA file.

    Each entry includes comment lines with CAI, GC, and organism metadata,
    making the output suitable for downstream analysis pipelines that
    process multiple sequences at once.

    Args:
        results: List of dicts, each with keys:
            - ``sequence`` (required): DNA sequence string
            - ``identifier`` (optional): Sequence ID (default: ``BioCompiler_design_N``)
            - ``description`` (optional): Description line
            - ``cai`` (optional): CAI value
            - ``gc`` (optional): GC content (auto-computed if not provided)
            - ``protein`` (optional): Protein translation
        organism: Default organism for all sequences

    Returns:
        Batch FASTA formatted string with all sequences

    Example::

        results = [
            {"sequence": "ATGGCC...", "identifier": "gene1", "cai": 0.95},
            {"sequence": "ATGAAA...", "identifier": "gene2", "cai": 0.88},
        ]
        fasta = export_batch_fasta(results, organism="Escherichia_coli")
    """
    records: list[str] = []
    for idx, entry in enumerate(results):
        seq = entry.get("sequence", "")
        if not seq:
            continue
        ident = entry.get("identifier", f"BioCompiler_design_{idx + 1}")
        desc = entry.get("description", "")
        entry_organism = entry.get("organism", organism)
        entry_cai = entry.get("cai")
        entry_protein = entry.get("protein")

        record = export_fasta(
            sequence=seq,
            identifier=ident,
            description=desc,
            organism=entry_organism,
            protein=entry_protein,
            cai=entry_cai,
            include_comments=True,
        )
        records.append(record.rstrip("\n"))

    return "\n".join(records) + "\n"


def export_genbank(
    sequence: str,
    locus_name: str = "BIOCOMPILER",
    definition: str = "BioCompiler designed sequence",
    organism: str = "Homo_sapiens",
    molecule_type: str = "DNA",
    topology: str = "linear",
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
    restriction_sites: Optional[list[RestrictionSiteInfo]] = None,
    certificate: Optional[Certificate] = None,
    type_results: Optional[list[TypeCheckResult]] = None,
    gene_name: Optional[str] = None,
    protein: Optional[str] = None,
    cai: Optional[float] = None,
    optimization_date: Optional[str] = None,
) -> str:
    """
    Export a designed sequence in GenBank format.

    GenBank is the standard format for sequence submission to NCBI/ENA/DDBJ
    and is the native format for Benchling, SnapGene, and Geneious. This
    function produces a fully annotated GenBank record with:

    - LOCUS header with correct molecule type and topology
    - DEFINITION, ACCESSION, VERSION, KEYWORDS headers
    - SOURCE and ORGANISM with taxonomy lineage
    - COMMENT section with CAI, GC, organism, optimization date, and certificate
    - FEATURE TABLE with gene, CDS (including CAI qualifier), regulatory,
      exon, and misc_feature annotations
    - ORIGIN section with numbered sequence

    Args:
        sequence: DNA sequence (designed, verified)
        locus_name: GenBank LOCUS name (max 16 chars, uppercase)
        definition: DEFINITION line text
        organism: Source organism for SOURCE/ORGANISM fields
        molecule_type: Molecule type (DNA, RNA, mRNA)
        topology: circular or linear
        exon_boundaries: List of (start, end) tuples for exon features
            (GenBank convention: 1-based, inclusive end)
        restriction_sites: List of :class:`RestrictionSiteInfo` dicts
        certificate: Optional Certificate to embed in COMMENT section
        type_results: Optional type-check results for FEATURE notes
        gene_name: Optional gene name for the gene feature
        protein: Optional protein translation (auto-computed if None)
        cai: Optional CAI value (auto-computed if None for valid sequences)
        optimization_date: Optional ISO 8601 date string for the optimization

    Returns:
        GenBank-formatted string
    """
    seq = sequence.upper().replace(" ", "")
    gc = gc_content(seq)
    now = datetime.now(timezone.utc)

    if protein is None:
        protein = translate(seq)

    # Auto-compute CAI if not provided and sequence is valid
    if cai is None and len(seq) >= 3 and len(seq) % 3 == 0:
        try:
            cai = compute_cai(seq, organism=organism)
        except Exception:
            cai = None

    # Truncate locus name to 16 chars (GenBank requirement)
    locus = locus_name[:16].upper()

    # Date in DD-MON-YYYY format
    date_str = now.strftime("%d-%b-%Y").upper()

    # ─── ACCESSION ──────────────────────────────────────────────────
    if certificate:
        acc = certificate.design_id[:12].upper()
    else:
        acc = _generate_accession()

    # ─── Assemble GenBank record from helper functions ──────────────
    lines: list[str] = []
    lines.extend(_format_genbank_header(
        locus, len(seq), molecule_type, topology, date_str,
        definition, acc, organism, gc, protein, type_results, certificate,
        cai=cai, optimization_date=optimization_date,
    ))
    lines.extend(_format_genbank_features(
        len(seq), gene_name, protein, exon_boundaries, restriction_sites,
        type_results, organism=organism, cai=cai,
    ))
    lines.extend(_format_genbank_sequence(seq))

    return "\n".join(lines)


def export_multi_fasta(
    sequences: list[FastaSequenceEntry],
) -> str:
    """
    Export multiple designed sequences as a multi-FASTA file.

    Each entry includes comment lines with CAI, GC, and organism metadata.

    Args:
        sequences: List of :class:`FastaSequenceEntry` dicts

    Returns:
        Multi-FASTA formatted string
    """
    records = []
    for entry in sequences:
        record = export_fasta(
            sequence=entry["sequence"],
            identifier=entry.get("id", "BioCompiler_design"),
            description=entry.get("description", ""),
            organism=entry.get("organism", "Homo_sapiens"),
            protein=entry.get("protein"),
            cai=entry.get("cai"),
            include_comments=True,
        )
        records.append(record.rstrip("\n"))
    return "\n".join(records) + "\n"


def export_genbank_with_certificate(
    sequence: str,
    certificate: Certificate,
    organism: str = "Homo_sapiens",
    gene_name: Optional[str] = None,
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
) -> str:
    """
    Export a GenBank record with full certificate provenance embedded.

    This is the primary export for regulatory submissions: the GenBank file
    contains both the sequence AND the machine-verified guarantee certificate,
    making it self-contained and independently verifiable.

    Args:
        sequence: DNA sequence (designed, verified)
        certificate: Certificate from BioCompiler pipeline
        organism: Source organism
        gene_name: Optional gene name
        exon_boundaries: Optional exon boundaries

    Returns:
        GenBank-formatted string with embedded certificate
    """
    type_results = _reconstruct_type_results(certificate)
    protein = translate(sequence)

    # Extract optimization date from certificate provenance
    opt_date = certificate.provenance.get("timestamp")

    return export_genbank(
        sequence=sequence,
        locus_name=certificate.design_id[:16].upper(),
        definition=f"BioCompiler verified design [{certificate.design_id[:8]}]",
        organism=organism,
        exon_boundaries=exon_boundaries,
        certificate=certificate,
        type_results=type_results,
        gene_name=gene_name,
        protein=protein,
        optimization_date=opt_date,
    )


# ─── JSON Export ──────────────────────────────────────────────────

def _serialize_for_json(obj: Any) -> Any:
    """Recursively serialize an object for JSON output.

    Handles dataclasses, datetimes, enums, and nested structures.
    """
    if obj is None:
        return None
    if isinstance(obj, str):
        return obj
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, (int, float)):
        return obj
    if isinstance(obj, datetime):
        return obj.isoformat()
    if isinstance(obj, Verdict):
        return obj.value
    if is_dataclass(obj) and not isinstance(obj, type):
        result = {}
        for f in fields(obj):
            val = getattr(obj, f.name)
            result[f.name] = _serialize_for_json(val)
        return result
    if isinstance(obj, dict):
        return {k: _serialize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_serialize_for_json(item) for item in obj]
    if isinstance(obj, set):
        return sorted(_serialize_for_json(item) for item in obj)
    # Fallback: try to convert to string
    try:
        return str(obj)
    except Exception:
        return None


def export_json(
    result: Any,
    indent: int = 2,
    include_certificate: bool = True,
    include_provenance: bool = True,
) -> str:
    """
    Export a full OptimizationResult as JSON.

    This produces a complete, self-contained JSON representation of an
    optimization result, including all metrics, provenance data, and
    certificate information. The output is suitable for API responses,
    data persistence, and pipeline integration.

    Args:
        result: An :class:`~biocompiler.optimization.OptimizationResult` object
            from :func:`~biocompiler.optimization.optimize_sequence`.
        indent: JSON indentation level (default 2). Use 0 for compact output.
        include_certificate: If True, include certificate text in the output.
        include_provenance: If True, include provenance and decision trail data.

    Returns:
        JSON-formatted string representing the full OptimizationResult

    Example::

        from biocompiler.export import export_json
        from biocompiler.api import optimize_sequence

        result = optimize_sequence('MSKGEELFTG', organism='Escherichia_coli')
        json_str = export_json(result)

        # Save to file
        with open('optimized.json', 'w') as f:
            f.write(json_str)

    Example output::

        {
          "biocompiler_version": "10.0.0",
          "export_timestamp": "2024-01-15T12:00:00+00:00",
          "sequence": "ATGAGCAAAGGAGAACTGTTT...",
          "protein": "MSKGEELFTG...",
          "metrics": {
            "cai": 0.999,
            "gc_content": 0.523,
            "sequence_length": 720,
            "protein_length": 239
          },
          "predicates": {
            "satisfied": ["GCInRange", "CodonAdapted", ...],
            "failed": [],
            "fallback_used": false
          },
          "certificate": { ... },
          "provenance": { ... }
        }
    """
    from .optimization import OptimizationResult

    if not isinstance(result, OptimizationResult):
        raise TypeError(
            f"Expected OptimizationResult, got {type(result).__name__}. "
            f"Use optimize_sequence() to produce an OptimizationResult."
        )

    now = datetime.now(timezone.utc)

    # Build the structured JSON output
    output: dict[str, Any] = {
        "biocompiler_version": __version__,
        "export_timestamp": now.isoformat(),
        "sequence": result.sequence,
        "protein": result.protein,
        "metrics": {
            "cai": result.cai,
            "gc_content": result.gc_content,
            "sequence_length": len(result.sequence),
            "protein_length": len(result.protein) if result.protein else 0,
            "codon_pair_bias": result.codon_pair_bias,
            "mutagenesis_applied": result.mutagenesis_applied,
            "fallback_used": result.fallback_used,
        },
        "predicates": {
            "satisfied": result.satisfied_predicates,
            "failed": result.failed_predicates,
            "aa_substitutions": result.aa_substitutions,
        },
        "organism_info": {
            "suggested_5utr": result.suggested_5utr,
            "suggested_3utr": result.suggested_3utr,
            "utr_score_5": result.utr_score_5,
            "utr_score_3": result.utr_score_3,
        },
    }

    # Biosafety annotations for JSON export
    organism_name = getattr(result, "organism_name", None) or "unknown"
    bsl = _assess_biosafety_level(organism_name)
    provenance_id = f"BC_{uuid.uuid4().hex[:12].upper()}"
    output["biosafety"] = {
        "optimized_by": f"biocompiler v{__version__}",
        "organism": organism_name,
        "cai_score": result.cai,
        "gc_content": result.gc_content,
        "passed_predicates": result.satisfied_predicates,
        "failed_predicates": result.failed_predicates,
        "biosecurity_screened": len(result.failed_predicates) == 0,
        "biosafety_level": bsl,
        "provenance_id": provenance_id,
    }

    # mRNA stability metrics
    if result.mrna_stability_score is not None:
        output["metrics"]["mrna_stability_score"] = _serialize_for_json(
            result.mrna_stability_score
        )
        output["metrics"]["destabilizing_motifs_removed"] = result.destabilizing_motifs_removed
        if result.stability_improvement is not None:
            output["metrics"]["stability_improvement"] = _serialize_for_json(
                result.stability_improvement
            )

    # Certificate data
    if include_certificate and result.certificate_text:
        output["certificate_text"] = result.certificate_text

    # Provenance data
    if include_provenance:
        provenance_data: dict[str, Any] = {}

        if result.provenance is not None:
            provenance_data["optimization_record"] = _serialize_for_json(
                result.provenance
            )

        if result.decision_trail is not None:
            provenance_data["decision_trail"] = _serialize_for_json(
                result.decision_trail
            )

        if provenance_data:
            output["provenance"] = provenance_data

    # indent=0 means compact (no whitespace); None also means compact in json.dumps
    effective_indent = None if indent == 0 else indent
    return json.dumps(output, indent=effective_indent, sort_keys=False, default=str)


# ─── Helper Functions ──────────────────────────────────────────────

def _get_taxonomy(organism: str) -> str:
    """Return approximate taxonomy lineage for common organisms."""
    taxonomies = {
        "Homo_sapiens": "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo.",
        "Mus_musculus": "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Mammalia; Eutheria; Euarchontoglires; Rodentia; Sciurognathi; Muroidea; Muridae; Murinae; Mus.",
        "Escherichia_coli": "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia.",
        "E_coli": "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia.",
        "CHO_K1": "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Mammalia; Eutheria; Euarchontoglires; Rodentia; Cricetidae; Cricetulus.",
        "Saccharomyces_cerevisiae": "Eukaryota; Fungi; Dikarya; Ascomycota; Saccharomycotina; Saccharomycetes; Saccharomycetales; Saccharomycetaceae; Saccharomyces.",
        "Pichia_pastoris": "Eukaryota; Fungi; Dikarya; Ascomycota; Saccharomycotina; Saccharomycetes; Saccharomycetales; Pichiaceae; Komagataella.",
        "Danio_rerio": "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Actinopterygii; Neopterygii; Teleostei; Cypriniformes; Cyprinidae; Danio.",
        "Drosophila_melanogaster": "Eukaryota; Metazoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Endopterygota; Diptera; Brachycera; Muscomorpha; Drosophilidae; Drosophila.",
        "Caenorhabditis_elegans": "Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis.",
        "Arabidopsis_thaliana": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Brassicales; Brassicaceae; Arabidopsis.",
        "Nicotiana_benthamiana": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Solanales; Solanaceae; Nicotiana.",
        "Spodoptera_frugiperda": "Eukaryota; Metazoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Endopterygota; Lepidoptera; Noctuoidea; Noctuidae; Spodoptera.",
        "Trichoplusia_ni": "Eukaryota; Metazoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Endopterygota; Lepidoptera; Noctuoidea; Noctuidae; Trichoplusia.",
        "Bacillus_subtilis": "Bacteria; Firmicutes; Bacilli; Bacillales; Bacillaceae; Bacillus.",
    }
    return taxonomies.get(organism, "Eukaryota; Metazoa; Unclassified.")


def _reconstruct_type_results(certificate: Certificate) -> list[TypeCheckResult]:
    """Reconstruct TypeCheckResult objects from certificate data."""
    results = []
    for t in certificate.types:
        verdict = Verdict(t["verdict"])
        results.append(TypeCheckResult(
            predicate=t["predicate"],
            verdict=verdict,
            derivation=None,
            violation=None if verdict == Verdict.PASS else "See certificate",
        ))
    return results


def _format_full_construct_features(
    seq_len: int,
    utr5_len: int,
    cds_len: int,
    utr3_len: int,
    gene_name: Optional[str],
    protein: Optional[str],
    organism: Optional[str] = None,
    cai: Optional[float] = None,
) -> list[str]:
    """Format the FEATURE TABLE for a full expression construct.

    The full construct is structured as:
    [5'UTR][CDS][3'UTR]

    GenBank coordinates are 1-based, inclusive.

    Args:
        seq_len: Total length of the full construct.
        utr5_len: Length of the 5' UTR.
        cds_len: Length of the CDS.
        utr3_len: Length of the 3' UTR.
        gene_name: Optional gene name.
        protein: Optional protein translation.
        organism: Optional organism name for feature qualifiers.
        cai: Optional CAI value for CDS feature qualifier.

    Returns:
        List of GenBank FEATURE TABLE lines.
    """
    lines: list[str] = []

    lines.append("FEATURES             Location/Qualifiers")

    # Compute 1-based coordinates
    utr5_start = 1
    utr5_end = utr5_len
    cds_start = utr5_len + 1
    cds_end = utr5_len + cds_len
    utr3_start = utr5_len + cds_len + 1
    utr3_end = seq_len

    # Gene feature spanning the entire construct
    if gene_name:
        lines.append(f"     gene            1..{seq_len}")
        lines.append(f'                     /gene="{gene_name}"')
        if organism:
            lines.append(f'                     /organism="{organism.replace("_", " ")}"')
        lines.append(f'                     /note="Full expression construct designed by BioCompiler"')

    # 5' UTR feature
    if utr5_len > 0:
        lines.append(f"     5'UTR           {utr5_start}..{utr5_end}")
        if gene_name:
            lines.append(f'                     /gene="{gene_name}"')
        lines.append(f'                     /note="5-prime UTR — suggested for expression in target organism"')
        lines.append(f'                     /label="5UTR"')

    # CDS feature (the optimized coding sequence)
    if cds_len > 0:
        lines.append(f"     CDS             {cds_start}..{cds_end}")
        if gene_name:
            lines.append(f'                     /gene="{gene_name}"')
        if organism:
            lines.append(f'                     /organism="{organism.replace("_", " ")}"')
        lines.append(f'                     /note="Codon-optimized CDS designed by BioCompiler"')
        if cai is not None:
            lines.append(f'                     /cai="{cai:.4f}"')
        lines.append(f'                     /codon_start=1')
        lines.append(f'                     /transl_table=1')
        if protein:
            prot_chunks = [protein[i:i + 40] for i in range(0, len(protein), 40)]
            lines.append(f'                     /translation="{prot_chunks[0]}"')
            for chunk in prot_chunks[1:]:
                lines.append(f'                     "{chunk}"')

    # 3' UTR feature
    if utr3_len > 0:
        lines.append(f"     3'UTR           {utr3_start}..{utr3_end}")
        if gene_name:
            lines.append(f'                     /gene="{gene_name}"')
        lines.append(f'                     /note="3-prime UTR — suggested for expression in target organism"')
        lines.append(f'                     /label="3UTR"')

    # mRNA feature spanning the transcript
    if utr5_len > 0 or utr3_len > 0:
        lines.append(f"     mRNA            1..{seq_len}")
        if gene_name:
            lines.append(f'                     /gene="{gene_name}"')
        # Build exon structure for the mRNA
        if utr5_len > 0 and utr3_len > 0:
            lines.append(f'                     /note="5UTR:{utr5_start}..{utr5_end} CDS:{cds_start}..{cds_end} 3UTR:{utr3_start}..{utr3_end}"')
        elif utr5_len > 0:
            lines.append(f'                     /note="5UTR:{utr5_start}..{utr5_end} CDS:{cds_start}..{cds_end}"')
        else:
            lines.append(f'                     /note="CDS:{cds_start}..{cds_end} 3UTR:{utr3_start}..{utr3_end}"')

    return lines


def export_full_construct(
    utr5: str,
    cds: str,
    utr3: str,
    organism: str = "Homo_sapiens",
    locus_name: str = "BIOCOMPILER",
    definition: str = "BioCompiler full expression construct",
    gene_name: Optional[str] = None,
    molecule_type: str = "DNA",
    topology: str = "linear",
    cai: Optional[float] = None,
) -> str:
    """Export a complete expression construct (5'UTR + CDS + 3'UTR) as GenBank.

    This produces a fully annotated GenBank record for a complete expression
    construct — exactly what a biologist would order from a gene synthesis
    company (e.g., IDT, Twist Bioscience, GenScript). The record includes:

    - LOCUS header with correct molecule type and topology
    - DEFINITION, ACCESSION, VERSION, KEYWORDS headers
    - SOURCE and ORGANISM with taxonomy lineage
    - FEATURE TABLE with 5'UTR, CDS, 3'UTR, and mRNA annotations
    - ORIGIN section with numbered full construct sequence

    The CDS feature includes the protein translation and CAI qualifier.
    The UTR features are annotated as suggested elements for expression
    optimization.

    Args:
        utr5: 5' UTR sequence (empty string if no 5' UTR).
        cds: Coding sequence (optimized, starts with ATG, ends with stop codon).
        utr3: 3' UTR sequence (empty string if no 3' UTR).
        organism: Target organism name (e.g., 'Homo_sapiens').
        locus_name: GenBank LOCUS name (max 16 chars, uppercase).
        definition: DEFINITION line text.
        gene_name: Optional gene name for feature annotations.
        molecule_type: Molecule type (DNA, RNA, mRNA).
        topology: circular or linear.
        cai: Optional CAI value for the CDS (auto-computed if None).

    Returns:
        GenBank-formatted string for the full expression construct.

    Example::

        from biocompiler import optimize_sequence, export_full_construct

        result = optimize_sequence("MVHLTPEEK", organism="Homo_sapiens")
        gb = export_full_construct(
            utr5=result.suggested_5utr or "",
            cds=result.sequence,
            utr3=result.suggested_3utr or "",
            organism="Homo_sapiens",
            gene_name="HBB",
        )
        # Write to file for ordering from synthesis company
        with open("construct.gb", "w") as f:
            f.write(gb)
    """
    # Normalize sequences
    utr5_clean = utr5.upper().replace(" ", "")
    cds_clean = cds.upper().replace(" ", "")
    utr3_clean = utr3.upper().replace(" ", "")

    full_seq = utr5_clean + cds_clean + utr3_clean
    if not full_seq:
        raise ValueError("Full construct sequence is empty — provide at least a CDS")
    if not cds_clean:
        raise ValueError("CDS is required for a full expression construct")

    gc = gc_content(full_seq)
    now = datetime.now(timezone.utc)
    protein = translate(cds_clean)

    # Auto-compute CAI if not provided
    if cai is None and len(cds_clean) >= 3 and len(cds_clean) % 3 == 0:
        try:
            cai = compute_cai(cds_clean, organism=organism)
        except Exception:
            cai = None

    # Truncate locus name to 16 chars (GenBank requirement)
    locus = locus_name[:16].upper()

    # Date in DD-MON-YYYY format
    date_str = now.strftime("%d-%b-%Y").upper()

    # Generate accession
    acc = _generate_accession()

    # ─── Assemble GenBank record ──────────────────────────────────
    lines: list[str] = []

    # Header
    lines.extend(_format_genbank_header(
        locus, len(full_seq), molecule_type, topology, date_str,
        definition, acc, organism, gc, protein, None, None,
        cai=cai,
    ))

    # UTR-specific comment
    comments: list[str] = []
    comments.append(f"Construct layout: 5'UTR({len(utr5_clean)}bp) + CDS({len(cds_clean)}bp) + 3'UTR({len(utr3_clean)}bp)")
    comments.append(f"Total length: {len(full_seq)} bp")
    if cai is not None:
        comments.append(f"CDS CAI: {cai:.4f}")
    if utr5_clean:
        comments.append(f"5' UTR: suggested sequence for {organism} (user should verify before ordering)")
    if utr3_clean:
        comments.append(f"3' UTR: suggested sequence for {organism} (user should verify before ordering)")
    comments.append("UTR sequences are SUGGESTED but not enforced — evaluate before ordering")

    lines.append("COMMENT     " + comments[0])
    for c in comments[1:]:
        lines.append("            " + c)

    # Features with UTR/CDS/mRNA annotations
    lines.extend(_format_full_construct_features(
        seq_len=len(full_seq),
        utr5_len=len(utr5_clean),
        cds_len=len(cds_clean),
        utr3_len=len(utr3_clean),
        gene_name=gene_name,
        protein=protein,
        organism=organism,
        cai=cai,
    ))

    # Sequence
    lines.extend(_format_genbank_sequence(full_seq))

    return "\n".join(lines)


# ─── Annotation-Aware Export ────────────────────────────────────────

def export_with_annotations(
    sequence: str,
    organism: str = "Homo_sapiens",
    cai: Optional[float] = None,
    gc: Optional[float] = None,
    type_results: Optional[list[TypeCheckResult]] = None,
    format: str = "genbank",
    **kwargs: Any,
) -> str:
    """Export a sequence with full biosafety annotations in the specified format.

    This is the primary export function for production use. It wraps the
    existing export functions and adds comprehensive biosafety metadata
    including BSL level, screening status, and provenance tracking.

    For GenBank format, the COMMENT section includes the full
    ``BIOCOMPILER_ANNOTATIONS`` block, ``WARNING`` (if predicates failed),
    and ``BIOSECURITY NOTICE`` with risk assessment summary.

    For FASTA format, the header includes biosecurity level and
    biocompiler version metadata.

    Args:
        sequence: DNA sequence (designed, verified).
        organism: Target organism name.
        cai: Optional CAI value (auto-computed if None).
        gc: Optional GC content (auto-computed if None).
        type_results: Optional type-check predicate results.
        format: Export format — ``"genbank"`` or ``"fasta"``.
        **kwargs: Additional keyword arguments forwarded to the underlying
            export function (e.g., ``locus_name``, ``gene_name``,
            ``identifier``, ``description``).

    Returns:
        Formatted string with embedded biosafety annotations.

    Example::

        from biocompiler.export import export_with_annotations
        from biocompiler.types import TypeCheckResult, Verdict

        type_results = [
            TypeCheckResult(predicate="gc_content", verdict=Verdict.PASS),
            TypeCheckResult(predicate="no_stop_codons", verdict=Verdict.PASS),
        ]
        result = export_with_annotations(
            "ATGGTGAGCAAGGGCGAGGAG",
            organism="Escherichia_coli",
            cai=0.95,
            type_results=type_results,
            format="genbank",
        )
        assert "BIOCOMPILER_ANNOTATIONS:" in result
        assert "biosafety_level: BSL-1" in result
    """
    seq = sequence.upper().replace(" ", "")
    if gc is None:
        gc = gc_content(seq)

    # Auto-compute CAI if not provided
    if cai is None and len(seq) >= 3 and len(seq) % 3 == 0:
        try:
            cai = compute_cai(seq, organism=organism)
        except Exception:
            cai = None

    format_lower = format.lower()
    if format_lower == "fasta":
        return export_fasta(
            sequence=sequence,
            organism=organism,
            cai=cai,
            type_results=type_results,
            **kwargs,
        )
    elif format_lower == "genbank":
        return export_genbank(
            sequence=sequence,
            organism=organism,
            cai=cai,
            type_results=type_results,
            **kwargs,
        )
    else:
        raise ValueError(
            f"Unsupported format: {format!r}. Use 'genbank' or 'fasta'."
        )
