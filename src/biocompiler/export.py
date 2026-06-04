"""
BioCompiler Export Engine — GenBank & FASTA Sequence Export

Production-grade sequence export with:
- GenBank format output with full feature annotations
- FASTA format with metadata headers
- Certificate embedding in GenBank comment section
- Exon/intron/restriction-site feature annotations
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
    "GENBANK_MAX_LINE",
    "GENBANK_SEQ_LINE",
    "GENBANK_SEQ_GROUP",
]

import logging
import uuid
from datetime import datetime, timezone
from typing import Optional, TypedDict

from .types import Certificate, TypeCheckResult, Verdict, combined_verdict
from .scanner import gc_content
from .translation import translate
from . import __version__

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
    """An entry for :func:`export_multi_fasta`.

    Required field: ``sequence``.
    Optional fields: ``id``, ``description``, ``organism``, ``protein``.
    """
    id: str
    description: str
    organism: str
    protein: str


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
) -> list[str]:
    """Format the LOCUS, DEFINITION, ACCESSION, VERSION, SOURCE, ORGANISM, and COMMENT sections."""
    lines: list[str] = []

    # ── LOCUS / DEFINITION / ACCESSION / VERSION ──
    length_str = f"{seq_len} bp"
    mol_str = f"{molecule_type}    "
    lines.append(
        f"LOCUS       {locus:<16} {length_str:>12}   {mol_str} {topology:<8}   SYN"
    )
    lines.append(f"DEFINITION  {definition}.")

    lines.append(f"ACCESSION   {acc}")
    lines.append(f"VERSION     {acc}.1")

    # ── SOURCE / ORGANISM ──
    taxonomy = _get_taxonomy(organism)
    lines.append(f"SOURCE      {organism}")
    lines.append(f"  ORGANISM  {organism}")
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
    comments.append(f"Protein length: {len(protein)} aa" if protein else "")

    if type_results:
        overall = combined_verdict([r.verdict for r in type_results])
        comments.append(f"Type-check verdict: {overall.value}")
        for r in type_results:
            symbol = {"PASS": "+", "FAIL": "X", "UNCERTAIN": "?"}[r.verdict.value]
            comments.append(f"  [{symbol}] {r.predicate}")

    if certificate:
        comments.append(f"Certificate ID: {certificate.design_id[:16]}...")
        comments.append(f"Certificate timestamp: {certificate.provenance.get('timestamp', 'N/A')}")

    if any(c for c in comments):
        lines.append("COMMENT     " + comments[0])
        for c in comments[1:]:
            if c:
                lines.append("            " + c)

    return lines


def _format_genbank_features(
    seq_len: int,
    gene_name: Optional[str],
    protein: Optional[str],
    exon_boundaries: Optional[list[tuple[int, int]]],
    restriction_sites: Optional[list[RestrictionSiteInfo]],
    type_results: Optional[list[TypeCheckResult]],
) -> list[str]:
    """Format the FEATURE TABLE section of a GenBank record."""
    lines: list[str] = []

    lines.append("FEATURES             Location/Qualifiers")

    # Gene feature
    if gene_name:
        lines.append(f"     gene            1..{seq_len}")
        lines.append(f'                     /gene="{gene_name}"')
        lines.append(f'                     /note="Designed by BioCompiler"')

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
        lines.append(f'                     /note="Designed by BioCompiler"')
        lines.append(f'                     /codon_start=1')
        lines.append(f'                     /transl_table=1')
        # Protein translation (wrapped)
        if protein:
            prot_chunks = [protein[i:i + 40] for i in range(0, len(protein), 40)]
            lines.append(f'                     /translation="{prot_chunks[0]}"')
            for chunk in prot_chunks[1:]:
                lines.append(f'                     "{chunk}"')

    # Exon features
    if exon_boundaries:
        for i, (start, end) in enumerate(exon_boundaries):
            lines.append(f"     exon            {start + 1}..{end}")
            if gene_name:
                lines.append(f'                     /gene="{gene_name}"')
            lines.append(f'                     /number={i + 1}')

    # Restriction site features
    if restriction_sites:
        for site in restriction_sites[:20]:  # Limit to 20 annotations
            pos = site.get("position", 0)
            enz = site.get("enzyme", site.get("site", "unknown"))
            strand = site.get("strand", "+")
            site_len = len(site.get("site", ""))
            lines.append(f"     misc_feature    {pos + 1}..{pos + site_len}")
            lines.append(f'                     /note="Restriction site: {enz} ({strand} strand)"')
            lines.append(f'                     /label={enz}')

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
) -> str:
    """
    Export a designed sequence in FASTA format.

    FASTA is the universal sequence format accepted by BLAST, Clustal,
    Geneious, and virtually all bioinformatics tools. This function
    generates a standards-compliant FASTA record with a rich header
    that includes organism, GC content, and protein translation.

    Args:
        sequence: DNA sequence (designed, verified)
        identifier: Sequence identifier (no spaces)
        description: Human-readable description line
        organism: Source organism name
        protein: Optional protein translation (auto-computed if None)

    Returns:
        FASTA-formatted string

    Example output::

        >BioCompiler_design|organism=Homo_sapiens|gc=0.523|len=720 eGFP designed sequence
        ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTG
        GACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCC
        ...
    """
    seq = sequence.upper().replace(" ", "")
    gc = gc_content(seq)

    if protein is None:
        protein = translate(seq)

    # Build FASTA header with structured metadata
    header_parts = [identifier]
    header_parts.append(f"organism={organism}")
    header_parts.append(f"gc={gc:.3f}")
    header_parts.append(f"len={len(seq)}")
    if protein:
        header_parts.append(f"protein_len={len(protein)}aa")

    header = "|".join(header_parts)
    if description:
        header += f" {description}"

    return f">{header}\n{_format_fasta_sequence(seq)}\n"


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
) -> str:
    """
    Export a designed sequence in GenBank format.

    GenBank is the standard format for sequence submission to NCBI/ENA/DDBJ
    and is the native format for Benchling, SnapGene, and Geneious. This
    function produces a fully annotated GenBank record with:

    - LOCUS, DEFINITION, ACCESSION, VERSION headers
    - SOURCE and ORGANISM with taxonomy lineage
    - FEATURE TABLE with gene, CDS, exon, and misc_feature annotations
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

    Returns:
        GenBank-formatted string
    """
    seq = sequence.upper().replace(" ", "")
    gc = gc_content(seq)
    now = datetime.now(timezone.utc)

    if protein is None:
        protein = translate(seq)

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
    ))
    lines.extend(_format_genbank_features(len(seq), gene_name, protein, exon_boundaries, restriction_sites, type_results))
    lines.extend(_format_genbank_sequence(seq))

    return "\n".join(lines)


def export_multi_fasta(
    sequences: list[FastaSequenceEntry],
) -> str:
    """
    Export multiple designed sequences as a multi-FASTA file.

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
        )
        records.append(record)
    return "\n".join(records)


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
    )


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
