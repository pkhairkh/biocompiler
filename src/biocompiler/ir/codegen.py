"""
BioCompiler IR Codegen — Backend Output Passes
==============================================

Emits standard biological file formats from IR objects:

  - IR_L0 → GenBank format (genomic DNA with annotations)
  - IR_L2 → FASTA format (mature mRNA)
  - IR_L3 → FASTA format (protein)
  - IR_L0 → SBOL3 (synthetic biology standard) — Phase 2

These are the compiler's "backend": once the IR pipeline has produced a
mature gene design (L0) or a translated protein (L3), these functions
serialise it into a string that downstream tools — Benchling, SnapGene,
BLAST, NCBI — can ingest directly.

Design rules
------------
* Pure functions: same IR → same string, no I/O, no side effects.
* The output is a ``str`` (not bytes) — callers can write it to a file
  with their own encoding policy.
* Coordinates in IR are **0-indexed half-open** ``[start, end)``;
  GenBank uses **1-indexed inclusive** ``start..end``.  Conversion is
  ``gb_start = ir_start + 1``, ``gb_end = ir_end`` (since ``ir_end`` is
  exclusive, it equals the 1-based inclusive end).
* FASTA wraps at 60 characters per line, per the de-facto standard.

Usage
-----
::

    from biocompiler.ir.codegen import to_genbank, to_fasta
    from biocompiler.ir import IR_L0_GenomicDNA, IRLevel, compile_gene

    gene = IR_L0_GenomicDNA(
        sequence="ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAG"
                  "GTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAA",
        regions=[], organism="Homo_sapiens", gene_name="HBB",
    )
    protein = compile_gene(gene, IRLevel.L3)

    print(to_genbank(gene))
    print(to_fasta_protein(protein))
"""

from __future__ import annotations

import re
import uuid
from datetime import datetime, timezone
from typing import Optional

from .types import (
    IR_L0_GenomicDNA,
    IR_L2_MatureMRNA,
    IR_L3_Polypeptide,
    GeneRegion,
    IRLevel,
)


__all__ = [
    "to_genbank",
    "to_fasta_dna",
    "to_fasta_rna",
    "to_fasta_protein",
    "to_fasta",
    "to_sbol3",
    "parse_genbank_sequence",
    "parse_genbank_features",
]


# ────────────────────────────────────────────────────────────────────
# Format constants
# ────────────────────────────────────────────────────────────────────

# FASTA: 60 residues per line — the de-facto standard (BLAST, NCBI).
FASTA_LINE_WIDTH = 60

# GenBank: 60 bases per sequence line, grouped in blocks of 10.
GENBANK_SEQ_LINE = 60
GENBANK_SEQ_GROUP = 10

# Map IR region_type strings → GenBank feature keys.
# Anything not in this map is emitted as ``misc_feature`` with a
# ``/label=`` qualifier carrying the original region_type, so the
# information is never lost.
_REGION_TYPE_TO_GB_KEY: dict[str, str] = {
    "exon": "exon",
    "intron": "intron",
    "5_utr": "5'UTR",
    "3_utr": "3'UTR",
    "promoter": "promoter",
    "terminator": "terminator",
    "cds": "CDS",
}


# ────────────────────────────────────────────────────────────────────
# Internal helpers
# ────────────────────────────────────────────────────────────────────

def _wrap_sequence(seq: str, width: int = FASTA_LINE_WIDTH) -> str:
    """Wrap a sequence string into fixed-width lines."""
    seq = seq.upper()
    if not seq:
        return ""
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def _gb_location(start: int, end: int) -> str:
    """Convert IR 0-based ``[start, end)`` → GenBank 1-based ``start..end``."""
    # IR end is exclusive; GenBank end is inclusive, so they coincide.
    return f"{start + 1}..{end}"


def _gb_date() -> str:
    """GenBank-style date: ``DD-MON-YYYY`` (e.g. ``01-JAN-2025``)."""
    return datetime.now(timezone.utc).strftime("%d-%b-%Y").upper()


def _pseudo_accession() -> str:
    """Generate a pseudo-accession ``BC`` + 8 hex chars (BC = BioCompiler)."""
    return f"BC{uuid.uuid4().hex[:8].upper()}"


def _organism_display(organism: str) -> str:
    """Pretty-print an organism key (``Homo_sapiens`` → ``Homo sapiens``)."""
    return organism.replace("_", " ") if organism else "unknown organism"


# ────────────────────────────────────────────────────────────────────
# GenBank codegen
# ────────────────────────────────────────────────────────────────────

def _build_genbank_features(ir_l0: IR_L0_GenomicDNA) -> list[str]:
    """Build the FEATURES table lines from IR-L0 regions.

    Emits (in order):

    1. A ``gene`` feature spanning the whole sequence.
    2. A ``CDS`` feature — ``join()`` of every ``exon``/``cds`` region
       (or the whole sequence if there are no region annotations, which
       is the prokaryotic / single-exon case).
    3. One ``exon`` feature per exon region, with a ``/number=`` qualifier.
    4. One feature per remaining region type (5'UTR, 3'UTR, promoter,
       terminator, intron, misc) with appropriate GenBank keys.
    """
    seq_len = len(ir_l0.sequence)
    gene_name = ir_l0.gene_name or "gene"
    lines: list[str] = []
    regions_sorted = sorted(ir_l0.regions, key=lambda r: r.start)

    lines.append("FEATURES             Location/Qualifiers")

    # ── gene feature: spans the whole sequence ──
    lines.append(f"     gene            {_gb_location(0, seq_len)}")
    lines.append(f'                     /gene="{gene_name}"')
    lines.append(f'                     /organism="{_organism_display(ir_l0.organism)}"')
    lines.append(f'                     /locus_tag="{gene_name.upper()}"')

    # ── CDS feature ──
    # Use explicit exon/cds regions if any; otherwise the whole sequence
    # (prokaryotic case — no introns).
    coding_regions = [
        r for r in regions_sorted if r.region_type in ("exon", "cds")
    ]
    if coding_regions:
        parts = [_gb_location(r.start, r.end) for r in coding_regions]
        cds_loc = f"join({','.join(parts)})" if len(parts) > 1 else parts[0]
    else:
        cds_loc = _gb_location(0, seq_len)

    lines.append(f"     CDS             {cds_loc}")
    lines.append(f'                     /gene="{gene_name}"')
    lines.append(f'                     /organism="{_organism_display(ir_l0.organism)}"')
    lines.append('                     /codon_start=1')
    lines.append('                     /transl_table=1')
    lines.append('                     /note="BioCompiler IR-L0 genomic DNA"')

    # ── exon features ──
    exon_idx = 0
    for r in regions_sorted:
        if r.region_type == "exon":
            exon_idx += 1
            lines.append(f"     exon            {_gb_location(r.start, r.end)}")
            lines.append(f'                     /gene="{gene_name}"')
            lines.append(f'                     /number={exon_idx}')

    # ── other region types (5'UTR, 3'UTR, promoter, intron, terminator, …) ──
    for r in regions_sorted:
        rtype = r.region_type
        if rtype in ("exon", "cds"):
            continue  # already covered by the CDS / exon features above
        gb_key = _REGION_TYPE_TO_GB_KEY.get(rtype, "misc_feature")
        lines.append(f"     {gb_key:<16}{_gb_location(r.start, r.end)}")
        lines.append(f'                     /gene="{gene_name}"')
        if rtype not in _REGION_TYPE_TO_GB_KEY:
            lines.append(f'                     /label="{rtype}"')
        else:
            lines.append(f'                     /note="{rtype}"')

    return lines


def _build_genbank_origin(seq: str) -> list[str]:
    """Build the ORIGIN section (numbered sequence) + terminator."""
    seq = seq.upper()
    lines = ["ORIGIN"]
    for i in range(0, len(seq), GENBANK_SEQ_LINE):
        chunk = seq[i:i + GENBANK_SEQ_LINE]
        groups = [
            chunk[j:j + GENBANK_SEQ_GROUP]
            for j in range(0, len(chunk), GENBANK_SEQ_GROUP)
        ]
        line_num = i + 1
        lines.append(f"{line_num:>9} {' '.join(groups)}")
    lines.append("//")
    return lines


def to_genbank(ir_l0: IR_L0_GenomicDNA) -> str:
    """Emit a GenBank-format string from an IR-L0 genomic DNA object.

    The output includes the standard GenBank record sections:

    * ``LOCUS``      — name, length, molecule type (DNA), topology, date
    * ``DEFINITION`` — gene + organism description
    * ``ACCESSION``  — pseudo-accession ``BC<hex>``
    * ``VERSION``    — accession.1
    * ``KEYWORDS``   — BioCompiler marker
    * ``SOURCE`` / ``ORGANISM`` — organism display name
    * ``FEATURES``   — ``gene``, ``CDS``, ``exon`` features (and 5'UTR /
      3'UTR / promoter / intron / terminator when present), with
      qualifiers (``/gene``, ``/organism``, ``/codon_start``,
      ``/transl_table``, ``/number``).
    * ``ORIGIN``     — numbered sequence (60 bases/line, groups of 10)
    * ``//``         — record terminator

    Coordinates are converted from IR's 0-based half-open ``[start, end)``
    to GenBank's 1-based inclusive ``start..end``.

    Args:
        ir_l0: the IR-L0 genomic DNA object to serialise.

    Returns:
        A GenBank-format string ending with ``//``.

    Example::

        >>> from biocompiler.ir import IR_L0_GenomicDNA
        >>> from biocompiler.ir.codegen import to_genbank
        >>> ir = IR_L0_GenomicDNA(
        ...     sequence="ATGGCTTAA", regions=[],
        ...     organism="e_coli", gene_name="test",
        ... )
        >>> gb = to_genbank(ir)
        >>> gb.startswith("LOCUS")
        True
        >>> gb.endswith("//\\n")
        True
    """
    seq = ir_l0.sequence.upper()
    seq_len = len(seq)
    gene_name = ir_l0.gene_name or "gene"
    organism_display = _organism_display(ir_l0.organism)
    accession = _pseudo_accession()
    date_str = _gb_date()

    # ── LOCUS ──
    # Field widths follow the GenBank standard:
    #   LOCUS  <name 16>  <len> bp   <mol_type 7>  <topology 8>  <div 3>  <date>
    locus_line = (
        f"LOCUS       {gene_name:<16} {seq_len:>11} bp    DNA    "
        f"linear   SYN {date_str}"
    )

    # ── DEFINITION / ACCESSION / VERSION / KEYWORDS ──
    definition = f"{gene_name} gene from {organism_display} (BioCompiler IR-L0)."
    # Wrap DEFINITION at 68 chars after the 12-char keyword prefix.
    wrapped_def: list[str] = []
    line = definition
    while line:
        wrapped_def.append(line[:68])
        line = line[68:]

    header_lines: list[str] = [
        locus_line,
        f"DEFINITION  {wrapped_def[0]}",
        *[f"            {ln}" for ln in wrapped_def[1:]],
        f"ACCESSION   {accession}",
        f"VERSION     {accession}.1",
        "KEYWORDS    BioCompiler; IR-L0; genomic DNA.",
        f"SOURCE      {organism_display}",
        f"  ORGANISM  {organism_display}",
        "            Eukaryota; BioCompiler synthetic gene.",
    ]

    features_lines = _build_genbank_features(ir_l0)
    origin_lines = _build_genbank_origin(seq)

    return "\n".join(header_lines + features_lines + origin_lines) + "\n"


# ────────────────────────────────────────────────────────────────────
# FASTA codegen
# ────────────────────────────────────────────────────────────────────

def _fasta_header(
    gene_name: Optional[str],
    organism: str,
    level: str,
    seq_len: int,
) -> str:
    """Build a FASTA header ``>gene|organism|level``.

    Falls back to ``BC_<hex>`` for the gene name when none is set, and
    ``unknown`` for an empty organism, so the header is never malformed.
    """
    name = gene_name or f"BC_{uuid.uuid4().hex[:6].upper()}"
    org = organism or "unknown"
    return f">{name}|{org}|{level}|len={seq_len}"


def to_fasta_dna(ir_l0: IR_L0_GenomicDNA) -> str:
    """Emit a FASTA-format string from an IR-L0 genomic DNA object.

    Header: ``>gene_name|organism|L0|len=<N>``
    Sequence: uppercase DNA, wrapped at 60 chars/line.

    Args:
        ir_l0: the IR-L0 genomic DNA object.

    Returns:
        A FASTA-format string.
    """
    seq = ir_l0.sequence.upper()
    header = _fasta_header(ir_l0.gene_name, ir_l0.organism, "L0", len(seq))
    body = _wrap_sequence(seq, FASTA_LINE_WIDTH)
    return f"{header}\n{body}\n" if body else f"{header}\n"


def to_fasta_rna(ir_l2: IR_L2_MatureMRNA) -> str:
    """Emit a FASTA-format string from an IR-L2 mature mRNA object.

    Header: ``>gene_name|organism|L2|len=<N>``
    Sequence: the **full mature mRNA** (5'UTR + CDS + 3'UTR), uppercase,
    wrapped at 60 chars/line.

    Args:
        ir_l2: the IR-L2 mature mRNA object.

    Returns:
        A FASTA-format string.
    """
    seq = ir_l2.sequence.upper()
    header = _fasta_header(ir_l2.gene_name, ir_l2.organism, "L2", len(seq))
    body = _wrap_sequence(seq, FASTA_LINE_WIDTH)
    return f"{header}\n{body}\n" if body else f"{header}\n"


def to_fasta_protein(ir_l3: IR_L3_Polypeptide) -> str:
    """Emit a FASTA-format string from an IR-L3 polypeptide object.

    Header: ``>gene_name|organism|L3|len=<N>``
    Sequence: single-letter amino acid codes (including ``*`` for stop),
    wrapped at 60 chars/line.

    Args:
        ir_l3: the IR-L3 polypeptide object.

    Returns:
        A FASTA-format string.

    Example::

        >>> from biocompiler.ir import IR_L3_Polypeptide
        >>> from biocompiler.ir.codegen import to_fasta_protein
        >>> p = IR_L3_Polypeptide(sequence="MAK*", organism="e_coli",
        ...                       gene_name="test")
        >>> to_fasta_protein(p)
        '>test|e_coli|L3|len=4\\nMAK*\\n'
    """
    seq = ir_l3.sequence.upper()
    header = _fasta_header(ir_l3.gene_name, ir_l3.organism, "L3", len(seq))
    body = _wrap_sequence(seq, FASTA_LINE_WIDTH)
    return f"{header}\n{body}\n" if body else f"{header}\n"


def to_fasta(ir) -> str:
    """Dispatch FASTA codegen based on the IR object's level.

    Convenience wrapper that picks :func:`to_fasta_dna`,
    :func:`to_fasta_rna`, or :func:`to_fasta_protein` based on the IR
    object's ``level`` property.

    Args:
        ir: an IR object with a ``level`` property
            (:class:`~biocompiler.ir.types.IR_L0_GenomicDNA`,
            :class:`~biocompiler.ir.types.IR_L2_MatureMRNA`, or
            :class:`~biocompiler.ir.types.IR_L3_Polypeptide`).

    Returns:
        A FASTA-format string.

    Raises:
        TypeError: if the IR object's level is not L0, L2, or L3.
    """
    level = getattr(ir, "level", None)
    if level == IRLevel.L0:
        return to_fasta_dna(ir)
    if level == IRLevel.L2:
        return to_fasta_rna(ir)
    if level == IRLevel.L3:
        return to_fasta_protein(ir)
    raise TypeError(
        f"to_fasta() supports IR levels L0, L2, L3 — got {level!r}"
    )


# ────────────────────────────────────────────────────────────────────
# SBOL3 codegen (RDF/Turtle)
# ────────────────────────────────────────────────────────────────────

# Map IR region_type strings → Sequence Ontology (SO) IDs used as
# ``sbol:role`` values on SequenceFeature objects.
_REGION_TYPE_TO_SO: dict[str, str] = {
    "exon": "0000147",        # SO:0000147 exon
    "intron": "0000188",      # SO:0000188 intron
    "5_utr": "0000204",       # SO:0000204 five_prime_UTR
    "3_utr": "0000205",       # SO:0000205 three_prime_UTR
    "promoter": "0000167",    # SO:0000167 promoter
    "terminator": "0000141",  # SO:0000141 terminator
    "cds": "0000316",         # SO:0000316 CDS
}

# SO:0000704 = gene (used for the whole-sequence gene feature).
_SO_GENE = "0000704"
# SO:0001411 = region (default for unrecognised region types).
_SO_MISC_REGION = "0001411"

# SBOL3 namespaces (as Turtle @prefix URIs).
_SBOL3_NS = "http://sbols.org/v3#"
_RDF_NS = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
_RDFS_NS = "http://www.w3.org/2000/01/rdf-schema#"
_BIOPAX_NS = "http://www.biopax.org/release/biopax-level3.owl#"
_SO_NS = "http://sequenceontology.org/resource/SO:"
# DNA encoding URI (IUPAC DNA alphabet).
_SBOL3_DNA_ENCODING = (
    "http://www.cheminfo.org/cheminfo/emolecules/encoding#DNA"
)

# Base URI for BioCompiler-emitted SBOL3 identities. Per the SBOL3 spec,
# identity URIs must be globally unique; we use a stable scheme based
# on ``gene_name`` + version so the same IR produces the same URI.
_SBOL3_BASE_URI = "http://biocompiler.org/gene"
_SBOL3_VERSION = "1"


def _sbol3_safe_display_id(gene_name: Optional[str]) -> str:
    """Sanitise a gene name into a valid SBOL3 ``displayId``.

    SBOL3 displayIds must match ``[a-zA-Z_][a-zA-Z0-9_]*``.  We replace
    any character outside that class with ``_`` and prefix ``gene_`` if
    the result would start with a digit.  Falls back to ``"gene"``
    when ``gene_name`` is ``None`` or empty.
    """
    raw = (gene_name or "").strip() or "gene"
    safe = re.sub(r"[^A-Za-z0-9_]", "_", raw)
    if not safe:
        safe = "gene"
    if not (safe[0].isalpha() or safe[0] == "_"):
        safe = "gene_" + safe
    return safe


def _escape_turtle_string(s: str) -> str:
    """Escape a string for use inside a Turtle ``"..."`` literal."""
    # Escape backslash and double-quote, then newlines/tabs.
    return (
        s.replace("\\", "\\\\")
         .replace('"', '\\"')
         .replace("\n", "\\n")
         .replace("\t", "\\t")
    )


def _build_sbol3_features(
    ir_l0: IR_L0_GenomicDNA,
    component_uri: str,
    sequence_uri: str,
) -> tuple[list[str], list[str]]:
    """Build the SequenceFeature + Range blocks for an IR-L0 object.

    Returns a tuple ``(feature_uris, feature_block_lines)``:

    * ``feature_uris`` — list of URIs to be referenced from the
      Component's ``sbol:hasFeature`` property.
    * ``feature_block_lines`` — Turtle lines describing each
      SequenceFeature and its Range location, to be appended after the
      Component and Sequence blocks.

    The following features are emitted:

    1. A **gene** feature spanning the whole sequence (role SO:0000704).
    2. One **CDS** feature per ``exon``/``cds`` region (role SO:0000316),
       *or* — in the prokaryotic case (no exon/cds regions) — a single
       CDS feature spanning the whole sequence.
    3. One feature per remaining region type (5'UTR, 3'UTR, promoter,
       intron, terminator, misc), with the appropriate SO role.
    """
    seq_len = len(ir_l0.sequence)
    regions_sorted = sorted(ir_l0.regions, key=lambda r: r.start)

    feature_uris: list[str] = []
    feature_lines: list[str] = []

    def _emit_feature(display_id: str, so_id: str, start_1based: int, end_1based: int,
                      slug: str) -> None:
        """Append a SequenceFeature + its Range location to the output."""
        f_uri = f"{component_uri}/feature/{slug}"
        l_uri = f"{f_uri}/location"
        feature_uris.append(f_uri)
        feature_lines.extend([
            f"<{f_uri}> a sbol:SequenceFeature ;",
            f'    sbol:displayId "{_escape_turtle_string(display_id)}" ;',
            f"    sbol:role so:{so_id} ;",
            f"    sbol:hasLocation <{l_uri}> .",
            f"<{l_uri}> a sbol:Range ;",
            f'    sbol:displayId "{_escape_turtle_string(display_id)}_location" ;',
            f"    sbol:sequence <{sequence_uri}> ;",
            f"    sbol:start {start_1based} ;",
            f"    sbol:end {end_1based} .",
            "",  # blank line between feature blocks
        ])

    # ── 1. Gene feature spanning the whole sequence ──
    _emit_feature("gene_feature", _SO_GENE, 1, seq_len, "gene")

    # ── 2. CDS feature(s) ──
    coding_regions = [
        r for r in regions_sorted if r.region_type in ("exon", "cds")
    ]
    if coding_regions:
        for i, r in enumerate(coding_regions, start=1):
            # IR 0-indexed half-open → SBOL3 1-indexed inclusive.
            _emit_feature(
                f"cds_feature_{i}", "0000316",
                r.start + 1, r.end,
                f"cds_{i}",
            )
    else:
        # Prokaryotic / no-region case: whole sequence is the CDS.
        _emit_feature("cds_feature", "0000316", 1, seq_len, "cds")

    # ── 3. Other region types (UTR, promoter, intron, terminator, misc) ──
    other_idx = 0
    for r in regions_sorted:
        rtype = r.region_type
        if rtype in ("exon", "cds"):
            continue
        other_idx += 1
        so_id = _REGION_TYPE_TO_SO.get(rtype, _SO_MISC_REGION)
        _emit_feature(
            f"{rtype}_feature_{other_idx}", so_id,
            r.start + 1, r.end,
            f"{rtype}_{other_idx}",
        )

    # Drop the trailing blank line so the document ends cleanly.
    while feature_lines and feature_lines[-1] == "":
        feature_lines.pop()

    return feature_uris, feature_lines


def to_sbol3(ir_l0: IR_L0_GenomicDNA) -> str:
    """Emit SBOL3 (Synthetic Biology Open Language 3) RDF/Turtle format.

    SBOL3 represents genetic designs as Component objects with Sequence,
    Feature, and Constraint sub-objects. This is the synthetic biology
    interchange standard.

    Output format: RDF/Turtle (SBOL3 native format).

    The emitted document contains:

    * ``@prefix`` declarations for sbol, rdf, rdfs, biopax, so, and the
      DNA encoding vocabulary.
    * A top-level **Component** (``sbol:Component``) with ``displayId``,
      ``name``, ``description``, ``type`` (``biopax:DnaRegion``),
      ``hasFeature`` references, and a ``hasSequence`` reference.
    * A top-level **Sequence** (``sbol:Sequence``) with ``displayId``,
      ``encoding`` (IUPAC DNA), and the uppercase ``elements`` string.
    * A **SequenceFeature** per region (gene, CDS, exon, intron, UTR,
      promoter, terminator) with a ``Range`` location giving 1-based
      inclusive ``start``/``end`` coordinates.

    Coordinate conversion: IR uses 0-based half-open ``[start, end)``;
    SBOL3 Range uses 1-based inclusive ``start..end``.  Conversion is
    ``sbol_start = ir_start + 1``, ``sbol_end = ir_end``.

    Args:
        ir_l0: the IR-L0 genomic DNA object to serialise.

    Returns:
        An SBOL3 Turtle-format string.

    Example::

        >>> from biocompiler.ir import IR_L0_GenomicDNA
        >>> from biocompiler.ir.codegen import to_sbol3
        >>> ir = IR_L0_GenomicDNA(
        ...     sequence="ATGGCTTAA", regions=[],
        ...     organism="e_coli", gene_name="test",
        ... )
        >>> doc = to_sbol3(ir)
        >>> "@prefix sbol:" in doc
        True
        >>> "sbol:Component" in doc
        True
        >>> "sbol:elements" in doc and "ATGGCTTAA" in doc
        True
    """
    seq = ir_l0.sequence.upper()
    seq_len = len(seq)
    gene_name = ir_l0.gene_name or "gene"
    safe_id = _sbol3_safe_display_id(ir_l0.gene_name)

    # Stable, unique URIs based on gene_name + version.
    component_uri = f"{_SBOL3_BASE_URI}/{safe_id}/{_SBOL3_VERSION}"
    sequence_uri = f"{component_uri}/sequence"

    # Build the feature blocks first (we need the URIs for hasFeature).
    feature_uris, feature_lines = _build_sbol3_features(
        ir_l0, component_uri, sequence_uri,
    )

    lines: list[str] = []

    # ── @prefix declarations ──
    lines.append(f"@prefix sbol: <{_SBOL3_NS}> .")
    lines.append(f"@prefix rdf: <{_RDF_NS}> .")
    lines.append(f"@prefix rdfs: <{_RDFS_NS}> .")
    lines.append(f"@prefix biopax: <{_BIOPAX_NS}> .")
    lines.append(f"@prefix so: <{_SO_NS}> .")
    lines.append(
        f"@prefix enc: <http://www.cheminfo.org/cheminfo/emolecules/encoding#> ."
    )
    lines.append("")

    # ── Component (top-level) ──
    lines.append(f"<{component_uri}> a sbol:Component ;")
    lines.append(f'    sbol:displayId "{_escape_turtle_string(safe_id)}" ;')
    lines.append(f'    sbol:name "{_escape_turtle_string(gene_name)}" ;')
    lines.append('    sbol:description "Gene design compiled by BioCompiler" ;')
    lines.append("    sbol:type biopax:DnaRegion ;")
    # Multi-valued hasFeature: comma-separated object list.
    feat_objs = ", ".join(f"<{u}>" for u in feature_uris)
    lines.append(f"    sbol:hasFeature {feat_objs} ;")
    lines.append(f"    sbol:hasSequence <{sequence_uri}> .")
    lines.append("")

    # ── Sequence (top-level) ──
    lines.append(f"<{sequence_uri}> a sbol:Sequence ;")
    lines.append(f'    sbol:displayId "{_escape_turtle_string(safe_id)}_sequence" ;')
    lines.append("    sbol:encoding enc:DNA ;")
    lines.append(f'    sbol:elements "{seq}" .')
    lines.append("")

    # ── Feature + Range blocks (referenced by the Component) ──
    if feature_lines:
        lines.extend(feature_lines)
        # Ensure document ends with a single trailing newline.
        if lines[-1] != "":
            lines.append("")

    return "\n".join(lines) + "\n"


# ────────────────────────────────────────────────────────────────────
# Minimal GenBank parser (for round-trip verification)
# ────────────────────────────────────────────────────────────────────

# Matches an ORIGIN-sequence data line: leading whitespace, a number,
# then space-separated groups of bases.  We capture only the bases.
_ORIGIN_LINE_RE = re.compile(
    r"^\s*\d+\s+([ACGTNacgtn\s]+?)\s*$"
)


def parse_genbank_sequence(genbank_text: str) -> str:
    """Extract the DNA sequence from a GenBank record's ORIGIN section.

    This is a **minimal** parser — it only recovers the ORIGIN sequence,
    not the full feature table.  It exists to support round-trip
    verification (``IR_L0 → to_genbank() → parse_genbank_sequence() →
    IR_L0``) without depending on Biopython.

    Args:
        genbank_text: a GenBank-format string (as produced by
            :func:`to_genbank`).

    Returns:
        The DNA sequence (uppercase, no whitespace), or an empty string
        if no ORIGIN section is present.

    Example::

        >>> from biocompiler.ir import IR_L0_GenomicDNA
        >>> from biocompiler.ir.codegen import to_genbank, parse_genbank_sequence
        >>> ir = IR_L0_GenomicDNA(
        ...     sequence="ATGGCTTAA", regions=[],
        ...     organism="e_coli", gene_name="test",
        ... )
        >>> gb = to_genbank(ir)
        >>> parse_genbank_sequence(gb)
        'ATGGCTTAA'
    """
    in_origin = False
    bases: list[str] = []
    for raw_line in genbank_text.splitlines():
        if raw_line.startswith("ORIGIN"):
            in_origin = True
            continue
        if not in_origin:
            continue
        if raw_line.startswith("//"):
            break
        m = _ORIGIN_LINE_RE.match(raw_line)
        if m:
            # Strip all whitespace and uppercase.
            bases.append(re.sub(r"\s+", "", m.group(1)).upper())
    return "".join(bases)


def parse_genbank_features(genbank_text: str) -> list[tuple[str, str]]:
    """Extract ``(feature_key, location)`` pairs from a GenBank FEATURES table.

    Minimal parser — recovers the feature key and location string for
    each feature line (the lines with 5 leading spaces and a feature
    key followed by a location).  Qualifiers are not parsed.

    Args:
        genbank_text: a GenBank-format string.

    Returns:
        A list of ``(feature_key, location)`` tuples in document order.
    """
    in_features = False
    features: list[tuple[str, str]] = []
    feat_line_re = re.compile(r"^     (\S+)\s+(\S.*)$")
    for raw_line in genbank_text.splitlines():
        if raw_line.startswith("FEATURES"):
            in_features = True
            continue
        if not in_features:
            continue
        # ORIGIN / // ends the features table.
        if raw_line.startswith("ORIGIN") or raw_line.startswith("//"):
            break
        # Feature lines have exactly 5 leading spaces, then the key.
        m = feat_line_re.match(raw_line)
        if m:
            key = m.group(1)
            location = m.group(2).strip()
            features.append((key, location))
    return features
