"""
BioCompiler Application — Export service.

Orchestrates FASTA, GenBank, and SBOL3 export operations.
"""

import json as _json
import logging
from typing import Optional

from biocompiler.export.core import export_fasta as _export_fasta, export_genbank as _export_genbank
from biocompiler.shared.types import Certificate

logger = logging.getLogger(__name__)


def export_fasta(
    sequence: str,
    identifier: str = "BioCompiler_design",
    description: str = "",
    organism: str = "Homo_sapiens",
) -> str:
    """Export a DNA sequence in FASTA format.

    Returns:
        FASTA-formatted string.
    """
    return _export_fasta(
        sequence=sequence,
        identifier=identifier,
        description=description,
        organism=organism,
    )


def export_genbank(
    sequence: str,
    locus_name: str = "BIOCOMPILER",
    definition: str = "BioCompiler designed sequence",
    organism: str = "Homo_sapiens",
    gene_name: Optional[str] = None,
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
    certificate: Optional[Certificate] = None,
) -> str:
    """Export a DNA sequence in GenBank format with optional certificate embedding.

    Returns:
        GenBank-formatted string.
    """
    return _export_genbank(
        sequence=sequence,
        locus_name=locus_name,
        definition=definition,
        organism=organism,
        gene_name=gene_name,
        exon_boundaries=exon_boundaries,
        certificate=certificate,
    )


def export_sbol3(
    sequence: str,
    organism: str,
    gene_name: str,
    fmt: str = "sbol3",
) -> str:
    """Export a DNA sequence in SBOL3 format.

    Args:
        sequence: DNA sequence (already validated/uppercased).
        organism: Source organism.
        gene_name: Gene name for the SBOL3 component.
        fmt: Output format — 'sbol3' (XML) or 'sbol3json' (JSON).

    Returns:
        SBOL3-formatted string.
    """
    seq = sequence.upper()

    if fmt == "sbol3json":
        sbol3_content = _json.dumps({
            "SBOL3": True,
            "type": "http://sbols.org/v3#Component",
            "identity": f"https://biocompiler.org/{gene_name}",
            "displayId": gene_name,
            "sequence": {
                "type": "http://sbols.org/v3#Sequence",
                "identity": f"https://biocompiler.org/{gene_name}_seq",
                "displayId": f"{gene_name}_seq",
                "elements": seq,
                "encoding": "http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html",
            },
            "organism": organism,
        }, indent=2)
    else:
        sbol3_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
         xmlns:sbol="http://sbols.org/v3#"
         xmlns:dcterms="http://purl.org/dc/terms/">
  <sbol:Component rdf:about="https://biocompiler.org/{gene_name}">
    <dcterms:title>{gene_name}</dcterms:title>
    <sbol:displayId>{gene_name}</sbol:displayId>
    <sbol:hasSequence rdf:resource="https://biocompiler.org/{gene_name}_seq"/>
    <sbol:organism>{organism}</sbol:organism>
  </sbol:Component>
  <sbol:Sequence rdf:about="https://biocompiler.org/{gene_name}_seq">
    <sbol:displayId>{gene_name}_seq</sbol:displayId>
    <sbol:elements>{seq}</sbol:elements>
    <sbol:encoding rdf:resource="http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html"/>
  </sbol:Sequence>
</rdf:RDF>"""

    return sbol3_content


def export_batch_item(
    sequence: str,
    fmt: str = "fasta",
    identifier: str = "BioCompiler_design",
    description: str = "",
    organism: str = "Homo_sapiens",
    locus_name: str = "BIOCOMPILER",
    definition: str = "BioCompiler designed sequence",
    gene_name: Optional[str] = None,
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
) -> str:
    """Export a single batch item.

    Returns:
        Exported content string.

    Raises:
        ValueError: If format is unsupported.
    """
    if fmt == "fasta":
        return export_fasta(
            sequence=sequence,
            identifier=identifier,
            description=description,
            organism=organism,
        )
    elif fmt == "genbank":
        return export_genbank(
            sequence=sequence,
            locus_name=locus_name,
            definition=definition,
            organism=organism,
            gene_name=gene_name,
            exon_boundaries=exon_boundaries,
        )
    else:
        raise ValueError(f"Unsupported export format: {fmt}")
