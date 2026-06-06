"""
BioCompiler SBOL3 Export — Synthetic Biology Open Language 3 RDF/XML Export

SBOL3 is the standard for exchanging genetic design information. This module
exports optimized sequences as valid SBOL3 documents in RDF/XML format,
enabling interoperability with SBOL-compatible tools (Benchling, Cello,
SBOLDesigner, etc.).

Features:
- Valid SBOL3 RDF/XML document generation
- Component (DNA sequence) with Sequence object
- Design metadata (organism, protein, CAI, GC content)
- Optional certificate annotation embedding
- Falls back to manual RDF/XML generation when sbol3 library is not installed

SBOL3 specification: https://sbolstandard.org/sbol3/
"""

from __future__ import annotations

import logging
import uuid
from datetime import datetime, timezone
from xml.sax.saxutils import escape as xml_escape
from typing import Optional

from .scanner import gc_content
from .translation import translate, compute_cai
from . import __version__

logger = logging.getLogger(__name__)

__all__ = [
    "export_sbol3",
]

# SBOL3 namespace URIs
_SBOL3_NS = "http://sbols.org/v3#"
_SBOL2_NS = "http://sbols.org/v2#"
_RDF_NS = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
_RDFS_NS = "http://www.w3.org/2000/01/rdf-schema#"
_DCT_NS = "http://purl.org/dc/terms/"
_BIOL_NS = "http://www.biopax.org/release/biopax-level3.owl#"
_BIOCOMPILER_NS = "https://biocompiler.dev/vocab/"

# SBOL3 object types
_COMPONENT_TYPE = f"{_SBOL3_NS}Component"
_SEQUENCE_TYPE = f"{_SBOL3_NS}Sequence"
_COLLECTION_TYPE = f"{_SBOL3_NS}Collection"

# SBOL3 encoding
_IUPAC_DNA_ENCODING = f"{_SBOL3_NS}iupacNucleicAcid"

# SBOL3 roles
_CDS_ROLE = "http://www.biopax.org/release/biopax-level3.owl#DnaRegion"
_ENGINEERED_REGION_ROLE = f"{_SBOL3_NS}engineeredRegion"


def _generate_sbol_identity(prefix: str, name: str) -> str:
    """Generate a unique SBOL3 identity URI.

    SBOL3 requires globally unique identity URIs for all top-level objects.
    We use a UUID-based scheme to guarantee uniqueness.
    """
    unique_id = uuid.uuid4().hex[:12]
    return f"{prefix}/{name}/{unique_id}"


def _build_rdf_xml(
    component_identity: str,
    component_display_id: str,
    component_name: str,
    component_description: str,
    sequence_identity: str,
    sequence_display_id: str,
    sequence_elements: str,
    sequence_length: int,
    organism: str,
    protein: str,
    metadata: dict | None = None,
    certificate: dict | None = None,
) -> str:
    """Build a minimal valid SBOL3 RDF/XML document manually.

    This is the fallback path when the ``sbol3`` Python library is not installed.
    It produces a valid SBOL3 document containing:

    - A Component (type: dna, role: CDS / engineered region)
    - A Sequence (IUPAC DNA encoding)
    - Metadata annotations (organism, protein, CAI, GC)
    - Optional certificate annotation

    The output conforms to the SBOL3 RDF/XML serialization.
    """
    now = datetime.now(timezone.utc).isoformat()

    # Build annotation elements
    annotation_elements = []

    # Core BioCompiler metadata — use namespace prefixes declared on the root element
    annotation_elements.append(
        f'    <biol:organism>{xml_escape(organism)}</biol:organism>'
    )
    annotation_elements.append(
        f'    <biol:protein>{xml_escape(protein)}</biol:protein>'
    )
    annotation_elements.append(
        f'    <bc:targetOrganism>{xml_escape(organism)}</bc:targetOrganism>'
    )
    annotation_elements.append(
        f'    <bc:proteinSequence>{xml_escape(protein)}</bc:proteinSequence>'
    )
    annotation_elements.append(
        f'    <bc:version>{xml_escape(__version__)}</bc:version>'
    )
    annotation_elements.append(
        f'    <bc:optimizationDate>{xml_escape(now)}</bc:optimizationDate>'
    )
    annotation_elements.append(
        f'    <bc:sequenceLength>{sequence_length}</bc:sequenceLength>'
    )

    # Additional metadata
    if metadata:
        for key, value in metadata.items():
            if value is not None:
                # Sanitize key to be a valid XML element name
                safe_key = xml_escape(key).replace(" ", "_")
                annotation_elements.append(
                    f'    <bc:{safe_key}>{xml_escape(str(value))}</bc:{safe_key}>'
                )

    # Certificate annotation
    if certificate:
        cert_design_id = certificate.get("design_id", "unknown")
        cert_version = certificate.get("version", "unknown")
        cert_status = certificate.get("provenance", {}).get("overall_status", "unknown")

        annotation_elements.append(
            f'    <bc:hasCertificate>true</bc:hasCertificate>'
        )
        annotation_elements.append(
            f'    <bc:certificateDesignId>{xml_escape(cert_design_id)}</bc:certificateDesignId>'
        )
        annotation_elements.append(
            f'    <bc:certificateVersion>{xml_escape(cert_version)}</bc:certificateVersion>'
        )
        annotation_elements.append(
            f'    <bc:certificateStatus>{xml_escape(cert_status)}</bc:certificateStatus>'
        )

        # Embed the full certificate as a structured annotation
        import json
        cert_json = json.dumps(certificate, default=str)
        annotation_elements.append(
            f'    <bc:certificateData>{xml_escape(cert_json)}</bc:certificateData>'
        )

    annotations_block = "\n".join(annotation_elements)

    rdf_xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="{_RDF_NS}"
         xmlns:rdfs="{_RDFS_NS}"
         xmlns:sbol="{_SBOL3_NS}"
         xmlns:dct="{_DCT_NS}"
         xmlns:biol="{_BIOL_NS}"
         xmlns:bc="{_BIOCOMPILER_NS}">

  <!-- SBOL3 Component: the designed DNA construct -->
  <sbol:Component rdf:about="{xml_escape(component_identity)}">
    <sbol:displayId>{xml_escape(component_display_id)}</sbol:displayId>
    <sbol:name>{xml_escape(component_name)}</sbol:name>
    <sbol:description>{xml_escape(component_description)}</sbol:description>
    <sbol:hasSequence rdf:resource="{xml_escape(sequence_identity)}"/>
    <sbol:type rdf:resource="{_BIOL_NS}DnaRegion"/>
    <sbol:role rdf:resource="{_ENGINEERED_REGION_ROLE}"/>
    <sbol:role rdf:resource="{_CDS_ROLE}"/>
    <dct:creator>BioCompiler v{xml_escape(__version__)}</dct:creator>
    <dct:created>{xml_escape(now)}</dct:created>
{annotations_block}
  </sbol:Component>

  <!-- SBOL3 Sequence: the IUPAC DNA encoding -->
  <sbol:Sequence rdf:about="{xml_escape(sequence_identity)}">
    <sbol:displayId>{xml_escape(sequence_display_id)}</sbol:displayId>
    <sbol:elements>{xml_escape(sequence_elements)}</sbol:elements>
    <sbol:encoding rdf:resource="{_IUPAC_DNA_ENCODING}"/>
  </sbol:Sequence>

</rdf:RDF>
"""
    return rdf_xml


def _try_sbol3_library(
    sequence: str,
    protein: str,
    organism: str,
    certificate: dict | None = None,
    metadata: dict | None = None,
) -> str | None:
    """Attempt to export using the ``sbol3`` Python library.

    Returns the RDF/XML string on success, or None if the library
    is not available.
    """
    try:
        import sbol3 as sbol3_lib  # type: ignore[import-untyped]
    except ImportError:
        return None

    try:
        # Set up SBOL3 document
        doc = sbol3_lib.Document()

        # Create the Sequence object
        seq_display_id = f"{organism}_optimized_seq"
        sbol_seq = sbol3_lib.Sequence(seq_display_id)
        sbol_seq.elements = sequence.upper()
        sbol_seq.encoding = sbol3_lib.IUPAC_DNA_ENCODING
        doc.add(sbol_seq)

        # Create the Component object
        comp_display_id = f"{organism}_optimized_component"
        sbol_comp = sbol3_lib.Component(comp_display_id, sbol3_lib.SBO_DNA)
        sbol_comp.sequences = [sbol_seq.identity]
        sbol_comp.name = f"BioCompiler optimized sequence for {organism}"
        sbol_comp.description = (
            f"Codon-optimized DNA sequence for {protein} "
            f"in {organism} (BioCompiler v{__version__})"
        )

        # Add roles
        sbol_comp.roles = [sbol3_lib.SBO_CDS, sbol3_lib.SBO_ENGINEERED_REGION]

        # Add annotations for metadata
        if metadata:
            for key, value in metadata.items():
                if value is not None:
                    annotation = sbol3_lib.TextProperty(
                        sbol_comp, f"{_BIOCOMPILER_NS}{key}", 0, 1
                    )
                    annotation.set(str(value))

        # Add certificate annotation
        if certificate:
            import json
            cert_annotation = sbol3_lib.TextProperty(
                sbol_comp, f"{_BIOCOMPILER_NS}certificateData", 0, 1
            )
            cert_annotation.set(json.dumps(certificate, default=str))

        doc.add(sbol_comp)

        # Write to RDF/XML string
        return doc.write_string(sbol3_lib.NTRIPLES)

    except Exception as e:
        logger.warning("sbol3 library export failed, falling back to manual RDF/XML: %s", e)
        return None


def export_sbol3(
    sequence: str,
    protein: str,
    organism: str,
    certificate: dict | None = None,
    metadata: dict | None = None,
) -> str:
    """Export a designed DNA sequence in SBOL3 RDF/XML format.

    SBOL3 (Synthetic Biology Open Language 3) is the standard for exchanging
    genetic design information. This function generates a valid SBOL3 document
    containing the optimized DNA sequence, associated metadata, and optionally
    a BioCompiler certificate as an annotation.

    The document includes:

    - **Component**: A DNA component representing the designed sequence,
      with type ``DnaRegion`` and roles ``CDS`` and ``engineeredRegion``.
    - **Sequence**: The IUPAC DNA encoding of the component.
    - **Annotations**: Organism, protein, CAI, GC content, optimization date,
      and BioCompiler version.
    - **Certificate** (optional): The full BioCompiler certificate embedded
      as a structured annotation, enabling downstream verification.

    If the ``sbol3`` Python library is installed, it is used for document
    generation. Otherwise, a minimal but valid RDF/XML document is generated
    manually.

    Args:
        sequence: DNA sequence string (A, C, G, T).
        protein: Target protein sequence (single-letter amino acid codes).
        organism: Target organism name (e.g., 'Escherichia_coli').
        certificate: Optional certificate dict to embed as an annotation.
            Should contain 'design_id', 'version', 'types', 'provenance'.
        metadata: Optional dict of additional metadata to include as
            annotations. Common keys: 'cai', 'gc_content', 'gene_name'.

    Returns:
        SBOL3 RDF/XML document as a string.

    Example::

        from biocompiler.export_sbol import export_sbol3
        from biocompiler import optimize_sequence

        result = optimize_sequence("MSKGEELFTG", organism="Escherichia_coli")
        sbol_doc = export_sbol3(
            sequence=result.sequence,
            protein=result.protein,
            organism="Escherichia_coli",
            certificate=result.certificate_text,
            metadata={"cai": result.cai, "gc_content": result.gc_content},
        )
        with open("design.xml", "w") as f:
            f.write(sbol_doc)
    """
    if not sequence:
        raise ValueError("Sequence must not be empty")
    if not protein:
        raise ValueError("Protein must not be empty")
    if not organism:
        raise ValueError("Organism must not be empty")

    # Normalize sequence
    seq = sequence.upper().replace(" ", "")
    gc = gc_content(seq)

    # Try to compute CAI for metadata
    cai_value = None
    if len(seq) >= 3 and len(seq) % 3 == 0:
        try:
            cai_value = compute_cai(seq, organism=organism)
        except Exception:
            pass

    # Build metadata dict (merge user-provided with computed values)
    effective_metadata: dict = {
        "gcContent": f"{gc:.4f}",
    }
    if cai_value is not None:
        effective_metadata["cai"] = f"{cai_value:.4f}"
    if metadata:
        effective_metadata.update(metadata)

    # Try sbol3 library first
    result = _try_sbol3_library(seq, protein, organism, certificate, effective_metadata)
    if result is not None:
        return result

    # Fallback: manual RDF/XML generation
    logger.info("sbol3 library not available; generating SBOL3 RDF/XML manually")

    # Generate SBOL3 identities
    prefix = f"https://biocompiler.dev/design"
    component_display_id = f"{organism}_optimized_component"
    sequence_display_id = f"{organism}_optimized_seq"
    component_identity = _generate_sbol_identity(prefix, component_display_id)
    sequence_identity = _generate_sbol_identity(prefix, sequence_display_id)

    component_name = f"BioCompiler optimized sequence for {organism}"
    component_description = (
        f"Codon-optimized DNA sequence for {protein} "
        f"in {organism} (BioCompiler v{__version__})"
    )

    return _build_rdf_xml(
        component_identity=component_identity,
        component_display_id=component_display_id,
        component_name=component_name,
        component_description=component_description,
        sequence_identity=sequence_identity,
        sequence_display_id=sequence_display_id,
        sequence_elements=seq,
        sequence_length=len(seq),
        organism=organism,
        protein=protein,
        metadata=effective_metadata,
        certificate=certificate,
    )
