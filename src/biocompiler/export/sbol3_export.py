"""
BioCompiler SBOL3 Export — Synthetic Biology Open Language 3 RDF/XML Export
===========================================================================

Pure-Python SBOL3 (Synthetic Biology Open Language 3) document generator.
No external SBOL library dependency required.

Generates valid SBOL3 RDF/XML documents containing:
- **Component**: DNA component representing the designed sequence, with type
  ``DnaRegion`` (BioPAX) and functional description roles (SO:0000316 CDS,
  SO:0000704 gene).
- **Sequence**: IUPAC DNA encoding linked to the component.
- **Constraint**: Structural constraints (e.g., sequence length) linking
  the component to its sequence.
- **Annotations**: Organism, protein name, CAI, GC content, optimization
  date, and BioCompiler version as Dublin Core / custom vocabulary terms.

SBOL3 specification: https://sbolstandard.org/sbol3/
"""

from __future__ import annotations

import json
import logging
import uuid
from datetime import datetime, timezone
from typing import Any, Optional
from xml.sax.saxutils import escape as xml_escape

from .. import __version__

logger = logging.getLogger(__name__)

__all__ = [
    "generate_sbol3_xml",
    "generate_sbol3_json",
    "export_sbol3",
    "_generate_identity",
    "BC_NS",
    "BIOPAX_NS",
    "DCT_NS",
    "DNAREGION_TYPE",
    "IUPAC_DNA_ENCODING",
    "RDF_NS",
    "SBOL3_NS",
    "SO_CDS",
    "SO_GENE",
]

# ─── SBOL3 Namespace Constants ────────────────────────────────────

SBOL3_NS = "http://sbols.org/v3#"
RDF_NS = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
RDFS_NS = "http://www.w3.org/2000/01/rdf-schema#"
DCT_NS = "http://purl.org/dc/terms/"
BIOPAX_NS = "http://www.biopax.org/release/biopax-level3.owl#"
PROV_NS = "http://www.w3.org/ns/prov#"
SO_NS = "http://sequenceontology.org/resource/SO:"
BC_NS = "https://biocompiler.dev/vocab/"

# SBOL3 encoding URIs
IUPAC_DNA_ENCODING = "http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html"
IUPAC_PROTEIN_ENCODING = "http://www.chem.qmul.ac.uk/iubmb/misc/aaseq.html"

# BioPAX type URIs
DNAREGION_TYPE = f"{BIOPAX_NS}DnaRegion"

# Sequence Ontology role URIs
SO_CDS = f"{SO_NS}0000316"
SO_GENE = f"{SO_NS}0000704"

# SBOL3 constraint types
MEETS_CONSTRAINT = f"{SBOL3_NS}meetsConstraint"


def _generate_identity(base_uri: str, display_id: str) -> str:
    """Generate a unique SBOL3 identity URI.

    SBOL3 requires globally unique identity URIs for all top-level objects.
    We use a UUID-based scheme to guarantee uniqueness.
    """
    uid = uuid.uuid4().hex[:12]
    return f"{base_uri}/{display_id}/{uid}"


def _format_timestamp() -> str:
    """Return the current UTC timestamp in ISO 8601 format."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def generate_sbol3_xml(
    sequence: str,
    protein_name: str,
    organism: str,
    cai: Optional[float] = None,
    gc_content: Optional[float] = None,
    satisfied_predicates: Optional[list[str]] = None,
    failed_predicates: Optional[list[str]] = None,
    fallback_used: bool = False,
    base_uri: str = "https://biocompiler.dev/design",
    certificate: Optional[dict[str, Any]] = None,
) -> str:
    """Generate an SBOL3-compliant RDF/XML document.

    Produces a valid SBOL3 document containing a Component (DnaRegion),
    Sequence (IUPAC DNA), and Constraint objects, along with functional
    description annotations and optimization result metadata.

    Args:
        sequence: DNA sequence string (A, C, G, T).
        protein_name: Name of the target protein (used as display ID).
        organism: Target organism name (e.g., 'Escherichia_coli').
        cai: Optional Codon Adaptation Index value.
        gc_content: Optional GC content fraction.
        satisfied_predicates: List of satisfied predicate names.
        failed_predicates: List of failed predicate names.
        fallback_used: Whether the optimizer used fallback.
        base_uri: Base URI for SBOL3 identity URIs.
        certificate: Optional certificate dict to embed as annotation.

    Returns:
        SBOL3 RDF/XML document as a string.

    Raises:
        ValueError: If sequence, protein_name, or organism is empty.
    """
    if not sequence or not sequence.strip():
        raise ValueError("Sequence must not be empty")
    if not protein_name or not protein_name.strip():
        raise ValueError("Protein name must not be empty")
    if not organism or not organism.strip():
        raise ValueError("Organism must not be empty")

    seq = sequence.upper().replace(" ", "")

    # Sanitise display IDs: replace spaces and special chars
    safe_protein = protein_name.replace(" ", "_").replace("(", "").replace(")", "")
    safe_organism = organism.replace(" ", "_")

    # Generate unique identities
    component_id = _generate_identity(base_uri, safe_protein)
    sequence_id = _generate_identity(base_uri, f"{safe_protein}_seq")
    constraint_id = _generate_identity(base_uri, f"{safe_protein}_constraint")

    now = _format_timestamp()

    # Build annotation elements
    annotations: list[str] = []

    # Dublin Core annotations
    annotations.append(
        f"    <dcterms:creator>BioCompiler v{xml_escape(__version__)}</dcterms:creator>"
    )
    annotations.append(
        f"    <dcterms:created>{xml_escape(now)}</dcterms:created>"
    )
    annotations.append(
        f'    <dcterms:title>{xml_escape(safe_protein)}</dcterms:title>'
    )

    # Functional description
    description = (
        f"Codon-optimized DNA sequence for {protein_name} "
        f"in {organism} (BioCompiler v{__version__})"
    )
    annotations.append(
        f"    <dcterms:description>{xml_escape(description)}</dcterms:description>"
    )

    # BioCompiler vocabulary annotations
    annotations.append(
        f"    <bc:targetOrganism>{xml_escape(organism)}</bc:targetOrganism>"
    )
    annotations.append(
        f"    <bc:proteinName>{xml_escape(protein_name)}</bc:proteinName>"
    )
    annotations.append(
        f"    <bc:biocompilerVersion>{xml_escape(__version__)}</bc:biocompilerVersion>"
    )
    annotations.append(
        f"    <bc:optimizationDate>{xml_escape(now)}</bc:optimizationDate>"
    )
    annotations.append(
        f"    <bc:sequenceLength>{len(seq)}</bc:sequenceLength>"
    )

    # Optimization result annotations
    if cai is not None:
        annotations.append(
            f'    <bc:cai>{cai:.4f}</bc:cai>'
        )
    if gc_content is not None:
        annotations.append(
            f'    <bc:gcContent>{gc_content:.4f}</bc:gcContent>'
        )
    if satisfied_predicates is not None:
        annotations.append(
            f'    <bc:satisfiedPredicates>{xml_escape(json.dumps(satisfied_predicates))}</bc:satisfiedPredicates>'
        )
    if failed_predicates is not None:
        annotations.append(
            f'    <bc:failedPredicates>{xml_escape(json.dumps(failed_predicates))}</bc:failedPredicates>'
        )
    if fallback_used:
        annotations.append(
            f'    <bc:fallbackUsed>true</bc:fallbackUsed>'
        )

    # Certificate annotation
    if certificate:
        annotations.append(
            '    <bc:hasCertificate>true</bc:hasCertificate>'
        )
        cert_design_id = certificate.get("design_id", "unknown")
        annotations.append(
            f'    <bc:certificateDesignId>{xml_escape(cert_design_id)}</bc:certificateDesignId>'
        )
        cert_json = json.dumps(certificate, default=str)
        annotations.append(
            f'    <bc:certificateData>{xml_escape(cert_json)}</bc:certificateData>'
        )

    annotations_block = "\n".join(annotations)

    # Build the full SBOL3 RDF/XML document
    rdf_xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="{RDF_NS}"
         xmlns:rdfs="{RDFS_NS}"
         xmlns:sbol="{SBOL3_NS}"
         xmlns:dcterms="{DCT_NS}"
         xmlns:biopax="{BIOPAX_NS}"
         xmlns:prov="{PROV_NS}"
         xmlns:so="{SO_NS}"
         xmlns:bc="{BC_NS}">

  <!-- SBOL3 Component: the designed DNA construct -->
  <sbol:Component rdf:about="{xml_escape(component_id)}">
    <sbol:displayId>{xml_escape(safe_protein)}</sbol:displayId>
    <sbol:name>{xml_escape(protein_name)}</sbol:name>
    <sbol:description>{xml_escape(description)}</sbol:description>
    <sbol:type rdf:resource="{DNAREGION_TYPE}"/>
    <sbol:role rdf:resource="{SO_GENE}"/>
    <sbol:role rdf:resource="{SO_CDS}"/>
    <sbol:hasSequence rdf:resource="{xml_escape(sequence_id)}"/>
    <sbol:hasConstraint rdf:resource="{xml_escape(constraint_id)}"/>
{annotations_block}
  </sbol:Component>

  <!-- SBOL3 Sequence: the IUPAC DNA encoding -->
  <sbol:Sequence rdf:about="{xml_escape(sequence_id)}">
    <sbol:displayId>{xml_escape(safe_protein)}_seq</sbol:displayId>
    <sbol:name>{xml_escape(protein_name)} sequence</sbol:name>
    <sbol:description>IUPAC DNA encoding of {xml_escape(protein_name)}</sbol:description>
    <sbol:elements>{xml_escape(seq)}</sbol:elements>
    <sbol:encoding rdf:resource="{IUPAC_DNA_ENCODING}"/>
  </sbol:Sequence>

  <!-- SBOL3 Constraint: structural constraint linking component to sequence -->
  <sbol:Constraint rdf:about="{xml_escape(constraint_id)}">
    <sbol:displayId>{xml_escape(safe_protein)}_constraint</sbol:displayId>
    <sbol:subject rdf:resource="{xml_escape(component_id)}"/>
    <sbol:object rdf:resource="{xml_escape(sequence_id)}"/>
    <sbol:restriction rdf:resource="{MEETS_CONSTRAINT}"/>
  </sbol:Constraint>

</rdf:RDF>
"""
    return rdf_xml


def generate_sbol3_json(
    sequence: str,
    protein_name: str,
    organism: str,
    cai: Optional[float] = None,
    gc_content: Optional[float] = None,
    satisfied_predicates: Optional[list[str]] = None,
    failed_predicates: Optional[list[str]] = None,
    fallback_used: bool = False,
    base_uri: str = "https://biocompiler.dev/design",
) -> str:
    """Generate an SBOL3-compliant JSON-LD document.

    Produces a simplified JSON representation of the SBOL3 data model
    for tools that prefer JSON over RDF/XML.

    Args:
        sequence: DNA sequence string (A, C, G, T).
        protein_name: Name of the target protein.
        organism: Target organism name.
        cai: Optional Codon Adaptation Index value.
        gc_content: Optional GC content fraction.
        satisfied_predicates: List of satisfied predicate names.
        failed_predicates: List of failed predicate names.
        fallback_used: Whether the optimizer used fallback.
        base_uri: Base URI for SBOL3 identity URIs.

    Returns:
        JSON-LD document as a string.

    Raises:
        ValueError: If sequence, protein_name, or organism is empty.
    """
    if not sequence or not sequence.strip():
        raise ValueError("Sequence must not be empty")
    if not protein_name or not protein_name.strip():
        raise ValueError("Protein name must not be empty")
    if not organism or not organism.strip():
        raise ValueError("Organism must not be empty")

    seq = sequence.upper().replace(" ", "")
    safe_protein = protein_name.replace(" ", "_")
    now = _format_timestamp()

    component_id = _generate_identity(base_uri, safe_protein)
    sequence_id = _generate_identity(base_uri, f"{safe_protein}_seq")
    constraint_id = _generate_identity(base_uri, f"{safe_protein}_constraint")

    doc: dict[str, Any] = {
        "@context": {
            "sbol": SBOL3_NS,
            "rdf": RDF_NS,
            "rdfs": RDFS_NS,
            "dcterms": DCT_NS,
            "biopax": BIOPAX_NS,
            "prov": PROV_NS,
            "so": SO_NS,
            "bc": BC_NS,
        },
        "@graph": [
            {
                "@id": component_id,
                "@type": "sbol:Component",
                "sbol:displayId": safe_protein,
                "sbol:name": protein_name,
                "sbol:description": (
                    f"Codon-optimized DNA sequence for {protein_name} "
                    f"in {organism} (BioCompiler v{__version__})"
                ),
                "sbol:type": {"@id": DNAREGION_TYPE},
                "sbol:role": [{"@id": SO_GENE}, {"@id": SO_CDS}],
                "sbol:hasSequence": {"@id": sequence_id},
                "sbol:hasConstraint": {"@id": constraint_id},
                "dcterms:creator": f"BioCompiler v{__version__}",
                "dcterms:created": now,
                "bc:targetOrganism": organism,
                "bc:proteinName": protein_name,
                "bc:sequenceLength": len(seq),
            },
            {
                "@id": sequence_id,
                "@type": "sbol:Sequence",
                "sbol:displayId": f"{safe_protein}_seq",
                "sbol:name": f"{protein_name} sequence",
                "sbol:elements": seq,
                "sbol:encoding": {"@id": IUPAC_DNA_ENCODING},
            },
            {
                "@id": constraint_id,
                "@type": "sbol:Constraint",
                "sbol:displayId": f"{safe_protein}_constraint",
                "sbol:subject": {"@id": component_id},
                "sbol:object": {"@id": sequence_id},
                "sbol:restriction": {"@id": MEETS_CONSTRAINT},
            },
        ],
    }

    # Add optional optimization metadata
    component_entry = doc["@graph"][0]
    if cai is not None:
        component_entry["bc:cai"] = round(cai, 4)
    if gc_content is not None:
        component_entry["bc:gcContent"] = round(gc_content, 4)
    if satisfied_predicates is not None:
        component_entry["bc:satisfiedPredicates"] = satisfied_predicates
    if failed_predicates is not None:
        component_entry["bc:failedPredicates"] = failed_predicates
    if fallback_used:
        component_entry["bc:fallbackUsed"] = True

    return json.dumps(doc, indent=2, sort_keys=False)


def export_sbol3(
    sequence: str,
    protein_name: str,
    organism: str,
    optimization_result: Optional[Any] = None,
    cai: Optional[float] = None,
    gc_content: Optional[float] = None,
    satisfied_predicates: Optional[list[str]] = None,
    failed_predicates: Optional[list[str]] = None,
    fallback_used: bool = False,
    format: str = "sbol3",
    base_uri: str = "https://biocompiler.dev/design",
    certificate: Optional[dict[str, Any]] = None,
) -> str:
    """Export a designed DNA sequence in SBOL3 format.

    This is the primary public API for SBOL3 export. It accepts either
    explicit parameters or an OptimizationResult object and generates
    a valid SBOL3 document.

    When an ``OptimizationResult`` is provided, its fields (CAI, GC content,
    satisfied/failed predicates, etc.) are automatically extracted and
    included as annotations.

    Args:
        sequence: DNA sequence string (A, C, G, T).
        protein_name: Name of the target protein.
        organism: Target organism name (e.g., 'Escherichia_coli').
        optimization_result: Optional OptimizationResult object.
            When provided, its metrics are extracted automatically.
        cai: Codon Adaptation Index (overridden by optimization_result if given).
        gc_content: GC content fraction (overridden by optimization_result if given).
        satisfied_predicates: List of satisfied predicates.
        failed_predicates: List of failed predicates.
        fallback_used: Whether the optimizer used fallback.
        format: Output format — ``"sbol3"`` (XML, default) or ``"sbol3json"`` (JSON-LD).
        base_uri: Base URI for SBOL3 identity URIs.
        certificate: Optional certificate dict to embed as annotation.

    Returns:
        SBOL3 document as a string (XML or JSON-LD).

    Raises:
        ValueError: If required fields are empty or format is invalid.
        TypeError: If optimization_result is not an OptimizationResult.
    """
    if format not in ("sbol3", "sbol3json"):
        raise ValueError(
            f"Unsupported SBOL format: {format!r}. Use 'sbol3' (XML) or 'sbol3json' (JSON-LD)."
        )

    # Extract fields from OptimizationResult if provided
    eff_cai = cai
    eff_gc = gc_content
    eff_satisfied = satisfied_predicates
    eff_failed = failed_predicates
    eff_fallback = fallback_used

    if optimization_result is not None:
        from ..optimization import OptimizationResult

        if not isinstance(optimization_result, OptimizationResult):
            raise TypeError(
                f"Expected OptimizationResult, got {type(optimization_result).__name__}. "
                f"Use optimize_sequence() to produce an OptimizationResult."
            )

        # OptimizationResult fields take precedence unless explicitly overridden
        if eff_cai is None:
            eff_cai = optimization_result.cai
        if eff_gc is None:
            eff_gc = optimization_result.gc_content
        if eff_satisfied is None:
            eff_satisfied = optimization_result.satisfied_predicates
        if eff_failed is None:
            eff_failed = optimization_result.failed_predicates
        if not eff_fallback:
            eff_fallback = optimization_result.fallback_used

    # Compute GC content from sequence if not provided
    if eff_gc is None and sequence:
        from ..scanner import gc_content as _gc_content
        eff_gc = _gc_content(sequence.upper().replace(" ", ""))

    if format == "sbol3json":
        return generate_sbol3_json(
            sequence=sequence,
            protein_name=protein_name,
            organism=organism,
            cai=eff_cai,
            gc_content=eff_gc,
            satisfied_predicates=eff_satisfied,
            failed_predicates=eff_failed,
            fallback_used=eff_fallback,
            base_uri=base_uri,
        )

    return generate_sbol3_xml(
        sequence=sequence,
        protein_name=protein_name,
        organism=organism,
        cai=eff_cai,
        gc_content=eff_gc,
        satisfied_predicates=eff_satisfied,
        failed_predicates=eff_failed,
        fallback_used=eff_fallback,
        base_uri=base_uri,
        certificate=certificate,
    )
