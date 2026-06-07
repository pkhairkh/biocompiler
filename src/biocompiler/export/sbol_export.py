"""
BioCompiler SBOL3 Export — Synthetic Biology Open Language v3
=============================================================

Export BioCompiler optimization results as SBOL3 (Synthetic Biology Open Language)
documents for interoperability with SBOL-compatible tools (Benchling, SynBioHub,
Pigeon, etc.).

SBOL3 format details:
- XML namespace: http://sbols.org/v3#
- Component types from SO (Sequence Ontology): SO:0000316 (CDS),
  SO:0000167 (promoter), SO:0000141 (terminator)
- Roles from SO: SO:0000316 (CDS), SO:0000167 (promoter), SO:0000141 (terminator)
- CAI/GC as SBOL Measure objects with OM (Ontology of Units of Measure) units
- Provenance metadata using PROV-O vocabulary

This module provides a pure-Python SBOL3 XML generator — no external SBOL
library dependency required.  The output conforms to the SBOL3 specification
and can be validated with sbol-validate or loaded into any SBOL3-compliant
tool.
"""

from __future__ import annotations

import json
import logging
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional
from xml.etree.ElementTree import Element, SubElement, ElementTree, indent, register_namespace

from .. import __version__

logger = logging.getLogger(__name__)

__all__ = [
    "SBOLComponent",
    "export_sbol",
    "export_sbol_collection",
    "SBOL3_NS",
    "RDF_NS",
    "SO_CDS",
    "SO_NS",
    "SO_PROMOTER",
    "SO_TERMINATOR",
    "_resolve_role_uri",
]

# ─── SBOL3 Namespace Constants ────────────────────────────────────

SBOL3_NS = "http://sbols.org/v3#"
RDF_NS = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
PROV_NS = "http://www.w3.org/ns/prov#"
SO_NS = "http://sequenceontology.org/resource/SO:"
OM_NS = "http://www.ontology-of-units-of-measure.org/resource/om-2/"
DCT_NS = "http://purl.org/dc/terms/"

# Register namespace prefixes so ElementTree uses human-readable names
# instead of auto-generated ns0, ns1, etc.  Must come after the
# constants are defined.
register_namespace("sbol3", SBOL3_NS)
register_namespace("rdf", RDF_NS)
register_namespace("prov", PROV_NS)
register_namespace("so", SO_NS)
register_namespace("om", OM_NS)
register_namespace("dct", DCT_NS)
register_namespace("dc", "http://purl.org/dc/elements/1.1/")

# Sequence Ontology identifiers
SO_CDS = f"{SO_NS}0000316"       # Coding Sequence
SO_PROMOTER = f"{SO_NS}0000167"  # Promoter
SO_TERMINATOR = f"{SO_NS}0000141"  # Terminator
SO_RBS = f"{SO_NS}0000139"       # Ribosome Binding Site
SO_GENE = f"{SO_NS}0000704"      # Gene
SO_DNA = f"{SO_NS}0000347"       # DNA molecule type

# OM (Ontology of Units of Measure) identifiers
OM_FRACTION = f"{OM_NS}fraction"
OM_DIMENSIONLESS = f"{OM_NS}dimensionless"

# Role mapping: BioCompiler concept → SO URI
_ROLE_MAP = {
    "CDS": SO_CDS,
    "cds": SO_CDS,
    "coding_sequence": SO_CDS,
    "promoter": SO_PROMOTER,
    "Promoter": SO_PROMOTER,
    "terminator": SO_TERMINATOR,
    "Terminator": SO_TERMINATOR,
    "RBS": SO_RBS,
    "rbs": SO_RBS,
    "ribosome_binding_site": SO_RBS,
    "gene": SO_GENE,
    "Gene": SO_GENE,
}


@dataclass
class SBOLComponent:
    """An SBOL3 Component definition.

    Represents a single biological component (e.g., a gene, CDS, promoter)
    in the SBOL3 data model.

    Attributes:
        identity: Unique URI identifying this component.
        display_id: Human-readable identifier (no spaces).
        component_type: Molecular type — ``"DNA"`` or ``"Protein"``.
        sequence: Nucleotide or amino acid sequence string.
        roles: SBOL role URIs (e.g., SO:0000316 for CDS).
        description: Free-text description.
    """

    identity: str
    display_id: str
    component_type: str  # "DNA", "Protein"
    sequence: str
    roles: list[str]  # SBOL roles like "CDS", "promoter", "terminator"
    description: str = ""

    def __post_init__(self) -> None:
        """Validate and normalise SBOLComponent fields."""
        self.component_type = self.component_type.upper()
        if self.component_type not in ("DNA", "PROTEIN"):
            raise ValueError(
                f"component_type must be 'DNA' or 'Protein', got {self.component_type!r}"
            )
        # Normalise display_id: replace spaces with underscores
        self.display_id = self.display_id.replace(" ", "_")
        if not self.display_id:
            self.display_id = f"component_{uuid.uuid4().hex[:8]}"


# ─── Internal Helpers ──────────────────────────────────────────────


def _make_identity(base_uri: str, display_id: str) -> str:
    """Construct a unique SBOL identity URI."""
    return f"{base_uri}/{display_id}"


def _resolve_role_uri(role: str) -> str:
    """Map a human-readable role name to its SO URI, or pass through if already a URI."""
    if role.startswith("http://") or role.startswith("https://"):
        return role
    return _ROLE_MAP.get(role, f"{SO_NS}0000000")  # Unknown maps to SO root


def _component_type_uri(comp_type: str) -> str:
    """Map component_type ('DNA'/'Protein') to an SBOL type URI."""
    if comp_type.upper() == "DNA":
        return SO_DNA
    # Protein — use a standard biomaterial type
    return f"{SO_NS}0000297"  # polypeptide


def _build_component_element(
    comp: SBOLComponent,
    base_uri: str,
    cai: float | None = None,
    gc: float | None = None,
    organism: str = "",
    provenance_id: str = "",
) -> Element:
    """Build an RDF/XML ``sbol:Component`` element for a single SBOLComponent.

    This includes:
    - Component with type and role annotations
    - Sequence object linked to the component
    - Measure objects for CAI and GC content
    - Provenance metadata (organism, BioCompiler version)
    """
    identity = comp.identity or _make_identity(base_uri, comp.display_id)

    # ── Component element ──
    comp_elem = Element(f"{{{SBOL3_NS}}}Component")
    comp_elem.set(f"{{{RDF_NS}}}about", identity)

    # displayId
    did = SubElement(comp_elem, f"{{{SBOL3_NS}}}displayId")
    did.text = comp.display_id

    # type (molecular type)
    type_elem = SubElement(comp_elem, f"{{{SBOL3_NS}}}type")
    type_elem.set(f"{{{RDF_NS}}}resource", _component_type_uri(comp.component_type))

    # roles
    for role_str in comp.roles:
        role_elem = SubElement(comp_elem, f"{{{SBOL3_NS}}}role")
        role_elem.set(f"{{{RDF_NS}}}resource", _resolve_role_uri(role_str))

    # description
    if comp.description:
        desc_elem = SubElement(comp_elem, f"{{{DCT_NS}}}description")
        desc_elem.text = comp.description

    # ── Sequence element (linked to component) ──
    if comp.sequence:
        seq_identity = f"{identity}/sequence"
        seq_elem = SubElement(comp_elem, f"{{{SBOL3_NS}}}hasSequence")
        seq_elem.set(f"{{{RDF_NS}}}resource", seq_identity)

        # We'll also create a separate Sequence element — that happens in the
        # parent builder, but we record the reference here.

    # ── Measure objects for CAI and GC ──
    measures_added = False
    if cai is not None:
        _add_measure(comp_elem, identity, "CAI", cai, OM_DIMENSIONLESS, "cai_measure")
        measures_added = True
    if gc is not None:
        _add_measure(comp_elem, identity, "GC_content", gc, OM_FRACTION, "gc_measure")
        measures_added = True

    # ── Provenance annotations ──
    if organism:
        org_elem = SubElement(comp_elem, f"{{{SBOL3_NS}}}wasDerivedFrom")
        org_elem.set(
            f"{{{RDF_NS}}}resource",
            f"http://identifiers.org/taxonomy/{organism.replace('_', ' ')}",
        )

    if provenance_id:
        prov_elem = SubElement(comp_elem, f"{{{PROV_NS}}}wasGeneratedBy")
        prov_elem.set(f"{{{RDF_NS}}}resource", f"{base_uri}/Activity/{provenance_id}")

    # BioCompiler version annotation
    bc_elem = SubElement(comp_elem, f"{{{DCT_NS}}}creator")
    bc_elem.text = f"BioCompiler v{__version__}"

    return comp_elem


def _add_measure(
    parent: Element,
    owner_identity: str,
    label: str,
    value: float,
    unit_uri: str,
    display_id: str,
) -> None:
    """Attach an SBOL Measure element to a parent Component element."""
    measure_identity = f"{owner_identity}/{display_id}"

    measure_wrapper = SubElement(parent, f"{{{SBOL3_NS}}}hasMeasure")

    measure_elem = SubElement(measure_wrapper, f"{{{SBOL3_NS}}}Measure")
    measure_elem.set(f"{{{RDF_NS}}}about", measure_identity)

    did = SubElement(measure_elem, f"{{{SBOL3_NS}}}displayId")
    did.text = display_id

    val = SubElement(measure_elem, f"{{{SBOL3_NS}}}value")
    val.set(f"{{{RDF_NS}}}datatype", "http://www.w3.org/2001/XMLSchema#decimal")
    val.text = f"{value:.6f}"

    unit = SubElement(measure_elem, f"{{{SBOL3_NS}}}unit")
    unit.set(f"{{{RDF_NS}}}resource", unit_uri)

    # Label for human readability
    name_elem = SubElement(measure_elem, f"{{{DCT_NS}}}title")
    name_elem.text = label


def _build_sequence_element(
    comp: SBOLComponent,
    base_uri: str,
) -> Element | None:
    """Build an RDF/XML ``sbol:Sequence`` element for a component's DNA sequence."""
    if not comp.sequence:
        return None

    identity = comp.identity or _make_identity(base_uri, comp.display_id)
    seq_identity = f"{identity}/sequence"

    seq_elem = Element(f"{{{SBOL3_NS}}}Sequence")
    seq_elem.set(f"{{{RDF_NS}}}about", seq_identity)

    did = SubElement(seq_elem, f"{{{SBOL3_NS}}}displayId")
    did.text = f"{comp.display_id}_seq"

    encoding = SubElement(seq_elem, f"{{{SBOL3_NS}}}encoding")
    if comp.component_type.upper() == "DNA":
        encoding.set(
            f"{{{RDF_NS}}}resource",
            "http://www.chem.qmul.ac.uk/iupac/DNA/",
        )
    else:
        encoding.set(
            f"{{{RDF_NS}}}resource",
            "http://www.chem.qmul.ac.uk/iupac/AA/",
        )

    elements = SubElement(seq_elem, f"{{{SBOL3_NS}}}elements")
    elements.text = comp.sequence.upper()

    return seq_elem


def _build_activity_element(
    base_uri: str,
    activity_id: str,
    organism: str = "",
) -> Element:
    """Build a PROV-O Activity element describing the BioCompiler optimization."""
    activity_identity = f"{base_uri}/Activity/{activity_id}"

    act_elem = Element(f"{{{PROV_NS}}}Activity")
    act_elem.set(f"{{{RDF_NS}}}about", activity_identity)

    did = SubElement(act_elem, f"{{{SBOL3_NS}}}displayId")
    did.text = activity_id

    started = SubElement(act_elem, f"{{{PROV_NS}}}startedAtTime")
    started.text = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    # Associate the plan (BioCompiler version)
    plan = SubElement(act_elem, f"{{{PROV_NS}}}hadPlan")
    plan.set(f"{{{RDF_NS}}}resource", f"{base_uri}/Plan/biocompiler_{__version__}")

    if organism:
        usage = SubElement(act_elem, f"{{{PROV_NS}}}used")
        usage.set(
            f"{{{RDF_NS}}}resource",
            f"http://identifiers.org/taxonomy/{organism.replace('_', ' ')}",
        )

    return act_elem


def _build_plan_element(base_uri: str) -> Element:
    """Build a PROV-O Plan element for the BioCompiler version."""
    plan_identity = f"{base_uri}/Plan/biocompiler_{__version__}"

    plan_elem = Element(f"{{{PROV_NS}}}Plan")
    plan_elem.set(f"{{{RDF_NS}}}about", plan_identity)

    did = SubElement(plan_elem, f"{{{SBOL3_NS}}}displayId")
    did.text = f"biocompiler_{__version__}"

    label = SubElement(plan_elem, f"{{{DCT_NS}}}title")
    label.text = f"BioCompiler v{__version__} codon optimization"

    return plan_elem


def _optimization_result_to_components(
    result: Any,
    gene_name: str = "",
) -> list[SBOLComponent]:
    """Convert an OptimizationResult into a list of SBOLComponent objects.

    Creates:
    - A top-level gene component
    - A CDS sub-component
    - Optionally promoter and terminator sub-components if present in metadata
    """
    from ..optimization import OptimizationResult

    if not isinstance(result, OptimizationResult):
        raise TypeError(
            f"Expected OptimizationResult, got {type(result).__name__}. "
            f"Use optimize_sequence() to produce an OptimizationResult."
        )

    components: list[SBOLComponent] = []
    display_id = gene_name or "optimized_gene"

    # Main gene component
    gene_comp = SBOLComponent(
        identity="",  # will be assigned by export_sbol
        display_id=display_id,
        component_type="DNA",
        sequence=result.sequence,
        roles=["gene"],
        description=f"BioCompiler-optimized gene: {display_id}",
    )
    components.append(gene_comp)

    # CDS sub-component
    cds_comp = SBOLComponent(
        identity="",
        display_id=f"{display_id}_CDS",
        component_type="DNA",
        sequence=result.sequence,
        roles=["CDS"],
        description=f"Coding sequence for {display_id} (CAI={result.cai:.4f}, GC={result.gc_content:.4f})",
    )
    components.append(cds_comp)

    return components


def _build_rdf_root(base_uri: str) -> Element:
    """Build the root RDF element.

    Namespace prefixes are registered via ``register_namespace`` at module
    load time, so ElementTree will emit human-readable prefixes (sbol3:,
    rdf:, prov:, etc.) automatically when serialising elements in those
    namespaces.
    """
    root = Element(f"{{{RDF_NS}}}RDF")
    return root


# ─── Public API ────────────────────────────────────────────────────


def export_sbol(
    optimization_result: Any,
    output_path: str,
    format: str = "sbol3",
    base_uri: str = "https://biocompiler.org/sbol3",
    gene_name: str = "",
    organism: str = "",
) -> str:
    """Export an optimization result as SBOL3 XML/JSON.

    Converts a BioCompiler :class:`~biocompiler.optimization.OptimizationResult`
    into an SBOL3 document containing:

    - **Component** objects for the gene and CDS, with SO role annotations
    - **Sequence** objects with the optimized DNA
    - **Measure** objects for CAI and GC content (with OM units)
    - **Activity** and **Plan** elements for provenance (PROV-O)

    The output is valid SBOL3 RDF/XML and can be loaded into SynBioHub,
    validated with ``sbol-validate``, or used with any SBOL3 library.

    Args:
        optimization_result: An :class:`~biocompiler.optimization.OptimizationResult`
            from :func:`~biocompiler.optimization.optimize_sequence`.
        output_path: File path to write the SBOL document.
        format: Output format — ``"sbol3"`` (XML, default) or ``"sbol3json"`` (JSON-LD).
        base_uri: Base URI for SBOL identity URIs. Defaults to
            ``"https://biocompiler.org/sbol3"``.
        gene_name: Optional gene name for display_id and annotations.
        organism: Target organism name (used for provenance and taxonomy annotation).

    Returns:
        The absolute path to the written file.

    Raises:
        TypeError: If *optimization_result* is not an OptimizationResult.
        ValueError: If *format* is not ``"sbol3"`` or ``"sbol3json"``.

    Example::

        from biocompiler.optimization import optimize_sequence
        from biocompiler.sbol_export import export_sbol

        result = optimize_sequence('MSKGEELFTG', organism='Escherichia_coli')
        path = export_sbol(result, 'gfp_sbol3.xml', gene_name='gfp',
                           organism='Escherichia_coli')
    """
    from ..optimization import OptimizationResult

    if not isinstance(optimization_result, OptimizationResult):
        raise TypeError(
            f"Expected OptimizationResult, got {type(optimization_result).__name__}. "
            f"Use optimize_sequence() to produce an OptimizationResult."
        )

    if format not in ("sbol3", "sbol3json"):
        raise ValueError(
            f"Unsupported SBOL format: {format!r}. Use 'sbol3' (XML) or 'sbol3json' (JSON-LD)."
        )

    # Determine organism from result metadata if not provided
    effective_organism = organism
    if not effective_organism:
        # Try to extract from provenance
        if optimization_result.provenance is not None:
            effective_organism = getattr(
                optimization_result.provenance, "organism", ""
            ) or ""

    # Generate provenance activity ID
    activity_id = f"optimization_{uuid.uuid4().hex[:8]}"

    # Convert OptimizationResult to SBOLComponent list
    components = _optimization_result_to_components(optimization_result, gene_name)

    # Assign identities
    for comp in components:
        if not comp.identity:
            comp.identity = _make_identity(base_uri, comp.display_id)

    # Build RDF document
    root = _build_rdf_root(base_uri)

    # Add Component and Sequence elements
    for comp in components:
        # Determine CAI/GC measures: attach to the gene component only
        cai_val = None
        gc_val = None
        if comp.roles and "gene" in comp.roles:
            cai_val = optimization_result.cai
            gc_val = optimization_result.gc_content

        comp_elem = _build_component_element(
            comp,
            base_uri,
            cai=cai_val,
            gc=gc_val,
            organism=effective_organism,
            provenance_id=activity_id,
        )
        root.append(comp_elem)

        # Sequence element
        seq_elem = _build_sequence_element(comp, base_uri)
        if seq_elem is not None:
            root.append(seq_elem)

    # Add Activity (provenance)
    act_elem = _build_activity_element(base_uri, activity_id, effective_organism)
    root.append(act_elem)

    # Add Plan
    plan_elem = _build_plan_element(base_uri)
    root.append(plan_elem)

    # Write output
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    if format == "sbol3":
        indent(root, space="  ")
        tree = ElementTree(root)
        tree.write(str(output), xml_declaration=True, encoding="utf-8")
    else:
        # JSON-LD format
        json_doc = _rdf_to_jsonld(root, base_uri)
        with open(str(output), "w", encoding="utf-8") as f:
            json.dump(json_doc, f, indent=2, sort_keys=True)

    logger.info("Exported SBOL3 document to %s (%d components)", output, len(components))
    return str(output.resolve())


def export_sbol_collection(
    results: list[Any],
    output_path: str,
    base_uri: str = "https://biocompiler.org/sbol3",
    collection_name: str = "BioCompiler_optimization_collection",
    organism: str = "",
) -> str:
    """Export multiple optimization results as an SBOL3 Collection.

    An SBOL3 Collection groups multiple Components together, useful for
    representing a library of optimized genes or a multi-gene construct.

    Each result is converted to SBOL Components (gene + CDS) with Measures
    for CAI and GC, and the Collection element links them all together.

    Args:
        results: List of :class:`~biocompiler.optimization.OptimizationResult` objects.
        output_path: File path to write the SBOL document.
        base_uri: Base URI for SBOL identity URIs.
        collection_name: Display ID for the Collection object.
        organism: Default target organism name.

    Returns:
        The absolute path to the written file.

    Example::

        from biocompiler.sbol_export import export_sbol_collection

        results = [optimize_sequence(p, organism='E_coli') for p in proteins]
        path = export_sbol_collection(results, 'gene_library.xml',
                                      collection_name='insulin_variants')
    """
    from ..optimization import OptimizationResult

    # Build RDF document
    root = _build_rdf_root(base_uri)

    # Collection element
    collection_identity = _make_identity(base_uri, collection_name)
    collection_elem = Element(f"{{{SBOL3_NS}}}Collection")
    collection_elem.set(f"{{{RDF_NS}}}about", collection_identity)

    did = SubElement(collection_elem, f"{{{SBOL3_NS}}}displayId")
    did.text = collection_name

    desc = SubElement(collection_elem, f"{{{DCT_NS}}}description")
    desc.text = (
        f"Collection of {len(results)} BioCompiler-optimized gene(s). "
        f"Generated by BioCompiler v{__version__}."
    )

    activity_id = f"batch_{uuid.uuid4().hex[:8]}"

    all_component_identities: list[str] = []

    for idx, result in enumerate(results):
        if not isinstance(result, OptimizationResult):
            logger.warning(
                "Skipping result at index %d: not an OptimizationResult", idx
            )
            continue

        gene_name = getattr(result, "gene_name", "") or f"gene_{idx + 1}"
        eff_organism = organism
        if not eff_organism and result.provenance is not None:
            eff_organism = getattr(result.provenance, "organism", "") or ""

        components = _optimization_result_to_components(result, gene_name)

        for comp in components:
            if not comp.identity:
                comp.identity = _make_identity(base_uri, comp.display_id)

            cai_val = None
            gc_val = None
            if comp.roles and "gene" in comp.roles:
                cai_val = result.cai
                gc_val = result.gc_content

            comp_elem = _build_component_element(
                comp,
                base_uri,
                cai=cai_val,
                gc=gc_val,
                organism=eff_organism,
                provenance_id=activity_id,
            )
            root.append(comp_elem)

            seq_elem = _build_sequence_element(comp, base_uri)
            if seq_elem is not None:
                root.append(seq_elem)

            # Only add gene-level components to the collection
            if comp.roles and "gene" in comp.roles:
                all_component_identities.append(comp.identity)

    # Add member references to collection
    for member_uri in all_component_identities:
        member_elem = SubElement(collection_elem, f"{{{SBOL3_NS}}}member")
        member_elem.set(f"{{{RDF_NS}}}resource", member_uri)

    # Add Activity and Plan
    act_elem = _build_activity_element(base_uri, activity_id, organism)
    root.append(act_elem)
    plan_elem = _build_plan_element(base_uri)
    root.append(plan_elem)

    # Collection must come first in the document
    root.insert(0, collection_elem)

    # Write output
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    indent(root, space="  ")
    tree = ElementTree(root)
    tree.write(str(output), xml_declaration=True, encoding="utf-8")

    logger.info(
        "Exported SBOL3 Collection to %s (%d members)",
        output, len(all_component_identities),
    )
    return str(output.resolve())


def _rdf_to_jsonld(root: Element, base_uri: str) -> dict[str, Any]:
    """Convert the RDF/XML Element tree to a JSON-LD document.

    This is a simplified conversion for the ``sbol3json`` format.
    For production use, a proper RDF library should be used.
    """
    components_data: list[dict[str, Any]] = []
    sequences_data: list[dict[str, Any]] = []
    activities_data: list[dict[str, Any]] = []

    for child in root:
        tag = child.tag
        about = child.get(f"{{{RDF_NS}}}about", "")

        if "Component" in tag:
            comp_dict: dict[str, Any] = {
                "@type": "sbol3:Component",
                "@id": about,
            }
            for sub in child:
                sub_tag = sub.tag.split("}")[-1] if "}" in sub.tag else sub.tag
                if sub_tag == "displayId":
                    comp_dict["sbol3:displayId"] = sub.text or ""
                elif sub_tag == "type":
                    comp_dict["sbol3:type"] = sub.get(f"{{{RDF_NS}}}resource", "")
                elif sub_tag == "role":
                    comp_dict.setdefault("sbol3:role", []).append(
                        sub.get(f"{{{RDF_NS}}}resource", "")
                    )
                elif sub_tag == "description":
                    comp_dict["dct:description"] = sub.text or ""
                elif sub_tag == "hasSequence":
                    comp_dict["sbol3:hasSequence"] = sub.get(
                        f"{{{RDF_NS}}}resource", ""
                    )
                elif sub_tag == "hasMeasure":
                    measures = comp_dict.setdefault("sbol3:hasMeasure", [])
                    measure_data: dict[str, Any] = {}
                    for msub in sub:
                        mtag = msub.tag.split("}")[-1] if "}" in msub.tag else msub.tag
                        if mtag == "Measure":
                            for msub2 in msub:
                                mtag2 = msub2.tag.split("}")[-1] if "}" in msub2.tag else msub2.tag
                                if mtag2 == "displayId":
                                    measure_data["displayId"] = msub2.text or ""
                                elif mtag2 == "value":
                                    measure_data["value"] = float(msub2.text or "0")
                                elif mtag2 == "unit":
                                    measure_data["unit"] = msub2.get(
                                        f"{{{RDF_NS}}}resource", ""
                                    )
                                elif mtag2 == "title":
                                    measure_data["title"] = msub2.text or ""
                    if measure_data:
                        measures.append(measure_data)
                elif sub_tag == "creator":
                    comp_dict["dct:creator"] = sub.text or ""
            components_data.append(comp_dict)

        elif "Sequence" in tag:
            seq_dict: dict[str, Any] = {
                "@type": "sbol3:Sequence",
                "@id": about,
            }
            for sub in child:
                sub_tag = sub.tag.split("}")[-1] if "}" in sub.tag else sub.tag
                if sub_tag == "displayId":
                    seq_dict["sbol3:displayId"] = sub.text or ""
                elif sub_tag == "elements":
                    seq_dict["sbol3:elements"] = sub.text or ""
                elif sub_tag == "encoding":
                    seq_dict["sbol3:encoding"] = sub.get(
                        f"{{{RDF_NS}}}resource", ""
                    )
            sequences_data.append(seq_dict)

        elif "Activity" in tag:
            act_dict: dict[str, Any] = {
                "@type": "prov:Activity",
                "@id": about,
            }
            for sub in child:
                sub_tag = sub.tag.split("}")[-1] if "}" in sub.tag else sub.tag
                if sub_tag == "displayId":
                    act_dict["sbol3:displayId"] = sub.text or ""
                elif sub_tag == "startedAtTime":
                    act_dict["prov:startedAtTime"] = sub.text or ""
            activities_data.append(act_dict)

    return {
        "@context": {
            "sbol3": SBOL3_NS,
            "prov": PROV_NS,
            "so": SO_NS,
            "om": OM_NS,
            "dct": DCT_NS,
            "rdf": RDF_NS,
        },
        "@id": base_uri,
        "components": components_data,
        "sequences": sequences_data,
        "activities": activities_data,
    }
