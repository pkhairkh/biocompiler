"""
BioCompiler SBOL3 Import — Synthetic Biology Open Language v3
=============================================================

Import SBOL3 documents and convert them to BioCompiler data structures
for re-optimization or further analysis.

Supports:
- SBOL3 RDF/XML format (the standard serialization)
- SBOL3 JSON-LD format (alternative serialization)
- Component and Sequence extraction
- Measure extraction (CAI, GC content)
- Role mapping from SO URIs to BioCompiler role names
- Conversion to GeneSpec objects for re-optimization

This module provides a pure-Python SBOL3 parser — no external SBOL
library dependency required.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Optional
from xml.etree.ElementTree import Element

from .sbol_export import (
    SBOL3_NS,
    RDF_NS,
    PROV_NS,
    SO_NS,
    OM_NS,
    DCT_NS,
    SO_CDS,
    SO_PROMOTER,
    SO_TERMINATOR,
    SO_RBS,
    SO_GENE,
    SBOLComponent,
)
from ..exceptions import FileFormatError

logger = logging.getLogger(__name__)

__all__ = [
    "import_sbol",
    "sbol_to_genespecs",
]

# Reverse mapping: SO URI → human-readable role name
_SO_ROLE_MAP: dict[str, str] = {
    SO_CDS: "CDS",
    SO_PROMOTER: "promoter",
    SO_TERMINATOR: "terminator",
    SO_RBS: "RBS",
    SO_GENE: "gene",
    f"{SO_NS}0000347": "DNA",
    f"{SO_NS}0000297": "Protein",
}


# ─── Internal Parsing Helpers ──────────────────────────────────────


def _resolve_input(path_or_text: str) -> str:
    """Resolve input: if path exists, read it; otherwise treat as text content."""
    p = Path(path_or_text)
    if p.exists():
        try:
            return p.read_text(encoding="utf-8")
        except (OSError, IOError) as e:
            raise FileFormatError(path_or_text, "SBOL3", f"Cannot read file: {e}") from e
    return path_or_text


def _detect_format(text: str) -> str:
    """Detect whether the input is XML or JSON."""
    stripped = text.strip()
    if stripped.startswith("{") or stripped.startswith("["):
        return "json"
    if stripped.startswith("<?xml") or stripped.startswith("<rdf:RDF"):
        return "xml"
    # Heuristic: if it starts with <, assume XML
    if stripped.startswith("<"):
        return "xml"
    raise FileFormatError(
        "<input>",
        "SBOL3",
        "Cannot detect format: expected XML or JSON",
    )


def _so_uri_to_role(uri: str) -> str:
    """Map an SO URI to a human-readable role name."""
    return _SO_ROLE_MAP.get(uri, "unknown")


def _parse_xml_components(root: Element) -> list[SBOLComponent]:
    """Parse SBOL3 RDF/XML and extract Component objects."""
    components: list[SBOLComponent] = []
    sequences: dict[str, str] = {}  # identity → sequence string
    component_seqs: dict[str, str] = {}  # component identity → seq identity

    # First pass: collect sequences
    for seq_elem in root.iter(f"{{{SBOL3_NS}}}Sequence"):
        identity = seq_elem.get(f"{{{RDF_NS}}}about", "")
        elements_elem = seq_elem.find(f"{{{SBOL3_NS}}}elements")
        seq_text = elements_elem.text if elements_elem is not None else ""
        if identity and seq_text:
            sequences[identity] = seq_text.strip().upper()

    # Second pass: collect components
    for comp_elem in root.iter(f"{{{SBOL3_NS}}}Component"):
        identity = comp_elem.get(f"{{{RDF_NS}}}about", "")

        # displayId
        did_elem = comp_elem.find(f"{{{SBOL3_NS}}}displayId")
        display_id = did_elem.text if did_elem is not None else ""

        # Component type
        type_elem = comp_elem.find(f"{{{SBOL3_NS}}}type")
        type_uri = type_elem.get(f"{{{RDF_NS}}}resource", "") if type_elem is not None else ""
        comp_type = "Protein" if "Protein" in type_uri or "0297" in type_uri else "DNA"

        # Roles
        roles: list[str] = []
        for role_elem in comp_elem.iter(f"{{{SBOL3_NS}}}role"):
            role_uri = role_elem.get(f"{{{RDF_NS}}}resource", "")
            if role_uri:
                roles.append(_so_uri_to_role(role_uri))

        # Description
        desc_elem = comp_elem.find(f"{{{DCT_NS}}}description")
        description = desc_elem.text if desc_elem is not None else ""

        # Sequence reference
        seq_ref_elem = comp_elem.find(f"{{{SBOL3_NS}}}hasSequence")
        if seq_ref_elem is not None:
            seq_ref = seq_ref_elem.get(f"{{{RDF_NS}}}resource", "")
            if seq_ref:
                component_seqs[identity] = seq_ref

        # Resolve sequence
        seq_identity = component_seqs.get(identity, "")
        sequence = sequences.get(seq_identity, "")

        # If no linked sequence found, try to find by convention
        if not sequence and identity:
            alt_seq_id = f"{identity}/sequence"
            sequence = sequences.get(alt_seq_id, "")

        comp = SBOLComponent(
            identity=identity,
            display_id=display_id,
            component_type=comp_type,
            sequence=sequence,
            roles=roles,
            description=description,
        )
        components.append(comp)

    return components


def _parse_json_components(data: dict[str, Any]) -> list[SBOLComponent]:
    """Parse SBOL3 JSON-LD and extract Component objects."""
    components: list[SBOLComponent] = []
    sequences: dict[str, str] = {}  # identity → sequence

    # Collect sequences
    for seq_data in data.get("sequences", []):
        seq_id = seq_data.get("@id", "")
        seq_elements = seq_data.get("sbol3:elements", "")
        if seq_id and seq_elements:
            sequences[seq_id] = seq_elements.strip().upper()

    # Collect components
    for comp_data in data.get("components", []):
        identity = comp_data.get("@id", "")
        display_id = comp_data.get("sbol3:displayId", "")

        # Type
        type_uri = comp_data.get("sbol3:type", "")
        comp_type = "Protein" if "Protein" in type_uri or "0297" in type_uri else "DNA"

        # Roles
        raw_roles = comp_data.get("sbol3:role", [])
        if isinstance(raw_roles, str):
            raw_roles = [raw_roles]
        roles = [_so_uri_to_role(r) for r in raw_roles]

        # Description
        description = comp_data.get("dct:description", "")

        # Sequence
        seq_ref = comp_data.get("sbol3:hasSequence", "")
        sequence = sequences.get(seq_ref, "")

        comp = SBOLComponent(
            identity=identity,
            display_id=display_id,
            component_type=comp_type,
            sequence=sequence,
            roles=roles,
            description=description,
        )
        components.append(comp)

    return components


def _extract_measures_xml(root: Element) -> dict[str, dict[str, float]]:
    """Extract Measure values (CAI, GC) from XML components.

    Returns a dict mapping component identity → {"CAI": float, "GC_content": float}.
    """
    measures: dict[str, dict[str, float]] = {}

    for comp_elem in root.iter(f"{{{SBOL3_NS}}}Component"):
        comp_identity = comp_elem.get(f"{{{RDF_NS}}}about", "")
        if not comp_identity:
            continue

        comp_measures: dict[str, float] = {}

        for measure_wrapper in comp_elem.iter(f"{{{SBOL3_NS}}}hasMeasure"):
            for measure_elem in measure_wrapper.iter(f"{{{SBOL3_NS}}}Measure"):
                label = ""
                value = 0.0

                title_elem = measure_elem.find(f"{{{DCT_NS}}}title")
                if title_elem is not None and title_elem.text:
                    label = title_elem.text

                value_elem = measure_elem.find(f"{{{SBOL3_NS}}}value")
                if value_elem is not None and value_elem.text:
                    try:
                        value = float(value_elem.text)
                    except ValueError:
                        continue

                if label:
                    comp_measures[label] = value

        if comp_measures:
            measures[comp_identity] = comp_measures

    return measures


def _extract_measures_json(data: dict[str, Any]) -> dict[str, dict[str, float]]:
    """Extract Measure values (CAI, GC) from JSON-LD components."""
    measures: dict[str, dict[str, float]] = {}

    for comp_data in data.get("components", []):
        comp_id = comp_data.get("@id", "")
        if not comp_id:
            continue

        comp_measures: dict[str, float] = {}
        raw_measures = comp_data.get("sbol3:hasMeasure", [])
        if isinstance(raw_measures, dict):
            raw_measures = [raw_measures]

        for m in raw_measures:
            label = m.get("title", "")
            value = m.get("value", 0.0)
            if label:
                try:
                    comp_measures[label] = float(value)
                except (ValueError, TypeError):
                    continue

        if comp_measures:
            measures[comp_id] = comp_measures

    return measures


# ─── Public API ────────────────────────────────────────────────────


def import_sbol(path: str) -> list[SBOLComponent]:
    """Parse an SBOL3 file and extract gene components.

    Reads an SBOL3 document (RDF/XML or JSON-LD format) and extracts
    all :class:`SBOLComponent` objects with their associated sequences,
    roles, and measure annotations.

    The function auto-detects the input format (XML or JSON) and handles
    both file paths and raw text content.

    Args:
        path: Path to an SBOL3 file, or raw SBOL3 text content.

    Returns:
        List of :class:`SBOLComponent` objects extracted from the document.

    Raises:
        FileFormatError: If the input cannot be parsed as SBOL3.

    Example::

        from biocompiler.sbol_import import import_sbol

        components = import_sbol('gfp_sbol3.xml')
        for comp in components:
            print(f"{comp.display_id}: {comp.roles}, {len(comp.sequence)} bp")
    """
    text = _resolve_input(path)
    fmt = _detect_format(text)

    if fmt == "xml":
        from xml.etree.ElementTree import fromstring
        try:
            root = fromstring(text)
        except Exception as e:
            raise FileFormatError(
                path if Path(path).exists() else "<text>",
                "SBOL3",
                f"XML parsing failed: {e}",
            ) from e

        components = _parse_xml_components(root)

    elif fmt == "json":
        try:
            data = json.loads(text)
        except json.JSONDecodeError as e:
            raise FileFormatError(
                path if Path(path).exists() else "<text>",
                "SBOL3",
                f"JSON parsing failed: {e}",
            ) from e

        components = _parse_json_components(data)
    else:
        raise FileFormatError(path, "SBOL3", f"Unsupported format: {fmt}")

    if not components:
        logger.warning("No SBOL3 Components found in input")

    logger.info("Imported %d SBOL3 Component(s)", len(components))
    return components


def sbol_to_genespecs(components: list[SBOLComponent]) -> list[Any]:
    """Convert SBOL components to GeneSpec objects for optimization.

    Filters for DNA components that have CDS or gene roles and a valid
    sequence, then translates the DNA to a protein sequence and creates
    :class:`~biocompiler.multigene.GeneSpec` objects suitable for
    re-optimization with :func:`~biocompiler.multigene.optimize_multigene`
    or :func:`~biocompiler.optimization.optimize_sequence`.

    Only components with:
    - ``component_type == "DNA"``
    - A non-empty sequence (multiple of 3 nucleotides)
    - At least one of: CDS, gene, or promoter role

    are converted.  Other components (e.g., pure promoter or terminator
    sequences without coding content) are skipped.

    Args:
        components: List of :class:`SBOLComponent` objects from
            :func:`import_sbol`.

    Returns:
        List of :class:`~biocompiler.multigene.GeneSpec` objects.

    Example::

        from biocompiler.sbol_import import import_sbol, sbol_to_genespecs

        components = import_sbol('gene_library.xml')
        gene_specs = sbol_to_genespecs(components)
        for spec in gene_specs:
            print(f"Gene: {spec.name}, Protein: {spec.protein[:20]}...")
    """
    from ..multigene import GeneSpec
    from ..translation import translate

    gene_specs: list[GeneSpec] = []

    for comp in components:
        # Filter: must be DNA
        if comp.component_type.upper() != "DNA":
            continue

        # Filter: must have a coding role
        coding_roles = {"CDS", "gene"}
        has_coding_role = bool(set(comp.roles) & coding_roles)
        if not has_coding_role:
            continue

        # Filter: must have a sequence that's a multiple of 3
        seq = comp.sequence.upper().strip()
        if not seq or len(seq) % 3 != 0:
            logger.debug(
                "Skipping component %s: sequence length %d not divisible by 3",
                comp.display_id, len(seq),
            )
            continue

        # Translate DNA to protein
        try:
            protein = translate(seq)
        except Exception as e:
            logger.warning(
                "Cannot translate sequence for %s: %s", comp.display_id, e
            )
            continue

        if not protein:
            logger.debug("Skipping component %s: empty protein translation", comp.display_id)
            continue

        # Extract promoter and terminator from roles
        promoter = ""
        terminator = ""
        for role in comp.roles:
            if role == "promoter":
                # Promoter sequence would need to be extracted from a separate
                # component; for now, leave empty
                pass
            elif role == "terminator":
                pass

        # Create GeneSpec
        try:
            spec = GeneSpec(
                protein=protein,
                name=comp.display_id,
                promoter=promoter,
                terminator=terminator,
            )
            gene_specs.append(spec)
        except Exception as e:
            logger.warning(
                "Cannot create GeneSpec for %s: %s", comp.display_id, e
            )
            continue

    logger.info(
        "Converted %d/%d SBOL components to GeneSpec objects",
        len(gene_specs), len(components),
    )
    return gene_specs
