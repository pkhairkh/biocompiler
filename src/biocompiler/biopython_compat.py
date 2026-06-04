"""
BioCompiler BioPython SeqRecord Interoperability

Provides seamless conversion between BioCompiler data structures and
BioPython SeqRecord/SeqFeature objects. BioPython is the standard library
for bioinformatics in Python — interop is essential for pipeline integration.

BioPython is OPTIONAL. All functions raise ImportError with a helpful message
if BioPython is not installed.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)

__all__ = [
    "to_seqrecord",
    "from_seqrecord",
    "optimize_to_seqrecord",
]


def _check_biopython() -> None:
    """Check that BioPython is installed, raise ImportError with helpful message if not."""
    try:
        import Bio  # noqa: F401
    except ImportError:
        raise ImportError(
            "BioPython is required for BioCompiler BioPython interop but is not installed. "
            "Install it with: pip install biopython>=1.80  "
            "or: pip install biocompiler[biopython]"
        )


def to_seqrecord(
    sequence: str,
    organism: str = "Homo_sapiens",
    gene_name: Optional[str] = None,
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
    type_results: Optional[list] = None,
    certificate: Optional[object] = None,
) -> "Bio.SeqRecord.SeqRecord":
    """
    Convert BioCompiler result data to a BioPython SeqRecord with features.

    Creates a fully annotated SeqRecord with:
    - Sequence set from the DNA string
    - Organism in annotations
    - Exon boundaries as SeqFeature(type="exon")
    - Type-check results as SeqFeature(type="misc_feature") with notes
    - Certificate embedded in record.annotations

    Args:
        sequence: DNA sequence string
        organism: Organism name (e.g. "Homo_sapiens")
        gene_name: Optional gene name for the record
        exon_boundaries: Optional list of (start, end) tuples for exon features
            (0-based, half-open intervals as used internally)
        type_results: Optional list of TypeCheckResult objects
        certificate: Optional Certificate object

    Returns:
        Bio.SeqRecord.SeqRecord with features and annotations

    Raises:
        ImportError: If BioPython is not installed
    """
    _check_biopython()

    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    seq_upper = sequence.upper()

    record = SeqRecord(
        Seq(seq_upper),
        id=gene_name or "BioCompiler_design",
        name=gene_name or "BioCompiler_design",
        description=f"BioCompiler designed sequence for {organism}",
    )

    # Annotations
    record.annotations["organism"] = organism
    record.annotations["topology"] = "linear"
    record.annotations["molecule_type"] = "DNA"

    from .scanner import gc_content as _gc_content
    record.annotations["gc_content"] = _gc_content(seq_upper)

    # Gene feature
    if gene_name:
        gene_feature = SeqFeature(
            FeatureLocation(0, len(seq_upper)),
            type="gene",
            qualifiers={"gene": [gene_name], "note": ["Designed by BioCompiler"]},
        )
        record.features.append(gene_feature)

    # CDS feature with exon join if boundaries provided
    if exon_boundaries and len(exon_boundaries) > 1:
        # Multi-exon CDS: use join
        from Bio.SeqFeature import CompoundLocation
        locations = [
            FeatureLocation(start, end)
            for start, end in exon_boundaries
        ]
        cds_location = CompoundLocation(locations)
    elif exon_boundaries and len(exon_boundaries) == 1:
        cds_location = FeatureLocation(exon_boundaries[0][0], exon_boundaries[0][1])
    else:
        cds_location = FeatureLocation(0, len(seq_upper))

    from .translation import translate as _translate
    protein = _translate(seq_upper)

    cds_qualifiers = {"note": ["Designed by BioCompiler"], "codon_start": ["1"], "transl_table": ["1"]}
    if gene_name:
        cds_qualifiers["gene"] = [gene_name]
    if protein:
        cds_qualifiers["translation"] = [protein]

    cds_feature = SeqFeature(
        cds_location,
        type="CDS",
        qualifiers=cds_qualifiers,
    )
    record.features.append(cds_feature)

    # Exon features
    if exon_boundaries:
        for i, (start, end) in enumerate(exon_boundaries):
            exon_qualifiers = {"number": [str(i + 1)]}
            if gene_name:
                exon_qualifiers["gene"] = [gene_name]
            exon_feature = SeqFeature(
                FeatureLocation(start, end),
                type="exon",
                qualifiers=exon_qualifiers,
            )
            record.features.append(exon_feature)

    # Type-check results as misc_feature annotations
    if type_results:
        for result in type_results:
            verdict = result.verdict.value if hasattr(result.verdict, "value") else str(result.verdict)
            notes = [
                f"BioCompiler type-check: {result.predicate}={verdict}",
            ]
            if result.violation:
                notes.append(f"Violation: {result.violation}")
            if result.knowledge_gap:
                notes.append(f"Knowledge gap: {result.knowledge_gap}")

            misc_feature = SeqFeature(
                FeatureLocation(0, len(seq_upper)),
                type="misc_feature",
                qualifiers={"note": notes},
            )
            record.features.append(misc_feature)

    # Embed certificate in annotations
    if certificate is not None:
        if hasattr(certificate, "to_dict"):
            cert_dict = certificate.to_dict()
        elif isinstance(certificate, dict):
            cert_dict = certificate
        else:
            cert_dict = {"error": "Certificate could not be serialized"}

        # Store certificate as structured annotation
        record.annotations["biocompiler_certificate"] = cert_dict
        record.annotations["biocompiler_design_id"] = cert_dict.get("design_id", "unknown")

    return record


def from_seqrecord(record: "Bio.SeqRecord.SeqRecord") -> dict:
    """
    Import a BioPython SeqRecord into BioCompiler's internal format.

    Extracts sequence, organism, gene name, exon boundaries, and certificate
    from a SeqRecord that was either created by BioCompiler or annotated
    with compatible features.

    Args:
        record: BioPython SeqRecord object

    Returns:
        Dictionary with keys:
            - sequence: str — DNA sequence
            - organism: str — organism name (from annotations)
            - gene_name: str or None — gene name (from gene/CDS features)
            - exon_boundaries: list[tuple[int, int]] — 0-based half-open exon intervals
            - protein: str or None — protein translation (from CDS feature)
            - certificate: dict or None — embedded certificate data
            - features: list[dict] — all features as simplified dicts

    Raises:
        ImportError: If BioPython is not installed
    """
    _check_biopython()

    result: dict = {
        "sequence": str(record.seq).upper(),
        "organism": record.annotations.get("organism", "Unknown"),
        "gene_name": None,
        "exon_boundaries": [],
        "protein": None,
        "certificate": None,
        "features": [],
    }

    # Extract gene name from gene or CDS features
    for feature in record.features:
        if feature.type == "gene":
            gene_qual = feature.qualifiers.get("gene")
            if gene_qual:
                result["gene_name"] = gene_qual[0] if isinstance(gene_qual, list) else gene_qual
                break

    if result["gene_name"] is None:
        for feature in record.features:
            if feature.type == "CDS":
                gene_qual = feature.qualifiers.get("gene")
                if gene_qual:
                    result["gene_name"] = gene_qual[0] if isinstance(gene_qual, list) else gene_qual
                break

    # Extract protein translation from CDS
    for feature in record.features:
        if feature.type == "CDS":
            translation_qual = feature.qualifiers.get("translation")
            if translation_qual:
                result["protein"] = translation_qual[0] if isinstance(translation_qual, list) else translation_qual
            break

    # Extract exon boundaries from exon features
    exon_features = [f for f in record.features if f.type == "exon"]
    if exon_features:
        # Sort by start position
        exon_features.sort(key=lambda f: int(f.location.start))
        result["exon_boundaries"] = [
            (int(f.location.start), int(f.location.end))
            for f in exon_features
        ]
    else:
        # Try to extract from CDS CompoundLocation (join)
        for feature in record.features:
            if feature.type == "CDS":
                if hasattr(feature.location, "parts"):
                    # CompoundLocation from multi-exon CDS
                    parts = sorted(feature.location.parts, key=lambda p: int(p.start))
                    result["exon_boundaries"] = [
                        (int(part.start), int(part.end))
                        for part in parts
                    ]
                break

    # Extract certificate from annotations
    cert_data = record.annotations.get("biocompiler_certificate")
    if cert_data is not None:
        result["certificate"] = cert_data

    # Collect all features as simplified dicts
    for feature in record.features:
        feat_dict = {
            "type": feature.type,
            "location": f"{feature.location}",
            "qualifiers": {
                k: (v[0] if isinstance(v, list) and len(v) == 1 else v)
                for k, v in feature.qualifiers.items()
            },
        }
        result["features"].append(feat_dict)

    # Compute GC content
    from .scanner import gc_content as _gc_content
    result["gc_content"] = _gc_content(result["sequence"])

    return result


def optimize_to_seqrecord(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.2,
    restriction_enzymes: Optional[list[str]] = None,
    gene_name: Optional[str] = None,
) -> "Bio.SeqRecord.SeqRecord":
    """
    One-shot: optimize a protein sequence and return a BioPython SeqRecord.

    Combines BioCompiler's optimize_sequence with to_seqrecord for
    a convenient pipeline integration point.

    Args:
        protein: Amino acid sequence (single-letter codes)
        organism: Target organism for codon optimization
        gc_lo: Minimum GC content fraction
        gc_hi: Maximum GC content fraction
        cai_threshold: Minimum CAI threshold
        restriction_enzymes: Optional list of restriction enzyme names to avoid
        gene_name: Optional gene name for the SeqRecord

    Returns:
        Bio.SeqRecord.SeqRecord with optimized sequence and features

    Raises:
        ImportError: If BioPython is not installed
        UnsupportedOrganismError: If organism is not supported
        InvalidProteinError: If protein contains invalid amino acids
    """
    from .optimization import optimize_sequence
    from .type_system import evaluate_all_predicates

    # Run optimization
    result = optimize_sequence(
        target_protein=protein,
        organism=organism,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        cai_threshold=cai_threshold,
        restriction_sites=restriction_enzymes,
    )

    # Run type-check on the optimized sequence
    type_results = evaluate_all_predicates(
        seq=result.sequence,
        known_exon_boundaries=[(0, len(result.sequence))],
        organism=organism,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        cai_threshold=cai_threshold,
        enzymes=restriction_enzymes or [],
    )

    # Convert to SeqRecord
    return to_seqrecord(
        sequence=result.sequence,
        organism=organism,
        gene_name=gene_name,
        exon_boundaries=[(0, len(result.sequence))],
        type_results=type_results,
    )
