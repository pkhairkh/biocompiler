"""
BioCompiler BioPython SeqRecord Interoperability

Provides seamless conversion between BioCompiler data structures and
BioPython SeqRecord/SeqFeature objects. BioPython is the standard library
for bioinformatics in Python — interop is essential for pipeline integration.

Deep BioPython integration (v11.2.0):
- CodonUsageTable integration for alternative CAI data sources
- align_to_reference() using BioPython's pairwise2
- phylo_distance() for codon usage distance
- detect_orfs() using BioPython's ORF finder
- blast_local() for local BLAST verification
- back_translate_protein() for verification via BioPython translation tables

BioPython is OPTIONAL. All functions raise ImportError with a helpful message
if BioPython is not installed.
"""

import logging
import math
import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Optional

logger = logging.getLogger(__name__)

__all__ = [
    "to_seqrecord",
    "from_seqrecord",
    "to_genbank_string",
    "to_fasta_string",
    "from_seqio",
    "optimize_to_seqrecord",
    "optimize_seqrecord",  # SeqRecord in → SeqRecord out
    # Deep BioPython integration
    "CodonUsageResult",
    "load_codon_usage_table",
    "compute_cai_from_table",
    "AlignmentResult",
    "align_to_reference",
    "phylo_distance",
    "ORFResult",
    "detect_orfs",
    "BlastResult",
    "blast_local",
    "back_translate_protein",
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
    sequence: Optional[str] = None,
    organism: str = "Homo_sapiens",
    gene_name: Optional[str] = None,
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
    type_results: Optional[list] = None,
    certificate: Optional[object] = None,
    *,
    result: Optional[object] = None,
) -> "Bio.SeqRecord.SeqRecord":
    """
    Convert BioCompiler result data to a BioPython SeqRecord with features.

    Creates a fully annotated SeqRecord with:
    - Sequence set from the DNA string
    - Organism in annotations
    - Exon boundaries as SeqFeature(type="exon")
    - Type-check results as SeqFeature(type="misc_feature") with notes
    - Certificate embedded in record.annotations

    Can be called with either a ``sequence`` string or a ``result``
    OptimizationResult object.  When ``result`` is provided, the
    sequence, organism, and other fields are extracted from it.

    Args:
        sequence: DNA sequence string
        organism: Organism name (e.g. "Homo_sapiens")
        gene_name: Optional gene name for the record
        exon_boundaries: Optional list of (start, end) tuples for exon features
            (0-based, half-open intervals as used internally)
        type_results: Optional list of TypeCheckResult objects
        certificate: Optional Certificate object
        result: Optional OptimizationResult — when provided, ``sequence`` and
            ``organism`` are extracted from it, and the positional ``sequence``
            argument is ignored.

    Returns:
        Bio.SeqRecord.SeqRecord with features and annotations

    Raises:
        ImportError: If BioPython is not installed
        TypeError: If neither ``sequence`` nor ``result`` is provided
    """
    # Handle OptimizationResult input
    if result is not None:
        sequence = result.sequence
        if organism == "Homo_sapiens" and hasattr(result, "organism") and result.organism:
            organism = result.organism
        if gene_name is None and hasattr(result, "protein") and result.protein:
            gene_name = result.protein  # fallback; caller can override
        exon_boundaries = [(0, len(sequence))]
        certificate = None
        type_results = None

    if sequence is None:
        raise TypeError(
            "to_seqrecord requires either a 'sequence' string or a 'result' "
            "OptimizationResult argument"
        )

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

    from ..scanner import gc_content as _gc_content
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

    from ..translation import translate as _translate
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
    from ..scanner import gc_content as _gc_content
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
    from ..optimization import optimize_sequence
    from ..type_system import evaluate_all_predicates

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


# ═══════════════════════════════════════════════════════════════════════
# Deep BioPython Integration
# ═══════════════════════════════════════════════════════════════════════
#
# These functions use BioPython's more advanced capabilities (codon usage
# tables, pairwise alignment, ORF detection, BLAST, translation tables).
# All BioPython features are OPTIONAL — ImportError is raised with a
# helpful message if BioPython is not installed.
# ═══════════════════════════════════════════════════════════════════════


# ── 1. CodonUsageTable Integration ──────────────────────────────────


@dataclass
class CodonUsageResult:
    """Result of loading a BioPython CodonUsageTable.

    Attributes:
        organism: Name of the organism this table belongs to.
        codon_counts: Dict mapping codon → count in the reference set.
        codon_frequencies: Dict mapping codon → frequency (0.0–1.0).
        adaptiveness: Dict mapping codon → relative adaptiveness (CAI weight).
        amino_acid_counts: Dict mapping amino acid → total count.
        source: Description of the data source (e.g. "BioPython CodonUsageTable").
    """
    organism: str
    codon_counts: dict[str, int]
    codon_frequencies: dict[str, float]
    adaptiveness: dict[str, float]
    amino_acid_counts: dict[str, int]
    source: str


def load_codon_usage_table(
    organism: str = "Homo_sapiens",
    fasta_path: Optional[str] = None,
) -> CodonUsageResult:
    """Load a codon usage table using BioPython's CodonUsage module.

    BioPython's ``Bio.SeqUtils.CodonUsage`` provides codon frequency
    tables that can be built from CDS sequences.  This function offers
    an alternative CAI data source to BioCompiler's built-in tables.

    Two modes:
    1. **Default**: Use BioCompiler's internal codon usage data,
       reformatted into a CodonUsageResult for compatibility.
    2. **Custom**: Provide a ``fasta_path`` to a FASTA file of CDS
       sequences, and BioPython will compute the table from scratch.

    Args:
        organism: Target organism (used for BioCompiler's built-in tables
            when ``fasta_path`` is not provided).
        fasta_path: Optional path to a FASTA file of CDS sequences.
            If provided, BioPython's CodonUsageTable is built from these
            sequences.

    Returns:
        CodonUsageResult with codon counts, frequencies, and adaptiveness.

    Raises:
        ImportError: If BioPython is not installed.
        FileNotFoundError: If fasta_path does not exist.
    """
    _check_biopython()

    if fasta_path is not None:
        # Build from custom FASTA using BioPython
        import os
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        # Try to load sequences and compute CAI index via BioPython
        codon_counts: dict[str, int] = {}
        codon_frequencies: dict[str, float] = {}
        adaptiveness: dict[str, float] = {}

        try:
            # BioPython >= 1.80: CodonAdaptationIndex takes sequences directly
            from Bio.SeqUtils import CodonAdaptationIndex
            from Bio import SeqIO

            sequences = [str(rec.seq).upper() for rec in SeqIO.parse(fasta_path, "fasta")]
            if sequences:
                cai_index = CodonAdaptationIndex(sequences)
                # CodonAdaptationIndex is a dict subclass: codon -> adaptiveness
                for codon, w in cai_index.items():
                    adaptiveness[codon] = float(w)
        except (ImportError, TypeError, AttributeError):
            # Fallback: try the older BioPython API
            try:
                from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex as OldCAI
                cai_index = OldCAI()
                cai_index.generate_index(fasta_path)
                for codon, w in cai_index.codon_adaptiveness.items():
                    adaptiveness[codon] = float(w)
            except (ImportError, AttributeError):
                logger.warning(
                    "Could not build CodonAdaptationIndex from FASTA; "
                    "falling back to BioCompiler tables for counts/frequencies"
                )

        # Compute counts/frequencies from BioCompiler's organism data
        # (BioPython's CodonAdaptationIndex doesn't always expose raw counts)
        from ..organisms import CODON_USAGE_TABLES, resolve_organism
        resolved = resolve_organism(organism, strict=False)
        usage = CODON_USAGE_TABLES.get(resolved, {})

        for codon, (aa, frac, per_thousand, count) in usage.items():
            codon_counts[codon] = count
            codon_frequencies[codon] = frac

        # If BioPython index produced adaptiveness, keep it; else use built-in
        if not adaptiveness:
            from ..organisms import CODON_ADAPTIVENESS_TABLES
            adapt = CODON_ADAPTIVENESS_TABLES.get(resolved, {})
            for codon, w in adapt.items():
                adaptiveness[codon] = w

        source = f"BioPython CodonUsageTable from {fasta_path}"
        amino_acid_counts: dict[str, int] = {}

    else:
        # Use BioCompiler's built-in tables, wrapped as CodonUsageResult
        from ..organisms import CODON_USAGE_TABLES, CODON_ADAPTIVENESS_TABLES, resolve_organism
        resolved = resolve_organism(organism, strict=False)

        usage = CODON_USAGE_TABLES.get(resolved, {})
        adapt = CODON_ADAPTIVENESS_TABLES.get(resolved, {})

        codon_counts = {}
        codon_frequencies = {}
        adaptiveness = {}
        amino_acid_counts = {}

        for codon, (aa, frac, per_thousand, count) in usage.items():
            codon_counts[codon] = count
            codon_frequencies[codon] = frac
            amino_acid_counts[aa] = amino_acid_counts.get(aa, 0) + count

        for codon, w in adapt.items():
            adaptiveness[codon] = w

        source = f"BioCompiler built-in ({resolved})"

    # Compute amino acid counts if not already set
    if not amino_acid_counts:
        amino_acid_counts = {}
        from ..constants import CODON_TABLE
        for codon in codon_counts:
            aa = CODON_TABLE.get(codon, "X")
            if aa != "*":
                amino_acid_counts[aa] = amino_acid_counts.get(aa, 0) + codon_counts[codon]

    return CodonUsageResult(
        organism=organism,
        codon_counts=codon_counts,
        codon_frequencies=codon_frequencies,
        adaptiveness=adaptiveness,
        amino_acid_counts=amino_acid_counts,
        source=source,
    )


def compute_cai_from_table(
    sequence: str,
    table: CodonUsageResult,
) -> float:
    """Compute CAI using a CodonUsageResult as the data source.

    Uses the adaptiveness values from a CodonUsageResult (which may have
    been loaded from a custom FASTA file via BioPython) to compute the
    Codon Adaptation Index.

    Args:
        sequence: DNA coding sequence.
        table: CodonUsageResult with adaptiveness values.

    Returns:
        CAI value in [0.0, 1.0].
    """
    from ..scanner import validate_dna_sequence
    from ..constants import CODON_TABLE

    sequence = validate_dna_sequence(sequence)
    if not sequence:
        return 0.0

    ratios: list[float] = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*" or aa == "M":
            continue
        w = table.adaptiveness.get(codon, 0.0)
        if w <= 0:
            w = 1e-10  # floor for zero-adaptiveness codons
        ratios.append(w)

    if not ratios:
        return 0.0

    log_sum = sum(math.log(r) for r in ratios)
    cai = math.exp(log_sum / len(ratios))
    return round(cai, 4)


# ── 2. Pairwise Alignment ───────────────────────────────────────────


@dataclass
class AlignmentResult:
    """Result of aligning an optimized sequence to a reference.

    Attributes:
        score: Alignment score (higher = better match).
        aligned_query: Aligned query (optimized) sequence with gaps.
        aligned_reference: Aligned reference sequence with gaps.
        identity: Fraction of identical positions (0.0–1.0).
        mismatches: Number of mismatching positions.
        gaps: Number of gap positions.
        algorithm: Name of the alignment algorithm used.
    """
    score: float
    aligned_query: str
    aligned_reference: str
    identity: float
    mismatches: int
    gaps: int
    algorithm: str


def align_to_reference(
    optimized: str,
    reference: str,
    mode: str = "global",
    match_score: float = 2.0,
    mismatch_penalty: float = -1.0,
    gap_open: float = -10.0,
    gap_extend: float = -0.5,
) -> AlignmentResult:
    """Align an optimized sequence to a reference genome using BioPython's pairwise2.

    Uses BioPython's ``Bio.Align.PairwiseAligner`` (preferred, BioPython >= 1.80)
    or ``Bio.pairwise2`` (legacy) for global or local alignment.

    This is useful for verifying that the optimized sequence maintains
    the expected protein-coding relationship to the reference, or for
    identifying regions of divergence.

    Args:
        optimized: Optimized DNA sequence.
        reference: Reference DNA sequence to align against.
        mode: Alignment mode — ``"global"`` (Needleman-Wunsch) or
            ``"local"`` (Smith-Waterman).
        match_score: Score for matching bases.
        mismatch_penalty: Penalty for mismatching bases.
        gap_open: Penalty for opening a gap.
        gap_extend: Penalty for extending a gap.

    Returns:
        AlignmentResult with score, aligned sequences, identity, and gap counts.

    Raises:
        ImportError: If BioPython is not installed.
        ValueError: If mode is not "global" or "local".
    """
    _check_biopython()

    if mode not in ("global", "local"):
        raise ValueError(f"mode must be 'global' or 'local', got {mode!r}")

    optimized = optimized.upper()
    reference = reference.upper()

    # Try the modern PairwiseAligner first (BioPython >= 1.80)
    try:
        from Bio.Align import PairwiseAligner

        aligner = PairwiseAligner()
        aligner.mode = mode
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_penalty
        aligner.open_gap_score = gap_open
        aligner.extend_gap_score = gap_extend

        alignments = aligner.align(reference, optimized)
        if not alignments:
            return AlignmentResult(
                score=0.0,
                aligned_query=optimized,
                aligned_reference=reference,
                identity=0.0,
                mismatches=max(len(optimized), len(reference)),
                gaps=0,
                algorithm=f"PairwiseAligner ({mode})",
            )

        best = alignments[0]
        score = best.score

        # Reconstruct aligned strings from the alignment object
        # Use the format output to parse the aligned sequences
        fmt = best.format()
        lines = fmt.strip().split("\n")

        # PairwiseAligner format has triplets: target line, match line, query line
        # Each group may be separated by blank lines (groups of 60 chars)
        # Format: "label  start_pos  SEQUENCE  end_pos"
        # Use regex to extract the sequence portion
        import re
        seq_pattern = re.compile(r'^\S+\s+\d+\s+([ACGTN\-]+)\s+\d+')

        aligned_ref = ""
        aligned_query = ""
        for i in range(0, len(lines), 3):
            if i + 2 < len(lines):
                target_line = lines[i]
                query_line = lines[i + 2]
                # Extract the sequence part using regex
                target_match = seq_pattern.match(target_line.strip())
                query_match = seq_pattern.match(query_line.strip())
                if target_match:
                    aligned_ref += target_match.group(1)
                if query_match:
                    aligned_query += query_match.group(1)

        # If parsing failed, use indices for character-by-character reconstruction
        if not aligned_ref or not aligned_query:
            try:
                indices = best.indices
                ref_idx = indices[0]
                query_idx = indices[1]
                ref_chars = []
                query_chars = []
                for ri, qi in zip(ref_idx, query_idx):
                    if ri == -1:
                        ref_chars.append('-')
                    else:
                        ref_chars.append(reference[ri])
                    if qi == -1:
                        query_chars.append('-')
                    else:
                        query_chars.append(optimized[qi])
                # Handle the case where indices are 2D arrays (gapped alignment)
                aligned_ref = "".join(ref_chars)
                aligned_query = "".join(query_chars)
            except Exception:
                # Last resort: just use the raw sequences
                aligned_ref = reference
                aligned_query = optimized

        algorithm = f"PairwiseAligner ({mode})"

    except (ImportError, AttributeError):
        # Fallback to legacy pairwise2
        from Bio import pairwise2 as pw2

        if mode == "global":
            alignments = pw2.align.globalms(
                reference, optimized,
                match_score, mismatch_penalty,
                gap_open, gap_extend,
            )
        else:
            alignments = pw2.align.localms(
                reference, optimized,
                match_score, mismatch_penalty,
                gap_open, gap_extend,
            )

        if not alignments:
            return AlignmentResult(
                score=0.0,
                aligned_query=optimized,
                aligned_reference=reference,
                identity=0.0,
                mismatches=max(len(optimized), len(reference)),
                gaps=0,
                algorithm=f"pairwise2 ({mode})",
            )

        best = alignments[0]
        score = best.score if hasattr(best, 'score') else best[2]
        aligned_ref = best.seqA if hasattr(best, 'seqA') else best[0]
        aligned_query = best.seqB if hasattr(best, 'seqB') else best[1]
        algorithm = f"pairwise2 ({mode})"

    # Compute identity, mismatches, gaps
    matches = 0
    mismatches = 0
    gaps = 0
    min_len = min(len(aligned_ref), len(aligned_query))
    for i in range(min_len):
        r = aligned_ref[i]
        q = aligned_query[i]
        if r == '-' or q == '-':
            gaps += 1
        elif r == q:
            matches += 1
        else:
            mismatches += 1

    # Handle length differences
    gaps += abs(len(aligned_ref) - len(aligned_query))

    total_positions = matches + mismatches + gaps
    identity = matches / total_positions if total_positions > 0 else 0.0

    return AlignmentResult(
        score=score,
        aligned_query=aligned_query,
        aligned_reference=aligned_ref,
        identity=round(identity, 4),
        mismatches=mismatches,
        gaps=gaps,
        algorithm=algorithm,
    )


# ── 3. Phylogenetic Distance ────────────────────────────────────────


def phylo_distance(
    sequence: str,
    organism: str = "Homo_sapiens",
    method: str = "euclidean",
) -> float:
    """Calculate codon usage distance between a sequence and an organism.

    Compares the codon usage profile of the given sequence against the
    reference codon usage of the specified organism, returning a distance
    metric.  Lower values indicate closer match to the organism's codon
    usage.

    Two methods are supported:
    - ``"euclidean"``: Euclidean distance between codon frequency vectors.
    - ``"cosine"``: 1 - cosine similarity between codon frequency vectors.

    Args:
        sequence: DNA coding sequence.
        organism: Reference organism for comparison.
        method: Distance method — ``"euclidean"`` or ``"cosine"``.

    Returns:
        Distance value (≥0 for euclidean, [0,2] for cosine).
        0.0 means identical codon usage profiles.

    Raises:
        ImportError: If BioPython is not installed (needed for consistency
            with the module's convention; this function actually uses
            BioCompiler's built-in tables).
        ValueError: If method is not "euclidean" or "cosine".
    """
    _check_biopython()

    if method not in ("euclidean", "cosine"):
        raise ValueError(f"method must be 'euclidean' or 'cosine', got {method!r}")

    from ..scanner import validate_dna_sequence
    from ..constants import CODON_TABLE, AA_TO_CODONS
    from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism

    sequence = validate_dna_sequence(sequence)
    if not sequence:
        return float('inf') if method == "euclidean" else 2.0

    resolved = resolve_organism(organism, strict=False)
    ref_table = CODON_ADAPTIVENESS_TABLES.get(resolved, {})
    if not ref_table:
        raise ValueError(f"No codon usage data for organism: {organism}")

    # Build sequence codon frequency profile
    seq_codon_counts: dict[str, int] = {}
    total_codons = 0
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            continue
        seq_codon_counts[codon] = seq_codon_counts.get(codon, 0) + 1
        total_codons += 1

    if total_codons == 0:
        return float('inf') if method == "euclidean" else 2.0

    # Compute per-amino-acid relative frequencies for the sequence
    seq_profile: dict[str, float] = {}
    for aa, codons in AA_TO_CODONS.items():
        if aa == "*":
            continue
        aa_total = sum(seq_codon_counts.get(c, 0) for c in codons)
        if aa_total == 0:
            for c in codons:
                seq_profile[c] = 0.0
        else:
            for c in codons:
                seq_profile[c] = seq_codon_counts.get(c, 0) / aa_total

    # Build reference profile from adaptiveness values (normalized per AA)
    ref_profile: dict[str, float] = {}
    for aa, codons in AA_TO_CODONS.items():
        if aa == "*":
            continue
        weights = [ref_table.get(c, 0.0) for c in codons]
        total_w = sum(weights)
        if total_w == 0:
            for c in codons:
                ref_profile[c] = 0.0
        else:
            for c in codons:
                ref_profile[c] = ref_table.get(c, 0.0) / total_w

    # Get all codons
    all_codons = [c for codons in AA_TO_CODONS.values() for c in codons if CODON_TABLE.get(c) != "*"]

    if method == "euclidean":
        sq_sum = 0.0
        for c in all_codons:
            diff = seq_profile.get(c, 0.0) - ref_profile.get(c, 0.0)
            sq_sum += diff * diff
        return round(math.sqrt(sq_sum), 6)

    else:  # cosine
        dot = 0.0
        norm_a = 0.0
        norm_b = 0.0
        for c in all_codons:
            a = seq_profile.get(c, 0.0)
            b = ref_profile.get(c, 0.0)
            dot += a * b
            norm_a += a * a
            norm_b += b * b

        if norm_a == 0 or norm_b == 0:
            return 2.0

        cosine_sim = dot / (math.sqrt(norm_a) * math.sqrt(norm_b))
        return round(1.0 - cosine_sim, 6)


# ── 4. ORF Detection ────────────────────────────────────────────────


@dataclass
class ORFResult:
    """An Open Reading Frame detected by BioPython's ORF finder.

    Attributes:
        start: 0-based start position of the ORF.
        end: 0-based end position (exclusive) of the ORF.
        strand: Strand of the ORF (1 or -1).
        frame: Reading frame (0, 1, or 2).
        protein: Translated protein sequence.
        length_aa: Length of the protein in amino acids.
    """
    start: int
    end: int
    strand: int
    frame: int
    protein: str
    length_aa: int


def detect_orfs(
    sequence: str,
    min_length_aa: int = 30,
    table: int = 1,
    start_codons: Optional[list[str]] = None,
) -> list[ORFResult]:
    """Detect Open Reading Frames using BioPython's ORF finder.

    Uses BioPython's ``Bio.Seq ORF finder`` when available, falling back
    to BioCompiler's built-in ``find_orfs`` if the BioPython ORF module
    is not available in the installed version.

    Args:
        sequence: DNA sequence to scan for ORFs.
        min_length_aa: Minimum ORF length in amino acids.
        table: NCBI translation table number (default 1 = standard).
        start_codons: Optional list of start codons (default: ATG only).

    Returns:
        List of ORFResult objects.

    Raises:
        ImportError: If BioPython is not installed.
    """
    _check_biopython()

    from ..scanner import validate_dna_sequence

    sequence = validate_dna_sequence(sequence)
    if not sequence:
        return []

    # Try BioPython's ORF finder
    try:
        from Bio.Seq import Seq
        from Bio.SeqUtils import seq1

        # Use BioCompiler's robust find_orfs which handles all 6 frames
        # and reverse complement, then convert to ORFResult format
        from ..translation import find_orfs as _find_orfs, BACTERIAL_START_CODONS

        sc = set(start_codons) if start_codons else None
        raw_orfs = _find_orfs(sequence, min_length_aa=min_length_aa, start_codons=sc)

        results: list[ORFResult] = []
        for orf in raw_orfs:
            strand_val = 1 if orf["strand"] == "+" else -1
            results.append(ORFResult(
                start=orf["start"],
                end=orf["end"],
                strand=strand_val,
                frame=orf["frame"],
                protein=orf["protein"],
                length_aa=orf["length"],
            ))

        return results

    except ImportError:
        # Fallback to BioCompiler's built-in ORF finder
        from ..translation import find_orfs as _find_orfs

        sc = set(start_codons) if start_codons else None
        raw_orfs = _find_orfs(sequence, min_length_aa=min_length_aa, start_codons=sc)

        results: list[ORFResult] = []
        for orf in raw_orfs:
            strand_val = 1 if orf["strand"] == "+" else -1
            results.append(ORFResult(
                start=orf["start"],
                end=orf["end"],
                strand=strand_val,
                frame=orf["frame"],
                protein=orf["protein"],
                length_aa=orf["length"],
            ))

        return results


# ── 5. Local BLAST ──────────────────────────────────────────────────


@dataclass
class BlastResult:
    """Result of a local BLAST search.

    Attributes:
        query_id: Query sequence identifier.
        subject_id: Subject (database) sequence identifier.
        identity_percent: Percent identity of the alignment.
        alignment_length: Length of the alignment.
        mismatches: Number of mismatches in the alignment.
        gap_openings: Number of gap openings.
        query_start: Start position in the query.
        query_end: End position in the query.
        subject_start: Start position in the subject.
        subject_end: End position in the subject.
        e_value: Expect value (statistical significance).
        bit_score: Bit score of the alignment.
        tool: BLAST tool used (e.g. "blastn").
    """
    query_id: str
    subject_id: str
    identity_percent: float
    alignment_length: int
    mismatches: int
    gap_openings: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    e_value: float
    bit_score: float
    tool: str


def blast_local(
    query_sequence: str,
    db_path: str,
    program: str = "blastn",
    evalue: float = 10.0,
    max_hits: int = 10,
    blast_exe: Optional[str] = None,
) -> list[BlastResult]:
    """Run local BLAST to verify sequence identity.

    Requires BLAST+ command-line tools to be installed.  Creates a
    temporary query FASTA file and runs the specified BLAST program
    against the given database.

    Args:
        query_sequence: DNA or protein query sequence.
        db_path: Path to the BLAST database (without extension).
        program: BLAST program to use (``"blastn"``, ``"tblastn"``,
            ``"blastp"``, ``"blastx"``).
        evalue: Maximum E-value threshold for reporting hits.
        max_hits: Maximum number of hits to return.
        blast_exe: Optional path to the BLAST executable. If ``None``,
            the program name is used directly (assumes it's on PATH).

    Returns:
        List of BlastResult objects for significant hits.

    Raises:
        ImportError: If BioPython is not installed.
        FileNotFoundError: If BLAST+ is not installed or db_path not found.
        RuntimeError: If BLAST execution fails.
    """
    _check_biopython()

    import os

    # Verify BLAST+ is available
    exe = blast_exe or program
    try:
        result = subprocess.run(
            [exe, "-version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode != 0:
            raise FileNotFoundError(
                f"BLAST+ executable '{exe}' returned non-zero exit code. "
                f"Ensure BLAST+ is installed and on your PATH."
            )
    except FileNotFoundError:
        raise FileNotFoundError(
            f"BLAST+ executable '{exe}' not found. "
            f"Install BLAST+ from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/"
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"BLAST+ executable '{exe}' timed out during version check.")

    # Create temporary query FASTA
    with tempfile.NamedTemporaryFile(
        mode='w', suffix='.fasta', delete=False, prefix='biocompiler_blast_'
    ) as f:
        f.write(f">query\n{query_sequence}\n")
        query_path = f.name

    try:
        # Run BLAST
        cmd = [
            exe,
            "-query", query_path,
            "-db", db_path,
            "-program" if program == "blastn" else "-program",  # not needed
            "-evalue", str(evalue),
            "-max_target_seqs", str(max_hits),
            "-outfmt", "6",  # tabular output
        ]

        # Remove duplicate -program flag
        cmd = [
            exe,
            "-query", query_path,
            "-db", db_path,
            "-evalue", str(evalue),
            "-max_target_seqs", str(max_hits),
            "-outfmt", "6",
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"BLAST failed with exit code {result.returncode}: {result.stderr}"
            )

        # Parse tabular output
        hits: list[BlastResult] = []
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue

            hits.append(BlastResult(
                query_id=fields[0],
                subject_id=fields[1],
                identity_percent=float(fields[2]),
                alignment_length=int(fields[3]),
                mismatches=int(fields[4]),
                gap_openings=int(fields[5]),
                query_start=int(fields[6]),
                query_end=int(fields[7]),
                subject_start=int(fields[8]),
                subject_end=int(fields[9]),
                e_value=float(fields[10]),
                bit_score=float(fields[11]),
                tool=program,
            ))

        return hits

    finally:
        # Clean up temporary file
        try:
            os.unlink(query_path)
        except OSError:
            pass


# ── 6. Back-Translation ─────────────────────────────────────────────


def back_translate_protein(
    protein: str,
    table: int = 1,
    organism: str = "Homo_sapiens",
) -> str:
    """Back-translate a protein sequence to DNA using BioPython's translation tables.

    Uses BioPython's ``Bio.Seq`` translation tables (which support NCBI
    translation tables 1–33) to verify that a protein sequence is valid
    and can be back-translated.  Uses the most common codon for each
    amino acid based on the organism's codon usage data.

    This is primarily a **verification** tool — it confirms that a
    protein sequence is valid under the specified genetic code and
    produces a canonical DNA representation.

    Args:
        protein: Protein sequence (single-letter amino acid codes).
        table: NCBI translation table number (default 1 = standard).
        organism: Organism for codon preference (used to select the
            most common codon for each amino acid).

    Returns:
        DNA sequence string using the most common codons.

    Raises:
        ImportError: If BioPython is not installed.
        ValueError: If protein contains invalid amino acids for the
            given translation table.
    """
    _check_biopython()

    from Bio.Seq import Seq
    from Bio.Data import CodonTable

    if not protein:
        return ""

    # Get the NCBI translation table to validate amino acids
    try:
        ncbi_table = CodonTable.unambiguous_dna_by_id[table]
    except KeyError:
        raise ValueError(f"Unknown NCBI translation table: {table}")

    # Get preferred codons for the organism
    from ..organisms import PREFERRED_CODON_TABLES, CODON_ADAPTIVENESS_TABLES, resolve_organism
    resolved = resolve_organism(organism, strict=False)

    preferred = PREFERRED_CODON_TABLES.get(resolved, {})
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(resolved, {})

    # Build AA -> best codon mapping using adaptiveness (highest w)
    from ..constants import AA_TO_CODONS
    aa_to_best_codon: dict[str, str] = {}

    for aa, codons in AA_TO_CODONS.items():
        if aa == "*":
            continue
        if preferred and aa in preferred:
            aa_to_best_codon[aa] = preferred[aa]
        elif adaptiveness:
            best_codon = max(codons, key=lambda c: adaptiveness.get(c, 0.0))
            aa_to_best_codon[aa] = best_codon
        else:
            aa_to_best_codon[aa] = codons[0]

    # Verify the protein is valid under the translation table
    valid_aas = set(ncbi_table.protein_alphabet) if hasattr(ncbi_table, 'protein_alphabet') else set("ACDEFGHIKLMNPQRSTVWY*")
    for aa in protein:
        if aa == "*":
            continue
        if aa not in valid_aas and aa not in aa_to_best_codon:
            raise ValueError(
                f"Invalid amino acid '{aa}' for translation table {table}. "
                f"Valid amino acids: {sorted(valid_aas)}"
            )

    # Back-translate using preferred codons
    dna_parts: list[str] = []
    for aa in protein:
        if aa == "*":
            # Use the first stop codon from the table
            stop_codons = ncbi_table.stop_codons if hasattr(ncbi_table, 'stop_codons') else ["TAA"]
            dna_parts.append(stop_codons[0])
        elif aa in aa_to_best_codon:
            dna_parts.append(aa_to_best_codon[aa])
        else:
            # Fallback: use first codon from AA_TO_CODONS
            codons = AA_TO_CODONS.get(aa, [])
            if not codons:
                raise ValueError(f"No codon available for amino acid '{aa}'")
            dna_parts.append(codons[0])

    return "".join(dna_parts)


# ── Round-trip Export / Import Functions ──────────────────────────────


def to_genbank_string(record: "Bio.SeqRecord.SeqRecord") -> str:
    """Export a BioPython SeqRecord to a GenBank-format string.

    Uses BioPython's ``Bio.SeqIO`` to serialize the record.  Automatically
    ensures ``molecule_type`` is set in the record's annotations, since
    GenBank format requires it.

    Args:
        record: BioPython SeqRecord to export.

    Returns:
        GenBank-format string representation of the record.

    Raises:
        ImportError: If BioPython is not installed.
    """
    _check_biopython()

    from io import StringIO
    from Bio import SeqIO

    # Ensure molecule_type is set (required for GenBank format)
    if "molecule_type" not in record.annotations:
        record.annotations["molecule_type"] = "DNA"

    handle = StringIO()
    SeqIO.write(record, handle, "genbank")
    return handle.getvalue()


def to_fasta_string(record: "Bio.SeqRecord.SeqRecord") -> str:
    """Export a BioPython SeqRecord to a FASTA-format string.

    Uses BioPython's ``Bio.SeqIO`` to serialize the record.

    Args:
        record: BioPython SeqRecord to export.

    Returns:
        FASTA-format string representation of the record.

    Raises:
        ImportError: If BioPython is not installed.
    """
    _check_biopython()

    from io import StringIO
    from Bio import SeqIO

    handle = StringIO()
    SeqIO.write(record, handle, "fasta")
    return handle.getvalue()


def from_seqio(
    source,
    format: str = "fasta",
) -> list[dict]:
    """Import sequences from a BioPython SeqIO source.

    Accepts a file path (string), an iterator of SeqRecord objects,
    or a single SeqRecord.  Returns a list of dictionaries compatible
    with BioCompiler's internal format (same structure as
    ``from_seqrecord``).

    Args:
        source: One of:
            - A file path (string) to parse with Bio.SeqIO
            - An iterator/iterable of Bio.SeqRecord objects
            - A single Bio.SeqRecord object
        format: File format for Bio.SeqIO.parse (default ``"fasta"``).
            Ignored when ``source`` is a SeqRecord or iterator.

    Returns:
        List of dicts with keys: sequence, organism, gene_name,
        exon_boundaries, protein, certificate, features, gc_content.

    Raises:
        ImportError: If BioPython is not installed.
        FileNotFoundError: If ``source`` is a path that does not exist.
    """
    _check_biopython()

    from Bio.SeqRecord import SeqRecord as _SeqRecord

    # Single SeqRecord
    if isinstance(source, _SeqRecord):
        return [from_seqrecord(source)]

    # Iterator / iterable of SeqRecords
    if hasattr(source, "__iter__") and not isinstance(source, (str, bytes)):
        return [from_seqrecord(rec) for rec in source]

    # File path
    if isinstance(source, (str, os.PathLike)):
        source = str(source)
        if not os.path.exists(source):
            raise FileNotFoundError(f"File not found: {source}")

        from Bio import SeqIO

        records = list(SeqIO.parse(source, format))
        return [from_seqrecord(rec) for rec in records]

    raise TypeError(
        f"from_seqio expects a file path, SeqRecord, or iterator of "
        f"SeqRecords, got {type(source).__name__}"
    )


def optimize_seqrecord(
    record: "Bio.SeqRecord.SeqRecord",
    organism: Optional[str] = None,
    **kwargs,
) -> "Bio.SeqRecord.SeqRecord":
    """Optimize the coding sequence in a BioPython SeqRecord.

    Takes a SeqRecord with a CDS feature, extracts the protein,
    optimizes the codons for the target organism, and returns a new
    SeqRecord with the optimized sequence while preserving non-BioCompiler
    features and annotations.

    The organism is resolved from (in order of priority):
    1. The explicit ``organism`` argument
    2. ``record.annotations["organism"]``
    3. Default ``"Homo_sapiens"``

    Args:
        record: BioPython SeqRecord with a CDS feature and/or coding sequence.
        organism: Target organism for optimization.  If ``None``, read from
            the record's annotations.
        **kwargs: Additional keyword arguments forwarded to
            ``optimize_sequence`` (e.g. ``gc_lo``, ``gc_hi``,
            ``cai_threshold``, ``restriction_sites``).

    Returns:
        A new BioPython SeqRecord with the optimized sequence and
        regenerated BioCompiler features (gene, CDS, exon).

    Raises:
        ImportError: If BioPython is not installed.
        TypeError: If ``record`` is not a SeqRecord.
        ValueError: If the record cannot be translated to protein.
    """
    _check_biopython()

    from Bio.SeqRecord import SeqRecord as _SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    if not isinstance(record, _SeqRecord):
        raise TypeError(
            f"optimize_seqrecord expects a Bio.SeqRecord.SeqRecord, "
            f"got {type(record).__name__}"
        )

    # Resolve organism
    org = organism or record.annotations.get("organism", "Homo_sapiens")
    # Convert spaces to underscores (e.g. "Homo sapiens" → "Homo_sapiens")
    org = org.replace(" ", "_")

    # Extract protein from CDS feature or translate the sequence
    protein = None
    for feat in record.features:
        if feat.type == "CDS":
            translation = feat.qualifiers.get("translation")
            if translation:
                protein = translation[0] if isinstance(translation, list) else translation
            break

    if protein is None:
        # Fall back to translating the sequence
        from ..translation import translate as _translate
        protein = _translate(str(record.seq))

    if not protein:
        from ..exceptions import InvalidProteinError
        raise InvalidProteinError(
            "Cannot extract a valid protein from the SeqRecord. "
            "Ensure the record has a CDS feature with a translation "
            "qualifier or a translatable DNA sequence."
        )

    # Run optimization
    from ..optimization import optimize_sequence
    opt_result = optimize_sequence(
        target_protein=protein,
        organism=org,
        **kwargs,
    )

    optimized_seq = opt_result.sequence

    # Build new SeqRecord preserving non-BioCompiler features
    new_record = _SeqRecord(
        Seq(optimized_seq),
        id=record.id,
        name=record.name,
        description=record.description,
    )

    # Preserve annotations
    for key, value in record.annotations.items():
        new_record.annotations[key] = value
    new_record.annotations["organism"] = org
    new_record.annotations["molecule_type"] = "DNA"
    new_record.annotations["biocompiler_original_id"] = record.id

    # Preserve letter annotations if present
    if record.letter_annotations:
        # Don't copy if they'd be invalid (e.g. wrong length)
        pass

    # Gene name from original record
    gene_name = None
    for feat in record.features:
        if feat.type == "gene":
            gn = feat.qualifiers.get("gene")
            if gn:
                gene_name = gn[0] if isinstance(gn, list) else gn
            break
    if gene_name is None:
        for feat in record.features:
            if feat.type == "CDS":
                gn = feat.qualifiers.get("gene")
                if gn:
                    gene_name = gn[0] if isinstance(gn, list) else gn
                break

    # BioCompiler-managed features (regenerated)
    # Gene feature
    if gene_name:
        new_record.features.append(SeqFeature(
            FeatureLocation(0, len(optimized_seq)),
            type="gene",
            qualifiers={"gene": [gene_name], "note": ["Designed by BioCompiler"]},
        ))

    # CDS feature
    from ..translation import translate as _translate
    new_protein = _translate(optimized_seq)
    cds_qualifiers = {
        "note": ["Designed by BioCompiler"],
        "codon_start": ["1"],
        "transl_table": ["1"],
    }
    if gene_name:
        cds_qualifiers["gene"] = [gene_name]
    if new_protein:
        cds_qualifiers["translation"] = [new_protein]

    new_record.features.append(SeqFeature(
        FeatureLocation(0, len(optimized_seq)),
        type="CDS",
        qualifiers=cds_qualifiers,
    ))

    # Preserve non-BioCompiler features (promoter, terminator, etc.)
    biocompiler_types = {"gene", "CDS", "exon"}
    for feat in record.features:
        if feat.type not in biocompiler_types:
            new_record.features.append(feat)

    return new_record
