"""
BioCompiler Sequence Annotation Enrichment

Automatically annotate biological features in optimized DNA sequences,
including ORFs, restriction sites, CpG islands, splice sites,
repetitive elements, and GC/AT-rich regions.  Produces fully-annotated
GenBank records suitable for submission to NCBI or import into
Benchling/SnapGene.

Feature types detected:
- CDS: Open reading frames
- promoter: Putative promoter elements (TATA box, GC box)
- RBS: Ribosome binding site (Shine-Dalgarno for prokaryotes)
- restriction_site: Known restriction enzyme recognition sites
- CpG_island: Regions of elevated CpG dinucleotide frequency
- splice_donor: GT donor sites with MaxEntScan scoring
- splice_acceptor: AG acceptor sites with MaxEntScan scoring
- simple_repeat: Dinucleotide/trinucleotide repeats
- GC_rich: Regions with high GC content (>70%)
- AT_rich: Regions with high AT content (>70%)
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass, field
from typing import Optional

from ..scanner import gc_content, validate_dna_sequence
from ..restriction_sites import RESTRICTION_SITES, get_eliminable_sites
from ..constants import CODON_TABLE, AA_TO_CODONS, reverse_complement, STOP_CODONS

logger = logging.getLogger(__name__)

__all__ = [
    "SequenceAnnotation",
    "annotate_sequence",
    "annotate_to_genbank",
    "_find_orfs",
    "_find_restriction_sites",
    "_find_cpg_islands",
    "_find_gc_at_rich_regions",
    "_find_rbs",
    "_find_simple_repeats",
    "_find_splice_sites",
]


@dataclass
class SequenceAnnotation:
    """A biological feature annotation on a DNA sequence.

    Attributes:
        feature_type: Type of feature (e.g. "CDS", "restriction_site",
            "CpG_island", "splice_donor", "splice_acceptor",
            "simple_repeat", "GC_rich", "AT_rich", "promoter", "RBS").
        start: 0-based start position (inclusive).
        end: 0-based end position (exclusive).
        strand: Strand direction: 1 (forward), -1 (reverse), 0 (both/neutral).
        qualifiers: Additional key-value metadata for the feature.
    """
    feature_type: str
    start: int
    end: int
    strand: int
    qualifiers: dict[str, str] = field(default_factory=dict)

    def __post_init__(self):
        if self.feature_type not in (
            "CDS", "promoter", "RBS", "restriction_site", "CpG_island",
            "splice_donor", "splice_acceptor", "simple_repeat",
            "GC_rich", "AT_rich",
        ):
            logger.debug("Non-standard feature_type: %s", self.feature_type)


# ── Annotation parameters ───────────────────────────────────────────

# CpG island detection parameters (per Gardiner-Garden & Frommer 1987)
CPG_WINDOW_SIZE: int = 200
CPG_MIN_GC: float = 0.50
CPG_MIN_OBS_EXP: float = 0.60

# Simple repeat detection
MIN_REPEAT_UNIT_LEN: int = 2
MAX_REPEAT_UNIT_LEN: int = 6
MIN_REPEAT_COPIES: int = 3

# GC/AT-rich region detection
GC_RICH_WINDOW: int = 50
GC_RICH_THRESHOLD: float = 0.70
AT_RICH_WINDOW: int = 50
AT_RICH_THRESHOLD: float = 0.70

# RBS detection (Shine-Dalgarno consensus)
SHINE_DALGARNO_PATTERNS: list[str] = ["AGGAGG", "AGGAG", "GGAGG", "AGGA", "GGAG"]
RBS_UPSTREAM_DIST: int = 20  # search window upstream of ATG


def _find_orfs(dna: str) -> list[SequenceAnnotation]:
    """Find open reading frames in all 6 frames."""
    from ..translation import find_orfs as _find_orfs_bio

    orfs = _find_orfs_bio(dna, min_length_aa=30)
    annotations: list[SequenceAnnotation] = []

    for orf in orfs:
        strand_val = 1 if orf["strand"] == "+" else -1
        protein = orf["protein"]
        qualifiers: dict[str, str] = {
            "translation": protein,
            "codon_start": "1",
            "transl_table": "1",
            "length_aa": str(orf["length"]),
            "frame": str(orf["frame"]),
        }
        annotations.append(SequenceAnnotation(
            feature_type="CDS",
            start=orf["start"],
            end=orf["end"],
            strand=strand_val,
            qualifiers=qualifiers,
        ))

    return annotations


def _find_restriction_sites(dna: str) -> list[SequenceAnnotation]:
    """Find restriction enzyme recognition sites in the sequence."""
    annotations: list[SequenceAnnotation] = []
    dna_upper = dna.upper()
    sites = get_eliminable_sites(min_length=4)  # include 4+ cutters for annotation

    for enzyme_name, site_seq in sites.items():
        if any(b not in "ACGT" for b in site_seq):
            # Skip IUPAC sites for simple annotation
            continue
        # Search forward strand
        start = 0
        while True:
            pos = dna_upper.find(site_seq, start)
            if pos == -1:
                break
            annotations.append(SequenceAnnotation(
                feature_type="restriction_site",
                start=pos,
                end=pos + len(site_seq),
                strand=1,
                qualifiers={
                    "enzyme": enzyme_name,
                    "site": site_seq,
                    "label": enzyme_name,
                },
            ))
            start = pos + 1

        # Search reverse complement
        rc_site = reverse_complement(site_seq)
        if rc_site != site_seq:  # skip palindromic sites (already found)
            start = 0
            while True:
                pos = dna_upper.find(rc_site, start)
                if pos == -1:
                    break
                annotations.append(SequenceAnnotation(
                    feature_type="restriction_site",
                    start=pos,
                    end=pos + len(rc_site),
                    strand=-1,
                    qualifiers={
                        "enzyme": enzyme_name,
                        "site": rc_site,
                        "label": enzyme_name,
                    },
                ))
                start = pos + 1

    return annotations


def _find_cpg_islands(dna: str) -> list[SequenceAnnotation]:
    """Find CpG islands using sliding window (Gardiner-Garden & Frommer 1987).

    A CpG island is defined as a region of at least 200 bp with:
    - GC content ≥ 50%
    - Obs/Exp CpG ratio ≥ 0.60
    """
    annotations: list[SequenceAnnotation] = []
    dna_upper = dna.upper()
    seq_len = len(dna_upper)

    if seq_len < CPG_WINDOW_SIZE:
        return annotations

    in_island = False
    island_start = 0

    for i in range(seq_len - CPG_WINDOW_SIZE + 1):
        window = dna_upper[i:i + CPG_WINDOW_SIZE]

        # GC content
        gc = (window.count('G') + window.count('C')) / len(window)

        # Obs/Exp CpG ratio
        c_count = window.count('C')
        g_count = window.count('G')
        cpg_count = 0
        for j in range(len(window) - 1):
            if window[j] == 'C' and window[j + 1] == 'G':
                cpg_count += 1

        expected_cpg = (c_count * g_count) / len(window) if len(window) > 0 else 0
        obs_exp = cpg_count / expected_cpg if expected_cpg > 0 else 0.0

        if gc >= CPG_MIN_GC and obs_exp >= CPG_MIN_OBS_EXP:
            if not in_island:
                island_start = i
                in_island = True
        else:
            if in_island:
                end = i + CPG_WINDOW_SIZE - 1
                length = end - island_start
                if length >= CPG_WINDOW_SIZE:
                    island_gc = gc_content(dna_upper[island_start:end])
                    annotations.append(SequenceAnnotation(
                        feature_type="CpG_island",
                        start=island_start,
                        end=end,
                        strand=0,
                        qualifiers={
                            "gc_content": f"{island_gc:.4f}",
                            "length": str(length),
                            "method": "Gardiner-Garden_1987",
                        },
                    ))
                in_island = False

    # Handle island that extends to end of sequence
    if in_island:
        end = seq_len
        length = end - island_start
        if length >= CPG_WINDOW_SIZE:
            island_gc = gc_content(dna_upper[island_start:end])
            annotations.append(SequenceAnnotation(
                feature_type="CpG_island",
                start=island_start,
                end=end,
                strand=0,
                qualifiers={
                    "gc_content": f"{island_gc:.4f}",
                    "length": str(length),
                    "method": "Gardiner-Garden_1987",
                },
            ))

    return annotations


def _find_splice_sites(dna: str) -> list[SequenceAnnotation]:
    """Find splice donor (GT) and acceptor (AG) sites with MaxEntScan scoring."""
    annotations: list[SequenceAnnotation] = []
    dna_upper = dna.upper()

    try:
        from ..maxentscan import score_donor, score_acceptor
        has_maxent = True
    except ImportError:
        has_maxent = False

    # Find donor sites (GT dinucleotides)
    for i in range(len(dna_upper) - 1):
        if dna_upper[i:i + 2] == "GT":
            if has_maxent:
                try:
                    score = score_donor(dna_upper, i)
                except Exception:
                    score = 0.0
            else:
                score = 0.0

            annotations.append(SequenceAnnotation(
                feature_type="splice_donor",
                start=i,
                end=i + 2,
                strand=1,
                qualifiers={
                    "consensus": "GT",
                    "score": f"{score:.2f}",
                },
            ))

    # Find acceptor sites (AG dinucleotides)
    for i in range(len(dna_upper) - 1):
        if dna_upper[i:i + 2] == "AG":
            if has_maxent:
                try:
                    score = score_acceptor(dna_upper, i)
                except Exception:
                    score = 0.0
            else:
                score = 0.0

            annotations.append(SequenceAnnotation(
                feature_type="splice_acceptor",
                start=i,
                end=i + 2,
                strand=1,
                qualifiers={
                    "consensus": "AG",
                    "score": f"{score:.2f}",
                },
            ))

    return annotations


def _find_simple_repeats(dna: str) -> list[SequenceAnnotation]:
    """Find simple dinucleotide/trinucleotide repeats."""
    annotations: list[SequenceAnnotation] = []
    dna_upper = dna.upper()

    for unit_len in range(MIN_REPEAT_UNIT_LEN, MAX_REPEAT_UNIT_LEN + 1):
        i = 0
        while i <= len(dna_upper) - unit_len * MIN_REPEAT_COPIES:
            unit = dna_upper[i:i + unit_len]
            if not all(b in "ACGT" for b in unit):
                i += 1
                continue

            # Count consecutive repeats
            copies = 1
            j = i + unit_len
            while j + unit_len <= len(dna_upper) and dna_upper[j:j + unit_len] == unit:
                copies += 1
                j += unit_len

            if copies >= MIN_REPEAT_COPIES:
                annotations.append(SequenceAnnotation(
                    feature_type="simple_repeat",
                    start=i,
                    end=j,
                    strand=0,
                    qualifiers={
                        "repeat_unit": unit,
                        "copies": str(copies),
                        "repeat_type": f"{unit_len}-mer",
                    },
                ))
                i = j  # skip past this repeat
            else:
                i += 1

    return annotations


def _find_gc_at_rich_regions(dna: str) -> list[SequenceAnnotation]:
    """Find GC-rich and AT-rich regions using a sliding window."""
    annotations: list[SequenceAnnotation] = []
    dna_upper = dna.upper()
    seq_len = len(dna_upper)

    # GC-rich regions
    if seq_len >= GC_RICH_WINDOW:
        in_gc_rich = False
        gc_start = 0
        for i in range(seq_len - GC_RICH_WINDOW + 1):
            window = dna_upper[i:i + GC_RICH_WINDOW]
            gc = (window.count('G') + window.count('C')) / GC_RICH_WINDOW
            if gc >= GC_RICH_THRESHOLD:
                if not in_gc_rich:
                    gc_start = i
                    in_gc_rich = True
            else:
                if in_gc_rich:
                    end = i + GC_RICH_WINDOW
                    annotations.append(SequenceAnnotation(
                        feature_type="GC_rich",
                        start=gc_start,
                        end=end,
                        strand=0,
                        qualifiers={
                            "gc_content": f"{gc_content(dna_upper[gc_start:end]):.4f}",
                            "window_size": str(GC_RICH_WINDOW),
                        },
                    ))
                    in_gc_rich = False

        if in_gc_rich:
            annotations.append(SequenceAnnotation(
                feature_type="GC_rich",
                start=gc_start,
                end=seq_len,
                strand=0,
                qualifiers={
                    "gc_content": f"{gc_content(dna_upper[gc_start:seq_len]):.4f}",
                    "window_size": str(GC_RICH_WINDOW),
                },
            ))

    # AT-rich regions
    if seq_len >= AT_RICH_WINDOW:
        in_at_rich = False
        at_start = 0
        for i in range(seq_len - AT_RICH_WINDOW + 1):
            window = dna_upper[i:i + AT_RICH_WINDOW]
            at = (window.count('A') + window.count('T')) / AT_RICH_WINDOW
            if at >= AT_RICH_THRESHOLD:
                if not in_at_rich:
                    at_start = i
                    in_at_rich = True
            else:
                if in_at_rich:
                    end = i + AT_RICH_WINDOW
                    annotations.append(SequenceAnnotation(
                        feature_type="AT_rich",
                        start=at_start,
                        end=end,
                        strand=0,
                        qualifiers={
                            "at_content": f"{1.0 - gc_content(dna_upper[at_start:end]):.4f}",
                            "window_size": str(AT_RICH_WINDOW),
                        },
                    ))
                    in_at_rich = False

        if in_at_rich:
            annotations.append(SequenceAnnotation(
                feature_type="AT_rich",
                start=at_start,
                end=seq_len,
                strand=0,
                qualifiers={
                    "at_content": f"{1.0 - gc_content(dna_upper[at_start:seq_len]):.4f}",
                    "window_size": str(AT_RICH_WINDOW),
                },
            ))

    return annotations


def _find_rbs(dna: str) -> list[SequenceAnnotation]:
    """Find putative ribosome binding sites (Shine-Dalgarno sequences)."""
    annotations: list[SequenceAnnotation] = []
    dna_upper = dna.upper()

    for i in range(len(dna_upper) - 2):
        codon = dna_upper[i:i + 3]
        if codon != "ATG":
            continue

        # Search upstream for RBS
        upstream_start = max(0, i - RBS_UPSTREAM_DIST)
        upstream = dna_upper[upstream_start:i]

        for pattern in SHINE_DALGARNO_PATTERNS:
            pos = upstream.rfind(pattern)
            if pos != -1:
                abs_start = upstream_start + pos
                annotations.append(SequenceAnnotation(
                    feature_type="RBS",
                    start=abs_start,
                    end=abs_start + len(pattern),
                    strand=1,
                    qualifiers={
                        "pattern": pattern,
                        "type": "Shine-Dalgarno",
                        "distance_to_atg": str(i - abs_start - len(pattern)),
                    },
                ))
                break  # only report the best (longest) match

    return annotations


def annotate_sequence(dna: str, organism: str = "") -> list[SequenceAnnotation]:
    """Automatically annotate biological features in an optimized sequence.

    Scans the DNA sequence for the following feature types:
    - **CDS**: Open reading frames (≥30 aa) in all 6 frames
    - **restriction_site**: Known restriction enzyme recognition sites
    - **CpG_island**: Regions of elevated CpG dinucleotide frequency
    - **splice_donor**: GT dinucleotide splice donor sites (with MaxEntScan scores)
    - **splice_acceptor**: AG dinucleotide splice acceptor sites
    - **simple_repeat**: Di/tri/tetra/penta/hexanucleotide repeats (≥3 copies)
    - **GC_rich**: Regions with >70% GC content
    - **AT_rich**: Regions with >70% AT content
    - **RBS**: Putative ribosome binding sites (Shine-Dalgarno)

    Args:
        dna: DNA sequence to annotate.
        organism: Optional organism name for context-dependent annotation.

    Returns:
        List of SequenceAnnotation objects, sorted by (start, feature_type).
    """
    dna = validate_dna_sequence(dna)
    if not dna:
        return []

    annotations: list[SequenceAnnotation] = []

    # 1. ORFs
    try:
        annotations.extend(_find_orfs(dna))
    except Exception:
        logger.warning("ORF detection failed", exc_info=True)

    # 2. Restriction sites
    try:
        annotations.extend(_find_restriction_sites(dna))
    except Exception:
        logger.warning("Restriction site detection failed", exc_info=True)

    # 3. CpG islands
    try:
        annotations.extend(_find_cpg_islands(dna))
    except Exception:
        logger.warning("CpG island detection failed", exc_info=True)

    # 4. Splice sites (eukaryotic feature)
    try:
        annotations.extend(_find_splice_sites(dna))
    except Exception:
        logger.warning("Splice site detection failed", exc_info=True)

    # 5. Simple repeats
    try:
        annotations.extend(_find_simple_repeats(dna))
    except Exception:
        logger.warning("Simple repeat detection failed", exc_info=True)

    # 6. GC/AT-rich regions
    try:
        annotations.extend(_find_gc_at_rich_regions(dna))
    except Exception:
        logger.warning("GC/AT-rich region detection failed", exc_info=True)

    # 7. RBS (prokaryotic feature)
    try:
        annotations.extend(_find_rbs(dna))
    except Exception:
        logger.warning("RBS detection failed", exc_info=True)

    # Sort by position then feature type
    annotations.sort(key=lambda a: (a.start, a.feature_type))

    logger.debug("Annotated %d features in %d bp sequence", len(annotations), len(dna))
    return annotations


def annotate_to_genbank(
    dna: str,
    name: str = "BioCompiler_design",
    organism: str = "Homo_sapiens",
    annotations: Optional[list[SequenceAnnotation]] = None,
) -> str:
    """Create a fully-annotated GenBank record from a DNA sequence.

    Automatically annotates the sequence using :func:`annotate_sequence`
    (if ``annotations`` is not provided), then creates a GenBank-format
    record with all features properly positioned and qualified.

    Args:
        dna: DNA sequence.
        name: Record name / locus name (max 16 chars).
        organism: Organism name for the SOURCE/ORGANISM fields.
        annotations: Optional pre-computed annotations. If ``None``,
            :func:`annotate_sequence` is called automatically.

    Returns:
        GenBank-formatted string with full feature annotations.
    """
    dna = validate_dna_sequence(dna)
    if not dna:
        return ""

    if annotations is None:
        annotations = annotate_sequence(dna, organism=organism)

    from .core import export_genbank
    from ..scanner import gc_content as _gc_content
    from ..translation import translate as _translate
    from ..constants import CODON_TABLE

    # Extract CDS annotations for the main feature
    cds_annots = [a for a in annotations if a.feature_type == "CDS"]
    exon_boundaries = None
    if cds_annots:
        # Sort by start position
        cds_annots.sort(key=lambda a: a.start)
        exon_boundaries = [(a.start, a.end) for a in cds_annots]

    # Build restriction site info for export
    restriction_sites = []
    for a in annotations:
        if a.feature_type == "restriction_site":
            restriction_sites.append({
                "enzyme": a.qualifiers.get("enzyme", "unknown"),
                "site": a.qualifiers.get("site", ""),
                "position": a.start,
                "strand": "+" if a.strand == 1 else "-",
            })

    # Translate protein
    protein = _translate(dna)

    # Build the base GenBank record
    locus = name[:16].upper()
    gb = export_genbank(
        sequence=dna,
        locus_name=locus,
        definition=f"BioCompiler annotated sequence for {organism}",
        organism=organism,
        exon_boundaries=exon_boundaries,
        restriction_sites=restriction_sites if restriction_sites else None,
        gene_name=name if name != "BioCompiler_design" else None,
        protein=protein,
    )

    # Add extra annotations not covered by export_genbank
    extra_features: list[str] = []
    for a in annotations:
        if a.feature_type in ("CDS", "exon", "gene", "restriction_site"):
            continue  # already handled by export_genbank

        # Convert 0-based [start, end) to 1-based [start+1, end]
        start_1 = a.start + 1
        end_1 = a.end
        location = f"{start_1}..{end_1}"

        if a.feature_type == "CpG_island":
            extra_features.append(f"     misc_feature    {location}")
            extra_features.append(f'                     /note="CpG island: GC={a.qualifiers.get("gc_content", "N/A")}"')
            extra_features.append(f'                     /label="CpG_island"')

        elif a.feature_type == "splice_donor":
            extra_features.append(f"     misc_feature    {location}")
            extra_features.append(f'                     /note="Splice donor GT: score={a.qualifiers.get("score", "N/A")}"')
            extra_features.append(f'                     /label="splice_donor"')

        elif a.feature_type == "splice_acceptor":
            extra_features.append(f"     misc_feature    {location}")
            extra_features.append(f'                     /note="Splice acceptor AG: score={a.qualifiers.get("score", "N/A")}"')
            extra_features.append(f'                     /label="splice_acceptor"')

        elif a.feature_type == "simple_repeat":
            extra_features.append(f"     misc_feature    {location}")
            extra_features.append(f'                     /note="Simple repeat: {a.qualifiers.get("repeat_unit", "?")}×{a.qualifiers.get("copies", "?")}"')
            extra_features.append(f'                     /label="simple_repeat"')

        elif a.feature_type == "GC_rich":
            extra_features.append(f"     misc_feature    {location}")
            extra_features.append(f'                     /note="GC-rich region: GC={a.qualifiers.get("gc_content", "N/A")}"')
            extra_features.append(f'                     /label="GC_rich"')

        elif a.feature_type == "AT_rich":
            extra_features.append(f"     misc_feature    {location}")
            extra_features.append(f'                     /note="AT-rich region: AT={a.qualifiers.get("at_content", "N/A")}"')
            extra_features.append(f'                     /label="AT_rich"')

        elif a.feature_type == "RBS":
            extra_features.append(f"     misc_feature    {location}")
            extra_features.append(f'                     /note="RBS: {a.qualifiers.get("pattern", "?")} (Shine-Dalgarno)"')
            extra_features.append(f'                     /label="RBS"')

    # Insert extra features before ORIGIN section
    if extra_features:
        origin_idx = gb.find("ORIGIN")
        if origin_idx > 0:
            insert_pos = gb.rfind("\n", 0, origin_idx)
            feature_block = "\n".join(extra_features) + "\n"
            gb = gb[:insert_pos] + "\n" + feature_block + gb[insert_pos:]

    return gb
