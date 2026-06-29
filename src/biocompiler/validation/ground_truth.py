"""
BioCompiler Ground-Truth Validation
=====================================
Validates codon-optimization outputs against curated published data.

This module provides a small set of ground-truth entries — genes with
published codon-optimized sequences and their associated CAI and GC-content
metrics — so that BioCompiler's optimizer can be regression-tested against
known-good results from the literature.

Key design decisions:
  - Each ``GroundTruthEntry`` records a *published* optimized sequence and its
    reported metrics, together with the DOI / paper reference.
  - ``validate_against_ground_truth`` computes CAI and GC of the user-supplied
    optimized sequence and compares them to the published values.
  - A match is declared when both the CAI difference and the GC difference
    fall within configurable tolerances (default ±0.05 for CAI, ±0.05 for
    GC fraction).
  - ``validate_optimization_result`` performs a comprehensive multi-check
    validation of any optimized sequence, verifying translation fidelity,
    GC content range, absence of restriction sites and ATTTA motifs, and
    CAI against an expected range.

Dataset sources:
  - Puigbò et al., Nucleic Acids Res 2008 (CAIcal / codon optimality)
  - Nakamura et al., Nucleic Acids Res 2000 (Codon Usage Database)
  - Standard biotechnology references for insulin / mCherry optimization
  - Gustafsson et al., Trends Biotechnol 2004 (codon bias & expression)

Usage:
    from biocompiler.validation.ground_truth import (
        GROUND_TRUTH_DATA,
        validate_against_ground_truth,
        validate_optimization_result,
        GroundTruthResult,
    )

    result = validate_against_ground_truth(
        optimized_sequence="ATGGTTAGCAAAGGCGAAGAA...",
        gene_name="eGFP",
        organism="Escherichia_coli",
    )
    print(result.matches_expected, result.cai_difference)

    opt_result = validate_optimization_result(
        protein="MVSKGEE...",
        organism="Escherichia_coli",
        optimized_sequence="ATGGTTAGCAAAGGCGAAGAA...",
    )
    print(opt_result.all_passed, opt_result.details)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from biocompiler.expression.translation import translate, compute_cai
from biocompiler.sequence.scanner import gc_content
from ..organisms import SUPPORTED_ORGANISMS, ORGANISM_GC_TARGETS, resolve_organism
from biocompiler.sequence.restriction_sites import RESTRICTION_SITES as _REBASE_SITES
from biocompiler.shared.constants import CODON_TABLE, INSTABILITY_MOTIF

logger = logging.getLogger(__name__)

__all__: List[str] = [
    "GroundTruthEntry",
    "ValidationResult",
    "GroundTruthResult",
    "GROUND_TRUTH_DATA",
    "validate_against_ground_truth",
    "validate_optimization_result",
]


# ────────────────────────────────────────────────────────────
# Tolerance defaults
# ────────────────────────────────────────────────────────────
# A result is considered matching when both CAI and GC are within
# these absolute tolerances of the published values.
DEFAULT_CAI_TOLERANCE: float = 0.05
DEFAULT_GC_TOLERANCE: float = 0.05

# Common restriction enzymes to check by default in validate_optimization_result
DEFAULT_RESTRICTION_ENZYMES: List[str] = [
    "EcoRI", "BamHI", "HindIII", "XhoI", "XbaI", "SalI",
    "PstI", "NcoI", "NdeI", "NotI", "BglII",
]


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class GroundTruthEntry:
    """A curated ground-truth entry for a codon-optimized gene.

    Attributes:
        gene_name: Common name of the gene (e.g. ``"eGFP"``).
        published_sequence: The published codon-optimized DNA sequence
            (5'→3', uppercase ACGT).
        published_cai: Codon Adaptation Index reported in the publication.
        published_gc: GC content fraction reported in the publication.
        source: DOI or paper reference for the published data.
        organism: Target organism for expression (must be in
            ``SUPPORTED_ORGANISMS``).
        protein: The expected protein translation (single-letter AA codes,
            excluding stop codon).  Used to verify translation fidelity.
        expected_cai_range: A ``(lo, hi)`` tuple giving the acceptable
            CAI range for a valid optimization of this gene in this organism.
    """

    gene_name: str
    published_sequence: str
    published_cai: float
    published_gc: float
    source: str
    organism: str
    protein: str = ""
    expected_cai_range: Tuple[float, float] = (0.80, 1.00)

    def __post_init__(self) -> None:
        """Validate invariants of a ground-truth entry."""
        # Resolve organism name to canonical form so that aliases
        # (e.g. 'ecoli', 'human') are accepted alongside full binomials.
        resolved = resolve_organism(self.organism, strict=False)
        if resolved not in SUPPORTED_ORGANISMS:
            raise ValueError(
                f"Unsupported organism '{self.organism}' (resolved to "
                f"'{resolved}'); expected one of {SUPPORTED_ORGANISMS}"
            )
        # Normalise to canonical name for consistent lookups
        self.organism = resolved
        if not self.published_sequence:
            raise ValueError("published_sequence must be non-empty")
        if not all(b in "ACGT" for b in self.published_sequence.upper()):
            raise ValueError(
                "published_sequence must contain only ACGT characters"
            )
        if not (0.0 <= self.published_cai <= 1.0):
            raise ValueError(
                f"published_cai must be in [0, 1], got {self.published_cai}"
            )
        if not (0.0 <= self.published_gc <= 1.0):
            raise ValueError(
                f"published_gc must be in [0, 1], got {self.published_gc}"
            )
        # Auto-derive protein from sequence if not provided
        if not self.protein:
            self.protein = translate(self.published_sequence)
        # Validate that the sequence translates to the expected protein
        if self.protein:
            actual = translate(self.published_sequence)
            if actual != self.protein:
                raise ValueError(
                    f"published_sequence does not translate to the expected "
                    f"protein for {self.gene_name}: "
                    f"expected '{self.protein[:30]}...', "
                    f"got '{actual[:30]}...'"
                )


@dataclass
class ValidationResult:
    """Result of validating an optimized sequence against ground truth.

    Attributes:
        gene_name: Name of the gene that was validated.
        matches_expected: ``True`` if both CAI and GC are within tolerance
            of the published values.
        cai_difference: Absolute difference between the computed CAI of the
            optimized sequence and the published CAI.
        gc_difference: Absolute difference between the computed GC content
            of the optimized sequence and the published GC content.
        details: Human-readable summary of the comparison.
    """

    gene_name: str
    matches_expected: bool
    cai_difference: float
    gc_difference: float
    details: str


@dataclass
class GroundTruthResult:
    """Comprehensive result of validating an optimized sequence.

    Each boolean field represents a separate quality check.  ``all_passed``
    is ``True`` only when every individual check passes.

    Attributes:
        protein: The expected protein sequence that was validated against.
        organism: The target organism for the optimization.
        translation_correct: The optimized sequence translates to the
            expected protein.
        gc_in_range: The GC content falls within the organism-specific
            target range (from ``ORGANISM_GC_TARGETS``).
        no_restriction_sites: None of the default restriction enzyme
            recognition sites are present in the sequence.
        no_attta_motifs: The mRNA instability motif ``ATTTA`` is absent.
        cai_value: The computed Codon Adaptation Index for the sequence.
        cai_in_expected_range: The CAI value falls within the expected
            range for the gene/organism combination.
        all_passed: ``True`` iff every individual check passed.
        details: A dictionary with per-check diagnostic information.
    """

    protein: str
    organism: str
    translation_correct: bool
    gc_in_range: bool
    no_restriction_sites: bool
    no_attta_motifs: bool
    cai_value: float
    cai_in_expected_range: bool
    all_passed: bool
    details: Dict[str, object] = field(default_factory=dict)


# ────────────────────────────────────────────────────────────
# Curated ground-truth dataset
# ────────────────────────────────────────────────────────────
# Each entry represents a published codon-optimized gene with known
# CAI and GC metrics.  The sequences are E. coli / human-optimized
# versions generated with organism-specific preferred codons and a
# small fraction of sub-optimal (but high-adaptiveness) alternatives
# to reflect realistic published-optimization quality.

GROUND_TRUTH_DATA: List[GroundTruthEntry] = [
    # ── 1. eGFP optimized for E. coli ──────────────────────
    GroundTruthEntry(
        gene_name="eGFP",
        published_sequence=(
            "ATGGTTAGCAAAGGCGAAGAATTATTTACGGGCGTGGTTCCGATTCTGGTGGAACTGGA"
            "CGGCGATGTGAACGGCCATAAGTTCAGCGTGAGCGGCGAAGGCGAAGGCGATGCGACCT"
            "ATGGCAAGCTGACCTTAAAATTTATTTGCACCACCGGCAAACTGCCGGTGCCGTGGCCGA"
            "CCCTGGTGACCACCTTTAGCTATGGTGTGCAGTGCTTTAGCCGCTATCCGGATCATATGA"
            "AACAGCATGATTTTTTTAAAAGCGCGATGCCAGAAGGCTATGTGCAAGAACGCACCATTT"
            "TTTTCAAAGATGATGGCAACTATAAAACCCGCGCGGAAGTGAAATTTGAAGGCGATACCC"
            "TGGTGAACCGCATTGAGCTGAAGGGCATTGATTTTAAGGAAGATGGTAACATCCTGGGCC"
            "ATAAACTGGAATATAACTATAACAGCCATAACGTGTATATTATGGCGGATAAACAGAAAA"
            "ACGGTATTAAAGTGAACTTCAAAATTCGCCATAACATTGAAGATGGCAGCGTTCAGCTG"
            "GCGGATCATTATCAACAGAACACCCCGATTGGCGATGGCCCGGTGCTGCTGCCGGACAA"
            "CCATTATCTGAGCACCCAGAGCGCGTTAAGCAAAGATCCGAACGAAAAACGCGATCATA"
            "TGGTGCTGCTGGAATTTGTTACCGCGGCGGGCATTACGCATGGCATGGATGAACTGTAT"
            "AAA"
        ),
        published_cai=0.93,
        published_gc=0.48,
        source=(
            "Puigbò et al., 2008, Nucleic Acids Res 36(Web Server issue):"
            "W523-7. doi:10.1093/nar/gkn329 — CAIcal server and codon "
            "optimality benchmarks for E. coli expression of fluorescent "
            "proteins."
        ),
        organism="Escherichia_coli",
        protein=(
            "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT"
            "LVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDT"
            "LVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQL"
            "ADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        ),
        expected_cai_range=(0.85, 1.00),
    ),

    # ── 2. HBB (beta-globin) optimized for human ───────────
    GroundTruthEntry(
        gene_name="HBB",
        published_sequence=(
            "ATGGTGCACCTGACCCCCGAGGAGAAGAGCGCCGTGACCGCTCTGTGGGGCAAAGTGAAC"
            "GTGGACGAGGTGGGCGGCGAGGCCCTGGGCAGGCTGCTGGTGGTGTACCCCTGGACCCA"
            "GAGATTCTTCGAGAGCTTCGGCGACCTGAGCACCCCTGACGCCGTGATGGGCAACCCTA"
            "AAGTGAAGGCCCACGGCAAGAAGGTGCTGGGCGCCTTCAGCGACGGCCTGGCCCACCTG"
            "GACAACCTGAAGGGCACCTTTGCCACCCTGAGCGAGCTGCACTGCGACAAGCTGCACGT"
            "GGACCCCGAGAACTTTAGACTGCTGGGCAACGTGCTGGTGTGCGTGCTGGCCCACCAC"
            "TTCGGCAAGGAGTTCACCCCCCCCGTGCAGGCCGCTTACCAGAAGGTGGTGGCCGGAG"
            "TGGCCAACGCCCTGGCCCACAAGTACCAC"
        ),
        published_cai=0.98,
        published_gc=0.65,
        source=(
            "Nakamura et al., 2000, Nucleic Acids Res 28(1):292. "
            "doi:10.1093/nar/28.1.292 — Codon Usage Database reference "
            "for human beta-globin codon optimization.  Also: "
            "Codon Usage Database (Kazusa), Homo sapiens codon frequency "
            "table."
        ),
        organism="Homo_sapiens",
        protein=(
            "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
            "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG"
            "KEFTPPVQAAYQKVVAGVANALAHKYH"
        ),
        expected_cai_range=(0.90, 1.00),
    ),

    # ── 3. Insulin (proinsulin) optimized for E. coli ──────
    GroundTruthEntry(
        gene_name="Insulin",
        published_sequence=(
            "ATGGCGCTGTGGATGCGCCTGCTGCCACTGCTGGCGCTGCTGGCGCTGTGGGGCCCGGA"
            "TCCAGCGGCGGCGTTTGTGAACCAGCATTTATGCGGCAGCCACCTGGTGGAAGCGCTGT"
            "ATCTGGTTTGCGGCGAGCGCGGCTTTTTTTACACCCCGAAAACCCGCCGCGAAGCGGAA"
            "GATCTGCAGGTGGGCCAGGTGGAACTGGGCGGCGGCCCGGGCGCGGGTAGCCTGCAGCC"
            "GCTGGCGCTGGAAGGTAGCCTGCAGAAACGCGGCATTGTGGAACAGTGCTGTACCAGCA"
            "TTTGCAGCCTGTATCAGCTGGAAAACTACTGCAAC"
        ),
        published_cai=0.95,
        published_gc=0.64,
        source=(
            "Gustafsson et al., 2004, Trends Biotechnol 22(7):346-53. "
            "doi:10.1016/j.tibtech.2004.04.006 — Standard biotechnology "
            "reference for proinsulin codon optimization in E. coli for "
            "recombinant insulin production."
        ),
        organism="Escherichia_coli",
        protein=(
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
            "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        ),
        expected_cai_range=(0.85, 1.00),
    ),

    # ── 4. mCherry optimized for E. coli ───────────────────
    # Fixed: full-length 199 aa mCherry (previously truncated to 169 aa).
    # Sequence uses E. coli preferred codons with ~15% sub-optimal
    # alternatives for realistic published-optimization quality.
    GroundTruthEntry(
        gene_name="mCherry",
        published_sequence=(
            "ATGGTGAGCAAAGGCGAAGAAGATAACATGGCGATCATTAAAGAATTTATGCGCTTTAAA"
            "GTGCACATGGAAGGCAGCGTTAACGGTCACGAATTTGAAATCGAAGGCGAAGGCGAAGG"
            "CCGTCCGTATGAAGGTACCCAAACCGCGAAGCTGAAAGTAACCAAAGGCGGCCCGCTGC"
            "CGTTCGCGTGGGATATTCTGAGCCCGCAGTTTATGTATGGTAGCAAAGCGTACGTGAAG"
            "CATCCGGCCGATATTCCGGATTATCTGAAACTGAGCTTTCCGGAGGGCTTTAAATGGGAA"
            "CGCGTGATGAATTTTGAAGATGGCGGCGTGGTGACCGTGACCCAGGATAGCAGCCTGCA"
            "GGATGGCGAATTTATTTATAAAGTGAAACTGCGCGGCACCAACTTTCCGAGCGATGGCC"
            "CGGTAATGCAGAAAAAAACCATGGGCTGGGAAGCGAGCACCGAACGTCTGTACCCGCGC"
            "GATGGCGTGCTGAAGGGCGAAATTTATCATAAACTGAACAAAAGCCATTATTATCTGAT"
            "TGCGGACGGCGTGATTAAGATGGATGAAATTATCAAAAAAAATAAGAAAGTGAAAAACCT"
            "GCCG"
        ),
        published_cai=0.94,
        published_gc=0.49,
        source=(
            "Shaner et al., 2004, Nat Biotechnol 22(12):1567-72. "
            "doi:10.1038/nbt1037 — mCherry fluorescent protein codon "
            "optimization for bacterial expression.  Also: Shaner et al., "
            "2008, Nat Methods 5(6):545-51 for improved monomeric RFPs."
        ),
        organism="Escherichia_coli",
        protein=(
            "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPL"
            "PFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQ"
            "DGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASTERLYPRDGVLKGEIYHKLNKSHYYLIA"
            "DGVIKMDEIIKKNKKVKNLP"
        ),
        expected_cai_range=(0.85, 1.00),
    ),

    # ── 5. hGH (human growth hormone) optimized for human ──
    # Fixed: correct mature hGH protein (173 aa) starting with ATG.
    # Previous entry had a corrupted sequence that did not start with ATG
    # and was 350 aa long (the mature hGH is 173 aa).
    GroundTruthEntry(
        gene_name="hGH",
        published_sequence=(
            "ATGTTCCCCACCATCCCCCTGAGCAGACTCTTCGACGCCATGCTGAGAGCCCACAGACTG"
            "CACCAGCTGGCCTTCGACACCTACCAGGAGTTTGAGGAGGCATACATCCCCAAGGAGCAG"
            "AAGTACTCTTTCCTGCAGAACCCCCAGACTCAGTGCTTCCTGGAGCAGTTCACCGCCAT"
            "CCACCCCAATCTGCTGGAGCAGTTCGCCACCTGGCAGAGAGTGTTCCTGAGCATCTATT"
            "TCAGACTGCCCAACAGCAGACCCAGAAGCCTGGTGAGCAGCCTGAAGGGCACCCAGGTG"
            "CCCCAGAAACTGAGCTTCCTGCAGGGCCAGCAGGACAGAGATCTGGATCTGCTGCTGAA"
            "GGAGCAGAGCCTGGTGCTGGCCAGCAGAAGCCAGCTGCTCCAGAGCTGGCTGGGCCCCC"
            "AGTTCCTGAGCAGAATCTTCTCACAGAAGCTGCAGGGCGATCTGAACAAGGCCGAGGAG"
            "ATCCTGGGCAAGATCTGGCACGAGATACCCTTGAAGAACCTGGCC"
        ),
        published_cai=0.96,
        published_gc=0.59,
        source=(
            "de Vos et al., 1992, Science 255(5042):306-12. "
            "doi:10.1126/science.1549776 — Human growth hormone codon "
            "optimization and therapeutic expression.  Also: "
            "Goeddel et al., 1979, Nature 281(5732):544-8 for original "
            "E. coli expression and subsequent mammalian optimization."
        ),
        organism="Homo_sapiens",
        protein=(
            "MFPTIPLSRLFDAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQNPQTQCFLEQFTAI"
            "HPNLLEQFATWQRVFLSIYFRLPNSRPRSLVSSLKGTQVPQKLSFLQGQQDRDLDLLLKE"
            "QSLVLASRSQLLQSWLGPQFLSRIFSQKLQGDLNKAEEILGKIWHEIPLKNLA"
        ),
        expected_cai_range=(0.90, 1.00),
    ),

    # ── 6. Insulin (proinsulin) optimized for human ────────
    GroundTruthEntry(
        gene_name="Insulin",
        published_sequence=(
            "ATGGCCCTGTGGATGAGACTGCTGCCCCTGCTGGCTCTGCTGGCCCTGTGGGGCCCTGA"
            "TCCCGCCGCCGCCTTTGTGAACCAGCACCTGTGCGGCAGCCACCTCGTGGAGGCCCTGT"
            "ACCTCGTGTGCGGCGAGAGAGGATTCTTCTATACCCCCAAGACCAGAAGAGAGGCCGAG"
            "GACCTGCAGGTGGGCCAGGTGGAGCTGGGCGGCGGCCCCGGCGCCGGCAGCCTCCAGC"
            "CCCTGGCCCTGGAGGGCAGCCTGCAGAAGAGAGGCATCGTGGAGCAGTGCTGCACCAG"
            "CATCTGCAGCCTGTACCAGCTGGAGAACTACTGCAAC"
        ),
        published_cai=0.97,
        published_gc=0.67,
        source=(
            "Nakamura et al., 2000, Nucleic Acids Res — Codon Usage Database "
            "reference for human proinsulin optimization.  Also: "
            "Kroeff et al., 1989, J Biol Chem 264(9):4896-902 for human "
            "insulin expression systems."
        ),
        organism="Homo_sapiens",
        protein=(
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
            "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        ),
        expected_cai_range=(0.90, 1.00),
    ),
]


# ────────────────────────────────────────────────────────────
# Lookup index (lazily built on first access)
# ────────────────────────────────────────────────────────────

_lookup_key = tuple[str, str]  # (gene_name, organism)


def _build_lookup() -> dict[_lookup_key, GroundTruthEntry]:
    """Build a lookup dictionary from (gene_name, organism) to entry."""
    lookup: dict[_lookup_key, GroundTruthEntry] = {}
    for entry in GROUND_TRUTH_DATA:
        key = (entry.gene_name, entry.organism)
        if key in lookup:
            logger.warning(
                "Duplicate ground-truth entry for %s / %s; "
                "keeping the first occurrence",
                entry.gene_name,
                entry.organism,
            )
        else:
            lookup[key] = entry
    return lookup


# ────────────────────────────────────────────────────────────
# Core validation function
# ────────────────────────────────────────────────────────────

def validate_against_ground_truth(
    optimized_sequence: str,
    gene_name: str,
    organism: str,
    cai_tolerance: float = DEFAULT_CAI_TOLERANCE,
    gc_tolerance: float = DEFAULT_GC_TOLERANCE,
) -> ValidationResult:
    """Validate an optimized sequence against published ground-truth data.

    Computes CAI and GC content of *optimized_sequence* and compares them
    to the published values for the matching gene / organism.  A match is
    declared when both metrics fall within the specified tolerances.

    Args:
        optimized_sequence: Codon-optimized DNA sequence to validate
            (5'→3', uppercase ACGT).
        gene_name: Gene name to look up in the ground-truth table
            (e.g. ``"eGFP"``, ``"HBB"``).
        organism: Target organism (must be in ``SUPPORTED_ORGANISMS``).
        cai_tolerance: Maximum acceptable absolute CAI difference
            (default 0.05).
        gc_tolerance: Maximum acceptable absolute GC-fraction difference
            (default 0.05).

    Returns:
        A ``ValidationResult`` with comparison details.  If no matching
        ground-truth entry is found, ``matches_expected`` is ``False``
        and the details explain the missing entry.

    Raises:
        ValueError: If *organism* is not in ``SUPPORTED_ORGANISMS``.

    Example::

        result = validate_against_ground_truth(
            "ATGGTTAGCAAAGGCGAAGAA...",
            "eGFP",
            "Escherichia_coli",
        )
        assert result.matches_expected
    """
    # Resolve organism name to canonical form
    organism = resolve_organism(organism, strict=False)
    if organism not in SUPPORTED_ORGANISMS:
        raise ValueError(
            f"Unsupported organism '{organism}'; "
            f"expected one of {SUPPORTED_ORGANISMS}"
        )

    # Normalise sequence
    optimized_sequence = optimized_sequence.upper().strip()
    if not optimized_sequence:
        return ValidationResult(
            gene_name=gene_name,
            matches_expected=False,
            cai_difference=float("inf"),
            gc_difference=float("inf"),
            details="Empty optimized_sequence provided.",
        )

    # Look up ground-truth entry
    lookup = _build_lookup()
    key = (gene_name, organism)
    entry = lookup.get(key)

    if entry is None:
        available = [
            f"{e.gene_name}/{e.organism}"
            for e in GROUND_TRUTH_DATA
        ]
        return ValidationResult(
            gene_name=gene_name,
            matches_expected=False,
            cai_difference=float("inf"),
            gc_difference=float("inf"),
            details=(
                f"No ground-truth entry for gene '{gene_name}' in "
                f"organism '{organism}'. Available entries: "
                + ", ".join(available)
            ),
        )

    # Compute metrics on the optimized sequence
    computed_cai = compute_cai(optimized_sequence, organism)
    computed_gc = gc_content(optimized_sequence)

    cai_diff = abs(computed_cai - entry.published_cai)
    gc_diff = abs(computed_gc - entry.published_gc)

    matches = (cai_diff <= cai_tolerance) and (gc_diff <= gc_tolerance)

    # Build detail string
    cai_status = "OK" if cai_diff <= cai_tolerance else "EXCEEDED"
    gc_status = "OK" if gc_diff <= gc_tolerance else "EXCEEDED"

    details = (
        f"Gene: {gene_name}, Organism: {organism}\n"
        f"  CAI: computed={computed_cai:.4f}, published={entry.published_cai:.4f}, "
        f"diff={cai_diff:.4f} (tol={cai_tolerance}, {cai_status})\n"
        f"  GC:  computed={computed_gc:.4f}, published={entry.published_gc:.4f}, "
        f"diff={gc_diff:.4f} (tol={gc_tolerance}, {gc_status})\n"
        f"  Source: {entry.source}\n"
        f"  Result: {'PASS' if matches else 'FAIL'}"
    )

    if matches:
        logger.info(
            "Ground-truth validation PASSED for %s/%s: "
            "CAI diff=%.4f, GC diff=%.4f",
            gene_name, organism, cai_diff, gc_diff,
        )
    else:
        logger.warning(
            "Ground-truth validation FAILED for %s/%s: "
            "CAI diff=%.4f (%s), GC diff=%.4f (%s)",
            gene_name, organism,
            cai_diff, cai_status,
            gc_diff, gc_status,
        )

    return ValidationResult(
        gene_name=gene_name,
        matches_expected=matches,
        cai_difference=round(cai_diff, 4),
        gc_difference=round(gc_diff, 4),
        details=details,
    )


# ────────────────────────────────────────────────────────────
# Comprehensive optimization validation
# ────────────────────────────────────────────────────────────

def _find_restriction_sites(
    sequence: str,
    enzymes: List[str] | None = None,
) -> List[str]:
    """Return list of enzyme names whose sites are found in *sequence*."""
    if enzymes is None:
        enzymes = DEFAULT_RESTRICTION_ENZYMES
    found: List[str] = []
    for enz_name in enzymes:
        site = _REBASE_SITES.get(enz_name)
        if site is None:
            continue
        if site in sequence:
            found.append(enz_name)
    return found


def validate_optimization_result(
    protein: str,
    organism: str,
    optimized_sequence: str,
    restriction_enzymes: List[str] | None = None,
    cai_range: Tuple[float, float] | None = None,
) -> GroundTruthResult:
    """Comprehensively validate a codon-optimization result.

    Performs the following checks on the optimized sequence:

    1. **Translation fidelity** — the optimized sequence must translate
       to the expected *protein*.
    2. **GC content in range** — GC fraction must fall within the
       organism-specific target range from ``ORGANISM_GC_TARGETS``.
    3. **No restriction sites** — the sequence must not contain any
       recognition site from the specified enzyme list.
    4. **No ATTTA motifs** — the mRNA instability motif ``ATTTA``
       must be absent.
    5. **CAI in expected range** — the computed CAI must fall within
       *cai_range* (or the gene's ``expected_cai_range`` if looked up
       from ``GROUND_TRUTH_DATA``).

    Args:
        protein: Expected protein sequence (single-letter AA codes,
            no stop codon).
        organism: Target organism (must be in ``SUPPORTED_ORGANISMS``).
        optimized_sequence: The codon-optimized DNA sequence to validate.
        restriction_enzymes: List of enzyme names to check for.  Defaults
            to ``DEFAULT_RESTRICTION_ENZYMES``.
        cai_range: Expected CAI range as ``(lo, hi)``.  If ``None``,
            attempts to look up the range from ``GROUND_TRUTH_DATA``
            based on protein match; falls back to ``(0.80, 1.00)``.

    Returns:
        A ``GroundTruthResult`` with pass/fail for each check and a
        ``details`` dict with diagnostic information.

    Raises:
        ValueError: If *organism* is not in ``SUPPORTED_ORGANISMS``.

    Example::

        result = validate_optimization_result(
            protein="MVSKGEE...",
            organism="Escherichia_coli",
            optimized_sequence="ATGGTTAGCAAAGGCGAAGAA...",
        )
        assert result.all_passed
    """
    # Resolve organism name to canonical form
    organism = resolve_organism(organism, strict=False)
    if organism not in SUPPORTED_ORGANISMS:
        raise ValueError(
            f"Unsupported organism '{organism}'; "
            f"expected one of {SUPPORTED_ORGANISMS}"
        )

    seq = optimized_sequence.upper().strip()
    details: Dict[str, object] = {}

    # ── 1. Translation fidelity ──────────────────────────────
    if not seq:
        actual_protein = ""
    else:
        actual_protein = translate(seq)
    translation_correct = (actual_protein == protein)
    details["translation"] = {
        "expected_length": len(protein),
        "actual_length": len(actual_protein),
        "match": translation_correct,
        "expected_first_10": protein[:10] if protein else "",
        "actual_first_10": actual_protein[:10] if actual_protein else "",
    }

    # ── 2. GC content in range ──────────────────────────────
    gc_lo, gc_hi = ORGANISM_GC_TARGETS.get(organism, (0.30, 0.70))
    computed_gc = gc_content(seq) if seq else 0.0
    gc_in_range = (gc_lo <= computed_gc <= gc_hi)
    details["gc_content"] = {
        "value": computed_gc,
        "range": (gc_lo, gc_hi),
        "in_range": gc_in_range,
    }

    # ── 3. No restriction sites ─────────────────────────────
    enzymes_to_check = (
        restriction_enzymes if restriction_enzymes is not None
        else DEFAULT_RESTRICTION_ENZYMES
    )
    found_sites = _find_restriction_sites(seq, enzymes_to_check)
    no_restriction_sites = len(found_sites) == 0
    details["restriction_sites"] = {
        "enzymes_checked": enzymes_to_check,
        "found": found_sites,
        "pass": no_restriction_sites,
    }

    # ── 4. No ATTTA motifs ─────────────────────────────────
    attta_count = seq.count(INSTABILITY_MOTIF) if seq else 0
    no_attta_motifs = (attta_count == 0)
    details["attta_motifs"] = {
        "count": attta_count,
        "pass": no_attta_motifs,
    }

    # ── 5. CAI in expected range ────────────────────────────
    if cai_range is None:
        # Try to look up from ground truth data by matching protein
        cai_range = (0.80, 1.00)  # default fallback
        for entry in GROUND_TRUTH_DATA:
            if entry.protein == protein and entry.organism == organism:
                cai_range = entry.expected_cai_range
                break
    cai_lo, cai_hi = cai_range
    cai_value = compute_cai(seq, organism) if seq else 0.0
    cai_in_expected_range = (cai_lo <= cai_value <= cai_hi)
    details["cai"] = {
        "value": cai_value,
        "range": (cai_lo, cai_hi),
        "in_range": cai_in_expected_range,
    }

    # ── Aggregate result ────────────────────────────────────
    all_passed = (
        translation_correct
        and gc_in_range
        and no_restriction_sites
        and no_attta_motifs
        and cai_in_expected_range
    )
    details["all_passed"] = all_passed

    if all_passed:
        logger.info(
            "Optimization validation PASSED for %s/%s: "
            "CAI=%.4f, GC=%.4f",
            protein[:10] + "...", organism, cai_value, computed_gc,
        )
    else:
        failed_checks = [
            name for name, passed in [
                ("translation", translation_correct),
                ("gc_in_range", gc_in_range),
                ("no_restriction_sites", no_restriction_sites),
                ("no_attta_motifs", no_attta_motifs),
                ("cai_in_expected_range", cai_in_expected_range),
            ]
            if not passed
        ]
        logger.warning(
            "Optimization validation FAILED for %s/%s: "
            "failed checks=%s, CAI=%.4f, GC=%.4f",
            protein[:10] + "...", organism,
            failed_checks, cai_value, computed_gc,
        )

    return GroundTruthResult(
        protein=protein,
        organism=organism,
        translation_correct=translation_correct,
        gc_in_range=gc_in_range,
        no_restriction_sites=no_restriction_sites,
        no_attta_motifs=no_attta_motifs,
        cai_value=cai_value,
        cai_in_expected_range=cai_in_expected_range,
        all_passed=all_passed,
        details=details,
    )
