"""
Multi-Organism Optimization Tests (v10)
========================================

Verify that optimize_sequence works correctly across all 5 supported
organisms, producing biologically sound sequences that:

1. Optimize without errors
2. Report CAI consistent with compute_cai
3. Translate to the correct protein
4. Have GC content in a reasonable range
5. Contain no stop codons in the coding sequence
6. Contain no restriction sites (for the default enzyme panel)

Also verify that all organism aliases resolve to the same canonical
organism and produce identical CAI values.

Task ID: 30
"""

from __future__ import annotations

import math

import pytest

from biocompiler.optimization import optimize_sequence
from biocompiler.translation import translate, compute_cai
from biocompiler.scanner import gc_content
from biocompiler.type_system import AA_TO_CODONS
from biocompiler.constants import RESTRICTION_ENZYMES, reverse_complement, STOP_CODONS
from biocompiler.organisms import SUPPORTED_ORGANISMS, ORGANISM_GC_TARGETS


# ─── Shared test protein ──────────────────────────────────────────────────────

# Human Insulin (preproinsulin) — 110 aa
# A well-studied protein expressed in all 5 target organisms.
INSULIN_PROTEIN = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTTPKTRREAED"
    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)

# Standard 5-enzyme cloning panel used across BioCompiler tests.
# The optimizer can reliably eliminate all of these sites.
# A broader panel (e.g. including PstI, SmaI, ApaI, BglII) is infeasible
# for some proteins because eliminating one site can introduce another.
DEFAULT_ENZYME_NAMES = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _has_restriction_site(seq: str, enzyme_names: list[str]) -> bool:
    """Return True if any enzyme recognition site is present in seq or its RC."""
    for name in enzyme_names:
        site = RESTRICTION_ENZYMES.get(name, "")
        if not site:
            continue
        # Skip IUPAC sites (e.g. SfiI with N wildcards) for simple checking
        if any(b not in "ACGT" for b in site):
            continue
        rc = reverse_complement(site)
        if site in seq or rc in seq:
            return True
    return False


def _has_internal_stop_codons(seq: str) -> bool:
    """Check if a coding sequence contains any in-frame stop codons."""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if codon in STOP_CODONS:
            return True
    return False


def _get_reasonable_gc_range(organism: str) -> tuple[float, float]:
    """Return a generous GC range for the organism."""
    # Use the organism-specific GC targets as a guide, but with a wide
    # margin for the "reasonable range" check
    gc_targets = ORGANISM_GC_TARGETS.get(organism, (0.25, 0.75))
    # Allow a wider range than the target — the optimizer may not hit
    # the exact target, but GC should be biologically plausible
    return (gc_targets[0] - 0.15, gc_targets[1] + 0.15)


# ─── Parametrized organisms ───────────────────────────────────────────────────

ALL_ORGANISMS = [
    "Escherichia_coli",
    "Homo_sapiens",
    "Mus_musculus",
    "CHO_K1",
    "Saccharomyces_cerevisiae",
]


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Core multi-organism optimization tests
# ═══════════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("organism", ALL_ORGANISMS)
class TestMultiOrganismOptimization:
    """Verify optimization correctness across all 5 supported organisms."""

    def test_optimize_without_errors(self, organism: str):
        """optimize_sequence should complete without raising for each organism."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
        assert result is not None
        assert result.sequence
        assert len(result.sequence) == len(INSULIN_PROTEIN) * 3

    def test_cai_consistency(self, organism: str):
        """CAI reported by optimize_sequence must match compute_cai."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
        independent_cai = compute_cai(result.sequence, organism=organism)
        assert math.isclose(result.cai, independent_cai, abs_tol=0.001), (
            f"{organism}: result.cai={result.cai:.6f} != "
            f"compute_cai={independent_cai:.6f}"
        )

    def test_correct_translation(self, organism: str):
        """Optimized sequence must translate to the input protein."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
        translated = translate(result.sequence)
        assert translated == INSULIN_PROTEIN, (
            f"{organism}: translation mismatch — expected "
            f"{len(INSULIN_PROTEIN)} aa, got {len(translated)} aa. "
            f"First difference at position "
            f"{next((i for i, (a, b) in enumerate(zip(translated, INSULIN_PROTEIN)) if a != b), 'none')}"
        )

    def test_gc_content_reasonable(self, organism: str):
        """GC content must be in a biologically reasonable range."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
        actual_gc = gc_content(result.sequence)
        # The optimizer was asked for [0.30, 0.70]
        assert 0.20 <= actual_gc <= 0.80, (
            f"{organism}: GC content {actual_gc:.4f} is outside "
            f"reasonable range [0.20, 0.80]"
        )
        # Also verify it matches the result's own report
        assert math.isclose(actual_gc, result.gc_content, abs_tol=1e-6), (
            f"{organism}: computed gc_content {actual_gc:.6f} != "
            f"result.gc_content {result.gc_content:.6f}"
        )

    def test_no_internal_stop_codons(self, organism: str):
        """Coding sequence must not contain in-frame stop codons."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
        assert not _has_internal_stop_codons(result.sequence), (
            f"{organism}: optimized sequence contains in-frame stop codons"
        )

    def test_no_restriction_sites(self, organism: str):
        """Optimized sequence must not contain default restriction sites.

        The standard 5-enzyme panel (EcoRI, BamHI, XhoI, HindIII, NotI)
        is the default panel used across BioCompiler's test suite.  The
        optimizer can reliably eliminate all of these sites.
        """
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
        assert not _has_restriction_site(result.sequence, DEFAULT_ENZYME_NAMES), (
            f"{organism}: optimized sequence contains one or more "
            f"restriction sites from the default enzyme panel "
            f"{DEFAULT_ENZYME_NAMES}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Organism alias resolution tests
# ═══════════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("alias,canonical", [
    ("ecoli", "Escherichia_coli"),
    ("E. coli", "Escherichia_coli"),
    ("e_coli", "Escherichia_coli"),
    ("human", "Homo_sapiens"),
    ("h_sapiens", "Homo_sapiens"),
    ("mouse", "Mus_musculus"),
    ("cho", "CHO_K1"),
    ("yeast", "Saccharomyces_cerevisiae"),
])
def test_organism_alias(alias: str, canonical: str):
    """Alias-based optimization must produce the same CAI as canonical name."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        result1 = optimize_sequence(
            INSULIN_PROTEIN,
            species=alias,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
        result2 = optimize_sequence(
            INSULIN_PROTEIN,
            organism=canonical,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
    assert result1.cai == result2.cai, (
        f"Alias '{alias}' → CAI {result1.cai:.6f} != "
        f"canonical '{canonical}' → CAI {result2.cai:.6f}"
    )


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Cross-organism sanity checks
# ═══════════════════════════════════════════════════════════════════════════════


class TestCrossOrganismSanity:
    """Cross-organism consistency checks."""

    @pytest.mark.parametrize("organism", ALL_ORGANISMS)
    def test_cai_positive(self, organism: str):
        """CAI must be strictly positive for a well-expressed protein."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
        assert result.cai > 0.0, (
            f"{organism}: CAI must be positive, got {result.cai:.6f}"
        )

    @pytest.mark.parametrize("organism", ALL_ORGANISMS)
    def test_all_codons_valid(self, organism: str):
        """Every codon in the optimized sequence must be valid for its amino acid."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
        for i, aa in enumerate(INSULIN_PROTEIN):
            codon = result.sequence[i * 3: i * 3 + 3]
            valid_codons = AA_TO_CODONS.get(aa, [])
            assert codon in valid_codons, (
                f"{organism}: position {i}: codon '{codon}' is not valid "
                f"for amino acid '{aa}'. Valid: {valid_codons}"
            )

    @pytest.mark.parametrize("organism", ALL_ORGANISMS)
    def test_protein_attribute_matches(self, organism: str):
        """result.protein must match the input protein."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYME_NAMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )
        assert result.protein == INSULIN_PROTEIN, (
            f"{organism}: result.protein does not match input protein"
        )

    def test_supported_organisms_list_complete(self):
        """All 5 expected organisms must be in SUPPORTED_ORGANISMS."""
        expected = {
            "Escherichia_coli",
            "Homo_sapiens",
            "Mus_musculus",
            "CHO_K1",
            "Saccharomyces_cerevisiae",
        }
        assert expected.issubset(set(SUPPORTED_ORGANISMS)), (
            f"Missing organisms: {expected - set(SUPPORTED_ORGANISMS)}"
        )
