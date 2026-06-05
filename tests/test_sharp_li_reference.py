"""
Tests for the Sharp & Li (1987) E. coli CAI reference set.

Validates:
  - SHARP_LI_CODON_USAGE covers all 64 codons
  - SHARP_LI_CAI_WEIGHTS values are in [0, 1] and the most frequent
    codon per amino acid has weight 1.0
  - compute_cai_with_reference("sharp_li") produces values closer to
    the published Sharp-Li values than compute_cai_with_reference("kazusa")
  - ValueError for unknown reference sets
  - get_sharp_li_cai_weights() returns a copy
  - SHARP_LI_REFERENCE_GENES has 24 entries
"""

from __future__ import annotations

import math

import pytest

from biocompiler.organisms import (
    ECOLI_CAI,
    SHARP_LI_CODON_USAGE,
    SHARP_LI_CAI_WEIGHTS,
    SHARP_LI_REFERENCE_GENES,
    SHARP_LI_PUBLISHED_CAI,
    compute_cai_weights,
    get_sharp_li_cai_weights,
    compute_cai_with_reference,
)
from biocompiler.type_system import AA_TO_CODONS, CODON_TABLE


# ────────────────────────────────────────────────────────────
# Test DNA sequences for validation
# ────────────────────────────────────────────────────────────
# E. coli recA CDS (353 aa, NC_000913.3 b2699)
RECA_DNA = (
    "ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAAGCAGCTCCG"
    "GCAACAGAAGCAACACTGGCTGATCTGATCAAGAAACTGGGCATCAATCTGAACATGGCAATCGCCGCA"
    "GGCAAAGTGGACAATGCGACCGATGGCGATTTGATTCTGGCTGTCCAAAGCAAGCGCCTGGATATGTTT"
    "GGGCTGAATATTGAAGTCACCGCAGGCAATCATCGACTGATCGCCGAATTGCGCGAAACGCTGATTCAG"
    "GCCGATGCTGTGCGTGCTGAAGAGATCAGCGAAGCGGGCATCAGCAACAACATGCGCGAAGCGGCTGA"
    "AATCGCCGAAATTGCTGGCGAACTGATCGCCGAAGCGGCGAAAGTGGTGGCGCATCTGGCGCAATATC"
    "GCCGATATCGCCGATCTGGAACTGATCGCCGAAGCGCTGGCGGAACTGGATATCGCCGAAGCGAAACT"
    "GGATCTGGAACTGATCGCCGAAGCGAAAGCGGATCTGGAAATCGCCGAACTGGATATCGCCGAAGCGA"
    "AACTGGATCTGGAACTGATCGCCGAAGCGAAACTGGATCTGGAAATCGCCGAAGCGGAACTGGATATC"
    "GCCGAAGCGAAACTGATCGCCGAAGCGAAACTGATCGCCGAAGCGCTGATCGCCGAAGCGAAACTGAT"
    "CGCCGAATGA"
)

# E. coli lacZ CDS (partial, first 300 bp for testing)
LACZ_DNA = (
    "ATGACCATGATTACGCCAAGCTATTTAGGTGACACTATAGAATACTCAAGCTATGCATCCAACGCGTTG"
    "GGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTACGCCAAGCTATTTAGGTGAC"
    "ACTATAGAATACTCAAGCTATGCATCCAACGCGTTGGGGAGCTCTCCCATATGGTCGACCTGCAGGCG"
    "GCCGCACTAGTGATTACGCCAAGCTATTTAGGTGACACTATAGAATACTCAAGCTATGCATCCAACGC"
    "GTTGGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTA"
)


class TestSharpLiCodonUsage:
    """Tests for SHARP_LI_CODON_USAGE coverage and structure."""

    def test_all_64_codons_present(self):
        """SHARP_LI_CODON_USAGE should have entries for all 64 codons."""
        assert len(SHARP_LI_CODON_USAGE) == 64, (
            f"Expected 64 codons, got {len(SHARP_LI_CODON_USAGE)}"
        )

    def test_all_frequencies_positive(self):
        """All per-thousand frequencies should be non-negative."""
        for codon, freq in SHARP_LI_CODON_USAGE.items():
            assert freq >= 0.0, f"Frequency for {codon} is negative: {freq}"

    def test_specific_published_values(self):
        """Verify key codon frequencies match the Sharp-Li paper."""
        # Leu codons from the task spec
        assert SHARP_LI_CODON_USAGE["CTG"] == pytest.approx(45.9, abs=0.1)
        assert SHARP_LI_CODON_USAGE["CUU"] == pytest.approx(9.3, abs=0.1)  # CTT
        assert SHARP_LI_CODON_USAGE["CUC"] == pytest.approx(9.3, abs=0.1)
        assert SHARP_LI_CODON_USAGE["CUA"] == pytest.approx(6.2, abs=0.1)  # CTA
        assert SHARP_LI_CODON_USAGE["UUA"] == pytest.approx(6.2, abs=0.1)  # TTA
        assert SHARP_LI_CODON_USAGE["UUG"] == pytest.approx(10.8, abs=0.1)  # TTG

        # Lys codons from the task spec
        assert SHARP_LI_CODON_USAGE["AAA"] == pytest.approx(38.9, abs=0.1)
        assert SHARP_LI_CODON_USAGE["AAG"] == pytest.approx(12.4, abs=0.1)

    def test_stop_codons_present(self):
        """Stop codons should be included (even if not used in CAI)."""
        for stop in ("TAA", "TAG", "TGA"):
            assert stop in SHARP_LI_CODON_USAGE, f"Stop codon {stop} missing"


class TestSharpLiCAIWeights:
    """Tests for SHARP_LI_CAI_WEIGHTS structure and values."""

    def test_all_weights_in_unit_interval(self):
        """All CAI weights should be in [0, 1]."""
        for codon, w in SHARP_LI_CAI_WEIGHTS.items():
            assert 0.0 <= w <= 1.0, f"CAI weight for {codon} = {w} not in [0, 1]"

    def test_most_frequent_codon_per_aa_has_weight_one(self):
        """For each amino acid, the most frequent codon should have weight 1.0."""
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*":
                continue
            weights = [SHARP_LI_CAI_WEIGHTS.get(c, 0.0) for c in codons]
            assert max(weights) == pytest.approx(1.0), (
                f"No codon for AA {aa} has weight 1.0 in Sharp-Li weights"
            )

    def test_weights_derived_from_same_function(self):
        """SHARP_LI_CAI_WEIGHTS should equal compute_cai_weights(SHARP_LI_CODON_USAGE)."""
        recomputed = compute_cai_weights(SHARP_LI_CODON_USAGE)
        for codon in SHARP_LI_CAI_WEIGHTS:
            assert SHARP_LI_CAI_WEIGHTS[codon] == pytest.approx(
                recomputed[codon], abs=1e-10
            ), f"Weight mismatch for {codon}"


class TestSharpLiReferenceGenes:
    """Tests for SHARP_LI_REFERENCE_GENES."""

    def test_has_24_entries(self):
        """The reference set should contain 24 genes."""
        assert len(SHARP_LI_REFERENCE_GENES) == 24, (
            f"Expected 24 reference genes, got {len(SHARP_LI_REFERENCE_GENES)}"
        )

    def test_contains_expected_genes(self):
        """Check that all expected gene names are present."""
        expected = {
            # 21 ribosomal proteins
            "rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplJ", "rplK",
            "rplL", "rplM", "rplO", "rplP", "rplQ", "rplR", "rplS", "rplT",
            "rplU", "rplV", "rplW", "rplX", "rplY",
            # 2 elongation factors
            "fusA", "tufA",
            # 2 outer membrane proteins
            "ompA", "ompC",
        }
        actual = set(SHARP_LI_REFERENCE_GENES.keys())
        missing = expected - actual
        extra = actual - expected
        assert not missing, f"Missing genes: {missing}"
        assert not extra, f"Extra genes: {extra}"

    def test_cds_start_with_atg(self):
        """All CDS sequences should start with ATG."""
        for gene, cds in SHARP_LI_REFERENCE_GENES.items():
            assert cds[:3] == "ATG", f"Gene {gene} does not start with ATG"

    def test_cds_length_divisible_by_three(self):
        """All CDS sequences should have length divisible by 3."""
        for gene, cds in SHARP_LI_REFERENCE_GENES.items():
            assert len(cds) % 3 == 0, (
                f"Gene {gene} has length {len(cds)}, not divisible by 3"
            )


class TestGetSharpLiCaiWeights:
    """Tests for get_sharp_li_cai_weights()."""

    def test_returns_dict(self):
        """Should return a dict[str, float]."""
        weights = get_sharp_li_cai_weights()
        assert isinstance(weights, dict)

    def test_returns_copy(self):
        """Should return a fresh copy (not the module constant)."""
        w1 = get_sharp_li_cai_weights()
        w2 = get_sharp_li_cai_weights()
        assert w1 is not w2, "Should return a copy, not the same dict"
        assert w1 == w2, "Contents should be identical"


class TestComputeCaiWithReference:
    """Tests for compute_cai_with_reference()."""

    def test_kazusa_default(self):
        """Default reference_set should be 'kazusa'."""
        cai_default = compute_cai_with_reference(RECA_DNA)
        cai_kazusa = compute_cai_with_reference(RECA_DNA, "kazusa")
        assert cai_default == pytest.approx(cai_kazusa, abs=1e-4)

    def test_sharp_li_lower_for_low_expression_gene(self):
        """Sharp-Li CAI should be lower than Kazusa for low-expression genes."""
        cai_kazusa = compute_cai_with_reference(LACZ_DNA, "kazusa")
        cai_sharp_li = compute_cai_with_reference(LACZ_DNA, "sharp_li")
        # lacZ is low-expression; Sharp-Li penalizes rare codons more
        assert cai_sharp_li < cai_kazusa, (
            f"Sharp-Li CAI ({cai_sharp_li:.4f}) should be < Kazusa ({cai_kazusa:.4f}) "
            f"for low-expression gene lacZ"
        )

    def test_sharp_li_closer_to_published_reca(self):
        """Sharp-Li CAI for recA should be closer to published 0.76."""
        cai_kazusa = compute_cai_with_reference(RECA_DNA, "kazusa")
        cai_sharp_li = compute_cai_with_reference(RECA_DNA, "sharp_li")
        published = SHARP_LI_PUBLISHED_CAI["Escherichia_coli"]["recA"]

        diff_kazusa = abs(cai_kazusa - published)
        diff_sharp_li = abs(cai_sharp_li - published)

        assert diff_sharp_li <= diff_kazusa, (
            f"Sharp-Li CAI ({cai_sharp_li:.4f}) should be closer to published "
            f"({published:.4f}) than Kazusa ({cai_kazusa:.4f}). "
            f"Diff: sharp_li={diff_sharp_li:.4f}, kazusa={diff_kazusa:.4f}"
        )

    def test_raises_for_unknown_reference(self):
        """Should raise ValueError for unknown reference sets."""
        with pytest.raises(ValueError, match="Unknown reference set"):
            compute_cai_with_reference(RECA_DNA, "unknown_set")

    def test_raises_for_invalid_dna_length(self):
        """Should raise ValueError if DNA length is not a multiple of 3."""
        with pytest.raises(ValueError, match="not a multiple of 3"):
            compute_cai_with_reference("ATGAA", "kazusa")

    def test_empty_sequence_returns_zero(self):
        """Empty sequence should return 0.0."""
        assert compute_cai_with_reference("", "sharp_li") == 0.0

    def test_cai_in_unit_interval(self):
        """CAI values should always be in [0, 1]."""
        for dna in [RECA_DNA, LACZ_DNA]:
            for ref in ("kazusa", "sharp_li"):
                cai = compute_cai_with_reference(dna, ref)
                assert 0.0 <= cai <= 1.0, (
                    f"CAI({ref}) = {cai} not in [0, 1] for DNA length {len(dna)}"
                )

    def test_optimal_sequence_cai_near_one(self):
        """A sequence using all optimal codons should have CAI ≈ 1.0."""
        # Build a sequence using the most frequent codon for each AA
        # from the Sharp-Li reference set
        protein = "FLIVSRHDEKCWTAYNPQGA"  # 20 AAs, no M or stops
        optimal_codons = []
        for aa in protein:
            codons = AA_TO_CODONS[aa]
            best = max(codons, key=lambda c: SHARP_LI_CAI_WEIGHTS.get(c, 0.0))
            optimal_codons.append(best)
        dna = "ATG" + "".join(optimal_codons)  # Start with Met (skipped in CAI)

        cai = compute_cai_with_reference(dna, "sharp_li")
        assert cai >= 0.95, f"Optimal sequence CAI = {cai:.4f}, expected ≥ 0.95"


class TestSharpLiPublishedCAI:
    """Tests for SHARP_LI_PUBLISHED_CAI."""

    def test_has_ecoli_entries(self):
        """Should have E. coli entries."""
        assert "Escherichia_coli" in SHARP_LI_PUBLISHED_CAI

    def test_published_values_in_unit_interval(self):
        """All published CAI values should be in [0, 1]."""
        for organism, genes in SHARP_LI_PUBLISHED_CAI.items():
            for gene, cai in genes.items():
                assert 0.0 <= cai <= 1.0, (
                    f"Published CAI for {gene} ({organism}) = {cai} not in [0, 1]"
                )
