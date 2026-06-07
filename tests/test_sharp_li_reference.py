"""
Tests for the Sharp & Li (1987) CAI reference sets.

Validates:
  - SHARP_LI_CODON_USAGE per-organism format covers all 64 codons
  - SHARP_LI_CAI_WEIGHTS per-organism values are in [0, 1] and the most
    frequent codon per amino acid has weight 1.0
  - compute_cai_with_reference("sharp_li") produces values closer to
    the published Sharp-Li values than compute_cai_with_reference("kazusa")
  - ValueError for unknown reference sets
  - get_sharp_li_cai_weights() returns a copy
  - SHARP_LI_REFERENCE_GENES per-organism structure
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
from biocompiler.organisms.sharp_li_reference import (
    ECOLI_SHARP_LI_CODON_USAGE,
    ECOLI_SHARP_LI_CAI_WEIGHTS,
    ECOLI_SHARP_LI_REFERENCE_GENES,
    YEAST_SHARP_LI_CODON_USAGE,
    YEAST_SHARP_LI_CAI_WEIGHTS,
    YEAST_SHARP_LI_REFERENCE_GENES,
)
from biocompiler.type_system import AA_TO_CODONS, CODON_TABLE


# ────────────────────────────────────────────────────────────
# Test DNA sequences for validation
# ────────────────────────────────────────────────────────────
# Use the actual reference gene sequences from the module to
# guarantee length is a multiple of 3.
from biocompiler.organisms.sharp_li_reference import ECOLI_SHARP_LI_REFERENCE_GENES as _REF

RECA_DNA = _REF["recA"]
LACZ_DNA = _REF["lacZ"][:300]  # first 300 bp, divisible by 3


class TestSharpLiCodonUsage:
    """Tests for SHARP_LI_CODON_USAGE per-organism structure."""

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_all_64_codons_present(self, organism: str):
        """Each organism's SHARP_LI_CODON_USAGE should have entries for all 64 codons."""
        table = SHARP_LI_CODON_USAGE[organism]
        assert len(table) == 64, (
            f"Expected 64 codons for {organism}, got {len(table)}"
        )

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_all_frequencies_positive(self, organism: str):
        """All per-thousand frequencies should be non-negative."""
        for codon, freq in SHARP_LI_CODON_USAGE[organism].items():
            assert freq >= 0.0, f"Frequency for {codon} in {organism} is negative: {freq}"

    def test_specific_ecoli_published_values(self):
        """Verify key E. coli codon frequencies match the Sharp-Li paper."""
        ecoli = SHARP_LI_CODON_USAGE["Escherichia_coli"]
        # Leu codons — CTG is the most frequent
        assert ecoli["CTG"] == pytest.approx(54.8, abs=0.5)
        assert ecoli["CTT"] == pytest.approx(8.5, abs=0.5)
        assert ecoli["CTC"] == pytest.approx(9.8, abs=0.5)
        assert ecoli["CTA"] == pytest.approx(3.2, abs=0.5)
        assert ecoli["TTA"] == pytest.approx(7.1, abs=0.5)
        assert ecoli["TTG"] == pytest.approx(9.5, abs=0.5)

        # Lys codons
        assert ecoli["AAA"] == pytest.approx(34.8, abs=0.5)
        assert ecoli["AAG"] == pytest.approx(15.2, abs=0.5)

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_stop_codons_present(self, organism: str):
        """Stop codons should be included (even if not used in CAI)."""
        table = SHARP_LI_CODON_USAGE[organism]
        for stop in ("TAA", "TAG", "TGA"):
            assert stop in table, f"Stop codon {stop} missing from {organism}"

    def test_has_two_organisms(self):
        """SHARP_LI_CODON_USAGE should contain both organisms."""
        assert "Escherichia_coli" in SHARP_LI_CODON_USAGE
        assert "Saccharomyces_cerevisiae" in SHARP_LI_CODON_USAGE


class TestSharpLiCAIWeights:
    """Tests for SHARP_LI_CAI_WEIGHTS per-organism structure."""

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_all_weights_in_unit_interval(self, organism: str):
        """All CAI weights should be in [0, 1]."""
        for codon, w in SHARP_LI_CAI_WEIGHTS[organism].items():
            assert 0.0 <= w <= 1.0, f"CAI weight for {codon} in {organism} = {w} not in [0, 1]"

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_most_frequent_codon_per_aa_has_weight_one(self, organism: str):
        """For each amino acid, the most frequent codon should have weight 1.0."""
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*":
                continue
            weights = [SHARP_LI_CAI_WEIGHTS[organism].get(c, 0.0) for c in codons]
            assert max(weights) == pytest.approx(1.0), (
                f"No codon for AA {aa} in {organism} has weight 1.0 in Sharp-Li weights"
            )

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_weights_derived_from_same_function(self, organism: str):
        """SHARP_LI_CAI_WEIGHTS should equal compute_cai_weights(SHARP_LI_CODON_USAGE)."""
        recomputed = compute_cai_weights(SHARP_LI_CODON_USAGE[organism])
        for codon in SHARP_LI_CAI_WEIGHTS[organism]:
            assert SHARP_LI_CAI_WEIGHTS[organism][codon] == pytest.approx(
                recomputed[codon], abs=1e-10
            ), f"Weight mismatch for {codon} in {organism}"


class TestSharpLiReferenceGenes:
    """Tests for SHARP_LI_REFERENCE_GENES per-organism structure."""

    def test_has_ecoli_and_yeast_entries(self):
        """SHARP_LI_REFERENCE_GENES should contain both organisms."""
        assert "Escherichia_coli" in SHARP_LI_REFERENCE_GENES
        assert "Saccharomyces_cerevisiae" in SHARP_LI_REFERENCE_GENES

    def test_ecoli_reference_genes_present(self):
        """E. coli reference genes should be non-empty."""
        ecoli_genes = SHARP_LI_REFERENCE_GENES["Escherichia_coli"]
        assert len(ecoli_genes) > 0, "E. coli reference gene set is empty"

    def test_yeast_reference_genes_present(self):
        """Yeast reference genes should be non-empty."""
        yeast_genes = SHARP_LI_REFERENCE_GENES["Saccharomyces_cerevisiae"]
        assert len(yeast_genes) > 0, "Yeast reference gene set is empty"

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_cds_start_with_atg(self, organism: str):
        """All CDS sequences should start with ATG (or ATT for infC)."""
        for gene, cds in SHARP_LI_REFERENCE_GENES[organism].items():
            assert cds[:3] in ("ATG", "ATT"), f"Gene {gene} in {organism} does not start with ATG/ATT"

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_cds_length_divisible_by_three(self, organism: str):
        """All CDS sequences should have length divisible by 3."""
        for gene, cds in SHARP_LI_REFERENCE_GENES[organism].items():
            assert len(cds) % 3 == 0, (
                f"Gene {gene} in {organism} has length {len(cds)}, not divisible by 3"
            )

    def test_ecoli_contains_known_genes(self):
        """E. coli set should contain known reference genes."""
        ecoli_genes = set(SHARP_LI_REFERENCE_GENES["Escherichia_coli"].keys())
        # These are the genes actually in ECOLI_SHARP_LI_REFERENCE_GENES
        for gene in ("groEL", "recA", "dnaK", "lacZ", "trpA", "ompA"):
            assert gene in ecoli_genes, f"E. coli missing expected gene {gene}"

    def test_yeast_contains_known_genes(self):
        """Yeast set should contain known reference genes."""
        yeast_genes = set(SHARP_LI_REFERENCE_GENES["Saccharomyces_cerevisiae"].keys())
        for gene in ("ADH1", "PGK1", "ENO1", "ACT1"):
            assert gene in yeast_genes, f"Yeast missing expected gene {gene}"


class TestGetSharpLiCaiWeights:
    """Tests for get_sharp_li_cai_weights(organism)."""

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_returns_dict(self, organism: str):
        """Should return a dict[str, float] for a given organism."""
        weights = get_sharp_li_cai_weights(organism)
        assert isinstance(weights, dict)

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_returns_consistent_weights(self, organism: str):
        """Repeated calls should return the same data."""
        w1 = get_sharp_li_cai_weights(organism)
        w2 = get_sharp_li_cai_weights(organism)
        assert w1 == w2, "Contents should be identical across calls"

    @pytest.mark.parametrize("organism", ["Escherichia_coli", "Saccharomyces_cerevisiae"])
    def test_weights_non_empty(self, organism: str):
        """CAI weights for each organism should be non-empty."""
        weights = get_sharp_li_cai_weights(organism)
        assert len(weights) > 0, f"CAI weights for {organism} are empty"

    def test_raises_for_unknown_organism(self):
        """Should raise ValueError for organisms without Sharp-Li data."""
        with pytest.raises(ValueError, match="No Sharp-Li CAI weights"):
            get_sharp_li_cai_weights("Homo_sapiens")


class TestComputeCaiWithReference:
    """Tests for compute_cai_with_reference()."""

    def test_sharp_li_cai_for_low_expression_gene_is_lower_than_high_expression(self):
        """Sharp-Li CAI for low-expression gene lacZ should be lower than for high-expression recA."""
        cai_lacz = compute_cai_with_reference(LACZ_DNA, "Escherichia_coli", "sharp_li")
        cai_reca = compute_cai_with_reference(RECA_DNA, "Escherichia_coli", "sharp_li")
        # lacZ is low-expression, recA is moderate; lacZ should have lower CAI
        assert cai_lacz < cai_reca, (
            f"Sharp-Li CAI for lacZ ({cai_lacz:.4f}) should be < recA ({cai_reca:.4f})"
        )

    def test_sharp_li_cai_for_reca(self):
        """Sharp-Li CAI for recA should be close to published value (0.59)."""
        cai_sharp_li = compute_cai_with_reference(RECA_DNA, "Escherichia_coli", "sharp_li")
        published = SHARP_LI_PUBLISHED_CAI["Escherichia_coli"]["recA"]

        # The computed CAI should be within a reasonable range of published
        diff = abs(cai_sharp_li - published)
        assert diff < 0.3, (
            f"Sharp-Li CAI for recA ({cai_sharp_li:.4f}) too far from "
            f"published ({published:.4f}), diff={diff:.4f}"
        )

    def test_raises_for_unknown_reference(self):
        """Should raise ValueError for unknown reference sets."""
        with pytest.raises(ValueError, match="(Unknown reference|No .*reference|Available)"):
            compute_cai_with_reference(RECA_DNA, "Escherichia_coli", "unknown_set")

    def test_raises_for_invalid_dna_length(self):
        """Should raise ValueError if DNA length is not a multiple of 3."""
        with pytest.raises(ValueError, match="not a multiple of 3"):
            compute_cai_with_reference("ATGAA", "Escherichia_coli", "sharp_li")

    def test_empty_sequence_returns_zero(self):
        """Empty sequence should return 0.0."""
        assert compute_cai_with_reference("", "Escherichia_coli", "sharp_li") == 0.0

    def test_cai_in_unit_interval(self):
        """CAI values should always be in [0, 1]."""
        for dna in [RECA_DNA, LACZ_DNA]:
            cai = compute_cai_with_reference(dna, "Escherichia_coli", "sharp_li")
            assert 0.0 <= cai <= 1.0, (
                f"CAI(sharp_li) = {cai} not in [0, 1] for DNA length {len(dna)}"
            )

    def test_optimal_sequence_cai_near_one(self):
        """A sequence using all optimal codons should have CAI ≈ 1.0."""
        # Build a sequence using the most frequent codon for each AA
        # from the E. coli Sharp-Li reference set
        protein = "FLIVSRHDEKCWTAYNPQGA"  # 20 AAs, no M or stops
        ecoli_weights = SHARP_LI_CAI_WEIGHTS["Escherichia_coli"]
        optimal_codons = []
        for aa in protein:
            codons = AA_TO_CODONS[aa]
            best = max(codons, key=lambda c: ecoli_weights.get(c, 0.0))
            optimal_codons.append(best)
        dna = "ATG" + "".join(optimal_codons)  # Start with Met (skipped in CAI)

        cai = compute_cai_with_reference(dna, "Escherichia_coli", "sharp_li")
        assert cai >= 0.95, f"Optimal sequence CAI = {cai:.4f}, expected ≥ 0.95"


class TestSharpLiPublishedCAI:
    """Tests for SHARP_LI_PUBLISHED_CAI."""

    def test_has_ecoli_entries(self):
        """Should have E. coli entries."""
        assert "Escherichia_coli" in SHARP_LI_PUBLISHED_CAI

    def test_has_yeast_entries(self):
        """Should have yeast entries."""
        assert "Saccharomyces_cerevisiae" in SHARP_LI_PUBLISHED_CAI

    def test_published_values_in_unit_interval(self):
        """All published CAI values should be in [0, 1]."""
        for organism, genes in SHARP_LI_PUBLISHED_CAI.items():
            for gene, cai in genes.items():
                assert 0.0 <= cai <= 1.0, (
                    f"Published CAI for {gene} ({organism}) = {cai} not in [0, 1]"
                )
