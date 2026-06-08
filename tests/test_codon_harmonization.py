"""
Tests for codon harmonization (HarmonizeRCA / Claassens method).

Tests cover:
- RCA computation with known values for E. coli
- Harmonization produces valid DNA (correct length, translates back)
- Harmonization score is higher for harmonized sequences than CAI-maximized
- harmonize_with_cai_fallback respects cai_weight
- Edge cases: single codon amino acids (M, W), stop codons
- Cross-organism pairs (E. coli → Human, Yeast → CHO, etc.)
- Pipeline integration (strategy='harmonize')
"""

import pytest

from biocompiler.optimizer.codon_harmonization import (
    compute_rca,
    harmonize_codons,
    harmonize_with_cai_fallback,
    compute_harmonization_score,
)
from biocompiler.organisms import (
    E_COLI_CODON_USAGE,
    HUMAN_CODON_USAGE,
    YEAST_CODON_USAGE,
    CHO_CODON_USAGE,
    MOUSE_CODON_USAGE,
    CODON_USAGE_TABLES,
    resolve_organism,
)
from biocompiler.translation import translate, compute_cai
from biocompiler.type_system import CODON_TABLE, AA_TO_CODONS


# ────────────────────────────────────────────────────────────
# Fixtures
# ────────────────────────────────────────────────────────────

@pytest.fixture
def ecoli_rca():
    """RCA values for E. coli."""
    return compute_rca(E_COLI_CODON_USAGE)


@pytest.fixture
def human_rca():
    """RCA values for Human."""
    return compute_rca(HUMAN_CODON_USAGE)


@pytest.fixture
def gfp_protein():
    """Green fluorescent protein (GFP) amino acid sequence (shortened)."""
    return "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"


@pytest.fixture
def insulin_protein():
    """Insulin B chain (short protein for quick tests)."""
    return "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"


# ────────────────────────────────────────────────────────────
# Test RCA Computation
# ────────────────────────────────────────────────────────────

class TestComputeRCA:
    """Tests for compute_rca()."""

    def test_rca_values_in_range(self, ecoli_rca):
        """All RCA values should be in [0.0, 1.0]."""
        for codon, rca_val in ecoli_rca.items():
            assert 0.0 <= rca_val <= 1.0, (
                f"RCA for {codon} = {rca_val} is out of [0, 1]"
            )

    def test_most_frequent_codon_has_rca_one(self, ecoli_rca):
        """Most frequent codon per amino acid should have RCA = 1.0."""
        # Group by amino acid and check
        aa_max_freq = {}
        for codon, (aa, _frac, freq, _count) in E_COLI_CODON_USAGE.items():
            if aa == "*":
                continue
            if aa not in aa_max_freq or freq > aa_max_freq[aa][1]:
                aa_max_freq[aa] = (codon, freq)

        for aa, (codon, _freq) in aa_max_freq.items():
            assert ecoli_rca[codon] == 1.0, (
                f"Most frequent codon {codon} for {aa} should have RCA=1.0, "
                f"got {ecoli_rca[codon]}"
            )

    def test_known_ecoli_rca_values(self, ecoli_rca):
        """Spot-check specific E. coli RCA values."""
        # Leucine: CTG is most frequent (per_thousand=54.8)
        # TTA per_thousand=7.1 → RCA = 7.1/54.8 ≈ 0.1296
        assert abs(ecoli_rca["TTA"] - 7.1 / 54.8) < 0.01, (
            f"RCA for TTA should be ~{7.1/54.8:.4f}, got {ecoli_rca['TTA']}"
        )

        # Phenylalanine: TTC is most frequent (per_thousand=22.1)
        # TTT per_thousand=15.2 → RCA = 15.2/22.1 ≈ 0.6878
        assert abs(ecoli_rca["TTT"] - 15.2 / 22.1) < 0.01, (
            f"RCA for TTT should be ~{15.2/22.1:.4f}, got {ecoli_rca['TTT']}"
        )

    def test_rca_excludes_stop_codons(self, ecoli_rca):
        """Stop codons should not appear in RCA output."""
        assert "TAA" not in ecoli_rca
        assert "TAG" not in ecoli_rca
        assert "TGA" not in ecoli_rca

    def test_single_codon_aas_have_rca_one(self, ecoli_rca):
        """Amino acids with only one codon (M, W) should have RCA = 1.0."""
        assert ecoli_rca["ATG"] == 1.0  # Methionine
        assert ecoli_rca["TGG"] == 1.0  # Tryptophan

    def test_human_rca_computation(self, human_rca):
        """Human RCA values should also be valid."""
        for codon, rca_val in human_rca.items():
            assert 0.0 <= rca_val <= 1.0
        # ATG should always be 1.0
        assert human_rca["ATG"] == 1.0

    def test_rca_with_custom_usage_table(self):
        """RCA computation should work with custom codon usage tables."""
        custom_usage = {
            "ATG": ("M", 1.0, 25.0, 250),
            "TTT": ("F", 0.3, 10.0, 100),
            "TTC": ("F", 0.7, 23.3, 233),
        }
        rca = compute_rca(custom_usage)
        assert rca["ATG"] == 1.0
        assert rca["TTC"] == 1.0  # Most frequent for F
        assert abs(rca["TTT"] - 10.0 / 23.3) < 0.01


# ────────────────────────────────────────────────────────────
# Test Codon Harmonization
# ────────────────────────────────────────────────────────────

class TestHarmonizeCodons:
    """Tests for harmonize_codons()."""

    def test_correct_length(self, gfp_protein):
        """Harmonized sequence length should be 3× protein length."""
        dna = harmonize_codons(gfp_protein, "Escherichia_coli", "Homo_sapiens")
        assert len(dna) == len(gfp_protein) * 3

    def test_translates_back(self, gfp_protein):
        """Harmonized DNA should translate back to the original protein."""
        dna = harmonize_codons(gfp_protein, "Escherichia_coli", "Homo_sapiens")
        translated = translate(dna, to_stop=True)
        assert translated == gfp_protein

    def test_only_valid_codons(self, gfp_protein):
        """All triplets in harmonized DNA should be valid codons."""
        dna = harmonize_codons(gfp_protein, "Escherichia_coli", "Homo_sapiens")
        for i in range(0, len(dna), 3):
            codon = dna[i:i + 3]
            assert codon in CODON_TABLE, f"Invalid codon: {codon}"

    def test_no_stop_codons_in_sequence(self, gfp_protein):
        """Harmonized DNA should not contain internal stop codons."""
        dna = harmonize_codons(gfp_protein, "Escherichia_coli", "Homo_sapiens")
        for i in range(0, len(dna), 3):
            codon = dna[i:i + 3]
            assert CODON_TABLE[codon] != "*", f"Internal stop codon at position {i}"

    def test_single_codon_amino_acids(self):
        """M and W should always use their only codons (ATG, TGG)."""
        protein = "MWMWMW"
        dna = harmonize_codons(protein, "Escherichia_coli", "Homo_sapiens")
        assert dna == "ATGTGGATGTGGATGTGG"

    def test_self_harmonization_uses_optimal_codons(self):
        """Self-harmonization (same source and target) should use
        high-frequency codons."""
        protein = "FFFF"  # Phenylalanine - 2 codons
        dna = harmonize_codons(protein, "Escherichia_coli", "Escherichia_coli")
        # Self-harmonization should produce all TTC (most frequent in E. coli)
        assert dna == "TTCTTCTTCTTC"

    def test_cross_organism_ecoli_to_human(self, insulin_protein):
        """E. coli → Human harmonization should produce valid DNA."""
        dna = harmonize_codons(
            insulin_protein,
            "Escherichia_coli",
            "Homo_sapiens",
        )
        assert len(dna) == len(insulin_protein) * 3
        assert translate(dna, to_stop=True) == insulin_protein

    def test_cross_organism_yeast_to_cho(self, insulin_protein):
        """Yeast → CHO harmonization should produce valid DNA."""
        dna = harmonize_codons(
            insulin_protein,
            "Saccharomyces_cerevisiae",
            "CHO_K1",
        )
        assert len(dna) == len(insulin_protein) * 3
        assert translate(dna, to_stop=True) == insulin_protein

    def test_short_protein(self):
        """Test with a very short protein."""
        protein = "MVHLTPEEK"
        dna = harmonize_codons(protein, "Escherichia_coli", "Homo_sapiens")
        assert len(dna) == 27
        assert translate(dna, to_stop=True) == protein

    def test_organism_alias_resolution(self, insulin_protein):
        """Short organism names should be resolved correctly."""
        # Both should work with short names
        dna1 = harmonize_codons(insulin_protein, "ecoli", "human")
        dna2 = harmonize_codons(insulin_protein, "Escherichia_coli", "Homo_sapiens")
        assert dna1 == dna2

    def test_invalid_organism_raises(self):
        """Invalid organism name should raise ValueError."""
        with pytest.raises(ValueError, match="No codon usage data"):
            harmonize_codons("MWW", "Nonexistent_organism", "Homo_sapiens")

    def test_custom_codon_usage_tables(self):
        """Custom codon usage tables should override lookup."""
        # Create custom source usage where TTT is very frequent for F
        custom_source = {
            "TTT": ("F", 0.9, 30.0, 300),
            "TTC": ("F", 0.1, 3.3, 33),
            "ATG": ("M", 1.0, 25.0, 250),
        }
        custom_target = {
            "TTT": ("F", 0.1, 5.0, 50),
            "TTC": ("F", 0.9, 45.0, 450),
            "ATG": ("M", 1.0, 25.0, 250),
        }
        protein = "MF"
        dna = harmonize_codons(
            protein,
            "Escherichia_coli",  # Will be overridden
            "Homo_sapiens",      # Will be overridden
            codon_usage_source=custom_source,
            codon_usage_target=custom_target,
        )
        # Source has TTT at RCA=1.0 (most frequent)
        # Target has TTT at RCA=5/45 ≈ 0.111, TTC at RCA=1.0
        # For F: source weighted RCA is high (~0.9), target TTT is low
        # So harmonization should pick TTC (closer to high RCA)
        # unless the weighted average calculation selects TTT
        assert translate(dna, to_stop=True) == "MF"


# ────────────────────────────────────────────────────────────
# Test Harmonize with CAI Fallback
# ────────────────────────────────────────────────────────────

class TestHarmonizeWithCAIFallback:
    """Tests for harmonize_with_cai_fallback()."""

    def test_pure_harmonization_same_as_harmonize_codons(self, insulin_protein):
        """cai_weight=0.0 should produce same result as harmonize_codons."""
        dna_harmonize = harmonize_codons(
            insulin_protein, "Escherichia_coli", "Homo_sapiens"
        )
        dna_fallback = harmonize_with_cai_fallback(
            insulin_protein, "Escherichia_coli", "Homo_sapiens",
            cai_weight=0.0,
        )
        assert dna_harmonize == dna_fallback

    def test_pure_cai_uses_optimal_codons(self):
        """cai_weight=1.0 should prefer CAI-optimal codons."""
        protein = "FFFF"  # Simple test with Phe
        dna = harmonize_with_cai_fallback(
            protein, "Escherichia_coli", "Escherichia_coli",
            cai_weight=1.0,
        )
        # With cai_weight=1.0, should select highest CAI codon
        # For E. coli F: TTC has highest adaptiveness
        assert all(dna[i:i+3] == "TTC" for i in range(0, len(dna), 3))

    def test_blended_weight_produces_valid_dna(self, gfp_protein):
        """cai_weight=0.5 should produce valid DNA."""
        dna = harmonize_with_cai_fallback(
            gfp_protein, "Escherichia_coli", "Homo_sapiens",
            cai_weight=0.5,
        )
        assert len(dna) == len(gfp_protein) * 3
        assert translate(dna, to_stop=True) == gfp_protein

    def test_invalid_cai_weight_raises(self):
        """cai_weight outside [0, 1] should raise ValueError."""
        with pytest.raises(ValueError, match="cai_weight must be in"):
            harmonize_with_cai_fallback("MW", "ecoli", "human", cai_weight=-0.1)
        with pytest.raises(ValueError, match="cai_weight must be in"):
            harmonize_with_cai_fallback("MW", "ecoli", "human", cai_weight=1.5)

    def test_different_weights_produce_different_sequences(self, insulin_protein):
        """Different cai_weights should generally produce different sequences."""
        dna_low = harmonize_with_cai_fallback(
            insulin_protein, "Escherichia_coli", "Homo_sapiens",
            cai_weight=0.0,
        )
        dna_mid = harmonize_with_cai_fallback(
            insulin_protein, "Escherichia_coli", "Homo_sapiens",
            cai_weight=0.5,
        )
        dna_high = harmonize_with_cai_fallback(
            insulin_protein, "Escherichia_coli", "Homo_sapiens",
            cai_weight=1.0,
        )
        # At least some of these should differ
        # (They might not all differ for every protein, but for insulin
        # with multiple codon choices, they should)
        sequences = {dna_low, dna_mid, dna_high}
        assert len(sequences) >= 2, (
            "Different cai_weights should produce at least 2 different sequences"
        )

    def test_higher_cai_weight_gives_higher_cai(self, insulin_protein):
        """Higher cai_weight should produce sequences with higher CAI."""
        dna_low = harmonize_with_cai_fallback(
            insulin_protein, "Escherichia_coli", "Homo_sapiens",
            cai_weight=0.0,
        )
        dna_high = harmonize_with_cai_fallback(
            insulin_protein, "Escherichia_coli", "Homo_sapiens",
            cai_weight=1.0,
        )
        cai_low = compute_cai(dna_low, organism="Homo_sapiens")
        cai_high = compute_cai(dna_high, organism="Homo_sapiens")
        assert cai_high >= cai_low, (
            f"CAI with weight=1.0 ({cai_high:.4f}) should be >= "
            f"CAI with weight=0.0 ({cai_low:.4f})"
        )


# ────────────────────────────────────────────────────────────
# Test Harmonization Score
# ────────────────────────────────────────────────────────────

class TestComputeHarmonizationScore:
    """Tests for compute_harmonization_score()."""

    def test_perfect_self_harmonization(self):
        """Self-harmonized (optimal codon) sequence should score high."""
        protein = "MVHLTPEEKSAVTALWGKVNV"
        # Use all preferred codons for the target organism
        dna = harmonize_codons(
            protein, "Escherichia_coli", "Escherichia_coli"
        )
        score = compute_harmonization_score(
            dna, "Escherichia_coli", "Escherichia_coli"
        )
        # Self-harmonization with optimal codons should score very high
        assert score > 0.8, f"Self-harmonization score should be > 0.8, got {score}"

    def test_harmonized_better_than_cai_maximized(self, insulin_protein):
        """Harmonized sequence should score higher than CAI-maximized sequence."""
        # Generate harmonized sequence
        dna_harmonized = harmonize_codons(
            insulin_protein, "Escherichia_coli", "Homo_sapiens"
        )

        # Generate CAI-maximized sequence (all optimal codons for target)
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        from biocompiler.type_system import AA_TO_CODONS as _AA2C

        cai_weights = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
        cai_codons = []
        for aa in insulin_protein:
            codons = _AA2C.get(aa, [])
            if codons:
                best = max(codons, key=lambda c: cai_weights.get(c, 0.0))
                cai_codons.append(best)
            else:
                cai_codons.append("ATG")
        dna_cai = "".join(cai_codons)

        score_harmonized = compute_harmonization_score(
            dna_harmonized, "Escherichia_coli", "Homo_sapiens"
        )
        score_cai = compute_harmonization_score(
            dna_cai, "Escherichia_coli", "Homo_sapiens"
        )

        assert score_harmonized >= score_cai, (
            f"Harmonized score ({score_harmonized:.4f}) should be >= "
            f"CAI-maximized score ({score_cai:.4f})"
        )

    def test_score_in_range(self, gfp_protein):
        """Harmonization score should be in [0.0, 1.0]."""
        dna = harmonize_codons(gfp_protein, "Escherichia_coli", "Homo_sapiens")
        score = compute_harmonization_score(
            dna, "Escherichia_coli", "Homo_sapiens"
        )
        assert 0.0 <= score <= 1.0

    def test_empty_sequence_returns_zero(self):
        """Empty sequence should return 0.0."""
        score = compute_harmonization_score("", "Escherichia_coli", "Homo_sapiens")
        assert score == 0.0

    def test_short_sequence(self):
        """Very short sequence should still produce a valid score."""
        dna = "ATGTGG"  # MW
        score = compute_harmonization_score(
            dna, "Escherichia_coli", "Homo_sapiens"
        )
        assert 0.0 <= score <= 1.0

    def test_invalid_organism_returns_zero(self):
        """Invalid organism should return 0.0 (not raise)."""
        score = compute_harmonization_score(
            "ATGTGG", "Nonexistent_organism", "Homo_sapiens"
        )
        assert score == 0.0


# ────────────────────────────────────────────────────────────
# Test Edge Cases
# ────────────────────────────────────────────────────────────

class TestEdgeCases:
    """Edge case tests for codon harmonization."""

    def test_single_amino_acid_protein(self):
        """Protein with single amino acid should work."""
        protein = "M"
        dna = harmonize_codons(protein, "Escherichia_coli", "Homo_sapiens")
        assert dna == "ATG"

    def test_all_same_amino_acid(self):
        """Protein with all same amino acid should work."""
        protein = "LLLLLL"
        dna = harmonize_codons(protein, "Escherichia_coli", "Homo_sapiens")
        assert len(dna) == 18
        translated = translate(dna, to_stop=True)
        assert translated == protein

    def test_methionine_only_codon(self):
        """Methionine (ATG) should always be the same regardless of organism."""
        protein = "MMMMM"
        dna = harmonize_codons(protein, "Escherichia_coli", "Homo_sapiens")
        assert dna == "ATGATGATGATGATG"

    def test_tryptophan_only_codon(self):
        """Tryptophan (TGG) should always be the same regardless of organism."""
        protein = "WWWWW"
        dna = harmonize_codons(protein, "Escherichia_coli", "Homo_sapiens")
        assert dna == "TGGTGGTGGTGGTGG"

    def test_all_20_amino_acids(self):
        """Protein with all 20 standard amino acids should work."""
        protein = "ACDEFGHIKLMNPQRSTVWY"
        dna = harmonize_codons(protein, "Escherichia_coli", "Homo_sapiens")
        assert len(dna) == 60
        assert translate(dna, to_stop=True) == protein

    def test_large_protein(self, gfp_protein):
        """Larger protein (GFP ~240 AA) should harmonize correctly."""
        dna = harmonize_codons(gfp_protein, "Escherichia_coli", "Homo_sapiens")
        assert translate(dna, to_stop=True) == gfp_protein

    def test_same_source_and_target(self, insulin_protein):
        """Same source and target should produce valid DNA."""
        dna = harmonize_codons(
            insulin_protein, "Escherichia_coli", "Escherichia_coli"
        )
        assert translate(dna, to_stop=True) == insulin_protein


# ────────────────────────────────────────────────────────────
# Test Cross-Organism Pairs
# ────────────────────────────────────────────────────────────

class TestCrossOrganismPairs:
    """Test harmonization across multiple organism pairs."""

    @pytest.mark.parametrize("source,target", [
        ("Escherichia_coli", "Homo_sapiens"),
        ("Homo_sapiens", "Escherichia_coli"),
        ("Saccharomyces_cerevisiae", "CHO_K1"),
        ("CHO_K1", "Saccharomyces_cerevisiae"),
        ("Mus_musculus", "Homo_sapiens"),
        ("Escherichia_coli", "Saccharomyces_cerevisiae"),
    ])
    def test_cross_organism_valid_dna(self, source, target, insulin_protein):
        """Cross-organism harmonization should produce valid DNA."""
        dna = harmonize_codons(insulin_protein, source, target)
        assert len(dna) == len(insulin_protein) * 3
        assert translate(dna, to_stop=True) == insulin_protein

    @pytest.mark.parametrize("source,target", [
        ("ecoli", "human"),
        ("human", "ecoli"),
        ("yeast", "cho"),
        ("mouse", "human"),
    ])
    def test_alias_resolution(self, source, target, insulin_protein):
        """Short organism names should work correctly."""
        dna = harmonize_codons(insulin_protein, source, target)
        assert len(dna) == len(insulin_protein) * 3
        assert translate(dna, to_stop=True) == insulin_protein


# ────────────────────────────────────────────────────────────
# Test Pipeline Integration
# ────────────────────────────────────────────────────────────

class TestPipelineIntegration:
    """Test the harmonize strategy through optimize_sequence()."""

    def test_harmonize_strategy_basic(self):
        """strategy='harmonize' should produce valid result."""
        from biocompiler import optimize_sequence
        result = optimize_sequence(
            "MVHLTPEEK",
            organism="Homo_sapiens",
            strategy="harmonize",
            source_organism="Escherichia_coli",
        )
        assert result.sequence
        assert len(result.sequence) == 27
        assert result.cai > 0.0
        assert result.gc_content > 0.0
        assert result.protein == "MVHLTPEEK"

    def test_harmonize_strategy_with_cai_weight(self):
        """harmonization_cai_weight should affect the result."""
        from biocompiler import optimize_sequence
        result_pure = optimize_sequence(
            "MVHLTPEEK",
            organism="Homo_sapiens",
            strategy="harmonize",
            source_organism="Escherichia_coli",
            harmonization_cai_weight=0.0,
        )
        result_blended = optimize_sequence(
            "MVHLTPEEK",
            organism="Homo_sapiens",
            strategy="harmonize",
            source_organism="Escherichia_coli",
            harmonization_cai_weight=0.5,
        )
        # Both should be valid
        assert result_pure.sequence
        assert result_blended.sequence
        assert translate(result_pure.sequence, to_stop=True) == "MVHLTPEEK"
        assert translate(result_blended.sequence, to_stop=True) == "MVHLTPEEK"

    def test_harmonize_strategy_defaults_source_to_target(self):
        """Without source_organism, should default to target organism."""
        from biocompiler import optimize_sequence
        result = optimize_sequence(
            "MVHLTPEEK",
            organism="Homo_sapiens",
            strategy="harmonize",
        )
        assert result.sequence
        assert translate(result.sequence, to_stop=True) == "MVHLTPEEK"

    def test_harmonize_strategy_produces_harmonization_score(self):
        """Result should include harmonization score as objective_score."""
        from biocompiler import optimize_sequence
        result = optimize_sequence(
            "MVHLTPEEK",
            organism="Homo_sapiens",
            strategy="harmonize",
            source_organism="Escherichia_coli",
        )
        # objective_score should be the harmonization score
        assert result.objective_score is not None
        assert 0.0 <= result.objective_score <= 1.0

    def test_harmonize_strategy_with_prokaryote_target(self):
        """Harmonize strategy should work with prokaryotic targets."""
        from biocompiler import optimize_sequence
        result = optimize_sequence(
            "MVHLTPEEK",
            organism="Escherichia_coli",
            strategy="harmonize",
            source_organism="Homo_sapiens",
        )
        assert result.sequence
        assert translate(result.sequence, to_stop=True) == "MVHLTPEEK"


# ────────────────────────────────────────────────────────────
# Test Objective Integration
# ────────────────────────────────────────────────────────────

class TestObjectiveIntegration:
    """Test the harmonization objective in OBJECTIVE_REGISTRY."""

    def test_harmonization_in_registry(self):
        """'harmonization' should be in OBJECTIVE_REGISTRY."""
        from biocompiler.optimizer.objectives import OBJECTIVE_REGISTRY
        assert "harmonization" in OBJECTIVE_REGISTRY

    def test_resolve_harmonization_objective(self):
        """resolve_objective('harmonization') should return a callable."""
        from biocompiler.optimizer.objectives import resolve_objective
        obj = resolve_objective("harmonization")
        assert callable(obj)

    def test_harmonization_objective_returns_float(self):
        """harmonization_objective should return a float."""
        from biocompiler.optimizer.objectives import harmonization_objective
        dna = "ATGTGGCCT"  # MWP
        score = harmonization_objective(dna, "MWP", "Escherichia_coli")
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_make_harmonization_objective(self):
        """make_harmonization_objective should create a configured objective."""
        from biocompiler.optimizer.objectives import make_harmonization_objective
        obj = make_harmonization_objective("Escherichia_coli")
        assert callable(obj)
        dna = "ATGTGGCCT"
        score = obj(dna, "MWP", "Homo_sapiens")
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_harmonization_objective_with_source(self):
        """harmonization_objective with source_organism should compute
        cross-organism score."""
        from biocompiler.optimizer.objectives import harmonization_objective
        protein = "MVHLTPEEK"
        dna_harm = harmonize_codons(protein, "Escherichia_coli", "Homo_sapiens")
        dna_cai = harmonize_with_cai_fallback(
            protein, "Escherichia_coli", "Homo_sapiens", cai_weight=1.0
        )
        score_harm = harmonization_objective(
            dna_harm, protein, "Homo_sapiens",
            source_organism="Escherichia_coli",
        )
        score_cai = harmonization_objective(
            dna_cai, protein, "Homo_sapiens",
            source_organism="Escherichia_coli",
        )
        # The harmonized sequence should have a better (or equal) score
        assert score_harm >= score_cai, (
            f"Harmonized score ({score_harm:.4f}) should be >= "
            f"CAI-maximized score ({score_cai:.4f})"
        )
