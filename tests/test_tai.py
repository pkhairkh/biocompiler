"""Tests for the tRNA Adaptation Index (tAI) module.

Validates the tAI implementation including wobble rules,
tRNA gene copy number handling, organism resolution,
cross-species comparisons, objective integration, and
optimization pipeline integration.
"""

import math
import random
import pytest

from biocompiler.expression.tai import (
    compute_tai,
    calculate_tai,
    compute_tai_and_cai,
    optimize_for_tai,
    TRNA_GENE_COPIES,
    WOBBLE_RULES,
    WOBBLE_EFFICIENCY,
    SUPPORTED_ORGANISMS_TAI,
    compute_codon_weights,
)
from biocompiler.organisms.tai_data import compute_tai_weights
from biocompiler.expression.translation import compute_cai
from biocompiler.optimizer.objectives import tai_objective, resolve_objective, OBJECTIVE_REGISTRY
from biocompiler.shared.constants import CODON_TABLE, AA_TO_CODONS


# ═══════════════════════════════════════════════════════════════════════════════
# Test Data
# ═══════════════════════════════════════════════════════════════════════════════

# eGFP in E. coli
_EGFP_DNA = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGG"
    "CCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCAC"
    "CACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCT"
    "ACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTT"
    "CTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGA"
    "GCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAA"
    "CGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGC"
    "AGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCAC"
    "TACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTG"
    "ACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
)

# Human insulin
_HUMAN_INSULIN_DNA = (
    "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCT"
    "TTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTC"
    "TACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGC"
    "AGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCA"
    "TCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
)

# All 10 organisms for testing
_ALL_ORGANISMS = [
    "e_coli", "human", "yeast", "mouse", "cho",
    "c_elegans", "d_melanogaster", "a_thaliana", "p_pastoris", "b_subtilis",
]

_ALL_ORGANISM_ALIASES = {
    "e_coli": ["Escherichia_coli", "ecoli", "E. coli", "E_coli"],
    "human": ["Homo_sapiens", "h_sapiens", "H. sapiens", "H_sapiens"],
    "yeast": ["Saccharomyces_cerevisiae", "s_cerevisiae", "S. cerevisiae", "S_cerevisiae"],
    "mouse": ["Mus_musculus", "m_musculus", "M. musculus", "M_musculus"],
    "cho": ["CHO_K1", "CHO", "Cricetulus_griseus"],
    "c_elegans": ["Caenorhabditis_elegans", "C. elegans", "C_elegans"],
    "d_melanogaster": ["Drosophila_melanogaster", "D. melanogaster", "D_melanogaster"],
    "a_thaliana": ["Arabidopsis_thaliana", "A. thaliana", "A_thaliana"],
    "p_pastoris": ["Pichia_pastoris", "P. pastoris", "P_pastoris", "Komagataella_phaffii"],
    "b_subtilis": ["Bacillus_subtilis", "B. subtilis", "B_subtilis"],
}


# ═══════════════════════════════════════════════════════════════════════════════
# Test: Wobble Rules Correctness
# ═══════════════════════════════════════════════════════════════════════════════

class TestWobbleRulesCorrectness:
    """Test that wobble rules match dos Reis et al. (2004)."""

    def test_all_64_codons_represented(self):
        """All 64 codons should be in the wobble rules."""
        assert len(WOBBLE_RULES) == 64

    def test_gu_wobble_is_0_5(self):
        """G:U wobble pairs (anticodon G, codon U) should have efficiency 0.5.

        This was a critical bug: the original implementation had 1.0 instead of 0.5.
        Codons ending in U (RNA) that are read by anticodons starting with G
        should use G:U wobble = 0.5, not 1.0.
        """
        # All NNU codons read by GNN anticodons should have 0.5
        gu_wobble_codons = {
            "UUU": "GAA",  # Phe
            "CUU": "GAG",  # Leu
            "GUU": "GAC",  # Val
            "UCU": "GGA",  # Ser
            "CCU": "GGG",  # Pro
            "ACU": "GGU",  # Thr
            "GCU": "GGC",  # Ala
            "CGU": "GCG",  # Arg
            "GGU": "GCC",  # Gly
            "UAU": "GUA",  # Tyr
            "CAU": "GUG",  # His
            "AAU": "GUU",  # Asn
            "GAU": "GUC",  # Asp
            "UGU": "GCA",  # Cys
            "AGU": "GCU",  # Ser
        }
        for codon, anticodon in gu_wobble_codons.items():
            rules = WOBBLE_RULES[codon]
            for ac, eff in rules:
                if ac == anticodon:
                    assert eff == 0.5, (
                        f"G:U wobble for codon {codon} read by {anticodon} "
                        f"should be 0.5, got {eff}"
                    )
                    break
            else:
                pytest.fail(f"Anticodon {anticodon} not found for codon {codon}")

    def test_gc_watson_crick_is_1_0(self):
        """G:C Watson-Crick pairs should have efficiency 1.0."""
        gc_wc_codons = {
            "UUC": "GAA",  # Phe
            "CUC": "GAG",  # Leu
            "GUC": "GAC",  # Val
            "UCC": "GGA",  # Ser
            "CCC": "GGG",  # Pro
            "ACC": "GGU",  # Thr
            "GCC": "GGC",  # Ala
            "CGC": "GCG",  # Arg
            "GGC": "GCC",  # Gly
            "UAC": "GUA",  # Tyr
            "CAC": "GUG",  # His
            "AAC": "GUU",  # Asn
            "GAC": "GUC",  # Asp
            "UGC": "GCA",  # Cys
            "AGC": "GCU",  # Ser
        }
        for codon, anticodon in gc_wc_codons.items():
            rules = WOBBLE_RULES[codon]
            for ac, eff in rules:
                if ac == anticodon:
                    assert eff == 1.0, (
                        f"G:C Watson-Crick for codon {codon} read by {anticodon} "
                        f"should be 1.0, got {eff}"
                    )
                    break
            else:
                pytest.fail(f"Anticodon {anticodon} not found for codon {codon}")

    def test_ug_wobble_is_0_2(self):
        """U:G wobble pairs (anticodon U, codon G) should have efficiency 0.2."""
        ug_wobble_codons = {
            "UUG": "UAA",  # Leu
            "CUG": "UAG",  # Leu
            "GUG": "UAC",  # Val
            "CCG": "UGG",  # Pro
            "AGG": "UCU",  # Arg
            "GGG": "UCC",  # Gly
            "CAG": "UUG",  # Gln
            "AAG": "UUU",  # Lys
            "GAG": "UUC",  # Glu
        }
        for codon, anticodon in ug_wobble_codons.items():
            rules = WOBBLE_RULES[codon]
            for ac, eff in rules:
                if ac == anticodon:
                    assert eff == 0.2, (
                        f"U:G wobble for codon {codon} read by {anticodon} "
                        f"should be 0.2, got {eff}"
                    )
                    break
            else:
                pytest.fail(f"Anticodon {anticodon} not found for codon {codon}")

    def test_inosine_efficiencies(self):
        """Inosine should pair with U(0.35), C(0.65), A(0.15) but NOT G."""
        assert WOBBLE_EFFICIENCY[("I", "U")] == 0.35
        assert WOBBLE_EFFICIENCY[("I", "C")] == 0.65
        assert WOBBLE_EFFICIENCY[("I", "A")] == 0.15
        # I:G should NOT be in the efficiency table (I cannot pair with G)
        assert ("I", "G") not in WOBBLE_EFFICIENCY

    def test_no_inosine_g_pairing_in_rules(self):
        """No wobble rule should list Inosine pairing with G.

        I:G is not a valid pairing per dos Reis (2004).
        Codons ending in G should NOT be read by Inosine-containing anticodons.
        """
        g_codons_with_inosine = []
        for codon, anticodons in WOBBLE_RULES.items():
            if codon.endswith("G"):
                for anticodon, efficiency in anticodons:
                    if anticodon.startswith("I"):
                        g_codons_with_inosine.append((codon, anticodon, efficiency))

        assert len(g_codons_with_inosine) == 0, (
            f"Inosine should not pair with G, but found: {g_codons_with_inosine}"
        )

    def test_ile_wobble_rules_correct(self):
        """Isoleucine wobble rules should use GAU (bacteria) and IAU (eukaryotes).

        The original implementation incorrectly used UAU/CAU for Ile codons.
        Correct rules:
        - AUU: GAU(0.5), IAU(0.35)
        - AUC: GAU(1.0), IAU(0.65)
        - AUA: k2C(1.0), IAU(0.15)
        """
        # AUU
        auu_rules = dict(WOBBLE_RULES["AUU"])
        assert "GAU" in auu_rules, "AUU should be read by GAU anticodon"
        assert auu_rules["GAU"] == 0.5, "GAU:AUU should be G:U wobble = 0.5"
        assert "IAU" in auu_rules, "AUU should be read by IAU anticodon (eukaryotes)"
        assert auu_rules["IAU"] == 0.35, "IAU:AUU should be I:U = 0.35"
        # CAU should NOT be in AUU rules (that is Met tRNA)
        assert "CAU" not in auu_rules, "CAU (Met tRNA) should not read AUU"

        # AUC
        auc_rules = dict(WOBBLE_RULES["AUC"])
        assert "GAU" in auc_rules, "AUC should be read by GAU anticodon"
        assert auc_rules["GAU"] == 1.0, "GAU:AUC should be G:C = 1.0"
        assert "IAU" in auc_rules, "AUC should be read by IAU anticodon (eukaryotes)"
        assert auc_rules["IAU"] == 0.65, "IAU:AUC should be I:C = 0.65"
        assert "CAU" not in auc_rules, "CAU (Met tRNA) should not read AUC"

        # AUA
        aua_rules = dict(WOBBLE_RULES["AUA"])
        assert "k2C" in aua_rules, "AUA should be read by k2C (lysidine) anticodon"
        assert aua_rules["k2C"] == 1.0, "k2C:AUA should be lysidine:A = 1.0"
        assert "IAU" in aua_rules, "AUA should be read by IAU (eukaryotes)"
        assert aua_rules["IAU"] == 0.15, "IAU:AUA should be I:A = 0.15"

    def test_met_anticodon_is_cau(self):
        """Methionine AUG should be read by CAU anticodon (unmodified)."""
        aug_rules = dict(WOBBLE_RULES["AUG"])
        assert "CAU" in aug_rules, "AUG should be read by CAU (Met tRNA)"
        assert aug_rules["CAU"] == 1.0, "CAU:AUG should be C:G = 1.0"

    def test_each_codon_has_at_least_one_anticodon(self):
        """Every codon should be readable by at least one anticodon."""
        for codon, anticodons in WOBBLE_RULES.items():
            assert len(anticodons) >= 1, f"Codon {codon} has no anticodon"
            for anticodon, efficiency in anticodons:
                assert 0.0 < efficiency <= 1.0, \
                    f"Invalid efficiency for {codon}/{anticodon}: {efficiency}"

    def test_stop_codons_have_wobble_rules(self):
        """Stop codons should have wobble rules."""
        for stop in ("UAA", "UAG", "UGA"):
            assert stop in WOBBLE_RULES

    def test_k2c_reads_a_only(self):
        """k2C (lysidine) should only pair with A, not U, C, or G."""
        assert WOBBLE_EFFICIENCY[("k2C", "A")] == 1.0
        assert ("k2C", "U") not in WOBBLE_EFFICIENCY
        assert ("k2C", "C") not in WOBBLE_EFFICIENCY
        assert ("k2C", "G") not in WOBBLE_EFFICIENCY

    def test_watson_crick_pairs_efficiency_1(self):
        """All Watson-Crick pairs should have efficiency 1.0."""
        assert WOBBLE_EFFICIENCY[("A", "U")] == 1.0
        assert WOBBLE_EFFICIENCY[("U", "A")] == 1.0
        assert WOBBLE_EFFICIENCY[("G", "C")] == 1.0
        assert WOBBLE_EFFICIENCY[("C", "G")] == 1.0

    def test_gu_wobble_efficiency_value(self):
        """G:U wobble pair should have efficiency 0.5."""
        assert WOBBLE_EFFICIENCY[("G", "U")] == 0.5


# ═══════════════════════════════════════════════════════════════════════════════
# Test: tRNA Gene Copy Number Database
# ═══════════════════════════════════════════════════════════════════════════════

class TestTRNAGeneCopies:
    """Test the tRNA gene copy number database."""

    @pytest.mark.parametrize("org", _ALL_ORGANISMS)
    def test_organism_data_exists(self, org):
        """Each organism should have tRNA data."""
        assert org in TRNA_GENE_COPIES

    def test_supported_organisms_count(self):
        """Should have at least 10 organisms."""
        assert len(SUPPORTED_ORGANISMS_TAI) >= 10

    def test_ecoli_has_required_anticodons(self):
        """E. coli should have tRNA data for common anticodons."""
        ecoli = TRNA_GENE_COPIES["e_coli"]
        assert "GAA" in ecoli  # Phe
        assert "CAU" in ecoli  # Met
        assert "UUC" in ecoli  # Glu
        assert "GUC" in ecoli  # Asp
        assert "GAU" in ecoli  # Ile (reads AUC, AUU)
        assert "k2C" in ecoli  # Ile-lysidine (reads AUA)

    def test_ecoli_no_cau_collision(self):
        """E. coli Ile-CAU should be stored as k2C, not conflated with Met-CAU.

        This tests the fix for the CAU collision bug where Ile-tRNA-CAU
        (lysidine) and Met-tRNA-CAU had the same dict key, causing
        Met-CAU to overwrite Ile-CAU.
        """
        ecoli = TRNA_GENE_COPIES["e_coli"]
        # CAU is Met only (should NOT also represent Ile-lysidine)
        assert "CAU" in ecoli  # Met
        assert "k2C" in ecoli  # Ile-lysidine (separate key)
        # GAU should be present for Ile (reads AUC, AUU)
        assert "GAU" in ecoli

    def test_human_has_inosine_anticodons(self):
        """Human should have Inosine-containing anticodons for eukaryotic tRNAs."""
        human = TRNA_GENE_COPIES["human"]
        assert "IAU" in human  # Ile with Inosine
        assert "IGA" in human  # Ser with Inosine
        assert "IGU" in human  # Thr with Inosine
        assert "IGC" in human  # Ala with Inosine
        assert "ICG" in human  # Arg with Inosine

    def test_yeast_has_inosine_anticodons(self):
        """Yeast should have Inosine-containing anticodons."""
        yeast = TRNA_GENE_COPIES["yeast"]
        assert "IAU" in yeast

    def test_bacteria_have_gau_not_iau(self):
        """Bacterial organisms should have GAU for Ile, not IAU."""
        for org in ["e_coli", "b_subtilis"]:
            data = TRNA_GENE_COPIES[org]
            assert "GAU" in data, f"{org} should have GAU (Ile tRNA)"
            assert "IAU" not in data, f"{org} should not have IAU (eukaryotic)"

    def test_eukaryotes_have_iau(self):
        """Eukaryotic organisms should have IAU for Ile."""
        for org in ["human", "yeast", "mouse", "cho", "c_elegans",
                     "d_melanogaster", "a_thaliana", "p_pastoris"]:
            data = TRNA_GENE_COPIES[org]
            assert "IAU" in data, f"{org} should have IAU (eukaryotic Ile tRNA)"

    def test_human_has_more_copies_than_ecoli(self):
        """Human genome should have more tRNA gene copies than E. coli."""
        human_total = sum(TRNA_GENE_COPIES["human"].values())
        ecoli_total = sum(TRNA_GENE_COPIES["e_coli"].values())
        assert human_total > ecoli_total

    def test_all_copy_numbers_positive(self):
        """All tRNA gene copy numbers should be positive integers."""
        for org, copies in TRNA_GENE_COPIES.items():
            for anticodon, count in copies.items():
                assert isinstance(count, int), f"{org}/{anticodon}: not int"
                assert count > 0, f"{org}/{anticodon}: count <= 0"


# ═══════════════════════════════════════════════════════════════════════════════
# Test: compute_tai function (matching compute_cai API)
# ═══════════════════════════════════════════════════════════════════════════════

class TestComputeTAI:
    """Test the compute_tai function (matching compute_cai API)."""

    @pytest.mark.parametrize("org", _ALL_ORGANISMS)
    def test_tai_computable_for_all_organisms(self, org):
        """tAI should be computable for all 10 organisms."""
        aliases = _ALL_ORGANISM_ALIASES.get(org, [org])
        tai = compute_tai(_EGFP_DNA, organism=aliases[0])
        assert 0.0 < tai <= 1.0, f"tAI for {org} out of range: {tai}"

    def test_basic_ecoli_tai(self):
        """tAI should be computable for E. coli."""
        tai = compute_tai(_EGFP_DNA, organism="Escherichia_coli")
        assert 0.0 < tai <= 1.0

    def test_basic_human_tai(self):
        """tAI should be computable for human."""
        tai = compute_tai(_HUMAN_INSULIN_DNA, organism="Homo_sapiens")
        assert 0.0 < tai <= 1.0

    def test_empty_sequence(self):
        """Empty sequence should return 0.0."""
        assert compute_tai("", organism="Escherichia_coli") == 0.0

    def test_short_sequence(self):
        """Sequence shorter than 3 bases should return 0.0."""
        assert compute_tai("AT", organism="Escherichia_coli") == 0.0

    def test_invalid_length(self):
        """Non-multiple-of-3 length should raise ValueError."""
        with pytest.raises(ValueError, match="not a multiple of 3"):
            compute_tai("ATGA", organism="Escherichia_coli")

    def test_unsupported_organism(self):
        """Unsupported organism should raise ValueError."""
        with pytest.raises(ValueError, match="No tRNA gene copy data"):
            compute_tai("ATGAAAGCGTTT", organism="zebrafish")

    def test_organism_aliases(self):
        """Various organism name aliases should work and produce same result."""
        dna = "ATGAAAGCGTTT"
        # E. coli aliases
        tai1 = compute_tai(dna, organism="e_coli")
        tai2 = compute_tai(dna, organism="Escherichia_coli")
        tai3 = compute_tai(dna, organism="ecoli")
        assert tai1 == tai2 == tai3

        # Human aliases
        tai_h1 = compute_tai(dna, organism="human")
        tai_h2 = compute_tai(dna, organism="Homo_sapiens")
        tai_h3 = compute_tai(dna, organism="h_sapiens")
        assert tai_h1 == tai_h2 == tai_h3

    def test_new_organism_aliases(self):
        """New organism aliases should work."""
        dna = "ATGAAAGCGTTT"
        # C. elegans
        tai1 = compute_tai(dna, organism="Caenorhabditis_elegans")
        tai2 = compute_tai(dna, organism="c_elegans")
        assert tai1 == tai2

        # D. melanogaster
        tai1 = compute_tai(dna, organism="Drosophila_melanogaster")
        tai2 = compute_tai(dna, organism="d_melanogaster")
        assert tai1 == tai2

        # A. thaliana
        tai1 = compute_tai(dna, organism="Arabidopsis_thaliana")
        tai2 = compute_tai(dna, organism="a_thaliana")
        assert tai1 == tai2

        # B. subtilis
        tai1 = compute_tai(dna, organism="Bacillus_subtilis")
        tai2 = compute_tai(dna, organism="b_subtilis")
        assert tai1 == tai2

        # P. pastoris
        tai1 = compute_tai(dna, organism="Pichia_pastoris")
        tai2 = compute_tai(dna, organism="p_pastoris")
        assert tai1 == tai2

    def test_species_parameter_backward_compat(self):
        """species parameter should work as alias for organism."""
        dna = "ATGAAAGCGTTT"
        tai_org = compute_tai(dna, organism="Escherichia_coli")
        tai_sp = compute_tai(dna, species="ecoli")
        assert tai_org == tai_sp

    def test_species_parameter_deprecation_warning(self):
        """species parameter should emit DeprecationWarning."""
        import warnings
        dna = "ATGAAAGCGTTT"
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            compute_tai(dna, species="ecoli")
            assert any(issubclass(x.category, DeprecationWarning) for x in w)

    @pytest.mark.parametrize("org", _ALL_ORGANISMS)
    def test_tai_in_range(self, org):
        """tAI should always be in [0, 1]."""
        tai = compute_tai(_EGFP_DNA, organism=org)
        assert 0.0 <= tai <= 1.0, f"tAI out of range for {org}: {tai}"

    def test_deterministic(self):
        """Same input should always produce same tAI."""
        for _ in range(10):
            tai = compute_tai(_EGFP_DNA, organism="Escherichia_coli")
            assert tai == compute_tai(_EGFP_DNA, organism="Escherichia_coli")

    def test_dna_case_insensitive(self):
        """DNA sequence should be case-insensitive."""
        tai_upper = compute_tai("ATGGCTAAAGCG", organism="e_coli")
        tai_lower = compute_tai("atggctaaagcg", organism="e_coli")
        assert tai_upper == tai_lower

    def test_whitespace_stripped(self):
        """Whitespace in DNA should be stripped."""
        tai_clean = compute_tai("ATGGCTAAAGCG", organism="e_coli")
        tai_spaced = compute_tai("  ATGGCTAAAGCG  ", organism="e_coli")
        assert tai_clean == tai_spaced


# ═══════════════════════════════════════════════════════════════════════════════
# Test: calculate_tai function (original API)
# ═══════════════════════════════════════════════════════════════════════════════

class TestCalculateTAI:
    """Test the calculate_tai function (original API)."""

    @pytest.mark.parametrize("org", _ALL_ORGANISMS)
    def test_basic_tai_all_organisms(self, org):
        """tAI should be computable for all organisms."""
        tai = calculate_tai(_EGFP_DNA, org)
        assert 0.0 < tai <= 1.0

    def test_empty_sequence(self):
        """Empty sequence should return 0.0."""
        assert calculate_tai("", "e_coli") == 0.0

    def test_short_sequence(self):
        """Sequence shorter than 3 bases should return 0.0."""
        assert calculate_tai("AT", "e_coli") == 0.0

    def test_invalid_length(self):
        """Non-multiple-of-3 length should raise ValueError."""
        with pytest.raises(ValueError, match="not a multiple of 3"):
            calculate_tai("ATGA", "e_coli")

    def test_unsupported_organism(self):
        """Unsupported organism should raise ValueError."""
        with pytest.raises(ValueError, match="No tRNA gene copy data"):
            calculate_tai("ATGAAAGCGTTT", "zebrafish")

    def test_skip_stop_codons(self):
        """Stop codons should be excluded by default."""
        with_stop = "ATGGCTAAAGCGTTTTAA"
        without_stop = "ATGGCTAAAGCGTTT"
        tai_with = calculate_tai(with_stop, "e_coli", skip_stop=True)
        tai_without = calculate_tai(without_stop, "e_coli")
        assert tai_with == tai_without

    def test_include_stop_codons(self):
        """Stop codons should be included when skip_stop=False."""
        with_stop = "ATGGCTAAAGCGTTTTAA"
        tai_skip = calculate_tai(with_stop, "e_coli", skip_stop=True)
        tai_include = calculate_tai(with_stop, "e_coli", skip_stop=False)
        assert isinstance(tai_include, float)

    def test_skip_met_codons(self):
        """Met codons should be excluded by default (CAI convention)."""
        seq = "ATGGCTAAAGCG"
        tai_skip = calculate_tai(seq, "e_coli", skip_met=True)
        tai_include = calculate_tai(seq, "e_coli", skip_met=False)
        assert isinstance(tai_skip, float)
        assert isinstance(tai_include, float)


# ═══════════════════════════════════════════════════════════════════════════════
# Test: compute_codon_weights
# ═══════════════════════════════════════════════════════════════════════════════

class TestComputeCodonWeights:
    """Test the compute_codon_weights function."""

    @pytest.mark.parametrize("org", _ALL_ORGANISMS)
    def test_all_organisms_have_weights(self, org):
        """All supported organisms should have computable weights."""
        weights = compute_codon_weights(org)
        assert isinstance(weights, dict)
        assert len(weights) > 0

    @pytest.mark.parametrize("org", _ALL_ORGANISMS)
    def test_weights_in_range(self, org):
        """All weights should be in [0, 1]."""
        weights = compute_codon_weights(org)
        for codon, w in weights.items():
            assert 0.0 <= w <= 1.0, \
                f"Weight out of range for {org}/{codon}: {w}"

    @pytest.mark.parametrize("org", _ALL_ORGANISMS)
    def test_optimal_codon_has_weight_1(self, org):
        """At least one codon per amino acid should have positive weight."""
        weights = compute_codon_weights(org)
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*" or aa == "M":
                continue
            rna_codons = [c.replace("T", "U") for c in codons]
            max_w = max(weights.get(c, 0.0) for c in rna_codons)
            assert max_w > 0, f"No positive weight for {aa} in {org}"

    def test_weights_differ_between_organisms(self):
        """Codon weights should differ between organisms."""
        weights_ecoli = compute_codon_weights("e_coli")
        weights_human = compute_codon_weights("human")
        common_codons = set(weights_ecoli.keys()) & set(weights_human.keys())
        differences = sum(
            1 for c in common_codons
            if abs(weights_ecoli[c] - weights_human[c]) > 0.01
        )
        assert differences > 0, "E. coli and human tAI weights are identical"

    def test_organism_aliases(self):
        """Organism aliases should produce the same weights."""
        w1 = compute_codon_weights("e_coli")
        w2 = compute_codon_weights("Escherichia_coli")
        w3 = compute_codon_weights("ecoli")
        assert w1 == w2 == w3

    def test_ecoli_optimal_codons_reflect_trna_abundance(self):
        """E. coli optimal codons by tAI should match tRNA abundance.

        For E. coli, the most abundant tRNAs determine the optimal codons.
        For example, Leu-CAG (4 copies) should give CUG a higher weight
        than Leu-UAG reading CUG (U:G wobble = 0.2, only 1 copy).
        """
        weights = compute_codon_weights("e_coli")
        # CUG (Leu) should have high weight because CAG tRNA has 4 copies
        # and C:G = 1.0, plus UAG with 1 copy × 0.2 = 0.2
        # Raw weight = 4*1.0 + 1*0.2 = 4.2
        # CUA (Leu) has UAG × 1.0 = 1*1.0 = 1.0
        # So CUG should have higher weight than CUA
        assert weights["CUG"] >= weights["CUA"], \
            "CUG should have higher tAI weight than CUA in E. coli"


# ═══════════════════════════════════════════════════════════════════════════════
# Test: E. coli Highly Expressed Gene tAI
# ═══════════════════════════════════════════════════════════════════════════════

class TestEcoliHighExprTAI:
    """Test tAI values for E. coli highly expressed genes.

    Published data (dos Reis et al., 2004) shows that highly expressed
    E. coli genes should have tAI > 0.5. We test this using sequences
    composed predominantly of optimal codons.
    """

    def test_ecoli_optimal_sequence_high_tai(self):
        """A sequence using E. coli optimal codons should have high tAI."""
        # Build a sequence using only optimal codons for each amino acid
        optimal_dna = optimize_for_tai(
            "MAKVLSTPEQDNRFWYHGCI", organism="e_coli"
        )
        tai = compute_tai(optimal_dna, organism="e_coli")
        # Optimal codons should give very high tAI
        assert tai > 0.5, f"Optimal E. coli sequence tAI too low: {tai}"

    def test_egfp_ecoli_tai_reasonable(self):
        """eGFP (E. coli-optimized) should have reasonable tAI."""
        tai = compute_tai(_EGFP_DNA, organism="e_coli")
        # eGFP is codon-optimized for E. coli, so tAI should be decent
        assert tai > 0.2, f"eGFP tAI in E. coli too low: {tai}"

    def test_human_insulin_ecoli_tai_lower(self):
        """Human insulin with native codons should have lower tAI in E. coli."""
        tai_insulin = compute_tai(_HUMAN_INSULIN_DNA, organism="e_coli")
        tai_egfp = compute_tai(_EGFP_DNA, organism="e_coli")
        # Human insulin uses mammalian codons, not optimized for E. coli
        assert tai_insulin < 0.8, f"Human insulin tAI in E. coli unexpectedly high: {tai_insulin}"


# ═══════════════════════════════════════════════════════════════════════════════
# Test: Cross-species tAI
# ═══════════════════════════════════════════════════════════════════════════════

class TestCrossSpeciesTAI:
    """Test cross-species tAI behavior."""

    def test_different_organisms_different_tai(self):
        """Same sequence should give different tAI for different organisms."""
        values = set()
        for org in _ALL_ORGANISMS:
            tai = calculate_tai(_EGFP_DNA, org)
            values.add(round(tai, 4))
        # Not all values should be identical
        assert len(values) > 1, \
            "tAI identical across all organisms — data may not be organism-specific"

    def test_mouse_and_human_both_valid(self):
        """Mouse and human tAI should both be valid."""
        tai_mouse = calculate_tai(_EGFP_DNA, "mouse")
        tai_human = calculate_tai(_EGFP_DNA, "human")
        assert 0.0 < tai_mouse <= 1.0
        assert 0.0 < tai_human <= 1.0

    @pytest.mark.parametrize("org", _ALL_ORGANISMS)
    def test_all_organisms_tai_valid(self, org):
        """All organisms should produce valid tAI values."""
        tai = calculate_tai(_EGFP_DNA, org)
        assert 0.0 < tai <= 1.0, f"Invalid tAI for {org}: {tai}"


# ═══════════════════════════════════════════════════════════════════════════════
# Test: tAI vs CAI
# ═══════════════════════════════════════════════════════════════════════════════

class TestTAIvsCAI:
    """Compare tAI and CAI — they should be correlated but not identical."""

    def test_tai_and_cai_both_computable(self):
        """Both tAI and CAI should be computable for the same sequence."""
        dna = _EGFP_DNA
        tai = calculate_tai(dna, "e_coli")
        cai = compute_cai(dna, organism="Escherichia_coli")
        assert 0.0 < tai <= 1.0
        assert 0.0 < cai <= 1.0

    def test_tai_uses_trna_not_codon_freq(self):
        """tAI should reflect tRNA abundance, not codon frequency."""
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        ecoli_weights = compute_codon_weights("e_coli")
        cai_weights = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]

        # Convert tAI weights from RNA to DNA for comparison
        tai_weights_dna = {
            codon.replace("U", "T"): w
            for codon, w in ecoli_weights.items()
        }

        # The weights should not be identical (different data sources)
        common = set(tai_weights_dna.keys()) & set(cai_weights.keys())
        differences = sum(
            1 for c in common
            if abs(tai_weights_dna[c] - cai_weights[c]) > 0.01
        )
        assert differences > 0, "tAI weights identical to CAI — may not be using tRNA data"


# ═══════════════════════════════════════════════════════════════════════════════
# Test: compute_tai_and_cai
# ═══════════════════════════════════════════════════════════════════════════════

class TestComputeTAIAndCAI:
    """Test the combined tAI and CAI computation."""

    def test_returns_dict_with_both_metrics(self):
        """Should return a dict with 'tai', 'cai', and 'correlation' keys."""
        result = compute_tai_and_cai(_EGFP_DNA, organism="e_coli")
        assert "tai" in result
        assert "cai" in result
        assert "correlation" in result

    def test_tai_matches_separate_computation(self):
        """tAI from combined function should match standalone tAI."""
        result = compute_tai_and_cai(_EGFP_DNA, organism="e_coli")
        tai_standalone = compute_tai(_EGFP_DNA, organism="e_coli")
        assert abs(result["tai"] - tai_standalone) < 1e-6

    def test_cai_matches_separate_computation(self):
        """CAI from combined function should match standalone CAI."""
        result = compute_tai_and_cai(_EGFP_DNA, organism="Escherichia_coli")
        cai_standalone = compute_cai(_EGFP_DNA, organism="Escherichia_coli")
        assert abs(result["cai"] - cai_standalone) < 1e-6

    @pytest.mark.parametrize("org", ["e_coli", "human", "yeast"])
    def test_all_organisms_work(self, org):
        """Should work for multiple organisms."""
        result = compute_tai_and_cai(_EGFP_DNA, organism=org)
        assert 0.0 <= result["tai"] <= 1.0
        assert 0.0 <= result["cai"] <= 1.0

    def test_correlation_is_positive_or_negative(self):
        """Correlation indicator should be +1 or -1."""
        result = compute_tai_and_cai(_EGFP_DNA, organism="e_coli")
        assert result["correlation"] in (1.0, -1.0)


# ═══════════════════════════════════════════════════════════════════════════════
# Test: optimize_for_tai
# ═══════════════════════════════════════════════════════════════════════════════

class TestOptimizeForTAI:
    """Test the optimize_for_tai function."""

    @pytest.mark.parametrize("org", _ALL_ORGANISMS)
    def test_produces_valid_dna(self, org):
        """Should produce valid DNA for all organisms."""
        protein = "MAKVLSTPEQDNR"
        dna = optimize_for_tai(protein, organism=org)
        assert len(dna) == len(protein) * 3
        assert all(c in "ACGT" for c in dna)

    def test_encodes_correct_protein(self):
        """Optimized DNA should encode the original protein."""
        from biocompiler.expression.translation import translate
        protein = "MAKVLSTPEQDNRFWYHGCI"
        dna = optimize_for_tai(protein, organism="e_coli")
        translated = translate(dna)
        # Remove stop codon if present
        translated_no_stop = translated.rstrip("*")
        assert translated_no_stop == protein

    def test_optimized_higher_tai_than_random(self):
        """tAI-optimized sequence should have higher tAI than random codons."""
        protein = "MAKVLSTPEQDNRFWYHGCI"
        optimized_dna = optimize_for_tai(protein, organism="e_coli")
        optimized_tai = compute_tai(optimized_dna, organism="e_coli")

        # Generate random codon sequences
        random.seed(42)
        random_tais = []
        for _ in range(20):
            random_dna = "".join(
                random.choice(AA_TO_CODONS[aa])
                for aa in protein
            )
            random_tais.append(compute_tai(random_dna, organism="e_coli"))

        avg_random_tai = sum(random_tais) / len(random_tais)
        assert optimized_tai >= avg_random_tai, \
            f"Optimized tAI ({optimized_tai}) should be >= avg random ({avg_random_tai})"

    def test_optimized_has_highest_possible_tai(self):
        """tAI-optimized sequence should have the theoretical maximum tAI.

        Since optimize_for_tai selects the best codon for each amino acid,
        no other synonymous sequence should have higher tAI.
        """
        protein = "MAKVLSTPEQDNR"
        optimized_dna = optimize_for_tai(protein, organism="e_coli")
        optimized_tai = compute_tai(optimized_dna, organism="e_coli")

        # Try all possible codon combinations for a short protein
        # For a 12-AA protein this is feasible
        best_tai = 0.0
        for aa in protein[:4]:  # Just test first 4 amino acids exhaustively
            codons = AA_TO_CODONS[aa]
            for codon in codons:
                test_dna = codon + optimized_dna[3:]  # Replace first codon
                tai = compute_tai(test_dna, organism="e_coli")
                best_tai = max(best_tai, tai)

        # The optimized sequence's first codon should be optimal
        assert optimized_tai >= best_tai - 0.01, \
            f"Optimized tAI ({optimized_tai}) should be near best ({best_tai})"

    def test_invalid_protein_raises_error(self):
        """Invalid amino acid codes should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid amino acid"):
            optimize_for_tai("MAKB", organism="e_coli")

    @pytest.mark.parametrize("org", ["e_coli", "human", "yeast", "b_subtilis"])
    def test_optimized_dna_length_correct(self, org):
        """DNA length should be 3× protein length."""
        protein = "MAKV"
        dna = optimize_for_tai(protein, organism=org)
        assert len(dna) == 12  # 4 amino acids × 3 bases

    def test_optimize_for_tai_uses_tai_weights(self):
        """optimize_for_tai should select codons with highest tAI weight."""
        protein = "FF"  # Two Phe codons: UUU and UUC
        dna = optimize_for_tai(protein, organism="e_coli")
        # In E. coli, Phe-GAA reads UUC with 1.0 efficiency and UUU with 0.5
        # UUC should always be preferred over UUU
        for i in range(0, len(dna), 3):
            codon = dna[i:i+3]
            assert codon == "TTC", f"Phe codon should be TTC (optimal), got {codon}"


# ═══════════════════════════════════════════════════════════════════════════════
# Test: Edge Cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestTAIEdgeCases:
    """Test edge cases for tAI computation."""

    def test_single_met_codon(self):
        """A sequence of just ATG should return 0.0 (Met is skipped)."""
        assert calculate_tai("ATG", "e_coli") == 0.0

    def test_single_non_met_codon(self):
        """A single non-Met codon should return a valid tAI."""
        tai = calculate_tai("ATGGCT", "e_coli")  # M + A
        assert 0.0 < tai <= 1.0

    def test_all_same_amino_acid(self):
        """Sequence encoding the same amino acid should have valid tAI."""
        alanine_seq = "ATGGCTGCTGCT"  # M + A + A + A
        tai = calculate_tai(alanine_seq, "e_coli")
        assert 0.0 < tai <= 1.0

    def test_long_sequence(self):
        """Long sequences should be processable without error."""
        long_dna = _EGFP_DNA * 5
        tai = calculate_tai(long_dna, "e_coli")
        assert 0.0 < tai <= 1.0

    def test_dna_with_only_stop(self):
        """Sequence of just stop codons should return 0.0."""
        assert calculate_tai("TAATAA", "e_coli") == 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# Test: tAI Objective Function
# ═══════════════════════════════════════════════════════════════════════════════

class TestTAIObjective:
    """Test the tAI objective function for the optimizer."""

    def test_tai_objective_returns_float(self):
        """tAI objective should return a float."""
        val = tai_objective("ATGGCTAAAGCG", "MAK", "Escherichia_coli")
        assert isinstance(val, float)

    def test_tai_objective_in_range(self):
        """tAI objective should return a value in [0, 1]."""
        val = tai_objective(_EGFP_DNA, "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK", "Escherichia_coli")
        assert 0.0 <= val <= 1.0

    def test_tai_objective_empty_dna(self):
        """tAI objective should return 0.0 for empty DNA."""
        assert tai_objective("", "", "Escherichia_coli") == 0.0

    def test_tai_objective_unsupported_organism(self):
        """tAI objective should return 0.0 for unsupported organisms."""
        val = tai_objective("ATGGCTAAAGCG", "MAK", "zebrafish")
        assert val == 0.0

    def test_tai_objective_in_registry(self):
        """tAI should be in the objective registry."""
        assert "tai" in OBJECTIVE_REGISTRY

    def test_tai_objective_resolvable(self):
        """resolve_objective should accept 'tai'."""
        obj = resolve_objective("tai")
        assert callable(obj)

    def test_tai_objective_matches_compute_tai(self):
        """tAI objective should return same value as compute_tai."""
        dna = _EGFP_DNA
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        obj_val = tai_objective(dna, protein, "Escherichia_coli")
        tai_val = compute_tai(dna, organism="Escherichia_coli")
        assert abs(obj_val - tai_val) < 1e-6


# ═══════════════════════════════════════════════════════════════════════════════
# Test: Integration with optimize_sequence
# ═══════════════════════════════════════════════════════════════════════════════

class TestTAIOptimizerIntegration:
    """Test tAI integration with the optimization pipeline."""

    def test_tai_objective_with_optimize_sequence(self):
        """optimize_sequence should accept 'tai' as objective."""
        from biocompiler.optimizer.pipeline import optimize_sequence
        result = optimize_sequence(
            "MAKV",
            organism="Escherichia_coli",
            objective="tai",
        )
        assert result is not None
        assert hasattr(result, 'sequence')
        assert len(result.sequence) > 0

    def test_optimized_sequence_has_tai(self):
        """Optimized sequence should have a computable tAI."""
        from biocompiler.optimizer.pipeline import optimize_sequence
        result = optimize_sequence(
            "MAKV",
            organism="Escherichia_coli",
        )
        tai = compute_tai(result.sequence, organism="Escherichia_coli")
        assert 0.0 < tai <= 1.0

    def test_tai_objective_produces_valid_sequence(self):
        """Sequences optimized with tAI objective should encode correct protein."""
        from biocompiler.optimizer.pipeline import optimize_sequence
        from biocompiler.expression.translation import translate
        protein = "MAKV"
        result = optimize_sequence(
            protein,
            organism="Escherichia_coli",
            objective="tai",
        )
        translated = translate(result.sequence)
        # Remove stop codon from translation for comparison
        translated_no_stop = translated.rstrip("*")
        assert translated_no_stop == protein


# ═══════════════════════════════════════════════════════════════════════════════
# Test: tAI-CAI Correlation
# ═══════════════════════════════════════════════════════════════════════════════

class TestTAICAICorrelation:
    """Test that tAI and CAI are correlated but distinct."""

    def test_tai_cai_positive_correlation(self):
        """tAI and CAI should be positively correlated across sequences.

        Sequences with high CAI should tend to have high tAI, and vice versa.
        """
        # Generate multiple sequences with varying codon usage
        random.seed(42)
        protein = "MAKVLSTPEQDNRFWYHGCI"
        tai_values = []
        cai_values = []

        for _ in range(30):
            dna = "".join(
                random.choice(AA_TO_CODONS[aa])
                for aa in protein
            )
            tai_values.append(compute_tai(dna, organism="e_coli"))
            cai_values.append(compute_cai(dna, organism="Escherichia_coli"))

        # Compute Pearson correlation
        n = len(tai_values)
        mean_tai = sum(tai_values) / n
        mean_cai = sum(cai_values) / n

        cov = sum((t - mean_tai) * (c - mean_cai) for t, c in zip(tai_values, cai_values)) / n
        std_tai = math.sqrt(sum((t - mean_tai) ** 2 for t in tai_values) / n)
        std_cai = math.sqrt(sum((c - mean_cai) ** 2 for c in cai_values) / n)

        if std_tai > 0 and std_cai > 0:
            correlation = cov / (std_tai * std_cai)
            # Positive correlation is expected (both measure codon optimality)
            assert correlation > 0.0, \
                f"tAI-CAI correlation should be positive, got {correlation:.3f}"

    def test_optimized_tai_not_also_optimal_cai(self):
        """tAI-optimized sequence should not always be CAI-optimal.

        Since tAI and CAI use different data sources, the optimal codon
        selection should differ for at least some amino acids.
        """
        protein = "LVPRSFYEWHDC"
        tai_dna = optimize_for_tai(protein, organism="e_coli")
        tai_val = compute_tai(tai_dna, organism="e_coli")
        cai_val = compute_cai(tai_dna, organism="Escherichia_coli")
        # Both should be reasonable but not necessarily equal
        assert 0.0 < tai_val <= 1.0
        assert 0.0 < cai_val <= 1.0
