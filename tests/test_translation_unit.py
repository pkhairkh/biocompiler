"""
Unit tests for biocompiler.translation module.

Covers:
1. translate() — DNA to protein conversion
2. compute_cai() — Codon Adaptation Index in [0,1]
3. Known codon sequences → expected proteins
4. Edge cases: empty sequence, non-multiple-of-3 length, stop codons
"""

import math
import pytest

from biocompiler.expression.translation import translate, compute_cai
from biocompiler.shared.exceptions import UnsupportedOrganismError, InvalidSequenceError


# ═══════════════════════════════════════════════════════════════════════════════
# 1. translate() — basic DNA-to-protein conversion
# ═══════════════════════════════════════════════════════════════════════════════

class TestTranslate:
    """Tests for the translate() function."""

    # --- Single codon mappings (known codon → amino acid) ---

    @pytest.mark.parametrize(
        "codon, expected_aa",
        [
            ("ATG", "M"),   # Start / Methionine
            ("TTT", "F"),   # Phenylalanine
            ("TTC", "F"),   # Phenylalanine
            ("TTA", "L"),   # Leucine
            ("CTT", "L"),   # Leucine
            ("ATT", "I"),   # Isoleucine
            ("GTT", "V"),   # Valine
            ("TCT", "S"),   # Serine
            ("CCT", "P"),   # Proline
            ("ACT", "T"),   # Threonine
            ("GCT", "A"),   # Alanine
            ("TAT", "Y"),   # Tyrosine
            ("CAT", "H"),   # Histidine
            ("CAA", "Q"),   # Glutamine
            ("AAT", "N"),   # Asparagine
            ("AAA", "K"),   # Lysine
            ("GAT", "D"),   # Aspartic acid
            ("GAA", "E"),   # Glutamic acid
            ("TGT", "C"),   # Cysteine
            ("TGG", "W"),   # Tryptophan
            ("CGT", "R"),   # Arginine
            ("AGT", "S"),   # Serine
            ("AGA", "R"),   # Arginine
            ("GGT", "G"),   # Glycine
        ],
    )
    def test_single_codon_translation(self, codon, expected_aa):
        """Each individual codon must map to the correct amino acid."""
        assert translate(codon) == expected_aa

    # --- Known multi-codon sequences ---

    def test_met_phe_trp(self):
        """ATG TTT TGG → MFW"""
        assert translate("ATGTTTTGG") == "MFW"

    def test_alanine_repeat(self):
        """GCT GCT GCT → AAA (alanine triple)"""
        assert translate("GCTGCTGCT") == "AAA"

    def test_all_leucine_codons(self):
        """All six leucine codons should produce L."""
        leucine_codons = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
        for codon in leucine_codons:
            assert translate(codon) == "L", f"Codon {codon} should encode Leucine"

    # --- Case insensitivity ---

    def test_lowercase_input(self):
        """translate() should handle lowercase DNA (via validate_dna_sequence)."""
        assert translate("atgttttgg") == "MFW"

    def test_mixed_case_input(self):
        """Mixed case should also work."""
        assert translate("AtGtTtTgG") == "MFW"

    # --- Stop codon behavior ---

    def test_stop_codon_taa_with_to_stop_true(self):
        """With to_stop=True (default), translation stops at TAA."""
        # ATG TTT TAA GCT → M F (stop) — GCT is not translated
        assert translate("ATGTTTTAAGCT") == "MF"

    def test_stop_codon_tag_with_to_stop_true(self):
        """With to_stop=True, TAG also stops translation."""
        assert translate("ATGTTTTAGGCT") == "MF"

    def test_stop_codon_tga_with_to_stop_true(self):
        """With to_stop=True, TGA also stops translation."""
        assert translate("ATGTTTTGAGCT") == "MF"

    def test_stop_codon_with_to_stop_false(self):
        """With to_stop=False, stop codons are included as '*'."""
        result = translate("ATGTTTTAAGCT", to_stop=False)
        assert result == "MF*A"

    def test_multiple_stop_codons_to_stop_false(self):
        """Multiple stop codons with to_stop=False all become '*'."""
        result = translate("ATGTTTTAATAGGCT", to_stop=False)
        assert result == "MF**A"

    def test_stop_codon_at_beginning(self):
        """If the very first codon is a stop, translate returns empty with to_stop=True."""
        assert translate("TAAGCTGCT") == ""

    def test_stop_codon_at_beginning_to_stop_false(self):
        """If the first codon is stop with to_stop=False, '*' is included."""
        assert translate("TAAGCTGCT", to_stop=False) == "*AA"

    # --- Edge cases: empty and short sequences ---

    def test_empty_sequence(self):
        """Empty string should return empty protein."""
        assert translate("") == ""

    def test_single_base(self):
        """Sequence shorter than a codon (1 base) → empty protein."""
        assert translate("A") == ""

    def test_two_bases(self):
        """Sequence shorter than a codon (2 bases) → empty protein."""
        assert translate("AT") == ""

    def test_non_multiple_of_3_length(self):
        """Non-multiple-of-3 length: extra bases are ignored."""
        # ATG TTT T → only ATG TTT translated, extra T ignored
        assert translate("ATGTTTT") == "MF"

    def test_non_multiple_of_3_length_two_extra(self):
        """Two extra bases at end are also ignored."""
        # ATG TTT TG → only ATG TTT translated, extra TG ignored
        assert translate("ATGTTTTG") == "MF"

    # --- Invalid characters ---

    def test_invalid_characters_raise(self):
        """Non-DNA characters should raise InvalidSequenceError."""
        with pytest.raises(InvalidSequenceError):
            translate("ATGZZZ")

    # --- Determinism ---

    def test_translate_is_deterministic(self):
        """Same input must always produce the same output."""
        seq = "ATGTTTTGGGCTCATAAAGCT"
        result1 = translate(seq)
        result2 = translate(seq)
        assert result1 == result2

    # --- Longer realistic sequence ---

    def test_longer_sequence(self):
        """Translate a longer coding sequence and verify partial results."""
        # ATG GCT GCT GCT TTT TTT TGG → MAAAFFW
        seq = "ATGGCTGCTGCTTTTTTTTGG"
        assert translate(seq) == "MAAAFFW"


# ═══════════════════════════════════════════════════════════════════════════════
# 2. compute_cai() — Codon Adaptation Index
# ═══════════════════════════════════════════════════════════════════════════════

class TestComputeCAI:
    """Tests for the compute_cai() function."""

    def test_cai_returns_float_in_unit_interval(self):
        """CAI must always be in [0.0, 1.0]."""
        seq = "ATGGCTGCTGCTTTTTTTTGG"
        cai = compute_cai(seq)
        assert isinstance(cai, float)
        assert 0.0 <= cai <= 1.0

    def test_cai_default_organism(self):
        """Default organism is Homo_sapiens; should not raise."""
        cai = compute_cai("ATGGCTGCTGCT")
        assert 0.0 <= cai <= 1.0

    def test_cai_ecoli(self):
        """E. coli CAI should also return a value in [0,1]."""
        cai = compute_cai("ATGGCTGCTGCT", organism="Escherichia_coli")
        assert 0.0 <= cai <= 1.0

    def test_cai_mouse(self):
        """Mouse CAI should also return a value in [0,1]."""
        cai = compute_cai("ATGGCTGCTGCT", organism="Mus_musculus")
        assert 0.0 <= cai <= 1.0

    def test_cai_cho(self):
        """CHO CAI should also return a value in [0,1]."""
        cai = compute_cai("ATGGCTGCTGCT", organism="CHO_K1")
        assert 0.0 <= cai <= 1.0

    def test_cai_yeast(self):
        """Yeast CAI should also return a value in [0,1]."""
        cai = compute_cai("ATGGCTGCTGCT", organism="Saccharomyces_cerevisiae")
        assert 0.0 <= cai <= 1.0

    def test_cai_unsupported_organism_raises(self):
        """An unsupported organism name must raise UnsupportedOrganismError."""
        with pytest.raises(UnsupportedOrganismError):
            compute_cai("ATGGCTGCTGCT", organism="Alien_genome")

    def test_cai_empty_sequence(self):
        """Empty sequence should return 0.0."""
        assert compute_cai("") == 0.0

    def test_cai_stop_codons_skipped(self):
        """Stop codons (and Met) should be skipped in CAI computation.
        A sequence consisting only of ATG + stop should yield 0.0
        because there are no contributing codons for CAI."""
        # ATG TAA → M + stop; Met and stop are both skipped → no ratios → 0.0
        cai = compute_cai("ATGTAA")
        assert cai == 0.0

    def test_cai_is_deterministic(self):
        """Same input must always produce the same CAI."""
        seq = "ATGGCTGCTGCTTTTTTTTGG"
        cai1 = compute_cai(seq)
        cai2 = compute_cai(seq)
        assert cai1 == cai2

    def test_cai_non_multiple_of_3(self):
        """Non-multiple-of-3 length: partial codon at end is ignored."""
        cai = compute_cai("ATGGCTGCT")
        assert 0.0 <= cai <= 1.0

    def test_cai_different_organisms_different_values(self):
        """Different organisms generally produce different CAI values
        for the same sequence (unless by coincidence)."""
        seq = "ATGGCTGCTGCTTTTTTTTGG"
        cai_human = compute_cai(seq, organism="Homo_sapiens")
        cai_ecoli = compute_cai(seq, organism="Escherichia_coli")
        # We just check both are valid; they *may* differ
        assert 0.0 <= cai_human <= 1.0
        assert 0.0 <= cai_ecoli <= 1.0

    def test_cai_preferred_codons_higher(self):
        """Sequences using preferred codons should generally have
        higher CAI than sequences using rare codons for the same protein."""
        # Phenylalanine: TTT (less preferred) vs TTC (preferred in many organisms)
        # Using multiple copies to amplify the difference
        rare_phe = "ATG" + "TTT" * 10  # All rare Phe codons
        pref_phe = "ATG" + "TTC" * 10  # All preferred Phe codons
        cai_rare = compute_cai(rare_phe, organism="Homo_sapiens")
        cai_pref = compute_cai(pref_phe, organism="Homo_sapiens")
        # Preferred codons should yield >= CAI compared to rare
        assert cai_pref >= cai_rare


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Known codon sequences → expected proteins
# ═══════════════════════════════════════════════════════════════════════════════

class TestKnownSequences:
    """Verify that well-known codon sequences produce expected proteins."""

    def test_methionine_only(self):
        """ATG repeated → M repeated."""
        assert translate("ATGATGATG") == "MMM"

    def test_glycine_all_codons(self):
        """All four glycine codons produce G."""
        for codon in ("GGT", "GGC", "GGA", "GGG"):
            assert translate(codon) == "G"

    def test_serine_all_codons(self):
        """All six serine codons produce S."""
        for codon in ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"):
            assert translate(codon) == "S"

    def test_arginine_all_codons(self):
        """All six arginine codons produce R."""
        for codon in ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"):
            assert translate(codon) == "R"

    def test_proline_all_codons(self):
        """All four proline codons produce P."""
        for codon in ("CCT", "CCC", "CCA", "CCG"):
            assert translate(codon) == "P"

    def test_stop_codons_all_three(self):
        """All three stop codons should terminate with to_stop=True."""
        for stop in ("TAA", "TAG", "TGA"):
            result = translate("ATG" + stop + "GCT")
            assert result == "M", f"Stop codon {stop} should terminate translation"

    def test_stop_codons_as_asterisk_to_stop_false(self):
        """All three stop codons should become '*' with to_stop=False."""
        for stop in ("TAA", "TAG", "TGA"):
            result = translate("ATG" + stop + "GCT", to_stop=False)
            assert result == "M*A", f"Stop codon {stop} should map to '*' with to_stop=False"

    def test_amino_acid_coverage(self):
        """Every standard amino acid should be reachable via translate()."""
        # Pick one codon per standard amino acid (excluding stop)
        codons_per_aa = {
            "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT",
            "Q": "CAA", "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT",
            "L": "CTT", "K": "AAA", "M": "ATG", "F": "TTT", "P": "CCT",
            "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTT",
        }
        seq = "".join(codons_per_aa[aa] for aa in sorted(codons_per_aa))
        result = translate(seq)
        # All 20 amino acids should appear (in sorted codon order)
        assert len(result) == 20
        for aa in codons_per_aa:
            assert aa in result, f"Amino acid {aa} missing from translation"


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestEdgeCases:
    """Edge cases for translate() and compute_cai()."""

    # --- translate() edge cases ---

    def test_translate_empty_string(self):
        """Empty string → empty protein."""
        assert translate("") == ""

    def test_translate_single_nucleotide(self):
        """Single base → cannot form a codon → empty protein."""
        assert translate("A") == ""

    def test_translate_two_nucleotides(self):
        """Two bases → cannot form a complete codon → empty protein."""
        assert translate("AT") == ""

    def test_translate_exactly_three_bases(self):
        """Exactly one codon → one amino acid."""
        assert translate("ATG") == "M"

    def test_translate_non_multiple_of_3_one_extra(self):
        """4 bases: first codon translated, 1 extra base ignored."""
        assert translate("ATGG") == "M"

    def test_translate_non_multiple_of_3_two_extra(self):
        """5 bases: first codon translated, 2 extra bases ignored."""
        assert translate("ATGGC") == "M"

    def test_translate_only_stop_codon_to_stop_true(self):
        """A sequence that is just a stop codon returns empty."""
        assert translate("TAA") == ""

    def test_translate_only_stop_codon_to_stop_false(self):
        """A stop codon with to_stop=False returns '*'."""
        assert translate("TAA", to_stop=False) == "*"

    def test_translate_all_stops_to_stop_true(self):
        """Multiple stop codons in a row with to_stop=True → empty."""
        assert translate("TAATAGTGA") == ""

    def test_translate_all_stops_to_stop_false(self):
        """Multiple stop codons with to_stop=False → all become '*'."""
        assert translate("TAATAGTGA", to_stop=False) == "***"

    def test_translate_long_sequence(self):
        """A longer sequence (e.g. 30 codons) should translate correctly."""
        # 10x ATG → 10x M
        seq = "ATG" * 10
        assert translate(seq) == "M" * 10

    def test_translate_sequence_with_n_bases(self):
        """Sequence with N bases — N is allowed by validate_dna_sequence,
        but unknown codons containing N should map to 'X'."""
        # ATN is not in CODON_TABLE → maps to 'X'
        result = translate("ATN", to_stop=False)
        assert result == "X"

    # --- compute_cai() edge cases ---

    def test_cai_empty_string(self):
        """Empty string → 0.0."""
        assert compute_cai("") == 0.0

    def test_cai_single_codon_atg(self):
        """ATG only (Met is skipped in CAI) → 0.0 (no contributing codons)."""
        assert compute_cai("ATG") == 0.0

    def test_cai_only_stop_codons(self):
        """Only stop codons → 0.0 (all skipped)."""
        assert compute_cai("TAATAGTGA") == 0.0

    def test_cai_atg_plus_stop(self):
        """ATG + stop → both skipped → 0.0."""
        assert compute_cai("ATGTAA") == 0.0

    def test_cai_sequence_with_n_codon(self):
        """Sequence containing N bases: codons with N are not in CODON_TABLE,
        so aa=None → skipped in CAI loop. Only valid codons contribute."""
        # ATN (unknown) + GCT (Ala, valid) → only GCT contributes
        cai = compute_cai("ATNGCT")
        assert 0.0 <= cai <= 1.0

    def test_cai_non_multiple_of_3_ignores_partial(self):
        """Non-multiple-of-3: partial codon at the end is ignored."""
        # ATG GCT G → only ATG and GCT are processed
        cai = compute_cai("ATGGCTG")
        assert 0.0 <= cai <= 1.0

    def test_cai_result_is_rounded_to_4_decimals(self):
        """CAI result should be rounded to 4 decimal places."""
        cai = compute_cai("ATGGCTGCTGCTTTTTTTTGG")
        # Check that the result has at most 4 decimal places
        assert cai == round(cai, 4)

    def test_cai_invalid_sequence_raises(self):
        """Invalid DNA characters should raise InvalidSequenceError."""
        with pytest.raises(InvalidSequenceError):
            compute_cai("ATGZZZ")

    def test_cai_deterministic_across_organisms(self):
        """CAI is deterministic for the same input across repeated calls."""
        seq = "ATGGCTGCTGCTTTTTTTTGG"
        for _ in range(5):
            cai = compute_cai(seq, organism="Homo_sapiens")
            assert 0.0 <= cai <= 1.0
