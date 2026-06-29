"""Unit tests for biocompiler.constants.

Covers:
1. CODON_TABLE has all 64 codons
2. AA_TO_CODONS maps all 20 standard amino acids
3. BLOSUM62 is a valid 20x20 matrix
4. HYDROPATHY has entries for all 20 AAs
5. RESTRICTION_ENZYMES has known enzymes
6. reverse_complement() works correctly
"""

from __future__ import annotations

import pytest

from biocompiler.shared.constants import (
    AA_TO_CODONS,
    BLOSUM62,
    CODON_TABLE,
    COMPLEMENT,
    HYDROPATHY,
    HYDROPHOBIC_AAS,
    RESTRICTION_ENZYMES,
    STANDARD_AAS,
    STANDARD_AAS_BLOSUM_ORDER,
    START_CODON,
    STOP_CODONS,
    reverse_complement,
)


# ═══════════════════════════════════════════════════════════════════════════
# 1. CODON_TABLE
# ═══════════════════════════════════════════════════════════════════════════

class TestCodonTable:
    """Verify CODON_TABLE is a complete, valid genetic code."""

    def test_has_64_codons(self):
        """The standard genetic code has exactly 64 codons."""
        assert len(CODON_TABLE) == 64

    def test_all_codons_are_three_bases(self):
        """Every key in CODON_TABLE must be a 3-character DNA string."""
        for codon in CODON_TABLE:
            assert len(codon) == 3, f"Codon {codon!r} is not 3 bases"
            assert all(b in "ACGT" for b in codon), f"Codon {codon!r} contains non-ACGT"

    def test_all_values_are_single_aa_or_stop(self):
        """Every value must be a single uppercase letter or '*' for stop."""
        for codon, aa in CODON_TABLE.items():
            assert len(aa) == 1, f"Codon {codon} maps to {aa!r} (not single char)"
            assert aa.isupper() or aa == "*", f"Codon {codon} maps to {aa!r} (not uppercase or stop)"

    def test_no_duplicate_codons(self):
        """All codon keys must be unique (guaranteed by dict, but verify)."""
        assert len(set(CODON_TABLE.keys())) == 64

    def test_covers_all_nucleotide_combinations(self):
        """All 4^3 = 64 possible triplet combinations are present."""
        import itertools
        expected = {"".join(combo) for combo in itertools.product("ACGT", repeat=3)}
        actual = set(CODON_TABLE.keys())
        assert actual == expected

    def test_stop_codons_match(self):
        """STOP_CODONS must be exactly the three codons mapping to '*'."""
        stops_in_table = {c for c, aa in CODON_TABLE.items() if aa == "*"}
        assert STOP_CODONS == stops_in_table
        assert len(STOP_CODONS) == 3

    def test_start_codon_is_atg(self):
        """START_CODON must be ATG (Methionine)."""
        assert START_CODON == "ATG"
        assert CODON_TABLE[START_CODON] == "M"

    def test_specific_codon_mappings(self):
        """Spot-check a few well-known codon assignments."""
        checks = {
            "ATG": "M",   # Methionine (start)
            "TGG": "W",   # Tryptophan (single codon)
            "TAA": "*",   # Stop (ochre)
            "TAG": "*",   # Stop (amber)
            "TGA": "*",   # Stop (opal)
            "TTT": "F",   # Phenylalanine
            "GGG": "G",   # Glycine
        }
        for codon, expected_aa in checks.items():
            assert CODON_TABLE[codon] == expected_aa, f"{codon} should map to {expected_aa}"


# ═══════════════════════════════════════════════════════════════════════════
# 2. AA_TO_CODONS
# ═══════════════════════════════════════════════════════════════════════════

class TestAAToCodons:
    """Verify AA_TO_CODONS covers all 20 standard amino acids."""

    def test_has_20_amino_acids(self):
        """AA_TO_CODONS must map exactly 20 standard amino acids."""
        assert len(AA_TO_CODONS) == 20

    def test_all_standard_aas_present(self):
        """All 20 standard amino acids must be keys."""
        standard = set(STANDARD_AAS)
        actual = set(AA_TO_CODONS.keys())
        assert actual == standard, f"Missing: {standard - actual}, Extra: {actual - standard}"

    def test_no_stop_codon_in_values(self):
        """Stop codons ('*') should not appear in AA_TO_CODONS."""
        assert "*" not in AA_TO_CODONS

    def test_all_codons_valid(self):
        """Every codon listed must be a valid 3-base DNA triplet."""
        for aa, codons in AA_TO_CODONS.items():
            for codon in codons:
                assert len(codon) == 3, f"Codon {codon!r} for {aa} is not 3 bases"
                assert all(b in "ACGT" for b in codon), f"Codon {codon!r} contains non-ACGT"

    def test_codons_consistent_with_codon_table(self):
        """AA_TO_CODONS must be the exact inverse of CODON_TABLE (excluding stops)."""
        # Forward: every codon listed under an AA must map to that AA in CODON_TABLE
        for aa, codons in AA_TO_CODONS.items():
            for codon in codons:
                assert CODON_TABLE[codon] == aa, f"CODON_TABLE[{codon}]={CODON_TABLE[codon]}, expected {aa}"

        # Reverse: every non-stop codon in CODON_TABLE must be listed under its AA
        for codon, aa in CODON_TABLE.items():
            if aa != "*":
                assert codon in AA_TO_CODONS[aa], f"Codon {codon} not listed under {aa}"

    def test_total_codons_sum_to_61(self):
        """61 sense codons (64 total minus 3 stop) distributed across 20 AAs."""
        total = sum(len(codons) for codons in AA_TO_CODONS.values())
        assert total == 61

    def test_known_codon_counts(self):
        """Spot-check degeneracy of specific amino acids."""
        known = {
            "M": 1,  # AUG only
            "W": 1,  # UGG only
            "L": 6,  # 6-fold degenerate
            "S": 6,  # 6-fold degenerate
            "R": 6,  # 6-fold degenerate
            "F": 2,  # 2-fold
            "K": 2,  # 2-fold
            "I": 3,  # 3-fold
        }
        for aa, expected_count in known.items():
            assert len(AA_TO_CODONS[aa]) == expected_count, (
                f"{aa} should have {expected_count} codons, got {len(AA_TO_CODONS[aa])}"
            )

    def test_each_codon_appears_exactly_once(self):
        """No codon should be listed under two different amino acids."""
        all_codons: list[str] = []
        for codons in AA_TO_CODONS.values():
            all_codons.extend(codons)
        assert len(all_codons) == len(set(all_codons)), "Duplicate codons across amino acids"


# ═══════════════════════════════════════════════════════════════════════════
# 3. BLOSUM62
# ═══════════════════════════════════════════════════════════════════════════

class TestBLOSUM62:
    """Verify BLOSUM62 is a valid 20x20 substitution matrix."""

    def test_is_20x20(self):
        """Matrix must have exactly 20 rows and 20 columns."""
        assert len(BLOSUM62) == 20
        for aa, row in BLOSUM62.items():
            assert len(row) == 20, f"Row for {aa} has {len(row)} columns, expected 20"

    def test_all_standard_aas_as_keys(self):
        """Row keys must be the 20 standard amino acids (BLOSUM index order)."""
        expected = set(STANDARD_AAS_BLOSUM_ORDER)
        actual = set(BLOSUM62.keys())
        assert actual == expected, f"Missing: {expected - actual}, Extra: {actual - expected}"

    def test_all_standard_aas_as_column_keys(self):
        """Each row's column keys must also be the 20 standard AAs."""
        expected = set(STANDARD_AAS_BLOSUM_ORDER)
        for aa, row in BLOSUM62.items():
            assert set(row.keys()) == expected, f"Row {aa} missing column AAs"

    def test_diagonal_scores_positive(self):
        """Diagonal (self-substitution) scores should all be positive."""
        for aa in BLOSUM62:
            assert BLOSUM62[aa][aa] > 0, f"Diagonal score for {aa} should be positive"

    def test_symmetric(self):
        """BLOSUM62 must be symmetric: B[i][j] == B[j][i]."""
        for a1 in BLOSUM62:
            for a2 in BLOSUM62:
                assert BLOSUM62[a1][a2] == BLOSUM62[a2][a1], (
                    f"BLOSUM62[{a1}][{a2}]={BLOSUM62[a1][a2]} != BLOSUM62[{a2}][{a1}]={BLOSUM62[a2][a1]}"
                )

    def test_scores_are_integers(self):
        """All scores in the BLOSUM62 matrix must be integers."""
        for a1, row in BLOSUM62.items():
            for a2, score in row.items():
                assert isinstance(score, int), f"BLOSUM62[{a1}][{a2}]={score!r} is not int"

    def test_score_range(self):
        """BLOSUM62 scores typically range from -4 to +11."""
        all_scores = [
            BLOSUM62[a1][a2]
            for a1 in BLOSUM62
            for a2 in BLOSUM62
        ]
        assert min(all_scores) >= -4, f"Min score {min(all_scores)} below expected range"
        assert max(all_scores) <= 11, f"Max score {max(all_scores)} above expected range"

    def test_known_values(self):
        """Spot-check a few well-known BLOSUM62 entries."""
        assert BLOSUM62["A"]["A"] == 4
        assert BLOSUM62["W"]["W"] == 11
        assert BLOSUM62["C"]["C"] == 9
        assert BLOSUM62["A"]["R"] == -1
        assert BLOSUM62["W"]["A"] == -3  # Symmetric with A-W


# ═══════════════════════════════════════════════════════════════════════════
# 4. HYDROPATHY
# ═══════════════════════════════════════════════════════════════════════════

class TestHydropathy:
    """Verify HYDROPATHY scale has entries for all 20 standard AAs."""

    def test_has_20_entries(self):
        """KYT_DOOLITTLE scale must cover all 20 standard amino acids."""
        assert len(HYDROPATHY) == 20

    def test_all_standard_aas_present(self):
        """Keys must be exactly the 20 standard amino acids."""
        standard = set(STANDARD_AAS)
        actual = set(HYDROPATHY.keys())
        assert actual == standard

    def test_values_are_numeric(self):
        """All hydropathy values must be floats (or int-coercible)."""
        for aa, val in HYDROPATHY.items():
            assert isinstance(val, (int, float)), f"HYDROPATHY[{aa}]={val!r} is not numeric"

    def test_value_range(self):
        """Kyte-Doolittle values range from about -4.5 to +4.5."""
        values = list(HYDROPATHY.values())
        assert min(values) >= -5.0, f"Min {min(values)} unexpectedly low"
        assert max(values) <= 5.0, f"Max {max(values)} unexpectedly high"

    def test_isoleucine_most_hydrophobic(self):
        """Isoleucine (I) should have the highest hydropathy score."""
        max_aa = max(HYDROPATHY, key=HYDROPATHY.get)  # type: ignore[arg-type]
        assert max_aa == "I", f"Expected I as most hydrophobic, got {max_aa}"

    def test_arginine_least_hydrophobic(self):
        """Arginine (R) should have the lowest (most negative) hydropathy score."""
        min_aa = min(HYDROPATHY, key=HYDROPATHY.get)  # type: ignore[arg-type]
        assert min_aa == "R", f"Expected R as least hydrophobic, got {min_aa}"

    def test_known_values(self):
        """Spot-check a few well-known Kyte-Doolittle values."""
        assert HYDROPATHY["I"] == pytest.approx(4.5)
        assert HYDROPATHY["A"] == pytest.approx(1.8)
        assert HYDROPATHY["R"] == pytest.approx(-4.5)
        assert HYDROPATHY["G"] == pytest.approx(-0.4)

    def test_hydrophobic_aas_consistent(self):
        """HYDROPHOBIC_AAS must contain exactly the AAs with hydropathy > 1.0."""
        expected = {aa for aa, val in HYDROPATHY.items() if val > 1.0}
        assert HYDROPHOBIC_AAS == expected


# ═══════════════════════════════════════════════════════════════════════════
# 5. RESTRICTION_ENZYMES
# ═══════════════════════════════════════════════════════════════════════════

class TestRestrictionEnzymes:
    """Verify RESTRICTION_ENZYMES contains known enzymes with valid sites."""

    def test_is_non_empty(self):
        """RESTRICTION_ENZYMES must contain at least one enzyme."""
        assert len(RESTRICTION_ENZYMES) > 0

    def test_has_common_enzymes(self):
        """Essential cloning enzymes must be present."""
        essential = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
        for enzyme in essential:
            assert enzyme in RESTRICTION_ENZYMES, f"Missing essential enzyme: {enzyme}"

    def test_enzyme_names_are_strings(self):
        """All enzyme names must be strings."""
        for name in RESTRICTION_ENZYMES:
            assert isinstance(name, str), f"Enzyme name {name!r} is not a string"

    def test_sites_are_strings(self):
        """All recognition sites must be strings."""
        for name, site in RESTRICTION_ENZYMES.items():
            assert isinstance(site, str), f"Site for {name} is {type(site)}, expected str"

    def test_sites_uppercase(self):
        """All recognition sites must be uppercase (standard convention)."""
        for name, site in RESTRICTION_ENZYMES.items():
            # Allow IUPAC codes (including N), all should be uppercase
            assert site == site.upper(), f"Site for {name} is not uppercase: {site!r}"

    def test_sites_are_palindromic_or_wildcard(self):
        """Most type II restriction sites are palindromic; some use IUPAC wildcards.

        We verify that standard (non-wildcard) Type IIP sites are palindromic,
        while wildcard-containing sites and Type IIS sites are exempted.
        Type IIS enzymes recognise non-palindromic sites and cut downstream.
        """
        # Type IIS enzymes have non-palindromic ACGT-only recognition sites
        TYPE_IIS_ENZYMES = {"BsaI", "BsmBI", "BbsI", "SapI"}
        for name, site in RESTRICTION_ENZYMES.items():
            # Sites with IUPAC ambiguity codes (e.g. N) are not strictly palindromic
            if any(b not in "ACGT" for b in site):
                continue  # Wildcard sites exempted
            if name in TYPE_IIS_ENZYMES:
                continue  # Type IIS sites are non-palindromic by design
            rc = reverse_complement(site)
            assert site == rc, f"Non-wildcard site {name}={site} is not palindromic (RC={rc})"

    def test_known_site_sequences(self):
        """Spot-check recognition sequences for well-known enzymes."""
        known = {
            "EcoRI": "GAATTC",
            "BamHI": "GGATCC",
            "XhoI": "CTCGAG",
            "HindIII": "AAGCTT",
            "NotI": "GCGGCCGC",
            "XbaI": "TCTAGA",
            "NdeI": "CATATG",
        }
        for enzyme, expected_site in known.items():
            assert RESTRICTION_ENZYMES[enzyme] == expected_site, (
                f"{enzyme} site mismatch: got {RESTRICTION_ENZYMES[enzyme]}, expected {expected_site}"
            )

    def test_site_lengths_reasonable(self):
        """Recognition sites should be between 4 and 16 bases.

        Wildcard-containing sites (e.g. SfiI with GGCCNNNNNGGCC, 13 chars)
        can be longer than typical 4-12 bp palindromic sites.
        """
        for name, site in RESTRICTION_ENZYMES.items():
            assert 4 <= len(site) <= 16, f"Site for {name} has unusual length {len(site)}: {site}"

    def test_at_least_20_enzymes(self):
        """Extended REBASE list should contain at least 20 enzymes."""
        assert len(RESTRICTION_ENZYMES) >= 20


# ═══════════════════════════════════════════════════════════════════════════
# 6. reverse_complement()
# ═══════════════════════════════════════════════════════════════════════════

class TestReverseComplement:
    """Verify reverse_complement() produces correct results."""

    def test_simple_sequence(self):
        """Basic reverse complement of a known sequence."""
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        """A palindromic sequence is its own reverse complement."""
        assert reverse_complement("GAATTC") == "GAATTC"  # EcoRI site

    def test_all_bases(self):
        """Each individual base complements correctly."""
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("C") == "G"
        assert reverse_complement("G") == "C"

    def test_lowercase_input(self):
        """Lowercase input should be complemented to lowercase."""
        assert reverse_complement("atcg") == "cgat"

    def test_mixed_case(self):
        """Mixed case input preserves case in output."""
        assert reverse_complement("AtCg") == "cGaT"

    def test_empty_string(self):
        """Empty sequence returns empty reverse complement."""
        assert reverse_complement("") == ""

    def test_long_sequence(self):
        """Test a longer, realistic DNA sequence."""
        seq = "ATGGTTTCTAAAGGTGAA"
        rc = reverse_complement(seq)
        # Double reverse complement should give back the original
        assert reverse_complement(rc) == seq

    def test_double_reverse_identity(self):
        """Reverse complement of reverse complement returns the original."""
        import itertools
        # Test all 256 possible 4-base sequences
        for combo in itertools.product("ACGT", repeat=4):
            seq = "".join(combo)
            assert reverse_complement(reverse_complement(seq)) == seq

    def test_iupac_ambiguity_codes(self):
        """IUPAC ambiguity codes should be complemented correctly."""
        assert reverse_complement("R") == "Y"   # purine -> pyrimidine
        assert reverse_complement("Y") == "R"   # pyrimidine -> purine
        assert reverse_complement("S") == "S"   # strong -> strong (self-complement)
        assert reverse_complement("W") == "W"   # weak -> weak (self-complement)
        assert reverse_complement("K") == "M"   # keto -> amino
        assert reverse_complement("M") == "K"   # amino -> keto
        assert reverse_complement("B") == "V"   # not-A -> not-T
        assert reverse_complement("V") == "B"   # not-T -> not-A
        assert reverse_complement("D") == "H"   # not-C -> not-G
        assert reverse_complement("H") == "D"   # not-G -> not-C
        assert reverse_complement("N") == "N"   # any -> any (self-complement)

    def test_unknown_base_raises_valueerror(self):
        """Unknown characters should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown base"):
            reverse_complement("Z")

    def test_unknown_base_error_message(self):
        """Error message should mention the bad character and supported bases."""
        with pytest.raises(ValueError, match="X") as exc_info:
            reverse_complement("ATXG")
        assert "Unknown base" in str(exc_info.value)
        assert "IUPAC" in str(exc_info.value)

    def test_restriction_site_palindromes(self):
        """Type IIP restriction sites (non-wildcard) are palindromic: RC == self.

        Type IIS enzymes (BsaI, BsmBI, BbsI, SapI) are excluded because
        they recognise non-palindromic sequences and cut downstream.
        """
        TYPE_IIS_ENZYMES = {"BsaI", "BsmBI", "BbsI", "SapI"}
        for name, site in RESTRICTION_ENZYMES.items():
            if any(b not in "ACGT" for b in site):
                continue  # Skip wildcard sites like SfiI
            if name in TYPE_IIS_ENZYMES:
                continue  # Type IIS sites are non-palindromic by design
            rc = reverse_complement(site)
            assert rc == site, f"{name} site {site} RC={rc} is not palindromic"

    def test_complement_dict_consistency(self):
        """COMPLEMENT dict must be the inverse of itself for standard DNA bases.

        RNA uracil (U/u) is excluded: U complements to A, but A complements
        to T (not U) because A-T is the DNA base pair. The gap character
        '-' is also excluded as it maps to itself but is not a base.
        """
        SKIP_BASES = {"U", "u", "-"}
        for base, comp in COMPLEMENT.items():
            if base in SKIP_BASES:
                continue
            assert COMPLEMENT[comp] == base, (
                f"COMPLEMENT[{base!r}]={comp!r}, but COMPLEMENT[{comp!r}]={COMPLEMENT[comp]!r}"
            )
