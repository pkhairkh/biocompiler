"""
Tests for biocompiler.iupac — IUPAC ambiguous base support.

Covers:
  - IUPAC_DNA dictionary correctness
  - is_ambiguous() single-base check
  - has_ambiguous() sequence-level check
  - validate_iupac_sequence() validation
  - resolve_ambiguous() with all strategies
  - expand_ambiguous() combinatorial expansion
  - Integration: validate_dna_sequence with allow_iupac flag
  - Integration: VALID_IUPAC_BASES in constants
"""

import sys
import importlib
import pytest


# Import directly from the submodules to avoid the top-level __init__.py
# circular import chain (type_system <-> sliding_gc).
from biocompiler.iupac import (
    IUPAC_DNA,
    is_ambiguous,
    has_ambiguous,
    validate_iupac_sequence,
    resolve_ambiguous,
    expand_ambiguous,
    VALID_IUPAC_BASES,
)
from biocompiler.constants import VALID_IUPAC_BASES as CONSTANTS_VALID_IUPAC_BASES


# ==============================================================================
# IUPAC_DNA dictionary
# ==============================================================================

class TestIUPACDNA:
    """Test the IUPAC_DNA dictionary is complete and correct."""

    def test_concrete_bases_are_singletons(self):
        for base in "ACGT":
            assert IUPAC_DNA[base] == {base}

    def test_purine(self):
        assert IUPAC_DNA["R"] == {"A", "G"}

    def test_pyrimidine(self):
        assert IUPAC_DNA["Y"] == {"C", "T"}

    def test_strong(self):
        assert IUPAC_DNA["S"] == {"G", "C"}

    def test_weak(self):
        assert IUPAC_DNA["W"] == {"A", "T"}

    def test_keto(self):
        assert IUPAC_DNA["K"] == {"G", "T"}

    def test_amino(self):
        assert IUPAC_DNA["M"] == {"A", "C"}

    def test_not_a(self):
        assert IUPAC_DNA["B"] == {"C", "G", "T"}

    def test_not_c(self):
        assert IUPAC_DNA["D"] == {"A", "G", "T"}

    def test_not_g(self):
        assert IUPAC_DNA["H"] == {"A", "C", "T"}

    def test_not_t(self):
        assert IUPAC_DNA["V"] == {"A", "C", "G"}

    def test_any(self):
        assert IUPAC_DNA["N"] == {"A", "C", "G", "T"}

    def test_all_15_codes_present(self):
        assert len(IUPAC_DNA) == 15

    def test_all_values_are_sets_of_acgt(self):
        for code, bases in IUPAC_DNA.items():
            assert isinstance(bases, set)
            assert bases.issubset({"A", "C", "G", "T"}), f"{code} has invalid bases: {bases}"


# ==============================================================================
# is_ambiguous
# ==============================================================================

class TestIsAmbiguous:
    def test_concrete_bases_are_not_ambiguous(self):
        for base in "ACGT":
            assert not is_ambiguous(base), f"{base} should not be ambiguous"

    def test_n_is_ambiguous(self):
        assert is_ambiguous("N")

    def test_r_is_ambiguous(self):
        assert is_ambiguous("R")

    def test_all_ambiguous_codes(self):
        for code in "RYSWKMBDHVN":
            assert is_ambiguous(code), f"{code} should be ambiguous"

    def test_case_insensitive(self):
        assert is_ambiguous("r")
        assert is_ambiguous("n")

    def test_invalid_base_is_not_ambiguous(self):
        assert not is_ambiguous("X")
        assert not is_ambiguous("Z")
        assert not is_ambiguous("1")


# ==============================================================================
# has_ambiguous
# ==============================================================================

class TestHasAmbiguous:
    def test_concrete_sequence(self):
        assert not has_ambiguous("ACGTACGT")

    def test_single_n(self):
        assert has_ambiguous("ACGNACGT")

    def test_all_ambiguous(self):
        assert has_ambiguous("RYSWKMBD")

    def test_empty_string(self):
        assert not has_ambiguous("")

    def test_case_insensitive(self):
        assert has_ambiguous("acgn")
        assert not has_ambiguous("acgt")


# ==============================================================================
# validate_iupac_sequence
# ==============================================================================

class TestValidateIUPACSequence:
    def test_valid_concrete(self):
        assert validate_iupac_sequence("ACGT") == "ACGT"

    def test_valid_mixed_iupac(self):
        assert validate_iupac_sequence("ATGRYSW") == "ATGRYSW"

    def test_case_insensitive(self):
        assert validate_iupac_sequence("atgr") == "ATGR"

    def test_strips_whitespace(self):
        assert validate_iupac_sequence("  ACGT  ") == "ACGT"

    def test_invalid_char_raises(self):
        with pytest.raises(ValueError, match="Invalid DNA bases"):
            validate_iupac_sequence("ACGTX")

    def test_numbers_raise(self):
        with pytest.raises(ValueError):
            validate_iupac_sequence("ACG123")

    def test_all_valid_iupac_codes(self):
        seq = "ACGTRYSWKMBDHVN"
        assert validate_iupac_sequence(seq) == seq


# ==============================================================================
# resolve_ambiguous
# ==============================================================================

class TestResolveAmbiguous:
    def test_no_ambiguous_bases(self):
        assert resolve_ambiguous("ACGT") == "ACGT"

    def test_most_common_strategy(self):
        # R = A or G, default freq A=0.295 > G=0.205, so R -> A
        result = resolve_ambiguous("R", strategy="most_common")
        assert result == "A"

    def test_most_common_y_resolves_to_t(self):
        # Y = C or T, default freq T=0.295 > C=0.205, so Y -> T
        result = resolve_ambiguous("Y", strategy="most_common")
        assert result == "T"

    def test_most_common_n_resolves_to_a_or_t(self):
        # N = any, A=0.295, T=0.295 (tied); max() picks A (alphabetical)
        result = resolve_ambiguous("N", strategy="most_common")
        assert result in ("A", "T")

    def test_custom_base_freq(self):
        # With G highest frequency, R should resolve to G
        freq = {"A": 0.1, "C": 0.1, "G": 0.7, "T": 0.1}
        result = resolve_ambiguous("R", strategy="most_common", base_freq=freq)
        assert result == "G"

    def test_first_strategy(self):
        # R = {A, G}, first alphabetically is A
        assert resolve_ambiguous("R", strategy="first") == "A"
        assert resolve_ambiguous("Y", strategy="first") == "C"
        assert resolve_ambiguous("N", strategy="first") == "A"

    def test_gc_balanced_strategy(self):
        # For a sequence that's all GC-balanced targets, should pick to keep GC ~50%
        result = resolve_ambiguous("R", strategy="gc_balanced", gc_target=0.5)
        # R = {A, G}; A is 0% GC, G is 100% GC
        # With only one base, GC target 0.5 means both are equally far, picks first
        assert result in ("A", "G")

    def test_gc_balanced_for_long_sequence(self):
        # A sequence of all R with GC target 0.5 should alternate or lean toward G
        seq = "RRRRRRRR"
        result = resolve_ambiguous(seq, strategy="gc_balanced", gc_target=0.5)
        # Result should be all ACGT
        assert all(c in "ACGT" for c in result)
        assert len(result) == 8

    def test_cai_optimal_strategy_with_table(self):
        # Simple CAI table where ATG has high CAI and ATA has low
        cai_table = {"ATG": 0.9, "ATA": 0.1, "ATC": 0.5, "ATT": 0.3}
        # R at position 2 of "ATR" -> ATR could be ATA or ATG
        result = resolve_ambiguous("ATR", strategy="cai_optimal", cai_table=cai_table)
        assert result == "ATG"  # ATG has highest CAI

    def test_cai_optimal_without_table_falls_back(self):
        # Should fall back to most_common when no table provided
        result = resolve_ambiguous("R", strategy="cai_optimal")
        assert result == "A"  # Same as most_common

    def test_mixed_concrete_and_ambiguous(self):
        # ACGT are left unchanged, only RYSW are resolved
        result = resolve_ambiguous("ACGTRYSW", strategy="first")
        assert result.startswith("ACGT")
        assert len(result) == 8
        assert all(c in "ACGT" for c in result)

    def test_unknown_strategy_raises(self):
        with pytest.raises(ValueError, match="Unknown resolution strategy"):
            resolve_ambiguous("R", strategy="invalid_strategy")

    def test_invalid_sequence_raises(self):
        with pytest.raises(ValueError):
            resolve_ambiguous("ACGTX")


# ==============================================================================
# expand_ambiguous
# ==============================================================================

class TestExpandAmbiguous:
    def test_no_ambiguous_returns_single(self):
        result = expand_ambiguous("ACGT")
        assert result == ["ACGT"]

    def test_single_r_expands(self):
        result = expand_ambiguous("R")
        assert sorted(result) == ["A", "G"]

    def test_ry_expands(self):
        result = expand_ambiguous("RY")
        assert sorted(result) == ["AC", "AT", "GC", "GT"]

    def test_n_expands_to_four(self):
        result = expand_ambiguous("N")
        assert sorted(result) == ["A", "C", "G", "T"]

    def test_atr_expands(self):
        result = expand_ambiguous("ATR")
        assert sorted(result) == ["ATA", "ATG"]

    def test_preserves_concrete_bases(self):
        result = expand_ambiguous("ACGNACGT")
        # N expands to 4 possibilities, concrete bases stay
        assert len(result) == 4
        for seq in result:
            assert seq.startswith("ACG")
            assert seq.endswith("ACGT")

    def test_all_15_codes_expand(self):
        # Each code should expand correctly
        for code, bases in IUPAC_DNA.items():
            if len(bases) == 1:
                result = expand_ambiguous(code)
                assert result == [code]
            else:
                result = expand_ambiguous(code)
                assert len(result) == len(bases)
                assert set(result) == set(bases)

    def test_length_preserved(self):
        result = expand_ambiguous("ATR")
        for seq in result:
            assert len(seq) == 3

    def test_combinatorial_count(self):
        # R(2) * Y(2) * N(4) = 16
        result = expand_ambiguous("RYN")
        assert len(result) == 16

    def test_cap_on_large_expansion(self):
        # Create a sequence that expands to more than 4096 combos
        # NNNNNNN = 4^7 = 16384 > 4096
        result = expand_ambiguous("NNNNNNN")
        assert result == []  # Should return empty list due to cap

    def test_invalid_char_raises(self):
        with pytest.raises(ValueError):
            expand_ambiguous("ACGTX")


# ==============================================================================
# Integration: constants.VALID_IUPAC_BASES
# ==============================================================================

class TestConstantsIntegration:
    def test_valid_iupac_bases_in_constants(self):
        """VALID_IUPAC_BASES should be available from constants.py."""
        assert CONSTANTS_VALID_IUPAC_BASES is not None

    def test_valid_iupac_bases_contains_all_codes(self):
        """VALID_IUPAC_BASES should contain all 15 IUPAC DNA codes."""
        expected = set("ACGTRYSWKMBDHVN")
        assert CONSTANTS_VALID_IUPAC_BASES == expected

    def test_consistent_with_iupac_module(self):
        """VALID_IUPAC_BASES in constants and iupac module should agree."""
        assert CONSTANTS_VALID_IUPAC_BASES == VALID_IUPAC_BASES


# ==============================================================================
# Integration: scanner.validate_dna_sequence with allow_iupac
# ==============================================================================

class TestScannerIntegration:
    def test_default_rejects_iupac(self):
        """Without allow_iupac, only ACGTN should be accepted."""
        from biocompiler.scanner import validate_dna_sequence
        from biocompiler.exceptions import InvalidSequenceError

        # ACGTN should pass
        assert validate_dna_sequence("ACGTN") == "ACGTN"

        # R should fail without allow_iupac
        with pytest.raises(InvalidSequenceError):
            validate_dna_sequence("ACGR")

    def test_allow_iupac_accepts_all_codes(self):
        """With allow_iupac=True, all IUPAC codes should be accepted."""
        from biocompiler.scanner import validate_dna_sequence

        result = validate_dna_sequence("ACGTRYSWKMBDHVN", allow_iupac=True)
        assert result == "ACGTRYSWKMBDHVN"

    def test_allow_iupac_rejects_invalid(self):
        """Even with allow_iupac=True, invalid chars should be rejected."""
        from biocompiler.scanner import validate_dna_sequence
        from biocompiler.exceptions import InvalidSequenceError

        with pytest.raises(InvalidSequenceError):
            validate_dna_sequence("ACGTX", allow_iupac=True)


# ==============================================================================
# Integration: Top-level package imports
# ==============================================================================

class TestPackageImports:
    def test_iupac_functions_importable(self):
        """IUPAC functions should be importable from top-level package."""
        import biocompiler
        assert biocompiler.resolve_ambiguous is not None
        assert biocompiler.is_ambiguous is not None
        assert biocompiler.expand_ambiguous is not None
        assert biocompiler.has_ambiguous is not None
        assert biocompiler.validate_iupac_sequence is not None
        assert biocompiler.IUPAC_DNA is not None

    def test_iupac_in_all(self):
        """IUPAC functions should be in __all__."""
        import biocompiler
        for name in ("IUPAC_DNA", "resolve_ambiguous", "is_ambiguous",
                      "expand_ambiguous", "has_ambiguous", "validate_iupac_sequence"):
            assert name in biocompiler.__all__


# ==============================================================================
# Property-based invariants
# ==============================================================================

class TestInvariants:
    def test_resolve_produces_only_acgt(self):
        """Resolved sequences should contain only ACGT."""
        for code in "RYSWKMBDHVN":
            result = resolve_ambiguous(code)
            assert all(c in "ACGT" for c in result)

    def test_expand_results_are_subsets(self):
        """Each expanded sequence should be a possible concrete version."""
        seq = "ATRSWN"
        expanded = expand_ambiguous(seq)
        for concrete in expanded:
            for i, base in enumerate(concrete):
                assert base in IUPAC_DNA[seq[i]], \
                    f"Base {base} at position {i} not in IUPAC expansion of {seq[i]}"

    def test_resolve_does_not_change_length(self):
        """Resolving should not change sequence length."""
        for seq in ["ACGT", "RYSWKM", "ACGNRYSW", "NNNN"]:
            result = resolve_ambiguous(seq)
            assert len(result) == len(seq)

    def test_expand_preserves_length(self):
        """All expanded sequences should have same length as input."""
        for seq in ["ACGT", "RY", "ATR", "NNN"]:
            expanded = expand_ambiguous(seq)
            for concrete in expanded:
                assert len(concrete) == len(seq)

    def test_resolve_deterministic(self):
        """Same input should always produce same output."""
        seq = "ACGTRYSWKMBDHVN"
        results = [resolve_ambiguous(seq) for _ in range(10)]
        assert len(set(results)) == 1  # All identical

    def test_complement_of_iupac_code_is_valid(self):
        """The complement of each IUPAC code should also be a valid IUPAC code."""
        from biocompiler.constants import COMPLEMENT
        for code in IUPAC_DNA:
            complement = COMPLEMENT.get(code.upper())
            if complement is not None:
                assert complement.upper() in IUPAC_DNA, \
                    f"Complement of {code} ({complement}) not in IUPAC_DNA"
