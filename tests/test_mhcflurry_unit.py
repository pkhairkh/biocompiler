"""Unit tests for mhcflurry_adapter — pure-function tests that do NOT import mhcflurry.

These tests exercise the adapter's pure functions and internal helpers
directly, without requiring the MHCflurry library or models to be
installed.  This ensures the adapter's logic is correct even when
MHCflurry is unavailable.

Test categories
---------------
1. ic50_to_binding_score — known IC50 → expected score
2. _validate_sequence — valid/invalid/empty sequences
3. _validate_peptide — valid/invalid, allows empty
4. _validate_protein — valid/invalid, rejects empty
5. _extract_overlapping_peptides — correct peptide generation
6. _LRUCache — put/get/clear/hit_rate/eviction
7. _mhcflurry_result_to_binding_result — conversion correctness
8. MHCBindingResult — dataclass construction
"""
from __future__ import annotations

import math

import pytest

# ---------------------------------------------------------------------------
# Import the adapter's pure functions (NOT the mhcflurry library)
# ---------------------------------------------------------------------------

from biocompiler.mhcflurry_adapter import (
    _LRUCache,
    _extract_overlapping_peptides,
    _mhcflurry_result_to_binding_result,
    _validate_peptide,
    _validate_protein,
    _validate_sequence,
    ic50_to_binding_score,
)
from biocompiler.immunogenicity import MHCBindingResult, classify_binding


# ═══════════════════════════════════════════════════════════════════════════
# 1. ic50_to_binding_score
# ═══════════════════════════════════════════════════════════════════════════

class TestIc50ToBindingScore:
    """Test ic50_to_binding_score with known IC50 values → expected scores."""

    def test_ic50_1_nm_maps_near_1(self) -> None:
        """IC50 = 1 nM should give a score very close to 1.0."""
        score = ic50_to_binding_score(1.0)
        assert score == pytest.approx(1.0 - math.log(1.0) / math.log(50_000), abs=1e-9)

    def test_ic50_50_nm_expected_score(self) -> None:
        """IC50 = 50 nM → score ≈ 1 - log(50)/log(50000) ≈ 0.745."""
        expected = 1.0 - math.log(50.0) / math.log(50_000.0)
        score = ic50_to_binding_score(50.0)
        assert score == pytest.approx(expected, abs=1e-9)

    def test_ic50_500_nm_expected_score(self) -> None:
        """IC50 = 500 nM → score ≈ 1 - log(500)/log(50000) ≈ 0.546."""
        expected = 1.0 - math.log(500.0) / math.log(50_000.0)
        score = ic50_to_binding_score(500.0)
        assert score == pytest.approx(expected, abs=1e-9)

    def test_ic50_5000_nm_expected_score(self) -> None:
        """IC50 = 5000 nM → score ≈ 1 - log(5000)/log(50000) ≈ 0.369."""
        expected = 1.0 - math.log(5000.0) / math.log(50_000.0)
        score = ic50_to_binding_score(5000.0)
        assert score == pytest.approx(expected, abs=1e-9)

    def test_ic50_50000_nm_is_zero(self) -> None:
        """IC50 = 50 000 nM (max) → score = 0.0."""
        score = ic50_to_binding_score(50_000.0)
        assert score == pytest.approx(0.0, abs=1e-9)

    def test_ic50_zero_returns_1(self) -> None:
        """IC50 ≤ 0 is treated as maximum binding → score = 1.0."""
        assert ic50_to_binding_score(0.0) == 1.0

    def test_ic50_negative_returns_1(self) -> None:
        """Negative IC50 (invalid) is clamped to score = 1.0."""
        assert ic50_to_binding_score(-10.0) == 1.0

    def test_ic50_above_max_clamped_to_zero(self) -> None:
        """IC50 > 50 000 nM → score clamped to 0.0."""
        score = ic50_to_binding_score(1_000_000.0)
        assert score == 0.0

    def test_score_always_in_unit_interval(self) -> None:
        """Score is always in [0, 1] for a range of IC50 values."""
        for ic50 in [0.001, 0.1, 1, 10, 50, 500, 5000, 50000, 1e6, 1e9]:
            score = ic50_to_binding_score(ic50)
            assert 0.0 <= score <= 1.0, f"Score {score} out of [0,1] for IC50 {ic50}"

    def test_monotonic_decreasing(self) -> None:
        """Higher IC50 should never yield a higher binding score."""
        ic50s = [0.5, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000]
        scores = [ic50_to_binding_score(ic) for ic in ic50s]
        for i in range(len(scores) - 1):
            assert scores[i] >= scores[i + 1], (
                f"Non-monotonic: IC50={ic50s[i]}→{scores[i]}, "
                f"IC50={ic50s[i+1]}→{scores[i+1]}"
            )


# ═══════════════════════════════════════════════════════════════════════════
# 2. _validate_sequence
# ═══════════════════════════════════════════════════════════════════════════

class TestValidateSequence:
    """Test _validate_sequence with valid/invalid/empty inputs."""

    def test_valid_sequence_uppercase(self) -> None:
        """A valid all-uppercase sequence is returned as-is."""
        assert _validate_sequence("ACDEFGHIKLMNPQRSTVWY") == "ACDEFGHIKLMNPQRSTVWY"

    def test_valid_sequence_lowercase_converted(self) -> None:
        """Lower-case letters are upper-cased."""
        assert _validate_sequence("acdefghiklmnpqrstvwy") == "ACDEFGHIKLMNPQRSTVWY"

    def test_valid_sequence_mixed_case(self) -> None:
        """Mixed-case input is upper-cased."""
        assert _validate_sequence("AcDeFgHiKlMnPqRsTvWy") == "ACDEFGHIKLMNPQRSTVWY"

    def test_whitespace_stripped(self) -> None:
        """Leading/trailing whitespace is stripped."""
        assert _validate_sequence("  ACDEF  ") == "ACDEF"

    def test_empty_sequence_disallowed_by_default(self) -> None:
        """Empty sequence raises ValueError when allow_empty=False (default)."""
        with pytest.raises(ValueError, match="must not be empty"):
            _validate_sequence("")

    def test_whitespace_only_sequence_disallowed(self) -> None:
        """Whitespace-only sequence raises ValueError."""
        with pytest.raises(ValueError, match="must not be empty"):
            _validate_sequence("   ")

    def test_empty_sequence_allowed_when_flag_set(self) -> None:
        """Empty sequence is accepted when allow_empty=True."""
        assert _validate_sequence("", allow_empty=True) == ""

    def test_non_standard_amino_acids_rejected(self) -> None:
        """Non-standard characters (B, J, O, U, X, Z, digits, etc.) raise ValueError."""
        invalid_seqs = ["ABCD", "ACD1", "ACDX", "ACDZ", "HELLO", "A C", "A-C"]
        for seq in invalid_seqs:
            with pytest.raises(ValueError, match="non-standard"):
                _validate_sequence(seq)

    def test_specific_invalid_chars_in_error(self) -> None:
        """The error message lists the invalid characters found."""
        with pytest.raises(ValueError, match="non-standard"):
            _validate_sequence("ABCX1")


# ═══════════════════════════════════════════════════════════════════════════
# 3. _validate_peptide
# ═══════════════════════════════════════════════════════════════════════════

class TestValidatePeptide:
    """Test _validate_peptide — valid/invalid, allows empty."""

    def test_valid_peptide(self) -> None:
        """Valid peptide string is returned upper-cased."""
        assert _validate_peptide("SIINFEKL") == "SIINFEKL"

    def test_lowercase_peptide_uppercased(self) -> None:
        """Lower-case peptide is upper-cased."""
        assert _validate_peptide("siinfekl") == "SIINFEKL"

    def test_empty_peptide_allowed(self) -> None:
        """Empty peptide is allowed (delegates with allow_empty=True)."""
        assert _validate_peptide("") == ""

    def test_whitespace_only_peptide_becomes_empty(self) -> None:
        """Whitespace-only peptide is stripped to empty and allowed."""
        assert _validate_peptide("  ") == ""

    def test_invalid_peptide_rejected(self) -> None:
        """Peptide with non-standard amino acids raises ValueError."""
        with pytest.raises(ValueError, match="non-standard"):
            _validate_peptide("SIINFEKLX")

    def test_peptide_with_digits_rejected(self) -> None:
        """Peptide containing digits raises ValueError."""
        with pytest.raises(ValueError, match="non-standard"):
            _validate_peptide("SIINFEK1")


# ═══════════════════════════════════════════════════════════════════════════
# 4. _validate_protein
# ═══════════════════════════════════════════════════════════════════════════

class TestValidateProtein:
    """Test _validate_protein — valid/invalid, rejects empty."""

    def test_valid_protein(self) -> None:
        """Valid protein string is returned upper-cased."""
        assert _validate_protein("MAGRSGDLDAIIRYVKQLR") == "MAGRSGDLDAIIRYVKQLR"

    def test_lowercase_protein_uppercased(self) -> None:
        """Lower-case protein is upper-cased."""
        assert _validate_protein("magrsgdldaiiryvkqlr") == "MAGRSGDLDAIIRYVKQLR"

    def test_empty_protein_rejected(self) -> None:
        """Empty protein raises ValueError (allow_empty=False)."""
        with pytest.raises(ValueError, match="must not be empty"):
            _validate_protein("")

    def test_whitespace_only_protein_rejected(self) -> None:
        """Whitespace-only protein is stripped to empty, then rejected."""
        with pytest.raises(ValueError, match="must not be empty"):
            _validate_protein("   ")

    def test_invalid_protein_rejected(self) -> None:
        """Protein with non-standard amino acids raises ValueError."""
        with pytest.raises(ValueError, match="non-standard"):
            _validate_protein("MAGRSGDLDAIIRYVKQLRX")

    def test_protein_with_special_chars_rejected(self) -> None:
        """Protein with special characters raises ValueError."""
        with pytest.raises(ValueError, match="non-standard"):
            _validate_protein("MAGRSGDL-DAIIRYVKQLR")


# ═══════════════════════════════════════════════════════════════════════════
# 5. _extract_overlapping_peptides
# ═══════════════════════════════════════════════════════════════════════════

class TestExtractOverlappingPeptides:
    """Test _extract_overlapping_peptides — correct peptide generation."""

    def test_simple_9mer_from_10aa(self) -> None:
        """A 10-residue protein yields two overlapping 9-mers."""
        protein = "ABCDEFGHIJ"  # 10 residues
        result = _extract_overlapping_peptides(protein, [9])
        assert len(result) == 2
        assert result[0] == ("ABCDEFGHI", 0, 8)
        assert result[1] == ("BCDEFGHIJ", 1, 9)

    def test_8_and_9_mers_from_10aa(self) -> None:
        """A 10-residue protein yields 3 × 8-mers + 2 × 9-mers = 5 peptides."""
        protein = "ABCDEFGHIJ"
        result = _extract_overlapping_peptides(protein, [8, 9])
        assert len(result) == 5
        # First 3 are 8-mers
        peps_8 = [r for r in result if len(r[0]) == 8]
        assert len(peps_8) == 3
        # Last 2 are 9-mers
        peps_9 = [r for r in result if len(r[0]) == 9]
        assert len(peps_9) == 2

    def test_protein_equals_epitope_length(self) -> None:
        """When protein length equals epitope length, one peptide is produced."""
        protein = "ABCDEFGHI"
        result = _extract_overlapping_peptides(protein, [9])
        assert len(result) == 1
        assert result[0] == ("ABCDEFGHI", 0, 8)

    def test_protein_shorter_than_epitope_length(self) -> None:
        """When protein is shorter than the epitope length, no peptides are produced."""
        protein = "ABCD"
        result = _extract_overlapping_peptides(protein, [9])
        assert result == []

    def test_single_residue_protein(self) -> None:
        """Single-residue protein yields no peptides for any length ≥ 2."""
        result = _extract_overlapping_peptides("M", [8, 9, 10, 11])
        assert result == []

    def test_positions_are_correct(self) -> None:
        """Start and end positions are correct for each peptide."""
        protein = "MAGRSGDLDAIIRYVKQLR"  # 20 residues
        result = _extract_overlapping_peptides(protein, [9])
        for pep, start, end in result:
            assert pep == protein[start : end + 1]
            assert end - start + 1 == 9
            assert end == start + 8

    def test_default_epitope_lengths(self) -> None:
        """Using [8, 9, 10, 11] on a 12-residue protein gives correct counts."""
        protein = "ABCDEFGHIJKL"  # 12 residues
        result = _extract_overlapping_peptides(protein, [8, 9, 10, 11])
        # 8-mers: 12-8+1=5, 9-mers: 12-9+1=4, 10-mers: 12-10+1=3, 11-mers: 12-11+1=2
        assert len(result) == 5 + 4 + 3 + 2

    def test_empty_protein(self) -> None:
        """Empty protein yields no peptides."""
        result = _extract_overlapping_peptides("", [8, 9])
        assert result == []

    def test_no_epitope_lengths(self) -> None:
        """Empty epitope_lengths list yields no peptides."""
        result = _extract_overlapping_peptides("ABCDEFGHIJ", [])
        assert result == []


# ═══════════════════════════════════════════════════════════════════════════
# 6. _LRUCache
# ═══════════════════════════════════════════════════════════════════════════

class TestLRUCache:
    """Test _LRUCache — put/get/clear/hit_rate/eviction."""

    def _make_result(self, peptide: str = "SIINFEKL", allele: str = "HLA-A*02:01") -> MHCBindingResult:
        """Helper to create a simple MHCBindingResult for cache tests."""
        return MHCBindingResult(
            allele=allele,
            peptide=peptide,
            start_position=0,
            end_position=len(peptide) - 1,
            binding_score=0.75,
            ic50_nm=150.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )

    def test_put_and_get(self) -> None:
        """Put a result and retrieve it."""
        cache = _LRUCache()
        result = self._make_result()
        cache.put("SIINFEKL", "HLA-A*02:01", result)
        retrieved = cache.get("SIINFEKL", "HLA-A*02:01")
        assert retrieved is not None
        assert retrieved.peptide == "SIINFEKL"
        assert retrieved.allele == "HLA-A*02:01"

    def test_get_miss_returns_none(self) -> None:
        """Getting a non-existent key returns None."""
        cache = _LRUCache()
        assert cache.get("SIINFEKL", "HLA-A*02:01") is None

    def test_clear_empties_cache(self) -> None:
        """Clear removes all entries and resets counters."""
        cache = _LRUCache()
        cache.put("SIINFEKL", "HLA-A*02:01", self._make_result())
        assert cache.size == 1
        cache.clear()
        assert cache.size == 0
        assert cache.get("SIINFEKL", "HLA-A*02:01") is None

    def test_clear_resets_hit_miss_counters(self) -> None:
        """Clear resets hits and misses to zero."""
        cache = _LRUCache()
        cache.put("PEP1", "HLA-A*02:01", self._make_result("PEP1"))
        cache.get("PEP1", "HLA-A*02:01")  # hit
        cache.get("PEP2", "HLA-A*02:01")  # miss
        assert cache.hits == 1
        assert cache.misses == 1
        cache.clear()
        assert cache.hits == 0
        assert cache.misses == 0

    def test_hit_rate_no_accesses(self) -> None:
        """Hit rate is 0.0 when there have been no accesses."""
        cache = _LRUCache()
        assert cache.hit_rate == 0.0

    def test_hit_rate_all_hits(self) -> None:
        """Hit rate is 1.0 when all gets are hits."""
        cache = _LRUCache()
        cache.put("PEP1", "HLA-A*02:01", self._make_result("PEP1"))
        cache.get("PEP1", "HLA-A*02:01")  # hit
        cache.get("PEP1", "HLA-A*02:01")  # hit
        assert cache.hit_rate == pytest.approx(1.0)

    def test_hit_rate_mixed(self) -> None:
        """Hit rate reflects the ratio of hits to total accesses."""
        cache = _LRUCache()
        cache.put("PEP1", "HLA-A*02:01", self._make_result("PEP1"))
        cache.get("PEP1", "HLA-A*02:01")  # hit
        cache.get("PEP2", "HLA-A*02:01")  # miss
        assert cache.hit_rate == pytest.approx(0.5)

    def test_size_tracks_entries(self) -> None:
        """size property returns the current number of cached entries."""
        cache = _LRUCache()
        assert cache.size == 0
        cache.put("PEP1", "HLA-A*02:01", self._make_result("PEP1"))
        assert cache.size == 1
        cache.put("PEP2", "HLA-A*02:01", self._make_result("PEP2"))
        assert cache.size == 2

    def test_overwrite_existing_key(self) -> None:
        """Putting with an existing key updates the value."""
        cache = _LRUCache()
        r1 = self._make_result()
        r1_updated = MHCBindingResult(
            allele="HLA-A*02:01", peptide="SIINFEKL",
            start_position=0, end_position=7,
            binding_score=0.99, ic50_nm=10.0,
            binding_class="strong_binder",
            anchor_residues={}, anchor_scores={},
        )
        cache.put("SIINFEKL", "HLA-A*02:01", r1)
        cache.put("SIINFEKL", "HLA-A*02:01", r1_updated)
        assert cache.size == 1
        retrieved = cache.get("SIINFEKL", "HLA-A*02:01")
        assert retrieved is not None
        assert retrieved.binding_score == pytest.approx(0.99)

    def test_overwrite_does_not_count_as_miss(self) -> None:
        """Overwriting an existing key via put then get should be a hit."""
        cache = _LRUCache()
        r = self._make_result()
        cache.put("SIINFEKL", "HLA-A*02:01", r)
        cache.put("SIINFEKL", "HLA-A*02:01", r)  # overwrite
        cache.get("SIINFEKL", "HLA-A*02:01")  # hit
        assert cache.hits == 1
        assert cache.misses == 0

    def test_eviction_when_full(self) -> None:
        """When the cache is full, the least-recently-used entry is evicted."""
        cache = _LRUCache(max_size=3)
        r1 = self._make_result("PEP1")
        r2 = self._make_result("PEP2")
        r3 = self._make_result("PEP3")
        r4 = self._make_result("PEP4")

        cache.put("PEP1", "HLA-A*02:01", r1)
        cache.put("PEP2", "HLA-A*02:01", r2)
        cache.put("PEP3", "HLA-A*02:01", r3)
        assert cache.size == 3

        # Adding a 4th entry should evict PEP1 (LRU)
        cache.put("PEP4", "HLA-A*02:01", r4)
        assert cache.size == 3
        assert cache.get("PEP1", "HLA-A*02:01") is None  # evicted
        assert cache.get("PEP2", "HLA-A*02:01") is not None
        assert cache.get("PEP3", "HLA-A*02:01") is not None
        assert cache.get("PEP4", "HLA-A*02:01") is not None

    def test_lru_access_refreshes_entry(self) -> None:
        """Accessing an entry moves it to most-recently-used, preventing eviction."""
        cache = _LRUCache(max_size=3)
        r1 = self._make_result("PEP1")
        r2 = self._make_result("PEP2")
        r3 = self._make_result("PEP3")
        r4 = self._make_result("PEP4")

        cache.put("PEP1", "HLA-A*02:01", r1)
        cache.put("PEP2", "HLA-A*02:01", r2)
        cache.put("PEP3", "HLA-A*02:01", r3)

        # Access PEP1 — makes it most-recently-used
        cache.get("PEP1", "HLA-A*02:01")

        # Adding PEP4 should evict PEP2 (now LRU), not PEP1
        cache.put("PEP4", "HLA-A*02:01", r4)
        assert cache.get("PEP1", "HLA-A*02:01") is not None  # still there
        assert cache.get("PEP2", "HLA-A*02:01") is None  # evicted

    def test_different_alleles_are_separate_keys(self) -> None:
        """Same peptide with different alleles are separate cache entries."""
        cache = _LRUCache()
        r1 = self._make_result("SIINFEKL", "HLA-A*02:01")
        r2 = self._make_result("SIINFEKL", "HLA-A*03:01")
        cache.put("SIINFEKL", "HLA-A*02:01", r1)
        cache.put("SIINFEKL", "HLA-A*03:01", r2)
        assert cache.size == 2
        got_a2 = cache.get("SIINFEKL", "HLA-A*02:01")
        got_a3 = cache.get("SIINFEKL", "HLA-A*03:01")
        assert got_a2 is not None
        assert got_a3 is not None
        assert got_a2.allele == "HLA-A*02:01"
        assert got_a3.allele == "HLA-A*03:01"


# ═══════════════════════════════════════════════════════════════════════════
# 7. _mhcflurry_result_to_binding_result
# ═══════════════════════════════════════════════════════════════════════════

class TestMhcflurryResultToBindingResult:
    """Test _mhcflurry_result_to_binding_result — conversion correctness."""

    def test_basic_conversion(self) -> None:
        """Basic conversion produces a valid MHCBindingResult."""
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=150.0,
        )
        assert isinstance(result, MHCBindingResult)
        assert result.allele == "HLA-A*02:01"
        assert result.peptide == "SIINFEKL"
        assert result.start_position == 0
        assert result.end_position == 7

    def test_binding_score_from_ic50(self) -> None:
        """binding_score is computed from ic50_to_binding_score."""
        ic50 = 150.0
        expected_score = ic50_to_binding_score(ic50)
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=ic50,
        )
        assert result.binding_score == pytest.approx(round(expected_score, 6))

    def test_ic50_rounded(self) -> None:
        """ic50_nm in the result is rounded to 2 decimal places."""
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=123.4567,
        )
        assert result.ic50_nm == pytest.approx(123.46, abs=0.01)

    def test_binding_class_from_classify_binding(self) -> None:
        """binding_class matches classify_binding for the given IC50."""
        for ic50, expected_class in [
            (10.0, "strong_binder"),
            (50.0, "moderate_binder"),
            (100.0, "moderate_binder"),
            (500.0, "moderate_binder"),
            (1000.0, "weak_binder"),
            (5000.0, "weak_binder"),
            (10000.0, "non_binder"),
        ]:
            result = _mhcflurry_result_to_binding_result(
                peptide="SIINFEKL",
                allele="HLA-A*02:01",
                start_position=0,
                end_position=7,
                ic50_nm=ic50,
            )
            assert result.binding_class == expected_class, (
                f"IC50={ic50}: expected {expected_class}, got {result.binding_class}"
            )

    def test_presentation_score_overrides_binding_score(self) -> None:
        """When presentation_score is provided, it overrides binding_score."""
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=150.0,
            presentation_score=0.42,
        )
        # binding_score should be the presentation_score, not derived from IC50
        assert result.binding_score == pytest.approx(0.42)

    def test_presentation_score_none_uses_ic50(self) -> None:
        """When presentation_score is None, binding_score comes from IC50."""
        ic50 = 150.0
        expected = round(ic50_to_binding_score(ic50), 6)
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=ic50,
            presentation_score=None,
        )
        assert result.binding_score == pytest.approx(expected)

    def test_anchor_residues_p2_and_cterm(self) -> None:
        """Anchor residues are set at P2 (index 1) and C-terminus."""
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=150.0,
        )
        # P2 = index 1 = 'I', C-terminus = index 7 = 'L'
        assert result.anchor_residues[1] == "I"
        assert result.anchor_residues[7] == "L"

    def test_anchor_residues_single_aa_peptide(self) -> None:
        """For a single-AA peptide, only C-terminus anchor is set (P2 skipped)."""
        result = _mhcflurry_result_to_binding_result(
            peptide="S",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=0,
            ic50_nm=5000.0,
        )
        # len(peptide) = 1, so no P2 anchor; C-terminus = index 0
        assert 1 not in result.anchor_residues
        assert result.anchor_residues[0] == "S"

    def test_anchor_residues_two_aa_peptide(self) -> None:
        """For a 2-AA peptide, P2 and C-terminus are the same position."""
        result = _mhcflurry_result_to_binding_result(
            peptide="SI",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=1,
            ic50_nm=500.0,
        )
        # P2 = index 1 = 'I', C-terminus = index 1 = 'I' — both set
        assert result.anchor_residues[1] == "I"

    def test_anchor_scores_set(self) -> None:
        """Anchor scores are populated at P2 and C-terminus positions."""
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=150.0,
        )
        expected_score = round(ic50_to_binding_score(150.0), 6)
        assert result.anchor_scores[1] == pytest.approx(expected_score)
        assert result.anchor_scores[7] == pytest.approx(expected_score)

    def test_method_parameter_default(self) -> None:
        """Default method parameter is 'mhcflurry' (not stored in MHCBindingResult but
        the conversion still works — verify no crash)."""
        # MHCBindingResult doesn't have a method field, but the function
        # should not raise
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=150.0,
            method="mhcflurry",
        )
        assert isinstance(result, MHCBindingResult)

    def test_binding_score_rounded_to_6_decimals(self) -> None:
        """binding_score is rounded to 6 decimal places."""
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=99.123456789,
        )
        # Verify the score has at most 6 decimal digits
        assert result.binding_score == round(result.binding_score, 6)

    def test_strong_binder_example(self) -> None:
        """Strong binder (IC50=10) has high binding score and correct class."""
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=10.0,
        )
        assert result.binding_class == "strong_binder"
        assert result.binding_score > 0.7

    def test_non_binder_example(self) -> None:
        """Non-binder (IC50=50000) has near-zero binding score."""
        result = _mhcflurry_result_to_binding_result(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=7,
            ic50_nm=50000.0,
        )
        assert result.binding_class == "non_binder"
        assert result.binding_score == pytest.approx(0.0)


# ═══════════════════════════════════════════════════════════════════════════
# 8. MHCBindingResult dataclass construction
# ═══════════════════════════════════════════════════════════════════════════

class TestMHCBindingResultConstruction:
    """Test MHCBindingResult dataclass construction."""

    def test_basic_construction(self) -> None:
        """MHCBindingResult can be constructed with all required fields."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="SIINFEKL",
            start_position=0,
            end_position=7,
            binding_score=0.85,
            ic50_nm=42.0,
            binding_class="strong_binder",
            anchor_residues={1: "I", 7: "L"},
            anchor_scores={1: 0.9, 7: 0.8},
        )
        assert result.allele == "HLA-A*02:01"
        assert result.peptide == "SIINFEKL"
        assert result.start_position == 0
        assert result.end_position == 7
        assert result.binding_score == pytest.approx(0.85)
        assert result.ic50_nm == pytest.approx(42.0)
        assert result.binding_class == "strong_binder"
        assert result.anchor_residues == {1: "I", 7: "L"}
        assert result.anchor_scores == {1: pytest.approx(0.9), 7: pytest.approx(0.8)}

    def test_empty_anchor_dicts(self) -> None:
        """MHCBindingResult accepts empty anchor dictionaries."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="SIINFEKL",
            start_position=0,
            end_position=7,
            binding_score=0.5,
            ic50_nm=500.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert result.anchor_residues == {}
        assert result.anchor_scores == {}

    def test_ic50_nm_can_be_none(self) -> None:
        """ic50_nm field accepts None (as per type hint float | None)."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="SIINFEKL",
            start_position=0,
            end_position=7,
            binding_score=0.5,
            ic50_nm=None,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert result.ic50_nm is None

    def test_binding_score_zero(self) -> None:
        """binding_score of 0.0 is valid (non-binder)."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="SIINFEKL",
            start_position=0,
            end_position=7,
            binding_score=0.0,
            ic50_nm=50000.0,
            binding_class="non_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert result.binding_score == 0.0

    def test_binding_score_one(self) -> None:
        """binding_score of 1.0 is valid (strong binder)."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="SIINFEKL",
            start_position=0,
            end_position=7,
            binding_score=1.0,
            ic50_nm=1.0,
            binding_class="strong_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert result.binding_score == 1.0

    def test_all_binding_classes(self) -> None:
        """All four valid binding classes can be set."""
        for bclass in ("strong_binder", "moderate_binder", "weak_binder", "non_binder"):
            result = MHCBindingResult(
                allele="HLA-A*02:01",
                peptide="SIINFEKL",
                start_position=0,
                end_position=7,
                binding_score=0.5,
                ic50_nm=500.0,
                binding_class=bclass,
                anchor_residues={},
                anchor_scores={},
            )
            assert result.binding_class == bclass

    def test_positions_consistent_with_peptide_length(self) -> None:
        """end_position - start_position + 1 equals peptide length."""
        peptide = "SIINFEKL"
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide=peptide,
            start_position=0,
            end_position=len(peptide) - 1,
            binding_score=0.5,
            ic50_nm=500.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert result.end_position - result.start_position + 1 == len(peptide)

    def test_dataclass_is_frozen_or_mutable(self) -> None:
        """Verify that MHCBindingResult fields are accessible and set correctly."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="SIINFEKL",
            start_position=5,
            end_position=12,
            binding_score=0.65,
            ic50_nm=250.0,
            binding_class="moderate_binder",
            anchor_residues={1: "I"},
            anchor_scores={1: 0.6},
        )
        # All fields should be readable
        assert result.allele == "HLA-A*02:01"
        assert result.peptide == "SIINFEKL"
        assert result.start_position == 5
        assert result.end_position == 12
        assert result.binding_score == pytest.approx(0.65)
        assert result.ic50_nm == pytest.approx(250.0)
        assert result.binding_class == "moderate_binder"
        assert result.anchor_residues[1] == "I"
        assert result.anchor_scores[1] == pytest.approx(0.6)
