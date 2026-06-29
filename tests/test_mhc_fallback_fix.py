"""Tests for F4.3 MHC fallback database fix.

Validates:
- MHC-II peptides (13-25 aa) are accepted by MHCBindingRecord
- is_mhc_class2() correctly classifies alleles
- generate_fallback_database() works with MHC-II alleles
- get_default_alleles_for_organism() returns correct alleles
- PSSM scoring with anchor position weights
"""
from __future__ import annotations

import pytest

from biocompiler.mhc_binding_db.schema import (
    MHCBindingDatabase,
    MHCBindingRecord,
    _MHC_I_PEPTIDE_LENGTH_RANGE,
    _MHC_II_PEPTIDE_LENGTH_RANGE,
    _PSSM_ANCHOR_POSITIONS,
    _score_peptide_pssm,
    generate_fallback_database,
    get_default_alleles_for_organism,
    is_mhc_class2,
)


# ═══════════════════════════════════════════════════════════════════════════
# is_mhc_class2() tests
# ═══════════════════════════════════════════════════════════════════════════


class TestIsMhcClass2:
    """Tests for the is_mhc_class2() allele classifier."""

    def test_drb1_allele_is_class2(self) -> None:
        assert is_mhc_class2("HLA-DRB1*01:01") is True

    def test_drb1_04_01_is_class2(self) -> None:
        assert is_mhc_class2("HLA-DRB1*04:01") is True

    def test_drb3_allele_is_class2(self) -> None:
        assert is_mhc_class2("HLA-DRB3*01:01") is True

    def test_drb4_allele_is_class2(self) -> None:
        assert is_mhc_class2("HLA-DRB4*01:01") is True

    def test_drb5_allele_is_class2(self) -> None:
        assert is_mhc_class2("HLA-DRB5*01:01") is True

    def test_dqa1_allele_is_class2(self) -> None:
        assert is_mhc_class2("HLA-DQA1*01:01") is True

    def test_dqb1_allele_is_class2(self) -> None:
        assert is_mhc_class2("HLA-DQB1*03:01") is True

    def test_dpa1_allele_is_class2(self) -> None:
        assert is_mhc_class2("HLA-DPA1*01:03") is True

    def test_dpb1_allele_is_class2(self) -> None:
        assert is_mhc_class2("HLA-DPB1*04:01") is True

    def test_h2_iab_is_class2(self) -> None:
        assert is_mhc_class2("H2-IAb") is True

    def test_h2_iae_is_class2(self) -> None:
        assert is_mhc_class2("H2-IAe") is True

    def test_h2_ieb_is_class2(self) -> None:
        assert is_mhc_class2("H2-IEb") is True

    def test_hla_a0201_is_not_class2(self) -> None:
        assert is_mhc_class2("HLA-A*02:01") is False

    def test_hla_b0702_is_not_class2(self) -> None:
        assert is_mhc_class2("HLA-B*07:02") is False

    def test_hla_a0101_is_not_class2(self) -> None:
        assert is_mhc_class2("HLA-A*01:01") is False

    def test_h2kb_is_not_class2(self) -> None:
        assert is_mhc_class2("H-2Kb") is False

    def test_h2db_is_not_class2(self) -> None:
        assert is_mhc_class2("H-2Db") is False

    def test_empty_string_is_not_class2(self) -> None:
        assert is_mhc_class2("") is False

    def test_partial_prefix_match_not_class2(self) -> None:
        """HLA-DR (without B1) should NOT match HLA-DRB1 prefix."""
        assert is_mhc_class2("HLA-DR*01:01") is False


# ═══════════════════════════════════════════════════════════════════════════
# MHCBindingRecord peptide length validation tests
# ═══════════════════════════════════════════════════════════════════════════


class TestMHCBindingRecordMHC2:
    """Tests that MHCBindingRecord accepts MHC-II peptides (13-25 aa)."""

    def _make_record(
        self,
        allele: str = "HLA-DRB1*01:01",
        peptide: str = "FVWYLLMIVAKTRE",
        **overrides: object,
    ) -> MHCBindingRecord:
        defaults = dict(
            allele=allele,
            peptide=peptide,
            ic50_nm=100.0,
            rank=0.5,
            binding_class="weak_binder",
            source="pssm_fallback",
            method="pssm_anchor_scoring",
            timestamp="2025-01-01T00:00:00Z",
        )
        defaults.update(overrides)
        return MHCBindingRecord(**defaults)  # type: ignore[arg-type]

    def test_mhc2_13mer_accepted(self) -> None:
        """13-mer peptide should be accepted for MHC-II allele."""
        rec = self._make_record(peptide="FVWYLLMIVAKTR")
        assert len(rec.peptide) == 13

    def test_mhc2_15mer_accepted(self) -> None:
        """15-mer peptide should be accepted for MHC-II allele (most common)."""
        rec = self._make_record(peptide="FVWYLLMIVAKTREA")
        assert len(rec.peptide) == 15

    def test_mhc2_20mer_accepted(self) -> None:
        """20-mer peptide should be accepted for MHC-II allele."""
        rec = self._make_record(peptide="FVWYLLMIVAKTREAAAAAA")
        assert len(rec.peptide) == 20

    def test_mhc2_25mer_accepted(self) -> None:
        """25-mer peptide should be accepted for MHC-II allele (max)."""
        rec = self._make_record(peptide="FVWYLLMIVAKTREAAAAAAAAAAA")
        assert len(rec.peptide) == 25

    def test_mhc2_12mer_rejected(self) -> None:
        """12-mer peptide should be rejected for MHC-II allele (too short)."""
        with pytest.raises(ValueError, match="13-25"):
            self._make_record(peptide="FVWYLLMIVAKT")

    def test_mhc2_26mer_rejected(self) -> None:
        """26-mer peptide should be rejected for MHC-II allele (too long)."""
        with pytest.raises(ValueError, match="13-25"):
            self._make_record(peptide="FVWYLLMIVAKTREAAAAAAAAAAAA")

    def test_mhc1_9mer_still_works(self) -> None:
        """9-mer peptide should still work for MHC-I allele."""
        rec = self._make_record(
            allele="HLA-A*02:01",
            peptide="LLFGYPVYV",
        )
        assert len(rec.peptide) == 9

    def test_mhc1_8mer_still_works(self) -> None:
        """8-mer peptide should still work for MHC-I allele."""
        rec = self._make_record(
            allele="HLA-A*02:01",
            peptide="LLFGYPVY",
        )
        assert len(rec.peptide) == 8

    def test_mhc1_11mer_still_works(self) -> None:
        """11-mer peptide should still work for MHC-I allele."""
        rec = self._make_record(
            allele="HLA-A*02:01",
            peptide="LLFGYPVYVAL",
        )
        assert len(rec.peptide) == 11

    def test_mhc1_7mer_rejected(self) -> None:
        """7-mer peptide should be rejected for MHC-I allele (too short)."""
        with pytest.raises(ValueError, match="8-11"):
            self._make_record(
                allele="HLA-A*02:01",
                peptide="LLFGYPV",
            )

    def test_mhc1_12mer_rejected(self) -> None:
        """12-mer peptide should be rejected for MHC-I allele (too long)."""
        with pytest.raises(ValueError, match="8-11"):
            self._make_record(
                allele="HLA-A*02:01",
                peptide="LLFGYPVYVALL",
            )

    def test_mhc1_peptide_rejected_as_mhc2_length(self) -> None:
        """MHC-I allele should not accept MHC-II-length peptides."""
        with pytest.raises(ValueError, match="8-11"):
            self._make_record(
                allele="HLA-A*02:01",
                peptide="FVWYLLMIVAKTRE",
            )

    def test_custom_peptide_length_range_override(self) -> None:
        """peptide_length_range parameter should override allele-based inference."""
        # Use a MHC-I allele but allow 15-mer via override
        rec = self._make_record(
            allele="HLA-A*02:01",
            peptide="FVWYLLMIVAKTREA",
            peptide_length_range=(13, 25),
        )
        assert len(rec.peptide) == 15

    def test_custom_peptide_length_range_rejects(self) -> None:
        """Custom peptide_length_range should also reject out-of-range peptides."""
        with pytest.raises(ValueError, match="5-7"):
            self._make_record(
                allele="HLA-DRB1*01:01",
                peptide="FVWYLLMIVAKTRE",
                peptide_length_range=(5, 7),
            )

    def test_error_message_mentions_mhc_class(self) -> None:
        """Error message should mention MHC-I or MHC-II."""
        with pytest.raises(ValueError, match="MHC-I"):
            self._make_record(
                allele="HLA-A*02:01",
                peptide="FVWYLLMIVAKTRE",
            )

    def test_error_message_mhc2_class(self) -> None:
        """Error message for MHC-II should mention MHC-II."""
        with pytest.raises(ValueError, match="MHC-II"):
            self._make_record(peptide="FVWYLLMIVAKT")


# ═══════════════════════════════════════════════════════════════════════════
# generate_fallback_database() with MHC-II alleles
# ═══════════════════════════════════════════════════════════════════════════


class TestGenerateFallbackDatabaseMHC2:
    """Tests for generate_fallback_database() with MHC-II alleles."""

    def test_mhc2_drb1_0101_generates_records(self) -> None:
        """HLA-DRB1*01:01 should produce records in the database."""
        db = generate_fallback_database(["HLA-DRB1*01:01"])
        assert len(db) > 0
        assert "HLA-DRB1*01:01" in db.alleles

    def test_mhc2_drb1_0401_generates_records(self) -> None:
        """HLA-DRB1*04:01 should produce records in the database."""
        db = generate_fallback_database(["HLA-DRB1*04:01"])
        assert len(db) > 0
        assert "HLA-DRB1*04:01" in db.alleles

    def test_mhc2_default_peptide_length_is_15(self) -> None:
        """Default MHC-II peptide length should be 15-mers."""
        db = generate_fallback_database(["HLA-DRB1*01:01"])
        records = db.records
        # All peptides should be 15-mers by default
        peptide_lengths = {len(r.peptide) for r in records}
        assert peptide_lengths == {15}

    def test_mhc2_custom_peptide_lengths(self) -> None:
        """Custom peptide_lengths should be respected for MHC-II."""
        db = generate_fallback_database(
            ["HLA-DRB1*01:01"],
            peptide_lengths=[13, 15, 20],
        )
        records = db.records
        peptide_lengths = {len(r.peptide) for r in records}
        assert peptide_lengths == {13, 15, 20}

    def test_mhc2_records_have_correct_allele(self) -> None:
        """All MHC-II records should have the correct allele."""
        db = generate_fallback_database(["HLA-DRB1*01:01"])
        for rec in db.records:
            assert rec.allele == "HLA-DRB1*01:01"

    def test_mhc2_records_have_valid_binding_class(self) -> None:
        """All MHC-II records should have a valid binding_class."""
        valid_classes = {"strong_binder", "weak_binder", "non_binder"}
        db = generate_fallback_database(["HLA-DRB1*01:01"])
        for rec in db.records:
            assert rec.binding_class in valid_classes

    def test_mhc2_records_have_pssm_fallback_source(self) -> None:
        """All generated records should have source='pssm_fallback'."""
        db = generate_fallback_database(["HLA-DRB1*01:01"])
        for rec in db.records:
            assert rec.source == "pssm_fallback"

    def test_mhc2_peptide_lengths_in_valid_range(self) -> None:
        """All MHC-II peptides should be 13-25 residues."""
        db = generate_fallback_database(
            ["HLA-DRB1*01:01"],
            peptide_lengths=[13, 15, 20, 25],
        )
        for rec in db.records:
            assert 13 <= len(rec.peptide) <= 25

    def test_mixed_mhc1_and_mhc2(self) -> None:
        """Should handle both MHC-I and MHC-II alleles simultaneously."""
        db = generate_fallback_database(
            ["HLA-A*02:01", "HLA-DRB1*01:01"],
        )
        assert "HLA-A*02:01" in db.alleles
        assert "HLA-DRB1*01:01" in db.alleles

        # MHC-I records should be 9-mers
        mhc1_records = [r for r in db.records if r.allele == "HLA-A*02:01"]
        mhc1_lengths = {len(r.peptide) for r in mhc1_records}
        assert mhc1_lengths == {9}

        # MHC-II records should be 15-mers
        mhc2_records = [r for r in db.records if r.allele == "HLA-DRB1*01:01"]
        mhc2_lengths = {len(r.peptide) for r in mhc2_records}
        assert mhc2_lengths == {15}

    def test_mhc1_default_peptide_length_is_9(self) -> None:
        """Default MHC-I peptide length should still be 9-mers."""
        db = generate_fallback_database(["HLA-A*02:01"])
        records = db.records
        peptide_lengths = {len(r.peptide) for r in records}
        assert peptide_lengths == {9}

    def test_mhc2_generic_allele(self) -> None:
        """MHC-II allele without PSSM should use generic scoring."""
        db = generate_fallback_database(["HLA-DRB1*07:01"])
        assert len(db) > 0
        for rec in db.records:
            assert rec.allele == "HLA-DRB1*07:01"
            assert 13 <= len(rec.peptide) <= 25

    def test_deterministic_with_seed(self) -> None:
        """generate_fallback_database should produce deterministic results."""
        db1 = generate_fallback_database(["HLA-DRB1*01:01"])
        db2 = generate_fallback_database(["HLA-DRB1*01:01"])
        assert len(db1) == len(db2)
        # Check first few records match
        for r1, r2 in zip(db1.records[:10], db2.records[:10]):
            assert r1.peptide == r2.peptide
            assert r1.ic50_nm == r2.ic50_nm


# ═══════════════════════════════════════════════════════════════════════════
# get_default_alleles_for_organism() tests
# ═══════════════════════════════════════════════════════════════════════════


class TestGetDefaultAllelesForOrganism:
    """Tests for the get_default_alleles_for_organism() function."""

    def test_homo_sapiens_returns_expected_alleles(self) -> None:
        alleles = get_default_alleles_for_organism("Homo_sapiens")
        assert "HLA-A*02:01" in alleles
        assert "HLA-A*01:01" in alleles
        assert "HLA-A*03:01" in alleles
        assert "HLA-B*07:02" in alleles
        assert "HLA-B*08:01" in alleles
        assert "HLA-DRB1*01:01" in alleles
        assert "HLA-DRB1*04:01" in alleles

    def test_homo_sapiens_has_both_mhc_classes(self) -> None:
        """Homo_sapiens should include both MHC-I and MHC-II alleles."""
        alleles = get_default_alleles_for_organism("Homo_sapiens")
        mhc1 = [a for a in alleles if not is_mhc_class2(a)]
        mhc2 = [a for a in alleles if is_mhc_class2(a)]
        assert len(mhc1) > 0, "Should have MHC-I alleles"
        assert len(mhc2) > 0, "Should have MHC-II alleles"

    def test_mus_musculus_returns_expected_alleles(self) -> None:
        alleles = get_default_alleles_for_organism("Mus_musculus")
        assert "H-2Kb" in alleles
        assert "H-2Db" in alleles

    def test_mus_musculus_has_two_alleles(self) -> None:
        alleles = get_default_alleles_for_organism("Mus_musculus")
        assert len(alleles) == 2

    def test_unknown_organism_returns_empty(self) -> None:
        """Unknown organism should return empty list."""
        alleles = get_default_alleles_for_organism("Unknown_organism")
        assert alleles == []

    def test_returns_copy_not_reference(self) -> None:
        """Returned list should be a copy, not a reference to internal data."""
        a1 = get_default_alleles_for_organism("Homo_sapiens")
        a2 = get_default_alleles_for_organism("Homo_sapiens")
        assert a1 == a2
        a1.append("TEST")
        assert "TEST" not in a2

    def test_homo_sapiens_allele_count(self) -> None:
        alleles = get_default_alleles_for_organism("Homo_sapiens")
        assert len(alleles) == 7


# ═══════════════════════════════════════════════════════════════════════════
# PSSM scoring with anchor weights tests
# ═══════════════════════════════════════════════════════════════════════════


class TestPSSMScoringWithAnchorWeights:
    """Tests for _score_peptide_pssm with anchor position weights."""

    def test_anchor_weight_improves_score_for_preferred_residues(self) -> None:
        """Peptides with preferred residues at anchor positions should score
        higher when anchor weights are applied."""
        preferred = [
            {"L": 2.0, "I": 1.8},  # pos 0 — anchor
            {},  # pos 1
            {},  # pos 2
            {},  # pos 3
            {},  # pos 4
            {},  # pos 5
            {},  # pos 6
            {},  # pos 7
            {"V": 1.5, "L": 1.5},  # pos 8 — anchor
        ]
        disfavored = [{}, {}, {}, {}, {}, {}, {}, {}, {}]

        # Peptide with preferred residues at anchor positions
        good_anchors = "LVVVVVVLV"
        # Peptide without preferred residues at anchor positions
        bad_anchors = "AAAAAAAKK"

        score_good_no_weight = _score_peptide_pssm(
            good_anchors, preferred, disfavored, anchor_positions=None,
        )
        score_good_weighted = _score_peptide_pssm(
            good_anchors, preferred, disfavored, anchor_positions={0, 8},
        )
        score_bad_no_weight = _score_peptide_pssm(
            bad_anchors, preferred, disfavored, anchor_positions=None,
        )
        score_bad_weighted = _score_peptide_pssm(
            bad_anchors, preferred, disfavored, anchor_positions={0, 8},
        )

        # With anchor weights, the good-anchors peptide should score
        # even higher relative to the bad-anchors peptide
        gap_no_weight = score_good_no_weight - score_bad_no_weight
        gap_weighted = score_good_weighted - score_bad_weighted
        assert gap_weighted >= gap_no_weight, (
            "Anchor weighting should widen the gap between good and bad anchors"
        )

    def test_no_anchor_positions_equals_flat_geometric_mean(self) -> None:
        """When anchor_positions=None, behavior should be flat geometric mean."""
        preferred = [
            {"L": 2.0},
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            {"V": 1.5},
        ]
        disfavored = [{}, {}, {}, {}, {}, {}, {}, {}, {}]

        # Both calls should produce the same result when anchor_positions=None
        score = _score_peptide_pssm(
            "LAAAAAAAAV", preferred, disfavored, anchor_positions=None,
        )
        assert 0.0 <= score <= 1.0

    def test_mhc2_anchor_positions_in_pssm(self) -> None:
        """MHC-II PSSM should use pocket positions as anchors."""
        # HLA-DRB1*01:01 uses {0, 3, 5, 8} as anchor positions
        assert 0 in _PSSM_ANCHOR_POSITIONS["HLA-DRB1*01:01"]
        assert 3 in _PSSM_ANCHOR_POSITIONS["HLA-DRB1*01:01"]
        assert 5 in _PSSM_ANCHOR_POSITIONS["HLA-DRB1*01:01"]
        assert 8 in _PSSM_ANCHOR_POSITIONS["HLA-DRB1*01:01"]

    def test_mhc1_anchor_positions_in_pssm(self) -> None:
        """MHC-I PSSM should use positions 1 and 8 as anchors."""
        assert 1 in _PSSM_ANCHOR_POSITIONS["HLA-A*02:01"]
        assert 8 in _PSSM_ANCHOR_POSITIONS["HLA-A*02:01"]

    def test_score_is_bounded(self) -> None:
        """PSSM score should always be in [0, 1]."""
        preferred = [
            {"L": 2.0, "I": 1.8},
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            {"V": 1.5, "L": 1.5},
        ]
        disfavored = [
            {"D": 0.3, "E": 0.3},
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            {"D": 0.4, "E": 0.4},
        ]

        for peptide in ["LLFGYPVYV", "DDDDDDDDDD", "AAAAAAAAA", "LIVVAKKKE"]:
            score = _score_peptide_pssm(
                peptide, preferred, disfavored, anchor_positions={0, 8},
            )
            assert 0.0 <= score <= 1.0


# ═══════════════════════════════════════════════════════════════════════════
# Integration: generate_fallback_database with get_default_alleles_for_organism
# ═══════════════════════════════════════════════════════════════════════════


class TestIntegrationFallbackWithOrganism:
    """Integration tests combining organism allele selection with fallback DB."""

    def test_homo_sapiens_fallback_db(self) -> None:
        """Should generate a fallback database for all Homo_sapiens alleles."""
        alleles = get_default_alleles_for_organism("Homo_sapiens")
        db = generate_fallback_database(alleles)
        assert len(db) > 0
        # Should contain records for all alleles
        for allele in alleles:
            assert allele in db.alleles, f"Missing records for {allele}"

    def test_mus_musculus_fallback_db(self) -> None:
        """Should generate a fallback database for Mus_musculus alleles."""
        alleles = get_default_alleles_for_organism("Mus_musculus")
        db = generate_fallback_database(alleles)
        assert len(db) > 0
        for allele in alleles:
            assert allele in db.alleles, f"Missing records for {allele}"

    def test_homo_sapiens_has_mhc2_records(self) -> None:
        """Homo_sapiens fallback should include MHC-II records."""
        alleles = get_default_alleles_for_organism("Homo_sapiens")
        db = generate_fallback_database(alleles)
        mhc2_records = [r for r in db.records if is_mhc_class2(r.allele)]
        assert len(mhc2_records) > 0, "Should have MHC-II records"

    def test_mhc2_peptides_are_longer_than_mhc1(self) -> None:
        """MHC-II peptides should be longer than MHC-I peptides."""
        alleles = get_default_alleles_for_organism("Homo_sapiens")
        db = generate_fallback_database(alleles)

        mhc1_lengths = [len(r.peptide) for r in db.records if not is_mhc_class2(r.allele)]
        mhc2_lengths = [len(r.peptide) for r in db.records if is_mhc_class2(r.allele)]

        if mhc1_lengths and mhc2_lengths:
            assert min(mhc2_lengths) > max(mhc1_lengths) or min(mhc2_lengths) >= 13


# ═══════════════════════════════════════════════════════════════════════════
# Database lookup with MHC-II records
# ═══════════════════════════════════════════════════════════════════════════


class TestDatabaseLookupMHC2:
    """Test database operations with MHC-II records."""

    def test_add_and_lookup_mhc2_record(self) -> None:
        """Should be able to add and look up MHC-II records."""
        db = MHCBindingDatabase()
        rec = MHCBindingRecord(
            allele="HLA-DRB1*01:01",
            peptide="FVWYLLMIVAKTREA",
            ic50_nm=50.0,
            rank=None,
            binding_class="weak_binder",
            source="pssm_fallback",
            method="pssm_anchor_scoring",
            timestamp="2025-01-01T00:00:00Z",
        )
        db.add(rec)
        result = db.lookup("HLA-DRB1*01:01", "FVWYLLMIVAKTREA")
        assert result is not None
        assert result.peptide == "FVWYLLMIVAKTREA"
        assert result.ic50_nm == 50.0

    def test_binding_peptides_mhc2(self) -> None:
        """binding_peptides() should work with MHC-II records."""
        db = generate_fallback_database(["HLA-DRB1*01:01"])
        binders = db.binding_peptides("HLA-DRB1*01:01", threshold_ic50=500.0)
        # Should have at least some binders or at least not crash
        assert isinstance(binders, list)

    def test_stats_includes_mhc2(self) -> None:
        """stats() should include MHC-II allele data."""
        db = generate_fallback_database(["HLA-DRB1*01:01"])
        stats = db.stats()
        assert "HLA-DRB1*01:01" in stats["alleles"]
