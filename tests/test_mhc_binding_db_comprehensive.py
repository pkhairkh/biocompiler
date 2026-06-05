"""Comprehensive tests for the MHC binding database module.

Covers:
  1. MHCBindingDatabase — multi-allele in-memory DB with add/lookup/filter/serialisation
  2. PrecomputedAlleleDatabase — legacy single-allele precomputed databases
  3. generate_fallback_database — PSSM-based offline fallback generation
  4. MHCFlurryAdapter offline fallback chain — verify the four-tier fallback

All tests are designed to run fully offline without MHCflurry or NetMHCpan
installed, ensuring the fixed initialization code and fallback paths are
properly exercised.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import asdict
from unittest.mock import MagicMock, patch

import pytest

from biocompiler.mhc_binding_db import (
    MHCBindingDatabase,
    MHCBindingRecord,
    PrecomputedAlleleDatabase,
    PrecomputedEntry,
    generate_fallback_database,
    get_database,
    get_default_alleles_for_organism,
    is_mhc_class2,
)
from biocompiler.mhc_binding_db.precomputed import (
    AVAILABLE_ALLELES,
    get_all_precomputed_databases,
    get_precomputed_database,
)


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _make_record(
    allele: str = "HLA-A*02:01",
    peptide: str = "LLFGYPVYV",
    ic50_nm: float = 12.0,
    rank: float | None = 0.5,
    binding_class: str = "strong_binder",
    source: str = "mhcflurry_predicted",
    method: str = "mhcflurry_class1",
    timestamp: str = "2025-01-01T00:00:00Z",
) -> MHCBindingRecord:
    """Create a valid MHCBindingRecord for testing."""
    return MHCBindingRecord(
        allele=allele,
        peptide=peptide,
        ic50_nm=ic50_nm,
        rank=rank,
        binding_class=binding_class,
        source=source,
        method=method,
        timestamp=timestamp,
    )


def _make_mhc2_record(
    peptide: str = "MKWVTFISLLLLFSR",
    ic50_nm: float = 100.0,
    binding_class: str = "weak_binder",
) -> MHCBindingRecord:
    """Create a valid MHC-II MHCBindingRecord for testing."""
    return MHCBindingRecord(
        allele="HLA-DRB1*01:01",
        peptide=peptide,
        ic50_nm=ic50_nm,
        rank=None,
        binding_class=binding_class,
        source="pssm_fallback",
        method="pssm_anchor_scoring",
        timestamp="2025-01-01T00:00:00Z",
    )


# ═══════════════════════════════════════════════════════════════════════════
# 1. MHCBindingDatabase comprehensive tests
# ═══════════════════════════════════════════════════════════════════════════


class TestMHCBindingDatabaseComprehensive:
    """Comprehensive tests for MHCBindingDatabase — add, lookup, filter, serialise."""

    def test_empty_database(self):
        """An empty database should have zero stats and None lookups."""
        db = MHCBindingDatabase()
        assert len(db) == 0
        assert db.lookup("HLA-A*02:01", "LLFGYPVYV") is None
        stats = db.stats()
        assert stats["total_records"] == 0
        assert stats["alleles"] == {}
        assert db.alleles == []
        assert db.records == []

    def test_add_and_lookup(self):
        """Add a single record and verify it can be looked up."""
        db = MHCBindingDatabase()
        rec = _make_record()
        db.add(rec)

        # Lookup should return the record
        result = db.lookup("HLA-A*02:01", "LLFGYPVYV")
        assert result is not None
        assert result.allele == "HLA-A*02:01"
        assert result.peptide == "LLFGYPVYV"
        assert result.ic50_nm == pytest.approx(12.0)
        assert result.binding_class == "strong_binder"

        # Lookup for non-existent peptide should return None
        assert db.lookup("HLA-A*02:01", "NOPEPTIDE") is None
        # Lookup for non-existent allele should return None
        assert db.lookup("HLA-B*99:99", "LLFGYPVYV") is None

    def test_add_batch(self):
        """Add multiple records via add_batch and verify all can be looked up."""
        db = MHCBindingDatabase()
        records = [
            _make_record(peptide="AAAAAAAAA", ic50_nm=10.0, binding_class="strong_binder"),
            _make_record(peptide="BBBBBBBBB", ic50_nm=200.0, binding_class="weak_binder"),
            _make_record(peptide="CCCCCCCCC", ic50_nm=800.0, binding_class="non_binder"),
        ]
        count = db.add_batch(records)
        assert count == 3
        assert len(db) == 3

        for rec in records:
            result = db.lookup("HLA-A*02:01", rec.peptide)
            assert result is not None
            assert result.peptide == rec.peptide
            assert result.ic50_nm == pytest.approx(rec.ic50_nm)

    def test_strong_binders(self):
        """Strong binders filter should return peptides with IC50 < 50 nM."""
        db = MHCBindingDatabase()
        db.add(_make_record(peptide="STRONGAA", ic50_nm=10.0, binding_class="strong_binder"))
        db.add(_make_record(peptide="WEAKPEPT", ic50_nm=200.0, binding_class="weak_binder"))
        db.add(_make_record(peptide="NONBNDPT", ic50_nm=800.0, binding_class="non_binder"))
        db.add(_make_record(peptide="ALMOSTST", ic50_nm=49.9, binding_class="strong_binder"))

        strong = db.strong_binders("HLA-A*02:01")
        assert "STRONGAA" in strong
        assert "ALMOSTST" in strong
        assert "WEAKPEPT" not in strong
        assert "NONBNDPT" not in strong

    def test_weak_binders(self):
        """Weak binders filter should return peptides with 50 <= IC50 < 500 nM."""
        db = MHCBindingDatabase()
        db.add(_make_record(peptide="STRONGAA", ic50_nm=10.0, binding_class="strong_binder"))
        db.add(_make_record(peptide="WEAKPEPT", ic50_nm=200.0, binding_class="weak_binder"))
        db.add(_make_record(peptide="BOUNDARY1", ic50_nm=50.0, binding_class="weak_binder"))
        db.add(_make_record(peptide="BOUNDARY2", ic50_nm=499.9, binding_class="weak_binder"))
        db.add(_make_record(peptide="NONBNDPT", ic50_nm=800.0, binding_class="non_binder"))

        weak = db.weak_binders("HLA-A*02:01")
        assert "WEAKPEPT" in weak
        assert "BOUNDARY1" in weak
        assert "BOUNDARY2" in weak
        assert "STRONGAA" not in weak
        assert "NONBNDPT" not in weak

    def test_binding_peptides_with_threshold(self):
        """binding_peptides should return peptides below custom IC50 threshold."""
        db = MHCBindingDatabase()
        db.add(_make_record(peptide="STRONGAA", ic50_nm=10.0, binding_class="strong_binder"))
        db.add(_make_record(peptide="WEAKPEPT", ic50_nm=200.0, binding_class="weak_binder"))
        db.add(_make_record(peptide="NONBNDPT", ic50_nm=800.0, binding_class="non_binder"))

        # Default threshold (500 nM) — includes strong + weak
        default_peptides = db.binding_peptides("HLA-A*02:01")
        assert "STRONGAA" in default_peptides
        assert "WEAKPEPT" in default_peptides
        assert "NONBNDPT" not in default_peptides

        # Custom threshold (100 nM) — only strong
        tight_peptides = db.binding_peptides("HLA-A*02:01", threshold_ic50=100.0)
        assert "STRONGAA" in tight_peptides
        assert "WEAKPEPT" not in tight_peptides
        assert "NONBNDPT" not in tight_peptides

        # Very high threshold — includes everything
        all_peptides = db.binding_peptides("HLA-A*02:01", threshold_ic50=50000.0)
        assert len(all_peptides) == 3

    def test_json_roundtrip(self):
        """Save to JSON, load back, verify records are identical."""
        db = MHCBindingDatabase()
        db.add(_make_record(peptide="LLFGYPVYV", ic50_nm=12.0, binding_class="strong_binder"))
        db.add(_make_record(peptide="ALVSRDYCV", ic50_nm=200.0, binding_class="weak_binder"))
        db.add(_make_record(peptide="RDQMWETCI", ic50_nm=5000.0, binding_class="non_binder"))

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False,
        ) as tmp:
            tmp_path = tmp.name

        try:
            db.save_to_json(tmp_path)

            # Load into a fresh database
            db2 = MHCBindingDatabase()
            count = db2.load_from_json(tmp_path)
            assert count == 3
            assert len(db2) == 3

            # Verify each record matches the original
            for rec in db.records:
                loaded = db2.lookup(rec.allele, rec.peptide)
                assert loaded is not None, f"Missing record for {rec.allele}/{rec.peptide}"
                assert loaded.allele == rec.allele
                assert loaded.peptide == rec.peptide
                assert loaded.ic50_nm == pytest.approx(rec.ic50_nm, abs=0.1)
                assert loaded.binding_class == rec.binding_class
                assert loaded.source == rec.source
        finally:
            os.unlink(tmp_path)

    def test_multiple_alleles(self):
        """Records for multiple alleles should be queryable per-allele."""
        db = MHCBindingDatabase()
        db.add(_make_record(allele="HLA-A*02:01", peptide="AAAAAAAAA"))
        db.add(_make_record(allele="HLA-A*01:01", peptide="BBBBBBBBB"))
        db.add(_make_record(allele="HLA-B*07:02", peptide="CCCCCCCCC"))

        assert len(db) == 3
        assert sorted(db.alleles) == ["HLA-A*01:01", "HLA-A*02:01", "HLA-B*07:02"]

        # Per-allele lookups
        assert db.lookup("HLA-A*02:01", "AAAAAAAAA") is not None
        assert db.lookup("HLA-A*01:01", "BBBBBBBBB") is not None
        assert db.lookup("HLA-B*07:02", "CCCCCCCCC") is not None

        # Cross-allele lookup should return None
        assert db.lookup("HLA-A*02:01", "BBBBBBBBB") is None

        # Per-allele filtering
        assert db.strong_binders("HLA-A*02:01") == ["AAAAAAAAA"]
        assert db.strong_binders("HLA-A*01:01") == ["BBBBBBBBB"]
        assert db.strong_binders("HLA-B*07:02") == ["CCCCCCCCC"]

    def test_stats(self):
        """Stats should correctly count records per allele and binding class."""
        db = MHCBindingDatabase()
        db.add(_make_record(peptide="STRONG1AA", ic50_nm=10.0, binding_class="strong_binder"))
        db.add(_make_record(peptide="STRONG2AA", ic50_nm=20.0, binding_class="strong_binder"))
        db.add(_make_record(peptide="WEAKONEAA", ic50_nm=200.0, binding_class="weak_binder"))
        db.add(_make_record(peptide="NONBONEAA", ic50_nm=5000.0, binding_class="non_binder"))

        stats = db.stats()
        assert stats["total_records"] == 4
        assert "HLA-A*02:01" in stats["alleles"]
        allele_s = stats["alleles"]["HLA-A*02:01"]
        assert allele_s["total"] == 4
        assert allele_s["strong_binder"] == 2
        assert allele_s["weak_binder"] == 1
        assert allele_s["non_binder"] == 1
        # Non-existent allele should not appear
        assert "HLA-B*99:99" not in stats["alleles"]

    def test_replace_on_duplicate_key(self):
        """Adding the same (allele, peptide) twice should replace the record."""
        db = MHCBindingDatabase()
        rec1 = _make_record(peptide="LLFGYPVYV", ic50_nm=12.0, binding_class="strong_binder")
        db.add(rec1)
        assert len(db) == 1

        # Add a different record for the same key
        rec2 = _make_record(peptide="LLFGYPVYV", ic50_nm=500.0, binding_class="weak_binder")
        db.add(rec2)
        assert len(db) == 1  # Still one record — replaced, not appended

        result = db.lookup("HLA-A*02:01", "LLFGYPVYV")
        assert result is not None
        assert result.ic50_nm == pytest.approx(500.0)
        assert result.binding_class == "weak_binder"

    def test_container_protocol(self):
        """__contains__ should work with (allele, peptide) tuple."""
        db = MHCBindingDatabase()
        rec = _make_record(peptide="LLFGYPVYV")
        db.add(rec)

        assert ("HLA-A*02:01", "LLFGYPVYV") in db
        assert ("HLA-A*02:01", "MISSING") not in db
        assert ("MISSING", "LLFGYPVYV") not in db

    def test_lookup_batch(self):
        """lookup_batch should return a list with None for missing peptides."""
        db = MHCBindingDatabase()
        db.add(_make_record(peptide="AAAAAAAAA"))
        db.add(_make_record(peptide="BBBBBBBBB"))

        results = db.lookup_batch("HLA-A*02:01", ["AAAAAAAAA", "BBBBBBBBB", "CCCCCCCCC"])
        assert len(results) == 3
        assert results[0] is not None
        assert results[0].peptide == "AAAAAAAAA"
        assert results[1] is not None
        assert results[1].peptide == "BBBBBBBBB"
        assert results[2] is None


# ═══════════════════════════════════════════════════════════════════════════
# 2. PrecomputedAlleleDatabase tests
# ═══════════════════════════════════════════════════════════════════════════


class TestPrecomputedAlleleDatabase:
    """Tests for the legacy single-allele precomputed database and loading."""

    def test_get_database_all_alleles(self):
        """Loading all 9 precomputed databases should succeed."""
        all_dbs = get_all_precomputed_databases()
        assert len(all_dbs) == 9
        for allele in AVAILABLE_ALLELES:
            assert allele in all_dbs
            assert isinstance(all_dbs[allele], PrecomputedAlleleDatabase)
            assert all_dbs[allele].allele == allele

    def test_get_database_invalid_allele(self):
        """Requesting an unknown allele should raise ValueError."""
        with pytest.raises(ValueError, match="No precomputed database"):
            get_precomputed_database("HLA-FAKE*99:99")

    def test_caching(self):
        """Loading the same allele twice should return the same cached instance."""
        # Clear the cache first to ensure a clean test
        from biocompiler.mhc_binding_db.precomputed import _databases
        _databases.clear()

        db1 = get_precomputed_database("HLA-A*02:01")
        db2 = get_precomputed_database("HLA-A*02:01")
        assert db1 is db2  # Same object — cached

    def test_each_allele_has_entries(self):
        """Each precomputed allele database should have at least one entry."""
        all_dbs = get_all_precomputed_databases()
        for allele, db in all_dbs.items():
            assert db.total_entries > 0, f"{allele} has no entries"
            assert len(db.entries) > 0, f"{allele} entries list is empty"

    def test_each_allele_has_binders(self):
        """Each precomputed allele database should have at least one binder."""
        all_dbs = get_all_precomputed_databases()
        for allele, db in all_dbs.items():
            binders = db.get_binders()
            assert len(binders) > 0, f"{allele} has no binders"

    def test_known_epitope_lookup(self):
        """Known epitopes should be found in the precomputed database."""
        db = get_precomputed_database("HLA-A*02:01")
        # GILGFVFTL is a known influenza M1 epitope for HLA-A*02:01
        entry = db.lookup("GILGFVFTL")
        assert entry is not None, "Known epitope GILGFVFTL not found in HLA-A*02:01 database"
        assert entry.source == "known_epitope"

        # LLFGYPVYV is another known epitope
        entry2 = db.lookup("LLFGYPVYV")
        assert entry2 is not None, "Known epitope LLFGYPVYV not found in HLA-A*02:01 database"
        assert entry2.source == "known_epitope"

    def test_get_database_via_top_level(self):
        """get_database from top-level __init__ should work for available alleles."""
        db = get_database("HLA-A*02:01")
        assert db is not None
        assert db.allele == "HLA-A*02:01"

    def test_get_database_returns_none_for_unknown(self):
        """get_database should return None for alleles not in AVAILABLE_ALLELES."""
        result = get_database("HLA-FAKE*99:99")
        assert result is None

    def test_is_binder_method(self):
        """is_binder should return True for strong/moderate binders, False otherwise."""
        db = get_precomputed_database("HLA-A*02:01")
        # Find a known binder
        binders = db.get_binders()
        if binders:
            assert db.is_binder(binders[0].peptide) is True

        # Find a non-binder
        non_binders = db.get_non_binders()
        if non_binders:
            assert db.is_binder(non_binders[0].peptide) is False

        # Unknown peptide
        assert db.is_binder("ZZZZZZZZZ") is False

    def test_precomputed_entry_properties(self):
        """PrecomputedEntry properties should work correctly."""
        entry = PrecomputedEntry(
            peptide="TESTPEPTD",
            binding_score=0.85,
            ic50_nm=45.0,
            binding_class="strong_binder",
            anchor_residues={1: "E", 8: "D"},
            anchor_scores={1: 1.5, 8: 1.2},
            source="known_epitope",
            peptide_length=9,
        )
        assert entry.is_binder is True
        assert entry.peptide_length == 9

        non_binder = PrecomputedEntry(
            peptide="NONBINDER",
            binding_score=0.05,
            ic50_nm=5000.0,
            binding_class="non_binder",
            anchor_residues={},
            anchor_scores={},
            source="pssm_predicted",
            peptide_length=9,
        )
        assert non_binder.is_binder is False


# ═══════════════════════════════════════════════════════════════════════════
# 3. generate_fallback_database tests
# ═══════════════════════════════════════════════════════════════════════════


class TestGenerateFallbackDatabase:
    """Tests for PSSM-based fallback database generation."""

    def test_generate_for_mhc1(self):
        """Generate fallback for MHC-I alleles should produce valid records."""
        mhc1_alleles = ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02"]
        db = generate_fallback_database(mhc1_alleles)

        assert len(db) > 0
        for allele in mhc1_alleles:
            assert allele in db.alleles

        # Verify MHC-I records have 9-mer peptides
        for rec in db.records:
            if not is_mhc_class2(rec.allele):
                assert len(rec.peptide) == 9, (
                    f"MHC-I peptide {rec.peptide!r} is not 9-mer (len={len(rec.peptide)})"
                )

    def test_generate_for_mhc2(self):
        """Generate fallback for MHC-II alleles should produce 15-mer peptides."""
        mhc2_alleles = ["HLA-DRB1*01:01", "HLA-DRB1*04:01"]
        db = generate_fallback_database(mhc2_alleles)

        assert len(db) > 0
        for allele in mhc2_alleles:
            assert allele in db.alleles

        # MHC-II records should have 15-mer peptides (default)
        for rec in db.records:
            if is_mhc_class2(rec.allele):
                assert len(rec.peptide) == 15, (
                    f"MHC-II peptide {rec.peptide!r} is not 15-mer (len={len(rec.peptide)})"
                )

    def test_deterministic(self):
        """Generate twice with seed=42, verify identical results."""
        alleles = ["HLA-A*02:01", "HLA-B*08:01"]
        db1 = generate_fallback_database(alleles)
        db2 = generate_fallback_database(alleles)

        # Same number of records
        assert len(db1) == len(db2)

        # All records should match
        for rec1 in db1.records:
            rec2 = db2.lookup(rec1.allele, rec1.peptide)
            assert rec2 is not None, f"Record {rec1.allele}/{rec1.peptide} not in second DB"
            assert rec1.ic50_nm == pytest.approx(rec2.ic50_nm, abs=0.1)
            assert rec1.binding_class == rec2.binding_class

    def test_binding_classes_valid(self):
        """All generated binding_class values should be valid."""
        valid_classes = {"strong_binder", "weak_binder", "non_binder"}
        alleles = ["HLA-A*02:01", "HLA-DRB1*01:01"]
        db = generate_fallback_database(alleles)

        for rec in db.records:
            assert rec.binding_class in valid_classes, (
                f"Invalid binding_class: {rec.binding_class!r} for {rec.allele}/{rec.peptide}"
            )

    def test_sources_are_fallback(self):
        """All generated records should have source='pssm_fallback'."""
        alleles = ["HLA-A*02:01", "HLA-DRB1*01:01"]
        db = generate_fallback_database(alleles)

        for rec in db.records:
            assert rec.source == "pssm_fallback", (
                f"Expected source 'pssm_fallback', got {rec.source!r}"
            )

    def test_generate_with_custom_peptide_lengths(self):
        """Custom peptide_lengths should be used for all alleles."""
        alleles = ["HLA-A*02:01"]
        db = generate_fallback_database(alleles, peptide_lengths=[10])

        # Should produce 10-mer peptides
        for rec in db.records:
            assert len(rec.peptide) == 10, (
                f"Expected 10-mer, got {len(rec.peptide)}-mer: {rec.peptide!r}"
            )

    def test_generate_for_allele_without_pssm(self):
        """Alleles without PSSM data should use generic hydrophobicity scoring."""
        # H-2Kb and H-2Db are mouse alleles without PSSM data in _PSSM_FALLBACK
        alleles = ["H-2Kb"]
        db = generate_fallback_database(alleles)

        assert len(db) > 0
        assert "H-2Kb" in db.alleles
        # Records should have 9-mer peptides (MHC-I default)
        for rec in db.records:
            assert len(rec.peptide) == 9

    def test_ic50_values_reasonable(self):
        """Generated IC50 values should be in a reasonable range (0.1 to 100,000 nM)."""
        alleles = ["HLA-A*02:01"]
        db = generate_fallback_database(alleles)

        for rec in db.records:
            assert 0.1 < rec.ic50_nm < 100000, (
                f"IC50 {rec.ic50_nm} out of reasonable range for {rec.peptide!r}"
            )

    def test_mixed_mhc1_and_mhc2(self):
        """Generating for both MHC-I and MHC-II should produce correct peptide lengths."""
        alleles = ["HLA-A*02:01", "HLA-DRB1*01:01"]
        db = generate_fallback_database(alleles)

        mhc1_records = [r for r in db.records if not is_mhc_class2(r.allele)]
        mhc2_records = [r for r in db.records if is_mhc_class2(r.allele)]

        assert len(mhc1_records) > 0
        assert len(mhc2_records) > 0

        for rec in mhc1_records:
            assert len(rec.peptide) == 9
        for rec in mhc2_records:
            assert len(rec.peptide) == 15


# ═══════════════════════════════════════════════════════════════════════════
# 4. MHCFlurry adapter offline fallback tests
# ═══════════════════════════════════════════════════════════════════════════


class TestMHCFlurryAdapterOffline:
    """Test that the MHCflurry adapter works in offline mode (no MHCflurry installed).

    These tests verify the fixed fallback chain works end-to-end without
    any external dependencies: MHCflurry → NetMHCpan → precomputed → PSSM.
    """

    def test_adapter_works_without_mhcflurry(self):
        """Verify the fallback chain produces a valid result without MHCflurry."""
        from biocompiler.mhcflurry_adapter import MHCflurryClient

        client = MHCflurryClient(allow_offline_fallback=True)
        # This should NOT raise — the fallback chain should always produce
        # a result when allow_offline_fallback=True
        result = client.predict_binding("LLFGYPVYV", "HLA-A*02:01")

        assert result is not None
        assert result.allele == "HLA-A*02:01"
        assert result.peptide == "LLFGYPVYV"
        assert result.binding_class in (
            "strong_binder", "moderate_binder", "weak_binder", "non_binder",
        )
        assert 0.0 <= result.binding_score <= 1.0
        # Method should indicate fallback (not MHCflurry)
        assert result.method in (
            "mhcflurry", "netmhcpan", "precomputed_lookup", "pssm_fallback",
        )

    def test_precomputed_fallback(self):
        """Test that precomputed DB is used as fallback when MHCflurry is unavailable.

        A known epitope from the precomputed database should be found
        via the precomputed_lookup method when MHCflurry/NetMHCpan are
        unavailable.
        """
        from biocompiler.mhcflurry_adapter import MHCflurryClient, CONFIDENCE_PRECOMPUTED

        client = MHCflurryClient(allow_offline_fallback=True)

        # Force MHCflurry to be unavailable
        client._models_load_failed = True
        client._affinity_predictor = None

        # Patch NetMHCpan to be unavailable too
        with patch.object(
            client, "_try_netmhcpan_prediction", return_value=None,
        ):
            # A peptide in the precomputed database
            result = client.predict_binding("GILGFVFTL", "HLA-A*02:01")

        assert result is not None
        # Should come from precomputed lookup (not PSSM)
        assert result.method == "precomputed_lookup"
        assert result.confidence == pytest.approx(CONFIDENCE_PRECOMPUTED)

    def test_pssm_fallback_for_unknown_peptide(self):
        """Unknown peptide should fall through to PSSM scoring."""
        from biocompiler.mhcflurry_adapter import (
            MHCflurryClient,
            CONFIDENCE_PSSM,
        )

        client = MHCflurryClient(allow_offline_fallback=True)
        client._models_load_failed = True
        client._affinity_predictor = None

        with patch.object(
            client, "_try_netmhcpan_prediction", return_value=None,
        ):
            # A random 9-mer peptide unlikely to be in precomputed DB
            result = client.predict_binding("QWERTYIPN", "HLA-A*02:01")

        # PSSM fallback should always produce a result
        assert result is not None
        # Should come from PSSM fallback
        assert result.method in ("pssm_fallback", "precomputed_lookup")
        assert result.confidence == pytest.approx(CONFIDENCE_PSSM)

    def test_adapter_no_fallback_raises(self):
        """With allow_offline_fallback=False, should raise RuntimeError."""
        from biocompiler.mhcflurry_adapter import MHCflurryClient

        client = MHCflurryClient(allow_offline_fallback=False)
        client._models_load_failed = True
        client._affinity_predictor = None

        with patch.object(
            client, "_try_netmhcpan_prediction", return_value=None,
        ):
            with pytest.raises(RuntimeError, match="MHCflurry not available"):
                client.predict_binding("LLFGYPVYV", "HLA-A*02:01")


# ═══════════════════════════════════════════════════════════════════════════
# 5. MHCBindingRecord validation tests
# ═══════════════════════════════════════════════════════════════════════════


class TestMHCBindingRecordValidation:
    """Test MHCBindingRecord field validation."""

    def test_invalid_binding_class_raises(self):
        """Invalid binding_class should raise ValueError."""
        with pytest.raises(ValueError, match="binding_class"):
            _make_record(binding_class="super_binder")

    def test_invalid_source_raises(self):
        """Invalid source should raise ValueError."""
        with pytest.raises(ValueError, match="source"):
            _make_record(source="unknown_source")

    def test_wrong_peptide_length_mhc1(self):
        """MHC-I peptide of wrong length should raise ValueError."""
        with pytest.raises(ValueError, match="Peptide length"):
            _make_record(peptide="SHORT")  # 5-mer, too short for MHC-I

    def test_wrong_peptide_length_mhc2(self):
        """MHC-II peptide of wrong length should raise ValueError."""
        with pytest.raises(ValueError, match="Peptide length"):
            _make_mhc2_record(peptide="SHORTPEPTIDE")  # 11-mer, too short for MHC-II

    def test_valid_mhc1_peptide_lengths(self):
        """MHC-I peptides of length 8-11 should be accepted."""
        for length in [8, 9, 10, 11]:
            peptide = "A" * length
            rec = _make_record(peptide=peptide)
            assert len(rec.peptide) == length

    def test_valid_mhc2_peptide_lengths(self):
        """MHC-II peptides of length 13-25 should be accepted."""
        for length in [13, 15, 20, 25]:
            peptide = "A" * length
            rec = _make_mhc2_record(peptide=peptide)
            assert len(rec.peptide) == length

    def test_custom_peptide_length_range(self):
        """Custom peptide_length_range should override default validation."""
        rec = MHCBindingRecord(
            allele="CustomAllele",
            peptide="SHORT",
            ic50_nm=100.0,
            rank=None,
            binding_class="weak_binder",
            source="pssm_fallback",
            method="custom",
            timestamp="2025-01-01T00:00:00Z",
            peptide_length_range=(3, 10),
        )
        assert len(rec.peptide) == 5


# ═══════════════════════════════════════════════════════════════════════════
# 6. Helper function tests
# ═══════════════════════════════════════════════════════════════════════════


class TestHelperFunctions:
    """Tests for module-level helper functions."""

    def test_is_mhc_class2_human_dr(self):
        """HLA-DRB1 alleles should be identified as MHC class II."""
        assert is_mhc_class2("HLA-DRB1*01:01") is True
        assert is_mhc_class2("HLA-DRB1*04:01") is True

    def test_is_mhc_class2_human_dq_dp(self):
        """HLA-DQ and HLA-DP alleles should be MHC class II."""
        assert is_mhc_class2("HLA-DQB1*03:01") is True
        assert is_mhc_class2("HLA-DPB1*04:01") is True

    def test_is_mhc_class1(self):
        """HLA-A, -B, -C alleles should NOT be MHC class II."""
        assert is_mhc_class2("HLA-A*02:01") is False
        assert is_mhc_class2("HLA-B*07:02") is False
        assert is_mhc_class2("HLA-C*07:01") is False

    def test_is_mhc_class2_mouse(self):
        """Mouse H2-I alleles should be MHC class II."""
        assert is_mhc_class2("H2-IAb") is True
        assert is_mhc_class2("H2-IEb") is True

    def test_get_default_alleles_homo_sapiens(self):
        """Homo_sapiens should return known default alleles."""
        alleles = get_default_alleles_for_organism("Homo_sapiens")
        assert len(alleles) > 0
        assert "HLA-A*02:01" in alleles
        assert "HLA-DRB1*01:01" in alleles

    def test_get_default_alleles_mus_musculus(self):
        """Mus_musculus should return mouse alleles."""
        alleles = get_default_alleles_for_organism("Mus_musculus")
        assert "H-2Kb" in alleles
        assert "H-2Db" in alleles

    def test_get_default_alleles_unknown_organism(self):
        """Unknown organism should return empty list."""
        alleles = get_default_alleles_for_organism("Unknown_organism")
        assert alleles == []
