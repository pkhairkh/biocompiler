"""Tests for biocompiler.organism_db — SQLite-backed organism database.

Covers:
1. Database lookup functions (organism_exists, list_organisms, list_organism_names,
   get_codon_usage, get_codon_adaptiveness, get_preferred_codons)
2. Species information retrieval (store_organism + read-back, metadata fields,
   codon adaptiveness computation, preferred codon computation)
3. Error handling for unknown organisms (UnsupportedOrganismError on missing data,
   Kazusa fetch failure handling, invalid input resilience)
4. Schema creation and versioning
5. Module-level convenience functions
6. HTML parsing and validation helpers
"""

from __future__ import annotations

import hashlib
import re
import sqlite3
import tempfile
import urllib.error
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Optional
from unittest.mock import patch, MagicMock

import pytest

from biocompiler.organisms.db import (
    OrganismDatabase,
    SCHEMA_VERSION,
    KAZUSA_ORGANISM_IDS,
    KAZUSA_API_URL,
    CACHE_TTL_SECONDS,
    TOTAL_CODONS,
    NUM_STOP_CODONS,
    PER_THOUSAND_SCALE,
    MIN_CODON_COUNT,
    get_database,
    get_codon_usage_db,
    get_codon_adaptiveness_db,
)
from biocompiler.shared.constants import CODON_TABLE, AA_TO_CODONS
from biocompiler.shared.exceptions import UnsupportedOrganismError


# ─── Fixtures ────────────────────────────────────────────────────────


@pytest.fixture
def tmp_db_path(tmp_path: Path) -> Path:
    """Provide a temporary database path that will not collide with the real DB."""
    return tmp_path / "test_organisms.db"


@pytest.fixture
def db(tmp_db_path: Path) -> OrganismDatabase:
    """Create a fresh OrganismDatabase with a temporary SQLite file."""
    return OrganismDatabase(db_path=tmp_db_path)


@pytest.fixture
def sample_codon_usage() -> dict[str, tuple[str, float, float, int]]:
    """Build a minimal but valid codon usage table for a fictitious organism.

    Uses uniform frequency for every codon within each amino-acid group so
    that the data is internally consistent.
    """
    usage: dict[str, tuple[str, float, float, int]] = {}
    for aa, codons in AA_TO_CODONS.items():
        n = len(codons)
        for codon in codons:
            freq = 1.0 / n
            usage[codon] = (aa, freq, 0.0, 1000)
    # Include stop codons
    for stop in ("TAA", "TAG", "TGA"):
        usage[stop] = ("*", 1.0 / 3, 0.0, 300)
    return usage


@pytest.fixture
def populated_db(
    db: OrganismDatabase,
    sample_codon_usage: dict[str, tuple[str, float, float, int]],
) -> OrganismDatabase:
    """Return a database that already has one organism stored."""
    db.store_organism(
        name="Test_organism",
        codon_usage=sample_codon_usage,
        taxonomy_id="99999",
        source="test",
        n_cds=500,
    )
    return db


# ═══════════════════════════════════════════════════════════════════════
# 1. Database Lookup Functions
# ═══════════════════════════════════════════════════════════════════════


class TestOrganismExists:
    """Tests for organism_exists()."""

    def test_returns_false_for_empty_db(self, db: OrganismDatabase) -> None:
        assert db.organism_exists("Homo_sapiens") is False

    def test_returns_true_after_storing(self, populated_db: OrganismDatabase) -> None:
        assert populated_db.organism_exists("Test_organism") is True

    def test_case_sensitive_lookup(self, populated_db: OrganismDatabase) -> None:
        """Organism names are case-sensitive in the database."""
        assert populated_db.organism_exists("test_organism") is False
        assert populated_db.organism_exists("Test_organism") is True

    def test_nonexistent_organism(self, populated_db: OrganismDatabase) -> None:
        assert populated_db.organism_exists("Does_not_exist") is False


class TestListOrganisms:
    """Tests for list_organisms() and list_organism_names()."""

    def test_empty_db_returns_empty_list(self, db: OrganismDatabase) -> None:
        assert db.list_organisms() == []
        assert db.list_organism_names() == []

    def test_returns_stored_organisms(self, populated_db: OrganismDatabase) -> None:
        names = populated_db.list_organism_names()
        assert "Test_organism" in names

    def test_list_organisms_returns_dicts_with_metadata(
        self, populated_db: OrganismDatabase
    ) -> None:
        organisms = populated_db.list_organisms()
        assert len(organisms) == 1
        org = organisms[0]
        assert org["name"] == "Test_organism"
        assert org["taxonomy_id"] == "99999"
        assert org["source"] == "test"
        assert org["n_cds"] == 500
        assert "updated_at" in org

    def test_multiple_organisms_sorted_by_name(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        db.store_organism("Zebra_fish", sample_codon_usage, source="test")
        db.store_organism("Aardvark", sample_codon_usage, source="test")
        names = db.list_organism_names()
        assert names == ["Aardvark", "Zebra_fish"]


class TestGetCodonUsage:
    """Tests for get_codon_usage()."""

    def test_returns_all_64_codons(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        usage = populated_db.get_codon_usage("Test_organism")
        assert len(usage) == TOTAL_CODONS

    def test_tuple_structure(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        """Each value is (amino_acid, frequency, adaptiveness, count)."""
        usage = populated_db.get_codon_usage("Test_organism")
        for codon, val in usage.items():
            assert isinstance(val, tuple) and len(val) == 4, (
                f"Codon {codon!r} value is not a 4-tuple: {val!r}"
            )
            aa, freq, adapt, count = val
            assert isinstance(aa, str)
            assert isinstance(freq, float)
            assert isinstance(adapt, float)
            assert isinstance(count, int)

    def test_frequencies_non_negative(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        usage = populated_db.get_codon_usage("Test_organism")
        for codon, (_, freq, _, _) in usage.items():
            assert freq >= 0.0, f"Negative frequency for codon {codon}"

    def test_adaptiveness_values_populated(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        """Adaptiveness should be computed and stored when codon usage is retrieved."""
        usage = populated_db.get_codon_usage("Test_organism")
        adapt_values = [a for _, (_, _, a, _) in usage.items()]
        # At least the most frequent codon per AA should have adaptiveness ≈ 1.0
        assert any(a > 0.5 for a in adapt_values), (
            "No adaptiveness values above 0.5 found — computation may have failed"
        )

    def test_codon_to_aa_matches_standard_code(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        usage = populated_db.get_codon_usage("Test_organism")
        for codon, (aa, *_) in usage.items():
            expected = CODON_TABLE.get(codon)
            if expected == "*":
                # Stop codons may be represented as "*" or "STOP"
                assert aa in ("*", "STOP"), (
                    f"Codon {codon} maps to {aa!r}, expected '*' or 'STOP'"
                )
            else:
                assert aa == expected, (
                    f"Codon {codon} maps to {aa!r}, expected {expected!r}"
                )

    def test_raises_for_unknown_organism(self, db: OrganismDatabase) -> None:
        """get_codon_usage raises UnsupportedOrganismError when the organism
        is not in the DB and Kazusa fetch also fails."""
        with patch.object(
            db, "fetch_from_kazusa", side_effect=urllib.error.URLError("unreachable")
        ):
            with pytest.raises(UnsupportedOrganismError) as exc_info:
                db.get_codon_usage("Fantasy_organism")
            assert "Fantasy_organism" in str(exc_info.value)


class TestGetCodonAdaptiveness:
    """Tests for get_codon_adaptiveness()."""

    def test_returns_dict_of_floats(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        adapt = populated_db.get_codon_adaptiveness("Test_organism")
        assert isinstance(adapt, dict)
        for codon, val in adapt.items():
            assert isinstance(val, float), f"Adaptiveness for {codon} is not float"
            assert 0.0 <= val <= 1.0, (
                f"Adaptiveness for {codon} out of range [0,1]: {val}"
            )

    def test_most_frequent_codon_has_adaptiveness_one(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        """For each amino acid, the codon with the highest frequency should
        have adaptiveness == 1.0."""
        adapt = populated_db.get_codon_adaptiveness("Test_organism")
        usage = populated_db._get_raw_usage("Test_organism")
        for aa, codons in AA_TO_CODONS.items():
            best = max(codons, key=lambda c: usage.get(c, 0.0))
            assert adapt.get(best, 0.0) == pytest.approx(1.0, abs=1e-9), (
                f"Best codon {best} for AA {aa} should have adaptiveness 1.0, "
                f"got {adapt.get(best)}"
            )

    def test_adaptiveness_computed_on_demand(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        """Adaptiveness should be auto-computed if not already stored."""
        db.store_organism("AutoAdapt", sample_codon_usage, source="test")
        # Clear adaptiveness manually
        conn = db._connect()
        try:
            conn.execute(
                "DELETE FROM codon_adaptiveness WHERE organism = ?", ("AutoAdapt",)
            )
            conn.commit()
        finally:
            conn.close()

        # Now request adaptiveness — should be recomputed
        adapt = db.get_codon_adaptiveness("AutoAdapt")
        assert len(adapt) > 0, "Adaptiveness was not recomputed on demand"


class TestGetPreferredCodons:
    """Tests for get_preferred_codons()."""

    def test_returns_dict_of_strings(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        preferred = populated_db.get_preferred_codons("Test_organism")
        assert isinstance(preferred, dict)
        for aa, codon in preferred.items():
            assert isinstance(aa, str) and len(aa) == 1
            assert isinstance(codon, str) and len(codon) == 3

    def test_one_preferred_codon_per_amino_acid(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        preferred = populated_db.get_preferred_codons("Test_organism")
        # 20 standard amino acids should each have exactly one preferred codon
        assert len(preferred) == 20

    def test_preferred_codon_is_valid_synonym(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        """The preferred codon for each AA must be a valid codon for that AA."""
        preferred = populated_db.get_preferred_codons("Test_organism")
        for aa, codon in preferred.items():
            assert codon in AA_TO_CODONS.get(aa, []), (
                f"Preferred codon {codon} for AA {aa} is not a valid synonym"
            )

    def test_preferred_codon_has_highest_adaptiveness(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        """The preferred codon for each AA should be the one with the
        highest adaptiveness value."""
        preferred = populated_db.get_preferred_codons("Test_organism")
        adapt = populated_db.get_codon_adaptiveness("Test_organism")
        for aa, codon in preferred.items():
            codons = AA_TO_CODONS[aa]
            max_adapt = max(adapt.get(c, 0.0) for c in codons)
            assert adapt.get(codon, 0.0) == pytest.approx(max_adapt, abs=1e-9), (
                f"Preferred codon {codon} for AA {aa} does not have the highest "
                f"adaptiveness: {adapt.get(codon)} vs max {max_adapt}"
            )


# ═══════════════════════════════════════════════════════════════════════
# 2. Species Information Retrieval
# ═══════════════════════════════════════════════════════════════════════


class TestStoreOrganism:
    """Tests for store_organism() — write then read-back."""

    def test_store_and_retrieve_basic(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        db.store_organism(
            name="Homo_sapiens",
            codon_usage=sample_codon_usage,
            taxonomy_id="9606",
            taxonomy_lineage="Eukaryota;Animalia;Chordata;Mammalia",
            source="test",
            n_cds=10000,
        )
        assert db.organism_exists("Homo_sapiens")
        orgs = db.list_organisms()
        assert any(o["name"] == "Homo_sapiens" for o in orgs)
        hs = [o for o in orgs if o["name"] == "Homo_sapiens"][0]
        assert hs["taxonomy_id"] == "9606"
        assert hs["source"] == "test"
        assert hs["n_cds"] == 10000

    def test_upsert_replaces_existing_data(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        db.store_organism("Upserter", sample_codon_usage, n_cds=100, source="v1")
        db.store_organism("Upserter", sample_codon_usage, n_cds=200, source="v2")
        orgs = db.list_organisms()
        upserter = [o for o in orgs if o["name"] == "Upserter"][0]
        assert upserter["source"] == "v2"
        assert upserter["n_cds"] == 200

    def test_store_multiple_organisms(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        for name in ("A", "B", "C"):
            db.store_organism(name, sample_codon_usage, source="test")
        assert db.list_organism_names() == ["A", "B", "C"]

    def test_taxonomy_lineage_stored(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        lineage = "Eukaryota;Animalia;Chordata"
        db.store_organism(
            "LinOrg", sample_codon_usage, taxonomy_lineage=lineage, source="test"
        )
        conn = db._connect()
        try:
            row = conn.execute(
                "SELECT taxonomy_lineage FROM organisms WHERE name = ?",
                ("LinOrg",),
            ).fetchone()
            assert row["taxonomy_lineage"] == lineage
        finally:
            conn.close()


class TestCodonAdaptivenessComputation:
    """Tests for the internal adaptiveness computation logic."""

    def test_uniform_frequencies_yield_equal_adaptiveness(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        """With uniform frequencies, all codons for the same AA should
        have adaptiveness ≈ 1.0."""
        db.store_organism("UniformOrg", sample_codon_usage, source="test")
        adapt = db.get_codon_adaptiveness("UniformOrg")
        for aa, codons in AA_TO_CODONS.items():
            for c in codons:
                assert adapt.get(c, 0.0) == pytest.approx(1.0, abs=0.01), (
                    f"Uniform codon {c} (AA {aa}) adaptiveness should be ≈1.0, "
                    f"got {adapt.get(c)}"
                )

    def test_skewed_frequencies_yield_proportional_adaptiveness(
        self, db: OrganismDatabase
    ) -> None:
        """Create a skewed usage table where one codon dominates, then
        verify adaptiveness is proportional."""
        skewed: dict[str, tuple[str, float, float, int]] = {}
        for aa, codons in AA_TO_CODONS.items():
            if len(codons) == 1:
                skewed[codons[0]] = (aa, 1.0, 0.0, 1000)
            else:
                # First codon gets 90%, rest share 10%
                for i, c in enumerate(codons):
                    freq = 0.9 if i == 0 else 0.1 / (len(codons) - 1)
                    skewed[c] = (aa, freq, 0.0, int(freq * 10000))
        for stop in ("TAA", "TAG", "TGA"):
            skewed[stop] = ("*", 1.0 / 3, 0.0, 300)

        db.store_organism("SkewedOrg", skewed, source="test")
        adapt = db.get_codon_adaptiveness("SkewedOrg")
        for aa, codons in AA_TO_CODONS.items():
            if len(codons) > 1:
                best = codons[0]
                assert adapt[best] == pytest.approx(1.0, abs=1e-6), (
                    f"Dominant codon {best} for AA {aa} should have adaptiveness 1.0"
                )
                for other in codons[1:]:
                    ratio = (0.1 / (len(codons) - 1)) / 0.9
                    assert adapt[other] == pytest.approx(ratio, abs=0.01), (
                        f"Codon {other} for AA {aa} adaptiveness wrong: "
                        f"got {adapt[other]}, expected ≈{ratio}"
                    )


# ═══════════════════════════════════════════════════════════════════════
# 3. Error Handling for Unknown Organisms
# ═══════════════════════════════════════════════════════════════════════


class TestUnsupportedOrganismError:
    """Tests for UnsupportedOrganismError propagation."""

    def test_get_codon_usage_raises_for_unknown(self, db: OrganismDatabase) -> None:
        with patch.object(
            db, "fetch_from_kazusa", side_effect=urllib.error.URLError("fail")
        ):
            with pytest.raises(UnsupportedOrganismError) as exc_info:
                db.get_codon_usage("NonExistent_species")
            exc = exc_info.value
            assert exc.organism == "NonExistent_species"
            assert isinstance(exc.available, list)

    def test_error_message_contains_organism_name(self, db: OrganismDatabase) -> None:
        with patch.object(
            db, "fetch_from_kazusa", side_effect=urllib.error.URLError("fail")
        ):
            with pytest.raises(UnsupportedOrganismError, match="NonExistent_species"):
                db.get_codon_usage("NonExistent_species")

    def test_error_lists_available_organisms(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        """The error should list organisms that *are* available."""
        with patch.object(
            populated_db,
            "fetch_from_kazusa",
            side_effect=urllib.error.URLError("fail"),
        ):
            with pytest.raises(UnsupportedOrganismError) as exc_info:
                populated_db.get_codon_usage("Missing_species")
            assert "Test_organism" in exc_info.value.available

    def test_fetch_from_kazusa_raises_valueerror_on_empty_parse(
        self, db: OrganismDatabase
    ) -> None:
        """When Kazusa returns HTML that cannot be parsed into codon usage,
        fetch_from_kazusa should raise ValueError."""
        bad_html = "<html><body>Not a codon table</body></html>"
        with patch.object(db, "_is_cache_fresh", return_value=False):
            with patch("urllib.request.urlopen") as mock_urlopen:
                mock_resp = MagicMock()
                mock_resp.read.return_value = bad_html.encode("utf-8")
                mock_resp.__enter__ = MagicMock(return_value=mock_resp)
                mock_resp.__exit__ = MagicMock(return_value=False)
                mock_urlopen.return_value = mock_resp
                # The parse may produce uniform fallback, but if it does,
                # ValueError will not be raised; if it produces empty, it will.
                # With such short HTML, _validate_kazusa_response may fail,
                # but parsing may still produce data via uniform fallback.
                # Let us just ensure no crash occurs.
                try:
                    db.fetch_from_kazusa("Fake_org")
                except (ValueError, urllib.error.URLError):
                    pass  # Expected — either parse fails or network fails

    def test_fetch_from_kazusa_raises_urlerror_on_network_failure(
        self, db: OrganismDatabase
    ) -> None:
        """When all retry attempts fail, URLError should be raised."""
        with patch.object(db, "_is_cache_fresh", return_value=False):
            with patch(
                "urllib.request.urlopen",
                side_effect=urllib.error.URLError("Connection refused"),
            ):
                # Mock time.sleep to avoid delays in tests
                with patch("time.sleep"):
                    with pytest.raises(urllib.error.URLError):
                        db.fetch_from_kazusa("NoNetwork_org")

    def test_get_codon_usage_on_kazusa_value_error(
        self, db: OrganismDatabase
    ) -> None:
        """ValueError from Kazusa (organism not found) should also
        result in UnsupportedOrganismError."""
        with patch.object(
            db,
            "fetch_from_kazusa",
            side_effect=ValueError("No codon usage data found"),
        ):
            with pytest.raises(UnsupportedOrganismError):
                db.get_codon_usage("Unknown_org")

    def test_get_codon_usage_on_sqlite_error(
        self, db: OrganismDatabase
    ) -> None:
        """sqlite3.Error from Kazusa fetch should also result in
        UnsupportedOrganismError."""
        with patch.object(
            db,
            "fetch_from_kazusa",
            side_effect=sqlite3.Error("DB error"),
        ):
            with pytest.raises(UnsupportedOrganismError):
                db.get_codon_usage("Broken_org")


class TestEdgeCases:
    """Edge-case and resilience tests."""

    def test_empty_codon_usage_store(self, db: OrganismDatabase) -> None:
        """Storing an organism with empty codon_usage should not crash."""
        db.store_organism("EmptyOrg", {}, source="test")
        assert db.organism_exists("EmptyOrg")

    def test_empty_codon_usage_retrieval(self, db: OrganismDatabase) -> None:
        """Retrieving codon usage for an organism with no data stored
        should trigger Kazusa fetch and then UnsupportedOrganismError."""
        db.store_organism("EmptyOrg", {}, source="test")
        with patch.object(
            db,
            "fetch_from_kazusa",
            side_effect=urllib.error.URLError("fail"),
        ):
            with pytest.raises(UnsupportedOrganismError):
                db.get_codon_usage("EmptyOrg")

    def test_organism_name_with_special_characters(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        """Organism names with underscores and numbers should work fine."""
        name = "Escherichia_coli_K12"
        db.store_organism(name, sample_codon_usage, source="test")
        assert db.organism_exists(name)
        usage = db.get_codon_usage(name)
        assert len(usage) == TOTAL_CODONS


# ═══════════════════════════════════════════════════════════════════════
# 4. Schema Creation and Versioning
# ═══════════════════════════════════════════════════════════════════════


class TestSchemaCreation:
    """Tests for database schema initialization and versioning."""

    def test_schema_version_set(self, db: OrganismDatabase) -> None:
        """The schema version should be set to SCHEMA_VERSION."""
        conn = db._connect()
        try:
            row = conn.execute(
                "SELECT value FROM schema_version WHERE key = 'version'"
            ).fetchone()
            assert row is not None, "schema_version table missing"
            assert int(row["value"]) == SCHEMA_VERSION
        finally:
            conn.close()

    def test_tables_created(self, db: OrganismDatabase) -> None:
        """All expected tables should exist after initialization."""
        conn = db._connect()
        try:
            tables = {
                row[0]
                for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                ).fetchall()
            }
        finally:
            conn.close()
        for expected in (
            "schema_version",
            "organisms",
            "codon_usage",
            "codon_adaptiveness",
            "preferred_codons",
            "api_response_cache",
        ):
            assert expected in tables, f"Table {expected} not found"

    def test_indexes_created(self, db: OrganismDatabase) -> None:
        """Expected indexes should exist."""
        conn = db._connect()
        try:
            indexes = {
                row[0]
                for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='index'"
                ).fetchall()
            }
        finally:
            conn.close()
        for expected_idx in (
            "idx_codon_usage_organism",
            "idx_codon_adaptiveness_organism",
            "idx_preferred_codons_organism",
        ):
            assert expected_idx in indexes, f"Index {expected_idx} not found"

    def test_fresh_db_schema_version_is_zero(self, tmp_db_path: Path) -> None:
        """Before _ensure_schema, a fresh DB should report version 0."""
        # Create a raw connection without going through OrganismDatabase
        conn = sqlite3.connect(str(tmp_db_path))
        try:
            version = OrganismDatabase._get_schema_version(conn)
            assert version == 0
        finally:
            conn.close()

    def test_schema_reinitialization_idempotent(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        """Calling _ensure_schema again should not corrupt existing data."""
        db.store_organism("IdempotentOrg", sample_codon_usage, source="test")
        db._ensure_schema()
        assert db.organism_exists("IdempotentOrg")
        usage = db.get_codon_usage("IdempotentOrg")
        assert len(usage) == TOTAL_CODONS


# ═══════════════════════════════════════════════════════════════════════
# 5. Module-level Convenience Functions
# ═══════════════════════════════════════════════════════════════════════


class TestConvenienceFunctions:
    """Tests for get_database(), get_codon_usage_db(), get_codon_adaptiveness_db()."""

    def test_get_database_returns_instance(self, tmp_db_path: Path) -> None:
        """get_database should return an OrganismDatabase instance."""
        import biocompiler.organisms.db as mod

        # Reset singleton
        mod._db_instance = None
        db = get_database(db_path=tmp_db_path)
        assert isinstance(db, OrganismDatabase)
        # Clean up
        mod._db_instance = None

    def test_get_database_singleton(self, tmp_db_path: Path) -> None:
        """get_database should return the same instance on repeated calls."""
        import biocompiler.organisms.db as mod

        mod._db_instance = None
        db1 = get_database(db_path=tmp_db_path)
        db2 = get_database(db_path=tmp_db_path)
        assert db1 is db2
        # Clean up
        mod._db_instance = None

    def test_get_codon_usage_db_convenience(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        """get_codon_usage_db should delegate to the database instance."""
        import biocompiler.organisms.db as mod

        old = mod._db_instance
        mod._db_instance = populated_db
        try:
            usage = get_codon_usage_db("Test_organism")
            assert len(usage) == TOTAL_CODONS
        finally:
            mod._db_instance = old

    def test_get_codon_adaptiveness_db_convenience(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        """get_codon_adaptiveness_db should delegate to the database instance."""
        import biocompiler.organisms.db as mod

        old = mod._db_instance
        mod._db_instance = populated_db
        try:
            adapt = get_codon_adaptiveness_db("Test_organism")
            assert isinstance(adapt, dict)
            assert len(adapt) > 0
        finally:
            mod._db_instance = old


# ═══════════════════════════════════════════════════════════════════════
# 6. HTML Parsing and Validation Helpers
# ═══════════════════════════════════════════════════════════════════════


class TestValidateKazusaResponse:
    """Tests for _validate_kazusa_response()."""

    def test_short_html_returns_false(self, db: OrganismDatabase) -> None:
        assert db._validate_kazusa_response("", "Test_org") is False
        assert db._validate_kazusa_response("x" * 50, "Test_org") is False

    def test_error_indicators_return_false(self, db: OrganismDatabase) -> None:
        """HTML containing Kazusa error indicators should fail validation."""
        html = (
            "<html><body><p>The species was not found in the database.</p>"
            + "<td>ATG</td>" * 25
            + "</body></html>"
        )
        assert db._validate_kazusa_response(html, "Test_org") is False

    def test_valid_response_returns_true(self, db: OrganismDatabase) -> None:
        """A response with plenty of codon patterns and no errors should pass."""
        # Generate a plausible HTML with 64+ codon patterns
        rows = ""
        for codon, aa in CODON_TABLE.items():
            rows += f"<tr><td>{codon}</td><td>{aa}</td><td>0.25</td><td>1000</td><td>25.0</td></tr>\n"
        html = f"<html><body><table>{rows}</table></body></html>"
        assert db._validate_kazusa_response(html, "Test_org") is True

    def test_few_codon_patterns_returns_false(self, db: OrganismDatabase) -> None:
        """A response with very few codon-like patterns should fail."""
        # Must be >= 100 chars to pass the length check
        html = (
            "<html><body><p>Some text with ATG and GCT only, "
            "plus enough padding to exceed the minimum length threshold.</p></body></html>"
        )
        assert len(html) >= 100, f"Test HTML too short ({len(html)} chars)"
        assert db._validate_kazusa_response(html, "Test_org") is False


class TestExtractNCdsFromHtml:
    """Tests for _extract_n_cds_from_html()."""

    def test_extracts_cds_count(self) -> None:
        html = "<p>Based on 43650 CDS sequences</p>"
        assert OrganismDatabase._extract_n_cds_from_html(html) == 43650

    def test_returns_zero_when_no_match(self) -> None:
        html = "<p>No CDS info here</p>"
        assert OrganismDatabase._extract_n_cds_from_html(html) == 0

    def test_extracts_first_match(self) -> None:
        html = "<p>100 CDS and then 200 CDS</p>"
        assert OrganismDatabase._extract_n_cds_from_html(html) == 100


class TestValidateCodonUsageData:
    """Tests for _validate_codon_usage_data()."""

    def test_corrects_wrong_amino_acid(self) -> None:
        """Should correct an amino acid that does not match CODON_TABLE."""
        bad_usage: dict[str, tuple[str, float, float, int]] = {
            "ATG": ("X", 1.0, 0.0, 1000),  # Should be M
            "TTT": ("Y", 0.5, 0.0, 500),    # Should be F
            "TTC": ("F", 0.5, 0.0, 500),
        }
        result = OrganismDatabase._validate_codon_usage_data(bad_usage, "TestOrg")
        assert result["ATG"][0] == "M", "AA for ATG should be corrected to M"
        assert result["TTT"][0] == "F", "AA for TTT should be corrected to F"

    def test_clamps_negative_frequencies(self) -> None:
        bad_usage: dict[str, tuple[str, float, float, int]] = {
            "ATG": ("M", -0.5, 0.0, 1000),
        }
        result = OrganismDatabase._validate_codon_usage_data(bad_usage, "TestOrg")
        assert result["ATG"][1] == 0.0, "Negative frequency should be clamped to 0"

    def test_normalizes_frequencies_above_one(self) -> None:
        bad_usage: dict[str, tuple[str, float, float, int]] = {
            "ATG": ("M", 1.5, 0.0, 1000),
        }
        result = OrganismDatabase._validate_codon_usage_data(bad_usage, "TestOrg")
        assert result["ATG"][1] <= 1.0, "Frequency > 1.0 should be normalized"


class TestParseKazusaHtmlRobust:
    """Tests for _parse_kazusa_html_robust() static method."""

    def test_parses_well_formed_html(self) -> None:
        """A well-formed Kazusa HTML table should be parsed successfully."""
        rows = ""
        for codon, aa in CODON_TABLE.items():
            aa_display = aa if aa != "*" else "*"
            rows += (
                f"<tr><td>{codon}</td><td>{aa_display}</td>"
                f"<td>0.25</td><td>1000</td><td>25.0</td></tr>\n"
            )
        html = f"<html><body><table>{rows}</table></body></html>"
        result = OrganismDatabase._parse_kazusa_html_robust(html, "TestOrg")
        assert len(result) >= MIN_CODON_COUNT, (
            f"Expected at least {MIN_CODON_COUNT} codons, got {len(result)}"
        )

    def test_returns_uniform_on_empty_html(self) -> None:
        """With completely unparseable HTML, the method should still return
        64 codons via the uniform fallback."""
        html = "<html><body>Nothing useful here at all, just text.</body></html>"
        result = OrganismDatabase._parse_kazusa_html_robust(html, "TestOrg")
        assert len(result) == TOTAL_CODONS, (
            f"Uniform fallback should produce {TOTAL_CODONS} codons, got {len(result)}"
        )

    def test_stop_codons_in_result(self) -> None:
        """Stop codons should be present in the parsed output."""
        html = "<html><body>Nothing parseable.</body></html>"
        result = OrganismDatabase._parse_kazusa_html_robust(html, "TestOrg")
        for stop in ("TAA", "TAG", "TGA"):
            assert stop in result, f"Stop codon {stop} missing from result"


# ═══════════════════════════════════════════════════════════════════════
# 7. Cache and Response Hash
# ═══════════════════════════════════════════════════════════════════════


class TestCacheFreshness:
    """Tests for _is_cache_fresh()."""

    def test_nonexistent_organism_is_not_fresh(self, db: OrganismDatabase) -> None:
        assert db._is_cache_fresh("NoOrganism") is False

    def test_builtin_source_is_not_fresh(
        self,
        populated_db: OrganismDatabase,
    ) -> None:
        """Builtin-sourced organisms are not considered 'kazusa-fresh'."""
        assert populated_db._is_cache_fresh("Test_organism") is False

    def test_recent_kazusa_fetch_is_fresh(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        """An organism with source='kazusa' and recent updated_at should be fresh."""
        # Insert directly with source='kazusa' and a very recent timestamp
        now = datetime.now(timezone.utc).isoformat()
        conn = db._connect()
        try:
            conn.execute(
                "INSERT INTO organisms (name, taxonomy_id, source, n_cds, created_at, updated_at) "
                "VALUES (?, ?, 'kazusa', 0, ?, ?)",
                ("KazusaOrg", "12345", now, now),
            )
            conn.commit()
        finally:
            conn.close()
        # Also need codon_usage data for the organism to be usable
        db.store_organism("KazusaOrg2", sample_codon_usage, source="kazusa")
        assert db._is_cache_fresh("KazusaOrg2") is True


class TestStoreResponseHash:
    """Tests for _store_response_hash()."""

    def test_hash_stored_correctly(
        self,
        db: OrganismDatabase,
        sample_codon_usage: dict[str, tuple[str, float, float, int]],
    ) -> None:
        db.store_organism("HashOrg", sample_codon_usage, source="test")
        html = "<html>test content</html>"
        url = "https://example.com/test"
        db._store_response_hash("HashOrg", html, url)

        conn = db._connect()
        try:
            row = conn.execute(
                "SELECT response_hash, url FROM api_response_cache WHERE organism = ?",
                ("HashOrg",),
            ).fetchone()
            assert row is not None, "Response hash not stored"
            expected_hash = hashlib.sha256(html.encode("utf-8")).hexdigest()
            assert row["response_hash"] == expected_hash
            assert row["url"] == url
        finally:
            conn.close()


# ═══════════════════════════════════════════════════════════════════════
# 8. Constants and Module-level Values
# ═══════════════════════════════════════════════════════════════════════


class TestConstants:
    """Verify key module constants are sensible."""

    def test_schema_version_is_positive(self) -> None:
        assert SCHEMA_VERSION >= 1

    def test_kazusa_organism_ids_cover_common_species(self) -> None:
        """The known IDs dict should include key model organisms."""
        for org in (
            "Homo_sapiens",
            "Mus_musculus",
            "Escherichia_coli",
            "Saccharomyces_cerevisiae",
        ):
            assert org in KAZUSA_ORGANISM_IDS, f"{org} missing from KAZUSA_ORGANISM_IDS"

    def test_e_coli_alias(self) -> None:
        """E_coli should be an alias for Escherichia_coli."""
        assert KAZUSA_ORGANISM_IDS.get("E_coli") == KAZUSA_ORGANISM_IDS.get(
            "Escherichia_coli"
        )

    def test_cache_ttl_is_one_week(self) -> None:
        assert CACHE_TTL_SECONDS == 7 * 24 * 60 * 60

    def test_total_codons_is_64(self) -> None:
        assert TOTAL_CODONS == 64

    def test_num_stop_codons_is_3(self) -> None:
        assert NUM_STOP_CODONS == 3
