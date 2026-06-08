"""
Tests for thread-safe Sharp-Li table access.

Validates:
  - get_sharp_li_table() returns correct data for known organisms
  - get_sharp_li_table() raises ValueError for unknown organisms
  - set_sharp_li_table() atomically updates the table
  - Concurrent access from multiple threads does not raise errors
  - _SharpLiState class provides lock-protected access to all 14 variables
"""

from __future__ import annotations

import threading
from typing import Dict

import pytest

from biocompiler.sharp_li_tables import (
    _SharpLiState,
    get_sharp_li_table,
    set_sharp_li_table,
)


# ────────────────────────────────────────────────────────────
# Basic functional tests
# ────────────────────────────────────────────────────────────


class TestGetSharpLiTable:
    """Tests for get_sharp_li_table()."""

    def test_ecoli_alias(self) -> None:
        """'ecoli' should resolve to E. coli CAI weights."""
        table = get_sharp_li_table("ecoli")
        assert isinstance(table, dict)
        assert len(table) > 0
        # CTG should be present (Leu codon, highly preferred in E. coli)
        assert "CTG" in table

    def test_ecoli_canonical(self) -> None:
        """'Escherichia_coli' should work as a key."""
        table = get_sharp_li_table("Escherichia_coli")
        assert isinstance(table, dict)
        assert len(table) > 0

    def test_yeast_alias(self) -> None:
        """'yeast' should resolve to S. cerevisiae CAI weights."""
        table = get_sharp_li_table("yeast")
        assert isinstance(table, dict)
        assert len(table) > 0
        # AGA should be present (Arg codon, preferred in yeast)
        assert "AGA" in table

    def test_yeast_canonical(self) -> None:
        """'Saccharomyces_cerevisiae' should work as a key."""
        table = get_sharp_li_table("Saccharomyces_cerevisiae")
        assert isinstance(table, dict)
        assert len(table) > 0

    def test_returns_copy(self) -> None:
        """Each call should return a fresh dict (not the internal state)."""
        t1 = get_sharp_li_table("ecoli")
        t2 = get_sharp_li_table("ecoli")
        assert t1 is not t2
        assert t1 == t2

    def test_unknown_organism_raises(self) -> None:
        """ValueError for unsupported organism."""
        with pytest.raises(ValueError, match="No Sharp-Li CAI weights"):
            get_sharp_li_table("nonexistent_organism")


class TestSetSharpLiTable:
    """Tests for set_sharp_li_table()."""

    def test_update_ecoli_table(self) -> None:
        """set_sharp_li_table should update the CAI weights."""
        original = get_sharp_li_table("ecoli")
        # Create a modified table
        modified = dict(original)
        modified["CTG"] = 0.5  # change one value

        try:
            set_sharp_li_table("ecoli", modified)
            updated = get_sharp_li_table("ecoli")
            assert updated["CTG"] == 0.5
        finally:
            # Restore original
            set_sharp_li_table("ecoli", original)

    def test_update_does_not_affect_input(self) -> None:
        """Mutating the input dict after set_sharp_li_table should
        not affect the stored state (defensive copy)."""
        original = get_sharp_li_table("ecoli")
        table = dict(original)
        table["CTG"] = 0.5

        try:
            set_sharp_li_table("ecoli", table)
            # Mutate input after setting
            table["CTG"] = 0.99
            stored = get_sharp_li_table("ecoli")
            assert stored["CTG"] == 0.5  # not 0.99
        finally:
            set_sharp_li_table("ecoli", original)


# ────────────────────────────────────────────────────────────
# _SharpLiState class tests
# ────────────────────────────────────────────────────────────


class TestSharpLiState:
    """Tests for the _SharpLiState class."""

    def test_all_14_accessors_return_copies(self) -> None:
        """Every accessor should return a fresh copy, not the internal state."""
        state = _SharpLiState()

        # Codon infrastructure (3)
        assert state.get_codon_table() is not state.get_codon_table()
        assert state.get_aa_to_codons() is not state.get_aa_to_codons()
        assert state.get_stop_codons() is not state.get_stop_codons()

        # E. coli (3)
        assert state.get_ecoli_reference_genes() is not state.get_ecoli_reference_genes()
        assert state.get_ecoli_codon_usage() is not state.get_ecoli_codon_usage()
        assert state.get_ecoli_cai_weights() is not state.get_ecoli_cai_weights()

        # Yeast (3)
        assert state.get_yeast_reference_genes() is not state.get_yeast_reference_genes()
        assert state.get_yeast_codon_usage() is not state.get_yeast_codon_usage()
        assert state.get_yeast_cai_weights() is not state.get_yeast_cai_weights()

        # Combined / published (5)
        assert state.get_published_cai() is not state.get_published_cai()
        assert state.get_reference_weights() is not state.get_reference_weights()
        assert state.get_reference_genes() is not state.get_reference_genes()
        assert state.get_codon_usage() is not state.get_codon_usage()
        assert state.get_cai_weights() is not state.get_cai_weights()

    def test_get_table_returns_copy(self) -> None:
        """get_table should return a copy."""
        state = _SharpLiState()
        t1 = state.get_table("ecoli")
        t2 = state.get_table("ecoli")
        assert t1 is not t2
        assert t1 == t2

    def test_get_table_unknown_returns_none(self) -> None:
        """get_table returns None for unknown organisms."""
        state = _SharpLiState()
        assert state.get_table("nonexistent") is None

    def test_set_table_updates_all_registries(self) -> None:
        """set_table should update cai_weights, codon_usage, and reference_weights."""
        state = _SharpLiState()
        original = state.get_table("Escherichia_coli")
        modified = dict(original)
        modified["CTG"] = 0.42

        state.set_table("Escherichia_coli", modified)

        # Check all three registries
        assert state.get_table("Escherichia_coli")["CTG"] == 0.42
        assert state.get_codon_usage()["Escherichia_coli"]["CTG"] == 0.42
        assert state.get_reference_weights()["Escherichia_coli"]["sharp_li"]["CTG"] == 0.42


# ────────────────────────────────────────────────────────────
# Thread-safety tests
# ────────────────────────────────────────────────────────────


class TestConcurrentAccess:
    """Concurrent access stress tests."""

    def test_concurrent_read_access(self) -> None:
        """Multiple threads reading concurrently should not raise errors."""
        errors: list[Exception] = []

        def access_table() -> None:
            try:
                for _ in range(100):
                    table = get_sharp_li_table("ecoli")
                    assert table is not None
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=access_table) for _ in range(10)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert len(errors) == 0, f"Errors during concurrent reads: {errors}"

    def test_concurrent_read_write_access(self) -> None:
        """Mix of reads and writes should not raise errors."""
        errors: list[Exception] = []
        original = get_sharp_li_table("ecoli")
        barrier = threading.Barrier(5)

        def reader() -> None:
            barrier.wait()
            try:
                for _ in range(50):
                    table = get_sharp_li_table("ecoli")
                    assert table is not None
            except Exception as e:
                errors.append(e)

        def writer() -> None:
            barrier.wait()
            try:
                for i in range(50):
                    modified = dict(original)
                    modified["CTG"] = 0.5 + (i % 10) * 0.01
                    set_sharp_li_table("ecoli", modified)
            except Exception as e:
                errors.append(e)
            finally:
                set_sharp_li_table("ecoli", original)

        threads = [
            threading.Thread(target=reader) for _ in range(4)
        ] + [threading.Thread(target=writer)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        # Restore original state
        set_sharp_li_table("ecoli", original)
        assert len(errors) == 0, f"Errors during concurrent read/write: {errors}"

    def test_concurrent_multi_organism_access(self) -> None:
        """Concurrent access to different organisms should be safe."""
        errors: list[Exception] = []

        def access_ecoli() -> None:
            try:
                for _ in range(100):
                    table = get_sharp_li_table("ecoli")
                    assert table is not None
            except Exception as e:
                errors.append(e)

        def access_yeast() -> None:
            try:
                for _ in range(100):
                    table = get_sharp_li_table("yeast")
                    assert table is not None
            except Exception as e:
                errors.append(e)

        threads = [
            threading.Thread(target=access_ecoli) for _ in range(5)
        ] + [
            threading.Thread(target=access_yeast) for _ in range(5)
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert len(errors) == 0, f"Errors during multi-organism access: {errors}"

    def test_concurrent_state_class_access(self) -> None:
        """Direct _SharpLiState access from multiple threads should be safe."""
        state = _SharpLiState()
        errors: list[Exception] = []

        def read_all_accessors() -> None:
            try:
                for _ in range(50):
                    state.get_codon_table()
                    state.get_aa_to_codons()
                    state.get_stop_codons()
                    state.get_ecoli_reference_genes()
                    state.get_ecoli_codon_usage()
                    state.get_ecoli_cai_weights()
                    state.get_yeast_reference_genes()
                    state.get_yeast_codon_usage()
                    state.get_yeast_cai_weights()
                    state.get_published_cai()
                    state.get_reference_weights()
                    state.get_reference_genes()
                    state.get_codon_usage()
                    state.get_cai_weights()
                    state.get_table("ecoli")
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=read_all_accessors) for _ in range(10)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert len(errors) == 0, f"Errors during state class concurrent access: {errors}"
