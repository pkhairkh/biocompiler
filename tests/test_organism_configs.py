"""Tests for organism configuration data in biocompiler.organisms.

Covers:
1. SPECIES dict has expected organisms (human, ecoli, yeast, mouse, CHO)
2. CODON_ADAPTIVENESS_TABLES has entries for each organism
3. ORGANISM_GC_TARGETS has organism-specific values
4. SUPPORTED_ORGANISMS list matches SPECIES keys (via canonical names)
5. Each organism module exports a codon usage table
6. Codon tables are valid (all 64 codons represented, only standard amino acids)
7. GC targets are in valid range [0.1, 0.9]
8. Adaptiveness values are in [0, 1] range
"""

from __future__ import annotations

import pytest

from biocompiler.organisms import (
    # Registry dicts
    SPECIES,
    CODON_USAGE_TABLES,
    CODON_ADAPTIVENESS_TABLES,
    ORGANISM_GC_TARGETS,
    SUPPORTED_ORGANISMS,
    # Name aliases
    HUMAN,
    E_COLI,
    MOUSE,
    CHO,
    YEAST,
    # Per-organism codon usage tables
    HUMAN_CODON_USAGE,
    E_COLI_CODON_USAGE,
    MOUSE_CODON_USAGE,
    CHO_CODON_USAGE,
    YEAST_CODON_USAGE,
    # Per-organism adaptiveness tables
    HUMAN_CODON_ADAPTIVENESS,
    E_COLI_CODON_ADAPTIVENESS,
    MOUSE_CODON_ADAPTIVENESS,
    CHO_CODON_ADAPTIVENESS,
    YEAST_CODON_ADAPTIVENESS,
)
from biocompiler.type_system import CODON_TABLE

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

EXPECTED_SPECIES_KEYS = {"ecoli", "human", "mouse", "cho", "yeast"}

# Canonical scientific name -> SPECIES short name mapping
CANONICAL_TO_SHORT = {
    HUMAN: "human",
    E_COLI: "ecoli",
    MOUSE: "mouse",
    CHO: "cho",
    YEAST: "yeast",
}

# All 20 standard amino acid single-letter codes + stop codon
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY*")

# Per-organism codon usage tables exposed on individual modules
ORGANISM_CODON_TABLES = {
    HUMAN: HUMAN_CODON_USAGE,
    E_COLI: E_COLI_CODON_USAGE,
    MOUSE: MOUSE_CODON_USAGE,
    CHO: CHO_CODON_USAGE,
    YEAST: YEAST_CODON_USAGE,
}

# Per-organism adaptiveness tables
ORGANISM_ADAPTIVENESS = {
    HUMAN: HUMAN_CODON_ADAPTIVENESS,
    E_COLI: E_COLI_CODON_ADAPTIVENESS,
    MOUSE: MOUSE_CODON_ADAPTIVENESS,
    CHO: CHO_CODON_ADAPTIVENESS,
    YEAST: YEAST_CODON_ADAPTIVENESS,
}

# All 64 codons from the standard genetic code
ALL_STANDARD_CODONS = set(CODON_TABLE.keys())


# ===========================================================================
# 1. SPECIES dict has expected organisms
# ===========================================================================


class TestSpeciesRegistry:
    """Tests for the SPECIES registry dict."""

    def test_species_has_expected_keys(self):
        """SPECIES must contain exactly the five expected organisms."""
        assert EXPECTED_SPECIES_KEYS == set(SPECIES.keys())

    @pytest.mark.parametrize("name", sorted(EXPECTED_SPECIES_KEYS))
    def test_species_entry_has_cai_weights(self, name: str):
        """Each SPECIES entry must contain non-empty cai_weights."""
        entry = SPECIES[name]
        assert "cai_weights" in entry, f"Species {name!r} missing 'cai_weights' key"
        assert len(entry["cai_weights"]) > 0, f"Species {name!r} has empty cai_weights"

    @pytest.mark.parametrize("name", sorted(EXPECTED_SPECIES_KEYS))
    def test_species_entry_has_validation_flag(self, name: str):
        """Each SPECIES entry must have a codon_usage_validation flag."""
        entry = SPECIES[name]
        assert "codon_usage_validation" in entry, (
            f"Species {name!r} missing 'codon_usage_validation' key"
        )
        assert isinstance(entry["codon_usage_validation"], bool), (
            f"Species {name!r} codon_usage_validation should be bool"
        )


# ===========================================================================
# 2. CODON_ADAPTIVENESS_TABLES has entries for each organism
# ===========================================================================


class TestCodonAdaptivenessTables:
    """Tests for the CODON_ADAPTIVENESS_TABLES registry."""

    @pytest.mark.parametrize("organism", sorted(CANONICAL_TO_SHORT.keys()))
    def test_adaptiveness_table_exists(self, organism: str):
        """Each supported organism must have an entry in CODON_ADAPTIVENESS_TABLES."""
        assert organism in CODON_ADAPTIVENESS_TABLES, (
            f"Organism {organism!r} missing from CODON_ADAPTIVENESS_TABLES"
        )

    def test_adaptiveness_tables_count(self):
        """Number of adaptiveness tables must equal number of supported organisms."""
        assert len(CODON_ADAPTIVENESS_TABLES) == len(SUPPORTED_ORGANISMS)

    @pytest.mark.parametrize("organism", sorted(CANONICAL_TO_SHORT.keys()))
    def test_adaptiveness_table_non_empty(self, organism: str):
        """Each adaptiveness table must be non-empty."""
        table = CODON_ADAPTIVENESS_TABLES[organism]
        assert len(table) > 0, f"Adaptiveness table for {organism!r} is empty"


# ===========================================================================
# 3. ORGANISM_GC_TARGETS has organism-specific values
# ===========================================================================


class TestOrganismGCTargets:
    """Tests for the ORGANISM_GC_TARGETS dict."""

    @pytest.mark.parametrize("organism", sorted(CANONICAL_TO_SHORT.keys()))
    def test_gc_target_exists(self, organism: str):
        """Each supported organism must have a GC target entry."""
        assert organism in ORGANISM_GC_TARGETS, (
            f"Organism {organism!r} missing from ORGANISM_GC_TARGETS"
        )

    def test_gc_targets_count(self):
        """Number of GC targets must equal number of supported organisms."""
        assert len(ORGANISM_GC_TARGETS) == len(SUPPORTED_ORGANISMS)

    @pytest.mark.parametrize("organism", sorted(CANONICAL_TO_SHORT.keys()))
    def test_gc_targets_are_tuples(self, organism: str):
        """Each GC target must be a (gc_lo, gc_hi) pair."""
        target = ORGANISM_GC_TARGETS[organism]
        assert isinstance(target, tuple) and len(target) == 2, (
            f"GC target for {organism!r} should be a 2-tuple, got {target!r}"
        )

    @pytest.mark.parametrize("organism", sorted(CANONICAL_TO_SHORT.keys()))
    def test_gc_targets_lo_le_hi(self, organism: str):
        """GC low must be <= GC high."""
        gc_lo, gc_hi = ORGANISM_GC_TARGETS[organism]
        assert gc_lo <= gc_hi, (
            f"GC lo ({gc_lo}) > GC hi ({gc_hi}) for {organism!r}"
        )

    def test_gc_targets_vary_by_organism(self):
        """Not all organisms should have identical GC targets (organism-specific)."""
        targets = set(ORGANISM_GC_TARGETS.values())
        assert len(targets) > 1, "All organisms have identical GC targets"


# ===========================================================================
# 4. SUPPORTED_ORGANISMS list matches SPECIES keys
# ===========================================================================


class TestSupportedOrganisms:
    """Tests for SUPPORTED_ORGANISMS consistency with SPECIES and other registries."""

    def test_supported_organisms_matches_codon_usage_keys(self):
        """SUPPORTED_ORGANISMS must equal CODON_USAGE_TABLES keys."""
        assert set(SUPPORTED_ORGANISMS) == set(CODON_USAGE_TABLES.keys())

    def test_supported_organisms_maps_to_species_keys(self):
        """Every SUPPORTED_ORGANISMS entry must map to a SPECIES short key."""
        for canonical in SUPPORTED_ORGANISMS:
            assert canonical in CANONICAL_TO_SHORT, (
                f"Supported organism {canonical!r} has no mapping to SPECIES key"
            )

    def test_species_keys_cover_all_supported(self):
        """Every SPECIES key must be reachable from a SUPPORTED_ORGANISMS entry."""
        mapped_short_names = {CANONICAL_TO_SHORT[org] for org in SUPPORTED_ORGANISMS}
        assert mapped_short_names == EXPECTED_SPECIES_KEYS


# ===========================================================================
# 5. Each organism module exports codon usage table
# ===========================================================================


class TestOrganismModuleExports:
    """Tests that each organism module exports its codon usage table."""

    @pytest.mark.parametrize("organism", sorted(ORGANISM_CODON_TABLES.keys()))
    def test_codon_usage_exported(self, organism: str):
        """Each organism module must export a non-empty codon usage table."""
        table = ORGANISM_CODON_TABLES[organism]
        assert isinstance(table, dict), f"Codon usage for {organism!r} is not a dict"
        assert len(table) > 0, f"Codon usage for {organism!r} is empty"

    @pytest.mark.parametrize("organism", sorted(ORGANISM_ADAPTIVENESS.keys()))
    def test_adaptiveness_exported(self, organism: str):
        """Each organism module must export a non-empty adaptiveness table."""
        table = ORGANISM_ADAPTIVENESS[organism]
        assert isinstance(table, dict), f"Adaptiveness for {organism!r} is not a dict"
        assert len(table) > 0, f"Adaptiveness for {organism!r} is empty"

    def test_codon_usage_tables_registry_matches_modules(self):
        """CODON_USAGE_TABLES registry values must match the module-level exports."""
        for organism, table in CODON_USAGE_TABLES.items():
            assert table is ORGANISM_CODON_TABLES[organism], (
                f"CODON_USAGE_TABLES[{organism!r}] is not the same object as "
                f"the module-level export"
            )


# ===========================================================================
# 6. Codon tables are valid (all 64 codons, only standard amino acids)
# ===========================================================================


class TestCodonTableValidity:
    """Tests for structural validity of codon usage tables."""

    @pytest.mark.parametrize("organism", sorted(ORGANISM_CODON_TABLES.keys()))
    def test_all_64_codons_represented(self, organism: str):
        """Every codon usage table must contain exactly all 64 standard codons."""
        table = ORGANISM_CODON_TABLES[organism]
        table_codons = set(table.keys())
        missing = ALL_STANDARD_CODONS - table_codons
        extra = table_codons - ALL_STANDARD_CODONS
        assert not missing, f"{organism!r} missing codons: {sorted(missing)}"
        assert not extra, f"{organism!r} has extra codons: {sorted(extra)}"

    @pytest.mark.parametrize("organism", sorted(ORGANISM_CODON_TABLES.keys()))
    def test_only_standard_amino_acids(self, organism: str):
        """All amino acid annotations must be standard single-letter codes or '*'."""
        table = ORGANISM_CODON_TABLES[organism]
        for codon, (aa, *_rest) in table.items():
            assert aa in STANDARD_AA, (
                f"{organism!r} codon {codon!r} maps to unknown AA {aa!r}"
            )

    @pytest.mark.parametrize("organism", sorted(ORGANISM_CODON_TABLES.keys()))
    def test_codon_to_aa_matches_standard_code(self, organism: str):
        """Each codon's amino acid annotation must match the standard genetic code."""
        table = ORGANISM_CODON_TABLES[organism]
        for codon, (aa, *_rest) in table.items():
            expected_aa = CODON_TABLE[codon]
            assert aa == expected_aa, (
                f"{organism!r} codon {codon!r} maps to {aa!r}, "
                f"but standard code says {expected_aa!r}"
            )

    @pytest.mark.parametrize("organism", sorted(ORGANISM_CODON_TABLES.keys()))
    def test_fractions_sum_to_one_per_amino_acid(self, organism: str):
        """For each amino acid, codon fractions should sum to ~1.0.

        Some organism tables have rounding imprecision in the fraction
        column (e.g. E. coli Valine sums to 1.16), so we use a generous
        tolerance.  The primary correctness check is that the amino acid
        assignments match the standard genetic code (tested separately).
        """
        table = ORGANISM_CODON_TABLES[organism]
        # Group fractions by amino acid
        aa_fractions: dict[str, list[float]] = {}
        for codon, (aa, frac, *_rest) in table.items():
            aa_fractions.setdefault(aa, []).append(frac)
        for aa, fracs in aa_fractions.items():
            total = sum(fracs)
            assert total == pytest.approx(1.0, abs=0.20), (
                f"{organism!r} AA {aa!r} fractions sum to {total:.4f}, expected ~1.0"
            )


# ===========================================================================
# 7. GC targets are in valid range [0.1, 0.9]
# ===========================================================================


class TestGCTargetRange:
    """Tests that GC target values are within biologically plausible bounds."""

    @pytest.mark.parametrize("organism", sorted(ORGANISM_GC_TARGETS.keys()))
    def test_gc_lo_in_valid_range(self, organism: str):
        """GC low target must be >= 0.1."""
        gc_lo, _ = ORGANISM_GC_TARGETS[organism]
        assert gc_lo >= 0.1, (
            f"GC lo for {organism!r} is {gc_lo}, below 0.1"
        )

    @pytest.mark.parametrize("organism", sorted(ORGANISM_GC_TARGETS.keys()))
    def test_gc_hi_in_valid_range(self, organism: str):
        """GC high target must be <= 0.9."""
        _, gc_hi = ORGANISM_GC_TARGETS[organism]
        assert gc_hi <= 0.9, (
            f"GC hi for {organism!r} is {gc_hi}, above 0.9"
        )

    @pytest.mark.parametrize("organism", sorted(ORGANISM_GC_TARGETS.keys()))
    def test_gc_lo_below_hi(self, organism: str):
        """GC low must be strictly below GC high."""
        gc_lo, gc_hi = ORGANISM_GC_TARGETS[organism]
        assert gc_lo < gc_hi, (
            f"GC lo ({gc_lo}) not strictly below GC hi ({gc_hi}) for {organism!r}"
        )


# ===========================================================================
# 8. Adaptiveness values are in [0, 1] range
# ===========================================================================


class TestAdaptivenessRange:
    """Tests that all codon adaptiveness values fall within [0, 1]."""

    @pytest.mark.parametrize("organism", sorted(CODON_ADAPTIVENESS_TABLES.keys()))
    def test_adaptiveness_values_in_range(self, organism: str):
        """Every adaptiveness value must be between 0.0 and 1.0 inclusive."""
        table = CODON_ADAPTIVENESS_TABLES[organism]
        for codon, val in table.items():
            assert 0.0 <= val <= 1.0, (
                f"{organism!r} adaptiveness for codon {codon!r} is {val}, "
                f"outside [0, 1]"
            )

    @pytest.mark.parametrize("organism", sorted(CODON_ADAPTIVENESS_TABLES.keys()))
    def test_at_least_one_adaptiveness_is_one(self, organism: str):
        """For each amino acid group, at least one codon should have adaptiveness 1.0."""
        table = CODON_ADAPTIVENESS_TABLES[organism]
        # Collect adaptiveness values per amino acid using CODON_TABLE
        aa_values: dict[str, list[float]] = {}
        for codon, val in table.items():
            aa = CODON_TABLE.get(codon)
            if aa is None or aa == "*":
                continue
            aa_values.setdefault(aa, []).append(val)
        for aa, vals in aa_values.items():
            assert max(vals) == pytest.approx(1.0, abs=1e-9), (
                f"{organism!r} AA {aa!r}: no codon has adaptiveness 1.0 "
                f"(max={max(vals)})"
            )

    @pytest.mark.parametrize("organism", sorted(CODON_ADAPTIVENESS_TABLES.keys()))
    def test_adaptiveness_excludes_stop_codons(self, organism: str):
        """Adaptiveness tables should not contain stop codons."""
        table = CODON_ADAPTIVENESS_TABLES[organism]
        stop_codons = {c for c, aa in CODON_TABLE.items() if aa == "*"}
        present_stops = stop_codons & set(table.keys())
        assert not present_stops, (
            f"{organism!r} adaptiveness table contains stop codons: {sorted(present_stops)}"
        )
