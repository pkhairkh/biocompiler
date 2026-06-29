"""Tests for the expanded DNA Chisel compatibility layer constraint mappings.

Covers:
1. CONSTRAINT_MAPPING registry — 10+ constraint types registered
2. build_constraint_spec — individual constraint building
3. translate_biocompiler_constraints — high-level constraint translation
4. Individual builder functions — enforce_translation, enforce_gc_content, etc.
5. Graceful degradation when DNA Chisel is not installed
6. New named constants for extended mappings

These tests work whether or not DNA Chisel is installed. When it is absent,
constraint builders return None, which we verify explicitly.
"""

from __future__ import annotations

import pytest
from unittest.mock import patch, MagicMock

from biocompiler.infrastructure.dna_chisel_compat import (
    CONSTRAINT_MAPPING,
    build_constraint_spec,
    translate_biocompiler_constraints,
    is_dna_chisel_available,
    _DNA_CHISEL_AVAILABLE,
    _DNA_CHISEL_CONSTRAINTS,
    _CONSTRAINT_BUILDERS,
    # New constants
    DEFAULT_HAIRPIN_STEM_SIZE,
    DEFAULT_HAIRPIN_BOOST,
    DEFAULT_KMER_UNIQUIFY_SIZE,
    DEFAULT_BACTERIAL_PROMOTER_LENGTH,
    GC_ENFORCEMENT_WINDOW,
    # Builder functions
    enforce_translation,
    enforce_gc_content,
    enforce_gc_content_local,
    avoid_pattern,
    avoid_bacterial_promoter,
    enforce_start_codon,
    enforce_stop_codon,
    uniquify_all_kmers,
    avoid_changes,
    enforce_sequence,
    avoid_hairpins,
)


# ---------------------------------------------------------------------------
# 1. CONSTRAINT_MAPPING registry
# ---------------------------------------------------------------------------

class TestConstraintMappingRegistry:
    """Tests for the CONSTRAINT_MAPPING registry."""

    def test_has_10_plus_entries(self):
        """CONSTRAINT_MAPPING has at least 10 entries."""
        assert len(CONSTRAINT_MAPPING) >= 10

    def test_contains_core_constraints(self):
        """Core constraints (EnforceTranslation, EnforceGCContent, AvoidPattern) are present."""
        assert "EnforceTranslation" in CONSTRAINT_MAPPING
        assert "EnforceGCContent" in CONSTRAINT_MAPPING
        assert "AvoidPattern" in CONSTRAINT_MAPPING

    def test_contains_local_gc(self):
        """EnforceGCContentLocal is registered."""
        assert "EnforceGCContentLocal" in CONSTRAINT_MAPPING

    def test_contains_bacterial_promoter(self):
        """AvoidBacterialPromoter is registered."""
        assert "AvoidBacterialPromoter" in CONSTRAINT_MAPPING

    def test_contains_start_codon(self):
        """EnforceStartCodon is registered."""
        assert "EnforceStartCodon" in CONSTRAINT_MAPPING

    def test_contains_stop_codon(self):
        """EnforceStopCodon is registered."""
        assert "EnforceStopCodon" in CONSTRAINT_MAPPING

    def test_contains_uniquify_kmers(self):
        """UniquifyAllKmers is registered."""
        assert "UniquifyAllKmers" in CONSTRAINT_MAPPING

    def test_contains_avoid_changes(self):
        """AvoidChanges is registered."""
        assert "AvoidChanges" in CONSTRAINT_MAPPING

    def test_contains_enforce_sequence(self):
        """EnforceSequence is registered."""
        assert "EnforceSequence" in CONSTRAINT_MAPPING

    def test_contains_avoid_hairpins(self):
        """AvoidHairpins is registered."""
        assert "AvoidHairpins" in CONSTRAINT_MAPPING

    def test_all_values_are_strings(self):
        """All values in CONSTRAINT_MAPPING are strings (builder function names)."""
        for key, val in CONSTRAINT_MAPPING.items():
            assert isinstance(val, str), f"Value for {key} is not a string"

    def test_all_builder_names_have_functions(self):
        """All builder function names in mapping have corresponding callables."""
        for key, builder_name in CONSTRAINT_MAPPING.items():
            assert builder_name in _CONSTRAINT_BUILDERS, (
                f"Builder {builder_name!r} for {key!r} not in _CONSTRAINT_BUILDERS"
            )


# ---------------------------------------------------------------------------
# 2. Named constants for extended mappings
# ---------------------------------------------------------------------------

class TestExtendedConstants:
    """Tests for new named constants added with the expansion."""

    def test_default_hairpin_stem_size_positive(self):
        """DEFAULT_HAIRPIN_STEM_SIZE is a positive integer."""
        assert isinstance(DEFAULT_HAIRPIN_STEM_SIZE, int)
        assert DEFAULT_HAIRPIN_STEM_SIZE > 0

    def test_default_hairpin_boost_positive(self):
        """DEFAULT_HAIRPIN_BOOST is a positive float."""
        assert isinstance(DEFAULT_HAIRPIN_BOOST, float)
        assert DEFAULT_HAIRPIN_BOOST > 0.0

    def test_default_kmer_uniquify_size_positive(self):
        """DEFAULT_KMER_UNIQUIFY_SIZE is a positive integer."""
        assert isinstance(DEFAULT_KMER_UNIQUIFY_SIZE, int)
        assert DEFAULT_KMER_UNIQUIFY_SIZE > 0

    def test_default_bacterial_promoter_length_positive(self):
        """DEFAULT_BACTERIAL_PROMOTER_LENGTH is a positive integer."""
        assert isinstance(DEFAULT_BACTERIAL_PROMOTER_LENGTH, int)
        assert DEFAULT_BACTERIAL_PROMOTER_LENGTH > 0

    def test_gc_enforcement_window_unchanged(self):
        """GC_ENFORCEMENT_WINDOW still exists and is positive."""
        assert isinstance(GC_ENFORCEMENT_WINDOW, int)
        assert GC_ENFORCEMENT_WINDOW > 0


# ---------------------------------------------------------------------------
# 3. build_constraint_spec
# ---------------------------------------------------------------------------

class TestBuildConstraintSpec:
    """Tests for the build_constraint_spec function."""

    @pytest.mark.skipif(
        not _DNA_CHISEL_AVAILABLE,
        reason="Requires DNA Chisel to be installed"
    )
    def test_returns_object_when_available(self):
        """Returns a DNA Chisel object when available."""
        result = build_constraint_spec("EnforceTranslation", protein="M")
        assert result is not None

    @pytest.mark.skipif(
        not _DNA_CHISEL_AVAILABLE,
        reason="Requires DNA Chisel to be installed"
    )
    def test_returns_none_for_unknown_type(self):
        """Returns None for unknown constraint type."""
        result = build_constraint_spec("NonExistentConstraint")
        assert result is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests behavior when DNA Chisel is NOT installed"
    )
    def test_returns_none_when_not_installed(self):
        """Returns None when DNA Chisel is not installed."""
        result = build_constraint_spec("EnforceTranslation", protein="M")
        assert result is None

    @pytest.mark.skipif(
        not _DNA_CHISEL_AVAILABLE,
        reason="Requires DNA Chisel to be installed"
    )
    def test_enforce_gc_content_builds(self):
        """EnforceGCContent can be built."""
        result = build_constraint_spec("EnforceGCContent", gc_lo=0.3, gc_hi=0.7)
        assert result is not None

    @pytest.mark.skipif(
        not _DNA_CHISEL_AVAILABLE,
        reason="Requires DNA Chisel to be installed"
    )
    def test_avoid_pattern_builds(self):
        """AvoidPattern can be built."""
        result = build_constraint_spec("AvoidPattern", pattern="GAATTC")
        assert result is not None


# ---------------------------------------------------------------------------
# 4. Individual builder functions (work without DNA Chisel installed)
# ---------------------------------------------------------------------------

class TestBuilderFunctionsWithoutChisel:
    """Tests for individual builder functions when DNA Chisel is not installed."""

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_enforce_translation_returns_none(self):
        """enforce_translation returns None when DNA Chisel unavailable."""
        assert enforce_translation(protein="M") is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_enforce_gc_content_returns_none(self):
        """enforce_gc_content returns None when DNA Chisel unavailable."""
        assert enforce_gc_content(gc_lo=0.3, gc_hi=0.7) is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_enforce_gc_content_local_returns_none(self):
        """enforce_gc_content_local returns None when DNA Chisel unavailable."""
        assert enforce_gc_content_local(gc_lo=0.3, gc_hi=0.7, window=50) is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_avoid_pattern_returns_none(self):
        """avoid_pattern returns None when DNA Chisel unavailable."""
        assert avoid_pattern(pattern="GAATTC") is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_avoid_bacterial_promoter_returns_none(self):
        """avoid_bacterial_promoter returns None when DNA Chisel unavailable."""
        assert avoid_bacterial_promoter() is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_enforce_start_codon_returns_none(self):
        """enforce_start_codon returns None when DNA Chisel unavailable."""
        assert enforce_start_codon() is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_enforce_stop_codon_returns_none(self):
        """enforce_stop_codon returns None when DNA Chisel unavailable."""
        assert enforce_stop_codon() is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_uniquify_all_kmers_returns_none(self):
        """uniquify_all_kmers returns None when DNA Chisel unavailable."""
        assert uniquify_all_kmers(kmer_size=9) is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_avoid_changes_returns_none(self):
        """avoid_changes returns None when DNA Chisel unavailable."""
        assert avoid_changes(zone="0-10") is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_enforce_sequence_returns_none(self):
        """enforce_sequence returns None when DNA Chisel unavailable."""
        assert enforce_sequence(sequence="ATGC") is None

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests graceful None returns when DNA Chisel is NOT installed"
    )
    def test_avoid_hairpins_returns_none(self):
        """avoid_hairpins returns None when DNA Chisel unavailable."""
        assert avoid_hairpins(stem_size=15) is None


# ---------------------------------------------------------------------------
# 5. Builder functions with DNA Chisel mocked as available
# ---------------------------------------------------------------------------

class TestBuilderFunctionsWithChiselMocked:
    """Tests for builder functions with DNA Chisel mocked as available."""

    def test_enforce_translation_builds_with_mock(self):
        """enforce_translation builds EnforceTranslation constraint."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"EnforceTranslation": mock_cls}, clear=False):
            result = enforce_translation(protein="MVSKGE")
            mock_cls.assert_called_once_with(translation="MVSKGE")
            assert result is not None

    def test_enforce_gc_content_builds_with_mock(self):
        """enforce_gc_content builds EnforceGCContent constraint."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"EnforceGCContent": mock_cls}, clear=False):
            result = enforce_gc_content(gc_lo=0.4, gc_hi=0.6)
            mock_cls.assert_called_once_with(mini=0.4, maxi=0.6)
            assert result is not None

    def test_enforce_gc_content_local_with_window(self):
        """enforce_gc_content_local adds window parameter."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"EnforceGCContent": mock_cls}, clear=False):
            result = enforce_gc_content_local(gc_lo=0.3, gc_hi=0.7, window=75)
            mock_cls.assert_called_once_with(mini=0.3, maxi=0.7, window=75)

    def test_avoid_pattern_builds_with_mock(self):
        """avoid_pattern builds AvoidPattern constraint."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"AvoidPattern": mock_cls}, clear=False):
            result = avoid_pattern(pattern="GAATTC")
            mock_cls.assert_called_once_with("GAATTC")

    def test_avoid_bacterial_promoter_with_mock(self):
        """avoid_bacterial_promoter builds with length parameter."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"AvoidBacterialPromoter": mock_cls}, clear=False):
            result = avoid_bacterial_promoter(length=40)
            mock_cls.assert_called_once_with(length=40)

    def test_enforce_start_codon_with_mock(self):
        """enforce_start_codon builds with start_codon parameter."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"EnforceStartCodon": mock_cls}, clear=False):
            result = enforce_start_codon(start_codon="ATG")
            mock_cls.assert_called_once_with(start_codon="ATG")

    def test_enforce_stop_codon_with_mock(self):
        """enforce_stop_codon builds with location parameter."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"EnforceStopCodon": mock_cls}, clear=False):
            result = enforce_stop_codon(location="end")
            mock_cls.assert_called_once_with(location="end")

    def test_uniquify_all_kmers_with_mock(self):
        """uniquify_all_kmers builds with k parameter (maps from kmer_size)."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"UniquifyAllKmers": mock_cls}, clear=False):
            result = uniquify_all_kmers(kmer_size=12)
            # DNA Chisel uses k=, not kmer_size=
            mock_cls.assert_called_once_with(k=12)

    def test_avoid_changes_with_zone(self):
        """avoid_changes builds with zone parameter (mapped to location)."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"AvoidChanges": mock_cls}, clear=False):
            result = avoid_changes(zone="0-50")
            # Zone "0-50" is parsed to location=(0, 50) for DNA Chisel
            mock_cls.assert_called_once_with(location=(0, 50))

    def test_avoid_changes_without_zone(self):
        """avoid_changes builds without zone parameter."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"AvoidChanges": mock_cls}, clear=False):
            result = avoid_changes()
            mock_cls.assert_called_once_with()

    def test_enforce_sequence_with_mock(self):
        """enforce_sequence builds with sequence parameter."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"EnforceSequence": mock_cls}, clear=False):
            result = enforce_sequence(sequence="ATGCGT")
            # EnforceSequence passes sequence as positional arg
            mock_cls.assert_called_once_with("ATGCGT")

    def test_avoid_hairpins_with_mock(self):
        """avoid_hairpins builds with stem_size and boost parameters."""
        mock_cls = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {"AvoidHairpins": mock_cls}, clear=False):
            result = avoid_hairpins(stem_size=20, boost=1.5)
            mock_cls.assert_called_once_with(stem_size=20, boost=1.5)

    def test_builder_returns_none_when_class_missing(self):
        """Builder returns None when the constraint class is not in _DNA_CHISEL_CONSTRAINTS."""
        # Remove a class temporarily
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {}, clear=True):
            result = enforce_translation(protein="M")
            assert result is None


# ---------------------------------------------------------------------------
# 6. translate_biocompiler_constraints
# ---------------------------------------------------------------------------

class TestTranslateBiocompilerConstraints:
    """Tests for the translate_biocompiler_constraints high-level function."""

    @pytest.mark.skipif(
        not _DNA_CHISEL_AVAILABLE,
        reason="Requires DNA Chisel to be installed"
    )
    def test_basic_translation_with_protein(self):
        """Basic translation with just protein and GC bounds."""
        specs = translate_biocompiler_constraints(
            protein="MVSKGE",
            gc_lo=0.3,
            gc_hi=0.7,
        )
        assert isinstance(specs, list)
        # Should have at least EnforceTranslation + EnforceGCContent
        assert len(specs) >= 2

    @pytest.mark.skipif(
        _DNA_CHISEL_AVAILABLE,
        reason="Tests empty list when DNA Chisel is NOT installed"
    )
    def test_returns_empty_when_not_installed(self):
        """Returns empty list when DNA Chisel is not installed."""
        specs = translate_biocompiler_constraints(protein="M")
        assert specs == []

    @pytest.mark.skipif(
        not _DNA_CHISEL_AVAILABLE,
        reason="Requires DNA Chisel to be installed"
    )
    def test_with_restriction_enzymes(self):
        """Adding restriction enzymes adds AvoidPattern specs."""
        specs_no_enz = translate_biocompiler_constraints(protein="MVSKGE", enforce_start=False, enforce_stop=False)
        specs_with_enz = translate_biocompiler_constraints(
            protein="MVSKGE",
            restriction_enzymes=["EcoRI"],
            enforce_start=False,
            enforce_stop=False,
        )
        assert len(specs_with_enz) > len(specs_no_enz)

    @pytest.mark.skipif(
        not _DNA_CHISEL_AVAILABLE,
        reason="Requires DNA Chisel to be installed"
    )
    def test_with_local_gc_window(self):
        """local_gc_window adds an extra EnforceGCContent with window."""
        specs_no_window = translate_biocompiler_constraints(
            protein="MVSKGE",
            enforce_start=False,
            enforce_stop=False,
        )
        specs_with_window = translate_biocompiler_constraints(
            protein="MVSKGE",
            local_gc_window=50,
            enforce_start=False,
            enforce_stop=False,
        )
        assert len(specs_with_window) > len(specs_no_window)

    @pytest.mark.skipif(
        not _DNA_CHISEL_AVAILABLE,
        reason="Requires DNA Chisel to be installed"
    )
    def test_all_flags_enabled(self):
        """Enabling all flags produces maximum constraints."""
        specs = translate_biocompiler_constraints(
            protein="MVSKGE",
            gc_lo=0.3,
            gc_hi=0.7,
            restriction_enzymes=["EcoRI"],
            local_gc_window=50,
            avoid_bacterial_promoters=True,
            enforce_start=True,
            enforce_stop=True,
            uniquify_kmers=None,  # Skip to avoid API version issues
            preserve_zones=["0-10"],
            enforce_sequence_str="ATGC",
            avoid_hairpins_flag=True,
        )
        assert isinstance(specs, list)
        # Should have many specs (exact count depends on which constraints are available)
        assert len(specs) >= 2

    def test_no_protein_skips_translation_constraint(self):
        """Empty protein skips EnforceTranslation."""
        # We mock as available to test the logic
        mock_et = MagicMock()
        mock_egc = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {
            "EnforceTranslation": mock_et,
            "EnforceGCContent": mock_egc,
            "AvoidPattern": MagicMock(),
        }, clear=True):
            with patch("biocompiler.infrastructure.dna_chisel_compat._DNA_CHISEL_AVAILABLE", True):
                specs = translate_biocompiler_constraints(
                    protein="",
                    gc_lo=0.3,
                    gc_hi=0.7,
                    enforce_start=False,
                    enforce_stop=False,
                )
                # EnforceTranslation should NOT have been called
                mock_et.assert_not_called()
                # But EnforceGCContent should have been called
                mock_egc.assert_called()

    def test_preserve_zones_multiple(self):
        """Multiple preserve zones create multiple AvoidChanges specs."""
        mock_ac = MagicMock()
        mock_egc = MagicMock()
        with patch.dict(_DNA_CHISEL_CONSTRAINTS, {
            "AvoidChanges": mock_ac,
            "EnforceGCContent": mock_egc,
        }, clear=True):
            with patch("biocompiler.infrastructure.dna_chisel_compat._DNA_CHISEL_AVAILABLE", True):
                specs = translate_biocompiler_constraints(
                    preserve_zones=["0-10", "20-30"],
                    enforce_start=False,
                    enforce_stop=False,
                )
                # AvoidChanges should have been called twice
                assert mock_ac.call_count == 2


# ---------------------------------------------------------------------------
# 7. Constraint mapping coverage
# ---------------------------------------------------------------------------

class TestConstraintMappingCoverage:
    """Verify that the expanded mapping covers 10+ distinct constraint types."""

    def test_mapping_has_all_11_entries(self):
        """CONSTRAINT_MAPPING has exactly 11 entries (3 original + 8 new)."""
        expected = {
            "EnforceTranslation",
            "EnforceGCContent",
            "EnforceGCContentLocal",
            "AvoidPattern",
            "AvoidBacterialPromoter",
            "EnforceStartCodon",
            "EnforceStopCodon",
            "UniquifyAllKmers",
            "AvoidChanges",
            "EnforceSequence",
            "AvoidHairpins",
        }
        assert set(CONSTRAINT_MAPPING.keys()) == expected

    def test_all_builders_registered(self):
        """All builder function names have corresponding callables."""
        expected_builders = {
            "enforce_translation",
            "enforce_gc_content",
            "enforce_gc_content_local",
            "avoid_pattern",
            "avoid_bacterial_promoter",
            "enforce_start_codon",
            "enforce_stop_codon",
            "uniquify_all_kmers",
            "avoid_changes",
            "enforce_sequence",
            "avoid_hairpins",
        }
        for name in expected_builders:
            assert name in _CONSTRAINT_BUILDERS, f"Missing builder: {name}"
            assert callable(_CONSTRAINT_BUILDERS[name]), f"Builder {name} is not callable"

    def test_distinct_builder_names(self):
        """All builder function names in CONSTRAINT_MAPPING are distinct."""
        values = list(CONSTRAINT_MAPPING.values())
        assert len(values) == len(set(values)), "Duplicate builder names found"
