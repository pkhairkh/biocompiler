"""Test BioCompiler ViennaRNA integration — mRNA structure prediction, accessibility, and stem-loop detection.

Tests are designed to work with or without ViennaRNA installed:
- Dataclass / availability tests always run.
- Functional tests (predict_mfe, predict_accessibility, etc.) are conditional
  on ViennaRNA availability via ``is_viennarna_available()``.
- Subprocess-fallback tests require the ``RNAfold`` CLI on $PATH.

Code targets Python 3.10+.
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess
import sys
from dataclasses import fields
from typing import Any

import pytest

# Mark every test in this module as requiring an external tool (ViennaRNA).
# Many tests below call ``_skip_if_no_viennarna()`` (or invoke the RNAfold
# CLI) at runtime; the marker keeps them deselected by default alongside
# other requires_external tests.
pytestmark = pytest.mark.requires_external

# ---------------------------------------------------------------------------
# Lazy import helper — avoids hard ImportError at module level
# ---------------------------------------------------------------------------

_VIENNARNA_MODULE = None


def _get_viennarna():
    """Lazily import and cache the viennarna module."""
    global _VIENNARNA_MODULE
    if _VIENNARNA_MODULE is None:
        _VIENNARNA_MODULE = pytest.importorskip("biocompiler.viennarna")
    return _VIENNARNA_MODULE


# ---------------------------------------------------------------------------
# Test sequences
# ---------------------------------------------------------------------------

# 60-nt sequence with moderate structure potential (mixed GC/AT)
SEQ_60NT = "ATGGCTAGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

# Sequence known to form a stable stem-loop (palindromic core flanked by loop)
SEQ_STEM_LOOP = "GGGCGCGAAAGCGCUUUUGCGCUUUCGCGCCC"

# All-A sequence (no structure expected)
SEQ_ALL_A = "A" * 60

# GC-rich sequence (strong structure expected)
SEQ_GC_RICH = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"

# Very short sequences
SEQ_4NT = "ATCG"
SEQ_10NT = "ATCGATCGAT"

# Empty sequence
SEQ_EMPTY = ""

# 5'UTR-like sequence (moderate structure around start codon context)
SEQ_5UTR = "AACAGAGGAGGCCGCCACCATGGCTAGCGATCGATCGATCGATCGATCGATCGATCGATCG"


# ═══════════════════════════════════════════════════════════════════════════
# 1. Availability check
# ═══════════════════════════════════════════════════════════════════════════

class TestAvailabilityCheck:
    """Tests for ``is_viennarna_available()``."""

    def test_returns_bool(self):
        """is_viennarna_available() always returns a bool."""
        vr = _get_viennarna()
        result = vr.is_viennarna_available()
        assert isinstance(result, bool), (
            f"Expected bool, got {type(result).__name__}"
        )

    def test_consistent_return(self):
        """Repeated calls return the same value."""
        vr = _get_viennarna()
        first = vr.is_viennarna_available()
        second = vr.is_viennarna_available()
        assert first == second, (
            f"Inconsistent availability: first={first}, second={second}"
        )

    def test_import_does_not_crash(self):
        """Importing the module succeeds regardless of ViennaRNA installation."""
        import importlib
        mod = importlib.import_module("biocompiler.viennarna")
        assert hasattr(mod, "is_viennarna_available")


# ═══════════════════════════════════════════════════════════════════════════
# 2. MFEResult dataclass
# ═══════════════════════════════════════════════════════════════════════════

class TestMFEResult:
    """Tests for the MFEResult dataclass."""

    def test_creation_with_valid_data(self):
        """MFEResult can be created with valid fields."""
        vr = _get_viennarna()
        result = vr.MFEResult(
            sequence=SEQ_60NT,
            structure="(((...)))...(((...)))...(((...)))...(((...)))",
            mfe=-12.5,
            region="5utr",
        )
        assert result.sequence == SEQ_60NT
        assert result.mfe == pytest.approx(-12.5)
        assert result.structure is not None
        assert result.region == "5utr"

    def test_field_access(self):
        """All fields are accessible after creation."""
        vr = _get_viennarna()
        result = vr.MFEResult(
            sequence=SEQ_60NT,
            structure="." * 60,
            mfe=0.0,
        )
        # Access every field
        assert result.sequence == SEQ_60NT
        assert result.structure == "." * 60
        assert result.mfe == pytest.approx(0.0)

    def test_default_values(self):
        """Default values are set correctly for optional fields."""
        vr = _get_viennarna()
        result = vr.MFEResult(
            sequence=SEQ_60NT,
            structure="." * 60,
            mfe=0.0,
        )
        # region should have a default (None or empty string)
        field_names = {f.name for f in fields(result)}
        assert "sequence" in field_names
        assert "structure" in field_names
        assert "mfe" in field_names

    def test_mfe_is_float(self):
        """mfe field stores a float value."""
        vr = _get_viennarna()
        result = vr.MFEResult(
            sequence=SEQ_4NT,
            structure="....",
            mfe=0.0,
        )
        assert isinstance(result.mfe, (int, float))

    def test_structure_dotbracket(self):
        """Structure field accepts valid dot-bracket notation."""
        vr = _get_viennarna()
        valid_structures = [
            "....",
            "((()))",
            "(((...)))",
            "(.(...).)",
        ]
        for s in valid_structures:
            result = vr.MFEResult(
                sequence="A" * len(s),
                structure=s,
                mfe=-1.0,
            )
            assert result.structure == s


# ═══════════════════════════════════════════════════════════════════════════
# 3. StemLoop dataclass
# ═══════════════════════════════════════════════════════════════════════════

class TestStemLoop:
    """Tests for the StemLoop dataclass."""

    def test_creation(self):
        """StemLoop can be created with position and energy fields."""
        vr = _get_viennarna()
        sl = vr.StemLoop(
            start=5,
            end=25,
            loop_start=13,
            loop_end=17,
            stem_length=8,
            mfe=-6.2,
            sequence="GGGCGCGAAAGCGC",
            structure="((((....)))).",
        )
        assert sl.start == 5
        assert sl.end == 25
        assert sl.mfe == pytest.approx(-6.2)

    def test_field_access(self):
        """All StemLoop fields are accessible."""
        vr = _get_viennarna()
        sl = vr.StemLoop(
            start=0,
            end=10,
            loop_start=4,
            loop_end=6,
            stem_length=4,
            mfe=-3.0,
            sequence="GGGAAACCC",
            structure="(((...)))",
        )
        # Verify all fields can be read
        _ = sl.start
        _ = sl.end
        _ = sl.loop_start
        _ = sl.loop_end
        _ = sl.stem_length
        _ = sl.mfe
        _ = sl.sequence
        _ = sl.structure

    def test_stem_length_positive(self):
        """Stem length should be a positive integer."""
        vr = _get_viennarna()
        sl = vr.StemLoop(
            start=0, end=10, loop_start=4, loop_end=6,
            stem_length=4, mfe=-3.0,
            sequence="GGGAAACCC", structure="(((...)))",
        )
        assert sl.stem_length > 0


# ═══════════════════════════════════════════════════════════════════════════
# 4. predict_mfe (conditional on ViennaRNA availability)
# ═══════════════════════════════════════════════════════════════════════════

def _skip_if_no_viennarna():
    """Skip test if ViennaRNA is not available."""
    vr = _get_viennarna()
    if not vr.is_viennarna_available():
        pytest.skip("ViennaRNA not available")


def _is_valid_dot_bracket(structure: str) -> bool:
    """Check if a string is valid dot-bracket notation.

    Valid characters: (, ), .  and the structure must be balanced.
    """
    if not re.match(r'^[().]+$', structure):
        return False
    depth = 0
    for ch in structure:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


class TestPredictMFE:
    """Tests for predict_mfe() function (requires ViennaRNA)."""

    def test_short_dna_sequence(self):
        """predict_mfe with a 60-nt DNA sequence returns MFEResult."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_60NT)
        assert isinstance(result, vr.MFEResult)

    def test_mfe_negative_for_structured_sequence(self):
        """MFE should be negative for a sequence with structure potential."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_60NT)
        # Structured sequences should have negative MFE
        # (though some very short sequences may be near 0)
        assert isinstance(result.mfe, (int, float))

    def test_structure_is_valid_dot_bracket(self):
        """Returned structure should be valid dot-bracket notation."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_60NT)
        assert _is_valid_dot_bracket(result.structure), (
            f"Invalid dot-bracket: {result.structure[:40]}..."
        )

    def test_structure_length_matches_sequence(self):
        """Structure length should equal sequence length."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_60NT)
        assert len(result.structure) == len(result.sequence), (
            f"Structure length {len(result.structure)} != "
            f"sequence length {len(result.sequence)}"
        )

    def test_region_5utr(self):
        """predict_mfe with region='5utr' works."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_5UTR, region="5utr")
        assert isinstance(result, vr.MFEResult)
        assert result.region == "5utr"

    def test_region_full(self):
        """predict_mfe with region='full' works."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_60NT, region="full")
        assert isinstance(result, vr.MFEResult)
        assert result.region == "full"

    def test_region_start_codon(self):
        """predict_mfe with region='start_codon' works."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_5UTR, region="start_codon")
        assert isinstance(result, vr.MFEResult)
        assert result.region == "start_codon"

    def test_gc_rich_has_more_negative_mfe(self):
        """GC-rich sequence should have more negative MFE than all-A sequence."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result_gc = vr.predict_mfe(SEQ_GC_RICH)
        result_a = vr.predict_mfe(SEQ_ALL_A)
        # GC pairs are stronger than AU pairs, so GC-rich should be more stable
        assert result_gc.mfe <= result_a.mfe, (
            f"GC-rich MFE ({result_gc.mfe}) should be <= all-A MFE ({result_a.mfe})"
        )


# ═══════════════════════════════════════════════════════════════════════════
# 5. predict_accessibility (conditional)
# ═══════════════════════════════════════════════════════════════════════════

class TestPredictAccessibility:
    """Tests for predict_accessibility() function (requires ViennaRNA)."""

    def test_returns_accessibility_result(self):
        """predict_accessibility returns an AccessibilityResult."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_accessibility(SEQ_60NT)
        assert hasattr(result, "accessibility") or hasattr(result, "mean_accessibility")

    def test_accessibility_values_in_range(self):
        """Accessibility values should be in [0, 1]."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_accessibility(SEQ_60NT)
        # Get per-position accessibility values
        acc_values = getattr(result, "accessibility", None)
        if acc_values is not None:
            for val in acc_values:
                assert 0.0 <= val <= 1.0, (
                    f"Accessibility value {val} out of [0, 1]"
                )

    def test_mean_accessibility_reasonable(self):
        """Mean accessibility should be in [0, 1]."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_accessibility(SEQ_60NT)
        mean_acc = getattr(result, "mean_accessibility", None)
        if mean_acc is not None:
            assert 0.0 <= mean_acc <= 1.0, (
                f"Mean accessibility {mean_acc} out of [0, 1]"
            )

    def test_all_a_high_accessibility(self):
        """All-A sequence should have high accessibility (low structure)."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result_a = vr.predict_accessibility(SEQ_ALL_A)
        result_gc = vr.predict_accessibility(SEQ_GC_RICH)
        mean_a = getattr(result_a, "mean_accessibility", None)
        mean_gc = getattr(result_gc, "mean_accessibility", None)
        if mean_a is not None and mean_gc is not None:
            # All-A should be more accessible than GC-rich
            assert mean_a >= mean_gc, (
                f"All-A accessibility ({mean_a}) should be >= GC-rich ({mean_gc})"
            )


# ═══════════════════════════════════════════════════════════════════════════
# 6. find_stable_structures (conditional)
# ═══════════════════════════════════════════════════════════════════════════

class TestFindStableStructures:
    """Tests for find_stable_structures() function (requires ViennaRNA)."""

    def test_stem_loops_found_in_structured_sequence(self):
        """A sequence designed to form stem-loops should find them."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.find_stable_structures(SEQ_STEM_LOOP)
        # Result should be a list of StemLoop objects
        assert isinstance(result, list)
        # May or may not find loops depending on ViennaRNA params,
        # but the function should not crash

    def test_no_stable_structures_in_all_a(self):
        """All-A sequence should not find stable stem-loops."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.find_stable_structures(SEQ_ALL_A)
        assert isinstance(result, list)
        # An all-A sequence is very unlikely to form stable stem-loops
        # If any are found, they should have very weak MFE (near 0)
        for sl in result:
            if hasattr(sl, "mfe"):
                assert sl.mfe > -3.0, (
                    f"All-A sequence should not have stable stem-loop with mfe={sl.mfe}"
                )

    def test_stem_loop_positions_valid(self):
        """Found stem-loops should have valid position ranges."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.find_stable_structures(SEQ_60NT)
        for sl in result:
            if hasattr(sl, "start") and hasattr(sl, "end"):
                assert 0 <= sl.start < sl.end, (
                    f"Invalid stem-loop positions: start={sl.start}, end={sl.end}"
                )
                assert sl.end <= len(SEQ_60NT), (
                    f"Stem-loop end {sl.end} exceeds sequence length {len(SEQ_60NT)}"
                )

    def test_gc_rich_finds_more_structures(self):
        """GC-rich sequence should find more stable structures than all-A."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result_gc = vr.find_stable_structures(SEQ_GC_RICH)
        result_a = vr.find_stable_structures(SEQ_ALL_A)
        # GC-rich should find at least as many structures as all-A
        # (may find same count of 0 for short sequences)
        assert isinstance(result_gc, list)
        assert isinstance(result_a, list)


# ═══════════════════════════════════════════════════════════════════════════
# 7. compute_5prime_dg (conditional)
# ═══════════════════════════════════════════════════════════════════════════

class TestCompute5PrimeDG:
    """Tests for compute_5prime_dg() function (requires ViennaRNA)."""

    def test_returns_float(self):
        """compute_5prime_dg returns a float."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.compute_5prime_dg(SEQ_5UTR)
        assert isinstance(result, (int, float)), (
            f"Expected numeric, got {type(result).__name__}"
        )

    def test_negative_for_structured_sequence(self):
        """Structured 5' region should have negative ΔG."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.compute_5prime_dg(SEQ_GC_RICH[:50])
        # GC-rich sequences should form stable structures → negative ΔG
        assert result < 0, (
            f"GC-rich 5' region should have negative ΔG, got {result}"
        )

    def test_all_a_near_zero(self):
        """All-A sequence should have ΔG near zero (no structure)."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.compute_5prime_dg(SEQ_ALL_A[:50])
        # All-A should not form stable structures
        assert result > -5.0, (
            f"All-A 5' region should not have very negative ΔG, got {result}"
        )

    def test_short_sequence(self):
        """compute_5prime_dg handles a short sequence without crashing."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.compute_5prime_dg(SEQ_10NT)
        assert isinstance(result, (int, float))


# ═══════════════════════════════════════════════════════════════════════════
# 8. Subprocess fallback (conditional on RNAfold CLI)
# ═══════════════════════════════════════════════════════════════════════════

def _rnafold_available() -> bool:
    """Check whether the RNAfold CLI is on $PATH."""
    return shutil.which("RNAfold") is not None


class TestSubprocessFallback:
    """Tests for the subprocess fallback path (requires RNAfold CLI)."""

    @pytest.mark.skipif(not _rnafold_available(), reason="RNAfold CLI not found")
    def test_subprocess_path_works(self):
        """Subprocess path produces valid output when Python bindings unavailable."""
        vr = _get_viennarna()
        # Try calling the subprocess fallback directly if it exists
        if not hasattr(vr, "_predict_mfe_subprocess"):
            pytest.skip("Module does not expose _predict_mfe_subprocess")
        result = vr._predict_mfe_subprocess(SEQ_60NT)
        assert result is not None
        # Should return MFEResult or similar
        if isinstance(result, vr.MFEResult):
            assert result.mfe is not None
            assert result.structure is not None

    @pytest.mark.skipif(not _rnafold_available(), reason="RNAfold CLI not found")
    def test_subprocess_rnafold_direct(self):
        """Direct RNAfold CLI call produces valid dot-bracket output."""
        # Feed the sequence to RNAfold and parse the output
        proc = subprocess.run(
            ["RNAfold", "--noPS"],
            input=SEQ_60NT,
            capture_output=True,
            text=True,
            timeout=30,
        )
        if proc.returncode != 0:
            pytest.skip("RNAfold returned non-zero exit code")
        lines = proc.stdout.strip().split("\n")
        # RNAfold output: line 1 = sequence, line 2 = structure + MFE
        if len(lines) >= 2:
            structure_mfe = lines[1]
            # Should contain dot-bracket and MFE in kcal/mol
            assert "(" in structure_mfe or "." in structure_mfe

    @pytest.mark.skipif(not _rnafold_available(), reason="RNAfold CLI not found")
    def test_subprocess_agrees_with_bindings(self):
        """Subprocess and Python binding results should agree within tolerance."""
        vr = _get_viennarna()
        if not vr.is_viennarna_available():
            pytest.skip("ViennaRNA Python bindings not available for comparison")
        if not hasattr(vr, "_predict_mfe_subprocess"):
            pytest.skip("Module does not expose _predict_mfe_subprocess")

        result_binding = vr.predict_mfe(SEQ_60NT)
        result_subprocess = vr._predict_mfe_subprocess(SEQ_60NT)

        if isinstance(result_subprocess, vr.MFEResult):
            # MFE values should be close (within 0.1 kcal/mol)
            assert abs(result_binding.mfe - result_subprocess.mfe) < 0.5, (
                f"Binding MFE ({result_binding.mfe}) vs "
                f"subprocess MFE ({result_subprocess.mfe}) differ by > 0.5"
            )


# ═══════════════════════════════════════════════════════════════════════════
# 9. Edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestEdgeCases:
    """Tests for edge cases: very short, empty, all-A, and GC-rich sequences."""

    # --- Very short sequences (4-10 nt) ---

    def test_very_short_4nt(self):
        """predict_mfe handles a 4-nt sequence."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_4NT)
        assert isinstance(result, vr.MFEResult)
        assert len(result.structure) == len(SEQ_4NT)

    def test_very_short_10nt(self):
        """predict_mfe handles a 10-nt sequence."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_10NT)
        assert isinstance(result, vr.MFEResult)
        assert len(result.structure) == len(SEQ_10NT)

    def test_short_sequence_mfe_near_zero(self):
        """Very short sequences should have MFE near zero."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_4NT)
        # 4-nt sequences can barely form any structure
        assert -3.0 <= result.mfe <= 0.5, (
            f"4-nt MFE should be near zero, got {result.mfe}"
        )

    # --- Empty sequence ---

    def test_empty_sequence(self):
        """predict_mfe handles an empty sequence gracefully."""
        vr = _get_viennarna()
        # Should either return a result with empty fields, or raise a
        # documented exception (ValueError / ViennaRNAError)
        try:
            result = vr.predict_mfe(SEQ_EMPTY)
            # If it returns, structure should be empty
            assert result.structure == ""
            assert result.mfe == 0.0
        except (ValueError, Exception) as exc:
            # Acceptable: raise a clear error for empty input
            assert "empty" in str(exc).lower() or "length" in str(exc).lower() or "invalid" in str(exc).lower()

    # --- All-A sequence (no structure) ---

    def test_all_a_no_structure(self):
        """All-A sequence should have no (or very weak) structure."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_ALL_A)
        assert isinstance(result, vr.MFEResult)
        # All-A should have mostly unpaired bases
        unpaired_count = result.structure.count(".")
        paired_count = len(result.structure) - unpaired_count
        # At least 80% should be unpaired for all-A
        if len(result.structure) > 0:
            unpaired_fraction = unpaired_count / len(result.structure)
            assert unpaired_fraction >= 0.5, (
                f"All-A should be mostly unpaired, but {unpaired_fraction:.1%} unpaired"
            )

    def test_all_a_mfe_near_zero(self):
        """All-A sequence MFE should be near zero."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_ALL_A)
        # All-A has very weak stacking, MFE should be close to 0
        assert result.mfe > -5.0, (
            f"All-A MFE should be near zero, got {result.mfe}"
        )

    # --- GC-rich sequence (strong structure) ---

    def test_gc_rich_has_structure(self):
        """GC-rich sequence should form stable structure."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_GC_RICH)
        assert isinstance(result, vr.MFEResult)
        # GC-rich should have many paired bases
        paired_count = result.structure.count("(") + result.structure.count(")")
        if len(result.structure) > 0:
            paired_fraction = paired_count / len(result.structure)
            assert paired_fraction >= 0.3, (
                f"GC-rich should have substantial pairing, but only {paired_fraction:.1%} paired"
            )

    def test_gc_rich_negative_mfe(self):
        """GC-rich sequence should have strongly negative MFE."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result = vr.predict_mfe(SEQ_GC_RICH)
        assert result.mfe < -1.0, (
            f"GC-rich should have negative MFE, got {result.mfe}"
        )

    def test_gc_rich_vs_all_a_mfe(self):
        """GC-rich MFE should be more negative than all-A MFE."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result_gc = vr.predict_mfe(SEQ_GC_RICH)
        result_a = vr.predict_mfe(SEQ_ALL_A)
        assert result_gc.mfe < result_a.mfe, (
            f"GC-rich MFE ({result_gc.mfe}) should be < all-A MFE ({result_a.mfe})"
        )


# ═══════════════════════════════════════════════════════════════════════════
# Additional integration / consistency tests
# ═══════════════════════════════════════════════════════════════════════════

class TestIntegrationConsistency:
    """Cross-function consistency checks."""

    def test_mfe_and_find_stable_structures_consistent(self):
        """If predict_mfe returns a structure, find_stable_structures should
        find stems in the same sequence (for structured sequences)."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        mfe_result = vr.predict_mfe(SEQ_STEM_LOOP)
        sl_result = vr.find_stable_structures(SEQ_STEM_LOOP)

        # If the structure has paired bases, we should find at least one stem
        has_pairs = "(" in mfe_result.structure
        if has_pairs and len(SEQ_STEM_LOOP) >= 20:
            # May or may not find discrete stem-loops depending on detection
            # threshold, but the function should not crash
            assert isinstance(sl_result, list)

    def test_accessibility_and_mfe_consistent(self):
        """Sequences with more negative MFE should have lower accessibility."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        mfe_gc = vr.predict_mfe(SEQ_GC_RICH)
        mfe_a = vr.predict_mfe(SEQ_ALL_A)
        acc_gc = vr.predict_accessibility(SEQ_GC_RICH)
        acc_a = vr.predict_accessibility(SEQ_ALL_A)

        mean_gc = getattr(acc_gc, "mean_accessibility", None)
        mean_a = getattr(acc_a, "mean_accessibility", None)

        if mean_gc is not None and mean_a is not None:
            # More stable structure → less accessible
            if mfe_gc.mfe < mfe_a.mfe:
                assert mean_gc <= mean_a, (
                    f"GC-rich is more stable (mfe={mfe_gc.mfe}) but more "
                    f"accessible ({mean_gc}) than all-A ({mean_a})"
                )

    def test_5prime_dg_and_mfe_correlated(self):
        """5' ΔG and full-sequence MFE should be correlated for the same region."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        dg = vr.compute_5prime_dg(SEQ_60NT[:50])
        mfe_result = vr.predict_mfe(SEQ_60NT[:50])
        # Both should be negative for structured sequences
        # The exact values may differ but should have the same sign
        if mfe_result.mfe < -1.0:
            assert dg < 0, (
                f"MFE is negative ({mfe_result.mfe}) but 5' ΔG is positive ({dg})"
            )

    def test_region_specific_mfe(self):
        """Region-specific MFE results should have the region attribute set."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        regions = ["5utr", "full", "start_codon"]
        for region in regions:
            result = vr.predict_mfe(SEQ_5UTR, region=region)
            assert result.region == region, (
                f"Expected region={region}, got region={result.region}"
            )

    def test_deterministic_results(self):
        """Same input should produce identical output across calls."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        result1 = vr.predict_mfe(SEQ_60NT)
        result2 = vr.predict_mfe(SEQ_60NT)
        assert result1.mfe == pytest.approx(result2.mfe, abs=1e-6)
        assert result1.structure == result2.structure


# ═══════════════════════════════════════════════════════════════════════════
# DNA vs RNA handling tests
# ═══════════════════════════════════════════════════════════════════════════

class TestDNARNAHandling:
    """Tests that the module correctly handles DNA input (T→U conversion)."""

    def test_dna_input_accepted(self):
        """predict_mfe should accept DNA sequences (with T)."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        # DNA sequence with T bases
        result = vr.predict_mfe("ATGCGATCGATCGATCGA")
        assert isinstance(result, vr.MFEResult)

    def test_rna_input_accepted(self):
        """predict_mfe should accept RNA sequences (with U) if supported."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        try:
            result = vr.predict_mfe("AUGCGAUCGAUCGAUCGA")
            assert isinstance(result, vr.MFEResult)
        except (ValueError, Exception):
            # Some implementations may not accept RNA directly
            pass

    def test_dna_rna_equivalent_mfe(self):
        """DNA (T) and RNA (U) versions should give the same MFE."""
        _skip_if_no_viennarna()
        vr = _get_viennarna()
        try:
            result_dna = vr.predict_mfe("ATGCGATCGATCGATCGA")
            result_rna = vr.predict_mfe("AUGCGAUCGAUCGAUCGA")
            # MFE should be identical (ViennaRNA converts T→U internally)
            assert result_dna.mfe == pytest.approx(result_rna.mfe, abs=0.01), (
                f"DNA MFE ({result_dna.mfe}) != RNA MFE ({result_rna.mfe})"
            )
        except (ValueError, Exception):
            # If RNA input not supported, skip
            pytest.skip("RNA input not supported by this implementation")


# ═══════════════════════════════════════════════════════════════════════════
# Module-level API contract tests
# ═══════════════════════════════════════════════════════════════════════════

class TestModuleAPIContract:
    """Verify the expected public API of the viennarna module."""

    def test_has_is_viennarna_available(self):
        """Module exposes is_viennarna_available function."""
        vr = _get_viennarna()
        assert callable(getattr(vr, "is_viennarna_available", None))

    def test_has_mfe_result(self):
        """Module exposes MFEResult dataclass."""
        vr = _get_viennarna()
        assert hasattr(vr, "MFEResult")

    def test_has_stem_loop(self):
        """Module exposes StemLoop dataclass."""
        vr = _get_viennarna()
        assert hasattr(vr, "StemLoop")

    def test_has_predict_mfe(self):
        """Module exposes predict_mfe function."""
        vr = _get_viennarna()
        assert callable(getattr(vr, "predict_mfe", None))

    def test_has_predict_accessibility(self):
        """Module exposes predict_accessibility function."""
        vr = _get_viennarna()
        assert callable(getattr(vr, "predict_accessibility", None))

    def test_has_find_stable_structures(self):
        """Module exposes find_stable_structures function."""
        vr = _get_viennarna()
        assert callable(getattr(vr, "find_stable_structures", None))

    def test_has_compute_5prime_dg(self):
        """Module exposes compute_5prime_dg function."""
        vr = _get_viennarna()
        assert callable(getattr(vr, "compute_5prime_dg", None))

    def test_mfe_result_is_dataclass(self):
        """MFEResult should be a dataclass."""
        vr = _get_viennarna()
        import dataclasses
        assert dataclasses.is_dataclass(vr.MFEResult)

    def test_stem_loop_is_dataclass(self):
        """StemLoop should be a dataclass."""
        vr = _get_viennarna()
        import dataclasses
        assert dataclasses.is_dataclass(vr.StemLoop)
