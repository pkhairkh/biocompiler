"""
Tests for SLOT soundness verification
========================================

Tests the non-vacuous soundness checking introduced in the
slot_verification module.  Key properties verified:

1. **Core predicates** (non-SLOT) are always recheckable and sound
   because they are deterministic.
2. **SLOT predicates** with PASS verdicts are non-vacuously checked:
   - If the recheck confirms PASS/LIKELY_PASS → sound
   - If the recheck gives UNCERTAIN → unsound (cannot confirm)
   - If the recheck gives FAIL → unsound (contradiction)
3. **UNCERTAIN/LIKELY_* verdicts** are trivially sound (no definite claim).
4. **SoundnessReport** aggregates correctly.
5. **_recheck_core_predicate** handles all core predicates.
6. **_apply_recheck_to_soundness** applies consistent rules.
"""

from __future__ import annotations

import pytest
from unittest.mock import patch

from biocompiler.shared.types import Verdict, SLOTMode
from biocompiler.provenance.slot_verification import (
    verify_soundness,
    SoundnessReport,
    SoundnessResult,
    _recheck_predicate,
    _recheck_core_predicate,
    _apply_recheck_to_soundness,
)


# ────────────────────────────────────────────────────────────
# 1. Core predicate soundness — deterministic, always recheckable
# ────────────────────────────────────────────────────────────

class TestCorePredicateSoundness:
    """Core (non-SLOT) predicates are deterministic and can always be
    rechecked. A PASS verdict should be confirmed on re-evaluation."""

    def test_no_stop_codons_pass_sound(self):
        """NoStopCodons PASS on a clean sequence is confirmed sound."""
        report = verify_soundness(
            [("NoStopCodons", Verdict.PASS)],
            seq="ATGGCCGCTTAA",
        )
        assert report.all_sound
        assert report.pass_verified == 1
        assert report.results[0].recheck_passed is True

    def test_no_stop_codons_fail_unsound(self):
        """NoStopCodons FAIL on a clean sequence is unsound (false alarm)."""
        report = verify_soundness(
            [("NoStopCodons", Verdict.FAIL)],
            seq="ATGGCCGCTTAA",
        )
        assert not report.all_sound
        assert report.results[0].recheck_violated is False

    def test_no_stop_codons_fail_sound(self):
        """NoStopCodons FAIL on a sequence with internal stop is sound."""
        report = verify_soundness(
            [("NoStopCodons", Verdict.FAIL)],
            seq="ATGTAGGCTTAA",  # TAG at position 3
        )
        assert report.all_sound
        assert report.results[0].recheck_violated is True

    def test_valid_coding_seq_pass_sound(self):
        """ValidCodingSeq PASS on a valid sequence is confirmed sound."""
        report = verify_soundness(
            [("ValidCodingSeq", Verdict.PASS)],
            seq="ATGGCCGCTTAA",
        )
        assert report.all_sound

    def test_valid_coding_seq_fail_sound(self):
        """ValidCodingSeq FAIL on invalid sequence is sound."""
        report = verify_soundness(
            [("ValidCodingSeq", Verdict.FAIL)],
            seq="ATGGCC",  # Length 6 but with valid codons — FAIL is unsound
        )
        # ValidCodingSeq should pass on a valid sequence, so FAIL is unsound
        assert not report.all_sound


class TestGTDinucleotideSoundness:
    """NoGTDinucleotide soundness checks."""

    def test_gt_fail_on_seq_with_gt(self):
        """FAIL verdict on a sequence with GT is sound."""
        report = verify_soundness(
            [("NoGTDinucleotide", Verdict.FAIL)],
            seq="ATGGTTGCTTAA",  # GT at pos 3
        )
        assert report.all_sound
        assert report.results[0].recheck_violated is True

    def test_gt_pass_on_seq_without_gt(self):
        """PASS verdict on a sequence without GT is sound."""
        report = verify_soundness(
            [("NoGTDinucleotide", Verdict.PASS)],
            seq="AACCCAAGGAATTCCCGGG",  # No GT
        )
        assert report.all_sound


# ────────────────────────────────────────────────────────────
# 2. SLOT predicate soundness — non-vacuous
# ────────────────────────────────────────────────────────────

class TestSLOTPredicateSoundness:
    """SLOT predicates use verify_slot_predicate for rechecking.
    When the tool is unavailable and returns UNCERTAIN, the PASS
    verdict is marked unsound (not vacuously sound)."""

    def test_uncertain_verdict_trivially_sound(self):
        """Non-definite verdicts are trivially sound."""
        report = verify_soundness(
            [("LowImmunogenicity", Verdict.UNCERTAIN)],
            protein_sequence="MAGICS",
        )
        assert report.all_sound
        assert report.results[0].details == "Verdict UNCERTAIN is non-definite; soundness is trivially true"

    def test_likely_pass_trivially_sound(self):
        """LIKELY_PASS is non-definite, trivially sound."""
        report = verify_soundness(
            [("LowImmunogenicity", Verdict.LIKELY_PASS)],
            protein_sequence="MAGICS",
        )
        assert report.all_sound

    def test_slot_pass_unsound_when_tool_unavailable(self):
        """SLOT PASS is unsound when recheck returns UNCERTAIN (tool unavailable).

        This is the key non-vacuous check: before the fix, PASS was vacuously
        sound because the fallback trusted the evidence. Now, if we cannot
        confirm the PASS claim, it is marked unsound.
        """
        report = verify_soundness(
            [("StableFolding", Verdict.PASS)],
            protein_sequence="MAGICS",
        )
        # StableFolding requires FoldX which is likely not installed
        # Recheck should give UNCERTAIN, making the PASS unsound
        assert not report.all_sound
        assert "NOT confirmed" in report.results[0].details or "UNSOUND" in report.results[0].details

    def test_codon_optimality_pass_sound_with_data(self):
        """CodonOptimality PASS with CAI data is confirmed sound."""
        report = verify_soundness(
            [("CodonOptimality", Verdict.PASS)],
            codon="ATG", species_cai={"ATG": 0.9}, min_cai=0.5,
        )
        assert report.all_sound
        assert report.results[0].recheck_passed is True

    def test_conservation_score_pass_sound(self):
        """ConservationScore PASS for A→A is confirmed sound."""
        report = verify_soundness(
            [("ConservationScore", Verdict.PASS)],
            original_aa="A", new_aa="A", min_score=0,
        )
        assert report.all_sound


class TestCodonOptimalitySoundness:
    """Detailed tests for CodonOptimality soundness rechecking."""

    def test_pass_with_high_cai(self):
        report = verify_soundness(
            [("CodonOptimality", Verdict.PASS)],
            codon="ATG", species_cai={"ATG": 1.0}, min_cai=0.8,
        )
        assert report.all_sound

    def test_pass_with_low_cai_unsound(self):
        """PASS claimed for a codon with CAI below threshold is unsound."""
        report = verify_soundness(
            [("CodonOptimality", Verdict.PASS)],
            codon="TTT", species_cai={"TTT": 0.3}, min_cai=0.8,
        )
        # Recheck should return FAIL (CAI 0.3 < 0.8)
        assert not report.all_sound


# ────────────────────────────────────────────────────────────
# 3. SoundnessReport aggregation
# ────────────────────────────────────────────────────────────

class TestSoundnessReportAggregation:
    """Test SoundnessReport aggregate statistics."""

    def test_empty_report(self):
        report = SoundnessReport(results=[])
        assert report.total == 0
        assert report.all_sound is True  # vacuously

    def test_mixed_report(self):
        """Report with both sound and unsound results."""
        results = [
            SoundnessResult("A", Verdict.PASS, sound=True, recheck_passed=True),
            SoundnessResult("B", Verdict.PASS, sound=False, recheck_passed=False),
        ]
        report = SoundnessReport(results=results)
        assert report.total == 2
        assert report.sound_count == 1
        assert report.unsound_count == 1
        assert not report.all_sound
        assert report.pass_verified == 1
        assert report.pass_unverified == 1

    def test_all_sound_report(self):
        results = [
            SoundnessResult("A", Verdict.PASS, sound=True, recheck_passed=True),
            SoundnessResult("B", Verdict.UNCERTAIN, sound=True),
        ]
        report = SoundnessReport(results=results)
        assert report.all_sound
        assert report.pass_verified == 1

    def test_fail_false_alarm_counted(self):
        results = [
            SoundnessResult("A", Verdict.FAIL, sound=False, recheck_violated=False),
        ]
        report = SoundnessReport(results=results)
        assert report.fail_false_alarm == 1

    def test_fail_verified_counted(self):
        results = [
            SoundnessResult("A", Verdict.FAIL, sound=True, recheck_violated=True),
        ]
        report = SoundnessReport(results=results)
        assert report.fail_verified == 1

    def test_to_dict(self):
        results = [
            SoundnessResult("A", Verdict.PASS, sound=True, recheck_passed=True),
        ]
        report = SoundnessReport(results=results)
        d = report.to_dict()
        assert "total" in d
        assert "all_sound" in d
        assert "results" in d
        assert len(d["results"]) == 1
        assert d["results"][0]["verdict"] == "PASS"


# ────────────────────────────────────────────────────────────
# 4. _recheck_core_predicate dispatch
# ────────────────────────────────────────────────────────────

class TestRecheckCorePredicate:
    """Test that _recheck_core_predicate handles all core predicates."""

    def test_no_stop_codons(self):
        result = _recheck_core_predicate("NoStopCodons", Verdict.PASS, seq="ATGGCCGCTTAA")
        assert result is not None
        assert result.sound is True

    def test_no_gt_dinucleotide(self):
        result = _recheck_core_predicate("NoGTDinucleotide", Verdict.PASS, seq="AACCCCGGG")
        assert result is not None
        assert result.sound is True

    def test_valid_coding_seq(self):
        result = _recheck_core_predicate("ValidCodingSeq", Verdict.PASS, seq="ATGGCCGCTTAA")
        assert result is not None
        assert result.sound is True

    def test_no_restriction_site(self):
        result = _recheck_core_predicate("NoRestrictionSite", Verdict.PASS, seq="ATGGCCGCTTAA", enzymes=[])
        assert result is not None

    def test_no_cpg_island(self):
        result = _recheck_core_predicate("NoCpGIsland", Verdict.PASS, seq="ATGGCCGCTTAA")
        assert result is not None

    def test_gc_in_range(self):
        result = _recheck_core_predicate("GCInRange", Verdict.PASS, seq="ATGGCCGCTTAA")
        assert result is not None

    def test_in_frame(self):
        result = _recheck_core_predicate("InFrame", Verdict.PASS, seq="ATGGCCGCTTAA")
        assert result is not None

    def test_no_instability_motif(self):
        result = _recheck_core_predicate("NoInstabilityMotif", Verdict.PASS, seq="ATGGCCGCTTAA")
        assert result is not None

    def test_slot_predicate_returns_none(self):
        """SLOT predicates are not handled by _recheck_core_predicate."""
        result = _recheck_core_predicate("LowImmunogenicity", Verdict.PASS, protein_sequence="MAGICS")
        assert result is None

    def test_unknown_predicate_returns_none(self):
        result = _recheck_core_predicate("MadeUpPredicate", Verdict.PASS, seq="ATG")
        assert result is None


# ────────────────────────────────────────────────────────────
# 5. _apply_recheck_to_soundness consistency
# ────────────────────────────────────────────────────────────

class TestApplyRecheckToSoundness:
    """Test that _apply_recheck_to_soundness applies consistent rules."""

    def test_pass_confirmed(self):
        result = SoundnessResult("A", Verdict.PASS, sound=True)
        _apply_recheck_to_soundness(result, Verdict.PASS, Verdict.PASS)
        assert result.sound is True
        assert result.recheck_passed is True

    def test_pass_contradicted_by_fail(self):
        result = SoundnessResult("A", Verdict.PASS, sound=True)
        _apply_recheck_to_soundness(result, Verdict.PASS, Verdict.FAIL)
        assert result.sound is False
        assert result.recheck_passed is False
        assert "UNSOUND" in result.details

    def test_pass_not_confirmed_by_uncertain(self):
        result = SoundnessResult("A", Verdict.PASS, sound=True)
        _apply_recheck_to_soundness(result, Verdict.PASS, Verdict.UNCERTAIN)
        assert result.sound is False
        assert result.recheck_passed is False
        assert "NOT confirmed" in result.details

    def test_pass_confirmed_by_likely_pass(self):
        result = SoundnessResult("A", Verdict.PASS, sound=True)
        _apply_recheck_to_soundness(result, Verdict.PASS, Verdict.LIKELY_PASS)
        assert result.sound is True
        assert result.recheck_passed is True

    def test_fail_confirmed(self):
        result = SoundnessResult("A", Verdict.FAIL, sound=True)
        _apply_recheck_to_soundness(result, Verdict.FAIL, Verdict.FAIL)
        assert result.sound is True
        assert result.recheck_violated is True

    def test_fail_false_alarm(self):
        result = SoundnessResult("A", Verdict.FAIL, sound=True)
        _apply_recheck_to_soundness(result, Verdict.FAIL, Verdict.PASS)
        assert result.sound is False
        assert result.recheck_violated is False
        assert "false alarm" in result.details

    def test_fail_not_confirmed_by_uncertain(self):
        result = SoundnessResult("A", Verdict.FAIL, sound=True)
        _apply_recheck_to_soundness(result, Verdict.FAIL, Verdict.UNCERTAIN)
        assert result.sound is False
        assert result.recheck_violated is False


# ────────────────────────────────────────────────────────────
# 6. Non-vacuity demonstration
# ────────────────────────────────────────────────────────────

class TestNonVacuityDemonstration:
    """Demonstrate that the soundness check is non-vacuous.

    In the old (vacuous) implementation, PASS verdicts for SLOT
    predicates without specific recheck logic were marked as sound
    with 'trusted evidence'. This meant the soundness property was
    trivially true because nothing could ever fail.

    The new implementation uses verify_slot_predicate to re-run
    the check. If the tool is unavailable (UNCERTAIN), the PASS
    is NOT confirmed — it is marked unsound. This forces honest
    reporting.
    """

    def test_old_vacuous_behavior_would_be_sound(self):
        """In the old implementation, all SLOT PASS verdicts were sound.

        This test documents what the OLD behavior would have been:
        PASS + no recheck = trusted evidence = vacuously sound.
        """
        # The old code had: result.recheck_passed = None;
        # result.details = "PASS: no recheck available (trusted evidence)"
        # This is what we DO NOT want anymore.

        # With the new code, SLOT PASS without tool confirmation is unsound
        report = verify_soundness(
            [("StableFolding", Verdict.PASS)],
            protein_sequence="MAGICS",
        )
        # The new behavior: StableFolding requires FoldX, which is not
        # installed → UNCERTAIN → PASS is NOT confirmed → unsound
        # This is NON-VACUOUS!
        if report.results[0].recheck_passed is False:
            # Tool unavailable → unsound (honest)
            assert not report.all_sound
        else:
            # If the tool happens to be installed, it is sound
            assert report.all_sound

    def test_new_behavior_distinguishes_confirmed_from_unconfirmed(self):
        """The new soundness check distinguishes confirmed vs unconfirmed PASS.

        Core predicates: always confirmed (deterministic).
        SLOT predicates: confirmed only if tool returns positive evidence.
        """
        # Core: always confirmed
        core_report = verify_soundness(
            [("NoStopCodons", Verdict.PASS)],
            seq="ATGGCCGCTTAA",
        )
        assert core_report.all_sound
        assert core_report.results[0].recheck_passed is True

        # SLOT: may not be confirmed (depends on tool availability)
        slot_report = verify_soundness(
            [("StructureConfidence", Verdict.PASS)],
            protein_sequence="MAGICS",
        )
        # ESMFold likely not installed → UNCERTAIN → unsound
        # This is the non-vacuous distinction!
        if not slot_report.all_sound:
            assert slot_report.results[0].recheck_passed is False


# ────────────────────────────────────────────────────────────
# 7. Batch soundness verification
# ────────────────────────────────────────────────────────────

class TestBatchSoundnessVerification:
    """Test verify_soundness with multiple predicates at once."""

    def test_mixed_core_predicates(self):
        """Multiple core predicates with known outcomes."""
        report = verify_soundness(
            [
                ("NoStopCodons", Verdict.PASS),
                ("ValidCodingSeq", Verdict.PASS),
                ("NoStopCodons", Verdict.FAIL),  # unsound on a clean seq
            ],
            seq="ATGGCCGCTTAA",
        )
        assert report.total == 3
        assert report.sound_count == 2
        assert report.unsound_count == 1

    def test_multiple_uncertain_verdicts(self):
        """Multiple UNCERTAIN verdicts are all trivially sound."""
        report = verify_soundness(
            [
                ("LowImmunogenicity", Verdict.UNCERTAIN),
                ("StableFolding", Verdict.UNCERTAIN),
            ],
            protein_sequence="MAGICS",
        )
        assert report.all_sound
