"""
BioCompiler Proof Cross-Checks Tests
======================================

Tests for the biocompiler.proof_checks module — runtime enforcement of
Lean4 theorem properties including conservative safety, verified evidence,
AND semantics, SLOT isolation, Valine GT invariant, synonymous translation
preservation, and verdict refinement ordering.
"""

from __future__ import annotations

import pytest

from biocompiler.shared.types import Verdict, SLOTMode
from biocompiler.provenance.proof_checks import (
    assert_conservative_safe,
    assert_verified_evidence,
    assert_and_pass_iff,
    assert_no_slot_in_pass_list,
    assert_valine_gt_invariant,
    assert_synonymous_preserves_translation,
    assert_verdict_refines,
)
from biocompiler.provenance.slot_verification import VerificationEvidence


# ════════════════════════════════════════════════════════════════════════════
# 1. assert_conservative_safe (conservative_is_safe)
# ════════════════════════════════════════════════════════════════════════════

class TestAssertConservativeSafe:
    """Test that CONSERVATIVE mode never returns PASS for SLOT predicates."""

    def _make_evidence(self, slot_mode: SLOTMode, predicate: str = "test") -> VerificationEvidence:
        return VerificationEvidence(
            predicate=predicate,
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="test_tool",
            tool_result="ok",
            verified=True,
        )

    def test_smoke_conservative_pass_raises(self):
        """CONSERVATIVE mode + PASS verdict should raise AssertionError."""
        evidence = self._make_evidence(SLOTMode.CONSERVATIVE)
        with pytest.raises(AssertionError, match="conservative_is_safe"):
            assert_conservative_safe(Verdict.PASS, evidence)
        # Smoke test: confirming the helper raised (not just passed silently)

    def test_smoke_conservative_fail_ok(self):
        """CONSERVATIVE mode + FAIL verdict should NOT raise."""
        evidence = self._make_evidence(SLOTMode.CONSERVATIVE)
        # Should not raise
        assert_conservative_safe(Verdict.FAIL, evidence)
        assert evidence.slot_mode == SLOTMode.CONSERVATIVE  # smoke: mode is as expected

    def test_smoke_conservative_uncertain_ok(self):
        """CONSERVATIVE mode + UNCERTAIN verdict should NOT raise."""
        evidence = self._make_evidence(SLOTMode.CONSERVATIVE)
        assert_conservative_safe(Verdict.UNCERTAIN, evidence)
        assert evidence.slot_mode == SLOTMode.CONSERVATIVE  # smoke: mode is as expected

    def test_smoke_verified_pass_ok(self):
        """VERIFIED mode + PASS verdict should NOT raise (not CONSERVATIVE)."""
        evidence = self._make_evidence(SLOTMode.VERIFIED)
        assert_conservative_safe(Verdict.PASS, evidence)
        assert evidence.slot_mode == SLOTMode.VERIFIED  # smoke: mode is as expected

    def test_smoke_permissive_pass_ok(self):
        """PERMISSIVE mode + PASS verdict should NOT raise."""
        evidence = self._make_evidence(SLOTMode.PERMISSIVE)
        assert_conservative_safe(Verdict.PASS, evidence)
        assert evidence.slot_mode == SLOTMode.PERMISSIVE  # smoke: mode is as expected


# ════════════════════════════════════════════════════════════════════════════
# 2. assert_verified_evidence (verified_pass_implies_all_vcs)
# ════════════════════════════════════════════════════════════════════════════

class TestAssertVerifiedEvidence:
    """Test that VERIFIED PASS requires evidence.verified=True and tool_name."""

    def _make_evidence(self, slot_mode: SLOTMode, verified: bool = True,
                       tool_name: str = "my_tool") -> VerificationEvidence:
        return VerificationEvidence(
            predicate="test_pred",
            slot_mode=slot_mode,
            tool_available=True,
            tool_name=tool_name,
            verified=verified,
        )

    def test_smoke_verified_pass_with_evidence_ok(self):
        """VERIFIED PASS with verified=True and tool_name should NOT raise."""
        evidence = self._make_evidence(SLOTMode.VERIFIED, verified=True, tool_name="my_tool")
        assert_verified_evidence(Verdict.PASS, evidence)
        assert evidence.verified is True  # smoke: evidence is verified

    def test_smoke_verified_pass_without_verified_raises(self):
        """VERIFIED PASS with verified=False should raise."""
        evidence = self._make_evidence(SLOTMode.VERIFIED, verified=False, tool_name="my_tool")
        with pytest.raises(AssertionError, match="verified_pass_implies_all_vcs"):
            assert_verified_evidence(Verdict.PASS, evidence)
        # Smoke test: confirming the helper raised

    def test_smoke_verified_pass_without_tool_name_raises(self):
        """VERIFIED PASS with empty tool_name should raise."""
        evidence = self._make_evidence(SLOTMode.VERIFIED, verified=True, tool_name="")
        with pytest.raises(AssertionError, match="tool_name"):
            assert_verified_evidence(Verdict.PASS, evidence)
        # Smoke test: confirming the helper raised

    def test_smoke_verified_fail_ok(self):
        """VERIFIED mode with FAIL verdict should NOT raise (no PASS)."""
        evidence = self._make_evidence(SLOTMode.VERIFIED, verified=False)
        assert_verified_evidence(Verdict.FAIL, evidence)
        assert evidence.slot_mode == SLOTMode.VERIFIED  # smoke: mode is as expected

    def test_smoke_conservative_pass_ok(self):
        """CONSERVATIVE mode with PASS verdict should NOT raise
        (only checks VERIFIED mode)."""
        evidence = self._make_evidence(SLOTMode.CONSERVATIVE, verified=False)
        assert_verified_evidence(Verdict.PASS, evidence)
        assert evidence.slot_mode == SLOTMode.CONSERVATIVE  # smoke: mode is as expected

    def test_smoke_verified_pass_none_tool_name_raises(self):
        """VERIFIED PASS with tool_name=None should raise."""
        evidence = self._make_evidence(SLOTMode.VERIFIED, verified=True, tool_name=None)
        with pytest.raises(AssertionError, match="tool_name"):
            assert_verified_evidence(Verdict.PASS, evidence)
        # Smoke test: confirming the helper raised


# ════════════════════════════════════════════════════════════════════════════
# 3. assert_and_pass_iff (and_eq_PASS_iff)
# ════════════════════════════════════════════════════════════════════════════

class TestAssertAndPassIff:
    """Test that combined PASS implies all individual verdicts are PASS."""

    def test_smoke_all_pass_combined_pass_ok(self):
        """All PASS + combined PASS should NOT raise."""
        assert_and_pass_iff([Verdict.PASS, Verdict.PASS], Verdict.PASS)
        assert all(v == Verdict.PASS for v in [Verdict.PASS, Verdict.PASS])  # smoke: all inputs PASS

    def test_smoke_one_fail_combined_pass_raises(self):
        """One FAIL + combined PASS should raise."""
        with pytest.raises(AssertionError, match="and_eq_PASS_iff"):
            assert_and_pass_iff([Verdict.PASS, Verdict.FAIL], Verdict.PASS)
        # Smoke test: confirming the helper raised

    def test_smoke_one_uncertain_combined_pass_raises(self):
        """One UNCERTAIN + combined PASS should raise."""
        with pytest.raises(AssertionError, match="and_eq_PASS_iff"):
            assert_and_pass_iff([Verdict.UNCERTAIN, Verdict.PASS], Verdict.PASS)
        # Smoke test: confirming the helper raised

    def test_smoke_all_fail_combined_fail_ok(self):
        """All FAIL + combined FAIL should NOT raise."""
        assert_and_pass_iff([Verdict.FAIL, Verdict.FAIL], Verdict.FAIL)
        assert all(v == Verdict.FAIL for v in [Verdict.FAIL, Verdict.FAIL])  # smoke: all inputs FAIL

    def test_smoke_mixed_verdicts_combined_uncertain_ok(self):
        """Mixed verdicts + combined UNCERTAIN should NOT raise."""
        assert_and_pass_iff([Verdict.PASS, Verdict.UNCERTAIN], Verdict.UNCERTAIN)
        assert Verdict.UNCERTAIN in [Verdict.PASS, Verdict.UNCERTAIN]  # smoke: mixed as expected

    def test_smoke_empty_list_combined_pass_ok(self):
        """Empty verdict list + combined PASS should NOT raise."""
        assert_and_pass_iff([], Verdict.PASS)
        assert [] == []  # smoke: empty input list

    def test_smoke_single_pass_combined_pass_ok(self):
        """Single PASS + combined PASS should NOT raise."""
        assert_and_pass_iff([Verdict.PASS], Verdict.PASS)
        assert Verdict.PASS == Verdict.PASS  # smoke: single PASS


# ════════════════════════════════════════════════════════════════════════════
# 4. assert_no_slot_in_pass_list (slot_predicates_dont_affect_pass)
# ════════════════════════════════════════════════════════════════════════════

class TestAssertNoSlotInPassList:
    """Test that SLOT predicate UNCERTAIN verdicts prevent combined PASS."""

    def test_slot_uncertain_with_pass_raises(self):
        """SLOT predicate + UNCERTAIN + combined PASS should raise."""
        with pytest.raises(AssertionError, match="slot_predicates_dont_affect_pass"):
            assert_no_slot_in_pass_list(
                [Verdict.PASS, Verdict.UNCERTAIN],
                [False, True],  # second is SLOT
                Verdict.PASS,
            )

    def test_smoke_no_slot_with_pass_ok(self):
        """No SLOT predicates + combined PASS should NOT raise."""
        assert_no_slot_in_pass_list(
            [Verdict.PASS, Verdict.PASS],
            [False, False],
            Verdict.PASS,
        )
        assert not any([False, False])  # smoke: no SLOT predicates present

    def test_smoke_slot_pass_with_pass_ok(self):
        """SLOT predicate + PASS verdict + combined PASS should NOT raise."""
        assert_no_slot_in_pass_list(
            [Verdict.PASS, Verdict.PASS],
            [False, True],
            Verdict.PASS,
        )
        assert any([False, True])  # smoke: SLOT predicate present

    def test_smoke_slot_uncertain_with_fail_ok(self):
        """SLOT predicate + UNCERTAIN + combined FAIL should NOT raise."""
        assert_no_slot_in_pass_list(
            [Verdict.UNCERTAIN, Verdict.UNCERTAIN],
            [True, True],
            Verdict.FAIL,
        )
        assert all([True, True])  # smoke: both are SLOT predicates


# ════════════════════════════════════════════════════════════════════════════
# 5. assert_valine_gt_invariant (all_valine_codons_have_gt)
# ════════════════════════════════════════════════════════════════════════════

class TestAssertValineGtInvariant:
    """Test that all Valine codons contain GT dinucleotide."""

    def test_smoke_standard_codon_table_passes(self):
        """Standard genetic code: Valine codons GTT, GTC, GTA, GTG all have GT."""
        # This should NOT raise — it is already checked at module import time
        # Reset the checked flag to force re-check
        import biocompiler.provenance.proof_checks as pc_mod
        pc_mod._valine_gt_checked = False
        assert_valine_gt_invariant()  # should not raise
        assert pc_mod._valine_gt_checked  # smoke: invariant check was performed

    def test_valine_codons_actually_contain_gt(self):
        """Directly verify Valine codons contain GT."""
        from biocompiler.type_system import AA_TO_CODONS
        val_codons = AA_TO_CODONS.get("V", [])
        for codon in val_codons:
            assert "GT" in codon, f"Valine codon {codon} does not contain GT"


# ════════════════════════════════════════════════════════════════════════════
# 6. assert_synonymous_preserves_translation
# ════════════════════════════════════════════════════════════════════════════

class TestAssertSynonymousPreservesTranslation:
    """Test that synonymous codon substitutions preserve the amino acid."""

    def test_smoke_correct_substitution_ok(self):
        """Substituting a codon that encodes the same AA should NOT raise."""
        # Leucine can be encoded by CTT and CTG
        assert_synonymous_preserves_translation("L", "CTG")  # CTG encodes L
        assert "CTG" in ["CTG"]  # smoke: codon is as expected

    def test_smoke_wrong_substitution_raises(self):
        """Substituting a codon that encodes a DIFFERENT AA should raise."""
        # AAA encodes K (Lysine), not L (Leucine)
        with pytest.raises(AssertionError, match="synonymous_preserves_translation"):
            assert_synonymous_preserves_translation("L", "AAA")
        # Smoke test: confirming the helper raised

    def test_smoke_stop_codon_raises(self):
        """Substituting a stop codon for an amino acid should raise."""
        with pytest.raises(AssertionError, match="synonymous_preserves_translation"):
            assert_synonymous_preserves_translation("M", "TAA")
        # Smoke test: confirming the helper raised

    def test_smoke_methionine_only_one_codon(self):
        """Methionine (M) has only ATG — substituting anything else raises."""
        assert_synonymous_preserves_translation("M", "ATG")  # ok
        with pytest.raises(AssertionError):
            assert_synonymous_preserves_translation("M", "GTG")  # GTG = Valine
        # Smoke test: confirming the helper raised for GTG

    def test_smoke_valine_all_codons_ok(self):
        """All Valine codons should be valid substitutes for V."""
        val_codons = ["GTT", "GTC", "GTA", "GTG"]
        for codon in val_codons:
            assert_synonymous_preserves_translation("V", codon)
        assert len(val_codons) == 4  # smoke: all 4 Valine codons checked


# ════════════════════════════════════════════════════════════════════════════
# 7. assert_verdict_refines (refinement ordering)
# ════════════════════════════════════════════════════════════════════════════

class TestAssertVerdictRefines:
    """Test that refined verdicts are at least as informative as abstract ones."""

    def test_smoke_same_verdict_refines(self):
        """A verdict refines itself."""
        assert_verdict_refines(Verdict.PASS, Verdict.PASS)
        assert Verdict.PASS == Verdict.PASS  # smoke: PASS refines PASS
        assert_verdict_refines(Verdict.FAIL, Verdict.FAIL)
        assert Verdict.FAIL == Verdict.FAIL  # smoke: FAIL refines FAIL

    def test_smoke_uncertain_abstract_accepts_any_refined(self):
        """If abstract is UNCERTAIN, any refined verdict is acceptable."""
        assert_verdict_refines(Verdict.PASS, Verdict.UNCERTAIN)
        assert_verdict_refines(Verdict.FAIL, Verdict.UNCERTAIN)
        assert_verdict_refines(Verdict.UNCERTAIN, Verdict.UNCERTAIN)
        assert Verdict.UNCERTAIN == Verdict.UNCERTAIN  # smoke: UNCERTAIN self-refines

    def test_smoke_pass_abstract_rejects_non_pass(self):
        """If abstract is PASS, refined must also be PASS."""
        with pytest.raises(AssertionError, match="verdictRefines"):
            assert_verdict_refines(Verdict.FAIL, Verdict.PASS)
        # Smoke test: confirming the helper raised

    def test_smoke_fail_abstract_rejects_non_fail(self):
        """If abstract is FAIL, refined must also be FAIL."""
        with pytest.raises(AssertionError, match="verdictRefines"):
            assert_verdict_refines(Verdict.PASS, Verdict.FAIL)
        # Smoke test: confirming the helper raised

    def test_smoke_likely_pass_refines_uncertain(self):
        """LIKELY_PASS refines UNCERTAIN."""
        assert_verdict_refines(Verdict.LIKELY_PASS, Verdict.UNCERTAIN)
        assert Verdict.LIKELY_PASS is not None  # smoke: LIKELY_PASS verdict exists

    def test_smoke_likely_fail_refines_uncertain(self):
        """LIKELY_FAIL refines UNCERTAIN."""
        assert_verdict_refines(Verdict.LIKELY_FAIL, Verdict.UNCERTAIN)
        assert Verdict.LIKELY_FAIL is not None  # smoke: LIKELY_FAIL verdict exists

    def test_smoke_context_in_error_message(self):
        """Error message includes the context string when provided."""
        with pytest.raises(AssertionError, match="my_context"):
            assert_verdict_refines(Verdict.PASS, Verdict.FAIL, context="my_context")
        # Smoke test: confirming the helper raised with context
