"""
BioCompiler Runtime Proof Cross-Checks
=======================================

Defensive assertions that verify properties proven in the Lean4 formal model
hold at runtime in the Python implementation. These are the runtime analogues
of Lean4 theorems — if any assertion fails, it indicates a proof-implementation
gap that needs investigation.

Lean4 theorems enforced at runtime:
  1. conservative_is_safe: CONSERVATIVE mode never returns PASS for SLOT predicates
  2. verified_pass_implies_all_vcs: VERIFIED PASS implies evidence.verified is True
  3. and_eq_PASS_iff: AND(a,b)=PASS iff a=PASS and b=PASS
  4. compositional_soundness: combined PASS implies all individual PASS
  5. slot_predicates_dont_affect_pass: UNCERTAIN in list prevents PASS
  6. all_valine_codons_have_gt: all Valine codons contain GT
  7. synonymous_preserves_translation: synonymous codons encode same AA

Usage:
    from biocompiler.proof_checks import assert_conservative_safe, assert_verified_evidence

Reference:
    proof/BioCompiler/*.lean
    docs/14-SLOT-Proof-Implementation-Gap.md
"""

from __future__ import annotations

import logging
from typing import Optional, TYPE_CHECKING

from .types import Verdict, SLOTMode, five_valued_and, combined_verdict, _VERDICT_ORDER
from .type_system import CODON_TABLE, AA_TO_CODONS

if TYPE_CHECKING:
    from .slot_verification import VerificationEvidence

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════
# 1. conservative_is_safe (SLOTVerification.lean)
# ═══════════════════════════════════════════════════════════════════════════

def assert_conservative_safe(verdict: Verdict, evidence: VerificationEvidence) -> None:
    """Runtime enforcement of Lean4 theorem: conservative_is_safe.

    In CONSERVATIVE mode, SLOT predicates must never return PASS.
    This is the foundational safety property from SLOTVerification.lean:
      "Conservative mode never returns PASS for SLOT predicates."

    Args:
        verdict: The verdict returned by a SLOT predicate evaluation.
        evidence: The accompanying evidence object.

    Raises:
        AssertionError: If CONSERVATIVE mode returns PASS for a SLOT predicate.
    """
    if evidence.slot_mode == SLOTMode.CONSERVATIVE:
        if verdict == Verdict.PASS:
            raise AssertionError(
                f"PROOF VIOLATION: conservative_is_safe — CONSERVATIVE mode "
                f"returned PASS for SLOT predicate '{evidence.predicate}'. "
                f"This violates the Lean4 theorem conservative_is_safe "
                f"(SLOTVerification.lean)."
            )


# ═══════════════════════════════════════════════════════════════════════════
# 2. verified_pass_implies_all_vcs (SLOTVerification.lean)
# ═══════════════════════════════════════════════════════════════════════════

def assert_verified_evidence(verdict: Verdict, evidence: VerificationEvidence) -> None:
    """Runtime enforcement of Lean4 theorem: verified_pass_implies_all_vcs.

    If VERIFIED mode returns PASS, then all verification conditions must hold.
    In Python, this means evidence.verified must be True.

    From SLOTVerification.lean:
      "If evaluateSLOT SLOTMode.verified vctx P seq ctx = PASS,
       then allVCHold vctx (slotVCs P) = true"

    Args:
        verdict: The verdict returned by a SLOT predicate evaluation.
        evidence: The accompanying evidence object.

    Raises:
        AssertionError: If VERIFIED PASS lacks proper evidence.
    """
    if evidence.slot_mode == SLOTMode.VERIFIED and verdict == Verdict.PASS:
        if not evidence.verified:
            raise AssertionError(
                f"PROOF VIOLATION: verified_pass_implies_all_vcs — VERIFIED mode "
                f"returned PASS for '{evidence.predicate}' but evidence.verified=False. "
                f"This violates the Lean4 theorem verified_pass_implies_all_vcs."
            )
        if evidence.tool_name is None or len(evidence.tool_name) == 0:
            raise AssertionError(
                f"PROOF VIOLATION: verified_pass_implies_all_vcs — VERIFIED mode "
                f"returned PASS for '{evidence.predicate}' but no tool_name documented. "
                f"Verification evidence must name the tool that produced the result."
            )


# ═══════════════════════════════════════════════════════════════════════════
# 3. and_eq_PASS_iff (ThreeValued.lean)
# ═══════════════════════════════════════════════════════════════════════════

def assert_and_pass_iff(verdicts: list[Verdict], combined: Verdict) -> None:
    """Runtime enforcement of Lean4 theorem: and_eq_PASS_iff.

    AND(a, b) = PASS iff a = PASS and b = PASS.
    Combined verdict PASS implies all individual verdicts are PASS.

    From ThreeValued.lean:
      "Verdict.and init hd = PASS ↔ init = PASS ∧ hd = PASS"

    Args:
        verdicts: List of individual verdicts.
        combined: The combined verdict computed from the list.

    Raises:
        AssertionError: If combined PASS doesn't imply all individual PASS.
    """
    if combined == Verdict.PASS:
        for i, v in enumerate(verdicts):
            if v != Verdict.PASS:
                raise AssertionError(
                    f"PROOF VIOLATION: and_eq_PASS_iff — combined verdict is PASS "
                    f"but element {i} is {v}. This violates the Lean4 theorem "
                    f"foldl_and_pass_implies_all_pass (ThreeValued.lean)."
                )


# ═══════════════════════════════════════════════════════════════════════════
# 4. slot_predicates_dont_affect_pass (Compositional.lean)
# ═══════════════════════════════════════════════════════════════════════════

def assert_no_slot_in_pass_list(
    verdicts: list[Verdict],
    is_slot: list[bool],
    combined: Verdict,
) -> None:
    """Runtime enforcement of Lean4 theorem: slot_predicates_dont_affect_pass.

    If any SLOT predicate is in the list and returns UNCERTAIN,
    the combined verdict cannot be PASS.

    From Compositional.lean:
      "If any SLOT-dependent predicate is in the list, evaluateAll
       cannot return PASS."

    Args:
        verdicts: List of individual verdicts.
        is_slot: Parallel list indicating which predicates are SLOT-dependent.
        combined: The combined verdict.

    Raises:
        AssertionError: If SLOT predicate UNCERTAIN verdicts don't prevent PASS.
    """
    if combined == Verdict.PASS:
        for i, (v, slot) in enumerate(zip(verdicts, is_slot)):
            if slot and v == Verdict.UNCERTAIN:
                raise AssertionError(
                    f"PROOF VIOLATION: slot_predicates_dont_affect_pass — "
                    f"combined verdict is PASS but element {i} is a SLOT predicate "
                    f"with verdict UNCERTAIN. This violates the Lean4 theorem "
                    f"slot_predicates_dont_affect_pass (Compositional.lean)."
                )


# ═══════════════════════════════════════════════════════════════════════════
# 5. all_valine_codons_have_gt (Mutagenesis.lean)
# ═══════════════════════════════════════════════════════════════════════════

_valine_gt_checked: bool = False


def assert_valine_gt_invariant() -> None:
    """Runtime enforcement of Lean4 theorem: all_valine_codons_have_gt.

    All four Valine codons must contain the GT dinucleotide.
    This is checked once at module import time.

    From Mutagenesis.lean:
      "Every Valine codon contains GT dinucleotide."

    Raises:
        AssertionError: If any Valine codon lacks GT.
    """
    global _valine_gt_checked
    if _valine_gt_checked:
        return
    _valine_gt_checked = True

    val_codons = AA_TO_CODONS.get("V", [])
    for codon in val_codons:
        if "GT" not in codon:
            raise AssertionError(
                f"PROOF VIOLATION: all_valine_codons_have_gt — "
                f"Valine codon '{codon}' does not contain GT. "
                f"This violates the Lean4 theorem all_valine_codons_have_gt "
                f"(Mutagenesis.lean)."
            )


# ═══════════════════════════════════════════════════════════════════════════
# 6. synonymous_preserves_translation (Mutagenesis.lean)
# ═══════════════════════════════════════════════════════════════════════════

def assert_synonymous_preserves_translation(original_aa: str, new_codon: str) -> None:
    """Runtime enforcement of Lean4 theorem: synonymous_preserves_translation.

    A synonymous codon substitution must preserve the amino acid.

    From Mutagenesis.lean:
      "codonToAA originalCodon = codonToAA newCodon"

    Args:
        original_aa: The amino acid encoded by the original codon.
        new_codon: The replacement codon.

    Raises:
        AssertionError: If the new codon doesn't encode the same amino acid.
    """
    new_aa = CODON_TABLE.get(new_codon)
    if new_aa != original_aa:
        raise AssertionError(
            f"PROOF VIOLATION: synonymous_preserves_translation — "
            f"substituting codon '{new_codon}' (encodes '{new_aa}') "
            f"for amino acid '{original_aa}' changes the protein. "
            f"This violates the Lean4 theorem synonymous_preserves_translation "
            f"(Mutagenesis.lean)."
        )


# ═══════════════════════════════════════════════════════════════════════════
# 7. Refinement ordering (Refinement.lean)
# ═══════════════════════════════════════════════════════════════════════════

def assert_verdict_refines(
    v_refined: Verdict,
    v_abstract: Verdict,
    context: str = "",
) -> None:
    """Runtime enforcement of Lean4 refinement ordering.

    v_refined must be at least as informative as v_abstract.
    From Refinement.lean:
      "verdictRefines v_refined v_abstract := v_abstract = UNCERTAIN ∨ v_refined = v_abstract"

    This is used to verify that VERIFIED/PERMISSIVE verdicts refine CONSERVATIVE verdicts.

    Args:
        v_refined: The more informative verdict (e.g., from VERIFIED mode).
        v_abstract: The less informative verdict (e.g., from CONSERVATIVE mode).
        context: Description of the check for error messages.

    Raises:
        AssertionError: If the refinement ordering is violated.
    """
    if v_abstract != Verdict.UNCERTAIN and v_refined != v_abstract:
        raise AssertionError(
            f"PROOF VIOLATION: verdictRefines — {context} "
            f"refined verdict {v_refined} does not refine abstract verdict "
            f"{v_abstract}. This violates the Lean4 theorem "
            f"verified_refines_conservative (Refinement.lean)."
        )


# ═══════════════════════════════════════════════════════════════════════════
# Initialization: Run invariant checks at import time
# ═══════════════════════════════════════════════════════════════════════════

# Check Valine GT invariant at import time
# This ensures the codon table hasn't been corrupted
try:
    assert_valine_gt_invariant()
    logger.debug("Lean4 proof cross-check: all_valine_codons_have_gt ✓")
except AssertionError:
    # Re-raise in development, log in production
    logger.error("Lean4 proof cross-check FAILED: all_valine_codons_have_gt")
    raise
