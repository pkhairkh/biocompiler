"""Test that immunogenicity predicates are honest about fabricated scores.

TIGHTEN-4
---------
The default MHC binding PSSMs in
:mod:`biocompiler.immunogenicity.core` are hand-crafted approximations
(the code comments literally say
``"guessed/approximate scores, NOT scores derived from experimental
binding data"``).  Until a real binding-data predictor (NetMHCpan or
MHCflurry) is installed, :func:`compute_immunogenicity` MUST:

1. Compute the guessed-PSSM scores for reference, but
2. Return ``verdict=Verdict.UNCERTAIN`` with ``reason="fabricated_scores"``
   and ``data_source="guessed_pssm"``,
3. NEVER return ``Verdict.PASS`` or ``Verdict.FAIL`` on the basis of
   fabricated scores, and
4. Carry a human-readable ``message`` that explains why the verdict is
   UNCERTAIN.

These tests pin that contract.
"""
from __future__ import annotations

import pytest

from biocompiler.immunogenicity.core import (
    ImmunogenicityResult,
    _check_real_binding_data_available,
    compute_immunogenicity,
)
from biocompiler.shared.types import Verdict

# Human hemoglobin subunit beta (HBB) — a real protein sequence
# used elsewhere in the test-suite.  30 aa N-terminus is enough to
# trigger the full T-cell / B-cell epitope pipeline.
HBB_FRAGMENT = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"

# A short peptide that still validates as a protein (>= 1 aa).
SHORT_PEPTIDE = "MVHLTPEEK"


class TestImmunogenicityHonesty:
    """Honesty-mode contract for :func:`compute_immunogenicity`."""

    def test_default_returns_uncertain(self):
        """Without real binding data, must return UNCERTAIN."""
        result = compute_immunogenicity(HBB_FRAGMENT, "human")
        assert result.verdict == Verdict.UNCERTAIN
        assert result.reason == "fabricated_scores"

    def test_never_pass_on_fabricated_data(self):
        """Must never return PASS on fabricated scores."""
        result = compute_immunogenicity(HBB_FRAGMENT, "human")
        assert result.verdict != Verdict.PASS

    def test_never_fail_on_fabricated_data(self):
        """Must never return FAIL on fabricated scores either.

        FAIL would be just as misleading as PASS — both imply the
        fabricated scores are trustworthy enough to issue a definite
        verdict.  UNCERTAIN is the only honest answer.
        """
        result = compute_immunogenicity(HBB_FRAGMENT, "human")
        assert result.verdict != Verdict.FAIL

    def test_never_likely_pass_or_likely_fail_on_fabricated_data(self):
        """Definite-ish verdicts (LIKELY_PASS / LIKELY_FAIL) are also
        forbidden on fabricated scores — only UNCERTAIN is allowed.
        """
        result = compute_immunogenicity(HBB_FRAGMENT, "human")
        assert result.verdict not in (Verdict.LIKELY_PASS, Verdict.LIKELY_FAIL)

    def test_scores_still_provided(self):
        """Guessed scores still in result for reference."""
        result = compute_immunogenicity(HBB_FRAGMENT, "human")
        # The result object must expose a `scores` attribute, and it
        # must be a non-empty dict so downstream tooling can still
        # inspect the (clearly-labelled) approximate values.
        assert hasattr(result, "scores")
        assert result.scores is not None
        assert isinstance(result.scores, dict)
        assert len(result.scores) > 0
        # The three canonical sub-scores must be present.
        for key in ("overall", "t_cell", "b_cell"):
            assert key in result.scores, f"missing score key: {key}"
            assert isinstance(result.scores[key], (int, float))

    def test_data_source_labelled_guessed_pssm(self):
        """``data_source`` must clearly say these are guessed PSSMs."""
        result = compute_immunogenicity(HBB_FRAGMENT, "human")
        assert result.data_source == "guessed_pssm"

    def test_clear_warning(self):
        """UNCERTAIN result must explain why."""
        result = compute_immunogenicity(SHORT_PEPTIDE, "human")
        msg = getattr(result, "message", "") or str(result)
        assert (
            "NOT real" in msg
            or "approximate" in msg.lower()
            or "UNCERTAIN" in msg
        ), f"message does not warn about fabricated scores: {msg!r}"

    def test_message_mentions_install_hint(self):
        """The UNCERTAIN message should point the user at a real
        predictor (NetMHCpan / MHCflurry) so they know how to upgrade.
        """
        result = compute_immunogenicity(HBB_FRAGMENT, "human")
        msg = getattr(result, "message", "") or ""
        msg_lower = msg.lower()
        assert "netmhcpan" in msg_lower or "mhcflurry" in msg_lower, (
            f"message should mention NetMHCpan or MHCflurry: {msg!r}"
        )

    def test_use_real_data_without_predictor_returns_uncertain(self):
        """``use_real_data=True`` but no real predictor installed -> UNCERTAIN.

        This is the same honesty contract, just signalled via a
        different ``reason`` ("no_real_predictor") so the caller can
        distinguish "I didn't ask for real data" from "I asked but
        none is available".
        """
        # If a real predictor happens to be installed in the test
        # environment, this test is a no-op (we can't fabricate the
        # absence of NetMHCpan / MHCflurry from here).  We still
        # assert the honesty contract holds either way.
        result = compute_immunogenicity(
            HBB_FRAGMENT, "human", use_real_data=True
        )
        if _check_real_binding_data_available():
            # Real predictor is installed — verdict may be PASS/FAIL/UNCERTAIN.
            # We only assert the result is well-formed.
            assert result.verdict in (
                Verdict.PASS,
                Verdict.LIKELY_PASS,
                Verdict.UNCERTAIN,
                Verdict.LIKELY_FAIL,
                Verdict.FAIL,
            )
        else:
            # No real predictor — must be UNCERTAIN, never PASS / FAIL.
            assert result.verdict == Verdict.UNCERTAIN
            assert result.reason == "no_real_predictor"
            assert result.verdict != Verdict.PASS
            assert result.verdict != Verdict.FAIL

    def test_result_is_immunogenicity_result_instance(self):
        """Sanity-check the return type so attribute access is safe."""
        result = compute_immunogenicity(HBB_FRAGMENT, "human")
        assert isinstance(result, ImmunogenicityResult)

    def test_verdict_field_present_on_result(self):
        """The honesty fields must be exposed on the result object."""
        result = compute_immunogenicity(HBB_FRAGMENT, "human")
        for attr in ("verdict", "reason", "message", "data_source", "scores"):
            assert hasattr(result, attr), f"result missing honesty field: {attr}"

    def test_default_verdict_constant_is_uncertain(self):
        """The ImmunogenicityResult dataclass default verdict must be
        UNCERTAIN (so any code path that forgets to set it explicitly
        still defaults to the honest answer).
        """
        # Build a bare result with only the required positional fields.
        # The honesty fields (verdict, reason, message, data_source,
        # scores) all have defaults — those defaults must encode the
        # honest "fabricated_scores" state, not PASS.
        bare = ImmunogenicityResult(
            sequence="",
            primary_score=0.0,
            classification="low",
        )
        assert bare.verdict == Verdict.UNCERTAIN
        assert bare.reason == "fabricated_scores"
        assert bare.data_source == "guessed_pssm"
        assert "NOT real" in bare.message or "approximate" in bare.message.lower()

    def test_self_protein_short_circuit_still_honest(self):
        """Self-protein auto-PASS is allowed because it does NOT rely
        on the fabricated PSSM scores — it's a biology-based shortcut
        (host organism tolerates its own proteins).  The result must
        clearly say ``data_source="self_protein"`` so the verdict is
        attributable to the self-protein rule, not to the PSSM scores.
        """
        result = compute_immunogenicity(
            HBB_FRAGMENT, "human", is_self_protein=True
        )
        assert result.verdict == Verdict.PASS
        assert result.reason == "self_protein_no_assessment"
        assert result.data_source == "self_protein"
        # Make sure the message does NOT claim the PSSMs are real.
        msg_lower = (result.message or "").lower()
        assert "guessed" not in msg_lower or "self" in msg_lower


if __name__ == "__main__":
    # Allow running this test file directly: `python -m pytest` is
    # preferred, but `python tests/test_immunogenicity_honest.py`
    # also works for quick smoke-tests.
    pytest.main([__file__, "-v"])
