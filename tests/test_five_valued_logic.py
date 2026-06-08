"""
Comprehensive tests for the five-valued logic module.

Tests cover:
- FiveValuedVerdict enum construction and values
- to_three_valued refinement mapping
- confidence_score for all verdicts
- combine_verdicts lattice combination
- verify_five_valued_soundness for all cases
- five_valued_and / five_valued_or operations
- refinement_is_sound bridge function
- Consistency with existing 3-valued Verdict enum
- Edge cases (empty lists, single verdicts)
"""

import pytest

from biocompiler.five_valued_logic import (
    FiveValuedVerdict,
    to_three_valued,
    confidence_score,
    combine_verdicts,
    verify_five_valued_soundness,
    five_valued_and,
    five_valued_or,
    refinement_is_sound,
)
from biocompiler.types import Verdict


# ═══════════════════════════════════════════════════════════════════════
# FiveValuedVerdict Enum Tests
# ═══════════════════════════════════════════════════════════════════════

class TestFiveValuedVerdictEnum:
    """Test FiveValuedVerdict enum construction and properties."""

    def test_all_five_values_exist(self):
        """All five verdict values are defined."""
        assert FiveValuedVerdict.PASS is not None
        assert FiveValuedVerdict.LIKELY_PASS is not None
        assert FiveValuedVerdict.UNCERTAIN is not None
        assert FiveValuedVerdict.LIKELY_FAIL is not None
        assert FiveValuedVerdict.FAIL is not None

    def test_five_distinct_values(self):
        """All five verdicts are distinct."""
        values = list(FiveValuedVerdict)
        assert len(values) == 5
        assert len(set(values)) == 5

    def test_string_values(self):
        """String values match enum names."""
        assert FiveValuedVerdict.PASS.value == "PASS"
        assert FiveValuedVerdict.LIKELY_PASS.value == "LIKELY_PASS"
        assert FiveValuedVerdict.UNCERTAIN.value == "UNCERTAIN"
        assert FiveValuedVerdict.LIKELY_FAIL.value == "LIKELY_FAIL"
        assert FiveValuedVerdict.FAIL.value == "FAIL"

    def test_is_str_enum(self):
        """FiveValuedVerdict is a str enum."""
        assert isinstance(FiveValuedVerdict.PASS, str)
        assert FiveValuedVerdict.PASS == "PASS"

    def test_iteration_order(self):
        """Iteration yields all five values."""
        values = list(FiveValuedVerdict)
        assert FiveValuedVerdict.PASS in values
        assert FiveValuedVerdict.LIKELY_PASS in values
        assert FiveValuedVerdict.UNCERTAIN in values
        assert FiveValuedVerdict.LIKELY_FAIL in values
        assert FiveValuedVerdict.FAIL in values


# ═══════════════════════════════════════════════════════════════════════
# to_three_valued Refinement Mapping Tests
# ═══════════════════════════════════════════════════════════════════════

class TestToThreeValued:
    """Test the refinement mapping from 5-valued to 3-valued logic."""

    def test_pass_maps_to_pass(self):
        assert to_three_valued(FiveValuedVerdict.PASS) == Verdict.PASS

    def test_likely_pass_maps_to_uncertain(self):
        assert to_three_valued(FiveValuedVerdict.LIKELY_PASS) == Verdict.UNCERTAIN

    def test_uncertain_maps_to_uncertain(self):
        assert to_three_valued(FiveValuedVerdict.UNCERTAIN) == Verdict.UNCERTAIN

    def test_likely_fail_maps_to_fail(self):
        assert to_three_valued(FiveValuedVerdict.LIKELY_FAIL) == Verdict.FAIL

    def test_fail_maps_to_fail(self):
        assert to_three_valued(FiveValuedVerdict.FAIL) == Verdict.FAIL

    def test_pass_is_only_source_of_three_valued_pass(self):
        """PASS is the ONLY 5-valued verdict that maps to 3-valued PASS."""
        for v in FiveValuedVerdict:
            if v != FiveValuedVerdict.PASS:
                assert to_three_valued(v) != Verdict.PASS, (
                    f"{v} should not map to Verdict.PASS"
                )

    def test_fail_and_likely_fail_are_only_sources_of_three_valued_fail(self):
        """Only FAIL and LIKELY_FAIL map to 3-valued FAIL."""
        fail_sources = [
            v for v in FiveValuedVerdict
            if to_three_valued(v) == Verdict.FAIL
        ]
        assert set(fail_sources) == {FiveValuedVerdict.FAIL, FiveValuedVerdict.LIKELY_FAIL}

    def test_uncertain_sources(self):
        """UNCERTAIN and LIKELY_PASS map to 3-valued UNCERTAIN."""
        uncertain_sources = [
            v for v in FiveValuedVerdict
            if to_three_valued(v) == Verdict.UNCERTAIN
        ]
        assert set(uncertain_sources) == {
            FiveValuedVerdict.UNCERTAIN,
            FiveValuedVerdict.LIKELY_PASS,
        }

    def test_refinement_is_surjective(self):
        """Every 3-valued verdict is reachable via refinement."""
        three_valued_results = {to_three_valued(v) for v in FiveValuedVerdict}
        assert Verdict.PASS in three_valued_results
        assert Verdict.UNCERTAIN in three_valued_results
        assert Verdict.FAIL in three_valued_results


# ═══════════════════════════════════════════════════════════════════════
# Confidence Score Tests
# ═══════════════════════════════════════════════════════════════════════

class TestConfidenceScore:
    """Test confidence scores for all verdicts."""

    def test_pass_confidence(self):
        assert confidence_score(FiveValuedVerdict.PASS) == 1.0

    def test_likely_pass_confidence(self):
        assert confidence_score(FiveValuedVerdict.LIKELY_PASS) == 0.75

    def test_uncertain_confidence(self):
        assert confidence_score(FiveValuedVerdict.UNCERTAIN) == 0.5

    def test_likely_fail_confidence(self):
        assert confidence_score(FiveValuedVerdict.LIKELY_FAIL) == 0.25

    def test_fail_confidence(self):
        assert confidence_score(FiveValuedVerdict.FAIL) == 0.0

    def test_confidence_is_monotonically_decreasing(self):
        """Confidence scores decrease with the ordering."""
        order = [
            FiveValuedVerdict.PASS,
            FiveValuedVerdict.LIKELY_PASS,
            FiveValuedVerdict.UNCERTAIN,
            FiveValuedVerdict.LIKELY_FAIL,
            FiveValuedVerdict.FAIL,
        ]
        for i in range(len(order) - 1):
            assert confidence_score(order[i]) > confidence_score(order[i + 1]), (
                f"confidence({order[i]}) should be > confidence({order[i+1]})"
            )

    def test_confidence_range(self):
        """All confidence scores are in [0.0, 1.0]."""
        for v in FiveValuedVerdict:
            score = confidence_score(v)
            assert 0.0 <= score <= 1.0

    def test_confidence_matches_verdict_confidence_property(self):
        """Confidence scores in FiveValuedVerdict match those in Verdict for shared values."""
        assert confidence_score(FiveValuedVerdict.PASS) == Verdict.PASS.confidence
        assert confidence_score(FiveValuedVerdict.UNCERTAIN) == Verdict.UNCERTAIN.confidence
        assert confidence_score(FiveValuedVerdict.FAIL) == Verdict.FAIL.confidence


# ═══════════════════════════════════════════════════════════════════════
# combine_verdicts Lattice Combination Tests
# ═══════════════════════════════════════════════════════════════════════

class TestCombineVerdicts:
    """Test the lattice combination of multiple verdicts."""

    def test_empty_list_returns_uncertain(self):
        assert combine_verdicts([]) == FiveValuedVerdict.UNCERTAIN

    def test_single_verdict(self):
        for v in FiveValuedVerdict:
            assert combine_verdicts([v]) == v

    def test_all_pass_yields_pass(self):
        assert combine_verdicts([
            FiveValuedVerdict.PASS,
            FiveValuedVerdict.PASS,
            FiveValuedVerdict.PASS,
        ]) == FiveValuedVerdict.PASS

    def test_any_fail_yields_fail(self):
        assert combine_verdicts([
            FiveValuedVerdict.PASS,
            FiveValuedVerdict.FAIL,
            FiveValuedVerdict.PASS,
        ]) == FiveValuedVerdict.FAIL

    def test_any_likely_fail_yields_likely_fail_or_worse(self):
        """LIKELY_FAIL is the weakest unless there's an actual FAIL."""
        result = combine_verdicts([
            FiveValuedVerdict.PASS,
            FiveValuedVerdict.LIKELY_FAIL,
        ])
        assert result in (FiveValuedVerdict.LIKELY_FAIL, FiveValuedVerdict.FAIL)
        assert result == FiveValuedVerdict.LIKELY_FAIL  # No actual FAIL present

    def test_pass_and_uncertain_yields_uncertain(self):
        assert combine_verdicts([
            FiveValuedVerdict.PASS,
            FiveValuedVerdict.UNCERTAIN,
        ]) == FiveValuedVerdict.UNCERTAIN

    def test_pass_and_likely_pass_yields_likely_pass(self):
        assert combine_verdicts([
            FiveValuedVerdict.PASS,
            FiveValuedVerdict.LIKELY_PASS,
        ]) == FiveValuedVerdict.LIKELY_PASS

    def test_likely_fail_and_fail_yields_fail(self):
        assert combine_verdicts([
            FiveValuedVerdict.LIKELY_FAIL,
            FiveValuedVerdict.FAIL,
        ]) == FiveValuedVerdict.FAIL

    def test_weakest_link_principle(self):
        """The combined verdict is the weakest (lowest-ordered) individual verdict."""
        # FAIL is weakest
        all_verdicts = list(FiveValuedVerdict)
        assert combine_verdicts(all_verdicts) == FiveValuedVerdict.FAIL

    def test_likely_pass_and_uncertain_yields_uncertain(self):
        assert combine_verdicts([
            FiveValuedVerdict.LIKELY_PASS,
            FiveValuedVerdict.UNCERTAIN,
        ]) == FiveValuedVerdict.UNCERTAIN

    def test_pairwise_combination_is_commutative(self):
        """combine_verdicts([a, b]) == combine_verdicts([b, a])."""
        for a in FiveValuedVerdict:
            for b in FiveValuedVerdict:
                assert combine_verdicts([a, b]) == combine_verdicts([b, a])

    def test_combination_preserves_three_valued_soundness(self):
        """Combining 5-valued verdicts then refining should be ≤
        refining each then combining in 3-valued logic."""
        for a in FiveValuedVerdict:
            for b in FiveValuedVerdict:
                combined = combine_verdicts([a, b])
                refined_combined = to_three_valued(combined)
                # The refinement of the combined result should match
                # the 3-valued combination of the refined individual results
                a3 = to_three_valued(a)
                b3 = to_three_valued(b)
                three_valued_combined = a3  # Start with first
                for v in [b3]:
                    # 3-valued AND: minimum ordering
                    from biocompiler.types import _VERDICT_ORDER
                    if _VERDICT_ORDER[v] < _VERDICT_ORDER[three_valued_combined]:
                        three_valued_combined = v
                assert refined_combined == three_valued_combined


# ═══════════════════════════════════════════════════════════════════════
# verify_five_valued_soundness Tests
# ═══════════════════════════════════════════════════════════════════════

class TestVerifyFiveValuedSoundness:
    """Test the soundness verification function."""

    def test_pass_with_true_condition_is_sound(self):
        assert verify_five_valued_soundness(FiveValuedVerdict.PASS, True) is True

    def test_pass_with_false_condition_is_unsound(self):
        assert verify_five_valued_soundness(FiveValuedVerdict.PASS, False) is False

    def test_likely_pass_with_true_condition_is_sound(self):
        assert verify_five_valued_soundness(FiveValuedVerdict.LIKELY_PASS, True) is True

    def test_likely_pass_with_false_condition_is_unsound(self):
        assert verify_five_valued_soundness(FiveValuedVerdict.LIKELY_PASS, False) is False

    def test_uncertain_always_sound(self):
        """UNCERTAIN makes no definite claim, so it's always sound."""
        assert verify_five_valued_soundness(FiveValuedVerdict.UNCERTAIN, True) is True
        assert verify_five_valued_soundness(FiveValuedVerdict.UNCERTAIN, False) is True

    def test_likely_fail_with_false_condition_is_sound(self):
        assert verify_five_valued_soundness(FiveValuedVerdict.LIKELY_FAIL, False) is True

    def test_likely_fail_with_true_condition_is_unsound(self):
        assert verify_five_valued_soundness(FiveValuedVerdict.LIKELY_FAIL, True) is False

    def test_fail_with_false_condition_is_sound(self):
        assert verify_five_valued_soundness(FiveValuedVerdict.FAIL, False) is True

    def test_fail_with_true_condition_is_unsound(self):
        assert verify_five_valued_soundness(FiveValuedVerdict.FAIL, True) is False

    @pytest.mark.parametrize("verdict,condition,expected", [
        (FiveValuedVerdict.PASS, True, True),
        (FiveValuedVerdict.PASS, False, False),
        (FiveValuedVerdict.LIKELY_PASS, True, True),
        (FiveValuedVerdict.LIKELY_PASS, False, False),
        (FiveValuedVerdict.UNCERTAIN, True, True),
        (FiveValuedVerdict.UNCERTAIN, False, True),
        (FiveValuedVerdict.LIKELY_FAIL, False, True),
        (FiveValuedVerdict.LIKELY_FAIL, True, False),
        (FiveValuedVerdict.FAIL, False, True),
        (FiveValuedVerdict.FAIL, True, False),
    ])
    def test_parametrized_soundness(self, verdict, condition, expected):
        """Parametrized test for all verdict/condition combinations."""
        assert verify_five_valued_soundness(verdict, condition) is expected


# ═══════════════════════════════════════════════════════════════════════
# five_valued_and / five_valued_or Tests
# ═══════════════════════════════════════════════════════════════════════

class TestFiveValuedAndOr:
    """Test 5-valued AND and OR operations."""

    def test_and_idempotent(self):
        for v in FiveValuedVerdict:
            assert five_valued_and(v, v) == v

    def test_or_idempotent(self):
        for v in FiveValuedVerdict:
            assert five_valued_or(v, v) == v

    def test_and_commutative(self):
        for a in FiveValuedVerdict:
            for b in FiveValuedVerdict:
                assert five_valued_and(a, b) == five_valued_and(b, a)

    def test_or_commutative(self):
        for a in FiveValuedVerdict:
            for b in FiveValuedVerdict:
                assert five_valued_or(a, b) == five_valued_or(b, a)

    def test_and_pass_is_identity(self):
        for v in FiveValuedVerdict:
            assert five_valued_and(FiveValuedVerdict.PASS, v) == v

    def test_or_fail_is_identity(self):
        for v in FiveValuedVerdict:
            assert five_valued_or(FiveValuedVerdict.FAIL, v) == v

    def test_and_fail_is_absorbing(self):
        for v in FiveValuedVerdict:
            assert five_valued_and(FiveValuedVerdict.FAIL, v) == FiveValuedVerdict.FAIL

    def test_or_pass_is_absorbing(self):
        for v in FiveValuedVerdict:
            assert five_valued_or(FiveValuedVerdict.PASS, v) == FiveValuedVerdict.PASS

    def test_and_returns_minimum(self):
        """AND returns the verdict with lower ordering."""
        order = [
            FiveValuedVerdict.FAIL,
            FiveValuedVerdict.LIKELY_FAIL,
            FiveValuedVerdict.UNCERTAIN,
            FiveValuedVerdict.LIKELY_PASS,
            FiveValuedVerdict.PASS,
        ]
        for a in order:
            for b in order:
                result = five_valued_and(a, b)
                # Result should be the lower-ordered of the two
                idx_a = order.index(a)
                idx_b = order.index(b)
                expected = order[min(idx_a, idx_b)]
                assert result == expected, (
                    f"five_valued_and({a}, {b}) = {result}, expected {expected}"
                )

    def test_or_returns_maximum(self):
        """OR returns the verdict with higher ordering."""
        order = [
            FiveValuedVerdict.FAIL,
            FiveValuedVerdict.LIKELY_FAIL,
            FiveValuedVerdict.UNCERTAIN,
            FiveValuedVerdict.LIKELY_PASS,
            FiveValuedVerdict.PASS,
        ]
        for a in order:
            for b in order:
                result = five_valued_or(a, b)
                idx_a = order.index(a)
                idx_b = order.index(b)
                expected = order[max(idx_a, idx_b)]
                assert result == expected, (
                    f"five_valued_or({a}, {b}) = {result}, expected {expected}"
                )

    def test_and_delegates_to_combine_verdicts(self):
        """five_valued_and is consistent with combine_verdicts for two elements."""
        for a in FiveValuedVerdict:
            for b in FiveValuedVerdict:
                assert five_valued_and(a, b) == combine_verdicts([a, b])


# ═══════════════════════════════════════════════════════════════════════
# refinement_is_sound Tests
# ═══════════════════════════════════════════════════════════════════════

class TestRefinementIsSound:
    """Test that the refinement to 3-valued logic preserves soundness."""

    def test_pass_refinement_sound_when_true(self):
        assert refinement_is_sound(FiveValuedVerdict.PASS, True) is True

    def test_pass_refinement_unsound_when_false(self):
        assert refinement_is_sound(FiveValuedVerdict.PASS, False) is False

    def test_likely_pass_refinement_always_sound(self):
        """LIKELY_PASS → UNCERTAIN: no definite claim in 3-valued model."""
        assert refinement_is_sound(FiveValuedVerdict.LIKELY_PASS, True) is True
        assert refinement_is_sound(FiveValuedVerdict.LIKELY_PASS, False) is True

    def test_uncertain_refinement_always_sound(self):
        """UNCERTAIN → UNCERTAIN: no definite claim."""
        assert refinement_is_sound(FiveValuedVerdict.UNCERTAIN, True) is True
        assert refinement_is_sound(FiveValuedVerdict.UNCERTAIN, False) is True

    def test_likely_fail_refinement_sound_when_false(self):
        """LIKELY_FAIL → FAIL: requires actual_condition = False."""
        assert refinement_is_sound(FiveValuedVerdict.LIKELY_FAIL, False) is True

    def test_likely_fail_refinement_unsound_when_true(self):
        assert refinement_is_sound(FiveValuedVerdict.LIKELY_FAIL, True) is False

    def test_fail_refinement_sound_when_false(self):
        assert refinement_is_sound(FiveValuedVerdict.FAIL, False) is True

    def test_fail_refinement_unsound_when_true(self):
        assert refinement_is_sound(FiveValuedVerdict.FAIL, True) is False


# ═══════════════════════════════════════════════════════════════════════
# Consistency with Existing Verdict Enum
# ═══════════════════════════════════════════════════════════════════════

class TestConsistencyWithThreeValued:
    """Test consistency between 5-valued and existing 3-valued Verdict."""

    def test_shared_values_match(self):
        """PASS/UNCERTAIN/FAIL values should have matching string representations."""
        assert FiveValuedVerdict.PASS.value == Verdict.PASS.value
        assert FiveValuedVerdict.UNCERTAIN.value == Verdict.UNCERTAIN.value
        assert FiveValuedVerdict.FAIL.value == Verdict.FAIL.value

    def test_refinement_is_identity_for_three_common_values(self):
        """to_three_valued is identity-like for PASS, UNCERTAIN, FAIL."""
        assert to_three_valued(FiveValuedVerdict.PASS) == Verdict.PASS
        assert to_three_valued(FiveValuedVerdict.UNCERTAIN) == Verdict.UNCERTAIN
        assert to_three_valued(FiveValuedVerdict.FAIL) == Verdict.FAIL

    def test_confidence_consistency(self):
        """Confidence scores are consistent with Verdict.confidence for shared values."""
        assert confidence_score(FiveValuedVerdict.PASS) == Verdict.PASS.confidence
        assert confidence_score(FiveValuedVerdict.UNCERTAIN) == Verdict.UNCERTAIN.confidence
        assert confidence_score(FiveValuedVerdict.FAIL) == Verdict.FAIL.confidence

    def test_likely_pass_confidence_between_pass_and_uncertain(self):
        """LIKELY_PASS confidence is between PASS and UNCERTAIN."""
        assert Verdict.PASS.confidence > confidence_score(FiveValuedVerdict.LIKELY_PASS)
        assert confidence_score(FiveValuedVerdict.LIKELY_PASS) > Verdict.UNCERTAIN.confidence

    def test_likely_fail_confidence_between_uncertain_and_fail(self):
        """LIKELY_FAIL confidence is between UNCERTAIN and FAIL."""
        assert Verdict.UNCERTAIN.confidence > confidence_score(FiveValuedVerdict.LIKELY_FAIL)
        assert confidence_score(FiveValuedVerdict.LIKELY_FAIL) > Verdict.FAIL.confidence


# ═══════════════════════════════════════════════════════════════════════
# Edge Cases
# ═══════════════════════════════════════════════════════════════════════

class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_combine_many_verdicts(self):
        """Combining many verdicts still yields the minimum."""
        verdicts = [FiveValuedVerdict.PASS] * 100 + [FiveValuedVerdict.UNCERTAIN]
        assert combine_verdicts(verdicts) == FiveValuedVerdict.UNCERTAIN

    def test_combine_all_same(self):
        """Combining identical verdicts yields that verdict."""
        for v in FiveValuedVerdict:
            assert combine_verdicts([v] * 10) == v

    def test_combine_single_fail_among_passes(self):
        """A single FAIL dominates all PASSes."""
        verdicts = [FiveValuedVerdict.PASS] * 99 + [FiveValuedVerdict.FAIL]
        assert combine_verdicts(verdicts) == FiveValuedVerdict.FAIL

    def test_confidence_sum_is_not_one(self):
        """Confidence scores don't sum to 1.0 (they're not probabilities)."""
        total = sum(confidence_score(v) for v in FiveValuedVerdict)
        assert total != 1.0

    def test_round_trip_identity_for_pass(self):
        """PASS is the only verdict where refinement → identity would hold."""
        v = FiveValuedVerdict.PASS
        v3 = to_three_valued(v)
        assert v3 == Verdict.PASS
        # Going back: PASS in 3-valued means PASS in 5-valued
        # (by toThreeValued_preserves_PASS)
