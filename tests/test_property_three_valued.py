"""
Property-based tests bridging Lean4 ThreeValued theorems and Python implementation.

The Lean4 proofs in proof/BioCompiler/ThreeValued.lean verify a simplified
three-valued model. These Hypothesis-based property tests verify that the
actual Python five-valued Kleene logic implementation satisfies the same
algebraic properties, closing the proof-implementation gap (Caveat 1).

Lean4 theorems covered:
  - and_PASS          → identity for AND
  - and_FAIL_left     → sticky FAIL (left)
  - and_FAIL_right    → sticky FAIL (right)
  - and_comm          → commutativity of AND
  - or_comm           → commutativity of OR
  - and_assoc         → associativity of AND
  - or_assoc          → associativity of OR
  - and_idem          → idempotence of AND
  - or_idem           → idempotence of OR
  - and_eq_PASS_iff   → AND(a,b)=PASS iff a=PASS and b=PASS
  - and_pass_pass     → both PASS implies AND = PASS
  - foldl_ne_pass_of_ne_pass  → combined_verdict from non-PASS never reaches PASS
  - foldl_and_pass_implies_all_pass → combined PASS implies all PASS
  - and_pass_all_pass → combined PASS implies no FAIL
  - and_monotone_left → monotonicity of AND
  - ordering_pass_highest / ordering_fail_lowest
  - (distributivity — lattice property, not in Lean4 but follows from min/max)
"""

import pytest
pytest.importorskip("hypothesis")
from hypothesis import given, strategies as st, assume, settings, Verbosity
from hypothesis.stateful import RuleBasedStateMachine, rule, invariant

from biocompiler.shared.types import (
    Verdict,
    five_valued_and,
    five_valued_or,
    combined_verdict,
)
from biocompiler.shared.types import _VERDICT_ORDER
from biocompiler.type_system import SpliceVerdict


# ────────────────────────────────────────────────────────────
# Strategies
# ────────────────────────────────────────────────────────────

verdict = st.sampled_from(list(Verdict))

# Three-valued subset (matching Lean4 model)
three_valued = st.sampled_from([Verdict.PASS, Verdict.UNCERTAIN, Verdict.FAIL])

verdict_lists = st.lists(verdict, min_size=0, max_size=20)
non_empty_verdict_lists = st.lists(verdict, min_size=1, max_size=20)

splice_verdict = st.sampled_from(list(SpliceVerdict))


# ────────────────────────────────────────────────────────────
# Kleene conjunction properties  (Lean4: and_PASS, and_FAIL_left/right)
# ────────────────────────────────────────────────────────────

class TestKleeneConjunction:
    """AND(PASS, x) = x  and  AND(FAIL, x) = FAIL  for all x."""

    @given(v=verdict)
    def test_and_pass_identity(self, v):
        """Lean4: and_PASS — Verdict.and PASS v = v"""
        assert five_valued_and(Verdict.PASS, v) == v

    @given(v=verdict)
    def test_and_pass_identity_right(self, v):
        """Lean4: and_PASS (right side) — Verdict.and v PASS = v"""
        assert five_valued_and(v, Verdict.PASS) == v

    @given(v=verdict)
    def test_and_fail_left(self, v):
        """Lean4: and_FAIL_left — Verdict.and FAIL v = FAIL"""
        assert five_valued_and(Verdict.FAIL, v) == Verdict.FAIL

    @given(v=verdict)
    def test_and_fail_right(self, v):
        """Lean4: and_FAIL_right — Verdict.and v FAIL = FAIL"""
        assert five_valued_and(v, Verdict.FAIL) == Verdict.FAIL

    @given(a=verdict, b=verdict)
    def test_and_pass_pass_yields_pass(self, a, b):
        """Lean4: and_pass_pass — if a=PASS and b=PASS then AND(a,b)=PASS"""
        assume(a == Verdict.PASS and b == Verdict.PASS)
        assert five_valued_and(a, b) == Verdict.PASS

    @given(v=verdict)
    def test_and_fail_fail_yields_fail(self, v):
        """AND(FAIL, FAIL) = FAIL (specific instance)"""
        assert five_valued_and(Verdict.FAIL, Verdict.FAIL) == Verdict.FAIL


# ────────────────────────────────────────────────────────────
# Kleene disjunction properties  (Lean4: or analogues)
# ────────────────────────────────────────────────────────────

class TestKleeneDisjunction:
    """OR(FAIL, x) = x  and  OR(PASS, x) = PASS  for all x."""

    @given(v=verdict)
    def test_or_fail_identity(self, v):
        """OR(FAIL, x) = x  (FAIL is identity for OR)"""
        assert five_valued_or(Verdict.FAIL, v) == v

    @given(v=verdict)
    def test_or_fail_identity_right(self, v):
        """OR(x, FAIL) = x  (FAIL is identity for OR)"""
        assert five_valued_or(v, Verdict.FAIL) == v

    @given(v=verdict)
    def test_or_pass_left(self, v):
        """OR(PASS, x) = PASS  (PASS annihilates OR)"""
        assert five_valued_or(Verdict.PASS, v) == Verdict.PASS

    @given(v=verdict)
    def test_or_pass_right(self, v):
        """OR(x, PASS) = PASS  (PASS annihilates OR)"""
        assert five_valued_or(v, Verdict.PASS) == Verdict.PASS

    @given(v=verdict)
    def test_or_fail_fail_yields_fail(self, v):
        """OR(FAIL, FAIL) = FAIL"""
        assert five_valued_or(Verdict.FAIL, Verdict.FAIL) == Verdict.FAIL


# ────────────────────────────────────────────────────────────
# Commutativity  (Lean4: and_comm, or_comm)
# ────────────────────────────────────────────────────────────

class TestCommutativity:
    @given(a=verdict, b=verdict)
    def test_and_commutative(self, a, b):
        """Lean4: and_comm — AND(a,b) = AND(b,a)"""
        assert five_valued_and(a, b) == five_valued_and(b, a)

    @given(a=verdict, b=verdict)
    def test_or_commutative(self, a, b):
        """Lean4: or_comm — OR(a,b) = OR(b,a)"""
        assert five_valued_or(a, b) == five_valued_or(b, a)


# ────────────────────────────────────────────────────────────
# Associativity  (Lean4: and_assoc, or_assoc)
# ────────────────────────────────────────────────────────────

class TestAssociativity:
    @given(a=verdict, b=verdict, c=verdict)
    def test_and_associative(self, a, b, c):
        """Lean4: and_assoc — AND(AND(a,b),c) = AND(a,AND(b,c))"""
        assert five_valued_and(five_valued_and(a, b), c) == five_valued_and(a, five_valued_and(b, c))

    @given(a=verdict, b=verdict, c=verdict)
    def test_or_associative(self, a, b, c):
        """Lean4: or_assoc — OR(OR(a,b),c) = OR(a,OR(b,c))"""
        assert five_valued_or(five_valued_or(a, b), c) == five_valued_or(a, five_valued_or(b, c))


# ────────────────────────────────────────────────────────────
# Distributivity  (lattice property: AND distributes over OR and vice versa)
# ────────────────────────────────────────────────────────────

class TestDistributivity:
    @given(a=verdict, b=verdict, c=verdict)
    def test_and_distributes_over_or(self, a, b, c):
        """AND(a, OR(b,c)) = OR(AND(a,b), AND(a,c))"""
        lhs = five_valued_and(a, five_valued_or(b, c))
        rhs = five_valued_or(five_valued_and(a, b), five_valued_and(a, c))
        assert lhs == rhs

    @given(a=verdict, b=verdict, c=verdict)
    def test_or_distributes_over_and(self, a, b, c):
        """OR(a, AND(b,c)) = AND(OR(a,b), OR(a,c))"""
        lhs = five_valued_or(a, five_valued_and(b, c))
        rhs = five_valued_and(five_valued_or(a, b), five_valued_or(a, c))
        assert lhs == rhs


# ────────────────────────────────────────────────────────────
# Idempotence  (Lean4: and_idem, or_idem)
# ────────────────────────────────────────────────────────────

class TestIdempotence:
    @given(v=verdict)
    def test_and_idempotent(self, v):
        """Lean4: and_idem — AND(v,v) = v"""
        assert five_valued_and(v, v) == v

    @given(v=verdict)
    def test_or_idempotent(self, v):
        """Lean4: or_idem — OR(v,v) = v"""
        assert five_valued_or(v, v) == v


# ────────────────────────────────────────────────────────────
# Identity elements  (Lean4: and_PASS; FAIL for OR)
# ────────────────────────────────────────────────────────────

class TestIdentity:
    @given(v=verdict)
    def test_and_pass_identity_element(self, v):
        """PASS is the identity element for AND: AND(PASS, x) = x"""
        assert five_valued_and(Verdict.PASS, v) == v

    @given(v=verdict)
    def test_or_fail_identity_element(self, v):
        """FAIL is the identity element for OR: OR(FAIL, x) = x"""
        assert five_valued_or(Verdict.FAIL, v) == v


# ────────────────────────────────────────────────────────────
# and_eq_PASS_iff  (Lean4: AND(a,b) = PASS ↔ a = PASS ∧ b = PASS)
# ────────────────────────────────────────────────────────────

class TestAndEqPassIff:
    @given(a=verdict, b=verdict)
    def test_and_eq_pass_iff(self, a, b):
        """Lean4: and_eq_PASS_iff — AND(a,b) = PASS iff a = PASS and b = PASS"""
        and_result = five_valued_and(a, b)
        # Forward: AND(a,b)=PASS => a=PASS and b=PASS
        if and_result == Verdict.PASS:
            assert a == Verdict.PASS and b == Verdict.PASS
        # Backward: a=PASS and b=PASS => AND(a,b)=PASS
        if a == Verdict.PASS and b == Verdict.PASS:
            assert and_result == Verdict.PASS


# ────────────────────────────────────────────────────────────
# combined_verdict properties  (Lean4: foldl theorems)
# ────────────────────────────────────────────────────────────

class TestCombinedVerdict:
    """Properties of combined_verdict (foldl with five_valued_and)."""

    @given(vs=verdict_lists)
    def test_combined_verdict_any_fail_yields_fail(self, vs):
        """If any element is FAIL, combined_verdict returns FAIL.
        This is the 'sticky FAIL' property (Lean4: and_FAIL_left/right)."""
        if Verdict.FAIL in vs:
            assert combined_verdict(vs) == Verdict.FAIL

    @given(vs=non_empty_verdict_lists)
    def test_combined_verdict_all_pass_yields_pass(self, vs):
        """If all elements are PASS, combined_verdict returns PASS.
        (Lean4: and_pass_pass)"""
        if all(v == Verdict.PASS for v in vs):
            assert combined_verdict(vs) == Verdict.PASS

    @given(vs=verdict_lists)
    def test_combined_verdict_pass_implies_all_pass(self, vs):
        """Lean4: foldl_and_pass_implies_all_pass — if combined = PASS,
        then every element must be PASS."""
        if combined_verdict(vs) == Verdict.PASS:
            for v in vs:
                assert v == Verdict.PASS

    @given(vs=verdict_lists)
    def test_combined_verdict_pass_implies_no_fail(self, vs):
        """Lean4: and_pass_all_pass — if combined = PASS,
        then no element is FAIL."""
        if combined_verdict(vs) == Verdict.PASS:
            assert Verdict.FAIL not in vs

    @given(vs=verdict_lists)
    def test_combined_empty_yields_uncertain(self, vs):
        """combined_verdict([]) = UNCERTAIN (convention for no evidence)."""
        assert combined_verdict([]) == Verdict.UNCERTAIN

    @given(vs=non_empty_verdict_lists)
    def test_combined_verdict_non_pass_init_never_pass(self, vs):
        """Lean4: foldl_ne_pass_of_ne_pass — folding from non-PASS
        never reaches PASS.

        If the first element is not PASS, the combined verdict cannot be PASS.
        """
        if vs[0] != Verdict.PASS:
            assert combined_verdict(vs) != Verdict.PASS

    @given(vs=verdict_lists)
    def test_combined_verdict_order_independent(self, vs):
        """Combined verdict is order-independent (commutativity + associativity)."""
        import random
        vs_shuffled = list(vs)
        random.shuffle(vs_shuffled)
        assert combined_verdict(vs) == combined_verdict(vs_shuffled)

    @given(vs=verdict_lists)
    def test_combined_verdict_worst_link(self, vs):
        """Combined verdict equals the minimum element by ordering.
        This is the 'weakest link' principle."""
        if not vs:
            return  # empty case handled separately
        min_verdict = min(vs, key=lambda v: _VERDICT_ORDER[v])
        assert combined_verdict(vs) == min_verdict


# ────────────────────────────────────────────────────────────
# Verdict ordering consistency  (Lean4: ordering, ordering_pass_highest, etc.)
# ────────────────────────────────────────────────────────────

class TestVerdictOrdering:
    """Verify the _VERDICT_ORDER mapping is consistent with the lattice."""

    @given(v=verdict)
    def test_pass_is_highest(self, v):
        """Lean4: ordering_pass_highest — PASS has the highest order."""
        assert _VERDICT_ORDER[Verdict.PASS] >= _VERDICT_ORDER[v]

    @given(v=verdict)
    def test_fail_is_lowest(self, v):
        """Lean4: ordering_fail_lowest — FAIL has the lowest order."""
        assert _VERDICT_ORDER[Verdict.FAIL] <= _VERDICT_ORDER[v]

    @given(v=verdict)
    def test_confidence_consistent_with_order(self, v):
        """Confidence score is monotonically consistent with _VERDICT_ORDER."""
        # Higher order => higher confidence
        # We check: PASS confidence > LIKELY_PASS > ... > FAIL
        for v1 in Verdict:
            for v2 in Verdict:
                if _VERDICT_ORDER[v1] > _VERDICT_ORDER[v2]:
                    assert v1.confidence > v2.confidence
                elif _VERDICT_ORDER[v1] == _VERDICT_ORDER[v2]:
                    assert v1.confidence == v2.confidence

    @given(a=verdict, b=verdict)
    def test_and_monotone_left(self, a, b):
        """Lean4: and_monotone_left — if order(a) <= order(b),
        then order(AND(a,v)) <= order(AND(b,v)) for any v."""
        v = Verdict.UNCERTAIN  # test with a fixed intermediate value
        if _VERDICT_ORDER[a] <= _VERDICT_ORDER[b]:
            assert _VERDICT_ORDER[five_valued_and(a, v)] <= _VERDICT_ORDER[five_valued_and(b, v)]

    @given(a=verdict, b=verdict, v=verdict)
    def test_and_monotone_left_general(self, a, b, v):
        """Lean4: and_monotone_left (full) — monotonicity for all v."""
        if _VERDICT_ORDER[a] <= _VERDICT_ORDER[b]:
            assert _VERDICT_ORDER[five_valued_and(a, v)] <= _VERDICT_ORDER[five_valued_and(b, v)]

    @given(a=verdict, b=verdict, v=verdict)
    def test_and_monotone_right(self, a, b, v):
        """Right monotonicity of AND (follows from commutativity)."""
        if _VERDICT_ORDER[a] <= _VERDICT_ORDER[b]:
            assert _VERDICT_ORDER[five_valued_and(v, a)] <= _VERDICT_ORDER[five_valued_and(v, b)]

    @given(v=verdict)
    def test_ordering_strictly_graded(self, v):
        """All five verdict values have distinct order values."""
        orders = [_VERDICT_ORDER[x] for x in Verdict]
        assert len(set(orders)) == len(list(Verdict))

    @given(v=verdict)
    def test_is_definite_consistent(self, v):
        """is_definite is True only for PASS and FAIL (endpoints)."""
        if v.is_definite:
            assert v in (Verdict.PASS, Verdict.FAIL)
        else:
            assert v in (Verdict.LIKELY_PASS, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL)


# ────────────────────────────────────────────────────────────
# Three-valued subset consistency  (Lean4 proofs are for 3 values)
# ────────────────────────────────────────────────────────────

class TestThreeValuedSubset:
    """Verify that the three-valued subset {PASS, UNCERTAIN, FAIL} satisfies
    all Lean4 theorems exactly as stated, when restricted to this subset."""

    @given(a=three_valued, b=three_valued)
    def test_three_valued_and_commutative(self, a, b):
        """Lean4: and_comm restricted to three values."""
        assert five_valued_and(a, b) == five_valued_and(b, a)

    @given(a=three_valued, b=three_valued, c=three_valued)
    def test_three_valued_and_associative(self, a, b, c):
        """Lean4: and_assoc restricted to three values."""
        assert five_valued_and(five_valued_and(a, b), c) == five_valued_and(a, five_valued_and(b, c))

    @given(a=three_valued, b=three_valued)
    def test_three_valued_or_commutative(self, a, b):
        """Lean4: or_comm restricted to three values."""
        assert five_valued_or(a, b) == five_valued_or(b, a)

    @given(a=three_valued, b=three_valued, c=three_valued)
    def test_three_valued_or_associative(self, a, b, c):
        """Lean4: or_assoc restricted to three values."""
        assert five_valued_or(five_valued_or(a, b), c) == five_valued_or(a, five_valued_or(b, c))

    @given(v=three_valued)
    def test_three_valued_and_idempotent(self, v):
        """Lean4: and_idem restricted to three values."""
        assert five_valued_and(v, v) == v

    @given(v=three_valued)
    def test_three_valued_or_idempotent(self, v):
        """Lean4: or_idem restricted to three values."""
        assert five_valued_or(v, v) == v

    @given(v=three_valued)
    def test_three_valued_and_pass_identity(self, v):
        """Lean4: and_PASS restricted to three values."""
        assert five_valued_and(Verdict.PASS, v) == v

    @given(v=three_valued)
    def test_three_valued_and_fail_absorbs(self, v):
        """Lean4: and_FAIL_left restricted to three values."""
        assert five_valued_and(Verdict.FAIL, v) == Verdict.FAIL

    @given(a=three_valued, b=three_valued)
    def test_three_valued_and_eq_pass_iff(self, a, b):
        """Lean4: and_eq_PASS_iff restricted to three values."""
        if five_valued_and(a, b) == Verdict.PASS:
            assert a == Verdict.PASS and b == Verdict.PASS
        if a == Verdict.PASS and b == Verdict.PASS:
            assert five_valued_and(a, b) == Verdict.PASS

    @given(a=three_valued, b=three_valued, c=three_valued)
    def test_three_valued_and_distributes_over_or(self, a, b, c):
        """Distributivity restricted to three values."""
        lhs = five_valued_and(a, five_valued_or(b, c))
        rhs = five_valued_or(five_valued_and(a, b), five_valued_and(a, c))
        assert lhs == rhs

    @given(a=three_valued, b=three_valued)
    def test_three_valued_and_pass_uncertain(self, a, b):
        """Lean4: and_pass_uncertain — AND(PASS, UNCERTAIN) = UNCERTAIN."""
        assert five_valued_and(Verdict.PASS, Verdict.UNCERTAIN) == Verdict.UNCERTAIN

    @given(a=three_valued, b=three_valued)
    def test_three_valued_and_uncertain_pass(self, a, b):
        """Lean4: and_uncertain_pass — AND(UNCERTAIN, PASS) = UNCERTAIN."""
        assert five_valued_and(Verdict.UNCERTAIN, Verdict.PASS) == Verdict.UNCERTAIN

    @given(vs=st.lists(three_valued, min_size=0, max_size=15))
    def test_three_valued_foldl_uncertain_ne_pass(self, vs):
        """Lean4: foldl_uncertain_ne_pass — foldl from UNCERTAIN never reaches PASS."""
        if vs and vs[0] == Verdict.UNCERTAIN:
            assert combined_verdict(vs) != Verdict.PASS

    @given(vs=st.lists(three_valued, min_size=1, max_size=15))
    def test_three_valued_foldl_fail_ne_pass(self, vs):
        """Lean4: foldl_fail_ne_pass — foldl from FAIL never reaches PASS."""
        if vs[0] == Verdict.FAIL:
            assert combined_verdict(vs) != Verdict.PASS


# ────────────────────────────────────────────────────────────
# SpliceVerdict dual-threshold classification monotonicity
# ────────────────────────────────────────────────────────────

def _classify_splice(score: float, low: float, high: float) -> SpliceVerdict:
    """Dual-threshold classification matching check_no_cryptic_splice logic."""
    if score < low:
        return SpliceVerdict.PASS
    elif score < high:
        return SpliceVerdict.UNCERTAIN
    else:
        return SpliceVerdict.FAIL


# SpliceVerdict ordering for monotonicity checks (higher = worse)
_SPLICE_ORDER = {
    SpliceVerdict.PASS: 0,
    SpliceVerdict.UNCERTAIN: 1,
    SpliceVerdict.FAIL: 2,
}


class TestSpliceVerdictMonotonicity:
    """Verify that SpliceVerdict dual-threshold classification is monotone:
    higher scores produce verdicts that are at least as 'bad'."""

    @given(
        score1=st.floats(min_value=-50.0, max_value=50.0, allow_nan=False, allow_infinity=False),
        score2=st.floats(min_value=-50.0, max_value=50.0, allow_nan=False, allow_infinity=False),
        low=st.floats(min_value=0.1, max_value=5.0, allow_nan=False, allow_infinity=False),
        high=st.floats(min_value=5.1, max_value=15.0, allow_nan=False, allow_infinity=False),
    )
    def test_splice_classification_monotone(self, score1, score2, low, high):
        """If score1 <= score2, then classify(score1) <= classify(score2)
        in the 'badness' ordering (PASS < UNCERTAIN < FAIL)."""
        assume(low < high)
        if score1 <= score2:
            v1 = _classify_splice(score1, low, high)
            v2 = _classify_splice(score2, low, high)
            assert _SPLICE_ORDER[v1] <= _SPLICE_ORDER[v2]

    @given(
        low=st.floats(min_value=0.1, max_value=5.0, allow_nan=False, allow_infinity=False),
        high=st.floats(min_value=5.1, max_value=15.0, allow_nan=False, allow_infinity=False),
    )
    def test_splice_below_low_is_pass(self, low, high):
        """Score below low threshold always yields PASS."""
        assume(low < high)
        score = low - 1.0  # guaranteed below low
        assert _classify_splice(score, low, high) == SpliceVerdict.PASS

    @given(
        low=st.floats(min_value=0.1, max_value=5.0, allow_nan=False, allow_infinity=False),
        high=st.floats(min_value=5.1, max_value=15.0, allow_nan=False, allow_infinity=False),
    )
    def test_splice_above_high_is_fail(self, low, high):
        """Score at or above high threshold always yields FAIL."""
        assume(low < high)
        score = high + 1.0  # guaranteed above high
        assert _classify_splice(score, low, high) == SpliceVerdict.FAIL

    @given(
        low=st.floats(min_value=0.1, max_value=5.0, allow_nan=False, allow_infinity=False),
        high=st.floats(min_value=5.1, max_value=15.0, allow_nan=False, allow_infinity=False),
    )
    def test_splice_between_thresholds_is_uncertain(self, low, high):
        """Score in [low, high) always yields UNCERTAIN."""
        assume(low < high)
        score = (low + high) / 2.0  # guaranteed in [low, high)
        assert _classify_splice(score, low, high) == SpliceVerdict.UNCERTAIN

    @given(
        low=st.floats(min_value=0.1, max_value=5.0, allow_nan=False, allow_infinity=False),
        high=st.floats(min_value=5.1, max_value=15.0, allow_nan=False, allow_infinity=False),
    )
    def test_splice_three_verdicts_exhaustive(self, low, high):
        """All three SpliceVerdict values are reachable for any valid thresholds."""
        assume(low < high)
        assert _classify_splice(low - 1.0, low, high) == SpliceVerdict.PASS
        assert _classify_splice((low + high) / 2, low, high) == SpliceVerdict.UNCERTAIN
        assert _classify_splice(high + 1.0, low, high) == SpliceVerdict.FAIL


# ────────────────────────────────────────────────────────────
# Absorption laws  (lattice property)
# ────────────────────────────────────────────────────────────

class TestAbsorption:
    """Absorption: AND(a, OR(a,b)) = a  and  OR(a, AND(a,b)) = a.
    These follow from the lattice structure and are critical for
    normal-form simplifications."""

    @given(a=verdict, b=verdict)
    def test_and_absorbs_or(self, a, b):
        """AND(a, OR(a,b)) = a"""
        assert five_valued_and(a, five_valued_or(a, b)) == a

    @given(a=verdict, b=verdict)
    def test_or_absorbs_and(self, a, b):
        """OR(a, AND(a,b)) = a"""
        assert five_valued_or(a, five_valued_and(a, b)) == a


# ────────────────────────────────────────────────────────────
# De Morgan-like consistency  (NOT is not implemented in Python,
# but we verify the dual relationship between AND and OR via ordering)
# ────────────────────────────────────────────────────────────

class TestLatticeConsistency:
    """Verify that AND and OR are consistent duals under the ordering."""

    @given(a=verdict, b=verdict)
    def test_and_result_order_le_min(self, a, b):
        """AND(a,b) has order <= min(order(a), order(b)).
        (For Kleene min, AND(a,b) = argmin, so order equals the min.)"""
        result = five_valued_and(a, b)
        assert _VERDICT_ORDER[result] == min(_VERDICT_ORDER[a], _VERDICT_ORDER[b])

    @given(a=verdict, b=verdict)
    def test_or_result_order_ge_max(self, a, b):
        """OR(a,b) has order >= max(order(a), order(b)).
        (For Kleene max, OR(a,b) = argmax, so order equals the max.)"""
        result = five_valued_or(a, b)
        assert _VERDICT_ORDER[result] == max(_VERDICT_ORDER[a], _VERDICT_ORDER[b])

    @given(a=verdict, b=verdict)
    def test_and_le_operands(self, a, b):
        """AND(a,b) <= a and AND(a,b) <= b in the ordering."""
        result = five_valued_and(a, b)
        assert _VERDICT_ORDER[result] <= _VERDICT_ORDER[a]
        assert _VERDICT_ORDER[result] <= _VERDICT_ORDER[b]

    @given(a=verdict, b=verdict)
    def test_or_ge_operands(self, a, b):
        """OR(a,b) >= a and OR(a,b) >= b in the ordering."""
        result = five_valued_or(a, b)
        assert _VERDICT_ORDER[result] >= _VERDICT_ORDER[a]
        assert _VERDICT_ORDER[result] >= _VERDICT_ORDER[b]
