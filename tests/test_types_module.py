"""
Tests for biocompiler.types — the canonical core types module.

Covers:
1. SLOTMode enum values and properties
2. Verdict enum values, confidence, and is_definite
3. five_valued_and / five_valued_or logic
4. combined_verdict function
5. PositionRange dataclass
6. Token dataclass (frozen, range property)
7. SpliceIsoform dataclass (frozen, repr)
8. TypeCheckResult dataclass (passed property, repr)
9. Certificate dataclass (to_dict, from_dict, validation)
10. Backward-compatible aliases (three_valued_and, three_valued_or)
11. __all__ exports
"""

from __future__ import annotations

import pytest

from biocompiler.types import (
    SLOTMode,
    Verdict,
    _VERDICT_ORDER,
    five_valued_and,
    five_valued_or,
    three_valued_and,
    three_valued_or,
    combined_verdict,
    PositionRange,
    Token,
    SpliceIsoform,
    TypeCheckResult,
    Certificate,
    __all__ as types_all,
)


# ═══════════════════════════════════════════════════════════════════════════
# 1. SLOTMode
# ═══════════════════════════════════════════════════════════════════════════

class TestSLOTMode:
    """Test the SLOTMode enum."""

    def test_conservative_value(self):
        assert SLOTMode.CONSERVATIVE.value == "conservative"

    def test_verified_value(self):
        assert SLOTMode.VERIFIED.value == "verified"

    def test_permissive_value(self):
        assert SLOTMode.PERMISSIVE.value == "permissive"

    def test_exactly_three_members(self):
        assert len(SLOTMode) == 3

    def test_all_values_distinct(self):
        values = {m.value for m in SLOTMode}
        assert len(values) == 3

    def test_members_are_enum_instances(self):
        for m in SLOTMode:
            assert isinstance(m, SLOTMode)

    def test_value_string_comparison(self):
        """SLOTMode members have .value that compares equal to the string."""
        assert SLOTMode.CONSERVATIVE.value == "conservative"
        assert SLOTMode.VERIFIED.value == "verified"
        assert SLOTMode.PERMISSIVE.value == "permissive"

    def test_iteration_order(self):
        """Iterating yields CONSERVATIVE, VERIFIED, PERMISSIVE (definition order)."""
        members = list(SLOTMode)
        assert members == [
            SLOTMode.CONSERVATIVE,
            SLOTMode.VERIFIED,
            SLOTMode.PERMISSIVE,
        ]


# ═══════════════════════════════════════════════════════════════════════════
# 2. Verdict
# ═══════════════════════════════════════════════════════════════════════════

class TestVerdictEnum:
    """Test the Verdict five-valued enum."""

    def test_pass_value(self):
        assert Verdict.PASS.value == "PASS"

    def test_likely_pass_value(self):
        assert Verdict.LIKELY_PASS.value == "LIKELY_PASS"

    def test_uncertain_value(self):
        assert Verdict.UNCERTAIN.value == "UNCERTAIN"

    def test_likely_fail_value(self):
        assert Verdict.LIKELY_FAIL.value == "LIKELY_FAIL"

    def test_fail_value(self):
        assert Verdict.FAIL.value == "FAIL"

    def test_exactly_five_members(self):
        assert len(Verdict) == 5

    def test_all_values_distinct(self):
        values = {v.value for v in Verdict}
        assert len(values) == 5

    def test_string_comparison(self):
        """Verdict inherits from str+Enum, so members compare equal to their values."""
        assert Verdict.PASS == "PASS"
        assert Verdict.FAIL == "FAIL"
        assert Verdict.UNCERTAIN == "UNCERTAIN"

    def test_iteration_order(self):
        members = list(Verdict)
        assert members == [
            Verdict.PASS,
            Verdict.LIKELY_PASS,
            Verdict.UNCERTAIN,
            Verdict.LIKELY_FAIL,
            Verdict.FAIL,
        ]


class TestVerdictConfidence:
    """Test the Verdict.confidence property."""

    @pytest.mark.parametrize(
        "verdict, expected",
        [
            (Verdict.PASS, 1.0),
            (Verdict.LIKELY_PASS, 0.75),
            (Verdict.UNCERTAIN, 0.5),
            (Verdict.LIKELY_FAIL, 0.25),
            (Verdict.FAIL, 0.0),
        ],
    )
    def test_confidence_values(self, verdict, expected):
        assert verdict.confidence == expected

    def test_confidence_monotonically_decreases(self):
        """Confidence decreases monotonically from PASS to FAIL."""
        members = list(Verdict)
        for i in range(len(members) - 1):
            assert members[i].confidence > members[i + 1].confidence

    def test_confidence_bounded(self):
        for v in Verdict:
            assert 0.0 <= v.confidence <= 1.0


class TestVerdictIsDefinite:
    """Test the Verdict.is_definite property."""

    @pytest.mark.parametrize(
        "verdict, expected",
        [
            (Verdict.PASS, True),
            (Verdict.LIKELY_PASS, False),
            (Verdict.UNCERTAIN, False),
            (Verdict.LIKELY_FAIL, False),
            (Verdict.FAIL, True),
        ],
    )
    def test_is_definite_values(self, verdict, expected):
        assert verdict.is_definite is expected

    def test_only_endpoints_are_definite(self):
        definite = [v for v in Verdict if v.is_definite]
        assert definite == [Verdict.PASS, Verdict.FAIL]


class TestVerdictOrder:
    """Test the _VERDICT_ORDER internal mapping."""

    def test_pass_highest(self):
        assert _VERDICT_ORDER[Verdict.PASS] == 4

    def test_likely_pass(self):
        assert _VERDICT_ORDER[Verdict.LIKELY_PASS] == 3

    def test_uncertain_middle(self):
        assert _VERDICT_ORDER[Verdict.UNCERTAIN] == 2

    def test_likely_fail(self):
        assert _VERDICT_ORDER[Verdict.LIKELY_FAIL] == 1

    def test_fail_lowest(self):
        assert _VERDICT_ORDER[Verdict.FAIL] == 0

    def test_all_orders_distinct(self):
        orders = list(_VERDICT_ORDER.values())
        assert len(set(orders)) == 5

    def test_consistent_with_confidence(self):
        """Higher order should mean higher confidence."""
        for v1 in Verdict:
            for v2 in Verdict:
                if _VERDICT_ORDER[v1] > _VERDICT_ORDER[v2]:
                    assert v1.confidence > v2.confidence


# ═══════════════════════════════════════════════════════════════════════════
# 3. five_valued_and / five_valued_or
# ═══════════════════════════════════════════════════════════════════════════

class TestFiveValuedAnd:
    """Test the Kleene-style conjunction."""

    def test_pass_identity(self):
        """AND(PASS, x) = x for all x."""
        for v in Verdict:
            assert five_valued_and(Verdict.PASS, v) == v

    def test_pass_identity_right(self):
        """AND(x, PASS) = x for all x."""
        for v in Verdict:
            assert five_valued_and(v, Verdict.PASS) == v

    def test_fail_absorbs(self):
        """AND(FAIL, x) = FAIL for all x."""
        for v in Verdict:
            assert five_valued_and(Verdict.FAIL, v) == Verdict.FAIL

    def test_fail_absorbs_right(self):
        """AND(x, FAIL) = FAIL for all x."""
        for v in Verdict:
            assert five_valued_and(v, Verdict.FAIL) == Verdict.FAIL

    def test_commutative(self):
        for a in Verdict:
            for b in Verdict:
                assert five_valued_and(a, b) == five_valued_and(b, a)

    def test_idempotent(self):
        for v in Verdict:
            assert five_valued_and(v, v) == v

    def test_returns_minimum(self):
        """AND returns the verdict with the lower order (weakest link)."""
        for a in Verdict:
            for b in Verdict:
                expected = min(a, b, key=lambda v: _VERDICT_ORDER[v])
                assert five_valued_and(a, b) == expected

    def test_specific_pairs(self):
        assert five_valued_and(Verdict.PASS, Verdict.UNCERTAIN) == Verdict.UNCERTAIN
        assert five_valued_and(Verdict.LIKELY_PASS, Verdict.UNCERTAIN) == Verdict.UNCERTAIN
        assert five_valued_and(Verdict.LIKELY_PASS, Verdict.LIKELY_FAIL) == Verdict.LIKELY_FAIL
        assert five_valued_and(Verdict.UNCERTAIN, Verdict.UNCERTAIN) == Verdict.UNCERTAIN
        assert five_valued_and(Verdict.PASS, Verdict.PASS) == Verdict.PASS
        assert five_valued_and(Verdict.FAIL, Verdict.UNCERTAIN) == Verdict.FAIL


class TestFiveValuedOr:
    """Test the Kleene-style disjunction."""

    def test_fail_identity(self):
        """OR(FAIL, x) = x for all x."""
        for v in Verdict:
            assert five_valued_or(Verdict.FAIL, v) == v

    def test_fail_identity_right(self):
        """OR(x, FAIL) = x for all x."""
        for v in Verdict:
            assert five_valued_or(v, Verdict.FAIL) == v

    def test_pass_absorbs(self):
        """OR(PASS, x) = PASS for all x."""
        for v in Verdict:
            assert five_valued_or(Verdict.PASS, v) == Verdict.PASS

    def test_pass_absorbs_right(self):
        """OR(x, PASS) = PASS for all x."""
        for v in Verdict:
            assert five_valued_or(v, Verdict.PASS) == Verdict.PASS

    def test_commutative(self):
        for a in Verdict:
            for b in Verdict:
                assert five_valued_or(a, b) == five_valued_or(b, a)

    def test_idempotent(self):
        for v in Verdict:
            assert five_valued_or(v, v) == v

    def test_returns_maximum(self):
        """OR returns the verdict with the higher order (strongest link)."""
        for a in Verdict:
            for b in Verdict:
                expected = max(a, b, key=lambda v: _VERDICT_ORDER[v])
                assert five_valued_or(a, b) == expected

    def test_specific_pairs(self):
        assert five_valued_or(Verdict.PASS, Verdict.UNCERTAIN) == Verdict.PASS
        assert five_valued_or(Verdict.FAIL, Verdict.UNCERTAIN) == Verdict.UNCERTAIN
        assert five_valued_or(Verdict.LIKELY_PASS, Verdict.LIKELY_FAIL) == Verdict.LIKELY_PASS
        assert five_valued_or(Verdict.UNCERTAIN, Verdict.UNCERTAIN) == Verdict.UNCERTAIN
        assert five_valued_or(Verdict.FAIL, Verdict.FAIL) == Verdict.FAIL
        assert five_valued_or(Verdict.PASS, Verdict.FAIL) == Verdict.PASS


class TestBackwardCompatibleAliases:
    """Test that three_valued_and/or are aliases for five_valued_and/or."""

    def test_three_valued_and_is_five_valued_and(self):
        assert three_valued_and is five_valued_and

    def test_three_valued_or_is_five_valued_or(self):
        assert three_valued_or is five_valued_or

    def test_alias_produces_same_results(self):
        for a in Verdict:
            for b in Verdict:
                assert three_valued_and(a, b) == five_valued_and(a, b)
                assert three_valued_or(a, b) == five_valued_or(a, b)


# ═══════════════════════════════════════════════════════════════════════════
# 4. combined_verdict
# ═══════════════════════════════════════════════════════════════════════════

class TestCombinedVerdict:
    """Test the combined_verdict function."""

    def test_empty_list_returns_uncertain(self):
        assert combined_verdict([]) == Verdict.UNCERTAIN

    def test_single_pass(self):
        assert combined_verdict([Verdict.PASS]) == Verdict.PASS

    def test_single_fail(self):
        assert combined_verdict([Verdict.FAIL]) == Verdict.FAIL

    def test_single_uncertain(self):
        assert combined_verdict([Verdict.UNCERTAIN]) == Verdict.UNCERTAIN

    def test_all_pass_yields_pass(self):
        assert combined_verdict([Verdict.PASS, Verdict.PASS, Verdict.PASS]) == Verdict.PASS

    def test_any_fail_yields_fail(self):
        assert combined_verdict([Verdict.PASS, Verdict.FAIL, Verdict.PASS]) == Verdict.FAIL

    def test_pass_and_uncertain_yields_uncertain(self):
        assert combined_verdict([Verdict.PASS, Verdict.UNCERTAIN]) == Verdict.UNCERTAIN

    def test_fail_dominates(self):
        """FAIL always drags the combined result to FAIL."""
        for v in Verdict:
            if v == Verdict.FAIL:
                continue
            assert combined_verdict([v, Verdict.FAIL]) == Verdict.FAIL
            assert combined_verdict([Verdict.FAIL, v]) == Verdict.FAIL

    def test_weakest_link(self):
        """Combined verdict is the minimum element by order."""
        verdicts = [Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN]
        assert combined_verdict(verdicts) == Verdict.UNCERTAIN

    def test_likely_fail_weakest(self):
        verdicts = [Verdict.PASS, Verdict.LIKELY_PASS, Verdict.LIKELY_FAIL]
        assert combined_verdict(verdicts) == Verdict.LIKELY_FAIL

    def test_all_likelihoods(self):
        """Full range: PASS, LIKELY_PASS, UNCERTAIN, LIKELY_FAIL, FAIL."""
        verdicts = list(Verdict)
        assert combined_verdict(verdicts) == Verdict.FAIL

    def test_order_independence(self):
        """Permuting the input should not change the result."""
        import itertools
        verdicts = [Verdict.PASS, Verdict.UNCERTAIN, Verdict.LIKELY_PASS]
        base = combined_verdict(verdicts)
        for perm in itertools.permutations(verdicts):
            assert combined_verdict(list(perm)) == base

    def test_combined_pass_implies_all_pass(self):
        """If combined result is PASS, every individual verdict must be PASS."""
        verdicts = [Verdict.PASS, Verdict.PASS]
        result = combined_verdict(verdicts)
        assert result == Verdict.PASS
        # Counter-example
        verdicts2 = [Verdict.PASS, Verdict.LIKELY_PASS]
        assert combined_verdict(verdicts2) != Verdict.PASS

    def test_associative_folding(self):
        """combined_verdict([a,b,c]) == AND(AND(a,b),c)."""
        a, b, c = Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN
        assert combined_verdict([a, b, c]) == five_valued_and(
            five_valued_and(a, b), c
        )

    def test_large_list(self):
        """Handles larger lists correctly."""
        verdicts = [Verdict.PASS] * 100 + [Verdict.UNCERTAIN]
        assert combined_verdict(verdicts) == Verdict.UNCERTAIN


# ═══════════════════════════════════════════════════════════════════════════
# 5. PositionRange
# ═══════════════════════════════════════════════════════════════════════════

class TestPositionRange:
    """Test the PositionRange frozen dataclass."""

    def test_construction(self):
        pr = PositionRange(start=5, end=10)
        assert pr.start == 5
        assert pr.end == 10

    def test_len(self):
        pr = PositionRange(start=3, end=8)
        assert len(pr) == 5

    def test_len_zero_width(self):
        pr = PositionRange(start=5, end=5)
        assert len(pr) == 0

    def test_overlaps_true(self):
        pr1 = PositionRange(start=0, end=10)
        pr2 = PositionRange(start=5, end=15)
        assert pr1.overlaps(pr2) is True
        assert pr2.overlaps(pr1) is True

    def test_overlaps_false_adjacent(self):
        """Adjacent half-open intervals don't overlap: [0,10) and [10,20)."""
        pr1 = PositionRange(start=0, end=10)
        pr2 = PositionRange(start=10, end=20)
        assert pr1.overlaps(pr2) is False

    def test_overlaps_false_disjoint(self):
        pr1 = PositionRange(start=0, end=5)
        pr2 = PositionRange(start=10, end=15)
        assert pr1.overlaps(pr2) is False

    def test_overlaps_self(self):
        pr = PositionRange(start=0, end=10)
        assert pr.overlaps(pr) is True

    def test_overlaps_contained(self):
        pr1 = PositionRange(start=0, end=20)
        pr2 = PositionRange(start=5, end=10)
        assert pr1.overlaps(pr2) is True

    def test_contains_position(self):
        pr = PositionRange(start=5, end=10)
        assert pr.contains(5) is True   # inclusive start
        assert pr.contains(9) is True   # just before end
        assert pr.contains(10) is False  # exclusive end
        assert pr.contains(4) is False   # before start
        assert pr.contains(0) is False

    def test_contains_at_boundaries(self):
        pr = PositionRange(start=0, end=1)
        assert pr.contains(0) is True
        assert pr.contains(1) is False

    def test_frozen(self):
        """PositionRange is frozen (immutable)."""
        pr = PositionRange(start=0, end=10)
        with pytest.raises(AttributeError):
            pr.start = 5  # type: ignore[misc]

    def test_equality(self):
        pr1 = PositionRange(start=5, end=10)
        pr2 = PositionRange(start=5, end=10)
        assert pr1 == pr2

    def test_hash(self):
        """Frozen dataclass should be hashable."""
        pr = PositionRange(start=5, end=10)
        assert isinstance(hash(pr), int)
        # Same values should produce same hash
        pr2 = PositionRange(start=5, end=10)
        assert hash(pr) == hash(pr2)


# ═══════════════════════════════════════════════════════════════════════════
# 6. Token
# ═══════════════════════════════════════════════════════════════════════════

class TestToken:
    """Test the Token frozen dataclass."""

    def test_construction_minimal(self):
        t = Token(position=10, element_type="promoter", match_sequence="TATAAA")
        assert t.position == 10
        assert t.element_type == "promoter"
        assert t.match_sequence == "TATAAA"
        assert t.score == 0.0
        assert t.frame is None
        assert t.strand == "+"

    def test_construction_full(self):
        t = Token(
            position=20,
            element_type="start_codon",
            match_sequence="ATG",
            score=0.95,
            frame=0,
            strand="-",
        )
        assert t.position == 20
        assert t.element_type == "start_codon"
        assert t.match_sequence == "ATG"
        assert t.score == 0.95
        assert t.frame == 0
        assert t.strand == "-"

    def test_range_property(self):
        t = Token(position=10, element_type="exon", match_sequence="ATGCGT")
        r = t.range
        assert isinstance(r, PositionRange)
        assert r.start == 10
        assert r.end == 16
        assert len(r) == 6

    def test_range_single_base(self):
        t = Token(position=5, element_type="base", match_sequence="A")
        r = t.range
        assert r.start == 5
        assert r.end == 6

    def test_range_position_zero(self):
        t = Token(position=0, element_type="start", match_sequence="ATG")
        r = t.range
        assert r.start == 0
        assert r.end == 3

    def test_frozen(self):
        t = Token(position=10, element_type="test", match_sequence="AAA")
        with pytest.raises(AttributeError):
            t.position = 20  # type: ignore[misc]

    def test_equality(self):
        t1 = Token(position=10, element_type="exon", match_sequence="ATG")
        t2 = Token(position=10, element_type="exon", match_sequence="ATG")
        assert t1 == t2

    def test_inequality(self):
        t1 = Token(position=10, element_type="exon", match_sequence="ATG")
        t2 = Token(position=11, element_type="exon", match_sequence="ATG")
        assert t1 != t2

    def test_hash(self):
        t = Token(position=10, element_type="exon", match_sequence="ATG")
        assert isinstance(hash(t), int)

    def test_default_score(self):
        t = Token(position=0, element_type="test", match_sequence="GCG")
        assert t.score == 0.0

    def test_negative_score(self):
        t = Token(position=0, element_type="test", match_sequence="GCG", score=-1.5)
        assert t.score == -1.5

    def test_frame_none_default(self):
        t = Token(position=0, element_type="test", match_sequence="GCG")
        assert t.frame is None

    def test_frame_values(self):
        for frame in (0, 1, 2):
            t = Token(position=0, element_type="test", match_sequence="GCG", frame=frame)
            assert t.frame == frame

    def test_strand_default(self):
        t = Token(position=0, element_type="test", match_sequence="GCG")
        assert t.strand == "+"

    def test_strand_negative(self):
        t = Token(position=0, element_type="test", match_sequence="GCG", strand="-")
        assert t.strand == "-"


# ═══════════════════════════════════════════════════════════════════════════
# 7. SpliceIsoform
# ═══════════════════════════════════════════════════════════════════════════

class TestSpliceIsoform:
    """Test the SpliceIsoform frozen dataclass."""

    def test_construction_minimal(self):
        si = SpliceIsoform(
            sequence="ATGCGTAA",
            exon_boundaries=[(0, 3), (4, 8)],
            parse_path=["exon1", "exon2"],
        )
        assert si.sequence == "ATGCGTAA"
        assert si.exon_boundaries == [(0, 3), (4, 8)]
        assert si.parse_path == ["exon1", "exon2"]
        assert si.score == 0.0

    def test_construction_with_score(self):
        si = SpliceIsoform(
            sequence="ATGCGT",
            exon_boundaries=[(0, 6)],
            parse_path=["single"],
            score=0.87,
        )
        assert si.score == 0.87

    def test_repr(self):
        si = SpliceIsoform(
            sequence="ATGCGTAA",
            exon_boundaries=[(0, 3), (4, 8)],
            parse_path=["e1", "e2"],
            score=0.5,
        )
        r = repr(si)
        assert "SpliceIsoform" in r
        assert "len=8" in r
        assert "exons=2" in r
        assert "score=0.50" in r
        assert "path=" in r

    def test_frozen(self):
        si = SpliceIsoform(
            sequence="ATG",
            exon_boundaries=[(0, 3)],
            parse_path=["e1"],
        )
        with pytest.raises(AttributeError):
            si.sequence = "GCT"  # type: ignore[misc]

    def test_equality(self):
        si1 = SpliceIsoform(
            sequence="ATG",
            exon_boundaries=[(0, 3)],
            parse_path=["e1"],
        )
        si2 = SpliceIsoform(
            sequence="ATG",
            exon_boundaries=[(0, 3)],
            parse_path=["e1"],
        )
        assert si1 == si2

    def test_hash_unhashable_fields(self):
        """SpliceIsoform contains list fields which are unhashable;
        even though the dataclass is frozen, hash raises TypeError."""
        si = SpliceIsoform(
            sequence="ATG",
            exon_boundaries=[(0, 3)],
            parse_path=["e1"],
        )
        with pytest.raises(TypeError, match="unhashable"):
            hash(si)

    def test_empty_exon_boundaries(self):
        si = SpliceIsoform(
            sequence="",
            exon_boundaries=[],
            parse_path=[],
        )
        assert si.exon_boundaries == []
        assert len(si.sequence) == 0

    def test_repr_multiple_exons(self):
        si = SpliceIsoform(
            sequence="A" * 100,
            exon_boundaries=[(0, 30), (30, 60), (60, 100)],
            parse_path=["e1", "e2", "e3"],
            score=0.123,
        )
        r = repr(si)
        assert "len=100" in r
        assert "exons=3" in r


# ═══════════════════════════════════════════════════════════════════════════
# 8. TypeCheckResult
# ═══════════════════════════════════════════════════════════════════════════

class TestTypeCheckResult:
    """Test the TypeCheckResult mutable dataclass."""

    def test_construction_minimal(self):
        tcr = TypeCheckResult(predicate="NoStopCodons", verdict=Verdict.PASS)
        assert tcr.predicate == "NoStopCodons"
        assert tcr.verdict == Verdict.PASS
        assert tcr.derivation is None
        assert tcr.violation is None
        assert tcr.knowledge_gap is None

    def test_construction_full(self):
        tcr = TypeCheckResult(
            predicate="GCContent",
            verdict=Verdict.FAIL,
            derivation=[{"step": 1, "detail": "computed GC"}],
            violation="GC=0.25 below 0.30",
            knowledge_gap=None,
        )
        assert tcr.predicate == "GCContent"
        assert tcr.verdict == Verdict.FAIL
        assert len(tcr.derivation) == 1
        assert tcr.violation == "GC=0.25 below 0.30"
        assert tcr.knowledge_gap is None

    def test_passed_true_for_pass(self):
        tcr = TypeCheckResult(predicate="P", verdict=Verdict.PASS)
        assert tcr.passed is True

    def test_passed_true_for_likely_pass(self):
        tcr = TypeCheckResult(predicate="P", verdict=Verdict.LIKELY_PASS)
        assert tcr.passed is True

    def test_passed_false_for_uncertain(self):
        tcr = TypeCheckResult(predicate="P", verdict=Verdict.UNCERTAIN)
        assert tcr.passed is False

    def test_passed_false_for_likely_fail(self):
        tcr = TypeCheckResult(predicate="P", verdict=Verdict.LIKELY_FAIL)
        assert tcr.passed is False

    def test_passed_false_for_fail(self):
        tcr = TypeCheckResult(predicate="P", verdict=Verdict.FAIL)
        assert tcr.passed is False

    def test_repr(self):
        tcr = TypeCheckResult(predicate="NoStopCodons", verdict=Verdict.PASS)
        r = repr(tcr)
        assert "TypeCheckResult" in r
        assert "NoStopCodons" in r
        assert "PASS" in r

    def test_repr_fail(self):
        tcr = TypeCheckResult(predicate="GCContent", verdict=Verdict.FAIL)
        r = repr(tcr)
        assert "FAIL" in r
        assert "GCContent" in r

    def test_not_frozen(self):
        """TypeCheckResult is mutable (not frozen)."""
        tcr = TypeCheckResult(predicate="P", verdict=Verdict.UNCERTAIN)
        tcr.verdict = Verdict.PASS
        assert tcr.verdict == Verdict.PASS

    def test_with_knowledge_gap(self):
        tcr = TypeCheckResult(
            predicate="NoCrypticSplice",
            verdict=Verdict.UNCERTAIN,
            knowledge_gap="SLOT predicate in CONSERVATIVE mode",
        )
        assert tcr.knowledge_gap is not None
        assert "SLOT" in tcr.knowledge_gap

    def test_with_derivation_list(self):
        derivation = [
            {"step": "check", "positions": [10, 25]},
            {"step": "score", "value": 3.5},
        ]
        tcr = TypeCheckResult(
            predicate="Splice",
            verdict=Verdict.LIKELY_FAIL,
            derivation=derivation,
        )
        assert tcr.derivation is not None
        assert len(tcr.derivation) == 2


# ═══════════════════════════════════════════════════════════════════════════
# 9. Certificate
# ═══════════════════════════════════════════════════════════════════════════

class TestCertificate:
    """Test the Certificate dataclass and its serialization."""

    @pytest.fixture
    def sample_cert(self):
        return Certificate(
            version="1.0",
            design_id="test-001",
            sequence="ATGCGTAA",
            types=[{"predicate": "NoStopCodons", "verdict": "PASS"}],
            provenance={"organism": "Homo_sapiens"},
        )

    def test_construction(self, sample_cert):
        assert sample_cert.version == "1.0"
        assert sample_cert.design_id == "test-001"
        assert sample_cert.sequence == "ATGCGTAA"
        assert len(sample_cert.types) == 1
        assert sample_cert.provenance["organism"] == "Homo_sapiens"

    def test_to_dict(self, sample_cert):
        d = sample_cert.to_dict()
        assert d["version"] == "1.0"
        assert d["design_id"] == "test-001"
        assert d["sequence"] == "ATGCGTAA"
        assert d["types"] == [{"predicate": "NoStopCodons", "verdict": "PASS"}]
        assert d["provenance"] == {"organism": "Homo_sapiens"}

    def test_to_dict_round_trip(self, sample_cert):
        d = sample_cert.to_dict()
        cert2 = Certificate.from_dict(d)
        assert cert2.version == sample_cert.version
        assert cert2.design_id == sample_cert.design_id
        assert cert2.sequence == sample_cert.sequence
        assert cert2.types == sample_cert.types
        assert cert2.provenance == sample_cert.provenance

    def test_from_dict_valid(self):
        data = {
            "version": "2.0",
            "design_id": "cert-42",
            "sequence": "ATG",
            "types": [],
            "provenance": {"tool": "biocompiler"},
        }
        cert = Certificate.from_dict(data)
        assert cert.version == "2.0"
        assert cert.design_id == "cert-42"

    def test_from_dict_missing_version(self):
        data = {
            "design_id": "cert-42",
            "sequence": "ATG",
            "types": [],
            "provenance": {},
        }
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict(data)

    def test_from_dict_missing_design_id(self):
        data = {
            "version": "1.0",
            "sequence": "ATG",
            "types": [],
            "provenance": {},
        }
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict(data)

    def test_from_dict_missing_sequence(self):
        data = {
            "version": "1.0",
            "design_id": "test",
            "types": [],
            "provenance": {},
        }
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict(data)

    def test_from_dict_missing_types(self):
        data = {
            "version": "1.0",
            "design_id": "test",
            "sequence": "ATG",
            "provenance": {},
        }
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict(data)

    def test_from_dict_missing_provenance(self):
        data = {
            "version": "1.0",
            "design_id": "test",
            "sequence": "ATG",
            "types": [],
        }
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict(data)

    def test_from_dict_multiple_missing_keys(self):
        data = {"version": "1.0"}
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict(data)

    def test_from_dict_empty_raises(self):
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict({})

    def test_from_dict_extra_keys_ignored(self):
        """Extra keys in the dict are silently ignored (no error)."""
        data = {
            "version": "1.0",
            "design_id": "test",
            "sequence": "ATG",
            "types": [],
            "provenance": {},
            "extra": "should not cause error",
        }
        cert = Certificate.from_dict(data)
        assert cert.version == "1.0"

    def test_to_dict_json_compatible_types(self, sample_cert):
        """All values in to_dict should be JSON-compatible (str, list, dict)."""
        d = sample_cert.to_dict()
        import json
        # This should not raise
        json.dumps(d)

    def test_not_frozen(self):
        """Certificate is mutable (not frozen)."""
        cert = Certificate(
            version="1.0",
            design_id="test",
            sequence="ATG",
            types=[],
            provenance={},
        )
        cert.version = "2.0"
        assert cert.version == "2.0"

    def test_types_field_can_be_complex(self):
        """Types field can contain nested dicts."""
        types_data = [
            {"predicate": "NoStopCodons", "verdict": "PASS", "details": "ok"},
            {"predicate": "GCContent", "verdict": "FAIL", "violation": "low"},
        ]
        cert = Certificate(
            version="1.0",
            design_id="test",
            sequence="ATG",
            types=types_data,
            provenance={},
        )
        assert len(cert.types) == 2
        d = cert.to_dict()
        assert d["types"] == types_data

    def test_provenance_can_contain_nested_data(self):
        provenance = {
            "organism": "Homo_sapiens",
            "slot_mode": "conservative",
            "timestamp": "2024-01-01",
            "nested": {"key": "value"},
        }
        cert = Certificate(
            version="1.0",
            design_id="test",
            sequence="ATG",
            types=[],
            provenance=provenance,
        )
        d = cert.to_dict()
        assert d["provenance"]["nested"]["key"] == "value"


# ═══════════════════════════════════════════════════════════════════════════
# 10. __all__ exports
# ═══════════════════════════════════════════════════════════════════════════

class TestModuleExports:
    """Verify that __all__ contains all expected public names."""

    EXPECTED_EXPORTS = [
        "SLOTMode",
        "Verdict",
        "five_valued_and",
        "five_valued_or",
        "three_valued_and",
        "three_valued_or",
        "combined_verdict",
        "PositionRange",
        "Token",
        "SpliceIsoform",
        "TypeCheckResult",
        "Certificate",
    ]

    def test_all_contains_expected(self):
        for name in self.EXPECTED_EXPORTS:
            assert name in types_all, f"{name} missing from __all__"

    def test_all_no_extra_unexpected(self):
        """All entries in __all__ should be in our expected list."""
        for name in types_all:
            assert name in self.EXPECTED_EXPORTS, f"Unexpected export: {name}"

    def test_all_count(self):
        assert len(types_all) == len(self.EXPECTED_EXPORTS)

    def test_all_names_importable(self):
        """Every name in __all__ should be importable from the module."""
        import biocompiler.types as types_mod
        for name in types_all:
            assert hasattr(types_mod, name), f"{name} in __all__ but not accessible"


# ═══════════════════════════════════════════════════════════════════════════
# 11. Cross-type integration tests
# ═══════════════════════════════════════════════════════════════════════════

class TestCrossTypeIntegration:
    """Test interactions between multiple types."""

    def test_token_range_in_position_range(self):
        """Token.range returns a PositionRange that correctly describes the token."""
        t = Token(position=10, element_type="exon", match_sequence="ATGCGT")
        r = t.range
        assert r.contains(10)
        assert r.contains(15)
        assert not r.contains(16)
        assert not r.contains(9)

    def test_type_check_result_with_certificate_types(self):
        """TypeCheckResult can feed into Certificate.types list."""
        tcr = TypeCheckResult(predicate="NoStopCodons", verdict=Verdict.PASS)
        cert = Certificate(
            version="1.0",
            design_id="test",
            sequence="ATG",
            types=[{"predicate": tcr.predicate, "verdict": tcr.verdict.value}],
            provenance={},
        )
        assert cert.types[0]["verdict"] == "PASS"

    def test_combined_verdict_feeds_type_check_results(self):
        """combined_verdict can be used with TypeCheckResult verdicts."""
        results = [
            TypeCheckResult(predicate="P1", verdict=Verdict.PASS),
            TypeCheckResult(predicate="P2", verdict=Verdict.UNCERTAIN),
            TypeCheckResult(predicate="P3", verdict=Verdict.LIKELY_PASS),
        ]
        verdicts = [r.verdict for r in results]
        combined = combined_verdict(verdicts)
        assert combined == Verdict.UNCERTAIN

    def test_position_range_from_token_overlaps(self):
        """Two overlapping tokens produce overlapping PositionRanges."""
        t1 = Token(position=10, element_type="a", match_sequence="ATGCGT")
        t2 = Token(position=13, element_type="b", match_sequence="CGTAAA")
        assert t1.range.overlaps(t2.range)

    def test_position_range_from_token_no_overlap(self):
        """Two non-overlapping tokens produce non-overlapping PositionRanges."""
        t1 = Token(position=0, element_type="a", match_sequence="ATG")
        t2 = Token(position=5, element_type="b", match_sequence="CGT")
        assert not t1.range.overlaps(t2.range)
