"""
Unit tests for biocompiler.provenance — DecisionRecord, ProvenanceTracker,
and OptimizationProvenance.

Covers:
1. DecisionRecord construction, serialization (to_dict / from_dict)
2. ProvenanceTracker record_decision, get_decisions_for_position,
   get_full_audit_trail
3. ProvenanceTracker to_dict / to_json roundtrip
4. OptimizationProvenance construction and serialization
5. Seed parameter for reproducibility
"""

import json
from datetime import datetime, timezone

import pytest

from biocompiler.provenance import (
    DecisionRecord,
    OptimizationProvenance,
    ProvenanceTracker,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _ts() -> str:
    """Return a consistent ISO 8601 timestamp for tests."""
    return "2025-07-01T12:00:00+00:00"


def _make_record(
    position: int = 0,
    decision_type: str = "codon_selected",
    chosen_value: str = "ATG",
    alternatives: list[str] | None = None,
    rationale: str = "Start codon",
    context: dict | None = None,
) -> DecisionRecord:
    """Factory for DecisionRecord instances with sensible defaults."""
    return DecisionRecord(
        timestamp=_ts(),
        decision_type=decision_type,
        position=position,
        chosen_value=chosen_value,
        alternatives_considered=alternatives or [],
        rationale=rationale,
        constraint_context=context or {},
    )


def _make_tracker_with_decisions(
    seed: int = 42,
    n: int = 5,
) -> ProvenanceTracker:
    """Return a tracker pre-loaded with *n* codon_selected decisions."""
    tracker = ProvenanceTracker(seed=seed)
    for i in range(n):
        tracker.record_decision(DecisionRecord(
            timestamp=_ts(),
            decision_type="codon_selected",
            position=i * 3,
            chosen_value=f"COD{i}",
            alternatives_considered=[f"ALT{i}a", f"ALT{i}b"],
            rationale=f"Position {i * 3} best codon",
            constraint_context={"cai": 0.8 + i * 0.02, "gc": 0.50},
        ))
    return tracker


# ═══════════════════════════════════════════════════════════════════════════════
# 1. DecisionRecord — construction & serialization
# ═══════════════════════════════════════════════════════════════════════════════

class TestDecisionRecordConstruction:
    """DecisionRecord dataclass: construction and field access."""

    def test_basic_construction(self):
        rec = _make_record()
        assert rec.timestamp == _ts()
        assert rec.decision_type == "codon_selected"
        assert rec.position == 0
        assert rec.chosen_value == "ATG"
        assert rec.alternatives_considered == []
        assert rec.rationale == "Start codon"
        assert rec.constraint_context == {}

    def test_all_fields_populated(self):
        rec = DecisionRecord(
            timestamp="2025-07-01T12:00:00+00:00",
            decision_type="mutation_applied",
            position=42,
            chosen_value="V42I",
            alternatives_considered=["V42L", "V42F"],
            rationale="Remove GT dinucleotide while preserving hydrophobicity",
            constraint_context={"blosum62": 3, "gc": 0.52, "max_donor_score": 2.1},
        )
        assert rec.decision_type == "mutation_applied"
        assert rec.position == 42
        assert rec.chosen_value == "V42I"
        assert len(rec.alternatives_considered) == 2
        assert rec.rationale.startswith("Remove GT")
        assert rec.constraint_context["blosum62"] == 3

    def test_decision_type_varieties(self):
        """Verify common decision_type values are accepted."""
        for dtype in ("codon_selected", "mutation_applied", "constraint_relaxed"):
            rec = _make_record(decision_type=dtype)
            assert rec.decision_type == dtype

    def test_alternatives_preserve_order(self):
        """alternatives_considered must preserve solver preference order."""
        alts = ["GTG", "GTA", "GTT", "GTC"]
        rec = _make_record(alternatives=alts)
        assert rec.alternatives_considered == alts

    def test_empty_alternatives_is_valid(self):
        rec = _make_record(alternatives=[])
        assert rec.alternatives_considered == []

    def test_constraint_context_arbitrary_keys(self):
        """constraint_context should accept any JSON-compatible structure."""
        ctx = {
            "cai": 0.92,
            "gc": 0.54,
            "active_predicates": ["GCInRange", "NoCrypticSplice"],
            "nested": {"a": 1, "b": [2, 3]},
        }
        rec = _make_record(context=ctx)
        assert rec.constraint_context["active_predicates"][0] == "GCInRange"
        assert rec.constraint_context["nested"]["b"] == [2, 3]

    def test_position_negative(self):
        """Negative positions are accepted (dataclass has no validation)."""
        rec = _make_record(position=-1)
        assert rec.position == -1

    def test_position_large(self):
        rec = _make_record(position=999_999)
        assert rec.position == 999_999


class TestDecisionRecordSerialization:
    """DecisionRecord to_dict / from_dict roundtrip and edge cases."""

    def test_to_dict_keys(self):
        rec = _make_record()
        d = rec.to_dict()
        expected_keys = {
            "timestamp", "decision_type", "position", "chosen_value",
            "alternatives_considered", "rationale", "constraint_context",
            "cai_impact", "codon_before", "codon_after",
        }
        assert set(d.keys()) == expected_keys

    def test_to_dict_values_match_fields(self):
        rec = DecisionRecord(
            timestamp="2025-07-01T12:00:00+00:00",
            decision_type="codon_selected",
            position=9,
            chosen_value="GTC",
            alternatives_considered=["GTG", "GTA"],
            rationale="Highest CAI",
            constraint_context={"cai": 0.91},
        )
        d = rec.to_dict()
        assert d["timestamp"] == "2025-07-01T12:00:00+00:00"
        assert d["decision_type"] == "codon_selected"
        assert d["position"] == 9
        assert d["chosen_value"] == "GTC"
        assert d["alternatives_considered"] == ["GTG", "GTA"]
        assert d["rationale"] == "Highest CAI"
        assert d["constraint_context"] == {"cai": 0.91}

    def test_to_dict_returns_plain_types(self):
        """to_dict should return plain list/dict, not dataclass references."""
        rec = _make_record(alternatives=["A", "B"], context={"x": 1})
        d = rec.to_dict()
        assert isinstance(d["alternatives_considered"], list)
        assert isinstance(d["constraint_context"], dict)

    def test_to_dict_is_json_serializable(self):
        rec = _make_record(
            alternatives=["GTG"],
            context={"cai": 0.88, "tags": ["best"]},
        )
        json_str = json.dumps(rec.to_dict())
        assert isinstance(json_str, str)
        parsed = json.loads(json_str)
        assert parsed["chosen_value"] == "ATG"

    def test_from_dict_roundtrip_identical(self):
        rec = DecisionRecord(
            timestamp="2025-07-01T12:00:00+00:00",
            decision_type="constraint_relaxed",
            position=15,
            chosen_value="gc_lo:0.28",
            alternatives_considered=["gc_lo:0.25"],
            rationale="Tight GC infeasible at 0.30",
            constraint_context={"gc": 0.29},
        )
        d = rec.to_dict()
        restored = DecisionRecord.from_dict(d)
        assert restored.timestamp == rec.timestamp
        assert restored.decision_type == rec.decision_type
        assert restored.position == rec.position
        assert restored.chosen_value == rec.chosen_value
        assert restored.alternatives_considered == rec.alternatives_considered
        assert restored.rationale == rec.rationale
        assert restored.constraint_context == rec.constraint_context

    def test_from_dict_coerces_position_to_int(self):
        """Position may arrive as a JSON number (float); from_dict should coerce."""
        d = {
            "timestamp": _ts(),
            "decision_type": "codon_selected",
            "position": 12.0,  # float from JSON
            "chosen_value": "GTC",
            "alternatives_considered": [],
            "rationale": "",
            "constraint_context": {},
        }
        rec = DecisionRecord.from_dict(d)
        assert isinstance(rec.position, int)
        assert rec.position == 12

    def test_from_dict_missing_single_key(self):
        with pytest.raises(ValueError, match="missing keys"):
            DecisionRecord.from_dict({
                "timestamp": _ts(),
                "decision_type": "codon_selected",
                # position missing
                "chosen_value": "ATG",
                "alternatives_considered": [],
                "rationale": "",
                "constraint_context": {},
            })

    def test_from_dict_missing_multiple_keys(self):
        with pytest.raises(ValueError, match="missing keys"):
            DecisionRecord.from_dict({"timestamp": _ts()})

    def test_from_dict_empty_dict(self):
        with pytest.raises(ValueError, match="missing keys"):
            DecisionRecord.from_dict({})

    def test_from_dict_extra_keys_ignored(self):
        """Extra keys in the dict should not cause errors."""
        d = _make_record().to_dict()
        d["extra_key"] = "ignored"
        restored = DecisionRecord.from_dict(d)
        assert restored.chosen_value == "ATG"

    def test_from_dict_preserves_alternatives_order(self):
        alts = ["GTG", "GTA", "GTT", "GTC"]
        d = _make_record(alternatives=alts).to_dict()
        restored = DecisionRecord.from_dict(d)
        assert restored.alternatives_considered == alts

    def test_from_dict_deep_copies_alternatives_and_context(self):
        """from_dict should create new list/dict instances, not share refs."""
        d = _make_record(alternatives=["A"], context={"k": "v"}).to_dict()
        restored = DecisionRecord.from_dict(d)
        # Mutating restored should not affect the original dict
        restored.alternatives_considered.append("B")
        restored.constraint_context["k2"] = "v2"
        assert len(d["alternatives_considered"]) == 1
        assert "k2" not in d["constraint_context"]


# ═══════════════════════════════════════════════════════════════════════════════
# 2. ProvenanceTracker — record_decision, querying
# ═══════════════════════════════════════════════════════════════════════════════

class TestProvenanceTrackerRecordDecision:
    """ProvenanceTracker.record_decision: storage and type enforcement."""

    def test_record_single_decision(self):
        tracker = ProvenanceTracker(seed=0)
        rec = _make_record(position=3, chosen_value="GTC")
        tracker.record_decision(rec)
        assert len(tracker) == 1

    def test_record_multiple_decisions(self):
        tracker = ProvenanceTracker(seed=0)
        for i in range(10):
            tracker.record_decision(_make_record(position=i))
        assert len(tracker) == 10

    def test_record_non_decision_record_raises_type_error(self):
        tracker = ProvenanceTracker(seed=0)
        with pytest.raises(TypeError, match="Expected DecisionRecord"):
            tracker.record_decision("not a record")  # type: ignore[arg-type]

    def test_record_dict_raises_type_error(self):
        tracker = ProvenanceTracker(seed=0)
        with pytest.raises(TypeError, match="Expected DecisionRecord"):
            tracker.record_decision({"position": 0})  # type: ignore[arg-type]

    def test_record_none_raises_type_error(self):
        tracker = ProvenanceTracker(seed=0)
        with pytest.raises(TypeError, match="Expected DecisionRecord"):
            tracker.record_decision(None)  # type: ignore[arg-type]


class TestProvenanceTrackerGetPosition:
    """ProvenanceTracker.get_decisions_for_position: querying by position."""

    def test_returns_matching_decisions(self):
        tracker = ProvenanceTracker(seed=0)
        tracker.record_decision(_make_record(position=6, chosen_value="GTC"))
        tracker.record_decision(_make_record(position=9, chosen_value="AAA"))
        result = tracker.get_decisions_for_position(6)
        assert len(result) == 1
        assert result[0].chosen_value == "GTC"

    def test_empty_for_unknown_position(self):
        tracker = ProvenanceTracker(seed=0)
        assert tracker.get_decisions_for_position(999) == []

    def test_multiple_decisions_same_position(self):
        """Multiple decisions at the same position are all returned."""
        tracker = ProvenanceTracker(seed=0)
        for val in ["GTC", "GTG"]:
            tracker.record_decision(_make_record(position=6, chosen_value=val))
        result = tracker.get_decisions_for_position(6)
        assert len(result) == 2
        assert [d.chosen_value for d in result] == ["GTC", "GTG"]

    def test_chronological_order_preserved(self):
        """Decisions at same position should be in insertion order."""
        tracker = ProvenanceTracker(seed=0)
        labels = ["first", "second", "third"]
        for label in labels:
            tracker.record_decision(_make_record(
                position=0, chosen_value=label,
            ))
        result = tracker.get_decisions_for_position(0)
        assert [d.chosen_value for d in result] == labels

    def test_returns_copy_not_reference(self):
        """Returned list should be a copy; modifying it does not affect tracker."""
        tracker = ProvenanceTracker(seed=0)
        tracker.record_decision(_make_record(position=0))
        result = tracker.get_decisions_for_position(0)
        result.clear()
        assert len(tracker.get_decisions_for_position(0)) == 1

    def test_position_zero(self):
        """Position 0 is valid and distinct from 'no position'."""
        tracker = ProvenanceTracker(seed=0)
        tracker.record_decision(_make_record(position=0, chosen_value="ATG"))
        assert len(tracker.get_decisions_for_position(0)) == 1
        assert tracker.get_decisions_for_position(0)[0].chosen_value == "ATG"

    def test_negative_position(self):
        tracker = ProvenanceTracker(seed=0)
        tracker.record_decision(_make_record(position=-5))
        assert len(tracker.get_decisions_for_position(-5)) == 1


class TestProvenanceTrackerFullAuditTrail:
    """ProvenanceTracker.get_full_audit_trail: complete ordered list."""

    def test_empty_tracker(self):
        tracker = ProvenanceTracker(seed=0)
        assert tracker.get_full_audit_trail() == []

    def test_single_decision(self):
        tracker = ProvenanceTracker(seed=0)
        tracker.record_decision(_make_record(position=3))
        trail = tracker.get_full_audit_trail()
        assert len(trail) == 1
        assert trail[0].position == 3

    def test_order_is_chronological(self):
        """Audit trail preserves insertion order."""
        tracker = ProvenanceTracker(seed=0)
        for i in range(5):
            tracker.record_decision(_make_record(position=i * 10))
        trail = tracker.get_full_audit_trail()
        positions = [d.position for d in trail]
        assert positions == [0, 10, 20, 30, 40]

    def test_returns_copy(self):
        """Modifying the returned list does not affect the tracker."""
        tracker = ProvenanceTracker(seed=0)
        tracker.record_decision(_make_record())
        trail = tracker.get_full_audit_trail()
        trail.clear()
        assert len(tracker.get_full_audit_trail()) == 1

    def test_len_matches_audit_trail(self):
        tracker = _make_tracker_with_decisions(n=7)
        assert len(tracker) == 7
        assert len(tracker.get_full_audit_trail()) == 7


# ═══════════════════════════════════════════════════════════════════════════════
# 3. ProvenanceTracker — to_dict / to_json roundtrip
# ═══════════════════════════════════════════════════════════════════════════════

class TestProvenanceTrackerSerialization:
    """ProvenanceTracker to_dict / to_json / from_dict / from_json roundtrip."""

    def test_to_dict_keys(self):
        tracker = ProvenanceTracker(seed=42)
        d = tracker.to_dict()
        assert set(d.keys()) == {"seed", "decision_count", "decisions", "optimization_records"}

    def test_to_dict_seed(self):
        tracker = ProvenanceTracker(seed=99)
        assert tracker.to_dict()["seed"] == 99

    def test_to_dict_decision_count(self):
        tracker = _make_tracker_with_decisions(n=5)
        assert tracker.to_dict()["decision_count"] == 5

    def test_to_dict_empty_decisions(self):
        tracker = ProvenanceTracker(seed=0)
        d = tracker.to_dict()
        assert d["decision_count"] == 0
        assert d["decisions"] == []

    def test_to_dict_decisions_are_serialized(self):
        tracker = ProvenanceTracker(seed=0)
        tracker.record_decision(_make_record(position=6, chosen_value="GTC"))
        d = tracker.to_dict()
        assert len(d["decisions"]) == 1
        assert d["decisions"][0]["chosen_value"] == "GTC"

    def test_to_dict_is_json_serializable(self):
        tracker = _make_tracker_with_decisions(n=3)
        json_str = json.dumps(tracker.to_dict())
        assert isinstance(json_str, str)
        parsed = json.loads(json_str)
        assert parsed["seed"] == 42
        assert len(parsed["decisions"]) == 3

    def test_to_json_returns_string(self):
        tracker = ProvenanceTracker(seed=0)
        assert isinstance(tracker.to_json(), str)

    def test_to_json_is_valid_json(self):
        tracker = _make_tracker_with_decisions(n=2)
        parsed = json.loads(tracker.to_json())
        assert parsed["seed"] == 42

    def test_to_json_sorted_keys(self):
        """to_json uses sort_keys=True, so keys are alphabetically ordered."""
        tracker = ProvenanceTracker(seed=1)
        tracker.record_decision(_make_record())
        json_str = tracker.to_json()
        # Verify it parses and keys are sorted by checking raw string
        parsed = json.loads(json_str)
        top_keys = list(parsed.keys())
        assert top_keys == sorted(top_keys)

    def test_to_json_indented(self):
        """to_json uses indent=2 for human-readable output."""
        tracker = ProvenanceTracker(seed=0)
        json_str = tracker.to_json()
        assert "\n" in json_str
        assert "  " in json_str  # 2-space indentation

    def test_from_dict_roundtrip(self):
        tracker = _make_tracker_with_decisions(seed=42, n=3)
        d = tracker.to_dict()
        restored = ProvenanceTracker.from_dict(d)
        assert restored.seed == 42
        assert len(restored) == 3
        # Verify decisions are properly reconstructed
        original_trail = tracker.get_full_audit_trail()
        restored_trail = restored.get_full_audit_trail()
        for orig, rest in zip(original_trail, restored_trail):
            assert rest.chosen_value == orig.chosen_value
            assert rest.position == orig.position
            assert rest.alternatives_considered == orig.alternatives_considered

    def test_from_dict_position_index_rebuilt(self):
        """After from_dict, get_decisions_for_position should work."""
        tracker = ProvenanceTracker(seed=10)
        tracker.record_decision(_make_record(position=9, chosen_value="AAA"))
        tracker.record_decision(_make_record(position=9, chosen_value="AAG"))
        tracker.record_decision(_make_record(position=12, chosen_value="GCT"))

        restored = ProvenanceTracker.from_dict(tracker.to_dict())
        pos_9 = restored.get_decisions_for_position(9)
        assert len(pos_9) == 2
        assert [d.chosen_value for d in pos_9] == ["AAA", "AAG"]
        pos_12 = restored.get_decisions_for_position(12)
        assert len(pos_12) == 1
        assert pos_12[0].chosen_value == "GCT"

    def test_from_dict_missing_seed_raises(self):
        with pytest.raises(ValueError, match="missing key 'seed'"):
            ProvenanceTracker.from_dict({"decisions": []})

    def test_from_dict_empty_decisions_key(self):
        """Missing 'decisions' key should default to empty tracker."""
        tracker = ProvenanceTracker.from_dict({"seed": 7})
        assert tracker.seed == 7
        assert len(tracker) == 0

    def test_from_json_roundtrip(self):
        tracker = _make_tracker_with_decisions(seed=55, n=4)
        json_str = tracker.to_json()
        restored = ProvenanceTracker.from_json(json_str)
        assert restored.seed == 55
        assert len(restored) == 4

    def test_from_json_roundtrip_with_decisions(self):
        tracker = ProvenanceTracker(seed=1)
        tracker.record_decision(DecisionRecord(
            timestamp="2025-07-01T12:00:00+00:00",
            decision_type="mutation_applied",
            position=42,
            chosen_value="V42I",
            alternatives_considered=["V42L"],
            rationale="Remove GT",
            constraint_context={"blosum62": 3},
        ))
        json_str = tracker.to_json()
        restored = ProvenanceTracker.from_json(json_str)
        decisions = restored.get_decisions_for_position(42)
        assert len(decisions) == 1
        assert decisions[0].chosen_value == "V42I"
        assert decisions[0].constraint_context["blosum62"] == 3

    def test_from_dict_with_malformed_decision_raises(self):
        """If a nested decision dict is missing keys, from_dict should raise."""
        d = {
            "seed": 0,
            "decision_count": 1,
            "decisions": [{"timestamp": _ts()}],  # missing most keys
        }
        with pytest.raises(ValueError, match="missing keys"):
            ProvenanceTracker.from_dict(d)

    def test_repr(self):
        tracker = ProvenanceTracker(seed=42)
        assert "seed=42" in repr(tracker)
        assert "decisions=0" in repr(tracker)
        tracker.record_decision(_make_record())
        assert "decisions=1" in repr(tracker)


# ═══════════════════════════════════════════════════════════════════════════════
# 4. OptimizationProvenance — construction & serialization
# ═══════════════════════════════════════════════════════════════════════════════

class TestOptimizationProvenanceConstruction:
    """OptimizationProvenance: construction and field access."""

    def _make_prov(self, **overrides) -> OptimizationProvenance:
        defaults = dict(
            input_protein="MVSKGE",
            organism="Homo_sapiens",
            solver_backend="greedy",
            config_snapshot={"gc_lo": 0.30, "gc_hi": 0.70},
            decisions=[_make_record(position=0, chosen_value="ATG")],
            final_sequence="ATGGTTTCTAAAGGTGAA",
            solve_time_seconds=1.5,
            constraints_active=["GCInRange", "NoCrypticSplice"],
        )
        defaults.update(overrides)
        return OptimizationProvenance(**defaults)

    def test_basic_construction(self):
        prov = self._make_prov()
        assert prov.input_protein == "MVSKGE"
        assert prov.organism == "Homo_sapiens"
        assert prov.solver_backend == "greedy"
        assert prov.solve_time_seconds == 1.5
        assert len(prov.decisions) == 1
        assert len(prov.constraints_active) == 2

    def test_all_solver_backends(self):
        for backend in ("greedy", "z3", "ortools"):
            prov = self._make_prov(solver_backend=backend)
            assert prov.solver_backend == backend

    def test_empty_decisions(self):
        prov = self._make_prov(decisions=[])
        assert prov.decisions == []

    def test_multiple_decisions(self):
        decisions = [
            _make_record(position=i * 3, chosen_value=f"COD{i}")
            for i in range(5)
        ]
        prov = self._make_prov(decisions=decisions)
        assert len(prov.decisions) == 5

    def test_config_snapshot_nested(self):
        config = {
            "gc_lo": 0.30,
            "gc_hi": 0.70,
            "restriction_enzymes": ["EcoRI", "BamHI"],
            "splice_threshold": 3.0,
            "nested": {"cai_threshold": 0.5},
        }
        prov = self._make_prov(config_snapshot=config)
        assert prov.config_snapshot["restriction_enzymes"] == ["EcoRI", "BamHI"]
        assert prov.config_snapshot["nested"]["cai_threshold"] == 0.5

    def test_zero_solve_time(self):
        prov = self._make_prov(solve_time_seconds=0.0)
        assert prov.solve_time_seconds == 0.0

    def test_repr(self):
        prov = self._make_prov()
        r = repr(prov)
        assert "Homo_sapiens" in r
        assert "greedy" in r
        assert "decisions=" in r


class TestOptimizationProvenanceSerialization:
    """OptimizationProvenance to_dict / to_json / from_dict / from_json."""

    def _make_prov(self) -> OptimizationProvenance:
        return OptimizationProvenance(
            input_protein="MVSKGE",
            organism="Homo_sapiens",
            solver_backend="z3",
            config_snapshot={"gc_lo": 0.30, "gc_hi": 0.70},
            decisions=[
                DecisionRecord(
                    timestamp="2025-07-01T12:00:00+00:00",
                    decision_type="codon_selected",
                    position=0,
                    chosen_value="ATG",
                    alternatives_considered=[],
                    rationale="Start codon",
                    constraint_context={},
                ),
                DecisionRecord(
                    timestamp="2025-07-01T12:00:01+00:00",
                    decision_type="codon_selected",
                    position=3,
                    chosen_value="GTT",
                    alternatives_considered=["GTC", "GTA", "GTG"],
                    rationale="Best CAI for Valine",
                    constraint_context={"cai": 0.88, "gc": 0.50},
                ),
            ],
            final_sequence="ATGGTTTCTAAAGGTGAA",
            solve_time_seconds=2.37,
            constraints_active=["GCInRange", "NoCrypticSplice", "NoRestrictionSite"],
        )

    def test_to_dict_keys(self):
        prov = self._make_prov()
        d = prov.to_dict()
        expected_keys = {
            "input_protein", "organism", "solver_backend",
            "config_snapshot", "decisions", "final_sequence",
            "solve_time_seconds", "constraints_active",
        }
        assert set(d.keys()) == expected_keys

    def test_to_dict_values(self):
        prov = self._make_prov()
        d = prov.to_dict()
        assert d["input_protein"] == "MVSKGE"
        assert d["organism"] == "Homo_sapiens"
        assert d["solver_backend"] == "z3"
        assert d["solve_time_seconds"] == 2.37
        assert len(d["decisions"]) == 2
        assert d["constraints_active"] == [
            "GCInRange", "NoCrypticSplice", "NoRestrictionSite"
        ]

    def test_to_dict_config_snapshot(self):
        prov = self._make_prov()
        d = prov.to_dict()
        assert d["config_snapshot"]["gc_lo"] == 0.30
        assert d["config_snapshot"]["gc_hi"] == 0.70

    def test_to_dict_decisions_serialized(self):
        prov = self._make_prov()
        d = prov.to_dict()
        assert d["decisions"][0]["chosen_value"] == "ATG"
        assert d["decisions"][1]["alternatives_considered"] == [
            "GTC", "GTA", "GTG"
        ]

    def test_to_dict_is_json_serializable(self):
        prov = self._make_prov()
        json_str = json.dumps(prov.to_dict())
        parsed = json.loads(json_str)
        assert parsed["organism"] == "Homo_sapiens"

    def test_from_dict_roundtrip(self):
        prov = self._make_prov()
        d = prov.to_dict()
        restored = OptimizationProvenance.from_dict(d)
        assert restored.input_protein == prov.input_protein
        assert restored.organism == prov.organism
        assert restored.solver_backend == prov.solver_backend
        assert restored.final_sequence == prov.final_sequence
        assert restored.solve_time_seconds == prov.solve_time_seconds
        assert restored.constraints_active == prov.constraints_active
        assert len(restored.decisions) == len(prov.decisions)
        assert restored.decisions[0].chosen_value == "ATG"
        assert restored.decisions[1].chosen_value == "GTT"

    def test_from_dict_coerces_solve_time(self):
        """solve_time_seconds may arrive as int from JSON; from_dict coerces."""
        d = self._make_prov().to_dict()
        d["solve_time_seconds"] = 3  # int instead of float
        restored = OptimizationProvenance.from_dict(d)
        assert isinstance(restored.solve_time_seconds, float)
        assert restored.solve_time_seconds == 3.0

    def test_from_dict_missing_keys_raises(self):
        with pytest.raises(ValueError, match="missing keys"):
            OptimizationProvenance.from_dict({"organism": "x"})

    def test_from_dict_missing_all_keys_raises(self):
        with pytest.raises(ValueError, match="missing keys"):
            OptimizationProvenance.from_dict({})

    def test_from_dict_partial_keys_raises(self):
        d = {
            "input_protein": "M",
            "organism": "E_coli",
            # missing: solver_backend, config_snapshot, decisions, etc.
        }
        with pytest.raises(ValueError, match="missing keys"):
            OptimizationProvenance.from_dict(d)

    def test_to_json_roundtrip(self):
        prov = self._make_prov()
        json_str = prov.to_json()
        restored = OptimizationProvenance.from_json(json_str)
        assert restored.input_protein == prov.input_protein
        assert restored.organism == prov.organism
        assert restored.solver_backend == prov.solver_backend

    def test_to_json_validates_parseable(self):
        prov = self._make_prov()
        json_str = prov.to_json()
        parsed = json.loads(json_str)
        assert "decisions" in parsed
        assert len(parsed["decisions"]) == 2

    def test_empty_decisions_roundtrip(self):
        prov = OptimizationProvenance(
            input_protein="M",
            organism="E_coli",
            solver_backend="greedy",
            config_snapshot={},
            decisions=[],
            final_sequence="ATG",
            solve_time_seconds=0.001,
            constraints_active=[],
        )
        d = prov.to_dict()
        restored = OptimizationProvenance.from_dict(d)
        assert restored.decisions == []
        assert restored.constraints_active == []

    def test_config_snapshot_roundtrip(self):
        config = {
            "gc_lo": 0.30,
            "gc_hi": 0.70,
            "restriction_enzymes": ["EcoRI", "BamHI"],
        }
        prov = OptimizationProvenance(
            input_protein="MK",
            organism="Homo_sapiens",
            solver_backend="ortools",
            config_snapshot=config,
            decisions=[],
            final_sequence="ATGAAA",
            solve_time_seconds=5.0,
            constraints_active=["GCInRange"],
        )
        restored = OptimizationProvenance.from_dict(prov.to_dict())
        assert restored.config_snapshot["restriction_enzymes"] == ["EcoRI", "BamHI"]


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Seed parameter for reproducibility
# ═══════════════════════════════════════════════════════════════════════════════

class TestSeedReproducibility:
    """Seed parameter ensures reproducibility: same seed + same input = same trail."""

    def test_default_seed_is_zero(self):
        tracker = ProvenanceTracker()
        assert tracker.seed == 0

    def test_custom_seed_stored(self):
        tracker = ProvenanceTracker(seed=12345)
        assert tracker.seed == 12345

    def test_negative_seed(self):
        tracker = ProvenanceTracker(seed=-1)
        assert tracker.seed == -1

    def test_large_seed(self):
        tracker = ProvenanceTracker(seed=2**31 - 1)
        assert tracker.seed == 2**31 - 1

    def test_zero_seed(self):
        tracker = ProvenanceTracker(seed=0)
        assert tracker.seed == 0

    def test_seed_preserved_in_to_dict(self):
        tracker = ProvenanceTracker(seed=42)
        assert tracker.to_dict()["seed"] == 42

    def test_seed_preserved_in_to_json(self):
        tracker = ProvenanceTracker(seed=42)
        parsed = json.loads(tracker.to_json())
        assert parsed["seed"] == 42

    def test_seed_restored_from_dict(self):
        tracker = ProvenanceTracker(seed=999)
        tracker.record_decision(_make_record())
        restored = ProvenanceTracker.from_dict(tracker.to_dict())
        assert restored.seed == 999

    def test_seed_restored_from_json(self):
        tracker = ProvenanceTracker(seed=777)
        tracker.record_decision(_make_record())
        restored = ProvenanceTracker.from_json(tracker.to_json())
        assert restored.seed == 777

    def test_same_seed_same_audit_trail_after_serialization(self):
        """Two trackers with the same seed and same decisions should serialize
        identically — demonstrating reproducibility."""
        decisions = [
            _make_record(position=i, chosen_value=f"COD{i}")
            for i in range(5)
        ]
        t1 = ProvenanceTracker(seed=42)
        t2 = ProvenanceTracker(seed=42)
        for d in decisions:
            t1.record_decision(d)
            t2.record_decision(d)
        assert t1.to_json() == t2.to_json()

    def test_different_seed_different_serialization(self):
        """Different seeds should produce different serialized output."""
        t1 = ProvenanceTracker(seed=1)
        t2 = ProvenanceTracker(seed=2)
        t1.record_decision(_make_record())
        t2.record_decision(_make_record())
        # JSON differs because seed field differs
        assert t1.to_json() != t2.to_json()

    def test_seed_in_optimization_provenance(self):
        """OptimizationProvenance should carry seed from tracker in config."""
        tracker = ProvenanceTracker(seed=42)
        prov = OptimizationProvenance(
            input_protein="MK",
            organism="Homo_sapiens",
            solver_backend="greedy",
            config_snapshot={"seed": tracker.seed},
            decisions=tracker.get_full_audit_trail(),
            final_sequence="ATGAAA",
            solve_time_seconds=0.5,
            constraints_active=["GCInRange"],
        )
        d = prov.to_dict()
        assert d["config_snapshot"]["seed"] == 42

    def test_full_reproducibility_workflow(self):
        """End-to-end: build tracker → serialize → restore → verify identical."""
        tracker = ProvenanceTracker(seed=42)
        # Simulate a full optimization decision trail
        protein = "MVSKGE"
        for i, aa in enumerate(protein):
            tracker.record_decision(DecisionRecord(
                timestamp=datetime.now(timezone.utc).isoformat(),
                decision_type="codon_selected",
                position=i * 3,
                chosen_value=f"COD_{aa}",
                alternatives_considered=[f"ALT_{aa}_1", f"ALT_{aa}_2"],
                rationale=f"Best codon for {aa} at position {i}",
                constraint_context={"cai": 0.9, "gc": 0.52},
            ))

        # Serialize to JSON and back
        json_str = tracker.to_json()
        restored = ProvenanceTracker.from_json(json_str)

        # Verify complete fidelity
        assert restored.seed == tracker.seed
        assert len(restored) == len(tracker)
        original_trail = tracker.get_full_audit_trail()
        restored_trail = restored.get_full_audit_trail()
        for orig, rest in zip(original_trail, restored_trail):
            assert rest.position == orig.position
            assert rest.chosen_value == orig.chosen_value
            assert rest.alternatives_considered == orig.alternatives_considered
            assert rest.rationale == orig.rationale
            assert rest.constraint_context == orig.constraint_context

        # Verify position-based queries work identically
        for i in range(len(protein)):
            pos = i * 3
            orig_decisions = tracker.get_decisions_for_position(pos)
            rest_decisions = restored.get_decisions_for_position(pos)
            assert len(orig_decisions) == len(rest_decisions)
            assert orig_decisions[0].chosen_value == rest_decisions[0].chosen_value

    def test_reproducibility_across_independent_runs(self):
        """Two independent trackers with same seed + same decisions = same output."""
        seed = 12345
        records = [
            _make_record(position=i, chosen_value=f"V{i}")
            for i in range(8)
        ]

        t1 = ProvenanceTracker(seed=seed)
        t2 = ProvenanceTracker(seed=seed)
        for rec in records:
            t1.record_decision(rec)
            t2.record_decision(rec)

        # Both to_dict and to_json must be identical
        assert t1.to_dict() == t2.to_dict()
        assert t1.to_json() == t2.to_json()


# ═══════════════════════════════════════════════════════════════════════════════
# Cross-cutting: integration between tracker and OptimizationProvenance
# ═══════════════════════════════════════════════════════════════════════════════

class TestTrackerProvenanceIntegration:
    """Verify that tracker and OptimizationProvenance work together correctly."""

    def test_tracker_decisions_feed_into_provenance(self):
        tracker = ProvenanceTracker(seed=42)
        for i in range(3):
            tracker.record_decision(_make_record(position=i * 3))

        prov = OptimizationProvenance(
            input_protein="MVL",
            organism="Homo_sapiens",
            solver_backend="greedy",
            config_snapshot={"seed": tracker.seed},
            decisions=tracker.get_full_audit_trail(),
            final_sequence="ATGGTGCTG",
            solve_time_seconds=0.3,
            constraints_active=["GCInRange"],
        )
        assert len(prov.decisions) == 3
        d = prov.to_dict()
        assert len(d["decisions"]) == 3

    def test_roundtrip_through_json(self):
        """Full workflow: tracker → provenance → JSON → restore → verify."""
        tracker = ProvenanceTracker(seed=42)
        tracker.record_decision(DecisionRecord(
            timestamp="2025-07-01T12:00:00+00:00",
            decision_type="codon_selected",
            position=0,
            chosen_value="ATG",
            alternatives_considered=[],
            rationale="Start codon, fixed",
            constraint_context={},
        ))

        prov = OptimizationProvenance(
            input_protein="M",
            organism="Homo_sapiens",
            solver_backend="z3",
            config_snapshot={"gc_lo": 0.3, "gc_hi": 0.7, "seed": 42},
            decisions=tracker.get_full_audit_trail(),
            final_sequence="ATG",
            solve_time_seconds=0.01,
            constraints_active=["NoStopCodons"],
        )

        json_str = prov.to_json()
        restored = OptimizationProvenance.from_json(json_str)
        assert restored.config_snapshot["seed"] == 42
        assert len(restored.decisions) == 1
        assert restored.decisions[0].chosen_value == "ATG"
