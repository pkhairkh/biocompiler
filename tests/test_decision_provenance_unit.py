"""
Unit tests for biocompiler.decision_provenance — CodonDecision, ConstraintDecision,
OptimizationDecisionTrail, DecisionProvenanceCollector.

These are pure unit tests (no optimization run required), covering:
- CodonDecision: construction, serialization (to_dict/from_dict), validation
- ConstraintDecision: construction, validation of action_taken, serialization
- OptimizationDecisionTrail: construction, to_dict/from_dict, to_json/from_json, repr
- DecisionProvenanceCollector: full lifecycle, error handling, summary generation
"""

from __future__ import annotations

import json
import os
import tempfile

import pytest

from biocompiler.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    OptimizationDecisionTrail,
    DecisionProvenanceCollector,
)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. CodonDecision
# ═══════════════════════════════════════════════════════════════════════════════

class TestCodonDecision:

    def test_construction(self):
        cd = CodonDecision(
            position=0,
            amino_acid="M",
            original_codon=None,
            chosen_codon="ATG",
            alternatives_considered=[],
            constraint_reason="maximize_cai",
            confidence=1.0,
        )
        assert cd.position == 0
        assert cd.amino_acid == "M"
        assert cd.chosen_codon == "ATG"
        assert cd.confidence == 1.0
        assert cd.cai_impact == 0.0  # default

    def test_cai_impact_default(self):
        cd = CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
        )
        assert cd.cai_impact == 0.0

    def test_cai_impact_explicit(self):
        cd = CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="gc_content", confidence=0.7,
            cai_impact=-0.005,
        )
        assert cd.cai_impact == -0.005

    def test_to_dict(self):
        cd = CodonDecision(
            position=2, amino_acid="L", original_codon="CTT",
            chosen_codon="CTG", alternatives_considered=[
                {"codon": "CTC", "cai_contribution": 0.8},
            ],
            constraint_reason="maximize_cai", confidence=0.95,
            cai_impact=-0.002,
        )
        d = cd.to_dict()
        assert d["position"] == 2
        assert d["amino_acid"] == "L"
        assert d["original_codon"] == "CTT"
        assert d["chosen_codon"] == "CTG"
        assert len(d["alternatives_considered"]) == 1
        assert d["constraint_reason"] == "maximize_cai"
        assert d["confidence"] == 0.95
        assert d["cai_impact"] == -0.002

    def test_from_dict(self):
        data = {
            "position": 3,
            "amino_acid": "V",
            "original_codon": "GTC",
            "chosen_codon": "GTG",
            "alternatives_considered": [
                {"codon": "GTT", "cai_contribution": 0.5},
            ],
            "constraint_reason": "avoid_restriction_site:EcoRI",
            "confidence": 0.8,
            "cai_impact": -0.01,
        }
        cd = CodonDecision.from_dict(data)
        assert cd.position == 3
        assert cd.amino_acid == "V"
        assert cd.chosen_codon == "GTG"
        assert cd.cai_impact == -0.01

    def test_from_dict_missing_keys_raises(self):
        data = {"position": 0}  # missing most keys
        with pytest.raises(ValueError, match="missing keys"):
            CodonDecision.from_dict(data)

    def test_from_dict_cai_impact_default(self):
        data = {
            "position": 0, "amino_acid": "M", "original_codon": None,
            "chosen_codon": "ATG", "alternatives_considered": [],
            "constraint_reason": "maximize_cai", "confidence": 1.0,
        }
        cd = CodonDecision.from_dict(data)
        assert cd.cai_impact == 0.0

    def test_roundtrip_to_dict_from_dict(self):
        cd = CodonDecision(
            position=5, amino_acid="G", original_codon="GGT",
            chosen_codon="GGC", alternatives_considered=[
                {"codon": "GGA", "cai_contribution": 0.3},
                {"codon": "GGG", "cai_contribution": 0.6},
            ],
            constraint_reason="gc_content", confidence=0.85,
            cai_impact=-0.003,
        )
        d = cd.to_dict()
        cd2 = CodonDecision.from_dict(d)
        assert cd2.position == cd.position
        assert cd2.amino_acid == cd.amino_acid
        assert cd2.chosen_codon == cd.chosen_codon
        assert cd2.cai_impact == cd.cai_impact
        assert cd2.confidence == cd.confidence


# ═══════════════════════════════════════════════════════════════════════════════
# 2. ConstraintDecision
# ═══════════════════════════════════════════════════════════════════════════════

class TestConstraintDecision:

    def test_construction_satisfied(self):
        cd = ConstraintDecision(
            constraint_name="GCInRange",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[3, 7],
            tradeoff_description="GC kept in range",
            impact_on_cai=-0.005,
        )
        assert cd.constraint_name == "GCInRange"
        assert cd.action_taken == "satisfied"
        assert cd.impact_on_cai == -0.005

    def test_invalid_action_raises(self):
        with pytest.raises(ValueError, match="Invalid action_taken"):
            ConstraintDecision(
                constraint_name="GCInRange",
                constraint_type="hard",
                action_taken="invalid_action",
                positions_affected=[],
                tradeoff_description="",
                impact_on_cai=0.0,
            )

    def test_all_valid_actions(self):
        for action in ("satisfied", "relaxed", "conflicted", "overridden"):
            cd = ConstraintDecision(
                constraint_name="Test",
                constraint_type="soft",
                action_taken=action,
                positions_affected=[],
                tradeoff_description="",
                impact_on_cai=0.0,
            )
            assert cd.action_taken == action

    def test_to_dict(self):
        cd = ConstraintDecision(
            constraint_name="NoCrypticSplice",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[5, 10],
            tradeoff_description="Eliminated cryptic donor at position 5",
            impact_on_cai=-0.02,
        )
        d = cd.to_dict()
        assert d["constraint_name"] == "NoCrypticSplice"
        assert d["action_taken"] == "satisfied"
        assert d["positions_affected"] == [5, 10]
        assert d["impact_on_cai"] == -0.02

    def test_from_dict(self):
        data = {
            "constraint_name": "NoRestrictionSite",
            "constraint_type": "hard",
            "action_taken": "relaxed",
            "positions_affected": [1, 2],
            "tradeoff_description": "Could not remove site",
            "impact_on_cai": -0.05,
        }
        cd = ConstraintDecision.from_dict(data)
        assert cd.constraint_name == "NoRestrictionSite"
        assert cd.action_taken == "relaxed"

    def test_from_dict_missing_keys_raises(self):
        with pytest.raises(ValueError, match="missing keys"):
            ConstraintDecision.from_dict({"constraint_name": "X"})

    def test_roundtrip_to_dict_from_dict(self):
        cd = ConstraintDecision(
            constraint_name="GCInRange",
            constraint_type="hard",
            action_taken="conflicted",
            positions_affected=[0, 1, 2],
            tradeoff_description="Tradeoff between GC and CAI",
            impact_on_cai=-0.1,
        )
        d = cd.to_dict()
        cd2 = ConstraintDecision.from_dict(d)
        assert cd2.constraint_name == cd.constraint_name
        assert cd2.action_taken == cd.action_taken
        assert cd2.impact_on_cai == cd.impact_on_cai


# ═══════════════════════════════════════════════════════════════════════════════
# 3. OptimizationDecisionTrail
# ═══════════════════════════════════════════════════════════════════════════════

def _make_trail():
    """Create a minimal OptimizationDecisionTrail for testing."""
    return OptimizationDecisionTrail(
        gene_name="TestGene",
        input_protein="MV",
        output_dna="ATGGTG",
        organism="Homo_sapiens",
        solver_backend="greedy",
        seed=42,
        total_cai=0.91,
        total_gc=0.54,
        codon_decisions=[
            CodonDecision(
                position=0, amino_acid="M", original_codon=None,
                chosen_codon="ATG", alternatives_considered=[],
                constraint_reason="maximize_cai", confidence=1.0,
            ),
            CodonDecision(
                position=1, amino_acid="V", original_codon=None,
                chosen_codon="GTG", alternatives_considered=[
                    {"codon": "GTC", "cai_contribution": 0.7},
                ],
                constraint_reason="maximize_cai", confidence=0.9,
                cai_impact=-0.002,
            ),
        ],
        constraint_decisions=[
            ConstraintDecision(
                constraint_name="GCInRange",
                constraint_type="hard",
                action_taken="satisfied",
                positions_affected=[0, 1],
                tradeoff_description="GC kept in range",
                impact_on_cai=-0.005,
            ),
        ],
        iteration_log=[
            {"step": 1, "action": "initial_assign", "score": 0.85},
        ],
        timestamp="2025-01-15T12:00:00+00:00",
        version="10.0.0",
    )


class TestOptimizationDecisionTrail:

    def test_construction(self):
        trail = _make_trail()
        assert trail.gene_name == "TestGene"
        assert trail.total_cai == 0.91
        assert len(trail.codon_decisions) == 2
        assert len(trail.constraint_decisions) == 1

    def test_to_dict(self):
        trail = _make_trail()
        d = trail.to_dict()
        assert d["gene_name"] == "TestGene"
        assert d["organism"] == "Homo_sapiens"
        assert len(d["codon_decisions"]) == 2
        assert len(d["constraint_decisions"]) == 1

    def test_to_json(self):
        trail = _make_trail()
        j = trail.to_json()
        parsed = json.loads(j)
        assert parsed["gene_name"] == "TestGene"

    def test_from_dict(self):
        trail = _make_trail()
        d = trail.to_dict()
        trail2 = OptimizationDecisionTrail.from_dict(d)
        assert trail2.gene_name == trail.gene_name
        assert trail2.total_cai == trail.total_cai
        assert len(trail2.codon_decisions) == len(trail.codon_decisions)

    def test_from_dict_missing_keys_raises(self):
        with pytest.raises(ValueError, match="missing keys"):
            OptimizationDecisionTrail.from_dict({"gene_name": "X"})

    def test_from_json(self):
        trail = _make_trail()
        j = trail.to_json()
        trail2 = OptimizationDecisionTrail.from_json(j)
        assert trail2.gene_name == trail.gene_name
        assert len(trail2.codon_decisions) == 2

    def test_repr(self):
        trail = _make_trail()
        r = repr(trail)
        assert "OptimizationDecisionTrail" in r
        assert "TestGene" in r
        assert "greedy" in r

    def test_roundtrip_preserves_data(self):
        trail = _make_trail()
        d = trail.to_dict()
        trail2 = OptimizationDecisionTrail.from_dict(d)
        d2 = trail2.to_dict()
        # Compare key fields
        assert d["gene_name"] == d2["gene_name"]
        assert d["total_cai"] == d2["total_cai"]
        assert d["total_gc"] == d2["total_gc"]
        assert len(d["codon_decisions"]) == len(d2["codon_decisions"])


# ═══════════════════════════════════════════════════════════════════════════════
# 4. DecisionProvenanceCollector
# ═══════════════════════════════════════════════════════════════════════════════

class TestDecisionProvenanceCollector:

    def test_initial_state(self):
        c = DecisionProvenanceCollector()
        assert repr(c) == "DecisionProvenanceCollector(state=idle, codon_decisions=0, constraint_decisions=0, iterations=0)"

    def test_full_lifecycle(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(
            protein="MV",
            organism="Homo_sapiens",
            constraints=["GCInRange"],
            gene_name="TestGene",
            solver_backend="greedy",
            seed=42,
        )
        c.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
        ))
        c.record_codon_decision(CodonDecision(
            position=1, amino_acid="V", original_codon=None,
            chosen_codon="GTG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=0.9,
        ))
        c.record_constraint_decision(ConstraintDecision(
            constraint_name="GCInRange",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[0, 1],
            tradeoff_description="GC in range",
            impact_on_cai=-0.005,
        ))
        c.record_iteration({"step": 1, "action": "assign", "score": 0.9})
        trail = c.finalize(output_dna="ATGGTG", cai=0.91, gc=0.54)
        assert isinstance(trail, OptimizationDecisionTrail)
        assert len(trail.codon_decisions) == 2
        assert len(trail.constraint_decisions) == 1
        assert len(trail.iteration_log) == 1

    def test_record_before_start_raises(self):
        c = DecisionProvenanceCollector()
        with pytest.raises(RuntimeError, match="start_optimization"):
            c.record_codon_decision(CodonDecision(
                position=0, amino_acid="M", original_codon=None,
                chosen_codon="ATG", alternatives_considered=[],
                constraint_reason="maximize_cai", confidence=1.0,
            ))

    def test_record_constraint_before_start_raises(self):
        c = DecisionProvenanceCollector()
        with pytest.raises(RuntimeError, match="start_optimization"):
            c.record_constraint_decision(ConstraintDecision(
                constraint_name="X", constraint_type="hard",
                action_taken="satisfied", positions_affected=[],
                tradeoff_description="", impact_on_cai=0.0,
            ))

    def test_record_iteration_before_start_raises(self):
        c = DecisionProvenanceCollector()
        with pytest.raises(RuntimeError, match="start_optimization"):
            c.record_iteration({"step": 1})

    def test_finalize_before_start_raises(self):
        c = DecisionProvenanceCollector()
        with pytest.raises(RuntimeError, match="start_optimization"):
            c.finalize(output_dna="ATG", cai=1.0, gc=0.5)

    def test_double_start_raises(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        with pytest.raises(RuntimeError, match="already in progress"):
            c.start_optimization(protein="LL", organism="Ecoli")

    def test_finalize_then_restart(self):
        """After finalize, a new optimization can be started."""
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        c.finalize(output_dna="ATGGTG", cai=0.9, gc=0.5)
        # Should be able to start a new one
        c.start_optimization(protein="LL", organism="Ecoli")
        c.record_codon_decision(CodonDecision(
            position=0, amino_acid="L", original_codon=None,
            chosen_codon="CTG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
        ))
        trail = c.finalize(output_dna="CTGCTG", cai=0.88, gc=0.67)
        assert len(trail.codon_decisions) == 1

    def test_record_wrong_type_raises(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        with pytest.raises(TypeError, match="Expected CodonDecision"):
            c.record_codon_decision("not a decision")

    def test_record_constraint_wrong_type_raises(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        with pytest.raises(TypeError, match="Expected ConstraintDecision"):
            c.record_constraint_decision("not a constraint")

    def test_to_json_before_finalize_raises(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        with pytest.raises(RuntimeError, match="finalize"):
            c.to_json("/tmp/test_trail.json")

    def test_to_json_writes_file(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(
            protein="MV", organism="Homo_sapiens",
            gene_name="TestGene", solver_backend="greedy",
        )
        c.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
        ))
        trail = c.finalize(output_dna="ATG", cai=1.0, gc=1.0)

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            filepath = f.name
        try:
            c.to_json(filepath)
            with open(filepath) as f:
                data = json.load(f)
            assert data["gene_name"] == "TestGene"
            assert len(data["codon_decisions"]) == 1
        finally:
            os.unlink(filepath)

    def test_summary_not_started(self):
        c = DecisionProvenanceCollector()
        assert "No optimization" in c.summary()

    def test_summary_active(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(
            protein="MV", organism="Ecoli", constraints=["GCInRange"],
            gene_name="TestGene",
        )
        summary = c.summary()
        assert "TestGene" in summary
        assert "Ecoli" in summary
        assert "IN PROGRESS" in summary

    def test_summary_finalized(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        c.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
        ))
        c.record_constraint_decision(ConstraintDecision(
            constraint_name="GCInRange",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[0],
            tradeoff_description="GC ok",
            impact_on_cai=-0.01,
        ))
        c.finalize(output_dna="ATG", cai=0.95, gc=0.5)
        summary = c.summary()
        assert "FINALIZED" in summary
        assert "GCInRange" in summary

    def test_repr_active(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        r = repr(c)
        assert "active" in r

    def test_repr_finalized(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        c.finalize(output_dna="ATG", cai=1.0, gc=1.0)
        r = repr(c)
        assert "finalized" in r

    def test_summary_with_iteration_log(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        c.record_iteration({"step": 1, "action": "assign", "score": 0.9})
        c.record_iteration({"step": 2, "action": "optimize", "score": 0.95})
        c.finalize(output_dna="ATG", cai=0.95, gc=0.5)
        summary = c.summary()
        assert "Iterations:" in summary
        assert "2" in summary

    def test_summary_low_confidence(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        c.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=0.3,
        ))
        c.finalize(output_dna="ATG", cai=1.0, gc=1.0)
        summary = c.summary()
        assert "Low-confidence" in summary

    def test_summary_cai_cost(self):
        c = DecisionProvenanceCollector()
        c.start_optimization(protein="MV", organism="Ecoli")
        c.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="gc_content", confidence=0.8,
            cai_impact=-0.05,
        ))
        c.finalize(output_dna="ATG", cai=0.95, gc=0.5)
        summary = c.summary()
        assert "CAI cost" in summary
