"""Comprehensive tests for the unified engine API.

Tests BaseEngineResult, MutationResult, BatchResult, EngineTimer,
EngineConfig, classify_score, validate_protein_sequence, and the
inheritance / property-alias integration across all engine result types.

Only this file is created/modified — no other project files are touched.
"""

from __future__ import annotations

import time

import pytest

# ── Unified base types ──────────────────────────────────────────────────────
from biocompiler.engines.base import (
    BaseEngineResult,
    BatchResult,
    EngineConfig,
    EngineResult,
    EngineTimer,
    MutationResult,
    classify_score,
    validate_protein_sequence,
)

# ── Engine-specific result types ────────────────────────────────────────────
from biocompiler.engines.esmfold import ESMFoldResult
from biocompiler.engines.foldx import FoldXResult
from biocompiler.engines.camsol import CamSolResult
from biocompiler.immunogenicity.core import ImmunogenicityResult
from biocompiler.immunogenicity.deimmunization import DeimmunizationResult
from biocompiler.engines.protein_design import DesignResult


# ═══════════════════════════════════════════════════════════════════════════
# 1. BaseEngineResult
# ═══════════════════════════════════════════════════════════════════════════

class TestBaseEngineResult:

    def test_construction_all_fields(self):
        r = BaseEngineResult(
            sequence="ACDE",
            primary_score=42.5,
            classification="good",
            success=True,
            error=None,
            execution_time_s=1.23,
            engine_name="test_engine",
            primary_score_label="test_metric",
        )
        assert r.sequence == "ACDE"
        assert r.primary_score == 42.5
        assert r.classification == "good"
        assert r.success is True
        assert r.error is None
        assert r.execution_time_s == 1.23
        assert r.engine_name == "test_engine"
        assert r.primary_score_label == "test_metric"

    def test_passed_property_returns_success(self):
        r_ok = BaseEngineResult(sequence="A", primary_score=1.0, classification="", success=True)
        r_fail = BaseEngineResult(sequence="A", primary_score=1.0, classification="", success=False)
        assert r_ok.passed is True
        assert r_fail.passed is False

    def test_default_values(self):
        r = BaseEngineResult(sequence="A", primary_score=0.0, classification="", success=True)
        assert r.error is None
        assert r.execution_time_s == 0.0
        assert r.engine_name == ""
        assert r.primary_score_label == "score"
        assert r.success is True

    def test_is_dataclass(self):
        """BaseEngineResult should be a dataclass (has __dataclass_fields__)."""
        assert hasattr(BaseEngineResult, "__dataclass_fields__")


# ═══════════════════════════════════════════════════════════════════════════
# 2. MutationResult
# ═══════════════════════════════════════════════════════════════════════════

class TestMutationResult:

    def test_construction_all_fields(self):
        m = MutationResult(
            position=5,
            original="A",
            mutant="G",
            delta_score=-1.2,
            score_type="ddg",
            engine="foldx",
            recommendation="stabilizing",
            description="A6G stabilizing mutation",
            details={"ddg_kcal": 1.2},
        )
        assert m.position == 5
        assert m.original == "A"
        assert m.mutant == "G"
        assert m.delta_score == -1.2
        assert m.score_type == "ddg"
        assert m.engine == "foldx"
        assert m.recommendation == "stabilizing"
        assert m.description == "A6G stabilizing mutation"
        assert m.details == {"ddg_kcal": 1.2}

    def test_score_property_alias(self):
        m = MutationResult(position=0, original="A", mutant="G",
                           delta_score=3.5, score_type="solubility", engine="camsol")
        assert m.score == 3.5
        assert m.score == m.delta_score

    def test_score_setter_alias(self):
        m = MutationResult(position=0, original="A", mutant="G",
                           delta_score=0.0, score_type="ddg", engine="foldx")
        m.score = 9.9
        assert m.delta_score == 9.9

    def test_original_aa_property_alias(self):
        m = MutationResult(position=0, original="L", mutant="P",
                           delta_score=-1.0, score_type="ddg", engine="foldx")
        assert m.original_aa == "L"

    def test_mutant_aa_property_alias(self):
        m = MutationResult(position=0, original="L", mutant="P",
                           delta_score=-1.0, score_type="ddg", engine="foldx")
        assert m.mutant_aa == "P"

    def test_str_representation(self):
        m = MutationResult(position=9, original="V", mutant="A",
                           delta_score=2.50, score_type="ddg", engine="foldx")
        s = str(m)
        # position+1 for 1-based display
        assert "V10A" in s
        assert "foldx" in s
        assert "ddg" in s

    def test_backward_compat_defaults(self):
        m = MutationResult(position=0, original="A", mutant="G",
                           delta_score=1.0, score_type="ddg", engine="foldx")
        assert m.recommendation == ""
        assert m.description == ""
        assert m.details == {}


# ═══════════════════════════════════════════════════════════════════════════
# 3. BatchResult
# ═══════════════════════════════════════════════════════════════════════════

class TestBatchResult:

    def test_auto_counting_successful_failed(self):
        r1 = BaseEngineResult(sequence="A", primary_score=1, classification="", success=True)
        r2 = BaseEngineResult(sequence="B", primary_score=2, classification="", success=True)
        r3 = BaseEngineResult(sequence="C", primary_score=0, classification="", success=False)
        br = BatchResult(results=[r1, r2, r3])
        assert br.successful == 2
        assert br.failed == 1

    def test_success_count_alias(self):
        r1 = BaseEngineResult(sequence="A", primary_score=1, classification="", success=True)
        br = BatchResult(results=[r1])
        assert br.success_count == br.successful

    def test_failure_count_alias(self):
        r1 = BaseEngineResult(sequence="A", primary_score=1, classification="", success=False)
        br = BatchResult(results=[r1])
        assert br.failure_count == br.failed

    def test_total_property(self):
        br = BatchResult(results=[1, 2, 3])  # generic T
        assert br.total == 3

    def test_empty_batch(self):
        br = BatchResult()
        assert br.results == []
        assert br.total == 0
        assert br.successful == 0
        assert br.failed == 0

    def test_explicit_counts_not_overridden(self):
        """If successful/failed are explicitly set, __post_init__ should not override."""
        r1 = BaseEngineResult(sequence="A", primary_score=1, classification="", success=True)
        br = BatchResult(results=[r1], successful=5, failed=2)
        assert br.successful == 5
        assert br.failed == 2

    def test_generic_type_parameter(self):
        """BatchResult should accept typed results."""
        esm1 = ESMFoldResult(protein="ACDE", primary_score=90.0, classification="Very high (experimental)")
        br: BatchResult[ESMFoldResult] = BatchResult(results=[esm1])
        assert br.total == 1
        assert isinstance(br.results[0], ESMFoldResult)


# ═══════════════════════════════════════════════════════════════════════════
# 4. EngineTimer
# ═══════════════════════════════════════════════════════════════════════════

class TestEngineTimer:

    def test_context_manager_usage(self):
        with EngineTimer() as timer:
            _ = sum(range(1000))
        assert timer.elapsed > 0

    def test_elapsed_positive_after_use(self):
        with EngineTimer() as timer:
            time.sleep(0.01)  # 10ms
        assert timer.elapsed >= 0.01

    def test_elapsed_zero_before_use(self):
        timer = EngineTimer()
        assert timer.elapsed == 0.0


# ═══════════════════════════════════════════════════════════════════════════
# 5. EngineConfig
# ═══════════════════════════════════════════════════════════════════════════

class TestEngineConfig:

    def test_default_values(self):
        cfg = EngineConfig()
        assert cfg.use_cache is True
        assert cfg.timeout_s == 300.0
        assert cfg.verbose is False
        assert cfg.max_workers == 4

    def test_custom_values(self):
        cfg = EngineConfig(use_cache=False, timeout_s=60.0, verbose=True, max_workers=8)
        assert cfg.use_cache is False
        assert cfg.timeout_s == 60.0
        assert cfg.verbose is True
        assert cfg.max_workers == 8


# ═══════════════════════════════════════════════════════════════════════════
# 6. classify_score
# ═══════════════════════════════════════════════════════════════════════════

class TestClassifyScore:

    def test_first_threshold_wins(self):
        result = classify_score(95.0, [(90, "high"), (70, "medium"), (50, "low")])
        assert result == "high"

    def test_middle_threshold(self):
        result = classify_score(75.0, [(90, "high"), (70, "medium"), (50, "low")])
        assert result == "medium"

    def test_last_threshold(self):
        result = classify_score(55.0, [(90, "high"), (70, "medium"), (50, "low")])
        assert result == "low"

    def test_below_all_thresholds_uses_fallback(self):
        result = classify_score(10.0, [(90, "high"), (70, "medium"), (50, "low")])
        assert result == "unknown"

    def test_custom_fallback(self):
        result = classify_score(10.0, [(90, "high")], fallback="very_low")
        assert result == "very_low"

    def test_exact_threshold_match(self):
        """Score exactly equal to threshold should match (>=)."""
        result = classify_score(90.0, [(90, "high"), (70, "medium")])
        assert result == "high"

    def test_empty_thresholds_returns_fallback(self):
        result = classify_score(50.0, [], fallback="none")
        assert result == "none"


# ═══════════════════════════════════════════════════════════════════════════
# 7. validate_protein_sequence
# ═══════════════════════════════════════════════════════════════════════════

class TestValidateProteinSequence:

    def test_valid_sequence_passes(self):
        result = validate_protein_sequence("ACDEFGHIKLMNPQRSTVWY")
        assert result == "ACDEFGHIKLMNPQRSTVWY"

    def test_valid_sequence_normalizes(self):
        result = validate_protein_sequence(" acd ", "test")
        assert result == "ACD"

    def test_invalid_characters_rejected(self):
        with pytest.raises(ValueError, match="invalid amino acids"):
            validate_protein_sequence("ACDX", "test")

    def test_empty_sequence_rejected(self):
        with pytest.raises(ValueError, match="must not be empty"):
            validate_protein_sequence("", "test")

    def test_whitespace_only_rejected(self):
        with pytest.raises(ValueError, match="must not be empty"):
            validate_protein_sequence("   ", "test")

    def test_engine_name_in_error(self):
        with pytest.raises(ValueError, match="MyEngine"):
            validate_protein_sequence("X", "MyEngine")


# ═══════════════════════════════════════════════════════════════════════════
# 8. Engine result type inheritance from BaseEngineResult
# ═══════════════════════════════════════════════════════════════════════════

class TestEngineResultInheritance:

    def test_esmfold_result_is_base_engine_result(self):
        r = ESMFoldResult(protein="ACDE", primary_score=90.0)
        assert isinstance(r, BaseEngineResult)

    def test_foldx_result_is_base_engine_result(self):
        r = FoldXResult(protein="ACDE", stability_kcal=-5.0)
        assert isinstance(r, BaseEngineResult)

    def test_camsol_result_is_base_engine_result(self):
        r = CamSolResult(sequence="ACDE", primary_score=1.5, classification="soluble", success=True)
        assert isinstance(r, BaseEngineResult)

    def test_immunogenicity_result_is_base_engine_result(self):
        r = ImmunogenicityResult(sequence="ACDE", primary_score=0.3, classification="low_risk")
        assert isinstance(r, BaseEngineResult)

    def test_deimmunization_result_is_base_engine_result(self):
        r = DeimmunizationResult(
            original_protein="ACDE", optimized_protein="ACDE",
            optimized_immunogenicity=0.1, success=True,
        )
        assert isinstance(r, BaseEngineResult)

    def test_design_result_is_base_engine_result(self):
        r = DesignResult(
            designed_protein="ACDE", stability_change=-2.0,
            success=True, constraints_satisfied=["min_stability"],
        )
        assert isinstance(r, BaseEngineResult)

    def test_all_satisfy_engine_result_protocol(self):
        """All engine result types should satisfy the EngineResult protocol."""
        results = [
            ESMFoldResult(protein="A", primary_score=50.0),
            FoldXResult(protein="A", stability_kcal=-5.0),
            CamSolResult(sequence="A", primary_score=1.0, classification="soluble", success=True),
            ImmunogenicityResult(sequence="A", primary_score=0.5, classification="medium_risk"),
            DeimmunizationResult(
                original_protein="A", optimized_protein="A",
                optimized_immunogenicity=0.1, success=True,
            ),
            DesignResult(
                designed_protein="A", stability_change=-1.0,
                success=True, constraints_satisfied=[],
            ),
        ]
        for r in results:
            assert isinstance(r, EngineResult), f"{type(r).__name__} does not satisfy EngineResult protocol"


# ═══════════════════════════════════════════════════════════════════════════
# 9. Property aliases on engine result types
# ═══════════════════════════════════════════════════════════════════════════

class TestESMFoldResultAliases:

    def test_plddt_alias_for_primary_score(self):
        r = ESMFoldResult(protein="ACDE", primary_score=88.5, classification="Confident")
        assert r.plddt == r.primary_score
        assert r.plddt == 88.5

    def test_plddt_setter(self):
        r = ESMFoldResult(protein="ACDE", primary_score=50.0)
        r.plddt = 92.0
        assert r.primary_score == 92.0

    def test_confidence_class_alias_for_classification(self):
        r = ESMFoldResult(protein="ACDE", primary_score=92.0, classification="Very high (experimental)")
        assert r.confidence_class == r.classification
        assert r.confidence_class == "Very high (experimental)"

    def test_confidence_class_setter(self):
        r = ESMFoldResult(protein="ACDE", primary_score=50.0)
        r.confidence_class = "Confident"
        assert r.classification == "Confident"

    def test_protein_alias_for_sequence(self):
        r = ESMFoldResult(protein="ACDE", primary_score=50.0)
        assert r.protein == r.sequence

    def test_protein_setter(self):
        r = ESMFoldResult(protein="ACDE", primary_score=50.0)
        r.protein = "FGHI"
        assert r.sequence == "FGHI"

    def test_mean_plddt_alias(self):
        r = ESMFoldResult(protein="ACDE", primary_score=75.0)
        assert r.mean_plddt == r.primary_score

    def test_passed_inherited(self):
        r = ESMFoldResult(protein="ACDE", primary_score=50.0, success=True)
        assert r.passed is True


class TestFoldXResultAliases:

    def test_ddg_alias_for_primary_score(self):
        r = FoldXResult(protein="ACDE", stability_kcal=-7.5)
        assert r.ddg == r.primary_score

    def test_stability_class_alias_for_classification(self):
        r = FoldXResult(protein="ACDE", stability_kcal=-7.5)
        assert r.stability_class == r.classification

    def test_stabilizing_mutations_alias_for_mutations(self):
        m = MutationResult(position=0, original="A", mutant="G",
                           delta_score=-1.0, score_type="ddg", engine="foldx")
        r = FoldXResult(protein="ACDE", stability_kcal=-5.0, mutations=[m])
        assert r.stabilizing_mutations == r.mutations
        assert len(r.stabilizing_mutations) == 1

    def test_passed_inherited(self):
        r = FoldXResult(protein="ACDE", stability_kcal=-5.0, success=True)
        assert r.passed is True


class TestCamSolResultAliases:

    def test_score_alias_for_primary_score(self):
        r = CamSolResult(sequence="ACDE", primary_score=1.8, classification="soluble", success=True)
        assert r.score == r.primary_score
        assert r.score == 1.8

    def test_solubility_class_alias_for_classification(self):
        r = CamSolResult(sequence="ACDE", primary_score=1.8, classification="highly_soluble", success=True)
        assert r.solubility_class == r.classification

    def test_solubility_class_setter(self):
        r = CamSolResult(sequence="ACDE", primary_score=1.8, classification="soluble", success=True)
        r.solubility_class = "insoluble"
        assert r.classification == "insoluble"

    def test_protein_alias_for_sequence(self):
        r = CamSolResult(sequence="ACDE", primary_score=1.0, classification="soluble", success=True)
        assert r.protein == r.sequence

    def test_passed_inherited(self):
        r = CamSolResult(sequence="ACDE", primary_score=1.0, classification="soluble", success=True)
        assert r.passed is True


class TestImmunogenicityResultAliases:

    def test_immunogenicity_score_alias_for_primary_score(self):
        r = ImmunogenicityResult(sequence="ACDE", primary_score=0.7, classification="high_risk")
        assert r.immunogenicity_score == r.primary_score

    def test_risk_class_alias_for_classification(self):
        r = ImmunogenicityResult(sequence="ACDE", primary_score=0.7, classification="high_risk")
        assert r.risk_class == r.classification

    def test_risk_class_setter(self):
        r = ImmunogenicityResult(sequence="ACDE", primary_score=0.7, classification="medium_risk")
        r.risk_class = "low_risk"
        assert r.classification == "low_risk"

    def test_immunogenicity_class_alias(self):
        r = ImmunogenicityResult(sequence="ACDE", primary_score=0.7, classification="high_risk")
        assert r.immunogenicity_class == r.classification

    def test_passed_inherited(self):
        r = ImmunogenicityResult(sequence="ACDE", primary_score=0.7, classification="high_risk", success=True)
        assert r.passed is True


class TestDeimmunizationResultAliases:

    def test_immunogenicity_score_alias(self):
        r = DeimmunizationResult(
            original_protein="ACDE", optimized_protein="ACDE",
            optimized_immunogenicity=0.2, success=True,
        )
        assert r.immunogenicity_score == r.optimized_immunogenicity

    def test_mutations_alias(self):
        r = DeimmunizationResult(
            original_protein="ACDE", optimized_protein="ACDE",
            optimized_immunogenicity=0.1, success=True,
            mutations_applied=[{"position": 0, "wildtype": "A", "mutant": "G"}],
        )
        assert r.mutations == r.mutations_applied

    def test_passed_inherited(self):
        r = DeimmunizationResult(
            original_protein="ACDE", optimized_protein="ACDE",
            optimized_immunogenicity=0.1, success=True,
        )
        assert r.passed is True


class TestDesignResultAliases:

    def test_designed_sequence_alias(self):
        r = DesignResult(
            designed_protein="ACDE", stability_change=-2.0,
            success=True, constraints_satisfied=[],
        )
        assert r.designed_sequence == r.designed_protein

    def test_passed_inherited(self):
        r = DesignResult(
            designed_protein="ACDE", stability_change=-2.0,
            success=True, constraints_satisfied=["min_stability"],
        )
        assert r.passed is True

    def test_classification_auto_computed(self):
        r_success = DesignResult(
            designed_protein="ACDE", stability_change=-2.0,
            success=True, constraints_satisfied=["min_stability"],
            constraints_violated=[],
        )
        assert r_success.classification == "design_success"

        r_partial = DesignResult(
            designed_protein="ACDE", stability_change=-2.0,
            success=True, constraints_satisfied=[],
            constraints_violated=["min_solubility"],
        )
        assert r_partial.classification == "design_partial"

        r_failed = DesignResult(
            designed_protein="ACDE", stability_change=0.0,
            success=False, constraints_satisfied=[],
            constraints_violated=["min_stability"],
        )
        assert r_failed.classification == "design_failed"
