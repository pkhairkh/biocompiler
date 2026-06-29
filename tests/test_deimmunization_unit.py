"""Unit tests for the BioCompiler Deimmunization Engine.

Covers four key areas:
  1. Basic deimmunization function exists and returns results
  2. Immunogenicity scoring with known sequences
  3. Mutation proposal reduces immunogenicity
  4. Constraint preservation during deimmunization
"""

from __future__ import annotations

import pytest

from biocompiler.immunogenicity.deimmunization import (
    DeimmunizationResult,
    EpitopeMutation,
    MutationImpactResult,
    ValidationReport,
    compute_mutation_impact,
    deimmunize,
    find_epitope_disrupting_mutations,
    rank_deimmunization_mutations,
    validate_deimmunized_protein,
)
from biocompiler.engines.base import MutationResult
from biocompiler.shared.exceptions import ImmunogenicityError


# ---------------------------------------------------------------------------
# Test proteins
# ---------------------------------------------------------------------------

# Short immunogenic peptide: repeated influenza M1 epitope (known HLA-A*02:01 binder)
# 18 AA — small enough for fast unit tests, large enough for 9-mer scanning
SHORT_IMMUNOGENIC = "GILGFVFTLGILGFVFTL"

# Human hemoglobin beta chain (well-characterised, moderate epitope load)
HBB_HUMAN = (
    "MVHLTPEEKSAVTALWGKVNVADIVGHALSDLHAKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
)


# ═══════════════════════════════════════════════════════════════════════════
# 1. Basic deimmunization function exists and returns results
# ═══════════════════════════════════════════════════════════════════════════

class TestDeimmunizeBasic:
    """Basic existence and return-type tests for the deimmunize() function."""

    def test_deimmunize_returns_result(self):
        """deimmunize() returns a DeimmunizationResult dataclass."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=1)
        assert isinstance(result, DeimmunizationResult)

    def test_deimmunize_result_has_required_fields(self):
        """Result must contain all documented domain-specific fields."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=1)
        assert hasattr(result, "original_protein")
        assert hasattr(result, "optimized_protein")
        assert hasattr(result, "mutations_applied")
        assert hasattr(result, "original_immunogenicity")
        assert hasattr(result, "optimized_immunogenicity")
        assert hasattr(result, "original_t_cell_epitopes")
        assert hasattr(result, "optimized_t_cell_epitopes")
        assert hasattr(result, "stability_preserved")
        assert hasattr(result, "iterations")
        assert hasattr(result, "method")
        assert hasattr(result, "success")

    def test_deimmunize_preserves_length(self):
        """Deimmunized protein must be the same length as the original."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=2)
        assert len(result.optimized_protein) == len(result.original_protein), (
            f"Length changed: {len(result.original_protein)} -> {len(result.optimized_protein)}"
        )

    def test_deimmunize_method_label(self):
        """Result method field should be 'iterative_epitope_disruption'."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=1)
        assert result.method == "iterative_epitope_disruption"

    def test_deimmunize_execution_time_recorded(self):
        """Result should have a non-negative execution time."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=1)
        assert result.execution_time_s >= 0.0

    def test_deimmunize_empty_sequence_raises(self):
        """Empty protein sequence must raise ImmunogenicityError."""
        with pytest.raises(ImmunogenicityError):
            deimmunize("")

    def test_deimmunize_invalid_sequence_raises(self):
        """Protein with non-standard amino acids must raise ImmunogenicityError."""
        with pytest.raises(ImmunogenicityError):
            deimmunize("MABCXYZ123")

    def test_deimmunize_zero_max_mutations(self):
        """With max_mutations=0, no mutations should be applied."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=0)
        assert result.iterations == 0
        assert len(result.mutations_applied) == 0
        # Optimized == original since no mutations were applied
        assert result.optimized_protein == result.original_protein

    def test_deimmunize_result_inherits_base_engine(self):
        """DeimmunizationResult should inherit from BaseEngineResult."""
        from biocompiler.engines.base import BaseEngineResult
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=1)
        assert isinstance(result, BaseEngineResult)

    def test_deimmunize_result_unified_api(self):
        """Unified API aliases should be accessible on result."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=1)
        # immunogenicity_score is a property alias for optimized_immunogenicity
        assert result.immunogenicity_score == result.optimized_immunogenicity
        # mutations is a property alias for mutations_applied
        assert result.mutations == result.mutations_applied


# ═══════════════════════════════════════════════════════════════════════════
# 2. Immunogenicity scoring with known sequences
# ═══════════════════════════════════════════════════════════════════════════

class TestImmunogenicityScoring:
    """Tests for immunogenicity scoring within the deimmunization module."""

    def test_score_in_range(self):
        """Immunogenicity score returned by deimmunize must be in [0, 1]."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=1)
        assert 0.0 <= result.original_immunogenicity <= 1.0, (
            f"Original score {result.original_immunogenicity} out of [0,1]"
        )
        assert 0.0 <= result.optimized_immunogenicity <= 1.0, (
            f"Optimized score {result.optimized_immunogenicity} out of [0,1]"
        )

    def test_score_range_hbb(self):
        """HBB protein immunogenicity score must be in [0, 1]."""
        result = deimmunize(HBB_HUMAN, max_mutations=0)
        assert 0.0 <= result.original_immunogenicity <= 1.0

    def test_compute_mutation_impact_returns_typed_dict(self):
        """compute_mutation_impact must return a MutationImpactResult dict."""
        impact = compute_mutation_impact(SHORT_IMMUNOGENIC, 5, "A")
        assert isinstance(impact, dict)
        assert "binding_impact" in impact
        assert "stability_impact" in impact
        assert "solubility_impact" in impact
        assert "blosum62" in impact

    def test_compute_mutation_impact_blosum62_correct(self):
        """BLOSUM62 field should match the known substitution score."""
        # BLOSUM62['G']['A'] = 0 (from the standard matrix)
        impact = compute_mutation_impact(SHORT_IMMUNOGENIC, 0, "A")  # G -> A
        assert impact["blosum62"] == 0, f"Expected BLOSUM62(G,A)=0, got {impact['blosum62']}"

    def test_estimate_ddg_blosum_heuristic(self):
        """BLOSUM62-based ddG heuristic (no FoldX) should be non-negative
        for a non-conservative substitution."""
        from biocompiler.immunogenicity.deimmunization import _estimate_ddg
        # Without protein context, the BLOSUM62 heuristic is used
        # W->P: BLOSUM62=-4, large hydropathy change, proline penalty
        ddg = _estimate_ddg("W", "P")
        assert ddg >= 0.0, f"Heuristic ddG for W->P should be >= 0, got {ddg}"

    def test_compute_mutation_impact_stability_is_numeric(self):
        """Stability impact (estimated ddG) should be a numeric value."""
        impact = compute_mutation_impact(SHORT_IMMUNOGENIC, 0, "P")  # G -> P (non-conservative)
        assert isinstance(impact["stability_impact"], (int, float)), (
            f"ddG should be numeric, got {type(impact['stability_impact'])}"
        )
        # FoldX delegation may return negative (stabilising) values,
        # so we just verify it is a reasonable magnitude
        assert abs(impact["stability_impact"]) < 100.0, (
            f"ddG {impact['stability_impact']} seems unreasonably large"
        )

    def test_compute_mutation_impact_empty_protein(self):
        """Empty protein should return zeroed impact result."""
        impact = compute_mutation_impact("", 0, "A")
        assert impact["binding_impact"] == []
        assert impact["stability_impact"] == 0.0
        assert impact["blosum62"] == 0

    def test_compute_mutation_impact_out_of_range(self):
        """Out-of-range position should return zeroed impact result."""
        impact = compute_mutation_impact(SHORT_IMMUNOGENIC, 999, "A")
        assert impact["binding_impact"] == []
        assert impact["stability_impact"] == 0.0

    def test_epitope_count_non_negative(self):
        """T-cell epitope counts should be non-negative integers."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=1)
        assert result.original_t_cell_epitopes >= 0
        assert result.optimized_t_cell_epitopes >= 0


# ═══════════════════════════════════════════════════════════════════════════
# 3. Test mutation proposal reduces immunogenicity
# ═══════════════════════════════════════════════════════════════════════════

class TestMutationReducesImmunogenicity:
    """Tests that proposed mutations actually reduce immunogenicity."""

    def test_deimmunize_reduces_or_maintains_score(self):
        """After deimmunization, the score should not increase significantly."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=2)
        assert result.optimized_immunogenicity <= result.original_immunogenicity + 0.05, (
            f"Immunogenicity increased: {result.original_immunogenicity:.4f} -> "
            f"{result.optimized_immunogenicity:.4f}"
        )

    def test_find_epitope_disrupting_mutations_returns_list(self):
        """find_epitope_disrupting_mutations returns a list of MutationResult."""
        mutations = find_epitope_disrupting_mutations(
            SHORT_IMMUNOGENIC, 0, 9, mhc_alleles=["HLA-A*02:01"]
        )
        assert isinstance(mutations, list)
        if mutations:
            assert isinstance(mutations[0], MutationResult)

    def test_find_epitope_mutations_have_positive_score(self):
        """All proposed mutations should have positive binding reduction (delta_score)."""
        mutations = find_epitope_disrupting_mutations(
            SHORT_IMMUNOGENIC, 0, 9, mhc_alleles=["HLA-A*02:01"]
        )
        for mut in mutations:
            assert mut.score > 0, (
                f"Mutation {mut.original}{mut.position}{mut.mutant} has non-positive score {mut.score}"
            )

    def test_find_epitope_mutations_sorted_descending(self):
        """Proposed mutations should be sorted by score descending."""
        mutations = find_epitope_disrupting_mutations(
            SHORT_IMMUNOGENIC, 0, 9, mhc_alleles=["HLA-A*02:01"]
        )
        if len(mutations) >= 2:
            for i in range(len(mutations) - 1):
                assert mutations[i].score >= mutations[i + 1].score, (
                    f"Mutations not sorted: {mutations[i].score} < {mutations[i+1].score}"
                )

    def test_rank_deimmunization_mutations_returns_list(self):
        """rank_deimmunization_mutations returns a list of MutationResult."""
        mutations = rank_deimmunization_mutations(
            SHORT_IMMUNOGENIC, mhc_alleles=["HLA-A*02:01"]
        )
        assert isinstance(mutations, list)
        if mutations:
            assert isinstance(mutations[0], MutationResult)

    def test_rank_mutations_deduplicated(self):
        """rank_deimmunization_mutations should not contain duplicate
        (position, mutant) pairs."""
        mutations = rank_deimmunization_mutations(
            SHORT_IMMUNOGENIC, mhc_alleles=["HLA-A*02:01"]
        )
        seen = set()
        for mut in mutations:
            key = (mut.position, mut.mutant)
            assert key not in seen, f"Duplicate mutation at position {mut.position}: {mut.mutant}"
            seen.add(key)

    def test_epitope_mutation_has_correct_engine_label(self):
        """EpitopeMutation objects should have engine='deimmunization'."""
        mutations = find_epitope_disrupting_mutations(
            SHORT_IMMUNOGENIC, 0, 9, mhc_alleles=["HLA-A*02:01"]
        )
        if mutations:
            epimut = mutations[0]
            # Could be EpitopeMutation subclass
            assert epimut.engine == "deimmunization" or hasattr(epimut, "to_mutation_result")

    def test_mutations_respect_blosum62_min_filter(self):
        """When blosum62_min > 0, all proposed mutations should have BLOSUM62 >= that."""
        blosum_min = 1
        mutations = find_epitope_disrupting_mutations(
            SHORT_IMMUNOGENIC, 0, 9,
            mhc_alleles=["HLA-A*02:01"],
            blosum62_min=blosum_min,
        )
        for mut in mutations:
            blosum_val = mut.details.get("blosum62", 0)
            assert blosum_val >= blosum_min, (
                f"Mutation {mut.original}{mut.position}{mut.mutant} has BLOSUM62={blosum_val} "
                f"< min={blosum_min}"
            )

    def test_deimmunize_with_target_score_zero(self):
        """Target score of 0.0 should push algorithm to apply as many
        mutations as possible up to max_mutations."""
        result = deimmunize(SHORT_IMMUNOGENIC, target_score=0.0, max_mutations=2)
        # Should attempt mutations (if any are available)
        assert isinstance(result, DeimmunizationResult)
        assert result.iterations <= 2


# ═══════════════════════════════════════════════════════════════════════════
# 4. Constraint preservation during deimmunization
# ═══════════════════════════════════════════════════════════════════════════

class TestConstraintPreservation:
    """Tests that deimmunization preserves structural and functional
    constraints (stability, solubility, BLOSUM62 conservation)."""

    def test_stability_preserved_default(self):
        """With default max_ddg=2.0, each mutation should have ddg < max_ddg."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=2, max_ddg=2.0)
        for mut in result.mutations_applied:
            assert mut["ddg"] < 2.0, (
                f"Mutation at pos {mut['position']} has ddg={mut['ddg']:.3f} >= max_ddg=2.0"
            )

    def test_cumulative_ddg_bounded(self):
        """Cumulative ddG across all mutations should be bounded."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=3, max_ddg=2.0)
        total_ddg = sum(mut["ddg"] for mut in result.mutations_applied)
        # The cumulative check in deimmunize uses max_ddg * _CUMULATIVE_DDG_MULTIPLIER (2)
        assert total_ddg < 2.0 * 2 + 0.01, (
            f"Cumulative ddG {total_ddg:.3f} exceeds threshold"
        )

    def test_preserve_positions_not_mutated(self):
        """Positions in preserve_positions should never be mutated."""
        preserve = [0, 1, 2]
        result = deimmunize(
            SHORT_IMMUNOGENIC, max_mutations=3, preserve_positions=preserve
        )
        for mut in result.mutations_applied:
            assert mut["position"] not in preserve, (
                f"Position {mut['position']} was mutated despite being in preserve_positions"
            )

    def test_blosum62_filter_enforced(self):
        """With blosum62_min=0, all mutations should have BLOSUM62 >= 0."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=2, blosum62_min=0)
        for mut in result.mutations_applied:
            assert mut["blosum62"] >= 0, (
                f"Mutation at pos {mut['position']} has BLOSUM62={mut['blosum62']} < 0"
            )

    def test_blosum62_strict_filter(self):
        """With blosum62_min=2, only very conservative mutations should be applied."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=2, blosum62_min=2)
        for mut in result.mutations_applied:
            assert mut["blosum62"] >= 2, (
                f"Mutation at pos {mut['position']} has BLOSUM62={mut['blosum62']} < 2"
            )

    def test_no_position_mutated_twice(self):
        """Each position should be mutated at most once."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=5)
        positions = [mut["position"] for mut in result.mutations_applied]
        assert len(positions) == len(set(positions)), (
            f"Duplicate positions in mutations: {positions}"
        )

    def test_max_mutations_respected(self):
        """Number of mutations should not exceed max_mutations."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=2)
        assert len(result.mutations_applied) <= 2, (
            f"Applied {len(result.mutations_applied)} mutations, max is 2"
        )

    def test_validate_deimmunized_protein_returns_report(self):
        """validate_deimmunized_protein returns a ValidationReport dict."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=2)
        report = validate_deimmunized_protein(
            result.optimized_protein, result.original_protein
        )
        assert isinstance(report, dict)
        assert "immunogenicity_reduced" in report
        assert "stability_preserved" in report
        assert "solubility_preserved" in report
        assert "all_mutations_conservative" in report
        assert "overall_valid" in report
        assert "mutations" in report

    def test_validate_report_scores_in_range(self):
        """Validation report immunogenicity scores should be in [0, 1]."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=1)
        report = validate_deimmunized_protein(
            result.optimized_protein, result.original_protein
        )
        assert 0.0 <= report["original_immunogenicity"] <= 1.0
        assert 0.0 <= report["optimized_immunogenicity"] <= 1.0

    def test_validate_report_mutation_details(self):
        """Each mutation in the validation report should have required keys."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=2)
        if not result.mutations_applied:
            pytest.skip("No mutations applied — cannot validate mutation details")
        report = validate_deimmunized_protein(
            result.optimized_protein, result.original_protein
        )
        if report["mutations"]:
            mut = report["mutations"][0]
            assert "position" in mut
            assert "wildtype" in mut
            assert "mutant" in mut
            assert "blosum62" in mut
            assert "ddg_estimate" in mut
            assert "conservative" in mut
            assert isinstance(mut["conservative"], bool)

    def test_validate_no_mutations_report(self):
        """Validating identical proteins (no mutations) should produce a
        report with empty mutations list."""
        report = validate_deimmunized_protein(SHORT_IMMUNOGENIC, SHORT_IMMUNOGENIC)
        assert report["mutations"] == []
        assert report["total_ddg"] == 0.0

    def test_validate_invalid_protein_raises(self):
        """Invalid protein sequence should raise ImmunogenicityError."""
        with pytest.raises(ImmunogenicityError):
            validate_deimmunized_protein("INVALID123", SHORT_IMMUNOGENIC)

    def test_max_ddg_tight_constraint(self):
        """With very tight max_ddg, mutations should still be bounded."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=2, max_ddg=0.5)
        for mut in result.mutations_applied:
            assert mut["ddg"] < 0.5, (
                f"Mutation at pos {mut['position']} has ddg={mut['ddg']:.3f} >= max_ddg=0.5"
            )

    def test_solubility_impact_recorded(self):
        """Each mutation in the result should have a solubility_impact field."""
        result = deimmunize(SHORT_IMMUNOGENIC, max_mutations=2)
        for mut in result.mutations_applied:
            assert "solubility_impact" in mut, (
                f"Mutation at pos {mut['position']} missing solubility_impact"
            )

    def test_epitope_mutation_to_mutation_result(self):
        """EpitopeMutation.to_mutation_result() should produce a valid MutationResult."""
        mutations = find_epitope_disrupting_mutations(
            SHORT_IMMUNOGENIC, 0, 9, mhc_alleles=["HLA-A*02:01"]
        )
        for mut in mutations:
            if isinstance(mut, EpitopeMutation) and hasattr(mut, "to_mutation_result"):
                converted = mut.to_mutation_result()
                assert isinstance(converted, MutationResult)
                assert converted.position == mut.position
                assert converted.original == mut.original
                assert converted.mutant == mut.mutant
                break  # Test at least one conversion
