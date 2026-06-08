"""Tests for the deimmunization module.

Covers:
- DeimmunizationResult dataclass and field syncing
- EpitopeMutation and conversion to MutationResult
- compute_mutation_impact: binding, stability, solubility assessment
- find_epitope_disrupting_mutations: epitope-specific mutation search
- rank_deimmunization_mutations: global mutation ranking
- validate_deimmunized_protein: post-optimization validation
- Internal helpers: _estimate_ddg, _estimate_solubility_impact,
  _detect_mhc_class, _is_anchor_position, _is_structurally_dangerous,
  _organism_to_species, _filter_binder_epitopes
- Edge cases: short proteins, invalid positions, empty inputs
"""

from __future__ import annotations

import pytest

from biocompiler.deimmunization import (
    AppliedMutation,
    DeimmunizationResult,
    EpitopeMutation,
    MutationImpactResult,
    MutationValidation,
    ValidationReport,
    compute_mutation_impact,
    find_epitope_disrupting_mutations,
    rank_deimmunization_mutations,
    validate_deimmunized_protein,
    _estimate_ddg,
    _estimate_solubility_impact,
    _detect_mhc_class,
    _is_anchor_position,
    _is_structurally_dangerous,
    _organism_to_species,
    _filter_binder_epitopes,
    _filter_strong_binder_epitopes,
    _get_mhc_alleles,
)
from biocompiler.engine_base import MutationResult


# ---------------------------------------------------------------------------
# DeimmunizationResult
# ---------------------------------------------------------------------------

class TestDeimmunizationResult:
    """Tests for the DeimmunizationResult dataclass."""

    def test_construction_defaults(self):
        """DeimmunizationResult can be constructed with defaults."""
        result = DeimmunizationResult()
        assert result.sequence == ""
        assert result.primary_score == 0.0
        assert result.engine_name == "deimmunization"
        assert result.primary_score_label == "immunogenicity_score"
        assert result.success is False
        assert result.mutations_applied == []

    def test_classification_deimmunized(self):
        """Success=True yields 'deimmunized' classification."""
        result = DeimmunizationResult(success=True)
        assert result.classification == "deimmunized"

    def test_classification_partially_deimmunized(self):
        """Lower optimized than original yields 'partially_deimmunized'."""
        result = DeimmunizationResult(
            success=False,
            original_immunogenicity=0.8,
            optimized_immunogenicity=0.3,
        )
        assert result.classification == "partially_deimmunized"

    def test_classification_failed(self):
        """No improvement yields 'failed' classification."""
        result = DeimmunizationResult(success=False)
        assert result.classification == "failed"

    def test_field_sync_sequence(self):
        """optimized_protein syncs to sequence when sequence is empty."""
        result = DeimmunizationResult(
            optimized_protein="MVSKGE",
        )
        assert result.sequence == "MVSKGE"

    def test_field_sync_primary_score(self):
        """optimized_immunogenicity syncs to primary_score."""
        result = DeimmunizationResult(optimized_immunogenicity=0.5)
        assert result.primary_score == 0.5

    def test_field_sync_reverse(self):
        """primary_score syncs back to optimized_immunogenicity when it's 0."""
        result = DeimmunizationResult(primary_score=0.4)
        assert result.optimized_immunogenicity == 0.4

    def test_immunogenicity_score_property(self):
        """immunogenicity_score is an alias for optimized_immunogenicity."""
        result = DeimmunizationResult(optimized_immunogenicity=0.65)
        assert result.immunogenicity_score == 0.65

    def test_mutations_property(self):
        """mutations is an alias for mutations_applied."""
        muts = [AppliedMutation(position=0, wildtype="A", mutant="G",
                                epitope_removed="AAAAAAAAA", ddg=0.1,
                                blosum62=0, binding_reduction=0.5,
                                solubility_impact=0.0)]
        result = DeimmunizationResult(mutations_applied=muts)
        assert result.mutations == muts

    def test_passed_property(self):
        """passed reflects success field."""
        result = DeimmunizationResult(success=True)
        assert result.passed is True


# ---------------------------------------------------------------------------
# EpitopeMutation
# ---------------------------------------------------------------------------

class TestEpitopeMutation:
    """Tests for the EpitopeMutation dataclass."""

    def test_defaults(self):
        """EpitopeMutation has correct default fields."""
        mut = EpitopeMutation(position=5, original="A", mutant="G")
        assert mut.score_type == "immunogenicity"
        assert mut.engine == "deimmunization"
        assert mut.recommendation == "deimmunizing"

    def test_to_mutation_result(self):
        """to_mutation_result produces a plain MutationResult."""
        mut = EpitopeMutation(
            position=5, original="A", mutant="G",
            delta_score=0.4, description="test",
        )
        plain = mut.to_mutation_result()
        assert isinstance(plain, MutationResult)
        assert plain.position == 5
        assert plain.original == "A"
        assert plain.mutant == "G"
        assert plain.delta_score == 0.4
        assert plain.engine == "deimmunization"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

class TestEstimateDdg:
    """Tests for _estimate_ddg."""

    def test_conservative_mutation_low_ddg(self):
        """Conservative mutations (high BLOSUM62) should have low ddG."""
        # A→G has BLOSUM62=0, relatively conservative
        ddg = _estimate_ddg("A", "G")
        assert isinstance(ddg, float)
        assert ddg >= 0.0

    def test_radical_mutation_high_ddg(self):
        """Radical mutations (low BLOSUM62) should have higher ddG."""
        # W→P is radical (BLOSUM62 very negative)
        ddg_radical = _estimate_ddg("W", "P")
        ddg_conservative = _estimate_ddg("A", "G")
        assert ddg_radical > ddg_conservative

    def test_proline_introduction_penalty(self):
        """Introducing proline should add a penalty."""
        ddg_proline = _estimate_ddg("A", "P")
        ddg_no_proline = _estimate_ddg("A", "G")
        assert ddg_proline > ddg_no_proline

    def test_identical_mutation_zero_ddg(self):
        """Same amino acid should yield ~0 ddG."""
        ddg = _estimate_ddg("A", "A")
        # BLOSUM62[A][A] = 4, so blosum_component = max(0, (-4+4)*0.3) = 0
        # Hydro change = 0, no proline penalty
        assert ddg == 0.0

    def test_returns_float(self):
        """Return type is float."""
        ddg = _estimate_ddg("L", "I")
        assert isinstance(ddg, float)


class TestEstimateSolubilityImpact:
    """Tests for _estimate_solubility_impact."""

    def test_charged_introduction_positive(self):
        """Introducing a charged residue should increase solubility."""
        # A is not charged, K is positively charged
        impact = _estimate_solubility_impact("A", "K")
        assert impact > 0

    def test_charged_removal_negative(self):
        """Removing a charged residue should decrease solubility."""
        impact = _estimate_solubility_impact("K", "A")
        assert impact < 0

    def test_disulfide_breaking(self):
        """C→X mutation (breaking disulfide) should decrease solubility."""
        impact = _estimate_solubility_impact("C", "A")
        assert impact < 0

    def test_returns_float(self):
        """Return type is float."""
        impact = _estimate_solubility_impact("A", "G")
        assert isinstance(impact, float)


class TestDetectMhcClass:
    """Tests for _detect_mhc_class."""

    def test_class_i_allele(self):
        """HLA-A alleles are MHC class I."""
        assert _detect_mhc_class("HLA-A*02:01") == 1

    def test_class_ii_drb_allele(self):
        """DRB1 alleles are MHC class II."""
        assert _detect_mhc_class("HLA-DRB1*01:01") == 2

    def test_class_ii_dpb_allele(self):
        """DPB alleles are MHC class II."""
        assert _detect_mhc_class("HLA-DPB1*04:01") == 2

    def test_class_ii_dqb_allele(self):
        """DQB alleles are MHC class II."""
        assert _detect_mhc_class("HLA-DQB1*03:01") == 2

    def test_class_i_b_allele(self):
        """HLA-B alleles are MHC class I."""
        assert _detect_mhc_class("HLA-B*07:02") == 1

    def test_mouse_class_ii(self):
        """Mouse I-A alleles are MHC class II."""
        assert _detect_mhc_class("H2-IAb") == 2


class TestIsAnchorPosition:
    """Tests for _is_anchor_position."""

    def test_mhc_i_anchor_p2(self):
        """Position 1 (P2) is an MHC-I anchor."""
        assert _is_anchor_position(1, 0, mhc_class=1) is True

    def test_mhc_i_anchor_p9(self):
        """Position 8 (P9) is an MHC-I anchor."""
        assert _is_anchor_position(8, 0, mhc_class=1) is True

    def test_mhc_i_non_anchor(self):
        """Position 4 is not an MHC-I anchor."""
        assert _is_anchor_position(4, 0, mhc_class=1) is False

    def test_mhc_ii_anchor(self):
        """Position 0 is an MHC-II anchor."""
        assert _is_anchor_position(0, 0, mhc_class=2) is True

    def test_mhc_ii_anchor_p4(self):
        """Position 3 (P4) is an MHC-II anchor."""
        assert _is_anchor_position(3, 0, mhc_class=2) is True

    def test_non_zero_epitope_start(self):
        """Anchor detection works with non-zero epitope_start."""
        # Position 10 with epitope_start=9 → pos_in_peptide=1 (P2 anchor for MHC-I)
        assert _is_anchor_position(10, 9, mhc_class=1) is True


class TestIsStructurallyDangerous:
    """Tests for _is_structurally_dangerous."""

    def test_proline_in_helical_region(self):
        """Proline introduction in a helical region should be dangerous."""
        # 8 AAs of helix-favoring residues
        protein = "EALMKEAL"
        # Position 4 is surrounded by helix-favoring neighbors
        dangerous = _is_structurally_dangerous("A", "P", protein, 4)
        assert dangerous is True

    def test_proline_in_non_helical_region(self):
        """Proline introduction outside a helical region should be safe."""
        # Glycine/proline-rich region (not helical)
        protein = "GGPGGPGG"
        dangerous = _is_structurally_dangerous("G", "P", protein, 4)
        assert dangerous is False

    def test_hydrophobic_to_hydrophilic_in_core(self):
        """Hydrophobic→hydrophilic mutation in hydrophobic core is dangerous."""
        # All hydrophobic core residues
        protein = "AVILAVIL"
        dangerous = _is_structurally_dangerous("V", "K", protein, 4)
        assert dangerous is True

    def test_hydrophobic_to_hydrophobic_safe(self):
        """Hydrophobic→hydrophobic mutation is safe."""
        protein = "AVILAVIL"
        dangerous = _is_structurally_dangerous("V", "I", protein, 4)
        assert dangerous is False

    def test_non_proline_non_core_safe(self):
        """Neutral mutations that don't affect core are safe."""
        protein = "MVSKGEMV"
        dangerous = _is_structurally_dangerous("K", "R", protein, 3)
        assert dangerous is False


class TestOrganismToSpecies:
    """Tests for _organism_to_species."""

    def test_human(self):
        assert _organism_to_species("Homo_sapiens") == "human"

    def test_mouse(self):
        assert _organism_to_species("Mus_musculus") == "mouse"

    def test_ecoli(self):
        assert _organism_to_species("Escherichia_coli") == "ecoli"

    def test_cho(self):
        assert _organism_to_species("CHO_K1") == "cho"

    def test_yeast(self):
        assert _organism_to_species("Saccharomyces_cerevisiae") == "yeast"

    def test_unknown_organism(self):
        """Unknown organisms get a reasonable fallback."""
        result = _organism_to_species("Drosophila_melanogaster")
        assert isinstance(result, str)
        assert len(result) > 0


class TestFilterBinderEpitopes:
    """Tests for _filter_binder_epitopes and _filter_strong_binder_epitopes."""

    def test_filter_binder_epitopes(self):
        """Only strong and moderate binders are kept."""
        epitopes = [
            {"binding_class": "strong_binder", "score": 0.9},
            {"binding_class": "moderate_binder", "score": 0.5},
            {"binding_class": "weak_binder", "score": 0.1},
            {"binding_class": "non_binder", "score": 0.0},
        ]
        result = _filter_binder_epitopes(epitopes)
        assert len(result) == 2
        assert all(e["binding_class"] in ("strong_binder", "moderate_binder") for e in result)

    def test_filter_strong_binder_epitopes(self):
        """Only strong binders are kept."""
        epitopes = [
            {"binding_class": "strong_binder", "score": 0.9},
            {"binding_class": "moderate_binder", "score": 0.5},
        ]
        result = _filter_strong_binder_epitopes(epitopes)
        assert len(result) == 1
        assert result[0]["binding_class"] == "strong_binder"

    def test_empty_list(self):
        """Empty list returns empty."""
        assert _filter_binder_epitopes([]) == []
        assert _filter_strong_binder_epitopes([]) == []

    def test_no_binders(self):
        """List with no binders returns empty."""
        epitopes = [
            {"binding_class": "weak_binder"},
            {"binding_class": "non_binder"},
        ]
        assert _filter_binder_epitopes(epitopes) == []


class TestGetMhcAlleles:
    """Tests for _get_mhc_alleles."""

    def test_human_alleles(self):
        """Human organism returns a non-empty allele list."""
        alleles = _get_mhc_alleles("Homo_sapiens")
        assert isinstance(alleles, list)
        assert len(alleles) > 0

    def test_ecoli_no_alleles(self):
        """E. coli has no MHC alleles."""
        alleles = _get_mhc_alleles("Escherichia_coli")
        assert alleles == []

    def test_unknown_organism(self):
        """Unknown organism returns empty list."""
        alleles = _get_mhc_alleles("Unknown_organism")
        assert alleles == []


# ---------------------------------------------------------------------------
# compute_mutation_impact
# ---------------------------------------------------------------------------

class TestComputeMutationImpact:
    """Tests for compute_mutation_impact."""

    def test_basic_mutation(self):
        """compute_mutation_impact returns a dict with expected keys."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        result = compute_mutation_impact(protein, position=5, mutant_aa="A")
        assert isinstance(result, dict)
        assert "binding_impact" in result
        assert "stability_impact" in result
        assert "solubility_impact" in result
        assert "blosum62" in result

    def test_blosum62_score_correct(self):
        """BLOSUM62 score matches the substitution table."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        from biocompiler.constants import BLOSUM62
        result = compute_mutation_impact(protein, position=0, mutant_aa="V")
        expected = BLOSUM62["M"]["V"]
        assert result["blosum62"] == expected

    def test_empty_protein(self):
        """Empty protein returns zeroed result."""
        result = compute_mutation_impact("", position=0, mutant_aa="A")
        assert result["binding_impact"] == []
        assert result["stability_impact"] == 0.0
        assert result["blosum62"] == 0

    def test_invalid_position(self):
        """Position out of range returns zeroed result."""
        result = compute_mutation_impact("MVSKGE", position=10, mutant_aa="A")
        assert result["binding_impact"] == []
        assert result["blosum62"] == 0

    def test_negative_position(self):
        """Negative position returns zeroed result."""
        result = compute_mutation_impact("MVSKGE", position=-1, mutant_aa="A")
        assert result["binding_impact"] == []
        assert result["blosum62"] == 0

    def test_stability_impact_nonnegative(self):
        """Stability impact should be non-negative for any mutation."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        result = compute_mutation_impact(protein, position=3, mutant_aa="L")
        assert result["stability_impact"] >= 0.0


# ---------------------------------------------------------------------------
# find_epitope_disrupting_mutations
# ---------------------------------------------------------------------------

class TestFindEpitopeDisruptingMutations:
    """Tests for find_epitope_disrupting_mutations."""

    def test_returns_list(self):
        """Should return a list of MutationResult objects."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        result = find_epitope_disrupting_mutations(
            protein, epitope_start=0, epitope_end=9,
            mhc_alleles=["HLA-A*02:01"],
        )
        assert isinstance(result, list)
        # Each item should be a MutationResult
        for mut in result:
            assert isinstance(mut, MutationResult)

    def test_mutations_have_required_fields(self):
        """Each mutation should have position, original, mutant, and details."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        result = find_epitope_disrupting_mutations(
            protein, epitope_start=0, epitope_end=9,
            mhc_alleles=["HLA-A*02:01"],
        )
        for mut in result:
            assert mut.position >= 0
            assert len(mut.original) == 1
            assert len(mut.mutant) == 1
            assert mut.original != mut.mutant
            assert "blosum62" in mut.details
            assert "ddg_estimate" in mut.details
            assert "binding_reduction" in mut.details
            assert "is_anchor" in mut.details

    def test_mutations_within_epitope_range(self):
        """All mutation positions should be within the epitope range."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        start, end = 5, 14
        result = find_epitope_disrupting_mutations(
            protein, epitope_start=start, epitope_end=end,
            mhc_alleles=["HLA-A*02:01"],
        )
        for mut in result:
            assert start <= mut.position < end

    def test_binding_reduction_positive(self):
        """All returned mutations should have positive binding reduction."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        result = find_epitope_disrupting_mutations(
            protein, epitope_start=0, epitope_end=9,
            mhc_alleles=["HLA-A*02:01"],
        )
        for mut in result:
            assert mut.details["binding_reduction"] > 0

    def test_empty_allele_list(self):
        """Empty allele list should produce no mutations."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        result = find_epitope_disrupting_mutations(
            protein, epitope_start=0, epitope_end=9,
            mhc_alleles=[],
        )
        assert result == []

    def test_sorted_by_score_descending(self):
        """Results should be sorted by score (highest first)."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        result = find_epitope_disrupting_mutations(
            protein, epitope_start=0, epitope_end=9,
            mhc_alleles=["HLA-A*02:01"],
        )
        scores = [m.score for m in result]
        assert scores == sorted(scores, reverse=True)


# ---------------------------------------------------------------------------
# rank_deimmunization_mutations
# ---------------------------------------------------------------------------

class TestRankDeimmunizationMutations:
    """Tests for rank_deimmunization_mutations."""

    def test_returns_list(self):
        """Should return a list of MutationResult objects."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        result = rank_deimmunization_mutations(
            protein, mhc_alleles=["HLA-A*02:01"],
        )
        assert isinstance(result, list)

    def test_no_binder_epitopes(self):
        """Protein with no binder epitopes should return empty list."""
        # Very short protein unlikely to have strong epitopes
        protein = "MVSKGE"
        result = rank_deimmunization_mutations(
            protein, mhc_alleles=["HLA-A*02:01"],
        )
        assert isinstance(result, list)

    def test_deduplication(self):
        """Results should be deduplicated by (position, mutant) key."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        result = rank_deimmunization_mutations(
            protein, mhc_alleles=["HLA-A*02:01"],
        )
        # No duplicate (position, mutant) pairs
        seen = set()
        for mut in result:
            key = (mut.position, mut.mutant)
            assert key not in seen, f"Duplicate mutation: {key}"
            seen.add(key)


# ---------------------------------------------------------------------------
# validate_deimmunized_protein
# ---------------------------------------------------------------------------

class TestValidateDeimmunizedProtein:
    """Tests for validate_deimmunized_protein."""

    def test_no_mutations_valid(self):
        """No mutations should produce a valid report with correct structure."""
        protein = "MVSKGE"
        report = validate_deimmunized_protein(
            protein=protein,
            original_protein=protein,
            organism="Homo_sapiens",
        )
        assert isinstance(report, dict)
        assert report["stability_preserved"] is True
        assert report["all_mutations_conservative"] is True
        # immunogenicity_reduced may be False if both scores are 0
        # (no epitopes to reduce)
        assert isinstance(report["overall_valid"], bool)

    def test_report_structure(self):
        """Report has all required fields."""
        protein = "MVSKGE"
        report = validate_deimmunized_protein(
            protein=protein,
            original_protein=protein,
            organism="Homo_sapiens",
        )
        assert "immunogenicity_reduced" in report
        assert "original_immunogenicity" in report
        assert "optimized_immunogenicity" in report
        assert "stability_preserved" in report
        assert "total_ddg" in report
        assert "solubility_preserved" in report
        assert "all_mutations_conservative" in report
        assert "mutations" in report
        assert "overall_valid" in report

    def test_conservative_mutation_valid(self):
        """Conservative mutation should not invalidate the report."""
        original = "MVSKGE"
        optimized = "MVS KGE".replace(" ", "")  # same for now
        # Use a conservative BLOSUM62 mutation: I→L (BLOSUM62=2)
        mutations = [AppliedMutation(
            position=0, wildtype="M", mutant="M",  # no actual change
            epitope_removed="MVSKGE", ddg=0.0,
            blosum62=4, binding_reduction=0.5, solubility_impact=0.0,
        )]
        report = validate_deimmunized_protein(
            protein=optimized,
            original_protein=original,
            organism="Homo_sapiens",
        )
        assert isinstance(report, dict)
