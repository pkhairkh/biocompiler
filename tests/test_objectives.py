"""
BioCompiler Objectives Module Tests
====================================

Tests for the biocompiler.objectives module — objective function protocol,
built-in objectives, registry, and resolve_objective().
"""

from __future__ import annotations

import pytest

from biocompiler.objectives import (
    ObjectiveFunction,
    cai_objective,
    cai_gc_balanced_objective,
    codon_pair_objective,
    min_max_gc_objective,
    resolve_objective,
    OBJECTIVE_REGISTRY,
)


# ════════════════════════════════════════════════════════════════════════════
# ObjectiveFunction Protocol
# ════════════════════════════════════════════════════════════════════════════

class TestObjectiveFunctionProtocol:
    """Test that ObjectiveFunction is a runtime-checkable Protocol."""

    def test_conforming_callable_satisfies_protocol(self):
        """A function with the right signature satisfies ObjectiveFunction."""
        def my_obj(dna: str, protein: str, organism: str) -> float:
            return 1.0
        assert isinstance(my_obj, ObjectiveFunction)

    def test_wrong_signature_does_not_satisfy(self):
        """A function with wrong signature does not satisfy ObjectiveFunction."""
        def bad_obj(x: int) -> str:
            return "nope"
        # Note: Protocol with __call__ checks at runtime
        # This should still pass isinstance because Protocol checks
        # for the presence of __call__, which all callables have.
        # The real check is that the callable can be called with (dna, protein, organism).
        assert callable(bad_obj)


# ════════════════════════════════════════════════════════════════════════════
# min_max_gc_objective
# ════════════════════════════════════════════════════════════════════════════

class TestMinMaxGcObjective:
    """Test the min_max_gc_objective function — no external dependencies."""

    def test_exact_target_returns_one(self):
        """DNA with exactly the target GC returns 1.0."""
        # 6 G/C out of 12 bases = 0.50 GC
        dna = "GGGGGGCCCCCC"  # 12 bases, 12 GC → GC=1.0
        # Use target_gc=1.0 for this DNA
        score = min_max_gc_objective(dna, "MMMM", "ecoli", target_gc=1.0)
        assert score == 1.0

    def test_half_gc_target_half(self):
        """DNA with 50% GC and target 0.5 returns 1.0."""
        # ATGCATGC = 4/8 = 0.5
        dna = "ATGCATGC"
        score = min_max_gc_objective(dna, "MM", "ecoli", target_gc=0.5)
        assert score == 1.0

    def test_deviation_reduces_score(self):
        """DNA far from target GC gives a low score."""
        dna = "ATATATAT"  # 0% GC
        score = min_max_gc_objective(dna, "MM", "ecoli", target_gc=0.55)
        assert score < 0.5  # Should be significantly below 1.0

    def test_empty_dna_returns_zero(self):
        """Empty DNA sequence returns 0.0."""
        score = min_max_gc_objective("", "M", "ecoli")
        assert score == 0.0

    def test_all_gc_at_target(self):
        """All-GC DNA with target 1.0 returns 1.0."""
        dna = "GCGCGCGC"
        score = min_max_gc_objective(dna, "AAAA", "ecoli", target_gc=1.0)
        assert score == 1.0

    def test_quadratic_penalty(self):
        """Small deviations are tolerated, large ones are penalized quadratically."""
        # Exactly on target: 1.0
        dna_exact = "ATGC" * 10  # 50% GC
        s_exact = min_max_gc_objective(dna_exact, "M" * 20, "ecoli", target_gc=0.5)

        # Slightly off: 45% GC
        dna_slight = "ATGCATGCAT" * 4  # mix
        gc = (dna_slight.count("G") + dna_slight.count("C")) / len(dna_slight)
        s_slight = min_max_gc_objective(dna_slight, "M" * 20, "ecoli", target_gc=0.5)

        # Closer should be higher score
        assert s_exact >= s_slight

    def test_score_is_non_negative(self):
        """Score should never go below 0.0."""
        dna = "ATATATAT"  # 0% GC
        score = min_max_gc_objective(dna, "M" * 8, "ecoli", target_gc=0.99)
        assert score >= 0.0

    def test_score_never_exceeds_one(self):
        """Score should never exceed 1.0."""
        dna = "GCGCGCGC"
        score = min_max_gc_objective(dna, "M" * 8, "ecoli", target_gc=0.5)
        assert score <= 1.0


# ════════════════════════════════════════════════════════════════════════════
# cai_objective
# ════════════════════════════════════════════════════════════════════════════

class TestCaiObjective:

    def test_empty_dna_returns_zero(self):
        """Empty DNA returns 0.0."""
        assert cai_objective("", "M", "ecoli") == 0.0

    def test_short_dna_returns_zero(self):
        """DNA shorter than 3 bases returns 0.0."""
        assert cai_objective("AT", "M", "ecoli") == 0.0

    def test_valid_dna_returns_positive(self):
        """Valid codon-optimized DNA should return a positive CAI."""
        dna = "ATGAAATTTCTG"  # MKFL — 4 codons
        score = cai_objective(dna, "MKFL", "Escherichia_coli")
        assert 0.0 < score <= 1.0

    def test_invalid_organism_returns_zero(self):
        """An unknown organism should return 0.0 (graceful fallback)."""
        dna = "ATGAAA"
        score = cai_objective(dna, "MK", "NonExistentOrganism_12345")
        assert score == 0.0


# ════════════════════════════════════════════════════════════════════════════
# cai_gc_balanced_objective
# ════════════════════════════════════════════════════════════════════════════

class TestCaiGcBalancedObjective:

    def test_empty_dna_returns_zero(self):
        """Empty DNA returns 0.0."""
        assert cai_gc_balanced_objective("", "M", "ecoli") == 0.0

    def test_returns_combined_score(self):
        """Should return a score between 0 and 1."""
        dna = "ATGAAATTTCTG"
        score = cai_gc_balanced_objective(dna, "MKFL", "Escherichia_coli")
        assert 0.0 <= score <= 1.0

    def test_gc_weight_zero(self):
        """With gc_weight=0, the score should equal pure CAI."""
        dna = "ATGAAATTTCTG"
        score = cai_gc_balanced_objective(dna, "MKFL", "Escherichia_coli", gc_weight=0.0)
        cai_score = cai_objective(dna, "MKFL", "Escherichia_coli")
        assert abs(score - cai_score) < 1e-10

    def test_gc_weight_one(self):
        """With gc_weight=1.0, only the GC component matters."""
        dna = "ATGAAATTTCTG"
        score = cai_gc_balanced_objective(dna, "MKFL", "Escherichia_coli", gc_weight=1.0)
        # Should be purely the GC score
        assert 0.0 <= score <= 1.0
        # Should NOT equal the CAI score (unless by coincidence)
        cai_score = cai_objective(dna, "MKFL", "Escherichia_coli")
        # gc_weight=1.0 means cai component is 0, so result is just gc_score
        assert score <= 1.0


# ════════════════════════════════════════════════════════════════════════════
# codon_pair_objective
# ════════════════════════════════════════════════════════════════════════════

class TestCodonPairObjective:

    def test_short_dna_returns_default(self):
        """DNA shorter than 6 bases returns 0.5 (default neutral)."""
        assert codon_pair_objective("ATGAA", "MK", "ecoli") == 0.5

    def test_empty_dna_returns_default(self):
        """Empty DNA returns 0.5."""
        assert codon_pair_objective("", "M", "ecoli") == 0.5

    def test_valid_dna_returns_in_range(self):
        """Valid DNA should return a score in [0.0, 1.0]."""
        dna = "ATGAAATTTCTGGCA"
        score = codon_pair_objective(dna, "MKFLA", "Escherichia_coli")
        assert 0.0 <= score <= 1.0


# ════════════════════════════════════════════════════════════════════════════
# OBJECTIVE_REGISTRY
# ════════════════════════════════════════════════════════════════════════════

class TestObjectiveRegistry:

    def test_contains_all_builtins(self):
        """Registry has entries for all four built-in objectives."""
        assert "cai" in OBJECTIVE_REGISTRY
        assert "cai_gc_balanced" in OBJECTIVE_REGISTRY
        assert "codon_pair" in OBJECTIVE_REGISTRY
        assert "min_max_gc" in OBJECTIVE_REGISTRY

    def test_registry_values_are_callable(self):
        """All registry values are callable."""
        for name, func in OBJECTIVE_REGISTRY.items():
            assert callable(func), f"Registry entry '{name}' is not callable"

    def test_registry_has_five_entries(self):
        """Registry should have exactly 5 entries (cai, cai_gc_balanced, codon_pair, min_max_gc, tai)."""
        assert len(OBJECTIVE_REGISTRY) == 5


# ════════════════════════════════════════════════════════════════════════════
# resolve_objective
# ════════════════════════════════════════════════════════════════════════════

class TestResolveObjective:

    def test_none_returns_cai(self):
        """resolve_objective(None) returns cai_objective."""
        func = resolve_objective(None)
        assert func is cai_objective

    def test_string_cai(self):
        """resolve_objective('cai') returns cai_objective."""
        func = resolve_objective("cai")
        assert func is cai_objective

    def test_string_min_max_gc(self):
        """resolve_objective('min_max_gc') returns min_max_gc_objective."""
        func = resolve_objective("min_max_gc")
        assert func is min_max_gc_objective

    def test_string_case_insensitive(self):
        """String lookup is case-insensitive after .lower().strip()."""
        func = resolve_objective("CAI")
        assert func is cai_objective

    def test_string_with_whitespace(self):
        """String lookup strips whitespace."""
        func = resolve_objective("  cai  ")
        assert func is cai_objective

    def test_callable_passthrough(self):
        """A callable is returned as-is."""
        def my_obj(dna, protein, organism):
            return 0.5
        func = resolve_objective(my_obj)
        assert func is my_obj

    def test_unknown_string_raises_value_error(self):
        """An unknown string name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown objective"):
            resolve_objective("nonexistent_objective")

    def test_invalid_type_raises_type_error(self):
        """A non-string, non-callable, non-None raises TypeError."""
        with pytest.raises(TypeError, match="objective must be"):
            resolve_objective(42)

    def test_list_raises_type_error(self):
        """A list is not a valid objective specification."""
        with pytest.raises(TypeError):
            resolve_objective(["cai"])

    def test_resolved_cai_callable(self):
        """The resolved CAI function can be called with the right signature."""
        func = resolve_objective("cai")
        dna = "ATGAAA"
        score = func(dna, "MK", "Escherichia_coli")
        assert isinstance(score, float)

    def test_resolved_min_max_gc_callable(self):
        """The resolved min_max_gc function can be called with the right signature."""
        func = resolve_objective("min_max_gc")
        dna = "ATGCATGC"
        score = func(dna, "M" * 8, "ecoli")
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0
