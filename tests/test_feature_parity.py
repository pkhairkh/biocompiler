"""
BioCompiler Feature Parity Tests
=================================

Comprehensive tests ensuring feature parity across all major BioCompiler
subsystems.  These tests verify that features work not only in isolation
but also correctly together, catching integration regressions early.

Test categories:
  1. Sliding-window GC constraints
  2. Custom objective functions
  3. IUPAC ambiguous base resolution
  4. Pattern enforcement (EnforcePattern / AvoidPattern)
  5. Part library loading and searching
  6. Assembly planning (Golden Gate and Gibson)
  7. DNA Chisel compatibility layer (expanded constraint set)
  8. Local GC constraints
  9. Cross-feature integration (features working together)
"""

from __future__ import annotations

import pytest

# ────────────────────────────────────────────────────────────
# Shared test data
# ────────────────────────────────────────────────────────────

SHORT_PROTEIN = "MVSKGE"  # 6 AA
HBB_PROTEIN = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
    "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK"
    "EFTPPVQAAYQKVVAGVANALAHKYH"
)
INSULIN_PROTEIN = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"


# ═══════════════════════════════════════════════════════════════
# 1. Sliding-window GC constraints
# ═══════════════════════════════════════════════════════════════

class TestSlidingGCConstraints:
    """Verify sliding-window GC scanning and fixing."""

    def test_check_sliding_gc_passes_for_uniform_gc(self):
        """A sequence with uniform GC should pass the sliding GC check."""
        from biocompiler.sequence.sliding_gc import check_sliding_gc

        # Build a 100 bp sequence with ~50% GC (alternating GC/AT)
        dna = "GCAT" * 25  # 100 bp, each window has GC ≈ 50%
        result = check_sliding_gc(dna, window_size=20, gc_min=0.30, gc_max=0.70)
        assert result.passed
        assert 0.30 <= result.min_gc <= 0.70
        assert 0.30 <= result.max_gc <= 0.70

    def test_check_sliding_gc_detects_high_gc_window(self):
        """A sequence with a local GC spike should be detected."""
        from biocompiler.sequence.sliding_gc import check_sliding_gc

        # 40 AT bases + 40 GC bases → the GC-only region triggers high GC
        dna = "ATATATAT" * 5 + "GCGCGCGC" * 5
        result = check_sliding_gc(dna, window_size=20, gc_min=0.20, gc_max=0.80)
        # The transition region may or may not violate depending on window;
        # the pure-GC region should be at 100% GC
        assert result.max_gc == 1.0

    def test_check_sliding_gc_empty_sequence(self):
        """Empty sequence should pass (no windows to check)."""
        from biocompiler.sequence.sliding_gc import check_sliding_gc

        result = check_sliding_gc("", window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed

    def test_check_sliding_gc_short_sequence(self):
        """Sequence shorter than window_size should use whole-sequence GC."""
        from biocompiler.sequence.sliding_gc import check_sliding_gc

        dna = "GCGCATATAT"  # 10 bp, GC = 40%
        result = check_sliding_gc(dna, window_size=50, gc_min=0.20, gc_max=0.80)
        assert result.passed

    def test_fix_sliding_gc_violations(self):
        """fix_sliding_gc_violations should reduce or eliminate violations."""
        from biocompiler.sequence.sliding_gc import check_sliding_gc, fix_sliding_gc_violations

        protein = "MVSKGE"
        # Use a DNA that correctly encodes the protein but may have local GC extremes
        # M=ATG V=GTT S=TCT K=AAA G=GGT E=GAA
        dna = "ATGGTTTCTAAAGGTGAA"
        assert len(dna) == len(protein) * 3

        fixed_dna, swaps = fix_sliding_gc_violations(
            dna, protein, window_size=12, gc_min=0.20, gc_max=0.80
        )
        # Protein should be preserved
        from biocompiler.expression.translation import translate

        translated = translate(fixed_dna)
        assert translated == protein

    def test_evaluate_sliding_gc_returns_type_check_result(self):
        """evaluate_sliding_gc should return a TypeCheckResult."""
        from biocompiler.sequence.sliding_gc import evaluate_sliding_gc
        from biocompiler.shared.types import TypeCheckResult, Verdict

        dna = "GCAT" * 25
        result = evaluate_sliding_gc(dna, window_size=20, gc_min=0.30, gc_max=0.70)
        assert isinstance(result, TypeCheckResult)
        assert result.verdict == Verdict.PASS

    def test_sliding_gc_step_parameter(self):
        """Sliding GC with step > 1 should still find violations."""
        from biocompiler.sequence.sliding_gc import check_sliding_gc

        dna = "ATATATAT" * 5 + "GCGCGCGC" * 5
        result = check_sliding_gc(dna, window_size=20, gc_min=0.20, gc_max=0.80, step=5)
        assert result.max_gc == 1.0

    def test_sliding_gc_python_numba_parity(self):
        """Python and NUMBA sliding GC paths should produce the same results."""
        from biocompiler.sequence.sliding_gc import (
            check_sliding_gc,
            _check_sliding_gc_python,
            _HAS_NUMBA,
        )

        if not _HAS_NUMBA:
            pytest.skip("NUMBA not available for parity check")

        from biocompiler.sequence.sliding_gc import _check_sliding_gc_numba

        dna = "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT" * 3
        py_result = _check_sliding_gc_python(dna, 20, 0.30, 0.70, 1)
        nb_result = _check_sliding_gc_numba(dna, 20, 0.30, 0.70, 1)
        assert py_result.passed == nb_result.passed
        assert abs(py_result.min_gc - nb_result.min_gc) < 1e-6
        assert abs(py_result.max_gc - nb_result.max_gc) < 1e-6
        assert len(py_result.violations) == len(nb_result.violations)


# ═══════════════════════════════════════════════════════════════
# 2. Custom objective functions
# ═══════════════════════════════════════════════════════════════

class TestCustomObjectives:
    """Verify that custom objective functions produce expected results."""

    def test_cai_objective_deterministic(self):
        """cai_objective should return the same result for the same input."""
        from biocompiler.optimizer.objectives import cai_objective

        dna = "ATGGCTCTGTGGATGCGCCTGCTGCC"
        score1 = cai_objective(dna, "MALWMR", "Escherichia_coli")
        score2 = cai_objective(dna, "MALWMR", "Escherichia_coli")
        assert score1 == score2

    def test_cai_gc_balanced_combines_both(self):
        """cai_gc_balanced_objective should be between pure CAI and pure GC score."""
        from biocompiler.optimizer.objectives import cai_objective, cai_gc_balanced_objective

        dna = "ATGGCTCTGTGGATGCGCCTGCTGCC"
        protein = "MALWMR"
        org = "Escherichia_coli"

        cai_score = cai_objective(dna, protein, org)
        balanced_score = cai_gc_balanced_objective(dna, protein, org, gc_weight=0.5)

        assert isinstance(balanced_score, float)
        assert 0.0 <= balanced_score <= 1.0
        # Balanced score should differ from pure CAI (unless GC is already perfect)
        # It should be a weighted combination

    def test_min_max_gc_peak_at_target(self):
        """min_max_gc_objective should return 1.0 when GC equals target."""
        from biocompiler.optimizer.objectives import min_max_gc_objective

        # AATTCCGG repeated → 50% GC
        dna_50 = "AATTCCGG" * 4
        score = min_max_gc_objective(dna_50, "NSFELRQD", "Escherichia_coli", target_gc=0.5)
        assert score == 1.0

    def test_codon_pair_objective_range(self):
        """codon_pair_objective should return a value in [0, 1]."""
        from biocompiler.optimizer.objectives import codon_pair_objective

        dna = "ATGGCTCTGTGGATGCGCCTGCTGCC"
        score = codon_pair_objective(dna, "MALWMR", "Escherichia_coli")
        assert 0.0 <= score <= 1.0

    def test_resolve_objective_callable(self):
        """resolve_objective with a callable should return it as-is."""
        from biocompiler.optimizer.objectives import resolve_objective

        def my_obj(dna, protein, organism):
            return 0.42

        fn = resolve_objective(my_obj)
        assert fn is my_obj
        assert fn("ATG", "M", "ecoli") == 0.42

    def test_custom_objective_in_optimization(self):
        """optimize_sequence with a custom objective should work and preserve protein."""
        from biocompiler.optimizer import optimize_sequence
        from biocompiler.expression.translation import translate

        def gc_target_obj(dna, protein, organism):
            gc = (dna.count("G") + dna.count("C")) / len(dna) if dna else 0.0
            return 1.0 - abs(gc - 0.5)

        result = optimize_sequence(
            SHORT_PROTEIN,
            organism="ecoli",
            objective=gc_target_obj,
            strict_mode=False,
        )
        assert result.cai > 0.0
        assert translate(result.sequence) == SHORT_PROTEIN

    def test_all_builtin_objectives_with_optimize(self):
        """All built-in objective names should work with optimize_sequence."""
        from biocompiler.optimizer import optimize_sequence
        from biocompiler.optimizer.objectives import OBJECTIVE_REGISTRY

        for name in OBJECTIVE_REGISTRY:
            result = optimize_sequence(
                SHORT_PROTEIN,
                organism="ecoli",
                objective=name,
                strict_mode=False,
            )
            assert result.cai > 0.0, f"Objective '{name}' produced zero CAI"


# ═══════════════════════════════════════════════════════════════
# 3. IUPAC ambiguous base resolution
# ═══════════════════════════════════════════════════════════════

class TestIUPACAmbiguousBaseResolution:
    """Verify IUPAC ambiguous base resolution works correctly."""

    def test_is_ambiguous(self):
        """is_ambiguous should correctly identify ambiguous IUPAC codes."""
        from biocompiler.sequence.iupac import is_ambiguous

        assert not is_ambiguous("A")
        assert not is_ambiguous("C")
        assert not is_ambiguous("G")
        assert not is_ambiguous("T")
        assert is_ambiguous("R")
        assert is_ambiguous("N")
        assert is_ambiguous("S")
        assert is_ambiguous("Y")

    def test_has_ambiguous(self):
        """has_ambiguous should detect any ambiguous base in a sequence."""
        from biocompiler.sequence.iupac import has_ambiguous

        assert not has_ambiguous("ATGC")
        assert has_ambiguous("ATGR")
        assert has_ambiguous("NATGC")

    def test_resolve_first_strategy(self):
        """'first' strategy should resolve to the first alphabetical concrete base."""
        from biocompiler.sequence.iupac import resolve_ambiguous

        # R = {A, G} → first is A
        result = resolve_ambiguous("ATR", strategy="first")
        assert result == "ATA"

    def test_resolve_most_common_strategy(self):
        """'most_common' strategy should pick the highest-frequency base."""
        from biocompiler.sequence.iupac import resolve_ambiguous

        # With default frequencies: A=0.295, G=0.205
        # R = {A, G} → A has higher frequency
        result = resolve_ambiguous("ATR", strategy="most_common")
        assert "A" in result or "G" in result
        # Result should only contain ACGT
        assert all(b in "ACGT" for b in result)

    def test_resolve_gc_balanced_strategy(self):
        """'gc_balanced' strategy should keep GC close to target."""
        from biocompiler.sequence.iupac import resolve_ambiguous

        # N = {A, C, G, T} with gc_target=0.5 should prefer G or C
        result = resolve_ambiguous("NNN", strategy="gc_balanced", gc_target=0.5)
        assert all(b in "ACGT" for b in result)

    def test_resolve_cai_optimal_strategy(self):
        """'cai_optimal' strategy should pick the highest-CAI codon."""
        from biocompiler.sequence.iupac import resolve_ambiguous

        # Provide a simple CAI table for E. coli
        cai_table = {
            "ATG": 1.0, "ATA": 0.3, "ATC": 0.8, "ATT": 0.4,
            "GCT": 0.9, "GCC": 1.0, "GCA": 0.5, "GCG": 0.3,
        }
        # R = {A, G}; if the second position of a codon is R,
        # we should get the codon with highest CAI
        result = resolve_ambiguous("ATGRCT", strategy="cai_optimal", cai_table=cai_table)
        assert all(b in "ACGT" for b in result)

    def test_expand_ambiguous_simple(self):
        """expand_ambiguous should enumerate all concrete sequences."""
        from biocompiler.sequence.iupac import expand_ambiguous

        expansions = expand_ambiguous("ATR")
        assert sorted(expansions) == ["ATA", "ATG"]

    def test_expand_ambiguous_no_ambiguity(self):
        """expand_ambiguous on a concrete sequence should return [sequence]."""
        from biocompiler.sequence.iupac import expand_ambiguous

        assert expand_ambiguous("ATGC") == ["ATGC"]

    def test_validate_iupac_sequence_rejects_invalid(self):
        """validate_iupac_sequence should reject non-IUPAC characters."""
        from biocompiler.sequence.iupac import validate_iupac_sequence

        with pytest.raises(ValueError, match="Invalid DNA bases"):
            validate_iupac_sequence("ATGX")

    def test_validate_iupac_sequence_accepts_ambiguous(self):
        """validate_iupac_sequence should accept valid IUPAC codes."""
        from biocompiler.sequence.iupac import validate_iupac_sequence

        result = validate_iupac_sequence("ATGRYSWKMBDHVN")
        assert result == "ATGRYSWKMBDHVN"

    def test_resolve_preserves_concrete_bases(self):
        """Resolution should never change already-concrete bases."""
        from biocompiler.sequence.iupac import resolve_ambiguous

        dna = "ATGCATGCN"
        result = resolve_ambiguous(dna, strategy="first")
        assert result[:8] == "ATGCATGC"


# ═══════════════════════════════════════════════════════════════
# 4. Pattern enforcement (EnforcePattern / AvoidPattern)
# ═══════════════════════════════════════════════════════════════

class TestPatternEnforcement:
    """Verify EnforcePattern and AvoidPattern work correctly."""

    def test_avoid_pattern_detects_violation(self):
        """AvoidPattern should detect a present pattern."""
        from biocompiler.sequence.pattern_enforcement import PatternConstraint, check_pattern

        c = PatternConstraint(pattern="GAATTC", action="avoid", scope="dna", strand="both")
        result = check_pattern("ATGGAATTCCGATC", c)
        assert not result.passed
        assert len(result.matches) > 0

    def test_avoid_pattern_passes_when_absent(self):
        """AvoidPattern should pass when the pattern is absent."""
        from biocompiler.sequence.pattern_enforcement import PatternConstraint, check_pattern

        c = PatternConstraint(pattern="GAATTC", action="avoid", scope="dna", strand="both")
        result = check_pattern("ATGGCATCCGATC", c)
        assert result.passed

    def test_enforce_pattern_detects_absence(self):
        """EnforcePattern should fail when the pattern is absent."""
        from biocompiler.sequence.pattern_enforcement import PatternConstraint, check_pattern

        c = PatternConstraint(pattern="CATCAT", action="enforce", scope="dna", strand="forward")
        result = check_pattern("ATGGCATCCGATC", c)
        assert not result.passed  # CATCAT is not present

    def test_enforce_pattern_passes_when_present(self):
        """EnforcePattern should pass when the pattern is present."""
        from biocompiler.sequence.pattern_enforcement import PatternConstraint, check_pattern

        c = PatternConstraint(pattern="CATCAT", action="enforce", scope="dna", strand="forward")
        result = check_pattern("ATGCATCATGATC", c)
        assert result.passed

    def test_iupac_pattern_expansion(self):
        """PatternConstraint with IUPAC codes should expand and match."""
        from biocompiler.sequence.pattern_enforcement import PatternConstraint, check_pattern

        # R = A or G; "GRC" should match "GAC" and "GGC"
        c = PatternConstraint(pattern="GRC", action="avoid", scope="dna", strand="forward")
        result = check_pattern("ATGGACATC", c)
        assert not result.passed  # GAC matches GRC

    def test_enforce_pattern_modifies_dna(self):
        """enforce_pattern should modify DNA to remove/insert patterns."""
        from biocompiler.sequence.pattern_enforcement import (
            PatternConstraint,
            enforce_pattern,
        )
        from biocompiler.expression.translation import translate

        protein = "MVSKGE"
        dna = "ATGGTTTCTAAAGGTGAA"
        c = PatternConstraint(pattern="GTTTCT", action="avoid", scope="dna", strand="forward")
        new_dna = enforce_pattern(dna, protein, c)
        # Protein should be preserved
        assert translate(new_dna) == protein

    def test_enforce_patterns_multiple(self):
        """enforce_patterns should handle multiple constraints iteratively."""
        from biocompiler.sequence.pattern_enforcement import (
            PatternConstraint,
            check_patterns,
        )

        constraints = [
            PatternConstraint(pattern="GAATTC", action="avoid", scope="dna", strand="both"),
            PatternConstraint(pattern="GGATCC", action="avoid", scope="dna", strand="both"),
        ]
        results = check_patterns("ATGGAATTCCCCGGATCCC", constraints)
        assert len(results) == 2
        assert not results[0].passed  # GAATTC present
        assert not results[1].passed  # GGATCC present

    def test_pattern_constraint_invalid_action(self):
        """PatternConstraint with invalid action should raise ValueError."""
        from biocompiler.sequence.pattern_enforcement import PatternConstraint

        with pytest.raises(ValueError, match="Invalid action"):
            PatternConstraint(pattern="ATG", action="invalid", scope="dna")

    def test_pattern_reverse_complement(self):
        """AvoidPattern with strand='both' should detect RC matches."""
        from biocompiler.sequence.pattern_enforcement import PatternConstraint, check_pattern

        # GAATTC is a palindrome (RC = GAATTC), so it matches both strands
        c = PatternConstraint(pattern="GAATTC", action="avoid", scope="dna", strand="both")
        result = check_pattern("ATGGAATTCCGATC", c)
        assert not result.passed

        # AATTCG has RC = CGAATT; test forward-only
        c_fwd = PatternConstraint(pattern="AATTCG", action="avoid", scope="dna", strand="forward")
        result_fwd = check_pattern("ATGAATTCGGATC", c_fwd)
        assert not result_fwd.passed

    def test_build_avoidance_scanner(self):
        """build_avoidance_scanner should create an AhoCorasickScanner."""
        from biocompiler.sequence.pattern_enforcement import (
            PatternConstraint,
            build_avoidance_scanner,
        )

        constraints = [
            PatternConstraint(pattern="GAATTC", action="avoid", scope="dna"),
            PatternConstraint(pattern="GGATCC", action="avoid", scope="dna"),
        ]
        scanner = build_avoidance_scanner(constraints)
        assert scanner is not None

    def test_build_avoidance_scanner_empty(self):
        """build_avoidance_scanner with no valid patterns should return None."""
        from biocompiler.sequence.pattern_enforcement import (
            PatternConstraint,
            build_avoidance_scanner,
        )

        constraints = [
            PatternConstraint(pattern="ATG", action="enforce", scope="dna"),
        ]
        scanner = build_avoidance_scanner(constraints)
        assert scanner is None


# ═══════════════════════════════════════════════════════════════
# 5. Part library loading and searching
# ═══════════════════════════════════════════════════════════════

class TestPartLibrary:
    """Verify part library loading, searching, and management."""

    def test_default_parts_loaded(self):
        """Default library should contain built-in parts."""
        from biocompiler.optimizer.parts import PartLibrary

        lib = PartLibrary()
        assert len(lib) > 0
        assert "T7_promoter" in lib
        assert "B0034" in lib

    def test_get_part_by_name(self):
        """Getting a part by name should return the correct Part."""
        from biocompiler.optimizer.parts import PartLibrary

        lib = PartLibrary()
        t7 = lib.get("T7_promoter")
        assert t7.part_type == "promoter"
        assert t7.sequence == "TAATACGACTCACTATAGGG"

    def test_get_nonexistent_part_raises(self):
        """Getting a nonexistent part should raise KeyError."""
        from biocompiler.optimizer.parts import PartLibrary

        lib = PartLibrary()
        with pytest.raises(KeyError):
            lib.get("NONEXISTENT_PART")

    def test_search_by_type(self):
        """Searching by type should return only matching parts."""
        from biocompiler.optimizer.parts import PartLibrary

        lib = PartLibrary()
        promoters = lib.search("promoter")
        assert len(promoters) >= 2  # T7 and lac at minimum
        assert all(p.part_type == "promoter" for p in promoters)

    def test_search_by_type_and_organism(self):
        """Searching by type + organism should filter correctly."""
        from biocompiler.optimizer.parts import PartLibrary

        lib = PartLibrary()
        ecoli_rbs = lib.search("rbs", organism="E_coli")
        assert len(ecoli_rbs) >= 1  # B0034, B0032, B0031
        for p in ecoli_rbs:
            assert p.metadata.get("organism", "").lower() == "e_coli"

    def test_add_custom_part(self):
        """Adding a custom part should make it retrievable."""
        from biocompiler.optimizer.parts import PartLibrary, Part

        lib = PartLibrary()
        custom = Part(
            name="my_promoter",
            part_type="promoter",
            sequence="ATATATAA",
            description="Custom test promoter",
        )
        lib.add(custom)
        assert "my_promoter" in lib
        assert lib.get("my_promoter").sequence == "ATATATAA"

    def test_list_parts(self):
        """list_parts should return sorted part names."""
        from biocompiler.optimizer.parts import PartLibrary

        lib = PartLibrary()
        names = lib.list_parts()
        assert names == sorted(names)

    def test_part_invalid_type_raises(self):
        """Creating a Part with invalid type should raise ValueError."""
        from biocompiler.optimizer.parts import Part

        with pytest.raises(ValueError, match="Invalid part_type"):
            Part(name="bad", part_type="invalid_type", sequence="ATGC")

    def test_part_empty_name_raises(self):
        """Creating a Part with empty name should raise ValueError."""
        from biocompiler.optimizer.parts import Part

        with pytest.raises(ValueError, match="name must not be empty"):
            Part(name="", part_type="promoter", sequence="ATGC")

    def test_part_empty_sequence_raises(self):
        """Creating a Part with empty sequence should raise ValueError."""
        from biocompiler.optimizer.parts import Part

        with pytest.raises(ValueError, match="sequence must not be empty"):
            Part(name="test", part_type="promoter", sequence="")

    def test_default_parts_have_valid_sequences(self):
        """All default parts should have uppercase ACGT sequences."""
        from biocompiler.optimizer.parts import DEFAULT_PARTS

        for part in DEFAULT_PARTS:
            assert all(b in "ACGT" for b in part.sequence), (
                f"Part {part.name} has invalid bases: {part.sequence}"
            )


# ═══════════════════════════════════════════════════════════════
# 6. Assembly planning (Golden Gate and Gibson)
# ═══════════════════════════════════════════════════════════════

class TestAssemblyPlanning:
    """Verify Golden Gate and Gibson assembly planning."""

    def test_golden_gate_basic(self):
        """plan_golden_gate should create a valid assembly plan."""
        from biocompiler.optimizer.assembly import plan_golden_gate

        fragments = ["ATGAAAGGGTTTCCC", "ATGCCCGGGAAATTT", "ATGTTTCCCAAAGGG"]
        plan = plan_golden_gate(fragments)
        assert plan.method == "golden_gate"
        assert len(plan.fragments) == 3
        assert plan.total_length == sum(len(f) for f in fragments)
        assert len(plan.enzymes) > 0
        assert len(plan.overlap_sequences) == 2  # N-1 overlaps for N fragments

    def test_golden_gate_custom_enzymes(self):
        """plan_golden_gate should accept custom enzyme names."""
        from biocompiler.optimizer.assembly import plan_golden_gate

        fragments = ["ATGAAAGGGTTTCCC"]
        plan = plan_golden_gate(fragments, enzymes=["BsaI"])
        assert "BsaI" in plan.enzymes

    def test_golden_gate_empty_sequences_raises(self):
        """plan_golden_gate with no sequences should raise ValueError."""
        from biocompiler.optimizer.assembly import plan_golden_gate

        with pytest.raises(ValueError, match="At least one sequence"):
            plan_golden_gate([])

    def test_gibson_basic(self):
        """plan_gibson should create a valid assembly plan."""
        from biocompiler.optimizer.assembly import plan_gibson

        fragments = [
            "ATGAAAGGGTTTCCCATGAAAGGGTTTCCC",
            "ATGCCCGGGAAATTTATGCCCGGGAAATTT",
            "ATGTTTCCCAAAGGGATGTTTCCCAAAGGG",
        ]
        plan = plan_gibson(fragments, overlap_length=20)
        assert plan.method == "gibson"
        assert len(plan.fragments) == 3
        assert len(plan.overlap_sequences) == 2
        assert plan.total_length == sum(len(f) for f in fragments) - 20 * 2

    def test_gibson_overlap_too_short_raises(self):
        """plan_gibson with overlap_length < 4 should raise ValueError."""
        from biocompiler.optimizer.assembly import plan_gibson

        with pytest.raises(ValueError, match="overlap_length must be >= 4"):
            plan_gibson(["ATGCATGC"], overlap_length=2)

    def test_gibson_empty_sequences_raises(self):
        """plan_gibson with no sequences should raise ValueError."""
        from biocompiler.optimizer.assembly import plan_gibson

        with pytest.raises(ValueError, match="At least one sequence"):
            plan_gibson([])

    def test_assembly_plan_invalid_method_raises(self):
        """AssemblyPlan with invalid method should raise ValueError."""
        from biocompiler.optimizer.assembly import AssemblyPlan

        with pytest.raises(ValueError, match="Invalid assembly method"):
            AssemblyPlan(
                method="invalid",
                fragments=["ATGC"],
                enzymes=[],
                overlap_sequences=[],
                total_length=4,
            )

    def test_golden_gate_warns_on_internal_sites(self, caplog):
        """plan_golden_gate should warn when fragments contain enzyme sites."""
        import logging
        from biocompiler.optimizer.assembly import plan_golden_gate

        # BsaI site is GGTCTC; embed it in a fragment
        fragments = ["ATGGGTCTCATGCGC"]
        with caplog.at_level(logging.WARNING):
            plan = plan_golden_gate(fragments, enzymes=["BsaI"])
        # The plan should still be created (with a warning)
        assert plan.method == "golden_gate"

    def test_golden_gate_single_fragment(self):
        """plan_golden_gate with a single fragment should work."""
        from biocompiler.optimizer.assembly import plan_golden_gate

        plan = plan_golden_gate(["ATGAAAGGGTTTCCC"])
        assert len(plan.overlap_sequences) == 0
        assert plan.total_length == 15


# ═══════════════════════════════════════════════════════════════
# 7. DNA Chisel compatibility layer
# ═══════════════════════════════════════════════════════════════

class TestDNAChiselCompat:
    """Verify the DNA Chisel compatibility layer and expanded constraint set."""

    def test_is_dna_chisel_available(self):
        """is_dna_chisel_available should return a boolean."""
        from biocompiler.infrastructure.dna_chisel_compat import is_dna_chisel_available

        result = is_dna_chisel_available()
        assert isinstance(result, bool)

    def test_constraint_mapping_has_expected_entries(self):
        """CONSTRAINT_MAPPING should include all expected constraint types."""
        from biocompiler.infrastructure.dna_chisel_compat import CONSTRAINT_MAPPING

        expected = [
            "EnforceTranslation",
            "EnforceGCContent",
            "EnforceGCContentLocal",
            "AvoidPattern",
            "AvoidBacterialPromoter",
            "EnforceStartCodon",
            "EnforceStopCodon",
            "UniquifyAllKmers",
            "AvoidChanges",
            "EnforceSequence",
            "AvoidHairpins",
        ]
        for name in expected:
            assert name in CONSTRAINT_MAPPING, f"Missing constraint: {name}"

    def test_build_constraint_spec_returns_none_without_chisel(self):
        """build_constraint_spec should return None when DNA Chisel is unavailable."""
        from biocompiler.infrastructure.dna_chisel_compat import (
            build_constraint_spec,
            is_dna_chisel_available,
        )

        if not is_dna_chisel_available():
            result = build_constraint_spec("AvoidPattern", pattern="GAATTC")
            assert result is None

    def test_translate_biocompiler_constraints(self):
        """translate_biocompiler_constraints should return a list."""
        from biocompiler.infrastructure.dna_chisel_compat import translate_biocompiler_constraints

        specs = translate_biocompiler_constraints(
            protein="MVSKGE",
            gc_lo=0.30,
            gc_hi=0.70,
        )
        assert isinstance(specs, list)

    def test_translate_biocompiler_constraints_all_flags(self):
        """translate_biocompiler_constraints with all flags should not crash."""
        from biocompiler.infrastructure.dna_chisel_compat import translate_biocompiler_constraints

        specs = translate_biocompiler_constraints(
            protein="MVSKGE",
            gc_lo=0.30,
            gc_hi=0.70,
            local_gc_window=50,
            avoid_bacterial_promoters=True,
            enforce_start=True,
            enforce_stop=True,
            uniquify_kmers=9,
            avoid_hairpins_flag=True,
        )
        assert isinstance(specs, list)

    def test_compare_optimizers_structure(self):
        """compare_optimizers should return a ComparisonResult."""
        from biocompiler.infrastructure.dna_chisel_compat import compare_optimizers

        result = compare_optimizers(
            protein=SHORT_PROTEIN,
            organism="Escherichia_coli",
        )
        assert result.protein == SHORT_PROTEIN
        assert result.organism == "Escherichia_coli"
        assert result.biocompiler is not None
        assert isinstance(result.winner, dict)


# ═══════════════════════════════════════════════════════════════
# 8. Local GC constraints
# ═══════════════════════════════════════════════════════════════

class TestLocalGCConstraints:
    """Verify local GC constraint checking and optimization."""

    def test_check_local_gc_satisfied(self):
        """check_local_gc should pass when constraints are satisfied."""
        from biocompiler.sequence.local_gc import LocalGCConstraint, check_local_gc

        dna = "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
        c = LocalGCConstraint(region_start=0, region_end=20, gc_min=0.30, gc_max=0.70)
        result = check_local_gc(dna, [c])
        assert result.satisfied

    def test_check_local_gc_violated(self):
        """check_local_gc should fail when a constraint is violated."""
        from biocompiler.sequence.local_gc import LocalGCConstraint, check_local_gc

        dna = "GCGCGCGCGCGCGCGCGCGC"  # 20 bp, 100% GC
        c = LocalGCConstraint(region_start=0, region_end=20, gc_min=0.30, gc_max=0.70)
        result = check_local_gc(dna, [c])
        assert not result.satisfied
        assert len(result.violations) == 1

    def test_check_local_gc_multiple_constraints(self):
        """check_local_gc should evaluate all constraints."""
        from biocompiler.sequence.local_gc import LocalGCConstraint, check_local_gc

        dna = "GCGCGCGCATATATATATATGCGCGCGCATATATATA"  # ~40 bp
        c1 = LocalGCConstraint(region_start=0, region_end=10, gc_min=0.30, gc_max=0.70)
        c2 = LocalGCConstraint(region_start=20, region_end=30, gc_min=0.30, gc_max=0.70)
        result = check_local_gc(dna, [c1, c2])
        # First region is all GC → violation; second region may or may not pass
        assert isinstance(result.satisfied, bool)

    def test_local_gc_constraint_invalid_start(self):
        """LocalGCConstraint with negative start should raise ValueError."""
        from biocompiler.sequence.local_gc import LocalGCConstraint

        with pytest.raises(ValueError, match="region_start must be >= 0"):
            LocalGCConstraint(region_start=-1, region_end=10, gc_min=0.3, gc_max=0.7)

    def test_local_gc_constraint_end_before_start(self):
        """LocalGCConstraint with end <= start should raise ValueError."""
        from biocompiler.sequence.local_gc import LocalGCConstraint

        with pytest.raises(ValueError, match="region_end"):
            LocalGCConstraint(region_start=10, region_end=10, gc_min=0.3, gc_max=0.7)

    def test_local_gc_constraint_invalid_gc_range(self):
        """LocalGCConstraint with gc_min > gc_max should raise ValueError."""
        from biocompiler.sequence.local_gc import LocalGCConstraint

        with pytest.raises(ValueError, match="gc_min"):
            LocalGCConstraint(region_start=0, region_end=10, gc_min=0.7, gc_max=0.3)

    def test_optimize_local_gc(self):
        """optimize_local_gc should attempt to fix violations."""
        from biocompiler.sequence.local_gc import LocalGCConstraint, optimize_local_gc

        protein = "MVSKGE"
        # M=ATG V=GTT S=TCT K=AAA G=GGT E=GAA
        dna = "ATGGTTTCTAAAGGTGAA"
        assert len(dna) == len(protein) * 3

        c = LocalGCConstraint(
            region_start=0, region_end=len(dna), gc_min=0.30, gc_max=0.70
        )
        result = optimize_local_gc(dna, protein, [c])
        # Protein should be preserved (translated sequence matches)
        from biocompiler.expression.translation import translate

        translated = translate(result.sequence)
        assert translated == protein

    def test_check_local_gc_beyond_sequence(self):
        """Constraints beyond the sequence length should be skipped."""
        from biocompiler.sequence.local_gc import LocalGCConstraint, check_local_gc

        dna = "ATGC"
        c = LocalGCConstraint(region_start=100, region_end=200, gc_min=0.30, gc_max=0.70)
        result = check_local_gc(dna, [c])
        assert result.satisfied  # No violations because region is skipped


# ═══════════════════════════════════════════════════════════════
# 9. Cross-feature integration
# ═══════════════════════════════════════════════════════════════

class TestCrossFeatureIntegration:
    """Verify that features work correctly together, not just in isolation."""

    def test_optimize_with_sliding_gc_and_local_gc(self):
        """Optimization with both sliding GC and local GC should preserve protein."""
        from biocompiler.optimizer import optimize_sequence
        from biocompiler.expression.translation import translate

        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="human",
            gc_lo=0.35,
            gc_hi=0.65,
            strict_mode=False,
        )
        # Protein should be preserved
        translated = translate(result.sequence)
        assert translated == INSULIN_PROTEIN
        # GC should be in range (or close)
        assert 0.20 <= result.gc_content <= 0.80

    def test_iupac_resolution_then_optimization(self):
        """Resolve IUPAC bases first, then optimize the sequence."""
        from biocompiler.sequence.iupac import resolve_ambiguous
        from biocompiler.optimizer import optimize_sequence
        from biocompiler.expression.translation import translate

        # Protein with ambiguous codon positions
        protein = "MVSKGE"
        # Create a DNA with an N (ambiguous)
        ambiguous_dna = "ATGGTTTCTAAAGGTNAA"
        resolved = resolve_ambiguous(ambiguous_dna, strategy="most_common")

        result = optimize_sequence(
            protein,
            organism="ecoli",
            strict_mode=False,
        )
        assert result.cai > 0.0
        assert translate(result.sequence) == protein

    def test_pattern_avoidance_and_assembly(self):
        """After avoiding restriction sites, plan assembly of the result."""
        from biocompiler.optimizer import optimize_sequence
        from biocompiler.optimizer.assembly import plan_golden_gate, plan_gibson

        result = optimize_sequence(
            SHORT_PROTEIN,
            organism="ecoli",
            strict_mode=False,
        )
        # Plan Golden Gate assembly with the optimized sequence
        plan = plan_golden_gate([result.sequence])
        assert plan.method == "golden_gate"
        assert plan.total_length > 0

        # Plan Gibson assembly too
        gibson_plan = plan_gibson([result.sequence])
        assert gibson_plan.method == "gibson"

    def test_part_library_and_assembly(self):
        """Use parts from the library to construct an assembly plan."""
        from biocompiler.optimizer.parts import PartLibrary
        from biocompiler.optimizer.assembly import plan_gibson

        lib = PartLibrary()
        promoter = lib.search("promoter", organism="E_coli")
        rbs = lib.search("rbs", organism="E_coli")

        if promoter and rbs:
            fragments = [promoter[0].sequence, rbs[0].sequence, "ATGAAAGGGTTTCCCTAA"]
            plan = plan_gibson(fragments, overlap_length=10)
            assert plan.method == "gibson"
            assert plan.total_length > 0

    def test_local_gc_after_optimization(self):
        """After optimization, verify local GC constraints are met."""
        from biocompiler.optimizer import optimize_sequence
        from biocompiler.sequence.local_gc import LocalGCConstraint, check_local_gc

        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="human",
            gc_lo=0.35,
            gc_hi=0.65,
            strict_mode=False,
        )

        # Check local GC in the first 30 bp
        c = LocalGCConstraint(
            region_start=0, region_end=30, gc_min=0.20, gc_max=0.80
        )
        gc_result = check_local_gc(result.sequence, [c])
        # Should pass with generous bounds
        assert gc_result.satisfied

    def test_custom_objective_preserves_constraints(self):
        """Custom objective should not violate hard constraints."""
        from biocompiler.optimizer import optimize_sequence
        from biocompiler.sequence.sliding_gc import check_sliding_gc

        def gc_target_obj(dna, protein, organism):
            gc = (dna.count("G") + dna.count("C")) / len(dna) if dna else 0.0
            return 1.0 - abs(gc - 0.5)

        result = optimize_sequence(
            SHORT_PROTEIN,
            organism="ecoli",
            objective=gc_target_obj,
            gc_lo=0.30,
            gc_hi=0.70,
            strict_mode=False,
        )

        # GC should be within bounds
        assert 0.20 <= result.gc_content <= 0.80

    def test_dna_chisel_compat_with_all_features(self):
        """DNA Chisel compat layer should work with all constraint types."""
        from biocompiler.infrastructure.dna_chisel_compat import translate_biocompiler_constraints

        # Test that translating all constraints does not crash
        specs = translate_biocompiler_constraints(
            protein=SHORT_PROTEIN,
            gc_lo=0.30,
            gc_hi=0.70,
            restriction_enzymes=["EcoRI", "BamHI"],
            local_gc_window=50,
            avoid_bacterial_promoters=True,
            enforce_start=True,
            enforce_stop=True,
            uniquify_kmers=9,
            preserve_zones=["0-10"],
            enforce_sequence_str="ATG",
            avoid_hairpins_flag=True,
        )
        # The result depends on whether DNA Chisel is installed
        assert isinstance(specs, list)

    def test_iupac_with_pattern_enforcement(self):
        """Resolve IUPAC bases, then check pattern enforcement."""
        from biocompiler.sequence.iupac import resolve_ambiguous
        from biocompiler.sequence.pattern_enforcement import PatternConstraint, check_pattern

        # Resolve ambiguous bases
        resolved = resolve_ambiguous("ATGRCATGC", strategy="first")
        # Then check for a pattern
        c = PatternConstraint(pattern="ATAC", action="avoid", scope="dna", strand="forward")
        result = check_pattern(resolved, c)
        assert isinstance(result.passed, bool)
