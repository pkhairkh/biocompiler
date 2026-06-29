"""
Property-based tests for the type_system module using Hypothesis.

Covers four core properties:
  1. evaluate_all_predicates always returns a list of TypeCheckResult
  2. Verdict enum values are always one of PASS/FAIL/UNCERTAIN/LIKELY_PASS/LIKELY_FAIL
  3. For any valid DNA sequence, type checking does not raise exceptions
  4. combined_verdict returns a valid Verdict
"""

import pytest
pytest.importorskip("hypothesis")
from hypothesis import given, strategies as st, assume, settings

from biocompiler.shared.types import Verdict, TypeCheckResult, combined_verdict
from biocompiler.type_system import (
    CODON_TABLE,
    AA_TO_CODONS,
    check_no_stop_codons,
    check_no_cpg_island,
    check_no_gt_dinucleotide,
    check_no_avoidable_gt,
    check_valid_coding_seq,
    check_no_cryptic_promoter,
    check_no_unexpected_tm_domain,
    check_mrna_secondary_structure,
    evaluate_all_predicates,
)


# ────────────────────────────────────────────────────────────
# Strategies
# ────────────────────────────────────────────────────────────

# Single DNA base
dna_base = st.sampled_from("ACGT")

# Short DNA sequence (any length, may not be valid coding)
dna_seq = st.text(alphabet="ACGT", min_size=0, max_size=300)

# Valid coding DNA: length divisible by 3, all codons from CODON_TABLE
valid_codons = st.sampled_from(list(CODON_TABLE.keys()))

valid_coding_seq = st.integers(min_value=1, max_value=50).flatmap(
    lambda n: st.lists(valid_codons, min_size=n, max_size=n).map(
        lambda codons: "".join(codons)
    )
)

# Strategy for a non-empty valid coding sequence (at least 1 codon)
valid_coding_seq_nonempty = st.integers(min_value=1, max_value=50).flatmap(
    lambda n: st.lists(valid_codons, min_size=n, max_size=n).map(
        lambda codons: "".join(codons)
    )
)

# Verdict strategy
verdict = st.sampled_from(list(Verdict))

# List of verdicts
verdict_lists = st.lists(verdict, min_size=0, max_size=30)


# ────────────────────────────────────────────────────────────
# Property 1: evaluate_all_predicates always returns a list of TypeCheckResult
# ────────────────────────────────────────────────────────────

class TestEvaluateAllPredicatesReturnType:
    """Property: evaluate_all_predicates always returns a list of TypeCheckResult."""

    @given(seq=valid_coding_seq)
    @settings(max_examples=30, deadline=5000)
    def test_returns_list_of_type_check_results(self, seq):
        """evaluate_all_predicates returns a list where every element is a TypeCheckResult."""
        results = evaluate_all_predicates(seq)
        assert isinstance(results, list), f"Expected list, got {type(results)}"
        for r in results:
            assert isinstance(r, TypeCheckResult), (
                f"Expected TypeCheckResult, got {type(r)}"
            )

    @given(seq=valid_coding_seq)
    @settings(max_examples=30, deadline=5000)
    def test_results_have_verdict_attribute(self, seq):
        """Every TypeCheckResult has a verdict attribute that is a Verdict."""
        results = evaluate_all_predicates(seq)
        for r in results:
            assert isinstance(r.verdict, Verdict), (
                f"Expected Verdict, got {type(r.verdict)} for predicate {r.predicate}"
            )

    @given(seq=valid_coding_seq)
    @settings(max_examples=30, deadline=5000)
    def test_results_have_predicate_name(self, seq):
        """Every TypeCheckResult has a non-empty predicate name."""
        results = evaluate_all_predicates(seq)
        for r in results:
            assert isinstance(r.predicate, str) and len(r.predicate) > 0, (
                f"Expected non-empty predicate string, got {r.predicate!r}"
            )

    @given(seq=valid_coding_seq)
    @settings(max_examples=20, deadline=5000)
    def test_returns_12_results(self, seq):
        """evaluate_all_predicates returns exactly 12 TypeCheckResults."""
        results = evaluate_all_predicates(seq)
        assert len(results) == 12, (
            f"Expected 12 results, got {len(results)}"
        )


# ────────────────────────────────────────────────────────────
# Property 2: Verdict enum values are always one of the five expected
# ────────────────────────────────────────────────────────────

class TestVerdictEnumValues:
    """Property: Verdict enum values are always one of
    PASS/FAIL/UNCERTAIN/LIKELY_PASS/LIKELY_FAIL."""

    VALID_VERDICT_NAMES = {"PASS", "FAIL", "UNCERTAIN", "LIKELY_PASS", "LIKELY_FAIL"}

    @given(v=verdict)
    def test_verdict_name_is_valid(self, v):
        """Every Verdict value's name is in the expected set."""
        assert v.name in self.VALID_VERDICT_NAMES, (
            f"Verdict {v.name} not in {self.VALID_VERDICT_NAMES}"
        )

    def test_verdict_enum_has_exactly_five_members(self):
        """The Verdict enum has exactly 5 members."""
        assert len(Verdict) == 5, f"Expected 5 Verdict members, got {len(Verdict)}"

    def test_all_five_verdicts_exist(self):
        """All five expected Verdict values exist in the enum."""
        for name in self.VALID_VERDICT_NAMES:
            assert hasattr(Verdict, name), f"Verdict.{name} is missing"

    @given(v=verdict)
    def test_verdict_is_instance_of_enum(self, v):
        """Every Verdict value is indeed an instance of the Verdict enum."""
        assert isinstance(v, Verdict)

    @given(v=verdict)
    def test_verdict_string_value_matches_name(self, v):
        """Since Verdict(str, Enum), the value equals the name."""
        assert v.value == v.name


# ────────────────────────────────────────────────────────────
# Property 3: For any valid DNA sequence, type checking does not raise exceptions
# ────────────────────────────────────────────────────────────

class TestTypeCheckingNoExceptions:
    """Property: For any valid DNA sequence, individual predicate checks
    and evaluate_all_predicates do not raise exceptions."""

    @given(seq=valid_coding_seq_nonempty)
    @settings(max_examples=50, deadline=5000)
    def test_check_no_stop_codons_no_exception(self, seq):
        """check_no_stop_codons never raises on valid coding sequences."""
        result = check_no_stop_codons(seq)
        assert result.predicate == "NoStopCodons"
        assert isinstance(result.passed, bool)

    @given(seq=valid_coding_seq_nonempty)
    @settings(max_examples=30, deadline=5000)
    def test_check_no_cpg_island_no_exception(self, seq):
        """check_no_cpg_island never raises on valid coding sequences."""
        result = check_no_cpg_island(seq)
        assert result.predicate == "NoCpGIsland"
        assert isinstance(result.passed, bool)

    @given(seq=valid_coding_seq_nonempty)
    @settings(max_examples=50, deadline=5000)
    def test_check_no_gt_dinucleotide_no_exception(self, seq):
        """check_no_gt_dinucleotide never raises on valid coding sequences."""
        result = check_no_gt_dinucleotide(seq)
        assert result.predicate == "NoGTDinucleotide"
        assert isinstance(result.passed, bool)

    @given(seq=valid_coding_seq_nonempty)
    @settings(max_examples=50, deadline=5000)
    def test_check_no_avoidable_gt_no_exception(self, seq):
        """check_no_avoidable_gt never raises on valid coding sequences."""
        result = check_no_avoidable_gt(seq)
        assert result.predicate == "NoGTDinucleotide"
        assert isinstance(result.passed, bool)

    @given(seq=valid_coding_seq_nonempty)
    @settings(max_examples=50, deadline=5000)
    def test_check_valid_coding_seq_no_exception(self, seq):
        """check_valid_coding_seq never raises on valid coding sequences."""
        result = check_valid_coding_seq(seq)
        assert result.predicate == "ValidCodingSeq"
        assert isinstance(result.passed, bool)

    @given(seq=valid_coding_seq_nonempty)
    @settings(max_examples=30, deadline=5000)
    def test_check_no_cryptic_promoter_no_exception(self, seq):
        """check_no_cryptic_promoter never raises on valid coding sequences."""
        result = check_no_cryptic_promoter(seq)
        assert result.predicate == "NoCrypticPromoter"
        assert isinstance(result.passed, bool)

    @given(seq=valid_coding_seq_nonempty)
    @settings(max_examples=30, deadline=5000)
    def test_check_no_unexpected_tm_domain_no_exception(self, seq):
        """check_no_unexpected_tm_domain never raises on valid coding sequences."""
        result = check_no_unexpected_tm_domain(seq)
        assert result.predicate == "NoUnexpectedTMDomain"
        assert isinstance(result.passed, bool)

    @given(seq=valid_coding_seq_nonempty)
    @settings(max_examples=30, deadline=5000)
    def test_check_mrna_secondary_structure_no_exception(self, seq):
        """check_mrna_secondary_structure never raises on valid coding sequences."""
        result = check_mrna_secondary_structure(seq, use_viennarna=False)
        assert result.predicate == "mRNASecondaryStructure"
        assert isinstance(result.passed, bool)

    @given(seq=valid_coding_seq_nonempty)
    @settings(max_examples=20, deadline=5000)
    def test_evaluate_all_predicates_no_exception(self, seq):
        """evaluate_all_predicates never raises on valid coding sequences."""
        results = evaluate_all_predicates(seq)
        assert isinstance(results, list)
        assert len(results) > 0

    # Also test with arbitrary DNA strings (not necessarily valid coding)
    @given(seq=dna_seq)
    @settings(max_examples=40, deadline=5000)
    def test_check_no_stop_codons_arbitrary_dna(self, seq):
        """check_no_stop_codons never raises on arbitrary DNA strings."""
        result = check_no_stop_codons(seq)
        assert isinstance(result.passed, bool)

    @given(seq=dna_seq)
    @settings(max_examples=30, deadline=5000)
    def test_check_no_cpg_island_arbitrary_dna(self, seq):
        """check_no_cpg_island never raises on arbitrary DNA strings."""
        # For very short sequences the default window=200 means no windows at all
        result = check_no_cpg_island(seq)
        assert isinstance(result.passed, bool)

    @given(seq=dna_seq)
    @settings(max_examples=40, deadline=5000)
    def test_check_valid_coding_seq_arbitrary_dna(self, seq):
        """check_valid_coding_seq never raises on arbitrary DNA strings."""
        result = check_valid_coding_seq(seq)
        assert isinstance(result.passed, bool)

    @given(seq=dna_seq)
    @settings(max_examples=30, deadline=5000)
    def test_check_no_cryptic_promoter_arbitrary_dna(self, seq):
        """check_no_cryptic_promoter never raises on arbitrary DNA strings."""
        result = check_no_cryptic_promoter(seq)
        assert isinstance(result.passed, bool)


# ────────────────────────────────────────────────────────────
# Property 4: combined_verdict returns a valid Verdict
# ────────────────────────────────────────────────────────────

class TestCombinedVerdictReturnsValidVerdict:
    """Property: combined_verdict always returns a valid Verdict enum member."""

    @given(vs=verdict_lists)
    def test_combined_verdict_returns_verdict(self, vs):
        """combined_verdict always returns a Verdict instance."""
        result = combined_verdict(vs)
        assert isinstance(result, Verdict), (
            f"Expected Verdict, got {type(result)}"
        )

    @given(vs=verdict_lists)
    def test_combined_verdict_name_in_valid_set(self, vs):
        """The returned Verdict's name is one of the five canonical names."""
        valid_names = {"PASS", "FAIL", "UNCERTAIN", "LIKELY_PASS", "LIKELY_FAIL"}
        result = combined_verdict(vs)
        assert result.name in valid_names, (
            f"Verdict name {result.name} not in {valid_names}"
        )

    @given(v=verdict)
    def test_single_verdict_returns_same(self, v):
        """combined_verdict([v]) == v for any single Verdict."""
        assert combined_verdict([v]) == v

    def test_empty_list_returns_uncertain(self):
        """combined_verdict([]) == UNCERTAIN (convention for no evidence)."""
        assert combined_verdict([]) == Verdict.UNCERTAIN

    @given(vs=verdict_lists)
    def test_any_fail_yields_fail(self, vs):
        """If any element is FAIL, combined_verdict returns FAIL."""
        if Verdict.FAIL in vs:
            assert combined_verdict(vs) == Verdict.FAIL

    @given(vs=st.lists(st.just(Verdict.PASS), min_size=1, max_size=20))
    def test_all_pass_yields_pass(self, vs):
        """If all elements are PASS, combined_verdict returns PASS."""
        assert combined_verdict(vs) == Verdict.PASS

    @given(vs=verdict_lists)
    def test_combined_verdict_worst_link(self, vs):
        """Combined verdict equals the minimum element by _VERDICT_ORDER."""
        if not vs:
            return  # empty case handled separately
        from biocompiler.shared.types import _VERDICT_ORDER
        min_verdict = min(vs, key=lambda v: _VERDICT_ORDER[v])
        assert combined_verdict(vs) == min_verdict


# ────────────────────────────────────────────────────────────
# Additional structural properties
# ────────────────────────────────────────────────────────────

class TestStructuralProperties:
    """Additional structural invariants for the type_system module."""

    def test_codon_table_has_64_entries(self):
        """The standard genetic code has exactly 64 codons."""
        assert len(CODON_TABLE) == 64, f"Expected 64 codons, got {len(CODON_TABLE)}"

    def test_all_codons_are_three_bases(self):
        """Every codon in CODON_TABLE is exactly 3 bases long."""
        for codon in CODON_TABLE:
            assert len(codon) == 3, f"Codon {codon!r} is not 3 bases"

    def test_all_codon_bases_are_acgt(self):
        """Every base in every codon is A, C, G, or T."""
        for codon in CODON_TABLE:
            for base in codon:
                assert base in "ACGT", f"Invalid base {base!r} in codon {codon!r}"

    def test_aa_to_codons_is_inverse_of_codon_table(self):
        """AA_TO_CODONS is a valid reverse mapping of CODON_TABLE."""
        for codon, aa in CODON_TABLE.items():
            assert codon in AA_TO_CODONS.get(aa, []), (
                f"Codon {codon} not found in AA_TO_CODONS[{aa!r}]"
            )
        for aa, codons in AA_TO_CODONS.items():
            for codon in codons:
                assert CODON_TABLE[codon] == aa, (
                    f"CODON_TABLE[{codon!r}] = {CODON_TABLE[codon]!r}, expected {aa!r}"
                )

    @given(codon=valid_codons)
    def test_valid_codon_in_table(self, codon):
        """Every codon generated by the valid_codons strategy is in CODON_TABLE."""
        assert codon in CODON_TABLE

    @given(seq=valid_coding_seq)
    def test_valid_coding_seq_passes_valid_coding_check(self, seq):
        """Sequences built from CODON_TABLE keys pass check_valid_coding_seq."""
        result = check_valid_coding_seq(seq)
        assert result.passed is True, (
            f"Valid coding sequence failed check: {result.details}"
        )
        assert result.verdict == Verdict.PASS
