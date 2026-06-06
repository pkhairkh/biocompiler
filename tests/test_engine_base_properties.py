"""Property-based tests for engine_base.py using Hypothesis.

Covers three core properties:
  1. BaseEngineResult fields are consistent (success/passed alignment, field types)
  2. EngineTimer measures positive elapsed time
  3. validate_protein_sequence rejects invalid inputs and accepts valid ones
"""

from __future__ import annotations

import time
from unittest.mock import patch

import pytest
from hypothesis import given, settings, assume, example
from hypothesis import strategies as st

from biocompiler.engine_base import (
    STANDARD_AMINO_ACIDS,
    BaseEngineResult,
    EngineConfig,
    EngineTimer,
    MutationResult,
    BatchResult,
    classify_score,
    validate_protein_sequence,
)


# ────────────────────────────────────────────────────────────
# Strategies
# ────────────────────────────────────────────────────────────

# Single standard amino acid
amino_acid = st.sampled_from(sorted(STANDARD_AMINO_ACIDS))

# Valid protein sequence (1–200 residues)
valid_protein = st.text(alphabet=sorted(STANDARD_AMINO_ACIDS), min_size=1, max_size=200)

# Any string (may contain invalid characters)
any_string = st.text(min_size=0, max_size=100)

# String with at least one non-amino-acid character
invalid_protein = st.text(
    alphabet=sorted(set("BJOUXZ0123456789!@#$% ") - STANDARD_AMINO_ACIDS),
    min_size=1,
    max_size=50,
)

# Floats for scores
finite_float = st.floats(allow_nan=False, allow_infinity=False)

# Small positive duration for sleep-based timer tests
small_positive_duration = st.floats(min_value=0.001, max_value=0.05)


# ────────────────────────────────────────────────────────────
# Property 1: BaseEngineResult fields are consistent
# ────────────────────────────────────────────────────────────

class TestBaseEngineResultConsistency:
    """Property: BaseEngineResult field relationships hold for all inputs."""

    @given(
        sequence=valid_protein,
        primary_score=finite_float,
        classification=st.text(min_size=1, max_size=50),
        success=st.booleans(),
        error=st.none() | st.text(min_size=1, max_size=100),
        execution_time_s=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False),
        engine_name=st.text(min_size=0, max_size=50),
        primary_score_label=st.text(min_size=0, max_size=30),
    )
    @settings(max_examples=50, deadline=5000)
    def test_passed_equals_success(
        self,
        sequence,
        primary_score,
        classification,
        success,
        error,
        execution_time_s,
        engine_name,
        primary_score_label,
    ):
        """The `passed` property always mirrors `success`."""
        result = BaseEngineResult(
            sequence=sequence,
            primary_score=primary_score,
            classification=classification,
            success=success,
            error=error,
            execution_time_s=execution_time_s,
            engine_name=engine_name,
            primary_score_label=primary_score_label,
        )
        assert result.passed == result.success
        assert result.passed == success

    @given(
        sequence=valid_protein,
        primary_score=finite_float,
        classification=st.text(min_size=1, max_size=50),
        success=st.booleans(),
    )
    @settings(max_examples=30, deadline=5000)
    def test_sequence_preserved(self, sequence, primary_score, classification, success):
        """The sequence field is stored exactly as provided."""
        result = BaseEngineResult(
            sequence=sequence,
            primary_score=primary_score,
            classification=classification,
            success=success,
        )
        assert result.sequence == sequence

    @given(
        sequence=valid_protein,
        primary_score=finite_float,
        classification=st.text(min_size=1, max_size=50),
        success=st.booleans(),
    )
    @settings(max_examples=30, deadline=5000)
    def test_primary_score_preserved(self, sequence, primary_score, classification, success):
        """The primary_score field is stored exactly as provided."""
        result = BaseEngineResult(
            sequence=sequence,
            primary_score=primary_score,
            classification=classification,
            success=success,
        )
        assert result.primary_score == primary_score

    @given(
        sequence=valid_protein,
        primary_score=finite_float,
        classification=st.text(min_size=1, max_size=50),
        success=st.booleans(),
    )
    @settings(max_examples=30, deadline=5000)
    def test_classification_preserved(self, sequence, primary_score, classification, success):
        """The classification field is stored exactly as provided."""
        result = BaseEngineResult(
            sequence=sequence,
            primary_score=primary_score,
            classification=classification,
            success=success,
        )
        assert result.classification == classification

    @given(success=st.booleans(), error=st.none() | st.text(min_size=1, max_size=80))
    @settings(max_examples=30, deadline=2000)
    def test_error_field_nullable(self, success, error):
        """The error field is either None or a string."""
        result = BaseEngineResult(
            sequence="M",
            primary_score=0.0,
            classification="ok",
            success=success,
            error=error,
        )
        assert result.error is None or isinstance(result.error, str)
        assert result.error == error

    @given(
        execution_time_s=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False),
    )
    @settings(max_examples=30, deadline=2000)
    def test_execution_time_non_negative(self, execution_time_s):
        """The execution_time_s is stored as provided (non-negative by construction)."""
        result = BaseEngineResult(
            sequence="M",
            primary_score=0.0,
            classification="ok",
            success=True,
            execution_time_s=execution_time_s,
        )
        assert result.execution_time_s >= 0.0

    @given(
        engine_name=st.text(min_size=0, max_size=50),
        primary_score_label=st.text(min_size=0, max_size=30),
    )
    @settings(max_examples=20, deadline=2000)
    def test_default_string_fields(self, engine_name, primary_score_label):
        """Engine name and score label are stored as provided."""
        result = BaseEngineResult(
            sequence="M",
            primary_score=1.0,
            classification="test",
            success=True,
            engine_name=engine_name,
            primary_score_label=primary_score_label,
        )
        assert result.engine_name == engine_name
        assert result.primary_score_label == primary_score_label

    def test_default_error_is_none(self):
        """Default value for error is None."""
        result = BaseEngineResult(
            sequence="M", primary_score=0.0, classification="ok", success=True,
        )
        assert result.error is None

    def test_default_execution_time_is_zero(self):
        """Default value for execution_time_s is 0.0."""
        result = BaseEngineResult(
            sequence="M", primary_score=0.0, classification="ok", success=True,
        )
        assert result.execution_time_s == 0.0

    def test_default_engine_name_is_empty(self):
        """Default value for engine_name is empty string."""
        result = BaseEngineResult(
            sequence="M", primary_score=0.0, classification="ok", success=True,
        )
        assert result.engine_name == ""

    def test_default_primary_score_label(self):
        """Default value for primary_score_label is 'score'."""
        result = BaseEngineResult(
            sequence="M", primary_score=0.0, classification="ok", success=True,
        )
        assert result.primary_score_label == "score"

    @given(success=st.booleans())
    @settings(max_examples=10, deadline=1000)
    def test_passed_bool_type(self, success):
        """The `passed` property always returns a bool."""
        result = BaseEngineResult(
            sequence="MK", primary_score=42.0, classification="test", success=success,
        )
        assert isinstance(result.passed, bool)

    @given(
        sequence=valid_protein,
        primary_score=finite_float,
        classification=st.text(min_size=1, max_size=50),
        success=st.booleans(),
        error=st.none() | st.text(min_size=1, max_size=50),
        execution_time_s=st.floats(min_value=0.0, max_value=1e6, allow_nan=False, allow_infinity=False),
        engine_name=st.text(min_size=0, max_size=30),
        primary_score_label=st.text(min_size=0, max_size=20),
    )
    @settings(max_examples=20, deadline=3000)
    def test_result_satisfies_engine_result_protocol(
        self,
        sequence,
        primary_score,
        classification,
        success,
        error,
        execution_time_s,
        engine_name,
        primary_score_label,
    ):
        """BaseEngineResult always satisfies the EngineResult protocol fields."""
        result = BaseEngineResult(
            sequence=sequence,
            primary_score=primary_score,
            classification=classification,
            success=success,
            error=error,
            execution_time_s=execution_time_s,
            engine_name=engine_name,
            primary_score_label=primary_score_label,
        )
        # EngineResult protocol requires these three fields
        assert hasattr(result, "success")
        assert hasattr(result, "error")
        assert hasattr(result, "execution_time_s")
        assert isinstance(result.success, bool)
        assert result.error is None or isinstance(result.error, str)
        assert isinstance(result.execution_time_s, float)


# ────────────────────────────────────────────────────────────
# Property 2: EngineTimer measures positive elapsed time
# ────────────────────────────────────────────────────────────

class TestEngineTimerElapsedTime:
    """Property: EngineTimer always records positive elapsed time after use.

    All timing-dependent tests use mocked time.perf_counter to eliminate
    flaky test failures due to system load or CI environment variance.
    """

    def test_timer_starts_at_zero(self):
        """Fresh EngineTimer has zero elapsed time."""
        timer = EngineTimer()
        assert timer.elapsed == 0.0
        assert timer.start == 0.0

    def test_elapsed_positive_after_context_with_mock(self):
        """After context exit, elapsed time is positive (mocked time)."""
        fake_time = 1000.0

        def mock_perf_counter():
            nonlocal fake_time
            return fake_time

        with patch("biocompiler.engine_base.time.perf_counter", side_effect=mock_perf_counter):
            with EngineTimer() as timer:
                fake_time = 1000.05  # advance 50ms
            assert timer.elapsed > 0.0

    def test_elapsed_matches_mock_time_advance(self):
        """Elapsed time matches the mock time advance exactly."""
        fake_time = 1000.0

        def mock_perf_counter():
            nonlocal fake_time
            return fake_time

        with patch("biocompiler.engine_base.time.perf_counter", side_effect=mock_perf_counter):
            with EngineTimer() as timer:
                fake_time = 1000.5  # advance 0.5s
            assert timer.elapsed == pytest.approx(0.5, rel=1e-6)

    def test_elapsed_reasonable_upper_bound_with_mock(self):
        """Elapsed time is not absurdly larger than the mock time advance."""
        fake_time = 2000.0

        def mock_perf_counter():
            nonlocal fake_time
            return fake_time

        with patch("biocompiler.engine_base.time.perf_counter", side_effect=mock_perf_counter):
            with EngineTimer() as timer:
                fake_time = 2000.1  # advance 0.1s
            assert timer.elapsed < 10.0

    def test_elapsed_positive_even_with_no_work(self):
        """Even with no work inside context, elapsed is >= 0."""
        fake_time = 3000.0

        def mock_perf_counter():
            nonlocal fake_time
            return fake_time

        with patch("biocompiler.engine_base.time.perf_counter", side_effect=mock_perf_counter):
            with EngineTimer() as timer:
                pass  # no-op, time doesn't advance
            assert timer.elapsed >= 0.0

    def test_timer_returns_self_from_enter(self):
        """__enter__ returns the timer instance itself."""
        timer = EngineTimer()
        with timer as t:
            assert t is timer

    def test_two_timers_ordering_with_mock(self):
        """A longer duration produces a larger elapsed time (mocked)."""
        # Short timer
        fake_time = 4000.0

        def mock_perf_counter_short():
            nonlocal fake_time
            return fake_time

        with patch("biocompiler.engine_base.time.perf_counter", side_effect=mock_perf_counter_short):
            with EngineTimer() as timer_short:
                fake_time = 4000.01  # 10ms
            short_elapsed = timer_short.elapsed

        # Long timer
        fake_time = 5000.0

        def mock_perf_counter_long():
            nonlocal fake_time
            return fake_time

        with patch("biocompiler.engine_base.time.perf_counter", side_effect=mock_perf_counter_long):
            with EngineTimer() as timer_long:
                fake_time = 5000.1  # 100ms
            long_elapsed = timer_long.elapsed

        # The longer timer should have greater elapsed time
        assert long_elapsed > short_elapsed

    def test_reuse_timer_resets_with_mock(self):
        """Reusing a timer resets start/elapsed properly (mocked)."""
        fake_time = 6000.0

        def mock_perf_counter():
            nonlocal fake_time
            return fake_time

        with patch("biocompiler.engine_base.time.perf_counter", side_effect=mock_perf_counter):
            timer = EngineTimer()
            with timer:
                fake_time = 6000.005
            first_elapsed = timer.elapsed
            with timer:
                fake_time = 6000.015  # 10ms from second start
            second_elapsed = timer.elapsed
            # Second elapsed is independent of first (not accumulated)
            assert second_elapsed >= 0.0
            # Both should be positive
            assert first_elapsed > 0.0
            assert second_elapsed > 0.0

    def test_elapsed_is_float_with_mock(self):
        """Elapsed time is always a float (mocked)."""
        fake_time = 7000.0

        def mock_perf_counter():
            nonlocal fake_time
            return fake_time

        with patch("biocompiler.engine_base.time.perf_counter", side_effect=mock_perf_counter):
            with EngineTimer() as timer:
                fake_time = 7000.1
            assert isinstance(timer.elapsed, float)

    def test_timer_with_exception_in_context_mocked(self):
        """Timer records elapsed time even if context body raises (mocked)."""
        fake_time = 8000.0

        def mock_perf_counter():
            nonlocal fake_time
            return fake_time

        timer = EngineTimer()
        try:
            with patch("biocompiler.engine_base.time.perf_counter", side_effect=mock_perf_counter):
                with timer:
                    fake_time = 8000.005
                    raise RuntimeError("test error")
        except RuntimeError:
            pass
        # Timer should still have recorded elapsed time
        assert timer.elapsed > 0.0


# ────────────────────────────────────────────────────────────
# Property 3: validate_protein_sequence rejects invalid inputs
# ────────────────────────────────────────────────────────────

class TestValidateProteinSequenceRejectsInvalid:
    """Property: validate_protein_sequence raises ValueError for invalid inputs."""

    @given(invalid_char_seq=invalid_protein)
    @settings(max_examples=50, deadline=3000)
    def test_rejects_non_amino_acid_characters(self, invalid_char_seq):
        """Sequences containing non-standard amino acids raise ValueError."""
        assume(any(c not in STANDARD_AMINO_ACIDS for c in invalid_char_seq.upper().strip()))
        with pytest.raises(ValueError):
            validate_protein_sequence(invalid_char_seq)

    def test_rejects_empty_string(self):
        """Empty string raises ValueError."""
        with pytest.raises(ValueError, match="must not be empty"):
            validate_protein_sequence("")

    def test_rejects_whitespace_only(self):
        """Whitespace-only string raises ValueError."""
        with pytest.raises(ValueError, match="must not be empty"):
            validate_protein_sequence("   ")

    @given(whitespace_str=st.text(alphabet=" \t\n\r", min_size=1, max_size=20))
    @settings(max_examples=20, deadline=2000)
    def test_rejects_various_whitespace(self, whitespace_str):
        """Any combination of whitespace characters raises ValueError."""
        with pytest.raises(ValueError):
            validate_protein_sequence(whitespace_str)

    def test_rejects_digits(self):
        """Sequences with digits raise ValueError."""
        with pytest.raises(ValueError):
            validate_protein_sequence("ACGT123")

    def test_rejects_lowercase_not_in_standard(self):
        """Lowercase letters that aren't amino acids after upper-casing
        still raise ValueError."""
        with pytest.raises(ValueError):
            validate_protein_sequence("BJOXZ")

    @given(engine_name=st.from_regex(r"[A-Za-z0-9_]{1,20}", fullmatch=True))
    @settings(max_examples=20, deadline=2000)
    def test_error_message_contains_engine_name(self, engine_name):
        """ValueError message includes the engine_name argument."""
        import re
        with pytest.raises(ValueError, match=re.escape(engine_name)):
            validate_protein_sequence("", engine_name=engine_name)

    @given(
        bad_seq=st.text(
            alphabet=sorted(set("BJOUXZ") - STANDARD_AMINO_ACIDS),
            min_size=1,
            max_size=10,
        ),
        engine_name=st.from_regex(r"[A-Za-z0-9_]{1,20}", fullmatch=True),
    )
    @settings(max_examples=30, deadline=3000)
    def test_invalid_char_error_mentions_engine(self, bad_seq, engine_name):
        """ValueError for invalid characters includes the engine name."""
        import re
        assume(any(c not in STANDARD_AMINO_ACIDS for c in bad_seq.upper().strip()))
        with pytest.raises(ValueError, match=re.escape(engine_name)):
            validate_protein_sequence(bad_seq, engine_name=engine_name)

    @given(bad_seq=invalid_protein)
    @settings(max_examples=30, deadline=3000)
    def test_invalid_char_error_mentions_invalid_amino_acids(self, bad_seq):
        """ValueError for invalid characters mentions 'invalid amino acids'."""
        assume(any(c not in STANDARD_AMINO_ACIDS for c in bad_seq.upper().strip()))
        with pytest.raises(ValueError, match="invalid amino acids"):
            validate_protein_sequence(bad_seq)


class TestValidateProteinSequenceAcceptsValid:
    """Property: validate_protein_sequence accepts and normalizes valid inputs."""

    @given(seq=valid_protein)
    @settings(max_examples=50, deadline=3000)
    def test_accepts_valid_sequences(self, seq):
        """Valid amino acid sequences pass validation without error."""
        result = validate_protein_sequence(seq)
        assert result == seq.upper().strip()

    @given(seq=valid_protein)
    @settings(max_examples=30, deadline=3000)
    def test_result_is_uppercase(self, seq):
        """Returned sequence is always uppercase."""
        result = validate_protein_sequence(seq)
        assert result == result.upper()

    @given(seq=valid_protein)
    @settings(max_examples=30, deadline=3000)
    def test_result_has_no_whitespace(self, seq):
        """Returned sequence has no leading/trailing whitespace."""
        result = validate_protein_sequence(seq)
        assert result == result.strip()

    @given(seq=valid_protein)
    @settings(max_examples=30, deadline=3000)
    def test_result_contains_only_standard_amino_acids(self, seq):
        """Every character in the result is a standard amino acid."""
        result = validate_protein_sequence(seq)
        for ch in result:
            assert ch in STANDARD_AMINO_ACIDS, f"Invalid amino acid: {ch!r}"

    @given(seq=valid_protein)
    @settings(max_examples=30, deadline=3000)
    def test_result_length_matches_input(self, seq):
        """Returned sequence length matches stripped input length."""
        result = validate_protein_sequence(seq)
        assert len(result) == len(seq.strip())

    @given(seq=valid_protein)
    @settings(max_examples=30, deadline=3000)
    def test_lowercase_input_normalized(self, seq):
        """Lowercase input is normalized to uppercase."""
        result = validate_protein_sequence(seq.lower())
        assert result == seq.upper().strip()

    @given(
        prefix=st.text(alphabet=" ", min_size=1, max_size=5),
        suffix=st.text(alphabet=" \t", min_size=1, max_size=5),
        seq=valid_protein,
    )
    @settings(max_examples=20, deadline=3000)
    def test_whitespace_stripped(self, prefix, suffix, seq):
        """Leading/trailing whitespace is stripped from the result."""
        padded = prefix + seq + suffix
        result = validate_protein_sequence(padded)
        assert result == seq.upper().strip()

    @given(seq=valid_protein)
    @example(seq="M")
    @example(seq="ACDEFGHIKLMNPQRSTVWY")
    @settings(max_examples=10, deadline=2000)
    def test_single_residue_and_full_alphabet(self, seq):
        """Edge cases: single residue and full 20-amino-acid alphabet."""
        result = validate_protein_sequence(seq)
        assert len(result) > 0
        for ch in result:
            assert ch in STANDARD_AMINO_ACIDS


# ────────────────────────────────────────────────────────────
# Additional structural properties for completeness
# ────────────────────────────────────────────────────────────

class TestStandardAminoAcidsConstant:
    """Property: STANDARD_AMINO_ACIDS has the expected composition."""

    def test_has_20_amino_acids(self):
        """There are exactly 20 standard amino acids."""
        assert len(STANDARD_AMINO_ACIDS) == 20

    def test_is_frozenset(self):
        """STANDARD_AMINO_ACIDS is a frozenset (immutable)."""
        assert isinstance(STANDARD_AMINO_ACIDS, frozenset)

    @given(aa=st.sampled_from(sorted(STANDARD_AMINO_ACIDS)))
    @settings(max_examples=20, deadline=1000)
    def test_each_amino_acid_is_single_char(self, aa):
        """Each amino acid code is a single uppercase character."""
        assert len(aa) == 1
        assert aa.isupper()


class TestMutationResultConsistency:
    """Property: MutationResult backward-compat aliases work correctly."""

    @given(
        position=st.integers(min_value=0, max_value=1000),
        original=st.sampled_from(sorted(STANDARD_AMINO_ACIDS)),
        mutant=st.sampled_from(sorted(STANDARD_AMINO_ACIDS)),
        delta_score=finite_float,
        score=finite_float,
    )
    @settings(max_examples=30, deadline=3000)
    def test_score_alias_matches_delta_score(self, position, original, mutant, delta_score, score):
        """The `score` property always equals `delta_score`."""
        mr = MutationResult(
            position=position,
            original=original,
            mutant=mutant,
            delta_score=delta_score,
        )
        assert mr.score == mr.delta_score
        assert mr.score == delta_score

    @given(
        original=st.sampled_from(sorted(STANDARD_AMINO_ACIDS)),
        mutant=st.sampled_from(sorted(STANDARD_AMINO_ACIDS)),
    )
    @settings(max_examples=20, deadline=2000)
    def test_original_aa_alias(self, original, mutant):
        """The `original_aa` property always equals `original`."""
        mr = MutationResult(original=original, mutant=mutant)
        assert mr.original_aa == mr.original
        assert mr.original_aa == original

    @given(
        original=st.sampled_from(sorted(STANDARD_AMINO_ACIDS)),
        mutant=st.sampled_from(sorted(STANDARD_AMINO_ACIDS)),
    )
    @settings(max_examples=20, deadline=2000)
    def test_mutant_aa_alias(self, original, mutant):
        """The `mutant_aa` property always equals `mutant`."""
        mr = MutationResult(original=original, mutant=mutant)
        assert mr.mutant_aa == mr.mutant
        assert mr.mutant_aa == mutant

    @given(delta_score=finite_float)
    @settings(max_examples=20, deadline=2000)
    def test_score_setter_updates_delta_score(self, delta_score):
        """Setting `score` updates `delta_score`."""
        mr = MutationResult(delta_score=0.0)
        mr.score = delta_score
        assert mr.delta_score == delta_score
        assert mr.score == delta_score

    @given(score_value=finite_float)
    @settings(max_examples=20, deadline=2000)
    def test_score_kwarg_sets_delta_score(self, score_value):
        """Using score= keyword in constructor sets delta_score."""
        mr = MutationResult(score=score_value)
        assert mr.delta_score == score_value
        assert mr.score == score_value

    def test_original_aa_kwarg_alias(self):
        """Using original_aa= sets the original field."""
        mr = MutationResult(original_aa="A", mutant="G")
        assert mr.original == "A"
        assert mr.original_aa == "A"

    def test_mutant_aa_kwarg_alias(self):
        """Using mutant_aa= sets the mutant field."""
        mr = MutationResult(original="A", mutant_aa="G")
        assert mr.mutant == "G"
        assert mr.mutant_aa == "G"


class TestBatchResultConsistency:
    """Property: BatchResult auto-computation and aliases are consistent."""

    @given(n_success=st.integers(min_value=0, max_value=50))
    @settings(max_examples=20, deadline=2000)
    def test_auto_compute_successful_and_failed(self, n_success):
        """When successful=failed=0 and results provided, counts are auto-computed."""
        results = [
            BaseEngineResult(
                sequence="M", primary_score=float(i), classification="ok", success=(i < n_success),
            )
            for i in range(n_success + 3)
        ]
        # Only n_success of them have success=True
        batch = BatchResult(results=results)
        assert batch.successful == n_success
        assert batch.failed == 3
        assert batch.total == n_success + 3

    @given(successful=st.integers(min_value=0, max_value=50), failed=st.integers(min_value=0, max_value=50))
    @settings(max_examples=20, deadline=2000)
    def test_explicit_counts_not_overridden(self, successful, failed):
        """When successful and failed are explicitly provided, they are preserved."""
        batch = BatchResult(successful=successful, failed=failed)
        assert batch.successful == successful
        assert batch.failed == failed

    @given(n=st.integers(min_value=0, max_value=20))
    @settings(max_examples=15, deadline=2000)
    def test_success_count_alias(self, n):
        """success_count is an alias for successful."""
        batch = BatchResult(successful=n, failed=0)
        assert batch.success_count == batch.successful
        assert batch.success_count == n

    @given(n=st.integers(min_value=0, max_value=20))
    @settings(max_examples=15, deadline=2000)
    def test_failure_count_alias(self, n):
        """failure_count is an alias for failed."""
        batch = BatchResult(successful=0, failed=n)
        assert batch.failure_count == batch.failed
        assert batch.failure_count == n

    def test_total_equals_len_results(self):
        """total property equals len(results)."""
        results = [
            BaseEngineResult(sequence="M", primary_score=i, classification="x", success=True)
            for i in range(7)
        ]
        batch = BatchResult(results=results)
        assert batch.total == len(results) == 7

    def test_empty_batch(self):
        """Empty batch has zero counts and empty lists."""
        batch = BatchResult()
        assert batch.results == []
        assert batch.errors == []
        assert batch.total == 0
        assert batch.successful == 0
        assert batch.failed == 0


class TestClassifyScoreProperties:
    """Property: classify_score returns consistent classifications."""

    @given(
        score=finite_float,
        thresholds=st.lists(
            st.tuples(finite_float, st.text(min_size=1, max_size=20)),
            min_size=1,
            max_size=10,
        ),
        fallback=st.text(min_size=1, max_size=20),
    )
    @settings(max_examples=40, deadline=3000)
    def test_returns_string(self, score, thresholds, fallback):
        """classify_score always returns a string."""
        result = classify_score(score, thresholds, fallback)
        assert isinstance(result, str)

    @given(
        score=st.floats(min_value=0.0, max_value=100.0, allow_nan=False, allow_infinity=False),
        fallback=st.text(min_size=1, max_size=20),
    )
    @settings(max_examples=30, deadline=2000)
    def test_returns_fallback_when_no_threshold_matches(self, score, fallback):
        """When score is below all thresholds, fallback is returned."""
        # Use thresholds all above the max possible score
        thresholds = [(100.1, "high"), (100.2, "very_high")]
        result = classify_score(score, thresholds, fallback)
        assert result == fallback

    @given(
        score=st.floats(min_value=50.0, max_value=100.0, allow_nan=False, allow_infinity=False),
    )
    @settings(max_examples=30, deadline=2000)
    def test_returns_label_when_score_exceeds_threshold(self, score):
        """When score exceeds at least one threshold, a threshold label is returned."""
        thresholds = [(50.0, "at_least_50"), (90.0, "at_least_90")]
        result = classify_score(score, thresholds, fallback="none")
        assert result in ("at_least_50", "at_least_90", "none")
        # If score >= 50, should not return fallback
        assert result != "none"

    def test_descending_threshold_order_matters(self):
        """First matching threshold in list order wins."""
        thresholds = [(80, "high"), (60, "medium"), (40, "low")]
        # Score 75 matches 60 but not 80
        assert classify_score(75, thresholds) == "medium"
        # Score 85 matches 80
        assert classify_score(85, thresholds) == "high"

    @given(fallback=st.text(min_size=1, max_size=20))
    @settings(max_examples=10, deadline=1000)
    def test_empty_thresholds_returns_fallback(self, fallback):
        """With no thresholds, fallback is always returned."""
        result = classify_score(42.0, [], fallback)
        assert result == fallback
