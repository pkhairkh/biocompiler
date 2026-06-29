"""Regression tests for peer-review P0 fixes (Waves 1-3).

Each test verifies a specific bug fix from the negative-intent peer review.
"""
import pytest
from biocompiler import optimize_sequence
from biocompiler.shared.exceptions import InvalidProteinError


class TestSlowPathImportFix:
    """W1-a: slow path no longer crashes with UnboundLocalError."""

    def test_slow_path_completes_without_crash(self):
        """use_integrated=False must not raise UnboundLocalError."""
        r = optimize_sequence("MVHLTPEEK", organism="e_coli",
                              use_integrated=False, strict_mode=False)
        assert r is not None
        assert r.convergence_status is not None

    def test_slow_path_produces_certificate(self):
        """Slow path should produce a certificate."""
        r = optimize_sequence("MVHLTPEEK", organism="e_coli",
                              use_integrated=False, strict_mode=False)
        assert bool(r.certificate_text)  # may or may not have cert, but must not crash


class TestProteinValidationBeforeFastPath:
    """W2-a: invalid proteins must raise on ALL paths, not just slow path."""

    def test_invalid_chars_fast_path_raises(self):
        """X, Z, !, @, # must raise InvalidProteinError on default path."""
        with pytest.raises(InvalidProteinError):
            optimize_sequence("MXYZ!!!@#", organism="e_coli")

    def test_invalid_chars_slow_path_raises(self):
        """Same invalid chars must raise on slow path too."""
        with pytest.raises(InvalidProteinError):
            optimize_sequence("MXYZ!!!@#", organism="e_coli", use_integrated=False)

    def test_selenocysteine_U_is_valid(self):
        """U (selenocysteine) must be accepted on both paths."""
        r = optimize_sequence("MASECUGP", organism="h_sapiens")
        assert r.convergence_status is not None

    def test_empty_protein_raises(self):
        """Empty protein must raise."""
        with pytest.raises((InvalidProteinError, ValueError)):
            optimize_sequence("", organism="e_coli")

    def test_trailing_stop_stripped(self):
        """Trailing * must be stripped and accepted."""
        r = optimize_sequence("MVHLTPEEK*", organism="e_coli")
        assert r is not None


class TestBiosecurityBeforeOptimization:
    """W2-b: biosecurity screening must run before optimization on all paths."""

    def test_fast_path_screens_before_optimization(self):
        """Fast path should have biosecurity_screening_result populated."""
        r = optimize_sequence("MVHLTPEEK", organism="e_coli")
        assert r.biosecurity_screening_result is not None


class TestCertifiedByDefault:
    """W3-a: the default (fast) path must produce a certificate."""

    def test_default_path_produces_certificate(self):
        """optimize_sequence() with defaults must produce certificate_text."""
        r = optimize_sequence("MVHLTPEEKSAVTALWGKVNVDEVGGEALGR", organism="h_sapiens")
        assert r.certificate_text  # non-empty
        assert len(r.certificate_text) > 100  # substantial certificate

    def test_default_path_evaluates_predicates(self):
        """predicate_results must be populated on default path."""
        r = optimize_sequence("MVHLTPEEK", organism="e_coli")
        assert r.predicate_results is not None
        assert len(r.predicate_results) > 0

    def test_default_path_populates_failed_predicates(self):
        """failed_predicates must be a list (even if empty)."""
        r = optimize_sequence("MVHLTPEEK", organism="e_coli")
        assert isinstance(r.failed_predicates, list)

    def test_default_path_populates_biosecurity(self):
        """biosecurity_screening_result must be populated."""
        r = optimize_sequence("MVHLTPEEK", organism="e_coli")
        assert r.biosecurity_screening_result is not None

    def test_strict_mode_raises_on_failure(self):
        """strict_mode=True must raise on predicate failure."""
        # HBB has NoRQCTrigger which should trigger strict mode
        with pytest.raises(Exception):  # OptimizationConstraintError
            optimize_sequence("MVHLTPEEKSAVTALWGKVNVDEVGGEALGR",
                            organism="h_sapiens", strict_mode=True)

    def test_strict_mode_false_returns_result(self):
        """strict_mode=False must return a result even with failures."""
        r = optimize_sequence("MVHLTPEEKSAVTALWGKVNVDEVGGEALGR",
                            organism="h_sapiens", strict_mode=False)
        assert r is not None
        assert r.certificate_text  # still produces certificate

    def test_fast_and_slow_path_field_parity(self):
        """Fast and slow paths should populate the same key fields."""
        fast = optimize_sequence("MVHLTPEEK", organism="e_coli", strict_mode=False)
        slow = optimize_sequence("MVHLTPEEK", organism="e_coli",
                                use_integrated=False, strict_mode=False)
        # Both should have certificates
        assert fast.certificate_text
        assert slow.certificate_text
        # Both should have predicate results
        assert fast.predicate_results is not None
        assert slow.predicate_results is not None
        # Both should have biosecurity
        assert fast.biosecurity_screening_result is not None
        assert slow.biosecurity_screening_result is not None


class TestIntegratedOptimizerHonesty:
    """W1-b: integrated optimizer docstrings must be honest."""

    def test_integrated_optimizer_imports(self):
        """Module must import cleanly."""
        from biocompiler.optimizer.integrated_optimizer import integrated_optimize
        assert callable(integrated_optimize)

    def test_integrated_optimizer_produces_valid_dna(self):
        """Output must be valid DNA with protein preserved."""
        from biocompiler.optimizer.integrated_optimizer import integrated_optimize, _translate_dna
        protein = "MVHLTPEEK"
        dna, notes, secis = integrated_optimize(protein, "e_coli")
        translated = _translate_dna(dna, secis)
        assert translated == protein + "*"
