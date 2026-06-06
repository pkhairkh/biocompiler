"""
BioCompiler Exception Hierarchy Tests
======================================

Tests for the biocompiler.exceptions module — every custom exception class,
its attributes, message formatting, and inheritance chain.
"""

from __future__ import annotations

import pytest

from biocompiler.exceptions import (
    BioCompilerError,
    InvalidSequenceError,
    CertificateGenerationError,
    CertificateVerificationError,
    UnknownPredicateError,
    OptimizationError,
    UnsupportedOrganismError,
    InvalidProteinError,
    FileFormatError,
    SplicingError,
    MutagenesisError,
    EngineError,
    ESMFoldError,
    FoldXError,
    CamSolError,
    ImmunogenicityError,
    OptimizationConstraintError,
    BiosecurityError,
    TranslationVerificationError,
)


# ════════════════════════════════════════════════════════════════════════════
# Base class
# ════════════════════════════════════════════════════════════════════════════

class TestBioCompilerError:
    """Test the root exception class."""

    def test_is_exception(self):
        """BioCompilerError inherits from Exception."""
        assert issubclass(BioCompilerError, Exception)

    def test_message(self):
        """BioCompilerError carries a message string."""
        err = BioCompilerError("test message")
        assert str(err) == "test message"

    def test_can_be_caught_as_exception(self):
        """BioCompilerError is catchable as a plain Exception."""
        with pytest.raises(Exception):
            raise BioCompilerError("boom")


# ════════════════════════════════════════════════════════════════════════════
# InvalidSequenceError
# ════════════════════════════════════════════════════════════════════════════

class TestInvalidSequenceError:
    """Test InvalidSequenceError with position context."""

    def test_inherits_from_base(self):
        assert issubclass(InvalidSequenceError, BioCompilerError)

    def test_stores_sequence_and_invalid_chars(self):
        seq = "ATXGZ"
        bad = {"X", "Z"}
        err = InvalidSequenceError(seq, bad)
        assert err.sequence == seq
        assert err.invalid_chars == bad

    def test_message_contains_positions(self):
        err = InvalidSequenceError("ATXGZ", {"X", "Z"})
        msg = str(err)
        assert "Invalid DNA bases" in msg
        assert "2" in msg  # position of X
        assert "4" in msg  # position of Z

    def test_message_truncates_many_positions(self):
        """When > 10 invalid positions, message says '... (N total)'."""
        seq = "X" * 20
        err = InvalidSequenceError(seq, {"X"})
        msg = str(err)
        assert "total" in msg

    def test_can_be_caught_as_biocompiler_error(self):
        with pytest.raises(BioCompilerError):
            raise InvalidSequenceError("AXG", {"X"})

    def test_empty_invalid_chars(self):
        """Edge case: no invalid chars (shouldn't normally happen)."""
        err = InvalidSequenceError("ATG", set())
        assert err.invalid_chars == set()


# ════════════════════════════════════════════════════════════════════════════
# CertificateGenerationError
# ════════════════════════════════════════════════════════════════════════════

class TestCertificateGenerationError:

    def test_inherits_from_base(self):
        assert issubclass(CertificateGenerationError, BioCompilerError)

    def test_stores_failures(self):
        class _FakeResult:
            predicate = "test_pred"
            verdict = type("V", (), {"value": "FAIL"})()
        failures = [_FakeResult(), _FakeResult()]
        err = CertificateGenerationError(failures)
        assert err.failures == failures

    def test_message_mentions_count(self):
        class _FakeResult:
            predicate = "test_pred"
            verdict = type("V", (), {"value": "FAIL"})()
        err = CertificateGenerationError([_FakeResult(), _FakeResult()])
        msg = str(err)
        assert "2 predicate(s) failed" in msg


# ════════════════════════════════════════════════════════════════════════════
# CertificateVerificationError
# ════════════════════════════════════════════════════════════════════════════

class TestCertificateVerificationError:

    def test_stores_reasons(self):
        reasons = ["sig mismatch", "missing key"]
        err = CertificateVerificationError(reasons)
        assert err.reasons == reasons

    def test_message_joins_reasons(self):
        err = CertificateVerificationError(["reason A", "reason B"])
        msg = str(err)
        assert "reason A" in msg
        assert "reason B" in msg


# ════════════════════════════════════════════════════════════════════════════
# UnknownPredicateError
# ════════════════════════════════════════════════════════════════════════════

class TestUnknownPredicateError:

    def test_stores_predicate_name(self):
        err = UnknownPredicateError("MyPredicate")
        assert err.predicate_name == "MyPredicate"

    def test_message_contains_name(self):
        err = UnknownPredicateError("FooPred")
        msg = str(err)
        assert "FooPred" in msg
        assert "Register it" in msg


# ════════════════════════════════════════════════════════════════════════════
# OptimizationError
# ════════════════════════════════════════════════════════════════════════════

class TestOptimizationError:

    def test_with_reason_only(self):
        err = OptimizationError("no solution found")
        assert str(err) == "Optimization failed: no solution found"
        assert err.unsat_core is None

    def test_with_unsat_core(self):
        core = ["gc_content", "restriction_site"]
        err = OptimizationError("conflict", unsat_core=core)
        assert err.unsat_core == core

    def test_inherits_from_base(self):
        assert issubclass(OptimizationError, BioCompilerError)


# ════════════════════════════════════════════════════════════════════════════
# UnsupportedOrganismError
# ════════════════════════════════════════════════════════════════════════════

class TestUnsupportedOrganismError:

    def test_stores_organism_and_available(self):
        avail = ["E_coli", "human"]
        err = UnsupportedOrganismError("Alien", avail)
        assert err.organism == "Alien"
        assert err.available == avail

    def test_message_lists_available(self):
        err = UnsupportedOrganismError("Alien", ["E_coli", "human"])
        msg = str(err)
        assert "Alien" in msg
        assert "E_coli" in msg


# ════════════════════════════════════════════════════════════════════════════
# InvalidProteinError
# ════════════════════════════════════════════════════════════════════════════

class TestInvalidProteinError:

    def test_stores_protein_and_invalid_chars(self):
        err = InvalidProteinError("MVLZ", {"Z"})
        assert err.protein == "MVLZ"
        assert err.invalid_chars == {"Z"}

    def test_message_mentions_amino_acid(self):
        err = InvalidProteinError("MX", {"X"})
        assert "amino acid" in str(err).lower()


# ════════════════════════════════════════════════════════════════════════════
# FileFormatError
# ════════════════════════════════════════════════════════════════════════════

class TestFileFormatError:

    def test_stores_path_format_reason(self):
        err = FileFormatError("/tmp/test.fasta", "FASTA", "bad header")
        assert err.path == "/tmp/test.fasta"
        assert err.format_name == "FASTA"
        assert err.reason == "bad header"

    def test_message_includes_all_fields(self):
        err = FileFormatError("/tmp/gb.gbk", "GenBank", "truncated")
        msg = str(err)
        assert "/tmp/gb.gbk" in msg
        assert "GenBank" in msg
        assert "truncated" in msg


# ════════════════════════════════════════════════════════════════════════════
# SplicingError
# ════════════════════════════════════════════════════════════════════════════

class TestSplicingError:

    def test_stores_reason(self):
        err = SplicingError("NDFST overflow")
        assert err.reason == "NDFST overflow"

    def test_message_format(self):
        err = SplicingError("timeout")
        assert "Splicing computation error" in str(err)


# ════════════════════════════════════════════════════════════════════════════
# MutagenesisError
# ════════════════════════════════════════════════════════════════════════════

class TestMutagenesisError:

    def test_stores_reason_and_substitutions(self):
        err = MutagenesisError("stuck", substitutions_applied=7)
        assert err.reason == "stuck"
        assert err.substitutions_applied == 7

    def test_default_substitutions_zero(self):
        err = MutagenesisError("fail")
        assert err.substitutions_applied == 0

    def test_message_includes_substitution_count(self):
        err = MutagenesisError("loop", substitutions_applied=3)
        msg = str(err)
        assert "3 substitutions" in msg


# ════════════════════════════════════════════════════════════════════════════
# EngineError base class
# ════════════════════════════════════════════════════════════════════════════

class TestEngineError:

    def test_inherits_from_biocompiler_error(self):
        assert issubclass(EngineError, BioCompilerError)

    def test_default_engine(self):
        err = EngineError("crash")
        assert err.engine == "unknown"
        assert err.reason == "crash"

    def test_custom_engine(self):
        err = EngineError("timeout", engine="CustomEngine")
        assert err.engine == "CustomEngine"
        assert "[CustomEngine]" in str(err)

    def test_str_uses_message(self):
        err = EngineError("fail", engine="Test")
        assert str(err) == "[Test] fail"


# ════════════════════════════════════════════════════════════════════════════
# ESMFoldError
# ════════════════════════════════════════════════════════════════════════════

class TestESMFoldError:

    def test_inherits_from_engine_error(self):
        assert issubclass(ESMFoldError, EngineError)

    def test_reason_only(self):
        err = ESMFoldError("API unreachable")
        # err.reason stores the raw reason; the formatted message includes prefix
        assert "API unreachable" in str(err)
        assert err.protein is None

    def test_with_protein(self):
        err = ESMFoldError("bad input", protein="MKFLILLF")
        assert err.protein == "MKFLILLF"
        msg = str(err)
        assert "8" in msg  # len("MKFLILLF") = 8

    def test_engine_is_esmfold(self):
        err = ESMFoldError("fail")
        assert err.engine == "ESMFold"


# ════════════════════════════════════════════════════════════════════════════
# FoldXError
# ════════════════════════════════════════════════════════════════════════════

class TestFoldXError:

    def test_inherits_from_engine_error(self):
        assert issubclass(FoldXError, EngineError)

    def test_with_command(self):
        err = FoldXError("timeout", command="FoldX --command")
        assert err.command == "FoldX --command"
        msg = str(err)
        assert "FoldX --command" in msg

    def test_without_command(self):
        err = FoldXError("missing install")
        assert err.command is None


# ════════════════════════════════════════════════════════════════════════════
# CamSolError
# ════════════════════════════════════════════════════════════════════════════

class TestCamSolError:

    def test_inherits_from_engine_error(self):
        assert issubclass(CamSolError, EngineError)

    def test_with_protein(self):
        err = CamSolError("bad protein", protein="ACDE")
        assert err.protein == "ACDE"
        assert "4" in str(err)  # len("ACDE") = 4

    def test_engine_is_camsol(self):
        err = CamSolError("fail")
        assert err.engine == "CamSol"


# ════════════════════════════════════════════════════════════════════════════
# ImmunogenicityError
# ════════════════════════════════════════════════════════════════════════════

class TestImmunogenicityError:

    def test_inherits_from_engine_error(self):
        assert issubclass(ImmunogenicityError, EngineError)

    def test_with_protein(self):
        err = ImmunogenicityError("epitope fail", protein="MKF")
        assert err.protein == "MKF"

    def test_engine_is_immunogenicity(self):
        err = ImmunogenicityError("fail")
        assert err.engine == "Immunogenicity"


# ════════════════════════════════════════════════════════════════════════════
# OptimizationConstraintError
# ════════════════════════════════════════════════════════════════════════════

class TestOptimizationConstraintError:

    def test_stores_failed_predicates(self):
        err = OptimizationConstraintError(["gc", "cai"])
        assert err.failed_predicates == ["gc", "cai"]
        assert err.partial_result is None

    def test_with_partial_result(self):
        partial = {"sequence": "ATG", "cai": 0.5}
        err = OptimizationConstraintError(["gc"], partial_result=partial)
        assert err.partial_result == partial

    def test_message_mentions_strict_mode(self):
        err = OptimizationConstraintError(["gc", "splice"])
        msg = str(err)
        assert "2 predicate(s)" in msg
        assert "strict mode" in msg
        assert "strict_mode=False" in msg


# ════════════════════════════════════════════════════════════════════════════
# BiosecurityError
# ════════════════════════════════════════════════════════════════════════════

class TestBiosecurityError:

    def test_legacy_form(self):
        err = BiosecurityError(
            "toxin detected",
            risk_level="high",
            flagged_categories=["toxin"],
            matches=[{"gene": "ctx"}],
        )
        assert err.reason == "toxin detected"
        assert err.risk_level == "high"
        assert err.flagged_categories == ["toxin"]
        assert err.report is None

    def test_legacy_form_defaults(self):
        err = BiosecurityError("danger")
        assert err.risk_level is None
        assert err.flagged_categories == []
        assert err.matches == []

    def test_report_form(self):
        class FakeReport:
            is_hazardous = True
            risk_level = "critical"
            flagged_categories = ["virulence"]
            matches = [{"hit": "x"}]
        err = BiosecurityError(FakeReport())
        assert err.report is not None
        assert err.risk_level == "critical"

    def test_message_includes_risk_level(self):
        err = BiosecurityError("hazard", risk_level="high")
        msg = str(err)
        assert "high" in msg
        assert "BLOCKED" in msg or "blocked" in msg


# ════════════════════════════════════════════════════════════════════════════
# TranslationVerificationError
# ════════════════════════════════════════════════════════════════════════════

class TestTranslationVerificationError:

    def test_rich_form_with_mismatches(self):
        mismatches = [
            {"position": 5, "expected": "M", "actual": "K", "codon_used": "AAA"},
        ]
        err = TranslationVerificationError(
            mismatches=mismatches,
            has_premature_stop=False,
            has_stop_codon=True,
            length_correct=True,
            translated_protein="MKFLIL",
            expected_protein="MKFLVL",
            dna_sequence="ATGAAATTTTTAATTTAA",
        )
        assert len(err.mismatches) == 1
        assert err.has_stop_codon is True
        assert err.length_correct is True
        assert err.translated_protein == "MKFLIL"
        assert err.expected_protein == "MKFLVL"

    def test_simple_form_with_reason(self):
        err = TranslationVerificationError(
            reason="mismatch at position 3",
        )
        msg = str(err)
        assert "mismatch at position 3" in msg

    def test_premature_stop_in_message(self):
        err = TranslationVerificationError(
            has_premature_stop=True,
            translated_protein="MKF",
            expected_protein="MKFLVL",
        )
        msg = str(err)
        assert "premature stop" in msg.lower()

    def test_length_mismatch_in_message(self):
        err = TranslationVerificationError(
            length_correct=False,
            dna_sequence="ATGAAA",
            expected_protein="MKFLVL",
        )
        msg = str(err)
        assert "DNA length" in msg

    def test_no_mismatches_no_premature_stop(self):
        err = TranslationVerificationError(
            translated_protein="MKF",
            expected_protein="MKF",
        )
        assert err.mismatches == []
        assert err.has_premature_stop is False


# ════════════════════════════════════════════════════════════════════════════
# Inheritance chain tests
# ════════════════════════════════════════════════════════════════════════════

class TestInheritanceChains:
    """Verify that specific exceptions can be caught at the right level."""

    def test_esmfold_caught_as_engine_error(self):
        with pytest.raises(EngineError):
            raise ESMFoldError("fail")

    def test_foldx_caught_as_engine_error(self):
        with pytest.raises(EngineError):
            raise FoldXError("fail")

    def test_camsol_caught_as_engine_error(self):
        with pytest.raises(EngineError):
            raise CamSolError("fail")

    def test_immunogenicity_caught_as_engine_error(self):
        with pytest.raises(EngineError):
            raise ImmunogenicityError("fail")

    def test_engine_error_caught_as_biocompiler_error(self):
        with pytest.raises(BioCompilerError):
            raise EngineError("fail")

    def test_all_engine_errors_caught_by_engine_error(self):
        """All four engine-specific errors should be catchable as EngineError."""
        for exc_cls in [ESMFoldError, FoldXError, CamSolError, ImmunogenicityError]:
            with pytest.raises(EngineError):
                raise exc_cls("test")
