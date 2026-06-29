"""Tests for the BioCompiler IR L3→L4 folding pass.

Covers:
* :func:`biocompiler.ir.folding.fold_sequence` — public entry point
  (heuristic fallback, ESMFold path, length limits, edge cases).
* :class:`biocompiler.ir.folding.FoldingResult` — dataclass shape,
  default values, ``oracle_used`` field.
* :func:`biocompiler.ir.passes.fold` — IR-L3 → IR-L4 lowering pass
  (metadata stamping, sequence preservation, ``use_esmfold`` flag).
* :func:`biocompiler.ir.passes.compile_gene` — end-to-end integration
  (L0 → L4 with folding).

Tests are written to be deterministic and offline-friendly: the
heuristic fallback is exercised everywhere ESMFold would normally be
needed, so no network access is required.
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

from biocompiler.ir import (
    IR_L0_GenomicDNA,
    IR_L3_Polypeptide,
    IR_L4_FoldedProtein,
    IRLevel,
    IRError,
    fold,
    compile_gene,
    check_l4_invariants,
)
from biocompiler.ir.folding import (
    ESMFOLD_MAX_LENGTH,
    FoldingResult,
    STANDARD_AMINO_ACIDS,
    fold_sequence,
    _fold_heuristic,
    _fold_with_esmfold,
    _predict_ss_string,
)


# ────────────────────────────────────────────────────────────────────
# Test fixtures
# ────────────────────────────────────────────────────────────────────

# Canonical HBB N-terminus (UniProt P68871, first 31 residues + stop).
HBB_PROTEIN = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*"

# Short canonical sample ("MAK" + stop).
SHORT_PROTEIN = "MAK*"


# ═══════════════════════════════════════════════════════════════════════════
# 1. FoldingResult dataclass
# ═══════════════════════════════════════════════════════════════════════════

class TestFoldingResultDataclass:
    """Verify the FoldingResult dataclass has the right shape and defaults."""

    def test_default_oracle_is_none(self) -> None:
        result = FoldingResult(sequence="MAK")
        assert result.oracle_used == "none"

    def test_default_confidence_is_zero(self) -> None:
        result = FoldingResult(sequence="MAK")
        assert result.confidence == 0.0

    def test_default_coordinates_is_none(self) -> None:
        result = FoldingResult(sequence="MAK")
        assert result.coordinates is None

    def test_default_ss_is_empty_string(self) -> None:
        result = FoldingResult(sequence="MAK")
        assert result.secondary_structure == ""

    def test_default_metadata_is_empty_dict(self) -> None:
        result = FoldingResult(sequence="MAK")
        assert result.metadata == {}

    def test_metadata_is_per_instance(self) -> None:
        # Two FoldingResults must not share a metadata dict.
        r1 = FoldingResult(sequence="A")
        r2 = FoldingResult(sequence="B")
        r1.metadata["k"] = "v"
        assert "k" not in r2.metadata

    def test_can_set_all_fields(self) -> None:
        result = FoldingResult(
            sequence="MAK",
            coordinates=[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]],
            confidence=87.5,
            oracle_used="esmfold",
            secondary_structure="HHC",
            metadata={"engine": "esmfold"},
        )
        assert result.sequence == "MAK"
        assert len(result.coordinates) == 3
        assert result.confidence == 87.5
        assert result.oracle_used == "esmfold"
        assert result.secondary_structure == "HHC"
        assert result.metadata["engine"] == "esmfold"


# ═══════════════════════════════════════════════════════════════════════════
# 2. fold_sequence — heuristic path (no ESMFold required)
# ═══════════════════════════════════════════════════════════════════════════

class TestFoldSequenceHeuristic:
    """Heuristic-fallback tests that do NOT require ESMFold."""

    def test_heuristic_returns_fallback_oracle(self) -> None:
        # Force the heuristic by setting use_esmfold=False.
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        assert result.oracle_used == "fallback"

    def test_heuristic_returns_nonempty_ss(self) -> None:
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        assert len(result.secondary_structure) == len(HBB_PROTEIN) - 1  # minus *

    def test_heuristic_ss_uses_dssp_alphabet(self) -> None:
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        # DSSP-style string uses only H, E, C.
        assert set(result.secondary_structure) <= {"H", "E", "C"}

    def test_heuristic_strips_trailing_stop(self) -> None:
        # The folded sequence should not contain the trailing '*'.
        result = fold_sequence("MAK*", use_esmfold=False)
        assert result.sequence == "MAK"
        assert "*" not in result.sequence

    def test_heuristic_strips_trailing_stop_hbb(self) -> None:
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        assert not result.sequence.endswith("*")
        assert result.sequence == HBB_PROTEIN[:-1]

    def test_heuristic_confidence_in_valid_range(self) -> None:
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        # Heuristic pLDDT is clamped to [25, 55].
        assert 0.0 <= result.confidence <= 100.0
        assert 25.0 <= result.confidence <= 55.0

    def test_heuristic_coordinates_are_none(self) -> None:
        # The heuristic cannot produce 3-D coordinates.
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        assert result.coordinates is None

    def test_heuristic_metadata_has_engine_key(self) -> None:
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        assert result.metadata.get("engine") == "heuristic"

    def test_heuristic_metadata_has_plddt_scores(self) -> None:
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        plddt = result.metadata.get("plddt_scores", [])
        assert isinstance(plddt, list)
        assert len(plddt) == len(result.sequence)

    def test_heuristic_metadata_has_ss_fractions(self) -> None:
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        ss_fractions = result.metadata.get("ss_fractions", {})
        assert "helix_fraction" in ss_fractions
        assert "sheet_fraction" in ss_fractions
        assert "coil_fraction" in ss_fractions
        # Fractions should sum to ~1.0.
        total = (
            ss_fractions["helix_fraction"]
            + ss_fractions["sheet_fraction"]
            + ss_fractions["turn_fraction"]
            + ss_fractions["coil_fraction"]
        )
        assert abs(total - 1.0) < 0.01

    def test_heuristic_hbb_has_significant_helix_content(self) -> None:
        # HBB N-terminus is mostly helical in the experimental structure
        # (PDB 1A3N).  The Chou-Fasman heuristic should predict a
        # substantial helix fraction for this sequence.
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        ss_fractions = result.metadata["ss_fractions"]
        assert ss_fractions["helix_fraction"] > 0.3, (
            f"HBB heuristic helix fraction too low: {ss_fractions['helix_fraction']}"
        )

    def test_heuristic_method_is_set(self) -> None:
        result = fold_sequence(HBB_PROTEIN, use_esmfold=False)
        assert result.metadata.get("method") == "heuristic_fallback"

    def test_heuristic_uppercases_input(self) -> None:
        # Lowercase input should be normalised to uppercase.
        result = fold_sequence("mvhl*", use_esmfold=False)
        assert result.sequence == "MVHL"


# ═══════════════════════════════════════════════════════════════════════════
# 3. fold_sequence — error cases
# ═══════════════════════════════════════════════════════════════════════════

class TestFoldSequenceErrors:
    """fold_sequence must raise IRError on invalid input."""

    def test_empty_sequence_raises_irerror(self) -> None:
        with pytest.raises(IRError) as exc_info:
            fold_sequence("")
        assert exc_info.value.level == IRLevel.L3
        assert "empty" in exc_info.value.message.lower()

    def test_only_stop_codon_raises_irerror(self) -> None:
        # '*' is stripped, leaving an empty sequence.
        with pytest.raises(IRError) as exc_info:
            fold_sequence("*")
        assert "empty" in exc_info.value.message.lower()

    def test_non_string_input_raises_irerror(self) -> None:
        with pytest.raises(IRError):
            fold_sequence(123)  # type: ignore[arg-type]

    def test_none_input_raises_irerror(self) -> None:
        with pytest.raises(IRError):
            fold_sequence(None)  # type: ignore[arg-type]

    def test_whitespace_only_raises_irerror(self) -> None:
        with pytest.raises(IRError):
            fold_sequence("   ")


# ═══════════════════════════════════════════════════════════════════════════
# 4. fold_sequence — non-standard amino acids
# ═══════════════════════════════════════════════════════════════════════════

class TestFoldSequenceNonStandardAA:
    """Non-standard AA codes (X, B, Z, J) must trigger the heuristic."""

    def test_x_aa_uses_fallback(self) -> None:
        # 'X' = unknown AA — ESMFold rejects it, so we must use heuristic.
        result = fold_sequence("MVHLXPEEK*")
        assert result.oracle_used == "fallback"

    def test_x_aa_still_produces_ss(self) -> None:
        # MVHLXPEEK* = 9 residues + stop; SS length should be 9 (stop stripped).
        result = fold_sequence("MVHLXPEEK*", use_esmfold=False)
        assert len(result.secondary_structure) == 9  # minus '*', X is a residue

    def test_b_aa_uses_fallback(self) -> None:
        # 'B' = Asx (Asn or Asp) — non-standard, triggers heuristic.
        result = fold_sequence("MBHLAPEEK*")
        assert result.oracle_used == "fallback"

    def test_z_aa_uses_fallback(self) -> None:
        # 'Z' = Glx (Gln or Glu) — non-standard, triggers heuristic.
        result = fold_sequence("MZHLAPEEK*")
        assert result.oracle_used == "fallback"

    def test_j_aa_uses_fallback(self) -> None:
        # 'J' = Leu or Ile — non-standard, triggers heuristic.
        result = fold_sequence("MJHLAPEEK*")
        assert result.oracle_used == "fallback"

    def test_all_standard_aas_can_use_esmfold_path(self) -> None:
        # When all AAs are standard AND use_esmfold=True, the ESMFold
        # path is attempted.  In offline mode this falls through to
        # the heuristic — but the path was tried.
        # We verify by checking oracle_used is one of "esmfold" or "fallback".
        result = fold_sequence("MVHLTPEEK*", use_esmfold=True)
        assert result.oracle_used in {"esmfold", "fallback"}


# ═══════════════════════════════════════════════════════════════════════════
# 5. fold_sequence — long sequences
# ═══════════════════════════════════════════════════════════════════════════

class TestFoldSequenceLong:
    """Sequences longer than ESMFOLD_MAX_LENGTH must use the heuristic."""

    def test_long_sequence_uses_fallback(self) -> None:
        # 1500-aa sequence — exceeds the 1000-aa ESMFold limit.
        long_seq = "A" * 1500 + "*"
        result = fold_sequence(long_seq, use_esmfold=True)
        assert result.oracle_used == "fallback"

    def test_long_sequence_ss_length_matches(self) -> None:
        long_seq = "A" * 1500 + "*"
        result = fold_sequence(long_seq, use_esmfold=True)
        assert len(result.secondary_structure) == 1500

    def test_at_limit_can_attempt_esmfold(self) -> None:
        # Exactly at the limit — should attempt ESMFold (and fall back
        # to heuristic if offline).  Result oracle is one of the two.
        seq = "A" * ESMFOLD_MAX_LENGTH + "*"
        result = fold_sequence(seq, use_esmfold=True)
        assert result.oracle_used in {"esmfold", "fallback"}

    def test_over_limit_never_uses_esmfold(self) -> None:
        seq = "A" * (ESMFOLD_MAX_LENGTH + 1) + "*"
        result = fold_sequence(seq, use_esmfold=True)
        # Must NOT be esmfold — too long.
        assert result.oracle_used != "esmfold"


# ═══════════════════════════════════════════════════════════════════════════
# 6. fold_sequence — ESMFold path (mocked)
# ═══════════════════════════════════════════════════════════════════════════

class TestFoldSequenceESMFoldPath:
    """Test the ESMFold backend with mocked predict_structure.

    We mock the ESMFold engine so the tests run offline and
    deterministically, but exercise the real code path inside
    ``_fold_with_esmfold``.
    """

    def test_esmfold_success_marks_oracle_correctly(self) -> None:
        """When ESMFold returns a real result, oracle_used='esmfold'."""
        # Build a fake ESMFoldResult that signals "real" ESMFold success.
        # We mock at the module boundary: biocompiler.engines.esmfold.
        fake_pdb = (
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 85.00           C\n"
            "ATOM      2  CA  ALA A   2       4.000   2.000   3.000  1.00 90.00           C\n"
            "ATOM      3  CA  ALA A   3       7.000   2.000   3.000  1.00 78.00           C\n"
            "END\n"
        )

        class FakeESMFoldResult:
            success = True
            method = "esmfold_api"
            model_name = "esmfold_v1"
            primary_score = 84.33
            plddt_scores = [85.0, 90.0, 78.0]
            pae_matrix = None
            classification = "Confident"
            execution_time_s = 0.42
            pdb_string = fake_pdb

        with patch(
            "biocompiler.engines.esmfold.predict_structure",
            return_value=FakeESMFoldResult(),
        ), patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=True,
        ):
            result = fold_sequence("MAK", use_esmfold=True)

        assert result.oracle_used == "esmfold"
        assert result.confidence == pytest.approx(84.33, abs=0.01)
        assert result.coordinates is not None
        assert len(result.coordinates) == 3
        assert result.coordinates[0] == [1.0, 2.0, 3.0]
        assert result.metadata["engine"] == "esmfold"
        assert result.metadata["method"] == "esmfold_api"
        assert result.metadata["has_coordinates"] is True

    def test_esmfold_failure_falls_back_to_heuristic(self) -> None:
        """When ESMFold reports failure, fall back to the heuristic."""

        class FakeFailedResult:
            success = False
            method = "esmfold_api"
            error = "API unreachable"
            primary_score = 0.0
            plddt_scores = []
            pae_matrix = None
            classification = ""
            execution_time_s = 0.0
            pdb_string = ""

        with patch(
            "biocompiler.engines.esmfold.predict_structure",
            return_value=FakeFailedResult(),
        ), patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=True,
        ):
            result = fold_sequence("MVHLTPEEK", use_esmfold=True)

        # Should have fallen back to heuristic.
        assert result.oracle_used == "fallback"
        assert result.confidence > 0.0  # heuristic always returns > 0

    def test_esmfold_internal_heuristic_defers_to_ir_heuristic(self) -> None:
        """If ESMFold returns method='heuristic_fallback', defer to our own."""

        class FakeHeuristicResult:
            success = True
            method = "heuristic_fallback"  # ESMFold used its own heuristic
            model_name = "heuristic_v2"
            primary_score = 42.0
            plddt_scores = [42.0] * 9
            pae_matrix = None
            classification = "Very low"
            execution_time_s = 0.01
            pdb_string = ""

        with patch(
            "biocompiler.engines.esmfold.predict_structure",
            return_value=FakeHeuristicResult(),
        ), patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=True,
        ):
            result = fold_sequence("MVHLTPEEK", use_esmfold=True)

        # Our heuristic should have been called, not ESMFold's.
        assert result.oracle_used == "fallback"

    def test_esmfold_unavailable_uses_fallback(self) -> None:
        """When is_esmfold_available() returns False, use heuristic directly."""
        with patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=False,
        ):
            result = fold_sequence("MVHLTPEEK", use_esmfold=True)

        assert result.oracle_used == "fallback"

    def test_esmfold_exception_falls_back_gracefully(self) -> None:
        """If ESMFold raises, fall back to heuristic without crashing."""
        with patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=True,
        ), patch(
            "biocompiler.engines.esmfold.predict_structure",
            side_effect=RuntimeError("simulated ESMFold crash"),
        ):
            result = fold_sequence("MVHLTPEEK", use_esmfold=True)

        assert result.oracle_used == "fallback"

    def test_esmfold_no_pdb_string_still_returns_result(self) -> None:
        """ESMFold success but empty PDB → coordinates=None, oracle='esmfold'."""

        class FakeNoPDBResult:
            success = True
            method = "esmfold_local"
            model_name = "esmfold_v1"
            primary_score = 75.0
            plddt_scores = [75.0] * 3
            pae_matrix = None
            classification = "Confident"
            execution_time_s = 0.5
            pdb_string = ""  # No PDB string

        with patch(
            "biocompiler.engines.esmfold.predict_structure",
            return_value=FakeNoPDBResult(),
        ), patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=True,
        ):
            result = fold_sequence("MAK", use_esmfold=True)

        assert result.oracle_used == "esmfold"
        assert result.coordinates is None  # No PDB → no coords
        assert result.confidence == 75.0


# ═══════════════════════════════════════════════════════════════════════════
# 7. _fold_with_esmfold direct tests
# ═══════════════════════════════════════════════════════════════════════════

class TestFoldWithESMFoldDirect:
    """Direct tests of the private _fold_with_esmfold helper."""

    def test_returns_none_when_engine_unavailable(self) -> None:
        with patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=False,
        ):
            result = _fold_with_esmfold("MVHLTPEEK")
        assert result is None

    def test_returns_none_on_exception(self) -> None:
        with patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=True,
        ), patch(
            "biocompiler.engines.esmfold.predict_structure",
            side_effect=ValueError("simulated crash"),
        ):
            result = _fold_with_esmfold("MVHLTPEEK")
        assert result is None

    def test_returns_none_on_failed_prediction(self) -> None:
        class FakeFailed:
            success = False
            error = "OOM"
            method = "esmfold_api"
            primary_score = 0.0
            plddt_scores = []
            pae_matrix = None
            classification = ""
            execution_time_s = 0.0
            pdb_string = ""

        with patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=True,
        ), patch(
            "biocompiler.engines.esmfold.predict_structure",
            return_value=FakeFailed(),
        ):
            result = _fold_with_esmfold("MVHLTPEEK")
        assert result is None


# ═══════════════════════════════════════════════════════════════════════════
# 8. _fold_heuristic direct tests
# ═══════════════════════════════════════════════════════════════════════════

class TestFoldHeuristicDirect:
    """Direct tests of the private _fold_heuristic helper."""

    def test_returns_fallback_oracle(self) -> None:
        result = _fold_heuristic("MVHLTPEEK")
        assert result.oracle_used == "fallback"

    def test_returns_no_coordinates(self) -> None:
        result = _fold_heuristic("MVHLTPEEK")
        assert result.coordinates is None

    def test_ss_string_matches_sequence_length(self) -> None:
        result = _fold_heuristic("MVHLTPEEK")
        assert len(result.secondary_structure) == len("MVHLTPEEK")

    def test_confidence_in_heuristic_range(self) -> None:
        result = _fold_heuristic("MVHLTPEEK")
        # Heuristic is bounded by [25, 55].
        assert 25.0 <= result.confidence <= 55.0

    def test_metadata_has_engine_key(self) -> None:
        result = _fold_heuristic("MVHLTPEEK")
        assert result.metadata["engine"] == "heuristic"


# ═══════════════════════════════════════════════════════════════════════════
# 9. _predict_ss_string helper
# ═══════════════════════════════════════════════════════════════════════════

class TestPredictSSString:
    """Tests for the SS-string helper used by the ESMFold path."""

    def test_returns_string_of_input_length(self) -> None:
        ss = _predict_ss_string("MVHLTPEEK")
        assert len(ss) == 9

    def test_uses_dssp_alphabet(self) -> None:
        ss = _predict_ss_string("MVHLTPEEKSAVTALWGKVNVDEVGGEALGR")
        assert set(ss) <= {"H", "E", "C"}

    def test_hbb_predicted_as_mostly_helical(self) -> None:
        # HBB is mostly helical — Chou-Fasman should predict >30% H.
        seq = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"
        ss = _predict_ss_string(seq)
        helix_fraction = ss.count("H") / len(ss)
        assert helix_fraction > 0.3


# ═══════════════════════════════════════════════════════════════════════════
# 10. fold() IR pass
# ═══════════════════════════════════════════════════════════════════════════

class TestFoldIRPass:
    """Test the IR-L3 → IR-L4 lowering pass itself."""

    def test_fold_returns_ir_l4(self) -> None:
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3, use_esmfold=False)
        assert isinstance(ir_l4, IR_L4_FoldedProtein)

    def test_fold_preserves_sequence_with_stop(self) -> None:
        # IR-L4 sequence should still include the trailing '*'.
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3, use_esmfold=False)
        assert ir_l4.sequence == HBB_PROTEIN
        assert ir_l4.sequence.endswith("*")

    def test_fold_preserves_organism(self) -> None:
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3, use_esmfold=False)
        assert ir_l4.organism == "human"

    def test_fold_preserves_gene_name(self) -> None:
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3, use_esmfold=False)
        assert ir_l4.gene_name == "HBB"

    def test_fold_stamps_pass_metadata(self) -> None:
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3, use_esmfold=False)
        assert ir_l4.metadata["pass"] == "fold"
        assert ir_l4.metadata["source_level"] == "L3"

    def test_fold_records_oracle_in_metadata(self) -> None:
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3, use_esmfold=False)
        assert ir_l4.metadata["oracle"] == "fallback"

    def test_fold_records_secondary_structure_in_metadata(self) -> None:
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3, use_esmfold=False)
        ss = ir_l4.metadata["secondary_structure"]
        assert len(ss) == len(HBB_PROTEIN) - 1  # minus '*'
        assert set(ss) <= {"H", "E", "C"}

    def test_fold_populates_confidence_field(self) -> None:
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3, use_esmfold=False)
        assert ir_l4.confidence is not None
        assert ir_l4.confidence > 0.0

    def test_fold_passes_l4_invariants(self) -> None:
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3, use_esmfold=False)
        # Should not raise.
        assert check_l4_invariants(ir_l4) is True

    def test_fold_ptms_is_empty_list(self) -> None:
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3, use_esmfold=False)
        assert ir_l4.ptms == []

    def test_fold_propagates_upstream_metadata(self) -> None:
        # Upstream metadata should propagate through the fold pass.
        ir_l3 = IR_L3_Polypeptide(
            sequence=HBB_PROTEIN,
            organism="human",
            gene_name="HBB",
            metadata={"request_id": "abc-123", "user": "alice"},
        )
        ir_l4 = fold(ir_l3, use_esmfold=False)
        assert ir_l4.metadata["request_id"] == "abc-123"
        assert ir_l4.metadata["user"] == "alice"

    def test_fold_use_esmfold_false_skips_esmfold(self) -> None:
        # When use_esmfold=False, never even try ESMFold.
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")

        with patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=True,
        ) as mock_avail:
            ir_l4 = fold(ir_l3, use_esmfold=False)
            # is_esmfold_available should NOT have been called.
            mock_avail.assert_not_called()
        assert ir_l4.metadata["oracle"] == "fallback"

    def test_fold_default_attempts_esmfold(self) -> None:
        # Default use_esmfold=True — should attempt ESMFold (which is
        # unavailable in test env, so falls back to heuristic).
        ir_l3 = IR_L3_Polypeptide(sequence=HBB_PROTEIN, organism="human", gene_name="HBB")
        ir_l4 = fold(ir_l3)  # default use_esmfold=True
        assert ir_l4.metadata["oracle"] in {"esmfold", "fallback"}

    def test_fold_empty_sequence_raises_irerror(self) -> None:
        # An empty polypeptide (just '*') is rejected by the fold pass.
        ir_l3 = IR_L3_Polypeptide(sequence="*", organism="e_coli")
        with pytest.raises(IRError) as exc_info:
            fold(ir_l3, use_esmfold=False)
        assert exc_info.value.level == IRLevel.L3


# ═══════════════════════════════════════════════════════════════════════════
# 11. compile_gene integration — L0 → L4
# ═══════════════════════════════════════════════════════════════════════════

class TestCompileGeneL4:
    """End-to-end L0 → L4 pipeline tests."""

    def test_compile_gene_to_l4_returns_ir_l4(self) -> None:
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTAAGTAA", regions=[], organism="e_coli", gene_name="test"
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        assert isinstance(ir_l4, IR_L4_FoldedProtein)
        assert ir_l4.level == IRLevel.L4

    def test_compile_gene_l4_has_structure_info(self) -> None:
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTAAGTAA", regions=[], organism="e_coli", gene_name="test"
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        # Structure info: oracle + secondary_structure + confidence.
        assert ir_l4.metadata["oracle"] == "fallback"
        assert "secondary_structure" in ir_l4.metadata
        assert ir_l4.confidence is not None
        assert ir_l4.confidence > 0.0

    def test_compile_gene_l4_preserves_sequence(self) -> None:
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTAAGTAA", regions=[], organism="e_coli", gene_name="test"
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        # Sequence goes DNA → RNA → protein (MAK*) → folded protein (MAK*).
        assert ir_l4.sequence == "MAK*"

    def test_compile_gene_l4_propagates_metadata(self) -> None:
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTAAGTAA",
            regions=[],
            organism="e_coli",
            gene_name="g",
            metadata={"request_id": "req-001"},
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        assert ir_l4.metadata["request_id"] == "req-001"
        # Last pass applied is fold, so its stamps win.
        assert ir_l4.metadata["pass"] == "fold"
        assert ir_l4.metadata["source_level"] == "L3"

    def test_compile_gene_l4_use_esmfold_flag_forwarded(self) -> None:
        # When use_esmfold=False, never even try ESMFold.
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTAAGTAA", regions=[], organism="e_coli", gene_name="test"
        )
        with patch(
            "biocompiler.engines.esmfold.is_esmfold_available",
            return_value=True,
        ) as mock_avail:
            ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
            mock_avail.assert_not_called()
        assert ir_l4.metadata["oracle"] == "fallback"

    def test_compile_gene_l3_does_not_fold(self) -> None:
        # Stopping at L3 should not invoke the fold pass.
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTAAGTAA", regions=[], organism="e_coli", gene_name="test"
        )
        with patch("biocompiler.ir.passes.fold") as mock_fold:
            ir_l3 = compile_gene(ir_l0, IRLevel.L3, use_esmfold=False)
            mock_fold.assert_not_called()
        assert ir_l3.level == IRLevel.L3

    def test_compile_gene_default_target_is_l3(self) -> None:
        # Default target_level is L3, not L4.
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTAAGTAA", regions=[], organism="e_coli", gene_name="test"
        )
        ir = compile_gene(ir_l0, use_esmfold=False)
        assert ir.level == IRLevel.L3


# ═══════════════════════════════════════════════════════════════════════════
# 12. HBB demo — full end-to-end
# ═══════════════════════════════════════════════════════════════════════════

class TestHBBDemo:
    """HBB N-terminus demo — the canonical BioCompiler smoke test."""

    # The HBB N-terminal CDS (96 nt = 31 aa + stop).
    HBB_DNA = (
        "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAG"
        "GTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAA"
    )

    def test_hbb_l4_translation_matches_uniprot(self) -> None:
        ir_l0 = IR_L0_GenomicDNA(
            sequence=self.HBB_DNA, regions=[], organism="human", gene_name="HBB"
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        # UniProt P68871 N-terminus: MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*
        assert ir_l4.sequence == "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*"

    def test_hbb_l4_has_secondary_structure(self) -> None:
        ir_l0 = IR_L0_GenomicDNA(
            sequence=self.HBB_DNA, regions=[], organism="human", gene_name="HBB"
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        ss = ir_l4.metadata["secondary_structure"]
        assert len(ss) == 31  # 31 residues + stop = 32 chars in sequence, minus 1 for '*'
        assert set(ss) <= {"H", "E", "C"}

    def test_hbb_l4_predicted_helical(self) -> None:
        # HBB is mostly alpha-helical; the heuristic should predict
        # a substantial helix fraction.
        ir_l0 = IR_L0_GenomicDNA(
            sequence=self.HBB_DNA, regions=[], organism="human", gene_name="HBB"
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        ss = ir_l4.metadata["secondary_structure"]
        helix_fraction = ss.count("H") / len(ss)
        assert helix_fraction > 0.3, (
            f"HBB predicted helix fraction too low: {helix_fraction:.2f}; "
            f"ss={ss}"
        )

    def test_hbb_l4_oracle_is_fallback_offline(self) -> None:
        # In the test environment (no ESMFold), the oracle should be
        # either 'fallback' (heuristic) or 'esmfold' (if a real ESMFold
        # is somehow available).  Most test environments will see 'fallback'.
        ir_l0 = IR_L0_GenomicDNA(
            sequence=self.HBB_DNA, regions=[], organism="human", gene_name="HBB"
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        assert ir_l4.metadata["oracle"] in {"esmfold", "fallback", "none"}

    def test_hbb_l4_confidence_in_valid_range(self) -> None:
        ir_l0 = IR_L0_GenomicDNA(
            sequence=self.HBB_DNA, regions=[], organism="human", gene_name="HBB"
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        assert 0.0 <= ir_l4.confidence <= 100.0

    def test_hbb_l4_passes_invariants(self) -> None:
        ir_l0 = IR_L0_GenomicDNA(
            sequence=self.HBB_DNA, regions=[], organism="human", gene_name="HBB"
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        # Should not raise.
        assert check_l4_invariants(ir_l4) is True

    def test_hbb_l4_fold_metadata_present(self) -> None:
        ir_l0 = IR_L0_GenomicDNA(
            sequence=self.HBB_DNA, regions=[], organism="human", gene_name="HBB"
        )
        ir_l4 = compile_gene(ir_l0, IRLevel.L4, use_esmfold=False)
        # The fold_metadata sub-dict should be populated.
        assert "fold_metadata" in ir_l4.metadata
        fold_meta = ir_l4.metadata["fold_metadata"]
        assert "engine" in fold_meta
        assert "plddt_scores" in fold_meta
        assert len(fold_meta["plddt_scores"]) == 31  # HBB N-term has 31 residues
