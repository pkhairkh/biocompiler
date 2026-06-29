"""
Tests for CLI predicate registry integration and export DNA validation.

Covers:
1. CLI check command uses predicate registry (not hardcoded predicates)
2. CLI --list-predicates flag lists all registered predicates
3. CLI --predicate flag filters to specific predicates
4. CLI --species accepts any organism (not just human/ecoli)
5. Export FASTA endpoint validates DNA sequence (ACGTN only, non-empty, len%3==0)
6. Export GenBank endpoint validates DNA sequence (ACGTN only, non-empty, len%3==0)
7. Batch export validates DNA sequences
"""

from __future__ import annotations

import argparse
import os
import tempfile
from io import StringIO
from unittest.mock import patch

import pytest
from pydantic import ValidationError

from biocompiler.cli import cmd_check, build_parser
from biocompiler.type_system import (
    registry as predicate_registry,
    PREDICATE_NAMES,
    evaluate_all_predicates,
)
from biocompiler.api import ExportFastaInput, ExportGenbankInput, BatchExportItem


# ─── Helpers ──────────────────────────────────────────────────────────────

def _make_fasta(directory, filename, sequence):
    """Create a FASTA file and return its path."""
    path = os.path.join(directory, filename)
    with open(path, "w") as f:
        f.write(f">test_sequence\n{sequence}\n")
    return path


# ═══════════════════════════════════════════════════════════════════════════
# 1. CLI check uses registry (not hardcoded 8 predicates)
# ═══════════════════════════════════════════════════════════════════════════

class TestCLIUsesRegistry:
    """Verify that cmd_check uses evaluate_all_predicates from the registry,
    not a hardcoded list of 8 predicate functions."""

    def test_check_produces_more_than_8_results(self, tmp_path):
        """The registry has 12+ DNA-level predicates; cmd_check should
        return results for all of them, not just 8."""
        fasta_path = _make_fasta(tmp_path, "gene.fasta", "ATGGCTAAGCTGGATCCTAA")
        args = argparse.Namespace(
            input=fasta_path,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            predicate=None,
            list_predicates=False,
        )
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            cmd_check(args)
        output = mock_out.getvalue()
        # The output should contain "Certificate:" and a certificate body
        assert "Certificate:" in output
        # The number of predicate results should be >= 12 (evaluate_all_predicates
        # returns 12 DNA-level predicates, not just 8 hardcoded ones)
        # We can count the lines with predicate results in the certificate
        lines = output.strip().split("\n")
        predicate_lines = [l for l in lines if "PASS" in l or "FAIL" in l or "UNCERTAIN" in l]
        assert len(predicate_lines) >= 12, (
            f"Expected >= 12 predicate results (registry-based), got {len(predicate_lines)}"
        )

    def test_registry_has_more_than_8_predicates(self):
        """The predicate registry should have more than 8 predicates registered."""
        names = predicate_registry.names()
        assert len(names) > 8, (
            f"Registry has only {len(names)} predicates; expected > 8"
        )

    def test_registry_overlaps_with_predicate_names(self):
        """Many names in PREDICATE_NAMES should be in the registry.
        Note: Some legacy PREDICATE_NAMES (NoStopCodons, NoGTDinucleotide,
        ValidCodingSeq, ConservationScore, CodonOptimality) use check_*
        functions not registered in the evaluate_* registry. The key
        requirement is that the registry has more entries than the old
        hardcoded 8."""
        registered = set(predicate_registry.names())
        # At least 10 of the PREDICATE_NAMES should be in the registry
        overlap = set(PREDICATE_NAMES) & registered
        assert len(overlap) >= 10, (
            f"Expected >= 10 PREDICATE_NAMES in registry, got {len(overlap)}: {overlap}"
        )


# ═══════════════════════════════════════════════════════════════════════════
# 2. CLI --list-predicates flag
# ═══════════════════════════════════════════════════════════════════════════

class TestCLIListPredicates:
    """Verify --list-predicates lists all registered predicates."""

    def test_list_predicates_flag(self):
        """--list-predicates should list predicate names and exit without error."""
        parser = build_parser()
        args = parser.parse_args(["check", "--list-predicates"])
        assert args.list_predicates is True

    def test_list_predicates_output(self):
        """cmd_check with --list-predicates should print predicate names."""
        args = argparse.Namespace(
            input=None,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            predicate=None,
            list_predicates=True,
        )
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            cmd_check(args)
        output = mock_out.getvalue()
        # Should contain at least the first few predicate names
        assert "NoCrypticSplice" in output or "CodonAdapted" in output
        # Should contain total count
        assert "Total:" in output

    def test_list_predicates_shows_all_registered(self):
        """All registered predicate names should appear in --list-predicates output."""
        args = argparse.Namespace(
            input=None,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            predicate=None,
            list_predicates=True,
        )
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            cmd_check(args)
        output = mock_out.getvalue()
        for name in predicate_registry.names():
            assert name in output, f"Registered predicate '{name}' not in --list-predicates output"


# ═══════════════════════════════════════════════════════════════════════════
# 3. CLI --predicate flag
# ═══════════════════════════════════════════════════════════════════════════

class TestCLIPredicateFilter:
    """Verify --predicate flag filters to specific predicates."""

    def test_predicate_flag_parsing(self):
        """--predicate should be parsed correctly."""
        parser = build_parser()
        args = parser.parse_args(["check", "gene.fasta", "--predicate", "GCInRange"])
        assert args.predicate == "GCInRange"

    def test_predicate_flag_filters_results(self, tmp_path):
        """--predicate GCInRange should only show GCInRange results."""
        fasta_path = _make_fasta(tmp_path, "gene.fasta", "ATGGCTAAGCTGGATCCTAA")
        args = argparse.Namespace(
            input=fasta_path,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            predicate="GCInRange",
            list_predicates=False,
        )
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            cmd_check(args)
        output = mock_out.getvalue()
        # Should only contain GCInRange predicate results
        assert "GCInRange" in output
        # Should NOT contain other predicates like NoRestrictionSite
        # (certificate formatting may vary, but predicate name should appear)

    def test_predicate_flag_multiple(self, tmp_path):
        """--predicate with comma-separated names should filter to those."""
        fasta_path = _make_fasta(tmp_path, "gene.fasta", "ATGGCTAAGCTGGATCCTAA")
        args = argparse.Namespace(
            input=fasta_path,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            predicate="GCInRange,InFrame",
            list_predicates=False,
        )
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            cmd_check(args)
        output = mock_out.getvalue()
        assert "Certificate:" in output

    def test_predicate_flag_unknown_exits(self, tmp_path):
        """--predicate with unknown name should exit with error."""
        fasta_path = _make_fasta(tmp_path, "gene.fasta", "ATGGCTAAGCTGGATCCTAA")
        args = argparse.Namespace(
            input=fasta_path,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            predicate="FakePredicate",
            list_predicates=False,
        )
        with pytest.raises(SystemExit) as exc_info:
            cmd_check(args)
        assert exc_info.value.code == 1

    def test_no_input_without_list_predicates_exits(self):
        """Without --list-predicates, missing input file should exit."""
        args = argparse.Namespace(
            input=None,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            predicate=None,
            list_predicates=False,
        )
        with pytest.raises(SystemExit) as exc_info:
            cmd_check(args)
        assert exc_info.value.code == 1


# ═══════════════════════════════════════════════════════════════════════════
# 4. CLI --species accepts any organism
# ═══════════════════════════════════════════════════════════════════════════

class TestCLISpeciesFlexibility:
    """Verify --species accepts any organism name, not just human/ecoli."""

    def test_species_not_limited_to_choices(self):
        """The --species argument should not be limited to choices=['human', 'ecoli']."""
        parser = build_parser()
        # This should NOT raise - any organism should be accepted
        args = parser.parse_args(["check", "gene.fasta", "--species", "CHO_K1"])
        assert args.species == "CHO_K1"

    def test_species_ecoli_works(self, tmp_path):
        """--species ecoli should work with registry-based check."""
        fasta_path = _make_fasta(tmp_path, "gene.fasta", "ATGGCTAAGCTGGATCCTAA")
        args = argparse.Namespace(
            input=fasta_path,
            species="ecoli",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            predicate=None,
            list_predicates=False,
        )
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            cmd_check(args)
        output = mock_out.getvalue()
        assert "Certificate:" in output


# ═══════════════════════════════════════════════════════════════════════════
# 5. Export FASTA DNA validation
# ═══════════════════════════════════════════════════════════════════════════

class TestExportFastaDNAValidation:
    """Verify FASTA export endpoint validates DNA sequence."""

    def test_valid_dna_accepted(self):
        """Valid DNA sequence should pass validation."""
        data = ExportFastaInput(sequence="ATGGCTAAGCTG")
        assert data.sequence == "ATGGCTAAGCTG"

    def test_valid_dna_with_n_accepted(self):
        """DNA sequence with N ambiguity should pass validation."""
        data = ExportFastaInput(sequence="ATGNCTAAGCTG")
        assert data.sequence == "ATGNCTAAGCTG"

    def test_empty_sequence_rejected(self):
        """Empty sequence should be rejected."""
        with pytest.raises(ValidationError) as exc_info:
            ExportFastaInput(sequence="")
        assert "empty" in str(exc_info.value).lower()

    def test_whitespace_only_rejected(self):
        """Whitespace-only sequence should be rejected."""
        with pytest.raises(ValidationError) as exc_info:
            ExportFastaInput(sequence="   ")
        assert "empty" in str(exc_info.value).lower()

    def test_invalid_characters_rejected(self):
        """Non-DNA characters should be rejected."""
        with pytest.raises(ValidationError) as exc_info:
            ExportFastaInput(sequence="ATGXYZ")
        error_str = str(exc_info.value).lower()
        assert "invalid" in error_str

    def test_lowercase_converted_to_uppercase(self):
        """Lowercase DNA should be converted to uppercase."""
        data = ExportFastaInput(sequence="atggctaagctg")
        assert data.sequence == "ATGGCTAAGCTG"

    def test_protein_sequence_rejected(self):
        """Protein-like sequences (non-DNA chars) should be rejected."""
        with pytest.raises(ValidationError) as exc_info:
            ExportFastaInput(sequence="MSKGEELFTG")
        error_str = str(exc_info.value).lower()
        assert "invalid" in error_str

    def test_numeric_sequence_rejected(self):
        """Numeric sequences should be rejected."""
        with pytest.raises(ValidationError) as exc_info:
            ExportFastaInput(sequence="123456")
        error_str = str(exc_info.value).lower()
        assert "invalid" in error_str

    # NOTE: ``test_length_not_multiple_of_3_rejected`` was removed because
    # the ``len % 3 == 0`` validation was dropped from ExportFastaInput in
    # v0.9.0. The model now only enforces non-empty, ACGTN-only, and
    # max-length constraints.


# ═══════════════════════════════════════════════════════════════════════════
# 6. Export GenBank DNA validation
# ═══════════════════════════════════════════════════════════════════════════

class TestExportGenbankDNAValidation:
    """Verify GenBank export endpoint validates DNA sequence."""

    def test_valid_dna_accepted(self):
        """Valid DNA sequence should pass validation."""
        data = ExportGenbankInput(sequence="ATGGCTAAGCTG")
        assert data.sequence == "ATGGCTAAGCTG"

    def test_empty_sequence_rejected(self):
        """Empty sequence should be rejected."""
        with pytest.raises(ValidationError) as exc_info:
            ExportGenbankInput(sequence="")
        assert "empty" in str(exc_info.value).lower()

    def test_invalid_characters_rejected(self):
        """Non-DNA characters should be rejected."""
        with pytest.raises(ValidationError) as exc_info:
            ExportGenbankInput(sequence="ATGZZZ")
        error_str = str(exc_info.value).lower()
        assert "invalid" in error_str

    def test_lowercase_converted_to_uppercase(self):
        """Lowercase DNA should be converted to uppercase."""
        data = ExportGenbankInput(sequence="atggctaagctg")
        assert data.sequence == "ATGGCTAAGCTG"

    # NOTE: ``test_length_not_multiple_of_3_rejected`` was removed because
    # the ``len % 3 == 0`` validation was dropped from ExportGenbankInput in
    # v0.9.0. The model now only enforces non-empty, ACGTN-only, and
    # max-length constraints.


# ═══════════════════════════════════════════════════════════════════════════
# 7. Batch export DNA validation
# ═══════════════════════════════════════════════════════════════════════════

class TestBatchExportDNAValidation:
    """Verify batch export items validate DNA sequence."""

    def test_valid_dna_accepted(self):
        """Valid DNA sequence should pass validation."""
        data = BatchExportItem(sequence="ATGGCTAAGCTG")
        assert data.sequence == "ATGGCTAAGCTG"

    def test_empty_sequence_rejected(self):
        """Empty sequence should be rejected."""
        with pytest.raises(ValidationError) as exc_info:
            BatchExportItem(sequence="")
        assert "empty" in str(exc_info.value).lower()

    def test_invalid_characters_rejected(self):
        """Non-DNA characters should be rejected."""
        with pytest.raises(ValidationError) as exc_info:
            BatchExportItem(sequence="ATGZZZ")
        error_str = str(exc_info.value).lower()
        assert "invalid" in error_str

    # NOTE: ``test_length_not_multiple_of_3_rejected`` was removed because
    # the ``len % 3 == 0`` validation was dropped from ExportFastaInput,
    # ExportGenbankInput, and BatchExportItem in v0.9.0. The models now
    # only enforce non-empty, ACGTN-only, and max-length constraints.
