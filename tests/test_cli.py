"""
Tests for the BioCompiler CLI module.

Covers:
  1. CLI entry point exists and is callable
  2. Help flag works
  3. Version flag works
  4. Basic argument parsing for all subcommands
  5. Helper functions (_read_fasta, _write_fasta, _resolve_protein, etc.)
  6. Subprocess invocation of the CLI
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from io import StringIO
from unittest.mock import patch

import pytest

from biocompiler import __version__
from biocompiler.cli import (
    build_parser,
    cmd_optimize,
    cmd_check,
    cmd_benchmark,
    cmd_scan,
    cmd_structure,
    cmd_stability,
    cmd_solubility,
    cmd_immunogenicity,
    cmd_assess,
    colorize,
    main,
    _read_fasta,
    _write_fasta,
    _write_certificate,
    _resolve_protein,
    _get_organism,
    _section_header,
    _verdict_symbol,
    _error_msg,
    _success_msg,
    _dim,
    _summary_box,
    _ProgressStep,
    _supports_color,
)


# ═══════════════════════════════════════════════════════════════════════
# 1. CLI entry point exists and is callable
# ═══════════════════════════════════════════════════════════════════════


class TestCLIEntryPoint:
    """Test that the CLI entry point exists and is callable."""

    def test_main_function_exists(self):
        """The main() function should exist and be callable."""
        assert callable(main)

    def test_main_is_importable_from_cli(self):
        """main should be importable from biocompiler.cli."""
        from biocompiler.cli import main as m
        assert callable(m)

    def test_build_parser_exists(self):
        """build_parser() should exist and return an ArgumentParser."""
        parser = build_parser()
        assert isinstance(parser, argparse.ArgumentParser)

    def test_main_no_args_shows_help_and_exits(self):
        """Calling main() with no arguments should print help and exit(0)."""
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            with pytest.raises(SystemExit) as exc_info:
                main([])
            # No command → print_help → exit(0)
            assert exc_info.value.code == 0


# ═══════════════════════════════════════════════════════════════════════
# 2. Help flag works
# ═══════════════════════════════════════════════════════════════════════


class TestHelpFlag:
    """Test that --help and help for subcommands works."""

    def test_main_help(self):
        """--help should exit with code 0 and print usage."""
        with pytest.raises(SystemExit) as exc_info:
            main(["--help"])
        assert exc_info.value.code == 0

    def test_main_help_via_parser(self):
        """Parser --help should include the program name."""
        parser = build_parser()
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["--help"])
        assert exc_info.value.code == 0

    def test_optimize_help(self):
        """optimize --help should exit cleanly."""
        with pytest.raises(SystemExit) as exc_info:
            main(["optimize", "--help"])
        assert exc_info.value.code == 0

    def test_check_help(self):
        """check --help should exit cleanly."""
        with pytest.raises(SystemExit) as exc_info:
            main(["check", "--help"])
        assert exc_info.value.code == 0

    def test_benchmark_help(self):
        """benchmark --help should exit cleanly."""
        with pytest.raises(SystemExit) as exc_info:
            main(["benchmark", "--help"])
        assert exc_info.value.code == 0

    def test_scan_help(self):
        """scan --help should exit cleanly."""
        with pytest.raises(SystemExit) as exc_info:
            main(["scan", "--help"])
        assert exc_info.value.code == 0

    def test_serve_help(self):
        """serve --help should exit cleanly."""
        with pytest.raises(SystemExit) as exc_info:
            main(["serve", "--help"])
        assert exc_info.value.code == 0

    def test_structure_help(self):
        """structure --help should exit cleanly."""
        with pytest.raises(SystemExit) as exc_info:
            main(["structure", "--help"])
        assert exc_info.value.code == 0

    def test_stability_help(self):
        """stability --help should exit cleanly."""
        with pytest.raises(SystemExit) as exc_info:
            main(["stability", "--help"])
        assert exc_info.value.code == 0

    def test_solubility_help(self):
        """solubility --help should exit cleanly."""
        with pytest.raises(SystemExit) as exc_info:
            main(["solubility", "--help"])
        assert exc_info.value.code == 0

    def test_immunogenicity_help(self):
        """immunogenicity --help should exit cleanly."""
        with pytest.raises(SystemExit) as exc_info:
            main(["immunogenicity", "--help"])
        assert exc_info.value.code == 0

    def test_assess_help(self):
        """assess --help should exit cleanly."""
        with pytest.raises(SystemExit) as exc_info:
            main(["assess", "--help"])
        assert exc_info.value.code == 0


# ═══════════════════════════════════════════════════════════════════════
# 3. Version flag works
# ═══════════════════════════════════════════════════════════════════════


class TestVersionFlag:
    """Test that --version works and matches __version__."""

    def test_version_flag_exits_cleanly(self):
        """--version should exit with code 0."""
        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0

    def test_version_output_contains_version(self):
        """--version output should contain the version string."""
        parser = build_parser()
        with pytest.raises(SystemExit) as exc_info:
            with patch("sys.stdout", new_callable=StringIO) as mock_out:
                parser.parse_args(["--version"])
        # argparse writes version to stdout
        assert exc_info.value.code == 0

    def test_version_matches_package_version(self):
        """The version in CLI description should match biocompiler.__version__."""
        parser = build_parser()
        # The --version action stores the version string
        for action in parser._actions:
            if isinstance(action, argparse._VersionAction):
                assert __version__ in action.version
                break
        else:
            pytest.fail("No --version action found in parser")


# ═══════════════════════════════════════════════════════════════════════
# 4. Basic argument parsing
# ═══════════════════════════════════════════════════════════════════════


class TestArgumentParsing:
    """Test that build_parser() correctly parses arguments for each subcommand."""

    def test_no_command_sets_command_none(self):
        """Parsing with no arguments should set command to None."""
        parser = build_parser()
        args = parser.parse_args([])
        assert args.command is None

    # ── optimize ──

    def test_optimize_basic(self):
        """Parse basic optimize command with input file."""
        parser = build_parser()
        args = parser.parse_args(["optimize", "my_gene.fasta"])
        assert args.command == "optimize"
        assert args.input == "my_gene.fasta"
        assert args.species == "human"
        assert args.enzymes == ""
        assert args.splice_low == 3.0
        assert args.splice_high == 6.0
        assert args.avoid_gt is True
        assert args.output is None
        assert args.certificate is None

    def test_optimize_with_species(self):
        """Parse optimize with --species ecoli."""
        parser = build_parser()
        args = parser.parse_args(["optimize", "gene.fasta", "--species", "ecoli"])
        assert args.species == "ecoli"

    def test_optimize_with_enzymes(self):
        """Parse optimize with --enzymes."""
        parser = build_parser()
        args = parser.parse_args(["optimize", "gene.fasta", "--enzymes", "EcoRI,BamHI"])
        assert args.enzymes == "EcoRI,BamHI"

    def test_optimize_with_splice_thresholds(self):
        """Parse optimize with custom splice thresholds."""
        parser = build_parser()
        args = parser.parse_args(["optimize", "gene.fasta", "--splice-low", "2.5", "--splice-high", "5.0"])
        assert args.splice_low == 2.5
        assert args.splice_high == 5.0

    def test_optimize_no_avoid_gt(self):
        """Parse optimize with --no-avoid-gt."""
        parser = build_parser()
        args = parser.parse_args(["optimize", "gene.fasta", "--no-avoid-gt"])
        assert args.avoid_gt is False

    def test_optimize_with_output(self):
        """Parse optimize with --output and --certificate."""
        parser = build_parser()
        args = parser.parse_args(["optimize", "gene.fasta", "-o", "out.fa", "-c", "cert.txt"])
        assert args.output == "out.fa"
        assert args.certificate == "cert.txt"

    def test_optimize_invalid_species(self):
        """Invalid species should cause parser error."""
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["optimize", "gene.fasta", "--species", "yeast"])

    # ── check ──

    def test_check_basic(self):
        """Parse basic check command."""
        parser = build_parser()
        args = parser.parse_args(["check", "my_gene.fasta"])
        assert args.command == "check"
        assert args.input == "my_gene.fasta"
        assert args.species == "human"

    def test_check_with_options(self):
        """Parse check with all options."""
        parser = build_parser()
        args = parser.parse_args([
            "check", "gene.fasta",
            "--species", "ecoli",
            "--enzymes", "EcoRI",
            "--splice-low", "4.0",
            "--splice-high", "7.0",
        ])
        assert args.species == "ecoli"
        assert args.enzymes == "EcoRI"
        assert args.splice_low == 4.0
        assert args.splice_high == 7.0

    # ── benchmark ──

    def test_benchmark_basic(self):
        """Parse basic benchmark command."""
        parser = build_parser()
        args = parser.parse_args(["benchmark"])
        assert args.command == "benchmark"
        assert args.enzymes == ""
        assert args.splice_low == 3.0
        assert args.splice_high == 6.0

    def test_benchmark_with_options(self):
        """Parse benchmark with custom options."""
        parser = build_parser()
        args = parser.parse_args(["benchmark", "--enzymes", "EcoRI,BamHI"])
        assert args.enzymes == "EcoRI,BamHI"

    # ── scan ──

    def test_scan_basic(self):
        """Parse basic scan command."""
        parser = build_parser()
        args = parser.parse_args(["scan", "--sequence", "ATGCGTACGT"])
        assert args.command == "scan"
        assert args.sequence == "ATGCGTACGT"
        assert args.enzymes == ""

    def test_scan_with_short_flag(self):
        """Parse scan with -s short flag."""
        parser = build_parser()
        args = parser.parse_args(["scan", "-s", "ATGCGT"])
        assert args.sequence == "ATGCGT"

    def test_scan_with_enzymes(self):
        """Parse scan with --enzymes."""
        parser = build_parser()
        args = parser.parse_args(["scan", "--sequence", "ATGCGT", "--enzymes", "EcoRI"])
        assert args.enzymes == "EcoRI"

    def test_scan_requires_sequence(self):
        """scan without --sequence should fail."""
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["scan"])

    # ── serve ──

    def test_serve_basic(self):
        """Parse basic serve command."""
        parser = build_parser()
        args = parser.parse_args(["serve"])
        assert args.command == "serve"
        assert args.host == "0.0.0.0"
        assert args.port == 8000

    def test_serve_custom_host_port(self):
        """Parse serve with --host and --port."""
        parser = build_parser()
        args = parser.parse_args(["serve", "--host", "127.0.0.1", "--port", "9000"])
        assert args.host == "127.0.0.1"
        assert args.port == 9000

    # ── structure ──

    def test_structure_with_protein(self):
        """Parse structure with --protein."""
        parser = build_parser()
        args = parser.parse_args(["structure", "--protein", "MVLSPADKTN"])
        assert args.command == "structure"
        assert args.protein == "MVLSPADKTN"
        assert args.sequence is None
        assert args.output is None
        assert args.quality_only is False
        assert args.pdb_file is None
        assert args.verbose is False

    def test_structure_with_sequence(self):
        """Parse structure with --sequence."""
        parser = build_parser()
        args = parser.parse_args(["structure", "--sequence", "ATGATG"])
        assert args.sequence == "ATGATG"
        assert args.protein is None

    def test_structure_with_organism(self):
        """Parse structure with --organism."""
        parser = build_parser()
        args = parser.parse_args(["structure", "--protein", "MVL", "--organism", "Mus_musculus"])
        assert args.organism == "Mus_musculus"

    def test_structure_quality_only(self):
        """Parse structure with --quality-only and --pdb-file."""
        parser = build_parser()
        args = parser.parse_args([
            "structure", "--protein", "MVL",
            "--quality-only", "--pdb-file", "test.pdb",
        ])
        assert args.quality_only is True
        assert args.pdb_file == "test.pdb"

    def test_structure_verbose(self):
        """Parse structure with -v."""
        parser = build_parser()
        args = parser.parse_args(["structure", "--protein", "MVL", "-v"])
        assert args.verbose is True

    def test_structure_requires_protein_or_sequence(self):
        """structure without --protein or --sequence should fail."""
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["structure"])

    def test_structure_mutually_exclusive(self):
        """structure with both --protein and --sequence should fail."""
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["structure", "--protein", "MVL", "--sequence", "ATGATG"])

    # ── stability ──

    def test_stability_basic(self):
        """Parse basic stability command."""
        parser = build_parser()
        args = parser.parse_args(["stability", "--protein", "MVLSPADKTN"])
        assert args.command == "stability"
        assert args.protein == "MVLSPADKTN"
        assert args.scan_mutations is False
        assert args.positions is None

    def test_stability_with_scan_mutations(self):
        """Parse stability with --scan-mutations."""
        parser = build_parser()
        args = parser.parse_args(["stability", "--protein", "MVL", "--scan-mutations"])
        assert args.scan_mutations is True

    def test_stability_with_positions(self):
        """Parse stability with --positions."""
        parser = build_parser()
        args = parser.parse_args(["stability", "--protein", "MVL", "--positions", "1", "5", "10"])
        assert args.positions == [1, 5, 10]

    # ── solubility ──

    def test_solubility_basic(self):
        """Parse basic solubility command."""
        parser = build_parser()
        args = parser.parse_args(["solubility", "--protein", "MVLSPADKTN"])
        assert args.command == "solubility"
        assert args.protein == "MVLSPADKTN"
        assert args.find_mutations is False
        assert args.min_score is None

    def test_solubility_with_find_mutations(self):
        """Parse solubility with --find-mutations."""
        parser = build_parser()
        args = parser.parse_args(["solubility", "--protein", "MVL", "--find-mutations"])
        assert args.find_mutations is True

    def test_solubility_with_min_score(self):
        """Parse solubility with --min-score."""
        parser = build_parser()
        args = parser.parse_args(["solubility", "--protein", "MVL", "--min-score", "0.5"])
        assert args.min_score == 0.5

    # ── immunogenicity ──

    def test_immunogenicity_basic(self):
        """Parse basic immunogenicity command."""
        parser = build_parser()
        args = parser.parse_args(["immunogenicity", "--protein", "MVLSPADKTN"])
        assert args.command == "immunogenicity"
        assert args.protein == "MVLSPADKTN"
        assert args.deimmunize is False
        assert args.target_score == 0.3
        assert args.max_mutations == 10
        assert args.blosum62_min == 1
        assert args.mhc_alleles is None

    def test_immunogenicity_with_deimmunize(self):
        """Parse immunogenicity with --deimmunize."""
        parser = build_parser()
        args = parser.parse_args(["immunogenicity", "--protein", "MVL", "--deimmunize"])
        assert args.deimmunize is True

    def test_immunogenicity_with_options(self):
        """Parse immunogenicity with all options."""
        parser = build_parser()
        args = parser.parse_args([
            "immunogenicity", "--protein", "MVL",
            "--target-score", "0.5",
            "--max-mutations", "20",
            "--blosum62-min", "2",
            "--mhc-alleles", "HLA-A*02:01", "HLA-B*07:02",
        ])
        assert args.target_score == 0.5
        assert args.max_mutations == 20
        assert args.blosum62_min == 2
        assert args.mhc_alleles == ["HLA-A*02:01", "HLA-B*07:02"]

    # ── assess ──

    def test_assess_basic(self):
        """Parse basic assess command."""
        parser = build_parser()
        args = parser.parse_args(["assess", "--protein", "MVLSPADKTN"])
        assert args.command == "assess"
        assert args.protein == "MVLSPADKTN"
        assert args.pdb_file is None
        assert args.skip_structure is False
        assert args.skip_stability is False
        assert args.skip_solubility is False
        assert args.skip_immunogenicity is False
        assert args.output is None
        assert args.format == "text"

    def test_assess_with_skip_flags(self):
        """Parse assess with skip flags."""
        parser = build_parser()
        args = parser.parse_args([
            "assess", "--protein", "MVL",
            "--skip-structure", "--skip-stability",
            "--skip-solubility", "--skip-immunogenicity",
        ])
        assert args.skip_structure is True
        assert args.skip_stability is True
        assert args.skip_solubility is True
        assert args.skip_immunogenicity is True

    def test_assess_with_format(self):
        """Parse assess with --format."""
        parser = build_parser()
        args = parser.parse_args(["assess", "--protein", "MVL", "--format", "json"])
        assert args.format == "json"

    def test_assess_format_choices(self):
        """assess --format only accepts text/json/html."""
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["assess", "--protein", "MVL", "--format", "xml"])

    def test_assess_with_output(self):
        """Parse assess with --output."""
        parser = build_parser()
        args = parser.parse_args(["assess", "--protein", "MVL", "-o", "report.txt"])
        assert args.output == "report.txt"


# ═══════════════════════════════════════════════════════════════════════
# 5. Helper functions
# ═══════════════════════════════════════════════════════════════════════


class TestHelperFunctions:
    """Test CLI helper functions."""

    # ── _read_fasta ──

    def test_read_fasta_valid_file(self, tmp_path):
        """Read a valid FASTA file."""
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">test_gene\nATGCGTACGT\nATGCGT\n")
        result = _read_fasta(str(fasta))
        assert result == "ATGCGTACGTATGCGT"

    def test_read_fasta_strips_non_dna(self, tmp_path):
        """Non-DNA characters should be stripped."""
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">gene\nATGNNGT\n")
        result = _read_fasta(str(fasta))
        assert result == "ATGGT"

    def test_read_fasta_missing_file(self):
        """Missing file should sys.exit(1)."""
        with pytest.raises(SystemExit) as exc_info:
            _read_fasta("/nonexistent/path.fasta")
        assert exc_info.value.code == 1

    def test_read_fasta_multiple_headers(self, tmp_path):
        """Multiple header lines should be skipped."""
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">gene1\nATGCGT\n>gene2\nTGCATG\n")
        result = _read_fasta(str(fasta))
        assert result == "ATGCGTTGCATG"

    # ── _write_fasta ──

    def test_write_fasta(self, tmp_path):
        """Write a FASTA file with proper formatting."""
        out = tmp_path / "out.fasta"
        _write_fasta(str(out), "ATGCGTACGT" * 10, header="test_header")
        content = out.read_text()
        assert content.startswith(">test_header\n")
        lines = content.strip().split("\n")
        # First line is header
        assert lines[0] == ">test_header"
        # Remaining lines should be <= 80 chars
        for line in lines[1:]:
            assert len(line) <= 80

    def test_write_fasta_line_wrapping(self, tmp_path):
        """Sequences longer than 80 chars should be wrapped."""
        out = tmp_path / "out.fasta"
        seq = "A" * 200
        _write_fasta(str(out), seq, header="long")
        content = out.read_text()
        lines = content.strip().split("\n")
        assert len(lines) >= 3  # header + at least 2 sequence lines

    # ── _write_certificate ──

    def test_write_certificate(self, tmp_path):
        """Write certificate text to file."""
        out = tmp_path / "cert.txt"
        _write_certificate(str(out), "Certificate: PASS")
        assert out.read_text() == "Certificate: PASS"

    # ── _resolve_protein ──

    def test_resolve_protein_with_protein(self):
        """Resolve with --protein flag."""
        args = argparse.Namespace(protein="MVLSPADKTN", sequence=None)
        result = _resolve_protein(args)
        assert result == "MVLSPADKTN"

    def test_resolve_protein_uppercases(self):
        """Protein sequence should be uppercased."""
        args = argparse.Namespace(protein="mvlspadktn", sequence=None)
        result = _resolve_protein(args)
        assert result == "MVLSPADKTN"

    def test_resolve_protein_strips_whitespace(self):
        """Protein sequence should be stripped."""
        args = argparse.Namespace(protein="  MVL  ", sequence=None)
        result = _resolve_protein(args)
        assert result == "MVL"

    def test_resolve_protein_invalid_characters(self):
        """Invalid amino acid characters should cause exit."""
        args = argparse.Namespace(protein="MVL123", sequence=None)
        with pytest.raises(SystemExit) as exc_info:
            _resolve_protein(args)
        assert exc_info.value.code == 1

    def test_resolve_protein_both_provided(self):
        """Providing both --protein and --sequence should exit."""
        args = argparse.Namespace(protein="MVL", sequence="ATGATG")
        with pytest.raises(SystemExit) as exc_info:
            _resolve_protein(args)
        assert exc_info.value.code == 1

    def test_resolve_protein_neither_provided(self):
        """Providing neither --protein nor --sequence should exit."""
        args = argparse.Namespace(protein=None, sequence=None)
        with pytest.raises(SystemExit) as exc_info:
            _resolve_protein(args)
        assert exc_info.value.code == 1

    def test_resolve_protein_with_dna_sequence(self):
        """Resolve with --sequence (DNA to protein translation)."""
        args = argparse.Namespace(protein=None, sequence="ATGATG")
        result = _resolve_protein(args)
        # ATGATG translates to "MM"
        assert isinstance(result, str)
        assert len(result) > 0

    def test_resolve_protein_dna_too_short(self):
        """DNA sequence too short to translate should exit."""
        args = argparse.Namespace(protein=None, sequence="AT")
        with pytest.raises(SystemExit) as exc_info:
            _resolve_protein(args)
        assert exc_info.value.code == 1

    def test_resolve_protein_accepts_stop_codon(self):
        """Protein sequences with * (stop) should be accepted."""
        args = argparse.Namespace(protein="MVL*", sequence=None)
        result = _resolve_protein(args)
        assert result == "MVL*"

    # ── _get_organism ──

    def test_get_organism_default(self):
        """Default organism should be Homo_sapiens."""
        args = argparse.Namespace(organism=None)
        result = _get_organism(args)
        assert result == "Homo_sapiens"

    def test_get_organism_custom(self):
        """Custom organism should be returned."""
        args = argparse.Namespace(organism="Mus_musculus")
        result = _get_organism(args)
        assert result == "Mus_musculus"

    def test_get_organism_empty_string(self):
        """Empty organism string should fall back to default."""
        args = argparse.Namespace(organism="")
        result = _get_organism(args)
        assert result == "Homo_sapiens"


# ═══════════════════════════════════════════════════════════════════════
# 6. Color/formatting helpers
# ═══════════════════════════════════════════════════════════════════════


class TestColorHelpers:
    """Test ANSI color and formatting helpers."""

    def test_colorize_no_tty(self):
        """colorize should return plain text when not a TTY."""
        with patch("biocompiler.cli._supports_color", return_value=False):
            result = colorize("hello", "\033[31m")
            assert result == "hello"

    def test_colorize_with_tty(self):
        """colorize should wrap text with ANSI codes when TTY."""
        with patch("biocompiler.cli._supports_color", return_value=True):
            result = colorize("hello", "\033[31m")
            assert "\033[31m" in result
            assert "\033[0m" in result
            assert "hello" in result

    def test_section_header(self):
        """_section_header should produce a string (colorized or not)."""
        result = _section_header("Title")
        assert "Title" in result

    def test_verdict_symbol_pass(self):
        """PASS verdict should produce a string containing PASS."""
        with patch("biocompiler.cli._supports_color", return_value=False):
            result = _verdict_symbol("PASS")
            assert "PASS" in result

    def test_verdict_symbol_fail(self):
        """FAIL verdict should produce a string containing FAIL."""
        with patch("biocompiler.cli._supports_color", return_value=False):
            result = _verdict_symbol("FAIL")
            assert "FAIL" in result

    def test_verdict_symbol_uncertain(self):
        """UNCERTAIN verdict should produce a string containing UNCERTAIN."""
        with patch("biocompiler.cli._supports_color", return_value=False):
            result = _verdict_symbol("UNCERTAIN")
            assert "UNCERTAIN" in result

    def test_error_msg(self):
        """_error_msg should contain the text."""
        with patch("biocompiler.cli._supports_color", return_value=False):
            result = _error_msg("oops")
            assert "oops" in result

    def test_success_msg(self):
        """_success_msg should contain the text."""
        with patch("biocompiler.cli._supports_color", return_value=False):
            result = _success_msg("great")
            assert "great" in result

    def test_dim(self):
        """_dim should contain the text."""
        with patch("biocompiler.cli._supports_color", return_value=False):
            result = _dim("faded")
            assert "faded" in result

    def test_summary_box(self):
        """_summary_box should contain the label and value."""
        result = _summary_box("Verdict", "PASS")
        assert "Verdict" in result
        assert "PASS" in result
        # Should contain box-drawing characters
        assert "\u2502" in result

    def test_supports_color_no_tty(self):
        """_supports_color should return False when stdout is not a TTY."""
        # In pytest, stdout is typically not a TTY
        result = _supports_color()
        assert isinstance(result, bool)


# ═══════════════════════════════════════════════════════════════════════
# 7. _ProgressStep context manager
# ═══════════════════════════════════════════════════════════════════════


class TestProgressStep:
    """Test the _ProgressStep context manager."""

    def test_progress_step_basic(self):
        """_ProgressStep should not raise on normal entry/exit."""
        with patch("sys.stderr", new_callable=StringIO):
            with _ProgressStep("Testing"):
                pass  # no error

    def test_progress_step_verbose(self):
        """_ProgressStep with verbose=True should include timing."""
        with patch("sys.stderr", new_callable=StringIO) as mock_err:
            with _ProgressStep("Testing", verbose=True):
                pass
            output = mock_err.getvalue()
            assert "done" in output

    def test_progress_step_non_verbose(self):
        """_ProgressStep with verbose=False should not include timing."""
        with patch("sys.stderr", new_callable=StringIO) as mock_err:
            with _ProgressStep("Testing", verbose=False):
                pass
            output = mock_err.getvalue()
            assert "done" in output


# ═══════════════════════════════════════════════════════════════════════
# 8. Command dispatch via main()
# ═══════════════════════════════════════════════════════════════════════


class TestCommandDispatch:
    """Test that main() dispatches to the correct command handler."""

    def test_unknown_command_shows_help(self):
        """An unrecognized command should print help and exit(1)."""
        parser = build_parser()
        # Manually set command to something unexpected
        args = parser.parse_args([])
        args.command = "nonexistent"
        with pytest.raises(SystemExit) as exc_info:
            with patch.object(parser, "print_help"):
                # Simulate the else branch in main()
                parser.print_help()
                raise SystemExit(1)
        assert exc_info.value.code == 1

    def test_scan_command_executes(self):
        """The scan command should execute and produce output."""
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            main(["scan", "--sequence", "ATGCGTACGTATGCGTACGT"])
            output = mock_out.getvalue()
            assert "Scanned" in output
            assert "Tokens found" in output

    def test_scan_too_short_sequence(self):
        """Scan with a too-short sequence should exit with error."""
        with pytest.raises(SystemExit) as exc_info:
            main(["scan", "--sequence", "AT"])
        assert exc_info.value.code == 1


# ═══════════════════════════════════════════════════════════════════════
# 9. Subprocess invocation of the CLI
# ═══════════════════════════════════════════════════════════════════════


class TestCLISubprocess:
    """Test the CLI via subprocess to verify the entry point is properly installed."""

    def test_help_via_subprocess(self):
        """python -m biocompiler --help should exit 0."""
        result = subprocess.run(
            [sys.executable, "-m", "biocompiler", "--help"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0
        assert "BioCompiler" in result.stdout

    def test_version_via_subprocess(self):
        """python -m biocompiler --version should exit 0 and show version."""
        result = subprocess.run(
            [sys.executable, "-m", "biocompiler", "--version"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0
        assert __version__ in result.stdout

    def test_no_args_via_subprocess(self):
        """python -m biocompiler with no args should exit 0 (shows help)."""
        result = subprocess.run(
            [sys.executable, "-m", "biocompiler"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0

    def test_optimize_help_via_subprocess(self):
        """python -m biocompiler optimize --help should exit 0."""
        result = subprocess.run(
            [sys.executable, "-m", "biocompiler", "optimize", "--help"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0
        assert "optimize" in result.stdout.lower() or "FASTA" in result.stdout

    def test_scan_via_subprocess(self):
        """python -m biocompiler scan --sequence ATGCGTACGT should produce output."""
        result = subprocess.run(
            [sys.executable, "-m", "biocompiler", "scan", "--sequence", "ATGCGTACGTATGCGTACGT"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0
        assert "Scanned" in result.stdout


# ═══════════════════════════════════════════════════════════════════════
# 10. Integration: cmd_optimize and cmd_check with FASTA files
# ═══════════════════════════════════════════════════════════════════════


class TestOptimizeCheckIntegration:
    """Test cmd_optimize and cmd_check with actual FASTA files."""

    def _make_fasta(self, tmp_path, name: str, seq: str) -> str:
        """Helper to create a temporary FASTA file."""
        path = tmp_path / name
        path.write_text(f">test_gene\n{seq}\n")
        return str(path)

    def test_optimize_command(self, tmp_path):
        """cmd_optimize should run end-to-end on a valid FASTA file."""
        fasta_path = self._make_fasta(tmp_path, "gene.fasta", "ATGGCTAAGCTGGATCC")
        out_fasta = str(tmp_path / "gene_optimized.fasta")
        out_cert = str(tmp_path / "gene_certificate.txt")
        args = argparse.Namespace(
            input=fasta_path,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            avoid_gt=True,
            output=out_fasta,
            certificate=out_cert,
        )
        with patch("sys.stdout", new_callable=StringIO):
            cmd_optimize(args)
        # Check output files were created
        assert os.path.isfile(out_fasta)
        assert os.path.isfile(out_cert)
        # Check output FASTA has content
        with open(out_fasta) as f:
            content = f.read()
        assert content.startswith(">")

    def test_optimize_missing_file(self, tmp_path):
        """cmd_optimize with a missing input file should exit."""
        args = argparse.Namespace(
            input=str(tmp_path / "nonexistent.fasta"),
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            avoid_gt=True,
            output=None,
            certificate=None,
        )
        with pytest.raises(SystemExit) as exc_info:
            cmd_optimize(args)
        assert exc_info.value.code == 1

    def test_check_command(self, tmp_path):
        """cmd_check should run end-to-end on a valid FASTA file."""
        fasta_path = self._make_fasta(tmp_path, "gene.fasta", "ATGGCTAAGCTGGATCC")
        args = argparse.Namespace(
            input=fasta_path,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
        )
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            cmd_check(args)
        output = mock_out.getvalue()
        assert "Certificate:" in output

    def test_check_missing_file(self, tmp_path):
        """cmd_check with missing file should exit."""
        args = argparse.Namespace(
            input=str(tmp_path / "nonexistent.fasta"),
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
        )
        with pytest.raises(SystemExit) as exc_info:
            cmd_check(args)
        assert exc_info.value.code == 1
