"""Tests for organism-aware constraint selection in the API and CLI (F1.4).

Validates that:
1. Auto mode correctly detects E. coli as prokaryote
2. Auto mode correctly detects human as eukaryote
3. Explicit domain overrides auto-detection
4. CLI option parsing works correctly
5. The resolve_organism_domain helper handles all cases
6. The ProteinInput model validates organism_domain correctly
7. The OptimizeResponse includes the resolved domain
"""

from __future__ import annotations

import argparse

import pytest

from biocompiler.api import (
    ProteinInput,
    OptimizeResponse,
    BatchOptimizeItem,
    resolve_organism_domain,
)
from biocompiler.cli import build_parser


# ─── resolve_organism_domain helper ─────────────────────────────────


class TestResolveOrganismDomain:
    """Tests for the resolve_organism_domain helper function."""

    def test_auto_ecoli_returns_prokaryote(self) -> None:
        """Auto mode should detect E. coli as prokaryote."""
        result = resolve_organism_domain("Escherichia_coli", "auto")
        assert result == "prokaryote"

    def test_auto_ecoli_alias_returns_prokaryote(self) -> None:
        """Auto mode should detect 'ecoli' alias as prokaryote."""
        result = resolve_organism_domain("ecoli", "auto")
        assert result == "prokaryote"

    def test_auto_ecoli_short_alias_returns_prokaryote(self) -> None:
        """Auto mode should detect 'E_coli' alias as prokaryote."""
        result = resolve_organism_domain("E_coli", "auto")
        assert result == "prokaryote"

    def test_auto_human_returns_eukaryote(self) -> None:
        """Auto mode should detect human as eukaryote."""
        result = resolve_organism_domain("Homo_sapiens", "auto")
        assert result == "eukaryote"

    def test_auto_human_alias_returns_eukaryote(self) -> None:
        """Auto mode should detect 'human' alias as eukaryote."""
        result = resolve_organism_domain("human", "auto")
        assert result == "eukaryote"

    def test_auto_mouse_returns_eukaryote(self) -> None:
        """Auto mode should detect mouse as eukaryote."""
        result = resolve_organism_domain("Mus_musculus", "auto")
        assert result == "eukaryote"

    def test_auto_yeast_returns_eukaryote(self) -> None:
        """Auto mode should detect yeast as eukaryote."""
        result = resolve_organism_domain("Saccharomyces_cerevisiae", "auto")
        assert result == "eukaryote"

    def test_auto_cho_returns_eukaryote(self) -> None:
        """Auto mode should detect CHO as eukaryote."""
        result = resolve_organism_domain("CHO_K1", "auto")
        assert result == "eukaryote"

    def test_auto_unknown_defaults_to_eukaryote(self) -> None:
        """Auto mode should default to eukaryote for unknown organisms (safe default)."""
        result = resolve_organism_domain("Unknown_organism_42", "auto")
        assert result == "eukaryote"


class TestExplicitDomainOverride:
    """Tests that explicit domain settings override auto-detection."""

    def test_explicit_eukaryote_for_ecoli(self) -> None:
        """Setting 'eukaryote' for E. coli should override auto-detection."""
        result = resolve_organism_domain("Escherichia_coli", "eukaryote")
        assert result == "eukaryote"

    def test_explicit_prokaryote_for_human(self) -> None:
        """Setting 'prokaryote' for human should override auto-detection."""
        result = resolve_organism_domain("Homo_sapiens", "prokaryote")
        assert result == "prokaryote"

    def test_explicit_eukaryote_for_human(self) -> None:
        """Setting 'eukaryote' for human should confirm auto-detection."""
        result = resolve_organism_domain("Homo_sapiens", "eukaryote")
        assert result == "eukaryote"

    def test_explicit_prokaryote_for_ecoli(self) -> None:
        """Setting 'prokaryote' for E. coli should confirm auto-detection."""
        result = resolve_organism_domain("Escherichia_coli", "prokaryote")
        assert result == "prokaryote"


class TestResolveOrganismDomainValidation:
    """Tests for invalid inputs to resolve_organism_domain."""

    def test_invalid_domain_raises_value_error(self) -> None:
        """An invalid organism_domain value should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid organism_domain"):
            resolve_organism_domain("Homo_sapiens", "archaea")

    def test_empty_domain_raises_value_error(self) -> None:
        """An empty organism_domain should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid organism_domain"):
            resolve_organism_domain("Homo_sapiens", "")

    def test_numeric_domain_raises_value_error(self) -> None:
        """A numeric organism_domain string should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid organism_domain"):
            resolve_organism_domain("Homo_sapiens", "42")


# ─── ProteinInput model validation ─────────────────────────────────


class TestProteinInputOrganismDomain:
    """Tests for the organism_domain field on ProteinInput."""

    def test_default_organism_domain_is_auto(self) -> None:
        """Default organism_domain should be 'auto'."""
        inp = ProteinInput(protein="MVLSPADKTN", organism="Homo_sapiens")
        assert inp.organism_domain == "auto"

    def test_valid_eukaryote_domain(self) -> None:
        """Setting organism_domain to 'eukaryote' should be accepted."""
        inp = ProteinInput(
            protein="MVLSPADKTN",
            organism="Escherichia_coli",
            organism_domain="eukaryote",
        )
        assert inp.organism_domain == "eukaryote"

    def test_valid_prokaryote_domain(self) -> None:
        """Setting organism_domain to 'prokaryote' should be accepted."""
        inp = ProteinInput(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            organism_domain="prokaryote",
        )
        assert inp.organism_domain == "prokaryote"

    def test_valid_auto_domain(self) -> None:
        """Setting organism_domain to 'auto' should be accepted."""
        inp = ProteinInput(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            organism_domain="auto",
        )
        assert inp.organism_domain == "auto"

    def test_invalid_domain_raises_validation_error(self) -> None:
        """An invalid organism_domain should cause validation error."""
        from pydantic import ValidationError
        with pytest.raises(ValidationError, match="Invalid organism_domain"):
            ProteinInput(
                protein="MVLSPADKTN",
                organism="Homo_sapiens",
                organism_domain="archaea",
            )

    def test_domain_is_case_insensitive(self) -> None:
        """organism_domain should be case-insensitive."""
        inp = ProteinInput(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            organism_domain="EUKARYOTE",
        )
        assert inp.organism_domain == "eukaryote"

    def test_domain_prokaryote_case_insensitive(self) -> None:
        """organism_domain='PROKARYOTE' should normalize to 'prokaryote'."""
        inp = ProteinInput(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            organism_domain="PROKARYOTE",
        )
        assert inp.organism_domain == "prokaryote"


# ─── BatchOptimizeItem model validation ─────────────────────────────


class TestBatchOptimizeItemOrganismDomain:
    """Tests for the organism_domain field on BatchOptimizeItem."""

    def test_default_organism_domain_is_auto(self) -> None:
        """Default organism_domain on BatchOptimizeItem should be 'auto'."""
        item = BatchOptimizeItem(protein="MVLSPADKTN", organism="Homo_sapiens")
        assert item.organism_domain == "auto"

    def test_explicit_domain(self) -> None:
        """Explicit organism_domain should be stored."""
        item = BatchOptimizeItem(
            protein="MVLSPADKTN",
            organism="Escherichia_coli",
            organism_domain="prokaryote",
        )
        assert item.organism_domain == "prokaryote"

    def test_invalid_domain_raises_validation_error(self) -> None:
        """Invalid organism_domain on BatchOptimizeItem should raise validation error."""
        from pydantic import ValidationError
        with pytest.raises(ValidationError, match="Invalid organism_domain"):
            BatchOptimizeItem(
                protein="MVLSPADKTN",
                organism="Homo_sapiens",
                organism_domain="bacteria",
            )


# ─── OptimizeResponse includes domain ──────────────────────────────


class TestOptimizeResponseOrganismDomain:
    """Tests for the organism_domain field on OptimizeResponse."""

    def test_response_includes_organism_domain(self) -> None:
        """OptimizeResponse should include organism_domain field."""
        resp = OptimizeResponse(
            sequence="ATGGCT",
            protein="MA",
            cai=0.8,
            gc_content=0.5,
            satisfied_predicates=["GCInRange"],
            failed_predicates=[],
            fallback_used=False,
            organism_domain="prokaryote",
        )
        assert resp.organism_domain == "prokaryote"

    def test_response_default_domain_is_eukaryote(self) -> None:
        """OptimizeResponse default for organism_domain should be 'eukaryote'."""
        resp = OptimizeResponse(
            sequence="ATGGCT",
            protein="MA",
            cai=0.8,
            gc_content=0.5,
            satisfied_predicates=["GCInRange"],
            failed_predicates=[],
            fallback_used=False,
        )
        assert resp.organism_domain == "eukaryote"


# ─── CLI option parsing ────────────────────────────────────────────


class TestCLIOrganismDomainOption:
    """Tests for --organism-domain CLI option parsing."""

    def test_default_is_auto(self) -> None:
        """Default --organism-domain should be 'auto'."""
        parser = build_parser()
        args = parser.parse_args(["optimize", "gene.fasta"])
        assert args.organism_domain == "auto"

    def test_explicit_auto(self) -> None:
        """--organism-domain auto should parse correctly."""
        parser = build_parser()
        args = parser.parse_args([
            "optimize", "gene.fasta", "--organism-domain", "auto"
        ])
        assert args.organism_domain == "auto"

    def test_explicit_eukaryote(self) -> None:
        """--organism-domain eukaryote should parse correctly."""
        parser = build_parser()
        args = parser.parse_args([
            "optimize", "gene.fasta", "--organism-domain", "eukaryote"
        ])
        assert args.organism_domain == "eukaryote"

    def test_explicit_prokaryote(self) -> None:
        """--organism-domain prokaryote should parse correctly."""
        parser = build_parser()
        args = parser.parse_args([
            "optimize", "gene.fasta", "--organism-domain", "prokaryote"
        ])
        assert args.organism_domain == "prokaryote"

    def test_invalid_choice_rejected(self) -> None:
        """--organism-domain with invalid value should cause parser error."""
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([
                "optimize", "gene.fasta", "--organism-domain", "archaea"
            ])

    def test_organism_domain_with_species_ecoli(self) -> None:
        """--organism-domain prokaryote with --species ecoli should parse."""
        parser = build_parser()
        args = parser.parse_args([
            "optimize", "gene.fasta",
            "--species", "ecoli",
            "--organism-domain", "prokaryote",
        ])
        assert args.species == "ecoli"
        assert args.organism_domain == "prokaryote"

    def test_organism_domain_eukaryote_overrides_ecoli(self) -> None:
        """--organism-domain eukaryote with --species ecoli should parse (override)."""
        parser = build_parser()
        args = parser.parse_args([
            "optimize", "gene.fasta",
            "--species", "ecoli",
            "--organism-domain", "eukaryote",
        ])
        assert args.species == "ecoli"
        assert args.organism_domain == "eukaryote"


# ─── End-to-end domain resolution from CLI args ────────────────────


class TestDomainResolutionFromCLIArgs:
    """Test that the full resolution pipeline works from CLI args to domain."""

    def test_ecoli_auto_resolves_prokaryote(self) -> None:
        """CLI: --species ecoli with default --organism-domain resolves to prokaryote."""
        parser = build_parser()
        args = parser.parse_args(["optimize", "gene.fasta", "--species", "ecoli"])
        resolved = resolve_organism_domain(args.species, args.organism_domain)
        assert resolved == "prokaryote"

    def test_human_auto_resolves_eukaryote(self) -> None:
        """CLI: --species human with default --organism-domain resolves to eukaryote."""
        parser = build_parser()
        args = parser.parse_args(["optimize", "gene.fasta", "--species", "human"])
        resolved = resolve_organism_domain(args.species, args.organism_domain)
        assert resolved == "eukaryote"

    def test_ecoli_with_eukaryote_override(self) -> None:
        """CLI: --species ecoli with --organism-domain eukaryote resolves to eukaryote."""
        parser = build_parser()
        args = parser.parse_args([
            "optimize", "gene.fasta",
            "--species", "ecoli",
            "--organism-domain", "eukaryote",
        ])
        resolved = resolve_organism_domain(args.species, args.organism_domain)
        assert resolved == "eukaryote"

    def test_human_with_prokaryote_override(self) -> None:
        """CLI: --species human with --organism-domain prokaryote resolves to prokaryote."""
        parser = build_parser()
        args = parser.parse_args([
            "optimize", "gene.fasta",
            "--species", "human",
            "--organism-domain", "prokaryote",
        ])
        resolved = resolve_organism_domain(args.species, args.organism_domain)
        assert resolved == "prokaryote"
