"""
End-to-end integration tests for BioCompiler.

These tests exercise the full pipeline from protein input through optimization,
constraint checking, certificate generation, and sequence export.  They verify
that the modules work together correctly as a cohesive system.

Test coverage:
  1. Optimize a short protein → get DNA sequence → verify translation
  2. Optimize with constraints → verify constraints satisfied
  3. Export to FASTA → verify format
  4. Export to GenBank → verify format
  5. Certificate generation → verify fields
  6. Full pipeline: protein input → optimization → constraint checking → export
"""

from __future__ import annotations

import hashlib
import re

import pytest

from biocompiler.optimizer import optimize_sequence, OptimizationResult
from biocompiler.export.core import export_fasta, export_genbank, export_genbank_with_certificate
from biocompiler.provenance.certificate import generate_certificate, verify_certificate
from biocompiler.expression.translation import translate, compute_cai
from biocompiler.sequence.scanner import gc_content
from biocompiler.shared.types import Certificate, TypeCheckResult, Verdict, SLOTMode
from biocompiler.type_system import (
    evaluate_all_predicates,
    check_no_restriction_site,
    check_no_stop_codons,
    check_valid_coding_seq,
    check_no_cpg_island,
)


# ─── Shared constants ───────────────────────────────────────────────────────

SHORT_PROTEIN = "MVSKGE"  # 6 AA — first 6 residues of eGFP
MEDIUM_PROTEIN = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"  # 60 AA
DEFAULT_ORGANISM = "Homo_sapiens"
DEFAULT_ENZYMES = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]


def _find_verdict(type_results, canonical_name):
    """Look up a verdict by canonical predicate name, handling
    parameterized names like 'GCInRange(0.3, 0.7)'."""
    for r in type_results:
        if r.predicate == canonical_name or r.predicate.startswith(canonical_name + "("):
            return r.verdict
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Optimize a short protein → get DNA sequence → verify translation
# ═══════════════════════════════════════════════════════════════════════════════

class TestOptimizeAndTranslate:
    """End-to-end: protein → optimization → DNA → re-translate → match."""

    def test_optimize_short_protein_translates_back(self):
        """The DNA produced by optimize_sequence should translate to the
        original protein when re-translated."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert result.sequence, "Optimized sequence must not be empty"

        # Re-translate and verify round-trip
        translated = translate(result.sequence)
        assert translated == SHORT_PROTEIN, (
            f"Translation mismatch: expected '{SHORT_PROTEIN}', got '{translated}'"
        )

    def test_optimized_sequence_length(self):
        """Optimized DNA length must equal 3 × protein length."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        assert len(result.sequence) == len(SHORT_PROTEIN) * 3, (
            f"Expected {len(SHORT_PROTEIN) * 3} bp, got {len(result.sequence)}"
        )

    def test_optimized_sequence_valid_dna(self):
        """All bases in the optimized sequence must be A/C/G/T."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        assert set(result.sequence) <= {"A", "C", "G", "T"}, (
            f"Invalid bases found: {set(result.sequence) - {'A', 'C', 'G', 'T'}}"
        )

    def test_optimize_medium_protein_translates_back(self):
        """Round-trip translation test with a longer protein (60 AA)."""
        result = optimize_sequence(MEDIUM_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        translated = translate(result.sequence)
        assert translated == MEDIUM_PROTEIN

    def test_optimize_ecoli_translates_back(self):
        """Optimization for E. coli should also preserve the protein."""
        result = optimize_sequence(SHORT_PROTEIN, organism="Escherichia_coli", strict_mode=False)
        translated = translate(result.sequence)
        assert translated == SHORT_PROTEIN

    def test_cai_is_nonzero(self):
        """The CAI score should be > 0 for any valid optimization."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        assert result.cai > 0.0, f"CAI should be positive, got {result.cai}"

    def test_gc_content_in_range(self):
        """GC content should be reported and within valid bounds."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        assert 0.0 <= result.gc_content <= 1.0


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Optimize with constraints → verify constraints satisfied
# ═══════════════════════════════════════════════════════════════════════════════

class TestOptimizeWithConstraints:
    """End-to-end: optimize with specific constraints and verify they hold."""

    def test_no_internal_stop_codons(self):
        """Optimized sequence must not contain internal stop codons."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        check = check_no_stop_codons(result.sequence)
        assert check.passed, f"Internal stop codons found: {check.details}"

    def test_valid_coding_sequence(self):
        """All codons must be valid and sequence length divisible by 3."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        check = check_valid_coding_seq(result.sequence)
        assert check.passed, f"Invalid coding sequence: {check.details}"

    def test_no_default_restriction_sites(self):
        """Optimized sequence should be free of the default enzyme sites."""
        result = optimize_sequence(
            SHORT_PROTEIN,
            organism=DEFAULT_ORGANISM,
            enzymes=DEFAULT_ENZYMES,
            strict_mode=False,
        )
        check = check_no_restriction_site(result.sequence, DEFAULT_ENZYMES)
        assert check.passed, f"Restriction sites found: {check.details}"

    def test_gc_within_specified_range(self):
        """GC content should fall within the specified [gc_lo, gc_hi] range."""
        gc_lo, gc_hi = 0.30, 0.70
        result = optimize_sequence(
            MEDIUM_PROTEIN,
            organism=DEFAULT_ORGANISM,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            strict_mode=False,
        )
        actual_gc = gc_content(result.sequence)
        assert gc_lo <= actual_gc <= gc_hi, (
            f"GC content {actual_gc:.4f} outside [{gc_lo}, {gc_hi}]"
        )

    def test_tight_gc_range_medium_protein(self):
        """Optimize with a tighter GC range and verify the optimizer attempts
        to comply.  Some proteins may not be able to reach very tight GC ranges
        due to amino acid composition constraints (e.g., amino acids with only
        low-GC codons).  We use a range that is achievable for this protein."""
        gc_lo, gc_hi = 0.30, 0.60
        result = optimize_sequence(
            MEDIUM_PROTEIN,
            organism=DEFAULT_ORGANISM,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            strict_mode=False,
        )
        actual_gc = gc_content(result.sequence)
        assert gc_lo <= actual_gc <= gc_hi, (
            f"GC content {actual_gc:.4f} outside range [{gc_lo}, {gc_hi}]"
        )

    def test_no_cpg_island(self):
        """Optimized sequence should not contain CpG islands."""
        result = optimize_sequence(
            SHORT_PROTEIN,
            organism=DEFAULT_ORGANISM,
            strict_mode=False,
        )
        check = check_no_cpg_island(result.sequence)
        assert check.passed, f"CpG island detected: {check.details}"

    def test_constraint_satisfaction_via_evaluate_all(self):
        """Run evaluate_all_predicates on the optimized sequence and verify
        that the hard-constraint predicates (InFrame, NoRestrictionSite)
        all pass."""
        result = optimize_sequence(
            MEDIUM_PROTEIN,
            organism=DEFAULT_ORGANISM,
            enzymes=DEFAULT_ENZYMES,
            gc_lo=0.30,
            gc_hi=0.70,
            strict_mode=False,
        )
        type_results = evaluate_all_predicates(
            result.sequence,
            organism=DEFAULT_ORGANISM,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
        )

        # Hard constraints must PASS (use prefix match for parameterized names)
        assert _find_verdict(type_results, "InFrame") == Verdict.PASS, (
            f"InFrame: {_find_verdict(type_results, 'InFrame')}"
        )
        assert _find_verdict(type_results, "NoRestrictionSite") == Verdict.PASS, (
            f"NoRestrictionSite: {_find_verdict(type_results, 'NoRestrictionSite')}"
        )
        assert _find_verdict(type_results, "GCInRange") == Verdict.PASS, (
            f"GCInRange: {_find_verdict(type_results, 'GCInRange')}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Export to FASTA → verify format
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportFastaE2E:
    """End-to-end: optimize → export FASTA → verify format compliance."""

    def test_fasta_header_starts_with_greater_than(self):
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        fasta = export_fasta(
            result.sequence,
            identifier="test_gene",
            organism=DEFAULT_ORGANISM,
            protein=SHORT_PROTEIN,
        )
        assert fasta.startswith(">"), "FASTA must start with '>'"

    def test_fasta_header_contains_organism(self):
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        fasta = export_fasta(
            result.sequence,
            identifier="test_gene",
            organism=DEFAULT_ORGANISM,
        )
        assert f"organism={DEFAULT_ORGANISM}" in fasta

    def test_fasta_header_contains_gc_and_length(self):
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        fasta = export_fasta(result.sequence, identifier="test_gene")
        assert "gc=" in fasta
        assert "len=" in fasta

    def test_fasta_sequence_preserves_optimized_dna(self):
        """All bases from the optimized DNA must appear in the FASTA body."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        fasta = export_fasta(result.sequence, identifier="test_gene")
        # Extract sequence lines (skip header)
        lines = fasta.strip().split("\n")
        seq_body = "".join(lines[1:]).replace(" ", "")
        assert seq_body == result.sequence, "FASTA sequence must match optimized DNA"

    def test_fasta_wrapping_at_60_chars(self):
        """Sequence lines should be at most 60 characters each."""
        result = optimize_sequence(MEDIUM_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        fasta = export_fasta(result.sequence, identifier="test_gene")
        lines = fasta.strip().split("\n")
        seq_lines = lines[1:]  # skip header
        for line in seq_lines:
            assert len(line) <= 60, f"FASTA line exceeds 60 chars: '{line}'"

    def test_fasta_translation_matches_protein(self):
        """Auto-computed protein translation in header should match."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        fasta = export_fasta(result.sequence, identifier="test_gene")
        assert "protein_len=" in fasta
        assert f"{len(SHORT_PROTEIN)}aa" in fasta

    def test_fasta_ends_with_newline(self):
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        fasta = export_fasta(result.sequence, identifier="test_gene")
        assert fasta.endswith("\n"), "FASTA output should end with a newline"


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Export to GenBank → verify format
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportGenbankE2E:
    """End-to-end: optimize → export GenBank → verify format compliance."""

    def test_genbank_has_required_sections(self):
        """GenBank output must contain LOCUS, FEATURES, ORIGIN, and //."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        gb = export_genbank(
            result.sequence,
            locus_name="TESTGENE",
            organism=DEFAULT_ORGANISM,
            gene_name="test_gene",
            protein=SHORT_PROTEIN,
        )
        assert "LOCUS" in gb, "Missing LOCUS section"
        assert "FEATURES" in gb, "Missing FEATURES section"
        assert "ORIGIN" in gb, "Missing ORIGIN section"
        assert gb.rstrip().endswith("//"), "Missing GenBank terminator //"

    def test_genbank_section_ordering(self):
        """Sections must appear in correct order: LOCUS < FEATURES < ORIGIN < //."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        gb = export_genbank(result.sequence, locus_name="TESTGENE")
        assert gb.index("LOCUS") < gb.index("FEATURES") < gb.index("ORIGIN")
        assert gb.index("ORIGIN") < gb.index("//")

    def test_genbank_locus_contains_bp(self):
        """LOCUS line must report the correct sequence length in bp."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        gb = export_genbank(result.sequence, locus_name="TESTGENE")
        expected_bp = f"{len(result.sequence)} bp"
        assert expected_bp in gb, f"Expected '{expected_bp}' in GenBank output"

    def test_genbank_cds_translation(self):
        """CDS feature must contain a /translation qualifier."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        gb = export_genbank(
            result.sequence,
            gene_name="test_gene",
            protein=SHORT_PROTEIN,
        )
        assert "/translation=" in gb, "Missing /translation qualifier in CDS"
        assert SHORT_PROTEIN in gb, "Protein sequence not found in GenBank output"

    def test_genbank_preserves_sequence_content(self):
        """All bases from the optimized DNA must appear in the ORIGIN section."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        gb = export_genbank(result.sequence, locus_name="TESTGENE")
        # Extract ORIGIN section, strip numbering and spaces
        origin_start = gb.index("ORIGIN")
        terminator = gb.index("//")
        origin_section = gb[origin_start:terminator]
        bases = re.sub(r'\d+', '', origin_section.replace("ORIGIN", ""))
        bases = bases.replace(" ", "").replace("\n", "")
        assert bases == result.sequence, "GenBank ORIGIN bases must match optimized DNA"

    def test_genbank_with_gene_name(self):
        """Gene name should appear in both gene and CDS features."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        gb = export_genbank(
            result.sequence,
            gene_name="eGFP",
            protein=SHORT_PROTEIN,
        )
        assert '/gene="eGFP"' in gb

    def test_genbank_accession_present(self):
        """ACCESSION line must be present."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        gb = export_genbank(result.sequence, locus_name="TESTGENE")
        assert "ACCESSION" in gb

    def test_genbank_source_organism(self):
        """SOURCE line must contain the organism name."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        gb = export_genbank(result.sequence, organism=DEFAULT_ORGANISM)
        assert "SOURCE" in gb
        assert DEFAULT_ORGANISM in gb

    def test_genbank_comment_includes_gc(self):
        """COMMENT section should include GC content."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        gb = export_genbank(result.sequence)
        assert "GC content:" in gb


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Certificate generation → verify fields
# ═══════════════════════════════════════════════════════════════════════════════

class TestCertificateGenerationE2E:
    """End-to-end: optimize → evaluate predicates → generate certificate → verify."""

    def test_generate_certificate_from_optimized_sequence(self):
        """Certificate should contain all required fields."""
        result = optimize_sequence(
            MEDIUM_PROTEIN,
            organism=DEFAULT_ORGANISM,
            enzymes=DEFAULT_ENZYMES,
            strict_mode=False,
        )
        type_results = evaluate_all_predicates(
            result.sequence,
            organism=DEFAULT_ORGANISM,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
        )
        input_params = {
            "gene": "test",
            "organism": DEFAULT_ORGANISM,
            "gc_lo": 0.30,
            "gc_hi": 0.70,
            "enzymes": DEFAULT_ENZYMES,
            "exon_boundaries": [(0, len(result.sequence))],
        }
        cert = generate_certificate(
            sequence=result.sequence,
            type_results=type_results,
            input_params=input_params,
        )
        assert isinstance(cert, Certificate)
        assert cert.version, "Certificate version must not be empty"
        assert cert.sequence == result.sequence
        assert cert.design_id, "Certificate design_id must not be empty"
        assert len(cert.types) > 0, "Certificate must contain at least one type result"
        assert "tool" in cert.provenance
        assert "version" in cert.provenance
        assert "timestamp" in cert.provenance
        assert "input_hash" in cert.provenance

    def test_certificate_design_id_matches_sequence_hash(self):
        """design_id should be the SHA-256 hash of the sequence."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        type_results = evaluate_all_predicates(result.sequence, organism=DEFAULT_ORGANISM)
        input_params = {"organism": DEFAULT_ORGANISM}
        cert = generate_certificate(
            sequence=result.sequence,
            type_results=type_results,
            input_params=input_params,
        )
        expected_hash = hashlib.sha256(result.sequence.encode()).hexdigest()
        assert cert.design_id == expected_hash, (
            f"design_id mismatch: cert={cert.design_id[:16]}... "
            f"expected={expected_hash[:16]}..."
        )

    def test_certificate_types_contain_verdicts(self):
        """Each type entry must have 'predicate' and 'verdict' keys."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        type_results = evaluate_all_predicates(result.sequence, organism=DEFAULT_ORGANISM)
        input_params = {"organism": DEFAULT_ORGANISM}
        cert = generate_certificate(
            sequence=result.sequence,
            type_results=type_results,
            input_params=input_params,
        )
        for t in cert.types:
            assert "predicate" in t, f"Missing 'predicate' in type entry: {t}"
            assert "verdict" in t, f"Missing 'verdict' in type entry: {t}"
            assert t["verdict"] in ("PASS", "LIKELY_PASS", "UNCERTAIN", "LIKELY_FAIL", "FAIL"), (
                f"Invalid verdict: {t['verdict']}"
            )

    def test_certificate_provenance_has_parameters(self):
        """Provenance should contain the full input parameters."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        type_results = evaluate_all_predicates(result.sequence, organism=DEFAULT_ORGANISM)
        input_params = {
            "organism": DEFAULT_ORGANISM,
            "gc_lo": 0.30,
            "gc_hi": 0.70,
        }
        cert = generate_certificate(
            sequence=result.sequence,
            type_results=type_results,
            input_params=input_params,
        )
        params = cert.provenance.get("parameters", {})
        assert params.get("organism") == DEFAULT_ORGANISM
        assert params.get("gc_lo") == 0.30
        assert params.get("gc_hi") == 0.70

    def test_certificate_serialization_round_trip(self):
        """Certificate should survive to_dict → from_dict round-trip."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        type_results = evaluate_all_predicates(result.sequence, organism=DEFAULT_ORGANISM)
        cert = generate_certificate(
            sequence=result.sequence,
            type_results=type_results,
            input_params={"organism": DEFAULT_ORGANISM},
        )
        cert_dict = cert.to_dict()
        restored = Certificate.from_dict(cert_dict)
        assert restored.version == cert.version
        assert restored.design_id == cert.design_id
        assert restored.sequence == cert.sequence
        assert len(restored.types) == len(cert.types)

    def test_certificate_verify_matches_design_id(self):
        """Independent verification should confirm design_id matches sequence."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        type_results = evaluate_all_predicates(result.sequence, organism=DEFAULT_ORGANISM)
        cert = generate_certificate(
            sequence=result.sequence,
            type_results=type_results,
            input_params={"organism": DEFAULT_ORGANISM},
        )
        cert_dict = cert.to_dict()
        # Verify certificate — at minimum, design_id must match
        computed_hash = hashlib.sha256(result.sequence.encode()).hexdigest()
        assert cert_dict["design_id"] == computed_hash, (
            "Certificate design_id does not match SHA-256 of sequence"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Full pipeline: protein input → optimization → constraint checking → export
# ═══════════════════════════════════════════════════════════════════════════════

class TestFullPipeline:
    """End-to-end: the complete pipeline from protein to export with certificate."""

    def test_full_pipeline_short_protein(self):
        """Complete pipeline: optimize → evaluate → certificate → export FASTA + GenBank."""
        # Step 1: Optimize
        result = optimize_sequence(
            SHORT_PROTEIN,
            organism=DEFAULT_ORGANISM,
            enzymes=DEFAULT_ENZYMES,
            gc_lo=0.30,
            gc_hi=0.70,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert result.sequence

        # Step 2: Verify translation round-trip
        translated = translate(result.sequence)
        assert translated == SHORT_PROTEIN

        # Step 3: Evaluate all predicates
        type_results = evaluate_all_predicates(
            result.sequence,
            organism=DEFAULT_ORGANISM,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
        )
        assert len(type_results) > 0, "evaluate_all_predicates returned no results"

        # Verify hard constraints pass
        assert _find_verdict(type_results, "InFrame") == Verdict.PASS
        assert _find_verdict(type_results, "NoRestrictionSite") == Verdict.PASS
        assert _find_verdict(type_results, "GCInRange") == Verdict.PASS

        # Step 4: Generate certificate
        input_params = {
            "gene": "test_short",
            "organism": DEFAULT_ORGANISM,
            "gc_lo": 0.30,
            "gc_hi": 0.70,
            "enzymes": DEFAULT_ENZYMES,
            "exon_boundaries": [(0, len(result.sequence))],
        }
        cert = generate_certificate(
            sequence=result.sequence,
            type_results=type_results,
            input_params=input_params,
        )
        assert isinstance(cert, Certificate)
        assert cert.design_id

        # Step 5: Export to FASTA
        fasta = export_fasta(
            result.sequence,
            identifier="pipeline_test",
            description="Full pipeline test",
            organism=DEFAULT_ORGANISM,
            protein=SHORT_PROTEIN,
        )
        assert fasta.startswith(">")
        assert "organism=Homo_sapiens" in fasta
        assert "pipeline_test" in fasta
        # Verify FASTA sequence body matches optimized DNA
        fasta_lines = fasta.strip().split("\n")
        fasta_seq = "".join(fasta_lines[1:])
        assert fasta_seq == result.sequence

        # Step 6: Export to GenBank (plain)
        gb = export_genbank(
            result.sequence,
            locus_name="PIPELINETEST",
            organism=DEFAULT_ORGANISM,
            gene_name="test_gene",
            protein=SHORT_PROTEIN,
            certificate=cert,
            type_results=type_results,
        )
        assert "LOCUS" in gb
        assert "FEATURES" in gb
        assert "ORIGIN" in gb
        assert gb.rstrip().endswith("//")
        assert "Certificate ID:" in gb

        # Step 7: Export to GenBank with certificate
        gb_cert = export_genbank_with_certificate(
            result.sequence,
            certificate=cert,
            organism=DEFAULT_ORGANISM,
            gene_name="test_gene",
        )
        assert "LOCUS" in gb_cert
        assert "Certificate ID:" in gb_cert

    def test_full_pipeline_medium_protein(self):
        """Full pipeline with a 60-AA protein to catch scaling issues."""
        # Optimize
        result = optimize_sequence(
            MEDIUM_PROTEIN,
            organism=DEFAULT_ORGANISM,
            enzymes=DEFAULT_ENZYMES,
            gc_lo=0.30,
            gc_hi=0.70,
            strict_mode=False,
        )
        # Translate
        assert translate(result.sequence) == MEDIUM_PROTEIN
        # Evaluate
        type_results = evaluate_all_predicates(
            result.sequence,
            organism=DEFAULT_ORGANISM,
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
        )
        assert _find_verdict(type_results, "InFrame") == Verdict.PASS
        assert _find_verdict(type_results, "NoRestrictionSite") == Verdict.PASS
        # Certificate
        cert = generate_certificate(
            sequence=result.sequence,
            type_results=type_results,
            input_params={
                "organism": DEFAULT_ORGANISM,
                "gc_lo": 0.30,
                "gc_hi": 0.70,
                "enzymes": DEFAULT_ENZYMES,
            },
        )
        assert isinstance(cert, Certificate)
        # FASTA export
        fasta = export_fasta(result.sequence, identifier="med_test")
        assert fasta.startswith(">")
        # GenBank export
        gb = export_genbank(result.sequence, locus_name="MEDTEST", protein=MEDIUM_PROTEIN)
        assert "LOCUS" in gb and "ORIGIN" in gb and gb.rstrip().endswith("//")

    def test_full_pipeline_ecoli(self):
        """Full pipeline targeting E. coli instead of human."""
        result = optimize_sequence(
            SHORT_PROTEIN,
            organism="Escherichia_coli",
            enzymes=DEFAULT_ENZYMES,
            strict_mode=False,
        )
        assert translate(result.sequence) == SHORT_PROTEIN
        type_results = evaluate_all_predicates(
            result.sequence,
            organism="Escherichia_coli",
            enzymes=DEFAULT_ENZYMES,
        )
        assert _find_verdict(type_results, "InFrame") == Verdict.PASS
        assert _find_verdict(type_results, "NoRestrictionSite") == Verdict.PASS

        cert = generate_certificate(
            sequence=result.sequence,
            type_results=type_results,
            input_params={"organism": "Escherichia_coli", "enzymes": DEFAULT_ENZYMES},
        )
        assert isinstance(cert, Certificate)

        # FASTA should reference E. coli
        fasta = export_fasta(
            result.sequence,
            identifier="ecoli_test",
            organism="Escherichia_coli",
        )
        assert "organism=Escherichia_coli" in fasta

        # GenBank should reference E. coli taxonomy
        gb = export_genbank(
            result.sequence,
            organism="Escherichia_coli",
            protein=SHORT_PROTEIN,
        )
        assert "Bacteria" in gb

    def test_pipeline_gc_constraint_enforced_in_output(self):
        """Verify that the GC constraint is actually enforced in the exported
        sequences, not just reported as passing."""
        gc_lo, gc_hi = 0.30, 0.70
        result = optimize_sequence(
            MEDIUM_PROTEIN,
            organism=DEFAULT_ORGANISM,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            strict_mode=False,
        )
        actual_gc = gc_content(result.sequence)
        assert gc_lo <= actual_gc <= gc_hi, (
            f"GC {actual_gc:.4f} not in [{gc_lo}, {gc_hi}]"
        )
        # Also verify in FASTA header
        fasta = export_fasta(result.sequence, identifier="gc_test")
        for part in fasta.split("\n")[0].split("|"):
            if part.startswith("gc="):
                reported_gc = float(part.split("=")[1])
                assert abs(reported_gc - actual_gc) < 0.01, (
                    f"FASTA gc={reported_gc} does not match actual {actual_gc:.4f}"
                )
        # Also verify in GenBank COMMENT
        gb = export_genbank(result.sequence)
        for line in gb.split("\n"):
            if "GC content:" in line:
                # Extract the GC value from the comment
                gc_match = re.search(r"GC content:\s*([\d.]+)", line)
                assert gc_match, f"Could not parse GC from: {line}"
                gb_gc = float(gc_match.group(1))
                assert abs(gb_gc - actual_gc) < 0.01, (
                    f"GenBank GC={gb_gc} does not match actual {actual_gc:.4f}"
                )

    def test_pipeline_certificate_embedded_in_genbank(self):
        """When a certificate is provided, it should be embedded in the
        GenBank COMMENT section."""
        result = optimize_sequence(SHORT_PROTEIN, organism=DEFAULT_ORGANISM, strict_mode=False)
        type_results = evaluate_all_predicates(result.sequence, organism=DEFAULT_ORGANISM)
        cert = generate_certificate(
            sequence=result.sequence,
            type_results=type_results,
            input_params={"organism": DEFAULT_ORGANISM},
        )
        gb = export_genbank_with_certificate(
            result.sequence,
            certificate=cert,
            organism=DEFAULT_ORGANISM,
        )
        assert "Certificate ID:" in gb
        assert "Certificate timestamp:" in gb
        assert "Type-check verdict:" in gb
