"""Test BioCompiler Certificate Generation — GOLD/SILVER/BRONZE computation."""

import pytest
from biocompiler.certificate import compute_certificate, format_certificate
from biocompiler.type_system import CertLevel, PredicateResult, SpliceVerdict


class TestComputeCertificate:
    """Tests for certificate level computation."""

    def test_gold_certificate(self):
        """All predicates satisfied with no mutagenesis → GOLD."""
        results = [
            PredicateResult("NoStopCodons", True, details="No internal stop codons"),
            PredicateResult("NoCrypticSplice", True, verdict=SpliceVerdict.PASS, details="No GT dinucleotides found"),
            PredicateResult("NoCpGIsland", True, details="Worst CpG Obs/Exp ratio 0.300 <= 0.6"),
            PredicateResult("NoRestrictionSite", True, details="No restriction sites found"),
            PredicateResult("NoGTDinucleotide", True, details="No GT dinucleotides found"),
            PredicateResult("ValidCodingSeq", True, details="All codons valid"),
            PredicateResult("ConservationScore", True, details="All AA conservation scores >= -1"),
            PredicateResult("CodonOptimality", True, details="Worst CAI: GCT=0.7244, min=0.0"),
        ]
        cert = compute_certificate(results)
        assert cert == CertLevel.GOLD

    def test_silver_certificate_unavoidable(self):
        """All predicates passed but some have unavoidable constraints → SILVER."""
        results = [
            PredicateResult("NoStopCodons", True, details="No internal stop codons"),
            PredicateResult("NoCrypticSplice", True, verdict=SpliceVerdict.PASS, details="No GT dinucleotides found"),
            PredicateResult("NoCpGIsland", True, details="Worst CpG Obs/Exp ratio 0.300 <= 0.6"),
            PredicateResult("NoRestrictionSite", True, details="No restriction sites found"),
            PredicateResult("NoGTDinucleotide", True, details="All 2 GT dinucleotides are unavoidable"),
            PredicateResult("ValidCodingSeq", True, details="All codons valid"),
            PredicateResult("ConservationScore", True, details="All AA conservation scores >= -1"),
            PredicateResult("CodonOptimality", True, details="Worst CAI: GCT=0.7244, min=0.0"),
        ]
        cert = compute_certificate(results)
        assert cert == CertLevel.SILVER

    def test_silver_certificate_mutagenesis(self):
        """All predicates passed but some involved mutagenesis → SILVER."""
        results = [
            PredicateResult("NoStopCodons", True, details="No internal stop codons"),
            PredicateResult("NoCrypticSplice", True, verdict=SpliceVerdict.PASS, details="No GT dinucleotides found"),
            PredicateResult("NoCpGIsland", True, details="Worst CpG Obs/Exp ratio 0.300 <= 0.6"),
            PredicateResult("NoRestrictionSite", True, details="No restriction sites found"),
            PredicateResult("NoGTDinucleotide", True, details="No GT dinucleotides found mutagenesis applied: pos 3:V→I"),
            PredicateResult("ValidCodingSeq", True, details="All codons valid"),
            PredicateResult("ConservationScore", True, details="All AA conservation scores >= -1"),
            PredicateResult("CodonOptimality", True, details="Worst CAI: GCT=0.7244, min=0.0"),
        ]
        cert = compute_certificate(results)
        assert cert == CertLevel.SILVER

    def test_bronze_certificate(self):
        """Unsatisfied predicates → BRONZE."""
        results = [
            PredicateResult("NoStopCodons", True, details="No internal stop codons"),
            PredicateResult("NoCrypticSplice", True, verdict=SpliceVerdict.PASS, details="No GT dinucleotides found"),
            PredicateResult("NoCpGIsland", True, details="Worst CpG Obs/Exp ratio 0.300 <= 0.6"),
            PredicateResult("NoRestrictionSite", True, details="No restriction sites found"),
            PredicateResult("NoGTDinucleotide", False, details="Avoidable GT dinucleotides at [3]"),
            PredicateResult("ValidCodingSeq", True, details="All codons valid"),
            PredicateResult("ConservationScore", True, details="All AA conservation scores >= -1"),
            PredicateResult("CodonOptimality", True, details="Worst CAI: GCT=0.7244, min=0.0"),
        ]
        cert = compute_certificate(results)
        assert cert == CertLevel.BRONZE


class TestFormatCertificate:
    """Tests for certificate formatting."""

    def test_format_certificate_contains_level(self):
        """Certificate text should contain the level name."""
        results = [
            PredicateResult("NoStopCodons", True, details="No internal stop codons"),
        ]
        text = format_certificate(results, "ATGGCTTAA", "ecoli")
        assert "GOLD" in text

    def test_format_certificate_contains_predicate(self):
        """Certificate text should list predicate results."""
        results = [
            PredicateResult("NoStopCodons", True, details="No internal stop codons"),
            PredicateResult("ValidCodingSeq", False, details="Sequence length 8 not divisible by 3"),
        ]
        text = format_certificate(results, "ATGGCTGC", "ecoli")
        assert "NoStopCodons" in text
        assert "ValidCodingSeq" in text
        assert "PASS" in text
        assert "FAIL" in text
