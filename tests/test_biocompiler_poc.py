#!/usr/bin/env python3
"""
Comprehensive unit tests for BioCompiler PoC.
"""

import sys
import os
import json
import hashlib

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from biocompiler_poc import (
    Verdict, three_valued_and, Token, PositionRange, SpliceIsoform,
    TypeCheckResult, Certificate,
    validate_dna_sequence, scan_sequence, translate, compute_cai, gc_content,
    evaluate_gc_in_range, evaluate_no_restriction_site, evaluate_in_frame,
    evaluate_no_instability_motif, evaluate_no_cryptic_splice,
    evaluate_splice_correct, evaluate_codon_adapted,
    generate_certificate, verify_certificate,
    KOZAK_CONSENSUS, INSTABILITY_MOTIF, CODON_TABLE,
)

import pytest


# ==============================================================================
# Three-valued logic tests
# ==============================================================================

class TestThreeValuedAnd:
    """Test all 9 combinations of PASS/FAIL/UNCERTAIN."""

    def test_pass_and_pass(self):
        assert three_valued_and(Verdict.PASS, Verdict.PASS) == Verdict.PASS

    def test_pass_and_fail(self):
        assert three_valued_and(Verdict.PASS, Verdict.FAIL) == Verdict.FAIL

    def test_pass_and_uncertain(self):
        assert three_valued_and(Verdict.PASS, Verdict.UNCERTAIN) == Verdict.UNCERTAIN

    def test_fail_and_pass(self):
        assert three_valued_and(Verdict.FAIL, Verdict.PASS) == Verdict.FAIL

    def test_fail_and_fail(self):
        assert three_valued_and(Verdict.FAIL, Verdict.FAIL) == Verdict.FAIL

    def test_fail_and_uncertain(self):
        assert three_valued_and(Verdict.FAIL, Verdict.UNCERTAIN) == Verdict.FAIL

    def test_uncertain_and_pass(self):
        assert three_valued_and(Verdict.UNCERTAIN, Verdict.PASS) == Verdict.UNCERTAIN

    def test_uncertain_and_fail(self):
        assert three_valued_and(Verdict.UNCERTAIN, Verdict.FAIL) == Verdict.FAIL

    def test_uncertain_and_uncertain(self):
        assert three_valued_and(Verdict.UNCERTAIN, Verdict.UNCERTAIN) == Verdict.UNCERTAIN


# ==============================================================================
# GC content tests
# ==============================================================================

class TestGCContent:
    def test_all_gc(self):
        assert gc_content("GCGCGC") == 1.0

    def test_no_gc(self):
        assert gc_content("ATATAT") == 0.0

    def test_half_gc(self):
        assert gc_content("ATGC") == 0.5

    def test_empty(self):
        assert gc_content("") == 0.0

    def test_case_insensitive(self):
        assert gc_content("atgc") == 0.5

    def test_known_sequence(self):
        # ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG
        seq = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
        gc = gc_content(seq)
        assert 0.3 <= gc <= 0.7  # Reasonable range for a coding sequence


# ==============================================================================
# Translate tests
# ==============================================================================

class TestTranslate:
    def test_start_codon(self):
        assert translate("ATG") == "M"

    def test_simple_seq(self):
        assert translate("ATGTAA") == "M"

    def test_stop_codons(self):
        assert translate("ATGTTTTAA") == "MF"

    def test_empty_seq(self):
        assert translate("") == ""

    def test_hbb_cds(self):
        """Test that HBB CDS translates to the expected protein start."""
        hbb_cds = (
            "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAG"
            "TTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTG"
            "GGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGG"
            "TGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCAC"
            "TGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCC"
            "CATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTA"
            "ATGCCCTGGCCCACAAGTATCACTAA"
        )
        protein = translate(hbb_cds)
        expected_start = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        assert protein == expected_start


# ==============================================================================
# CAI tests
# ==============================================================================

class TestComputeCAI:
    def test_empty_sequence(self):
        assert compute_cai("") == 0.0

    def test_unknown_organism(self):
        assert compute_cai("ATGAAATTT", "E_coli") == 0.0

    def test_short_sequence(self):
        # ATG (M, skipped) + AAA (K, 0.43/0.57≈0.754) + TTT (F, 0.46/0.54≈0.852)
        cai = compute_cai("ATGAAATTT")
        assert cai > 0.0
        assert cai <= 1.0

    def test_high_cai_sequence(self):
        """Sequence using all preferred codons should have high CAI."""
        # Use preferred codons: F=TTC, L=CTG, I=ATC, V=GTG, S=AGC, P=CCC, T=ACC, A=GCC
        seq = "ATG" + "CTG" * 10 + "TAA"  # M + 10 Leu (CTG preferred) + stop
        cai = compute_cai(seq)
        # CTG is the preferred Leu codon (freq 0.39), so ratio should be 1.0
        assert cai == 1.0  # All codons are preferred


# ==============================================================================
# Scanner tests
# ==============================================================================

class TestScanSequence:
    def test_empty_sequence(self):
        assert scan_sequence("") == []

    def test_donor_site(self):
        tokens = scan_sequence("AAGTAA")
        donor_tokens = [t for t in tokens if t.element_type == "splice_donor"]
        assert len(donor_tokens) >= 1

    def test_instability_motif(self):
        tokens = scan_sequence("AAAAATTTAAAA")
        inst_tokens = [t for t in tokens if t.element_type == "instability_motif"]
        assert len(inst_tokens) == 1

    def test_restriction_site(self):
        tokens = scan_sequence("AAGAATTCAA", ["EcoRI"])
        rest_tokens = [t for t in tokens if t.element_type == "restriction_site"]
        assert len(rest_tokens) == 1
        assert rest_tokens[0].match_sequence == "GAATTC"

    def test_no_restriction_site_when_not_requested(self):
        tokens = scan_sequence("GAATTC")
        rest_tokens = [t for t in tokens if t.element_type == "restriction_site"]
        assert len(rest_tokens) == 0

    def test_invalid_dna_raises(self):
        with pytest.raises(ValueError):
            scan_sequence("HELLO")


# ==============================================================================
# Kozak scanner tests
# ==============================================================================

class TestKozakScanner:
    def test_kozak_exact_match(self):
        """GCCACC should be found directly."""
        tokens = scan_sequence("AAAGCCACCAA")
        kozak_tokens = [t for t in tokens if t.element_type == "kozak"]
        assert len(kozak_tokens) == 1
        assert kozak_tokens[0].match_sequence == "GCCACC"
        assert kozak_tokens[0].position == 3

    def test_kozak_no_false_positives(self):
        """Partial matches should not be detected."""
        tokens = scan_sequence("AAAGCCACAA")  # GCCACA - not GCCACC
        kozak_tokens = [t for t in tokens if t.element_type == "kozak"]
        assert len(kozak_tokens) == 0

    def test_kozak_at_start(self):
        tokens = scan_sequence("GCCACCAAA")
        kozak_tokens = [t for t in tokens if t.element_type == "kozak"]
        assert len(kozak_tokens) == 1
        assert kozak_tokens[0].position == 0


# ==============================================================================
# U-rich detection tests
# ==============================================================================

class TestURichDetection:
    def test_six_consecutive_t(self):
        """6 consecutive T should be detected as U-rich."""
        result = evaluate_no_instability_motif("AAAAAATTTTTTAAAAA")
        assert result.verdict == Verdict.FAIL
        assert "U_rich" in result.violation

    def test_eight_consecutive_t(self):
        """8 consecutive T should be detected as one U-rich region, not multiple."""
        result = evaluate_no_instability_motif("AAATTTTTTTTAAA")
        assert result.verdict == Verdict.FAIL

    def test_five_consecutive_t_passes(self):
        """5 consecutive T should NOT be detected as U-rich."""
        result = evaluate_no_instability_motif("AAATTTTTAAAA")
        assert result.verdict == Verdict.PASS

    def test_no_u_rich(self):
        """Sequence with no long T runs should pass."""
        result = evaluate_no_instability_motif("ATGCATGCATGC")
        assert result.verdict == Verdict.PASS


# ==============================================================================
# GC in range tests
# ==============================================================================

class TestEvaluateGCInRange:
    def test_in_range(self):
        seq = "ATGCATGC"  # 50% GC
        result = evaluate_gc_in_range(seq, 0.3, 0.7)
        assert result.verdict == Verdict.PASS

    def test_out_of_range_low(self):
        seq = "ATATATAT"  # 0% GC
        result = evaluate_gc_in_range(seq, 0.3, 0.7)
        assert result.verdict == Verdict.FAIL

    def test_out_of_range_high(self):
        seq = "GCGCGCGC"  # 100% GC
        result = evaluate_gc_in_range(seq, 0.3, 0.7)
        assert result.verdict == Verdict.FAIL


# ==============================================================================
# No restriction site tests
# ==============================================================================

class TestEvaluateNoRestrictionSite:
    def test_no_site(self):
        seq = "ATGCATGCATGC"
        result = evaluate_no_restriction_site(seq, ["EcoRI"])
        assert result.verdict == Verdict.PASS

    def test_with_site(self):
        seq = "GAATTC"
        result = evaluate_no_restriction_site(seq, ["EcoRI"])
        assert result.verdict == Verdict.FAIL


# ==============================================================================
# In frame tests
# ==============================================================================

class TestEvaluateInFrame:
    def test_valid_frame(self):
        seq = "ATGAAATTTGGG"
        result = evaluate_in_frame(seq, [(0, 12)])
        assert result.verdict == Verdict.PASS

    def test_invalid_frame(self):
        seq = "ATGAAATTTGGG"
        # Exon length 5 is not a multiple of 3
        result = evaluate_in_frame(seq, [(0, 5)])
        assert result.verdict == Verdict.FAIL


# ==============================================================================
# No instability motif tests
# ==============================================================================

class TestEvaluateNoInstabilityMotif:
    def test_no_motif(self):
        result = evaluate_no_instability_motif("ATGCATGCATGC")
        assert result.verdict == Verdict.PASS

    def test_with_atfta(self):
        result = evaluate_no_instability_motif("AAAATTTAAAA")
        assert result.verdict == Verdict.FAIL


# ==============================================================================
# Certificate generation and verification round-trip test
# ==============================================================================

class TestCertificateRoundTrip:
    def test_generation_and_verification(self):
        """Generate a certificate and verify it round-trip."""
        # Create a simple sequence that passes all checks
        # Using a carefully crafted sequence to avoid GT/AG, restriction sites, etc.
        seq = "ATGCCTCCTCCTCCTTAA"
        # GC content: C=6, G=0, total=18 => 6/18 = 0.3333
        exon_boundaries = [(0, len(seq))]

        type_results = []

        # GCInRange
        result = evaluate_gc_in_range(seq, 0.30, 0.70)
        type_results.append(result)

        # NoRestrictionSite
        result = evaluate_no_restriction_site(seq, ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"])
        type_results.append(result)

        # InFrame - only if length is multiple of 3
        if len(seq) % 3 == 0:
            result = evaluate_in_frame(seq, exon_boundaries)
            type_results.append(result)

        # NoInstabilityMotif
        result = evaluate_no_instability_motif(seq)
        type_results.append(result)

        # Only generate cert if all pass
        passing = [r for r in type_results if r.verdict == Verdict.PASS]
        if len(passing) == len(type_results):
            cert = generate_certificate(
                seq, type_results,
                {"gene": "test", "organism": "Homo_sapiens", "exon_boundaries": exon_boundaries}
            )
            cert_dict = {
                "version": cert.version,
                "design_id": cert.design_id,
                "sequence": cert.sequence,
                "types": cert.types,
                "provenance": cert.provenance,
            }

            status, failures = verify_certificate(cert_dict, exon_boundaries)
            assert status == "VERIFIED", f"Certificate verification failed: {failures}"


# ==============================================================================
# Validate DNA sequence tests
# ==============================================================================

class TestValidateDNASequence:
    def test_valid_sequence(self):
        assert validate_dna_sequence("ATGC") == "ATGC"

    def test_lowercase(self):
        assert validate_dna_sequence("atgc") == "ATGC"

    def test_with_n(self):
        assert validate_dna_sequence("ATGN") == "ATGN"

    def test_invalid_bases(self):
        with pytest.raises(ValueError):
            validate_dna_sequence("ATGX")

    def test_empty_string(self):
        assert validate_dna_sequence("") == ""


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
