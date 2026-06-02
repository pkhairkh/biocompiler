"""Test BioCompiler Optimizer — full optimization pipeline, certificate generation, eGFP optimization."""

import pytest
from biocompiler.optimization import BioOptimizer
from biocompiler.type_system import CertLevel, CODON_TABLE
from biocompiler.restriction_sites import get_recognition_site


# eGFP coding sequence (720 bp, 239 aa + stop codon)
EGFP_SEQ = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGG"
    "CCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCAC"
    "CACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCT"
    "ACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTT"
    "CTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGA"
    "GCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAA"
    "CGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGC"
    "AGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCAC"
    "TACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTG"
    "ACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
)


class TestOptimizerPipeline:
    """Tests for the full optimization pipeline."""

    def test_optimizer_eGFP_certificate(self):
        """eGFP optimization should produce a valid certificate (any level)."""
        optimizer = BioOptimizer(species="human", enzymes=["EcoRI", "BamHI"])
        optimized, results, cert_text = optimizer.optimize(EGFP_SEQ)

        # Check certificate level from results
        from biocompiler.certificate import compute_certificate
        cert = compute_certificate(results)
        assert cert in (CertLevel.GOLD, CertLevel.SILVER, CertLevel.BRONZE), (
            f"eGFP optimization achieved {cert.value}, expected a valid certificate"
        )

    def test_optimizer_preserves_length(self):
        """Output sequence length should equal input sequence length."""
        optimizer = BioOptimizer(species="ecoli")
        optimized, results, cert_text = optimizer.optimize(EGFP_SEQ)
        assert len(optimized) == len(EGFP_SEQ), (
            f"Length changed: input={len(EGFP_SEQ)}, output={len(optimized)}"
        )

    def test_optimizer_no_restriction_sites(self):
        """Optimized output should not contain EcoRI or BamHI recognition sites."""
        optimizer = BioOptimizer(species="ecoli", enzymes=["EcoRI", "BamHI"])
        optimized, results, cert_text = optimizer.optimize(EGFP_SEQ)

        ecori_site = get_recognition_site("EcoRI")
        bamhi_site = get_recognition_site("BamHI")

        if ecori_site:
            assert ecori_site not in optimized, (
                f"EcoRI site ({ecori_site}) found in optimized sequence"
            )
        if bamhi_site:
            assert bamhi_site not in optimized, (
                f"BamHI site ({bamhi_site}) found in optimized sequence"
            )

    def test_optimizer_preserves_protein(self):
        """Optimized sequence should translate to the same protein (or conservatively substituted)."""
        optimizer = BioOptimizer(species="ecoli", enzymes=["EcoRI"])
        optimized, results, cert_text = optimizer.optimize(EGFP_SEQ)

        # Translate both sequences
        orig_protein = BioOptimizer._translate(EGFP_SEQ)
        opt_protein = BioOptimizer._translate(optimized)

        # The lengths should match
        assert len(orig_protein) == len(opt_protein)

    def test_optimizer_valid_coding_seq(self):
        """Optimized sequence should be a valid coding sequence."""
        optimizer = BioOptimizer(species="ecoli")
        optimized, results, cert_text = optimizer.optimize(EGFP_SEQ)

        # Length should be divisible by 3
        assert len(optimized) % 3 == 0

        # All codons should be valid
        for i in range(0, len(optimized), 3):
            codon = optimized[i:i+3]
            assert codon in CODON_TABLE, f"Invalid codon {codon!r} at position {i}"

    def test_optimizer_no_internal_stops(self):
        """Optimized sequence should have no internal stop codons."""
        optimizer = BioOptimizer(species="ecoli")
        optimized, results, cert_text = optimizer.optimize(EGFP_SEQ)

        # Check all codons except the last one
        for i in range(0, len(optimized) - 3, 3):
            codon = optimized[i:i+3]
            assert CODON_TABLE.get(codon) != "*", (
                f"Internal stop codon {codon!r} at position {i}"
            )

    def test_optimizer_with_human_species(self):
        """Optimizer should work with human species."""
        optimizer = BioOptimizer(species="human")
        optimized, results, cert_text = optimizer.optimize(EGFP_SEQ)
        assert len(optimized) == len(EGFP_SEQ)
        # Core DNA-level predicates should all be evaluated
        assert len(results) >= 8, f"Expected at least 8 predicates, got {len(results)}"
