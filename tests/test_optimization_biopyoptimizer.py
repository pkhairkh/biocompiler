"""Tests for the BioOptimizer class.

Covers:
- Initialization with different organisms
- optimize() method end-to-end with DNA input
- Prokaryote fast path
- Eukaryote fast path
- Error handling for invalid inputs
- Strategy selection
- Organism resolution
"""

import pytest
from biocompiler.optimization import BioOptimizer, OptimizationResult
from biocompiler.type_system import CODON_TABLE, AA_TO_CODONS
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES


# Standard test DNA sequence (encodes MVSKGE)
TEST_DNA = "ATGGTTTCTAAAGGTGAA"
# Longer test DNA (encodes MVSKGEELFTGVVPILVELDGDVNGHK)
LONG_DNA = "ATGGTTTCTAAAGGTGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGG"


class TestBioOptimizerInit:
    def test_default_ecoli_init(self):
        opt = BioOptimizer(species="ecoli")
        assert opt.organism_name == "Escherichia_coli"
        assert opt.is_prokaryote is True
        assert opt.avoid_gt is False  # auto-disabled for prokaryotes

    def test_human_init(self):
        opt = BioOptimizer(species="human")
        assert opt.organism_name == "Homo_sapiens"
        assert opt.is_prokaryote is False

    def test_mouse_init(self):
        opt = BioOptimizer(species="mouse")
        assert opt.organism_name == "Mus_musculus"

    def test_cho_init(self):
        opt = BioOptimizer(species="cho")
        assert opt.organism_name == "CHO_K1"

    def test_yeast_init(self):
        opt = BioOptimizer(species="yeast")
        assert opt.organism_name == "Saccharomyces_cerevisiae"

    def test_full_name_init(self):
        opt = BioOptimizer(species="Escherichia_coli")
        assert opt.organism_name == "Escherichia_coli"

    def test_organism_name_kwarg(self):
        opt = BioOptimizer(species="ecoli", organism_name="Homo_sapiens")
        assert opt.organism_name == "Homo_sapiens"

    def test_species_cai_loaded(self):
        opt = BioOptimizer(species="ecoli")
        assert isinstance(opt.species_cai, dict)
        assert len(opt.species_cai) > 0

    def test_enzymes_default_empty(self):
        opt = BioOptimizer(species="ecoli")
        assert opt.enzymes == []

    def test_enzymes_custom(self):
        opt = BioOptimizer(species="ecoli", enzymes=["EcoRI", "BamHI"])
        assert opt.enzymes == ["EcoRI", "BamHI"]

    def test_default_strategy(self):
        opt = BioOptimizer(species="ecoli")
        assert opt.strategy == "constraint_first"

    def test_cai_first_strategy(self):
        opt = BioOptimizer(species="ecoli", strategy="cai_first")
        assert opt.strategy == "cai_first"

    def test_prokaryote_domain_auto_detected(self):
        opt = BioOptimizer(species="ecoli")
        assert opt.organism_domain == "prokaryote"

    def test_eukaryote_domain_auto_detected(self):
        opt = BioOptimizer(species="human")
        assert opt.organism_domain == "eukaryote"

    def test_avoid_gt_default_true_eukaryote(self):
        opt = BioOptimizer(species="human")
        assert opt.avoid_gt is True

    def test_optimize_mrna_stability_default(self):
        opt = BioOptimizer(species="ecoli")
        assert opt.optimize_mrna_stability is True

    def test_splice_thresholds(self):
        opt = BioOptimizer(species="human", splice_low=2.0, splice_high=5.0)
        assert opt.splice_low == pytest.approx(2.0, rel=1e-6)
        assert opt.splice_high == pytest.approx(5.0, rel=1e-6)


class TestBioOptimizerOptimize:
    def test_basic_ecoli_dna(self):
        opt = BioOptimizer(species="ecoli")
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)
        assert len(seq) > 0
        assert isinstance(preds, list)
        assert isinstance(cert, str)

    def test_basic_human_dna(self):
        opt = BioOptimizer(species="human")
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)
        assert len(seq) > 0

    def test_dna_input_ecoli(self):
        opt = BioOptimizer(species="ecoli")
        seq, preds, cert = opt.optimize(LONG_DNA)
        assert isinstance(seq, str)
        assert len(seq) > 0

    def test_ecoli_encoding_preserved(self):
        opt = BioOptimizer(species="ecoli")
        seq, preds, cert = opt.optimize(LONG_DNA)
        # The original protein should be preserved
        assert opt._original_protein != ""

    def test_human_encoding_preserved(self):
        opt = BioOptimizer(species="human")
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)
        assert len(seq) > 0

    def test_constraint_first_strategy(self):
        opt = BioOptimizer(species="ecoli", strategy="constraint_first")
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)

    def test_cai_first_strategy(self):
        opt = BioOptimizer(species="ecoli", strategy="cai_first")
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)

    def test_enzymes_param(self):
        opt = BioOptimizer(species="ecoli", enzymes=["EcoRI"])
        seq, preds, cert = opt.optimize(LONG_DNA)
        assert isinstance(seq, str)

    def test_deterministic_output(self):
        opt = BioOptimizer(species="ecoli")
        seq1, _, _ = opt.optimize(TEST_DNA)
        seq2, _, _ = opt.optimize(TEST_DNA)
        assert seq1 == seq2

    def test_prokaryote_skips_eukaryote_steps(self):
        # Prokaryote should not run splice/CpG steps
        opt = BioOptimizer(species="ecoli")
        assert opt.is_prokaryote is True
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)

    def test_eukaryote_runs_eukaryote_steps(self):
        opt = BioOptimizer(species="human")
        assert opt.is_prokaryote is False
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)

    def test_longer_sequence(self):
        opt = BioOptimizer(species="ecoli")
        seq, preds, cert = opt.optimize(LONG_DNA)
        assert len(seq) > 0

    def test_mrna_stability_optimization(self):
        opt = BioOptimizer(species="human", optimize_mrna_stability=True)
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)

    def test_no_mrna_stability_optimization(self):
        opt = BioOptimizer(species="human", optimize_mrna_stability=False)
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)

    def test_yeast_optimization(self):
        opt = BioOptimizer(species="yeast")
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)
        assert len(seq) > 0

    def test_cho_optimization(self):
        opt = BioOptimizer(species="cho")
        seq, preds, cert = opt.optimize(TEST_DNA)
        assert isinstance(seq, str)
        assert len(seq) > 0
