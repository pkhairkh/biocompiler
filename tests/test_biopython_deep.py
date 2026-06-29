"""Tests for biocompiler.biopython_compat — Deep BioPython integration.

Covered areas:
1. CodonUsageTable integration (load_codon_usage_table, compute_cai_from_table)
2. Pairwise alignment (align_to_reference)
3. Phylogenetic distance (phylo_distance)
4. ORF detection (detect_orfs)
5. Local BLAST (blast_local) — mocked since BLAST+ is rarely installed
6. Back-translation (back_translate_protein)
7. Graceful degradation when BioPython is not installed

All BioPython-dependent tests use ``pytest.importorskip("Bio")`` so the
entire suite is skipped cleanly when BioPython is absent.
"""

from __future__ import annotations

import sys
import types
from unittest import mock

import pytest

# Mark every test in this module as requiring an external tool (BioPython).
# Functional tests below use ``pytest.importorskip("Bio")``; the marker
# keeps them deselected by default alongside other requires_external tests.
pytestmark = pytest.mark.requires_external


# ═══════════════════════════════════════════════════════════════════════
# 0. Module import / new public names
# ═══════════════════════════════════════════════════════════════════════


class TestDeepModuleImport:
    """Verify the new deep-integration names are exposed."""

    def test_module_importable(self):
        """biocompiler.biopython_compat can still be imported without BioPython."""
        import biocompiler.shared.biopython_compat as mod
        assert mod is not None

    def test_new_public_names_exist(self):
        """New deep-integration functions are in __all__."""
        import biocompiler.shared.biopython_compat as mod

        expected = [
            "CodonUsageResult", "load_codon_usage_table", "compute_cai_from_table",
            "AlignmentResult", "align_to_reference",
            "phylo_distance",
            "ORFResult", "detect_orfs",
            "BlastResult", "blast_local",
            "back_translate_protein",
        ]
        for name in expected:
            assert name in mod.__all__, f"{name} not in __all__"

    def test_new_functions_are_callable(self):
        """New functions are callable."""
        import biocompiler.shared.biopython_compat as mod

        for name in ["load_codon_usage_table", "compute_cai_from_table",
                      "align_to_reference", "phylo_distance",
                      "detect_orfs", "blast_local", "back_translate_protein"]:
            assert callable(getattr(mod, name)), f"{name} is not callable"

    def test_new_dataclasses_importable(self):
        """New dataclass types can be accessed."""
        import biocompiler.shared.biopython_compat as mod

        assert mod.CodonUsageResult is not None
        assert mod.AlignmentResult is not None
        assert mod.ORFResult is not None
        assert mod.BlastResult is not None


class TestDeepGracefulDegradation:
    """All new deep-integration functions must raise ImportError when
    BioPython is unavailable."""

    @staticmethod
    def _make_biopython_unavailable():
        """Return a context manager that makes ``import Bio`` fail."""
        import builtins
        real_import = builtins.__import__

        def _fake_import(name, *args, **kwargs):
            if name == "Bio" or name.startswith("Bio."):
                raise ImportError(f"No module named '{name}'")
            return real_import(name, *args, **kwargs)

        return mock.patch("builtins.__import__", _fake_import)

    def test_load_codon_usage_table_raises_without_biopython(self):
        import biocompiler.shared.biopython_compat as mod
        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                mod.load_codon_usage_table()

    def test_align_to_reference_raises_without_biopython(self):
        import biocompiler.shared.biopython_compat as mod
        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                mod.align_to_reference("ATG", "ATG")

    def test_phylo_distance_raises_without_biopython(self):
        import biocompiler.shared.biopython_compat as mod
        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                mod.phylo_distance("ATGCGT")

    def test_detect_orfs_raises_without_biopython(self):
        import biocompiler.shared.biopython_compat as mod
        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                mod.detect_orfs("ATGCGT")

    def test_blast_local_raises_without_biopython(self):
        import biocompiler.shared.biopython_compat as mod
        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                mod.blast_local("ATGCGT", "/fake/db")

    def test_back_translate_protein_raises_without_biopython(self):
        import biocompiler.shared.biopython_compat as mod
        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                mod.back_translate_protein("MVSKGE")


# ═══════════════════════════════════════════════════════════════════════
# Functional tests (require BioPython)
# ═══════════════════════════════════════════════════════════════════════

Bio = pytest.importorskip("Bio", reason="BioPython not installed — skipping deep integration tests")


class TestCodonUsageTable:
    """Test CodonUsageTable integration."""

    def test_load_default_table(self):
        """load_codon_usage_table returns a CodonUsageResult for default organism."""
        from biocompiler.infrastructure.biopython_compat import load_codon_usage_table, CodonUsageResult

        result = load_codon_usage_table()
        assert isinstance(result, CodonUsageResult)
        assert result.organism == "Homo_sapiens"
        assert len(result.codon_counts) > 0
        assert len(result.adaptiveness) > 0

    def test_load_ecoli_table(self):
        """load_codon_usage_table works for E. coli."""
        from biocompiler.infrastructure.biopython_compat import load_codon_usage_table

        result = load_codon_usage_table(organism="Escherichia_coli")
        assert "ATG" in result.codon_counts
        assert result.adaptiveness.get("ATG", 0) > 0

    def test_codon_counts_positive(self):
        """Codon counts are positive integers."""
        from biocompiler.infrastructure.biopython_compat import load_codon_usage_table

        result = load_codon_usage_table()
        for codon, count in result.codon_counts.items():
            assert isinstance(count, int)
            assert count >= 0

    def test_adaptiveness_range(self):
        """Adaptiveness values are in [0, 1]."""
        from biocompiler.infrastructure.biopython_compat import load_codon_usage_table

        result = load_codon_usage_table()
        for codon, w in result.adaptiveness.items():
            assert 0.0 <= w <= 1.0, f"Adaptiveness for {codon} = {w} out of range"

    def test_source_description(self):
        """Source is descriptive."""
        from biocompiler.infrastructure.biopython_compat import load_codon_usage_table

        result = load_codon_usage_table()
        assert "BioCompiler" in result.source or "BioPython" in result.source

    def test_nonexistent_fasta_raises(self):
        """Providing a non-existent FASTA path raises FileNotFoundError."""
        from biocompiler.infrastructure.biopython_compat import load_codon_usage_table

        with pytest.raises(FileNotFoundError):
            load_codon_usage_table(fasta_path="/nonexistent/path.fasta")


class TestComputeCAIFromTable:
    """Test compute_cai_from_table function."""

    def test_cai_from_default_table(self, sample_dna):
        """CAI computed from table matches expected range."""
        from biocompiler.infrastructure.biopython_compat import load_codon_usage_table, compute_cai_from_table

        table = load_codon_usage_table()
        cai = compute_cai_from_table(sample_dna, table)
        assert isinstance(cai, float)
        assert 0.0 <= cai <= 1.0

    def test_cai_empty_sequence(self):
        """Empty sequence returns 0.0."""
        from biocompiler.infrastructure.biopython_compat import load_codon_usage_table, compute_cai_from_table

        table = load_codon_usage_table()
        assert compute_cai_from_table("", table) == 0.0


class TestAlignToReference:
    """Test align_to_reference function."""

    def test_identical_sequences_global(self):
        """Global alignment of identical sequences has identity 1.0."""
        from biocompiler.infrastructure.biopython_compat import align_to_reference

        seq = "ATGGTTTCTAAAGGTGAA"
        result = align_to_reference(seq, seq, mode="global")
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.gaps == 0
        assert result.score > 0

    def test_identical_sequences_local(self):
        """Local alignment of identical sequences has identity 1.0."""
        from biocompiler.infrastructure.biopython_compat import align_to_reference

        seq = "ATGGTTTCTAAAGGTGAA"
        result = align_to_reference(seq, seq, mode="local")
        assert result.identity == 1.0

    def test_different_sequences(self):
        """Alignment of different sequences has identity < 1.0."""
        from biocompiler.infrastructure.biopython_compat import align_to_reference

        result = align_to_reference("ATGCGTCGA", "ATTCGTCGA", mode="global")
        assert result.identity < 1.0
        assert result.mismatches >= 1

    def test_invalid_mode_raises(self):
        """Invalid mode raises ValueError."""
        from biocompiler.infrastructure.biopython_compat import align_to_reference

        with pytest.raises(ValueError, match="mode"):
            align_to_reference("ATG", "ATG", mode="semi-global")

    def test_result_has_algorithm(self):
        """Result includes the algorithm name."""
        from biocompiler.infrastructure.biopython_compat import align_to_reference

        result = align_to_reference("ATGCGT", "ATGCGT")
        assert "global" in result.algorithm.lower() or "PairwiseAligner" in result.algorithm


class TestPhyloDistance:
    """Test phylo_distance function."""

    def test_distance_nonnegative(self, sample_dna):
        """Phylogenetic distance is non-negative."""
        from biocompiler.infrastructure.biopython_compat import phylo_distance

        dist = phylo_distance(sample_dna)
        assert dist >= 0.0

    def test_euclidean_vs_cosine(self, sample_dna):
        """Both methods return valid distances."""
        from biocompiler.infrastructure.biopython_compat import phylo_distance

        euc = phylo_distance(sample_dna, method="euclidean")
        cos = phylo_distance(sample_dna, method="cosine")
        assert euc >= 0.0
        assert 0.0 <= cos <= 2.0

    def test_different_organisms_different_distances(self, sample_dna):
        """Different organisms may yield different distances."""
        from biocompiler.infrastructure.biopython_compat import phylo_distance

        dist_human = phylo_distance(sample_dna, organism="Homo_sapiens")
        dist_ecoli = phylo_distance(sample_dna, organism="Escherichia_coli")
        # They may or may not differ for this short sequence, but both should be valid
        assert dist_human >= 0.0
        assert dist_ecoli >= 0.0

    def test_empty_sequence_large_distance(self):
        """Empty sequence returns large distance."""
        from biocompiler.infrastructure.biopython_compat import phylo_distance

        euc = phylo_distance("", method="euclidean")
        assert euc == float('inf')

    def test_invalid_method_raises(self, sample_dna):
        """Invalid method raises ValueError."""
        from biocompiler.infrastructure.biopython_compat import phylo_distance

        with pytest.raises(ValueError, match="method"):
            phylo_distance(sample_dna, method="manhattan")


class TestDetectORFs:
    """Test detect_orfs function."""

    def test_detect_orfs_in_coding_sequence(self):
        """detect_orfs finds ORFs in a known coding sequence."""
        from biocompiler.infrastructure.biopython_compat import detect_orfs

        # eGFP is 720bp, 239 aa — should find at least the main ORF
        egfp = (
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
        orfs = detect_orfs(egfp, min_length_aa=50)
        assert len(orfs) >= 1
        # The main ORF should start at position 0
        main_orfs = [o for o in orfs if o.start == 0]
        assert len(main_orfs) >= 1
        assert main_orfs[0].length_aa >= 200

    def test_orf_result_fields(self):
        """ORFResult has all expected fields."""
        from biocompiler.infrastructure.biopython_compat import detect_orfs, ORFResult

        seq = "ATG" + "GCT" * 40 + "TAA"  # 123 bp, 40 aa ORF
        orfs = detect_orfs(seq, min_length_aa=30)
        assert len(orfs) >= 1
        orf = orfs[0]
        assert isinstance(orf, ORFResult)
        assert orf.start >= 0
        assert orf.end > orf.start
        assert orf.strand in (1, -1)
        assert orf.frame in (0, 1, 2)
        assert len(orf.protein) > 0
        assert orf.length_aa >= 30

    def test_no_orfs_in_short_sequence(self):
        """Very short sequences have no ORFs above threshold."""
        from biocompiler.infrastructure.biopython_compat import detect_orfs

        orfs = detect_orfs("ATGCGT", min_length_aa=30)
        assert len(orfs) == 0

    def test_empty_sequence_no_orfs(self):
        """Empty sequence returns no ORFs."""
        from biocompiler.infrastructure.biopython_compat import detect_orfs

        orfs = detect_orfs("")
        assert len(orfs) == 0


class TestBlastLocal:
    """Test blast_local function — mostly mocking since BLAST+ is rarely available."""

    def test_blast_not_available_raises(self):
        """blast_local raises FileNotFoundError when BLAST+ is not on PATH."""
        from biocompiler.infrastructure.biopython_compat import blast_local

        with pytest.raises(FileNotFoundError, match="BLAST"):
            blast_local("ATGCGT", "/nonexistent/db", blast_exe="/nonexistent/blastn")


class TestBackTranslateProtein:
    """Test back_translate_protein function."""

    def test_back_translate_simple(self):
        """Back-translation of a simple protein produces valid DNA."""
        from biocompiler.infrastructure.biopython_compat import back_translate_protein

        dna = back_translate_protein("MVSKGE")
        assert len(dna) == 18  # 6 amino acids × 3 bases
        assert all(b in "ACGT" for b in dna)

    def test_back_translate_empty(self):
        """Empty protein returns empty string."""
        from biocompiler.infrastructure.biopython_compat import back_translate_protein

        assert back_translate_protein("") == ""

    def test_back_translate_preserves_translation(self):
        """Back-translated DNA translates back to the same protein."""
        from biocompiler.infrastructure.biopython_compat import back_translate_protein
        from biocompiler.expression.translation import translate

        protein = "MVSKGE"
        dna = back_translate_protein(protein)
        retranslated = translate(dna)
        assert retranslated == protein

    def test_back_translate_organism_specific(self):
        """Different organisms may produce different DNA."""
        from biocompiler.infrastructure.biopython_compat import back_translate_protein

        dna_human = back_translate_protein("MVSKGE", organism="Homo_sapiens")
        dna_ecoli = back_translate_protein("MVSKGE", organism="Escherichia_coli")
        # Both should be valid 18-bp DNA sequences
        assert len(dna_human) == 18
        assert len(dna_ecoli) == 18

    def test_back_translate_with_stop(self):
        """Stop codon (*) is handled correctly."""
        from biocompiler.infrastructure.biopython_compat import back_translate_protein

        dna = back_translate_protein("MVSKGE*")
        assert len(dna) == 21  # 7 × 3
        # The last codon should be a stop codon
        from biocompiler.shared.constants import CODON_TABLE
        last_codon = dna[-3:]
        assert CODON_TABLE.get(last_codon) == "*"

    def test_invalid_translation_table_raises(self):
        """Invalid NCBI table number raises ValueError."""
        from biocompiler.infrastructure.biopython_compat import back_translate_protein

        with pytest.raises(ValueError, match="translation table"):
            back_translate_protein("M", table=999)


class TestDataclassInstantiation:
    """Verify dataclass constructors work correctly."""

    def test_codon_usage_result(self):
        from biocompiler.infrastructure.biopython_compat import CodonUsageResult

        r = CodonUsageResult(
            organism="test",
            codon_counts={"ATG": 100},
            codon_frequencies={"ATG": 1.0},
            adaptiveness={"ATG": 1.0},
            amino_acid_counts={"M": 100},
            source="test",
        )
        assert r.organism == "test"

    def test_alignment_result(self):
        from biocompiler.infrastructure.biopython_compat import AlignmentResult

        r = AlignmentResult(
            score=10.0,
            aligned_query="ATG",
            aligned_reference="ATG",
            identity=1.0,
            mismatches=0,
            gaps=0,
            algorithm="test",
        )
        assert r.score == 10.0

    def test_orf_result(self):
        from biocompiler.infrastructure.biopython_compat import ORFResult

        r = ORFResult(start=0, end=100, strand=1, frame=0, protein="M", length_aa=1)
        assert r.start == 0

    def test_blast_result(self):
        from biocompiler.infrastructure.biopython_compat import BlastResult

        r = BlastResult(
            query_id="q", subject_id="s",
            identity_percent=99.0, alignment_length=100,
            mismatches=1, gap_openings=0,
            query_start=1, query_end=100,
            subject_start=1, subject_end=100,
            e_value=1e-50, bit_score=150.0,
            tool="blastn",
        )
        assert r.identity_percent == 99.0
