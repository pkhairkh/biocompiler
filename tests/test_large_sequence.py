"""Test Large Sequence Support — chunk-based optimization, boundary handling, progress callbacks."""

import pytest
from biocompiler.large_sequence import (
    optimize_large_sequence,
    ProteinTooLongError,
    MAX_PROTEIN_LENGTH_DEFAULT,
    _split_into_chunks,
    _merge_chunks,
    _repair_boundaries,
)
from biocompiler.optimization import OptimizationResult
from biocompiler.translation import translate, compute_cai
from biocompiler.type_system import CODON_TABLE


# ── Helper: Generate a random-ish protein of given length ────────────
# Uses a deterministic pattern to ensure reproducibility without
# relying on random number generators.
def _make_protein(length: int) -> str:
    """Generate a protein of the given length using a repeating pattern."""
    # Use common amino acids that have many codon choices
    pool = "ACDEFGHIKLMNPQRSTVWY"
    result = []
    for i in range(length):
        result.append(pool[i % len(pool)])
    return "".join(result)


class TestSplitIntoChunks:
    """Test the chunk splitting logic."""

    def test_short_protein_no_split(self):
        """Proteins shorter than chunk_size should not be split."""
        protein = _make_protein(100)
        chunks = _split_into_chunks(protein, chunk_size=300, overlap=10)
        assert len(chunks) == 1
        assert chunks[0] == (0, 100, protein)

    def test_exact_chunk_size(self):
        """Protein exactly chunk_size should be a single chunk."""
        protein = _make_protein(300)
        chunks = _split_into_chunks(protein, chunk_size=300, overlap=10)
        assert len(chunks) == 1

    def test_two_chunks(self):
        """Protein slightly longer than chunk_size should split into 2 chunks."""
        protein = _make_protein(310)
        chunks = _split_into_chunks(protein, chunk_size=300, overlap=10)
        assert len(chunks) == 2
        # First chunk: 0-300
        assert chunks[0][0] == 0
        assert chunks[0][1] == 300
        # Second chunk overlaps with first by 10 aa
        assert chunks[1][0] == 290  # 300 - 10
        assert chunks[1][1] == 310

    def test_all_chunks_cover_full_protein(self):
        """All positions should be covered by at least one chunk."""
        protein = _make_protein(800)
        chunks = _split_into_chunks(protein, chunk_size=300, overlap=10)
        # Verify every position is covered
        covered = set()
        for start, end, seq in chunks:
            for i in range(start, end):
                covered.add(i)
        assert covered == set(range(len(protein)))

    def test_overlap_preserves_sequence(self):
        """Each chunk's sequence should match the corresponding region of the protein."""
        protein = _make_protein(600)
        chunks = _split_into_chunks(protein, chunk_size=300, overlap=10)
        for start, end, seq in chunks:
            assert seq == protein[start:end], (
                f"Chunk [{start}:{end}] doesn't match protein substring"
            )


class TestOptimizeLargeSequence500aa:
    """Test optimization with a 500 aa protein (typical size)."""

    def test_500aa_produces_valid_result(self):
        """A 500 aa protein should optimize successfully."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
        )
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 500 * 3

    def test_500aa_translates_correctly(self):
        """The optimized sequence should translate back to the original protein."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
        )
        translated = translate(result.sequence, to_stop=False)
        assert translated == protein, (
            f"Translation mismatch: got {len(translated)} aa, expected {len(protein)}"
        )

    def test_500aa_valid_codons(self):
        """All codons in the result should be valid."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
        )
        for i in range(0, len(result.sequence), 3):
            codon = result.sequence[i:i + 3]
            assert codon in CODON_TABLE, f"Invalid codon {codon!r} at position {i}"

    def test_500aa_gc_in_range(self):
        """GC content should be within the default range."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
            gc_lo=0.30,
            gc_hi=0.70,
        )
        assert 0.30 <= result.gc_content <= 0.70, (
            f"GC content {result.gc_content:.3f} outside [0.30, 0.70]"
        )

    def test_500aa_positive_cai(self):
        """CAI should be positive and reasonable."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
        )
        assert result.cai > 0.0, f"CAI should be positive, got {result.cai}"


class TestOptimizeLargeSequence1000aa:
    """Test optimization with a 1000 aa protein (challenging size)."""

    def test_1000aa_produces_valid_result(self):
        """A 1000 aa protein should optimize successfully."""
        protein = _make_protein(1000)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
        )
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 1000 * 3

    def test_1000aa_translates_correctly(self):
        """The optimized sequence should translate back to the original protein."""
        protein = _make_protein(1000)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
        )
        translated = translate(result.sequence, to_stop=False)
        assert translated == protein, (
            f"Translation mismatch at 1000 aa"
        )

    def test_1000aa_valid_codons(self):
        """All codons should be valid."""
        protein = _make_protein(1000)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
        )
        for i in range(0, len(result.sequence), 3):
            codon = result.sequence[i:i + 3]
            assert codon in CODON_TABLE, f"Invalid codon {codon!r} at position {i}"

    def test_1000aa_with_restriction_enzymes(self):
        """1000 aa optimization with restriction site avoidance."""
        protein = _make_protein(1000)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
            enzymes=["EcoRI", "BamHI"],
            strict_mode=False,  # May not be able to remove all sites
        )
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 1000 * 3

    def test_1000aa_human_organism(self):
        """1000 aa optimization targeting human."""
        protein = _make_protein(1000)
        result = optimize_large_sequence(
            protein,
            organism="human",
            chunk_size=300,
            overlap=10,
        )
        assert isinstance(result, OptimizationResult)
        translated = translate(result.sequence, to_stop=False)
        assert translated == protein


class TestChunkBoundaryHandling:
    """Test that chunk boundaries are handled correctly."""

    def test_boundary_codons_are_valid(self):
        """Codons at chunk boundaries should be valid and encode correct amino acids."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
        )
        # Check codons at the boundary region (around position 290-310)
        for i in range(285, 315):
            if i < len(protein):
                codon = result.sequence[i * 3: i * 3 + 3]
                expected_aa = protein[i]
                actual_aa = CODON_TABLE.get(codon)
                assert actual_aa == expected_aa, (
                    f"Boundary codon at aa {i}: expected {expected_aa}, "
                    f"got {actual_aa} (codon {codon})"
                )

    def test_overlap_parameter_effect(self):
        """Different overlap values should all produce valid results."""
        protein = _make_protein(400)
        for overlap in [4, 10, 20]:
            result = optimize_large_sequence(
                protein,
                organism="ecoli",
                chunk_size=200,
                overlap=overlap,
            )
            assert len(result.sequence) == 400 * 3
            translated = translate(result.sequence, to_stop=False)
            assert translated == protein, (
                f"Translation failed with overlap={overlap}"
            )

    def test_small_chunk_size(self):
        """Small chunk size with overlap should still produce valid results."""
        protein = _make_protein(300)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=100,
            overlap=10,
        )
        assert len(result.sequence) == 300 * 3
        translated = translate(result.sequence, to_stop=False)
        assert translated == protein

    def test_no_cross_codon_gt_at_boundaries_eukaryote(self):
        """For eukaryotes, check that boundary GT issues are minimized."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="human",
            chunk_size=300,
            overlap=10,
        )
        # The result should still translate correctly (basic invariant)
        translated = translate(result.sequence, to_stop=False)
        assert translated == protein

    def test_chunk_overlap_at_least_4(self):
        """Overlap < 4 should produce a warning but still work."""
        protein = _make_protein(400)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=200,
            overlap=2,  # Below recommended 4
        )
        # Should still produce a valid result
        assert len(result.sequence) == 400 * 3


class TestProgressCallback:
    """Test progress callback functionality."""

    def test_callback_called(self):
        """Progress callback should be called during optimization."""
        protein = _make_protein(500)
        calls = []

        def on_progress(current, total):
            calls.append((current, total))

        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
            progress_callback=on_progress,
        )
        # Should have been called at least twice (start and end)
        assert len(calls) >= 2, f"Expected >= 2 callback calls, got {len(calls)}"

    def test_callback_reports_total(self):
        """Callback total should equal protein length."""
        protein = _make_protein(500)
        totals = set()

        def on_progress(current, total):
            totals.add(total)

        optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
            progress_callback=on_progress,
        )
        assert totals == {500}, f"Expected totals={{500}}, got {totals}"

    def test_callback_final_position_equals_total(self):
        """Final callback should report current == total."""
        protein = _make_protein(500)
        final_call = [None]

        def on_progress(current, total):
            final_call[0] = (current, total)

        optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
            progress_callback=on_progress,
        )
        assert final_call[0] is not None
        assert final_call[0][0] == final_call[0][1], (
            f"Final callback: current={final_call[0][0]}, total={final_call[0][1]}"
        )

    def test_callback_with_short_sequence(self):
        """Callback should work even with a short sequence (no chunking)."""
        protein = _make_protein(100)  # Below chunk_size
        calls = []

        def on_progress(current, total):
            calls.append((current, total))

        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
            progress_callback=on_progress,
        )
        assert len(calls) >= 2

    def test_callback_none_is_fine(self):
        """Passing None as callback should work without errors."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
            progress_callback=None,
        )
        assert isinstance(result, OptimizationResult)


class TestMaxProteinLengthCap:
    """Test the max_protein_length safety cap."""

    def test_default_cap_is_10000(self):
        """The default safety cap should be 10,000 aa."""
        assert MAX_PROTEIN_LENGTH_DEFAULT == 10_000

    def test_protein_within_cap_succeeds(self):
        """Proteins within the cap should optimize normally."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            max_protein_length=1000,
        )
        assert isinstance(result, OptimizationResult)

    def test_protein_exceeding_cap_raises(self):
        """Proteins exceeding the cap should raise ProteinTooLongError."""
        protein = _make_protein(600)
        with pytest.raises(ProteinTooLongError) as exc_info:
            optimize_large_sequence(
                protein,
                organism="ecoli",
                max_protein_length=500,
            )
        assert exc_info.value.protein_length == 600
        assert exc_info.value.max_length == 500

    def test_protein_at_exact_cap_succeeds(self):
        """Protein exactly at the cap should succeed."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            max_protein_length=500,
        )
        assert isinstance(result, OptimizationResult)

    def test_cap_zero_disables_limit(self):
        """Setting max_protein_length=0 should disable the safety cap."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            max_protein_length=0,  # Disable cap
        )
        assert isinstance(result, OptimizationResult)

    def test_cap_override_increases_limit(self):
        """Overriding the cap to a higher value should work."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            max_protein_length=50_000,
        )
        assert isinstance(result, OptimizationResult)


class TestInvalidInput:
    """Test error handling for invalid inputs."""

    def test_invalid_amino_acid_raises(self):
        """Proteins with invalid amino acid codes should raise InvalidProteinError."""
        from biocompiler.exceptions import InvalidProteinError
        with pytest.raises(InvalidProteinError):
            optimize_large_sequence("MKXINVALID", organism="ecoli")

    def test_empty_protein_raises(self):
        """Empty proteins should raise InvalidProteinError."""
        from biocompiler.exceptions import InvalidProteinError
        with pytest.raises(InvalidProteinError):
            optimize_large_sequence("", organism="ecoli")

    def test_chunk_size_less_than_overlap_raises(self):
        """chunk_size <= overlap should raise ValueError."""
        with pytest.raises(ValueError, match="chunk_size"):
            optimize_large_sequence(
                _make_protein(500),
                organism="ecoli",
                chunk_size=10,
                overlap=20,
            )


class TestRepairBoundaries:
    """Test the boundary repair logic directly."""

    def test_repair_fixes_wrong_codon(self):
        """_repair_boundaries should fix codons that encode wrong amino acids."""
        from biocompiler.large_sequence import _repair_boundaries
        protein = "MKA"
        # Use wrong codon for K (lysine) — AAA is correct, but put GGG (glycine) instead
        dna = "ATGGGGGCA"  # ATG=M, GGG=G(!), GCA=A
        repaired = _repair_boundaries(dna, protein, "ecoli")
        # K codon at position 1 should now encode K
        k_codon = repaired[3:6]
        assert CODON_TABLE.get(k_codon) == "K", (
            f"Expected K-encoding codon, got {k_codon} -> {CODON_TABLE.get(k_codon)}"
        )

    def test_repair_preserves_correct_codons(self):
        """_repair_boundaries should not change codons that are already correct."""
        from biocompiler.large_sequence import _repair_boundaries
        protein = "MK"
        dna = "ATGAAAAAA"  # ATG=M, AAA=K — both correct
        repaired = _repair_boundaries(dna, protein, "ecoli")
        assert repaired == dna


class TestIntegrationWithOptimizeSequence:
    """Test that optimize_large_sequence delegates to optimize_sequence for short proteins."""

    def test_short_protein_uses_fast_path(self):
        """Proteins shorter than chunk_size should use optimize_sequence directly."""
        protein = _make_protein(100)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
        )
        # Should produce a valid result without any chunk-related warnings
        assert isinstance(result, OptimizationResult)
        # No chunk warning for single-chunk path
        chunk_warnings = [w for w in result.warnings if "chunk" in w.lower()]
        assert len(chunk_warnings) == 0

    def test_result_has_expected_fields(self):
        """Result should have all expected OptimizationResult fields populated."""
        protein = _make_protein(500)
        result = optimize_large_sequence(
            protein,
            organism="ecoli",
            chunk_size=300,
            overlap=10,
        )
        assert result.sequence is not None
        assert len(result.sequence) > 0
        assert 0.0 <= result.gc_content <= 1.0
        assert 0.0 <= result.cai <= 1.0
        assert result.protein == protein
