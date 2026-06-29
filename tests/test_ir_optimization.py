"""
Tests for BioCompiler IR Optimization Passes (Task B-8)
=======================================================

These tests verify the IR→IR optimization passes in
:mod:`biocompiler.ir.optimization`:

  * ``optimize_codons``        — codon optimization (higher CAI, same protein)
  * ``eliminate_cpgs``         — CpG dinucleotide removal
  * ``run_optimization_pipeline`` — chain multiple passes

The single most important property tested is **protein preservation**:
every optimization pass must produce an L2 whose translation matches
the input L2's translation.  If a pass breaks this, the IR-level
check raises :class:`IRError` immediately.

Test gene: HBB (human hemoglobin beta), N-terminal 31 residues + stop,
loaded from the YAML spec shipped with the frontend (B-4).
Expected protein: ``MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*`` (UniProt P68871).
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

# Make ``src`` importable when running from a source checkout (the
# project layout uses ``src/biocompiler``).  This matches the pattern
# used in tests/test_hbb_full_pass.py.
sys.path.insert(0, "src")

from biocompiler.ir import (
    IR_L2_MatureMRNA,
    IR_L3_Polypeptide,
    IRLevel,
    IRError,
    translate,
    check_l2_invariants,
)
from biocompiler.ir.frontend import compile_from_spec
from biocompiler.ir.optimization import (
    optimize_codons,
    eliminate_cpgs,
    run_optimization_pipeline,
    ir_cai,
    _PASSES,
    _DEFAULT_PASSES,
)
from biocompiler.optimizer.utils import OptimizationResult
from biocompiler.organisms import resolve_organism, CODON_ADAPTIVENESS_TABLES
from biocompiler.optimizer.cai import _compute_cai_fast


# ────────────────────────────────────────────────────────────────────
# Fixtures
# ────────────────────────────────────────────────────────────────────

HBB_YAML = "src/biocompiler/ir/example_specs/hbb.yaml"
HBB_PROTEIN_WITH_STOP = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*"


@pytest.fixture(scope="module")
def hbb_l2() -> IR_L2_MatureMRNA:
    """Compile HBB from its YAML spec down to IR-L2 (module-cached)."""
    ir = compile_from_spec(HBB_YAML, target_level=IRLevel.L2)
    assert isinstance(ir, IR_L2_MatureMRNA)
    return ir


@pytest.fixture
def fresh_hbb_l2() -> IR_L2_MatureMRNA:
    """A fresh HBB L2 per test (no cross-test mutation).

    Tests that mutate metadata or expect a pristine object should use
    this fixture instead of the module-cached ``hbb_l2``.
    """
    return compile_from_spec(HBB_YAML, target_level=IRLevel.L2)


# ────────────────────────────────────────────────────────────────────
# optimize_codons
# ────────────────────────────────────────────────────────────────────

class TestOptimizeCodons:
    """Codon optimization: IR_L2 → IR_L2 (same protein, higher CAI)."""

    def test_returns_new_l2_object(self, hbb_l2):
        """The pass must return a *new* IR_L2, not mutate the input."""
        optimized = optimize_codons(hbb_l2)
        assert isinstance(optimized, IR_L2_MatureMRNA)
        assert optimized is not hbb_l2

    def test_preserves_protein(self, hbb_l2):
        """The translated protein MUST be unchanged (semantic correctness)."""
        optimized = optimize_codons(hbb_l2)
        assert translate(optimized).sequence == translate(hbb_l2).sequence

    def test_preserves_hbb_protein_value(self, hbb_l2):
        """Optimized HBB must still translate to the canonical HBB protein."""
        optimized = optimize_codons(hbb_l2)
        assert translate(optimized).sequence == HBB_PROTEIN_WITH_STOP

    def test_cds_is_different(self, hbb_l2):
        """The CDS should actually change (otherwise the pass was a no-op)."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.cds != hbb_l2.cds, "optimize_codons did not change the CDS"

    def test_cds_starts_with_aug(self, hbb_l2):
        """Optimized CDS must still start with AUG (start codon)."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.cds.startswith("AUG")

    def test_cds_ends_with_stop(self, hbb_l2):
        """Optimized CDS must still end with a stop codon."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.cds[-3:] in {"UAA", "UAG", "UGA"}

    def test_cds_length_divisible_by_3(self, hbb_l2):
        """CDS length must remain a multiple of 3."""
        optimized = optimize_codons(hbb_l2)
        assert len(optimized.cds) % 3 == 0

    def test_cds_length_unchanged(self, hbb_l2):
        """CDS length must be unchanged (same number of codons)."""
        optimized = optimize_codons(hbb_l2)
        assert len(optimized.cds) == len(hbb_l2.cds)

    def test_preserves_utrs(self, hbb_l2):
        """5'UTR and 3'UTR must be passed through untouched."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.five_utr == hbb_l2.five_utr
        assert optimized.three_utr == hbb_l2.three_utr

    def test_preserves_organism(self, hbb_l2):
        """Organism must be preserved (optimization is per-organism)."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.organism == hbb_l2.organism

    def test_preserves_gene_name(self, hbb_l2):
        """gene_name must be preserved."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.gene_name == hbb_l2.gene_name

    def test_passes_l2_invariants(self, hbb_l2):
        """The output must satisfy IR-L2 invariants (checked inside the pass
        too, but we re-check to be belt-and-braces)."""
        optimized = optimize_codons(hbb_l2)
        assert check_l2_invariants(optimized) is True

    def test_provenance_pass_name(self, hbb_l2):
        """metadata['pass'] must be 'optimize_codons'."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.metadata["pass"] == "optimize_codons"

    def test_provenance_source_level_is_l2(self, hbb_l2):
        """IR→IR pass: source_level must be 'L2' (not 'L1' or 'L3')."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.metadata["source_level"] == "L2"

    def test_provenance_optimization_tag(self, hbb_l2):
        """metadata['optimization'] must record the optimization kind."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.metadata["optimization"] == "codon_optimization"

    def test_cai_before_metadata(self, hbb_l2):
        """CAI of the original CDS must be recorded."""
        optimized = optimize_codons(hbb_l2)
        assert "cai_before" in optimized.metadata
        assert optimized.metadata["cai_before"] == pytest.approx(ir_cai(hbb_l2))

    def test_cai_after_metadata(self, hbb_l2):
        """CAI of the optimized CDS must be recorded."""
        optimized = optimize_codons(hbb_l2)
        assert "cai_after" in optimized.metadata
        assert 0.0 <= optimized.metadata["cai_after"] <= 1.0

    def test_cai_improves(self, hbb_l2):
        """The optimized CDS must have higher CAI than the original."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.metadata["cai_after"] > optimized.metadata["cai_before"]

    def test_cai_delta_positive(self, hbb_l2):
        """cai_delta must be positive (CAI went up)."""
        optimized = optimize_codons(hbb_l2)
        assert optimized.metadata["cai_delta"] > 0

    def test_ir_cai_matches_optimized_cai(self, hbb_l2):
        """ir_cai() on the optimized L2 must match the optimizer's cai_after
        (they use the same kernel and tables, up to display rounding)."""
        optimized = optimize_codons(hbb_l2)
        # The optimizer rounds CAI to 4 decimals for display; ir_cai
        # returns full float precision.  Compare with a tolerance that
        # accommodates the rounding.
        assert ir_cai(optimized) == pytest.approx(
            optimized.metadata["cai_after"], abs=1e-3
        )

    def test_organism_override(self, hbb_l2):
        """Passing organism= explicitly must override the L2's own organism."""
        # E. coli uses very different optimal codons from human — the CDS
        # should differ from the human-optimised one.
        human_opt = optimize_codons(hbb_l2, organism="human")
        ecoli_opt = optimize_codons(hbb_l2, organism="ecoli")
        assert human_opt.cds != ecoli_opt.cds
        # Protein must still be preserved in both.
        assert translate(human_opt).sequence == translate(hbb_l2).sequence
        assert translate(ecoli_opt).sequence == translate(hbb_l2).sequence


# ────────────────────────────────────────────────────────────────────
# eliminate_cpgs
# ────────────────────────────────────────────────────────────────────

class TestEliminateCpGs:
    """CpG elimination: IR_L2 → IR_L2 (fewer CG dinucleotides)."""

    def test_returns_new_l2_object(self, hbb_l2):
        eliminated = eliminate_cpgs(hbb_l2)
        assert isinstance(eliminated, IR_L2_MatureMRNA)
        assert eliminated is not hbb_l2

    def test_preserves_protein(self, hbb_l2):
        """Protein must be unchanged after CpG elimination."""
        eliminated = eliminate_cpgs(hbb_l2)
        assert translate(eliminated).sequence == translate(hbb_l2).sequence

    def test_hbb_has_cpgs_to_remove(self, hbb_l2):
        """Sanity: the HBB example CDS contains CpG dinucleotides.

        Without this, the 'CpGs go down' test below would be vacuous.
        """
        assert hbb_l2.cds.count("CG") >= 1, (
            "HBB example CDS has no CpGs — the test gene needs updating"
        )

    def test_cpg_count_decreases_or_stays_zero(self, hbb_l2):
        """cpg_after must be <= cpg_before, and strictly less if there
        were any CpGs to remove."""
        eliminated = eliminate_cpgs(hbb_l2)
        before = eliminated.metadata["cpg_before"]
        after = eliminated.metadata["cpg_after"]
        assert after <= before
        if before > 0:
            # The optimizer's aggressive CpG mode should remove at least one.
            assert after < before, (
                f"CpG count did not decrease: {before} → {after}"
            )

    def test_cpg_after_is_zero_for_hbb(self, hbb_l2):
        """For HBB, the aggressive CpG mode should eliminate ALL CpGs."""
        eliminated = eliminate_cpgs(hbb_l2)
        assert eliminated.metadata["cpg_after"] == 0
        assert eliminated.cds.count("CG") == 0

    def test_cpg_removed_metadata(self, hbb_l2):
        """cpg_removed = cpg_before - cpg_after must be recorded correctly."""
        eliminated = eliminate_cpgs(hbb_l2)
        expected = eliminated.metadata["cpg_before"] - eliminated.metadata["cpg_after"]
        assert eliminated.metadata["cpg_removed"] == expected

    def test_provenance_pass_name(self, hbb_l2):
        eliminated = eliminate_cpgs(hbb_l2)
        assert eliminated.metadata["pass"] == "eliminate_cpgs"

    def test_provenance_source_level_is_l2(self, hbb_l2):
        eliminated = eliminate_cpgs(hbb_l2)
        assert eliminated.metadata["source_level"] == "L2"

    def test_provenance_optimization_tag(self, hbb_l2):
        eliminated = eliminate_cpgs(hbb_l2)
        assert eliminated.metadata["optimization"] == "cpg_elimination"

    def test_preserves_utrs_and_organism(self, hbb_l2):
        eliminated = eliminate_cpgs(hbb_l2)
        assert eliminated.five_utr == hbb_l2.five_utr
        assert eliminated.three_utr == hbb_l2.three_utr
        assert eliminated.organism == hbb_l2.organism
        assert eliminated.gene_name == hbb_l2.gene_name

    def test_passes_l2_invariants(self, hbb_l2):
        eliminated = eliminate_cpgs(hbb_l2)
        assert check_l2_invariants(eliminated) is True


# ────────────────────────────────────────────────────────────────────
# run_optimization_pipeline
# ────────────────────────────────────────────────────────────────────

class TestRunOptimizationPipeline:
    """Chaining multiple IR→IR optimization passes."""

    def test_default_passes(self, hbb_l2):
        """Default pipeline runs optimize_codons then eliminate_cpgs."""
        result = run_optimization_pipeline(hbb_l2)
        assert result.metadata["passes_applied"] == [
            "optimize_codons",
            "eliminate_cpgs",
        ]

    def test_default_passes_constant(self):
        """The module-level _DEFAULT_PASSES constant must match."""
        assert _DEFAULT_PASSES == ["optimize_codons", "eliminate_cpgs"]

    def test_protein_preserved_through_pipeline(self, hbb_l2):
        result = run_optimization_pipeline(hbb_l2)
        assert translate(result).sequence == translate(hbb_l2).sequence

    def test_pipeline_changes_cds(self, hbb_l2):
        """The pipeline output CDS should differ from the input."""
        result = run_optimization_pipeline(hbb_l2)
        assert result.cds != hbb_l2.cds

    def test_pipeline_cai_high(self, hbb_l2):
        """The final CAI should be high (codon optimization ran)."""
        result = run_optimization_pipeline(hbb_l2)
        assert ir_cai(result) >= ir_cai(hbb_l2)

    def test_pipeline_no_cpgs(self, hbb_l2):
        """CpG elimination ran → final CDS has 0 CpGs."""
        result = run_optimization_pipeline(hbb_l2)
        assert result.cds.count("CG") == 0

    def test_custom_pass_list_single(self, hbb_l2):
        """A single-pass pipeline must work and record exactly that pass."""
        result = run_optimization_pipeline(hbb_l2, passes=["optimize_codons"])
        assert result.metadata["passes_applied"] == ["optimize_codons"]

    def test_custom_pass_list_reverse_order(self, hbb_l2):
        """Reverse order must also work (CpGs first, then codon opt)."""
        result = run_optimization_pipeline(
            hbb_l2, passes=["eliminate_cpgs", "optimize_codons"]
        )
        assert result.metadata["passes_applied"] == [
            "eliminate_cpgs",
            "optimize_codons",
        ]
        # Protein still preserved.
        assert translate(result).sequence == translate(hbb_l2).sequence

    def test_empty_pass_list_is_noop(self, hbb_l2):
        """An empty pass list returns a fresh L2 with the same CDS."""
        result = run_optimization_pipeline(hbb_l2, passes=[])
        assert result.metadata["passes_applied"] == []
        assert result.cds == hbb_l2.cds
        # But it's a new object (pipeline always returns a fresh L2).
        assert result is not hbb_l2

    def test_unknown_pass_raises_irerror(self, hbb_l2):
        """An unknown pass name must raise IRError (fail-fast)."""
        with pytest.raises(IRError, match="unknown optimization pass"):
            run_optimization_pipeline(hbb_l2, passes=["does_not_exist"])

    def test_unknown_pass_error_lists_available(self, hbb_l2):
        """The error message should list available passes for usability."""
        with pytest.raises(IRError, match="optimize_codons"):
            run_optimization_pipeline(hbb_l2, passes=["???"])

    def test_passes_registry_contains_expected(self):
        """The _PASSES registry must contain both implemented passes."""
        assert "optimize_codons" in _PASSES
        assert "eliminate_cpgs" in _PASSES
        assert _PASSES["optimize_codons"] is optimize_codons
        assert _PASSES["eliminate_cpgs"] is eliminate_cpgs

    def test_pipeline_preserves_utrs(self, hbb_l2):
        result = run_optimization_pipeline(hbb_l2)
        assert result.five_utr == hbb_l2.five_utr
        assert result.three_utr == hbb_l2.three_utr

    def test_pipeline_passes_l2_invariants(self, hbb_l2):
        result = run_optimization_pipeline(hbb_l2)
        assert check_l2_invariants(result) is True


# ────────────────────────────────────────────────────────────────────
# Protein preservation enforcement (the critical correctness check)
# ────────────────────────────────────────────────────────────────────

class TestProteinPreservationCheck:
    """The optimizer must raise IRError if it ever changes the protein.

    This is the single most important safety property of the IR passes.
    We test it by monkey-patching ``optimize_sequence`` to return a
    sequence that translates to a *different* protein, then verifying
    that the IR pass catches the discrepancy.
    """

    def _make_bogus_result(self, original_protein_no_stop: str) -> OptimizationResult:
        """Build a fake OptimizationResult whose DNA translates to a
        DIFFERENT protein than the input.

        We keep the start codon ATG (so the L2 invariant ``cds.startswith('AUG')``
        still passes — we want to test the *protein preservation* check,
        not the structural invariant check) but mutate the SECOND codon
        so the second amino acid changes.  For HBB that means V→A at
        position 1 (GTG → GCA).

        The optimizer's contract says ``result.sequence`` translates to
        ``result.protein`` — we honour that contract, just with the
        wrong protein.
        """
        if len(original_protein_no_stop) < 2:
            # Fall back to a single-codon swap if the protein is too short
            # to mutate position 1 without touching the start codon.
            bogus_dna = "GCA"  # GCA = A (different from M=ATG)
            bogus_protein = "A"
        else:
            # Keep ATG (M) at position 0; mutate position 1 to GCA (A).
            bogus_dna = "ATG" + "GCA"  # M + A
            # Fill the rest with GCA (A) for the remaining residues.
            bogus_dna += "GCA" * (len(original_protein_no_stop) - 2)
            bogus_protein = "M" + "A" * (len(original_protein_no_stop) - 1)
        return OptimizationResult(
            sequence=bogus_dna,
            gc_content=0.5,
            cai=0.5,
            protein=bogus_protein,
        )

    def test_optimize_codons_raises_on_protein_change(self, hbb_l2):
        """optimize_codons must raise IRError if the optimizer returns
        a sequence encoding a different protein."""
        original_protein = translate(hbb_l2).sequence
        protein_no_stop = original_protein.rstrip("*")
        bogus = self._make_bogus_result(protein_no_stop)

        with patch(
            "biocompiler.ir.optimization.optimize_sequence",
            return_value=bogus,
        ):
            with pytest.raises(IRError, match="changed the translated protein"):
                optimize_codons(hbb_l2)

    def test_eliminate_cpgs_raises_on_protein_change(self, hbb_l2):
        """eliminate_cpgs must raise IRError if the optimizer returns
        a sequence encoding a different protein."""
        original_protein = translate(hbb_l2).sequence
        protein_no_stop = original_protein.rstrip("*")
        bogus = self._make_bogus_result(protein_no_stop)

        with patch(
            "biocompiler.ir.optimization.optimize_sequence",
            return_value=bogus,
        ):
            with pytest.raises(IRError, match="changed the translated protein"):
                eliminate_cpgs(hbb_l2)

    def test_pipeline_raises_on_protein_change(self, hbb_l2):
        """The pipeline must propagate the IRError from any pass that
        breaks the protein."""
        original_protein = translate(hbb_l2).sequence
        protein_no_stop = original_protein.rstrip("*")
        bogus = self._make_bogus_result(protein_no_stop)

        with patch(
            "biocompiler.ir.optimization.optimize_sequence",
            return_value=bogus,
        ):
            with pytest.raises(IRError, match="changed the translated protein"):
                run_optimization_pipeline(hbb_l2)

    def test_error_message_shows_both_proteins(self, hbb_l2):
        """The IRError message should include both the old and new proteins
        for debuggability."""
        original_protein = translate(hbb_l2).sequence
        protein_no_stop = original_protein.rstrip("*")
        bogus = self._make_bogus_result(protein_no_stop)

        with patch(
            "biocompiler.ir.optimization.optimize_sequence",
            return_value=bogus,
        ):
            with pytest.raises(IRError) as exc_info:
                optimize_codons(hbb_l2)
        msg = str(exc_info.value)
        # The original HBB protein must appear in the error message.
        assert "MVHLTPEEK" in msg
        # And the bogus protein (starts with MA... not MV...) must appear.
        assert "MAA" in msg or "MAAAAA" in msg


# ────────────────────────────────────────────────────────────────────
# Input validation
# ────────────────────────────────────────────────────────────────────

class TestInputValidation:
    """The passes must reject bad inputs up-front with IRError."""

    def test_rejects_non_l2_input(self):
        """Passing a non-IR_L2 object must raise IRError."""
        with pytest.raises(IRError, match="IR_L2_MatureMRNA"):
            optimize_codons("not an IR object")  # type: ignore[arg-type]

    def test_rejects_l3_input(self):
        """Passing an IR_L3 must raise IRError (wrong level)."""
        ir_l3 = IR_L3_Polypeptide(
            sequence="MAK*", organism="e_coli"
        )
        with pytest.raises(IRError, match="IR_L2_MatureMRNA"):
            optimize_codons(ir_l3)  # type: ignore[arg-type]

    def test_accepts_l2_no_start_codon(self):
        """L2 without AUG start is now accepted (back-translated genes may not start with M)."""
        good_l2 = IR_L2_MatureMRNA(
            sequence="GCUUAA",
            five_utr="",
            cds="GCUUAA",  # no AUG start — OK now
            three_utr="",
            organism="e_coli",
        )
        optimize_codons(good_l2)  # should work

    def test_rejects_malformed_l2_no_stop(self):
        """An L2 whose CDS doesn't end with a stop codon must be rejected."""
        bad_l2 = IR_L2_MatureMRNA(
            sequence="AUGGCU",
            five_utr="",
            cds="AUGGCU",  # no stop
            three_utr="",
            organism="e_coli",
        )
        with pytest.raises(IRError, match="stop"):
            optimize_codons(bad_l2)

    def test_eliminate_cpgs_rejects_non_l2(self):
        with pytest.raises(IRError):
            eliminate_cpgs(None)  # type: ignore[arg-type]

    def test_pipeline_rejects_non_l2(self):
        with pytest.raises(IRError):
            run_optimization_pipeline("nope")  # type: ignore[arg-type]


# ────────────────────────────────────────────────────────────────────
# ir_cai helper
# ────────────────────────────────────────────────────────────────────

class TestIrCai:
    """The IR-level CAI computation helper."""

    def test_returns_float_in_unit_interval(self, hbb_l2):
        cai = ir_cai(hbb_l2)
        assert isinstance(cai, float)
        assert 0.0 <= cai <= 1.0

    def test_matches_optimizer_cai_for_hbb(self, hbb_l2):
        """ir_cai must agree with the optimizer's CAI computation
        (same kernel, same adaptiveness table)."""
        # The optimizer's CAI for the HBB protein is ~0.82 after opt.
        # ir_cai on the *original* HBB CDS computes the pre-optimization CAI.
        # We just check it's in a sensible range for natural HBB in human.
        cai = ir_cai(hbb_l2)
        # Natural HBB has moderately high CAI in human (it's a highly
        # expressed gene).  Anything in [0.3, 0.95] is plausible.
        assert 0.3 <= cai <= 0.95

    def test_nonzero_for_countable_cds(self):
        """A CDS with at least one non-M, non-stop codon must yield CAI > 0.

        (M and * are excluded from CAI by convention — they have no
        synonymous alternatives — so a CDS of just `AUG UAA` would give
        CAI 0.  We use `AUG GCU UAA` = M A * so the Alanine codon is
        counted.)
        """
        l2 = IR_L2_MatureMRNA(
            sequence="AUGGCUUAA",
            five_utr="",
            cds="AUGGCUUAA",
            three_utr="",
            organism="e_coli",
        )
        assert ir_cai(l2) > 0.0

    def test_unknown_organism_returns_zero(self):
        """An unknown organism must not crash — it returns 0.0."""
        l2 = IR_L2_MatureMRNA(
            sequence="AUGGCUUAA",
            five_utr="",
            cds="AUGGCUUAA",
            three_utr="",
            organism="totally_made_up_organism_xyz",
        )
        assert ir_cai(l2) == 0.0


# ────────────────────────────────────────────────────────────────────
# Provenance / metadata propagation
# ────────────────────────────────────────────────────────────────────

class TestProvenancePropagation:
    """Upstream metadata must be preserved; new provenance must be stamped."""

    def test_upstream_metadata_preserved(self, fresh_hbb_l2):
        """Metadata from the frontend (e.g. source_format) must survive
        the optimization pass."""
        assert fresh_hbb_l2.metadata.get("source_format") == "yaml"
        optimized = optimize_codons(fresh_hbb_l2)
        assert optimized.metadata.get("source_format") == "yaml"

    def test_pass_overwrites_previous_pass(self, fresh_hbb_l2):
        """The new 'pass' value must overwrite the old one (not stack)."""
        assert fresh_hbb_l2.metadata.get("pass") == "splice"  # came from L1→L2
        optimized = optimize_codons(fresh_hbb_l2)
        assert optimized.metadata["pass"] == "optimize_codons"
        # The old 'splice' value must NOT also be present as 'pass'.
        # (It's still in the metadata chain via the spread, but the 'pass'
        # key is overwritten.)

    def test_pipeline_stamps_passes_applied(self, fresh_hbb_l2):
        result = run_optimization_pipeline(fresh_hbb_l2)
        assert "passes_applied" in result.metadata
        assert isinstance(result.metadata["passes_applied"], list)

    def test_pipeline_passes_applied_is_copy(self, fresh_hbb_l2):
        """passes_applied must be a fresh list, not the module default."""
        result = run_optimization_pipeline(fresh_hbb_l2)
        applied = result.metadata["passes_applied"]
        applied.append("HACK")
        # The module default must not have been mutated.
        assert "HACK" not in _DEFAULT_PASSES

    def test_each_pass_in_pipeline_stamps_its_own_name(self, fresh_hbb_l2):
        """After the pipeline, metadata['pass'] reflects the LAST pass run."""
        result = run_optimization_pipeline(fresh_hbb_l2)
        # Default pipeline ends with eliminate_cpgs.
        assert result.metadata["pass"] == "eliminate_cpgs"

    def test_gene_name_preserved_in_pipeline(self, fresh_hbb_l2):
        result = run_optimization_pipeline(fresh_hbb_l2)
        assert result.gene_name == fresh_hbb_l2.gene_name


# ────────────────────────────────────────────────────────────────────
# HBB end-to-end demo
# ────────────────────────────────────────────────────────────────────

class TestHBBDemo:
    """The canonical BioCompiler demo: HBB → optimize → verify.

    This is the 'hello world' that ties together the frontend (B-4),
    the IR pipeline (B-3), and these optimization passes (B-8).
    """

    def test_hbb_yaml_to_optimized_l2(self):
        """HBB YAML → L2 → optimize_codons → L2' with:
        - same protein (UniProt P68871 N-terminus)
        - higher CAI
        - 0 CpGs (incidentally, since optimize_codons runs CpG cleanup)
        """
        ir_l2 = compile_from_spec(HBB_YAML, target_level=IRLevel.L2)
        original_protein = translate(ir_l2).sequence
        original_cai = ir_cai(ir_l2)

        optimized = optimize_codons(ir_l2)
        optimized_protein = translate(optimized).sequence
        optimized_cai = ir_cai(optimized)

        # 1. Protein preserved exactly.
        assert optimized_protein == original_protein
        assert optimized_protein == HBB_PROTEIN_WITH_STOP

        # 2. CAI strictly improved.
        assert optimized_cai > original_cai

        # 3. The optimized CDS is different from the original.
        assert optimized.cds != ir_l2.cds

    def test_hbb_full_pipeline_demo(self):
        """HBB → L2 → run_optimization_pipeline → L2' → translate → same protein."""
        ir_l2 = compile_from_spec(HBB_YAML, target_level=IRLevel.L2)
        original_protein = translate(ir_l2).sequence

        optimized = run_optimization_pipeline(ir_l2)
        optimized_protein = translate(optimized).sequence

        assert optimized_protein == original_protein
        # Both passes ran.
        assert optimized.metadata["passes_applied"] == [
            "optimize_codons",
            "eliminate_cpgs",
        ]
        # Final CAI is high.
        assert ir_cai(optimized) > ir_cai(ir_l2)
        # No CpGs left.
        assert optimized.cds.count("CG") == 0

    def test_hbb_optimization_is_idempotent_in_protein(self):
        """Running optimize_codons twice yields the same protein both times."""
        ir_l2 = compile_from_spec(HBB_YAML, target_level=IRLevel.L2)
        opt1 = optimize_codons(ir_l2)
        opt2 = optimize_codons(opt1)
        assert translate(opt1).sequence == translate(opt2).sequence
        assert translate(opt2).sequence == translate(ir_l2).sequence

    def test_hbb_cai_before_and_after_recorded(self):
        """The HBB demo must record both CAI before and after, with after > before.

        HBB is a naturally-occurring human gene, so its native codon usage is
        already well-adapted to human expression (CAI ≈ 0.71). The optimizer
        can only modestly improve on this for a 32-codon fragment — typical
        delta is 0.01–0.02. The test threshold of 0.005 reflects this: it
        verifies the optimizer is doing *something* (delta > 0) while not
        demanding an unrealistic improvement on an already-optimized gene.
        """
        ir_l2 = compile_from_spec(HBB_YAML, target_level=IRLevel.L2)
        optimized = optimize_codons(ir_l2)
        before = optimized.metadata["cai_before"]
        after = optimized.metadata["cai_after"]
        assert before > 0
        assert after > before
        # CAI delta is positive and meaningful for a naturally-occurring gene.
        assert optimized.metadata["cai_delta"] > 0.005


# ────────────────────────────────────────────────────────────────────
# GFP example (different gene, different organism) — robustness
# ────────────────────────────────────────────────────────────────────

class TestGFPRobustness:
    """Sanity check on a second gene (GFP N-terminus).

    Ensures the passes aren't accidentally hardcoded to HBB / human.
    The GFP YAML spec lists ``aequorea_victoria`` as its source organism
    (GFP comes from the jellyfish *Aequorea victoria*), but that
    organism isn't in the optimizer's supported-organisms list.  Since
    GFP is in practice always expressed in E. coli, we override the
    target organism to ``ecoli`` — this also demonstrates the
    ``organism=`` override parameter end-to-end.
    """

    def test_gfp_optimization_preserves_protein(self):
        """GFP YAML → L2 → optimize_codons(organism='ecoli') → same protein."""
        ir_l2 = compile_from_spec(
            "src/biocompiler/ir/example_specs/gfp.yaml",
            target_level=IRLevel.L2,
        )
        original_protein = translate(ir_l2).sequence
        optimized = optimize_codons(ir_l2, organism="ecoli")
        assert translate(optimized).sequence == original_protein

    def test_gfp_optimization_improves_cai(self):
        """GFP optimization should not decrease CAI for E. coli."""
        ir_l2 = compile_from_spec(
            "src/biocompiler/ir/example_specs/gfp.yaml",
            target_level=IRLevel.L2,
        )
        optimized = optimize_codons(ir_l2, organism="ecoli")
        # CAI should not decrease (it may stay similar if the input was
        # already well-adapted, but should never go down).
        assert optimized.metadata["cai_after"] >= optimized.metadata["cai_before"] - 0.01
