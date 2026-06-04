"""Tests for ConstraintSpec.check() alignment with HardConstraint.check().

Ensures that ConstraintSpec.check() and the corresponding HardConstraint.check()
ALWAYS agree on the same input sequences.  This is the core invariant that
allows post-solve verification via ConstraintSpec to trust the same logic
that the optimizer used.

Task F1.7: Constraint enforcement fix.
"""

from __future__ import annotations

import pytest

from biocompiler.solver.types import ConstraintSpec, ConstraintType, ConstraintStrictness
from biocompiler.solver.constraints import (
    GCRangeConstraint,
    NoCrypticSpliceConstraint,
    NoCpGIslandConstraint,
    NoRestrictionSiteConstraint,
    NoATTTAMotifConstraint,
    NoTRunConstraint,
    TranslationConstraint,
)


# ── Helper: build a ConstraintSpec from a HardConstraint ──────────────

def _spec_from_gc(gc_lo: float = 0.30, gc_hi: float = 0.70) -> ConstraintSpec:
    return ConstraintSpec(
        ctype=ConstraintType.GC_CONTENT,
        name="GCRangeConstraint",
        params={"gc_lo": gc_lo, "gc_hi": gc_hi},
    )


def _spec_from_cpg(window: int = 200, threshold: float = 0.6,
                    organism: str = "") -> ConstraintSpec:
    params: dict = {"window": window, "threshold": threshold}
    if organism:
        params["organism"] = organism
    return ConstraintSpec(
        ctype=ConstraintType.NO_CPG,
        name="NoCpGIslandConstraint",
        params=params,
    )


def _spec_from_splice(threshold: float = 3.0,
                       organism: str = "") -> ConstraintSpec:
    params: dict = {"threshold": threshold}
    if organism:
        params["organism"] = organism
    return ConstraintSpec(
        ctype=ConstraintType.NO_CRYPTIC_SPLICE,
        name="NoCrypticSpliceConstraint",
        params=params,
    )


def _spec_from_restriction(sites: list[str]) -> ConstraintSpec:
    return ConstraintSpec(
        ctype=ConstraintType.RESTRICTION_SITE,
        name="NoRestrictionSiteConstraint",
        params={"sites": sites},
    )


def _spec_from_instability(motif: str = "ATTTA") -> ConstraintSpec:
    return ConstraintSpec(
        ctype=ConstraintType.NO_INSTABILITY_MOTIF,
        name="NoATTTAMotifConstraint",
        params={"motif": motif},
    )


def _spec_from_trun(max_run: int = 5) -> ConstraintSpec:
    return ConstraintSpec(
        ctype=ConstraintType.MRNA_STABILITY,
        name="NoTRunConstraint",
        params={"max_run": max_run},
    )


def _spec_from_translation(protein: str) -> ConstraintSpec:
    return ConstraintSpec(
        ctype=ConstraintType.AMINO_ACID_IDENTITY,
        name="TranslationConstraint",
        params={"protein": protein},
    )


# ═══════════════════════════════════════════════════════════════════════
# GC Content — ConstraintSpec vs GCRangeConstraint
# ═══════════════════════════════════════════════════════════════════════

class TestGCContentAlignment:
    """ConstraintSpec.check() must agree with GCRangeConstraint.check()."""

    @pytest.mark.parametrize("seq,gc_lo,gc_hi", [
        ("ATGGTTTCTAAAGGTGAA", 0.20, 0.50),   # GC≈0.28 — inside
        ("GCGCGCGCATATATAT", 0.40, 0.60),      # GC=0.50 — inside
        ("ATATATATATATATAT", 0.30, 0.70),       # GC=0.00 — outside (too low)
        ("GCGCGCGCGCGCGCGC", 0.30, 0.70),       # GC=1.00 — outside (too high)
        ("ATGCATGC", 0.40, 0.60),               # GC=0.50 — inside
        ("ATGCATGC", 0.60, 0.80),               # GC=0.50 — outside (too low)
    ])
    def test_gc_spec_matches_constraint(self, seq: str, gc_lo: float, gc_hi: float):
        hc = GCRangeConstraint(gc_lo=gc_lo, gc_hi=gc_hi)
        spec = _spec_from_gc(gc_lo=gc_lo, gc_hi=gc_hi)
        assert spec.check(seq) == hc.check(seq), (
            f"Discrepancy: spec.check({seq!r})={spec.check(seq)}, "
            f"hc.check({seq!r})={hc.check(seq)} for gc=[{gc_lo},{gc_hi}]"
        )

    def test_gc_empty_sequence(self):
        hc = GCRangeConstraint(gc_lo=0.30, gc_hi=0.70)
        spec = _spec_from_gc()
        assert spec.check("") == hc.check("")

    def test_gc_boundary_values(self):
        """GC exactly at lower boundary should be satisfied (inclusive)."""
        # 4 G/C out of 10 = 0.40
        seq = "ATGCATGCAT"
        hc = GCRangeConstraint(gc_lo=0.40, gc_hi=0.50)
        spec = _spec_from_gc(gc_lo=0.40, gc_hi=0.50)
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is True


# ═══════════════════════════════════════════════════════════════════════
# CpG Island — ConstraintSpec vs NoCpGIslandConstraint
# ═══════════════════════════════════════════════════════════════════════

class TestCpGIslandAlignment:
    """ConstraintSpec.check() must agree with NoCpGIslandConstraint.check()."""

    @pytest.mark.parametrize("window,threshold", [
        (200, 0.6),
        (100, 0.6),
        (200, 0.8),
        (50, 0.65),
    ])
    def test_cpg_no_island_sequence(self, window: int, threshold: float):
        """Sequence with no CpG island — both should return True."""
        # AT-rich sequence, no CG dinucleotides
        seq = "ATATATATAT" * 30  # 300 bp
        hc = NoCpGIslandConstraint(window=window, threshold=threshold)
        spec = _spec_from_cpg(window=window, threshold=threshold)
        assert spec.check(seq) == hc.check(seq)

    def test_cpg_with_island(self):
        """Sequence with CpG island — both should detect violation."""
        # Construct a window-sized sequence rich in CG dinucleotides
        # CGCGCG... repeated gives c_count ~ window/2, g_count ~ window/2
        # cg_count ~ window/2, expected ~ window/4, obs_exp ~ 2.0
        window = 50
        threshold = 0.6
        seq = "CGCGCGCGCG" * 5  # 50 bp, all CG
        hc = NoCpGIslandConstraint(window=window, threshold=threshold)
        spec = _spec_from_cpg(window=window, threshold=threshold)
        assert spec.check(seq) == hc.check(seq)
        # This should be a violation (obs/exp >> 0.6)
        assert spec.check(seq) is False

    def test_cpg_short_sequence(self):
        """Sequence shorter than window — both should return True."""
        seq = "ATGCGT"
        hc = NoCpGIslandConstraint(window=200, threshold=0.6)
        spec = _spec_from_cpg()
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is True

    def test_cpg_moderate_sequence(self):
        """Mixed GC content — both should agree."""
        seq = (
            "ATGGTTTCTAAAGGTGAA" * 5    # 90 bp, some CG
            + "CGCGATCGATCGATCGATCG" * 5   # 100 bp, more CG
        )
        hc = NoCpGIslandConstraint(window=100, threshold=0.6)
        spec = _spec_from_cpg(window=100, threshold=0.6)
        assert spec.check(seq) == hc.check(seq)


# ═══════════════════════════════════════════════════════════════════════
# Cryptic Splice — ConstraintSpec vs NoCrypticSpliceConstraint
# ═══════════════════════════════════════════════════════════════════════

class TestCrypticSpliceAlignment:
    """ConstraintSpec.check() must agree with NoCrypticSpliceConstraint.check()."""

    def test_splice_no_gt_ag(self):
        """Sequence without GT or AG dinucleotides — both should return True."""
        seq = "ATCCCTTTAAACCC"
        hc = NoCrypticSpliceConstraint(threshold=3.0)
        spec = _spec_from_splice(threshold=3.0)
        assert spec.check(seq) == hc.check(seq)

    def test_splice_with_dinucleotides(self):
        """Sequence with GT/AG but below threshold — both should agree."""
        # A typical coding sequence will have GT/AG but most won't be
        # strong splice sites
        seq = "ATGGTTTCTAAAGGTGAAATGCATGCTAGCTAG"
        hc = NoCrypticSpliceConstraint(threshold=3.0)
        spec = _spec_from_splice(threshold=3.0)
        assert spec.check(seq) == hc.check(seq)

    def test_splice_very_high_threshold(self):
        """With an extremely high threshold, everything passes."""
        seq = "ATGGTTTCTAAAGGTGAA"
        hc = NoCrypticSpliceConstraint(threshold=100.0)
        spec = _spec_from_splice(threshold=100.0)
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is True

    def test_splice_very_low_threshold(self):
        """With a very low threshold, GT/AG sites may fail."""
        seq = "ATGGTTTCTAAAGGTGAA"
        hc = NoCrypticSpliceConstraint(threshold=0.01)
        spec = _spec_from_splice(threshold=0.01)
        assert spec.check(seq) == hc.check(seq)

    def test_splice_sequence_passed_is_uppercased(self):
        """Verify that score_donor receives uppercased sequence (the F1.7 fix).

        If the original (non-uppercased) sequence were passed to
        score_donor, results could differ.  We test with a mixed-case
        sequence that has GT dinucleotides in lowercase.
        """
        # This sequence has 'gt' in lowercase — if we don't upper before
        # scanning for GT, it would be missed.  But the spec check should
        # match the HardConstraint which also uppercases.
        seq_upper = "ATGGTTTCTAAAGGTGAA"
        seq_mixed = "ATGgttTCTAAAGGTGAA"
        hc = NoCrypticSpliceConstraint(threshold=3.0)
        spec = _spec_from_splice(threshold=3.0)
        # HardConstraint uppercases internally
        hc_result = hc.check(seq_mixed)
        # ConstraintSpec should also handle mixed case correctly
        spec_result = spec.check(seq_mixed)
        assert spec_result == hc_result, (
            f"Mixed-case handling mismatch: spec={spec_result}, hc={hc_result}"
        )


# ═══════════════════════════════════════════════════════════════════════
# Restriction Site — ConstraintSpec vs NoRestrictionSiteConstraint
# ═══════════════════════════════════════════════════════════════════════

class TestRestrictionSiteAlignment:
    """ConstraintSpec.check() must agree with NoRestrictionSiteConstraint.check()."""

    def test_restriction_no_site(self):
        """Sequence without restriction sites — both should return True."""
        seq = "ATGGTTTCTAAAGGTGAA"
        sites = ["GAATTC"]  # EcoRI
        hc = NoRestrictionSiteConstraint(sites=sites)
        spec = _spec_from_restriction(sites=sites)
        assert spec.check(seq) == hc.check(seq)

    def test_restriction_site_present(self):
        """Sequence with EcoRI site — both should return False."""
        seq = "ATGGAATTCGTTTCTAAA"
        sites = ["GAATTC"]
        hc = NoRestrictionSiteConstraint(sites=sites)
        spec = _spec_from_restriction(sites=sites)
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is False

    def test_restriction_rc_site(self):
        """Sequence with reverse complement of restriction site."""
        seq = "ATGGAATTCTCTAAAGGT"
        sites = ["GAATTC"]  # RC = GAATTC (palindrome)
        hc = NoRestrictionSiteConstraint(sites=sites)
        spec = _spec_from_restriction(sites=sites)
        assert spec.check(seq) == hc.check(seq)

    def test_restriction_multiple_sites(self):
        """Multiple restriction sites to avoid."""
        seq = "ATGAAGCTTGATCC"  # Has AAGCTT (HindIII) inside
        sites = ["GAATTC", "GGATCC", "AAGCTT"]
        hc = NoRestrictionSiteConstraint(sites=sites)
        spec = _spec_from_restriction(sites=sites)
        assert spec.check(seq) == hc.check(seq)

    def test_restriction_empty_sites(self):
        """Empty sites list — both should return True."""
        seq = "ATGGTTTCTAAAGGTGAA"
        hc = NoRestrictionSiteConstraint(sites=[])
        spec = _spec_from_restriction(sites=[])
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is True


# ═══════════════════════════════════════════════════════════════════════
# Instability Motif — ConstraintSpec vs NoATTTAMotifConstraint
# ═══════════════════════════════════════════════════════════════════════

class TestInstabilityMotifAlignment:
    """ConstraintSpec.check() must agree with NoATTTAMotifConstraint.check()."""

    def test_motif_absent(self):
        """Sequence without ATTTA — both should return True."""
        seq = "ATGGTTTCTAAAGGTGAA"
        hc = NoATTTAMotifConstraint()
        spec = _spec_from_instability()
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is True

    def test_motif_present(self):
        """Sequence with ATTTA — both should return False."""
        seq = "ATGATTTAGTTTCTAAA"
        hc = NoATTTAMotifConstraint()
        spec = _spec_from_instability()
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is False

    def test_motif_at_end(self):
        """ATTTA at the end of the sequence."""
        seq = "ATGGTTTCTAAAATTTA"
        hc = NoATTTAMotifConstraint()
        spec = _spec_from_instability()
        assert spec.check(seq) == hc.check(seq)

    def test_motif_custom(self):
        """Custom motif parameter."""
        seq = "ATGTTATTTATGAA"  # Contains TTATTTAT
        hc = NoATTTAMotifConstraint()  # Checks for ATTTA
        spec = _spec_from_instability(motif="ATTTA")
        assert spec.check(seq) == hc.check(seq)


# ═══════════════════════════════════════════════════════════════════════
# T-Run — ConstraintSpec vs NoTRunConstraint
# ═══════════════════════════════════════════════════════════════════════

class TestTRunAlignment:
    """ConstraintSpec.check() must agree with NoTRunConstraint.check()."""

    def test_trun_short(self):
        """No T-run exceeds max — both should return True."""
        seq = "ATGGTTTCTAAAGGTGAA"
        hc = NoTRunConstraint(max_run=5)
        spec = _spec_from_trun(max_run=5)
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is True

    def test_trun_long(self):
        """T-run of 6 exceeds max_run=5 — both should return False."""
        seq = "ATGTTTTTTGAA"  # 6 consecutive T's
        hc = NoTRunConstraint(max_run=5)
        spec = _spec_from_trun(max_run=5)
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is False

    def test_trun_exactly_at_limit(self):
        """T-run exactly at max_run — should still be satisfied."""
        seq = "ATGTTTTTGAA"  # 5 consecutive T's
        hc = NoTRunConstraint(max_run=5)
        spec = _spec_from_trun(max_run=5)
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is True

    def test_trun_custom_max_run(self):
        """Custom max_run value."""
        seq = "ATGTTTTGAA"  # 4 T's
        hc = NoTRunConstraint(max_run=3)
        spec = _spec_from_trun(max_run=3)
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is False


# ═══════════════════════════════════════════════════════════════════════
# Translation — ConstraintSpec vs TranslationConstraint
# ═══════════════════════════════════════════════════════════════════════

class TestTranslationAlignment:
    """ConstraintSpec.check() must agree with TranslationConstraint.check()."""

    def test_translation_correct(self):
        """Correctly translated sequence — both should return True."""
        protein = "MVSKGE"
        seq = "ATGGTTTCTAAAGGTGAA"  # Correct translation
        hc = TranslationConstraint(protein)
        spec = _spec_from_translation(protein)
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is True

    def test_translation_wrong_codon(self):
        """Wrong codon at position — both should return False."""
        protein = "MVSKGE"
        seq = "ATGAAATCTAAAGGTGAA"  # AAA at pos 1 → K, not V
        hc = TranslationConstraint(protein)
        spec = _spec_from_translation(protein)
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is False

    def test_translation_wrong_length(self):
        """Wrong length sequence — both should return False."""
        protein = "MVSKGE"
        seq = "ATGGTTTCTAAA"  # Too short
        hc = TranslationConstraint(protein)
        spec = _spec_from_translation(protein)
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is False

    def test_translation_empty_protein_empty_seq(self):
        """Empty protein with empty sequence — both should return True."""
        # Note: With empty protein, TranslationConstraint checks
        # len(seq) == 0 which is vacuously true for empty seq.
        seq = ""
        hc = TranslationConstraint("")
        spec = _spec_from_translation("")
        assert spec.check(seq) == hc.check(seq)
        assert spec.check(seq) is True


# ═══════════════════════════════════════════════════════════════════════
# Organism-awareness for NO_CRYPTIC_SPLICE and NO_CPG
# ═══════════════════════════════════════════════════════════════════════

class TestOrganismAwareness:
    """Prokaryotic organisms should skip eukaryote-specific constraints."""

    def test_cryptic_splice_prokaryote_returns_true(self):
        """For prokaryotic organisms, NO_CRYPTIC_SPLICE should always be
        satisfied (True) regardless of sequence content."""
        # Sequence with GT/AG that might be splice sites
        seq = "ATGGTTTCTAAAGGTGAA"
        spec = _spec_from_splice(threshold=0.01, organism="E_coli_K12")
        assert spec.check(seq) is True

    def test_cryptic_splice_eukaryote_checks_normally(self):
        """For eukaryotic organisms, NO_CRYPTIC_SPLICE should still check."""
        seq = "ATGGTTTCTAAAGGTGAA"
        spec = _spec_from_splice(threshold=0.01, organism="Homo_sapiens")
        # At threshold=0.01, GT/AG sites may be flagged
        hc = NoCrypticSpliceConstraint(threshold=0.01)
        assert spec.check(seq) == hc.check(seq)

    def test_cpg_prokaryote_returns_true(self):
        """For prokaryotic organisms, NO_CPG should always be satisfied."""
        # CG-rich sequence
        seq = "CGCGCGCGCG" * 5
        spec = _spec_from_cpg(window=50, threshold=0.6, organism="E_coli_K12")
        assert spec.check(seq) is True

    def test_cpg_eukaryote_checks_normally(self):
        """For eukaryotic organisms, NO_CPG should still check."""
        seq = "CGCGCGCGCG" * 5
        spec = _spec_from_cpg(window=50, threshold=0.6, organism="Homo_sapiens")
        hc = NoCpGIslandConstraint(window=50, threshold=0.6)
        assert spec.check(seq) == hc.check(seq)

    def test_no_organism_param_checks_normally(self):
        """Without organism param, constraint should check normally."""
        seq = "ATGGTTTCTAAAGGTGAA"
        spec_no_org = _spec_from_splice(threshold=0.01)
        hc = NoCrypticSpliceConstraint(threshold=0.01)
        assert spec_no_org.check(seq) == hc.check(seq)

    def test_prokaryote_aliases(self):
        """Legacy aliases for prokaryotic organisms should also work."""
        seq = "ATGGTTTCTAAAGGTGAA"
        for organism in ["E_coli_K12", "E_coli_BL21", "ecoli", "Escherichia_coli"]:
            spec = _spec_from_splice(threshold=0.01, organism=organism)
            assert spec.check(seq) is True, (
                f"Expected True for prokaryotic organism {organism!r}"
            )

    def test_eukaryote_aliases(self):
        """Legacy aliases for eukaryotic organisms should still check."""
        seq = "ATGGTTTCTAAAGGTGAA"
        for organism in ["Homo_sapiens", "human", "Mus_musculus", "mouse"]:
            spec = _spec_from_splice(threshold=0.01, organism=organism)
            hc = NoCrypticSpliceConstraint(threshold=0.01)
            assert spec.check(seq) == hc.check(seq), (
                f"Discrepancy for eukaryotic organism {organism!r}"
            )

    def test_cpg_prokaryote_via_ecoli_bl21(self):
        """E. coli BL21 (prokaryote) should skip CpG check."""
        seq = "CGCGCGCGCG" * 5
        spec = _spec_from_cpg(window=50, threshold=0.6, organism="E_coli_BL21")
        assert spec.check(seq) is True


# ═══════════════════════════════════════════════════════════════════════
# Soft / always-satisfied constraints
# ═══════════════════════════════════════════════════════════════════════

class TestSoftConstraints:
    """Soft constraints (CODON_USAGE, CODON_PAIR_BIAS) should always be
    satisfied via ConstraintSpec.check()."""

    def test_codon_usage_always_satisfied(self):
        spec = ConstraintSpec(
            ctype=ConstraintType.CODON_USAGE,
            name="MaximizeCAI",
        )
        assert spec.check("ATGGTTTCTAAAGGTGAA") is True

    def test_codon_pair_bias_always_satisfied(self):
        spec = ConstraintSpec(
            ctype=ConstraintType.CODON_PAIR_BIAS,
            name="MinimizeCodonPairBias",
        )
        assert spec.check("ATGGTTTCTAAAGGTGAA") is True

    def test_mhc_binding_always_satisfied(self):
        spec = ConstraintSpec(
            ctype=ConstraintType.MHC_BINDING,
            name="MHCBinding",
        )
        assert spec.check("ATGGTTTCTAAAGGTGAA") is True

    def test_tcell_epitope_always_satisfied(self):
        spec = ConstraintSpec(
            ctype=ConstraintType.TCELL_EPITOPE,
            name="TCellEpitope",
        )
        assert spec.check("ATGGTTTCTAAAGGTGAA") is True

    def test_protein_stability_always_satisfied(self):
        spec = ConstraintSpec(
            ctype=ConstraintType.PROTEIN_STABILITY,
            name="ProteinStability",
        )
        assert spec.check("ATGGTTTCTAAAGGTGAA") is True


# ═══════════════════════════════════════════════════════════════════════
# Edge cases
# ═══════════════════════════════════════════════════════════════════════

class TestEdgeCases:
    """Edge cases for ConstraintSpec.check()."""

    def test_empty_sequence(self):
        """Empty sequence should always be satisfied."""
        for ctype in [ConstraintType.GC_CONTENT, ConstraintType.NO_CPG,
                       ConstraintType.NO_CRYPTIC_SPLICE, ConstraintType.RESTRICTION_SITE,
                       ConstraintType.NO_INSTABILITY_MOTIF, ConstraintType.MRNA_STABILITY]:
            spec = ConstraintSpec(ctype=ctype, name=f"test_{ctype.value}")
            assert spec.check("") is True

    def test_unknown_constraint_type(self):
        """Unknown constraint types should default to satisfied."""
        spec = ConstraintSpec(ctype=ConstraintType.CUSTOM, name="UnknownConstraint")
        assert spec.check("ATGGTTTCTAAAGGTGAA") is True

    def test_no_gt_dinucleotide(self):
        """NO_GT_DINUCLEOTIDE type — no GT in sequence = satisfied."""
        seq_no_gt = "ATCCCATCCC"
        seq_with_gt = "ATGGTTCCC"
        spec = ConstraintSpec(
            ctype=ConstraintType.NO_GT_DINUCLEOTIDE,
            name="NoGTDinucleotide",
        )
        assert spec.check(seq_no_gt) is True
        assert spec.check(seq_with_gt) is False
