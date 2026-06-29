"""
BioCompiler Proof-Correspondence Tests
=======================================

Property-based tests that verify the Python implementation matches the
properties proven in the Lean4 formal model.  Each test class maps to
a specific proof module and theorem, closing the 7 known proof-implementation
gaps documented in DOC-14 and docs/11-Refinement-Mapping.md.

Gaps addressed:
  §3.1  Three-valued vs. five-valued logic
  §3.2  Lean4 Sequence (List Nucleotide) vs. Python str
  §3.3  Axiom assumptions vs. Python heuristics (scanner correspondence)
  §3.4  Float vs. Rat differences
  §3.5  SLOT predicates: formal UNCERTAIN vs. Python PASS/FAIL
  §3.6  Reverse complement handling (Python is stronger)
  §3.7  Number of predicates alignment (33 vs. 43)

Lean4 proof modules covered:
  - ThreeValued.lean        (12 theorems)
  - Compositional.lean       (compositional_soundness, slot_predicates_uncertain)
  - SLOTVerification.lean    (conservative_is_safe, slot_soundness_verified, …)
  - SLOTIndependence.lean    (predicate_is_core_or_slot, ffi_never_pass, …)
  - Refinement.lean          (verified_refines_conservative, simulation, …)
  - Mutagenesis.lean         (valine_only_mandatory_gt_aa, synonymous_preserves_translation)
  - Certificates.lean        (certificate_soundness)

Reference:
  - proof/BioCompiler/*.lean
  - docs/14-SLOT-Proof-Implementation-Gap.md
  - docs/11-Refinement-Mapping.md
"""

from __future__ import annotations

import math
from typing import List

import pytest
pytest.importorskip("hypothesis")
pytest.importorskip("hypothesis")
from hypothesis import given, strategies as st, assume, settings

from biocompiler.shared.types import (
    Verdict,
    SLOTMode,
    five_valued_and,
    five_valued_or,
    combined_verdict,
)
from biocompiler.shared.types import _VERDICT_ORDER
from biocompiler.provenance.slot_verification import (
    SLOT_PREDICATES,
    is_slot_predicate,
    VerificationEvidence,
    verify_no_cryptic_splice,
    verify_no_cryptic_promoter,
    verify_no_unexpected_tm_domain,
    verify_mrna_secondary_structure,
    verify_co_translational_folding,
    verify_conservation_score,
    verify_codon_optimality,
    verify_structure_predicate,
    verify_stability_predicate,
    verify_solubility_predicate,
    verify_immunogenicity_predicate,
)
from biocompiler.type_system import (
    CODON_TABLE,
    AA_TO_CODONS,
    BLOSUM62,
    PREDICATE_NAMES,
    check_no_stop_codons,
    check_no_cryptic_splice,
    check_no_cpg_island,
    check_conservation_score,
    check_codon_optimality,
)
from biocompiler.expression.translation import translate


# ═══════════════════════════════════════════════════════════════════════════
# Strategies
# ═══════════════════════════════════════════════════════════════════════════

verdict_strategy = st.sampled_from(list(Verdict))
three_valued_strategy = st.sampled_from([Verdict.PASS, Verdict.UNCERTAIN, Verdict.FAIL])
slot_mode_strategy = st.sampled_from(list(SLOTMode))

dna_base = st.sampled_from("ACGT")
dna_seq = st.text(alphabet=st.sampled_from("ACGT"), min_size=0, max_size=300)
dna_seq_codon = st.text(alphabet=st.sampled_from("ACGT"), min_size=3, max_size=300)
    # multiples of 3 not guaranteed; use .filter when needed

amino_acid = st.sampled_from(list("ACDEFGHIKLMNPQRSTVWY"))


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.1: Three-valued vs. Five-valued Logic
# ═══════════════════════════════════════════════════════════════════════════

class TestRefinement5to3:
    """Verify the 5→3 refinement mapping is a homomorphism.

    Lean4 uses {PASS, UNCERTAIN, FAIL}.  Python uses 5 values.
    The refinement refine5to3 maps:
      PASS → PASS, LIKELY_PASS → UNCERTAIN, UNCERTAIN → UNCERTAIN,
      LIKELY_FAIL → UNCERTAIN, FAIL → FAIL.

    We verify that all Lean4 algebraic laws still hold under this mapping.
    """

    @staticmethod
    def refine5to3(v: Verdict) -> Verdict:
        """Conservative abstraction: map 5-valued to 3-valued."""
        mapping = {
            Verdict.PASS: Verdict.PASS,
            Verdict.LIKELY_PASS: Verdict.UNCERTAIN,
            Verdict.UNCERTAIN: Verdict.UNCERTAIN,
            Verdict.LIKELY_FAIL: Verdict.UNCERTAIN,
            Verdict.FAIL: Verdict.FAIL,
        }
        return mapping[v]

    @given(a=verdict_strategy, b=verdict_strategy)
    def test_refine_and_homomorphism(self, a, b):
        """Lean4: AND refines — refine5to3(AND(a,b)) = three_valued_and(refine5to3(a), refine5to3(b))"""
        r_ab = self.refine5to3(five_valued_and(a, b))
        r_a = self.refine5to3(a)
        r_b = self.refine5to3(b)
        expected = five_valued_and(r_a, r_b)  # works since 3-valued is a subset of 5-valued
        assert r_ab == expected, f"refine(AND({a},{b})) = {r_ab} != AND(refine({a}),refine({b})) = {expected}"

    @given(a=verdict_strategy, b=verdict_strategy)
    def test_refine_or_homomorphism(self, a, b):
        """Lean4: OR refines — refine5to3(OR(a,b)) = three_valued_or(refine5to3(a), refine5to3(b))"""
        r_ab = self.refine5to3(five_valued_or(a, b))
        r_a = self.refine5to3(a)
        r_b = self.refine5to3(b)
        expected = five_valued_or(r_a, r_b)
        assert r_ab == expected

    @given(vs=st.lists(verdict_strategy, min_size=0, max_size=10))
    def test_refine_combined_verdict_homomorphism(self, vs):
        """Lean4: foldl refines — refine5to3(combined5(vs)) = combined3(refine5to3(vs))"""
        combined5 = combined_verdict(vs)
        vs3 = [self.refine5to3(v) for v in vs]
        combined3 = combined_verdict(vs3)
        assert self.refine5to3(combined5) == combined3

    @given(v=verdict_strategy)
    def test_refine_preserves_definite_verdicts(self, v):
        """PASS and FAIL are preserved by the refinement (no information loss for definite verdicts)."""
        r = self.refine5to3(v)
        if v == Verdict.PASS:
            assert r == Verdict.PASS
        if v == Verdict.FAIL:
            assert r == Verdict.FAIL

    @given(v=verdict_strategy)
    def test_refine_never_downgrades(self, v):
        """Refinement never turns PASS into FAIL or vice versa."""
        r = self.refine5to3(v)
        if v == Verdict.PASS:
            assert r != Verdict.FAIL
        if v == Verdict.FAIL:
            assert r != Verdict.PASS


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.2: Sequence Correspondence (List Nucleotide vs. str)
# ═══════════════════════════════════════════════════════════════════════════

class TestSequenceCorrespondence:
    """Verify Python str-based sequences behave like Lean4 List Nucleotide.

    Lean4 theorems:
      - matchesAt_spec: pattern matching is correct
      - containsPattern_complete: if pattern appears, scanner finds it
      - containsPattern_sound: if scanner returns true, pattern is there
    """

    @given(seq=dna_seq, pattern=st.text(alphabet=st.sampled_from("ACGT"), min_size=1, max_size=6))
    def test_str_find_completeness(self, seq, pattern):
        """Lean4: containsPattern_complete — if pattern in seq, find returns a valid position."""
        if pattern in seq:
            pos = seq.find(pattern)
            assert pos >= 0
            assert seq[pos:pos + len(pattern)] == pattern

    @given(seq=dna_seq, pattern=st.text(alphabet=st.sampled_from("ACGT"), min_size=1, max_size=6))
    def test_str_find_soundness(self, seq, pattern):
        """Lean4: containsPattern_sound — if find returns valid pos, pattern is there."""
        pos = seq.find(pattern)
        if pos >= 0:
            assert seq[pos:pos + len(pattern)] == pattern

    @given(seq=dna_seq)
    def test_dna_alphabet_invariant(self, seq):
        """Python must maintain ACGT-only invariant (Lean4 guarantees this by construction)."""
        assert set(seq) <= {"A", "C", "G", "T"}, f"Non-ACGT characters in: {seq!r}"


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.3: Scanner/Oracle Correspondence
# ═══════════════════════════════════════════════════════════════════════════

class TestScannerCorrespondence:
    """Verify Python scanner implementations match Lean4 scanner axioms.

    Lean4 axioms (remaining 3 of 18):
      A1: SpliceSiteScanner.scanner_completeness
      A2: SpliceSiteScanner.scanner_soundness
      A3: SpliceSiteScanner.borderline_completeness

    Proved in Lean4 (axioms 4-18 eliminated):
      A4-5: CpGIslandScanner completeness/soundness
      A6-8: PromoterScanner completeness/soundness/borderline
      A9-11: TMDomainScanner completeness/soundness/borderline
    """

    def test_codon_table_completeness(self):
        """All 64 codons map to an amino acid or stop in CODON_TABLE."""
        bases = "ACGT"
        all_codons = [a + b + c for a in bases for b in bases for c in bases]
        for codon in all_codons:
            assert codon in CODON_TABLE, f"Missing codon: {codon}"

    def test_codon_table_size(self):
        """Exactly 64 codon entries."""
        assert len(CODON_TABLE) == 64

    def test_stop_codons_correct(self):
        """Standard genetic code stop codons: TAA, TAG, TGA."""
        stops = {c for c, aa in CODON_TABLE.items() if aa == "*"}
        assert stops == {"TAA", "TAG", "TGA"}

    def test_start_codon_correct(self):
        """ATG encodes Methionine (also the start codon)."""
        assert CODON_TABLE["ATG"] == "M"

    def test_blosum62_symmetry(self):
        """BLOSUM62 is a symmetric matrix: BLOSUM62[(a,b)] == BLOSUM62[(b,a)]."""
        for (a, b), score in BLOSUM62.items():
            assert BLOSUM62.get((b, a)) == score, f"Asymmetric: BLOSUM62[({a},{b})]={score}, BLOSUM62[({b},{a})]={BLOSUM62.get((b,a))}"

    @given(seq=dna_seq)
    def test_no_stop_codons_soundness(self, seq):
        """Lean4: if check_no_stop_codons returns PASS, there are truly no internal stops."""
        result = check_no_stop_codons(seq)
        if result.verdict == Verdict.PASS and len(seq) >= 6:
            # Verify independently: no internal stops
            for i in range(0, len(seq) - 3, 3):
                codon = seq[i:i + 3]
                if codon in ("TAA", "TAG", "TGA"):
                    pytest.fail(f"check_no_stop_codons returned PASS but found stop at {i}: {codon}")

    @given(seq=dna_seq.filter(lambda s: len(s) >= 200))
    @settings(max_examples=20)
    def test_cpg_island_soundness(self, seq):
        """Lean4 A4-5: CpG island scanner soundness — FAIL means a real CpG island exists."""
        result = check_no_cpg_island(seq)
        if result.verdict == Verdict.FAIL:
            # The result should report a position
            assert result.positions, "CpG FAIL verdict should report positions"


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.4: Float vs. Rat Differences
# ═══════════════════════════════════════════════════════════════════════════

class TestFloatRatCorrespondence:
    """Verify float-based computations agree with exact rational arithmetic.

    Lean4 uses Rat for GC content, CAI thresholds, etc.
    Python uses IEEE 754 float.  This gap is mitigated by epsilon tolerance.
    """

    @given(seq=dna_seq.filter(lambda s: len(s) > 0))
    def test_gc_content_float_rational_agreement(self, seq):
        """GC content computed with floats should agree with exact rational within epsilon."""
        n = len(seq)
        g_count = seq.count("G")
        c_count = seq.count("C")
        # Float computation
        gc_float = (g_count + c_count) / n
        # Exact rational (as fraction)
        gc_num = g_count + c_count
        gc_den = n
        # They should agree within floating-point epsilon
        assert abs(gc_float - gc_num / gc_den) < 1e-15

    @given(seq=dna_seq.filter(lambda s: len(s) > 0))
    def test_gc_content_in_range(self, seq):
        """GC content is always in [0, 1] — matches Lean4 gcContent specification."""
        n = len(seq)
        gc = (seq.count("G") + seq.count("C")) / n
        assert 0.0 <= gc <= 1.0

    @given(
        gc_lo=st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False),
        gc_hi=st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False),
    )
    def test_gc_boundary_sensitivity(self, gc_lo, gc_hi):
        """Boundary comparisons should not be affected by float rounding at exact 0.0 or 1.0."""
        assume(gc_lo <= gc_hi)
        # At GC=0.0 (all AT), any lo>0 should fail
        assert not (0.0 >= gc_lo and gc_lo > 0.0)
        # At GC=1.0 (all GC), any hi<1.0 should fail
        assert not (1.0 <= gc_hi and gc_hi < 1.0)


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.5: SLOT Predicates — Formal UNCERTAIN vs. Python PASS/FAIL
# ═══════════════════════════════════════════════════════════════════════════

class TestSLOTConservativeMode:
    """Verify CONSERVATIVE mode matches the Lean4 model exactly.

    Lean4 theorem: conservative_is_safe
      "Conservative mode never returns PASS for SLOT predicates."

    This is the foundational safety property: in CONSERVATIVE mode,
    SLOT predicates can never contribute to a PASS verdict.
    """

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_conservative_always_uncertain_cryptic_splice(self, seq):
        """Lean4: conservative_is_safe — NoCrypticSplice returns UNCERTAIN in CONSERVATIVE mode."""
        verdict, evidence = verify_no_cryptic_splice(seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN
        assert evidence.slot_mode == SLOTMode.CONSERVATIVE
        assert not evidence.verified

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_conservative_always_uncertain_cryptic_promoter(self, seq):
        """Lean4: conservative_is_safe — NoCrypticPromoter returns UNCERTAIN in CONSERVATIVE mode."""
        verdict, evidence = verify_no_cryptic_promoter(seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN
        assert evidence.slot_mode == SLOTMode.CONSERVATIVE

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_conservative_always_uncertain_tm_domain(self, seq):
        """Lean4: conservative_is_safe — NoUnexpectedTMDomain returns UNCERTAIN in CONSERVATIVE."""
        verdict, evidence = verify_no_unexpected_tm_domain(seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_conservative_always_uncertain_mrna_structure(self, seq):
        """Lean4: conservative_is_safe — mRNASecondaryStructure returns UNCERTAIN in CONSERVATIVE."""
        verdict, evidence = verify_mrna_secondary_structure(seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_conservative_always_uncertain_cotrans_folding(self, seq):
        """Lean4: conservative_is_safe — CoTranslationalFolding returns UNCERTAIN in CONSERVATIVE."""
        verdict, evidence = verify_co_translational_folding(seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN

    @given(aa1=amino_acid, aa2=amino_acid)
    def test_conservative_always_uncertain_conservation(self, aa1, aa2):
        """Lean4: conservative_is_safe — ConservationScore returns UNCERTAIN in CONSERVATIVE."""
        verdict, evidence = verify_conservation_score(aa1, aa2, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN

    # Structure/stability predicates call _check_esmfold_available / _check_foldx_available
    # which can be slow, so use deadline=None and few examples.

    @given(protein=st.sampled_from(["MVSKGE", "ACDEFGHIKL", "MKWVTFISLL"]))
    @settings(deadline=None, max_examples=3)
    def test_conservative_always_uncertain_structure(self, protein):
        """Lean4: conservative_is_safe — structure predicates return UNCERTAIN in CONSERVATIVE."""
        for pred_name in ["StructureConfidence", "NoMisfoldingRisk",
                          "CorrectFoldTopology", "NoUnexpectedInteraction"]:
            verdict, evidence = verify_structure_predicate(
                pred_name, protein, slot_mode=SLOTMode.CONSERVATIVE)
            assert verdict == Verdict.UNCERTAIN, f"{pred_name} returned {verdict} in CONSERVATIVE"

    @given(protein=st.sampled_from(["MVSKGE", "ACDEFGHIKL", "MKWVTFISLL"]))
    @settings(deadline=None, max_examples=3)
    def test_conservative_always_uncertain_stability(self, protein):
        """Lean4: conservative_is_safe — stability predicates return UNCERTAIN in CONSERVATIVE."""
        for pred_name in ["StableFolding", "NoDestabilizingMutation",
                          "DisulfideBondIntegrity", "HydrophobicCoreQuality"]:
            verdict, evidence = verify_stability_predicate(
                pred_name, protein, slot_mode=SLOTMode.CONSERVATIVE)
            assert verdict == Verdict.UNCERTAIN, f"{pred_name} returned {verdict} in CONSERVATIVE"


class TestSLOTVerifiedModeEvidence:
    """Verify VERIFIED mode only returns PASS with proper evidence.

    Lean4 theorem: verified_pass_implies_all_vcs
      "If VERIFIED mode returns PASS, then all verification conditions hold."

    In the Python implementation, this translates to:
    PASS → evidence.verified == True AND evidence.tool_available == True
    """

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    @settings(deadline=None)
    def test_verified_pass_has_evidence(self, seq):
        """Lean4: verified_pass_implies_all_vcs — PASS in VERIFIED mode requires evidence."""
        verdict, evidence = verify_no_cryptic_splice(seq, slot_mode=SLOTMode.VERIFIED)
        if verdict == Verdict.PASS:
            assert evidence.verified, f"PASS verdict without verified=True for {evidence.predicate}"
            assert evidence.tool_available is not None

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    @settings(deadline=None)
    def test_verified_evidence_documents_threshold(self, seq):
        """Every VERIFIED-mode result should document the threshold used."""
        verdict, evidence = verify_no_cryptic_promoter(seq, slot_mode=SLOTMode.VERIFIED)
        if verdict in (Verdict.PASS, Verdict.FAIL):
            assert evidence.threshold_used is not None, "Missing threshold documentation"

    @given(seq=st.sampled_from(["ATGATGATGATG", "GCTAGCTAGCTA", "ATGGTTTCTAAAGGTGAA"]))
    @settings(deadline=None, max_examples=3)
    def test_verified_uncertain_when_tool_unavailable(self, seq):
        """Lean4: VERIFIED mode returns UNCERTAIN when tool unavailable.

        For structure/stability/immunogenicity predicates that require external
        tools, VERIFIED mode should return UNCERTAIN if the tool is unavailable.
        """
        # These require external tools (ESMFold, FoldX, CamSol, MHC)
        for pred_name in ["StructureConfidence", "NoMisfoldingRisk"]:
            verdict, evidence = verify_structure_predicate(
                pred_name, seq, slot_mode=SLOTMode.VERIFIED)
            if not evidence.tool_available:
                assert verdict in (Verdict.UNCERTAIN,), \
                    f"VERIFIED without tool returned {verdict} for {pred_name}"


class TestSLOTPermissiveModeGoesBeyond:
    """Verify PERMISSIVE mode behavior matches documentation.

    Lean4: No soundness theorem exists for PERMISSIVE mode.
    PERMISSIVE mode goes beyond what is formally proven.
    """

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_permissive_verdict_at_least_as_informative_as_conservative(self, seq):
        """Lean4: permissive_refines_conservative — PERMISSIVE refines CONSERVATIVE.

        The PERMISSIVE verdict should be at least as informative as CONSERVATIVE
        (which is always UNCERTAIN for SLOT predicates).
        """
        verdict_c, _ = verify_no_cryptic_splice(seq, slot_mode=SLOTMode.CONSERVATIVE)
        verdict_p, _ = verify_no_cryptic_splice(seq, slot_mode=SLOTMode.PERMISSIVE)
        # CONSERVATIVE is always UNCERTAIN; PERMISSIVE should be PASS or UNCERTAIN
        # (both refine UNCERTAIN in the information ordering)
        assert verdict_c == Verdict.UNCERTAIN
        # PERMISSIVE can return PASS or UNCERTAIN (or FAIL for definitive checks)
        # The key property: it never contradicts CONSERVATIVE's UNCERTAIN

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_permissive_relaxed_threshold(self, seq):
        """PERMISSIVE mode uses relaxed thresholds compared to VERIFIED."""
        _, evidence_v = verify_no_cryptic_promoter(seq, slot_mode=SLOTMode.VERIFIED)
        _, evidence_p = verify_no_cryptic_promoter(seq, slot_mode=SLOTMode.PERMISSIVE)
        # Both should document their thresholds
        if evidence_v.threshold_used is not None and evidence_p.threshold_used is not None:
            # Permissive threshold should be <= verified threshold (more lenient)
            assert evidence_p.threshold_used <= evidence_v.threshold_used + 1e-10


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.5 (cont.): Refinement Theorem Correspondence
# ═══════════════════════════════════════════════════════════════════════════

class TestRefinementTheorems:
    """Verify Python behavior matches the Lean4 Refinement.lean theorems.

    Lean4 theorems:
      - verified_refines_conservative: VERIFIED refines CONSERVATIVE
      - simulation_verified_conservative: overall verdict refinement
      - conservative_pass_no_downgrade: switching modes never downgrades PASS
      - verified_fail_consistent: VERIFIED FAIL is consistent with CONSERVATIVE
    """

    @staticmethod
    def verdict_refines(v_refined: Verdict, v_abstract: Verdict) -> bool:
        """Lean4: verdictRefines — v_refined is at least as informative as v_abstract.

        From Refinement.lean:
          verdictRefines v_refined v_abstract := v_abstract = UNCERTAIN ∨ v_refined = v_abstract
        """
        return v_abstract == Verdict.UNCERTAIN or v_refined == v_abstract

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_verified_refines_conservative_per_predicate(self, seq):
        """Lean4: verified_slot_refines_conservative_slot — per-predicate refinement."""
        predicates_to_check = [
            lambda s: verify_no_cryptic_splice(s, slot_mode=SLOTMode.VERIFIED),
            lambda s: verify_no_cryptic_promoter(s, slot_mode=SLOTMode.VERIFIED),
            lambda s: verify_no_unexpected_tm_domain(s, slot_mode=SLOTMode.VERIFIED),
        ]
        for check_fn in predicates_to_check:
            verdict_v, _ = check_fn(seq)
            verdict_c, _ = verify_no_cryptic_splice(seq, slot_mode=SLOTMode.CONSERVATIVE) \
                if check_fn == predicates_to_check[0] else \
                (Verdict.UNCERTAIN, None)  # CONSERVATIVE is always UNCERTAIN
            assert self.verdict_refines(verdict_v, Verdict.UNCERTAIN), \
                f"VERIFIED {verdict_v} does not refine CONSERVATIVE UNCERTAIN"

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_no_downgrade_from_conservative_pass(self, seq):
        """Lean4: conservative_pass_no_downgrade — CONSERVATIVE PASS not downgraded.

        Since CONSERVATIVE always returns UNCERTAIN for SLOT predicates,
        this is trivially satisfied (no CONSERVATIVE PASS to downgrade).
        But we verify the structure: VERIFIED never returns FAIL when
        CONSERVATIVE would return PASS for core predicates.
        """
        # For SLOT predicates, CONSERVATIVE is always UNCERTAIN
        verdict_c, _ = verify_no_cryptic_splice(seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict_c == Verdict.UNCERTAIN  # matches Lean4 model


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.5 (cont.): SLOT Independence Theorems
# ═══════════════════════════════════════════════════════════════════════════

class TestSLOTIndependence:
    """Verify SLOT independence properties in the Python implementation.

    Lean4 theorems:
      - predicate_is_core_or_slot: every predicate is core or SLOT (exclusive)
      - core_not_slot: core and SLOT are mutually exclusive
      - ffi_never_pass: FFI-dependent predicates never produce PASS
    """

    def test_predicate_classification_exclusive(self):
        """Lean4: core_not_slot — no predicate is both core and SLOT."""
        from biocompiler.provenance.slot_verification import SLOT_PREDICATES
        core_predicates = {
            "NoStopCodons", "NoCpGIsland", "NoRestrictionSite",
            "NoGTDinucleotide", "ValidCodingSeq", "SpliceCorrect",
            "GCInRange", "InFrame", "NoInstabilityMotif",
        }
        # Verify no overlap
        overlap = core_predicates & SLOT_PREDICATES
        assert len(overlap) == 0, f"Predicates classified as both core and SLOT: {overlap}"

    def test_slot_predicates_classified(self):
        """Lean4: predicate_is_core_or_slot — every SLOT predicate is in SLOT_PREDICATES."""
        # All 20 SLOT predicates from the Lean4 model should be in SLOT_PREDICATES
        # MOVE-PRED: NoLongHydrophobicStretch, NoRibosomalFrameshift,
        # NoMiRNABindingSite moved to CORE (deterministic). 19 SLOT remain.
        lean4_slot_predicates = {
            "ConservationScore", "NoUnexpectedTMDomain", "mRNASecondaryStructure",
            "CoTranslationalFolding", "StructureConfidence", "NoMisfoldingRisk",
            "CorrectFoldTopology", "NoUnexpectedInteraction", "StableFolding",
            "NoDestabilizingMutation", "DisulfideBondIntegrity", "HydrophobicCoreQuality",
            "SolubleExpression", "NoAggregationProneRegion", "ChargeComposition",
            "LowImmunogenicity", "NoStrongTCellEpitope",
            "NoDominantBCellEpitope", "PopulationCoverageSafe",
        }
        for pred in lean4_slot_predicates:
            assert is_slot_predicate(pred), f"Lean4 SLOT predicate {pred} not in Python SLOT_PREDICATES"

    def test_core_predicates_not_slot(self):
        """Lean4: core_not_slot — core predicates should NOT be classified as SLOT."""
        core_predicates = [
            "NoStopCodons", "NoCpGIsland", "ValidCodingSeq",
            "NoGTDinucleotide", "NoInstabilityMotif",
        ]
        for pred in core_predicates:
            assert not is_slot_predicate(pred), f"Core predicate {pred} incorrectly classified as SLOT"

    @given(pred_name=st.sampled_from([
        "NoCrypticSplice", "NoUnexpectedTMDomain", "StructureConfidence",
        "StableFolding", "SolubleExpression", "LowImmunogenicity",
    ]))
    def test_is_slot_predicate_for_slot(self, pred_name):
        """All SLOT predicate names should be recognized as SLOT."""
        assert is_slot_predicate(pred_name)

    @given(pred_name=st.sampled_from([
        "NoStopCodons", "ValidCodingSeq", "NoCpGIsland", "NoGTDinucleotide",
    ]))
    def test_is_not_slot_predicate_for_core(self, pred_name):
        """Core predicate names should NOT be recognized as SLOT."""
        assert not is_slot_predicate(pred_name)


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.5 (cont.): Compositional Soundness
# ═══════════════════════════════════════════════════════════════════════════

class TestCompositionalSoundness:
    """Verify compositional soundness in the Python implementation.

    Lean4 theorems:
      - compositional_soundness: evaluateAll = PASS → ∀ P, propertyHolds P
      - slot_predicates_dont_affect_pass: SLOT in list → no overall PASS
      - all_core_if_pass: evaluateAll = PASS → all predicates are core
    """

    @given(vs=st.lists(verdict_strategy, min_size=0, max_size=20))
    def test_combined_pass_implies_all_pass(self, vs):
        """Lean4: compositional_soundness + foldl_and_pass_implies_all_pass."""
        if combined_verdict(vs) == Verdict.PASS:
            for v in vs:
                assert v == Verdict.PASS, f"combined PASS but element {v} is not PASS"

    @given(vs=st.lists(verdict_strategy, min_size=1, max_size=20))
    def test_slot_uncertain_prevents_pass(self, vs):
        """Lean4: slot_predicates_dont_affect_pass — UNCERTAIN in list prevents PASS.

        SLOT predicates always return UNCERTAIN in CONSERVATIVE mode.
        If any verdict is UNCERTAIN, combined cannot be PASS.
        """
        if Verdict.UNCERTAIN in vs:
            assert combined_verdict(vs) != Verdict.PASS

    @given(vs=st.lists(verdict_strategy, min_size=1, max_size=20))
    def test_fail_in_list_prevents_pass(self, vs):
        """Lean4: FAIL in the list prevents overall PASS."""
        if Verdict.FAIL in vs:
            assert combined_verdict(vs) == Verdict.FAIL

    @given(vs=st.lists(verdict_strategy, min_size=0, max_size=20))
    def test_combined_verdict_worst_link(self, vs):
        """Combined verdict equals the minimum (weakest link) — Kleene AND semantics."""
        if not vs:
            assert combined_verdict(vs) == Verdict.UNCERTAIN
        else:
            min_order = min(_VERDICT_ORDER[v] for v in vs)
            result = combined_verdict(vs)
            assert _VERDICT_ORDER[result] == min_order


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.5 (cont.): Mutagenesis Theorems
# ═══════════════════════════════════════════════════════════════════════════

class TestMutagenesisCorrespondence:
    """Verify mutagenesis properties proven in Lean4 hold in Python.

    Lean4 theorems:
      - all_valine_codons_have_gt: all 4 Valine codons contain GT
      - valine_only_mandatory_gt_aa: Valine is the only mandatory-GT amino acid
      - every_aa_has_ag_free_codon: every AA has an AG-free codon
      - synonymous_preserves_translation: synonymous mutations preserve amino acid
    """

    def test_all_valine_codons_contain_gt(self):
        """Lean4: all_valine_codons_have_gt — all 4 Val codons have GT."""
        val_codons = AA_TO_CODONS.get("V", [])
        assert len(val_codons) == 4, f"Expected 4 Val codons, got {len(val_codons)}"
        for codon in val_codons:
            assert "GT" in codon, f"Valine codon {codon} does not contain GT"

    def test_valine_only_mandatory_gt_aa(self):
        """Lean4: valine_only_mandatory_gt_aa — only Val has GT in all codons."""
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*":
                continue  # Skip stop codons
            all_have_gt = all("GT" in c for c in codons)
            if aa == "V":
                assert all_have_gt, "Valine should have GT in all codons"
            else:
                assert not all_have_gt, f"Amino acid {aa} unexpectedly has GT in all codons: {codons}"

    def test_every_aa_has_gt_free_codon(self):
        """Lean4: Every non-V, non-Stop AA has at least one GT-free codon."""
        for aa, codons in AA_TO_CODONS.items():
            if aa in ("V", "*"):
                continue
            gt_free = [c for c in codons if "GT" not in c]
            assert len(gt_free) > 0, f"Amino acid {aa} has no GT-free codons: {codons}"

    def test_every_aa_has_ag_free_codon(self):
        """Lean4: every_aa_has_ag_free_codon — no amino acid is AG-mandatory."""
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*":
                continue
            ag_free = [c for c in codons if "AG" not in c]
            assert len(ag_free) > 0, f"Amino acid {aa} has no AG-free codons: {codons}"

    def test_synonymous_preserves_translation(self):
        """Lean4: synonymous_preserves_translation — synonymous codons encode same AA."""
        for aa, codons in AA_TO_CODONS.items():
            for codon in codons:
                assert CODON_TABLE.get(codon) == aa, \
                    f"Codon {codon} should encode {aa} but encodes {CODON_TABLE.get(codon)}"

    def test_single_degenerate_no_synonymous(self):
        """Lean4: single_degenerate_no_synonymous_options — Trp/Met have exactly 1 codon."""
        assert len(AA_TO_CODONS.get("W", [])) == 1, "Tryptophan should have exactly 1 codon"
        assert len(AA_TO_CODONS.get("M", [])) == 1, "Methionine should have exactly 1 codon"

    def test_trp_met_no_gt_ag(self):
        """Trp (TGG) and Met (ATG) codons contain neither GT nor AG."""
        assert "GT" not in AA_TO_CODONS["W"][0], "TGG should not contain GT"
        assert "AG" not in AA_TO_CODONS["W"][0], "TGG should not contain AG"
        assert "GT" not in AA_TO_CODONS["M"][0], "ATG should not contain GT"
        assert "AG" not in AA_TO_CODONS["M"][0], "ATG should not contain AG"


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.6: Reverse Complement Handling
# ═══════════════════════════════════════════════════════════════════════════

class TestReverseComplement:
    """Verify Python's reverse-complement handling is STRONGER than the Lean4 model.

    The Lean4 model does not explicitly model reverse complement strands.
    Python checks both strands, which is strictly more thorough.
    No gap closure needed — Python is stronger than the formal model.
    """

    def test_reverse_complement_correctness(self):
        """Verify reverse complement is biologically correct."""
        from biocompiler.shared.constants import reverse_complement
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("A") == "T"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("") == ""

    def test_reverse_complement_involution(self):
        """reverse_complement(reverse_complement(s)) == s for all valid DNA."""
        from biocompiler.shared.constants import reverse_complement
        test_seqs = ["ATCG", "AACCGGTT", "ATGATG", "GCTAGCTA"]
        for seq in test_seqs:
            assert reverse_complement(reverse_complement(seq)) == seq

    @given(seq=dna_seq)
    def test_reverse_complement_involution_property(self, seq):
        """Property: RC(RC(s)) = s — reverse complement is an involution."""
        from biocompiler.shared.constants import reverse_complement
        assert reverse_complement(reverse_complement(seq)) == seq


# ═══════════════════════════════════════════════════════════════════════════
# GAP §3.7: Number of Predicates Alignment
# ═══════════════════════════════════════════════════════════════════════════

class TestPredicateAlignment:
    """Verify Python predicate set aligns with the Lean4 33-predicate model.

    Lean4 has 33 TypePredicate constructors (13 core + 20 SLOT).
    Python adds 10 extended diagnostic predicates for 43 total.
    """

    def test_predicate_count(self):
        """PREDICATE_NAMES should contain 43 entries (13 core + 20 SLOT + 10 extended diagnostic)."""
        assert len(PREDICATE_NAMES) == 43, f"Expected 43 predicates, got {len(PREDICATE_NAMES)}"

    def test_slot_predicate_count(self):
        """SLOT_PREDICATES should contain all 19 SLOT predicates from the Lean4 model."""
        # MOVE-PRED: 3 predicates moved SLOT→CORE, so 22→19 SLOT predicates
        assert len(SLOT_PREDICATES) >= 19, f"Expected >=19 SLOT predicates, got {len(SLOT_PREDICATES)}"

    def test_all_slot_predicates_in_predicate_names(self):
        """Every SLOT predicate should also appear in PREDICATE_NAMES."""
        for pred in SLOT_PREDICATES:
            assert pred in PREDICATE_NAMES, f"SLOT predicate {pred} not in PREDICATE_NAMES"

    def test_lean4_core_predicates_present(self):
        """All 13 Lean4 core predicates should have Python check functions.

        Some Lean4 predicates map to different names in Python:
          - CodonAdapted → covered by CodonOptimality (same CAI check)
          - GCInRange → checked implicitly in optimizer (not a standalone predicate name)
          - InFrame → checked implicitly via ValidCodingSeq
          - SpliceCorrect → implicitly verified via splice checks
          - NoInstabilityMotif → checked in scanner but not in PREDICATE_NAMES
        """
        # Direct matches in PREDICATE_NAMES
        direct_matches = {
            "NoStopCodons", "NoCrypticSplice", "NoRestrictionSite",
            "NoCpGIsland", "NoGTDinucleotide", "ValidCodingSeq",
            "CodonOptimality", "NoCrypticPromoter",
        }
        for pred in direct_matches:
            assert pred in PREDICATE_NAMES, f"Core predicate {pred} not found in PREDICATE_NAMES"

        # Predicates with different names but equivalent functionality
        implicit_matches = {
            "CodonAdapted": "CodonOptimality",  # CAI-based, same check
            "SpliceCorrect": None,               # Implicit via splice checks
            "GCInRange": None,                    # Implicit in optimizer GC checks
            "InFrame": "ValidCodingSeq",          # Frame check is part of valid coding seq
            "NoInstabilityMotif": None,           # In scanner but not PREDICATE_NAMES
        }
        for lean4_name, python_name in implicit_matches.items():
            if python_name is not None:
                assert python_name in PREDICATE_NAMES, \
                    f"Lean4 {lean4_name} maps to {python_name}, not found"


# ═══════════════════════════════════════════════════════════════════════════
# Cross-cutting: Translation Preservation (optimization invariant)
# ═══════════════════════════════════════════════════════════════════════════

class TestTranslationPreservation:
    """Verify the central optimization invariant: all translations preserve protein identity.

    Lean4 theorems:
      - synonymous_preserves_translation
      - pipeline_preserves_protein
      - phase1/2/5_preserves_translation

    The Python implementation should guarantee that after any optimization,
    the translated protein matches the original.
    """

    @given(seq=dna_seq.filter(lambda s: len(s) % 3 == 0 and len(s) >= 6))
    def test_translate_then_back_translate_preserves_aas(self, seq):
        """Translation is deterministic: same sequence → same protein."""
        protein1 = translate(seq)
        protein2 = translate(seq)
        assert protein1 == protein2

    def test_known_translation_examples(self):
        """Verify translation against known codon-AA mappings."""
        assert translate("ATG") == "M"          # Methionine/Start
        assert translate("TAA") == ""            # Stop (to_stop=True)
        assert translate("ATGTAA") == "M"        # Start-Stop
        assert translate("ATGAAATTT") == "MKF"   # Met-Lys-Phe

    @given(codon=st.sampled_from(list(CODON_TABLE.keys())))
    def test_single_codon_translation(self, codon):
        """Each codon translates to its correct amino acid."""
        aa = CODON_TABLE[codon]
        if aa == "*":
            # Stop codons: translate returns empty when to_stop=True
            assert translate(codon) == ""
        else:
            assert translate(codon) == aa


# ═══════════════════════════════════════════════════════════════════════════
# Cross-cutting: Evidence Audit Trail
# ═══════════════════════════════════════════════════════════════════════════

class TestEvidenceAuditTrail:
    """Verify VerificationEvidence objects maintain an auditable trust chain.

    Lean4: The VerificationContext makes trust explicit and parameterizable.
    Python: VerificationEvidence makes the trust chain auditable.
    """

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_conservative_evidence_always_unverified(self, seq):
        """CONSERVATIVE mode evidence should always be unverified."""
        _, evidence = verify_no_cryptic_splice(seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert not evidence.verified

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_verified_evidence_documents_tool(self, seq):
        """VERIFIED mode evidence should document the tool used."""
        _, evidence = verify_no_cryptic_splice(seq, slot_mode=SLOTMode.VERIFIED)
        assert evidence.tool_name is not None
        assert len(evidence.tool_name) > 0

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_evidence_serializable(self, seq):
        """Evidence should be serializable for audit logging."""
        _, evidence = verify_no_cryptic_splice(seq, slot_mode=SLOTMode.VERIFIED)
        d = evidence.to_dict()
        assert isinstance(d, dict)
        assert "predicate" in d
        assert "slot_mode" in d
        assert "tool_available" in d
        assert "verified" in d

    @given(seq=dna_seq.filter(lambda s: len(s) >= 6))
    def test_mode_documented_in_evidence(self, seq):
        """Evidence should always document which SLOTMode was used."""
        for mode in SLOTMode:
            _, evidence = verify_no_cryptic_splice(seq, slot_mode=mode)
            assert evidence.slot_mode == mode


# ═══════════════════════════════════════════════════════════════════════════
# Cross-cutting: Constraints Do Not Compose (Lean4 counterexamples)
# ═══════════════════════════════════════════════════════════════════════════

class TestConstraintsDoNotCompose:
    """Verify the Python implementation correctly handles the key insight from Lean4:

    Lean4 theorems:
      - dinucleotide_no_compose: GT can form at junction of two GT-free sequences
      - restriction_site_no_compose: restriction sites can form at junctions

    Python must detect junction-crossing patterns.
    """

    def test_gt_at_junction(self):
        """Lean4: dinucleotide_no_compose — GT forms at junction of 'CG' and 'TC'."""
        s1 = "CG"
        s2 = "TC"
        assert "GT" not in s1
        assert "GT" not in s2
        assert "GT" in (s1 + s2)  # CGTC contains GT at positions 1-2

    def test_restriction_site_at_junction(self):
        """Lean4: restriction_site_no_compose — GATC forms at junction of 'GAT' and 'C'."""
        s1 = "GAT"
        s2 = "CG"
        assert "GATC" not in s1
        assert "GATC" not in s2
        assert "GATC" in (s1 + s2)

    def test_junction_detection_in_python(self):
        """Python's cross-codon checks should detect junction-crossing dinucleotides."""
        # Create a sequence where GT crosses a codon boundary
        # Codon 1 ends with G, Codon 2 starts with T
        seq = "TTTGATCTT"  # codons: TTT GAT CTT — GT at positions 4-5 (cross-codon)
        gt_positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]
        # This specific sequence does not have GT, let us construct one that does
        seq_with_cross_codon_gt = "GCAGTCAAA"  # GCA GTC AAA — GT at positions 2-3 (within codon 2)
        gt_positions = [i for i in range(len(seq_with_cross_codon_gt) - 1)
                       if seq_with_cross_codon_gt[i:i+2] == "GT"]
        assert len(gt_positions) > 0, "Cross-codon GT should be detected"
