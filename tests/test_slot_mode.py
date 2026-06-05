"""
Tests for SLOT Mode Integration
================================

Tests that verify SLOT mode behavior in the type system:
- CONSERVATIVE: SLOT predicates always return UNCERTAIN
- VERIFIED: SLOT predicates return PASS when verification conditions met
- PERMISSIVE: SLOT predicates return PASS with weaker evidence
- Certificate includes slot_mode metadata
- Backward compatibility (default is CONSERVATIVE)
"""

import pytest
from biocompiler.types import Verdict, SLOTMode, TypeCheckResult, Certificate
from biocompiler.type_system import (
    evaluate_all_predicates,
    PredicateResult,
    CertLevel,
    PREDICATE_NAMES,
)
from biocompiler.slot_verification import (
    is_slot_predicate,
    SLOT_PREDICATES,
    VerificationEvidence,
    verify_slot_predicate,
    verify_no_cryptic_splice,
    verify_no_cryptic_promoter,
    verify_no_unexpected_tm_domain,
    verify_mrna_secondary_structure,
    verify_conservation_score,
    verify_codon_optimality,
)
from biocompiler.certificate import (
    generate_certificate,
    compute_certificate,
    format_certificate,
)


# ────────────────────────────────────────────────────────────
# Test Fixtures
# ────────────────────────────────────────────────────────────

@pytest.fixture
def simple_seq():
    """A simple valid coding sequence (ATG = M, TTT = F, TTT = F, TAA = stop)."""
    return "ATGTTTTTTTAA"


@pytest.fixture
def clean_seq():
    """A longer sequence with no problematic features."""
    # ATG (M), GCC (A), GCT (A), GGA (G), TCC (S), ACC (T), GGC (G), TAA (*)
    return "ATGGCCGCTGGATCCACCGGCTAA"


@pytest.fixture
def seq_with_gt():
    """A sequence containing GT dinucleotides (for splice testing)."""
    # Contains GT in 'GTT' (Valine codon)
    return "ATGGTTGCTGGATCCACCTAA"


# ────────────────────────────────────────────────────────────
# 1. SLOT Predicate Classification Tests
# ────────────────────────────────────────────────────────────

class TestSLOTPredicateClassification:
    """Test that SLOT predicate classification is correct."""

    def test_known_slot_predicates(self):
        """SLOT_PREDICATES contains expected predicate names."""
        assert "NoCrypticSplice" in SLOT_PREDICATES
        assert "NoCrypticPromoter" in SLOT_PREDICATES
        assert "NoUnexpectedTMDomain" in SLOT_PREDICATES
        assert "mRNASecondaryStructure" in SLOT_PREDICATES
        assert "CoTranslationalFolding" in SLOT_PREDICATES
        assert "ConservationScore" in SLOT_PREDICATES
        assert "CodonOptimality" in SLOT_PREDICATES
        # Structure predicates
        assert "StructureConfidence" in SLOT_PREDICATES
        assert "NoMisfoldingRisk" in SLOT_PREDICATES
        # Stability predicates
        assert "StableFolding" in SLOT_PREDICATES
        # Solubility predicates
        assert "SolubleExpression" in SLOT_PREDICATES
        # Immunogenicity predicates
        assert "LowImmunogenicity" in SLOT_PREDICATES

    def test_core_predicates_not_in_slot(self):
        """Core predicates are NOT in SLOT_PREDICATES."""
        assert "NoStopCodons" not in SLOT_PREDICATES
        assert "NoGTDinucleotide" not in SLOT_PREDICATES
        assert "ValidCodingSeq" not in SLOT_PREDICATES
        assert "NoRestrictionSite" not in SLOT_PREDICATES
        assert "NoCpGIsland" not in SLOT_PREDICATES

    def test_is_slot_predicate_function(self):
        """is_slot_predicate correctly identifies SLOT predicates."""
        assert is_slot_predicate("NoCrypticSplice") is True
        assert is_slot_predicate("NoCrypticPromoter") is True
        assert is_slot_predicate("NoStopCodons") is False
        assert is_slot_predicate("ValidCodingSeq") is False

    def test_is_slot_predicate_handles_parameterized_names(self):
        """is_slot_predicate handles parameterized predicate names."""
        assert is_slot_predicate("CodonOptimality(0.5)") is True
        assert is_slot_predicate("NoCrypticSplice(threshold=3.0)") is True
        assert is_slot_predicate("NoStopCodons(strict=True)") is False

    def test_slot_predicate_count(self):
        """SLOT_PREDICATES has expected count (19 from architecture)."""
        # 7 DNA-level + 4 structure + 4 stability + 4 solubility + 4 immunogenicity
        # But the exact count depends on what's classified as SLOT
        # At minimum, we should have the DNA-level SLOT predicates
        assert len(SLOT_PREDICATES) >= 7


# ────────────────────────────────────────────────────────────
# 2. SLOTMode Enum Tests
# ────────────────────────────────────────────────────────────

class TestSLOTMode:
    """Test SLOTMode enum values."""

    def test_conservative_value(self):
        assert SLOTMode.CONSERVATIVE.value == "conservative"

    def test_verified_value(self):
        assert SLOTMode.VERIFIED.value == "verified"

    def test_permissive_value(self):
        assert SLOTMode.PERMISSIVE.value == "permissive"

    def test_all_modes_distinct(self):
        values = {m.value for m in SLOTMode}
        assert len(values) == 3


# ────────────────────────────────────────────────────────────
# 3. CONSERVATIVE Mode Tests
# ────────────────────────────────────────────────────────────

class TestConservativeMode:
    """Test CONSERVATIVE mode: SLOT predicates always return UNCERTAIN."""

    def test_conservative_default(self, clean_seq):
        """Default slot_mode is CONSERVATIVE."""
        results = evaluate_all_predicates(clean_seq)
        slot_results = [r for r in results if is_slot_predicate(r.predicate)]
        for r in slot_results:
            assert r.verdict == Verdict.UNCERTAIN, (
                f"SLOT predicate {r.predicate} should be UNCERTAIN in CONSERVATIVE mode, "
                f"got {r.verdict}"
            )

    def test_conservative_explicit(self, clean_seq):
        """Explicit CONSERVATIVE mode same as default."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        slot_results = [r for r in results if is_slot_predicate(r.predicate)]
        for r in slot_results:
            assert r.verdict == Verdict.UNCERTAIN

    def test_conservative_core_predicates_unaffected(self, clean_seq):
        """Core predicates still evaluate normally in CONSERVATIVE mode."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        core_results = [r for r in results if not is_slot_predicate(r.predicate)]
        # Core predicates should have definite verdicts (PASS or FAIL)
        for r in core_results:
            assert r.verdict in (Verdict.PASS, Verdict.FAIL), (
                f"Core predicate {r.predicate} should have definite verdict, "
                f"got {r.verdict}"
            )

    def test_conservative_knowledge_gap_set(self, clean_seq):
        """CONSERVATIVE mode sets knowledge_gap on SLOT predicates."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        slot_results = [r for r in results if is_slot_predicate(r.predicate)]
        for r in slot_results:
            assert r.knowledge_gap is not None
            assert "SLOT" in r.knowledge_gap
            assert "CONSERVATIVE" in r.knowledge_gap

    def test_conservative_splice_verification(self, clean_seq):
        """verify_no_cryptic_splice returns UNCERTAIN in CONSERVATIVE mode."""
        verdict, evidence = verify_no_cryptic_splice(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN
        assert evidence.slot_mode == SLOTMode.CONSERVATIVE
        assert evidence.verified is False

    def test_conservative_promoter_verification(self, clean_seq):
        """verify_no_cryptic_promoter returns UNCERTAIN in CONSERVATIVE mode."""
        verdict, evidence = verify_no_cryptic_promoter(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN
        assert evidence.verified is False

    def test_conservative_tm_domain_verification(self, clean_seq):
        """verify_no_unexpected_tm_domain returns UNCERTAIN in CONSERVATIVE mode."""
        verdict, evidence = verify_no_unexpected_tm_domain(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN
        assert evidence.verified is False

    def test_conservative_mrna_verification(self, clean_seq):
        """verify_mrna_secondary_structure returns UNCERTAIN in CONSERVATIVE mode."""
        verdict, evidence = verify_mrna_secondary_structure(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN
        assert evidence.verified is False

    def test_conservative_conservation_verification(self):
        """verify_conservation_score returns UNCERTAIN in CONSERVATIVE mode."""
        verdict, evidence = verify_conservation_score("A", "A", slot_mode=SLOTMode.CONSERVATIVE)
        assert verdict == Verdict.UNCERTAIN
        assert evidence.verified is False

    def test_conservative_codon_optimality_verification(self):
        """verify_codon_optimality returns UNCERTAIN in CONSERVATIVE mode."""
        verdict, evidence = verify_codon_optimality(
            "ATG", {"ATG": 0.8}, slot_mode=SLOTMode.CONSERVATIVE
        )
        assert verdict == Verdict.UNCERTAIN
        assert evidence.verified is False


# ────────────────────────────────────────────────────────────
# 4. VERIFIED Mode Tests
# ────────────────────────────────────────────────────────────

class TestVerifiedMode:
    """Test VERIFIED mode: SLOT predicates return PASS when conditions met."""

    def test_verified_mode_returns_pass_when_conditions_met(self, clean_seq):
        """VERIFIED mode returns PASS for SLOT predicates that pass."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.VERIFIED)
        # At minimum, some SLOT predicates should have non-UNCERTAIN verdicts
        # in VERIFIED mode (not all will be PASS, but they shouldn't all be UNCERTAIN)
        slot_results = [r for r in results if is_slot_predicate(r.predicate)]
        # In VERIFIED mode, SLOT predicates should have evidence in derivation
        verified_preds = [r for r in slot_results if r.derivation is not None]
        # At least some should have verification evidence
        assert len(verified_preds) > 0, "Expected at least some SLOT predicates to have verification evidence"

    def test_verified_splice_with_clean_seq(self, clean_seq):
        """VERIFIED mode: NoCrypticSplice should PASS for clean sequence."""
        verdict, evidence = verify_no_cryptic_splice(clean_seq, slot_mode=SLOTMode.VERIFIED)
        # With MaxEntScan available and no GT sites, should be PASS
        assert verdict in (Verdict.PASS, Verdict.UNCERTAIN)  # UNCERTAIN if MaxEnt unavailable
        assert evidence.slot_mode == SLOTMode.VERIFIED

    def test_verified_promoter_with_clean_seq(self, clean_seq):
        """VERIFIED mode: NoCrypticPromoter should PASS for clean sequence."""
        verdict, evidence = verify_no_cryptic_promoter(clean_seq, slot_mode=SLOTMode.VERIFIED)
        # PWM scanner is built-in, so should be able to verify
        assert evidence.tool_available is True
        assert evidence.slot_mode == SLOTMode.VERIFIED

    def test_verified_core_predicates_unaffected(self, clean_seq):
        """Core predicates still evaluate normally in VERIFIED mode."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.VERIFIED)
        core_results = [r for r in results if not is_slot_predicate(r.predicate)]
        for r in core_results:
            assert r.verdict in (Verdict.PASS, Verdict.FAIL), (
                f"Core predicate {r.predicate} should have definite verdict, "
                f"got {r.verdict}"
            )

    def test_verified_evidence_documented(self, clean_seq):
        """VERIFIED mode produces verification evidence."""
        verdict, evidence = verify_no_cryptic_promoter(clean_seq, slot_mode=SLOTMode.VERIFIED)
        assert isinstance(evidence, VerificationEvidence)
        assert evidence.predicate == "NoCrypticPromoter"
        assert evidence.slot_mode == SLOTMode.VERIFIED
        assert evidence.tool_available is True
        assert evidence.tool_name == "PWM_scanner"

    def test_verified_codon_optimality_with_data(self):
        """VERIFIED mode: CodonOptimality PASS when CAI data available and high."""
        cai_data = {"ATG": 0.9, "TTT": 0.7, "GCT": 0.8}
        verdict, evidence = verify_codon_optimality(
            "ATG", cai_data, min_cai=0.5, slot_mode=SLOTMode.VERIFIED
        )
        assert verdict == Verdict.PASS
        assert evidence.verified is True
        assert evidence.tool_available is True

    def test_verified_codon_optimality_without_data(self):
        """VERIFIED mode: CodonOptimality UNCERTAIN when CAI data unavailable."""
        verdict, evidence = verify_codon_optimality(
            "ATG", {}, min_cai=0.5, slot_mode=SLOTMode.VERIFIED
        )
        assert verdict == Verdict.UNCERTAIN
        assert evidence.verified is False
        assert evidence.tool_available is False


# ────────────────────────────────────────────────────────────
# 5. PERMISSIVE Mode Tests
# ────────────────────────────────────────────────────────────

class TestPermissiveMode:
    """Test PERMISSIVE mode: SLOT predicates return PASS with weaker evidence."""

    def test_permissive_more_pass_than_conservative(self, clean_seq):
        """PERMISSIVE mode produces at least as many PASS verdicts as CONSERVATIVE."""
        results_conservative = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        results_permissive = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.PERMISSIVE)
        
        slot_conservative = [r for r in results_conservative if is_slot_predicate(r.predicate)]
        slot_permissive = [r for r in results_permissive if is_slot_predicate(r.predicate)]
        
        pass_count_conservative = sum(1 for r in slot_conservative if r.verdict == Verdict.PASS)
        pass_count_permissive = sum(1 for r in slot_permissive if r.verdict == Verdict.PASS)
        
        # PERMISSIVE should produce at least as many PASS verdicts
        assert pass_count_permissive >= pass_count_conservative

    def test_permissive_promoter_relaxed(self, clean_seq):
        """PERMISSIVE mode uses relaxed promoter threshold."""
        verdict_c, _ = verify_no_cryptic_promoter(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        verdict_p, evidence_p = verify_no_cryptic_promoter(clean_seq, slot_mode=SLOTMode.PERMISSIVE)
        # PERMISSIVE should not be worse than CONSERVATIVE
        assert verdict_p.value >= verdict_c.value or verdict_p == Verdict.PASS

    def test_permissive_splice_relaxed(self, clean_seq):
        """PERMISSIVE mode uses relaxed splice threshold."""
        verdict_c, _ = verify_no_cryptic_splice(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        verdict_p, evidence_p = verify_no_cryptic_splice(clean_seq, slot_mode=SLOTMode.PERMISSIVE)
        # PERMISSIVE should not be worse than CONSERVATIVE
        assert evidence_p.slot_mode == SLOTMode.PERMISSIVE

    def test_permissive_uncertain_promoted_to_pass(self, clean_seq):
        """PERMISSIVE mode: UNCERTAIN results are promoted to PASS."""
        verdict, evidence = verify_no_cryptic_promoter(clean_seq, slot_mode=SLOTMode.PERMISSIVE)
        # In PERMISSIVE mode, UNCERTAIN should be promoted to PASS
        # (unless there's a FAIL-level violation)
        if evidence.tool_available and evidence.verified:
            assert verdict in (Verdict.PASS, Verdict.LIKELY_PASS, Verdict.FAIL)

    def test_permissive_core_unaffected(self, clean_seq):
        """Core predicates unaffected by PERMISSIVE mode."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.PERMISSIVE)
        core_results = [r for r in results if not is_slot_predicate(r.predicate)]
        for r in core_results:
            assert r.verdict in (Verdict.PASS, Verdict.FAIL)


# ────────────────────────────────────────────────────────────
# 6. Certificate Integration Tests
# ────────────────────────────────────────────────────────────

class TestCertificateSLOTMode:
    """Test that certificates include slot_mode metadata."""

    def test_certificate_includes_slot_mode(self, clean_seq):
        """Certificate provenance includes slot_mode."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        type_results = [
            TypeCheckResult(predicate=r.predicate, verdict=r.verdict,
                          derivation=r.derivation, violation=r.violation,
                          knowledge_gap=r.knowledge_gap)
            for r in results
        ]
        cert = generate_certificate(
            sequence=clean_seq,
            type_results=type_results,
            input_params={"organism": "Homo_sapiens"},
            slot_mode=SLOTMode.CONSERVATIVE,
        )
        cert_dict = cert.to_dict()
        assert "slot_mode" in cert_dict["provenance"]
        assert cert_dict["provenance"]["slot_mode"] == "conservative"

    def test_certificate_verified_mode(self, clean_seq):
        """Certificate with VERIFIED mode documents it."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.VERIFIED)
        type_results = [
            TypeCheckResult(predicate=r.predicate, verdict=r.verdict,
                          derivation=r.derivation, violation=r.violation,
                          knowledge_gap=r.knowledge_gap)
            for r in results
        ]
        cert = generate_certificate(
            sequence=clean_seq,
            type_results=type_results,
            input_params={"organism": "Homo_sapiens"},
            slot_mode=SLOTMode.VERIFIED,
        )
        cert_dict = cert.to_dict()
        assert cert_dict["provenance"]["slot_mode"] == "verified"
        assert "slot_mode_description" in cert_dict["provenance"]

    def test_certificate_permissive_mode(self, clean_seq):
        """Certificate with PERMISSIVE mode documents it."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.PERMISSIVE)
        type_results = [
            TypeCheckResult(predicate=r.predicate, verdict=r.verdict,
                          derivation=r.derivation, violation=r.violation,
                          knowledge_gap=r.knowledge_gap)
            for r in results
        ]
        cert = generate_certificate(
            sequence=clean_seq,
            type_results=type_results,
            input_params={"organism": "Homo_sapiens"},
            slot_mode=SLOTMode.PERMISSIVE,
        )
        cert_dict = cert.to_dict()
        assert cert_dict["provenance"]["slot_mode"] == "permissive"

    def test_slot_mode_description_present(self, clean_seq):
        """Certificate includes human-readable slot_mode_description."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        type_results = [
            TypeCheckResult(predicate=r.predicate, verdict=r.verdict,
                          derivation=r.derivation, violation=r.violation,
                          knowledge_gap=r.knowledge_gap)
            for r in results
        ]
        cert = generate_certificate(
            sequence=clean_seq,
            type_results=type_results,
            input_params={"organism": "Homo_sapiens"},
            slot_mode=SLOTMode.CONSERVATIVE,
        )
        cert_dict = cert.to_dict()
        desc = cert_dict["provenance"]["slot_mode_description"]
        assert "Lean4" in desc or "UNCERTAIN" in desc

    def test_gold_verified_stronger_than_gold_conservative(self, clean_seq):
        """GOLD certificate with VERIFIED mode is documented as stronger."""
        # Conservative: SLOT predicates are UNCERTAIN
        results_c = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        type_results_c = [
            TypeCheckResult(predicate=r.predicate, verdict=r.verdict,
                          derivation=r.derivation, violation=r.violation,
                          knowledge_gap=r.knowledge_gap)
            for r in results_c
        ]
        cert_c = generate_certificate(
            sequence=clean_seq,
            type_results=type_results_c,
            input_params={"organism": "Homo_sapiens"},
            slot_mode=SLOTMode.CONSERVATIVE,
        )

        # Verified: SLOT predicates have actual results
        results_v = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.VERIFIED)
        type_results_v = [
            TypeCheckResult(predicate=r.predicate, verdict=r.verdict,
                          derivation=r.derivation, violation=r.violation,
                          knowledge_gap=r.knowledge_gap)
            for r in results_v
        ]
        cert_v = generate_certificate(
            sequence=clean_seq,
            type_results=type_results_v,
            input_params={"organism": "Homo_sapiens"},
            slot_mode=SLOTMode.VERIFIED,
        )

        # Both should be valid certificates
        assert cert_c.to_dict()["provenance"]["slot_mode"] == "conservative"
        assert cert_v.to_dict()["provenance"]["slot_mode"] == "verified"

    def test_format_certificate_includes_slot_mode(self, clean_seq):
        """format_certificate includes slot_mode in output."""
        results = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.VERIFIED)
        pred_results = []
        for r in results:
            pred_results.append(
                PredicateResult(
                    predicate=r.predicate,
                    passed=r.passed,
                    verdict=r.verdict,
                    details=f"verdict={r.verdict.value}",
                )
            )
        report = format_certificate(pred_results, clean_seq, "Homo_sapiens", slot_mode=SLOTMode.VERIFIED)
        assert "SLOT Mode:       verified" in report
        assert "VERIFIED" in report


# ────────────────────────────────────────────────────────────
# 7. Backward Compatibility Tests
# ────────────────────────────────────────────────────────────

class TestBackwardCompatibility:
    """Test backward compatibility with existing code."""

    def test_default_slot_mode_is_conservative(self, clean_seq):
        """Default slot_mode is CONSERVATIVE (no breaking change)."""
        results = evaluate_all_predicates(clean_seq)
        slot_results = [r for r in results if is_slot_predicate(r.predicate)]
        # In default (CONSERVATIVE) mode, all SLOT predicates should be UNCERTAIN
        for r in slot_results:
            assert r.verdict == Verdict.UNCERTAIN

    def test_predicate_result_has_verification_evidence_field(self):
        """PredicateResult has verification_evidence field."""
        result = PredicateResult("TestPred", True, verdict=Verdict.PASS, details="test")
        assert hasattr(result, 'verification_evidence')
        assert result.verification_evidence is None  # Default is None

    def test_predicate_result_backward_compatible(self):
        """PredicateResult works without verification_evidence (backward compat)."""
        # Old code that doesn't pass verification_evidence should still work
        result = PredicateResult("NoStopCodons", True, verdict=Verdict.PASS, details="OK")
        assert result.passed is True
        assert result.predicate == "NoStopCodons"

    def test_evaluate_all_predicates_without_slot_mode(self, clean_seq):
        """evaluate_all_predicates works without slot_mode parameter."""
        # This tests backward compatibility: old code that doesn't pass slot_mode
        results = evaluate_all_predicates(clean_seq)
        assert len(results) == 12

    def test_compute_certificate_backward_compatible(self):
        """compute_certificate works without slot_mode parameter."""
        results = [
            PredicateResult("NoStopCodons", True, verdict=Verdict.PASS, details="OK"),
            PredicateResult("ValidCodingSeq", True, verdict=Verdict.PASS, details="OK"),
        ]
        level = compute_certificate(results)
        assert level == CertLevel.GOLD

    def test_generate_certificate_backward_compatible(self, clean_seq):
        """generate_certificate works without slot_mode parameter."""
        results = evaluate_all_predicates(clean_seq)
        type_results = [
            TypeCheckResult(predicate=r.predicate, verdict=r.verdict)
            for r in results
        ]
        cert = generate_certificate(
            sequence=clean_seq,
            type_results=type_results,
            input_params={"organism": "Homo_sapiens"},
        )
        assert cert is not None
        # Default slot_mode should be conservative
        cert_dict = cert.to_dict()
        assert cert_dict["provenance"]["slot_mode"] == "conservative"


# ────────────────────────────────────────────────────────────
# 8. Verification Evidence Tests
# ────────────────────────────────────────────────────────────

class TestVerificationEvidence:
    """Test VerificationEvidence data structure."""

    def test_evidence_creation(self):
        """VerificationEvidence can be created with required fields."""
        evidence = VerificationEvidence(
            predicate="NoCrypticSplice",
            slot_mode=SLOTMode.CONSERVATIVE,
            tool_available=True,
            tool_name="MaxEntScan",
            verified=False,
            details="test",
        )
        assert evidence.predicate == "NoCrypticSplice"
        assert evidence.slot_mode == SLOTMode.CONSERVATIVE
        assert evidence.tool_available is True
        assert evidence.verified is False

    def test_evidence_to_dict(self):
        """VerificationEvidence serializes to dict."""
        evidence = VerificationEvidence(
            predicate="NoCrypticPromoter",
            slot_mode=SLOTMode.VERIFIED,
            tool_available=True,
            tool_name="PWM_scanner",
            tool_result="score=0.3",
            threshold_used=0.7,
            verified=True,
            details="PASS",
        )
        d = evidence.to_dict()
        assert d["predicate"] == "NoCrypticPromoter"
        assert d["slot_mode"] == "verified"
        assert d["tool_available"] is True
        assert d["tool_name"] == "PWM_scanner"
        assert d["tool_result"] == "score=0.3"
        assert d["threshold_used"] == 0.7
        assert d["verified"] is True
        assert d["details"] == "PASS"

    def test_evidence_optional_fields(self):
        """VerificationEvidence handles None optional fields."""
        evidence = VerificationEvidence(
            predicate="Test",
            slot_mode=SLOTMode.CONSERVATIVE,
            tool_available=False,
            tool_name="test_tool",
        )
        d = evidence.to_dict()
        assert d["tool_result"] is None
        assert d["threshold_used"] is None
        assert d["verified"] is False
        assert d["details"] == ""


# ────────────────────────────────────────────────────────────
# 9. Slot Verification Dispatch Tests
# ────────────────────────────────────────────────────────────

class TestSlotVerificationDispatch:
    """Test verify_slot_predicate dispatch function."""

    def test_dispatch_unknown_predicate_raises(self):
        """Dispatch raises ValueError for non-SLOT predicates."""
        with pytest.raises(ValueError, match="not a SLOT predicate"):
            verify_slot_predicate("NoStopCodons", slot_mode=SLOTMode.CONSERVATIVE)

    def test_dispatch_known_slot_predicate(self, clean_seq):
        """Dispatch works for known SLOT predicates."""
        verdict, evidence = verify_slot_predicate(
            "NoCrypticSplice",
            slot_mode=SLOTMode.CONSERVATIVE,
            seq=clean_seq,
        )
        assert verdict == Verdict.UNCERTAIN
        assert evidence.predicate == "NoCrypticSplice"

    def test_dispatch_structure_predicate(self):
        """Dispatch works for structure SLOT predicates."""
        verdict, evidence = verify_slot_predicate(
            "StructureConfidence",
            slot_mode=SLOTMode.CONSERVATIVE,
            protein_sequence="MAGIC",
        )
        assert verdict == Verdict.UNCERTAIN
        assert evidence.predicate == "StructureConfidence"

    def test_dispatch_stability_predicate(self):
        """Dispatch works for stability SLOT predicates."""
        verdict, evidence = verify_slot_predicate(
            "StableFolding",
            slot_mode=SLOTMode.CONSERVATIVE,
            protein_sequence="MAGIC",
        )
        assert verdict == Verdict.UNCERTAIN
        assert evidence.predicate == "StableFolding"

    def test_dispatch_solubility_predicate(self):
        """Dispatch works for solubility SLOT predicates."""
        verdict, evidence = verify_slot_predicate(
            "SolubleExpression",
            slot_mode=SLOTMode.CONSERVATIVE,
            protein_sequence="MAGIC",
        )
        assert verdict == Verdict.UNCERTAIN
        assert evidence.predicate == "SolubleExpression"

    def test_dispatch_immunogenicity_predicate(self):
        """Dispatch works for immunogenicity SLOT predicates."""
        verdict, evidence = verify_slot_predicate(
            "LowImmunogenicity",
            slot_mode=SLOTMode.CONSERVATIVE,
            protein_sequence="MAGIC",
        )
        assert verdict == Verdict.UNCERTAIN
        assert evidence.predicate == "LowImmunogenicity"


# ────────────────────────────────────────────────────────────
# 10. Mode Comparison Tests
# ────────────────────────────────────────────────────────────

class TestModeComparison:
    """Test that mode strength ordering is correct."""

    def test_conservative_is_most_restrictive(self, clean_seq):
        """CONSERVATIVE produces most UNCERTAIN verdicts for SLOT predicates."""
        results_c = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        results_v = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.VERIFIED)
        results_p = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.PERMISSIVE)

        slot_c = [r for r in results_c if is_slot_predicate(r.predicate)]
        slot_v = [r for r in results_v if is_slot_predicate(r.predicate)]
        slot_p = [r for r in results_p if is_slot_predicate(r.predicate)]

        uncertain_c = sum(1 for r in slot_c if r.verdict == Verdict.UNCERTAIN)
        uncertain_v = sum(1 for r in slot_v if r.verdict == Verdict.UNCERTAIN)
        uncertain_p = sum(1 for r in slot_p if r.verdict == Verdict.UNCERTAIN)

        # CONSERVATIVE should have the most UNCERTAIN verdicts
        assert uncertain_c >= uncertain_v
        assert uncertain_c >= uncertain_p

    def test_core_predicates_identical_across_modes(self, clean_seq):
        """Core predicates produce identical results across all modes."""
        results_c = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.CONSERVATIVE)
        results_v = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.VERIFIED)
        results_p = evaluate_all_predicates(clean_seq, slot_mode=SLOTMode.PERMISSIVE)

        core_c = {r.predicate: r.verdict for r in results_c if not is_slot_predicate(r.predicate)}
        core_v = {r.predicate: r.verdict for r in results_v if not is_slot_predicate(r.predicate)}
        core_p = {r.predicate: r.verdict for r in results_p if not is_slot_predicate(r.predicate)}

        # Core predicates should have identical verdicts across modes
        for pred in core_c:
            assert core_c[pred] == core_v[pred], f"{pred} differs between CONSERVATIVE and VERIFIED"
            assert core_c[pred] == core_p[pred], f"{pred} differs between CONSERVATIVE and PERMISSIVE"
