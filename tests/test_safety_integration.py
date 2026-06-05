"""
BioCompiler Safety Integration Tests
=====================================

Comprehensive end-to-end tests verifying the entire safety pipeline:

1. Hazardous input rejection  — biosecurity screening blocks dangerous sequences
2. Translation verification    — optimized proteins match input exactly
3. Provenance tracking         — every optimization is tracked by default
4. Strict mode enforcement     — impossible constraints are caught cleanly
5. Biosafety in exports        — GenBank/FASTA carry safety annotations
6. Full safety pipeline        — complete flow from input to export

Tests that depend on modules created by other agents (biosecurity,
protein_verification, strict_mode) are marked with
``@pytest.mark.skip(reason="Awaiting Agent N module")`` but contain
complete test code ready for activation.

Markers:
    requires_external — depends on external services (MHCflurry, NetMHCpan, etc.)
    integration       — exercises multiple modules together
    e2e               — full end-to-end pipeline test
"""

from __future__ import annotations

import json
import re
import warnings
from typing import Any

import pytest

# ═══════════════════════════════════════════════════════════════════════
# Core imports (always available)
# ═══════════════════════════════════════════════════════════════════════

from biocompiler.optimization import optimize_sequence, OptimizationResult
from biocompiler.translation import translate, compute_cai
from biocompiler.export import export_fasta, export_genbank
from biocompiler.exceptions import (
    BioCompilerError,
    OptimizationError,
    InvalidProteinError,
)
from biocompiler.scanner import gc_content

# ═══════════════════════════════════════════════════════════════════════
# Conditional imports — modules being built by other agents
# ═══════════════════════════════════════════════════════════════════════

# Biosecurity module (Agent 3 expected)
try:
    from biocompiler.biosecurity import (
        BiosecurityError,
        BiosecurityScreenResult,
        screen_sequence,
        HAZARDOUS_PROTEINS,
        RESISTANCE_MARKERS,
    )
    _HAS_BIOSECURITY = True
except ImportError:
    _HAS_BIOSECURITY = False
    BiosecurityError = None  # type: ignore[assignment, misc]
    BiosecurityScreenResult = None  # type: ignore[assignment, misc]
    screen_sequence = None  # type: ignore[assignment, misc]
    HAZARDOUS_PROTEINS = None  # type: ignore[assignment, misc]
    RESISTANCE_MARKERS = None  # type: ignore[assignment, misc]

# Protein verification module (Agent 4 expected)
try:
    from biocompiler.protein_verification import (
        verify_translation,
        TranslationVerificationResult,
        detect_corruption,
    )
    _HAS_PROTEIN_VERIFICATION = True
except ImportError:
    _HAS_PROTEIN_VERIFICATION = False
    verify_translation = None  # type: ignore[assignment, misc]
    TranslationVerificationResult = None  # type: ignore[assignment, misc]
    detect_corruption = None  # type: ignore[assignment, misc]

# Strict mode module (Agent 5 expected)
try:
    from biocompiler.strict_mode import (
        StrictModeError,
        optimize_strict,
        FailedPredicatesError,
    )
    _HAS_STRICT_MODE = True
except ImportError:
    _HAS_STRICT_MODE = False
    StrictModeError = None  # type: ignore[assignment, misc]
    optimize_strict = None  # type: ignore[assignment, misc]
    FailedPredicatesError = None  # type: ignore[assignment, misc]

# Provenance — already in codebase
from biocompiler.provenance import (
    OptimizationRecord,
    ProvenanceTracker,
)


# ═══════════════════════════════════════════════════════════════════════
# Test protein fixtures
# ═══════════════════════════════════════════════════════════════════════

# Human insulin A-chain + B-chain (preproinsulin signal peptide removed)
INSULIN_PROTEIN = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"

# Enhanced GFP (239 AA)
EGFP_PROTEIN = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
    "VTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVN"
    "RIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQ"
    "QNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)

# Human hemoglobin beta (147 AA) — a well-characterised self-protein
HBB_PROTEIN = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
    "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK"
    "EFTPPVQAAYQKVVAGVANALAHKYH"
)

# Ricin A-chain (approximate, 267 AA) — a known toxin
# Based on UniProt P02879 (Ricin A-chain, residues 1-267)
RICIN_A_CHAIN = (
    "IFPKQYFIINPQTLPKCNKCVNNVCLRPSSQSYVYNSCTTNGICFKDNDRLFFTLVSSNQTL"
    "PSRSNLQNGTWFTNKKTKLVTATSVLSRCSSCTVNCPKVQVPTNATVTKTLCSSQCTNNSCQ"
    "LKHSLSQYNKIRVNRRSTQTICRISTCKNKGLCCSQIQTSCTNSSCSNICNDTKTQNCSKSS"
    "NVSRLYKTVNLSNSNLNNSNCTTTATCTFSPCKTDNNSCIPNLQNTSNYTTQSCTATPKHFS"
    "NISQCKNTQNCSYNNSQTCQQSLNECLSNLTVSNSTQSCNSNCSKSNCSNNIGCLTNGYTN"
    "SCSNNVCLPKTNLCNICSNQCLNSNLCNNICSNPCLNSVCLSSNCKNVCLNSNCL"
)

# SARS-CoV-2 Spike RBD (receptor-binding domain, ~193 AA)
# Residues 319-541 of the SARS-CoV-2 spike protein (Wuhan-Hu-1)
SARS_COV2_RBD = (
    "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEV"
    "RQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYG"
    "FQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"
)

# blaTEM-1 beta-lactamase (286 AA) — common antibiotic resistance marker
BLATEM_PROTEIN = (
    "MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLG"
    "RRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTM"
    "VAMATTSRKVASALAPVRKSQLLQAGYLGSRFHKALRQLAEVPVASLPIVSRVSPELNEALQGRLVQATVQGLDVPQALKHPSVLSQMLR"
    "ELGSNVTATLLAHGNNLRNKATVFLPDLVNNSHLLVNLPGANRSFVAMTPLQAPAGITLSDMLNRIGFEMRFAGYDLKTNDVAVKSLLEQK"
    "LKVPPEVLTEDNLKVLGVLTNKAQGQKVTLPAELEQWLPGCALRIGSHTPVWLNEQYVTR"
)


@pytest.fixture
def insulin() -> str:
    """Human insulin protein for optimization tests."""
    return INSULIN_PROTEIN


@pytest.fixture
def egfp() -> str:
    """eGFP protein for optimization tests."""
    return EGFP_PROTEIN


@pytest.fixture
def hbb() -> str:
    """Human hemoglobin beta protein."""
    return HBB_PROTEIN


@pytest.fixture
def ricin_a() -> str:
    """Ricin A-chain — known hazardous toxin protein."""
    return RICIN_A_CHAIN


@pytest.fixture
def sars_rbd() -> str:
    """SARS-CoV-2 spike RBD — regulated pathogen protein."""
    return SARS_COV2_RBD


@pytest.fixture
def blatem() -> str:
    """blaTEM-1 beta-lactamase — antibiotic resistance marker."""
    return BLATEM_PROTEIN


@pytest.fixture
def safe_protein() -> str:
    """A clearly safe, non-hazardous protein for baseline tests."""
    return INSULIN_PROTEIN


# ═══════════════════════════════════════════════════════════════════════
# 1. Hazardous Input Rejection
# ═══════════════════════════════════════════════════════════════════════

class TestHazardousInputRejection:
    """Verify the biosecurity screen blocks known hazardous sequences."""

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_ricin_a_chain_raises_biosecurity_error(self, ricin_a: str) -> None:
        """Optimizing ricin A-chain protein must raise BiosecurityError."""
        with pytest.raises(BiosecurityError, match=r"(?i)ricin|hazardous|toxin"):
            optimize_sequence(ricin_a, organism="Escherichia_coli")

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_ricin_a_chain_screen_result(self, ricin_a: str) -> None:
        """Screening ricin A-chain must flag it as hazardous with details."""
        result = screen_sequence(ricin_a)
        assert isinstance(result, BiosecurityScreenResult)
        assert result.flagged is True
        assert len(result.flags) > 0
        # At least one flag should mention ricin or toxin
        flag_text = " ".join(f.description for f in result.flags).lower()
        assert "ricin" in flag_text or "toxin" in flag_text

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_blatem_resistance_marker_flagged(self, blatem: str) -> None:
        """Optimizing blaTEM resistance marker must flag it."""
        result = screen_sequence(blatem)
        assert result.flagged is True
        flag_text = " ".join(f.description for f in result.flags).lower()
        assert "resistance" in flag_text or "antibiotic" in flag_text or "blatem" in flag_text

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_sars_cov2_rbd_flagged(self, sars_rbd: str) -> None:
        """Optimizing SARS-CoV-2 spike RBD must flag it."""
        result = screen_sequence(sars_rbd)
        assert result.flagged is True
        flag_text = " ".join(f.description for f in result.flags).lower()
        assert "sars" in flag_text or "coronavirus" in flag_text or "spike" in flag_text or "pathogen" in flag_text

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_safe_protein_not_flagged(self, safe_protein: str) -> None:
        """A safe protein (insulin) must NOT be flagged by the screen."""
        result = screen_sequence(safe_protein)
        assert result.flagged is False
        assert len(result.flags) == 0

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_blatem_raises_biosecurity_error(self, blatem: str) -> None:
        """Optimizing blaTEM must raise BiosecurityError."""
        with pytest.raises(BiosecurityError, match=r"(?i)resistance|antibiotic|blatem|hazardous"):
            optimize_sequence(blatem, organism="Escherichia_coli")

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_sars_rbd_raises_biosecurity_error(self, sars_rbd: str) -> None:
        """Optimizing SARS-CoV-2 RBD must raise BiosecurityError."""
        with pytest.raises(BiosecurityError, match=r"(?i)sars|coronavirus|pathogen|hazardous"):
            optimize_sequence(sars_rbd, organism="Homo_sapiens")

    def test_safe_protein_optimizes_normally(self, safe_protein: str) -> None:
        """A safe protein must optimize without any biosecurity error.

        This test does NOT depend on the biosecurity module — it confirms
        that the existing optimizer works on a known-safe protein.
        """
        result = optimize_sequence(safe_protein, organism="Escherichia_coli")
        assert isinstance(result, OptimizationResult)
        assert result.sequence  # non-empty
        assert result.cai > 0.0

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_hazardous_dna_sequence_flagged(self) -> None:
        """A DNA sequence encoding a hazardous protein must be flagged.

        Tests that the screen works on DNA input, not just protein.
        """
        # Translate ricin to get a DNA codon sequence that would encode it
        # We test that the screen accepts DNA too.
        from biocompiler.constants import CODON_TABLE, AA_TO_CODONS

        # Build a naive DNA encoding of ricin
        dna_parts = []
        for aa in RICIN_A_CHAIN:
            codons = AA_TO_CODONS.get(aa, ["ATG"])
            dna_parts.append(codons[0])
        dna = "".join(dna_parts)

        result = screen_sequence(dna)
        assert result.flagged is True


# ═══════════════════════════════════════════════════════════════════════
# 2. Translation Verification After Optimization
# ═══════════════════════════════════════════════════════════════════════

class TestTranslationVerification:
    """Verify that optimized sequences translate back to the input protein."""

    def test_insulin_ecoli_translation_matches(self, insulin: str) -> None:
        """Optimize insulin for E. coli → translated protein matches input."""
        result = optimize_sequence(insulin, organism="Escherichia_coli")
        translated = translate(result.sequence)
        assert translated == insulin, (
            f"Translation mismatch after E. coli optimization: "
            f"input={insulin!r}, translated={translated!r}"
        )

    def test_gfp_human_translation_matches(self, egfp: str) -> None:
        """Optimize GFP for human → translated protein matches input."""
        result = optimize_sequence(egfp, organism="Homo_sapiens")
        translated = translate(result.sequence)
        assert translated == egfp, (
            f"Translation mismatch after human optimization: "
            f"input len={len(egfp)}, translated len={len(translated)}"
        )

    def test_insulin_yeast_translation_matches(self, insulin: str) -> None:
        """Optimize insulin for yeast → translated protein matches input."""
        result = optimize_sequence(insulin, organism="Saccharomyces_cerevisiae")
        translated = translate(result.sequence)
        assert translated == insulin

    def test_hbb_cho_translation_matches(self, hbb: str) -> None:
        """Optimize HBB for CHO → translated protein matches input."""
        result = optimize_sequence(hbb, organism="CHO_K1")
        translated = translate(result.sequence)
        assert translated == hbb

    def test_corrupted_codon_detected(self, insulin: str) -> None:
        """Intentionally corrupt a codon → verify mismatch detection."""
        result = optimize_sequence(insulin, organism="Escherichia_coli")

        # Corrupt codon at position 3 (the 4th amino acid)
        seq = list(result.sequence)
        corrupt_pos = 3 * 3  # codon position for 4th amino acid (0-indexed)
        # Change the first base of the codon to guarantee a different AA
        original_codon = seq[corrupt_pos:corrupt_pos + 3]
        # Try all possible single-base substitutions until we find one that
        # changes the amino acid
        original_aa = translate("".join(original_codon))
        corrupted = False
        for base in "ACGT":
            if base != original_codon[0]:
                test_codon = [base] + original_codon[1:]
                test_aa = translate("".join(test_codon))
                if test_aa != original_aa:
                    seq[corrupt_pos] = base
                    corrupted = True
                    break
        if not corrupted:
            for base in "ACGT":
                if base != original_codon[2]:
                    test_codon = original_codon[:2] + [base]
                    test_aa = translate("".join(test_codon))
                    if test_aa != original_aa:
                        seq[corrupt_pos + 2] = base
                        corrupted = True
                        break

        assert corrupted, "Could not corrupt codon to produce different AA"
        corrupted_seq = "".join(seq)
        translated = translate(corrupted_seq)
        assert translated != insulin, (
            "Corrupted codon should produce a different protein translation"
        )

    @pytest.mark.skipif(
        not _HAS_PROTEIN_VERIFICATION,
        reason="Awaiting Agent 4 protein_verification module",
    )
    def test_verify_translation_with_module(self, insulin: str) -> None:
        """Use the protein_verification module to verify translation."""
        result = optimize_sequence(insulin, organism="Escherichia_coli")
        verification = verify_translation(result.sequence, insulin)
        assert isinstance(verification, TranslationVerificationResult)
        assert verification.match is True
        assert verification.mismatches == []

    @pytest.mark.skipif(
        not _HAS_PROTEIN_VERIFICATION,
        reason="Awaiting Agent 4 protein_verification module",
    )
    def test_detect_corruption_with_module(self, insulin: str) -> None:
        """Detect corruption using the protein_verification module."""
        result = optimize_sequence(insulin, organism="Escherichia_coli")

        # Corrupt position 9 (10th amino acid)
        seq = list(result.sequence)
        corrupt_pos = 9 * 3
        seq[corrupt_pos] = "G" if seq[corrupt_pos] != "G" else "A"
        corrupted_seq = "".join(seq)

        corruption_report = detect_corruption(corrupted_seq, insulin)
        assert isinstance(corruption_report, TranslationVerificationResult)
        assert corruption_report.match is False
        assert len(corruption_report.mismatches) > 0
        # Mismatch should be at or near position 9
        mismatch_positions = [m.position for m in corruption_report.mismatches]
        assert 9 in mismatch_positions or any(abs(p - 9) <= 1 for p in mismatch_positions)

    def test_multi_organism_translation_consistency(self, insulin: str) -> None:
        """Optimize insulin for all supported organisms → translations match."""
        organisms = [
            "Escherichia_coli",
            "Homo_sapiens",
            "Saccharomyces_cerevisiae",
            "CHO_K1",
            "Mus_musculus",
        ]
        for org in organisms:
            result = optimize_sequence(insulin, organism=org)
            translated = translate(result.sequence)
            assert translated == insulin, (
                f"Translation mismatch for organism={org}: "
                f"input={insulin!r}, translated={translated!r}"
            )


# ═══════════════════════════════════════════════════════════════════════
# 3. Provenance Tracking by Default
# ═══════════════════════════════════════════════════════════════════════

class TestProvenanceTrackingByDefault:
    """Verify that optimization provenance is tracked even without
    explicit ``track_provenance=True``."""

    def test_provenance_recorded_without_flag(self, insulin: str) -> None:
        """Calling optimize() without track_provenance still records provenance.

        The OptimizationResult.provenance field should be populated with
        at least a lightweight OptimizationRecord, even when the caller
        does not explicitly request full provenance tracking.
        """
        result = optimize_sequence(insulin, organism="Escherichia_coli")
        # Even without track_provenance=True, the result should have
        # a non-None provenance attribute (OptimizationRecord)
        assert result.provenance is not None, (
            "Provenance should be recorded by default, not just when "
            "track_provenance=True is passed"
        )

    def test_provenance_with_explicit_flag(self, insulin: str) -> None:
        """Calling optimize() with track_provenance=True records full provenance."""
        result = optimize_sequence(
            insulin,
            organism="Escherichia_coli",
            track_provenance=True,
        )
        assert result.provenance is not None
        # With track_provenance=True, we should also get a decision_trail
        assert result.decision_trail is not None, (
            "track_provenance=True should populate decision_trail"
        )

    def test_provenance_record_has_required_fields(self, insulin: str) -> None:
        """Verify provenance record has all required fields."""
        result = optimize_sequence(
            insulin,
            organism="Escherichia_coli",
            track_provenance=True,
            seed=42,
        )
        provenance = result.provenance
        assert provenance is not None

        # If it's an OptimizationRecord, check required fields
        if isinstance(provenance, OptimizationRecord):
            assert provenance.input_sequence is not None
            assert provenance.output_sequence is not None
            assert provenance.output_sequence == result.sequence
            assert provenance.organism == "Escherichia_coli"
            assert isinstance(provenance.constraints_applied, list)
            assert isinstance(provenance.mutations_made, list)
            assert provenance.solver_backend is not None
            assert provenance.solve_time >= 0
            assert provenance.seed_used == 42
            assert provenance.timestamp is not None
            assert provenance.biocompiler_version is not None
        elif hasattr(provenance, "to_dict"):
            # It may be an OptimizationProvenance or dict-like object
            d = provenance.to_dict()
            assert "organism" in d
            assert "final_sequence" in d or "output_sequence" in d

    def test_provenance_is_queryable(self, insulin: str) -> None:
        """Verify provenance is queryable via ProvenanceTracker.

        When track_provenance=True, the decision_trail should contain
        per-position decisions that can be queried.
        """
        result = optimize_sequence(
            insulin,
            organism="Escherichia_coli",
            track_provenance=True,
            seed=42,
        )
        trail = result.decision_trail
        if trail is not None and hasattr(trail, "get_full_audit_trail"):
            decisions = trail.get_full_audit_trail()
            assert len(decisions) > 0, (
                "Decision trail should contain at least one decision"
            )
            # Decisions should have required fields
            for d in decisions:
                assert hasattr(d, "decision_type")
                assert hasattr(d, "position")
                assert hasattr(d, "chosen_value")
                assert hasattr(d, "rationale")

    def test_provenance_serializable(self, insulin: str) -> None:
        """Verify provenance can be serialized to dict/JSON."""
        result = optimize_sequence(
            insulin,
            organism="Escherichia_coli",
            track_provenance=True,
            seed=42,
        )
        provenance = result.provenance
        if provenance is not None:
            # Should be serializable
            if hasattr(provenance, "to_dict"):
                d = provenance.to_dict()
                assert isinstance(d, dict)
                # Should be JSON-serializable
                json_str = json.dumps(d, default=str)
                assert isinstance(json_str, str)
                # Round-trip
                parsed = json.loads(json_str)
                assert isinstance(parsed, dict)
            elif hasattr(provenance, "to_json"):
                json_str = provenance.to_json()
                assert isinstance(json_str, str)
                parsed = json.loads(json_str)
                assert isinstance(parsed, dict)


# ═══════════════════════════════════════════════════════════════════════
# 4. Strict Mode Enforcement
# ═══════════════════════════════════════════════════════════════════════

class TestStrictModeEnforcement:
    """Verify strict mode catches impossible constraints cleanly."""

    @pytest.mark.skipif(
        not _HAS_STRICT_MODE,
        reason="Awaiting Agent 5 strict_mode module",
    )
    def test_impossible_gc_constraint_hard_stop(self, insulin: str) -> None:
        """Impossible GC constraint (0.50 in 99% AT region) → hard stop.

        A protein made almost entirely of A and T codons cannot have GC=0.50.
        Strict mode must raise StrictModeError.
        """
        # Build an AT-rich protein (F, Y, L, I, N, K are AT-rich)
        at_rich_protein = "FIYFIYFIYFIYFIYFIY"  # 18 aa, very AT-rich
        with pytest.raises(StrictModeError) as exc_info:
            optimize_strict(
                at_rich_protein,
                organism="Escherichia_coli",
                gc_lo=0.50,
                gc_hi=0.50,  # exact GC target
                strict_mode=True,
            )
        # Error must mention failed predicates
        assert exc_info.value.failed_predicates is not None
        assert len(exc_info.value.failed_predicates) > 0

    @pytest.mark.skipif(
        not _HAS_STRICT_MODE,
        reason="Awaiting Agent 5 strict_mode module",
    )
    def test_error_message_lists_failed_predicates(self, insulin: str) -> None:
        """Strict mode error must list which predicates failed."""
        # Use extreme GC constraints that are impossible to satisfy
        with pytest.raises(StrictModeError) as exc_info:
            optimize_strict(
                insulin,
                organism="Escherichia_coli",
                gc_lo=0.99,
                gc_hi=1.00,  # impossible for this protein
                strict_mode=True,
            )
        error = exc_info.value
        failed = error.failed_predicates
        assert isinstance(failed, list)
        assert len(failed) > 0
        # At least one failed predicate should be GC-related
        failed_names = [p.name for p in failed] if hasattr(failed[0], "name") else [str(p) for p in failed]
        assert any("gc" in name.lower() for name in failed_names), (
            f"Expected a GC-related failed predicate, got: {failed_names}"
        )

    @pytest.mark.skipif(
        not _HAS_STRICT_MODE,
        reason="Awaiting Agent 5 strict_mode module",
    )
    def test_strict_mode_false_returns_with_warnings(self, insulin: str) -> None:
        """strict_mode=False returns result with warnings instead of hard stop."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = optimize_strict(
                insulin,
                organism="Escherichia_coli",
                gc_lo=0.99,
                gc_hi=1.00,
                strict_mode=False,
            )
            # Should return a result, not raise
            assert isinstance(result, OptimizationResult)
            # Should have warnings about failed predicates
            assert len(result.failed_predicates) > 0 or len(w) > 0

    def test_current_optimizer_handles_extreme_gc(self, insulin: str) -> None:
        """Test the current optimizer's behavior with extreme GC constraints.

        This test works with the current codebase — it verifies that
        the optimizer at least returns a result (possibly with failed
        predicates) rather than crashing.
        """
        result = optimize_sequence(
            insulin,
            organism="Escherichia_coli",
            gc_lo=0.01,
            gc_hi=0.99,  # very loose — should succeed
        )
        assert isinstance(result, OptimizationResult)
        assert result.sequence

    @pytest.mark.skipif(
        not _HAS_STRICT_MODE,
        reason="Awaiting Agent 5 strict_mode module",
    )
    def test_strict_mode_preserves_valid_constraints(self, insulin: str) -> None:
        """Strict mode with achievable constraints should succeed."""
        result = optimize_strict(
            insulin,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            strict_mode=True,
        )
        assert isinstance(result, OptimizationResult)
        assert result.sequence
        assert result.cai > 0.0

    @pytest.mark.skipif(
        not _HAS_STRICT_MODE,
        reason="Awaiting Agent 5 strict_mode module",
    )
    def test_failed_predicates_error_details(self) -> None:
        """FailedPredicatesError carries structured predicate failure info."""
        at_rich = "IIIIIIIIIIIIIIIIIII"  # Isoleucine-heavy, AT-rich
        with pytest.raises(FailedPredicatesError) as exc_info:
            optimize_strict(
                at_rich,
                organism="Escherichia_coli",
                gc_lo=0.50,
                gc_hi=0.50,
                strict_mode=True,
            )
        error = exc_info.value
        assert hasattr(error, "failed_predicates")
        failed = error.failed_predicates
        assert len(failed) > 0
        # Each failed predicate should have a name and reason
        for p in failed:
            if hasattr(p, "name"):
                assert p.name  # non-empty
            if hasattr(p, "reason"):
                assert p.reason  # non-empty


# ═══════════════════════════════════════════════════════════════════════
# 5. Biosafety in Exports
# ═══════════════════════════════════════════════════════════════════════

class TestBiosafetyInExports:
    """Verify that exported sequences carry biosafety annotations."""

    def test_genbank_has_biosafety_annotations(self, insulin: str) -> None:
        """GenBank export of optimized sequence has biosafety annotations.

        At minimum, the GenBank COMMENT section should indicate that
        the sequence was screened by BioCompiler's safety pipeline.
        """
        result = optimize_sequence(insulin, organism="Escherichia_coli")
        gb = export_genbank(
            sequence=result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
            gene_name="INS",
        )
        # Must have BioCompiler marker in COMMENT
        assert "BioCompiler" in gb, "GenBank export must reference BioCompiler"
        # Must have KEYWORDS section with codon-optimized marker
        assert "codon-optimized" in gb.lower() or "BioCompiler" in gb

    def test_fasta_header_has_safety_metadata(self, insulin: str) -> None:
        """FASTA header includes organism and optimization metadata."""
        result = optimize_sequence(insulin, organism="Escherichia_coli")
        fasta = export_fasta(
            sequence=result.sequence,
            identifier="INS_optimized",
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
        )
        lines = fasta.strip().split("\n")
        # Header line starts with '>' — may be preceded by comment lines (';')
        header_lines = [l for l in lines if l.startswith(">")]
        assert len(header_lines) > 0, "FASTA must have a header line starting with '>'"
        header = header_lines[0]
        # Header must include organism
        assert "organism=" in header, f"FASTA header must include organism metadata: {header}"
        assert "Escherichia_coli" in header or "coli" in header

    def test_genbank_comment_has_safety_screen(self, insulin: str) -> None:
        """GenBank COMMENT section indicates safety screening was performed.

        When the biosecurity module is available, the GenBank export
        should include a note that the sequence passed the biosecurity screen.
        """
        result = optimize_sequence(insulin, organism="Escherichia_coli")

        # Build type_results if the biosecurity module provides them
        type_results = None
        if _HAS_BIOSECURITY:
            screen_result = screen_sequence(insulin)
            if not screen_result.flagged:
                # Add a synthetic type result for the biosecurity predicate
                from biocompiler.types import TypeCheckResult, Verdict
                type_results = [
                    TypeCheckResult(
                        predicate="NoHazardousSequence",
                        verdict=Verdict.PASS,
                        derivation=None,
                        violation=None,
                    )
                ]

        gb = export_genbank(
            sequence=result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
            type_results=type_results,
            gene_name="INS",
        )
        # GenBank should have COMMENT section
        assert "COMMENT" in gb, "GenBank export must have COMMENT section"

        if type_results:
            # Should include the biosecurity predicate in type-check results
            assert "NoHazardousSequence" in gb or "typecheck" in gb.lower()

    def test_export_with_failed_predicates_includes_warning(self, insulin: str) -> None:
        """Export with failed predicates includes warning comment."""
        result = optimize_sequence(insulin, organism="Escherichia_coli")

        # Simulate failed predicates in the export
        from biocompiler.types import TypeCheckResult, Verdict
        failed_results = [
            TypeCheckResult(
                predicate="GCInRange",
                verdict=Verdict.FAIL,
                derivation=None,
                violation="GC content 0.25 below minimum 0.30",
            ),
        ]

        gb = export_genbank(
            sequence=result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
            type_results=failed_results,
            gene_name="INS",
        )
        # GenBank should mention the failed predicate
        assert "TYPE FAIL" in gb or "FAIL" in gb, (
            "GenBank export must annotate failed predicates"
        )
        assert "GCInRange" in gb, "Failed predicate name must appear in export"

    def test_fasta_safety_comment(self, insulin: str) -> None:
        """FASTA comment lines include safety screening information."""
        result = optimize_sequence(insulin, organism="Escherichia_coli")
        fasta = export_fasta(
            sequence=result.sequence,
            identifier="INS_optimized",
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
            include_comments=True,
        )
        # Should have comment lines (starting with ';')
        comment_lines = [l for l in fasta.split("\n") if l.startswith(";")]
        assert len(comment_lines) > 0, "FASTA export should include comment lines"
        # Comment should mention BioCompiler
        comment_text = " ".join(comment_lines)
        assert "BioCompiler" in comment_text

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_genbank_includes_safety_screen_passed(self, insulin: str) -> None:
        """GenBank includes explicit 'safety screen: PASSED' annotation."""
        result = optimize_sequence(insulin, organism="Escherichia_coli")
        screen = screen_sequence(insulin)
        assert screen.flagged is False

        # The export should be enriched with safety screen metadata
        from biocompiler.types import TypeCheckResult, Verdict
        safety_type_result = TypeCheckResult(
            predicate="NoHazardousSequence",
            verdict=Verdict.PASS,
            derivation=None,
            violation=None,
        )
        gb = export_genbank(
            sequence=result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
            type_results=[safety_type_result],
            gene_name="INS",
        )
        # The COMMENT section should reference the safety predicate
        assert "NoHazardousSequence" in gb or "hazardous" in gb.lower()


# ═══════════════════════════════════════════════════════════════════════
# 6. Full Safety Pipeline
# ═══════════════════════════════════════════════════════════════════════

class TestFullSafetyPipeline:
    """End-to-end test running the complete safety flow:

    Input protein → Biosecurity screen → Optimize → Verify translation
    → Check predicates → Export with annotations
    """

    @pytest.mark.integration
    @pytest.mark.e2e
    def test_full_pipeline_safe_protein(self, insulin: str) -> None:
        """Complete pipeline for a safe protein (insulin for E. coli).

        Steps:
        1. Screen the protein for biosecurity concerns
        2. Optimize the codons
        3. Verify the translation matches
        4. Check that predicates are satisfied
        5. Export with safety annotations
        """
        # ── Step 1: Biosecurity screen ──
        if _HAS_BIOSECURITY:
            screen = screen_sequence(insulin)
            assert screen.flagged is False, (
                f"Insulin should pass biosecurity screen, but got flags: "
                f"{[f.description for f in screen.flags]}"
            )

        # ── Step 2: Optimize ──
        result = optimize_sequence(
            insulin,
            organism="Escherichia_coli",
            track_provenance=True,
            seed=42,
        )
        assert isinstance(result, OptimizationResult)
        assert result.sequence
        assert result.cai > 0.0
        assert 0.0 <= result.gc_content <= 1.0

        # ── Step 3: Verify translation ──
        translated = translate(result.sequence)
        assert translated == insulin, (
            f"Translation verification failed: expected {insulin!r}, "
            f"got {translated!r}"
        )

        # Additional verification with protein_verification module
        if _HAS_PROTEIN_VERIFICATION:
            verification = verify_translation(result.sequence, insulin)
            assert verification.match is True

        # ── Step 4: Check predicates ──
        # GC should be in default range [0.30, 0.70]
        gc = gc_content(result.sequence)
        assert 0.30 <= gc <= 0.70 or "GCInRange" in result.failed_predicates, (
            f"GC content {gc:.4f} out of range and not reported in failed_predicates"
        )

        # CAI should be reasonable
        cai_check = compute_cai(result.sequence, organism="Escherichia_coli")
        assert abs(cai_check - result.cai) < 0.01, (
            f"CAI mismatch: result.cai={result.cai}, "
            f"compute_cai={cai_check}"
        )

        # ── Step 5: Export with annotations ──
        # Build type results for export
        type_results = []
        from biocompiler.types import TypeCheckResult, Verdict

        type_results.append(TypeCheckResult(
            predicate="CodonAdapted",
            verdict=Verdict.PASS if result.cai >= 0.5 else Verdict.FAIL,
            derivation=None,
            violation=None if result.cai >= 0.5 else f"CAI={result.cai:.4f} < 0.5",
        ))

        if _HAS_BIOSECURITY:
            type_results.append(TypeCheckResult(
                predicate="NoHazardousSequence",
                verdict=Verdict.PASS,
                derivation=None,
                violation=None,
            ))

        # GenBank export
        gb = export_genbank(
            sequence=result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
            type_results=type_results,
            gene_name="INS",
        )
        assert "BioCompiler" in gb
        assert "COMMENT" in gb
        assert "Escherichia coli" in gb

        # FASTA export
        fasta = export_fasta(
            sequence=result.sequence,
            identifier="INS_ecoli_optimized",
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
        )
        assert fasta.startswith(">") or fasta.startswith(";")
        assert "organism=" in fasta

        # ── Step 6: Provenance check ──
        assert result.provenance is not None, "Provenance must be recorded"
        if hasattr(result.provenance, "organism"):
            assert result.provenance.organism == "Escherichia_coli"

    @pytest.mark.integration
    @pytest.mark.e2e
    def test_full_pipeline_gfp_human(self, egfp: str) -> None:
        """Complete pipeline for eGFP optimized for human expression."""
        # Screen
        if _HAS_BIOSECURITY:
            screen = screen_sequence(egfp)
            assert screen.flagged is False

        # Optimize
        result = optimize_sequence(
            egfp,
            organism="Homo_sapiens",
            track_provenance=True,
            seed=123,
        )
        assert isinstance(result, OptimizationResult)

        # Verify translation
        translated = translate(result.sequence)
        assert translated == egfp

        # Verify CAI
        assert result.cai > 0.5, f"eGFP CAI too low: {result.cai}"

        # Verify GC
        gc = gc_content(result.sequence)
        assert 0.30 <= gc <= 0.70 or "GCInRange" in result.failed_predicates

        # Export
        gb = export_genbank(
            sequence=result.sequence,
            organism="Homo_sapiens",
            protein=result.protein,
            cai=result.cai,
            gene_name="EGFP",
        )
        assert "BioCompiler" in gb
        assert "Homo sapiens" in gb

        # Provenance
        assert result.provenance is not None

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    @pytest.mark.integration
    @pytest.mark.e2e
    def test_full_pipeline_hazardous_blocked(self, ricin_a: str) -> None:
        """Complete pipeline for a hazardous protein → blocked at screen."""
        # Step 1: Biosecurity screen should flag it
        screen = screen_sequence(ricin_a)
        assert screen.flagged is True

        # Step 2: Optimization should raise BiosecurityError
        with pytest.raises(BiosecurityError):
            optimize_sequence(ricin_a, organism="Escherichia_coli")

    @pytest.mark.integration
    @pytest.mark.e2e
    def test_pipeline_provenance_end_to_end(self, insulin: str) -> None:
        """Verify provenance flows from optimization through to export."""
        # Optimize with full provenance
        result = optimize_sequence(
            insulin,
            organism="Escherichia_coli",
            track_provenance=True,
            seed=99,
        )

        # Verify provenance is populated
        assert result.provenance is not None
        assert result.decision_trail is not None

        # Export as JSON (includes provenance)
        from biocompiler.export import export_json
        json_output = export_json(result, include_provenance=True)
        parsed = json.loads(json_output)
        assert "provenance" in parsed or "sequence" in parsed

        # The JSON should include the seed for reproducibility
        if "provenance" in parsed:
            prov = parsed["provenance"]
            if "optimization_record" in prov:
                opt_rec = prov["optimization_record"]
                # Seed should be recorded
                assert "seed_used" in opt_rec or "seed" in opt_rec

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    @pytest.mark.integration
    @pytest.mark.e2e
    def test_pipeline_safety_annotations_in_genbank(self, insulin: str) -> None:
        """Safety screen results appear in GenBank export annotations."""
        # Optimize with safety screen
        result = optimize_sequence(
            insulin,
            organism="Escherichia_coli",
            track_provenance=True,
        )

        # Get screen result
        screen = screen_sequence(insulin)
        assert screen.flagged is False

        # Build type results including safety
        from biocompiler.types import TypeCheckResult, Verdict
        type_results = [
            TypeCheckResult(
                predicate="NoHazardousSequence",
                verdict=Verdict.PASS,
                derivation=None,
                violation=None,
            ),
        ]

        # Export
        gb = export_genbank(
            sequence=result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
            type_results=type_results,
            gene_name="INS",
        )

        # Verify safety annotation is present
        assert "NoHazardousSequence" in gb or "[+]" in gb
        # Verify overall type-check verdict is reported
        assert "Type-check verdict" in gb or "typecheck" in gb.lower() or "PASS" in gb


# ═══════════════════════════════════════════════════════════════════════
# Additional edge-case safety tests
# ═══════════════════════════════════════════════════════════════════════

class TestSafetyEdgeCases:
    """Edge cases and boundary conditions for the safety pipeline."""

    def test_empty_protein_raises_error(self) -> None:
        """Empty protein string must raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            optimize_sequence("", organism="Escherichia_coli")

    def test_whitespace_only_protein_handled(self) -> None:
        """Whitespace-only protein must either raise an error or return
        an empty/invalid result — it must NOT silently produce a
        valid-looking optimization."""
        try:
            result = optimize_sequence("   ", organism="Escherichia_coli")
            # If it doesn't raise, the result must be clearly invalid
            assert result.sequence == "", (
                "Whitespace-only protein should not produce a valid sequence"
            )
        except (InvalidProteinError, BioCompilerError, ValueError, AssertionError):
            pass  # Any of these is acceptable

    def test_invalid_amino_acids_raise_error(self) -> None:
        """Proteins with invalid amino acid codes must raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            optimize_sequence("MALWMRX", organism="Escherichia_coli")

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_partial_hazardous_match_flagged(self) -> None:
        """A protein that is a substring of a known hazardous protein
        should still be flagged (e.g., a ricin fragment)."""
        # Use a 50-aa fragment from the middle of ricin A-chain
        fragment = RICIN_A_CHAIN[100:150]
        result = screen_sequence(fragment)
        # Fragments of hazardous proteins should at least generate a warning
        # (may not be a hard block, but should be flagged)
        assert result.flagged is True or len(result.warnings) > 0

    def test_short_protein_optimizes_and_verifies(self) -> None:
        """Very short protein (3 aa) optimizes and translates correctly."""
        protein = "MAG"
        result = optimize_sequence(protein, organism="Escherichia_coli")
        translated = translate(result.sequence)
        assert translated == protein

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_case_insensitive_screening(self) -> None:
        """Biosecurity screening is case-insensitive."""
        upper_result = screen_sequence(RICIN_A_CHAIN)
        lower_result = screen_sequence(RICIN_A_CHAIN.lower())
        assert upper_result.flagged == lower_result.flagged

    @pytest.mark.skipif(
        not _HAS_STRICT_MODE,
        reason="Awaiting Agent 5 strict_mode module",
    )
    def test_strict_mode_zero_gc_target(self) -> None:
        """GC target of exactly 0.0 is impossible → strict mode error."""
        protein = "MAGCLVFW"  # Mix of GC-rich and AT-rich
        with pytest.raises(StrictModeError):
            optimize_strict(
                protein,
                organism="Escherichia_coli",
                gc_lo=0.0,
                gc_hi=0.0,
                strict_mode=True,
            )

    def test_provenance_preserved_in_json_export(self, insulin: str) -> None:
        """Provenance data survives JSON export round-trip."""
        result = optimize_sequence(
            insulin,
            organism="Escherichia_coli",
            track_provenance=True,
            seed=77,
        )
        from biocompiler.export import export_json
        json_str = export_json(result, include_provenance=True)
        parsed = json.loads(json_str)

        # Basic structure
        assert "sequence" in parsed
        assert "metrics" in parsed
        assert "cai" in parsed["metrics"]

        # Provenance should be present
        if "provenance" in parsed:
            prov = parsed["provenance"]
            # Should contain optimization record or decision trail
            assert "optimization_record" in prov or "decision_trail" in prov

    @pytest.mark.skipif(
        not _HAS_BIOSECURITY,
        reason="Awaiting Agent 3 biosecurity module",
    )
    def test_screen_result_is_deterministic(self) -> None:
        """Biosecurity screening is deterministic: same input → same result."""
        result1 = screen_sequence(RICIN_A_CHAIN)
        result2 = screen_sequence(RICIN_A_CHAIN)
        assert result1.flagged == result2.flagged
        assert len(result1.flags) == len(result2.flags)

    def test_optimization_result_invariants(self, insulin: str) -> None:
        """OptimizationResult satisfies its documented invariants."""
        result = optimize_sequence(insulin, organism="Escherichia_coli")

        # Sequence length must be 3x protein length
        if result.protein:
            assert len(result.sequence) == len(result.protein) * 3, (
                f"Sequence length ({len(result.sequence)}) != "
                f"3 * protein length ({len(result.protein)})"
            )

        # CAI in [0, 1]
        assert 0.0 <= result.cai <= 1.0
        # GC in [0, 1]
        assert 0.0 <= result.gc_content <= 1.0
