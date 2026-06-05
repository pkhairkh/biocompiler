"""
Agent 50: Comprehensive Biosecurity Screening Tests
=====================================================

Expanded test suite for biocompiler.biosecurity covering:
- Biosecurity screening runs before every optimization
- Known toxins are blocked (ricin, botulinum)
- Antibiotic resistance markers are flagged
- Hard-stop mode works correctly
- Risk levels are appropriate (GFP should be "none" or "low")
- Biosecurity annotations in exports
- Valid research sequences are not broken by the screen
"""

import warnings

import pytest

from biocompiler.biosecurity import (
    BiosecurityReport,
    HazardMatch,
    HAZARD_SIGNATURE_COUNT,
    screen_hazardous_sequence,
    check_biosecurity_before_optimize,
    sig_risk_for_match,
    _build_recommendations,
    _max_risk,
    _HAZARD_SIGNATURES,
    _PROTEIN_SIGNATURES,
    _DNA_SIGNATURES,
)
from biocompiler.exceptions import BiosecurityError, BioCompilerError
from biocompiler.optimization import optimize_sequence, OptimizationResult
from biocompiler.export import export_genbank, export_fasta


# ─── Standard protein fixtures ─────────────────────────────────────

# Safe / benign proteins
INSULIN_PROTEIN = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
EGFP_FULL = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
    "VTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVN"
    "RIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQ"
    "QNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)
HBB_PROTEIN = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
    "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK"
    "EFTPPVQAAYQKVVAGVANALAHKYH"
)

# Hazardous proteins — select agent toxins
RICIN_A_MOTIF = "MISRDNIRVGLPIISTNKYEDKQL"  # Contains NIRVGLPIIS
BOTULINUM_MOTIF = "AAAHEHETQSNLRDLAAAA"  # Contains HETQSNLRDL
SHIGA_MOTIF = "AAANVRVNIPFSRAAAA"  # Contains NVRVNIPFSR
ANTHRAX_LF_MOTIF = "AAAHETHFGVVSYAAAA"  # Contains HETHFGVVSY
DIPHTHERIA_MOTIF = "AAAGVADVIQEINAAA"  # Contains GVADVIQEIN

# Antibiotic resistance markers
BLATEM_PROTEIN_MOTIF = "AAAHPETLALKFGAAA"  # Contains HPETLALKFG
NPTII_MOTIF = "AAARPMTIHGSGSAAA"  # Contains RPMTIHGSGS
NDM1_MOTIF = "MKHHHHHDGHFGGPVKVTLEAP"  # Contains HHHHDGHFGGP
VANA_MOTIF = "AAAHGLSSAVPGLAAA"  # Contains HGLSSAVPGL

# Viral surface proteins
INFLUENZA_HA_MOTIF = "AAAGLFGAIAGFIAAAA"  # Contains GLFGAIAGFI
SARS2_RBD_MOTIF = "AAAVYYHKNNKSWAAAA"  # Contains VYYHKNNKSW
HIV_V3_MOTIF = "AAAGPGRAFYTIGAAAA"  # Contains GPGRAFYTIG

# Antibiotic resistance DNA
BLATEM_DNA = "ATGAGTATTCAACATTTCCGTG"
NDM1_DNA = "ATGGAATTGCCCAATATTATG"
NPTII_DNA = "ATGATTGAACAAGATGGATTG"


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 1: Screening Runs Before Optimization
# ═══════════════════════════════════════════════════════════════════


class TestScreeningBeforeOptimization:
    """Verify that biosecurity screening is correctly invoked via the
    check_biosecurity_before_optimize gate and blocks hazardous inputs.

    The biosecurity screen is a gate that should be called before
    optimization.  It raises BiosecurityError for critical/high risk
    sequences and returns a BiosecurityReport for safe ones.
    """

    def test_safe_protein_passes_optimization(self):
        """Safe protein (insulin) passes the biosecurity gate and optimization."""
        # Step 1: Pass the biosecurity gate
        report = check_biosecurity_before_optimize(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        assert isinstance(report, BiosecurityReport)
        assert report.risk_level == "none"
        assert not report.is_hazardous

        # Step 2: Optimization can proceed
        result = optimize_sequence(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        assert isinstance(result, OptimizationResult)
        assert result.sequence
        assert result.cai > 0.0

    def test_hazardous_protein_blocked_at_gate(self):
        """Hazardous protein (ricin motif) is blocked at the biosecurity gate."""
        with pytest.raises(BiosecurityError) as exc_info:
            check_biosecurity_before_optimize(
                RICIN_A_MOTIF, organism="Escherichia_coli",
            )
        assert exc_info.value.report.risk_level in ("high", "critical")

    def test_botulinum_blocked_at_gate(self):
        """Botulinum toxin motif is blocked at the biosecurity gate."""
        with pytest.raises(BiosecurityError):
            check_biosecurity_before_optimize(
                BOTULINUM_MOTIF, organism="Escherichia_coli",
            )

    def test_check_biosecurity_returns_report_for_safe(self):
        """check_biosecurity_before_optimize returns a report for safe proteins."""
        report = check_biosecurity_before_optimize(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        assert isinstance(report, BiosecurityReport)
        assert report.risk_level == "none"
        assert not report.is_hazardous

    def test_check_biosecurity_raises_for_critical(self):
        """check_biosecurity_before_optimize raises for critical hazards."""
        with pytest.raises(BiosecurityError) as exc_info:
            check_biosecurity_before_optimize(
                RICIN_A_MOTIF, organism="Escherichia_coli",
            )
        assert exc_info.value.report.risk_level == "critical"

    def test_screen_function_always_returns_report(self):
        """screen_hazardous_sequence always returns a BiosecurityReport."""
        # Safe protein
        report = screen_hazardous_sequence(INSULIN_PROTEIN)
        assert isinstance(report, BiosecurityReport)

        # Hazardous protein
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        assert isinstance(report, BiosecurityReport)

        # Empty protein
        report = screen_hazardous_sequence("")
        assert isinstance(report, BiosecurityReport)

    def test_safe_protein_screen_then_optimize_workflow(self):
        """Complete workflow: screen → optimize → verify for safe protein.

        This demonstrates the recommended pattern: always call the
        biosecurity gate before optimizing.
        """
        protein = EGFP_FULL
        organism = "Homo_sapiens"

        # Step 1: Screen
        report = check_biosecurity_before_optimize(protein, organism=organism)
        assert report.risk_level in ("none", "low")

        # Step 2: Optimize (only reached if screen passed)
        result = optimize_sequence(protein, organism=organism)
        assert isinstance(result, OptimizationResult)
        assert result.sequence

    def test_hazardous_protein_workflow_blocked(self):
        """Complete workflow: screen blocks hazardous protein, optimization
        is never reached."""
        protein = RICIN_A_MOTIF
        organism = "Escherichia_coli"

        # Step 1: Screen raises error
        with pytest.raises(BiosecurityError):
            check_biosecurity_before_optimize(protein, organism=organism)

        # Step 2: Optimization should NOT be called — but even if
        # someone bypasses the gate, the screen catches it independently
        independent_report = screen_hazardous_sequence(protein)
        assert independent_report.risk_level == "critical"


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 2: Known Toxins Blocked
# ═══════════════════════════════════════════════════════════════════


class TestKnownToxinsBlocked:
    """Verify that known select agent toxins are detected and blocked."""

    def test_ricin_a_chain_catalytic_detected(self):
        """Ricin A-chain catalytic motif (NIRVGLPIIS) is detected."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        assert report.is_hazardous
        assert report.risk_level == "critical"
        assert "select_agent" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "ricin_A_chain_catalytic" in names

    def test_ricin_b_chain_lectin_detected(self):
        """Ricin B-chain lectin motif (QNRWIIRYVQ) is detected."""
        protein = "AAAQNRWIIRYVQAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "ricin_B_chain_lectin" in names

    def test_botulinum_zinc_protease_detected(self):
        """Botulinum zinc protease motif is detected."""
        report = screen_hazardous_sequence(BOTULINUM_MOTIF)
        assert report.is_hazardous
        assert "select_agent" in report.flagged_categories
        assert any("botulinum" in m.name for m in report.matches)

    def test_shiga_toxin_detected(self):
        """Shiga toxin A subunit motif is detected."""
        report = screen_hazardous_sequence(SHIGA_MOTIF)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "shiga_toxin_A_subunit" in names

    def test_anthrax_lf_detected(self):
        """Anthrax lethal factor motif is detected."""
        report = screen_hazardous_sequence(ANTHRAX_LF_MOTIF)
        assert report.is_hazardous
        assert any("anthrax_LF" in m.name for m in report.matches)

    def test_diphtheria_toxin_detected(self):
        """Diphtheria toxin ADR motif is detected."""
        report = screen_hazardous_sequence(DIPHTHERIA_MOTIF)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "diphtheria_toxin_ADR" in names

    def test_abrin_detected(self):
        """Abrin A-chain motif is detected."""
        protein = "AAANVRVGIPISRAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "abrin_A_chain" in names

    def test_tetanus_toxin_detected(self):
        """Tetanus toxin zinc protease motif is detected."""
        protein = "AAAHEIKSNIASKAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "tetanus_toxin_zinc" in names

    def test_cholera_toxin_detected(self):
        """Cholera toxin motif is detected."""
        protein = "AAARYVHHVSGQNAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "cholera_toxin_A1" in names

    def test_seb_superantigen_detected(self):
        """Staphylococcal enterotoxin B (SEB) is detected."""
        protein = "AAAVVPDLKDKSKAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "SEB_superantigen" in names

    def test_all_toxins_are_critical_or_high_risk(self):
        """All select agent toxin matches should be critical or high risk."""
        toxin_proteins = [
            RICIN_A_MOTIF,
            BOTULINUM_MOTIF,
            SHIGA_MOTIF,
            ANTHRAX_LF_MOTIF,
        ]
        for protein in toxin_proteins:
            report = screen_hazardous_sequence(protein)
            assert report.risk_level in ("high", "critical"), (
                f"Toxin protein should have high/critical risk, "
                f"got {report.risk_level}"
            )


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 3: Antibiotic Resistance Markers Flagged
# ═══════════════════════════════════════════════════════════════════


class TestAntibioticResistanceFlagged:
    """Verify antibiotic resistance markers are flagged correctly."""

    def test_blatem_protein_flagged(self):
        """blaTEM protein motif is flagged."""
        report = screen_hazardous_sequence(BLATEM_PROTEIN_MOTIF)
        assert report.is_hazardous
        assert "antibiotic_resistance" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "blaTEM_protein" in names

    def test_blatem_dna_flagged(self):
        """blaTEM DNA motif is flagged."""
        report = screen_hazardous_sequence("", BLATEM_DNA)
        assert report.is_hazardous
        assert "antibiotic_resistance" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "blaTEM_dna" in names

    def test_nptII_protein_flagged(self):
        """nptII protein motif is flagged."""
        report = screen_hazardous_sequence(NPTII_MOTIF)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "nptII_protein" in names

    def test_nptII_dna_flagged(self):
        """nptII DNA motif is flagged."""
        report = screen_hazardous_sequence("", NPTII_DNA)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "nptII_dna" in names

    def test_ndm1_protein_flagged_critical(self):
        """NDM-1 protein motif is flagged as critical."""
        report = screen_hazardous_sequence(NDM1_MOTIF)
        assert report.is_hazardous
        assert report.risk_level == "critical"
        names = [m.name for m in report.matches]
        assert "ndm1_protein" in names

    def test_ndm1_dna_flagged(self):
        """NDM-1 DNA motif is flagged."""
        report = screen_hazardous_sequence("", NDM1_DNA)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "ndm1_dna" in names

    def test_vana_protein_flagged(self):
        """vanA protein motif is flagged."""
        report = screen_hazardous_sequence(VANA_MOTIF)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "vanA_protein" in names

    def test_cat_protein_flagged(self):
        """cat (chloramphenicol acetyltransferase) protein is flagged."""
        protein = "AAAFHRGVCTNKAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "cat_protein" in names

    def test_both_protein_and_dna_detected_simultaneously(self):
        """Both protein and DNA resistance markers can be detected at once."""
        report = screen_hazardous_sequence(
            BLATEM_PROTEIN_MOTIF, BLATEM_DNA,
        )
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "blaTEM_protein" in names
        assert "blaTEM_dna" in names

    def test_resistance_category_in_recommendations(self):
        """Recommendations mention antibiotic resistance when detected."""
        report = screen_hazardous_sequence(BLATEM_PROTEIN_MOTIF)
        rec_text = " ".join(report.recommendations)
        assert "resistance" in rec_text.lower() or "antibiotic" in rec_text.lower()


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 4: Hard-Stop Mode
# ═══════════════════════════════════════════════════════════════════


class TestHardStopMode:
    """Verify that the hard-stop mode works correctly for high/critical risks."""

    def test_critical_raises_biosecurity_error(self):
        """Critical risk raises BiosecurityError (hard stop)."""
        with pytest.raises(BiosecurityError) as exc_info:
            check_biosecurity_before_optimize(RICIN_A_MOTIF)
        assert exc_info.value.report.risk_level == "critical"

    def test_high_raises_biosecurity_error(self):
        """High risk raises BiosecurityError (hard stop)."""
        # Influenza HA fusion peptide is high risk
        with pytest.raises(BiosecurityError):
            check_biosecurity_before_optimize(INFLUENZA_HA_MOTIF)

    def test_error_contains_full_report(self):
        """BiosecurityError carries the complete BiosecurityReport."""
        with pytest.raises(BiosecurityError) as exc_info:
            check_biosecurity_before_optimize(RICIN_A_MOTIF)
        error = exc_info.value
        assert hasattr(error, "report")
        assert isinstance(error.report, BiosecurityReport)
        assert error.report.is_hazardous is True
        assert len(error.report.matches) > 0
        assert error.report.risk_level == "critical"

    def test_error_is_biocompiler_error_subclass(self):
        """BiosecurityError is a subclass of BioCompilerError."""
        with pytest.raises(BioCompilerError):
            check_biosecurity_before_optimize(RICIN_A_MOTIF)

    def test_error_is_exception_subclass(self):
        """BiosecurityError is a subclass of Exception."""
        with pytest.raises(Exception):
            check_biosecurity_before_optimize(RICIN_A_MOTIF)

    def test_error_string_contains_risk_level(self):
        """Error string mentions the risk level."""
        with pytest.raises(BiosecurityError) as exc_info:
            check_biosecurity_before_optimize(RICIN_A_MOTIF)
        error_str = str(exc_info.value)
        assert "critical" in error_str

    def test_error_string_mentions_blocked(self):
        """Error string mentions 'blocked' or similar."""
        with pytest.raises(BiosecurityError) as exc_info:
            check_biosecurity_before_optimize(RICIN_A_MOTIF)
        error_str = str(exc_info.value).lower()
        assert "block" in error_str or "biosecurity" in error_str

    def test_medium_does_not_raise_error(self):
        """Medium risk does NOT raise BiosecurityError — only warns."""
        protein = NPTII_MOTIF
        report = screen_hazardous_sequence(protein)
        if report.risk_level == "medium":
            result = check_biosecurity_before_optimize(protein)
            assert isinstance(result, BiosecurityReport)

    def test_low_does_not_raise_error(self):
        """Low risk does NOT raise BiosecurityError."""
        # VEGF heparin binding — low risk
        protein = "AAAAPMAVPPPKKAAA"
        report = screen_hazardous_sequence(protein)
        if report.risk_level == "low":
            result = check_biosecurity_before_optimize(protein)
            assert isinstance(result, BiosecurityReport)

    def test_none_does_not_raise_error(self):
        """No risk does NOT raise BiosecurityError."""
        result = check_biosecurity_before_optimize(INSULIN_PROTEIN)
        assert isinstance(result, BiosecurityReport)
        assert result.risk_level == "none"


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 5: Risk Levels Appropriate
# ═══════════════════════════════════════════════════════════════════


class TestRiskLevelsAppropriate:
    """Verify that risk levels are appropriate for different sequence types."""

    def test_gfp_is_none_or_low_risk(self):
        """GFP (standard reporter protein) should be "none" or "low" risk."""
        report = screen_hazardous_sequence(EGFP_FULL)
        assert report.risk_level in ("none", "low"), (
            f"GFP should have 'none' or 'low' risk, got '{report.risk_level}'"
        )
        assert not report.is_hazardous or report.risk_level == "low"

    def test_insulin_is_none_risk(self):
        """Insulin should have "none" risk."""
        report = screen_hazardous_sequence(INSULIN_PROTEIN)
        assert report.risk_level == "none"
        assert not report.is_hazardous

    def test_hbb_is_none_or_low_risk(self):
        """Human hemoglobin beta should be "none" or "low" risk."""
        report = screen_hazardous_sequence(HBB_PROTEIN)
        assert report.risk_level in ("none", "low"), (
            f"HBB should have 'none' or 'low' risk, got '{report.risk_level}'"
        )

    def test_ricin_is_critical(self):
        """Ricin A-chain should be "critical" risk."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        assert report.risk_level == "critical"

    def test_botulinum_is_critical(self):
        """Botulinum should be "critical" risk."""
        report = screen_hazardous_sequence(BOTULINUM_MOTIF)
        assert report.risk_level == "critical"

    def test_ndm1_is_critical(self):
        """NDM-1 should be "critical" risk."""
        report = screen_hazardous_sequence(NDM1_MOTIF)
        assert report.risk_level == "critical"

    def test_blatem_protein_risk_level(self):
        """blaTEM protein should be "high" risk."""
        report = screen_hazardous_sequence(BLATEM_PROTEIN_MOTIF)
        assert report.risk_level in ("high", "critical")

    def test_viral_surface_at_least_medium(self):
        """Viral surface proteins should be at least "medium" risk."""
        report = screen_hazardous_sequence(INFLUENZA_HA_MOTIF)
        assert report.risk_level in ("medium", "high", "critical")

    def test_empty_protein_is_none_risk(self):
        """Empty protein should have "none" risk."""
        report = screen_hazardous_sequence("")
        assert report.risk_level == "none"
        assert not report.is_hazardous

    def test_short_random_peptide_is_none_or_low(self):
        """Short random peptide should typically be "none" or "low" risk."""
        protein = "MALWMRLLPL"
        report = screen_hazardous_sequence(protein)
        assert report.risk_level in ("none", "low"), (
            f"Short random peptide should be 'none' or 'low' risk, "
            f"got '{report.risk_level}'"
        )

    def test_highest_risk_wins_when_multiple_matches(self):
        """When multiple categories match, the highest risk level wins."""
        # Protein with both VEGF motif (low) and ricin motif (critical)
        protein = "APMAVPPPPKNIRVGLPIIS"
        report = screen_hazardous_sequence(protein)
        assert report.risk_level == "critical"


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 6: Biosecurity Annotations in Exports
# ═══════════════════════════════════════════════════════════════════


class TestBiosecurityAnnotationsInExports:
    """Verify that biosecurity information appears correctly in exports."""

    def test_genbank_has_biocompiler_annotations(self):
        """GenBank export includes BIOCOMPILER_ANNOTATIONS section."""
        result = optimize_sequence(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        gb = export_genbank(
            result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
        )
        assert "BIOCOMPILER_ANNOTATIONS:" in gb

    def test_genbank_has_biosecurity_screened(self):
        """GenBank export includes biosecurity_screened field."""
        result = optimize_sequence(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        gb = export_genbank(
            result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
        )
        assert "biosecurity_screened:" in gb

    def test_genbank_has_biosafety_level(self):
        """GenBank export includes biosafety_level field."""
        result = optimize_sequence(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        gb = export_genbank(
            result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
        )
        assert "biosafety_level:" in gb

    def test_genbank_ecoli_is_bsl1(self):
        """E. coli constructs are annotated as BSL-1."""
        result = optimize_sequence(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        gb = export_genbank(
            result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
        )
        assert "BSL-1" in gb

    def test_genbank_human_is_bsl2(self):
        """Human constructs are annotated as BSL-2."""
        result = optimize_sequence(
            INSULIN_PROTEIN, organism="Homo_sapiens",
            strict_mode=False,  # Avoid GC constraint failures for short seqs
        )
        gb = export_genbank(
            result.sequence,
            organism="Homo_sapiens",
            protein=result.protein,
            cai=result.cai,
        )
        assert "BSL-2" in gb

    def test_fasta_has_biosecurity_metadata(self):
        """FASTA header includes biosecurity metadata."""
        result = optimize_sequence(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        fasta = export_fasta(
            result.sequence,
            identifier="INS",
            organism="Escherichia_coli",
            cai=result.cai,
        )
        assert "biosecurity=" in fasta

    def test_genbank_has_biosecurity_notice(self):
        """GenBank export includes BIOSECURITY NOTICE."""
        result = optimize_sequence(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        gb = export_genbank(
            result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
        )
        assert "BIOSECURITY NOTICE:" in gb

    def test_genbank_notice_mentions_bsl1_for_ecoli(self):
        """Biosecurity notice for E. coli mentions BSL-1."""
        result = optimize_sequence(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        gb = export_genbank(
            result.sequence,
            organism="Escherichia_coli",
            protein=result.protein,
            cai=result.cai,
        )
        assert "BSL-1" in gb


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 7: Valid Research Sequences Not Broken
# ═══════════════════════════════════════════════════════════════════


class TestValidResearchNotBroken:
    """Verify that the biosecurity screen does not break valid research
    sequences — common lab proteins should pass without false positives."""

    def test_gfp_passes_screen(self):
        """GFP passes biosecurity screening without false positives."""
        report = screen_hazardous_sequence(EGFP_FULL)
        assert not report.is_hazardous or report.risk_level in ("none", "low")
        if report.risk_level == "none":
            assert len(report.matches) == 0

    def test_insulin_passes_screen(self):
        """Insulin passes biosecurity screening."""
        report = screen_hazardous_sequence(INSULIN_PROTEIN)
        assert report.risk_level == "none"
        assert not report.is_hazardous
        assert len(report.matches) == 0

    def test_hbb_passes_screen(self):
        """Human hemoglobin beta passes biosecurity screening."""
        report = screen_hazardous_sequence(HBB_PROTEIN)
        assert report.risk_level in ("none", "low")

    def test_common_lab_enzymes_pass(self):
        """Common lab enzyme sequences pass screening."""
        # T4 DNA ligase (partial)
        protein = "MRGKLLFFVLTSLLQSVQATKA"
        report = screen_hazardous_sequence(protein)
        assert report.risk_level in ("none", "low", "medium"), (
            f"Common lab enzyme should pass screening, got {report.risk_level}"
        )

    def test_optimized_sequence_passes_screen(self):
        """Optimized DNA for safe protein passes biosecurity screening."""
        result = optimize_sequence(
            INSULIN_PROTEIN, organism="Escherichia_coli",
        )
        # Screen the protein
        report = screen_hazardous_sequence(result.protein)
        assert report.risk_level == "none"

    def test_multiple_safe_proteins_all_pass(self):
        """Multiple well-known safe proteins all pass screening."""
        safe_proteins = [
            ("insulin", INSULIN_PROTEIN),
            ("GFP", EGFP_FULL),
            ("HBB", HBB_PROTEIN),
        ]
        for name, protein in safe_proteins:
            report = screen_hazardous_sequence(protein)
            assert report.risk_level in ("none", "low"), (
                f"{name} should be safe, got risk_level={report.risk_level}"
            )

    def test_case_insensitivity_safe_protein(self):
        """Case-insensitive screening works for safe proteins too."""
        report_lower = screen_hazardous_sequence(INSULIN_PROTEIN.lower())
        report_upper = screen_hazardous_sequence(INSULIN_PROTEIN.upper())
        assert report_lower.risk_level == report_upper.risk_level

    def test_short_peptide_not_false_positive(self):
        """Short peptides that aren't toxins should not trigger false positives."""
        protein = "MALWMR"
        report = screen_hazardous_sequence(protein)
        assert report.risk_level == "none"

    def test_research_antibody_sequence_passes(self):
        """A typical antibody variable region should not be flagged."""
        # Antibody heavy chain variable region (partial, generic)
        protein = "QVQLVQSGAEVKKPGASVKVSCKASGYTFT"
        report = screen_hazardous_sequence(protein)
        assert report.risk_level in ("none", "low", "medium"), (
            f"Antibody sequence should not be flagged as hazardous, "
            f"got {report.risk_level}"
        )


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 8: Hazard Signature Database Integrity
# ═══════════════════════════════════════════════════════════════════


class TestHazardDatabaseIntegrity:
    """Verify the hazard signature database is comprehensive and well-formed."""

    def test_at_least_50_signatures(self):
        """Hazard database has at least 50 signatures."""
        assert HAZARD_SIGNATURE_COUNT >= 50

    def test_protein_signatures_present(self):
        """Protein signatures exist in the database."""
        assert len(_PROTEIN_SIGNATURES) > 0

    def test_dna_signatures_present(self):
        """DNA signatures exist in the database."""
        assert len(_DNA_SIGNATURES) > 0

    def test_all_signatures_have_required_fields(self):
        """Every signature has all required fields."""
        required_fields = {"category", "name", "motif", "confidence", "risk", "type", "source"}
        for sig in _HAZARD_SIGNATURES:
            missing = required_fields - set(sig.keys())
            assert not missing, (
                f"Signature {sig.get('name', 'unknown')} missing fields: {missing}"
            )

    def test_all_categories_valid(self):
        """All signatures have valid categories."""
        valid_categories = {"select_agent", "viral_surface", "antibiotic_resistance", "oncogene"}
        for sig in _HAZARD_SIGNATURES:
            assert sig["category"] in valid_categories, (
                f"Invalid category: {sig['category']} for {sig['name']}"
            )

    def test_all_risk_levels_valid(self):
        """All signatures have valid risk levels."""
        valid_risks = {"low", "medium", "high", "critical"}
        for sig in _HAZARD_SIGNATURES:
            assert sig["risk"] in valid_risks, (
                f"Invalid risk: {sig['risk']} for {sig['name']}"
            )

    def test_all_confidence_values_in_range(self):
        """All signature confidence values are in [0.0, 1.0]."""
        for sig in _HAZARD_SIGNATURES:
            assert 0.0 <= sig["confidence"] <= 1.0, (
                f"Confidence out of range for {sig['name']}: {sig['confidence']}"
            )

    def test_all_types_are_protein_or_dna(self):
        """All signatures have type 'protein' or 'dna'."""
        for sig in _HAZARD_SIGNATURES:
            assert sig["type"] in ("protein", "dna"), (
                f"Invalid type: {sig['type']} for {sig['name']}"
            )

    def test_all_motifs_non_empty(self):
        """All signature motifs are non-empty."""
        for sig in _HAZARD_SIGNATURES:
            assert len(sig["motif"]) > 0, (
                f"Empty motif for {sig['name']}"
            )

    def test_select_agent_signatures_are_medium_or_above(self):
        """Select agent signatures should be medium, high, or critical risk.

        Note: Some select agent entries (e.g., T-2 mycotoxin target RPL3)
        are classified as 'medium' risk because they target host proteins
        rather than being direct toxin signatures.
        """
        for sig in _HAZARD_SIGNATURES:
            if sig["category"] == "select_agent":
                assert sig["risk"] in ("medium", "high", "critical"), (
                    f"Select agent {sig['name']} should be medium/high/critical, "
                    f"got {sig['risk']}"
                )


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 9: Report Structure and Recommendations
# ═══════════════════════════════════════════════════════════════════


class TestReportStructureAndRecommendations:
    """Verify BiosecurityReport structure and recommendation quality."""

    def test_report_fields_present(self):
        """BiosecurityReport has all required fields."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        assert hasattr(report, "is_hazardous")
        assert hasattr(report, "risk_level")
        assert hasattr(report, "flagged_categories")
        assert hasattr(report, "matches")
        assert hasattr(report, "recommendations")

    def test_hazard_match_fields_present(self):
        """Each HazardMatch has all required fields."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        for match in report.matches:
            assert hasattr(match, "category")
            assert hasattr(match, "name")
            assert hasattr(match, "position")
            assert hasattr(match, "matched_sequence")
            assert hasattr(match, "confidence")
            assert hasattr(match, "source")

    def test_match_confidence_in_range(self):
        """All match confidence values are in [0.0, 1.0]."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        for match in report.matches:
            assert 0.0 <= match.confidence <= 1.0

    def test_match_positions_non_negative(self):
        """All match positions are non-negative."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        for match in report.matches:
            assert match.position >= 0

    def test_flagged_categories_sorted(self):
        """Flagged categories are sorted."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        assert report.flagged_categories == sorted(report.flagged_categories)

    def test_recommendations_for_critical_hazard(self):
        """Critical hazard has actionable recommendations."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        assert len(report.recommendations) > 0
        rec_text = " ".join(report.recommendations).lower()
        assert "select agent" in rec_text or "toxin" in rec_text or "critical" in rec_text

    def test_recommendations_for_safe_sequence(self):
        """Safe sequence has 'no concerns' recommendation."""
        report = screen_hazardous_sequence(INSULIN_PROTEIN)
        assert len(report.recommendations) > 0
        assert "No biosecurity concerns" in report.recommendations[0]

    def test_recommendations_mention_biosafety_officer(self):
        """Recommendations for critical risk mention biosafety officer."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        rec_text = " ".join(report.recommendations).lower()
        assert "biosafety officer" in rec_text or "ibo" in rec_text or "42 cfr" in rec_text

    def test_recommendations_for_antibiotic_resistance(self):
        """Recommendations for antibiotic resistance mention NIH guidelines."""
        report = screen_hazardous_sequence(BLATEM_PROTEIN_MOTIF)
        rec_text = " ".join(report.recommendations).lower()
        assert "resistance" in rec_text or "antibiotic" in rec_text


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 10: Risk Level Helper Functions
# ═══════════════════════════════════════════════════════════════════


class TestRiskLevelHelpers:
    """Test internal risk level helper functions."""

    def test_max_risk_empty(self):
        """_max_risk with no arguments returns 'none'."""
        assert _max_risk() == "none"

    def test_max_risk_single(self):
        """_max_risk with single level returns that level."""
        assert _max_risk("low") == "low"
        assert _max_risk("critical") == "critical"

    def test_max_risk_multiple(self):
        """_max_risk returns the highest risk level."""
        assert _max_risk("low", "medium") == "medium"
        assert _max_risk("low", "critical") == "critical"
        assert _max_risk("none", "low", "medium", "high", "critical") == "critical"

    def test_max_risk_ordering(self):
        """Risk level ordering: none < low < medium < high < critical."""
        assert _max_risk("none", "low") == "low"
        assert _max_risk("low", "medium") == "medium"
        assert _max_risk("medium", "high") == "high"
        assert _max_risk("high", "critical") == "critical"

    def test_sig_risk_for_match(self):
        """sig_risk_for_match returns correct risk for known matches."""
        # Create a match with a known signature
        match = HazardMatch(
            category="select_agent",
            name="ricin_A_chain_catalytic",
            position=0,
            matched_sequence="NIRVGLPIIS",
            confidence=0.95,
            source="test",
        )
        risk = sig_risk_for_match(match)
        assert risk == "critical"

    def test_sig_risk_for_unknown_match(self):
        """sig_risk_for_match returns category default for unknown match."""
        match = HazardMatch(
            category="select_agent",
            name="nonexistent_toxin",
            position=0,
            matched_sequence="AAAAAAAAAA",
            confidence=0.5,
            source="test",
        )
        risk = sig_risk_for_match(match)
        assert risk == "critical"  # Fallback for select_agent


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 11: Viral Surface Protein Detection
# ═══════════════════════════════════════════════════════════════════


class TestViralSurfaceProteinDetection:
    """Verify viral surface protein detection across multiple pathogens."""

    def test_influenza_ha_detected(self):
        """Influenza HA fusion peptide is detected."""
        report = screen_hazardous_sequence(INFLUENZA_HA_MOTIF)
        assert "viral_surface" in report.flagged_categories
        assert any("influenza_HA" in m.name for m in report.matches)

    def test_sars2_spike_rbd_detected(self):
        """SARS-CoV-2 spike RBD is detected."""
        report = screen_hazardous_sequence(SARS2_RBD_MOTIF)
        assert "viral_surface" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "SARS2_spike_RBD" in names

    def test_hiv_env_v3_detected(self):
        """HIV-1 Env V3 loop is detected."""
        report = screen_hazardous_sequence(HIV_V3_MOTIF)
        assert "viral_surface" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "HIV_env_V3_loop" in names

    def test_ebola_gp_detected(self):
        """Ebola GP1 receptor binding motif is detected."""
        protein = "AAAWIPVQNQCGPAAAA"
        report = screen_hazardous_sequence(protein)
        assert "viral_surface" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "ebola_GP1_receptor" in names

    def test_variola_envelope_detected(self):
        """Variola (smallpox) envelope protein is detected."""
        protein = "AAAYDDVVRVYKVAAA"
        report = screen_hazardous_sequence(protein)
        assert "viral_surface" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "variola_envelope" in names

    def test_viral_surface_recommendations_mention_durc(self):
        """Viral surface protein recommendations mention DURC policies."""
        report = screen_hazardous_sequence(INFLUENZA_HA_MOTIF)
        rec_text = " ".join(report.recommendations).lower()
        assert "dual-use" in rec_text or "durc" in rec_text or "viral" in rec_text or "vaccine" in rec_text


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 12: Oncogene Detection
# ═══════════════════════════════════════════════════════════════════


class TestOncogeneDetection:
    """Verify oncogene and growth factor detection."""

    def test_ras_gtp_binding_detected(self):
        """RAS GTP-binding motif is detected."""
        protein = "AAALVGNKCDLPSAAA"
        report = screen_hazardous_sequence(protein)
        assert "oncogene" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "RAS_GTP_binding" in names

    def test_braf_activation_detected(self):
        """BRAF activation segment motif is detected."""
        protein = "AAAIGDFGLATVKAAA"
        report = screen_hazardous_sequence(protein)
        assert "oncogene" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "BRAF_activation" in names

    def test_egfr_kinase_detected(self):
        """EGFR kinase domain motif is detected."""
        protein = "AAAIKHRDLAARNAAA"
        report = screen_hazardous_sequence(protein)
        assert "oncogene" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "EGFR_kinase" in names

    def test_oncogene_recommendations_mention_research(self):
        """Oncogene recommendations mention legitimate research."""
        protein = "AAALVGNKCDLPSAAA"
        report = screen_hazardous_sequence(protein)
        rec_text = " ".join(report.recommendations).lower()
        assert "research" in rec_text or "oncogene" in rec_text


# ═══════════════════════════════════════════════════════════════════
# Agent 50 — Test Suite 13: Case Insensitivity and Edge Cases
# ═══════════════════════════════════════════════════════════════════


class TestCaseInsensitivityAndEdgeCases:
    """Test case insensitivity and edge cases for biosecurity screening."""

    def test_lowercase_protein_detected(self):
        """Lowercase protein sequence is detected correctly."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF.lower())
        assert report.is_hazardous
        assert any("ricin" in m.name for m in report.matches)

    def test_mixed_case_protein_detected(self):
        """Mixed-case protein sequence is detected correctly."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF.title())
        assert report.is_hazardous

    def test_lowercase_dna_detected(self):
        """Lowercase DNA sequence is detected correctly."""
        report = screen_hazardous_sequence("", BLATEM_DNA.lower())
        assert report.is_hazardous

    def test_multiple_occurrences_detected(self):
        """Multiple occurrences of the same motif are all detected."""
        protein = "NIRVGLPIISNIRVGLPIIS"
        report = screen_hazardous_sequence(protein)
        ricin_matches = [m for m in report.matches if m.name == "ricin_A_chain_catalytic"]
        assert len(ricin_matches) >= 2

    def test_motif_at_start_of_protein(self):
        """Motif at the very start of a protein is detected."""
        protein = "NIRVGLPIISAAAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        ricin_matches = [m for m in report.matches if m.name == "ricin_A_chain_catalytic"]
        assert len(ricin_matches) >= 1
        assert ricin_matches[0].position == 0

    def test_motif_at_end_of_protein(self):
        """Motif at the very end of a protein is detected."""
        protein = "AAAAAANIRVGLPIIS"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous

    def test_dna_only_screening(self):
        """Screening with only DNA (empty protein) works."""
        report = screen_hazardous_sequence("", BLATEM_DNA)
        assert report.is_hazardous
        assert "antibiotic_resistance" in report.flagged_categories

    def test_medium_risk_emits_warning(self):
        """Medium risk check emits a UserWarning."""
        protein = NPTII_MOTIF
        report = screen_hazardous_sequence(protein)
        if report.risk_level == "medium":
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = check_biosecurity_before_optimize(protein)
                assert len(w) >= 1
                assert any("Biosecurity" in str(warning.message) for warning in w)

    def test_medium_risk_does_not_raise(self):
        """Medium risk does not raise BiosecurityError."""
        protein = NPTII_MOTIF
        report = screen_hazardous_sequence(protein)
        if report.risk_level == "medium":
            result = check_biosecurity_before_optimize(protein)
            assert isinstance(result, BiosecurityReport)
            assert result.risk_level == "medium"

    def test_low_risk_no_warning_no_error(self):
        """Low risk does not warn or raise (only logs)."""
        protein = "AAAAPMAVPPPKKAAA"
        report = screen_hazardous_sequence(protein)
        if report.risk_level == "low":
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = check_biosecurity_before_optimize(protein)
                # Low risk should not emit a UserWarning (only logs)
                biosec_warnings = [
                    warning for warning in w
                    if "Biosecurity" in str(warning.message)
                ]
                assert len(biosec_warnings) == 0
                assert isinstance(result, BiosecurityReport)

    def test_is_hazardous_false_for_none(self):
        """is_hazardous is False for 'none' risk level."""
        report = screen_hazardous_sequence(INSULIN_PROTEIN)
        assert report.is_hazardous is False

    def test_is_hazardous_true_for_medium(self):
        """is_hazardous is True for 'medium' risk level."""
        protein = NPTII_MOTIF
        report = screen_hazardous_sequence(protein)
        if report.risk_level == "medium":
            assert report.is_hazardous is True

    def test_is_hazardous_true_for_critical(self):
        """is_hazardous is True for 'critical' risk level."""
        report = screen_hazardous_sequence(RICIN_A_MOTIF)
        assert report.is_hazardous is True
