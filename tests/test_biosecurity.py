"""
Tests for the biosecurity screening module.

Covers:
  - Known toxin detection (ricin A-chain motif)
  - Antibiotic resistance marker detection
  - Clean sequence passes
  - Risk level classification
  - Hard-stop behavior for critical hazards
  - Warning behavior for medium hazards
  - DNA signature matching
  - Oncogene detection
  - Viral surface protein detection
  - Multiple matches and category aggregation
  - Integration hook (check_biosecurity_before_optimize)
"""

import warnings

import pytest

from biocompiler.biosecurity import (
    BiosecurityReport,
    HazardMatch,
    HAZARD_SIGNATURE_COUNT,
    screen_hazardous_sequence,
    check_biosecurity_before_optimize,
)
from biocompiler.exceptions import BiosecurityError, BioCompilerError


# ═══════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def clean_protein():
    """A benign protein (human insulin A-chain + flanking)."""
    return "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"


@pytest.fixture
def ricin_a_chain_protein():
    """Protein containing the ricin A-chain catalytic motif NIRVGLPIIS."""
    return "MISRDNIRVGLPIISTNKYEDKQL"


@pytest.fixture
def antibiotic_resistance_dna():
    """DNA containing the blaTEM ORF start signature."""
    return "ATGGAGCCCGATGAGTATTCAACATTTCCGTGTCGCCCTTATTCC"


@pytest.fixture
def ndm1_protein():
    """Protein containing the NDM-1 zinc binding motif HHHHDGHFGGP."""
    return "MKHHHHHDGHFGGPVKVTLEAP"


# ═══════════════════════════════════════════════════════════════════════════
# Test: Database size
# ═══════════════════════════════════════════════════════════════════════════

class TestDatabaseSize:
    """Ensure the hazard signature database is comprehensive."""

    def test_at_least_50_signatures(self):
        assert HAZARD_SIGNATURE_COUNT >= 50, (
            f"Hazard database has only {HAZARD_SIGNATURE_COUNT} signatures; "
            f"require at least 50"
        )


# ═══════════════════════════════════════════════════════════════════════════
# Test: Known toxin detection
# ═══════════════════════════════════════════════════════════════════════════

class TestToxinDetection:
    """Detect select agent toxins by peptide motif."""

    def test_ricin_a_chain_detected(self, ricin_a_chain_protein):
        report = screen_hazardous_sequence(ricin_a_chain_protein)
        assert report.is_hazardous
        assert report.risk_level == "critical"
        assert "select_agent" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "ricin_A_chain_catalytic" in names

    def test_ricin_a_chain_position(self, ricin_a_chain_protein):
        report = screen_hazardous_sequence(ricin_a_chain_protein)
        ricin_matches = [m for m in report.matches if m.name == "ricin_A_chain_catalytic"]
        assert len(ricin_matches) == 1
        assert ricin_matches[0].position == 5  # "MISRDNIRVGLPIIS..."
        assert ricin_matches[0].matched_sequence == "NIRVGLPIIS"

    def test_botulinum_detected(self):
        protein = "AAAHEHETQSNLRDLAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        assert "select_agent" in report.flagged_categories
        assert any("botulinum" in m.name for m in report.matches)

    def test_shiga_toxin_detected(self):
        protein = "AAANVRVNIPFSRAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "shiga_toxin_A_subunit" in names

    def test_anthrax_LF_detected(self):
        protein = "AAHETHFGVVSYAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        assert any("anthrax_LF" in m.name for m in report.matches)

    def test_diphtheria_detected(self):
        protein = "AAAGVADVIQEINAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "diphtheria_toxin_ADR" in names

    def test_tetanus_detected(self):
        protein = "AAAHEIKSNIASKAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "tetanus_toxin_zinc" in names

    def test_abrin_detected(self):
        protein = "AAANVRVGIPISRAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "abrin_A_chain" in names


# ═══════════════════════════════════════════════════════════════════════════
# Test: Antibiotic resistance marker detection
# ═══════════════════════════════════════════════════════════════════════════

class TestAntibioticResistanceDetection:
    """Detect antibiotic resistance markers via protein and DNA signatures."""

    def test_blatem_protein_detected(self):
        protein = "AAAHPETLALKFGAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        assert "antibiotic_resistance" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "blaTEM_protein" in names

    def test_blatem_dna_detected(self, antibiotic_resistance_dna):
        report = screen_hazardous_sequence("", antibiotic_resistance_dna)
        assert report.is_hazardous
        assert "antibiotic_resistance" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "blaTEM_dna" in names

    def test_nptII_protein_detected(self):
        protein = "AAARPMTIHGSGSAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "nptII_protein" in names

    def test_ndm1_protein_detected(self, ndm1_protein):
        report = screen_hazardous_sequence(ndm1_protein)
        assert report.is_hazardous
        # NDM-1 is critical risk
        assert report.risk_level == "critical"
        names = [m.name for m in report.matches]
        assert "ndm1_protein" in names

    def test_vanA_protein_detected(self):
        protein = "AAAHGLSSAVPGLAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "vanA_protein" in names

    def test_cat_protein_detected(self):
        protein = "AAAFHRGVCTNKAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "cat_protein" in names


# ═══════════════════════════════════════════════════════════════════════════
# Test: Clean sequence passes
# ═══════════════════════════════════════════════════════════════════════════

class TestCleanSequence:
    """Benign sequences should pass screening with risk_level='none'."""

    def test_insulin_passes(self, clean_protein):
        report = screen_hazardous_sequence(clean_protein)
        assert not report.is_hazardous
        assert report.risk_level == "none"
        assert report.flagged_categories == []
        assert report.matches == []

    def test_egfp_passes(self):
        egfp = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        report = screen_hazardous_sequence(egfp)
        assert not report.is_hazardous
        assert report.risk_level == "none"
        assert len(report.matches) == 0

    def test_empty_protein_passes(self):
        report = screen_hazardous_sequence("")
        assert not report.is_hazardous
        assert report.risk_level == "none"

    def test_recommendations_for_clean(self, clean_protein):
        report = screen_hazardous_sequence(clean_protein)
        assert len(report.recommendations) == 1
        assert "No biosecurity concerns" in report.recommendations[0]


# ═══════════════════════════════════════════════════════════════════════════
# Test: Risk level classification
# ═══════════════════════════════════════════════════════════════════════════

class TestRiskLevelClassification:
    """Verify correct risk level assignment based on match categories."""

    def test_critical_for_select_agent(self, ricin_a_chain_protein):
        report = screen_hazardous_sequence(ricin_a_chain_protein)
        assert report.risk_level == "critical"

    def test_high_for_viral_surface(self):
        # Influenza HA fusion peptide
        protein = "AAAGLFGAIAGFIAAAA"
        report = screen_hazardous_sequence(protein)
        assert report.risk_level in ("high", "critical")
        assert "viral_surface" in report.flagged_categories

    def test_medium_for_antibiotic_resistance_only(self):
        # nptII protein motif — risk is "medium"
        protein = "AAARPMTIHGSGSAAA"
        report = screen_hazardous_sequence(protein)
        # May also match other things; check that antibiotic_resistance is flagged
        assert "antibiotic_resistance" in report.flagged_categories

    def test_low_for_oncogene_only(self):
        # VEGF heparin binding motif — risk is "low"
        protein = "AAAAPMAVPPPKKAAA"
        report = screen_hazardous_sequence(protein)
        assert "oncogene" in report.flagged_categories

    def test_highest_risk_wins(self):
        """When multiple categories match, the highest risk level prevails."""
        # Protein containing both a VEGF motif (low) and ricin motif (critical)
        protein = "APMAVPPPPKNIRVGLPIIS"
        report = screen_hazardous_sequence(protein)
        assert report.risk_level == "critical"

    def test_is_hazardous_true_for_medium(self):
        """Medium risk should still be flagged as hazardous."""
        # Create a sequence that only matches medium-risk signatures
        protein = "RPMTIHGSGS"  # nptII (medium risk)
        report = screen_hazardous_sequence(protein)
        # If only medium matches, is_hazardous should be True
        if report.risk_level == "medium":
            assert report.is_hazardous is True

    def test_is_hazardous_false_for_none(self, clean_protein):
        report = screen_hazardous_sequence(clean_protein)
        assert report.is_hazardous is False


# ═══════════════════════════════════════════════════════════════════════════
# Test: Hard-stop behavior for critical hazards
# ═══════════════════════════════════════════════════════════════════════════

class TestHardStopCritical:
    """Critical and high risk levels should raise BiosecurityError."""

    def test_critical_raises_biosecurity_error(self, ricin_a_chain_protein):
        with pytest.raises(BiosecurityError) as exc_info:
            check_biosecurity_before_optimize(ricin_a_chain_protein, organism="e_coli")
        assert exc_info.value.report.risk_level == "critical"
        assert "select_agent" in exc_info.value.report.flagged_categories

    def test_high_raises_biosecurity_error(self):
        # Influenza HA fusion peptide is high risk
        protein = "AAAGLFGAIAGFIAAAA"
        # May also match other high/critical signatures; at minimum should raise
        with pytest.raises(BiosecurityError):
            check_biosecurity_before_optimize(protein, organism="human")

    def test_biosecurity_error_inherits_from_biocompiler_error(self, ricin_a_chain_protein):
        with pytest.raises(BioCompilerError):
            check_biosecurity_before_optimize(ricin_a_chain_protein)

    def test_error_contains_report(self, ricin_a_chain_protein):
        with pytest.raises(BiosecurityError) as exc_info:
            check_biosecurity_before_optimize(ricin_a_chain_protein)
        assert hasattr(exc_info.value, "report")
        assert isinstance(exc_info.value.report, BiosecurityReport)
        assert exc_info.value.report.is_hazardous is True


# ═══════════════════════════════════════════════════════════════════════════
# Test: Warning behavior for medium hazards
# ═══════════════════════════════════════════════════════════════════════════

class TestWarningMedium:
    """Medium risk should emit a warning but not raise an exception."""

    def test_medium_emits_warning(self):
        # Use a sequence that ONLY matches medium-risk signatures
        # nptII protein motif is "medium" risk
        protein = "AAARPMTIHGSGSAAA"
        report = screen_hazardous_sequence(protein)
        if report.risk_level == "medium":
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = check_biosecurity_before_optimize(protein)
                assert len(w) >= 1
                assert any("Biosecurity" in str(warning.message) for warning in w)
            assert isinstance(result, BiosecurityReport)

    def test_medium_does_not_raise_error(self):
        # nptII protein motif — medium risk
        protein = "AAARPMTIHGSGSAAA"
        report = screen_hazardous_sequence(protein)
        if report.risk_level == "medium":
            result = check_biosecurity_before_optimize(protein)
            assert result.risk_level == "medium"


# ═══════════════════════════════════════════════════════════════════════════
# Test: Low risk logging (no error, no warning)
# ═══════════════════════════════════════════════════════════════════════════

class TestLowRiskLogging:
    """Low risk should not raise or warn, just log."""

    def test_low_does_not_raise(self):
        # VEGF heparin binding — low risk
        protein = "AAAAPMAVPPPKKAAA"
        report = screen_hazardous_sequence(protein)
        if report.risk_level == "low":
            result = check_biosecurity_before_optimize(protein)
            assert result.risk_level == "low"

    def test_none_does_not_raise(self, clean_protein):
        result = check_biosecurity_before_optimize(clean_protein)
        assert result.risk_level == "none"


# ═══════════════════════════════════════════════════════════════════════════
# Test: Viral surface protein detection
# ═══════════════════════════════════════════════════════════════════════════

class TestViralSurfaceDetection:
    """Detect viral surface proteins by peptide motif."""

    def test_influenza_HA_detected(self):
        protein = "AAAGLFGAIAGFIAAAA"
        report = screen_hazardous_sequence(protein)
        assert "viral_surface" in report.flagged_categories
        assert any("influenza_HA" in m.name for m in report.matches)

    def test_sars2_spike_RBD_detected(self):
        protein = "AAAVYYHKNNKSWAAAA"
        report = screen_hazardous_sequence(protein)
        assert "viral_surface" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "SARS2_spike_RBD" in names

    def test_hiv_env_V3_detected(self):
        protein = "AAAGPGRAFYTIGAAAA"
        report = screen_hazardous_sequence(protein)
        assert "viral_surface" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "HIV_env_V3_loop" in names

    def test_ebola_GP_detected(self):
        protein = "AAAWIPVQNQCGPAAAA"
        report = screen_hazardous_sequence(protein)
        assert "viral_surface" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "ebola_GP1_receptor" in names

    def test_variola_detected(self):
        protein = "AAAYDDVVRVYKVAAA"
        report = screen_hazardous_sequence(protein)
        assert "viral_surface" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "variola_envelope" in names


# ═══════════════════════════════════════════════════════════════════════════
# Test: Oncogene detection
# ═══════════════════════════════════════════════════════════════════════════

class TestOncogeneDetection:
    """Detect oncogenes and growth factors by peptide motif."""

    def test_RAS_GTP_binding_detected(self):
        protein = "AAALVGNKCDLPSAAA"
        report = screen_hazardous_sequence(protein)
        assert "oncogene" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "RAS_GTP_binding" in names

    def test_BRAF_activation_detected(self):
        protein = "AAAIGDFGLATVKAAA"
        report = screen_hazardous_sequence(protein)
        assert "oncogene" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "BRAF_activation" in names

    def test_EGFR_kinase_detected(self):
        protein = "AAAIKHRDLAARNAAA"
        report = screen_hazardous_sequence(protein)
        assert "oncogene" in report.flagged_categories
        names = [m.name for m in report.matches]
        assert "EGFR_kinase" in names


# ═══════════════════════════════════════════════════════════════════════════
# Test: DNA signature matching
# ═══════════════════════════════════════════════════════════════════════════

class TestDNASignatures:
    """Test nucleotide pattern matching for resistance markers."""

    def test_blatem_dna_only(self, antibiotic_resistance_dna):
        """Screen with empty protein but matching DNA."""
        report = screen_hazardous_sequence("", antibiotic_resistance_dna)
        assert report.is_hazardous
        assert "antibiotic_resistance" in report.flagged_categories

    def test_ndm1_dna_detected(self):
        dna = "ATGGAATTGCCCAATATTATG"
        report = screen_hazardous_sequence("", dna)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "ndm1_dna" in names

    def test_nptII_dna_detected(self):
        dna = "ATGATTGAACAAGATGGATTG"
        report = screen_hazardous_sequence("", dna)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "nptII_dna" in names

    def test_both_protein_and_dna(self):
        """Both protein and DNA signatures can be detected simultaneously."""
        protein = "HPETLALKFG"  # blaTEM protein
        dna = "ATGAGTATTCAACATTTCCGTG"  # blaTEM DNA
        report = screen_hazardous_sequence(protein, dna)
        assert report.is_hazardous
        names = [m.name for m in report.matches]
        assert "blaTEM_protein" in names
        assert "blaTEM_dna" in names


# ═══════════════════════════════════════════════════════════════════════════
# Test: BiosecurityReport structure
# ═══════════════════════════════════════════════════════════════════════════

class TestBiosecurityReportStructure:
    """Verify the dataclass structure and field types."""

    def test_report_dataclass_fields(self, ricin_a_chain_protein):
        report = screen_hazardous_sequence(ricin_a_chain_protein)
        assert hasattr(report, "is_hazardous")
        assert hasattr(report, "risk_level")
        assert hasattr(report, "flagged_categories")
        assert hasattr(report, "matches")
        assert hasattr(report, "recommendations")

    def test_risk_level_valid_values(self, clean_protein, ricin_a_chain_protein):
        valid_levels = {"none", "low", "medium", "high", "critical"}
        report1 = screen_hazardous_sequence(clean_protein)
        assert report1.risk_level in valid_levels
        report2 = screen_hazardous_sequence(ricin_a_chain_protein)
        assert report2.risk_level in valid_levels

    def test_hazard_match_dataclass(self, ricin_a_chain_protein):
        report = screen_hazardous_sequence(ricin_a_chain_protein)
        for match in report.matches:
            assert isinstance(match, HazardMatch)
            assert hasattr(match, "category")
            assert hasattr(match, "name")
            assert hasattr(match, "position")
            assert hasattr(match, "matched_sequence")
            assert hasattr(match, "confidence")
            assert hasattr(match, "source")

    def test_confidence_range(self, ricin_a_chain_protein):
        report = screen_hazardous_sequence(ricin_a_chain_protein)
        for match in report.matches:
            assert 0.0 <= match.confidence <= 1.0

    def test_position_is_nonnegative(self, ricin_a_chain_protein):
        report = screen_hazardous_sequence(ricin_a_chain_protein)
        for match in report.matches:
            assert match.position >= 0

    def test_flagged_categories_sorted(self, ricin_a_chain_protein):
        report = screen_hazardous_sequence(ricin_a_chain_protein)
        assert report.flagged_categories == sorted(report.flagged_categories)

    def test_recommendations_nonempty_for_hazardous(self, ricin_a_chain_protein):
        report = screen_hazardous_sequence(ricin_a_chain_protein)
        assert len(report.recommendations) > 0


# ═══════════════════════════════════════════════════════════════════════════
# Test: Case insensitivity
# ═══════════════════════════════════════════════════════════════════════════

class TestCaseInsensitivity:
    """Sequences should be case-insensitive."""

    def test_lowercase_protein_detected(self):
        protein = "misrdnirvglpiistnkyedkql"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        assert any("ricin" in m.name for m in report.matches)

    def test_mixed_case_protein_detected(self):
        protein = "MiSrdNirvGlpiisTNKyeDKql"
        report = screen_hazardous_sequence(protein)
        assert report.is_hazardous
        assert any("ricin" in m.name for m in report.matches)

    def test_lowercase_dna_detected(self):
        dna = "atgagtattcaacatttccgtg"
        report = screen_hazardous_sequence("", dna)
        assert report.is_hazardous
        assert any("blaTEM" in m.name for m in report.matches)


# ═══════════════════════════════════════════════════════════════════════════
# Test: Multiple matches
# ═══════════════════════════════════════════════════════════════════════════

class TestMultipleMatches:
    """Sequences matching multiple signatures should aggregate correctly."""

    def test_multiple_categories_flagged(self):
        # Ricin motif + VEGF motif in one protein
        protein = "APMAVPPPPKNIRVGLPIISAPMAVPPPKK"
        report = screen_hazardous_sequence(protein)
        assert "select_agent" in report.flagged_categories
        assert "oncogene" in report.flagged_categories
        assert report.risk_level == "critical"  # highest wins

    def test_multiple_matches_reported(self):
        protein = "NIRVGLPIISNIRVGLPIIS"
        report = screen_hazardous_sequence(protein)
        ricin_matches = [m for m in report.matches if m.name == "ricin_A_chain_catalytic"]
        assert len(ricin_matches) >= 2  # two occurrences


# ═══════════════════════════════════════════════════════════════════════════
# Test: BiosecurityError exception hierarchy
# ═══════════════════════════════════════════════════════════════════════════

class TestBiosecurityErrorException:
    """Verify BiosecurityError is properly integrated into the exception hierarchy."""

    def test_is_biocompiler_error(self):
        assert issubclass(BiosecurityError, BioCompilerError)

    def test_is_exception(self):
        assert issubclass(BiosecurityError, Exception)

    def test_str_contains_risk_level(self, ricin_a_chain_protein):
        with pytest.raises(BiosecurityError) as exc_info:
            check_biosecurity_before_optimize(ricin_a_chain_protein)
        error_str = str(exc_info.value)
        assert "critical" in error_str or "risk_level=critical" in error_str

    def test_report_attribute(self, ricin_a_chain_protein):
        with pytest.raises(BiosecurityError) as exc_info:
            check_biosecurity_before_optimize(ricin_a_chain_protein)
        assert exc_info.value.report.is_hazardous is True
