"""
BioCompiler CamSol IDP Sensitivity Tests
=========================================
Tests for the Urry hydrophobicity scale integration and IDP sensitivity
improvement. Validates that:

  1. The Urry scale is auto-selected for IDP sequences
  2. Known IDPs (alpha-synuclein, tau, p53 TAD) score as aggregation-prone
  3. Sensitivity improves from 33% (3/9) to >= 67% (6/9)
  4. Soluble proteins are not adversely affected by the new scale
  5. The IDP prediction heuristic works correctly
"""

import pytest

from biocompiler.camsol import (
    CAMSOL_HYDROPATHY,
    URRY_HYDROPATHY,
    HydrophobicityScale,
    CamSolResult,
    classify_solubility,
    clear_cache,
    compute_intrinsic_solubility,
    predict_idp,
    select_hydropathy_scale,
)
from biocompiler.validation.camsol_benchmark import (
    BENCHMARK_DATASET,
    compute_enhanced_benchmark_score,
)


# ────────────────────────────────────────────────────────────
# Known IDP sequences for testing
# ────────────────────────────────────────────────────────────

# Alpha-synuclein (140 aa) — Parkinson's disease amyloid
ALPHA_SYNUCLEIN = (
    "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGA"
    "VVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
)

# Tau K18 fragment (155 aa) — tauopathy amyloid core
TAU_K18 = (
    "VQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETH"
    "KLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
)

# p53 Transactivation Domain (residues 1-73) — intrinsically disordered
P53_TAD = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"

# Huntingtin Exon-1 (17Q) — polyQ expansion aggregates
HUNTINGTIN_EX1 = (
    "MATLEKLMKAFESLKSFQQQQQQQQQQQQQQQQQQPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"
    "PPPPPPPPPPPPPPPPPPPP"
)

# Amyloid-beta 42 — prototypical amyloid
AB42 = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"

# Human Prion Protein (mature 23-231) — prion disease
PRION = (
    "QVYYRPVDQYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRG"
    "SSMVLFSSPPVILLISFLIFLIVG"
)

# Human IAPP (37 aa) — type 2 diabetes amyloid
IAPP = "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"

# Prothymosin alpha — highly disordered, very acidic
PROTHYMOSIN_ALPHA = (
    "MSDAAVDTSSEITTKDLKEKKEVVEEAENGRDAPANGNNENEENGEQEADNEVDEEEEEGGEEEEEEEEEG"
    "DGEEEDGDEEAEIAKEAEGEEDEEEEGDEEAGDEEG"
)

# c-Myc Transactivation Domain — well-known IDP
CMYC_TAD = (
    "MPLNVSFNRDADLQKLSVEQFEMRLQEIQQSIQLRQEEENRQSSSSESIIPNVSKRTRSSSSLQSIKDEI"
    "LNRLLQHEMLNSVEDNRLNHRV"
)

# Well-known soluble/globular proteins for specificity testing
THIOREDOXIN = (
    "MSTKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIKAPILDEIADEYQGKLTVAKLNIDQNPGTAP"
    "KYGIRGIPTLLFKNGEVAATKVGALSKGQLKEFLDANLAGS"
)

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)


# ────────────────────────────────────────────────────────────
# Test: Urry scale data integrity
# ────────────────────────────────────────────────────────────

class TestUrryScaleIntegrity:
    """Verify the Urry hydrophobicity scale is well-formed."""

    def test_all_20_standard_aas_present(self):
        """URRY_HYDROPATHY should contain all 20 standard amino acids."""
        standard_aas = set("ACDEFGHIKLMNPQRSTVWY")
        assert set(URRY_HYDROPATHY.keys()) == standard_aas

    def test_hydrophobic_residues_negative(self):
        """Hydrophobic residues should have negative values (CamSol convention)."""
        hydrophobic = {"I", "L", "V", "F", "W", "Y", "M"}
        for aa in hydrophobic:
            assert URRY_HYDROPATHY[aa] < 0, (
                f"Hydrophobic residue {aa} should be negative in Urry scale, "
                f"got {URRY_HYDROPATHY[aa]}"
            )

    def test_charged_residues_positive(self):
        """Charged residues should have positive values (CamSol convention)."""
        charged = {"K", "R", "D", "E"}
        for aa in charged:
            assert URRY_HYDROPATHY[aa] > 0, (
                f"Charged residue {aa} should be positive in Urry scale, "
                f"got {URRY_HYDROPATHY[aa]}"
            )

    def test_urry_more_extreme_than_ww(self):
        """Urry scale should have more extreme values for hydrophobic residues."""
        for aa in ["I", "L", "V", "F", "W", "Y"]:
            assert URRY_HYDROPATHY[aa] < CAMSOL_HYDROPATHY[aa], (
                f"Urry value for {aa} ({URRY_HYDROPATHY[aa]}) should be more "
                f"negative than Wimley-White ({CAMSOL_HYDROPATHY[aa]})"
            )

    def test_urry_less_extreme_charged(self):
        """Urry scale should have less extreme positive values for charged residues."""
        for aa in ["K", "E", "D"]:
            assert URRY_HYDROPATHY[aa] < CAMSOL_HYDROPATHY[aa], (
                f"Urry value for {aa} ({URRY_HYDROPATHY[aa]}) should be less "
                f"positive than Wimley-White ({CAMSOL_HYDROPATHY[aa]})"
            )


# ────────────────────────────────────────────────────────────
# Test: HydrophobicityScale enum
# ────────────────────────────────────────────────────────────

class TestHydrophobicityScaleEnum:
    """Test the HydrophobicityScale enum."""

    def test_enum_values(self):
        """Enum should have WIMLEY_WHITE and URRY members."""
        assert HydrophobicityScale.WIMLEY_WHITE.value == "wimley_white"
        assert HydrophobicityScale.URRY.value == "urry"

    def test_enum_is_string(self):
        """Enum members should be strings."""
        assert isinstance(HydrophobicityScale.WIMLEY_WHITE, str)
        assert isinstance(HydrophobicityScale.URRY, str)


# ────────────────────────────────────────────────────────────
# Test: IDP prediction heuristic
# ────────────────────────────────────────────────────────────

class TestIDPPrediction:
    """Test the IDP prediction heuristic (predict_idp)."""

    def test_alpha_synuclein_predicted_idp(self):
        """Alpha-synuclein should be predicted as IDP."""
        assert predict_idp(ALPHA_SYNUCLEIN), (
            "Alpha-synuclein should be predicted as IDP"
        )

    def test_tau_k18_predicted_idp(self):
        """Tau K18 fragment should be predicted as IDP."""
        assert predict_idp(TAU_K18), (
            "Tau K18 fragment should be predicted as IDP"
        )

    def test_huntingtin_predicted_idp(self):
        """Huntingtin exon-1 should be predicted as IDP."""
        assert predict_idp(HUNTINGTIN_EX1), (
            "Huntingtin exon-1 should be predicted as IDP"
        )

    def test_p53_tad_predicted_idp(self):
        """p53 TAD should be predicted as IDP."""
        assert predict_idp(P53_TAD), (
            "p53 transactivation domain should be predicted as IDP"
        )

    def test_prothymosin_predicted_idp(self):
        """Prothymosin alpha should be predicted as IDP (highly acidic)."""
        assert predict_idp(PROTHYMOSIN_ALPHA), (
            "Prothymosin alpha should be predicted as IDP"
        )

    def test_thioredoxin_not_predicted_idp(self):
        """E. coli Thioredoxin should NOT be predicted as IDP."""
        assert not predict_idp(THIOREDOXIN), (
            "Thioredoxin should NOT be predicted as IDP"
        )

    def test_ubiquitin_not_predicted_idp(self):
        """Ubiquitin should NOT be predicted as IDP."""
        assert not predict_idp(UBIQUITIN), (
            "Ubiquitin should NOT be predicted as IDP"
        )

    def test_short_sequence_not_idp(self):
        """Sequences shorter than 10 residues should not be predicted as IDP."""
        assert not predict_idp("KEKEKEKEK")


# ────────────────────────────────────────────────────────────
# Test: Scale auto-selection
# ────────────────────────────────────────────────────────────

class TestScaleAutoSelection:
    """Test that the Urry scale is auto-selected for IDP sequences."""

    def test_auto_selects_urry_for_alpha_synuclein(self):
        """Auto mode should select Urry scale for alpha-synuclein."""
        _, scale_name = select_hydropathy_scale(ALPHA_SYNUCLEIN, "auto")
        assert scale_name == "urry", (
            f"Expected Urry scale for alpha-synuclein, got {scale_name}"
        )

    def test_auto_selects_urry_for_tau(self):
        """Auto mode should select Urry scale for tau K18."""
        _, scale_name = select_hydropathy_scale(TAU_K18, "auto")
        assert scale_name == "urry", (
            f"Expected Urry scale for tau K18, got {scale_name}"
        )

    def test_auto_selects_urry_for_huntingtin(self):
        """Auto mode should select Urry scale for huntingtin."""
        _, scale_name = select_hydropathy_scale(HUNTINGTIN_EX1, "auto")
        assert scale_name == "urry", (
            f"Expected Urry scale for huntingtin, got {scale_name}"
        )

    def test_auto_selects_urry_for_p53_tad(self):
        """Auto mode should select Urry scale for p53 TAD."""
        _, scale_name = select_hydropathy_scale(P53_TAD, "auto")
        assert scale_name == "urry", (
            f"Expected Urry scale for p53 TAD, got {scale_name}"
        )

    def test_auto_selects_ww_for_thioredoxin(self):
        """Auto mode should select Wimley-White scale for thioredoxin."""
        _, scale_name = select_hydropathy_scale(THIOREDOXIN, "auto")
        assert scale_name == "wimley_white", (
            f"Expected Wimley-White scale for thioredoxin, got {scale_name}"
        )

    def test_auto_selects_ww_for_ubiquitin(self):
        """Auto mode should select Wimley-White scale for ubiquitin."""
        _, scale_name = select_hydropathy_scale(UBIQUITIN, "auto")
        assert scale_name == "wimley_white", (
            f"Expected Wimley-White scale for ubiquitin, got {scale_name}"
        )

    def test_explicit_urry_selection(self):
        """Explicit 'urry' scale should always return Urry."""
        hydropathy, scale_name = select_hydropathy_scale(THIOREDOXIN, "urry")
        assert scale_name == "urry"
        assert hydropathy is URRY_HYDROPATHY

    def test_explicit_ww_selection(self):
        """Explicit 'wimley_white' scale should always return Wimley-White."""
        hydropathy, scale_name = select_hydropathy_scale(ALPHA_SYNUCLEIN, "wimley_white")
        assert scale_name == "wimley_white"
        assert hydropathy is CAMSOL_HYDROPATHY

    def test_enum_selection(self):
        """HydrophobicityScale enum should be accepted."""
        _, scale_name = select_hydropathy_scale(ALPHA_SYNUCLEIN, HydrophobicityScale.URRY)
        assert scale_name == "urry"
        _, scale_name = select_hydropathy_scale(ALPHA_SYNUCLEIN, HydrophobicityScale.WIMLEY_WHITE)
        assert scale_name == "wimley_white"

    def test_invalid_scale_raises(self):
        """Invalid scale name should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown hydrophobicity scale"):
            select_hydropathy_scale(ALPHA_SYNUCLEIN, "invalid_scale")


# ────────────────────────────────────────────────────────────
# Test: CamSolResult includes scale used
# ────────────────────────────────────────────────────────────

class TestCamSolResultScaleField:
    """Test that CamSolResult records which scale was used."""

    def test_auto_result_has_scale_field(self):
        """Result from auto mode should have hydrophobicity_scale_used."""
        clear_cache()
        result = compute_intrinsic_solubility(ALPHA_SYNUCLEIN)
        assert hasattr(result, "hydrophobicity_scale_used")
        assert result.hydrophobicity_scale_used in ("wimley_white", "urry")

    def test_idp_result_uses_urry(self):
        """IDP result should show Urry scale was used."""
        clear_cache()
        result = compute_intrinsic_solubility(ALPHA_SYNUCLEIN)
        assert result.hydrophobicity_scale_used == "urry", (
            f"Alpha-synuclein should use Urry scale, got {result.hydrophobicity_scale_used}"
        )

    def test_globular_result_uses_ww(self):
        """Globular protein result should show Wimley-White scale was used."""
        clear_cache()
        result = compute_intrinsic_solubility(THIOREDOXIN)
        assert result.hydrophobicity_scale_used == "wimley_white", (
            f"Thioredoxin should use Wimley-White scale, got {result.hydrophobicity_scale_used}"
        )

    def test_explicit_urry_overrides_auto(self):
        """Explicit 'urry' should override auto-selection."""
        clear_cache()
        result = compute_intrinsic_solubility(THIOREDOXIN, hydrophobicity_scale="urry")
        assert result.hydrophobicity_scale_used == "urry"

    def test_explicit_ww_overrides_auto(self):
        """Explicit 'wimley_white' should override auto-selection."""
        clear_cache()
        result = compute_intrinsic_solubility(ALPHA_SYNUCLEIN, hydrophobicity_scale="wimley_white")
        assert result.hydrophobicity_scale_used == "wimley_white"


# ────────────────────────────────────────────────────────────
# Test: IDP sensitivity improvement
# ────────────────────────────────────────────────────────────

class TestIDPSensitivity:
    """Test that IDP sensitivity improves with Urry scale auto-selection.

    The 9 IDP test proteins are:
    1. Alpha-synuclein
    2. Tau K18 fragment
    3. p53 Transactivation Domain
    4. Huntingtin Exon-1
    5. Amyloid-beta 42
    6. Human Prion Protein
    7. Human IAPP
    8. Prothymosin alpha
    9. c-Myc Transactivation Domain
    """

    IDP_PROTEINS = [
        ("Alpha-synuclein", ALPHA_SYNUCLEIN),
        ("Tau K18", TAU_K18),
        ("p53 TAD", P53_TAD),
        ("Huntingtin Exon-1", HUNTINGTIN_EX1),
        ("Amyloid-beta 42", AB42),
        ("Prion Protein", PRION),
        ("IAPP", IAPP),
        ("Prothymosin alpha", PROTHYMOSIN_ALPHA),
        ("c-Myc TAD", CMYC_TAD),
    ]

    @pytest.fixture(autouse=True)
    def _clear(self):
        """Clear cache before each test."""
        clear_cache()
        yield
        clear_cache()

    def test_idp_sensitivity_with_auto_scale(self):
        """IDP sensitivity with auto scale should be >= 67% (6/9)."""
        correct = 0
        total = len(self.IDP_PROTEINS)
        details = []

        for name, seq in self.IDP_PROTEINS:
            result = compute_intrinsic_solubility(seq)
            enhanced = compute_enhanced_benchmark_score(result)
            is_aggregation_prone = enhanced < 0
            if is_aggregation_prone:
                correct += 1
            details.append(
                f"{name}: score={result.intrinsic_score:.4f}, "
                f"enhanced={enhanced:.4f}, scale={result.hydrophobicity_scale_used}, "
                f"{'correct' if is_aggregation_prone else 'MISSED'}"
            )

        sensitivity = correct / total
        assert correct >= 6, (
            f"IDP sensitivity {sensitivity:.1%} ({correct}/{total}) is below "
            f"67% target (need at least 6/9).\nDetails:\n" + "\n".join(details)
        )

    def test_idp_sensitivity_improvement_over_ww_only(self):
        """Auto-scale sensitivity should be strictly better than Wimley-White only."""
        auto_correct = 0
        ww_correct = 0

        for name, seq in self.IDP_PROTEINS:
            # Auto scale
            clear_cache()
            auto_result = compute_intrinsic_solubility(seq, hydrophobicity_scale="auto")
            auto_enhanced = compute_enhanced_benchmark_score(auto_result)
            if auto_enhanced < 0:
                auto_correct += 1

            # Wimley-White only
            clear_cache()
            ww_result = compute_intrinsic_solubility(seq, hydrophobicity_scale="wimley_white")
            ww_enhanced = compute_enhanced_benchmark_score(ww_result)
            if ww_enhanced < 0:
                ww_correct += 1

        assert auto_correct > ww_correct, (
            f"Auto-scale ({auto_correct}/9) should improve over "
            f"Wimley-White only ({ww_correct}/9)"
        )

    def test_alpha_synuclein_correctly_predicted(self):
        """Alpha-synuclein should be predicted as aggregation-prone with Urry scale."""
        result = compute_intrinsic_solubility(ALPHA_SYNUCLEIN)
        enhanced = compute_enhanced_benchmark_score(result)
        assert enhanced < 0, (
            f"Alpha-synuclein enhanced score {enhanced:.4f} should be negative "
            f"(aggregation-prone). Scale used: {result.hydrophobicity_scale_used}"
        )

    def test_tau_correctly_predicted(self):
        """Tau K18 should be predicted as aggregation-prone with Urry scale."""
        result = compute_intrinsic_solubility(TAU_K18)
        enhanced = compute_enhanced_benchmark_score(result)
        assert enhanced < 0, (
            f"Tau K18 enhanced score {enhanced:.4f} should be negative "
            f"(aggregation-prone). Scale used: {result.hydrophobicity_scale_used}"
        )

    def test_p53_tad_correctly_predicted(self):
        """p53 TAD should be predicted as aggregation-prone with Urry scale."""
        result = compute_intrinsic_solubility(P53_TAD)
        enhanced = compute_enhanced_benchmark_score(result)
        assert enhanced < 0, (
            f"p53 TAD enhanced score {enhanced:.4f} should be negative "
            f"(aggregation-prone). Scale used: {result.hydrophobicity_scale_used}"
        )

    def test_huntingtin_urry_scale_selected(self):
        """Huntingtin exon-1 should auto-select Urry scale (high disorder content).

        Note: Huntingtin's aggregation is driven by polyglutamine (polyQ) expansion,
        which uses a 'polar zipper' mechanism rather than hydrophobic interactions.
        Neither the Wimley-White nor the Urry scale captures this mechanism because
        glutamine is classified as polar/soluble in both scales. The Urry scale is
        still selected (correctly identifying the sequence as IDP-like), but the
        overall solubility score may remain positive.
        """
        result = compute_intrinsic_solubility(HUNTINGTIN_EX1)
        assert result.hydrophobicity_scale_used == "urry", (
            f"Huntingtin should auto-select Urry scale, got {result.hydrophobicity_scale_used}"
        )
        # Verify that Urry scale gives a lower score than WW for huntingtin
        clear_cache()
        ww_result = compute_intrinsic_solubility(HUNTINGTIN_EX1, hydrophobicity_scale="wimley_white")
        assert result.intrinsic_score < ww_result.intrinsic_score, (
            f"Urry score ({result.intrinsic_score:.4f}) should be lower than "
            f"WW score ({ww_result.intrinsic_score:.4f}) for huntingtin"
        )


# ────────────────────────────────────────────────────────────
# Test: Specificity preserved for soluble proteins
# ────────────────────────────────────────────────────────────

class TestSpecificityPreserved:
    """Test that the Urry scale does not degrade specificity for soluble proteins."""

    SOLUBLE_PROTEINS = [
        ("Thioredoxin", THIOREDOXIN),
        ("Ubiquitin", UBIQUITIN),
    ]

    @pytest.fixture(autouse=True)
    def _clear(self):
        clear_cache()
        yield
        clear_cache()

    def test_soluble_proteins_still_positive(self):
        """Known soluble proteins should still have positive enhanced scores with auto scale."""
        for name, seq in self.SOLUBLE_PROTEINS:
            result = compute_intrinsic_solubility(seq)
            enhanced = compute_enhanced_benchmark_score(result)
            assert enhanced > 0, (
                f"{name} enhanced score {enhanced:.4f} should be positive "
                f"(soluble). Scale used: {result.hydrophobicity_scale_used}"
            )

    def test_benchmark_high_solubility_specificity(self):
        """All benchmark high-solubility proteins should still have positive enhanced scores."""
        high_entries = [
            (name, seq)
            for name, _uid, seq, known, _ref in BENCHMARK_DATASET
            if known == "high"
        ]
        failures = []
        for name, seq in high_entries:
            clear_cache()
            result = compute_intrinsic_solubility(seq)
            enhanced = compute_enhanced_benchmark_score(result)
            if enhanced <= 0:
                failures.append(f"{name}: enhanced={enhanced:.4f}")

        assert len(failures) == 0, (
            f"{len(failures)} benchmark high-solubility proteins have non-positive "
            f"enhanced scores:\n" + "\n".join(failures)
        )


# ────────────────────────────────────────────────────────────
# Test: Urry scale gives lower scores for IDPs than Wimley-White
# ────────────────────────────────────────────────────────────

class TestUrryScaleEffect:
    """Test that Urry scale produces lower (more aggregation-prone) scores for IDPs."""

    @pytest.fixture(autouse=True)
    def _clear(self):
        clear_cache()
        yield
        clear_cache()

    def test_urry_lower_than_ww_for_alpha_synuclein(self):
        """Urry scale should give a lower score for alpha-synuclein than Wimley-White."""
        clear_cache()
        ww_result = compute_intrinsic_solubility(ALPHA_SYNUCLEIN, hydrophobicity_scale="wimley_white")
        clear_cache()
        urry_result = compute_intrinsic_solubility(ALPHA_SYNUCLEIN, hydrophobicity_scale="urry")
        assert urry_result.intrinsic_score < ww_result.intrinsic_score, (
            f"Urry score ({urry_result.intrinsic_score:.4f}) should be lower than "
            f"Wimley-White ({ww_result.intrinsic_score:.4f}) for alpha-synuclein"
        )

    def test_urry_lower_than_ww_for_tau(self):
        """Urry scale should give a lower score for tau K18 than Wimley-White."""
        clear_cache()
        ww_result = compute_intrinsic_solubility(TAU_K18, hydrophobicity_scale="wimley_white")
        clear_cache()
        urry_result = compute_intrinsic_solubility(TAU_K18, hydrophobicity_scale="urry")
        assert urry_result.intrinsic_score < ww_result.intrinsic_score, (
            f"Urry score ({urry_result.intrinsic_score:.4f}) should be lower than "
            f"Wimley-White ({ww_result.intrinsic_score:.4f}) for tau K18"
        )

    def test_urry_lower_than_ww_for_huntingtin(self):
        """Urry scale should give a lower score for huntingtin than Wimley-White."""
        clear_cache()
        ww_result = compute_intrinsic_solubility(HUNTINGTIN_EX1, hydrophobicity_scale="wimley_white")
        clear_cache()
        urry_result = compute_intrinsic_solubility(HUNTINGTIN_EX1, hydrophobicity_scale="urry")
        assert urry_result.intrinsic_score < ww_result.intrinsic_score, (
            f"Urry score ({urry_result.intrinsic_score:.4f}) should be lower than "
            f"Wimley-White ({ww_result.intrinsic_score:.4f}) for huntingtin"
        )
