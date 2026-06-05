"""
BioCompiler v12 Integration Smoke Tests
=========================================

End-to-end tests verifying the full v12 pipeline:
  biosecurity screen → optimize → verify translation → strict mode check → export

Tests cover all 5 organisms (e_coli, human, yeast, mouse, cho) and
3 reference proteins (insulin, GFP, HBB).

CAI expected ranges are based on published values and the v10+ CAI
table unification (organism-specific CODON_ADAPTIVENESS_TABLES).
"""

import pytest
import warnings

import biocompiler
from biocompiler import (
    __version__,
    optimize_sequence,
    compute_cai,
    translate,
    export_fasta,
    export_genbank,
    BiosecurityError,
    OptimizationConstraintError,
    TranslationVerificationError,
    UnsupportedOrganismError,
)
from biocompiler.biosecurity import (
    screen_hazardous_sequence,
    check_biosecurity_before_optimize,
)
from biocompiler.protein_verification import verify_translation, verify_and_raise
from biocompiler.organism_config import (
    get_organism_config,
    is_eukaryotic_organism,
    get_constraint_profile,
    ORGANISM_CONFIGS,
)
from biocompiler.exceptions import (
    BioCompilerError,
    EngineError,
    ESMFoldError,
    FoldXError,
    CamSolError,
    ImmunogenicityError,
    BiosecurityError,
    TranslationVerificationError,
    OptimizationConstraintError,
    InvalidSequenceError,
    CertificateGenerationError,
    CertificateVerificationError,
    UnknownPredicateError,
    OptimizationError,
    UnsupportedOrganismError,
    InvalidProteinError,
    FileFormatError,
    SplicingError,
    MutagenesisError,
)


# ──────────────────────────────────────────────────────────────────────
# Reference protein sequences
# ──────────────────────────────────────────────────────────────────────

INSULIN = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKTREEDFPQVEKLEGGGPRGADDVRVLMEDCGLNVEFLPSAFFLSRDLRAEVDGPKVRDSFVKNIYIVDSQVTLPTEEPQVDKLQGQVLNLPVDNHNMTSIFSVQKRLGKLNLC"

GFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

HBB = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"


# ──────────────────────────────────────────────────────────────────────
# Organism keys accepted by optimize_sequence
# ──────────────────────────────────────────────────────────────────────

ORGANISMS = ["e_coli", "human", "yeast", "mouse", "cho"]

# Canonical organism names for compute_cai
ORGANISM_CANONICAL = {
    "e_coli": "Escherichia_coli",
    "human": "Homo_sapiens",
    "yeast": "Saccharomyces_cerevisiae",
    "mouse": "Mus_musculus",
    "cho": "CHO_K1",
}

# Expected CAI ranges per organism/protein.
# These are conservative ranges; the optimizer should achieve at least
# the lower bound. Upper bound is not enforced strictly — just a sanity
# check that CAI doesn't exceed 1.0.
EXPECTED_CAI_RANGES = {
    "e_coli": {
        "insulin": (0.50, 1.0),
        "gfp":     (0.70, 1.0),
        "hbb":     (0.50, 1.0),
    },
    "human": {
        "insulin": (0.50, 1.0),
        "gfp":     (0.50, 1.0),
        "hbb":     (0.60, 1.0),
    },
    "yeast": {
        "insulin": (0.40, 1.0),
        "gfp":     (0.50, 1.0),
        "hbb":     (0.40, 1.0),
    },
    "mouse": {
        "insulin": (0.40, 1.0),
        "gfp":     (0.50, 1.0),
        "hbb":     (0.40, 1.0),
    },
    "cho": {
        "insulin": (0.40, 1.0),
        "gfp":     (0.50, 1.0),
        "hbb":     (0.40, 1.0),
    },
}

PROTEINS = {
    "insulin": INSULIN,
    "gfp": GFP,
    "hbb": HBB,
}


# ═══════════════════════════════════════════════════════════════════════
# Test: Package version
# ═══════════════════════════════════════════════════════════════════════

class TestVersion:
    """Verify the package version is 12.0.0."""

    def test_version_is_12_0_0(self):
        assert __version__ == "12.0.0", f"Expected version 12.0.0, got {__version__}"

    def test_version_string_format(self):
        parts = __version__.split(".")
        assert len(parts) == 3, f"Version should be semver, got {__version__}"
        assert all(p.isdigit() for p in parts), f"Version parts should be numeric, got {__version__}"


# ═══════════════════════════════════════════════════════════════════════
# Test: All new modules import correctly
# ═══════════════════════════════════════════════════════════════════════

class TestModuleImports:
    """Verify that all new v12 modules import correctly."""

    def test_biosecurity_module(self):
        from biocompiler import biosecurity
        assert hasattr(biosecurity, "screen_hazardous_sequence")
        assert hasattr(biosecurity, "check_biosecurity_before_optimize")
        assert hasattr(biosecurity, "BiosecurityReport")
        assert hasattr(biosecurity, "HazardMatch")

    def test_protein_verification_module(self):
        from biocompiler import protein_verification
        assert hasattr(protein_verification, "verify_translation")
        assert hasattr(protein_verification, "verify_and_raise")
        assert hasattr(protein_verification, "VerificationResult")
        assert hasattr(protein_verification, "PositionMismatch")

    def test_organism_config_module(self):
        from biocompiler import organism_config
        assert hasattr(organism_config, "OrganismConfig")
        assert hasattr(organism_config, "ORGANISM_CONFIGS")
        assert hasattr(organism_config, "get_organism_config")
        assert hasattr(organism_config, "is_eukaryotic_organism")

    def test_exceptions_module(self):
        from biocompiler import exceptions
        assert hasattr(exceptions, "BioCompilerError")
        assert hasattr(exceptions, "EngineError")
        assert hasattr(exceptions, "ESMFoldError")
        assert hasattr(exceptions, "FoldXError")
        assert hasattr(exceptions, "CamSolError")
        assert hasattr(exceptions, "ImmunogenicityError")
        assert hasattr(exceptions, "BiosecurityError")
        assert hasattr(exceptions, "TranslationVerificationError")
        assert hasattr(exceptions, "OptimizationConstraintError")
        assert hasattr(exceptions, "InvalidSequenceError")
        assert hasattr(exceptions, "CertificateGenerationError")
        assert hasattr(exceptions, "CertificateVerificationError")
        assert hasattr(exceptions, "UnknownPredicateError")
        assert hasattr(exceptions, "OptimizationError")
        assert hasattr(exceptions, "UnsupportedOrganismError")
        assert hasattr(exceptions, "InvalidProteinError")
        assert hasattr(exceptions, "FileFormatError")
        assert hasattr(exceptions, "SplicingError")
        assert hasattr(exceptions, "MutagenesisError")

    def test_organisms_subpackage(self):
        from biocompiler import organisms
        assert hasattr(organisms, "OrganismDatabase")
        assert hasattr(organisms, "CODON_ADAPTIVENESS_TABLES")
        assert hasattr(organisms, "SUPPORTED_ORGANISMS")

    def test_solver_subpackage(self):
        from biocompiler import solver
        assert hasattr(solver, "CSPSolver")

    def test_provenance_module(self):
        from biocompiler import provenance
        assert hasattr(provenance, "DecisionRecord")
        assert hasattr(provenance, "ProvenanceTracker")

    def test_incremental_module(self):
        from biocompiler import incremental
        assert hasattr(incremental, "IncrementalSequenceState")

    def test_hybrid_optimizer_module(self):
        from biocompiler import hybrid_optimizer
        assert hasattr(hybrid_optimizer, "HybridOptimizer")

    def test_objectives_module(self):
        from biocompiler import objectives
        assert hasattr(objectives, "cai_objective")

    def test_export_module(self):
        from biocompiler import export
        assert hasattr(export, "export_fasta")
        assert hasattr(export, "export_genbank")

    def test_translation_module(self):
        from biocompiler import translation
        assert hasattr(translation, "translate")
        assert hasattr(translation, "compute_cai")


# ═══════════════════════════════════════════════════════════════════════
# Test: Exception hierarchy
# ═══════════════════════════════════════════════════════════════════════

class TestExceptionHierarchy:
    """Verify exception inheritance and all new exceptions exist."""

    def test_engine_errors_inherit_from_engine_error(self):
        assert issubclass(ESMFoldError, EngineError)
        assert issubclass(FoldXError, EngineError)
        assert issubclass(CamSolError, EngineError)
        assert issubclass(ImmunogenicityError, EngineError)

    def test_all_errors_inherit_from_biocompiler_error(self):
        for exc_class in [
            BioCompilerError, EngineError, ESMFoldError, FoldXError,
            CamSolError, ImmunogenicityError, BiosecurityError,
            TranslationVerificationError, OptimizationConstraintError,
            InvalidSequenceError, CertificateGenerationError,
            CertificateVerificationError, UnknownPredicateError,
            OptimizationError, UnsupportedOrganismError,
            InvalidProteinError, FileFormatError, SplicingError,
            MutagenesisError,
        ]:
            assert issubclass(exc_class, BioCompilerError), (
                f"{exc_class.__name__} should inherit from BioCompilerError"
            )

    def test_biosecurity_error_instantiation(self):
        err = BiosecurityError("test hazard", risk_level="high")
        assert err.risk_level == "high"
        assert "test hazard" in str(err)

    def test_translation_verification_error_instantiation(self):
        err = TranslationVerificationError(
            reason="mismatch found",
            mismatches=[],
            translated_protein="M",
            expected_protein="V",
        )
        assert "mismatch found" in str(err)

    def test_optimization_constraint_error_instantiation(self):
        err = OptimizationConstraintError(
            failed_predicates=["no_cryptic_splice", "gc_in_range"],
            partial_result=None,
        )
        assert err.failed_predicates == ["no_cryptic_splice", "gc_in_range"]
        assert "strict mode" in str(err)


# ═══════════════════════════════════════════════════════════════════════
# Test: Biosecurity screening
# ═══════════════════════════════════════════════════════════════════════

class TestBiosecurityScreening:
    """Verify biosecurity screening for safe proteins."""

    @pytest.mark.parametrize("protein_name,protein", [
        ("insulin", INSULIN),
        ("gfp", GFP),
        ("hbb", HBB),
    ])
    def test_safe_proteins_pass_screening(self, protein_name, protein):
        report = screen_hazardous_sequence(protein)
        # Our reference proteins should not be flagged as hazardous
        # (they are therapeutic/research proteins, not toxins)
        # Note: Some may trigger low-risk oncogene hits (e.g. VEGF-like motifs)
        # but should not be high/critical
        assert report.risk_level in ("none", "low", "medium"), (
            f"{protein_name} should not be high/critical risk, "
            f"got risk_level={report.risk_level}, "
            f"matches={[m.name for m in report.matches]}"
        )

    def test_biosecurity_gate_does_not_block_safe_proteins(self):
        """check_biosecurity_before_optimize should not raise for safe proteins."""
        for name, protein in PROTEINS.items():
            # Should not raise BiosecurityError
            report = check_biosecurity_before_optimize(protein, organism="Homo_sapiens")
            assert report is not None

    def test_hazardous_protein_raises_biosecurity_error(self):
        """A protein matching a select agent toxin should raise BiosecurityError."""
        # Ricin A-chain active site motif
        hazardous_protein = "NIRVGLPIIS" * 5  # repeat to ensure detection
        with pytest.raises(BiosecurityError):
            check_biosecurity_before_optimize(hazardous_protein, organism="Homo_sapiens")


# ═══════════════════════════════════════════════════════════════════════
# Test: Full pipeline per organism × protein
# ═══════════════════════════════════════════════════════════════════════

class TestFullPipeline:
    """Test the full v12 pipeline for each organism × protein combination."""

    @pytest.mark.parametrize("organism", ORGANISMS)
    @pytest.mark.parametrize("protein_name", ["insulin", "gfp", "hbb"])
    def test_full_pipeline_no_crash(self, organism, protein_name):
        """Verify the full pipeline does not crash for any combination."""
        protein = PROTEINS[protein_name]

        # Step 1: Biosecurity screen
        report = screen_hazardous_sequence(protein)

        # Step 2: Optimize (with strict_mode=False to avoid raising on
        # marginal cases — we verify strict mode separately)
        result = optimize_sequence(
            protein,
            organism=organism,
            strict_mode=False,
        )

        # Step 3: Verify translation
        verification = verify_translation(result.sequence, protein)
        assert verification.is_valid, (
            f"Translation verification failed for {protein_name}/{organism}: "
            f"mismatches={len(verification.mismatches)}, "
            f"premature_stop={verification.has_premature_stop}, "
            f"length_correct={verification.length_correct}"
        )

        # Step 4: Check CAI is in expected range
        cai_lo, cai_hi = EXPECTED_CAI_RANGES[organism][protein_name]
        assert cai_lo <= result.cai <= cai_hi, (
            f"CAI out of range for {protein_name}/{organism}: "
            f"got {result.cai:.4f}, expected [{cai_lo}, {cai_hi}]"
        )

        # Step 5: Export should work
        fasta = export_fasta(result.sequence, organism=organism, protein=protein)
        # FASTA may start with comment lines (;) before the > header
        assert ">" in fasta
        assert len(fasta) > len(result.sequence)

    @pytest.mark.parametrize("organism", ORGANISMS)
    @pytest.mark.parametrize("protein_name", ["insulin", "gfp", "hbb"])
    def test_cai_value_in_expected_range(self, organism, protein_name):
        """Verify CAI values are in expected ranges after optimization."""
        protein = PROTEINS[protein_name]
        result = optimize_sequence(
            protein,
            organism=organism,
            strict_mode=False,
        )
        cai_lo, cai_hi = EXPECTED_CAI_RANGES[organism][protein_name]
        assert cai_lo <= result.cai <= cai_hi, (
            f"CAI={result.cai:.4f} out of [{cai_lo}, {cai_hi}] "
            f"for {protein_name}/{organism}"
        )

    @pytest.mark.parametrize("organism", ORGANISMS)
    @pytest.mark.parametrize("protein_name", ["insulin", "gfp", "hbb"])
    def test_translation_verification_passes(self, organism, protein_name):
        """Verify that optimized DNA translates back to the input protein."""
        protein = PROTEINS[protein_name]
        result = optimize_sequence(
            protein,
            organism=organism,
            strict_mode=False,
        )
        verification = verify_translation(result.sequence, protein)
        assert verification.is_valid, (
            f"Translation verification failed: {protein_name}/{organism}, "
            f"mismatches={verification.mismatches}"
        )

    @pytest.mark.parametrize("organism", ORGANISMS)
    @pytest.mark.parametrize("protein_name", ["insulin", "gfp", "hbb"])
    def test_verify_and_raise_does_not_raise(self, organism, protein_name):
        """verify_and_raise should succeed for properly optimized sequences."""
        protein = PROTEINS[protein_name]
        result = optimize_sequence(
            protein,
            organism=organism,
            strict_mode=False,
        )
        # Should not raise TranslationVerificationError
        verification = verify_and_raise(result.sequence, protein, organism=organism)
        assert verification.is_valid


# ═══════════════════════════════════════════════════════════════════════
# Test: Strict mode
# ═══════════════════════════════════════════════════════════════════════

class TestStrictMode:
    """Verify strict mode behavior."""

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_strict_mode_returns_valid_result_or_raises(self, organism):
        """In strict mode, optimize_sequence either returns a valid result
        or raises OptimizationConstraintError."""
        protein = GFP  # Use GFP as a well-behaved test case
        try:
            result = optimize_sequence(
                protein,
                organism=organism,
                strict_mode=True,
            )
            # If it returns, the result should have no failed predicates
            # (or at most a few marginal ones)
            assert result.cai > 0.0
            assert 0.0 <= result.gc_content <= 1.0
        except OptimizationConstraintError as e:
            # If it raises, the error should list the failed predicates
            assert len(e.failed_predicates) > 0
            assert e.partial_result is not None


# ═══════════════════════════════════════════════════════════════════════
# Test: Organism configuration
# ═══════════════════════════════════════════════════════════════════════

class TestOrganismConfiguration:
    """Verify organism configuration for all 5 organisms."""

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_organism_config_exists(self, organism):
        config = get_organism_config(organism)
        assert config is not None
        assert config.name  # Should have a human-readable name

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_organism_domain_classification(self, organism):
        if organism == "e_coli":
            assert not is_eukaryotic_organism(organism), (
                "E. coli should be classified as prokaryote"
            )
        else:
            assert is_eukaryotic_organism(organism), (
                f"{organism} should be classified as eukaryote"
            )

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_organism_constraint_profile(self, organism):
        profile = get_constraint_profile(organism)
        assert "cai" in profile
        assert "gc_content" in profile
        # E. coli should not have splice avoidance or CpG avoidance
        if organism == "e_coli":
            assert profile.get("splice_avoidance") is False
            assert profile.get("cpg_avoidance") is False
        # Eukaryotes should have splice avoidance and CpG avoidance
        elif organism in ("human", "mouse", "cho"):
            assert profile.get("splice_avoidance") is True
            assert profile.get("cpg_avoidance") is True

    def test_all_five_organisms_in_configs(self):
        """Verify all 5 target organisms are in ORGANISM_CONFIGS."""
        for org in ORGANISMS:
            config = get_organism_config(org)
            assert config is not None, f"Missing config for organism: {org}"


# ═══════════════════════════════════════════════════════════════════════
# Test: Export formats
# ═══════════════════════════════════════════════════════════════════════

class TestExportFormats:
    """Verify export functionality works for optimized sequences."""

    @pytest.mark.parametrize("organism", ORGANISMS[:2])  # e_coli and human only for speed
    def test_fasta_export(self, organism):
        protein = GFP
        result = optimize_sequence(protein, organism=organism, strict_mode=False)
        fasta = export_fasta(
            result.sequence,
            identifier="test_gfp",
            organism=organism,
            protein=protein,
            cai=result.cai,
        )
        assert ">test_gfp" in fasta
        assert "organism=" in fasta
        assert "CAI=" in fasta
        assert "GC=" in fasta

    @pytest.mark.parametrize("organism", ORGANISMS[:2])  # e_coli and human only for speed
    def test_genbank_export(self, organism):
        protein = GFP
        result = optimize_sequence(protein, organism=organism, strict_mode=False)
        gb = export_genbank(
            result.sequence,
            locus_name="TESTGFP",
            organism=organism,
            gene_name="GFP",
            protein=protein,
            cai=result.cai,
        )
        assert "LOCUS" in gb
        assert "ORIGIN" in gb
        assert "//" in gb
        assert "GFP" in gb


# ═══════════════════════════════════════════════════════════════════════
# Test: CAI computation consistency
# ═══════════════════════════════════════════════════════════════════════

class TestCAIConsistency:
    """Verify that compute_cai agrees with the optimizer's CAI value."""

    @pytest.mark.parametrize("organism", ORGANISMS[:2])
    @pytest.mark.parametrize("protein_name", ["gfp"])
    def test_cai_matches_optimizer(self, organism, protein_name):
        protein = PROTEINS[protein_name]
        result = optimize_sequence(protein, organism=organism, strict_mode=False)
        canonical = ORGANISM_CANONICAL[organism]
        cai_from_compute = compute_cai(result.sequence, organism=canonical)
        # Allow small rounding difference
        assert abs(result.cai - cai_from_compute) < 0.01, (
            f"CAI mismatch: optimizer={result.cai:.4f}, "
            f"compute_cai={cai_from_compute:.4f}"
        )


# ═══════════════════════════════════════════════════════════════════════
# Test: Deterministic behavior
# ═══════════════════════════════════════════════════════════════════════

class TestDeterminism:
    """Verify that the pipeline is deterministic."""

    def test_same_input_same_output(self):
        """Running the same optimization twice should produce identical results."""
        for organism in ["e_coli", "human"]:
            result1 = optimize_sequence(GFP, organism=organism, strict_mode=False)
            result2 = optimize_sequence(GFP, organism=organism, strict_mode=False)
            assert result1.sequence == result2.sequence, (
                f"Non-deterministic output for {organism}"
            )
            assert result1.cai == result2.cai


# ═══════════════════════════════════════════════════════════════════════
# Test: Edge cases
# ═══════════════════════════════════════════════════════════════════════

class TestEdgeCases:
    """Verify edge case handling."""

    def test_unsupported_organism_raises(self):
        with pytest.raises(UnsupportedOrganismError):
            optimize_sequence(GFP, organism="totally_fake_organism_xyz", strict_mode=False)

    def test_empty_protein_raises(self):
        with pytest.raises((InvalidProteinError, ValueError, AssertionError)):
            optimize_sequence("", organism="e_coli", strict_mode=False)

    def test_short_protein_works(self):
        """A very short protein (e.g., 10 aa) should still optimize."""
        short_protein = "MSKGEELFTG"
        result = optimize_sequence(short_protein, organism="e_coli", strict_mode=False)
        assert result.sequence  # Non-empty result
        assert len(result.sequence) == len(short_protein) * 3
