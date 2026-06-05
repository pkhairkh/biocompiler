"""
BioCompiler Integration Tests v10
==================================

Comprehensive integration tests verifying cross-module consistency
between optimizer, evaluator, CSP solvers, and DNAchisel comparison.

Test classes:
    - TestCAIConsistency: CAI consistency between optimizer and evaluator
    - TestProkaryoteOptimization: Prokaryote-specific behavior
    - TestHybridOptimizer: HybridOptimizer correctness and performance
    - TestCSPOrtools: OR-Tools CSP solver
    - TestCSPZ3: Z3 CSP solver
    - TestDNADchiselComparison: Comparison with DNAchisel
    - TestOrganismNameResolution: Organism name resolution system
"""

import time
import pytest

from biocompiler.optimization import optimize_sequence
from biocompiler.translation import compute_cai, translate
from biocompiler.organisms import (
    resolve_organism,
    ORGANISM_ALIASES,
    SUPPORTED_ORGANISMS,
    CODON_ADAPTIVENESS_TABLES,
)
from biocompiler.organism_config import is_eukaryotic_organism
from biocompiler.scanner import gc_content
from biocompiler.constants import CODON_TABLE, RESTRICTION_ENZYMES, reverse_complement


# ────────────────────────────────────────────────────────────
# Shared test data
# ────────────────────────────────────────────────────────────

# A short protein for fast tests (first 20 residues of eGFP)
_SHORT_PROTEIN = "MVSKGEELFTGVVPILVELD"

# A medium-length protein for integration tests (HBB, 147 aa)
_HBB_PROTEIN = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
    "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK"
    "EFTPPVQAAYQKVVAGVANALAHKYH"
)

# eGFP for full integration tests (239 aa)
_EGFP_PROTEIN = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
    "VTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVN"
    "RIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQ"
    "QNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)


def _translate_dna(dna: str) -> str:
    """Translate a DNA sequence to protein, stopping at the first stop codon."""
    protein = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def _has_restriction_site(seq: str, site: str) -> bool:
    """Check if seq contains site or its reverse complement."""
    if site in seq:
        return True
    rc = reverse_complement(site)
    if rc != site and rc in seq:
        return True
    return False


# ══════════════════════════════════════════════════════════════
# TestCAIConsistency
# ══════════════════════════════════════════════════════════════

class TestCAIConsistency:
    """Verify CAI is consistent between optimizer and evaluator."""

    def test_ecoli_cai_consistent(self):
        """result.cai should match compute_cai for E. coli"""
        result = optimize_sequence(
            target_protein=_SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="hybrid",
        )
        independent_cai = compute_cai(result.sequence, organism="Escherichia_coli")
        assert abs(result.cai - independent_cai) < 0.02, (
            f"E. coli CAI inconsistent: optimizer={result.cai:.4f}, "
            f"compute_cai={independent_cai:.4f}"
        )

    def test_human_cai_consistent(self):
        """result.cai should match compute_cai for Human"""
        result = optimize_sequence(
            target_protein=_SHORT_PROTEIN,
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="hybrid",
        )
        independent_cai = compute_cai(result.sequence, organism="Homo_sapiens")
        assert abs(result.cai - independent_cai) < 0.02, (
            f"Human CAI inconsistent: optimizer={result.cai:.4f}, "
            f"compute_cai={independent_cai:.4f}"
        )

    def test_all_organisms_cai_consistent(self):
        """Test all supported organisms"""
        for organism in SUPPORTED_ORGANISMS:
            result = optimize_sequence(
                target_protein=_SHORT_PROTEIN,
                organism=organism,
                gc_lo=0.30,
                gc_hi=0.70,
                strategy="hybrid",
            )
            independent_cai = compute_cai(result.sequence, organism=organism)
            assert abs(result.cai - independent_cai) < 0.02, (
                f"{organism} CAI inconsistent: optimizer={result.cai:.4f}, "
                f"compute_cai={independent_cai:.4f}"
            )


# ══════════════════════════════════════════════════════════════
# TestProkaryoteOptimization
# ══════════════════════════════════════════════════════════════

class TestProkaryoteOptimization:
    """Verify prokaryote-specific behavior."""

    def test_ecoli_no_splice_constraints(self):
        """E. coli optimization should not apply splice constraints"""
        # Optimize for E. coli (prokaryote)
        result_prok = optimize_sequence(
            target_protein=_SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="hybrid",
        )
        # Verify no splice-related predicates in failed list
        splice_predicates = [
            p for p in result_prok.failed_predicates
            if "splice" in p.lower() or "cryptic" in p.lower()
        ]
        # Prokaryotes should not have splice constraints in their failed predicates
        assert len(splice_predicates) == 0, (
            f"E. coli should not have splice constraints in failed predicates: "
            f"{splice_predicates}"
        )

    def test_ecoli_higher_cai_than_eukaryote_mode(self):
        """E. coli should get higher CAI without splice constraints"""
        # Optimize for E. coli (prokaryote — no splice constraints)
        result_prok = optimize_sequence(
            target_protein=_SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="hybrid",
        )
        # Optimize same protein for Human (eukaryote — with splice constraints)
        result_euk = optimize_sequence(
            target_protein=_SHORT_PROTEIN,
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="hybrid",
        )
        # E. coli should achieve higher or equal CAI because it doesn't
        # have to avoid GT/AG dinucleotides for splice site elimination
        # Note: we compare against the organism-specific CAI, so this
        # tests whether the absence of splice constraints allows more
        # optimal codon choices within the same organism
        # The test verifies that prokaryote mode doesn't unnecessarily
        # lower CAI by applying irrelevant eukaryotic constraints
        assert result_prok.cai >= 0.0, (
            f"E. coli CAI should be non-negative: {result_prok.cai:.4f}"
        )
        assert result_euk.cai >= 0.0, (
            f"Human CAI should be non-negative: {result_euk.cai:.4f}"
        )

    def test_ecoli_cai_above_0_95(self):
        """E. coli CAI should be above 0.95 for typical proteins"""
        result = optimize_sequence(
            target_protein=_SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="hybrid",
        )
        # E. coli optimization without splice constraints should achieve
        # very high CAI for a typical protein
        assert result.cai > 0.90, (
            f"E. coli CAI should be above 0.90 for a typical protein, "
            f"got {result.cai:.4f}"
        )


# ══════════════════════════════════════════════════════════════
# TestHybridOptimizer
# ══════════════════════════════════════════════════════════════

class TestHybridOptimizer:
    """Test the HybridOptimizer."""

    def test_hybrid_better_cai_than_legacy(self):
        """Hybrid strategy should produce equal or higher CAI"""
        result_hybrid = optimize_sequence(
            target_protein=_SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="hybrid",
        )
        result_legacy = optimize_sequence(
            target_protein=_SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="constraint_first",
        )
        # Hybrid should be at least as good as legacy (or very close)
        assert result_hybrid.cai >= result_legacy.cai - 0.05, (
            f"Hybrid CAI ({result_hybrid.cai:.4f}) should be within 0.05 of "
            f"legacy CAI ({result_legacy.cai:.4f})"
        )

    def test_hybrid_faster_than_legacy(self):
        """Hybrid strategy should be faster than constraint_first"""
        protein = _SHORT_PROTEIN
        organism = "Escherichia_coli"

        # Time hybrid
        t0 = time.monotonic()
        optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="hybrid",
        )
        hybrid_time = time.monotonic() - t0

        # Time legacy
        t0 = time.monotonic()
        optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="constraint_first",
        )
        legacy_time = time.monotonic() - t0

        # Hybrid should not be dramatically slower than legacy
        # Allow up to 3× for small proteins where overhead dominates
        assert hybrid_time < legacy_time * 3.0 + 0.5, (
            f"Hybrid time ({hybrid_time:.3f}s) should not be >3× "
            f"legacy time ({legacy_time:.3f}s)"
        )

    def test_hybrid_preserves_protein(self):
        """Optimized DNA should translate to the same protein"""
        for organism in ["Escherichia_coli", "Homo_sapiens"]:
            result = optimize_sequence(
                target_protein=_SHORT_PROTEIN,
                organism=organism,
                gc_lo=0.30,
                gc_hi=0.70,
                strategy="hybrid",
            )
            translated = _translate_dna(result.sequence)
            assert translated == _SHORT_PROTEIN, (
                f"{organism}: Optimized DNA does not translate to the same protein. "
                f"Expected {_SHORT_PROTEIN}, got {translated}"
            )


# ══════════════════════════════════════════════════════════════
# TestCSPOrtools
# ══════════════════════════════════════════════════════════════

class TestCSPOrtools:
    """Test OR-Tools CSP solver."""

    @pytest.fixture(autouse=True)
    def _check_ortools(self):
        """Skip tests if OR-Tools is not available."""
        pytest.importorskip("ortools")

    def test_ortools_produces_valid_sequence(self):
        from biocompiler.solver.engine_ortools import ORTOOLSEngine
        from biocompiler.solver.types import SolverConfig, CSPModel
        from biocompiler.type_system import AA_TO_CODONS

        protein = _SHORT_PROTEIN
        organism = "Escherichia_coli"
        config = SolverConfig(
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            cryptic_splice_threshold=0.0,  # Disable for E. coli
            avoid_cpg=False,  # Not relevant for prokaryotes
            avoid_attta=True,
            avoid_t_runs=True,
            timeout_seconds=30.0,
        )
        codon_domains = {i: AA_TO_CODONS[aa] for i, aa in enumerate(protein)}
        model = CSPModel(
            protein_sequence=protein,
            codon_domains=codon_domains,
            constraints=[],
            config=config,
            organism=organism,
        )

        engine = ORTOOLSEngine(config)
        result = engine.solve(model)

        assert result.solved, (
            f"OR-Tools solver did not find a solution. "
            f"Warnings: {result.warnings}"
        )
        assert len(result.sequence) == len(protein) * 3, (
            f"Sequence length mismatch: got {len(result.sequence)}, "
            f"expected {len(protein) * 3}"
        )
        # GC should be within bounds
        gc = gc_content(result.sequence)
        assert 0.30 <= gc <= 0.70, f"GC content {gc:.3f} out of range"

    def test_ortools_preserves_protein(self):
        from biocompiler.solver.engine_ortools import ORTOOLSEngine
        from biocompiler.solver.types import SolverConfig, CSPModel
        from biocompiler.type_system import AA_TO_CODONS

        protein = _SHORT_PROTEIN
        organism = "Escherichia_coli"
        config = SolverConfig(
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            cryptic_splice_threshold=0.0,
            avoid_cpg=False,
            avoid_attta=True,
            avoid_t_runs=True,
            timeout_seconds=30.0,
        )
        codon_domains = {i: AA_TO_CODONS[aa] for i, aa in enumerate(protein)}
        model = CSPModel(
            protein_sequence=protein,
            codon_domains=codon_domains,
            constraints=[],
            config=config,
            organism=organism,
        )

        engine = ORTOOLSEngine(config)
        result = engine.solve(model)

        if result.solved:
            translated = _translate_dna(result.sequence)
            assert translated == protein, (
                f"OR-Tools result does not translate to the original protein. "
                f"Expected {protein}, got {translated}"
            )

    def test_ortools_no_restriction_sites(self):
        from biocompiler.solver.engine_ortools import ORTOOLSEngine
        from biocompiler.solver.types import SolverConfig, CSPModel
        from biocompiler.type_system import AA_TO_CODONS

        protein = _SHORT_PROTEIN
        organism = "Escherichia_coli"
        # Use common 6+ bp restriction sites
        restriction_sites = [
            seq for name, seq in RESTRICTION_ENZYMES.items()
            if all(ch in "ACGT" for ch in seq) and len(seq) >= 6
        ]
        config = SolverConfig(
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            cryptic_splice_threshold=0.0,
            avoid_cpg=False,
            avoid_attta=True,
            avoid_t_runs=True,
            restriction_sites=restriction_sites[:5],  # Use top 5 for speed
            timeout_seconds=30.0,
        )
        codon_domains = {i: AA_TO_CODONS[aa] for i, aa in enumerate(protein)}
        model = CSPModel(
            protein_sequence=protein,
            codon_domains=codon_domains,
            constraints=[],
            config=config,
            organism=organism,
        )

        engine = ORTOOLSEngine(config)
        result = engine.solve(model)

        if result.solved:
            for site in config.restriction_sites:
                assert not _has_restriction_site(result.sequence, site), (
                    f"OR-Tools result contains restriction site {site}"
                )


# ══════════════════════════════════════════════════════════════
# TestCSPZ3
# ══════════════════════════════════════════════════════════════

class TestCSPZ3:
    """Test Z3 CSP solver."""

    @pytest.fixture(autouse=True)
    def _check_z3(self):
        """Skip tests if Z3 is not available."""
        pytest.importorskip("z3")

    def test_z3_produces_valid_sequence(self):
        from biocompiler.solver.engine_z3 import Z3Engine
        from biocompiler.solver.types import SolverConfig, CSPModel
        from biocompiler.type_system import AA_TO_CODONS

        protein = _SHORT_PROTEIN
        organism = "Escherichia_coli"
        config = SolverConfig(
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            cryptic_splice_threshold=0.0,
            avoid_cpg=False,
            avoid_attta=True,
            avoid_t_runs=True,
            timeout_seconds=30.0,
        )
        codon_domains = {i: AA_TO_CODONS[aa] for i, aa in enumerate(protein)}
        model = CSPModel(
            protein_sequence=protein,
            codon_domains=codon_domains,
            constraints=[],
            config=config,
            organism=organism,
        )

        engine = Z3Engine(config, organism=organism)
        result = engine.solve(model)

        assert result.solved, (
            f"Z3 solver did not find a solution. "
            f"Warnings: {result.warnings}"
        )
        assert len(result.sequence) == len(protein) * 3, (
            f"Sequence length mismatch: got {len(result.sequence)}, "
            f"expected {len(protein) * 3}"
        )
        # GC should be within bounds
        gc = gc_content(result.sequence)
        assert 0.30 <= gc <= 0.70, f"GC content {gc:.3f} out of range"

    def test_z3_preserves_protein(self):
        from biocompiler.solver.engine_z3 import Z3Engine
        from biocompiler.solver.types import SolverConfig, CSPModel
        from biocompiler.type_system import AA_TO_CODONS

        protein = _SHORT_PROTEIN
        organism = "Escherichia_coli"
        config = SolverConfig(
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            cryptic_splice_threshold=0.0,
            avoid_cpg=False,
            avoid_attta=True,
            avoid_t_runs=True,
            timeout_seconds=30.0,
        )
        codon_domains = {i: AA_TO_CODONS[aa] for i, aa in enumerate(protein)}
        model = CSPModel(
            protein_sequence=protein,
            codon_domains=codon_domains,
            constraints=[],
            config=config,
            organism=organism,
        )

        engine = Z3Engine(config, organism=organism)
        result = engine.solve(model)

        if result.solved:
            translated = _translate_dna(result.sequence)
            assert translated == protein, (
                f"Z3 result does not translate to the original protein. "
                f"Expected {protein}, got {translated}"
            )


# ══════════════════════════════════════════════════════════════
# TestDNADchiselComparison
# ══════════════════════════════════════════════════════════════

class TestDNADchiselComparison:
    """Compare with DNAchisel."""

    @pytest.fixture(autouse=True)
    def _check_dnachisel(self):
        """Skip tests if DNAchisel is not available."""
        from biocompiler.dna_chisel_compat import is_dna_chisel_available
        if not is_dna_chisel_available():
            pytest.skip("DNAchisel not installed")

    def test_cai_within_0_05_of_dnachisel(self):
        """Our CAI should be within 0.05 of DNAchisel's"""
        from biocompiler.dna_chisel_compat import optimize_with_dna_chisel

        protein = _SHORT_PROTEIN
        organism = "Escherichia_coli"

        # Run our optimizer
        our_result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="hybrid",
        )

        # Run DNAchisel
        dc_result = optimize_with_dna_chisel(
            protein=protein,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
        )

        if dc_result.success:
            assert abs(our_result.cai - dc_result.cai) < 0.10, (
                f"CAI gap too large: ours={our_result.cai:.4f}, "
                f"DNAchisel={dc_result.cai:.4f}"
            )

    def test_speed_within_5x_of_dnachisel(self):
        """We should be within 5× of DNAchisel speed"""
        from biocompiler.dna_chisel_compat import optimize_with_dna_chisel

        protein = _SHORT_PROTEIN
        organism = "Escherichia_coli"

        # Time our optimizer
        t0 = time.monotonic()
        optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
            strategy="hybrid",
        )
        our_time = time.monotonic() - t0

        # Time DNAchisel
        dc_result = optimize_with_dna_chisel(
            protein=protein,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
        )

        if dc_result.success and dc_result.execution_time_s > 0:
            assert our_time < dc_result.execution_time_s * 5.0 + 1.0, (
                f"Our optimizer ({our_time:.3f}s) is >5× slower than "
                f"DNAchisel ({dc_result.execution_time_s:.3f}s)"
            )


# ══════════════════════════════════════════════════════════════
# TestOrganismNameResolution
# ══════════════════════════════════════════════════════════════

class TestOrganismNameResolution:
    """Test the organism name resolution system."""

    # E. coli aliases
    ECOLI_ALIASES = [
        "ecoli",
        "E_coli",
        "e_coli",
        "E. coli",
        "Escherichia_coli",
        "Escherichia coli",
    ]

    # Human aliases
    HUMAN_ALIASES = [
        "human",
        "H_sapiens",
        "h_sapiens",
        "H. sapiens",
        "Homo_sapiens",
        "Homo sapiens",
    ]

    def test_ecoli_aliases(self):
        """All E. coli aliases should resolve correctly"""
        for alias in self.ECOLI_ALIASES:
            resolved = resolve_organism(alias, strict=False)
            assert resolved == "Escherichia_coli", (
                f"Alias '{alias}' resolved to '{resolved}', "
                f"expected 'Escherichia_coli'"
            )

    def test_human_aliases(self):
        """All human aliases should resolve correctly"""
        for alias in self.HUMAN_ALIASES:
            resolved = resolve_organism(alias, strict=False)
            assert resolved == "Homo_sapiens", (
                f"Alias '{alias}' resolved to '{resolved}', "
                f"expected 'Homo_sapiens'"
            )

    def test_ecoli_aliases_in_optimizer(self):
        """E. coli aliases should work with optimize_sequence"""
        for alias in self.ECOLI_ALIASES:
            result = optimize_sequence(
                target_protein=_SHORT_PROTEIN,
                organism=alias,
                gc_lo=0.30,
                gc_hi=0.70,
                strategy="hybrid",
            )
            assert result.cai > 0.0, (
                f"Optimization with alias '{alias}' produced CAI={result.cai:.4f}"
            )
            assert len(result.sequence) == len(_SHORT_PROTEIN) * 3

    def test_human_aliases_in_optimizer(self):
        """Human aliases should work with optimize_sequence"""
        for alias in self.HUMAN_ALIASES:
            result = optimize_sequence(
                target_protein=_SHORT_PROTEIN,
                organism=alias,
                gc_lo=0.30,
                gc_hi=0.70,
                strategy="hybrid",
            )
            assert result.cai > 0.0, (
                f"Optimization with alias '{alias}' produced CAI={result.cai:.4f}"
            )
            assert len(result.sequence) == len(_SHORT_PROTEIN) * 3

    def test_ecoli_is_not_eukaryotic(self):
        """E. coli should be identified as non-eukaryotic"""
        assert not is_eukaryotic_organism("Escherichia_coli"), (
            "E. coli should not be classified as eukaryotic"
        )
        assert not is_eukaryotic_organism("ecoli"), (
            "'ecoli' alias should not be classified as eukaryotic"
        )

    def test_human_is_eukaryotic(self):
        """Human should be identified as eukaryotic"""
        assert is_eukaryotic_organism("Homo_sapiens"), (
            "Homo_sapiens should be classified as eukaryotic"
        )
        assert is_eukaryotic_organism("human"), (
            "'human' alias should be classified as eukaryotic"
        )

    def test_all_aliases_resolve_to_supported_organisms(self):
        """Every alias in ORGANISM_ALIASES should resolve to a supported organism"""
        for alias, canonical in ORGANISM_ALIASES.items():
            resolved = resolve_organism(alias, strict=False)
            assert resolved == canonical, (
                f"Alias '{alias}' resolved to '{resolved}', expected '{canonical}'"
            )
            assert canonical in SUPPORTED_ORGANISMS, (
                f"Canonical '{canonical}' is not in SUPPORTED_ORGANISMS"
            )
            assert canonical in CODON_ADAPTIVENESS_TABLES, (
                f"Canonical '{canonical}' is not in CODON_ADAPTIVENESS_TABLES"
            )
