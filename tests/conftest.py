"""BioCompiler test suite — shared fixtures and pytest configuration.

This conftest.py centralizes commonly-used test data that was previously
duplicated across individual test files.  Import fixtures here instead of
redefining protein sequences, DNA sequences, organism names, solver configs,
immunogenicity test data, or optimization parameters in each test module.

Fixtures are organised into the following groups:
    1. Standard protein sequences
    2. DNA sequences
    3. Organisms
    4. Solver components
    5. Immunogenicity test data
    6. Optimization helpers
    7. Standard amino-acid data (duplicated in 43+ test files)
    8. Solubility test proteins
    9. Immunogenicity test proteins
    10. Structure test data
    11. Temp directories and I/O helpers
    12. Custom pytest markers
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

# Set a very high rate limit for tests so the API rate limiter doesn't
# interfere with test execution. This must be set before any biocompiler
# module is imported (the rate limit is read at import time).
os.environ.setdefault("BIOCOMPILER_RATE_LIMIT", "1000000")

import pytest


# ═══════════════════════════════════════════════════════════════════════════
# 1. Standard protein sequences
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def sample_protein() -> str:
    """Very short protein (6 AA) for quick unit tests.

    Corresponds to the first 6 residues of eGFP: M-V-S-K-G-E.
    """
    return "MVSKGE"


@pytest.fixture
def hbb_protein() -> str:
    """Human hemoglobin beta chain (HBB), 147 AA.

    A well-characterised human self-protein with low immunogenicity
    in humans.  Used in integration tests and the HBB full-pass suite.
    """
    return (
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK"
        "EFTPPVQAAYQKVVAGVANALAHKYH"
    )


@pytest.fixture
def egfp_protein() -> str:
    """Enhanced green fluorescent protein (eGFP), 239 AA.

    A standard reporter protein used across molecular biology.
    Frequently used for optimizer integration tests.
    """
    return (
        "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
        "VTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVN"
        "RIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQ"
        "QNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
    )


@pytest.fixture
def insulin_protein() -> str:
    """Human insulin B-chain (first 30 residues) for immunogenicity tests.

    A small, well-studied therapeutic protein commonly used in
    biopharmaceutical immunogenicity benchmarking.
    """
    return "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"


# ═══════════════════════════════════════════════════════════════════════════
# 2. DNA sequences
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def sample_dna() -> str:
    """DNA sequence encoding sample_protein ('MVSKGE').

    Uses the most common human codons for each amino acid:
    M=ATG  V=GTT  S=TCT  K=AAA  G=GGT  E=GAA
    """
    return "ATGGTTTCTAAAGGTGAA"


@pytest.fixture
def egfp_dna() -> str:
    """eGFP coding sequence (720 bp, 239 aa + stop codon).

    The canonical eGFP CDS used by the optimizer and CSP solver tests.
    """
    return (
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGG"
        "CCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCAC"
        "CACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCT"
        "ACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTT"
        "CTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGA"
        "GCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAA"
        "CGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGC"
        "AGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCAC"
        "TACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTG"
        "ACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
    )


# ═══════════════════════════════════════════════════════════════════════════
# 3. Organisms
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def default_organism() -> str:
    """Default organism name used throughout the test suite."""
    return "Homo_sapiens"


@pytest.fixture
def supported_organisms() -> list[str]:
    """List of all organisms natively supported by BioCompiler.

    Each entry is a valid key for the organisms module / SolverConfig.
    """
    return [
        "Homo_sapiens",
        "Mus_musculus",
        "Escherichia_coli",
        "Saccharomyces_cerevisiae",
        "CHO_K1",
        "Drosophila_melanogaster",
        "Caenorhabditis_elegans",
        "Danio_rerio",
        "Arabidopsis_thaliana",
        "Pichia_pastoris",
    ]


# ═══════════════════════════════════════════════════════════════════════════
# 4. Solver components
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def default_solver_config():
    """Default SolverConfig with sensible GC bounds and Homo_sapiens organism.

    Returns a ``biocompiler.solver.types.SolverConfig`` instance.  If the
    solver types module is not available the fixture is skipped.
    """
    mod = pytest.importorskip("biocompiler.solver.types")
    return mod.SolverConfig(gc_lo=0.30, gc_hi=0.70)


@pytest.fixture
def sample_csp_model():
    """A minimal CSPModel for the 6-AA sample_protein ('MVSKGE').

    Uses the default SolverConfig and builds codon domains from
    ``biocompiler.type_system.AA_TO_CODONS``.  Skipped if the solver
    types module is unavailable.
    """
    mod = pytest.importorskip("biocompiler.solver.types")
    from biocompiler.type_system import AA_TO_CODONS

    protein = "MVSKGE"
    codon_domains = {i: AA_TO_CODONS[aa] for i, aa in enumerate(protein)}
    config = mod.SolverConfig(gc_lo=0.30, gc_hi=0.70)
    return mod.CSPModel(
        protein_sequence=protein,
        codon_domains=codon_domains,
        constraints=[],
        config=config,
    )


# ═══════════════════════════════════════════════════════════════════════════
# 5. Immunogenicity test data
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def sample_alleles() -> list[str]:
    """Common MHC class I alleles used in immunogenicity tests.

    HLA-A*02:01 is the most frequent allele in many populations and
    is the best-studied for binding prediction benchmarks.
    """
    return [
        "HLA-A*02:01",
        "HLA-A*01:01",
        "HLA-B*07:02",
        "HLA-B*08:01",
    ]


@pytest.fixture
def sample_peptides() -> dict[str, str]:
    """Well-known MHC-binding peptides for validation tests.

    Each key is a descriptive name; the value is the peptide sequence.

    Includes:
      - GILGFVFTL: Influenza M1 epitope, strong HLA-A*02:01 binder
      - SIINFEKL:  OVA epitope, 8-mer, classic model peptide
      - ELAGIGILTV: Melanoma MART-1 epitope, HLA-A*02:01 binder
    """
    return {
        "influenza_m1": "GILGFVFTL",
        "ova": "SIINFEKL",
        "mart1": "ELAGIGILTV",
    }


# ═══════════════════════════════════════════════════════════════════════════
# 6. Optimization helpers
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def default_enzymes() -> list[str]:
    """Default list of restriction enzymes for optimization tests.

    These are the five most commonly-avoided cloning enzymes.
    """
    return ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]


@pytest.fixture
def sample_optimization_result():
    """A sample OptimizationResult for unit-testing downstream consumers.

    Uses the sample_protein ('MVSKGE') with reasonable default metrics.
    Skipped if the optimization module is unavailable.
    """
    mod = pytest.importorskip("biocompiler.optimizer")
    return mod.OptimizationResult(
        sequence="ATGGTTTCTAAAGGTGAA",
        gc_content=0.333,
        cai=0.78,
        failed_predicates=[],
        predicate_results=[],
        certificate_text="",
        protein="MVSKGE",
        fallback_used=False,
        satisfied_predicates=[],
        aa_substitutions=[],
        mutagenesis_applied=False,
    )


# ═══════════════════════════════════════════════════════════════════════════
# 7. Standard amino-acid data
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def standard_aas() -> str:
    """The 20 standard amino acids as a single-letter string.

    "ACDEFGHIKLMNPQRSTVWY" — duplicated in 43+ test files.
    Use this fixture in new tests instead of redefining it locally.
    """
    return "ACDEFGHIKLMNPQRSTVWY"


@pytest.fixture
def standard_aas_set() -> set[str]:
    """The 20 standard amino acids as a set of single-letter codes.

    Convenient for ``assert aa in standard_aas_set`` membership checks.
    """
    return set("ACDEFGHIKLMNPQRSTVWY")


@pytest.fixture
def standard_aas_three_letter() -> list[str]:
    """The 20 standard amino acids as three-letter codes.

    Ordered to match the canonical single-letter ordering
    A→ALA, C→CYS, D→ASP, E→GLU, F→PHE, G→GLY, H→HIS, I→ILE,
    K→LYS, L→LEU, M→MET, N→ASN, P→PRO, Q→GLN, R→ARG, S→SER,
    T→THR, V→VAL, W→TRP, Y→TYR.
    """
    return [
        "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
        "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER",
        "THR", "VAL", "TRP", "TYR",
    ]


# ═══════════════════════════════════════════════════════════════════════════
# 8. Solubility test proteins
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def soluble_protein() -> str:
    """Highly soluble protein (charged, hydrophilic).

    Duplicated in test_camsol.py, test_camsol_idp.py, test_foldx.py,
    test_immuno_predicates_ext.py, and test_extended_predicates.py.
    """
    return "MSEKKDKKEKEKKDEKKDEEKKDESKKDEEKKDEEKKDESKKDEEKKDEEKK"


@pytest.fixture
def insoluble_protein() -> str:
    """Hydrophobic / aggregation-prone protein.

    Duplicated in test_camsol.py, test_camsol_idp.py, test_foldx.py,
    test_immuno_predicates_ext.py, and test_extended_predicates.py.
    """
    return "MVVVIIVVVLLLFLLLLFFFFWWWAAAIIIMMM"


@pytest.fixture
def balanced_protein() -> str:
    """Protein with average amino-acid composition.

    Duplicated in test_camsol.py, test_camsol_idp.py, test_foldx.py,
    test_immuno_predicates_ext.py, and test_extended_predicates.py.
    """
    return "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAVDILSKKGDVQVIK"


# ═══════════════════════════════════════════════════════════════════════════
# 9. Immunogenicity test proteins
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def foreign_protein() -> str:
    """Foreign-like sequence with high predicted immunogenicity.

    Based on serum albumin — a non-self protein that scores high on
    immunogenicity heuristics.  Duplicated in test_immunogenicity.py,
    test_immunogenicity_unit.py, test_immuno_predicates_ext.py,
    test_deimmunization_unit.py, and test_mhcflurry_integration.py.
    """
    return (
        "MKWVTFISLLLLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDK"
        "SLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAP"
        "ELLYYANKY"
    )


@pytest.fixture
def short_peptide() -> str:
    """Known HLA-A*02:01 binder from influenza M1 protein (9-mer).

    "GILGFVFTL" — one of the best-characterised T-cell epitopes.
    Duplicated in test_immunogenicity.py, test_immunogenicity_unit.py,
    test_immuno_predicates_ext.py, and test_mhcflurry_integration.py.
    """
    return "GILGFVFTL"


# ═══════════════════════════════════════════════════════════════════════════
# 10. Structure test data
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def mini_pdb() -> str:
    """Minimal PDB string with 3 CA atoms (M-A-G).

    B-factors represent pLDDT scores: 85.0, 90.0, 78.0.
    Duplicated in test_esmfold.py, test_esmfold_unit.py,
    test_protein_analysis.py, and test_extended_predicates.py.
    """
    return (
        "ATOM      1  CA  MET A   1       1.000   2.000   3.000  1.00 85.00           C\n"
        "ATOM      2  CA  ALA A   2       4.000   2.000   3.000  1.00 90.00           C\n"
        "ATOM      3  CA  GLY A   3       7.000   2.000   3.000  1.00 78.00           C\n"
        "END\n"
    )


# ═══════════════════════════════════════════════════════════════════════════
# 11. Temp directories and I/O helpers
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def tmp_output_dir(tmp_path: Path) -> Path:
    """A temporary directory for output files (FASTA, GenBank, reports).

    Uses pytest's built-in tmp_path under the hood so the directory is
    automatically cleaned up after the test session.  Prefer this over
    calling tempfile.mkdtemp() manually.
    """
    out = tmp_path / "output"
    out.mkdir()
    return out


@pytest.fixture
def sample_organism_config():
    """A sample OrganismConfig for Homo_sapiens.

    Skipped if the organism_config module is unavailable.
    """
    mod = pytest.importorskip("biocompiler.organisms.config")
    return mod.get_organism_config("Homo_sapiens")


# ═══════════════════════════════════════════════════════════════════════════
# 12. Custom pytest markers
# ═══════════════════════════════════════════════════════════════════════════

def pytest_configure(config):
    """Register custom pytest markers used across the test suite.

    Markers:
        slow       — Tests that are computationally expensive or require
                     network access (e.g. MHCflurry model downloads,
                     full eGFP optimization).  Run with ``-m slow`` or
                     skip with ``-m "not slow"``.

        integration — End-to-end integration tests that exercise multiple
                      modules together.  May require external dependencies
                      (OR-Tools, Z3, MHCflurry, ViennaRNA).

        benchmark   — Benchmark tests comparing against other tools.

        e2e         — Full end-to-end pipeline tests.

        requires_mhcflurry  — Tests requiring the mhcflurry package.

        requires_netmhcpan  — Tests requiring the NetMHCpan tool.

        requires_dnachisel  — Tests requiring the DNAchisel package.

        requires_external   — Tests needing any external tool.
    """
    config.addinivalue_line("markers", "slow: computationally expensive or network-dependent tests")
    config.addinivalue_line("markers", "integration: end-to-end integration tests requiring multiple modules")
    config.addinivalue_line("markers", "benchmark: benchmark tests comparing against other tools")
    config.addinivalue_line("markers", "e2e: full end-to-end pipeline tests")
    config.addinivalue_line("markers", "requires_mhcflurry: tests requiring mhcflurry")
    config.addinivalue_line("markers", "requires_netmhcpan: tests requiring netmhcpan")
    config.addinivalue_line("markers", "requires_dnachisel: tests requiring DNAchisel")
    config.addinivalue_line("markers", "requires_external: tests needing any external tool")
    config.addinivalue_line("markers", "deprecation: tests verifying deprecated API warning behavior")


def pytest_collection_modifyitems(config, items):
    """Add the 'deprecation' marker to the deprecation test suite automatically."""
    for item in items:
        if "test_deprecation_suite" in str(item.fspath):
            item.add_marker(pytest.mark.deprecation)
