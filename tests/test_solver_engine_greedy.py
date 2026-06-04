"""Tests for solver.engine_greedy.GreedyEngine.

Covers:
1. GreedyEngine construction (config, seed, defaults)
2. solve() returns SolverResult with correct fields
3. Deterministic results with seed
4. Fallback behavior on invalid input (unknown AA, unknown organism, empty domain)
5. Constraint satisfaction for basic cases (translation, GC, CAI)
"""

from __future__ import annotations

import math
from unittest.mock import patch

import pytest

from biocompiler.solver.engine_greedy import GreedyEngine
from biocompiler.solver.types import (
    CSPModel,
    SolverBackend,
    SolverConfig,
    SolverResult,
)
from biocompiler.constants import AA_TO_CODONS, CODON_TABLE
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES


# ── Fixtures ────────────────────────────────────────────────────────

@pytest.fixture
def default_config() -> SolverConfig:
    """Default SolverConfig for greedy engine tests."""
    return SolverConfig(gc_lo=0.30, gc_hi=0.70)


@pytest.fixture
def short_protein() -> str:
    """Short protein for quick tests: MVSKGE (6 AA)."""
    return "MVSKGE"


@pytest.fixture
def short_model(short_protein: str, default_config: SolverConfig) -> CSPModel:
    """Minimal CSPModel for the 6-AA short_protein."""
    codon_domains = {i: AA_TO_CODONS.get(aa, []) for i, aa in enumerate(short_protein)}
    return CSPModel(
        protein_sequence=short_protein,
        codon_domains=codon_domains,
        constraints=[],
        config=default_config,
    )


@pytest.fixture
def hbb_protein() -> str:
    """Human hemoglobin beta chain (147 AA) — realistic-length protein."""
    return (
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK"
        "EFTPPVQAAYQKVVAGVANALAHKYH"
    )


# ════════════════════════════════════════════════════════════════════
# 1. GreedyEngine construction
# ════════════════════════════════════════════════════════════════════

class TestGreedyEngineConstruction:
    """Tests for GreedyEngine.__init__."""

    def test_init_with_default_seed(self, default_config: SolverConfig):
        """Engine should be constructable with just a config (seed defaults to 0)."""
        engine = GreedyEngine(default_config)
        assert engine.config is default_config

    def test_init_with_explicit_seed(self, default_config: SolverConfig):
        """Engine should accept an explicit seed for deterministic tie-breaking."""
        engine = GreedyEngine(default_config, seed=42)
        assert engine.config is default_config

    def test_init_with_none_seed_uses_default(self, default_config: SolverConfig):
        """Passing seed=None should use the internal default (0)."""
        engine_none = GreedyEngine(default_config, seed=None)
        engine_zero = GreedyEngine(default_config, seed=0)
        # Both should produce identical RNG state → identical results
        model = CSPModel(
            protein_sequence="LLL",
            codon_domains={i: AA_TO_CODONS["L"] for i in range(3)},
            constraints=[],
            config=default_config,
        )
        assert engine_none.solve(model).sequence == engine_zero.solve(model).sequence

    def test_config_stored_as_attribute(self, default_config: SolverConfig):
        """The config passed at construction should be accessible as engine.config."""
        engine = GreedyEngine(default_config)
        assert engine.config.gc_lo == 0.30
        assert engine.config.gc_hi == 0.70


# ════════════════════════════════════════════════════════════════════
# 2. solve() returns SolverResult
# ════════════════════════════════════════════════════════════════════

class TestSolveReturnsSolverResult:
    """Tests that solve() returns a properly populated SolverResult."""

    def test_returns_solver_result_instance(self, default_config: SolverConfig, short_model: CSPModel):
        """solve() must return a SolverResult."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert isinstance(result, SolverResult)

    def test_solved_is_true(self, default_config: SolverConfig, short_model: CSPModel):
        """Greedy engine always reports solved=True (it always produces output)."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert result.solved is True

    def test_backend_used_is_greedy_fallback(self, default_config: SolverConfig, short_model: CSPModel):
        """The backend_used must be GREEDY_FALLBACK."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert result.backend_used == SolverBackend.GREEDY_FALLBACK

    def test_fallback_used_is_true(self, default_config: SolverConfig, short_model: CSPModel):
        """Greedy results always have fallback_used=True."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert result.fallback_used is True

    def test_sequence_length_is_three_times_protein_length(
        self, default_config: SolverConfig, short_protein: str, short_model: CSPModel
    ):
        """Output DNA sequence length must equal 3 * len(protein)."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert len(result.sequence) == 3 * len(short_protein)

    def test_protein_field_matches_input(self, default_config: SolverConfig, short_model: CSPModel):
        """Result.protein must match the input protein sequence (uppercased)."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert result.protein == short_model.protein_sequence.upper()

    def test_gc_content_is_populated(self, default_config: SolverConfig, short_model: CSPModel):
        """Result.gc_content must be a float in [0, 1]."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert isinstance(result.gc_content, float)
        assert 0.0 <= result.gc_content <= 1.0

    def test_cai_is_populated(self, default_config: SolverConfig, short_model: CSPModel):
        """Result.cai must be a float in [0, 1]."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert isinstance(result.cai, float)
        assert 0.0 <= result.cai <= 1.0

    def test_solve_time_seconds_non_negative(self, default_config: SolverConfig, short_model: CSPModel):
        """Solve time should be a non-negative float."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert result.solve_time_seconds >= 0.0

    def test_warnings_contain_fallback_message(self, default_config: SolverConfig, short_model: CSPModel):
        """Warnings should indicate greedy fallback was used."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert any("greedy" in w.lower() or "fallback" in w.lower() for w in result.warnings)

    def test_sequence_contains_only_valid_bases(self, default_config: SolverConfig, short_model: CSPModel):
        """Output sequence must only contain A, C, G, T characters."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert set(result.sequence).issubset({"A", "C", "G", "T"})

    def test_organism_field_populated(self, default_config: SolverConfig, short_model: CSPModel):
        """Result.organism should reflect the organism used (Homo_sapiens by default)."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        # Default fallback is Homo_sapiens
        assert result.organism == "Homo_sapiens"


# ════════════════════════════════════════════════════════════════════
# 3. Deterministic results with seed
# ════════════════════════════════════════════════════════════════════

class TestDeterministicResults:
    """Tests that the engine produces deterministic, reproducible results with seed."""

    def test_same_seed_same_result(self, default_config: SolverConfig, short_model: CSPModel):
        """Two engines with the same seed must produce identical sequences."""
        engine_a = GreedyEngine(default_config, seed=123)
        engine_b = GreedyEngine(default_config, seed=123)
        result_a = engine_a.solve(short_model)
        result_b = engine_b.solve(short_model)
        assert result_a.sequence == result_b.sequence

    def test_different_seeds_may_differ(self, default_config: SolverConfig):
        """Engines with different seeds may produce different sequences when ties exist.

        Leucine (L) has 6 codons, several with similar CAI. Running on a
        leucine-rich protein maximises the chance of observing tie-breaking.
        """
        protein = "LLLLLLLLLL"  # 10 leucines
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={i: AA_TO_CODONS["L"] for i in range(len(protein))},
            constraints=[],
            config=default_config,
        )
        result_a = GreedyEngine(default_config, seed=0).solve(model)
        result_b = GreedyEngine(default_config, seed=99).solve(model)
        # They *may* differ — but both must still translate to LLLLLLLLLL
        assert len(result_a.sequence) == 30
        assert len(result_b.sequence) == 30

    def test_default_seed_is_deterministic(self, default_config: SolverConfig, short_model: CSPModel):
        """Engines constructed without explicit seed should produce identical results."""
        engine_a = GreedyEngine(default_config)
        engine_b = GreedyEngine(default_config)
        result_a = engine_a.solve(short_model)
        result_b = engine_b.solve(short_model)
        assert result_a.sequence == result_b.sequence

    def test_repeated_solve_same_engine_same_result(self, default_config: SolverConfig, short_model: CSPModel):
        """Calling solve() multiple times on the same engine should produce identical results."""
        engine = GreedyEngine(default_config, seed=7)
        result_1 = engine.solve(short_model)
        result_2 = engine.solve(short_model)
        assert result_1.sequence == result_2.sequence

    def test_deterministic_with_hbb_protein(self, default_config: SolverConfig, hbb_protein: str):
        """Determinism on a realistic-length protein (147 AA)."""
        model = CSPModel(
            protein_sequence=hbb_protein,
            codon_domains={i: AA_TO_CODONS.get(aa, []) for i, aa in enumerate(hbb_protein)},
            constraints=[],
            config=default_config,
        )
        result_a = GreedyEngine(default_config, seed=42).solve(model)
        result_b = GreedyEngine(default_config, seed=42).solve(model)
        assert result_a.sequence == result_b.sequence


# ════════════════════════════════════════════════════════════════════
# 4. Fallback behavior on invalid input
# ════════════════════════════════════════════════════════════════════

class TestFallbackBehavior:
    """Tests for graceful degradation on invalid/unusual inputs."""

    def test_unknown_amino_acid_produces_placeholder(
        self, default_config: SolverConfig
    ):
        """Unknown amino acid (not in AA_TO_CODONS) should produce 'NNN' placeholder."""
        protein = "MX"  # X is not a standard AA
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={0: AA_TO_CODONS["M"], 1: []},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        # M should produce ATG; X should produce NNN
        assert result.sequence.startswith("ATG")
        assert result.sequence[3:] == "NNN"

    def test_unknown_organism_falls_back_to_homo_sapiens(
        self, default_config: SolverConfig
    ):
        """If the organism has no adaptiveness table, fallback organisms are tried."""
        protein = "MK"
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={0: AA_TO_CODONS["M"], 1: AA_TO_CODONS["K"]},
            constraints=[],
            config=default_config,
        )
        # Set a completely unknown organism
        model.config._organism = "UnknownOrganism_xyz"
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        # Should still produce a result (falling back to Homo_sapiens)
        assert result.solved is True
        assert len(result.sequence) == 6

    def test_all_unknown_amino_acids_produces_all_nnn(
        self, default_config: SolverConfig
    ):
        """Protein with entirely unknown amino acids should produce all NNN codons."""
        protein = "ZZZ"
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={i: [] for i in range(3)},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert result.sequence == "NNNNNNNNN"
        assert result.cai == 0.0  # No valid codons → CAI = 0

    def test_empty_protein_produces_empty_sequence(
        self, default_config: SolverConfig
    ):
        """An empty protein sequence should produce an empty DNA sequence."""
        model = CSPModel(
            protein_sequence="",
            codon_domains={},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert result.sequence == ""
        assert result.gc_content == 0.0

    def test_single_amino_acid_protein(self, default_config: SolverConfig):
        """Single amino acid protein should produce exactly one codon."""
        model = CSPModel(
            protein_sequence="M",
            codon_domains={0: AA_TO_CODONS["M"]},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert result.sequence == "ATG"
        assert len(result.sequence) == 3

    def test_mixed_valid_and_invalid_amino_acids(self, default_config: SolverConfig):
        """Protein with a mix of valid and invalid AAs should handle both."""
        protein = "MXK"  # M=valid, X=invalid, K=valid
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={0: AA_TO_CODONS["M"], 1: [], 2: AA_TO_CODONS["K"]},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert result.sequence[:3] == "ATG"  # M → ATG
        assert result.sequence[3:6] == "NNN"  # X → NNN
        assert len(result.sequence) == 9


# ════════════════════════════════════════════════════════════════════
# 5. Constraint satisfaction for basic cases
# ════════════════════════════════════════════════════════════════════

class TestConstraintSatisfaction:
    """Tests that the greedy engine satisfies basic biological constraints."""

    def test_translation_constraint_all_codons_valid(
        self, default_config: SolverConfig, short_protein: str, short_model: CSPModel
    ):
        """Every codon in the output should translate to the expected amino acid."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        for i, expected_aa in enumerate(short_protein):
            codon = result.sequence[i * 3 : i * 3 + 3]
            translated_aa = CODON_TABLE.get(codon)
            assert translated_aa == expected_aa, (
                f"Codon {codon} at position {i} translates to {translated_aa}, "
                f"expected {expected_aa}"
            )

    def test_translation_constraint_hbb(
        self, default_config: SolverConfig, hbb_protein: str
    ):
        """Translation check on a realistic-length protein (147 AA)."""
        model = CSPModel(
            protein_sequence=hbb_protein,
            codon_domains={i: AA_TO_CODONS.get(aa, []) for i, aa in enumerate(hbb_protein)},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        for i, expected_aa in enumerate(hbb_protein):
            codon = result.sequence[i * 3 : i * 3 + 3]
            translated_aa = CODON_TABLE.get(codon)
            assert translated_aa == expected_aa, (
                f"Codon {codon} at position {i} translates to {translated_aa}, "
                f"expected {expected_aa}"
            )

    def test_gc_content_computation_is_accurate(
        self, default_config: SolverConfig, short_model: CSPModel
    ):
        """Reported GC content should match manual calculation."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        manual_gc = sum(1 for b in result.sequence if b in "GC") / len(result.sequence)
        assert math.isclose(result.gc_content, manual_gc, abs_tol=1e-6), (
            f"Reported GC {result.gc_content} != manual {manual_gc}"
        )

    def test_gc_content_accurate_hbb(
        self, default_config: SolverConfig, hbb_protein: str
    ):
        """GC content accuracy on a realistic-length protein."""
        model = CSPModel(
            protein_sequence=hbb_protein,
            codon_domains={i: AA_TO_CODONS.get(aa, []) for i, aa in enumerate(hbb_protein)},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        manual_gc = sum(1 for b in result.sequence if b in "GC") / len(result.sequence)
        assert math.isclose(result.gc_content, manual_gc, abs_tol=1e-6)

    def test_cai_positive_for_homo_sapiens(
        self, default_config: SolverConfig, short_model: CSPModel
    ):
        """CAI should be positive for a valid protein with Homo_sapiens codon table."""
        engine = GreedyEngine(default_config)
        result = engine.solve(short_model)
        assert result.cai > 0.0, "CAI should be positive for a valid protein"

    def test_cai_high_for_greedy(self, default_config: SolverConfig, short_protein: str):
        """Greedy engine selects highest-CAI codons, so CAI should be relatively high."""
        model = CSPModel(
            protein_sequence=short_protein,
            codon_domains={i: AA_TO_CODONS.get(aa, []) for i, aa in enumerate(short_protein)},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        # Greedy picks the best codon per position, so CAI should be >= 0.5
        assert result.cai >= 0.5, f"CAI {result.cai:.3f} is unexpectedly low for greedy selection"

    def test_methionine_always_atg(self, default_config: SolverConfig):
        """Methionine has only one codon (ATG), so it must always be ATG."""
        protein = "MMM"
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={i: AA_TO_CODONS["M"] for i in range(3)},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        for i in range(3):
            codon = result.sequence[i * 3 : i * 3 + 3]
            assert codon == "ATG", f"Position {i}: expected ATG, got {codon}"

    def test_tryptophan_always_tgg(self, default_config: SolverConfig):
        """Tryptophan has only one codon (TGG), so it must always be TGG."""
        protein = "WWWW"
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={i: AA_TO_CODONS["W"] for i in range(4)},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        for i in range(4):
            codon = result.sequence[i * 3 : i * 3 + 3]
            assert codon == "TGG", f"Position {i}: expected TGG, got {codon}"

    def test_valine_codons_all_contain_gt(self, default_config: SolverConfig):
        """All Valine codons contain GT — greedy must pick one of them."""
        protein = "VVV"
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={i: AA_TO_CODONS["V"] for i in range(3)},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        for i in range(3):
            codon = result.sequence[i * 3 : i * 3 + 3]
            assert "GT" in codon, f"Valine codon {codon} should contain GT"
            assert CODON_TABLE.get(codon) == "V"

    def test_escherichia_coli_organism(self, default_config: SolverConfig):
        """Engine should work correctly with E. coli organism."""
        protein = "MK"
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={0: AA_TO_CODONS["M"], 1: AA_TO_CODONS["K"]},
            constraints=[],
            config=default_config,
        )
        model.config._organism = "Escherichia_coli"
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert result.solved is True
        assert result.organism == "Escherichia_coli"
        assert len(result.sequence) == 6

    def test_pick_best_codon_selects_highest_cai(self, default_config: SolverConfig):
        """_pick_best_codon should select the codon with highest adaptiveness."""
        engine = GreedyEngine(default_config, seed=0)
        adaptiveness = {"AAA": 0.2, "AAG": 0.9}
        result = engine._pick_best_codon(["AAA", "AAG"], adaptiveness)
        assert result == "AAG"

    def test_pick_best_codon_tie_breaks_with_seed(self, default_config: SolverConfig):
        """When multiple codons share the same CAI, seed-based tie-breaking is used."""
        engine = GreedyEngine(default_config, seed=0)
        adaptiveness = {"CTT": 0.5, "CTC": 0.5}  # Both have same score
        result = engine._pick_best_codon(["CTT", "CTC"], adaptiveness)
        # Should deterministically pick one
        assert result in ("CTT", "CTC")

    def test_pick_best_codon_empty_domain(self, default_config: SolverConfig):
        """_pick_best_codon with empty domain should return 'NNN'."""
        engine = GreedyEngine(default_config, seed=0)
        result = engine._pick_best_codon([], {})
        assert result == "NNN"

    def test_compute_gc_empty_string(self):
        """_compute_gc on empty string should return 0.0."""
        assert GreedyEngine._compute_gc("") == 0.0

    def test_compute_gc_all_gc(self):
        """_compute_gc on all-GC sequence should return 1.0."""
        assert GreedyEngine._compute_gc("GCGCGC") == 1.0

    def test_compute_gc_all_at(self):
        """_compute_gc on all-AT sequence should return 0.0."""
        assert GreedyEngine._compute_gc("ATATAT") == 0.0

    def test_compute_gc_mixed(self):
        """_compute_gc on mixed sequence should return correct fraction."""
        # ATGC → 2/4 = 0.5
        assert GreedyEngine._compute_gc("ATGC") == 0.5

    def test_compute_cai_unknown_organism(self):
        """_compute_cai with unknown organism should return 0.0."""
        result = GreedyEngine._compute_cai("ATGAAA", "MK", "NonExistent_organism")
        assert result == 0.0

    def test_compute_cai_valid_organism(self):
        """_compute_cai with a known organism should return a positive value."""
        result = GreedyEngine._compute_cai("ATGAAA", "MK", "Homo_sapiens")
        assert isinstance(result, float)
        assert result > 0.0


# ════════════════════════════════════════════════════════════════════
# 6. Edge cases & robustness
# ════════════════════════════════════════════════════════════════════

class TestEdgeCases:
    """Edge case tests for the GreedyEngine."""

    def test_case_insensitive_protein(self, default_config: SolverConfig):
        """Lowercase protein input should work (engine uppercases internally)."""
        protein = "mk"
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={0: AA_TO_CODONS["M"], 1: AA_TO_CODONS["K"]},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert result.solved is True
        assert len(result.sequence) == 6
        assert result.sequence.startswith("ATG")

    def test_all_same_amino_acid(self, default_config: SolverConfig):
        """Protein of all the same amino acid should produce valid output."""
        protein = "AAAAAAAAAA"  # 10 Alanines
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={i: AA_TO_CODONS["A"] for i in range(len(protein))},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert len(result.sequence) == 30
        for i in range(10):
            codon = result.sequence[i * 3 : i * 3 + 3]
            assert CODON_TABLE.get(codon) == "A", f"Codon {codon} doesn't encode Alanine"

    def test_protein_with_stop_codon_marker(self, default_config: SolverConfig):
        """Protein containing stop codon marker (*) should use NNN placeholder."""
        protein = "M*"  # * is not a standard amino acid
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={0: AA_TO_CODONS["M"], 1: []},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert result.sequence[:3] == "ATG"
        assert result.sequence[3:6] == "NNN"

    def test_model_without_organism_attribute(self, default_config: SolverConfig):
        """SolverConfig without _organism attribute should fall back gracefully."""
        protein = "MK"
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={0: AA_TO_CODONS["M"], 1: AA_TO_CODONS["K"]},
            constraints=[],
            config=default_config,
        )
        # Ensure no _organism attribute
        assert not hasattr(model.config, "_organism")
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert result.solved is True
        # Should default to Homo_sapiens
        assert result.organism == "Homo_sapiens"

    def test_long_protein_does_not_crash(self, default_config: SolverConfig, hbb_protein: str):
        """Realistic-length protein should solve without errors."""
        model = CSPModel(
            protein_sequence=hbb_protein,
            codon_domains={i: AA_TO_CODONS.get(aa, []) for i, aa in enumerate(hbb_protein)},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert result.solved is True
        assert len(result.sequence) == 3 * len(hbb_protein)
        assert result.solve_time_seconds >= 0.0

    def test_sequence_only_valid_dna_for_hbb(self, default_config: SolverConfig, hbb_protein: str):
        """All codons for HBB should be valid DNA (no NNN placeholders)."""
        model = CSPModel(
            protein_sequence=hbb_protein,
            codon_domains={i: AA_TO_CODONS.get(aa, []) for i, aa in enumerate(hbb_protein)},
            constraints=[],
            config=default_config,
        )
        engine = GreedyEngine(default_config)
        result = engine.solve(model)
        assert "N" not in result.sequence, "HBB protein should have no unknown codons"
