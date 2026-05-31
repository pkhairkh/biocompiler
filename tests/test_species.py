"""Test BioCompiler Species Data — CAI computation and species registry."""

import pytest
from biocompiler.species import (
    compute_cai_weights, ECOLI_CODON_USAGE, ECOLI_CAI,
    HUMAN_CODON_USAGE, HUMAN_CAI, SPECIES,
)
from biocompiler.type_system import AA_TO_CODONS


class TestCAIComputation:
    """Tests for CAI weight computation."""

    def test_ecoli_cai_highest_codon(self):
        """CTG should have CAI=1.0 for Leucine in E.coli (highest usage)."""
        # Leu codons: TTA, TTG, CTT, CTC, CTA, CTG
        # E.coli usage: TTA=7.6, TTG=11.0, CTT=10.5, CTC=10.5, CTA=3.9, CTG=51.0
        # CTG has the highest usage, so its CAI weight should be 1.0
        assert ECOLI_CAI["CTG"] == pytest.approx(1.0)

    def test_human_cai_range(self):
        """All CAI values for Human should be between 0 and 1 (inclusive)."""
        for codon, cai in HUMAN_CAI.items():
            assert 0.0 <= cai <= 1.0, f"CAI({codon})={cai} is out of range [0, 1]"

    def test_ecoli_cai_range(self):
        """All CAI values for E.coli should be between 0 and 1 (inclusive)."""
        for codon, cai in ECOLI_CAI.items():
            assert 0.0 <= cai <= 1.0, f"CAI({codon})={cai} is out of range [0, 1]"


class TestSpeciesRegistry:
    """Tests for the species registry."""

    def test_species_registry(self):
        """Both ecoli and human should be in the registry."""
        assert "ecoli" in SPECIES
        assert "human" in SPECIES

    def test_species_registry_has_cai_data(self):
        """Each species in registry should have non-empty CAI data."""
        for name, cai_data in SPECIES.items():
            assert len(cai_data) > 0, f"Species {name} has empty CAI data"

    def test_ecoli_highest_cai_per_aa(self):
        """For each amino acid, at least one codon should have CAI=1.0."""
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*":
                continue
            cai_values = [ECOLI_CAI.get(c, 0.0) for c in codons]
            assert max(cai_values) == pytest.approx(1.0), (
                f"No codon for AA {aa} has CAI=1.0 in E.coli"
            )
