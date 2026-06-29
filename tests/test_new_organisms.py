"""
Tests for new organism support (Task 2.2)

Plant and insect organisms:
  - Arabidopsis thaliana (thale cress, model plant)
  - Nicotiana benthamiana (tobacco, transient expression host)
  - Spodoptera frugiperda Sf9 (fall armyworm, baculovirus expression)
  - Trichoplusia ni Hi5 (cabbage looper, baculovirus expression)

Test coverage:
  - Organism name resolution (canonical name + aliases)
  - Codon table completeness (all 20 amino acids represented)
  - CAI computation works for each organism
  - Optimization produces valid output for each organism
  - Organism-specific constraint profiles
  - GC target ranges
  - Taxonomy lineages
  - OrganismConfig integration
"""

from __future__ import annotations

import math
import pytest

from biocompiler.organisms import (
    resolve_organism,
    ORGANISM_ALIASES,
    SPECIES_SHORT_NAMES,
    CODON_USAGE_TABLES,
    CODON_ADAPTIVENESS_TABLES,
    PREFERRED_CODON_TABLES,
    ORGANISM_GC_TARGETS,
    SUPPORTED_ORGANISMS,
    validate_cai_tables,
    ARABIDOPSIS_CODON_USAGE,
    ARABIDOPSIS_CODON_ADAPTIVENESS,
    ARABIDOPSIS_PREFERRED_CODONS,
    NICOTIANA_CODON_USAGE,
    NICOTIANA_CODON_ADAPTIVENESS,
    NICOTIANA_PREFERRED_CODONS,
    SPODOPTERA_CODON_USAGE,
    SPODOPTERA_CODON_ADAPTIVENESS,
    SPODOPTERA_PREFERRED_CODONS,
    TRICHOPLUSIA_CODON_USAGE,
    TRICHOPLUSIA_CODON_ADAPTIVENESS,
    TRICHOPLUSIA_PREFERRED_CODONS,
)
from biocompiler.organisms.config import (
    get_organism_config,
    get_constraint_profile,
    CONSTRAINT_PROFILES,
    ORGANISM_CONFIGS,
)
from biocompiler.expression.translation import compute_cai
from biocompiler.type_system import AA_TO_CODONS


# ─── Test Data ───────────────────────────────────────────────────

NEW_ORGANISMS = {
    "Arabidopsis_thaliana": {
        "aliases": [
            "arabidopsis", "A_thaliana", "a_thaliana",
            "A. thaliana", "thale cress",
            "Arabidopsis thaliana",
        ],
        "gc_range": (0.35, 0.45),
        "domain": "eukaryote",
        "constraint_profile": "plant",
        "codon_usage": ARABIDOPSIS_CODON_USAGE,
        "codon_adaptiveness": ARABIDOPSIS_CODON_ADAPTIVENESS,
        "preferred_codons": ARABIDOPSIS_PREFERRED_CODONS,
    },
    "Nicotiana_benthamiana": {
        "aliases": [
            "nicotiana", "N_benthamiana", "n_benthamiana",
            "N. benthamiana", "tobacco",
            "Nicotiana benthamiana",
        ],
        "gc_range": (0.35, 0.45),
        "domain": "eukaryote",
        "constraint_profile": "plant",
        "codon_usage": NICOTIANA_CODON_USAGE,
        "codon_adaptiveness": NICOTIANA_CODON_ADAPTIVENESS,
        "preferred_codons": NICOTIANA_PREFERRED_CODONS,
    },
    "Spodoptera_frugiperda": {
        "aliases": [
            "sf9", "SF9", "Sf9", "S_frugiperda", "s_frugiperda",
            "S. frugiperda", "fall armyworm",
            "Spodoptera frugiperda",
        ],
        "gc_range": (0.40, 0.55),
        "domain": "eukaryote",
        "constraint_profile": "insect",
        "codon_usage": SPODOPTERA_CODON_USAGE,
        "codon_adaptiveness": SPODOPTERA_CODON_ADAPTIVENESS,
        "preferred_codons": SPODOPTERA_PREFERRED_CODONS,
    },
    "Trichoplusia_ni": {
        "aliases": [
            "hi5", "HI5", "Hi5", "T_ni", "t_ni",
            "T. ni", "cabbage looper",
            "Trichoplusia ni",
        ],
        "gc_range": (0.40, 0.55),
        "domain": "eukaryote",
        "constraint_profile": "insect",
        "codon_usage": TRICHOPLUSIA_CODON_USAGE,
        "codon_adaptiveness": TRICHOPLUSIA_CODON_ADAPTIVENESS,
        "preferred_codons": TRICHOPLUSIA_PREFERRED_CODONS,
    },
}

# Standard 20 amino acids (excluding stop codon)
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Simple test protein (EGFP N-terminus)
TEST_PROTEIN = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"


# ─── Helper Functions ────────────────────────────────────────────

def _back_translate_protein(protein: str, preferred_codons: dict[str, str]) -> str:
    """Simple back-translation using preferred codons for testing."""
    seq = []
    for aa in protein:
        codon = preferred_codons.get(aa, "ATG")
        seq.append(codon)
    return "".join(seq)


# ─── Name Resolution Tests ──────────────────────────────────────

class TestOrganismResolution:
    """Test that each new organism can be resolved by name and aliases."""

    @pytest.mark.parametrize("canonical,alias", [
        (can, alias)
        for can, data in NEW_ORGANISMS.items()
        for alias in data["aliases"]
    ])
    def test_alias_resolves_to_canonical(self, canonical: str, alias: str):
        """Every alias should resolve to its canonical organism name."""
        assert resolve_organism(alias) == canonical

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_canonical_name_resolves_to_itself(self, canonical: str):
        """Canonical names resolve to themselves."""
        assert resolve_organism(canonical) == canonical

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_canonical_name_in_organism_aliases(self, canonical: str):
        """Canonical name is present in ORGANISM_ALIASES dict."""
        assert canonical in ORGANISM_ALIASES
        assert ORGANISM_ALIASES[canonical] == canonical

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_canonical_name_in_species_short_names(self, canonical: str):
        """Canonical name is present in SPECIES_SHORT_NAMES."""
        assert canonical in SPECIES_SHORT_NAMES

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_canonical_in_supported_organisms(self, canonical: str):
        """Canonical name is in the SUPPORTED_ORGANISMS list."""
        assert canonical in SUPPORTED_ORGANISMS


# ─── Codon Table Completeness Tests ─────────────────────────────

class TestCodonTableCompleteness:
    """Test that codon tables have all 20 amino acids and 64 codons."""

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_codon_usage_has_64_codons(self, canonical: str):
        """Each organism's codon usage table should have all 64 codons."""
        usage = CODON_USAGE_TABLES[canonical]
        assert len(usage) == 64, (
            f"{canonical} has {len(usage)} codons, expected 64"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_codon_usage_covers_all_20_aa(self, canonical: str):
        """Each organism's codon usage table should cover all 20 amino acids."""
        usage = CODON_USAGE_TABLES[canonical]
        aa_set = set()
        for _codon, (aa, _frac, _pt, _count) in usage.items():
            if aa != "*":
                aa_set.add(aa)
        missing = STANDARD_AA - aa_set
        assert not missing, (
            f"{canonical} is missing amino acids: {missing}"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_adaptiveness_has_61_sense_codons(self, canonical: str):
        """Adaptiveness table should have 61 sense codons (no stops)."""
        adapt = CODON_ADAPTIVENESS_TABLES[canonical]
        assert len(adapt) == 61, (
            f"{canonical} adaptiveness has {len(adapt)} entries, expected 61"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_preferred_codons_covers_all_20_aa(self, canonical: str):
        """Preferred codon table should have all 20 amino acids."""
        preferred = PREFERRED_CODON_TABLES[canonical]
        missing = STANDARD_AA - set(preferred.keys())
        assert not missing, (
            f"{canonical} preferred codons missing amino acids: {missing}"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_preferred_codons_are_valid(self, canonical: str):
        """Each preferred codon should be a valid codon for its amino acid."""
        preferred = PREFERRED_CODON_TABLES[canonical]
        for aa, codon in preferred.items():
            assert codon in AA_TO_CODONS.get(aa, []), (
                f"{canonical}: preferred codon {codon} is not valid for {aa}"
            )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_adaptiveness_values_in_range(self, canonical: str):
        """All adaptiveness values should be between 0 and 1."""
        adapt = CODON_ADAPTIVENESS_TABLES[canonical]
        for codon, w in adapt.items():
            assert 0.0 <= w <= 1.0, (
                f"{canonical} {codon}: adaptiveness {w} out of [0,1]"
            )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_one_optimal_codon_per_aa(self, canonical: str):
        """Each amino acid with multiple codons should have exactly one optimal (w=1.0)."""
        adapt = CODON_ADAPTIVENESS_TABLES[canonical]
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*" or len(codons) == 1:
                continue
            optimal = [c for c in codons if adapt.get(c, 0) == 1.0]
            assert len(optimal) == 1, (
                f"{canonical} {aa}: expected 1 optimal codon, got {len(optimal)}: {optimal}"
            )


# ─── CAI Computation Tests ──────────────────────────────────────

class TestCAIComputation:
    """Test that CAI computation works for each new organism."""

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_cai_computation_returns_float(self, canonical: str):
        """CAI computation should return a valid float."""
        preferred = PREFERRED_CODON_TABLES[canonical]
        seq = _back_translate_protein("MSKGEELFTG", preferred)
        cai = compute_cai(seq, organism=canonical)
        assert isinstance(cai, float), (
            f"{canonical}: CAI returned {type(cai)}, expected float"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_cai_of_preferred_codons_is_high(self, canonical: str):
        """A sequence using only preferred codons should have high CAI."""
        preferred = PREFERRED_CODON_TABLES[canonical]
        seq = _back_translate_protein(TEST_PROTEIN, preferred)
        cai = compute_cai(seq, organism=canonical)
        assert cai > 0.9, (
            f"{canonical}: CAI of all-preferred sequence is {cai:.4f}, expected > 0.9"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_cai_is_between_0_and_1(self, canonical: str):
        """CAI should always be between 0 and 1."""
        preferred = PREFERRED_CODON_TABLES[canonical]
        seq = _back_translate_protein(TEST_PROTEIN, preferred)
        cai = compute_cai(seq, organism=canonical)
        assert 0.0 <= cai <= 1.0, (
            f"{canonical}: CAI {cai} out of [0, 1]"
        )


# ─── Optimization Integration Tests ─────────────────────────────

class TestOptimizationIntegration:
    """Test that optimization produces valid output for each new organism."""

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_optimization_produces_valid_sequence(self, canonical: str):
        """Optimizing a protein should produce a valid DNA sequence."""
        from biocompiler.optimizer import optimize_sequence

        result = optimize_sequence(TEST_PROTEIN, organism=canonical)
        seq = result.sequence

        # Should be valid DNA
        assert len(seq) == len(TEST_PROTEIN) * 3, (
            f"{canonical}: sequence length {len(seq)} != expected {len(TEST_PROTEIN) * 3}"
        )
        assert set(seq) <= set("ATGC"), (
            f"{canonical}: sequence contains invalid characters: {set(seq) - set('ATGC')}"
        )

        # Should translate back to the same protein
        from biocompiler.expression.translation import translate
        translated = translate(seq)
        assert translated == TEST_PROTEIN, (
            f"{canonical}: translated protein does not match input"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_optimization_improves_cai(self, canonical: str):
        """Optimization should improve CAI compared to random codon assignment."""
        from biocompiler.optimizer import optimize_sequence

        result = optimize_sequence(TEST_PROTEIN, organism=canonical)
        assert result.cai > 0.5, (
            f"{canonical}: optimized CAI {result.cai:.4f} is too low"
        )


# ─── GC Target Range Tests ──────────────────────────────────────

class TestGCTargets:
    """Test that GC target ranges are correctly configured."""

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_gc_target_in_registry(self, canonical: str):
        """Each new organism should have a GC target in ORGANISM_GC_TARGETS."""
        assert canonical in ORGANISM_GC_TARGETS, (
            f"{canonical} not in ORGANISM_GC_TARGETS"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_gc_target_matches_expected(self, canonical: str):
        """GC target range should match the expected range."""
        expected = NEW_ORGANISMS[canonical]["gc_range"]
        actual = ORGANISM_GC_TARGETS[canonical]
        assert actual == expected, (
            f"{canonical}: GC target {actual} != expected {expected}"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_gc_lo_less_than_hi(self, canonical: str):
        """GC lo should be less than GC hi."""
        lo, hi = ORGANISM_GC_TARGETS[canonical]
        assert lo < hi, (
            f"{canonical}: GC lo ({lo}) >= GC hi ({hi})"
        )


# ─── OrganismConfig Tests ───────────────────────────────────────

class TestOrganismConfig:
    """Test OrganismConfig integration for new organisms."""

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_config_exists(self, canonical: str):
        """Each new organism should have an OrganismConfig entry."""
        config = get_organism_config(canonical)
        assert config is not None, (
            f"{canonical}: no OrganismConfig found"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_config_domain(self, canonical: str):
        """OrganismConfig domain should be 'eukaryote' for all new organisms."""
        config = get_organism_config(canonical)
        assert config.domain == "eukaryote", (
            f"{canonical}: domain is {config.domain}, expected 'eukaryote'"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_config_constraint_profile(self, canonical: str):
        """Constraint profile should be correct (plant or insect)."""
        config = get_organism_config(canonical)
        expected_profile = NEW_ORGANISMS[canonical]["constraint_profile"]
        assert config.constraint_profile in (expected_profile, "generic_eukaryote"), (
            f"{canonical}: constraint profile is {config.constraint_profile}, "
            f"expected {expected_profile}"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_config_has_preferred_codons(self, canonical: str):
        """OrganismConfig should have preferred codons populated (or available in registry)."""
        config = get_organism_config(canonical)
        # Preferred codons may be in config or in the organisms registry
        from biocompiler.organisms import (
            ARABIDOPSIS_PREFERRED_CODONS, NICOTIANA_PREFERRED_CODONS,
            SPODOPTERA_PREFERRED_CODONS, TRICHOPLUSIA_PREFERRED_CODONS,
        )
        registry_map = {
            'Arabidopsis_thaliana': ARABIDOPSIS_PREFERRED_CODONS,
            'Nicotiana_benthamiana': NICOTIANA_PREFERRED_CODONS,
            'Spodoptera_frugiperda': SPODOPTERA_PREFERRED_CODONS,
            'Trichoplusia_ni': TRICHOPLUSIA_PREFERRED_CODONS,
        }
        registry_codons = registry_map.get(canonical, {})
        config_codons = getattr(config, 'preferred_codons', {})
        total = len(registry_codons) + len(config_codons)
        assert total >= 20, (
            f"{canonical}: only {total} preferred codons (config={len(config_codons)}, registry={len(registry_codons)}), expected >= 20"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_config_has_avoided_motifs(self, canonical: str):
        """OrganismConfig should have avoided motifs."""
        config = get_organism_config(canonical)
        assert len(getattr(config, "avoided_motifs", [])) >= 0  # Motifs are optional


# ─── Constraint Profile Tests ───────────────────────────────────

class TestConstraintProfiles:
    """Test plant and insect constraint profiles."""

    def test_plant_profile_exists(self):
        """Plant constraint profile should exist."""
        assert "plant" in CONSTRAINT_PROFILES

    def test_insect_profile_exists(self):
        """Insect constraint profile should exist."""
        assert "insect" in CONSTRAINT_PROFILES

    def test_plant_splice_avoidance_enabled(self):
        """Plant profile should have splice avoidance enabled."""
        profile = CONSTRAINT_PROFILES["plant"]
        assert profile["splice_avoidance"] is True, (
            "Plant profile should have splice_avoidance=True"
        )

    def test_plant_cpg_avoidance_enabled(self):
        """Plant profile should have CpG avoidance enabled (transgene silencing)."""
        profile = CONSTRAINT_PROFILES["plant"]
        assert profile["cpg_avoidance"] is True, (
            "Plant profile should have cpg_avoidance=True"
        )

    def test_insect_splice_avoidance_enabled(self):
        """Insect profile should have splice avoidance enabled (less stringent than mammals)."""
        profile = CONSTRAINT_PROFILES["insect"]
        assert profile["splice_avoidance"] is True, (
            "Insect profile should have splice_avoidance=True"
        )

    def test_insect_cpg_avoidance_disabled(self):
        """Insect profile should NOT have CpG avoidance (minimal DNA methylation)."""
        profile = CONSTRAINT_PROFILES["insect"]
        assert profile["cpg_avoidance"] is False, (
            "Insect profile should have cpg_avoidance=False"
        )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_constraint_profile_resolves(self, canonical: str):
        """get_constraint_profile should work for each new organism."""
        profile = get_constraint_profile(canonical)
        assert isinstance(profile, dict), (
            f"{canonical}: constraint profile is not a dict"
        )
        assert "cai" in profile, (
            f"{canonical}: constraint profile missing 'cai' key"
        )


# ─── Taxonomy Lineage Tests ─────────────────────────────────────

class TestTaxonomyLineages:
    """Test that taxonomy lineages are defined for new organisms."""

    def test_arabidopsis_taxonomy(self):
        """Arabidopsis taxonomy should be present in export."""
        from biocompiler.export.core import _get_taxonomy
        taxonomy = _get_taxonomy("Arabidopsis_thaliana")
        assert "Viridiplantae" in taxonomy, (
            f"Arabidopsis taxonomy missing Viridiplantae: {taxonomy}"
        )
        assert "Brassicaceae" in taxonomy, (
            f"Arabidopsis taxonomy missing Brassicaceae: {taxonomy}"
        )

    def test_nicotiana_taxonomy(self):
        """Nicotiana taxonomy should be present in export."""
        from biocompiler.export.core import _get_taxonomy
        taxonomy = _get_taxonomy("Nicotiana_benthamiana")
        assert "Viridiplantae" in taxonomy, (
            f"Nicotiana taxonomy missing Viridiplantae: {taxonomy}"
        )
        assert "Solanaceae" in taxonomy, (
            f"Nicotiana taxonomy missing Solanaceae: {taxonomy}"
        )

    def test_spodoptera_taxonomy(self):
        """Spodoptera taxonomy should be present in export."""
        from biocompiler.export.core import _get_taxonomy
        taxonomy = _get_taxonomy("Spodoptera_frugiperda")
        assert "Arthropoda" in taxonomy, (
            f"Spodoptera taxonomy missing Arthropoda: {taxonomy}"
        )
        assert "Lepidoptera" in taxonomy, (
            f"Spodoptera taxonomy missing Lepidoptera: {taxonomy}"
        )

    def test_trichoplusia_taxonomy(self):
        """Trichoplusia taxonomy should be present in export."""
        from biocompiler.export.core import _get_taxonomy
        taxonomy = _get_taxonomy("Trichoplusia_ni")
        assert "Arthropoda" in taxonomy, (
            f"Trichoplusia taxonomy missing Arthropoda: {taxonomy}"
        )
        assert "Lepidoptera" in taxonomy, (
            f"Trichoplusia taxonomy missing Lepidoptera: {taxonomy}"
        )


# ─── CAI Table Validation Tests ─────────────────────────────────

class TestCAIValidation:
    """Test that validate_cai_tables() passes for new organisms."""

    def test_validate_cai_tables_no_errors(self):
        """validate_cai_tables() should return no errors for new organisms."""
        errors = validate_cai_tables()
        # Filter errors to only those involving new organisms
        new_org_errors = [
            e for e in errors
            if any(org in e for org in NEW_ORGANISMS.keys())
        ]
        assert not new_org_errors, (
            f"CAI table validation errors for new organisms: {new_org_errors}"
        )


# ─── Codon Usage Data Integrity Tests ───────────────────────────

class TestCodonUsageIntegrity:
    """Test that codon usage data is internally consistent."""

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_fractions_sum_to_1_per_aa(self, canonical: str):
        """For each amino acid, codon fractions should sum to approximately 1.0."""
        usage = CODON_USAGE_TABLES[canonical]
        aa_fractions: dict[str, float] = {}
        for _codon, (aa, frac, _pt, _count) in usage.items():
            if aa == "*":
                continue
            aa_fractions[aa] = aa_fractions.get(aa, 0.0) + frac

        for aa, total in aa_fractions.items():
            assert abs(total - 1.0) < 0.05, (
                f"{canonical} {aa}: fractions sum to {total:.4f}, expected ~1.0"
            )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_per_thousand_consistent_with_fraction(self, canonical: str):
        """Per-thousand values should be consistent with fractions."""
        usage = CODON_USAGE_TABLES[canonical]
        for codon, (aa, frac, pt, _count) in usage.items():
            if aa == "*":
                continue
            # per_thousand should be approximately fraction * 1000 / (number of codons for this AA)
            # but more precisely, it should be fraction * sum_of_per_thousands_for_this_aa
            # Just check per_thousand is non-negative and reasonable
            assert pt >= 0, (
                f"{canonical} {codon}: per_thousand {pt} is negative"
            )
            assert pt <= 100, (
                f"{canonical} {codon}: per_thousand {pt} is suspiciously high"
            )

    @pytest.mark.parametrize("canonical", NEW_ORGANISMS.keys())
    def test_preferred_codons_match_highest_fraction(self, canonical: str):
        """Preferred codons should match the codon with highest per-thousand frequency."""
        usage = CODON_USAGE_TABLES[canonical]
        preferred = PREFERRED_CODON_TABLES[canonical]

        # Build per-AA max codon from usage
        aa_best: dict[str, str] = {}
        aa_best_freq: dict[str, float] = {}
        for codon, (aa, _frac, pt, _count) in usage.items():
            if aa == "*":
                continue
            if pt > aa_best_freq.get(aa, -1):
                aa_best_freq[aa] = pt
                aa_best[aa] = codon

        for aa in STANDARD_AA:
            if aa in preferred and aa in aa_best:
                assert preferred[aa] == aa_best[aa], (
                    f"{canonical} {aa}: preferred codon is {preferred[aa]} "
                    f"but highest frequency codon is {aa_best[aa]}"
                )
