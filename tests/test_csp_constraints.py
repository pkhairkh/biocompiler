"""
Comprehensive pytest tests for biocompiler.solver.constraints.

Covers:
1. CSP Model building — variables, domains, constraint registration
2. Constraint check functions — restriction sites, GC range, cryptic splice,
   CpG islands, ATTTA motif, T-runs
3. Helper functions — codon_gc_count, codon_contains_gt, codon_contains_cpg,
   compute_gc_from_codons
4. Variable domain correctness — AA_TO_CODONS agreement, stop codon exclusion,
   single-codon AAs
5. Edge cases — short proteins, single-codon-only proteins, tight GC bounds
"""

from __future__ import annotations

import pytest

from biocompiler.constants import (
    AA_TO_CODONS,
    CODON_TABLE,
    STOP_CODONS,
    RESTRICTION_ENZYMES,
    INSTABILITY_MOTIF,
)


def _import_constraints():
    """Import constraints module; skip tests if not yet implemented."""
    try:
        from biocompiler.solver import constraints
        return constraints
    except ImportError as exc:
        pytest.skip(f"solver.constraints not yet available: {exc}")


# ── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def sample_protein() -> str:
    """First 20 AA of human alpha-globin."""
    return "MVLSPADKTNVKAAWGKVGA"


@pytest.fixture
def default_config() -> dict:
    return {
        "organism": "Homo_sapiens",
        "gc_lo": 0.30, "gc_hi": 0.70,
        "restriction_sites": list(RESTRICTION_ENZYMES.values()),
        "cryptic_splice_threshold": 3.0,
        "avoid_cpg": True, "avoid_attta": True, "max_t_run": 5,
    }


@pytest.fixture
def tight_gc_config() -> dict:
    return {
        "organism": "Homo_sapiens",
        "gc_lo": 0.49, "gc_hi": 0.51,
        "restriction_sites": [],
        "cryptic_splice_threshold": 3.0,
        "avoid_cpg": False, "avoid_attta": False, "max_t_run": 6,
    }


@pytest.fixture
def no_restriction_config() -> dict:
    return {
        "organism": "Homo_sapiens",
        "gc_lo": 0.20, "gc_hi": 0.80,
        "restriction_sites": [],
        "cryptic_splice_threshold": 3.0,
        "avoid_cpg": False, "avoid_attta": False, "max_t_run": 6,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 1. CSP Model Building
# ═══════════════════════════════════════════════════════════════════════════════

class TestCSPModelBuilding:

    def test_correct_number_of_variables(self, sample_protein, default_config):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        assert len(model.variables) == len(sample_protein)

    def test_variable_domains_match_aa_to_codons(self, sample_protein, default_config):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        for i, aa in enumerate(sample_protein):
            assert set(model.variables[i].domain) == set(AA_TO_CODONS[aa]), (
                f"Position {i} AA={aa}: domain mismatch"
            )

    def test_hard_constraints_registered(self, sample_protein, default_config):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        assert len(model.hard_constraints) > 0

    def test_soft_constraints_registered(self, sample_protein, default_config):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        assert len(model.soft_constraints) > 0

    def test_model_stores_protein_and_organism(self, sample_protein, default_config):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        assert model.protein == sample_protein
        assert model.organism == "Homo_sapiens"

    def test_no_restriction_sites_config_omits_rs_constraints(self, sample_protein, no_restriction_config):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", no_restriction_config)
        rs_names = [n for n in (c_.name for c_ in model.hard_constraints)
                    if "restriction" in n.lower() or "rs_" in n.lower()]
        assert len(rs_names) == 0

    def test_build_model_deterministic(self, sample_protein, default_config):
        c = _import_constraints()
        m1 = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        m2 = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        assert len(m1.variables) == len(m2.variables)
        for v1, v2 in zip(m1.variables, m2.variables):
            assert set(v1.domain) == set(v2.domain)
        assert len(m1.hard_constraints) == len(m2.hard_constraints)
        assert len(m1.soft_constraints) == len(m2.soft_constraints)


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Constraint Check Functions
# ═══════════════════════════════════════════════════════════════════════════════

class TestNoRestrictionSiteConstraint:

    def test_gaattc_detected(self):
        c = _import_constraints()
        checker = c.NoRestrictionSiteConstraint(sites=["GAATTC"])
        # GAA (Glu) + TTC (Phe) → GAATTC (EcoRI)
        assert checker.is_violated(["GAA", "TTC"])

    def test_clean_sequence_passes(self):
        c = _import_constraints()
        checker = c.NoRestrictionSiteConstraint(sites=["GAATTC", "GGATCC"])
        assert not checker.is_violated(["AAA", "AAA"])

    def test_bamhi_detected(self):
        c = _import_constraints()
        checker = c.NoRestrictionSiteConstraint(sites=["GGATCC"])
        # GGA (Gly) + TCC (Ser) → GGATCC (BamHI)
        assert checker.is_violated(["GGA", "TCC"])

    def test_cross_codon_boundary(self):
        c = _import_constraints()
        checker = c.NoRestrictionSiteConstraint(sites=["AAGCTT"])  # HindIII
        # AAG (Lys) + CTT (Leu) → AAGCTT
        assert checker.is_violated(["AAG", "CTT"])


class TestGCRangeConstraint:

    def test_above_upper_bound(self):
        c = _import_constraints()
        checker = c.GCRangeConstraint(gc_lo=0.30, gc_hi=0.50)
        assert checker.is_violated(["GCG", "GCG"])  # 100% GC > 50%

    def test_below_lower_bound(self):
        c = _import_constraints()
        checker = c.GCRangeConstraint(gc_lo=0.50, gc_hi=0.70)
        assert checker.is_violated(["AAA", "AAA"])  # 0% GC < 50%

    def test_within_bounds(self):
        c = _import_constraints()
        checker = c.GCRangeConstraint(gc_lo=0.30, gc_hi=0.70)
        assert not checker.is_violated(["ATG", "GCG"])  # 4/6 ≈ 67%

    def test_exact_lower_bound(self):
        c = _import_constraints()
        checker = c.GCRangeConstraint(gc_lo=0.50, gc_hi=0.70)
        assert not checker.is_violated(["GCG", "AAA"])  # 3/6 = 50%

    def test_exact_upper_bound(self):
        c = _import_constraints()
        checker = c.GCRangeConstraint(gc_lo=0.30, gc_hi=0.50)
        assert not checker.is_violated(["GCG", "AAA"])  # 3/6 = 50%


class TestNoCrypticSpliceConstraint:

    def test_gt_codon_considered(self):
        """Constraint should evaluate without error for GT-containing codons."""
        c = _import_constraints()
        checker = c.NoCrypticSpliceConstraint(threshold=3.0)
        result = checker.is_violated(["GTT"])  # Valine: contains GT
        assert isinstance(result, bool)

    def test_no_gt_ag_passes(self):
        c = _import_constraints()
        checker = c.NoCrypticSpliceConstraint(threshold=3.0)
        assert not checker.is_violated(["AAA", "AAA"])

    def test_cross_codon_gt_evaluated(self):
        """GT at codon boundary (CAG|TTT) should be considered."""
        c = _import_constraints()
        checker = c.NoCrypticSpliceConstraint(threshold=3.0)
        result = checker.is_violated(["CAG", "TTT"])
        assert isinstance(result, bool)


class TestNoCpGIslandConstraint:

    def test_cg_within_codon_flagged(self):
        c = _import_constraints()
        checker = c.NoCpGIslandConstraint()
        assert checker.is_violated(["CGT", "CGC"])

    def test_no_cg_passes(self):
        c = _import_constraints()
        checker = c.NoCpGIslandConstraint()
        assert not checker.is_violated(["AAA", "AAA"])

    def test_cross_codon_cg_detected(self):
        """AAC|GTT → boundary CG."""
        c = _import_constraints()
        checker = c.NoCpGIslandConstraint()
        assert checker.is_violated(["AAC", "GTT"])


class TestNoATTTAMotifConstraint:

    def test_attta_detected(self):
        c = _import_constraints()
        checker = c.NoATTTAMotifConstraint()
        # AAT|TTA → "AATTTA" contains "ATTTA" at pos 1
        assert checker.is_violated(["AAT", "TTA"])

    def test_no_attta_passes(self):
        c = _import_constraints()
        checker = c.NoATTTAMotifConstraint()
        assert not checker.is_violated(["AAA", "CCC"])

    def test_attta_across_three_codons(self):
        c = _import_constraints()
        checker = c.NoATTTAMotifConstraint()
        # GAT|TTA|AAA → "GATTTAAAA" contains ATTTA at pos 1
        assert checker.is_violated(["GAT", "TTA", "AAA"])


class TestNoTRunConstraint:

    def test_six_ts_flagged(self):
        c = _import_constraints()
        checker = c.NoTRunConstraint(max_run=5)
        assert checker.is_violated(["TTT", "TTT"])  # 6 T's > 5

    def test_no_long_t_run_passes(self):
        c = _import_constraints()
        checker = c.NoTRunConstraint(max_run=5)
        assert not checker.is_violated(["AAA", "GCC"])

    def test_nine_ts_flagged(self):
        c = _import_constraints()
        checker = c.NoTRunConstraint(max_run=5)
        assert checker.is_violated(["TTT", "TTT", "TTT"])


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Helper Functions
# ═══════════════════════════════════════════════════════════════════════════════

class TestCodonGCCount:

    @pytest.mark.parametrize("codon,expected", [
        ("GCG", 3), ("ATA", 0), ("ATG", 1), ("TTC", 1), ("ACG", 2), ("GCC", 3),
    ])
    def test_gc_counts(self, codon, expected):
        c = _import_constraints()
        assert c.codon_gc_count(codon) == expected

    def test_all_gc_codons_return_three(self):
        c = _import_constraints()
        for codon in ["GGG", "GGC", "GCG", "GCC", "CGG", "CGC", "CCG", "CCC"]:
            assert c.codon_gc_count(codon) == 3

    def test_all_at_codons_return_zero(self):
        c = _import_constraints()
        for codon in ["AAA", "AAT", "ATA", "ATT", "TAA", "TAT", "TTA", "TTT"]:
            assert c.codon_gc_count(codon) == 0


class TestCodonContainsGT:

    @pytest.mark.parametrize("codon,expected", [
        ("GTT", True), ("ATG", False), ("AGT", True), ("GTG", True),
        ("TGG", False), ("TGT", True),
    ])
    def test_gt_detection(self, codon, expected):
        c = _import_constraints()
        assert c.codon_contains_gt(codon) is expected

    def test_all_valine_codons_contain_gt(self):
        c = _import_constraints()
        for codon in AA_TO_CODONS["V"]:
            assert c.codon_contains_gt(codon) is True


class TestCodonContainsCpG:

    @pytest.mark.parametrize("codon,expected", [
        ("CGT", True),   # CG at pos 0-1
        ("ACG", False),  # CG at pos 1-2, not canonical within-codon CpG start
        ("CGC", True),   # CG at pos 0-1
        ("CGG", True),   # CG at pos 0-1
        ("AAA", False),
    ])
    def test_cpg_detection(self, codon, expected):
        c = _import_constraints()
        assert c.codon_contains_cpg(codon) is expected


class TestComputeGCFromCodons:

    def test_gcg_ata_atg(self):
        """GCG(3) + ATA(0) + ATG(1) = 4/9."""
        c = _import_constraints()
        assert abs(c.compute_gc_from_codons(["GCG", "ATA", "ATG"]) - 4 / 9) < 1e-9

    def test_all_gc(self):
        c = _import_constraints()
        assert abs(c.compute_gc_from_codons(["GCG", "GCC", "CGC"]) - 1.0) < 1e-9

    def test_all_at(self):
        c = _import_constraints()
        assert abs(c.compute_gc_from_codons(["AAA", "ATA", "TTT"]) - 0.0) < 1e-9

    def test_single_codon(self):
        c = _import_constraints()
        assert abs(c.compute_gc_from_codons(["ATG"]) - 1 / 3) < 1e-9

    def test_empty_list(self):
        c = _import_constraints()
        assert c.compute_gc_from_codons([]) == 0.0

    def test_mixed_five_codons(self):
        """GCG(3)+ATA(0)+TTT(0)+GCC(3)+AAT(0) = 6/15 = 0.4."""
        c = _import_constraints()
        assert abs(c.compute_gc_from_codons(["GCG", "ATA", "TTT", "GCC", "AAT"]) - 0.4) < 1e-9


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Variable Domain Correctness
# ═══════════════════════════════════════════════════════════════════════════════

class TestVariableDomainCorrectness:

    def test_all_aas_domain_matches_constants(self):
        c = _import_constraints()
        for aa, expected in AA_TO_CODONS.items():
            assert set(c.get_domain_for_aa(aa)) == set(expected), f"AA={aa} mismatch"

    def test_stop_codons_not_in_any_domain(self):
        c = _import_constraints()
        for aa in AA_TO_CODONS:
            for codon in c.get_domain_for_aa(aa):
                assert codon not in STOP_CODONS, f"Stop codon {codon} in domain for {aa}"

    def test_methionine_single_codon(self):
        c = _import_constraints()
        assert set(c.get_domain_for_aa("M")) == {"ATG"}

    def test_tryptophan_single_codon(self):
        c = _import_constraints()
        assert set(c.get_domain_for_aa("W")) == {"TGG"}

    @pytest.mark.parametrize("aa,expected_count", [("L", 6), ("S", 6), ("R", 6)])
    def test_six_codon_aas(self, aa, expected_count):
        c = _import_constraints()
        assert len(c.get_domain_for_aa(aa)) == expected_count

    def test_all_domains_non_empty(self):
        c = _import_constraints()
        for aa in AA_TO_CODONS:
            assert len(c.get_domain_for_aa(aa)) > 0, f"AA={aa} has empty domain"

    def test_all_domain_codons_valid(self):
        """Every codon in every domain must be in CODON_TABLE mapping to the right AA."""
        c = _import_constraints()
        for aa in AA_TO_CODONS:
            for codon in c.get_domain_for_aa(aa):
                assert CODON_TABLE[codon] == aa, f"{codon} maps to {CODON_TABLE[codon]}, expected {aa}"

    def test_valine_all_contain_gt(self):
        """All Valine codons contain GT — critical for splice constraint reasoning."""
        c = _import_constraints()
        for codon in c.get_domain_for_aa("V"):
            assert "GT" in codon


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Edge Cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestEdgeCases:

    def test_single_aa_protein(self, default_config):
        c = _import_constraints()
        model = c.build_csp_model("M", "Homo_sapiens", default_config)
        assert len(model.variables) == 1
        assert set(model.variables[0].domain) == {"ATG"}

    def test_two_aa_protein(self, default_config):
        c = _import_constraints()
        model = c.build_csp_model("MA", "Homo_sapiens", default_config)
        assert len(model.variables) == 2

    def test_only_single_codon_aas(self, default_config):
        """Protein using only M and W should have singleton domains."""
        c = _import_constraints()
        model = c.build_csp_model("MWMWMW", "Homo_sapiens", default_config)
        for var in model.variables:
            assert len(var.domain) == 1

    def test_tight_gc_bounds(self, tight_gc_config):
        c = _import_constraints()
        protein = "ACDEFGHIKLMNPQRSTVWY"
        model = c.build_csp_model(protein, "Homo_sapiens", tight_gc_config)
        assert len(model.variables) == len(protein)
        gc_names = [c_.name for c_ in model.hard_constraints if "gc" in c_.name.lower()]
        assert len(gc_names) > 0

    def test_no_restriction_config_builds(self, no_restriction_config, sample_protein):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", no_restriction_config)
        assert len(model.variables) == len(sample_protein)
        rs_names = [c_.name for c_ in model.hard_constraints + model.soft_constraints
                    if "restriction" in c_.name.lower() or "rs_" in c_.name.lower()]
        assert len(rs_names) == 0

    def test_long_protein_100aa(self, default_config):
        c = _import_constraints()
        protein = "MVLSPADKTNVKAAWGKVGA" * 5  # 100 AAs
        model = c.build_csp_model(protein, "Homo_sapiens", default_config)
        assert len(model.variables) == 100

    def test_all_twenty_aas(self, default_config):
        c = _import_constraints()
        protein = "ACDEFGHIKLMNPQRSTVWY"
        model = c.build_csp_model(protein, "Homo_sapiens", default_config)
        assert len(model.variables) == 20
        for i, aa in enumerate(protein):
            assert set(model.variables[i].domain) == set(AA_TO_CODONS[aa])

    def test_gc_constraint_single_codon(self):
        c = _import_constraints()
        assert not c.GCRangeConstraint(0.0, 1.0).is_violated(["ATG"])

    def test_gc_constraint_impossible_bounds(self):
        c = _import_constraints()
        assert c.GCRangeConstraint(0.50, 0.70).is_violated(["AAA", "AAT", "ATA"])

    def test_cpg_constraint_single_codon(self):
        c = _import_constraints()
        assert c.NoCpGIslandConstraint().is_violated(["CGT"])

    def test_motif_constraint_single_codon(self):
        c = _import_constraints()
        assert not c.NoATTTAMotifConstraint().is_violated(["AAA"])

    def test_t_run_constraint_single_codon(self):
        c = _import_constraints()
        assert not c.NoTRunConstraint(max_run=5).is_violated(["TTT"])

    def test_specific_restriction_site_config(self):
        c = _import_constraints()
        config = {
            "organism": "Homo_sapiens", "gc_lo": 0.30, "gc_hi": 0.70,
            "restriction_sites": ["GAATTC"], "cryptic_splice_threshold": 3.0,
            "avoid_cpg": False, "avoid_attta": False, "max_t_run": 6,
        }
        model = c.build_csp_model("MVLSPADKTNVKAAWGKVGA", "Homo_sapiens", config)
        assert len(model.variables) == 20


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Model + Constraint Integration
# ═══════════════════════════════════════════════════════════════════════════════

class TestModelConstraintIntegration:

    def test_restriction_constraint_present(self, default_config):
        c = _import_constraints()
        model = c.build_csp_model("EF", "Homo_sapiens", default_config)  # GAA+TTC=GAATTC
        rs = [c_ for c_ in model.hard_constraints if "restriction" in c_.name.lower()]
        assert len(rs) > 0

    def test_gc_constraint_present(self, default_config, sample_protein):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        gc = [c_ for c_ in model.hard_constraints if "gc" in c_.name.lower()]
        assert len(gc) > 0

    def test_cpg_constraint_when_enabled(self, default_config, sample_protein):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        cpg = [c_ for c_ in model.soft_constraints + model.hard_constraints
               if "cpg" in c_.name.lower()]
        assert len(cpg) > 0

    def test_attta_constraint_when_enabled(self, default_config, sample_protein):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        attta = [c_ for c_ in model.soft_constraints + model.hard_constraints
                 if "attta" in c_.name.lower() or "instability" in c_.name.lower()]
        assert len(attta) > 0

    def test_t_run_constraint_registered(self, default_config, sample_protein):
        c = _import_constraints()
        model = c.build_csp_model(sample_protein, "Homo_sapiens", default_config)
        trun = [c_ for c_ in model.soft_constraints + model.hard_constraints
                if "t_run" in c_.name.lower() or "trun" in c_.name.lower()]
        assert len(trun) > 0


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Constants Consistency
# ═══════════════════════════════════════════════════════════════════════════════

class TestConstantsConsistency:

    def test_solver_domains_match_constants(self):
        c = _import_constraints()
        for aa, expected in AA_TO_CODONS.items():
            assert set(c.get_domain_for_aa(aa)) == set(expected)

    def test_no_stop_codons_in_domains(self):
        c = _import_constraints()
        for aa in AA_TO_CODONS:
            for codon in c.get_domain_for_aa(aa):
                assert codon not in STOP_CODONS

    def test_restriction_enzymes_usable(self):
        c = _import_constraints()
        checker = c.NoRestrictionSiteConstraint(sites=["GAATTC"])
        assert checker is not None

    def test_instability_motif_is_attta(self):
        assert INSTABILITY_MOTIF == "ATTTA"
