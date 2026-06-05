"""
Tests for extended gene sets in biocompiler.benchmarking.gene_sets.

Covers:
  1. All new gene sets have valid protein sequences (single-letter AA codes)
  2. Each gene entry has all required fields (name, protein, organism, expected_cai_range, category)
  3. Protein sequences contain only valid amino acid codes
  4. Gene set sizes are as expected
  5. Organism consistency within sets
  6. expected_cai_range tuples are well-formed
"""

import pytest

from biocompiler.benchmarking.gene_sets import (
    E_COLI_EXTENDED,
    HUMAN_SIGNALING,
    HUMAN_THERAPEUTIC,
    MOUSE_MODEL,
    YEAST_INDUSTRIAL,
    get_all_gene_sets,
)

# Standard 20 amino acid single-letter codes
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Required fields for extended gene set entries
REQUIRED_FIELDS = {"name", "protein", "organism", "expected_cai_range", "category"}

# All new gene sets for parametrized tests
NEW_GENE_SETS = {
    "E_COLI_EXTENDED": E_COLI_EXTENDED,
    "HUMAN_THERAPEUTIC": HUMAN_THERAPEUTIC,
    "HUMAN_SIGNALING": HUMAN_SIGNALING,
    "YEAST_INDUSTRIAL": YEAST_INDUSTRIAL,
    "MOUSE_MODEL": MOUSE_MODEL,
}


# ---------------------------------------------------------------------------
# 1. Gene set size tests
# ---------------------------------------------------------------------------

class TestGeneSetSizes:
    """Tests that each gene set has the expected number of entries."""

    def test_ecoli_extended_has_20_genes(self):
        assert len(E_COLI_EXTENDED) == 20

    def test_human_therapeutic_has_5_genes(self):
        assert len(HUMAN_THERAPEUTIC) == 5

    def test_human_signaling_has_4_genes(self):
        assert len(HUMAN_SIGNALING) == 4

    def test_yeast_industrial_has_3_genes(self):
        assert len(YEAST_INDUSTRIAL) == 3

    def test_mouse_model_has_6_genes(self):
        assert len(MOUSE_MODEL) == 6


# ---------------------------------------------------------------------------
# 2. Required fields tests
# ---------------------------------------------------------------------------

class TestRequiredFields:
    """Tests that every gene entry has all required fields."""

    @pytest.mark.parametrize(
        "set_name,gene_set",
        list(NEW_GENE_SETS.items()),
    )
    def test_all_required_fields_present(self, set_name, gene_set):
        """Each gene entry must contain all required fields."""
        for gene_key, entry in gene_set.items():
            missing = REQUIRED_FIELDS - set(entry.keys())
            assert not missing, (
                f"{set_name}[{gene_key}] missing fields: {missing}"
            )

    @pytest.mark.parametrize(
        "set_name,gene_set",
        list(NEW_GENE_SETS.items()),
    )
    def test_name_matches_key(self, set_name, gene_set):
        """The 'name' field should match the dictionary key."""
        for gene_key, entry in gene_set.items():
            assert entry["name"] == gene_key, (
                f"{set_name}[{gene_key}]: name field '{entry['name']}' != key '{gene_key}'"
            )


# ---------------------------------------------------------------------------
# 3. Valid protein sequence tests
# ---------------------------------------------------------------------------

class TestValidProteinSequences:
    """Tests that all protein sequences contain only valid amino acid codes."""

    @pytest.mark.parametrize(
        "set_name,gene_set",
        list(NEW_GENE_SETS.items()),
    )
    def test_proteins_contain_only_valid_aa(self, set_name, gene_set):
        """Every character in every protein sequence must be a valid AA code."""
        for gene_key, entry in gene_set.items():
            seq = entry["protein"]
            invalid = set(seq) - VALID_AA
            assert not invalid, (
                f"{set_name}[{gene_key}] has invalid amino acids: {invalid}"
            )

    @pytest.mark.parametrize(
        "set_name,gene_set",
        list(NEW_GENE_SETS.items()),
    )
    def test_proteins_are_non_empty(self, set_name, gene_set):
        """No protein sequence should be empty."""
        for gene_key, entry in gene_set.items():
            assert len(entry["protein"]) > 0, (
                f"{set_name}[{gene_key}] has empty protein sequence"
            )

    @pytest.mark.parametrize(
        "set_name,gene_set",
        list(NEW_GENE_SETS.items()),
    )
    def test_protein_is_string(self, set_name, gene_set):
        """The protein field must be a string."""
        for gene_key, entry in gene_set.items():
            assert isinstance(entry["protein"], str), (
                f"{set_name}[{gene_key}] protein is not a string: {type(entry['protein'])}"
            )


# ---------------------------------------------------------------------------
# 4. Organism consistency tests
# ---------------------------------------------------------------------------

class TestOrganismConsistency:
    """Tests that organisms are consistent and correctly specified."""

    def test_ecoli_organisms(self):
        """All E. coli genes should have Escherichia coli organism."""
        for gene_key, entry in E_COLI_EXTENDED.items():
            assert entry["organism"] == "Escherichia coli", (
                f"lacZ[{gene_key}] organism is '{entry['organism']}'"
            )

    def test_human_therapeutic_organisms(self):
        """All human therapeutic genes should be from Homo sapiens."""
        for gene_key, entry in HUMAN_THERAPEUTIC.items():
            assert entry["organism"] == "Homo sapiens", (
                f"HUMAN_THERAPEUTIC[{gene_key}] organism is '{entry['organism']}'"
            )

    def test_human_signaling_organisms(self):
        """All human signaling genes should be from Homo sapiens."""
        for gene_key, entry in HUMAN_SIGNALING.items():
            assert entry["organism"] == "Homo sapiens", (
                f"HUMAN_SIGNALING[{gene_key}] organism is '{entry['organism']}'"
            )

    def test_yeast_industrial_organisms(self):
        """All yeast industrial genes should be from S. cerevisiae."""
        for gene_key, entry in YEAST_INDUSTRIAL.items():
            assert entry["organism"] == "Saccharomyces cerevisiae", (
                f"YEAST_INDUSTRIAL[{gene_key}] organism is '{entry['organism']}'"
            )

    def test_mouse_model_organisms(self):
        """All mouse model genes should be from Mus musculus."""
        for gene_key, entry in MOUSE_MODEL.items():
            assert entry["organism"] == "Mus musculus", (
                f"MOUSE_MODEL[{gene_key}] organism is '{entry['organism']}'"
            )


# ---------------------------------------------------------------------------
# 5. expected_cai_range tests
# ---------------------------------------------------------------------------

class TestExpectedCaiRange:
    """Tests that expected_cai_range tuples are well-formed."""

    @pytest.mark.parametrize(
        "set_name,gene_set",
        list(NEW_GENE_SETS.items()),
    )
    def test_cai_range_is_tuple_of_two(self, set_name, gene_set):
        """expected_cai_range should be a tuple of exactly 2 floats."""
        for gene_key, entry in gene_set.items():
            cai_range = entry["expected_cai_range"]
            assert isinstance(cai_range, tuple), (
                f"{set_name}[{gene_key}] expected_cai_range is not a tuple"
            )
            assert len(cai_range) == 2, (
                f"{set_name}[{gene_key}] expected_cai_range length != 2"
            )

    @pytest.mark.parametrize(
        "set_name,gene_set",
        list(NEW_GENE_SETS.items()),
    )
    def test_cai_range_values_in_bounds(self, set_name, gene_set):
        """expected_cai_range values should be in [0, 1]."""
        for gene_key, entry in gene_set.items():
            lo, hi = entry["expected_cai_range"]
            assert 0.0 <= lo <= 1.0, (
                f"{set_name}[{gene_key}] CAI range low {lo} out of [0,1]"
            )
            assert 0.0 <= hi <= 1.0, (
                f"{set_name}[{gene_key}] CAI range high {hi} out of [0,1]"
            )
            assert lo < hi, (
                f"{set_name}[{gene_key}] CAI range low {lo} >= high {hi}"
            )


# ---------------------------------------------------------------------------
# 6. Category tests
# ---------------------------------------------------------------------------

class TestCategories:
    """Tests that categories are valid and consistent."""

    def test_ecoli_extended_categories(self):
        """E. coli genes should have known functional categories."""
        valid_categories = {
            "housekeeping", "regulatory", "dna_repair", "dna_replication",
            "transcription", "chaperone", "translation", "resistance",
            "transport", "membrane", "metabolic",
        }
        for gene_key, entry in E_COLI_EXTENDED.items():
            assert entry["category"] in valid_categories, (
                f"E_COLI_EXTENDED[{gene_key}] unknown category: {entry['category']}"
            )

    def test_human_therapeutic_category(self):
        """All human therapeutic genes should have 'therapeutic' category."""
        for gene_key, entry in HUMAN_THERAPEUTIC.items():
            assert entry["category"] == "therapeutic", (
                f"HUMAN_THERAPEUTIC[{gene_key}] category is '{entry['category']}'"
            )

    def test_human_signaling_category(self):
        """All human signaling genes should have 'signaling' category."""
        for gene_key, entry in HUMAN_SIGNALING.items():
            assert entry["category"] == "signaling", (
                f"HUMAN_SIGNALING[{gene_key}] category is '{entry['category']}'"
            )

    def test_yeast_industrial_category(self):
        """All yeast industrial genes should have 'industrial' category."""
        for gene_key, entry in YEAST_INDUSTRIAL.items():
            assert entry["category"] == "industrial", (
                f"YEAST_INDUSTRIAL[{gene_key}] category is '{entry['category']}'"
            )

    def test_mouse_model_category(self):
        """All mouse model genes should have 'model_organism' category."""
        for gene_key, entry in MOUSE_MODEL.items():
            assert entry["category"] == "model_organism", (
                f"MOUSE_MODEL[{gene_key}] category is '{entry['category']}'"
            )


# ---------------------------------------------------------------------------
# 7. Specific gene presence tests
# ---------------------------------------------------------------------------

class TestSpecificGenePresence:
    """Tests that specific required genes are present in each set."""

    def test_ecoli_has_lacz(self):
        assert "lacZ" in E_COLI_EXTENDED

    def test_ecoli_has_tufa(self):
        assert "tufA" in E_COLI_EXTENDED

    def test_ecoli_has_groel(self):
        assert "groEL" in E_COLI_EXTENDED

    def test_human_therapeutic_has_insulin(self):
        assert "INS" in HUMAN_THERAPEUTIC

    def test_human_therapeutic_has_hgh(self):
        assert "GH1" in HUMAN_THERAPEUTIC

    def test_human_therapeutic_has_ifna2(self):
        assert "IFNA2" in HUMAN_THERAPEUTIC

    def test_human_therapeutic_has_epo(self):
        assert "EPO" in HUMAN_THERAPEUTIC

    def test_human_therapeutic_has_gcsf(self):
        assert "CSF3" in HUMAN_THERAPEUTIC

    def test_human_signaling_has_egfr(self):
        assert "EGFR" in HUMAN_SIGNALING

    def test_human_signaling_has_her2(self):
        assert "ERBB2" in HUMAN_SIGNALING

    def test_human_signaling_has_vegf(self):
        assert "VEGFA" in HUMAN_SIGNALING

    def test_human_signaling_has_pd1(self):
        assert "PDCD1" in HUMAN_SIGNALING

    def test_yeast_has_invertase(self):
        assert "SUC2" in YEAST_INDUSTRIAL

    def test_yeast_has_amylase(self):
        assert "AMY1" in YEAST_INDUSTRIAL

    def test_yeast_has_cellulase(self):
        assert "CEL3A" in YEAST_INDUSTRIAL

    def test_mouse_has_albumin(self):
        assert "ALB_MOUSE" in MOUSE_MODEL

    def test_mouse_has_p53(self):
        assert "TP53_MOUSE" in MOUSE_MODEL

    def test_mouse_has_hbb(self):
        assert "HBB_MOUSE" in MOUSE_MODEL


# ---------------------------------------------------------------------------
# 8. get_all_gene_sets integration test
# ---------------------------------------------------------------------------

class TestGetAllGeneSets:
    """Tests for the get_all_gene_sets aggregation function."""

    def test_merged_contains_new_sets(self):
        """The merged dict should contain genes from all new sets."""
        merged = get_all_gene_sets()
        # Check at least one gene from each new set
        assert "lacZ" in merged
        assert "INS" in merged
        assert "EGFR" in merged
        assert "SUC2" in merged
        assert "ALB_MOUSE" in merged

    def test_merged_total_count(self):
        """The merged dict should have genes from all sets combined."""
        merged = get_all_gene_sets()
        # Legacy: ~26 genes, Extended: 20+5+4+3+6 = 38 genes
        # Some overlap in keys (e.g. INS exists in both legacy and new)
        assert len(merged) >= 38
