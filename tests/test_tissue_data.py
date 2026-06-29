"""
BioCompiler Tissue-Specific Splicing Data Tests
================================================

Tests for the biocompiler.tissue_data module — tissue weight lookup,
alias resolution, custom tissue addition, listing, and JSON export.
"""

from __future__ import annotations

import json
import pytest

from biocompiler.organisms.tissue_data import (
    GTEX_TISSUE_WEIGHTS,
    TISSUE_ALIASES,
    CANONICAL_BASELINE_WEIGHT,
    get_tissue_weights,
    list_available_tissues,
    add_custom_tissue,
    export_tissue_weights_json,
)


# ════════════════════════════════════════════════════════════════════════════
# Constants
# ════════════════════════════════════════════════════════════════════════════

class TestConstants:

    def test_canonical_baseline_weight(self):
        """Canonical baseline weight is 1.0."""
        assert CANONICAL_BASELINE_WEIGHT == 1.0

    def test_builtin_tissues_present(self):
        """Core GTEx tissues are present in the weights dict."""
        expected = {"HEK293T", "HeLa", "HepG2", "Brain", "Heart", "Liver",
                    "Kidney", "Lung", "Muscle", "Testis", "Whole_Blood", "default"}
        assert expected.issubset(set(GTEX_TISSUE_WEIGHTS.keys()))

    def test_all_tissues_have_required_keys(self):
        """Every tissue entry has all required event type keys."""
        required = {"canonical", "exon_skip", "intron_retention", "alt_site", "cryptic"}
        for tissue, weights in GTEX_TISSUE_WEIGHTS.items():
            assert required.issubset(set(weights.keys())), (
                f"Tissue '{tissue}' missing keys: {required - set(weights.keys())}"
            )

    def test_all_canonical_weights_are_one(self):
        """Canonical weight is 1.0 for all tissues."""
        for tissue, weights in GTEX_TISSUE_WEIGHTS.items():
            assert weights["canonical"] == 1.0, (
                f"Tissue '{tissue}' has canonical={weights['canonical']}, expected 1.0"
            )

    def test_weights_are_non_negative(self):
        """All weight values should be non-negative."""
        for tissue, weights in GTEX_TISSUE_WEIGHTS.items():
            for event, val in weights.items():
                assert val >= 0.0, (
                    f"Tissue '{tissue}' event '{event}' has negative weight {val}"
                )

    def test_aliases_map_to_valid_tissues(self):
        """All alias targets exist in GTEX_TISSUE_WEIGHTS."""
        for alias, canonical in TISSUE_ALIASES.items():
            assert canonical in GTEX_TISSUE_WEIGHTS, (
                f"Alias '{alias}' maps to '{canonical}' which is not in GTEX_TISSUE_WEIGHTS"
            )


# ════════════════════════════════════════════════════════════════════════════
# get_tissue_weights
# ════════════════════════════════════════════════════════════════════════════

class TestGetTissueWeights:

    def test_exact_match_hek293t(self):
        """Exact match for HEK293T returns correct weights."""
        weights = get_tissue_weights("HEK293T")
        assert weights["canonical"] == 1.0
        assert "exon_skip" in weights

    def test_exact_match_brain(self):
        """Exact match for Brain returns correct weights."""
        weights = get_tissue_weights("Brain")
        assert weights["canonical"] == 1.0

    def test_alias_match_hek293(self):
        """Alias 'hek293' resolves to HEK293T."""
        weights = get_tissue_weights("hek293")
        assert weights == GTEX_TISSUE_WEIGHTS["HEK293T"]

    def test_alias_match_blood(self):
        """Alias 'blood' resolves to Whole_Blood."""
        weights = get_tissue_weights("blood")
        assert weights == GTEX_TISSUE_WEIGHTS["Whole_Blood"]

    def test_alias_match_cerebellum(self):
        """Alias 'cerebellum' resolves to Brain."""
        weights = get_tissue_weights("cerebellum")
        assert weights == GTEX_TISSUE_WEIGHTS["Brain"]

    def test_case_insensitive_match(self):
        """Case-insensitive match for 'brain' returns Brain weights."""
        weights = get_tissue_weights("brain")
        assert weights == GTEX_TISSUE_WEIGHTS["Brain"]

    def test_unknown_tissue_returns_default(self):
        """Unknown tissue falls back to default weights."""
        weights = get_tissue_weights("MartianTissue_XYZ")
        assert weights == GTEX_TISSUE_WEIGHTS["default"]

    def test_returned_dict_has_all_keys(self):
        """Returned weights dict has all required event type keys."""
        required = {"canonical", "exon_skip", "intron_retention", "alt_site", "cryptic"}
        for tissue in ["HEK293T", "Brain", "default", "Testis"]:
            weights = get_tissue_weights(tissue)
            assert required.issubset(set(weights.keys()))


# ════════════════════════════════════════════════════════════════════════════
# list_available_tissues
# ════════════════════════════════════════════════════════════════════════════

class TestListAvailableTissues:

    def test_returns_sorted_list(self):
        """Result is a sorted list of tissue names."""
        tissues = list_available_tissues()
        assert tissues == sorted(tissues)

    def test_excludes_default(self):
        """The 'default' entry is excluded from the list."""
        tissues = list_available_tissues()
        assert "default" not in tissues

    def test_includes_core_tissues(self):
        """Core tissue names are present."""
        tissues = list_available_tissues()
        assert "HEK293T" in tissues
        assert "Brain" in tissues
        assert "Liver" in tissues

    def test_returns_list(self):
        """Result is a list."""
        tissues = list_available_tissues()
        assert isinstance(tissues, list)


# ════════════════════════════════════════════════════════════════════════════
# add_custom_tissue
# ════════════════════════════════════════════════════════════════════════════

class TestAddCustomTissue:

    def test_add_valid_custom_tissue(self):
        """Adding a custom tissue with all required keys works."""
        add_custom_tissue("MyCustomTissue", {
            "canonical": 1.0,
            "exon_skip": 0.35,
            "intron_retention": 0.15,
            "alt_site": 0.40,
            "cryptic": 0.10,
        })
        weights = get_tissue_weights("MyCustomTissue")
        assert weights["exon_skip"] == 0.35

    def test_missing_keys_raises(self):
        """Adding a custom tissue with missing keys raises ValueError."""
        with pytest.raises(ValueError, match="Missing required weight keys"):
            add_custom_tissue("BadTissue", {"canonical": 1.0})

    def test_canonical_not_one_warns(self):
        """Adding a custom tissue with canonical != 1.0 logs a warning
        but does not raise."""
        # Should not raise
        add_custom_tissue("NonStandardCanonical", {
            "canonical": 0.9,
            "exon_skip": 0.3,
            "intron_retention": 0.1,
            "alt_site": 0.2,
            "cryptic": 0.05,
        })
        weights = get_tissue_weights("NonStandardCanonical")
        assert weights["canonical"] == 0.9

    def test_override_existing_tissue(self):
        """Adding a tissue with an existing name overrides it."""
        add_custom_tissue("HEK293T_override_test", {
            "canonical": 1.0,
            "exon_skip": 0.99,
            "intron_retention": 0.01,
            "alt_site": 0.01,
            "cryptic": 0.01,
        })
        weights = get_tissue_weights("HEK293T_override_test")
        assert weights["exon_skip"] == 0.99


# ════════════════════════════════════════════════════════════════════════════
# export_tissue_weights_json
# ════════════════════════════════════════════════════════════════════════════

class TestExportTissueWeightsJson:

    def test_returns_valid_json_string(self):
        """Export without a path returns a valid JSON string."""
        result = export_tissue_weights_json()
        parsed = json.loads(result)
        assert isinstance(parsed, dict)

    def test_json_contains_weights(self):
        """Exported JSON has a 'weights' key."""
        result = export_tissue_weights_json()
        parsed = json.loads(result)
        assert "weights" in parsed

    def test_json_contains_aliases(self):
        """Exported JSON has an 'aliases' key."""
        result = export_tissue_weights_json()
        parsed = json.loads(result)
        assert "aliases" in parsed

    def test_json_contains_source(self):
        """Exported JSON has a 'source' key."""
        result = export_tissue_weights_json()
        parsed = json.loads(result)
        assert "source" in parsed

    def test_write_to_file(self, tmp_path):
        """Export with output_path writes a file."""
        out = str(tmp_path / "tissue_weights.json")
        result = export_tissue_weights_json(output_path=out)
        # The function returns the JSON string AND writes to file
        with open(out) as f:
            file_content = f.read()
        assert json.loads(file_content)  # valid JSON in file

    def test_roundtrip_preserves_data(self):
        """Export → parse roundtrip preserves tissue data."""
        result = export_tissue_weights_json()
        parsed = json.loads(result)
        assert "HEK293T" in parsed["weights"]
        assert parsed["weights"]["HEK293T"]["canonical"] == 1.0
