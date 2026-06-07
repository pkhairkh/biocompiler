"""Test BioCompiler Grammar Loader — YAML configuration loading and validation."""

import textwrap
from pathlib import Path

import pytest
import yaml

from biocompiler.grammar_loader import (
    _DEFAULT_CAI_THRESHOLD,
    _DEFAULT_CELLULAR_CONTEXT,
    _DEFAULT_CRYPTIC_SPLICE_THRESHOLD,
    _DEFAULT_ENZYMES,
    _DEFAULT_EXON_BOUNDARIES,
    _DEFAULT_GC_HI,
    _DEFAULT_GC_LO,
    _DEFAULT_ORGANISM,
    _DEFAULT_UNCERTAIN_LO,
    grammar_to_predicate_params,
    list_builtin_grammars,
    load_builtin_grammar,
    load_grammar,
)


# ────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────

def _write_grammar(tmp_path: Path, content: str, name: str = "test_grammar.yaml") -> Path:
    """Write a grammar YAML string to a temp file and return the path."""
    p = tmp_path / name
    p.write_text(textwrap.dedent(content), encoding="utf-8")
    return p


def _minimal_grammar() -> dict:
    """Return a minimally valid grammar dict (for unit-test use)."""
    return {
        "gene": {"name": "TestGene", "organism": "Homo_sapiens"},
        "type_system": {
            "predicates": [
                {"name": "SpliceCorrect"},
                {"name": "GCInRange", "lo": 0.30, "hi": 0.70},
            ]
        },
    }


# ────────────────────────────────────────────────────────────
# 1. Grammar loading from YAML files
# ────────────────────────────────────────────────────────────

class TestLoadGrammar:
    """Tests for load_grammar() — loading grammar from a YAML file path."""

    def test_load_valid_grammar(self, tmp_path):
        """A well-formed grammar YAML file loads without error."""
        content = """\
        gene:
          name: "TestGene"
          organism: "Homo_sapiens"
          expression_organism: "Mus_musculus"
        type_system:
          predicates:
            - name: "SpliceCorrect"
            - name: "GCInRange"
              lo: 0.35
              hi: 0.65
        """
        path = _write_grammar(tmp_path, content)
        grammar = load_grammar(path)
        assert isinstance(grammar, dict)
        assert grammar["gene"]["name"] == "TestGene"
        assert grammar["type_system"]["predicates"][0]["name"] == "SpliceCorrect"

    def test_load_grammar_with_all_sections(self, tmp_path):
        """Grammar with gene, pre_mrna, exons, introns, ndfst, type_system loads."""
        content = """\
        gene:
          name: "HBB"
          organism: "Homo_sapiens"
        pre_mrna:
          length_bp: 1608
        exons:
          - id: "exon1"
            start: 1
            end: 142
          - id: "exon2"
            start: 273
            end: 495
        introns:
          - id: "IVS1"
            start: 143
            end: 272
        ndfst:
          cell_type: "HEK293T"
          splice_consensus:
            donor: "GT"
            acceptor: "AG"
        type_system:
          predicates:
            - name: "SpliceCorrect"
            - name: "NoCrypticSplice"
              cryptic_threshold: 3.0
            - name: "CodonAdapted"
              threshold: 0.70
            - name: "GCInRange"
              lo: 0.40
              hi: 0.60
            - name: "NoRestrictionSite"
              enzyme_sites:
                - "EcoRI"
                - "BamHI"
            - name: "InFrame"
            - name: "NoInstabilityMotif"
            - name: "NoCpGIsland"
        """
        path = _write_grammar(tmp_path, content)
        grammar = load_grammar(path)
        assert len(grammar["exons"]) == 2
        assert len(grammar["introns"]) == 1
        assert grammar["ndfst"]["cell_type"] == "HEK293T"
        assert len(grammar["type_system"]["predicates"]) == 8

    def test_load_grammar_accepts_pathlib_path(self, tmp_path):
        """load_grammar() accepts a pathlib.Path argument."""
        content = """\
        gene:
          name: "GFP"
          organism: "Aequorea_victoria"
        type_system:
          predicates:
            - name: "SpliceCorrect"
        """
        path = _write_grammar(tmp_path, content)
        result = load_grammar(path)  # Already a Path object
        assert result["gene"]["name"] == "GFP"

    def test_load_grammar_accepts_string_path(self, tmp_path):
        """load_grammar() accepts a plain string path."""
        content = """\
        gene:
          name: "GFP"
          organism: "Aequorea_victoria"
        type_system:
          predicates:
            - name: "SpliceCorrect"
        """
        path = _write_grammar(tmp_path, content)
        result = load_grammar(str(path))  # Pass as string
        assert result["gene"]["name"] == "GFP"

    def test_load_grammar_preserves_predicate_params(self, tmp_path):
        """Predicate-specific parameters (thresholds, lo/hi) are preserved."""
        content = """\
        gene:
          name: "Test"
          organism: "Homo_sapiens"
        type_system:
          predicates:
            - name: "GCInRange"
              lo: 0.45
              hi: 0.55
            - name: "CodonAdapted"
              threshold: 0.80
            - name: "NoCrypticSplice"
              cryptic_threshold: 5.0
              uncertain_lo: 2.0
        """
        path = _write_grammar(tmp_path, content)
        grammar = load_grammar(path)
        preds = grammar["type_system"]["predicates"]
        gc = next(p for p in preds if p["name"] == "GCInRange")
        assert gc["lo"] == 0.45
        assert gc["hi"] == 0.55
        cai = next(p for p in preds if p["name"] == "CodonAdapted")
        assert cai["threshold"] == 0.80


# ────────────────────────────────────────────────────────────
# 2. Built-in grammar loading
# ────────────────────────────────────────────────────────────

class TestBuiltinGrammars:
    """Tests for load_builtin_grammar() and list_builtin_grammars()."""

    def test_list_builtin_grammars_returns_list(self):
        """list_builtin_grammars() returns a list of strings."""
        result = list_builtin_grammars()
        assert isinstance(result, list)
        for name in result:
            assert isinstance(name, str)

    def test_list_builtin_grammars_includes_known(self):
        """Built-in grammars include the shipped YAML files."""
        names = list_builtin_grammars()
        # The repo ships egfp_hek293t and hbb_hek293t
        assert "egfp_hek293t" in names
        assert "hbb_hek293t" in names

    def test_load_builtin_egfp(self):
        """Load the built-in EGFP grammar by name."""
        grammar = load_builtin_grammar("egfp_hek293t")
        assert grammar["gene"]["name"] == "EGFP"
        assert "type_system" in grammar
        assert "predicates" in grammar["type_system"]

    def test_load_builtin_hbb(self):
        """Load the built-in HBB grammar by name."""
        grammar = load_builtin_grammar("hbb_hek293t")
        assert grammar["gene"]["name"] == "HBB"
        assert grammar["gene"]["organism"] == "Homo_sapiens"

    def test_load_builtin_unknown_raises(self):
        """Loading a non-existent built-in grammar raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="No built-in grammar named"):
            load_builtin_grammar("nonexistent_gene")

    def test_load_builtin_unknown_mentions_available(self):
        """Error message for unknown grammar mentions available grammars."""
        with pytest.raises(FileNotFoundError, match="Available grammars"):
            load_builtin_grammar("does_not_exist")


# ────────────────────────────────────────────────────────────
# 3. Invalid grammar handling
# ────────────────────────────────────────────────────────────

class TestInvalidGrammar:
    """Tests for error handling when grammars are malformed."""

    def test_file_not_found(self):
        """Non-existent file path raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="Grammar file not found"):
            load_grammar("/nonexistent/path/grammar.yaml")

    def test_empty_yaml_file(self, tmp_path):
        """An empty YAML file (loads as None) raises ValueError."""
        path = _write_grammar(tmp_path, "")
        with pytest.raises(ValueError, match="YAML mapping"):
            load_grammar(path)

    def test_yaml_list_instead_of_dict(self, tmp_path):
        """A YAML file containing a list raises ValueError."""
        content = """\
        - item1
        - item2
        """
        path = _write_grammar(tmp_path, content)
        with pytest.raises(ValueError, match="YAML mapping"):
            load_grammar(path)

    def test_yaml_scalar_instead_of_dict(self, tmp_path):
        """A YAML file containing a scalar raises ValueError."""
        path = _write_grammar(tmp_path, "just a string\n")
        with pytest.raises(ValueError, match="YAML mapping"):
            load_grammar(path)

    def test_missing_gene_section(self, tmp_path):
        """Grammar without a 'gene' key raises ValueError."""
        content = """\
        type_system:
          predicates:
            - name: "SpliceCorrect"
        """
        path = _write_grammar(tmp_path, content)
        with pytest.raises(ValueError, match="Missing 'gene' section"):
            load_grammar(path)

    def test_missing_type_system_section(self, tmp_path):
        """Grammar without a 'type_system' key raises ValueError."""
        content = """\
        gene:
          name: "TestGene"
          organism: "Homo_sapiens"
        """
        path = _write_grammar(tmp_path, content)
        with pytest.raises(ValueError, match="Missing 'type_system' section"):
            load_grammar(path)

    def test_gene_not_a_mapping(self, tmp_path):
        """Gene section that is not a dict raises ValueError."""
        content = """\
        gene: "not_a_dict"
        type_system:
          predicates:
            - name: "SpliceCorrect"
        """
        path = _write_grammar(tmp_path, content)
        with pytest.raises(ValueError, match="'gene' must be a mapping"):
            load_grammar(path)

    def test_type_system_not_a_mapping(self, tmp_path):
        """type_system section that is not a dict raises ValueError."""
        content = """\
        gene:
          name: "TestGene"
        type_system: "not_a_dict"
        """
        path = _write_grammar(tmp_path, content)
        with pytest.raises(ValueError, match="'type_system' must be a mapping"):
            load_grammar(path)

    def test_predicates_not_a_list(self, tmp_path):
        """type_system.predicates that is not a list raises ValueError."""
        content = """\
        gene:
          name: "TestGene"
        type_system:
          predicates: "not_a_list"
        """
        path = _write_grammar(tmp_path, content)
        with pytest.raises(ValueError, match="'type_system.predicates' must be a list"):
            load_grammar(path)

    def test_predicate_not_a_mapping(self, tmp_path):
        """A predicate that is not a dict raises ValueError."""
        content = """\
        gene:
          name: "TestGene"
        type_system:
          predicates:
            - "just_a_string"
        """
        path = _write_grammar(tmp_path, content)
        with pytest.raises(ValueError, match="Predicate 0 must be a mapping"):
            load_grammar(path)

    def test_predicate_missing_name(self, tmp_path):
        """A predicate without a 'name' key raises ValueError."""
        content = """\
        gene:
          name: "TestGene"
        type_system:
          predicates:
            - threshold: 0.5
        """
        path = _write_grammar(tmp_path, content)
        with pytest.raises(ValueError, match="Predicate 0 missing 'name' key"):
            load_grammar(path)

    def test_unknown_predicate_name_warns(self, tmp_path, caplog):
        """An unrecognized predicate name logs a warning but does not raise."""
        content = """\
        gene:
          name: "TestGene"
        type_system:
          predicates:
            - name: "SpliceCorrect"
            - name: "CustomPredicate"
        """
        path = _write_grammar(tmp_path, content)
        import logging
        with caplog.at_level(logging.WARNING, logger="biocompiler.grammar_loader"):
            grammar = load_grammar(path)
        assert grammar is not None  # Should succeed (just warns)
        assert "CustomPredicate" in caplog.text

    def test_multiple_validation_errors(self, tmp_path):
        """Multiple validation issues are detected and all reported."""
        content = """\
        gene: "not_a_dict"
        type_system: "also_not_a_dict"
        """
        path = _write_grammar(tmp_path, content)
        with pytest.raises(ValueError) as exc_info:
            load_grammar(path)
        msg = str(exc_info.value)
        assert "'gene' must be a mapping" in msg
        assert "'type_system' must be a mapping" in msg

    def test_missing_gene_and_type_system_both_reported(self, tmp_path):
        """When both gene and type_system are missing, both errors appear."""
        content = """\
        ndfst:
          cell_type: "HEK293T"
        """
        path = _write_grammar(tmp_path, content)
        with pytest.raises(ValueError) as exc_info:
            load_grammar(path)
        msg = str(exc_info.value)
        assert "Missing 'gene' section" in msg
        assert "Missing 'type_system' section" in msg


# ────────────────────────────────────────────────────────────
# 4. Default grammar behavior (grammar_to_predicate_params)
# ────────────────────────────────────────────────────────────

class TestGrammarToPredicateParams:
    """Tests for grammar_to_predicate_params() — extracting parameters with defaults."""

    def test_minimal_grammar_gets_defaults(self):
        """A grammar with no predicate params returns all defaults."""
        grammar = _minimal_grammar()
        params = grammar_to_predicate_params(grammar)
        assert params["organism"] == "Homo_sapiens"
        assert params["gc_lo"] == _DEFAULT_GC_LO
        assert params["gc_hi"] == _DEFAULT_GC_HI
        assert params["cai_threshold"] == _DEFAULT_CAI_THRESHOLD
        assert params["cryptic_splice_threshold"] == _DEFAULT_CRYPTIC_SPLICE_THRESHOLD
        assert params["uncertain_lo"] == _DEFAULT_UNCERTAIN_LO
        assert params["enzymes"] == list(_DEFAULT_ENZYMES)
        assert params["exon_boundaries"] == list(_DEFAULT_EXON_BOUNDARIES)
        assert params["cellular_context"] == _DEFAULT_CELLULAR_CONTEXT

    def test_expression_organism_overrides_organism(self):
        """expression_organism takes precedence over organism in gene section."""
        grammar = {
            "gene": {"name": "GFP", "organism": "Aequorea_victoria", "expression_organism": "Mus_musculus"},
            "type_system": {"predicates": []},
        }
        params = grammar_to_predicate_params(grammar)
        assert params["organism"] == "Mus_musculus"

    def test_organism_fallback_when_no_expression_organism(self):
        """When expression_organism is absent, organism is used."""
        grammar = {
            "gene": {"name": "GFP", "organism": "Aequorea_victoria"},
            "type_system": {"predicates": []},
        }
        params = grammar_to_predicate_params(grammar)
        assert params["organism"] == "Aequorea_victoria"

    def test_default_organism_when_both_absent(self):
        """When both organism and expression_organism are absent, default is used."""
        grammar = {"gene": {"name": "GFP"}, "type_system": {"predicates": []}}
        params = grammar_to_predicate_params(grammar)
        assert params["organism"] == _DEFAULT_ORGANISM

    def test_gc_in_range_params(self):
        """GCInRange predicate params are extracted from grammar."""
        grammar = {
            "gene": {"name": "T"},
            "type_system": {"predicates": [{"name": "GCInRange", "lo": 0.40, "hi": 0.60}]},
        }
        params = grammar_to_predicate_params(grammar)
        assert params["gc_lo"] == 0.40
        assert params["gc_hi"] == 0.60

    def test_codon_adapted_params(self):
        """CodonAdapted predicate params are extracted from grammar."""
        grammar = {
            "gene": {"name": "T"},
            "type_system": {"predicates": [{"name": "CodonAdapted", "threshold": 0.80}]},
        }
        params = grammar_to_predicate_params(grammar)
        assert params["cai_threshold"] == 0.80

    def test_no_cryptic_splice_params(self):
        """NoCrypticSplice predicate params are extracted from grammar."""
        grammar = {
            "gene": {"name": "T"},
            "type_system": {
                "predicates": [
                    {"name": "NoCrypticSplice", "cryptic_threshold": 5.0, "uncertain_lo": 2.5}
                ]
            },
        }
        params = grammar_to_predicate_params(grammar)
        assert params["cryptic_splice_threshold"] == 5.0
        assert params["uncertain_lo"] == 2.5

    def test_no_restriction_site_params(self):
        """NoRestrictionSite predicate params are extracted from grammar."""
        grammar = {
            "gene": {"name": "T"},
            "type_system": {
                "predicates": [{"name": "NoRestrictionSite", "enzyme_sites": ["EcoRI", "XhoI"]}]
            },
        }
        params = grammar_to_predicate_params(grammar)
        assert params["enzymes"] == ["EcoRI", "XhoI"]

    def test_cellular_context_from_ndfst(self):
        """cellular_context comes from ndfst.cell_type."""
        grammar = {
            "gene": {"name": "T"},
            "ndfst": {"cell_type": "HeLa"},
            "type_system": {"predicates": []},
        }
        params = grammar_to_predicate_params(grammar)
        assert params["cellular_context"] == "HeLa"

    def test_cellular_context_default(self):
        """Default cellular_context is used when ndfst has no cell_type."""
        grammar = {"gene": {"name": "T"}, "type_system": {"predicates": []}}
        params = grammar_to_predicate_params(grammar)
        assert params["cellular_context"] == _DEFAULT_CELLULAR_CONTEXT

    def test_exon_boundaries_from_exons(self):
        """Exon boundaries are built from exons with 1-based to 0-based conversion."""
        grammar = {
            "gene": {"name": "T"},
            "exons": [
                {"id": "e1", "start": 1, "end": 142},
                {"id": "e2", "start": 273, "end": 495},
            ],
            "type_system": {"predicates": []},
        }
        params = grammar_to_predicate_params(grammar)
        # start is converted: 1 → 0, 273 → 272; end stays as-is
        assert params["exon_boundaries"] == [(0, 142), (272, 495)]

    def test_exon_boundaries_default_when_no_exons(self):
        """Default exon_boundaries when grammar has no exons section."""
        grammar = {"gene": {"name": "T"}, "type_system": {"predicates": []}}
        params = grammar_to_predicate_params(grammar)
        assert params["exon_boundaries"] == list(_DEFAULT_EXON_BOUNDARIES)

    def test_exons_without_start_end_ignored(self):
        """Exon entries without start/end are silently skipped."""
        grammar = {
            "gene": {"name": "T"},
            "exons": [
                {"id": "e1", "start": 1, "end": 100},
                {"id": "e2", "description": "incomplete"},
            ],
            "type_system": {"predicates": []},
        }
        params = grammar_to_predicate_params(grammar)
        assert params["exon_boundaries"] == [(0, 100)]

    def test_default_enzymes_is_copy(self):
        """Returned default enzymes list is a copy, not the module constant."""
        grammar = {"gene": {"name": "T"}, "type_system": {"predicates": []}}
        params1 = grammar_to_predicate_params(grammar)
        params2 = grammar_to_predicate_params(grammar)
        params1["enzymes"].append("NewEnzyme")
        # params2 should not be affected
        assert "NewEnzyme" not in params2["enzymes"]


# ────────────────────────────────────────────────────────────
# 5. Integration: full grammar round-trip
# ────────────────────────────────────────────────────────────

class TestGrammarRoundTrip:
    """Integration tests: load a full grammar, then extract params."""

    def test_egfp_builtin_round_trip(self):
        """Load built-in EGFP grammar and extract params without error."""
        grammar = load_builtin_grammar("egfp_hek293t")
        params = grammar_to_predicate_params(grammar)
        assert params["organism"] == "Homo_sapiens"
        assert params["cellular_context"] == "HEK293T"
        # EGFP has a single exon 1..720 → boundary (0, 720)
        assert (0, 720) in params["exon_boundaries"]
        # Grammar explicitly sets these values
        assert params["gc_lo"] == 0.30
        assert params["gc_hi"] == 0.70
        assert params["cai_threshold"] == 0.70
        assert params["cryptic_splice_threshold"] == 3.0

    def test_hbb_builtin_round_trip(self):
        """Load built-in HBB grammar and extract params without error."""
        grammar = load_builtin_grammar("hbb_hek293t")
        params = grammar_to_predicate_params(grammar)
        assert params["organism"] == "Homo_sapiens"
        assert params["cellular_context"] == "HEK293T"
        # HBB has three exons
        assert len(params["exon_boundaries"]) == 3
        # First exon: start=1 → 0-based start=0, end=142
        assert params["exon_boundaries"][0] == (0, 142)
        assert params["gc_lo"] == 0.40
        assert params["gc_hi"] == 0.60

    def test_custom_grammar_round_trip(self, tmp_path):
        """Load a custom grammar from file and extract params."""
        content = """\
        gene:
          name: "CustomGene"
          expression_organism: "Saccharomyces_cerevisiae"
        ndfst:
          cell_type: "Yeast"
        exons:
          - id: "exon1"
            start: 1
            end: 500
        type_system:
          predicates:
            - name: "SpliceCorrect"
            - name: "GCInRange"
              lo: 0.35
              hi: 0.55
            - name: "CodonAdapted"
              threshold: 0.60
            - name: "NoRestrictionSite"
              enzyme_sites:
                - "EcoRI"
        """
        path = _write_grammar(tmp_path, content)
        grammar = load_grammar(path)
        params = grammar_to_predicate_params(grammar)
        assert params["organism"] == "Saccharomyces_cerevisiae"
        assert params["cellular_context"] == "Yeast"
        assert params["exon_boundaries"] == [(0, 500)]
        assert params["gc_lo"] == 0.35
        assert params["gc_hi"] == 0.55
        assert params["cai_threshold"] == 0.60
        assert params["enzymes"] == ["EcoRI"]


# ────────────────────────────────────────────────────────────
# 6. YAML import guard
# ────────────────────────────────────────────────────────────

class TestYAMLImportGuard:
    """Test that missing PyYAML is properly detected."""

    @pytest.mark.skip(reason="Shim module doesn't support monkeypatching internal 'yaml' attribute")
    def test_check_yaml_available_raises_when_missing(self, monkeypatch):
        """_check_yaml_available raises ImportError if yaml is None."""
        from biocompiler import grammar_loader as gl
        monkeypatch.setattr(gl, "yaml", None)
        with pytest.raises(ImportError, match="PyYAML is required"):
            gl._check_yaml_available()

    @pytest.mark.skip(reason="Shim module doesn't support monkeypatching internal 'yaml' attribute")
    def test_load_grammar_raises_when_yaml_missing(self, monkeypatch, tmp_path):
        """load_grammar raises ImportError when PyYAML is not available."""
        from biocompiler import grammar_loader as gl
        monkeypatch.setattr(gl, "yaml", None)
        path = tmp_path / "dummy.yaml"
        path.write_text("key: value\n")
        with pytest.raises(ImportError, match="PyYAML is required"):
            gl.load_grammar(str(path))
