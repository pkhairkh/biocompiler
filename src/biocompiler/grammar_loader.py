"""
BioCompiler Grammar Loader — YAML Configuration Support

Loads YAML grammar files that define gene-specific NDFST rules,
type system parameters, and organism-specific settings.
"""

import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore


def _check_yaml_available():
    """Check that PyYAML is available."""
    if yaml is None:
        raise ImportError(
            "PyYAML is required for YAML grammar loading. "
            "Install it with: pip install pyyaml"
        )


def load_grammar(path: str | Path) -> dict[str, Any]:
    """
    Load a BioCompiler grammar from a YAML file.

    The grammar defines:
    - gene: name, organism, reference
    - pre_mrna: sequence source and length
    - exons/introns: structure
    - ndfst: cell type, splice consensus, expected isoforms
    - type_system: predicate parameters (thresholds, organisms, etc.)

    Args:
        path: path to the YAML grammar file

    Returns:
        Parsed grammar as a dict

    Raises:
        ImportError: if PyYAML is not installed
        FileNotFoundError: if the grammar file doesn't exist
        ValueError: if the grammar file is malformed
    """
    _check_yaml_available()

    grammar_path = Path(path)
    if not grammar_path.exists():
        raise FileNotFoundError(f"Grammar file not found: {path}")

    with open(grammar_path) as f:
        grammar = yaml.safe_load(f)

    if not isinstance(grammar, dict):
        raise ValueError(f"Grammar file must contain a YAML mapping, got {type(grammar)}")

    # Validate required top-level keys
    _validate_grammar(grammar, path)
    return grammar


def _validate_grammar(grammar: dict, path: str):
    """Validate that a grammar has required structure."""
    issues = []

    if "gene" not in grammar:
        issues.append("Missing 'gene' section")
    elif not isinstance(grammar["gene"], dict):
        issues.append("'gene' must be a mapping")

    if "type_system" not in grammar:
        issues.append("Missing 'type_system' section")
    elif not isinstance(grammar["type_system"], dict):
        issues.append("'type_system' must be a mapping")

    # Validate predicates
    type_sys = grammar.get("type_system", {})
    predicates = type_sys.get("predicates", [])
    if not isinstance(predicates, list):
        issues.append("'type_system.predicates' must be a list")
    else:
        known_predicates = {
            "SpliceCorrect", "NoCrypticSplice", "CodonAdapted",
            "GCInRange", "NoRestrictionSite", "InFrame",
            "NoInstabilityMotif", "NoCpGIsland",
        }
        for i, pred in enumerate(predicates):
            if not isinstance(pred, dict):
                issues.append(f"Predicate {i} must be a mapping")
            elif "name" not in pred:
                issues.append(f"Predicate {i} missing 'name' key")
            elif pred["name"] not in known_predicates:
                logger.warning(
                    "Predicate '%s' (index %d) is not a built-in predicate. "
                    "Make sure it is registered in the predicate registry.",
                    pred["name"], i,
                )

    if issues:
        raise ValueError(f"Grammar validation errors in {path}: {'; '.join(issues)}")


def grammar_to_predicate_params(grammar: dict) -> dict[str, Any]:
    """
    Extract type-checking parameters from a grammar definition.

    Converts the YAML grammar's type_system.predicates section into
    a format suitable for evaluate_all_predicates().

    Args:
        grammar: parsed grammar dict from load_grammar()

    Returns:
        Dict of parameters for type checking:
        - organism, gc_lo, gc_hi, cai_threshold
        - exon_boundaries, cellular_context, enzymes
        - cryptic_splice_threshold
    """
    gene = grammar.get("gene", {})
    type_sys = grammar.get("type_system", {})
    ndfst = grammar.get("ndfst", {})
    exons = grammar.get("exons", [])
    predicates = type_sys.get("predicates", [])

    params: dict[str, Any] = {
        "organism": gene.get("expression_organism", gene.get("organism", "Homo_sapiens")),
        "cellular_context": ndfst.get("cell_type", "HEK293T"),
    }

    # Build exon boundaries from grammar
    if exons:
        params["exon_boundaries"] = [
            (e["start"] - 1, e["end"])  # Convert 1-based to 0-based
            for e in exons
            if "start" in e and "end" in e
        ]

    # Extract predicate-specific parameters
    for pred in predicates:
        name = pred.get("name", "")
        if name == "GCInRange":
            params.setdefault("gc_lo", pred.get("lo", 0.30))
            params.setdefault("gc_hi", pred.get("hi", 0.70))
        elif name == "CodonAdapted":
            params.setdefault("cai_threshold", pred.get("threshold", 0.5))
        elif name == "NoCrypticSplice":
            params.setdefault("cryptic_splice_threshold", pred.get("cryptic_threshold", 3.0))
        elif name == "NoRestrictionSite":
            params.setdefault("enzymes", pred.get("enzyme_sites", []))

    # Set defaults for any missing params
    params.setdefault("gc_lo", 0.30)
    params.setdefault("gc_hi", 0.70)
    params.setdefault("cai_threshold", 0.5)
    params.setdefault("cryptic_splice_threshold", 3.0)
    params.setdefault("enzymes", ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"])
    params.setdefault("exon_boundaries", [(0, 0)])

    return params
