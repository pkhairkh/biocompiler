"""
BioCompiler Grammar Loader — YAML Configuration Support

Loads YAML grammar files that define gene-specific NDFST rules,
type system parameters, and organism-specific settings.

Built-in grammars are shipped inside the ``biocompiler.grammars`` package
directory and can be loaded by name via :func:`load_builtin_grammar`.
User-supplied grammar files can be loaded by path via :func:`load_grammar`.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from biocompiler.sequence.maxentscan import (
    CRYPTIC_SPLICE_THRESHOLD as _CRYPTIC_SPLICE_THRESHOLD_FROM_MAXENT,
)

__all__ = [
    "load_grammar",
    "grammar_to_predicate_params",
    "list_builtin_grammars",
    "load_builtin_grammar",
]

logger = logging.getLogger(__name__)

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore[assignment]
    logger.debug("PyYAML not installed; YAML grammar loading will be unavailable.")

# Directory where built-in grammar YAML files live
_GRAMMARS_DIR: Path = Path(__file__).resolve().parent.parent / "grammars"

# ── 1-based to 0-based coordinate offset ─────────────────────────────────
_ONE_BASED_OFFSET: int = 1

# ── Default parameter constants (used when grammar omits values) ──────────
_DEFAULT_GC_LO: float = 0.30
_DEFAULT_GC_HI: float = 0.70
_DEFAULT_CAI_THRESHOLD: float = 0.5
_DEFAULT_CRYPTIC_SPLICE_THRESHOLD: float = _CRYPTIC_SPLICE_THRESHOLD_FROM_MAXENT
# Derived from maxentscan.CRYPTIC_SPLICE_THRESHOLD (single source of truth).
# Previously hardcoded to 3.0 (PWM-era value); now automatically tracks
# the Markov-model-calibrated threshold.
_DEFAULT_UNCERTAIN_LO: float = 1.5
_DEFAULT_ENZYMES: list[str] = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
_DEFAULT_EXON_BOUNDARIES: list[tuple[int, int]] = [(0, 0)]
_DEFAULT_ORGANISM: str = "Homo_sapiens"
_DEFAULT_CELLULAR_CONTEXT: str = "HEK293T"


def _check_yaml_available() -> None:
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
        FileNotFoundError: if the grammar file does not exist
        ValueError: if the grammar file is malformed
    """
    _check_yaml_available()

    grammar_path = Path(path)
    if not grammar_path.exists():
        raise FileNotFoundError(f"Grammar file not found: {path}")

    with open(grammar_path, encoding="utf-8") as f:
        grammar = yaml.safe_load(f)

    if not isinstance(grammar, dict):
        raise ValueError(f"Grammar file must contain a YAML mapping, got {type(grammar)}")

    # Validate required top-level keys
    _validate_grammar(grammar, path)
    return grammar


def _validate_grammar(grammar: dict[str, Any], path: str | Path) -> None:
    """Validate that a grammar has required structure."""
    issues: list[str] = []

    if "gene" not in grammar:
        issues.append("Missing 'gene' section")
    elif not isinstance(grammar["gene"], dict):
        issues.append("'gene' must be a mapping")

    type_sys_is_valid: bool = True
    if "type_system" not in grammar:
        issues.append("Missing 'type_system' section")
        type_sys_is_valid = False
    elif not isinstance(grammar["type_system"], dict):
        issues.append("'type_system' must be a mapping")
        type_sys_is_valid = False

    # Validate predicates only when type_system is a proper dict
    if type_sys_is_valid:
        type_sys: dict[str, Any] = grammar["type_system"]
        predicates: Any = type_sys.get("predicates", [])
        if not isinstance(predicates, list):
            issues.append("'type_system.predicates' must be a list")
        else:
            known_predicates: set[str] = {
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


def grammar_to_predicate_params(grammar: dict[str, Any]) -> dict[str, Any]:
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
    gene: Any = grammar.get("gene", {})
    type_sys: Any = grammar.get("type_system", {})
    ndfst: Any = grammar.get("ndfst", {})
    exons: Any = grammar.get("exons", [])

    # Guard against non-dict sections (defensive — validation should catch this)
    if not isinstance(gene, dict):
        logger.warning("gene is not a mapping in grammar; treating as empty.")
        gene = {}
    if not isinstance(type_sys, dict):
        logger.warning("type_system is not a mapping in grammar; treating as empty.")
        type_sys = {}
    if not isinstance(ndfst, dict):
        logger.warning("ndfst is not a mapping in grammar; treating as empty.")
        ndfst = {}
    predicates: Any = type_sys.get("predicates", [])

    params: dict[str, Any] = {
        "organism": gene.get("expression_organism", gene.get("organism", _DEFAULT_ORGANISM)),
        "cellular_context": ndfst.get("cell_type", _DEFAULT_CELLULAR_CONTEXT),
    }

    # Build exon boundaries from grammar (convert 1-based → 0-based)
    if exons:
        params["exon_boundaries"] = [
            (e["start"] - _ONE_BASED_OFFSET, e["end"])
            for e in exons
            if "start" in e and "end" in e
        ]

    # Extract predicate-specific parameters
    for pred in predicates:
        name = pred.get("name", "")
        if name == "GCInRange":
            params.setdefault("gc_lo", pred.get("lo", _DEFAULT_GC_LO))
            params.setdefault("gc_hi", pred.get("hi", _DEFAULT_GC_HI))
        elif name == "CodonAdapted":
            params.setdefault("cai_threshold", pred.get("threshold", _DEFAULT_CAI_THRESHOLD))
        elif name == "NoCrypticSplice":
            params.setdefault("cryptic_splice_threshold", pred.get("cryptic_threshold", _DEFAULT_CRYPTIC_SPLICE_THRESHOLD))
            params.setdefault("uncertain_lo", pred.get("uncertain_lo", _DEFAULT_UNCERTAIN_LO))
        elif name == "NoRestrictionSite":
            params.setdefault("enzymes", pred.get("enzyme_sites", []))

    # Set defaults for any missing params
    params.setdefault("gc_lo", _DEFAULT_GC_LO)
    params.setdefault("gc_hi", _DEFAULT_GC_HI)
    params.setdefault("cai_threshold", _DEFAULT_CAI_THRESHOLD)
    params.setdefault("cryptic_splice_threshold", _DEFAULT_CRYPTIC_SPLICE_THRESHOLD)
    params.setdefault("uncertain_lo", _DEFAULT_UNCERTAIN_LO)
    params.setdefault("enzymes", list(_DEFAULT_ENZYMES))
    params.setdefault("exon_boundaries", list(_DEFAULT_EXON_BOUNDARIES))

    return params


def list_builtin_grammars() -> list[str]:
    """
    List the names of built-in grammar files shipped with BioCompiler.

    Built-in grammars are YAML files in the ``biocompiler.grammars`` package
    directory. They can be loaded by name using :func:`load_builtin_grammar`.

    Returns:
        List of grammar names (without ``.yaml`` extension)
    """
    if not _GRAMMARS_DIR.exists():
        return []
    return sorted(p.stem for p in _GRAMMARS_DIR.glob("*.yaml"))


def load_builtin_grammar(name: str) -> dict[str, Any]:
    """
    Load a built-in grammar by name.

    Built-in grammars are YAML files shipped inside the ``biocompiler.grammars``
    package directory. Use :func:`list_builtin_grammars` to discover available
    names.

    Args:
        name: Grammar name (without ``.yaml`` extension), e.g. ``"egfp_hek293t"``

    Returns:
        Parsed grammar as a dict

    Raises:
        ImportError: if PyYAML is not installed
        FileNotFoundError: if no built-in grammar with that name exists
        ValueError: if the grammar file is malformed
    """
    path = _GRAMMARS_DIR / f"{name}.yaml"
    if not path.exists():
        available = list_builtin_grammars()
        raise FileNotFoundError(
            f"No built-in grammar named '{name}'. "
            f"Available grammars: {available}"
        )
    return load_grammar(path)
