"""
BioCompiler Structure Subpackage
=================================
Consolidated structure analysis, quality assessment, type-system
predicates, and assessment reporting.

Previously these lived in four separate flat modules:
    biocompiler.structure
    biocompiler.structure_quality
    biocompiler.structure_predicates
    biocompiler.structure_report

All public names are re-exported here so that external code can do::

    from biocompiler.structure import ProteinStructure, compute_structure_quality
"""

# ── Submodule references ────────────────────────────────────────
from . import parser, quality, predicates, report

# ── Parser (PDB parsing, data models, dihedrals) ────────────────
from .parser import (
    Atom,
    Residue,
    Chain,
    ProteinStructure,
    parse_pdb,
    parse_pdb_file,
    compute_dihedral,
    compute_ramachandran,
    secondary_structure_estimate,
    THREE_TO_ONE,
    ONE_TO_THREE,
)

# ── Quality (structure quality assessment) ──────────────────────
from .quality import (
    StructureQualityReport,
    assess_plddt,
    assess_ramachandran,
    compute_clash_score,
    compute_packing_density,
    compute_exposed_hydrophobic,
    compute_structure_quality,
    find_low_confidence_regions,
    compute_sasa_approximation,
    KYTE_DOOLITTLE,
    VDW_RADII,
    MAX_SASA,
)

# ── Predicates (type-system structure predicates) ──────────────
from .predicates import (
    evaluate_structure_confidence,
    evaluate_no_misfolding_risk,
    evaluate_correct_fold_topology,
    evaluate_no_unexpected_interaction,
    expected_radius_of_gyration,
    compute_secondary_structure_fractions,
    find_surface_charge_patches,
    find_unstructured_regions,
)

# ── Report (assessment reporting & visualization) ──────────────
# NOTE: These symbols live in biocompiler.structure.report (NOT
# biocompiler.provenance.report). A previous version imported them
# from the wrong module, raising ImportError and silently degrading
# the entire structure subsystem to None fallbacks in the top-level
# biocompiler.__init__ (see historical audit notes).
from .report import (
    ProteinAssessmentReport,
    assess_protein,
    format_assessment_text,
    format_assessment_json,
    format_assessment_html,
    generate_recommendations,
    compute_overall_verdict,
    plot_plddt_bar_svg,
    plot_solubility_profile_svg,
)

__all__ = [
    # parser
    "Atom",
    "Residue",
    "Chain",
    "ProteinStructure",
    "parse_pdb",
    "parse_pdb_file",
    "compute_dihedral",
    "compute_ramachandran",
    "secondary_structure_estimate",
    "THREE_TO_ONE",
    "ONE_TO_THREE",
    # quality
    "StructureQualityReport",
    "assess_plddt",
    "assess_ramachandran",
    "compute_clash_score",
    "compute_packing_density",
    "compute_exposed_hydrophobic",
    "compute_structure_quality",
    "find_low_confidence_regions",
    "compute_sasa_approximation",
    "KYTE_DOOLITTLE",
    "VDW_RADII",
    "MAX_SASA",
    # predicates
    "evaluate_structure_confidence",
    "evaluate_no_misfolding_risk",
    "evaluate_correct_fold_topology",
    "evaluate_no_unexpected_interaction",
    "expected_radius_of_gyration",
    "compute_secondary_structure_fractions",
    "find_surface_charge_patches",
    "find_unstructured_regions",
    # report
    "ProteinAssessmentReport",
    "assess_protein",
    "format_assessment_text",
    "format_assessment_json",
    "format_assessment_html",
    "generate_recommendations",
    "compute_overall_verdict",
    "plot_plddt_bar_svg",
    "plot_solubility_profile_svg",
]
