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
)

# ── Predicates (type-system structure predicates) ──────────────
from .predicates import (
    evaluate_structure_confidence,
    evaluate_no_misfolding_risk,
    evaluate_correct_fold_topology,
    evaluate_no_unexpected_interaction,
)

# ── Report (assessment reporting & visualization) ──────────────
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
    # predicates
    "evaluate_structure_confidence",
    "evaluate_no_misfolding_risk",
    "evaluate_correct_fold_topology",
    "evaluate_no_unexpected_interaction",
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
