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

Also includes ligand_binding_v2 for protein-ligand binding analysis with
optional RDKit, AutoDock Vina, and ProLIF integration, prolif_integration
for comprehensive protein-ligand interaction fingerprint analysis,
vina_wrapper for dedicated AutoDock Vina molecular docking, and
rdkit_integration for RDKit-based cheminformatics (SMILES parsing,
molecular descriptors, conformer generation, pharmacophore extraction,
RMSD, conformer clustering).

All public names are re-exported here so that external code can do::

    from biocompiler.structure import ProteinStructure, compute_structure_quality
    from biocompiler.structure import detect_binding_sites, parse_smiles_features_rdkit
"""

# ── Submodule references ────────────────────────────────────────
from . import parser, quality, predicates, report, ligand_binding_v2, prolif_integration, vina_wrapper, rdkit_integration

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

# ── Ligand Binding v2 (protein-ligand analysis) ──────────────
from .ligand_binding_v2 import (
    BindingSite,
    PharmacophoreFeature,
    LigandInfo,
    DockingResult,
    detect_binding_sites,
    score_binding_site,
    decompose_per_residue_energy,
    parse_smiles_features_rdkit,
    generate_3d_conformer,
    dock_ligand_vina,
    compute_interaction_fingerprint,
    HAS_RDKIT,
    HAS_VINA,
    HAS_PROLIF,
)

# ── ProLIF Integration (interaction fingerprints) ──────────────
from .prolif_integration import (
    InteractionFingerprint,
    InteractionReport,
    compute_interaction_fingerprint as prolif_compute_interaction_fingerprint,
    compute_interaction_fingerprint_from_mols,
    detect_interaction_types,
    estimate_binding_affinity,
    compare_interaction_patterns,
    is_prolif_available,
    SUPPORTED_INTERACTION_TYPES,
)

# ── Vina Wrapper (AutoDock Vina molecular docking) ─────────────
from .vina_wrapper import (
    VinaResult,
    DockingConfig,
    dock_smiles,
    dock_pdbqt,
    smiles_to_pdbqt,
    pdb_to_pdbqt,
    compute_box_from_receptor,
    is_vina_available,
    score_ligand_binding,
)

# ── RDKit Integration (cheminformatics for ligand analysis) ─────
from .rdkit_integration import (
    ConformerInfo,
    PharmacophoreFeature3D,
    is_rdkit_available,
    parse_smiles,
    compute_molecular_descriptors,
    generate_conformers,
    compute_partial_charges,
    compute_rotatable_bonds,
    extract_pharmacophore_features_3d,
    compute_rmsd,
    cluster_poses_by_rmsd,
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
    # ligand_binding_v2
    "BindingSite",
    "PharmacophoreFeature",
    "LigandInfo",
    "DockingResult",
    "detect_binding_sites",
    "score_binding_site",
    "decompose_per_residue_energy",
    "parse_smiles_features_rdkit",
    "generate_3d_conformer",
    "dock_ligand_vina",
    "compute_interaction_fingerprint",
    "HAS_RDKIT",
    "HAS_VINA",
    "HAS_PROLIF",
    # prolif_integration
    "InteractionFingerprint",
    "InteractionReport",
    "prolif_compute_interaction_fingerprint",
    "compute_interaction_fingerprint_from_mols",
    "detect_interaction_types",
    "estimate_binding_affinity",
    "compare_interaction_patterns",
    "is_prolif_available",
    "SUPPORTED_INTERACTION_TYPES",
    # vina_wrapper
    "VinaResult",
    "DockingConfig",
    "dock_smiles",
    "dock_pdbqt",
    "smiles_to_pdbqt",
    "pdb_to_pdbqt",
    "compute_box_from_receptor",
    "is_vina_available",
    "score_ligand_binding",
    # rdkit_integration
    "ConformerInfo",
    "PharmacophoreFeature3D",
    "is_rdkit_available",
    "parse_smiles",
    "compute_molecular_descriptors",
    "generate_conformers",
    "compute_partial_charges",
    "compute_rotatable_bonds",
    "extract_pharmacophore_features_3d",
    "compute_rmsd",
    "cluster_poses_by_rmsd",
]
