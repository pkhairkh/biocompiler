"""Precomputed MHC binding data for HLA-DRB1*01:01.

MHC-II presents 13-25 aa peptides with a 9-mer core binding region.
Peptides are 15-mers scored by scanning all possible
9-mer core registers within each peptide.

Anchor residues: position 1 (F/Y/W), position 4 (L/I/V), position 6 (S/T), position 9 (V/L)
Known epitopes: PKYVKQNTLKLAT, YVKQNTLKLAT
Data source: IEDB and SYFPEITHI databases + PSSM predictions

Binding scores are derived from PSSM-based prediction with expected
AUC-ROC of 0.60-0.75. Do not replace experimental validation or
NetMHCpan predictions where available.

References
----------
- Southwood et al., J Immunol 1998; 160:3363 (DR binding motifs)
- Rammensee et al., Immunogenetics 1995 (SYFPEITHI)
"""
from __future__ import annotations

from .. import PrecomputedAlleleDatabase, PrecomputedEntry

ALLELE: str = "HLA-DRB1*01:01"
MHC_CLASS: str = "II"
PEPTIDE_LENGTH: int = 15
ANCHOR_DESCRIPTION: str = "position 1 (F/Y/W), position 4 (L/I/V), position 6 (S/T), position 9 (V/L)"
KNOWN_EPITOPES: list[str] = ['PKYVKQNTLKLAT', 'YVKQNTLKLAT']

# ═══════════════════════════════════════════════════════════════════════════
# Position-Specific Scoring Matrix (PSSM) — core 9-mer
# ═══════════════════════════════════════════════════════════════════════════

PSSM: list[dict[str, float]] = [
    {"A": 0.9, "C": 0.9, "D": 0.4, "E": 0.4, "F": 1.8, "G": 0.9, "H": 0.9, "I": 1.3, "K": 0.5, "L": 1.4, "M": 1.2, "N": 0.9, "P": 0.9, "Q": 0.9, "R": 0.5, "S": 0.9, "T": 0.9, "V": 1.3, "W": 1.6, "Y": 1.7},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 0.8, "C": 0.8, "D": 0.4, "E": 0.4, "F": 1.2, "G": 0.8, "H": 0.8, "I": 1.6, "K": 0.5, "L": 1.7, "M": 1.3, "N": 0.8, "P": 0.5, "Q": 0.8, "R": 0.5, "S": 0.8, "T": 0.8, "V": 1.5, "W": 0.8, "Y": 0.8},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.2, "C": 0.8, "D": 0.8, "E": 0.8, "F": 0.5, "G": 0.8, "H": 0.8, "I": 0.8, "K": 0.8, "L": 0.8, "M": 0.8, "N": 1.1, "P": 0.8, "Q": 0.8, "R": 0.6, "S": 1.7, "T": 1.6, "V": 0.8, "W": 0.5, "Y": 0.6},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.1, "C": 0.8, "D": 0.4, "E": 0.4, "F": 0.8, "G": 0.8, "H": 0.8, "I": 1.4, "K": 0.5, "L": 1.6, "M": 1.2, "N": 0.8, "P": 0.4, "Q": 0.8, "R": 0.5, "S": 0.8, "T": 0.8, "V": 1.7, "W": 0.8, "Y": 0.8},
]


def _build_entries() -> list[PrecomputedEntry]:
    """Construct all precomputed entries for HLA-DRB1*01:01."""
    entries: list[PrecomputedEntry] = []

    # --- Binders (strong + moderate + weak) ---
    entries.append(PrecomputedEntry(
        peptide="LIVLDSGLVFFYRTH",
        binding_score=0.884956,
        ic50_nm=54.53,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 3: "L", 5: "S", 8: "V"},
        anchor_scores={0: 1.4, 3: 1.7, 5: 1.7, 8: 1.7},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KNGYDVINFLFSLEV",
        binding_score=0.852953,
        ic50_nm=65.56,
        binding_class="moderate_binder",
        anchor_residues={6: "I", 9: "L", 11: "S", 14: "V"},
        anchor_scores={6: 1.3, 9: 1.7, 11: 1.7, 14: 1.7},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NKKAQNLPWLPSIEI",
        binding_score=0.802896,
        ic50_nm=87.45,
        binding_class="moderate_binder",
        anchor_residues={6: "L", 9: "L", 11: "S", 14: "I"},
        anchor_scores={6: 1.4, 9: 1.7, 11: 1.7, 14: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="QPFAAWYYIFASVSL",
        binding_score=0.794516,
        ic50_nm=91.77,
        binding_class="moderate_binder",
        anchor_residues={6: "Y", 9: "F", 11: "S", 14: "L"},
        anchor_scores={6: 1.7, 9: 1.2, 11: 1.7, 14: 1.6},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="QHPQFIHICAKSLFV",
        binding_score=0.793113,
        ic50_nm=92.52,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "I", 9: "A", 12: "L"},
        anchor_scores={4: 1.8, 7: 1.6, 9: 1.2, 12: 1.6},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="FVPIGTCGMEVFEEY",
        binding_score=0.793113,
        ic50_nm=92.52,
        binding_class="moderate_binder",
        anchor_residues={0: "F", 3: "I", 5: "T", 8: "M"},
        anchor_scores={0: 1.8, 3: 1.6, 5: 1.6, 8: 1.2},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="LYKRLHNVLVKPRDI",
        binding_score=0.783884,
        ic50_nm=97.56,
        binding_class="moderate_binder",
        anchor_residues={1: "Y", 4: "L", 6: "N", 9: "V"},
        anchor_scores={1: 1.7, 4: 1.7, 6: 1.1, 9: 1.7},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HMEQRCYCEMLTDYI",
        binding_score=0.749269,
        ic50_nm=119.08,
        binding_class="moderate_binder",
        anchor_residues={6: "Y", 9: "M", 11: "T", 14: "I"},
        anchor_scores={6: 1.7, 9: 1.3, 11: 1.6, 14: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WQNINNDQLINVMSY",
        binding_score=0.71335,
        ic50_nm=146.43,
        binding_class="moderate_binder",
        anchor_residues={0: "W", 3: "I", 5: "N", 8: "L"},
        anchor_scores={0: 1.6, 3: 1.6, 5: 1.1, 8: 1.6},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HFVMLMASQMDNGAR",
        binding_score=0.705039,
        ic50_nm=153.6,
        binding_class="moderate_binder",
        anchor_residues={1: "F", 4: "L", 6: "A", 9: "M"},
        anchor_scores={1: 1.8, 4: 1.7, 6: 1.2, 9: 1.2},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="INHHPPLTTLETNKA",
        binding_score=0.68639,
        ic50_nm=171.01,
        binding_class="moderate_binder",
        anchor_residues={6: "L", 9: "L", 11: "T", 14: "A"},
        anchor_scores={6: 1.4, 9: 1.7, 11: 1.6, 14: 1.1},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="YMFTAWTICSFEVML",
        binding_score=0.684017,
        ic50_nm=173.36,
        binding_class="moderate_binder",
        anchor_residues={4: "A", 7: "I", 9: "S", 12: "V"},
        anchor_scores={4: 0.9, 7: 1.6, 9: 1.7, 12: 1.7},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="FRYWWNIAPVLTGIM",
        binding_score=0.646301,
        ic50_nm=215.4,
        binding_class="moderate_binder",
        anchor_residues={6: "I", 9: "V", 11: "T", 14: "M"},
        anchor_scores={6: 1.3, 9: 1.5, 11: 1.6, 14: 1.2},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KFNFNLRSTYIYPEM",
        binding_score=0.63664,
        ic50_nm=227.72,
        binding_class="moderate_binder",
        anchor_residues={2: "N", 5: "L", 7: "S", 10: "I"},
        anchor_scores={2: 0.9, 5: 1.7, 7: 1.7, 10: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WARNFYKIICAKSAI",
        binding_score=0.631775,
        ic50_nm=234.19,
        binding_class="moderate_binder",
        anchor_residues={5: "Y", 8: "I", 10: "A", 13: "A"},
        anchor_scores={5: 1.7, 8: 1.6, 10: 1.2, 13: 1.1},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PKYVKQNTLKLAT",
        binding_score=0.621252,
        ic50_nm=248.81,
        binding_class="moderate_binder",
        anchor_residues={2: "Y", 5: "Q", 7: "T", 10: "L"},
        anchor_scores={2: 1.7, 5: 0.8, 7: 1.6, 10: 1.6},
        source="known_epitope",
        peptide_length=13,
    ))
    entries.append(PrecomputedEntry(
        peptide="YVKQNTLKLAT",
        binding_score=0.621252,
        ic50_nm=248.81,
        binding_class="moderate_binder",
        anchor_residues={0: "Y", 3: "Q", 5: "T", 8: "L"},
        anchor_scores={0: 1.7, 3: 0.8, 5: 1.6, 8: 1.6},
        source="known_epitope",
        peptide_length=11,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKGWRTFDNMRADSM",
        binding_score=0.610206,
        ic50_nm=265.15,
        binding_class="moderate_binder",
        anchor_residues={6: "F", 9: "M", 11: "A", 14: "M"},
        anchor_scores={6: 1.8, 9: 1.3, 11: 1.2, 14: 1.2},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRVSVFQSFGNADID",
        binding_score=0.605885,
        ic50_nm=271.82,
        binding_class="moderate_binder",
        anchor_residues={5: "F", 8: "F", 10: "N", 13: "I"},
        anchor_scores={5: 1.8, 8: 1.2, 10: 1.1, 13: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="LYKLCQIIVVRDYWK",
        binding_score=0.596813,
        ic50_nm=286.4,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 3: "L", 5: "Q", 8: "V"},
        anchor_scores={0: 1.4, 3: 1.7, 5: 0.8, 8: 1.7},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="YYMNMNAPRMMRRKS",
        binding_score=0.591229,
        ic50_nm=295.75,
        binding_class="moderate_binder",
        anchor_residues={1: "Y", 4: "M", 6: "A", 9: "M"},
        anchor_scores={1: 1.7, 4: 1.3, 6: 1.2, 9: 1.2},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="ACIMDLCTVEQLCEP",
        binding_score=0.553433,
        ic50_nm=367.64,
        binding_class="moderate_binder",
        anchor_residues={2: "I", 5: "L", 7: "T", 10: "Q"},
        anchor_scores={2: 1.3, 5: 1.7, 7: 1.6, 10: 0.8},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="VKLPKVEEFSARWIP",
        binding_score=0.529845,
        ic50_nm=421.1,
        binding_class="moderate_binder",
        anchor_residues={5: "V", 8: "F", 10: "A", 13: "I"},
        anchor_scores={5: 1.3, 8: 1.2, 10: 1.2, 13: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="ALDWIKFYRVTKMCM",
        binding_score=0.526491,
        ic50_nm=429.31,
        binding_class="moderate_binder",
        anchor_residues={6: "F", 9: "V", 11: "K", 14: "M"},
        anchor_scores={6: 1.8, 9: 1.5, 11: 0.8, 14: 1.2},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PYPQLWFWWVECRQM",
        binding_score=0.526491,
        ic50_nm=429.31,
        binding_class="moderate_binder",
        anchor_residues={6: "F", 9: "V", 11: "C", 14: "M"},
        anchor_scores={6: 1.8, 9: 1.5, 11: 0.8, 14: 1.2},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="IKYVTAEAAAEFNWT",
        binding_score=0.524383,
        ic50_nm=434.55,
        binding_class="moderate_binder",
        anchor_residues={0: "I", 3: "V", 5: "A", 8: "A"},
        anchor_scores={0: 1.3, 3: 1.5, 5: 1.2, 8: 1.1},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="IWLPTVKNTNAGWTH",
        binding_score=0.520496,
        ic50_nm=444.38,
        binding_class="moderate_binder",
        anchor_residues={2: "L", 5: "V", 7: "N", 10: "A"},
        anchor_scores={2: 1.4, 5: 1.5, 7: 1.1, 10: 1.1},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="ILVVMTANWADEQCD",
        binding_score=0.515147,
        ic50_nm=458.28,
        binding_class="moderate_binder",
        anchor_residues={0: "I", 3: "V", 5: "T", 8: "W"},
        anchor_scores={0: 1.3, 3: 1.5, 5: 1.6, 8: 0.8},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="TMIANFQNNIIERQA",
        binding_score=0.503831,
        ic50_nm=489.13,
        binding_class="moderate_binder",
        anchor_residues={2: "I", 5: "F", 7: "N", 10: "I"},
        anchor_scores={2: 1.3, 5: 1.2, 7: 1.1, 10: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="VRYGRQQMMNNQAVV",
        binding_score=0.476868,
        ic50_nm=571.25,
        binding_class="weak_binder",
        anchor_residues={5: "Q", 8: "M", 10: "N", 13: "V"},
        anchor_scores={5: 0.9, 8: 1.3, 10: 1.1, 13: 1.7},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="AIKWLKNDHCGRASQ",
        binding_score=0.444344,
        ic50_nm=688.87,
        binding_class="weak_binder",
        anchor_residues={1: "I", 4: "L", 6: "N", 9: "C"},
        anchor_scores={1: 1.3, 4: 1.7, 6: 1.1, 9: 0.8},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="YTTLVTKPPSQEMML",
        binding_score=0.430959,
        ic50_nm=744.05,
        binding_class="weak_binder",
        anchor_residues={0: "Y", 3: "L", 5: "T", 8: "P"},
        anchor_scores={0: 1.7, 3: 1.7, 5: 1.6, 8: 0.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="VETWRYPAPISMTME",
        binding_score=0.414154,
        ic50_nm=819.63,
        binding_class="weak_binder",
        anchor_residues={5: "Y", 8: "P", 10: "S", 13: "M"},
        anchor_scores={5: 1.7, 8: 0.5, 10: 1.7, 13: 1.2},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DLMYCWMCCSPEHLA",
        binding_score=0.399759,
        ic50_nm=890.43,
        binding_class="weak_binder",
        anchor_residues={5: "W", 8: "C", 10: "P", 13: "L"},
        anchor_scores={5: 1.6, 8: 0.8, 10: 0.8, 13: 1.6},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="QSVELWYQIQMITSM",
        binding_score=0.399759,
        ic50_nm=890.43,
        binding_class="weak_binder",
        anchor_residues={5: "W", 8: "I", 10: "M", 13: "S"},
        anchor_scores={5: 1.6, 8: 1.6, 10: 0.8, 13: 0.8},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KWGPLWCLYRMMADT",
        binding_score=0.389286,
        ic50_nm=945.77,
        binding_class="weak_binder",
        anchor_residues={4: "L", 7: "L", 9: "R", 12: "A"},
        anchor_scores={4: 1.4, 7: 1.7, 9: 0.6, 12: 1.1},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="MEVYYWQMPNCIPHW",
        binding_score=0.381756,
        ic50_nm=987.66,
        binding_class="weak_binder",
        anchor_residues={3: "Y", 6: "Q", 8: "P", 11: "I"},
        anchor_scores={3: 1.7, 6: 0.8, 8: 0.8, 11: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="MTKPFDTVTWNDAVY",
        binding_score=0.375615,
        ic50_nm=1023.2,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "V", 9: "W", 12: "A"},
        anchor_scores={4: 1.8, 7: 1.5, 9: 0.5, 12: 1.1},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="QQSNVLHDTNMCEAG",
        binding_score=0.372984,
        ic50_nm=1038.82,
        binding_class="weak_binder",
        anchor_residues={2: "S", 5: "L", 7: "D", 10: "M"},
        anchor_scores={2: 0.9, 5: 1.7, 7: 0.8, 10: 1.2},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="CLFRIQPGCNEPAWN",
        binding_score=0.367209,
        ic50_nm=1073.93,
        binding_class="weak_binder",
        anchor_residues={1: "L", 4: "I", 6: "P", 9: "N"},
        anchor_scores={1: 1.4, 4: 1.6, 6: 0.8, 9: 0.8},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NVLEWNGMYRLCANW",
        binding_score=0.367209,
        ic50_nm=1073.93,
        binding_class="weak_binder",
        anchor_residues={2: "L", 5: "N", 7: "M", 10: "L"},
        anchor_scores={2: 1.4, 5: 0.8, 7: 0.8, 10: 1.6},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="CVVRNCICMNLNAHW",
        binding_score=0.349939,
        ic50_nm=1186.19,
        binding_class="weak_binder",
        anchor_residues={2: "V", 5: "C", 7: "C", 10: "L"},
        anchor_scores={2: 1.3, 5: 0.8, 7: 0.8, 10: 1.6},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="AHSFTWDHYFYLSVY",
        binding_score=0.345506,
        ic50_nm=1216.85,
        binding_class="weak_binder",
        anchor_residues={5: "W", 8: "Y", 10: "Y", 13: "V"},
        anchor_scores={5: 1.6, 8: 0.8, 10: 0.6, 13: 1.7},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PNSNFPWTGPMQAHR",
        binding_score=0.338763,
        ic50_nm=1265.0,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "T", 9: "P", 12: "A"},
        anchor_scores={4: 1.8, 7: 0.8, 9: 0.8, 12: 1.1},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GLNSFEYARMLYDLF",
        binding_score=0.328433,
        ic50_nm=1342.51,
        binding_class="weak_binder",
        anchor_residues={1: "L", 4: "F", 6: "Y", 9: "M"},
        anchor_scores={1: 1.4, 4: 1.2, 6: 0.6, 9: 1.2},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="SWQWWGEDNDRVATF",
        binding_score=0.3261,
        ic50_nm=1360.66,
        binding_class="weak_binder",
        anchor_residues={3: "W", 6: "E", 8: "N", 11: "V"},
        anchor_scores={3: 1.6, 6: 0.4, 8: 1.1, 11: 1.7},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="IMHDDIQHSMQMCCA",
        binding_score=0.294274,
        ic50_nm=1634.23,
        binding_class="weak_binder",
        anchor_residues={6: "Q", 9: "M", 11: "M", 14: "A"},
        anchor_scores={6: 0.9, 9: 1.3, 11: 0.8, 14: 1.1},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RPCVDSDKFADDDFV",
        binding_score=0.292364,
        ic50_nm=1652.3,
        binding_class="weak_binder",
        anchor_residues={0: "R", 3: "V", 5: "S", 8: "F"},
        anchor_scores={0: 0.5, 3: 1.5, 5: 1.7, 8: 0.8},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="SMDKLEKMKEIHRIE",
        binding_score=0.228891,
        ic50_nm=2381.07,
        binding_class="weak_binder",
        anchor_residues={4: "L", 7: "M", 9: "E", 12: "R"},
        anchor_scores={4: 1.4, 7: 1.3, 9: 0.8, 12: 0.5},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RLEKPHASGPWFNVG",
        binding_score=0.200051,
        ic50_nm=2811.07,
        binding_class="weak_binder",
        anchor_residues={5: "H", 8: "G", 10: "W", 13: "V"},
        anchor_scores={5: 0.9, 8: 0.8, 10: 0.5, 13: 1.7},
        source="pssm_predicted",
        peptide_length=15,
    ))

    # --- Non-binders ---
    entries.append(PrecomputedEntry(
        peptide="GKEEPENDQKNQEQD",
        binding_score=0.099264,
        ic50_nm=5021.58,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPPEPENENKKNDGK",
        binding_score=0.099264,
        ic50_nm=5021.58,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NKKPDKEDDQENRKH",
        binding_score=0.099264,
        ic50_nm=5021.58,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKENEHHNDRNPGDP",
        binding_score=0.099264,
        ic50_nm=5021.58,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EGPGRENHREDEGPR",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GRHDGEKKEHPQENE",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="REKNNPDKHGDKPEN",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GRNPERRQGQPDKHD",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GNDNDRKERGHKRKG",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KHRNKQRHRPQPNPG",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KGRGEKEHGQDQDKD",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GKEDRQRNEQRNGRK",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RNQEPKDKGKKGDNP",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="QRQNGKERPHGNDPR",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PKQPRRKNNKPEGNR",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPRRRRRQEPRPGNK",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DRKQQNDEEGDPKPQ",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRGQGKKHGHRRPKD",
        binding_score=0.08976,
        ic50_nm=5303.95,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="QEDKRNNQRNPQDDD",
        binding_score=0.086531,
        ic50_nm=5403.46,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RNEDDGNGDPHGQNK",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RPKPNPERGDRGPED",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PDPEPHQEDDDERED",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NGKDGPGEEPERDGH",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GDQHRHEREEEQHHQ",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NNRKHGPQEDKKEEN",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="REGRQDDPNRHGPNQ",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PRKHGQDEERKQHNK",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PHEQRDPQEPNRPKN",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PNKNEENDPPDEPQR",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NDNDDENEHDPDPPP",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NHRKEKPKKGPRHNK",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DGDHDRERHNERDDH",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GQPDQDKDRPENQKQ",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="QHHDEHGREGKGGHE",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KGEEKKGPQEEEHKN",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GDDKGGRHDKHQDRR",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DRNHPNNKRPDPPKD",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDKHNQKDEEEPHED",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKKKPQGERPKHNDR",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKENENEDERPNDKH",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KGKRPDNDRDNQRNQ",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PKNHDPQRPPEPRED",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NDHDDDHHGRHRRDE",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RQNEDEHRKGPRGDE",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RPQDNEGPKPGGKDK",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="ENHGPPHHEDDEDQR",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PGRDPDEDNKGKGQE",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GHEEDDKGGDKHPKK",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="REEGDEPNRQGHQRE",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RNERGPRKDKPRDQD",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GPDKDKHNDQDRPDK",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RDDNNNDHKGQPEPH",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GEDKGRGEPDPKNQR",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="QNNDKERKQPGHGHG",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NPDEDDPERGRHRQE",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RNRDHKEERPQHKEG",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DNRDHPRQNEDDEEK",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DGEEGGGEHPERQDQ",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DRQRKDKGQQGNDGP",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PHHQDRQRREHDRQQ",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DENKRQKPPHEGQRP",
        binding_score=0.079918,
        ic50_nm=5613.12,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="QRKKNGKQRRRNKPH",
        binding_score=0.075818,
        ic50_nm=5747.18,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKDGENEKRNNGKEE",
        binding_score=0.075818,
        ic50_nm=5747.18,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RERGRKQNRRGPPER",
        binding_score=0.075818,
        ic50_nm=5747.18,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NNPEPRDQHRRPKGR",
        binding_score=0.074225,
        ic50_nm=5800.13,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDDPKRQKRRHKPRE",
        binding_score=0.074225,
        ic50_nm=5800.13,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DRPHNKDRQRPEQDH",
        binding_score=0.074225,
        ic50_nm=5800.13,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RNDDNDRKKRHGDGE",
        binding_score=0.074225,
        ic50_nm=5800.13,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NGQDGKRRRRPKQHE",
        binding_score=0.074225,
        ic50_nm=5800.13,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GGDPEPEQPDDPQRD",
        binding_score=0.069715,
        ic50_nm=5952.68,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PEPKDPKQRGPENKK",
        binding_score=0.069715,
        ic50_nm=5952.68,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PNEKDDKHDEGNQPN",
        binding_score=0.069715,
        ic50_nm=5952.68,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NQGEDKKNKDPPHGH",
        binding_score=0.069715,
        ic50_nm=5952.68,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RDGEEREHHHRRNRG",
        binding_score=0.069715,
        ic50_nm=5952.68,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RDQKERKNPEDNDDR",
        binding_score=0.066942,
        ic50_nm=6048.47,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PNHPKEKPKEDEGKD",
        binding_score=0.059469,
        ic50_nm=6314.33,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DHKHDEKPKKPKDHH",
        binding_score=0.059469,
        ic50_nm=6314.33,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDNNEKEQPRKDGNR",
        binding_score=0.059469,
        ic50_nm=6314.33,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DGPNEKKKDKRPKNH",
        binding_score=0.059469,
        ic50_nm=6314.33,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RQHNEDKQERERNPR",
        binding_score=0.059469,
        ic50_nm=6314.33,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKKEDEQRKRNPRDR",
        binding_score=0.059469,
        ic50_nm=6314.33,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HDDPRKKPRRGDKRE",
        binding_score=0.059469,
        ic50_nm=6314.33,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NNGEKRQKDKDENRD",
        binding_score=0.059469,
        ic50_nm=6314.33,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DQKKERDDPKGNNNQ",
        binding_score=0.057783,
        ic50_nm=6375.93,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KEEENEQPEEERPEG",
        binding_score=0.056427,
        ic50_nm=6425.89,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRDQDPERRQKQEPP",
        binding_score=0.056427,
        ic50_nm=6425.89,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DGPPRERRRRQPQNR",
        binding_score=0.056427,
        ic50_nm=6425.89,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHDDKRHKQEDEREK",
        binding_score=0.056427,
        ic50_nm=6425.89,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NPQDDEREKPPQRGQ",
        binding_score=0.050943,
        ic50_nm=6631.96,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PDHEPRRKPGPHDRR",
        binding_score=0.050943,
        ic50_nm=6631.96,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDHKGKHPGPPKERE",
        binding_score=0.050943,
        ic50_nm=6631.96,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDRRDRKGQNGHDDD",
        binding_score=0.050943,
        ic50_nm=6631.96,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EEQRDKRDERPPNHN",
        binding_score=0.050943,
        ic50_nm=6631.96,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKPEEKRKRQDRERE",
        binding_score=0.043959,
        ic50_nm=6904.03,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GEDDRRNPPPQKEED",
        binding_score=0.043959,
        ic50_nm=6904.03,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HPGRKEEGDEDQDPK",
        binding_score=0.043959,
        ic50_nm=6904.03,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DEQRGDRRRKDPEND",
        binding_score=0.043959,
        ic50_nm=6904.03,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EQKNEPKHDRRPRER",
        binding_score=0.043959,
        ic50_nm=6904.03,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DGDEERGHGPDRRER",
        binding_score=0.039972,
        ic50_nm=7064.32,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKEEDEEDQHPRNDP",
        binding_score=0.036845,
        ic50_nm=7192.62,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))

    return entries


# Lazy singleton — built on first access
_database: PrecomputedAlleleDatabase | None = None


def get_database() -> PrecomputedAlleleDatabase:
    """Return the precomputed binding database for HLA-DRB1*01:01."""
    global _database
    if _database is None:
        _database = PrecomputedAlleleDatabase(
            allele=ALLELE,
            peptide_length=PEPTIDE_LENGTH,
            anchor_description=ANCHOR_DESCRIPTION,
            known_epitopes=KNOWN_EPITOPES,
            entries=_build_entries(),
        )
    return _database