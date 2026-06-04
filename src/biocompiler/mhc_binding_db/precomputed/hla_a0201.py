"""Precomputed MHC binding data for HLA-A*02:01.

HLA-A*02:01 is the most studied and most common MHC-I allele worldwide.
It has well-known binding motifs with primary anchor residues at
position 2 (L/M/V/I) and position 9 (V/L/I/A).

MHC-I presents 9-mer peptides.

Anchor residues: position 2 (L/M/V/I), position 9 (V/L/I/A)
Known epitopes: GILGFVFTL, LLFGYPVYV, YVLDHLIVV, MLGEQLFKA, ILDKVLVHL, SVYDFFVWL, KVAELVHFL
Data source: IEDB and SYFPEITHI databases + PSSM predictions

Binding scores are derived from PSSM-based prediction (SYFPEITHI
Rammensee et al. 1999) with expected AUC-ROC of 0.60-0.75.
Do not replace experimental validation or NetMHCpan predictions
where available.
"""

from __future__ import annotations

from .. import MHCBindingDatabase, PrecomputedEntry

ALLELE: str = "HLA-A*02:01"
MHC_CLASS: str = "I"
PEPTIDE_LENGTH: int = 9
ANCHOR_DESCRIPTION: str = "position 2 (L/M/V/I), position 9 (V/L/I/A)"
KNOWN_EPITOPES: list[str] = ['GILGFVFTL', 'LLFGYPVYV', 'YVLDHLIVV', 'MLGEQLFKA', 'ILDKVLVHL', 'SVYDFFVWL', 'KVAELVHFL']

# ============================================================================
# Position-Specific Scoring Matrix (PSSM) -- 9-mer
# Based on SYFPEITHI (Rammensee et al. 1999) binding motif data
# Anchor positions: P2 (L/M/V/I), P9 (V/L/I/A)
#
# Position 2: L(+9), M(+7), V(+5), I(+3), A(+1), T(0),
#             D/E/K/R/P/G strongly disfavored
# Position 9: V(+9), L(+7), I(+5), A(+2), T(+1),
#             D/E/K/R/P/G strongly disfavored
# Non-anchor positions: default=1.0 (permissive)
# ============================================================================

PSSM: list[dict[str, float]] = [
    {"A": 1.1, "C": 1.0, "D": 0.5, "E": 0.5, "F": 1.1, "G": 1.0, "H": 1.0, "I": 1.2, "K": 0.5, "L": 1.2, "M": 1.2, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 0.5, "S": 1.0, "T": 1.0, "V": 1.2, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 0.6, "D": 0.3, "E": 0.3, "F": 0.6, "G": 0.4, "H": 0.6, "I": 1.4, "K": 0.3, "L": 2.0, "M": 1.8, "N": 0.6, "P": 0.3, "Q": 0.6, "R": 0.3, "S": 0.6, "T": 0.9, "V": 1.6, "W": 0.6, "Y": 0.6},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.2, "C": 0.6, "D": 0.3, "E": 0.3, "F": 0.6, "G": 0.4, "H": 0.6, "I": 1.6, "K": 0.3, "L": 1.8, "M": 0.6, "N": 0.6, "P": 0.3, "Q": 0.6, "R": 0.3, "S": 0.6, "T": 0.9, "V": 2.0, "W": 0.6, "Y": 0.6},
]


def _build_entries() -> list[PrecomputedEntry]:
    """Construct all precomputed entries for HLA-A*02:01."""
    entries: list[PrecomputedEntry] = []

    # --- Binders (strong + moderate + weak) ---
    entries.append(PrecomputedEntry(
        peptide="LLFGYPVYV",
        binding_score=1.0,
        ic50_nm=28.12,
        binding_class="strong_binder",
        anchor_residues={0: "L", 1: "L", 8: "V"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 2.0},
        source="known_epitope",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ILRLQVSFV",
        binding_score=1.0,
        ic50_nm=28.12,
        binding_class="strong_binder",
        anchor_residues={0: "I", 1: "L", 8: "V"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ALVSRDYCV",
        binding_score=0.953029,
        ic50_nm=36.85,
        binding_class="strong_binder",
        anchor_residues={0: "A", 1: "L", 8: "V"},
        anchor_scores={0: 1.1, 1: 2.0, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ILDKVLVHL",
        binding_score=0.943324,
        ic50_nm=38.97,
        binding_class="strong_binder",
        anchor_residues={0: "I", 1: "L", 8: "L"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 1.8},
        source="known_epitope",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IMKIIFEVV",
        binding_score=0.943324,
        ic50_nm=38.97,
        binding_class="strong_binder",
        anchor_residues={0: "I", 1: "M", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MLPKRNGNL",
        binding_score=0.943324,
        ic50_nm=38.97,
        binding_class="strong_binder",
        anchor_residues={0: "M", 1: "L", 8: "L"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MMPNRTREV",
        binding_score=0.943324,
        ic50_nm=38.97,
        binding_class="strong_binder",
        anchor_residues={0: "M", 1: "M", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VMVAEDFVV",
        binding_score=0.943324,
        ic50_nm=38.97,
        binding_class="strong_binder",
        anchor_residues={0: "V", 1: "M", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LMCNHIENV",
        binding_score=0.943324,
        ic50_nm=38.97,
        binding_class="strong_binder",
        anchor_residues={0: "L", 1: "M", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IMTEPWHKV",
        binding_score=0.943324,
        ic50_nm=38.97,
        binding_class="strong_binder",
        anchor_residues={0: "I", 1: "M", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MMTYLEHLV",
        binding_score=0.943324,
        ic50_nm=38.97,
        binding_class="strong_binder",
        anchor_residues={0: "M", 1: "M", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IMRDLIKML",
        binding_score=0.888933,
        ic50_nm=53.29,
        binding_class="moderate_binder",
        anchor_residues={0: "I", 1: "M", 8: "L"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LMILLRDIL",
        binding_score=0.888933,
        ic50_nm=53.29,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "M", 8: "L"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MMRGLACML",
        binding_score=0.888933,
        ic50_nm=53.29,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "M", 8: "L"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ILQHEEQNI",
        binding_score=0.882667,
        ic50_nm=55.25,
        binding_class="moderate_binder",
        anchor_residues={0: "I", 1: "L", 8: "I"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VVADIGQSV",
        binding_score=0.882667,
        ic50_nm=55.25,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "V", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.6, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LVHTKFNDV",
        binding_score=0.882667,
        ic50_nm=55.25,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "V", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.6, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ILGQAGMQI",
        binding_score=0.882667,
        ic50_nm=55.25,
        binding_class="moderate_binder",
        anchor_residues={0: "I", 1: "L", 8: "I"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VLIELYEWI",
        binding_score=0.882667,
        ic50_nm=55.25,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "L", 8: "I"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LLPIQCMSI",
        binding_score=0.882667,
        ic50_nm=55.25,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "L", 8: "I"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IVMAGKCFV",
        binding_score=0.882667,
        ic50_nm=55.25,
        binding_class="moderate_binder",
        anchor_residues={0: "I", 1: "V", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.6, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MLMKNFSDI",
        binding_score=0.882667,
        ic50_nm=55.25,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "L", 8: "I"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="CLSTTGCTL",
        binding_score=0.850609,
        ic50_nm=66.45,
        binding_class="moderate_binder",
        anchor_residues={0: "C", 1: "L", 8: "L"},
        anchor_scores={0: 1.0, 1: 2.0, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HMIKPFLRV",
        binding_score=0.850609,
        ic50_nm=66.45,
        binding_class="moderate_binder",
        anchor_residues={0: "H", 1: "M", 8: "V"},
        anchor_scores={0: 1.0, 1: 1.8, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AVMPKQKDV",
        binding_score=0.839638,
        ic50_nm=70.78,
        binding_class="moderate_binder",
        anchor_residues={0: "A", 1: "V", 8: "V"},
        anchor_scores={0: 1.1, 1: 1.6, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VMKDHWMHI",
        binding_score=0.830753,
        ic50_nm=74.49,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "M", 8: "I"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VVPWWADQL",
        binding_score=0.830753,
        ic50_nm=74.49,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "V", 8: "L"},
        anchor_scores={0: 1.2, 1: 1.6, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IMIHAYFII",
        binding_score=0.830753,
        ic50_nm=74.49,
        binding_class="moderate_binder",
        anchor_residues={0: "I", 1: "M", 8: "I"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MVDPAKVEL",
        binding_score=0.830753,
        ic50_nm=74.49,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "V", 8: "L"},
        anchor_scores={0: 1.2, 1: 1.6, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VIGRPKVIV",
        binding_score=0.817239,
        ic50_nm=80.52,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "I", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LIHHVRFQV",
        binding_score=0.817239,
        ic50_nm=80.52,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "I", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VIQYYYITV",
        binding_score=0.817239,
        ic50_nm=80.52,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "I", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VIDYPYISV",
        binding_score=0.817239,
        ic50_nm=80.52,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "I", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MISEIPWNV",
        binding_score=0.817239,
        ic50_nm=80.52,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "I", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LIYQEEVHV",
        binding_score=0.817239,
        ic50_nm=80.52,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "I", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YVLDHLIVV",
        binding_score=0.794194,
        ic50_nm=91.94,
        binding_class="moderate_binder",
        anchor_residues={0: "Y", 1: "V", 8: "V"},
        anchor_scores={0: 1.0, 1: 1.6, 8: 2.0},
        source="known_epitope",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LVISSRADI",
        binding_score=0.77526,
        ic50_nm=102.53,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "V", 8: "I"},
        anchor_scores={0: 1.2, 1: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VVCWVVMEI",
        binding_score=0.77526,
        ic50_nm=102.53,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "V", 8: "I"},
        anchor_scores={0: 1.2, 1: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IVNEWTHFI",
        binding_score=0.77526,
        ic50_nm=102.53,
        binding_class="moderate_binder",
        anchor_residues={0: "I", 1: "V", 8: "I"},
        anchor_scores={0: 1.2, 1: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VIVQIRWKL",
        binding_score=0.76804,
        ic50_nm=106.88,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "I", 8: "L"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MIYYEGLEL",
        binding_score=0.76804,
        ic50_nm=106.88,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "I", 8: "L"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VIVGDRNWL",
        binding_score=0.76804,
        ic50_nm=106.88,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "I", 8: "L"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MLGEQLFKA",
        binding_score=0.745963,
        ic50_nm=121.37,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "L", 8: "A"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 1.2},
        source="known_epitope",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SVYDFFVWL",
        binding_score=0.745963,
        ic50_nm=121.37,
        binding_class="moderate_binder",
        anchor_residues={0: "S", 1: "V", 8: "L"},
        anchor_scores={0: 1.0, 1: 1.6, 8: 1.8},
        source="known_epitope",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MLFSEAWLA",
        binding_score=0.745963,
        ic50_nm=121.37,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "L", 8: "A"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LLLIYNIHA",
        binding_score=0.745963,
        ic50_nm=121.37,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "L", 8: "A"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AVEFKEEVI",
        binding_score=0.735942,
        ic50_nm=128.57,
        binding_class="moderate_binder",
        anchor_residues={0: "A", 1: "V", 8: "I"},
        anchor_scores={0: 1.1, 1: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AVDDMYFDI",
        binding_score=0.735942,
        ic50_nm=128.57,
        binding_class="moderate_binder",
        anchor_residues={0: "A", 1: "V", 8: "I"},
        anchor_scores={0: 1.1, 1: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SIYVMYIDV",
        binding_score=0.733416,
        ic50_nm=130.46,
        binding_class="moderate_binder",
        anchor_residues={0: "S", 1: "I", 8: "V"},
        anchor_scores={0: 1.0, 1: 1.4, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="CIPYCAHLV",
        binding_score=0.733416,
        ic50_nm=130.46,
        binding_class="moderate_binder",
        anchor_residues={0: "C", 1: "I", 8: "V"},
        anchor_scores={0: 1.0, 1: 1.4, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FICISYDRL",
        binding_score=0.728974,
        ic50_nm=133.83,
        binding_class="moderate_binder",
        anchor_residues={0: "F", 1: "I", 8: "L"},
        anchor_scores={0: 1.1, 1: 1.4, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MITPVMAEI",
        binding_score=0.715491,
        ic50_nm=144.63,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "I", 8: "I"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IIEQIGTRI",
        binding_score=0.715491,
        ic50_nm=144.63,
        binding_class="moderate_binder",
        anchor_residues={0: "I", 1: "I", 8: "I"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MIETRAFQI",
        binding_score=0.715491,
        ic50_nm=144.63,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "I", 8: "I"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MIRQKHTEI",
        binding_score=0.715491,
        ic50_nm=144.63,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "I", 8: "I"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ALNLMWYDA",
        binding_score=0.707675,
        ic50_nm=151.29,
        binding_class="moderate_binder",
        anchor_residues={0: "A", 1: "L", 8: "A"},
        anchor_scores={0: 1.1, 1: 2.0, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VMIDCMDTA",
        binding_score=0.699776,
        ic50_nm=158.33,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "M", 8: "A"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MMMELTLQA",
        binding_score=0.699776,
        ic50_nm=158.33,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "M", 8: "A"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MMSHNKMKA",
        binding_score=0.699776,
        ic50_nm=158.33,
        binding_class="moderate_binder",
        anchor_residues={0: "M", 1: "M", 8: "A"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LMMNRKLKA",
        binding_score=0.699776,
        ic50_nm=158.33,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "M", 8: "A"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IMMNAGFWA",
        binding_score=0.699776,
        ic50_nm=158.33,
        binding_class="moderate_binder",
        anchor_residues={0: "I", 1: "M", 8: "A"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GILGFVFTL",
        binding_score=0.687767,
        ic50_nm=169.66,
        binding_class="moderate_binder",
        anchor_residues={0: "G", 1: "I", 8: "L"},
        anchor_scores={0: 1.0, 1: 1.4, 8: 1.8},
        source="known_epitope",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="WIQVYYVAL",
        binding_score=0.687767,
        ic50_nm=169.66,
        binding_class="moderate_binder",
        anchor_residues={0: "W", 1: "I", 8: "L"},
        anchor_scores={0: 1.0, 1: 1.4, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LAYIGLVAV",
        binding_score=0.667299,
        ic50_nm=190.88,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "A", 8: "V"},
        anchor_scores={0: 1.2, 1: 1.0, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LVKCDQKCA",
        binding_score=0.650495,
        ic50_nm=210.26,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "V", 8: "A"},
        anchor_scores={0: 1.2, 1: 1.6, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VTIYCTLRV",
        binding_score=0.624517,
        ic50_nm=244.18,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "T", 8: "V"},
        anchor_scores={0: 1.2, 1: 0.9, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AVWLPVALA",
        binding_score=0.615638,
        ic50_nm=256.98,
        binding_class="moderate_binder",
        anchor_residues={0: "A", 1: "V", 8: "A"},
        anchor_scores={0: 1.1, 1: 1.6, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AVHTNYSRA",
        binding_score=0.615638,
        ic50_nm=256.98,
        binding_class="moderate_binder",
        anchor_residues={0: "A", 1: "V", 8: "A"},
        anchor_scores={0: 1.1, 1: 1.6, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VIWKMIDKA",
        binding_score=0.59753,
        ic50_nm=285.22,
        binding_class="moderate_binder",
        anchor_residues={0: "V", 1: "I", 8: "A"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LIHHKDGIA",
        binding_score=0.59753,
        ic50_nm=285.22,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 1: "I", 8: "A"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FIFDCGLYA",
        binding_score=0.564626,
        ic50_nm=344.69,
        binding_class="moderate_binder",
        anchor_residues={0: "F", 1: "I", 8: "A"},
        anchor_scores={0: 1.1, 1: 1.4, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RMFSLTKQL",
        binding_score=0.517172,
        ic50_nm=452.97,
        binding_class="moderate_binder",
        anchor_residues={0: "R", 1: "M", 8: "L"},
        anchor_scores={0: 0.5, 1: 1.8, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KVAELVHFL",
        binding_score=0.477037,
        ic50_nm=570.7,
        binding_class="weak_binder",
        anchor_residues={0: "K", 1: "V", 8: "L"},
        anchor_scores={0: 0.5, 1: 1.6, 8: 1.8},
        source="known_epitope",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ILEDPENNY",
        binding_score=0.477037,
        ic50_nm=570.7,
        binding_class="weak_binder",
        anchor_residues={0: "I", 1: "L", 8: "Y"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VLQQMVRQC",
        binding_score=0.477037,
        ic50_nm=570.7,
        binding_class="weak_binder",
        anchor_residues={0: "V", 1: "L", 8: "C"},
        anchor_scores={0: 1.2, 1: 2.0, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IQWPNIFTL",
        binding_score=0.442927,
        ic50_nm=694.52,
        binding_class="weak_binder",
        anchor_residues={0: "I", 1: "Q", 8: "L"},
        anchor_scores={0: 1.2, 1: 0.6, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LMSGDLTYM",
        binding_score=0.442927,
        ic50_nm=694.52,
        binding_class="weak_binder",
        anchor_residues={0: "L", 1: "M", 8: "M"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IMWHDVHWH",
        binding_score=0.442927,
        ic50_nm=694.52,
        binding_class="weak_binder",
        anchor_residues={0: "I", 1: "M", 8: "H"},
        anchor_scores={0: 1.2, 1: 1.8, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DIQMTKALL",
        binding_score=0.434087,
        ic50_nm=730.77,
        binding_class="weak_binder",
        anchor_residues={0: "D", 1: "I", 8: "L"},
        anchor_scores={0: 0.5, 1: 1.4, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="QYDPPYRTV",
        binding_score=0.419051,
        ic50_nm=796.84,
        binding_class="weak_binder",
        anchor_residues={0: "Q", 1: "Y", 8: "V"},
        anchor_scores={0: 1.0, 1: 0.6, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LSYQKCNHI",
        binding_score=0.406735,
        ic50_nm=855.39,
        binding_class="weak_binder",
        anchor_residues={0: "L", 1: "S", 8: "I"},
        anchor_scores={0: 1.2, 1: 0.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AWNFDLMQI",
        binding_score=0.381275,
        ic50_nm=990.4,
        binding_class="weak_binder",
        anchor_residues={0: "A", 1: "W", 8: "I"},
        anchor_scores={0: 1.1, 1: 0.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ANCPKHERI",
        binding_score=0.381275,
        ic50_nm=990.4,
        binding_class="weak_binder",
        anchor_residues={0: "A", 1: "N", 8: "I"},
        anchor_scores={0: 1.1, 1: 0.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MIQMMEGMQ",
        binding_score=0.368098,
        ic50_nm=1068.45,
        binding_class="weak_binder",
        anchor_residues={0: "M", 1: "I", 8: "Q"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LYTQVRGSA",
        binding_score=0.326541,
        ic50_nm=1357.21,
        binding_class="weak_binder",
        anchor_residues={0: "L", 1: "Y", 8: "A"},
        anchor_scores={0: 1.2, 1: 0.6, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VCRDMKMEA",
        binding_score=0.326541,
        ic50_nm=1357.21,
        binding_class="weak_binder",
        anchor_residues={0: "V", 1: "C", 8: "A"},
        anchor_scores={0: 1.2, 1: 0.6, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GISLMERDF",
        binding_score=0.319289,
        ic50_nm=1415.07,
        binding_class="weak_binder",
        anchor_residues={0: "G", 1: "I", 8: "F"},
        anchor_scores={0: 1.0, 1: 1.4, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AQSEQNRFA",
        binding_score=0.304475,
        ic50_nm=1541.03,
        binding_class="weak_binder",
        anchor_residues={0: "A", 1: "Q", 8: "A"},
        anchor_scores={0: 1.1, 1: 0.6, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="APSANLPQV",
        binding_score=0.261367,
        ic50_nm=1975.07,
        binding_class="weak_binder",
        anchor_residues={0: "A", 1: "P", 8: "V"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 2.0},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VKMELQGRL",
        binding_score=0.257256,
        ic50_nm=2022.36,
        binding_class="weak_binder",
        anchor_residues={0: "V", 1: "K", 8: "L"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AMSEEIVFP",
        binding_score=0.238317,
        ic50_nm=2255.31,
        binding_class="weak_binder",
        anchor_residues={0: "A", 1: "M", 8: "P"},
        anchor_scores={0: 1.1, 1: 1.8, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IVKFFVKGE",
        binding_score=0.231835,
        ic50_nm=2341.06,
        binding_class="weak_binder",
        anchor_residues={0: "I", 1: "V", 8: "E"},
        anchor_scores={0: 1.2, 1: 1.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AVWIWYCYD",
        binding_score=0.214104,
        ic50_nm=2592.63,
        binding_class="weak_binder",
        anchor_residues={0: "A", 1: "V", 8: "D"},
        anchor_scores={0: 1.1, 1: 1.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VIKQSAPMP",
        binding_score=0.20498,
        ic50_nm=2732.44,
        binding_class="weak_binder",
        anchor_residues={0: "V", 1: "I", 8: "P"},
        anchor_scores={0: 1.2, 1: 1.4, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LFFSVGKTY",
        binding_score=0.176475,
        ic50_nm=3219.68,
        binding_class="weak_binder",
        anchor_residues={0: "L", 1: "F", 8: "Y"},
        anchor_scores={0: 1.2, 1: 0.6, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MHRNLIIAH",
        binding_score=0.176475,
        ic50_nm=3219.68,
        binding_class="weak_binder",
        anchor_residues={0: "M", 1: "H", 8: "H"},
        anchor_scores={0: 1.2, 1: 0.6, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MPMSVCYDA",
        binding_score=0.176475,
        ic50_nm=3219.68,
        binding_class="weak_binder",
        anchor_residues={0: "M", 1: "P", 8: "A"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AELTYHFNA",
        binding_score=0.161524,
        ic50_nm=3509.06,
        binding_class="weak_binder",
        anchor_residues={0: "A", 1: "E", 8: "A"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 1.2},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KICGSTRKG",
        binding_score=0.109669,
        ic50_nm=4729.63,
        binding_class="weak_binder",
        anchor_residues={0: "K", 1: "I", 8: "G"},
        anchor_scores={0: 0.5, 1: 1.4, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDPSDWCFL",
        binding_score=0.10491,
        ic50_nm=4860.99,
        binding_class="weak_binder",
        anchor_residues={0: "E", 1: "D", 8: "L"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))

    # --- Non-binders ---
    entries.append(PrecomputedEntry(
        peptide="RDQMWETCI",
        binding_score=0.090342,
        ic50_nm=5286.21,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "D", 8: "I"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RVWCGMAHR",
        binding_score=0.090342,
        ic50_nm=5286.21,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "V", 8: "R"},
        anchor_scores={0: 0.5, 1: 1.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RVYDEINGP",
        binding_score=0.090342,
        ic50_nm=5286.21,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "V", 8: "P"},
        anchor_scores={0: 0.5, 1: 1.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KEFQCKFFI",
        binding_score=0.090342,
        ic50_nm=5286.21,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "E", 8: "I"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EVYYYLALD",
        binding_score=0.090342,
        ic50_nm=5286.21,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "V", 8: "D"},
        anchor_scores={0: 0.5, 1: 1.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DVMPWDLQR",
        binding_score=0.090342,
        ic50_nm=5286.21,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "V", 8: "R"},
        anchor_scores={0: 0.5, 1: 1.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YGTFRGGWW",
        binding_score=0.090342,
        ic50_nm=5286.21,
        binding_class="non_binder",
        anchor_residues={0: "Y", 1: "G", 8: "W"},
        anchor_scores={0: 1.0, 1: 0.4, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VPSHMLAHH",
        binding_score=0.078377,
        ic50_nm=5663.14,
        binding_class="non_binder",
        anchor_residues={0: "V", 1: "P", 8: "H"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LCKIVLGRD",
        binding_score=0.078377,
        ic50_nm=5663.14,
        binding_class="non_binder",
        anchor_residues={0: "L", 1: "C", 8: "D"},
        anchor_scores={0: 1.2, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IWGPLLFGR",
        binding_score=0.078377,
        ic50_nm=5663.14,
        binding_class="non_binder",
        anchor_residues={0: "I", 1: "W", 8: "R"},
        anchor_scores={0: 1.2, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RIRNAQCPD",
        binding_score=0.075344,
        ic50_nm=5762.88,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "I", 8: "D"},
        anchor_scores={0: 0.5, 1: 1.4, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AKPIRKMLW",
        binding_score=0.069231,
        ic50_nm=5969.28,
        binding_class="non_binder",
        anchor_residues={0: "A", 1: "K", 8: "W"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AEHRLDQSF",
        binding_score=0.069231,
        ic50_nm=5969.28,
        binding_class="non_binder",
        anchor_residues={0: "A", 1: "E", 8: "F"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FDFMYCYLS",
        binding_score=0.069231,
        ic50_nm=5969.28,
        binding_class="non_binder",
        anchor_residues={0: "F", 1: "D", 8: "S"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GNNLPQMYD",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "G", 1: "N", 8: "D"},
        anchor_scores={0: 1.0, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="QDRPMGSSN",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "Q", 1: "D", 8: "N"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="WDFGYYMWC",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "W", 1: "D", 8: "C"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HYTPTCCCE",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "H", 1: "Y", 8: "E"},
        anchor_scores={0: 1.0, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HWVFGMTRR",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "H", 1: "W", 8: "R"},
        anchor_scores={0: 1.0, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GFQNPRPPK",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "G", 1: "F", 8: "K"},
        anchor_scores={0: 1.0, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GPTGVGMVQ",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "G", 1: "P", 8: "Q"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="PSFWKQPRR",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "P", 1: "S", 8: "R"},
        anchor_scores={0: 1.0, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NCSWGNGKR",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "N", 1: "C", 8: "R"},
        anchor_scores={0: 1.0, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="WCQFILKMK",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "W", 1: "C", 8: "K"},
        anchor_scores={0: 1.0, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SEKQPPCWW",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "S", 1: "E", 8: "W"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SFHALQETD",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "S", 1: "F", 8: "D"},
        anchor_scores={0: 1.0, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KYELEESHY",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "Y", 8: "Y"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SRGEFCGEQ",
        binding_score=0.059955,
        ic50_nm=6296.68,
        binding_class="non_binder",
        anchor_residues={0: "S", 1: "R", 8: "Q"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDDESYVAT",
        binding_score=0.036442,
        ic50_nm=7209.32,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "D", 8: "T"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.9},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ETLIWMWYP",
        binding_score=0.036442,
        ic50_nm=7209.32,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "T", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.9, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KTNVYYLLE",
        binding_score=0.036442,
        ic50_nm=7209.32,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "T", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.9, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AGHHCICRP",
        binding_score=0.034878,
        ic50_nm=7274.55,
        binding_class="non_binder",
        anchor_residues={0: "A", 1: "G", 8: "P"},
        anchor_scores={0: 1.1, 1: 0.4, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FDGIYKHEG",
        binding_score=0.034878,
        ic50_nm=7274.55,
        binding_class="non_binder",
        anchor_residues={0: "F", 1: "D", 8: "G"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KGRQSRHMY",
        binding_score=0.028662,
        ic50_nm=7539.54,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "G", 8: "Y"},
        anchor_scores={0: 0.5, 1: 0.4, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="QGQSPNVPE",
        binding_score=0.028662,
        ic50_nm=7539.54,
        binding_class="non_binder",
        anchor_residues={0: "Q", 1: "G", 8: "E"},
        anchor_scores={0: 1.0, 1: 0.4, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NEVYHHLRG",
        binding_score=0.028662,
        ic50_nm=7539.54,
        binding_class="non_binder",
        anchor_residues={0: "N", 1: "E", 8: "G"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VPLKYFMDP",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "V", 1: "P", 8: "P"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MRMFFCLSD",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "M", 1: "R", 8: "D"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VDWLSANMK",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "V", 1: "D", 8: "K"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IDRNFKHED",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "I", 1: "D", 8: "D"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IDPVMPFDD",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "I", 1: "D", 8: "D"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IPAPNSVSD",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "I", 1: "P", 8: "D"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MPAFWFETK",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "M", 1: "P", 8: "K"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MRAKHIDNP",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "M", 1: "R", 8: "P"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VDMYYSAQE",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "V", 1: "D", 8: "E"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LDSHEAHGP",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "L", 1: "D", 8: "P"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VKNAMTHAD",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "V", 1: "K", 8: "D"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LRVHRMNCK",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "L", 1: "R", 8: "K"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IKPMLRFFR",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "I", 1: "K", 8: "R"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VKVDLDKEE",
        binding_score=0.022565,
        ic50_nm=7808.87,
        binding_class="non_binder",
        anchor_residues={0: "V", 1: "K", 8: "E"},
        anchor_scores={0: 1.2, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ADHFVYTQP",
        binding_score=0.018124,
        ic50_nm=8011.08,
        binding_class="non_binder",
        anchor_residues={0: "A", 1: "D", 8: "P"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FPYWIICTP",
        binding_score=0.018124,
        ic50_nm=8011.08,
        binding_class="non_binder",
        anchor_residues={0: "F", 1: "P", 8: "P"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FDTWAYGRP",
        binding_score=0.018124,
        ic50_nm=8011.08,
        binding_class="non_binder",
        anchor_residues={0: "F", 1: "D", 8: "P"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FEACIFHPR",
        binding_score=0.018124,
        ic50_nm=8011.08,
        binding_class="non_binder",
        anchor_residues={0: "F", 1: "E", 8: "R"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FKCAWYPVR",
        binding_score=0.018124,
        ic50_nm=8011.08,
        binding_class="non_binder",
        anchor_residues={0: "F", 1: "K", 8: "R"},
        anchor_scores={0: 1.1, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KQEIQWPTE",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "Q", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EEETTHNNF",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "E", 8: "F"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KWCNTDLRD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "W", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YETPRWVCP",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "Y", 1: "E", 8: "P"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RFCREMDTR",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "F", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ECMYFYPFK",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "C", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NDANDNIEK",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "N", 1: "D", 8: "K"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TKVSKCHLK",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "T", 1: "K", 8: "K"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GRNRCNYQP",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "G", 1: "R", 8: "P"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SPPVEKKRE",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "S", 1: "P", 8: "E"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NKVWGFQTK",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "N", 1: "K", 8: "K"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KNTLDKHVE",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "N", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RDRIHWNCC",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "D", 8: "C"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YKFMYMGPR",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "Y", 1: "K", 8: "R"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="WEICWSGTD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "W", 1: "E", 8: "D"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RSMSNVKLK",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "S", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKAVGQSSH",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "K", 8: "H"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="CPLWYEMLR",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "C", 1: "P", 8: "R"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="CERLGPTLE",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "C", 1: "E", 8: "E"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HRYDTRTND",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "H", 1: "R", 8: "D"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KCRFTQRWD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "C", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HPDCQNDVP",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "H", 1: "P", 8: "P"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EYCRYCICP",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "Y", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKCVWGQVS",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "K", 8: "S"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DRVHCKPYY",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "R", 8: "Y"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SDGNGFVSD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "S", 1: "D", 8: "D"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRQFIDFAH",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "R", 8: "H"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SKAVVQAAP",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "S", 1: "K", 8: "P"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DDDNKISYY",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "D", 8: "Y"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EEPMRMQYF",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "E", 8: "F"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="PDRYSWPVD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "P", 1: "D", 8: "D"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RFTPQYFWP",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "F", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKNDREIHW",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "K", 8: "W"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KNRSMWADR",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "N", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKEAHSDFW",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "K", 8: "W"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ECITQNRDR",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "C", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KERKEDGLY",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "E", 8: "Y"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ESQVSCDKP",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "S", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DFERGGHHE",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "F", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DWTQEQFFD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "W", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EENQKPSWY",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "E", 8: "Y"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPPMICSKN",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "P", 8: "N"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EEFLRMKQY",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "E", 8: "Y"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KQPNRKHRD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "Q", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKLHFTICY",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "K", 8: "Y"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KNTLGVGLE",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "N", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DYDSHTEPE",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "Y", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HPWGSDFYD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "H", 1: "P", 8: "D"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NEGPAIIFD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "N", 1: "E", 8: "D"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DDLGCVSHS",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "D", 8: "S"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YKAGTSSGD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "Y", 1: "K", 8: "D"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DPERRWTFS",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "P", 8: "S"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SPGWGFEPR",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "S", 1: "P", 8: "R"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DEQATEAPF",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "E", 8: "F"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EFCLHVCNE",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "F", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPCPPQMVC",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "P", 8: "C"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="REIGTDTGQ",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "E", 8: "Q"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKQGAMIIQ",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "K", 8: "Q"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RCMETKKWD",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "C", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.6, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TKPITPNNE",
        binding_score=0.013857,
        ic50_nm=8210.26,
        binding_class="non_binder",
        anchor_residues={0: "T", 1: "K", 8: "E"},
        anchor_scores={0: 1.0, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KGKDNKDNG",
        binding_score=0.009422,
        ic50_nm=8422.56,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "G", 8: "G"},
        anchor_scores={0: 0.5, 1: 0.4, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPCQKVKFG",
        binding_score=0.002281,
        ic50_nm=8776.02,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "P", 8: "G"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KGVIFPVLR",
        binding_score=0.002281,
        ic50_nm=8776.02,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "G", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.4, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DGIDTKHFK",
        binding_score=0.002281,
        ic50_nm=8776.02,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "G", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.4, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDQGFPVIP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "D", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRWQWPQLK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "R", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDSRWWRWK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "D", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RDDLVPMFE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "D", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERCSNVETP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "R", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DPQGDTNDP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "P", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERSGFDRCK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "R", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRYGVLMED",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "R", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRGAFYRCE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "R", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDNRDWFTE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "D", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRTIEHFID",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "R", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPHWQYASD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "P", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RDNTVEMIP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "D", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKIGMWMHD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "K", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKKNPNCPR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "K", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKKYYEFEP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "K", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERTHVMSTK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "R", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DPIVCGQGK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "P", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DEHGYPQTE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "E", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERVPTFIAP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "R", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKPIVMLPP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "K", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EESMHNMQR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "E", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKQLFHMIP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "K", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPETFIHYP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "P", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RDIILDCDP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "D", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPYYPMAKR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "P", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDYLIDDKD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "D", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRLRNKEFD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "R", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="REWTGPGFP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "E", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERSWIMPKR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "R", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKEKITTVR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "K", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKNASELQE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "K", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DKNSERRMD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "K", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDPAYHRWR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "D", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EENNRYGLE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "E", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RDHPVTSQE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "D", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="REKVNVCRD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "E", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKYFHHCWK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "K", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="REERAEGRR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "E", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KESDIPDEE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "E", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="REHDEFEWD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "E", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DEWQNVGSD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "E", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKLTQQPDE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "K", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKLGEAFAE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "K", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDAMQKTDR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "D", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DPAICSNPE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "P", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRHFWEFEE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "R", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDCLKGAPR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "D", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRADHRFDP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "R", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERPQGNHGK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "R", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERIKVCIWK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "R", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDNWWFWHR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "D", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DKYMYIWED",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "K", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERDCLAMKD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "R", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKVVKLLID",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "K", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RPKDVPGNE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "P", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRCFEDNVP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "R", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DRVSWISLD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "R", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKEAINSYK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "K", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDINILQMR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "D", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDWRFMGDK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "D", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KEWLGDCLP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "E", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DPVYCYCEK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "P", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDFDTDEKE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "D", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RPKDNMGEP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "P", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RDPNYAFSK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "D", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EPWWIMTTR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "P", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RPIVPPIDD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "R", 1: "P", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DEEIAMNLP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "E", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKYHHEPHP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "K", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRYYICMPE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "R", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DRTRSCLHD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "R", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDQEHVQDK",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "D", 8: "K"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DDWIKQPRR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "D", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DPCDNTELD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "P", 8: "D"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERYDHMSRR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "R", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRNLSVEGR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "K", 1: "R", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKDDATAWR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "K", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DRFNPIWLE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "D", 1: "R", 8: "E"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKCISIDKP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "K", 8: "P"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EESSHWGQR",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={0: "E", 1: "E", 8: "R"},
        anchor_scores={0: 0.5, 1: 0.3, 8: 0.3},
        source="pssm_predicted",
        peptide_length=9,
    ))

    return entries


# Lazy singleton -- built on first access
_database: MHCBindingDatabase | None = None


def get_database() -> MHCBindingDatabase:
    """Return the precomputed binding database for HLA-A*02:01."""
    global _database
    if _database is None:
        _database = MHCBindingDatabase(
            allele=ALLELE,
            peptide_length=PEPTIDE_LENGTH,
            anchor_description=ANCHOR_DESCRIPTION,
            known_epitopes=KNOWN_EPITOPES,
            entries=_build_entries(),
        )
    return _database


# ============================================================================
# Dict export for quick lookup
# ============================================================================

HLA_A0201_BINDERS: dict[str, PrecomputedEntry] = {
    entry.peptide: entry for entry in _build_entries()
}


def get_hla_a0201_database() -> MHCBindingDatabase:
    """Load and return the HLA-A*02:01 precomputed binding database.

    This is an alias for :func:`get_database` that provides a
    more explicit, allele-specific name.

    Returns
    -------
    MHCBindingDatabase
        The precomputed binding database for HLA-A*02:01.
    """
    return get_database()
