"""Precomputed MHC binding data for H-2Kb.

MHC-I presents 8-mer peptides.

Anchor residues: position 5 (F/Y), position 8 (V/L/I)
Known epitopes: SIINFEKL, FAPGNYPAL, SSYSYSSSY
Data source: IEDB and SYFPEITHI databases + PSSM predictions

Binding scores are derived from PSSM-based prediction with expected
AUC-ROC of 0.60-0.75. Do not replace experimental validation or
NetMHCpan predictions where available.

References
----------
- Falk et al., Nature 1991; 351:290 (H-2Kb motif)
- Rammensee et al., Immunogenetics 1995 (SYFPEITHI)
"""
from __future__ import annotations

from .. import PrecomputedAlleleDatabase, PrecomputedEntry

ALLELE: str = "H-2Kb"
MHC_CLASS: str = "I"
PEPTIDE_LENGTH: int = 8
ANCHOR_DESCRIPTION: str = "position 5 (F/Y), position 8 (V/L/I)"
KNOWN_EPITOPES: list[str] = ['SIINFEKL', 'FAPGNYPAL', 'SSYSYSSSY']

# ═══════════════════════════════════════════════════════════════════════════
# Position-Specific Scoring Matrix (PSSM) — core 9-mer
# ═══════════════════════════════════════════════════════════════════════════

PSSM: list[dict[str, float]] = [
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.4, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.3, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 0.8, "C": 0.8, "D": 0.4, "E": 0.4, "F": 1.8, "G": 0.8, "H": 0.8, "I": 0.8, "K": 0.4, "L": 0.8, "M": 0.8, "N": 0.8, "P": 0.8, "Q": 0.8, "R": 0.4, "S": 0.8, "T": 0.8, "V": 0.8, "W": 0.8, "Y": 1.7},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.2, "E": 1.3, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 0.8, "C": 0.8, "D": 0.4, "E": 0.4, "F": 0.8, "G": 0.8, "H": 0.8, "I": 1.4, "K": 0.5, "L": 1.6, "M": 1.2, "N": 0.8, "P": 0.4, "Q": 0.8, "R": 0.5, "S": 0.8, "T": 0.8, "V": 1.7, "W": 0.8, "Y": 0.8},
]


def _build_entries() -> list[PrecomputedEntry]:
    """Construct all precomputed entries for H-2Kb."""
    entries: list[PrecomputedEntry] = []

    # --- Binders (strong + moderate + weak) ---
    entries.append(PrecomputedEntry(
        peptide="HITYYVHV",
        binding_score=0.793424,
        ic50_nm=92.35,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "V"},
        anchor_scores={4: 1.7, 7: 1.7},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="SIINFEKL",
        binding_score=0.791356,
        ic50_nm=93.46,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "L"},
        anchor_scores={4: 1.8, 7: 1.6},
        source="known_epitope",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="EAKIYFEV",
        binding_score=0.749992,
        ic50_nm=118.58,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "V"},
        anchor_scores={4: 1.7, 7: 1.7},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="CILDYEPI",
        binding_score=0.683027,
        ic50_nm=174.35,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "I"},
        anchor_scores={4: 1.7, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="VLEFYEEI",
        binding_score=0.643725,
        ic50_nm=218.62,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "I"},
        anchor_scores={4: 1.7, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="NQDMYMEI",
        binding_score=0.643725,
        ic50_nm=218.62,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "I"},
        anchor_scores={4: 1.7, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="GTELFYHV",
        binding_score=0.637995,
        ic50_nm=225.95,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "V"},
        anchor_scores={4: 1.8, 7: 1.7},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="VHLPFRTV",
        binding_score=0.637995,
        ic50_nm=225.95,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "V"},
        anchor_scores={4: 1.8, 7: 1.7},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRLQFSFV",
        binding_score=0.637995,
        ic50_nm=225.95,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "V"},
        anchor_scores={4: 1.8, 7: 1.7},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="YDRQFHPV",
        binding_score=0.637995,
        ic50_nm=225.95,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "V"},
        anchor_scores={4: 1.8, 7: 1.7},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="YKCRYEPV",
        binding_score=0.608887,
        ic50_nm=267.17,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "V"},
        anchor_scores={4: 1.7, 7: 1.7},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="DQQMYRQV",
        binding_score=0.608887,
        ic50_nm=267.17,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "V"},
        anchor_scores={4: 1.7, 7: 1.7},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="RMDAFYWL",
        binding_score=0.607151,
        ic50_nm=269.85,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "L"},
        anchor_scores={4: 1.8, 7: 1.6},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="VFHQFGYL",
        binding_score=0.607151,
        ic50_nm=269.85,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "L"},
        anchor_scores={4: 1.8, 7: 1.6},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="IAHPYKDI",
        binding_score=0.602971,
        ic50_nm=276.42,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "I"},
        anchor_scores={4: 1.7, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="VAPWYADI",
        binding_score=0.602971,
        ic50_nm=276.42,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "I"},
        anchor_scores={4: 1.7, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="CECPYEIL",
        binding_score=0.578975,
        ic50_nm=317.37,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "L"},
        anchor_scores={4: 1.7, 7: 1.6},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="IKGEYCSL",
        binding_score=0.578975,
        ic50_nm=317.37,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "L"},
        anchor_scores={4: 1.7, 7: 1.6},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="MQYTYPWL",
        binding_score=0.578975,
        ic50_nm=317.37,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "L"},
        anchor_scores={4: 1.7, 7: 1.6},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WSISYHEM",
        binding_score=0.566772,
        ic50_nm=340.46,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "M"},
        anchor_scores={4: 1.7, 7: 1.2},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="GQMKFHMI",
        binding_score=0.542672,
        ic50_nm=391.13,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "I"},
        anchor_scores={4: 1.8, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WWTMFHMI",
        binding_score=0.542672,
        ic50_nm=391.13,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "I"},
        anchor_scores={4: 1.8, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="FKGRFQVI",
        binding_score=0.542672,
        ic50_nm=391.13,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "I"},
        anchor_scores={4: 1.8, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="QWYMFRRI",
        binding_score=0.542672,
        ic50_nm=391.13,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "I"},
        anchor_scores={4: 1.8, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="SRKIFTSI",
        binding_score=0.542672,
        ic50_nm=391.13,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "I"},
        anchor_scores={4: 1.8, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WPSAFLPI",
        binding_score=0.542672,
        ic50_nm=391.13,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "I"},
        anchor_scores={4: 1.8, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDSAFCNI",
        binding_score=0.542672,
        ic50_nm=391.13,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "I"},
        anchor_scores={4: 1.8, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WVIWFADC",
        binding_score=0.490864,
        ic50_nm=527.04,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "C"},
        anchor_scores={4: 1.8, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="YMSAFNLM",
        binding_score=0.473941,
        ic50_nm=580.96,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "M"},
        anchor_scores={4: 1.8, 7: 1.2},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="DYDIFEWM",
        binding_score=0.473941,
        ic50_nm=580.96,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "M"},
        anchor_scores={4: 1.8, 7: 1.2},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="TYMDFLIM",
        binding_score=0.473941,
        ic50_nm=580.96,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "M"},
        anchor_scores={4: 1.8, 7: 1.2},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="EQIGFRCM",
        binding_score=0.473941,
        ic50_nm=580.96,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "M"},
        anchor_scores={4: 1.8, 7: 1.2},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRIDFVEC",
        binding_score=0.415363,
        ic50_nm=813.94,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "C"},
        anchor_scores={4: 1.8, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WVFQFCLN",
        binding_score=0.415363,
        ic50_nm=813.94,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "N"},
        anchor_scores={4: 1.8, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="EVKMYLQG",
        binding_score=0.393318,
        ic50_nm=924.07,
        binding_class="weak_binder",
        anchor_residues={4: "Y", 7: "G"},
        anchor_scores={4: 1.7, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="QWPNFFTS",
        binding_score=0.320272,
        ic50_nm=1407.09,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "S"},
        anchor_scores={4: 1.8, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="IGQSFHPC",
        binding_score=0.320272,
        ic50_nm=1407.09,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "C"},
        anchor_scores={4: 1.8, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="MCCWFTTG",
        binding_score=0.320272,
        ic50_nm=1407.09,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "G"},
        anchor_scores={4: 1.8, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="IKMMFDFF",
        binding_score=0.320272,
        ic50_nm=1407.09,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "F"},
        anchor_scores={4: 1.8, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="HKDYYVIG",
        binding_score=0.301565,
        ic50_nm=1567.07,
        binding_class="weak_binder",
        anchor_residues={4: "Y", 7: "G"},
        anchor_scores={4: 1.7, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="NDINLGRV",
        binding_score=0.301565,
        ic50_nm=1567.07,
        binding_class="weak_binder",
        anchor_residues={4: "L", 7: "V"},
        anchor_scores={4: 0.8, 7: 1.7},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="DVFFYVGK",
        binding_score=0.239322,
        ic50_nm=2242.3,
        binding_class="weak_binder",
        anchor_residues={4: "Y", 7: "K"},
        anchor_scores={4: 1.7, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="VIMCFCMP",
        binding_score=0.21453,
        ic50_nm=2586.27,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "P"},
        anchor_scores={4: 1.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WMHSYRFK",
        binding_score=0.172769,
        ic50_nm=3289.11,
        binding_class="weak_binder",
        anchor_residues={4: "Y", 7: "K"},
        anchor_scores={4: 1.7, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="TVYMRYTL",
        binding_score=0.167903,
        ic50_nm=3382.53,
        binding_class="weak_binder",
        anchor_residues={4: "R", 7: "L"},
        anchor_scores={4: 0.4, 7: 1.6},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="FGQYDPPV",
        binding_score=0.126003,
        ic50_nm=4305.19,
        binding_class="weak_binder",
        anchor_residues={4: "D", 7: "V"},
        anchor_scores={4: 0.4, 7: 1.7},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="AEKGYKCE",
        binding_score=0.126003,
        ic50_nm=4305.19,
        binding_class="weak_binder",
        anchor_residues={4: "Y", 7: "E"},
        anchor_scores={4: 1.7, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="HKTSKCDI",
        binding_score=0.123762,
        ic50_nm=4361.1,
        binding_class="weak_binder",
        anchor_residues={4: "K", 7: "I"},
        anchor_scores={4: 0.4, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="DTIKSHVF",
        binding_score=0.114767,
        ic50_nm=4592.86,
        binding_class="weak_binder",
        anchor_residues={4: "S", 7: "F"},
        anchor_scores={4: 0.8, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="IGGQAGMQ",
        binding_score=0.114767,
        ic50_nm=4592.86,
        binding_class="weak_binder",
        anchor_residues={4: "A", 7: "Q"},
        anchor_scores={4: 0.8, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))

    # --- Non-binders ---
    entries.append(PrecomputedEntry(
        peptide="LWPPKDWI",
        binding_score=0.092129,
        ic50_nm=5232.13,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "I"},
        anchor_scores={4: 0.4, 7: 1.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="LIYNIHYK",
        binding_score=0.092129,
        ic50_nm=5232.13,
        binding_class="non_binder",
        anchor_residues={4: "I", 7: "K"},
        anchor_scores={4: 0.8, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="AVRQCHTK",
        binding_score=0.08078,
        ic50_nm=5585.34,
        binding_class="non_binder",
        anchor_residues={4: "C", 7: "K"},
        anchor_scores={4: 0.8, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="AFYRCFDK",
        binding_score=0.069466,
        ic50_nm=5961.22,
        binding_class="non_binder",
        anchor_residues={4: "C", 7: "K"},
        anchor_scores={4: 0.8, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="EISADVMY",
        binding_score=0.06048,
        ic50_nm=6277.71,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "Y"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="VCWVRMEQ",
        binding_score=0.051597,
        ic50_nm=6607.07,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "Q"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="NVNEKWIQ",
        binding_score=0.051597,
        ic50_nm=6607.07,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "Q"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="TVNDDCQA",
        binding_score=0.051597,
        ic50_nm=6607.07,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "A"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="MEVGHHSK",
        binding_score=0.047211,
        ic50_nm=6775.98,
        binding_class="non_binder",
        anchor_residues={4: "H", 7: "K"},
        anchor_scores={4: 0.8, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="SDSNQMMK",
        binding_score=0.047211,
        ic50_nm=6775.98,
        binding_class="non_binder",
        anchor_residues={4: "Q", 7: "K"},
        anchor_scores={4: 0.8, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="VFPRNVQK",
        binding_score=0.047211,
        ic50_nm=6775.98,
        binding_class="non_binder",
        anchor_residues={4: "N", 7: "K"},
        anchor_scores={4: 0.8, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="TLRACSPK",
        binding_score=0.047211,
        ic50_nm=6775.98,
        binding_class="non_binder",
        anchor_residues={4: "C", 7: "K"},
        anchor_scores={4: 0.8, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="MEDMLLRK",
        binding_score=0.047211,
        ic50_nm=6775.98,
        binding_class="non_binder",
        anchor_residues={4: "L", 7: "K"},
        anchor_scores={4: 0.8, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="SEQNRFQG",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "G"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="TYLEDLIN",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "N"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="SRMGCKSD",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "C", 7: "D"},
        anchor_scores={4: 0.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="RYSLCILE",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "C", 7: "E"},
        anchor_scores={4: 0.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="KWHQEVIF",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "F"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="KTVSRDYC",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "C"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="CHKVELRS",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "S"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="TFDKQMTP",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "Q", 7: "P"},
        anchor_scores={4: 0.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="PPFSDFTW",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "W"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="DMVPMSVP",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "M", 7: "P"},
        anchor_scores={4: 0.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="RLVGRRNW",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "W"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="MQGHKVNT",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "T"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="ESYQECNH",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "H"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="IWMAGKCD",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "G", 7: "D"},
        anchor_scores={4: 0.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="SRNTDERT",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "T"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="PYTLDTYQ",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "Q"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="YYVAKKAG",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "G"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="DTHPKRMG",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "G"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="STLCKPYC",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "C"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="MFTDETNA",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "A"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="TCMYEYPF",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "F"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="HNWRELPT",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "T"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="YFKCSNVD",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "S", 7: "D"},
        anchor_scores={4: 0.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="LDGKETFQ",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "Q"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="MFCNVMGE",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "V", 7: "E"},
        anchor_scores={4: 0.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="FDRCLHCP",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "L", 7: "P"},
        anchor_scores={4: 0.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="NTKDQDQP",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "Q", 7: "P"},
        anchor_scores={4: 0.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="NPMKREAG",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "G"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKKREYLS",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "S"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WDFMERRT",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "T"},
        anchor_scores={4: 0.4, 7: 0.8},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WCADCWKD",
        binding_score=0.026261,
        ic50_nm=7644.47,
        binding_class="non_binder",
        anchor_residues={4: "C", 7: "D"},
        anchor_scores={4: 0.8, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="MIMGKTEE",
        binding_score=0.019368,
        ic50_nm=7953.89,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "E"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="YITPRRLR",
        binding_score=0.016831,
        ic50_nm=8070.93,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WQWPKLER",
        binding_score=0.01255,
        ic50_nm=8272.28,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="SVLLDIPK",
        binding_score=0.01255,
        ic50_nm=8272.28,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="CYEHDKDK",
        binding_score=0.008665,
        ic50_nm=8459.36,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="KSLMERDK",
        binding_score=0.008665,
        ic50_nm=8459.36,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="PQMYRMDR",
        binding_score=0.008665,
        ic50_nm=8459.36,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WNRVRYDK",
        binding_score=0.008665,
        ic50_nm=8459.36,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="LVAVDDIE",
        binding_score=0.003563,
        ic50_nm=8711.48,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "E"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="YEWCKVQK",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="KYVSEQWK",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="SRADKIPR",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="FKCCDLFK",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="DCGLKWLR",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="KWCGRTRK",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="LPVCRDMK",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="FSLTEQSR",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="RETRKFQK",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="SYRQDHTR",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="NLDDDYTR",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="HMYSKFCR",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="FMDWDNLK",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="CQNIEYLK",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="TNDTRTTK",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="MKENEPRK",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "K"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WCYFRNNR",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="GNTIEHFR",
        binding_score=0.002565,
        ic50_nm=8761.69,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "R"},
        anchor_scores={4: 0.4, 7: 0.5},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="KNFSKDDD",
        binding_score=0.001704,
        ic50_nm=8805.24,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="DKRIRWYP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="GSHNKMKD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="WWSFRVSP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="IDQEERGD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="LNNQEITD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="YPYIKWFD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="RLHPKEIP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="FECLRCWP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="VNDPKKVP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="YAKFDVKP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="VHQRDQMP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="CLYLREWD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="FKLMRASE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "E"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="RNDQDIQP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="GSRWKRWD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="SRLKRCND",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="DTGCERRP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="NAQCDTND",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="RYADRKHP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="ILECRQYD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="TWIFRQAD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="NKANDNIE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "E"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="NRCNRQKE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "E"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="PFNTKEME",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 7: "E"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="TRIHRNCD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 7: "D"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="EQFKEPNP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="DPTRDKYP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))
    entries.append(PrecomputedEntry(
        peptide="VNFHDTPP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 7: "P"},
        anchor_scores={4: 0.4, 7: 0.4},
        source="pssm_predicted",
        peptide_length=8,
    ))

    return entries


# Lazy singleton — built on first access
_database: PrecomputedAlleleDatabase | None = None


def get_database() -> PrecomputedAlleleDatabase:
    """Return the precomputed binding database for H-2Kb."""
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