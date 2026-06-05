"""Precomputed MHC binding data for H-2Db.

MHC-I presents 9-mer peptides.

Anchor residues: position 5 (N/M), position 9 (M/V/I/L)
Known epitopes: ASNENMETM, NMAVMDFAV, LTVQVARVK
Data source: IEDB and SYFPEITHI databases + PSSM predictions

Binding scores are derived from PSSM-based prediction with expected
AUC-ROC of 0.60-0.75. Do not replace experimental validation or
NetMHCpan predictions where available.

References
----------
- Falk et al., Nature 1991; 351:290 (H-2Db motif)
- Rammensee et al., Immunogenetics 1995 (SYFPEITHI)
"""
from __future__ import annotations

from .. import PrecomputedAlleleDatabase, PrecomputedEntry

ALLELE: str = "H-2Db"
MHC_CLASS: str = "I"
PEPTIDE_LENGTH: int = 9
ANCHOR_DESCRIPTION: str = "position 5 (N/M), position 9 (M/V/I/L)"
KNOWN_EPITOPES: list[str] = ['ASNENMETM', 'NMAVMDFAV', 'LTVQVARVK']

# ═══════════════════════════════════════════════════════════════════════════
# Position-Specific Scoring Matrix (PSSM) — core 9-mer
# ═══════════════════════════════════════════════════════════════════════════

PSSM: list[dict[str, float]] = [
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.3, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.2, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 0.8, "C": 0.8, "D": 0.5, "E": 0.5, "F": 0.8, "G": 0.8, "H": 0.8, "I": 0.8, "K": 0.5, "L": 0.8, "M": 1.6, "N": 1.7, "P": 0.8, "Q": 0.8, "R": 0.5, "S": 0.8, "T": 0.8, "V": 0.8, "W": 0.8, "Y": 0.8},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 0.8, "C": 0.8, "D": 0.4, "E": 0.4, "F": 0.8, "G": 0.8, "H": 0.8, "I": 1.5, "K": 0.5, "L": 1.4, "M": 1.8, "N": 0.8, "P": 0.4, "Q": 0.8, "R": 0.5, "S": 0.8, "T": 0.8, "V": 1.6, "W": 0.8, "Y": 0.8},
]


def _build_entries() -> list[PrecomputedEntry]:
    """Construct all precomputed entries for H-2Db."""
    entries: list[PrecomputedEntry] = []

    # --- Binders (strong + moderate + weak) ---
    entries.append(PrecomputedEntry(
        peptide="ASNENMETM",
        binding_score=1.0,
        ic50_nm=28.12,
        binding_class="strong_binder",
        anchor_residues={4: "N", 8: "M"},
        anchor_scores={4: 1.7, 8: 1.8},
        source="known_epitope",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KQSKNLMAM",
        binding_score=0.938338,
        ic50_nm=40.1,
        binding_class="strong_binder",
        anchor_residues={4: "N", 8: "M"},
        anchor_scores={4: 1.7, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MNNHNPHLI",
        binding_score=0.863151,
        ic50_nm=61.82,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "I"},
        anchor_scores={4: 1.7, 8: 1.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YVTLNLGGM",
        binding_score=0.807062,
        ic50_nm=85.38,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "M"},
        anchor_scores={4: 1.7, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GTADNWQIM",
        binding_score=0.807062,
        ic50_nm=85.38,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "M"},
        anchor_scores={4: 1.7, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LCYDMKFHM",
        binding_score=0.76615,
        ic50_nm=108.05,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "M"},
        anchor_scores={4: 1.6, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LHQPMFWRM",
        binding_score=0.76615,
        ic50_nm=108.05,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "M"},
        anchor_scores={4: 1.6, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="QIRYMGTTM",
        binding_score=0.76615,
        ic50_nm=108.05,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "M"},
        anchor_scores={4: 1.6, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AYWNNKRRV",
        binding_score=0.728795,
        ic50_nm=133.97,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "V"},
        anchor_scores={4: 1.7, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NMAVMDFAV",
        binding_score=0.690442,
        ic50_nm=167.07,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "V"},
        anchor_scores={4: 1.6, 8: 1.6},
        source="known_epitope",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YWEVMIDAV",
        binding_score=0.690442,
        ic50_nm=167.07,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "V"},
        anchor_scores={4: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VLCLMTVYV",
        binding_score=0.690442,
        ic50_nm=167.07,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "V"},
        anchor_scores={4: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TKDIMHWRV",
        binding_score=0.690442,
        ic50_nm=167.07,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "V"},
        anchor_scores={4: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NFFPMDWGV",
        binding_score=0.690442,
        ic50_nm=167.07,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "V"},
        anchor_scores={4: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RMVNMPWGV",
        binding_score=0.690442,
        ic50_nm=167.07,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "V"},
        anchor_scores={4: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FMQDMINMV",
        binding_score=0.690442,
        ic50_nm=167.07,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "V"},
        anchor_scores={4: 1.6, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NMRKNGQWI",
        binding_score=0.68801,
        ic50_nm=169.42,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "I"},
        anchor_scores={4: 1.7, 8: 1.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LKLNNSTAI",
        binding_score=0.68801,
        ic50_nm=169.42,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "I"},
        anchor_scores={4: 1.7, 8: 1.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SEKVNIAMI",
        binding_score=0.68801,
        ic50_nm=169.42,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "I"},
        anchor_scores={4: 1.7, 8: 1.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TPWCMRPWI",
        binding_score=0.651024,
        ic50_nm=209.62,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "I"},
        anchor_scores={4: 1.6, 8: 1.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VWTAMVVWI",
        binding_score=0.651024,
        ic50_nm=209.62,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "I"},
        anchor_scores={4: 1.6, 8: 1.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HARCMNAEI",
        binding_score=0.651024,
        ic50_nm=209.62,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "I"},
        anchor_scores={4: 1.6, 8: 1.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TVYKNCSNL",
        binding_score=0.646018,
        ic50_nm=215.75,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "L"},
        anchor_scores={4: 1.7, 8: 1.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDVFNMIAL",
        binding_score=0.646018,
        ic50_nm=215.75,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "L"},
        anchor_scores={4: 1.7, 8: 1.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YWMPNKPIL",
        binding_score=0.646018,
        ic50_nm=215.75,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "L"},
        anchor_scores={4: 1.7, 8: 1.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="QMIANIPEL",
        binding_score=0.646018,
        ic50_nm=215.75,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "L"},
        anchor_scores={4: 1.7, 8: 1.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IEKMNAVFL",
        binding_score=0.646018,
        ic50_nm=215.75,
        binding_class="moderate_binder",
        anchor_residues={4: "N", 8: "L"},
        anchor_scores={4: 1.7, 8: 1.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="CGMKMMEQL",
        binding_score=0.610463,
        ic50_nm=264.75,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "L"},
        anchor_scores={4: 1.6, 8: 1.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DWTQMHTML",
        binding_score=0.610463,
        ic50_nm=264.75,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "L"},
        anchor_scores={4: 1.6, 8: 1.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IWGWMSFIL",
        binding_score=0.610463,
        ic50_nm=264.75,
        binding_class="moderate_binder",
        anchor_residues={4: "M", 8: "L"},
        anchor_scores={4: 1.6, 8: 1.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TRNFNCMLW",
        binding_score=0.483305,
        ic50_nm=550.47,
        binding_class="weak_binder",
        anchor_residues={4: "N", 8: "W"},
        anchor_scores={4: 1.7, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NFSVMIEGY",
        binding_score=0.416114,
        ic50_nm=810.43,
        binding_class="weak_binder",
        anchor_residues={4: "M", 8: "Y"},
        anchor_scores={4: 1.6, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LKDFAVGLM",
        binding_score=0.387293,
        ic50_nm=956.68,
        binding_class="weak_binder",
        anchor_residues={4: "A", 8: "M"},
        anchor_scores={4: 0.8, 8: 1.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SAFANKSEG",
        binding_score=0.362795,
        ic50_nm=1101.57,
        binding_class="weak_binder",
        anchor_residues={4: "N", 8: "G"},
        anchor_scores={4: 1.7, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="PTVANWCRW",
        binding_score=0.362795,
        ic50_nm=1101.57,
        binding_class="weak_binder",
        anchor_residues={4: "N", 8: "W"},
        anchor_scores={4: 1.7, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GEPPNLEGQ",
        binding_score=0.362795,
        ic50_nm=1101.57,
        binding_class="weak_binder",
        anchor_residues={4: "N", 8: "Q"},
        anchor_scores={4: 1.7, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VSCSNMTTS",
        binding_score=0.362795,
        ic50_nm=1101.57,
        binding_class="weak_binder",
        anchor_residues={4: "N", 8: "S"},
        anchor_scores={4: 1.7, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NRILMTSKH",
        binding_score=0.337847,
        ic50_nm=1271.69,
        binding_class="weak_binder",
        anchor_residues={4: "M", 8: "H"},
        anchor_scores={4: 1.6, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RIVFHVCMV",
        binding_score=0.337847,
        ic50_nm=1271.69,
        binding_class="weak_binder",
        anchor_residues={4: "H", 8: "V"},
        anchor_scores={4: 0.8, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KSVFSSRRV",
        binding_score=0.337847,
        ic50_nm=1271.69,
        binding_class="weak_binder",
        anchor_residues={4: "S", 8: "V"},
        anchor_scores={4: 0.8, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HPGKYANCL",
        binding_score=0.286564,
        ic50_nm=1708.41,
        binding_class="weak_binder",
        anchor_residues={4: "Y", 8: "L"},
        anchor_scores={4: 0.8, 8: 1.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="CYQVNMRTK",
        binding_score=0.195864,
        ic50_nm=2879.65,
        binding_class="weak_binder",
        anchor_residues={4: "N", 8: "K"},
        anchor_scores={4: 1.7, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DFPYMDRPK",
        binding_score=0.178556,
        ic50_nm=3181.33,
        binding_class="weak_binder",
        anchor_residues={4: "M", 8: "K"},
        anchor_scores={4: 1.6, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HMCYEIWRV",
        binding_score=0.178556,
        ic50_nm=3181.33,
        binding_class="weak_binder",
        anchor_residues={4: "E", 8: "V"},
        anchor_scores={4: 0.5, 8: 1.6},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IICCDYNWL",
        binding_score=0.143616,
        ic50_nm=3890.1,
        binding_class="weak_binder",
        anchor_residues={4: "D", 8: "L"},
        anchor_scores={4: 0.5, 8: 1.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FRKSNYNKE",
        binding_score=0.136595,
        ic50_nm=4050.53,
        binding_class="weak_binder",
        anchor_residues={4: "N", 8: "E"},
        anchor_scores={4: 1.7, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRPYMPFSE",
        binding_score=0.12255,
        ic50_nm=4391.61,
        binding_class="weak_binder",
        anchor_residues={4: "M", 8: "E"},
        anchor_scores={4: 1.6, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LHQRSVEYG",
        binding_score=0.12255,
        ic50_nm=4391.61,
        binding_class="weak_binder",
        anchor_residues={4: "S", 8: "G"},
        anchor_scores={4: 0.8, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LTWQQIPKA",
        binding_score=0.12255,
        ic50_nm=4391.61,
        binding_class="weak_binder",
        anchor_residues={4: "Q", 8: "A"},
        anchor_scores={4: 0.8, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YVKYPHCSA",
        binding_score=0.12255,
        ic50_nm=4391.61,
        binding_class="weak_binder",
        anchor_residues={4: "P", 8: "A"},
        anchor_scores={4: 0.8, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))

    # --- Non-binders ---
    entries.append(PrecomputedEntry(
        peptide="PYSERLPGY",
        binding_score=0.067202,
        ic50_nm=6039.43,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "Y"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LTVQVARVK",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "V", 8: "K"},
        anchor_scores={4: 0.8, 8: 0.5},
        source="known_epitope",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKKFDPERC",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "C"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HDWTECKNF",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "F"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ATQCKHIDH",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "H"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FRIMKNPVG",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "G"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TQYIGAQCR",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "G", 8: "R"},
        anchor_scores={4: 0.8, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SKCSEWQAT",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "T"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HLDSEYWIG",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "G"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FEKIFSTSR",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "F", 8: "R"},
        anchor_scores={4: 0.8, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TFIIKFPWF",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "F"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MIVYDMSLW",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "W"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YQAFDEKGQ",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "Q"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HSMCKLEMT",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "T"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TSWRRFGFN",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "N"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HYDYDMMKA",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "A"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TQWHESGAT",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "T"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="VQKHATRNR",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "A", 8: "R"},
        anchor_scores={4: 0.8, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TCRDELKVY",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "Y"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IIEWEETVG",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "G"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ITYKEVAHH",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "H"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="PIQDWDSTK",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "W", 8: "K"},
        anchor_scores={4: 0.8, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IVHGDRSKH",
        binding_score=0.041269,
        ic50_nm=7011.77,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "H"},
        anchor_scores={4: 0.5, 8: 0.8},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NTAQWNLFE",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "W", 8: "E"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HSRMFYCTD",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "F", 8: "D"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DHIHLHQCE",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "L", 8: "E"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HVLQLEDED",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "L", 8: "D"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YVCGYMCRD",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "Y", 8: "D"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GYFFAYDLP",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "A", 8: "P"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LRCLPERIE",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "P", 8: "E"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SYWRQQYKE",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "Q", 8: "E"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YKWPFWDSD",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "F", 8: "D"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DVCHGLPRP",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "G", 8: "P"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HPPLVCSFE",
        binding_score=0.018505,
        ic50_nm=7993.5,
        binding_class="non_binder",
        anchor_residues={4: "V", 8: "E"},
        anchor_scores={4: 0.8, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LFSHKMWKR",
        binding_score=0.013673,
        ic50_nm=8218.98,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HQNSEIQMP",
        binding_score=0.005634,
        ic50_nm=8608.26,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EDNRDYWRE",
        binding_score=0.005634,
        ic50_nm=8608.26,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IKRTEWDSK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="IQDDEYYAR",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="PEKEKCGLK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KAPHRIQER",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HWLEENYAR",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DGDHDAVCR",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MVPHEWWSR",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="WHLCESWAK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="FGKKRNLTK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERMKEPTTK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KHKDEPKFR",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ANLEKMRLK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YGCLEHIPR",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RQHIKMHAK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DIRKKNGNK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GEFARMKEK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LIMMDCMHR",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KEFPDCMGK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RSQGDVWRK",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "K"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DNIVKHPQR",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="QEREEAATR",
        binding_score=0.004058,
        ic50_nm=8686.72,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "R"},
        anchor_scores={4: 0.5, 8: 0.5},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YESRDVPPP",
        binding_score=0.002696,
        ic50_nm=8755.05,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ATSDDHGRE",
        binding_score=0.002696,
        ic50_nm=8755.05,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="PGSKKVDKE",
        binding_score=0.002696,
        ic50_nm=8755.05,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DCIAEIRTP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="LAIKDGGHE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="QEKWECHRP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TYTEEERMP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHVYDARCE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="DDISREETP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EHCSDNNHP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YDGAEKTFE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GMQIKHKKD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "D"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NRQNKWFQP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TECHRKNMP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="PCVHKYHTE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AAYHDFHFP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ERLHDSDCP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="AFQCKWDKD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "D"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NTARDCTCE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KGGEDMVFP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MNHNDWALE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="PHQSKTCCP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="HSTTKMGDD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "D"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ETGYELFQD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "D"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GWQDDHLME",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="PDICKGHCE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="WFTPRHFDP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="RQKSESCVE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="YAANDTKFP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="NNCKRTNVE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KCVQKNWNE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="WVFFRKVIE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="GFDSKVWFD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "D"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MRKLKIEGD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "D"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="KQHGKSKPP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="MKGMRLPKE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="REFLKRCTE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "K", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="SSYFESIVP",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "E", 8: "P"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="TGWHRFLPE",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "E"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="ASQHRKGYD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "D"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EQIHRHNID",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "R", 8: "D"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))
    entries.append(PrecomputedEntry(
        peptide="EKTCDRNPD",
        binding_score=0.0,
        ic50_nm=8892.01,
        binding_class="non_binder",
        anchor_residues={4: "D", 8: "D"},
        anchor_scores={4: 0.5, 8: 0.4},
        source="pssm_predicted",
        peptide_length=9,
    ))

    return entries


# Lazy singleton — built on first access
_database: PrecomputedAlleleDatabase | None = None


def get_database() -> PrecomputedAlleleDatabase:
    """Return the precomputed binding database for H-2Db."""
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