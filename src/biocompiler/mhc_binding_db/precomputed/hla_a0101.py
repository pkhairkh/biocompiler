"""Precomputed MHC binding data for HLA-A*01:01.

Anchor residues: position 2 (T/S), position 9 (Y)
Known epitopes: YLDVSSNYI, ATLGFFSQY, SLYSGFFYI
Data source: IEDB and SYFPEITHI databases + PSSM predictions

Binding scores are derived from PSSM-based prediction with expected
AUC-ROC of 0.60-0.75. Do not replace experimental validation or
NetMHCpan predictions where available.
"""

from __future__ import annotations

from .. import PrecomputedAlleleDatabase, PrecomputedEntry

ALLELE: str = "HLA-A*01:01"
PEPTIDE_LENGTH: int = 9
ANCHOR_DESCRIPTION: str = "position 2 (T/S), position 9 (Y)"
KNOWN_EPITOPES: list[str] = ['YLDVSSNYI', 'ATLGFFSQY', 'SLYSGFFYI']


def _build_entries() -> list[PrecomputedEntry]:
    """Construct all precomputed entries for HLA-A*01:01."""
    entries: list[PrecomputedEntry] = []

    # --- Binders (strong + moderate) ---
    entries.append(PrecomputedEntry(
        peptide="ATLGFFSQY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="known_epitope",
    ))
    entries.append(PrecomputedEntry(
        peptide="STSWCEMNY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATWKIWKKY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATHPVADYY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STKLWPIWY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATAVKEHTY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STDDKVHSY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATSVDAEAY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STTAVINFY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATFKLDTMY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATIWREMHY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STGASEPRY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATMVIVSNY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATFYSRKLY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STYGCFLKY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STPWIWGTY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STYDRLVWY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STYGIVMFY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STDYGGVIY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STKRDWYAY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STQRRLQFY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STVFINNDY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STSRPGWPY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STKIAHLYY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STYATVCPY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATWVIDIMY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STKWGYGDY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STWFWRMEY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATVEKEFMY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATHHTQVQY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STWWRTMCY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STEKKYTRY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATPKKMAIY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STGGMPRTY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STMPPVSFY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STFPENVIY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATSWHAQEY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STFRQDWYY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STKNYMSPY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATGWMEQGY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATEPPVKQY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STFKPKYSY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STTCIAGLY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STYDCKNEY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STEETHQCY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATQQSWYVY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STNGCGINY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STWQNGQIY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATHKRQDFY",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'T', 8: 'Y'},
        anchor_scores={1: 1.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGHHYWYDY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'G', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SCGLKSWAY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'C', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATRINDFRQ",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'Q'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STIRLPLIT",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'T'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATWINPVLT",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'T'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATHNYAASL",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'L'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SNLADSGPY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'N', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AKGGVFKKY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'K', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SNPTQHLFY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'N', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AYLWIWYKY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'Y', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SGVPMKSYY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'G', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AYDTLHLKY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'Y', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AWCDPSLLY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SYRVKGNMY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'Y', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AYQIRSLCY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'Y', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SRWPPLCLY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AKHPWNFDY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'K', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APIYGTKYY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANEKTIQSY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'N', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATVDADGWA",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'A'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AQKWPFEWY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'Q', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGNMCFGVY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'G', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SRCNPQFHY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SAKFPKSCY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AHWRTTPIY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'H', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SNHRWRMTY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'N', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACSTTIHYY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'C', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APKRMHFNY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SWCNFWPVY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SNAQKDCTY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'N', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AHMKKAPAY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'H', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SPNPHFQLY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STRQHELEH",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'H'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATIFCKWVG",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'G'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SYCIYVMIY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'Y', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SRDATPPWY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STFCVYLNP",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'P'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ARQNNVCIY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SMQGQDGEY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'M', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ARNAKYEQY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SHCECLTRY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'H', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AQQVLLVDY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'Q', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STSQTKIYT",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'T'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SATKSCRWY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STHLYCTPC",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'C'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SRSRLLPGY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STQVCVTTV",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'V'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AWEFQEEIY",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'Y'},
        anchor_scores={1: 0.8, 8: 1.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATRKWMGEC",
        binding_score=0.749645,
        ic50_nm=50.3,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'C'},
        anchor_scores={1: 1.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SLYSGFFYI",
        binding_score=0.55,
        ic50_nm=400.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'L', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.8},
        source="known_epitope",
    ))
    entries.append(PrecomputedEntry(
        peptide="YLDVSSNYI",
        binding_score=0.55,
        ic50_nm=400.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'L', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.8},
        source="known_epitope",
    ))

    # --- Non-binders (weak + non) ---
    entries.append(PrecomputedEntry(
        peptide="YMKYGHTVG",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'G'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NNKMQEYLT",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'N', 8: 'T'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QYYFMERES",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'S'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HMSNNFVLM",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'M'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MGTAKLTIH",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'H'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VNDAYCLRQ",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'N', 8: 'Q'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PMMETQMKP",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'P'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CCEWLMTPH",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'C', 8: 'H'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WSENTCGWW",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'S', 8: 'W'},
        anchor_scores={1: 1.6, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LYEQECLLS",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'S'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PYETEQSVM",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'M'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NAKTPCEIA",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'A'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ESRKWVTIR",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'S', 8: 'R'},
        anchor_scores={1: 1.6, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IAPAFEPWV",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'V'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CHLNKQGQI",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'H', 8: 'I'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QRRQWGQRT",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'R', 8: 'T'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KPGMNRSVS",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'P', 8: 'S'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGWKPMCSA",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'A'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GCPISSVRQ",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'C', 8: 'Q'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NYKEHSSHA",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'A'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KQIQSTGWG",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Q', 8: 'G'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YYPYERPWI",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'I'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LYKDGSLHG",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'G'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NKHMHMFSN",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'K', 8: 'N'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GNMLWDPEV",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'N', 8: 'V'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LKYGNPSQM",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'K', 8: 'M'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QKQHLHTYL",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GNWRGEDVM",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'N', 8: 'M'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GWFQLHMSL",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'W', 8: 'L'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KCRCRKFIW",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'C', 8: 'W'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VHKLRGVEC",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'H', 8: 'C'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VHGSQTMMG",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'H', 8: 'G'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TYVSPHLCL",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'L'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FHATHDSQA",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'H', 8: 'A'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="INYVFPLKC",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'N', 8: 'C'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ECHVAYTIS",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'C', 8: 'S'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YSQFTRNDK",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'S', 8: 'K'},
        anchor_scores={1: 1.6, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FRIQQMIPP",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'R', 8: 'P'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TASMYHRCT",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'T'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YGAKSMDRW",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'W'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CPFCINSCI",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'P', 8: 'I'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QCVCRMFLC",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'C', 8: 'C'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RSMVCHMVS",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'S', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LWHSHTKNT",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'W', 8: 'T'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MNFWRIYHL",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'N', 8: 'L'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PQWQVNYGV",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Q', 8: 'V'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DPDYEGDKC",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'P', 8: 'C'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FHVNMEHDQ",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'H', 8: 'Q'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ICENQMSKT",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'C', 8: 'T'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WSTQALRRI",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'S', 8: 'I'},
        anchor_scores={1: 1.6, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IGQAQQQYW",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'W'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QHCHYQFRL",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'H', 8: 'L'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YYILRDVVC",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'C'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGISRDATN",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'N'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CAFRVFNCN",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'N'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NYQEGYDGS",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'S'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EYPYDYVYS",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'S'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VYPNMIPPI",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'I'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GNDFNHWRQ",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'N', 8: 'Q'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DCRIHSNWC",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'C', 8: 'C'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LYLTSVIKN",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'N'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WKGTAPGIF",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'K', 8: 'F'},
        anchor_scores={1: 0.8, 8: 1.6},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IAEGNAVWM",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'M'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EYYHMWLKN",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'N'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CPSIHYPKG",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'P', 8: 'G'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TWKNMGDMN",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'W', 8: 'N'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NAEMCDNPH",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'H'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FGYEYAAWI",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MYRFWYDFM",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Y', 8: 'M'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VNKPGHEEQ",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'N', 8: 'Q'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CPQIDEHLI",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'P', 8: 'I'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QRRTISGIL",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NNDQLIKPC",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'N', 8: 'C'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KQYALKEAT",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'Q', 8: 'T'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WSMFRLFEN",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'S', 8: 'N'},
        anchor_scores={1: 1.6, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HCFPAERSA",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'C', 8: 'A'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WSIVMQYCN",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'S', 8: 'N'},
        anchor_scores={1: 1.6, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DAHTINCIS",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'S'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LHTKYDCIP",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'H', 8: 'P'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CWECNFWWL",
        binding_score=0.495302,
        ic50_nm=523.4,
        binding_class="weak_binder",
        anchor_residues={1: 'W', 8: 'L'},
        anchor_scores={1: 0.8, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVGWCQITK",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVDYQGTEK",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIAICTTME",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIIWGKWTE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVTFAPGVR",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVYKKNWKR",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ALEDAEDWE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIYPYPSID",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SLFTGGAFD",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SFLCAWQKE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVYTFDYWK",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SLACGYKHD",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SFFTPRWLD",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVPQFTKPK",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ALALIDNME",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ALEIGDKCE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVKQEYWFR",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIMMPRQNE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SFSWHRDIE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AFYGLSSID",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SLVETQGVD",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AFTSTYIND",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SLELFSEFE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SFQKCMWWE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVKKNTGAR",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ALMDPLMCD",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AFYQDIQDE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVTEFNAGK",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ALWGHTQGE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SFWSTMYYE",
        binding_score=0.226689,
        ic50_nm=6212.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LIDYDNSWD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RHALKGCQE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'E'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RRCDGVLTD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'R', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WGNGGVYMD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WHYQLKVKE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'E'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RVALQYQSN",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'N'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KFRVLMYNE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VVQQKSTRK",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TFCHMRQDD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WVRIWGAQN",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'N'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WYSTHAQWD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'Y', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TVKTEAAEK",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FLTMPLGCD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WKLDEMLKE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'K', 8: 'E'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WVTILVGGL",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'L'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PVELLTRDK",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WGCPKAKWE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'E'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GISPAGADD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IFPFNMLRE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LVRSNYEKR",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WGQMVLEMD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HITIPGISE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RVEDGDWSP",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'P'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WCFFQSDWD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RAQHWWYGD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GLVNVGYSE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RVNANTRNV",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'V'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KIDIIKMVD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LFFEVMDGD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WRTSGVHCD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'R', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QIIMAKMND",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FFPKVSMHE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GVNKHKAGR",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WVGICVHKW",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'W'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WVYPYGWTH",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'H'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WWGDNNHRD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'W', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NIMNTDDPD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GLGHNDANE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DVTQYDWTR",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CVTFSQYCK",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QFLHRQFNE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GFAVAVEWE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WVRYRMCIM",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'M'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RATHKYIWD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EFPHIFECE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RCENYVEFD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ILYHCICWE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ILKDTMYHE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YFEDMWQSD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EIEIPWMAD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PLRLFSAKE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MVKKRQLHR",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WVVDFAFDV",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'V'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LLRANHKDD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PVYVPKSDK",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WNRRFIMVE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'E'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WVAPCGGNI",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'I'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PLLLDSWKD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WRGPGTIVE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'R', 8: 'E'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RYPSSIFLE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'Y', 8: 'E'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NINKTFVSE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GFFWSWRSE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EIMYDRTQE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YFGTWQTLE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WVASICHHW",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'W'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RQKNHKQFD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'Q', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NVQHCLPGR",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VVEFGVHPK",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LFSQWKRFE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RMKRRYMPD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NLCVMNWLE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GLYCQWGVD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ILCETYCDD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HIMFRHKLD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MVDRGNCQK",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FVKTQWKLK",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RHRHEIYAE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'E'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CVENGDRIR",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KIQLPWRFE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RVPAESYRI",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'I'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RCMKHWKQE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'E'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WYGDQPEWD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'Y', 8: 'D'},
        anchor_scores={1: 0.8, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WVMVCGRRG",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'G'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MFACYEKYD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FFNMGIDNE",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HVDHHEQER",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.5, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RVLKDNAMP",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'P'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HFTVCKSYD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CLACEKCKD",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.4, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WVEVDNAYM",
        binding_score=0.204228,
        ic50_nm=7639.9,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'M'},
        anchor_scores={1: 0.5, 8: 0.8},
        source="pssm_predicted",
    ))

    return entries


# Lazy singleton — built on first access
_database: PrecomputedAlleleDatabase | None = None


def get_database() -> PrecomputedAlleleDatabase:
    """Return the precomputed binding database for HLA-A*01:01."""
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
