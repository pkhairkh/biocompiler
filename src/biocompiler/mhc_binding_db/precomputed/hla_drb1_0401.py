"""Precomputed MHC binding data for HLA-DRB1*04:01.

MHC-II presents 13-25 aa peptides with a 9-mer core binding region.
Peptides are 15-mers scored by scanning all possible
9-mer core registers within each peptide.

Anchor residues: position 1 (F/Y/W), position 4 (D/E), position 6 (A/S/G), position 9 (L/I/V)
Known epitopes: GQYRAEMFDSDV, PVSKMRMATPLLMQA
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

ALLELE: str = "HLA-DRB1*04:01"
MHC_CLASS: str = "II"
PEPTIDE_LENGTH: int = 15
ANCHOR_DESCRIPTION: str = "position 1 (F/Y/W), position 4 (D/E), position 6 (A/S/G), position 9 (L/I/V)"
KNOWN_EPITOPES: list[str] = ['GQYRAEMFDSDV', 'PVSKMRMATPLLMQA']

# ═══════════════════════════════════════════════════════════════════════════
# Position-Specific Scoring Matrix (PSSM) — core 9-mer
# ═══════════════════════════════════════════════════════════════════════════

PSSM: list[dict[str, float]] = [
    {"A": 0.9, "C": 0.9, "D": 0.4, "E": 0.4, "F": 1.8, "G": 0.9, "H": 0.9, "I": 0.9, "K": 0.5, "L": 1.3, "M": 0.9, "N": 0.9, "P": 0.9, "Q": 0.9, "R": 0.5, "S": 0.9, "T": 0.9, "V": 0.9, "W": 1.6, "Y": 1.7},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 0.8, "C": 0.8, "D": 1.8, "E": 1.6, "F": 0.8, "G": 0.8, "H": 0.8, "I": 0.8, "K": 0.4, "L": 0.8, "M": 0.8, "N": 0.8, "P": 0.8, "Q": 0.8, "R": 0.4, "S": 0.8, "T": 0.8, "V": 0.8, "W": 0.5, "Y": 0.8},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.6, "C": 0.9, "D": 0.9, "E": 0.9, "F": 0.6, "G": 1.3, "H": 0.9, "I": 0.9, "K": 0.9, "L": 0.9, "M": 0.9, "N": 1.2, "P": 0.9, "Q": 0.9, "R": 0.9, "S": 1.4, "T": 0.9, "V": 0.9, "W": 0.5, "Y": 0.6},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 1.0, "C": 1.0, "D": 1.0, "E": 1.0, "F": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "K": 1.0, "L": 1.0, "M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0, "R": 1.0, "S": 1.0, "T": 1.0, "V": 1.0, "W": 1.0, "Y": 1.0},
    {"A": 0.9, "C": 0.9, "D": 0.5, "E": 0.5, "F": 1.3, "G": 0.9, "H": 0.9, "I": 1.3, "K": 0.5, "L": 1.4, "M": 0.9, "N": 0.9, "P": 0.9, "Q": 0.9, "R": 0.9, "S": 0.9, "T": 0.9, "V": 1.3, "W": 0.9, "Y": 0.9},
]


def _build_entries() -> list[PrecomputedEntry]:
    """Construct all precomputed entries for HLA-DRB1*04:01."""
    entries: list[PrecomputedEntry] = []

    # --- Binders (strong + moderate + weak) ---
    entries.append(PrecomputedEntry(
        peptide="YVSDTAQFIFKVTNI",
        binding_score=0.93511,
        ic50_nm=40.85,
        binding_class="strong_binder",
        anchor_residues={0: "Y", 3: "D", 5: "A", 8: "I"},
        anchor_scores={0: 1.7, 3: 1.8, 5: 1.6, 8: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DHGNFSLDDSYMFED",
        binding_score=0.898776,
        ic50_nm=50.36,
        binding_class="moderate_binder",
        anchor_residues={4: "F", 7: "D", 9: "S", 12: "F"},
        anchor_scores={4: 1.8, 7: 1.8, 9: 1.4, 12: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="SFFAVDYGGHIIGFA",
        binding_score=0.864492,
        ic50_nm=61.34,
        binding_class="moderate_binder",
        anchor_residues={2: "F", 5: "D", 7: "G", 10: "I"},
        anchor_scores={2: 1.8, 5: 1.8, 7: 1.3, 10: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="ICWFFMDKGEHVFIG",
        binding_score=0.864492,
        ic50_nm=61.34,
        binding_class="moderate_binder",
        anchor_residues={3: "F", 6: "D", 8: "G", 11: "V"},
        anchor_scores={3: 1.8, 6: 1.8, 8: 1.3, 11: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="LGIDMANGIMRSAWW",
        binding_score=0.811929,
        ic50_nm=83.02,
        binding_class="moderate_binder",
        anchor_residues={0: "L", 3: "D", 5: "A", 8: "I"},
        anchor_scores={0: 1.3, 3: 1.8, 5: 1.6, 8: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="VCFKIEMNYTLTFHE",
        binding_score=0.809342,
        ic50_nm=84.26,
        binding_class="moderate_binder",
        anchor_residues={2: "F", 5: "E", 7: "N", 10: "L"},
        anchor_scores={2: 1.8, 5: 1.6, 7: 1.2, 10: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HSLILIYEVASLLQS",
        binding_score=0.793026,
        ic50_nm=92.56,
        binding_class="moderate_binder",
        anchor_residues={4: "L", 7: "E", 9: "A", 12: "L"},
        anchor_scores={4: 1.3, 7: 1.6, 9: 1.6, 12: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GGSYKWEKGSGVSPS",
        binding_score=0.787254,
        ic50_nm=95.69,
        binding_class="moderate_binder",
        anchor_residues={3: "Y", 6: "E", 8: "G", 11: "V"},
        anchor_scores={3: 1.7, 6: 1.6, 8: 1.3, 11: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="MMTLMADCSHSIKIV",
        binding_score=0.755116,
        ic50_nm=115.14,
        binding_class="moderate_binder",
        anchor_residues={3: "L", 6: "D", 8: "S", 11: "I"},
        anchor_scores={3: 1.3, 6: 1.8, 8: 1.4, 11: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="ESCYSWELNLYVFTH",
        binding_score=0.753595,
        ic50_nm=116.15,
        binding_class="moderate_binder",
        anchor_residues={3: "Y", 6: "E", 8: "N", 11: "V"},
        anchor_scores={3: 1.7, 6: 1.6, 8: 1.2, 11: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="CVDGHNWIFDLQYDL",
        binding_score=0.690755,
        ic50_nm=166.77,
        binding_class="moderate_binder",
        anchor_residues={6: "W", 9: "D", 11: "Q", 14: "L"},
        anchor_scores={6: 1.6, 9: 1.8, 11: 0.9, 14: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GYENDPQNEIPYVPL",
        binding_score=0.685515,
        ic50_nm=171.87,
        binding_class="moderate_binder",
        anchor_residues={1: "Y", 4: "D", 6: "Q", 9: "I"},
        anchor_scores={1: 1.7, 4: 1.8, 6: 0.9, 9: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="AYLTYEVGDLFHCKS",
        binding_score=0.67844,
        ic50_nm=179.02,
        binding_class="moderate_binder",
        anchor_residues={2: "L", 5: "E", 7: "G", 10: "F"},
        anchor_scores={2: 1.3, 5: 1.6, 7: 1.3, 10: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WTHCLEYEHGFTIYP",
        binding_score=0.67844,
        ic50_nm=179.02,
        binding_class="moderate_binder",
        anchor_residues={4: "L", 7: "E", 9: "G", 12: "I"},
        anchor_scores={4: 1.3, 7: 1.6, 9: 1.3, 12: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="SQWMRYMDEEMHTLM",
        binding_score=0.668734,
        ic50_nm=189.31,
        binding_class="moderate_binder",
        anchor_residues={5: "Y", 8: "E", 10: "M", 13: "L"},
        anchor_scores={5: 1.7, 8: 1.6, 10: 0.9, 13: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="QTVEGPTLPDPANDI",
        binding_score=0.6623,
        ic50_nm=196.45,
        binding_class="moderate_binder",
        anchor_residues={6: "T", 9: "D", 11: "A", 14: "I"},
        anchor_scores={6: 0.9, 9: 1.8, 11: 1.6, 14: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="LSGPPKYWMMWAYTL",
        binding_score=0.624904,
        ic50_nm=243.64,
        binding_class="moderate_binder",
        anchor_residues={6: "Y", 9: "M", 11: "A", 14: "L"},
        anchor_scores={6: 1.7, 9: 0.8, 11: 1.6, 14: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="SELRSDGSSIQQMQI",
        binding_score=0.613098,
        ic50_nm=260.77,
        binding_class="moderate_binder",
        anchor_residues={2: "L", 5: "D", 7: "S", 10: "Q"},
        anchor_scores={2: 1.3, 5: 1.8, 7: 1.4, 10: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GWGVEGNKPNCHEHA",
        binding_score=0.590335,
        ic50_nm=297.28,
        binding_class="moderate_binder",
        anchor_residues={1: "W", 4: "E", 6: "N", 9: "N"},
        anchor_scores={1: 1.6, 4: 1.6, 6: 1.2, 9: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HTYAVWNNEDNCCPV",
        binding_score=0.590335,
        ic50_nm=297.28,
        binding_class="moderate_binder",
        anchor_residues={5: "W", 8: "E", 10: "N", 13: "P"},
        anchor_scores={5: 1.6, 8: 1.6, 10: 1.2, 13: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="ALYADSGFNMQMAVK",
        binding_score=0.586912,
        ic50_nm=303.19,
        binding_class="moderate_binder",
        anchor_residues={1: "L", 4: "D", 6: "G", 9: "M"},
        anchor_scores={1: 1.3, 4: 1.8, 6: 1.3, 9: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="YNLNSSTKFCQQRMS",
        binding_score=0.552531,
        ic50_nm=369.55,
        binding_class="moderate_binder",
        anchor_residues={0: "Y", 3: "N", 5: "S", 8: "F"},
        anchor_scores={0: 1.7, 3: 0.8, 5: 1.4, 8: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="CWIYQYETDFADCDE",
        binding_score=0.548849,
        ic50_nm=377.46,
        binding_class="moderate_binder",
        anchor_residues={5: "Y", 8: "D", 10: "A", 13: "D"},
        anchor_scores={5: 1.7, 8: 1.8, 10: 1.6, 13: 0.5},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DLLAVECVMCVLTRY",
        binding_score=0.546891,
        ic50_nm=381.75,
        binding_class="moderate_binder",
        anchor_residues={2: "L", 5: "E", 7: "V", 10: "V"},
        anchor_scores={2: 1.3, 5: 1.6, 7: 0.9, 10: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EGLWWWGDNYWNLII",
        binding_score=0.544925,
        ic50_nm=386.09,
        binding_class="moderate_binder",
        anchor_residues={4: "W", 7: "D", 9: "Y", 12: "L"},
        anchor_scores={4: 1.6, 7: 1.8, 9: 0.6, 12: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="DLQMHFAPCLHMRYM",
        binding_score=0.532538,
        ic50_nm=414.62,
        binding_class="moderate_binder",
        anchor_residues={1: "L", 4: "H", 6: "A", 9: "L"},
        anchor_scores={1: 1.3, 4: 0.8, 6: 1.6, 9: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KDCPHWAKDYAMCKT",
        binding_score=0.528948,
        ic50_nm=423.28,
        binding_class="moderate_binder",
        anchor_residues={5: "W", 8: "D", 10: "A", 13: "K"},
        anchor_scores={5: 1.6, 8: 1.8, 10: 1.6, 13: 0.5},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHEYYTPLCGCVVER",
        binding_score=0.528159,
        ic50_nm=425.21,
        binding_class="moderate_binder",
        anchor_residues={4: "Y", 7: "L", 9: "G", 12: "V"},
        anchor_scores={4: 1.7, 7: 0.8, 9: 1.3, 12: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="AIHEINLQDSSDKTH",
        binding_score=0.490672,
        ic50_nm=527.62,
        binding_class="weak_binder",
        anchor_residues={5: "N", 8: "D", 10: "S", 13: "T"},
        anchor_scores={5: 0.9, 8: 1.8, 10: 1.4, 13: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="FLNLKWSASMNHKIC",
        binding_score=0.483929,
        ic50_nm=548.5,
        binding_class="weak_binder",
        anchor_residues={5: "W", 8: "S", 10: "N", 13: "I"},
        anchor_scores={5: 1.6, 8: 0.8, 10: 1.2, 13: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EAMYDLRCVFRQGDL",
        binding_score=0.468202,
        ic50_nm=600.48,
        binding_class="weak_binder",
        anchor_residues={1: "A", 4: "D", 6: "R", 9: "F"},
        anchor_scores={1: 0.9, 4: 1.8, 6: 0.9, 9: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="TQCPELLIGMYSKGF",
        binding_score=0.467792,
        ic50_nm=601.89,
        binding_class="weak_binder",
        anchor_residues={6: "L", 9: "M", 11: "S", 14: "F"},
        anchor_scores={6: 1.3, 9: 0.8, 11: 1.4, 14: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="CDVNTFCKEDGMPKS",
        binding_score=0.464504,
        ic50_nm=613.39,
        binding_class="weak_binder",
        anchor_residues={5: "F", 8: "E", 10: "G", 13: "K"},
        anchor_scores={5: 1.8, 8: 1.6, 10: 1.3, 13: 0.5},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WKSFCNTETKEELRM",
        binding_score=0.45529,
        ic50_nm=646.81,
        binding_class="weak_binder",
        anchor_residues={4: "C", 7: "E", 9: "K", 12: "L"},
        anchor_scores={4: 0.9, 7: 1.6, 9: 0.9, 12: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="CIRIVDASAPVMSEV",
        binding_score=0.425993,
        ic50_nm=765.63,
        binding_class="weak_binder",
        anchor_residues={2: "R", 5: "D", 7: "S", 10: "V"},
        anchor_scores={2: 0.5, 5: 1.8, 7: 1.4, 10: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PIHHKWNECSSYLRV",
        binding_score=0.421666,
        ic50_nm=784.94,
        binding_class="weak_binder",
        anchor_residues={5: "W", 8: "C", 10: "S", 13: "R"},
        anchor_scores={5: 1.6, 8: 0.8, 10: 1.4, 13: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WNDHMLCALRWKNIP",
        binding_score=0.421666,
        ic50_nm=784.94,
        binding_class="weak_binder",
        anchor_residues={0: "W", 3: "H", 5: "L", 8: "L"},
        anchor_scores={0: 1.6, 3: 0.8, 5: 0.9, 8: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NNGPVLWVRSLQCAL",
        binding_score=0.421666,
        ic50_nm=784.94,
        binding_class="weak_binder",
        anchor_residues={6: "W", 9: "S", 11: "Q", 14: "L"},
        anchor_scores={6: 1.6, 9: 0.8, 11: 0.9, 14: 1.4},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="EFDTRQAMANEAMHV",
        binding_score=0.401385,
        ic50_nm=882.14,
        binding_class="weak_binder",
        anchor_residues={6: "A", 9: "N", 11: "A", 14: "V"},
        anchor_scores={6: 0.9, 9: 0.8, 11: 1.6, 14: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="TPIGRAVLVPWHLPT",
        binding_score=0.401385,
        ic50_nm=882.14,
        binding_class="weak_binder",
        anchor_residues={0: "T", 3: "G", 5: "A", 8: "V"},
        anchor_scores={0: 0.9, 3: 0.8, 5: 1.6, 8: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NHETTHWRDLVMYCT",
        binding_score=0.366848,
        ic50_nm=1076.17,
        binding_class="weak_binder",
        anchor_residues={5: "H", 8: "D", 10: "V", 13: "C"},
        anchor_scores={5: 0.9, 8: 1.8, 10: 0.9, 13: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="IKSMPTNDDCQLQMG",
        binding_score=0.366848,
        ic50_nm=1076.17,
        binding_class="weak_binder",
        anchor_residues={4: "P", 7: "D", 9: "C", 12: "Q"},
        anchor_scores={4: 0.9, 7: 1.8, 9: 0.9, 12: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="ELCNPLCMNSRWVFC",
        binding_score=0.3665,
        ic50_nm=1078.32,
        binding_class="weak_binder",
        anchor_residues={4: "P", 7: "M", 9: "S", 12: "V"},
        anchor_scores={4: 0.9, 7: 0.8, 9: 1.4, 12: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="YCKLNMDIYCQCNWH",
        binding_score=0.35736,
        ic50_nm=1136.58,
        binding_class="weak_binder",
        anchor_residues={3: "L", 6: "D", 8: "Y", 11: "C"},
        anchor_scores={3: 1.3, 6: 1.8, 8: 0.6, 11: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="CTALLNWTQCQCVVY",
        binding_score=0.348035,
        ic50_nm=1199.26,
        binding_class="weak_binder",
        anchor_residues={4: "L", 7: "T", 9: "C", 12: "V"},
        anchor_scores={4: 1.3, 7: 0.8, 9: 0.9, 12: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="IQGSFTKWSGANTSC",
        binding_score=0.313792,
        ic50_nm=1460.56,
        binding_class="weak_binder",
        anchor_residues={4: "F", 7: "W", 9: "G", 12: "T"},
        anchor_scores={4: 1.8, 7: 0.5, 9: 1.3, 12: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RPREIHMNVTDLAQH",
        binding_score=0.287583,
        ic50_nm=1698.41,
        binding_class="weak_binder",
        anchor_residues={0: "R", 3: "E", 5: "H", 8: "V"},
        anchor_scores={0: 0.5, 3: 1.6, 5: 0.9, 8: 1.3},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="VLARAQLWQQCPFHE",
        binding_score=0.265378,
        ic50_nm=1929.99,
        binding_class="weak_binder",
        anchor_residues={1: "L", 4: "A", 6: "L", 9: "Q"},
        anchor_scores={1: 1.3, 4: 0.8, 6: 0.9, 9: 0.9},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="GQYRAEMFDSDV",
        binding_score=0.258892,
        ic50_nm=2003.41,
        binding_class="weak_binder",
        anchor_residues={2: "Y", 5: "E", 7: "F", 10: "D"},
        anchor_scores={2: 1.7, 5: 1.6, 7: 0.6, 10: 0.5},
        source="known_epitope",
        peptide_length=12,
    ))
    entries.append(PrecomputedEntry(
        peptide="PVSKMRMATPLLMQA",
        binding_score=0.256508,
        ic50_nm=2031.1,
        binding_class="weak_binder",
        anchor_residues={2: "S", 5: "R", 7: "A", 10: "L"},
        anchor_scores={2: 0.9, 5: 0.4, 7: 1.6, 10: 1.4},
        source="known_epitope",
        peptide_length=15,
    ))

    # --- Non-binders ---
    entries.append(PrecomputedEntry(
        peptide="KKHWRRRRHRPKWKP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHHKRRKHPRWKKHR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KHPWKKRPWKPRKWH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HRPRPRKKWKWKKKP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HPHRRKWRKRWWPWH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HKRRWRKKWWHWKPW",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHWRKKKWHKWPHPK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPKRKWRWRWKPPKP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PRPRRKWKWRPKWHK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRPWRRRKWKWHWHK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHKRRRKRKWKHKPK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PRRKPRWKRRWWKWP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRHRRKKKKKRWHWK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PKPRRKKKHWWHRPW",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HRKKWKRWPWWPKRR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKHKKKPKWKWRWPP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHRRKHRKKRKPPHR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HKRKHKKRKKPRWKR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KHHKRRRHRWKPRKW",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PPPKRRKRRKRKRPH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PPPRRKRKWPWWHHK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHRRRKRHWWKPHPR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KHRRRKPHWWPWKWH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PPRRHWWRRKWRPRK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHHRKKPKPRWKWWK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PKRRPKRKKKRHWPW",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PWPKKRWRRRPPPHK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRWHWRRRKRKWKHR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WRPRKKKRKRPKHRK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KHPKKRPKWKRPPWW",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHPRRKRPWWPRWRP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HPKKRRRWRWWKWWR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HRHRRKPWWWKHKRK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RPRKRHKKRRKHKPH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKKKRPKKKKKKPRH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHPRKRWWRRPWKWK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PRRRKRRPRRHPKKW",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRKPRWKWRPWPWRK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PHHKKRRWKRHWHHP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HRKKKRKRRKWRHHH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PRPRRKKRWKWWPKK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RHRRKKPKWRRRPWH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PPHKKKPPWKHPKHP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKKHKRRPKWHPKPK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRHRPKRRHKHRRKK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PKPRHRPRWKKKHRH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HRPRWRKKWWKHKKH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WPHKRRRKKWHPRKP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HRKRPKWKPRWPRPK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPHWRRRRKHHKHHK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PRHRRRKWWRPRHWR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKKHKWKWKKHWHKK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PPWKKRKRKRKPKRK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KHRWRKKRWHPWKPR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRHPRKRKKRKPRPH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KHPWKKKKRKWKHHR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HPPRKPPWKWKKKWK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRKPKHKKKHKPKKK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RHPRKKPWWRWKHKH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KHHPRRRPKWWWHPW",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PHPKRRRKKRWKWKK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKKKRHRRRRKHPWH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PHPKKRKRRWHHPPW",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHPRRKWRKKHHRPK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRHWWRKKKWRKPWH",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HRPKRRRRKWHRKHR",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PKRKKRKPHWRKKKK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KPRRRWKPKRRKKKW",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKPKRRKPRWWKRKP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="NNKKKRKNKWPWKRK",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RPPWKKKWWWRWWRP",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RNPRKRRPKKWWKRW",
        binding_score=0.099682,
        ic50_nm=5009.51,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WRKKRRRKKKKPRPH",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WKWRKKRWKPWHWRK",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PRKRWRRKKKRHKHH",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKPRPRWWKKPWKKR",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KWRRKRWHPRWRKWK",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PKWRRKRWKWRKKWK",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WPPRKWRWKKWWRPP",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKRKWRHKWKKWKHR",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRKKWWRRWRWKKKW",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRRWWHRKRWKKWKK",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRRWRRKPWWPWPHK",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKPRKWKKRWKPPKP",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="WPKKKRRRKKKPPHW",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HKHWKWRWRWWKKPR",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HHHWKWRWKKHKRKH",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KNRRRKWKWRRPPPK",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKWKKKRRNWKNRKR",
        binding_score=0.098298,
        ic50_nm=5049.57,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HKRKPKRWKKHRKKK",
        binding_score=0.063721,
        ic50_nm=6161.67,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRWKKKRWRKKWRRH",
        binding_score=0.063721,
        ic50_nm=6161.67,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRRRRRRWRKKKRKR",
        binding_score=0.063721,
        ic50_nm=6161.67,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PKRKKRKWKRKHRPR",
        binding_score=0.063721,
        ic50_nm=6161.67,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="HKHRKKRWKKHPHRR",
        binding_score=0.063721,
        ic50_nm=6161.67,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KRKKRRPHWWWKKWK",
        binding_score=0.063721,
        ic50_nm=6161.67,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RKRRPHWWRKHWKKK",
        binding_score=0.063721,
        ic50_nm=6161.67,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="PRPKKRRWKKHRHKK",
        binding_score=0.063721,
        ic50_nm=6161.67,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RRHKKRRWPKHRKKW",
        binding_score=0.054069,
        ic50_nm=6513.71,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="KKKKKKPRKRKKRPK",
        binding_score=0.054069,
        ic50_nm=6513.71,
        binding_class="non_binder",
        anchor_residues={},
        anchor_scores={},
        source="pssm_predicted",
        peptide_length=15,
    ))
    entries.append(PrecomputedEntry(
        peptide="RNKKKKRRWKNKRKK",
        binding_score=0.0462,
        ic50_nm=6815.52,
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
    """Return the precomputed binding database for HLA-DRB1*04:01."""
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