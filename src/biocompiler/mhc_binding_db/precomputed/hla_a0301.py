"""Precomputed MHC binding data for HLA-A*03:01.

Anchor residues: position 2 (V/L/I), position 9 (K/R)
Known epitopes: RLIYLCRLV, RLRAEAQVK, KLVWAGIGL
Data source: IEDB and SYFPEITHI databases + PSSM predictions

Binding scores are derived from PSSM-based prediction with expected
AUC-ROC of 0.60-0.75. Do not replace experimental validation or
NetMHCpan predictions where available.
"""

from __future__ import annotations

from .. import PrecomputedAlleleDatabase, PrecomputedEntry

ALLELE: str = "HLA-A*03:01"
PEPTIDE_LENGTH: int = 9
ANCHOR_DESCRIPTION: str = "position 2 (V/L/I), position 9 (K/R)"
KNOWN_EPITOPES: list[str] = ['RLIYLCRLV', 'RLRAEAQVK', 'KLVWAGIGL']


def _build_entries() -> list[PrecomputedEntry]:
    """Construct all precomputed entries for HLA-A*03:01."""
    entries: list[PrecomputedEntry] = []

    # --- Binders (strong + moderate) ---
    entries.append(PrecomputedEntry(
        peptide="SIGDWFSYK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SINARDQFK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIWQHDHQK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SITTDSWPK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SITMKWCCK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVAVFPKWK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVDHMWQQK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIWVIDEGK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIQHNWKIK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVWHFITFK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVDMFGYAK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVYEVLMDK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIHFAFAPK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVTCQMSFK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIGGIKDLK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIQCIITFK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIWSFARIK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVLPMCFVK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIHWATEYK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIGKEGTPK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIVAMSVHK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVDMWVCMK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVQVFCVCK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVCNNVHYK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVRSIDVAK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVEVLETQK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVSNPPGMK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AITLTKFPK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVQEWKFHK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVFDSSQFK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVDTLVTGK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVCNQFTLK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIKIAHLYK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIWDDQIKK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIFTVDFQK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AISCVAMTK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SICAIYDCK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIHGYCECK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIQDLYVFK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVTKIHWEK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVFNARLIK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIVEWGANK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVVAQWMTK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SVMAASKVK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'V', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIWLELSNK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SIGMMNPHK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIIDYAVMK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIYPGVSCK",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 1.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RLRAEAQVK",
        binding_score=0.927316,
        ic50_nm=9.8,
        binding_class="strong_binder",
        anchor_residues={1: 'L', 8: 'K'},
        anchor_scores={1: 1.6, 8: 2.0},
        source="known_epitope",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATTKEHPEK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SWKKCHHIK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AWCDPSLLK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AAPAYSGFK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ARYCCHCGK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STPNMCHLK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACITMLRWK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'C', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ASWGKHVGK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'S', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SSPFTGNKK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'S', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SFGVREDQK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'F', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SHFPKHNIK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'H', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STSYMCDPK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SYYIQVYEK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'Y', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SAYAALYEK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACMPKHIFK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'C', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APKWTTGKK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AKSDLGSFK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'K', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AYHYDQESK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'Y', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ARMINYKTK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SKPVEQFAK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'K', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SWQMNFYFK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AAQYPGQTK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SAWTPRWMK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ARHNNGSTK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SGWRLNYIK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGEHYMYAK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATECWDHLK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SAFWREDQK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APSNMAMRK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGVPGWTEK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ARWSLRSKK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SWNSNTDEK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AWSLGKYPK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SRGPFYFTK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SPYTPEEMK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SKQSKLIIK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'K', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AAPIINNRK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ARRENANLK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'R', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SKQIPIGLK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'K', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AKPDFFKSK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'K', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SSDFEAGCK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'S', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SWYTNTLAK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SCAANFCNK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'C', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACITGYNNK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'C', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AYDAAIQQK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'Y', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AAYNTLCIK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ASGHSECVK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'S', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="STLPDEHTK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SCERTVEIK",
        binding_score=0.732352,
        ic50_nm=59.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'C', 8: 'K'},
        anchor_scores={1: 0.8, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RLIYLCRLV",
        binding_score=0.593314,
        ic50_nm=212.2,
        binding_class="moderate_binder",
        anchor_residues={1: 'L', 8: 'V'},
        anchor_scores={1: 1.6, 8: 0.7},
        source="known_epitope",
    ))
    entries.append(PrecomputedEntry(
        peptide="KLVWAGIGL",
        binding_score=0.593314,
        ic50_nm=212.2,
        binding_class="moderate_binder",
        anchor_residues={1: 'L', 8: 'L'},
        anchor_scores={1: 1.6, 8: 0.7},
        source="known_epitope",
    ))

    # --- Non-binders (weak + non) ---
    entries.append(PrecomputedEntry(
        peptide="VDACEGWVK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MLRCWTGIT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EDTHDVMDK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QLTLCSCVS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NDFGNHSQK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FLEFHRYFT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TDATPCLKK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KELCGGPYK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MDCIDQRAK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TLSVQCDNT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DLLNSADMT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GLEVIHYMS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DEALDKMDK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HMFVPTRLT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RMHVGEQES",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FMALDMRKS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DMPLITCTS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DDQHCVFHK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CDHHCHLFK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TLFQKIHNS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QLEKPRMRT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YDNVDQAQK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NMQHDLTFS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MEAVAVSGK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RMYHTHQVT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EDSQFDMTK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RDDATQFIK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EDMQVGNVK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IEKDIRNWK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HDFLCEIGK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FMGNYDSRS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LDEVMNPNK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LEDSCITLK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PDLCHMDWK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RMYYYQNFT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QDNAQFNMK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TMVKPWHYS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WDIWLNCCK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DMAFDSSVT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MMNGINDCT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YLWSNGNKS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PENTPMGKK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CDRLNLQYK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PLFPGLYVS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FLIFWIIIT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CLSYWGCRS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QEKQLGTKK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HLSPDEKRT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FDDMPTQTK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PDHYSGGGK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NMPTCLLDT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VDTHKRCAK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HMNFCWRYT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LLDRRWAES",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ILNHETIRT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GDDIPLQAK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IMMKQNTWT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KDGSMKIPK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HLTCTTFYS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DEMFYWGPK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VEFDWQRIK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NLEKDCTYT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TEEYWVTYK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EDRFKYNVK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NETICRYIK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FEPPRINEK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IESSQIIDK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IESTFIEHK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LMWWFIEIS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GLAQIVRLS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DLRQGPAET",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NDQRMACLK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VEIGVSSNK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YDDASGACK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LMDTSGCHS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'M', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GEYWLRRIK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'E', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NLEPEHPNS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IDAYCPWFK",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'D', 8: 'K'},
        anchor_scores={1: 0.4, 8: 2.0},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WLKKQMFLT",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'T'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TLWDKEHFS",
        binding_score=0.49425,
        ic50_nm=528.4,
        binding_class="weak_binder",
        anchor_residues={1: 'L', 8: 'S'},
        anchor_scores={1: 1.6, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEWVWIDGM",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADSKVVGMA",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEMWYLPKG",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'G'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADENINVKC",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'C'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEVNMHWTQ",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEYSHQERW",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SESKINKWG",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'G'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEWGMCKFP",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'P'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADNATMFNA",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEGGYYWYQ",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AERNTQMQV",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'V'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADYPRPWGN",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'N'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEYMQLVDV",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'V'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDIMAMPLY",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Y'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADDQFHATV",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'V'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEQDMNQWP",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'P'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADQWPYSPL",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'L'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SERPIFDHV",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'V'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDSYNKWIF",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'F'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AETREYREQ",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADEHIKCCW",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SETADSLLL",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'L'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADADPIDSM",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEWAYSDNL",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'L'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDYDWTNFN",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'N'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDLMDNIYI",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEWRHDENY",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Y'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADCLHWYHI",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEFAMHSMQ",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDTMRVDIL",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'L'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADLPEGGKQ",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDLRQMSHP",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'P'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADNLCHSPL",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'L'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADTISACMQ",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADTWNEMAY",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Y'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEAWGFLYP",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'P'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDRCWWDIA",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEEGWTWEI",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEYSWTPWQ",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADRDSDIPN",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'N'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADRDPCHWF",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'F'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDHGQFRIA",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEDGKMSYC",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'C'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDPGGRGAY",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Y'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEFDVIWVM",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADGSYSLLM",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEWKIDWMF",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'F'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDPMQKFPL",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'L'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDRRQPGHQ",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEFCQPLKP",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'P'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADSYRYPKW",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AENHFNRAM",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADEAGFKNM",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SERLMCQWC",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'C'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEEQPNVYM",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADKHRPSIQ",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADGEQENFI",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADTQVWTRC",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'C'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SESWCTWFG",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'G'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEWAYLINI",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SESVDEIYP",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'P'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDDFQPPHI",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDMDTDHEW",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AELWTDDCF",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'F'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEWSEGMEF",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'F'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SECWPCMGW",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADNWHAFYG",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'G'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDDSTYGNM",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDVRGIKCL",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'L'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEIGGICFV",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'V'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADQRKIYHW",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDQEPELTI",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEQYGDYFI",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AEDVQKAPG",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'G'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADMLWNRFY",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Y'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADRWQLVLN",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'N'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADECEYKMN",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'N'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AETDAINPV",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'V'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDNFFTTVQ",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDLNSPHPC",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'C'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AENMHVFCW",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SELLAGSSI",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEGAHGWVW",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADFQTPAEY",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Y'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDEIAKNFL",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'L'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDEYKAPMG",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'G'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEYRDTIRM",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADHVWGRVA",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDTGFSGPL",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'L'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEHNTEQNY",
        binding_score=0.232567,
        ic50_nm=5884.8,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Y'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PDWGIQYEL",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'L'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YESGGLACW",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IDTPRWLNW",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GEQTWRITV",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'V'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KEWWVRGSA",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IEHFRTWEQ",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WDCHPTMII",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NDEMKQGDA",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WDMHCNSHC",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'C'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MEKHHNQKW",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'W'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FEEAHHCHI",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="REVIVFAIQ",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Q'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LEVQWYRRM",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QEHVYMCRA",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HEPHVSGGV",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'V'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IEVGINSDA",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TDVQWHGHF",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'F'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VEHPIDARN",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'N'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NEQEHKKSA",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PDDNCSIII",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'I'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VDHLPMRCG",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'G'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PDLGMNNLA",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'A'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HENQIMYPM",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VDDEAHDIC",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'C'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VEGIEIMDP",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'P'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HDSHGTSKG",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'G'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WDREFTQCP",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'P'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MDYGDFDTP",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'P'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HDHHSKRDC",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'C'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TEGTFKENP",
        binding_score=0.207929,
        ic50_nm=7383.9,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'P'},
        anchor_scores={1: 0.4, 8: 0.7},
        source="pssm_predicted",
    ))

    return entries


# Lazy singleton — built on first access
_database: PrecomputedAlleleDatabase | None = None


def get_database() -> PrecomputedAlleleDatabase:
    """Return the precomputed binding database for HLA-A*03:01."""
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
