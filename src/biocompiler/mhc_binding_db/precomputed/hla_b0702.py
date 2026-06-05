"""Precomputed MHC binding data for HLA-B*07:02.

Anchor residues: position 2 (P), position 9 (L)
Known epitopes: RPPIFIRRL, SPRWYFYYL
Data source: IEDB and SYFPEITHI databases + PSSM predictions

Binding scores are derived from PSSM-based prediction with expected
AUC-ROC of 0.60-0.75. Do not replace experimental validation or
NetMHCpan predictions where available.
"""

from __future__ import annotations

from .. import PrecomputedAlleleDatabase, PrecomputedEntry

ALLELE: str = "HLA-B*07:02"
PEPTIDE_LENGTH: int = 9
ANCHOR_DESCRIPTION: str = "position 2 (P), position 9 (L)"
KNOWN_EPITOPES: list[str] = ['RPPIFIRRL', 'SPRWYFYYL']


def _build_entries() -> list[PrecomputedEntry]:
    """Construct all precomputed entries for HLA-B*07:02."""
    entries: list[PrecomputedEntry] = []

    # --- Binders (strong + moderate) ---
    entries.append(PrecomputedEntry(
        peptide="APITLKKLL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APVWEIVFL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPWHIYYIL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPWDLLKNL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPFIMLSFL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APIVRGTVL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPSSRSICL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APQDYAESL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APASHSFPL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APIDQDESL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPGAMYIAL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPHKWEAPL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPPSATSPL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APEHFRTYL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPIMCNETL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APWGVKMFL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APCMQQKEL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APHEHLCFL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPKISRDKL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPNAWTYKL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APDWGHRVL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPHCPAIQL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APTLNSQML",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APQQAMWWL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APACCGFNL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPSDDDMKL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APRCDTHFL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPWAWDHRL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APVDYLADL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APVHILEVL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPGFIMCSL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPHYIDPKL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPGMATTGL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APWEAKNFL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APKRMHFNL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APRRPCEFL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPEYQIYYL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPVFEKDRL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APSPLRSCL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPAMNNNPL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APGCYWLDL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APERAQWCL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPLTWYICL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPGTYRFTL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APHRFDMEL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APACERKTI",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'I'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPRCKTSCL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="APQTLWQSL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PPIDRKVWL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RPPIFIRRL",
        binding_score=0.965804,
        ic50_nm=6.9,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="known_epitope",
    ))
    entries.append(PrecomputedEntry(
        peptide="SPRWYFYYL",
        binding_score=0.965804,
        ic50_nm=6.9,
        binding_class="strong_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 2.0, 8: 1.5},
        source="known_epitope",
    ))
    entries.append(PrecomputedEntry(
        peptide="LPFGCTGMS",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'S'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WPSELPIYS",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'S'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EPHQAWTNT",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'T'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KPFLCWRCM",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'M'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IPSYCESLA",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'A'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WPWAGHAMS",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'S'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VPRWLVHMN",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'N'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NPCSFIDSN",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'N'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NPTGSFMPQ",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'Q'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GPWRSVTDN",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'N'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IPSCVWRGA",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'A'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QPYPPGQSG",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'G'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HPIRDSLMM",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'M'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RPVFHYRTT",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'T'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MPIPNPSRT",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'T'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SPGRAYESQ",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'Q'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CPNYCWWIG",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'G'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SPVHFHWMH",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'H'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DPYQVHFLC",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'C'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EPEHGAMTM",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'M'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QPWVYPLIG",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'G'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YPLLQELNP",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'P'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MPPCAPQEF",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'F'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FPGKKWEDY",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'Y'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LPMRDVTQP",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'P'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IPISMEVRW",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'W'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NPMTTCFIP",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'P'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EPPVAAIPM",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'M'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YPGAFRTDT",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'T'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EPFQRFSVM",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'M'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FPEKPLDGG",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'G'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EPLDSYVKH",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'H'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FPAQHSESP",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'P'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GPAHPCQYF",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'F'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IPELNYDYH",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'H'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TPVKGVIYC",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'C'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EPHGMCTWQ",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'Q'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EPGQDKTEM",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'M'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GPTVPTHEW",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'W'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NPIQKMRCF",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'F'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CPMARHMHM",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'M'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KPRTGCSWM",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'M'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QPTIMEPDF",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'F'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WPEYYITWP",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'P'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TPKSVLDNF",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'F'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EPNHAKYQQ",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'Q'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NPMFNSGKY",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'Y'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EPFCIPFMF",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'F'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FPYNKGVCW",
        binding_score=0.749114,
        ic50_nm=50.5,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'W'},
        anchor_scores={1: 2.0, 8: 0.8},
        source="pssm_predicted",
    ))

    # --- Non-binders (weak + non) ---
    entries.append(PrecomputedEntry(
        peptide="IANMENKID",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MACYYFNGD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FAPVQQLLR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TATMWAQLK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KAKLPNHDK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KASDPHFMR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IAETMWMHK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EADLFACMD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EASPRKGIK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KAHGDRHTK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GARQYMVSR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EAFFYHGFD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SAYFSHEED",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LAYPHCMAE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DAHWSVPRR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TAHTLNLTD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YALTNNGYE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WAAFVSEMK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FAWDYKNRE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RAWWEYWGK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CAWGCYTNK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FACQQYDWE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HAVPPKGRR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HAHMGGTLE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YAYSSIDCK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TANEYWTGE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GAMINRDPD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FAVKICAWR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FAWADLHGR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TANNHWPSR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FACRGCNHK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MAAHLNLME",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FAHHLHYSD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FAEIDITHD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NAIHSVRDD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EARITMSEK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LAEWTPGSR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DAERFNKTR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WANECAMVK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GADKSLRFD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TAKDKDMRK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VAMYHVTCK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SAVGIFPER",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EACPNAAWE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GAGHRWWLE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SAEVGTTSR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HAHFHRILK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RAWEPIKYE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YAVWPSHVD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RANERLACK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KANFHSHKD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VATKSRWWD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TAKQGMKDD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HAGYMLGYR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FARYQCEKD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TAYKPPIHD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YANDDKCHR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KAYYCPTAD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QATYDWRKK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DAYFTKCRK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DACGWVSTR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FAHVGWFMR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RAFNGPRKR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KASRYWYYE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IANVFMWIR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LAFNEHPSD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YAAQPRIKR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HANDWFIPR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QAVRHQWKK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DAKRSVKCR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GALPIYNMR",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'R'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IAYLEYWWE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FAKLKGQEK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RATHVNMTK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TAAVHLCTE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VAYHSPFLE",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'E'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GAFTYARFK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GASKEYHDK",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'K'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DAPTAHVLD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SAVGHYTTD",
        binding_score=0.494808,
        ic50_nm=525.7,
        binding_class="weak_binder",
        anchor_residues={1: 'A', 8: 'D'},
        anchor_scores={1: 1.8, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PMRSEVTVE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANEAMRVHR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PFWPPSQMD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PYGDDGMNE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Y', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PFCQMCNRE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVCCNYGWD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PFGLRWDGK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PVARWPPFR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PIHFADEPD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANHAPDVEK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGMIAFRLD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ASFEFNMND",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PVCLFIYRD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AMNPTYDCR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PIECPMNLD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGPYHICDD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PMNQDPPER",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PMLFSHKFK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AFFEHKWFD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AHYWRTCGD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PHEAWEFTE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PCHFVCLDE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIYPYPSID",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACFAWEWGR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AHCHCTNLK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANGLGRHFD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PVGPTKPND",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANSNMMSIR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PQFANHSWK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Q', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AQNPCSCFR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Q', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIQILRWYD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AQLMRIMND",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Q', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PIDQWAWVK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PNVHQDATE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACPCIGPGD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PCNYHVYIK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANNKPAAED",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PHNRWWHSE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ASRYDRPLK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AMVMDEYNE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PHETQGFGD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AMNRLCHLE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AYNWENMKD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Y', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACTMQWQCR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACDGQKETD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AYCPRHKYR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Y', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PSVIRFRID",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIWIKAVDE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGKSVEQLE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGTHRKPID",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PFMDIMWDE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AQTSCYIGK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Q', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVGCTSCRE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PNQSCTISR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANLETFSDR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AYFLQWCIE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Y', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PYATALFER",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Y', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PMKVNYVRE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PIGGDMFIR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANITEGNWE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGFYGSEIE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PTPCMEWSK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'T', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PIPPDTLTK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ALHPAMGFD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGVNLDNIK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ASYETRSAR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PTNDTYNYK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'T', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ASTQWMCWE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGTCIEPKK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AMFNINLND",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PQWATPAHD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Q', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AIYLHDWNK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AQIGWAEED",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Q', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AFWGQYGTK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ASSYSAPRD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PMWSESGDK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANKEGEWLK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AHVILRGRD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ASVSSVCLR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PFSPACSHE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACLCTTCHK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AHVDAVNVD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATVFDDLCE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'T', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PISGWIIFD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'I', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PYQDDVHDK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Y', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGINMHYWR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PTKGEIFGE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'T', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGYSQMETR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PSMWTFTCK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANLGHTQDR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGIDDFRVE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AMPMHTTME",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PSHWNVTGE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AHAVGTKED",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PMKMIPPTR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'M', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATMKQIDHK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'T', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PNKVSKIFD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PVCRCFYMD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PNGTHPCSR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PNLMVQDGD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AYNTEYRGR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Y', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ASHIHNKCD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATWYCGCRD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'T', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACQVSHCER",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGHTCTLDD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PNWFRCQRE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PFKTLSNNR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PCIIISAKK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGQQNRSRR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AFQETYGTD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'F', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACPFSTDLK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACCDVIWED",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'C', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PQTCIVEFR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'Q', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AVSKCVFNR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'V', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANQHGINDK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'N', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGDYNGGAD",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'D'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PHYFSYIER",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'H', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATTDGDPQR",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'T', 8: 'R'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PSCFKSPGK",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'S', 8: 'K'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PLQGEICQE",
        binding_score=0.248017,
        ic50_nm=5104.2,
        binding_class="non_binder",
        anchor_residues={1: 'L', 8: 'E'},
        anchor_scores={1: 0.7, 8: 0.4},
        source="pssm_predicted",
    ))

    return entries


# Lazy singleton — built on first access
_database: PrecomputedAlleleDatabase | None = None


def get_database() -> PrecomputedAlleleDatabase:
    """Return the precomputed binding database for HLA-B*07:02."""
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
