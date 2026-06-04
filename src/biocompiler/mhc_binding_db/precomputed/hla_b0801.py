"""Precomputed MHC binding data for HLA-B*08:01.

Anchor residues: position 3 (K/R), position 9 (L/I)
Known epitopes: RAKFKQLL, ELRSRYWAI, GPKVKRWWL
Data source: IEDB and SYFPEITHI databases + PSSM predictions

Binding scores are derived from PSSM-based prediction with expected
AUC-ROC of 0.60-0.75. Do not replace experimental validation or
NetMHCpan predictions where available.
"""

from __future__ import annotations

from .. import PrecomputedAlleleDatabase, PrecomputedEntry

ALLELE: str = "HLA-B*08:01"
PEPTIDE_LENGTH: int = 9
ANCHOR_DESCRIPTION: str = "position 3 (K/R), position 9 (L/I)"
KNOWN_EPITOPES: list[str] = ['RAKFKQLL', 'ELRSRYWAI', 'GPKVKRWWL']


def _build_entries() -> list[PrecomputedEntry]:
    """Construct all precomputed entries for HLA-B*08:01."""
    entries: list[PrecomputedEntry] = []

    # --- Binders (strong + moderate) ---
    entries.append(PrecomputedEntry(
        peptide="WRKKRYFNL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NRRWQNYYL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TKKFSLGAL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AKYIILGNL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MKRAGRIKL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FKRRPQLML",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CKNRKEFYL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YRKMTATNL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ERRYMHTYL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AKRDLGSFL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CKKFAWEEL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IKRRKKYWL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LRRGGWSNL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ERRCNGLML",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NRKHTLGEL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TRRAIFATL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WRRMHIGLL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GRNWFINDL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YKRFHPMIL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SRDATPPWL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FRRQPEVML",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CKKEDPSDL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TRKTYKDAL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IKRVTENDL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ERWHEKQFL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GRPGYATGL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LRKNTCHCL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TRRWKQRIL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GRRVTRNTL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GRKWKHDNL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NRKLIRCDL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IKKVCLYVL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KRRCCCLML",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AKPCNGYVL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DRKAQQSRL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SRRWLRYGL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MRRTKFRQL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PRKDWFVVL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ARRQIKVNL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SRGAKIWHL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ARNAKYEQL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RRKGMILYL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LRTCSTICL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HKKMLFFTL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GKKTPFHPL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TRIEATVIL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DRKLHRPAL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'R', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KKKAPQDNL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LKLSFFGYL",
        binding_score=1.0,
        ic50_nm=5.0,
        binding_class="strong_binder",
        anchor_residues={1: 'K', 8: 'L'},
        anchor_scores={1: 1.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ATKFRINYL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AWREWPVCL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TIQHQWAML",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'I', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IFRWNYDAL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'F', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LQRIIYYWL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'Q', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LMVLAWKHL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'M', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GTRVGHWRL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DLRPETYDL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'L', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DHSYDEMPL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'H', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CFKCVPYTL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'F', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DNKSMLKYL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'N', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MYRDSMHIL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'Y', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VAKIWCPML",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IIKWHGFYL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'I', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VVRVVIPDL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'V', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TFPFIQCTL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'F', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SSTMDHDDL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'S', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VLRDPHKKL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'L', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FMKAAYNTL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'M', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MFRSHREDL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'F', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WCRDCTETL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'C', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KLAQWMHNL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'L', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HTRLIPQAL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HVNGVWIEL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'V', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LVRLKDHYL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'V', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LAHQGSMPL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MMRNPPYCL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'M', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MMKELGRCL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'M', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AWRLIACFL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ANNNTCWDL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'N', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QMRFCLTGL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'M', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PAFRTTYPL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'A', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LYRAPQLLL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'Y', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KIRWIKNSL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'I', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QWNYACWNL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NSRWYDQSL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'S', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LIRHKSPKL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'I', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TQKMPPVTL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'Q', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QNFTSWCML",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'N', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SWKRDWYAL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QHKAFASPL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'H', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CHCELMGFL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'H', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NLKFAIYVL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'L', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GLKVYNQWL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'L', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DTSDIPKIL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'T', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ACRYWLFLL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'C', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FWKRGWWFL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'W', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NIKHDIVSL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'I', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EIHKFHRKL",
        binding_score=0.705374,
        ic50_nm=75.6,
        binding_class="moderate_binder",
        anchor_residues={1: 'I', 8: 'L'},
        anchor_scores={1: 0.8, 8: 1.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ELRSRYWAI",
        binding_score=0.656083,
        ic50_nm=119.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'L', 8: 'I'},
        anchor_scores={1: 0.8, 8: 1.3},
        source="known_epitope",
    ))
    entries.append(PrecomputedEntry(
        peptide="GPKVKRWWL",
        binding_score=0.55,
        ic50_nm=400.0,
        binding_class="moderate_binder",
        anchor_residues={1: 'P', 8: 'L'},
        anchor_scores={1: 0.4, 8: 1.5},
        source="known_epitope",
    ))

    # --- Non-binders (weak + non) ---
    entries.append(PrecomputedEntry(
        peptide="TGRGFEPYI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EGRPTFIAI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FGRDKTIHI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YGKFDFNHI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TGLGNKLYI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IGKRPYWGI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGKEYHHMI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGRLQKKWI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGRYVCKDI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KGRVMGKRI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VGRQSFPGI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EGKKACCFI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FGGSVWMPV",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'V'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QGHCPSNII",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WGVLRIEII",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGKVAQEKV",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'V'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGHKRRKMI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGRYRHCHI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KGKPYQYKI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QGCEGRNCI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NGSHICCQI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SGKILCDTI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YGKIWMIAI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IGKEWGWVI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGGEGKANI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YGDIFTHYI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FGRTPNFCI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGRRKFNQI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RGMRWTSWV",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'V'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IGRCYGSYI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EGKRMQEVI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TGTVRNTSI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGHQDIFYI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGKIHDYFI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGFCCPVTI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NGCMIWCNI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MGYVARDSI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TGKVTMLFI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGRAGWYCI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SGTWLEMMI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGIQMWLTV",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'V'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MGRREAAFI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VGRITAEEI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VGPQKPGDI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGMGEPVYV",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'V'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SGKIMAYRI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGNFNHTSV",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'V'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QGRGEKVEI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SGTNWPHQI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGKQEVCII",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IGKLSNEQI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KGKANHITI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WGRGDNRLV",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'V'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGMDTEEYI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGKPAANYI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGKETQQRI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VGKAKSIAI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NGVNYASQI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TGRWVDGFI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CGTGMELCI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IGRLLYRDI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VGESVDCNI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGMDKHQIV",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'V'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGKNDECVI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGKWTSQQI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EGRTGVGMI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGKAIRINI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGFQNRSPI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGVFDSHNI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QGRIAEKSI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGRLVSCLI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MGRVMCRQI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WGRLTHYNI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RGRLIGCNI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LGKMFAEGI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QGREGFLGI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TGVGHAASI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TGCKEVYDI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGRRTFQDI",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'I'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LGSYTIITV",
        binding_score=0.499595,
        ic50_nm=503.1,
        binding_class="weak_binder",
        anchor_residues={1: 'G', 8: 'V'},
        anchor_scores={1: 0.5, 8: 1.3},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGNDKAGSK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MGCANHHDK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RGMNDVHFK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KGRTYQSEK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RGQNYNGEK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CGCCQEMYK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LGCLRAKYK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RGRVGPSAK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VGQDDYTEK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QGKWMTFDK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HGSWVESMK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGILICFIK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WGAHWYADK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EGANWNHYK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LGEELMWWK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SGPTVFYIK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KGLPGVHWK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LGKVMHDDK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NGHYMCNWK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGWYVHHWK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGAGWWCNK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EGMQNMGMK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QGKTYFFPK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YGSVEYKAK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TGAKPRNAK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MGLRTCTEK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGQVYGCTK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGKMLRAIK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IGKNPSMVK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGKYDIYCK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGKDLKIAK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VGKVDYNYK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KGHYRVYIK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RGKQIDEWK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MGRVYCIRK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QGKVHDFIK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGYGNSIRK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="CGPLASVKK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="AGCWEMKRK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGGKLLVIK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GGRQFWVMK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WGTFGKLNK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QGAPNPHHK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WGPSWITQK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGHLGWNIK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DGTTHAVDK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MGKDTRRHK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KGRLTKRQK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PGTSGHPEK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NGAKQKEGK",
        binding_score=0.205575,
        ic50_nm=7545.7,
        binding_class="non_binder",
        anchor_residues={1: 'G', 8: 'K'},
        anchor_scores={1: 0.5, 8: 0.5},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RENEDSSVC",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'C'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PDYKRNLCT",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'T'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FEKEDKPRM",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WELPHQDFG",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'G'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GERCVEEPQ",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Q'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DEWQACNGG",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'G'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LEFFICEFT",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'T'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EEWASFCAT",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'T'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DDWRNMRWY",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Y'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDNYYIEGW",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'W'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MEDSFNHWM",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SEKNPMADN",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'N'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MEHDQLERN",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'N'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADSLMWKHY",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Y'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HENSGDTYP",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'P'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VEREIQDDN",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'N'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="REWAALYKH",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'H'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HETNSRAKR",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'R'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KEQMLKPAN",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'N'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TEWAACYGH",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'H'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QENTCCKYY",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Y'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HEGINDFKP",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'P'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="ADKLKCRGM",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'M'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FESVEKIAH",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'H'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SERAEWFGG",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'G'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DECMHIEFP",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'P'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TEQWMDIWS",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'S'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="LEKNCTMPA",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'A'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="REWFAIAVF",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'F'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RENYPQMTF",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'F'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DDFALQCRP",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'P'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EERFKKDFG",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'G'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VDTHKIPRS",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'S'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MEKLQDKWS",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'S'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDYKTASRA",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'A'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QEDPYRMMA",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'A'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TDHTVLMPM",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'M'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QDSLHPSAM",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'M'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GELKFKKCH",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'H'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PEDRPHCTR",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'R'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="REQHFDAAM",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDQYNSINQ",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'Q'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PDKDRMCWA",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'A'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PDHNNLTRW",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'W'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="HERRNYAPG",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'G'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="IECGYMHTY",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Y'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KEDYFKYQT",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'T'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KDRHMSYQS",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'S'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NEVHSMKHM",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NEDCNMFWQ",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Q'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDGQLDLRH",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'H'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="YETSSYCQQ",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Q'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KDIMDKATC",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'C'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="KDIWTVPRR",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'R'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="MDVPCGEFA",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'A'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="SDVRRFVLT",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'T'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QEQSHEDFN",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'N'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="EENCGQKIR",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'R'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RDRLIKDCP",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'P'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="QDARLFKET",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'T'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="RETITHTPC",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'C'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="DEQFAFSKR",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'R'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="PEKEHQYSW",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'W'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="WELSHSFNA",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'A'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="GERFPAIPA",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'A'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VDRNASLGC",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'C'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="TERYEESDM",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'M'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="NDRLYEKNN",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'N'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="VDVFYQHYG",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'D', 8: 'G'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))
    entries.append(PrecomputedEntry(
        peptide="FERRDTMWY",
        binding_score=0.193696,
        ic50_nm=8418.1,
        binding_class="non_binder",
        anchor_residues={1: 'E', 8: 'Y'},
        anchor_scores={1: 0.3, 8: 0.8},
        source="pssm_predicted",
    ))

    return entries


# Lazy singleton — built on first access
_database: PrecomputedAlleleDatabase | None = None


def get_database() -> PrecomputedAlleleDatabase:
    """Return the precomputed binding database for HLA-B*08:01."""
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
