"""
Legacy pathogen signature database and motif-to-pathogen mapping.

Each entry is a tuple: (sequence, pathogen_name, risk_level, description)
These are longer signature sequences (N-terminal / signal peptides) used for
exact substring matching in the check_biosecurity_before_optimize pipeline.
"""

from __future__ import annotations

_PATHOGEN_SIGNATURES: list[tuple[str, str, str, str]] = [
    # Bacillus anthracis — lethal factor (LF) N-terminal signal peptide
    (
        "MEFKLRILVVSVATLFVSSGYSQHGVRNEQYADLAKR",
        "Bacillus_anthracis",
        "CRITICAL",
        "Anthrax lethal factor signal peptide",
    ),
    # Clostridium botulinum — neurotoxin N-terminal signal peptide
    (
        "MFVKLSFVKILIFQSSQHGVR",
        "Clostridium_botulinum",
        "CRITICAL",
        "Botulinum neurotoxin signal peptide",
    ),
    # Ricinus communis — ricin A-chain N-terminal region
    (
        "AITNLFGRRLDKVKDTSINL",
        "Ricinus_communis",
        "CRITICAL",
        "Ricin A-chain catalytic region",
    ),
    # Shigella dysenteriae — Shiga toxin A subunit
    (
        "MYTILFSILLVSQHGVR",
        "Shigella_dysenteriae",
        "CRITICAL",
        "Shiga toxin A subunit signal peptide",
    ),
    # Vibrio cholerae — cholera toxin A subunit signal peptide
    (
        "MIKLCSFVSILLFFSSFSSG",
        "Vibrio_cholerae",
        "HIGH",
        "Cholera toxin A subunit signal peptide",
    ),
    # Corynebacterium diphtheriae — diphtheria toxin signal peptide
    (
        "MGKKVAVLLLASSVFAHGVR",
        "Corynebacterium_diphtheriae",
        "CRITICAL",
        "Diphtheria toxin signal peptide",
    ),
    # Yersinia pestis — V antigen signal peptide
    (
        "MKKISSVVVMTTMTFSSQHG",
        "Yersinia_pestis",
        "HIGH",
        "Plague V antigen signal peptide",
    ),
]

# Mapping from motif-based HazardMatch names to pathogen organism names
# Used by check_biosecurity_before_optimize to merge motif findings
# into the BiosecurityScreeningResult format.
_MOTIF_TO_PATHOGEN: dict[str, str] = {
    "ricin_A_chain_catalytic": "Ricinus_communis",
    "ricin_A_chain_rRNA": "Ricinus_communis",
    "ricin_A_chain_active": "Ricinus_communis",
    "ricin_B_chain_lectin": "Ricinus_communis",
    "abrin_A_chain": "Abrus_precatorius",
    "abrin_A_chain_rRNA": "Abrus_precatorius",
    "botulinum_zinc_protease": "Clostridium_botulinum",
    "botulinum_receptor": "Clostridium_botulinum",
    "botulinum_light_chain": "Clostridium_botulinum",
    "botulinum_zinc_HELIH": "Clostridium_botulinum",
    "botulinum_catalytic_EIDH": "Clostridium_botulinum",
    "botulinum_zinc_dna": "Clostridium_botulinum",
    "botulinum_catalytic_dna": "Clostridium_botulinum",
    "botulinum_receptor_dna": "Clostridium_botulinum",
    "shiga_toxin_A_subunit": "Shigella_dysenteriae",
    "shiga_toxin_B_subunit": "Shigella_dysenteriae",
    "diphtheria_toxin_ADR": "Corynebacterium_diphtheriae",
    "diphtheria_toxin_catalytic": "Corynebacterium_diphtheriae",
    "tetanus_toxin_zinc": "Clostridium_tetani",
    "tetanus_toxin_light": "Clostridium_tetani",
    "cholera_toxin_A1": "Vibrio_cholerae",
    "cholera_toxin_NAD": "Vibrio_cholerae",
    "anthrax_EF_cyclase": "Bacillus_anthracis",
    "anthrax_EF_calmodulin": "Bacillus_anthracis",
    "anthrax_LF_protease": "Bacillus_anthracis",
    "anthrax_LF_substrate": "Bacillus_anthracis",
    "anthrax_PA_pore": "Bacillus_anthracis",
    "SEB_superantigen": "Staphylococcus_aureus",
    "T2_mycotoxin_target": "Fusarium_spp",
    "influenza_HA_fusion": "Influenza_virus",
    "influenza_HA_receptor": "Influenza_virus",
    "influenza_HA_cleavage": "Influenza_virus",
    "influenza_NA_active": "Influenza_virus",
    "influenza_NA_framework": "Influenza_virus",
    "SARS2_spike_RBD": "SARS_CoV_2",
    "SARS2_spike_fusion": "SARS_CoV_2",
    "SARS2_spike_furin": "SARS_CoV_2",
    "SARS2_spike_heptad": "SARS_CoV_2",
    "HIV_env_V3_loop": "HIV_1",
    "HIV_env_gp41_fusion": "HIV_1",
    "HIV_env_CD4_binding": "HIV_1",
    "ebola_GP1_receptor": "Ebolavirus",
    "ebola_GP_fusion": "Ebolavirus",
    "ebola_GP_mucin": "Ebolavirus",
    "variola_envelope": "Variola_virus",
}
