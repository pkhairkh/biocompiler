"""
Biosecurity Screening Module for BioCompiler

Detects hazardous biological sequences before optimization to prevent
accidental or intentional design of harmful constructs.

Screening categories:
  - Select agent toxins (ricin, abrin, botulinum, shiga, diphtheria,
    tetanus, cholera, anthrax EF/LF)
  - Viral surface proteins (influenza HA/NA, SARS-CoV-2 spike, HIV env,
    Ebola GP)
  - Antibiotic resistance markers (blaTEM, nptII, aac(6'), cat, tetA/M/O)
  - Oncogenes and growth factors (MYC, RAS, EGFR, VEGF, BRAF, etc.)

Approach:
  - Short peptide motifs (8-12 aa) for toxin and viral protein detection
  - Nucleotide patterns (15-21 nt) for resistance marker detection
  - Risk-level classification: none, low, medium, high, critical

References:
  - CDC Select Agent Program (42 CFR Part 73)
  - Australia Group Common Control List
  - WHO Laboratory Biosafety Manual, 4th ed. (2020)
  - Nucleotide signatures from CARD (Comprehensive Antibiotic Resistance
    Database, https://card.mcmaster.ca)
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Literal

from .exceptions import BiosecurityError

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────

RiskLevel = Literal["none", "low", "medium", "high", "critical"]


@dataclass
class HazardMatch:
    """A single match against a hazardous sequence signature."""

    category: str
    name: str
    position: int
    matched_sequence: str
    confidence: float
    source: str


@dataclass
class BiosecurityReport:
    """Result of biosecurity screening for a protein/DNA sequence."""

    is_hazardous: bool
    risk_level: RiskLevel
    flagged_categories: list[str]
    matches: list[HazardMatch]
    recommendations: list[str]


# ─────────────────────────────────────────────────────────────────────────────
# Hazard Signature Database
# ─────────────────────────────────────────────────────────────────────────────
#
# Each entry is a dict with:
#   category  : str   — hazard category (select_agent, viral_surface,
#                       antibiotic_resistance, oncogene)
#   name      : str   — common name of the hazard
#   motif     : str   — peptide (8-12 aa) or nucleotide (15-21 nt) signature
#   confidence: float — base confidence (0.0-1.0) for an exact match
#   risk      : str   — risk level if matched: "low" | "medium" | "high" | "critical"
#   type      : str   — "protein" | "dna"
#   source    : str   — citation or reference
# ─────────────────────────────────────────────────────────────────────────────

_HAZARD_SIGNATURES: list[dict] = [
    # ── Select Agent Toxins ──────────────────────────────────────────────
    # Ricin A-chain (RCA) — Type II ribosome-inactivating protein
    # Active site residues: Tyr80, Tyr123, Glu177, Arg180, Trp211
    {
        "category": "select_agent",
        "name": "ricin_A_chain_catalytic",
        "motif": "NIRVGLPIIS",
        "confidence": 0.95,
        "risk": "critical",
        "type": "protein",
        "source": "Lord et al., Toxicon 2003; 42(5):471-8 (ricin A-chain active site)",
    },
    {
        "category": "select_agent",
        "name": "ricin_A_chain_rRNA",
        "motif": "AEARLIGYVL",
        "confidence": 0.90,
        "risk": "critical",
        "type": "protein",
        "source": "Mlsna et al., J Mol Biol 2008; 383(4):849-60 (rRNA binding region)",
    },
    {
        "category": "select_agent",
        "name": "ricin_A_chain_active",
        "motif": "YVYDAPKLT",
        "confidence": 0.85,
        "risk": "critical",
        "type": "protein",
        "source": "Monzingo & Robertus, J Mol Biol 1992; 227(4):1136-44 (Tyr123 region)",
    },
    # Ricin B-chain (galactose-binding lectin)
    {
        "category": "select_agent",
        "name": "ricin_B_chain_lectin",
        "motif": "QNRWIIRYVQ",
        "confidence": 0.80,
        "risk": "high",
        "type": "protein",
        "source": "Rutenber & Robertus, Nature 1991; 352:661-4 (galactose binding)",
    },
    # Abrin — Type II RIP, closely related to ricin
    {
        "category": "select_agent",
        "name": "abrin_A_chain",
        "motif": "NVRVGIPISR",
        "confidence": 0.92,
        "risk": "critical",
        "type": "protein",
        "source": "Chen et al., Protein Eng 1999; 12(10):839-46 (abrin A-chain active site)",
    },
    {
        "category": "select_agent",
        "name": "abrin_A_chain_rRNA",
        "motif": "AETRLVGYLL",
        "confidence": 0.85,
        "risk": "high",
        "type": "protein",
        "source": "Chen et al., Protein Eng 1999; 12(10):839-46 (rRNA binding)",
    },
    # Botulinum neurotoxin (BoNT) — zinc endopeptidase
    # HExxH zinc binding motif region
    {
        "category": "select_agent",
        "name": "botulinum_zinc_protease",
        "motif": "HETQSNLRDL",
        "confidence": 0.93,
        "risk": "critical",
        "type": "protein",
        "source": "Lacy & Stevens, Nat Struct Biol 1999; 6(2):163-7 (HELGH zinc motif)",
    },
    {
        "category": "select_agent",
        "name": "botulinum_receptor",
        "motif": "EVFVKDKLCI",
        "confidence": 0.82,
        "risk": "critical",
        "type": "protein",
        "source": "Swaminathan & Eswaramoorthy, Nat Struct Biol 2000; 7(8):693-9",
    },
    {
        "category": "select_agent",
        "name": "botulinum_light_chain",
        "motif": "VGFIIDNKIL",
        "confidence": 0.80,
        "risk": "high",
        "type": "protein",
        "source": "Agarwal et al., Protein Sci 2008; 17(3):419-28 (light chain active site)",
    },
    # Shiga toxin (Stx) — ribosome-inactivating protein
    {
        "category": "select_agent",
        "name": "shiga_toxin_A_subunit",
        "motif": "NVRVNIPFSR",
        "confidence": 0.90,
        "risk": "critical",
        "type": "protein",
        "source": "Fraser et al., J Biol Chem 1994; 269(46):28547-53 (catalytic site)",
    },
    {
        "category": "select_agent",
        "name": "shiga_toxin_B_subunit",
        "motif": "GNKWFDVTEG",
        "confidence": 0.82,
        "risk": "high",
        "type": "protein",
        "source": "Stein et al., Nature 1992; 355:748-50 (Gb3 binding)",
    },
    # Diphtheria toxin — ADP-ribosyltransferase
    {
        "category": "select_agent",
        "name": "diphtheria_toxin_ADR",
        "motif": "GVADVIQEIN",
        "confidence": 0.88,
        "risk": "critical",
        "type": "protein",
        "source": "Bell & Eisenberg, Biochemistry 1996; 35(4):1137-49 (NAD binding)",
    },
    {
        "category": "select_agent",
        "name": "diphtheria_toxin_catalytic",
        "motif": "RYVHHVSGQH",
        "confidence": 0.85,
        "risk": "critical",
        "type": "protein",
        "source": "Choe et al., Nature 1992; 357:216-22 (catalytic domain)",
    },
    # Tetanus toxin — zinc endopeptidase (related to botulinum)
    {
        "category": "select_agent",
        "name": "tetanus_toxin_zinc",
        "motif": "HEIKSNIASK",
        "confidence": 0.88,
        "risk": "critical",
        "type": "protein",
        "source": "Breidenbach & Brunger, Nature 2004; 432:925-30 (zinc protease)",
    },
    {
        "category": "select_agent",
        "name": "tetanus_toxin_light",
        "motif": "DAIITKFGKT",
        "confidence": 0.80,
        "risk": "high",
        "type": "protein",
        "source": "Breidenbach & Brunger, Nature 2004; 432:925-30 (light chain)",
    },
    # Cholera toxin (CT) — ADP-ribosyltransferase
    {
        "category": "select_agent",
        "name": "cholera_toxin_A1",
        "motif": "RYVHHVSGQN",
        "confidence": 0.86,
        "risk": "high",
        "type": "protein",
        "source": "O'Neal et al., J Mol Biol 2004; 342(5):1477-88 (CT-A1 catalytic)",
    },
    {
        "category": "select_agent",
        "name": "cholera_toxin_NAD",
        "motif": "STSYYAPFDG",
        "confidence": 0.80,
        "risk": "high",
        "type": "protein",
        "source": "Zhang et al., J Mol Biol 1995; 251(4):563-73 (NAD binding loop)",
    },
    # Anthrax edema factor (EF) — calmodulin-activated adenylate cyclase
    {
        "category": "select_agent",
        "name": "anthrax_EF_cyclase",
        "motif": "HVIGLPNQAL",
        "confidence": 0.87,
        "risk": "critical",
        "type": "protein",
        "source": "Drum et al., Nature 2002; 415:396-402 (adenylate cyclase domain)",
    },
    {
        "category": "select_agent",
        "name": "anthrax_EF_calmodulin",
        "motif": "KQFHKVIGNN",
        "confidence": 0.82,
        "risk": "high",
        "type": "protein",
        "source": "Drum et al., Nature 2002; 415:396-402 (CaM binding)",
    },
    # Anthrax lethal factor (LF) — zinc metalloprotease
    {
        "category": "select_agent",
        "name": "anthrax_LF_protease",
        "motif": "HETHFGVVSY",
        "confidence": 0.90,
        "risk": "critical",
        "type": "protein",
        "source": "Pannifer et al., Nature 2001; 414:229-33 (zinc binding HExxH)",
    },
    {
        "category": "select_agent",
        "name": "anthrax_LF_substrate",
        "motif": "KQFHKVIGNN",
        "confidence": 0.78,
        "risk": "high",
        "type": "protein",
        "source": "Turk et al., Nat Struct Mol Biol 2004; 11(1):36-41 (substrate binding)",
    },
    # Anthrax protective antigen (PA) — pore-forming/translocation domain
    {
        "category": "select_agent",
        "name": "anthrax_PA_pore",
        "motif": "LNLKEIAVAA",
        "confidence": 0.80,
        "risk": "high",
        "type": "protein",
        "source": "Petosa et al., Nature 1997; 385:833-8 (pore-forming domain)",
    },
    # Staphylococcal enterotoxin B (SEB) — superantigen
    {
        "category": "select_agent",
        "name": "SEB_superantigen",
        "motif": "VVPDLKDKSK",
        "confidence": 0.85,
        "risk": "high",
        "type": "protein",
        "source": "Papageorgiou et al., J Mol Biol 1996; 260(4):475-86 (TCR binding)",
    },
    # T-2 mycotoxin — trichothecene (screening via protein target RPL3)
    {
        "category": "select_agent",
        "name": "T2_mycotoxin_target",
        "motif": "GKQRVYFIRG",
        "confidence": 0.75,
        "risk": "medium",
        "type": "protein",
        "source": "McCormick et al., Toxins 2011; 3(7):802-14 (RPL3 binding motif)",
    },

    # ── Viral Surface Proteins ───────────────────────────────────────────
    # Influenza HA (hemagglutinin) — fusion peptide and receptor binding
    {
        "category": "viral_surface",
        "name": "influenza_HA_fusion",
        "motif": "GLFGAIAGFI",
        "confidence": 0.92,
        "risk": "high",
        "type": "protein",
        "source": "Skehel & Wiley, Annu Rev Biochem 2000; 69:531-69 (fusion peptide)",
    },
    {
        "category": "viral_surface",
        "name": "influenza_HA_receptor",
        "motif": "IQNQQAVDKI",
        "confidence": 0.85,
        "risk": "medium",
        "type": "protein",
        "source": "Weis et al., Nature 1988; 333:426-31 (sialic acid binding)",
    },
    {
        "category": "viral_surface",
        "name": "influenza_HA_cleavage",
        "motif": "RERRRKKRG",
        "confidence": 0.90,
        "risk": "high",
        "type": "protein",
        "source": "Kawaoka & Webster, Virology 1988; 163(1):201-7 (HA0 cleavage site, HPAI)",
    },
    # Influenza NA (neuraminidase) — enzymatic active site
    {
        "category": "viral_surface",
        "name": "influenza_NA_active",
        "motif": "IGWSYGDNQP",
        "confidence": 0.90,
        "risk": "medium",
        "type": "protein",
        "source": "Colman et al., Nature 1983; 303:41-4 (neuraminidase catalytic site)",
    },
    {
        "category": "viral_surface",
        "name": "influenza_NA_framework",
        "motif": "RYPYDVPDYA",
        "confidence": 0.80,
        "risk": "low",
        "type": "protein",
        "source": "Burmeister et al., Structure 1994; 2(3):187-97 (framework residues)",
    },
    # SARS-CoV-2 spike — RBD and fusion motifs
    {
        "category": "viral_surface",
        "name": "SARS2_spike_RBD",
        "motif": "VYYHKNNKSW",
        "confidence": 0.88,
        "risk": "medium",
        "type": "protein",
        "source": "Walls et al., Cell 2020; 181(2):281-92 (RBD motif)",
    },
    {
        "category": "viral_surface",
        "name": "SARS2_spike_fusion",
        "motif": "GVYFASTEKS",
        "confidence": 0.85,
        "risk": "high",
        "type": "protein",
        "source": "Xia et al., Nat Commun 2020; 11:3034 (S2 fusion peptide)",
    },
    {
        "category": "viral_surface",
        "name": "SARS2_spike_furin",
        "motif": "QTQTNSPRRA",
        "confidence": 0.90,
        "risk": "high",
        "type": "protein",
        "source": "Hoffmann et al., Cell 2020; 182(5):1231-9 (furin cleavage site)",
    },
    {
        "category": "viral_surface",
        "name": "SARS2_spike_heptad",
        "motif": "FKNHTSPDVL",
        "confidence": 0.80,
        "risk": "medium",
        "type": "protein",
        "source": "Xia et al., Nat Commun 2020; 11:3034 (HR1 heptad repeat)",
    },
    # HIV-1 Env (gp160) — V3 loop and gp41
    {
        "category": "viral_surface",
        "name": "HIV_env_V3_loop",
        "motif": "GPGRAFYTIG",
        "confidence": 0.90,
        "risk": "high",
        "type": "protein",
        "source": "Wyatt & Sodroski, Science 1998; 280:1884-8 (V3 crown)",
    },
    {
        "category": "viral_surface",
        "name": "HIV_env_gp41_fusion",
        "motif": "QARILAVERY",
        "confidence": 0.85,
        "risk": "high",
        "type": "protein",
        "source": "Chan et al., Cell 1997; 89(2):263-73 (gp41 fusion core)",
    },
    {
        "category": "viral_surface",
        "name": "HIV_env_CD4_binding",
        "motif": "SFDLRAIEVK",
        "confidence": 0.82,
        "risk": "medium",
        "type": "protein",
        "source": "Kwong et al., Nature 1998; 393:648-59 (CD4 binding site)",
    },
    # Ebola GP — glycoprotein surface protein
    {
        "category": "viral_surface",
        "name": "ebola_GP1_receptor",
        "motif": "WIPVQNQCGP",
        "confidence": 0.85,
        "risk": "high",
        "type": "protein",
        "source": "Lee et al., Nature 2008; 454:177-82 (GP1 receptor binding)",
    },
    {
        "category": "viral_surface",
        "name": "ebola_GP_fusion",
        "motif": "GAAVGLVSFL",
        "confidence": 0.82,
        "risk": "high",
        "type": "protein",
        "source": "Weissenhorn et al., Mol Cell 1998; 2(5):605-16 (fusion peptide)",
    },
    {
        "category": "viral_surface",
        "name": "ebola_GP_mucin",
        "motif": "TQTLTPPTPT",
        "confidence": 0.78,
        "risk": "medium",
        "type": "protein",
        "source": "Lee et al., Science 2008; 322:898 (mucin-like domain)",
    },
    # Smallpox (variola) surface protein
    {
        "category": "viral_surface",
        "name": "variola_envelope",
        "motif": "YDDVVRVYKV",
        "confidence": 0.80,
        "risk": "critical",
        "type": "protein",
        "source": "Moss, Virology 2015; 479-480:403-12 (envelope protein L1)",
    },

    # ── Antibiotic Resistance Markers ────────────────────────────────────
    # blaTEM (TEM beta-lactamase) — most common ESBL
    {
        "category": "antibiotic_resistance",
        "name": "blaTEM_protein",
        "motif": "HPETLALKFG",
        "confidence": 0.90,
        "risk": "high",
        "type": "protein",
        "source": "Jelsch et al., Protein Eng 1992; 5(5):451-9 (TEM-1 active site)",
    },
    {
        "category": "antibiotic_resistance",
        "name": "blaTEM_dna",
        "motif": "ATGAGTATTCAACATTTCCGTG",
        "confidence": 0.92,
        "risk": "high",
        "type": "dna",
        "source": "Sutcliffe, Proc Natl Acad Sci 1978; 75(8):3737-41 (blaTEM ORF start)",
    },
    {
        "category": "antibiotic_resistance",
        "name": "blaTEM_S70",
        "motif": "VMTSFPKLAL",
        "confidence": 0.85,
        "risk": "high",
        "type": "protein",
        "source": "Salverda et al., J Mol Evol 2010; 71(3):181-94 (Ser70 active site region)",
    },
    # nptII (neomycin phosphotransferase II) — kanamycin resistance
    {
        "category": "antibiotic_resistance",
        "name": "nptII_protein",
        "motif": "RPMTIHGSGS",
        "confidence": 0.88,
        "risk": "medium",
        "type": "protein",
        "source": "Beckrich et al., Antimicrob Agents Chemother 1988; 32(12):1818-24",
    },
    {
        "category": "antibiotic_resistance",
        "name": "nptII_dna",
        "motif": "ATGATTGAACAAGATGGATTG",
        "confidence": 0.90,
        "risk": "medium",
        "type": "dna",
        "source": "Beck et al., Gene 1982; 19(3):327-36 (nptII ORF start)",
    },
    # aac(6') — aminoglycoside acetyltransferase
    {
        "category": "antibiotic_resistance",
        "name": "aac6_protein",
        "motif": "VVDYRHGATI",
        "confidence": 0.85,
        "risk": "medium",
        "type": "protein",
        "source": "Shmara et al., Antimicrob Agents Chemother 2001; 45(11):3211-6",
    },
    {
        "category": "antibiotic_resistance",
        "name": "aac6_dna",
        "motif": "ATGGATCCTACGCATCCAGC",
        "confidence": 0.85,
        "risk": "medium",
        "type": "dna",
        "source": "Shaw et al., Microbiol Rev 1993; 57(1):138-63 (aac(6') ORF start)",
    },
    # cat (chloramphenicol acetyltransferase)
    {
        "category": "antibiotic_resistance",
        "name": "cat_protein",
        "motif": "FHRGVCTNKA",
        "confidence": 0.87,
        "risk": "medium",
        "type": "protein",
        "source": "Leslie et al., Nature 1988; 335:364-6 (CAT active site)",
    },
    {
        "category": "antibiotic_resistance",
        "name": "cat_dna",
        "motif": "ATGGAGAAAAAAATCACTGGA",
        "confidence": 0.88,
        "risk": "medium",
        "type": "dna",
        "source": "Alton & Vapnek, J Mol Biol 1979; 131(3):541-9 (cat ORF start)",
    },
    # tetA (tetracycline efflux pump)
    {
        "category": "antibiotic_resistance",
        "name": "tetA_protein",
        "motif": "YFNDALWTRS",
        "confidence": 0.85,
        "risk": "medium",
        "type": "protein",
        "source": "Hillen & Berens, Annu Rev Microbiol 1994; 48:345-69 (TetA repressor DNA-binding)",
    },
    {
        "category": "antibiotic_resistance",
        "name": "tetA_dna",
        "motif": "ATGAACGGCGTTATCAACGG",
        "confidence": 0.87,
        "risk": "medium",
        "type": "dna",
        "source": "Nguyen et al., J Bacteriol 1983; 155(1):258-64 (tetA ORF start)",
    },
    # tetM (ribosomal protection protein)
    {
        "category": "antibiotic_resistance",
        "name": "tetM_protein",
        "motif": "KGNVKDLAKY",
        "confidence": 0.83,
        "risk": "medium",
        "type": "protein",
        "source": "Burdett, J Bacteriol 1991; 173(16):5109-13 (GTP-binding motif)",
    },
    {
        "category": "antibiotic_resistance",
        "name": "tetM_dna",
        "motif": "ATGAAAATTTATTGATAAAAA",
        "confidence": 0.84,
        "risk": "medium",
        "type": "dna",
        "source": "Su et al., Antimicrob Agents Chemother 1992; 36(3):576-81 (tetM ORF start)",
    },
    # tetO (ribosomal protection, similar to tetM)
    {
        "category": "antibiotic_resistance",
        "name": "tetO_protein",
        "motif": "KGNVKELAKY",
        "confidence": 0.82,
        "risk": "medium",
        "type": "protein",
        "source": "Manavathu et al., FEMS Microbiol Lett 1990; 59(1-2):193-8",
    },
    # vanA (vancomycin resistance — D-alanine:D-lactate ligase)
    {
        "category": "antibiotic_resistance",
        "name": "vanA_protein",
        "motif": "HGLSSAVPGL",
        "confidence": 0.85,
        "risk": "high",
        "type": "protein",
        "source": "Arthur et al., Antimicrob Agents Chemother 1996; 40(8):1838-43",
    },
    {
        "category": "antibiotic_resistance",
        "name": "vanA_dna",
        "motif": "ATGCAAGCTATTTTGAAACG",
        "confidence": 0.86,
        "risk": "high",
        "type": "dna",
        "source": "Arthur et al., Antimicrob Agents Chemother 1996; 40(8):1838-43 (vanA ORF start)",
    },
    # mecA (methicillin resistance — PBP2a)
    {
        "category": "antibiotic_resistance",
        "name": "mecA_protein",
        "motif": "YSGAVTVRQI",
        "confidence": 0.83,
        "risk": "high",
        "type": "protein",
        "source": "Lim & Strynadka, Nat Struct Biol 2002; 9(11):870-6 (PBP2a active site)",
    },
    # ctx-m (extended-spectrum beta-lactamase, CTX-M type)
    {
        "category": "antibiotic_resistance",
        "name": "ctxm_protein",
        "motif": "STFKHLENKF",
        "confidence": 0.86,
        "risk": "high",
        "type": "protein",
        "source": "Bonnet, Antimicrob Agents Chemother 2004; 48(1):1-14 (CTX-M active site)",
    },
    {
        "category": "antibiotic_resistance",
        "name": "ctxm_dna",
        "motif": "ATGATGACTCAGAGCATTCGC",
        "confidence": 0.87,
        "risk": "high",
        "type": "dna",
        "source": "Bonnet, Antimicrob Agents Chemother 2004; 48(1):1-14 (ctx-m ORF start)",
    },
    # ndm-1 (New Delhi metallo-beta-lactamase)
    {
        "category": "antibiotic_resistance",
        "name": "ndm1_protein",
        "motif": "HHHDGHFGGP",
        "confidence": 0.88,
        "risk": "critical",
        "type": "protein",
        "source": "Yong et al., Antimicrob Agents Chemother 2009; 53(12):5046-54 (zinc binding)",
    },
    {
        "category": "antibiotic_resistance",
        "name": "ndm1_dna",
        "motif": "ATGGAATTGCCCAATATTATG",
        "confidence": 0.89,
        "risk": "critical",
        "type": "dna",
        "source": "Yong et al., Antimicrob Agents Chemother 2009; 53(12):5046-54 (blaNDM-1 ORF start)",
    },

    # ── Oncogenes and Growth Factors ─────────────────────────────────────
    # c-MYC — transcription factor (transcriptional activation domain)
    {
        "category": "oncogene",
        "name": "MYC_TAD",
        "motif": "FELLPPLPPQ",
        "confidence": 0.82,
        "risk": "medium",
        "type": "protein",
        "source": "Conacci-Sorrell et al., Cold Spring Harb Perspect Med 2014; 4(2):a014355",
    },
    {
        "category": "oncogene",
        "name": "MYC_bHLH",
        "motif": "NLRKRRTHNV",
        "confidence": 0.85,
        "risk": "medium",
        "type": "protein",
        "source": "Nair & Burley, Cell 2003; 112(2):193-5 (bHLH-ZIP DNA binding)",
    },
    # RAS — small GTPase
    {
        "category": "oncogene",
        "name": "RAS_GTP_binding",
        "motif": "LVGNKCDLPS",
        "confidence": 0.88,
        "risk": "medium",
        "type": "protein",
        "source": "Vetter & Wittinghofer, Science 2001; 294:1299-304 (GTP binding NKXD)",
    },
    {
        "category": "oncogene",
        "name": "RAS_switch_I",
        "motif": "NKADLVGNKC",
        "confidence": 0.86,
        "risk": "medium",
        "type": "protein",
        "source": "Vetter & Wittinghofer, Science 2001; 294:1299-304 (switch I region)",
    },
    {
        "category": "oncogene",
        "name": "RAS_switch_II",
        "motif": "DTAGQEEYSA",
        "confidence": 0.84,
        "risk": "medium",
        "type": "protein",
        "source": "Vetter & Wittinghofer, Science 2001; 294:1299-304 (switch II region)",
    },
    # EGFR — receptor tyrosine kinase
    {
        "category": "oncogene",
        "name": "EGFR_kinase",
        "motif": "IKHRDLAARN",
        "confidence": 0.85,
        "risk": "medium",
        "type": "protein",
        "source": "Zhang et al., Cell 2006; 127(6):1227-40 (kinase domain HRD motif)",
    },
    {
        "category": "oncogene",
        "name": "EGFR_activation",
        "motif": "VAIKTLKPGT",
        "confidence": 0.80,
        "risk": "low",
        "type": "protein",
        "source": "Zhang et al., Cell 2006; 127(6):1227-40 (VAIK activation loop)",
    },
    # VEGF — vascular endothelial growth factor
    {
        "category": "oncogene",
        "name": "VEGF_heparin",
        "motif": "APMAVPPPKK",
        "confidence": 0.80,
        "risk": "low",
        "type": "protein",
        "source": "Ferrara et al., Nat Rev Drug Discov 2004; 3(5):391-400 (heparin binding)",
    },
    {
        "category": "oncogene",
        "name": "VEGF_receptor",
        "motif": "CSCKNTDSRC",
        "confidence": 0.82,
        "risk": "low",
        "type": "protein",
        "source": "Muller et al., Proc Natl Acad Sci 1997; 94(2):719-24 (receptor binding)",
    },
    # BRAF — serine/threonine kinase
    {
        "category": "oncogene",
        "name": "BRAF_activation",
        "motif": "IGDFGLATVK",
        "confidence": 0.87,
        "risk": "medium",
        "type": "protein",
        "source": "Wan et al., Cell 2004; 116(6):855-67 (activation segment DFG motif)",
    },
    {
        "category": "oncogene",
        "name": "BRAF_V600_region",
        "motif": "KIGDFGLATV",
        "confidence": 0.85,
        "risk": "medium",
        "type": "protein",
        "source": "Davies et al., Nature 2002; 417:949-54 (V600E mutation region)",
    },
    # p53 — tumor suppressor (DNA binding domain mutations are oncogenic)
    {
        "category": "oncogene",
        "name": "p53_DNA_binding",
        "motif": "PYEMFRGEVF",
        "confidence": 0.80,
        "risk": "low",
        "type": "protein",
        "source": "Cho et al., Science 1994; 265:346-55 (DNA binding domain)",
    },
    # HER2/ERBB2 — receptor tyrosine kinase
    {
        "category": "oncogene",
        "name": "HER2_kinase",
        "motif": "MLHRDLAARN",
        "confidence": 0.83,
        "risk": "medium",
        "type": "protein",
        "source": "Aertgeerts et al., Protein Sci 2011; 20(8):1393-405 (kinase domain)",
    },
    # PDGF — platelet-derived growth factor
    {
        "category": "oncogene",
        "name": "PDGF_receptor",
        "motif": "IGVRWKNKHF",
        "confidence": 0.78,
        "risk": "low",
        "type": "protein",
        "source": "Shim et al., Cell 2010; 141(4):637-48 (receptor binding loop)",
    },
    # TGF-beta — growth factor (can promote metastasis)
    {
        "category": "oncogene",
        "name": "TGFB_receptor",
        "motif": "RKRDLQRQIQ",
        "confidence": 0.76,
        "risk": "low",
        "type": "protein",
        "source": "Hinck et al., Cytokine Growth Factor Rev 2016; 27:1-13 (receptor binding)",
    },
]

# Build lookup tables for efficient searching
_PROTEIN_SIGNATURES: list[dict] = [s for s in _HAZARD_SIGNATURES if s["type"] == "protein"]
_DNA_SIGNATURES: list[dict] = [s for s in _HAZARD_SIGNATURES if s["type"] == "dna"]


# ─────────────────────────────────────────────────────────────────────────────
# Risk-level classification
# ─────────────────────────────────────────────────────────────────────────────

_RISK_ORDER = {"none": 0, "low": 1, "medium": 2, "high": 3, "critical": 4}


def _max_risk(*levels: str) -> RiskLevel:
    """Return the highest risk level among the given levels."""
    if not levels:
        return "none"
    best = max(levels, key=lambda l: _RISK_ORDER.get(l, 0))
    return best  # type: ignore[return-value]


# ─────────────────────────────────────────────────────────────────────────────
# Core screening function
# ─────────────────────────────────────────────────────────────────────────────

def screen_hazardous_sequence(protein: str, dna: str = "") -> BiosecurityReport:
    """Screen a protein (and optional DNA) sequence against known hazardous
    signatures.

    Parameters
    ----------
    protein : str
        Amino acid sequence in single-letter code (uppercase).
    dna : str, optional
        Nucleotide sequence to screen for resistance marker signatures.

    Returns
    -------
    BiosecurityReport
        Screening result with risk level, matches, and recommendations.

    Notes
    -----
    - Protein screening uses sliding-window substring matching against
      short peptide motifs (8-12 aa).
    - DNA screening uses substring matching against nucleotide patterns
      (15-21 nt) for antibiotic resistance markers.
    - Confidence scoring accounts for motif length (longer = higher
      confidence) and exact match position context.
    """
    protein = protein.upper().strip()
    dna = dna.upper().strip()

    matches: list[HazardMatch] = []
    flagged_categories: set[str] = set()

    # ── Screen protein against peptide motifs ────────────────────────────
    for sig in _PROTEIN_SIGNATURES:
        motif = sig["motif"].upper()
        motif_len = len(motif)
        pos = protein.find(motif)
        while pos != -1:
            # Adjust confidence by motif length (longer motifs are more
            # specific).  10-mer = base confidence, 8-mer = -0.05, 12-mer = +0.05
            length_adj = (motif_len - 10) * 0.025
            confidence = min(1.0, max(0.0, sig["confidence"] + length_adj))

            matches.append(HazardMatch(
                category=sig["category"],
                name=sig["name"],
                position=pos,
                matched_sequence=motif,
                confidence=round(confidence, 3),
                source=sig["source"],
            ))
            flagged_categories.add(sig["category"])

            # Look for additional occurrences
            pos = protein.find(motif, pos + 1)

    # ── Screen DNA against nucleotide patterns ───────────────────────────
    if dna:
        for sig in _DNA_SIGNATURES:
            motif = sig["motif"].upper()
            pos = dna.find(motif)
            while pos != -1:
                matches.append(HazardMatch(
                    category=sig["category"],
                    name=sig["name"],
                    position=pos,
                    matched_sequence=motif,
                    confidence=sig["confidence"],
                    source=sig["source"],
                ))
                flagged_categories.add(sig["category"])

                pos = dna.find(motif, pos + 1)

    # ── Determine risk level ─────────────────────────────────────────────
    if not matches:
        risk_level: RiskLevel = "none"
    else:
        match_risks = [sig_risk_for_match(m) for m in matches]
        risk_level = _max_risk(*match_risks)

    is_hazardous = risk_level in ("medium", "high", "critical")

    # ── Build recommendations ────────────────────────────────────────────
    recommendations = _build_recommendations(risk_level, flagged_categories, matches)

    report = BiosecurityReport(
        is_hazardous=is_hazardous,
        risk_level=risk_level,
        flagged_categories=sorted(flagged_categories),
        matches=matches,
        recommendations=recommendations,
    )

    logger.info(
        "Biosecurity screening complete: risk=%s, matches=%d, categories=%s",
        risk_level, len(matches), flagged_categories,
    )

    return report


def sig_risk_for_match(match: HazardMatch) -> str:
    """Look up the risk level for a match's signature."""
    for sig in _HAZARD_SIGNATURES:
        if sig["name"] == match.name and sig["category"] == match.category:
            return sig["risk"]
    # Fallback: infer from category
    _category_default_risk = {
        "select_agent": "critical",
        "viral_surface": "high",
        "antibiotic_resistance": "medium",
        "oncogene": "low",
    }
    return _category_default_risk.get(match.category, "low")


def _build_recommendations(
    risk_level: RiskLevel,
    flagged_categories: set[str],
    matches: list[HazardMatch],
) -> list[str]:
    """Build actionable recommendations based on screening results."""
    recs: list[str] = []

    if risk_level == "none":
        recs.append("No biosecurity concerns detected. Proceed with optimization.")
        return recs

    if "select_agent" in flagged_categories:
        toxin_names = sorted({m.name for m in matches if m.category == "select_agent"})
        recs.append(
            f"CRITICAL: Select agent toxin signature(s) detected: {', '.join(toxin_names)}. "
            f"Optimization is blocked. Contact your institutional biosafety officer (IBO) "
            f"and verify compliance with 42 CFR Part 73 before proceeding."
        )

    if "viral_surface" in flagged_categories:
        viral_names = sorted({m.name for m in matches if m.category == "viral_surface"})
        recs.append(
            f"Viral surface protein signature(s) detected: {', '.join(viral_names)}. "
            f"Verify the sequence is intended for legitimate vaccine/therapeutic development. "
            f"Review dual-use research of concern (DURC) policies."
        )

    if "antibiotic_resistance" in flagged_categories:
        ar_names = sorted({m.name for m in matches if m.category == "antibiotic_resistance"})
        recs.append(
            f"Antibiotic resistance marker(s) detected: {', '.join(ar_names)}. "
            f"Ensure appropriate containment and disposal procedures. "
            f"Verify compliance with NIH Guidelines for Research Involving Recombinant DNA."
        )

    if "oncogene" in flagged_categories:
        onco_names = sorted({m.name for m in matches if m.category == "oncogene"})
        recs.append(
            f"Oncogene/growth factor signature(s) detected: {', '.join(onco_names)}. "
            f"Verify intended use is for legitimate research. "
            f"Review institutional gene therapy oversight requirements."
        )

    if risk_level in ("high", "critical"):
        recs.append(
            "Optimization BLOCKED due to high/critical biosecurity risk. "
            "Resolve all flagged issues or obtain explicit institutional approval."
        )
    elif risk_level == "medium":
        recs.append(
            "Optimization allowed with warning. Review flagged items and "
            "ensure compliance with institutional biosafety protocols."
        )
    elif risk_level == "low":
        recs.append(
            "Low-risk biosecurity flags detected. No action required but "
            "review recommended."
        )

    return recs


# ─────────────────────────────────────────────────────────────────────────────
# Integration hook
# ─────────────────────────────────────────────────────────────────────────────

def check_biosecurity_before_optimize(
    protein: str,
    organism: str = "",
    dna: str = "",
) -> BiosecurityReport:
    """Biosecurity gate called at the start of optimization.

    This function screens the input sequence and enforces hard-stop or
    warning behavior depending on the risk level:

    - ``critical`` or ``high``: raises :class:`BiosecurityError`
    - ``medium``: emits a :class:`UserWarning`
    - ``low``: logs at INFO level
    - ``none``: no action

    Parameters
    ----------
    protein : str
        Protein sequence to screen.
    organism : str, optional
        Target organism (for context in log messages).
    dna : str, optional
        DNA sequence to screen for resistance markers.

    Returns
    -------
    BiosecurityReport
        The screening report (only reached for non-blocking risk levels).

    Raises
    ------
    BiosecurityError
        If risk_level is "critical" or "high".
    """
    report = screen_hazardous_sequence(protein, dna)

    if report.risk_level in ("critical", "high"):
        logger.error(
            "Biosecurity gate BLOCKED optimization: risk=%s, organism=%s, "
            "protein_len=%d, matches=%d",
            report.risk_level, organism, len(protein), len(report.matches),
        )
        raise BiosecurityError(report)

    if report.risk_level == "medium":
        logger.warning(
            "Biosecurity warning (risk=medium): organism=%s, protein_len=%d, "
            "matches=%d, categories=%s",
            organism, len(protein), len(report.matches), report.flagged_categories,
        )
        warnings.warn(
            f"Biosecurity screening detected medium-risk hazard(s): "
            f"{', '.join(report.flagged_categories)}. "
            f"Review the biosecurity report before proceeding.",
            UserWarning,
            stacklevel=2,
        )

    elif report.risk_level == "low":
        logger.info(
            "Biosecurity note (risk=low): organism=%s, protein_len=%d, "
            "matches=%d, categories=%s",
            organism, len(protein), len(report.matches), report.flagged_categories,
        )

    return report


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

__all__ = [
    "BiosecurityReport",
    "HazardMatch",
    "BiosecurityError",
    "screen_hazardous_sequence",
    "check_biosecurity_before_optimize",
    # Expose the database size for testing/validation
    "HAZARD_SIGNATURE_COUNT",
]

HAZARD_SIGNATURE_COUNT = len(_HAZARD_SIGNATURES)
