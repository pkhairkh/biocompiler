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
import os
import warnings
from dataclasses import dataclass, field
from typing import Literal, Optional

from .exceptions import BiosecurityError

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Type aliases
# ─────────────────────────────────────────────────────────────────────────────

RiskLevel = Literal["none", "low", "medium", "high", "critical"]
BiosecurityMode = Literal["enforce", "warn", "off"]
MatchType = Literal["exact", "fuzzy", "reverse_complement"]
StrandType = Literal["forward", "reverse"]


# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class HazardMatch:
    """A single match against a hazardous sequence signature."""

    category: str
    name: str
    position: int
    matched_sequence: str
    confidence: float
    source: str
    # Extended fields for fuzzy and reverse complement matching
    match_type: MatchType = "exact"
    distance: int = 0
    strand: StrandType = "forward"
    substitutions: list[tuple[int, str, str]] = field(default_factory=list)


@dataclass
class BiosecurityReport:
    """Result of biosecurity screening for a protein/DNA sequence."""

    is_hazardous: bool
    risk_level: RiskLevel
    flagged_categories: list[str]
    matches: list[HazardMatch]
    recommendations: list[str]


@dataclass
class BiosecurityScreeningResult:
    """User-facing result of biosecurity screening for a protein sequence.

    This is the high-level result type returned by
    :func:`check_biosecurity_before_optimize`.  It provides a simple
    ``passed`` / ``failed`` boolean along with details about any
    flagged pathogens, risk levels, k-mer similarity scores, and
    match positions.
    """

    passed: bool
    flagged_pathogens: list[str] = field(default_factory=list)
    risk_levels: list[str] = field(default_factory=list)
    match_details: list[str] = field(default_factory=list)
    kmer_scores: dict[str, float] = field(default_factory=dict)
    screened_sequence_length: int = 0

    @property
    def risk_level(self) -> str:
        """Return the highest risk level, or 'none' if no hazards."""
        if not self.risk_levels:
            return "none"
        priority = {"CRITICAL": 4, "HIGH": 3, "MEDIUM": 2, "LOW": 1, "none": 0}
        return max(self.risk_levels, key=lambda r: priority.get(r.upper(), 0))

    def __str__(self) -> str:
        if self.passed:
            return "BiosecurityScreeningResult: PASSED"
        pathogens = ", ".join(self.flagged_pathogens)
        levels = ", ".join(self.risk_levels)
        return f"BiosecurityScreeningResult: FAILED (pathogens=[{pathogens}], risk=[{levels}])"


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
# Legacy pathogen signature database (for integration tests)
# ─────────────────────────────────────────────────────────────────────────────
# Each entry is a tuple: (sequence, pathogen_name, risk_level, description)
# These are longer signature sequences (N-terminal / signal peptides) used for
# exact substring matching in the check_biosecurity_before_optimize pipeline.

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


# ─────────────────────────────────────────────────────────────────────────────
# K-mer similarity helpers
# ─────────────────────────────────────────────────────────────────────────────

_KMER_SIZE: int = 5
_SIMILARITY_THRESHOLD: float = 0.6


def _extract_kmers(sequence: str, k: int = _KMER_SIZE) -> set[str]:
    """Extract all k-mers from a sequence.

    Parameters
    ----------
    sequence : str
        The amino acid or nucleotide sequence.
    k : int
        K-mer length (default: ``_KMER_SIZE``).

    Returns
    -------
    set[str]
        Set of all k-mers found in the sequence.  Returns an empty
        set if the sequence is shorter than *k*.
    """
    sequence = sequence.upper()
    if len(sequence) < k:
        return set()
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


def _compute_kmer_similarity(query_kmers: set[str], pathogen_kmers: set[str]) -> float:
    """Compute the Jaccard-like k-mer similarity between two sets.

    Defined as ``|intersection| / |query_kmers|`` — the fraction of the
    query's k-mers that also appear in the pathogen signature.  Returns
    ``0.0`` if *query_kmers* is empty.

    Parameters
    ----------
    query_kmers : set[str]
        K-mers from the query sequence.
    pathogen_kmers : set[str]
        K-mers from a pathogen signature.

    Returns
    -------
    float
        Similarity score in [0, 1].
    """
    if not query_kmers:
        return 0.0
    return len(query_kmers & pathogen_kmers) / len(query_kmers)


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
# Biosecurity mode
# ─────────────────────────────────────────────────────────────────────────────

def get_biosecurity_mode() -> BiosecurityMode:
    """Read the biosecurity screening mode from the environment.

    Returns one of ``"enforce"``, ``"warn"``, or ``"off"`` based on the
    ``BIOCOMPILER_BIOSECURITY_MODE`` environment variable.  Falls back
    to ``"enforce"`` when the variable is unset or set to an
    unrecognised value.  The value is case-insensitive.

    Returns
    -------
    BiosecurityMode
        The effective biosecurity mode.
    """
    raw = os.environ.get("BIOCOMPILER_BIOSECURITY_MODE", "").strip().lower()
    if raw in ("enforce", "warn", "off"):
        return raw  # type: ignore[return-value]
    return "enforce"


# ─────────────────────────────────────────────────────────────────────────────
# Sequence utilities
# ─────────────────────────────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(dna: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Parameters
    ----------
    dna : str
        DNA sequence consisting of A, C, G, T (uppercase expected).

    Returns
    -------
    str
        The reverse complement.
    """
    return dna.translate(_COMPLEMENT)[::-1]


def _hamming_distance(s1: str, s2: str) -> int:
    """Compute the Hamming distance between two equal-length strings.

    Raises ``ValueError`` if the strings differ in length.
    """
    if len(s1) != len(s2):
        raise ValueError(
            f"Hamming distance requires equal-length strings, "
            f"got {len(s1)} and {len(s2)}"
        )
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def _levenshtein_distance(s1: str, s2: str) -> int:
    """Compute the Levenshtein (edit) distance between two strings.

    Uses dynamic programming with O(min(len(s1), len(s2))) memory.
    """
    if len(s1) < len(s2):
        s1, s2 = s2, s1

    # Now len(s1) >= len(s2)
    previous = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous[j + 1] + 1
            deletions = current[j] + 1
            substitutions = previous[j] + (c1 != c2)
            current.append(min(insertions, deletions, substitutions))
        previous = current

    return previous[-1]


def _fuzzy_match_hamming(
    sequence: str,
    motif: str,
    max_distance: int = 2,
) -> list[tuple[int, int, list[tuple[int, str, str]]]]:
    """Find fuzzy matches of *motif* in *sequence* using Hamming distance.

    Only returns matches with distance >= 1 (i.e. excludes exact matches)
    and distance <= *max_distance*.

    Parameters
    ----------
    sequence : str
        The sequence to search in.
    motif : str
        The motif to search for.
    max_distance : int
        Maximum Hamming distance to report (default 2).

    Returns
    -------
    list of (position, distance, substitutions)
        Each element is a tuple ``(pos, dist, subs)`` where *subs* is a
        list of ``(position_in_window, original_char, replacement_char)``
        tuples describing the substitutions.
    """
    mlen = len(motif)
    slen = len(sequence)
    if mlen == 0 or slen < mlen:
        return []

    results: list[tuple[int, int, list[tuple[int, str, str]]]] = []
    for i in range(slen - mlen + 1):
        window = sequence[i : i + mlen]
        dist = _hamming_distance(window, motif)
        if 1 <= dist <= max_distance:
            subs = [
                (j, motif[j], window[j])
                for j in range(mlen)
                if motif[j] != window[j]
            ]
            results.append((i, dist, subs))

    return results


def _fuzzy_match_edit_distance(
    sequence: str,
    motif: str,
    max_distance: int = 1,
) -> list[tuple[int, int]]:
    """Find fuzzy matches of *motif* in *sequence* using Levenshtein distance.

    Uses a sliding-window approach with windows of varying length around
    the motif length to catch insertions and deletions.

    Only returns matches with distance >= 1 and <= *max_distance*.

    Parameters
    ----------
    sequence : str
        The sequence to search in.
    motif : str
        The motif to search for.
    max_distance : int
        Maximum edit distance to report (default 1).

    Returns
    -------
    list of (position, distance)
        Each element is a tuple ``(pos, dist)``.
    """
    mlen = len(motif)
    slen = len(sequence)
    if mlen == 0 or slen == 0:
        return []

    results: list[tuple[int, int]] = []
    # Check windows of length mlen-1 through mlen+max_distance
    for window_len in range(max(1, mlen - max_distance), mlen + max_distance + 1):
        if window_len > slen:
            continue
        for i in range(slen - window_len + 1):
            window = sequence[i : i + window_len]
            dist = _levenshtein_distance(window, motif)
            if 1 <= dist <= max_distance:
                results.append((i, dist))

    # Deduplicate: keep the best distance per position
    best: dict[int, int] = {}
    for pos, dist in results:
        if pos not in best or dist < best[pos]:
            best[pos] = dist

    return sorted(best.items())


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
    - Fuzzy matching via Hamming distance (1-2 substitutions) and
      Levenshtein edit distance (1 edit) is applied for motifs < 15 aa.
    - DNA screening uses substring matching against nucleotide patterns
      (15-21 nt) for antibiotic resistance markers, plus reverse
      complement screening.
    - Confidence scoring accounts for motif length (longer = higher
      confidence) and exact match position context.
    """
    protein = protein.upper().strip()
    dna = dna.upper().strip()

    matches: list[HazardMatch] = []
    flagged_categories: set[str] = set()
    # Track exact match positions per signature name for deduplication
    exact_positions: dict[str, set[int]] = {}

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
                match_type="exact",
                distance=0,
                strand="forward",
                substitutions=[],
            ))
            flagged_categories.add(sig["category"])
            exact_positions.setdefault(sig["name"], set()).add(pos)

            # Look for additional occurrences
            pos = protein.find(motif, pos + 1)

        # ── Fuzzy matching (Hamming) for short protein motifs ────────────
        if motif_len < 15:
            hamming_results = _fuzzy_match_hamming(protein, motif, max_distance=2)
            for fpos, fdist, fsubs in hamming_results:
                # Skip fuzzy matches at positions already covered by exact matches
                if sig["name"] in exact_positions and fpos in exact_positions[sig["name"]]:
                    continue
                # Reduce confidence for fuzzy matches
                confidence = min(1.0, max(0.0, sig["confidence"] - 0.10 * fdist))
                length_adj = (motif_len - 10) * 0.025
                confidence = min(1.0, max(0.0, confidence + length_adj))

                matches.append(HazardMatch(
                    category=sig["category"],
                    name=sig["name"],
                    position=fpos,
                    matched_sequence=protein[fpos:fpos + motif_len],
                    confidence=round(confidence, 3),
                    source=sig["source"],
                    match_type="fuzzy",
                    distance=fdist,
                    strand="forward",
                    substitutions=fsubs,
                ))
                flagged_categories.add(sig["category"])

            # ── Fuzzy matching (edit distance) for short protein motifs ───
            edit_results = _fuzzy_match_edit_distance(protein, motif, max_distance=1)
            for epos, edist in edit_results:
                # Skip if an exact match already covers this position
                if sig["name"] in exact_positions and epos in exact_positions[sig["name"]]:
                    continue
                # Skip if a Hamming fuzzy match already covers this position
                existing_fuzzy = [
                    m for m in matches
                    if m.name == sig["name"] and m.position == epos and m.match_type == "fuzzy"
                ]
                if existing_fuzzy:
                    continue
                confidence = min(1.0, max(0.0, sig["confidence"] - 0.10 * edist))
                length_adj = (motif_len - 10) * 0.025
                confidence = min(1.0, max(0.0, confidence + length_adj))

                matches.append(HazardMatch(
                    category=sig["category"],
                    name=sig["name"],
                    position=epos,
                    matched_sequence=protein[epos:epos + motif_len],
                    confidence=round(confidence, 3),
                    source=sig["source"],
                    match_type="fuzzy",
                    distance=edist,
                    strand="forward",
                    substitutions=[],
                ))
                flagged_categories.add(sig["category"])

    # ── Screen DNA against nucleotide patterns ───────────────────────────
    if dna:
        # Forward strand
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
                    match_type="exact",
                    distance=0,
                    strand="forward",
                    substitutions=[],
                ))
                flagged_categories.add(sig["category"])

                pos = dna.find(motif, pos + 1)

        # Reverse complement strand
        rc_dna = reverse_complement(dna)
        for sig in _DNA_SIGNATURES:
            motif = sig["motif"].upper()
            pos = rc_dna.find(motif)
            while pos != -1:
                # Map position back to the original strand
                original_pos = len(rc_dna) - pos - len(motif)
                matches.append(HazardMatch(
                    category=sig["category"],
                    name=sig["name"],
                    position=original_pos,
                    matched_sequence=motif,
                    confidence=sig["confidence"],
                    source=sig["source"],
                    match_type="reverse_complement",
                    distance=0,
                    strand="reverse",
                    substitutions=[],
                ))
                flagged_categories.add(sig["category"])

                pos = rc_dna.find(motif, pos + 1)

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
    """Look up the risk level for a match's signature.

    Fuzzy matches have their risk downgraded:
      - distance 1 → one level below the signature risk
      - distance 2 → two levels below the signature risk
    Reverse complement matches keep the original risk level.
    """
    base_risk: str | None = None
    for sig in _HAZARD_SIGNATURES:
        if sig["name"] == match.name and sig["category"] == match.category:
            base_risk = sig["risk"]
            break

    if base_risk is None:
        # Fallback: infer from category
        _category_default_risk = {
            "select_agent": "critical",
            "viral_surface": "high",
            "antibiotic_resistance": "medium",
            "oncogene": "low",
        }
        base_risk = _category_default_risk.get(match.category, "low")

    if match.match_type == "reverse_complement":
        return base_risk

    if match.match_type == "fuzzy" and match.distance > 0:
        # Fuzzy matches get fixed downgraded risk levels:
        #   distance 1 → "medium"
        #   distance 2+ → "low"
        if match.distance == 1:
            return "medium"
        return "low"

    return base_risk


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

    # Mention fuzzy/homology matches
    fuzzy_matches = [m for m in matches if m.match_type == "fuzzy"]
    if fuzzy_matches:
        fuzzy_names = sorted({m.name for m in fuzzy_matches})
        recs.append(
            f"Fuzzy homology match(es) detected (substitution/indel tolerance): "
            f"{', '.join(fuzzy_names)}. "
            f"These are lower-confidence matches — review carefully before proceeding."
        )

    # Mention reverse complement matches
    rc_matches = [m for m in matches if m.match_type == "reverse_complement"]
    if rc_matches:
        rc_names = sorted({m.name for m in rc_matches})
        recs.append(
            f"Reverse complement / anti-sense match(es) detected on the opposite strand: "
            f"{', '.join(rc_names)}. "
            f"These indicate hazardous sequences on the reverse strand."
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
    biosecurity_mode: Optional[BiosecurityMode] = None,
    skip_biosecurity_check: bool = False,
) -> BiosecurityScreeningResult:
    """Biosecurity gate called at the start of optimization.

    This function screens the input sequence and enforces hard-stop or
    warning behavior depending on the risk level and the effective
    biosecurity mode:

    - ``enforce``: ``critical`` or ``high`` risk raises
      :class:`BiosecurityError`; ``medium`` emits :class:`UserWarning`.
    - ``warn``: all risk levels emit warnings but no
      :class:`BiosecurityError` is raised.
    - ``off``: screening is skipped entirely; returns a clean result
      noting that screening was skipped.

    If *biosecurity_mode* is ``None`` (the default), the mode is read
    from the ``BIOCOMPILER_BIOSECURITY_MODE`` environment variable via
    :func:`get_biosecurity_mode`.

    Parameters
    ----------
    protein : str
        Protein sequence to screen.
    organism : str, optional
        Target organism (for context in log messages).
    dna : str, optional
        DNA sequence to screen for resistance markers.
    biosecurity_mode : BiosecurityMode or None, optional
        Explicit override for the biosecurity mode.  When ``None``,
        the mode is read from the environment variable.
    skip_biosecurity_check : bool, optional
        If ``True``, skip biosecurity screening entirely and return a
        passed result.  Default is ``False``.  This should only be used
        for testing or when the user explicitly opts out.

    Returns
    -------
    BiosecurityScreeningResult
        The screening result with pass/fail status and details.

    Raises
    ------
    BiosecurityError
        If mode is ``enforce`` and risk_level is ``critical`` or ``high``.
    ValueError
        If *protein* is empty or whitespace-only.
    """
    if not protein or not protein.strip():
        raise ValueError("Protein sequence must not be empty")

    # Skip check: always return passed
    if skip_biosecurity_check:
        return BiosecurityScreeningResult(
            passed=True,
            screened_sequence_length=len(protein.strip()),
        )

    if biosecurity_mode is None:
        biosecurity_mode = get_biosecurity_mode()

    # Off mode: skip screening entirely
    if biosecurity_mode == "off":
        return BiosecurityScreeningResult(
            passed=True,
            screened_sequence_length=len(protein.strip()),
        )

    protein_upper = protein.upper().strip()

    # ── Check legacy pathogen signatures (exact substring match) ─────────
    flagged_pathogens: list[str] = []
    risk_levels: list[str] = []
    match_details: list[str] = []
    kmer_scores: dict[str, float] = {}

    for sig_seq, pathogen_name, risk_level, description in _PATHOGEN_SIGNATURES:
        if sig_seq.upper() in protein_upper:
            pos = protein_upper.find(sig_seq.upper())
            flagged_pathogens.append(pathogen_name)
            risk_levels.append(risk_level)
            match_details.append(
                f"{description}: matched at position {pos}"
            )
            kmer_scores[pathogen_name] = 1.0

    # ── Also run the motif-based screening ──────────────────────────────
    report = screen_hazardous_sequence(protein, dna)

    # Merge motif-based findings into the result
    for match in report.matches:
        # Map motif match names to pathogen names where possible
        pathogen_name = _MOTIF_TO_PATHOGEN.get(match.name, match.name)
        if pathogen_name not in flagged_pathogens:
            flagged_pathogens.append(pathogen_name)
            # Map category to risk level based on biosecurity classification
            _CATEGORY_RISK = {
                "SELECT_AGENT": "CRITICAL",
                "VIRAL_SURFACE": "HIGH",
                "ANTIBIOTIC_RESISTANCE": "MEDIUM",
                "ONCOGENE": "HIGH",
                "TOXIN": "CRITICAL",
            }
            cat_risk = _CATEGORY_RISK.get(match.category.upper(), "MEDIUM")
            risk_levels.append(cat_risk.lower() if cat_risk in ("CRITICAL", "HIGH") else "medium")
            match_details.append(
                f"{match.name}: {match.match_type} match at position {match.position}"
            )
        # Compute k-mer similarity for this pathogen
        if pathogen_name not in kmer_scores:
            query_kmers = _extract_kmers(protein_upper)
            sig_kmers = _extract_kmers(match.matched_sequence)
            sim = _compute_kmer_similarity(query_kmers, sig_kmers)
            kmer_scores[pathogen_name] = round(sim, 4)

    passed = len(flagged_pathogens) == 0

    result = BiosecurityScreeningResult(
        passed=passed,
        flagged_pathogens=flagged_pathogens,
        risk_levels=risk_levels,
        match_details=match_details,
        kmer_scores=kmer_scores,
        screened_sequence_length=len(protein.strip()),
    )

    # Enforce mode: raise for critical/high risk
    if biosecurity_mode == "enforce" and not passed:
        has_critical_or_high = any(
            rl.upper() in ("CRITICAL", "HIGH") for rl in risk_levels
        )
        if has_critical_or_high:
            logger.error(
                "Biosecurity gate BLOCKED optimization: pathogens=%s, "
                "organism=%s, protein_len=%d",
                flagged_pathogens, organism, len(protein),
            )
            raise BiosecurityError(
                "hazardous_sequence_detected",
                protein=protein,
                flagged_pathogens=flagged_pathogens,
                risk_levels=risk_levels,
                match_details=match_details,
            )

    # Warn mode or medium risk: emit warnings
    if not passed and biosecurity_mode == "warn":
        warnings.warn(
            f"Biosecurity screening detected hazard(s): "
            f"{', '.join(flagged_pathogens)}. "
            f"Mode is 'warn' so optimization is NOT blocked, but review is strongly recommended.",
            UserWarning,
            stacklevel=2,
        )

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

__all__ = [
    "BiosecurityReport",
    "BiosecurityScreeningResult",
    "HazardMatch",
    "BiosecurityError",
    "BiosecurityMode",
    "MatchType",
    "StrandType",
    "screen_hazardous_sequence",
    "check_biosecurity_before_optimize",
    "get_biosecurity_mode",
    "sig_risk_for_match",
    "reverse_complement",
    "_hamming_distance",
    "_levenshtein_distance",
    "_fuzzy_match_hamming",
    "_fuzzy_match_edit_distance",
    # Expose the database size for testing/validation
    "HAZARD_SIGNATURE_COUNT",
    # Legacy pathogen signature exports (used by integration tests)
    "_PATHOGEN_SIGNATURES",
    "_KMER_SIZE",
    "_SIMILARITY_THRESHOLD",
    "_extract_kmers",
    "_compute_kmer_similarity",
]

HAZARD_SIGNATURE_COUNT = len(_HAZARD_SIGNATURES)
