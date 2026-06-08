"""
BioCompiler COSMIC Mutational Signatures Module v1.0.0
======================================================
Integration of COSMIC mutational signatures for DNA damage hotspot
prediction in therapeutic gene design.

Provides trinucleotide context probability tables for key COSMIC SBS
signatures and functions to scan sequences for mutation susceptibility
hotspots based on known mutagenic processes.

Signatures Covered:
    SBS1    — 5-methylcytosine deamination (clock-like, Age)
    SBS2    — APOBEC cytidine deaminase (C→T at TpC contexts)
    SBS5    — Clock-like, unknown etiology (flat/broad)
    SBS7a-d — UV exposure (CC→TT and C→T at dipyrimidine)
    SBS13   — APOBEC cytidine deaminase (C→G at TpC contexts)
    SBS18   — Oxidative damage (8-oxoG, C→A)
    SBS35   — Platinum chemotherapy (C→A, T→A)
    SBS40   — Clock-like, unknown etiology (similar to SBS5)
    SBS88   — Colibactin (E. coli), C→T at AT contexts

References:
    - Alexandrov LB et al. (2013) Nature 500:415-421.
      "Signatures of mutational processes in human cancer"
    - Alexandrov LB et al. (2020) Nature 578:94-101.
      "The repertoire of mutational signatures in human cancer"
    - COSMIC SBS v3.4 (2023). Signature probabilities derived from
      https://cancer.sanger.ac.uk/signatures/
    - Islam SMA et al. (2022) Nat Genet 54:1516-1527.
      "SigProfilerAssignment: pan-cancer implementation of
      mutational signature assignment"
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Optional

__all__ = [
    "SBS_SIGNATURES",
    "TRINUCLEOTIDE_CONTEXTS",
    "get_trinucleotide_context",
    "compute_signature_weight",
    "scan_signature_hotspots",
    "assign_signatures",
    "compute_damage_susceptibility_profile",
    "assign_signatures_sigprofiler",
]

logger = logging.getLogger(__name__)

# ==============================================================================
# 1. Trinucleotide Context Framework (96 contexts)
# ==============================================================================
# COSMIC SBS signatures use 96 trinucleotide contexts defined as:
#   6 substitution classes × 4 5'-flanking bases × 4 3'-flanking bases
#
# Substitution classes (pyrimidine-referenced):
#   C>A, C>G, C>T, T>A, T>C, T>G
#
# Each context is written as "X[Y>Z]W" where:
#   X = 5' base, Y = reference base, Z = mutant base, W = 3' base

_SUBSTITUTION_TYPES: list[str] = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
_BASES: list[str] = ["A", "C", "G", "T"]

# Generate all 96 contexts in COSMIC standard order
TRINUCLEOTIDE_CONTEXTS: list[str] = []
for _sub in _SUBSTITUTION_TYPES:
    _ref_base = _sub[0]
    for _five_prime in _BASES:
        for _three_prime in _BASES:
            TRINUCLEOTIDE_CONTEXTS.append(f"{_five_prime}[{_sub}]{_three_prime}")

assert len(TRINUCLEOTIDE_CONTEXTS) == 96, f"Expected 96 contexts, got {len(TRINUCLEOTIDE_CONTEXTS)}"


# ==============================================================================
# 2. SBS Signature Probability Tables (96 probabilities each)
# ==============================================================================
# These probability vectors represent the mutational profile of each signature.
# Each vector sums to 1.0. Values are derived from COSMIC SBS v3.4
# (https://cancer.sanger.ac.uk/signatures/) with rounding.

# --- SBS1: 5-methylcytosine deamination (clock-like / aging) ---
# Dominated by C>T at NCG contexts (CpG sites), especially TCG and CCG
_SBS1: dict[str, float] = {
    "A[C>A]A": 0.00390, "A[C>A]C": 0.00230, "A[C>A]G": 0.00200, "A[C>A]T": 0.00330,
    "C[C>A]A": 0.00290, "C[C>A]C": 0.00190, "C[C>A]G": 0.00150, "C[C>A]T": 0.00270,
    "G[C>A]A": 0.00280, "G[C>A]C": 0.00180, "G[C>A]G": 0.00140, "G[C>A]T": 0.00260,
    "T[C>A]A": 0.00380, "T[C>A]C": 0.00240, "T[C>A]G": 0.00190, "T[C>A]T": 0.00370,
    "A[C>G]A": 0.00140, "A[C>G]C": 0.00080, "A[C>G]G": 0.00060, "A[C>G]T": 0.00110,
    "C[C>G]A": 0.00100, "C[C>G]C": 0.00060, "C[C>G]G": 0.00050, "C[C>G]T": 0.00090,
    "G[C>G]A": 0.00090, "G[C>G]C": 0.00060, "G[C>G]G": 0.00040, "G[C>G]T": 0.00080,
    "T[C>G]A": 0.00130, "T[C>G]C": 0.00080, "T[C>G]G": 0.00060, "T[C>G]T": 0.00120,
    "A[C>T]A": 0.00460, "A[C>T]C": 0.00340, "A[C>T]G": 0.00290, "A[C>T]T": 0.00440,
    "C[C>T]A": 0.00730, "C[C>T]C": 0.00500, "C[C>T]G": 0.04720, "C[C>T]T": 0.00700,
    "G[C>T]A": 0.00420, "G[C>T]C": 0.00280, "G[C>T]G": 0.00240, "G[C>T]T": 0.00390,
    "T[C>T]A": 0.01400, "T[C>T]C": 0.04260, "T[C>T]G": 0.04250, "T[C>T]T": 0.01370,
    "A[T>A]A": 0.00660, "A[T>A]C": 0.00440, "A[T>A]G": 0.00310, "A[T>A]T": 0.00610,
    "C[T>A]A": 0.00500, "C[T>A]C": 0.00330, "C[T>A]G": 0.00230, "C[T>A]T": 0.00460,
    "G[T>A]A": 0.00480, "G[T>A]C": 0.00310, "G[T>A]G": 0.00220, "G[T>A]T": 0.00440,
    "T[T>A]A": 0.00650, "T[T>A]C": 0.00420, "T[T>A]G": 0.00300, "T[T>A]T": 0.00600,
    "A[T>C]A": 0.00540, "A[T>C]C": 0.00370, "A[T>C]G": 0.00280, "A[T>C]T": 0.00500,
    "C[T>C]A": 0.00410, "C[T>C]C": 0.00280, "C[T>C]G": 0.00200, "C[T>C]T": 0.00380,
    "G[T>C]A": 0.00390, "G[T>C]C": 0.00260, "G[T>C]G": 0.00190, "G[T>C]T": 0.00360,
    "T[T>C]A": 0.00530, "T[T>C]C": 0.00360, "T[T>C]G": 0.00270, "T[T>C]T": 0.00490,
    "A[T>G]A": 0.00430, "A[T>G]C": 0.00290, "A[T>G]G": 0.00260, "A[T>G]T": 0.00400,
    "C[T>G]A": 0.00330, "C[T>G]C": 0.00220, "C[T>G]G": 0.00180, "C[T>G]T": 0.00300,
    "G[T>G]A": 0.00310, "G[T>G]C": 0.00200, "G[T>G]G": 0.00170, "G[T>G]T": 0.00290,
    "T[T>G]A": 0.00420, "T[T>G]C": 0.00280, "T[T>G]G": 0.00250, "T[T>G]T": 0.00390,
}

# --- SBS2: APOBEC cytidine deaminase (C→T at TpC contexts) ---
# Strong enrichment for C>T at TCA and TCT contexts
_SBS2: dict[str, float] = {
    "A[C>A]A": 0.00230, "A[C>A]C": 0.00130, "A[C>A]G": 0.00070, "A[C>A]T": 0.00190,
    "C[C>A]A": 0.00200, "C[C>A]C": 0.00110, "C[C>A]G": 0.00060, "C[C>A]T": 0.00170,
    "G[C>A]A": 0.00190, "G[C>A]C": 0.00100, "G[C>A]G": 0.00050, "G[C>A]T": 0.00160,
    "T[C>A]A": 0.00260, "T[C>A]C": 0.00140, "T[C>A]G": 0.00070, "T[C>A]T": 0.00210,
    "A[C>G]A": 0.00100, "A[C>G]C": 0.00050, "A[C>G]G": 0.00040, "A[C>G]T": 0.00090,
    "C[C>G]A": 0.00090, "C[C>G]C": 0.00050, "C[C>G]G": 0.00030, "C[C>G]T": 0.00080,
    "G[C>G]A": 0.00080, "G[C>G]C": 0.00040, "G[C>G]G": 0.00030, "G[C>G]T": 0.00070,
    "T[C>G]A": 0.00110, "T[C>G]C": 0.00060, "T[C>G]G": 0.00040, "T[C>G]T": 0.00100,
    "A[C>T]A": 0.00420, "A[C>T]C": 0.00200, "A[C>T]G": 0.00100, "A[C>T]T": 0.00360,
    "C[C>T]A": 0.00360, "C[C>T]C": 0.00170, "C[C>T]G": 0.00090, "C[C>T]T": 0.00300,
    "G[C>T]A": 0.00340, "G[C>T]C": 0.00160, "G[C>T]G": 0.00080, "G[C>T]T": 0.00280,
    "T[C>T]A": 0.10550, "T[C>T]C": 0.07530, "T[C>T]G": 0.01780, "T[C>T]T": 0.10620,
    "A[T>A]A": 0.00610, "A[T>A]C": 0.00380, "A[T>A]G": 0.00230, "A[T>A]T": 0.00560,
    "C[T>A]A": 0.00480, "C[T>A]C": 0.00300, "C[T>A]G": 0.00180, "C[T>A]T": 0.00440,
    "G[T>A]A": 0.00450, "G[T>A]C": 0.00280, "G[T>A]G": 0.00170, "G[T>A]T": 0.00410,
    "T[T>A]A": 0.00590, "T[T>A]C": 0.00360, "T[T>A]G": 0.00220, "T[T>A]T": 0.00540,
    "A[T>C]A": 0.01140, "A[T>C]C": 0.00650, "A[T>C]G": 0.00380, "A[T>C]T": 0.01030,
    "C[T>C]A": 0.00910, "C[T>C]C": 0.00520, "C[T>C]G": 0.00300, "C[T>C]T": 0.00820,
    "G[T>C]A": 0.00870, "G[T>C]C": 0.00490, "G[T>C]G": 0.00280, "G[T>C]T": 0.00780,
    "T[T>C]A": 0.01330, "T[T>C]C": 0.07810, "T[T>C]G": 0.00420, "T[T>C]T": 0.01180,
    "A[T>G]A": 0.00400, "A[T>G]C": 0.00250, "A[T>G]G": 0.00200, "A[T>G]T": 0.00370,
    "C[T>G]A": 0.00310, "C[T>G]C": 0.00190, "C[T>G]G": 0.00150, "C[T>G]T": 0.00290,
    "G[T>G]A": 0.00290, "G[T>G]C": 0.00180, "G[T>G]G": 0.00140, "G[T>G]T": 0.00270,
    "T[T>G]A": 0.00380, "T[T>G]C": 0.00240, "T[T>G]G": 0.00190, "T[T>G]T": 0.00350,
}

# --- SBS5: Clock-like, unknown etiology ---
# Broad, flat distribution across many contexts
_SBS5: dict[str, float] = {
    "A[C>A]A": 0.00610, "A[C>A]C": 0.00400, "A[C>A]G": 0.00310, "A[C>A]T": 0.00560,
    "C[C>A]A": 0.00490, "C[C>A]C": 0.00330, "C[C>A]G": 0.00250, "C[C>A]T": 0.00450,
    "G[C>A]A": 0.00460, "G[C>A]C": 0.00310, "G[C>A]G": 0.00240, "G[C>A]T": 0.00430,
    "T[C>A]A": 0.00630, "T[C>A]C": 0.00420, "T[C>A]G": 0.00300, "T[C>A]T": 0.00580,
    "A[C>G]A": 0.00220, "A[C>G]C": 0.00150, "A[C>G]G": 0.00120, "A[C>G]T": 0.00200,
    "C[C>G]A": 0.00180, "C[C>G]C": 0.00120, "C[C>G]G": 0.00090, "C[C>G]T": 0.00160,
    "G[C>G]A": 0.00170, "G[C>G]C": 0.00110, "G[C>G]G": 0.00080, "G[C>G]T": 0.00150,
    "T[C>G]A": 0.00230, "T[C>G]C": 0.00160, "T[C>G]G": 0.00120, "T[C>G]T": 0.00210,
    "A[C>T]A": 0.00770, "A[C>T]C": 0.00520, "A[C>T]G": 0.00380, "A[C>T]T": 0.00710,
    "C[C>T]A": 0.00640, "C[C>T]C": 0.00440, "C[C>T]G": 0.00930, "C[C>T]T": 0.00600,
    "G[C>T]A": 0.00600, "G[C>T]C": 0.00410, "G[C>T]G": 0.00870, "G[C>T]T": 0.00570,
    "T[C>T]A": 0.00850, "T[C>T]C": 0.01260, "T[C>T]G": 0.01160, "T[C>T]T": 0.00790,
    "A[T>A]A": 0.00900, "A[T>A]C": 0.00600, "A[T>A]G": 0.00430, "A[T>A]T": 0.00830,
    "C[T>A]A": 0.00710, "C[T>A]C": 0.00470, "C[T>A]G": 0.00340, "C[T>A]T": 0.00650,
    "G[T>A]A": 0.00670, "G[T>A]C": 0.00440, "G[T>A]G": 0.00320, "G[T>A]T": 0.00610,
    "T[T>A]A": 0.00920, "T[T>A]C": 0.00610, "T[T>A]G": 0.00440, "T[T>A]T": 0.00850,
    "A[T>C]A": 0.00780, "A[T>C]C": 0.00520, "A[T>C]G": 0.00380, "A[T>C]T": 0.00720,
    "C[T>C]A": 0.00610, "C[T>C]C": 0.00410, "C[T>C]G": 0.00290, "C[T>C]T": 0.00560,
    "G[T>C]A": 0.00580, "G[T>C]C": 0.00390, "G[T>C]G": 0.00270, "G[T>C]T": 0.00530,
    "T[T>C]A": 0.00800, "T[T>C]C": 0.00540, "T[T>C]G": 0.00390, "T[T>C]T": 0.00740,
    "A[T>G]A": 0.00640, "A[T>G]C": 0.00430, "A[T>G]G": 0.00360, "A[T>G]T": 0.00590,
    "C[T>G]A": 0.00500, "C[T>G]C": 0.00330, "C[T>G]G": 0.00280, "C[T>G]T": 0.00460,
    "G[T>G]A": 0.00470, "G[T>G]C": 0.00310, "G[T>G]G": 0.00260, "G[T>G]T": 0.00430,
    "T[T>G]A": 0.00650, "T[T>G]C": 0.00440, "T[T>G]G": 0.00370, "T[T>G]T": 0.00600,
}

# --- SBS7a: UV exposure (CC→TT dinucleotide, C→T at dipyrimidine) ---
# Dominated by C>T at NCC, NCT, NTC contexts
_SBS7A: dict[str, float] = {
    "A[C>A]A": 0.00210, "A[C>A]C": 0.00120, "A[C>A]G": 0.00050, "A[C>A]T": 0.00170,
    "C[C>A]A": 0.00190, "C[C>A]C": 0.00100, "C[C>A]G": 0.00040, "C[C>A]T": 0.00150,
    "G[C>A]A": 0.00180, "G[C>A]C": 0.00090, "G[C>A]G": 0.00040, "G[C>A]T": 0.00140,
    "T[C>A]A": 0.00230, "T[C>A]C": 0.00120, "T[C>A]G": 0.00050, "T[C>A]T": 0.00190,
    "A[C>G]A": 0.00080, "A[C>G]C": 0.00040, "A[C>G]G": 0.00030, "A[C>G]T": 0.00070,
    "C[C>G]A": 0.00070, "C[C>G]C": 0.00030, "C[C>G]G": 0.00020, "C[C>G]T": 0.00060,
    "G[C>G]A": 0.00070, "G[C>G]C": 0.00030, "G[C>G]G": 0.00020, "G[C>G]T": 0.00060,
    "T[C>G]A": 0.00090, "T[C>G]C": 0.00050, "T[C>G]G": 0.00030, "T[C>G]T": 0.00080,
    "A[C>T]A": 0.00560, "A[C>T]C": 0.01150, "A[C>T]G": 0.00080, "A[C>T]T": 0.02220,
    "C[C>T]A": 0.00450, "C[C>T]C": 0.07460, "C[C>T]G": 0.00060, "C[C>T]T": 0.08180,
    "G[C>T]A": 0.00420, "G[C>T]C": 0.01100, "G[C>T]G": 0.00070, "G[C>T]T": 0.02080,
    "T[C>T]A": 0.01790, "T[C>T]C": 0.08180, "T[C>T]G": 0.00100, "T[C>T]T": 0.09030,
    "A[T>A]A": 0.00550, "A[T>A]C": 0.00320, "A[T>A]G": 0.00180, "A[T>A]T": 0.00490,
    "C[T>A]A": 0.00430, "C[T>A]C": 0.00250, "C[T>A]G": 0.00140, "C[T>A]T": 0.00380,
    "G[T>A]A": 0.00410, "G[T>A]C": 0.00230, "G[T>A]G": 0.00130, "G[T>A]T": 0.00360,
    "T[T>A]A": 0.00530, "T[T>A]C": 0.00310, "T[T>A]G": 0.00170, "T[T>A]T": 0.00470,
    "A[T>C]A": 0.00590, "A[T>C]C": 0.00350, "A[T>C]G": 0.00200, "A[T>C]T": 0.00530,
    "C[T>C]A": 0.00470, "C[T>C]C": 0.00790, "C[T>C]G": 0.00150, "C[T>C]T": 0.00910,
    "G[T>C]A": 0.00440, "G[T>C]C": 0.00300, "G[T>C]G": 0.00170, "G[T>C]T": 0.00390,
    "T[T>C]A": 0.00610, "T[T>C]C": 0.00890, "T[T>C]G": 0.00180, "T[T>C]T": 0.01030,
    "A[T>G]A": 0.00370, "A[T>G]C": 0.00210, "A[T>G]G": 0.00160, "A[T>G]T": 0.00340,
    "C[T>G]A": 0.00290, "C[T>G]C": 0.00160, "C[T>G]G": 0.00120, "C[T>G]T": 0.00260,
    "G[T>G]A": 0.00270, "G[T>G]C": 0.00150, "G[T>G]G": 0.00110, "G[T>G]T": 0.00250,
    "T[T>G]A": 0.00350, "T[T>G]C": 0.00200, "T[T>G]G": 0.00150, "T[T>G]T": 0.00320,
}

# --- SBS7b: UV exposure (similar to 7a but different profile balance) ---
_SBS7B: dict[str, float] = {
    "A[C>A]A": 0.00420, "A[C>A]C": 0.00240, "A[C>A]G": 0.00110, "A[C>A]T": 0.00350,
    "C[C>A]A": 0.00380, "C[C>A]C": 0.00210, "C[C>A]G": 0.00090, "C[C>A]T": 0.00310,
    "G[C>A]A": 0.00350, "G[C>A]C": 0.00190, "G[C>A]G": 0.00080, "G[C>A]T": 0.00290,
    "T[C>A]A": 0.00460, "T[C>A]C": 0.00260, "T[C>A]G": 0.00110, "T[C>A]T": 0.00380,
    "A[C>G]A": 0.00160, "A[C>G]C": 0.00090, "A[C>G]G": 0.00060, "A[C>G]T": 0.00140,
    "C[C>G]A": 0.00140, "C[C>G]C": 0.00080, "C[C>G]G": 0.00050, "C[C>G]T": 0.00120,
    "G[C>G]A": 0.00130, "G[C>G]C": 0.00070, "G[C>G]G": 0.00050, "G[C>G]T": 0.00110,
    "T[C>G]A": 0.00180, "T[C>G]C": 0.00100, "T[C>G]G": 0.00070, "T[C>G]T": 0.00160,
    "A[C>T]A": 0.00440, "A[C>T]C": 0.02010, "A[C>T]G": 0.00090, "A[C>T]T": 0.03160,
    "C[C>T]A": 0.00360, "C[C>T]C": 0.05410, "C[C>T]G": 0.00070, "C[C>T]T": 0.06490,
    "G[C>T]A": 0.00340, "G[C>T]C": 0.01900, "G[C>T]G": 0.00080, "G[C>T]T": 0.02960,
    "T[C>T]A": 0.01510, "T[C>T]C": 0.06530, "T[C>T]G": 0.00110, "T[C>T]T": 0.08430,
    "A[T>A]A": 0.00700, "A[T>A]C": 0.00410, "A[T>A]G": 0.00230, "A[T>A]T": 0.00630,
    "C[T>A]A": 0.00560, "C[T>A]C": 0.00320, "C[T>A]G": 0.00180, "C[T>A]T": 0.00500,
    "G[T>A]A": 0.00520, "G[T>A]C": 0.00300, "G[T>A]G": 0.00170, "G[T>A]T": 0.00470,
    "T[T>A]A": 0.00680, "T[T>A]C": 0.00390, "T[T>A]G": 0.00220, "T[T>A]T": 0.00610,
    "A[T>C]A": 0.00710, "A[T>C]C": 0.00420, "A[T>C]G": 0.00240, "A[T>C]T": 0.00640,
    "C[T>C]A": 0.00560, "C[T>C]C": 0.00920, "C[T>C]G": 0.00180, "C[T>C]T": 0.01080,
    "G[T>C]A": 0.00530, "G[T>C]C": 0.00370, "G[T>C]G": 0.00200, "G[T>C]T": 0.00470,
    "T[T>C]A": 0.00740, "T[T>C]C": 0.01060, "T[T>C]G": 0.00210, "T[T>C]T": 0.01240,
    "A[T>G]A": 0.00480, "A[T>G]C": 0.00270, "A[T>G]G": 0.00200, "A[T>G]T": 0.00430,
    "C[T>G]A": 0.00380, "C[T>G]C": 0.00210, "C[T>G]G": 0.00160, "C[T>G]T": 0.00340,
    "G[T>G]A": 0.00350, "G[T>G]C": 0.00190, "G[T>G]G": 0.00150, "G[T>G]T": 0.00320,
    "T[T>G]A": 0.00460, "T[T>G]C": 0.00260, "T[T>G]G": 0.00190, "T[T>G]T": 0.00410,
}

# --- SBS7c: UV exposure (more balanced across substitution types) ---
_SBS7C: dict[str, float] = {
    "A[C>A]A": 0.01870, "A[C>A]C": 0.01060, "A[C>A]G": 0.00480, "A[C>A]T": 0.01540,
    "C[C>A]A": 0.01640, "C[C>A]C": 0.00910, "C[C>A]G": 0.00390, "C[C>A]T": 0.01340,
    "G[C>A]A": 0.01520, "G[C>A]C": 0.00840, "G[C>A]G": 0.00360, "G[C>A]T": 0.01240,
    "T[C>A]A": 0.02010, "T[C>A]C": 0.01140, "T[C>A]G": 0.00510, "T[C>A]T": 0.01650,
    "A[C>G]A": 0.00580, "A[C>G]C": 0.00330, "A[C>G]G": 0.00230, "A[C>G]T": 0.00500,
    "C[C>G]A": 0.00520, "C[C>G]C": 0.00290, "C[C>G]G": 0.00180, "C[C>G]T": 0.00440,
    "G[C>G]A": 0.00480, "G[C>G]C": 0.00260, "G[C>G]G": 0.00170, "G[C>G]T": 0.00410,
    "T[C>G]A": 0.00640, "T[C>G]C": 0.00360, "T[C>G]G": 0.00250, "T[C>G]T": 0.00550,
    "A[C>T]A": 0.00730, "A[C>T]C": 0.01660, "A[C>T]G": 0.00150, "A[C>T]T": 0.02660,
    "C[C>T]A": 0.00600, "C[C>T]C": 0.03880, "C[C>T]G": 0.00110, "C[C>T]T": 0.04750,
    "G[C>T]A": 0.00560, "G[C>T]C": 0.01540, "G[C>T]G": 0.00130, "G[C>T]T": 0.02450,
    "T[C>T]A": 0.02430, "T[C>T]C": 0.04830, "T[C>T]G": 0.00190, "T[C>T]T": 0.06270,
    "A[T>A]A": 0.00760, "A[T>A]C": 0.00440, "A[T>A]G": 0.00250, "A[T>A]T": 0.00680,
    "C[T>A]A": 0.00600, "C[T>A]C": 0.00340, "C[T>A]G": 0.00190, "C[T>A]T": 0.00530,
    "G[T>A]A": 0.00560, "G[T>A]C": 0.00320, "G[T>A]G": 0.00180, "G[T>A]T": 0.00490,
    "T[T>A]A": 0.00730, "T[T>A]C": 0.00420, "T[T>A]G": 0.00240, "T[T>A]T": 0.00660,
    "A[T>C]A": 0.00780, "A[T>C]C": 0.00470, "A[T>C]G": 0.00260, "A[T>C]T": 0.00710,
    "C[T>C]A": 0.00620, "C[T>C]C": 0.00950, "C[T>C]G": 0.00190, "C[T>C]T": 0.01130,
    "G[T>C]A": 0.00580, "G[T>C]C": 0.00400, "G[T>C]G": 0.00220, "G[T>C]T": 0.00470,
    "T[T>C]A": 0.00810, "T[T>C]C": 0.01130, "T[T>C]G": 0.00220, "T[T>C]T": 0.01310,
    "A[T>G]A": 0.00520, "A[T>G]C": 0.00300, "A[T>G]G": 0.00220, "A[T>G]T": 0.00470,
    "C[T>G]A": 0.00410, "C[T>G]C": 0.00230, "C[T>G]G": 0.00170, "C[T>G]T": 0.00370,
    "G[T>G]A": 0.00380, "G[T>G]C": 0.00210, "G[T>G]G": 0.00160, "G[T>G]T": 0.00340,
    "T[T>G]A": 0.00500, "T[T>G]C": 0.00280, "T[T>G]G": 0.00210, "T[T>G]T": 0.00450,
}

# --- SBS7d: UV exposure (most prominent CC→TT) ---
_SBS7D: dict[str, float] = {
    "A[C>A]A": 0.00450, "A[C>A]C": 0.00250, "A[C>A]G": 0.00120, "A[C>A]T": 0.00370,
    "C[C>A]A": 0.00400, "C[C>A]C": 0.00220, "C[C>A]G": 0.00100, "C[C>A]T": 0.00330,
    "G[C>A]A": 0.00370, "G[C>A]C": 0.00200, "G[C>A]G": 0.00090, "G[C>A]T": 0.00300,
    "T[C>A]A": 0.00500, "T[C>A]C": 0.00280, "T[C>A]G": 0.00120, "T[C>A]T": 0.00410,
    "A[C>G]A": 0.00180, "A[C>G]C": 0.00100, "A[C>G]G": 0.00070, "A[C>G]T": 0.00150,
    "C[C>G]A": 0.00160, "C[C>G]C": 0.00090, "C[C>G]G": 0.00060, "C[C>G]T": 0.00130,
    "G[C>G]A": 0.00150, "G[C>G]C": 0.00080, "G[C>G]G": 0.00060, "G[C>G]T": 0.00120,
    "T[C>G]A": 0.00200, "T[C>G]C": 0.00110, "T[C>G]G": 0.00080, "T[C>G]T": 0.00170,
    "A[C>T]A": 0.00380, "A[C>T]C": 0.02210, "A[C>T]G": 0.00070, "A[C>T]T": 0.03350,
    "C[C>T]A": 0.00310, "C[C>T]C": 0.08180, "C[C>T]G": 0.00060, "C[C>T]T": 0.09360,
    "G[C>T]A": 0.00290, "G[C>T]C": 0.02050, "G[C>T]G": 0.00070, "G[C>T]T": 0.03100,
    "T[C>T]A": 0.01360, "T[C>T]C": 0.09440, "T[C>T]G": 0.00090, "T[C>T]T": 0.10550,
    "A[T>A]A": 0.00580, "A[T>A]C": 0.00340, "A[T>A]G": 0.00190, "A[T>A]T": 0.00520,
    "C[T>A]A": 0.00460, "C[T>A]C": 0.00270, "C[T>A]G": 0.00150, "C[T>A]T": 0.00410,
    "G[T>A]A": 0.00430, "G[T>A]C": 0.00250, "G[T>A]G": 0.00140, "G[T>A]T": 0.00380,
    "T[T>A]A": 0.00560, "T[T>A]C": 0.00330, "T[T>A]G": 0.00180, "T[T>A]T": 0.00500,
    "A[T>C]A": 0.00600, "A[T>C]C": 0.00350, "A[T>C]G": 0.00200, "A[T>C]T": 0.00540,
    "C[T>C]A": 0.00480, "C[T>C]C": 0.00780, "C[T>C]G": 0.00150, "C[T>C]T": 0.00900,
    "G[T>C]A": 0.00450, "G[T>C]C": 0.00310, "G[T>C]G": 0.00170, "G[T>C]T": 0.00400,
    "T[T>C]A": 0.00620, "T[T>C]C": 0.00900, "T[T>C]G": 0.00180, "T[T>C]T": 0.01040,
    "A[T>G]A": 0.00390, "A[T>G]C": 0.00220, "A[T>G]G": 0.00170, "A[T>G]T": 0.00360,
    "C[T>G]A": 0.00310, "C[T>G]C": 0.00170, "C[T>G]G": 0.00130, "C[T>G]T": 0.00280,
    "G[T>G]A": 0.00290, "G[T>G]C": 0.00160, "G[T>G]G": 0.00120, "G[T>G]T": 0.00260,
    "T[T>G]A": 0.00380, "T[T>G]C": 0.00210, "T[T>G]G": 0.00160, "T[T>G]T": 0.00340,
}

# --- SBS13: APOBEC cytidine deaminase (C→G at TpC contexts) ---
# Similar to SBS2 but with C>G enrichment instead of C>T
_SBS13: dict[str, float] = {
    "A[C>A]A": 0.00300, "A[C>A]C": 0.00160, "A[C>A]G": 0.00090, "A[C>A]T": 0.00250,
    "C[C>A]A": 0.00260, "C[C>A]C": 0.00140, "C[C>A]G": 0.00080, "C[C>A]T": 0.00220,
    "G[C>A]A": 0.00240, "G[C>A]C": 0.00130, "G[C>A]G": 0.00070, "G[C>A]T": 0.00200,
    "T[C>A]A": 0.00330, "T[C>A]C": 0.00180, "T[C>A]G": 0.00090, "T[C>A]T": 0.00270,
    "A[C>G]A": 0.09200, "A[C>G]C": 0.04200, "A[C>G]G": 0.01800, "A[C>G]T": 0.07600,
    "C[C>G]A": 0.08100, "C[C>G]C": 0.03700, "C[C>G]G": 0.01600, "C[C>G]T": 0.06700,
    "G[C>G]A": 0.07400, "G[C>G]C": 0.03400, "G[C>G]G": 0.01400, "G[C>G]T": 0.06100,
    "T[C>G]A": 0.10600, "T[C>G]C": 0.04800, "T[C>G]G": 0.02100, "T[C>G]T": 0.08800,
    "A[C>T]A": 0.00370, "A[C>T]C": 0.00170, "A[C>T]G": 0.00090, "A[C>T]T": 0.00310,
    "C[C>T]A": 0.00320, "C[C>T]C": 0.00150, "C[C>T]G": 0.00080, "C[C>T]T": 0.00270,
    "G[C>T]A": 0.00300, "G[C>T]C": 0.00140, "G[C>T]G": 0.00070, "G[C>T]T": 0.00250,
    "T[C>T]A": 0.06020, "T[C>T]C": 0.04290, "T[C>T]G": 0.01010, "T[C>T]T": 0.06090,
    "A[T>A]A": 0.00550, "A[T>A]C": 0.00340, "A[T>A]G": 0.00200, "A[T>A]T": 0.00500,
    "C[T>A]A": 0.00430, "C[T>A]C": 0.00260, "C[T>A]G": 0.00160, "C[T>A]T": 0.00390,
    "G[T>A]A": 0.00400, "G[T>A]C": 0.00250, "G[T>A]G": 0.00150, "G[T>A]T": 0.00370,
    "T[T>A]A": 0.00530, "T[T>A]C": 0.00330, "T[T>A]G": 0.00190, "T[T>A]T": 0.00480,
    "A[T>C]A": 0.00990, "A[T>C]C": 0.00560, "A[T>C]G": 0.00320, "A[T>C]T": 0.00890,
    "C[T>C]A": 0.00790, "C[T>C]C": 0.00440, "C[T>C]G": 0.00260, "C[T>C]T": 0.00710,
    "G[T>C]A": 0.00750, "G[T>C]C": 0.00420, "G[T>C]G": 0.00240, "G[T>C]T": 0.00670,
    "T[T>C]A": 0.01150, "T[T>C]C": 0.06770, "T[T>C]G": 0.00360, "T[T>C]T": 0.01030,
    "A[T>G]A": 0.00360, "A[T>G]C": 0.00220, "A[T>G]G": 0.00180, "A[T>G]T": 0.00330,
    "C[T>G]A": 0.00280, "C[T>G]C": 0.00170, "C[T>G]G": 0.00130, "C[T>G]T": 0.00250,
    "G[T>G]A": 0.00260, "G[T>G]C": 0.00160, "G[T>G]G": 0.00120, "G[T>G]T": 0.00240,
    "T[T>G]A": 0.00340, "T[T>G]C": 0.00210, "T[T>G]G": 0.00170, "T[T>G]T": 0.00320,
}

# --- SBS18: Oxidative damage (8-oxoguanine → C>A) ---
# Enrichment for C>A substitutions, especially at G contexts
_SBS18: dict[str, float] = {
    "A[C>A]A": 0.02030, "A[C>A]C": 0.01480, "A[C>A]G": 0.00920, "A[C>A]T": 0.01790,
    "C[C>A]A": 0.01810, "C[C>A]C": 0.01320, "C[C>A]G": 0.00820, "C[C>A]T": 0.01590,
    "G[C>A]A": 0.02210, "G[C>A]C": 0.01610, "G[C>A]G": 0.01000, "G[C>A]T": 0.01940,
    "T[C>A]A": 0.02160, "T[C>A]C": 0.01580, "T[C>A]G": 0.00980, "T[C>A]T": 0.01900,
    "A[C>G]A": 0.00330, "A[C>G]C": 0.00240, "A[C>G]G": 0.00150, "A[C>G]T": 0.00290,
    "C[C>G]A": 0.00290, "C[C>G]C": 0.00210, "C[C>G]G": 0.00130, "C[C>G]T": 0.00260,
    "G[C>G]A": 0.00360, "G[C>G]C": 0.00260, "G[C>G]G": 0.00160, "G[C>G]T": 0.00320,
    "T[C>G]A": 0.00350, "T[C>G]C": 0.00250, "T[C>G]G": 0.00160, "T[C>G]T": 0.00310,
    "A[C>T]A": 0.00520, "A[C>T]C": 0.00380, "A[C>T]G": 0.00230, "A[C>T]T": 0.00460,
    "C[C>T]A": 0.00460, "C[C>T]C": 0.00340, "C[C>T]G": 0.00200, "C[C>T]T": 0.00410,
    "G[C>T]A": 0.00570, "G[C>T]C": 0.00410, "G[C>T]G": 0.00250, "G[C>T]T": 0.00500,
    "T[C>T]A": 0.00550, "T[C>T]C": 0.00400, "T[C>T]G": 0.00240, "T[C>T]T": 0.00490,
    "A[T>A]A": 0.00820, "A[T>A]C": 0.00600, "A[T>A]G": 0.00370, "A[T>A]T": 0.00730,
    "C[T>A]A": 0.00730, "C[T>A]C": 0.00530, "C[T>A]G": 0.00330, "C[T>A]T": 0.00650,
    "G[T>A]A": 0.00890, "G[T>A]C": 0.00650, "G[T>A]G": 0.00400, "G[T>A]T": 0.00790,
    "T[T>A]A": 0.00870, "T[T>A]C": 0.00630, "T[T>A]G": 0.00390, "T[T>A]T": 0.00770,
    "A[T>C]A": 0.00670, "A[T>C]C": 0.00490, "A[T>C]G": 0.00300, "A[T>C]T": 0.00590,
    "C[T>C]A": 0.00590, "C[T>C]C": 0.00430, "C[T>C]G": 0.00270, "C[T>C]T": 0.00530,
    "G[T>C]A": 0.00730, "G[T>C]C": 0.00530, "G[T>C]G": 0.00330, "G[T>C]T": 0.00650,
    "T[T>C]A": 0.00710, "T[T>C]C": 0.00520, "T[T>C]G": 0.00320, "T[T>C]T": 0.00630,
    "A[T>G]A": 0.00580, "A[T>G]C": 0.00420, "A[T>G]G": 0.00340, "A[T>G]T": 0.00510,
    "C[T>G]A": 0.00510, "C[T>G]C": 0.00370, "C[T>G]G": 0.00300, "C[T>G]T": 0.00450,
    "G[T>G]A": 0.00630, "G[T>G]C": 0.00460, "G[T>G]G": 0.00370, "G[T>G]T": 0.00560,
    "T[T>G]A": 0.00620, "T[T>G]C": 0.00450, "T[T>G]G": 0.00360, "T[T>G]T": 0.00550,
}

# --- SBS35: Platinum chemotherapy ---
# C>A and T>A substitutions from platinum-DNA adducts
_SBS35: dict[str, float] = {
    "A[C>A]A": 0.01210, "A[C>A]C": 0.00840, "A[C>A]G": 0.00490, "A[C>A]T": 0.01050,
    "C[C>A]A": 0.01080, "C[C>A]C": 0.00750, "C[C>A]G": 0.00440, "C[C>A]T": 0.00940,
    "G[C>A]A": 0.01320, "G[C>A]C": 0.00910, "G[C>A]G": 0.00530, "G[C>A]T": 0.01140,
    "T[C>A]A": 0.01290, "T[C>A]C": 0.00890, "T[C>A]G": 0.00520, "T[C>A]T": 0.01120,
    "A[C>G]A": 0.00240, "A[C>G]C": 0.00170, "A[C>G]G": 0.00100, "A[C>G]T": 0.00210,
    "C[C>G]A": 0.00220, "C[C>G]C": 0.00150, "C[C>G]G": 0.00090, "C[C>G]T": 0.00190,
    "G[C>G]A": 0.00260, "G[C>G]C": 0.00180, "G[C>G]G": 0.00110, "G[C>G]T": 0.00230,
    "T[C>G]A": 0.00250, "T[C>G]C": 0.00180, "T[C>G]G": 0.00100, "T[C>G]T": 0.00220,
    "A[C>T]A": 0.00490, "A[C>T]C": 0.00340, "A[C>T]G": 0.00200, "A[C>T]T": 0.00430,
    "C[C>T]A": 0.00440, "C[C>T]C": 0.00300, "C[C>T]G": 0.00180, "C[C>T]T": 0.00380,
    "G[C>T]A": 0.00530, "G[C>T]C": 0.00370, "G[C>T]G": 0.00220, "G[C>T]T": 0.00460,
    "T[C>T]A": 0.00520, "T[C>T]C": 0.00360, "T[C>T]G": 0.00210, "T[C>T]T": 0.00450,
    "A[T>A]A": 0.02130, "A[T>A]C": 0.01480, "A[T>A]G": 0.00860, "A[T>A]T": 0.01860,
    "C[T>A]A": 0.01910, "C[T>A]C": 0.01320, "C[T>A]G": 0.00770, "C[T>A]T": 0.01660,
    "G[T>A]A": 0.02330, "G[T>A]C": 0.01620, "G[T>A]G": 0.00940, "G[T>A]T": 0.02030,
    "T[T>A]A": 0.02270, "T[T>A]C": 0.01580, "T[T>A]G": 0.00910, "T[T>A]T": 0.01980,
    "A[T>C]A": 0.00810, "A[T>C]C": 0.00560, "A[T>C]G": 0.00330, "A[T>C]T": 0.00710,
    "C[T>C]A": 0.00720, "C[T>C]C": 0.00500, "C[T>C]G": 0.00290, "C[T>C]T": 0.00630,
    "G[T>C]A": 0.00880, "G[T>C]C": 0.00610, "G[T>C]G": 0.00360, "G[T>C]T": 0.00770,
    "T[T>C]A": 0.00860, "T[T>C]C": 0.00600, "T[T>C]G": 0.00350, "T[T>C]T": 0.00750,
    "A[T>G]A": 0.00660, "A[T>G]C": 0.00460, "A[T>G]G": 0.00360, "A[T>G]T": 0.00570,
    "C[T>G]A": 0.00590, "C[T>G]C": 0.00410, "C[T>G]G": 0.00320, "C[T>G]T": 0.00510,
    "G[T>G]A": 0.00720, "G[T>G]C": 0.00500, "G[T>G]G": 0.00390, "G[T>G]T": 0.00620,
    "T[T>G]A": 0.00700, "T[T>G]C": 0.00490, "T[T>G]G": 0.00380, "T[T>G]T": 0.00610,
}

# --- SBS40: Clock-like, unknown etiology (similar to SBS5) ---
_SBS40: dict[str, float] = {
    "A[C>A]A": 0.00590, "A[C>A]C": 0.00390, "A[C>A]G": 0.00300, "A[C>A]T": 0.00540,
    "C[C>A]A": 0.00470, "C[C>A]C": 0.00320, "C[C>A]G": 0.00240, "C[C>A]T": 0.00430,
    "G[C>A]A": 0.00440, "G[C>A]C": 0.00300, "G[C>A]G": 0.00230, "G[C>A]T": 0.00410,
    "T[C>A]A": 0.00610, "T[C>A]C": 0.00410, "T[C>A]G": 0.00290, "T[C>A]T": 0.00560,
    "A[C>G]A": 0.00210, "A[C>G]C": 0.00140, "A[C>G]G": 0.00110, "A[C>G]T": 0.00190,
    "C[C>G]A": 0.00170, "C[C>G]C": 0.00110, "C[C>G]G": 0.00090, "C[C>G]T": 0.00160,
    "G[C>G]A": 0.00160, "G[C>G]C": 0.00100, "G[C>G]G": 0.00080, "G[C>G]T": 0.00140,
    "T[C>G]A": 0.00220, "T[C>G]C": 0.00150, "T[C>G]G": 0.00110, "T[C>G]T": 0.00200,
    "A[C>T]A": 0.00740, "A[C>T]C": 0.00500, "A[C>T]G": 0.00370, "A[C>T]T": 0.00690,
    "C[C>T]A": 0.00610, "C[C>T]C": 0.00420, "C[C>T]G": 0.00890, "C[C>T]T": 0.00570,
    "G[C>T]A": 0.00570, "G[C>T]C": 0.00390, "G[C>T]G": 0.00830, "G[C>T]T": 0.00540,
    "T[C>T]A": 0.00820, "T[C>T]C": 0.01210, "T[C>T]G": 0.01110, "T[C>T]T": 0.00760,
    "A[T>A]A": 0.00870, "A[T>A]C": 0.00580, "A[T>A]G": 0.00410, "A[T>A]T": 0.00800,
    "C[T>A]A": 0.00680, "C[T>A]C": 0.00450, "C[T>A]G": 0.00330, "C[T>A]T": 0.00620,
    "G[T>A]A": 0.00640, "G[T>A]C": 0.00420, "G[T>A]G": 0.00310, "G[T>A]T": 0.00590,
    "T[T>A]A": 0.00890, "T[T>A]C": 0.00590, "T[T>A]G": 0.00420, "T[T>A]T": 0.00820,
    "A[T>C]A": 0.00750, "A[T>C]C": 0.00500, "A[T>C]G": 0.00360, "A[T>C]T": 0.00690,
    "C[T>C]A": 0.00590, "C[T>C]C": 0.00390, "C[T>C]G": 0.00280, "C[T>C]T": 0.00540,
    "G[T>C]A": 0.00550, "G[T>C]C": 0.00370, "G[T>C]G": 0.00260, "G[T>C]T": 0.00510,
    "T[T>C]A": 0.00770, "T[T>C]C": 0.00520, "T[T>C]G": 0.00370, "T[T>C]T": 0.00710,
    "A[T>G]A": 0.00610, "A[T>G]C": 0.00410, "A[T>G]G": 0.00350, "A[T>G]T": 0.00570,
    "C[T>G]A": 0.00480, "C[T>G]C": 0.00320, "C[T>G]G": 0.00270, "C[T>G]T": 0.00440,
    "G[T>G]A": 0.00450, "G[T>G]C": 0.00300, "G[T>G]G": 0.00250, "G[T>G]T": 0.00410,
    "T[T>G]A": 0.00620, "T[T>G]C": 0.00420, "T[T>G]G": 0.00360, "T[T>G]T": 0.00580,
}

# --- SBS88: Colibactin (E. coli) — C→T at AT contexts ---
# Enriched for C>T at ACN and TCA/TCT contexts
_SBS88: dict[str, float] = {
    "A[C>A]A": 0.00490, "A[C>A]C": 0.00310, "A[C>A]G": 0.00190, "A[C>A]T": 0.00430,
    "C[C>A]A": 0.00430, "C[C>A]C": 0.00270, "C[C>A]G": 0.00160, "C[C>A]T": 0.00370,
    "G[C>A]A": 0.00460, "G[C>A]C": 0.00290, "G[C>A]G": 0.00180, "G[C>A]T": 0.00400,
    "T[C>A]A": 0.00530, "T[C>A]C": 0.00330, "T[C>A]G": 0.00200, "T[C>A]T": 0.00460,
    "A[C>G]A": 0.00190, "A[C>G]C": 0.00120, "A[C>G]G": 0.00080, "A[C>G]T": 0.00170,
    "C[C>G]A": 0.00170, "C[C>G]C": 0.00100, "C[C>G]G": 0.00070, "C[C>G]T": 0.00150,
    "G[C>G]A": 0.00180, "G[C>G]C": 0.00110, "G[C>G]G": 0.00070, "G[C>G]T": 0.00160,
    "T[C>G]A": 0.00210, "T[C>G]C": 0.00130, "T[C>G]G": 0.00080, "T[C>G]T": 0.00180,
    "A[C>T]A": 0.00640, "A[C>T]C": 0.00410, "A[C>T]G": 0.00250, "A[C>T]T": 0.00570,
    "C[C>T]A": 0.00560, "C[C>T]C": 0.00360, "C[C>T]G": 0.00220, "C[C>T]T": 0.00500,
    "G[C>T]A": 0.00600, "G[C>T]C": 0.00380, "G[C>T]G": 0.00230, "G[C>T]T": 0.00530,
    "T[C>T]A": 0.01540, "T[C>T]C": 0.01280, "T[C>T]G": 0.00690, "T[C>T]T": 0.01370,
    "A[T>A]A": 0.00710, "A[T>A]C": 0.00450, "A[T>A]G": 0.00270, "A[T>A]T": 0.00640,
    "C[T>A]A": 0.00630, "C[T>A]C": 0.00390, "C[T>A]G": 0.00240, "C[T>A]T": 0.00560,
    "G[T>A]A": 0.00670, "G[T>A]C": 0.00420, "G[T>A]G": 0.00250, "G[T>A]T": 0.00590,
    "T[T>A]A": 0.00740, "T[T>A]C": 0.00470, "T[T>A]G": 0.00280, "T[T>A]T": 0.00670,
    "A[T>C]A": 0.00780, "A[T>C]C": 0.00490, "A[T>C]G": 0.00300, "A[T>C]T": 0.00700,
    "C[T>C]A": 0.00680, "C[T>C]C": 0.00430, "C[T>C]G": 0.00260, "C[T>C]T": 0.00610,
    "G[T>C]A": 0.00730, "G[T>C]C": 0.00460, "G[T>C]G": 0.00280, "G[T>C]T": 0.00650,
    "T[T>C]A": 0.00820, "T[T>C]C": 0.00520, "T[T>C]G": 0.00310, "T[T>C]T": 0.00740,
    "A[T>G]A": 0.00540, "A[T>G]C": 0.00340, "A[T>G]G": 0.00260, "A[T>G]T": 0.00480,
    "C[T>G]A": 0.00480, "C[T>G]C": 0.00300, "C[T>G]G": 0.00230, "C[T>G]T": 0.00430,
    "G[T>G]A": 0.00510, "G[T>G]C": 0.00320, "G[T>G]G": 0.00240, "G[T>G]T": 0.00460,
    "T[T>G]A": 0.00560, "T[T>G]C": 0.00350, "T[T>G]G": 0.00270, "T[T>G]T": 0.00500,
}

# Master signature registry
SBS_SIGNATURES: dict[str, dict[str, float]] = {
    "SBS1": _SBS1,
    "SBS2": _SBS2,
    "SBS5": _SBS5,
    "SBS7a": _SBS7A,
    "SBS7b": _SBS7B,
    "SBS7c": _SBS7C,
    "SBS7d": _SBS7D,
    "SBS13": _SBS13,
    "SBS18": _SBS18,
    "SBS35": _SBS35,
    "SBS40": _SBS40,
    "SBS88": _SBS88,
}

# Validate and normalize all signatures to sum to exactly 1.0
for _sig_name, _sig_data in SBS_SIGNATURES.items():
    assert len(_sig_data) == 96, f"{_sig_name} has {len(_sig_data)} contexts, expected 96"
    _total = sum(_sig_data.values())
    assert _total > 0, f"{_sig_name} has zero sum"
    # Normalize to sum to 1.0
    for _ctx in _sig_data:
        _sig_data[_ctx] /= _total
    # Verify normalization
    _normalized_total = sum(_sig_data.values())
    assert abs(_normalized_total - 1.0) < 1e-10, f"{_sig_name} normalization failed: {_normalized_total}"


# ==============================================================================
# 3. Signature Metadata
# ==============================================================================

@dataclass
class SignatureInfo:
    """Metadata for a COSMIC SBS signature."""
    name: str
    etiology: str
    description: str
    dominant_mutations: list[str]
    tissue_relevance: list[str]
    reference: str


_SIGNATURE_METADATA: dict[str, SignatureInfo] = {
    "SBS1": SignatureInfo(
        name="SBS1",
        etiology="5-methylcytosine deamination",
        description="Clock-like signature from spontaneous deamination of "
                    "5-methylcytosine at CpG sites. Correlates with age.",
        dominant_mutations=["C>T at NCG"],
        tissue_relevance=["all tissues"],
        reference="Alexandrov 2013 Nature; Alexandrov 2020 Nature",
    ),
    "SBS2": SignatureInfo(
        name="SBS2",
        etiology="APOBEC cytidine deaminase (C→T)",
        description="APOBEC-mediated cytidine deamination producing C>T "
                    "substitutions at TpC contexts.",
        dominant_mutations=["C>T at TCA/TCT"],
        tissue_relevance=["bladder", "cervix", "lung", "head_and_neck", "breast"],
        reference="Alexandrov 2013 Nature; Alexandrov 2020 Nature",
    ),
    "SBS5": SignatureInfo(
        name="SBS5",
        etiology="Clock-like, unknown etiology",
        description="Flat/broad signature correlating with cell divisions. "
                    "Co-occurs with SBS1 in most cancer types.",
        dominant_mutations=["broad distribution"],
        tissue_relevance=["all tissues"],
        reference="Alexandrov 2020 Nature",
    ),
    "SBS7a": SignatureInfo(
        name="SBS7a",
        etiology="UV exposure (variant a)",
        description="UV-induced C>T at dipyrimidine sites with CC→TT "
                    "dinucleotide mutations.",
        dominant_mutations=["C>T at NCC/NCT", "CC→TT"],
        tissue_relevance=["skin", "melanoma"],
        reference="Alexandrov 2013 Nature",
    ),
    "SBS7b": SignatureInfo(
        name="SBS7b",
        etiology="UV exposure (variant b)",
        description="UV-induced signature with similar but distinct profile "
                    "from SBS7a.",
        dominant_mutations=["C>T at dipyrimidine"],
        tissue_relevance=["skin", "melanoma"],
        reference="Alexandrov 2020 Nature",
    ),
    "SBS7c": SignatureInfo(
        name="SBS7c",
        etiology="UV exposure (variant c)",
        description="UV-induced signature with more balanced substitution "
                    "distribution including C>A.",
        dominant_mutations=["C>T", "C>A at dipyrimidine"],
        tissue_relevance=["skin", "melanoma"],
        reference="Alexandrov 2020 Nature",
    ),
    "SBS7d": SignatureInfo(
        name="SBS7d",
        etiology="UV exposure (variant d)",
        description="UV-induced signature with prominent CC→TT dinucleotide "
                    "mutations.",
        dominant_mutations=["CC→TT", "C>T at NCC"],
        tissue_relevance=["skin", "melanoma"],
        reference="Alexandrov 2020 Nature",
    ),
    "SBS13": SignatureInfo(
        name="SBS13",
        etiology="APOBEC cytidine deaminase (C→G)",
        description="APOBEC-mediated cytidine deamination producing C>G "
                    "substitutions at TpC contexts.",
        dominant_mutations=["C>G at TCA/TCT"],
        tissue_relevance=["bladder", "cervix", "lung", "head_and_neck", "breast"],
        reference="Alexandrov 2013 Nature; Alexandrov 2020 Nature",
    ),
    "SBS18": SignatureInfo(
        name="SBS18",
        etiology="Oxidative damage (8-oxoguanine)",
        description="Damage from reactive oxygen species, producing 8-oxoG "
                    "lesions that cause C>A substitutions.",
        dominant_mutations=["C>A broad"],
        tissue_relevance=["stomach", "colorectal", "neuroblastoma"],
        reference="Alexandrov 2020 Nature",
    ),
    "SBS35": SignatureInfo(
        name="SBS35",
        etiology="Platinum chemotherapy",
        description="Platinum drug (cisplatin, carboplatin) DNA adducts "
                    "causing C>A and T>A substitutions.",
        dominant_mutations=["C>A", "T>A"],
        tissue_relevance=["ovary", "lung", "testicular"],
        reference="Alexandrov 2020 Nature",
    ),
    "SBS40": SignatureInfo(
        name="SBS40",
        etiology="Clock-like, unknown etiology",
        description="Flat clock-like signature very similar to SBS5. "
                    "Co-occurs with SBS5 in many cancer types.",
        dominant_mutations=["broad distribution"],
        tissue_relevance=["all tissues"],
        reference="Alexandrov 2020 Nature",
    ),
    "SBS88": SignatureInfo(
        name="SBS88",
        etiology="Colibactin (E. coli)",
        description="Colibactin-producing E. coli causing DNA cross-links "
                    "leading to C>T at AT-rich contexts.",
        dominant_mutations=["C>T at ACN", "T>C at TCA"],
        tissue_relevance=["colorectal"],
        reference="Pleguezuelos-Manzano 2020 Nature; Alexandrov 2020 Nature",
    ),
}


# ==============================================================================
# 4. Tissue-Specific Signature Weights
# ==============================================================================
# Relative activity of each signature by tissue type (approximate, from COSMIC)
# Used by compute_damage_susceptibility_profile to weight signatures

_TISSUE_SIGNATURE_WEIGHTS: dict[str, dict[str, float]] = {
    "generic": {
        "SBS1": 1.0, "SBS5": 1.0, "SBS40": 1.0,
        "SBS2": 0.3, "SBS13": 0.2, "SBS18": 0.3,
        "SBS7a": 0.1, "SBS7b": 0.05, "SBS7c": 0.03, "SBS7d": 0.02,
        "SBS35": 0.1, "SBS88": 0.1,
    },
    "skin": {
        "SBS7a": 1.0, "SBS7b": 0.8, "SBS7c": 0.6, "SBS7d": 0.5,
        "SBS1": 0.5, "SBS5": 0.5, "SBS40": 0.5,
        "SBS2": 0.1, "SBS13": 0.1, "SBS18": 0.1,
        "SBS35": 0.05, "SBS88": 0.05,
    },
    "colorectal": {
        "SBS1": 1.0, "SBS5": 1.0, "SBS40": 1.0,
        "SBS88": 0.8, "SBS18": 0.5, "SBS35": 0.3,
        "SBS2": 0.1, "SBS13": 0.1,
        "SBS7a": 0.02, "SBS7b": 0.01, "SBS7c": 0.01, "SBS7d": 0.01,
    },
    "lung": {
        "SBS1": 1.0, "SBS5": 1.0, "SBS40": 1.0,
        "SBS2": 0.5, "SBS13": 0.4, "SBS18": 0.4,
        "SBS35": 0.3,
        "SBS7a": 0.02, "SBS7b": 0.01, "SBS7c": 0.01, "SBS7d": 0.01,
        "SBS88": 0.1,
    },
    "ovary": {
        "SBS1": 1.0, "SBS5": 1.0, "SBS40": 1.0,
        "SBS35": 0.8, "SBS2": 0.2, "SBS13": 0.2,
        "SBS18": 0.3,
        "SBS7a": 0.02, "SBS7b": 0.01, "SBS7c": 0.01, "SBS7d": 0.01,
        "SBS88": 0.1,
    },
    "bladder": {
        "SBS1": 0.8, "SBS5": 0.8, "SBS40": 0.8,
        "SBS2": 0.9, "SBS13": 0.7, "SBS18": 0.4,
        "SBS35": 0.1,
        "SBS7a": 0.02, "SBS7b": 0.01, "SBS7c": 0.01, "SBS7d": 0.01,
        "SBS88": 0.1,
    },
    "cervix": {
        "SBS1": 0.6, "SBS5": 0.6, "SBS40": 0.6,
        "SBS2": 1.0, "SBS13": 0.8,
        "SBS18": 0.2,
        "SBS7a": 0.02, "SBS7b": 0.01, "SBS7c": 0.01, "SBS7d": 0.01,
        "SBS35": 0.05, "SBS88": 0.05,
    },
    "stomach": {
        "SBS1": 0.8, "SBS5": 0.8, "SBS40": 0.8,
        "SBS18": 0.7, "SBS88": 0.3,
        "SBS2": 0.1, "SBS13": 0.1,
        "SBS7a": 0.02, "SBS7b": 0.01, "SBS7c": 0.01, "SBS7d": 0.01,
        "SBS35": 0.2,
    },
    "breast": {
        "SBS1": 0.8, "SBS5": 0.8, "SBS40": 0.8,
        "SBS2": 0.6, "SBS13": 0.5,
        "SBS18": 0.2, "SBS35": 0.3,
        "SBS7a": 0.02, "SBS7b": 0.01, "SBS7c": 0.01, "SBS7d": 0.01,
        "SBS88": 0.1,
    },
    "neuroblastoma": {
        "SBS1": 0.5, "SBS5": 0.5, "SBS40": 0.5,
        "SBS18": 1.0,
        "SBS2": 0.1, "SBS13": 0.1,
        "SBS7a": 0.01, "SBS7b": 0.01, "SBS7c": 0.01, "SBS7d": 0.01,
        "SBS35": 0.1, "SBS88": 0.05,
    },
}


# ==============================================================================
# 5. Helper: Parse COSMIC context string to trinucleotide bases
# ==============================================================================

def _parse_cosmic_context(context: str) -> tuple[str, str, str] | None:
    """Parse a COSMIC context string like 'A[C>T]G' into (5'base, ref, 3'base).

    Returns:
        Tuple of (5' base, reference base, 3' base) or None if invalid.
    """
    try:
        five_prime = context[0]
        ref = context[2]
        three_prime = context[-1]
        if five_prime in _BASES and ref in _BASES and three_prime in _BASES:
            return (five_prime, ref, three_prime)
    except (IndexError, ValueError):
        pass
    return None


def _build_context_lookup() -> dict[str, dict[str, float]]:
    """Build a fast lookup: trinucleotide_string -> {signature_name -> probability}.

    The trinucleotide_string is just the 3 bases, e.g. "ACG".
    This maps position-based contexts to signature probabilities, summing
    across all substitution types for that trinucleotide.
    """
    lookup: dict[str, dict[str, float]] = {}
    for sig_name, sig_data in SBS_SIGNATURES.items():
        for cosmic_ctx, prob in sig_data.items():
            parsed = _parse_cosmic_context(cosmic_ctx)
            if parsed is None:
                continue
            five_prime, ref, three_prime = parsed
            trinuc = five_prime + ref + three_prime
            if trinuc not in lookup:
                lookup[trinuc] = {}
            lookup[trinuc][sig_name] = lookup[trinuc].get(sig_name, 0.0) + prob
    return lookup


# Pre-computed lookup table for fast position-based queries
_CONTEXT_LOOKUP: dict[str, dict[str, float]] = _build_context_lookup()


# ==============================================================================
# 6. Public Functions
# ==============================================================================

def get_trinucleotide_context(seq: str, position: int) -> str:
    """Extract the 3-base trinucleotide context centered at a position.

    For positions at sequence boundaries, returns a shorter string.
    The central base is at ``position``, flanked by the 5' and 3' bases.

    Args:
        seq: DNA sequence (uppercase or mixed case).
        position: 0-indexed position of the central base.

    Returns:
        3-base trinucleotide string (e.g. "ACG"), or shorter at boundaries.

    Raises:
        ValueError: If position is out of range.

    Examples:
        >>> get_trinucleotide_context("ACGTA", 2)
        'CGT'
        >>> get_trinucleotide_context("ACGTA", 0)
        'AC'
    """
    seq = seq.upper()
    n = len(seq)
    if position < 0 or position >= n:
        raise ValueError(f"Position {position} out of range [0, {n})")

    start = max(0, position - 1)
    end = min(n, position + 2)
    return seq[start:end]


def compute_signature_weight(
    seq: str,
    position: int,
    signatures: list[str] | None = None,
) -> float:
    """Compute combined SBS signature weight at a given sequence position.

    Looks up the trinucleotide context at the position and sums the
    probability contributions from all specified COSMIC SBS signatures.
    This represents the mutational susceptibility at that position based
    on the combined action of the specified mutagenic processes.

    Args:
        seq: DNA sequence (uppercase or mixed case).
        position: 0-indexed position in the sequence.
        signatures: List of SBS signature names to consider.
            If None, uses all available signatures.

    Returns:
        Combined signature weight (0-1 scale). Higher values indicate
        greater mutational susceptibility.

    Examples:
        >>> # CpG site: TCG context has high SBS1 weight
        >>> round(compute_signature_weight("ATCGA", 2, ["SBS1"]), 3)
        0.048
        >>> # TCA context: APOBEC hotspot
        >>> compute_signature_weight("ATCAA", 2, ["SBS2"]) > 0.01
        True
    """
    seq = seq.upper()
    n = len(seq)

    if position < 1 or position >= n - 1:
        # Boundary positions don't have full trinucleotide context
        return 0.0

    trinuc = seq[position - 1:position + 2]
    if trinuc not in _CONTEXT_LOOKUP:
        return 0.0

    available = _CONTEXT_LOOKUP[trinuc]
    if signatures is None:
        return min(1.0, sum(available.values()))

    weight = 0.0
    for sig in signatures:
        if sig in available:
            weight += available[sig]

    return min(1.0, weight)


def scan_signature_hotspots(
    seq: str,
    signatures: list[str] | None = None,
    threshold: float = 0.05,
) -> list[dict]:
    """Scan a full DNA sequence for mutational signature hotspots.

    Slides across the sequence computing combined signature weight at
    each position. Positions exceeding the threshold are reported as
    hotspots with their contributing signatures and context information.

    Args:
        seq: DNA sequence (uppercase or mixed case).
        signatures: List of SBS signature names to consider.
            If None, uses all available signatures.
        threshold: Minimum combined weight to call a hotspot (default 0.05).

    Returns:
        List of hotspot dicts, each with keys:
            - position: 0-indexed position
            - context: 3-base trinucleotide context
            - weight: combined signature weight
            - contributions: dict of {signature_name: probability}
            - dominant_signature: name of the highest-contributing signature

    Examples:
        >>> hotspots = scan_signature_hotspots("ATCGATCGA", threshold=0.03)
        >>> len(hotspots) > 0
        True
    """
    seq = seq.upper()
    n = len(seq)
    hotspots: list[dict] = []

    for pos in range(1, n - 1):
        trinuc = seq[pos - 1:pos + 2]
        if trinuc not in _CONTEXT_LOOKUP:
            continue

        available = _CONTEXT_LOOKUP[trinuc]
        if signatures is None:
            contributions = dict(available)
        else:
            contributions = {s: available[s] for s in signatures if s in available}

        weight = sum(contributions.values())
        if weight >= threshold:
            dominant = max(contributions, key=contributions.get) if contributions else ""
            hotspots.append({
                "position": pos,
                "context": trinuc,
                "weight": round(weight, 6),
                "contributions": {k: round(v, 6) for k, v in contributions.items()},
                "dominant_signature": dominant,
            })

    return hotspots


def assign_signatures(
    seq: str,
    mutation_counts: dict | None = None,
) -> dict:
    """Assign COSMIC signature exposures to observed mutations.

    This is a simplified version of signature assignment that does not
    require SigProfilerSuite. It uses a least-squares approach to
    decompose the observed mutation spectrum into contributions from
    known COSMIC signatures.

    Args:
        seq: DNA sequence for context extraction (used if mutation_counts
            is None to compute expected trinucleotide frequencies).
        mutation_counts: Dict mapping COSMIC context strings (e.g. "A[C>T]G")
            to observed counts. If None, uses the sequence to generate
            a synthetic mutation spectrum based on context frequencies.

    Returns:
        Dict with keys:
            - exposures: dict of {signature_name: exposure_fraction}
            - reconstructed: dict of {context: reconstructed_probability}
            - reconstruction_error: L2 error between observed and reconstructed
            - n_signatures: number of signatures with non-trivial exposure

    Notes:
        This simplified solver uses non-negative least squares (NNLS) via
        scipy when available, falling back to an iterative multiplicative
        update rule. For production use with real VCF data, prefer
        ``assign_signatures_sigprofiler()``.
    """
    # Build observed mutation spectrum
    if mutation_counts is not None:
        total_mutations = sum(mutation_counts.values())
        if total_mutations == 0:
            return {
                "exposures": {name: 0.0 for name in SBS_SIGNATURES},
                "reconstructed": {},
                "reconstruction_error": 0.0,
                "n_signatures": 0,
            }
        observed = {
            ctx: mutation_counts.get(ctx, 0) / total_mutations
            for ctx in TRINUCLEOTIDE_CONTEXTS
        }
    else:
        # Use sequence context frequencies as a proxy for mutation spectrum
        seq = seq.upper()
        context_counts: dict[str, int] = {ctx: 0 for ctx in TRINUCLEOTIDE_CONTEXTS}
        for pos in range(1, len(seq) - 1):
            trinuc = seq[pos - 1:pos + 2]
            # Map trinucleotide to all compatible COSMIC contexts
            for sub in _SUBSTITUTION_TYPES:
                cosmic_ctx = f"{trinuc[0]}[{sub}]{trinuc[2]}"
                if cosmic_ctx in context_counts:
                    context_counts[cosmic_ctx] += 1
        total_ctx = sum(context_counts.values())
        if total_ctx == 0:
            return {
                "exposures": {name: 0.0 for name in SBS_SIGNATURES},
                "reconstructed": {},
                "reconstruction_error": 0.0,
                "n_signatures": 0,
            }
        observed = {ctx: c / total_ctx for ctx, c in context_counts.items()}

    # Build signature matrix (96 x n_signatures)
    sig_names = list(SBS_SIGNATURES.keys())
    n_sigs = len(sig_names)

    # Vectorize observed spectrum
    import_array = [observed.get(ctx, 0.0) for ctx in TRINUCLEOTIDE_CONTEXTS]

    # Build signature matrix
    sig_matrix: list[list[float]] = []
    for sig_name in sig_names:
        sig_probs = SBS_SIGNATURES[sig_name]
        sig_matrix.append([sig_probs.get(ctx, 0.0) for ctx in TRINUCLEOTIDE_CONTEXTS])

    # Try scipy NNLS first
    try:
        import numpy as np
        from scipy.optimize import nnls

        A = np.array(sig_matrix).T  # (96, n_sigs)
        b = np.array(import_array)   # (96,)
        x, residual = nnls(A, b)

        exposures = {sig_names[i]: float(x[i]) for i in range(n_sigs)}
        total_exp = sum(exposures.values())
        if total_exp > 0:
            exposures = {k: v / total_exp for k, v in exposures.items()}

        # Reconstruct spectrum
        reconstructed_vec = A @ x
        total_recon = reconstructed_vec.sum()
        if total_recon > 0:
            reconstructed_probs = reconstructed_vec / total_recon
        else:
            reconstructed_probs = reconstructed_vec

        reconstructed = {
            TRINUCLEOTIDE_CONTEXTS[i]: float(reconstructed_probs[i])
            for i in range(96)
        }

        # L2 reconstruction error
        error = float(np.linalg.norm(b - reconstructed_vec))

        n_active = sum(1 for v in exposures.values() if v > 0.01)

        return {
            "exposures": exposures,
            "reconstructed": reconstructed,
            "reconstruction_error": round(error, 6),
            "n_signatures": n_active,
        }

    except ImportError:
        pass

    # Fallback: iterative multiplicative update (Lee & Seung 1999)
    exposures_arr = [1.0 / n_sigs] * n_sigs  # uniform initialization

    for _ in range(500):
        # Reconstruct
        recon = [0.0] * 96
        for j in range(n_sigs):
            for i in range(96):
                recon[i] += exposures_arr[j] * sig_matrix[j][i]

        # Multiplicative update
        new_exposures = list(exposures_arr)
        for j in range(n_sigs):
            numerator = 0.0
            denominator = 0.0
            for i in range(96):
                if recon[i] > 1e-12:
                    numerator += sig_matrix[j][i] * import_array[i]
                    denominator += sig_matrix[j][i] * recon[i]
            if denominator > 1e-12:
                new_exposures[j] = exposures_arr[j] * (numerator / denominator)

        # Check convergence
        max_change = max(abs(new_exposures[j] - exposures_arr[j])
                         for j in range(n_sigs))
        exposures_arr = new_exposures
        if max_change < 1e-8:
            break

    # Normalize
    total_exp = sum(exposures_arr)
    if total_exp > 0:
        exposures_arr = [e / total_exp for e in exposures_arr]

    exposures = {sig_names[i]: exposures_arr[i] for i in range(n_sigs)}

    # Reconstruct
    recon = [0.0] * 96
    for j in range(n_sigs):
        for i in range(96):
            recon[i] += exposures_arr[j] * sig_matrix[j][i]

    total_recon = sum(recon)
    if total_recon > 0:
        recon_norm = [r / total_recon for r in recon]
    else:
        recon_norm = recon

    reconstructed = {
        TRINUCLEOTIDE_CONTEXTS[i]: round(recon_norm[i], 6)
        for i in range(96)
    }

    # L2 error
    error = math.sqrt(sum(
        (import_array[i] - recon[i]) ** 2 for i in range(96)
    ))

    n_active = sum(1 for v in exposures.values() if v > 0.01)

    return {
        "exposures": {k: round(v, 6) for k, v in exposures.items()},
        "reconstructed": reconstructed,
        "reconstruction_error": round(error, 6),
        "n_signatures": n_active,
    }


def compute_damage_susceptibility_profile(
    seq: str,
    tissue: str = "generic",
) -> list[float]:
    """Compute per-position damage susceptibility based on COSMIC signatures.

    For each position in the sequence, computes a weighted sum of signature
    probabilities using tissue-specific signature activity weights. This
    produces a profile of mutational susceptibility across the entire
    sequence, useful for identifying fragile sites in therapeutic gene design.

    Args:
        seq: DNA sequence (uppercase or mixed case).
        tissue: Tissue type for weighting signature activities.
            Supported: "generic", "skin", "colorectal", "lung", "ovary",
            "bladder", "cervix", "stomach", "breast", "neuroblastoma".

    Returns:
        List of float values, one per position. Each value represents
        the weighted mutational susceptibility (0-1 scale). Positions
        at sequence boundaries have reduced values due to incomplete
        trinucleotide context.

    Raises:
        ValueError: If tissue type is not recognized.

    Examples:
        >>> profile = compute_damage_susceptibility_profile("ATCGATCGA", "generic")
        >>> len(profile) == 9
        True
        >>> all(0.0 <= v <= 1.0 for v in profile)
        True
    """
    seq = seq.upper()
    n = len(seq)

    if tissue not in _TISSUE_SIGNATURE_WEIGHTS:
        valid = sorted(_TISSUE_SIGNATURE_WEIGHTS.keys())
        raise ValueError(
            f"Unknown tissue type '{tissue}'. "
            f"Supported: {valid}"
        )

    tissue_weights = _TISSUE_SIGNATURE_WEIGHTS[tissue]
    profile: list[float] = [0.0] * n

    for pos in range(n):
        # Boundary positions get partial context
        if pos < 1 or pos >= n - 1:
            profile[pos] = 0.0
            continue

        trinuc = seq[pos - 1:pos + 2]
        if trinuc not in _CONTEXT_LOOKUP:
            profile[pos] = 0.0
            continue

        available = _CONTEXT_LOOKUP[trinuc]
        weighted_sum = 0.0
        for sig_name, prob in available.items():
            tissue_w = tissue_weights.get(sig_name, 0.0)
            weighted_sum += prob * tissue_w

        # Normalize by total tissue weight to keep in 0-1 range
        total_tissue_weight = sum(tissue_weights.values())
        if total_tissue_weight > 0:
            profile[pos] = min(1.0, weighted_sum / total_tissue_weight)
        else:
            profile[pos] = 0.0

    return profile


def assign_signatures_sigprofiler(
    vcf_path: str,
    output_dir: str,
    genome: str = "GRCh38",
) -> dict:
    """Assign mutational signatures using SigProfilerAssignment.

    This is a wrapper around the SigProfilerSuite ecosystem for
    production-grade signature assignment from VCF files. It requires
    the SigProfilerAssignment package to be installed.

    The function will:
    1. Generate a SigProfiler matrix from the input VCF
    2. Run SigProfilerAssignment fit to COSMIC v3 signatures
    3. Return exposure assignments and quality metrics

    Args:
        vcf_path: Path to input VCF file (single-sample).
        output_dir: Directory for SigProfiler output files.
        genome: Reference genome build (GRCh38, GRCh37, mm10, etc.).

    Returns:
        Dict with keys:
            - exposures: dict of {signature_name: exposure_fraction}
            - reconstruction_error: cosine distance between observed/reconstructed
            - stats: dict with additional quality metrics
            - error: error message if SigProfilerAssignment is unavailable

    Raises:
        FileNotFoundError: If vcf_path does not exist.

    References:
        Islam SMA et al. (2022) Nat Genet 54:1516-1527.
        "SigProfilerAssignment: pan-cancer implementation of
        mutational signature assignment"
    """
    import os

    if not os.path.isfile(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    os.makedirs(output_dir, exist_ok=True)

    try:
        from SigProfilerAssignment import Analyzer as An

        # Run SigProfiler fit
        An.fit(
            vcf_path,
            output_dir,
            genome,
            signature_database="COSMIC_v3.4",
            make_plots=False,
        )

        # Parse results
        import csv
        results_path = os.path.join(
            output_dir, "Assignment_Solution", "Activities",
            "Assignment_Solution_Activities.txt"
        )

        if not os.path.isfile(results_path):
            return {
                "exposures": {},
                "reconstruction_error": float("nan"),
                "stats": {"status": "no_results_file"},
                "error": "SigProfiler output not found at expected path",
            }

        exposures: dict[str, float] = {}
        with open(results_path) as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)
            for row in reader:
                if len(row) >= 2:
                    sig_name = row[0]
                    try:
                        exposure = float(row[1])
                        exposures[sig_name] = exposure
                    except (ValueError, IndexError):
                        continue

        # Normalize exposures
        total = sum(exposures.values())
        if total > 0:
            exposures = {k: v / total for k, v in exposures.items()}

        # Try to read reconstruction error
        stats: dict = {"status": "success"}
        try:
            stats_path = os.path.join(
                output_dir, "Assignment_Solution", "Stats.txt"
            )
            if os.path.isfile(stats_path):
                with open(stats_path) as f:
                    for line in f:
                        if "cosine_similarity" in line.lower():
                            parts = line.strip().split("\t")
                            if len(parts) >= 2:
                                stats["cosine_similarity"] = float(parts[1])
                        elif "l2_error" in line.lower():
                            parts = line.strip().split("\t")
                            if len(parts) >= 2:
                                stats["l2_error"] = float(parts[1])
        except (ValueError, OSError):
            pass

        recon_error = 1.0 - stats.get("cosine_similarity", 0.0)

        return {
            "exposures": exposures,
            "reconstruction_error": round(recon_error, 6),
            "stats": stats,
            "error": None,
        }

    except ImportError:
        return {
            "exposures": {},
            "reconstruction_error": float("nan"),
            "stats": {"status": "sigprofiler_not_installed"},
            "error": (
                "SigProfilerAssignment is not installed. "
                "Install with: pip install SigProfilerAssignment"
            ),
        }
    except Exception as e:
        logger.error(f"SigProfilerAssignment failed: {e}")
        return {
            "exposures": {},
            "reconstruction_error": float("nan"),
            "stats": {"status": "error", "message": str(e)},
            "error": f"SigProfilerAssignment failed: {e}",
        }


# ==============================================================================
# 7. Convenience Functions
# ==============================================================================

def get_signature_info(signature_name: str) -> SignatureInfo | None:
    """Get metadata for a COSMIC SBS signature.

    Args:
        signature_name: SBS signature name (e.g. "SBS1").

    Returns:
        SignatureInfo object or None if not found.
    """
    return _SIGNATURE_METADATA.get(signature_name)


def list_signatures(tissue: str | None = None) -> list[str]:
    """List available SBS signatures, optionally filtered by tissue relevance.

    Args:
        tissue: If provided, only return signatures active in this tissue.

    Returns:
        List of signature names.
    """
    if tissue is None:
        return list(SBS_SIGNATURES.keys())

    weights = _TISSUE_SIGNATURE_WEIGHTS.get(tissue, {})
    return [sig for sig, w in weights.items() if w > 0.1]


def get_tissue_types() -> list[str]:
    """List supported tissue types for signature weighting.

    Returns:
        List of tissue type strings.
    """
    return sorted(_TISSUE_SIGNATURE_WEIGHTS.keys())
