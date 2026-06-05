"""
BioCompiler mRNA Stability Prediction Models
=============================================

Organism-specific models for scoring mRNA stability based on
codon optimality and cis-regulatory sequence motifs.

mRNA stability is a critical determinant of heterologous protein expression.
Two complementary mechanisms are modelled:

1. **Codon optimality** — Optimal (high-adaptiveness) codons promote
   efficient translation and mRNA stability; non-optimal codons lead to
   ribosome stalling and accelerated decay.  The Codon Adaptation Index
   (CAI) is used as a primary proxy for stability (Sharp & Li 1987;
   Presnyak et al. 2015).

2. **Cis-regulatory motifs** — AU-rich elements (AREs), RNase recognition
   sites, and other decay-promoting motifs can severely limit yield even
   when codon usage is fully optimized.

Three organism models are provided:
  - *E. coli*:  RNase E sites, stem-loop stabilizers; stable if CAI > 0.7
    (Belasco & Higgins 1988; Carpousis 2010).
  - *Human / mammalian*:  ARE subtypes (ATTTA = instability motif);
    stable if CAI > 0.7 AND no ATTTA motifs (Shaw & Kamen 1986;
    Chen & Shyu 1995; Barreau et al. 2005).
  - *S. cerevisiae*:  PGK1 stabilizing element, rapid-decay elements;
    similar thresholds to human (Muhlrad & Parker 1992).

Usage::

    from biocompiler.mrna_stability import (
        score_mrna_stability,
        compute_mrna_half_life_score,
        predict_mrna_stability,
        suggest_mutations_for_stability,
    )

    # Predict stability category
    dna = "ATGGCCGCGATTTAGCGGCGCCTAA"
    category = predict_mrna_stability(dna, organism="Homo_sapiens")
    print(f"Stability: {category}")  # STABLE, MODERATE, or UNSTABLE

    # Get a half-life score (0-1)
    score = compute_mrna_half_life_score(dna, organism="Escherichia_coli")
    print(f"Half-life score: {score:.3f}")

    # Detailed motif-based scoring
    result = score_mrna_stability(dna, organism="Homo_sapiens")
    print(f"Stability score: {result.overall_score:.3f} (risk: {result.risk_level})")
    print(f"  Stabilizing motifs: {result.stabilizing_count}")
    print(f"  Destabilizing motifs: {result.destabilizing_count}")

    # Get suggestions for improving stability via synonymous codon changes
    suggestions = suggest_mutations_for_stability(dna, organism="Homo_sapiens")
    for s in suggestions:
        print(f"  Position {s['position']}: {s['original_codon']} -> "
              f"{s['suggested_codon']} ({s['amino_acid']}) "
              f"removes {s['motif_removed']}")

Public API
----------
- ``STABILITY_MOTIFS`` — per-organism motif catalogue
- ``MRNAStabilityScore`` — dataclass returned by the scorer
- ``score_mrna_stability()`` — scan a DNA sequence and return a composite score
- ``compute_mrna_half_life_score()`` — CAI-based half-life proxy score (0-1)
- ``predict_mrna_stability()`` — return STABLE / MODERATE / UNSTABLE category
- ``suggest_mutations_for_stability()`` — synonymous-codon fixes

Notes
-----
All motif patterns are expressed in **DNA** (T, not U).  The input
``dna`` argument is the coding strand (5'→3').  UTR positions are
approximated from the CDS boundaries: the first ~10 % of the sequence
is treated as the 5' UTR proxy and the last ~15 % as the 3' UTR proxy.
Callers who know the exact UTR boundaries should pre-slice and pass
only the relevant region.

References:
  Sharp, P.M. & Li, W.-H. (1987). "The codon Adaptation Index - a measure
  of directional synonymous codon usage bias." *Nucleic Acids Research*
  15:1281–1295.
  Presnyak, V. et al. (2015). "Codon optimality is a major determinant of
  mRNA stability." *Cell* 160:1111–1124.
  Shaw, G. & Kamen, R. (1986). "A conserved AU sequence from the 3'
  untranslated region of GM-CSF mRNA mediates selective mRNA degradation."
  *Cell* 46:659–667.
  Chen, C.-Y.A. & Shyu, A.-B. (1995). "AU-rich elements: characterization
  and importance in mRNA degradation." *Trends in Biochemical Sciences*
  20:465–470.
  Belasco, J.G. & Higgins, C.F. (1988). "Mechanisms of mRNA decay in
  bacteria." *Microbiological Reviews* 52:411–425.
  Barreau, C., Paillard, L., & Osborne, H.B. (2005). "AU-rich elements
  and associated factors: are there unifying principles?" *Nucleic Acids
  Research* 33:7138–7150.
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass, field
from typing import Any, Literal

from .constants import AA_TO_CODONS, CODON_TABLE

logger = logging.getLogger(__name__)

__all__ = [
    "STABILITY_MOTIFS",
    "MRNAStabilityScore",
    "score_mrna_stability",
    "compute_mrna_half_life_score",
    "predict_mrna_stability",
    "suggest_mutations_for_stability",
]

# Type alias for stability category
StabilityCategory = Literal["STABLE", "MODERATE", "UNSTABLE"]

# ────────────────────────────────────────────────────────────
# 1. Motif catalogue
# ────────────────────────────────────────────────────────────

# Each motif entry is a tuple: (pattern, description, source_paper)
# Patterns are DNA regexes (T not U).  They are matched against the
# upper-cased input sequence with re.finditer.

_ECOLI_MOTIFS: dict[str, list[tuple[str, str, str]]] = {
    "stabilizing": [
        # 5' stem-loop: GC-rich palindrome that can form a hairpin.
        # Simplified as runs of G/C flanked by complementary pairs.
        (
            r"(?:GC){3,}",
            "5' stem-loop structure (GC-rich hairpin proxy)",
            "Belasco & Higgins 1988",
        ),
        # Rho-independent terminator: GC hairpin followed by U-tract.
        # In DNA: GC-hairpin then T-tract (>=4 T's within 10 nt).
        (
            r"(?:GC){2,}.{0,8}(?:CG){2,}.{0,6}T{4,}",
            "Rho-independent terminator signal (GC-hairpin + T-tract)",
            "Carpousis 2010",
        ),
    ],
    "destabilizing": [
        # RNase E prefers single-stranded AU-rich regions.
        (
            r"[AT]{5,}",
            "RNase E recognition site (AT-rich single-stranded region)",
            "Carpousis 2010",
        ),
        # Consensus RNase E cleavage motif (AU-rich with positional bias).
        (
            r"A[AT]T[AT]A[AT]T",
            "RNase E cleavage consensus (AT-rich motif)",
            "Belasco & Higgins 1988",
        ),
        # Endonuclease cleavage site (A/U-rich with preference for A).
        (
            r"ATTA[AT]",
            "Endonuclease cleavage site (AUUA[AU] motif)",
            "Belasco & Higgins 1988",
        ),
    ],
}

_HUMAN_MOTIFS: dict[str, list[tuple[str, str, str]]] = {
    "stabilizing": [
        # GU-rich elements (GRE) — stabilizing in 3' UTR of some
        # short-lived transcripts when bound by CUG-BP1.
        (
            r"(?:GT){2,}",
            "GU-rich element (GRE)",
            "Barreau et al. 2005",
        ),
        # Stabilizing ARE variants — non-canonical ARE that recruit
        # stabilizing factors (e.g., HuR binding).
        (
            r"ATTTA.{1,10}ATTTA",
            "Stabilizing ARE variant (spaced AUUUA repeats - HuR binding)",
            "Chen & Shyu 1995",
        ),
        # PUF protein binding sites — typically UGUR motifs in 3' UTR.
        (
            r"TGT[ATCG]TAT",
            "PUF protein binding site (UGUR motif variant)",
            "Barreau et al. 2005",
        ),
    ],
    "destabilizing": [
        # ---- ARE Class I: scattered AUUUA in U-rich context ----
        (
            r"ATTTA",
            "ARE Class I (scattered AUUUA in U-rich context)",
            "Shaw & Kamen 1986; Chen & Shyu 1995",
        ),
        # ---- ARE Class II: overlapping AUUUA repeats ----
        (
            r"ATTTATTTA",
            "ARE Class II (overlapping AUUUA repeats - rapid decay)",
            "Chen & Shyu 1995",
        ),
        # ---- ARE Class III: U-rich without AUUUA ----
        (
            r"T{6,}",
            "ARE Class III (U-rich element without AUUUA)",
            "Chen & Shyu 1995; Barreau et al. 2005",
        ),
        # GU-UG repeats — associated with mRNA decay.
        (
            r"(?:GTTG){2,}",
            "GU-UG repeat (decay-promoting)",
            "Barreau et al. 2005",
        ),
    ],
}

_YEAST_MOTIFS: dict[str, list[tuple[str, str, str]]] = {
    "stabilizing": [
        # PGK1 stabilizing element — GC-rich region in the coding
        # sequence that protects against deadenylation-dependent decay.
        (
            r"(?:GC){3,}",
            "PGK1 stabilizing element (GC-rich coding region proxy)",
            "Muhlrad & Parker 1992",
        ),
        # C-type cyclin stabilizer — short motif found in CLN3 mRNA.
        (
            r"CCGCAC",
            "C-type cyclin stabilizer (CLN3-derived motif)",
            "Muhlrad & Parker 1992",
        ),
    ],
    "destabilizing": [
        # AUUUA-like motif (yeast variant — often in 3' UTR).
        (
            r"ATTTA",
            "AUUUA-like destabilizing motif",
            "Muhlrad & Parker 1992",
        ),
        # Rapid decay element — AU-rich stretch in 3' UTR.
        (
            r"[AT]{6,}",
            "Rapid decay element (AT-rich stretch)",
            "Muhlrad & Parker 1992",
        ),
    ],
}

# Public registry: organism name -> motif dict
STABILITY_MOTIFS: dict[str, dict[str, list[tuple[str, str, str]]]] = {
    "Escherichia_coli": _ECOLI_MOTIFS,
    "E_coli": _ECOLI_MOTIFS,
    "Homo_sapiens": _HUMAN_MOTIFS,
    "Mus_musculus": _HUMAN_MOTIFS,  # mammalian motifs shared
    "CHO_K1": _HUMAN_MOTIFS,
    "Saccharomyces_cerevisiae": _YEAST_MOTIFS,
}

# ────────────────────────────────────────────────────────────
# 2. Result dataclass
# ────────────────────────────────────────────────────────────


@dataclass
class MRNAStabilityScore:
    """Structured result of mRNA stability motif scanning.

    Attributes:
        overall_score: Stability score in [0, 1].  Higher is more stable.
            0.5 is neutral; <0.3 is high risk; >0.7 is low risk.
        stabilizing_count: Number of stabilizing motif hits found.
        destabilizing_count: Number of destabilizing motif hits found.
        motif_details: Per-hit information: position, matched text,
            effect (``"stabilizing"``/``"destabilizing"``), description.
        risk_level: ``"low"`` / ``"medium"`` / ``"high"`` bucketing
            based on the overall score.
    """

    overall_score: float
    stabilizing_count: int
    destabilizing_count: int
    motif_details: list[dict[str, Any]] = field(default_factory=list)
    risk_level: str = "medium"

    def __post_init__(self) -> None:
        # Clamp score
        self.overall_score = max(0.0, min(1.0, self.overall_score))
        # Derive risk level (consistent with predict_mrna_stability thresholds)
        # score > 0.7 → low risk (STABLE)
        # score 0.4–0.7 → medium risk (MODERATE)
        # score < 0.4 → high risk (UNSTABLE)
        if self.overall_score > 0.7:
            self.risk_level = "low"
        elif self.overall_score >= 0.4:
            self.risk_level = "medium"
        else:
            self.risk_level = "high"


# ────────────────────────────────────────────────────────────
# 3. Position helpers
# ────────────────────────────────────────────────────────────


def _fractional_position(pos: int, seq_len: int) -> str:
    """Classify a position as 5' UTR / CDS / 3' UTR proxy."""
    frac = pos / max(seq_len, 1)
    if frac < 0.10:
        return "5'UTR"
    elif frac > 0.85:
        return "3'UTR"
    return "CDS"


# Position-dependent weight multipliers.
# Destabilizing motifs in 3' UTR are especially harmful; stabilizing
# motifs in 5' UTR are especially helpful.

_POSITION_WEIGHTS: dict[str, dict[str, float]] = {
    "stabilizing": {
        "5'UTR": 1.5,
        "CDS": 1.0,
        "3'UTR": 1.2,
    },
    "destabilizing": {
        "5'UTR": 0.8,
        "CDS": 1.0,
        "3'UTR": 1.8,
    },
}


# ────────────────────────────────────────────────────────────
# 4. score_mrna_stability
# ────────────────────────────────────────────────────────────


# ────────────────────────────────────────────────────────────
# Precompiled regex patterns (Optimization: avoid re.compile per call)
# ────────────────────────────────────────────────────────────

_COMPILED_MOTIFS: dict[str, dict[str, list[tuple[Any, str, str]]]] = {}


def _get_compiled_motifs(organism: str) -> dict[str, list[tuple[Any, str, str]]]:
    """Get or build precompiled motif patterns for the given organism.

    Caches compiled regex patterns so they are compiled only once per
    organism, not once per call to score_mrna_stability.
    """
    if organism in _COMPILED_MOTIFS:
        return _COMPILED_MOTIFS[organism]

    raw = STABILITY_MOTIFS.get(organism, STABILITY_MOTIFS.get("Homo_sapiens", {}))
    compiled: dict[str, list[tuple[Any, str, str]]] = {"stabilizing": [], "destabilizing": []}

    for effect in ("stabilizing", "destabilizing"):
        for pattern, description, source in raw.get(effect, []):
            try:
                compiled_pattern = re.compile(pattern, re.IGNORECASE)
            except re.error as exc:
                logger.debug("Skipping invalid regex '%s': %s", pattern, exc)
                continue
            compiled[effect].append((compiled_pattern, description, source))

    _COMPILED_MOTIFS[organism] = compiled
    return compiled


def score_mrna_stability(dna: str, organism: str) -> MRNAStabilityScore:
    """Scan *dna* for stability motifs and return a composite score.

    The algorithm combines two complementary models:

    **Codon optimality (CAI-based)**: primary component
      Computes the Codon Adaptation Index (CAI) as a proxy for mRNA
      stability.  Optimal codons (high adaptiveness) promote efficient
      translation and mRNA stability; non-optimal codons accelerate
      decay (Presnyak et al. 2015).  CAI > 0.8 → likely stable;
      CAI < 0.5 → likely unstable.

    **Motif-based scoring**: modifier component
      1. Look up the organism in ``STABILITY_MOTIFS`` (falls back to
         *Homo_sapiens* with a warning).
      2. Scan for each motif using ``re.finditer`` (case-insensitive).
      3. Weight each hit by its position (5' UTR / CDS / 3' UTR proxy).
      4. Destabilizing motifs reduce the CAI-based score; stabilizing
         motifs boost it slightly.

    The final composite score is:

        ``overall = CAI - beta * sum(destabilizing_weights)
                         + alpha * sum(stabilizing_weights)``

    This integration ensures that sequences with good codon optimality
    receive high stability scores, while sequences with destabilizing
    motifs are appropriately penalised.

    Args:
        dna: DNA coding-strand sequence (5'->3').
        organism: Organism name (e.g. ``"Homo_sapiens"``).

    Returns:
        ``MRNAStabilityScore`` with score, counts, and per-motif details.
    """
    dna = dna.upper()
    seq_len = len(dna)

    if seq_len == 0:
        return MRNAStabilityScore(
            overall_score=0.5,
            stabilizing_count=0,
            destabilizing_count=0,
        )

    # ---- Codon optimality (CAI-based) component ----
    cai_score = _compute_cai(dna, organism)

    # ---- Motif-based component ----
    # Resolve motif database (use precompiled patterns for performance)
    compiled_motifs = _get_compiled_motifs(organism)

    details: list[dict[str, Any]] = []
    stab_weighted_sum = 0.0
    destab_weighted_sum = 0.0

    for effect in ("stabilizing", "destabilizing"):
        for compiled, description, source in compiled_motifs.get(effect, []):
            for m in compiled.finditer(dna):
                pos = m.start()
                region = _fractional_position(pos, seq_len)
                weight = _POSITION_WEIGHTS[effect][region]

                detail = {
                    "position": pos,
                    "region": region,
                    "motif": m.group(),
                    "effect": effect,
                    "description": description,
                    "source": source,
                    "weight": weight,
                }
                details.append(detail)

                if effect == "stabilizing":
                    stab_weighted_sum += weight
                else:
                    destab_weighted_sum += weight

    stab_count = sum(1 for d in details if d["effect"] == "stabilizing")
    destab_count = sum(1 for d in details if d["effect"] == "destabilizing")

    # -- Motif scoring formula --
    # Base score is 0.5 (neutral).  Each stabilizing weighted hit adds
    # alpha per hit; each destabilizing weighted hit subtracts beta per hit.
    _ALPHA = 0.04  # stabilizing contribution per weighted hit
    _BETA = 0.05   # destabilizing penalty per weighted hit

    motif_score = 0.5 + _ALPHA * stab_weighted_sum - _BETA * destab_weighted_sum

    # -- Density penalty --
    # If destabilizing motifs exceed 3 per kb, apply an extra penalty.
    destab_per_kb = destab_count / (seq_len / 1000.0) if seq_len > 0 else 0.0
    _DESTAB_DENSITY_THRESHOLD = 3.0
    _DESTAB_DENSITY_PENALTY = 0.03  # per extra hit/kb

    if destab_per_kb > _DESTAB_DENSITY_THRESHOLD:
        excess = destab_per_kb - _DESTAB_DENSITY_THRESHOLD
        motif_score -= _DESTAB_DENSITY_PENALTY * excess

    # ---- Composite score: CAI primary + motif adjustments ----
    # CAI serves as the primary stability proxy (Presnyak et al. 2015).
    # Motif scores act as modifiers: destabilizing motifs reduce the score,
    # stabilizing motifs boost it slightly.
    score = cai_score - _BETA * destab_weighted_sum + _ALPHA * stab_weighted_sum

    # ---- ATTTA one-level downgrade for consistency with predict_mrna_stability ----
    # When ATTTA motifs are present for organisms that check them,
    # ensure the score is downgraded by one category level so that:
    #   score > 0.7 → STABLE, 0.4–0.7 → MODERATE, < 0.4 → UNSTABLE
    # matches what predict_mrna_stability would return.
    predict_thresholds = _ORGANISM_STABILITY_THRESHOLDS.get(
        organism,
        _ORGANISM_STABILITY_THRESHOLDS.get(
            "Homo_sapiens", _DEFAULT_STABILITY_THRESHOLD
        ),
    )
    check_motifs = predict_thresholds[2]

    if check_motifs and _has_instability_motif(dna):
        if score > 0.7:
            # STABLE → MODERATE: shift score down by 0.3
            score -= 0.3
        elif score >= 0.4:
            # MODERATE → UNSTABLE: shift score down by 0.4
            score -= 0.4

    return MRNAStabilityScore(
        overall_score=score,
        stabilizing_count=stab_count,
        destabilizing_count=destab_count,
        motif_details=details,
    )


# ────────────────────────────────────────────────────────────
# 5. Codon optimality model (CAI-based stability scoring)
# ────────────────────────────────────────────────────────────

# Codon optimality coefficients: derived from CAI adaptiveness tables.
# These represent the "codon optimality" concept from Presnyak et al. (2015):
# optimal codons → stable mRNA, non-optimal codons → unstable mRNA.
#
# The values are the organism-specific CAI adaptiveness weights, which
# range from 0 (rare/non-optimal) to 1.0 (most frequent/optimal).
# They are loaded lazily from the organisms module.

_CODON_OPTIMALITY: dict[str, dict[str, float]] = {}


def _get_codon_optimality(organism: str) -> dict[str, float]:
    """Get codon optimality coefficients for the given organism.

    Returns the CAI adaptiveness table for the organism, which serves
    as a proxy for codon optimality.  Values are cached after first load.
    """
    if organism in _CODON_OPTIMALITY:
        return _CODON_OPTIMALITY[organism]

    try:
        from .organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism

        canonical = resolve_organism(organism)
        table = CODON_ADAPTIVENESS_TABLES.get(canonical, {})
    except (ImportError, AttributeError):
        logger.debug("Could not load CAI adaptiveness for organism '%s'", organism)
        table = {}

    _CODON_OPTIMALITY[organism] = table
    return table


# ────────────────────────────────────────────────────────────
# 6. Organism-specific stability thresholds
# ────────────────────────────────────────────────────────────

# Each entry: (stable_threshold, unstable_threshold, check_instability_motifs)
# - stable_threshold: CAI/half-life score above this → STABLE
# - unstable_threshold: CAI/half-life score below this → UNSTABLE
# - check_instability_motifs: if True, ATTTA motifs override CAI to UNSTABLE
#
# Rationale:
#   E. coli: mRNA half-life typically 2-8 min; CAI > 0.7 predicts stability
#     (Ikemura 1985; Sharp & Li 1987)
#   Human: mRNA half-life varies widely; CAI > 0.7 AND no ATTTA instability
#     motifs predict stability (Shaw & Kamen 1986; Presnyak et al. 2015)
#   Yeast: similar to human; CAI > 0.7 AND no ATTTA (Muhlrad & Parker 1992)

_ORGANISM_STABILITY_THRESHOLDS: dict[str, tuple[float, float, bool]] = {
    "Escherichia_coli": (0.7, 0.5, False),  # E. coli: no ATTTA check needed
    "E_coli": (0.7, 0.5, False),            # alias
    "Homo_sapiens": (0.7, 0.5, True),        # Human: check ATTTA
    "Mus_musculus": (0.7, 0.5, True),        # Mouse: same as human
    "CHO_K1": (0.7, 0.5, True),             # CHO: same as human
    "Saccharomyces_cerevisiae": (0.7, 0.5, True),  # Yeast: check ATTTA
}

# Default thresholds for unrecognised organisms
_DEFAULT_STABILITY_THRESHOLD: tuple[float, float, bool] = (0.7, 0.5, True)

# ATTTA instability motif (DNA form of AUUUA)
_INSTABILITY_MOTIF_ATTTA = "ATTTA"


def _has_instability_motif(dna: str) -> bool:
    """Check if the DNA sequence contains the ATTTA instability motif."""
    return _INSTABILITY_MOTIF_ATTTA in dna.upper()


def _compute_cai(dna: str, organism: str) -> float:
    """Compute CAI for a DNA sequence using the organism's adaptiveness table.

    This is a lightweight internal CAI computation that doesn't depend on
    the translation module, avoiding circular imports.

    Returns:
        CAI value in [0, 1]. Returns 0.0 for empty or invalid sequences.
    """
    dna = dna.upper()
    if len(dna) < 3:
        return 0.0

    adaptiveness = _get_codon_optimality(organism)
    if not adaptiveness:
        return 0.0

    _EPSILON = 1e-10
    log_sum = 0.0
    count = 0

    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa = CODON_TABLE.get(codon)
        # Skip Met (ATG, only codon) and stop codons — Sharp & Li (1987)
        if aa is None or aa == "*" or aa == "M":
            continue
        w = adaptiveness.get(codon, 0.0)
        if w <= 0:
            w = _EPSILON
        log_sum += math.log(w)
        count += 1

    if count == 0:
        return 0.0

    return math.exp(log_sum / count)


def compute_mrna_half_life_score(dna_sequence: str, organism: str) -> float:
    """Compute a normalised mRNA half-life score using codon optimality.

    This function predicts mRNA stability primarily through codon optimality,
    using CAI as a proxy.  The score reflects the expected half-life of the
    mRNA relative to organism-specific baselines:

    - **E. coli**: baseline half-life 2-8 min; optimal codons extend half-life
      up to ~20 min (BP; Belasco 1993).
    - **Human**: half-life varies from minutes to hours; optimal codons and
      absence of AREs promote stability.
    - **Yeast**: similar to human, with additional sensitivity to AT-rich
      decay elements.

    The model combines:

      1. **CAI score** (primary): directly reflects codon optimality.
         CAI > 0.8 → likely stable; CAI < 0.5 → likely unstable.
      2. **Motif penalties**: ATTTA instability motifs reduce the score.
      3. **GC content**: extreme GC (<30% or >70%) penalises stability.

    Args:
        dna_sequence: DNA coding-strand sequence (5'->3').
        organism: Organism name (e.g. ``"Homo_sapiens"``, ``"Escherichia_coli"``).

    Returns:
        Float in [0, 1] where:
          - > 0.7 → likely stable (long half-life)
          - 0.5 – 0.7 → moderate stability
          - < 0.5 → likely unstable (short half-life)
    """
    dna = dna_sequence.upper()
    seq_len = len(dna)

    if seq_len < 3:
        return 0.5  # insufficient data → neutral

    # --- 1. CAI-based codon optimality score (primary) ---
    cai = _compute_cai(dna, organism)

    # CAI directly serves as the base stability score.
    # CAI > 0.8 → likely stable; CAI < 0.5 → likely unstable
    # (Presnyak et al. 2015: codon optimality is a major determinant
    #  of mRNA stability).
    base_score = cai

    # --- 2. Instability motif penalty ---
    motif_penalty = 0.0
    thresholds = _ORGANISM_STABILITY_THRESHOLDS.get(
        organism,
        _ORGANISM_STABILITY_THRESHOLDS.get(
            "Homo_sapiens", _DEFAULT_STABILITY_THRESHOLD
        ),
    )
    check_motifs = thresholds[2]

    if check_motifs:
        # Count ATTTA motifs — each one reduces stability
        atttta_count = dna.count(_INSTABILITY_MOTIF_ATTTA)
        if atttta_count > 0:
            # Each ATTTA motif reduces score by 0.05 (capped at 0.30)
            motif_penalty = min(0.30, atttta_count * 0.05)

    # For E. coli, check RNase E sites instead
    if organism in ("Escherichia_coli", "E_coli"):
        # AT-rich runs of 5+ are RNase E targets
        at_rich_runs = len(re.findall(r"[AT]{5,}", dna))
        if at_rich_runs > 0:
            motif_penalty = min(0.25, at_rich_runs * 0.04)

    # --- 3. GC content adjustment ---
    gc_contribution = 0.0
    if seq_len > 0:
        gc = (dna.count("G") + dna.count("C")) / seq_len
        # Optimal GC range: 40-60%. Outside this range, penalise.
        if gc < 0.30:
            gc_contribution = -0.10 * (0.30 - gc) / 0.30  # up to -0.10
        elif gc > 0.70:
            gc_contribution = -0.10 * (gc - 0.70) / 0.30  # up to -0.10
        else:
            gc_contribution = 0.0  # within optimal range

    # --- Composite score ---
    # CAI directly serves as the stability proxy (Presnyak et al. 2015).
    # Motif penalties and GC adjustments are subtracted from the CAI base.
    # This ensures CAI > 0.7 → score > 0.7 → STABLE when no motifs present.
    score = base_score - motif_penalty + gc_contribution

    return max(0.0, min(1.0, score))


def predict_mrna_stability(
    dna_sequence: str,
    organism: str,
) -> StabilityCategory:
    """Predict mRNA stability category for a DNA sequence.

    Uses the Codon Adaptation Index (CAI) directly to classify mRNA
    stability into one of three categories, with an ATTTA motif
    penalty that downgrades by at most one level:

      - **STABLE**: CAI >= 0.8 and no ATTTA motifs (or organism
        doesn't check them).  The mRNA is expected to have a long
        half-life suitable for high expression.
      - **MODERATE**: CAI 0.5–0.8, or CAI >= 0.8 with ATTTA motifs
        (one-level downgrade from STABLE).
      - **UNSTABLE**: CAI < 0.5, or CAI 0.5–0.8 with ATTTA motifs
        (one-level downgrade from MODERATE).

    Organism-specific rules:

      - *E. coli*: CAI-based only (no ATTTA check needed)
      - *Human / mammalian / yeast*: ATTTA motif presence
        downgrades the CAI-based category by one level

    Args:
        dna_sequence: DNA coding-strand sequence (5'->3').
        organism: Organism name (e.g. ``"Homo_sapiens"``).

    Returns:
        One of ``"STABLE"``, ``"MODERATE"``, or ``"UNSTABLE"``.
    """
    dna = dna_sequence.upper()
    cai = _compute_cai(dna, organism)

    # Step 1: Map CAI to stability category
    if cai >= 0.8:
        category: StabilityCategory = "STABLE"
    elif cai >= 0.5:
        category = "MODERATE"
    else:
        category = "UNSTABLE"

    # Step 2: ATTTA penalty — downgrade by one level only
    thresholds = _ORGANISM_STABILITY_THRESHOLDS.get(
        organism,
        _ORGANISM_STABILITY_THRESHOLDS.get(
            "Homo_sapiens", _DEFAULT_STABILITY_THRESHOLD
        ),
    )
    check_motifs = thresholds[2]

    if check_motifs and _has_instability_motif(dna):
        if category == "STABLE":
            category = "MODERATE"
        elif category == "MODERATE":
            category = "UNSTABLE"

    return category


# ────────────────────────────────────────────────────────────
# 7. suggest_mutations_for_stability
# ────────────────────────────────────────────────────────────

# Build reverse lookup: codon -> amino acid
_CODON_TO_AA: dict[str, str] = {}
for _codon, _aa in CODON_TABLE.items():
    if _aa != "*":
        _CODON_TO_AA[_codon] = _aa


def suggest_mutations_for_stability(
    dna: str,
    organism: str,
) -> list[dict[str, Any]]:
    """Identify destabilizing motifs and propose synonymous codon changes.

    For each destabilizing motif hit found by ``score_mrna_stability``,
    the algorithm:
      1. Determines the codons overlapping the motif.
      2. For each overlapping codon, looks up alternative (synonymous)
         codons that encode the same amino acid.
      3. Selects the alternative whose substitution breaks the motif
         pattern (i.e., the motif regex no longer matches the modified
         sequence).
      4. Prefers higher-adaptiveness codons (organism CAI rank) as
         tie-breaker.

    Args:
        dna: DNA coding-strand sequence (5'->3').
        organism: Organism name.

    Returns:
        List of dicts, each with keys:
          - ``position``: 0-based position of the suggested codon change
          - ``original_codon``: the current codon
          - ``suggested_codon``: the replacement codon
          - ``amino_acid``: the amino acid (unchanged)
          - ``motif_removed``: the destabilizing motif that would be broken
          - ``motif_position``: start position of the removed motif
    """
    dna = dna.upper()
    seq_len = len(dna)

    if seq_len < 3:
        return []

    # Get the stability report to find destabilizing hits
    report = score_mrna_stability(dna, organism)
    destab_hits = [d for d in report.motif_details if d["effect"] == "destabilizing"]

    if not destab_hits:
        return []

    # Try to load organism-specific CAI weights for codon ranking
    cai_weights: dict[str, float] = {}
    try:
        from .organisms import CODON_ADAPTIVENESS_TABLES

        org_key = organism
        if org_key not in CODON_ADAPTIVENESS_TABLES:
            # Try common aliases
            if org_key in ("E_coli",):
                org_key = "Escherichia_coli"
            elif org_key in ("human",):
                org_key = "Homo_sapiens"
            elif org_key in ("mouse",):
                org_key = "Mus_musculus"
            elif org_key in ("yeast",):
                org_key = "Saccharomyces_cerevisiae"
            elif org_key in ("cho", "CHO_K1"):
                org_key = "CHO_K1"
        cai_weights = CODON_ADAPTIVENESS_TABLES.get(org_key, {})
    except (ImportError, AttributeError):
        logger.debug("Could not load CAI weights for organism '%s'", organism)

    suggestions: list[dict[str, Any]] = []
    seen_positions: set[int] = set()  # avoid duplicate suggestions

    for hit in destab_hits:
        motif_start = hit["position"]
        motif_end = motif_start + len(hit["motif"])
        # Re-derive the regex pattern from the motif database
        motif_pattern = _find_pattern_for_description(organism, hit["description"])

        # Find codons overlapping the motif
        # Codons are on frame boundaries starting from position 0
        first_codon_start = (motif_start // 3) * 3
        last_codon_start = ((motif_end - 1) // 3) * 3

        for codon_start in range(first_codon_start, last_codon_start + 1, 3):
            if codon_start in seen_positions:
                continue
            if codon_start + 3 > seq_len:
                continue

            original_codon = dna[codon_start : codon_start + 3]
            aa = _CODON_TO_AA.get(original_codon)
            if aa is None:
                continue

            # Get synonymous codons
            synonyms = AA_TO_CODONS.get(aa, [])
            if len(synonyms) <= 1:
                continue  # no alternative

            # Try each synonym and check if it breaks the motif
            best_candidate: dict[str, Any] | None = None
            best_cai = -1.0

            for syn in synonyms:
                if syn == original_codon:
                    continue

                # Construct the modified sequence
                modified = dna[:codon_start] + syn + dna[codon_start + 3 :]

                # Check if the motif still matches near this position
                if motif_pattern is not None:
                    try:
                        search_start = max(0, motif_start - 3)
                        search_end = min(len(modified), motif_end + 3)
                        if re.search(
                            motif_pattern,
                            modified[search_start:search_end],
                            re.IGNORECASE,
                        ):
                            # Motif still present — this synonym does not help
                            continue
                    except re.error:
                        continue
                else:
                    # No pattern to test — check if the motif text is altered
                    modified_motif = modified[motif_start : min(motif_end, len(modified))]
                    if modified_motif == hit["motif"]:
                        continue

                # Motif broken — rank by CAI weight
                cai = cai_weights.get(syn, 0.5)
                if cai > best_cai:
                    best_cai = cai
                    best_candidate = {
                        "position": codon_start,
                        "original_codon": original_codon,
                        "suggested_codon": syn,
                        "amino_acid": aa,
                        "motif_removed": hit["motif"],
                        "motif_position": motif_start,
                    }

            if best_candidate is not None:
                suggestions.append(best_candidate)
                seen_positions.add(codon_start)

    return suggestions


def _find_pattern_for_description(
    organism: str,
    description: str,
) -> str | None:
    """Look up the regex pattern for a given motif description and organism."""
    motifs = STABILITY_MOTIFS.get(organism, STABILITY_MOTIFS.get("Homo_sapiens", {}))
    for effect in ("stabilizing", "destabilizing"):
        for pattern, desc, _source in motifs.get(effect, []):
            if desc == description:
                return pattern
    return None
