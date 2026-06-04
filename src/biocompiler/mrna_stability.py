"""
BioCompiler mRNA Stability Motif Models
========================================

Organism-specific models for scoring mRNA stability based on
cis-regulatory sequence motifs in the coding and UTR regions.

mRNA stability is a critical determinant of heterologous protein expression.
Unstable transcripts — often caused by AU-rich elements (AREs), RNase
recognition sites, or other decay-promoting motifs — can severely limit
yield even when codon usage is fully optimized.  This module scans DNA
sequences for organism-specific stabilizing and destabilizing motifs,
computes a composite stability score, and proposes synonymous codon
substitutions that disrupt destabilizing motifs while preserving amino
acid identity and CAI.

Three organism models are provided:
  - *E. coli*:  RNase E sites, stem-loop stabilizers, Rho-independent
    terminators (Belasco & Higgins 1988; Carpousis 2010).
  - *Human / mammalian*:  ARE subtypes (Shaw & Kamen 1986), GRE,
    PUF-binding sites (Chen & Shyu 1995; Barreau et al. 2005).
  - *S. cerevisiae*:  PGK1 stabilizing element, rapid-decay elements
    (Muhlrad & Parker 1992).

Usage::

    from biocompiler.mrna_stability import score_mrna_stability, suggest_mutations_for_stability

    # Score mRNA stability of a coding sequence
    dna = "ATGGCCGCGATTTAGCGGCGCCTAA"
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
- ``score_mrna_stability()`` — scan a DNA sequence and return a score
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
import re
from dataclasses import dataclass, field
from typing import Any

from .constants import AA_TO_CODONS, CODON_TABLE

logger = logging.getLogger(__name__)

__all__ = [
    "STABILITY_MOTIFS",
    "MRNAStabilityScore",
    "score_mrna_stability",
    "suggest_mutations_for_stability",
]

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
        # Derive risk level
        if self.overall_score >= 0.7:
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


def score_mrna_stability(dna: str, organism: str) -> MRNAStabilityScore:
    """Scan *dna* for stability motifs and return a composite score.

    The algorithm:
      1. Look up the organism in ``STABILITY_MOTIFS`` (falls back to
         *Homo_sapiens* with a warning).
      2. Scan for each motif using ``re.finditer`` (case-insensitive).
      3. Weight each hit by its position (5' UTR / CDS / 3' UTR proxy).
      4. Compute a score:

         ``score = 0.5 + alpha * sum(stabilizing_weights)
                        - beta  * sum(destabilizing_weights)``

         where alpha and beta are calibrated so that a typical CDS
         without extreme motif density lands near 0.5.

      5. Density penalty: if destabilizing hits per kb exceeds a
         threshold, an additional penalty is applied.

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

    # Resolve motif database
    motifs = STABILITY_MOTIFS.get(organism)
    if motifs is None:
        logger.warning(
            "No stability motifs for organism '%s'; "
            "falling back to Homo_sapiens",
            organism,
        )
        motifs = STABILITY_MOTIFS["Homo_sapiens"]

    details: list[dict[str, Any]] = []
    stab_weighted_sum = 0.0
    destab_weighted_sum = 0.0

    for effect in ("stabilizing", "destabilizing"):
        for pattern, description, source in motifs.get(effect, []):
            try:
                compiled = re.compile(pattern, re.IGNORECASE)
            except re.error as exc:
                logger.debug("Skipping invalid regex '%s': %s", pattern, exc)
                continue

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

    # -- Scoring formula --
    # Base score is 0.5 (neutral).  Each stabilizing weighted hit adds
    # alpha per hit; each destabilizing weighted hit subtracts beta per hit.
    # alpha, beta are calibrated so that moderate motif counts stay in [0.3, 0.7].
    _ALPHA = 0.04  # stabilizing contribution per weighted hit
    _BETA = 0.05   # destabilizing penalty per weighted hit

    score = 0.5 + _ALPHA * stab_weighted_sum - _BETA * destab_weighted_sum

    # -- Density penalty --
    # If destabilizing motifs exceed 3 per kb, apply an extra penalty.
    destab_per_kb = destab_count / (seq_len / 1000.0) if seq_len > 0 else 0.0
    _DESTAB_DENSITY_THRESHOLD = 3.0
    _DESTAB_DENSITY_PENALTY = 0.03  # per extra hit/kb

    if destab_per_kb > _DESTAB_DENSITY_THRESHOLD:
        excess = destab_per_kb - _DESTAB_DENSITY_THRESHOLD
        score -= _DESTAB_DENSITY_PENALTY * excess

    return MRNAStabilityScore(
        overall_score=score,
        stabilizing_count=stab_count,
        destabilizing_count=destab_count,
        motif_details=details,
    )


# ────────────────────────────────────────────────────────────
# 5. suggest_mutations_for_stability
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
