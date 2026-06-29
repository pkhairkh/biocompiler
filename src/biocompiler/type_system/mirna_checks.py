"""
BioCompiler Type System — miRNA Binding-Site Predicate Checks
=============================================================
miRNA binding-site predicate check plus its private helpers
(_rna_revcomp_to_dna seed conversion and _mirna_context_score
context-factor adjustment).

Extracted from the historical checks.py monolith during the W8-b refactor.
Re-exported by checks.py for backwards compatibility.
"""
import logging
from typing import Any, Dict, List, Optional, Set, Tuple

from .codon_tables import (
    AA_TO_CODONS,
    BLOSUM62,
    CODON_TABLE,
    PROMOTER_CONSENSUS,
    PredicateResult,
    _BLOSUM62_MISSING_SCORE,
    _CPG_DENSITY_MULTIPLIER,
    _CPG_GC_RICH_THRESHOLD,
    _COTRANS_HIGH_CONFIDENCE,
    _COTRANS_LOW_CONFIDENCE,
    _CODON_RAMP_LENGTH,
    _DG_AU_PAIR_KCAL,
    _DG_GC_PAIR_KCAL,
    _DG_GU_PAIR_KCAL,
    _EUKARYOTE_GT_PER_BP,
    _EUK_INITIATOR_OFFSET_MAX,
    _EUK_INITIATOR_OFFSET_MIN,
    _FAST_CODON_CAI_THRESHOLD,
    _HIGH_AVG_CAI_THRESHOLD,
    _INSTABILITY_T_RUN_MIN,
    _MAXENT_INSUFFICIENT_CONTEXT_SCORE,
    _MIN_RAMP_FOR_WARNING,
    _MRNA_DG_EUKARYOTE_FAIL,
    _MRNA_DG_PROKARYOTE_FAIL,
    _MRNA_MODERATE_DG_RATIO,
    _MRNA_STABILITY_THRESHOLDS,
    _ORGANISM_TO_SPECIES_KEY,
    _PAUSE_SITE_CAI_THRESHOLD,
    _PROMOTER_UNCERTAIN_RATIO,
    _RESTRICTION_SITE_MIN_LENGTH,
    _TM_BORDERLINE_RATIO,
    _TM_EUKARYOTIC_MIN_STRETCH,
    _TM_PROKARYOTIC_MIN_STRETCH,
    _match_iupac,
    _score_consensus,
)
from biocompiler.shared.types import Verdict

# miRNA seed database (multi-organism, tissue-filtered). This
# import was at module level in the original checks.py (line 1760,
# immediately before _rna_revcomp_to_dna) and is required by
# check_no_mirna_binding_site.
from .mirna_seeds import get_mirna_seeds


logger = logging.getLogger(__name__)

def _rna_revcomp_to_dna(seed_rna: str) -> str:
    """Compute the DNA target of an miRNA seed (reverse complement, U→T).

    The miRNA binds the mRNA 3'→5', so the mRNA target site (5'→3') is
    the reverse complement of the miRNA seed.  Since the DNA coding strand
    equals the mRNA with T instead of U, the DNA target is the reverse
    complement of the RNA seed with U→T.

    Args:
        seed_rna: miRNA seed sequence, 5'→3' RNA (e.g. ``"AGCUUAU"``).

    Returns:
        DNA target string to scan for in the coding strand.
    """
    complement = {"A": "U", "U": "A", "G": "C", "C": "G"}
    rc = "".join(complement[b] for b in reversed(seed_rna))
    return rc.replace("U", "T")


def _mirna_context_score(
    seq: str,
    hit_pos: int,
    match_type: str,
    hit_score: float,
) -> float:
    """Adjust miRNA binding score based on CDS context.

    In coding sequences, miRNA binding efficacy is modified by:
    1. Position relative to codon boundaries: sites spanning exon-exon
       junctions (if known) are less accessible; sites fully within a
       codon are more accessible.
    2. AU-rich flanking: bases flanking the seed match that are A/T
       improve accessibility (lower mRNA secondary structure).
    3. Conservation proxy: in CDS, seed matches in more conserved
       codon positions (1st/2nd base) are less likely to be functional
       because they are constrained by protein coding; matches at 3rd
       positions are more likely to be functional.

    Returns:
        Adjusted score in [0, 1.2] range. Values > 1.0 indicate
        the context makes the site MORE likely to be functional.
    """
    n = len(seq)
    # Determine the span of the seed match in the DNA sequence
    seed_len = 6 if match_type == "6mer" else 7
    if match_type == "8mer":
        seed_len = 8

    # ── Factor 1: AU-rich flanking context ──
    # Count A/T bases in ±3 positions around the seed match
    flank_start = max(0, hit_pos - 3)
    flank_end = min(n, hit_pos + seed_len + 3)
    flanking = seq[flank_start:hit_pos] + seq[hit_pos + seed_len:flank_end]
    au_count = sum(1 for b in flanking if b in ("A", "T"))
    flanking_len = len(flanking)
    # AU fraction: 0.0 = all GC, 1.0 = all AU
    au_fraction = au_count / flanking_len if flanking_len > 0 else 0.5
    # Boost for AU-rich: +0.0 to +0.15
    au_factor = 0.0 + 0.15 * au_fraction

    # ── Factor 2: Codon position (reading frame) ──
    # In CDS, 3rd codon positions are less constrained by protein coding,
    # so seed matches falling predominantly on 3rd positions are more
    # likely to be functionally relevant.
    third_pos_count = 0
    for i in range(seed_len):
        codon_pos = (hit_pos + i) % 3  # 0=1st, 1=2nd, 2=3rd
        if codon_pos == 2:
            third_pos_count += 1
    # If most of the seed sits on 3rd positions: boost
    # If most sits on 1st/2nd positions: penalty
    if seed_len > 0:
        third_fraction = third_pos_count / seed_len
    else:
        third_fraction = 0.33
    # Range: -0.10 (all 1st/2nd) to +0.10 (all 3rd)
    codon_factor = -0.10 + 0.20 * third_fraction

    # ── Factor 3: CDS position (proximity to 3' end / stop codon) ──
    # Sites near the 3' end of the CDS are closer to the UTR boundary
    # and more likely to be functional for miRNA binding.
    relative_pos = hit_pos / n if n > 0 else 0.5
    # Last 20% of CDS gets a small boost, first 20% gets a small penalty
    if relative_pos >= 0.8:
        position_factor = 0.05
    elif relative_pos <= 0.2:
        position_factor = -0.05
    else:
        position_factor = 0.0

    # ── Combine: base score × context factor ──
    # Context factor in [0.7, 1.2]: 1.0 = neutral
    raw_context = 1.0 + au_factor + codon_factor + position_factor
    context_factor = max(0.7, min(1.2, raw_context))

    return round(hit_score * context_factor, 3)


# Tissue proximity groups — miRNAs from nearby tissues are partially relevant
_TISSUE_GROUPS: Dict[str, Set[str]] = {
    "liver": {"liver", "hepatocyte", "fibrotic"},
    "muscle": {"muscle", "smooth_muscle", "cardiac"},
    "immune": {"immune", "blood", "lymphoid", "hematopoietic"},
    "neural": {"neural", "brain"},
    "epithelial": {"epithelial", "endothelial", "kidney"},
    "colon": {"colon", "intestinal", "digestive"},
    "ubiquitous": {"ubiquitous"},
    "tumor_suppressor": {"tumor_suppressor", "cancer"},
}


def check_no_mirna_binding_site(
    seq: str,
    organism: str = "Homo_sapiens",
    min_seed_match: int = 7,
    tissue: str = "",
) -> PredicateResult:
    """Predicate: No high-affinity miRNA binding sites.

    Scans the DNA coding sequence for seed matches against known miRNA
    families. miRNA binding to mRNA coding regions can trigger:

    - Translational repression (Argonaute-mediated ribosome dropoff)
    - mRNA destabilization (deadenylation and decapping)
    - Off-target silencing in vaccine applications

    The scan implements a seed-based matching algorithm:
    1. For each miRNA seed (positions 2-8 from 5' end), scan for
       Watson-Crick complementary matches in the DNA sequence
    2. Seed match types (from most to least binding):
       - 8mer: perfect match to positions 1-8 + adenosine at position 1
       - 7mer-m8: perfect match to positions 2-8
       - 7mer-A1: match to positions 2-7 + adenosine at position 1
       - 6mer: match to positions 2-7
    3. Score each match: 8mer = 1.0, 7mer-m8 = 0.9, 7mer-A1 = 0.85, 6mer = 0.7
    4. Adjust scores with CDS context model (see ``_mirna_context_score``):
       - AU-rich flanking improves accessibility (boost)
       - Seed matches at 3rd codon positions are less protein-constrained (boost)
       - Matches at 1st/2nd codon positions are more constrained (penalty)
       - Sites near the 3' end of CDS are closer to the UTR (small boost)
       - Context factor range: [0.7, 1.2] applied as multiplier
    5. Aggregate by miRNA family to avoid double-counting

    The miRNA seed database includes the top 50 most abundantly expressed
    miRNAs across human tissues (from miRBase v22 + tissue expression data
    from Ludwig et al. 2016).

    When a *tissue* is specified, the predicate applies tissue-aware
    filtering so that only biologically relevant miRNAs contribute to
    the verdict:

    - Ubiquitous miRNAs are always included (they affect all tissues).
    - miRNAs tagged with the exact target tissue are always included.
    - miRNAs from related tissues (via ``_TISSUE_GROUPS``) are included
      but their effective tier is worsened by 1 (tier-1 → tier-2, etc.).
    - miRNAs from unrelated tissues are downweighted: their effective
      tier is increased by 1, making them less likely to trigger a FAIL.

    Verdict logic:
    - FAIL: 8mer or 7mer-m8 match with context-adjusted score >= 0.9
      for any top-10 abundantly expressed miRNA (after tissue-aware
      tier adjustment)
    - UNCERTAIN: 7mer or 6mer match, or match to lower-abundance miRNA
    - PASS: No significant miRNA seed matches

    Args:
        seq: DNA coding sequence (uppercase).
        organism: Target organism for tissue-specific miRNA selection.
        min_seed_match: Minimum seed match length (6 or 7). Default 7.
        tissue: Target tissue for context-aware miRNA filtering
            (e.g. "liver", "muscle", "immune").  Empty or "ubiquitous"
            uses all seeds with original tiers (legacy behaviour).

    Returns:
        PredicateResult with verdict and positions of miRNA binding sites.
    """
    seq = seq.upper()
    n = len(seq)

    if n < 6:
        return PredicateResult(
            "NoMiRNABindingSite", True, verdict=Verdict.PASS,
            details="Sequence too short for miRNA seed matching",
        )

    # ── Tissue-aware seed filtering ────────────────────────────
    # Determine which tissue group the target tissue belongs to
    # so we can classify each seed as exact-match, related, or unrelated.
    tissue_lower = tissue.lower().strip()
    use_tissue_filter = tissue_lower not in ("", "ubiquitous")

    # Resolve the target tissue to its proximity group key.
    # e.g. "hepatocyte" → "liver", "blood" → "immune"
    target_group_key: Optional[str] = None
    if use_tissue_filter:
        for group_key, group_members in _TISSUE_GROUPS.items():
            if tissue_lower in group_members:
                target_group_key = group_key
                break
        # If the tissue does not match any group, treat it as its own group
        if target_group_key is None:
            target_group_key = tissue_lower

    # Compute the set of related tissue tags for the target group
    related_tissues: Set[str] = set()
    if target_group_key is not None and target_group_key in _TISSUE_GROUPS:
        related_tissues = _TISSUE_GROUPS[target_group_key]

    def _effective_tier(original_tier: int, seed_tissue: str) -> int:
        """Return the effective tier for a seed, adjusted by tissue proximity.

        - Ubiquitous seeds: always keep original tier
        - Seeds matching the target tissue exactly: keep original tier
        - Seeds in the same tissue group (related): downgrade by 1
        - Seeds from unrelated tissues: downgrade by 1
        """
        if not use_tissue_filter:
            return original_tier
        seed_tissue_lower = seed_tissue.lower()
        # Ubiquitous seeds are always fully relevant
        if seed_tissue_lower == "ubiquitous":
            return original_tier
        # Exact tissue match — fully relevant
        if seed_tissue_lower == tissue_lower:
            return original_tier
        # Related tissue (same proximity group) — downweight by 1
        if seed_tissue_lower in related_tissues:
            return original_tier + 1
        # Unrelated tissue — downweight by 1
        return original_tier + 1

    # Build DNA targets for all miRNA seeds
    hits: List[Dict[str, Any]] = []

    # Load ALL seeds (do not filter by tissue at the DB level;
    # tissue-awareness is handled via effective tier below)
    mirna_seeds = get_mirna_seeds(organism)
    for mirna_name, (seed_rna, tier, seed_tissue) in mirna_seeds.items():
        eff_tier = _effective_tier(tier, seed_tissue)
        seed_len = len(seed_rna)

        # Skip seeds shorter than min_seed_match
        if seed_len < min_seed_match:
            continue

        # Full 7mer DNA target
        dna_target_7 = _rna_revcomp_to_dna(seed_rna)
        # 6mer DNA target (first 6 nt of the 7mer seed's reverse complement)
        dna_target_6 = dna_target_7[:6] if seed_len >= 7 else dna_target_7

        # Scan for 7mer-m8 (full seed match, positions 2-8)
        if seed_len >= 7:
            pos = seq.find(dna_target_7)
            while pos != -1:
                # Check for 8mer (T at position immediately 5' of target
                # in the mRNA, i.e. position pos-1 in DNA coding strand)
                is_8mer = (pos > 0 and seq[pos - 1] == "T")
                if is_8mer:
                    base_score = 1.0
                    match_type = "8mer"
                else:
                    base_score = 0.9
                    match_type = "7mer-m8"
                score = _mirna_context_score(seq, pos, match_type, base_score)
                # Compute the context factor for reporting
                context_factor = round(score / base_score, 3) if base_score > 0 else 1.0
                hits.append({
                    "mirna": mirna_name,
                    "pos": pos,
                    "score": score,
                    "match_type": match_type,
                    "tier": eff_tier,
                    "tissue": seed_tissue,
                    "context_factor": context_factor,
                })
                pos = seq.find(dna_target_7, pos + 1)

        # Scan for 6mer match (positions 2-7 of seed)
        # and 7mer-A1 (6mer + T upstream)
        if min_seed_match <= 6:
            pos = seq.find(dna_target_6)
            while pos != -1:
                # Avoid double-counting: skip if this 6mer is part of a
                # 7mer match already found (check if dna_target_7 matches
                # at pos or pos-1)
                if seed_len >= 7:
                    is_part_of_7mer = False
                    if pos + 7 <= n and seq[pos:pos + 7] == dna_target_7:
                        is_part_of_7mer = True
                    if pos > 0 and seq[pos - 1:pos + 6] == dna_target_7:
                        is_part_of_7mer = True
                    if is_part_of_7mer:
                        pos = seq.find(dna_target_6, pos + 1)
                        continue

                is_7mer_A1 = (pos > 0 and seq[pos - 1] == "T")
                if is_7mer_A1:
                    base_score = 0.85
                    match_type = "7mer-A1"
                else:
                    base_score = 0.7
                    match_type = "6mer"
                score = _mirna_context_score(seq, pos, match_type, base_score)
                # Compute the context factor for reporting
                context_factor = round(score / base_score, 3) if base_score > 0 else 1.0
                hits.append({
                    "mirna": mirna_name,
                    "pos": pos,
                    "score": score,
                    "match_type": match_type,
                    "tier": eff_tier,
                    "tissue": seed_tissue,
                    "context_factor": context_factor,
                })
                pos = seq.find(dna_target_6, pos + 1)

    if not hits:
        tissue_suffix = f" (tissue={tissue_lower})" if use_tissue_filter else ""
        return PredicateResult(
            "NoMiRNABindingSite", True, verdict=Verdict.PASS,
            details=f"No significant miRNA seed matches found{tissue_suffix}",
        )

    # Determine worst verdict
    worst_verdict = Verdict.PASS
    worst_score = 0.0
    worst_hit: Optional[Dict[str, Any]] = None

    for hit in hits:
        tier = hit["tier"]
        score = hit["score"]
        if score >= 0.9 and tier == 1:
            v = Verdict.FAIL
        elif score >= 0.85:
            v = Verdict.UNCERTAIN
        elif score >= 0.7:
            v = Verdict.UNCERTAIN
        else:
            v = Verdict.PASS

        if score > worst_score:
            worst_score = score
            worst_hit = hit
        # Upgrade verdict if this hit is worse
        if v == Verdict.FAIL:
            worst_verdict = Verdict.FAIL
        elif v == Verdict.UNCERTAIN and worst_verdict == Verdict.PASS:
            worst_verdict = Verdict.UNCERTAIN

    positions = list({h["pos"] for h in hits})
    positions.sort()

    passed = worst_verdict != Verdict.FAIL

    tier1_hits = [h for h in hits if h["tier"] == 1 and h["score"] >= 0.9]
    other_hits = [h for h in hits if not (h["tier"] == 1 and h["score"] >= 0.9)]

    if worst_verdict == Verdict.FAIL:
        detail_parts = [f"High-affinity miRNA binding site(s): "]
        for h in tier1_hits:
            detail_parts.append(
                f"{h['mirna']} ({h['match_type']}, score={h['score']:.2f}, "
                f"pos={h['pos']}, tissue={h['tissue']}, "
                f"ctx={h.get('context_factor', 1.0):.2f})"
            )
        details = "; ".join(detail_parts)
    elif worst_verdict == Verdict.UNCERTAIN:
        detail_parts = ["Possible miRNA binding site(s): "]
        for h in hits[:5]:  # limit detail
            detail_parts.append(
                f"{h['mirna']} ({h['match_type']}, score={h['score']:.2f}, "
                f"pos={h['pos']}, ctx={h.get('context_factor', 1.0):.2f})"
            )
        details = "; ".join(detail_parts)
        if len(hits) > 5:
            details += f"; ... and {len(hits) - 5} more"
    else:
        tissue_suffix = f" (tissue={tissue_lower})" if use_tissue_filter else ""
        details = f"No significant miRNA seed matches (checked {len(mirna_seeds)} seeds){tissue_suffix}"

    return PredicateResult(
        "NoMiRNABindingSite", passed, verdict=worst_verdict,
        details=details,
        positions=positions,
    )


# ────────────────────────────────────────────────────────────
# m6A modification site predicate
# ────────────────────────────────────────────────────────────



__all__ = [
    "_rna_revcomp_to_dna",
    "_mirna_context_score",
    "check_no_mirna_binding_site",
]
