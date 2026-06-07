"""
BioCompiler Type System — High-Level Evaluate API
==================================================
All evaluate_* functions that return TypeCheckResult objects.
These are the public-facing API for predicate evaluation.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

from .codon_tables import (
    AA_TO_CODONS,
    BLOSUM62,
    CODON_TABLE,
    PredicateResult,
    _BLOSUM62_MISSING_SCORE,
    _CODON_RAMP_LENGTH,
    _CPG_DENSITY_MULTIPLIER,
    _CPG_GC_RICH_THRESHOLD,
    _DG_AU_PAIR_KCAL,
    _DG_GC_PAIR_KCAL,
    _DG_GU_PAIR_KCAL,
    _EUKARYOTE_GT_PER_BP,
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
    _RESTRICTION_SITE_MIN_LENGTH,
    _TM_BORDERLINE_RATIO,
)
from .checks import (
    check_co_translational_folding,
    check_codon_optimality,
    check_conservation_score,
    check_mrna_secondary_structure,
    check_mrna_stability,
    check_no_avoidable_gt,
    check_no_cpg_island,
    check_no_cryptic_promoter,
    check_no_cryptic_splice,
    check_no_gt_dinucleotide,
    check_no_gt_dinucleotide_soft,
    check_no_restriction_site,
    check_no_stop_codons,
    check_no_unexpected_tm_domain,
    check_valid_coding_seq,
    _compute_codon_ramp_score,
    _is_prokaryotic_organism,
    _resolve_species_cai,
    find_cross_codon_cg,
    find_cross_codon_gt,
    find_cross_codon_restriction,
)
from ..types import Verdict, SLOTMode, TypeCheckResult
from ..sliding_gc import evaluate_sliding_gc, check_sliding_gc, SlidingGCResult, WindowViolation

logger = logging.getLogger(__name__)


def evaluate_gc_in_range(seq: str, gc_lo: float = 0.30, gc_hi: float = 0.70) -> TypeCheckResult:
    """Evaluate whether GC content falls within the specified range."""
    if not seq:
        return TypeCheckResult(
            predicate=f"GCInRange({gc_lo}, {gc_hi})",
            verdict=Verdict.UNCERTAIN,
            violation="Empty sequence",
        )
    seq = seq.upper()
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    if gc_lo <= gc <= gc_hi:
        return TypeCheckResult(
            predicate=f"GCInRange({gc_lo}, {gc_hi})",
            verdict=Verdict.PASS,
        )
    direction = "below" if gc < gc_lo else "above"
    return TypeCheckResult(
        predicate=f"GCInRange({gc_lo}, {gc_hi})",
        verdict=Verdict.FAIL,
        violation=f"GC content {gc:.3f} is {direction} range [{gc_lo}, {gc_hi}]",
    )


def evaluate_no_cryptic_splice(
    seq: str,
    boundaries: List[Tuple[int, int]] | None = None,
    cryptic_threshold: float = 3.0,
    uncertain_lo: float = 0.0,
    organism: str = "",
) -> TypeCheckResult:
    """Evaluate whether the sequence contains cryptic splice sites."""
    # Skip for prokaryotic organisms
    if organism and _is_prokaryotic_organism(organism):
        return TypeCheckResult(
            predicate="NoCrypticSplice",
            verdict=Verdict.PASS,
        )

    # For eukaryotes, use a higher default threshold matching the optimizer
    effective_threshold = cryptic_threshold
    if organism and not _is_prokaryotic_organism(organism):
        effective_threshold = max(cryptic_threshold, 8.0)

    from ..maxentscan import score_donor, score_acceptor

    seq = seq.upper()
    worst_score = _MAXENT_INSUFFICIENT_CONTEXT_SCORE
    worst_pos = -1
    worst_verdict = Verdict.PASS

    # Scan donor sites (GT dinucleotides)
    for i in range(len(seq) - 1):
        if seq[i:i+2] == "GT":
            score = score_donor(seq, i)
            if score <= _MAXENT_INSUFFICIENT_CONTEXT_SCORE:
                score = 0.0

            if score >= effective_threshold:
                v = Verdict.FAIL
            elif uncertain_lo > 0 and score >= uncertain_lo:
                v = Verdict.UNCERTAIN
            else:
                v = Verdict.PASS

            if score > worst_score:
                worst_score = score
                worst_pos = i
                worst_verdict = v

    # Scan acceptor sites (AG dinucleotides)
    for i in range(len(seq) - 1):
        if seq[i:i+2] == "AG":
            score = score_acceptor(seq, i)
            if score <= _MAXENT_INSUFFICIENT_CONTEXT_SCORE:
                score = 0.0

            if score >= effective_threshold:
                v = Verdict.FAIL
            elif uncertain_lo > 0 and score >= uncertain_lo:
                v = Verdict.UNCERTAIN
            else:
                v = Verdict.PASS

            if score > worst_score:
                worst_score = score
                worst_pos = i
                worst_verdict = v

    return TypeCheckResult(
        predicate="NoCrypticSplice",
        verdict=worst_verdict,
        violation=(
            f"Worst cryptic splice score {worst_score:.2f} at pos {worst_pos}"
            if worst_verdict != Verdict.PASS else None
        ),
    )


def evaluate_splice_correct(
    seq: str,
    boundaries: List[Tuple[int, int]] | None = None,
    cellular_context: str = "HEK293T",
) -> TypeCheckResult:
    """Evaluate splice correctness — whether known exon boundaries are respected."""
    seq = seq.upper()

    # Guard: no boundaries or empty boundaries -> no splice correction needed
    if boundaries is None or not boundaries or len(boundaries) < 2:
        return TypeCheckResult(
            predicate="SpliceCorrect",
            verdict=Verdict.PASS,
        )

    # Check canonical splice signals at each intron boundary
    n = len(seq)
    for i in range(len(boundaries) - 1):
        try:
            intron_start = boundaries[i][1]
            intron_end = boundaries[i + 1][0]
        except (IndexError, TypeError):
            continue

        if intron_start >= intron_end:
            continue

        # Check donor (GT) at intron start
        if intron_start + 2 <= n and intron_start >= 0:
            donor = seq[intron_start:intron_start + 2]
            if donor != "GT":
                return TypeCheckResult(
                    predicate="SpliceCorrect",
                    verdict=Verdict.FAIL,
                    violation=f"Non-canonical donor {donor} at pos {intron_start}",
                )

        # Check acceptor (AG) at intron end
        if intron_end - 2 >= 0 and intron_end <= n:
            acceptor = seq[intron_end - 2:intron_end]
            if acceptor != "AG":
                return TypeCheckResult(
                    predicate="SpliceCorrect",
                    verdict=Verdict.FAIL,
                    violation=f"Non-canonical acceptor {acceptor} at pos {intron_end - 2}",
                )

    return TypeCheckResult(
        predicate="SpliceCorrect",
        verdict=Verdict.PASS,
    )


def evaluate_codon_adapted(
    seq: str,
    organism: str = "Homo_sapiens",
    threshold: float = 0.5,
) -> TypeCheckResult:
    """Evaluate whether codon adaptation index meets the threshold."""
    from ..translation import compute_cai

    cai = compute_cai(seq, organism)
    if cai >= threshold:
        return TypeCheckResult(
            predicate=f"CodonAdapted({organism}, {threshold})",
            verdict=Verdict.PASS,
        )
    return TypeCheckResult(
        predicate=f"CodonAdapted({organism}, {threshold})",
        verdict=Verdict.FAIL,
        violation=f"CAI {cai:.4f} is below threshold {threshold}",
    )


def evaluate_no_restriction_site(
    seq: str,
    enzymes: List[str] | None = None,
    enzyme_set: List[str] | None = None,
    min_site_length: int = _RESTRICTION_SITE_MIN_LENGTH,
) -> TypeCheckResult:
    """Evaluate whether the sequence contains restriction enzyme recognition sites."""
    from ..restriction_sites import get_recognition_site

    effective_enzymes = enzymes or enzyme_set or []
    seq = seq.upper()
    violations = []

    for enzyme in effective_enzymes:
        site = get_recognition_site(enzyme)
        if site is None:
            continue
        if len(site) < min_site_length:
            continue
        # Check for IUPAC patterns
        has_iupac = any(b not in "ACGT" for b in site.upper())
        if has_iupac:
            for i in range(len(seq) - len(site) + 1):
                window = seq[i:i + len(site)]
                match = True
                for s_base, p_base in zip(window, site.upper()):
                    from ..constants import IUPAC_EXPAND
                    allowed = IUPAC_EXPAND.get(p_base, p_base)
                    if s_base not in allowed:
                        match = False
                        break
                if match:
                    violations.append((i, enzyme))
        else:
            pos = seq.find(site)
            while pos != -1:
                violations.append((pos, enzyme))
                pos = seq.find(site, pos + 1)

    if violations:
        details = ", ".join(f"{e}@{p}" for p, e in violations)
        return TypeCheckResult(
            predicate="NoRestrictionSite",
            verdict=Verdict.FAIL,
            violation=f"Restriction sites found: {details}",
        )
    return TypeCheckResult(
        predicate="NoRestrictionSite",
        verdict=Verdict.PASS,
    )


def evaluate_in_frame(
    seq: str,
    boundaries: List[Tuple[int, int]] | None = None,
    exon_boundaries: List[Tuple[int, int]] | None = None,
) -> TypeCheckResult:
    """Evaluate whether the coding sequence is in-frame (valid ORF)."""
    effective_boundaries = boundaries or exon_boundaries or [(0, len(seq))]
    seq = seq.upper()

    for start, end in effective_boundaries:
        exon_seq = seq[start:end]
        if len(exon_seq) % 3 != 0:
            return TypeCheckResult(
                predicate="InFrame",
                verdict=Verdict.FAIL,
                violation=f"Exon [{start}, {end}) length {len(exon_seq)} not divisible by 3",
            )
        for i in range(0, len(exon_seq), 3):
            codon = exon_seq[i:i+3]
            if codon not in CODON_TABLE:
                return TypeCheckResult(
                    predicate="InFrame",
                    verdict=Verdict.FAIL,
                    violation=f"Invalid codon '{codon}' at position {start + i}",
                )

    return TypeCheckResult(
        predicate="InFrame",
        verdict=Verdict.PASS,
    )


def evaluate_no_instability_motif(seq: str) -> TypeCheckResult:
    """Evaluate whether the sequence contains mRNA instability motifs."""
    seq = seq.upper()
    positions = []

    # Check for ATTTA (AUUUA in mRNA)
    for i in range(len(seq) - 4):
        if seq[i:i+5] == "ATTTA":
            positions.append(i)

    # Check for U-rich regions (6+ consecutive T's in DNA)
    i = 0
    while i < len(seq):
        if seq[i] == "T":
            run_start = i
            while i < len(seq) and seq[i] == "T":
                i += 1
            run_len = i - run_start
            if run_len >= _INSTABILITY_T_RUN_MIN:
                positions.append(run_start)
        else:
            i += 1

    if positions:
        return TypeCheckResult(
            predicate="NoInstabilityMotif",
            verdict=Verdict.FAIL,
            violation=f"Instability motifs at positions {positions}",
        )
    return TypeCheckResult(
        predicate="NoInstabilityMotif",
        verdict=Verdict.PASS,
    )


def evaluate_no_unexpected_tm_domain(
    seq: str,
    is_cytosolic: bool = True,
    window_size: int = 19,
    threshold: float = 0.68,
    organism: str = "",
) -> TypeCheckResult:
    """Evaluate whether a cytosolic protein has gained unexpected TM domains."""
    result = check_no_unexpected_tm_domain(
        seq, is_cytosolic=is_cytosolic, window_size=window_size,
        threshold=threshold, organism=organism,
    )

    if not is_cytosolic:
        return TypeCheckResult(
            predicate="NoUnexpectedTMDomain",
            verdict=Verdict.PASS,
        )

    if result.verdict == Verdict.FAIL:
        return TypeCheckResult(
            predicate=f"NoUnexpectedTMDomain({is_cytosolic}, {threshold})",
            verdict=Verdict.FAIL,
            violation=result.details,
        )
    elif result.verdict == Verdict.UNCERTAIN:
        return TypeCheckResult(
            predicate=f"NoUnexpectedTMDomain({is_cytosolic}, {threshold})",
            verdict=Verdict.UNCERTAIN,
            violation=result.details,
        )
    else:
        return TypeCheckResult(
            predicate=f"NoUnexpectedTMDomain({is_cytosolic}, {threshold})",
            verdict=Verdict.LIKELY_PASS,
        )


def evaluate_mrna_secondary_structure(
    seq: str,
    window_start: int = 0,
    window_end: int = 50,
    dg_threshold: float = -15.0,
    organism: str = "",
) -> TypeCheckResult:
    """Evaluate mRNA secondary structure stability around the RBS/start codon."""
    # Apply organism-specific dG threshold
    effective_threshold = dg_threshold
    if organism:
        if _is_prokaryotic_organism(organism):
            effective_threshold = _MRNA_DG_PROKARYOTE_FAIL
        else:
            effective_threshold = _MRNA_DG_EUKARYOTE_FAIL

    seq = seq.upper()
    effective_end = min(window_end, len(seq))
    window_seq = seq[window_start:effective_end]

    if len(window_seq) < 4:
        return TypeCheckResult(
            predicate="mRNASecondaryStructure",
            verdict=Verdict.LIKELY_PASS,
        )

    # Convert DNA to RNA for pairing analysis (T -> U)
    rna = window_seq.replace("T", "U")

    # Count potential base pairs using simplified hairpin model
    gc_pairs = 0
    au_pairs = 0
    gu_pairs = 0

    half = len(rna) // 2
    first_half = rna[:half]
    second_half = rna[half:2 * half]

    for i in range(min(len(first_half), len(second_half))):
        j = len(second_half) - 1 - i
        if j < 0:
            break
        base_5 = first_half[i]
        base_3 = second_half[j]

        if (base_5 == "G" and base_3 == "C") or (base_5 == "C" and base_3 == "G"):
            gc_pairs += 1
        elif (base_5 == "A" and base_3 == "U") or (base_5 == "U" and base_3 == "A"):
            au_pairs += 1
        elif (base_5 == "G" and base_3 == "U") or (base_5 == "U" and base_3 == "G"):
            gu_pairs += 1

    # Simplified nearest-neighbor dG estimate
    dg = _DG_GC_PAIR_KCAL * gc_pairs + _DG_AU_PAIR_KCAL * au_pairs + _DG_GU_PAIR_KCAL * gu_pairs

    if dg <= effective_threshold:
        verdict = Verdict.FAIL
    elif dg <= effective_threshold * _MRNA_MODERATE_DG_RATIO:
        verdict = Verdict.UNCERTAIN
    else:
        verdict = Verdict.LIKELY_PASS

    violation = None
    if verdict == Verdict.FAIL:
        violation = (
            f"Strong mRNA secondary structure: dG={dg:.1f} kcal/mol "
            f"<= {effective_threshold}"
        )
    elif verdict == Verdict.UNCERTAIN:
        # Refine using GC content of the window
        window_gc = (window_seq.count('G') + window_seq.count('C')) / len(window_seq) if window_seq else 0.5
        if window_gc < 0.5:
            verdict = Verdict.LIKELY_PASS
            violation = (
                f"Moderate mRNA secondary structure: dG={dg:.1f} kcal/mol, "
                f"but AT-rich window (GC={window_gc:.0%}) weakens structure"
            )
        else:
            verdict = Verdict.LIKELY_FAIL
            violation = (
                f"Moderate mRNA secondary structure: dG={dg:.1f} kcal/mol "
                f"<= {effective_threshold * _MRNA_MODERATE_DG_RATIO:.1f}"
            )

    return TypeCheckResult(
        predicate=f"mRNASecondaryStructure({window_start}, {window_end}, {effective_threshold})",
        verdict=verdict,
        violation=violation,
    )


def evaluate_no_cryptic_promoter(
    seq: str,
    organism: str = "E_coli",
    threshold: float = 0.7,
) -> TypeCheckResult:
    """Evaluate whether the sequence contains cryptic promoter sites."""
    result = check_no_cryptic_promoter(seq, organism, threshold)

    if result.verdict == Verdict.FAIL:
        return TypeCheckResult(
            predicate=f"NoCrypticPromoter({organism}, {threshold})",
            verdict=Verdict.FAIL,
            violation=result.details,
        )
    elif result.verdict == Verdict.UNCERTAIN:
        return TypeCheckResult(
            predicate=f"NoCrypticPromoter({organism}, {threshold})",
            verdict=Verdict.UNCERTAIN,
            violation=result.details,
        )
    else:
        return TypeCheckResult(
            predicate=f"NoCrypticPromoter({organism}, {threshold})",
            verdict=Verdict.LIKELY_PASS,
        )


def evaluate_no_cpg_island(
    seq: str,
    window: int = 200,
    threshold: float = 0.6,
    organism: str = "",
) -> TypeCheckResult:
    """Evaluate whether the sequence contains CpG islands."""
    result = check_no_cpg_island(seq, window=window, threshold=threshold, organism=organism)

    if result.verdict == Verdict.FAIL:
        return TypeCheckResult(
            predicate="NoCpGIsland",
            verdict=Verdict.FAIL,
            violation=result.details,
        )
    return TypeCheckResult(
        predicate="NoCpGIsland",
        verdict=Verdict.PASS,
    )


def analyze_codon_at_position(
    seq: str,
    position: int,
    organism: str = "Homo_sapiens",
) -> Dict[str, Any]:
    """Analyze the codon at a given position for optimality and alternatives."""
    codon_start = (position // 3) * 3
    if codon_start + 3 > len(seq):
        return {"codon": "N/A", "amino_acid": "N/A", "cai": 0.0, "alternatives": [], "position": codon_start}

    seq = seq.upper()
    codon = seq[codon_start:codon_start + 3]
    aa = CODON_TABLE.get(codon, "?")
    species_cai = _resolve_species_cai(organism)
    cai = species_cai.get(codon, 0.0)

    alternatives = []
    for alt in AA_TO_CODONS.get(aa, []):
        if alt != codon:
            alternatives.append({
                "codon": alt,
                "cai": species_cai.get(alt, 0.0),
            })

    return {
        "codon": codon,
        "amino_acid": aa,
        "cai": cai,
        "alternatives": sorted(alternatives, key=lambda x: x["cai"], reverse=True),
        "position": codon_start,
    }


def evaluate_co_translational_folding(
    seq: str,
    organism: str = "Homo_sapiens",
    domain_boundaries: List[int] | None = None,
    min_pause_cai: float = 0.3,
    structure_confidence: float | None = None,
    plddt_score: float | None = None,
) -> TypeCheckResult:
    """Evaluate co-translational folding preservation after codon optimization."""
    seq = seq.upper()

    # Resolve organism name to flat codon->CAI-weight dict
    species_cai = _resolve_species_cai(organism)

    # Run the low-level predicate check
    result = check_co_translational_folding(
        seq, species_cai, domain_boundaries, min_pause_cai,
        structure_confidence=structure_confidence,
        plddt_score=plddt_score,
    )

    # Map the PredicateResult verdict to a TypeCheckResult verdict using
    # the five-valued logic, incorporating additional context about the
    # severity of the findings.
    ramp_info = _compute_codon_ramp_score(seq, species_cai)
    ramp_score = ramp_info["ramp_score"]
    speed_disruptions = ramp_info["speed_disruptions"]
    pause_sites = ramp_info["pause_sites"]

    ramp_length = min(_CODON_RAMP_LENGTH, len(seq) // 3)
    ramp_all_fast = ramp_length >= _MIN_RAMP_FOR_WARNING and all(
        species_cai.get(seq[i * 3:(i + 1) * 3], 0.0) > _FAST_CODON_CAI_THRESHOLD
        for i in range(ramp_length)
    )

    # Count domain boundary disruptions
    domain_disrupted = 0
    if domain_boundaries:
        num_codons = len(seq) // 3
        for bp in domain_boundaries:
            if 0 <= bp < num_codons:
                codon = seq[bp * 3:(bp + 1) * 3]
                if species_cai.get(codon, 0.0) > _FAST_CODON_CAI_THRESHOLD:
                    domain_disrupted += 1

    # Determine the TypeCheckResult verdict
    if domain_boundaries and domain_disrupted > 0 and ramp_all_fast:
        verdict = Verdict.FAIL
        violation = (
            f"Ramp destroyed (avg CAI={ramp_score:.3f}) and {domain_disrupted} "
            f"domain boundary(ies) disrupted by fast codons"
        )
    elif domain_boundaries and domain_disrupted >= 2:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"{domain_disrupted} domain boundary(ies) have fast codons "
            f"(CAI > 0.7) where pause sites are needed"
        )
    elif domain_boundaries and domain_disrupted == 1:
        verdict = Verdict.UNCERTAIN
        violation = (
            f"1 domain boundary has a fast codon (CAI > 0.7) where "
            f"a pause site may be needed"
        )
    elif ramp_all_fast and ramp_length >= _MIN_RAMP_FOR_WARNING:
        verdict = Verdict.UNCERTAIN
        violation = (
            f"Codon ramp (first {ramp_length} codons) is entirely fast "
            f"(avg CAI={ramp_score:.3f}) — ribosome jam risk"
        )
    elif speed_disruptions:
        if len(speed_disruptions) <= 2:
            verdict = Verdict.LIKELY_PASS
            violation = (
                f"{len(speed_disruptions)} potential pause site(s) replaced by "
                f"fast codons — minor concern"
            )
        else:
            verdict = Verdict.UNCERTAIN
            violation = (
                f"{len(speed_disruptions)} potential pause site(s) replaced by "
                f"fast codons — may affect co-translational folding"
            )
    elif pause_sites:
        verdict = Verdict.PASS
        violation = None
    else:
        verdict = Verdict.LIKELY_PASS
        violation = None

    return TypeCheckResult(
        predicate=f"CoTranslationalFolding({organism}, {min_pause_cai})",
        verdict=verdict,
        violation=violation,
    )


def evaluate_no_stop_codons(seq: str) -> TypeCheckResult:
    """Evaluate whether the DNA sequence contains internal stop codons."""
    result = check_no_stop_codons(seq)
    if result.passed:
        return TypeCheckResult(
            predicate="NoStopCodons",
            verdict=Verdict.PASS,
        )
    return TypeCheckResult(
        predicate="NoStopCodons",
        verdict=Verdict.FAIL,
        violation=result.details,
    )


def evaluate_no_gt_dinucleotide(
    seq: str,
    organism: str = "",
    max_gt_count: int | None = None,
) -> TypeCheckResult:
    """Evaluate whether the sequence contains avoidable GT dinucleotides."""
    result = check_no_gt_dinucleotide_soft(
        seq, organism=organism, max_gt_count=max_gt_count,
    )

    if result.verdict == Verdict.PASS:
        return TypeCheckResult(
            predicate="NoGTDinucleotide",
            verdict=Verdict.PASS,
        )
    elif result.verdict == Verdict.LIKELY_FAIL:
        return TypeCheckResult(
            predicate="NoGTDinucleotide",
            verdict=Verdict.LIKELY_FAIL,
            violation=result.details,
        )
    else:
        return TypeCheckResult(
            predicate="NoGTDinucleotide",
            verdict=Verdict.FAIL,
            violation=result.details,
        )


def evaluate_valid_coding_seq(seq: str) -> TypeCheckResult:
    """Evaluate whether the DNA sequence is a valid coding sequence."""
    result = check_valid_coding_seq(seq)
    if result.passed:
        return TypeCheckResult(
            predicate="ValidCodingSeq",
            verdict=Verdict.PASS,
        )
    return TypeCheckResult(
        predicate="ValidCodingSeq",
        verdict=Verdict.FAIL,
        violation=result.details,
    )


def evaluate_conservation_score(
    seq: str,
    protein: str = "",
    min_score: int = 0,
) -> TypeCheckResult:
    """Evaluate whether the BLOSUM62 conservation score meets the threshold."""
    if not protein:
        return TypeCheckResult(
            predicate="ConservationScore",
            verdict=Verdict.PASS,
        )
    result = check_conservation_score(seq, protein, min_score=min_score)
    if result.passed:
        return TypeCheckResult(
            predicate="ConservationScore",
            verdict=Verdict.PASS,
        )
    return TypeCheckResult(
        predicate="ConservationScore",
        verdict=Verdict.FAIL,
        violation=result.details,
    )


def evaluate_codon_optimality(
    seq: str,
    organism: str = "Homo_sapiens",
    threshold: float = 0.5,
) -> TypeCheckResult:
    """Evaluate whether the codon adaptation index (CAI) meets the threshold."""
    result = check_codon_optimality(seq, organism, min_cai=threshold)
    if result.passed:
        return TypeCheckResult(
            predicate="CodonOptimality",
            verdict=Verdict.PASS,
        )
    return TypeCheckResult(
        predicate="CodonOptimality",
        verdict=Verdict.FAIL,
        violation=result.details,
    )


def evaluate_mrna_stability(
    seq: str,
    organism: str = "Homo_sapiens",
    threshold: float | None = None,
) -> TypeCheckResult:
    """Evaluate mRNA stability for a coding sequence."""
    result = check_mrna_stability(seq, organism, threshold)

    if result.verdict == Verdict.PASS:
        return TypeCheckResult(
            predicate=f"MRNAStability({organism})",
            verdict=Verdict.LIKELY_PASS if result.details.startswith("mRNA stability score") else Verdict.PASS,
        )
    elif result.verdict == Verdict.UNCERTAIN:
        return TypeCheckResult(
            predicate=f"MRNAStability({organism})",
            verdict=Verdict.UNCERTAIN,
            violation=result.details,
        )
    else:
        return TypeCheckResult(
            predicate=f"MRNAStability({organism})",
            verdict=Verdict.FAIL,
            violation=result.details,
        )


def evaluate_all_predicates(
    seq: str,
    boundaries: List[Tuple[int, int]] | None = None,
    known_exon_boundaries: List[Tuple[int, int]] | None = None,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    enzymes: List[str] | None = None,
    cryptic_threshold: float = 3.0,
    uncertain_lo: float = 0.0,
    cai_threshold: float = 0.5,
    promoter_threshold: float = 0.7,
    is_cytosolic: bool = True,
    tm_threshold: float = 0.68,
    mrna_window: int = 50,
    mrna_dg_threshold: float = -15.0,
    folding_threshold: float = 0.3,
    domain_boundaries: List[int] | None = None,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> List[TypeCheckResult]:
    """Evaluate all 12 type predicates against a sequence."""
    # Backward compatibility: accept known_exon_boundaries as alias for boundaries
    if boundaries is None and known_exon_boundaries is not None:
        boundaries = known_exon_boundaries

    # Track whether the caller supplied real exon boundaries.
    has_real_boundaries = boundaries is not None and len(boundaries) >= 2

    if boundaries is None:
        boundaries = [(0, len(seq))]

    if enzymes is None:
        enzymes = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
    # Map organism name to promoter organism key
    promoter_organism = "E_coli" if organism in ("E_coli", "ecoli", "Escherichia_coli") else "eukaryote"

    # Import slot_verification lazily to avoid circular imports
    from ..slot_verification import is_slot_predicate, verify_slot_predicate

    # Core (non-SLOT) predicates: always evaluate normally
    if has_real_boundaries:
        splice_result = evaluate_splice_correct(seq, boundaries)
    else:
        splice_result = TypeCheckResult(
            predicate="SpliceCorrect",
            verdict=Verdict.PASS,
        )

    results: List[TypeCheckResult] = [
        splice_result,
        evaluate_gc_in_range(seq, gc_lo, gc_hi),
        evaluate_no_restriction_site(seq, enzymes),
        evaluate_in_frame(seq, boundaries),
        evaluate_no_instability_motif(seq),
        evaluate_no_cpg_island(seq, organism=organism),
    ]

    # SLOT predicates: behavior depends on slot_mode
    slot_predicates = [
        ("NoCrypticSplice", lambda: evaluate_no_cryptic_splice(seq, boundaries, cryptic_threshold, uncertain_lo, organism=organism)),
        ("CodonAdapted", lambda: evaluate_codon_adapted(seq, organism, cai_threshold)),
        ("NoCrypticPromoter", lambda: evaluate_no_cryptic_promoter(seq, promoter_organism, promoter_threshold)),
        ("NoUnexpectedTMDomain", lambda: evaluate_no_unexpected_tm_domain(seq, is_cytosolic, 19, tm_threshold)),
        ("mRNASecondaryStructure", lambda: evaluate_mrna_secondary_structure(seq, 0, mrna_window, mrna_dg_threshold)),
        ("CoTranslationalFolding", lambda: evaluate_co_translational_folding(seq, organism, domain_boundaries, folding_threshold)),
    ]

    for pred_name, eval_fn in slot_predicates:
        if slot_mode == SLOTMode.CONSERVATIVE and is_slot_predicate(pred_name):
            result = eval_fn()
            result.verdict = Verdict.UNCERTAIN
            result.knowledge_gap = f"SLOT predicate: {pred_name} returns UNCERTAIN in CONSERVATIVE mode"
            results.append(result)
        elif slot_mode in (SLOTMode.VERIFIED, SLOTMode.PERMISSIVE) and is_slot_predicate(pred_name):
            verdict, evidence = verify_slot_predicate(
                pred_name,
                slot_mode=slot_mode,
                seq=seq,
                low_thresh=uncertain_lo,
                high_thresh=cryptic_threshold,
                organism=promoter_organism,
                threshold=promoter_threshold,
                is_cytosolic=is_cytosolic,
                window_end=mrna_window,
                dg_threshold=mrna_dg_threshold,
                domain_boundaries=domain_boundaries,
                min_pause_cai=folding_threshold,
            )
            result = eval_fn()
            result.verdict = verdict
            result.knowledge_gap = f"SLOT predicate: {pred_name} evaluated in {slot_mode.value} mode"
            if evidence.verified:
                result.derivation = [{"evidence": evidence.to_dict()}]
            results.append(result)
        else:
            results.append(eval_fn())

    # Reorder to match the expected predicate order
    name_to_result = {r.predicate: r for r in results}
    ordered = []
    for canonical in [
        "NoCrypticSplice", "SpliceCorrect", "GCInRange", "CodonAdapted",
        "NoRestrictionSite", "InFrame", "NoInstabilityMotif", "NoCpGIsland",
        "NoCrypticPromoter", "NoUnexpectedTMDomain", "mRNASecondaryStructure",
        "CoTranslationalFolding",
    ]:
        if canonical in name_to_result:
            ordered.append(name_to_result[canonical])
        else:
            for name, result in name_to_result.items():
                if name.startswith(canonical):
                    ordered.append(result)
                    break

    if len(ordered) != 12:
        ordered = results

    return ordered
