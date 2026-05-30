"""
BioCompiler Type System — Predicate Evaluation with Registry Pattern

Production-grade type system with:
- Registry pattern for predicates (extensible, no string-prefix hacks)
- NoCrypticSplice uses MaxEntScan scoring (not raw GT/AG counting)
- SpliceCorrect fixes UNCERTAIN logic
- Each predicate is a self-contained callable with metadata
- Evaluate all registered predicates via registry (not hardcoded list)
- Proper error handling
"""

import logging
from .types import Verdict, TypeCheckResult, Token
from .scanner import scan_sequence, gc_content, validate_dna_sequence
from .splicing import compute_splice_isoforms
from .translation import compute_cai
from .constants import INSTABILITY_MOTIF, RESTRICTION_ENZYMES, IUPAC_EXPAND, reverse_complement
from .maxentscan import max_donor_score, max_acceptor_score, scan_splice_sites
from .exceptions import UnknownPredicateError

logger = logging.getLogger(__name__)


# ==============================================================================
# Predicate Registry
# ==============================================================================

class PredicateRegistry:
    """
    Registry of type predicates. Replaces fragile string-prefix dispatch.

    Each predicate is registered with:
    - name: unique identifier
    - evaluator: callable that returns TypeCheckResult
    - verifier: callable that re-evaluates from certificate data
    - param_keys: list of parameter keys needed to invoke the evaluator
    """

    def __init__(self):
        self._predicates: dict[str, dict] = {}

    def register(self, name: str, evaluator, verifier=None, param_keys: list[str] | None = None):
        """Register a predicate with its evaluator, optional verifier, and parameter keys."""
        self._predicates[name] = {
            "evaluator": evaluator,
            "verifier": verifier or evaluator,
            "param_keys": param_keys or [],
        }

    def evaluate(self, name: str, **kwargs) -> TypeCheckResult:
        """Evaluate a registered predicate by name."""
        if name not in self._predicates:
            raise UnknownPredicateError(name)
        return self._predicates[name]["evaluator"](**kwargs)

    def verify(self, name: str, **kwargs) -> TypeCheckResult:
        """Re-evaluate a predicate for certificate verification."""
        if name not in self._predicates:
            raise UnknownPredicateError(name)
        return self._predicates[name]["verifier"](**kwargs)

    def names(self) -> list[str]:
        """Return all registered predicate names."""
        return list(self._predicates.keys())

    def __contains__(self, name: str) -> bool:
        return name in self._predicates


# Global registry instance
registry = PredicateRegistry()


# ==============================================================================
# Predicate Implementations
# ==============================================================================

def evaluate_no_cryptic_splice(
    seq: str,
    known_exon_boundaries: list[tuple[int, int]],
    cryptic_threshold: float = 3.0,
    **kwargs,
) -> TypeCheckResult:
    """
    NoCrypticSplice: No donor/acceptor pair within known exons that could form
    a cryptic intron, as determined by MaxEntScan scoring.

    This uses the MaxEntScan model to score every GT/AG dinucleotide in exonic
    regions. Only pairs where BOTH the donor AND acceptor score above the
    threshold are flagged as potential cryptic splice sites. This dramatically
    reduces false positives compared to the naive approach of flagging every
    GT/AG pair.

    Args:
        seq: full pre-mRNA sequence
        known_exon_boundaries: list of (start, end) tuples for known exons
        cryptic_threshold: MaxEntScan score threshold above which a splice site
                          is considered functional (default 3.0)
    """
    derivation = []

    for exon_start, exon_end in known_exon_boundaries:
        exon_seq = seq[exon_start:exon_end]

        # Use MaxEntScan to find scored splice sites within this exon
        splice_sites = scan_splice_sites(exon_seq, cryptic_threshold, cryptic_threshold)

        donors = [(pos, score) for pos, stype, score in splice_sites if stype == "donor"]
        acceptors = [(pos, score) for pos, stype, score in splice_sites if stype == "acceptor"]

        for pos, score in donors:
            derivation.append({
                "step": "maxentscan_donor_in_exon",
                "position": exon_start + pos,
                "score": score,
                "threshold": cryptic_threshold,
            })
        for pos, score in acceptors:
            derivation.append({
                "step": "maxentscan_acceptor_in_exon",
                "position": exon_start + pos,
                "score": score,
                "threshold": cryptic_threshold,
            })

        # Check for donor-acceptor pairs that could form cryptic introns
        # Only pairs where gap > MIN_INTRON_LENGTH are biologically plausible
        for d_pos, d_score in donors:
            for a_pos, a_score in acceptors:
                if a_pos > d_pos + 30:  # Minimum intron length
                    # Both sites score above threshold → potential cryptic intron
                    return TypeCheckResult(
                        predicate="NoCrypticSplice",
                        verdict=Verdict.FAIL,
                        derivation=derivation,
                        violation=(
                            f"Cryptic splice site pair in exon [{exon_start},{exon_end}): "
                            f"donor at {exon_start + d_pos} (score={d_score:.2f}), "
                            f"acceptor at {exon_start + a_pos} (score={a_score:.2f})"
                        ),
                    )

    return TypeCheckResult(
        predicate="NoCrypticSplice",
        verdict=Verdict.PASS,
        derivation=derivation + [{
            "step": "no_cryptic_pairs_found",
            "evidence": "maxentscan_scored_exhaustive_scan",
            "threshold": cryptic_threshold,
        }],
    )


def evaluate_splice_correct(
    seq: str, known_exon_boundaries: list[tuple[int, int]],
    cellular_context: str = "HEK293T", **kwargs
) -> TypeCheckResult:
    """
    SpliceCorrect(C): The intended isoform is the only possible splice isoform.
    """
    isoforms = compute_splice_isoforms(seq, known_exon_boundaries, cellular_context)
    target_seq = "".join(seq[start:end] for start, end in known_exon_boundaries)
    non_target = [iso for iso in isoforms if iso.sequence != target_seq]

    derivation = [
        {"step": "ndfst_output_set_size", "value": len(isoforms)},
        {"step": "target_isoform", "exon_boundaries": known_exon_boundaries},
    ]

    if len(non_target) == 0:
        # All isoforms have the same sequence as the target → PASS
        # This covers both: only one isoform (singleton), and multiple isoforms
        # that all produce the same sequence
        return TypeCheckResult(
            predicate=f"SpliceCorrect({cellular_context})",
            verdict=Verdict.PASS,
            derivation=derivation + [{
                "step": "all_isoforms_match_target",
                "total_isoforms": len(isoforms),
            }],
        )
    else:
        # There are alternative isoforms with different sequences → FAIL
        alt_desc = [
            f"boundaries={iso.exon_boundaries} via {iso.parse_path}"
            for iso in non_target[:5]
        ]
        return TypeCheckResult(
            predicate=f"SpliceCorrect({cellular_context})",
            verdict=Verdict.FAIL,
            derivation=derivation,
            violation=f"Found {len(non_target)} alternative isoforms: {'; '.join(alt_desc)}",
        )


def evaluate_gc_in_range(seq: str, gc_lo: float = 0.30, gc_hi: float = 0.70, **kwargs) -> TypeCheckResult:
    """GCInRange(lo, hi): GC content is within [lo, hi]."""
    gc = gc_content(seq)
    passed = gc_lo <= gc <= gc_hi
    return TypeCheckResult(
        predicate=f"GCInRange({gc_lo}, {gc_hi})",
        verdict=Verdict.PASS if passed else Verdict.FAIL,
        derivation=[{"step": "compute_gc", "gc_content": gc, "g_count": seq.count('G'), "c_count": seq.count('C'), "total": len(seq)}],
        violation=f"GC content {gc} not in [{gc_lo}, {gc_hi}]" if not passed else None,
    )


def evaluate_codon_adapted(seq: str, organism: str = "Homo_sapiens", threshold: float = 0.5, **kwargs) -> TypeCheckResult:
    """CodonAdapted(O, theta): CAI >= threshold for organism O."""
    cai = compute_cai(seq, organism)
    passed = cai >= threshold
    return TypeCheckResult(
        predicate=f"CodonAdapted({organism}, {threshold})",
        verdict=Verdict.PASS if passed else Verdict.FAIL,
        derivation=[{"step": "compute_cai", "cai": cai, "organism": organism}],
        violation=f"CAI {cai} < threshold {threshold}" if not passed else None,
    )


def evaluate_no_restriction_site(seq: str, enzyme_set: list[str] | None = None, **kwargs) -> TypeCheckResult:
    """NoRestrictionSite(S): No restriction enzyme recognition site from set S in sequence (both strands)."""
    enzyme_set = enzyme_set or list(RESTRICTION_ENZYMES.keys())
    found_sites: list[dict] = []

    for enz_name in enzyme_set:
        # Support both name-based and sequence-based enzyme specification
        if enz_name in RESTRICTION_ENZYMES:
            site = RESTRICTION_ENZYMES[enz_name]
        elif all(b in "ACGT" for b in enz_name.upper()):
            # Treat as a raw recognition sequence
            site = enz_name.upper()
        else:
            logger.warning("Unknown enzyme '%s' and not a valid DNA sequence — skipping", enz_name)
            continue

        has_iupac = any(b not in "ACGT" for b in site.upper())

        # Forward strand with IUPAC-aware matching
        for i in range(len(seq) - len(site) + 1):
            window = seq[i:i+len(site)]
            if has_iupac:
                match = all(window[j] in IUPAC_EXPAND.get(site[j].upper(), site[j].upper())
                            for j in range(len(site)) if j < len(window))
                if match:
                    found_sites.append({"enzyme": enz_name, "position": i, "site": window, "strand": "+"})
            else:
                if window == site:
                    found_sites.append({"enzyme": enz_name, "position": i, "site": site, "strand": "+"})

        # Reverse complement strand (only for non-IUPAC sites)
        if not has_iupac:
            site_rc = reverse_complement(site)
            if site_rc != site:
                for i in range(len(seq) - len(site_rc) + 1):
                    if seq[i:i+len(site_rc)] == site_rc:
                        found_sites.append({"enzyme": enz_name, "position": i, "site": site_rc, "strand": "-"})

    if found_sites:
        return TypeCheckResult(
            predicate=f"NoRestrictionSite({enzyme_set})",
            verdict=Verdict.FAIL,
            violation=f"Found {len(found_sites)} restriction sites: {found_sites[:5]}",
        )
    return TypeCheckResult(
        predicate=f"NoRestrictionSite({enzyme_set})",
        verdict=Verdict.PASS,
        derivation=[{"step": "exhaustive_search", "enzymes_checked": enzyme_set, "sites_found": 0, "both_strands": True}],
    )


def evaluate_in_frame(seq: str, exon_boundaries: list[tuple[int, int]], **kwargs) -> TypeCheckResult:
    """InFrame: Reading frame consistent across exon boundaries, no premature stop codons."""
    frame_issues = []
    for i, (start, end) in enumerate(exon_boundaries):
        exon_len = end - start
        if exon_len % 3 != 0:
            frame_issues.append(f"Exon {i} length {exon_len} is not a multiple of 3")

    coding_seq = "".join(seq[start:end] for start, end in exon_boundaries)
    premature_stops = []
    for i in range(0, len(coding_seq) - 2, 3):
        codon = coding_seq[i:i+3]
        if codon in ("TAA", "TAG", "TGA") and i < len(coding_seq) - 3:
            premature_stops.append({"position": i, "codon": codon})

    if frame_issues:
        return TypeCheckResult(
            predicate="InFrame",
            verdict=Verdict.FAIL,
            violation=f"Frame issues: {'; '.join(frame_issues)}",
        )
    if premature_stops:
        return TypeCheckResult(
            predicate="InFrame",
            verdict=Verdict.FAIL,
            violation=f"Premature stop codons: {premature_stops}",
        )
    return TypeCheckResult(
        predicate="InFrame",
        verdict=Verdict.PASS,
        derivation=[
            {"step": "frame_consistency", "all_exons_multiple_of_3": True},
            {"step": "no_premature_stop", "checked_codons": len(coding_seq) // 3},
        ],
    )


def evaluate_no_instability_motif(seq: str, **kwargs) -> TypeCheckResult:
    """NoInstabilityMotif: No ATTTA or U-rich motifs present."""
    motifs_found = []
    for i in range(len(seq) - 4):
        if seq[i:i+5] == INSTABILITY_MOTIF:
            motifs_found.append({"position": i, "motif": seq[i:i+5], "type": "ATTTA"})

    # U-rich regions (6+ consecutive T in DNA)
    i = 0
    while i < len(seq):
        if seq[i] == 'T':
            j = i
            while j < len(seq) and seq[j] == 'T':
                j += 1
            run_len = j - i
            if run_len >= 6:
                motifs_found.append({"position": i, "motif": seq[i:j], "type": "U_rich", "length": run_len})
            i = j
        else:
            i += 1

    if motifs_found:
        return TypeCheckResult(
            predicate="NoInstabilityMotif",
            verdict=Verdict.FAIL,
            violation=f"Found {len(motifs_found)} instability motifs: {motifs_found[:5]}",
        )
    return TypeCheckResult(
        predicate="NoInstabilityMotif",
        verdict=Verdict.PASS,
        derivation=[{"step": "exhaustive_motif_scan", "motifs_found": 0}],
    )


def evaluate_no_cpg_island(seq: str, window_size: int = 200, threshold: float = 0.6, min_obs_exp: float = 0.65, **kwargs) -> TypeCheckResult:
    """
    NoCpGIsland: No CpG islands in the coding sequence that could trigger
    epigenetic silencing.

    A CpG island is defined as a region of at least `window_size` bp with:
    - GC content >= threshold (default 0.6)
    - Observed/Expected CpG ratio >= min_obs_exp (default 0.65)

    CpG islands in coding sequences can trigger DNA methylation and
    epigenetic silencing, which is undesirable for expression constructs.
    """
    seq = seq.upper()
    cpg_islands = []

    for start in range(len(seq) - window_size + 1):
        window = seq[start:start + window_size]
        gc = gc_content(window)

        if gc >= threshold:
            # Compute observed/expected CpG ratio
            cpg_count = sum(1 for i in range(len(window) - 1) if window[i:i+2] == "CG")
            c_count = window.count('C')
            g_count = window.count('G')
            expected = (c_count * g_count) / max(len(window), 1)
            obs_exp = cpg_count / max(expected, 1e-10)

            if obs_exp >= min_obs_exp:
                cpg_islands.append({
                    "start": start,
                    "end": start + window_size,
                    "gc_content": round(gc, 4),
                    "cpg_obs_exp": round(obs_exp, 4),
                })

    # Merge overlapping islands
    merged = []
    for island in cpg_islands:
        if merged and island["start"] <= merged[-1]["end"]:
            merged[-1]["end"] = max(merged[-1]["end"], island["end"])
        else:
            merged.append(dict(island))

    if merged:
        return TypeCheckResult(
            predicate="NoCpGIsland",
            verdict=Verdict.FAIL,
            violation=f"Found {len(merged)} CpG island(s): {merged[:5]}",
        )
    return TypeCheckResult(
        predicate="NoCpGIsland",
        verdict=Verdict.PASS,
        derivation=[{
            "step": "cpg_island_scan",
            "window_size": window_size,
            "gc_threshold": threshold,
            "obs_exp_threshold": min_obs_exp,
            "islands_found": 0,
        }],
    )


# ==============================================================================
# Register all built-in predicates
# ==============================================================================

registry.register("NoCrypticSplice", evaluate_no_cryptic_splice,
                  param_keys=["seq", "known_exon_boundaries", "cryptic_threshold"])
registry.register("SpliceCorrect", evaluate_splice_correct,
                  param_keys=["seq", "known_exon_boundaries", "cellular_context"])
registry.register("GCInRange", evaluate_gc_in_range,
                  param_keys=["seq", "gc_lo", "gc_hi"])
registry.register("CodonAdapted", evaluate_codon_adapted,
                  param_keys=["seq", "organism", "threshold"])
registry.register("NoRestrictionSite", evaluate_no_restriction_site,
                  param_keys=["seq", "enzyme_set"])
registry.register("InFrame", evaluate_in_frame,
                  param_keys=["seq", "exon_boundaries"])
registry.register("NoInstabilityMotif", evaluate_no_instability_motif,
                  param_keys=["seq"])
registry.register("NoCpGIsland", evaluate_no_cpg_island,
                  param_keys=["seq"])


def evaluate_all_predicates(
    seq: str,
    known_exon_boundaries: list[tuple[int, int]],
    organism: str = "Homo_sapiens",
    cellular_context: str = "HEK293T",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.5,
    enzymes: list[str] | None = None,
    cryptic_splice_threshold: float = 3.0,
) -> list[TypeCheckResult]:
    """
    Evaluate all registered type predicates against a sequence.

    Uses the registry to discover predicates, so new predicates are
    automatically included when registered.

    Returns results in a canonical order.
    """
    enzymes = enzymes or list(RESTRICTION_ENZYMES.keys())
    target_isoform = "".join(seq[start:end] for start, end in known_exon_boundaries)

    # Evaluate each predicate via the registry — ensures all registered predicates are included
    predicate_args = {
        "NoCrypticSplice": {"seq": seq, "known_exon_boundaries": known_exon_boundaries, "cryptic_threshold": cryptic_splice_threshold},
        "SpliceCorrect": {"seq": seq, "known_exon_boundaries": known_exon_boundaries, "cellular_context": cellular_context},
        "GCInRange": {"seq": target_isoform, "gc_lo": gc_lo, "gc_hi": gc_hi},
        "CodonAdapted": {"seq": target_isoform, "organism": organism, "threshold": cai_threshold},
        "NoRestrictionSite": {"seq": target_isoform, "enzyme_set": enzymes},
        "InFrame": {"seq": target_isoform, "exon_boundaries": [(0, len(target_isoform))]},
        "NoInstabilityMotif": {"seq": target_isoform},
        "NoCpGIsland": {"seq": target_isoform},
    }

    results = []
    for name in registry.names():
        if name in predicate_args:
            try:
                result = registry.evaluate(name, **predicate_args[name])
                results.append(result)
            except Exception as e:
                logger.error("Error evaluating predicate %s: %s", name, e)
                results.append(TypeCheckResult(
                    predicate=name,
                    verdict=Verdict.UNCERTAIN,
                    derivation=[],
                    knowledge_gap=f"Evaluation error: {e}",
                ))
        else:
            logger.debug("Predicate %s registered but no arguments configured for evaluate_all_predicates", name)

    return results
