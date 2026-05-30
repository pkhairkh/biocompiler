"""
BioCompiler Type System — Predicate Evaluation with Registry Pattern

FIXES from toy model:
- Registry pattern for predicates (extensible, no string-prefix hacks)
- Each predicate is a self-contained callable with metadata
- Easy to add new predicates without modifying verification code
- Proper error handling
"""

import logging
from .types import Verdict, TypeCheckResult, Token
from .scanner import scan_sequence, gc_content, validate_dna_sequence
from .splicing import compute_splice_isoforms
from .translation import compute_cai
from .constants import INSTABILITY_MOTIF, RESTRICTION_ENZYMES, reverse_complement
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
    """

    def __init__(self):
        self._predicates: dict[str, dict] = {}

    def register(self, name: str, evaluator, verifier=None):
        """Register a predicate with its evaluator and optional verifier."""
        self._predicates[name] = {
            "evaluator": evaluator,
            "verifier": verifier or evaluator,
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

def evaluate_no_cryptic_splice(seq: str, known_exon_boundaries: list[tuple[int, int]], **kwargs) -> TypeCheckResult:
    """
    NoCrypticSplice: No donor/acceptor pair within known exons that could form a cryptic intron.
    """
    derivation = []
    for exon_start, exon_end in known_exon_boundaries:
        exon_seq = seq[exon_start:exon_end]
        tokens = scan_sequence(exon_seq)
        donors = [t for t in tokens if t.element_type == "splice_donor"]
        acceptors = [t for t in tokens if t.element_type == "splice_acceptor"]

        for d in donors:
            derivation.append({
                "step": "scan_donor_in_exon",
                "position": exon_start + d.position,
                "sequence": d.match_sequence,
                "score": d.score,
            })
        for a in acceptors:
            derivation.append({
                "step": "scan_acceptor_in_exon",
                "position": exon_start + a.position,
                "sequence": a.match_sequence,
                "score": a.score,
            })

        # Check for donor-acceptor pairs forming potential cryptic intron
        for d in donors:
            for a in acceptors:
                if a.position > d.position + 30:
                    return TypeCheckResult(
                        predicate="NoCrypticSplice",
                        verdict=Verdict.FAIL,
                        derivation=derivation,
                        violation=(
                            f"Cryptic splice site pair in exon [{exon_start},{exon_end}): "
                            f"donor at {exon_start + d.position}, acceptor at {exon_start + a.position}"
                        ),
                    )

    return TypeCheckResult(
        predicate="NoCrypticSplice",
        verdict=Verdict.PASS,
        derivation=derivation + [{"step": "no_cryptic_pairs_found", "evidence": "exhaustive_scan"}],
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

    if len(non_target) == 0 and len(isoforms) == 1:
        return TypeCheckResult(
            predicate=f"SpliceCorrect({cellular_context})",
            verdict=Verdict.PASS,
            derivation=derivation + [{"step": "singleton_isoform_set"}],
        )
    elif non_target:
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
    else:
        return TypeCheckResult(
            predicate=f"SpliceCorrect({cellular_context})",
            verdict=Verdict.UNCERTAIN,
            derivation=derivation,
            knowledge_gap="Multiple isoforms but cannot determine which will occur in vivo",
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
        if enz_name not in RESTRICTION_ENZYMES:
            continue
        site = RESTRICTION_ENZYMES[enz_name]
        site_rc = reverse_complement(site)

        # Forward strand
        for i in range(len(seq) - len(site) + 1):
            if seq[i:i+len(site)] == site:
                found_sites.append({"enzyme": enz_name, "position": i, "site": site, "strand": "+"})

        # Reverse complement strand (skip palindromes)
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


# ==============================================================================
# Register all built-in predicates
# ==============================================================================

registry.register("NoCrypticSplice", evaluate_no_cryptic_splice)
registry.register("SpliceCorrect", evaluate_splice_correct)
registry.register("GCInRange", evaluate_gc_in_range)
registry.register("CodonAdapted", evaluate_codon_adapted)
registry.register("NoRestrictionSite", evaluate_no_restriction_site)
registry.register("InFrame", evaluate_in_frame)
registry.register("NoInstabilityMotif", evaluate_no_instability_motif)


def evaluate_all_predicates(
    seq: str,
    known_exon_boundaries: list[tuple[int, int]],
    organism: str = "Homo_sapiens",
    cellular_context: str = "HEK293T",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.5,
    enzymes: list[str] | None = None,
) -> list[TypeCheckResult]:
    """
    Evaluate all registered type predicates against a sequence.

    Returns results in a canonical order.
    """
    enzymes = enzymes or list(RESTRICTION_ENZYMES.keys())
    target_isoform = "".join(seq[start:end] for start, end in known_exon_boundaries)

    results = [
        evaluate_no_cryptic_splice(seq, known_exon_boundaries),
        evaluate_splice_correct(seq, known_exon_boundaries, cellular_context),
        evaluate_gc_in_range(target_isoform, gc_lo, gc_hi),
        evaluate_codon_adapted(target_isoform, organism, cai_threshold),
        evaluate_no_restriction_site(target_isoform, enzymes),
        evaluate_in_frame(target_isoform, [(0, len(target_isoform))]),
        evaluate_no_instability_motif(target_isoform),
    ]
    return results
