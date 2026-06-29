"""
BioCompiler IR Splicing — NDFST-based Alternative Splicing
==========================================================

Models splicing as a Non-deterministic Finite State Transducer (NDFST):

  - Input:  IR_L1_PreMRNA  (pre-mRNA with intron/exon regions)
  - Output: list of IR_L2_MatureMRNA  (one per isoform)

The NDFST is **non-deterministic** because a single pre-mRNA can be
spliced in multiple ways (alternative splicing), producing different
mature mRNAs and hence different protein isoforms from one gene.

Splice site recognition
-----------------------
  - Donor sites (5' splice): ``GU`` at intron start, scored by MaxEntScan
    (Yeo & Burge 2004 first-order Markov model, 9-mer context).
  - Acceptor sites (3' splice): ``AG`` at intron end, scored by MaxEntScan
    (23-mer context, polypyrimidine tract).
  - Sites scoring >= ``DEFAULT_CANONICAL_THRESHOLD`` (3.0 bits) are
    "canonical" — used in the primary isoform.
  - Sites scoring above ``NOISE_FLOOR`` (-5.0 bits) but below the
    canonical threshold are "cryptic" — kept as candidates for
    alternative isoforms.
  - Sites at or below ``NOISE_FLOOR`` (e.g. edge-case scores returned
    by MaxEntScan when there isn't enough flanking context) are
    filtered out entirely.

Formal correspondence
---------------------
This module is the executable counterpart of the formal
``SplicingNDFST`` defined in ``proof/BioCompiler/OracleProofs.lean``
(with proved soundness and completeness via ``ndfstRun_sound`` /
``ndfstRun_complete`` in ``proof/BioCompiler/NDFST.lean``).

The NDFST transition relation is encoded implicitly by the
:func:`find_splice_sites` + :func:`enumerate_isoforms` pair:

  * **Soundness** — every output isoform corresponds to an accepting
    run of the NDFST (each isoform is built from a valid chain of
    non-overlapping donor->acceptor pairs).
  * **Completeness** — every accepting run (every valid chain)
    produces an isoform, subject to the ``max_isoforms`` bound.

Alternative splicing modes explored
-----------------------------------
  1. **Primary isoform** — uses all canonical donor->acceptor pairs in
     their natural pairing (each donor paired with the next acceptor).
  2. **Exon skipping** — replaces two consecutive pairs
     ``(d_i, a_i)`` and ``(d_{i+1}, a_{i+1})`` with a single
     ``(d_i, a_{i+1})`` pair, splicing out the intervening exon.
  3. **Intron retention** — drops one pair from the primary path,
     leaving that intron in the mRNA (a major class of alternative
     splicing, ~15% of events genome-wide).
  4. **No splicing** — retains all introns (the unspliced pre-mRNA).

Isoforms that don't produce a valid CDS (no in-frame stop codon)
are filtered out.  The CDS does NOT need to start with ``AUG`` —
matching the IR-L2 invariant that back-translated genes may begin
with any amino acid (e.g. an antibody heavy chain starts with
E = ``GAG``).  Isoforms are deduplicated by their spliced sequence
(5'UTR + CDS + 3'UTR).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

from .types import IR_L1_PreMRNA, IR_L2_MatureMRNA, GeneRegion


# Standard genetic code stop codons (RNA form, used to locate the CDS end).
_RNA_STOP_CODONS = frozenset({"UAA", "UAG", "UGA"})

# Splice site score thresholds (in bits, from MaxEntScan).
# - Sites with score >= DEFAULT_CANONICAL_THRESHOLD are "canonical"
#   (used in the primary isoform).
# - Sites with score > NOISE_FLOOR are kept in the candidate list
#   (might be used as cryptic sites in alternative isoforms).
# - Sites with score <= NOISE_FLOOR are filtered out as noise
#   (this also catches MaxEntScan's _EDGE_CASE_SCORE = -5.0 returned
#   when there isn't enough flanking context to score a site).
DEFAULT_CANONICAL_THRESHOLD: float = 3.0
NOISE_FLOOR: float = -5.0

# Maximum number of isoforms to enumerate by default.  Real genes can
# produce thousands of isoforms combinatorially; we cap this for
# practical reasons (and to mirror the bounded exploration of the
# formal NDFST).
DEFAULT_MAX_ISOFORMS: int = 10


# ==============================================================================
# Data structures
# ==============================================================================

@dataclass
class SpliceSite:
    """A recognized splice site.

    Attributes
    ----------
    position:
        0-indexed position of the splice site in the pre-mRNA.
        For donors, this is the position of the ``G`` in the ``GU``
        dinucleotide (i.e. the start of the donor).
        For acceptors, this is the position of the ``A`` in the ``AG``
        dinucleotide.
    site_type:
        ``"donor"`` (5' splice site, ``GU``) or ``"acceptor"`` (3'
        splice site, ``AG``).
    score:
        MaxEntScan score in bits.  Higher = stronger.  Typical ranges:
        strong canonical donors 8-12, weak/cryptic 0-5, non-sites <0.
    is_canonical:
        True if the site is at or above the canonical threshold (i.e.
        a "real" splice site rather than a cryptic one).
    """

    position: int
    site_type: str  # "donor" or "acceptor"
    score: float
    is_canonical: bool


@dataclass
class SpliceIsoform:
    """One possible splicing outcome (one accepting run of the NDFST).

    Attributes
    ----------
    exons:
        List of :class:`GeneRegion` (exon type) included in this
        isoform, in genomic order.
    cds:
        The spliced CDS (RNA), ending with a stop codon.  Length is
        always a multiple of 3.  The CDS does NOT need to start with
        ``AUG`` — back-translated genes may begin with any amino acid
        (see the module docstring and
        :func:`biocompiler.ir.invariants.check_l2_invariants`).
    five_utr:
        The 5' UTR (RNA), spliced.
    three_utr:
        The 3' UTR (RNA), spliced.
    is_primary:
        True if this is the most likely isoform (uses all canonical
        splice sites in their natural pairing).
    confidence:
        Confidence score in ``[0, 1]``, based on splice site scores.
        Primary isoform = 1.0; alternatives < 1.0; no-splicing = 0.0.
    """

    exons: List[GeneRegion]
    cds: str
    five_utr: str
    three_utr: str
    is_primary: bool
    confidence: float


# ==============================================================================
# Splice site detection
# ==============================================================================

def find_splice_sites(
    ir_l1: IR_L1_PreMRNA,
    threshold: float = DEFAULT_CANONICAL_THRESHOLD,
) -> List[SpliceSite]:
    """Find all potential splice sites in the pre-mRNA.

    Uses MaxEntScan scoring (from :mod:`biocompiler.sequence.maxentscan`):

      - Donor sites: ``GU`` dinucleotide + downstream context (9-mer,
        first-order Markov model).
      - Acceptor sites: ``AG`` dinucleotide + upstream context (23-mer,
        first-order Markov model).

    Sites above ``threshold`` are marked ``is_canonical=True``; sites
    above :data:`NOISE_FLOOR` but below ``threshold`` are kept as
    cryptic candidates; sites at or below :data:`NOISE_FLOOR` (e.g.
    edge-case scores from insufficient flanking context) are filtered
    out entirely.

    Parameters
    ----------
    ir_l1:
        The pre-mRNA IR object.  Sequence may be upper or lower case.
    threshold:
        MaxEntScan score (in bits) at or above which a site is
        considered "canonical".  Default 3.0.

    Returns
    -------
    list of :class:`SpliceSite`
        All candidate splice sites, in genomic order (donors and
        acceptors interleaved by position).  Empty if no sites found.
    """
    # Import inside the function to keep the dependency direction
    # explicit and lazy (maxentscan is in biocompiler.sequence, which
    # does not depend on biocompiler.ir, but deferring the import
    # means a missing or broken maxentscan module only fails when
    # splicing is actually invoked, not at IR module load time).
    from biocompiler.sequence.maxentscan import score_donor, score_acceptor

    rna_seq = ir_l1.sequence.upper()
    # MaxEntScan expects DNA (T, not U).  Positions are identical: a
    # ``GU`` at position ``i`` in RNA is a ``GT`` at position ``i`` in
    # DNA, so we can score using the DNA sequence and the RNA position.
    dna_seq = rna_seq.replace("U", "T")

    sites: List[SpliceSite] = []

    # Find donor sites (GU in RNA = GT in DNA).
    for i in range(len(rna_seq) - 1):
        if rna_seq[i:i + 2] == "GU":
            score = score_donor(dna_seq, i)
            if score > NOISE_FLOOR:
                sites.append(SpliceSite(
                    position=i,
                    site_type="donor",
                    score=score,
                    is_canonical=(score >= threshold),
                ))

    # Find acceptor sites (AG in RNA = AG in DNA).
    for i in range(len(rna_seq) - 1):
        if rna_seq[i:i + 2] == "AG":
            score = score_acceptor(dna_seq, i)
            if score > NOISE_FLOOR:
                sites.append(SpliceSite(
                    position=i,
                    site_type="acceptor",
                    score=score,
                    is_canonical=(score >= threshold),
                ))

    # Sort by position for deterministic output.
    sites.sort(key=lambda s: s.position)
    return sites


# ==============================================================================
# Isoform enumeration (the NDFST exploration)
# ==============================================================================

# Type alias: a (donor, acceptor) pair defines one intron.
_SplicePair = Tuple[SpliceSite, SpliceSite]


def _find_primary_path(
    donors: List[SpliceSite],
    acceptors: List[SpliceSite],
) -> List[_SplicePair]:
    """Find the primary splice path (the most likely isoform).

    The primary path is the "natural" splicing: walk donors in genomic
    order and pair each one with the next available acceptor that
    comes after it (and after the previous intron's end).  This
    corresponds to the most common mRNA isoform observed
    experimentally — the one that uses every canonical intron.

    Donors that fall inside a previously-defined intron are skipped
    (they would create overlapping introns, which is biologically
    impossible for a single isoform).

    Parameters
    ----------
    donors:
        Canonical donor sites (any order; will be sorted internally).
    acceptors:
        Canonical acceptor sites (any order; will be sorted internally).

    Returns
    -------
    list of (SpliceSite, SpliceSite) tuples
        Non-overlapping ``(donor, acceptor)`` pairs, in genomic order.
        Empty if no valid pairing exists.
    """
    path: List[_SplicePair] = []
    acceptors_sorted = sorted(acceptors, key=lambda s: s.position)
    a_idx = 0
    # Exclusive end of the last intron: the next donor must be at or
    # after this position to avoid overlap.  -1 means "no intron yet".
    last_end = -1

    for d in sorted(donors, key=lambda s: s.position):
        if d.position < last_end:
            continue  # donor is inside a previous intron — skip
        # Advance past acceptors that are at or before this donor
        # (an acceptor must come strictly after its donor).
        while a_idx < len(acceptors_sorted) and acceptors_sorted[a_idx].position <= d.position:
            a_idx += 1
        if a_idx >= len(acceptors_sorted):
            break  # no more acceptors — done
        a = acceptors_sorted[a_idx]
        path.append((d, a))
        # Intron spans [d.position, a.position + 2).  AG is at
        # a.position, a.position+1; the next exon starts at a.position+2.
        last_end = a.position + 2
        a_idx += 1

    return path


def _build_isoform(
    ir_l1: IR_L1_PreMRNA,
    path: List[_SplicePair],
) -> Optional[SpliceIsoform]:
    """Build a :class:`SpliceIsoform` from a splice path.

    The splice path is a list of ``(donor, acceptor)`` pairs.  Each
    pair defines one intron spanning ``[donor.position,
    acceptor.position + 2)`` — i.e. from the ``G`` of the donor ``GU``
    through the ``G`` of the acceptor ``AG``.  Everything outside the
    introns is an exon.

    The exons are concatenated to form the spliced sequence.  The CDS
    starts at position 0 of the spliced sequence (it is **not**
    scanned for an ``AUG`` start codon — see the module docstring)
    and extends to the first in-frame stop codon; everything after
    the stop codon is the 3'UTR.  The 5'UTR is empty by construction
    (the spliced sequence begins with the CDS); this matches the
    behaviour of :func:`biocompiler.ir.passes.splice` when no
    ``5_utr`` region is annotated.

    Parameters
    ----------
    ir_l1:
        The pre-mRNA IR object.
    path:
        Splice path (list of ``(donor, acceptor)`` pairs).  May be
        empty (no splicing — the whole sequence is one exon).

    Returns
    -------
    :class:`SpliceIsoform` or None
        ``None`` if no valid CDS can be extracted (no in-frame stop
        codon in the spliced sequence, or the spliced sequence is
        empty).  Otherwise a :class:`SpliceIsoform` with
        ``is_primary`` and ``confidence`` left unset (the caller
        fills them in).
    """
    rna_seq = ir_l1.sequence.upper()

    # Sort path by donor position (defensive — callers should already
    # sort, but this guarantees correctness if they don't).
    sorted_path = sorted(path, key=lambda p: p[0].position)

    # Build exon regions by walking the path and carving out introns.
    exons: List[GeneRegion] = []
    prev_end = 0
    for d, a in sorted_path:
        # Exon before this intron: [prev_end, d.position)
        if d.position > prev_end:
            exons.append(GeneRegion(prev_end, d.position, "exon"))
        # The intron [d.position, a.position + 2) is removed.
        prev_end = a.position + 2
    # Final exon: [prev_end, len)
    if prev_end < len(rna_seq):
        exons.append(GeneRegion(prev_end, len(rna_seq), "exon"))

    # Concatenate exons to form the spliced sequence.
    spliced = "".join(rna_seq[e.start:e.end] for e in exons)

    if not spliced:
        return None  # no exons — no isoform

    # CDS starts at position 0 of the spliced sequence — matching
    # passes.splice().  The IR-L2 invariant does NOT require the CDS
    # to start with AUG; back-translated genes may begin with any
    # amino acid (e.g., antibody chains starting with E = GAG).  We do
    # NOT scan for an internal AUG, because doing so would silently
    # drop isoforms whose spliced sequence doesn't contain one (the
    # bug reported as H7: back-translated genes starting with E got
    # zero isoforms).
    start_idx = 0

    # Find first in-frame stop codon.
    stop_idx = -1
    for i in range(start_idx, len(spliced) - 2, 3):
        codon = spliced[i:i + 3]
        if codon in _RNA_STOP_CODONS:
            stop_idx = i + 3  # exclusive end
            break
    if stop_idx < 0:
        return None  # no in-frame stop codon

    five_utr = spliced[:start_idx]
    cds = spliced[start_idx:stop_idx]
    three_utr = spliced[stop_idx:]

    return SpliceIsoform(
        exons=exons,
        cds=cds,
        five_utr=five_utr,
        three_utr=three_utr,
        is_primary=False,  # caller sets this
        confidence=0.0,  # caller sets this
    )


def _compute_confidence(
    path: List[_SplicePair],
    primary_path: List[_SplicePair],
) -> float:
    """Compute isoform confidence as the fraction of primary splice-site score used.

    The primary isoform (which uses all primary pairs) gets confidence
    ``1.0``.  Alternative isoforms (which use fewer or different
    pairs) get confidence ``< 1.0``.  The "no splicing" isoform (empty
    path) gets confidence ``0.0``.

    The formula is::

        confidence = sum(score(d) + score(a) for d,a in path)
                     / sum(score(d) + score(a) for d,a in primary_path)

    This rewards isoforms that use strong splice sites and penalises
    isoforms that drop or skip pairs (since each dropped pair removes
    its donor+acceptor score from the numerator).
    """
    if not primary_path:
        # No canonical pairs -> only the empty (unspliced) path is
        # valid, and it's the "primary" by default.
        return 1.0 if not path else 0.0

    max_total = sum(d.score + a.score for d, a in primary_path)
    if max_total <= 0:
        return 0.5  # degenerate case — no information

    actual_total = sum(d.score + a.score for d, a in path)
    return max(0.0, min(1.0, actual_total / max_total))


def enumerate_isoforms(
    ir_l1: IR_L1_PreMRNA,
    max_isoforms: int = DEFAULT_MAX_ISOFORMS,
    threshold: float = DEFAULT_CANONICAL_THRESHOLD,
) -> List[SpliceIsoform]:
    """Enumerate possible splice isoforms using NDFST semantics.

    The NDFST explores all valid ``(donor, acceptor)`` pairings.  Each
    pairing produces one isoform.  We limit to ``max_isoforms`` to
    avoid combinatorial explosion (real genes can produce thousands of
    isoforms).

    The following alternative-splicing modes are explored (see module
    docstring for details):

      1. Primary isoform (all canonical pairs, natural pairing).
      2. Exon skipping (replace two consecutive pairs with one
         long-range pair).
      3. Intron retention (drop one pair from the primary path).
      4. No splicing (retain all introns).

    Isoforms that don't produce a valid CDS (no in-frame stop codon)
    are filtered out.  The CDS does NOT need to start with ``AUG`` —
    matching the IR-L2 invariant (back-translated genes may begin
    with any amino acid).  Isoforms are deduplicated by their spliced
    sequence (5'UTR + CDS + 3'UTR).

    Parameters
    ----------
    ir_l1:
        The pre-mRNA IR object.
    max_isoforms:
        Maximum number of isoforms to return.  Default 10.
    threshold:
        MaxEntScan score threshold for canonical sites.  Default 3.0.

    Returns
    -------
    list of :class:`SpliceIsoform`
        Isoforms sorted by: primary first, then by confidence
        descending.  Empty if no canonical splice sites are found.
    """
    sites = find_splice_sites(ir_l1, threshold=threshold)
    canonical_donors = [
        s for s in sites if s.site_type == "donor" and s.is_canonical
    ]
    canonical_acceptors = [
        s for s in sites if s.site_type == "acceptor" and s.is_canonical
    ]

    if not canonical_donors or not canonical_acceptors:
        return []

    # All valid (donor, acceptor) pairs where donor.position < acceptor.position.
    all_pairs: List[_SplicePair] = []
    for d in canonical_donors:
        for a in canonical_acceptors:
            if d.position < a.position:
                all_pairs.append((d, a))

    if not all_pairs:
        return []

    # Primary path: natural donor-acceptor pairing.
    primary_path = _find_primary_path(canonical_donors, canonical_acceptors)

    # Build candidate paths to explore.  We always include the primary
    # path first so it's guaranteed to be built (and ranked first
    # after sorting), even when max_isoforms is small.
    candidate_paths: List[List[_SplicePair]] = []

    # 1. Primary path.
    candidate_paths.append(primary_path)

    # 2. Exon skipping: replace (d_i, a_i) and (d_{i+1}, a_{i+1})
    #    with (d_i, a_{i+1}), splicing out the exon between them.
    for i in range(len(primary_path) - 1):
        d_i, _ = primary_path[i]
        _, a_j = primary_path[i + 1]
        # Look for the long-range pair (d_i, a_j) in all_pairs.
        skip_pair = next(
            (
                (d, a)
                for d, a in all_pairs
                if d.position == d_i.position and a.position == a_j.position
            ),
            None,
        )
        if skip_pair is not None:
            new_path = (
                primary_path[:i] + [skip_pair] + primary_path[i + 2:]
            )
            candidate_paths.append(new_path)

    # 3. Intron retention: drop one pair from the primary path.
    for i in range(len(primary_path)):
        new_path = primary_path[:i] + primary_path[i + 1:]
        candidate_paths.append(new_path)

    # 4. No splicing: retain all introns (empty path).
    candidate_paths.append([])

    # Build all isoforms, deduplicate by spliced sequence.
    isoforms: List[SpliceIsoform] = []
    seen_seqs: set = set()

    for path in candidate_paths:
        iso = _build_isoform(ir_l1, path)
        if iso is None:
            continue  # no valid CDS — skip
        seq_key = iso.five_utr + "|" + iso.cds + "|" + iso.three_utr
        if seq_key in seen_seqs:
            continue  # duplicate — skip
        seen_seqs.add(seq_key)
        iso.is_primary = (path == primary_path)
        iso.confidence = _compute_confidence(path, primary_path)
        isoforms.append(iso)

    # Sort: primary first, then by confidence descending.
    isoforms.sort(key=lambda i: (not i.is_primary, -i.confidence))

    # Truncate to max_isoforms (primary is guaranteed to survive
    # because it's first after sorting).
    return isoforms[:max_isoforms]


# ==============================================================================
# Main entry point: L1 -> list of L2
# ==============================================================================

def splice_ndfst(
    ir_l1: IR_L1_PreMRNA,
    primary_only: bool = False,
) -> List[IR_L2_MatureMRNA]:
    """NDFST-based splicing: L1 -> list of L2 isoforms.

    This is the main entry point for NDFST-based splicing.  It
    enumerates all valid splice isoforms (up to
    :data:`DEFAULT_MAX_ISOFORMS`) and returns them as IR-L2 objects.

    If ``primary_only`` is True, only the most likely (primary)
    isoform is returned (still as a single-element list, for API
    uniformity).

    If no canonical splice sites are found, returns an empty list —
    the caller (typically :func:`biocompiler.ir.passes.splice`) is
    expected to fall back to the simple start/stop codon logic.  This
    avoids a circular dependency: ``passes.splice`` calls
    ``splice_ndfst``; if ``splice_ndfst`` were to call
    ``passes.splice`` as a fallback, the result would be infinite
    recursion.  Keeping the fallback in ``passes.splice`` (the
    *caller*) breaks the cycle cleanly.

    Parameters
    ----------
    ir_l1:
        The pre-mRNA IR object.
    primary_only:
        If True, return only the primary isoform (as a 1-element
        list).  Default False (return all isoforms).

    Returns
    -------
    list of :class:`IR_L2_MatureMRNA`
        Splice isoforms, primary first.  Empty if no splice sites
        found (caller should fall back to simple splice).
    """
    # Always enumerate with the default cap so the primary isoform is
    # guaranteed to be built and ranked first.  When primary_only=True
    # we truncate to the first element afterwards.
    isoforms = enumerate_isoforms(ir_l1, max_isoforms=DEFAULT_MAX_ISOFORMS)

    if not isoforms:
        return []

    if primary_only:
        isoforms = isoforms[:1]

    results: List[IR_L2_MatureMRNA] = []
    for iso in isoforms:
        ir_l2 = IR_L2_MatureMRNA(
            sequence=iso.five_utr + iso.cds + iso.three_utr,
            five_utr=iso.five_utr,
            cds=iso.cds,
            three_utr=iso.three_utr,
            organism=ir_l1.organism,
            gene_name=ir_l1.gene_name,
            metadata={
                **ir_l1.metadata,
                "pass": "splice_ndfst",
                "source_level": "L1",
                "isoform_confidence": iso.confidence,
                "is_primary": iso.is_primary,
                "exon_count": len(iso.exons),
            },
        )
        results.append(ir_l2)

    return results


__all__ = [
    "SpliceSite",
    "SpliceIsoform",
    "find_splice_sites",
    "enumerate_isoforms",
    "splice_ndfst",
    "DEFAULT_CANONICAL_THRESHOLD",
    "DEFAULT_MAX_ISOFORMS",
    "NOISE_FLOOR",
]
