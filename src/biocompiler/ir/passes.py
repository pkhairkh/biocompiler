"""
BioCompiler IR — Lowering Passes
================================

Each pass transforms one IR level to the next, matching the Central
Dogma:

    L0 → L1   ``transcribe``   DNA → pre-mRNA (T → U, regions preserved)
    L1 → L2   ``splice``       remove introns, assemble 5'UTR + CDS + 3'UTR
    L2 → L3   ``translate``    codon → amino acid (standard genetic code)
    L3 → L4   ``fold``         protein folding oracle (ESMFold/AlphaFold — Phase 2)

Design rules
------------
* Each pass is a **pure function**: same input → same output, no side
  effects, no I/O.  This makes passes trivially testable and composable.
* Each pass validates its **input** (raising :class:`IRError` if
  invariants are violated) and constructs only **well-formed** output
  (raising :class:`IRError` if it cannot).
* Each pass stamps ``metadata["pass"]`` and ``metadata["source_level"]``
  so the resulting IR object remembers how it was produced — this is
  the provenance trail that makes the IR auditable.

Phase 1 scope
-------------
The ``transcribe``, ``splice``, and ``translate`` passes are fully
implemented as pure functions.  The ``fold`` pass delegates to the
folding oracle in :mod:`biocompiler.ir.folding` — it tries ESMFold
first (real 3-D structure prediction) and falls back to a pure-Python
Chou-Fasman heuristic when ESMFold is unavailable.  The resulting
:class:`IR_L4_FoldedProtein` records which oracle was used in
``metadata["oracle"]``.
"""

from __future__ import annotations

from typing import Union

from .types import (
    IR_L0_GenomicDNA,
    IR_L1_PreMRNA,
    IR_L2_MatureMRNA,
    IR_L3_Polypeptide,
    IR_L4_FoldedProtein,
    IRError,
    IRLevel,
    GeneRegion,
)
from ..type_system.codon_tables import CODON_TABLE

# Re-export the union of all IR object types for type hints.
IRObject = Union[
    IR_L0_GenomicDNA,
    IR_L1_PreMRNA,
    IR_L2_MatureMRNA,
    IR_L3_Polypeptide,
    IR_L4_FoldedProtein,
]

# Standard genetic code stop codons (RNA form, used by ``splice``).
_RNA_STOP_CODONS = frozenset({"UAA", "UAG", "UGA"})


def transcribe(ir_l0: IR_L0_GenomicDNA) -> IR_L1_PreMRNA:
    """L0 → L1: Transcription.

    DNA → pre-mRNA: replace every ``T`` with ``U`` and copy all region
    annotations unchanged (transcription does not move sequence
    positions).  This is essentially identity at the sequence level.

    Raises
    ------
    IRError
        If the input sequence is empty or contains bases outside
        ``{A, C, G, T, N}``.
    """
    # --- validate input -------------------------------------------------
    if not ir_l0.sequence:
        raise IRError(IRLevel.L0, "empty sequence")
    valid_bases = set("ACGTN")
    seq_upper = ir_l0.sequence.upper()
    if not set(seq_upper) <= valid_bases:
        bad = set(seq_upper) - valid_bases
        raise IRError(IRLevel.L0, f"invalid bases in sequence: {sorted(bad)}")

    # --- transcribe: DNA → RNA (T → U) ---------------------------------
    rna_sequence = seq_upper.replace("T", "U")

    # --- copy regions (same positions in pre-mRNA) ---------------------
    regions = [
        GeneRegion(r.start, r.end, r.region_type, dict(r.metadata))
        for r in ir_l0.regions
    ]

    return IR_L1_PreMRNA(
        sequence=rna_sequence,
        regions=regions,
        organism=ir_l0.organism,
        gene_name=ir_l0.gene_name,
        secis_positions=list(ir_l0.secis_positions),  # propagate SECIS positions
        # NOTE: spread the upstream metadata FIRST so the new "pass" /
        # "source_level" keys below take precedence over any inherited
        # ones (otherwise the previous pass's "pass" value would survive
        # and the provenance trail would be wrong).
        metadata={
            **ir_l0.metadata,
            "pass": "transcribe",
            "source_level": "L0",
        },
    )


def splice(ir_l1: IR_L1_PreMRNA) -> IR_L2_MatureMRNA:
    """L1 → L2: Splicing.

    Remove introns and assemble the three functional sub-regions of the
    mature mRNA: 5'UTR, CDS, 3'UTR.

    Two input shapes are supported:

    1. **No region annotations** (``ir_l1.regions == []``): the whole
       sequence is treated as the CDS.  The CDS starts at position 0
       (it is **not** scanned for an ``AUG`` start codon — see the
       "Start-codon policy" note below) and extends to the first
       in-frame stop codon; everything after the stop codon is the
       3'UTR.  No 5'UTR is extracted (there is no annotation to mark
       where it would end and the CDS would begin).
    2. **Region annotations present**: ``exon`` and ``cds`` regions are
       concatenated (sorted by ``start``) to form the spliced
       transcript; ``intron`` and other non-exon regions are dropped.
       As in case (1), the CDS starts at position 0 of the spliced
       transcript (no ``AUG`` scan) and extends to the first in-frame
       stop codon; everything after the stop codon is the 3'UTR.  A
       separate ``5_utr`` region annotation, if present, is extracted
       as the 5'UTR (concatenated in genomic order if multiple
       ``5_utr`` regions exist); otherwise the 5'UTR is empty.

    Start-codon policy
    ------------------
    The CDS does **not** need to start with ``AUG``.  This matches the
    IR-L2 invariant (see
    :func:`biocompiler.ir.invariants.check_l2_invariants`):
    back-translated genes may begin with any amino acid — e.g. an
    antibody heavy chain starts with E = ``GAG``, not M = ``AUG``.
    Scanning for an internal ``AUG`` would silently truncate such
    genes, so the splicing pass never does so.  Callers that need
    start-codon semantics should annotate a ``5_utr`` region ending at
    the ``AUG``; the CDS will then begin at the first base of the
    first ``exon``/``cds`` region (which the caller is expected to
    place immediately after the 5'UTR).

    SECIS-aware stop-codon scan
    ---------------------------
    A ``UGA`` codon at a position listed in ``ir_l1.secis_positions``
    is interpreted as selenocysteine (not stop) and is skipped by the
    stop-codon scan, so selenoprotein CDSs terminate at the next
    non-SECIS stop codon.

    Raises
    ------
    IRError
        If no in-frame stop codon is found, or if the resulting CDS
        violates IR-L2 invariants (must end with a stop codon and be
        a multiple of 3 long).  The CDS does **not** need to start
        with ``AUG`` — see the "Start-codon policy" note above.
    """
    rna = ir_l1.sequence.upper()

    if not ir_l1.regions:
        # Case 1: whole sequence is the input to the start/stop scan.
        spliced = rna
    else:
        # Case 2: concatenate exons (and explicit 'cds' regions) to
        # form the spliced transcript.  Introns are dropped silently.
        exons = sorted(
            (
                r
                for r in ir_l1.regions
                if r.region_type in ("exon", "cds")
            ),
            key=lambda r: r.start,
        )
        if not exons:
            raise IRError(
                IRLevel.L1, "no exon/cds regions found in region annotations"
            )
        spliced = "".join(rna[r.start:r.end] for r in exons)

    # --- locate CDS boundaries ----------------------------------------
    # When no regions: CDS starts at position 0 (the whole sequence IS the
    # CDS). We do NOT search for an internal AUG, because back-translated
    # genes may not start with M (e.g., antibody chains start with E = GAG).
    # When regions present: exons are already concatenated; CDS starts at
    # position 0 of the spliced transcript.
    #
    # The CDS extends to the first in-frame stop codon (scanning from 0
    # in steps of 3).
    start_idx = 0

    # --- locate first in-frame stop codon ------------------------------
    # SECIS-aware: a UGA codon at a SECIS position is selenocysteine, not stop.
    # Skip such codons when scanning for the terminal stop codon.
    secis_set = set(ir_l1.secis_positions)
    stop_idx = -1
    codon_idx = 0
    for i in range(start_idx, len(spliced) - 2, 3):
        codon = spliced[i:i + 3]
        if codon in _RNA_STOP_CODONS:
            # Check if this UGA is a SECIS position (selenocysteine, not stop)
            if codon == "UGA" and codon_idx in secis_set:
                codon_idx += 1
                continue  # skip — this is selenocysteine, not stop
            stop_idx = i + 3  # exclusive end
            break
        codon_idx += 1
    if stop_idx < 0:
        raise IRError(
            IRLevel.L1, "no stop codon (UAA/UAG/UGA) found in-frame"
        )

    # Extract the 5'UTR from an explicit ``5_utr`` region annotation,
    # if one is present.  When no 5'UTR region is annotated, the
    # 5'UTR is empty by construction: the spliced transcript begins
    # with the CDS (the first ``exon``/``cds`` region), so there is no
    # leading sequence to attribute to the 5'UTR.  Callers that need
    # start-codon semantics must annotate a ``5_utr`` region ending at
    # the AUG; see the splice() docstring's "Start-codon policy".
    five_utr = ""
    if ir_l1.regions:
        five_utr_regions = sorted(
            (r for r in ir_l1.regions if r.region_type == "5_utr"),
            key=lambda r: r.start,
        )
        if five_utr_regions:
            five_utr = "".join(rna[r.start:r.end] for r in five_utr_regions)
    cds = spliced[start_idx:stop_idx]
    three_utr = spliced[stop_idx:]

    # --- validate IR-L2 invariants before constructing ----------------
    # CDS does NOT need to start with AUG — back-translated genes may
    # start with any amino acid. It just needs to end with a stop codon
    # and be a multiple of 3.
    if not cds:
        raise IRError(IRLevel.L2, "empty CDS")
    if cds[-3:] not in _RNA_STOP_CODONS:
        raise IRError(
            IRLevel.L2, f"CDS doesn't end with stop codon: {cds[-3:]!r}"
        )
    if len(cds) % 3 != 0:
        raise IRError(
            IRLevel.L2, f"CDS length not multiple of 3: {len(cds)}"
        )

    return IR_L2_MatureMRNA(
        sequence=five_utr + cds + three_utr,
        five_utr=five_utr,
        cds=cds,
        three_utr=three_utr,
        organism=ir_l1.organism,
        gene_name=ir_l1.gene_name,
        secis_positions=list(ir_l1.secis_positions),  # propagate SECIS positions
        metadata={
            **ir_l1.metadata,
            "pass": "splice",
            "source_level": "L1",
        },
    )


def translate(ir_l2: IR_L2_MatureMRNA) -> IR_L3_Polypeptide:
    """L2 → L3: Translation.

    Translate the CDS codon-by-codon using the standard genetic code
    (:data:`biocompiler.type_system.codon_tables.CODON_TABLE`).

    The CDS is expected to satisfy IR-L2 invariants (start with AUG,
    end with a stop codon, length divisible by 3).  Each codon is
    translated to its single-letter amino acid; unknown codons (e.g.
    containing ``N``) become ``"X"``.  The trailing stop codon is
    translated to ``"*"``.

    Raises
    ------
    IRError
        If the CDS is empty, has length not divisible by 3, or does
        not start with ``AUG``.
    """
    cds = ir_l2.cds.upper()

    if not cds:
        raise IRError(IRLevel.L2, "empty CDS")
    if len(cds) % 3 != 0:
        raise IRError(
            IRLevel.L2,
            f"CDS length not multiple of 3: {len(cds)}",
        )
    # CODON_TABLE is keyed by DNA codons (T, not U), so convert.
    # CDS does NOT need to start with AUG — back-translated genes may
    # start with any amino acid.
    #
    # Selenocysteine (U): the codon UGA (TGA in DNA) is normally a stop
    # codon, but in selenoproteins a SECIS element in the 3'UTR causes
    # the ribosome to insert selenocysteine instead. The secis_positions
    # field lists CDS codon indices where UGA should be recoded to U.
    secis_set = set(ir_l2.secis_positions)
    protein_chars: list[str] = []
    codon_idx = 0
    for i in range(0, len(cds), 3):
        rna_codon = cds[i:i + 3]
        dna_codon = rna_codon.replace("U", "T")
        if dna_codon == "TGA" and codon_idx in secis_set:
            protein_chars.append("U")  # selenocysteine (Sec)
        else:
            protein_chars.append(CODON_TABLE.get(dna_codon, "X"))
        codon_idx += 1

    protein_seq = "".join(protein_chars)

    return IR_L3_Polypeptide(
        sequence=protein_seq,
        organism=ir_l2.organism,
        gene_name=ir_l2.gene_name,
        metadata={
            **ir_l2.metadata,
            "pass": "translate",
            "source_level": "L2",
        },
    )


def fold(
    ir_l3: IR_L3_Polypeptide,
    use_esmfold: bool = True,
) -> IR_L4_FoldedProtein:
    """L3 → L4: Protein folding via ESMFold oracle (with heuristic fallback).

    This is an ORACLE-DEPENDENT pass — the structure prediction comes
    from an external ML model (ESMFold), not a proven computation.
    The :class:`IR_L4_FoldedProtein` object records which oracle was
    used (``metadata["oracle"]``) and the achieved confidence
    (``confidence`` field, mean pLDDT on a 0–100 scale).

    Strategy:
      1. Try ESMFold (ESM Atlas API or locally-installed ``esm``
         package) — produces real 3-D coordinates + per-residue pLDDT.
      2. On any ESMFold failure (offline, OOM on long sequences,
         non-standard amino acids, ...) fall back to the pure-Python
         Chou-Fasman heuristic — produces a secondary-structure string
         and a calibrated low-confidence pLDDT estimate, but no 3-D
         coordinates.

    The trailing ``"*"`` stop codon on the polypeptide sequence is
    stripped before folding (it is not a residue).  The original
    sequence (with ``"*"``) is preserved in ``IR_L4_FoldedProtein.sequence``
    so round-tripping through L3 → L4 → L3 is lossless.

    Parameters
    ----------
    ir_l3:
        The IR-L3 polypeptide to fold.
    use_esmfold:
        If True (default), attempt ESMFold first.  If False, skip
        straight to the heuristic — useful for offline tests, fast
        CI runs, or when the caller knows the ESM Atlas API is
        unreachable.

    Returns
    -------
    IR_L4_FoldedProtein
        Always non-``None``.  The ``metadata["oracle"]`` field records
        which backend produced the result (``"esmfold"`` /
        ``"fallback"`` / ``"none"``).

    Raises
    ------
    IRError
        If the polypeptide sequence is empty (after stripping the
        trailing ``"*"``), or if :func:`fold_sequence` rejects it.
    """
    from .folding import fold_sequence

    result = fold_sequence(ir_l3.sequence, use_esmfold=use_esmfold)

    return IR_L4_FoldedProtein(
        sequence=ir_l3.sequence,
        coordinates=result.coordinates,
        confidence=result.confidence,
        ptms=[],  # PTMs are Phase 3 — not yet implemented.
        organism=ir_l3.organism,
        gene_name=ir_l3.gene_name,
        metadata={
            **ir_l3.metadata,
            "pass": "fold",
            "source_level": "L3",
            "oracle": result.oracle_used,
            "secondary_structure": result.secondary_structure,
            "fold_metadata": result.metadata,
        },
    )


def compile_gene(
    ir_l0: IR_L0_GenomicDNA,
    target_level: IRLevel = IRLevel.L3,
    use_esmfold: bool = True,
) -> IRObject:
    """Run the full compilation pipeline from L0 to ``target_level``.

    The pipeline applies lowering passes in order, stopping once the
    target level is reached::

        L0 ──transcribe──▶ L1 ──splice──▶ L2 ──translate──▶ L3 ──fold──▶ L4

    Comparison of levels uses :attr:`IRLevel.order` (a numerical index),
    **not** the string ``.value`` of the enum member — the latter is
    alphabetically ordered and would yield wrong pipeline skips.

    Parameters
    ----------
    ir_l0:
        The source-level IR-L0 object.
    target_level:
        The IR level to stop at.  Defaults to :attr:`IRLevel.L3`
        (polypeptide), the deepest level with a deterministic (non-
        oracle) lowering pass.
    use_esmfold:
        Forwarded to :func:`fold` if ``target_level`` is
        :attr:`IRLevel.L4`.  If True (default), the fold pass attempts
        ESMFold first and falls back to the heuristic on failure.  If
        False, the fold pass uses the heuristic directly — useful for
        offline tests or fast CI runs.

    Returns
    -------
    IRObject
        An IR object of the requested level (L0 through L4).

    Raises
    ------
    IRError
        If any lowering pass rejects its input.
    ValueError
        If ``target_level`` is not a known :class:`IRLevel` member.
    """
    target_idx = target_level.order  # numerical 0..4
    if target_idx < 0:
        raise ValueError(f"unknown target level: {target_level!r}")

    current: IRObject = ir_l0
    if target_idx >= IRLevel.L1.order:
        current = transcribe(current)  # type: ignore[arg-type]
    if target_idx >= IRLevel.L2.order:
        current = splice(current)  # type: ignore[arg-type]
    if target_idx >= IRLevel.L3.order:
        current = translate(current)  # type: ignore[arg-type]
    if target_idx >= IRLevel.L4.order:
        current = fold(current, use_esmfold=use_esmfold)  # type: ignore[arg-type]
    return current


__all__ = [
    "transcribe",
    "splice",
    "translate",
    "fold",
    "compile_gene",
    "IRObject",
]


def splice_all_isoforms(ir_l1: IR_L1_PreMRNA) -> list[IR_L2_MatureMRNA]:
    """L1 → list of L2: All splice isoforms (NDFST non-determinism).
    
    Returns all possible mature mRNAs from alternative splicing.
    Uses the NDFST-based splicing module.
    """
    from .splicing import splice_ndfst
    return splice_ndfst(ir_l1, primary_only=False)
