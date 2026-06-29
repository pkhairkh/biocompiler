"""
BioCompiler IR Optimization Passes
==================================

These passes transform IR objects while preserving semantics:

  - ``optimize_codons``        IR_L2 → IR_L2 (different CDS, same protein)
  - ``eliminate_cpgs``         IR_L2 → IR_L2 (remove CpG dinucleotides)
  - ``run_optimization_pipeline``  chain multiple IR_L2 → IR_L2 passes

Each pass:

  1. Checks input invariants (the L2 must already be well-formed).
  2. Translates the L2 CDS back to its protein (semantic identity).
  3. Delegates to the existing ``optimize_sequence`` API — the IR pass
     is a *thin adapter* that wraps the optimizer in IR-level typing
     and provenance, **not** a re-implementation of codon optimization.
  4. Re-attaches the original stop codon (the optimizer returns only
     the protein-coding region, without stop) and converts back to RNA.
  5. Builds a new :class:`IR_L2_MatureMRNA` with the optimized CDS,
     preserving the 5'/3' UTRs, organism, gene name, and upstream
     metadata.
  6. Checks output invariants.
  7. **CRITICAL — protein preservation check**: re-translates both the
     original and the optimized L2 and raises :class:`IRError` if the
     proteins differ.  This is the semantic correctness guarantee —
     analogous to LLVM's ``-verify`` after each pass.

The key correctness property: every optimization pass preserves the
translated protein.  A pass that breaks the protein is a bug, and the
check turns that bug into a loud, immediate failure rather than silent
corruption of the user's gene design.

Design note — IR-level vs. optimizer-level
------------------------------------------
The existing :func:`biocompiler.optimizer.optimize_sequence` operates
on raw strings (protein in, DNA out) and returns an
:class:`OptimizationResult` dataclass with metrics.  It is *unaware*
of the IR — it doesn't know about 5'/3' UTRs, IR-L2 invariants, or
provenance trails.

These IR passes are the bridge: they take a typed IR object, extract
what the optimizer needs (the protein), call the optimizer, and wrap
the result back into a typed IR object that satisfies L2 invariants
and carries a provenance stamp.  This is exactly how LLVM's
``InstCombine`` wraps low-level peephole rewrites: the rewrite logic
is shared, the IR-level pass adds typing and verification.

This mirrors the compiler analogy:

    LLVM:        IR ─ InstCombine ─▶ IR   (same semantics, better code)
    BioCompiler: IR_L2 ─ optimize_codons ─▶ IR_L2  (same protein, better codons)
"""

from __future__ import annotations

from typing import Optional

from .types import IR_L2_MatureMRNA, IRLevel, IRError
from .invariants import check_l2_invariants
from .passes import translate
from ..optimizer.pipeline_core import optimize_sequence
from ..organisms import resolve_organism, CODON_ADAPTIVENESS_TABLES
from ..optimizer.cai import _compute_cai_fast


__all__ = [
    "optimize_codons",
    "eliminate_cpgs",
    "run_optimization_pipeline",
    "ir_cai",
]


# ────────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────────

def ir_cai(ir_l2: IR_L2_MatureMRNA) -> float:
    """Compute the Codon Adaptation Index (CAI) of an IR-L2's CDS.

    Uses the same adaptiveness tables and the same CAI kernel as the
    optimizer, so the value is directly comparable to
    :attr:`OptimizationResult.cai`.  This is the IR-level "objective
    function" that optimization passes aim to improve.

    Returns 0.0 if the organism is unknown or the CDS is empty.
    """
    if not ir_l2.cds:
        return 0.0
    try:
        canonical = resolve_organism(ir_l2.organism, strict=False)
    except Exception:
        return 0.0
    weights = CODON_ADAPTIVENESS_TABLES.get(canonical)
    if weights is None:
        return 0.0
    # IR L2 uses RNA (U); the CAI kernel uses DNA (T).
    dna_cds = ir_l2.cds.upper().replace("U", "T")
    return _compute_cai_fast(dna_cds, weights)


def _validate_input(ir_l2: IR_L2_MatureMRNA) -> None:
    """Reject non-L2 or malformed-L2 inputs before optimizing.

    Raises :class:`IRError` (level L2) with a descriptive message.
    """
    if not isinstance(ir_l2, IR_L2_MatureMRNA):
        raise IRError(
            IRLevel.L2,
            f"optimization passes require IR_L2_MatureMRNA, got {type(ir_l2).__name__}",
        )
    # check_l2_invariants raises IRError on failure — let it propagate.
    check_l2_invariants(ir_l2)


def _protein_for_optimizer(ir_l2: IR_L2_MatureMRNA) -> str:
    """Translate L2 → protein and strip the trailing stop (``*``).

    :func:`optimize_sequence` expects "Amino acid sequence (1-letter
    codes, no stop)" — see its docstring — so the trailing ``*`` from
    :func:`biocompiler.ir.passes.translate` must be removed.  The stop
    codon is preserved separately (see :func:`_re_attach_stop`) so the
    optimized CDS still ends with the *same* stop codon, satisfying the
    IR-L2 invariant that ``cds`` ends with ``UAA``/``UAG``/``UGA``.
    """
    protein = translate(ir_l2).sequence
    if protein.endswith("*"):
        return protein[:-1]
    return protein


def _re_attach_stop(optimized_dna: str, original_l2: IR_L2_MatureMRNA) -> str:
    """Append the original stop codon (as DNA) to the optimized CDS.

    The optimizer's ``result.sequence`` contains only the
    protein-coding region (no stop), so we re-attach the original stop
    codon to keep the IR-L2 invariant (``cds`` ends with a stop codon)
    intact.  Preserving the *original* stop codon choice (rather than
    always using ``TAA``) is a small but real semantic preservation:
    some organisms prefer one stop over another for translational
    termination efficiency.
    """
    # IR L2 CDS is RNA; convert the last 3 nt (the stop codon) to DNA.
    original_stop_dna = original_l2.cds[-3:].upper().replace("U", "T")
    return optimized_dna.upper() + original_stop_dna


def _build_optimized_l2(
    ir_l2: IR_L2_MatureMRNA,
    optimized_dna_with_stop: str,
    *,
    pass_name: str,
    extra_metadata: dict,
) -> IR_L2_MatureMRNA:
    """Construct a new IR_L2 from the optimized DNA CDS.

    Preserves the 5'UTR, 3'UTR, organism, and gene_name from the input.
    Stamps provenance: ``metadata["pass"]`` and ``metadata["source_level"]``
    (matching the convention used by the lowering passes in
    :mod:`biocompiler.ir.passes`), plus pass-specific metrics from
    ``extra_metadata``.
    """
    optimized_cds_rna = optimized_dna_with_stop.replace("T", "U")
    return IR_L2_MatureMRNA(
        sequence=ir_l2.five_utr + optimized_cds_rna + ir_l2.three_utr,
        five_utr=ir_l2.five_utr,
        cds=optimized_cds_rna,
        three_utr=ir_l2.three_utr,
        organism=ir_l2.organism,
        gene_name=ir_l2.gene_name,
        metadata={
            **ir_l2.metadata,
            "pass": pass_name,
            "source_level": "L2",  # IR→IR pass: source and target are both L2
            **extra_metadata,
        },
    )


def _assert_protein_preserved(
    original: IR_L2_MatureMRNA, optimized: IR_L2_MatureMRNA, pass_name: str
) -> None:
    """Re-translate both L2 objects and verify the proteins match.

    This is the **semantic correctness check** — the single most
    important invariant of an IR optimization pass.  If an optimization
    pass changes the protein, it has corrupted the user's gene design;
    we raise :class:`IRError` immediately rather than let the bad
    sequence flow further down the pipeline.
    """
    old_protein = translate(original).sequence
    new_protein = translate(optimized).sequence
    if old_protein != new_protein:
        raise IRError(
            IRLevel.L2,
            f"optimization pass '{pass_name}' changed the translated protein: "
            f"{old_protein!r} → {new_protein!r}",
        )


# ────────────────────────────────────────────────────────────────────
# Public IR optimization passes
# ────────────────────────────────────────────────────────────────────

def optimize_codons(
    ir_l2: IR_L2_MatureMRNA,
    organism: Optional[str] = None,
) -> IR_L2_MatureMRNA:
    """IR optimization pass: codon optimization.

    Transforms the CDS to use optimal codons for the target organism
    while **preserving the translated protein** (semantic correctness).

    This is the compiler's "optimization pass" — analogous to LLVM's
    ``-O2``.  Just as ``-O2`` rewrites instructions to faster
    equivalents that compute the same value, ``optimize_codons``
    rewrites codons to higher-CAI equivalents that encode the same
    amino acid.

    Parameters
    ----------
    ir_l2:
        A well-formed :class:`IR_L2_MatureMRNA` (invariants are checked).
    organism:
        Optional target organism override.  If ``None``, the L2's own
        ``organism`` field is used (so optimising a human HBB L2
        defaults to human codon usage).

    Returns
    -------
    IR_L2_MatureMRNA
        A new IR-L2 with the same protein, same 5'/3' UTRs, same
        organism, same gene_name, but a codon-optimized CDS.  The
        metadata carries:

        * ``pass``         = ``"optimize_codons"``
        * ``source_level`` = ``"L2"`` (IR→IR pass)
        * ``optimization`` = ``"codon_optimization"``
        * ``cai_before``   — CAI of the input CDS
        * ``cai_after``    — CAI of the optimized CDS (from the optimizer)
        * ``cai_delta``    — improvement
        * ``cpg_mode``     — CpG-elimination mode forwarded to
          :func:`biocompiler.optimizer.optimize_sequence`.  Always
          ``"aggressive"`` (see "CpG side-effect" note below).
        * ``cpg_before``   — number of ``CG`` dinucleotides in the
          input CDS.
        * ``cpg_after``    — number of ``CG`` dinucleotides in the
          optimized CDS.
        * ``cpg_removed``  — net reduction in CpG count.

    CpG side-effect
    ----------------
    This pass forwards ``cpg_mode="aggressive"`` to the underlying
    optimizer (the default of :func:`optimize_sequence`).
    "aggressive" CpG elimination runs *as part of* codon optimization,
    not just in the dedicated :func:`eliminate_cpgs` pass.  Removing
    ``CG`` dinucleotides can force the optimizer away from the
    highest-CAI synonymous codon for some amino acids (notably Arg,
    where the highest-CAI codon in many organisms is ``CGN``), so the
    recorded ``cai_after`` may be lower than a pure CAI-maximising
    optimisation would achieve.  The ``cpg_mode`` / ``cpg_before`` /
    ``cpg_after`` / ``cpg_removed`` metadata keys make this side
    effect explicit so downstream consumers (and tests) can attribute
    CAI suppression to CpG elimination rather than to a bug.

    Raises
    ------
    IRError
        If the input is not a valid IR-L2, or if the optimization
        changed the translated protein (semantic corruption).
    """
    _validate_input(ir_l2)
    org = organism or ir_l2.organism
    protein_no_stop = _protein_for_optimizer(ir_l2)

    # Run the existing optimizer — it takes a protein, returns a DNA
    # sequence (no stop codon) and metrics including CAI.
    #
    # CpG side-effect: optimize_sequence() defaults to
    # cpg_mode="aggressive", which removes CG dinucleotides as part of
    # the optimisation.  This is desirable for gene-therapy applications
    # (CpGs trigger TLR9) but it also suppresses CAI: the highest-CAI
    # codon for several amino acids (notably Arg = CGN) is incompatible
    # with CpG elimination.  We pass cpg_mode EXPLICITLY here (rather
    # than relying on the optimizer's default) so the side effect is
    # visible at the call site, and we record it in the output metadata
    # so consumers can attribute CAI suppression to CpG elimination
    # rather than to a bug.  See the optimize_codons() docstring's
    # "CpG side-effect" section.
    cpg_mode_for_pass = "aggressive"
    result = optimize_sequence(
        target_protein=protein_no_stop,
        organism=org,
        cpg_mode=cpg_mode_for_pass,
        strict_mode=False,
    )
    optimized_dna = result.sequence

    # Re-attach the original stop codon (the optimizer omits it).
    optimized_dna_with_stop = _re_attach_stop(optimized_dna, ir_l2)

    cai_before = ir_cai(ir_l2)
    cai_after = float(getattr(result, "cai", 0.0))

    # Count CpGs in the input and output CDS so the side effect is
    # visible in metadata.  IR L2 uses RNA (U); "CG" is identical in
    # RNA and DNA (C and G are unchanged by T↔U transcription).
    cpg_before = ir_l2.cds.upper().count("CG")
    cpg_after = optimized_dna_with_stop.count("CG")

    new_l2 = _build_optimized_l2(
        ir_l2,
        optimized_dna_with_stop,
        pass_name="optimize_codons",
        extra_metadata={
            "optimization": "codon_optimization",
            "cai_before": cai_before,
            "cai_after": cai_after,
            "cai_delta": cai_after - cai_before,
            "gc_content": float(getattr(result, "gc_content", 0.0)),
            "cpg_mode": cpg_mode_for_pass,
            "cpg_before": cpg_before,
            "cpg_after": cpg_after,
            "cpg_removed": cpg_before - cpg_after,
        },
    )

    # Check structural invariants on the new IR-L2.
    check_l2_invariants(new_l2)
    # Check semantic invariants: protein MUST be preserved.
    _assert_protein_preserved(ir_l2, new_l2, "optimize_codons")
    return new_l2


def eliminate_cpgs(
    ir_l2: IR_L2_MatureMRNA,
    organism: Optional[str] = None,
) -> IR_L2_MatureMRNA:
    """IR optimization pass: eliminate CpG dinucleotides.

    Removes ``CG`` dinucleotides from the CDS.  This is important for
    gene therapy applications, where CpG dinucleotides can trigger
    Toll-like receptor 9 (TLR9) immune responses that degrade the
    therapeutic transcript.

    PRESERVES the translated protein (semantic correctness).

    Parameters
    ----------
    ir_l2:
        A well-formed :class:`IR_L2_MatureMRNA`.
    organism:
        Optional target organism override.

    Returns
    -------
    IR_L2_MatureMRNA
        A new IR-L2 with fewer (ideally zero) ``CG`` dinucleotides in
        the CDS.  The metadata carries:

        * ``pass``         = ``"eliminate_cpgs"``
        * ``source_level`` = ``"L2"``
        * ``optimization`` = ``"cpg_elimination"``
        * ``cpg_before``   — number of ``CG`` dinucleotides in the input CDS
        * ``cpg_after``    — number of ``CG`` dinucleotides in the output CDS
        * ``cpg_removed``  — net reduction

    Raises
    ------
    IRError
        If the input is not a valid IR-L2, or if the optimization
        changed the translated protein.
    """
    _validate_input(ir_l2)
    org = organism or ir_l2.organism
    protein_no_stop = _protein_for_optimizer(ir_l2)

    # The optimizer's default cpg_mode is already "aggressive"; we pass
    # it explicitly to make the intent of THIS IR pass self-documenting.
    result = optimize_sequence(
        target_protein=protein_no_stop,
        organism=org,
        cpg_mode="aggressive",
        strict_mode=False,
    )
    optimized_dna = result.sequence
    optimized_dna_with_stop = _re_attach_stop(optimized_dna, ir_l2)

    # Count CpGs in both the original and the optimized CDS.
    # IR L2 uses RNA (U); CpG dinucleotides are CG in both alphabets
    # (C and G are unchanged by T↔U transcription).
    old_cpg = ir_l2.cds.upper().count("CG")
    # Count in the DNA form (with stop re-attached) — same result.
    new_cpg = optimized_dna_with_stop.count("CG")

    new_l2 = _build_optimized_l2(
        ir_l2,
        optimized_dna_with_stop,
        pass_name="eliminate_cpgs",
        extra_metadata={
            "optimization": "cpg_elimination",
            "cpg_before": old_cpg,
            "cpg_after": new_cpg,
            "cpg_removed": old_cpg - new_cpg,
            "cai_after": float(getattr(result, "cai", 0.0)),
        },
    )

    check_l2_invariants(new_l2)
    _assert_protein_preserved(ir_l2, new_l2, "eliminate_cpgs")
    return new_l2


# ────────────────────────────────────────────────────────────────────
# Pipeline runner
# ────────────────────────────────────────────────────────────────────

# Registry of available IR optimization passes, keyed by the names
# users pass to :func:`run_optimization_pipeline`.  This is deliberately
# a module-level dict (not a class) so new passes can be added by simply
# assigning a new entry — no inheritance, no registration ceremony.
_PASSES = {
    "optimize_codons": optimize_codons,
    "eliminate_cpgs": eliminate_cpgs,
}

# Default pass sequence when the caller doesn't specify one.  Mirrors
# the typical gene-synthesis workflow: first maximise expression (CAI),
# then remove immunostimulatory motifs (CpGs).  The order matters —
# CpG elimination runs *after* codon optimization so it can repair any
# CpGs that the high-CAI codon choices happened to introduce.
_DEFAULT_PASSES: list[str] = ["optimize_codons", "eliminate_cpgs"]


def run_optimization_pipeline(
    ir_l2: IR_L2_MatureMRNA,
    passes: Optional[list[str]] = None,
) -> IR_L2_MatureMRNA:
    """Run a sequence of IR optimization passes on an IR-L2.

    Each pass is an IR_L2 → IR_L2 transformation that preserves the
    translated protein.  Passes are applied left-to-right; the output
    of pass *i* is the input to pass *i+1*.  This mirrors the way LLVM
    applies a pass manager's pipeline: ``-O2`` is internally a list of
    passes run in sequence.

    Parameters
    ----------
    ir_l2:
        A well-formed :class:`IR_L2_MatureMRNA`.
    passes:
        List of pass names to run, in order.  ``None`` (default) uses
        :data:`_DEFAULT_PASSES` = ``["optimize_codons", "eliminate_cpgs"]``.

    Returns
    -------
    IR_L2_MatureMRNA
        The final IR-L2 after all passes have been applied.  The
        metadata carries the list of passes that ran, in order, under
        ``metadata["passes_applied"]``.

    Raises
    ------
    IRError
        If the input is not a valid IR-L2, if any pass name is unknown,
        or if any pass changes the translated protein.

    Examples
    --------
    >>> from biocompiler.ir.frontend import compile_from_spec
    >>> from biocompiler.ir.optimization import run_optimization_pipeline
    >>> ir_l2 = compile_from_spec("gene.yaml", target_level=IRLevel.L2)
    >>> optimized = run_optimization_pipeline(ir_l2)
    >>> optimized.metadata["passes_applied"]
    ['optimize_codons', 'eliminate_cpgs']
    """
    if passes is None:
        passes = list(_DEFAULT_PASSES)

    # Validate pass names up-front so a typo fails fast, before any
    # optimization work is done (matching the fail-fast principle).
    for name in passes:
        if name not in _PASSES:
            raise IRError(
                IRLevel.L2,
                f"unknown optimization pass: {name!r}. "
                f"Available: {sorted(_PASSES.keys())}",
            )

    _validate_input(ir_l2)
    current = ir_l2
    applied: list[str] = []
    for pass_name in passes:
        current = _PASSES[pass_name](current)
        applied.append(pass_name)

    # Stamp the pipeline-level provenance on top of the last pass's
    # metadata.  We do this on a *copy* of the final L2 to avoid
    # mutating the object returned by the last pass in-place (dataclasses
    # are not frozen, but our convention is "passes return fresh objects").
    final = IR_L2_MatureMRNA(
        sequence=current.sequence,
        five_utr=current.five_utr,
        cds=current.cds,
        three_utr=current.three_utr,
        organism=current.organism,
        gene_name=current.gene_name,
        metadata={
            **current.metadata,
            "passes_applied": list(applied),
        },
    )
    # Final belt-and-braces invariant check on the pipeline output.
    check_l2_invariants(final)
    _assert_protein_preserved(ir_l2, final, "run_optimization_pipeline")
    return final
