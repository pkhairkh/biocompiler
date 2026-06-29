"""
BioOptimizer pipeline: the high-level ``optimize_sequence()`` API.

Historical note (v0.9.0 cleanup): this module previously hosted a large
multi-strategy slow-path dispatch (greedy / hybrid / harmonize /
constraint_first / cai_first / CSP solver) plus the ``BioOptimizer``
class with its mixin-based optimization pipeline.  As of v0.9.0 all of
that has been removed in favour of the single-pass integrated optimizer
(``integrated_optimizer.integrated_optimize``).  The legacy strategies
are still accepted as keyword arguments to ``optimize_sequence`` for
backward compatibility, but they now emit a ``DeprecationWarning`` and
redirect to the integrated optimizer.

The ``BioOptimizer`` class is retained as a deprecated stub so that
``from biocompiler import BioOptimizer`` continues to work; instantiating
it emits a ``DeprecationWarning`` and its ``optimize()`` /
``_evaluate_all_predicates()`` methods raise ``NotImplementedError``.

The default fast path runs entirely through:

    ``optimize_sequence(use_integrated=True)``
      -> ``integrated_optimize()``
      -> ``certify_fast_path_result()``
      -> ``_attach_ir_pipeline()``

and is the only supported execution path.
"""

from typing import Any, Callable, List, Optional

import logging
import signal
import warnings

from ..type_system import AA_TO_CODONS
from ..organisms import resolve_organism
from ..shared.exceptions import (
    InvalidProteinError,
    OptimizationConstraintError,
)

from .constraints import _organism_to_species_key
from .utils import (
    OptimizationResult,
    DEFAULT_MAX_ITERATIONS,
    OPTIMIZER_TIMEOUT_SECONDS,
)
from .objectives import resolve_objective as _resolve_objective
from .pipeline_certification import certify_fast_path_result
from .pipeline_batch import batch_optimize  # re-exported for backward compat

logger = logging.getLogger(__name__)


# ── IR pipeline attachment (TIGHTEN-1) ───────────────────────────────
def _attach_ir_pipeline(
    result: OptimizationResult,
    target_protein: str,
    organism: str,
    kwargs: dict,
    log: logging.Logger,
) -> None:
    """Run the BioCompiler IR pipeline on ``result.sequence`` and attach
    the IR objects + codegen outputs to ``result`` in-place.

    This is the TIGHTEN-1 integration: every ``optimize_sequence()`` call
    runs its optimized DNA through the IR (L0->L1->L2->L3) so that:

    * ``result.ir_verified`` is ``True`` iff the IR-translated polypeptide
      matches the input protein (the "compiler" correctness check).
    * ``result.ir_l0`` / ``result.ir_l3`` carry the typed IR objects.
    * ``result.genbank`` / ``result.fasta`` / ``result.sbol3`` carry
      standard codegen outputs ready for downstream tools (Benchling,
      SnapGene, BLAST, NCBI).

    Best-effort: any IR failure is logged at DEBUG level and
    ``result.ir_verified`` stays ``False`` (its dataclass default).  The
    optimization itself is never broken by an IR issue.

    The optimizer's output DNA typically lacks a terminal stop codon
    (the optimizer returns exactly ``len(protein) * 3`` nucleotides).
    The IR ``splice`` pass requires an in-frame stop codon to mark the
    end of the CDS, so we append ``TAA`` when the optimized sequence
    does not already end with a stop codon.  This added stop codon is
    purely for IR/codegen purposes -- it is NOT part of
    ``result.sequence``.
    """
    from biocompiler.ir import IR_L0_GenomicDNA, IRLevel
    from biocompiler.ir.passes import compile_gene
    from biocompiler.ir.codegen import to_genbank, to_fasta, to_sbol3

    try:
        optimized_seq = result.sequence
        # Append a stop codon if needed so the splice pass can locate
        # the CDS boundary.
        _seq_for_ir = optimized_seq
        if not _seq_for_ir[-3:].upper() in {"TAA", "TAG", "TGA"}:
            _seq_for_ir = _seq_for_ir + "TAA"

        ir_l0 = IR_L0_GenomicDNA(
            sequence=_seq_for_ir,
            regions=[],
            organism=organism,
            gene_name=kwargs.get("gene_name"),
            secis_positions=kwargs.get('secis_positions', []),
        )
        ir_l3 = compile_gene(ir_l0, IRLevel.L3)

        # Verify protein preservation: the IR-translated polypeptide
        # must equal the input protein plus a trailing stop ('*').
        expected_protein = (
            target_protein + "*"
            if not target_protein.endswith("*")
            else target_protein
        )
        ir_verified = (ir_l3.sequence == expected_protein)

        # Attach IR objects + codegen outputs to the result.
        result.ir_l0 = ir_l0
        result.ir_l3 = ir_l3
        result.ir_verified = ir_verified
        result.genbank = to_genbank(ir_l0)
        result.fasta = to_fasta(ir_l3)
        result.sbol3 = to_sbol3(ir_l0)

        if not ir_verified:
            log.warning(
                "IR verification FAILED: IR-translated protein %r does not "
                "match expected %r (input protein was %r)",
                ir_l3.sequence, expected_protein, target_protein,
            )
    except Exception as _ir_exc:
        # IR verification is best-effort -- don't fail the optimization
        # if the IR has an issue.  Log the exception for diagnosis and
        # leave ir_verified=False (its dataclass default).
        log.debug(
            "IR pipeline integration failed: %s: %s",
            type(_ir_exc).__name__, _ir_exc, exc_info=True,
        )
        result.ir_verified = False


class OptimizerTimeout(Exception):
    """Raised when the optimizer exceeds its wall-clock time budget.

    Caught by ``optimize_sequence()`` which then returns a partial
    :class:`OptimizationResult` with ``timed_out=True`` instead of
    letting the process hang or OOM.  (Task TIGHTEN-2)
    """
    pass


def _timeout_handler(signum, frame):
    """SIGALRM handler — raises :class:`OptimizerTimeout`.

    Installed via :func:`signal.signal` around the body of
    :func:`optimize_sequence`.  When the alarm fires, Python raises
    this exception in the main thread, unwinding the (possibly deep)
    optimization call stack back to the ``try/except`` in
    ``optimize_sequence``.
    """
    raise OptimizerTimeout(
        "Optimizer exceeded its time budget (SIGALRM fired)"
    )


def _build_timeout_result(
    target_protein: str,
    organism: str,
    optimized_seq: Optional[str] = None,
    timeout_seconds: int = 30,
    extra_warnings: Optional[List[str]] = None,
) -> "OptimizationResult":
    """Build a partial :class:`OptimizationResult` on timeout.

    The returned result has ``timed_out=True`` and
    ``convergence_status="timeout"``.  The sequence is the best
    partial output available; if none was produced yet, a
    back-translated fallback is used so the length invariant
    (``len(sequence) == len(protein) * 3``) always holds.
    """
    warnings_list: List[str] = list(extra_warnings or [])
    warnings_list.append(
        f"Optimizer timed out after exceeding the {timeout_seconds}s "
        f"time budget; result is partial."
    )

    seq = optimized_seq
    # Simple back-translation fallback.  The previous
    # ``pipeline_backtranslate`` helper was removed in v0.9.0 along
    # with the rest of the slow-path stack; a plain ``AA_TO_CODONS``
    # first-codon translation is sufficient for the timeout fallback
    # (which is itself best-effort).
    if not seq or len(seq) != len(target_protein) * 3:
        seq = "".join(
            AA_TO_CODONS.get(aa, ["GCG"])[0] for aa in target_protein
        )

    if len(seq) != len(target_protein) * 3:
        seq = "".join(
            AA_TO_CODONS.get(aa, ["GCG"])[0] for aa in target_protein
        )

    gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
    try:
        from biocompiler.expression.translation import compute_cai
        cai_val = compute_cai(seq, organism)
    except Exception:
        cai_val = 0.0
    cai_val = max(0.0, min(1.0, float(cai_val)))

    return OptimizationResult(
        sequence=seq,
        gc_content=round(gc, 4),
        cai=cai_val,
        protein=target_protein,
        timed_out=True,
        convergence_status="timeout",
        warnings=warnings_list,
    )


def optimize_sequence(
    target_protein: str | None = None,
    organism: str = "Homo_sapiens",
    species: str | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.5,
    enzymes: list[str] | None = None,
    strategy: str = "hybrid",
    source_organism: str | None = None,
    harmonization_cai_weight: float = 0.5,
    use_csp_solver: bool = False,
    optimize_mrna_stability: bool = True,
    strict_mode: bool = False,
    seed: int | None = None,
    include_utr: bool = True,
    consider_codon_pair_bias: bool = False,
    track_provenance: bool = True,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
    objective: str | Callable[..., float] | None = None,
    gc_window_size: int = 50,
    gc_window_min: float | None = None,
    gc_window_max: float | None = None,
    biosecurity_mode: str | None = None,
    skip_biosecurity_check: bool = False,
    cpg_mode: str = "aggressive",
    timeout_seconds: int = OPTIMIZER_TIMEOUT_SECONDS,
    use_integrated: bool = True,
    use_context_aware_gt: bool = False,
    **kwargs: Any,
) -> OptimizationResult:
    """Optimize a protein sequence for expression in the target organism.

    This is the high-level convenience API.  It takes a protein sequence
    and returns an :class:`OptimizationResult` with the optimized DNA
    sequence and quality metrics.

    As of v0.9.0 the integrated optimizer is the only supported
    execution path.  The legacy ``use_integrated=False``,
    ``use_csp_solver=True``, and non-default ``strategy`` values
    (``"harmonize"``, ``"constraint_first"``, ``"cai_first"``) are
    accepted for backward compatibility but emit a
    :class:`DeprecationWarning` and redirect to the integrated
    optimizer.  They will become an error in v1.0.0.

    Organism Specification:

        The target organism can be specified using **either** the
        ``organism`` parameter **or** the ``species`` parameter.  Both
        accept the same set of names — short aliases, abbreviated
        binomials, display names, or full canonical names — and both
        map to the same internal representation via
        :func:`~biocompiler.organisms.resolve_organism`.

        +-------------------------+------------------------+
        | Example                 | Resolved to            |
        +=========================+========================+
        | ``species='ecoli'``     | ``'Escherichia_coli'`` |
        | ``organism='ecoli'``    | ``'Escherichia_coli'`` |
        | ``species='E. coli'``   | ``'Escherichia_coli'`` |
        | ``organism='Escherichia_coli'`` | same           |
        | ``species='human'``     | ``'Homo_sapiens'``     |
        | ``species='h_sapiens'`` | ``'Homo_sapiens'``     |
        +-------------------------+------------------------+

        If both ``species`` and ``organism`` are provided, ``species``
        takes precedence and a :class:`DeprecationWarning` is emitted
        recommending the use of ``organism`` only (the canonical
        parameter).  If neither matches the other after resolution, a
        warning is logged.

    Args:
        target_protein: Amino acid sequence (1-letter codes, no stop).
            Can also be passed as the first positional argument.
        organism: Target organism.  Accepts canonical binomials
            (e.g., ``'Homo_sapiens'``, ``'Escherichia_coli'``),
            short keys (``'ecoli'``, ``'human'``), abbreviated
            binomials (``'E_coli'``, ``'h_sapiens'``), or display
            names (``'E. coli'``).  All forms are resolved via
            :func:`~biocompiler.organisms.resolve_organism`.
        species: Alias for ``organism``.  Accepts the same values.
            If provided **together with** ``organism``, ``species``
            takes precedence and a deprecation warning is emitted.
            Prefer using ``organism`` in new code; ``species`` is
            retained for backward compatibility.
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        gc_window_size: Sliding window size (in nucleotides) for local
            GC content checking.  Set to 0 to disable sliding-window GC
            constraints.  Default: 50.
        gc_window_min: Minimum acceptable GC fraction for each sliding
            window.  If None (default), uses ``gc_lo``.
        gc_window_max: Maximum acceptable GC fraction for each sliding
            window.  If None (default), uses ``gc_hi``.
        cai_threshold: Minimum CAI score for the CodonAdapted predicate.
        enzymes: List of restriction enzyme names to avoid.
        strategy: Optimization strategy.  **DEPRECATED in v0.9.0** —
            the only supported strategy is the integrated optimizer
            (formerly ``"hybrid"``).  The values ``"constraint_first"``,
            ``"cai_first"``, and ``"harmonize"`` are accepted for
            backward compatibility but emit a ``DeprecationWarning`` and
            fall through to the integrated optimizer.  Will become an
            error in v1.0.0.
        source_organism: Source organism for codon harmonization.
            **DEPRECATED** — only meaningful for the removed
            ``strategy='harmonize'`` path; now ignored.
        harmonization_cai_weight: Blend weight for CAI fallback in
            harmonization.  **DEPRECATED** — only meaningful for the
            removed ``strategy='harmonize'`` path; now ignored.
        use_csp_solver: If True, try CSP/SMT solver before greedy.
            **DEPRECATED in v0.9.0** — accepted for backward
            compatibility but emits a ``DeprecationWarning`` and
            redirects to the integrated optimizer.
        optimize_mrna_stability: If True, run a soft mRNA stability
            improvement pass after CAI optimization.  Passed through to
            the integrated optimizer.
        seed: Deterministic seed for reproducibility.  Currently unused
            as the integrated optimizer is fully deterministic, but
            reserved for future randomized strategies.
        include_utr: If True, generate organism-appropriate 5' and 3' UTR
            suggestions and include them in the result.
        consider_codon_pair_bias: If True, run a soft codon pair bias
            (CPB) improvement pass after the CAI optimization pass.
        track_provenance: If True, collect decision-level provenance.
        strict_mode: If True (default), raise
            :class:`OptimizationConstraintError` when any predicates
            fail instead of returning a result with
            ``failed_predicates``.
        max_iterations: Maximum number of optimization iterations
            (advisory; the integrated optimizer is single-pass).
        objective: Custom optimization objective.
        **kwargs: Additional arguments (e.g., splice_low, splice_high,
            restriction_sites, organism_domain).
        biosecurity_mode: Biosecurity screening mode. One of ``"enforce"``
            (default — raises :class:`BiosecurityError` on high/critical
            risk), ``"warn"`` (logs warning but proceeds), or ``"off"``
            (skips screening entirely, for testing only).  If ``None``
            (default), the ``BIOCOMPILER_BIOSECURITY_MODE`` environment
            variable is consulted, falling back to ``"enforce"``.
        use_integrated: If True (default), use the single-pass integrated
            optimizer.  If False, **DEPRECATED** — emit a
            ``DeprecationWarning`` and use the integrated optimizer
            anyway.

    Returns:
        OptimizationResult with optimized sequence and metrics.

    Raises:
        InvalidProteinError: if the protein contains invalid amino acid codes.
        UnsupportedOrganismError: if the organism is not supported.
        BiosecurityError: if ``biosecurity_mode="enforce"`` and the
            sequence matches known hazardous signatures (high/critical risk).
        OptimizationConstraintError: if ``strict_mode=True`` and the
            optimized sequence has one or more failed predicates.
    """
    # ── Resolve custom objective ──────────────────────────────────────
    # Validated early so an unknown objective name fails fast; the
    # resolved callable is unused on the integrated fast path.
    _resolve_objective(objective)

    # Set deterministic seed if provided (reserved for future randomized steps).
    # Use a per-call random.Random instance instead of random.seed() to avoid
    # polluting global random state (thread-safety / reproducibility concern).
    if seed is not None:
        import random as _random_mod
        _random_mod.Random(seed)

    # Start timing for provenance
    import time as _time
    _start_time = _time.monotonic()

    # Handle positional arg: optimize_sequence("MVHLTPEEK", organism="...")
    if target_protein is None and len(kwargs.get("_args", [])) > 0:
        target_protein = kwargs["_args"][0]

    # ── Organism resolution ────────────────────────────────────────
    # Support both 'species' and 'organism' parameters.  Both accept
    # the same set of aliases and resolve to the same canonical name.
    # 'species' is retained for backward compatibility but is not the
    # preferred parameter — use 'organism' in new code.
    _resolved_organism: str | None = None

    # Legacy: species passed via kwargs (old calling convention)
    _species_from_kwargs = kwargs.pop("species", None)
    if _species_from_kwargs is not None and species is None:
        species = _species_from_kwargs

    if species is not None:
        _resolved_organism = resolve_organism(species, strict=False)
        # Emit deprecation warning if both species and organism are given
        if organism != "Homo_sapiens":
            _resolved_organism_from_explicit = resolve_organism(organism, strict=False)
            if _resolved_organism != _resolved_organism_from_explicit:
                warnings.warn(
                    f"Both 'species={species!r}' and 'organism={organism!r}' "
                    f"were provided but resolve to different organisms "
                    f"({_resolved_organism!r} vs {_resolved_organism_from_explicit!r}). "
                    f"Using 'species' ({_resolved_organism!r}). "
                    f"Prefer using only 'organism' in new code.",
                    DeprecationWarning,
                    stacklevel=2,
                )
            else:
                warnings.warn(
                    f"Both 'species' and 'organism' were provided. "
                    f"Prefer using only 'organism' in new code; "
                    f"'species' is retained for backward compatibility.",
                    DeprecationWarning,
                    stacklevel=2,
                )
        else:
            # species was given but organism was left at default —
            # no conflict, just use species. Still emit a gentle hint.
            warnings.warn(
                f"The 'species' parameter is deprecated in favor of 'organism'. "
                f"Use organism='{_resolved_organism}' instead of "
                f"species='{species}'. Both accept the same aliases.",
                DeprecationWarning,
                stacklevel=2,
            )

    if _resolved_organism is not None:
        organism = _resolved_organism
    else:
        organism = resolve_organism(organism, strict=False)

    # ── Protein normalization & validation ──────────────────────────────
    # Normalize (strip whitespace, uppercase, drop trailing stop-codon
    # marker '*') and validate the amino-acid alphabet BEFORE entering
    # the integrated fast path.
    #
    # This MUST run before the fast path: previously the validation lived
    # *after* the `if use_integrated and target_protein:` block, so an
    # invalid input like 'MXYZ!!!@#' took the fast path, was silently
    # translated by the integrated optimizer (X/Z/!/etc -> Alanine), and
    # the raw invalid string was stored verbatim as result.protein — a
    # silent data-corruption bug. (Task W2-a.)
    if target_protein is None:
        target_protein = ""
    target_protein = target_protein.strip().upper().rstrip("*")

    if not target_protein:
        raise InvalidProteinError("", {"<empty>"})

    # Allow the 20 standard AAs plus U (selenocysteine) and * (stop).
    # Trailing '*' is stripped above; '*' is kept in the valid set as a
    # safety net so an internal stop marker doesn't surface as a confusing
    # "invalid amino acid" error from a downstream component.
    valid_aas = set("ACDEFGHIKLMNPQRSTVWYU*")
    invalid = set(target_protein) - valid_aas
    if invalid:
        raise InvalidProteinError(target_protein, invalid)

    # ── Context-aware GT avoidance (Task W2-A1) ────────────────────────
    # If the caller requested context-aware GT avoidance, delegate to the
    # two-pass optimizer instead of running the legacy single-pass
    # integrated optimizer.  The two-pass optimizer:
    #   (1) calls this function (optimize_sequence) for a pass-1 baseline,
    #   (2) scans the result for high-scoring cryptic splice donors
    #       (MaxEntScan score > threshold), and
    #   (3) repairs only those positions by swapping to the next-best
    #       codon that doesn't create a new cryptic site.
    # This recovers CAI lost to the legacy global-GT-avoidance logic
    # (which discards otherwise-optimal valine GTN codons) while still
    # preventing functional cryptic splice sites.
    #
    # This is an ADD-ONLY change: when ``use_context_aware_gt=False``
    # (the default), the existing code path below runs unchanged.
    if use_context_aware_gt:
        from .two_pass import optimize_two_pass
        # Pass through all user-supplied parameters so the pass-1
        # baseline optimization inside optimize_two_pass respects
        # gc_lo / gc_hi / strict_mode / enzymes / etc.  The two-pass
        # optimizer accepts **kwargs and forwards them to its internal
        # optimize_sequence() call.  Biosecurity screening, organism
        # resolution, and protein validation have already run above; the
        # two-pass optimizer's internal optimize_sequence() call will
        # re-run biosecurity screening (idempotent) for defense in depth.
        return optimize_two_pass(
            target_protein,
            organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cai_threshold=cai_threshold,
            enzymes=enzymes,
            strict_mode=strict_mode,
            seed=seed,
            track_provenance=track_provenance,
            include_utr=include_utr,
            max_iterations=max_iterations,
            objective=objective,
            gc_window_size=gc_window_size,
            gc_window_min=gc_window_min,
            gc_window_max=gc_window_max,
            biosecurity_mode=biosecurity_mode,
            skip_biosecurity_check=skip_biosecurity_check,
            cpg_mode=cpg_mode,
            timeout_seconds=timeout_seconds,
            **kwargs,
        )

    # ── Biosecurity screening (BEFORE any optimization work) ────────────
    # Screen the input protein for known hazardous signatures BEFORE any
    # optimization begins (README invariant: "input sequences screened
    # against pathogen/toxin signature databases before any optimization
    # work"). In "enforce" mode (default), a high/critical finding raises
    # BiosecurityError and aborts before the optimizer is invoked. (W2-b.)
    from ..biosecurity import check_biosecurity_before_optimize
    _biosecurity_result = check_biosecurity_before_optimize(
        target_protein,
        organism=organism,
        biosecurity_mode=biosecurity_mode,
        skip_biosecurity_check=skip_biosecurity_check,
    )

    # ── Slow path / non-default strategies (DEPRECATED) ──────────────────
    # As of v0.9.0, the integrated optimizer is the only supported path.
    # The legacy slow-path strategies (hybrid, harmonize, constraint_first,
    # cai_first) and CSP solver have been removed. use_integrated=False and
    # use_csp_solver=True now redirect to the integrated optimizer with a
    # DeprecationWarning.  The `strategy="hybrid"` default and
    # `strategy=None` remain compatible with the integrated path — the
    # integrated optimizer implements the hybrid strategy in spirit.
    if (not use_integrated) or use_csp_solver or (strategy not in ("hybrid", None)):
        warnings.warn(
            f"The legacy slow path (use_integrated={use_integrated}, "
            f"strategy={strategy!r}, use_csp_solver={use_csp_solver}) is deprecated "
            f"and has been removed. Redirecting to the integrated optimizer. "
            f"This will become an error in v1.0.0.",
            DeprecationWarning,
            stacklevel=2,
        )

    # ── Optimizer timeout (Task TIGHTEN-2) ────────────────────────
    # Wrap the integrated fast path in a SIGALRM-based timeout so the
    # optimizer NEVER hangs — it either completes or returns a partial
    # result with timed_out=True.  signal.alarm is the safety net; the
    # integrated optimizer itself is single-pass and fast.
    _timeout_old_handler = None
    _timeout_alarm_set = False
    if timeout_seconds and timeout_seconds > 0 and hasattr(signal, 'SIGALRM'):
        try:
            _timeout_old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
            signal.alarm(int(timeout_seconds))
            _timeout_alarm_set = True
        except (ValueError, OSError):
            # signal.alarm only works in the main thread on Unix;
            # in worker threads / Windows we fall back to no-timeout mode.
            logger.debug(
                "Could not install SIGALRM timeout (not main thread / "
                "unsupported platform); integrated optimizer will run "
                "without an outer wall-clock cap."
            )

    try:
        # ── Integrated optimizer fast path (only supported path) ─────────
        # Single-pass constraint-solving optimizer.  10-17x faster than
        # the removed sequential slow path.  Falls back to a partial
        # timeout result if SIGALRM fires; all other exceptions
        # propagate to the caller (no slow-path fallback anymore).
        from .integrated_optimizer import integrated_optimize
        _int_dna, _int_notes, _secis_pos = integrated_optimize(
            protein=target_protein,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            # Default the enzyme list the SAME way every other call site
            # in this file does (`enzymes or ["EcoRI", "BamHI", "XhoI",
            # "HindIII", "NotI"]`).  Without this default, the integrated
            # optimizer would not avoid restriction sites when the caller
            # omits `enzymes`, while certify_fast_path_result WOULD check
            # against the default enzyme list -- producing a
            # NoRestrictionSite failure (e.g. EcoRI at pos 363 on HBB)
            # that the optimizer never had a chance to avoid. (Task W5-b.)
            enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            cpg_mode=cpg_mode,
            # NOTE: strict_mode is enforced downstream by
            # certify_fast_path_result(strict_mode=strict_mode, ...) which
            # raises OptimizationConstraintError on predicate failure.
            # The integrated_optimize() function itself no longer accepts
            # a (previously dead) strict_mode parameter. (Task W4 / H2.)
        )

        # NOTE: Biosecurity screening now runs BEFORE integrated_optimize()
        # — see the "Biosecurity screening" block above the fast path.
        # (Task W2-b: previously called here, AFTER optimization, which
        # violated the README "screened before any optimization work"
        # invariant.)

        # Build result.
        # NOTE: OptimizationResult is imported at module top (see the
        # `from .utils import ...` block at the top of this file). Do NOT
        # re-import it locally here: a local import would make Python treat
        # `OptimizationResult` as a function-local name for the ENTIRE
        # optimize_sequence() scope, breaking the timeout fallback path in
        # _build_timeout_result. (Task W1-a.)
        from biocompiler.sequence.scanner import gc_content as _gc_fn
        from biocompiler.expression.translation import compute_cai as _cai_fn
        _int_cai = _cai_fn(_int_dna, organism) if _int_dna else 0.0
        _int_gc = _gc_fn(_int_dna) if _int_dna else 0.0

        # Strip stop codon — OptimizationResult expects len(seq) == len(protein)*3
        _int_dna_no_stop = _int_dna[:-3] if _int_dna.endswith(("TAA", "TAG", "TGA")) else _int_dna
        _int_result = OptimizationResult(
            sequence=_int_dna_no_stop,
            cai=_int_cai,
            gc_content=_int_gc,
            protein=target_protein,
            biosecurity_screening_result=_biosecurity_result,
        )
        _int_result.timed_out = False
        # The integrated optimizer always converges in a single forward
        # pass + bounded cleanup (deterministic, no iteration cap), so we
        # report "converged" rather than a bespoke "integrated" status.
        # This aligns the fast path with the ConvergenceTracker contract
        # expected by tests/test_convergence.py (Task W4 / H10).
        _int_result.convergence_status = "converged"
        _int_result.iterations_used = 1  # single-pass integrated optimizer

        # ── Compute organism-specific constraint flags ──
        # so the certification helper applies the right
        # prokaryote/eukaryote-aware predicate set.
        _int_organism_domain = kwargs.get("organism_domain", "auto")
        if _int_organism_domain not in ("auto", "eukaryote", "prokaryote"):
            _int_organism_domain = "auto"
        if _int_organism_domain == "auto":
            from biocompiler.organisms.config import is_eukaryotic_organism
            _int_is_prok = not is_eukaryotic_organism(organism)
        elif _int_organism_domain == "prokaryote":
            _int_is_prok = True
        else:
            _int_is_prok = False
        _int_avoid_gt = kwargs.get("avoid_gt", not _int_is_prok)
        _int_splice_low = kwargs.get("splice_low", 3.0 if not _int_is_prok else 999.0)
        _int_species_key = _organism_to_species_key(organism)

        # ── Post-optimization certification (post-optimization certification) ──
        # Evaluate all 20 predicates on the optimized DNA, generate the
        # certificate, populate provenance / decision_trail / UTR fields,
        # and enforce strict_mode.  MUST run BEFORE _attach_ir_pipeline
        # so the strict_mode raise (if any) skips the expensive IR pass.
        certify_fast_path_result(
            result=_int_result,
            target_protein=target_protein,
            organism=organism,
            species_key=_int_species_key,
            gc_lo=gc_lo, gc_hi=gc_hi,
            cai_threshold=cai_threshold,
            enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            effective_splice_low=_int_splice_low,
            effective_avoid_gt=_int_avoid_gt,
            is_prokaryote=_int_is_prok,
            gc_window_size=gc_window_size,
            gc_window_min=gc_window_min,
            gc_window_max=gc_window_max,
            cpg_mode=cpg_mode,
            strict_mode=strict_mode,
            track_provenance=track_provenance,
            include_utr=include_utr,
            seed=seed,
            start_time=_start_time,
            tissue=kwargs.get("tissue", ""),
            enable_mutagenesis=kwargs.get("enable_mutagenesis", False),
            logger=logger,
        )
        # If strict_mode=True and predicates failed,
        # certify_fast_path_result raises OptimizationConstraintError
        # (with partial_result=_int_result, certificate already
        # populated) — caught by the dedicated except clause below so it
        # propagates to the caller.

        # Surface the integrated optimizer's cleanup notes (e.g., "GT
        # count = 7 left after cleanup") to the caller as warnings. (W2-c
        # spec §3.5 rule #3.)
        if _int_notes:
            _int_result.warnings = list(_int_result.warnings) + list(_int_notes)

        # Attach IR pipeline + codegen.
        # Pass SECIS positions through kwargs so the IR pipeline can use them.
        _int_kwargs = dict(kwargs)
        _int_kwargs['secis_positions'] = _secis_pos
        _attach_ir_pipeline(_int_result, target_protein, organism, _int_kwargs, logger)
        logger.info("Integrated optimizer: %dbp CAI=%.3f GC=%.3f IR=%s",
                    len(_int_dna), _int_cai, _int_gc, _int_result.ir_verified)
        return _int_result

    except OptimizationConstraintError:
        # strict_mode contract — propagate the
        # raise instead of swallowing it.  The partial_result on the
        # exception already carries a populated certificate +
        # predicate_results (per W2-c spec §4.1).
        raise
    except OptimizerTimeout:
        # The optimizer exceeded its time budget.  Return a partial
        # result with timed_out=True rather than letting the process
        # hang.
        logger.warning(
            "Optimizer exceeded %ds time budget for protein of length "
            "%d; returning partial result.",
            timeout_seconds, len(target_protein),
        )
        return _build_timeout_result(
            target_protein=target_protein,
            organism=organism,
            timeout_seconds=timeout_seconds,
        )
    finally:
        if _timeout_alarm_set:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, _timeout_old_handler)


# batch_optimize is in pipeline_batch.py — re-exported at top of file.


class BioOptimizer:
    """DEPRECATED stub retained for backward import compatibility.

    As of v0.9.0 the multi-strategy BioOptimizer pipeline (greedy /
    hybrid / constraint_first / cai_first / CSP solver) has been removed
    in favour of the single-pass integrated optimizer.  This class is
    kept only so that ``from biocompiler import BioOptimizer`` continues
    to succeed; instantiating it emits a :class:`DeprecationWarning` and
    its ``optimize()`` / ``_evaluate_all_predicates()`` methods raise
    :class:`NotImplementedError`.  New code should call
    :func:`optimize_sequence` directly.

    The strategy-mixin base classes (``CAIFirstStrategyMixin``,
    ``ConstraintFirstStrategyMixin``, ``CrossCodonMixin``) and the
    ``_step_*`` / ``_optimize_*`` method bodies that previously
    composed the slow-path pipeline have all been removed along with
    the slow-path stack (greedy, hybrid_*, strategy_*, pipeline_cross_codon,
    incremental, performance, pipeline_cai/cpg/gc/rs/postprocessing,
    pipeline_backtranslate, pipeline_validation, blast_avoidance, etc.).
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        warnings.warn(
            "BioOptimizer is deprecated as of v0.9.0 and its multi-strategy "
            "pipeline has been removed. Use biocompiler.optimize_sequence() "
            "instead. This will become an error in v1.0.0.",
            DeprecationWarning,
            stacklevel=2,
        )
        # Store common attributes for legacy code that reads them after
        # construction (e.g. cli_services.py, pipeline_batch.py).  These
        # are best-effort: callers should migrate to optimize_sequence().
        self.species: str = kwargs.get("species", "ecoli")
        self.enzymes: List[str] = list(kwargs.get("enzymes") or [])
        self.organism_name: str = kwargs.get("organism_name") or self.species
        self.strategy: str = kwargs.get("strategy", "constraint_first")
        self._applied_mutagenesis: List[dict] = []
        self._original_protein: str = ""

    def optimize(self, *args: Any, **kwargs: Any):
        raise NotImplementedError(
            "BioOptimizer.optimize() has been removed in v0.9.0 along with "
            "the multi-strategy slow path. Use "
            "biocompiler.optimize_sequence() instead."
        )

    def _evaluate_all_predicates(self, *args: Any, **kwargs: Any):
        raise NotImplementedError(
            "BioOptimizer._evaluate_all_predicates() has been removed in "
            "v0.9.0 along with the CrossCodonMixin. Use "
            "biocompiler.optimize_sequence() instead."
        )
