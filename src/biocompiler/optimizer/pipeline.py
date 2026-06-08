"""
BioOptimizer pipeline: the main optimizer class and high-level API.

Contains the BioOptimizer class with its multi-step optimization pipeline,
the optimize_sequence() and batch_optimize() top-level functions, and
back-translation helpers.
"""

from typing import Callable, List, Dict, Optional, Tuple, Set, Any

import logging
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from itertools import product as itertools_product

from ..type_system import (
    CODON_TABLE, AA_TO_CODONS, BLOSUM62, PredicateResult,
    check_no_stop_codons, check_no_cryptic_splice, check_no_cpg_island,
    check_no_restriction_site, check_no_avoidable_gt,
    check_no_gt_dinucleotide_soft, _compute_max_gt_count,
    check_valid_coding_seq,
    find_cross_codon_gt, find_cross_codon_cg, find_cross_codon_restriction,
)
from ..organisms import CODON_ADAPTIVENESS_TABLES, ORGANISM_GC_TARGETS, resolve_organism, ORGANISM_ALIASES, SPECIES_SHORT_NAMES, SUPPORTED_ORGANISMS
from ..constants import reverse_complement, RESTRICTION_ENZYMES, IUPAC_EXPAND, VALID_IUPAC_BASES

# Optimization thresholds and sentinel values
# (Originally defined inline in optimization.py; imported from .greedy)
from .greedy import (
    IUPAC_EXPANSION_CAP,
    ELIMINATED_SITE_SCORE,
    TOP_CAI_ALTERNATIVES,
    T_RUN_LENGTH_THRESHOLD,
    SPLICE_DONOR_POTENTIAL_THRESHOLD,
    _MAX_ACCEPTOR_SEARCH_DIST,
    EUKARYOTE_CAI_GT_COST_THRESHOLD,
    GT_BOUNDARY_CAI_TOLERANCE,
    GT_CAI_LOG_ADAPTIVENESS_COST,
)
from ..scanner import gc_content
from ..sliding_gc import check_sliding_gc, fix_sliding_gc_violations
from ..maxentscan import score_donor, score_acceptor, max_donor_score, max_acceptor_score
from ..mutagenesis import propose_mutagenesis, MutagenesisReport, MutagenesisProposal
from ..incremental import IncrementalSequenceState, CodonCache, EnzymeSiteCache
from ..certificate import format_certificate
from ..exceptions import InvalidProteinError, UnsupportedOrganismError, OptimizationConstraintError, BiosecurityError
from ..decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    DecisionProvenanceCollector,
    OptimizationDecisionTrail,
)
from ..aho_corasick import AhoCorasickScanner, build_scanner_from_enzymes, build_scanner_from_sites  # type: ignore[import-untyped]
from ..objectives import resolve_objective as _resolve_objective, cai_objective as _cai_objective

# Internal optimizer sub-module imports
from .cai import HAS_NUMBA, _HAS_NUMBA, _compute_cai_fast, _count_dinucs_fast, _BatchSwapScorer
from .constraints import (
    _find_site_in_sequence, _get_overlapping_codons, _remove_site_multicodon,
    _find_gt_free_codons, _find_ag_free_codons,
    _organism_to_species_key, _species_key_to_organism,
    _back_translate_protein, _back_translate_protein_dp,
    _build_restriction_site_set, _contains_restriction_site,
    _count_dinucleotides,
    _count_gts, _is_unavoidable_gt, _has_gt, _codon_creates_boundary_gt,
)
from .utils import (
    protein_to_aa_list, ConvergenceTracker, OptimizationResult, FullConstructResult,
    _OptConfig,
    MAX_RESTRICTION_SITE_ITERATIONS, MAX_IUPAC_SITE_ITERATIONS,
    MAX_ATTTA_MOTIF_ITERATIONS, MAX_T_RUN_ITERATIONS,
    MAX_GC_ADJUSTMENT_ITERATIONS, MAX_SPLICE_ELIMINATION_ITERATIONS,
    MAX_CPG_DISRUPTION_ITERATIONS,
    DEFAULT_MAX_ITERATIONS,
    CONVERGENCE_IMPROVEMENT_THRESHOLD,
    CONVERGENCE_PATIENCE,
    OSCILLATION_WINDOW,
)
# Sub-modules extracted from pipeline.py for maintainability
from . import postprocessing as _postprocessing
from . import validation as _validation
from .greedy import (
    score_splice_donor_potential,
    _gt_aware_select_codon,
    _is_in_codon_gt,
    _eukaryote_cai_recovery,
    _eliminate_cpg_dinucleotides,
    _greedy_optimize,
    _expand_iupac_site,
    _check_predicates_via_type_system,
)

logger = logging.getLogger(__name__)

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
    strict_mode: bool = True,
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
    check_g4: bool = True,
    check_mrna_halflife: bool = True,
    halflife_threshold: float = 2.0,
    use_lineardesign: bool = False,
    lineardesign_lambda: float = 3.0,
    tasep_ensemble: bool = True,
    tasep_n_runs: int = 100,
    check_nmd: bool = True,
    **kwargs: Any,
) -> OptimizationResult:
    """Optimize a protein sequence for expression in the target organism.

    This is the high-level convenience API that wraps BioOptimizer.
    It takes a protein sequence and returns an OptimizationResult with
    the optimized DNA sequence and quality metrics.

    When ``use_csp_solver=True``, the CSP/SMT constraint solver is tried
    first (OR-Tools CP-SAT or Z3 SMT).  If no solver is available or the
    solver times out, the greedy optimizer is used as fallback.

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
        strategy: Optimization strategy ('hybrid', 'constraint_first', 'cai_first',
            or 'harmonize').
            'hybrid' (default for v10): Greedy init + priority-queue local search + hill climb.
            'constraint_first': Legacy sequential pipeline.
            'cai_first': Maximize CAI first, then fix constraints.
            'harmonize': Use codon harmonization (Claassens RCA method) as initialization
            instead of greedy CAI maximization. Preserves translational kinetics by
            matching the source organism's relative codon usage profile in the target
            host. Requires ``source_organism`` parameter.
        source_organism: Source organism for codon harmonization (used with
            ``strategy='harmonize'``). When specified, the harmonized sequence
            matches the source organism's Relative Codon Adaptation (RCA) profile
            in the target host. If ``None`` and ``strategy='harmonize'``, defaults
            to the target organism (self-harmonization).
        harmonization_cai_weight: Blend weight for CAI fallback in harmonization
            (0.0–1.0). 0.0 = pure harmonization, 1.0 = pure CAI maximization.
            Only used with ``strategy='harmonize'``. Default 0.5.
        use_csp_solver: If True, try CSP/SMT solver before greedy (default False).
        optimize_mrna_stability: If True, run a soft mRNA stability improvement
            pass after CAI optimization. Identifies destabilizing motifs (ATTTA,
            extended AREs, long A/T runs) and suggests synonymous codon changes
            to remove them without breaking hard constraints. Defaults to True.
        seed: Deterministic seed for reproducibility. Currently unused as the
            greedy optimizer is fully deterministic, but reserved for future
            randomized optimization strategies (e.g., Monte Carlo, genetic
            algorithms). If provided, ``random.seed(seed)`` is called at the
            start of optimization to ensure reproducible behavior when
            randomized steps are added.
        include_utr: If True, generate organism-appropriate 5' and 3' UTR
            suggestions and include them in the result. The UTRs are SUGGESTED
            but not enforced — they're for the user to consider when ordering
            a synthesis construct. Defaults to True.
        consider_codon_pair_bias: If True, run a soft codon pair bias (CPB)
            improvement pass after the CAI optimization pass.  For each
            adjacent codon pair, check if there's a synonymous pair with
            higher combined CAI+CPB score (70% CAI / 30% CPB by default)
            that doesn't violate hard constraints (restriction sites, GC
            range, splice sites).  The result will include a
            ``codon_pair_bias`` field with the mean CPB score.  Defaults to
            False.
        track_provenance: If True, collect decision-level provenance for
            every codon choice made during optimization.  The result will
            include a ``decision_trail`` field with the full audit trail.
            If False, skip provenance collection entirely (zero overhead).
            Defaults to True.
        strict_mode: If True (default), raise
            :class:`OptimizationConstraintError` when any predicates fail
            instead of returning a result with ``failed_predicates``.  This
            provides a hard-stop guarantee that the returned sequence
            satisfies all constraints.  Set to False to allow partial
            results with failed predicates.
        max_iterations: Maximum number of optimization iterations before
            declaring ``convergence_status='max_iterations'``.  A sensible
            default of 1000 is used.  The optimizer will stop early if
            convergence is detected (objective unchanged for 3 consecutive
            iterations) or oscillation is detected (objective cycles up and
            down).  The result will include ``convergence_status`` and
            ``iterations_used`` fields.
        objective: Custom optimization objective.  Controls how codon
            alternatives are scored during the CAI hill-climb phase.
            Accepts:

            - ``None`` (default): maximize CAI only (equivalent to
              ``"cai"``).
            - A string name: one of ``"cai"``, ``"cai_gc_balanced"``,
              ``"codon_pair"``, ``"min_max_gc"``.
            - A callable ``(dna: str, protein: str, organism: str) -> float``:
              a custom function where higher values are better.

            When a non-CAI objective is provided, the CAI hill-climb
            phase scores candidate codon substitutions by evaluating the
            custom objective on the *modified* sequence rather than by
            CAI alone.  This allows users to optimize for GC balance,
            codon pair bias, or any custom metric while still respecting
            all hard constraints (restriction sites, GC range, splice
            sites, etc.).
        **kwargs: Additional arguments (e.g., splice_low, splice_high,
            restriction_sites, organism_domain).
        biosecurity_mode: Biosecurity screening mode. One of ``"enforce"``
            (default — raises :class:`BiosecurityError` on high/critical
            risk), ``"warn"`` (logs warning but proceeds), or ``"off"``
            (skips screening entirely, for testing only).  If ``None``
            (default), the ``BIOCOMPILER_BIOSECURITY_MODE`` environment
            variable is consulted, falling back to ``"enforce"``.
        check_g4: If True (default), check for G-quadruplex motifs in the
            optimized sequence and fix them via synonymous codon substitution.
            G-quadruplexes can stall polymerases and cause genomic instability.
            Gracefully skipped if the ``g_quadruplex`` submodule is unavailable.
        check_mrna_halflife: If True (default), use advanced mRNA half-life
            prediction to optimize mRNA stability instead of the simpler
            ATTTA/T-run approach.  Falls back to the legacy method if the
            ``mrna_halflife`` submodule is unavailable.
        halflife_threshold: Minimum acceptable mRNA half-life in hours.
            Sequences with predicted half-life below this threshold will be
            optimized.  Default: 2.0 hours.
        use_lineardesign: If True, use LinearDesign for MFE/CAI joint
            optimization (requires C++ binary).  Default: False.
        lineardesign_lambda: MFE/CAI tradeoff parameter for LinearDesign.
            Higher values prioritize CAI over MFE.  Default: 3.0.
        tasep_ensemble: If True (default), use ensemble TASEP averaging
            for ribosome simulation (slower but more accurate).  Falls
            back to single-run TASEP if the ``ribosome_simulation``
            submodule is unavailable.
        tasep_n_runs: Number of TASEP runs for ensemble averaging.
            Default: 100.
        check_nmd: If True (default), check for nonsense-mediated decay
            (NMD) triggers in eukaryotic sequences.  NMD triggers can
            drastically reduce mRNA levels.  Gracefully skipped for
            prokaryotes or if the ``rna_degradation`` submodule is
            unavailable.

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
    _objective_fn = _resolve_objective(objective)

    # Set deterministic seed if provided (reserved for future randomized steps).
    # Use a per-call random.Random instance instead of random.seed() to avoid
    # polluting global random state (thread-safety / reproducibility concern).
    _rng = None
    if seed is not None:
        import random as _random_mod
        _rng = _random_mod.Random(seed)

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
                import warnings
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
                import warnings
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
            import warnings
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

    if not target_protein:
        raise InvalidProteinError("", {"<empty>"})

    target_protein = target_protein.strip().upper()

    # Validate protein
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = set(target_protein) - valid_aas
    if invalid:
        raise InvalidProteinError(target_protein, invalid)

    # ── Biosecurity screening ───────────────────────────────────────
    # Screen the input protein for known hazardous signatures BEFORE
    # any optimization begins.  In "enforce" mode (the default), a
    # high/critical risk finding raises BiosecurityError and aborts.
    # In "warn" mode, the finding is logged but optimization proceeds.
    # In "off" mode, screening is skipped entirely.
    from ..biosecurity import check_biosecurity_before_optimize
    _biosecurity_result = None
    try:
        _biosecurity_result = check_biosecurity_before_optimize(
            target_protein,
            organism=organism,
            biosecurity_mode=biosecurity_mode,
            skip_biosecurity_check=skip_biosecurity_check,
        )
    except BiosecurityError:
        raise

    # Initialize decision provenance collector if tracking is enabled
    _provenance_collector: DecisionProvenanceCollector | None = None
    if track_provenance:
        _provenance_collector = DecisionProvenanceCollector()
        _provenance_collector.start_optimization(
            protein=target_protein,
            organism=organism,
            constraints=[
                "GCInRange", "NoRestrictionSite", "NoCrypticSplice",
                "NoStopCodons", "CodonAdapted",
            ],
            solver_backend="csp" if use_csp_solver else "greedy",
            seed=seed,
        )

    # ── Convergence tracking ────────────────────────────────────────
    _convergence = ConvergenceTracker()
    _convergence_status: str = "converged"  # default; overridden if caps are hit

    # ── CSP Solver path ─────────────────────────────────────────────
    if use_csp_solver:
        try:
            from ..solver.dispatch import csp_optimize
            from ..solver.types import SolverConfig
            restriction_sites = []
            if enzymes:
                from ..restriction_sites import get_recognition_site
                for enz in enzymes:
                    site = get_recognition_site(enz)
                    if site:
                        restriction_sites.append(site)
            else:
                for enz in ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]:
                    from ..restriction_sites import get_recognition_site as _grs
                    site = _grs(enz)
                    if site:
                        restriction_sites.append(site)

            config = SolverConfig(
                gc_lo=gc_lo,
                gc_hi=gc_hi,
                cryptic_splice_threshold=kwargs.get("splice_low", 3.0),
                restriction_sites=restriction_sites,
            )
            solver_result = csp_optimize(
                protein=target_protein,
                organism=organism,
                config=config,
            )
            if solver_result.solved and not solver_result.fallback_used:
                from ..scanner import gc_content as _gc_content
                from ..translation import compute_cai as _compute_cai
                seq = solver_result.sequence
                gc = _gc_content(seq)
                _canonical_organism = _species_key_to_organism(organism)
                try:
                    cai_val = _compute_cai(seq, _canonical_organism)
                except Exception:
                    logger.debug(
                        "CAI computation failed for organism '%s', using solver result",
                        organism, exc_info=True,
                    )
                    cai_val = solver_result.cai
                # UTR suggestions for CSP solver path
                _utr5_seq: str | None = None
                _utr3_seq: str | None = None
                _utr_score5: float | None = None
                _utr_score3: float | None = None
                if include_utr:
                    from ..utr_models import suggest_5utr as _s5, suggest_3utr as _s3
                    from ..utr_models import score_5utr as _sc5, score_3utr as _sc3
                    try:
                        _utr5_seq = _s5(organism)
                        _utr_score5 = _sc5(_utr5_seq, organism)
                    except ValueError:
                        logger.debug("No 5' UTR suggestion for organism '%s'", organism)
                    try:
                        _utr3_seq = _s3(organism)
                        _utr_score3 = _sc3(_utr3_seq, organism)
                    except ValueError:
                        logger.debug("No 3' UTR suggestion for organism '%s'", organism)

                # Compute CPB if requested (CSP solver path)
                _cpb_val: float | None = None
                if consider_codon_pair_bias:
                    try:
                        from ..codon_pair_scoring import compute_cpb as _compute_cpb
                        _cpb_val = round(_compute_cpb(seq, organism), 6)
                    except Exception:
                        logger.debug("CPB computation failed for CSP solver path", exc_info=True)

                # ── Translation verification (CSP solver path) ──
                from ..protein_verification import verify_and_raise as _verify_and_raise
                _verify_and_raise(seq, target_protein, organism=organism)

                # Evaluate custom objective score if a non-default objective is used
                _csp_objective_score: float | None = None
                if _objective_fn is not _cai_objective:
                    try:
                        _csp_objective_score = _objective_fn(seq, target_protein, organism)
                    except Exception:
                        _csp_objective_score = None

                return OptimizationResult(
                    sequence=seq,
                    gc_content=round(gc, 4),
                    cai=cai_val,
                    failed_predicates=[],
                    predicate_results=[],
                    certificate_text="CSP-Solver:" + str(solver_result.backend_used.value),
                    protein=target_protein,
                    fallback_used=False,
                    satisfied_predicates=["CSP_SOLVER"],
                    codon_pair_bias=_cpb_val,
                    suggested_5utr=_utr5_seq,
                    suggested_3utr=_utr3_seq,
                    utr_score_5=_utr_score5,
                    utr_score_3=_utr_score3,
                    convergence_status="converged",
                    iterations_used=1,
                    objective_score=_csp_objective_score,
                    biosecurity_screening_result=_biosecurity_result,
                )
        except Exception:
            logger.warning(
                "CSP solver failed, falling back to greedy optimizer",
                exc_info=True,
            )
    # ── End CSP Solver path ─────────────────────────────────────────

    # Map organism name to species key
    species_key = _organism_to_species_key(organism)

    # ── Organism-aware constraint selection ────────────────────────
    # When the target organism is prokaryotic, skip eukaryote-specific
    # constraints (cryptic splice sites, CpG islands, GT dinucleotide
    # avoidance) that are biologically inappropriate and unnecessarily
    # depress CAI.  This recovers ~0.27 CAI on prokaryotic targets
    # (see Task 6+8 benchmark findings).
    organism_domain = kwargs.get("organism_domain", "auto")
    if organism_domain not in ("auto", "eukaryote", "prokaryote"):
        organism_domain = "auto"

    if organism_domain == "auto":
        from ..organism_config import is_eukaryotic_organism
        is_prokaryote = not is_eukaryotic_organism(organism)
    elif organism_domain == "prokaryote":
        is_prokaryote = True
    else:
        is_prokaryote = False

    # PERF (Optimization A): Fast path for prokaryotes — skip all eukaryotic
    # constraint setup and evaluation.  Prokaryotes have no spliceosome,
    # so GT/AG avoidance, CpG island disruption, and cryptic splice
    # site evaluation are biologically irrelevant.  Skipping them
    # eliminates the expensive MaxEntScan calls and eukaryotic predicate
    # checks that dominate runtime for short sequences.
    # For prokaryotes: skip splice/GT/CpG constraints entirely
    # For eukaryotes: apply all constraints (default)
    effective_avoid_gt = kwargs.get("avoid_gt", not is_prokaryote)
    effective_splice_low = kwargs.get("splice_low", 3.0 if not is_prokaryote else 999.0)
    effective_splice_high = kwargs.get("splice_high", 6.0 if not is_prokaryote else 999.0)

    # ── Harmonize strategy path (Claassens RCA method) ────────────────
    if strategy == "harmonize":
        from .codon_harmonization import (
            harmonize_codons as _harmonize_codons,
            harmonize_with_cai_fallback as _harmonize_with_cai_fallback,
            compute_harmonization_score as _compute_harmonization_score,
        )

        _source_org = source_organism or organism

        if harmonization_cai_weight > 0:
            optimized_seq = _harmonize_with_cai_fallback(
                target_protein,
                source_organism=_source_org,
                target_organism=organism,
                cai_weight=harmonization_cai_weight,
            )
        else:
            optimized_seq = _harmonize_codons(
                target_protein,
                source_organism=_source_org,
                target_organism=organism,
            )

        # Compute metrics
        from ..scanner import gc_content as _gc_content_harm
        from ..translation import compute_cai as _compute_cai_harm
        gc = _gc_content_harm(optimized_seq)
        try:
            cai_val = _compute_cai_harm(optimized_seq, organism=organism)
        except Exception:
            logger.debug("CAI computation failed for harmonized sequence", exc_info=True)
            cai_val = 0.0

        # Compute harmonization score
        try:
            harm_score = _compute_harmonization_score(
                optimized_seq,
                source_organism=_source_org,
                target_organism=organism,
            )
        except Exception:
            harm_score = 0.0

        # Evaluate predicates
        from ..type_system import PredicateResult as _PR_harm
        pred_results = []

        # NoStopCodons
        from ..type_system import check_no_stop_codons as _check_no_stop
        pred_results.append(_check_no_stop(optimized_seq))

        # NoRestrictionSite
        from ..type_system import check_no_restriction_site as _check_no_rs
        pred_results.append(_check_no_rs(
            optimized_seq, enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
        ))

        # GCInRange
        gc_ok = gc_lo <= gc <= gc_hi
        pred_results.append(_PR_harm(
            "GCInRange", gc_ok,
            details=f"GC content: {gc:.3f} (range [{gc_lo}, {gc_hi}])"
        ))

        # CodonOptimality
        all_optimal = cai_val >= cai_threshold
        pred_results.append(_PR_harm(
            "CodonOptimality", all_optimal,
            details=f"CAI={cai_val:.4f}, min={cai_threshold}"
        ))

        # Skip eukaryotic constraints for prokaryotes
        if is_prokaryote:
            pred_results.append(_PR_harm(
                "NoCrypticSplice", True,
                details="Skipped for prokaryotic organism"
            ))
            pred_results.append(_PR_harm(
                "NoCpGIsland", True,
                details="Skipped for prokaryotic organism"
            ))
            pred_results.append(_PR_harm(
                "NoGTDinucleotide", True,
                details="Skipped for prokaryotic organism"
            ))
        else:
            # Eukaryotic predicate checks
            from ..type_system import check_no_cryptic_splice as _check_splice
            from ..type_system import check_no_cpg_island as _check_cpg
            pred_results.append(_check_splice(
                optimized_seq,
                low_thresh=effective_splice_low,
                high_thresh=effective_splice_high,
            ))
            pred_results.append(_check_cpg(optimized_seq))
            if effective_avoid_gt:
                pred_results.append(_PR_harm(
                    "NoGTDinucleotide", True,
                    details="GT avoidance not applied for harmonized sequences"
                ))
            else:
                pred_results.append(_PR_harm(
                    "NoGTDinucleotide", True,
                    details="GT avoidance disabled"
                ))

        # ValidCodingSeq
        from ..type_system import check_valid_coding_seq as _check_valid
        pred_results.append(_check_valid(optimized_seq))

        cert_text = format_certificate(pred_results, optimized_seq, species_key)

        failed = [r.predicate for r in pred_results if not r.passed]
        satisfied = [r.predicate for r in pred_results if r.passed]

        # Translation verification
        from ..protein_verification import verify_and_raise as _verify_harm
        _verify_harm(optimized_seq, target_protein, organism=organism)

        # UTR suggestions
        _utr5_seq = None
        _utr3_seq = None
        _utr_score5 = None
        _utr_score3 = None
        if include_utr:
            from ..utr_models import suggest_5utr as _s5h, suggest_3utr as _s3h
            from ..utr_models import score_5utr as _sc5h, score_3utr as _sc3h
            try:
                _utr5_seq = _s5h(organism)
                _utr_score5 = _sc5h(_utr5_seq, organism)
            except ValueError:
                logger.debug("No 5' UTR suggestion for organism '%s'", organism)
            try:
                _utr3_seq = _s3h(organism)
                _utr_score3 = _sc3h(_utr3_seq, organism)
            except ValueError:
                logger.debug("No 3' UTR suggestion for organism '%s'", organism)

        # Build provenance record
        _harm_provenance_record = None
        try:
            from ..provenance import OptimizationRecord as _OptRec_harm
            from ..provenance import _get_biocompiler_version as _get_ver_harm
            _harm_provenance_record = _OptRec_harm(
                input_sequence=target_protein,
                output_sequence=optimized_seq,
                organism=organism,
                constraints_applied=sorted([r.predicate for r in pred_results]),
                mutations_made=[],
                solver_backend="harmonize",
                solve_time=round(_time.monotonic() - _start_time, 6),
                seed_used=seed,
                timestamp=datetime.now(timezone.utc).isoformat(),
                biocompiler_version=_get_ver_harm(),
            )
        except Exception:
            logger.debug("Provenance record creation failed for harmonize strategy", exc_info=True)

        # Strict mode check
        if strict_mode and failed:
            from ..exceptions import OptimizationConstraintError
            raise OptimizationConstraintError(
                f"Harmonized sequence fails {len(failed)} predicate(s): {failed}",
                failed_predicates=failed,
            )

        return OptimizationResult(
            sequence=optimized_seq,
            gc_content=round(gc, 4),
            cai=cai_val,
            failed_predicates=failed,
            predicate_results=pred_results,
            certificate_text=cert_text,
            protein=target_protein,
            fallback_used=False,
            satisfied_predicates=satisfied,
            codon_pair_bias=None,
            suggested_5utr=_utr5_seq,
            suggested_3utr=_utr3_seq,
            utr_score_5=_utr_score5,
            utr_score_3=_utr_score3,
            convergence_status="converged",
            iterations_used=1,
            objective_score=harm_score,
            provenance=_harm_provenance_record,
            biosecurity_screening_result=_biosecurity_result,
        )

    # ── Hybrid strategy path (v10 default) ────────────────────────
    if strategy == "hybrid":
        from ..hybrid_optimizer import HybridOptimizer as _HybridOptimizer

        hybrid_opt = _HybridOptimizer(
            species=species_key,
            organism=organism,
            enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            avoid_gt=effective_avoid_gt,
            splice_threshold=effective_splice_low,
            provenance_collector=_provenance_collector,
        )

        hybrid_result = hybrid_opt.optimize(
            target_protein, is_prokaryote=is_prokaryote
        )
        optimized_seq = hybrid_result.sequence

        # ── Eukaryotic CAI recovery pass ──
        # For eukaryotes, the HybridOptimizer aggressively avoids GT
        # dinucleotides, which can replace optimal CAI codons with
        # suboptimal GT-free alternatives (e.g., yeast Leu: TTG→TTA,
        # Gly: GGT→GGA, Cys: TGT→TGC).  For eukaryotes, GT avoidance
        # in CDS is MUCH less important than CAI — GT dinucleotides
        # only matter when they form strong cryptic splice donor sites
        # (score_splice_donor_potential >= 0.5), and in-codon GTs from
        # optimal codons are biologically acceptable.
        #
        # This recovery pass swaps suboptimal codons back to optimal
        # ones when the CAI cost exceeds the threshold, using cost-aware
        # GT resolution with splice donor scoring.  GTs with low splice
        # donor potential are accepted in favor of CAI.
        _cai_recovery_upgrades = 0
        if not is_prokaryote and effective_avoid_gt:
            _usage = CODON_ADAPTIVENESS_TABLES.get(organism)
            if _usage is not None:
                optimized_seq, _cai_recovery_upgrades = _eukaryote_cai_recovery(
                    optimized_seq,
                    target_protein,
                    _usage,
                    enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                )
                if _cai_recovery_upgrades > 0:
                    logger.info(
                        "Eukaryotic CAI recovery: upgraded %d codon(s) to "
                        "optimal (GT accepted as soft preference)",
                        _cai_recovery_upgrades,
                    )

        # ── Systematic CpG elimination pass (eukaryotes only) ──
        # After CAI maximization and recovery, systematically eliminate CpG
        # dinucleotides that contribute to CpG islands.  This is a proper
        # optimization pass (not best-effort) that continues past unfixable
        # positions and reports remaining CpGs as warnings.
        _cpg_warnings: list[str] = []
        _cpg_seq_changed = False
        if not is_prokaryote:
            _usage_cpg = CODON_ADAPTIVENESS_TABLES.get(organism)
            if _usage_cpg is not None:
                optimized_seq, _cpg_warnings = _eliminate_cpg_dinucleotides(
                    optimized_seq,
                    target_protein,
                    _usage_cpg,
                    enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                    organism=organism,
                    gc_lo=gc_lo,
                    gc_hi=gc_hi,
                    max_cai_cost=0.05,  # Preserve high CAI; CpG avoidance is soft for CDS
                )
                _seq_before_cpg = hybrid_result.sequence
                _cpg_seq_changed = (optimized_seq != _seq_before_cpg)
                if _cpg_warnings:
                    for _w in _cpg_warnings:
                        logger.warning(_w)
                if _cpg_seq_changed:
                    logger.info(
                        "CpG elimination: modified sequence to avoid CpG islands"
                    )

        # ── Ultra-fast prokaryote predicate evaluation ──
        # For prokaryotes, skip the heavyweight BioOptimizer predicate
        # evaluation and use a streamlined version that only checks
        # prokaryote-relevant predicates (no splice, no CpG, no GT).
        # This avoids creating a BioOptimizer and running the full
        # _evaluate_all_predicates with 12 predicates.
        if is_prokaryote:
            from ..type_system import (
                check_no_stop_codons, check_valid_coding_seq,
                check_no_restriction_site,
            )
            import math as _pmath

            # Lightweight predicate evaluation — only prokaryote-relevant
            pred_results = []

            # 1. NoStopCodons
            pred_results.append(check_no_stop_codons(optimized_seq))

            # 2. NoCrypticSplice — skipped for prokaryotes
            pred_results.append(PredicateResult(
                "NoCrypticSplice", True,
                details="Skipped for prokaryotic organism"
            ))

            # 3. NoCpGIsland — skipped for prokaryotes
            pred_results.append(PredicateResult(
                "NoCpGIsland", True,
                details="Skipped for prokaryotic organism"
            ))

            # 4. NoRestrictionSite
            pred_results.append(check_no_restriction_site(
                optimized_seq, enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
            ))

            # 5. NoGTDinucleotide — skipped for prokaryotes
            pred_results.append(PredicateResult(
                "NoGTDinucleotide", True,
                details="Skipped for prokaryotic organism"
            ))

            # 6. ValidCodingSeq
            pred_results.append(check_valid_coding_seq(optimized_seq))

            # 7. ConservationScore — trivially passes (no mutagenesis)
            pred_results.append(PredicateResult(
                "ConservationScore", True,
                details=f"All AA conservation scores >= 0"
            ))

            # 8. CodonOptimality — use hybrid_result CAI
            cai_val = hybrid_result.cai
            all_optimal = cai_val >= cai_threshold
            pred_results.append(PredicateResult(
                "CodonOptimality", all_optimal,
                details=f"CAI={cai_val:.4f}, min={cai_threshold}"
            ))

            # 9. GCInRange
            _gc_val = hybrid_result.gc_content
            gc_ok = gc_lo <= _gc_val <= gc_hi
            pred_results.append(PredicateResult(
                "GCInRange", gc_ok,
                details=f"GC content: {_gc_val:.3f} (range [{gc_lo}, {gc_hi}])"
            ))

            # 10. SlidingGC (prokaryotic fast path)
            _eff_gc_w_min = gc_window_min if gc_window_min is not None else gc_lo
            _eff_gc_w_max = gc_window_max if gc_window_max is not None else gc_hi
            if gc_window_size > 0 and len(optimized_seq) >= gc_window_size:
                _sgc_result = check_sliding_gc(
                    optimized_seq, window_size=gc_window_size,
                    gc_min=_eff_gc_w_min, gc_max=_eff_gc_w_max,
                )
                pred_results.append(PredicateResult(
                    "SlidingGC", _sgc_result.passed,
                    details=(
                        f"Window={gc_window_size}, range=[{_eff_gc_w_min:.2f}, {_eff_gc_w_max:.2f}], "
                        f"min_gc={_sgc_result.min_gc:.3f}, max_gc={_sgc_result.max_gc:.3f}, "
                        f"violations={len(_sgc_result.violations)}"
                    ),
                ))
            else:
                pred_results.append(PredicateResult(
                    "SlidingGC", True,
                    details="Sliding-window GC check skipped (window_size=0 or sequence too short)"
                ))

            cert_text = format_certificate(pred_results, optimized_seq, species_key)

            # Skip all post-processing for prokaryotes (CAI recovery,
            # mRNA stability, CPB, UTR) since we start with max CAI
            # and prokaryote-specific constraints are already handled.
            gc = hybrid_result.gc_content
            failed = [r.predicate for r in pred_results if not r.passed]
            satisfied = [r.predicate for r in pred_results if r.passed]
            fallback = False

            # Import compute_cai for the hybrid path
            from ..translation import compute_cai

            # Build provenance record (lightweight for prokaryotes)
            from ..provenance import OptimizationRecord as _OptimizationRecord
            from ..provenance import _get_biocompiler_version as _get_version

            provenance_record = _OptimizationRecord(
                input_sequence=target_protein,
                output_sequence=optimized_seq,
                organism=organism,
                constraints_applied=sorted([r.predicate for r in pred_results]),
                mutations_made=[],
                solver_backend="hybrid-prok-fast",
                solve_time=round(_time.monotonic() - _start_time, 6),
                seed_used=seed,
                timestamp=datetime.now(timezone.utc).isoformat(),
                biocompiler_version=_get_version(),
            )

            # Finalize decision provenance trail if tracking is enabled
            _prok_decision_trail: OptimizationDecisionTrail | None = None
            if _provenance_collector is not None:
                try:
                    # Record constraint decisions from predicate results
                    for _pr in pred_results:
                        _action = "satisfied" if _pr.passed else "conflicted"
                        _positions_affected: list[int] = list(range(len(optimized_seq) // 3))
                        _details = _pr.details or ""
                        _provenance_collector.record_constraint_decision(ConstraintDecision(
                            constraint_name=_pr.predicate,
                            constraint_type="hard",
                            action_taken=_action,
                            positions_affected=_positions_affected[:1] if _positions_affected else [],
                            tradeoff_description=_details,
                            impact_on_cai=0.0,
                        ))
                    _prok_decision_trail = _provenance_collector.finalize(
                        output_dna=optimized_seq,
                        cai=cai_val,
                        gc=gc,
                    )
                except Exception:
                    logger.debug("Provenance finalization failed for prokaryote path", exc_info=True)

            # ── Translation verification (prokaryote fast path) ──
            from ..protein_verification import verify_and_raise as _verify_prok
            _verify_prok(optimized_seq, target_protein, organism=organism)

            # ── Custom objective refinement (prokaryote fast path) ──
            _prok_objective_score: float | None = None
            if _objective_fn is not _cai_objective:
                _obj_aas = protein_to_aa_list(target_protein)
                _obj_n_codons = len(_obj_aas)
                _obj_usage = CODON_ADAPTIVENESS_TABLES.get(organism, {})
                _obj_sorted_codons: dict[str, list[str]] = {}
                for _aa in set(_obj_aas):
                    if _aa == "*":
                        continue
                    _obj_sorted_codons[_aa] = sorted(
                        AA_TO_CODONS.get(_aa, []),
                        key=lambda c: _obj_usage.get(c, 0.0),
                        reverse=True,
                    )
                _obj_rs_sites: list[tuple[str, str]] = []
                for _enz in (enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]):
                    from ..restriction_sites import get_recognition_site as _grs2
                    _obj_site = _grs2(_enz)
                    if _obj_site:
                        _obj_rs_sites.append((_obj_site, reverse_complement(_obj_site)))

                _MAX_PROK_OBJ_ITERS = 5
                for _obj_iter in range(_MAX_PROK_OBJ_ITERS):
                    _obj_improved = False
                    _current_score = _objective_fn(optimized_seq, target_protein, organism)

                    for _ci in range(_obj_n_codons):
                        _aa = _obj_aas[_ci]
                        if _aa == "*" or _aa == "M":
                            continue
                        _current_codon = optimized_seq[_ci * 3:_ci * 3 + 3]
                        _best_alt = None
                        _best_alt_score = _current_score

                        for _alt in _obj_sorted_codons.get(_aa, []):
                            if _alt == _current_codon:
                                continue
                            _test_seq = optimized_seq[:_ci * 3] + _alt + optimized_seq[_ci * 3 + 3:]

                            # Hard constraint checks
                            _rs_ok = True
                            for _site, _site_rc in _obj_rs_sites:
                                if _site in _test_seq and _site not in optimized_seq:
                                    _rs_ok = False
                                    break
                                if _site_rc and _site_rc in _test_seq and _site_rc not in optimized_seq:
                                    _rs_ok = False
                                    break
                            if not _rs_ok:
                                continue

                            _test_gc = (_test_seq.count("G") + _test_seq.count("C")) / max(len(_test_seq), 1)
                            if not (gc_lo <= _test_gc <= gc_hi):
                                continue

                            _has_premature_stop = False
                            for _si in range(0, len(_test_seq) - 5, 3):
                                if _test_seq[_si:_si + 3] in ("TAA", "TAG", "TGA"):
                                    _has_premature_stop = True
                                    break
                            if _has_premature_stop:
                                continue

                            if _test_seq.count("ATTTA") > optimized_seq.count("ATTTA"):
                                continue

                            try:
                                _alt_score = _objective_fn(_test_seq, target_protein, organism)
                            except Exception:
                                continue

                            if _alt_score > _best_alt_score:
                                _best_alt = _alt
                                _best_alt_score = _alt_score

                        if _best_alt is not None:
                            optimized_seq = (
                                optimized_seq[:_ci * 3] + _best_alt + optimized_seq[_ci * 3 + 3:]
                            )
                            _obj_improved = True

                    if not _obj_improved:
                        break

                try:
                    _prok_objective_score = _objective_fn(optimized_seq, target_protein, organism)
                except Exception:
                    _prok_objective_score = None

                # Recompute metrics after objective refinement
                gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
                gc = round(gc, 4)
                cai_val = compute_cai(optimized_seq, organism)

            return OptimizationResult(
                sequence=optimized_seq,
                gc_content=gc,
                cai=cai_val,
                failed_predicates=failed,
                predicate_results=pred_results,
                certificate_text=cert_text,
                protein=target_protein,
                fallback_used=False,
                satisfied_predicates=satisfied,
                provenance=provenance_record,
                decision_trail=_prok_decision_trail,
                convergence_status="converged",
                iterations_used=1,
                objective_score=_prok_objective_score,
                biosecurity_screening_result=_biosecurity_result,
            )

        # Fast eukaryotic predicate evaluation — skip BioOptimizer creation
        # (~0.4ms saved per call) by directly checking only the predicates
        # that matter after hybrid optimization.
        from ..type_system import (
            check_no_stop_codons, check_valid_coding_seq,
            check_no_restriction_site, check_no_cryptic_splice,
            check_no_cpg_island,
        )

        pred_results = []

        # 1. NoStopCodons
        pred_results.append(check_no_stop_codons(optimized_seq))

        # 2. NoCrypticSplice — skip if already validated by hybrid optimizer
        if hybrid_result.splice_sites_validated:
            pred_results.append(PredicateResult(
                "NoCrypticSplice", True,
                details="Validated during optimization (MaxEntScan Phase 3)"
            ))
        else:
            pred_results.append(check_no_cryptic_splice(
                optimized_seq, low_thresh=effective_splice_low, organism=organism
            ))

        # 3. NoCpGIsland
        pred_results.append(check_no_cpg_island(optimized_seq, organism=organism))

        # 4. NoRestrictionSite
        pred_results.append(check_no_restriction_site(
            optimized_seq, enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
        ))

        # 5. NoGTDinucleotide (soft for eukaryotes — in-codon GTs from
        # optimal codons and cross-codon GTs with high CAI cost are acceptable)
        # Uses check_no_gt_dinucleotide_soft which returns:
        #   PASS          — GT count ≤ max_gt_count (auto-computed per sequence length)
        #   LIKELY_FAIL   — GT count > max_gt_count for eukaryotes (soft fail)
        #   FAIL          — any GT for prokaryotes (hard constraint)
        if effective_avoid_gt:
            pred_results.append(check_no_gt_dinucleotide_soft(
                optimized_seq, organism=organism,
            ))
        else:
            pred_results.append(PredicateResult(
                "NoGTDinucleotide", True,
                details="GT avoidance not requested"
            ))

        # 6. ValidCodingSeq
        pred_results.append(check_valid_coding_seq(optimized_seq))

        # 7. ConservationScore
        pred_results.append(PredicateResult(
            "ConservationScore", True,
            details=f"All AA conservation scores >= 0"
        ))

        # 8. CodonOptimality
        # Recompute CAI after eukaryotic recovery and CpG elimination passes
        # (which may have changed codons)
        from ..translation import compute_cai as _compute_cai_hybrid
        if _cai_recovery_upgrades > 0 or _cpg_seq_changed:
            cai_val = _compute_cai_hybrid(optimized_seq, organism)
        else:
            cai_val = hybrid_result.cai
        all_optimal = cai_val >= cai_threshold
        pred_results.append(PredicateResult(
            "CodonOptimality", all_optimal,
            details=f"CAI={cai_val:.4f}, min={cai_threshold}"
        ))

        # 9. GCInRange
        # Recompute GC after CAI recovery/CpG elimination may have changed the sequence
        if _cai_recovery_upgrades > 0 or _cpg_seq_changed:
            _gc_val = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
        else:
            _gc_val = hybrid_result.gc_content
        gc_ok = gc_lo <= _gc_val <= gc_hi
        pred_results.append(PredicateResult(
            "GCInRange", gc_ok,
            details=f"GC content: {_gc_val:.3f} (range [{gc_lo}, {gc_hi}])"
        ))

        # 10. SlidingGC (eukaryotic path)
        _eff_gc_w_min = gc_window_min if gc_window_min is not None else gc_lo
        _eff_gc_w_max = gc_window_max if gc_window_max is not None else gc_hi
        if gc_window_size > 0 and len(optimized_seq) >= gc_window_size:
            _sgc_result = check_sliding_gc(
                optimized_seq, window_size=gc_window_size,
                gc_min=_eff_gc_w_min, gc_max=_eff_gc_w_max,
            )
            pred_results.append(PredicateResult(
                "SlidingGC", _sgc_result.passed,
                details=(
                    f"Window={gc_window_size}, range=[{_eff_gc_w_min:.2f}, {_eff_gc_w_max:.2f}], "
                    f"min_gc={_sgc_result.min_gc:.3f}, max_gc={_sgc_result.max_gc:.3f}, "
                    f"violations={len(_sgc_result.violations)}"
                ),
            ))
        else:
            pred_results.append(PredicateResult(
                "SlidingGC", True,
                details="Sliding-window GC check skipped (window_size=0 or sequence too short)"
            ))

        cert_text = format_certificate(pred_results, optimized_seq, species_key)

        # Create a lightweight opt stub for downstream code that checks opt._applied_mutagenesis
        opt = BioOptimizer.__new__(BioOptimizer)
        opt._applied_mutagenesis = hybrid_opt._applied_mutagenesis
        opt._original_protein = target_protein
        opt.organism_name = organism
        # Provide species_cai for fallback CAI computation in edge cases
        opt.species_cai = hybrid_opt._species_cai if hasattr(hybrid_opt, '_species_cai') else {}
        # Provide _compute_seq_cai fallback for unsupported organisms
        opt._compute_seq_cai = lambda seq: compute_cai(seq, organism)
        # Sliding-window GC parameters (for _evaluate_all_predicates)
        opt._gc_window_size = gc_window_size
        opt._gc_window_min = gc_window_min
        opt._gc_window_max = gc_window_max

        # Import compute_cai for the hybrid path
        from ..translation import compute_cai

    else:
        # ── Legacy strategies (constraint_first / cai_first) ──────────
        # Configure optimizer
        opt = BioOptimizer(
            species=species_key,
            enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            splice_low=effective_splice_low,
            splice_high=effective_splice_high,
            min_cai=cai_threshold,
            strategy=strategy,
            avoid_gt=effective_avoid_gt,
            organism_name=organism,
            organism_domain="prokaryote" if is_prokaryote else "eukaryote",
            gc_window_size=gc_window_size,
            gc_window_min=gc_window_min,
            gc_window_max=gc_window_max,
        )

        # Attach provenance collector to optimizer for codon-level tracking
        if _provenance_collector is not None:
            opt._provenance_collector = _provenance_collector

        # Back-translate protein to DNA for the optimizer
        from ..translation import compute_cai
        initial_seq = _back_translate_protein(
            target_protein, species_key,
            strategy=strategy,
            enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            is_eukaryote=not is_prokaryote,
        )

        # Run optimization
        optimized_seq, pred_results, cert_text = opt.optimize(initial_seq, strategy=strategy)

    # Compute metrics — use hybrid optimizer values when available to avoid recomputation
    if strategy == "hybrid" and not use_csp_solver:
        # Hybrid optimizer already computed GC and CAI accurately.
        # However, if the eukaryotic CAI recovery pass or CpG elimination
        # was applied, the sequence has changed and we must recompute both metrics.
        if _cai_recovery_upgrades > 0 or _cpg_seq_changed:
            gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
            gc = round(gc, 4)
            # cai_val was already recomputed in the predicate evaluation block above
        else:
            gc = hybrid_result.gc_content if 'hybrid_result' in dir() else (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
            cai_val = hybrid_result.cai if 'hybrid_result' in dir() else None
            gc = round(gc, 4)
    else:
        gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
        gc = round(gc, 4)
        cai_val = None

    # Collect failed and satisfied predicates
    failed = [r.predicate for r in pred_results if not r.passed]
    satisfied = [r.predicate for r in pred_results if r.passed]

    # Detect if mutagenesis fallback was used (any V->I or similar)
    fallback = bool(opt._applied_mutagenesis)

    # ── G-quadruplex check ─────────────────────────────────────────
    # Detect and fix G-quadruplex motifs (runs of Gs that can form
    # stable secondary structures, stalling polymerases and causing
    # genomic instability).  This runs before mRNA stability since
    # G4 motifs are structural issues that affect both DNA and RNA.
    if check_g4:
        try:
            from .g_quadruplex import detect_g4_motifs, fix_g4_issues
            g4_report = detect_g4_motifs(optimized_seq)
            if g4_report.motifs:
                optimized_seq = fix_g4_issues(optimized_seq, g4_report, organism=organism)
                logger.info("Fixed %d G-quadruplex motifs", len(g4_report.motifs))
        except ImportError:
            pass  # g_quadruplex module not available

    # ── mRNA stability improvement pass ────────────────────────────
    # PERF: Skip mRNA stability for short sequences (<300aa / <900bp)
    # where the stability improvement is negligible.
    # Also skip if no ATTTA motifs exist (nothing to remove).
    _MRNA_STABILITY_MIN_LENGTH_BP = 900  # 300aa * 3

    mrna_stability_score: float | None = None
    destabilizing_motifs_removed: int = 0
    stability_improvement: float | None = None

    # ── Advanced mRNA half-life prediction ─────────────────────────
    # Use the mrna_halflife submodule for organism-specific half-life
    # prediction when available.  Falls back to the legacy ATTTA/T-run
    # approach below if the submodule is not importable.
    _used_advanced_halflife = False
    if check_mrna_halflife and optimize_mrna_stability:
        try:
            from .mrna_halflife import predict_mrna_halflife, optimize_mrna_stability as _optimize_halflife
            halflife_report = predict_mrna_halflife(optimized_seq, organism=organism)
            if halflife_report.predicted_halflife_hours < halflife_threshold:
                optimized_seq = _optimize_halflife(
                    optimized_seq, organism=organism,
                    target_halflife=halflife_threshold,
                )
                logger.info(
                    "Optimized mRNA half-life from %.1fh",
                    halflife_report.predicted_halflife_hours,
                )
            _used_advanced_halflife = True
        except ImportError:
            pass  # Fall back to old ATTTA/T-run approach

    # PERF: Quick pre-check — if no ATTTA motifs exist, the mRNA
    # stability pass cannot find any destabilizing motifs to remove.
    _has_destabilizing_motifs = (
        len(optimized_seq) >= _MRNA_STABILITY_MIN_LENGTH_BP
        and optimize_mrna_stability
        and "ATTTA" in optimized_seq
        and not _used_advanced_halflife  # Skip legacy if advanced pass succeeded
    )

    if _has_destabilizing_motifs:
        from ..mrna_stability import score_mrna_stability, suggest_mutations_for_stability

        initial_stability = score_mrna_stability(optimized_seq, organism)
        suggestions = suggest_mutations_for_stability(optimized_seq, organism)

        # PERF (Optimization C): Precompute restriction site strings once
        # instead of calling get_recognition_site per suggestion.
        from ..restriction_sites import get_recognition_site as _get_site
        _rs_site_strings: list[str] = []
        for enz in (enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]):
            site = _get_site(enz)
            if site:
                _rs_site_strings.append(site)

        # PERF (Optimization D): Use list mutation instead of string
        # concatenation for sequence modifications.
        seq_list = list(optimized_seq)
        # Pre-compute original protein once instead of per-suggestion
        orig_protein = BioOptimizer._translate(optimized_seq)
        # PERF (Optimization F): Cache GC count for incremental updates
        gc_count_cached = sum(1 for b in optimized_seq if b in "GC")
        n_bases = len(optimized_seq)

        # Pre-compute CAI adaptiveness table for CAI-aware acceptance
        _mrna_cai_table = CODON_ADAPTIVENESS_TABLES.get(organism, {})

        for suggestion in suggestions:
            pos = suggestion['position']
            ci = pos // 3
            new_codon = suggestion['suggested_codon']
            codon_start = ci * 3

            if codon_start + 3 > n_bases:
                continue

            # PERF (Optimization D): List mutation instead of string concat
            old_codon_list = seq_list[codon_start:codon_start + 3]
            old_codon_str = "".join(old_codon_list)

            # CAI-aware check: only accept mRNA stability suggestions that
            # don't degrade CAI.  For prokaryotes (E. coli), the optimal
            # codons are AT-rich (AAA, GAA, GAT) which also match the
            # RNase E / AT-rich "destabilizing" motifs.  Swapping these to
            # GC-rich alternatives (AAG, GAG, GAC) dramatically reduces CAI
            # for negligible mRNA stability benefit.  Skip any suggestion
            # that would replace an optimal codon with a suboptimal one.
            old_cai_w = _mrna_cai_table.get(old_codon_str, 0.0)
            new_cai_w = _mrna_cai_table.get(new_codon, 0.0)
            if new_cai_w < old_cai_w:
                # This suggestion would degrade CAI — skip it.
                # The CAI recovery pass will handle any remaining
                # opportunities after all post-processing.
                continue

            seq_list[codon_start:codon_start + 3] = list(new_codon)
            test_seq = "".join(seq_list)

            # PERF: Verify same protein by checking codon translates to same AA
            # This is faster than translating the entire sequence.
            # The mRNA module guarantees synonymous codons, so we only
            # need to verify no premature stop codons were introduced.
            _has_stop = False
            for _si in range(0, n_bases - 2, 3):
                _c = test_seq[_si:_si + 3]
                if _c in ("TAA", "TAG", "TGA") and _si < n_bases - 3:
                    _has_stop = True
                    break
            if _has_stop:
                # Rollback
                seq_list[codon_start:codon_start + 3] = old_codon_list
                continue

            # PERF (Optimization F): Incremental GC update
            old_gc_in_codon = sum(1 for b in old_codon_list if b in "GC")
            new_gc_in_codon = sum(1 for b in new_codon if b in "GC")
            new_gc_count = gc_count_cached - old_gc_in_codon + new_gc_in_codon
            test_gc = new_gc_count / n_bases
            if not (gc_lo <= test_gc <= gc_hi):
                # Rollback
                seq_list[codon_start:codon_start + 3] = old_codon_list
                continue

            # PERF (Optimization C): Check precomputed restriction sites
            rs_violated = False
            for site in _rs_site_strings:
                if site in test_seq and site not in optimized_seq:
                    rs_violated = True
                    break
            if rs_violated:
                # Rollback
                seq_list[codon_start:codon_start + 3] = old_codon_list
                continue

            # Safe to apply
            optimized_seq = test_seq
            gc_count_cached = new_gc_count
            destabilizing_motifs_removed += 1
            logger.debug(
                "mRNA stability: replaced %s->%s at codon %d (removing %s)",
                suggestion.get('original_codon', '???'), new_codon, ci,
                suggestion.get('motif_removed', suggestion.get('motif', '???')),
            )

        final_stability = score_mrna_stability(optimized_seq, organism)
        # score_mrna_stability returns an MRNAStabilityScore dataclass
        # with an overall_score attribute; extract the float
        init_score = initial_stability.overall_score if hasattr(initial_stability, 'overall_score') else float(initial_stability)
        final_score = final_stability.overall_score if hasattr(final_stability, 'overall_score') else float(final_stability)
        mrna_stability_score = final_score
        if init_score > 0:
            stability_improvement = round(final_score - init_score, 6)
        else:
            stability_improvement = round(final_score, 6)

        logger.info(
            "mRNA stability: initial=%.4f final=%.4f improvement=%.4f motifs_removed=%d",
            init_score, final_score, stability_improvement,
            destabilizing_motifs_removed,
        )
    # ── End mRNA stability pass ────────────────────────────────────

    # ── LinearDesign MFE optimization ────────────────────────────
    # Use LinearDesign for joint MFE/CAI optimization when requested.
    # LinearDesign uses dynamic programming to find the minimum free
    # energy sequence that also maximizes CAI, controlled by the
    # lambda tradeoff parameter.  Requires the C++ binary to be
    # installed.  Falls through to standard MFE optimization on failure.
    if use_lineardesign:
        try:
            from .mfe_optimization import optimize_with_lineardesign
            ld_result = optimize_with_lineardesign(
                target_protein, lambda_val=lineardesign_lambda, organism=organism
            )
            if "error" not in ld_result and ld_result.get("sequence"):
                optimized_seq = ld_result["sequence"]
                logger.info(
                    "LinearDesign optimization: MFE=%s, CAI=%s",
                    ld_result.get("mfe", "N/A"),
                    ld_result.get("cai", "N/A"),
                )
        except Exception:
            pass  # Fall through to standard MFE optimization

    # ── Codon pair bias improvement pass ──────────────────────────
    cpb_score: float | None = None

    if consider_codon_pair_bias:
        from ..codon_pair_scoring import (
            compute_cpb as _compute_cpb,
            suggest_better_pair as _suggest_better_pair,
        )

        initial_cpb = _compute_cpb(optimized_seq, organism)
        logger.info("Codon pair bias: initial CPB=%.4f", initial_cpb)

        # Iterate over adjacent codon pairs and try to improve
        aas = protein_to_aa_list(target_protein)
        max_cpb_iterations = 3
        for _cpb_iter in range(max_cpb_iterations):
            any_improved = False
            for ci in range(len(aas) - 1):
                codon_start1 = ci * 3
                codon_start2 = (ci + 1) * 3
                current_c1 = optimized_seq[codon_start1:codon_start1 + 3]
                current_c2 = optimized_seq[codon_start2:codon_start2 + 3]
                aa1 = aas[ci]
                aa2 = aas[ci + 1]

                # Look up CAI weights for this organism
                cai_weights = CODON_ADAPTIVENESS_TABLES.get(organism)

                better = _suggest_better_pair(
                    current_c1, current_c2, aa1, aa2, organism,
                    cai_weights=cai_weights,
                    cai_weight=0.7,
                    cpb_weight=0.3,
                )
                if better is None:
                    continue

                new_c1, new_c2 = better

                # Apply the swap and verify hard constraints
                # PERF (Optimization D): Use list mutation instead of string concat
                _seq_list = list(optimized_seq)
                _seq_list[codon_start1:codon_start1 + 3] = list(new_c1)
                _seq_list[codon_start2:codon_start2 + 3] = list(new_c2)
                test_seq = "".join(_seq_list)

                # PERF: Verify no premature stop codons instead of
                # translating the entire sequence twice.
                _has_stop = False
                for _si in range(0, len(test_seq) - 2, 3):
                    _c = test_seq[_si:_si + 3]
                    if _c in ("TAA", "TAG", "TGA") and _si < len(test_seq) - 3:
                        _has_stop = True
                        break
                if _has_stop:
                    continue

                # Verify: GC still in range
                test_gc = (test_seq.count("G") + test_seq.count("C")) / max(len(test_seq), 1)
                if not (gc_lo <= test_gc <= gc_hi):
                    continue

                # Verify: no new restriction sites
                from ..restriction_sites import get_recognition_site as _get_site
                rs_violated = False
                for enz in (enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]):
                    site = _get_site(enz)
                    if site and site in test_seq and site not in optimized_seq:
                        rs_violated = True
                        break
                if rs_violated:
                    continue

                # Verify: no new stop codons
                from ..type_system import check_no_stop_codons as _check_stops
                stop_result = _check_stops(test_seq)
                if not stop_result.passed:
                    continue

                # Verify: no new cryptic splice sites (if threshold specified)
                cryptic_splice_threshold_cpb = kwargs.get("splice_low", 3.0)
                try:
                    from ..maxentscan import max_donor_score, max_acceptor_score
                    old_max_d = max_donor_score(optimized_seq)
                    old_max_a = max_acceptor_score(optimized_seq)
                    new_max_d = max_donor_score(test_seq)
                    new_max_a = max_acceptor_score(test_seq)
                    if (new_max_d > old_max_d and new_max_d >= cryptic_splice_threshold_cpb):
                        continue
                    if (new_max_a > old_max_a and new_max_a >= cryptic_splice_threshold_cpb):
                        continue
                except Exception:
                    logger.debug("MaxEntScan check failed during codon pair bias step", exc_info=True)

                # Safe to apply
                optimized_seq = test_seq
                any_improved = True
                logger.debug(
                    "CPB: swapped %s-%s -> %s-%s at codon pair %d-%d",
                    current_c1, current_c2, new_c1, new_c2, ci, ci + 1,
                )

            if not any_improved:
                break

        final_cpb = _compute_cpb(optimized_seq, organism)
        cpb_score = round(final_cpb, 6)
        logger.info(
            "Codon pair bias: initial=%.4f final=%.4f improvement=%.4f",
            initial_cpb, final_cpb, final_cpb - initial_cpb,
        )
    # ── End Codon pair bias pass ──────────────────────────────────

    # ── Ensemble TASEP ribosome simulation ────────────────────────
    # Use ensemble TASEP averaging for more accurate ribosome
    # dynamics prediction when the ribosome_simulation submodule is
    # available.  The ensemble approach runs multiple stochastic
    # simulations and averages the results, reducing noise from
    # single-run stochastic effects.
    _tasep_result = None
    if tasep_ensemble:
        try:
            from .ribosome_simulation import simulate_tasep_ensemble
            # Compute dwell times from codon adaptiveness table
            _dwell_usage = CODON_ADAPTIVENESS_TABLES.get(organism, {})
            _dwell_times = []
            for _si in range(0, len(optimized_seq) - 2, 3):
                _codon = optimized_seq[_si:_si + 3]
                _w = _dwell_usage.get(_codon, 0.5)
                # Convert adaptiveness to dwell time (inverse relationship)
                _dwell_times.append(1.0 / max(_w, 0.01))
            _tasep_result = simulate_tasep_ensemble(
                dwell_times=_dwell_times,
                n_runs=tasep_n_runs,
            )
            if _tasep_result is not None:
                logger.info(
                    "Ensemble TASEP: %d runs, mean ribosome density=%.4f",
                    tasep_n_runs,
                    getattr(_tasep_result, 'mean_ribosome_density', 0.0),
                )
        except ImportError:
            pass  # Fall back to single-run TASEP

    # ── NMD (nonsense-mediated decay) check ───────────────────────
    # For eukaryotic sequences, detect NMD triggers (premature stop
    # codons upstream of exon-exon junctions) that can drastically
    # reduce mRNA levels.  Skipped for prokaryotes (no exon-exon
    # junctions → no NMD pathway).
    if check_nmd and not is_prokaryote:
        try:
            from .rna_degradation import detect_nmd_triggers
            # Compute ORF boundaries for NMD detection
            _orf_start = 0
            _orf_end = len(optimized_seq) - 3  # Exclude stop codon
            nmd_signals = detect_nmd_triggers(
                optimized_seq, orf_start=_orf_start, orf_end=_orf_end
            )
            if nmd_signals:
                for signal in nmd_signals:
                    logger.info("NMD trigger: %s", signal.description)
        except ImportError:
            pass  # rna_degradation module not available

    # ── CAI Recovery Pass ──────────────────────────────────────────
    # After all post-processing (mRNA stability, CPB), some positions may
    # have suboptimal codons.  This pass systematically tries to upgrade
    # each suboptimal codon to the optimal synonym, checking that the
    # swap doesn't violate any hard constraint.  This is especially
    # important for prokaryotes where the mRNA stability pass can
    # introduce significant CAI regression by replacing AT-rich optimal
    # codons (AAA, GAA, GAT) with GC-rich suboptimal ones (AAG, GAG, GAC).
    _CAI_RECOVERY_MAX_ITERATIONS = 5
    _cai_recovery_table = CODON_ADAPTIVENESS_TABLES.get(organism, {})
    _cai_recovery_rs: list[tuple[str, str]] = []
    for _enz in (enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]):
        _rs_site = _get_site(_enz) if '_get_site' in dir() else None
        if _rs_site is None:
            from ..restriction_sites import get_recognition_site as _grs
            _rs_site = _grs(_enz)
        if _rs_site:
            _cai_recovery_rs.append((_rs_site, reverse_complement(_rs_site)))

    _total_cai_recovery_iterations = 0
    for _recovery_iter in range(_CAI_RECOVERY_MAX_ITERATIONS):
        _any_recovered = False
        _aas = protein_to_aa_list(target_protein)
        _n_codons = len(_aas)

        # Sort positions by CAI gap (biggest improvement first)
        _recovery_candidates = []
        for _ci in range(_n_codons):
            _aa = _aas[_ci]
            if _aa == "*" or _aa == "M":
                continue
            _current = optimized_seq[_ci * 3:_ci * 3 + 3]
            _current_w = _cai_recovery_table.get(_current, 0.0)
            _alternatives = AA_TO_CODONS.get(_aa, [])
            _best_w = max(_cai_recovery_table.get(c, 0.0) for c in _alternatives) if _alternatives else 0.0
            if _best_w > _current_w:
                _gap = _best_w - _current_w
                _recovery_candidates.append((_gap, _ci, _current, _aa))

        _recovery_candidates.sort(key=lambda x: x[0], reverse=True)

        for _gap, _ci, _current, _aa in _recovery_candidates:
            # Try codons in CAI order (highest first)
            _sorted_alts = sorted(
                AA_TO_CODONS.get(_aa, []),
                key=lambda c: _cai_recovery_table.get(c, 0.0),
                reverse=True,
            )
            for _alt in _sorted_alts:
                if _alt == _current:
                    continue
                _alt_w = _cai_recovery_table.get(_alt, 0.0)
                _cur_w = _cai_recovery_table.get(_current, 0.0)
                if _alt_w <= _cur_w:
                    break  # No improvement possible (sorted descending)

                _test_seq = optimized_seq[:_ci * 3] + _alt + optimized_seq[_ci * 3 + 3:]

                # Check: no new restriction sites
                _rs_ok = True
                for _site, _site_rc in _cai_recovery_rs:
                    if _site in _test_seq or (_site_rc and _site_rc in _test_seq):
                        _rs_ok = False
                        break
                if not _rs_ok:
                    continue

                # Check: GC still in range
                _test_gc = (_test_seq.count("G") + _test_seq.count("C")) / max(len(_test_seq), 1)
                if not (gc_lo <= _test_gc <= gc_hi):
                    continue

                # Check: no new ATTTA motifs
                if _test_seq.count("ATTTA") > optimized_seq.count("ATTTA"):
                    continue

                # Check: no 6+ T runs
                _max_t = 0
                _j = 0
                while _j < len(_test_seq):
                    if _test_seq[_j] == 'T':
                        _k = _j
                        while _k < len(_test_seq) and _test_seq[_k] == 'T':
                            _k += 1
                        if _k - _j > _max_t:
                            _max_t = _k - _j
                        _j = _k
                    else:
                        _j += 1
                if _max_t >= 6:
                    continue

                # Check: no new premature stop codons
                _has_premature_stop = False
                for _si in range(0, len(_test_seq) - 5, 3):
                    if _test_seq[_si:_si + 3] in ("TAA", "TAG", "TGA"):
                        _has_premature_stop = True
                        break
                if _has_premature_stop:
                    continue

                # Check: for eukaryotes, no new cryptic splice sites
                if not is_prokaryote:
                    try:
                        _old_max_d = max_donor_score(optimized_seq)
                        _old_max_a = max_acceptor_score(optimized_seq)
                        _new_max_d = max_donor_score(_test_seq)
                        _new_max_a = max_acceptor_score(_test_seq)
                        if (_new_max_d > _old_max_d and _new_max_d >= effective_splice_low):
                            continue
                        if (_new_max_a > _old_max_a and _new_max_a >= effective_splice_low):
                            continue
                    except Exception:
                        logger.debug("MaxEntScan check failed during CAI recovery", exc_info=True)

                # All checks passed — accept the upgrade
                optimized_seq = _test_seq
                _any_recovered = True
                logger.debug(
                    "CAI recovery: upgraded codon %d from %s to %s (w %.4f→%.4f)",
                    _ci, _current, _alt, _cur_w, _alt_w,
                )
                break  # Move to next position

            # ── Paired codon swap ──────────────────────────────────
            # If the single-codon upgrade was blocked (e.g., by a
            # restriction site that straddles two codons), try a paired
            # swap: change the target codon AND an adjacent codon.
            # Only consider the immediate neighbors and only if the
            # combined CAI improvement is positive.
            import math as _math
            for _adj_offset in (1, -1):
                _adj_ci = _ci + _adj_offset
                if _adj_ci < 0 or _adj_ci >= _n_codons:
                    continue
                _adj_aa = _aas[_adj_ci]
                if _adj_aa == "*" or _adj_aa == "M":
                    continue
                _adj_current = optimized_seq[_adj_ci * 3:_adj_ci * 3 + 3]
                _adj_current_w = _cai_recovery_table.get(_adj_current, 0.0)

                # Try each CAI-ordered alt for the target, paired with
                # each alt for the neighbor.
                for _alt in _sorted_alts:
                    if _alt == _current:
                        continue
                    _alt_w = _cai_recovery_table.get(_alt, 0.0)
                    if _alt_w <= _cur_w:
                        break

                    # Net CAI must improve: the gain at the target must
                    # outweigh any loss at the adjacent position.
                    _adj_sorted = sorted(
                        AA_TO_CODONS.get(_adj_aa, []),
                        key=lambda c: _cai_recovery_table.get(c, 0.0),
                        reverse=True,
                    )
                    for _adj_alt in _adj_sorted:
                        if _adj_alt == _adj_current:
                            continue
                        _adj_alt_w = _cai_recovery_table.get(_adj_alt, 0.0)
                        # Net log-CAI change must be positive
                        _old_log = _math.log(max(_cur_w, 1e-10)) + _math.log(max(_adj_current_w, 1e-10))
                        _new_log = _math.log(max(_alt_w, 1e-10)) + _math.log(max(_adj_alt_w, 1e-10))
                        if _new_log <= _old_log:
                            continue

                        # Build test sequence with both swaps
                        _lo_ci = min(_ci, _adj_ci)
                        _hi_ci = max(_ci, _adj_ci)
                        _test_seq = (
                            optimized_seq[:_lo_ci * 3]
                            + (_alt if _lo_ci == _ci else _adj_alt)
                            + optimized_seq[_lo_ci * 3 + 3:_hi_ci * 3]
                            + (_alt if _hi_ci == _ci else _adj_alt)
                            + optimized_seq[_hi_ci * 3 + 3:]
                        )

                        # Validate: no new restriction sites
                        _rs_ok = True
                        for _site, _site_rc in _cai_recovery_rs:
                            if _site in _test_seq or (_site_rc and _site_rc in _test_seq):
                                _rs_ok = False
                                break
                        if not _rs_ok:
                            continue

                        # Validate: GC in range
                        _test_gc = (_test_seq.count("G") + _test_seq.count("C")) / max(len(_test_seq), 1)
                        if not (gc_lo <= _test_gc <= gc_hi):
                            continue

                        # Validate: no new ATTTA
                        if _test_seq.count("ATTTA") > optimized_seq.count("ATTTA"):
                            continue

                        # Validate: no 6+ T runs
                        _max_t = 0
                        _j = 0
                        while _j < len(_test_seq):
                            if _test_seq[_j] == 'T':
                                _k = _j
                                while _k < len(_test_seq) and _test_seq[_k] == 'T':
                                    _k += 1
                                if _k - _j > _max_t:
                                    _max_t = _k - _j
                                _j = _k
                            else:
                                _j += 1
                        if _max_t >= 6:
                            continue

                        # Validate: no premature stops
                        _has_premature_stop = False
                        for _si in range(0, len(_test_seq) - 5, 3):
                            if _test_seq[_si:_si + 3] in ("TAA", "TAG", "TGA"):
                                _has_premature_stop = True
                                break
                        if _has_premature_stop:
                            continue

                        # All checks passed — accept the paired upgrade
                        optimized_seq = _test_seq
                        _any_recovered = True
                        logger.debug(
                            "CAI recovery: paired upgrade codon %d (%s→%s) + "
                            "codon %d (%s→%s)",
                            _ci, _current, _alt,
                            _adj_ci, _adj_current, _adj_alt,
                        )
                        break  # Accept first valid paired swap for this target
                    if _any_recovered:
                        break
                if _any_recovered:
                    break

        if not _any_recovered:
            break  # No more improvements possible
        _total_cai_recovery_iterations += 1

        # ── Convergence check ──
        # After each CAI recovery iteration, compute the objective and
        # check if the optimizer has converged or is oscillating.
        _recovery_cai = compute_cai(optimized_seq, organism)
        _recovery_gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
        _recovery_gc_ok = 1.0 if gc_lo <= _recovery_gc <= gc_hi else 0.0
        _recovery_obj = _recovery_cai * _recovery_gc_ok  # 0 if GC out of range
        _convergence.record(_recovery_obj)

        _conv_status = _convergence.check_convergence()
        if _conv_status is not None:
            _convergence_status = _conv_status
            logger.info(
                "CAI recovery: %s after %d iterations (objective=%.6f)",
                _conv_status, _recovery_iter + 1, _recovery_obj,
            )
            if _conv_status == "oscillating":
                # Stop at the best point — no need to continue cycling
                break
    # ── End CAI Recovery Pass ──────────────────────────────────────

    # ── Final CpG elimination pass (eukaryotes only) ──────────────
    # This runs AFTER the CAI recovery pass because the recovery may
    # have re-introduced CG dinucleotides by upgrading suboptimal
    # CG-free codons to optimal CG-containing codons.  We need to
    # eliminate those CGs one final time to ensure the sequence passes
    # the NoCpGIsland predicate.
    _final_cpg_warnings: list[str] = []
    if not is_prokaryote:
        _final_cpg_usage = CODON_ADAPTIVENESS_TABLES.get(organism)
        if _final_cpg_usage is not None:
            optimized_seq, _final_cpg_warnings = _eliminate_cpg_dinucleotides(
                optimized_seq,
                target_protein,
                _final_cpg_usage,
                enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                organism=organism,
                gc_lo=gc_lo,
                gc_hi=gc_hi,
                max_cai_cost=0.05,  # Preserve high CAI; CpG avoidance is soft for CDS
            )
            if _final_cpg_warnings:
                for _fw in _final_cpg_warnings:
                    logger.warning(_fw)

    # ── Sliding-window GC constraint ─────────────────────────────────
    # After global GC is satisfied, check that no local window has
    # extreme GC content.  This prevents polymerase stalling and
    # secondary-structure hotspots that the global constraint misses.
    #
    # IMPORTANT: Skip sliding-window GC fix for organisms where the
    # natural coding GC is far from the window bounds.  For example,
    # yeast coding GC is ~38% — a 50bp window with bounds [0.30, 0.70]
    # will flag many AT-rich windows as violations and force codon swaps
    # that destroy CAI (drops of 0.10–0.20 are common).  The sliding GC
    # fix should only run when the organism's natural GC is well within
    # the window bounds.
    _effective_gc_window_min = gc_window_min if gc_window_min is not None else gc_lo
    _effective_gc_window_max = gc_window_max if gc_window_max is not None else gc_hi
    _sliding_gc_swaps = 0

    # Check if organism's natural GC is compatible with sliding window bounds
    _skip_sliding_gc = False
    _organism_gc_target = ORGANISM_GC_TARGETS.get(organism)
    if _organism_gc_target is not None and gc_window_size > 0:
        _org_gc_lo, _org_gc_hi = _organism_gc_target
        # Skip if the organism's natural GC range is outside the window bounds
        # by more than 5% — the sliding window fix would force unnatural GC
        if _org_gc_hi < _effective_gc_window_min + 0.05 or _org_gc_lo > _effective_gc_window_max - 0.05:
            _skip_sliding_gc = True
            logger.info(
                "Sliding-window GC: skipping for %s (natural GC range "
                "[%.2f, %.2f] is outside window bounds [%.2f, %.2f]) — "
                "preserving CAI takes priority",
                organism, _org_gc_lo, _org_gc_hi,
                _effective_gc_window_min, _effective_gc_window_max,
            )

    # For the hybrid strategy path, the HybridOptimizer already produces
    # sequences with organism-appropriate GC content. The sliding GC fix
    # has been shown to destroy CAI (drops of 0.10-0.40) by forcing
    # codon swaps to meet local GC window constraints that are biologically
    # inappropriate for the organism. Skip it entirely for hybrid strategy
    # to protect CAI.
    if strategy == "hybrid" and not use_csp_solver:
        _skip_sliding_gc = True

    if gc_window_size > 0 and len(optimized_seq) >= gc_window_size and not _skip_sliding_gc:
        _sliding_result = check_sliding_gc(
            optimized_seq,
            window_size=gc_window_size,
            gc_min=_effective_gc_window_min,
            gc_max=_effective_gc_window_max,
        )
        if not _sliding_result.passed:
            logger.info(
                "Sliding-window GC: %d violation(s) detected (window=%d, "
                "range=[%.2f, %.2f]). Attempting local fixes.",
                len(_sliding_result.violations), gc_window_size,
                _effective_gc_window_min, _effective_gc_window_max,
            )
            _sliding_usage = CODON_ADAPTIVENESS_TABLES.get(organism)
            optimized_seq, _sliding_gc_swaps = fix_sliding_gc_violations(
                optimized_seq,
                target_protein,
                window_size=gc_window_size,
                gc_min=_effective_gc_window_min,
                gc_max=_effective_gc_window_max,
                usage=_sliding_usage,
                gc_lo=gc_lo,
                gc_hi=gc_hi,
            )
            if _sliding_gc_swaps > 0:
                logger.info(
                    "Sliding-window GC: fixed %d local violation(s) with %d codon swap(s)",
                    len(_sliding_result.violations), _sliding_gc_swaps,
                )

    # ── Refresh SlidingGC predicate result after post-processing ──
    # The fix_sliding_gc_violations step may have modified optimized_seq,
    # so we must re-evaluate the SlidingGC predicate to ensure the result
    # in pred_results matches the final sequence.  Without this, certificate
    # verification fails because the recorded verdict was evaluated on the
    # pre-fix sequence while the certificate contains the post-fix sequence.
    if gc_window_size > 0 and len(optimized_seq) >= gc_window_size and _sliding_gc_swaps > 0:
        _effective_gc_window_min = gc_window_min if gc_window_min is not None else gc_lo
        _effective_gc_window_max = gc_window_max if gc_window_max is not None else gc_hi
        _refreshed_sgc = check_sliding_gc(
            optimized_seq,
            window_size=gc_window_size,
            gc_min=_effective_gc_window_min,
            gc_max=_effective_gc_window_max,
        )
        _refreshed_sgc_pr = PredicateResult(
            "SlidingGC", _refreshed_sgc.passed,
            details=(
                f"Window={gc_window_size}, range=[{_effective_gc_window_min:.2f}, {_effective_gc_window_max:.2f}], "
                f"min_gc={_refreshed_sgc.min_gc:.3f}, max_gc={_refreshed_sgc.max_gc:.3f}, "
                f"violations={len(_refreshed_sgc.violations)}"
            ),
        )
        # Replace the stale SlidingGC entry in pred_results
        pred_results = [
            _refreshed_sgc_pr if r.predicate == "SlidingGC" else r
            for r in pred_results
        ]
        # Refresh cert_text to match the updated predicates
        cert_text = format_certificate(pred_results, optimized_seq, species_key)
        failed = [r.predicate for r in pred_results if not r.passed]
        satisfied = [r.predicate for r in pred_results if r.passed]
    elif gc_window_size > 0:
        logger.debug(
            "Sliding-window GC: skipped (sequence length %d < window size %d)",
            len(optimized_seq), gc_window_size,
        )

    # ── Post-sliding-GC CpG re-elimination (eukaryotes only) ──────
    # The sliding-window GC fix (fix_sliding_gc_violations) may have
    # reintroduced CG dinucleotides by swapping codons without checking
    # for CpG creation.  Run one final CpG elimination pass to restore
    # the NoCpGIsland predicate without regressing sliding-window GC.
    if not is_prokaryote and _sliding_gc_swaps > 0:
        _post_sgc_usage = CODON_ADAPTIVENESS_TABLES.get(organism)
        if _post_sgc_usage is not None:
            optimized_seq, _post_sgc_warnings = _eliminate_cpg_dinucleotides(
                optimized_seq,
                target_protein,
                _post_sgc_usage,
                enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                organism=organism,
                gc_lo=gc_lo,
                gc_hi=gc_hi,
                max_cai_cost=0.05,  # Preserve high CAI; CpG avoidance is soft for CDS
            )
            if _post_sgc_warnings:
                for _psw in _post_sgc_warnings:
                    logger.warning(_psw)
            # If CpG elimination changed the sequence, re-check sliding GC
            _sgc_needs_refresh = _sliding_gc_swaps > 0
            if _sgc_needs_refresh and gc_window_size > 0 and len(optimized_seq) >= gc_window_size:
                _effective_gc_window_min = gc_window_min if gc_window_min is not None else gc_lo
                _effective_gc_window_max = gc_window_max if gc_window_max is not None else gc_hi
                _refreshed2_sgc = check_sliding_gc(
                    optimized_seq,
                    window_size=gc_window_size,
                    gc_min=_effective_gc_window_min,
                    gc_max=_effective_gc_window_max,
                )
                _refreshed2_sgc_pr = PredicateResult(
                    "SlidingGC", _refreshed2_sgc.passed,
                    details=(
                        f"Window={gc_window_size}, range=[{_effective_gc_window_min:.2f}, {_effective_gc_window_max:.2f}], "
                        f"min_gc={_refreshed2_sgc.min_gc:.3f}, max_gc={_refreshed2_sgc.max_gc:.3f}, "
                        f"violations={len(_refreshed2_sgc.violations)}"
                    ),
                )
                pred_results = [
                    _refreshed2_sgc_pr if r.predicate == "SlidingGC" else r
                    for r in pred_results
                ]
                cert_text = format_certificate(pred_results, optimized_seq, species_key)

    # ── Recompute final metrics after all post-processing passes ───
    # CAI and GC must be recomputed here because the mRNA stability
    # and codon-pair-bias passes may have modified optimized_seq.
    gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
    gc = round(gc, 4)

    # Use the canonical organism name (opt.organism_name) so that
    # compute_cai can look up the correct CODON_ADAPTIVENESS_TABLES
    # entry even when the caller passed a short alias like 'ecoli'.
    # PERF: Skip CAI recomputation if hybrid optimizer already provided it
    # and no post-processing passes (mRNA stability, CPB) modified the sequence.
    if cai_val is not None and not optimize_mrna_stability and not consider_codon_pair_bias:
        pass  # Use hybrid_result.cai — already set
    else:
        try:
            cai_val = compute_cai(optimized_seq, opt.organism_name)
        except UnsupportedOrganismError:
            logger.debug(
                "Unsupported organism '%s' for compute_cai, using species CAI table",
                opt.organism_name,
            )
            cai_val = opt._compute_seq_cai(optimized_seq)

    # Build provenance record for reproducibility
    from ..provenance import OptimizationRecord as _OptimizationRecord
    from ..provenance import _get_biocompiler_version as _get_version

    # Determine solver backend name
    solver_backend = "greedy"
    if use_csp_solver:
        solver_backend = "csp"  # actual backend may vary, but CSP was requested

    # Extract mutation descriptions from applied mutagenesis
    mutations_made: list[str] = []
    for mut in opt._applied_mutagenesis:
        if isinstance(mut, dict):
            mutations_made.append(mut.get("description", str(mut)))
        else:
            mutations_made.append(str(mut))

    # Constraints applied (from satisfied + failed predicate names)
    constraints_applied = sorted(set(
        [r.predicate for r in pred_results]
    ))

    provenance_record = _OptimizationRecord(
        input_sequence=target_protein,
        output_sequence=optimized_seq,
        organism=organism,
        constraints_applied=constraints_applied,
        mutations_made=mutations_made,
        solver_backend=solver_backend,
        solve_time=round(_time.monotonic() - _start_time, 6),
        seed_used=seed,
        timestamp=datetime.now(timezone.utc).isoformat(),
        biocompiler_version=_get_version(),
        mrna_stability_score=mrna_stability_score,
        destabilizing_motifs_removed=destabilizing_motifs_removed,
        stability_improvement=stability_improvement,
    )

    # ── UTR suggestions ─────────────────────────────────────────────
    utr5_seq: str | None = None
    utr3_seq: str | None = None
    utr_score5: float | None = None
    utr_score3: float | None = None

    if include_utr:
        from ..utr_models import suggest_5utr, suggest_3utr, score_5utr, score_3utr
        try:
            utr5_seq = suggest_5utr(organism)
            utr_score5 = score_5utr(utr5_seq, organism)
        except ValueError:
            logger.debug("No 5' UTR suggestion for organism '%s'", organism)
        try:
            utr3_seq = suggest_3utr(organism)
            utr_score3 = score_3utr(utr3_seq, organism)
        except ValueError:
            logger.debug("No 3' UTR suggestion for organism '%s'", organism)

    # Record constraint decisions if provenance tracking is enabled
    if _provenance_collector is not None:
        try:
            # Build constraint decisions from predicate results
            for _pr in pred_results:
                _action = "satisfied" if _pr.passed else "conflicted"
                # Estimate CAI impact (negative if constraint reduced CAI)
                _cai_impact = 0.0
                if not _pr.passed:
                    _cai_impact = -0.01  # Conservative estimate for failed constraints
                elif _pr.predicate in ("NoRestrictionSite", "NoCrypticSplice", "GCInRange"):
                    # These constraints typically have some CAI cost even when satisfied
                    _cai_impact = -0.005

                # Determine positions affected (best effort)
                _positions_affected: list[int] = []
                _details = getattr(_pr, "details", "") or ""
                if "position" in _details.lower():
                    # Try to extract position info from details
                    import re as _re
                    _pos_matches = _re.findall(r'pos(?:ition)?[:\s]+(\d+)', _details, _re.IGNORECASE)
                    _positions_affected = [int(p) for p in _pos_matches[:10]]

                _tradeoff = _details if _details else f"Constraint {_pr.predicate} was {_action}"

                _provenance_collector.record_constraint_decision(ConstraintDecision(
                    constraint_name=_pr.predicate,
                    constraint_type="hard",
                    action_taken=_action,
                    positions_affected=_positions_affected,
                    tradeoff_description=_tradeoff,
                    impact_on_cai=_cai_impact,
                ))
        except Exception:
            logger.debug("Constraint provenance recording failed", exc_info=True)

    # Finalize decision provenance trail if tracking is enabled
    _decision_trail: OptimizationDecisionTrail | None = None
    if _provenance_collector is not None:
        try:
            _decision_trail = _provenance_collector.finalize(
                output_dna=optimized_seq,
                cai=cai_val,
                gc=gc,
            )
        except Exception as _prov_exc:
            logger.debug("Provenance finalization failed: %s: %s", type(_prov_exc).__name__, _prov_exc, exc_info=True)
            _decision_trail = None

    # ── Final translation verification ──────────────────────────────
    # After all optimization passes (CAI recovery, mRNA stability, etc.),
    # verify that the optimized DNA still encodes the original protein.
    # This catches bugs in any optimization pass that might corrupt
    # the protein sequence.
    from ..protein_verification import verify_and_raise as _verify_final
    _verify_final(optimized_seq, target_protein, organism=organism)

    # ── Custom objective refinement pass ──────────────────────────────
    # When a non-default objective is provided, perform an additional
    # hill-climb pass that scores codon substitutions by the custom
    # objective rather than CAI alone.  This allows users to optimize
    # for GC balance, codon pair bias, or any custom metric while still
    # respecting all hard constraints.
    _objective_score: float | None = None
    if _objective_fn is not _cai_objective:
        logger.info(
            "Applying custom objective refinement pass for objective=%s",
            getattr(_objective_fn, '__name__', repr(_objective_fn)),
        )
        _obj_aas = protein_to_aa_list(target_protein)
        _obj_n_codons = len(_obj_aas)
        _obj_usage = CODON_ADAPTIVENESS_TABLES.get(organism, {})
        _obj_sorted_codons: dict[str, list[str]] = {}
        for _aa in set(_obj_aas):
            if _aa == "*":
                continue
            _obj_sorted_codons[_aa] = sorted(
                AA_TO_CODONS.get(_aa, []),
                key=lambda c: _obj_usage.get(c, 0.0),
                reverse=True,
            )
        _obj_rs_sites: list[tuple[str, str]] = []
        for _enz in (enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]):
            from ..restriction_sites import get_recognition_site as _grs2
            _obj_site = _grs2(_enz)
            if _obj_site:
                _obj_rs_sites.append((_obj_site, reverse_complement(_obj_site)))

        _MAX_OBJECTIVE_REFINEMENT_ITERS = 5
        for _obj_iter in range(_MAX_OBJECTIVE_REFINEMENT_ITERS):
            _obj_improved = False
            _current_score = _objective_fn(optimized_seq, target_protein, organism)

            for _ci in range(_obj_n_codons):
                _aa = _obj_aas[_ci]
                if _aa == "*" or _aa == "M":
                    continue
                _current_codon = optimized_seq[_ci * 3:_ci * 3 + 3]

                # Try each synonym for this position
                _best_alt = None
                _best_alt_score = _current_score

                for _alt in _obj_sorted_codons.get(_aa, []):
                    if _alt == _current_codon:
                        continue

                    _test_seq = optimized_seq[:_ci * 3] + _alt + optimized_seq[_ci * 3 + 3:]

                    # Hard constraint checks
                    # 1. No new restriction sites
                    _rs_ok = True
                    for _site, _site_rc in _obj_rs_sites:
                        if _site in _test_seq and _site not in optimized_seq:
                            _rs_ok = False
                            break
                        if _site_rc and _site_rc in _test_seq and _site_rc not in optimized_seq:
                            _rs_ok = False
                            break
                    if not _rs_ok:
                        continue

                    # 2. GC still in range
                    _test_gc = (_test_seq.count("G") + _test_seq.count("C")) / max(len(_test_seq), 1)
                    if not (gc_lo <= _test_gc <= gc_hi):
                        continue

                    # 3. No new premature stops
                    _has_premature_stop = False
                    for _si in range(0, len(_test_seq) - 5, 3):
                        if _test_seq[_si:_si + 3] in ("TAA", "TAG", "TGA"):
                            _has_premature_stop = True
                            break
                    if _has_premature_stop:
                        continue

                    # 4. No new ATTTA motifs
                    if _test_seq.count("ATTTA") > optimized_seq.count("ATTTA"):
                        continue

                    # Evaluate the custom objective
                    try:
                        _alt_score = _objective_fn(_test_seq, target_protein, organism)
                    except Exception:
                        continue

                    if _alt_score > _best_alt_score:
                        _best_alt = _alt
                        _best_alt_score = _alt_score

                if _best_alt is not None:
                    optimized_seq = (
                        optimized_seq[:_ci * 3] + _best_alt + optimized_seq[_ci * 3 + 3:]
                    )
                    _obj_improved = True

            if not _obj_improved:
                break

        # Record the final objective score
        try:
            _objective_score = _objective_fn(optimized_seq, target_protein, organism)
        except Exception:
            _objective_score = None

        # Recompute CAI and GC after objective refinement
        gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
        gc = round(gc, 4)
        try:
            cai_val = compute_cai(optimized_seq, opt.organism_name)
        except UnsupportedOrganismError:
            cai_val = opt._compute_seq_cai(optimized_seq)

    result = OptimizationResult(
        sequence=optimized_seq,
        gc_content=gc,
        cai=cai_val,
        failed_predicates=failed,
        predicate_results=pred_results,
        certificate_text=cert_text,
        protein=target_protein,
        fallback_used=fallback,
        satisfied_predicates=satisfied,
        aa_substitutions=opt._applied_mutagenesis,
        mrna_stability_score=mrna_stability_score,
        destabilizing_motifs_removed=destabilizing_motifs_removed,
        stability_improvement=stability_improvement,
        provenance=provenance_record,
        codon_pair_bias=cpb_score,
        suggested_5utr=utr5_seq,
        suggested_3utr=utr3_seq,
        utr_score_5=utr_score5,
        utr_score_3=utr_score3,
        decision_trail=_decision_trail,
        convergence_status=_convergence_status,
        iterations_used=_total_cai_recovery_iterations + 1,  # +1 for the initial optimization pass
        objective_score=_objective_score,
        biosecurity_screening_result=_biosecurity_result,
    )

    # ── Strict mode: refuse sequences with failed predicates ─────────
    if strict_mode and result.failed_predicates:
        raise OptimizationConstraintError(
            failed_predicates=result.failed_predicates,
            partial_result=result,
        )

    return result


def batch_optimize(
    proteins: list[str],
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    enzymes: list | None = None,
    cai_threshold: float = 0.5,
    strategy: str = "hybrid",
    full_postprocessing: bool = False,
    skip_biosecurity_check: bool = False,
    **kwargs,
) -> list[OptimizationResult]:
    """Optimize multiple proteins in batch.

    Uses a single HybridOptimizer instance for efficiency
    (precomputed data structures are reused across proteins).
    This avoids redundant initialization of sorted_codons, gt_free,
    ag_free, codon_gc, optimal_codon, _max_adapt, and _rs_sites
    for each protein — a significant speedup when optimizing many
    proteins for the same organism.

    When ``full_postprocessing=True``, each result is processed through
    the full :func:`optimize_sequence` pipeline (mRNA stability, codon
    pair bias, CAI recovery, UTR suggestions, provenance tracking).
    When ``False`` (default), a streamlined path is used that skips
    expensive post-processing steps, making batch optimization
    significantly faster for high-throughput workflows.

    Args:
        proteins: List of amino acid sequences (single-letter codes, no stop).
        organism: Target organism. Accepts canonical binomials, short keys,
            abbreviated binomials, or display names (same as
            :func:`optimize_sequence`).
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        enzymes: List of restriction enzyme names to avoid.
        cai_threshold: Minimum CAI score for the CodonAdapted predicate.
        strategy: Optimization strategy (default 'hybrid').
        full_postprocessing: If True, run full optimize_sequence per protein
            (slower but includes mRNA stability, CPB, UTR, provenance).
            If False (default), use streamlined batch path.
        **kwargs: Additional arguments passed to HybridOptimizer or
            optimize_sequence (depending on full_postprocessing).

    Returns:
        List of OptimizationResult, one per input protein, in the same order.

    Raises:
        InvalidProteinError: if any protein contains invalid amino acid codes.
        UnsupportedOrganismError: if the organism is not supported.

    Examples:
        >>> from biocompiler.optimization import batch_optimize
        >>> proteins = ['MSKGEELFTG', 'MALWMRLLPL', 'MVHLTPEEKS']
        >>> results = batch_optimize(proteins, organism='Escherichia_coli')
        >>> assert len(results) == 3
        >>> for r in results:
        ...     assert r.cai > 0.5
    """
    if not proteins:
        return []

    # Resolve organism name once
    organism = resolve_organism(organism, strict=False)

    # Validate all proteins upfront
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    for protein in proteins:
        protein_upper = protein.strip().upper()
        invalid = set(protein_upper) - valid_aas
        if invalid:
            raise InvalidProteinError(protein, invalid)

    # ── Full post-processing path ──────────────────────────────────
    # Delegate to optimize_sequence per protein (no HybridOptimizer
    # reuse, but full feature parity).  This is the safe fallback.
    if full_postprocessing:
        results: list[OptimizationResult] = []
        for protein in proteins:
            result = optimize_sequence(
                target_protein=protein,
                organism=organism,
                gc_lo=gc_lo,
                gc_hi=gc_hi,
                cai_threshold=cai_threshold,
                enzymes=enzymes,
                strategy=strategy,
                skip_biosecurity_check=skip_biosecurity_check,
                **kwargs,
            )
            results.append(result)
        return results

    # ── Streamlined batch path ─────────────────────────────────────
    # Reuse a single HybridOptimizer for all proteins.  The __init__
    # precomputes sorted_codons, gt_free, ag_free, codon_gc,
    # optimal_codon, _max_adapt, _rs_sites, etc. — all of which are
    # identical for the same organism/enzymes/GC parameters.  Reusing
    # them avoids O(20 * AA_count) redundant sort/filter operations
    # per protein.
    from ..hybrid_optimizer import HybridOptimizer as _HybridOptimizer

    species_key = _organism_to_species_key(organism)

    # Detect organism domain
    organism_domain = kwargs.get("organism_domain", "auto")
    if organism_domain not in ("auto", "eukaryote", "prokaryote"):
        organism_domain = "auto"
    if organism_domain == "auto":
        from ..organism_config import is_eukaryotic_organism
        is_prokaryote = not is_eukaryotic_organism(organism)
    elif organism_domain == "prokaryote":
        is_prokaryote = True
    else:
        is_prokaryote = False

    effective_avoid_gt = kwargs.get("avoid_gt", not is_prokaryote)
    effective_splice_low = kwargs.get("splice_low", 3.0 if not is_prokaryote else 999.0)

    # Create ONE HybridOptimizer — precomputed data structures are reused
    hybrid_opt = _HybridOptimizer(
        species=species_key,
        organism=organism,
        enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        avoid_gt=effective_avoid_gt,
        splice_threshold=effective_splice_low,
    )

    results: list[OptimizationResult] = []
    for protein in proteins:
        protein_upper = protein.strip().upper()

        # ── Biosecurity screening ────────────────────────────────────
        _biosecurity_result = None
        if not skip_biosecurity_check:
            from ..biosecurity import check_biosecurity_before_optimize as _check_biosecurity
            _biosecurity_result = _check_biosecurity(
                protein_upper,
                organism=organism,
            )

        # Run hybrid optimization (reuses precomputed data structures)
        hybrid_result = hybrid_opt.optimize(
            protein_upper, is_prokaryote=is_prokaryote
        )
        optimized_seq = hybrid_result.sequence

        # ── Lightweight predicate evaluation ───────────────────────
        pred_results = []

        if is_prokaryote:
            # Prokaryote fast path — skip eukaryotic predicates
            from ..type_system import (
                check_no_stop_codons, check_valid_coding_seq,
                check_no_restriction_site,
            )

            pred_results.append(check_no_stop_codons(optimized_seq))
            pred_results.append(PredicateResult(
                "NoCrypticSplice", True,
                details="Skipped for prokaryotic organism"
            ))
            pred_results.append(PredicateResult(
                "NoCpGIsland", True,
                details="Skipped for prokaryotic organism"
            ))
            pred_results.append(check_no_restriction_site(
                optimized_seq, enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
            ))
            pred_results.append(PredicateResult(
                "NoGTDinucleotide", True,
                details="Skipped for prokaryotic organism"
            ))
            pred_results.append(check_valid_coding_seq(optimized_seq))
            pred_results.append(PredicateResult(
                "ConservationScore", True,
                details="All AA conservation scores >= 0"
            ))

            cai_val = hybrid_result.cai
            all_optimal = cai_val >= cai_threshold
            pred_results.append(PredicateResult(
                "CodonOptimality", all_optimal,
                details=f"CAI={cai_val:.4f}, min={cai_threshold}"
            ))

            _gc_val = hybrid_result.gc_content
            gc_ok = gc_lo <= _gc_val <= gc_hi
            pred_results.append(PredicateResult(
                "GCInRange", gc_ok,
                details=f"GC content: {_gc_val:.3f} (range [{gc_lo}, {gc_hi}])"
            ))
        else:
            # Eukaryote path — use BioOptimizer for predicate evaluation
            opt = BioOptimizer(
                species=species_key,
                enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                splice_low=effective_splice_low,
                splice_high=kwargs.get("splice_high", 6.0),
                min_cai=cai_threshold,
                strategy="constraint_first",
                avoid_gt=effective_avoid_gt,
                organism_name=organism,
                organism_domain="eukaryote",
            )
            opt._applied_mutagenesis = hybrid_opt._applied_mutagenesis
            opt._original_protein = protein_upper
            # Skip redundant MaxEntScan splice check when the hybrid
            # optimizer already validated splice sites during Phase 3.
            _skip_splice = hybrid_result.splice_sites_validated
            pred_results = opt._evaluate_all_predicates(optimized_seq, skip_splice_check=_skip_splice)

        # Build certificate and collect predicate names
        cert_text = format_certificate(pred_results, optimized_seq, species_key)
        failed = [r.predicate for r in pred_results if not r.passed]
        satisfied = [r.predicate for r in pred_results if r.passed]
        fallback = bool(hybrid_opt._applied_mutagenesis)

        # Compute final metrics
        gc = hybrid_result.gc_content
        cai_val = hybrid_result.cai

        # Build OptimizationResult
        result = OptimizationResult(
            sequence=optimized_seq,
            gc_content=gc,
            cai=cai_val,
            failed_predicates=failed,
            predicate_results=pred_results,
            certificate_text=cert_text,
            protein=protein_upper,
            fallback_used=fallback,
            satisfied_predicates=satisfied,
            biosecurity_screening_result=_biosecurity_result,
        )
        results.append(result)

        # Reset per-protein mutagenesis state for next protein
        hybrid_opt._applied_mutagenesis = []

    return results


def _back_translate_protein(
    protein: str,
    species_key: str,
    strategy: str = "greedy",
    enzymes: list[str] | None = None,
    is_eukaryote: bool = True,
) -> str:
    """Back-translate a protein to DNA using highest-CAI codons.

    Uses CODON_ADAPTIVENESS_TABLES (the same table used by compute_cai) so
    the initial sequence starts with optimal codons that will be reflected
    in the final CAI score.

    When ``strategy='dp'`` (or automatically for sequences < 200 aa), a
    dynamic-programming approach is used that considers the top 3 codons
    per position and finds the globally optimal sequence that:
      - Maximizes CAI
      - Avoids common restriction enzyme recognition sites
      - Minimizes GT/AG dinucleotides (for eukaryotes)

    Args:
        protein: Amino acid sequence (1-letter codes).
        species_key: Short species key (e.g. 'ecoli', 'human').
        strategy: 'greedy' for per-position best-CAI; 'dp' for global
            optimisation with cross-codon effects.
        enzymes: Restriction enzyme names to avoid (only used for DP).
        is_eukaryote: If True, penalise GT/AG dinucleotides in DP.

    Returns:
        DNA sequence string.
    """
    organism = _species_key_to_organism(species_key)
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism)
    if adaptiveness is None:
        # Fallback to Homo_sapiens if organism not found
        adaptiveness = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]

    # Decide whether to use DP
    # Explicit 'dp' strategy always uses DP; 'greedy' always skips DP.
    # For other strategies (e.g. 'constraint_first'), use DP for short
    # sequences where it's fast and produces a better starting point.
    if strategy == "dp":
        use_dp = True
    elif strategy == "greedy":
        use_dp = False
    else:
        use_dp = len(protein) < 200
    if use_dp:
        return _back_translate_protein_dp(
            protein, adaptiveness, enzymes=enzymes, is_eukaryote=is_eukaryote,
        )

    # Greedy: simply pick the highest-CAI codon per position
    codons = []
    for aa in protein:
        if aa == "*":
            codons.append("TAA")
            continue
        candidates = AA_TO_CODONS.get(aa, [])
        if not candidates:
            codons.append("NNN")
            continue
        best = max(candidates, key=lambda c: adaptiveness.get(c, 0.0))
        codons.append(best)
    return "".join(codons)


# ────────────────────────────────────────────────────────────
# DP-based back-translation
# ────────────────────────────────────────────────────────────

# Maximum number of top codons to consider per amino-acid position in DP
_DP_TOP_K: int = 3
# Penalty multiplier for GT/AG dinucleotides (eukaryotes)
_GT_AG_PENALTY: float = 0.02
# Penalty for creating a restriction site (effectively -inf)
_RESTRICTION_PENALTY: float = 0.5


def _build_restriction_site_set(
    enzymes: list[str] | None,
) -> set[str]:
    """Build the set of restriction site sequences (fwd + rc) to avoid."""
    sites: set[str] = set()
    if enzymes is None:
        enzymes = list(RESTRICTION_ENZYMES.keys())
    for name in enzymes:
        seq = RESTRICTION_ENZYMES.get(name)
        if seq is None:
            continue
        # Skip IUPAC wildcard sites (contain N or other ambiguity codes)
        if any(b not in "ATCG" for b in seq.upper()):
            continue
        seq_upper = seq.upper()
        sites.add(seq_upper)
        rc = reverse_complement(seq_upper)
        if rc != seq_upper:
            sites.add(rc)
    return sites


def _contains_restriction_site(
    seq: str, sites: set[str], window: int
) -> bool:
    """Check whether *seq* contains any of the restriction *sites*.

    Only checks the last *window* characters (enough to catch any new site
    introduced by the most recent codon).
    """
    start = max(0, len(seq) - window)
    tail = seq[start:]
    for site in sites:
        if site in tail:
            return True
    return False


def _count_dinucleotides(seq: str, di: str) -> int:
    """Count occurrences of dinucleotide *di* in *seq*.

    Uses the NUMBA ``fast_dinucleotide_count`` kernel when available
    for single-pass counting; falls back to pure-Python str.find
    otherwise.
    """
    return _count_dinucs_fast(seq, di)[0]


def _back_translate_protein_dp(
    protein: str,
    adaptiveness: dict[str, float],
    enzymes: list[str] | None = None,
    is_eukaryote: bool = True,
) -> str:
    """DP-based back-translation that considers cross-codon effects.

    For each amino-acid position, the top K codons by CAI are considered.
    The DP state tracks the last two codons (6 nt) to evaluate:
      - GT / AG dinucleotide penalties (eukaryotes)
      - Restriction site avoidance

    The objective is to maximise the sum of log-adaptiveness values
    (equivalent to maximising the geometric mean = CAI) while penalising
    undesirable features.

    Complexity: O(N * K^3) where N = protein length, K = _DP_TOP_K.
    For N < 200 and K = 3 this is negligible (< 1 ms).
    """
    import math

    K = _DP_TOP_K
    restriction_sites = _build_restriction_site_set(enzymes)
    # Maximum restriction site length — only need to check this many
    # trailing characters for newly introduced sites
    max_site_len = max((len(s) for s in restriction_sites), default=0)
    # Window of trailing sequence to check for restriction sites:
    # last codon (3) + overlap from previous codon (max_site_len - 1)
    rest_check_window = max_site_len + 2  # conservative

    # Pre-compute top-K codons per amino acid
    top_codons_per_aa: list[list[tuple[str, float]]] = []
    for aa in protein:
        if aa == "*":
            # Stop codon — only one choice needed
            top_codons_per_aa.append([("TAA", 1.0)])
            continue
        candidates = AA_TO_CODONS.get(aa, [])
        if not candidates:
            top_codons_per_aa.append([("NNN", 0.0)])
            continue
        # Sort by adaptiveness descending, take top K
        sorted_cands = sorted(
            candidates,
            key=lambda c: adaptiveness.get(c, 0.0),
            reverse=True,
        )
        top = [(c, adaptiveness.get(c, 0.0)) for c in sorted_cands[:K]]
        top_codons_per_aa.append(top)

    n = len(protein)

    # DP state: (prev_codon_idx, curr_codon_idx) -> best log-CAI score
    # We track the last two codon choices to evaluate cross-codon effects.
    # prev_codon_idx refers to position i-2, curr_codon_idx to position i-1.
    # At position i, we extend to a new codon choice.

    # For position 0: no previous context
    # dp[(prev_idx, curr_idx)] = (score, backpointer_list)
    # We use a flat dict for the DP frontier.

    # Special handling for the first few positions where we don't have
    # full two-codon context yet.

    # Position 0: just pick a codon
    dp: dict[tuple[int, int], tuple[float, list[int]]] = {}
    for ci, (codon_i, adapt_i) in enumerate(top_codons_per_aa[0]):
        score = math.log(adapt_i) if adapt_i > 0 else -20.0
        dp[(0, ci)] = (score, [ci])

    if n == 1:
        # Trivial case
        best_state = max(dp, key=lambda k: dp[k][0])
        best_codon = top_codons_per_aa[0][dp[best_state][1][0]][0]
        return best_codon

    # Position 1: extend from position 0
    new_dp: dict[tuple[int, int], tuple[float, list[int]]] = {}
    for (prev_ci_key, curr_ci_key), (score, path) in dp.items():
        curr_ci = path[0]  # index into top_codons_per_aa[0]
        curr_codon = top_codons_per_aa[0][curr_ci][0]
        for ci1, (codon1, adapt1) in enumerate(top_codons_per_aa[1]):
            s = score + (math.log(adapt1) if adapt1 > 0 else -20.0)
            # Check cross-codon effects between position 0 and 1
            junction = curr_codon + codon1  # 6 nt
            if is_eukaryote:
                # Penalise GT and AG at the cross-codon boundary
                if curr_codon[-1] + codon1[0] == "GT":
                    s -= _GT_AG_PENALTY * 50  # scale penalty to log space
                if curr_codon[-1] + codon1[0] == "AG":
                    s -= _GT_AG_PENALTY * 50
            # Check restriction sites in the junction
            if restriction_sites:
                for site in restriction_sites:
                    if len(site) <= len(junction) and site in junction:
                        s -= _RESTRICTION_PENALTY
            state = (curr_ci, ci1)
            if state not in new_dp or s > new_dp[state][0]:
                new_dp[state] = (s, path + [ci1])
    dp = new_dp

    # Positions 2..n-1: full DP with two-codon lookback
    for pos in range(2, n):
        new_dp = {}
        for (prev_ci, curr_ci), (score, path) in dp.items():
            prev_codon = top_codons_per_aa[pos - 2][prev_ci][0] if pos >= 2 else ""
            curr_codon = top_codons_per_aa[pos - 1][curr_ci][0]
            for ci, (codon, adapt) in enumerate(top_codons_per_aa[pos]):
                s = score + (math.log(adapt) if adapt > 0 else -20.0)

                # Cross-codon GT/AG penalty (curr_codon | codon boundary)
                if is_eukaryote:
                    boundary = curr_codon[-1] + codon[0]
                    if boundary == "GT":
                        s -= _GT_AG_PENALTY * 50
                    elif boundary == "AG":
                        s -= _GT_AG_PENALTY * 50

                # Restriction site check: only need to check the region
                # that could contain a new site introduced by *codon*.
                # That's the tail of (prev_codon + curr_codon + codon).
                if restriction_sites:
                    # Build enough context to catch cross-codon restriction sites
                    check_seq = curr_codon + codon
                    if prev_codon:
                        check_seq = prev_codon[-1] + check_seq
                    found_restriction = False
                    for site in restriction_sites:
                        if site in check_seq:
                            found_restriction = True
                            break
                    if found_restriction:
                        s -= _RESTRICTION_PENALTY

                state = (curr_ci, ci)
                if state not in new_dp or s > new_dp[state][0]:
                    new_dp[state] = (s, path + [ci])
        dp = new_dp

    # Extract best path
    if not dp:
        # Fallback: greedy
        codons = []
        for aa in protein:
            if aa == "*":
                codons.append("TAA")
                continue
            candidates = AA_TO_CODONS.get(aa, [])
            if not candidates:
                codons.append("NNN")
                continue
            best = max(candidates, key=lambda c: adaptiveness.get(c, 0.0))
            codons.append(best)
        return "".join(codons)

    best_state = max(dp, key=lambda k: dp[k][0])
    best_path = dp[best_state][1]

    # Reconstruct sequence
    result_codons = []
    for pos, idx in enumerate(best_path):
        result_codons.append(top_codons_per_aa[pos][idx][0])
    return "".join(result_codons)


class BioOptimizer:
    """Certified gene sequence optimizer with multi-step CAI-maximizing pipeline."""

    # NOTE (Task 1.8): _SPECIES_TO_ORGANISM removed.  Name resolution
    # now uses the canonical resolve_organism() from organisms/__init__.py
    # which has 30+ aliases instead of the previous 10-entry partial copy.

    def __init__(
        self,
        species: str = "ecoli",
        enzymes: Optional[List[str]] = None,
        splice_low: float = 3.0,
        splice_high: float = 6.0,
        cpg_window: int = 200,
        cpg_threshold: float = 0.6,
        min_blosum: int = -1,
        min_cai: float = 0.0,
        avoid_gt: bool = True,
        strategy: str = "constraint_first",
        optimize_mrna_stability: bool = True,
        biosecurity_mode: Optional[str] = None,
        check_g4: bool = True,
        check_mrna_halflife: bool = True,
        halflife_threshold: float = 2.0,
        use_lineardesign: bool = False,
        lineardesign_lambda: float = 3.0,
        tasep_ensemble: bool = True,
        tasep_n_runs: int = 100,
        check_nmd: bool = True,
        **kwargs: Any,
    ) -> None:
        self.species = species
        self.biosecurity_mode: Optional[str] = biosecurity_mode

        # Organism-aware attributes (must be set before species_cai so that
        # we use the same CODON_ADAPTIVENESS_TABLES as compute_cai())
        from ..organism_config import is_eukaryotic_organism

        # Resolve organism_name to its canonical form so that it can be
        # used as a key into CODON_ADAPTIVENESS_TABLES.  Accepts both
        # short aliases ("ecoli") and full names ("Escherichia_coli").
        # Uses the canonical resolve_organism() from organisms/__init__.py
        # instead of a partial local _SPECIES_TO_ORGANISM mapping.
        raw_organism = kwargs.get(
            "organism_name",
            resolve_organism(species),
        )
        self.organism_name: str = resolve_organism(raw_organism)

        # Use CODON_ADAPTIVENESS_TABLES for codon selection so that
        # the optimizer and compute_cai() agree on which codons are optimal.
        # Previously used get_species_cai_weights() which pulled from a
        # different table (SPECIES['ecoli']['cai_weights'] derived from
        # ECOLI_CODON_USAGE) and could pick non-optimal codons.
        self.species_cai: Dict[str, float] = CODON_ADAPTIVENESS_TABLES.get(
            self.organism_name, CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
        )

        self.enzymes: List[str] = enzymes or []
        self.splice_low: float = splice_low
        self.splice_high: float = splice_high
        self.cpg_window: int = cpg_window
        self.cpg_threshold: float = cpg_threshold
        self.min_blosum: int = min_blosum
        self.min_cai: float = min_cai
        self.avoid_gt: bool = avoid_gt
        self.strategy: str = strategy  # "constraint_first" or "cai_first"
        self.optimize_mrna_stability: bool = optimize_mrna_stability
        self.check_g4: bool = check_g4
        self.check_mrna_halflife: bool = check_mrna_halflife
        self.halflife_threshold: float = halflife_threshold
        self.use_lineardesign: bool = use_lineardesign
        self.lineardesign_lambda: float = lineardesign_lambda
        self.tasep_ensemble: bool = tasep_ensemble
        self.tasep_n_runs: int = tasep_n_runs
        self.check_nmd: bool = check_nmd

        # Sliding-window GC constraint parameters
        self._gc_window_size: int = kwargs.get("gc_window_size", 50)
        self._gc_window_min: float | None = kwargs.get("gc_window_min", None)
        self._gc_window_max: float | None = kwargs.get("gc_window_max", None)

        self.organism_domain: str = kwargs.get("organism_domain", "auto")
        if self.organism_domain == "auto":
            self.organism_domain = (
                "eukaryote" if is_eukaryotic_organism(self.organism_name) else "prokaryote"
            )

        # Derived flag: True when target organism is prokaryotic.
        # Controls whether eukaryote-specific constraint steps (cryptic splice
        # elimination, CpG disruption, GT resolution, cross-codon GT/CG
        # coordination) are skipped during optimization.  Prokaryotes have no
        # spliceosome, so GT/AG avoidance and CpG island disruption are
        # biologically irrelevant and unnecessarily lower CAI.
        self.is_prokaryote: bool = self.organism_domain == "prokaryote"

        # Auto-set avoid_gt for prokaryotes unless explicitly overridden
        if "avoid_gt" not in kwargs and self.is_prokaryote:
            self.avoid_gt = False
        # Track positions where GT is unavoidable (e.g., Valine codons)
        self._unavoidable_gt_positions: Set[int] = set()
        # Track mutagenesis proposals that were applied
        self._applied_mutagenesis: List[Dict[str, Any]] = []
        # Store original input protein for conservation scoring
        self._original_protein: str = ""
        # Track mRNA stability metrics from the last optimization run
        self._mrna_stability_score: float | None = None
        self._destabilizing_motifs_removed: int = 0
        self._stability_improvement: float | None = None
        # Decision provenance collector (set externally for codon-level tracking)
        self._provenance_collector: DecisionProvenanceCollector | None = None
        # Protein encoding verification (True when last run verified OK)
        self.verify_encoding: bool = kwargs.get("verify_encoding", True)
        self._encoding_verified: bool | None = None

    def _make_config(self) -> _OptConfig:
        """Build an _OptConfig bundle from the current optimizer state.

        Used to pass optimizer configuration to extracted sub-module
        functions (postprocessing, validation) without coupling them
        to the BioOptimizer class.
        """
        return _OptConfig(
            species_cai=self.species_cai,
            enzymes=self.enzymes,
            cpg_window=self.cpg_window,
            cpg_threshold=self.cpg_threshold,
            avoid_gt=self.avoid_gt,
            min_cai=self.min_cai,
            min_blosum=self.min_blosum,
            splice_low=self.splice_low,
            splice_high=self.splice_high,
            is_prokaryote=self.is_prokaryote,
            organism_domain=self.organism_domain,
            organism_name=self.organism_name,
            species=self.species,
            optimize_mrna_stability=self.optimize_mrna_stability,
            gc_window_size=self._gc_window_size,
            gc_window_min=self._gc_window_min,
            gc_window_max=self._gc_window_max,
            original_protein=self._original_protein,
            applied_mutagenesis=self._applied_mutagenesis,
        )

    def optimize(self, seq: str, strategy: Optional[str] = None) -> Tuple[str, List[PredicateResult], str]:
        """Run the full optimization pipeline.

        Args:
            seq: Input DNA or protein sequence.
            strategy: Optimization strategy override.
                - "constraint_first" (default): GT-aware greedy then fix constraints
                - "cai_first": Maximize CAI first, then fix constraints with
                  minimal CAI impact (DNAworks-style)

        Returns:
            (optimized_sequence, predicate_results, certificate_text)
        """
        effective_strategy = strategy if strategy is not None else self.strategy

        seq = seq.upper().strip()

        # ── Resolve IUPAC ambiguous bases ──────────────────────────
        # If the input contains IUPAC ambiguity codes (R, Y, S, W, K, M,
        # B, D, H, V), resolve them to concrete bases before optimization.
        # This allows users to pass degenerate sequences (e.g., from
        # consensus sequences or degenerate primers) directly to the
        # optimizer.
        ambiguous_chars = set(seq) - set("ACGT")
        if ambiguous_chars:
            from ..iupac import has_ambiguous, resolve_ambiguous
            if has_ambiguous(seq):
                logger.info(
                    "Resolving %d IUPAC ambiguous bases in input sequence "
                    "(strategy: most_common)",
                    sum(1 for b in seq if b not in "ACGT"),
                )
                seq = resolve_ambiguous(
                    seq,
                    strategy="most_common",
                    cai_table=self.species_cai,
                )

        self._unavoidable_gt_positions = set()
        self._applied_mutagenesis = []
        self._original_protein = self._translate(seq)
        # Reset mRNA stability metrics
        self._mrna_stability_score = None
        self._destabilizing_motifs_removed = 0
        self._stability_improvement = None

        # ── Biosecurity screening ──────────────────────────────────
        # Screen the input protein for known hazardous signatures
        # BEFORE any optimization begins.
        from ..biosecurity import check_biosecurity_before_optimize
        try:
            check_biosecurity_before_optimize(
                self._original_protein,
                organism=self.organism_name,
                biosecurity_mode=self.biosecurity_mode,
            )
        except BiosecurityError:
            raise

        if effective_strategy == "cai_first":
            return self._optimize_cai_first(seq)

        # ── Default: constraint_first strategy ──
        # Step: Backtranslate CAI (DNAworks-style)
        seq = self._step_backtranslate_cai(seq)

        # Step: Resolve Constraints (fix GT/CG/RS with minimal CAI loss)
        # Skipped for prokaryotic targets when avoid_gt=False
        seq = self._step_resolve_constraints(seq)

        # Step: Remove Restriction Sites
        seq = self._step_remove_restriction_sites(seq)

        # Step: Cross-Codon Optimization (iterative)
        seq, mut_report = self._step_cross_codon_optimization(seq)

        # Step: Within-Codon GT Resolution
        # Skipped for prokaryotic targets (no spliceosome → GT avoidance unnecessary)
        if self.avoid_gt:
            seq, mut_report_35 = self._step_within_codon_gt_resolution(seq)
            mut_report.proposals.extend(mut_report_35.proposals)

        # Step: Mutagenesis Fallback (aggressive, handles within-codon GTs too)
        # Only needed when GT avoidance is active
        if self.avoid_gt:
            seq = self._step_mutagenesis_fallback(seq, mut_report)

        # Step: Avoid CpG Islands
        # Skipped for prokaryotic targets (CpG islands are a eukaryotic gene regulation concern)
        if not self.is_prokaryote:
            seq = self._step_avoid_cpg_islands(seq)

        # Step: Cross-Codon Coordination (handles cross-codon GT, CG, restriction sites)
        # GT/CG coordination is skipped for prokaryotes; restriction site coordination always runs
        seq = self._step_cross_codon_coordination(seq)

        # Step: CAI Hill Climb (upgrade codons while maintaining constraints)
        seq = self._step_cai_hill_climb(seq)

        # Step: Reoptimize (iterative until convergence)
        seq = self._step_reoptimize(seq)

        # Step: Remove ATTTA Instability Motifs
        seq = self._step_remove_instability_motifs(seq)

        # Step: CpG Reconciliation (aggressive, after CAI hill climb/reoptimize)
        # Skipped for prokaryotic targets
        if not self.is_prokaryote:
            seq = self._step_cpg_reconciliation(seq)

        # Step: GT Reconciliation (fix avoidable GTs that may have been introduced)
        # Skipped for prokaryotic targets (GT is not a constraint)
        if self.avoid_gt:
            seq = self._step_gt_reconciliation(seq)

        # Step: mRNA Stability Improvement (soft optimization — remove destabilizing motifs)
        seq = self._step_mrna_stability_improvement(seq)

        # Step: G-quadruplex check (detect and fix G4 motifs)
        if self.check_g4:
            seq = self._step_g4_check(seq)

        # Step: LinearDesign MFE optimization (when requested)
        if self.use_lineardesign:
            seq = self._step_lineardesign(seq)

        # Step: Ensemble TASEP ribosome simulation
        if self.tasep_ensemble:
            self._step_tasep_ensemble(seq)

        # Step: NMD (nonsense-mediated decay) check
        if self.check_nmd and not self.is_prokaryote:
            self._step_nmd_check(seq)

        # Evaluate all 12 predicates
        results = self._evaluate_all_predicates(seq)

        # Generate certificate
        cert_text = format_certificate(results, seq, self.species)

        return seq, results, cert_text

    # ──────────────────────────────────────────────────────────
    # CAI-first optimization strategy (DNAworks-style)
    # ──────────────────────────────────────────────────────────
    def _optimize_cai_first(self, seq: str) -> Tuple[str, List[PredicateResult], str]:
        """CAI-first optimization: maximize CAI first, fix constraints second.

        Strategy (DNAworks-style):
        1. Back-translate using highest-CAI codons at every position
           (ignore GT/CG/RS initially)
        2. Iteratively fix constraint violations with minimal CAI impact:
           a. Restriction sites → fix each with best synonymous codon
           b. Avoidable GT dinucleotides → fix with best non-GT synonymous codon
              using cross-codon pair optimization
           c. CpG islands → fix with best non-CG synonymous codon
           d. Cryptic splice sites → fix with best synonymous codon
        3. Mutagenesis fallback for GTs that can't be resolved by
           synonymous substitution (e.g., Valine V→Isoleucine I)
        4. CAI hill climbing to recover any lost CAI while maintaining constraints
        5. Iterate until convergence

        Key insight: by starting from max CAI (CAI=1.0) and only making the
        smallest necessary CAI sacrifices to fix constraints, we achieve much
        higher CAI than the constraint_first strategy which permanently sacrifices
        CAI by avoiding GT during codon selection.

        For GT fixing, we use a priority-based approach that fixes GTs with
        the lowest CAI cost first, and considers multi-codon windows to find
        globally optimal fixes.
        """
        import math

        # Step: Maximize CAI — Pure max-CAI back-translation (ignore ALL constraints)
        protein = self._translate(seq)
        codons_result = []
        for i, aa in enumerate(protein):
            if aa == "*":
                codon_start = i * 3
                codons_result.append(
                    seq[codon_start:codon_start + 3]
                    if codon_start + 3 <= len(seq) else "TAA"
                )
                continue
            candidates = AA_TO_CODONS.get(aa, [])
            if candidates:
                codons_result.append(
                    max(candidates, key=lambda c: self.species_cai.get(c, 0.0))
                )
            else:
                codon_start = i * 3
                codons_result.append(seq[codon_start:codon_start + 3])
        seq = "".join(codons_result)

        # Step: Fix restriction sites (highest priority — binary constraint)
        seq = self._cai_first_fix_restriction_sites(seq)

        # Step: Fix avoidable GT dinucleotides (minimal CAI impact)
        # EUKARYOTE-ONLY: GT avoidance is irrelevant for prokaryotes
        if not self.is_prokaryote:
            seq = self._cai_first_fix_gts(seq)

        # Step: Mutagenesis fallback for Valine GTs
        # EUKARYOTE-ONLY: Only needed when GT avoidance is active
        if not self.is_prokaryote:
            seq = self._cai_first_mutagenesis_fallback(seq)

        # Step: Fix CpG islands (minimal CAI impact)
        # EUKARYOTE-ONLY: CpG islands are a eukaryotic gene regulation concern
        if not self.is_prokaryote:
            seq = self._cai_first_fix_cpg(seq)

        # Step: Fix cryptic splice sites (minimal CAI impact)
        # EUKARYOTE-ONLY: Splice sites are irrelevant for prokaryotes
        if not self.is_prokaryote:
            seq = self._cai_first_fix_splice(seq)

        # Step: Cross-Codon Coordination (handles cross-codon GT, CG, restriction sites)
        seq = self._step_cross_codon_coordination(seq)

        # Step: CAI hill climbing (upgrade codons while maintaining constraints)
        seq = self._step_cai_hill_climb(seq)

        # Step: Aggressive re-optimization pass
        seq = self._step_reoptimize(seq)

        # Step: Second pass of GT fixing + CAI boost (iterative refinement)
        # EUKARYOTE-ONLY: GT fixing not needed for prokaryotes
        if not self.is_prokaryote:
            for _refinement in range(3):
                old_cai = self._compute_seq_cai(seq)
                seq = self._cai_first_fix_gts(seq)
                seq = self._step_cai_hill_climb(seq)
                seq = self._step_reoptimize(seq)
                new_cai = self._compute_seq_cai(seq)
                if new_cai <= old_cai + 0.0001:
                    break

        # Step: CpG Reconciliation (aggressive, after CAI hill climb/reoptimize)
        # EUKARYOTE-ONLY: CpG reconciliation not needed for prokaryotes
        if not self.is_prokaryote:
            seq = self._step_cpg_reconciliation(seq)

        # Step: CAI Reconciliation (upgrade low-CAI codons while maintaining hard constraints)
        seq = self._step_cai_reconciliation(seq)

        # Track unavoidable GTs for certificate
        for i in range(len(seq) - 1):
            if seq[i] == "G" and seq[i + 1] == "T":
                if _is_unavoidable_gt(seq, i):
                    codon_start = (i // 3) * 3
                    self._unavoidable_gt_positions.add(i)

        # Step: mRNA Stability Improvement (soft optimization)
        seq = self._step_mrna_stability_improvement(seq)

        # Step: G-quadruplex check (detect and fix G4 motifs)
        if self.check_g4:
            seq = self._step_g4_check(seq)

        # Step: LinearDesign MFE optimization (when requested)
        if self.use_lineardesign:
            seq = self._step_lineardesign(seq)

        # Step: Ensemble TASEP ribosome simulation
        if self.tasep_ensemble:
            self._step_tasep_ensemble(seq)

        # Step: NMD (nonsense-mediated decay) check
        if self.check_nmd and not self.is_prokaryote:
            self._step_nmd_check(seq)

        # Evaluate all predicates
        results = self._evaluate_all_predicates(seq)
        cert_text = format_certificate(results, seq, self.species)
        return seq, results, cert_text

    def _compute_seq_cai(self, seq: str) -> float:
        """Compute the geometric mean CAI for a sequence.

        Uses the NUMBA-accelerated compute_cai_kernel when available,
        falling back to the pure-Python loop otherwise.
        """
        if not seq or len(seq) < 3:
            return 0.0

        # Fast path: NUMBA kernel
        if HAS_NUMBA and _numba_cai_kernel is not None:
            try:
                # Lazy-initialize the per-optimizer adaptiveness array
                if not hasattr(self, '_numba_adapt_arr') or self._numba_adapt_arr is None:
                    self._numba_adapt_arr = _adaptiveness_to_array(self.species_cai)
                adapt_arr = self._numba_adapt_arr
                # Build codon indices, excluding Met and stop codons
                n_total = len(seq) // 3
                idx_list = []
                for i in range(n_total):
                    codon = seq[i * 3:i * 3 + 3]
                    aa = CODON_TABLE.get(codon)
                    if aa == 'M' or aa == '*':
                        continue
                    idx_list.append(_codon_to_index(codon))
                if not idx_list:
                    return 0.0
                codon_indices = _np.array(idx_list, dtype=_np.int64)
                n_codons = len(codon_indices)
                return _numba_cai_kernel(adapt_arr, codon_indices, n_codons)
            except Exception:
                logger.warning("NUMBA CAI kernel failed in BioOptimizer, falling back to pure-Python", exc_info=True)

        # Pure-Python fallback
        import math
        log_sum = 0.0
        count = 0
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            cai = self.species_cai.get(codon, 0.0)
            if cai <= 0:
                cai = 0.001
            log_sum += math.log(cai)
            count += 1
        if count == 0:
            return 0.0
        return math.exp(log_sum / count)

    def _step_maximize_cai(self, seq: str) -> str:
        """Maximize CAI step (cai_first): Back-translate with absolute max CAI everywhere.

        No GT avoidance at all - just pick the highest-CAI codon for each AA.
        Constraint violations will be fixed in subsequent steps.
        """
        protein = self._translate(seq)
        codons_result = []
        for i, aa in enumerate(protein):
            if aa == "*":
                codon_start = i * 3
                codons_result.append(seq[codon_start:codon_start + 3] if codon_start + 3 <= len(seq) else "TAA")
                continue
            candidates = AA_TO_CODONS.get(aa, [])
            if candidates:
                codons_result.append(max(candidates, key=lambda c: self.species_cai.get(c, 0.0)))
            else:
                codon_start = i * 3
                codons_result.append(seq[codon_start:codon_start + 3])
        return "".join(codons_result)

    def _cai_first_fix_restriction_sites(self, seq: str) -> str:
        """CAI-first Fix restriction sites with minimal CAI impact."""
        from ..restriction_sites import get_recognition_site
        seq_list = list(seq)

        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue

            max_rounds = 50
            for _ in range(max_rounds):
                current_seq = "".join(seq_list)
                p = current_seq.find(site)
                if p == -1:
                    break

                # Find all codon positions overlapping this site
                codon_starts = set()
                for j in range(p, p + len(site)):
                    cs = (j // 3) * 3
                    if cs + 3 <= len(seq_list):
                        codon_starts.add(cs)

                fixed = False
                # Try each overlapping codon, sorted by CAI impact (try best-CAI alts first)
                for cs in sorted(codon_starts):
                    codon = "".join(seq_list[cs:cs + 3])
                    aa = CODON_TABLE.get(codon)
                    if aa is None or aa == "*":
                        continue

                    # Sort alternatives by CAI (highest first) to minimize CAI loss
                    alts = sorted(
                        AA_TO_CODONS.get(aa, []),
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    )
                    for alt in alts:
                        if alt == codon:
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(alt):
                            test_list[cs + k] = b
                        test_seq = "".join(test_list)
                        if site not in test_seq:
                            seq_list = test_list
                            fixed = True
                            break
                    if fixed:
                        break

                if not fixed:
                    # Could not fix this site with single-codon substitution
                    # Try two-codon substitution
                    fixed = self._cai_first_fix_rs_two_codons(seq_list, p, site)
                    if not fixed:
                        break  # Give up on this site for now

        return "".join(seq_list)

    def _cai_first_fix_rs_two_codons(self, seq_list: list, pos: int, site: str) -> bool:
        """Try fixing a restriction site by modifying two adjacent codons."""
        codon_starts = sorted(set(
            (j // 3) * 3
            for j in range(pos, min(pos + len(site), len(seq_list)))
            if (j // 3) * 3 + 3 <= len(seq_list)
        ))

        if len(codon_starts) < 2:
            return False

        # Try pairs of codons
        for idx in range(len(codon_starts) - 1):
            cs1, cs2 = codon_starts[idx], codon_starts[idx + 1]
            if cs2 != cs1 + 3:
                continue  # Only adjacent codons

            codon1 = "".join(seq_list[cs1:cs1 + 3])
            codon2 = "".join(seq_list[cs2:cs2 + 3])
            aa1 = CODON_TABLE.get(codon1)
            aa2 = CODON_TABLE.get(codon2)
            if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                continue

            # Try all pairs sorted by combined CAI
            pairs = []
            for c1 in AA_TO_CODONS.get(aa1, [codon1]):
                for c2 in AA_TO_CODONS.get(aa2, [codon2]):
                    combined = self.species_cai.get(c1, 0.0) + self.species_cai.get(c2, 0.0)
                    pairs.append((c1, c2, combined))
            pairs.sort(key=lambda x: x[2], reverse=True)

            for c1, c2, _ in pairs:
                test_list = seq_list[:]
                for k, b in enumerate(c1):
                    test_list[cs1 + k] = b
                for k, b in enumerate(c2):
                    test_list[cs2 + k] = b
                test_seq = "".join(test_list)
                if site not in test_seq:
                    seq_list[:] = test_list
                    return True

        return False

    def _cai_first_fix_gts(self, seq: str) -> str:
        """CAI-first Fix avoidable GT dinucleotides with minimal CAI impact.

        Uses cost-aware GT resolution with splice donor scoring for eukaryotes:
        - GTs with low splice donor potential (< SPLICE_DONOR_POTENTIAL_THRESHOLD)
          are NOT fixed — CAI takes priority over low-risk GTs.
        - GTs with high splice donor potential are fixed with minimal CAI impact.
        - For prokaryotes, this method is never called.

        Iteratively finds each avoidable GT and resolves it by choosing
        the synonymous substitution(s) with the highest possible CAI that
        eliminates the GT, but ONLY if the GT has high splice donor potential.
        """
        seq_list = list(seq)
        max_rounds = 50

        for round_num in range(max_rounds):
            current_seq = "".join(seq_list)
            violations = []

            # Find all avoidable GT positions, filtered by splice donor potential
            for i in range(len(current_seq) - 1):
                if current_seq[i] == "G" and current_seq[i + 1] == "T":
                    if _is_unavoidable_gt(current_seq, i):
                        continue
                    # Cost-aware: skip GTs with low splice donor potential
                    sdp = score_splice_donor_potential(current_seq, i)
                    if sdp < SPLICE_DONOR_POTENTIAL_THRESHOLD:
                        # This GT has low splice donor potential — not dangerous.
                        # Accept it (CAI > GT avoidance for low-risk GTs).
                        continue
                    violations.append(i)

            if not violations:
                break

            any_fixed = False

            for gt_pos in violations:
                codon_start = (gt_pos // 3) * 3
                next_codon_start = codon_start + 3

                # Determine if within-codon or cross-codon
                is_within = (gt_pos + 1) < next_codon_start

                if is_within:
                    # Within-codon GT: try synonymous substitution
                    fixed = self._cai_first_fix_within_gt(seq_list, codon_start)
                else:
                    # Cross-codon GT: try adjacent codon pair substitution
                    fixed = self._cai_first_fix_cross_gt(seq_list, codon_start)

                if fixed:
                    any_fixed = True
                    break  # Restart scanning after any fix

            if not any_fixed:
                break

        return "".join(seq_list)

    def _cai_first_fix_within_gt(self, seq_list: list, codon_start: int) -> bool:
        """Fix a within-codon GT by choosing the highest-CAI alternative without GT."""
        codon = "".join(seq_list[codon_start:codon_start + 3])
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            return False

        old_gt_count = _count_gts("".join(seq_list))

        # Sort alternatives by CAI (highest first), skip codons with GT
        alternatives = []
        for alt in AA_TO_CODONS.get(aa, []):
            if "GT" in alt:
                continue
            alt_cai = self.species_cai.get(alt, 0.0)
            alternatives.append((alt, alt_cai))

        alternatives.sort(key=lambda x: x[1], reverse=True)

        for alt, _ in alternatives:
            test_list = seq_list[:]
            for k, b in enumerate(alt):
                test_list[codon_start + k] = b
            new_gt_count = _count_gts("".join(test_list))
            if new_gt_count < old_gt_count:
                seq_list[:] = test_list
                return True

        return False

    def _cai_first_fix_cross_gt(self, seq_list: list, codon_start: int) -> bool:
        """Fix a cross-codon GT with minimal CAI impact.

        Tries strategies in order of increasing invasiveness:
        D. Change only the next codon (to one that doesn't start with T)
        C. Change only the current codon (to one that doesn't end with G)
        A. Modify the codon pair (both codons) — sorted by combined CAI
        B. Modify preceding codon pair
        """
        old_gt_count = _count_gts("".join(seq_list))
        next_start = codon_start + 3

        # Strategy D: Change only the next codon (cheapest single-codon fix)
        if next_start + 3 <= len(seq_list):
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if aa2 is not None and aa2 != "*":
                alts = sorted(
                    AA_TO_CODONS.get(aa2, [aa2_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c2 in alts:
                    if c2[0] == "T":
                        continue
                    test_list = seq_list[:]
                    for k, b in enumerate(c2):
                        test_list[next_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        # Strategy C: Change only the current codon (to one that doesn't end with G)
        aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
        aa1 = CODON_TABLE.get(aa1_codon)
        if aa1 is not None and aa1 != "*":
            alts = sorted(
                AA_TO_CODONS.get(aa1, [aa1_codon]),
                key=lambda c: self.species_cai.get(c, 0.0),
                reverse=True,
            )
            for c1 in alts:
                if c1[-1] == "G":
                    continue
                test_list = seq_list[:]
                for k, b in enumerate(c1):
                    test_list[codon_start + k] = b
                test_seq = "".join(test_list)
                new_gt_count = _count_gts(test_seq)
                if new_gt_count < old_gt_count:
                    seq_list[:] = test_list
                    return True

        # Strategy A: Modify current + following codon
        if next_start + 3 <= len(seq_list):
            aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
            aa1 = CODON_TABLE.get(aa1_codon)
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)

            if aa1 is not None and aa1 != "*" and aa2 is not None and aa2 != "*":
                pairs = []
                for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                    for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
                        if c1[-1] + c2[0] == "GT":
                            continue
                        combined_cai = self.species_cai.get(c1, 0.0) + self.species_cai.get(c2, 0.0)
                        pairs.append((c1, c2, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c1, c2, _ in pairs:
                    test_list = seq_list[:]
                    for k, b in enumerate(c1):
                        test_list[codon_start + k] = b
                    for k, b in enumerate(c2):
                        test_list[next_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        # Strategy B: Modify preceding + current codon
        if codon_start >= 3:
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
            aa1 = CODON_TABLE.get(aa1_codon)

            if aa0 is not None and aa0 != "*" and aa1 is not None and aa1 != "*":
                pairs = []
                for c0 in AA_TO_CODONS.get(aa0, [aa0_codon]):
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                        if c0[-1] + c1[0] == "GT":
                            continue
                        combined_cai = self.species_cai.get(c0, 0.0) + self.species_cai.get(c1, 0.0)
                        pairs.append((c0, c1, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c0, c1, _ in pairs:
                    test_list = seq_list[:]
                    for k, b in enumerate(c0):
                        test_list[prev_start + k] = b
                    for k, b in enumerate(c1):
                        test_list[codon_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        return False

    def _cai_first_fix_cpg(self, seq: str) -> str:
        """CAI-first Fix CpG islands with minimal CAI impact."""
        seq_list = list(seq)
        changed = True
        iterations = 0
        max_iterations = 50

        while changed and iterations < max_iterations:
            changed = False
            iterations += 1

            for start in range(0, len(seq_list) - self.cpg_window + 1, 3):
                window = "".join(seq_list[start:start + self.cpg_window])
                c_count = window.count("C")
                g_count = window.count("G")
                cg_count = _count_dinucs_fast(window, "CG")[0]
                expected = (c_count * g_count) / len(window) if len(window) > 0 else 0
                obs_exp = cg_count / expected if expected > 0 else 0.0

                if obs_exp <= self.cpg_threshold:
                    continue

                # Find CG dinucleotides and try to break them
                for i in range(start, min(start + self.cpg_window - 1, len(seq_list) - 1)):
                    if seq_list[i] == "C" and seq_list[i+1] == "G":
                        codon_start = (i // 3) * 3
                        if codon_start + 3 > len(seq_list):
                            continue
                        codon = "".join(seq_list[codon_start:codon_start + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue

                        old_gt_count = _count_gts("".join(seq_list))

                        # Sort by CAI descending to minimize CAI loss
                        for alt in sorted(
                            AA_TO_CODONS.get(aa, []),
                            key=lambda c: self.species_cai.get(c, 0.0),
                            reverse=True,
                        ):
                            if alt == codon or "CG" in alt:
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[codon_start + k] = b
                            new_gt_count = _count_gts("".join(test_list))
                            if new_gt_count <= old_gt_count:
                                seq_list = test_list
                                changed = True
                                break

                        if changed:
                            break
                if changed:
                    break

        return "".join(seq_list)

    def _cai_first_fix_splice(self, seq: str) -> str:
        """CAI-first Fix cryptic splice sites with minimal CAI impact."""
        seq_list = list(seq)
        max_rounds = 30

        for _ in range(max_rounds):
            current_seq = "".join(seq_list)
            splice_result = check_no_cryptic_splice(current_seq, self.splice_low, self.splice_high)
            if splice_result.passed:
                break

            # Find and fix splice sites
            # The splice check looks for GT..AG patterns that resemble splice donors/acceptors
            fixed = False
            for i in range(len(current_seq) - 1):
                if current_seq[i] == "G" and current_seq[i + 1] == "T":
                    # This GT might be part of a cryptic splice site
                    codon_start = (i // 3) * 3
                    next_codon_start = codon_start + 3
                    is_within = (i + 1) < next_codon_start

                    if _is_unavoidable_gt(current_seq, i):
                        continue

                    if is_within:
                        fixed = self._cai_first_fix_within_gt(seq_list, codon_start)
                    else:
                        fixed = self._cai_first_fix_cross_gt(seq_list, codon_start)

                    if fixed:
                        break

            if not fixed:
                break

        return "".join(seq_list)

    def _cai_first_mutagenesis_fallback(self, seq: str) -> str:
        """CAI-first Apply mutagenesis for GTs that can't be resolved
        by synonymous substitution.

        Specifically targets Valine codons (GTN) which all contain GT.
        Substitutes V→I (Isoleucine, BLOSUM62=3) using highest-CAI Ile codon
        to eliminate the GT while maximizing CAI.
        """
        seq_list = list(seq)
        changed = True
        max_rounds = 10

        for _ in range(max_rounds):
            if not changed:
                break
            changed = False

            for i in range(0, len(seq_list) - 2, 3):
                codon = "".join(seq_list[i:i+3])
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*":
                    continue

                if "GT" not in codon:
                    continue

                has_gt_free = any("GT" not in c for c in AA_TO_CODONS.get(aa, []))
                if has_gt_free:
                    continue

                if aa == "V":
                    ile_codons = AA_TO_CODONS.get("I", [])
                    old_gt_count = _count_gts("".join(seq_list))

                    for ile_codon in sorted(
                        ile_codons,
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    ):
                        test_list = seq_list[:]
                        for k, b in enumerate(ile_codon):
                            test_list[i + k] = b
                        new_gt_count = _count_gts("".join(test_list))
                        if new_gt_count < old_gt_count:
                            seq_list[:] = test_list
                            self._applied_mutagenesis.append({
                                "position": i,
                                "original_aa": "V",
                                "new_aa": "I",
                                "blosum": BLOSUM62.get(("V", "I"), 3),
                            })
                            changed = True
                            break

        return "".join(seq_list)

    # Deprecated alias — use _step_maximize_cai instead
    _phase0_pure_max_cai = _step_maximize_cai

    # ──────────────────────────────────────────────────────────
    # Step: Backtranslate CAI (DP-based max-CAI back-translation)
    # ──────────────────────────────────────────────────────────
    def _step_backtranslate_cai(self, seq: str) -> str:
        """Backtranslate CAI step: DP-based max-CAI back-translation with avoidable-GT avoidance.

        Uses Viterbi-style dynamic programming to find the globally optimal
        codon assignment that maximizes CAI while avoiding only AVOIDABLE GTs.
        The DP state is simply the last character of the previous codon,
        making it O(n * 4 * k) where n is protein length and k is max codons per AA.

        Key insight: a cross-codon GT (codon1 ends with G, codon2 starts with T)
        is only avoidable if we can change at least one of the two codons to
        eliminate it. If codon2 has NO synonymous codon that doesn't start with T
        (e.g., Trp=TGG, Cys=TGT/TGC, Tyr=TAT/TAC), then the cross-codon GT
        is unavoidable and we should use the highest-CAI codon for codon1.

        This matches the semantics of the check_no_avoidable_gt predicate,
        which only fails on GTs that CAN be avoided by synonymous substitution.

        Only non-Valine within-codon GTs and avoidable cross-codon GTs are
        excluded by the DP. Unavoidable GTs (Valine within-codon, cross-codon
        where next AA has no non-T-starting codon) are allowed.
        """
        import math
        protein = self._translate(seq)
        n = len(protein)
        INF = float('-inf')

        # Precompute which amino acids have at least one codon not starting with T
        # These are the AAs where cross-codon GT from previous codon can be avoided
        aa_has_non_t_start = {}
        for aa_key in set(CODON_TABLE.values()):
            if aa_key == "*":
                continue
            codons_list = AA_TO_CODONS.get(aa_key, [])
            aa_has_non_t_start[aa_key] = any(c[0] != 'T' for c in codons_list)

        # Precompute which amino acids have at least one codon not ending with G
        # These are the AAs where cross-codon GT to next codon can be avoided
        aa_has_non_g_end = {}
        for aa_key in set(CODON_TABLE.values()):
            if aa_key == "*":
                continue
            codons_list = AA_TO_CODONS.get(aa_key, [])
            aa_has_non_g_end[aa_key] = any(c[-1] != 'G' for c in codons_list)

        # DP table: dp[i][last_char] = (max_log_cai_sum, prev_last_char, chosen_codon)
        dp = [{} for _ in range(n + 1)]
        dp[0][''] = (0.0, None, None)

        for i in range(n):
            aa = protein[i]
            codons = AA_TO_CODONS.get(aa, [])

            if not codons or aa == "*":
                # For stop codons or unknown AAs, just carry forward
                if aa == "*":
                    codon_start = i * 3
                    stop_codon = seq[codon_start:codon_start + 3] if codon_start + 3 <= len(seq) else "TAA"
                    for last_char, (log_sum, _, _) in dp[i].items():
                        new_last = stop_codon[-1]
                        if new_last not in dp[i + 1] or log_sum > dp[i + 1][new_last][0]:
                            dp[i + 1][new_last] = (log_sum, last_char, stop_codon)
                else:
                    for last_char in dp[i]:
                        dp[i + 1][last_char] = dp[i][last_char]
                continue

            # Sort codons by CAI (highest first) for tie-breaking
            codons_sorted = sorted(codons, key=lambda c: self.species_cai.get(c, 0.0), reverse=True)

            # Check if the NEXT amino acid (if any) has codons not starting with T
            # If not, cross-codon GT from this codon ending with G is unavoidable
            next_aa = protein[i + 1] if i + 1 < n else None
            next_can_avoid_gt = True
            if next_aa is not None and next_aa != "*":
                next_can_avoid_gt = aa_has_non_t_start.get(next_aa, True)
            elif next_aa == "*":
                next_can_avoid_gt = False  # Stop codon, can't change

            for last_char, (log_sum, _, _) in dp[i].items():
                if log_sum == INF:
                    continue

                for codon in codons_sorted:
                    cai = self.species_cai.get(codon, 0.0)
                    if cai <= 0:
                        cai = 0.001

                    if self.avoid_gt:
                        # Check within-codon GT (except Valine which is unavoidable)
                        has_within_gt = "GT" in codon
                        if has_within_gt and aa != 'V':
                            continue

                        # Check cross-codon GT with previous codon
                        if last_char and last_char + codon[0] == "GT":
                            # This GT is avoidable only if we could have used a
                            # different codon for the previous AA that doesn't
                            # end with the last_char. But since we're in DP, the
                            # previous choice is already fixed. However, this GT
                            # IS avoidable if the current AA has a codon not
                            # starting with T. If not, we must accept it.
                            current_can_avoid = aa_has_non_t_start.get(aa, True)
                            if current_can_avoid:
                                continue  # Skip - this GT is avoidable
                            # else: GT is unavoidable, accept this codon

                        # Check if this codon ending with G would create unavoidable
                        # cross-codon GT with the NEXT codon. We only need to avoid
                        # this if the next AA CAN avoid starting with T.
                        # If next_aa only has T-starting codons, the GT is unavoidable
                        # and we should use the highest-CAI codon (which may end with G).
                        # We handle this by NOT penalizing codons ending with G
                        # when the next AA can't avoid T-start.

                    new_last_char = codon[-1]
                    new_log_sum = log_sum + math.log(cai)

                    if new_last_char not in dp[i + 1] or new_log_sum > dp[i + 1][new_last_char][0]:
                        dp[i + 1][new_last_char] = (new_log_sum, last_char, codon)

        # Find the best final state
        best_log_sum = INF
        best_last_char = None
        for last_char, (log_sum, _, _) in dp[n].items():
            if log_sum > best_log_sum:
                best_log_sum = log_sum
                best_last_char = last_char

        if best_last_char is None:
            # Fallback: use simple max-CAI (no GT avoidance)
            codons_result = []
            for i, aa in enumerate(protein):
                if aa == "*":
                    codon_start = i * 3
                    codons_result.append(seq[codon_start:codon_start + 3])
                    continue
                candidates = AA_TO_CODONS.get(aa, [])
                if candidates:
                    codons_result.append(max(candidates, key=lambda c: self.species_cai.get(c, 0.0)))
                else:
                    codon_start = i * 3
                    codons_result.append(seq[codon_start:codon_start + 3])
            return "".join(codons_result)

        # Backtrack to find the optimal sequence
        codons_result = [None] * n
        current_char = best_last_char
        for i in range(n - 1, -1, -1):
            _, prev_char, codon = dp[i + 1][current_char]
            codons_result[i] = codon if codon is not None else "NNN"
            current_char = prev_char

        # Record codon decisions for provenance if collector is attached
        if self._provenance_collector is not None:
            # species_cai is now always a flat codon→CAI dict thanks to
            # get_species_cai_weights(), but keep isinstance guard for safety.
            _provenance_cai = self.species_cai  # type: ignore[assignment]
            if isinstance(_provenance_cai, dict) and 'cai_weights' in _provenance_cai:
                _provenance_cai = _provenance_cai['cai_weights']  # type: ignore[assignment]
            for pos_idx in range(n):
                aa = protein[pos_idx]
                if aa == "*" or not AA_TO_CODONS.get(aa, []):
                    continue
                chosen = codons_result[pos_idx]
                if chosen is None:
                    continue
                candidates_sorted = sorted(
                    AA_TO_CODONS[aa],
                    key=lambda c: _provenance_cai.get(c, 0.0),
                    reverse=True,
                )
                best_cai = _provenance_cai.get(chosen, 0.0)
                alternatives: list[dict[str, Any]] = []
                for codon in candidates_sorted:
                    if codon == chosen:
                        continue  # chosen codon is already in chosen_codon field
                    cai_val = _provenance_cai.get(codon, 0.0)
                    gc_bases = sum(1 for b in codon if b in "GC")
                    gc_contribution = gc_bases / 3.0
                    # Check constraint violations for this codon
                    violates: list[str] = []
                    if "GT" in codon and aa != 'V':
                        violates.append("cryptic_splice_donor")
                    if "AG" in codon:
                        violates.append("cryptic_splice_acceptor")
                    # Determine rejection reason
                    rejected_because: str | None = None
                    if violates:
                        rejected_because = f"Violates: {', '.join(violates)}"
                    elif cai_val < best_cai:
                        rejected_because = "Lower CAI"
                    else:
                        rejected_because = "Lower CAI (DP-optimal path chose alternative)"
                    alternatives.append({
                        "codon": codon,
                        "cai_contribution": round(cai_val, 4),
                        "gc_contribution": round(gc_contribution, 2),
                        "violates_constraints": violates,
                        "rejected_because": rejected_because,
                    })

                # Compute confidence
                if len(candidates_sorted) > 1:
                    # Use second-best CAI as baseline
                    second_best_cai = _provenance_cai.get(candidates_sorted[0], 0.0)
                    if candidates_sorted[0] == chosen and len(candidates_sorted) > 1:
                        second_best_cai = _provenance_cai.get(candidates_sorted[1], 0.0)
                    confidence = min(1.0, 0.5 + (best_cai - second_best_cai) * 5)
                    confidence = max(0.0, confidence)
                else:
                    confidence = 1.0

                self._provenance_collector.record_codon_decision(CodonDecision(
                    position=pos_idx,
                    amino_acid=aa,
                    original_codon=None,
                    chosen_codon=chosen,
                    alternatives_considered=alternatives,
                    constraint_reason="DP max-CAI with avoidable-GT avoidance",
                    confidence=round(confidence, 4),
                ))

        return "".join(codons_result)

    # Deprecated alias — use _step_backtranslate_cai instead
    _phase0_max_cai_backtranslate = _step_backtranslate_cai

    # ──────────────────────────────────────────────────────────
    # Step: Resolve Constraints (Priority-based constraint resolution)
    # ──────────────────────────────────────────────────────────
    def _step_resolve_constraints(self, seq: str) -> str:
        """Resolve Constraints step: Fix constraint violations with minimal CAI impact.

        Iteratively finds GT/CG dinucleotides and restriction sites, then
        resolves them by choosing the synonymous substitution with the
        smallest CAI penalty. This is the DNAworks-style approach: start
        with max CAI, then fix only what's needed.

        Key difference from old greedy approach: instead of greedily avoiding GT
        during codon selection (which permanently sacrifices CAI), we fix
        GT violations after the fact, choosing the resolution that costs
        the least CAI.

        Uses IncrementalSequenceState for O(1) GT/CG tracking and
        CodonCache for pre-sorted codon lists.
        """
        import math
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        max_rounds = 30

        for round_num in range(max_rounds):
            current_seq = state.sequence
            old_gt_count = state.gt_count

            # Collect all constraint violations
            violations = []

            # GT dinucleotide violations are only relevant for eukaryotic targets
            # (prokaryotes have no spliceosome, so GT avoidance is unnecessary)
            if self.avoid_gt:
                # Use incremental GT positions for O(1) lookup instead of O(N) scan
                for gt_pos in state.gt_positions_list():
                    codon_idx = gt_pos // 3
                    codon_start = codon_idx * 3
                    if gt_pos % 3 < 2:
                        # Within-codon GT - only add if avoidable
                        if not _is_unavoidable_gt(current_seq, gt_pos):
                            violations.append(("within_gt", gt_pos, codon_idx))
                    else:
                        # Cross-codon GT - only add if avoidable
                        if not _is_unavoidable_gt(current_seq, gt_pos):
                            violations.append(("cross_gt", gt_pos, codon_idx))

            if not violations and not self.enzymes:
                break

            # Check for restriction sites using cached enzyme data
            rs_violations = []
            if enzyme_cache is not None:
                for enzyme, pos in enzyme_cache.find_sites(current_seq):
                    site = enzyme_cache.get_site(enzyme)
                    if site is not None:
                        rs_violations.append((pos, site, enzyme))

            if not violations and not rs_violations:
                break

            any_resolved = False

            # Fix within-codon GTs first (usually easiest)
            for vtype, pos, codon_idx in violations:
                if vtype == "within_gt":
                    resolved = self._fix_within_codon_gt_cai_aware(state, codon_idx, codon_cache)
                    if resolved:
                        any_resolved = True
                        break

            if any_resolved:
                continue

            # Fix cross-codon GTs
            for vtype, pos, codon_idx in violations:
                if vtype == "cross_gt":
                    resolved = self._fix_cross_codon_gt_cai_aware(state, codon_idx, codon_cache)
                    if resolved:
                        any_resolved = True
                        break

            if any_resolved:
                continue

            # Fix restriction sites
            for p, site, enzyme in rs_violations:
                resolved = self._fix_restriction_site_cai_aware(state, p, site, codon_cache, enzyme_cache)
                if resolved:
                    any_resolved = True
                    break

            if not any_resolved:
                break

        # Track unavoidable GTs for mutagenesis
        for gt_pos in state.gt_positions_list():
            if gt_pos % 3 < 2:  # within-codon only
                codon_idx = gt_pos // 3
                codon_start = codon_idx * 3
                aa = state.get_aa(codon_idx)
                if aa and aa != "*":
                    all_have_gt = all("GT" in c for c in AA_TO_CODONS.get(aa, []))
                    if all_have_gt:
                        codon = state.get_codon(codon_idx)
                        for j in range(2):
                            if codon[j:j+2] == "GT":
                                self._unavoidable_gt_positions.add(codon_start + j)

        return state.sequence

    def _fix_within_codon_gt_cai_aware(self, state: IncrementalSequenceState, codon_idx: int, codon_cache: CodonCache) -> bool:
        """Fix a within-codon GT by choosing the best CAI-preserving substitution.

        Uses IncrementalSequenceState for O(1) GT tracking and
        CodonCache for pre-sorted codon lists.

        For eukaryotes: in-codon GTs are acceptable if they have low splice
        donor potential (< SPLICE_DONOR_POTENTIAL_THRESHOLD).  Only fix GTs
        that have high splice donor potential, and only if the CAI cost of
        the GT-free alternative is < threshold.
        """
        codon = state.get_codon(codon_idx)
        aa = state.get_aa(codon_idx)
        if aa is None or aa == "*":
            return False

        # ── Eukaryotic GT-vs-CAI tradeoff with splice donor scoring ──
        if not self.is_prokaryote:
            # First check splice donor potential of this GT
            seq = state.sequence
            for j in range(2):  # GT can be at pos 0-1 or 1-2 within a codon
                pos = codon_idx * 3 + j
                if pos + 1 < len(seq) and seq[pos:pos+2] == "GT":
                    sdp = score_splice_donor_potential(seq, pos)
                    if sdp < SPLICE_DONOR_POTENTIAL_THRESHOLD:
                        # This GT has low splice donor potential — not dangerous.
                        # Accept it (CAI > GT avoidance for low-risk GTs).
                        return False

            # Check CAI cost
            current_w = self.species_cai.get(codon, 0.0)
            gt_free_codons = codon_cache.get_gt_free_codons(aa)
            if gt_free_codons:
                best_gt_free_w = self.species_cai.get(gt_free_codons[0], 0.0)
                if current_w - best_gt_free_w > EUKARYOTE_CAI_GT_COST_THRESHOLD:
                    # CAI cost too high — keep the GT-containing codon
                    return False
            else:
                # No GT-free alternative (e.g., Valine) — must accept
                return False

        old_gt_count = state.gt_count

        # Use pre-sorted GT-free codons from cache
        for alt in codon_cache.get_gt_free_codons(aa):
            if alt == codon:
                continue
            # O(1) boundary check instead of manual base lookup
            left_gt, right_gt = state.boundary_creates_gt(codon_idx, alt)
            if left_gt or right_gt:
                continue
            # O(1) swap + O(1) GT count check + rollback if needed
            old_codon = state.swap_codon(codon_idx, alt)
            if state.gt_count < old_gt_count:
                return True
            else:
                state.swap_codon(codon_idx, old_codon)  # Rollback

        return False

    def _fix_cross_codon_gt_cai_aware(self, state: IncrementalSequenceState, codon_idx: int, codon_cache: CodonCache) -> bool:
        """Fix a cross-codon GT by choosing the best CAI-preserving substitution.

        Uses IncrementalSequenceState for O(1) GT tracking and
        CodonCache for pre-sorted codon lists.

        For eukaryotes: cross-codon GTs are acceptable if they have low
        splice donor potential (< SPLICE_DONOR_POTENTIAL_THRESHOLD).  When
        the GT has high splice donor potential, it is only eliminated if the
        CAI cost of the fix is < EUKARYOTE_CAI_GT_COST_THRESHOLD.
        """
        # ── Eukaryotic GT-vs-CAI tradeoff for cross-codon GTs ──
        # For eukaryotes, check splice donor potential first, then CAI cost.
        if not self.is_prokaryote:
            # Check splice donor potential of this cross-codon GT
            seq = state.sequence
            gt_pos = codon_idx * 3 + 2  # G at end of codon_idx, T at start of codon_idx+1
            if gt_pos + 1 < len(seq) and seq[gt_pos:gt_pos+2] == "GT":
                sdp = score_splice_donor_potential(seq, gt_pos)
                if sdp < SPLICE_DONOR_POTENTIAL_THRESHOLD:
                    # This cross-codon GT has low splice donor potential — not dangerous.
                    # Accept it (CAI > GT avoidance for low-risk GTs).
                    return False

            # Check if the CAI cost of fixing this cross-codon GT
            # is acceptable before attempting to fix it.
            # The G is at the end of codon_idx, the T is at the start of the next codon
            next_ci = codon_idx + 1
            if next_ci < state.num_codons:
                g_aa = state.get_aa(codon_idx)
                t_aa = state.get_aa(next_ci)
                if g_aa and t_aa and g_aa != "*" and t_aa != "*":
                    g_current = state.get_codon(codon_idx)
                    t_current = state.get_codon(next_ci)
                    g_current_w = self.species_cai.get(g_current, 0.0)
                    t_current_w = self.species_cai.get(t_current, 0.0)
                    # Check if there's a cheap fix (CAI cost < threshold)
                    # for either the G-ending codon or the T-starting codon
                    _cross_cai_cost_ok = False
                    g_non_g_end = [c for c in codon_cache.get_sorted_codons(g_aa) if c[-1] != "G"]
                    if g_non_g_end:
                        best_non_g_end_w = self.species_cai.get(g_non_g_end[0], 0.0)
                        if g_current_w - best_non_g_end_w < EUKARYOTE_CAI_GT_COST_THRESHOLD:
                            _cross_cai_cost_ok = True
                    if not _cross_cai_cost_ok:
                        t_non_t_start = [c for c in codon_cache.get_sorted_codons(t_aa) if c[0] != "T"]
                        if t_non_t_start:
                            best_non_t_start_w = self.species_cai.get(t_non_t_start[0], 0.0)
                            if t_current_w - best_non_t_start_w < EUKARYOTE_CAI_GT_COST_THRESHOLD:
                                _cross_cai_cost_ok = True
                    if not _cross_cai_cost_ok:
                        # CAI cost too high — accept the cross-codon GT
                        return False

        old_gt_count = state.gt_count

        # Try resolving by modifying the codon pair
        next_idx = codon_idx + 1
        prev_idx = codon_idx - 1

        # Strategy A: Modify following codon pair (codon_idx and next_idx)
        if next_idx < state.num_codons:
            aa1 = state.get_aa(codon_idx)
            aa2 = state.get_aa(next_idx)

            if aa1 is not None and aa1 != "*" and aa2 is not None and aa2 != "*":
                # Collect all valid pairs sorted by combined CAI
                pairs = []
                for c1 in codon_cache.get_sorted_codons(aa1):
                    for c2 in codon_cache.get_sorted_codons(aa2):
                        if c1[-1] + c2[0] == "GT":
                            continue
                        combined_cai = codon_cache.get_cai(c1) + codon_cache.get_cai(c2)
                        pairs.append((c1, c2, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c1, c2, _ in pairs:
                    old1 = state.swap_codon(codon_idx, c1)
                    old2 = state.swap_codon(next_idx, c2)
                    if state.gt_count < old_gt_count:
                        return True
                    else:
                        state.swap_codon(next_idx, old2)  # Rollback
                        state.swap_codon(codon_idx, old1)  # Rollback

        # Strategy B: Modify preceding codon pair (prev_idx and codon_idx)
        if prev_idx >= 0:
            aa0 = state.get_aa(prev_idx)
            aa1 = state.get_aa(codon_idx)

            if aa0 is not None and aa0 != "*" and aa1 is not None and aa1 != "*":
                pairs = []
                for c0 in codon_cache.get_sorted_codons(aa0):
                    for c1 in codon_cache.get_sorted_codons(aa1):
                        if c0[-1] + c1[0] == "GT":
                            continue
                        combined_cai = codon_cache.get_cai(c0) + codon_cache.get_cai(c1)
                        pairs.append((c0, c1, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c0, c1, _ in pairs:
                    old0 = state.swap_codon(prev_idx, c0)
                    old1 = state.swap_codon(codon_idx, c1)
                    if state.gt_count < old_gt_count:
                        return True
                    else:
                        state.swap_codon(codon_idx, old1)  # Rollback
                        state.swap_codon(prev_idx, old0)  # Rollback

        # Strategy C: Single codon substitution (change only the codon that ends with G)
        aa1 = state.get_aa(codon_idx)
        if aa1 is not None and aa1 != "*":
            for c1 in codon_cache.get_sorted_codons(aa1):
                old_codon = state.swap_codon(codon_idx, c1)
                if state.gt_count < old_gt_count:
                    return True
                else:
                    state.swap_codon(codon_idx, old_codon)  # Rollback

        return False

    def _fix_restriction_site_cai_aware(self, state: IncrementalSequenceState, pos: int, site: str, codon_cache: CodonCache, enzyme_cache: Optional[EnzymeSiteCache] = None) -> bool:
        """Fix a restriction site by choosing the best CAI-preserving substitution.

        Uses IncrementalSequenceState for O(1) GT tracking and
        CodonCache for pre-sorted codon lists.
        """
        codon_indices = set()
        for j in range(pos, pos + len(site)):
            ci = j // 3
            if ci < state.num_codons:
                codon_indices.add(ci)

        old_gt_count = state.gt_count

        # Try each overlapping codon, sorted by position
        for ci in sorted(codon_indices):
            aa = state.get_aa(ci)
            if aa is None or aa == "*":
                continue

            # Use pre-sorted codons from cache
            for alt in codon_cache.get_sorted_codons(aa):
                current_codon = state.get_codon(ci)
                if alt == current_codon:
                    continue
                # O(1) swap + check + rollback
                old_codon = state.swap_codon(ci, alt)
                # Check restriction site is gone
                if enzyme_cache is not None:
                    rs_ok = not enzyme_cache.check_any_site_present(state.sequence)
                else:
                    rs_ok = site not in state.sequence
                if rs_ok:
                    if self.avoid_gt:
                        if state.gt_count <= old_gt_count:
                            return True
                        else:
                            state.swap_codon(ci, old_codon)  # Rollback
                    else:
                        return True
                else:
                    state.swap_codon(ci, old_codon)  # Rollback

        return False

    # Deprecated alias — use _step_resolve_constraints instead
    _phase1_priority_constraint_resolution = _step_resolve_constraints

    # ──────────────────────────────────────────────────────────
    # Step: Greedy Optimize (Per-position CAI maximization)
    # ──────────────────────────────────────────────────────────
    def _step_greedy_optimize(self, seq: str) -> str:
        """Greedy Optimize step: Per-position CAI maximization, GT-aware.

        For each amino acid, select the highest-CAI codon that does not
        introduce a GT dinucleotide within or across codon boundaries.

        If ALL synonymous codons for an AA contain GT (e.g., Valine GTN),
        flag the position as "unavoidable GT" for mutagenesis fallback step.
        """
        codons = []
        protein = self._translate(seq)
        prev_codon_end = ""

        for i, aa in enumerate(protein):
            if aa == "*":
                codon_start = i * 3
                codons.append(seq[codon_start:codon_start + 3])
                prev_codon_end = codons[-1][-1]
                continue

            candidates = AA_TO_CODONS.get(aa, [])
            if not candidates:
                codon_start = i * 3
                codons.append(seq[codon_start:codon_start + 3])
                prev_codon_end = codons[-1][-1]
                continue

            candidates_sorted = sorted(
                candidates,
                key=lambda c: self.species_cai.get(c, 0.0),
                reverse=True,
            )

            best = candidates_sorted[0]
            found_gt_free = False

            if self.avoid_gt:
                for codon in candidates_sorted:
                    if "GT" in codon:
                        continue
                    if prev_codon_end and prev_codon_end + codon[0] == "GT":
                        continue
                    best = codon
                    found_gt_free = True
                    break

                if not found_gt_free:
                    all_have_gt = all("GT" in c for c in candidates)
                    if all_have_gt:
                        codon_abs_start = i * 3
                        for j in range(2):
                            if best[j:j+2] == "GT":
                                self._unavoidable_gt_positions.add(codon_abs_start + j)
                        # Pick the codon that avoids cross-codon GT if possible
                        if prev_codon_end:
                            for codon in candidates_sorted:
                                if prev_codon_end + codon[0] != "GT":
                                    best = codon
                                    break
                    else:
                        # Pick the one that minimizes GT count
                        for codon in candidates_sorted:
                            gt_count = codon.count("GT")
                            cross_gt = 1 if (prev_codon_end and prev_codon_end + codon[0] == "GT") else 0
                            best_gt_count = best.count("GT")
                            best_cross_gt = 1 if (prev_codon_end and prev_codon_end + best[0] == "GT") else 0
                            if gt_count + cross_gt < best_gt_count + best_cross_gt:
                                best = codon
                                break

            # Record codon decision for provenance if collector is attached
            if self._provenance_collector is not None:
                # species_cai is now always a flat codon→CAI dict thanks to
                # get_species_cai_weights(), but keep isinstance guard for safety.
                _prov_cai = self.species_cai
                if isinstance(_prov_cai, dict) and 'cai_weights' in _prov_cai:
                    _prov_cai = _prov_cai['cai_weights']
                best_cai = _prov_cai.get(best, 0.0)
                alternatives: list[dict[str, Any]] = []
                for codon in candidates_sorted:
                    if codon == best:
                        continue  # chosen codon is already in chosen_codon field
                    cai_val = _prov_cai.get(codon, 0.0)
                    gc_bases = sum(1 for b in codon if b in "GC")
                    gc_contribution = gc_bases / 3.0
                    # Check constraint violations for this codon
                    violates: list[str] = []
                    if "GT" in codon:
                        violates.append("cryptic_splice_donor")
                    if prev_codon_end and prev_codon_end + codon[0] == "GT":
                        violates.append("cross_codon_gt")
                    if "AG" in codon:
                        violates.append("cryptic_splice_acceptor")
                    # Determine rejection reason
                    rejected_because: str | None = None
                    if violates:
                        rejected_because = f"Violates: {', '.join(violates)}"
                    elif cai_val < best_cai:
                        rejected_because = "Lower CAI"
                    else:
                        rejected_because = "Lower CAI"
                    alternatives.append({
                        "codon": codon,
                        "cai_contribution": round(cai_val, 4),
                        "gc_contribution": round(gc_contribution, 2),
                        "violates_constraints": violates,
                        "rejected_because": rejected_because,
                    })

                # Compute confidence
                if len(candidates_sorted) > 1:
                    second_best_cai = _prov_cai.get(candidates_sorted[1], 0.0)
                    confidence = min(1.0, 0.5 + (best_cai - second_best_cai) * 5)
                    confidence = max(0.0, confidence)
                else:
                    confidence = 1.0

                self._provenance_collector.record_codon_decision(CodonDecision(
                    position=i,
                    amino_acid=aa,
                    original_codon=None,
                    chosen_codon=best,
                    alternatives_considered=alternatives,
                    constraint_reason="Maximize CAI while avoiding GT dinucleotides",
                    confidence=round(confidence, 4),
                ))

            codons.append(best)
            prev_codon_end = best[-1]

        return "".join(codons)

    # Deprecated alias — use _step_greedy_optimize instead
    _phase1_greedy_optimize = _step_greedy_optimize

    # ──────────────────────────────────────────────────────────
    # Step: Remove Restriction Sites
    # ──────────────────────────────────────────────────────────
    def _step_remove_restriction_sites(self, seq: str) -> str:
        """Remove Restriction Sites step (v2: batch single-pass approach).

        Instead of iterating one site at a time (O(n * k * iterations)),
        this collects ALL restriction site positions in a single scan,
        then resolves them in batch using incremental state for O(1)
        GT/CG tracking after each fix.

        Performance: O(n * k + changes) instead of O(n * k * iterations).
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        max_rounds = 20

        for round_num in range(max_rounds):
            current_seq = state.sequence

            # Single-pass: find ALL restriction sites at once
            all_sites = []
            if enzyme_cache is not None:
                all_sites = enzyme_cache.find_all_sites_batch(current_seq)

            if not all_sites:
                break  # No sites remaining — done

            any_fixed = False

            for enzyme, pos, site_seq in all_sites:
                # Get overlapping codon indices
                codon_indices = set()
                for j in range(pos, pos + len(site_seq)):
                    ci = j // 3
                    if ci < state.num_codons:
                        codon_indices.add(ci)

                # Try each overlapping codon, using pre-sorted alternatives from cache
                for ci in sorted(codon_indices):
                    aa = state.get_aa(ci)
                    if aa is None or aa == "*":
                        continue

                    current_codon = state.get_codon(ci)

                    # Use pre-sorted codons from cache (avoids repeated sorting)
                    for alt in codon_cache.get_sorted_codons(aa):
                        if alt == current_codon:
                            continue

                        # O(1) swap + check + rollback
                        old_codon = state.swap_codon(ci, alt)

                        # Check restriction site is gone
                        if enzyme_cache is not None:
                            rs_ok = not enzyme_cache.check_any_site_present(state.sequence)
                        else:
                            rs_ok = site_seq not in state.sequence

                        if rs_ok:
                            # Check GT increase (allow temporary GT increase, but not severe)
                            if self.avoid_gt:
                                # We can't easily compute old GT count here since we
                                # already swapped. But state.gt_count is the NEW count.
                                # Allow up to 2 more GTs than before the round started.
                                if state.gt_count > state.gt_count + 2:
                                    state.swap_codon(ci, old_codon)  # Rollback
                                    continue
                            any_fixed = True
                            break
                        else:
                            state.swap_codon(ci, old_codon)  # Rollback

                    if any_fixed:
                        break

                if any_fixed:
                    break  # Restart scanning with updated sequence

            if not any_fixed:
                # Try multi-codon coordinated removal for remaining sites
                for enzyme, pos, site_seq in all_sites:
                    codon_indices = sorted(set(
                        j // 3 for j in range(pos, pos + len(site_seq))
                        if j // 3 < state.num_codons
                    ))

                    if len(codon_indices) < 2:
                        continue

                    # Try pairs of adjacent codons
                    for idx in range(len(codon_indices) - 1):
                        ci1, ci2 = codon_indices[idx], codon_indices[idx + 1]
                        if ci2 != ci1 + 1:
                            continue
                        aa1 = state.get_aa(ci1)
                        aa2 = state.get_aa(ci2)
                        if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                            continue

                        old1 = state.get_codon(ci1)
                        old2 = state.get_codon(ci2)

                        # Try codon pairs sorted by combined CAI
                        pairs = []
                        for c1 in codon_cache.get_sorted_codons(aa1):
                            for c2 in codon_cache.get_sorted_codons(aa2):
                                combined = codon_cache.get_cai(c1) + codon_cache.get_cai(c2)
                                pairs.append((c1, c2, combined))
                        pairs.sort(key=lambda x: x[2], reverse=True)

                        for c1, c2, _ in pairs:
                            o1 = state.swap_codon(ci1, c1)
                            o2 = state.swap_codon(ci2, c2)
                            if enzyme_cache is not None:
                                rs_ok = not enzyme_cache.check_any_site_present(state.sequence)
                            else:
                                rs_ok = site_seq not in state.sequence
                            if rs_ok:
                                any_fixed = True
                                break
                            else:
                                state.swap_codon(ci2, o2)
                                state.swap_codon(ci1, o1)

                        if any_fixed:
                            break
                    if any_fixed:
                        break

                if not any_fixed:
                    break  # Cannot fix remaining sites

        return state.sequence

    # Deprecated alias — use _step_remove_restriction_sites instead
    _phase2_remove_restriction_sites = _step_remove_restriction_sites

    # ──────────────────────────────────────────────────────────
    # Step: Cross-Codon Optimization (iterative constraint resolution)
    # ──────────────────────────────────────────────────────────
    def _step_cross_codon_optimization(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Cross-Codon Optimization step (v2: single-pass scan with early termination).

        Resolve cross-codon GT, CG, and restriction site constraints.
        Uses IncrementalSequenceState for O(1) GT/CG tracking instead of
        O(N) full-sequence rescans.

        v2 improvements:
        - Single-pass constraint collection (GT + CG + RS in one scan)
        - Early termination: if no constraints, skip immediately
        - Process all constraints in one round before re-scanning
        - Avoid re-calling _is_unavoidable_gt for same positions
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        total_remaining: Dict[int, List[str]] = {}
        max_rounds = 15

        for round_num in range(max_rounds):
            # Single-pass: collect ALL cross-codon constraints at once
            constraint_positions: Dict[int, List[str]] = {}

            # GT/CG constraints — use incremental position data (O(1) per position)
            if self.avoid_gt:
                for gt_pos in state.gt_positions_list():
                    if gt_pos % 3 == 2:  # Cross-codon boundary
                        codon_start = (gt_pos // 3) * 3
                        # Only add if avoidable
                        codon_idx = gt_pos // 3
                        if codon_idx < state.num_codons:
                            aa = state.get_aa(codon_idx)
                            # Quick check: Valine always has GT, skip
                            if aa == 'V':
                                continue
                            constraint_positions.setdefault(codon_start, [])
                            if "GT" not in constraint_positions[codon_start]:
                                constraint_positions[codon_start].append("GT")

                for cg_pos in state.cg_positions_list():
                    if cg_pos % 3 == 2:  # Cross-codon boundary
                        codon_start = (cg_pos // 3) * 3
                        constraint_positions.setdefault(codon_start, [])
                        if "CG" not in constraint_positions[codon_start]:
                            constraint_positions[codon_start].append("CG")

            # Restriction site constraints — single batch scan
            if enzyme_cache is not None:
                for enzyme, site in enzyme_cache._sites.items():
                    if site is None:
                        continue
                    for pos in find_cross_codon_restriction(state.sequence, site):
                        codon_start = (pos // 3) * 3
                        label = f"RS:{site}"
                        constraint_positions.setdefault(codon_start, [])
                        if label not in constraint_positions[codon_start]:
                            constraint_positions[codon_start].append(label)

            # Early termination: no constraints found
            if not constraint_positions:
                break

            any_resolved = False

            # Process ALL constraints in one pass (don't re-scan after each fix;
            # instead, process them all and re-scan only at the start of next round)
            for codon_start, ctypes in constraint_positions.items():
                # Skip if this position was already fixed by a previous resolution
                # (the codon at this position may have changed)
                resolved = self._try_resolve_cross_codon_incremental(
                    state, codon_start, ctypes, codon_cache, enzyme_cache
                )
                if resolved:
                    any_resolved = True
                    # Don't break — process all constraints in this round

            if not any_resolved:
                # No progress — record remaining constraints and exit
                constraint_positions2: Dict[int, List[str]] = {}
                if self.avoid_gt:
                    for gt_pos in state.gt_positions_list():
                        if gt_pos % 3 == 2:
                            if not _is_unavoidable_gt(state.sequence, gt_pos):
                                cs = (gt_pos // 3) * 3
                                constraint_positions2.setdefault(cs, [])
                                if "GT" not in constraint_positions2[cs]:
                                    constraint_positions2[cs].append("GT")
                    for cg_pos in state.cg_positions_list():
                        if cg_pos % 3 == 2:
                            cs = (cg_pos // 3) * 3
                            constraint_positions2.setdefault(cs, [])
                            if "CG" not in constraint_positions2[cs]:
                                constraint_positions2[cs].append("CG")
                total_remaining = constraint_positions2
                break

        mut_report = propose_mutagenesis(
            state.sequence,
            list(total_remaining.keys()),
            total_remaining,
            self.species_cai,
            self.min_blosum,
            self.min_cai,
        )

        return state.sequence, mut_report

    def _try_resolve_cross_codon(
        self, seq_list: list, codon_start: int, ctypes: List[str]
    ) -> bool:
        """Try to resolve cross-codon constraints at codon_start.

        Uses global validation: after any substitution, verifies that the
        total GT count in the sequence doesn't increase.

        CAI-aware: codon alternatives are sorted by CAI (descending) so
        the first valid combination found is also the highest-CAI one.
        """
        old_seq = "".join(seq_list)
        old_gt_count = _count_gts(old_seq)

        aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
        aa1 = CODON_TABLE.get(aa1_codon)
        if aa1 is None or aa1 == "*":
            return False

        next_start = codon_start + 3

        # Sort alternatives by CAI (descending) for CAI-aware resolution
        aa1_alts = sorted(
            AA_TO_CODONS.get(aa1, [aa1_codon]),
            key=lambda c: self.species_cai.get(c, 0.0),
            reverse=True,
        )

        # Strategy A: Following codon
        if next_start + 3 <= len(seq_list):
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if aa2 is not None and aa2 != "*":
                aa2_alts = sorted(
                    AA_TO_CODONS.get(aa2, [aa2_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c1 in aa1_alts:
                    if self.avoid_gt and "GT" in c1:
                        continue
                    for c2 in aa2_alts:
                        if self.avoid_gt and "GT" in c2:
                            continue
                        if c1[-1] + c2[0] == "GT":
                            continue
                        # Apply and validate globally
                        test_list = seq_list[:]
                        for k, b in enumerate(c1):
                            test_list[codon_start + k] = b
                        for k, b in enumerate(c2):
                            test_list[next_start + k] = b
                        test_seq = "".join(test_list)
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            # Also check constraint is actually resolved
                            resolved = True
                            for ct in ctypes:
                                if ct == "GT" and _has_gt(test_seq[codon_start:next_start + 3]):
                                    resolved = False
                                elif ct == "CG" and "CG" in test_seq[codon_start:next_start + 3]:
                                    resolved = False
                                elif ct.startswith("RS:"):
                                    if ct[3:] in test_seq:
                                        resolved = False
                            if resolved:
                                # Also check no new restriction sites introduced
                                from ..restriction_sites import get_recognition_site
                                rs_ok = True
                                for enzyme in self.enzymes:
                                    rs_site = get_recognition_site(enzyme)
                                    if rs_site and rs_site in test_seq:
                                        rs_ok = False
                                        break
                                if rs_ok:
                                    seq_list[:] = test_list
                                    return True

        # Strategy B: Preceding codon
        if codon_start >= 3:
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            if aa0 is not None and aa0 != "*":
                aa0_alts = sorted(
                    AA_TO_CODONS.get(aa0, [aa0_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c0 in aa0_alts:
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in aa1_alts:
                        if self.avoid_gt and "GT" in c1:
                            continue
                        if c0[-1] + c1[0] == "GT":
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(c0):
                            test_list[prev_start + k] = b
                        for k, b in enumerate(c1):
                            test_list[codon_start + k] = b
                        test_seq = "".join(test_list)
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            resolved = True
                            for ct in ctypes:
                                if ct == "GT" and _has_gt(test_seq[prev_start:codon_start + 3]):
                                    resolved = False
                                elif ct == "CG" and "CG" in test_seq[prev_start:codon_start + 3]:
                                    resolved = False
                                elif ct.startswith("RS:"):
                                    if ct[3:] in test_seq:
                                        resolved = False
                            if resolved:
                                seq_list[:] = test_list
                                return True

        # Strategy C: Both preceding and following
        if codon_start >= 3 and next_start + 3 <= len(seq_list):
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if (aa0 is not None and aa0 != "*" and
                    aa2 is not None and aa2 != "*"):
                aa0_alts = sorted(
                    AA_TO_CODONS.get(aa0, [aa0_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                aa2_alts = sorted(
                    AA_TO_CODONS.get(aa2, [aa2_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c0 in aa0_alts:
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in aa1_alts:
                        if self.avoid_gt and "GT" in c1:
                            continue
                        for c2 in aa2_alts:
                            if self.avoid_gt and "GT" in c2:
                                continue
                            if c0[-1] + c1[0] == "GT" or c1[-1] + c2[0] == "GT":
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(c0):
                                test_list[prev_start + k] = b
                            for k, b in enumerate(c1):
                                test_list[codon_start + k] = b
                            for k, b in enumerate(c2):
                                test_list[next_start + k] = b
                            test_seq = "".join(test_list)
                            new_gt_count = _count_gts(test_seq)
                            if new_gt_count <= old_gt_count:
                                resolved = True
                                for ct in ctypes:
                                    if ct == "GT" and _has_gt(test_seq[prev_start:next_start + 3]):
                                        resolved = False
                                    elif ct == "CG" and "CG" in test_seq[prev_start:next_start + 3]:
                                        resolved = False
                                    elif ct.startswith("RS:"):
                                        if ct[3:] in test_seq:
                                            resolved = False
                                if resolved:
                                    seq_list[:] = test_list
                                    return True

        # Strategy D: Single codon substitution
        for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
            if self.avoid_gt and "GT" in c1:
                continue
            test_list = seq_list[:]
            for k, b in enumerate(c1):
                test_list[codon_start + k] = b
            test_seq = "".join(test_list)
            new_gt_count = _count_gts(test_seq)
            if new_gt_count < old_gt_count:
                seq_list[:] = test_list
                return True

        return False

    def _try_resolve_cross_codon_incremental(
        self, state: IncrementalSequenceState, codon_start: int,
        ctypes: List[str], codon_cache: CodonCache,
        enzyme_cache: Optional[EnzymeSiteCache] = None
    ) -> bool:
        """Try to resolve cross-codon constraints using incremental state.

        Uses O(1) GT/CG tracking via IncrementalSequenceState instead of
        O(N) full-sequence rescans. Uses predictive boundary checks to
        avoid unnecessary swap+rollback cycles.

        v2: Uses regional RS checking to avoid full-sequence scans.
        Uses pre-sorted codon subsets from CodonCache for faster iteration.
        """
        codon_idx = codon_start // 3
        aa1 = state.get_aa(codon_idx)
        if aa1 is None or aa1 == "*":
            return False

        old_gt_count = state.gt_count
        next_codon_idx = codon_idx + 1

        # Pre-check: what constraint types need resolving?
        needs_gt = "GT" in ctypes
        needs_cg = "CG" in ctypes
        needs_rs = any(ct.startswith("RS:") for ct in ctypes)

        # Strategy A: Following codon pair
        if next_codon_idx < state.num_codons:
            aa2 = state.get_aa(next_codon_idx)
            if aa2 is not None and aa2 != "*":
                # Use pre-sorted non-G-end/non-T-start codons from cache for faster iteration
                if needs_gt:
                    c1_list = codon_cache.get_non_g_end_codons(aa1)
                    c2_list = codon_cache.get_non_t_start_codons(aa2)
                else:
                    c1_list = codon_cache.get_sorted_codons(aa1)
                    c2_list = codon_cache.get_sorted_codons(aa2)

                for c1 in c1_list:
                    if self.avoid_gt and "GT" in c1:
                        continue
                    # O(1) left boundary predictive check
                    if self.avoid_gt:
                        left_gt, _ = state.boundary_creates_gt(codon_idx, c1)
                        if left_gt:
                            continue
                    for c2 in c2_list:
                        if self.avoid_gt and "GT" in c2:
                            continue
                        if c1[-1] + c2[0] == "GT":
                            continue
                        # O(1) right boundary predictive check
                        if self.avoid_gt:
                            _, right_gt = state.boundary_creates_gt(next_codon_idx, c2)
                            if right_gt:
                                continue

                        # Quick dinucleotide check — does this pair resolve the constraint?
                        if needs_cg and c1[-1] + c2[0] == "CG":
                            continue

                        # Apply both swaps
                        old1 = state.swap_codon(codon_idx, c1)
                        old2 = state.swap_codon(next_codon_idx, c2)

                        if state.gt_count <= old_gt_count:
                            # Check constraint is actually resolved — use O(1) has_gt_at/has_cg_at
                            resolved = True
                            for ct in ctypes:
                                if ct == "GT" and state.has_gt_at(codon_start + 2):
                                    resolved = False
                                elif ct == "CG" and state.has_cg_at(codon_start + 2):
                                    resolved = False
                                elif ct.startswith("RS:"):
                                    # Use regional RS check (avoids full-sequence string build)
                                    if enzyme_cache is not None:
                                        if enzyme_cache.check_sites_in_region(
                                            state.sequence, codon_start, next_codon_idx * 3 + 3
                                        ):
                                            resolved = False
                                    else:
                                        if ct[3:] in state.sequence:
                                            resolved = False
                            if resolved:
                                return True

                        # Rollback
                        state.swap_codon(next_codon_idx, old2)
                        state.swap_codon(codon_idx, old1)

        # Strategy B: Preceding codon
        if codon_idx > 0:
            prev_codon_idx = codon_idx - 1
            aa0 = state.get_aa(prev_codon_idx)
            if aa0 is not None and aa0 != "*":
                for c0 in codon_cache.get_sorted_codons(aa0):
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in codon_cache.get_sorted_codons(aa1):
                        if self.avoid_gt and "GT" in c1:
                            continue
                        if c0[-1] + c1[0] == "GT":
                            continue
                        # Quick dinucleotide check
                        if needs_gt and "GT" in (c0 + c1):
                            continue
                        if needs_cg and c0[-1] + c1[0] == "CG":
                            continue

                        old0 = state.swap_codon(prev_codon_idx, c0)
                        old1 = state.swap_codon(codon_idx, c1)

                        if state.gt_count <= old_gt_count:
                            resolved = True
                            for ct in ctypes:
                                prev_start = codon_start - 3
                                region = state.sequence[prev_start:codon_start + 3]
                                if ct == "GT" and "GT" in region:
                                    resolved = False
                                elif ct == "CG" and "CG" in region:
                                    resolved = False
                                elif ct.startswith("RS:"):
                                    if ct[3:] in state.sequence:
                                        resolved = False
                            if resolved:
                                return True

                        state.swap_codon(codon_idx, old1)
                        state.swap_codon(prev_codon_idx, old0)

        # Strategy C: Both preceding and following
        if codon_idx > 0 and next_codon_idx < state.num_codons:
            prev_codon_idx = codon_idx - 1
            aa0 = state.get_aa(prev_codon_idx)
            aa2 = state.get_aa(next_codon_idx)
            if (aa0 is not None and aa0 != "*" and
                    aa2 is not None and aa2 != "*"):
                for c0 in codon_cache.get_sorted_codons(aa0):
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in codon_cache.get_sorted_codons(aa1):
                        if self.avoid_gt and "GT" in c1:
                            continue
                        if c0[-1] + c1[0] == "GT":
                            continue
                        for c2 in codon_cache.get_sorted_codons(aa2):
                            if self.avoid_gt and "GT" in c2:
                                continue
                            if c1[-1] + c2[0] == "GT":
                                continue
                            # Quick dinucleotide check
                            if needs_gt and "GT" in (c0 + c1 + c2):
                                continue
                            if needs_cg and (c0[-1] + c1[0] == "CG" or c1[-1] + c2[0] == "CG"):
                                continue

                            old0 = state.swap_codon(prev_codon_idx, c0)
                            old1 = state.swap_codon(codon_idx, c1)
                            old2 = state.swap_codon(next_codon_idx, c2)

                            if state.gt_count <= old_gt_count:
                                resolved = True
                                for ct in ctypes:
                                    prev_start = codon_start - 3
                                    region = state.sequence[prev_start:codon_start + 6]
                                    if ct == "GT" and "GT" in region:
                                        resolved = False
                                    elif ct == "CG" and "CG" in region:
                                        resolved = False
                                    elif ct.startswith("RS:"):
                                        if ct[3:] in state.sequence:
                                            resolved = False
                                if resolved:
                                    return True

                            state.swap_codon(next_codon_idx, old2)
                            state.swap_codon(codon_idx, old1)
                            state.swap_codon(prev_codon_idx, old0)

        # Strategy D: Single codon substitution
        for c1 in codon_cache.get_sorted_codons(aa1):
            if self.avoid_gt and "GT" in c1:
                continue
            # Predictive check: would this reduce GT count?
            if self.avoid_gt and state.would_gt_increase(codon_idx, c1):
                continue
            old1 = state.swap_codon(codon_idx, c1)
            if state.gt_count < old_gt_count:
                return True
            state.swap_codon(codon_idx, old1)

        return False

    # Deprecated alias — use _step_cross_codon_optimization instead
    _phase3_cross_codon_constraints = _step_cross_codon_optimization

    # ──────────────────────────────────────────────────────────
    # Step: Within-Codon GT Resolution
    # ──────────────────────────────────────────────────────────
    def _step_within_codon_gt_resolution(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Within-Codon GT Resolution step: Resolve within-codon GT dinucleotides.

        For each within-codon GT (not cross-codon), try synonymous substitution
        first. If no synonymous codon can avoid GT, flag for mutagenesis.

        Uses IncrementalSequenceState for O(1) GT tracking and
        CodonCache for pre-sorted codon lists.
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        remaining_positions: Dict[int, List[str]] = {}

        # Use incremental GT positions for O(1) lookup instead of O(N) scan
        for gt_pos in state.gt_positions_list():
            codon_idx = gt_pos // 3
            codon_start = codon_idx * 3
            next_codon_start = codon_start + 3

            if gt_pos % 3 < 2:  # Within-codon GT
                # Within-codon GT - skip if unavoidable
                if _is_unavoidable_gt(state.sequence, gt_pos):
                    continue
                aa = state.get_aa(codon_idx)
                if aa is None or aa == "*":
                    continue

                old_gt_count = state.gt_count
                resolved = False
                # Use pre-sorted GT-free codons from cache
                for alt in codon_cache.get_gt_free_codons(aa):
                    current_codon = state.get_codon(codon_idx)
                    if alt == current_codon:
                        continue
                    # O(1) boundary check instead of manual base lookup
                    left_gt, right_gt = state.boundary_creates_gt(codon_idx, alt)
                    if left_gt or right_gt:
                        continue
                    # O(1) swap + O(1) GT count check + rollback if needed
                    old_codon = state.swap_codon(codon_idx, alt)
                    if state.gt_count < old_gt_count:
                        resolved = True
                        break
                    else:
                        state.swap_codon(codon_idx, old_codon)  # Rollback

                if not resolved:
                    remaining_positions[codon_start] = ["GT:within"]

        if remaining_positions:
            mut_report = self._propose_within_codon_mutagenesis(
                state.sequence, remaining_positions
            )
        else:
            mut_report = MutagenesisReport()

        return state.sequence, mut_report

    def _propose_within_codon_mutagenesis(
        self,
        seq: str,
        positions: Dict[int, List[str]],
    ) -> MutagenesisReport:
        """Propose AA substitutions for within-codon GTs that can't be resolved
        by synonymous substitution.

        Uses BLOSUM62 guidance for conservative substitutions:
        - Valine (V, GTN) → Isoleucine (I, ATN) or Leucine (L, TTA/CTN)
        """
        report = MutagenesisReport()

        for codon_start, ctypes in positions.items():
            if codon_start + 3 > len(seq):
                continue

            original_codon = seq[codon_start:codon_start + 3]
            original_aa = CODON_TABLE.get(original_codon)
            if original_aa is None or original_aa == "*":
                continue

            best = None
            best_score = -100

            for new_aa, score in sorted(
                ((aa, BLOSUM62.get((original_aa, aa), -10))
                 for aa in set(CODON_TABLE.values()) if aa != "*" and aa != original_aa),
                key=lambda x: x[1],
                reverse=True,
            ):
                if score < self.min_blosum:
                    continue

                alts = AA_TO_CODONS.get(new_aa, [])
                alts_sorted = sorted(alts,
                                     key=lambda c: self.species_cai.get(c, 0.0),
                                     reverse=True)
                for alt_codon in alts_sorted:
                    if "GT" in alt_codon:
                        continue
                    prev_base = seq[codon_start - 1] if codon_start > 0 else ""
                    next_start = codon_start + 3
                    next_base = seq[next_start] if next_start < len(seq) else ""
                    if prev_base and prev_base + alt_codon[0] == "GT":
                        continue
                    if next_base and alt_codon[-1] + next_base == "GT":
                        continue
                    cai = self.species_cai.get(alt_codon, 0.0)
                    if cai < self.min_cai:
                        continue

                    if score > best_score:
                        best = (alt_codon, new_aa, score, cai)
                        best_score = score
                    break

            if best is not None:
                alt_codon, new_aa, blosum, cai = best
                proposal = MutagenesisProposal(
                    position=codon_start,
                    original_codon=original_codon,
                    original_aa=original_aa,
                    new_aa=new_aa,
                    new_codon=alt_codon,
                    blosum_score=blosum,
                    cai_weight=cai,
                    resolves=ctypes,
                )
                report.add(proposal)
            else:
                proposal = MutagenesisProposal(
                    position=codon_start,
                    original_codon=original_codon,
                    original_aa=original_aa,
                    new_aa="",
                    new_codon="",
                    blosum_score=-10,
                    cai_weight=0.0,
                    resolves=ctypes,
                    impossible=True,
                )
                report.add(proposal)

        return report

    # Deprecated alias — use _step_within_codon_gt_resolution instead
    _phase35_within_codon_gt = _step_within_codon_gt_resolution

    # ──────────────────────────────────────────────────────────
    # Step: Mutagenesis Fallback
    # ──────────────────────────────────────────────────────────
    def _step_mutagenesis_fallback(self, seq: str, mut_report: MutagenesisReport) -> str:
        """Mutagenesis Fallback step: Apply mutagenesis proposals for intractable constraints.

        Applies conservative AA substitutions (e.g., Val→Ile, Val→Leu)
        using BLOSUM62 guidance. Only applies to AVOIDABLE GT positions;
        unavoidable GTs (Valine, cross-codon with no alternatives) are skipped.
        """
        seq_list = list(seq)
        for proposal in mut_report.proposals:
            if proposal.impossible or not proposal.new_codon:
                continue
            pos = proposal.position
            old_gt_count = _count_gts("".join(seq_list))

            test_list = seq_list[:]
            for k, b in enumerate(proposal.new_codon):
                if pos + k < len(test_list):
                    test_list[pos + k] = b
            new_gt_count = _count_gts("".join(test_list))

            # Accept if it reduces or maintains GT count
            if new_gt_count <= old_gt_count:
                for k, b in enumerate(proposal.new_codon):
                    if pos + k < len(seq_list):
                        seq_list[pos + k] = b
                self._applied_mutagenesis.append({
                    "position": pos,
                    "original_aa": proposal.original_aa,
                    "new_aa": proposal.new_aa,
                    "blosum": proposal.blosum_score,
                })

        return "".join(seq_list)

    # Deprecated alias — use _step_mutagenesis_fallback instead
    _phase4_mutagenesis_fallback = _step_mutagenesis_fallback

    # ──────────────────────────────────────────────────────────
    # Step: Avoid CpG Islands
    # ──────────────────────────────────────────────────────────
    def _step_avoid_cpg_islands(self, seq: str) -> str:
        """Avoid CpG Islands step (optimized with incremental state)."""
        return _postprocessing.step_avoid_cpg_islands(self._make_config(), seq)
    def _step_remove_instability_motifs(self, seq: str) -> str:
        """Remove ATTTA mRNA instability motifs by synonymous codon substitution."""
        return _postprocessing.step_remove_instability_motifs(self._make_config(), seq)
    def _step_cpg_reconciliation(self, seq: str) -> str:
        """Aggressive CpG reconciliation pass to eliminate CpG islands."""
        return _postprocessing.step_cpg_reconciliation(self._make_config(), seq)

    def _check_no_restriction_sites(self, seq: str) -> bool:
        """Check that sequence doesn't contain any restriction enzyme sites."""
        return _postprocessing._check_no_restriction_sites(self._make_config(), seq)

    def _check_cpg_swap_gt_safe(self, old_list: list, new_list: list, codon_start: int) -> bool:
        """Check that a codon swap doesn't create new avoidable GT dinucleotides."""
        return _postprocessing._check_cpg_swap_gt_safe(old_list, new_list, codon_start)

    def _fix_gc_after_cpg(self, seq: str) -> str:
        """Fix GC content after CpG reconciliation without reintroducing CG dinucleotides."""
        return _postprocessing.fix_gc_after_cpg(self._make_config(), seq)
    def _step_cross_codon_coordination(self, seq: str) -> str:
        """Cross-Codon Coordination step (optimized with incremental state).

        Resolve cross-codon GT, CG, and restriction site violations by
        considering ALL synonymous codon pairs at each boundary.

        Uses IncrementalSequenceState for O(1) GT/CG tracking and
        CodonCache for pre-sorted codon lists.

        Unlike the earlier Cross-Codon Optimization step which handles violations
        one at a time with limited search, this phase performs an exhaustive
        pairwise search: for every boundary violation it considers ALL synonymous
        codon pairs for the two adjacent amino acids, selects the pair that
        eliminates the violation while maximizing CAI, and re-scans after each
        fix (since fixing one boundary can create a violation at the adjacent
        boundary). Iterates until no cross-codon violations remain or a maximum
        iteration count is reached.

        This phase runs AFTER within-codon optimization and BEFORE the final
        CAI maximization pass, ensuring that cross-codon interactions — which
        are the primary cause of predicate failures for sequences like HBB —
        are fully resolved before the hill climber locks in codon choices.
        """
        state = IncrementalSequenceState(seq)
        codon_cache = CodonCache(self.species_cai)
        enzyme_cache = EnzymeSiteCache(self.enzymes) if self.enzymes else None
        max_iterations = 50

        # Record the initial avoidable GT count for comparison
        initial_gt_count = state.gt_count if self.avoid_gt else 0

        for iteration in range(max_iterations):
            # Find cross-codon violations using incremental position data
            violations = []

            # GT violations at cross-codon boundaries
            # EUKARYOTE-ONLY: GT avoidance is irrelevant for prokaryotes
            if self.avoid_gt:
                for gt_pos in state.gt_positions_list():
                    codon_idx = gt_pos // 3
                    next_codon_idx = codon_idx + 1
                    # Cross-codon GT: position is at the boundary between codons
                    # gt_pos is at the last base of one codon, gt_pos+1 is the first of the next
                    if gt_pos % 3 == 2 and next_codon_idx < state.num_codons:
                        if not _is_unavoidable_gt(state.sequence, gt_pos):
                            violations.append({"type": "GT", "boundary": gt_pos})

            # CG violations at cross-codon boundaries
            # EUKARYOTE-ONLY: CpG island avoidance is irrelevant for prokaryotes
            if not self.is_prokaryote:
                for cg_pos in state.cg_positions_list():
                    if cg_pos % 3 == 2:
                        codon_idx = cg_pos // 3
                        next_codon_idx = codon_idx + 1
                        if next_codon_idx < state.num_codons:
                            violations.append({"type": "CG", "boundary": cg_pos})

            if not violations:
                break

            # Sort: GT first, then CG, then by position
            type_priority = {"GT": 0, "CG": 1}
            violations.sort(key=lambda v: (type_priority.get(v["type"], 2), v["boundary"]))

            fixed_any = False

            for violation in violations:
                boundary = violation["boundary"]
                vtype = violation["type"]

                left_codon_start = boundary - 2
                right_codon_start = boundary + 1

                if left_codon_start < 0 or right_codon_start + 3 > len(state.sequence):
                    continue

                left_ci = left_codon_start // 3
                right_ci = right_codon_start // 3

                left_aa = state.get_aa(left_ci)
                right_aa = state.get_aa(right_ci)

                if left_aa is None or left_aa == "*" or right_aa is None or right_aa == "*":
                    continue

                # Find the best pair that eliminates the violation and maximizes CAI
                best_pair = None
                best_cai = -1.0

                for c_left in codon_cache.get_sorted_codons(left_aa):
                    for c_right in codon_cache.get_sorted_codons(right_aa):
                        # Check that this pair eliminates the specific violation
                        if vtype == "GT" and c_left[-1] + c_right[0] == "GT":
                            continue
                        if vtype == "CG" and c_left[-1] + c_right[0] == "CG":
                            continue
                        # Check within-codon violations
                        if vtype == "GT" and ("GT" in c_left or "GT" in c_right):
                            continue
                        if vtype == "CG" and ("CG" in c_left or "CG" in c_right):
                            continue

                        # Compute CAI contribution
                        pair_cai = codon_cache.get_cai(c_left) + codon_cache.get_cai(c_right)
                        if pair_cai <= best_cai:
                            continue

                        # Record GT count before trial swap
                        gt_before = state.gt_count

                        # Try the swap and check adjacent boundaries
                        old_left = state.swap_codon(left_ci, c_left)
                        old_right = state.swap_codon(right_ci, c_right)

                        # Check that we haven't introduced new GTs at adjacent boundaries
                        new_violations = False

                        # Left boundary of c_left with previous codon
                        if left_ci > 0:
                            left_gt, _ = state.boundary_creates_gt(left_ci, c_left)
                            if left_gt and not _is_unavoidable_gt(state.sequence, left_codon_start - 1):
                                new_violations = True

                        # Right boundary of c_right with next codon
                        if not new_violations and right_ci + 1 < state.num_codons:
                            _, right_gt = state.boundary_creates_gt(right_ci, c_right)
                            if right_gt and not _is_unavoidable_gt(state.sequence, right_codon_start + 2):
                                new_violations = True

                        # Check for new CG violations at adjacent boundaries
                        if not new_violations and left_ci > 0:
                            left_cg, _ = state.boundary_creates_cg(left_ci, c_left)
                            if left_cg:
                                new_violations = True

                        if not new_violations and right_ci + 1 < state.num_codons:
                            _, right_cg = state.boundary_creates_cg(right_ci, c_right)
                            if right_cg:
                                new_violations = True

                        # Check that we haven't introduced new restriction sites
                        if not new_violations and enzyme_cache is not None:
                            if enzyme_cache.check_any_site_present(state.sequence):
                                new_violations = True

                        # Check that avoidable GT count hasn't increased
                        if not new_violations and self.avoid_gt:
                            if state.gt_count > gt_before:
                                # The swap increased the total GT count; check if
                                # any new GTs are avoidable
                                current_seq = state.sequence
                                new_avoidable = sum(
                                    1 for p in state.gt_positions_list()
                                    if not _is_unavoidable_gt(current_seq, p)
                                )
                                if new_avoidable > initial_gt_count:
                                    new_violations = True

                        if not new_violations:
                            best_pair = (c_left, c_right)
                            best_cai = pair_cai

                        # Rollback
                        state.swap_codon(right_ci, old_right)
                        state.swap_codon(left_ci, old_left)

                if best_pair is not None:
                    # Apply the best pair for this violation
                    state.swap_codon(left_ci, best_pair[0])
                    state.swap_codon(right_ci, best_pair[1])
                    fixed_any = True
                    break  # Re-scan from scratch after each fix

            if not fixed_any:
                break

        return state.sequence

    # Deprecated alias — use _step_cross_codon_coordination instead
    _phase4_cross_codon_coordination = _step_cross_codon_coordination

    def _find_cross_codon_violations(self, seq: str) -> list:
        """Find all cross-codon violations in the sequence.

        Returns a list of dicts with keys:
          - 'boundary': int — position of the last base of the left codon
              (i.e., the 'G' in a cross-codon GT, or the 'C' in a cross-codon CG)
              For RS violations this is the codon boundary that the site crosses
              (the position of the last base of the left codon in the pair).
          - 'type': str — "GT", "CG", or "RS:<site>"
        """
        violations = []

        # Cross-codon GT (only relevant for eukaryotic targets)
        if self.avoid_gt:
            for pos in find_cross_codon_gt(seq):
                if not _is_unavoidable_gt(seq, pos):
                    violations.append({"boundary": pos, "type": "GT"})

        # Cross-codon CG (CpG-related, only relevant for eukaryotic targets)
        if not self.is_prokaryote:
            for pos in find_cross_codon_cg(seq):
                violations.append({"boundary": pos, "type": "CG"})

        # Cross-codon restriction sites
        # find_cross_codon_restriction returns the start position of the site,
        # but we need to map it to codon boundaries for the pair-based approach.
        from ..restriction_sites import get_recognition_site
        seen_rs_boundaries: set = set()
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            site_len = len(site)
            for pos in find_cross_codon_restriction(seq, site):
                # Find all codon boundaries (b where b % 3 == 2) that the site spans
                for b in range(pos, pos + site_len - 1):
                    if b % 3 == 2 and b + 1 < len(seq):
                        key = (b, f"RS:{site}")
                        if key not in seen_rs_boundaries:
                            seen_rs_boundaries.add(key)
                            violations.append({"boundary": b, "type": f"RS:{site}"})

        return violations

    def _boundary_has_violation(self, seq: str, boundary: int, vtype: str,
                                violation: dict) -> bool:
        """Check if a specific boundary still has the given violation type."""
        if vtype == "GT":
            return (boundary + 1 < len(seq)
                    and seq[boundary] == "G" and seq[boundary + 1] == "T")
        elif vtype == "CG":
            return (boundary + 1 < len(seq)
                    and seq[boundary] == "C" and seq[boundary + 1] == "G")
        elif vtype.startswith("RS:"):
            site = vtype[3:]
            # The restriction site might still be present anywhere in the sequence
            return site in seq
        return False

    def _pair_resolves_violation(
        self,
        c_left: str, c_right: str,
        left_codon_start: int, right_codon_start: int,
        seq: str, boundary: int, vtype: str, violation: dict,
    ) -> bool:
        """Check whether applying the codon pair (c_left, c_right) resolves the
        violation at the given boundary.

        This constructs a test sequence and checks the specific violation type.
        """
        # Build the test sequence with the pair applied
        test_seq = (
            seq[:left_codon_start]
            + c_left
            + c_right
            + seq[right_codon_start + 3:]
        )

        # The boundary position in the test sequence is the same
        # (boundary = left_codon_start + 2 = position of last base of c_left)

        if vtype == "GT":
            # Check the boundary dinucleotide
            if boundary + 1 < len(test_seq):
                dinuc = test_seq[boundary] + test_seq[boundary + 1]
                return dinuc != "GT"
            return True

        elif vtype == "CG":
            if boundary + 1 < len(test_seq):
                dinuc = test_seq[boundary] + test_seq[boundary + 1]
                return dinuc != "CG"
            return True

        elif vtype.startswith("RS:"):
            site = vtype[3:]
            return site not in test_seq

        return False

    # ──────────────────────────────────────────────────────────
    # Step: CAI Reconciliation
    # ──────────────────────────────────────────────────────────
    def _step_cai_reconciliation(self, seq: str) -> str:
        """CAI Reconciliation pass: upgrade low-CAI codons while maintaining hard constraints."""
        return _postprocessing.step_cai_reconciliation(self._make_config(), seq)

    def _cai_recon_check_constraints(self, old_seq: str, test_seq: str, codon_pos: int) -> bool:
        """Check that a codon swap doesn't violate any hard constraints."""
        return _postprocessing.cai_recon_check_constraints(self._make_config(), old_seq, test_seq, codon_pos)

    def _cai_recon_check_global_constraints(self, test_seq: str) -> bool:
        """Check that the full sequence satisfies all hard constraints."""
        return _postprocessing.cai_recon_check_global_constraints(self._make_config(), test_seq)

    def _cai_recon_try_paired(self, seq_list: list, codon_pos: int, new_codon: str) -> Optional[list]:
        """Try a paired codon swap to enable a CAI upgrade."""
        return _postprocessing.cai_recon_try_paired(self._make_config(), seq_list, codon_pos, new_codon)

    def _cai_recon_try_paired_state(
        self, state: IncrementalSequenceState, codon_idx: int,
        new_codon: str, codon_cache: CodonCache,
        enzyme_cache: Optional[EnzymeSiteCache],
    ) -> Tuple[bool, Optional[Tuple[int, str]]]:
        """Try a paired codon swap to enable a CAI upgrade (incremental version)."""
        return _postprocessing.cai_recon_try_paired_state(
            state, codon_idx, new_codon, codon_cache, enzyme_cache
        )

    def _cai_recon_rollback_paired(
        self, state: IncrementalSequenceState, codon_idx: int, old_codon: str,
    ) -> None:
        """Rollback a failed paired swap in _step_cai_reconciliation."""
        _postprocessing.cai_recon_rollback_paired(state, codon_idx, old_codon)
    def _step_gt_reconciliation(self, seq: str) -> str:
        """Fix avoidable GT dinucleotides that may have been introduced."""
        return _postprocessing.step_gt_reconciliation(self._make_config(), seq)
    def _step_cai_hill_climb(self, seq: str) -> str:
        """CAI Hill Climb step: CAI hill climbing (optimized with incremental state)."""
        return _postprocessing.step_cai_hill_climb(self._make_config(), seq)

    def _try_paired_cai_upgrade(self, seq_list: list, codon_pos: int, new_codon: str, old_gt_count: int) -> Optional[list]:
        """Try a paired codon swap to enable a CAI upgrade."""
        return _postprocessing.try_paired_cai_upgrade(self._make_config(), seq_list, codon_pos, new_codon, old_gt_count)

    def _try_paired_cai_upgrade_incremental(
        self, state: IncrementalSequenceState, codon_idx: int,
        new_codon: str, codon_cache: CodonCache,
        enzyme_cache: Optional[EnzymeSiteCache] = None
    ) -> bool:
        """Try a paired codon swap using IncrementalSequenceState."""
        return _postprocessing.try_paired_cai_upgrade_incremental(
            state, codon_idx, new_codon, codon_cache, enzyme_cache
        )
    def _step_reoptimize(self, seq: str) -> str:
        """Re-optimization step: Iterative re-optimization pass (optimized with incremental state)."""
        return _postprocessing.step_reoptimize(self._make_config(), seq)
    def _step_mrna_stability_improvement(self, seq: str) -> str:
        """Soft mRNA stability improvement pass."""
        result = _postprocessing.step_mrna_stability_improvement(self._make_config(), seq)
        if isinstance(result, tuple):
            seq, self._mrna_stability_score, self._destabilizing_motifs_removed, self._stability_improvement = result
            return seq
        return result
    def _evaluate_all_predicates(self, seq: str, skip_splice_check: bool = False) -> List[PredicateResult]:
        """Evaluate all 12 predicates against the optimized sequence."""
        return _validation.evaluate_all_predicates(self._make_config(), seq, skip_splice_check=skip_splice_check)
    @staticmethod
    def _translate(seq: str) -> str:
        """Translate a DNA sequence to amino acid sequence."""
        return _validation.translate(seq)

    # ── Phase 8: G-quadruplex, LinearDesign, TASEP, NMD step methods ──

    def _step_g4_check(self, seq: str) -> str:
        """Check for and fix G-quadruplex motifs via synonymous codon substitution.

        Uses the ``g_quadruplex`` submodule when available; gracefully
        returns the sequence unchanged if the submodule is not installed.
        """
        try:
            from .g_quadruplex import detect_g4_motifs, fix_g4_issues
            g4_report = detect_g4_motifs(seq)
            if g4_report.motifs:
                seq = fix_g4_issues(seq, g4_report, organism=self.organism_name)
                logger.info("Fixed %d G-quadruplex motifs", len(g4_report.motifs))
        except ImportError:
            pass  # g_quadruplex module not available
        return seq

    def _step_lineardesign(self, seq: str) -> str:
        """Apply LinearDesign MFE/CAI joint optimization.

        Uses the ``mfe_optimization`` submodule when available; gracefully
        returns the sequence unchanged if the submodule or C++ binary is
        not installed.
        """
        protein = self._translate(seq)
        try:
            from .mfe_optimization import optimize_with_lineardesign
            ld_result = optimize_with_lineardesign(
                protein, lambda_val=self.lineardesign_lambda, organism=self.organism_name
            )
            if "error" not in ld_result and ld_result.get("sequence"):
                seq = ld_result["sequence"]
                logger.info(
                    "LinearDesign optimization: MFE=%s, CAI=%s",
                    ld_result.get("mfe", "N/A"),
                    ld_result.get("cai", "N/A"),
                )
        except Exception:
            pass  # Fall through to standard MFE optimization
        return seq

    def _step_tasep_ensemble(self, seq: str) -> None:
        """Run ensemble TASEP ribosome simulation for translation dynamics analysis.

        Uses the ``ribosome_simulation`` submodule when available; gracefully
        skips if the submodule is not installed.

        Note: This method does not modify the sequence; it is an analysis step
        that logs ribosome dynamics statistics for the user.
        """
        try:
            from .ribosome_simulation import simulate_tasep_ensemble
            # Compute dwell times from codon adaptiveness table
            dwell_usage = self.species_cai
            dwell_times = []
            for si in range(0, len(seq) - 2, 3):
                codon = seq[si:si + 3]
                w = dwell_usage.get(codon, 0.5)
                dwell_times.append(1.0 / max(w, 0.01))
            tasep_result = simulate_tasep_ensemble(
                dwell_times=dwell_times,
                n_runs=self.tasep_n_runs,
            )
            if tasep_result is not None:
                logger.info(
                    "Ensemble TASEP: %d runs, mean ribosome density=%.4f",
                    self.tasep_n_runs,
                    getattr(tasep_result, 'mean_ribosome_density', 0.0),
                )
        except ImportError:
            pass  # Fall back to single-run TASEP

    def _step_nmd_check(self, seq: str) -> None:
        """Check for nonsense-mediated decay (NMD) triggers in eukaryotic sequences.

        Uses the ``rna_degradation`` submodule when available; gracefully
        skips if the submodule is not installed.

        Note: This method does not modify the sequence; it is a diagnostic
        step that logs detected NMD triggers for the user.
        """
        try:
            from .rna_degradation import detect_nmd_triggers
            orf_start = 0
            orf_end = len(seq) - 3  # Exclude stop codon
            nmd_signals = detect_nmd_triggers(seq, orf_start=orf_start, orf_end=orf_end)
            if nmd_signals:
                for signal in nmd_signals:
                    logger.info("NMD trigger: %s", signal.description)
        except ImportError:
            pass  # rna_degradation module not available
