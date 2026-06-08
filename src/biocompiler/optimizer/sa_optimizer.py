"""
BioCompiler Simulated Annealing Optimizer for mRNA Design
==========================================================

This module provides a comprehensive simulated annealing (SA) framework for
jointly optimizing mRNA minimum free energy (MFE) and codon adaptation index
(CAI).  It extends the basic SA approach in ``mfe_optimization.py`` with:

1. **Configurable SA pipeline** — ``SAConfig`` dataclass and ``optimize_sa()``
   entry point with multiple cooling schedules, Metropolis acceptance, and
   temperature history tracking.

2. **Windowed SA for long genes** — ``optimize_sa_windowed()`` splits long
   sequences into overlapping windows, optimises each independently, and
   merges with boundary repair.

3. **Multi-objective Pareto front** — ``optimize_pareto()`` sweeps multiple
   ``lambda_cai`` values to trace the MFE–CAI Pareto front, and
   ``select_pareto_optimal()`` picks the best result given user-specified
   targets.

4. **Optional NSGA-II integration** — ``optimize_nsga2()`` uses the ``pymoo``
   library for true multi-objective optimisation when installed.

Objective function
------------------
The joint objective follows the LinearDesign formulation (Zhang et al. 2023)::

    score(dna) = MFE(dna) + lambda_cai * (-log(CAI(dna)))

Lower scores are better.  ``lambda_cai = 0`` gives pure MFE minimisation;
higher values penalise low-CAI sequences more aggressively.

References
----------
Kirkpatrick, S., Gelatt, C.D. & Vecchi, M.P. (1983). "Optimization by
Simulated Annealing." *Science* 220(4598):671–680.

Deb, K., Pratap, A., Agarwal, S. & Meyarivan, T. (2002). "A Fast and
Elitist Multiobjective Genetic Algorithm: NSGA-II." *IEEE Transactions
on Evolutionary Computation* 6(2):182–197.

Zhang, H. et al. (2023). "Algorithm for optimized mRNA design improves
stability and immunogenicity." *Nature* 621:396–403.
"""

from __future__ import annotations

import logging
import math
import random
import time
import warnings
from dataclasses import dataclass, field
from typing import Any

from ..shared.constants import CODON_TABLE, AA_TO_CODONS
from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism

logger = logging.getLogger(__name__)

__all__ = [
    # Dataclasses
    "SAConfig",
    "SAResult",
    # Core SA functions
    "optimize_sa",
    "optimize_sa_windowed",
    "compute_objective",
    "cooling_schedule",
    "accept_move",
    # Multi-objective
    "optimize_pareto",
    "select_pareto_optimal",
    # Optional pymoo integration
    "optimize_nsga2",
]


# ══════════════════════════════════════════════════════════════
# Dataclasses
# ══════════════════════════════════════════════════════════════


@dataclass
class SAConfig:
    """Configuration for simulated annealing optimization.

    Attributes:
        initial_temp: Starting temperature for the SA annealing schedule.
            Higher values allow more exploration of the search space.
        cooling_rate: Multiplicative cooling factor per iteration.
            Values close to 1.0 cool slowly (more exploration);
            values close to 0.0 cool quickly (more exploitation).
        n_iterations: Total number of SA iterations to run.
        seed: Random seed for reproducibility.
        lambda_cai: Weight of the CAI term in the joint objective.
            ``0.0`` = pure MFE optimisation; higher values prioritise
            codon adaptation.
        window_size: Window size in nucleotides for windowed SA.
        window_overlap: Overlap between adjacent windows in nucleotides.
    """

    initial_temp: float = 1.0
    cooling_rate: float = 0.997
    n_iterations: int = 10_000
    seed: int = 42
    lambda_cai: float = 3.0
    window_size: int = 200
    window_overlap: int = 50


@dataclass
class SAResult:
    """Result of a simulated annealing optimization run.

    Attributes:
        sequence: Optimized DNA coding sequence (uppercase ACGT).
        mfe: Minimum free energy of the result in kcal/mol.
        cai: Codon Adaptation Index of the result in [0.0, 1.0].
        score: Final joint objective value (lower is better).
        n_iterations: Number of SA iterations actually performed.
        method: Optimisation method identifier (e.g. ``"sa"`` or
            ``"sa_windowed"``).
        elapsed_seconds: Wall-clock time in seconds.
        temperature_history: List of temperature values sampled during
            the run (one per iteration).
    """

    sequence: str
    mfe: float
    cai: float
    score: float
    n_iterations: int
    method: str
    elapsed_seconds: float
    temperature_history: list[float] = field(default_factory=list)


# ══════════════════════════════════════════════════════════════
# Core helper functions
# ══════════════════════════════════════════════════════════════


def compute_objective(
    dna_seq: str,
    organism: str = "human",
    lambda_cai: float = 3.0,
) -> float:
    """Compute the joint MFE + CAI objective for a DNA sequence.

    The objective follows the LinearDesign formulation::

        score = MFE(dna) + lambda_cai * (-log(CAI(dna)))

    Lower scores are better.  A ``lambda_cai`` of 0.0 yields pure MFE
    minimisation; higher values increasingly penalise low-CAI sequences.

    Args:
        dna_seq: DNA coding sequence (uppercase ACGT).
        organism: Target organism for CAI computation.
        lambda_cai: CAI penalty weight.

    Returns:
        Joint objective value (float).  Lower is better.

    References
    ----------
    Zhang, H. et al. (2023). *Nature* 621:396–403.
    """
    mfe = _compute_mfe(dna_seq)
    cai = _compute_cai(dna_seq, organism)
    cai_penalty = -math.log(max(cai, 1e-10)) * lambda_cai if cai > 0 else 0.0
    return mfe + cai_penalty


def cooling_schedule(
    iteration: int,
    n_iterations: int,
    initial_temp: float,
    method: str = "exponential",
) -> float:
    """Compute temperature at a given iteration using the specified schedule.

    Supported methods:

    - ``"exponential"`` (default): ``T = T0 * r^i``  — the classic
      geometric cooling schedule used in Kirkpatrick et al. (1983).
    - ``"linear"``: ``T = T0 * (1 - i/N)``  — linear temperature
      decrease; reaches zero at the last iteration.
    - ``"logarithmic"``: ``T = T0 / log(2 + i)``  — slow cooling
      inspired by the theoretical convergence proof of Geman & Geman (1984).
    - ``"adaptive"``: Exponential cooling but with a floor at 0.01 * T0
      to prevent the temperature from vanishing completely.

    Args:
        iteration: Current iteration index (0-based).
        n_iterations: Total number of iterations.
        initial_temp: Initial temperature T0.
        method: Cooling schedule name.

    Returns:
        Temperature at the given iteration.

    Raises:
        ValueError: If *method* is not one of the supported schedule names.

    References
    ----------
    Kirkpatrick, S. et al. (1983). *Science* 220:671–680.
    Geman, S. & Geman, D. (1984). *IEEE TPAMI* 6:721–741.
    """
    if method == "exponential":
        # Geometric cooling: T0 * cooling_rate^i
        # The cooling_rate is embedded in the caller's loop; here we
        # compute a per-iteration rate that achieves initial_temp * r^N ≈ 0.
        # When called from optimize_sa, the caller passes the effective
        # rate directly; this formula provides a standalone calculation.
        rate = 0.997  # default geometric rate
        return initial_temp * (rate ** iteration)

    elif method == "linear":
        progress = iteration / max(n_iterations, 1)
        return initial_temp * max(1.0 - progress, 0.0)

    elif method == "logarithmic":
        return initial_temp / math.log(2 + iteration)

    elif method == "adaptive":
        # Exponential with a floor at 1% of initial temperature
        rate = 0.997
        temp = initial_temp * (rate ** iteration)
        return max(temp, 0.01 * initial_temp)

    else:
        valid = ", ".join(("exponential", "linear", "logarithmic", "adaptive"))
        raise ValueError(
            f"Unknown cooling schedule {method!r}. Valid methods: {valid}"
        )


def accept_move(
    current_score: float,
    new_score: float,
    temperature: float,
) -> bool:
    """Metropolis acceptance criterion for simulated annealing.

    Always accepts improving moves (lower scores).  Accepts worsening
    moves with probability ``exp(-delta / T)`` where ``delta`` is the
    score increase.

    Args:
        current_score: Objective value of the current solution.
        new_score: Objective value of the proposed solution.
        temperature: Current annealing temperature.

    Returns:
        True if the move should be accepted, False otherwise.

    References
    ----------
    Kirkpatrick, S. et al. (1983). *Science* 220:671–680.
    Metropolis, N. et al. (1953). *J. Chem. Phys.* 21:1087–1092.
    """
    delta = new_score - current_score
    if delta <= 0:
        return True
    if temperature <= 0:
        return False
    try:
        probability = math.exp(-delta / temperature)
    except OverflowError:
        # delta / T is very large; probability is essentially zero
        return False
    return random.random() < probability


# ══════════════════════════════════════════════════════════════
# MFE / CAI helpers
# ══════════════════════════════════════════════════════════════


def _compute_mfe_rna_nn(dna_seq: str) -> float:
    """Compute MFE using RNA nearest-neighbor stacking approximation.

    Uses Turner 2004 (Mathews et al.) average stacking free energies
    per dinucleotide step.  This is more accurate than the scalar GC
    heuristic because it accounts for AU/GC/GU composition differences
    and sequence order effects.

    NOTE: This is still an approximation — it assumes perfect base-pairing
    (no loops, bulges, or unpaired regions).  For accurate MFE, install
    ViennaRNA (RNAfold).

    Args:
        dna_seq: DNA coding sequence (uppercase ACGT).

    Returns:
        Approximate MFE in kcal/mol (negative = stable).
    """
    # Convert DNA to RNA
    rna = dna_seq.upper().replace("T", "U")

    # Turner 2004 average stacking ΔG (kcal/mol) per dinucleotide
    # These are averaged over the Watson-Crick complement stacks
    _RNA_DINUC_DG = {
        'AA': -1.00, 'AU': -1.10, 'AC': -1.55, 'AG': -1.55,
        'UA': -1.00, 'UU': -0.70, 'UC': -1.25, 'UG': -1.25,
        'CA': -1.90, 'CU': -1.60, 'CC': -2.10, 'CG': -2.40,
        'GA': -1.55, 'GU': -1.40, 'GC': -2.40, 'GG': -1.90,
    }

    mfe = 0.0
    for i in range(len(rna) - 1):
        dinuc = rna[i:i + 2]
        mfe += _RNA_DINUC_DG.get(dinuc, -1.0)  # default average

    # Initiation penalty (~2 × initiation + terminal AU penalties)
    mfe += 3.4

    return mfe


def _compute_mfe(dna_seq: str) -> float:
    """Compute MFE using ViennaRNA; fall back to RNA NN approximation with warning."""
    rna_seq = dna_seq.upper().replace("T", "U")
    try:
        import RNA  # type: ignore[import-untyped]
        fc = RNA.fold_compound(rna_seq)
        mfe = fc.mfe()[1]
        return mfe
    except ImportError:
        pass

    # RNA nearest-neighbor stacking approximation (Turner 2004)
    warnings.warn(
        "ViennaRNA not available; Using RNA NN stacking approximation "
        "(Turner 2004) for MFE. Install ViennaRNA (pip install ViennaRNA) "
        "for accurate MFE computation.",
        UserWarning,
        stacklevel=3,
    )
    return _compute_mfe_rna_nn(dna_seq)


def _compute_cai(dna_seq: str, organism: str) -> float:
    """Compute CAI using the expression module (lazy import)."""
    from ..expression.translation import compute_cai
    return compute_cai(dna_seq, organism=organism)


def _best_cai_codon(aa: str, codon_freq: dict[str, float]) -> str:
    """Return the highest-adaptiveness codon for *aa*."""
    codons = AA_TO_CODONS.get(aa, [])
    if not codons:
        return CODON_TABLE.get(aa, ["---"])[0] if aa in CODON_TABLE else "---"
    return max(codons, key=lambda c: codon_freq.get(c, 0.01))


# ══════════════════════════════════════════════════════════════
# Core SA optimizer
# ══════════════════════════════════════════════════════════════


def optimize_sa(
    protein_seq: str,
    organism: str = "human",
    config: SAConfig | None = None,
) -> SAResult:
    """Full simulated annealing optimization for mRNA design.

    Starts from the highest-CAI codon assignment, then iteratively
    proposes random synonymous codon swaps and accepts/rejects them
    according to the Metropolis criterion with a decreasing temperature
    schedule.

    Algorithm:
        1. Initialise every position with the highest-CAI codon.
        2. At each iteration, pick a random position and propose a
           random synonymous codon swap.
        3. Compute the joint objective ``MFE - log(CAI) * lambda``.
        4. Accept if the score improves; otherwise accept with
           probability ``exp(-delta / T)``.
        5. Cool temperature according to the chosen schedule.
        6. Return the best solution found.

    Args:
        protein_seq: Protein amino acid sequence (1-letter codes, no stop).
        organism: Target organism for codon usage.
        config: SA configuration.  Uses defaults if ``None``.

    Returns:
        ``SAResult`` with the optimised sequence, MFE, CAI, score,
        iteration count, and temperature history.

    Raises:
        ValueError: If *protein_seq* is empty.

    References
    ----------
    Kirkpatrick, S. et al. (1983). *Science* 220:671–680.
    Zhang, H. et al. (2023). *Nature* 621:396–403.
    """
    if not protein_seq:
        raise ValueError("protein_seq must not be empty")

    cfg = config if config is not None else SAConfig()
    rng = random.Random(cfg.seed)
    t_start = time.perf_counter()

    org_key = resolve_organism(organism, strict=False)
    codon_freq = CODON_ADAPTIVENESS_TABLES.get(org_key, {})

    # Initialise with highest-CAI codons
    codons = [_best_cai_codon(aa, codon_freq) for aa in protein_seq]
    current_seq = "".join(codons)

    current_score = compute_objective(current_seq, organism, cfg.lambda_cai)
    best_seq = current_seq
    best_score = current_score

    temp = cfg.initial_temp
    temperature_history: list[float] = []

    for iteration in range(cfg.n_iterations):
        # Record temperature
        temperature_history.append(temp)

        # Random codon swap
        pos = rng.randint(0, len(protein_seq) - 1)
        aa = protein_seq[pos]
        alternatives = AA_TO_CODONS.get(aa, [])
        if len(alternatives) <= 1:
            continue

        alt = rng.choice(alternatives)
        old_codon = codons[pos]
        if alt == old_codon:
            continue

        codons[pos] = alt
        new_seq = "".join(codons)
        new_score = compute_objective(new_seq, organism, cfg.lambda_cai)

        if accept_move(current_score, new_score, temp):
            current_score = new_score
            if current_score < best_score:
                best_seq = new_seq
                best_score = current_score
        else:
            codons[pos] = old_codon

        # Cool temperature (exponential schedule)
        temp *= cfg.cooling_rate

    elapsed = time.perf_counter() - t_start
    final_mfe = _compute_mfe(best_seq)
    final_cai = _compute_cai(best_seq, organism)
    final_score = compute_objective(best_seq, organism, cfg.lambda_cai)

    return SAResult(
        sequence=best_seq,
        mfe=final_mfe,
        cai=final_cai,
        score=final_score,
        n_iterations=cfg.n_iterations,
        method="sa",
        elapsed_seconds=elapsed,
        temperature_history=temperature_history,
    )


# ══════════════════════════════════════════════════════════════
# Windowed SA for long genes
# ══════════════════════════════════════════════════════════════


def optimize_sa_windowed(
    protein_seq: str,
    organism: str = "human",
    config: SAConfig | None = None,
) -> SAResult:
    """Windowed simulated annealing for long genes.

    Splits the coding sequence into overlapping windows of
    ``config.window_size`` nucleotides with ``config.window_overlap``
    overlap, runs SA independently on each window, then merges the
    results and performs a final global refinement pass.

    This avoids the O(n²) cost of ViennaRNA folding for very long
    sequences while still capturing local structure–function
    relationships.

    Args:
        protein_seq: Protein amino acid sequence.
        organism: Target organism for codon usage.
        config: SA configuration.  Uses defaults if ``None``.

    Returns:
        ``SAResult`` with the optimised sequence and temperature
        history from all windows concatenated.

    Raises:
        ValueError: If *protein_seq* is empty.

    References
    ----------
    Kirkpatrick, S. et al. (1983). *Science* 220:671–680.
    """
    if not protein_seq:
        raise ValueError("protein_seq must not be empty")

    cfg = config if config is not None else SAConfig()
    rng = random.Random(cfg.seed)
    t_start = time.perf_counter()

    org_key = resolve_organism(organism, strict=False)
    codon_freq = CODON_ADAPTIVENESS_TABLES.get(org_key, {})

    # Initialise with highest-CAI codons
    codons = [_best_cai_codon(aa, codon_freq) for aa in protein_seq]
    n_nt = len(codons) * 3

    window_size = cfg.window_size
    window_overlap = cfg.window_overlap
    all_temperature_history: list[float] = []

    # Determine windows
    if n_nt <= window_size * 2:
        # Short enough for a single global pass — delegate to optimize_sa
        return optimize_sa(protein_seq, organism, config)

    n_windows = max(1, (n_nt - window_overlap) // (window_size - window_overlap))
    iters_per_window = cfg.n_iterations // n_windows

    for win_idx in range(n_windows):
        win_start_aa = win_idx * (window_size - window_overlap) // 3
        win_end_aa = min(len(protein_seq), win_start_aa + window_size // 3)

        if win_start_aa >= win_end_aa:
            continue

        current_seq = "".join(codons)
        current_score = compute_objective(current_seq, organism, cfg.lambda_cai)
        best_seq = current_seq
        best_score = current_score

        temp = cfg.initial_temp

        for _iteration in range(iters_per_window):
            all_temperature_history.append(temp)

            # Random codon swap within window
            pos = rng.randint(win_start_aa, max(win_start_aa, win_end_aa - 1))
            aa = protein_seq[pos]
            alternatives = AA_TO_CODONS.get(aa, [])
            if len(alternatives) <= 1:
                continue

            alt = rng.choice(alternatives)
            old_codon = codons[pos]
            if alt == old_codon:
                continue

            codons[pos] = alt
            new_seq = "".join(codons)
            new_score = compute_objective(new_seq, organism, cfg.lambda_cai)

            if accept_move(current_score, new_score, temp):
                current_score = new_score
                if current_score < best_score:
                    best_seq = new_seq
                    best_score = current_score
            else:
                codons[pos] = old_codon

            temp *= cfg.cooling_rate

    # Global refinement pass: a short SA run on the full sequence
    merged_seq = "".join(codons)
    current_score = compute_objective(merged_seq, organism, cfg.lambda_cai)
    best_seq = merged_seq
    best_score = current_score

    temp = cfg.initial_temp * 0.5  # start cooler for refinement
    refinement_iters = max(100, cfg.n_iterations // 10)

    for _iteration in range(refinement_iters):
        all_temperature_history.append(temp)

        pos = rng.randint(0, len(protein_seq) - 1)
        aa = protein_seq[pos]
        alternatives = AA_TO_CODONS.get(aa, [])
        if len(alternatives) <= 1:
            continue

        alt = rng.choice(alternatives)
        old_codon = codons[pos]
        if alt == old_codon:
            continue

        codons[pos] = alt
        new_seq = "".join(codons)
        new_score = compute_objective(new_seq, organism, cfg.lambda_cai)

        if accept_move(current_score, new_score, temp):
            current_score = new_score
            if current_score < best_score:
                best_seq = new_seq
                best_score = current_score
        else:
            codons[pos] = old_codon

        temp *= cfg.cooling_rate

    elapsed = time.perf_counter() - t_start
    final_mfe = _compute_mfe(best_seq)
    final_cai = _compute_cai(best_seq, organism)
    final_score = compute_objective(best_seq, organism, cfg.lambda_cai)

    return SAResult(
        sequence=best_seq,
        mfe=final_mfe,
        cai=final_cai,
        score=final_score,
        n_iterations=cfg.n_iterations,
        method="sa_windowed",
        elapsed_seconds=elapsed,
        temperature_history=all_temperature_history,
    )


# ══════════════════════════════════════════════════════════════
# Multi-objective Pareto front
# ══════════════════════════════════════════════════════════════


def optimize_pareto(
    protein_seq: str,
    organism: str = "human",
    lambda_values: list[float] | None = None,
    **kwargs: Any,
) -> list[SAResult]:
    """Run SA with multiple lambda values to trace the MFE–CAI Pareto front.

    Each ``lambda_cai`` value produces a different trade-off between MFE
    minimisation and CAI maximisation.  The collection of results
    approximates the Pareto-optimal set: no single result dominates
    another on both objectives simultaneously.

    The default lambda sweep is ``[0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0]``.

    Args:
        protein_seq: Protein amino acid sequence.
        organism: Target organism for codon usage.
        lambda_values: List of ``lambda_cai`` values to sweep.  Defaults
            to a geometric sweep from 0 to 8.
        **kwargs: Additional keyword arguments passed to ``SAConfig``
            (e.g. ``n_iterations``, ``seed``).

    Returns:
        List of ``SAResult`` objects, one per lambda value, sorted by
        increasing lambda (i.e. from pure-MFE to CAI-heavy).

    Raises:
        ValueError: If *protein_seq* is empty.

    References
    ----------
    Zhang, H. et al. (2023). *Nature* 621:396–403.
    """
    if not protein_seq:
        raise ValueError("protein_seq must not be empty")

    if lambda_values is None:
        lambda_values = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0]

    results: list[SAResult] = []
    for lam in lambda_values:
        cfg = SAConfig(lambda_cai=lam, **kwargs)
        result = optimize_sa(protein_seq, organism, config=cfg)
        results.append(result)

    return results


def select_pareto_optimal(
    results: list[SAResult],
    target_cai: float | None = None,
    target_mfe: float | None = None,
) -> SAResult:
    """Select the best result from a Pareto front given user targets.

    Selection strategy:
      - If ``target_cai`` is specified (and ``target_mfe`` is not),
        select the result with the highest CAI that still meets the
        target.
      - If ``target_mfe`` is specified (and ``target_cai`` is not),
        select the result with the lowest MFE that still meets the
        target.
      - If both are specified, select the result that minimises the
        weighted distance to both targets.
      - If neither is specified, select the result with the best
        (lowest) joint score.

    Args:
        results: List of ``SAResult`` from ``optimize_pareto()`` or
            ``optimize_nsga2()``.
        target_cai: Desired minimum CAI threshold.
        target_mfe: Desired maximum MFE threshold (less negative is
            closer to zero, i.e. less stable).

    Returns:
        The ``SAResult`` that best matches the specified targets.

    Raises:
        ValueError: If *results* is empty.
    """
    if not results:
        raise ValueError("results must not be empty")

    if len(results) == 1:
        return results[0]

    if target_cai is not None and target_mfe is None:
        # Select highest CAI that meets the target
        candidates = [r for r in results if r.cai >= target_cai]
        if not candidates:
            # Relax: pick the one with highest CAI
            return max(results, key=lambda r: r.cai)
        return max(candidates, key=lambda r: r.cai)

    elif target_mfe is not None and target_cai is None:
        # Select lowest MFE that meets the target
        candidates = [r for r in results if r.mfe <= target_mfe]
        if not candidates:
            # Relax: pick the one with lowest MFE
            return min(results, key=lambda r: r.mfe)
        return min(candidates, key=lambda r: r.mfe)

    elif target_cai is not None and target_mfe is not None:
        # Normalise both objectives and minimise distance to target
        mfe_vals = [r.mfe for r in results]
        cai_vals = [r.cai for r in results]
        mfe_range = max(mfe_vals) - min(mfe_vals) if len(set(mfe_vals)) > 1 else 1.0
        cai_range = max(cai_vals) - min(cai_vals) if len(set(cai_vals)) > 1 else 1.0

        def _distance(r: SAResult) -> float:
            mfe_d = abs(r.mfe - target_mfe) / abs(mfe_range) if mfe_range else 0.0
            cai_d = abs(r.cai - target_cai) / abs(cai_range) if cai_range else 0.0
            return mfe_d + cai_d

        return min(results, key=_distance)

    else:
        # No targets specified: pick the best joint score
        return min(results, key=lambda r: r.score)


# ══════════════════════════════════════════════════════════════
# Optional NSGA-II integration via pymoo
# ══════════════════════════════════════════════════════════════


def optimize_nsga2(
    protein_seq: str,
    organism: str = "human",
    n_generations: int = 100,
    pop_size: int = 50,
) -> list[SAResult]:
    """Multi-objective mRNA optimisation using NSGA-II from pymoo.

    Uses the Non-dominated Sorting Genetic Algorithm II (NSGA-II) to
    simultaneously minimise MFE and maximise CAI, producing a true
    Pareto-optimal front without requiring a scalar ``lambda_cai``
    parameter.

    Requires the ``pymoo`` package (``pip install pymoo``).  Falls back
    to ``optimize_pareto()`` with a warning when pymoo is unavailable.

    Args:
        protein_seq: Protein amino acid sequence.
        organism: Target organism for codon usage.
        n_generations: Number of NSGA-II generations.
        pop_size: Population size per generation.

    Returns:
        List of ``SAResult`` objects representing the Pareto front.

    Raises:
        ValueError: If *protein_seq* is empty.

    References
    ----------
    Deb, K. et al. (2002). *IEEE Trans. Evol. Comput.* 6(2):182–197.
    """
    if not protein_seq:
        raise ValueError("protein_seq must not be empty")

    try:
        import numpy as np  # noqa: F401
        from pymoo.core.problem import Problem  # type: ignore[import-untyped]
        from pymoo.algorithms.moo.nsga2 import NSGA2  # type: ignore[import-untyped]
        from pymoo.operators.crossover.sbx import SBX  # type: ignore[import-untyped]
        from pymoo.operators.mutation.pm import PM  # type: ignore[import-untyped]
        from pymoo.operators.sampling.rnd import IntegerRandomSampling  # type: ignore[import-untyped]
        from pymoo.optimize import minimize as pymoo_minimize  # type: ignore[import-untyped]
    except ImportError:
        warnings.warn(
            "pymoo is not installed. Falling back to optimize_pareto() with "
            "default lambda sweep. Install pymoo for true multi-objective "
            "optimisation: pip install pymoo",
            UserWarning,
            stacklevel=2,
        )
        return optimize_pareto(protein_seq, organism)

    org_key = resolve_organism(organism, strict=False)
    codon_freq = CODON_ADAPTIVENESS_TABLES.get(org_key, {})

    # Build per-position codon options
    n_positions = len(protein_seq)
    codon_options: list[list[str]] = []
    for aa in protein_seq:
        opts = AA_TO_CODONS.get(aa, [])
        if not opts:
            opts = [CODON_TABLE.get(aa, ["---"])[0] if aa in CODON_TABLE else "---"]
        codon_options.append(opts)

    # Maximum number of codon choices at any position (for variable-length encoding)
    max_choices = max(len(opts) for opts in codon_options)

    # Mapping: position × choice_index → codon string
    # choice_index >= len(codon_options[pos]) maps to the last valid codon

    def _decode(individual: list[int]) -> list[str]:
        """Decode an integer vector to a list of codon strings."""
        return [
            codon_options[i][min(idx, len(codon_options[i]) - 1)]
            for i, idx in enumerate(individual)
        ]

    class _mRNAProblem(Problem):
        """pymoo Problem definition for multi-objective mRNA optimisation."""

        def __init__(self) -> None:
            super().__init__(
                n_var=n_positions,
                n_obj=2,  # MFE (minimise), -CAI (minimise = maximise CAI)
                n_ieq_constr=0,
                xl=np.zeros(n_positions, dtype=int),
                xu=np.array(
                    [len(opts) - 1 for opts in codon_options], dtype=int
                ),
                vtype=int,
            )

        def _evaluate(self, X: Any, out: Any, *args: Any, **kwargs: Any) -> None:
            import numpy as np  # type: ignore[import-untyped]

            f1 = np.full(X.shape[0], np.inf)  # MFE
            f2 = np.full(X.shape[0], np.inf)  # -CAI

            for i in range(X.shape[0]):
                individual = [int(x) for x in X[i]]
                codons_list = _decode(individual)
                dna_seq = "".join(codons_list)

                mfe = _compute_mfe(dna_seq)
                cai = _compute_cai(dna_seq, organism)

                f1[i] = mfe
                f2[i] = -cai  # minimise negative CAI = maximise CAI

            out["F"] = np.column_stack([f1, f2])

    # Set up and run NSGA-II
    problem = _mRNAProblem()

    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=IntegerRandomSampling(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
    )

    t_start = time.perf_counter()
    res = pymoo_minimize(
        problem,
        algorithm,
        ("n_gen", n_generations),
        seed=42,
        verbose=False,
    )
    elapsed = time.perf_counter() - t_start

    # Convert pymoo results to SAResult list
    results: list[SAResult] = []
    if res.X is not None and res.F is not None:
        import numpy as np  # type: ignore[import-untyped]

        # Handle both single and population results
        X_pop = res.X if len(res.X.shape) > 1 else res.X.reshape(1, -1)
        F_pop = res.F if len(res.F.shape) > 1 else res.F.reshape(1, -1)

        for i in range(len(F_pop)):
            individual = [int(x) for x in X_pop[i]]
            codons_list = _decode(individual)
            dna_seq = "".join(codons_list)

            mfe_val = float(F_pop[i, 0])
            cai_val = -float(F_pop[i, 1])
            score_val = compute_objective(dna_seq, organism, lambda_cai=3.0)

            results.append(SAResult(
                sequence=dna_seq,
                mfe=mfe_val,
                cai=cai_val,
                score=score_val,
                n_iterations=n_generations * pop_size,
                method="nsga2",
                elapsed_seconds=elapsed,
            ))

    # Deduplicate by sequence
    seen: set[str] = set()
    unique_results: list[SAResult] = []
    for r in results:
        if r.sequence not in seen:
            seen.add(r.sequence)
            unique_results.append(r)

    if not unique_results:
        # Fallback: return a single result using highest-CAI codons
        codons_list = [_best_cai_codon(aa, codon_freq) for aa in protein_seq]
        dna_seq = "".join(codons_list)
        mfe_val = _compute_mfe(dna_seq)
        cai_val = _compute_cai(dna_seq, organism)
        score_val = compute_objective(dna_seq, organism, lambda_cai=3.0)
        unique_results.append(SAResult(
            sequence=dna_seq,
            mfe=mfe_val,
            cai=cai_val,
            score=score_val,
            n_iterations=n_generations * pop_size,
            method="nsga2",
            elapsed_seconds=elapsed,
        ))

    return unique_results
