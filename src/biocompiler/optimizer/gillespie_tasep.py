"""
BioCompiler Gillespie TASEP — Standalone Ribosome Traffic Simulation
=====================================================================

A self-contained, dataclass-driven Gillespie TASEP implementation for
simulating ribosome traffic on mRNA transcripts.  This module is designed
to be imported directly without pulling in the full optimizer pipeline.

The Gillespie Direct Method (Gillespie 1977) is preferred over
discrete-time TASEP because it samples exact reaction times from
exponential distributions, avoiding the lattice-synchronization artefacts
of discrete-time updates and naturally handling rare events (stalls,
collisions) without fine time-step tuning.

Algorithm Overview
------------------
1. Compute all possible reaction rates (initiation + per-ribosome
   elongation hops) at the current state.
2. Draw waiting time Δt ~ Exp(Σ rates).
3. Select a reaction proportionally to its rate.
4. Execute the reaction, enforcing the simple exclusion rule (two
   ribosomes cannot overlap within *ribosome_footprint* codons).
5. Advance the clock by Δt and repeat until *max_time*.

Key Design Decisions
--------------------
- **Ribosome positions are stored in codon coordinates** (0-indexed).
- **Exclusion is checked via pairwise distance**; a ribosome at position
  *p* blocks positions *[p − footprint + 1, p + footprint − 1]* for
  other ribosomes.
- **Termination** at the last codon uses a fast rate (2 × elongation_rate)
  following Shah et al. (2013).
- **Stall-site detection** uses a 3× mean-density threshold, consistent
  with ribosome-profiling heuristics.

References
----------
- MacDonald, C.T. & Gibbs, J.H. (1968) "Kinetics of biopolymerization
  on nucleic acid templates." *Biopolymers* 6(1):1-25.
- Ciandrini, L., Stansfield, I. & Romano, M.C. (2013) "Ribosome
  traffic on mRNAs maps to gene orientation: roles of the 5' UTR and
  the coding sequence." *PLoS Computational Biology* 9(8):e1002866.
- Gillespie, D.T. (1977) "Exact stochastic simulation of coupled
  chemical reactions." *The Journal of Chemical Physics* 81(25):2340.
- Shah, P., Ding, Y., Niemczyk, M., Kudla, G. & Plotkin, J.B. (2013)
  "Rate-limiting steps in yeast protein translation." *Cell*
  153(7):1589-1601.
"""

from __future__ import annotations

import logging
import math
import random
from dataclasses import dataclass, field
from typing import Any

__all__ = [
    # Dataclasses
    "RibosomeState",
    "TASEPResult",
    # Core simulation
    "simulate_gillespie",
    # Ensemble runner
    "run_ensemble",
    # Initiation rate estimation
    "estimate_initiation_rate",
    # Dwell time computation
    "compute_dwell_times",
    # Optional Pinetree integration
    "simulate_pinetree",
]

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Dataclasses
# ═══════════════════════════════════════════════════════════════════════════════


@dataclass
class RibosomeState:
    """Snapshot of ribosome positions on the mRNA at a given time.

    Attributes:
        positions: Codon positions (0-indexed) of each ribosome on the
            transcript, sorted in ascending order.
        time: Simulation clock value (seconds) at this snapshot.
    """

    positions: list[int] = field(default_factory=list)
    time: float = 0.0


@dataclass
class TASEPResult:
    """Result of a single Gillespie TASEP simulation run.

    Attributes:
        codon_density: Per-codon occupancy probability (dimensionless,
            0.0–1.0).  Length equals the number of codons in the ORF.
        stall_sites: Indices of codons whose density exceeds 3× the
            mean, indicating putative ribosome stalls.
        collision_count: Total number of blocked elongation attempts
            (i.e. a ribosome wanted to hop forward but the exclusion
            rule prevented it).
        mean_ribosome_count: Average number of ribosomes on the mRNA
            over the course of the simulation.
        time_elapsed: Total simulation time (seconds) that elapsed.
        method: Simulation method identifier (always ``"gillespie"``).
    """

    codon_density: list[float] = field(default_factory=list)
    stall_sites: list[int] = field(default_factory=list)
    collision_count: int = 0
    mean_ribosome_count: float = 0.0
    time_elapsed: float = 0.0
    method: str = "gillespie"


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Core Gillespie TASEP Simulation
# ═══════════════════════════════════════════════════════════════════════════════


def simulate_gillespie(
    dwell_times: list[float],
    initiation_rate: float = 0.1,
    elongation_rate: float = 10.0,
    ribosome_footprint: int = 10,
    max_time: float = 1000.0,
    seed: int = 42,
) -> TASEPResult:
    """Simulate ribosome traffic using the continuous-time Gillespie TASEP.

    Implements the Gillespie Direct Method (Gillespie 1977) applied to
    the Totally Asymmetric Simple Exclusion Process (MacDonald & Gibbs
    1968) with codon-specific elongation rates derived from dwell times.

    Args:
        dwell_times: Per-codon dwell times in **seconds**.  A larger
            value means a slower codon (e.g. rare-codon stall).  Length
            must equal the number of codons in the ORF.
        initiation_rate: Rate at which new ribosomes load onto the 5'
            end of the mRNA (ribosomes/s).  Default 0.1.
        elongation_rate: Baseline elongation rate (codons/s) used to
            scale the per-codon hop rates.  Default 10.0.
        ribosome_footprint: Number of codons occupied by a single
            ribosome (exclusion radius).  Default 10.
        max_time: Maximum simulation wall-clock in seconds.  Default
            1000.0.
        seed: Random seed for reproducibility.  Default 42.

    Returns:
        A :class:`TASEPResult` dataclass with codon densities, stall
        sites, collision count, mean ribosome count, and elapsed time.

    Raises:
        ValueError: If *dwell_times* is empty.
    """
    if not dwell_times:
        raise ValueError("dwell_times must be a non-empty list of per-codon dwell times")

    rng = random.Random(seed)
    n_codons = len(dwell_times)

    # Convert dwell times → hop rates.  Longer dwell → slower hop.
    # Floor at 0.001 s to avoid division-by-zero.
    hop_rates = [1.0 / max(dt, 0.001) for dt in dwell_times]

    # Normalise so that the fastest codon has rate = elongation_rate.
    max_rate = max(hop_rates)
    if max_rate > 0:
        hop_rates = [r / max_rate * elongation_rate for r in hop_rates]

    # ── State ──────────────────────────────────────────────────────────────
    ribosomes: list[int] = []  # positions in codon coordinates, ascending
    time = 0.0

    # ── Accumulators ───────────────────────────────────────────────────────
    codon_occupancy = [0] * n_codons
    collision_count = 0
    sample_count = 0
    ribosome_count_sum = 0

    # ── Main loop ──────────────────────────────────────────────────────────
    while time < max_time:
        rates: list[float] = []
        events: list[tuple[str, int]] = []  # (event_type, ribosome_index)

        # --- Initiation event ---
        # A new ribosome can load if position 0 is free (no ribosome
        # within `ribosome_footprint` of the start).
        can_initiate = True
        for pos in ribosomes:
            if pos < ribosome_footprint:
                can_initiate = False
                break
        if can_initiate:
            rates.append(initiation_rate)
            events.append(("initiate", -1))

        # --- Elongation / termination events ---
        for idx, pos in enumerate(ribosomes):
            if pos < n_codons - 1:
                # Attempt to hop forward by 1 codon.
                next_pos = pos + 1
                can_hop = True
                for other_idx, other_pos in enumerate(ribosomes):
                    if other_idx != idx and abs(other_pos - next_pos) < ribosome_footprint:
                        can_hop = False
                        collision_count += 1
                        break
                if can_hop:
                    rates.append(hop_rates[min(pos, n_codons - 1)])
                    events.append(("hop", idx))
            else:
                # At last codon → termination (fast off-rate).
                rates.append(elongation_rate * 2.0)
                events.append(("terminate", idx))

        if not rates:
            break

        # 2. Draw waiting time from exponential(total_rate).
        total_rate = sum(rates)
        if total_rate <= 0:
            break
        dt = rng.expovariate(total_rate)
        time += dt

        if time > max_time:
            break

        # 3. Select reaction proportionally.
        r = rng.random() * total_rate
        cumsum = 0.0
        selected = 0
        for i, rate in enumerate(rates):
            cumsum += rate
            if r <= cumsum:
                selected = i
                break

        # 4. Execute selected reaction.
        event_type, event_idx = events[selected]

        if event_type == "initiate":
            ribosomes.insert(0, 0)
        elif event_type == "hop":
            ribosomes[event_idx] += 1
        elif event_type == "terminate":
            ribosomes.pop(event_idx)

        # 5. Record occupancy snapshot.
        sample_count += 1
        ribosome_count_sum += len(ribosomes)
        for pos in ribosomes:
            if 0 <= pos < n_codons:
                codon_occupancy[pos] += 1

    # ── Compute final statistics ───────────────────────────────────────────
    if sample_count > 0:
        density = [o / sample_count for o in codon_occupancy]
    else:
        density = [0.0] * n_codons

    mean_density = sum(density) / max(1, len(density))
    stall_sites = [i for i, d in enumerate(density) if d > 3 * mean_density and mean_density > 0]

    return TASEPResult(
        codon_density=density,
        stall_sites=stall_sites,
        collision_count=collision_count,
        mean_ribosome_count=ribosome_count_sum / max(1, sample_count),
        time_elapsed=time,
        method="gillespie",
    )


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Ensemble Runner
# ═══════════════════════════════════════════════════════════════════════════════


def run_ensemble(
    dwell_times: list[float],
    n_runs: int = 100,
    burnin_fraction: float = 0.1,
    confidence: float = 0.95,
    **kwargs: Any,
) -> dict[str, Any]:
    """Run multiple independent Gillespie TASEP simulations and compute
    ensemble statistics with confidence intervals.

    Each run uses a different random seed (base_seed + run_index × 1000).
    A configurable fraction of initial runs is discarded as burn-in.

    Args:
        dwell_times: Per-codon dwell times (seconds).
        n_runs: Number of independent simulation runs.  Default 100.
        burnin_fraction: Fraction of runs to discard from the start.
            Default 0.1 (10 %).
        confidence: Confidence level for the interval (e.g. 0.95).
            Default 0.95.
        **kwargs: Forwarded to :func:`simulate_gillespie`
            (e.g. ``initiation_rate``, ``elongation_rate``).

    Returns:
        A dict with the following keys:

        - ``mean_density`` (list[float]): Per-codon mean occupancy.
        - ``sem_density`` (list[float]): Standard error of the mean.
        - ``ci_lower`` / ``ci_upper`` (list[float]): Confidence-interval
          bounds at the requested level.
        - ``stall_sites`` (list[int]): Consensus stall sites (present
          in >50 % of post-burn-in runs).
        - ``mean_ribosome_count`` (float): Ensemble-average ribosome
          count.
        - ``mean_collision_count`` (float): Ensemble-average collisions.
        - ``n_runs_kept`` (int): Number of runs after burn-in discard.
        - ``confidence`` (float): Requested confidence level.
        - ``method`` (str): ``"gillespie_ensemble"``.
    """
    # Try to import scipy for accurate t-values; fall back otherwise.
    try:
        from scipy import stats as scipy_stats  # type: ignore[import-untyped]

        _has_scipy = True
    except ImportError:
        _has_scipy = False
        logger.warning("scipy not available; using approximate t-values for CI")

    # Try to import numpy for fast array operations; fall back otherwise.
    try:
        import numpy as np  # type: ignore[import-untyped]

        _has_numpy = True
    except ImportError:
        _has_numpy = False
        logger.info("numpy not available; using pure-Python statistics")

    base_seed = kwargs.pop("seed", 42)

    all_densities: list[list[float]] = []
    all_stalls: list[list[int]] = []
    all_collisions: list[int] = []
    all_ribosome_counts: list[float] = []

    for run_idx in range(n_runs):
        result = simulate_gillespie(
            dwell_times=dwell_times,
            seed=base_seed + run_idx * 1000,
            **kwargs,
        )
        all_densities.append(result.codon_density)
        all_stalls.append(result.stall_sites)
        all_collisions.append(result.collision_count)
        all_ribosome_counts.append(result.mean_ribosome_count)

    # Discard burn-in runs.
    n_burnin = int(n_runs * burnin_fraction)
    all_densities = all_densities[n_burnin:]
    all_stalls = all_stalls[n_burnin:]
    all_collisions = all_collisions[n_burnin:]
    all_ribosome_counts = all_ribosome_counts[n_burnin:]

    n_codons = len(dwell_times)
    n_kept = len(all_densities)

    # ── Mean ± SEM ─────────────────────────────────────────────────────────
    if _has_numpy:
        import numpy as _np  # type: ignore[import-untyped]

        density_matrix = _np.array(all_densities)  # shape (n_kept, n_codons)
        mean_density = density_matrix.mean(axis=0).tolist()
        sem_density = (density_matrix.std(axis=0, ddof=1) / math.sqrt(n_kept)).tolist() if n_kept > 1 else density_matrix.std(axis=0, ddof=0).tolist()
    else:
        mean_density = [0.0] * n_codons
        sem_density = [0.0] * n_codons
        for j in range(n_codons):
            col = [all_densities[i][j] for i in range(n_kept)]
            m = sum(col) / max(1, n_kept)
            mean_density[j] = m
            if n_kept > 1:
                variance = sum((x - m) ** 2 for x in col) / (n_kept - 1)
                sem_density[j] = math.sqrt(variance) / math.sqrt(n_kept)

    # ── Confidence interval ────────────────────────────────────────────────
    if _has_scipy and n_kept > 1:
        t_val = float(scipy_stats.t.ppf((1 + confidence) / 2, max(1, n_kept - 1)))
    else:
        # Approximate t-value.
        if n_kept > 30:
            t_val = 1.96  # Normal approximation
        elif n_kept > 1:
            t_val = 2.0 + 2.0 / max(1, n_kept)
        else:
            t_val = 1.96

    ci_lower = [m - t_val * s for m, s in zip(mean_density, sem_density)]
    ci_upper = [m + t_val * s for m, s in zip(mean_density, sem_density)]

    # ── Consensus stall sites ──────────────────────────────────────────────
    stall_counts = [0] * n_codons
    for stalls in all_stalls:
        for s in stalls:
            if 0 <= s < n_codons:
                stall_counts[s] += 1
    consensus_stalls = [i for i, c in enumerate(stall_counts) if c > n_kept * 0.5]

    # ── Aggregate statistics ───────────────────────────────────────────────
    mean_rib_count = sum(all_ribosome_counts) / max(1, n_kept)
    mean_collisions = sum(all_collisions) / max(1, n_kept)

    return {
        "mean_density": mean_density,
        "sem_density": sem_density,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "stall_sites": consensus_stalls,
        "mean_ribosome_count": mean_rib_count,
        "mean_collision_count": mean_collisions,
        "n_runs_kept": n_kept,
        "confidence": confidence,
        "method": "gillespie_ensemble",
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Initiation Rate Estimation
# ═══════════════════════════════════════════════════════════════════════════════

# Baseline initiation rates (ribosomes/s) from Shah et al. (2013) and
# Ciandrini et al. (2013).
_BASE_INITIATION_RATE = 0.1  # ribosomes/s (yeast / mammalian baseline)

# Kozak consensus strength multipliers (Noderer & Aldred 2024 heuristic).
# Strong Kozak (≥6/9 match): 1.5×, moderate (4-5/9): 1.0×, weak (<4/9): 0.5×
_KOZAK_MULTIPLIERS = {
    "strong": 1.5,
    "moderate": 1.0,
    "weak": 0.5,
}

# Shine-Dalgarno strength multipliers for prokaryotes.
_SD_MULTIPLIERS = {
    "strong": 2.0,   # exact complement (AGGAGG), ΔG ≈ −13 kcal/mol
    "moderate": 1.0,  # partial complement
    "weak": 0.4,      # no / very weak complement
}

# uORF penalty: each upstream open reading frame reduces initiation at
# the main ORF by ~30 % (Calvo et al. 2009 PNAS 106:7507).
_UORF_PENALTY_FACTOR = 0.70


def estimate_initiation_rate(
    kozak_score: float | None = None,
    sd_strength: float | None = None,
    uorf_count: int = 0,
) -> float:
    """Estimate the translation initiation rate from sequence features.

    The initiation rate is a key determinant of ribosome density on an
    mRNA.  This function combines three known regulators:

    1. **Kozak sequence strength** (eukaryotes): A stronger Kozak
       context increases the probability that the 43S pre-initiation
       complex successfully recruits the mRNA.
    2. **Shine-Dalgarno strength** (prokaryotes): Stronger mRNA-16S
       rRNA complementarity increases ribosome binding.
    3. **Upstream ORF count**: Each uORF competes for scanning
       ribosomes, reducing the flux that reaches the main ORF.

    At least one of *kozak_score* or *sd_strength* should be provided;
    otherwise the baseline rate is returned.

    Args:
        kozak_score: Kozak consensus match score in [0, 1].  Values
            ≥0.67 are "strong", 0.44–0.56 are "moderate", <0.44 are
            "weak".  Pass ``None`` to skip Kozak adjustment.
        sd_strength: Shine-Dalgarno binding strength in [0, 1].
            Values ≥0.7 are "strong", 0.3–0.7 are "moderate",
            <0.3 are "weak".  Pass ``None`` to skip SD adjustment.
        uorf_count: Number of upstream open reading frames.  Default 0.

    Returns:
        Estimated initiation rate in ribosomes/s.

    References:
        - Shah et al. (2013) *Cell* 153:1589-1601.
        - Calvo et al. (2009) *PNAS* 106:7507-7512.
        - Noderer & Aldred (2024) Quantitative models of translation
          initiation.
    """
    rate = _BASE_INITIATION_RATE

    # ── Kozak adjustment (eukaryotic) ──────────────────────────────────────
    if kozak_score is not None:
        if kozak_score >= 0.67:
            kozak_mult = _KOZAK_MULTIPLIERS["strong"]
        elif kozak_score >= 0.44:
            kozak_mult = _KOZAK_MULTIPLIERS["moderate"]
        else:
            kozak_mult = _KOZAK_MULTIPLIERS["weak"]
        rate *= kozak_mult

    # ── Shine-Dalgarno adjustment (prokaryotic) ───────────────────────────
    if sd_strength is not None:
        if sd_strength >= 0.7:
            sd_mult = _SD_MULTIPLIERS["strong"]
        elif sd_strength >= 0.3:
            sd_mult = _SD_MULTIPLIERS["moderate"]
        else:
            sd_mult = _SD_MULTIPLIERS["weak"]
        rate *= sd_mult

    # ── uORF penalty ──────────────────────────────────────────────────────
    if uorf_count > 0:
        rate *= _UORF_PENALTY_FACTOR ** uorf_count

    return rate


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Dwell-Time Computation
# ═══════════════════════════════════════════════════════════════════════════════

# Organism-specific baseline elongation rates (aa/s).
# Sources: Shah et al. 2013 (yeast), Ciandrini et al. 2013 (human/eukaryotes),
# Li et al. 2012 (E. coli).
_ELONGATION_RATES: dict[str, float] = {
    "human": 5.6,
    "mouse": 5.6,
    "cho": 5.6,
    "yeast": 9.3,
    "p_pastoris": 9.3,
    "e_coli": 15.0,
    "b_subtilis": 14.0,
    "c_elegans": 5.0,
    "d_melanogaster": 5.0,
    "a_thaliana": 5.0,
}

# Default organism key.
_DEFAULT_ORGANISM = "human"

# Arrhenius temperature correction factor (Q10 ≈ 2 for enzymatic rates).
_Q10 = 2.0
_REF_TEMP = 37.0  # °C


def compute_dwell_times(
    dna_seq: str,
    organism: str = "human",
    temperature: float = 37.0,
) -> list[float]:
    """Compute per-codon dwell times from a DNA coding sequence.

    Dwell times are inversely proportional to the tRNA Adaptation Index
    (tAI) weight for each codon.  Rare codons (low tAI) receive longer
    dwell times; optimal codons (high tAI) receive shorter dwell times.

    A temperature correction is applied using a Q10 = 2 Arrhenius model:
    ``rate(T) = rate(T_ref) × Q10^((T − T_ref) / 10)``.

    Args:
        dna_seq: DNA coding sequence (length must be a multiple of 3).
            Should not include the stop codon.
        organism: Organism key for tAI lookup.  Supported values:
            ``"human"``, ``"mouse"``, ``"cho"``, ``"yeast"``,
            ``"p_pastoris"``, ``"e_coli"``, ``"b_subtilis"``,
            ``"c_elegans"``, ``"d_melanogaster"``, ``"a_thaliana"``.
            Default ``"human"``.
        temperature: Incubation temperature in °C.  Default 37.0.

    Returns:
        List of per-codon dwell times in **seconds**.

    Raises:
        ValueError: If *dna_seq* length is not a multiple of 3 or
            contains invalid codons.
    """
    from ..organisms.tai_data import compute_tai_weights, SUPPORTED_ORGANISMS_TAI  # type: ignore[import-untyped]

    seq = dna_seq.upper().strip()
    if len(seq) == 0:
        raise ValueError("dna_seq must not be empty")
    if len(seq) % 3 != 0:
        raise ValueError(f"dna_seq length ({len(seq)}) is not a multiple of 3")

    n_codons = len(seq) // 3

    # Resolve organism key.
    org_key = organism.lower().replace(" ", "_")
    # Map common aliases.
    _ALIASES: dict[str, str] = {
        "homo_sapiens": "human",
        "h_sapiens": "human",
        "mus_musculus": "mouse",
        "m_musculus": "mouse",
        "s_cerevisiae": "yeast",
        "saccharomyces_cerevisiae": "yeast",
        "e_coli": "e_coli",
        "escherichia_coli": "e_coli",
        "cricetulus_griseus": "cho",
        "cho_k1": "cho",
        "c_elegans": "c_elegans",
        "caenorhabditis_elegans": "c_elegans",
        "d_melanogaster": "d_melanogaster",
        "drosophila_melanogaster": "d_melanogaster",
        "a_thaliana": "a_thaliana",
        "arabidopsis_thaliana": "a_thaliana",
        "b_subtilis": "b_subtilis",
        "bacillus_subtilis": "b_subtilis",
        "p_pastoris": "p_pastoris",
        "pichia_pastoris": "p_pastoris",
        "komagataella_phaffii": "p_pastoris",
    }
    org_key = _ALIASES.get(org_key, org_key)

    # Get tAI weights for this organism.
    if org_key in SUPPORTED_ORGANISMS_TAI:
        tai_weights = compute_tai_weights(org_key)
    else:
        logger.warning(
            "Organism '%s' not in tAI database; falling back to '%s'. "
            "Supported: %s",
            organism,
            _DEFAULT_ORGANISM,
            SUPPORTED_ORGANISMS_TAI,
        )
        tai_weights = compute_tai_weights(_DEFAULT_ORGANISM)
        org_key = _DEFAULT_ORGANISM

    # Get baseline elongation rate for this organism.
    base_rate = _ELONGATION_RATES.get(org_key, _ELONGATION_RATES[_DEFAULT_ORGANISM])

    # Temperature correction: rate(T) = rate(T_ref) * Q10^((T - T_ref)/10)
    temp_factor = _Q10 ** ((temperature - _REF_TEMP) / 10.0)
    effective_rate = base_rate * temp_factor

    # Compute dwell times.
    dwell_times: list[float] = []
    for i in range(n_codons):
        dna_codon = seq[i * 3 : i * 3 + 3]
        rna_codon = dna_codon.replace("T", "U")

        # tAI weight for this codon; default to 0.5 if not found.
        w = tai_weights.get(rna_codon, 0.5)
        if w <= 0:
            w = 0.01  # floor to avoid infinite dwell times

        # Dwell time = 1 / (effective_rate * tAI_weight).
        # Low tAI → long dwell (slow codon).
        dwell = 1.0 / (effective_rate * w)
        dwell_times.append(dwell)

    return dwell_times


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Optional Pinetree Integration
# ═══════════════════════════════════════════════════════════════════════════════


def simulate_pinetree(
    dna_seq: str,
    organism: str = "human",
    **kwargs: Any,
) -> dict[str, Any]:
    """Simulate ribosome traffic using the *pinetree* package.

    `pinetree <https://github.com/benjaminjack/pinetree>`_ is a
    stochastic simulator for co-translational processes that includes
    dynamic tRNA charging and ribosome recycling.  It is particularly
    useful for modelling time-varying tRNA availability, which the
    standard Gillespie TASEP assumes to be constant.

    When *pinetree* is not installed, this function returns a
    descriptive error dict instead of raising.

    Args:
        dna_seq: DNA coding sequence.
        organism: Organism key for tAI-based dwell times.  Default
            ``"human"``.
        **kwargs: Forwarded to :func:`simulate_gillespie` as fallback
            parameters if pinetree is unavailable (``initiation_rate``,
            ``elongation_rate``, ``ribosome_footprint``, ``max_time``,
            ``seed``).

    Returns:
        A dict with simulation results.  If pinetree is available, the
        dict contains pinetree-native output fields.  If not, the dict
        contains the key ``"error"`` plus the Gillespie TASEP fallback
        result.

    References:
        - Jack, B., Shi, M. & Bhatt, A. (2019) "Pinetree: a
          stochastic simulator of co-transcriptional translation."
          *Bioinformatics* 35(22):4739-4740.
    """
    # Compute dwell times from the DNA sequence for the fallback path.
    dwell_times = compute_dwell_times(dna_seq, organism=organism)

    try:
        import pinetree  # type: ignore[import-untyped]

        # Pinetree simulation setup.
        seq_upper = dna_seq.upper()
        n_codons = len(seq_upper) // 3
        initiation_rate = kwargs.get("initiation_rate", 0.1)
        elongation_rate = kwargs.get("elongation_rate", 10.0)
        max_time = kwargs.get("max_time", 1000.0)
        seed = kwargs.get("seed", 42)

        # Build the pinetree simulator.
        sim = pinetree.Simulator(seed=seed)

        # Register the transcript with per-codon elongation rates.
        # Pinetree uses a per-transcript model; we convert dwell times
        # to per-position rates.
        rates = [1.0 / max(dt, 0.001) for dt in dwell_times]
        max_r = max(rates) if rates else 1.0
        rates = [r / max_r * elongation_rate for r in rates]

        # Pinetree expects a gene with start/stop and per-codon rates.
        # We create a simple transcript specification.
        transcript = sim.add_gene(
            name="transcript",
            start=0,
            end=len(seq_upper),
            initiation_rate=initiation_rate,
        )

        # Set per-codon elongation rates.
        for i, rate in enumerate(rates):
            transcript.set_elongation_rate(pos=i * 3, rate=rate)

        # Run the simulation.
        sim.run(time=max_time)

        # Extract ribosome density from pinetree output.
        ribosome_counts = sim.ribosome_counts.get("transcript", [])
        time_elapsed = sim.time

        # Compute per-codon density from pinetree ribosome positions.
        n = n_codons
        codon_density = [0.0] * n
        if ribosome_counts:
            for pos in ribosome_counts:
                codon_idx = pos // 3
                if 0 <= codon_idx < n:
                    codon_density[codon_idx] += 1
            max_density = max(codon_density) if codon_density else 1.0
            if max_density > 0:
                codon_density = [d / max_density for d in codon_density]

        mean_density = sum(codon_density) / max(1, len(codon_density))
        stall_sites = [i for i, d in enumerate(codon_density) if d > 3 * mean_density and mean_density > 0]

        return {
            "codon_density": codon_density,
            "stall_sites": stall_sites,
            "time_elapsed": time_elapsed,
            "method": "pinetree",
            "ribosome_count": len(ribosome_counts) if ribosome_counts else 0,
        }

    except ImportError:
        logger.info(
            "pinetree not installed; falling back to Gillespie TASEP. "
            "Install with: pip install pinetree"
        )
    except Exception as exc:
        logger.warning(
            "pinetree simulation failed (%s); falling back to Gillespie TASEP",
            exc,
        )

    # ── Fallback: Gillespie TASEP ──────────────────────────────────────────
    result = simulate_gillespie(dwell_times=dwell_times, **kwargs)
    return {
        "codon_density": result.codon_density,
        "stall_sites": result.stall_sites,
        "collision_count": result.collision_count,
        "mean_ribosome_count": result.mean_ribosome_count,
        "time_elapsed": result.time_elapsed,
        "method": "gillespie_fallback",
        "error": "pinetree not available; used Gillespie TASEP fallback",
    }
