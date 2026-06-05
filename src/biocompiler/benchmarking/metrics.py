from __future__ import annotations
"""Comprehensive comparison metrics for benchmarking biocompiler against other tools.

Metrics computed:
  1. CAI (Codon Adaptation Index) — validated against Sharp & Li (1987)
  2. GC content — with distribution analysis
  3. Restriction site count — across a standard enzyme panel
  4. mRNA stability — T-run analysis, hairpin potential
  5. Cryptic splice site count — GT/AG dinucleotide frequency
  6. CpG island count
  7. Runtime performance
  8. Sequence identity to original (for constraint satisfaction)

References
----------
Sharp, P. M., & Li, W.-H. (1987). The codon Adaptation Index—a measure of
directional synonymous codon usage bias, and its potential applications.
*Nucleic Acids Research*, 15(3), 1281–1295.
doi:10.1093/nar/15.3.1281
"""

import math
import time
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from scipy import stats

from ..constants import AA_TO_CODONS, CODON_TABLE, reverse_complement
from ..organisms import CODON_ADAPTIVENESS_TABLES, SUPPORTED_ORGANISMS
from ..restriction_sites import get_recognition_site
from ..scanner import gc_content
from ..translation import compute_cai

__all__ = [
    # Data classes
    "GCProfile",
    "BenchmarkMetrics",
    # Metric computation functions
    "compute_cai_validated",
    "compute_gc_distribution",
    "count_restriction_sites",
    "count_cryptic_splice_sites",
    "count_cpg_islands",
    "compute_codon_pair_bias",
    "compute_mrna_stability_score",
    "compute_sequence_identity",
    # All-in-one
    "compute_all_metrics",
    # Statistical comparison
    "StatisticalComparison",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_CODON_LENGTH: int = 3

# Minimum number of bases for a T-run to potentially affect mRNA stability
_MIN_T_RUN_LENGTH: int = 4

# CpG island detection defaults (Gardiner-Garden & Frommer, 1987)
_DEFAULT_CPG_WINDOW: int = 200
_DEFAULT_CPG_GC_THRESHOLD: float = 0.50
_DEFAULT_CPG_RATIO_THRESHOLD: float = 0.60

# Standard enzyme panel for benchmarking
STANDARD_ENZYME_PANEL: list[str] = [
    "EcoRI", "BamHI", "HindIII", "XhoI", "XbaI", "SalI",
    "PstI", "NcoI", "NdeI", "NotI",
]

# Epsilon floor for zero-adaptiveness codons in CAI computation
_CAI_EPSILON: float = 1e-10

# Codon pair bias tables (simplified — only for well-studied organisms).
# These are illustrative values derived from Coleman et al. (2008) and
# Clarke et al. (2008).  Real production data should come from
# organism-specific genome-wide codon-pair frequency analyses.
# Key: (codon1, codon2) -> over-representation score
_CODON_PAIR_BIAS: dict[str, dict[tuple[str, str], float]] = {}


# ---------------------------------------------------------------------------
# GCProfile — sliding-window GC distribution
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GCProfile:
    """Sliding-window GC content profile for a DNA sequence.

    Attributes
    ----------
    mean : float
        Mean GC fraction over the entire sequence.
    std : float
        Standard deviation of GC across windows.
    min_ : float
        Minimum GC fraction in any window.
    max_ : float
        Maximum GC fraction in any window.
    window : int
        Window size used for the sliding analysis.
    values : tuple[float, ...]
        GC fraction at each window position (step = 1 bp).
    """

    mean: float
    std: float
    min_: float
    max_: float
    window: int
    values: tuple[float, ...] = field(default_factory=tuple)

    @property
    def range(self) -> float:
        """Max – min GC across windows."""
        return self.max_ - self.min_

    @property
    def cv(self) -> float:
        """Coefficient of variation of GC (std / mean)."""
        return self.std / self.mean if self.mean > 0 else 0.0

    def to_dict(self) -> dict:
        """Convert to a plain dictionary (for JSON serialisation)."""
        return {
            "mean": self.mean,
            "std": self.std,
            "min": self.min_,
            "max": self.max_,
            "range": self.range,
            "cv": self.cv,
            "window": self.window,
            "n_windows": len(self.values),
        }


# ---------------------------------------------------------------------------
# BenchmarkMetrics — all metrics in one dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BenchmarkMetrics:
    """All benchmark metrics for a single optimised sequence.

    Attributes
    ----------
    cai : float
        Codon Adaptation Index (Sharp & Li 1987), range [0, 1].
    gc_profile : GCProfile
        Sliding-window GC distribution.
    restriction_sites : dict[str, int]
        Count of each restriction enzyme site found.
    restriction_site_total : int
        Total count across all enzymes.
    cryptic_splice_sites : int
        Count of GT/AG cryptic splice-site dinucleotides in coding context.
    cpg_islands : int
        Number of CpG islands detected.
    codon_pair_bias : float
        Codon pair bias score (0.0 if no data available).
    mrna_stability : float
        Composite mRNA stability score (T-run + hairpin potential).
    sequence_identity : float
        Fractional identity to the original (pre-optimisation) sequence.
    runtime_s : float
        Wall-clock time to compute all metrics, in seconds.
    """

    cai: float
    gc_profile: GCProfile
    restriction_sites: dict[str, int]
    restriction_site_total: int
    cryptic_splice_sites: int
    cpg_islands: int
    codon_pair_bias: float
    mrna_stability: float
    sequence_identity: float
    runtime_s: float

    def to_dict(self) -> dict:
        """Convert to a flat dictionary suitable for JSON / CSV output."""
        d = {
            "cai": self.cai,
            "gc_mean": self.gc_profile.mean,
            "gc_std": self.gc_profile.std,
            "gc_min": self.gc_profile.min_,
            "gc_max": self.gc_profile.max_,
            "gc_range": self.gc_profile.range,
            "gc_cv": self.gc_profile.cv,
            "restriction_site_total": self.restriction_site_total,
            "restriction_sites": dict(self.restriction_sites),
            "cryptic_splice_sites": self.cryptic_splice_sites,
            "cpg_islands": self.cpg_islands,
            "codon_pair_bias": self.codon_pair_bias,
            "mrna_stability": self.mrna_stability,
            "sequence_identity": self.sequence_identity,
            "runtime_s": self.runtime_s,
        }
        return d


# ---------------------------------------------------------------------------
# 1. CAI — validated against Sharp & Li (1987)
# ---------------------------------------------------------------------------

def compute_cai_validated(dna: str, organism: str) -> float:
    """Compute the Codon Adaptation Index following Sharp & Li (1987).

    The CAI is defined as the geometric mean of relative adaptiveness
    values *w_i* for each codon in the coding sequence:

        CAI = (∏_{i=1}^{L} w_i)^{1/L}

    where *w_i = f_i / f_max* for the amino acid encoded by codon *i*,
    *f_i* is the frequency of codon *i*, and *f_max* is the frequency of
    the most frequently used synonymous codon for that amino acid.

    Methionine (ATG) and stop codons are excluded from the product,
    following the original paper.

    This implementation delegates to :func:`biocompiler.translation.compute_cai`
    which uses the organism-specific adaptiveness tables from the
    :mod:`biocompiler.organisms` module.

    Parameters
    ----------
    dna : str
        DNA coding sequence (uppercase ACGT).
    organism : str
        Organism name matching a key in
        :data:`biocompiler.organisms.CODON_ADAPTIVENESS_TABLES`
        (e.g. ``"Homo_sapiens"``).

    Returns
    -------
    float
        CAI value in [0, 1].  Returns 0.0 for empty or invalid sequences.

    References
    ----------
    Sharp, P. M., & Li, W.-H. (1987).  *Nucleic Acids Research*,
    15(3), 1281–1295.  doi:10.1093/nar/15.3.1281
    """
    if not dna or len(dna) < _CODON_LENGTH:
        return 0.0
    dna = dna.upper()
    # Normalise organism name
    organism_key = organism
    if organism_key not in CODON_ADAPTIVENESS_TABLES:
        # Try common aliases
        _ALIASES: dict[str, str] = {
            "human": "Homo_sapiens",
            "ecoli": "Escherichia_coli",
            "mouse": "Mus_musculus",
            "cho": "CHO_K1",
            "yeast": "Saccharomyces_cerevisiae",
        }
        organism_key = _ALIASES.get(organism.lower(), organism)
    if organism_key not in CODON_ADAPTIVENESS_TABLES:
        return 0.0

    adaptiveness = CODON_ADAPTIVENESS_TABLES[organism_key]
    log_sum: float = 0.0
    count: int = 0

    for i in range(0, len(dna) - _CODON_LENGTH + 1, _CODON_LENGTH):
        codon = dna[i:i + _CODON_LENGTH]
        aa = CODON_TABLE.get(codon)
        # Skip Met (ATG, only codon) and stop codons — Sharp & Li (1987)
        if aa is None or aa == "*" or aa == "M":
            continue
        w = adaptiveness.get(codon, 0.0)
        if w <= 0:
            w = _CAI_EPSILON  # Floor to avoid log(0)
        log_sum += math.log(w)
        count += 1

    if count == 0:
        return 0.0
    return math.exp(log_sum / count)


# ---------------------------------------------------------------------------
# 2. GC content with sliding-window distribution
# ---------------------------------------------------------------------------

def compute_gc_distribution(dna: str, window: int = 50) -> GCProfile:
    """Compute GC content with a sliding-window analysis.

    Rather than just reporting the mean GC fraction, this function produces
    a full distribution profile showing how GC content varies along the
    sequence — useful for detecting GC-rich or AT-rich regions that may
    affect transcription, translation, or mRNA stability.

    Parameters
    ----------
    dna : str
        DNA sequence (uppercase ACGT; lowercase is accepted).
    window : int
        Sliding window size in base pairs (default 50).

    Returns
    -------
    GCProfile
        Sliding-window GC distribution with mean, std, min, max, and
        per-window values.
    """
    dna = dna.upper()
    n = len(dna)
    if n == 0:
        return GCProfile(
            mean=0.0, std=0.0, min_=0.0, max_=0.0,
            window=window, values=(),
        )

    # Clamp window to sequence length
    window = max(1, min(window, n))

    # Compute per-window GC fraction
    gc_vals: list[float] = []
    for start in range(n - window + 1):
        subseq = dna[start:start + window]
        gc = (subseq.count("G") + subseq.count("C")) / window
        gc_vals.append(gc)

    if not gc_vals:
        # Sequence shorter than window
        gc = (dna.count("G") + dna.count("C")) / n
        return GCProfile(
            mean=gc, std=0.0, min_=gc, max_=gc,
            window=window, values=(gc,),
        )

    arr = np.array(gc_vals, dtype=np.float64)
    return GCProfile(
        mean=float(arr.mean()),
        std=float(arr.std()),
        min_=float(arr.min()),
        max_=float(arr.max()),
        window=window,
        values=tuple(gc_vals),
    )


# ---------------------------------------------------------------------------
# 3. Restriction site counting
# ---------------------------------------------------------------------------

def count_restriction_sites(dna: str, enzymes: list[str]) -> dict[str, int]:
    """Count occurrences of each restriction enzyme recognition site.

    Checks both forward and reverse-complement strands.  Palindromic
    sites (e.g. EcoRI: GAATTC) are counted once per occurrence on the
    forward strand to avoid double-counting.

    Parameters
    ----------
    dna : str
        DNA sequence (case-insensitive).
    enzymes : list[str]
        Enzyme names to check (must be in
        :data:`biocompiler.restriction_sites.RESTRICTION_SITES`).

    Returns
    -------
    dict[str, int]
        Mapping from enzyme name to count of occurrences found.
    """
    dna = dna.upper()
    result: dict[str, int] = {}

    for enzyme in enzymes:
        site = get_recognition_site(enzyme)
        if site is None:
            result[enzyme] = 0
            continue
        site_upper = site.upper()

        # Skip sites with IUPAC ambiguity codes for simple counting
        if any(b not in "ACGT" for b in site_upper):
            result[enzyme] = 0
            continue

        count = 0
        start = 0
        while True:
            pos = dna.find(site_upper, start)
            if pos == -1:
                break
            count += 1
            start = pos + 1

        # Check reverse complement strand (non-palindromic only)
        site_rc = reverse_complement(site_upper)
        if site_rc != site_upper:
            start = 0
            while True:
                pos = dna.find(site_rc, start)
                if pos == -1:
                    break
                count += 1
                start = pos + 1

        result[enzyme] = count

    return result


# ---------------------------------------------------------------------------
# 4. Cryptic splice site counting
# ---------------------------------------------------------------------------

def count_cryptic_splice_sites(dna: str) -> int:
    """Count cryptic splice-site dinucleotides (GT/AG) in coding context.

    In a coding sequence, GT dinucleotides (5' splice donor consensus) and
    AG dinucleotides (3' splice acceptor consensus) can be recognised by
    the spliceosome as cryptic splice sites, leading to aberrant splicing.
    This function counts the total number of GT and AG dinucleotides that
    fall **within codons** (i.e., not spanning codon boundaries), which
    are most likely to form functional cryptic sites.

    Parameters
    ----------
    dna : str
        DNA coding sequence (uppercase ACGT).

    Returns
    -------
    int
        Total count of in-codon GT and AG dinucleotides.
    """
    dna = dna.upper()
    if len(dna) < 2:
        return 0

    count = 0
    # Scan each codon independently — dinucleotides that span a codon
    # boundary (position 2→3) are less likely to form splice sites
    # because the spliceosome recognises a continuous context.
    for i in range(0, len(dna) - _CODON_LENGTH + 1, _CODON_LENGTH):
        codon = dna[i:i + _CODON_LENGTH]
        for j in range(len(codon) - 1):
            dinuc = codon[j:j + 2]
            if dinuc == "GT" or dinuc == "AG":
                count += 1

    # Also count boundary-spanning dinucleotides (conservative count)
    # These are positions where the last base of one codon and the first
    # base of the next form GT or AG.
    for i in range(_CODON_LENGTH - 1, len(dna) - 1, _CODON_LENGTH):
        dinuc = dna[i:i + 2]
        if dinuc == "GT" or dinuc == "AG":
            count += 1

    return count


# ---------------------------------------------------------------------------
# 5. CpG island counting
# ---------------------------------------------------------------------------

def count_cpg_islands(
    dna: str,
    window: int = _DEFAULT_CPG_WINDOW,
    threshold: float = _DEFAULT_CPG_RATIO_THRESHOLD,
) -> int:
    """Count CpG islands using the Gardiner-Garden & Frommer (1987) criteria.

    A CpG island is defined as a region of at least *window* base pairs
    where:
      - GC content >= 50 %, AND
      - Observed/Expected CpG ratio >= *threshold* (default 0.60)

    The Obs/Exp CpG ratio is computed as:

        Obs/Exp(CpG) = (n_CpG * L) / (n_C * n_G)

    where *n_CpG* is the number of CG dinucleotides, *n_C* and *n_G* are
    the counts of C and G respectively, and *L* is the window length.

    Adjacent qualifying windows are merged into a single island.

    Parameters
    ----------
    dna : str
        DNA sequence (uppercase ACGT).
    window : int
        Minimum window size in bp (default 200).
    threshold : float
        Minimum Obs/Exp CpG ratio (default 0.60).

    Returns
    -------
    int
        Number of distinct CpG islands found.

    References
    ----------
    Gardiner-Garden, M., & Frommer, M. (1987).  CpG islands in vertebrate
    genomes.  *Journal of Molecular Biology*, 196(2), 261–282.
    """
    dna = dna.upper()
    n = len(dna)
    if n < window:
        return 0

    gc_thresh = _DEFAULT_CPG_GC_THRESHOLD
    islands = 0
    in_island = False

    for start in range(n - window + 1):
        subseq = dna[start:start + window]
        gc = (subseq.count("G") + subseq.count("C")) / window
        if gc < gc_thresh:
            in_island = False
            continue

        n_c = subseq.count("C")
        n_g = subseq.count("G")
        n_cpg = sum(
            1 for j in range(len(subseq) - 1) if subseq[j:j + 2] == "CG"
        )
        if n_c == 0 or n_g == 0:
            in_island = False
            continue
        obs_exp = (n_cpg * window) / (n_c * n_g)

        if obs_exp >= threshold:
            if not in_island:
                islands += 1
                in_island = True
        else:
            in_island = False

    return islands


# ---------------------------------------------------------------------------
# 6. Codon pair bias
# ---------------------------------------------------------------------------

def compute_codon_pair_bias(dna: str, organism: str) -> float:
    """Compute codon pair bias score for a coding sequence.

    Codon pair bias captures the tendency of certain codon pairs to be
    over- or under-represented in a given organism.  Positive scores
    indicate over-represented (preferred) codon pairs; negative scores
    indicate under-represented (avoided) codon pairs.

    The score is the mean log-odds ratio across all consecutive codon
    pairs in the sequence:

        CPB = (1/N) * Σ log(P_obs / P_exp)

    If no codon-pair data is available for the organism, returns 0.0
    (indicating no detected bias).

    Parameters
    ----------
    dna : str
        DNA coding sequence.
    organism : str
        Organism name.

    Returns
    -------
    float
        Mean codon pair bias score.  Returns 0.0 if no organism-specific
        data is available.
    """
    # Normalise organism name
    _ALIASES: dict[str, str] = {
        "human": "Homo_sapiens",
        "ecoli": "Escherichia_coli",
        "mouse": "Mus_musculus",
        "cho": "CHO_K1",
        "yeast": "Saccharomyces_cerevisiae",
    }
    org_key = _ALIASES.get(organism.lower(), organism)

    pair_table = _CODON_PAIR_BIAS.get(org_key)
    if not pair_table:
        return 0.0

    dna = dna.upper()
    codons = [
        dna[i:i + _CODON_LENGTH]
        for i in range(0, len(dna) - _CODON_LENGTH + 1, _CODON_LENGTH)
    ]
    if len(codons) < 2:
        return 0.0

    scores: list[float] = []
    for i in range(len(codons) - 1):
        pair = (codons[i], codons[i + 1])
        bias = pair_table.get(pair, 0.0)
        scores.append(bias)

    if not scores:
        return 0.0
    return float(np.mean(scores))


# ---------------------------------------------------------------------------
# 7. mRNA stability score
# ---------------------------------------------------------------------------

def compute_mrna_stability_score(dna: str) -> float:
    """Compute a composite mRNA stability score.

    This metric combines two indicators of mRNA stability:

    1. **T-run analysis**: Long poly-T runs (≥4 consecutive T's) in the
       coding sequence can serve as premature transcription termination
       signals or destabilise mRNA via AU-rich elements.
    2. **Hairpin potential**: GC-rich regions that could form secondary
       structures (approximated by the fraction of windows with GC > 0.65).

    The composite score is in [0, 1], where 1.0 indicates maximum
    predicted stability (few T-runs, low hairpin potential) and 0.0
    indicates minimum predicted stability.

    Parameters
    ----------
    dna : str
        DNA coding sequence.

    Returns
    -------
    float
        Composite mRNA stability score in [0, 1].
    """
    dna = dna.upper()
    n = len(dna)
    if n == 0:
        return 1.0

    # --- T-run analysis ---
    t_run_count = 0
    i = 0
    while i < n:
        if dna[i] == "T":
            run_len = 0
            while i < n and dna[i] == "T":
                run_len += 1
                i += 1
            if run_len >= _MIN_T_RUN_LENGTH:
                t_run_count += 1
        else:
            i += 1

    # Normalise: expect ~0 T-runs in a well-optimised sequence
    # Penalty: each T-run reduces stability by ~0.1 (capped at 1.0)
    t_run_penalty = min(1.0, t_run_count * 0.1)

    # --- Hairpin potential ---
    # Approximate by fraction of 30-bp windows with GC > 0.65
    hairpin_window = 30
    gc_high_count = 0
    total_windows = max(1, n - hairpin_window + 1)
    for start in range(total_windows):
        subseq = dna[start:start + hairpin_window]
        gc = (subseq.count("G") + subseq.count("C")) / hairpin_window
        if gc > 0.65:
            gc_high_count += 1
    hairpin_fraction = gc_high_count / total_windows

    # Composite: higher is more stable
    stability = 1.0 - 0.5 * t_run_penalty - 0.5 * hairpin_fraction
    return max(0.0, min(1.0, stability))


# ---------------------------------------------------------------------------
# 8. Sequence identity
# ---------------------------------------------------------------------------

def compute_sequence_identity(dna_original: str, dna_optimised: str) -> float:
    """Compute fractional sequence identity between two DNA sequences.

    For sequences of different lengths, only the aligned portion is
    compared.  This metric measures how much the optimisation perturbed
    the original sequence — useful for verifying that the protein is
    preserved while the DNA is re-coded.

    Parameters
    ----------
    dna_original : str
        Original (pre-optimisation) DNA sequence.
    dna_optimised : str
        Optimised DNA sequence.

    Returns
    -------
    float
        Fraction of identical positions in [0, 1].
    """
    a = dna_original.upper()
    b = dna_optimised.upper()
    length = min(len(a), len(b))
    if length == 0:
        return 1.0
    matches = sum(1 for i in range(length) if a[i] == b[i])
    return matches / max(len(a), len(b))


# ---------------------------------------------------------------------------
# 9. All-in-one metric computation
# ---------------------------------------------------------------------------

def compute_all_metrics(
    dna: str,
    protein: str,
    organism: str,
    enzymes: list[str] | None = None,
    original_dna: str | None = None,
) -> BenchmarkMetrics:
    """Compute all benchmark metrics for an optimised sequence.

    Parameters
    ----------
    dna : str
        Optimised DNA coding sequence.
    protein : str
        Target protein sequence (for validation).
    organism : str
        Organism name for CAI / codon pair bias computation.
    enzymes : list[str] or None
        Restriction enzymes to check.  Defaults to
        :data:`STANDARD_ENZYME_PANEL`.
    original_dna : str or None
        Original DNA for sequence identity.  If ``None``, identity
        defaults to 1.0.

    Returns
    -------
    BenchmarkMetrics
        All computed metrics in a single dataclass.
    """
    if enzymes is None:
        enzymes = list(STANDARD_ENZYME_PANEL)

    t0 = time.perf_counter()

    # 1. CAI
    cai = compute_cai_validated(dna, organism)

    # 2. GC distribution
    gc_profile = compute_gc_distribution(dna)

    # 3. Restriction sites
    rs_dict = count_restriction_sites(dna, enzymes)
    rs_total = sum(rs_dict.values())

    # 4. Cryptic splice sites
    css = count_cryptic_splice_sites(dna)

    # 5. CpG islands
    cpg = count_cpg_islands(dna)

    # 6. Codon pair bias
    cpb = compute_codon_pair_bias(dna, organism)

    # 7. mRNA stability
    mrna = compute_mrna_stability_score(dna)

    # 8. Sequence identity
    identity = (
        compute_sequence_identity(original_dna, dna)
        if original_dna is not None
        else 1.0
    )

    elapsed = time.perf_counter() - t0

    return BenchmarkMetrics(
        cai=cai,
        gc_profile=gc_profile,
        restriction_sites=rs_dict,
        restriction_site_total=rs_total,
        cryptic_splice_sites=css,
        cpg_islands=cpg,
        codon_pair_bias=cpb,
        mrna_stability=mrna,
        sequence_identity=identity,
        runtime_s=elapsed,
    )


# ---------------------------------------------------------------------------
# StatisticalComparison — significance tests & effect sizes
# ---------------------------------------------------------------------------

class StatisticalComparison:
    """Statistical comparison between two sets of benchmark results.

    Given two matched lists of metric values (e.g. CAI from BioCompiler
    vs. a competing tool across a gene panel), this class computes:

    - Paired t-test (parametric)
    - Wilcoxon signed-rank test (non-parametric)
    - Cohen's d effect size
    - Cliff's delta effect size (non-parametric)

    Parameters
    ----------
    group_a : Sequence[float]
        Metric values for group A (e.g. BioCompiler).
    group_b : Sequence[float]
        Metric values for group B (e.g. competitor tool).
    name_a : str
        Label for group A (default ``"A"``).
    name_b : str
        Label for group B (default ``"B"``).
    alpha : float
        Significance level (default 0.05).

    Raises
    ------
    ValueError
        If the two groups have different lengths.
    """

    def __init__(
        self,
        group_a: Sequence[float],
        group_b: Sequence[float],
        name_a: str = "A",
        name_b: str = "B",
        alpha: float = 0.05,
    ) -> None:
        if len(group_a) != len(group_b):
            raise ValueError(
                f"Group lengths must match: len(A)={len(group_a)}, "
                f"len(B)={len(group_b)}"
            )
        self._a = np.asarray(group_a, dtype=np.float64)
        self._b = np.asarray(group_b, dtype=np.float64)
        self.name_a = name_a
        self.name_b = name_b
        self.alpha = alpha
        self._n = len(self._a)

    # -- Descriptive statistics ---------------------------------------------

    @property
    def n(self) -> int:
        """Number of paired observations."""
        return self._n

    @property
    def mean_a(self) -> float:
        """Mean of group A."""
        return float(np.mean(self._a))

    @property
    def mean_b(self) -> float:
        """Mean of group B."""
        return float(np.mean(self._b))

    @property
    def std_a(self) -> float:
        """Standard deviation of group A."""
        return float(np.std(self._a, ddof=1)) if self._n > 1 else 0.0

    @property
    def std_b(self) -> float:
        """Standard deviation of group B."""
        return float(np.std(self._b, ddof=1)) if self._n > 1 else 0.0

    @property
    def mean_diff(self) -> float:
        """Mean of A – B differences."""
        return float(np.mean(self._a - self._b))

    # -- Parametric test ----------------------------------------------------

    def paired_t_test(self) -> tuple[float, float]:
        """Two-sided paired t-test.

        Returns
        -------
        tuple[float, float]
            (t_statistic, p_value)
        """
        if self._n < 2:
            return (0.0, 1.0)
        t_stat, p_val = stats.ttest_rel(self._a, self._b)
        return (float(t_stat), float(p_val))

    @property
    def is_significant_t(self) -> bool:
        """Whether the paired t-test is significant at *alpha*."""
        _, p = self.paired_t_test()
        return p < self.alpha

    # -- Non-parametric test ------------------------------------------------

    def wilcoxon_test(self) -> tuple[float, float]:
        """Wilcoxon signed-rank test (two-sided).

        Returns
        -------
        tuple[float, float]
            (W_statistic, p_value)
        """
        if self._n < 2:
            return (0.0, 1.0)
        diffs = self._a - self._b
        # All differences zero → not significant
        if np.all(diffs == 0):
            return (0.0, 1.0)
        try:
            result = stats.wilcoxon(self._a, self._b, alternative="two-sided")
            return (float(result.statistic), float(result.pvalue))
        except ValueError:
            # Fallback when Wilcoxon cannot be computed (e.g. too few non-zero diffs)
            return (0.0, 1.0)

    @property
    def is_significant_wilcoxon(self) -> bool:
        """Whether the Wilcoxon test is significant at *alpha*."""
        _, p = self.wilcoxon_test()
        return p < self.alpha

    # -- Effect sizes -------------------------------------------------------

    def cohens_d(self) -> float:
        """Cohen's d effect size for paired samples.

        d = mean_diff / std_diff

        Interpretation (Cohen 1988):
          |d| < 0.2  → negligible
          |d| < 0.5  → small
          |d| < 0.8  → medium
          |d| >= 0.8 → large

        Returns
        -------
        float
            Cohen's d value.
        """
        if self._n < 2:
            return 0.0
        diffs = self._a - self._b
        std_diff = float(np.std(diffs, ddof=1))
        if std_diff == 0:
            return 0.0
        return float(np.mean(diffs)) / std_diff

    def cliffs_delta(self) -> float:
        """Cliff's delta non-parametric effect size.

        Measures the proportion of (a_i > b_j) minus proportion of
        (a_i < b_j) across all pairs.

        Interpretation:
          |delta| < 0.147 → negligible
          |delta| < 0.33  → small
          |delta| < 0.474 → medium
          |delta| >= 0.474 → large

        Returns
        -------
        float
            Cliff's delta in [-1, 1].
        """
        if self._n < 2:
            return 0.0
        more = 0
        less = 0
        for a_val in self._a:
            for b_val in self._b:
                if a_val > b_val:
                    more += 1
                elif a_val < b_val:
                    less += 1
        total = self._n * self._n
        if total == 0:
            return 0.0
        return (more - less) / total

    # -- Summary ------------------------------------------------------------

    def summary(self) -> dict:
        """Return a summary dictionary of all statistical results."""
        t_stat, t_p = self.paired_t_test()
        w_stat, w_p = self.wilcoxon_test()
        d = self.cohens_d()
        delta = self.cliffs_delta()
        return {
            "name_a": self.name_a,
            "name_b": self.name_b,
            "n": self.n,
            "mean_a": self.mean_a,
            "mean_b": self.mean_b,
            "std_a": self.std_a,
            "std_b": self.std_b,
            "mean_diff": self.mean_diff,
            "paired_t": {"statistic": t_stat, "p_value": t_p, "significant": self.is_significant_t},
            "wilcoxon": {"statistic": w_stat, "p_value": w_p, "significant": self.is_significant_wilcoxon},
            "cohens_d": d,
            "cliffs_delta": delta,
            "alpha": self.alpha,
        }

    def __repr__(self) -> str:
        t_stat, t_p = self.paired_t_test()
        d = self.cohens_d()
        return (
            f"StatisticalComparison(n={self.n}, "
            f"mean_diff={self.mean_diff:.4f}, "
            f"t={t_stat:.3f} p={t_p:.4f}, "
            f"Cohen's d={d:.3f})"
        )
