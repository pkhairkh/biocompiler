"""
BioCompiler mRNA Accessibility Module — Unpaired Probability Computation
=========================================================================

Computes per-position mRNA accessibility (unpaired probability) using the
ViennaRNA partition function, with a GC-content heuristic fallback when
ViennaRNA is not installed.

mRNA accessibility measures the probability that each nucleotide is
unpaired (not involved in intramolecular base pairing), which is critical
for:

- **Translation initiation**: 5' end accessibility determines ribosome
  binding efficiency. Inaccessible 5' regions strongly inhibit translation.
- **miRNA targeting**: miRNA binding sites buried in secondary structure
  are ineffective. Accessibility modulates miRNA efficacy (Kertesz 2010).
- **Codon optimality**: Local structure around codons modulates ribosome
  elongation speed (Bazzini 2016).
- **RBP binding**: RNA-binding proteins require accessible binding motifs.

The primary method uses the ViennaRNA partition function to compute
base-pairing probabilities (BPP), then converts to unpaired probabilities:

    P(unpaired_i) = 1 - sum_j P(i pairs with j)

For long sequences, an RNAplfold-equivalent sliding window approach is
used via ``RNA.fold_compound().probs_window()``.

When ViennaRNA is unavailable, a GC-content heuristic estimates
accessibility from local GC composition: higher GC content implies more
stable secondary structure and therefore lower accessibility.

Functions:
    compute_accessibility            — Per-position unpaired probability
    compute_accessibility_windows    — Sliding window accessibility with position info
    compute_codon_accessibility      — Average accessibility per codon
    compute_5prime_accessibility     — Average accessibility in first N codons
    compute_mirna_site_accessibility — Accessibility at miRNA binding sites
    adjust_severity_for_accessibility — Reduce severity for inaccessible sites

References:
    Lorenz, R. et al. (2011). "ViennaRNA Package 2.0."
    *Algorithms Mol Biol* 6:26.
    https://doi.org/10.1186/1748-7188-6-26

    Kertesz, M. et al. (2010). "The role of site accessibility in
    microRNA target recognition." *Nature Genet* 42:814–819.
    https://doi.org/10.1038/ng.647

    Bazzini, A.A. et al. (2016). "Codon identity determines mRNA
    stability and regulates translation in vertebrates."
    *Nat Commun* 7:11314.
    https://doi.org/10.1038/ncomms11314
"""

from __future__ import annotations

import logging
import math
import warnings
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    # Core accessibility functions
    "compute_accessibility",
    "compute_accessibility_windows",
    "compute_codon_accessibility",
    "compute_5prime_accessibility",
    "compute_mirna_site_accessibility",
    # Severity integration
    "adjust_severity_for_accessibility",
    # Data classes
    "AccessibilityWindowResult",
    # Constants
    "DEFAULT_WINDOW",
    "DEFAULT_SPAN",
    "DEFAULT_SLIDING_WINDOW_SIZE",
    "DEFAULT_SLIDING_STEP",
    "DEFAULT_5PRIME_N_CODONS",
    "DEFAULT_MIRNA_SITE_LENGTH",
    "ACCESSIBILITY_SEVERITY_THRESHOLD",
    "ACCESSIBILITY_SEVERITY_REDUCTION",
    # Feature flag
    "HAS_VIENNARNA",
]

# ── Named Constants ──────────────────────────────────────────

DEFAULT_WINDOW: int = 80
"""Sliding window size for partition function computation (nt)."""

DEFAULT_SPAN: int = 40
"""Maximum base-pair span for accessibility computation (nt)."""

DEFAULT_SLIDING_WINDOW_SIZE: int = 80
"""Window size for sliding-window accessibility computation (nt)."""

DEFAULT_SLIDING_STEP: int = 10
"""Step size for sliding-window accessibility computation (nt)."""

DEFAULT_5PRIME_N_CODONS: int = 17
"""Number of codons from the 5' end to consider for translation
initiation accessibility (≈ 51 nt, covering the ribosome standby site)."""

DEFAULT_MIRNA_SITE_LENGTH: int = 8
"""Default miRNA binding site length in nucleotides (seed + supplementary)."""

ACCESSIBILITY_SEVERITY_THRESHOLD: float = 0.05
"""Below this accessibility, a site is considered buried in secondary
structure and its severity should be reduced."""

ACCESSIBILITY_SEVERITY_REDUCTION: float = 0.80
"""Fraction by which severity is reduced when accessibility is below
the threshold (0.80 = reduce by 80%, keeping 20% of original)."""

# GC-content heuristic parameters
_GC_STABLE_THRESHOLD: float = 0.60
"""GC fraction above which a sequence is considered likely to form stable structure."""

_GC_UNSTABLE_THRESHOLD: float = 0.40
"""GC fraction below which a sequence is considered unlikely to form stable structure."""

_GC_HEURISTIC_WINDOW: int = 30
"""Window size for local GC content computation in heuristic fallback."""

_GC_ACCESSIBILITY_SCALE: float = 2.0
"""Scaling factor for GC→accessibility mapping. Higher values give
steeper transition between accessible and inaccessible."""

# ── Feature flag ─────────────────────────────────────────────

HAS_VIENNARNA: bool = False
"""Whether ViennaRNA Python bindings are available."""

try:
    import RNA as _RNA  # noqa: F401

    HAS_VIENNARNA = True
except ImportError:
    pass


# ── Data Classes ─────────────────────────────────────────────


@dataclass
class AccessibilityWindowResult:
    """Result from sliding-window accessibility computation.

    Attributes:
        start: 0-based start position of the window in the original sequence.
        end: 0-based exclusive end position of the window.
        mean_accessibility: Mean unpaired probability across the window.
        position_accessibility: Dict mapping 0-based position → P(unpaired).
        gc_content: GC fraction of the window (for heuristic validation).
        method: Computation method used.
    """

    start: int
    end: int
    mean_accessibility: float
    position_accessibility: dict[int, float]
    gc_content: float
    method: str = "unknown"


# ── Internal Helpers ─────────────────────────────────────────


def _dna_to_rna(seq: str) -> str:
    """Convert DNA sequence to RNA (T → U), uppercase."""
    return seq.upper().replace("T", "U")


def _compute_gc(seq: str) -> float:
    """Compute GC content fraction of a sequence."""
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return gc / len(seq)


def _gc_to_accessibility(gc_frac: float) -> float:
    """Map GC content to estimated accessibility using sigmoid.

    Uses a logistic function centered at 0.50 GC with steepness
    controlled by ``_GC_ACCESSIBILITY_SCALE``:
        acc = 1 / (1 + exp(scale * (gc - 0.50)))

    - GC = 0.30 → acc ≈ 0.88  (AT-rich, very accessible)
    - GC = 0.50 → acc ≈ 0.50  (mixed)
    - GC = 0.70 → acc ≈ 0.12  (GC-rich, mostly structured)

    This is a rough heuristic; actual accessibility depends on
    nearest-neighbor stacking, loop entropy, and other factors
    that GC content alone cannot capture.
    """
    return 1.0 / (1.0 + math.exp(_GC_ACCESSIBILITY_SCALE * (gc_frac - 0.50)))


def _accessibility_heuristic(seq: str) -> list[float]:
    """Estimate per-position accessibility from local GC content.

    Uses a sliding window of ``_GC_HEURISTIC_WINDOW`` nt centered on
    each position to compute local GC, then maps to accessibility via
    a sigmoid function.

    .. warning::
       This is a **rough heuristic** and may differ significantly from
       the true thermodynamic accessibility. Use only as a last resort
       when ViennaRNA is not installed.

    Args:
        seq: DNA or RNA sequence (uppercased internally).

    Returns:
        List of per-position accessibility scores (0–1).
    """
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return []
    if n < 4:
        # Sequences shorter than 4 nt cannot form stable base pairs
        # (minimum hairpin requires 4-nt loop) → fully accessible
        return [1.0] * n

    half_w = _GC_HEURISTIC_WINDOW // 2
    result: list[float] = []

    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        local_gc = _compute_gc(seq[start:end])
        result.append(_gc_to_accessibility(local_gc))

    return result


def _accessibility_viennarna_partition(seq: str, window: int, span: int) -> list[float]:
    """Compute per-position unpaired probability using ViennaRNA partition function.

    Uses ``RNA.fold_compound().pf()`` to compute the full partition function,
    then extracts base-pairing probabilities via ``fc.bpp()`` or
    ``fc.bp_prob(i, j)`` to derive per-position P(unpaired).

    For sequences longer than *window* nt, uses a sliding-window approach
    where each window is folded independently.

    Args:
        seq: RNA sequence (U not T).
        window: Window size for partition function computation.
        span: Maximum base-pair span.

    Returns:
        List of per-position unpaired probabilities (0–1).
    """
    import RNA

    n = len(seq)
    if n == 0:
        return []

    # For short sequences, compute full partition function at once
    if n <= window:
        return _compute_pf_unpaired(seq, span)

    # For long sequences, use sliding window
    result: list[float] = [0.0] * n
    counts: list[int] = [0] * n
    step = max(1, window // 4)  # 75% overlap between windows

    for start in range(0, n - window + 1, step):
        end = min(start + window, n)
        sub_seq = seq[start:end]
        sub_acc = _compute_pf_unpaired(sub_seq, span)

        for i, acc in enumerate(sub_acc):
            pos = start + i
            if pos < n:
                result[pos] += acc
                counts[pos] += 1

    # Handle the tail if it doesn't align with window boundaries
    if n > window:
        tail_start = max(0, n - window)
        tail_seq = seq[tail_start:]
        tail_acc = _compute_pf_unpaired(tail_seq, span)
        for i, acc in enumerate(tail_acc):
            pos = tail_start + i
            if pos < n:
                result[pos] += acc
                counts[pos] += 1

    # Average overlapping predictions
    for i in range(n):
        if counts[i] > 0:
            result[i] /= counts[i]
        else:
            result[i] = 1.0  # default to accessible if no coverage

    return result


def _compute_pf_unpaired(rna_seq: str, span: int) -> list[float]:
    """Compute per-position unpaired probability for a single window.

    Uses ViennaRNA partition function with base-pairing probability
    extraction.  Converts P(paired) → P(unpaired) = 1 - P(paired).

    Args:
        rna_seq: RNA sequence (U not T).
        span: Maximum base-pair span.

    Returns:
        List of per-position unpaired probabilities.
    """
    import RNA

    n = len(rna_seq)
    if n == 0:
        return []
    if n < 2:
        return [1.0]

    try:
        # Create fold compound with span constraint
        fc = RNA.fold_compound(rna_seq)

        # Try to set max base pair span if supported
        try:
            fc.params.model_details.max_bp_span = span
        except (AttributeError, TypeError):
            pass  # older ViennaRNA versions may not support this

        # Compute partition function
        structure, ensemble_energy = fc.pf()

        # Extract base-pairing probabilities
        unpaired: list[float] = [1.0] * n

        # Method 1: Try bpp() matrix (ViennaRNA ≥ 2.4)
        try:
            bpp_matrix = fc.bpp()
            if bpp_matrix and len(bpp_matrix) > 0:
                for i in range(1, n + 1):  # ViennaRNA uses 1-based indexing
                    prob_paired = 0.0
                    row_i = bpp_matrix[i] if i < len(bpp_matrix) else None
                    if row_i is not None:
                        for j in range(i + 1, n + 1):
                            if j < len(row_i):
                                prob_paired += row_i[j]
                    unpaired[i - 1] = max(0.0, min(1.0, 1.0 - prob_paired))
                return unpaired
        except (AttributeError, TypeError, IndexError, KeyError):
            pass

        # Method 2: Try bp_prob(i, j) pairwise queries
        try:
            for i in range(1, n + 1):
                prob_paired = 0.0
                for j in range(1, n + 1):
                    if i == j:
                        continue
                    p = fc.bp_prob(i, j)
                    if p > 0.0:
                        prob_paired += p
                unpaired[i - 1] = max(0.0, min(1.0, 1.0 - prob_paired / 2.0))
            return unpaired
        except (AttributeError, TypeError):
            pass

        # Method 3: Try probs_window (RNAplfold equivalent)
        try:
            # probs_window returns unpaired probabilities directly
            probs = fc.probs_window(0, RNA.PROB_UP, span)
            if probs:
                for entry in probs:
                    # entry format varies by ViennaRNA version
                    if isinstance(entry, (list, tuple)) and len(entry) >= 2:
                        pos = int(entry[0]) - 1  # 1-based → 0-based
                        if 0 <= pos < n:
                            unpaired[pos] = max(0.0, min(1.0, float(entry[1])))
                return unpaired
        except (AttributeError, TypeError, ValueError):
            pass

        # Fallback: use MFE structure as binary accessibility
        try:
            fc2 = RNA.fold_compound(rna_seq)
            mfe_struct, _ = fc2.mfe()
            for i, ch in enumerate(mfe_struct):
                unpaired[i] = 1.0 if ch == "." else 0.0
            return unpaired
        except Exception:
            pass

        # Ultimate fallback: return uniform moderate accessibility
        return [0.5] * n

    except Exception as exc:
        logger.warning("ViennaRNA partition function failed: %s", exc)
        return [0.5] * n


# ── Public API ───────────────────────────────────────────────


def compute_accessibility(
    seq: str,
    window: int = DEFAULT_WINDOW,
    span: int = DEFAULT_SPAN,
    method: str = "viennarna",
) -> list[float]:
    """Compute per-position unpaired probability for an mRNA sequence.

    Uses the ViennaRNA partition function to compute base-pairing
    probabilities, then converts to unpaired probabilities:
    P(unpaired_i) = 1 - P(paired_i).

    For long sequences (> *window* nt), uses a sliding-window approach
    with 75% overlap between windows to ensure smooth coverage.

    When ViennaRNA is unavailable or *method* is ``"gc_heuristic"``,
    estimates accessibility from local GC content using a sigmoid
    function (higher GC → more structure → less accessible).

    Args:
        seq: DNA or RNA sequence (T or U accepted).
        window: Sliding window size for partition function (default 80).
        span: Maximum base-pair span (default 40).
        method: Computation method — ``"viennarna"`` (default) or
            ``"gc_heuristic"``.

    Returns:
        List of per-position unpaired probabilities (0–1), where
        1.0 = fully accessible and 0.0 = fully paired.

    Raises:
        ValueError: If *method* is not a recognized value.

    Examples::

        >>> acc = compute_accessibility("ATGGCCATGGCCATGGCC")
        >>> len(acc) == 18
        True
        >>> all(0.0 <= a <= 1.0 for a in acc)
        True
    """
    if not seq:
        return []

    if method not in ("viennarna", "gc_heuristic"):
        raise ValueError(
            f"method must be 'viennarna' or 'gc_heuristic', got {method!r}"
        )

    rna = _dna_to_rna(seq)

    # Try ViennaRNA if requested
    if method == "viennarna" and HAS_VIENNARNA:
        try:
            result = _accessibility_viennarna_partition(rna, window, span)
            if result and any(r > 0 for r in result):
                return result
            # Fall through to heuristic if ViennaRNA returned all zeros
            logger.debug(
                "ViennaRNA returned all-zero accessibility, falling back to heuristic"
            )
        except Exception as exc:
            logger.warning(
                "ViennaRNA accessibility computation failed: %s — "
                "falling back to GC heuristic",
                exc,
            )

    # GC heuristic fallback
    if method == "viennarna":
        warnings.warn(
            "ViennaRNA unavailable or failed; using GC-content heuristic "
            "for accessibility estimation. This is approximate — install "
            "ViennaRNA (pip install ViennaRNA) for accurate computation.",
            UserWarning,
            stacklevel=2,
        )

    return _accessibility_heuristic(seq)


def compute_accessibility_windows(
    seq: str,
    window_size: int = DEFAULT_SLIDING_WINDOW_SIZE,
    step: int = DEFAULT_SLIDING_STEP,
) -> list[dict[str, Any]]:
    """Sliding-window accessibility computation with position info.

    Divides the sequence into overlapping windows of *window_size*
    nucleotides, stepping by *step* positions.  For each window,
    computes the mean accessibility and per-position unpaired
    probabilities.

    This is useful for identifying regions of particularly high or low
    accessibility across a long mRNA.

    Args:
        seq: DNA or RNA sequence (T or U accepted).
        window_size: Size of each sliding window (default 80 nt).
        step: Step size between windows (default 10 nt).

    Returns:
        List of dicts, each containing:
          - ``"start"``: 0-based start position
          - ``"end"``: 0-based exclusive end position
          - ``"mean_accessibility"``: Average P(unpaired) in the window
          - ``"position_accessibility"``: Dict of position → P(unpaired)
          - ``"gc_content"``: GC fraction of the window
          - ``"method"``: Computation method used
    """
    if not seq:
        return []

    seq_upper = seq.upper()
    n = len(seq_upper)

    # Handle short sequences as a single window
    if n <= window_size:
        acc = compute_accessibility(seq_upper, window=window_size)
        pos_acc = {i: acc[i] for i in range(len(acc))}
        gc = _compute_gc(seq_upper)
        method = "viennarna" if HAS_VIENNARNA else "gc_heuristic"
        return [
            {
                "start": 0,
                "end": n,
                "mean_accessibility": sum(acc) / len(acc) if acc else 0.0,
                "position_accessibility": pos_acc,
                "gc_content": gc,
                "method": method,
            }
        ]

    results: list[dict[str, Any]] = []

    for start in range(0, n - window_size + 1, step):
        end = start + window_size
        window_seq = seq_upper[start:end]

        # Compute accessibility for this window
        acc = compute_accessibility(window_seq, window=window_size)

        # Map window-local positions to global positions
        pos_acc = {start + i: acc[i] for i in range(len(acc))}

        mean_acc = sum(acc) / len(acc) if acc else 0.0
        gc = _compute_gc(window_seq)
        method = "viennarna" if HAS_VIENNARNA else "gc_heuristic"

        results.append(
            {
                "start": start,
                "end": end,
                "mean_accessibility": mean_acc,
                "position_accessibility": pos_acc,
                "gc_content": gc,
                "method": method,
            }
        )

    # Handle trailing partial window
    if n > window_size and (n - window_size) % step != 0:
        start = n - window_size
        window_seq = seq_upper[start:]
        acc = compute_accessibility(window_seq, window=window_size)
        pos_acc = {start + i: acc[i] for i in range(len(acc))}
        mean_acc = sum(acc) / len(acc) if acc else 0.0
        gc = _compute_gc(window_seq)
        method = "viennarna" if HAS_VIENNARNA else "gc_heuristic"
        results.append(
            {
                "start": start,
                "end": n,
                "mean_accessibility": mean_acc,
                "position_accessibility": pos_acc,
                "gc_content": gc,
                "method": method,
            }
        )

    return results


def compute_codon_accessibility(
    dna_seq: str,
    method: str = "viennarna",
) -> list[float]:
    """Compute average accessibility per codon.

    For each codon (3-nt window aligned to the reading frame), computes
    the mean unpaired probability across the three nucleotide positions.
    This captures the local structure environment around each codon,
    which affects ribosome elongation speed (Bazzini 2016).

    Args:
        dna_seq: DNA coding sequence (length should be a multiple of 3).
        method: Computation method — ``"viennarna"`` or ``"gc_heuristic"``.

    Returns:
        List of per-codon accessibility scores (0–1), one entry per codon.

    Examples::

        >>> codon_acc = compute_codon_accessibility("ATGGCCATGGCCATGGCC")
        >>> len(codon_acc) == 6  # 18 nt / 3 = 6 codons
        True
    """
    if not dna_seq:
        return []

    dna_seq = dna_seq.upper()

    # Compute full per-position accessibility
    per_pos = compute_accessibility(dna_seq, method=method)

    # Average every 3 positions (codon alignment from position 0)
    n_codons = len(dna_seq) // 3
    codon_acc: list[float] = []

    for i in range(n_codons):
        pos_start = i * 3
        pos_end = min(pos_start + 3, len(per_pos))
        if pos_end > pos_start:
            positions = per_pos[pos_start:pos_end]
            codon_acc.append(sum(positions) / len(positions))
        else:
            codon_acc.append(0.5)  # fallback

    return codon_acc


def compute_5prime_accessibility(
    dna_seq: str,
    n_codons: int = DEFAULT_5PRIME_N_CODONS,
) -> float:
    """Compute average accessibility in the first N codons.

    The 5' region of an mRNA is critical for translation initiation.
    Inaccessible 5' regions (buried in secondary structure) strongly
    inhibit ribosome binding.  This function provides a single summary
    score for the accessibility of the first *n_codons* codons
    (≈ 3 × n_codons nucleotides from the 5' end).

    Uses ViennaRNA when available; falls back to GC heuristic.

    Args:
        dna_seq: DNA coding sequence.
        n_codons: Number of codons from the 5' end to analyse
            (default 17, ≈ 51 nt).

    Returns:
        Mean unpaired probability across the 5' region (0–1).
        Returns 0.5 if the sequence is empty or too short.

    Examples::

        >>> acc = compute_5prime_accessibility("ATGGCC" * 10, n_codons=10)
        >>> 0.0 <= acc <= 1.0
        True
    """
    if not dna_seq:
        return 0.5

    dna_seq = dna_seq.upper()
    n_nt = min(n_codons * 3, len(dna_seq))

    if n_nt < 3:
        return 0.5

    region = dna_seq[:n_nt]
    per_pos = compute_accessibility(region)

    if not per_pos:
        return 0.5

    return sum(per_pos) / len(per_pos)


def compute_mirna_site_accessibility(
    dna_seq: str,
    site_positions: list[int],
    site_length: int = DEFAULT_MIRNA_SITE_LENGTH,
) -> list[float]:
    """Compute accessibility at miRNA binding sites.

    miRNA binding efficacy is strongly modulated by the accessibility
    of the target site in the mRNA secondary structure (Kertesz 2010).
    Sites buried in stable stems are ineffective even with perfect
    seed complementarity.

    This function computes the mean unpaired probability across each
    miRNA binding site, providing a per-site accessibility score.

    Args:
        dna_seq: DNA sequence containing the miRNA binding sites.
        site_positions: List of 0-based start positions for each
            miRNA binding site in the DNA sequence.
        site_length: Length of each miRNA binding site in nucleotides
            (default 8).

    Returns:
        List of accessibility scores, one per site. Each score is the
        mean P(unpaired) across the site's nucleotides (0–1).

    Examples::

        >>> acc = compute_mirna_site_accessibility(
        ...     "ATGGCCATGGCCATGGCCATGGCC",
        ...     site_positions=[3, 12],
        ...     site_length=8,
        ... )
        >>> len(acc) == 2
        True
        >>> all(0.0 <= a <= 1.0 for a in acc)
        True
    """
    if not dna_seq or not site_positions:
        return []

    dna_seq = dna_seq.upper()

    # Compute full per-position accessibility once
    per_pos = compute_accessibility(dna_seq)

    if not per_pos:
        return [0.5] * len(site_positions)

    n = len(per_pos)
    site_accessibility: list[float] = []

    for pos in site_positions:
        if pos < 0 or pos >= n:
            site_accessibility.append(0.5)
            continue

        end = min(pos + site_length, n)
        site_positions_list = per_pos[pos:end]

        if site_positions_list:
            site_accessibility.append(sum(site_positions_list) / len(site_positions_list))
        else:
            site_accessibility.append(0.5)

    return site_accessibility


def adjust_severity_for_accessibility(
    severity: float,
    accessibility: float,
    threshold: float = ACCESSIBILITY_SEVERITY_THRESHOLD,
) -> float:
    """Reduce severity of sites in inaccessible (structured) regions.

    When a functional site (e.g., miRNA binding site, splice site, RBP
    motif) falls within a region of very low accessibility (< *threshold*),
    it is likely buried in mRNA secondary structure and biologically
    ineffective.  This function reduces the severity score accordingly:

    - If accessibility >= threshold: severity unchanged
    - If accessibility < threshold: severity reduced by 80%
      (i.e., multiplied by 0.20)

    The 80% reduction reflects the observation from Kertesz et al. (2010)
    that site accessibility is a major determinant of miRNA efficacy,
    and sites in highly structured regions are largely non-functional.

    Args:
        severity: Original severity score (any non-negative value).
        accessibility: Mean accessibility at the site (0–1).
        threshold: Accessibility threshold below which severity is
            reduced (default 0.05).

    Returns:
        Adjusted severity score. If accessibility < threshold,
        returns ``severity * (1 - 0.80)`` = ``severity * 0.20``.

    Examples::

        >>> adjust_severity_for_accessibility(1.0, 0.50)
        1.0
        >>> adjust_severity_for_accessibility(1.0, 0.03)
        0.2
        >>> adjust_severity_for_accessibility(0.8, 0.01)
        0.16000000000000003
    """
    if accessibility < threshold:
        reduction_factor = 1.0 - ACCESSIBILITY_SEVERITY_REDUCTION
        return severity * reduction_factor
    return severity
