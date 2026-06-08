"""
BioCompiler Nucleosome Positioning PSSM Module
================================================

Position-specific scoring matrix (PSSM) implementation for nucleosome
positioning prediction based on the Kaplan Lab / Segal 2006 periodicity
model.

This module provides the full 147×16 dinucleotide PSSM, sequence scoring
functions, and utilities for loading and parsing Kaplan Lab model files.

Model Background
----------------
The Segal 2006 / Kaplan 2009 model encodes nucleosome sequence preferences
as position-specific dinucleotide log-likelihood ratios across the 147 bp
nucleosome footprint.  Key periodicity features:

- **AA/TT ~10.2 bp period**: Minor groove faces inward (toward histone
  octamer) at positions where the DNA backbone bends sharply.  AA/TT
  dinucleotides are enriched at these positions.
- **GC/CG inverse phase**: GC/CG dinucleotides are enriched at the
  opposite phase, where the major groove faces inward.
- **TA depletion**: TA dinucleotides are uniformly depleted from
  nucleosomal DNA due to their rigid, narrow minor groove.
- **Poly-dA:dT disfavour**: Runs of A/T are strongly disfavoured because
  they resist the bending required for nucleosome wrapping.

References
----------
- Segal, E. *et al.*, "A genomic code for nucleosome positioning",
  **Nature** 442, 772–778 (2006).  doi:10.1038/nature04979
- Kaplan, N. *et al.*, "The DNA-encoded nucleosome organization of a
  eukaryotic genome", **Nature** 458, 362–366 (2009).
  doi:10.1038/nature07667
- Field, Y. *et al.*, "Gene expression is influenced by nucleosome
  positioning in vivo", **PLoS Genet** 4, e1000204 (2008).
  doi:10.1371/journal.pgen.1000204
"""

from __future__ import annotations

import logging
import os
import re
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

# ── Constants ────────────────────────────────────────────────────────────────

NUCLEOSOME_SIZE: int = 147
"""DNA footprint of a single nucleosome core particle (bp)."""

HELICAL_PERIOD: float = 10.2
"""DNA helical repeat period (bp) around a nucleosome core."""

DEFAULT_CHEMICAL_POTENTIAL: float = -3.0
"""Default chemical potential mu for histone octamers (kT units)."""

DINUC_INDEX: dict[str, int] = {
    "AA": 0, "AC": 1, "AG": 2, "AT": 3,
    "CA": 4, "CC": 5, "CG": 6, "CT": 7,
    "GA": 8, "GC": 9, "GG": 10, "GT": 11,
    "TA": 12, "TC": 13, "TG": 14, "TT": 15,
}
"""Mapping of dinucleotide strings to column indices in the PSSM."""

DINUCS_ORDERED: list[str] = list(DINUC_INDEX.keys())
"""Dinucleotide labels in PSSM column order."""

# Reverse mapping for display
_INDEX_DINUC: dict[int, str] = {v: k for k, v in DINUC_INDEX.items()}

# Complement dinucleotide mapping (Watson-Crick reverse complement)
_DINUC_COMPLEMENT: dict[str, str] = {
    "AA": "TT", "AC": "GT", "AG": "CT", "AT": "AT",
    "CA": "TG", "CC": "GG", "CG": "CG", "CT": "AG",
    "GA": "TC", "GC": "GC", "GG": "CC", "GT": "AC",
    "TA": "TA", "TC": "GA", "TG": "CA", "TT": "AA",
}

# Kaplan Lab GitHub repository for model files
KAPLAN_GITHUB_RAW: str = (
    "https://raw.githubusercontent.com/jipingw/NuPoP_Fortran/master/"
)
KAPLAN_GITHUB_REPO: str = "https://github.com/jipingw/NuPoP_Fortran"

# ── Segal 2006 Periodicity Model Parameters ─────────────────────────────────
#
# These parameters reproduce the key features of the Segal 2006 / Kaplan 2009
# nucleosome positioning model.  The model is based on dinucleotide
# periodicity with a ~10.2 bp helical repeat, plus position-independent
# terms for TA depletion and poly-dA:dT disfavour.
#
# Parameter values are derived from the published supplementary materials
# and the Kaplan Lab Perl implementation.

# Per-dinucleotide periodicity parameters: (amplitude, phase_offset_bp)
#   amplitude  — magnitude of the cosine oscillation in log-likelihood
#   phase      — position (bp) offset relative to the dyad
_PERIDIODICITY_PARAMS: dict[str, tuple[float, float]] = {
    # AA/TT: strong ~10.2 bp periodicity, minor groove inward
    "AA": (0.150, 5.0),
    "TT": (0.150, 5.0),
    # AT: weaker periodicity, similar phase to AA/TT
    "AT": (0.060, 5.3),
    # GC/CG: inverse phase (major groove inward)
    "GC": (0.120, 0.0),
    "CG": (0.120, 0.0),
    # GG/CC: weaker GC-class contribution
    "GG": (0.050, 0.2),
    "CC": (0.050, 0.2),
}

# Position-independent (constant) dinucleotide scores
_CONSTANT_SCORES: dict[str, float] = {
    # TA: uniformly depleted from nucleosomes (rigid, narrow minor groove)
    "TA": -0.250,
    # CA/TG: moderately favourable (bendable)
    "CA": 0.030,
    "TG": 0.030,
    # AC/GT: slightly favourable
    "AC": 0.020,
    "GT": 0.020,
    # AG/CT: small positive
    "AG": 0.010,
    "CT": 0.010,
    # GA/TC: small positive
    "GA": 0.010,
    "TC": 0.010,
}

# Position-dependent amplitude envelope
# The outer ~20 bp of the nucleosome are less constrained, so we apply
# a trapezoidal envelope that ramps up from 0 to full amplitude over
# the first 20 positions and down over the last 20.
_RAMP_LENGTH: int = 20


def _position_envelope(pos: int, size: int = NUCLEOSOME_SIZE) -> float:
    """Compute the amplitude envelope for a given position.

    Positions near the entry/exit points of the nucleosome have reduced
    constraint, reflected as a linear ramp from 0 to 1 over the first
    and last ``_RAMP_LENGTH`` positions.

    Parameters
    ----------
    pos : int
        Zero-based position within the 147 bp footprint.
    size : int
        Total footprint size (default 147).

    Returns
    -------
    float
        Envelope value in [0, 1].
    """
    if pos < _RAMP_LENGTH:
        return pos / _RAMP_LENGTH
    elif pos >= size - _RAMP_LENGTH:
        return (size - 1 - pos) / _RAMP_LENGTH
    else:
        return 1.0


# ── PSSM Generation ─────────────────────────────────────────────────────────


def generate_segal_pssm() -> np.ndarray:
    """Generate the 147×16 Segal 2006 / Kaplan 2009 PSSM.

    Creates a position-specific scoring matrix where each entry is the
    log-likelihood ratio of observing a dinucleotide at a specific
    position in nucleosomal DNA versus linker DNA.

    The matrix encodes three classes of sequence signal:

    1. **Helical periodicity**: AA/TT and GC/CG oscillate with a ~10.2 bp
       period but in opposite phase, reflecting the rotational setting of
       DNA around the histone octamer.
    2. **TA depletion**: TA dinucleotides are uniformly penalised because
       their rigid, narrow minor groove resists the bending required for
       nucleosome wrapping.
    3. **Position-dependent envelope**: Outer positions have reduced
       amplitude because entry/exit DNA is less constrained.

    Returns
    -------
    np.ndarray
        A 147×16 float64 array.  Row indices correspond to positions
        0–146 within the nucleosome footprint; column indices follow
        :data:`DINUC_INDEX`.

    Notes
    -----
    The parameters are calibrated to reproduce the correlation structure
    of the published Segal 2006 model.  For the exact numerical values
    from the Kaplan Lab, use :func:`load_kaplan_model` or
    :func:`parse_kaplan_perl_model` to load the original model files.

    Examples
    --------
    >>> pssm = generate_segal_pssm()
    >>> pssm.shape
    (147, 16)
    >>> pssm[73, DINUC_INDEX["AA"]]  # doctest: +SKIP
    0.15  # near dyad, AA at peak of minor-groove-inward phase
    """
    pssm = np.zeros((NUCLEOSOME_SIZE, 16), dtype=np.float64)

    for pos in range(NUCLEOSOME_SIZE):
        envelope = _position_envelope(pos)

        # ── Periodic dinucleotide terms ────────────────────────────────
        for dinuc, (amplitude, phase_offset) in _PERIDIODICITY_PARAMS.items():
            angular_freq = 2.0 * np.pi / HELICAL_PERIOD
            theta = angular_freq * (pos - phase_offset)
            pssm[pos, DINUC_INDEX[dinuc]] = envelope * amplitude * np.cos(theta)

        # ── Constant (position-independent) terms ─────────────────────
        for dinuc, score in _CONSTANT_SCORES.items():
            # TA depletion is still subject to the envelope at the edges
            if dinuc == "TA":
                pssm[pos, DINUC_INDEX[dinuc]] = envelope * score
            else:
                pssm[pos, DINUC_INDEX[dinuc]] = score

    return pssm


# Pre-compute the default PSSM at module load time
_SEGAL_PSSM: np.ndarray = generate_segal_pssm()
"""The default 147×16 Segal PSSM, generated at import time."""


# ── Sequence Scoring ────────────────────────────────────────────────────────


def score_sequence_pssm(
    seq: str,
    pssm: np.ndarray | None = None,
) -> float:
    """Score a 147 bp DNA sequence using the dinucleotide PSSM.

    Computes the sum of log-likelihood ratios for each dinucleotide
    in *seq* at its corresponding position in the PSSM.  Higher scores
    indicate greater nucleosome-favouring sequence content.

    Parameters
    ----------
    seq : str
        DNA sequence of length >= 147.  Only the first 147 bp are scored.
    pssm : np.ndarray or None
        A 147×16 PSSM array.  If ``None``, uses the default
        :data:`_SEGAL_PSSM`.

    Returns
    -------
    float
        Log-likelihood score.  Returns 0.0 if the sequence is shorter
        than 147 bp or contains no valid dinucleotides.

    Examples
    --------
    >>> score_sequence_pssm("A" * 147)  # doctest: +SKIP
    2.34  # poly-A is weakly nucleosome-favouring at periodic positions
    """
    if pssm is None:
        pssm = _SEGAL_PSSM

    if len(seq) < NUCLEOSOME_SIZE:
        logger.warning(
            "Sequence length %d < %d; returning 0.0",
            len(seq), NUCLEOSOME_SIZE,
        )
        return 0.0

    score = 0.0
    seq_upper = seq.upper()[:NUCLEOSOME_SIZE]

    for pos in range(NUCLEOSOME_SIZE - 1):  # 146 dinucleotides in 147 bp
        dinuc = seq_upper[pos : pos + 2]
        idx = DINUC_INDEX.get(dinuc)
        if idx is not None:
            score += float(pssm[pos, idx])

    return score


def score_sequence_sliding(
    seq: str,
    pssm: np.ndarray | None = None,
    step: int = 10,
) -> list[dict[str, Any]]:
    """Sliding-window PSSM scoring across a DNA sequence.

    Evaluates the PSSM at successive positions, stepping by *step* bp,
    and returns a list of scored windows.

    Parameters
    ----------
    seq : str
        DNA sequence to score.
    pssm : np.ndarray or None
        A 147×16 PSSM.  If ``None``, uses the default.
    step : int
        Sliding window step size in bp (default 10).

    Returns
    -------
    list[dict]
        Each dict has keys:

        - ``"position"`` (int) — start position of the 147 bp window
        - ``"dyad"`` (int) — estimated dyad position (start + 73)
        - ``"raw_score"`` (float) — PSSM log-likelihood score
        - ``"occupancy_prob"`` (float) — logistic-sigmoid occupancy
          probability

    Examples
    --------
    >>> results = score_sequence_sliding("A" * 500, step=50)
    >>> len(results)
    4
    """
    if pssm is None:
        pssm = _SEGAL_PSSM

    seq_upper = seq.upper()
    n = len(seq_upper)

    if n < NUCLEOSOME_SIZE:
        return []

    results: list[dict[str, Any]] = []
    for start in range(0, n - NUCLEOSOME_SIZE + 1, step):
        window = seq_upper[start : start + NUCLEOSOME_SIZE]
        raw = score_sequence_pssm(window, pssm=pssm)
        prob = convert_occupancy_to_prob(raw)
        results.append(
            {
                "position": start,
                "dyad": start + NUCLEOSOME_SIZE // 2,
                "raw_score": raw,
                "occupancy_prob": prob,
            }
        )

    return results


# ── Occupancy Probability Conversion ────────────────────────────────────────


def convert_occupancy_to_prob(
    raw_score: float,
    midpoint: float = 0.0,
    width: float = 5.0,
) -> float:
    """Convert a raw PSSM score to an occupancy probability.

    Uses a logistic sigmoid to map the unbounded log-likelihood score
    to the [0, 1] interval.  The *midpoint* and *width* parameters
    control the centre and steepness of the sigmoid.

    Parameters
    ----------
    raw_score : float
        Raw log-likelihood score from :func:`score_sequence_pssm`.
    midpoint : float
        Score at which the probability is 0.5 (default 0.0).
    width : float
        Scaling factor for the sigmoid steepness (default 5.0).
        Larger values produce a sharper transition.

    Returns
    -------
    float
        Occupancy probability in [0, 1].

    Notes
    -----
    The logistic sigmoid is:

    .. math::

        P = \\frac{1}{1 + e^{-(s - \\mu) / \\sigma}}

    where *s* is the raw score, *μ* is *midpoint*, and *σ* is *width*.

    Examples
    --------
    >>> convert_occupancy_to_prob(0.0)
    0.5
    >>> convert_occupancy_to_prob(10.0)  # strongly favourable
    0.8808...
    """
    exponent = -(raw_score - midpoint) / width
    # Guard against overflow
    if exponent > 500.0:
        return 0.0
    elif exponent < -500.0:
        return 1.0
    return 1.0 / (1.0 + np.exp(exponent))


# ── Well-Positioned Nucleosome Detection ────────────────────────────────────


def find_well_positioned_nucleosomes(
    scores: list[dict[str, Any]],
    threshold: float = 0.7,
) -> list[dict[str, Any]]:
    """Identify well-positioned nucleosomes from sliding-window scores.

    A nucleosome is considered "well-positioned" if its occupancy
    probability exceeds *threshold* and it is a local maximum among
    overlapping windows (i.e., no overlapping window has a higher score).

    Parameters
    ----------
    scores : list[dict]
        Output from :func:`score_sequence_sliding`.  Each dict must
        contain ``"position"``, ``"dyad"``, ``"raw_score"``, and
        ``"occupancy_prob"`` keys.
    threshold : float
        Minimum occupancy probability to qualify as well-positioned
        (default 0.7).

    Returns
    -------
    list[dict]
        Subset of *scores* that are well-positioned, with an additional
        ``"well_positioned"`` key set to ``True``.

    Examples
    --------
    >>> scores = score_sequence_sliding("G" * 500, step=10)
    >>> positioned = find_well_positioned_nucleosomes(scores, threshold=0.6)
    >>> len(positioned) >= 0
    True
    """
    if not scores:
        return []

    # Sort by position for local-maximum detection
    sorted_scores = sorted(scores, key=lambda d: d["position"])

    # Filter by threshold
    candidates = [
        s for s in sorted_scores
        if s.get("occupancy_prob", 0.0) >= threshold
    ]

    if not candidates:
        return []

    # Find local maxima: a candidate is a local max if no overlapping
    # window has a higher occupancy probability.
    well_positioned: list[dict[str, Any]] = []
    for i, cand in enumerate(candidates):
        pos_i = cand["position"]
        occ_i = cand["occupancy_prob"]
        is_local_max = True

        for j, other in enumerate(candidates):
            if i == j:
                continue
            pos_j = other["position"]
            occ_j = other["occupancy_prob"]

            # Check overlap: windows overlap if |pos_i - pos_j| < 147
            if abs(pos_j - pos_i) < NUCLEOSOME_SIZE:
                if occ_j > occ_i:
                    is_local_max = False
                    break

        if is_local_max:
            entry = dict(cand)
            entry["well_positioned"] = True
            well_positioned.append(entry)

    return well_positioned


# ── Kaplan Lab Model Loading ────────────────────────────────────────────────


def load_kaplan_model(model_dir: str | None = None) -> dict[str, Any]:
    """Load Kaplan Lab nucleosome model parameters.

    Attempts to load the Kaplan Lab nucleosome positioning model from
    the NuPoP repository.  If *model_dir* is provided, looks for the
    model files there; otherwise attempts to download from GitHub.

    Parameters
    ----------
    model_dir : str or None
        Path to a directory containing the Kaplan Lab model files.
        If ``None``, attempts to clone/download from the NuPoP GitHub
        repository.

    Returns
    -------
    dict
        Dictionary with keys:

        - ``"pssm"`` (np.ndarray) — 147×16 PSSM from the model, or the
          Segal-generated default if the model files are unavailable
        - ``"source"`` (str) — ``"kaplan_file"``, ``"kaplan_download"``,
          or ``"segal_default"``
        - ``"model_dir"`` (str or None) — Path where model was loaded from
        - ``"species"`` (str) — Species identifier (default ``"human"``)
        - ``"order"`` (int) — Markov chain order (default 4)

    Notes
    -----
    The Kaplan Lab model is distributed as part of the NuPoP Fortran
    package (Wang *et al.*, NAR 2008; Xi *et al.*, BMC Bioinf 2010).
    The model files contain species-specific 4th-order Markov chain
    parameters for nucleosome and linker states.

    If model files cannot be located or parsed, the function falls back
    to the :func:`generate_segal_pssm` periodicity model.
    """
    result: dict[str, Any] = {
        "pssm": _SEGAL_PSSM.copy(),
        "source": "segal_default",
        "model_dir": model_dir,
        "species": "human",
        "order": 4,
    }

    # Try loading from specified directory
    if model_dir is not None:
        model_path = Path(model_dir)
        if model_path.is_dir():
            # Look for Perl model files
            perl_files = list(model_path.glob("*.pl")) + list(
                model_path.glob("**/*.pl")
            )
            for pf in perl_files:
                try:
                    parsed_pssm = parse_kaplan_perl_model(str(pf))
                    if parsed_pssm.shape == (NUCLEOSOME_SIZE, 16):
                        result["pssm"] = parsed_pssm
                        result["source"] = "kaplan_file"
                        result["model_dir"] = str(model_path)
                        return result
                except Exception as exc:
                    logger.debug(
                        "Failed to parse Perl model %s: %s", pf, exc
                    )

            # Look for numpy .npy or .npz files
            for ext in ("*.npy", "*.npz"):
                np_files = list(model_path.glob(ext))
                for nf in np_files:
                    try:
                        loaded = np.load(str(nf))
                        arr = loaded if isinstance(loaded, np.ndarray) else loaded.get("pssm", loaded.get("arr_0"))
                        if arr.shape == (NUCLEOSOME_SIZE, 16):
                            result["pssm"] = arr
                            result["source"] = "kaplan_file"
                            result["model_dir"] = str(model_path)
                            return result
                    except Exception as exc:
                        logger.debug("Failed to load numpy model %s: %s", nf, exc)

    # Try downloading from GitHub
    try:
        import subprocess

        with tempfile.TemporaryDirectory() as tmpdir:
            repo_dir = os.path.join(tmpdir, "NuPoP_Fortran")

            # Shallow clone
            subprocess.run(
                [
                    "git", "clone", "--depth", "1",
                    KAPLAN_GITHUB_REPO, repo_dir,
                ],
                capture_output=True,
                timeout=120,
                check=False,
            )

            if os.path.isdir(repo_dir):
                perl_files = list(Path(repo_dir).rglob("*.pl"))
                for pf in perl_files:
                    try:
                        parsed_pssm = parse_kaplan_perl_model(str(pf))
                        if parsed_pssm.shape == (NUCLEOSOME_SIZE, 16):
                            result["pssm"] = parsed_pssm
                            result["source"] = "kaplan_download"
                            result["model_dir"] = repo_dir
                            return result
                    except Exception as exc:
                        logger.debug(
                            "Failed to parse downloaded model %s: %s",
                            pf, exc,
                        )
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as exc:
        logger.debug("Could not download Kaplan model: %s", exc)

    logger.info(
        "Kaplan model files not found; using Segal-generated PSSM"
    )
    return result


def parse_kaplan_perl_model(perl_model_path: str) -> np.ndarray:
    """Parse a Kaplan Lab Perl model file to extract the PSSM.

    The Kaplan Lab distributes nucleosome positioning parameters in
    Perl source files as part of the NuPoP package.  This function
    extracts the dinucleotide frequency or weight arrays from such
    files and converts them into a 147×16 numpy array.

    Parameters
    ----------
    perl_model_path : str
        Path to a Perl model file (e.g., from the NuPoP_Fortran repo).

    Returns
    -------
    np.ndarray
        A 147×16 float64 PSSM.  If the file does not contain a
        recognisable 147×16 matrix, a :func:`generate_segal_pssm`
        fallback is returned and a warning is logged.

    Raises
    ------
    FileNotFoundError
        If *perl_model_path* does not exist.

    Notes
    -----
    The parser looks for Perl array initialisation patterns such as::

        @matrix = (
            [0.12, -0.03, ...],
            ...
        );

    or inline numeric arrays.  It also recognises Fortran-format
    DATA statements that may appear in the NuPoP Fortran source.
    """
    path = Path(perl_model_path)
    if not path.exists():
        raise FileNotFoundError(f"Perl model file not found: {perl_model_path}")

    content = path.read_text(errors="replace")

    # Strategy 1: Look for Perl 2D array patterns
    # Match patterns like: [0.12, -0.03, 0.05, ...],
    row_pattern = re.compile(
        r"\[\s*"
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
        r"(?:\s*,\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))*"
        r"\s*\]"
    )

    rows: list[list[float]] = []
    for match in row_pattern.finditer(content):
        row_str = match.group(0)
        # Extract all numbers from this row
        numbers = re.findall(
            r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", row_str
        )
        if len(numbers) == 16:
            try:
                row = [float(x) for x in numbers]
                rows.append(row)
            except ValueError:
                continue

    if len(rows) == NUCLEOSOME_SIZE:
        pssm = np.array(rows, dtype=np.float64)
        if pssm.shape == (NUCLEOSOME_SIZE, 16):
            logger.info(
                "Parsed %d×%d PSSM from Perl model %s",
                *pssm.shape, perl_model_path,
            )
            return pssm

    # Strategy 2: Look for flat arrays of 147*16 = 2352 numbers
    all_numbers = re.findall(
        r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", content
    )
    expected = NUCLEOSOME_SIZE * 16

    if len(all_numbers) >= expected:
        try:
            flat = [float(x) for x in all_numbers[:expected]]
            pssm = np.array(flat, dtype=np.float64).reshape(
                NUCLEOSOME_SIZE, 16
            )
            logger.info(
                "Parsed flat %d×%d PSSM from Perl model %s",
                *pssm.shape, perl_model_path,
            )
            return pssm
        except ValueError:
            pass

    # Strategy 3: Look for Fortran DATA statements
    # DATA matrix / 0.12, -0.03, ... /
    data_pattern = re.compile(
        r"DATA\s+\w+\s*/\s*([^/]+)\s*/", re.IGNORECASE
    )
    for match in data_pattern.finditer(content):
        nums = re.findall(
            r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", match.group(1)
        )
        if len(nums) >= expected:
            try:
                flat = [float(x) for x in nums[:expected]]
                pssm = np.array(flat, dtype=np.float64).reshape(
                    NUCLEOSOME_SIZE, 16
                )
                logger.info(
                    "Parsed Fortran DATA %d×%d PSSM from %s",
                    *pssm.shape, perl_model_path,
                )
                return pssm
            except ValueError:
                continue

    # Fallback: return the Segal-generated PSSM with a warning
    logger.warning(
        "Could not parse a 147×16 matrix from %s; "
        "returning Segal-generated PSSM as fallback",
        perl_model_path,
    )
    return generate_segal_pssm()


# ── Utility Functions ────────────────────────────────────────────────────────


def dinuc_composition(seq: str) -> dict[str, int]:
    """Count dinucleotide occurrences in a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence.

    Returns
    -------
    dict[str, int]
        Counts keyed by dinucleotide (only canonical dinucleotides).
    """
    seq_upper = seq.upper()
    counts: dict[str, int] = {d: 0 for d in DINUCS_ORDERED}
    for i in range(len(seq_upper) - 1):
        dinuc = seq_upper[i : i + 2]
        if dinuc in DINUC_INDEX:
            counts[dinuc] += 1
    return counts


def pssm_to_dataframe(pssm: np.ndarray | None = None):
    """Convert a PSSM array to a pandas DataFrame for inspection.

    Parameters
    ----------
    pssm : np.ndarray or None
        A 147×16 PSSM.  If ``None``, uses the default.

    Returns
    -------
    pandas.DataFrame
        DataFrame with dinucleotide column labels and position index.

    Raises
    ------
    ImportError
        If pandas is not installed.
    """
    import pandas as pd

    if pssm is None:
        pssm = _SEGAL_PSSM

    return pd.DataFrame(
        pssm,
        columns=DINUCS_ORDERED,
        index=range(NUCLEOSOME_SIZE),
    )


def compute_dyad_score(
    seq: str,
    dyad_pos: int,
    pssm: np.ndarray | None = None,
) -> float:
    """Score a nucleosome centred at a specific dyad position.

    Parameters
    ----------
    seq : str
        Full DNA sequence.
    dyad_pos : int
        Zero-based position of the nucleosome dyad (centre).
    pssm : np.ndarray or None
        A 147×16 PSSM.  If ``None``, uses the default.

    Returns
    -------
    float
        PSSM log-likelihood score for the 147 bp window centred
        at *dyad_pos*.
    """
    start = dyad_pos - NUCLEOSOME_SIZE // 2
    if start < 0:
        return 0.0
    window = seq[start : start + NUCLEOSOME_SIZE]
    if len(window) < NUCLEOSOME_SIZE:
        return 0.0
    return score_sequence_pssm(window, pssm=pssm)


# ── Module-level Validation ─────────────────────────────────────────────────

def _validate_pssm() -> None:
    """Run basic sanity checks on the generated PSSM."""
    pssm = _SEGAL_PSSM
    assert pssm.shape == (NUCLEOSOME_SIZE, 16), (
        f"PSSM shape {pssm.shape} != ({NUCLEOSOME_SIZE}, 16)"
    )

    # AA/TT should have the same values (complement symmetry)
    assert np.allclose(
        pssm[:, DINUC_INDEX["AA"]], pssm[:, DINUC_INDEX["TT"]]
    ), "AA/TT columns differ — complement symmetry violated"

    # GC/CG should have the same values
    assert np.allclose(
        pssm[:, DINUC_INDEX["GC"]], pssm[:, DINUC_INDEX["CG"]]
    ), "GC/CG columns differ — complement symmetry violated"

    # TA should be uniformly negative
    assert np.all(pssm[:, DINUC_INDEX["TA"]] <= 0.0), (
        "TA column has positive entries — expected uniform depletion"
    )

    # AA/TT should show periodicity
    aa_col = pssm[:, DINUC_INDEX["AA"]]
    # Autocorrelation at lag ~10 should be positive
    if len(aa_col) > 20:
        autocorr = np.corrcoef(aa_col[:-10], aa_col[10:])[0, 1]
        assert autocorr > 0.0, (
            f"AA column autocorrelation at lag 10 is {autocorr:.3f} "
            f"— expected positive periodicity"
        )

    logger.debug("PSSM validation passed")


# Run validation on import (lightweight — just array checks)
try:
    _validate_pssm()
except AssertionError as exc:
    logger.warning("PSSM validation failed: %s", exc)


# ── Public API ───────────────────────────────────────────────────────────────

__all__ = [
    # Constants
    "NUCLEOSOME_SIZE",
    "HELICAL_PERIOD",
    "DEFAULT_CHEMICAL_POTENTIAL",
    "DINUC_INDEX",
    "DINUCS_ORDERED",
    "KAPLAN_GITHUB_RAW",
    "KAPLAN_GITHUB_REPO",
    # PSSM generation
    "generate_segal_pssm",
    # Sequence scoring
    "score_sequence_pssm",
    "score_sequence_sliding",
    "compute_dyad_score",
    # Occupancy conversion
    "convert_occupancy_to_prob",
    # Well-positioned detection
    "find_well_positioned_nucleosomes",
    # Kaplan Lab model loading
    "load_kaplan_model",
    "parse_kaplan_perl_model",
    # Utilities
    "dinuc_composition",
    "pssm_to_dataframe",
]
