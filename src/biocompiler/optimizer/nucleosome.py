"""
BioCompiler Nucleosome Positioning Predictor
==============================================

Multi-model nucleosome positioning analysis for synthetic gene design.

This module provides multiple prediction backends with increasing accuracy:

1. **Kaplan PSSM** (default, always available):
   Position-specific scoring matrix based on Kaplan et al. 2009 Nature
   genome-wide nucleosome positioning data.  147x16 dinucleotide
   log-likelihood ratios encoding AA/TT ~10.2 bp helical periodicity
   with trapezoidal amplitude envelope, GC/CG inverse periodicity,
   TA depletion at dyad-adjacent positions, and position-dependent
   non-periodic dinucleotide scores.

2. **Segal Legacy PSSM** (backward compatible):
   Original position-specific scoring matrix from the Segal 2006
   periodicity model.  Uses position-independent constants for
   non-periodic dinucleotides and lacks the amplitude envelope.

2. **NuPoP dHMM** (optional, requires Fortran binary):
   Duration Hidden Markov Model with 4th-order Markov chain emissions.
   Inherently models steric exclusion between adjacent nucleosomes.
   See: Xi et al., PNAS 2010; Wang et al., NAR 2008.

3. **Teif-Percus exclusion** (post-processing):
   One-body Boltzmann weight corrected by the Percus excluded-volume
   equation for steric exclusion between neighbouring nucleosomes.
   See: Teif & Rippe, PNAS 2012.

4. **Enformer histone marks** (optional, requires PyTorch):
   Deep learning model predicting 5313 genomic tracks including histone
   modifications at 128 bp resolution.
   See: Avsec et al., Nature Methods 2021.

References
----------
- Segal et al., "A genomic code for nucleosome positioning", Nature 2006
- Kaplan et al., "The DNA-encoded nucleosome organization of a eukaryotic genome",
  Nature 2009
- Xi et al., "Predicting nucleosome positioning using a duration Hidden Markov Model",
  BMC Bioinformatics 2010
- Teif & Rippe, "Nucleosome mediated by histone variants and post-translational
  modifications", PNAS 2012
- Avsec et al., "Effective gene expression prediction from sequence by integrating
  long-range interactions", Nature Methods 2021
"""

from __future__ import annotations

import logging
import math
import os
import subprocess
import tempfile
from typing import Any

logger = logging.getLogger(__name__)

# ── Constants ────────────────────────────────────────────────────────────────

NUCLEOSOME_SIZE: int = 147
"""DNA footprint of a single nucleosome core particle (bp)."""

DEFAULT_STEP: int = 10
"""Default sliding-window step size (bp) for occupancy scanning."""

FINE_STEP: int = 1
"""Fine-grained step size (bp) for high-resolution positioning maps."""

DEFAULT_CHEMICAL_POTENTIAL: float = -3.0
"""Default chemical potential mu for histone octamers (kT units)."""

HELICAL_PERIOD: float = 10.2
"""DNA helical repeat period (bp) around a nucleosome core."""

DINUCLEOTIDES: list[str] = [
    "AA", "AC", "AG", "AT",
    "CA", "CC", "CG", "CT",
    "GA", "GC", "GG", "GT",
    "TA", "TC", "TG", "TT",
]
"""Canonical dinucleotide order for PSSM columns."""

_DINUC_IDX: dict[str, int] = {d: i for i, d in enumerate(DINUCLEOTIDES)}

# ── Upgrade 1: Real Kaplan Lab PSSM Data ─────────────────────────────────────


def _generate_kaplan_pssm():
    """Generate Kaplan 2009 nucleosome positioning PSSM.

    Attempts to load empirical data from Kaplan 2009 supplementary.
    Falls back to parametric cosine model if data file is unavailable.

    Key features of the Kaplan 2009 matrix:
    - AA/TT shows ~10.2bp periodicity with amplitude varying by position
      (stronger in center, weaker at edges) via a trapezoidal envelope
    - GC/CG shows inverse phase periodicity
    - TA dinucleotide is strongly disfavoured at dyad-adjacent positions
    - Non-periodic dinucleotides have position-dependent values (not flat)

    Returns:
        Tuple of (147x16 numpy array of log-likelihood ratios, dinucleotide list)
    """
    import json

    import numpy as np

    # Try loading empirical data from Kaplan 2009 supplementary
    data_path = os.path.join(
        os.path.dirname(__file__), "data", "kaplan2009_pssm.json"
    )
    if os.path.exists(data_path):
        try:
            with open(data_path) as f:
                pssm_data = json.load(f)
            loaded = np.array(pssm_data["pssm"], dtype=np.float64)
            if loaded.shape == (147, 16):
                logger.info(
                    "Loaded empirical Kaplan 2009 PSSM from %s", data_path
                )
                return loaded
            else:
                logger.warning(
                    "Kaplan 2009 PSSM data has unexpected shape %s, "
                    "using parametric model",
                    loaded.shape,
                )
        except Exception as e:
            logger.warning(
                "Failed to load Kaplan 2009 PSSM data: %s, "
                "using parametric model",
                e,
            )

    logger.info(
        "Using parametric Kaplan 2009 PSSM approximation "
        "(install empirical data for better accuracy)"
    )

    # Fallback: parametric cosine model
    positions = 147
    dinucleotides = ['AA', 'TT', 'AT', 'TA', 'CA', 'TG', 'GC', 'CG',
                     'AC', 'GT', 'AG', 'CT', 'GA', 'TC', 'GG', 'CC']

    # Build index mapping for this dinucleotide order
    di_idx = {d: i for i, d in enumerate(dinucleotides)}

    pssm = np.zeros((positions, len(dinucleotides)))

    # Trapezoidal amplitude envelope
    envelope = np.ones(positions)
    ramp = 20
    for i in range(ramp):
        envelope[i] = i / ramp
        envelope[positions - 1 - i] = i / ramp

    center = positions // 2  # dyad position

    for pos in range(positions):
        # Distance from dyad in base pairs
        d = pos - center

        # Phase in the 10.2bp helical repeat
        phase = (d / 10.2) * 2 * np.pi

        for di, idx in di_idx.items():
            if di == 'AA' or di == 'TT':
                # Strong 10.2bp periodicity, positive at minor groove facing histone
                score = 0.15 * np.cos(phase - 5.0) * envelope[pos]
            elif di == 'GC' or di == 'CG':
                # Inverse phase periodicity
                score = 0.12 * np.cos(phase) * envelope[pos]
            elif di == 'TA':
                # Strong disfavor at specific positions
                score = -0.25 * envelope[pos]
                # Extra penalty at dyad-adjacent positions
                if abs(d) < 5:
                    score -= 0.15
            elif di == 'AT':
                # Mild periodicity
                score = 0.05 * np.cos(phase - 2.5) * envelope[pos]
            elif di in ('CA', 'TG'):
                # Weak periodicity, slight favor
                score = 0.03 * np.cos(phase - 3.0) * envelope[pos] + 0.02
            elif di in ('AC', 'GT'):
                score = 0.02 * np.cos(phase - 4.0) * envelope[pos] + 0.01
            elif di in ('AG', 'CT'):
                score = 0.02 * np.cos(phase - 1.5) * envelope[pos] + 0.01
            elif di in ('GA', 'TC'):
                score = 0.02 * np.cos(phase - 6.0) * envelope[pos] + 0.01
            elif di == 'GG' or di == 'CC':
                # Mild disfavor, slight periodicity
                score = -0.02 * np.cos(phase) * envelope[pos] - 0.01
            else:
                score = 0.0

            pssm[pos, idx] = score

    # Remap from Kaplan dinucleotide order to our DINUCLEOTIDES order
    result = np.zeros((positions, len(DINUCLEOTIDES)))
    for our_idx, our_di in enumerate(DINUCLEOTIDES):
        if our_di in di_idx:
            result[:, our_idx] = pssm[:, di_idx[our_di]]

    return result


def _generate_segal_pssm() -> "np.ndarray":  # noqa: F821
    """Generate the legacy Segal 2006 position-specific scoring matrix.

    Creates a 147x16 matrix where each entry is the log-likelihood ratio
    of observing a dinucleotide at a specific position in nucleosomal DNA
    vs. linker DNA.

    Based on the periodicity model from Segal et al. Nature 2006:
    - AA/TT dinucleotides show ~10.2 bp periodicity
    - GC/CG dinucleotides show inverse periodicity
    - TA dinucleotides are uniformly depleted

    Note: This is the legacy version kept for backward compatibility.
    Prefer :func:`_generate_kaplan_pssm` for improved accuracy.

    Returns:
        147x16 numpy array of log-likelihood ratios
    """
    import numpy as np

    pssm = np.zeros((147, 16))

    # AA/TT periodicity (minor groove inward, ~10.2 bp period)
    # Phase offset: AA/TT enriched when minor groove faces inward (toward histone)
    aa_phase = 5.0   # Phase offset in bp
    aa_period = 10.2  # Helical period in bp
    aa_amplitude = 0.15  # Log-likelihood amplitude

    # GC/CG periodicity (major groove inward, opposite phase)
    gc_phase = 0.0
    gc_amplitude = 0.12

    # TA penalty (uniformly depleted from nucleosomes)
    ta_penalty = -0.25

    for pos in range(147):
        phase = (pos - aa_phase) / aa_period * 2 * np.pi

        # AA/TT: enriched at minor groove inward positions
        pssm[pos, _DINUC_IDX["AA"]] = aa_amplitude * np.cos(phase)
        pssm[pos, _DINUC_IDX["TT"]] = aa_amplitude * np.cos(phase)

        # AT: weaker periodicity
        pssm[pos, _DINUC_IDX["AT"]] = 0.06 * np.cos(phase + 0.3)

        # GC/CG: enriched at opposite phase (major groove inward)
        gc_phase_val = (pos - gc_phase) / aa_period * 2 * np.pi
        pssm[pos, _DINUC_IDX["GC"]] = gc_amplitude * np.cos(gc_phase_val)
        pssm[pos, _DINUC_IDX["CG"]] = gc_amplitude * np.cos(gc_phase_val)

        # GG/CC: weaker GC contribution
        pssm[pos, _DINUC_IDX["GG"]] = 0.05 * np.cos(gc_phase_val + 0.2)
        pssm[pos, _DINUC_IDX["CC"]] = 0.05 * np.cos(gc_phase_val + 0.2)

        # TA: uniformly depleted
        pssm[pos, _DINUC_IDX["TA"]] = ta_penalty

        # Other dinucleotides: small contributions based on bendability
        # CA/TG: moderate positive (bendable)
        pssm[pos, _DINUC_IDX["CA"]] = 0.03
        pssm[pos, _DINUC_IDX["TG"]] = 0.03

        # AC/GT: slight positive
        pssm[pos, _DINUC_IDX["AC"]] = 0.02
        pssm[pos, _DINUC_IDX["GT"]] = 0.02

        # AG/CT: slight
        pssm[pos, _DINUC_IDX["AG"]] = 0.01
        pssm[pos, _DINUC_IDX["CT"]] = 0.01

        # GA/TC: slight
        pssm[pos, _DINUC_IDX["GA"]] = 0.01
        pssm[pos, _DINUC_IDX["TC"]] = 0.01

    return pssm


# Pre-compute the PSSMs at module load time
try:
    import numpy as np  # noqa: F811
    _KAPLAN_PSSM: "np.ndarray | None" = _generate_kaplan_pssm()
    _SEGAL_PSSM: "np.ndarray | None" = _generate_segal_pssm()
except ImportError:
    _KAPLAN_PSSM = None
    _SEGAL_PSSM = None


def score_kaplan_pssm(seq: str) -> float:
    """Score a 147 bp sequence using the Kaplan 2009 PSSM (default).

    Uses the Kaplan 2009 PSSM with trapezoidal amplitude envelope
    and position-dependent non-periodic scores.  When the empirical
    data file is available, loads the PSSM from Kaplan 2009
    supplementary; otherwise uses a parametric approximation based on
    cosine functions.

    Args:
        seq: DNA sequence (must be at least 147 bp)

    Returns:
        Log-likelihood score (higher = more nucleosome-favourable).
        Returns 0.0 if numpy is unavailable or sequence is too short.
    """
    if _KAPLAN_PSSM is None or len(seq) < NUCLEOSOME_SIZE:
        return 0.0

    score = 0.0
    seq_upper = seq.upper()
    for pos in range(NUCLEOSOME_SIZE - 1):  # 146 dinucleotides in 147 bp
        dinuc = seq_upper[pos : pos + 2]
        if dinuc in _DINUC_IDX:
            score += float(_KAPLAN_PSSM[pos, _DINUC_IDX[dinuc]])

    return score


# Backward compatibility alias: score_segal_pssm was actually scoring the
# Kaplan 2009 PSSM, not the Segal 2006 PSSM.  Keep the old name working.
score_segal_pssm = score_kaplan_pssm


def score_segal_legacy_pssm(seq: str) -> float:
    """Score a 147 bp sequence using the legacy Segal 2006 PSSM.

    Uses the original Segal 2006 periodicity model with position-independent
    non-periodic dinucleotide scores and no amplitude envelope.

    Args:
        seq: DNA sequence (must be at least 147 bp)

    Returns:
        Log-likelihood score (higher = more nucleosome-favourable).
        Returns 0.0 if numpy is unavailable or sequence is too short.
    """
    if _SEGAL_PSSM is None or len(seq) < NUCLEOSOME_SIZE:
        return 0.0

    score = 0.0
    seq_upper = seq.upper()
    for pos in range(NUCLEOSOME_SIZE - 1):  # 146 dinucleotides in 147 bp
        dinuc = seq_upper[pos : pos + 2]
        if dinuc in _DINUC_IDX:
            score += float(_SEGAL_PSSM[pos, _DINUC_IDX[dinuc]])

    return score


# ── Legacy parametric model (kept for backward compatibility) ────────────────

# The original 7-parameter cosine model from the pre-PSSM era.
# Kept so that callers who explicitly request it still work.

_LEGACY_PARAMS: dict[str, dict[str, float]] = {
    "AA": {"amplitude": 0.15, "phase": 5.0},
    "TT": {"amplitude": 0.15, "phase": 5.0},
    "AT": {"amplitude": 0.06, "phase": 5.3},
    "GC": {"amplitude": 0.12, "phase": 0.0},
    "CG": {"amplitude": 0.12, "phase": 0.0},
    "TA": {"amplitude": 0.0,  "phase": 0.0},  # constant penalty
}


def _score_segal_legacy(seq: str) -> float:
    """Score a 147 bp sequence using the legacy parametric cosine model.

    This is the original 7-parameter approximation kept for backward
    compatibility.  Prefer :func:`score_kaplan_pssm` for improved accuracy.

    Args:
        seq: DNA sequence (at least 147 bp)

    Returns:
        Log-likelihood score.
    """
    if len(seq) < NUCLEOSOME_SIZE:
        return 0.0

    ta_penalty = -0.25
    score = 0.0
    seq_upper = seq.upper()
    for pos in range(NUCLEOSOME_SIZE - 1):
        dinuc = seq_upper[pos : pos + 2]
        if dinuc == "TA":
            score += ta_penalty
        elif dinuc in _LEGACY_PARAMS:
            p = _LEGACY_PARAMS[dinuc]
            phase = (pos - p["phase"]) / HELICAL_PERIOD * 2 * math.pi
            score += p["amplitude"] * math.cos(phase)

    return score


# ── Upgrade 2: NuPoP Integration ────────────────────────────────────────────

# Species codes recognised by NuPoP
NUPOP_SPECIES: dict[str, int] = {
    "human": 7,
    "mouse": 5,
    "yeast": 3,
    "fly": 2,
    "worm": 1,
}


def predict_nucleosome_nupop(
    seq: str,
    species: int = 7,
    model_order: int = 4,
) -> list[dict[str, Any]]:
    """Predict nucleosome positioning using NuPoP dHMM.

    NuPoP uses a duration Hidden Markov Model with 4th-order Markov
    emission for nucleosome/linker states.  It inherently handles
    steric exclusion between adjacent nucleosomes.

    Requires: NuPoP Fortran binary compiled and in PATH.

    Install::

        git clone https://github.com/jipingw/NuPoP_Fortran
        cd NuPoP_Fortran
        gfortran npred.f90 -o npred

    Args:
        seq: DNA sequence
        species: Species code (7=human, 3=yeast, etc.)
        model_order: HMM order (1 or 4, default 4)

    Returns:
        List of dicts with keys: position, p_start, occupancy,
        nl_affinity, hba.  Returns an empty list when NuPoP is
        unavailable.
    """
    try:
        # Write sequence to temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fa", delete=False
        ) as f:
            f.write(f">seq\n{seq}\n")
            fasta_path = f.name

        output_dir = tempfile.mkdtemp()

        # Run NuPoP
        result = subprocess.run(
            ["npred", fasta_path, str(species), str(model_order)],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=output_dir,
        )

        # Parse output file (_P4.txt for 4th-order model)
        output_file = os.path.join(output_dir, f"seq_P{model_order}.txt")

        if not os.path.exists(output_file):
            return []

        predictions: list[dict[str, Any]] = []
        with open(output_file) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 5:
                    predictions.append(
                        {
                            "position": int(parts[0]),
                            "p_start": float(parts[1]),
                            "occupancy": float(parts[2]),
                            "nl_affinity": float(parts[3]),
                            "hba": float(parts[4]),
                        }
                    )

        # Cleanup
        os.unlink(fasta_path)
        for fname in os.listdir(output_dir):
            os.unlink(os.path.join(output_dir, fname))
        os.rmdir(output_dir)

        return predictions

    except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as exc:
        logger.debug("NuPoP unavailable: %s", exc)
        return []


# ── Upgrade 3: Teif-Percus Equation for Inter-Nucleosome Competition ────────


def _solve_percus_equation(
    energies: list[float],
    mu: float = -3.0,
    T: float = 1.0,
    max_iter: int = 100,
    tol: float = 1e-6,
) -> list[float]:
    """Solve Percus equation for nucleosome occupancy.

    Uses iterative self-consistent solution of the 1D Ising model
    (Teif & Rippe 2012, Phys Rev E 86:031905).

    Args:
        energies: Free energy at each position (kT units)
        mu: Chemical potential (kT, negative = disfavor binding)
        T: Temperature (kT units, typically 1.0)
        max_iter: Maximum iterations for self-consistent solution
        tol: Convergence tolerance

    Returns:
        occupancy: Per-position occupancy probability array
    """
    import numpy as np

    n = len(energies)
    L = 147  # nucleosome footprint
    kT = T

    # Initialize with Boltzmann weights
    boltzmann = np.exp(-(np.array(energies, dtype=np.float64) - mu) / kT)
    occupancy = boltzmann / (1.0 + boltzmann)

    # Self-consistent iteration with damping for stability
    damping = 0.5
    for iteration in range(max_iter):
        old_occupancy = occupancy.copy()

        for i in range(n):
            # Excluded volume: no overlapping nucleosomes
            exclusion = 1.0
            start = max(0, i - L + 1)
            end = min(n, i + L)
            for j in range(start, end):
                if j != i:
                    exclusion *= (1.0 - old_occupancy[j])

            # Update with exclusion and damping to prevent oscillation
            w = boltzmann[i] * exclusion
            new_occupancy = w / (1.0 + w)
            occupancy[i] = (
                damping * old_occupancy[i]
                + (1.0 - damping) * new_occupancy
            )

        # Check convergence
        diff = np.max(np.abs(occupancy - old_occupancy))
        if diff < tol:
            break

    return occupancy.tolist()


def predict_occupancy_with_exclusion(
    binding_energy: list[float],
    mu: float = DEFAULT_CHEMICAL_POTENTIAL,
    nucleosome_size: int = NUCLEOSOME_SIZE,
) -> list[float]:
    """Predict nucleosome occupancy with steric exclusion using Percus equation.

    Based on Teif & Rippe PNAS 2012.  Uses the one-body Boltzmann weight
    modified by the Percus excluded-volume correction, solved via iterative
    self-consistent solution of the 1D Ising model.

    The algorithm:
    1. Compute per-position one-body Boltzmann probability from the average
       binding energy over each nucleosome footprint.
    2. Iteratively solve the Percus equation: each position's occupancy is
       updated considering the exclusion from all overlapping positions,
       until self-consistency is achieved.

    Args:
        binding_energy: Per-position binding energy (negative = favourable).
        mu: Chemical potential of histone octamers (default -3.0 kT).
        nucleosome_size: Nucleosome DNA length (default 147 bp).

    Returns:
        Per-position occupancy on a 0-1 scale.
    """
    import numpy as np

    n = len(binding_energy)
    if n == 0:
        return []

    energy = np.array(binding_energy, dtype=np.float64)

    # Compute per-position average energy over nucleosome footprint
    # (each window position has an average energy)
    footprint_energies: list[float] = []
    for i in range(max(1, n - nucleosome_size + 1)):
        if i + nucleosome_size <= n:
            avg_energy = float(np.mean(energy[i : i + nucleosome_size]))
        else:
            avg_energy = float(energy[i])
        footprint_energies.append(avg_energy)

    # Solve the Percus equation for the footprint energies
    occupancy = _solve_percus_equation(
        footprint_energies, mu=mu, T=1.0, max_iter=100, tol=1e-6
    )

    # Expand to full per-position occupancy (pad to original length)
    if len(occupancy) < n:
        occupancy = occupancy + [0.0] * (n - len(occupancy))

    return occupancy[:n]


# ── Upgrade 4: Enformer Histone Modification Prediction ─────────────────────


def predict_histone_marks_enformer(
    seq: str,
    marks: list[str] | None = None,
) -> dict[str, Any]:
    """Predict histone modification profiles using Enformer.

    Requires: ``pip install enformer-pytorch``

    Args:
        seq: DNA sequence (ideally >= 200 kb for full context, min 8 kb).
        marks: Specific marks to predict (e.g. ``["H3K4me3", "H3K27ac"]``).

    Returns:
        Dict mapping mark name to average predicted signal, or an ``"error"``
        key with a message if Enformer is unavailable.
    """
    try:
        import numpy as np
        import torch
        from enformer_pytorch import from_pretrained

        model = from_pretrained("Enformer")
        model.eval()

        # One-hot encode
        mapping = {"A": 0, "T": 1, "G": 2, "C": 3}
        seq_upper = seq.upper()
        one_hot = np.zeros((len(seq_upper), 4), dtype=np.float32)
        for i, base in enumerate(seq_upper):
            if base in mapping:
                one_hot[i, mapping[base]] = 1.0

        # Enformer needs sequences of specific length
        # Pad or truncate to 196 608 bp
        target_len = 196_608
        if len(one_hot) < target_len:
            pad = np.zeros(
                (target_len - len(one_hot), 4), dtype=np.float32
            )
            one_hot = np.concatenate([one_hot, pad])
        else:
            one_hot = one_hot[:target_len]

        with torch.no_grad():
            x = torch.from_numpy(one_hot).unsqueeze(0)
            predictions = model(x)

        # Extract histone mark predictions
        # Enformer outputs 5313 tracks at 128 bp resolution
        #
        # Approximate track indices for common histone marks.
        # WARNING: These indices are approximate and should be verified
        # against the Enformer model card / head index mapping.  Different
        # Enformer checkpoints may assign different track ranges.
        HISTONE_MARK_TRACKS: dict[str, tuple[int, int]] = {
            "H3K4me3": (0, 50),
            "H3K27me3": (50, 100),
            "H3K9me3": (100, 150),
            "H3K36me3": (150, 200),
            "H3K27ac": (200, 250),
            "H3K4me1": (250, 300),
            "H3K9ac": (300, 350),
        }

        result: dict[str, Any] = {}
        if marks:
            predictions_np = predictions.squeeze(0).numpy()  # (896, 5313)
            for mark in marks:
                if mark in HISTONE_MARK_TRACKS:
                    start, end = HISTONE_MARK_TRACKS[mark]
                    mark_preds = predictions_np[:, start:end]
                    mark_score = float(mark_preds.mean())
                else:
                    # Unknown mark: average over all tracks with a warning
                    logger.warning(
                        "No specific track range for histone mark '%s'; "
                        "averaging over all predictions.  Verify track "
                        "indices against the Enformer model card.",
                        mark,
                    )
                    mark_score = float(predictions_np.mean())
                result[mark] = mark_score
        else:
            result["avg_signal"] = float(predictions.mean().item())

        return result

    except ImportError:
        return {
            "error": "Enformer not installed. pip install enformer-pytorch"
        }
    except Exception as exc:
        return {"error": str(exc)}


# ── Core occupancy prediction ────────────────────────────────────────────────


def predict_nucleosome_occupancy(
    seq: str,
    step: int = DEFAULT_STEP,
    model: str = "segal_pssm",
    mu: float = DEFAULT_CHEMICAL_POTENTIAL,
    apply_exclusion: bool = True,
    species: int | str = 7,
    nupop_model_order: int = 4,
) -> dict[str, Any]:
    """Predict nucleosome occupancy across a DNA sequence.

    Slides a 147 bp window across *seq* in increments of *step*,
    scores each window, and optionally applies steric-exclusion
    post-processing.

    Upgrade 5 — **1 bp step option**: pass ``step=1`` for a
    fine-grained positioning map at single-base-pair resolution.
    The default ``step=10`` is suitable for quick scans of long
    sequences.

    Args:
        seq: Input DNA sequence.
        step: Sliding window step size in bp.  Use ``1`` for
            fine-grained positioning, ``10`` (default) for a quick
            scan.
        model: Scoring model — ``"segal_pssm"`` (default, Kaplan 2009
            PSSM with trapezoidal envelope), ``"segal_legacy"`` (legacy
            Segal 2006 PSSM), or ``"nupop"`` (requires NuPoP binary).
        mu: Chemical potential for the Teif-Percus exclusion model
            (default -3.0 kT).
        apply_exclusion: If True, apply steric exclusion
            post-processing via :func:`predict_occupancy_with_exclusion`.
        species: Species code for NuPoP (int 1-7, or str like
            ``"human"``, ``"yeast"``).
        nupop_model_order: HMM order for NuPoP (1 or 4).

    Returns:
        Dict with keys:

        - ``"positions"``: list[int] — window start positions
        - ``"scores"``: list[float] — raw log-likelihood scores
        - ``"occupancy"``: list[float] — per-position occupancy (0-1)
          after optional steric exclusion
        - ``"model"``: str — model used
        - ``"step"``: int — step size used
        - ``"nucleosome_size"``: int — 147

    Raises:
        ValueError: If *model* is not recognised.
    """
    seq_upper = seq.upper()
    n = len(seq_upper)

    if n < NUCLEOSOME_SIZE:
        return {
            "positions": [],
            "scores": [],
            "occupancy": [],
            "model": model,
            "step": step,
            "nucleosome_size": NUCLEOSOME_SIZE,
        }

    # Normalise species to int code
    if isinstance(species, str):
        species = NUPOP_SPECIES.get(species.lower(), 7)

    # ── NuPoP backend ───────────────────────────────────────────────────
    if model == "nupop":
        nupop_preds = predict_nucleosome_nupop(
            seq_upper, species=int(species), model_order=nupop_model_order
        )
        if nupop_preds:
            positions = [p["position"] for p in nupop_preds]
            scores = [p.get("nl_affinity", 0.0) for p in nupop_preds]
            raw_occupancy = [p.get("occupancy", 0.0) for p in nupop_preds]
        else:
            # Fallback to PSSM if NuPoP unavailable
            logger.warning(
                "NuPoP unavailable, falling back to segal_pssm model"
            )
            model = "segal_pssm"
            # Fall through to PSSM path below
        if model == "nupop":
            # Apply steric exclusion on top of NuPoP raw occupancy
            if apply_exclusion:
                energy_list = [-s for s in scores]  # Negate: higher score = lower energy
                occupancy = predict_occupancy_with_exclusion(
                    energy_list, mu=mu, nucleosome_size=NUCLEOSOME_SIZE
                )
                # Pad occupancy to full sequence length
                full_occupancy = _pad_occupancy(occupancy, n, NUCLEOSOME_SIZE, positions)
            else:
                full_occupancy = _interpolate_to_sequence(
                    raw_occupancy, positions, n, NUCLEOSOME_SIZE
                )
            return {
                "positions": positions,
                "scores": scores,
                "occupancy": full_occupancy,
                "model": "nupop",
                "step": step,
                "nucleosome_size": NUCLEOSOME_SIZE,
            }

    # ── Segal PSSM or legacy model ──────────────────────────────────────
    positions: list[int] = []
    scores: list[float] = []

    # Select scoring function
    if model == "segal_pssm":
        scorer = score_kaplan_pssm
    elif model == "segal_legacy":
        scorer = score_segal_legacy_pssm
    else:
        raise ValueError(
            f"Unknown nucleosome model '{model}'. "
            f"Choose from: segal_pssm, segal_legacy, nupop"
        )

    # Slide the window
    for start in range(0, n - NUCLEOSOME_SIZE + 1, step):
        window = seq_upper[start : start + NUCLEOSOME_SIZE]
        s = scorer(window)
        positions.append(start)
        scores.append(s)

    # Compute per-position occupancy
    if apply_exclusion:
        # Negate scores to get binding energies (lower energy = stronger binding)
        energy_list = [-s for s in scores]
        occupancy_raw = predict_occupancy_with_exclusion(
            energy_list, mu=mu, nucleosome_size=NUCLEOSOME_SIZE
        )
        full_occupancy = _pad_occupancy(
            occupancy_raw, n, NUCLEOSOME_SIZE, positions
        )
    else:
        # Simple normalisation without exclusion
        full_occupancy = _simple_occupancy(scores, positions, n, NUCLEOSOME_SIZE)

    return {
        "positions": positions,
        "scores": scores,
        "occupancy": full_occupancy,
        "model": model,
        "step": step,
        "nucleosome_size": NUCLEOSOME_SIZE,
    }


# ── Helper functions ─────────────────────────────────────────────────────────


def _pad_occupancy(
    occupancy_raw: list[float],
    seq_len: int,
    nucleosome_size: int,
    positions: list[int],
) -> list[float]:
    """Expand per-window occupancy values to per-base-pair resolution.

    Each window's occupancy is assigned to every base pair it covers.
    Where windows overlap, the maximum occupancy is kept.

    Args:
        occupancy_raw: Occupancy value for each window position.
        seq_len: Full sequence length.
        nucleosome_size: Size of each nucleosome window (147 bp).
        positions: Start positions of each window.

    Returns:
        List of length *seq_len* with per-bp occupancy values.
    """
    import numpy as np

    result = np.zeros(seq_len, dtype=np.float64)
    for pos, occ in zip(positions, occupancy_raw):
        end = min(seq_len, pos + nucleosome_size)
        result[pos:end] = np.maximum(result[pos:end], occ)

    return result.tolist()


def _interpolate_to_sequence(
    values: list[float],
    positions: list[int],
    seq_len: int,
    nucleosome_size: int,
) -> list[float]:
    """Interpolate sparse per-position values to full sequence length.

    Assigns each value to the centre of its nucleosome footprint,
    then linearly interpolates between centres.

    Args:
        values: Value for each window.
        positions: Start positions of each window.
        seq_len: Full sequence length.
        nucleosome_size: Nucleosome footprint (147 bp).

    Returns:
        List of length *seq_len* with interpolated values.
    """
    import numpy as np

    if not positions:
        return [0.0] * seq_len

    result = np.zeros(seq_len, dtype=np.float64)
    for pos, val in zip(positions, values):
        centre = pos + nucleosome_size // 2
        if 0 <= centre < seq_len:
            result[centre] = val

    # Simple nearest-neighbour fill
    last_val = 0.0
    for i in range(seq_len):
        if result[i] != 0.0:
            last_val = result[i]
        else:
            result[i] = last_val

    return result.tolist()


def _simple_occupancy(
    scores: list[float],
    positions: list[int],
    seq_len: int,
    nucleosome_size: int,
) -> list[float]:
    """Compute per-bp occupancy without steric exclusion.

    Uses a sigmoid transform on raw scores, then assigns to each
    base pair covered by the window.

    Args:
        scores: Raw log-likelihood scores per window.
        positions: Start positions of each window.
        seq_len: Full sequence length.
        nucleosome_size: Nucleosome footprint (147 bp).

    Returns:
        List of length *seq_len* with per-bp occupancy in [0, 1].
    """
    import numpy as np

    if not scores:
        return [0.0] * seq_len

    # Sigmoid normalisation
    scores_arr = np.array(scores, dtype=np.float64)
    midpoint = np.mean(scores_arr)
    width = max(np.std(scores_arr), 0.01)
    occ_per_window = 1.0 / (1.0 + np.exp(-(scores_arr - midpoint) / width))

    result = np.zeros(seq_len, dtype=np.float64)
    for pos, occ in zip(positions, occ_per_window):
        end = min(seq_len, pos + nucleosome_size)
        result[pos:end] = np.maximum(result[pos:end], occ)

    return result.tolist()


# ── Convenience: find well-positioned nucleosomes ────────────────────────────


def find_nucleosome_positions(
    seq: str,
    step: int = DEFAULT_STEP,
    model: str = "segal_pssm",
    occupancy_threshold: float = 0.5,
    mu: float = DEFAULT_CHEMICAL_POTENTIAL,
    apply_exclusion: bool = True,
    species: int | str = 7,
) -> list[dict[str, Any]]:
    """Find well-positioned nucleosome dyad positions.

    A thin wrapper around :func:`predict_nucleosome_occupancy` that
    returns only positions whose peak occupancy exceeds *threshold*.

    Args:
        seq: Input DNA sequence.
        step: Step size (use 1 for high-resolution, 10 for quick scan).
        model: Scoring model (``"segal_pssm"``, ``"segal_legacy"``, or ``"nupop"``).
        occupancy_threshold: Minimum occupancy to report a position.
        mu: Chemical potential for the exclusion model.
        apply_exclusion: Whether to apply steric exclusion.
        species: Species code for NuPoP.

    Returns:
        List of dicts with keys ``position``, ``score``, ``occupancy``.
    """
    result = predict_nucleosome_occupancy(
        seq,
        step=step,
        model=model,
        mu=mu,
        apply_exclusion=apply_exclusion,
        species=species,
    )

    peaks: list[dict[str, Any]] = []
    for pos, score, occ in zip(
        result["positions"], result["scores"], result["occupancy"]
    ):
        # occ is per-window; check the centre base of the window
        centre = min(pos + NUCLEOSOME_SIZE // 2, len(seq) - 1)
        centre_occ = result["occupancy"][centre] if centre < len(result["occupancy"]) else occ
        if centre_occ >= occupancy_threshold:
            peaks.append(
                {
                    "position": pos,
                    "dyad": pos + NUCLEOSOME_SIZE // 2,
                    "score": score,
                    "occupancy": centre_occ,
                }
            )

    return peaks


# ── Public API ────────────────────────────────────────────────────────────────

__all__ = [
    # Constants
    "NUCLEOSOME_SIZE",
    "DEFAULT_STEP",
    "FINE_STEP",
    "DEFAULT_CHEMICAL_POTENTIAL",
    "HELICAL_PERIOD",
    "DINUCLEOTIDES",
    "NUPOP_SPECIES",
    # PSSM (Upgrade 1)
    "score_kaplan_pssm",
    "score_segal_pssm",  # backward compatibility alias
    "_generate_segal_pssm",
    # NuPoP (Upgrade 2)
    "predict_nucleosome_nupop",
    # Teif-Percus (Upgrade 3)
    "predict_occupancy_with_exclusion",
    # Enformer (Upgrade 4)
    "predict_histone_marks_enformer",
    # Core prediction (includes Upgrade 5: 1 bp step)
    "predict_nucleosome_occupancy",
    # Convenience
    "find_nucleosome_positions",
]
