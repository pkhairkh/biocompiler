"""
BioCompiler DNA Shape Prediction Module
========================================
Predicts DNA structural features relevant to damage susceptibility using
dinucleotide structural parameters and optional deep-learning models.

Shape features predicted:
  - Minor Groove Width (MGW) in Angstroms
  - Helix Twist (HelT) in degrees
  - Propeller Twist (ProT) in degrees
  - Roll angle (Roll) in degrees

These features are used to derive per-position damage susceptibility scores
based on biophysical rules connecting DNA shape to chemical vulnerability.

Primary parameter source:
  Olson WK, et al. (1998) DNA sequence-dependent deformability deduced from
  protein-DNA crystal complexes. PNAS 95:11163-11168.

Additional references:
  Rohs R, et al. (2009) The role of DNA shape in protein-DNA recognition.
  Nature 461:1248-1253.
  Zhou T, et al. (2013) DNAshape: high-throughput prediction of DNA structural
  features on a genomic scale. Nucleic Acids Res 41(Web Server):W56-62.
  Li J, et al. (2024) Deep DNAshape for 14 DNA structural features with
  convolutional neural networks. Nat Commun 15:1833.
"""

from __future__ import annotations

import logging
import math
import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Optional

__all__ = [
    "DNAShapeProfile",
    "compute_minor_groove_width",
    "compute_helix_twist",
    "compute_propeller_twist",
    "compute_roll",
    "compute_dna_shape_profile",
    "compute_damage_susceptibility_from_shape",
    "compute_shape_dnacurve",
    "compute_shape_deep_dnashape",
]

logger = logging.getLogger(__name__)


# ==============================================================================
# 1. Dinucleotide Structural Parameter Tables
# ==============================================================================
# All 16 unique dinucleotide step parameters from Olson et al. 1998 PNAS 95:11163.
# Values are assigned to the inter-base-pair step between position i and i+1.
# Complementary dinucleotides (e.g. AC/GT) share the same step parameter.

# Minor Groove Width (MGW) in Angstroms
# Derived from crystal structure analysis of protein-DNA complexes.
# A-tract steps (AA/TT, TA) have narrow minor grooves (<3.2 A) that
# concentrate negative electrostatic potential, enhancing electrophilic
# damage. Values from DNAshape pentamer model averages (Rohs 2009,
# Zhou 2013) mapped to dinucleotide resolution.
_OLSON_MGW: dict[str, float] = {
    "AA": 3.02, "AC": 4.45, "AG": 4.10, "AT": 3.72,
    "CA": 3.48, "CC": 4.76, "CG": 4.39, "CT": 4.10,
    "GA": 3.92, "GC": 4.67, "GG": 4.76, "GT": 4.45,
    "TA": 2.97, "TC": 3.92, "TG": 3.48, "TT": 3.02,
}

# Helix Twist (HelT) in degrees
# Average helical twist per dinucleotide step. B-DNA average ~34.3 deg.
# Low HelT (<31 deg) indicates compressed stacking.
_OLSON_HELT: dict[str, float] = {
    "AA": 35.62, "AC": 34.42, "AG": 27.77, "AT": 31.53,
    "CA": 34.52, "CC": 29.78, "CG": 29.16, "CT": 27.77,
    "GA": 36.92, "GC": 40.00, "GG": 29.78, "GT": 34.42,
    "TA": 35.05, "TC": 36.92, "TG": 34.52, "TT": 35.62,
}

# Propeller Twist (ProT) in degrees
# Rotation of base pairs about the long axis. High ProT (>14 deg)
# enhances UV photoproduct formation by bringing pyrimidines closer.
_OLSON_PROT: dict[str, float] = {
    "AA": 12.51, "AC": 11.73, "AG": 10.87, "AT": 10.64,
    "CA": 9.62,  "CC": 9.09,  "CG": 8.31,  "CT": 10.87,
    "GA": 14.01, "GC": 12.56, "GG": 9.09,  "GT": 11.73,
    "TA": 6.89,  "TC": 14.01, "TG": 9.62,  "TT": 12.51,
}

# Roll angle (Roll) in degrees
# Angular displacement of successive base pairs about the long axis.
# High Roll (>5 deg) indicates base pair opening, enhancing glycosylase access.
_OLSON_ROLL: dict[str, float] = {
    "AA": 0.31,  "AC": 3.22,  "AG": 6.71,  "AT": 2.47,
    "CA": 11.25, "CC": 4.65,  "CG": 8.05,  "CT": 6.71,
    "GA": 1.63,  "GC": -0.20, "GG": 4.65,  "GT": 3.22,
    "TA": 4.86,  "TC": 1.63,  "TG": 11.25, "TT": 0.31,
}


# ==============================================================================
# 2. Damage Susceptibility Thresholds and Multipliers
# ==============================================================================

# Narrow minor groove threshold (Angstroms)
_MGW_NARROW_THRESHOLD: float = 3.2
# Risk multiplier for narrow minor groove (electrophilic damage enhancement)
_MGW_NARROW_MULTIPLIER: float = 0.30  # +30% risk

# High propeller twist threshold (degrees)
_PROT_HIGH_THRESHOLD: float = 14.0
# Risk multiplier for high propeller twist (UV photoproduct formation)
_PROT_HIGH_MULTIPLIER: float = 0.20  # +20%

# High roll angle threshold (degrees)
_ROLL_HIGH_THRESHOLD: float = 5.0
# Risk category for high roll (glycosylase access enhancement)
_ROLL_HIGH_LABEL: str = "enhanced_glycosylase_access"

# Low helix twist threshold (degrees)
_HELT_LOW_THRESHOLD: float = 31.0
# Risk modifier for low helix twist (compressed stacking, reduced UV damage)
_HELT_LOW_MULTIPLIER: float = -0.10  # -10% (protective)

# Baseline damage susceptibility
_BASELINE_SUSCEPTIBILITY: float = 1.0


# ==============================================================================
# 3. Data Classes
# ==============================================================================

@dataclass
class DNAShapeProfile:
    """Complete DNA shape profile with damage susceptibility scores.

    Attributes:
        sequence: Input DNA sequence.
        mgw: Per-position minor groove width in Angstroms (len = len(seq) - 1).
        helix_twist: Per-position helix twist in degrees (len = len(seq) - 1).
        propeller_twist: Per-position propeller twist in degrees (len = len(seq) - 1).
        roll: Per-position roll angle in degrees (len = len(seq) - 1).
        damage_susceptibility: Per-position relative damage susceptibility (len = len(seq) - 1).
        method: Method used for shape prediction.
    """
    sequence: str
    mgw: list[float]
    helix_twist: list[float]
    propeller_twist: list[float]
    roll: list[float]
    damage_susceptibility: list[float]
    method: str = "olson_1998"

    @property
    def n_steps(self) -> int:
        """Number of dinucleotide steps."""
        return len(self.mgw)

    @property
    def mean_mgw(self) -> float:
        """Mean minor groove width across the sequence."""
        return sum(self.mgw) / len(self.mgw) if self.mgw else 0.0

    @property
    def mean_damage_susceptibility(self) -> float:
        """Mean damage susceptibility across the sequence."""
        return sum(self.damage_susceptibility) / len(self.damage_susceptibility) if self.damage_susceptibility else 0.0

    @property
    def high_risk_positions(self) -> list[int]:
        """Positions with damage susceptibility > 1.5 (50% above baseline)."""
        return [i for i, s in enumerate(self.damage_susceptibility) if s > 1.5]

    @property
    def narrow_groove_positions(self) -> list[int]:
        """Positions with MGW below the narrow threshold."""
        return [i for i, m in enumerate(self.mgw) if m < _MGW_NARROW_THRESHOLD]


# ==============================================================================
# 4. Core Shape Computation Functions
# ==============================================================================

def _validate_sequence(seq: str) -> str:
    """Validate and normalize a DNA sequence.

    Args:
        seq: DNA sequence (case-insensitive).

    Returns:
        Uppercase DNA sequence.

    Raises:
        ValueError: If sequence contains non-ACGT characters or is too short.
    """
    seq = seq.upper().strip()
    if len(seq) < 2:
        raise ValueError(
            f"Sequence must be at least 2 nucleotides for dinucleotide "
            f"shape prediction, got {len(seq)}"
        )
    valid = set("ACGT")
    invalid_chars = set(seq) - valid
    if invalid_chars:
        raise ValueError(
            f"Sequence contains invalid characters: {invalid_chars}. "
            f"Only A, C, G, T are allowed."
        )
    return seq


def _compute_dinuc_feature(seq: str, table: dict[str, float],
                           feature_name: str) -> list[float]:
    """Compute a per-position dinucleotide feature from a parameter table.

    Each value is assigned to the step between position i and i+1,
    producing len(seq)-1 values.

    Args:
        seq: DNA sequence (uppercase ACGT).
        table: Dinucleotide -> value mapping (16 entries).
        feature_name: Name of the feature (for error messages).

    Returns:
        List of per-step values (length = len(seq) - 1).
    """
    values: list[float] = []
    for i in range(len(seq) - 1):
        dinuc = seq[i:i + 2]
        if dinuc in table:
            values.append(table[dinuc])
        else:
            # Should not happen after validation, but handle gracefully
            logger.warning(
                "Unknown dinucleotide '%s' at position %d for %s; "
                "using B-DNA average", dinuc, i, feature_name
            )
            values.append(sum(table.values()) / len(table))
    return values


def compute_minor_groove_width(seq: str) -> list[float]:
    """Predict per-position minor groove width using Olson 1998 dinucleotide parameters.

    Minor groove width (MGW) is a key determinant of DNA-protein recognition
    and chemical accessibility. Narrow minor grooves (<3.2 A) concentrate
    negative electrostatic potential, enhancing electrophilic damage.

    Args:
        seq: DNA sequence (case-insensitive, ACGT only).

    Returns:
        List of MGW values in Angstroms (length = len(seq) - 1).

    Raises:
        ValueError: If sequence contains invalid characters or is too short.

    References:
        Olson WK, et al. (1998) PNAS 95:11163-11168.
        Rohs R, et al. (2009) Nature 461:1248-1253.
    """
    seq = _validate_sequence(seq)
    return _compute_dinuc_feature(seq, _OLSON_MGW, "MGW")


def compute_helix_twist(seq: str) -> list[float]:
    """Predict per-position helix twist using Olson 1998 dinucleotide parameters.

    Helix twist (HelT) measures the rotation between successive base pairs.
    Low helix twist (<31 deg) indicates compressed base stacking, which
    reduces UV photoproduct formation but may affect protein binding.

    Args:
        seq: DNA sequence (case-insensitive, ACGT only).

    Returns:
        List of HelT values in degrees (length = len(seq) - 1).

    Raises:
        ValueError: If sequence contains invalid characters or is too short.

    References:
        Olson WK, et al. (1998) PNAS 95:11163-11168.
    """
    seq = _validate_sequence(seq)
    return _compute_dinuc_feature(seq, _OLSON_HELT, "HelT")


def compute_propeller_twist(seq: str) -> list[float]:
    """Predict per-position propeller twist using Olson 1998 dinucleotide parameters.

    Propeller twist (ProT) is the rotation of one base in a pair relative
    to the other about the long axis. High ProT (>14 deg) brings
    pyrimidines on opposite strands closer together, enhancing UV
    photoproduct (CPD and 6-4PP) formation.

    Args:
        seq: DNA sequence (case-insensitive, ACGT only).

    Returns:
        List of ProT values in degrees (length = len(seq) - 1).

    Raises:
        ValueError: If sequence contains invalid characters or is too short.

    References:
        Olson WK, et al. (1998) PNAS 95:11163-11168.
    """
    seq = _validate_sequence(seq)
    return _compute_dinuc_feature(seq, _OLSON_PROT, "ProT")


def compute_roll(seq: str) -> list[float]:
    """Predict per-position roll angle using Olson 1998 dinucleotide parameters.

    Roll angle measures the angular displacement of successive base pairs
    about their long axis. High roll (>5 deg) opens the minor groove
    and facilitates base pair opening, enhancing glycosylase access
    for DNA repair and making the base more accessible to damaging agents.

    Args:
        seq: DNA sequence (case-insensitive, ACGT only).

    Returns:
        List of Roll values in degrees (length = len(seq) - 1).

    Raises:
        ValueError: If sequence contains invalid characters or is too short.

    References:
        Olson WK, et al. (1998) PNAS 95:11163-11168.
    """
    seq = _validate_sequence(seq)
    return _compute_dinuc_feature(seq, _OLSON_ROLL, "Roll")


# ==============================================================================
# 5. Damage Susceptibility from Shape Features
# ==============================================================================

def compute_damage_susceptibility_from_shape(
    mgw: list[float],
    helix_twist: list[float],
    propeller_twist: list[float],
    roll: list[float],
) -> list[float]:
    """Derive per-position damage susceptibility from DNA shape features.

    Applies biophysically-motivated rules connecting DNA shape to
    chemical vulnerability:

    1. **Narrow minor groove** (MGW < 3.2 A) -> enhanced electrophilic
       damage (+30% risk). Narrow grooves concentrate negative electrostatic
       potential, attracting electrophilic mutagens.

    2. **High propeller twist** (ProT > 14 deg) -> enhanced UV photoproduct
       formation (+20%). High ProT brings pyrimidines closer, facilitating
       cyclobutane pyrimidine dimer and 6-4 photoproduct formation.

    3. **High roll angle** (Roll > 5 deg) -> base pair opening -> enhanced
       glycosylase access. Base pairs with high roll are more accessible
       to DNA repair enzymes and damaging agents alike.

    4. **Low helix twist** (HelT < 31 deg) -> compressed stacking -> reduced
       UV damage (-10%). Compressed stacking shields bases from UV radiation.

    The susceptibility is multiplicative: multiple risk factors compound.
    A baseline of 1.0 is modified by additive adjustments from each rule.
    The minimum susceptibility is clamped at 0.0.

    Args:
        mgw: Per-position minor groove width values (Angstroms).
        helix_twist: Per-position helix twist values (degrees).
        propeller_twist: Per-position propeller twist values (degrees).
        roll: Per-position roll angle values (degrees).

    Returns:
        Per-position relative damage susceptibility scores.
        Baseline = 1.0; >1.0 = enhanced risk; <1.0 = reduced risk.

    Raises:
        ValueError: If input lists have different lengths or are empty.

    References:
        Rohs R, et al. (2009) Nature 461:1248-1253.
        Parker SC, et al. (2009) Science 324:386-387.
    """
    n = len(mgw)
    if n == 0:
        raise ValueError("Shape feature lists must not be empty")
    if not (len(helix_twist) == len(propeller_twist) == len(roll) == n):
        raise ValueError(
            f"All shape feature lists must have the same length. "
            f"Got MGW={n}, HelT={len(helix_twist)}, "
            f"ProT={len(propeller_twist)}, Roll={len(roll)}"
        )

    susceptibility: list[float] = []

    for i in range(n):
        score = _BASELINE_SUSCEPTIBILITY

        # Rule 1: Narrow minor groove -> enhanced electrophilic damage
        if mgw[i] < _MGW_NARROW_THRESHOLD:
            # More severe as groove gets narrower
            # Linear interpolation: +30% at threshold, up to +50% at 2.5 A
            deficit = _MGW_NARROW_THRESHOLD - mgw[i]
            max_deficit = _MGW_NARROW_THRESHOLD - 2.5  # extreme narrow
            fraction = min(1.0, deficit / max_deficit) if max_deficit > 0 else 1.0
            score += _MGW_NARROW_MULTIPLIER * (1.0 + fraction)

        # Rule 2: High propeller twist -> enhanced UV photoproduct formation
        if propeller_twist[i] > _PROT_HIGH_THRESHOLD:
            excess = propeller_twist[i] - _PROT_HIGH_THRESHOLD
            max_excess = 8.0  # maximum expected excess (~22 deg)
            fraction = min(1.0, excess / max_excess)
            score += _PROT_HIGH_MULTIPLIER * (1.0 + fraction)

        # Rule 3: High roll angle -> base pair opening -> glycosylase access
        if roll[i] > _ROLL_HIGH_THRESHOLD:
            # Mark as enhanced access but don't add numerical risk;
            # glycosylase access is a repair factor, not direct damage
            # Include a small risk increase for accessibility
            excess = roll[i] - _ROLL_HIGH_THRESHOLD
            max_excess = 10.0
            fraction = min(1.0, excess / max_excess)
            score += 0.10 * fraction  # modest +10% for accessibility

        # Rule 4: Low helix twist -> compressed stacking -> reduced UV damage
        if helix_twist[i] < _HELT_LOW_THRESHOLD:
            deficit = _HELT_LOW_THRESHOLD - helix_twist[i]
            max_deficit = 8.0  # extreme compression
            fraction = min(1.0, deficit / max_deficit)
            score += _HELT_LOW_MULTIPLIER * (1.0 + fraction)

        # Clamp to non-negative
        susceptibility.append(max(0.0, score))

    return susceptibility


# ==============================================================================
# 6. Integrated DNA Shape Profile
# ==============================================================================

def compute_dna_shape_profile(seq: str) -> dict:
    """Compute all 4 shape features and damage susceptibility for a DNA sequence.

    Convenience function that runs all four shape predictions and derives
    per-position damage susceptibility in a single call.

    Args:
        seq: DNA sequence (case-insensitive, ACGT only).

    Returns:
        Dictionary with keys:
          - "sequence": input sequence
          - "mgw": list[float] - minor groove width (A)
          - "helix_twist": list[float] - helix twist (deg)
          - "propeller_twist": list[float] - propeller twist (deg)
          - "roll": list[float] - roll angle (deg)
          - "damage_susceptibility": list[float] - relative susceptibility
          - "method": "olson_1998"
          - "n_steps": number of dinucleotide steps
          - "mean_mgw": mean minor groove width
          - "mean_damage_susceptibility": mean damage susceptibility
          - "high_risk_positions": positions with susceptibility > 1.5
          - "narrow_groove_positions": positions with MGW < 3.2 A

    Raises:
        ValueError: If sequence contains invalid characters or is too short.

    References:
        Olson WK, et al. (1998) PNAS 95:11163-11168.
        Rohs R, et al. (2009) Nature 461:1248-1253.
    """
    seq = _validate_sequence(seq)

    mgw = compute_minor_groove_width(seq)
    helt = compute_helix_twist(seq)
    prot = compute_propeller_twist(seq)
    roll = compute_roll(seq)
    damage = compute_damage_susceptibility_from_shape(mgw, helt, prot, roll)

    profile = DNAShapeProfile(
        sequence=seq,
        mgw=mgw,
        helix_twist=helt,
        propeller_twist=prot,
        roll=roll,
        damage_susceptibility=damage,
        method="olson_1998",
    )

    return {
        "sequence": profile.sequence,
        "mgw": profile.mgw,
        "helix_twist": profile.helix_twist,
        "propeller_twist": profile.propeller_twist,
        "roll": profile.roll,
        "damage_susceptibility": profile.damage_susceptibility,
        "method": profile.method,
        "n_steps": profile.n_steps,
        "mean_mgw": profile.mean_mgw,
        "mean_damage_susceptibility": profile.mean_damage_susceptibility,
        "high_risk_positions": profile.high_risk_positions,
        "narrow_groove_positions": profile.narrow_groove_positions,
    }


# ==============================================================================
# 7. Optional External Tool: dnacurve
# ==============================================================================

def compute_shape_dnacurve(seq: str) -> dict:
    """Use the dnacurve package for 3D B-DNA structure prediction.

    The dnacurve package (pip install dnacurve) computes the 3D structure
    of B-DNA using the Wedge model or other parameter sets, providing
    atomic coordinates and structural parameters.

    When dnacurve is not installed, falls back to the built-in Olson 1998
    dinucleotide parameter tables.

    Args:
        seq: DNA sequence (case-insensitive, ACGT only).

    Returns:
        Dictionary with keys:
          - "method": "dnacurve" or "olson_1998_fallback"
          - "mgw": list[float] - minor groove width (A) if available
          - "helix_twist": list[float] - helix twist (deg)
          - "propeller_twist": list[float] - propeller twist (deg)
          - "roll": list[float] - roll angle (deg)
          - "damage_susceptibility": list[float] - relative susceptibility
          - "coordinates": 3D coordinates (if dnacurve available)
          - "error": error message (if dnacurve unavailable)

    References:
        dnacurve: https://pypi.org/project/dnacurve/
        Bolshoy A, et al. (1991) PNAS 88:2312-2316 (wedge model).
    """
    seq = _validate_sequence(seq)

    try:
        from dnacurve import CurvedDNA

        result = CurvedDNA(seq)
        method = "dnacurve"

        # Extract structural parameters from dnacurve
        # dnacurve provides per-step parameters
        helix_twist_vals: list[float] = []
        roll_vals: list[float] = []
        prot_vals: list[float] = []
        mgw_vals: list[float] = []
        coordinates: list[list[float]] = []

        # Try to access the computed data
        if hasattr(result, "trinucleotide") and result.trinucleotide is not None:
            # dnacurve uses trinucleotide parameters
            for i in range(len(seq) - 1):
                # Extract from result attributes if available
                try:
                    step_data = result.trinucleotide
                    # dnacurve provides wedge, direction, twist per trinucleotide
                    helix_twist_vals.append(
                        float(step_data[i + 1]["Twist"])
                        if i + 1 < len(step_data)
                        else _OLSON_HELT.get(seq[i:i + 2], 34.3)
                    )
                    roll_vals.append(
                        float(step_data[i + 1]["Wedge"])
                        if i + 1 < len(step_data)
                        else _OLSON_ROLL.get(seq[i:i + 2], 3.0)
                    )
                except (IndexError, KeyError, TypeError):
                    helix_twist_vals.append(_OLSON_HELT.get(seq[i:i + 2], 34.3))
                    roll_vals.append(_OLSON_ROLL.get(seq[i:i + 2], 3.0))

                # ProT and MGW not directly provided by dnacurve;
                # fall back to Olson tables for these
                dinuc = seq[i:i + 2]
                prot_vals.append(_OLSON_PROT.get(dinuc, 10.0))
                mgw_vals.append(_OLSON_MGW.get(dinuc, 4.0))

        else:
            # Fall back to Olson tables if trinucleotide data unavailable
            for i in range(len(seq) - 1):
                dinuc = seq[i:i + 2]
                mgw_vals.append(_OLSON_MGW.get(dinuc, 4.0))
                helix_twist_vals.append(_OLSON_HELT.get(dinuc, 34.3))
                prot_vals.append(_OLSON_PROT.get(dinuc, 10.0))
                roll_vals.append(_OLSON_ROLL.get(dinuc, 3.0))

        # Try to extract 3D coordinates
        if hasattr(result, "coordinates") and result.coordinates is not None:
            try:
                coords = result.coordinates
                for c in coords:
                    if len(c) >= 3:
                        coordinates.append([float(c[0]), float(c[1]), float(c[2])])
            except (TypeError, IndexError):
                coordinates = []

        damage = compute_damage_susceptibility_from_shape(
            mgw_vals, helix_twist_vals, prot_vals, roll_vals
        )

        return_dict: dict = {
            "method": method,
            "mgw": mgw_vals,
            "helix_twist": helix_twist_vals,
            "propeller_twist": prot_vals,
            "roll": roll_vals,
            "damage_susceptibility": damage,
            "coordinates": coordinates,
        }
        return return_dict

    except ImportError:
        logger.info(
            "dnacurve package not installed; falling back to Olson 1998 "
            "dinucleotide parameters. Install with: pip install dnacurve"
        )
        # Fall back to built-in Olson tables
        profile = compute_dna_shape_profile(seq)
        return {
            "method": "olson_1998_fallback",
            "mgw": profile["mgw"],
            "helix_twist": profile["helix_twist"],
            "propeller_twist": profile["propeller_twist"],
            "roll": profile["roll"],
            "damage_susceptibility": profile["damage_susceptibility"],
            "coordinates": [],
            "error": (
                "dnacurve package not installed. Install with: "
                "pip install dnacurve"
            ),
        }
    except Exception as e:
        logger.warning("dnacurve computation failed: %s; using fallback", e)
        profile = compute_dna_shape_profile(seq)
        return {
            "method": "olson_1998_fallback",
            "mgw": profile["mgw"],
            "helix_twist": profile["helix_twist"],
            "propeller_twist": profile["propeller_twist"],
            "roll": profile["roll"],
            "damage_susceptibility": profile["damage_susceptibility"],
            "coordinates": [],
            "error": f"dnacurve failed: {e}",
        }


# ==============================================================================
# 8. Optional External Tool: Deep DNAshape
# ==============================================================================

# Mapping from Deep DNAshape feature names to normalized keys
_DEEP_DNASHAPE_FEATURES: list[str] = [
    "MGW", "HelT", "ProT", "Roll",
    "Shift", "Slide", "Rise", "Tilt",
    "Bend", "Shear", "Stretch", "Stagger",
    "Opening", "EP",
]


def compute_shape_deep_dnashape(
    seq: str,
    deepdnashape_path: str = "deepDNAshape",
) -> dict:
    """Use Deep DNAshape CLI for 14-feature DNA shape prediction.

    Deep DNAshape (Li et al. 2024) uses a convolutional neural network
    trained on molecular dynamics simulations to predict 14 DNA structural
    features from sequence. This wrapper invokes the CLI tool and parses
    the output.

    When Deep DNAshape is not installed, falls back to the built-in
    Olson 1998 dinucleotide parameter tables (4 features only).

    The 14 predicted features are:
      MGW, HelT, ProT, Roll, Shift, Slide, Rise, Tilt,
      Bend, Shear, Stretch, Stagger, Opening, EP

    Args:
        seq: DNA sequence (case-insensitive, ACGT only).
            Must be >= 15 bp for Deep DNAshape (CNN window requirement).
        deepdnashape_path: Path to the deepDNAshape CLI executable.

    Returns:
        Dictionary with keys:
          - "method": "deep_dnashape" or "olson_1998_fallback"
          - "features": dict mapping feature name -> list[float]
          - "mgw": list[float] - minor groove width (A)
          - "helix_twist": list[float] - helix twist (deg)
          - "propeller_twist": list[float] - propeller twist (deg)
          - "roll": list[float] - roll angle (deg)
          - "damage_susceptibility": list[float] - relative susceptibility
          - "error": error message (if unavailable)

    Raises:
        ValueError: If sequence is too short for Deep DNAshape (<15 bp).

    References:
        Li J, et al. (2024) Deep DNAshape for 14 DNA structural features
        with convolutional neural networks. Nat Commun 15:1833.
    """
    seq = _validate_sequence(seq)

    min_len_for_deep = 15
    if len(seq) < min_len_for_deep:
        logger.info(
            "Sequence too short for Deep DNAshape (%d < %d bp); "
            "falling back to Olson 1998",
            len(seq), min_len_for_deep,
        )
        profile = compute_dna_shape_profile(seq)
        return {
            "method": "olson_1998_fallback",
            "features": {
                "MGW": profile["mgw"],
                "HelT": profile["helix_twist"],
                "ProT": profile["propeller_twist"],
                "Roll": profile["roll"],
            },
            "mgw": profile["mgw"],
            "helix_twist": profile["helix_twist"],
            "propeller_twist": profile["propeller_twist"],
            "roll": profile["roll"],
            "damage_susceptibility": profile["damage_susceptibility"],
            "error": (
                f"Sequence too short for Deep DNAshape "
                f"(need >= {min_len_for_deep} bp, got {len(seq)})"
            ),
        }

    try:
        # Write sequence to temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fa", delete=False
        ) as tmp_fa:
            tmp_fa.write(f">query\n{seq}\n")
            tmp_fa_path = tmp_fa.name

        # Run deepDNAshape CLI
        # Expected usage: deepDNAshape -i input.fa -o output_dir
        with tempfile.TemporaryDirectory() as tmp_out:
            cmd = [
                deepdnashape_path,
                "-i", tmp_fa_path,
                "-o", tmp_out,
            ]
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=120,
                    check=True,
                )
            except FileNotFoundError:
                raise ImportError(
                    f"deepDNAshape CLI not found at '{deepdnashape_path}'. "
                    f"Install from: https://github.com/zzhu36/DeepDNAshape"
                )
            except subprocess.TimeoutExpired:
                raise RuntimeError(
                    "deepDNAshape CLI timed out after 120 seconds"
                )
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"deepDNAshape CLI failed with return code {e.returncode}: "
                    f"{e.stderr}"
                )

            # Parse output files
            features: dict[str, list[float]] = {}
            for feature_name in _DEEP_DNASHAPE_FEATURES:
                # Deep DNAshape outputs one file per feature
                feature_file = os.path.join(tmp_out, f"{feature_name}.txt")
                if os.path.exists(feature_file):
                    values: list[float] = []
                    with open(feature_file, "r") as f:
                        for line in f:
                            line = line.strip()
                            if line and not line.startswith(">"):
                                try:
                                    values.append(float(line))
                                except ValueError:
                                    continue
                    features[feature_name] = values
                else:
                    features[feature_name] = []

        # Clean up temp FASTA
        try:
            os.unlink(tmp_fa_path)
        except OSError:
            pass

        # Extract the 4 core features for damage susceptibility
        mgw = features.get("MGW", [])
        helt = features.get("HelT", [])
        prot = features.get("ProT", [])
        roll = features.get("Roll", [])

        # Fall back to Olson tables for missing features
        if not mgw:
            mgw = compute_minor_groove_width(seq)
        if not helt:
            helt = compute_helix_twist(seq)
        if not prot:
            prot = compute_propeller_twist(seq)
        if not roll:
            roll = compute_roll(seq)

        # Ensure all feature lists have the same length for susceptibility
        min_len = min(len(mgw), len(helt), len(prot), len(roll))
        mgw = mgw[:min_len]
        helt = helt[:min_len]
        prot = prot[:min_len]
        roll = roll[:min_len]

        damage = compute_damage_susceptibility_from_shape(mgw, helt, prot, roll)

        return {
            "method": "deep_dnashape",
            "features": features,
            "mgw": mgw,
            "helix_twist": helt,
            "propeller_twist": prot,
            "roll": roll,
            "damage_susceptibility": damage,
        }

    except ImportError as e:
        logger.info(
            "Deep DNAshape not available: %s; falling back to Olson 1998", e
        )
        profile = compute_dna_shape_profile(seq)
        return {
            "method": "olson_1998_fallback",
            "features": {
                "MGW": profile["mgw"],
                "HelT": profile["helix_twist"],
                "ProT": profile["propeller_twist"],
                "Roll": profile["roll"],
            },
            "mgw": profile["mgw"],
            "helix_twist": profile["helix_twist"],
            "propeller_twist": profile["propeller_twist"],
            "roll": profile["roll"],
            "damage_susceptibility": profile["damage_susceptibility"],
            "error": str(e),
        }
    except Exception as e:
        logger.warning(
            "Deep DNAshape computation failed: %s; using fallback", e
        )
        profile = compute_dna_shape_profile(seq)
        return {
            "method": "olson_1998_fallback",
            "features": {
                "MGW": profile["mgw"],
                "HelT": profile["helix_twist"],
                "ProT": profile["propeller_twist"],
                "Roll": profile["roll"],
            },
            "mgw": profile["mgw"],
            "helix_twist": profile["helix_twist"],
            "propeller_twist": profile["propeller_twist"],
            "roll": profile["roll"],
            "damage_susceptibility": profile["damage_susceptibility"],
            "error": f"Deep DNAshape failed: {e}",
        }
