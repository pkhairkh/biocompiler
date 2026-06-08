"""
Thermal stability and melting temperature (Tm) calculations for DNA oligonucleotides.

This module provides state-of-the-art Tm prediction using:

- **SantaLucia & Hicks 2004** nearest-neighbor parameters (DNA_NN4) — the
  community-standard thermodynamic table that supersedes the 1998 set with
  corrected initiation parameters.
- **Self-complementary sequence detection** with the correct symmetry-corrected
  Tm formula (ΔS penalty + ln(C_T) instead of ln(C_T/4)).
- **Owczarzy et al. 2008** magnesium correction with a unified mixed-regime
  (0.22 ≤ R < 6.0) implementation using Eq. 16 with monovalent-adjusted
  coefficients.
- **Ka-based dNTP chelation** using the quadratic solution for Mg-dNTP
  binding equilibrium (K_a = 3×10⁴ M⁻¹), replacing the simple 1:1 model.
- **Proper K⁺/Tris⁺ aggregation** per Owczarzy 2008 (K⁺ and Na⁺ are
  interchangeable; Tris⁺ counts fully toward the monovalent pool).

References
----------
SantaLucia J Jr, Hicks D.  The thermodynamics of DNA structural motifs.
Annu Rev Biophys Biomol Struct. 2004;33:415-40.

Owczarzy R, Moreira BG, You Y, et al.  Predicting stability of DNA
duplexes in solutions containing magnesium and monovalent cations.
Biochemistry. 2008;47(19):5336-53.

Owczarzy R, You Y, Moreira BG, et al.  Effects of sodium ions on DNA
duplex oligonucleotides.  Biophys Chem. 2004;111(3):197-213.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

__all__ = [
    # Core Tm functions
    "calculate_tm",
    "calculate_tm_nearest_neighbor",
    "calculate_tm_wallace",
    # Self-complementarity
    "is_self_complementary",
    # Owczarzy salt corrections
    "owczarzy_2004_na_correction",
    "owczarzy_2008_mg_correction",
    # Free Mg calculation
    "compute_free_mg",
    # Thermodynamic computation
    "compute_dh_ds",
    "compute_tm_from_thermo",
    # Data classes
    "ThermoResult",
    "TmResult",
    # Constants
    "R_CAL",
    "DNA_NN4",
    "DNTP_MG_KA",
]

# ═══════════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════════

R_CAL = 1.987  # Gas constant in cal/(mol·K)

# dNTP-Mg binding constant (M⁻¹) — used for Ka-based chelation model.
# From Ahsen et al. (2001) and von Ahsen et al. (2001), Clin Chem.
DNTP_MG_KA = 3.0e4

# Complementarity map
_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "G": "C", "C": "G"}

# ═══════════════════════════════════════════════════════════════════════════
# Upgrade 1: SantaLucia & Hicks 2004 Nearest-Neighbor Parameters (DNA_NN4)
# ═══════════════════════════════════════════════════════════════════════════
#
# The 10 Watson-Crick NN pairs are IDENTICAL between SantaLucia 1998 (DNA_NN3)
# and SantaLucia & Hicks 2004 (DNA_NN4).  Only the initiation parameters
# changed:
#
#   Parameter        DNA_NN3 (1998)      DNA_NN4 (2004)
#   ─────────────    ──────────────      ──────────────
#   Initiation       (0.0,  0.0)         (0.2, -5.7)    ← general
#   Init A/T ΔH      2.3 kcal/mol        2.3 kcal/mol  (unchanged)
#   Init A/T ΔS      4.1 cal/(mol·K)     6.9 cal/(mol·K)  ← shifted
#   Init G/C ΔH      0.1 kcal/mol        0.0 kcal/mol  ← zeroed
#   Init G/C ΔS     -2.8 cal/(mol·K)     0.0 cal/(mol·K)  ← zeroed
#
# Reference: SantaLucia J Jr & Hicks D, Annu Rev Biophys Biomol Struct
#            2004;33:415-40, Table 2.

DNA_NN4: dict[str, tuple[float, float]] = {
    # ── Watson-Crick nearest-neighbor pairs (ΔH kcal/mol, ΔS cal/(mol·K)) ──
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

# ── Initiation parameters — SantaLucia & Hicks 2004 (DNA_NN4) ─────────
# General initiation: (0.2 kcal/mol, -5.7 cal/(mol·K))
_INIT_GENERAL_DH = 0.2    # kcal/mol
_INIT_GENERAL_DS = -5.7   # cal/(mol·K)

# Per-terminal-A/T initiation: (2.3 kcal/mol, 6.9 cal/(mol·K))
_INIT_AT_DH = 2.3   # kcal/mol
_INIT_AT_DS = 6.9    # cal/(mol·K)  — was 4.1 in 1998

# Per-terminal-G/C initiation: (0, 0) — was (0.1, -2.8) in 1998
_INIT_GC_DH = 0.0   # kcal/mol
_INIT_GC_DS = 0.0    # cal/(mol·K)


# ═══════════════════════════════════════════════════════════════════════════
# Upgrade 2: Self-complementary sequence detection
# ═══════════════════════════════════════════════════════════════════════════

def is_self_complementary(seq: str) -> bool:
    """Check if a DNA sequence is self-complementary.

    A sequence is self-complementary if it equals its own reverse complement.
    Examples: ATCGAT, AATATT

    Args:
        seq: DNA sequence

    Returns:
        True if self-complementary
    """
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    rev_comp = "".join(comp.get(b, "N") for b in reversed(seq.upper()))
    return seq.upper() == rev_comp


# ═══════════════════════════════════════════════════════════════════════════
# Data classes
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class ThermoResult:
    """Result of nearest-neighbor thermodynamic computation.

    Attributes:
        dH: Total enthalpy in kcal/mol.
        dS: Total entropy in cal/(mol·K).
        Nbp: Number of base pairs in the duplex.
        f_GC: GC fraction of the sequence.
        is_self_comp: Whether the sequence is self-complementary.
    """
    dH: float
    dS: float
    Nbp: int
    f_GC: float
    is_self_comp: bool


@dataclass
class TmResult:
    """Result of a melting temperature calculation.

    Attributes:
        Tm: Melting temperature in °C.
        Tm_1M: Melting temperature at 1 M NaCl in °C (before salt correction).
        Tm_1M_kelvin: Melting temperature at 1 M NaCl in Kelvin.
        dH: Total enthalpy in kcal/mol.
        dS: Total entropy in cal/(mol·K).
        Nbp: Number of base pairs.
        f_GC: GC fraction.
        is_self_comp: Whether the sequence is self-complementary.
        salt_correction: Salt correction applied in °C (Owczarzy method).
        method: Description of the Tm method used.
    """
    Tm: float
    Tm_1M: float
    Tm_1M_kelvin: float
    dH: float
    dS: float
    Nbp: int
    f_GC: float
    is_self_comp: bool
    salt_correction: float
    method: str


# ═══════════════════════════════════════════════════════════════════════════
# Core thermodynamic computation
# ═══════════════════════════════════════════════════════════════════════════

def compute_dh_ds(seq: str) -> ThermoResult:
    """Compute ΔH and ΔS for a DNA duplex using SantaLucia 2004 NN parameters.

    Uses DNA_NN4 initiation parameters:
    - General initiation: (0.2, -5.7) kcal/mol, cal/(mol·K)
    - Per-terminal A/T:  (2.3, 6.9) kcal/mol, cal/(mol·K)
    - Per-terminal G/C:  (0.0, 0.0) kcal/mol, cal/(mol·K)

    Args:
        seq: DNA sequence (5'→3', DNA bases only).

    Returns:
        ThermoResult with dH, dS, Nbp, f_GC, is_self_comp.
    """
    seq = seq.upper().replace("U", "T")
    Nbp = len(seq)
    if Nbp == 0:
        return ThermoResult(dH=0.0, dS=0.0, Nbp=0, f_GC=0.0,
                            is_self_comp=False)

    # GC fraction
    gc_count = seq.count("G") + seq.count("C")
    f_GC = gc_count / Nbp if Nbp > 0 else 0.0

    # Check self-complementarity
    self_comp = is_self_complementary(seq)

    # ── Initiation parameters (SantaLucia 2004) ────────────────────────
    # General initiation (applies to all sequences)
    dH = _INIT_GENERAL_DH   # kcal/mol
    dS = _INIT_GENERAL_DS   # cal/(mol·K)

    # Terminal base initiation
    first_base = seq[0]
    last_base = seq[-1]

    for terminal in (first_base, last_base):
        if terminal in ("G", "C"):
            dH += _INIT_GC_DH
            dS += _INIT_GC_DS
        else:  # A or T
            dH += _INIT_AT_DH
            dS += _INIT_AT_DS

    # ── Sum nearest-neighbor parameters ────────────────────────────────
    for i in range(Nbp - 1):
        dinuc = seq[i:i + 2]
        params = DNA_NN4.get(dinuc)
        if params is not None:
            dH += params[0]
            dS += params[1]

    return ThermoResult(dH=dH, dS=dS, Nbp=Nbp, f_GC=f_GC,
                        is_self_comp=self_comp)


def compute_tm_from_thermo(
    dH: float,
    dS: float,
    ct: float,
    is_self_comp: bool = False,
) -> float:
    """Compute Tm in °C from thermodynamic parameters.

    Uses the correct formula based on self-complementarity:

    - **Non-self-complementary**:
      Tm = ΔH / (ΔS + R × ln(C_T / 4)) − 273.15

    - **Self-complementary**:
      Tm = ΔH / (ΔS − 1.4 + R × ln(C_T)) − 273.15
      (includes the −1.4 cal/(mol·K) symmetry correction)

    Args:
        dH: Total enthalpy in kcal/mol.
        dS: Total entropy in cal/(mol·K).
        ct: Total strand concentration in M (for non-self-comp, this is
            the total concentration of both strands; for self-comp, it is
            the total concentration of the single self-complementary strand).
        is_self_comp: Whether the sequence is self-complementary.

    Returns:
        Melting temperature in °C.
    """
    dH_cal = dH * 1000.0  # Convert kcal to cal

    if is_self_comp:
        # Symmetry correction: ΔS += −1.4 cal/(mol·K)
        dS_corrected = dS - 1.4
        # Self-complementary: use ln(C_T) not ln(C_T/4)
        denominator = dS_corrected + R_CAL * math.log(ct) if ct > 0 else dS_corrected
    else:
        # Non-self-complementary: use ln(C_T/4)
        if ct / 4.0 > 0:
            denominator = dS + R_CAL * math.log(ct / 4.0)
        else:
            denominator = dS

    if abs(denominator) < 1e-10:
        return 0.0

    return (dH_cal / denominator) - 273.15


# ═══════════════════════════════════════════════════════════════════════════
# Upgrade 4: Ka-based dNTP chelation
# ═══════════════════════════════════════════════════════════════════════════

def compute_free_mg(
    total_Mg: float,
    total_dNTP: float,
) -> float:
    """Compute free Mg²⁺ concentration using Ka-based chelation model.

    Solves the quadratic equilibrium equation for Mg²⁺ binding to dNTP:

        K_a = [Mg·dNTP] / ([Mg_free] × [dNTP_free])

    The quadratic is derived from mass balance:
        [Mg_total]    = [Mg_free] + [Mg·dNTP]
        [dNTP_total]  = [dNTP_free] + [Mg·dNTP]

    This replaces the simple 1:1 chelation model (free_Mg = total_Mg − dNTP)
    with the full equilibrium solution, which is more accurate when dNTP
    concentrations are similar to Mg²⁺ concentrations.

    Args:
        total_Mg: Total Mg²⁺ concentration in M.
        total_dNTP: Total dNTP concentration in M (sum of all four dNTPs).

    Returns:
        Free Mg²⁺ concentration in M.
    """
    if total_dNTP <= 0:
        return total_Mg
    if total_Mg <= 0:
        return 0.0

    Ka = DNTP_MG_KA  # 3.0e4 M⁻¹

    # Solve the quadratic equilibrium for [Mg_free]:
    #   Ka·[Mg_free]² + (Ka·dNTP − Ka·Mg + 1)·[Mg_free] − Mg = 0
    # Solution (positive root):
    #   [Mg_free] = (−(Ka·dNTP − Ka·Mg + 1)
    #               + sqrt((Ka·dNTP − Ka·Mg + 1)² + 4·Ka·Mg)) / (2·Ka)
    inner = Ka * total_dNTP - Ka * total_Mg + 1
    discriminant = inner ** 2 + 4 * Ka * total_Mg

    if discriminant >= 0:
        free_Mg = (-inner + math.sqrt(discriminant)) / (2 * Ka)
        # Clamp to physical range
        free_Mg = max(0.0, min(free_Mg, total_Mg))
    else:
        # Fallback to simple 1:1 model if numerical issues arise
        free_Mg = max(0.0, total_Mg - total_dNTP)

    return free_Mg


# ═══════════════════════════════════════════════════════════════════════════
# Upgrade 5: Proper K⁺/Tris⁺ aggregation
# ═══════════════════════════════════════════════════════════════════════════

def compute_monovalent_conc(
    sodium: float = 0.0,
    potassium: float = 0.0,
    tris: float = 0.0,
) -> float:
    """Compute total monovalent cation concentration.

    Per Owczarzy 2008:
    - Na⁺ and K⁺ are interchangeable (Fig 2 of Owczarzy et al. 2008)
    - Tris⁺ counts fully as a monovalent cation

    Args:
        sodium: Na⁺ concentration in M.
        potassium: K⁺ concentration in M.
        tris: Tris⁺ concentration in M.

    Returns:
        Total monovalent cation concentration in M.
    """
    return sodium + potassium + tris


# ═══════════════════════════════════════════════════════════════════════════
# Owczarzy 2004 Na⁺ correction
# ═══════════════════════════════════════════════════════════════════════════

def owczarzy_2004_na_correction(
    Tm_1M_kelvin: float,
    Nbp: int,
    f_GC: float,
    sodium: float,
) -> float:
    """Apply Owczarzy et al. 2004 sodium correction to Tm.

    Corrects Tm from 1 M NaCl conditions to the specified monovalent
    cation concentration using Owczarzy 2004 Eq 3:

        1/Tm = 1/Tm(1M) + (4.29×f_GC − 3.95)×10⁻⁵ × ln[Mon⁺]
               + 9.87×10⁻⁶ × (ln[Mon⁺])²

    This formulation correctly returns Tm(1M) when [Mon⁺] = 1 M
    (ln[1] = 0 → correction = 0).

    Reference: Owczarzy R, You Y, Moreira BG, et al. Effects of sodium
    ions on DNA duplex oligonucleotides. Biophys Chem. 2004;111(3):197-213.

    Args:
        Tm_1M_kelvin: Tm at 1 M NaCl in Kelvin.
        Nbp: Number of base pairs in the duplex (unused in Eq 3, kept
            for API compatibility).
        f_GC: GC fraction of the sequence.
        sodium: Monovalent cation concentration in M.

    Returns:
        Corrected Tm in Kelvin.
    """
    if sodium <= 0:
        return Tm_1M_kelvin

    ln_mon = math.log(sodium)

    # Owczarzy 2004 Eq 3 — simplified monovalent correction
    inv_Tm = (1.0 / Tm_1M_kelvin
              + (4.29e-5 * f_GC - 3.95e-5) * ln_mon
              + 9.87e-6 * ln_mon ** 2)

    return 1.0 / inv_Tm if inv_Tm > 0 else Tm_1M_kelvin


# ═══════════════════════════════════════════════════════════════════════════
# Upgrade 3: Unified Owczarzy 2008 Mg²⁺ correction
# ═══════════════════════════════════════════════════════════════════════════
#
# Previously there were TWO different implementations of the Owczarzy 2008
# mixed regime (0.22 ≤ R < 6.0): one using an additive formula and one using
# linear interpolation.  Both were incorrect.
#
# The correct approach uses Eq 16 from Owczarzy et al. 2008, with
# monovalent-adjusted coefficients for the mixed regime.
#
# Reference: Owczarzy R, Moreira BG, You Y, et al. Predicting stability
# of DNA duplexes in solutions containing magnesium and monovalent cations.
# Biochemistry. 2008;47(19):5336-53.

def owczarzy_2008_mg_correction(
    Tm_1M_kelvin: float,
    Nbp: int,
    f_GC: float,
    mg_conc: float,
    monovalent_conc: float = 0.0,
) -> float:
    """Apply Owczarzy et al. 2008 magnesium correction to Tm.

    Implements the full Owczarzy 2008 algorithm with three regimes:

    1. **R < 0.22** (monovalent-dominated): Use Owczarzy 2004 Na⁺ correction.
    2. **0.22 ≤ R < 6.0** (mixed): Use Eq 16 with monovalent-adjusted
       coefficients (a_adj, d_adj, g_adj).
    3. **R ≥ 6.0** (magnesium-dominated): Use Eq 16 with base coefficients.

    where R = sqrt(Mg²⁺) / Mon⁺.

    Reference: Owczarzy R, Moreira BG, You Y, et al. Predicting stability
    of DNA duplexes in solutions containing magnesium and monovalent cations.
    Biochemistry. 2008;47(19):5336-53.

    Args:
        Tm_1M_kelvin: Tm at 1 M NaCl in Kelvin.
        Nbp: Number of base pairs in the duplex.
        f_GC: GC fraction of the sequence.
        mg_conc: Free Mg²⁺ concentration in M.
        monovalent_conc: Total monovalent cation concentration in M.

    Returns:
        Corrected Tm in Kelvin.
    """
    if mg_conc <= 0 and monovalent_conc <= 0:
        return Tm_1M_kelvin

    # ── No Mg²⁺: pure monovalent correction ────────────────────────────
    if mg_conc <= 0:
        return owczarzy_2004_na_correction(
            Tm_1M_kelvin, Nbp, f_GC, monovalent_conc)

    # ── Compute R = sqrt([Mg²⁺]) / [Mon⁺] ─────────────────────────────
    sqrt_Mg = math.sqrt(mg_conc)

    if monovalent_conc > 0:
        R = sqrt_Mg / monovalent_conc
    else:
        # No monovalent ions — treat as Mg-dominated regime
        R = float("inf")

    ln_Mg = math.log(mg_conc)

    # ── Regime 1: R < 0.22 — monovalent-dominated ──────────────────────
    if R < 0.22:
        # Use Owczarzy 2004 Na⁺ correction
        return owczarzy_2004_na_correction(
            Tm_1M_kelvin, Nbp, f_GC, monovalent_conc)

    # ── Base coefficients for Eq 16 (Owczarzy 2008 Table 4) ────────────
    a_base = 3.92e-5
    b_base = -9.11e-6
    c_base = 6.26e-5
    e_base = -4.82e-4
    f_base = 5.25e-4

    # ── Regime 2: 0.22 ≤ R < 6.0 — mixed regime ───────────────────────
    # Adjusted coefficients per Owczarzy 2008 Eq 16 with monovalent
    # ion competition adjustments.
    if R < 6.0:
        mon = monovalent_conc  # [Mon⁺] in molar
        ln_mon = math.log(mon) if mon > 0 else 0.0

        a_adj = 3.92e-5 * (0.843 - 0.352 * math.sqrt(mon) * ln_mon)
        d_adj = 1.42e-5 * (1.279 - 4.03e-3 * ln_mon
                            - 8.03e-3 * ln_mon ** 2)
        g_adj = 8.31e-5 * (0.486 - 0.258 * ln_mon
                            + 5.25e-3 * ln_mon ** 3)

        # b, c, e, f keep their base values
        inv_Tm_mg = (1.0 / Tm_1M_kelvin
                     + (a_adj + b_base * ln_Mg
                        + f_GC * (c_base + d_adj * ln_Mg)
                        + (1.0 / (2 * (Nbp - 1)))
                        * (e_base + f_base * ln_Mg + g_adj * ln_Mg ** 2)))
    else:
        # ── Regime 3: R ≥ 6.0 — magnesium-dominated ───────────────────
        # Use base coefficients (no monovalent adjustment)
        d_base = 1.42e-5
        g_base = 8.31e-5

        inv_Tm_mg = (1.0 / Tm_1M_kelvin
                     + (a_base + b_base * ln_Mg
                        + f_GC * (c_base + d_base * ln_Mg)
                        + (1.0 / (2 * (Nbp - 1)))
                        * (e_base + f_base * ln_Mg + g_base * ln_Mg ** 2)))

    if inv_Tm_mg > 0:
        return 1.0 / inv_Tm_mg
    else:
        return Tm_1M_kelvin


# ═══════════════════════════════════════════════════════════════════════════
# Main Tm calculation functions
# ═══════════════════════════════════════════════════════════════════════════

def calculate_tm(
    sequence: str,
    sodium: float = 0.05,
    potassium: float = 0.0,
    tris: float = 0.0,
    mg_concentration: float = 0.0,
    dntp_concentration: float = 0.0,
    primer_concentration: float = 5e-7,
) -> float:
    """Calculate melting temperature using the nearest-neighbor method.

    Uses SantaLucia & Hicks 2004 (DNA_NN4) thermodynamic parameters with:

    - Automatic self-complementary detection and corrected Tm formula
    - Owczarzy 2004/2008 salt corrections (unified mixed regime)
    - Ka-based dNTP-Mg²⁺ chelation model
    - Proper K⁺/Tris⁺ aggregation into monovalent pool

    For sequences shorter than 4 bp, returns 0.0.
    For sequences ≤ 8 bp, uses the Wallace rule as an approximation.

    Args:
        sequence: Primer sequence (5'→3', DNA bases only).
        sodium: Na⁺ concentration in M (default 50 mM).
        potassium: K⁺ concentration in M (default 0). K⁺ is interchangeable
            with Na⁺ per Owczarzy 2008.
        tris: Tris⁺ concentration in M (default 0). Tris⁺ counts fully
            toward the monovalent pool per Owczarzy 2008.
        mg_concentration: Total Mg²⁺ concentration in M (default 0).
        dntp_concentration: Total dNTP concentration in M (default 0).
            Free Mg²⁺ is computed using the Ka-based chelation model.
        primer_concentration: Total strand concentration in M (default 0.5 µM).

    Returns:
        Melting temperature in °C.
    """
    sequence = sequence.upper().replace("U", "T")

    # Short sequences: use simple formula
    if len(sequence) < 4:
        return 0.0

    if len(sequence) <= 8:
        return calculate_tm_wallace(sequence)

    # ── Compute thermodynamics ─────────────────────────────────────────
    thermo = compute_dh_ds(sequence)

    # ── Compute total monovalent concentration ─────────────────────────
    # Upgrade 5: Proper K⁺/Tris⁺ aggregation
    monovalent_conc = compute_monovalent_conc(sodium, potassium, tris)

    # ── Compute free Mg²⁺ with Ka-based chelation ─────────────────────
    # Upgrade 4: Ka-based dNTP chelation
    free_mg = compute_free_mg(mg_concentration, dntp_concentration)

    # ── Compute Tm at 1 M NaCl ────────────────────────────────────────
    ct = primer_concentration if primer_concentration > 0 else 5e-7
    Tm_1M = compute_tm_from_thermo(thermo.dH, thermo.dS, ct,
                                   thermo.is_self_comp)
    Tm_1M_kelvin = Tm_1M + 273.15

    # ── Apply salt correction ──────────────────────────────────────────
    if free_mg > 0:
        # Use Owczarzy 2008 Mg²⁺ correction (handles all three regimes)
        Tm_kelvin = owczarzy_2008_mg_correction(
            Tm_1M_kelvin, thermo.Nbp, thermo.f_GC,
            free_mg, monovalent_conc)
    elif monovalent_conc > 0:
        # Use Owczarzy 2004 Na⁺ correction
        Tm_kelvin = owczarzy_2004_na_correction(
            Tm_1M_kelvin, thermo.Nbp, thermo.f_GC, monovalent_conc)
    else:
        Tm_kelvin = Tm_1M_kelvin

    return round(Tm_kelvin - 273.15, 1)


def calculate_tm_detailed(
    sequence: str,
    sodium: float = 0.05,
    potassium: float = 0.0,
    tris: float = 0.0,
    mg_concentration: float = 0.0,
    dntp_concentration: float = 0.0,
    primer_concentration: float = 5e-7,
) -> TmResult:
    """Calculate melting temperature with detailed results.

    Same as :func:`calculate_tm` but returns a :class:`TmResult` object
    with all intermediate values for inspection.

    Args:
        sequence: Primer sequence (5'→3', DNA bases only).
        sodium: Na⁺ concentration in M (default 50 mM).
        potassium: K⁺ concentration in M (default 0).
        tris: Tris⁺ concentration in M (default 0).
        mg_concentration: Total Mg²⁺ concentration in M (default 0).
        dntp_concentration: Total dNTP concentration in M (default 0).
        primer_concentration: Total strand concentration in M (default 0.5 µM).

    Returns:
        TmResult with Tm and all intermediate values.
    """
    sequence = sequence.upper().replace("U", "T")

    thermo = compute_dh_ds(sequence)
    monovalent_conc = compute_monovalent_conc(sodium, potassium, tris)
    free_mg = compute_free_mg(mg_concentration, dntp_concentration)

    ct = primer_concentration if primer_concentration > 0 else 5e-7
    Tm_1M = compute_tm_from_thermo(thermo.dH, thermo.dS, ct,
                                   thermo.is_self_comp)
    Tm_1M_kelvin = Tm_1M + 273.15

    if free_mg > 0:
        Tm_kelvin = owczarzy_2008_mg_correction(
            Tm_1M_kelvin, thermo.Nbp, thermo.f_GC,
            free_mg, monovalent_conc)
        method = "Owczarzy2008"
    elif monovalent_conc > 0:
        Tm_kelvin = owczarzy_2004_na_correction(
            Tm_1M_kelvin, thermo.Nbp, thermo.f_GC, monovalent_conc)
        method = "Owczarzy2004"
    else:
        Tm_kelvin = Tm_1M_kelvin
        method = "SantaLucia2004_1M"

    salt_correction = (Tm_kelvin - 273.15) - Tm_1M
    Tm = Tm_kelvin - 273.15

    return TmResult(
        Tm=round(Tm, 1),
        Tm_1M=round(Tm_1M, 2),
        Tm_1M_kelvin=round(Tm_1M_kelvin, 4),
        dH=thermo.dH,
        dS=thermo.dS,
        Nbp=thermo.Nbp,
        f_GC=thermo.f_GC,
        is_self_comp=thermo.is_self_comp,
        salt_correction=round(salt_correction, 2),
        method=method,
    )


def calculate_tm_nearest_neighbor(
    sequence: str,
    na_conc: float = 0.05,
    mg_conc: float = 0.0,
    primer_conc: float = 5e-7,
) -> float:
    """Calculate Tm using nearest-neighbor thermodynamics (SantaLucia 2004).

    Backward-compatible interface matching the original primer_design.py API.
    For sequences shorter than 14 bp, falls back to the Wallace rule.
    For longer sequences, uses the full nearest-neighbor model with
    salt correction.

    Args:
        sequence: Primer sequence (5'→3', DNA bases only).
        na_conc: Na⁺ concentration in M (default 50 mM).
        mg_conc: Mg²⁺ concentration in M (default 0).
        primer_conc: Primer concentration in M (default 0.5 µM).

    Returns:
        Melting temperature in °C.
    """
    sequence = sequence.upper().replace("U", "T")
    if len(sequence) < 14:
        return calculate_tm_wallace(sequence)
    return calculate_tm(
        sequence,
        sodium=na_conc,
        mg_concentration=mg_conc,
        primer_concentration=primer_conc,
    )


def calculate_tm_wallace(sequence: str) -> float:
    """Calculate Tm using the simple Wallace rule: 2×(A+T) + 4×(G+C).

    Suitable only for short oligonucleotides (< 14 bp) in ideal conditions.
    For longer sequences, use :func:`calculate_tm_nearest_neighbor` instead.

    Args:
        sequence: Primer sequence (5'→3', DNA bases only).

    Returns:
        Melting temperature in °C.
    """
    sequence = sequence.upper().replace("U", "T")
    gc = sequence.count("G") + sequence.count("C")
    at = sequence.count("A") + sequence.count("T")
    return float(2 * at + 4 * gc)


# ═══════════════════════════════════════════════════════════════════════════
# Internal helpers
# ═══════════════════════════════════════════════════════════════════════════

def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def _compute_gc_content(seq: str) -> float:
    """Compute GC fraction of a DNA sequence."""
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)
