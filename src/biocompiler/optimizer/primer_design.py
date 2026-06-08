"""
Primer design constraints for DNA sequence optimization.

Ensures that the optimized sequence is compatible with PCR primer design
by controlling melting temperature, avoiding secondary structures,
and preventing primer dimer formation.

Based on IDT primer design guidelines and the SantaLucia (1998)
nearest-neighbor thermodynamic model for Tm calculation.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from ..type_system import CODON_TABLE, AA_TO_CODONS

logger = logging.getLogger(__name__)

__all__ = [
    "PrimerCheckResult",
    "PrimerDesignResult",
    "PrimerConstraintResult",
    "calculate_tm",
    "calculate_tm_nearest_neighbor",
    "calculate_tm_wallace",
    "check_gc_clamp",
    "check_self_complementarity",
    "check_heterodimer",
    "check_primer_pair",
    "design_primers",
    "evaluate_primer_constraint",
]

# ── Nearest-neighbor thermodynamic parameters (SantaLucia 1998) ───────
# ΔH in kcal/mol, ΔS in cal/(mol·K)

_NN_PARAMS: dict[str, tuple[float, float]] = {
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

# Initiation parameters
_INIT_GC_DH = 0.1   # kcal/mol (per G/C terminal)
_INIT_GC_DS = -2.8  # cal/(mol·K)
_INIT_AT_DH = 2.3   # kcal/mol (per A/T terminal)
_INIT_AT_DS = 4.1    # cal/(mol·K)

# Complementarity map (for self-complementarity / heterodimer checks)
_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "G": "C", "C": "G"}


@dataclass
class PrimerCheckResult:
    """Result of validating a primer pair.

    Attributes:
        valid: Whether the primer pair passes all checks.
        forward_tm: Melting temperature of the forward primer (°C).
        reverse_tm: Melting temperature of the reverse primer (°C).
        tm_diff: Absolute Tm difference between primers (°C).
        gc_clamp_forward: Whether the forward primer has a GC clamp.
        gc_clamp_reverse: Whether the reverse primer has a GC clamp.
        self_comp_forward: Self-complementary positions in forward primer.
        self_comp_reverse: Self-complementary positions in reverse primer.
        heterodimer_positions: Complementary positions between primers.
        issues: List of issues found.
    """

    valid: bool
    forward_tm: float
    reverse_tm: float
    tm_diff: float
    gc_clamp_forward: bool
    gc_clamp_reverse: bool
    self_comp_forward: list[tuple] = field(default_factory=list)
    self_comp_reverse: list[tuple] = field(default_factory=list)
    heterodimer_positions: list[tuple] = field(default_factory=list)
    issues: list[str] = field(default_factory=list)


@dataclass
class PrimerDesignResult:
    """Result of designing primers for a target region.

    Attributes:
        forward_primer: Forward primer sequence (5'→3').
        reverse_primer: Reverse primer sequence (5'→3', already reverse-complemented).
        forward_tm: Melting temperature of the forward primer (°C).
        reverse_tm: Melting temperature of the reverse primer (°C).
        forward_gc: GC content of the forward primer (0.0-1.0).
        reverse_gc: GC content of the reverse primer (0.0-1.0).
        product_length: Length of the PCR product in bp.
        issues: List of issues found during design.
    """

    forward_primer: str
    reverse_primer: str
    forward_tm: float
    reverse_tm: float
    forward_gc: float = 0.0
    reverse_gc: float = 0.0
    product_length: int = 0
    issues: list[str] = field(default_factory=list)


@dataclass
class PrimerConstraintResult:
    """Result of evaluating primer compatibility constraints on a sequence.

    Attributes:
        satisfied: Whether the sequence satisfies primer design constraints.
        forward_valid: Whether the forward primer region is valid.
        reverse_valid: Whether the reverse primer region is valid.
        issues: List of issues found.
    """

    satisfied: bool
    forward_valid: bool
    reverse_valid: bool
    issues: list[str] = field(default_factory=list)


def calculate_tm(
    sequence: str,
    na_concentration: float = 0.05,
    mg_concentration: float = 0.0,
    dntp_concentration: float = 0.0,
    primer_concentration: float = 0.0000005,
) -> float:
    """Calculate melting temperature using the nearest-neighbor method (SantaLucia 1998).

    Args:
        sequence: Primer sequence (5'→3', DNA bases only).
        na_concentration: Na⁺ concentration in M (default 50 mM).
        mg_concentration: Mg²⁺ concentration in M (default 0).
        dntp_concentration: dNTP concentration in M (default 0).
        primer_concentration: Primer concentration in M (default 0.5 µM).

    Returns:
        Melting temperature in °C.
    """
    sequence = sequence.upper().replace("U", "T")

    # Short sequences: use simple formula
    if len(sequence) < 4:
        return 0.0

    if len(sequence) <= 8:
        # For very short oligos, use Wallace rule as approximation
        gc = sequence.count("G") + sequence.count("C")
        at = sequence.count("A") + sequence.count("T")
        return float(2 * at + 4 * gc)

    # Nearest-neighbor method
    dh_total = 0.0  # kcal/mol
    ds_total = 0.0  # cal/(mol·K)

    # Add initiation parameters
    first_base = sequence[0]
    last_base = sequence[-1]

    for terminal in (first_base, last_base):
        if terminal in ("G", "C"):
            dh_total += _INIT_GC_DH
            ds_total += _INIT_GC_DS
        else:  # A or T
            dh_total += _INIT_AT_DH
            ds_total += _INIT_AT_DS

    # Sum nearest-neighbor parameters
    for i in range(len(sequence) - 1):
        dinuc = sequence[i:i + 2]
        params = _NN_PARAMS.get(dinuc)
        if params:
            dh_total += params[0]
            ds_total += params[1]

    # Salt correction (Owczarzy et al. 2004, simplified)
    # For monovalent cations
    if mg_concentration > 0 and na_concentration > 0:
        # Simplified: use effective monovalent concentration
        ratio = mg_concentration ** 0.5 / na_concentration
        if ratio < 0.22:
            # Use Na+ correction only
            salt_correction = 12.0 * (na_concentration ** 0.5 - 0.05) / (0.05 ** 0.5)
        else:
            # Use Mg2+ correction
            salt_correction = 4.29 * (mg_concentration ** 0.5) - 3.95
    elif na_concentration > 0:
        salt_correction = 12.0 * (na_concentration ** 0.5 - 0.05) / (0.05 ** 0.5)
    else:
        salt_correction = 0.0

    # Convert ΔH to cal/mol for consistent units with ΔS
    dh_cal = dh_total * 1000.0  # kcal to cal

    # Tm = ΔH / (ΔS + R * ln(Ct/4)) - 273.15 + salt_correction
    R = 1.987  # cal/(mol·K)
    ct = primer_concentration
    if ct <= 0:
        ct = 0.5e-6

    denominator = ds_total + R * (1.0 if ct / 4.0 <= 0 else __import__("math").log(ct / 4.0))
    if abs(denominator) < 1e-10:
        return 0.0

    tm = (dh_cal / denominator) - 273.15 + salt_correction
    return round(tm, 1)


def check_gc_clamp(
    sequence: str,
    min_gc_3prime: int = 1,
    window: int = 5,
) -> bool:
    """Check if the last bases have at least ``min_gc_3prime`` G/C bases (GC clamp).

    A GC clamp at the 3' end of a primer ensures stable binding during PCR
    initiation.  IDT recommends at least 1 G/C in the last 5 bases.

    Args:
        sequence: Primer sequence (5'→3').
        min_gc_3prime: Minimum number of G/C bases required in the 3' window.
        window: Number of bases at the 3' end to check.

    Returns:
        True if the GC clamp requirement is satisfied.
    """
    sequence = sequence.upper()
    if len(sequence) < window:
        window = len(sequence)
    if window == 0:
        return True

    tail = sequence[-window:]
    gc_count = tail.count("G") + tail.count("C")
    return gc_count >= min_gc_3prime


def check_self_complementarity(
    sequence: str,
    max_complement: int = 4,
) -> list[tuple[int, int]]:
    """Check for self-complementary regions that could form hairpins.

    Scans for complementary subsequences within the primer that could
    base-pair with each other, forming secondary structures.

    Args:
        sequence: Primer sequence (5'→3').
        max_complement: Maximum length of complementary stretch to flag.

    Returns:
        List of (start1, start2) tuples indicating complementary regions.
    """
    sequence = sequence.upper()
    comp = "".join(_COMPLEMENT.get(b, "N") for b in sequence)
    positions: list[tuple[int, int]] = []

    n = len(sequence)
    for length in range(max_complement, n // 2 + 1):
        for i in range(n - length + 1):
            subseq = sequence[i:i + length]
            # Check if this subsequence appears in the complement
            # (which means it can base-pair with another part of the primer)
            comp_strand = "".join(_COMPLEMENT.get(b, "N") for b in subseq)
            # Reverse the complement strand (hairpin formation)
            comp_rev = comp_strand[::-1]
            for j in range(i + length, n - length + 1):
                if sequence[j:j + length] == comp_rev:
                    positions.append((i, j))

    # Deduplicate — keep only the longest non-overlapping pairs
    return _deduplicate_complement_positions(positions)


def check_heterodimer(
    seq1: str,
    seq2: str,
    max_complement: int = 4,
) -> list[tuple[int, int, int, int]]:
    """Check for complementarity between two primers (primer dimer formation).

    Args:
        seq1: First primer sequence (5'→3').
        seq2: Second primer sequence (5'→3').
        max_complement: Maximum length of complementary stretch to flag.

    Returns:
        List of (start1, end1, start2, end2) tuples for complementary regions.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Reverse complement of seq2 for alignment checking
    seq2_rc = "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq2))

    positions: list[tuple[int, int, int, int]] = []

    n1 = len(seq1)
    n2_rc = len(seq2_rc)

    for length in range(max_complement, min(n1, n2_rc) + 1):
        for i in range(n1 - length + 1):
            subseq = seq1[i:i + length]
            for j in range(n2_rc - length + 1):
                if seq2_rc[j:j + length] == subseq:
                    # Map j back to original seq2 position
                    orig_j = n2_rc - j - length
                    positions.append((i, i + length, orig_j, orig_j + length))

    return positions


def check_primer_pair(
    forward: str,
    reverse: str,
    min_tm: float = 55.0,
    max_tm: float = 65.0,
    max_tm_diff: float = 5.0,
    max_length: int = 30,
    min_length: int = 18,
) -> PrimerCheckResult:
    """Comprehensive primer pair validation.

    Validates:
    - Primer length within bounds
    - Melting temperature within range
    - Tm difference between primers
    - GC clamp at 3' end
    - Self-complementarity (hairpin potential)
    - Heterodimer formation

    Args:
        forward: Forward primer sequence (5'→3').
        reverse: Reverse primer sequence (5'→3').
        min_tm: Minimum acceptable Tm (°C).
        max_tm: Maximum acceptable Tm (°C).
        max_tm_diff: Maximum Tm difference between primers (°C).
        max_length: Maximum primer length.
        min_length: Minimum primer length.

    Returns:
        PrimerCheckResult with validation details.
    """
    forward = forward.upper()
    reverse = reverse.upper()
    issues: list[str] = []

    # Length checks
    for name, seq in [("Forward", forward), ("Reverse", reverse)]:
        if len(seq) < min_length:
            issues.append(f"{name} primer too short: {len(seq)} < {min_length}")
        if len(seq) > max_length:
            issues.append(f"{name} primer too long: {len(seq)} > {max_length}")

    # Tm checks
    fwd_tm = calculate_tm(forward)
    rev_tm = calculate_tm(reverse)
    tm_diff = abs(fwd_tm - rev_tm)

    if fwd_tm < min_tm:
        issues.append(f"Forward Tm too low: {fwd_tm:.1f}°C < {min_tm}°C")
    if fwd_tm > max_tm:
        issues.append(f"Forward Tm too high: {fwd_tm:.1f}°C > {max_tm}°C")
    if rev_tm < min_tm:
        issues.append(f"Reverse Tm too low: {rev_tm:.1f}°C < {min_tm}°C")
    if rev_tm > max_tm:
        issues.append(f"Reverse Tm too high: {rev_tm:.1f}°C > {max_tm}°C")
    if tm_diff > max_tm_diff:
        issues.append(f"Tm difference too large: {tm_diff:.1f}°C > {max_tm_diff}°C")

    # GC clamp checks
    gc_clamp_fwd = check_gc_clamp(forward)
    gc_clamp_rev = check_gc_clamp(reverse)
    if not gc_clamp_fwd:
        issues.append("Forward primer lacks GC clamp (need ≥1 G/C in last 5 bases)")
    if not gc_clamp_rev:
        issues.append("Reverse primer lacks GC clamp (need ≥1 G/C in last 5 bases)")

    # Self-complementarity checks
    self_comp_fwd = check_self_complementarity(forward)
    self_comp_rev = check_self_complementarity(reverse)
    if self_comp_fwd:
        issues.append(f"Forward primer has self-complementary regions: {len(self_comp_fwd)} found")
    if self_comp_rev:
        issues.append(f"Reverse primer has self-complementary regions: {len(self_comp_rev)} found")

    # Heterodimer check
    hetero = check_heterodimer(forward, reverse)
    if hetero:
        issues.append(f"Primer dimer detected: {len(hetero)} complementary regions")

    return PrimerCheckResult(
        valid=len(issues) == 0,
        forward_tm=fwd_tm,
        reverse_tm=rev_tm,
        tm_diff=tm_diff,
        gc_clamp_forward=gc_clamp_fwd,
        gc_clamp_reverse=gc_clamp_rev,
        self_comp_forward=self_comp_fwd,
        self_comp_reverse=self_comp_rev,
        heterodimer_positions=hetero,
        issues=issues,
    )


def design_primers(
    sequence: str,
    target_region: tuple[int, int],
    organism: str = "",
    min_tm: float = 58.0,
    max_tm: float = 62.0,
    primer_length: int = 20,
) -> PrimerDesignResult:
    """Design primer pair for a target region.

    Designs forward and reverse primers flanking the target region,
    adjusting length to achieve the desired melting temperature.

    Args:
        sequence: Template DNA sequence.
        target_region: (start, end) of the target region (0-based).
        organism: Target organism (for codon context).
        min_tm: Minimum acceptable Tm (°C).
        max_tm: Maximum acceptable Tm (°C).
        primer_length: Starting primer length (adjusted for Tm).

    Returns:
        PrimerDesignResult with designed primer pair.
    """
    sequence = sequence.upper()
    issues: list[str] = []

    target_start, target_end = target_region

    # Design forward primer (upstream of target)
    fwd_start = max(0, target_start - primer_length)
    fwd_end = target_start
    fwd_primer = sequence[fwd_start:fwd_end]

    # Adjust forward primer length for Tm
    fwd_tm = calculate_tm(fwd_primer)
    for length_adj in range(0, 15):
        trial_length = primer_length + length_adj
        trial_start = max(0, target_start - trial_length)
        trial_primer = sequence[trial_start:target_start]
        if len(trial_primer) < trial_length:
            break
        trial_tm = calculate_tm(trial_primer)
        if min_tm <= trial_tm <= max_tm:
            fwd_primer = trial_primer
            fwd_start = trial_start
            fwd_tm = trial_tm
            break
    else:
        issues.append(f"Could not achieve target Tm for forward primer (Tm={fwd_tm:.1f}°C)")

    # Design reverse primer (downstream of target)
    rev_start = target_end
    rev_end = min(len(sequence), target_end + primer_length)
    rev_template = sequence[rev_start:rev_end]

    # Reverse complement for the actual primer
    rev_primer = _reverse_complement(rev_template)
    rev_tm = calculate_tm(rev_primer)

    for length_adj in range(0, 15):
        trial_length = primer_length + length_adj
        trial_end = min(len(sequence), target_end + trial_length)
        trial_template = sequence[target_end:trial_end]
        trial_primer = _reverse_complement(trial_template)
        if len(trial_primer) < trial_length:
            break
        trial_tm = calculate_tm(trial_primer)
        if min_tm <= trial_tm <= max_tm:
            rev_primer = trial_primer
            rev_start = target_end
            rev_end = trial_end
            rev_tm = trial_tm
            break
    else:
        issues.append(f"Could not achieve target Tm for reverse primer (Tm={rev_tm:.1f}°C)")

    product_length = rev_end - fwd_start
    forward_gc = _compute_gc_content(fwd_primer)
    reverse_gc = _compute_gc_content(rev_primer)

    return PrimerDesignResult(
        forward_primer=fwd_primer,
        reverse_primer=rev_primer,
        forward_tm=fwd_tm,
        reverse_tm=rev_tm,
        forward_gc=forward_gc,
        reverse_gc=reverse_gc,
        product_length=product_length,
        issues=issues,
    )


def evaluate_primer_constraint(
    sequence: str,
    region_start: int,
    region_end: int,
    min_tm: float = 55.0,
    max_tm: float = 65.0,
) -> PrimerConstraintResult:
    """Evaluate if the sequence at the given region satisfies primer design constraints.

    Designs primers flanking the region and checks if they meet
    Tm, GC clamp, and secondary structure requirements.

    Args:
        sequence: Template DNA sequence.
        region_start: Start of the region (0-based).
        region_end: End of the region (0-based, exclusive).
        min_tm: Minimum acceptable Tm (°C).
        max_tm: Maximum acceptable Tm (°C).

    Returns:
        PrimerConstraintResult with evaluation details.
    """
    result = design_primers(
        sequence, (region_start, region_end),
        min_tm=min_tm, max_tm=max_tm,
    )

    check = check_primer_pair(
        result.forward_primer,
        result.reverse_primer,
        min_tm=min_tm, max_tm=max_tm,
    )

    forward_valid = (min_tm <= result.forward_tm <= max_tm
                     and check.gc_clamp_forward
                     and len(check.self_comp_forward) == 0)
    reverse_valid = (min_tm <= result.reverse_tm <= max_tm
                     and check.gc_clamp_reverse
                     and len(check.self_comp_reverse) == 0)

    return PrimerConstraintResult(
        satisfied=forward_valid and reverse_valid,
        forward_valid=forward_valid,
        reverse_valid=reverse_valid,
        issues=result.issues + check.issues,
    )


# ── Internal helpers ─────────────────────────────────────────────────


def calculate_tm_nearest_neighbor(
    sequence: str,
    na_conc: float = 0.05,
    mg_conc: float = 0.0,
    primer_conc: float = 5e-7,
) -> float:
    """Calculate Tm using nearest-neighbor thermodynamics (SantaLucia 1998).

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
        na_concentration=na_conc,
        mg_concentration=mg_conc,
        primer_concentration=primer_conc,
    )


def calculate_tm_wallace(sequence: str) -> float:
    """Calculate Tm using the simple Wallace rule: 2*(A+T) + 4*(G+C).

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


def _compute_gc_content(seq: str) -> float:
    """Compute GC fraction of a DNA sequence."""
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def _deduplicate_complement_positions(
    positions: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Deduplicate complementary position pairs, keeping the longest."""
    if not positions:
        return []

    # Sort by distance between positions (longer hairpins first)
    positions.sort(key=lambda p: abs(p[1] - p[0]), reverse=True)

    deduped: list[tuple[int, int]] = []
    used: set[int] = set()

    for p in positions:
        if p[0] not in used and p[1] not in used:
            deduped.append(p)
            used.add(p[0])
            used.add(p[1])

    return sorted(deduped)
