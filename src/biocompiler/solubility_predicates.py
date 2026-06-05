"""
BioCompiler Solubility Predicates v9.2.0
=========================================
Type system predicates for protein solubility assessment.

Uses CamSol-style intrinsic solubility scoring and biophysical heuristics
to predict whether a protein will express in soluble form or aggregate.

Predicates:
  - SolubleExpression: protein predicted soluble via CamSol scoring
  - NoAggregationProneRegion: no long aggregation-prone regions
  - ChargeComposition: adequate charged residue fraction and safe pI
  - NoLongHydrophobicStretch: no excessively long hydrophobic stretches

Helper functions:
  - compute_approximate_pI: estimate isoelectric point via bisection
  - compute_net_charge: net charge at a given pH (Henderson-Hasselbalch)
  - find_hydrophobic_stretches: detect maximal consecutive hydrophobic runs
"""

from __future__ import annotations

from .type_system import Verdict, TypeCheckResult

__all__ = [
    "PKA_VALUES",
    "compute_net_charge",
    "compute_approximate_pI",
    "find_hydrophobic_stretches",
    "evaluate_soluble_expression",
    "evaluate_no_aggregation_prone_region",
    "evaluate_charge_composition",
    "evaluate_no_long_hydrophobic_stretch",
]

# ────────────────────────────────────────────────────────────
# pKa values for ionizable groups (charge calculation)
# ────────────────────────────────────────────────────────────
PKA_VALUES: dict[str, float] = {
    "N_term": 9.69,   # N-terminal amino group
    "C_term": 2.34,   # C-terminal carboxyl group
    "K": 10.54,       # Lysine side chain
    "R": 12.48,       # Arginine side chain
    "H": 6.04,        # Histidine side chain
    "D": 3.90,        # Aspartic acid side chain
    "E": 4.07,        # Glutamic acid side chain
    "C": 8.28,        # Cysteine side chain
    "Y": 10.07,       # Tyrosine side chain
}

# ────────────────────────────────────────────────────────────
# CamSol intrinsic solubility profile (per-residue)
# Based on amino acid physicochemical properties.
# Positive values → soluble propensity; negative → aggregation-prone.
# Scaled so that overall score × 3.0 aligns with CamSol thresholds.
# ────────────────────────────────────────────────────────────
_CAMSOL_INTRINSIC: dict[str, float] = {
    "A":  0.10, "R":  1.30, "N":  0.60, "D":  1.20, "C": -0.50,
    "Q":  0.50, "E":  1.40, "G":  0.00, "H":  0.30, "I": -1.50,
    "L": -1.30, "K":  1.50, "M": -0.80, "F": -1.60, "P":  0.20,
    "S":  0.50, "T":  0.10, "W": -1.70, "Y": -0.90, "V": -1.20,
}

# Scaling factor: converts mean smoothed profile to CamSol-score range
# where > 1.5 = highly soluble, < -1.0 = insoluble.
_CAMSOL_SCALE = 3.0

# Default smoothing window for CamSol profile
_CAMSOL_WINDOW = 7

# Default hydrophobic residue set (AILMFWV)
_DEFAULT_HYDROPHOBIC: set[str] = set("AILMFWV")

# Membrane protein residue set (strongly hydrophobic transmembrane indicators)
_MEMBRANE_HYDROPHOBIC: set[str] = set("AILMFWV")

# Minimum consecutive hydrophobic residues to be flagged as aggregation-prone
_AGG_MIN_CONSECUTIVE_HYDROPHOBIC: int = 6

# Maximum aggregation score cutoff for FAIL verdict
_AGG_FAIL_SCORE_CUTOFF: float = 1.5

# N-terminal aggregation penalty factor (N-term aggregation is worse)
_AGG_NTERM_PENALTY: float = 1.5

# Hydrophobic fraction thresholds for solubility
_HYDROPHOBIC_FRACTION_UNCERTAIN_LO: float = 0.35  # > this → UNCERTAIN zone
_HYDROPHOBIC_FRACTION_FAIL: float = 0.55           # > this → FAIL zone

# Net charge thresholds at pH 7
_NET_CHARGE_PASS: float = 2.0
_NET_CHARGE_WARN: float = 1.0

# Charge/length ratio thresholds
_CHARGE_RATIO_LO: float = -0.15
_CHARGE_RATIO_HI: float = 0.15

# Organism-specific threshold adjustments
_ORGANISM_THRESHOLDS: dict[str, dict] = {
    "ecoli": {
        "net_charge_min_offset": 1.0,    # E. coli needs more positive charge
        "charge_ratio_lo": -0.20,        # E. coli prefers slightly negative
        "charge_ratio_hi": 0.10,
    },
    "mammalian": {
        "net_charge_min_offset": -0.5,    # Mammalian needs less positive charge
        "charge_ratio_lo": -0.10,         # Mammalian prefers more neutral
        "charge_ratio_hi": 0.15,
    },
}

# ────────────────────────────────────────────────────────────
# Verdict threshold constants
# ────────────────────────────────────────────────────────────

# CamSol overall score thresholds for solubility verdicts
_CAMSOL_HIGHLY_SOLUBLE: float = 1.5   # score > this → PASS
_CAMSOL_MARGINAL: float = -1.0       # score < this → LIKELY_FAIL

# Aggregation-prone region length thresholds (residues)
_AGG_BORDERLINE_MAX: int = 7    # max_region_length+1 … 7 → LIKELY_PASS
_AGG_UNCERTAIN_MAX: int = 10   # 8 … 10 → UNCERTAIN
_AGG_LIKELY_FAIL_MAX: int = 15 # 11 … 15 → LIKELY_FAIL; >15 → FAIL

# Hydrophobic stretch excess thresholds (residues beyond max_stretch)
_HYDRO_EXCESS_BORDERLINE: int = 3  # excess 1…3 → LIKELY_PASS
_HYDRO_EXCESS_UNCERTAIN: int = 6   # excess 4…6 → UNCERTAIN; >6 → FAIL

# Bisection parameters for approximate pI computation
_PH_MIN: float = 0.0
_PH_MAX: float = 14.0
_BISECTION_ITERATIONS: int = 100


# ────────────────────────────────────────────────────────────
# Internal CamSol helpers
# ────────────────────────────────────────────────────────────

def _camsol_smoothed_profile(protein: str, window: int = _CAMSOL_WINDOW) -> list[float]:
    """Compute per-residue CamSol smoothed solubility profile.

    Applies a sliding-window average over the raw intrinsic scores,
    then scales by ``_CAMSOL_SCALE`` so the values align with CamSol
    threshold conventions (overall > 1.5 = highly soluble).

    Args:
        protein: Upper-cased amino-acid sequence.
        window: Smoothing window width (odd number preferred).

    Returns:
        List of smoothed, scaled per-residue scores.  Length == len(protein).
    """
    n = len(protein)
    if n == 0:
        return []

    raw = [_CAMSOL_INTRINSIC.get(aa, 0.0) for aa in protein]
    half_w = window // 2

    smoothed: list[float] = []
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        avg = sum(raw[start:end]) / (end - start)
        smoothed.append(avg * _CAMSOL_SCALE)

    return smoothed


def _camsol_overall_score(protein: str, window: int = _CAMSOL_WINDOW) -> float:
    """Compute overall CamSol intrinsic solubility score.

    Args:
        protein: Upper-cased amino-acid sequence.
        window: Smoothing window width.

    Returns:
        Overall CamSol score (typical range -3 to +3).
        > 1.5 highly soluble, 0–1.5 soluble, -1–0 marginal, < -1 insoluble.
    """
    profile = _camsol_smoothed_profile(protein, window)
    if not profile:
        return 0.0
    return sum(profile) / len(profile)


def _find_aggregation_regions(
    protein: str,
    window: int = _CAMSOL_WINDOW,
    score_threshold: float = -1.0,
    min_consecutive: int = _AGG_MIN_CONSECUTIVE_HYDROPHOBIC,
) -> list[tuple[int, int, float]]:
    """Find aggregation-prone regions from CamSol per-residue scores.

    A region is a maximal consecutive run of residues whose smoothed
    CamSol score is below *score_threshold* AND has at least
    *min_consecutive* residues in the run.  Short runs are filtered
    out to reduce false positives from transient hydrophobic patches.

    Args:
        protein: Upper-cased amino-acid sequence.
        window: Smoothing window width.
        score_threshold: Per-residue score below which a residue is
            considered aggregation-prone.
        min_consecutive: Minimum number of consecutive low-scoring
            residues to qualify as an aggregation-prone region.

    Returns:
        List of (start, end, avg_score) tuples (end is exclusive).
    """
    profile = _camsol_smoothed_profile(protein, window)
    if not profile:
        return []

    # First find all maximal runs below threshold
    raw_regions: list[tuple[int, int, float]] = []
    i = 0
    n = len(profile)
    while i < n:
        if profile[i] < score_threshold:
            start = i
            while i < n and profile[i] < score_threshold:
                i += 1
            end = i
            raw_regions.append((start, end))
        else:
            i += 1

    # Filter to regions meeting the minimum consecutive residue requirement
    regions: list[tuple[int, int, float]] = []
    for start, end in raw_regions:
        length = end - start
        if length >= min_consecutive:
            avg_score = sum(profile[start:end]) / length
            regions.append((start, end, round(avg_score, 3)))

    return regions


# ────────────────────────────────────────────────────────────
# Public helper functions
# ────────────────────────────────────────────────────────────

def compute_net_charge(protein: str, pH: float) -> float:
    """Compute the net charge of a protein at a given pH.

    Uses the Henderson-Hasselbalch equation for each ionizable group.

    Positive contributions (protonated at low pH):
        N-terminus, K, R, H  →  count × 1 / (1 + 10^(pH − pKa))

    Negative contributions (deprotonated at high pH):
        C-terminus, D, E     →  −count × 1 / (1 + 10^(pKa − pH))

    Args:
        protein: Amino-acid sequence (single-letter codes).
        pH: pH value at which to compute net charge.

    Returns:
        Net charge (float).  Positive = basic; negative = acidic.
    """
    if not protein:
        return 0.0

    protein = protein.upper()
    charge = 0.0

    # Positive groups
    # N-terminal amino group
    charge += 1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["N_term"]))

    for aa in protein:
        if aa == "K":
            charge += 1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["K"]))
        elif aa == "R":
            charge += 1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["R"]))
        elif aa == "H":
            charge += 1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["H"]))

    # Negative groups
    # C-terminal carboxyl group
    charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["C_term"] - pH))

    for aa in protein:
        if aa == "D":
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["D"] - pH))
        elif aa == "E":
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["E"] - pH))
        elif aa == "C":
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["C"] - pH))
        elif aa == "Y":
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["Y"] - pH))

    return charge


def compute_approximate_pI(protein: str) -> float:
    """Compute approximate isoelectric point (pI) of a protein.

    The pI is the pH at which the net charge is zero.  Uses bisection
    over [0, 14] with 100 iterations (~1e-30 precision) to find the
    crossover.

    Args:
        protein: Amino-acid sequence (single-letter codes).

    Returns:
        Approximate pI value (float, in range [0, 14]).
        Returns 7.0 for empty input.
    """
    if not protein:
        return 7.0

    lo, hi = _PH_MIN, _PH_MAX
    for _ in range(_BISECTION_ITERATIONS):
        mid = (lo + hi) / 2.0
        charge = compute_net_charge(protein, mid)
        if charge > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def find_hydrophobic_stretches(
    protein: str,
    hydrophobic: set[str] | None = None,
) -> list[tuple[int, int]]:
    """Find all maximal consecutive hydrophobic stretches.

    A stretch is a maximal run of residues where every residue belongs
    to the *hydrophobic* set.

    Args:
        protein: Amino-acid sequence (single-letter codes).
        hydrophobic: Set of single-letter hydrophobic residue codes.
            Defaults to ``{'A', 'I', 'L', 'M', 'F', 'W', 'V'}``.

    Returns:
        List of (start, end) tuples (end exclusive) for every maximal
        hydrophobic stretch of length ≥ 1.
    """
    if not protein:
        return []

    protein = protein.upper()
    hydro = hydrophobic if hydrophobic is not None else _DEFAULT_HYDROPHOBIC

    stretches: list[tuple[int, int]] = []
    start: int | None = None

    for i, aa in enumerate(protein):
        if aa in hydro:
            if start is None:
                start = i
        else:
            if start is not None:
                stretches.append((start, i))
                start = None

    # Handle stretch extending to the C-terminus
    if start is not None:
        stretches.append((start, len(protein)))

    return stretches


# ────────────────────────────────────────────────────────────
# Organism key helper
# ────────────────────────────────────────────────────────────

def _organism_key(organism: str) -> str:
    """Normalize organism name to a lookup key for threshold tables."""
    org_lower = organism.lower().strip()
    if any(tok in org_lower for tok in ("ecoli", "e.coli", "e. coli", "escherichia")):
        return "ecoli"
    if any(tok in org_lower for tok in ("mammal", "human", "mouse", "cho", "hek", "hepg2")):
        return "mammalian"
    return "default"


def _get_org_threshold(organism: str, param: str, default: float) -> float:
    """Retrieve an organism-specific threshold parameter, falling back to *default*."""
    key = _organism_key(organism)
    if key in _ORGANISM_THRESHOLDS and param in _ORGANISM_THRESHOLDS[key]:
        return _ORGANISM_THRESHOLDS[key][param]
    return default


# ────────────────────────────────────────────────────────────
# Membrane protein heuristic
# ────────────────────────────────────────────────────────────

def _is_likely_membrane_protein(protein: str) -> bool:
    """Heuristic: detect if a protein is likely a membrane protein.

    Returns True if there are ≥2 hydrophobic stretches of ≥19 residues
    (typical transmembrane helix length), suggesting a membrane protein
    with multiple transmembrane domains.
    """
    stretches = find_hydrophobic_stretches(protein, hydrophobic=_MEMBRANE_HYDROPHOBIC)
    long_stretches = [s for s in stretches if (s[1] - s[0]) >= 19]
    return len(long_stretches) >= 2


# Signal peptide detection constants
_SIGNAL_PEPTIDE_MAX_LENGTH: int = 30   # Typical signal peptide: 15-30 residues
_SIGNAL_PEPTIDE_MIN_HYDRO: int = 7     # Minimum hydrophobic residues in core
_SIGNAL_PEPTIDE_NTERM_WINDOW: int = 30 # Only look in first 30 residues

def _detect_signal_peptide(protein: str) -> tuple[int, int] | None:
    """Detect a predicted N-terminal signal peptide.

    Signal peptides are short (15-30 residue) N-terminal sequences with a
    hydrophobic core.  This heuristic checks for a hydrophobic stretch of
    ≥7 residues starting within the first 30 residues of the protein.

    Returns:
        (start, end) of the predicted signal peptide (end exclusive),
        or None if no signal peptide is detected.
    """
    if not protein:
        return None

    # Only examine the N-terminal region
    nterm = protein[:_SIGNAL_PEPTIDE_NTERM_WINDOW]
    stretches = find_hydrophobic_stretches(nterm, hydrophobic=_DEFAULT_HYDROPHOBIC)

    for start, end in stretches:
        length = end - start
        if length >= _SIGNAL_PEPTIDE_MIN_HYDRO and start < _SIGNAL_PEPTIDE_MAX_LENGTH:
            # Extend to typical signal peptide boundaries (include charged n-region)
            sp_start = max(0, start - 3)
            sp_end = min(end + 5, _SIGNAL_PEPTIDE_MAX_LENGTH, len(protein))
            return (sp_start, sp_end)

    return None


# ────────────────────────────────────────────────────────────
# Predicate 1: Soluble Expression
# ────────────────────────────────────────────────────────────

def evaluate_soluble_expression(
    sequence: str,
    protein: str,
    organism: str,
    min_solubility_score: float = 0.0,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check if protein is predicted to be soluble using CamSol scoring,
    hydrophobic fraction, and net charge at pH 7.

    Verdict logic considers three signals with organism-specific thresholds:

    1. **CamSol intrinsic score** (primary):
       - Score > 1.5   → PASS
       - Score 0.0–1.5 → LIKELY_PASS
       - Score -1.0–0.0 → UNCERTAIN (moderate / WARN zone)
       - Score < -1.0   → LIKELY_FAIL

    2. **Hydrophobic fraction** (PASS if ≤ 0.35):
       - ≤ 0.35            → no penalty
       - 0.35 to 0.55      → UNCERTAIN (downgrade to UNCERTAIN at most)
       - > 0.55            → additional FAIL signal

    3. **Net charge at pH 7** (PASS if |net_charge| > 2):
       - |net_charge| > 2   → no penalty
       - 1 < |net_charge| ≤ 2 → WARN (UNCERTAIN)
       - |net_charge| ≤ 1    → additional FAIL signal

    Organism-specific adjustments:
       - E. coli: needs more positive charge (offset +1.0)
       - Mammalian: less positive charge needed (offset -0.5)

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (single-letter codes).
        organism: Target organism name.
        min_solubility_score: Minimum acceptable CamSol score (default 0.0).
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict and CamSol score derivation.
    """
    protein = protein.upper()

    if not protein:
        return TypeCheckResult(
            predicate="SolubleExpression",
            verdict=Verdict.FAIL,
            violation="Empty protein sequence",
        )

    n = len(protein)

    # --- CamSol intrinsic score (primary signal) ---
    camsol_score = _camsol_overall_score(protein)

    # --- Hydrophobic fraction ---
    hydrophobic_count = sum(1 for aa in protein if aa in _DEFAULT_HYDROPHOBIC)
    hydrophobic_fraction = hydrophobic_count / n if n > 0 else 0.0

    # --- Net charge at pH 7 ---
    net_charge = compute_net_charge(protein, 7.0)
    abs_net_charge = abs(net_charge)

    # Compute charged fraction for downstream adjustments
    pos_count = sum(1 for aa in protein if aa in {"K", "R", "H"})
    neg_count = sum(1 for aa in protein if aa in {"D", "E"})
    charged_fraction = (pos_count + neg_count) / n if n > 0 else 0.0

    # --- Organism-specific thresholds ---
    org_charge_offset = _get_org_threshold(organism, "net_charge_min_offset", 0.0)
    effective_net_charge_min = _NET_CHARGE_PASS + org_charge_offset

    # --- Determine verdict from CamSol score (primary) ---
    if camsol_score > _CAMSOL_HIGHLY_SOLUBLE:
        verdict = Verdict.PASS
        violation = None
    elif camsol_score >= 0.0:
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif camsol_score >= _CAMSOL_MARGINAL:
        # Marginal CamSol score: use charge composition to refine.
        # If the protein has adequate net charge, it's more likely soluble
        # even with a borderline CamSol score → LIKELY_PASS.
        if abs_net_charge > _NET_CHARGE_WARN:
            verdict = Verdict.LIKELY_PASS
            violation = None
        else:
            verdict = Verdict.LIKELY_FAIL
            violation = (
                f"Marginal solubility: CamSol score {camsol_score:.3f} "
                f"in [{_CAMSOL_MARGINAL}, 0.0) with low net charge ({net_charge:.1f})"
            )
    else:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Insoluble protein: CamSol score {camsol_score:.3f} < {_CAMSOL_MARGINAL}"
        )

    # --- Hydrophobic fraction adjustments ---
    hydro_warnings: list[str] = []
    if hydrophobic_fraction > _HYDROPHOBIC_FRACTION_FAIL:
        # Strong failure signal (hydrophobic_fraction > 0.55) — downgrade by up to 2 levels
        # Use LIKELY_FAIL instead of UNCERTAIN since high hydrophobicity
        # is a meaningful signal, not ambiguous.
        if verdict == Verdict.PASS:
            verdict = Verdict.LIKELY_FAIL
        elif verdict == Verdict.LIKELY_PASS:
            verdict = Verdict.LIKELY_FAIL
        hydro_warnings.append(
            f"High hydrophobic fraction ({hydrophobic_fraction:.1%} > "
            f"{_HYDROPHOBIC_FRACTION_FAIL:.0%})"
        )
    elif hydrophobic_fraction > _HYDROPHOBIC_FRACTION_UNCERTAIN_LO:
        # Elevated hydrophobic zone (0.35 < hydrophobic_fraction <= 0.55)
        # Use LIKELY_FAIL instead of UNCERTAIN since elevated hydrophobicity
        # is a meaningful signal. However, if charged fraction is high
        # (>25%), the charges may compensate → only downgrade by 1 level.
        if hydrophobic_fraction <= 0.45 and charged_fraction > 0.25:
            # Moderate hydrophobicity with good charge compensation
            if verdict == Verdict.PASS:
                verdict = Verdict.LIKELY_PASS
            elif verdict == Verdict.LIKELY_PASS:
                verdict = Verdict.LIKELY_FAIL
        else:
            if verdict == Verdict.PASS:
                verdict = Verdict.LIKELY_FAIL
            elif verdict == Verdict.LIKELY_PASS:
                verdict = Verdict.LIKELY_FAIL
        hydro_warnings.append(
            f"Elevated hydrophobic fraction ({hydrophobic_fraction:.1%} > "
            f"{_HYDROPHOBIC_FRACTION_UNCERTAIN_LO:.0%})"
        )

    # --- Net charge adjustments ---
    # Note: net charge alone is not a reliable solubility indicator when the
    # protein has many charged residues (high charged_fraction). A protein
    # with balanced K/D pairs (net charge ~0) is perfectly soluble. Only
    # downgrade when BOTH net charge is low AND charged fraction is low.
    charge_warnings: list[str] = []
    if abs_net_charge <= _NET_CHARGE_WARN and charged_fraction < 0.15:
        # Very low charge AND few charged residues — meaningful signal
        if verdict == Verdict.PASS:
            verdict = Verdict.LIKELY_FAIL
        elif verdict == Verdict.LIKELY_PASS:
            verdict = Verdict.LIKELY_FAIL
        charge_warnings.append(
            f"Very low net charge at pH 7 (|{net_charge:.1f}| ≤ {_NET_CHARGE_WARN}) "
            f"with low charged fraction ({charged_fraction:.1%})"
        )
    elif abs_net_charge <= effective_net_charge_min and charged_fraction < 0.15:
        # Moderate charge AND few charged residues — downgrade by 1 level
        if verdict == Verdict.PASS:
            verdict = Verdict.LIKELY_PASS
        elif verdict == Verdict.LIKELY_PASS:
            verdict = Verdict.LIKELY_FAIL
        charge_warnings.append(
            f"Low net charge at pH 7 (|{net_charge:.1f}| ≤ {effective_net_charge_min:.1f}) "
            f"with low charged fraction ({charged_fraction:.1%})"
        )

    # If the score is below the user-specified minimum, that's a stronger
    # failure signal
    if camsol_score < min_solubility_score and verdict in (
        Verdict.PASS, Verdict.LIKELY_PASS,
    ):
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"CamSol score {camsol_score:.3f} is below the minimum "
            f"acceptable score {min_solubility_score}"
        )

    # Combine violation messages
    all_warnings = hydro_warnings + charge_warnings
    if all_warnings and violation:
        violation = f"{violation}; {'; '.join(all_warnings)}"
    elif all_warnings:
        violation = "; ".join(all_warnings)

    # Identify aggregation-prone regions for derivation
    agg_regions = _find_aggregation_regions(protein)

    # Build derivation
    derivation: list[dict] = [
        {"step": "camsol_intrinsic_score", "value": round(camsol_score, 3)},
        {"step": "min_solubility_score", "value": min_solubility_score},
        {"step": "hydrophobic_fraction", "value": round(hydrophobic_fraction, 4)},
        {"step": "hydrophobic_fraction_pass_threshold", "value": _HYDROPHOBIC_FRACTION_UNCERTAIN_LO},
        {"step": "net_charge_pH7", "value": round(net_charge, 2)},
        {"step": "net_charge_pass_threshold", "value": effective_net_charge_min},
        {"step": "organism", "value": organism},
    ]
    if agg_regions:
        derivation.append({
            "step": "aggregation_prone_regions",
            "value": [
                {"start": s, "end": e, "avg_score": sc}
                for s, e, sc in agg_regions
            ],
        })
    else:
        derivation.append({
            "step": "aggregation_prone_regions",
            "value": [],
        })

    # Knowledge gap when no PDB structure
    knowledge_gap: str | None = None
    if pdb_string is None:
        knowledge_gap = (
            "No PDB structure provided; solubility estimated from "
            "intrinsic sequence properties only.  Structural correction "
            "would improve accuracy."
        )

    return TypeCheckResult(
        predicate=f"SolubleExpression(min={min_solubility_score})",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )


# ────────────────────────────────────────────────────────────
# Predicate 2: No Aggregation-Prone Region
# ────────────────────────────────────────────────────────────

def evaluate_no_aggregation_prone_region(
    sequence: str,
    protein: str,
    organism: str,
    max_region_length: int = 5,
    score_threshold: float = -1.0,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check for aggregation-prone regions using CamSol per-residue scoring.

    A region is a consecutive run of ≥6 residues whose smoothed CamSol score
    is below *score_threshold*.  Verdict depends on the weighted worst region,
    with position-weighted aggregation scores (N-terminal regions penalised
    more heavily).

    Membrane proteins (detected heuristically via ≥2 long hydrophobic
    stretches of ≥19 residues) automatically PASS, since they naturally
    contain aggregation-prone transmembrane domains.

    Verdict logic (position-weighted max aggregation score):
    - No qualifying regions (≥6 consecutive)  → PASS
    - max_weighted_score ≤ score_threshold    → PASS (mild)
    - max_weighted_score ≤ 1.0                → LIKELY_PASS (borderline)
    - max_weighted_score ≤ 1.5                → UNCERTAIN
    - max_weighted_score ≤ 2.5                → LIKELY_FAIL
    - max_weighted_score > 2.5                → FAIL

    Only FAIL if the maximum aggregation score > 1.5 (not the previous 1.0).

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (single-letter codes).
        organism: Target organism name.
        max_region_length: Maximum acceptable region length (default 5).
        score_threshold: Per-residue CamSol score threshold (default -1.0).
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict and aggregation-prone region details.
    """
    protein = protein.upper()

    if not protein:
        return TypeCheckResult(
            predicate="NoAggregationProneRegion",
            verdict=Verdict.PASS,
        )

    # --- Membrane protein auto-PASS ---
    if _is_likely_membrane_protein(protein):
        return TypeCheckResult(
            predicate=f"NoAggregationProneRegion(max={max_region_length}, "
                      f"threshold={score_threshold})",
            verdict=Verdict.PASS,
            derivation=[
                {"step": "membrane_protein_detected", "value": True},
                {"step": "aggregation_prone_regions", "value": []},
                {"step": "note", "value": "Auto-PASS: membrane proteins naturally "
                 "contain aggregation-prone transmembrane domains"},
            ],
        )

    # --- Signal peptide detection ---
    # If a signal peptide is predicted, exclude it from aggregation scoring
    # since membrane-associated proteins naturally have N-terminal hydrophobic
    # stretches that are not aggregation-prone in vivo.
    signal_peptide = _detect_signal_peptide(protein)
    sp_info: dict | None = None
    if signal_peptide is not None:
        sp_start, sp_end = signal_peptide
        sp_info = {"start": sp_start, "end": sp_end, "length": sp_end - sp_start}

    agg_regions = _find_aggregation_regions(
        protein,
        score_threshold=score_threshold,
        min_consecutive=_AGG_MIN_CONSECUTIVE_HYDROPHOBIC,
    )

    # Filter out aggregation regions that overlap with the predicted signal peptide
    if signal_peptide is not None:
        sp_start, sp_end = signal_peptide
        filtered_regions: list[tuple[int, int, float]] = []
        for start, end, avg_score in agg_regions:
            # Keep region only if it does NOT overlap with the signal peptide
            if end <= sp_start or start >= sp_end:
                filtered_regions.append((start, end, avg_score))
        agg_regions = filtered_regions

    if not agg_regions:
        derivation_base: list[dict] = [
            {"step": "aggregation_prone_regions", "value": []},
            {"step": "longest_region", "value": 0},
            {"step": "membrane_protein_detected", "value": False},
        ]
        if sp_info is not None:
            derivation_base.append({"step": "signal_peptide_excluded", "value": sp_info})
        return TypeCheckResult(
            predicate=f"NoAggregationProneRegion(max={max_region_length}, "
                      f"threshold={score_threshold})",
            verdict=Verdict.PASS,
            derivation=derivation_base,
        )

    # --- Position-weighted aggregation scores ---
    # N-terminal aggregation is worse than C-terminal.
    # Weight = 1.0 + (1.0 - position_fraction) * penalty_factor
    #   where position_fraction = region_midpoint / protein_length
    n = len(protein)
    worst_weighted_score = 0.0
    worst_region_info: dict = {}

    for start, end, avg_score in agg_regions:
        length = end - start
        midpoint = (start + end) / 2.0
        position_fraction = midpoint / n if n > 0 else 0.5
        # N-terminal (position_fraction near 0) gets higher weight
        position_weight = 1.0 + (1.0 - position_fraction) * (_AGG_NTERM_PENALTY - 1.0)
        # Weighted aggregation score: magnitude of avg_score × length_factor × position_weight
        length_factor = min(length / 10.0, 2.0)  # cap at 2× for very long regions
        weighted_score = abs(avg_score) * length_factor * position_weight

        if weighted_score > worst_weighted_score:
            worst_weighted_score = weighted_score
            worst_region_info = {
                "start": start,
                "end": end,
                "length": length,
                "avg_score": avg_score,
                "position_weight": round(position_weight, 3),
                "weighted_score": round(weighted_score, 3),
            }

    # --- Determine verdict based on weighted max aggregation score ---
    # Only FAIL if max aggregation score > _AGG_FAIL_SCORE_CUTOFF (1.5)
    if worst_weighted_score <= 0.5:
        verdict = Verdict.PASS
        violation = None
    elif worst_weighted_score <= 1.0:
        verdict = Verdict.LIKELY_PASS
        violation = (
            f"Borderline aggregation-prone region: weighted score "
            f"{worst_weighted_score:.3f} (mild)"
        )
    elif worst_weighted_score <= _AGG_FAIL_SCORE_CUTOFF:
        # Moderate aggregation score — LIKELY_FAIL (not UNCERTAIN)
        # since the score provides meaningful evidence
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Aggregation-prone region detected: weighted score "
            f"{worst_weighted_score:.3f} (moderate)"
        )
    elif worst_weighted_score <= 2.5:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Significant aggregation-prone region: weighted score "
            f"{worst_weighted_score:.3f} > {_AGG_FAIL_SCORE_CUTOFF}"
        )
    else:
        verdict = Verdict.FAIL
        violation = (
            f"Very significant aggregation-prone region: weighted score "
            f"{worst_weighted_score:.3f} >> {_AGG_FAIL_SCORE_CUTOFF}"
        )

    # Build derivation: list all aggregation-prone regions with weighting
    derivation: list[dict] = [
        {
            "step": "aggregation_prone_regions",
            "value": [
                {
                    "start": s, "end": e, "length": e - s, "avg_score": sc,
                    "position_weight": round(
                        1.0 + (1.0 - ((s + e) / 2.0) / n) * (_AGG_NTERM_PENALTY - 1.0), 3
                    ) if n > 0 else 1.0,
                }
                for s, e, sc in agg_regions
            ],
        },
        {"step": "worst_weighted_score", "value": round(worst_weighted_score, 3)},
        {"step": "worst_region", "value": worst_region_info},
        {"step": "fail_score_cutoff", "value": _AGG_FAIL_SCORE_CUTOFF},
        {"step": "min_consecutive_hydrophobic", "value": _AGG_MIN_CONSECUTIVE_HYDROPHOBIC},
        {"step": "membrane_protein_detected", "value": False},
    ]
    if sp_info is not None:
        derivation.append({"step": "signal_peptide_excluded", "value": sp_info})

    knowledge_gap: str | None = None
    if pdb_string is None:
        knowledge_gap = (
            "Aggregation-prone regions identified from intrinsic "
            "sequence only.  Structural accessibility correction "
            "would refine the prediction."
        )

    return TypeCheckResult(
        predicate=f"NoAggregationProneRegion(max={max_region_length}, "
                  f"threshold={score_threshold})",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )


# ────────────────────────────────────────────────────────────
# Predicate 3: Charge Composition
# ────────────────────────────────────────────────────────────

def evaluate_charge_composition(
    sequence: str,
    protein: str,
    organism: str,
    min_charged_fraction: float = 0.10,
    max_pI: float = 9.0,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check charge composition for solubility.

    Computes three metrics:
    1. **Charged fraction** — fraction of K, R, H, D, E residues.
       Below *min_charged_fraction* → LIKELY_FAIL (too few charges
       for solubility).
    2. **Isoelectric point (pI)** — pH at which net charge is zero.
       Above *max_pI* → UNCERTAIN (protein may precipitate near its pI
       in typical buffers).
    3. **Net charge / length ratio** — PASS if between -0.15 and +0.15.
       Values outside this range indicate extreme charge imbalance.

    Organism-specific charge preferences:
       - E. coli: prefers slightly negative net charge/length ratio
         (wider negative bound: -0.20)
       - Mammalian: prefers more neutral ratio
         (narrower negative bound: -0.10)

    Uses standard pKa values: Asp=3.90, Glu=4.07, His=6.04, Cys=8.28,
    Tyr=10.07, Lys=10.54, Arg=12.48.

    Both OK → PASS.

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (single-letter codes).
        organism: Target organism name.
        min_charged_fraction: Minimum fraction of charged residues (default 0.10).
        max_pI: Maximum acceptable isoelectric point (default 9.0).
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict, charged fraction, pI, net charge
        ratio, and residue counts in the derivation.
    """
    protein = protein.upper()

    if not protein:
        return TypeCheckResult(
            predicate="ChargeComposition",
            verdict=Verdict.FAIL,
            violation="Empty protein sequence",
        )

    n = len(protein)

    # Count charged residues
    pos_count = sum(1 for aa in protein if aa in {"K", "R", "H"})
    neg_count = sum(1 for aa in protein if aa in {"D", "E"})
    charged_count = pos_count + neg_count
    charged_fraction = charged_count / n

    # Compute isoelectric point
    pI = compute_approximate_pI(protein)

    # Compute net charge / length ratio at pH 7
    net_charge = compute_net_charge(protein, 7.0)
    charge_ratio = net_charge / n if n > 0 else 0.0

    # Organism-specific charge ratio bounds
    org_ratio_lo = _get_org_threshold(organism, "charge_ratio_lo", _CHARGE_RATIO_LO)
    org_ratio_hi = _get_org_threshold(organism, "charge_ratio_hi", _CHARGE_RATIO_HI)

    # Evaluate conditions
    low_charge = charged_fraction < min_charged_fraction
    high_pI = pI > max_pI
    extreme_charge_ratio = charge_ratio < org_ratio_lo or charge_ratio > org_ratio_hi

    if low_charge and high_pI and extreme_charge_ratio:
        verdict = Verdict.FAIL
        violation = (
            f"Low charged fraction ({charged_fraction:.1%} < {min_charged_fraction:.0%}), "
            f"high pI ({pI:.2f} > {max_pI}), and extreme charge ratio "
            f"({charge_ratio:.3f} outside [{org_ratio_lo}, {org_ratio_hi}])"
        )
    elif low_charge and high_pI:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Low charged fraction ({charged_fraction:.1%} < {min_charged_fraction:.0%}) "
            f"and high pI ({pI:.2f} > {max_pI})"
        )
    elif low_charge and extreme_charge_ratio:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Low charged fraction ({charged_fraction:.1%} < {min_charged_fraction:.0%}) "
            f"and extreme charge ratio ({charge_ratio:.3f} outside "
            f"[{org_ratio_lo}, {org_ratio_hi}])"
        )
    elif low_charge:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Low charged fraction ({charged_fraction:.1%} < {min_charged_fraction:.0%}): "
            f"insufficient surface charges for solubility"
        )
    elif high_pI and extreme_charge_ratio:
        # Both high pI and extreme charge ratio → LIKELY_FAIL (not UNCERTAIN)
        # But if charged fraction is high (>20%), the protein has enough
        # charges to be soluble despite the imbalance → just UNCERTAIN
        if charged_fraction > 0.20:
            verdict = Verdict.UNCERTAIN
            violation = (
                f"High pI ({pI:.2f} > {max_pI}) and extreme "
                f"charge ratio ({charge_ratio:.3f}), but high charged fraction "
                f"({charged_fraction:.1%}) suggests solubility is maintained"
            )
        else:
            verdict = Verdict.LIKELY_FAIL
            violation = (
                f"High isoelectric point (pI={pI:.2f} > {max_pI}) and extreme "
                f"charge ratio ({charge_ratio:.3f} outside [{org_ratio_lo}, {org_ratio_hi}])"
            )
    elif high_pI:
        # High isoelectric point alone is ambiguous — some proteins
        # naturally have high pI. Keep as UNCERTAIN.
        verdict = Verdict.UNCERTAIN
        violation = (
            f"High isoelectric point (pI={pI:.2f} > {max_pI}): "
            f"protein may precipitate near its pI in typical buffers"
        )
    elif extreme_charge_ratio:
        # Extreme charge ratio alone is ambiguous — acidic/basic proteins
        # naturally have extreme ratios. If charged fraction is high (>20%),
        # the protein is clearly charged enough to be soluble → PASS.
        if charged_fraction > 0.20:
            verdict = Verdict.PASS
            violation = None
        else:
            verdict = Verdict.UNCERTAIN
            violation = (
                f"Extreme charge ratio ({charge_ratio:.3f} outside "
                f"[{org_ratio_lo}, {org_ratio_hi}]) with low charged fraction "
                f"({charged_fraction:.1%}): charge imbalance may affect solubility"
            )
    else:
        verdict = Verdict.PASS
        violation = None

    derivation: list[dict] = [
        {"step": "charged_fraction", "value": round(charged_fraction, 4)},
        {"step": "min_charged_fraction", "value": min_charged_fraction},
        {"step": "isoelectric_point", "value": round(pI, 2)},
        {"step": "max_pI", "value": max_pI},
        {"step": "net_charge_pH7", "value": round(net_charge, 2)},
        {"step": "charge_ratio", "value": round(charge_ratio, 4)},
        {"step": "charge_ratio_lo", "value": org_ratio_lo},
        {"step": "charge_ratio_hi", "value": org_ratio_hi},
        {"step": "positive_residues", "value": pos_count},
        {"step": "negative_residues", "value": neg_count},
        {"step": "total_charged", "value": charged_count},
        {"step": "protein_length", "value": n},
        {"step": "organism", "value": organism},
    ]

    knowledge_gap: str | None = None
    if pdb_string is None:
        knowledge_gap = (
            "Charge composition assessed from sequence alone.  Surface "
            "accessibility of charged residues (from structure) would "
            "improve solubility prediction."
        )

    return TypeCheckResult(
        predicate=f"ChargeComposition(min_charged={min_charged_fraction}, "
                  f"max_pI={max_pI})",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )


# ────────────────────────────────────────────────────────────
# Predicate 4: No Long Hydrophobic Stretch
# ────────────────────────────────────────────────────────────

def evaluate_no_long_hydrophobic_stretch(
    sequence: str,
    protein: str,
    organism: str,
    max_stretch: int = 7,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check for long consecutive hydrophobic stretches (AILMFWV).

    Stretches longer than *max_stretch* are aggregation-prone.

    Verdict logic (based on the longest stretch found):
    - No stretch > *max_stretch*                     → PASS
    - Longest stretch *max_stretch*+1 to *max_stretch*+3 → LIKELY_PASS
    - Longest stretch *max_stretch*+4 to *max_stretch*+6 → UNCERTAIN
    - Longest stretch > *max_stretch*+6              → FAIL

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (single-letter codes).
        organism: Target organism name.
        max_stretch: Maximum acceptable hydrophobic stretch length (default 7).
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict and positions of long hydrophobic
        stretches in the derivation.
    """
    protein = protein.upper()

    if not protein:
        return TypeCheckResult(
            predicate="NoLongHydrophobicStretch",
            verdict=Verdict.PASS,
        )

    stretches = find_hydrophobic_stretches(protein)

    # Filter to stretches longer than max_stretch
    long_stretches = [(s, e) for s, e in stretches if (e - s) > max_stretch]

    if not long_stretches:
        # All stretches are within the allowed length
        max_found = max((e - s) for s, e in stretches) if stretches else 0
        return TypeCheckResult(
            predicate=f"NoLongHydrophobicStretch(max={max_stretch})",
            verdict=Verdict.PASS,
            derivation=[
                {"step": "max_stretch_found", "value": max_found},
                {"step": "max_stretch_allowed", "value": max_stretch},
                {"step": "long_stretches", "value": []},
            ],
        )

    # Find the longest stretch
    longest = max(e - s for s, e in long_stretches)
    excess = longest - max_stretch

    # Determine verdict based on how much the longest stretch exceeds the limit
    if excess <= _HYDRO_EXCESS_BORDERLINE:
        verdict = Verdict.LIKELY_PASS
        violation = (
            f"Hydrophobic stretch of {longest} residues slightly exceeds "
            f"limit of {max_stretch} (borderline)"
        )
    elif excess <= _HYDRO_EXCESS_UNCERTAIN:
        # Moderate excess — LIKELY_FAIL (not UNCERTAIN) since the
        # quantitative excess provides meaningful evidence
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Hydrophobic stretch of {longest} residues exceeds limit "
            f"of {max_stretch} by {excess}"
        )
    else:
        verdict = Verdict.FAIL
        violation = (
            f"Very long hydrophobic stretch of {longest} residues exceeds "
            f"limit of {max_stretch} by {excess} (aggregation-prone)"
        )

    # Build derivation: list all long hydrophobic stretches
    derivation: list[dict] = [
        {
            "step": "long_stretches",
            "value": [
                {"start": s, "end": e, "length": e - s}
                for s, e in long_stretches
            ],
        },
        {"step": "longest_stretch", "value": longest},
        {"step": "max_stretch_allowed", "value": max_stretch},
    ]

    knowledge_gap: str | None = None
    if pdb_string is None:
        knowledge_gap = (
            "Hydrophobic stretches assessed from sequence alone.  "
            "Structural context (buried vs. exposed) would clarify "
            "whether stretches contribute to aggregation."
        )

    return TypeCheckResult(
        predicate=f"NoLongHydrophobicStretch(max={max_stretch})",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )
