"""
BioCompiler ViennaRNA Wrapper — mRNA Secondary Structure Prediction
====================================================================

Wraps the ViennaRNA library for thermodynamic folding of mRNA sequences,
replacing the toy hairpin model in ``type_system.py`` with real
nearest-neighbor free-energy calculations.

Expected ViennaRNA version: ≥ 2.5.0
(https://www.tbi.univie.ac.at/RNA/)

The module is **always importable**.  When ViennaRNA is not installed,
public functions return sentinel results with ``success=False`` and
``method="unavailable"`` so callers can degrade gracefully.

Two backends are attempted in order:

1. **Python bindings** — ``import RNA`` (fastest, most features).
2. **CLI fallback** — subprocess call to ``RNAfold --noPS``
   (broader compatibility, MFE only).

References:
    Lorenz et al., Monatshefte für Chemie 2011; 142:345–349 (ViennaRNA 2.0)
    Lorenz et al., Bioinformatics 2011; 27:1827–1828 (RNAfold / RNA.pf_fold)
"""

from __future__ import annotations

import logging
import re
import subprocess
from dataclasses import dataclass, field
from typing import Literal, TypedDict

logger = logging.getLogger(__name__)

__all__ = [
    # Data classes
    "StemLoop", "MFEResult", "AccessibilityResult", "MRNAStructureResult",
    # Core prediction functions
    "is_viennarna_available", "predict_mfe", "predict_accessibility",
    "find_stable_structures", "compute_5prime_dg",
    "check_mrna_structure_viennarna", "predict_mfe_batch",
    "predict_accessibility_batch",
    # Organism / heuristic helpers (Issue 1 & 3)
    "get_organism_dg_threshold", "estimate_dg_from_gc",
    "find_most_stable_region",
    # Public constants needed by callers
    "DEFAULT_DG_THRESHOLD", "CDS_DG_THRESHOLD", "ORGANISM_DG_THRESHOLDS",
    "REGION_FULL", "REGION_5UTR", "REGION_START_CODON", "REGION_CDS",
    "SLIDING_WINDOW_SIZE", "SLIDING_WINDOW_STEP",
    "DEFAULT_5PRIME_WINDOW",
    "DEFAULT_FOLD_TIMEOUT_SECONDS",
    "DEFAULT_FULL_LENGTH_CUTOFF",
    "DEFAULT_OVERLAP_THRESHOLD",
    "DEFAULT_STEP",
    "DEFAULT_WINDOW_SIZE",
    "EXPECTED_VIENNARNA_VERSION",
    "MAX_FOLD_TIMEOUT_SECONDS",
    "NEAREST_NEIGHBOR_AU",
    "NEAREST_NEIGHBOR_GC",
    "NEAREST_NEIGHBOR_GU",
]

# ── Constants ──────────────────────────────────────────────

EXPECTED_VIENNARNA_VERSION = (2, 5, 0)
DEFAULT_DG_THRESHOLD: float = -15.0
DEFAULT_WINDOW_SIZE: int = 80
DEFAULT_STEP: int = 20
DEFAULT_5PRIME_WINDOW: int = 50
REGION_FULL = "full"
REGION_5UTR = "5utr"
REGION_START_CODON = "start_codon"
REGION_CDS = "cds"

# Named constants replacing raw magic numbers
DEFAULT_FULL_LENGTH_CUTOFF: int = 100
DEFAULT_FOLD_TIMEOUT_SECONDS: int = 5
MAX_FOLD_TIMEOUT_SECONDS: int = 30
DEFAULT_OVERLAP_THRESHOLD: float = 0.5

# Keep in sync with solver/constraints.py NEAREST_NEIGHBOR_* constants
NEAREST_NEIGHBOR_GC: float = -1.5
NEAREST_NEIGHBOR_AU: float = -0.5
NEAREST_NEIGHBOR_GU: float = -0.3

# ── Issue 1: Organism-specific ΔG thresholds ──────────────
# Different organisms have different tolerances for mRNA secondary
# structure.  Prokaryotes (especially E. coli) are sensitive to
# structure at the RBS; eukaryotes tolerate more structure in the
# 5'UTR but have a different threshold.  Coding regions (CDS)
# normally have secondary structure and should use a much more
# relaxed threshold.

ORGANISM_DG_THRESHOLDS: dict[str, float] = {
    # E. coli: flag if ΔG < -15 kcal/mol in RBS region (first 50 nt)
    "E_coli": -15.0,
    "E_coli_K12": -15.0,
    "E_coli_BL21": -15.0,
    # Human: flag if ΔG < -25 kcal/mol in 5'UTR
    "Homo_sapiens": -25.0,
    "human": -25.0,
    # Other eukaryotes: moderate threshold
    "Saccharomyces_cerevisiae": -20.0,
    "yeast": -20.0,
    "Mus_musculus": -25.0,
    "mouse": -25.0,
    "CHO_K1": -25.0,
    "cho": -25.0,
}

# CDS regions normally contain secondary structure; use a much more
# relaxed (more negative) threshold so that normal coding-sequence
# structure is not flagged as problematic.
CDS_DG_THRESHOLD: float = -50.0

# ── Issue 2: Region-aware check parameters ────────────────
# Prokaryotes: only check RBS/start codon region (first ~50 nt)
PROKARYOTE_CHECK_WINDOW: int = 50
# Eukaryotes: check 5'UTR region (typically first ~50 nt of UTR)
EUKARYOTE_CHECK_WINDOW: int = 50

# ── Issue 3: GC-content heuristic fallback parameters ─────
GC_STABLE_THRESHOLD: float = 0.60   # GC% > 60% → likely stable
GC_UNSTABLE_THRESHOLD: float = 0.40  # GC% < 40% → likely unstable
# Rough ΔG estimate per nucleotide for GC-rich sequences (kcal/mol/nt)
GC_DG_PER_NT: float = -0.4
# Rough ΔG estimate per nucleotide for AT-rich sequences (kcal/mol/nt)
AT_DG_PER_NT: float = -0.1

# ── Issue 4: Sliding window parameters for batch prediction ─
SLIDING_WINDOW_SIZE: int = 50
SLIDING_WINDOW_STEP: int = 10

# ── Data classes ───────────────────────────────────────────

@dataclass
class StemLoop:
    """A single stem-loop (hairpin) region in an RNA secondary structure.

    Attributes:
        start: 0-based start position in the original full sequence.
        end:   0-based exclusive end position.
        structure: Dot-bracket notation for this region only.
        mfe:   Minimum free energy (kcal/mol). More negative = more stable.
    """
    start: int
    end: int
    structure: str
    mfe: float

    def __repr__(self) -> str:
        return (f"StemLoop(start={self.start}, end={self.end}, "
                f"mfe={self.mfe:.1f}, structure={self.structure!r})")


@dataclass
class MFEResult:
    """Result of a minimum free-energy (MFE) folding calculation.

    Attributes:
        structure: Dot-bracket notation for the folded sequence.
        mfe:       Minimum free energy in kcal/mol.
        sequence:  The RNA sequence that was folded (U not T).
        base_pairing_probs: Position → P(paired). Populated only when
                            partition-function computation succeeds.
        stem_loops: Identified stem-loops with ΔG below threshold.
        success:    Whether the computation succeeded.
        method:     Backend used — ``"viennarna_python"``,
                    ``"viennarna_cli"``, ``"trivial"``, or ``"unavailable"``.
        error:      Error message if *success* is False; None otherwise.
    """
    structure: str = ""
    mfe: float = 0.0
    sequence: str = ""
    base_pairing_probs: dict[int, float] = field(default_factory=dict)
    stem_loops: list[StemLoop] = field(default_factory=list)
    success: bool = False
    method: str = "unavailable"
    error: str | None = None


@dataclass
class AccessibilityResult:
    """Per-region RNA accessibility (fraction unpaired) from the
    partition function.

    Attributes:
        region:              Region label (``"5utr"``, ``"start_codon"``,
                             ``"full"``, or custom).
        mean_accessibility:  Average P(unpaired) across all positions (0–1).
        position_accessibility: Position → P(unpaired).
        success:             Whether the computation succeeded.
        method:              Backend used.
        error:               Error message if *success* is False.
    """
    region: str = ""
    mean_accessibility: float = 0.0
    position_accessibility: dict[int, float] = field(default_factory=dict)
    success: bool = False
    method: str = "unavailable"
    error: str | None = None


# ── Internal helpers ───────────────────────────────────────

def _dna_to_rna(dna: str) -> str:
    """Convert DNA to RNA (T → U)."""
    return dna.upper().replace("T", "U")


def _extract_region(dna: str, region: str) -> str:
    """Extract a functional region from a DNA sequence."""
    dna = dna.upper()
    if region == REGION_FULL:
        return dna
    if region == REGION_5UTR:
        return dna[:DEFAULT_5PRIME_WINDOW]
    if region == REGION_START_CODON:
        return dna[:DEFAULT_FULL_LENGTH_CUTOFF]
    logger.debug("Unknown region %r, using full sequence", region)
    return dna


def _identify_stem_loops(structure: str) -> list[tuple[int, int, str]]:
    """Identify stem-loop regions from a dot-bracket string.

    Returns list of (start, end, sub_structure) tuples.
    """
    n = len(structure)
    visited: set[int] = set()
    loops: list[tuple[int, int, str]] = []
    i = 0
    while i < n:
        if structure[i] == "(" and i not in visited:
            stem_start = i
            j = i
            while j < n and structure[j] == "(":
                j += 1
            stem_open_end = j
            # Find loop (dots between arms)
            loop_end = stem_open_end
            while loop_end < n and structure[loop_end] == ".":
                loop_end += 1
            # Find closing arm
            close_count = 0
            k = loop_end
            while k < n and structure[k] == ")":
                close_count += 1
                k += 1
            open_count = stem_open_end - stem_start
            if open_count >= 1 and close_count >= 1 and loop_end > stem_open_end:
                paired = min(open_count, close_count)
                actual_end = min(stem_start + paired + (loop_end - stem_open_end) + paired, n)
                sub = structure[stem_start:actual_end]
                if sub and sub[0] == "(" and sub[-1] == ")":
                    loops.append((stem_start, actual_end, sub))
                    for p in range(stem_start, actual_end):
                        visited.add(p)
                    i = actual_end
                    continue
        i += 1
    return loops


def _parse_version(version_str: str) -> tuple[int, int, int] | None:
    """Parse a version string like '2.5.1' into (2, 5, 1)."""
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", version_str)
    if m:
        return (int(m.group(1)), int(m.group(2)),
                int(m.group(3)) if m.group(3) else 0)
    return None


# ── Availability check ─────────────────────────────────────

def is_viennarna_available() -> bool:
    """Check if ViennaRNA is available for RNA folding.

    Returns True if either the ViennaRNA Python bindings (``import RNA``)
    can be imported, or the ``RNAfold`` CLI is found on ``$PATH``.
    """
    try:
        import RNA  # noqa: F401
        return True
    except ImportError:
        logger.debug("ViennaRNA Python bindings not available, trying CLI")
    try:
        r = subprocess.run(["RNAfold", "--version"],
                           capture_output=True, text=True, timeout=DEFAULT_FOLD_TIMEOUT_SECONDS)
        return r.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return False


def _get_viennarna_version() -> tuple[int, int, int] | None:
    """Return ViennaRNA version as (major, minor, patch) or None."""
    try:
        import RNA
        v = getattr(RNA, "__version__", "")
        if not v and hasattr(RNA, "fold_version"):
            v = RNA.fold_version()
        return _parse_version(v)
    except ImportError:
        logger.debug("ViennaRNA Python bindings not available for version check, trying CLI")
    try:
        r = subprocess.run(["RNAfold", "--version"],
                           capture_output=True, text=True, timeout=DEFAULT_FOLD_TIMEOUT_SECONDS)
        if r.returncode == 0:
            return _parse_version(r.stdout + r.stderr)
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        logger.debug("RNAfold CLI not found or timed out for version check")
    return None


# ── Backend: ViennaRNA Python bindings ─────────────────────

def _predict_mfe_python(rna_sequence: str) -> tuple[str, float] | None:
    """Fold using ViennaRNA Python bindings. Returns (structure, mfe) or None."""
    try:
        import RNA
    except ImportError:
        return None
    try:
        fc = RNA.fold_compound(rna_sequence)
        structure, mfe = fc.mfe()
        return (structure, mfe)
    except Exception as exc:
        logger.warning("ViennaRNA RNA.fold() failed: %s", exc)
        return None


def _predict_pf_python(rna_sequence: str) -> tuple[str, float, dict[int, float]] | None:
    """Compute partition function and base-pairing probabilities.

    Returns (structure, ensemble_energy, pairing_probs) or None.
    """
    try:
        import RNA
    except ImportError:
        return None
    try:
        fc = RNA.fold_compound(rna_sequence)
        structure, energy = fc.pf()
        pairing_probs: dict[int, float] = {}
        n = len(rna_sequence)
        for i in range(1, n + 1):  # ViennaRNA uses 1-based indexing
            prob_paired = 0.0
            for j in range(i + 1, n + 1):
                p = fc.bp_prob(i, j)
                if p > 0.0:
                    prob_paired += p
            pairing_probs[i - 1] = min(prob_paired, 1.0)
        return (structure, energy, pairing_probs)
    except Exception as exc:
        logger.warning("ViennaRNA partition function failed: %s", exc)
        return None


# ── Backend: RNAfold CLI fallback ──────────────────────────

_RNAFOLD_RE = re.compile(
    r"^([.(){}<>]+)\s+\(\s*([-]?\d+\.\d+)\s*\)", re.MULTILINE
)


def _predict_mfe_cli(rna_sequence: str) -> tuple[str, float] | None:
    """Fold using ``RNAfold --noPS`` CLI. Returns (structure, mfe) or None."""
    try:
        proc = subprocess.run(
            ["RNAfold", "--noPS"], input=rna_sequence,
            capture_output=True, text=True, timeout=MAX_FOLD_TIMEOUT_SECONDS,
        )
    except FileNotFoundError:
        return None
    except (subprocess.TimeoutExpired, OSError) as exc:
        logger.warning("RNAfold CLI failed: %s", exc)
        return None
    if proc.returncode != 0:
        logger.warning("RNAfold exit %d: %s", proc.returncode, proc.stderr[:200])
        return None
    lines = [ln.strip() for ln in proc.stdout.strip().splitlines() if ln.strip()]
    if len(lines) < 2:
        return None
    m = _RNAFOLD_RE.search(lines[1])
    if not m:
        logger.warning("Cannot parse RNAfold output: %r", lines[1][:200])
        return None
    return (m.group(1), float(m.group(2)))


# ── Unified folding with automatic backend selection ───────

def _fold_mfe(rna_sequence: str) -> MFEResult:
    """Fold an RNA sequence, trying Python bindings then CLI fallback."""
    # Strategy 1: Python bindings
    result = _predict_mfe_python(rna_sequence)
    if result is not None:
        structure, mfe = result
        pairing_probs: dict[int, float] = {}
        pf_result = _predict_pf_python(rna_sequence)
        if pf_result is not None:
            pairing_probs = pf_result[2]
        return MFEResult(structure=structure, mfe=mfe, sequence=rna_sequence,
                         base_pairing_probs=pairing_probs, success=True,
                         method="viennarna_python")

    # Strategy 2: CLI fallback
    result = _predict_mfe_cli(rna_sequence)
    if result is not None:
        structure, mfe = result
        return MFEResult(structure=structure, mfe=mfe, sequence=rna_sequence,
                         success=True, method="viennarna_cli")

    # Strategy 3: Unavailable
    return MFEResult(sequence=rna_sequence, success=False, method="unavailable",
                     error="ViennaRNA not available: Python bindings not importable "
                           "and RNAfold CLI not found")


# ── Public API: predict_mfe ────────────────────────────────

def predict_mfe(
    dna_sequence: str,
    region: str = REGION_FULL,
    dg_threshold: float = DEFAULT_DG_THRESHOLD,
) -> MFEResult:
    """Predict the minimum free-energy secondary structure of an mRNA.

    Converts DNA → RNA, extracts the requested region, folds using
    ViennaRNA, and identifies stem-loops with ΔG below *dg_threshold*.

    Region mapping:
      - ``"full"``         — fold the entire sequence
      - ``"5utr"``         — fold the first 50 nt
      - ``"start_codon"``  — fold the first 100 nt
      - Any other string   — fold the full sequence

    Args:
        dna_sequence: DNA sequence (T, not U).
        region:       Functional region to analyse (default ``"full"``).
        dg_threshold: ΔG threshold for stem-loop identification
                      (default -15.0 kcal/mol).

    Returns:
        :class:`MFEResult` with structure, ΔG, pairing probabilities,
        and identified stem-loops.
    """
    if not dna_sequence:
        return MFEResult(success=False, error="Empty DNA sequence provided")

    rna = _dna_to_rna(_extract_region(dna_sequence, region))

    if len(rna) < 4:
        return MFEResult(structure="." * len(rna), mfe=0.0, sequence=rna,
                         success=True, method="trivial")

    result = _fold_mfe(rna)
    if not result.success:
        return result

    # Identify stem-loops by re-folding each detected hairpin region
    for start, end, _ in _identify_stem_loops(result.structure):
        sub_rna = rna[start:end]
        if len(sub_rna) < 4:
            continue
        sub_result = _fold_mfe(sub_rna)
        if sub_result.success and sub_result.mfe < dg_threshold:
            result.stem_loops.append(
                StemLoop(start=start, end=end,
                         structure=sub_result.structure, mfe=sub_result.mfe))

    return result


# ── Public API: predict_accessibility ──────────────────────

def predict_accessibility(
    dna_sequence: str,
    region: str = REGION_5UTR,
) -> AccessibilityResult:
    """Compute RNA accessibility (fraction unpaired) for a region.

    Uses the ViennaRNA partition function to compute per-position
    P(unpaired) = 1 − P(paired).  The 5′ UTR region (first ~50 nt)
    is the default target as it is most relevant for ribosome binding.

    Args:
        dna_sequence: DNA sequence (T, not U).
        region:       Functional region (default ``"5utr"``).

    Returns:
        :class:`AccessibilityResult` with mean and per-position accessibility.
    """
    if not dna_sequence:
        return AccessibilityResult(region=region, success=False,
                                   error="Empty DNA sequence provided")

    rna = _dna_to_rna(_extract_region(dna_sequence, region))

    if len(rna) < 4:
        pos_acc = {i: 1.0 for i in range(len(rna))}
        return AccessibilityResult(region=region, mean_accessibility=1.0,
                                   position_accessibility=pos_acc,
                                   success=True, method="trivial")

    # Try Python bindings (partition function)
    pf_result = _predict_pf_python(rna)
    if pf_result is not None:
        pos_acc = {i: max(0.0, min(1.0, 1.0 - pf_result[2].get(i, 0.0)))
                   for i in range(len(rna))}
        mean = sum(pos_acc.values()) / len(pos_acc) if pos_acc else 0.0
        return AccessibilityResult(region=region, mean_accessibility=mean,
                                   position_accessibility=pos_acc,
                                   success=True, method="viennarna_python")

    # CLI fallback: binary approximation from MFE structure
    cli_result = _predict_mfe_cli(rna)
    if cli_result is not None:
        structure, _ = cli_result
        pos_acc = {i: (1.0 if ch == "." else 0.0)
                   for i, ch in enumerate(structure)}
        mean = sum(pos_acc.values()) / len(pos_acc) if pos_acc else 0.0
        return AccessibilityResult(region=region, mean_accessibility=mean,
                                   position_accessibility=pos_acc,
                                   success=True, method="viennarna_cli")

    return AccessibilityResult(region=region, success=False, method="unavailable",
                               error="ViennaRNA not available: cannot compute accessibility")


# ── Public API: find_stable_structures ─────────────────────

def find_stable_structures(
    dna_sequence: str,
    dg_threshold: float = DEFAULT_DG_THRESHOLD,
    window_size: int = DEFAULT_WINDOW_SIZE,
    step: int = DEFAULT_STEP,
) -> list[StemLoop]:
    """Scan an mRNA for locally stable secondary structures.

    Slides a window of *window_size* nucleotides across the sequence,
    folding each window independently.  Returns stem-loops whose ΔG is
    below (more negative than) *dg_threshold*.  Overlapping hits are
    merged (≥ 50% overlap), keeping the most stable.

    Args:
        dna_sequence: DNA sequence (T, not U).
        dg_threshold: ΔG cutoff in kcal/mol (default -15.0).
        window_size:  Sliding window size in nt (default 80).
        step:         Step size for the sliding window in nt (default 20).

    Returns:
        List of :class:`StemLoop` objects sorted by start position.
    """
    if not dna_sequence:
        return []

    rna = _dna_to_rna(dna_sequence.upper())
    n = len(rna)

    if n < window_size:
        result = _fold_mfe(rna)
        if result.success and result.mfe < dg_threshold:
            return [StemLoop(start=0, end=n,
                             structure=result.structure, mfe=result.mfe)]
        return []

    candidates: list[StemLoop] = []
    for start in range(0, n - window_size + 1, step):
        result = _fold_mfe(rna[start:start + window_size])
        if result.success and result.mfe < dg_threshold:
            candidates.append(StemLoop(start=start, end=start + window_size,
                                       structure=result.structure, mfe=result.mfe))

    # Merge overlapping candidates (≥ 50% overlap → keep more stable)
    if not candidates:
        return []
    candidates.sort(key=lambda s: s.start)
    merged: list[StemLoop] = [candidates[0]]
    for cur in candidates[1:]:
        prev = merged[-1]
        overlap = max(0, min(prev.end, cur.end) - max(prev.start, cur.start))
        shorter = min(prev.end - prev.start, cur.end - cur.start)
        if shorter > 0 and overlap / shorter >= DEFAULT_OVERLAP_THRESHOLD:
            if cur.mfe < prev.mfe:
                merged[-1] = cur
        else:
            merged.append(cur)
    return merged


# ── Public API: compute_5prime_dg ──────────────────────────

def compute_5prime_dg(
    dna_sequence: str,
    window: int = DEFAULT_5PRIME_WINDOW,
) -> float:
    """Quick ΔG computation for the 5′ region of an mRNA.

    Folds the first *window* nucleotides and returns the MFE.  This is
    the region most relevant for ribosome binding; the CSP solver uses
    this value as a soft constraint to penalize sequences with very
    stable 5′ structures.

    Args:
        dna_sequence: DNA sequence (T, not U).
        window:       Number of nucleotides from the 5′ end to fold (default 50).

    Returns:
        ΔG in kcal/mol.  Returns 0.0 if ViennaRNA is not available.
    """
    if not dna_sequence:
        return 0.0
    rna = _dna_to_rna(dna_sequence.upper()[:window])
    if len(rna) < 4:
        return 0.0
    result = _fold_mfe(rna)
    return result.mfe if result.success else 0.0


# ── Integration helper for type_system.py ──────────────────

class MRNAStructureResult(TypedDict):
    """Return type for check_mrna_structure_viennarna.

    Attributes:
        dg: Free energy of the folded structure (kcal/mol).
            More negative = more stable. 0.0 for trivial/unavailable.
        method: Backend used — ``"viennarna_python"``, ``"viennarna_cli"``,
            ``"gc_heuristic_fallback"``, ``"toy_hairpin_fallback"``, or ``"trivial"``.
        structure: Dot-bracket notation for the folded region.
            Empty string when ViennaRNA is unavailable.
        viennarna_used: True when the result came from ViennaRNA
            (Python bindings or CLI); False for toy/trivial fallback.
        has_hairpin: True when at least one stem-loop with ΔG below
            *dg_threshold* was detected in the analysis window.
        hairpin_positions: List of ``(start, end)`` tuples for each
            detected hairpin, 0-based within the *window_start:window_end*
            subsequence.  Empty when no hairpins are found.
        region: The functional region that was checked (``"rbs"``,
            ``"5utr"``, ``"cds"``, or ``"full"``).
        organism: The organism key used for threshold selection.
        dg_threshold: The ΔG threshold that was applied.
        most_stable_window_start: Start position of the most stable
            sliding window (0-based in the original sequence).
            None when no sliding window scan was performed.
    """
    dg: float
    method: str
    structure: str
    viennarna_used: bool
    has_hairpin: bool
    hairpin_positions: list[tuple[int, int]]
    region: str
    organism: str
    dg_threshold: float
    most_stable_window_start: int | None


def get_organism_dg_threshold(
    organism: str,
    region: str = REGION_FULL,
) -> float:
    """Return the organism-specific ΔG threshold for flagging mRNA structure.

    Organism-specific thresholds:
      - E. coli: only flag if ΔG < -15 kcal/mol in RBS region
      - Human: only flag if ΔG < -25 kcal/mol in 5'UTR region
      - CDS (coding) regions: use a much more relaxed threshold (-50)
        because secondary structure in the CDS is normal and expected.

    Args:
        organism: Organism key (e.g. ``"E_coli_K12"``, ``"Homo_sapiens"``).
        region: Functional region being checked (``"5utr"``, ``"start_codon"``,
                ``"cds"``, ``"full"``).

    Returns:
        ΔG threshold in kcal/mol.
    """
    # CDS regions normally contain secondary structure — much more relaxed
    if region == REGION_CDS:
        return CDS_DG_THRESHOLD

    # Look up organism-specific threshold
    threshold = ORGANISM_DG_THRESHOLDS.get(organism)
    if threshold is not None:
        return threshold

    # Try to resolve via organism_config for unknown keys
    try:
        from ..organisms.config import get_organism_config
        config = get_organism_config(organism)
        if not config.is_eukaryote:
            return -15.0  # prokaryote default
        return -25.0  # eukaryote default
    except Exception:
        pass

    # Ultimate fallback: use the default threshold
    return DEFAULT_DG_THRESHOLD


def _determine_check_region(
    organism: str,
) -> tuple[str, int]:
    """Determine which region to check based on organism type.

    Issue 2: Region awareness — not all mRNA secondary structure is bad.
    Only check the RBS/start codon region for prokaryotes, and the 5'UTR
    for eukaryotes.  Structure in the middle of the CDS is expected and
    not problematic.

    Args:
        organism: Organism key.

    Returns:
        (region_label, window_size) tuple.
    """
    try:
        from ..organisms.config import get_organism_config
        config = get_organism_config(organism)
        if not config.is_eukaryote:
            # Prokaryote: check RBS/start codon region (first ~50 nt)
            return ("rbs", PROKARYOTE_CHECK_WINDOW)
        else:
            # Eukaryote: check 5'UTR region
            return (REGION_5UTR, EUKARYOTE_CHECK_WINDOW)
    except Exception:
        pass

    # Fallback: use a simple heuristic based on the organism key
    org_lower = organism.lower()
    if any(k in org_lower for k in ("ecoli", "e_coli", "escherichia")):
        return ("rbs", PROKARYOTE_CHECK_WINDOW)
    return (REGION_5UTR, EUKARYOTE_CHECK_WINDOW)


def estimate_dg_from_gc(
    dna_sequence: str,
    region: str = REGION_FULL,
) -> float:
    """Estimate ΔG from GC content when ViennaRNA is unavailable.

    .. warning::
       **APPROXIMATE** — This is a rough heuristic, NOT a thermodynamic
       calculation.  It assumes a simple linear relationship between GC
       content and folding free energy, which is a major simplification.
       Actual ΔG depends on nearest-neighbour stacking, loop entropy,
       and other factors that GC content alone cannot capture.  Use only
       as a last resort when ViennaRNA is not installed.

    Issue 3: Fallback heuristic.  Instead of returning UNCERTAIN when
    ViennaRNA is not installed, estimate ΔG from the GC content and
    sequence length:
      - GC% > 60% → likely stable structure (more negative ΔG)
      - GC% < 40% → likely unstable structure (ΔG close to 0)
      - 40% ≤ GC% ≤ 60% → intermediate estimate

    The estimate uses a linear interpolation between AT-rich and
    GC-rich per-nucleotide ΔG contributions, scaled by sequence length.
    This is approximate and should not be relied upon for precise
    thermodynamic predictions.

    Args:
        dna_sequence: DNA sequence (T, not U).
        region: Functional region to analyse (default ``"full"``).

    Returns:
        Estimated ΔG in kcal/mol.  More negative = more stable.
        **This value is approximate and may differ significantly from
        the true thermodynamic ΔG.**
    """
    dna = _extract_region(dna_sequence.upper(), region)
    if not dna:
        return 0.0

    n = len(dna)
    if n < 4:
        return 0.0

    gc_count = dna.count("G") + dna.count("C")
    gc_frac = gc_count / n

    if gc_frac >= GC_STABLE_THRESHOLD:
        # GC-rich → likely stable structure
        dg_per_nt = GC_DG_PER_NT
    elif gc_frac <= GC_UNSTABLE_THRESHOLD:
        # AT-rich → likely unstable
        dg_per_nt = AT_DG_PER_NT
    else:
        # Linear interpolation between unstable and stable
        t = (gc_frac - GC_UNSTABLE_THRESHOLD) / (
            GC_STABLE_THRESHOLD - GC_UNSTABLE_THRESHOLD
        )
        dg_per_nt = AT_DG_PER_NT + t * (GC_DG_PER_NT - AT_DG_PER_NT)

    # Scale by window length; only the region that can form a hairpin
    # contributes.  A rough estimate: half the window can pair.
    effective_pairing_length = n // 2
    estimated_dg = dg_per_nt * effective_pairing_length

    # For CDS regions, apply a relaxation factor since structure in
    # coding regions is normal and expected
    if region == REGION_CDS:
        estimated_dg *= 0.5  # Halve the estimated stability for CDS

    return estimated_dg


def find_most_stable_region(
    dna_sequence: str,
    window_size: int = SLIDING_WINDOW_SIZE,
    step: int = SLIDING_WINDOW_STEP,
) -> tuple[int, float, str]:
    """Find the most stable region using a sliding window approach.

    Issue 4: Batch prediction.  For long sequences, fold the entire
    mRNA using a sliding window and extract regional ΔG.  Uses
    50nt windows with 10nt step to find the most stable region,
    and only reports that.

    For sequences shorter than *window_size* (default 50 nt), the
    full sequence is folded as a single window rather than applying
    a sliding window.

    Args:
        dna_sequence: DNA sequence (T, not U).
        window_size: Sliding window size in nt (default 50).
        step:        Step size in nt (default 10).

    Returns:
        (start_position, mfe, method) tuple for the most stable window.
        If ViennaRNA is unavailable, falls back to GC heuristic.
        Returns (0, 0.0, "trivial") for empty/too-short sequences.
    """
    rna = _dna_to_rna(dna_sequence.upper())
    n = len(rna)

    if n < 4:
        return (0, 0.0, "trivial")

    best_start = 0
    best_dg = 0.0
    best_method = "trivial"

    # If sequence is shorter than the sliding window, fold the full
    # sequence as a single window (Issue 2: handle sequences < 50nt).
    if n < window_size:
        result = _fold_mfe(rna)
        if result.success:
            return (0, result.mfe, result.method)
        # ViennaRNA unavailable — use GC heuristic on full sequence
        est = estimate_dg_from_gc(dna_sequence)
        return (0, est, "gc_heuristic_fallback")

    # Sequence fits in exactly one window — fold it
    if n == window_size:
        result = _fold_mfe(rna)
        if result.success:
            return (0, result.mfe, result.method)
        est = estimate_dg_from_gc(dna_sequence)
        return (0, est, "gc_heuristic_fallback")

    # Slide the window across the sequence (n > window_size guaranteed here)
    for start in range(0, n - window_size + 1, step):
        sub_rna = rna[start:start + window_size]
        result = _fold_mfe(sub_rna)
        if result.success:
            if result.mfe < best_dg:
                best_dg = result.mfe
                best_start = start
                best_method = result.method
        else:
            # ViennaRNA unavailable for this window — use GC heuristic
            est = estimate_dg_from_gc(dna_sequence[start:start + window_size])
            if est < best_dg:
                best_dg = est
                best_start = start
                best_method = "gc_heuristic_fallback"

    # If all windows failed or nothing stable found, try GC heuristic on full seq
    if best_method == "trivial":
        est = estimate_dg_from_gc(dna_sequence)
        return (0, est, "gc_heuristic_fallback")

    return (best_start, best_dg, best_method)


def check_mrna_structure_viennarna(
    dna_sequence: str,
    window_start: int = 0,
    window_end: int = 50,
    dg_threshold: float = -15.0,
    organism: str = "",
    region: str = "",
) -> MRNAStructureResult:
    """ViennaRNA-backed mRNA secondary structure check for predicate 11.

    Drop-in replacement for the toy hairpin model in
    ``type_system.check_mrna_secondary_structure``.  When ViennaRNA is
    not available, uses a GC-content heuristic instead of returning
    UNCERTAIN.

    Issue 1 (Threshold calibration): Uses organism-specific thresholds:
      - E. coli: only flag if ΔG < -15 kcal/mol in RBS region
      - Human: only flag if ΔG < -25 kcal/mol in 5'UTR region
      - CDS regions: much more relaxed threshold (-50 kcal/mol)

    Issue 2 (Region awareness): Only checks the RBS/start codon region
    for prokaryotes and the 5'UTR for eukaryotes.  Structure in the
    middle of the CDS is expected and not problematic.

    Issue 3 (GC heuristic fallback): When ViennaRNA is not installed,
    estimates ΔG from GC content instead of returning UNCERTAIN.

    Issue 4 (Batch/sliding window): For long sequences, uses a sliding
    window approach (50nt windows, 10nt step) to find the most stable
    region and only reports that.

    Args:
        dna_sequence: DNA sequence (uppercase, T not U).
        window_start: Start position of the analysis window (default 0).
        window_end:   End position of the analysis window (default 50).
        dg_threshold: ΔG threshold for FAIL verdict (default -15.0).
                      Overridden by organism-specific threshold when
                      *organism* is provided.
        organism:     Organism key (e.g. ``"E_coli_K12"``, ``"Homo_sapiens"``).
                      When provided, enables organism-specific thresholds
                      and region-aware checking.
        region:       Functional region (``"5utr"``, ``"start_codon"``,
                      ``"cds"``, ``"full"``).  When empty, auto-detected
                      from *organism* if provided.

    Returns:
        :class:`MRNAStructureResult` with keys ``dg``, ``method``,
        ``structure``, ``viennarna_used``, ``has_hairpin``,
        ``hairpin_positions``, ``region``, ``organism``,
        ``dg_threshold``, ``most_stable_window_start``.
    """
    # ── Issue 2: Region awareness ────────────────────────────
    # Determine which region to check based on organism
    check_region = region
    check_window = window_end - window_start  # caller-specified window size
    if not check_region and organism:
        check_region, check_window = _determine_check_region(organism)
    elif not check_region:
        check_region = REGION_START_CODON  # backward-compatible default

    # ── Issue 1: Organism-specific ΔG threshold ──────────────
    # Use organism-specific threshold if organism is provided.
    # When the caller supplied default window bounds (0, 50) and an
    # organism is specified, also apply the organism-specific check
    # window size so the correct subsequence is analysed.
    effective_threshold = dg_threshold
    if organism:
        org_threshold = get_organism_dg_threshold(organism, check_region)
        effective_threshold = org_threshold
        # Apply organism-specific check window when caller uses defaults
        if window_start == 0 and window_end == 50:
            window_end = check_window

    # Determine the effective window to check
    # For CDS regions, we still scan but with a much more relaxed threshold
    if check_region == REGION_CDS:
        # CDS: use sliding window but relaxed threshold
        effective_threshold = CDS_DG_THRESHOLD

    # Extract the subsequence to analyse
    rna = _dna_to_rna(dna_sequence[window_start:window_end])
    if len(rna) < 4:
        return {"dg": 0.0, "method": "trivial",
                "structure": "." * len(rna), "viennarna_used": False,
                "has_hairpin": False, "hairpin_positions": [],
                "region": check_region, "organism": organism,
                "dg_threshold": effective_threshold,
                "most_stable_window_start": None}

    # ── Issue 4: Sliding window batch prediction ─────────────
    # For long sequences, use sliding window to find most stable region.
    # Issue 2: For sequences shorter than SLIDING_WINDOW_SIZE, fold the
    # full sequence directly instead of applying a sliding window.
    most_stable_start: int | None = None
    if len(rna) > SLIDING_WINDOW_SIZE:
        best_start, best_dg, best_method = find_most_stable_region(
            dna_sequence[window_start:window_end],
            window_size=SLIDING_WINDOW_SIZE,
            step=SLIDING_WINDOW_STEP,
        )
        most_stable_start = best_start

        # Re-fold the most stable window to get full details.
        # Use min() to handle sequences where the best window is at the
        # end and is shorter than SLIDING_WINDOW_SIZE.
        best_end = min(best_start + SLIDING_WINDOW_SIZE, len(rna))
        best_rna = rna[best_start:best_end]
        if len(best_rna) >= 4:
            result = _fold_mfe(best_rna)
            if result.success:
                hairpin_positions: list[tuple[int, int]] = []
                for s, e, _ in _identify_stem_loops(result.structure):
                    sub_rna = best_rna[s:e]
                    if len(sub_rna) < 4:
                        continue
                    sub_result = _fold_mfe(sub_rna)
                    if sub_result.success and sub_result.mfe < effective_threshold:
                        hairpin_positions.append((s, e))
                return {"dg": result.mfe, "method": result.method,
                        "structure": result.structure, "viennarna_used": True,
                        "has_hairpin": len(hairpin_positions) > 0,
                        "hairpin_positions": hairpin_positions,
                        "region": check_region, "organism": organism,
                        "dg_threshold": effective_threshold,
                        "most_stable_window_start": most_stable_start}
            # ── Issue 3: GC heuristic fallback ─────────────────
            # ViennaRNA unavailable for this window — use GC heuristic
            est_dg = estimate_dg_from_gc(
                dna_sequence[window_start + best_start:window_start + best_start + SLIDING_WINDOW_SIZE],
                region=check_region,
            )
            return {"dg": est_dg, "method": "gc_heuristic_fallback",
                    "structure": "", "viennarna_used": False,
                    "has_hairpin": est_dg < effective_threshold,
                    "hairpin_positions": [],
                    "region": check_region, "organism": organism,
                    "dg_threshold": effective_threshold,
                    "most_stable_window_start": most_stable_start,
                    "verdict_suggestion": "UNCERTAIN"}

    # ── Standard path: fold the window directly ─────────────
    result = _fold_mfe(rna)
    if result.success:
        # Detect hairpins in the folded result
        hairpin_positions = []
        for start, end, _ in _identify_stem_loops(result.structure):
            sub_rna = rna[start:end]
            if len(sub_rna) < 4:
                continue
            sub_result = _fold_mfe(sub_rna)
            if sub_result.success and sub_result.mfe < effective_threshold:
                hairpin_positions.append((start, end))
        return {"dg": result.mfe, "method": result.method,
                "structure": result.structure, "viennarna_used": True,
                "has_hairpin": len(hairpin_positions) > 0,
                "hairpin_positions": hairpin_positions,
                "region": check_region, "organism": organism,
                "dg_threshold": effective_threshold,
                "most_stable_window_start": most_stable_start}

    # ── Issue 3: GC heuristic fallback when ViennaRNA unavailable ──
    # Instead of returning UNCERTAIN or the toy hairpin model,
    # estimate ΔG from GC content.
    est_dg = estimate_dg_from_gc(dna_sequence[window_start:window_end], region=check_region)
    return {"dg": est_dg, "method": "gc_heuristic_fallback",
            "structure": "", "viennarna_used": False,
            "has_hairpin": est_dg < effective_threshold,
            "hairpin_positions": [],
            "region": check_region, "organism": organism,
            "dg_threshold": effective_threshold,
            "most_stable_window_start": most_stable_start,
            "verdict_suggestion": "UNCERTAIN"}


# ── Batch processing ───────────────────────────────────────

def predict_mfe_batch(
    dna_sequences: list[str],
    region: str = REGION_FULL,
    dg_threshold: float = DEFAULT_DG_THRESHOLD,
) -> list[MFEResult]:
    """Fold multiple DNA sequences in batch.  Individual failures are isolated."""
    results: list[MFEResult] = []
    for seq in dna_sequences:
        try:
            results.append(predict_mfe(seq, region=region, dg_threshold=dg_threshold))
        except Exception as exc:
            results.append(MFEResult(sequence=_dna_to_rna(seq) if seq else "",
                                     success=False, error=f"Batch fold failed: {exc}"))
    return results


def predict_accessibility_batch(
    dna_sequences: list[str],
    region: str = REGION_5UTR,
) -> list[AccessibilityResult]:
    """Compute accessibility for multiple sequences in batch."""
    results: list[AccessibilityResult] = []
    for seq in dna_sequences:
        try:
            results.append(predict_accessibility(seq, region=region))
        except Exception as exc:
            results.append(AccessibilityResult(region=region, success=False,
                                                error=f"Batch accessibility failed: {exc}"))
    return results


# ── Module-level version check (informational) ────────────

def _log_version_info() -> None:
    """Log ViennaRNA version info at module import time."""
    version = _get_viennarna_version()
    if version is None:
        logger.info("ViennaRNA not available — mRNA folding will use fallback")
        return
    v_str = ".".join(str(v) for v in version)
    exp_str = ".".join(str(v) for v in EXPECTED_VIENNARNA_VERSION)
    if version >= EXPECTED_VIENNARNA_VERSION:
        logger.info("ViennaRNA %s detected (≥ %s)", v_str, exp_str)
    else:
        logger.warning("ViennaRNA %s detected but ≥ %s expected", v_str, exp_str)

try:
    _log_version_info()
except Exception:
    logger.warning("ViennaRNA operation failed", exc_info=True)
