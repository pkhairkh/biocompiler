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

logger = logging.getLogger(__name__)

__all__ = [
    "StemLoop", "MFEResult", "AccessibilityResult",
    "is_viennarna_available", "predict_mfe", "predict_accessibility",
    "find_stable_structures", "compute_5prime_dg",
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
        return dna[:50]
    if region == REGION_START_CODON:
        return dna[:100]
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
        pass
    try:
        r = subprocess.run(["RNAfold", "--version"],
                           capture_output=True, text=True, timeout=5)
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
        pass
    try:
        r = subprocess.run(["RNAfold", "--version"],
                           capture_output=True, text=True, timeout=5)
        if r.returncode == 0:
            return _parse_version(r.stdout + r.stderr)
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        pass
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
            capture_output=True, text=True, timeout=30,
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
        if shorter > 0 and overlap / shorter >= 0.5:
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

def check_mrna_structure_viennarna(
    dna_sequence: str,
    window_start: int = 0,
    window_end: int = 50,
    dg_threshold: float = -15.0,
) -> dict:
    """ViennaRNA-backed mRNA secondary structure check for predicate 11.

    Drop-in replacement for the toy hairpin model in
    ``type_system.check_mrna_secondary_structure``.  Falls back to the
    toy model if ViennaRNA is not available.

    Args:
        dna_sequence: DNA sequence (uppercase, T not U).
        window_start: Start position of the analysis window (default 0).
        window_end:   End position of the analysis window (default 50).
        dg_threshold: ΔG threshold for FAIL verdict (default -15.0).

    Returns:
        Dict with keys ``dg``, ``method``, ``structure``, ``viennarna_used``.
    """
    rna = _dna_to_rna(dna_sequence[window_start:window_end])
    if len(rna) < 4:
        return {"dg": 0.0, "method": "trivial",
                "structure": "." * len(rna), "viennarna_used": False}

    result = _fold_mfe(rna)
    if result.success:
        return {"dg": result.mfe, "method": result.method,
                "structure": result.structure, "viennarna_used": True}

    # Fallback: toy hairpin model (mirrors type_system.py logic)
    gc = au = gu = 0
    half = len(rna) // 2
    first, second = rna[:half], rna[half:2 * half]
    for i in range(min(len(first), len(second))):
        b5, b3 = first[i], second[len(second) - 1 - i]
        if (b5 + b3) in ("GC", "CG"):
            gc += 1
        elif (b5 + b3) in ("AU", "UA"):
            au += 1
        elif (b5 + b3) in ("GU", "UG"):
            gu += 1
    return {"dg": -1.5 * gc - 0.5 * au - 0.3 * gu,
            "method": "toy_hairpin_fallback", "structure": "",
            "viennarna_used": False}


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
    pass
