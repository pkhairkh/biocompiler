"""
Performance optimizations for BioCompiler's optimization pipeline.

Target: achieve <10ms per gene for prokaryotic optimization and
<50ms for eukaryotic optimization, approaching DNAchisel's speed
while maintaining BioCompiler's constraint guarantees.

Key strategies:
1. Skip unnecessary constraint passes when constraints are already satisfied
2. Batch constraint checking instead of per-violation checks
3. Fast-path detection for simple optimization problems
4. Lazy initialization of expensive data structures
5. Reduced GC recalculation overhead
6. Warm NUMBA kernel cache

Architecture:
    - should_skip_constraint: Determine if a constraint check can be skipped
    - batch_detect_violations: Detect multiple violation types in a single pass
    - estimate_optimization_complexity: Classify problems as simple/moderate/complex
    - get_fast_path_config: Return optimizer configuration tuned for complexity
    - warm_numba_cache: Pre-compile all NUMBA kernels
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .hybrid_types import HybridResult, Violation

logger = logging.getLogger(__name__)


# ── Complexity thresholds ──────────────────────────────────────────────
_SIMPLE_PROTEIN_LEN = 200       # aa — short proteins
_MODERATE_PROTEIN_LEN = 500    # aa — medium proteins
_SIMPLE_GC_RANGE = 0.30        # wide GC range (0.30-0.70) = simple
_TIGHT_GC_RANGE = 0.10         # tight GC range (<0.10 width) = complex
_SIMPLE_ENZYME_COUNT = 2       # 0-2 enzymes = simple
_MANY_ENZYMES = 5              # 5+ enzymes = complex

# ── CAI thresholds for early termination ───────────────────────────────
_PROKARYOTE_CAI_SKIP_HILLCLIMB = 0.95   # Skip hill climbing if CAI > this
_EUKARYOTE_CAI_REDUCE_HILLCLIMB = 0.90  # Reduce hill climbing if CAI > this

# ── Violation types that are eukaryotic-only ───────────────────────────
_EUKARYOTIC_ONLY_VIOLATIONS = frozenset({
    "cpg_island",
    "cryptic_splice_donor",
    "cryptic_splice_acceptor",
    "avoidable_gt",
})

# ── Organism data cache ────────────────────────────────────────────────
# Avoids re-resolving organism names and reloading tables on every call.
_organism_data_cache: dict[str, dict] = {}


def should_skip_constraint(
    violation_type: str,
    is_prokaryote: bool,
    previous_result: "HybridResult | None" = None,
) -> bool:
    """Determine if a constraint check can be skipped.

    Skips eukaryotic-only constraints for prokaryotes, and skips
    constraints that were already satisfied in a previous optimization
    pass (when previous_result is provided and shows no issues).

    Args:
        violation_type: The constraint type to check (e.g., "cpg_island").
        is_prokaryote: Whether the target organism is prokaryotic.
        previous_result: An optional previous HybridResult to check for
            already-satisfied constraints.

    Returns:
        True if the constraint check should be skipped.
    """
    # Prokaryotes never need eukaryotic constraint checks
    if is_prokaryote and violation_type in _EUKARYOTIC_ONLY_VIOLATIONS:
        return True

    # If we have a previous result with no warnings and high CAI,
    # many constraint re-checks are unnecessary
    if previous_result is not None:
        # If previous result had no violations and high CAI, skip re-check
        if (previous_result.violations_fixed == 0
                and previous_result.cai >= _PROKARYOTE_CAI_SKIP_HILLCLIMB):
            # Soft constraints can be skipped when hard constraints are met
            if violation_type in ("cpg_island", "atttta_motif", "t_run"):
                return True

    return False


def batch_detect_violations(
    sequence: str,
    violation_types: list[str],
    organism_data: dict,
) -> list["Violation"]:
    """Detect multiple violation types in a single pass through the sequence.

    Instead of scanning the sequence separately for each violation type
    (O(n) per type), this function scans the sequence once and checks
    for all requested violation types at each position (O(n) total).

    Currently supports:
    - "restriction_site": Check for restriction enzyme recognition sites
    - "atttta_motif": Check for ATTTA mRNA instability motifs
    - "t_run": Check for 6+ consecutive T nucleotides
    - "gt_dinucleotide": Check for GT dinucleotides (eukaryotic)
    - "ag_dinucleotide": Check for AG dinucleotides (eukaryotic)
    - "cpg_dinucleotide": Check for CpG dinucleotides (eukaryotic)

    Args:
        sequence: DNA sequence string (uppercase ACGT).
        violation_types: List of violation type strings to detect.
        organism_data: Dictionary with organism-specific data, including:
            - "enzymes": List of restriction enzyme names
            - "rs_sites": List of (site, reverse_complement) tuples
            - "is_prokaryote": Whether the organism is prokaryotic

    Returns:
        List of Violation objects for all detected violations.
    """
    from .hybrid_types import Violation, SEVERITY_WEIGHTS

    violations: list[Violation] = []
    n = len(sequence)
    n_codons = n // 3

    # Build set of checks to perform
    check_attta = "atttta_motif" in violation_types
    check_gt = "gt_dinucleotide" in violation_types
    check_ag = "ag_dinucleotide" in violation_types
    check_cpg = "cpg_dinucleotide" in violation_types
    check_trun = "t_run" in violation_types
    check_rs = "restriction_site" in violation_types

    # Skip eukaryotic checks for prokaryotes
    is_prokaryote = organism_data.get("is_prokaryote", False)
    if is_prokaryote:
        check_gt = False
        check_ag = False
        check_cpg = False

    # Pre-process restriction sites
    rs_sites: list[tuple[str, str]] = organism_data.get("rs_sites", [])

    # Single pass through the sequence
    t_run_start = -1
    t_run_len = 0

    for i in range(n):
        b = sequence[i]

        # ── T-run tracking (consecutive T nucleotides) ──
        if check_trun:
            if b == 'T':
                if t_run_len == 0:
                    t_run_start = i
                t_run_len += 1
            else:
                if t_run_len >= 6:
                    codon_indices = list(range(
                        t_run_start // 3,
                        min(n_codons, (t_run_start + t_run_len - 1) // 3 + 1)
                    ))
                    violations.append(Violation(
                        violation_type="t_run",
                        position=t_run_start,
                        severity=SEVERITY_WEIGHTS["t_run"] * (t_run_len - 5),
                        codon_indices=codon_indices,
                        details=f"T-run of {t_run_len} at position {t_run_start}",
                    ))
                t_run_len = 0

        # ── Dinucleotide checks (GT, AG, CG) at position i ──
        if i < n - 1:
            b_next = sequence[i + 1]

            if check_gt and b == 'G' and b_next == 'T':
                codon_idx = i // 3
                next_codon_start = codon_idx * 3 + 3
                is_within = (i + 1) < next_codon_start
                if is_within:
                    involved = [codon_idx]
                else:
                    involved = (
                        [codon_idx, codon_idx + 1]
                        if codon_idx + 1 < n_codons
                        else [codon_idx]
                    )
                violations.append(Violation(
                    violation_type="gt_dinucleotide",
                    position=i,
                    severity=SEVERITY_WEIGHTS["avoidable_gt"],
                    codon_indices=involved,
                    details=f"GT at pos {i} ({'within' if is_within else 'cross'}-codon)",
                ))

            if check_ag and b == 'A' and b_next == 'G':
                codon_idx = i // 3
                next_codon_start = codon_idx * 3 + 3
                is_within = (i + 1) < next_codon_start
                if is_within:
                    involved = [codon_idx]
                else:
                    involved = (
                        [codon_idx, codon_idx + 1]
                        if codon_idx + 1 < n_codons
                        else [codon_idx]
                    )
                violations.append(Violation(
                    violation_type="ag_dinucleotide",
                    position=i,
                    severity=SEVERITY_WEIGHTS["cryptic_splice_acceptor"],
                    codon_indices=involved,
                    details=f"AG at pos {i} ({'within' if is_within else 'cross'}-codon)",
                ))

            if check_cpg and b == 'C' and b_next == 'G':
                codon_idx = i // 3
                violations.append(Violation(
                    violation_type="cpg_dinucleotide",
                    position=i,
                    severity=SEVERITY_WEIGHTS["cpg_island"],
                    codon_indices=[codon_idx],
                    details=f"CpG at pos {i}",
                ))

        # ── ATTTA motif check ──
        if check_attta and b == 'A' and i + 4 < n:
            if sequence[i:i + 5] == "ATTTA":
                codon_indices = list(range(
                    max(0, i // 3 - 1),
                    min(n_codons, (i + 4) // 3 + 2)
                ))
                violations.append(Violation(
                    violation_type="atttta_motif",
                    position=i,
                    severity=SEVERITY_WEIGHTS["atttta_motif"],
                    codon_indices=codon_indices,
                    details=f"ATTTA at position {i}",
                ))

    # Handle T-run at end of sequence
    if check_trun and t_run_len >= 6:
        codon_indices = list(range(
            t_run_start // 3,
            min(n_codons, (t_run_start + t_run_len - 1) // 3 + 1)
        ))
        violations.append(Violation(
            violation_type="t_run",
            position=t_run_start,
            severity=SEVERITY_WEIGHTS["t_run"] * (t_run_len - 5),
            codon_indices=codon_indices,
            details=f"T-run of {t_run_len} at position {t_run_start}",
        ))

    # ── Restriction site check (separate pass for correctness) ──
    if check_rs and rs_sites:
        for site, site_rc in rs_sites:
            pos = 0
            while True:
                p = sequence.find(site, pos)
                if p == -1:
                    break
                codon_indices = list(range(
                    max(0, p // 3),
                    min(n_codons, (p + len(site) - 1) // 3 + 1)
                ))
                violations.append(Violation(
                    violation_type="restriction_site",
                    position=p,
                    severity=SEVERITY_WEIGHTS["restriction_site"],
                    codon_indices=codon_indices,
                    details=f"Site {site} at position {p}",
                ))
                pos = p + 1

            if site_rc:
                pos = 0
                while True:
                    p = sequence.find(site_rc, pos)
                    if p == -1:
                        break
                    # Avoid duplicate if site == site_rc (palindromic)
                    if site_rc == site:
                        break
                    codon_indices = list(range(
                        max(0, p // 3),
                        min(n_codons, (p + len(site_rc) - 1) // 3 + 1)
                    ))
                    violations.append(Violation(
                        violation_type="restriction_site",
                        position=p,
                        severity=SEVERITY_WEIGHTS["restriction_site"],
                        codon_indices=codon_indices,
                        details=f"Site {site_rc} at position {p}",
                    ))
                    pos = p + 1

    return violations


def estimate_optimization_complexity(
    protein: str,
    organism: str,
    enzymes: list[str] | None,
    gc_lo: float,
    gc_hi: float,
) -> str:
    """Classify optimization as "simple", "moderate", or "complex".

    Classification criteria:
    - Simple: short protein, no/few enzymes, wide GC range
    - Moderate: medium protein, some enzymes, moderate GC range
    - Complex: long protein, many enzymes, tight GC range

    Args:
        protein: Amino acid sequence (single-letter codes).
        organism: Target organism name.
        enzymes: List of restriction enzyme names to avoid.
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.

    Returns:
        One of "simple", "moderate", or "complex".
    """
    protein_len = len(protein)
    n_enzymes = len(enzymes) if enzymes else 0
    gc_range = gc_hi - gc_lo

    score = 0

    # Protein length scoring
    if protein_len <= _SIMPLE_PROTEIN_LEN:
        score += 0  # simple
    elif protein_len <= _MODERATE_PROTEIN_LEN:
        score += 1  # moderate
    else:
        score += 2  # complex

    # Enzyme count scoring
    if n_enzymes <= _SIMPLE_ENZYME_COUNT:
        score += 0  # simple
    elif n_enzymes < _MANY_ENZYMES:
        score += 1  # moderate
    else:
        score += 2  # complex

    # GC range scoring
    if gc_range >= _SIMPLE_GC_RANGE * 2:
        score += 0  # very wide range = simple
    elif gc_range >= _SIMPLE_GC_RANGE:
        score += 0  # standard range = simple
    elif gc_range > _TIGHT_GC_RANGE:
        score += 1  # moderate range
    else:
        score += 2  # very tight range = complex

    if score <= 1:
        return "simple"
    elif score <= 3:
        return "moderate"
    else:
        return "complex"


@dataclass
class FastPathConfig:
    """Optimizer configuration tuned for a specific complexity level."""
    max_local_search_iterations: int
    max_hill_climb_iterations: int
    skip_hill_climbing: bool
    hill_climb_passes: int
    micro_optimization_passes: int
    attta_max_iterations: int
    trun_max_iterations: int
    rs_max_iterations: int
    cpg_max_iterations: int
    early_termination: bool
    convergence_window: int


# Pre-built configs for each complexity level
_SIMPLE_PROKARYOTE_CONFIG = FastPathConfig(
    max_local_search_iterations=10,
    max_hill_climb_iterations=3,
    skip_hill_climbing=True,      # Skip if CAI already > 0.95
    hill_climb_passes=1,
    micro_optimization_passes=2,
    attta_max_iterations=20,
    trun_max_iterations=20,
    rs_max_iterations=30,
    cpg_max_iterations=0,         # No CpG for prokaryotes
    early_termination=True,
    convergence_window=2,
)

_MODERATE_PROKARYOTE_CONFIG = FastPathConfig(
    max_local_search_iterations=20,
    max_hill_climb_iterations=5,
    skip_hill_climbing=True,
    hill_climb_passes=2,
    micro_optimization_passes=3,
    attta_max_iterations=50,
    trun_max_iterations=50,
    rs_max_iterations=50,
    cpg_max_iterations=0,
    early_termination=True,
    convergence_window=3,
)

_COMPLEX_PROKARYOTE_CONFIG = FastPathConfig(
    max_local_search_iterations=30,
    max_hill_climb_iterations=5,
    skip_hill_climbing=False,
    hill_climb_passes=3,
    micro_optimization_passes=5,
    attta_max_iterations=100,
    trun_max_iterations=100,
    rs_max_iterations=100,
    cpg_max_iterations=0,
    early_termination=True,
    convergence_window=3,
)

_SIMPLE_EUKARYOTE_CONFIG = FastPathConfig(
    max_local_search_iterations=15,
    max_hill_climb_iterations=3,
    skip_hill_climbing=False,
    hill_climb_passes=2,
    micro_optimization_passes=2,
    attta_max_iterations=30,
    trun_max_iterations=30,
    rs_max_iterations=30,
    cpg_max_iterations=15,
    early_termination=True,
    convergence_window=2,
)

_MODERATE_EUKARYOTE_CONFIG = FastPathConfig(
    max_local_search_iterations=25,
    max_hill_climb_iterations=5,
    skip_hill_climbing=False,
    hill_climb_passes=3,
    micro_optimization_passes=3,
    attta_max_iterations=50,
    trun_max_iterations=50,
    rs_max_iterations=50,
    cpg_max_iterations=20,
    early_termination=True,
    convergence_window=3,
)

_COMPLEX_EUKARYOTE_CONFIG = FastPathConfig(
    max_local_search_iterations=50,
    max_hill_climb_iterations=10,
    skip_hill_climbing=False,
    hill_climb_passes=5,
    micro_optimization_passes=5,
    attta_max_iterations=100,
    trun_max_iterations=100,
    rs_max_iterations=100,
    cpg_max_iterations=30,
    early_termination=True,
    convergence_window=3,
)


def get_fast_path_config(complexity: str, is_prokaryote: bool) -> FastPathConfig:
    """Return optimizer configuration tuned for the estimated complexity level.

    Simple optimizations should skip hill climbing, use fewer iterations,
    and enable early termination. Complex optimizations use more iterations
    and disable early skip heuristics.

    Args:
        complexity: One of "simple", "moderate", or "complex" (from
            estimate_optimization_complexity).
        is_prokaryote: Whether the target organism is prokaryotic.

    Returns:
        FastPathConfig with appropriate iteration limits and flags.
    """
    if is_prokaryote:
        if complexity == "simple":
            return _SIMPLE_PROKARYOTE_CONFIG
        elif complexity == "moderate":
            return _MODERATE_PROKARYOTE_CONFIG
        else:
            return _COMPLEX_PROKARYOTE_CONFIG
    else:
        if complexity == "simple":
            return _SIMPLE_EUKARYOTE_CONFIG
        elif complexity == "moderate":
            return _MODERATE_EUKARYOTE_CONFIG
        else:
            return _COMPLEX_EUKARYOTE_CONFIG


def warm_numba_cache() -> None:
    """Pre-compile all NUMBA kernels so first optimization call isn't slow.

    This function triggers NUMBA JIT compilation for all kernels used
    during optimization. It should be called once at startup or during
    application initialization to avoid ~1-2s latency on the first
    optimization call.

    If NUMBA is not available, this function is a no-op.
    """
    try:
        from .numba_kernels import (
            HAS_NUMBA,
            count_gc,
            count_dinucleotides,
            compute_cai_kernel,
            scan_restriction_sites,
            find_all_dinucleotide_positions,
            compute_cai_incremental,
            batch_codon_swap_score,
            fast_gc_window,
            fast_dinucleotide_count,
            count_gc_parallel,
            scan_restriction_sites_multi,
            seq_to_bytes,
        )

        if not HAS_NUMBA:
            logger.debug("NUMBA not available; warm_numba_cache is a no-op")
            return

        import numpy as np

        # Create tiny dummy inputs for compilation
        dummy_seq = np.array([65, 84, 71, 67], dtype=np.uint8)  # b'ATGC'
        dummy_dinuc = np.array([71, 84], dtype=np.uint8)  # b'GT'
        dummy_adapt = np.array([0.5, 1.0], dtype=np.float64)
        dummy_indices = np.array([0, 1], dtype=np.int64)

        # Warm up each kernel
        try:
            count_gc(dummy_seq)
            count_dinucleotides(dummy_seq, dummy_dinuc)
            compute_cai_kernel(dummy_adapt, dummy_indices, 2)
            scan_restriction_sites(dummy_seq, dummy_dinuc, 2)
            find_all_dinucleotide_positions(dummy_seq, dummy_dinuc)
            compute_cai_incremental(-1.3862943611198906, 2, 0.5, 1.0)
            batch_codon_swap_score(
                dummy_adapt, dummy_indices, 2, 0,
                np.array([0, 1], dtype=np.int64), 2, -1.3862943611198906
            )
            fast_gc_window(dummy_seq, 2)
            fast_dinucleotide_count(
                dummy_seq,
                np.array([[71, 84], [67, 71]], dtype=np.uint8), 2
            )
            count_gc_parallel(dummy_seq)
            scan_restriction_sites_multi(
                dummy_seq, dummy_dinuc,
                np.array([0], dtype=np.int64),
                np.array([2], dtype=np.int64), 1
            )
            # Also warm seq_to_bytes conversion
            seq_to_bytes("ATGC")
        except Exception as e:
            logger.debug("NUMBA warmup failed (non-fatal): %s", e)

        logger.debug("NUMBA kernel cache warmed successfully")

    except ImportError:
        logger.debug("NUMBA kernels not importable; warm_numba_cache is a no-op")


def get_organism_data(organism: str) -> dict:
    """Get cached organism data, resolving the organism name only once.

    This avoids repeated calls to resolve_organism() and
    CODON_ADAPTIVENESS_TABLES lookups which can be expensive when
    optimizing many genes in sequence.

    Args:
        organism: Organism name (canonical or alias).

    Returns:
        Dictionary with organism data including:
        - "species_cai": Codon adaptiveness table
        - "gc_target": Target GC fraction
        - "is_prokaryote": Whether the organism is prokaryotic
    """
    if organism in _organism_data_cache:
        return _organism_data_cache[organism]

    try:
        from ..organisms import (
            CODON_ADAPTIVENESS_TABLES,
            ORGANISM_GC_TARGETS,
            resolve_organism,
        )
        from ..organism_config import is_eukaryotic_organism

        canonical = resolve_organism(organism, strict=False)
        species_cai = dict(CODON_ADAPTIVENESS_TABLES.get(
            canonical, CODON_ADAPTIVENESS_TABLES.get("Escherichia_coli", {})
        ))
        gc_target = ORGANISM_GC_TARGETS.get(canonical, 0.50)
        is_prok = not is_eukaryotic_organism(canonical)

        data = {
            "species_cai": species_cai,
            "gc_target": gc_target,
            "is_prokaryote": is_prok,
            "canonical_name": canonical,
        }
    except ImportError:
        data = {
            "species_cai": {},
            "gc_target": 0.50,
            "is_prokaryote": False,
            "canonical_name": organism,
        }

    _organism_data_cache[organism] = data
    return data


def clear_caches() -> None:
    """Clear all performance caches.

    This should be called when organism configurations change or
    between test runs to avoid stale cache data.
    """
    _organism_data_cache.clear()


def should_skip_mrna_stability(sequence: str) -> bool:
    """Determine if mRNA stability pass can be skipped.

    The mRNA stability pass is unnecessary when there are no ATTTA
    motifs or other destabilizing motifs in the sequence.

    Args:
        sequence: Optimized DNA sequence.

    Returns:
        True if the mRNA stability pass can be skipped.
    """
    # Quick check for the most common destabilizing motif
    if "ATTTA" not in sequence:
        return True
    return False


def should_skip_cpg_elimination(is_prokaryote: bool, sequence: str) -> bool:
    """Determine if CpG elimination pass can be skipped.

    CpG elimination is never needed for prokaryotes and can be
    skipped for eukaryotic sequences that have no CpG dinucleotides.

    Args:
        is_prokaryote: Whether the target organism is prokaryotic.
        sequence: Optimized DNA sequence.

    Returns:
        True if the CpG elimination pass can be skipped.
    """
    if is_prokaryote:
        return True
    # Quick scan for CpG dinucleotides
    if "CG" not in sequence:
        return True
    return False


def should_skip_utr_suggestions(include_utr: bool) -> bool:
    """Determine if UTR suggestion pass can be skipped.

    Args:
        include_utr: Whether UTR suggestions were requested.

    Returns:
        True if UTR suggestions should be skipped.
    """
    return not include_utr
