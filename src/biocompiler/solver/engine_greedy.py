"""Greedy Fallback Engine for BioCompiler CSP Solver.

Minimal deterministic greedy codon optimization engine.  Used when neither
OR-Tools nor Z3 is available.  Selects the highest-CAI codon for each amino
acid position, then applies lightweight post-hoc constraint fixes.

This is NOT the main greedy optimizer (which is in ``optimization.py`` and
runs 8+ sequential fix-up passes).  This is a minimal fallback for the
solver module that picks the best codon per-position and applies
organism-aware constraint fixes.

Architecture:
    - For each amino acid, select the synonymous codon with the highest
      codon adaptiveness value for the target organism.
    - Apply lightweight post-hoc constraint fixes:
      * ATTTA motif removal (all organisms)
      * T-run breaking (all organisms)
      * GC content adjustment (all organisms)
      * Cryptic splice site avoidance (eukaryotes only)
      * CpG island avoidance (eukaryotes only)
    - For prokaryotic targets, eukaryote-specific constraints (splice
      sites, CpG islands) are automatically skipped with clear logging.
    - The result always has ``fallback_used=True`` so the caller can
      distinguish it from a proper CSP solution.
"""

from __future__ import annotations

import logging
import math
import random
import time
from typing import Optional

from .types import (
    SolverBackend,
    SolverConfig,
    SolverResult,
    CSPModel,
    SolverBackendProtocol,
)
from ..constants import AA_TO_CODONS, reverse_complement
from ..organisms import CODON_ADAPTIVENESS_TABLES
from ..organism_config import is_eukaryotic_organism

logger = logging.getLogger(__name__)

__all__ = [
    "GreedyEngine",
]

# ── Module-level constants ─────────────────────────────────────────────
_CODON_LENGTH: int = 3
"""Number of nucleotides in a codon."""

_FALLBACK_ORGANISMS: list[str] = ["Homo_sapiens", "Escherichia_coli"]
"""Organisms tried in order when the requested one has no adaptiveness table."""

_UNKNOWN_CODON_PLACEHOLDER: str = "NNN"
"""Placeholder codon when no synonymous codons exist for an amino acid."""

_MAX_ATTTA_ITERATIONS: int = 50
_MAX_T_RUN_ITERATIONS: int = 50
_MAX_GC_ITERATIONS: int = 100
_MAX_SPLICE_ITERATIONS: int = 50
_MAX_CPG_ITERATIONS: int = 50
_T_RUN_THRESHOLD: int = 6


class GreedyEngine(SolverBackendProtocol):
    """Organism-aware greedy codon optimization engine.

    Always available (no external dependencies).  Produces a valid
    translation by selecting the highest-CAI synonymous codon for each
    amino acid position, then applies lightweight post-hoc constraint
    fixes.  Eukaryote-specific constraints (cryptic splice site avoidance,
    CpG island avoidance) are automatically skipped for prokaryotic
    organisms with clear logging.

    Constraints applied for ALL organisms:
      - ATTTA instability motif removal
      - T-run (6+ consecutive T) breaking
      - GC content adjustment into configured range

    Constraints applied for EUKARYOTIC organisms only:
      - Cryptic splice site avoidance (GT/AG dinucleotide disruption)
      - CpG island avoidance (CG dinucleotide disruption)

    This engine exists so that ``CSPSolver._get_engine()`` always has a
    working fallback instead of crashing with ``ImportError`` when
    both OR-Tools and Z3 are unavailable.

    When multiple codons share the same highest CAI score for a given
    amino acid, the tie is broken using a seeded ``random.Random``
    instance for reproducibility.

    Args:
        config: Solver configuration (used for organism info and timeouts).
        seed: Optional deterministic seed for tie-breaking among codons
            with equal CAI scores.  When ``None``, a default seed of 0 is
            used so results are still reproducible.
    """

    def __init__(self, config: SolverConfig, seed: Optional[int] = None) -> None:
        self.config: SolverConfig = config
        self._rng: random.Random = random.Random(seed if seed is not None else 0)

    def solve(self, model: CSPModel) -> SolverResult:
        """Solve the CSP model using greedy per-position codon selection."""
        start_time: float = time.monotonic()
        protein: str = model.protein_sequence.upper()

        # Determine the target organism from the config or model
        target_organism: str = getattr(model.config, "organism", None) or _FALLBACK_ORGANISMS[0]
        if not hasattr(model.config, "organism"):
            logger.debug(
                "Config has no 'organism' attribute; defaulting to %s",
                target_organism,
            )

        adaptiveness: dict[str, float] = CODON_ADAPTIVENESS_TABLES.get(target_organism, {})

        if not adaptiveness:
            for alias in [target_organism] + _FALLBACK_ORGANISMS:
                adaptiveness = CODON_ADAPTIVENESS_TABLES.get(alias, {})
                if adaptiveness:
                    logger.info(
                        "Organism %r not found in adaptiveness tables; using %r instead",
                        target_organism, alias,
                    )
                    target_organism = alias
                    break
            else:
                logger.warning(
                    "No codon adaptiveness table available for %r or any fallback organism; "
                    "CAI will be 0.0", target_organism,
                )

        codons: list[str] = []
        aas: list[str] = list(protein)
        for aa in aas:
            domain: list[str] = AA_TO_CODONS.get(aa, [])
            if not domain:
                logger.warning(
                    "No codons found for amino acid '%s'; using %s placeholder",
                    aa, _UNKNOWN_CODON_PLACEHOLDER,
                )
                codons.append(_UNKNOWN_CODON_PLACEHOLDER)
                continue
            best: str = self._pick_best_codon(domain, adaptiveness)
            codons.append(best)

        sequence: str = "".join(codons)

        # Build sorted codon lists for constraint fix-up (CAI-descending)
        sorted_codons: dict[str, list[str]] = {}
        for aa in set(aas):
            aa_codons = AA_TO_CODONS.get(aa, [])
            sorted_codons[aa] = sorted(
                aa_codons,
                key=lambda c: adaptiveness.get(c, 0.0),
                reverse=True,
            )

        warnings: list[str] = ["Used greedy fallback (lightweight constraint fix-up)"]

        # Determine organism domain for constraint gating
        is_eukaryote: bool = is_eukaryotic_organism(target_organism)

        # ── Constraints applied for ALL organisms ────────────────────
        sequence, attta_warnings = self._fix_attta_motifs(
            sequence, aas, sorted_codons, adaptiveness,
        )
        warnings.extend(attta_warnings)

        sequence, trun_warnings = self._fix_t_runs(
            sequence, aas, sorted_codons, adaptiveness,
        )
        warnings.extend(trun_warnings)

        gc_lo = self.config.gc_lo
        gc_hi = self.config.gc_hi
        sequence, gc_warnings = self._adjust_gc_content(
            sequence, aas, sorted_codons, adaptiveness, gc_lo, gc_hi,
        )
        warnings.extend(gc_warnings)

        # ── Constraints applied for EUKARYOTIC organisms only ─────────
        if is_eukaryote:
            splice_threshold = self.config.cryptic_splice_threshold
            sequence, splice_warnings = self._avoid_splice_sites(
                sequence, aas, sorted_codons, adaptiveness, splice_threshold,
            )
            warnings.extend(splice_warnings)

            if self.config.avoid_cpg:
                sequence, cpg_warnings = self._avoid_cpg_islands(
                    sequence, aas, sorted_codons, adaptiveness,
                )
                warnings.extend(cpg_warnings)
        else:
            logger.info(
                "Organism %r is prokaryotic — skipping eukaryote-specific "
                "constraints (cryptic splice site avoidance, CpG island "
                "avoidance) for compatibility", target_organism,
            )
            warnings.append(
                f"Skipped eukaryote-specific constraints (splice/CpG) "
                f"for prokaryotic organism '{target_organism}'"
            )

        # Compute basic metrics
        gc_val: float = self._compute_gc(sequence)
        cai_val: float = self._compute_cai(sequence, protein, target_organism)

        try:
            result = SolverResult(
                sequence=sequence,
                solved=True,
                backend_used=SolverBackend.GREEDY_FALLBACK,
                protein=protein,
                organism=target_organism,
                cai=cai_val,
                gc_content=gc_val,
                solve_time_seconds=time.monotonic() - start_time,
                fallback_used=True,
                warnings=warnings,
            )
        except Exception:
            logger.exception(
                "Failed to construct SolverResult for protein %r with organism %r",
                protein, target_organism,
            )
            raise

        return result

    # ── Internal helpers ──────────────────────────────────────

    def _pick_best_codon(
        self, domain: list[str], adaptiveness: dict[str, float],
    ) -> str:
        if not domain:
            return _UNKNOWN_CODON_PLACEHOLDER
        scored: list[tuple[str, float]] = [
            (codon, adaptiveness.get(codon, 0.0)) for codon in domain
        ]
        max_score: float = max(s for _, s in scored)
        top_codons: list[str] = [codon for codon, score in scored if score == max_score]
        if len(top_codons) == 1:
            return top_codons[0]
        return self._rng.choice(top_codons)

    # ── Lightweight constraint fix-up methods ────────────────

    def _swap_codon(self, sequence: str, codon_idx: int, new_codon: str) -> str:
        start = codon_idx * _CODON_LENGTH
        return sequence[:start] + new_codon + sequence[start + _CODON_LENGTH:]

    def _fix_attta_motifs(
        self, sequence: str, aas: list[str], sorted_codons: dict[str, list[str]],
        adaptiveness: dict[str, float],
    ) -> tuple[str, list[str]]:
        warnings: list[str] = []
        for _ in range(_MAX_ATTTA_ITERATIONS):
            pos = sequence.find("ATTTA")
            if pos == -1:
                break
            codon_idx = pos // 3
            fixed = False
            for ci in range(max(0, codon_idx - 1), min(len(aas), codon_idx + 2)):
                aa = aas[ci]
                current = sequence[ci * _CODON_LENGTH:ci * _CODON_LENGTH + _CODON_LENGTH]
                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    test = self._swap_codon(sequence, ci, alt)
                    if "ATTTA" not in test:
                        sequence = test
                        fixed = True
                        break
                if fixed:
                    break
            if not fixed:
                warnings.append("ATTTA motif: cannot remove at position %d" % pos)
                break
        return sequence, warnings

    def _fix_t_runs(
        self, sequence: str, aas: list[str], sorted_codons: dict[str, list[str]],
        adaptiveness: dict[str, float],
    ) -> tuple[str, list[str]]:
        warnings: list[str] = []
        for _ in range(_MAX_T_RUN_ITERATIONS):
            max_run, max_pos = 0, -1
            i = 0
            while i < len(sequence):
                if sequence[i] == "T":
                    j = i
                    while j < len(sequence) and sequence[j] == "T":
                        j += 1
                    if j - i > max_run:
                        max_run, max_pos = j - i, i
                    i = j
                else:
                    i += 1
            if max_run < _T_RUN_THRESHOLD:
                break
            codon_idx = (max_pos + max_run // 2) // 3
            if codon_idx >= len(aas):
                break
            aa = aas[codon_idx]
            current = sequence[codon_idx * _CODON_LENGTH:codon_idx * _CODON_LENGTH + _CODON_LENGTH]
            fixed = False
            for alt in sorted_codons.get(aa, []):
                if alt == current:
                    continue
                test = self._swap_codon(sequence, codon_idx, alt)
                if not any(
                    test[k:k + _T_RUN_THRESHOLD] == "T" * _T_RUN_THRESHOLD
                    for k in range(len(test) - _T_RUN_THRESHOLD + 1)
                ):
                    sequence = test
                    fixed = True
                    break
            if not fixed:
                warnings.append("Consecutive T run: cannot fix at position %d" % max_pos)
                break
        return sequence, warnings

    def _adjust_gc_content(
        self, sequence: str, aas: list[str], sorted_codons: dict[str, list[str]],
        adaptiveness: dict[str, float], gc_lo: float, gc_hi: float,
    ) -> tuple[str, list[str]]:
        warnings: list[str] = []
        if not sequence:
            return sequence, warnings
        n_bases = len(sequence)
        gc_count = sum(1 for b in sequence.upper() if b in "GC")
        gc_val = gc_count / n_bases
        if gc_lo <= gc_val <= gc_hi:
            return sequence, warnings
        target = gc_lo if gc_val < gc_lo else gc_hi
        for _ in range(_MAX_GC_ITERATIONS):
            gc_val = gc_count / n_bases
            if gc_lo <= gc_val <= gc_hi:
                break
            best_alt = None
            best_ci = -1
            best_diff = abs(gc_val - target)
            for ci in range(len(aas)):
                aa = aas[ci]
                current = sequence[ci * _CODON_LENGTH:ci * _CODON_LENGTH + _CODON_LENGTH]
                current_gc = sum(1 for b in current if b in "GC")
                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    alt_gc = sum(1 for b in alt if b in "GC")
                    new_gc_count = gc_count - current_gc + alt_gc
                    new_frac = new_gc_count / n_bases
                    diff = abs(new_frac - target)
                    if diff < best_diff:
                        best_diff = diff
                        best_alt = alt
                        best_ci = ci
            if best_alt is None:
                break
            old_codon = sequence[best_ci * _CODON_LENGTH:best_ci * _CODON_LENGTH + _CODON_LENGTH]
            gc_count = gc_count - sum(1 for b in old_codon if b in "GC") + sum(1 for b in best_alt if b in "GC")
            sequence = self._swap_codon(sequence, best_ci, best_alt)
        gc_val = gc_count / n_bases
        if not (gc_lo <= gc_val <= gc_hi):
            warnings.append(
                f"GC adjustment: could not reach [{gc_lo:.2f}, {gc_hi:.2f}], "
                f"current GC={gc_val:.3f}"
            )
        return sequence, warnings

    def _avoid_splice_sites(
        self, sequence: str, aas: list[str], sorted_codons: dict[str, list[str]],
        adaptiveness: dict[str, float], threshold: float,
    ) -> tuple[str, list[str]]:
        """Disrupt GT/AG dinucleotides (eukaryotes only)."""
        warnings: list[str] = []
        if threshold <= 0:
            return sequence, warnings
        for _ in range(_MAX_SPLICE_ITERATIONS):
            fixed_any = False
            for ci in range(len(aas)):
                codon_start = ci * _CODON_LENGTH
                codon_end = codon_start + _CODON_LENGTH
                codon = sequence[codon_start:codon_end]
                if "GT" in codon:
                    aa = aas[ci]
                    gt_free = [c for c in sorted_codons.get(aa, []) if "GT" not in c]
                    if gt_free:
                        sequence = self._swap_codon(sequence, ci, gt_free[0])
                        fixed_any = True
                        break
                if "AG" in codon:
                    aa = aas[ci]
                    ag_free = [c for c in sorted_codons.get(aa, []) if "AG" not in c]
                    if ag_free:
                        sequence = self._swap_codon(sequence, ci, ag_free[0])
                        fixed_any = True
                        break
            if not fixed_any:
                break
        return sequence, warnings

    def _avoid_cpg_islands(
        self, sequence: str, aas: list[str], sorted_codons: dict[str, list[str]],
        adaptiveness: dict[str, float],
    ) -> tuple[str, list[str]]:
        """Disrupt CG dinucleotides (eukaryotes only)."""
        warnings: list[str] = []
        for _ in range(_MAX_CPG_ITERATIONS):
            cpg_positions = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
            if not cpg_positions:
                break
            fixed = False
            for pos in cpg_positions:
                left_ci = pos // 3
                right_ci = (pos + 1) // 3
                for ci in ([left_ci, right_ci] if left_ci != right_ci else [left_ci]):
                    if ci < 0 or ci >= len(aas):
                        continue
                    aa = aas[ci]
                    current = sequence[ci * _CODON_LENGTH:ci * _CODON_LENGTH + _CODON_LENGTH]
                    for alt in sorted_codons.get(aa, []):
                        if alt == current:
                            continue
                        test = self._swap_codon(sequence, ci, alt)
                        if test[pos:pos+2] != "CG":
                            new_cpg = sum(1 for i in range(len(test) - 1) if test[i:i+2] == "CG")
                            old_cpg = sum(1 for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG")
                            if new_cpg < old_cpg:
                                sequence = test
                                fixed = True
                                break
                    if fixed:
                        break
                if fixed:
                    break
            if not fixed:
                break
        return sequence, warnings

    # ── Utility methods ───────────────────────────────────────

    @staticmethod
    def _compute_gc(sequence: str) -> float:
        if not sequence:
            return 0.0
        gc: int = sum(1 for b in sequence.upper() if b in "GC")
        return gc / len(sequence)

    @staticmethod
    def _compute_cai(sequence: str, protein: str, organism: str) -> float:
        adaptiveness: dict[str, float] = CODON_ADAPTIVENESS_TABLES.get(organism, {})
        if not adaptiveness:
            return 0.0
        codons: list[str] = [
            sequence[i : i + _CODON_LENGTH] for i in range(0, len(sequence), _CODON_LENGTH)
        ]
        log_cai_sum: float = 0.0
        n: int = 0
        for codon, aa in zip(codons, protein):
            w: float = adaptiveness.get(codon, 0.0)
            if w > 0:
                log_cai_sum += math.log(w)
                n += 1
        if n == 0:
            return 0.0
        return math.exp(log_cai_sum / n)
