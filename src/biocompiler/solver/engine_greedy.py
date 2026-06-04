"""Greedy Fallback Engine for BioCompiler CSP Solver.

Minimal deterministic greedy codon optimization engine.  Used when neither
OR-Tools nor Z3 is available.  Selects the highest-CAI codon for each amino
acid position without global constraint satisfaction.

This is NOT the main greedy optimizer (which is in ``optimization.py`` and
runs 8+ sequential fix-up passes).  This is a minimal fallback for the
solver module that simply picks the best codon per-position.

Architecture:
    - For each amino acid, select the synonymous codon with the highest
      codon adaptiveness value for the target organism.
    - No constraint checking or post-hoc repair — just a fast baseline.
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
from ..constants import AA_TO_CODONS
from ..organisms import CODON_ADAPTIVENESS_TABLES

logger = logging.getLogger(__name__)

__all__ = [
    "GreedyEngine",
]

# ── Module-level constants ─────────────────────────────────────────────
_CODON_LENGTH: int = 3
"""Number of nucleotides in a codon."""

_FALLBACK_ORGANISMS: list[str] = ["Homo_sapiens", "E_coli"]
"""Organisms tried in order when the requested one has no adaptiveness table."""

_UNKNOWN_CODON_PLACEHOLDER: str = "NNN"
"""Placeholder codon when no synonymous codons exist for an amino acid."""


class GreedyEngine(SolverBackendProtocol):
    """Minimal greedy codon optimization engine.

    Always available (no external dependencies).  Produces a valid
    translation by selecting the highest-CAI synonymous codon for each
    amino acid position.  Does **not** enforce constraints like GC range,
    restriction site avoidance, or splice site avoidance.

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
        """Initialize the greedy engine with solver configuration.

        Args:
            config: Solver configuration (used for organism info and timeouts).
            seed: Optional deterministic seed for tie-breaking among codons
                with equal CAI scores.  When ``None``, a default seed of 0 is
                used so results are still reproducible.
        """
        self.config: SolverConfig = config
        self._rng: random.Random = random.Random(seed if seed is not None else 0)

    def solve(self, model: CSPModel) -> SolverResult:
        """Solve the CSP model using greedy per-position codon selection.

        Parameters
        ----------
        model : CSPModel
            The constraint model (used for protein sequence and organism).

        Returns
        -------
        SolverResult
            A solved result with ``fallback_used=True`` and
            ``backend_used=SolverBackend.GREEDY_FALLBACK``.
        """
        start_time: float = time.monotonic()
        protein: str = model.protein_sequence.upper()

        # Determine the target organism from the config or model
        target_organism: str = getattr(model.config, "_organism", None) or _FALLBACK_ORGANISMS[0]
        if not hasattr(model.config, "_organism"):
            logger.debug(
                "Config has no '_organism' attribute; defaulting to %s",
                target_organism,
            )

        adaptiveness: dict[str, float] = CODON_ADAPTIVENESS_TABLES.get(target_organism, {})

        if not adaptiveness:
            # Try common aliases
            for alias in [target_organism] + _FALLBACK_ORGANISMS:
                adaptiveness = CODON_ADAPTIVENESS_TABLES.get(alias, {})
                if adaptiveness:
                    logger.info(
                        "Organism %r not found in adaptiveness tables; using %r instead",
                        target_organism,
                        alias,
                    )
                    target_organism = alias
                    break
            else:
                logger.warning(
                    "No codon adaptiveness table available for %r or any fallback organism; "
                    "CAI will be 0.0",
                    target_organism,
                )

        codons: list[str] = []
        for aa in protein:
            domain: list[str] = AA_TO_CODONS.get(aa, [])
            if not domain:
                logger.warning(
                    "No codons found for amino acid '%s'; using %s placeholder",
                    aa,
                    _UNKNOWN_CODON_PLACEHOLDER,
                )
                codons.append(_UNKNOWN_CODON_PLACEHOLDER)
                continue
            # Pick highest-CAI codon for this amino acid, with seeded tie-breaking
            best: str = self._pick_best_codon(domain, adaptiveness)
            codons.append(best)

        sequence: str = "".join(codons)

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
                warnings=["Used greedy fallback (no constraint satisfaction)"],
            )
        except Exception:
            logger.exception(
                "Failed to construct SolverResult for protein %r with organism %r",
                protein,
                target_organism,
            )
            raise

        return result

    # ── Internal helpers ──────────────────────────────────────

    def _pick_best_codon(
        self,
        domain: list[str],
        adaptiveness: dict[str, float],
    ) -> str:
        """Select the highest-CAI codon from *domain*, breaking ties with the RNG.

        When multiple codons share the same maximum adaptiveness value, the
        seeded random number generator deterministically picks one, ensuring
        reproducibility across runs with the same seed.
        """
        if not domain:
            return _UNKNOWN_CODON_PLACEHOLDER

        # Compute CAI scores for all candidates
        scored: list[tuple[str, float]] = [
            (codon, adaptiveness.get(codon, 0.0)) for codon in domain
        ]
        max_score: float = max(s for _, s in scored)
        top_codons: list[str] = [codon for codon, score in scored if score == max_score]

        if len(top_codons) == 1:
            return top_codons[0]

        # Tie-break using seeded RNG for reproducibility
        return self._rng.choice(top_codons)

    # ── Utility methods ───────────────────────────────────────

    @staticmethod
    def _compute_gc(sequence: str) -> float:
        """Compute GC content fraction of a DNA sequence."""
        if not sequence:
            return 0.0
        gc: int = sum(1 for b in sequence.upper() if b in "GC")
        return gc / len(sequence)

    @staticmethod
    def _compute_cai(sequence: str, protein: str, organism: str) -> float:
        """Compute Codon Adaptation Index for a sequence.

        Uses the organism's codon adaptiveness table.  Returns 0.0 if
        the organism table is unavailable.
        """
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
