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
import time
from typing import Optional

from .types import (
    SolverBackend,
    SolverConfig,
    SolverResult,
    CSPModel,
)
from ..constants import AA_TO_CODONS
from ..organisms import CODON_ADAPTIVENESS_TABLES

logger = logging.getLogger(__name__)


class GreedyEngine:
    """Minimal greedy codon optimization engine.

    Always available (no external dependencies).  Produces a valid
    translation by selecting the highest-CAI synonymous codon for each
    amino acid position.  Does **not** enforce constraints like GC range,
    restriction site avoidance, or splice site avoidance.

    This engine exists so that ``CSPSolver._get_engine()`` always has a
    working fallback instead of crashing with ``ImportError`` when
    both OR-Tools and Z3 are unavailable.

    Args:
        config: Solver configuration (used for organism info and timeouts).
    """

    def __init__(self, config: SolverConfig) -> None:
        self.config = config

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
        start_time = time.monotonic()
        protein = model.protein_sequence.upper()
        organism = model.config.backend.value if hasattr(model.config, "organism") else "Homo_sapiens"

        # Try to get the organism from the config or model
        target_organism = getattr(model.config, "_organism", None) or "Homo_sapiens"

        adaptiveness = CODON_ADAPTIVENESS_TABLES.get(target_organism, {})

        if not adaptiveness:
            # Try common aliases
            for alias in [target_organism, "Homo_sapiens", "E_coli"]:
                adaptiveness = CODON_ADAPTIVENESS_TABLES.get(alias, {})
                if adaptiveness:
                    target_organism = alias
                    break

        codons = []
        for aa in protein:
            domain = AA_TO_CODONS.get(aa, [])
            if not domain:
                logger.warning("No codons found for amino acid '%s'; using NNN placeholder", aa)
                codons.append("NNN")
                continue
            # Pick highest-CAI codon for this amino acid
            best = max(domain, key=lambda c: adaptiveness.get(c, 0.0))
            codons.append(best)

        sequence = "".join(codons)

        # Compute basic metrics
        gc_val = self._compute_gc(sequence)
        cai_val = self._compute_cai(sequence, protein, target_organism)

        return SolverResult(
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

    # ── Utility methods ───────────────────────────────────────

    @staticmethod
    def _compute_gc(sequence: str) -> float:
        """Compute GC content fraction of a DNA sequence."""
        if not sequence:
            return 0.0
        gc = sum(1 for b in sequence.upper() if b in "GC")
        return gc / len(sequence)

    @staticmethod
    def _compute_cai(sequence: str, protein: str, organism: str) -> float:
        """Compute Codon Adaptation Index for a sequence.

        Uses the organism's codon adaptiveness table.  Returns 0.0 if
        the organism table is unavailable.
        """
        import math

        adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})
        if not adaptiveness:
            return 0.0

        codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
        log_cai_sum = 0.0
        n = 0
        for codon, aa in zip(codons, protein):
            w = adaptiveness.get(codon, 0.0)
            if w > 0:
                log_cai_sum += math.log(w)
                n += 1

        if n == 0:
            return 0.0
        return math.exp(log_cai_sum / n)
