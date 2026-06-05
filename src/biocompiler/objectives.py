"""
BioCompiler Custom Objective Functions
=======================================

User-defined objective functions for gene optimization.

The default optimizer maximizes Codon Adaptation Index (CAI).  This module
provides a protocol and built-in objective functions that allow users to
optimize for other criteria (GC content, codon pair bias, or a combination).

Architecture:

- ``ObjectiveFunction`` is a ``Protocol`` that any callable satisfying the
  signature ``(dna: str, protein: str, organism: str) -> float`` conforms to.
  Higher return values are better.
- Built-in objectives (``cai_objective``, ``cai_gc_balanced_objective``,
  ``codon_pair_objective``, ``min_max_gc_objective``) cover the most common
  use cases.
- The ``resolve_objective()`` function accepts either a callable or a string
  name (``"cai"``, ``"cai_gc_balanced"``, ``"codon_pair"``, ``"min_max_gc"``)
  and returns a callable conforming to the protocol.
- The optimizer's codon selection uses the custom objective for *scoring*
  alternatives: instead of always picking the highest-CAI codon, it picks
  the codon that maximises the objective when substituted into the current
  sequence context.

Usage::

    from biocompiler.objectives import cai_gc_balanced_objective
    result = optimize_sequence("MVHLTPEEK", organism="human",
                               objective=cai_gc_balanced_objective)

    # Or by name:
    result = optimize_sequence("MVHLTPEEK", organism="human",
                               objective="min_max_gc")

    # Or a custom callable:
    def my_obj(dna, protein, organism):
        gc = (dna.count("G") + dna.count("C")) / len(dna)
        return -abs(gc - 0.55)

    result = optimize_sequence("MVHLTPEEK", organism="human",
                               objective=my_obj)
"""

from __future__ import annotations

import math
from typing import Protocol, Callable, Union, runtime_checkable

__all__ = [
    "ObjectiveFunction",
    "cai_objective",
    "cai_gc_balanced_objective",
    "codon_pair_objective",
    "min_max_gc_objective",
    "resolve_objective",
    "OBJECTIVE_REGISTRY",
]


# ────────────────────────────────────────────────────────────
# Protocol
# ────────────────────────────────────────────────────────────

@runtime_checkable
class ObjectiveFunction(Protocol):
    """Protocol for custom optimization objectives.

    Any callable that accepts ``(dna, protein, organism)`` and returns
    a ``float`` satisfies this protocol.  Higher values are better.

    The ``dna`` argument is the full coding sequence (uppercase, length
    is a multiple of 3).  ``protein`` is the amino acid sequence
    (1-letter codes, no stop codon).  ``organism`` is the resolved
    canonical organism name (e.g. ``"Escherichia_coli"``).
    """

    def __call__(self, dna: str, protein: str, organism: str) -> float:
        """Evaluate the objective. Higher is better."""
        ...


# ────────────────────────────────────────────────────────────
# Built-in objectives
# ────────────────────────────────────────────────────────────

def cai_objective(dna: str, protein: str, organism: str) -> float:
    """Maximize Codon Adaptation Index (default objective).

    This is the standard CAI maximization that the optimizer has always
    used.  Including it as an explicit objective makes it easy to swap
    in alternative objectives or compose it with others.

    Args:
        dna: Full coding DNA sequence (uppercase, length % 3 == 0).
        protein: Amino acid sequence (1-letter codes, no stop).
        organism: Resolved canonical organism name.

    Returns:
        CAI value in [0.0, 1.0]. Returns 0.0 for empty/invalid sequences.
    """
    if not dna or len(dna) < 3:
        return 0.0
    from .translation import compute_cai
    try:
        return compute_cai(dna, organism=organism)
    except Exception:
        return 0.0


def cai_gc_balanced_objective(
    dna: str,
    protein: str,
    organism: str,
    gc_weight: float = 0.3,
) -> float:
    """Balance CAI with GC content near the organism target.

    The combined score is::

        (1 - gc_weight) * CAI + gc_weight * gc_score

    where ``gc_score`` is 1.0 when GC equals the organism's target GC,
    falling off linearly toward 0.0 as GC deviates from target.

    Args:
        dna: Full coding DNA sequence.
        protein: Amino acid sequence.
        organism: Resolved canonical organism name.
        gc_weight: Weight for the GC term (0–1). Default 0.3 means
            70% CAI / 30% GC balance.

    Returns:
        Combined score in [0.0, 1.0].
    """
    if not dna:
        return 0.0

    from .translation import compute_cai
    from .organisms import ORGANISM_GC_TARGETS

    try:
        cai_val = compute_cai(dna, organism=organism)
    except Exception:
        cai_val = 0.0

    gc = (dna.count("G") + dna.count("C")) / len(dna)
    gc_lo, gc_hi = ORGANISM_GC_TARGETS.get(organism, (0.40, 0.60))
    target_gc = (gc_lo + gc_hi) / 2.0
    max_deviation = max(abs(1.0 - target_gc), abs(target_gc - 0.0))
    deviation = abs(gc - target_gc)
    gc_score = max(0.0, 1.0 - deviation / max(max_deviation, 1e-9))

    return (1.0 - gc_weight) * cai_val + gc_weight * gc_score


def codon_pair_objective(dna: str, protein: str, organism: str) -> float:
    """Maximize codon pair score (avoid rare pairs).

    Returns the mean codon pair bias (CPB) score for the sequence.
    Positive CPB indicates over-represented (favoured) codon pairs;
    negative CPB indicates under-represented (disfavoured) pairs.

    This objective normalises the raw CPB to approximately [0.0, 1.0]
    so it can be compared with other objectives.  The mapping is::

        normalised = 0.5 + cpb / (2.0 * max_abs)

    where ``max_abs`` is the maximum absolute CPB in the organism's
    data table (clamped to 0.3 to avoid division by tiny numbers).

    Args:
        dna: Full coding DNA sequence.
        protein: Amino acid sequence.
        organism: Resolved canonical organism name.

    Returns:
        Normalised codon pair score, approximately in [0.0, 1.0].
    """
    if not dna or len(dna) < 6:
        return 0.5

    from .codon_pair_scoring import compute_cpb, get_codon_pair_data

    try:
        cpb = compute_cpb(dna, organism)
    except Exception:
        return 0.5

    # Normalise: published CPB typically ranges ~[-0.3, +0.3]
    # Use the organism data to determine the actual range
    data = get_codon_pair_data(organism)
    if data:
        max_abs = max(abs(v) for v in data.values())
        max_abs = max(max_abs, 0.3)  # floor to avoid amplifying noise
    else:
        max_abs = 0.3

    return max(0.0, min(1.0, 0.5 + cpb / (2.0 * max_abs)))


def min_max_gc_objective(
    dna: str,
    protein: str,
    organism: str,
    target_gc: float = 0.55,
) -> float:
    """Minimize deviation from target GC.

    Returns 1.0 when GC equals the target, falling off quadratically
    as GC deviates.  The quadratic penalty makes the optimizer strongly
    prefer sequences close to the target GC.

    Args:
        dna: Full coding DNA sequence.
        protein: Amino acid sequence.
        organism: Resolved canonical organism name.
        target_gc: Desired GC fraction. Default 0.55.

    Returns:
        Score in [0.0, 1.0].  1.0 = exactly on target.
    """
    if not dna:
        return 0.0

    gc = (dna.count("G") + dna.count("C")) / len(dna)
    deviation = abs(gc - target_gc)
    # Quadratic penalty: small deviations are tolerated, large ones are not
    return max(0.0, 1.0 - deviation * deviation * 4.0)


# ────────────────────────────────────────────────────────────
# Registry & resolution
# ────────────────────────────────────────────────────────────

OBJECTIVE_REGISTRY: dict[str, Callable[..., float]] = {
    "cai": cai_objective,
    "cai_gc_balanced": cai_gc_balanced_objective,
    "codon_pair": codon_pair_objective,
    "min_max_gc": min_max_gc_objective,
}
"""Built-in objective name → callable mapping.

Keys:
    - ``"cai"``: :func:`cai_objective` (default)
    - ``"cai_gc_balanced"``: :func:`cai_gc_balanced_objective`
    - ``"codon_pair"``: :func:`codon_pair_objective`
    - ``"min_max_gc"``: :func:`min_max_gc_objective`
"""


def resolve_objective(
    objective: Union[str, Callable[..., float], None] = None,
) -> Callable[..., float]:
    """Resolve an objective specification to a callable.

    Args:
        objective: One of:
            - ``None`` (default): returns :func:`cai_objective`.
            - A string name: looked up in :data:`OBJECTIVE_REGISTRY`.
            - A callable: returned as-is.

    Returns:
        A callable ``(dna, protein, organism) -> float``.

    Raises:
        ValueError: If ``objective`` is a string not in the registry.
        TypeError: If ``objective`` is neither a string, callable, nor None.
    """
    if objective is None:
        return cai_objective

    if isinstance(objective, str):
        key = objective.lower().strip()
        if key not in OBJECTIVE_REGISTRY:
            valid = ", ".join(sorted(OBJECTIVE_REGISTRY.keys()))
            raise ValueError(
                f"Unknown objective name {objective!r}. "
                f"Valid names: {valid}"
            )
        return OBJECTIVE_REGISTRY[key]

    if callable(objective):
        return objective

    raise TypeError(
        f"objective must be a string name, callable, or None; "
        f"got {type(objective).__name__}"
    )
