"""
BioCompiler Five-Valued Logic — Refinement Mapping and Soundness Bridge

This module bridges the gap between the runtime's 5-valued logic
(PASS / FAIL / UNCERTAIN / LIKELY_PASS / LIKELY_FAIL) and the Lean4
formal proofs which use 3-valued logic (PASS / FAIL / UNCERTAIN).

Refinement Mapping:
  The 5-valued verdicts are refined into 3-valued ones as follows:

    5-valued         →  3-valued       Rationale
    ───────────────────────────────────────────────────────────────────
    PASS             →  PASS           Definite pass (formally verified)
    LIKELY_PASS      →  UNCERTAIN      High confidence but not formally
                                       verified; conservatively mapped
    UNCERTAIN        →  UNCERTAIN      Insufficient evidence for pass/fail
    LIKELY_FAIL      →  FAIL           High confidence of failure;
                                       mapped to FAIL for soundness
    FAIL             →  FAIL           Definite failure (formally verified)

  Soundness guarantee: For any 5-valued verdict v,
    if to_three_valued(v) = PASS in the 3-valued model,
    then v must be PASS in the 5-valued model.
    This means the 3-valued PASS verdict is preserved exactly.

  The LIKELY_FAIL → FAIL mapping is conservative: we treat "probably
  failing" as "definitely failing" for soundness. This ensures that
  the Lean4 soundness proof (which only knows about PASS/FAIL/UNCERTAIN)
  remains valid when applied to 5-valued verdicts.

  The LIKELY_PASS → UNCERTAIN mapping is also conservative: we treat
  "probably passing" as "uncertain" rather than "passing". This
  ensures that the 3-valued model never claims PASS unless it is
  definitely PASS.

Confidence Score Ordering:
  PASS       → 1.0   (formally verified)
  LIKELY_PASS → 0.75  (high confidence, e.g., tool available + good result)
  UNCERTAIN  → 0.5   (insufficient evidence)
  LIKELY_FAIL → 0.25  (high confidence of failure)
  FAIL       → 0.0   (formally verified failure)

Lattice Combination:
  The combine_verdicts function implements a meet (greatest lower bound)
  operation on the 5-valued lattice:
    PASS > LIKELY_PASS > UNCERTAIN > LIKELY_FAIL > FAIL

  This matches the Kleene-style conjunction where AND takes the minimum.
"""

from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .types import Verdict


class FiveValuedVerdict(str, Enum):
    """Five-valued logic for type-check verdicts.

    Ordering: PASS > LIKELY_PASS > UNCERTAIN > LIKELY_FAIL > FAIL

    This extends the 3-valued logic (PASS/FAIL/UNCERTAIN) used in the
    Lean4 formal proofs with two intermediate verdicts:
    - LIKELY_PASS: High confidence of passing but not formally verified
    - LIKELY_FAIL: High confidence of failing but not formally verified
    """
    PASS = "PASS"
    LIKELY_PASS = "LIKELY_PASS"
    UNCERTAIN = "UNCERTAIN"
    LIKELY_FAIL = "LIKELY_FAIL"
    FAIL = "FAIL"


# Internal ordering: PASS > LIKELY_PASS > UNCERTAIN > LIKELY_FAIL > FAIL
_FIVE_VALUED_ORDER: dict[FiveValuedVerdict, int] = {
    FiveValuedVerdict.PASS: 4,
    FiveValuedVerdict.LIKELY_PASS: 3,
    FiveValuedVerdict.UNCERTAIN: 2,
    FiveValuedVerdict.LIKELY_FAIL: 1,
    FiveValuedVerdict.FAIL: 0,
}


def to_three_valued(v: FiveValuedVerdict) -> "Verdict":
    """Refine a 5-valued verdict to the 3-valued Lean4 model.

    Mapping:
      PASS        → PASS       (exact)
      LIKELY_PASS → UNCERTAIN  (conservative: not formally verified)
      UNCERTAIN   → UNCERTAIN  (exact)
      LIKELY_FAIL → FAIL       (conservative: treat as failure)
      FAIL        → FAIL       (exact)

    This refinement is sound: if to_three_valued(v) = PASS, then v
    must be PASS. The 3-valued PASS is preserved exactly.

    Args:
        v: A FiveValuedVerdict value.

    Returns:
        The corresponding Verdict in the 3-valued Lean4 model.
    """
    from .types import Verdict

    _REFINEMENT_MAP: dict[FiveValuedVerdict, Verdict] = {
        FiveValuedVerdict.PASS: Verdict.PASS,
        FiveValuedVerdict.LIKELY_PASS: Verdict.UNCERTAIN,
        FiveValuedVerdict.UNCERTAIN: Verdict.UNCERTAIN,
        FiveValuedVerdict.LIKELY_FAIL: Verdict.FAIL,
        FiveValuedVerdict.FAIL: Verdict.FAIL,
    }
    return _REFINEMENT_MAP[v]


def confidence_score(v: FiveValuedVerdict) -> float:
    """Compute a confidence score for a 5-valued verdict.

    Scores:
      PASS        → 1.0   (formally verified)
      LIKELY_PASS → 0.75  (high confidence)
      UNCERTAIN   → 0.5   (insufficient evidence)
      LIKELY_FAIL → 0.25  (high confidence of failure)
      FAIL        → 0.0   (formally verified failure)

    Args:
        v: A FiveValuedVerdict value.

    Returns:
        Float confidence score in [0.0, 1.0].
    """
    _CONFIDENCE_MAP: dict[FiveValuedVerdict, float] = {
        FiveValuedVerdict.PASS: 1.0,
        FiveValuedVerdict.LIKELY_PASS: 0.75,
        FiveValuedVerdict.UNCERTAIN: 0.5,
        FiveValuedVerdict.LIKELY_FAIL: 0.25,
        FiveValuedVerdict.FAIL: 0.0,
    }
    return _CONFIDENCE_MAP[v]


def combine_verdicts(verdicts: list[FiveValuedVerdict]) -> FiveValuedVerdict:
    """Combine multiple 5-valued verdicts using lattice meet (greatest lower bound).

    The combination follows the Kleene-style conjunction: AND takes the
    minimum in the ordering PASS > LIKELY_PASS > UNCERTAIN > LIKELY_FAIL > FAIL.
    This is the weakest-link principle: the overall verdict is determined
    by the worst individual verdict.

    An empty list returns UNCERTAIN (no evidence either way).

    Args:
        verdicts: List of FiveValuedVerdict values to combine.

    Returns:
        The combined FiveValuedVerdict.
    """
    if not verdicts:
        return FiveValuedVerdict.UNCERTAIN

    result = verdicts[0]
    for v in verdicts[1:]:
        if _FIVE_VALUED_ORDER[v] < _FIVE_VALUED_ORDER[result]:
            result = v
    return result


def verify_five_valued_soundness(
    verdict: FiveValuedVerdict,
    actual_condition: bool,
) -> bool:
    """Check whether a 5-valued verdict is sound with respect to the actual condition.

    A verdict is sound if it does not claim PASS when the condition is
    false, and does not claim FAIL when the condition is true.

    Soundness conditions:
      - PASS requires actual_condition = True
      - LIKELY_PASS requires actual_condition = True
      - UNCERTAIN is always sound (makes no definite claim)
      - LIKELY_FAIL requires actual_condition = False
      - FAIL requires actual_condition = False

    Note: UNCERTAIN is always sound because it makes no definite claim.
    This is consistent with the 3-valued model where UNCERTAIN means
    "insufficient evidence to claim PASS or FAIL".

    Args:
        verdict: The FiveValuedVerdict to check.
        actual_condition: The ground truth boolean condition.

    Returns:
        True if the verdict is sound, False otherwise.
    """
    if verdict == FiveValuedVerdict.PASS:
        return actual_condition is True
    elif verdict == FiveValuedVerdict.LIKELY_PASS:
        return actual_condition is True
    elif verdict == FiveValuedVerdict.UNCERTAIN:
        return True  # UNCERTAIN makes no definite claim
    elif verdict == FiveValuedVerdict.LIKELY_FAIL:
        return actual_condition is False
    elif verdict == FiveValuedVerdict.FAIL:
        return actual_condition is False
    return False  # Unreachable for valid enum values


def five_valued_and(a: FiveValuedVerdict, b: FiveValuedVerdict) -> FiveValuedVerdict:
    """Conjunction in five-valued logic (Kleene-style). AND takes the minimum."""
    if _FIVE_VALUED_ORDER[a] <= _FIVE_VALUED_ORDER[b]:
        return a
    return b


def five_valued_or(a: FiveValuedVerdict, b: FiveValuedVerdict) -> FiveValuedVerdict:
    """Disjunction in five-valued logic (Kleene-style). OR takes the maximum."""
    if _FIVE_VALUED_ORDER[a] >= _FIVE_VALUED_ORDER[b]:
        return a
    return b


def refinement_is_sound(verdict: FiveValuedVerdict, actual_condition: bool) -> bool:
    """Check that the refinement to 3-valued logic preserves soundness.

    For any 5-valued verdict v and actual condition c:
      - If to_three_valued(v) = PASS, then v = PASS and c = True
      - If to_three_valued(v) = FAIL, then c = False (conservative for LIKELY_FAIL)
      - If to_three_valued(v) = UNCERTAIN, soundness is trivially preserved

    This function verifies that the refinement does not introduce unsoundness.
    """
    three_v = to_three_valued(verdict)
    from .types import Verdict

    if three_v == Verdict.PASS:
        # PASS in 3-valued model requires actual_condition = True
        return actual_condition is True
    elif three_v == Verdict.FAIL:
        # FAIL in 3-valued model is conservative: LIKELY_FAIL → FAIL
        # Sound if actual_condition is False, or if verdict is LIKELY_FAIL
        # and we accept conservative failure
        if verdict == FiveValuedVerdict.LIKELY_FAIL:
            # LIKELY_FAIL mapped to FAIL: soundness requires actual_condition = False
            return actual_condition is False
        return actual_condition is False
    else:
        # UNCERTAIN: always sound (makes no definite claim)
        return True


__all__ = [
    "FiveValuedVerdict",
    "to_three_valued",
    "confidence_score",
    "combine_verdicts",
    "verify_five_valued_soundness",
    "five_valued_and",
    "five_valued_or",
    "refinement_is_sound",
]
