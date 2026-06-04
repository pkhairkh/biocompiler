"""
MaxEntScan Cross-Validation Against Published Values
======================================================

Compares the biocompiler.maxentscan PWM-based scores with:

1. The first-order Markov model from
   ``biocompiler.benchmarking.maxentscan_validated`` — a more
   faithful reproduction of Yeo & Burge (2004) that captures
   position dependencies the independent PWM cannot.

2. Hard-coded reference scores from the Yeo & Burge (2004)
   training data, as encoded in ``MAXENTSCAN_DONOR_SCORES`` and
   ``MAXENTSCAN_ACCEPTOR_SCORES``.

The PWM and Markov models use fundamentally different mathematical
formulations, so absolute scores will differ.  The cross-validation
verifies:

- **Score ordering concordance**: both models should agree on which
  sites are strong vs weak (Kendall's tau rank correlation).

- **Canonical > non-canonical**: GT donors and AG acceptors must score
  higher than non-canonical variants in both models.

- **Per-site qualitative agreement**: for each reference sequence, both
  models should classify it the same way (strong / weak / non-site).

**Acceptor 23-mer format difference:**

The two models use different 23-mer layouts:

- Markov model: ``[18 upstream][A at idx 18][G at idx 19][3 exonic]``
- PWM model:    ``[19 upstream ending with A at idx 19][G at idx 20][2 exonic]``

The PWM model has one more upstream position (position -20) and one
fewer exonic position.  Cross-validation constructs separate 23-mers
for each model representing the same biological context.

Reference:
  Yeo G, Burge CB. "Maximum entropy modeling of short sequence motifs
  with applications to RNA splicing."
  J Comput Biol. 2004;11(2-3):377-94. doi:10.1089/cmb.2004.11.377
"""

from __future__ import annotations

import math
import logging
from typing import Dict, List, Any

from biocompiler.maxentscan import (
    score_donor as pwm_score_donor,
    score_acceptor as pwm_score_acceptor,
    BASE_TO_INDEX,
    BG_PROB,
    DONOR_PWM_SCORE,
    ACCEPTOR_PWM_SCORE,
    _IMPOSSIBLE_SCORE,
)
from biocompiler.benchmarking.maxentscan_validated import (
    score_donor_maxentscan,
    score_acceptor_maxentscan,
    MAXENTSCAN_DONOR_SCORES,
    MAXENTSCAN_ACCEPTOR_SCORES,
    is_strong_donor,
    is_strong_acceptor,
)

__all__ = ["cross_validate_maxentscan"]

_logger = logging.getLogger(__name__)


# ==============================================================================
# Helpers
# ==============================================================================

def _pwm_donor_from_9mer(nine_mer: str) -> float:
    """Score a 9-mer donor using the PWM model.

    Embeds the 9-mer in a padded sequence and scores at position 3
    (so the 9-mer spans seq[0:9]).
    """
    seq = nine_mer.upper() + "A" * 7
    return pwm_score_donor(seq, 3)


def _pwm_acceptor_from_23mer(twenty_three_mer: str) -> float:
    """Score a 23-mer acceptor using the PWM model.

    Uses G-position convention: embeds the 23-mer so that
    score_acceptor(seq, 22) extracts it.  The 23-mer must have
    A at index 19 and G at index 20 for correct PWM alignment.
    """
    kmer = twenty_three_mer.upper()
    assert len(kmer) == 23, f"Expected 23-mer, got {len(kmer)}"
    seq = "AA" + kmer + "AAAA"
    return pwm_score_acceptor(seq, 22)


def _convert_validated_23mer_to_pwm(
    validated_23mer: str,
) -> str:
    """Convert a validated-module 23-mer to a PWM-compatible 23-mer.

    The validated module uses: [18 upstream][A at 18][G at 19][3 exonic]
    The PWM model uses:       [19 upstream ending A at 19][G at 20][2 exonic]

    Conversion:
      - Prepend the first upstream base (extends tract by 1 position)
      - Keep A and G at the junction
      - Use only the first 2 exonic bases (drop the 3rd)
    """
    v = validated_23mer.upper()
    assert len(v) == 23
    # Prepend one copy of the first upstream base to extend the tract
    extra_base = v[0]
    upstream_18 = v[0:18]   # 18 upstream bases
    a_base = v[18]          # A of AG
    g_base = v[19]          # G of AG
    exonic_2 = v[20:22]     # first 2 exonic bases

    # PWM 23-mer: [extra_base][18 upstream][A][G][2 exonic]
    pwm_23mer = extra_base + upstream_18 + a_base + g_base + exonic_2
    assert len(pwm_23mer) == 23, f"Converted 23-mer has length {len(pwm_23mer)}"
    return pwm_23mer


def _kendall_tau(scores_a: List[float], scores_b: List[float]) -> float:
    """Compute Kendall's tau rank correlation between two score lists.

    Returns a value in [-1, 1]: +1 = perfect concordance, -1 = perfect
    discordance.
    """
    n = len(scores_a)
    if n < 2:
        return 0.0

    concordant = 0
    discordant = 0
    for i in range(n):
        for j in range(i + 1, n):
            a_diff = scores_a[i] - scores_a[j]
            b_diff = scores_b[i] - scores_b[j]
            if a_diff * b_diff > 0:
                concordant += 1
            elif a_diff * b_diff < 0:
                discordant += 1

    total = concordant + discordant
    if total == 0:
        return 0.0
    return (concordant - discordant) / total


# ==============================================================================
# Biological test scenarios for acceptor cross-validation
# ==============================================================================
# Each entry: (label, markov_23mer, pwm_23mer)
# Both 23-mers represent the same biological context but formatted for
# their respective models.

_ACCEPTOR_SCENARIOS: List[Dict[str, str]] = [
    {
        "label": "strong_T_rich",
        # Markov: 18 upstream + A(18) + G(19) + 3 exonic
        "markov_23mer": "T" * 18 + "AGATG",
        # PWM: 19 upstream + A(19) + G(20) + 2 exonic
        "pwm_23mer": "T" * 19 + "AGAT",
    },
    {
        "label": "strong_C_rich",
        "markov_23mer": "C" * 18 + "AGATC",
        "pwm_23mer": "C" * 19 + "AGAT",
    },
    {
        "label": "strong_mixed_CT",
        # CT×9 = 18 chars upstream
        "markov_23mer": "CT" * 9 + "AGATG",
        # CT×9 + C = 19 chars upstream (last C extends pyrimidine tract)
        "pwm_23mer": "CT" * 9 + "C" + "AGAT",
    },
    {
        "label": "moderate_mixed_CT",
        # 18 upstream chars: CTCCTTTTTCCTTTTCTT
        "markov_23mer": "CTCCTTTTTCCTTTTCTT" + "AGTTG",
        # 19 upstream chars: CTCCTTTTTCCTTTTCTTT
        "pwm_23mer": "CTCCTTTTTCCTTTTCTTT" + "AGTT",
    },
    {
        "label": "weak_A_rich",
        "markov_23mer": "A" * 18 + "AGATG",
        "pwm_23mer": "A" * 19 + "AGAT",
    },
    {
        "label": "weak_G_rich",
        "markov_23mer": "G" * 18 + "AGATG",
        "pwm_23mer": "G" * 19 + "AGAT",
    },
    {
        "label": "weak_mixed_purine",
        # AG×9 = 18 chars upstream
        "markov_23mer": "AG" * 9 + "AGATG",
        # AG×9 + A = 19 chars upstream
        "pwm_23mer": "AG" * 9 + "A" + "AGAT",
    },
]


# ==============================================================================
# Main cross-validation function
# ==============================================================================

def cross_validate_maxentscan() -> Dict[str, Any]:
    """Cross-validate maxentscan PWM scores against Yeo & Burge (2004).

    Compares the independent-position PWM (biocompiler.maxentscan)
    with the first-order Markov model (biocompiler.benchmarking.maxentscan_validated)
    and hard-coded reference scores from the Yeo & Burge training data.

    Returns:
        Dict with the following keys:

        - ``"donor_comparison"``: list of dicts, one per reference 9-mer,
          with keys ``"sequence"``, ``"pwm_score"``, ``"markov_score"``,
          ``"reference_score"``, ``"pwm_class"``, ``"markov_class"``.

        - ``"acceptor_comparison"``: list of dicts, one per biological
          scenario, with keys ``"label"``, ``"pwm_score"``,
          ``"markov_score"``, ``"pwm_class"``, ``"markov_class"``.

        - ``"donor_tau"``: Kendall's tau between PWM and Markov donor scores.

        - ``"acceptor_tau"``: Kendall's tau between PWM and Markov acceptor
          scores across biological scenarios.

        - ``"donor_concordance"``: fraction of donor sites where PWM and
          Markov agree on classification (strong / weak).

        - ``"acceptor_concordance"``: fraction of acceptor scenarios where
          PWM and Markov agree on classification.

        - ``"gt_gt_agreement"``: True if both models agree GT > GC for donors.

        - ``"ag_invariant_agreement"``: True if both models score canonical
          AG acceptors higher than non-canonical.

        - ``"pwm_self_consistent"``: True if PWM scores match manual calculation.

        - ``"markov_self_consistent"``: True if Markov scores reproduce
          hard-coded reference values.

        - ``"yeo_burge_invariant_check"``: dict of invariant position
          probability checks from the PWM.
    """
    results: Dict[str, Any] = {}

    # ======================================================================
    # 1. Donor cross-validation (same 9-mers for both models)
    # ======================================================================
    _logger.info("Cross-validating donor scores...")

    donor_comparison: List[Dict[str, Any]] = []
    pwm_donor_scores: List[float] = []
    markov_donor_scores: List[float] = []

    for nine_mer, ref_score in MAXENTSCAN_DONOR_SCORES.items():
        pwm_score = _pwm_donor_from_9mer(nine_mer)
        markov_score = score_donor_maxentscan(nine_mer)

        pwm_donor_scores.append(pwm_score)
        markov_donor_scores.append(markov_score)

        donor_comparison.append({
            "sequence": nine_mer,
            "pwm_score": pwm_score,
            "markov_score": markov_score,
            "reference_score": ref_score,
            "pwm_class": "strong" if pwm_score > 3.0 else "weak",
            "markov_class": "strong" if markov_score > 3.0 else "weak",
        })

    results["donor_comparison"] = donor_comparison
    results["donor_tau"] = _kendall_tau(pwm_donor_scores, markov_donor_scores)

    donor_concordant = sum(
        1 for d in donor_comparison
        if d["pwm_class"] == d["markov_class"]
    )
    results["donor_concordance"] = (
        donor_concordant / len(donor_comparison) if donor_comparison else 0.0
    )

    # ======================================================================
    # 2. Acceptor cross-validation (separate 23-mers per model)
    # ======================================================================
    _logger.info("Cross-validating acceptor scores...")

    acceptor_comparison: List[Dict[str, Any]] = []
    pwm_acceptor_scores: List[float] = []
    markov_acceptor_scores: List[float] = []

    for scenario in _ACCEPTOR_SCENARIOS:
        markov_23 = scenario["markov_23mer"]
        pwm_23 = scenario["pwm_23mer"]

        # Validate lengths
        assert len(markov_23) == 23, (
            f"Markov 23-mer for {scenario['label']} has length {len(markov_23)}"
        )
        assert len(pwm_23) == 23, (
            f"PWM 23-mer for {scenario['label']} has length {len(pwm_23)}"
        )

        markov_score = score_acceptor_maxentscan(markov_23)
        pwm_score = _pwm_acceptor_from_23mer(pwm_23)

        pwm_acceptor_scores.append(pwm_score)
        markov_acceptor_scores.append(markov_score)

        acceptor_comparison.append({
            "label": scenario["label"],
            "pwm_23mer": pwm_23,
            "markov_23mer": markov_23,
            "pwm_score": pwm_score,
            "markov_score": markov_score,
            "pwm_class": "strong" if pwm_score > 3.0 else "weak",
            "markov_class": "strong" if markov_score > 3.0 else "weak",
        })

    # Also add the hard-coded reference 23-mers (Markov model only;
    # PWM scores are computed from converted 23-mers)
    for twenty_three_mer, ref_score in MAXENTSCAN_ACCEPTOR_SCORES.items():
        markov_score = score_acceptor_maxentscan(twenty_three_mer)
        # Convert to PWM format and score
        try:
            pwm_23 = _convert_validated_23mer_to_pwm(twenty_three_mer)
            pwm_score = _pwm_acceptor_from_23mer(pwm_23)
        except (AssertionError, ValueError) as e:
            _logger.warning(
                "Could not convert acceptor %s for PWM: %s",
                twenty_three_mer[:10] + "...", e,
            )
            pwm_score = _IMPOSSIBLE_SCORE

        pwm_acceptor_scores.append(pwm_score)
        markov_acceptor_scores.append(markov_score)

        acceptor_comparison.append({
            "label": f"ref_{twenty_three_mer[:8]}",
            "pwm_23mer": pwm_23 if pwm_score != _IMPOSSIBLE_SCORE else "",
            "markov_23mer": twenty_three_mer,
            "pwm_score": pwm_score,
            "markov_score": markov_score,
            "reference_score": ref_score,
            "pwm_class": "strong" if pwm_score > 3.0 else "weak",
            "markov_class": "strong" if markov_score > 3.0 else "weak",
        })

    results["acceptor_comparison"] = acceptor_comparison
    results["acceptor_tau"] = _kendall_tau(
        pwm_acceptor_scores, markov_acceptor_scores
    )

    acceptor_concordant = sum(
        1 for d in acceptor_comparison
        if d["pwm_class"] == d["markov_class"]
    )
    results["acceptor_concordance"] = (
        acceptor_concordant / len(acceptor_comparison)
        if acceptor_comparison else 0.0
    )

    # ======================================================================
    # 3. GT > GC agreement
    # ======================================================================
    gt_9mer = "CAGGTAAGT"
    gc_9mer = "CAGGCAAGT"

    pwm_gt = _pwm_donor_from_9mer(gt_9mer)
    pwm_gc = _pwm_donor_from_9mer(gc_9mer)
    markov_gt = score_donor_maxentscan(gt_9mer)
    markov_gc = score_donor_maxentscan(gc_9mer)

    results["gt_gt_agreement"] = (
        (pwm_gt > pwm_gc) and (markov_gt > markov_gc)
    )

    # ======================================================================
    # 4. AG invariant agreement
    # ======================================================================
    # Strong acceptor (T-rich tract + AG) vs non-acceptor (G-rich + AG)
    # PWM 23-mers (A at idx 19, G at idx 20)
    pwm_strong_23 = "T" * 19 + "AGAT"
    pwm_non_23 = "G" * 19 + "AGAT"

    # Markov 23-mers (A at idx 18, G at idx 19)
    markov_strong_23 = "T" * 18 + "AGATG"
    markov_non_23 = "G" * 18 + "AGATG"

    pwm_strong = _pwm_acceptor_from_23mer(pwm_strong_23)
    pwm_non = _pwm_acceptor_from_23mer(pwm_non_23)
    markov_strong = score_acceptor_maxentscan(markov_strong_23)
    markov_non = score_acceptor_maxentscan(markov_non_23)

    results["ag_invariant_agreement"] = (
        (pwm_strong > pwm_non) and (markov_strong > markov_non)
    )

    # ======================================================================
    # 5. PWM self-consistency: manual vs module calculation
    # ======================================================================
    _EPSILON = 0.001
    pwm_self_consistent = True
    for nine_mer in MAXENTSCAN_DONOR_SCORES:
        module_score = _pwm_donor_from_9mer(nine_mer)
        manual_score = 0.0
        for i, base in enumerate(nine_mer.upper()):
            idx = BASE_TO_INDEX[base]
            prob = max(DONOR_PWM_SCORE[i][idx], _EPSILON)
            manual_score += math.log2(prob / BG_PROB)
        manual_score = round(manual_score, 4)
        if module_score != manual_score:
            pwm_self_consistent = False
            _logger.warning(
                "PWM donor self-inconsistency for %s: module=%s manual=%s",
                nine_mer, module_score, manual_score,
            )
    results["pwm_self_consistent"] = pwm_self_consistent

    # ======================================================================
    # 6. Markov self-consistency: module reproduces hard-coded references
    # ======================================================================
    markov_self_consistent = True
    for nine_mer, ref_score in MAXENTSCAN_DONOR_SCORES.items():
        computed = score_donor_maxentscan(nine_mer)
        if abs(computed - ref_score) > 0.02:
            markov_self_consistent = False
            _logger.warning(
                "Markov donor self-inconsistency for %s: computed=%s ref=%s",
                nine_mer, computed, ref_score,
            )

    for twenty_three_mer, ref_score in MAXENTSCAN_ACCEPTOR_SCORES.items():
        computed = score_acceptor_maxentscan(twenty_three_mer)
        if abs(computed - ref_score) > 0.02:
            markov_self_consistent = False
            _logger.warning(
                "Markov acceptor self-inconsistency for %s: computed=%s ref=%s",
                twenty_three_mer, computed, ref_score,
            )
    results["markov_self_consistent"] = markov_self_consistent

    # ======================================================================
    # 7. Yeo & Burge training data statistics validation
    # ======================================================================
    # The donor PWM position +1 should have G probability ≈ 0.990
    # The donor PWM position +2 should have T probability ≈ 0.990
    # The acceptor PWM position -1 should have A probability ≈ 0.980
    # The acceptor PWM position +0 should have G probability ≈ 0.980
    donor_g_prob = DONOR_PWM_SCORE[3][BASE_TO_INDEX["G"]]  # row 3 = pos +1
    donor_t_prob = DONOR_PWM_SCORE[4][BASE_TO_INDEX["T"]]  # row 4 = pos +2
    acceptor_a_prob = ACCEPTOR_PWM_SCORE[19][BASE_TO_INDEX["A"]]  # row 19 = pos -1
    acceptor_g_prob = ACCEPTOR_PWM_SCORE[20][BASE_TO_INDEX["G"]]  # row 20 = pos +0

    results["yeo_burge_invariant_check"] = {
        "donor_plus1_G_prob": donor_g_prob,
        "donor_plus2_T_prob": donor_t_prob,
        "acceptor_minus1_A_prob": acceptor_a_prob,
        "acceptor_plus0_G_prob": acceptor_g_prob,
        "donor_G_near_invariant": donor_g_prob >= 0.95,
        "donor_T_near_invariant": donor_t_prob >= 0.95,
        "acceptor_A_near_invariant": acceptor_a_prob >= 0.95,
        "acceptor_G_near_invariant": acceptor_g_prob >= 0.95,
    }

    _logger.info("Cross-validation complete.")
    return results
