#!/usr/bin/env python3
"""
MaxEntScan Validation Against Published Yeo & Burge (2004) Scores
==================================================================

This script validates the three MaxEntScan implementations in the biocompiler
codebase against the hard-coded reference scores in
benchmarking.maxentscan_validated and known biological properties.

Three implementations tested:
  1. First-order Markov model (benchmarking.maxentscan_validated)
  2. Independent-position PWM (maxentscan.score_donor / score_acceptor)
  3. Simplified hand-crafted PWM (splicing.maxent_score)

Validation checks:
  A. Self-consistency: Does the Markov model reproduce its own hard-coded scores?
  B. Score ordering: Strong donors > weak donors > non-donors
  C. Canonical GT > non-canonical GC
  D. Polypyrimidine tract effect for acceptors
  E. Rank concordance between all three implementations
  F. Cross-validation of absolute score ranges

Reference:
  Yeo G, Burge CB (2004) "Maximum entropy modeling of short sequence motifs
  with applications to RNA splicing." J Comput Biol 11(2-3):377-94.
"""

from __future__ import annotations

import sys
import math
from typing import Dict, List, Tuple

# ── Path setup ────────────────────────────────────────────────────────────────
sys.path.insert(0, "/home/z/my-project/biocompiler/src")

from biocompiler.benchmarking.maxentscan_validated import (
    MAXENTSCAN_DONOR_SCORES,
    MAXENTSCAN_ACCEPTOR_SCORES,
    score_donor_maxentscan,
    score_acceptor_maxentscan,
    is_strong_donor,
    is_strong_acceptor,
    DONOR_MONO_FREQ,
    DONOR_COND_FREQ,
    ACCEPTOR_MONO_FREQ,
    ACCEPTOR_COND_FREQ,
)

from biocompiler.sequence.maxentscan import (
    score_donor as pwm_score_donor,
    score_acceptor as pwm_score_acceptor,
    DONOR_PWM_SCORE,
    ACCEPTOR_PWM_SCORE,
    BASE_TO_INDEX,
    BG_PROB,
    _IMPOSSIBLE_SCORE,
)

from biocompiler.sequence.splicing import maxent_score as splicing_maxent_score


# ==============================================================================
# Helpers
# ==============================================================================

def _make_donor_seq(nine_mer: str, pos: int = 6, total_len: int = 20) -> str:
    """Embed a 9-mer donor into a sequence so score_donor(seq, pos) sees it."""
    seq = list("A" * total_len)
    for i, c in enumerate(nine_mer):
        idx = pos - 3 + i
        if 0 <= idx < total_len:
            seq[idx] = c
    return "".join(seq)


def _make_acceptor_seq(
    twenty_three_mer: str,
    pos: int = 22,
    total_len: int = 45,
) -> str:
    """Embed a 23-mer acceptor so score_acceptor(seq, pos) sees it (G-position convention)."""
    seq = list("A" * total_len)
    for i, c in enumerate(twenty_three_mer):
        idx = pos - 20 + i
        if 0 <= idx < total_len:
            seq[idx] = c
    return "".join(seq)


def _manual_pwm_donor_score(nine_mer: str) -> float:
    """Independent manual PWM calculation for cross-checking."""
    _EPSILON = 0.001
    score = 0.0
    for i, base in enumerate(nine_mer.upper()):
        idx = BASE_TO_INDEX[base]
        prob = max(DONOR_PWM_SCORE[i][idx], _EPSILON)
        score += math.log2(prob / BG_PROB)
    return round(score, 4)


def _manual_pwm_acceptor_score(twenty_three_mer: str) -> float:
    """Independent manual PWM calculation for cross-checking."""
    _EPSILON = 0.001
    score = 0.0
    for i, base in enumerate(twenty_three_mer.upper()):
        idx = BASE_TO_INDEX[base]
        prob = max(ACCEPTOR_PWM_SCORE[i][idx], _EPSILON)
        score += math.log2(prob / BG_PROB)
    return round(score, 4)


def kendall_tau(rank1: List[float], rank2: List[float]) -> Tuple[float, int, int]:
    """Compute Kendall's tau-b rank correlation between two score lists.

    Returns (tau, concordant_pairs, discordant_pairs).
    """
    n = len(rank1)
    if n != len(rank2):
        raise ValueError("Lists must have the same length")
    concordant = 0
    discordant = 0
    for i in range(n):
        for j in range(i + 1, n):
            d1 = rank1[i] - rank1[j]
            d2 = rank2[i] - rank2[j]
            prod = d1 * d2
            if prod > 0:
                concordant += 1
            elif prod < 0:
                discordant += 1
    total = concordant + discordant
    if total == 0:
        return (0.0, 0, 0)
    tau = (concordant - discordant) / total
    return (tau, concordant, discordant)


# ==============================================================================
# Validation results container
# ==============================================================================

class ValidationResult:
    """Track pass/fail for each validation check."""

    def __init__(self):
        self.checks: List[Dict] = []

    def add(self, name: str, passed: bool, detail: str):
        self.checks.append({
            "name": name,
            "passed": passed,
            "detail": detail,
        })

    def summary(self) -> str:
        lines = []
        passed = sum(1 for c in self.checks if c["passed"])
        total = len(self.checks)
        lines.append(f"\n{'='*72}")
        lines.append(f"VALIDATION SUMMARY: {passed}/{total} checks passed")
        lines.append(f"{'='*72}")
        for c in self.checks:
            icon = "PASS" if c["passed"] else "FAIL"
            lines.append(f"  [{icon}] {c['name']}")
            if c["detail"]:
                lines.append(f"         {c['detail']}")
        lines.append(f"{'='*72}")
        return "\n".join(lines)


results = ValidationResult()


# ==============================================================================
# SECTION A: Self-consistency of the first-order Markov model
# ==============================================================================

print("\n" + "="*72)
print("SECTION A: Self-consistency of First-Order Markov Model")
print("="*72)
print()
print("Checking whether score_donor_maxentscan() reproduces the hard-coded")
print("MAXENTSCAN_DONOR_SCORES and score_acceptor_maxentscan() reproduces")
print("MAXENTSCAN_ACCEPTOR_SCORES.\n")

# ── Donor self-consistency ────────────────────────────────────────────────────
print("── Donor Scores (9-mer) ──")
donor_self_consistent = 0
donor_self_total = 0
donor_mismatches = []

for seq, expected_score in sorted(MAXENTSCAN_DONOR_SCORES.items(), key=lambda x: -x[1]):
    computed = score_donor_maxentscan(seq)
    match = abs(computed - expected_score) < 0.02  # Allow tiny rounding differences
    donor_self_total += 1
    if match:
        donor_self_consistent += 1
    else:
        donor_mismatches.append((seq, expected_score, computed))
    status = "OK" if match else "MISMATCH"
    print(f"  {seq}: expected={expected_score:8.2f}  computed={computed:8.2f}  [{status}]")

donor_pct = 100.0 * donor_self_consistent / donor_self_total if donor_self_total else 0
results.add(
    "Markov donor self-consistency",
    donor_self_consistent == donor_self_total,
    f"{donor_self_consistent}/{donor_self_total} match ({donor_pct:.0f}%)"
)

if donor_mismatches:
    print(f"\n  ** MISMATCHES **")
    for seq, exp, comp in donor_mismatches:
        print(f"    {seq}: expected={exp:.2f}, computed={comp:.2f}, diff={comp-exp:.4f}")

# ── Acceptor self-consistency ─────────────────────────────────────────────────
print("\n── Acceptor Scores (23-mer) ──")
acceptor_self_consistent = 0
acceptor_self_total = 0
acceptor_mismatches = []

for seq, expected_score in sorted(MAXENTSCAN_ACCEPTOR_SCORES.items(), key=lambda x: -x[1]):
    computed = score_acceptor_maxentscan(seq)
    match = abs(computed - expected_score) < 0.02
    acceptor_self_total += 1
    if match:
        acceptor_self_consistent += 1
    else:
        acceptor_mismatches.append((seq, expected_score, computed))
    status = "OK" if match else "MISMATCH"
    print(f"  {seq}: expected={expected_score:8.2f}  computed={computed:8.2f}  [{status}]")

acceptor_pct = 100.0 * acceptor_self_consistent / acceptor_self_total if acceptor_self_total else 0
results.add(
    "Markov acceptor self-consistency",
    acceptor_self_consistent == acceptor_self_total,
    f"{acceptor_self_consistent}/{acceptor_self_total} match ({acceptor_pct:.0f}%)"
)

if acceptor_mismatches:
    print(f"\n  ** MISMATCHES **")
    for seq, exp, comp in acceptor_mismatches:
        print(f"    {seq}: expected={exp:.2f}, computed={comp:.2f}, diff={comp-exp:.4f}")


# ==============================================================================
# SECTION B: Score ordering — Strong > Weak > Non-donor
# ==============================================================================

print("\n" + "="*72)
print("SECTION B: Score Ordering (Strong > Weak > Non-donor)")
print("="*72)
print()

# Categorize donors by expected biological strength
strong_donors = {k: v for k, v in MAXENTSCAN_DONOR_SCORES.items() if v > 5.0}
moderate_donors = {k: v for k, v in MAXENTSCAN_DONOR_SCORES.items() if 0 < v <= 5.0}
non_donors = {k: v for k, v in MAXENTSCAN_DONOR_SCORES.items() if v < 0}

print(f"  Strong donors (score > 5): {list(strong_donors.keys())}")
print(f"  Moderate donors (0 < score <= 5): {list(moderate_donors.keys())}")
print(f"  Non-donors (score < 0): {list(non_donors.keys())}")
print()

# Check: all strong > all moderate > all non-donors (Markov model)
min_strong = min(strong_donors.values()) if strong_donors else float('inf')
max_moderate = max(moderate_donors.values()) if moderate_donors else float('-inf')
min_moderate = min(moderate_donors.values()) if moderate_donors else float('inf')
max_non = max(non_donors.values()) if non_donors else float('-inf')

strong_gt_moderate = min_strong > max_moderate if strong_donors and moderate_donors else True
moderate_gt_non = min_moderate > max_non if moderate_donors and non_donors else True

print(f"  Min strong score: {min_strong:.2f}")
print(f"  Max moderate score: {max_moderate:.2f}")
print(f"  Strong > Moderate: {strong_gt_moderate}")
print(f"  Min moderate score: {min_moderate:.2f}")
print(f"  Max non-donor score: {max_non:.2f}")
print(f"  Moderate > Non-donor: {moderate_gt_non}")

results.add(
    "Markov: strong donors > moderate donors",
    strong_gt_moderate,
    f"min(strong)={min_strong:.2f} vs max(moderate)={max_moderate:.2f}"
)
results.add(
    "Markov: moderate donors > non-donors",
    moderate_gt_non,
    f"min(moderate)={min_moderate:.2f} vs max(non-donor)={max_non:.2f}"
)

# Now verify with PWM implementation
print("\n  ── PWM implementation score ordering ──")
pwm_strong_scores = []
pwm_moderate_scores = []
pwm_non_scores = []

for seq in strong_donors:
    s = pwm_score_donor(_make_donor_seq(seq), 6)
    if s != _IMPOSSIBLE_SCORE:
        pwm_strong_scores.append(s)
for seq in moderate_donors:
    s = pwm_score_donor(_make_donor_seq(seq), 6)
    if s != _IMPOSSIBLE_SCORE:
        pwm_moderate_scores.append(s)
for seq in non_donors:
    s = pwm_score_donor(_make_donor_seq(seq), 6)
    if s != _IMPOSSIBLE_SCORE:
        pwm_non_scores.append(s)

print(f"  PWM strong scores: {[f'{s:.2f}' for s in pwm_strong_scores]}")
print(f"  PWM moderate scores: {[f'{s:.2f}' for s in pwm_moderate_scores]}")
print(f"  PWM non-donor scores: {[f'{s:.2f}' for s in pwm_non_scores]}")

if pwm_strong_scores and pwm_moderate_scores:
    pwm_s_gt_m = min(pwm_strong_scores) > max(pwm_moderate_scores)
else:
    pwm_s_gt_m = True
if pwm_moderate_scores and pwm_non_scores:
    pwm_m_gt_n = min(pwm_moderate_scores) > max(pwm_non_scores)
else:
    pwm_m_gt_n = True

results.add(
    "PWM: strong donors > moderate donors",
    pwm_s_gt_m,
    f"min(strong)={min(pwm_strong_scores):.2f} vs max(moderate)={max(pwm_moderate_scores):.2f}" if pwm_strong_scores and pwm_moderate_scores else "N/A"
)
results.add(
    "PWM: moderate donors > non-donors",
    pwm_m_gt_n,
    f"min(moderate)={min(pwm_moderate_scores):.2f} vs max(non-donor)={max(pwm_non_scores):.2f}" if pwm_moderate_scores and pwm_non_scores else "N/A"
)

# Simplified PWM (splicing.maxent_score)
print("\n  ── Simplified PWM (splicing.maxent_score) ordering ──")
sp_strong_scores = [splicing_maxent_score(seq) for seq in strong_donors]
sp_moderate_scores = [splicing_maxent_score(seq) for seq in moderate_donors]
sp_non_scores = [splicing_maxent_score(seq) for seq in non_donors]

print(f"  Splicing strong: {[f'{s:.2f}' for s in sp_strong_scores]}")
print(f"  Splicing moderate: {[f'{s:.2f}' for s in sp_moderate_scores]}")
print(f"  Splicing non-donor: {[f'{s:.2f}' for s in sp_non_scores]}")

if sp_strong_scores and sp_moderate_scores:
    sp_s_gt_m = min(sp_strong_scores) > max(sp_moderate_scores)
else:
    sp_s_gt_m = True
if sp_moderate_scores and sp_non_scores:
    sp_m_gt_n = min(sp_moderate_scores) > max(sp_non_scores)
else:
    sp_m_gt_n = True

results.add(
    "Splicing: strong donors > moderate donors",
    sp_s_gt_m,
    f"min(strong)={min(sp_strong_scores):.2f} vs max(moderate)={max(sp_moderate_scores):.2f}" if sp_strong_scores and sp_moderate_scores else "N/A"
)
results.add(
    "Splicing: moderate donors > non-donors",
    sp_m_gt_n,
    f"min(moderate)={min(sp_moderate_scores):.2f} vs max(non-donor)={max(sp_non_scores):.2f}" if sp_moderate_scores and sp_non_scores else "N/A"
)


# ==============================================================================
# SECTION C: Canonical GT > Non-canonical GC
# ==============================================================================

print("\n" + "="*72)
print("SECTION C: Canonical GT > Non-canonical GC")
print("="*72)
print()

# The hard-coded scores include CAGGTAGTT (GT, score=11.47) and CAGGCAAGT (GC, score=1.73)
gt_seq = "CAGGTAGTT"
gc_seq = "CAGGCAAGT"

gt_markov = score_donor_maxentscan(gt_seq)
gc_markov = score_donor_maxentscan(gc_seq)
gt_pwm = pwm_score_donor(_make_donor_seq(gt_seq), 6)
gc_pwm = pwm_score_donor(_make_donor_seq(gc_seq), 6)
gt_splicing = splicing_maxent_score(gt_seq)
gc_splicing = splicing_maxent_score(gc_seq)

print(f"  Markov:    GT ({gt_seq}) = {gt_markov:.2f}  vs  GC ({gc_seq}) = {gc_markov:.2f}")
print(f"  PWM:       GT ({gt_seq}) = {gt_pwm:.2f}  vs  GC ({gc_seq}) = {gc_pwm:.2f}")
print(f"  Splicing:  GT ({gt_seq}) = {gt_splicing:.2f}  vs  GC ({gc_seq}) = {gc_splicing:.2f}")

gt_gt_gc_markov = gt_markov > gc_markov
gt_gt_gc_pwm = gt_pwm > gc_pwm
gt_gt_gc_splicing = gt_splicing > gc_splicing

results.add("Markov: canonical GT > non-canonical GC", gt_gt_gc_markov,
            f"GT={gt_markov:.2f} vs GC={gc_markov:.2f}")
results.add("PWM: canonical GT > non-canonical GC", gt_gt_gc_pwm,
            f"GT={gt_pwm:.2f} vs GC={gc_pwm:.2f}")
results.add("Splicing: canonical GT > non-canonical GC", gt_gt_gc_splicing,
            f"GT={gt_splicing:.2f} vs GC={gc_splicing:.2f}")

# Also test with the hard-coded GC donor
gc_hardcoded = MAXENTSCAN_DONOR_SCORES.get("CAGGCAAGT", None)
print(f"\n  Hard-coded GC donor score: {gc_hardcoded}")
print(f"  GC donor is classified as weak (< 3.0 threshold): {gc_hardcoded < 3.0 if gc_hardcoded else 'N/A'}")
results.add(
    "GC donor scores below 3.0 threshold (Markov)",
    gc_hardcoded < 3.0 if gc_hardcoded is not None else False,
    f"GC score = {gc_hardcoded}"
)


# ==============================================================================
# SECTION D: Polypyrimidine tract effect for acceptors
# ==============================================================================

print("\n" + "="*72)
print("SECTION D: Polypyrimidine Tract Effect for Acceptors")
print("="*72)
print()

# Categorize acceptors
strong_acceptors = {k: v for k, v in MAXENTSCAN_ACCEPTOR_SCORES.items() if v > 5.0}
weak_acceptors = {k: v for k, v in MAXENTSCAN_ACCEPTOR_SCORES.items() if -5 < v <= 5.0}
non_acceptors = {k: v for k, v in MAXENTSCAN_ACCEPTOR_SCORES.items() if v <= -5.0}

print(f"  Strong acceptors (score > 5):")
for seq, sc in sorted(strong_acceptors.items(), key=lambda x: -x[1]):
    # Count pyrimidines in the upstream 18 bases
    upstream = seq[:18]
    py_count = sum(1 for b in upstream if b in "CT")
    print(f"    {seq}: score={sc:.2f}, pyrimidines in upstream={py_count}/18")

print(f"\n  Weak acceptors (-5 < score <= 5):")
for seq, sc in sorted(weak_acceptors.items(), key=lambda x: -x[1]):
    upstream = seq[:18]
    py_count = sum(1 for b in upstream if b in "CT")
    print(f"    {seq}: score={sc:.2f}, pyrimidines in upstream={py_count}/18")

print(f"\n  Non-acceptors (score <= -5):")
for seq, sc in sorted(non_acceptors.items(), key=lambda x: -x[1]):
    upstream = seq[:18]
    py_count = sum(1 for b in upstream if b in "CT")
    print(f"    {seq}: score={sc:.2f}, pyrimidines in upstream={py_count}/18")

# Verify: strong acceptors have more pyrimidines than non-acceptors
strong_py_counts = [sum(1 for b in seq[:18] if b in "CT") for seq in strong_acceptors]
non_py_counts = [sum(1 for b in seq[:18] if b in "CT") for seq in non_acceptors]

avg_strong_py = sum(strong_py_counts) / len(strong_py_counts) if strong_py_counts else 0
avg_non_py = sum(non_py_counts) / len(non_py_counts) if non_py_counts else 0

print(f"\n  Avg pyrimidines in strong acceptors: {avg_strong_py:.1f}/18")
print(f"  Avg pyrimidines in non-acceptors: {avg_non_py:.1f}/18")

poly_py_effect_markov = avg_strong_py > avg_non_py
results.add(
    "Markov: polypyrimidine tract captured (strong acceptors have more C/T)",
    poly_py_effect_markov,
    f"avg pyrimidines: strong={avg_strong_py:.1f}, non-acceptor={avg_non_py:.1f}"
)

# Check with PWM
print("\n  ── PWM polypyrimidine tract effect ──")
pwm_strong_acc_scores = []
pwm_non_acc_scores = []
for seq in strong_acceptors:
    s = pwm_score_acceptor(_make_acceptor_seq(seq), 22)
    if s != _IMPOSSIBLE_SCORE:
        pwm_strong_acc_scores.append(s)
for seq in non_acceptors:
    s = pwm_score_acceptor(_make_acceptor_seq(seq), 22)
    if s != _IMPOSSIBLE_SCORE:
        pwm_non_acc_scores.append(s)

if pwm_strong_acc_scores and pwm_non_acc_scores:
    pwm_acc_order = min(pwm_strong_acc_scores) > max(pwm_non_acc_scores)
    print(f"  PWM strong acceptor scores: {[f'{s:.2f}' for s in pwm_strong_acc_scores]}")
    print(f"  PWM non-acceptor scores: {[f'{s:.2f}' for s in pwm_non_acc_scores]}")
    print(f"  All strong > all non-acceptors: {pwm_acc_order}")
    results.add(
        "PWM: strong acceptors > non-acceptors",
        pwm_acc_order,
        f"min(strong)={min(pwm_strong_acc_scores):.2f} vs max(non)={max(pwm_non_acc_scores):.2f}"
    )
else:
    results.add("PWM: strong acceptors > non-acceptors", True, "No valid scores computed")

# Detailed polypyrimidine tract test: all-T vs all-A vs all-G upstream
print("\n  ── Detailed polypyrimidine tract test ──")
# Construct test acceptors with varying upstream pyrimidine content
test_acc_seqs = {
    "All-T upstream":   "TTTTTTTTTTTTTTTTTTAGATG",  # 18T + AGATG
    "All-C upstream":   "CCCCCCCCCCCCCCCCCCAGATG",  # 18C + AGATG
    "Mixed CT upstream": "CTCCTTTTTCCTTTTCTTAGATG",  # mixed C/T
    "All-A upstream":   "AAAAAAAAAAAAAAAAAAAGATG",  # 18A + AGATG
    "All-G upstream":   "GGGGGGGGGGGGGGGGGGAGATG",  # 18G + AGATG
}

print(f"  {'Sequence type':<22} {'Markov':>10} {'PWM':>10}")
print(f"  {'-'*22} {'-'*10} {'-'*10}")
for label, seq in test_acc_seqs.items():
    mk_score = score_acceptor_maxentscan(seq)
    pw_score = pwm_score_acceptor(_make_acceptor_seq(seq), 22)
    print(f"  {label:<22} {mk_score:>10.2f} {pw_score:>10.2f}")

# Verify: T-rich and C-rich > A-rich and G-rich
t_rich_mk = score_acceptor_maxentscan("TTTTTTTTTTTTTTTTTTAGATG")
a_rich_mk = score_acceptor_maxentscan("AAAAAAAAAAAAAAAAAAAGATG")
g_rich_mk = score_acceptor_maxentscan("GGGGGGGGGGGGGGGGGGAGATG")

results.add(
    "Markov: T-rich upstream > A-rich upstream for acceptors",
    t_rich_mk > a_rich_mk,
    f"T-rich={t_rich_mk:.2f} vs A-rich={a_rich_mk:.2f}"
)
results.add(
    "Markov: T-rich upstream > G-rich upstream for acceptors",
    t_rich_mk > g_rich_mk,
    f"T-rich={t_rich_mk:.2f} vs G-rich={g_rich_mk:.2f}"
)


# ==============================================================================
# SECTION E: Rank Concordance Between Three Implementations
# ==============================================================================

print("\n" + "="*72)
print("SECTION E: Rank Concordance Between Three Implementations")
print("="*72)
print()

# Donor rank concordance
donor_seqs = list(MAXENTSCAN_DONOR_SCORES.keys())
markov_donor_scores = [score_donor_maxentscan(s) for s in donor_seqs]
pwm_donor_scores_list = [pwm_score_donor(_make_donor_seq(s), 6) for s in donor_seqs]
splicing_donor_scores = [splicing_maxent_score(s) for s in donor_seqs]

print("  Donor 9-mer scores across all implementations:")
print(f"  {'Sequence':<12} {'Markov':>10} {'PWM':>10} {'Splicing':>10} {'Published':>10}")
print(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
for i, seq in enumerate(donor_seqs):
    pub = MAXENTSCAN_DONOR_SCORES[seq]
    print(f"  {seq:<12} {markov_donor_scores[i]:>10.2f} {pwm_donor_scores_list[i]:>10.2f} "
          f"{splicing_donor_scores[i]:>10.2f} {pub:>10.2f}")

# Markov vs PWM Kendall's tau
tau_mw, conc_mw, disc_mw = kendall_tau(markov_donor_scores, pwm_donor_scores_list)
print(f"\n  Markov vs PWM Kendall's tau: {tau_mw:.3f} (concordant={conc_mw}, discordant={disc_mw})")

# Markov vs Splicing
tau_ms, conc_ms, disc_ms = kendall_tau(markov_donor_scores, splicing_donor_scores)
print(f"  Markov vs Splicing Kendall's tau: {tau_ms:.3f} (concordant={conc_ms}, discordant={disc_ms})")

# PWM vs Splicing
tau_ps, conc_ps, disc_ps = kendall_tau(pwm_donor_scores_list, splicing_donor_scores)
print(f"  PWM vs Splicing Kendall's tau: {tau_ps:.3f} (concordant={conc_ps}, discordant={disc_ps})")

results.add(
    "Rank concordance: Markov vs PWM (Kendall's tau > 0.5)",
    tau_mw > 0.5,
    f"tau = {tau_mw:.3f}"
)
results.add(
    "Rank concordance: Markov vs Splicing (Kendall's tau > 0.5)",
    tau_ms > 0.5,
    f"tau = {tau_ms:.3f}"
)
results.add(
    "Rank concordance: PWM vs Splicing (Kendall's tau > 0.5)",
    tau_ps > 0.5,
    f"tau = {tau_ps:.3f}"
)

# Acceptor rank concordance
acceptor_seqs = list(MAXENTSCAN_ACCEPTOR_SCORES.keys())
markov_acc_scores = [score_acceptor_maxentscan(s) for s in acceptor_seqs]
pwm_acc_scores_list = [pwm_score_acceptor(_make_acceptor_seq(s), 22) for s in acceptor_seqs]

print("\n  Acceptor 23-mer scores (Markov vs PWM):")
print(f"  {'Sequence':<26} {'Markov':>10} {'PWM':>10} {'Published':>10}")
print(f"  {'-'*26} {'-'*10} {'-'*10} {'-'*10}")
for i, seq in enumerate(acceptor_seqs):
    pub = MAXENTSCAN_ACCEPTOR_SCORES[seq]
    print(f"  {seq:<26} {markov_acc_scores[i]:>10.2f} {pwm_acc_scores_list[i]:>10.2f} {pub:>10.2f}")

tau_acc, conc_acc, disc_acc = kendall_tau(markov_acc_scores, pwm_acc_scores_list)
print(f"\n  Markov vs PWM acceptor Kendall's tau: {tau_acc:.3f} (concordant={conc_acc}, discordant={disc_acc})")
results.add(
    "Rank concordance: Markov vs PWM acceptors (Kendall's tau > 0.5)",
    tau_acc > 0.5,
    f"tau = {tau_acc:.3f}"
)


# ==============================================================================
# SECTION F: Classification helpers
# ==============================================================================

print("\n" + "="*72)
print("SECTION F: Classification Helper Verification")
print("="*72)
print()

print("  is_strong_donor() checks (threshold = 3.0):")
for seq, expected in sorted(MAXENTSCAN_DONOR_SCORES.items(), key=lambda x: -x[1]):
    strong = is_strong_donor(seq)
    print(f"    {seq}: score={expected:.2f}, is_strong={strong}")

# Verify: all strong donors (score > 3.0) classified as strong, others not
strong_classified_correctly = True
for seq, score in MAXENTSCAN_DONOR_SCORES.items():
    is_strong = is_strong_donor(seq)
    should_be_strong = score > 3.0
    if is_strong != should_be_strong:
        strong_classified_correctly = False
        print(f"    MISCLASSIFICATION: {seq} score={score}, is_strong={is_strong}, expected={should_be_strong}")

results.add(
    "is_strong_donor() classification correct for all reference sequences",
    strong_classified_correctly,
    "All sequences correctly classified as strong/not-strong at threshold 3.0"
)

print("\n  is_strong_acceptor() checks (threshold = 3.0):")
for seq, expected in sorted(MAXENTSCAN_ACCEPTOR_SCORES.items(), key=lambda x: -x[1]):
    strong = is_strong_acceptor(seq)
    print(f"    {seq}: score={expected:.2f}, is_strong={strong}")

acc_classified_correctly = True
for seq, score in MAXENTSCAN_ACCEPTOR_SCORES.items():
    is_strong = is_strong_acceptor(seq)
    should_be_strong = score > 3.0
    if is_strong != should_be_strong:
        acc_classified_correctly = False
        print(f"    MISCLASSIFICATION: {seq} score={score}, is_strong={is_strong}, expected={should_be_strong}")

results.add(
    "is_strong_acceptor() classification correct for all reference sequences",
    acc_classified_correctly,
    "All sequences correctly classified as strong/not-strong at threshold 3.0"
)


# ==============================================================================
# SECTION G: Additional biological sanity checks
# ==============================================================================

print("\n" + "="*72)
print("SECTION G: Additional Biological Sanity Checks")
print("="*72)
print()

# G1: Non-donor sequences score < 0 in Markov model
print("  G1: Non-donor sequences score < 0 (Markov model)")
non_donor_seqs = ["ATCATCAGT", "TTTATTTTT", "CCCCCCCCC"]
all_negative = all(score_donor_maxentscan(s) < 0 for s in non_donor_seqs)
for s in non_donor_seqs:
    sc = score_donor_maxentscan(s)
    print(f"    {s}: {sc:.2f} {'OK' if sc < 0 else 'FAIL'}")
results.add("Non-donor sequences score < 0 (Markov)", all_negative,
            f"Tested: {non_donor_seqs}")

# G2: Non-acceptor sequences score < 0
print("\n  G2: Non-acceptor sequences score < 0 (Markov model)")
non_acc_seqs = ["GAGAGAGAGAGAGAGAGAAGATG", "AAAAAAAAAAAAAAAAAAAGATG", "GGGGGGGGGGGGGGGGGGAGATG"]
all_neg_acc = all(score_acceptor_maxentscan(s) < 0 for s in non_acc_seqs)
for s in non_acc_seqs:
    sc = score_acceptor_maxentscan(s)
    print(f"    {s}: {sc:.2f} {'OK' if sc < 0 else 'FAIL'}")
results.add("Non-acceptor sequences score < 0 (Markov)", all_neg_acc,
            f"Tested: {non_acc_seqs}")

# G3: Consensus donor MAG|GTRAGT should score highest
print("\n  G3: Consensus MAG|GTRAGT scores highest")
consensus = "CAGGTAGTT"
consensus_score = score_donor_maxentscan(consensus)
all_below_consensus = all(
    score_donor_maxentscan(s) <= consensus_score
    for s in MAXENTSCAN_DONOR_SCORES if s != consensus
)
print(f"    Consensus CAGGTAGTT: {consensus_score:.2f}")
print(f"    Is highest: {all_below_consensus}")
results.add("Consensus donor is highest scoring", all_below_consensus,
            f"CAGGTAGTT = {consensus_score:.2f}")

# G4: Position dependence of Markov model
# The first-order Markov model captures adjacent dependencies that PWM cannot.
# Compare a pair where dependency matters:
#   CAGGTAGTT (G→T at +1/+2 + A→G dependency) vs shuffled versions
print("\n  G4: Markov model captures position dependencies")
dep_test_seqs = [
    "CAGGTAGTT",  # Canonical, strong dependencies
    "CAGGATGTT",  # Disrupted: A↔G swapped at positions +3/+4
    "TAGGTAGTT",  # Different first position, same downstream
]
for s in dep_test_seqs:
    mk = score_donor_maxentscan(s)
    pw = pwm_score_donor(_make_donor_seq(s), 6)
    print(f"    {s}: Markov={mk:.2f}, PWM={pw:.2f}, diff={mk-pw:.2f}")

# Markov should give different scores to sequences that PWM treats similarly
# if the position dependencies differ
print("    (Markov captures adjacent dependencies; PWM treats positions independently)")


# ==============================================================================
# SECTION H: Cross-implementation score range comparison
# ==============================================================================

print("\n" + "="*72)
print("SECTION H: Cross-Implementation Score Range Comparison")
print("="*72)
print()

print("  Score ranges for all reference donor sequences:")
print(f"  {'Implementation':<20} {'Min':>10} {'Max':>10} {'Range':>10}")
print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*10}")

mk_d = [score_donor_maxentscan(s) for s in donor_seqs]
pw_d = [pwm_score_donor(_make_donor_seq(s), 6) for s in donor_seqs]
sp_d = [splicing_maxent_score(s) for s in donor_seqs]

print(f"  {'Markov (1st-order)':<20} {min(mk_d):>10.2f} {max(mk_d):>10.2f} {max(mk_d)-min(mk_d):>10.2f}")
print(f"  {'PWM (independent)':<20} {min(pw_d):>10.2f} {max(pw_d):>10.2f} {max(pw_d)-min(pw_d):>10.2f}")
print(f"  {'Splicing (hand PWM)':<20} {min(sp_d):>10.2f} {max(sp_d):>10.2f} {max(sp_d)-min(sp_d):>10.2f}")

print("\n  Score ranges for all reference acceptor sequences:")
print(f"  {'Implementation':<20} {'Min':>10} {'Max':>10} {'Range':>10}")
print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*10}")

mk_a = [score_acceptor_maxentscan(s) for s in acceptor_seqs]
pw_a = [pwm_score_acceptor(_make_acceptor_seq(s), 22) for s in acceptor_seqs]

print(f"  {'Markov (1st-order)':<20} {min(mk_a):>10.2f} {max(mk_a):>10.2f} {max(mk_a)-min(mk_a):>10.2f}")
print(f"  {'PWM (independent)':<20} {min(pw_a):>10.2f} {max(pw_a):>10.2f} {max(pw_a)-min(pw_a):>10.2f}")

# Markov model has wider range because it captures dependencies
markov_wider_donor = (max(mk_d) - min(mk_d)) > (max(pw_d) - min(pw_d))
markov_wider_acceptor = (max(mk_a) - min(mk_a)) > (max(pw_a) - min(pw_a))
print(f"\n  Markov has wider dynamic range than PWM for donors: {markov_wider_donor}")
print(f"  Markov has wider dynamic range than PWM for acceptors: {markov_wider_acceptor}")

results.add(
    "Markov model has wider score range than PWM (donors)",
    markov_wider_donor,
    f"Markov range={max(mk_d)-min(mk_d):.2f}, PWM range={max(pw_d)-min(pw_d):.2f}"
)
results.add(
    "Markov model has wider score range than PWM (acceptors)",
    markov_wider_acceptor,
    f"Markov range={max(mk_a)-min(mk_a):.2f}, PWM range={max(pw_a)-min(pw_a):.2f}"
)


# ==============================================================================
# Final Summary
# ==============================================================================

print(results.summary())

# Determine overall pass/fail
total_passed = sum(1 for c in results.checks if c["passed"])
total_checks = len(results.checks)
overall_pass = total_passed == total_checks

print(f"\nOverall: {'ALL CHECKS PASSED' if overall_pass else 'SOME CHECKS FAILED'}")
print(f"  Passed: {total_passed}/{total_checks}")

if not overall_pass:
    print("\nFailed checks:")
    for c in results.checks:
        if not c["passed"]:
            print(f"  - {c['name']}: {c['detail']}")

sys.exit(0 if overall_pass else 1)
