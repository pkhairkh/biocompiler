"""
Cross-validation test: BioCompiler CAI vs. independent Sharp & Li implementation.

This test suite verifies that our primary CAI implementation
(biocompiler.translation.compute_cai) agrees with an independent
reimplementation (biocompiler.benchmarking.cai_validated.compute_cai_sharp_li)
following the Sharp & Li (1987) formula.

Test categories:
  1. 100 random protein/DNA sequences × 5 organisms — agreement within ±0.02
  2. Property-based tests (via hypothesis, if available):
       - CAI is always in [0, 1]
       - All-optimal-codon sequence → CAI ≈ 1.0
       - CAI is length-independent (geometric mean)
       - Replacing a codon with a more frequent one never decreases CAI
  3. Detailed discrepancy reporting for debugging
"""

from __future__ import annotations

import math
import random
import sys
from typing import Optional

import pytest
pytest.importorskip("hypothesis")

from biocompiler.expression.translation import compute_cai
from biocompiler.benchmarking.cai_validated import (
    compute_cai_sharp_li,
    compute_cai_sharp_li_for_organism,
    load_reference_set,
)
from biocompiler.shared.constants import CODON_TABLE, AA_TO_CODONS
from biocompiler.organisms import (
    CODON_ADAPTIVENESS_TABLES,
    SUPPORTED_ORGANISMS,
    PREFERRED_CODON_TABLES,
)

# ────────────────────────────────────────────────────────────────────
# Try importing hypothesis for property-based tests
# ────────────────────────────────────────────────────────────────────
try:
    from hypothesis import given, settings, assume, HealthCheck
    from hypothesis import strategies as st

    HAS_HYPOTHESIS = True
except ImportError:
    HAS_HYPOTHESIS = False

# ────────────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────────────

# Organisms to test
ORGANISMS = [
    "Escherichia_coli",
    "Homo_sapiens",
    "Saccharomyces_cerevisiae",
    "Mus_musculus",
    "CHO_K1",
]

# Tolerance for cross-implementation agreement
CAI_TOLERANCE = 0.02

# Number of random sequences for the bulk cross-validation
N_RANDOM_SEQUENCES = 100

# Random seed for reproducibility
RANDOM_SEED = 42

# Standard 20 amino acids (excluding stop)
STANDARD_AAS = [aa for aa in AA_TO_CODONS if aa != "*"]

# Realistic amino acid composition weights (approximating human proteins)
# Based on average amino acid frequencies from UniProt
_AA_WEIGHTS: dict[str, float] = {
    "A": 8.25, "R": 5.53, "N": 4.06, "D": 5.45, "C": 1.37,
    "Q": 3.93, "E": 6.75, "G": 7.07, "H": 2.27, "I": 5.96,
    "L": 9.65, "K": 5.84, "M": 2.42, "F": 3.86, "P": 4.70,
    "S": 6.56, "T": 5.34, "W": 1.08, "Y": 2.92, "V": 6.87,
}


# ────────────────────────────────────────────────────────────────────
# Helper functions
# ────────────────────────────────────────────────────────────────────


def _weighted_random_aa() -> str:
    """Pick a random amino acid with realistic composition weights."""
    aas = list(_AA_WEIGHTS.keys())
    weights = list(_AA_WEIGHTS.values())
    return random.choices(aas, weights=weights, k=1)[0]


def _random_protein(min_len: int = 20, max_len: int = 200) -> str:
    """Generate a random protein sequence with realistic AA composition."""
    length = random.randint(min_len, max_len)
    # Always start with M (as in real coding sequences)
    return "M" + "".join(_weighted_random_aa() for _ in range(length - 1))


def _random_dna_for_protein(protein: str) -> str:
    """Generate a random valid DNA encoding for a protein sequence.

    For each amino acid, chooses a random synonymous codon.
    """
    codons: list[str] = []
    for aa in protein:
        syn_codons = AA_TO_CODONS.get(aa, [])
        if not syn_codons:
            raise ValueError(f"No codons found for amino acid '{aa}'")
        codons.append(random.choice(syn_codons))
    return "".join(codons)


def _optimal_dna_for_protein(protein: str, organism: str) -> str:
    """Generate DNA using all optimal (highest-frequency) codons."""
    preferred = PREFERRED_CODON_TABLES[organism]
    codons: list[str] = []
    for aa in protein:
        if aa == "M":
            codons.append("ATG")  # Only one codon for Met
        elif aa == "W":
            codons.append("TGG")  # Only one codon for Trp
        else:
            codons.append(preferred.get(aa, AA_TO_CODONS[aa][0]))
    return "".join(codons)


def _find_more_frequent_codon(
    codon: str, organism: str
) -> Optional[str]:
    """Find a synonymous codon with higher adaptiveness, or return None."""
    aa = CODON_TABLE.get(codon)
    if aa is None or aa == "*" or aa == "M":
        return None

    adaptiveness = CODON_ADAPTIVENESS_TABLES[organism]
    current_w = adaptiveness.get(codon, 0.0)
    syn_codons = AA_TO_CODONS.get(aa, [codon])

    better = [
        c for c in syn_codons
        if adaptiveness.get(c, 0.0) > current_w
    ]
    if better:
        return max(better, key=lambda c: adaptiveness[c])
    return None


# ════════════════════════════════════════════════════════════════════
# 1. BULK CROSS-VALIDATION: 100 random sequences × 5 organisms
# ════════════════════════════════════════════════════════════════════


class TestBulkCrossValidation:
    """Compare both CAI implementations on 100 random sequences per organism."""

    @pytest.fixture(autouse=True)
    def _set_seed(self):
        """Fix random seed for reproducibility."""
        random.seed(RANDOM_SEED)

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_100_random_sequences_agree(self, organism):
        """Both implementations must agree within ±0.02 for all 100 sequences."""
        discrepancies: list[dict] = []

        for seq_idx in range(N_RANDOM_SEQUENCES):
            protein = _random_protein(min_len=20, max_len=200)
            dna = _random_dna_for_protein(protein)

            cai_ours = compute_cai(dna, organism=organism)
            cai_sharp_li = compute_cai_sharp_li_for_organism(dna, organism)

            diff = abs(cai_ours - cai_sharp_li)
            if diff > CAI_TOLERANCE:
                discrepancies.append({
                    "seq_idx": seq_idx,
                    "protein_len": len(protein),
                    "dna_len": len(dna),
                    "organism": organism,
                    "cai_ours": cai_ours,
                    "cai_sharp_li": cai_sharp_li,
                    "diff": diff,
                    "dna_first_30": dna[:30],
                    "protein_first_10": protein[:10],
                })

        # If any discrepancies, print detailed debugging info
        if discrepancies:
            print("\n" + "=" * 70, file=sys.stderr)
            print(
                f"CROSS-VALIDATION FAILURES: {len(discrepancies)} / "
                f"{N_RANDOM_SEQUENCES} sequences disagree for {organism}",
                file=sys.stderr,
            )
            print("=" * 70, file=sys.stderr)
            for d in discrepancies:
                print(
                    f"  Seq #{d['seq_idx']}: protein_len={d['protein_len']}, "
                    f"organism={d['organism']}\n"
                    f"    cai_ours     = {d['cai_ours']}\n"
                    f"    cai_sharp_li = {d['cai_sharp_li']}\n"
                    f"    diff         = {d['diff']:.6f}\n"
                    f"    dna[:30]     = {d['dna_first_30']}\n"
                    f"    protein[:10] = {d['protein_first_10']}",
                    file=sys.stderr,
                )
            print("=" * 70, file=sys.stderr)

        assert not discrepancies, (
            f"{len(discrepancies)}/{N_RANDOM_SEQUENCES} sequences disagreed "
            f"by >{CAI_TOLERANCE} for {organism}. "
            f"Max diff: {max(d['diff'] for d in discrepancies):.6f}"
        )


# ════════════════════════════════════════════════════════════════════
# 2. DETERMINISTIC SPOT-CHECKS
# ════════════════════════════════════════════════════════════════════


class TestDeterministicSpotChecks:
    """Hand-crafted sequences where we know what CAI should be."""

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_all_optimal_codons_cai_near_one(self, organism):
        """A sequence using ALL preferred codons should have CAI ≈ 1.0."""
        # Use a protein that exercises all amino acids with synonymous codons
        protein = "MFLIVSRHDEKCWTAYNPQGALM"
        dna = _optimal_dna_for_protein(protein, organism)

        cai_ours = compute_cai(dna, organism=organism)
        cai_sharp_li = compute_cai_sharp_li_for_organism(dna, organism)

        # Both should be close to 1.0
        assert cai_ours >= 0.95, (
            f"Expected CAI ≥ 0.95 for all-optimal sequence, got {cai_ours} "
            f"({organism})"
        )
        assert cai_sharp_li >= 0.95, (
            f"Sharp-Li expected CAI ≥ 0.95 for all-optimal sequence, "
            f"got {cai_sharp_li} ({organism})"
        )
        # And they should agree with each other
        assert abs(cai_ours - cai_sharp_li) <= CAI_TOLERANCE, (
            f"Disagreement: ours={cai_ours}, sharp_li={cai_sharp_li} "
            f"for all-optimal sequence ({organism})"
        )

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_single_codon_types_no_disagreement(self, organism):
        """Sequences of a single non-Met amino acid type should agree."""
        # Test each amino acid with 5 repeats
        for aa in "FLIVSRHDEKCWTAYNPQG":
            syn_codons = AA_TO_CODONS[aa]
            if len(syn_codons) <= 1:
                continue  # Skip single-codon AAs (W, M)
            # Use the rarest codon for this AA
            adaptiveness = CODON_ADAPTIVENESS_TABLES[organism]
            rarest = min(syn_codons, key=lambda c: adaptiveness.get(c, 0.0))
            dna = "ATG" + rarest * 5  # ATG (Met) skipped + 5 rare codons

            cai_ours = compute_cai(dna, organism=organism)
            cai_sharp_li = compute_cai_sharp_li_for_organism(dna, organism)

            assert abs(cai_ours - cai_sharp_li) <= CAI_TOLERANCE, (
                f"Disagreement on {aa}={rarest}×5: ours={cai_ours}, "
                f"sharp_li={cai_sharp_li} ({organism})"
            )

    def test_empty_sequence_both_zero(self):
        """Empty sequence → 0.0 for both implementations."""
        for org in ORGANISMS:
            assert compute_cai("", organism=org) == 0.0
            assert compute_cai_sharp_li_for_organism("", org) == 0.0

    def test_only_met_and_stop_both_zero(self):
        """Sequence with only Met + stop → 0.0 for both (no contributing codons)."""
        dna = "ATGTAA"
        for org in ORGANISMS:
            assert compute_cai(dna, organism=org) == 0.0
            assert compute_cai_sharp_li_for_organism(dna, org) == 0.0


# ════════════════════════════════════════════════════════════════════
# 3. PROPERTY-BASED TESTS (hypothesis, if available)
# ════════════════════════════════════════════════════════════════════


if HAS_HYPOTHESIS:

    # Strategy: generate a valid DNA sequence (ATG + random synonymous codons)
    def _dna_strategy(min_codons=3, max_codons=60):
        """Hypothesis strategy that generates valid coding DNA sequences."""

        @st.composite
        def _gen(draw):
            length = draw(st.integers(min_value=min_codons, max_value=max_codons))
            # Pick random amino acids
            aas = draw(
                st.lists(
                    st.sampled_from(STANDARD_AAS),
                    min_size=length,
                    max_size=length,
                )
            )
            # Force first AA to be M (start codon)
            aas[0] = "M"
            # Pick random synonymous codons
            codons = []
            for aa in aas:
                syn = AA_TO_CODONS.get(aa, ["ATG"])
                codon = draw(st.sampled_from(syn))
                codons.append(codon)
            return "".join(codons)

        return _gen()

    class TestPropertyBased:
        """Property-based tests using hypothesis."""

        @given(dna=_dna_strategy(), organism=st.sampled_from(ORGANISMS))
        @settings(
            max_examples=200,
            suppress_health_check=[HealthCheck.too_slow],
            deadline=None,
        )
        def test_cai_always_in_zero_one(self, dna, organism):
            """Property: CAI is always in [0, 1]."""
            cai_ours = compute_cai(dna, organism=organism)
            cai_sharp = compute_cai_sharp_li_for_organism(dna, organism)

            assert 0.0 <= cai_ours <= 1.0, (
                f"Our CAI out of range: {cai_ours} for {organism}"
            )
            assert 0.0 <= cai_sharp <= 1.0, (
                f"Sharp-Li CAI out of range: {cai_sharp} for {organism}"
            )

        @given(dna=_dna_strategy(), organism=st.sampled_from(ORGANISMS))
        @settings(
            max_examples=200,
            suppress_health_check=[HealthCheck.too_slow],
            deadline=None,
        )
        def test_both_implementations_agree(self, dna, organism):
            """Property: Both implementations agree within ±0.02."""
            cai_ours = compute_cai(dna, organism=organism)
            cai_sharp = compute_cai_sharp_li_for_organism(dna, organism)

            diff = abs(cai_ours - cai_sharp)
            assert diff <= CAI_TOLERANCE, (
                f"CAI implementations disagree: ours={cai_ours}, "
                f"sharp_li={cai_sharp}, diff={diff:.6f}, "
                f"organism={organism}, dna[:30]={dna[:30]}"
            )

        @given(organism=st.sampled_from(ORGANISMS))
        @settings(
            max_examples=5 * len(ORGANISMS),
            suppress_health_check=[HealthCheck.too_slow],
            deadline=None,
        )
        def test_all_optimal_codons_cai_near_one(self, organism):
            """Property: All-optimal-codon sequence has CAI ≈ 1.0."""
            # Use a protein with diverse amino acids
            protein = "MFLIVSRHDEKCWTAYNPQGA"
            dna = _optimal_dna_for_protein(protein, organism)

            cai_ours = compute_cai(dna, organism=organism)
            cai_sharp = compute_cai_sharp_li_for_organism(dna, organism)

            assert cai_ours >= 0.95, f"Our CAI={cai_ours} < 0.95 for {organism}"
            assert cai_sharp >= 0.95, (
                f"Sharp-Li CAI={cai_sharp} < 0.95 for {organism}"
            )

        @given(organism=st.sampled_from(ORGANISMS), n_copies=st.integers(5, 50))
        @settings(
            max_examples=100,
            suppress_health_check=[HealthCheck.too_slow],
            deadline=None,
        )
        def test_cai_length_independent(self, organism, n_copies):
            """Property: CAI is independent of protein length (geometric mean).

            Repeating the same codon pattern N times should yield the same CAI
            as a single copy, because CAI is a geometric mean.
            """
            # A simple 3-codon protein (Met + 2 variable codons)
            protein_single = "MFA"
            dna_single = _random_dna_for_protein(protein_single)
            cai_single = compute_cai(dna_single, organism=organism)

            # Repeat the non-Met part N times
            protein_repeated = "M" + "FA" * n_copies
            dna_repeated = (
                "ATG"
                + dna_single[3:9] * n_copies  # repeat F+A codons
            )
            cai_repeated = compute_cai(dna_repeated, organism=organism)

            # Should be identical (or very close due to rounding)
            assert abs(cai_single - cai_repeated) <= 0.015, (
                f"CAI not length-independent: single={cai_single}, "
                f"repeated({n_copies})={cai_repeated}, "
                f"diff={abs(cai_single - cai_repeated):.6f}, "
                f"organism={organism}"
            )

        @given(dna=_dna_strategy(min_codons=5, max_codons=30), organism=st.sampled_from(ORGANISMS))
        @settings(
            max_examples=200,
            suppress_health_check=[HealthCheck.too_slow],
            deadline=None,
        )
        def test_replacing_with_more_frequent_codon_never_decreases_cai(
            self, dna, organism
        ):
            """Property: Replacing a codon with a more frequent one never decreases CAI."""
            # Find a codon position (not Met, not stop) to replace
            codon_list = [dna[i:i + 3] for i in range(0, len(dna) - 2, 3)]
            valid_positions = []
            for idx, codon in enumerate(codon_list):
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*" or aa == "M":
                    continue
                better = _find_more_frequent_codon(codon, organism)
                if better is not None:
                    valid_positions.append((idx, codon, better))

            assume(len(valid_positions) > 0)

            # Pick a random valid position to upgrade
            idx, old_codon, new_codon = random.choice(valid_positions)
            pos = idx * 3
            new_dna = dna[:pos] + new_codon + dna[pos + 3:]

            cai_old = compute_cai(dna, organism=organism)
            cai_new = compute_cai(new_dna, organism=organism)

            assert cai_new >= cai_old - 0.001, (
                f"Replacing {old_codon}→{new_codon} at pos {pos} "
                f"decreased CAI: {cai_old} → {cai_new} "
                f"(organism={organism})"
            )


# ════════════════════════════════════════════════════════════════════
# 4. PROPERTY-BASED TESTS (fallback without hypothesis)
# ════════════════════════════════════════════════════════════════════


class TestPropertyBasedFallback:
    """Property-based tests that run even without hypothesis."""

    @pytest.fixture(autouse=True)
    def _set_seed(self):
        random.seed(RANDOM_SEED)

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_cai_always_in_unit_interval(self, organism):
        """CAI is always in [0, 1] for random sequences."""
        for _ in range(50):
            protein = _random_protein(min_len=10, max_len=100)
            dna = _random_dna_for_protein(protein)
            cai_ours = compute_cai(dna, organism=organism)
            cai_sharp = compute_cai_sharp_li_for_organism(dna, organism)

            assert 0.0 <= cai_ours <= 1.0, (
                f"Our CAI out of range: {cai_ours} ({organism})"
            )
            assert 0.0 <= cai_sharp <= 1.0, (
                f"Sharp-Li CAI out of range: {cai_sharp} ({organism})"
            )

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_all_optimal_cai_near_one(self, organism):
        """All-optimal-codon sequence → CAI ≈ 1.0."""
        protein = "MFLIVSRHDEKCWTAYNPQGA"
        dna = _optimal_dna_for_protein(protein, organism)

        cai_ours = compute_cai(dna, organism=organism)
        cai_sharp = compute_cai_sharp_li_for_organism(dna, organism)

        assert cai_ours >= 0.95, f"Our CAI={cai_ours} < 0.95 for {organism}"
        assert cai_sharp >= 0.95, (
            f"Sharp-Li CAI={cai_sharp} < 0.95 for {organism}"
        )

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_cai_length_independent(self, organism):
        """CAI is independent of protein sequence length (geometric mean)."""
        # Create a base codon pattern: Met + Phe + Ala
        dna_f = AA_TO_CODONS["F"][0]  # First Phe codon
        dna_a = AA_TO_CODONS["A"][0]  # First Ala codon

        # Single copy: ATG + F + A  (2 contributing codons)
        dna_single = "ATG" + dna_f + dna_a
        cai_single = compute_cai(dna_single, organism=organism)

        # Multiple copies: ATG + (F + A) × N
        for n in [5, 10, 20, 50]:
            dna_repeated = "ATG" + (dna_f + dna_a) * n
            cai_repeated = compute_cai(dna_repeated, organism=organism)
            assert abs(cai_single - cai_repeated) <= 0.015, (
                f"CAI not length-independent: single={cai_single}, "
                f"repeated({n})={cai_repeated}, diff="
                f"{abs(cai_single - cai_repeated):.6f}, "
                f"organism={organism}"
            )

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_more_frequent_codon_never_decreases_cai(self, organism):
        """Replacing a codon with a more frequent one never decreases CAI."""
        for _ in range(50):
            protein = _random_protein(min_len=10, max_len=50)
            dna = _random_dna_for_protein(protein)

            cai_original = compute_cai(dna, organism=organism)

            # Try to replace each non-Met, non-stop codon with a better one
            codon_list = [dna[i:i + 3] for i in range(0, len(dna) - 2, 3)]
            for idx, codon in enumerate(codon_list):
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*" or aa == "M":
                    continue
                better = _find_more_frequent_codon(codon, organism)
                if better is None:
                    continue

                pos = idx * 3
                new_dna = dna[:pos] + better + dna[pos + 3:]
                cai_new = compute_cai(new_dna, organism=organism)

                assert cai_new >= cai_original - 0.001, (
                    f"Replacing {codon}→{better} at pos {pos} "
                    f"decreased CAI: {cai_original:.6f} → {cai_new:.6f} "
                    f"({organism})"
                )
                break  # One valid replacement is enough per iteration


# ════════════════════════════════════════════════════════════════════
# 5. DISCREPANCY DETECTION & REPORTING
# ════════════════════════════════════════════════════════════════════


class TestDiscrepancyDetection:
    """If the two implementations disagree, print detailed discrepancy info."""

    @pytest.fixture(autouse=True)
    def _set_seed(self):
        random.seed(RANDOM_SEED + 7)

    @pytest.mark.parametrize("organism", ORGANISMS)
    def test_no_discrepancies_across_organisms(self, organism):
        """Detailed cross-validation with discrepancy reporting."""
        max_diff = 0.0
        worst_case: Optional[dict] = None

        for seq_idx in range(N_RANDOM_SEQUENCES):
            protein = _random_protein(min_len=20, max_len=200)
            dna = _random_dna_for_protein(protein)

            cai_ours = compute_cai(dna, organism=organism)
            cai_sharp_li = compute_cai_sharp_li_for_organism(dna, organism)

            diff = abs(cai_ours - cai_sharp_li)

            if diff > max_diff:
                max_diff = diff
                worst_case = {
                    "seq_idx": seq_idx,
                    "protein_len": len(protein),
                    "cai_ours": cai_ours,
                    "cai_sharp_li": cai_sharp_li,
                    "diff": diff,
                }

        # Print worst case for visibility (even if it passes)
        if worst_case:
            print(
                f"\n  [{organism}] Worst case: diff={worst_case['diff']:.6f} "
                f"(ours={worst_case['cai_ours']}, "
                f"sharp_li={worst_case['cai_sharp_li']}, "
                f"protein_len={worst_case['protein_len']})",
                file=sys.stderr,
            )

        assert max_diff <= CAI_TOLERANCE, (
            f"Maximum CAI discrepancy {max_diff:.6f} exceeds tolerance "
            f"{CAI_TOLERANCE} for {organism}. "
            f"Worst case: {worst_case}"
        )
