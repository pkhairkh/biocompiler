"""Unit tests for individual organism codon usage modules.

Covers:
1. Each organism module (e_coli, human, yeast, mouse, cho) exports codon_usage dict
2. Codon usage fractions sum to ~1.0 per amino acid
3. All 64 codons represented in each table
4. GC content of codon usage matches organism expectations
"""

from __future__ import annotations

import itertools
from typing import NamedTuple

import pytest

from biocompiler.organisms.escherichia import E_COLI_CODON_USAGE
from biocompiler.organisms.human import HUMAN_CODON_USAGE
from biocompiler.organisms.mouse import MOUSE_CODON_USAGE
from biocompiler.organisms.yeast import YEAST_CODON_USAGE
from biocompiler.organisms.cho import CHO_CODON_USAGE
from biocompiler.organisms import ORGANISM_GC_TARGETS
from biocompiler.type_system import CODON_TABLE

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# All 64 standard codons from the standard genetic code
ALL_STANDARD_CODONS: set[str] = set(CODON_TABLE.keys())
assert len(ALL_STANDARD_CODONS) == 64

# Standard amino acid single-letter codes + stop
STANDARD_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY*")

# Organic bases for generating all possible codons (validation cross-check)
BASES = "ACGT"


class OrganismSpec(NamedTuple):
    """Metadata for a single organism under test."""

    name: str  # human-readable name
    canonical: str  # canonical scientific name used in ORGANISM_GC_TARGETS
    codon_usage: dict[str, tuple[str, float, float, int]]


ORGANISMS: list[OrganismSpec] = [
    OrganismSpec("E. coli", "Escherichia_coli", E_COLI_CODON_USAGE),
    OrganismSpec("Human", "Homo_sapiens", HUMAN_CODON_USAGE),
    OrganismSpec("Yeast", "Saccharomyces_cerevisiae", YEAST_CODON_USAGE),
    OrganismSpec("Mouse", "Mus_musculus", MOUSE_CODON_USAGE),
    OrganismSpec("CHO", "CHO_K1", CHO_CODON_USAGE),
]

# Map from canonical name → OrganismSpec for easy lookup
ORGANISM_BY_CANONICAL: dict[str, OrganismSpec] = {
    o.canonical: o for o in ORGANISMS
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _gc_count(codon: str) -> int:
    """Return the number of G/C bases in a codon (0–3)."""
    return sum(1 for base in codon if base in "GC")


def _compute_gc_from_usage(
    usage: dict[str, tuple[str, float, float, int]],
) -> float:
    """Compute overall GC content from a codon usage table.

    Uses the *per-thousand* frequency (index 2 of the tuple) as the weight
    for each codon, which reflects the actual abundance of that codon in
    the organism's coding sequences.

    Returns:
        GC fraction in [0, 1].
    """
    total_weighted_gc = 0.0
    total_bases = 0.0
    for codon, (_aa, _frac, per_thousand, _count) in usage.items():
        gc = _gc_count(codon)
        total_weighted_gc += gc * per_thousand
        total_bases += 3.0 * per_thousand
    if total_bases == 0:
        return 0.0
    return total_weighted_gc / total_bases


def _group_fractions_by_aa(
    usage: dict[str, tuple[str, float, float, int]],
) -> dict[str, list[float]]:
    """Group codon fractions by amino acid."""
    aa_fractions: dict[str, list[float]] = {}
    for _codon, (aa, frac, *_rest) in usage.items():
        aa_fractions.setdefault(aa, []).append(frac)
    return aa_fractions


# ===========================================================================
# 1. Each organism module exports codon_usage dict
# ===========================================================================


class TestCodonUsageExport:
    """Each organism module must export a properly structured codon_usage dict."""

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_codon_usage_is_dict(self, org: OrganismSpec) -> None:
        """The codon_usage attribute must be a dict."""
        assert isinstance(org.codon_usage, dict), (
            f"{org.name} codon_usage is not a dict"
        )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_codon_usage_non_empty(self, org: OrganismSpec) -> None:
        """The codon_usage dict must be non-empty."""
        assert len(org.codon_usage) > 0, f"{org.name} codon_usage is empty"

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_codon_usage_tuple_values(self, org: OrganismSpec) -> None:
        """Every value in the codon_usage dict must be a 4-tuple
        (amino_acid, fraction, per_thousand, count)."""
        for codon, value in org.codon_usage.items():
            assert isinstance(value, tuple) and len(value) == 4, (
                f"{org.name} codon {codon!r} value is not a 4-tuple: {value!r}"
            )
            aa, frac, per_thousand, count = value
            assert isinstance(aa, str) and len(aa) == 1, (
                f"{org.name} codon {codon!r}: amino acid {aa!r} is not a single char"
            )
            assert isinstance(frac, float), (
                f"{org.name} codon {codon!r}: fraction {frac!r} is not float"
            )
            assert isinstance(per_thousand, float), (
                f"{org.name} codon {codon!r}: per_thousand {per_thousand!r} is not float"
            )
            assert isinstance(count, int), (
                f"{org.name} codon {codon!r}: count {count!r} is not int"
            )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_codon_usage_keys_are_valid_codons(self, org: OrganismSpec) -> None:
        """Every key in the codon_usage dict must be a valid 3-base codon."""
        for codon in org.codon_usage:
            assert len(codon) == 3, (
                f"{org.name} key {codon!r} is not 3 bases long"
            )
            assert all(base in BASES for base in codon), (
                f"{org.name} key {codon!r} contains invalid bases"
            )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_amino_acids_are_standard(self, org: OrganismSpec) -> None:
        """All amino acid annotations must be standard single-letter codes or '*'."""
        for codon, (aa, *_) in org.codon_usage.items():
            assert aa in STANDARD_AA, (
                f"{org.name} codon {codon!r} maps to unknown AA {aa!r}"
            )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_fractions_non_negative(self, org: OrganismSpec) -> None:
        """All codon fractions must be >= 0."""
        for codon, (aa, frac, *_) in org.codon_usage.items():
            assert frac >= 0.0, (
                f"{org.name} codon {codon!r} ({aa}) has negative fraction {frac}"
            )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_per_thousand_non_negative(self, org: OrganismSpec) -> None:
        """All per-thousand values must be >= 0."""
        for codon, (aa, _frac, per_thousand, _) in org.codon_usage.items():
            assert per_thousand >= 0.0, (
                f"{org.name} codon {codon!r} ({aa}) has negative per_thousand {per_thousand}"
            )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_counts_non_negative(self, org: OrganismSpec) -> None:
        """All count values must be >= 0."""
        for codon, (aa, _frac, _per_thousand, count) in org.codon_usage.items():
            assert count >= 0, (
                f"{org.name} codon {codon!r} ({aa}) has negative count {count}"
            )


# ===========================================================================
# 2. Codon usage fractions sum to ~1.0 per amino acid
# ===========================================================================


class TestCodonFractionsSum:
    """For each amino acid, codon fractions should sum to approximately 1.0."""

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_fractions_sum_approx_one_per_aa(self, org: OrganismSpec) -> None:
        """For each amino acid, the sum of fractions across synonymous codons
        must be approximately 1.0.

        Some rounding imprecision exists in source data (especially E. coli
        Valine), so we use a generous tolerance of ±0.20.
        """
        aa_fractions = _group_fractions_by_aa(org.codon_usage)
        for aa, fracs in aa_fractions.items():
            total = sum(fracs)
            assert total == pytest.approx(1.0, abs=0.20), (
                f"{org.name} AA {aa!r} fractions sum to {total:.4f}, "
                f"expected ~1.0 (fractions: {fracs})"
            )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_single_codon_aas_have_fraction_one(self, org: OrganismSpec) -> None:
        """Amino acids with only one codon (Methionine, Tryptophan) must
        have fraction exactly 1.0."""
        aa_fractions = _group_fractions_by_aa(org.codon_usage)
        for aa in ("M", "W"):
            assert aa in aa_fractions, (
                f"{org.name} is missing amino acid {aa!r}"
            )
            fracs = aa_fractions[aa]
            assert len(fracs) == 1, (
                f"{org.name} AA {aa!r} should have exactly 1 codon, "
                f"found {len(fracs)}"
            )
            assert fracs[0] == pytest.approx(1.0, abs=0.01), (
                f"{org.name} AA {aa!r} single-codon fraction should be 1.0, "
                f"got {fracs[0]}"
            )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_stop_codon_fractions_sum_approx_one(self, org: OrganismSpec) -> None:
        """Stop codon fractions should also sum to approximately 1.0."""
        aa_fractions = _group_fractions_by_aa(org.codon_usage)
        assert "*" in aa_fractions, f"{org.name} is missing stop codons"
        total = sum(aa_fractions["*"])
        assert total == pytest.approx(1.0, abs=0.20), (
            f"{org.name} stop codon fractions sum to {total:.4f}, expected ~1.0"
        )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_no_fraction_exceeds_one(self, org: OrganismSpec) -> None:
        """No single codon fraction should exceed 1.0."""
        for codon, (aa, frac, *_) in org.codon_usage.items():
            assert frac <= 1.01, (  # tiny float tolerance
                f"{org.name} codon {codon!r} ({aa}) fraction {frac} exceeds 1.0"
            )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_no_zero_fraction_for_multi_codon_aa(self, org: OrganismSpec) -> None:
        """For amino acids with multiple codons, no fraction should be exactly
        0.0 (that would mean the codon is never used, which is biologically
        unlikely for standard coding sequences)."""
        aa_fractions = _group_fractions_by_aa(org.codon_usage)
        for aa, fracs in aa_fractions.items():
            if len(fracs) > 1:
                for i, frac in enumerate(fracs):
                    assert frac > 0.0, (
                        f"{org.name} AA {aa!r} codon index {i} has fraction 0.0"
                    )


# ===========================================================================
# 3. All 64 codons represented in each table
# ===========================================================================


class TestAllCodonsRepresented:
    """Every organism codon usage table must contain all 64 standard codons."""

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_exact_64_codons(self, org: OrganismSpec) -> None:
        """Each codon usage table must contain exactly 64 entries."""
        assert len(org.codon_usage) == 64, (
            f"{org.name} has {len(org.codon_usage)} codons, expected 64"
        )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_all_standard_codons_present(self, org: OrganismSpec) -> None:
        """Every standard codon must be present in the table."""
        table_codons = set(org.codon_usage.keys())
        missing = ALL_STANDARD_CODONS - table_codons
        assert not missing, (
            f"{org.name} missing codons: {sorted(missing)}"
        )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_no_extra_codons(self, org: OrganismSpec) -> None:
        """No non-standard codons should be present in the table."""
        table_codons = set(org.codon_usage.keys())
        extra = table_codons - ALL_STANDARD_CODONS
        assert not extra, (
            f"{org.name} has extra codons: {sorted(extra)}"
        )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_codon_to_aa_matches_standard_code(self, org: OrganismSpec) -> None:
        """Each codon's amino acid annotation must match the standard genetic code."""
        for codon, (aa, *_) in org.codon_usage.items():
            expected_aa = CODON_TABLE[codon]
            assert aa == expected_aa, (
                f"{org.name} codon {codon!r} maps to {aa!r}, "
                f"but standard code says {expected_aa!r}"
            )

    def test_all_64_codons_generated_independently(self) -> None:
        """Cross-check: generate all 64 codons from base permutations and
        verify they match CODON_TABLE keys."""
        generated = {a + b + c for a, b, c in itertools.product(BASES, repeat=3)}
        assert generated == ALL_STANDARD_CODONS
        assert len(generated) == 64


# ===========================================================================
# 4. GC content of codon usage matches organism expectations
# ===========================================================================


class TestGCContent:
    """GC content computed from codon usage tables should match organism
    expectations (i.e., fall within the ORGANISM_GC_TARGETS range)."""

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_gc_content_within_target_range(self, org: OrganismSpec) -> None:
        """Overall GC content computed from per-thousand frequencies should
        fall within the organism's GC target range."""
        gc = _compute_gc_from_usage(org.codon_usage)
        gc_lo, gc_hi = ORGANISM_GC_TARGETS[org.canonical]
        assert gc_lo <= gc <= gc_hi, (
            f"{org.name} GC content {gc:.4f} outside target range "
            f"[{gc_lo}, {gc_hi}]"
        )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_gc_content_is_biologically_plausible(self, org: OrganismSpec) -> None:
        """GC content must be in the biologically plausible range [0.20, 0.75]."""
        gc = _compute_gc_from_usage(org.codon_usage)
        assert 0.20 <= gc <= 0.75, (
            f"{org.name} GC content {gc:.4f} is outside biologically "
            f"plausible range [0.20, 0.75]"
        )

    def test_ecoli_gc_higher_than_yeast(self) -> None:
        """E. coli coding GC (~50.8%) should be higher than yeast (~38.3%)."""
        ecoli_gc = _compute_gc_from_usage(E_COLI_CODON_USAGE)
        yeast_gc = _compute_gc_from_usage(YEAST_CODON_USAGE)
        assert ecoli_gc > yeast_gc, (
            f"E. coli GC ({ecoli_gc:.4f}) should be > yeast GC ({yeast_gc:.4f})"
        )

    def test_mammalian_gc_similar(self) -> None:
        """Human, mouse, and CHO are all mammalian; their GC contents
        should be within 10 percentage points of each other."""
        human_gc = _compute_gc_from_usage(HUMAN_CODON_USAGE)
        mouse_gc = _compute_gc_from_usage(MOUSE_CODON_USAGE)
        cho_gc = _compute_gc_from_usage(CHO_CODON_USAGE)
        max_diff = max(human_gc, mouse_gc, cho_gc) - min(human_gc, mouse_gc, cho_gc)
        assert max_diff < 0.10, (
            f"Mammalian GC spread too large: human={human_gc:.4f}, "
            f"mouse={mouse_gc:.4f}, CHO={cho_gc:.4f}, "
            f"max_diff={max_diff:.4f}"
        )

    def test_yeast_gc_lowest(self) -> None:
        """Yeast should have the lowest coding GC among all organisms."""
        yeast_gc = _compute_gc_from_usage(YEAST_CODON_USAGE)
        for org in ORGANISMS:
            if org.name == "Yeast":
                continue
            other_gc = _compute_gc_from_usage(org.codon_usage)
            assert other_gc > yeast_gc, (
                f"Yeast GC ({yeast_gc:.4f}) should be lower than "
                f"{org.name} GC ({other_gc:.4f})"
            )

    @pytest.mark.parametrize("org", ORGANISMS, ids=lambda o: o.name)
    def test_gc_computed_from_count_matches_per_thousand(self, org: OrganismSpec) -> None:
        """GC computed from raw counts should closely match GC from
        per-thousand values (they should be proportional)."""
        # Compute GC from counts
        total_gc_bases = 0
        total_bases = 0
        for codon, (_aa, _frac, _per_thousand, count) in org.codon_usage.items():
            gc = _gc_count(codon)
            total_gc_bases += gc * count
            total_bases += 3 * count
        gc_from_counts = total_gc_bases / total_bases if total_bases > 0 else 0.0

        gc_from_per_thousand = _compute_gc_from_usage(org.codon_usage)

        assert gc_from_counts == pytest.approx(gc_from_per_thousand, abs=0.02), (
            f"{org.name} GC from counts ({gc_from_counts:.4f}) disagrees "
            f"with GC from per_thousand ({gc_from_per_thousand:.4f})"
        )
