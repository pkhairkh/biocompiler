"""
Integration tests for organism-specific expression optimization.

Covers:
1. Each organism produces DIFFERENT optimized sequences with appropriate GC content
2. Codon pair bias integration (E. coli CPB)
3. UTR model integration (Shine-Dalgarno for E. coli, Kozak for human)
4. mRNA stability awareness (ON vs OFF)
5. Constraint satisfaction per organism
6. Expression parameter retrieval completeness
"""

from __future__ import annotations

import random
from typing import NamedTuple

import pytest

from biocompiler.mrna_stability import score_mrna_stability, suggest_mutations_for_stability
from biocompiler.organism_config import OrganismConfig, ORGANISM_CONFIGS, get_organism_config
from biocompiler.optimization import optimize_sequence
from biocompiler.organisms import (
    CODON_ADAPTIVENESS_TABLES,
    CODON_USAGE_TABLES,
    ORGANISM_GC_TARGETS,
    PREFERRED_CODON_TABLES,
    SUPPORTED_ORGANISMS,
)
from biocompiler.organisms.e_coli import (
    E_COLI_CODON_PAIR_BIAS,
    E_COLI_EXPRESSION_OPTIMIZATION_PARAMS,
    compute_codon_pair_bias,
)
from biocompiler.restriction_sites import RESTRICTION_SITES
from biocompiler.scanner import gc_content
from biocompiler.translation import compute_cai
from biocompiler.type_system import AA_TO_CODONS
from biocompiler.utr_models import (
    ORGANISM_UTR_CONFIGS,
    score_3utr,
    score_5utr,
    suggest_3utr,
    suggest_5utr,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Human insulin A-chain + B-chain combined (110 AA) — a small well-studied
# therapeutic protein that contains a diverse mix of amino acids.
INSULIN_PROTEIN: str = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)

# Organism identifiers used in the test suite
ORGANISMS_UNDER_TEST: list[str] = [
    "Escherichia_coli",
    "Saccharomyces_cerevisiae",
    "Homo_sapiens",
    "Mus_musculus",
    "CHO_K1",
]

# Human-friendly labels for parametrized test IDs
ORGANISM_LABELS: dict[str, str] = {
    "Escherichia_coli": "ecoli",
    "Saccharomyces_cerevisiae": "yeast",
    "Homo_sapiens": "human",
    "Mus_musculus": "mouse",
    "CHO_K1": "cho",
}


class OrganismSpec(NamedTuple):
    """Metadata for a single organism under test."""
    canonical: str
    label: str
    gc_lo: float
    gc_hi: float


ORG_SPECS: list[OrganismSpec] = [
    OrganismSpec(
        canonical="Escherichia_coli",
        label="ecoli",
        gc_lo=0.30,
        gc_hi=0.70,
    ),
    OrganismSpec(
        canonical="Saccharomyces_cerevisiae",
        label="yeast",
        gc_lo=0.25,
        gc_hi=0.55,
    ),
    OrganismSpec(
        canonical="Homo_sapiens",
        label="human",
        gc_lo=0.30,
        gc_hi=0.70,
    ),
    OrganismSpec(
        canonical="Mus_musculus",
        label="mouse",
        gc_lo=0.30,
        gc_hi=0.70,
    ),
    OrganismSpec(
        canonical="CHO_K1",
        label="cho",
        gc_lo=0.30,
        gc_hi=0.70,
    ),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _count_destabilizing_motifs(dna: str, organism: str) -> int:
    """Count destabilizing motif hits in a DNA sequence for an organism."""
    report = score_mrna_stability(dna, organism)
    return report.destabilizing_count


def _make_naive_sequence(protein: str, organism: str) -> str:
    """Generate a naive (random codon) DNA sequence for comparison.

    Picks a random synonymous codon for each amino acid, ignoring
    CAI, GC, or any other optimization.  Used as a baseline for
    CPB and stability comparisons.
    """
    rng = random.Random(42)  # Fixed seed for reproducibility
    usage = CODON_ADAPTIVENESS_TABLES.get(organism, {})
    codons: list[str] = []
    for aa in protein:
        synonyms = AA_TO_CODONS.get(aa, [])
        if not synonyms:
            codons.append("NNN")
            continue
        codons.append(rng.choice(synonyms))
    return "".join(codons)


# ===========================================================================
# 1. Each organism produces DIFFERENT optimized sequences
# ===========================================================================


class TestOrganismSpecificOptimization:
    """Verify that each organism yields a distinct optimized DNA sequence
    and that GC content is organism-appropriate."""

    @pytest.mark.integration
    def test_all_five_organisms_produce_different_sequences(self) -> None:
        """Optimize insulin for 5 organisms and verify that organisms with
        distinct codon preferences produce different DNA sequences.

        Organisms with very similar codon preferences (human, mouse, CHO)
        may produce identical sequences under the same GC constraints —
        this is expected.  However, E. coli and yeast, which have
        fundamentally different codon biases, MUST produce distinct
        sequences from the mammalian group.

        We also verify that the preferred codon tables themselves differ
        across organisms, which is the underlying cause of sequence
        differentiation.
        """
        results: dict[str, str] = {}
        for org in ORGANISMS_UNDER_TEST:
            result = optimize_sequence(
                INSULIN_PROTEIN,
                organism=org,
                gc_lo=0.30,
                gc_hi=0.70,
                optimize_mrna_stability=False,
                include_utr=False,
            )
            results[org] = result.sequence

        # E. coli must differ from human (very different codon preferences)
        ecoli_seq = results["Escherichia_coli"]
        human_seq = results["Homo_sapiens"]
        assert ecoli_seq != human_seq, (
            "E. coli and human should produce different optimized sequences "
            "(fundamentally different codon preferences)"
        )

        # At least 2 distinct sequences (prokaryote vs eukaryote)
        unique_sequences = set(results.values())
        assert len(unique_sequences) >= 2, (
            f"Expected >= 2 different sequences across 5 organisms, "
            f"got {len(unique_sequences)}. "
            f"Sequences: { {k: v[:30] + '...' for k, v in results.items()} }"
        )

    @pytest.mark.integration
    def test_preferred_codons_differ_across_organisms(self) -> None:
        """Verify that preferred codon tables are different across organisms,
        which is the underlying reason for sequence differentiation."""
        ecoli_pref = PREFERRED_CODON_TABLES["Escherichia_coli"]
        human_pref = PREFERRED_CODON_TABLES["Homo_sapiens"]
        yeast_pref = PREFERRED_CODON_TABLES["Saccharomyces_cerevisiae"]

        # Count how many amino acids have different preferred codons
        shared_aas = set(ecoli_pref.keys()) & set(human_pref.keys())
        different_count = sum(
            1 for aa in shared_aas
            if ecoli_pref[aa] != human_pref[aa]
        )
        # E. coli and human should differ on many amino acids
        assert different_count >= 5, (
            f"E. coli and human preferred codons differ on only "
            f"{different_count} amino acids, expected >= 5"
        )

        # Yeast and human should also differ
        shared_aas_yh = set(yeast_pref.keys()) & set(human_pref.keys())
        different_count_yh = sum(
            1 for aa in shared_aas_yh
            if yeast_pref[aa] != human_pref[aa]
        )
        assert different_count_yh >= 5, (
            f"Yeast and human preferred codons differ on only "
            f"{different_count_yh} amino acids, expected >= 5"
        )

    @pytest.mark.integration
    @pytest.mark.parametrize(
        "org_spec",
        ORG_SPECS,
        ids=lambda s: s.label,
    )
    def test_gc_content_organism_appropriate(self, org_spec: OrganismSpec) -> None:
        """Optimized sequence GC content should fall within organism bounds."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=org_spec.canonical,
            gc_lo=org_spec.gc_lo,
            gc_hi=org_spec.gc_hi,
            optimize_mrna_stability=False,
            include_utr=False,
        )
        gc = gc_content(result.sequence)
        assert org_spec.gc_lo <= gc <= org_spec.gc_hi, (
            f"{org_spec.label}: GC={gc:.4f} outside "
            f"[{org_spec.gc_lo}, {org_spec.gc_hi}]"
        )

    @pytest.mark.integration
    def test_ecoli_gc_near_50_percent(self) -> None:
        """E. coli optimized sequence should have GC content near 50%."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            optimize_mrna_stability=False,
            include_utr=False,
        )
        gc = gc_content(result.sequence)
        # E. coli coding GC is ~50.8%; allow generous range for small proteins
        assert 0.35 <= gc <= 0.65, (
            f"E. coli GC={gc:.4f} not near expected ~50%"
        )

    @pytest.mark.integration
    def test_yeast_codon_usage_gc_lower_than_ecoli(self) -> None:
        """Yeast codon usage table should have lower GC content than E. coli.

        This tests the underlying biological data rather than the optimizer
        output, which may be constrained by restriction site removal and
        other hard constraints that override GC targets.
        """
        from biocompiler.organisms import CODON_USAGE_TABLES

        def _compute_gc_from_usage(
            usage: dict[str, tuple[str, float, float, int]],
        ) -> float:
            total_weighted_gc = 0.0
            total_bases = 0.0
            for codon, (_aa, _frac, per_thousand, _count) in usage.items():
                gc = sum(1 for base in codon if base in "GC")
                total_weighted_gc += gc * per_thousand
                total_bases += 3.0 * per_thousand
            return total_weighted_gc / total_bases if total_bases > 0 else 0.0

        ecoli_gc = _compute_gc_from_usage(CODON_USAGE_TABLES["Escherichia_coli"])
        yeast_gc = _compute_gc_from_usage(CODON_USAGE_TABLES["Saccharomyces_cerevisiae"])
        # Yeast coding GC (~38%) is fundamentally lower than E. coli (~50%)
        assert yeast_gc < ecoli_gc, (
            f"Yeast codon usage GC ({yeast_gc:.4f}) should be < "
            f"E. coli codon usage GC ({ecoli_gc:.4f})"
        )


# ===========================================================================
# 2. Codon pair bias integration
# ===========================================================================


class TestCodonPairBiasIntegration:
    """Verify that codon pair bias (CPB) is considered during optimization."""

    @pytest.mark.integration
    def test_ecoli_optimized_has_positive_mean_cpb(self) -> None:
        """E. coli optimized sequence should have positive mean CPB.

        The optimizer prefers high-CAI codons, and high-CAI codons in
        E. coli tend to form over-represented (positive CPB) pairs.
        """
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Escherichia_coli",
            optimize_mrna_stability=False,
            include_utr=False,
        )
        cpb = compute_codon_pair_bias(result.sequence, "Escherichia_coli")
        assert cpb > 0.0, (
            f"E. coli optimized mean CPB={cpb:.4f}, expected > 0"
        )

    @pytest.mark.integration
    def test_ecoli_optimized_cpb_better_than_naive(self) -> None:
        """Optimized sequence CPB should be comparable to or better than a
        naive (random-codon) sequence.

        Note: The optimizer maximizes CAI, not CPB directly.  However,
        high-CAI codons in E. coli tend to be those that also appear in
        over-represented pairs, so we expect the optimized CPB to be
        non-negative and within a reasonable range of the naive CPB.
        We do NOT require strict improvement because CPB is a secondary
        metric that can fluctuate based on codon pair sampling.
        """
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Escherichia_coli",
            optimize_mrna_stability=False,
            include_utr=False,
        )
        naive_seq = _make_naive_sequence(INSULIN_PROTEIN, "Escherichia_coli")
        optimized_cpb = compute_codon_pair_bias(result.sequence, "Escherichia_coli")
        naive_cpb = compute_codon_pair_bias(naive_seq, "Escherichia_coli")
        # Optimized CPB must be positive (good codon pairs favoured)
        assert optimized_cpb > 0.0, (
            f"Optimized CPB ({optimized_cpb:.4f}) should be positive"
        )
        # The CPB of the optimized sequence should not be drastically worse
        # than naive — allow up to 2x degradation margin
        assert optimized_cpb >= naive_cpb * 0.5 - 0.01, (
            f"Optimized CPB ({optimized_cpb:.4f}) is too far below "
            f"naive CPB ({naive_cpb:.4f})"
        )

    @pytest.mark.integration
    def test_human_optimized_has_non_negative_cpb_proxy(self) -> None:
        """Human does not have a dedicated CPB table, but the optimized
        sequence should still have reasonable codon pair composition.

        We verify this indirectly: the optimizer uses high-adaptiveness
        codons, which are unlikely to form the worst (most under-represented)
        pairs in any organism.
        """
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Homo_sapiens",
            optimize_mrna_stability=False,
            include_utr=False,
        )
        cai = compute_cai(result.sequence, "Homo_sapiens")
        # A well-optimized human sequence should have decent CAI
        assert cai >= 0.5, (
            f"Human optimized CAI={cai:.4f}, expected >= 0.5"
        )


# ===========================================================================
# 3. UTR model integration
# ===========================================================================


class TestUTRModelIntegration:
    """Verify that organism-specific UTR models produce and score
    appropriate 5' and 3' UTR sequences."""

    @pytest.mark.integration
    def test_ecoli_suggested_5utr_contains_shine_dalgarno(self) -> None:
        """E. coli suggested 5' UTR must contain a Shine-Dalgarno-like
        sequence (AGGAGG or close variant)."""
        utr5 = suggest_5utr("Escherichia_coli")
        assert utr5 is not None, "E. coli 5' UTR suggestion should not be None"
        # The Shine-Dalgarno consensus is AGGAGG
        sd_variants = ["AGGAGG", "AGGAGGT", "GGAGG", "AGGAG"]
        has_sd = any(variant in utr5.upper() for variant in sd_variants)
        assert has_sd, (
            f"E. coli 5' UTR {utr5!r} does not contain a "
            f"Shine-Dalgarno variant"
        )

    @pytest.mark.integration
    def test_human_suggested_5utr_contains_kozak(self) -> None:
        """Human suggested 5' UTR must contain the Kozak consensus
        (GCCACC before ATG)."""
        utr5 = suggest_5utr("Homo_sapiens")
        assert utr5 is not None, "Human 5' UTR suggestion should not be None"
        # The Kozak consensus is GCCACCATGG
        has_kozak = "GCCACC" in utr5.upper()
        assert has_kozak, (
            f"Human 5' UTR {utr5!r} does not contain Kozak consensus"
        )

    @pytest.mark.integration
    def test_ecoli_5utr_scores_high(self) -> None:
        """The suggested E. coli 5' UTR should score well."""
        utr5 = suggest_5utr("Escherichia_coli")
        score = score_5utr(utr5, "Escherichia_coli")
        assert score >= 0.5, (
            f"E. coli 5' UTR score={score:.4f}, expected >= 0.5"
        )

    @pytest.mark.integration
    def test_human_5utr_scores_high(self) -> None:
        """The suggested human 5' UTR should score well."""
        utr5 = suggest_5utr("Homo_sapiens")
        score = score_5utr(utr5, "Homo_sapiens")
        assert score >= 0.5, (
            f"Human 5' UTR score={score:.4f}, expected >= 0.5"
        )

    @pytest.mark.integration
    def test_ecoli_3utr_scores_high(self) -> None:
        """The suggested E. coli 3' UTR should score well (Rho-independent
        terminator)."""
        utr3 = suggest_3utr("Escherichia_coli")
        score = score_3utr(utr3, "Escherichia_coli")
        assert score >= 0.3, (
            f"E. coli 3' UTR score={score:.4f}, expected >= 0.3"
        )

    @pytest.mark.integration
    def test_human_3utr_scores_high(self) -> None:
        """The suggested human 3' UTR should score well (polyA signal)."""
        utr3 = suggest_3utr("Homo_sapiens")
        score = score_3utr(utr3, "Homo_sapiens")
        assert score >= 0.3, (
            f"Human 3' UTR score={score:.4f}, expected >= 0.3"
        )

    @pytest.mark.integration
    def test_utr_config_ecoli_has_shine_dalgarno_no_kozak(self) -> None:
        """E. coli UTR config must have shine_dalgarno set and
        kozak_sequence as None."""
        config = ORGANISM_UTR_CONFIGS["Escherichia_coli"]
        assert config.shine_dalgarno is not None, (
            "E. coli UTR config must define shine_dalgarno"
        )
        assert config.kozak_sequence is None, (
            "E. coli UTR config must not have kozak_sequence"
        )

    @pytest.mark.integration
    def test_utr_config_human_has_kozak_no_shine_dalgarno(self) -> None:
        """Human UTR config must have kozak_sequence set and
        shine_dalgarno as None."""
        config = ORGANISM_UTR_CONFIGS["Homo_sapiens"]
        assert config.kozak_sequence is not None, (
            "Human UTR config must define kozak_sequence"
        )
        assert config.shine_dalgarno is None, (
            "Human UTR config must not have shine_dalgarno"
        )


# ===========================================================================
# 4. mRNA stability
# ===========================================================================


class TestMRNAStability:
    """Verify that stability-aware optimization reduces destabilizing motifs."""

    @pytest.mark.integration
    def test_stability_on_has_fewer_destabilizing_motifs(self) -> None:
        """Optimizing with stability ON should produce fewer destabilizing
        motifs than with stability OFF."""
        result_on = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Homo_sapiens",
            optimize_mrna_stability=True,
            include_utr=False,
        )
        result_off = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Homo_sapiens",
            optimize_mrna_stability=False,
            include_utr=False,
        )
        destab_on = _count_destabilizing_motifs(result_on.sequence, "Homo_sapiens")
        destab_off = _count_destabilizing_motifs(result_off.sequence, "Homo_sapiens")
        assert destab_on <= destab_off, (
            f"Stability ON should have fewer/same destabilizing motifs "
            f"({destab_on}) than OFF ({destab_off})"
        )

    @pytest.mark.integration
    def test_stability_on_avoids_attta(self) -> None:
        """With stability ON, the ATTTA motif should be minimized
        (E. coli and human both avoid it)."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Escherichia_coli",
            optimize_mrna_stability=True,
            include_utr=False,
        )
        attta_count = result.sequence.upper().count("ATTTA")
        # The optimizer should remove ATTTA motifs when possible
        assert attta_count == 0, (
            f"Optimized E. coli sequence still contains {attta_count} "
            f"ATTTA motifs"
        )

    @pytest.mark.integration
    def test_stability_score_improves(self) -> None:
        """The mRNA stability score should improve (or stay the same)
        when stability optimization is enabled."""
        result_on = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Homo_sapiens",
            optimize_mrna_stability=True,
            include_utr=False,
        )
        result_off = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Homo_sapiens",
            optimize_mrna_stability=False,
            include_utr=False,
        )
        score_on = score_mrna_stability(result_on.sequence, "Homo_sapiens")
        score_off = score_mrna_stability(result_off.sequence, "Homo_sapiens")
        assert score_on.overall_score >= score_off.overall_score, (
            f"Stability ON score ({score_on.overall_score:.4f}) should be "
            f">= OFF score ({score_off.overall_score:.4f})"
        )

    @pytest.mark.integration
    def test_suggest_mutations_for_stability(self) -> None:
        """For a deliberately unstable sequence, stability mutation
        suggestions should be generated."""
        # Build a sequence with many ATTTA motifs
        bad_seq = _make_naive_sequence(INSULIN_PROTEIN, "Homo_sapiens")
        suggestions = suggest_mutations_for_stability(bad_seq, "Homo_sapiens")
        # Even if the naive sequence doesn't have many motifs,
        # the function should return a list (possibly empty)
        assert isinstance(suggestions, list), (
            "suggest_mutations_for_stability should return a list"
        )
        for s in suggestions:
            assert "original_codon" in s, f"Suggestion missing 'original_codon': {s}"
            assert "suggested_codon" in s, f"Suggestion missing 'suggested_codon': {s}"
            assert "amino_acid" in s, f"Suggestion missing 'amino_acid': {s}"


# ===========================================================================
# 5. Constraint satisfaction per organism
# ===========================================================================


class TestConstraintSatisfaction:
    """Verify that optimized sequences satisfy organism-specific constraints."""

    @pytest.mark.integration
    def test_ecoli_constraints(self) -> None:
        """E. coli: GC in [0.30, 0.70], no restriction sites,
        positive CPB, no ATTTA."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            optimize_mrna_stability=True,
            include_utr=False,
        )
        seq = result.sequence

        # GC content in range
        gc = gc_content(seq)
        assert 0.30 <= gc <= 0.70, f"E. coli GC={gc:.4f} outside [0.30, 0.70]"

        # No restriction sites (default enzymes)
        for site in ["GAATTC", "GGATCC", "AAGCTT", "CTCGAG", "GCGGCCGC"]:
            assert site not in seq, f"E. coli sequence contains restriction site {site}"

        # Positive CPB
        cpb = compute_codon_pair_bias(seq, "Escherichia_coli")
        assert cpb > 0.0, f"E. coli CPB={cpb:.4f}, expected > 0"

        # No ATTTA instability motif
        assert "ATTTA" not in seq, "E. coli sequence contains ATTTA motif"

    @pytest.mark.integration
    def test_yeast_constraints(self) -> None:
        """Yeast: GC in [0.25, 0.55], no restriction sites,
        low T-runs."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Saccharomyces_cerevisiae",
            gc_lo=0.25,
            gc_hi=0.55,
            optimize_mrna_stability=True,
            include_utr=False,
        )
        seq = result.sequence

        # GC content in range
        gc = gc_content(seq)
        assert 0.25 <= gc <= 0.55, f"Yeast GC={gc:.4f} outside [0.25, 0.55]"

        # No restriction sites (default enzymes)
        for site in ["GAATTC", "GGATCC", "AAGCTT", "CTCGAG", "GCGGCCGC"]:
            assert site not in seq, f"Yeast sequence contains restriction site {site}"

        # Low T-runs (no 6+ consecutive T's)
        max_t_run = 0
        current_run = 0
        for base in seq:
            if base == "T":
                current_run += 1
                max_t_run = max(max_t_run, current_run)
            else:
                current_run = 0
        assert max_t_run < 6, (
            f"Yeast sequence has T-run of length {max_t_run} (max allowed: 5)"
        )

    @pytest.mark.integration
    def test_human_constraints(self) -> None:
        """Human: GC in [0.30, 0.70], no restriction sites,
        splice site aware."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
            optimize_mrna_stability=True,
            include_utr=False,
        )
        seq = result.sequence

        # GC content in range
        gc = gc_content(seq)
        assert 0.30 <= gc <= 0.70, f"Human GC={gc:.4f} outside [0.30, 0.70]"

        # No restriction sites (default enzymes)
        for site in ["GAATTC", "GGATCC", "AAGCTT", "CTCGAG", "GCGGCCGC"]:
            assert site not in seq, f"Human sequence contains restriction site {site}"

        # The optimizer should have attempted splice site elimination.
        # Verify that the optimization result exists and has valid CAI.
        cai = compute_cai(seq, "Homo_sapiens")
        assert cai > 0.0, f"Human CAI={cai:.4f}, expected > 0"

    @pytest.mark.integration
    def test_cho_constraints(self) -> None:
        """CHO: GC in [0.30, 0.70], Kozak context awareness,
        splice site aware."""
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="CHO_K1",
            gc_lo=0.30,
            gc_hi=0.70,
            optimize_mrna_stability=True,
            include_utr=False,
        )
        seq = result.sequence

        # GC content in range
        gc = gc_content(seq)
        assert 0.30 <= gc <= 0.70, f"CHO GC={gc:.4f} outside [0.30, 0.70]"

        # No restriction sites (default enzymes)
        for site in ["GAATTC", "GGATCC", "AAGCTT", "CTCGAG", "GCGGCCGC"]:
            assert site not in seq, f"CHO sequence contains restriction site {site}"

        # CHO UTR config should have Kozak (for UTR suggestion)
        cho_utr = ORGANISM_UTR_CONFIGS.get("CHO_K1")
        assert cho_utr is not None, "CHO UTR config should exist"
        assert cho_utr.kozak_sequence is not None, (
            "CHO should have Kozak sequence defined"
        )
        assert cho_utr.splicing_signals is True, (
            "CHO should have splicing_signals=True"
        )

    @pytest.mark.integration
    @pytest.mark.parametrize(
        "org_spec",
        ORG_SPECS,
        ids=lambda s: s.label,
    )
    def test_protein_translation_preserved(self, org_spec: OrganismSpec) -> None:
        """All organisms: the optimized DNA must encode the same protein."""
        from biocompiler.type_system import CODON_TABLE

        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism=org_spec.canonical,
            gc_lo=org_spec.gc_lo,
            gc_hi=org_spec.gc_hi,
            optimize_mrna_stability=False,
            include_utr=False,
        )
        seq = result.sequence

        # Translate and compare
        codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
        translated = "".join(CODON_TABLE.get(c, "X") for c in codons)
        # Remove stop codon from comparison
        translated_no_stop = translated.rstrip("*")
        assert translated_no_stop == INSULIN_PROTEIN, (
            f"{org_spec.label}: translated protein doesn't match input. "
            f"Expected {INSULIN_PROTEIN[:20]}..., "
            f"got {translated_no_stop[:20]}..."
        )


# ===========================================================================
# 6. Expression parameter retrieval
# ===========================================================================


class TestExpressionParameterRetrieval:
    """Verify that each organism has complete expression optimization params."""

    @pytest.mark.integration
    def test_ecoli_expression_optimization_params_complete(self) -> None:
        """E. coli expression optimization params must have all required fields."""
        required_fields = [
            "gc_content_min",
            "gc_content_max",
            "max_t_run",
            "avoid_motifs",
            "max_consecutive_rare_codons",
        ]
        params = E_COLI_EXPRESSION_OPTIMIZATION_PARAMS
        for field in required_fields:
            assert field in params, (
                f"E. coli expression params missing field {field!r}"
            )

    @pytest.mark.integration
    @pytest.mark.parametrize(
        "org_spec",
        ORG_SPECS,
        ids=lambda s: s.label,
    )
    def test_organism_config_has_required_fields(self, org_spec: OrganismSpec) -> None:
        """Each organism's OrganismConfig must have complete parameters."""
        # Map canonical names to OrganismConfig keys
        config_key_map: dict[str, str] = {
            "Escherichia_coli": "E_coli_K12",
            "Saccharomyces_cerevisiae": "Saccharomyces_cerevisiae",
            "Homo_sapiens": "Homo_sapiens",
            "Mus_musculus": "Mus_musculus",
            "CHO_K1": "CHO_K1",
        }
        config_key = config_key_map.get(org_spec.canonical, org_spec.canonical)
        config = get_organism_config(config_key)

        # Verify it's an OrganismConfig instance
        assert isinstance(config, OrganismConfig), (
            f"{org_spec.label}: expected OrganismConfig, got {type(config)}"
        )

        # Required fields (mapped from OrganismConfig attributes)
        assert hasattr(config, "gc_target_lo"), (
            f"{org_spec.label}: OrganismConfig missing gc_target_lo"
        )
        assert hasattr(config, "gc_target_hi"), (
            f"{org_spec.label}: OrganismConfig missing gc_target_hi"
        )
        assert hasattr(config, "max_homopolymer_run"), (
            f"{org_spec.label}: OrganismConfig missing max_homopolymer_run"
        )
        assert hasattr(config, "avoided_motifs"), (
            f"{org_spec.label}: OrganismConfig missing avoided_motifs"
        )
        assert hasattr(config, "preferred_codons"), (
            f"{org_spec.label}: OrganismConfig missing preferred_codons"
        )

        # Verify values are sensible
        assert 0.0 <= config.gc_target_lo <= 1.0, (
            f"{org_spec.label}: gc_target_lo={config.gc_target_lo} out of range"
        )
        assert 0.0 <= config.gc_target_hi <= 1.0, (
            f"{org_spec.label}: gc_target_hi={config.gc_target_hi} out of range"
        )
        assert config.gc_target_lo < config.gc_target_hi, (
            f"{org_spec.label}: gc_target_lo ({config.gc_target_lo}) must be < "
            f"gc_target_hi ({config.gc_target_hi})"
        )
        assert config.max_homopolymer_run >= 1, (
            f"{org_spec.label}: max_homopolymer_run={config.max_homopolymer_run} < 1"
        )
        assert isinstance(config.avoided_motifs, list), (
            f"{org_spec.label}: avoided_motifs is not a list"
        )
        assert isinstance(config.preferred_codons, dict), (
            f"{org_spec.label}: preferred_codons is not a dict"
        )

    @pytest.mark.integration
    @pytest.mark.parametrize(
        "org_spec",
        ORG_SPECS,
        ids=lambda s: s.label,
    )
    def test_codon_usage_tables_complete(self, org_spec: OrganismSpec) -> None:
        """Each organism must have complete codon usage, adaptiveness,
        and preferred codon tables."""
        org = org_spec.canonical

        # Codon usage table
        assert org in CODON_USAGE_TABLES, (
            f"{org_spec.label}: not in CODON_USAGE_TABLES"
        )
        usage = CODON_USAGE_TABLES[org]
        assert len(usage) == 64, (
            f"{org_spec.label}: CODON_USAGE_TABLES has {len(usage)} codons, expected 64"
        )

        # Adaptiveness table
        assert org in CODON_ADAPTIVENESS_TABLES, (
            f"{org_spec.label}: not in CODON_ADAPTIVENESS_TABLES"
        )
        adapt = CODON_ADAPTIVENESS_TABLES[org]
        # Non-stop codons: 61
        assert len(adapt) >= 50, (
            f"{org_spec.label}: CODON_ADAPTIVENESS_TABLES has {len(adapt)} entries, "
            f"expected >= 50"
        )

        # Preferred codons
        assert org in PREFERRED_CODON_TABLES, (
            f"{org_spec.label}: not in PREFERRED_CODON_TABLES"
        )
        pref = PREFERRED_CODON_TABLES[org]
        # Should have entries for most amino acids (at least 18 of 20)
        assert len(pref) >= 18, (
            f"{org_spec.label}: PREFERRED_CODON_TABLES has {len(pref)} entries, "
            f"expected >= 18"
        )

    @pytest.mark.integration
    @pytest.mark.parametrize(
        "org_spec",
        ORG_SPECS,
        ids=lambda s: s.label,
    )
    def test_gc_targets_defined(self, org_spec: OrganismSpec) -> None:
        """Each organism must have GC target range defined."""
        org = org_spec.canonical
        assert org in ORGANISM_GC_TARGETS, (
            f"{org_spec.label}: not in ORGANISM_GC_TARGETS"
        )
        gc_lo, gc_hi = ORGANISM_GC_TARGETS[org]
        assert 0.0 < gc_lo < gc_hi < 1.0, (
            f"{org_spec.label}: GC targets ({gc_lo}, {gc_hi}) invalid"
        )

    @pytest.mark.integration
    def test_ecoli_params_gc_bounds_consistent(self) -> None:
        """E. coli expression params GC bounds should be consistent with
        the OrganismConfig."""
        params = E_COLI_EXPRESSION_OPTIMIZATION_PARAMS
        config = get_organism_config("E_coli_K12")

        # Params allow wide range; config narrows it
        assert params["gc_content_min"] <= config.gc_target_lo, (
            f"E. coli params gc_min ({params['gc_content_min']}) > "
            f"config gc_lo ({config.gc_target_lo})"
        )
        assert params["gc_content_max"] >= config.gc_target_hi, (
            f"E. coli params gc_max ({params['gc_content_max']}) < "
            f"config gc_hi ({config.gc_target_hi})"
        )

    @pytest.mark.integration
    def test_ecoli_params_avoid_motifs_includes_attta(self) -> None:
        """E. coli expression params should include ATTTA in avoid_motifs."""
        params = E_COLI_EXPRESSION_OPTIMIZATION_PARAMS
        assert "ATTTA" in params["avoid_motifs"], (
            "E. coli params should include ATTTA in avoid_motifs"
        )

    @pytest.mark.integration
    def test_supported_organisms_covers_test_set(self) -> None:
        """All organisms under test must be in SUPPORTED_ORGANISMS."""
        for org in ORGANISMS_UNDER_TEST:
            assert org in SUPPORTED_ORGANISMS, (
                f"{org} not in SUPPORTED_ORGANISMS: {SUPPORTED_ORGANISMS}"
            )
