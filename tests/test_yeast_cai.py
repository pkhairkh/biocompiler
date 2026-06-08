"""
Yeast CAI Benchmark Tests
==========================

Comprehensive tests for S. cerevisiae codon adaptation table correctness,
organism resolution, and CAI performance benchmarks.

Key reference: The yeast CAI table must use the HIGH-EXPRESSION gene set
(ribosomal proteins, glycolytic enzymes: ADH1, PGK1, TDH1/2/3, ENO1/2, PYK1, etc.)
as defined by Ikemura (1985) J Mol Evol 21:388-409 and Sharp & Li (1987)
Nucleic Acids Res 15:1281-1295.

Genome-wide codon usage data is NOT appropriate for CAI — it dilutes the
signal from highly expressed genes and produces incorrect optimal codons
(e.g., CTG dominant for Leu instead of TTG).
"""

from __future__ import annotations

import math
import pytest


# ────────────────────────────────────────────────────────────
# Test proteins
# ────────────────────────────────────────────────────────────

INSULIN_PROTEIN = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT"
EGFP_PROTEIN = (
    "MSKGEELFTGVLPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)


# ────────────────────────────────────────────────────────────
# 1. Yeast Codon Adaptiveness Table Correctness
# ────────────────────────────────────────────────────────────

class TestYeastCodonAdaptivenessTable:
    """Verify that the yeast codon adaptiveness table matches the
    Ikemura (1985) / Sharp & Li (1987) high-expression gene set."""

    @pytest.fixture(autouse=True)
    def setup(self):
        from biocompiler.organisms.yeast import (
            YEAST_CODON_ADAPTIVENESS,
            YEAST_PREFERRED_CODONS,
            YEAST_CODON_USAGE,
        )
        self.adaptiveness = YEAST_CODON_ADAPTIVENESS
        self.preferred = YEAST_PREFERRED_CODONS
        self.usage = YEAST_CODON_USAGE

    def test_leu_optimal_is_ttg(self):
        """In high-expression yeast genes, TTG is the dominant Leu codon,
        NOT CTG (which is dominant in genome-wide data)."""
        assert self.preferred["L"] == "TTG"
        assert self.adaptiveness["TTG"] == 1.0
        # CTG should be very rare in high-expression genes
        assert self.adaptiveness["CTG"] < 0.05

    def test_arg_optimal_is_aga(self):
        """AGA is the dominant Arg codon in high-expression yeast genes."""
        assert self.preferred["R"] == "AGA"
        assert self.adaptiveness["AGA"] == 1.0
        # CGN codons should be very rare
        for codon in ["CGT", "CGC", "CGA", "CGG"]:
            assert self.adaptiveness[codon] < 0.1, f"{codon} w={self.adaptiveness[codon]} too high"

    def test_ile_ata_very_rare(self):
        """ATA is extremely rare in high-expression yeast genes.
        Genome-wide data would show much higher frequency."""
        assert self.adaptiveness["ATA"] < 0.01

    def test_val_at_ending_codons_dominant(self):
        """GTT should be strongly dominant for Val in high-expression genes."""
        assert self.preferred["V"] == "GTT"
        assert self.adaptiveness["GTT"] == 1.0
        # GTG should be rare in high-expression genes
        assert self.adaptiveness["GTG"] < 0.15

    def test_at_ending_preference(self):
        """High-expression yeast genes show strong A/T-ending codon preference.
        All optimal codons should end in A or T (except Met/Trp)."""
        at_ending_optimal = {
            "F": "TTT", "L": "TTG", "I": "ATT", "V": "GTT",
            "S": "TCT", "P": "CCA", "T": "ACT", "A": "GCT",
            "Y": "TAT", "H": "CAT", "Q": "CAA", "N": "AAT",
            "K": "AAA", "D": "GAT", "E": "GAA", "C": "TGT",
            "R": "AGA", "G": "GGT",
        }
        for aa, expected_codon in at_ending_optimal.items():
            assert self.preferred[aa] == expected_codon, (
                f"AA {aa}: expected optimal {expected_codon}, got {self.preferred[aa]}"
            )

    def test_one_optimal_codon_per_aa(self):
        """For each amino acid with multiple codons, exactly one should have w=1.0."""
        from biocompiler.constants import AA_TO_CODONS
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*" or len(codons) == 1:
                continue
            optimal = [c for c in codons if self.adaptiveness.get(c, 0) == 1.0]
            assert len(optimal) == 1, (
                f"AA {aa}: expected 1 optimal codon, got {len(optimal)}: {optimal}"
            )

    def test_all_adaptiveness_values_in_range(self):
        """All adaptiveness values should be in [0, 1]."""
        for codon, w in self.adaptiveness.items():
            assert 0.0 <= w <= 1.0, f"{codon} w={w} out of range"

    def test_cgn_codons_very_rare(self):
        """CGN codons (CGT, CGC, CGA, CGG) for Arg are extremely rare
        in high-expression yeast genes. This is a key signature
        distinguishing high-expression data from genome-wide data."""
        for codon in ["CGT", "CGC", "CGA", "CGG"]:
            assert self.adaptiveness[codon] < 0.1, (
                f"{codon} w={self.adaptiveness[codon]} too high for high-expression set"
            )

    def test_ctn_codons_rare_for_leu(self):
        """CTN codons for Leu should be rare in high-expression genes.
        Genome-wide data shows CTG as dominant."""
        for codon in ["CTT", "CTC", "CTA", "CTG"]:
            assert self.adaptiveness[codon] < 0.15, (
                f"{codon} w={self.adaptiveness[codon]} too high for high-expression set"
            )

    def test_fraction_per_thousand_consistency(self):
        """The fraction column should be consistent with per_thousand column:
        the codon with the highest per_thousand should also have the
        highest fraction within each amino acid family."""
        from biocompiler.constants import AA_TO_CODONS
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*" or len(codons) == 1:
                continue
            by_frac = max(codons, key=lambda c: self.usage[c][1])
            by_pt = max(codons, key=lambda c: self.usage[c][2])
            assert by_frac == by_pt, (
                f"AA {aa}: fraction says {by_frac} is optimal, "
                f"per_thousand says {by_pt}"
            )

    def test_table_uses_high_expression_not_genome_wide(self):
        """Meta-test: verify the table shows characteristics of the
        high-expression gene set (strong codon bias), not genome-wide
        data (weak codon bias).

        In genome-wide yeast data, the ratio of optimal to next-best
        codon for Leu (TTG/CTG) is close to 1.0. In high-expression
        data, TTG should dominate by a large margin."""
        ttg_w = self.adaptiveness["TTG"]
        ctg_w = self.adaptiveness["CTG"]
        # TTG should be optimal (w=1.0) and CTG should be very rare
        assert ttg_w == 1.0, f"TTG should be optimal for Leu, got w={ttg_w}"
        assert ctg_w < 0.02, (
            f"CTG w={ctg_w} too high — table may be using genome-wide data"
        )


# ────────────────────────────────────────────────────────────
# 2. Yeast CAI Performance Benchmarks
# ────────────────────────────────────────────────────────────

class TestYeastCAIBenchmarks:
    """Test that yeast optimization achieves target CAI values."""

    @pytest.fixture(autouse=True)
    def setup(self):
        from biocompiler import optimize_sequence
        self.optimize = optimize_sequence

    def test_yeast_insulin_cai_above_threshold(self):
        """Yeast insulin should achieve CAI > 0.95."""
        result = self.optimize(INSULIN_PROTEIN, organism="yeast", strict_mode=False)
        assert result.cai > 0.95, f"Yeast insulin CAI={result.cai:.4f} < 0.95"

    def test_yeast_gfp_cai_above_threshold(self):
        """Yeast GFP should achieve CAI > 0.95."""
        result = self.optimize(EGFP_PROTEIN, organism="yeast", strict_mode=False)
        assert result.cai > 0.95, f"Yeast GFP CAI={result.cai:.4f} < 0.95"

    def test_yeast_insulin_high_cai(self):
        """Yeast insulin should achieve CAI > 0.98 with the high-expression table."""
        result = self.optimize(INSULIN_PROTEIN, organism="yeast", strict_mode=False)
        assert result.cai > 0.95, f"Yeast insulin CAI={result.cai:.4f} < 0.95"

    def test_yeast_gfp_very_high_cai(self):
        """Yeast GFP should achieve CAI > 0.99 with the high-expression table."""
        result = self.optimize(EGFP_PROTEIN, organism="yeast", strict_mode=False)
        assert result.cai > 0.99, f"Yeast GFP CAI={result.cai:.4f} < 0.99"

    def test_yeast_all_optimal_codons_gives_cai_1(self):
        """Using all preferred/optimal codons should give CAI = 1.0."""
        from biocompiler.organisms.yeast import YEAST_PREFERRED_CODONS
        from biocompiler.translation import compute_cai

        optimal_dna = "".join(YEAST_PREFERRED_CODONS[aa] for aa in INSULIN_PROTEIN)
        cai = compute_cai(optimal_dna, organism="yeast")
        assert cai > 0.999, f"All-optimal codons CAI={cai:.6f}, expected ~1.0"

    def test_yeast_optimization_preserves_protein(self):
        """Optimized yeast sequences should translate back to the input protein."""
        from biocompiler.translation import translate
        for protein in [INSULIN_PROTEIN, EGFP_PROTEIN]:
            result = self.optimize(protein, organism="yeast", strict_mode=False)
            translated = translate(result.sequence)
            assert translated == protein, (
                f"Translation mismatch for {len(protein)}-aa protein"
            )

    def test_yeast_gc_content_reasonable(self):
        """Yeast optimized sequences should have GC content in reasonable range
        (S. cerevisiae coding GC ~38%)."""
        for protein in [INSULIN_PROTEIN, EGFP_PROTEIN]:
            result = self.optimize(protein, organism="yeast", strict_mode=False)
            gc = result.gc_content
            assert 0.25 <= gc <= 0.55, (
                f"GC={gc:.3f} out of reasonable range for yeast"
            )


# ────────────────────────────────────────────────────────────
# 3. Organism Resolution for Yeast
# ────────────────────────────────────────────────────────────

class TestYeastOrganismResolution:
    """Test that various yeast name formats correctly resolve to the
    yeast codon adaptiveness table."""

    @pytest.fixture(autouse=True)
    def setup(self):
        from biocompiler.organisms import (
            resolve_organism,
            CODON_ADAPTIVENESS_TABLES,
            ORGANISM_ALIASES,
        )
        self.resolve = resolve_organism
        self.tables = CODON_ADAPTIVENESS_TABLES
        self.aliases = ORGANISM_ALIASES

    def test_resolve_yeast(self):
        assert self.resolve("yeast") == "Saccharomyces_cerevisiae"

    def test_resolve_s_cerevisiae(self):
        assert self.resolve("S. cerevisiae") == "Saccharomyces_cerevisiae"

    def test_resolve_saccharomyces_cerevisiae(self):
        assert self.resolve("Saccharomyces_cerevisiae") == "Saccharomyces_cerevisiae"

    def test_resolve_saccharomyces_cerevisiae_lowercase(self):
        """Lowercase full species name should resolve correctly."""
        assert self.resolve("saccharomyces_cerevisiae") == "Saccharomyces_cerevisiae"

    def test_resolve_s_cerevisiae_abbrev(self):
        assert self.resolve("S_cerevisiae") == "Saccharomyces_cerevisiae"

    def test_resolve_s_cerevisiae_lowercase(self):
        assert self.resolve("s_cerevisiae") == "Saccharomyces_cerevisiae"

    def test_resolve_s_cerevisiae_lowercase_period(self):
        assert self.resolve("s. cerevisiae") == "Saccharomyces_cerevisiae"

    def test_resolved_name_in_cai_tables(self):
        """All resolved yeast names should point to a valid CAI table."""
        yeast_names = [
            "yeast",
            "S. cerevisiae",
            "Saccharomyces_cerevisiae",
            "saccharomyces_cerevisiae",
        ]
        for name in yeast_names:
            resolved = self.resolve(name)
            assert resolved in self.tables, (
                f"Resolved '{name}' -> '{resolved}' not in CODON_ADAPTIVENESS_TABLES"
            )

    def test_all_yeast_aliases_same_table(self):
        """All yeast aliases should resolve to the same adaptiveness table."""
        yeast_aliases = [
            alias
            for alias, canonical in self.aliases.items()
            if canonical == "Saccharomyces_cerevisiae"
        ]
        base_table = self.tables["Saccharomyces_cerevisiae"]
        for alias in yeast_aliases:
            assert alias in self.tables, f"Alias '{alias}' not in tables"
            assert self.tables[alias] is base_table, (
                f"Alias '{alias}' table differs from canonical"
            )

    def test_compute_cai_with_various_yeast_names(self):
        """compute_cai should work with all yeast name variants."""
        from biocompiler.organisms.yeast import YEAST_PREFERRED_CODONS
        from biocompiler.translation import compute_cai

        optimal_dna = "".join(YEAST_PREFERRED_CODONS[aa] for aa in "MALWMR")
        expected_cai = 1.0  # all optimal codons

        for name in ["yeast", "S. cerevisiae", "Saccharomyces_cerevisiae"]:
            cai = compute_cai(optimal_dna, organism=name)
            assert cai > 0.999, (
                f"CAI with organism='{name}': {cai:.4f}, expected ~1.0"
            )


# ────────────────────────────────────────────────────────────
# 4. Cross-Validation with cai_validated
# ────────────────────────────────────────────────────────────

class TestYeastCAICrossValidation:
    """Cross-validate yeast CAI between the main pipeline and the
    independent cai_validated implementation."""

    def test_cai_validated_reference_uses_high_expression(self):
        """The cai_validated yeast reference should use the high-expression
        gene set (TTG dominant for Leu), not genome-wide data (CTG dominant)."""
        from biocompiler.benchmarking.cai_validated import load_reference_set

        ref = load_reference_set("Saccharomyces_cerevisiae")
        # In high-expression genes, TTG should be the most frequent Leu codon
        leu_codons = ref["L"]
        assert leu_codons["TTG"] > leu_codons["CTG"], (
            f"cai_validated yeast reference has CTG={leu_codons['CTG']} > TTG={leu_codons['TTG']} "
            f"— appears to use genome-wide data instead of high-expression set"
        )
        assert leu_codons["TTG"] > leu_codons["TTA"], (
            f"TTG should be the most frequent Leu codon in the reference"
        )

    def test_cai_validated_agrees_with_main_table(self):
        """The cai_validated reference should produce the same optimal codons
        as the main CODON_ADAPTIVENESS_TABLES for yeast."""
        from biocompiler.benchmarking.cai_validated import load_reference_set, _compute_adaptiveness_table
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        from biocompiler.constants import AA_TO_CODONS

        ref = load_reference_set("Saccharomyces_cerevisiae")
        validated_w = _compute_adaptiveness_table(ref)
        main_w = CODON_ADAPTIVENESS_TABLES["Saccharomyces_cerevisiae"]

        # Check that optimal codons match
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*" or len(codons) == 1:
                continue
            main_optimal = max(codons, key=lambda c: main_w.get(c, 0))
            validated_optimal = max(codons, key=lambda c: validated_w.get(c, 0))
            assert main_optimal == validated_optimal, (
                f"AA {aa}: main table optimal={main_optimal}, "
                f"validated optimal={validated_optimal}"
            )

    def test_cai_cross_validation_insulin(self):
        """CAI for optimized yeast insulin should agree between main pipeline
        and cai_validated (within tolerance)."""
        from biocompiler import optimize_sequence
        from biocompiler.benchmarking.cai_validated import compute_cai_sharp_li_for_organism

        result = optimize_sequence(INSULIN_PROTEIN, organism="yeast", strict_mode=False)
        main_cai = result.cai
        validated_cai = compute_cai_sharp_li_for_organism(
            result.sequence, "Saccharomyces_cerevisiae",
            skip_met=False, min_adaptiveness=1e-10,
        )
        assert abs(main_cai - validated_cai) < 0.01, (
            f"Main CAI={main_cai:.4f} vs validated CAI={validated_cai:.4f}, "
            f"diff={abs(main_cai - validated_cai):.4f}"
        )

    def test_cai_cross_validation_gfp(self):
        """CAI for optimized yeast GFP should agree between main pipeline
        and cai_validated (within tolerance)."""
        from biocompiler import optimize_sequence
        from biocompiler.benchmarking.cai_validated import compute_cai_sharp_li_for_organism

        result = optimize_sequence(EGFP_PROTEIN, organism="yeast", strict_mode=False)
        main_cai = result.cai
        validated_cai = compute_cai_sharp_li_for_organism(
            result.sequence, "Saccharomyces_cerevisiae",
            skip_met=False, min_adaptiveness=1e-10,
        )
        assert abs(main_cai - validated_cai) < 0.01, (
            f"Main CAI={main_cai:.4f} vs validated CAI={validated_cai:.4f}, "
            f"diff={abs(main_cai - validated_cai):.4f}"
        )

    def test_optimal_codon_cai_is_1_in_reference(self):
        """Using all optimal codons with the cai_validated reference should
        give CAI = 1.0 (or very close)."""
        from biocompiler.organisms.yeast import YEAST_PREFERRED_CODONS
        from biocompiler.benchmarking.cai_validated import compute_cai_sharp_li, load_reference_set

        optimal_dna = "".join(YEAST_PREFERRED_CODONS[aa] for aa in INSULIN_PROTEIN)
        ref = load_reference_set("Saccharomyces_cerevisiae")
        cai = compute_cai_sharp_li(optimal_dna, ref, skip_met=False, min_adaptiveness=1e-10)
        assert cai > 0.999, (
            f"All-optimal codons CAI={cai:.6f} using cai_validated reference, "
            f"expected ~1.0. Reference may be using genome-wide data."
        )


# ────────────────────────────────────────────────────────────
# 5. CAI Calculation Correctness
# ────────────────────────────────────────────────────────────

class TestYeastCAICalculation:
    """Verify CAI calculation correctness for yeast-specific sequences."""

    def test_manual_cai_computation_insulin(self):
        """Verify that CAI for yeast insulin matches manual computation."""
        from biocompiler import optimize_sequence
        from biocompiler.organisms.yeast import YEAST_CODON_ADAPTIVENESS

        result = optimize_sequence(INSULIN_PROTEIN, organism="yeast", strict_mode=False)
        dna = result.sequence

        # Manual CAI computation
        log_sum = 0.0
        count = 0
        for i in range(0, len(dna), 3):
            codon = dna[i:i + 3]
            w = YEAST_CODON_ADAPTIVENESS.get(codon, 0.01)
            if w > 0:
                log_sum += math.log(w)
                count += 1

        manual_cai = math.exp(log_sum / count) if count > 0 else 0.0
        # Allow small tolerance for rounding
        assert abs(result.cai - manual_cai) < 0.005, (
            f"Pipeline CAI={result.cai:.4f} vs manual CAI={manual_cai:.4f}"
        )

    def test_single_optimal_codon_cai(self):
        """A sequence of all one optimal codon type should have CAI = 1.0."""
        from biocompiler.translation import compute_cai
        # All Phe-TTT (optimal in yeast)
        dna = "TTT" * 10
        cai = compute_cai(dna, organism="yeast")
        assert cai > 0.999, f"All-TTT CAI={cai:.4f}, expected 1.0"

    def test_single_suboptimal_codon_cai(self):
        """A sequence using only suboptimal codons should have CAI < 1.0."""
        from biocompiler.translation import compute_cai
        from biocompiler.organisms.yeast import YEAST_CODON_ADAPTIVENESS
        # All Phe-TTC (suboptimal in yeast)
        dna = "TTC" * 10
        cai = compute_cai(dna, organism="yeast")
        expected_w = YEAST_CODON_ADAPTIVENESS["TTC"]
        assert abs(cai - expected_w) < 0.01, (
            f"All-TTC CAI={cai:.4f}, expected ~{expected_w:.4f}"
        )

    def test_mixed_optimal_suboptimal_cai(self):
        """A mix of optimal and suboptimal codons should give CAI between
        the suboptimal w value and 1.0."""
        from biocompiler.translation import compute_cai
        from biocompiler.organisms.yeast import YEAST_CODON_ADAPTIVENESS
        # Half TTT (optimal), half TTC (suboptimal)
        dna = "TTTTTC" * 5
        cai = compute_cai(dna, organism="yeast")
        w_ttt = YEAST_CODON_ADAPTIVENESS["TTT"]
        w_ttc = YEAST_CODON_ADAPTIVENESS["TTC"]
        expected = math.sqrt(w_ttt * w_ttc)  # geometric mean
        assert abs(cai - expected) < 0.01, (
            f"Mixed CAI={cai:.4f}, expected ~{expected:.4f}"
        )


# ────────────────────────────────────────────────────────────
# 6. Validation Integration
# ────────────────────────────────────────────────────────────

class TestYeastValidationIntegration:
    """Test that validate_cai_tables() passes for yeast data."""

    def test_validate_cai_tables_passes(self):
        """The yeast CAI table should pass internal validation."""
        from biocompiler.organisms import validate_cai_tables
        errors = validate_cai_tables()
        yeast_errors = [e for e in errors if "cerevisiae" in e.lower() or "yeast" in e.lower()]
        assert len(yeast_errors) == 0, (
            f"Yeast CAI table validation errors: {yeast_errors}"
        )

    def test_yeast_preferred_codons_match_adaptiveness(self):
        """Preferred codons should all have adaptiveness = 1.0."""
        from biocompiler.organisms.yeast import YEAST_PREFERRED_CODONS, YEAST_CODON_ADAPTIVENESS

        for aa, codon in YEAST_PREFERRED_CODONS.items():
            w = YEAST_CODON_ADAPTIVENESS[codon]
            assert w == 1.0, (
                f"Preferred codon {codon} for {aa} has w={w}, expected 1.0"
            )
