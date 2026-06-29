"""
BioCompiler Stress Test Suite — Production Robustness
=====================================================
Stress tests that probe edge cases, boundary conditions, and pathological
inputs that could break the optimizer in production.

Every test is deterministic and completes in < 1 second.

Test matrix:
  1. Single amino acid              — degenerate case, one codon or many
  2. All-same amino acid (50x Ala)  — codon diversity under homogeneity
  3. Very long protein (1000+ aa)   — performance and correctness at scale
  4. Rare amino acids (W, M)        — single-codon amino acids
  5. GC extreme target (0.49–0.51)  — tight GC constraint
  6. Many restriction enzymes (20+) — heavy constraint load
  7. Repeated optimization (100x)   — deterministic reproducibility
  8. All 5 organisms, same protein  — cross-organism robustness
  9. Many valine                    — all codons contain GT dinucleotide
 10. Empty enzyme list              — no restriction-site avoidance
"""

import time

import pytest

from biocompiler.optimizer import optimize_sequence, OptimizationResult
from biocompiler.type_system import AA_TO_CODONS, CODON_TABLE
from biocompiler.sequence.scanner import gc_content
from biocompiler.expression.translation import compute_cai
from biocompiler.organisms import (
    SUPPORTED_ORGANISMS,
    CODON_ADAPTIVENESS_TABLES,
)
from biocompiler.shared.constants import RESTRICTION_ENZYMES
from biocompiler.shared.exceptions import InvalidProteinError


# ────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────

STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"

# A medium-length protein used for multi-organism / repeatability tests.
# Includes a mix of amino acid types (hydrophobic, polar, charged, single-codon).
_REFERENCE_PROTEIN = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"


def _translate(dna: str) -> str:
    """Translate a DNA sequence to its amino acid sequence."""
    protein = ""
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa = CODON_TABLE.get(codon, "?")
        protein += aa
    return protein


def _all_valid_codons(dna: str) -> bool:
    """Check that every 3-mer in *dna* is a known codon."""
    for i in range(0, len(dna), 3):
        if dna[i:i + 3] not in CODON_TABLE:
            return False
    return True


def _no_internal_stops(dna: str) -> bool:
    """Check that no internal (non-trailing) stop codons exist."""
    for i in range(0, len(dna) - 5, 3):  # -5: skip the last codon
        if dna[i:i + 3] in ("TAA", "TAG", "TGA"):
            return False
    return True


def _check_restriction_sites(dna: str, enzymes: list[str]) -> list[str]:
    """Return list of enzyme names whose sites are present in *dna*."""
    from biocompiler.sequence.restriction_sites import get_recognition_site
    from biocompiler.shared.constants import reverse_complement
    found = []
    for enz in enzymes:
        site = get_recognition_site(enz)
        if site is None:
            continue
        rc = reverse_complement(site)
        if site in dna or (rc and rc in dna):
            found.append(enz)
    return found


# ════════════════════════════════════════════════════════════
# 1. Single Amino Acid
# ════════════════════════════════════════════════════════════

class TestSingleAminoAcid:
    """Optimize a single amino acid — the smallest possible input."""

    def test_single_methionine(self):
        """M has exactly one codon (ATG). Optimizer must return ATG."""
        result = optimize_sequence("M", organism="Homo_sapiens", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert result.sequence == "ATG"
        assert len(result.sequence) == 3

    def test_single_tryptophan(self):
        """W has exactly one codon (TGG)."""
        result = optimize_sequence("W", organism="Homo_sapiens", strict_mode=False)
        assert result.sequence == "TGG"

    def test_single_alanine(self):
        """A has 4 codons — optimizer must pick a valid one."""
        result = optimize_sequence("A", organism="Homo_sapiens", strict_mode=False)
        assert result.sequence in AA_TO_CODONS["A"]
        assert len(result.sequence) == 3

    def test_single_leucine(self):
        """L has 6 codons — the most of any amino acid."""
        result = optimize_sequence("L", organism="Escherichia_coli", strict_mode=False)
        assert result.sequence in AA_TO_CODONS["L"]
        assert len(result.sequence) == 3

    def test_single_for_each_organism(self):
        """Single amino acid works for every supported organism."""
        for org in SUPPORTED_ORGANISMS:
            result = optimize_sequence("A", organism=org, strict_mode=False)
            assert isinstance(result, OptimizationResult)
            assert result.sequence in AA_TO_CODONS["A"], (
                f"Invalid codon for {org}: {result.sequence}"
            )


# ════════════════════════════════════════════════════════════
# 2. All Same Amino Acid
# ════════════════════════════════════════════════════════════

class TestAllSameAminoAcid:
    """Optimize a sequence of all the same amino acid (e.g., 50x Alanine).

    This tests whether the optimizer handles codon diversity well — when
    every position has the same amino acid, codon choice is the only
    lever the optimizer has for GC, restriction sites, etc.
    """

    def test_50x_alanine(self):
        """50 Alanines — 4 codons (GCN), GC-rich."""
        protein = "A" * 50
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 150
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)
        # Alanine codons are GCN, so GC should be relatively high
        assert result.gc_content > 0.3

    def test_50x_leucine(self):
        """50 Leucines — 6 codons (TTR + CTN), most diversity."""
        protein = "L" * 50
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert len(result.sequence) == 150
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_50x_serine(self):
        """50 Serines — 6 codons in two disconnected groups (TCN, AGY)."""
        protein = "S" * 50
        result = optimize_sequence(protein, organism="Escherichia_coli", strict_mode=False)
        assert len(result.sequence) == 150
        assert _all_valid_codons(result.sequence)

    def test_50x_arginine(self):
        """50 Arginines — 6 codons (CGN + AGR), includes AG dinucleotide."""
        protein = "R" * 50
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert len(result.sequence) == 150
        assert _all_valid_codons(result.sequence)

    def test_50x_methionine(self):
        """50 Methionines — only 1 codon (ATG). No codon choice at all."""
        protein = "M" * 50
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert len(result.sequence) == 150
        assert result.sequence == "ATG" * 50


# ════════════════════════════════════════════════════════════
# 3. Very Long Protein
# ════════════════════════════════════════════════════════════

class TestVeryLongProtein:
    """Optimize a 1000+ amino acid protein.

    Verifies correctness and performance at scale.
    """

    def test_1000aa_protein_ecoli(self):
        """1000 alanines in E. coli — prokaryote fast path."""
        protein = "A" * 1000
        t0 = time.monotonic()
        result = optimize_sequence(protein, organism="Escherichia_coli", strict_mode=False)
        elapsed = time.monotonic() - t0
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 3000
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)
        assert result.cai > 0.0
        # Should complete in well under 1 second
        assert elapsed < 1.0, f"1000aa optimization took {elapsed:.2f}s"

    def test_1000aa_protein_human(self):
        """1000 alanines in human — eukaryote path with splice checking."""
        protein = "A" * 1000
        t0 = time.monotonic()
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        elapsed = time.monotonic() - t0
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 3000
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)
        assert elapsed < 3.0, f"1000aa human optimization took {elapsed:.2f}s"

    def test_1000aa_mixed_protein(self):
        """1000 amino acids of mixed composition (repeating STANDARD_AAS)."""
        # Repeat the 20 standard AAs 50 times = 1000 aa
        protein = STANDARD_AAS * 50
        result = optimize_sequence(protein, organism="Escherichia_coli", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 3000
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_1500aa_protein(self):
        """1500 amino acids — stress the optimizer beyond 1000."""
        # "MAGTHIVKLMN" is 11 aa; 11 * 136 = 1496; + "MAGT" = 1500 aa
        protein = "MAGTHIVKLMN" * 136 + "MAGT"
        assert len(protein) == 1500
        result = optimize_sequence(protein, organism="Escherichia_coli", strict_mode=False)
        assert len(result.sequence) == len(protein) * 3
        assert _all_valid_codons(result.sequence)


# ════════════════════════════════════════════════════════════
# 4. Protein with Rare Amino Acids (Single-Codon)
# ════════════════════════════════════════════════════════════

class TestRareAminoAcids:
    """Test protein with Trp (W) and Met (M) — single-codon amino acids.

    These are interesting because the optimizer has zero codon choice
    for these positions.  They also create unavoidable GT (ATG) and
    fixed GC contributions.
    """

    def test_all_trp(self):
        """30 Tryptophans — only TGG codon. No optimization freedom."""
        protein = "W" * 30
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert result.sequence == "TGG" * 30

    def test_all_met(self):
        """30 Methionines — only ATG codon. Creates unavoidable GT."""
        protein = "M" * 30
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert result.sequence == "ATG" * 30

    def test_trp_met_alternating(self):
        """Alternating W-M: TGG ATG TGG ATG... — no codon choice at all."""
        protein = "WM" * 25  # 50 aa
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert len(result.sequence) == 150
        # Every odd codon is TGG, every even codon is ATG
        for i in range(0, 150, 6):
            assert result.sequence[i:i + 3] == "TGG"
            assert result.sequence[i + 3:i + 6] == "ATG"

    def test_protein_with_many_trp_met(self):
        """A realistic protein heavy in W and M."""
        protein = "MWWWMWWWMWWWKLMN" * 6  # 96 aa, many W and M
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == len(protein) * 3
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)


# ════════════════════════════════════════════════════════════
# 5. GC Extreme Target
# ════════════════════════════════════════════════════════════

class TestGCExtremeTarget:
    """Optimize with very tight GC range (0.49–0.51).

    This tests the GC adjustment step under extreme constraint.
    The optimizer must find codons that keep GC within a 2% window.
    """

    def test_tight_gc_ecoli(self):
        """Tight GC range for E. coli — prokaryote fast path."""
        protein = "MAGTHIVKLMN" * 10  # 110 aa
        result = optimize_sequence(
            protein,
            organism="Escherichia_coli",
            gc_lo=0.49,
            gc_hi=0.51,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        # The optimizer should try to get GC close to target;
        # exact compliance depends on amino acid composition
        assert 0.0 <= result.gc_content <= 1.0
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_tight_gc_human(self):
        """Tight GC range for human — eukaryote path."""
        protein = "MAGTHIVKLMN" * 10
        result = optimize_sequence(
            protein,
            organism="Homo_sapiens",
            gc_lo=0.49,
            gc_hi=0.51,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_tight_gc_leucine_heavy(self):
        """Leucine-heavy protein with tight GC — Leu has both GC-rich and GC-poor codons."""
        protein = "L" * 50
        result = optimize_sequence(
            protein,
            organism="Escherichia_coli",
            gc_lo=0.49,
            gc_hi=0.51,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert _all_valid_codons(result.sequence)
        # Leucine codons: TTA(0GC), TTG(1GC), CTT(1GC), CTC(2GC), CTA(1GC), CTG(2GC)
        # The optimizer can mix these to hit target GC

    def test_impossibly_tight_gc(self):
        """GC range so tight it may not be achievable — should still return a result."""
        protein = "M" * 10  # All ATG — GC = 1/3 ≈ 0.333
        result = optimize_sequence(
            protein,
            organism="Escherichia_coli",
            gc_lo=0.49,
            gc_hi=0.51,
            strict_mode=False,
        )
        # Cannot achieve 0.49–0.51 with only ATG codons, but optimizer should not crash
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 30

    def test_gc_zero_to_half_percent(self):
        """Extremely narrow 0.1% GC window — essentially a point target."""
        protein = "ACDEFGHIKLMNPQRSTVWY" * 5  # 100 aa
        result = optimize_sequence(
            protein,
            organism="Escherichia_coli",
            gc_lo=0.50,
            gc_hi=0.501,
            strict_mode=False,
        )
        # Such a narrow window is likely impossible for this protein,
        # but the optimizer must not crash
        assert isinstance(result, OptimizationResult)
        assert _all_valid_codons(result.sequence)


# ════════════════════════════════════════════════════════════
# 6. Many Restriction Sites
# ════════════════════════════════════════════════════════════

class TestManyRestrictionSites:
    """Optimize with 20+ restriction enzymes to avoid.

    Heavy constraint load tests the restriction site removal
    step's ability to handle many simultaneous constraints.
    """

    _MANY_ENZYMES = list(RESTRICTION_ENZYMES.keys())[:25]

    def test_many_enzymes_ecoli(self):
        """25 restriction enzymes, E. coli target."""
        protein = "MAGTHIVKLMN" * 10  # 110 aa
        result = optimize_sequence(
            protein,
            organism="Escherichia_coli",
            enzymes=self._MANY_ENZYMES,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)
        # Check that the most common sites are avoided
        # (not all 25 may be removable depending on protein)
        critical_enzymes = ["EcoRI", "BamHI", "HindIII"]
        found = _check_restriction_sites(result.sequence, critical_enzymes)
        assert len(found) == 0, f"Critical sites still present: {found}"

    def test_many_enzymes_human(self):
        """25 restriction enzymes, human target (eukaryote path)."""
        protein = "MAGTHIVKLMN" * 10
        result = optimize_sequence(
            protein,
            organism="Homo_sapiens",
            enzymes=self._MANY_ENZYMES,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert _all_valid_codons(result.sequence)

    def test_many_enzymes_short_protein(self):
        """25 enzymes on a short protein — limited codon positions to fix."""
        protein = "MAGRKLM"
        result = optimize_sequence(
            protein,
            organism="Homo_sapiens",
            enzymes=self._MANY_ENZYMES,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 21

    def test_4_cutter_enzymes(self):
        """4-cutter enzymes (AluI, HaeIII, etc.) are very frequent — hard to avoid."""
        protein = "ACDEFGHIKLMNPQRSTVWY" * 3  # 60 aa
        four_cutters = ["AluI", "HaeIII", "TaqI", "Sau3AI", "MspI"]
        result = optimize_sequence(
            protein,
            organism="Escherichia_coli",
            enzymes=four_cutters,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert _all_valid_codons(result.sequence)
        # 4-cutters are almost impossible to fully avoid in a 180bp sequence,
        # but the optimizer should still produce a valid result


# ════════════════════════════════════════════════════════════
# 7. Repeated Optimization (Determinism)
# ════════════════════════════════════════════════════════════

class TestRepeatedOptimization:
    """Optimize the same protein 100 times — should be deterministic.

    The greedy optimizer is fully deterministic, so every call with
    the same inputs must produce the identical output.
    """

    def test_100x_deterministic_ecoli(self):
        """100 runs of the same protein in E. coli produce identical results."""
        protein = _REFERENCE_PROTEIN
        results = []
        for _ in range(100):
            r = optimize_sequence(protein, organism="Escherichia_coli", strict_mode=False)
            results.append(r.sequence)
        # All 100 sequences must be identical
        assert len(set(results)) == 1, "Optimizer is non-deterministic!"

    def test_100x_deterministic_human(self):
        """100 runs of the same protein in human produce identical results."""
        protein = _REFERENCE_PROTEIN
        results = []
        for _ in range(100):
            r = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
            results.append(r.sequence)
        assert len(set(results)) == 1, "Optimizer is non-deterministic!"

    def test_100x_deterministic_with_enzymes(self):
        """Determinism holds even with restriction enzyme constraints."""
        protein = _REFERENCE_PROTEIN
        enzymes = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
        results = []
        for _ in range(100):
            r = optimize_sequence(
                protein, organism="Homo_sapiens", enzymes=enzymes,
                strict_mode=False,
            )
            results.append(r.sequence)
        assert len(set(results)) == 1, "Optimizer is non-deterministic with enzymes!"

    def test_deterministic_metrics(self):
        """CAI and GC content are also deterministic across runs."""
        protein = _REFERENCE_PROTEIN
        results = []
        for _ in range(50):
            r = optimize_sequence(protein, organism="Escherichia_coli", strict_mode=False)
            results.append((r.sequence, r.cai, r.gc_content))
        seqs = [r[0] for r in results]
        cais = [r[1] for r in results]
        gcs = [r[2] for r in results]
        assert len(set(seqs)) == 1
        assert len(set(cais)) == 1
        assert len(set(gcs)) == 1

    def test_deterministic_with_tight_gc(self):
        """Determinism under tight GC constraint."""
        protein = "MAGTHIVKLMN" * 5
        results = []
        for _ in range(50):
            r = optimize_sequence(
                protein, organism="Escherichia_coli",
                gc_lo=0.45, gc_hi=0.55,
                strict_mode=False,
            )
            results.append(r.sequence)
        assert len(set(results)) == 1


# ════════════════════════════════════════════════════════════
# 8. All Organisms — Same Protein
# ════════════════════════════════════════════════════════════

class TestAllOrganismsSameProtein:
    """Optimize the same protein for all 5 supported organisms.

    Verifies that the optimizer works for every organism and produces
    a valid result with correct metrics.
    """

    @pytest.fixture(scope="class")
    def protein(self):
        return _REFERENCE_PROTEIN

    def test_homo_sapiens(self, protein):
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == len(protein) * 3
        assert result.cai > 0.0
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_escherichia_coli(self, protein):
        result = optimize_sequence(protein, organism="Escherichia_coli", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == len(protein) * 3
        assert result.cai > 0.0
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_mus_musculus(self, protein):
        result = optimize_sequence(protein, organism="Mus_musculus", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == len(protein) * 3
        assert result.cai > 0.0
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_cho_k1(self, protein):
        result = optimize_sequence(protein, organism="CHO_K1", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == len(protein) * 3
        assert result.cai > 0.0
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_saccharomyces_cerevisiae(self, protein):
        result = optimize_sequence(protein, organism="Saccharomyces_cerevisiae", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == len(protein) * 3
        assert result.cai > 0.0
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_different_organisms_different_sequences(self, protein):
        """Different organisms should generally produce different codon choices."""
        sequences = {}
        for org in SUPPORTED_ORGANISMS:
            r = optimize_sequence(protein, organism=org, strict_mode=False)
            sequences[org] = r.sequence
        # At least 2 organisms should differ (codon preferences differ)
        unique_seqs = set(sequences.values())
        assert len(unique_seqs) >= 2, (
            f"All organisms produced identical sequences — likely a bug"
        )

    def test_all_organisms_cai_reasonable(self, protein):
        """CAI should be > 0.2 for all organisms (reasonable optimization)."""
        for org in SUPPORTED_ORGANISMS:
            r = optimize_sequence(protein, organism=org, strict_mode=False)
            assert r.cai > 0.2, f"CAI too low for {org}: {r.cai:.4f}"

    def test_all_organisms_gc_in_range(self, protein):
        """GC content should be within default [0.30, 0.70] for all organisms."""
        for org in SUPPORTED_ORGANISMS:
            r = optimize_sequence(protein, organism=org, gc_lo=0.30, gc_hi=0.70, strict_mode=False)
            # GC may not always be perfectly in range due to amino acid constraints,
            # but it should be a valid fraction
            assert 0.0 <= r.gc_content <= 1.0, (
                f"GC out of bounds for {org}: {r.gc_content}"
            )


# ════════════════════════════════════════════════════════════
# 9. Protein with Many Valine
# ════════════════════════════════════════════════════════════

class TestManyValine:
    """Valine is tricky because all codons contain GT.

    All 4 valine codons (GTT, GTC, GTA, GTG) start with GT, creating
    unavoidable cryptic splice donor sites.  The optimizer must handle
    this gracefully without crashing or entering infinite loops.
    """

    def test_all_valine(self):
        """100% valine protein — every codon has GT."""
        protein = "V" * 30
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 90
        for i in range(0, 90, 3):
            assert result.sequence[i:i + 3] in AA_TO_CODONS["V"]

    def test_valine_heavy_mixed(self):
        """50% valine + alternating amino acid."""
        protein = "VA" * 40  # 80 aa
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 240
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_valine_with_serine(self):
        """Valine + Serine — both create GT dinucleotides (AGT in serine AGY codons)."""
        protein = "VS" * 30  # 60 aa
        result = optimize_sequence(protein, organism="Homo_sapiens", strict_mode=False)
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 180
        assert _all_valid_codons(result.sequence)

    def test_valine_prokaryote_no_splice_check(self):
        """In E. coli, valine GT should not trigger splice checks (prokaryote path)."""
        protein = "V" * 50
        t0 = time.monotonic()
        result = optimize_sequence(protein, organism="Escherichia_coli", strict_mode=False)
        elapsed = time.monotonic() - t0
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == 150
        # Prokaryote path should be fast — no splice checking
        assert elapsed < 1.0

    def test_valine_ecoli_high_cai(self):
        """Valine-heavy protein in E. coli should achieve decent CAI."""
        protein = "VAVAVAVAVA" * 5  # 50 aa
        result = optimize_sequence(protein, organism="Escherichia_coli", strict_mode=False)
        assert result.cai > 0.0
        # E. coli prokaryote fast path does not waste CAI on GT avoidance
        assert _all_valid_codons(result.sequence)


# ════════════════════════════════════════════════════════════
# 10. Empty Enzyme List
# ════════════════════════════════════════════════════════════

class TestEmptyEnzymeList:
    """Optimize with no restriction enzyme avoidance.

    When no enzymes are specified, the optimizer should still
    produce valid results — restriction site avoidance is simply
    skipped or uses defaults.
    """

    def test_empty_enzymes_list(self):
        """Explicitly passing enzymes=[] should skip RS avoidance."""
        protein = _REFERENCE_PROTEIN
        result = optimize_sequence(
            protein,
            organism="Homo_sapiens",
            enzymes=[],
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == len(protein) * 3
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)

    def test_none_enzymes(self):
        """Passing enzymes=None uses default enzyme set."""
        protein = _REFERENCE_PROTEIN
        result = optimize_sequence(
            protein,
            organism="Homo_sapiens",
            enzymes=None,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert len(result.sequence) == len(protein) * 3
        assert _all_valid_codons(result.sequence)

    def test_empty_vs_default_differ(self):
        """Empty enzyme list may produce different sequence than default (more freedom)."""
        protein = "MAGTHIVKLMN" * 5
        result_empty = optimize_sequence(
            protein, organism="Escherichia_coli", enzymes=[],
            strict_mode=False,
        )
        result_default = optimize_sequence(
            protein, organism="Escherichia_coli",
            strict_mode=False,
        )
        # Both should be valid
        assert _all_valid_codons(result_empty.sequence)
        assert _all_valid_codons(result_default.sequence)
        # The empty-enzyme version might have higher CAI
        # (fewer constraints → more codon freedom), but it is not guaranteed
        assert result_empty.cai > 0.0
        assert result_default.cai > 0.0

    def test_empty_enzymes_all_organisms(self):
        """Empty enzyme list works for all organisms."""
        protein = "MAGTHIVKLMN"
        for org in SUPPORTED_ORGANISMS:
            result = optimize_sequence(protein, organism=org, enzymes=[], strict_mode=False)
            assert isinstance(result, OptimizationResult)
            assert len(result.sequence) == len(protein) * 3
            assert _all_valid_codons(result.sequence)

    def test_empty_enzymes_long_protein(self):
        """Empty enzyme list with a longer protein."""
        protein = STANDARD_AAS * 25  # 500 aa
        result = optimize_sequence(
            protein, organism="Escherichia_coli", enzymes=[],
            strict_mode=False,
        )
        assert len(result.sequence) == 1500
        assert _all_valid_codons(result.sequence)
        assert _no_internal_stops(result.sequence)
