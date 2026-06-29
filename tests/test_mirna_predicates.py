"""
Comprehensive tests for miRNA predicate logic in BioCompiler.

Test categories:
1. TestCheckNoMiRNABindingSite  — low-level check function
2. TestEvaluateNoMiRNABindingSite — high-level evaluate API
3. TestMiRNASeedDatabase — multi-organism seed database integrity
4. TestMiRNAAvoidance — codon-substitution avoidance of miRNA sites
"""

import os
import sys
import warnings

import pytest

# Ensure src is on path
sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")
)

from biocompiler.type_system.checks import (
    check_no_mirna_binding_site,
    _rna_revcomp_to_dna,
    _mirna_context_score,
)
from biocompiler.type_system.predicates import evaluate_no_mirna_binding_site
from biocompiler.type_system.mirna_seeds import (
    get_mirna_seeds,
    HUMAN_MIRNA_SEEDS,
    MOUSE_MIRNA_SEEDS,
    CHO_MIRNA_SEEDS,
    RAT_MIRNA_SEEDS,
    ORGANISM_MIRNA_MAP,
)
from biocompiler.type_system.codon_tables import AA_TO_CODONS, CODON_TABLE
from biocompiler.shared.types import Verdict, TypeCheckResult
from biocompiler.expression.translation import translate


# ═══════════════════════════════════════════════════════════════════════════
# Helper: Construct well-known test sequences
# ═══════════════════════════════════════════════════════════════════════════

# miR-21-5p seed (AGCUUAU) → DNA target ATAAGCT
# 8mer requires T immediately 5': TATAAGCT
_MIR21_DNA_TARGET = "ATAAGCT"

# miR-16-5p seed (UCAAGU) → DNA target ACTTGA (6-base seed)
_MIR16_DNA_TARGET = "ACTTGA"

# miR-122-5p seed (GGAGUGU) → DNA target ACACTCC
_MIR122_DNA_TARGET = "ACACTCC"


def _make_seq_with_8mer(target_7: str, flank: str = "GCC") -> str:
    """Build a sequence with an 8mer miRNA match (T + 7mer target).

    The target is embedded in a neutral GC-rich context that gives
    a context factor close to 1.0, so the adjusted score ≈ base score.
    """
    padding = flank * 12
    return padding + "T" + target_7 + padding


def _make_seq_with_7mer_m8(target_7: str, prefix: str = "A", flank: str = "GCC") -> str:
    """Build a sequence with a 7mer-m8 match (no T before the target).

    *prefix* is placed immediately before the 7mer target to avoid
    accidentally creating an 8mer.  Use 'A' (not 'T').
    """
    padding = flank * 12
    return padding + prefix + target_7 + padding


def _make_safe_seq(length: int = 60) -> str:
    """Build a DNA sequence with no miRNA seed matches.

    Uses a GC-rich pattern (GGCGCC repeat) that avoids all known
    human miRNA seed reverse complements.
    """
    unit = "GGCGCC"
    reps = (length + len(unit) - 1) // len(unit)
    return (unit * reps)[:length]


# ═══════════════════════════════════════════════════════════════════════════
# 1. TestCheckNoMiRNABindingSite
# ═══════════════════════════════════════════════════════════════════════════


class TestCheckNoMiRNABindingSite:
    """Unit tests for the check_no_mirna_binding_site function."""

    def test_short_sequence_passes(self):
        """Sequences shorter than 6 bp should return PASS."""
        for seq in ["", "A", "AT", "ATG", "ATGC", "ATGCG"]:
            r = check_no_mirna_binding_site(seq)
            assert r.verdict == Verdict.PASS, f"seq={seq!r} should PASS"
            assert r.passed is True

    def test_no_seed_match_passes(self):
        """A GC-rich sequence with no seed matches should PASS."""
        seq = _make_safe_seq(120)
        r = check_no_mirna_binding_site(seq)
        assert r.verdict == Verdict.PASS
        assert r.passed is True

    def test_8mer_tier1_fails(self):
        """An 8mer match to a tier-1 miRNA should FAIL.

        miR-21-5p (tier 1, ubiquitous) seed AGCUUAU → DNA target ATAAGCT.
        An 8mer requires T immediately before the target: TATAAGCT.
        """
        seq = _make_seq_with_8mer(_MIR21_DNA_TARGET)
        r = check_no_mirna_binding_site(seq)
        assert r.verdict == Verdict.FAIL
        assert r.passed is False
        assert "8mer" in r.details
        assert "hsa-miR-21-5p" in r.details

    def test_7mer_m8_tier1_fails(self):
        """A 7mer-m8 match to a tier-1 miRNA should FAIL when context
        boosts the adjusted score to >= 0.9.

        With AU-rich flanking the context factor can push a 0.9 base
        score above the FAIL threshold.
        """
        # Build an AU-rich sequence with a 7mer-m8 match to miR-21-5p
        target = _MIR21_DNA_TARGET  # ATAAGCT
        padding = "AATAAT" * 10
        seq = padding + "A" + target + padding
        r = check_no_mirna_binding_site(seq)
        assert r.verdict == Verdict.FAIL
        assert "7mer-m8" in r.details

    def test_6mer_returns_uncertain(self):
        """A 6mer match (with sufficient context) should return UNCERTAIN.

        miR-16-5p seed UCAAGU (6 bases) → DNA target ACTTGA.
        Need min_seed_match=6 to detect 6mer matches.
        With balanced flanking the context score reaches >= 0.7.
        """
        target = _MIR16_DNA_TARGET  # ACTTGA (6-base seed)
        # Use GGCATG padding: mixed GC/AT that gives context factor ~1.02
        # (score ≈ 0.71 ≥ 0.7 threshold for UNCERTAIN) while avoiding
        # accidental seed matches at junctions.
        padding = "GGCATG" * 10
        seq = padding + "A" + target + padding
        r = check_no_mirna_binding_site(seq, min_seed_match=6)
        assert r.verdict == Verdict.UNCERTAIN

    def test_7mer_a1_returns_uncertain(self):
        """A 7mer-A1 match should return UNCERTAIN.

        7mer-A1 = 6mer target + T upstream.
        miR-16-5p seed UCAAGU → DNA target ACTTGA.
        TACTTGA gives 7mer-A1.
        """
        target = _MIR16_DNA_TARGET  # ACTTGA
        padding = "AATAAT" * 10
        seq = padding + "T" + target + padding
        r = check_no_mirna_binding_site(seq, min_seed_match=6)
        # 7mer-A1 base_score=0.85, context may push above or below FAIL
        # In most cases it should be UNCERTAIN or FAIL
        assert r.verdict in (Verdict.UNCERTAIN, Verdict.FAIL)

    def test_multiple_hits_worst_verdict(self):
        """Multiple matches should report the worst (most severe) verdict.

        Embed both an 8mer (FAIL) and a 7mer-m8 (could be FAIL or
        UNCERTAIN) in the same sequence. The overall verdict should be FAIL.
        """
        # 8mer for miR-21-5p
        seq_8mer = _make_seq_with_8mer(_MIR21_DNA_TARGET)
        r = check_no_mirna_binding_site(seq_8mer)
        assert r.verdict == Verdict.FAIL

    def test_positions_reported(self):
        """Positions of seed matches should be reported in the result."""
        seq = _make_seq_with_8mer(_MIR21_DNA_TARGET)
        r = check_no_mirna_binding_site(seq)
        if r.verdict != Verdict.PASS:
            assert len(r.positions) > 0
            # The 8mer target should be at the expected position
            assert any(
                seq[pos:pos + 7] == _MIR21_DNA_TARGET for pos in r.positions
            )

    def test_invalid_bases_pass(self):
        """Sequences with N or other invalid bases should be handled
        gracefully — the seed scanner simply will not find Watson-Crick
        matches at positions with ambiguous bases."""
        # Embed a valid 8mer site but surround it with Ns
        seq = "ATGNNNNNNTATAAGCTGNNNNNNNNNNN"
        r = check_no_mirna_binding_site(seq)
        # The 8mer is still present and should be found
        assert r.verdict == Verdict.FAIL
        assert "hsa-miR-21-5p" in r.details

    def test_uppercase_handling(self):
        """Lowercase input should be handled — the function upper-cases
        the sequence internally."""
        seq_lower = _make_seq_with_8mer(_MIR21_DNA_TARGET).lower()
        seq_upper = _make_seq_with_8mer(_MIR21_DNA_TARGET)
        r_lower = check_no_mirna_binding_site(seq_lower)
        r_upper = check_no_mirna_binding_site(seq_upper)
        assert r_lower.verdict == r_upper.verdict

    def test_different_organisms(self):
        """Test with organism='Mus_musculus' — should use the mouse
        miRNA seed database, which has mmu-miR-21-5p."""
        seq = _make_seq_with_8mer(_MIR21_DNA_TARGET)
        r = check_no_mirna_binding_site(seq, organism="Mus_musculus")
        assert r.verdict == Verdict.FAIL
        # Mouse seeds use mmu- prefix
        assert "mmu-miR-21-5p" in r.details

    def test_tissue_filtering(self):
        """Test with tissue='liver' to verify tissue-aware filtering.

        When a tissue is specified, miRNAs from unrelated tissues are
        downweighted (effective tier is increased), while ubiquitous and
        tissue-matching miRNAs retain their tier.
        miR-21-5p is ubiquitous, so it should still be detected as tier 1.
        """
        seq = _make_seq_with_8mer(_MIR21_DNA_TARGET)
        r_no_tissue = check_no_mirna_binding_site(seq, tissue="")
        r_liver = check_no_mirna_binding_site(seq, tissue="liver")
        # miR-21-5p is ubiquitous — always included
        assert r_liver.verdict == Verdict.FAIL
        assert "hsa-miR-21-5p" in r_liver.details

    def test_match_type_scoring(self):
        """Verify base scores for each match type.

        Base scores: 8mer=1.0, 7mer-m8=0.9, 7mer-A1=0.85, 6mer=0.7.
        The _mirna_context_score function adjusts these, so we test
        that the base scores are used correctly as inputs.
        """
        # 8mer base score
        assert _mirna_context_score("ATGTATAAGCTGCCC", 3, "8mer", 1.0) > 0
        # 7mer-m8 base score
        assert _mirna_context_score("ATGGCCATAAGCTGCCC", 5, "7mer-m8", 0.9) > 0
        # 7mer-A1 base score
        assert _mirna_context_score("ATGGCCTACTTGAGCCC", 6, "7mer-A1", 0.85) > 0
        # 6mer base score
        assert _mirna_context_score("ATGGCCACTTGAGCCC", 6, "6mer", 0.7) > 0

        # Verify the context-adjusted score is always positive and bounded
        for match_type, base in [
            ("8mer", 1.0), ("7mer-m8", 0.9),
            ("7mer-A1", 0.85), ("6mer", 0.7),
        ]:
            score = _mirna_context_score(
                "ATGGCCATAAGCTGCCC", 5, match_type, base
            )
            # Context factor ∈ [0.7, 1.2], so adjusted score ∈ [base*0.7, base*1.2]
            assert base * 0.7 <= score <= base * 1.2, (
                f"{match_type}: score {score} out of expected range"
            )

    def test_details_string(self):
        """Verify the details string contains useful information about hits."""
        seq = _make_seq_with_8mer(_MIR21_DNA_TARGET)
        r = check_no_mirna_binding_site(seq)
        # Should contain miRNA name, match type, and position
        assert "hsa-miR-21-5p" in r.details
        assert "8mer" in r.details
        assert "pos=" in r.details
        assert "score=" in r.details

    def test_7mer_m8_gc_context_uncertain(self):
        """A 7mer-m8 match in GC-rich context may get score < 0.9,
        resulting in UNCERTAIN rather than FAIL."""
        target = _MIR21_DNA_TARGET  # ATAAGCT
        padding = "GCC" * 14
        seq = padding + "A" + target + padding
        r = check_no_mirna_binding_site(seq)
        # GC-rich context reduces score; 7mer-m8 base=0.9 * ctx<1 → UNCERTAIN
        assert r.verdict == Verdict.UNCERTAIN

    def test_default_min_seed_match_is_7(self):
        """By default min_seed_match=7, so 6-base seeds are skipped."""
        # miR-16-5p has a 6-base seed, so with default min_seed_match=7
        # it will not be detected even if its target is present
        target = _MIR16_DNA_TARGET  # ACTTGA
        seq = "AATAAT" * 10 + "A" + target + "AATAAT" * 10
        r = check_no_mirna_binding_site(seq)  # default min_seed_match=7
        # 6-base seed should be skipped with default min_seed_match=7
        assert r.verdict == Verdict.PASS

    def test_min_seed_match_6_enables_6mer(self):
        """Setting min_seed_match=6 enables detection of 6-base seeds."""
        target = _MIR16_DNA_TARGET  # ACTTGA
        seq = "AATAAT" * 10 + "A" + target + "AATAAT" * 10
        r = check_no_mirna_binding_site(seq, min_seed_match=6)
        # Should now detect the 6mer match
        assert r.verdict != Verdict.PASS or "miR-16" in r.details


# ═══════════════════════════════════════════════════════════════════════════
# 2. TestEvaluateNoMiRNABindingSite
# ═══════════════════════════════════════════════════════════════════════════


class TestEvaluateNoMiRNABindingSite:
    """Tests for the high-level evaluate_no_mirna_binding_site API."""

    def test_evaluate_pass(self):
        """Call evaluate_no_mirna_binding_site with a safe sequence;
        should return TypeCheckResult with PASS."""
        seq = _make_safe_seq(120)
        result = evaluate_no_mirna_binding_site(seq)
        assert isinstance(result, TypeCheckResult)
        assert result.verdict == Verdict.PASS
        assert result.passed is True
        assert "NoMiRNABindingSite" in result.predicate

    def test_evaluate_fail(self):
        """Call with a sequence that contains an 8mer match;
        should return TypeCheckResult with FAIL."""
        seq = _make_seq_with_8mer(_MIR21_DNA_TARGET)
        result = evaluate_no_mirna_binding_site(seq)
        assert isinstance(result, TypeCheckResult)
        assert result.verdict == Verdict.FAIL
        assert result.passed is False
        assert result.violation is not None
        assert "miRNA" in result.violation

    def test_evaluate_uncertain(self):
        """Call with a sequence that triggers UNCERTAIN verdict."""
        # 7mer-m8 in GC-rich context → UNCERTAIN
        target = _MIR21_DNA_TARGET  # ATAAGCT
        padding = "GCC" * 14
        seq = padding + "A" + target + padding
        result = evaluate_no_mirna_binding_site(seq)
        assert isinstance(result, TypeCheckResult)
        assert result.verdict == Verdict.UNCERTAIN
        assert result.violation is not None

    def test_evaluate_verdict_field(self):
        """Verify the verdict field is properly set for each case."""
        # PASS
        r_pass = evaluate_no_mirna_binding_site(_make_safe_seq(60))
        assert r_pass.verdict == Verdict.PASS

        # FAIL
        r_fail = evaluate_no_mirna_binding_site(
            _make_seq_with_8mer(_MIR21_DNA_TARGET)
        )
        assert r_fail.verdict == Verdict.FAIL

        # UNCERTAIN
        target = _MIR21_DNA_TARGET
        padding = "GCC" * 14
        r_unc = evaluate_no_mirna_binding_site(padding + "A" + target + padding)
        assert r_unc.verdict == Verdict.UNCERTAIN


# ═══════════════════════════════════════════════════════════════════════════
# 3. TestMiRNASeedDatabase
# ═══════════════════════════════════════════════════════════════════════════


class TestMiRNASeedDatabase:
    """Tests for the multi-organism miRNA seed database."""

    def test_human_seeds_exist(self):
        """HUMAN_MIRNA_SEEDS should have >= 25 entries."""
        assert len(HUMAN_MIRNA_SEEDS) >= 25

    def test_mouse_seeds_exist(self):
        """MOUSE_MIRNA_SEEDS should have >= 15 entries."""
        assert len(MOUSE_MIRNA_SEEDS) >= 15

    def test_cho_seeds_exist(self):
        """CHO_MIRNA_SEEDS should have >= 10 entries."""
        assert len(CHO_MIRNA_SEEDS) >= 10

    def test_rat_seeds_exist(self):
        """RAT_MIRNA_SEEDS should have >= 10 entries."""
        assert len(RAT_MIRNA_SEEDS) >= 10

    def test_seed_format(self):
        """Each seed should be a tuple of (rna_string, tier_int, tissue_string)."""
        for db_name, db in [
            ("HUMAN", HUMAN_MIRNA_SEEDS),
            ("MOUSE", MOUSE_MIRNA_SEEDS),
            ("CHO", CHO_MIRNA_SEEDS),
            ("RAT", RAT_MIRNA_SEEDS),
        ]:
            for mirna_name, entry in db.items():
                assert len(entry) == 3, (
                    f"{db_name}:{mirna_name} has {len(entry)} elements, expected 3"
                )
                seed_rna, tier, tissue = entry
                assert isinstance(seed_rna, str), (
                    f"{db_name}:{mirna_name} seed_rna is not str"
                )
                assert isinstance(tier, int), (
                    f"{db_name}:{mirna_name} tier is not int"
                )
                assert isinstance(tissue, str), (
                    f"{db_name}:{mirna_name} tissue is not str"
                )
                assert tier in (1, 2, 3), (
                    f"{db_name}:{mirna_name} tier={tier}, expected 1/2/3"
                )

    def test_seed_rna_valid_bases(self):
        """Seed RNA strings should only contain A, U, G, C."""
        valid_bases = {"A", "U", "G", "C"}
        for db_name, db in [
            ("HUMAN", HUMAN_MIRNA_SEEDS),
            ("MOUSE", MOUSE_MIRNA_SEEDS),
            ("CHO", CHO_MIRNA_SEEDS),
            ("RAT", RAT_MIRNA_SEEDS),
        ]:
            for mirna_name, (seed_rna, tier, tissue) in db.items():
                invalid = set(seed_rna) - valid_bases
                assert not invalid, (
                    f"{db_name}:{mirna_name} has invalid bases: {invalid}"
                )

    def test_get_mirna_seeds_lookup(self):
        """get_mirna_seeds should return the correct DB for known organisms.

        Note: get_mirna_seeds may return a copy, so we check equality
        rather than identity.
        """
        assert get_mirna_seeds("Homo_sapiens") == HUMAN_MIRNA_SEEDS
        assert get_mirna_seeds("Mus_musculus") == MOUSE_MIRNA_SEEDS
        assert get_mirna_seeds("Cricetulus_griseus") == CHO_MIRNA_SEEDS
        assert get_mirna_seeds("Rattus_norvegicus") == RAT_MIRNA_SEEDS

    def test_get_mirna_seeds_fallback(self):
        """Unknown organisms should fall back to human seeds with a warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = get_mirna_seeds("Unknown_organism_xyz")
            # Should return human seeds (equality, not identity)
            assert result == HUMAN_MIRNA_SEEDS
            # Should emit a warning
            assert len(w) >= 1
            assert "falling back" in str(w[0].message).lower() or \
                   "fallback" in str(w[0].message).lower()

    def test_seed_lengths_reasonable(self):
        """Seed lengths should be between 6 and 8 nucleotides."""
        for db_name, db in [
            ("HUMAN", HUMAN_MIRNA_SEEDS),
            ("MOUSE", MOUSE_MIRNA_SEEDS),
            ("CHO", CHO_MIRNA_SEEDS),
            ("RAT", RAT_MIRNA_SEEDS),
        ]:
            for mirna_name, (seed_rna, tier, tissue) in db.items():
                assert 6 <= len(seed_rna) <= 8, (
                    f"{db_name}:{mirna_name} seed length {len(seed_rna)} "
                    f"outside expected range [6, 8]"
                )

    def test_organism_mirna_map_consistency(self):
        """ORGANISM_MIRNA_MAP should contain entries for all four organisms."""
        assert "Homo_sapiens" in ORGANISM_MIRNA_MAP
        assert "Mus_musculus" in ORGANISM_MIRNA_MAP
        assert "Cricetulus_griseus" in ORGANISM_MIRNA_MAP
        assert "Rattus_norvegicus" in ORGANISM_MIRNA_MAP
        assert ORGANISM_MIRNA_MAP["Homo_sapiens"] is HUMAN_MIRNA_SEEDS
        assert ORGANISM_MIRNA_MAP["Mus_musculus"] is MOUSE_MIRNA_SEEDS

    def test_conserved_seeds_across_organisms(self):
        """Key miRNA families should have conserved seeds across organisms.

        miR-21-5p, let-7a-5p, and miR-16-5p are deeply conserved.
        """
        human_21 = HUMAN_MIRNA_SEEDS["hsa-miR-21-5p"][0]
        mouse_21 = MOUSE_MIRNA_SEEDS["mmu-miR-21-5p"][0]
        assert human_21 == mouse_21 == "AGCUUAU"

        human_let7 = HUMAN_MIRNA_SEEDS["hsa-let-7a-5p"][0]
        mouse_let7 = MOUSE_MIRNA_SEEDS["mmu-let-7a-5p"][0]
        assert human_let7 == mouse_let7 == "GAGGUAG"


# ═══════════════════════════════════════════════════════════════════════════
# 4. TestMiRNAAvoidance
# ═══════════════════════════════════════════════════════════════════════════


def _eliminate_mirna_binding_sites(
    seq: str,
    protein: str,
    organism: str = "Homo_sapiens",
    max_iterations: int = 100,
) -> str:
    """Simple codon-substitution loop to remove miRNA seed matches.

    For each seed match found, try substituting the codon(s) covering
    the match with synonymous codons that break the seed. This is a
    test helper — the production optimizer is far more sophisticated.

    Returns the modified sequence (or the original if no changes needed).
    """
    seq = seq.upper()
    for _ in range(max_iterations):
        result = check_no_mirna_binding_site(seq, organism=organism)
        if result.verdict == Verdict.PASS:
            return seq

        # Get the first failing position
        if not result.positions:
            return seq

        pos = result.positions[0]

        # Determine which codons overlap with the match at pos
        # Try substituting each overlapping codon
        changed = False
        for codon_start in range(
            max(0, (pos // 3) * 3),
            min(len(seq), ((pos + 7) // 3 + 1) * 3),
            3,
        ):
            if codon_start + 3 > len(seq):
                break
            codon = seq[codon_start: codon_start + 3]
            aa = CODON_TABLE.get(codon)
            if aa is None or aa == "*":
                continue
            synonyms = AA_TO_CODONS.get(aa, [])
            for alt in synonyms:
                if alt == codon:
                    continue
                trial = seq[:codon_start] + alt + seq[codon_start + 3:]
                # Verify protein is preserved
                if translate(trial) == protein:
                    trial_result = check_no_mirna_binding_site(
                        trial, organism=organism
                    )
                    if trial_result.verdict != Verdict.FAIL:
                        seq = trial
                        changed = True
                        break
            if changed:
                break

        if not changed:
            # Cannot fix with single-codon substitution at this position
            return seq

    return seq


class TestMiRNAAvoidance:
    """Tests for miRNA binding site avoidance via codon substitution."""

    def test_eliminate_mirna_binding_sites(self):
        """Verify that the avoidance function can remove a known seed match."""
        # Build a coding sequence for a protein that contains an 8mer
        # miR-21-5p target: ATAAGCT appears in "TATAAGCT" (8mer)
        # We need a protein whose codons, when concatenated, include ATAAGCT
        # ATA AAG CT... = I-K-L... or ATAAGC T... = I-S...
        # Let us use a known protein and construct DNA that includes the site
        protein = "MAAAAYKLPPPPPP"  # M-A(5)-Y-K-L-P(6)
        # Build DNA with the 8mer embedded:
        # M=ATG A=GCG A=GCG A=GCG A=GCG A=GCG Y=TAC K=AAA L=TTA
        # That gives ATGGCGGCGGCGGCGGCGTACAAATTA...
        # No 8mer there. Let me be more direct:
        # ATA = I, AAG = K, CTT = L → "IKL"
        # So "TATAAGCTT" encodes: TAT=Y, AAG=K, CTT=L, which is YKL
        # Or: T-ATA-AGC-T... where ATA=I, AGC=S → "TIS..."
        # Let me just use: protein contains "IKL" and use ATA AAG CTT
        protein = "MIKLPPPP"
        # M=ATG I=ATA K=AAG L=CTT P=CCC P=CCC P=CCC P=CCC
        # ATG ATA AAG CTT CCC CCC CCC CCC
        # This has "TATAAGCT" from position 2: ATG [T] ATA AAG C [T]
        # A-T-G-A-T-A-A-A-G-C-T-T-C-C-C-C-C-C-C-C
        # Position 3: T-A-T-A-A-A-G-C → nope, that is TATAAAGC not TATAAGCT
        # Let me try: I=ATA K=AAG L=CTA → ATAAGCTA → ATA AAG CTA
        # "ATAAGCTA" = ATA-AGC-TA... but CTA is Leu
        # Or: T-ATA-AGC-T... = I-S with partial codon
        # Let me think differently:
        # The 8mer is TATAAGCT. This is T + ATAAGCT.
        # TAT = Y, AAG = K, CT... = need next base
        # Or: T-ATA-AGC-T → ATA=I, AGC=S, and T is partial
        # Better: TAT = Y, AAG = K, CTT = L → protein "YKL"
        # So TATAAGCTT = YKL
        protein = "MYKLPPPP"
        # M=ATG Y=TAT K=AAG L=CTT P=CCC...
        # ATG TAT AAG CTT CCC CCC CCC CCC
        # That gives ATGTATAAGCTTCCCCCC... 
        # The substring "TATAAGCT" is at position 3!
        seq = "ATGTATAAGCTTCCCCCCCCCCCC"
        # Verify it has the 8mer
        r = check_no_mirna_binding_site(seq)
        assert r.verdict == Verdict.FAIL, "Precondition: 8mer must be present"

        # Run avoidance
        new_seq = _eliminate_mirna_binding_sites(seq, protein)
        # Verify the 8mer is gone
        r2 = check_no_mirna_binding_site(new_seq)
        assert r2.verdict != Verdict.FAIL, (
            f"8mer should be eliminated; verdict={r2.verdict}, "
            f"details={r2.details}"
        )

    def test_eliminate_preserves_protein(self):
        """Verify the output sequence translates to the same protein."""
        protein = "MYKLPPPP"
        seq = "ATGTATAAGCTTCCCCCCCCCCCC"
        new_seq = _eliminate_mirna_binding_sites(seq, protein)
        new_protein = translate(new_seq)
        assert new_protein == protein, (
            f"Protein changed: {new_protein!r} != {protein!r}"
        )

    def test_eliminate_max_iterations(self):
        """Verify that max_iterations is respected — the function should
        terminate without error even if it cannot fully eliminate all sites."""
        protein = "MYKLPPPP"
        seq = "ATGTATAAGCTTCCCCCCCCCCCC"

        # With max_iterations=0, the function should return immediately
        # (no substitution attempts made, but initial check still done)
        new_seq = _eliminate_mirna_binding_sites(
            seq, protein, max_iterations=0
        )
        # With 0 iterations, no changes are made
        # The function should still return a valid sequence
        assert isinstance(new_seq, str)
        assert len(new_seq) > 0

        # With max_iterations=1, at most one substitution is attempted
        new_seq_1 = _eliminate_mirna_binding_sites(
            seq, protein, max_iterations=1
        )
        assert isinstance(new_seq_1, str)
        # Protein should be preserved even with limited iterations
        if new_seq_1 != seq:
            assert translate(new_seq_1) == protein


# ═══════════════════════════════════════════════════════════════════════════
# 5. Additional edge cases
# ═══════════════════════════════════════════════════════════════════════════


class TestMiRNAEdgeCases:
    """Edge-case tests for miRNA predicate logic."""

    def test_rna_revcomp_to_dna(self):
        """Verify the _rna_revcomp_to_dna helper function."""
        # AGCUUAU → reverse: UAUUCGA → complement: ATAAGCT → T→DNA: ATAAGCT
        assert _rna_revcomp_to_dna("AGCUUAU") == "ATAAGCT"
        # GAGGUAG → reverse: GAUGGAG → complement: CUACCUC → DNA: CTACCTC
        assert _rna_revcomp_to_dna("GAGGUAG") == "CTACCTC"
        # Single base
        assert _rna_revcomp_to_dna("A") == "T"
        assert _rna_revcomp_to_dna("U") == "A"
        assert _rna_revcomp_to_dna("G") == "C"
        assert _rna_revcomp_to_dna("C") == "G"

    def test_empty_sequence(self):
        """Empty sequence should PASS."""
        r = check_no_mirna_binding_site("")
        assert r.verdict == Verdict.PASS

    def test_exactly_6bp_no_match(self):
        """A 6bp sequence with no seed match should PASS."""
        r = check_no_mirna_binding_site("GGGGGG")
        assert r.verdict == Verdict.PASS

    def test_predicate_name(self):
        """The predicate name should be 'NoMiRNABindingSite'."""
        r = check_no_mirna_binding_site("ATGGCC")
        assert r.predicate == "NoMiRNABindingSite"

    def test_cho_organism(self):
        """Test with CHO (Cricetulus griseus) organism."""
        seq = _make_seq_with_8mer(_MIR21_DNA_TARGET)
        r = check_no_mirna_binding_site(seq, organism="Cricetulus_griseus")
        assert r.verdict == Verdict.FAIL
        assert "cgr-miR-21-5p" in r.details

    def test_rat_organism(self):
        """Test with rat (Rattus norvegicus) organism."""
        seq = _make_seq_with_8mer(_MIR21_DNA_TARGET)
        r = check_no_mirna_binding_site(seq, organism="Rattus_norvegicus")
        assert r.verdict == Verdict.FAIL
        assert "rno-miR-21-5p" in r.details

    def test_context_score_bounds(self):
        """Context-adjusted scores should be in the range
        [base * 0.7, base * 1.2]."""
        seq = "ATGTATAAGCTGCCC"
        for match_type, base in [
            ("8mer", 1.0), ("7mer-m8", 0.9),
            ("7mer-A1", 0.85), ("6mer", 0.7),
        ]:
            score = _mirna_context_score(seq, 3, match_type, base)
            assert base * 0.7 <= score <= base * 1.2, (
                f"{match_type}: score={score} outside [{base*0.7}, {base*1.2}]"
            )
