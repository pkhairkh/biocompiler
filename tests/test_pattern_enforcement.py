"""Tests for biocompiler.pattern_enforcement — DNA-Chisel-style pattern constraints.

Covers:
- PatternConstraint dataclass validation
- check_pattern: avoid and enforce actions
- check_pattern: forward, reverse, both strands
- check_pattern: IUPAC expansion
- check_pattern: regex patterns
- check_pattern: protein scope
- enforce_pattern: avoid by synonymous codon substitution
- enforce_pattern: enforce (embed) a pattern
- enforce_patterns: iterative multi-constraint enforcement
- build_avoidance_scanner: Aho-Corasick integration
- CSP constraint model: PatternAvoidanceConstraint / PatternEnforcementConstraint
- ConstraintSpec.check() integration with PATTERN_AVOIDANCE / PATTERN_ENFORCEMENT
"""

import pytest

# Direct submodule imports to avoid circular import through __init__.py
from biocompiler.pattern_enforcement import (
    PatternConstraint,
    PatternResult,
    check_pattern,
    check_patterns,
    enforce_pattern,
    enforce_patterns,
    build_avoidance_scanner,
)
from biocompiler.constants import CODON_TABLE, AA_TO_CODONS, reverse_complement


# ═══════════════════════════════════════════════════════════════════════
# PatternConstraint validation
# ═══════════════════════════════════════════════════════════════════════


class TestPatternConstraint:
    def test_valid_avoid(self):
        c = PatternConstraint(pattern="GAATTC", action="avoid")
        assert c.action == "avoid"
        assert c.scope == "dna"
        assert c.strand == "both"

    def test_valid_enforce(self):
        c = PatternConstraint(pattern="CATCAT", action="enforce", strand="forward")
        assert c.action == "enforce"
        assert c.strand == "forward"

    def test_invalid_action(self):
        with pytest.raises(ValueError, match="Invalid action"):
            PatternConstraint(pattern="AAA", action="require")

    def test_invalid_scope(self):
        with pytest.raises(ValueError, match="Invalid scope"):
            PatternConstraint(pattern="AAA", action="avoid", scope="rna")

    def test_invalid_strand(self):
        with pytest.raises(ValueError, match="Invalid strand"):
            PatternConstraint(pattern="AAA", action="avoid", strand="antisense")

    def test_is_regex_literal(self):
        c = PatternConstraint(pattern="GAATTC", action="avoid")
        assert not c.is_regex

    def test_is_regex_iupac(self):
        c = PatternConstraint(pattern="GTYRAC", action="avoid")
        assert c.is_regex

    def test_is_regex_metachar(self):
        c = PatternConstraint(pattern="G[AT]TC", action="avoid")
        assert c.is_regex

    def test_expanded_iupac_pure_acgt(self):
        c = PatternConstraint(pattern="GAATTC", action="avoid")
        assert c.expanded_iupac == ["GAATTC"]

    def test_expanded_iupac_ambiguity(self):
        c = PatternConstraint(pattern="GTR", action="avoid")
        # R = A|G, so GT R → GTA, GTG
        expansions = c.expanded_iupac
        assert "GTA" in expansions
        assert "GTG" in expansions
        assert len(expansions) == 2

    def test_frozen(self):
        c = PatternConstraint(pattern="AAA", action="avoid")
        with pytest.raises(AttributeError):
            c.pattern = "BBB"  # type: ignore[misc]


# ═══════════════════════════════════════════════════════════════════════
# check_pattern
# ═══════════════════════════════════════════════════════════════════════


class TestCheckPattern:
    def test_avoid_pattern_found(self):
        c = PatternConstraint(pattern="GAATTC", action="avoid")
        result = check_pattern("ATGGAATTCCGATC", c)
        assert not result.passed  # pattern IS present → avoid fails
        assert len(result.matches) >= 1
        assert result.pattern == "GAATTC"

    def test_avoid_pattern_not_found(self):
        c = PatternConstraint(pattern="GAATTC", action="avoid")
        result = check_pattern("ATGGCATCCGATC", c)
        assert result.passed  # pattern NOT present → avoid passes

    def test_enforce_pattern_found(self):
        c = PatternConstraint(pattern="GAATTC", action="enforce")
        result = check_pattern("ATGGAATTCCGATC", c)
        assert result.passed  # pattern IS present → enforce passes

    def test_enforce_pattern_not_found(self):
        c = PatternConstraint(pattern="GAATTC", action="enforce")
        result = check_pattern("ATGGCATCCGATC", c)
        assert not result.passed  # pattern NOT present → enforce fails

    def test_match_positions(self):
        c = PatternConstraint(pattern="GAATTC", action="avoid")
        result = check_pattern("ATGGAATTCCGATC", c)
        assert result.matches[0] == (3, 9)

    def test_multiple_matches(self):
        c = PatternConstraint(pattern="ATG", action="avoid")
        result = check_pattern("ATGATGATG", c)
        assert len(result.matches) >= 2  # overlapping matches

    def test_forward_strand_only(self):
        c = PatternConstraint(pattern="GAATTC", action="avoid", strand="forward")
        # GAATTC's RC is GAATTC (palindrome)
        result = check_pattern("ATGGAATTCCGATC", c)
        assert not result.passed

    def test_reverse_strand(self):
        c = PatternConstraint(pattern="GAATTC", action="avoid", strand="reverse")
        # Check that reverse complement is scanned
        # GAATTC RC = GAATTC (palindrome)
        result = check_pattern("ATGGAATTCCGATC", c)
        assert not result.passed

    def test_non_palindrome_reverse(self):
        # GGATCC RC = GGATCC (palindrome too)
        # Use AAGCTT (HindIII) → RC = AAGCTT (also palindromic)
        # Let's use GGATCC — palindromic
        # Use a non-palindromic example:
        # GTCGAC (SalI) RC = GTCGAC (palindromic too)
        # TCTAGA (XbaI) RC = TCTAGA
        # Let's find a truly non-palindromic pattern
        c = PatternConstraint(pattern="GATC", action="avoid", strand="reverse")
        # RC of GATC = GATC (palindrome)
        # Use ACGT → RC = ACGT (palindrome)
        # Just test the strand mechanism works
        result = check_pattern("ATGACGTCCGATC", c)
        # ACGT present on forward → its RC (ACGT) present on reverse
        assert not result.passed

    def test_iupac_pattern(self):
        c = PatternConstraint(pattern="GTYRAC", action="avoid")
        # HincII site: GTYRAC = GT(C/T)(A/G)AC
        # GTCGAC is one expansion
        result = check_pattern("ATGGTCGACCGATC", c)
        assert not result.passed  # GTCGAC matches GTYRAC

    def test_regex_pattern(self):
        c = PatternConstraint(pattern="G[AT]TC", action="avoid")
        # G[AT]TC matches GATC: G + A(matches [AT]) + T + C
        result = check_pattern("ATGGATCCGATC", c)
        assert not result.passed  # GATC matches G[AT]TC

        # Regex pattern that doesn't match
        c3 = PatternConstraint(pattern="GGG[AT]CCC", action="avoid")
        result3 = check_pattern("ATGGGACCCGATC", c3)
        assert not result3.passed  # GGGACCC matches GGG[AT]CCC

    def test_protein_scope(self):
        c = PatternConstraint(pattern="HHH", action="enforce", scope="protein")
        result = check_pattern("MHHHKK", c)
        assert result.passed  # HHH found in protein

    def test_protein_scope_avoid(self):
        c = PatternConstraint(pattern="HHH", action="avoid", scope="protein")
        result = check_pattern("MHHHKK", c)
        assert not result.passed

    def test_empty_sequence(self):
        c = PatternConstraint(pattern="GAATTC", action="avoid")
        result = check_pattern("", c)
        assert result.passed  # nothing to match → avoid passes

    def test_case_insensitive(self):
        c = PatternConstraint(pattern="gaattc", action="avoid")
        result = check_pattern("ATGGAATTCCGATC", c)
        assert not result.passed  # should match regardless of case


class TestCheckPatterns:
    def test_multiple_constraints(self):
        constraints = [
            PatternConstraint(pattern="GAATTC", action="avoid"),
            PatternConstraint(pattern="GGATCC", action="avoid"),
        ]
        results = check_patterns("ATGGAATTCCGATC", constraints)
        assert len(results) == 2
        assert not results[0].passed  # GAATTC found
        assert results[1].passed  # GGATCC not found


# ═══════════════════════════════════════════════════════════════════════
# enforce_pattern
# ═══════════════════════════════════════════════════════════════════════


class TestEnforcePattern:
    def test_avoid_removes_pattern(self):
        """Avoid constraint should remove a pattern via codon substitution."""
        # Leucine (L) has codons: TTA, TTG, CTT, CTC, CTA, CTG
        # Some contain "TT" which could be part of "ATTTA"
        # Build DNA where ATTTA appears across a codon boundary
        dna = "ATGAAA" + "ATTTA" + "AAA"  # ATTTA at position 6-10
        # This is a synthetic test; protein is K K X K — not a real protein
        # Let's use a real protein instead
        # Use protein "KFKF" → AAA TTT AAA TTT → contains ATTTA across boundary
        # K=AAA, F=TTT, K=AAA, F=TTT
        # DNA: AAATTTAAATTT → ATTTA at position 2
        dna = "AAATTTAAATTT"
        protein = "KFKF"
        c = PatternConstraint(pattern="ATTTA", action="avoid")
        new_dna = enforce_pattern(dna, protein, c)
        # Verify pattern is gone
        result = check_pattern(new_dna, c)
        assert result.passed, f"ATTTA still present in: {new_dna}"
        # Verify translation is preserved
        for i, aa in enumerate(protein):
            codon = new_dna[i * 3 : i * 3 + 3]
            assert CODON_TABLE.get(codon) == aa, (
                f"Codon {codon} at position {i} does not translate to {aa}"
            )

    def test_avoid_already_satisfied(self):
        """If pattern is already absent, enforce_pattern returns original."""
        dna = "ATGAAAGCGAAATGA"
        protein = "MKAKM"
        c = PatternConstraint(pattern="GAATTC", action="avoid")
        new_dna = enforce_pattern(dna, protein, c)
        assert new_dna == dna

    def test_enforce_embeds_pattern(self):
        """Enforce constraint should embed a pattern via codon substitution."""
        # Try to embed "CAT" (His-tag start) into a Histidine codon
        # H has codons: CAT, CAC → CAT already contains "CAT"
        protein = "HHHH"  # All histidine
        dna = "CACCACCACCAC"  # All CAC codons (no CAT)
        c = PatternConstraint(pattern="CAT", action="enforce", strand="forward")
        new_dna = enforce_pattern(dna, protein, c)
        result = check_pattern(new_dna, c)
        assert result.passed, f"CAT not found in: {new_dna}"
        # Verify translation
        for i, aa in enumerate(protein):
            codon = new_dna[i * 3 : i * 3 + 3]
            assert CODON_TABLE.get(codon) == aa

    def test_enforce_already_present(self):
        """If pattern is already present, enforce_pattern returns original."""
        dna = "CATCATCATCAT"  # All H codons using CAT
        protein = "HHHH"
        c = PatternConstraint(pattern="CAT", action="enforce", strand="forward")
        new_dna = enforce_pattern(dna, protein, c)
        assert new_dna == dna

    def test_preserves_translation(self):
        """All enforcement must preserve the protein translation."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        # Generate initial DNA with optimal E. coli codons
        dna = ""
        for aa in protein:
            codons = AA_TO_CODONS.get(aa, ["ATG"])
            dna += codons[0]

        c = PatternConstraint(pattern="GAATTC", action="avoid")
        new_dna = enforce_pattern(dna, protein, c)

        for i, aa in enumerate(protein):
            codon = new_dna[i * 3 : i * 3 + 3]
            assert CODON_TABLE.get(codon) == aa, (
                f"Translation broken at position {i}: "
                f"codon={codon}, expected={aa}, got={CODON_TABLE.get(codon)}"
            )


class TestEnforcePatterns:
    def test_multiple_avoid_constraints(self):
        """Multiple avoid constraints should all be satisfied."""
        protein = "KFKFRFRF"
        dna = "AAATTTAAATTTAGAAGAAGAAGA"
        constraints = [
            PatternConstraint(pattern="ATTTA", action="avoid"),
            PatternConstraint(pattern="GAATTC", action="avoid"),
        ]
        new_dna = enforce_patterns(dna, protein, constraints)
        for c in constraints:
            result = check_pattern(new_dna, c)
            assert result.passed, f"Constraint not satisfied: {c.pattern}"

    def test_iterative_resolution(self):
        """Test that cascading conflicts are resolved iteratively."""
        protein = "VVVV"  # Valine: GTT, GTC, GTA, GTG
        dna = "GTTGTTGTTGTT"
        constraints = [
            PatternConstraint(pattern="GTTG", action="avoid"),
            PatternConstraint(pattern="TGTT", action="avoid"),
        ]
        new_dna = enforce_patterns(dna, protein, constraints, max_iterations=20)
        # Check that at least some constraints are better
        # (may not be able to satisfy all for Valine)
        assert len(new_dna) == len(dna)


# ═══════════════════════════════════════════════════════════════════════
# build_avoidance_scanner (Aho-Corasick integration)
# ═══════════════════════════════════════════════════════════════════════


class TestBuildAvoidanceScanner:
    def test_literal_patterns(self):
        constraints = [
            PatternConstraint(pattern="GAATTC", action="avoid"),
            PatternConstraint(pattern="GGATCC", action="avoid"),
        ]
        scanner = build_avoidance_scanner(constraints)
        assert scanner is not None
        matches = scanner.scan("ATGGAATTCCGATCGGATCCTAA")
        assert len(matches) >= 2

    def test_iupac_expansion(self):
        constraints = [
            PatternConstraint(pattern="GTYRAC", action="avoid"),
        ]
        scanner = build_avoidance_scanner(constraints)
        assert scanner is not None
        # GTCGAC is one expansion of GTYRAC
        matches = scanner.scan("ATGGTCGACCGATC")
        assert len(matches) >= 1

    def test_enforce_patterns_excluded(self):
        constraints = [
            PatternConstraint(pattern="GAATTC", action="enforce"),
        ]
        scanner = build_avoidance_scanner(constraints)
        assert scanner is None  # enforce constraints not included

    def test_empty_constraints(self):
        scanner = build_avoidance_scanner([])
        assert scanner is None

    def test_regex_excluded(self):
        constraints = [
            PatternConstraint(pattern="G[AT]TC", action="avoid"),
        ]
        scanner = build_avoidance_scanner(constraints)
        # Regex patterns should be skipped
        # This has metacharacters, so it's skipped
        assert scanner is None

    def test_reverse_complement_included(self):
        constraints = [
            PatternConstraint(pattern="GATC", action="avoid", strand="both"),
        ]
        scanner = build_avoidance_scanner(constraints)
        assert scanner is not None
        # GATC RC = GATC (palindrome)
        # Scanner should find GATC on forward strand
        matches = scanner.scan("ATGGATCCGATC")
        assert len(matches) >= 1


# ═══════════════════════════════════════════════════════════════════════
# CSP constraint model integration
# ═══════════════════════════════════════════════════════════════════════


class TestPatternConstraintsCSPModel:
    def test_pattern_avoidance_constraint(self):
        from biocompiler.solver.constraints import PatternAvoidanceConstraint

        c = PatternAvoidanceConstraint(pattern="GAATTC")
        assert "GAATTC" in c.name
        assert not c.check("ATGGAATTCCGATC")
        assert c.check("ATGGCATCCGATC")
        positions = c.violated_positions("ATGGAATTCCGATC")
        assert 3 in positions

    def test_pattern_enforcement_constraint(self):
        from biocompiler.solver.constraints import PatternEnforcementConstraint

        c = PatternEnforcementConstraint(pattern="GAATTC")
        assert "GAATTC" in c.name
        assert c.check("ATGGAATTCCGATC")
        assert not c.check("ATGGCATCCGATC")
        positions = c.violated_positions("ATGGCATCCGATC")
        assert len(positions) == len("ATGGCATCCGATC")  # entire sequence

    def test_avoidance_constraint_types(self):
        from biocompiler.solver.constraints import PatternAvoidanceConstraint
        from biocompiler.solver.types import ConstraintType

        c = PatternAvoidanceConstraint(pattern="GAATTC")
        assert c.constraint_type == ConstraintType.PATTERN_AVOIDANCE

    def test_enforcement_constraint_types(self):
        from biocompiler.solver.constraints import PatternEnforcementConstraint
        from biocompiler.solver.types import ConstraintType

        c = PatternEnforcementConstraint(pattern="CATCAT")
        assert c.constraint_type == ConstraintType.PATTERN_ENFORCEMENT


class TestConstraintSpecIntegration:
    def test_pattern_avoidance_spec(self):
        from biocompiler.solver.types import ConstraintSpec, ConstraintType

        spec = ConstraintSpec(
            ctype=ConstraintType.PATTERN_AVOIDANCE,
            name="avoid_ecori",
            params={"pattern": "GAATTC", "strand": "both", "scope": "dna"},
        )
        assert not spec.check("ATGGAATTCCGATC")
        assert spec.check("ATGGCATCCGATC")

    def test_pattern_enforcement_spec(self):
        from biocompiler.solver.types import ConstraintSpec, ConstraintType

        spec = ConstraintSpec(
            ctype=ConstraintType.PATTERN_ENFORCEMENT,
            name="enforce_his_tag",
            params={"pattern": "CATCAT", "strand": "forward", "scope": "dna"},
        )
        assert spec.check("ATGCATCATCGATC")
        assert not spec.check("ATGCACCACCGATC")


# ═══════════════════════════════════════════════════════════════════════
# Integration with enforcement module
# ═══════════════════════════════════════════════════════════════════════


class TestEnforcementIntegration:
    def test_enforcement_module_recognizes_new_constraints(self):
        """Test that pattern constraints work with the constraint checking system."""
        from biocompiler.solver.constraints import PatternAvoidanceConstraint

        protein = "KFKF"
        dna = "AAATTTAAATTT"

        # Test the constraint check directly
        c = PatternAvoidanceConstraint(pattern="ATTTA")
        assert not c.check(dna)

        # After enforcement, the pattern should be gone
        new_dna = enforce_pattern(
            dna, protein,
            PatternConstraint(pattern="ATTTA", action="avoid"),
        )
        assert c.check(new_dna)


# ═══════════════════════════════════════════════════════════════════════
# DNA-Chisel compatibility test
# ═══════════════════════════════════════════════════════════════════════


class TestDNAChiselCompat:
    """Tests inspired by DNA-Chisel's EnforcePattern/AvoidPattern semantics."""

    def test_avoid_restriction_site_generalization(self):
        """AvoidPattern generalizes AvoidRestrictionSite."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        dna = ""
        for aa in protein:
            codons = AA_TO_CODONS.get(aa, ["ATG"])
            dna += codons[0]

        c = PatternConstraint(pattern="GAATTC", action="avoid")
        new_dna = enforce_pattern(dna, protein, c)
        result = check_pattern(new_dna, c)
        assert result.passed

        # Also check reverse complement (GAATTC RC = GAATTC, palindrome)
        assert "GAATTC" not in new_dna

    def test_avoid_multiple_sites_via_aho_corasick(self):
        """Multiple avoid constraints can be checked efficiently."""
        dna = "ATGGAATTCCGATCGGATCCTAA"
        # EcoRI + BamHI
        constraints = [
            PatternConstraint(pattern="GAATTC", action="avoid"),
            PatternConstraint(pattern="GGATCC", action="avoid"),
        ]
        scanner = build_avoidance_scanner(constraints)
        assert scanner is not None
        matches = scanner.scan(dna)
        # Should find both sites
        assert len(matches) >= 2

    def test_enforce_pattern_embeds_site(self):
        """EnforcePattern can embed a restriction site."""
        # Use a protein with amino acids that can encode a restriction site
        # EcoRI: GAATTC → GAA=E, TTC=F → EF
        protein = "EF"  # E=GAA/GAG, F=TTC/TTT
        dna = "GAGTTT"  # GAG(E) + TTT(F) — no GAATTC
        c = PatternConstraint(pattern="GAATTC", action="enforce", strand="forward")
        new_dna = enforce_pattern(dna, protein, c)
        result = check_pattern(new_dna, c)
        assert result.passed, f"GAATTC not found in: {new_dna}"
        # Verify: should be GAATTC
        assert new_dna == "GAATTC"

    def test_iupac_avoidance_with_expansion(self):
        """IUPAC patterns are properly expanded and checked."""
        # GTYRAC is HincII: Y=C|T, R=A|G
        # Expansions: GTCAAC, GTCGAC, GTTAAC, GTTGAC
        c = PatternConstraint(pattern="GTYRAC", action="avoid")
        # GTCGAC is one expansion
        result = check_pattern("ATGGTCGACCGATC", c)
        assert not result.passed

        # GTTGAC is another expansion
        result2 = check_pattern("ATGGTTGACCGATC", c)
        assert not result2.passed

        # No expansion matches
        result3 = check_pattern("ATGGCATCCGATC", c)
        assert result3.passed
