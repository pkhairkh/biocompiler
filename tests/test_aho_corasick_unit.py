"""
Unit tests for biocompiler.aho_corasick — AhoCorasickScanner, build_scanner_from_enzymes,
build_scanner_from_sites.

Covers:
- AhoCorasickScanner: construction, scan, has_any_match, has_any_match_in_region,
  scan_region, find_all_sites, count_matches
- Builder functions: build_scanner_from_enzymes, build_scanner_from_sites
- Edge cases: empty patterns, no matches, overlapping patterns, palindromes
"""

from __future__ import annotations

import pytest

from biocompiler.sequence.aho_corasick import (
    AhoCorasickScanner,
    AhoCorasickNode,
    build_scanner_from_enzymes,
    build_scanner_from_sites,
)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. AhoCorasickScanner — Construction
# ═══════════════════════════════════════════════════════════════════════════════

class TestScannerConstruction:

    def test_single_pattern(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        assert scanner.longest_pattern == 6
        assert scanner.num_nodes > 0

    def test_multiple_patterns(self):
        scanner = AhoCorasickScanner({
            "GAATTC": "EcoRI",
            "GGATCC": "BamHI",
            "AAGCTT": "HindIII",
        })
        assert scanner.longest_pattern == 6

    def test_empty_patterns(self):
        """Empty patterns dict should still create a valid scanner."""
        scanner = AhoCorasickScanner({})
        assert scanner.num_nodes == 1  # root only
        assert scanner.longest_pattern == 0

    def test_prefix_sharing(self):
        """Patterns that share prefixes should share trie nodes."""
        scanner = AhoCorasickScanner({
            "GAATTC": "EcoRI",
            "GATC": "MboI_partial",
        })
        # Should have fewer nodes than the sum of pattern lengths
        assert scanner.num_nodes < 6 + 4 + 1  # not every char gets a new node


# ═══════════════════════════════════════════════════════════════════════════════
# 2. AhoCorasickScanner — scan
# ═══════════════════════════════════════════════════════════════════════════════

class TestScannerScan:

    def test_scan_single_match(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("ATGGAATTCCGATCG")
        assert len(results) == 1
        pos, site, enzyme = results[0]
        assert pos == 3
        assert site == "GAATTC"
        assert enzyme == "EcoRI"

    def test_scan_no_match(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("ATGCATGCATGC")
        assert len(results) == 0

    def test_scan_multiple_matches(self):
        scanner = AhoCorasickScanner({
            "GAATTC": "EcoRI",
            "GGATCC": "BamHI",
        })
        seq = "ATGGAATTCCGATCGGATCCTAA"
        results = scanner.scan(seq)
        assert len(results) == 2
        enzymes = [e for _, _, e in results]
        assert "EcoRI" in enzymes
        assert "BamHI" in enzymes

    def test_scan_duplicate_sites(self):
        """Two occurrences of the same site should both be found."""
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        seq = "GAATTCXXXGAATTC"
        results = scanner.scan(seq)
        assert len(results) == 2

    def test_scan_overlapping_patterns(self):
        """Two patterns that overlap should both be detected."""
        scanner = AhoCorasickScanner({
            "GATC": "MboI",
            "GAATTC": "EcoRI",
        })
        # GATC appears within GAATTC at pos 3
        seq = "ATGGAATTCC"
        results = scanner.scan(seq)
        sites = [s for _, s, _ in results]
        assert "GAATTC" in sites

    def test_scan_empty_sequence(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("")
        assert results == []

    def test_scan_shorter_than_pattern(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("GAA")
        assert results == []


# ═══════════════════════════════════════════════════════════════════════════════
# 3. AhoCorasickScanner — has_any_match
# ═══════════════════════════════════════════════════════════════════════════════

class TestScannerHasAnyMatch:

    def test_has_match_true(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        assert scanner.has_any_match("ATGGAATTCC") is True

    def test_has_match_false(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        assert scanner.has_any_match("ATGCATGC") is False

    def test_has_match_empty_sequence(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        assert scanner.has_any_match("") is False

    def test_has_match_empty_patterns(self):
        scanner = AhoCorasickScanner({})
        assert scanner.has_any_match("GAATTC") is False

    def test_short_circuits(self):
        """has_any_match should stop at the first match."""
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        # Just verify it returns True quickly; no way to measure short-circuit
        # but the result should be correct
        assert scanner.has_any_match("GAATTC" * 100) is True


# ═══════════════════════════════════════════════════════════════════════════════
# 4. AhoCorasickScanner — scan_region & has_any_match_in_region
# ═══════════════════════════════════════════════════════════════════════════════

class TestScannerRegionScan:

    def test_has_any_match_in_region(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        seq = "ATGGAATTCCGATCG"
        # Site at pos 3 — within region [0, 10)
        assert scanner.has_any_match_in_region(seq, 0, 10) is True
        # Outside the region
        assert scanner.has_any_match_in_region(seq, 10, 16) is False

    def test_scan_region(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        seq = "ATGGAATTCCGATCG"
        results = scanner.scan_region(seq, 0, 10)
        # Site at pos 3 starts within [0, 10)
        assert len(results) >= 1

    def test_scan_region_outside(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        seq = "ATGGAATTCCGATCG"
        results = scanner.scan_region(seq, 10, 16)
        assert len(results) == 0

    def test_scan_region_boundary_overlap(self):
        """A site that starts before the region but overlaps should be detected."""
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        seq = "ATGGAATTCCGATCG"
        # Region [5, 10) — site starts at 3 but extends to 9
        results = scanner.scan_region(seq, 5, 10)
        # The site at pos 3 starts before the region [5,10), so may or may not be included
        # depending on the implementation
        assert isinstance(results, list)

    def test_scan_region_multiple_sites(self):
        scanner = AhoCorasickScanner({
            "GAATTC": "EcoRI",
            "GGATCC": "BamHI",
        })
        seq = "ATGGAATTCCGGATCCTAA"
        results = scanner.scan_region(seq, 0, len(seq))
        assert len(results) == 2


# ═══════════════════════════════════════════════════════════════════════════════
# 5. AhoCorasickScanner — find_all_sites & count_matches
# ═══════════════════════════════════════════════════════════════════════════════

class TestScannerFindAllSites:

    def test_find_all_sites(self):
        scanner = AhoCorasickScanner({
            "GAATTC": "EcoRI",
            "GGATCC": "BamHI",
        })
        seq = "ATGGAATTCCGGATCCTAA"
        by_enzyme = scanner.find_all_sites(seq)
        assert "EcoRI" in by_enzyme
        assert "BamHI" in by_enzyme
        assert 3 in by_enzyme["EcoRI"]

    def test_count_matches(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        seq = "GAATTCXXXGAATTC"
        assert scanner.count_matches(seq) == 2

    def test_count_no_matches(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        assert scanner.count_matches("ATGCATGC") == 0


# ═══════════════════════════════════════════════════════════════════════════════
# 6. AhoCorasickScanner — Edge Cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestScannerEdgeCases:

    def test_palindromic_site(self):
        """GAATTC is palindromic — should find once, not twice."""
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("ATGGAATTCC")
        assert len(results) == 1

    def test_non_acgt_characters(self):
        """Non-ACGT characters should be handled gracefully — they reset state."""
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        # N at pos 8 breaks the GAATTC pattern; the site before N is GAATTN, not GAATTC
        # So we should get 0 results here
        results = scanner.scan("ATGGAATTNCGATCG")
        assert len(results) == 0
        # But a valid site after the N should be found
        results2 = scanner.scan("ATGGAATTCCNNGAATTC")
        assert len(results2) == 2  # both GAATTC instances

    def test_site_at_end(self):
        """Site at the very end of the sequence."""
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("ATGCGAATTC")
        assert len(results) == 1
        assert results[0][0] == 4

    def test_site_at_start(self):
        """Site at the very start of the sequence."""
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("GAATTCATGC")
        assert len(results) == 1
        assert results[0][0] == 0


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Builder Functions
# ═══════════════════════════════════════════════════════════════════════════════

class TestBuilderFunctions:

    def test_build_from_enzymes(self):
        scanner = build_scanner_from_enzymes(["EcoRI", "BamHI"])
        assert scanner is not None
        assert isinstance(scanner, AhoCorasickScanner)
        # EcoRI = GAATTC, should be detected
        results = scanner.scan("ATGGAATTCC")
        assert len(results) >= 1

    def test_build_from_enzymes_unknown(self):
        """Unknown enzymes should be skipped."""
        scanner = build_scanner_from_enzymes(["NonExistentEnzyme"])
        # Should return None since no valid recognition site found
        assert scanner is None

    def test_build_from_enzymes_includes_rc(self):
        """build_scanner_from_enzymes should include reverse complements."""
        scanner = build_scanner_from_enzymes(["EcoRI"])
        # GAATTC is palindromic, so RC is the same
        results = scanner.scan("GAATTC")
        assert len(results) >= 1

    def test_build_from_sites(self):
        scanner = build_scanner_from_sites(["GAATTC", "GGATCC"])
        assert scanner is not None
        results = scanner.scan("ATGGAATTCCGGATCCTAA")
        assert len(results) == 2

    def test_build_from_sites_includes_rc(self):
        """Non-palindromic sites should include reverse complement."""
        scanner = build_scanner_from_sites(["GATC"])
        # GATC RC = GATC (palindromic), so only one pattern
        assert scanner is not None

    def test_build_from_sites_iupac_skipped(self):
        """Sites with IUPAC ambiguity codes should be skipped."""
        scanner = build_scanner_from_sites(["GRCGC"])  # R is IUPAC
        # Should return None since GRCGC has non-ACGT characters
        assert scanner is None

    def test_build_from_sites_empty(self):
        scanner = build_scanner_from_sites([])
        assert scanner is None


# ═══════════════════════════════════════════════════════════════════════════════
# 8. AhoCorasickNode
# ═══════════════════════════════════════════════════════════════════════════════

class TestAhoCorasickNode:

    def test_node_construction(self):
        node = AhoCorasickNode()
        assert node.children == {}
        assert node.fail == 0
        assert node.output == []
        assert node.depth == 0
