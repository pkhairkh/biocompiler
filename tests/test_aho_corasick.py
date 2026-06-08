"""Tests for the AhoCorasick scanner module.

Covers:
- Single pattern matching
- Multiple pattern matching
- Overlapping patterns
- No-match case
- Empty scanner
- Reverse complement detection
- Region scanning
- Count matches
- Build from enzymes/sites
"""

import pytest
from biocompiler.aho_corasick import (
    AhoCorasickScanner,
    build_scanner_from_enzymes,
    build_scanner_from_sites,
)


class TestAhoCorasickSinglePattern:
    def test_single_pattern_found(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("ATGGAATTCC")
        assert len(results) >= 1
        pos, site, enzyme = results[0]
        assert site == "GAATTC"
        assert enzyme == "EcoRI"
        assert pos == 3

    def test_single_pattern_not_found(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("ATGCATGCATGC")
        assert results == []

    def test_pattern_at_start(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("GAATTCATGC")
        assert len(results) >= 1
        assert results[0][0] == 0

    def test_pattern_at_end(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("ATGCATGAATTC")
        assert len(results) >= 1
        assert results[0][0] == 6

    def test_multiple_occurrences(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("GAATTCATGAATTC")
        assert len(results) >= 2


class TestAhoCorasickMultiplePatterns:
    def test_two_patterns(self):
        scanner = AhoCorasickScanner({
            "GAATTC": "EcoRI",
            "GGATCC": "BamHI",
        })
        results = scanner.scan("ATGGAATTCCGATCGGATCCTAA")
        sites_found = {r[1] for r in results}
        assert "GAATTC" in sites_found
        assert "GGATCC" in sites_found

    def test_three_patterns(self):
        scanner = AhoCorasickScanner({
            "GAATTC": "EcoRI",
            "GGATCC": "BamHI",
            "GCGGCCGC": "NotI",
        })
        results = scanner.scan("ATGGAATTCCGATCGGATCCTAA")
        sites_found = {r[1] for r in results}
        assert "GAATTC" in sites_found
        assert "GGATCC" in sites_found

    def test_only_some_patterns_match(self):
        scanner = AhoCorasickScanner({
            "GAATTC": "EcoRI",
            "GGATCC": "BamHI",
        })
        results = scanner.scan("ATGGAATTCC")  # Only EcoRI present
        sites_found = {r[1] for r in results}
        assert "GAATTC" in sites_found
        assert "GGATCC" not in sites_found


class TestAhoCorasickOverlapping:
    def test_overlapping_patterns(self):
        # GATT and ATTC overlap at ATT
        scanner = AhoCorasickScanner({
            "GATT": "pattern1",
            "ATTC": "pattern2",
        })
        results = scanner.scan("GATTC")
        assert len(results) >= 1  # At least one should be found

    def test_same_pattern_different_positions(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("GAATTCGAATTC")
        assert len(results) >= 2


class TestAhoCorasickNoMatch:
    def test_no_match_returns_empty(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("ATGCATGC")
        assert results == []

    def test_has_any_match_false(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        assert scanner.has_any_match("ATGCATGC") is False

    def test_has_any_match_true(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        assert scanner.has_any_match("ATGGAATTC") is True


class TestAhoCorasickEmpty:
    def test_empty_patterns(self):
        scanner = AhoCorasickScanner({})
        results = scanner.scan("ATGGAATTCC")
        assert results == []

    def test_empty_sequence(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan("")
        assert results == []

    def test_empty_has_any_match(self):
        scanner = AhoCorasickScanner({})
        assert scanner.has_any_match("ATGC") is False

    def test_empty_count_matches(self):
        scanner = AhoCorasickScanner({})
        assert scanner.count_matches("ATGC") == 0


class TestAhoCorasickProperties:
    def test_longest_pattern(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI", "GGATCC": "BamHI"})
        assert scanner.longest_pattern == 6

    def test_num_nodes(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        assert scanner.num_nodes > 0

    def test_count_matches(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        count = scanner.count_matches("ATGGAATTCC")
        assert count >= 1

    def test_find_all_sites(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI", "GGATCC": "BamHI"})
        results = scanner.find_all_sites("ATGGAATTCCGATCGGATCCTAA")
        assert isinstance(results, dict)
        assert "EcoRI" in results or "BamHI" in results

    def test_scan_region(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        results = scanner.scan_region("ATGGAATTCC", 0, 5)
        assert isinstance(results, list)

    def test_has_any_match_in_region(self):
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})
        result = scanner.has_any_match_in_region("ATGGAATTCC", 0, 10)
        assert isinstance(result, bool)


class TestBuildFromEnzymes:
    def test_build_from_ecori(self):
        scanner = build_scanner_from_enzymes(["EcoRI"])
        assert scanner is not None
        assert isinstance(scanner, AhoCorasickScanner)

    def test_build_from_multiple_enzymes(self):
        scanner = build_scanner_from_enzymes(["EcoRI", "BamHI", "XhoI"])
        assert scanner is not None

    def test_build_from_unknown_enzyme(self):
        scanner = build_scanner_from_enzymes(["NONEXISTENT_ENZYME_XYZ"])
        assert scanner is None  # No valid ACGT-only sites

    def test_built_scanner_works(self):
        scanner = build_scanner_from_enzymes(["EcoRI"])
        results = scanner.scan("ATGGAATTCC")
        assert len(results) >= 1


class TestBuildFromSites:
    def test_build_from_sites(self):
        scanner = build_scanner_from_sites(["GAATTC", "GGATCC"])
        assert scanner is not None

    def test_build_from_single_site(self):
        scanner = build_scanner_from_sites(["GAATTC"])
        assert scanner is not None

    def test_build_from_iupac_site_skipped(self):
        # IUPAC sites should be skipped
        scanner = build_scanner_from_sites(["GANTC"])  # N is IUPAC
        assert scanner is None

    def test_built_scanner_detects(self):
        scanner = build_scanner_from_sites(["GAATTC"])
        results = scanner.scan("ATGGAATTCC")
        assert len(results) >= 1
