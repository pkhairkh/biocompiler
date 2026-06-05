"""Aho-Corasick multi-pattern string matching for restriction site detection.

Scans a DNA sequence for all restriction enzyme recognition sites simultaneously,
including reverse complements. Complexity: O(L + M) where L = sequence length,
M = total number of matches.

This replaces the per-enzyme O(N * L * site_len) scanning approach with a single
O(L + M) pass, where N = number of enzymes, L = sequence length.

The automaton is built once (in HybridOptimizer.__init__ or on first use) and
reused for every scan, making repeated constraint checks extremely fast.

Performance note:
    The scanner uses a precomputed flat-array transition table (delta) for
    O(1) state transitions per character. This avoids the dict-lookup overhead
    of a typical Python Aho-Corasick implementation and makes the scan loop
    as fast as possible in pure Python.

Usage:
    from biocompiler.aho_corasick import AhoCorasickScanner

    # Build scanner with all enzyme sites + reverse complements
    patterns = {
        "GAATTC": "EcoRI",
        "GGATCC": "BamHI",
    }
    scanner = AhoCorasickScanner(patterns)

    # Scan full sequence
    matches = scanner.scan("ATGGAATTCCGATCGGATCCTAA")
    # [(3, "GAATTC", "EcoRI"), (15, "GGATCC", "BamHI")]

    # Check if ANY site is present (boolean, stops at first match)
    has_site = scanner.has_any_match("ATGGAATTCCGATCGGATCCTAA")
    # True

    # Check a local region only (for incremental constraint checking)
    local_matches = scanner.scan_region("ATGGAATTCCGATCGGATCCTAA", 0, 10)
    # [(3, "GAATTC", "EcoRI")]
"""

from __future__ import annotations

from collections import deque
from typing import Optional

# ── NUMBA integration ──────────────────────────────────────────────
try:
    from .numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        scan_restriction_sites as _numba_scan_rs,
        seq_to_bytes as _seq_to_bytes,
    )
except ImportError:
    _HAS_NUMBA = False

HAS_NUMBA: bool = _HAS_NUMBA

__all__ = [
    "AhoCorasickNode",
    "AhoCorasickScanner",
    "build_scanner_from_enzymes",
    "build_scanner_from_sites",
]

# DNA alphabet encoding: A=0, C=1, G=2, T=3
_CHAR_TO_IDX: dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}
_ALPHABET_SIZE = 4
_IDX_TO_CHAR = "ACGT"


class AhoCorasickNode:
    """Node in the Aho-Corasick trie/automaton (used during construction only).

    Uses __slots__ for memory efficiency since we may create many nodes
    (one per distinct prefix across all patterns).

    After construction, the transition table is flattened into a list
    for faster scanning.
    """

    __slots__ = ("children", "fail", "output", "depth")

    def __init__(self) -> None:
        self.children: dict[str, int] = {}
        self.fail: int = 0
        self.output: list[tuple[str, str]] = []  # (site_string, enzyme_name)
        self.depth: int = 0


class AhoCorasickScanner:
    """Multi-pattern DNA scanner using Aho-Corasick algorithm.

    Builds a finite automaton from a set of patterns (enzyme recognition sites)
    that can scan any DNA sequence in O(L + M) time where L = sequence length
    and M = total number of pattern matches.

    The automaton is built once and reused for all scans. The transition table
    is stored as a flat list for fast integer-indexed access during scanning.

    Complexity:
        Build: O(total_pattern_length * alphabet_size) for trie + failure links
        Scan: O(L + M) where L = len(sequence), M = total matches

    Args:
        patterns: Mapping of site_string -> enzyme_name. All site strings
            must be uppercase DNA (A, C, G, T only). Reverse complements
            should be included as separate entries if both strands need
            to be detected.

    Examples:
        >>> scanner = AhoCorasickScanner({"GAATTC": "EcoRI", "GGATCC": "BamHI"})
        >>> scanner.scan("ATGGAATTCCGATCGGATCCTAA")
        [(3, 'GAATTC', 'EcoRI'), (15, 'GGATCC', 'BamHI')]

        >>> scanner.has_any_match("ATGCATGCATGC")
        False

        >>> scanner.has_any_match_in_region("ATGGAATTCC", 0, 10)
        True
    """

    def __init__(self, patterns: dict[str, str]) -> None:
        self._patterns: dict[str, str] = dict(patterns)  # defensive copy
        self._longest_pattern: int = 0
        self._num_nodes: int = 0

        # Build the automaton using node-based construction
        nodes: list[AhoCorasickNode] = []
        self._build_trie(nodes)
        self._build_failure_links(nodes)

        # Flatten the transition table into a list for fast scanning.
        # delta[state * ALPHABET_SIZE + char_idx] = next_state
        self._delta: list[int] = [0] * (self._num_nodes * _ALPHABET_SIZE)
        # output_table[state] = list of (site_string, enzyme_name) or None
        self._output_table: list[list[tuple[str, str]] | None] = [None] * self._num_nodes

        for state in range(self._num_nodes):
            for ch, child_idx in nodes[state].children.items():
                ci = _CHAR_TO_IDX[ch]
                self._delta[state * _ALPHABET_SIZE + ci] = child_idx
            # Fill in remaining transitions from failure links (already computed
            # in _build_failure_links which sets up the children dict completely)
            # Actually, we need to compute the full delta including failure transitions
            pass

        # Recompute delta with full failure-link transitions (the nodes' children
        # only contain trie edges, not the full transition function)
        self._compute_full_delta(nodes)

        # Build output table
        for state in range(self._num_nodes):
            if nodes[state].output:
                self._output_table[state] = nodes[state].output

    def _build_trie(self, nodes: list[AhoCorasickNode]) -> None:
        """Phase 1: Build the trie from all patterns."""
        root = AhoCorasickNode()
        root.depth = 0
        nodes.append(root)

        for site, enzyme in self._patterns.items():
            current = 0  # start at root
            for ch in site:
                if ch not in nodes[current].children:
                    new_idx = len(nodes)
                    new_node = AhoCorasickNode()
                    new_node.depth = nodes[current].depth + 1
                    nodes.append(new_node)
                    nodes[current].children[ch] = new_idx
                current = nodes[current].children[ch]
            # Record that this pattern ends at the current node
            nodes[current].output.append((site, enzyme))
            if len(site) > self._longest_pattern:
                self._longest_pattern = len(site)

        self._num_nodes = len(nodes)

    def _build_failure_links(self, nodes: list[AhoCorasickNode]) -> None:
        """Phase 2: Compute failure links using BFS."""
        queue: deque[int] = deque()

        # Initialize: children of root have failure link = root
        for ch, child_idx in nodes[0].children.items():
            nodes[child_idx].fail = 0
            queue.append(child_idx)

        while queue:
            current = queue.popleft()
            current_node = nodes[current]

            for ch, child_idx in current_node.children.items():
                fail = current_node.fail
                while fail != 0 and ch not in nodes[fail].children:
                    fail = nodes[fail].fail

                if ch in nodes[fail].children and nodes[fail].children[ch] != child_idx:
                    nodes[child_idx].fail = nodes[fail].children[ch]
                else:
                    nodes[child_idx].fail = 0

                # Merge output from failure link
                fail_node = nodes[nodes[child_idx].fail]
                if fail_node.output:
                    nodes[child_idx].output = (
                        nodes[child_idx].output + fail_node.output
                    )

                queue.append(child_idx)

    def _compute_full_delta(self, nodes: list[AhoCorasickNode]) -> None:
        """Compute the full transition table including failure-link transitions.

        For each state and each character, delta[state][char] gives the
        next state, following failure links as needed. This is the key
        optimization: instead of following failure links at scan time
        (which requires a while loop), we precompute all transitions.
        """
        # Use BFS to compute the full delta table
        for state in range(self._num_nodes):
            for ci in range(_ALPHABET_SIZE):
                ch = _IDX_TO_CHAR[ci]
                # Follow trie edges first, then failure links
                s = state
                while s != 0 and ch not in nodes[s].children:
                    s = nodes[s].fail
                if ch in nodes[s].children:
                    self._delta[state * _ALPHABET_SIZE + ci] = nodes[s].children[ch]
                else:
                    self._delta[state * _ALPHABET_SIZE + ci] = 0

    @property
    def longest_pattern(self) -> int:
        """Length of the longest pattern in the automaton."""
        return self._longest_pattern

    @property
    def num_nodes(self) -> int:
        """Number of nodes in the automaton."""
        return self._num_nodes

    def scan(self, sequence: str) -> list[tuple[int, str, str]]:
        """Scan sequence for all patterns simultaneously.

        Single-pass O(L + M) scan where L = len(sequence) and M = total
        number of matches. Returns all positions where any pattern matches.

        Args:
            sequence: DNA sequence string (uppercase ACGT).

        Returns:
            List of (position, site_string, enzyme_name) tuples, sorted by
            position. If multiple patterns match at the same position, they
            are all included.
        """
        results: list[tuple[int, str, str]] = []
        state = 0
        delta = self._delta
        output_table = self._output_table
        alpha = _ALPHABET_SIZE
        char_to_idx = _CHAR_TO_IDX.get

        for i in range(len(sequence)):
            ci = char_to_idx(sequence[i])
            if ci is None:
                state = 0
                continue
            state = delta[state * alpha + ci]
            output = output_table[state]
            if output is not None:
                for site, enzyme in output:
                    pos = i - len(site) + 1
                    results.append((pos, site, enzyme))

        return results

    def has_any_match(self, sequence: str) -> bool:
        """Check if ANY pattern matches in the sequence.

        Short-circuits on the first match found. Much faster than scan()
        when you only need a boolean answer.

        Args:
            sequence: DNA sequence string (uppercase ACGT).

        Returns:
            True if any pattern is found, False otherwise.
        """
        state = 0
        delta = self._delta
        output_table = self._output_table
        alpha = _ALPHABET_SIZE
        char_to_idx = _CHAR_TO_IDX.get

        for i in range(len(sequence)):
            ci = char_to_idx(sequence[i])
            if ci is None:
                state = 0
                continue
            state = delta[state * alpha + ci]
            if output_table[state] is not None:
                return True

        return False

    def has_any_match_in_region(self, sequence: str, start: int, end: int) -> bool:
        """Check if any pattern matches within a region of the sequence.

        Extracts the local region and scans it from the initial state,
        avoiding the O(start) fast-forward overhead. The region is extended
        slightly to catch patterns that overlap with the boundaries.

        Args:
            sequence: Full DNA sequence string (uppercase ACGT).
            start: Start position of region (inclusive).
            end: End position of region (exclusive).

        Returns:
            True if any pattern is found in the region, False otherwise.
        """
        # Extend the region to catch patterns that overlap boundaries
        region_start = max(0, start - self._longest_pattern + 1)
        region_end = min(len(sequence), end + self._longest_pattern - 1)
        region = sequence[region_start:region_end]

        # Scan from initial state — no fast-forward needed
        state = 0
        delta = self._delta
        output_table = self._output_table
        alpha = _ALPHABET_SIZE
        char_to_idx = _CHAR_TO_IDX.get

        for i in range(len(region)):
            ci = char_to_idx(region[i])
            if ci is None:
                state = 0
                continue
            state = delta[state * alpha + ci]
            if output_table[state] is not None:
                return True

        return False

    def scan_region(self, sequence: str, start: int, end: int) -> list[tuple[int, str, str]]:
        """Scan a region of the sequence for all patterns.

        Detects patterns that START within [start, end), even if they extend
        past end. This is useful for incremental constraint checking after
        a codon swap — only the affected region needs to be rescanned.

        Args:
            sequence: Full DNA sequence string (uppercase ACGT).
            start: Start position of region (inclusive).
            end: End position of region (exclusive).

        Returns:
            List of (position, site_string, enzyme_name) tuples where
            position is within [start, end).
        """
        results: list[tuple[int, str, str]] = []
        state = 0
        delta = self._delta
        output_table = self._output_table
        alpha = _ALPHABET_SIZE
        char_to_idx = _CHAR_TO_IDX.get
        longest = self._longest_pattern
        seq_len = len(sequence)

        # Need to start scanning from (start - longest_pattern + 1) to catch
        # patterns that overlap with the start of our region
        scan_start = max(0, start - longest + 1)
        scan_end = min(seq_len, end + longest - 1)

        # Fast-forward state to scan_start position
        for i in range(scan_start):
            ci = char_to_idx(sequence[i])
            if ci is None:
                state = 0
                continue
            state = delta[state * alpha + ci]

        # Scan the region of interest
        for i in range(scan_start, scan_end):
            ci = char_to_idx(sequence[i])
            if ci is None:
                state = 0
                continue
            state = delta[state * alpha + ci]
            output = output_table[state]
            if output is not None:
                for site, enzyme in output:
                    pos = i - len(site) + 1
                    if start <= pos < end:
                        results.append((pos, site, enzyme))

        return results

    def find_all_sites(self, sequence: str) -> dict[str, list[int]]:
        """Scan sequence and return matches grouped by enzyme name.

        Convenience method that groups the flat results from scan() into
        a dictionary mapping enzyme_name -> sorted list of positions.

        Args:
            sequence: DNA sequence string (uppercase ACGT).

        Returns:
            Dictionary mapping enzyme_name to sorted list of positions
            where the site was found.
        """
        results: dict[str, list[int]] = {}
        for pos, site, enzyme in self.scan(sequence):
            results.setdefault(enzyme, []).append(pos)
        for enzyme in results:
            results[enzyme].sort()
        return results

    def count_matches(self, sequence: str) -> int:
        """Count total number of pattern matches in sequence.

        Faster than len(scanner.scan(sequence)) because it doesn't
        build the result list.

        Args:
            sequence: DNA sequence string (uppercase ACGT).

        Returns:
            Total number of pattern matches.
        """
        count = 0
        state = 0
        delta = self._delta
        output_table = self._output_table
        alpha = _ALPHABET_SIZE
        char_to_idx = _CHAR_TO_IDX.get

        for i in range(len(sequence)):
            ci = char_to_idx(sequence[i])
            if ci is None:
                state = 0
                continue
            state = delta[state * alpha + ci]
            output = output_table[state]
            if output is not None:
                count += len(output)

        return count


def build_scanner_from_enzymes(
    enzyme_names: list[str],
) -> Optional[AhoCorasickScanner]:
    """Build an AhoCorasickScanner from a list of enzyme names.

    Convenience function that resolves enzyme names to recognition sites,
    includes reverse complements, and filters out IUPAC-only sites.

    Args:
        enzyme_names: List of restriction enzyme names (e.g. ["EcoRI", "BamHI"]).

    Returns:
        AhoCorasickScanner ready to scan, or None if no concrete (ACGT-only)
        sites were found.
    """
    from .restriction_sites import get_recognition_site
    from .constants import reverse_complement

    patterns: dict[str, str] = {}

    for enzyme in enzyme_names:
        site = get_recognition_site(enzyme)
        if site is None:
            continue
        site_upper = site.upper()
        # Skip sites with IUPAC ambiguity codes (they need regex matching)
        if any(b not in "ACGT" for b in site_upper):
            continue
        patterns[site_upper] = enzyme
        rc = reverse_complement(site_upper)
        if rc != site_upper:
            patterns[rc] = enzyme

    if not patterns:
        return None

    return AhoCorasickScanner(patterns)


def build_scanner_from_sites(
    sites: list[str],
) -> Optional[AhoCorasickScanner]:
    """Build an AhoCorasickScanner from a list of site sequences.

    Convenience function for BioOptimizer's _greedy_optimize which works
    with site sequences directly rather than enzyme names. Includes reverse
    complements automatically and filters out IUPAC-only sites.

    Args:
        sites: List of recognition site sequences (e.g. ["GAATTC", "GGATCC"]).

    Returns:
        AhoCorasickScanner ready to scan, or None if no concrete (ACGT-only)
        sites were found.
    """
    from .constants import reverse_complement

    patterns: dict[str, str] = {}

    for site in sites:
        site_upper = site.upper()
        # Skip sites with IUPAC ambiguity codes (they need regex matching)
        if any(b not in "ACGT" for b in site_upper):
            continue
        # Use the site string as both the pattern and the label
        patterns[site_upper] = site_upper
        rc = reverse_complement(site_upper)
        if rc != site_upper:
            patterns[rc] = site_upper

    if not patterns:
        return None

    return AhoCorasickScanner(patterns)
