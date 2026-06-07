"""BioCompiler Arbitrary Pattern Enforcement (EnforcePattern / AvoidPattern).

Generalizes restriction-site avoidance into a DNA-Chisel-style pattern
constraint system.  Two complementary actions are supported:

- **EnforcePattern**: The pattern *must* appear in the sequence.
  Useful for embedding restriction sites, affinity tags, or other
  sequence motifs that the design requires.

- **AvoidPattern**: The pattern *must not* appear in the sequence.
  Generalization of ``AvoidRestrictionSite`` — works with any regex
  or literal pattern, not just enzyme recognition sequences.

Both actions preserve the protein translation (the codon table
constraint is never violated).

Pattern matching uses:
- **regex** for patterns containing IUPAC ambiguity codes or ``[...]``
  character classes;
- **Aho-Corasick** for multi-literal *avoid* constraints (O(L+M)
  scan);
- **simple string search** for single literal patterns.

Usage::

    from biocompiler.pattern_enforcement import (
        PatternConstraint, PatternResult, check_pattern, enforce_pattern,
    )

    # Avoid a restriction site (generalized)
    c = PatternConstraint(pattern="GAATTC", action="avoid", scope="dna", strand="both")
    result = check_pattern("ATGGAATTCCGATC", c)
    assert not result.passed  # GAATTC is present → violation

    # Enforce that a His-tag motif appears
    c = PatternConstraint(pattern="CATCAT", action="enforce", scope="dna", strand="forward")
    dna = enforce_pattern(dna, protein, c, codon_table)
"""

from __future__ import annotations

import re
import logging
from dataclasses import dataclass, field
from typing import Optional

from biocompiler.shared.constants import CODON_TABLE, AA_TO_CODONS, reverse_complement, IUPAC_EXPAND
from .aho_corasick import AhoCorasickScanner

__all__ = [
    "PatternConstraint",
    "PatternResult",
    "check_pattern",
    "check_patterns",
    "enforce_pattern",
    "enforce_patterns",
    "build_avoidance_scanner",
]

logger = logging.getLogger(__name__)


# ─── Data structures ──────────────────────────────────────────────────────────


@dataclass(frozen=True)
class PatternConstraint:
    """A pattern-based constraint on a DNA or protein sequence.

    Attributes:
        pattern: The pattern to match.  Can be a literal DNA string
            (e.g. ``"GAATTC"``) or a regex (e.g. ``"G[AT]TC"``).
            IUPAC ambiguity codes are also supported (e.g. ``"GTYRAC"``).
        action: ``"enforce"`` (pattern MUST appear) or ``"avoid"``
            (pattern MUST NOT appear).
        scope: ``"dna"`` (match against the DNA sequence) or
            ``"protein"`` (match against the translated protein).
        strand: ``"both"`` (check forward + reverse complement),
            ``"forward"`` (forward strand only), or ``"reverse"``
            (reverse complement only).
    """

    pattern: str
    action: str  # "enforce" or "avoid"
    scope: str = "dna"  # "dna" or "protein"
    strand: str = "both"  # "both", "forward", "reverse"

    def __post_init__(self) -> None:
        if self.action not in ("enforce", "avoid"):
            raise ValueError(
                f"Invalid action '{self.action}'. Must be 'enforce' or 'avoid'."
            )
        if self.scope not in ("dna", "protein"):
            raise ValueError(
                f"Invalid scope '{self.scope}'. Must be 'dna' or 'protein'."
            )
        if self.strand not in ("both", "forward", "reverse"):
            raise ValueError(
                f"Invalid strand '{self.strand}'. Must be 'both', 'forward', or 'reverse'."
            )

    @property
    def is_regex(self) -> bool:
        """True if the pattern contains regex metacharacters or IUPAC codes."""
        meta_chars = set(r"\.^$*+?{}[]|()")
        if any(c in meta_chars for c in self.pattern):
            return True
        # IUPAC ambiguity codes (non-ACGT) make it regex-like
        if self.scope == "dna" and any(b not in "ACGTacgt" for b in self.pattern):
            return True
        return False

    @property
    def expanded_iupac(self) -> list[str]:
        """Expand IUPAC ambiguity codes into all concrete ACGT sequences.

        Returns a singleton list [pattern] if the pattern is already pure ACGT.
        Returns an empty list if the pattern contains regex metacharacters
        (use regex matching instead of literal expansion).
        """
        pat_upper = self.pattern.upper()
        # If the pattern contains regex metacharacters, don't expand
        # (use regex matching instead)
        meta_chars = set(r"\\.^$*+?{}[]|()")
        if any(c in meta_chars for c in self.pattern):
            return []
        if all(b in "ACGT" for b in pat_upper):
            return [pat_upper]
        expansions: list[str] = [""]
        for base in pat_upper:
            options = IUPAC_EXPAND.get(base)
            if options is None:
                raise ValueError(
                    f"Unknown IUPAC code '{base}' in pattern '{self.pattern}'."
                )
            expansions = [prefix + opt for prefix in expansions for opt in options]
        return sorted(expansions)


@dataclass
class PatternResult:
    """Result of checking a pattern constraint against a sequence.

    Attributes:
        passed: True if the constraint is satisfied.
            - For ``"avoid"``: True when the pattern is NOT found.
            - For ``"enforce"``: True when the pattern IS found.
        matches: List of ``(start, end)`` tuples for each match.
        pattern: The pattern that was checked.
    """

    passed: bool
    matches: list[tuple[int, int]] = field(default_factory=list)
    pattern: str = ""


# ─── Internal helpers ─────────────────────────────────────────────────────────


def _find_matches_literal(sequence: str, pattern: str) -> list[tuple[int, int]]:
    """Find all literal substring matches, returning (start, end) tuples."""
    results: list[tuple[int, int]] = []
    seq_upper = sequence.upper()
    pat_upper = pattern.upper()
    start = 0
    while True:
        pos = seq_upper.find(pat_upper, start)
        if pos == -1:
            break
        results.append((pos, pos + len(pat_upper)))
        start = pos + 1  # overlapping matches allowed
    return results


def _find_matches_regex(sequence: str, pattern: str) -> list[tuple[int, int]]:
    """Find all regex matches, returning (start, end) tuples."""
    try:
        compiled = re.compile(pattern, re.IGNORECASE)
    except re.error as e:
        raise ValueError(f"Invalid regex pattern '{pattern}': {e}") from None
    results: list[tuple[int, int]] = []
    for m in compiled.finditer(sequence):
        results.append((m.start(), m.end()))
    return results


def _find_matches(
    sequence: str, constraint: PatternConstraint
) -> list[tuple[int, int]]:
    """Find all matches of a pattern on the appropriate strand(s).

    For IUPAC patterns that expand to a small number of concrete
    sequences, we use literal matching for each expansion rather than
    converting to a complex regex.
    """
    all_matches: list[tuple[int, int]] = []

    if constraint.scope == "protein":
        # Protein-scope patterns only match forward
        if constraint.is_regex:
            all_matches.extend(_find_matches_regex(sequence, constraint.pattern))
        else:
            all_matches.extend(_find_matches_literal(sequence, constraint.pattern))
        return all_matches

    # DNA-scope: handle strand and IUPAC expansion
    expansions = constraint.expanded_iupac

    def _match_on_seq(seq: str) -> list[tuple[int, int]]:
        """Match all expansions on a single strand."""
        matches: list[tuple[int, int]] = []
        has_regex_meta = any(
            c in constraint.pattern for c in r"\.^$*+?{}[]|()"
        )
        if len(expansions) == 1 and not has_regex_meta:
            # Simple literal match (fast path)
            matches.extend(_find_matches_literal(seq, expansions[0]))
        elif expansions and len(expansions) <= 16 and not has_regex_meta:
            # IUPAC expansion: do literal match for each concrete sequence
            for exp in expansions:
                matches.extend(_find_matches_literal(seq, exp))
        else:
            # Complex regex or empty expansion (regex metacharacters)
            matches.extend(_find_matches_regex(seq, constraint.pattern))
        return matches

    if constraint.strand in ("forward", "both"):
        all_matches.extend(_match_on_seq(sequence))

    if constraint.strand in ("reverse", "both"):
        rc_seq = reverse_complement(sequence)
        rc_matches = _match_on_seq(rc_seq)
        # Convert RC positions back to forward coordinates
        seq_len = len(sequence)
        for rc_start, rc_end in rc_matches:
            fwd_start = seq_len - rc_end
            fwd_end = seq_len - rc_start
            all_matches.append((fwd_start, fwd_end))

    # Deduplicate overlapping matches
    all_matches.sort()
    return all_matches


# ─── Public API: checking ─────────────────────────────────────────────────────


def check_pattern(
    sequence: str, constraint: PatternConstraint
) -> PatternResult:
    """Check if a sequence satisfies a pattern constraint.

    Args:
        sequence: DNA or protein sequence to check.
        constraint: The pattern constraint to evaluate.

    Returns:
        A :class:`PatternResult` indicating whether the constraint
        passed and where the matches are.
    """
    matches = _find_matches(sequence, constraint)

    if constraint.action == "avoid":
        passed = len(matches) == 0
    else:  # enforce
        passed = len(matches) > 0

    return PatternResult(
        passed=passed,
        matches=matches,
        pattern=constraint.pattern,
    )


def check_patterns(
    sequence: str,
    constraints: list[PatternConstraint],
) -> list[PatternResult]:
    """Check multiple pattern constraints against a sequence.

    Args:
        sequence: DNA or protein sequence to check.
        constraints: List of pattern constraints.

    Returns:
        List of :class:`PatternResult`, one per constraint.
    """
    return [check_pattern(sequence, c) for c in constraints]


# ─── Public API: enforcement ──────────────────────────────────────────────────


def enforce_pattern(
    dna: str,
    protein: str,
    constraint: PatternConstraint,
    codon_table: dict[str, str] | None = None,
) -> str:
    """Modify DNA to enforce or avoid a pattern while preserving translation.

    For **avoid** constraints, synonymous codon substitutions are used
    to remove pattern matches without changing the protein. For **enforce**
    constraints, synonymous codon substitutions are used to *create* the
    required pattern at the first available position.

    Args:
        dna: Input DNA sequence (uppercase ACGT).
        protein: Target protein sequence (single-letter codes).
        constraint: The pattern constraint to enforce.
        codon_table: Codon → amino acid mapping. Defaults to
            :data:`~biocompiler.constants.CODON_TABLE`.

    Returns:
        Modified DNA sequence satisfying the constraint (if possible).
        Returns the original sequence unchanged if no valid modification
        exists.
    """
    if codon_table is None:
        codon_table = CODON_TABLE

    result = check_pattern(dna, constraint)
    if result.passed:
        return dna  # Already satisfied

    if constraint.action == "avoid":
        return _enforce_avoid(dna, protein, constraint, codon_table)
    else:  # enforce
        return _enforce_must_appear(dna, protein, constraint, codon_table)


def enforce_patterns(
    dna: str,
    protein: str,
    constraints: list[PatternConstraint],
    codon_table: dict[str, str] | None = None,
    max_iterations: int = 50,
) -> str:
    """Enforce multiple pattern constraints iteratively.

    Repeatedly applies :func:`enforce_pattern` for each violated
    constraint until all are satisfied or ``max_iterations`` is
    reached.  This handles cascading conflicts where fixing one
    constraint may break another.

    Args:
        dna: Input DNA sequence.
        protein: Target protein sequence.
        constraints: List of pattern constraints.
        codon_table: Codon → amino acid mapping.
        max_iterations: Maximum number of enforcement rounds.

    Returns:
        Modified DNA sequence satisfying all constraints (if possible).
    """
    if codon_table is None:
        codon_table = CODON_TABLE

    current_dna = dna
    for iteration in range(max_iterations):
        all_satisfied = True
        for constraint in constraints:
            result = check_pattern(current_dna, constraint)
            if not result.passed:
                new_dna = enforce_pattern(
                    current_dna, protein, constraint, codon_table
                )
                if new_dna != current_dna:
                    current_dna = new_dna
                    all_satisfied = False
                else:
                    logger.warning(
                        "Could not satisfy constraint %s (action=%s, pattern=%s) "
                        "at iteration %d",
                        constraint.action, constraint.action, constraint.pattern,
                        iteration,
                    )
        if all_satisfied:
            break

    return current_dna


# ─── Avoidance enforcement ────────────────────────────────────────────────────


def _enforce_avoid(
    dna: str,
    protein: str,
    constraint: PatternConstraint,
    codon_table: dict[str, str],
) -> str:
    """Remove all occurrences of a pattern by synonymous codon substitution.

    For each match position, identifies the codon(s) overlapping the
    match and tries alternative synonymous codons that do not contain
    the pattern locally.  Prefers codons that maintain the highest CAI.
    """
    seq_list = list(dna)
    n_codons = len(dna) // 3
    matches = _find_matches(dna, constraint)

    for match_start, match_end in matches:
        # Find codons overlapping this match
        first_codon = match_start // 3
        last_codon = min((match_end - 1) // 3, n_codons - 1)

        fixed = False
        for ci in range(first_codon, last_codon + 1):
            if ci < 0 or ci >= len(protein):
                continue
            aa = protein[ci]
            if aa == "*" or aa not in AA_TO_CODONS:
                continue

            current_codon = "".join(seq_list[ci * 3 : ci * 3 + 3])
            alternatives = AA_TO_CODONS[aa]

            for alt in alternatives:
                if alt == current_codon:
                    continue
                # Apply swap
                for j, base in enumerate(alt):
                    seq_list[ci * 3 + j] = base

                # Check if the pattern is gone locally
                new_dna = "".join(seq_list)
                local_start = max(0, ci * 3 - len(constraint.pattern))
                local_end = min(len(new_dna), ci * 3 + 3 + len(constraint.pattern))
                local_seq = new_dna[local_start:local_end]

                # Check forward and reverse complement as needed
                local_has_pattern = _has_pattern_local(
                    local_seq, constraint, new_dna, local_start, local_end
                )

                if not local_has_pattern:
                    fixed = True
                    break
                else:
                    # Revert
                    for j, base in enumerate(current_codon):
                        seq_list[ci * 3 + j] = base

            if fixed:
                break

        if not fixed:
            logger.debug(
                "Could not avoid pattern '%s' at position %d-%d",
                constraint.pattern, match_start, match_end,
            )

    return "".join(seq_list)


def _has_pattern_local(
    local_seq: str,
    constraint: PatternConstraint,
    full_dna: str,
    local_start: int,
    local_end: int,
) -> bool:
    """Check if a pattern appears in a local region, respecting strand settings."""
    if constraint.is_regex:
        try:
            return bool(re.search(constraint.pattern, local_seq, re.IGNORECASE))
        except re.error:
            return local_seq.upper().find(constraint.pattern.upper()) != -1

    # Literal / IUPAC check
    expansions = constraint.expanded_iupac
    for exp in expansions:
        if exp.upper() in local_seq.upper():
            return True

    # Check reverse complement if needed
    if constraint.strand in ("both", "reverse"):
        rc_local = reverse_complement(local_seq)
        for exp in expansions:
            if exp.upper() in rc_local.upper():
                return True

    return False


# ─── Enforce-must-appear ─────────────────────────────────────────────────────


def _enforce_must_appear(
    dna: str,
    protein: str,
    constraint: PatternConstraint,
    codon_table: dict[str, str],
) -> str:
    """Ensure a pattern appears by synonymous codon substitution.

    Strategy: scan each codon position and try to create the pattern
    by swapping to a synonymous codon that introduces the pattern.
    If the pattern spans a codon boundary, try swapping adjacent codon
    pairs.
    """
    pat_len = len(constraint.pattern)
    n_codons = len(dna) // 3
    seq_list = list(dna)

    # Strategy 1: try to embed pattern within a single codon swap
    for ci in range(len(protein)):
        aa = protein[ci]
        if aa == "*" or aa not in AA_TO_CODONS:
            continue

        current_codon = "".join(seq_list[ci * 3 : ci * 3 + 3])
        alternatives = AA_TO_CODONS[aa]

        for alt in alternatives:
            if alt == current_codon:
                continue
            # Apply swap
            for j, base in enumerate(alt):
                seq_list[ci * 3 + j] = base

            new_dna = "".join(seq_list)
            result = check_pattern(new_dna, constraint)
            if result.passed:
                return new_dna

            # Revert
            for j, base in enumerate(current_codon):
                seq_list[ci * 3 + j] = base

    # Strategy 2: try adjacent codon pair swaps for longer patterns
    for ci in range(len(protein) - 1):
        aa1 = protein[ci]
        aa2 = protein[ci + 1]
        if aa1 == "*" or aa2 == "*":
            continue
        if aa1 not in AA_TO_CODONS or aa2 not in AA_TO_CODONS:
            continue

        codon1 = "".join(seq_list[ci * 3 : ci * 3 + 3])
        codon2 = "".join(seq_list[(ci + 1) * 3 : (ci + 1) * 3 + 3])

        for alt1 in AA_TO_CODONS[aa1]:
            if alt1 == codon1:
                continue
            for alt2 in AA_TO_CODONS[aa2]:
                if alt2 == codon2:
                    continue
                # Apply both swaps
                for j, base in enumerate(alt1):
                    seq_list[ci * 3 + j] = base
                for j, base in enumerate(alt2):
                    seq_list[(ci + 1) * 3 + j] = base

                new_dna = "".join(seq_list)
                result = check_pattern(new_dna, constraint)
                if result.passed:
                    return new_dna

                # Revert both
                for j, base in enumerate(codon1):
                    seq_list[ci * 3 + j] = base
                for j, base in enumerate(codon2):
                    seq_list[(ci + 1) * 3 + j] = base

    logger.warning(
        "Could not enforce pattern '%s' (action=enforce) — no valid "
        "codon substitution found that introduces the pattern while "
        "preserving translation",
        constraint.pattern,
    )
    return dna


# ─── Aho-Corasick multi-pattern avoidance ─────────────────────────────────────


def build_avoidance_scanner(
    constraints: list[PatternConstraint],
) -> Optional[AhoCorasickScanner]:
    """Build an Aho-Corasick scanner from multiple avoid-pattern constraints.

    Only includes constraints with:
    - ``action == "avoid"``
    - ``scope == "dna"``
    - Pure-ACGT literal patterns (no regex/IUPAC)

    For IUPAC patterns, expands them into concrete ACGT sequences
    before adding to the scanner.

    Args:
        constraints: List of pattern constraints.

    Returns:
        An :class:`~biocompiler.aho_corasick.AhoCorasickScanner` ready
        to scan, or ``None`` if no concrete ACGT patterns were found.
    """
    patterns: dict[str, str] = {}

    for constraint in constraints:
        if constraint.action != "avoid":
            continue
        if constraint.scope != "dna":
            continue
        # Skip regex patterns (contain metacharacters)
        if any(c in constraint.pattern for c in r"\.^$*+?{}[]|()"):
            continue

        # Expand IUPAC and add concrete sequences
        try:
            expansions = constraint.expanded_iupac
        except ValueError:
            continue

        for exp in expansions:
            patterns[exp] = constraint.pattern
            # Add reverse complement if strand == "both"
            if constraint.strand in ("both", "reverse"):
                rc = reverse_complement(exp)
                if rc != exp:
                    patterns[rc] = constraint.pattern

    if not patterns:
        return None

    return AhoCorasickScanner(patterns)
