"""
BioCompiler MaxEntScan CSP Encoding
====================================

Encodes MaxEntScan splice scoring as CSP (Constraint Satisfaction Problem)
constraints, enabling the formal solver to avoid cryptic splice sites
instead of the greedy optimizer's iterative fix-up approach.

Core Challenge
--------------
MaxEntScan scores are continuous values computed from 9-mer (donor) and
23-mer (acceptor) sequence context, but CSP solvers work with discrete/
integer variables. This module bridges the gap by:

1. **Pre-computation**: Computing all possible splice scores for every
   codon choice BEFORE building the CSP model. This converts the
   continuous scoring function into a lookup table.

2. **Table constraints**: For OR-Tools, encoding "forbidden codon
   assignments" as table constraints — the solver simply cannot select
   a codon that would create a cryptic splice site.

3. **Quantization**: For soft constraints (where we want to minimize
   the worst splice score rather than forbid it), quantizing continuous
   scores into discrete bins for the objective function.

4. **Cross-codon encoding**: GT/AG dinucleotides that span codon
   boundaries require constraints on *pairs* of adjacent codon variables.
   This is the key advantage over the greedy optimizer, which struggles
   with these interactions.

Key Insight
-----------
The MaxEntScan score for a GT at position *p* depends on the 9-mer
context (3 bases before GT, GT itself, 6 bases after). This means
changing one codon can affect splice scores for GT dinucleotides up to
3 codons away. Similarly, the 23-mer acceptor context means one codon
change can affect acceptor scores up to 7 codons away.

For OR-Tools CP-SAT: encode as table constraints on tuples of adjacent
codons (arity-2 constraints for cross-boundary sites).

For Z3 SMT: encode as if-then-else chains on the codon string variables.

References
----------
Yeo & Burge (2004) "Maximum entropy modeling of short sequence motifs
with applications to RNA splicing" J Comp Biol 11(2-3):377-94.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from itertools import product
from typing import Optional

from ..constants import AA_TO_CODONS
from ..maxentscan import score_donor, score_acceptor
from .types import SolverConfig, SpliceConstraint, CrossCodonSpliceConstraint

logger = logging.getLogger(__name__)

# Sentinel value used by MaxEntScan when the scoring context is out of bounds.
_OUT_OF_BOUNDS_SCORE: float = -50.0


# ==============================================================================
# 1. Pre-compute splice scores
# ==============================================================================


def _build_sequence_with_codon(
    protein: str,
    codon_assignments: dict[int, str],
    position: int,
    codon: str,
) -> str:
    """Build a DNA sequence string with a specific codon at a specific position.

    Uses a reference context from codon_assignments for all other positions.
    Positions not in codon_assignments use the first (highest-CAI) codon.

    Args:
        protein: Target protein sequence.
        codon_assignments: Partial codon assignments {position: codon}.
        position: The position to override.
        codon: The codon to place at position.

    Returns:
        DNA sequence string of length len(protein) * 3.
    """
    seq_parts: list[str] = []
    assignments = dict(codon_assignments)
    assignments[position] = codon
    for i in range(len(protein)):
        if i in assignments:
            seq_parts.append(assignments[i])
        else:
            seq_parts.append(AA_TO_CODONS[protein[i]][0])
    return "".join(seq_parts)


def precompute_splice_scores(
    protein: str,
    config: SolverConfig,
) -> dict[tuple[int, str], dict]:
    """Pre-compute MaxEntScan splice scores for every codon choice at every position.

    For each codon position *i* and each possible codon choice *c_i*, this
    function computes the MaxEntScan donor and acceptor scores for all GT/AG
    dinucleotides that would be created or affected by this codon choice.

    This pre-computation is done BEFORE building the CSP model, converting
    the continuous MaxEntScan scoring function into a finite lookup table
    that the CSP solver can reason about.

    Algorithm
    ---------
    For each position i and codon c:
    1. Build a candidate DNA sequence with c at position i (highest-CAI
       codons at all other positions as default context).
    2. Scan the local context window for GT and AG dinucleotides whose
       scoring context (9-mer for donors, 23-mer for acceptors) overlaps
       with the bases of codon i.
    3. Score each GT/AG found with MaxEntScan.
    4. Record all scores where the scoring context overlaps with position i.

    Args:
        protein: Target protein sequence (single-letter amino acid codes).
        config: Solver configuration with splice thresholds.

    Returns:
        Dictionary mapping ``(position, codon)`` to a dict with keys:
        - ``"donor_scores"``: list of ``(base_position, score)`` for GT sites
          whose 9-mer context overlaps with this codon's bases.
        - ``"acceptor_scores"``: list of ``(base_position, score)`` for AG sites
          whose 23-mer context overlaps with this codon's bases.

    Example
    -------
    >>> config = SolverConfig()
    >>> data = precompute_splice_scores("MV", config)
    >>> (0, "ATG") in data  # position 0, codon ATG for Methionine
    True
    """
    protein = protein.upper()
    n_codons = len(protein)
    result: dict[tuple[int, str], dict] = {}

    # Build a default sequence (highest-CAI codon per position) for context
    default_assignments: dict[int, str] = {}
    for i, aa in enumerate(protein):
        if aa in AA_TO_CODONS:
            default_assignments[i] = AA_TO_CODONS[aa][0]

    for i in range(n_codons):
        aa = protein[i]
        if aa not in AA_TO_CODONS:
            logger.warning("Unknown amino acid '%s' at position %d, skipping", aa, i)
            continue

        for codon in AA_TO_CODONS[aa]:
            # Build sequence with this codon at position i
            seq = _build_sequence_with_codon(protein, default_assignments, i, codon)

            # Scan for donor sites (GT dinucleotides) in the affected range.
            # Donor 9-mer: bases [p-3, p+8). A GT at base p is affected by
            # codon i if codon i's bases [3i, 3i+3) overlap with [p-3, p+8).
            codon_start = i * 3
            codon_end = i * 3 + 3

            donor_scores: list[tuple[int, float]] = []
            gt_lo = max(0, codon_start - 8)
            gt_hi = min(len(seq) - 1, codon_end + 2)

            for p in range(gt_lo, gt_hi):
                if p + 1 < len(seq) and seq[p] == "G" and seq[p + 1] == "T":
                    s = score_donor(seq, p)
                    mer_start = p - 3
                    mer_end = p + 8
                    if mer_start < codon_end and mer_end > codon_start:
                        donor_scores.append((p, s))

            # Scan for acceptor sites (AG dinucleotides) in the affected range.
            # Acceptor 23-mer: bases [p-20, p+3), AG at base p.
            acceptor_scores: list[tuple[int, float]] = []
            ag_lo = max(0, codon_start - 22)
            ag_hi = min(len(seq) - 1, codon_end + 20)

            for p in range(ag_lo, ag_hi):
                if p + 1 < len(seq) and seq[p] == "A" and seq[p + 1] == "G":
                    s = score_acceptor(seq, p)
                    mer_start = p - 20
                    mer_end = p + 3
                    if mer_start < codon_end and mer_end > codon_start:
                        acceptor_scores.append((p, s))

            result[(i, codon)] = {
                "donor_scores": donor_scores,
                "acceptor_scores": acceptor_scores,
            }

    logger.debug(
        "Pre-computed splice scores for %d (position, codon) pairs",
        len(result),
    )
    return result


# ==============================================================================
# 2. Build splice constraint table
# ==============================================================================


def build_splice_constraint_table(
    protein: str,
    config: SolverConfig,
    score_data: dict[tuple[int, str], dict],
) -> list[SpliceConstraint]:
    """Build a table of forbidden codon assignments based on splice scores.

    A codon choice is FORBIDDEN if placing that codon at its position would
    create a cryptic splice site — i.e., a GT with MaxEntScan donor score
    above the donor threshold, or an AG with MaxEntScan acceptor score
    above the acceptor threshold.

    This is the hard constraint encoding: the solver must never assign
    a forbidden codon to a position. For soft constraint variants (penalize
    high scores rather than forbid them), use quantize_maxent_scores() and
    add the quantized scores to the objective function.

    Args:
        protein: Target protein sequence.
        config: Solver configuration with splice thresholds.
        score_data: Pre-computed splice scores from precompute_splice_scores().

    Returns:
        List of SpliceConstraint objects, one per forbidden (position, codon)
        assignment. Each constraint records which splice site would be
        created and its score. Deduplicated: only the worst score per
        (position, codon, site_type) triple is kept.

    Example
    -------
    >>> config = SolverConfig(cryptic_splice_threshold=3.0)
    >>> data = precompute_splice_scores("MV", config)
    >>> constraints = build_splice_constraint_table("MV", config, data)
    >>> # All constraints have score > threshold
    >>> all(c.score > c.threshold for c in constraints)
    True
    """
    protein = protein.upper()
    donor_thresh = config.effective_donor_threshold
    acceptor_thresh = config.effective_acceptor_threshold
    forbidden: list[SpliceConstraint] = []

    for (pos, codon), scores in score_data.items():
        # Check donor scores
        for base_pos, d_score in scores["donor_scores"]:
            if d_score > donor_thresh:
                forbidden.append(SpliceConstraint(
                    position=pos,
                    codon=codon,
                    site_type="donor",
                    score=d_score,
                    threshold=donor_thresh,
                ))

        # Check acceptor scores
        for base_pos, a_score in scores["acceptor_scores"]:
            if a_score > acceptor_thresh:
                forbidden.append(SpliceConstraint(
                    position=pos,
                    codon=codon,
                    site_type="acceptor",
                    score=a_score,
                    threshold=acceptor_thresh,
                ))

    # Deduplicate: same (pos, codon, site_type) may appear from multiple
    # GT/AG positions. Keep only the worst score per key.
    best_per_key: dict[tuple[int, str, str], SpliceConstraint] = {}
    for c in forbidden:
        key = (c.position, c.codon, c.site_type)
        if key not in best_per_key or c.score > best_per_key[key].score:
            best_per_key[key] = c

    result = list(best_per_key.values())
    logger.debug(
        "Built splice constraint table: %d forbidden (position, codon) assignments",
        len(result),
    )
    return result


# ==============================================================================
# 3. Quantize MaxEntScan scores
# ==============================================================================


def quantize_maxent_scores(
    scores: list[float],
    n_bins: int = 20,
    threshold: float = 3.0,
) -> list[int]:
    """Quantize continuous MaxEntScan scores into discrete bins for CSP encoding.

    Converts a list of continuous MaxEntScan scores (which can range from
    about -50 to ~12 for donors and ~14 for acceptors) into integer bin
    indices. This allows CSP solvers to reason about splice scores in
    their objective function without continuous variables.

    Binning Strategy
    ----------------
    - Bins are defined over [0, threshold] in equal-width intervals.
    - Scores < 0 map to bin 0 (very weak, no splice activity).
    - Scores >= threshold map to bin n_bins (very strong, cryptic site).
    - Intermediate scores map linearly: bin = floor(score / (threshold / n_bins))

    This non-linear mapping concentrates resolution near the threshold,
    which is exactly where the decision boundary matters most. Very weak
    and very strong sites don't need fine-grained distinction.

    Args:
        scores: List of continuous MaxEntScan scores.
        n_bins: Number of bins for quantization (default 20).
        threshold: The splice site threshold; scores above this are
            all mapped to the maximum bin (default 3.0).

    Returns:
        List of integer bin indices in [0, n_bins]. Bin 0 = very weak,
        bin n_bins = very strong (above threshold).

    Example
    -------
    >>> quantize_maxent_scores([0.0, 1.5, 3.0, 5.0], n_bins=10, threshold=3.0)
    [0, 5, 10, 10]
    >>> quantize_maxent_scores([-5.0, 0.3], n_bins=20, threshold=3.0)
    [0, 2]
    """
    if n_bins < 2:
        raise ValueError(f"Need at least 2 bins, got {n_bins}")
    if threshold <= 0:
        raise ValueError(f"Threshold must be positive, got {threshold}")

    bin_width = threshold / n_bins
    result: list[int] = []

    for score in scores:
        if score <= 0.0:
            result.append(0)
        elif score >= threshold:
            result.append(n_bins)
        else:
            bin_idx = int(score / bin_width)
            # Clamp to [0, n_bins]
            bin_idx = max(0, min(n_bins, bin_idx))
            result.append(bin_idx)

    return result


# ==============================================================================
# 4. Encode cross-codon splice context
# ==============================================================================


def encode_cross_codon_splice_context(
    protein: str,
    config: SolverConfig,
) -> list[tuple[int, str, str, bool, float]]:
    """Enumerate cross-codon GT/AG dinucleotides that span codon boundaries.

    Cross-codon splice sites occur when a GT or AG dinucleotide straddles
    the boundary between two adjacent codons. These are the hardest cases
    for the greedy optimizer because fixing one codon alone doesn't
    eliminate the dinucleotide — you need to change BOTH codons.

    This function enumerates all (codon_left, codon_right) pairs at each
    adjacent codon boundary that would create a GT or AG at the boundary.

    Boundary dinucleotides
    ----------------------
    At the boundary between codon i (bases 3i, 3i+1, 3i+2) and codon i+1
    (bases 3i+3, 3i+4, 3i+5), the cross-codon dinucleotide is formed by:
    - base at position 3i+2 (last base of codon i)
    - base at position 3i+3 (first base of codon i+1)

    So we check: does codon_left[-1] + codon_right[0] == "GT" or "AG"?

    For each such pair, we compute the MaxEntScan score in the context of
    a default surrounding sequence (highest-CAI codons at other positions).

    Args:
        protein: Target protein sequence.
        config: Solver configuration (used for thresholds in downstream filtering).

    Returns:
        List of tuples ``(position, codon_left, codon_right, is_donor, score)``:
        - position: the left codon position (boundary is between position and position+1)
        - codon_left: the 3-letter codon at the left position
        - codon_right: the 3-letter codon at the right position
        - is_donor: True if the boundary dinucleotide is GT (donor), False if AG (acceptor)
        - score: the MaxEntScan score for this splice site in context

    Example
    -------
    >>> config = SolverConfig()
    >>> cross = encode_cross_codon_splice_context("MV", config)
    >>> # Each entry has a cross-boundary GT or AG
    >>> all(len(t[1]) == 3 and len(t[2]) == 3 for t in cross)
    True
    """
    protein = protein.upper()
    n_codons = len(protein)
    result: list[tuple[int, str, str, bool, float]] = []

    # Build default sequence for context
    default_seq = "".join(
        AA_TO_CODONS[protein[i]][0] for i in range(n_codons)
    )

    for i in range(n_codons - 1):
        aa_left = protein[i]
        aa_right = protein[i + 1]

        if aa_left not in AA_TO_CODONS or aa_right not in AA_TO_CODONS:
            continue

        codons_left = AA_TO_CODONS[aa_left]
        codons_right = AA_TO_CODONS[aa_right]

        # Boundary base position: last base of codon i
        boundary_base = i * 3 + 2

        for cl, cr in product(codons_left, codons_right):
            # Check cross-boundary dinucleotide: last base of cl + first base of cr
            dinuc = cl[2] + cr[0]

            if dinuc == "GT":
                # Build the full sequence with this pair
                seq_list = list(default_seq)
                seq_list[i * 3: i * 3 + 3] = list(cl)
                seq_list[(i + 1) * 3: (i + 1) * 3 + 3] = list(cr)
                seq = "".join(seq_list)

                score = score_donor(seq, boundary_base)
                result.append((i, cl, cr, True, score))

            elif dinuc == "AG":
                seq_list = list(default_seq)
                seq_list[i * 3: i * 3 + 3] = list(cl)
                seq_list[(i + 1) * 3: (i + 1) * 3 + 3] = list(cr)
                seq = "".join(seq_list)

                score = score_acceptor(seq, boundary_base)
                result.append((i, cl, cr, False, score))

    # Filter out entries where context was out of bounds (score = sentinel)
    result = [t for t in result if t[4] > _OUT_OF_BOUNDS_SCORE]

    logger.debug(
        "Encoded %d cross-codon splice context entries across %d boundaries",
        len(result),
        n_codons - 1,
    )
    return result


# ==============================================================================
# 5. SpliceConstraintEncoder class
# ==============================================================================


@dataclass
class SpliceConstraintEncoder:
    """Encapsulates MaxEntScan splice scoring as CSP constraints.

    This class provides a high-level interface for encoding splice site
    avoidance into a CSP model. It caches precomputed scores for reuse
    across multiple solver invocations and supports both OR-Tools table
    constraints and Z3 if-then-else encoding.

    Usage
    -----
    >>> config = SolverConfig()
    >>> encoder = SpliceConstraintEncoder("MVLSPADKTN", config)
    >>> # Get forbidden codon assignments
    >>> forbidden = encoder.get_forbidden_assignments()
    >>> # Get cross-codon constraints
    >>> cross = encoder.get_cross_codon_constraints()
    >>> # Build OR-Tools model (requires ortools installed)
    >>> encoder.encode_or_tools_table_constraints(model, codon_vars, domains)

    Attributes:
        protein: Target protein sequence.
        config: Solver configuration.
        _score_cache: Cached pre-computed splice scores.
        _forbidden_cache: Cached forbidden assignments.
        _cross_codon_cache: Cached cross-codon raw data.
        _cross_codon_constraint_cache: Cached CrossCodonSpliceConstraint objects.
        _initialized: Whether initialization has been run.
    """

    protein: str
    config: SolverConfig
    _score_cache: dict[tuple[int, str], dict] = field(default_factory=dict, repr=False)
    _forbidden_cache: list[SpliceConstraint] = field(default_factory=list, repr=False)
    _cross_codon_cache: list[tuple[int, str, str, bool, float]] = field(
        default_factory=list, repr=False
    )
    _cross_codon_constraint_cache: list[CrossCodonSpliceConstraint] = field(
        default_factory=list, repr=False
    )
    _initialized: bool = field(default=False, repr=False)

    def initialize(self) -> None:
        """Run pre-computation and build all caches.

        This must be called before any other query method. It computes splice
        scores, builds the forbidden assignment table, and enumerates
        cross-codon constraints. Results are cached for reuse.

        The computation is deterministic: same protein + config always
        produces same caches.
        """
        if self._initialized:
            return

        protein = self.protein.upper()

        # Step 1: Pre-compute splice scores
        self._score_cache = precompute_splice_scores(protein, self.config)

        # Step 2: Build forbidden assignment table
        self._forbidden_cache = build_splice_constraint_table(
            protein, self.config, self._score_cache
        )

        # Step 3: Enumerate cross-codon constraints
        self._cross_codon_cache = encode_cross_codon_splice_context(
            protein, self.config
        )

        # Step 3b: Convert to CrossCodonSpliceConstraint objects (filtered by threshold)
        self._cross_codon_constraint_cache = self._build_cross_codon_constraints()

        self._initialized = True
        logger.info(
            "SpliceConstraintEncoder initialized: %d score entries, "
            "%d forbidden assignments, %d cross-codon constraints",
            len(self._score_cache),
            len(self._forbidden_cache),
            len(self._cross_codon_constraint_cache),
        )

    def _build_cross_codon_constraints(self) -> list[CrossCodonSpliceConstraint]:
        """Build CrossCodonSpliceConstraint objects from cross-codon data.

        Only includes pairs where the MaxEntScan score exceeds the
        relevant threshold (donor or acceptor).
        """
        constraints: list[CrossCodonSpliceConstraint] = []
        donor_thresh = self.config.effective_donor_threshold
        acceptor_thresh = self.config.effective_acceptor_threshold

        for pos, cl, cr, is_donor, score in self._cross_codon_cache:
            threshold = donor_thresh if is_donor else acceptor_thresh
            if score > threshold:
                constraints.append(CrossCodonSpliceConstraint(
                    position_left=pos,
                    position_right=pos + 1,
                    codon_left=cl,
                    codon_right=cr,
                    is_donor=is_donor,
                    score=score,
                ))
        return constraints

    # ────────────────────────────────────────────────────────────
    # Public query methods
    # ────────────────────────────────────────────────────────────

    def get_forbidden_assignments(self) -> list[SpliceConstraint]:
        """Return the list of forbidden (position, codon) assignments.

        Each SpliceConstraint records a codon that cannot be used at a
        position because it would create a cryptic splice site.

        Returns:
            List of SpliceConstraint objects.
        """
        self.initialize()
        return list(self._forbidden_cache)

    def get_cross_codon_constraints(self) -> list[CrossCodonSpliceConstraint]:
        """Return cross-codon constraints for adjacent codon pairs.

        These constraints forbid specific (codon_left, codon_right) pairs
        at adjacent positions where the boundary dinucleotide forms a
        cryptic splice site.

        Returns:
            List of CrossCodonSpliceConstraint objects.
        """
        self.initialize()
        return list(self._cross_codon_constraint_cache)

    def get_scores_for_position(
        self, position: int, codon: str
    ) -> Optional[dict]:
        """Get pre-computed splice scores for a specific (position, codon) pair.

        Args:
            position: 0-indexed codon position.
            codon: 3-letter DNA codon string.

        Returns:
            Dict with "donor_scores" and "acceptor_scores" lists, or None
            if the (position, codon) pair was not computed.
        """
        self.initialize()
        return self._score_cache.get((position, codon))

    def get_max_splice_score(
        self, position: int, codon: str
    ) -> tuple[float, float]:
        """Get the maximum donor and acceptor scores for a (position, codon) pair.

        Useful for building soft constraint objectives — minimize the
        worst-case splice score across all positions.

        Args:
            position: 0-indexed codon position.
            codon: 3-letter DNA codon string.

        Returns:
            (max_donor_score, max_acceptor_score) tuple. Returns
            (_OUT_OF_BOUNDS_SCORE, _OUT_OF_BOUNDS_SCORE) if no scores
            are available.
        """
        scores = self.get_scores_for_position(position, codon)
        if scores is None:
            return (_OUT_OF_BOUNDS_SCORE, _OUT_OF_BOUNDS_SCORE)

        max_donor = max((s for _, s in scores["donor_scores"]), default=_OUT_OF_BOUNDS_SCORE)
        max_acceptor = max((s for _, s in scores["acceptor_scores"]), default=_OUT_OF_BOUNDS_SCORE)
        return (max_donor, max_acceptor)

    def quantize_position_scores(
        self, position: int, codon: str
    ) -> tuple[int, int]:
        """Get quantized splice scores for a (position, codon) pair.

        Convenience method combining get_max_splice_score and
        quantize_maxent_scores.

        Args:
            position: 0-indexed codon position.
            codon: 3-letter DNA codon string.

        Returns:
            (donor_bin, acceptor_bin) tuple with integer bin indices.
        """
        max_d, max_a = self.get_max_splice_score(position, codon)
        donor_bins = quantize_maxent_scores(
            [max_d],
            n_bins=self.config.n_quantize_bins,
            threshold=self.config.effective_donor_threshold,
        )
        acceptor_bins = quantize_maxent_scores(
            [max_a],
            n_bins=self.config.n_quantize_bins,
            threshold=self.config.effective_acceptor_threshold,
        )
        return (donor_bins[0], acceptor_bins[0])

    # ────────────────────────────────────────────────────────────
    # OR-Tools encoding
    # ────────────────────────────────────────────────────────────

    def encode_or_tools_table_constraints(
        self,
        model: "cp_model.CpModel",
        codon_vars: list["cp_model.IntVar"],
        codon_domains: list[list[str]],
    ) -> None:
        """Encode splice constraints as OR-Tools table constraints.

        For within-codon constraints (a single codon creates a cryptic site),
        we add inequality constraints excluding forbidden codon indices from
        each variable's domain.

        For cross-codon constraints (adjacent codons create a boundary GT/AG),
        we add Boolean reified constraints that forbid the specific
        (codon_i, codon_i+1) value pair.

        Args:
            model: OR-Tools CpModel instance.
            codon_vars: List of OR-Tools IntVar variables, one per codon position.
                Each variable's domain is an index into the corresponding
                codon_domains entry.
            codon_domains: List of codon domain lists. codon_domains[i] is the
                list of allowed codon strings at position i. The IntVar value
                is an index into this list.

        Raises:
            ImportError: If ortools is not installed.
        """
        from ortools.sat.python import cp_model  # noqa: F811

        self.initialize()
        n_codons = len(self.protein)

        # Build a mapping: position -> codon -> domain index
        codon_to_idx: dict[int, dict[str, int]] = {}
        for pos in range(n_codons):
            codon_to_idx[pos] = {}
            for idx, codon in enumerate(codon_domains[pos]):
                codon_to_idx[pos][codon] = idx

        # Within-codon forbidden assignments:
        # Add var != forbidden_index for each forbidden codon.
        within_forbidden: dict[int, set[int]] = {}
        for fc in self._forbidden_cache:
            if fc.position not in within_forbidden:
                within_forbidden[fc.position] = set()
            if fc.codon in codon_to_idx.get(fc.position, {}):
                within_forbidden[fc.position].add(
                    codon_to_idx[fc.position][fc.codon]
                )

        for pos, forbidden_indices in within_forbidden.items():
            if pos < n_codons:
                all_indices = set(range(len(codon_domains[pos])))
                allowed = all_indices - forbidden_indices
                if allowed:
                    for fi in forbidden_indices:
                        model.Add(codon_vars[pos] != fi)
                else:
                    logger.warning(
                        "Position %d: ALL codons forbidden by splice constraints! "
                        "This protein position may require amino acid substitution.",
                        pos,
                    )

        # Cross-codon forbidden pairs:
        # For each CrossCodonSpliceConstraint, add a forbidden tuple constraint.
        cross_pairs: dict[tuple[int, int], list[tuple[int, int]]] = {}
        for cc in self._cross_codon_constraint_cache:
            key = (cc.position_left, cc.position_right)
            if key not in cross_pairs:
                cross_pairs[key] = []
            left_idx = codon_to_idx.get(cc.position_left, {}).get(cc.codon_left)
            right_idx = codon_to_idx.get(cc.position_right, {}).get(cc.codon_right)
            if left_idx is not None and right_idx is not None:
                cross_pairs[key].append((left_idx, right_idx))

        for (pos_l, pos_r), forbidden_tuples in cross_pairs.items():
            if pos_l < n_codons and pos_r < n_codons and forbidden_tuples:
                # Encode: NOT (vars[pos_l] == a AND vars[pos_r] == b)
                # Equivalently: vars[pos_l] != a OR vars[pos_r] != b
                for left_val, right_val in forbidden_tuples:
                    b1 = model.NewBoolVar(
                        f"splice_not_{pos_l}_{left_val}"
                    )
                    b2 = model.NewBoolVar(
                        f"splice_not_{pos_r}_{right_val}"
                    )
                    model.Add(codon_vars[pos_l] != left_val).OnlyEnforceIf(b1)
                    model.Add(codon_vars[pos_r] != right_val).OnlyEnforceIf(b2)
                    # At least one must be different
                    model.AddBoolOr([b1, b2])

    # ────────────────────────────────────────────────────────────
    # Z3 encoding
    # ────────────────────────────────────────────────────────────

    def encode_z3_constraints(
        self,
        solver: "z3.Solver",
        codon_vars: list["z3.ExprRef"],
        codon_domains: list[list[str]],
    ) -> None:
        """Encode splice constraints as Z3 if-then-else constraints.

        For each forbidden (position, codon) assignment, adds an assertion
        that the codon variable at that position must NOT equal the forbidden
        codon's index.

        For cross-codon constraints, adds assertions that adjacent codon
        variables must NOT simultaneously take on the forbidden pair values.

        Args:
            solver: Z3 Solver instance.
            codon_vars: List of Z3 integer variables, one per codon position.
            codon_domains: List of codon domain lists (same as OR-Tools).

        Raises:
            ImportError: If z3-solver is not installed.
        """
        import z3

        self.initialize()
        n_codons = len(self.protein)

        # Build index mapping
        codon_to_idx: dict[int, dict[str, int]] = {}
        for pos in range(n_codons):
            codon_to_idx[pos] = {}
            for idx, codon in enumerate(codon_domains[pos]):
                codon_to_idx[pos][codon] = idx

        # Within-codon: add NOT(var == forbidden_index) for each forbidden codon
        for fc in self._forbidden_cache:
            if fc.position < n_codons:
                idx = codon_to_idx.get(fc.position, {}).get(fc.codon)
                if idx is not None:
                    solver.add(
                        codon_vars[fc.position] != idx,
                        f"splice_{fc.site_type}_pos{fc.position}_codon{fc.codon}",
                    )

        # Cross-codon: add NOT(var_left == a AND var_right == b)
        for cc in self._cross_codon_constraint_cache:
            left_idx = codon_to_idx.get(cc.position_left, {}).get(cc.codon_left)
            right_idx = codon_to_idx.get(cc.position_right, {}).get(cc.codon_right)
            if left_idx is not None and right_idx is not None:
                site_label = "donor" if cc.is_donor else "acceptor"
                solver.add(
                    z3.Or(
                        codon_vars[cc.position_left] != left_idx,
                        codon_vars[cc.position_right] != right_idx,
                    ),
                    f"splice_cross_{site_label}_pos{cc.position_left}_{cc.position_right}",
                )

    # ────────────────────────────────────────────────────────────
    # Objective function helpers
    # ────────────────────────────────────────────────────────────

    def build_or_tools_splice_objective(
        self,
        model: "cp_model.CpModel",
        codon_vars: list["cp_model.IntVar"],
        codon_domains: list[list[str]],
    ) -> "cp_model.IntVar":
        """Build a soft constraint objective that minimizes the worst splice score.

        Creates auxiliary integer variables representing the quantized splice
        score at each position, and returns a variable representing the
        maximum across all positions. The solver should minimize this.

        Args:
            model: OR-Tools CpModel instance.
            codon_vars: List of OR-Tools IntVar codon choice variables.
            codon_domains: List of codon domain lists.

        Returns:
            An IntVar representing the maximum quantized splice score.
            The solver should minimize this to avoid cryptic splice sites.

        Raises:
            ImportError: If ortools is not installed.
        """
        from ortools.sat.python import cp_model  # noqa: F811

        self.initialize()
        n_codons = len(self.protein)

        # For each position, create an auxiliary variable for the splice score
        max_bin = self.config.n_quantize_bins
        score_vars = []

        for pos in range(n_codons):
            # Build a lookup: codon_domain_index -> (donor_bin, acceptor_bin)
            domain = codon_domains[pos]
            lookup: dict[int, tuple[int, int]] = {}
            for idx, codon in enumerate(domain):
                d_bin, a_bin = self.quantize_position_scores(pos, codon)
                lookup[idx] = (d_bin, a_bin)

            # Create the score variable and link it to the codon choice
            score_var = model.NewIntVar(0, max_bin, f"splice_score_{pos}")

            for idx, (d_bin, a_bin) in lookup.items():
                worst = max(d_bin, a_bin)
                b = model.NewBoolVar(f"splice_score_{pos}_is_{idx}")
                model.Add(codon_vars[pos] == idx).OnlyEnforceIf(b)
                model.Add(score_var == worst).OnlyEnforceIf(b)

            score_vars.append(score_var)

        # Max of all score_vars
        max_score = model.NewIntVar(0, max_bin, "max_splice_score")
        model.AddMaxEquality(max_score, score_vars)
        return max_score

    def clear_cache(self) -> None:
        """Clear all cached pre-computation results.

        Useful when the configuration changes and caches need to be rebuilt.
        """
        self._score_cache.clear()
        self._forbidden_cache.clear()
        self._cross_codon_cache.clear()
        self._cross_codon_constraint_cache.clear()
        self._initialized = False
