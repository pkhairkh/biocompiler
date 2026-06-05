"""
BioCompiler CSP Constraint Model
==================================

Formalizes all gene optimization constraints as a Constraint Satisfaction
Problem (CSP) model.  The key insight (DOC-10 S5): gene optimization is
CONSTRAINT SATISFACTION, not optimization.  We seek ANY solution satisfying
all hard constraints, then optimize CAI within the feasible region.

Constraint hierarchy:
    HARD (must hold — infeasible if violated):
        - TranslationConstraint      translate(codon_i) = protein[i]
        - NoRestrictionSiteConstraint  no forbidden substrings
        - GCRangeConstraint           gc_lo <= GC <= gc_hi
        - NoCrypticSpliceConstraint   MaxEnt score < threshold at all GT/AG
        - NoCpGIslandConstraint       CpG Obs/Exp ratio <= threshold
        - NoATTTAMotifConstraint      no ATTTA instability motifs
        - NoTRunConstraint            no 6+ consecutive T bases

    SOFT / objective (MAXIMIZE):
        - MaximizeCAI     sum_i log(adaptiveness(codon_i))
        - MinimizeCpG     minimize CG dinucleotide count
        - MinimizeMRNADG  minimize mRNA folding dG at 5' end

Usage::

    from biocompiler.solver.constraints import build_csp_model

    model = build_csp_model(protein="MVSKGE", organism="Homo_sapiens", config=cfg)
    # model.variables  -> [CodonVariable(...), ...]
    # model.hard_constraints -> [TranslationConstraint(...), ...]
    # model.soft_constraints -> [MaximizeCAI(...), ...]
"""

from __future__ import annotations

import math
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, ClassVar, Sequence

from .types import CodonVariable, ConstraintPriority, ConstraintSpec, SolverConfig, ConstraintStrictness, ConstraintType
from ..constants import AA_TO_CODONS, CODON_TABLE, RESTRICTION_ENZYMES, INSTABILITY_MOTIF
from ..organisms import (
    CODON_ADAPTIVENESS_TABLES,
    ORGANISM_GC_TARGETS,
    get_sharp_li_adaptiveness_tables,
    SUPPORTED_ORGANISMS,
)
from ..organism_config import is_eukaryotic_organism
from ..maxentscan import score_donor, score_acceptor

logger = logging.getLogger(__name__)


__all__ = [
    # Constants
    "DEFAULT_GC_LO",
    "DEFAULT_GC_HI",
    "DEFAULT_CPG_WINDOW",
    "DEFAULT_CPG_THRESHOLD",
    "DEFAULT_MAX_T_RUN",
    "CAI_LOG_EPSILON",
    "NEAREST_NEIGHBOR_GC",
    "NEAREST_NEIGHBOR_AU",
    "NEAREST_NEIGHBOR_GU",
    # Helper functions
    "codon_gc_count",
    "codon_contains_gt",
    "codon_contains_ag",
    "codon_contains_cpg",
    "compute_gc_from_codons",
    # Abstract base classes
    "HardConstraint",
    "SoftConstraint",
    # Hard constraints
    "TranslationConstraint",
    "NoRestrictionSiteConstraint",
    "GCRangeConstraint",
    "NoCrypticSpliceConstraint",
    "NoCpGIslandConstraint",
    "NoATTTAMotifConstraint",
    "NoTRunConstraint",
    # Soft constraints
    "MaximizeCAI",
    "MinimizeCpG",
    "MinimizeCodonPairBias",
    "MinimizeMRNADG",
    # Model
    "CSPModel",
    # Builder
    "build_csp_model",
]


# ==============================================================================
# Module-level named constants
# ==============================================================================

DEFAULT_GC_LO = 0.30
DEFAULT_GC_HI = 0.70
DEFAULT_CPG_WINDOW = 200
DEFAULT_CPG_THRESHOLD = 0.6
DEFAULT_MAX_T_RUN = 5
CAI_LOG_EPSILON = 1e-10
NEAREST_NEIGHBOR_GC = -1.5
NEAREST_NEIGHBOR_AU = -0.5
NEAREST_NEIGHBOR_GU = -0.3


# ==============================================================================
# Helper functions
# ==============================================================================

def codon_gc_count(codon: str) -> int:
    """Count G/C bases in a codon (0-3).

    Args:
        codon: A 3-character DNA string (e.g. ``"GGC"``).

    Returns:
        Number of bases that are G or C (0, 1, 2, or 3).

    Examples:
        >>> codon_gc_count("ATG")
        1
        >>> codon_gc_count("GGC")
        3
        >>> codon_gc_count("ATA")
        0
    """
    return sum(1 for b in codon if b in "GC")


def codon_contains_gt(codon: str) -> bool:
    """Check whether a codon contains the GT dinucleotide (5' splice donor).

    Only checks *within* the codon — cross-codon GT is handled separately
    at the constraint level.

    Args:
        codon: A 3-character DNA string.

    Returns:
        True if ``"GT"`` appears at positions 0-1 or 1-2 within the codon.

    Examples:
        >>> codon_contains_gt("GTT")
        True
        >>> codon_contains_gt("AGT")
        True
        >>> codon_contains_gt("GAT")
        False
    """
    return "GT" in codon


def codon_contains_ag(codon: str) -> bool:
    """Check whether a codon contains the AG dinucleotide (3' splice acceptor).

    Only checks *within* the codon — cross-codon AG is handled separately
    at the constraint level.

    Args:
        codon: A 3-character DNA string.

    Returns:
        True if ``"AG"`` appears at positions 0-1 or 1-2 within the codon.

    Examples:
        >>> codon_contains_ag("AGT")
        True
        >>> codon_contains_ag("CAG")
        True
        >>> codon_contains_ag("ACG")
        False
    """
    return "AG" in codon


def codon_contains_cpg(codon: str) -> bool:
    """Check whether a codon contains a CpG dinucleotide (CG) *within* the codon.

    Cross-codon CpG (C at position 2 of one codon, G at position 0 of the next)
    is not detected by this helper — it is handled at the constraint level by
    scanning the full sequence.

    Args:
        codon: A 3-character DNA string.

    Returns:
        True if ``"CG"`` appears at positions 0-1 or 1-2 within the codon.

    Examples:
        >>> codon_contains_cpg("CGT")
        True
        >>> codon_contains_cpg("ACG")
        True
        >>> codon_contains_cpg("CAG")
        False
    """
    return "CG" in codon


def compute_gc_from_codons(codons: Sequence[str]) -> float:
    """Compute GC content fraction from a sequence of codon strings.

    Args:
        codons: An iterable of 3-character DNA codon strings.

    Returns:
        GC fraction in [0.0, 1.0]. Returns 0.0 for empty input.

    Examples:
        >>> compute_gc_from_codons(["ATG", "GGC"])
        0.6666666666666666
        >>> compute_gc_from_codons(["AAA", "TTT"])
        0.0
    """
    seq = "".join(codons)
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return gc / len(seq)


def _sequence_from_codons(codons: Sequence[str]) -> str:
    """Concatenate codon strings into a full DNA sequence.

    Args:
        codons: Iterable of 3-character DNA codon strings.

    Returns:
        Concatenated DNA string.
    """
    return "".join(codons)


# ==============================================================================
# Abstract base classes for constraints
# ==============================================================================

class HardConstraint(ABC):
    """Abstract base for hard constraints (MUST hold).

    Every hard constraint must implement:
    - ``name``: unique string identifier
    - ``check(sequence)``: return True if the constraint holds
    - ``violated_positions(sequence)``: return positions where violated
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Unique identifier for this constraint."""
        ...

    @property
    def strictness(self) -> ConstraintStrictness:
        """Constraint classification (always HARD)."""
        return ConstraintStrictness.HARD

    @property
    def constraint_type(self) -> ConstraintType:
        """Constraint category — subclasses should override."""
        return ConstraintType.CUSTOM

    @abstractmethod
    def check(self, sequence: str) -> bool:
        """Return True if the constraint is satisfied for *sequence*.

        Args:
            sequence: Full DNA sequence to validate.

        Returns:
            True if the constraint holds, False otherwise.
        """
        ...

    @abstractmethod
    def violated_positions(self, sequence: str) -> list[int]:
        """Return nucleotide positions where the constraint is violated.

        Args:
            sequence: Full DNA sequence to validate.

        Returns:
            Sorted list of 0-based nucleotide positions involved in violations.
            Empty list if the constraint is satisfied.
        """
        ...

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name={self.name!r})"


class SoftConstraint(ABC):
    """Abstract base for soft constraints / objectives (PREFER to satisfy).

    Soft constraints contribute to the optimization objective.  They may be
    violated if they conflict with hard constraints.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Unique identifier for this soft constraint."""
        ...

    @property
    def strictness(self) -> ConstraintStrictness:
        """Constraint classification (always SOFT)."""
        return ConstraintStrictness.SOFT

    @property
    def constraint_type(self) -> ConstraintType:
        """Constraint category — subclasses should override."""
        return ConstraintType.CUSTOM

    @abstractmethod
    def check(self, sequence: str) -> bool:
        """Return True if the soft constraint is satisfied.

        Note: for optimization objectives like MaximizeCAI, this always
        returns True since there is no hard threshold — only a gradient
        to follow.
        """
        ...

    @abstractmethod
    def violated_positions(self, sequence: str) -> list[int]:
        """Return positions contributing to sub-optimality.

        For optimization objectives, this may return positions that could
        be improved (not necessarily "violated" in the hard sense).
        """
        ...

    @abstractmethod
    def score(self, sequence: str) -> float:
        """Compute the objective contribution (higher = better).

        Args:
            sequence: Full DNA sequence to score.

        Returns:
            Numerical score. For maximization objectives, higher is better.
            For minimization objectives, the convention is to return the
            *negated* cost so that the solver can always maximize.
        """
        ...

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name={self.name!r})"


# ==============================================================================
# Hard constraints
# ==============================================================================

class TranslationConstraint(HardConstraint):
    """Enforce that every codon translates to the correct amino acid.

    This constraint is implicitly enforced by the domain of each
    CodonVariable (only synonymous codons are in the domain), but
    we include it as an explicit constraint for:
    1. Completeness in the constraint model
    2. Verification of solutions (defense-in-depth)
    3. MUS diagnosis when domains are modified

    Attributes:
        protein: The target amino acid sequence.
    """

    def __init__(self, protein: str) -> None:
        self._protein = protein

    @property
    def name(self) -> str:
        return "TranslationConstraint"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.AMINO_ACID_IDENTITY

    @property
    def protein(self) -> str:
        return self._protein

    def check(self, sequence: str) -> bool:
        """Verify the sequence translates to the target protein.

        Checks that:
        1. len(sequence) == len(protein) * 3
        2. Every codon maps to the expected amino acid via CODON_TABLE
        """
        if len(sequence) != len(self._protein) * 3:
            return False
        for i, expected_aa in enumerate(self._protein):
            codon = sequence[i * 3 : i * 3 + 3]
            if CODON_TABLE.get(codon) != expected_aa:
                return False
        return True

    def violated_positions(self, sequence: str) -> list[int]:
        """Return codon start positions where translation mismatches."""
        violations: list[int] = []
        for i, expected_aa in enumerate(self._protein):
            codon_start = i * 3
            if codon_start + 3 > len(sequence):
                violations.append(codon_start)
                continue
            codon = sequence[codon_start : codon_start + 3]
            if CODON_TABLE.get(codon) != expected_aa:
                violations.append(codon_start)
        return violations


class NoRestrictionSiteConstraint(HardConstraint):
    """No restriction enzyme recognition site may appear in the sequence.

    Checks both the forward strand and reverse complement for each
    enzyme's recognition sequence.

    Attributes:
        sites: List of recognition sequences to avoid (e.g. ``["GAATTC"]``).
    """

    def __init__(self, sites: list[str]) -> None:
        self._sites = sites

    @property
    def name(self) -> str:
        return "NoRestrictionSiteConstraint"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.RESTRICTION_SITE

    @property
    def sites(self) -> list[str]:
        return list(self._sites)

    def check(self, sequence: str) -> bool:
        """Return True if no forbidden recognition site is present."""
        from ..constants import reverse_complement

        for site in self._sites:
            if site in sequence:
                return False
            rc = reverse_complement(site)
            if rc != site and rc in sequence:
                return False
        return True

    def violated_positions(self, sequence: str) -> list[int]:
        """Return all positions where a restriction site is found."""
        from ..constants import reverse_complement

        positions: list[int] = []
        for site in self._sites:
            pos = sequence.find(site)
            while pos != -1:
                positions.append(pos)
                pos = sequence.find(site, pos + 1)
            rc = reverse_complement(site)
            if rc != site:
                pos = sequence.find(rc)
                while pos != -1:
                    positions.append(pos)
                    pos = sequence.find(rc, pos + 1)
        return sorted(set(positions))


class GCRangeConstraint(HardConstraint):
    """Overall GC content must fall within [gc_lo, gc_hi].

    This is a global constraint on the entire sequence, not a per-window
    constraint.  The GC fraction is computed as the count of G and C bases
    divided by the total sequence length.

    Attributes:
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
    """

    def __init__(self, gc_lo: float = DEFAULT_GC_LO, gc_hi: float = DEFAULT_GC_HI) -> None:
        if not (0.0 <= gc_lo < gc_hi <= 1.0):
            raise ValueError(f"Invalid GC range: [{gc_lo}, {gc_hi}]")
        self._gc_lo = gc_lo
        self._gc_hi = gc_hi

    @property
    def name(self) -> str:
        return "GCRangeConstraint"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.GC_CONTENT

    @property
    def gc_lo(self) -> float:
        return self._gc_lo

    @property
    def gc_hi(self) -> float:
        return self._gc_hi

    @staticmethod
    def _gc_count(sequence: str) -> int:
        """Count G/C bases in a sequence."""
        return sum(1 for b in sequence if b in "GC")

    def _gc_fraction(self, sequence: str) -> tuple[int, float]:
        """Return (gc_count, gc_fraction) for the sequence.

        Computes the GC count and fraction once so that both ``check()``
        and ``violated_positions()`` can reuse the result without
        calling ``_gc_count`` twice.
        """
        gc = self._gc_count(sequence)
        return gc, gc / len(sequence)

    def check(self, sequence: str) -> bool:
        """Return True if GC content is within the allowed range."""
        if not sequence:
            return True
        _, frac = self._gc_fraction(sequence)
        return self._gc_lo <= frac <= self._gc_hi

    def violated_positions(self, sequence: str) -> list[int]:
        """Return positions of G/C bases if GC is out of range.

        If GC is too high, returns positions of all G/C bases.
        If GC is too low, returns positions of all A/T bases (the positions
        that *should* be G/C).
        If GC is in range, returns empty list.
        """
        if not sequence:
            return []
        gc, frac = self._gc_fraction(sequence)
        if self._gc_lo <= frac <= self._gc_hi:
            return []

        if frac > self._gc_hi:
            # Too much GC — report GC positions
            return [i for i, b in enumerate(sequence) if b in "GC"]
        else:
            # Too little GC — report AT positions (candidates for GC increase)
            return [i for i, b in enumerate(sequence) if b in "AT"]


class NoCrypticSpliceConstraint(HardConstraint):
    """No cryptic splice sites with MaxEnt score >= threshold.

    Scans every GT (donor) and AG (acceptor) dinucleotide in the
    sequence and checks that the MaxEntScan score is below the
    threshold.  Both within-codon and cross-codon dinucleotides
    are considered.

    Amino acids like Valine (V) have only GT-containing codons
    (GTT, GTC, GTA, GTG), making GT avoidance at those positions
    impossible.  When a ``protein`` sequence is provided, the
    constraint automatically skips GT dinucleotides that are within
    Valine codons — these are *unavoidable* and should not be
    reported as violations.

    Attributes:
        threshold: Maximum allowed MaxEntScan score. Sites scoring
            at or above this threshold are considered cryptic and
            violate the constraint.
        protein: The amino acid sequence.  When provided, Valine
            codon positions are excluded from GT scanning since
            GT within Valine codons is unavoidable.
    """

    # Amino acids whose codons ALL contain the GT dinucleotide
    _GT_UNAVOIDABLE_AAS: ClassVar[set[str]] = {"V"}

    def __init__(self, threshold: float = 3.0, protein: str = "") -> None:
        if threshold <= 0:
            raise ValueError(f"Threshold must be positive, got {threshold}")
        self._threshold = threshold
        self._protein = protein.upper() if protein else ""

    @property
    def name(self) -> str:
        return "NoCrypticSpliceConstraint"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.NO_CRYPTIC_SPLICE

    @property
    def threshold(self) -> float:
        return self._threshold

    @property
    def protein(self) -> str:
        return self._protein

    def _unavoidable_gt_positions(self, sequence: str) -> set[int]:
        """Compute nucleotide positions where GT is unavoidable (Valine codons).

        For each Valine position in the protein, the GT dinucleotide
        occupies nucleotide positions [i*3, i*3+1] within the codon.
        These positions cannot be avoided and should be excluded from
        violation reporting.

        Args:
            sequence: DNA sequence (used only for length validation).

        Returns:
            Set of nucleotide positions where GT is within an unavoidable
            Valine codon.
        """
        positions: set[int] = set()
        if not self._protein:
            return positions
        for i, aa in enumerate(self._protein):
            if aa in self._GT_UNAVOIDABLE_AAS:
                # Valine codons are always GTN, so GT is at positions
                # [i*3, i*3+1] within the codon
                nuc_start = i * 3
                positions.add(nuc_start)  # Position of the G in GT
        return positions

    def _scan_splice_sites(self, sequence: str, threshold: float) -> list[int]:
        """Scan for cryptic splice sites at or above *threshold*.

        Returns a list of positions (0-based) where GT donor or AG acceptor
        dinucleotides have MaxEntScan scores >= threshold.

        GT dinucleotides at Valine codon positions are excluded from
        scanning since they are unavoidable.

        Note: The canonical sliding-window GT/AG scanning implementation lives
        in :func:`biocompiler.maxentscan.scan_splice_sites`, which returns
        richer ``(position, site_type, score)`` tuples.  This method is a
        lightweight wrapper used by the constraint model.
        """
        sequence = sequence.upper()
        positions: list[int] = []
        unavoidable_gt = self._unavoidable_gt_positions(sequence)

        for i in range(len(sequence) - 1):
            if sequence[i : i + 2] == "GT":
                # Skip GT within unavoidable Valine codons
                if i in unavoidable_gt:
                    continue
                s = score_donor(sequence, i)
                if s >= threshold:
                    positions.append(i)

        for i in range(len(sequence) - 1):
            if sequence[i : i + 2] == "AG":
                s = score_acceptor(sequence, i)
                if s >= threshold:
                    positions.append(i)

        return positions

    def check(self, sequence: str) -> bool:
        """Return True if no cryptic splice site exceeds the threshold."""
        return not self._scan_splice_sites(sequence, self._threshold)

    def violated_positions(self, sequence: str) -> list[int]:
        """Return positions of cryptic splice sites above threshold."""
        return sorted(set(self._scan_splice_sites(sequence, self._threshold)))


class NoCpGIslandConstraint(HardConstraint):
    """No CpG island (Obs/Exp CG ratio > threshold in any sliding window).

    A CpG island is detected when the observed-to-expected ratio of CG
    dinucleotides in any sliding window exceeds the threshold.  The
    expected CG count is ``C_count * G_count / window_length``.

    CpG island avoidance is primarily relevant for mammalian expression
    systems where CpG methylation can lead to gene silencing.  For
    prokaryotic organisms, this constraint is automatically satisfied
    (check returns True, violated_positions returns []) since CpG
    islands have no known regulatory significance in prokaryotes.

    Attributes:
        window: Window size in nucleotides (default 200).
        threshold: Maximum allowed Obs/Exp CG ratio (default 0.6).
        organism: Target organism name.  If prokaryotic, the constraint
            is automatically satisfied.
    """

    def __init__(self, window: int = DEFAULT_CPG_WINDOW, threshold: float = DEFAULT_CPG_THRESHOLD, organism: str = "") -> None:
        if window <= 0:
            raise ValueError(f"Window must be positive, got {window}")
        if threshold <= 0:
            raise ValueError(f"Threshold must be positive, got {threshold}")
        self._window = window
        self._threshold = threshold
        self._organism = organism

    @property
    def name(self) -> str:
        return "NoCpGIslandConstraint"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.NO_CPG

    @property
    def window(self) -> int:
        return self._window

    @property
    def threshold(self) -> float:
        return self._threshold

    @property
    def organism(self) -> str:
        return self._organism

    def _is_prokaryotic(self) -> bool:
        """Return True if the configured organism is prokaryotic."""
        if not self._organism:
            return False
        return not is_eukaryotic_organism(self._organism)

    def _scan_cpg(self, sequence: str, window: int, threshold: float) -> list[int]:
        """Scan for CpG islands in sliding windows.

        Returns start positions of windows where the Obs/Exp CG ratio
        exceeds *threshold*.

        Note: The canonical CpG island sliding-window implementation lives
        in :func:`biocompiler.type_system.check_no_cpg_island` and
        :func:`biocompiler.type_system.evaluate_no_cpg_island`.  This method
        is a lightweight variant used by the constraint model.
        """
        sequence = sequence.upper()
        violating_starts: list[int] = []

        if len(sequence) < window:
            return violating_starts

        for start in range(len(sequence) - window + 1):
            window_seq = sequence[start : start + window]
            c_count = window_seq.count("C")
            g_count = window_seq.count("G")
            cg_count = sum(
                1 for i in range(len(window_seq) - 1) if window_seq[i : i + 2] == "CG"
            )
            expected = (c_count * g_count) / window if window > 0 else 0
            if expected > 0:
                obs_exp = cg_count / expected
                if obs_exp > threshold:
                    violating_starts.append(start)

        return violating_starts

    def check(self, sequence: str) -> bool:
        """Return True if no sliding window has CpG Obs/Exp above threshold.

        Automatically returns True for prokaryotic organisms.
        """
        if self._is_prokaryotic():
            logger.info(
                "NoCpGIslandConstraint.check() skipped for prokaryotic organism '%s'",
                self._organism,
            )
            return True
        return not self._scan_cpg(sequence, self._window, self._threshold)

    def violated_positions(self, sequence: str) -> list[int]:
        """Return start positions of windows with CpG island violations.

        Automatically returns [] for prokaryotic organisms.
        """
        if self._is_prokaryotic():
            return []
        return self._scan_cpg(sequence, self._window, self._threshold)


class NoATTTAMotifConstraint(HardConstraint):
    """No ATTTA instability motifs may appear in the sequence.

    ATTTA is an AU-rich element (ARE) that targets mRNA for rapid
    degradation.  Removing these motifs is critical for stable
    transgene expression in mammalian cells.

    This constraint scans for the literal substring ``"ATTTA"`` in
    the full sequence, including cross-codon occurrences.
    """

    @property
    def name(self) -> str:
        return "NoATTTAMotifConstraint"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.NO_INSTABILITY_MOTIF

    def check(self, sequence: str) -> bool:
        """Return True if no ATTTA motif is present."""
        return INSTABILITY_MOTIF not in sequence.upper()

    def violated_positions(self, sequence: str) -> list[int]:
        """Return start positions of all ATTTA occurrences."""
        sequence = sequence.upper()
        positions: list[int] = []
        pos = sequence.find(INSTABILITY_MOTIF)
        while pos != -1:
            positions.append(pos)
            pos = sequence.find(INSTABILITY_MOTIF, pos + 1)
        return positions


class NoTRunConstraint(HardConstraint):
    """No run of 6 or more consecutive T bases.

    Long poly-T runs can act as premature transcription termination
    signals (Pol II pause sites) in eukaryotic expression systems.
    This constraint flags any stretch of 6+ consecutive T characters.

    Attributes:
        max_run: Maximum allowed consecutive T count (default 5, so
            runs of 6+ are forbidden).
    """

    def __init__(self, max_run: int = DEFAULT_MAX_T_RUN) -> None:
        if max_run < 1:
            raise ValueError(f"max_run must be >= 1, got {max_run}")
        self._max_run = max_run

    @property
    def name(self) -> str:
        return "NoTRunConstraint"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.MRNA_STABILITY

    @property
    def max_run(self) -> int:
        return self._max_run

    def check(self, sequence: str) -> bool:
        """Return True if no T-run exceeds max_run."""
        sequence = sequence.upper()
        run_length = 0
        for base in sequence:
            if base == "T":
                run_length += 1
                if run_length > self._max_run:
                    return False
            else:
                run_length = 0
        return True

    def violated_positions(self, sequence: str) -> list[int]:
        """Return start positions of T-runs exceeding max_run."""
        sequence = sequence.upper()
        positions: list[int] = []
        i = 0
        while i < len(sequence):
            if sequence[i] == "T":
                run_start = i
                while i < len(sequence) and sequence[i] == "T":
                    i += 1
                run_length = i - run_start
                if run_length > self._max_run:
                    positions.append(run_start)
            else:
                i += 1
        return positions


# ==============================================================================
# Soft constraints / objectives
# ==============================================================================

class MaximizeCAI(SoftConstraint):
    """Maximize the Codon Adaptation Index (CAI) across all positions.

    CAI is the geometric mean of relative adaptiveness values:

        CAI = (prod_i w_i) ^ (1/N)

    where w_i is the relative adaptiveness of the codon chosen at
    position i, and N is the number of codons.

    For the CSP solver, we maximize the log-CAI (equivalent, but
    avoids numerical underflow):

        log(CAI) = (1/N) * sum_i log(w_i)

    Attributes:
        adaptiveness: Dict mapping codon strings to relative adaptiveness
            values (0.0-1.0) for the target organism.
        protein: Target amino acid sequence (needed to identify which
            codons are relevant at each position).
    """

    def __init__(self, adaptiveness: dict[str, float], protein: str) -> None:
        self._adaptiveness = adaptiveness
        self._protein = protein

    @property
    def name(self) -> str:
        return "MaximizeCAI"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.CODON_USAGE

    @property
    def adaptiveness(self) -> dict[str, float]:
        return dict(self._adaptiveness)

    def check(self, sequence: str) -> bool:
        """Always True — CAI is an optimization objective, not a threshold."""
        return True

    def violated_positions(self, sequence: str) -> list[int]:
        """Return positions where CAI could be improved (non-optimal codons).

        A codon is "sub-optimal" if another synonymous codon has higher
        adaptiveness.  These positions are candidates for improvement.
        """
        positions: list[int] = []
        for i, aa in enumerate(self._protein):
            codon = sequence[i * 3 : i * 3 + 3]
            current_w = self._adaptiveness.get(codon, 0.0)
            best_w = max(
                self._adaptiveness.get(c, 0.0) for c in AA_TO_CODONS.get(aa, [codon])
            )
            if current_w < best_w:
                positions.append(i * 3)
        return positions

    def score(self, sequence: str) -> float:
        """Compute log-CAI (higher = better).

        Returns:
            The sum of log(adaptiveness) across all codon positions.
            If any adaptiveness is 0, uses a small epsilon to avoid -inf.
        """
        total = 0.0
        for i in range(len(self._protein)):
            codon = sequence[i * 3 : i * 3 + 3]
            w = self._adaptiveness.get(codon, CAI_LOG_EPSILON)
            total += math.log(max(w, CAI_LOG_EPSILON))
        return total

    def cai(self, sequence: str) -> float:
        """Compute the actual CAI value (geometric mean of adaptiveness).

        Returns:
            CAI in [0.0, 1.0]. Returns 0.0 if any adaptiveness is zero.
        """
        n = len(self._protein)
        if n == 0:
            return 0.0
        log_product = self.score(sequence)
        return math.exp(log_product / n)


class MinimizeCpG(SoftConstraint):
    """Minimize the number of CpG (CG) dinucleotides in the sequence.

    CpG dinucleotides are methylation targets in mammalian cells and can
    lead to gene silencing.  While complete elimination may be infeasible
    (some amino acids like Arginine require CG-containing codons), we
    minimize the total count.

    For prokaryotic organisms, CpG minimization is irrelevant (no DNA
    methylation silencing mechanism), so the score always returns 0.0
    and violated_positions returns [].

    Both within-codon and cross-codon CG dinucleotides are counted
    (for eukaryotic organisms).

    Attributes:
        organism: Target organism name.  If prokaryotic, the soft
            constraint becomes a no-op (score always 0.0).
    """

    def __init__(self, organism: str = "") -> None:
        self._organism = organism

    @property
    def name(self) -> str:
        return "MinimizeCpG"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.NO_CPG

    @property
    def organism(self) -> str:
        return self._organism

    def _is_prokaryotic(self) -> bool:
        """Return True if the configured organism is prokaryotic."""
        if not self._organism:
            return False
        return not is_eukaryotic_organism(self._organism)

    def check(self, sequence: str) -> bool:
        """Always True — CpG count is an optimization objective."""
        return True

    def violated_positions(self, sequence: str) -> list[int]:
        """Return positions of all CG dinucleotides (candidates for removal).

        Returns [] for prokaryotic organisms since CpG is irrelevant.
        """
        if self._is_prokaryotic():
            return []
        sequence = sequence.upper()
        return [i for i in range(len(sequence) - 1) if sequence[i : i + 2] == "CG"]

    def score(self, sequence: str) -> float:
        """Return negated CG count (higher = fewer CG = better).

        The solver maximizes the objective, so we negate the count.
        For prokaryotic organisms, returns 0.0 (neutral objective).
        """
        if self._is_prokaryotic():
            return 0.0
        sequence = sequence.upper()
        cpg_count = sum(
            1 for i in range(len(sequence) - 1) if sequence[i : i + 2] == "CG"
        )
        return -float(cpg_count)

    def cpg_count(self, sequence: str) -> int:
        """Return the raw number of CG dinucleotides.

        Returns 0 for prokaryotic organisms since CpG is irrelevant.
        """
        if self._is_prokaryotic():
            return 0
        sequence = sequence.upper()
        return sum(
            1 for i in range(len(sequence) - 1) if sequence[i : i + 2] == "CG"
        )


class MinimizeCodonPairBias(SoftConstraint):
    """Maximize codon pair bias (CPB) across consecutive codon pairs.

    Codon pair bias measures the over/under-representation of consecutive
    codon pairs relative to expected frequency.  Over-represented pairs
    (positive CPB) are favoured for expression; under-represented pairs
    (negative CPB) are disfavoured.

    The score is the arithmetic mean of per-pair CPB scores across the
    sequence.  Higher values indicate the sequence uses over-represented
    (favoured) codon pairs.

    For the CSP solver, we return the mean CPB directly (higher = better).

    Attributes:
        organism: Target organism identifier (e.g. ``"Escherichia_coli"``).
        cpb_data: Dict mapping ``"{codon1}-{codon2}"`` to CPB score.
    """

    def __init__(self, organism: str) -> None:
        from ..codon_pair_scoring import get_codon_pair_data
        self._organism = organism
        self._cpb_data = get_codon_pair_data(organism)

    @property
    def name(self) -> str:
        return "MinimizeCodonPairBias"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.CODON_PAIR_BIAS

    @property
    def organism(self) -> str:
        return self._organism

    @property
    def cpb_data(self) -> dict[str, float]:
        return dict(self._cpb_data)

    def check(self, sequence: str) -> bool:
        """Always True -- CPB is an optimization objective, not a threshold."""
        return True

    def violated_positions(self, sequence: str) -> list[int]:
        """Return start positions of codon pairs with negative CPB.

        Codon pairs with negative CPB are candidates for improvement
        (substitution with synonymous codons that form a better pair).
        """
        sequence = sequence.upper()
        if len(sequence) < 6:
            return []
        positions: list[int] = []
        for i in range(0, len(sequence) - 3, 3):
            codon1 = sequence[i:i + 3]
            codon2 = sequence[i + 3:i + 6]
            pair_key = f"{codon1}-{codon2}"
            if self._cpb_data.get(pair_key, 0.0) < 0.0:
                positions.append(i)
        return positions

    def score(self, sequence: str) -> float:
        """Compute mean codon pair bias score (higher = better).

        Returns:
            The arithmetic mean of CPB scores for all consecutive
            codon pairs.  Pairs not in the CPB table default to 0.0.
            Returns 0.0 for sequences shorter than two codons.
        """
        sequence = sequence.upper()
        if len(sequence) < 6 or not self._cpb_data:
            return 0.0

        codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
        scores: list[float] = []
        for i in range(len(codons) - 1):
            pair_key = f"{codons[i]}-{codons[i + 1]}"
            scores.append(self._cpb_data.get(pair_key, 0.0))

        return sum(scores) / len(scores) if scores else 0.0

    def geometric_mean_score(self, sequence: str) -> float:
        """Compute geometric mean of 2^CPB across all codon pairs.

        For each codon pair, the contribution is 2^CPB (which equals
        observed/expected).  The geometric mean of these ratios indicates
        the overall fold-enrichment of codon pairs in the sequence.

        Returns:
            Geometric mean of 2^CPB values.  Returns 1.0 for sequences
            shorter than two codons or if no CPB data is available.
        """
        sequence = sequence.upper()
        if len(sequence) < 6 or not self._cpb_data:
            return 1.0

        codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
        log_values: list[float] = []
        for i in range(len(codons) - 1):
            pair_key = f"{codons[i]}-{codons[i + 1]}"
            cpb = self._cpb_data.get(pair_key, 0.0)
            # 2^CPB = observed/expected; log2(obs/exp) = CPB
            log_values.append(cpb * math.log(2))  # ln(2^CPB) = CPB * ln(2)

        if not log_values:
            return 1.0
        mean_log = sum(log_values) / len(log_values)
        return math.exp(mean_log)


class MinimizeMRNADG(SoftConstraint):
    """Minimize mRNA folding free energy (dG) at the 5' end.

    Strong secondary structure near the 5' end (around the RBS/start codon)
    can block ribosome binding and reduce translation efficiency.  This
    objective prefers sequences with weaker 5' structure (less negative dG).

    Uses ViennaRNA via ``compute_5prime_dg()`` when available for accurate
    thermodynamic folding.  Falls back to a simplified nearest-neighbor
    approximation when ViennaRNA is not installed.

    Attributes:
        window_start: Start of the analysis window (default 0).
        window_end: End of the analysis window (default 50).
    """

    def __init__(
        self,
        window_start: int = 0,
        window_end: int = 50,
    ) -> None:
        self._window_start = window_start
        self._window_end = window_end

    @property
    def name(self) -> str:
        return "MinimizeMRNADG"

    @property
    def constraint_type(self) -> ConstraintType:
        return ConstraintType.MRNA_STABILITY

    @property
    def window_start(self) -> int:
        return self._window_start

    @property
    def window_end(self) -> int:
        return self._window_end

    def check(self, sequence: str) -> bool:
        """Always True — mRNA dG is an optimization objective."""
        return True

    def violated_positions(self, sequence: str) -> list[int]:
        """Return positions in the 5' window that contribute to stable structure.

        Positions with G/C bases in the analysis window are reported as
        candidates for A/T substitution to reduce structure stability.
        """
        effective_end = min(self._window_end, len(sequence))
        window_seq = sequence[self._window_start : effective_end].upper()
        return [
            self._window_start + i
            for i, b in enumerate(window_seq)
            if b in "GC"
        ]

    def score(self, sequence: str) -> float:
        """Return negated |dG| approximation (higher = weaker structure = better).

        Uses ViennaRNA when available for accurate MFE computation.
        Falls back to a simplified nearest-neighbor dG estimate based on
        potential base pairs in the 5' window when ViennaRNA is not installed.
        """
        dg = self.compute_dg(sequence)
        # We want to MAXIMIZE, so return the negated |dG|
        # (less stable structure = higher score = better)
        return -abs(dg)

    def compute_dg(self, sequence: str) -> float:
        """Return the dG value (more negative = more stable structure).

        Uses ViennaRNA via ``compute_5prime_dg`` when available for accurate
        MFE computation.  Falls back to a simplified nearest-neighbor
        approximation when ViennaRNA is not installed.
        """
        # --- Try ViennaRNA first ---
        try:
            from ..viennarna import compute_5prime_dg, is_viennarna_available
            if is_viennarna_available():
                window_len = self._window_end - self._window_start
                # compute_5prime_dg folds the first `window` nt from the
                # given sequence; we slice our window out first.
                window_seq = sequence[self._window_start : self._window_end]
                return compute_5prime_dg(window_seq, window=window_len)
        except ImportError:
            logger.warning("ViennaRNA not installed, using fallback dG approximation", exc_info=True)
        except Exception:
            logger.warning("ViennaRNA computation failed, using fallback dG approximation", exc_info=True)

        # --- Fallback: simplified nearest-neighbor approximation ---
        sequence = sequence.upper()
        effective_end = min(self._window_end, len(sequence))
        window_seq = sequence[self._window_start : effective_end]

        if len(window_seq) < 4:
            return 0.0

        rna = window_seq.replace("T", "U")
        half = len(rna) // 2
        first_half = rna[:half]
        second_half = rna[half : 2 * half]

        gc_pairs = 0
        au_pairs = 0
        gu_pairs = 0

        for i in range(min(len(first_half), len(second_half))):
            j = len(second_half) - 1 - i
            if j < 0:
                break
            base_5 = first_half[i]
            base_3 = second_half[j]
            if (base_5 == "G" and base_3 == "C") or (base_5 == "C" and base_3 == "G"):
                gc_pairs += 1
            elif (base_5 == "A" and base_3 == "U") or (base_5 == "U" and base_3 == "A"):
                au_pairs += 1
            elif (base_5 == "G" and base_3 == "U") or (base_5 == "U" and base_3 == "G"):
                gu_pairs += 1

        return NEAREST_NEIGHBOR_GC * gc_pairs + NEAREST_NEIGHBOR_AU * au_pairs + NEAREST_NEIGHBOR_GU * gu_pairs


# ==============================================================================
# Helper: extract params from HardConstraint / SoftConstraint instances
# ==============================================================================

def _extract_constraint_params(constraint: HardConstraint | SoftConstraint) -> dict[str, Any]:
    """Extract check-relevant parameters from a constraint object.

    Returns a dict of parameters that :meth:`ConstraintSpec.check` can use
    to verify the constraint against a sequence without needing the original
    ``HardConstraint`` / ``SoftConstraint`` instance.
    """
    from .types import ConstraintType as CT

    ctype = constraint.constraint_type

    if ctype == CT.GC_CONTENT:
        # GCRangeConstraint has gc_lo / gc_hi properties
        return {
            "gc_lo": getattr(constraint, "gc_lo", 0.30),
            "gc_hi": getattr(constraint, "gc_hi", 0.70),
        }

    if ctype == CT.NO_CPG:
        # NoCpGIslandConstraint has window / threshold
        return {
            "window": getattr(constraint, "window", 200),
            "threshold": getattr(constraint, "threshold", 0.6),
        }

    if ctype in (CT.NO_CRYPTIC_SPLICE, CT.SPLICE_DONOR_AVOIDANCE):
        # NoCrypticSpliceConstraint has threshold and protein
        return {
            "threshold": getattr(constraint, "threshold", 3.0),
            "protein": getattr(constraint, "protein", ""),
            "organism": "",
        }

    if ctype == CT.RESTRICTION_SITE:
        # NoRestrictionSiteConstraint has sites
        return {
            "sites": getattr(constraint, "sites", []),
        }

    if ctype == CT.NO_INSTABILITY_MOTIF:
        # NoATTTAMotifConstraint uses constant INSTABILITY_MOTIF
        return {
            "motif": INSTABILITY_MOTIF,
        }

    if ctype == CT.MRNA_STABILITY:
        # NoTRunConstraint has max_run
        return {
            "max_run": getattr(constraint, "max_run", 5),
        }

    if ctype == CT.AMINO_ACID_IDENTITY:
        # TranslationConstraint has protein
        return {
            "protein": getattr(constraint, "protein", ""),
        }

    # For soft constraints (CODON_USAGE, etc.) and unknown types,
    # no params needed — ConstraintSpec.check() defaults to True.
    return {}


# Default enforcement priority for each well-known constraint type.
_CONSTRAINT_PRIORITY_MAP: dict[str, ConstraintPriority] = {
    "TranslationConstraint": ConstraintPriority.CRITICAL,
    "NoRestrictionSiteConstraint": ConstraintPriority.HIGH,
    "GCRangeConstraint": ConstraintPriority.MEDIUM,
    "NoCrypticSpliceConstraint": ConstraintPriority.MEDIUM,
    "NoCpGIslandConstraint": ConstraintPriority.MEDIUM,
    "NoATTTAMotifConstraint": ConstraintPriority.MEDIUM,
    "NoTRunConstraint": ConstraintPriority.LOW,
    "MaximizeCAI": ConstraintPriority.LOW,
    "MinimizeCpG": ConstraintPriority.LOW,
    "MinimizeCodonPairBias": ConstraintPriority.LOW,
    "MinimizeMRNADG": ConstraintPriority.LOW,
}


def _default_priority_for(name: str) -> ConstraintPriority:
    """Return the default enforcement priority for a constraint by name.

    Args:
        name: The constraint's ``name`` property (e.g. ``"GCRangeConstraint"``).

    Returns:
        The corresponding :class:`ConstraintPriority`, or MEDIUM if unknown.
    """
    return _CONSTRAINT_PRIORITY_MAP.get(name, ConstraintPriority.MEDIUM)


def _default_weight_for(name: str, config: SolverConfig) -> float:
    """Return the default scoring weight for a soft constraint by name.

    Uses the solver configuration's objective weights for the standard
    soft constraints (CAI, CpG, mRNA dG).  Returns 1.0 for unknown
    constraints.

    Args:
        name: The constraint's ``name`` property.
        config: Solver configuration providing objective weights.

    Returns:
        Weight multiplier for soft constraint scoring.
    """
    weight_map: dict[str, float] = {
        "MaximizeCAI": config.cai_weight,
        "MinimizeCpG": config.cpg_weight,
        "MinimizeMRNADG": config.mrna_dg_weight,
        "MinimizeCodonPairBias": config.codon_pair_bias_weight,
    }
    return weight_map.get(name, 1.0)


# ==============================================================================
# CSPModel dataclass
# ==============================================================================

@dataclass
class CSPModel:
    """Complete constraint satisfaction model for gene optimization.

    Encapsulates all decision variables, constraints, and configuration
    needed by the solver backends (OR-Tools, Z3, greedy).

    Construction is done via the factory function :func:`build_csp_model`,
    which creates variables from the protein sequence and registers all
    constraints based on the solver configuration.

    Attributes:
        variables: One CodonVariable per amino acid position.
        hard_constraints: Constraints that MUST be satisfied.
        soft_constraints: Optimization objectives (prefer to satisfy).
        protein: Target amino acid sequence.
        organism: Target organism name (e.g. ``"Homo_sapiens"``).
        config: Solver configuration (GC bounds, thresholds, weights).
    """

    variables: list[CodonVariable]
    hard_constraints: list[HardConstraint]
    soft_constraints: list[SoftConstraint]
    protein: str
    organism: str
    config: SolverConfig

    @property
    def num_variables(self) -> int:
        """Number of decision variables (codon positions)."""
        return len(self.variables)

    @property
    def num_hard_constraints(self) -> int:
        """Number of hard constraints."""
        return len(self.hard_constraints)

    @property
    def num_soft_constraints(self) -> int:
        """Number of soft constraints / objectives."""
        return len(self.soft_constraints)

    @property
    def protein_sequence(self) -> str:
        """Alias for :attr:`protein` — compatibility with ``types.CSPModel``."""
        return self.protein

    @property
    def codon_domains(self) -> dict[int, list[str]]:
        """Mapping from codon position to allowed codon domain.

        Derived from :attr:`variables`.  Provides compatibility with
        ``types.CSPModel`` which stores domains in this dict format.
        """
        return {var.position: list(var.domain) for var in self.variables}

    @property
    def constraints(self) -> list[ConstraintSpec]:
        """Constraint specs derived from hard and soft constraints.

        Converts the typed constraint objects (:class:`HardConstraint` /
        :class:`SoftConstraint`) into :class:`ConstraintSpec` instances so
        that code expecting ``types.CSPModel.constraints`` (a list of
        ``ConstraintSpec``) works transparently.

        Each spec carries an enforcement ``priority`` (CRITICAL for
        translation, HIGH for restriction sites, MEDIUM for most biological
        constraints, LOW for advisory/soft constraints) and a ``weight``
        derived from the solver configuration's objective weights.
        """
        specs: list[ConstraintSpec] = []
        for hc in self.hard_constraints:
            specs.append(ConstraintSpec(
                ctype=hc.constraint_type,
                name=hc.name,
                strictness=ConstraintStrictness.HARD,
                params=_extract_constraint_params(hc),
                priority=_default_priority_for(hc.name),
            ))
        for sc in self.soft_constraints:
            specs.append(ConstraintSpec(
                ctype=sc.constraint_type,
                name=sc.name,
                strictness=ConstraintStrictness.SOFT,
                params=_extract_constraint_params(sc),
                priority=_default_priority_for(sc.name),
                weight=_default_weight_for(sc.name, self.config),
            ))
        return specs

    def sequence_from_assignment(self, assignment: dict[int, str] | None = None) -> str:
        """Build a DNA sequence from variable assignments.

        Args:
            assignment: Dict mapping codon position index to codon string.
                If None, uses ``current_value`` from each variable.

        Returns:
            Concatenated DNA sequence.

        Raises:
            ValueError: If any variable has no assigned value.
        """
        codons: list[str] = []
        for var in self.variables:
            if assignment is not None:
                codon = assignment.get(var.position)
            else:
                codon = var.current_value
            if codon is None:
                raise ValueError(
                    f"CodonVariable at position {var.position} has no assignment"
                )
            codons.append(codon)
        return _sequence_from_codons(codons)

    def check_all_hard(self, sequence: str) -> bool:
        """Check whether all hard constraints are satisfied.

        Args:
            sequence: DNA sequence to validate.

        Returns:
            True if every hard constraint holds.
        """
        return all(c.check(sequence) for c in self.hard_constraints)

    def check_all_soft(self, sequence: str) -> dict[str, bool]:
        """Check each soft constraint.

        Args:
            sequence: DNA sequence to evaluate.

        Returns:
            Dict mapping constraint name to satisfaction status.
        """
        return {c.name: c.check(sequence) for c in self.soft_constraints}

    def objective_value(self, sequence: str) -> float:
        """Compute the weighted sum of all soft constraint scores.

        Uses the weights from the SolverConfig to combine individual
        soft constraint scores into a single objective value.

        Args:
            sequence: DNA sequence to score.

        Returns:
            Weighted objective value (higher = better).
        """
        total = 0.0
        weight_map: dict[str, float] = {
            "MaximizeCAI": self.config.cai_weight,
            "MinimizeCpG": self.config.cpg_weight,
            "MinimizeMRNADG": self.config.mrna_dg_weight,
            "MinimizeCodonPairBias": self.config.codon_pair_bias_weight,
        }
        for sc in self.soft_constraints:
            w = weight_map.get(sc.name, 0.0)
            total += w * sc.score(sequence)
        return total

    def hard_violations(self, sequence: str) -> dict[str, list[int]]:
        """Return all hard constraint violations.

        Args:
            sequence: DNA sequence to validate.

        Returns:
            Dict mapping constraint name to list of violated positions.
            Only constraints with violations are included.
        """
        result: dict[str, list[int]] = {}
        for c in self.hard_constraints:
            vp = c.violated_positions(sequence)
            if vp:
                result[c.name] = vp
        return result


# ==============================================================================
# Model builder
# ==============================================================================

def build_csp_model(
    protein: str,
    organism: str | None = None,
    config: SolverConfig | None = None,
) -> CSPModel:
    """Build a complete CSP model for gene optimization.

    This is the main factory function for the constraint model.  It:

    1. Creates one ``CodonVariable`` per amino acid position, with domain
       set to the synonymous codons for that amino acid (ordered by CAI).
    2. Registers all applicable hard constraints based on the config.
    3. Registers all soft constraints (optimization objectives).
    4. Returns the complete model ready for solving.

    Args:
        protein: Target amino acid sequence (e.g. ``"MVSKGE"``).
        organism: Target organism name.  If ``None``, falls back to
            ``config.organism``.  When both are provided, the explicit
            *organism* parameter takes precedence.
        config: Solver configuration. Uses defaults if not provided.

    Returns:
        A fully constructed :class:`CSPModel` ready for the solver.

    Raises:
        ValueError: If the protein contains invalid amino acid codes.
        ValueError: If the organism is not supported.
    """
    if config is None:
        config = SolverConfig()

    # Resolve organism: explicit parameter > config.organism > default
    if organism is None:
        organism = config.organism

    protein = protein.upper().strip()
    _validate_protein(protein)
    _validate_organism(organism)

    # Get organism-specific codon adaptiveness based on reference set
    if config.cai_reference_set == "sharp_li":
        sharp_li_tables = get_sharp_li_adaptiveness_tables()
        adaptiveness = sharp_li_tables.get(organism, {})
        logger.info("Using Sharp-Li (1987) CAI reference set for %s", organism)
    else:
        adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})
        logger.info("Using Kazusa CAI reference set for %s", organism)
    if not adaptiveness:
        logger.warning("No codon adaptiveness table for %s; using uniform weights", organism)
        adaptiveness = {codon: 1.0 for codon in CODON_TABLE}

    # --- Step 1: Create variables ---
    variables: list[CodonVariable] = []
    for i, aa in enumerate(protein):
        codons = AA_TO_CODONS.get(aa, [])
        if not codons:
            raise ValueError(f"No codons found for amino acid '{aa}' at position {i}")
        # Sort domain by adaptiveness (highest first) for solver guidance
        sorted_codons = sorted(
            codons,
            key=lambda c: adaptiveness.get(c, 0.0),
            reverse=True,
        )
        variables.append(CodonVariable(
            position=i,
            amino_acid=aa,
            domain=sorted_codons,
        ))

    # --- Determine whether to apply eukaryote-specific constraints ---
    # When auto_detect_organism_domain is True (default), the solver
    # automatically skips constraints that are only relevant for eukaryotes
    # (splice sites, CpG islands) when the target organism is prokaryotic.
    # This recovers ~0.27 CAI on prokaryotic targets by removing constraints
    # that unnecessarily restrict the codon search space.
    is_eukaryote: bool = True  # safe default
    if config.auto_detect_organism_domain:
        is_eukaryote = is_eukaryotic_organism(organism)
        if not is_eukaryote:
            logger.info(
                "Organism %s is prokaryotic; will skip eukaryote-only constraints",
                organism,
            )

    # --- Step 2: Register hard constraints ---
    hard_constraints: list[HardConstraint] = []

    # 2a. Translation (always present, domain-enforced but explicit)
    hard_constraints.append(TranslationConstraint(protein))

    # 2b. Restriction sites (from config or all known enzymes)
    sites = config.restriction_sites
    if not sites:
        # Default: avoid common 6+ bp restriction enzyme sites.
        # 4 bp sites (AluI, HaeIII, MboI, etc.) are excluded because
        # they appear too frequently to avoid in long sequences and
        # would make the solver model infeasible.
        sites = [
            seq for seq in RESTRICTION_ENZYMES.values()
            if all(b in "ACGT" for b in seq.upper()) and len(seq) >= 6
        ]
    # Filter to concrete sites (no IUPAC ambiguity codes)
    concrete_sites = [s for s in sites if all(b in "ACGT" for b in s.upper())]
    if concrete_sites:
        hard_constraints.append(NoRestrictionSiteConstraint(concrete_sites))

    # 2c. GC range
    hard_constraints.append(GCRangeConstraint(gc_lo=config.gc_lo, gc_hi=config.gc_hi))

    # 2d. Cryptic splice sites (eukaryote-only)
    if is_eukaryote:
        hard_constraints.append(
            NoCrypticSpliceConstraint(threshold=config.cryptic_splice_threshold, protein=protein)
        )
    else:
        logger.info(
            "Skipping eukaryote-only constraint NoCrypticSpliceConstraint "
            "for prokaryotic organism %s",
            organism,
        )

    # 2e. CpG islands (eukaryote-only)
    if is_eukaryote:
        if config.avoid_cpg:
            hard_constraints.append(NoCpGIslandConstraint(organism=organism))
    else:
        if config.avoid_cpg:
            logger.info(
                "Skipping eukaryote-only constraint NoCpGIslandConstraint "
                "for prokaryotic organism %s",
                organism,
            )

    # 2f. ATTTA instability motifs (relevant to both domains)
    if config.avoid_attta:
        hard_constraints.append(NoATTTAMotifConstraint())

    # 2g. T-runs (relevant to both domains)
    if config.avoid_t_runs:
        hard_constraints.append(NoTRunConstraint())

    # --- Step 3: Register soft constraints / objectives ---
    soft_constraints: list[SoftConstraint] = []

    # 3a. Maximize CAI
    soft_constraints.append(MaximizeCAI(adaptiveness=adaptiveness, protein=protein))

    # 3b. Minimize CpG (eukaryote-only objective)
    if is_eukaryote:
        if config.avoid_cpg:
            soft_constraints.append(MinimizeCpG(organism=organism))
    else:
        if config.avoid_cpg:
            logger.info(
                "Skipping eukaryote-only soft constraint MinimizeCpG "
                "for prokaryotic organism %s",
                organism,
            )

    # 3c. Minimize mRNA dG at 5' end
    soft_constraints.append(MinimizeMRNADG())

    # 3d. Minimize codon pair bias (optional -- off by default)
    if config.optimize_codon_pair_bias:
        cpb_constraint = MinimizeCodonPairBias(organism=organism)
        # Only add if CPB data is available for the organism
        if cpb_constraint.cpb_data:
            soft_constraints.append(cpb_constraint)
        else:
            logger.warning(
                "Codon pair bias optimization requested but no CPB data "
                "available for organism '%s'; skipping.", organism,
            )

    # --- Step 4: Build and return the model ---
    model = CSPModel(
        variables=variables,
        hard_constraints=hard_constraints,
        soft_constraints=soft_constraints,
        protein=protein,
        organism=organism,
        config=config,
    )

    logger.info(
        "Built CSP model: %d variables, %d hard constraints, %d soft constraints "
        "for %s (%s, %d aa)",
        model.num_variables,
        model.num_hard_constraints,
        model.num_soft_constraints,
        organism,
        protein[:20] + "..." if len(protein) > 20 else protein,
        len(protein),
    )

    return model


# ==============================================================================
# Internal validation helpers
# ==============================================================================

def _validate_protein(protein: str) -> None:
    """Validate that a protein string contains only standard amino acid codes.

    Args:
        protein: Uppercase protein string.

    Raises:
        ValueError: If any character is not a valid amino acid code.
    """
    valid_aas = set(AA_TO_CODONS.keys())
    for i, ch in enumerate(protein):
        if ch not in valid_aas:
            raise ValueError(
                f"Invalid amino acid '{ch}' at position {i} in protein. "
                f"Valid codes: {sorted(valid_aas)}"
            )


def _validate_organism(organism: str) -> None:
    """Validate that the organism is supported.

    Args:
        organism: Organism name string.

    Raises:
        ValueError: If the organism is not in SUPPORTED_ORGANISMS.
    """
    if organism not in SUPPORTED_ORGANISMS:
        raise ValueError(
            f"Unsupported organism '{organism}'. "
            f"Supported: {SUPPORTED_ORGANISMS}"
        )
