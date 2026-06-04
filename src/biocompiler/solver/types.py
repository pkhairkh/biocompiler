"""CSP Solver shared types and data structures.

Defines the core data model for constraint-based gene optimization:
- Variables: codon positions with finite domains (AA_TO_CODONS)
- Constraints: hard (must satisfy) and soft (optimize)
- Results: solver output with metadata

This module is the single source of truth for all solver-related types.
Other solver sub-modules (engines, constraints, etc.) import from here,
never the reverse.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional, Protocol


class SolverBackend(Enum):
    """Available solver backends for constraint satisfaction.

    ORTOOLS:
        Google OR-Tools CP-SAT solver. Primary backend.
        Supports automaton constraints for forbidden substrings (restriction
        sites, splice motifs) and table constraints for cross-codon coordination.
        Best for discrete, combinatorial problems.

    Z3:
        Microsoft Z3 SMT solver. Fallback backend.
        Rich theory support for continuous constraints (GC content ranges,
        weighted CAI optimization). Better for mixed integer-real problems.

    GREEDY_FALLBACK:
        Deterministic greedy algorithm. Used when neither OR-Tools nor Z3
        is installed. Sequential constraint resolution — no global optimality
        guarantee, but always available.
    """

    ORTOOLS = "ortools"
    Z3 = "z3"
    GREEDY_FALLBACK = "greedy"


class ConstraintStrictness(Enum):
    """Classification of constraint strictness.

    HARD:
        Must satisfy. If a hard constraint cannot be satisfied, the problem
        is infeasible and the solver reports a MUS (Minimal Unsatisfiable
        Subset) identifying the conflicting constraints.

    SOFT:
        Prefer to satisfy. Soft constraints contribute to the optimization
        objective but may be violated if they conflict with hard constraints.
        Violation severity is tracked in the result.
    """

    HARD = "hard"      # Must satisfy (infeasible if violated)
    SOFT = "soft"      # Prefer to satisfy (optimization objective)


class ConstraintType(str, Enum):
    """Enumeration of all supported CSP constraint categories.

    Each value corresponds to a biological design rule that the solver
    enforces. Used by the MUS module to classify and explain conflicts.
    """

    # Sequence composition constraints
    GC_CONTENT = "gc_content"
    CODON_USAGE = "codon_usage"
    NO_CPG = "no_cpg"
    NO_GT_DINUCLEOTIDE = "no_gt_dinucleotide"

    # Splicing constraints
    NO_CRYPTIC_SPLICE = "no_cryptic_splice"
    SPLICE_DONOR_AVOIDANCE = "splice_donor_avoidance"

    # Restriction site constraints
    RESTRICTION_SITE = "restriction_site"

    # mRNA structure constraints
    MRNA_STABILITY = "mrna_stability"
    NO_INSTABILITY_MOTIF = "no_instability_motif"

    # Immunogenicity constraints
    MHC_BINDING = "mhc_binding"
    TCELL_EPITOPE = "tcell_epitope"

    # Protein-level constraints
    AMINO_ACID_IDENTITY = "amino_acid_identity"
    PROTEIN_STABILITY = "protein_stability"

    # Custom / user-defined
    CUSTOM = "custom"


class SolverStatus(str, Enum):
    """Solver run outcome."""

    OPTIMAL = "optimal"
    SATISFIED = "satisfied"
    INFEASIBLE = "infeasible"
    UNKNOWN = "unknown"
    TIMEOUT = "timeout"
    ERROR = "error"


@dataclass
class CodonVariable:
    """A decision variable representing a codon position in the gene sequence.

    Each codon position is a finite-domain variable whose domain is the set
    of synonymous codons for its amino acid. The solver assigns exactly one
    codon from the domain to each variable.

    Attributes:
        position: 0-based index of this codon in the protein sequence.
        amino_acid: Single-letter amino acid code (e.g. 'M', 'V', 'K').
        domain: Ordered list of possible codons for this amino acid.
            The order reflects CAI preference (highest first).
        current_value: The assigned codon, or None if unassigned.
    """

    position: int
    amino_acid: str
    domain: list[str]  # Possible codons for this AA
    current_value: Optional[str] = None


@dataclass(frozen=True)
class ConstraintSpec:
    """A single constraint in the CSP model.

    Attributes:
        ctype: The constraint category classification.
        name: Human-readable name (e.g. "NoCrypticSplice_pos42").
        strictness: Whether this is a hard or soft constraint.
        params: Constraint-specific parameters (e.g. GC bounds, threshold).
        positions: Codon positions this constraint applies to (empty = global).
        priority: Suggestion priority for relaxation (lower = harder to relax).
    """

    ctype: ConstraintType
    name: str
    strictness: ConstraintStrictness = ConstraintStrictness.HARD
    params: dict[str, Any] = field(default_factory=dict)
    positions: list[int] = field(default_factory=list)
    priority: int = 5  # 1=critical, 10=easiest to relax

    def __hash__(self) -> int:
        return hash((self.ctype, self.name))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ConstraintSpec):
            return NotImplemented
        return self.ctype == other.ctype and self.name == other.name


@dataclass
class ConstraintViolation:
    """A constraint that is violated in a solution.

    Captures both the identity and severity of a violated constraint,
    enabling actionable diagnostics for the user.

    Attributes:
        constraint_name: Human-readable identifier (e.g. 'no_EcoRI_site').
        constraint_type: Whether this was a hard or soft constraint.
        description: Human-readable explanation of the violation.
        positions: Nucleotide positions involved in the violation.
        severity: How badly the constraint is violated, on a 0-1 scale.
            0.0 = marginally violated, 1.0 = severely violated.
    """

    constraint_name: str
    constraint_type: ConstraintStrictness
    description: str
    positions: list[int] = field(default_factory=list)
    severity: float = 0.0  # 0-1, how badly violated


@dataclass
class MUSReport:
    """Minimal Unsatisfiable Subset — the smallest set of conflicting constraints.

    When the solver determines that hard constraints cannot all be satisfied
    simultaneously, it computes a MUS: a minimal subset of constraints that
    is itself unsatisfiable. Removing any single constraint from the MUS
    makes the remaining constraints satisfiable.

    This enables targeted relaxation: the user can decide which constraint
    to relax based on the suggested relaxations.

    Attributes:
        mus_constraints: The minimal set of conflicting constraint specs.
        conflicting_constraints: Names of constraints in the MUS (legacy).
        all_constraints: The full constraint set (for context).
        explanation: Human-readable explanation of the conflict.
        suggested_relaxations: Ordered list of constraints to relax,
            from least impactful to most impactful.
        iterations: Number of solver calls during MUS extraction.
        solve_time_seconds: Total time for MUS computation.
        conflict_positions: Codon positions where conflicts manifest.
    """

    mus_constraints: list[ConstraintSpec] = field(default_factory=list)
    conflicting_constraints: list[str] = field(default_factory=list)
    all_constraints: list[ConstraintSpec] = field(default_factory=list)
    explanation: str = ""
    suggested_relaxations: list[str] = field(default_factory=list)
    iterations: int = 0
    solve_time_seconds: float = 0.0
    conflict_positions: list[int] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Populate conflicting_constraints from mus_constraints if empty."""
        if not self.conflicting_constraints and self.mus_constraints:
            self.conflicting_constraints = [c.name for c in self.mus_constraints]


@dataclass
class SolverConfig:
    """Configuration for the CSP solver.

    Controls which backend to use, timeout, constraint parameters, and
    objective weights. All parameters have sensible defaults for mammalian
    gene optimization.

    Attributes:
        backend: Which solver backend to use.
        timeout_seconds: Maximum solve time before falling back.
        max_codons: Maximum sequence length the solver handles.
            Sequences longer than this are split or fall back to greedy.
        gc_lo: Minimum acceptable GC content fraction (hard constraint).
        gc_hi: Maximum acceptable GC content fraction (hard constraint).
        cryptic_splice_threshold: MaxEnt score threshold for cryptic splice
            site detection. Sites scoring above this are forbidden.
            Used as default for both donor and acceptor thresholds unless
            donor_threshold or acceptor_threshold are explicitly set.
        donor_threshold: MaxEntScan score threshold specific to donor sites.
            If None, falls back to cryptic_splice_threshold.
        acceptor_threshold: MaxEntScan score threshold specific to acceptor sites.
            If None, falls back to cryptic_splice_threshold.
        n_quantize_bins: Number of bins for quantizing continuous MaxEntScan
            scores into discrete CSP domain values (default 20).
        restriction_sites: List of restriction enzyme recognition sequences
            to avoid (e.g. ['GAATTC', 'GGATCC']).
        avoid_cpg: Whether to minimize CpG dinucleotides.
        avoid_attta: Whether to avoid ATTTA instability motifs.
        avoid_t_runs: Whether to break runs of 6+ consecutive T bases.
        cai_weight: Weight for CAI maximization in the objective.
        cpg_weight: Weight for CpG minimization in the objective.
        mrna_dg_weight: Weight for mRNA structure stability (ViennaRNA).
        verbose: Whether to emit detailed solver diagnostics.
    """

    backend: SolverBackend = SolverBackend.ORTOOLS
    timeout_seconds: float = 60.0
    max_codons: int = 5000  # Max sequence length the solver handles
    gc_lo: float = 0.30
    gc_hi: float = 0.70
    cryptic_splice_threshold: float = 3.0
    donor_threshold: Optional[float] = None  # Falls back to cryptic_splice_threshold
    acceptor_threshold: Optional[float] = None  # Falls back to cryptic_splice_threshold
    n_quantize_bins: int = 20
    restriction_sites: list[str] = field(default_factory=list)
    avoid_cpg: bool = True
    avoid_attta: bool = True
    avoid_t_runs: bool = True
    cai_weight: float = 1.0
    cpg_weight: float = 0.5
    mrna_dg_weight: float = 0.3  # Weight for mRNA structure objective (ViennaRNA)
    verbose: bool = False

    @property
    def effective_donor_threshold(self) -> float:
        """Donor threshold, falling back to cryptic_splice_threshold if not set."""
        return self.donor_threshold if self.donor_threshold is not None else self.cryptic_splice_threshold

    @property
    def effective_acceptor_threshold(self) -> float:
        """Acceptor threshold, falling back to cryptic_splice_threshold if not set."""
        return self.acceptor_threshold if self.acceptor_threshold is not None else self.cryptic_splice_threshold


@dataclass
class SolverResult:
    """Result from the CSP solver.

    Contains the optimized DNA sequence, solution status, quality metrics,
    and diagnostic information. Always check ``solved`` before using the
    sequence — if False, the sequence may be a best-effort fallback.

    Attributes:
        sequence: The optimized DNA sequence (or best-effort if unsolved).
        solved: Whether a valid solution satisfying all hard constraints
            was found.
        backend_used: Which solver backend produced this result.
        cai: Codon Adaptation Index of the result sequence.
        gc_content: GC fraction of the result sequence.
        solve_time_seconds: Wall-clock time spent solving.
        violations: List of constraints violated in this solution.
        mus_report: If infeasible, the Minimal Unsatisfiable Subset.
        objective_value: Value of the optimization objective (higher = better).
        num_constraints: Total number of constraints posted to the solver.
        num_variables: Total number of decision variables (codon positions).
        fallback_used: Whether the greedy fallback was used because the
            primary solver was unavailable or timed out.
        warnings: Non-fatal diagnostic messages.
    """

    sequence: str
    solved: bool
    backend_used: SolverBackend
    cai: float = 0.0
    gc_content: float = 0.0
    solve_time_seconds: float = 0.0
    violations: list[ConstraintViolation] = field(default_factory=list)
    mus_report: Optional[MUSReport] = None
    objective_value: float = 0.0
    num_constraints: int = 0
    num_variables: int = 0
    fallback_used: bool = False
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Splice-specific constraint types (added by A4: maxent_encoding)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SpliceConstraint:
    """A constraint forbidding a specific codon assignment at a position.

    Encodes the fact that choosing ``codon`` at ``position`` would create a
    cryptic splice site (donor or acceptor) with MaxEntScan score exceeding
    the configured threshold.

    Attributes:
        position: 0-indexed codon position in the protein.
        codon: The 3-letter DNA codon that is forbidden.
        site_type: "donor" or "acceptor".
        score: The MaxEntScan score that would result from this codon choice.
        threshold: The configured threshold that this score exceeds.
    """

    position: int
    codon: str
    site_type: str  # "donor" or "acceptor"
    score: float
    threshold: float

    @property
    def margin(self) -> float:
        """How far above the threshold this score is (positive = violation)."""
        return self.score - self.threshold


@dataclass(frozen=True)
class CrossCodonSpliceConstraint:
    """A constraint involving two adjacent codon positions.

    Some GT/AG dinucleotides span codon boundaries. This constraint forbids
    specific (codon_left, codon_right) pairs that would create such
    cross-codon cryptic splice sites.

    Attributes:
        position_left: 0-indexed position of the left codon.
        position_right: 0-indexed position of the right codon.
        codon_left: The codon at the left position.
        codon_right: The codon at the right position.
        is_donor: True if the splice site is a donor (GT), False for acceptor (AG).
        score: The MaxEntScan score for this cross-codon splice site.
    """

    position_left: int
    position_right: int
    codon_left: str
    codon_right: str
    is_donor: bool
    score: float


# ---------------------------------------------------------------------------
# CSP Model
# ---------------------------------------------------------------------------

@dataclass
class CSPModel:
    """The full constraint satisfaction problem model.

    Attributes:
        protein_sequence: The amino acid sequence to codon-optimize.
        codon_domains: Mapping from position index to list of allowed codons.
        constraints: Ordered list of constraints to satisfy.
        config: Solver configuration.
    """

    protein_sequence: str
    codon_domains: dict[int, list[str]]
    constraints: list[ConstraintSpec]
    config: SolverConfig = field(default_factory=SolverConfig)

    @property
    def length(self) -> int:
        """Number of codon positions."""
        return len(self.protein_sequence)


# ---------------------------------------------------------------------------
# Solver backend protocol (for MUS dependency injection)
# ---------------------------------------------------------------------------

class SolverBackendProtocol(Protocol):
    """Protocol that any CSP/SMT solver backend must implement.

    This allows the MUS module to call the solver without depending
    on a specific implementation (Z3, OR-Tools, etc.).
    """

    def solve(self, model: CSPModel) -> SolverResult:
        """Solve the CSP model and return the result."""
        ...
