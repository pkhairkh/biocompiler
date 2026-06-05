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
from typing import Any, ClassVar, Optional, Protocol


__all__ = [
    # ── Enums ───────────────────────────────────────────────
    "SolverBackend",
    "ConstraintStrictness",
    "ConstraintType",
    "ConstraintPriority",
    "SolverStatus",

    # ── Data classes ────────────────────────────────────────
    "CodonVariable",
    "ConstraintSpec",
    "ConstraintViolation",
    "MUSReport",
    "SolverConfig",
    "SolverResult",
    "SpliceConstraint",
    "CrossCodonSpliceConstraint",
    "CSPModel",

    # ── Protocols ───────────────────────────────────────────
    "SolverBackendProtocol",
]


class SolverBackend(Enum):
    """Available solver backends for constraint satisfaction.

    NONE:
        No backend was used. Used for fallback results when all backends
        are unavailable or infeasible.

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

    NONE = "none"
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

    Note: Re-exported by solver/__init__.py as part of the solver package's
    public API.
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
    CODON_PAIR_BIAS = "codon_pair_bias"
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


class ConstraintPriority(Enum):
    """Enforcement priority for constraint scoring and violation ranking.

    Each priority level determines how severely a constraint violation
    impacts the composite solution score:

        CRITICAL: Must satisfy — violation yields a score of 0.0
            (infeasible solution).  Used for constraints like translation
            fidelity that cannot be compromised.
        HIGH: Severe penalty — violation heavily penalises the score.
            Used for constraints like restriction-site avoidance where
            violation has serious experimental consequences.
        MEDIUM: Moderate penalty — default priority.  Used for most
            biological constraints (GC range, splice sites, CpG islands).
        LOW: Minor penalty — violation has small score impact.
            Used for advisory constraints that are nice-to-have.

    The ``penalty_weight`` property returns a numeric multiplier used
    by :class:`ConstraintScorer` to compute weighted violation penalties.
    """

    CRITICAL = "critical"   # Must satisfy — violation = score 0.0
    HIGH = "high"           # Severe penalty
    MEDIUM = "medium"       # Moderate penalty (default)
    LOW = "low"             # Minor penalty

    @property
    def penalty_weight(self) -> float:
        """Numeric penalty multiplier for this priority level.

        Returns:
            Float multiplier used to weight violation penalties:
            - CRITICAL → 1000.0  (effectively infinite)
            - HIGH     → 10.0
            - MEDIUM   → 3.0
            - LOW      → 1.0
        """
        return {
            ConstraintPriority.CRITICAL: 1000.0,
            ConstraintPriority.HIGH: 10.0,
            ConstraintPriority.MEDIUM: 3.0,
            ConstraintPriority.LOW: 1.0,
        }[self]

    @property
    def rank(self) -> int:
        """Numeric rank for sorting (lower = more critical).

        Returns:
            Integer rank: CRITICAL=0, HIGH=1, MEDIUM=2, LOW=3.
        """
        return {
            ConstraintPriority.CRITICAL: 0,
            ConstraintPriority.HIGH: 1,
            ConstraintPriority.MEDIUM: 2,
            ConstraintPriority.LOW: 3,
        }[self]

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, ConstraintPriority):
            return NotImplemented
        return self.rank < other.rank

    def __le__(self, other: object) -> bool:
        if not isinstance(other, ConstraintPriority):
            return NotImplemented
        return self.rank <= other.rank

    def __gt__(self, other: object) -> bool:
        if not isinstance(other, ConstraintPriority):
            return NotImplemented
        return self.rank > other.rank

    def __ge__(self, other: object) -> bool:
        if not isinstance(other, ConstraintPriority):
            return NotImplemented
        return self.rank >= other.rank


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


def _is_eukaryotic(organism: str) -> bool:
    """Check whether *organism* is eukaryotic.

    Wraps :func:`biocompiler.organism_config.is_eukaryotic_organism`
    with a safe fallback (returns True for unknown organisms, matching
    the convention in :mod:`~biocompiler.organism_config`).

    This helper is used by :meth:`ConstraintSpec.check` to skip
    eukaryote-specific constraints (cryptic splice, CpG islands) when
    the target organism is prokaryotic.
    """
    try:
        from ..organism_config import is_eukaryotic_organism
        return is_eukaryotic_organism(organism)
    except Exception:
        # If organism config is unavailable, default to True (eukaryote)
        # to avoid false negatives.
        return True


@dataclass(frozen=True)
class ConstraintSpec:
    """A single constraint in the CSP model.

    Attributes:
        ctype: The constraint category classification.
        name: Human-readable name (e.g. "NoCrypticSplice_pos42").
        strictness: Whether this is a hard or soft constraint.
        params: Constraint-specific parameters (e.g. GC bounds, threshold).
        positions: Codon positions this constraint applies to (empty = global).
        priority: Enforcement priority for scoring — determines penalty weight
            when this constraint is violated.  CRITICAL constraints yield a
            composite score of 0.0 if violated; HIGH/MEDIUM/LOW impose
            progressively smaller penalties.  Defaults to MEDIUM.
        weight: Soft constraint weight multiplier for scoring.  Controls the
            relative contribution of this constraint to the overall solution
            score.  Defaults to 1.0.
    """

    ctype: ConstraintType
    name: str
    strictness: ConstraintStrictness = ConstraintStrictness.HARD
    params: dict[str, Any] = field(default_factory=dict)
    positions: list[int] = field(default_factory=list)
    priority: ConstraintPriority = ConstraintPriority.MEDIUM
    weight: float = 1.0

    def __hash__(self) -> int:
        # Intentional: only ctype + name used for hashing/equality.
        # Two constraints of the same type and name are considered identical
        # regardless of positions, params, or priority — this supports dedup
        # in sets and dicts (e.g. MUS constraint sets).
        return hash((self.ctype, self.name))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ConstraintSpec):
            return NotImplemented
        return self.ctype == other.ctype and self.name == other.name

    def check(self, sequence: str) -> bool:
        """Return True if the constraint is satisfied for *sequence*.

        Dispatches to the appropriate checking logic based on ``self.ctype``,
        using ``self.params`` for constraint-specific parameters (GC bounds,
        thresholds, forbidden sites, etc.).  Returns True for unknown
        constraint types (no violation).

        The checking logic is kept strictly consistent with the
        ``HardConstraint.check()`` implementations in
        :mod:`~biocompiler.solver.constraints` so that
        :meth:`ConstraintSpec.check` and the corresponding
        ``HardConstraint.check`` always agree on the same input.

        Parameters
        ----------
        sequence : str
            Full DNA sequence to validate.

        Returns
        -------
        bool
            True if the constraint is satisfied, False if violated.
        """
        if not sequence:
            return True

        seq = sequence.upper()

        # ── GC content ────────────────────────────────────────
        if self.ctype == ConstraintType.GC_CONTENT:
            gc_lo = self.params.get("gc_lo", 0.30)
            gc_hi = self.params.get("gc_hi", 0.70)
            gc_frac = sum(1 for b in seq if b in "GC") / len(seq)
            return gc_lo <= gc_frac <= gc_hi

        # ── CpG island avoidance ──────────────────────────────
        if self.ctype == ConstraintType.NO_CPG:
            # Prokaryotes lack DNA methylation → CpG islands are irrelevant.
            organism = self.params.get("organism", "")
            if organism and not _is_eukaryotic(organism):
                return True
            window = self.params.get("window", 200)
            threshold = self.params.get("threshold", 0.6)
            if len(seq) < window:
                return True
            for start in range(len(seq) - window + 1):
                w = seq[start : start + window]
                c_count = w.count("C")
                g_count = w.count("G")
                cg_count = sum(
                    1 for i in range(len(w) - 1) if w[i : i + 2] == "CG"
                )
                expected = (c_count * g_count) / window if window > 0 else 0
                if expected > 0 and cg_count / expected > threshold:
                    return False
            return True

        # ── Cryptic splice site avoidance ─────────────────────
        if self.ctype in (
            ConstraintType.NO_CRYPTIC_SPLICE,
            ConstraintType.SPLICE_DONOR_AVOIDANCE,
        ):
            # Prokaryotes lack spliceosomes → cryptic splice sites are
            # irrelevant.
            organism = self.params.get("organism", "")
            if organism and not _is_eukaryotic(organism):
                return True
            threshold = self.params.get("threshold", 3.0)
            protein = self.params.get("protein", "")

            # Compute unavoidable GT positions (Valine codons)
            unavoidable_gt: set[int] = set()
            if protein:
                for idx, aa in enumerate(protein):
                    if aa == "V":  # Valine codons all start with GT
                        unavoidable_gt.add(idx * 3)

            from ..maxentscan import score_donor, score_acceptor

            for i in range(len(seq) - 1):
                if seq[i : i + 2] == "GT":
                    # Skip GT within unavoidable Valine codons
                    if i in unavoidable_gt:
                        continue
                    if score_donor(seq, i) >= threshold:
                        return False
                if seq[i : i + 2] == "AG":
                    if score_acceptor(seq, i) >= threshold:
                        return False
            return True

        # ── Restriction site avoidance ────────────────────────
        if self.ctype == ConstraintType.RESTRICTION_SITE:
            sites = self.params.get("sites", [])
            from ..constants import reverse_complement

            for site in sites:
                site_upper = site.upper()
                if site_upper in seq:
                    return False
                rc = reverse_complement(site).upper()
                if rc != site_upper and rc in seq:
                    return False
            return True

        # ── Instability motif avoidance (ATTTA) ───────────────
        if self.ctype == ConstraintType.NO_INSTABILITY_MOTIF:
            motif = self.params.get("motif", "ATTTA")
            return motif.upper() not in seq

        # ── mRNA stability / T-run constraint ─────────────────
        if self.ctype == ConstraintType.MRNA_STABILITY:
            max_run = self.params.get("max_run", 5)
            run_length = 0
            for base in seq:
                if base == "T":
                    run_length += 1
                    if run_length > max_run:
                        return False
                else:
                    run_length = 0
            return True

        # ── Amino acid identity (translation fidelity) ────────
        if self.ctype == ConstraintType.AMINO_ACID_IDENTITY:
            protein = self.params.get("protein", "")
            if not protein:
                return True
            if len(sequence) != len(protein) * 3:
                return False
            from ..constants import CODON_TABLE

            for i, expected_aa in enumerate(protein):
                codon = sequence[i * 3 : i * 3 + 3]
                if CODON_TABLE.get(codon) != expected_aa:
                    return False
            return True

        # ── Codon usage (soft/optimization — always satisfied) ─
        if self.ctype == ConstraintType.CODON_USAGE:
            return True

        # ── Codon pair bias (soft/optimization — always satisfied) ─
        if self.ctype == ConstraintType.CODON_PAIR_BIAS:
            return True

        # ── GT dinucleotide avoidance ─────────────────────────
        if self.ctype == ConstraintType.NO_GT_DINUCLEOTIDE:
            return "GT" not in seq

        # ── Immunogenicity / MHC / T-cell (no sequence check) ─
        if self.ctype in (
            ConstraintType.MHC_BINDING,
            ConstraintType.TCELL_EPITOPE,
        ):
            # These require specialized immunogenicity models; default
            # to satisfied since we cannot check from sequence alone.
            return True

        # ── Protein stability (no sequence check) ─────────────
        if self.ctype == ConstraintType.PROTEIN_STABILITY:
            return True

        # ── Custom / unknown — assume satisfied ───────────────
        return True


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
        priority: Enforcement priority of the violated constraint.
            Used by :class:`ConstraintScorer` to rank violations by
            severity.  Defaults to MEDIUM.
        weight: Weight multiplier from the violated constraint's spec.
            Defaults to 1.0.
    """

    constraint_name: str
    constraint_type: ConstraintStrictness
    description: str
    positions: list[int] = field(default_factory=list)
    severity: float = 0.0  # 0-1, how badly violated
    priority: ConstraintPriority = ConstraintPriority.MEDIUM
    weight: float = 1.0

    def __hash__(self) -> int:
        # Consistent with ConstraintSpec: identity is name + type.
        # Two violations of the same constraint name and type are considered
        # identical regardless of positions or severity — supports dedup in
        # sets and dicts.
        return hash((self.constraint_name, self.constraint_type))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ConstraintViolation):
            return NotImplemented
        return (
            self.constraint_name == other.constraint_name
            and self.constraint_type == other.constraint_type
        )


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


# ---------------------------------------------------------------------------
# Organism-specific presets for SolverConfig.for_organism()
# ---------------------------------------------------------------------------

_ORGANISM_PRESETS: dict[str, dict[str, Any]] = {
    "E_coli_K12": {
        "cai_weight": 1.0, "cpg_weight": 0.0, "mrna_dg_weight": 0.2,
        "gc_lo": 0.45, "gc_hi": 0.55, "avoid_cpg": False,
        "cryptic_splice_threshold": 0.0,
    },
    "E_coli_BL21": {
        "cai_weight": 1.0, "cpg_weight": 0.0, "mrna_dg_weight": 0.2,
        "gc_lo": 0.45, "gc_hi": 0.55, "avoid_cpg": False,
        "cryptic_splice_threshold": 0.0,
    },
    "Homo_sapiens": {
        "cai_weight": 1.0, "cpg_weight": 0.5, "mrna_dg_weight": 0.3,
        "gc_lo": 0.40, "gc_hi": 0.60, "avoid_cpg": True,
        "cryptic_splice_threshold": 3.0,
    },
    "Mus_musculus": {
        "cai_weight": 1.0, "cpg_weight": 0.4, "mrna_dg_weight": 0.2,
        "gc_lo": 0.40, "gc_hi": 0.55, "avoid_cpg": True,
        "cryptic_splice_threshold": 3.0,
    },
    "CHO_K1": {
        "cai_weight": 1.0, "cpg_weight": 0.5, "mrna_dg_weight": 0.3,
        "gc_lo": 0.40, "gc_hi": 0.60, "avoid_cpg": True,
        "cryptic_splice_threshold": 3.0,
    },
    "Saccharomyces_cerevisiae": {
        "cai_weight": 1.0, "cpg_weight": 0.1, "mrna_dg_weight": 0.2,
        "gc_lo": 0.35, "gc_hi": 0.45, "avoid_cpg": False,
        "cryptic_splice_threshold": 3.0,
    },
}

_ORGANISM_ALIASES: dict[str, str] = {
    "E. coli": "E_coli_K12", "ecoli": "E_coli_K12",
    "E_coli": "E_coli_K12", "Escherichia_coli": "E_coli_K12",
    "human": "Homo_sapiens", "mouse": "Mus_musculus",
    "cho": "CHO_K1", "yeast": "Saccharomyces_cerevisiae",
}


@dataclass
class SolverConfig:
    """Configuration for the CSP solver.

    Controls which backend to use, timeout, constraint parameters, and
    objective weights. All parameters have sensible defaults for mammalian
    gene optimization.

    Abbreviation convention: short field names (gc_lo, gc_hi) are used for
    frequently-referenced numeric bounds to keep call-sites concise, while
    descriptive names (timeout_seconds, solve_time_seconds) are used for
    one-off or less frequently used fields. Readable aliases (gc_low,
    gc_high) are provided as @property accessors.

    Attributes:
        organism: Target organism name (e.g. ``"Homo_sapiens"``,
            ``"E_coli_K12"``).  Stored on the config so that every stage
            of the solver pipeline can make organism-aware decisions
            without requiring a separate parameter.  Defaults to
            ``"Homo_sapiens"``.
        auto_detect_organism_domain: When ``True`` (default), the solver
            automatically determines whether the target organism is
            eukaryotic and skips eukaryote-only constraints (cryptic
            splice sites, CpG islands) for prokaryotic targets.  Set to
            ``False`` to treat the organism as eukaryotic regardless.
        backend: Which solver backend to use.
        timeout_seconds: Maximum solve time before falling back.
        max_codons: Maximum sequence length the solver handles.
            Sequences longer than this are split or fall back to greedy.
        gc_lo: Minimum acceptable GC content fraction (hard constraint).
            Alias: gc_low
        gc_hi: Maximum acceptable GC content fraction (hard constraint).
            Alias: gc_high
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

    organism: str = "Homo_sapiens"  # Target organism for the solver pipeline
    auto_detect_organism_domain: bool = True  # Auto-skip eukaryote-only constraints for prokaryotes
    backend: SolverBackend = SolverBackend.ORTOOLS
    timeout_seconds: float = 60.0
    max_codons: int = 5000  # Max sequence length the solver handles
    gc_lo: float = 0.30  # Abbreviated; use gc_low property for readable access
    gc_hi: float = 0.70  # Abbreviated; use gc_high property for readable access
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
    optimize_codon_pair_bias: bool = False  # Enable codon pair bias soft constraint
    codon_pair_bias_weight: float = 0.2  # Weight for codon pair bias objective
    cai_reference_set: str = "kazusa"  # CAI reference set: "kazusa" (default) or "sharp_li"
    verbose: bool = False


    # Valid values for cai_reference_set
    _VALID_CAI_REFERENCE_SETS: ClassVar[tuple[str, ...]] = ("kazusa", "sharp_li")

    def __post_init__(self) -> None:
        """Validate cai_reference_set is one of the allowed values."""
        if self.cai_reference_set not in self._VALID_CAI_REFERENCE_SETS:
            raise ValueError(
                f"Invalid cai_reference_set '{self.cai_reference_set}'."
                f" Valid values: {list(self._VALID_CAI_REFERENCE_SETS)}"
            )

    @classmethod
    def for_organism(cls, organism: str) -> SolverConfig:
        """Return a pre-tuned SolverConfig for the given organism.

        Uses organism-specific scoring weights and constraint parameters.
        For prokaryotes, CpG-avoidance is irrelevant (no DNA methylation),
        so ``cpg_weight`` is set to 0.0.

        The *organism* argument accepts both canonical keys (e.g.
        ``"E_coli_K12"``) and legacy aliases (e.g. ``"ecoli"``,
        ``"human"``).  If the organism is not found, a default config
        is returned with a warning.
        """
        import logging as _logging
        _logger = _logging.getLogger(__name__)

        canonical = _ORGANISM_ALIASES.get(organism, organism)
        preset = _ORGANISM_PRESETS.get(canonical)
        if preset is None:
            _logger.warning(
                "No preset config for organism %r; using default SolverConfig",
                organism,
            )
            return cls()

        config = cls(**preset)
        _logger.info(
            "Created SolverConfig for organism %r: "
            "cai=%.2f, cpg=%.2f, mrna_dg=%.2f, gc=[%.2f,%.2f], avoid_cpg=%s",
            organism, config.cai_weight, config.cpg_weight,
            config.mrna_dg_weight, config.gc_lo, config.gc_hi,
            config.avoid_cpg,
        )
        return config

    @property
    def gc_low(self) -> float:
        """Readable alias for gc_lo."""
        return self.gc_lo

    @property
    def gc_high(self) -> float:
        """Readable alias for gc_hi."""
        return self.gc_hi

    @property
    def effective_donor_threshold(self) -> float:
        """Donor threshold, falling back to cryptic_splice_threshold if not set."""
        return self.donor_threshold if self.donor_threshold is not None else self.cryptic_splice_threshold

    @property
    def effective_acceptor_threshold(self) -> float:
        """Acceptor threshold, falling back to cryptic_splice_threshold if not set."""
        return self.acceptor_threshold if self.acceptor_threshold is not None else self.cryptic_splice_threshold

    @property
    def is_eukaryotic(self) -> bool:
        """Whether the target organism is eukaryotic.

        When ``auto_detect_organism_domain`` is ``True`` (the default),
        this property delegates to :func:`is_eukaryotic_organism` using
        the configured ``organism`` name to determine the domain of life.
        When ``auto_detect_organism_domain`` is ``False``, the property
        assumes the organism is eukaryotic (conservative default).

        Returns:
            ``True`` if the organism is eukaryotic, ``False`` otherwise.
        """
        if self.auto_detect_organism_domain:
            from ..organism_config import is_eukaryotic_organism
            return is_eukaryotic_organism(self.organism)
        # When auto-detect is off, assume eukaryote (conservative default)
        return True


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
        protein: The input amino acid sequence.
        organism: The target organism name.
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
        metadata: Additional key-value metadata (e.g. reason for fallback).
    """

    sequence: str
    solved: bool
    backend_used: SolverBackend
    protein: str = ""
    organism: str = ""
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
    metadata: dict[str, Any] = field(default_factory=dict)


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
        organism: Target organism name (e.g. ``"Homo_sapiens"``).
    """

    protein_sequence: str
    codon_domains: dict[int, list[str]]
    constraints: list[ConstraintSpec]
    config: SolverConfig = field(default_factory=SolverConfig)
    organism: str = ""

    @property
    def length(self) -> int:
        """Number of codon positions."""
        return len(self.protein_sequence)

    @property
    def protein(self) -> str:
        """Alias for :attr:`protein_sequence` — compatibility with other CSPModel variants."""
        return self.protein_sequence

    @property
    def hard_constraints(self) -> list[ConstraintSpec]:
        """Hard constraints derived from :attr:`constraints`.

        Provides compatibility with :class:`constraints.CSPModel` which
        stores hard and soft constraints separately.
        """
        return [c for c in self.constraints if c.strictness == ConstraintStrictness.HARD]

    @property
    def soft_constraints(self) -> list[ConstraintSpec]:
        """Soft constraints derived from :attr:`constraints`.

        Provides compatibility with :class:`constraints.CSPModel` which
        stores hard and soft constraints separately.
        """
        return [c for c in self.constraints if c.strictness == ConstraintStrictness.SOFT]


# ---------------------------------------------------------------------------
# Solver backend protocol (for MUS dependency injection)
# ---------------------------------------------------------------------------

class SolverBackendProtocol(Protocol):
    """Protocol that any CSP/SMT solver backend must implement.

    This allows the MUS module to call the solver without depending
    on a specific implementation (Z3, OR-Tools, etc.).

    Note: Re-exported by solver/__init__.py as part of the solver package's
    public API.
    """
    # Protocol for MUS dependency injection — used internally by solver/mus.py

    def solve(self, model: CSPModel) -> SolverResult:
        """Solve the CSP model and return the result."""
        ...
