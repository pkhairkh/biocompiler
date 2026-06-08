"""Type stubs for biocompiler.solver — CSP solver, config, result, and constraints."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional, Protocol, Union


# ────────────────────────────────────────────────────────────
# Enums
# ────────────────────────────────────────────────────────────

class SolverBackend(Enum):
    NONE: str
    ORTOOLS: str
    Z3: str
    GREEDY_FALLBACK: str


class ConstraintStrictness(Enum):
    HARD: str
    SOFT: str


class ConstraintType(str, Enum):
    GC_CONTENT: str
    CODON_USAGE: str
    CODON_PAIR_BIAS: str
    NO_CPG: str
    NO_GT_DINUCLEOTIDE: str
    NO_CRYPTIC_SPLICE: str
    SPLICE_DONOR_AVOIDANCE: str
    RESTRICTION_SITE: str
    PATTERN_ENFORCEMENT: str
    PATTERN_AVOIDANCE: str
    MRNA_STABILITY: str
    NO_INSTABILITY_MOTIF: str
    MHC_BINDING: str
    TCELL_EPITOPE: str
    AMINO_ACID_IDENTITY: str
    PROTEIN_STABILITY: str
    CUSTOM: str


class ConstraintPriority(Enum):
    CRITICAL: str
    HIGH: str
    MEDIUM: str
    LOW: str

    @property
    def penalty_weight(self) -> float: ...

    @property
    def rank(self) -> int: ...


class SolverStatus(str, Enum):
    OPTIMAL: str
    SATISFIED: str
    INFEASIBLE: str
    UNKNOWN: str
    TIMEOUT: str
    ERROR: str


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class CodonVariable:
    position: int
    amino_acid: str
    domain: list[str]
    current_value: Optional[str]


@dataclass(frozen=True)
class ConstraintSpec:
    ctype: ConstraintType
    name: str
    strictness: ConstraintStrictness
    params: dict[str, Any]
    positions: list[int]
    priority: ConstraintPriority
    weight: float

    def check(self, sequence: str) -> bool: ...


@dataclass
class ConstraintViolation:
    constraint_name: str
    constraint_type: ConstraintStrictness
    description: str
    positions: list[int]
    severity: float
    priority: ConstraintPriority
    weight: float


@dataclass
class MUSReport:
    mus_constraints: list[ConstraintSpec]
    conflicting_constraints: list[str]
    all_constraints: list[ConstraintSpec]
    explanation: str
    suggested_relaxations: list[str]
    iterations: int
    solve_time_seconds: float
    conflict_positions: list[int]


@dataclass
class InfeasibilityReport:
    is_infeasible: bool
    conflicting_constraints: list[str]
    unsat_core: list[str] | None
    suggested_relaxations: list[str]
    explanation: str
    detection_method: str


@dataclass
class SolverConfig:
    organism: str
    auto_detect_organism_domain: bool
    backend: SolverBackend
    timeout_seconds: float
    max_codons: int
    gc_lo: float
    gc_hi: float
    cryptic_splice_threshold: float
    donor_threshold: Optional[float]
    acceptor_threshold: Optional[float]
    n_quantize_bins: int
    restriction_sites: list[str]
    avoid_cpg: bool
    avoid_attta: bool
    avoid_t_runs: bool
    add_default_restriction_sites: bool
    cai_weight: float
    cpg_weight: float
    mrna_dg_weight: float
    optimize_codon_pair_bias: bool
    codon_pair_bias_weight: float
    cai_reference_set: str
    verbose: bool

    @classmethod
    def for_organism(cls, organism: str) -> SolverConfig: ...

    @property
    def gc_low(self) -> float: ...

    @property
    def gc_high(self) -> float: ...

    @property
    def effective_donor_threshold(self) -> float: ...

    @property
    def effective_acceptor_threshold(self) -> float: ...

    @property
    def is_eukaryotic(self) -> bool: ...


@dataclass
class SolverResult:
    sequence: str
    solved: bool
    backend_used: SolverBackend
    protein: str
    organism: str
    cai: float
    gc_content: float
    solve_time_seconds: float
    violations: list[ConstraintViolation]
    mus_report: Optional[MUSReport]
    infeasibility_report: Optional[InfeasibilityReport]
    objective_value: float
    num_constraints: int
    num_variables: int
    fallback_used: bool
    warnings: list[str]
    metadata: dict[str, Any]

    @property
    def backend(self) -> SolverBackend: ...


@dataclass(frozen=True)
class SpliceConstraint:
    position: int
    codon: str
    site_type: str
    score: float
    threshold: float

    @property
    def margin(self) -> float: ...


@dataclass(frozen=True)
class CrossCodonSpliceConstraint:
    position_left: int
    position_right: int
    codon_left: str
    codon_right: str
    site_type: str
    score: float
    threshold: float


class CSPModel:
    protein: str
    organism: str
    variables: list[CodonVariable]
    hard_constraints: list[ConstraintSpec]
    soft_constraints: list[ConstraintSpec]
    config: SolverConfig


# ────────────────────────────────────────────────────────────
# Protocol
# ────────────────────────────────────────────────────────────

class SolverBackendProtocol(Protocol):
    def solve(self, protein: str, organism: str) -> SolverResult: ...


EngineType = Union["ORTOOLSEngine", "Z3Engine", "GreedyEngine"]


# ────────────────────────────────────────────────────────────
# CSPSolver class
# ────────────────────────────────────────────────────────────

class CSPSolver:
    config: SolverConfig

    def __init__(self, config: Optional[SolverConfig] = ...) -> None: ...
    def solve(self, protein: str, organism: str = ..., **overrides: Any) -> SolverResult: ...


# ────────────────────────────────────────────────────────────
# Convenience functions
# ────────────────────────────────────────────────────────────

def solve(protein: str, organism: str = ..., config: Optional[SolverConfig] = ..., **overrides: Any) -> SolverResult: ...
def is_solver_available(backend: Optional[SolverBackend] = ...) -> bool: ...


# ────────────────────────────────────────────────────────────
# Dispatch functions
# ────────────────────────────────────────────────────────────

def solve_with_csp(protein: str, organism: str = ..., **kwargs: Any) -> SolverResult: ...
def solve_csp(protein: str, organism: str = ..., **kwargs: Any) -> SolverResult: ...
def get_csp_availability() -> bool: ...
def csp_optimize(protein: str, organism: str = ..., **kwargs: Any) -> SolverResult: ...
def validate_csp_solution(sequence: str, protein: str, organism: str = ..., **kwargs: Any) -> bool: ...


# ────────────────────────────────────────────────────────────
# Constraint model
# ────────────────────────────────────────────────────────────

class HardConstraint:
    def check(self, sequence: str) -> bool: ...

class SoftConstraint:
    def check(self, sequence: str) -> bool: ...

class TranslationConstraint(HardConstraint): ...
class NoRestrictionSiteConstraint(HardConstraint): ...
class GCRangeConstraint(HardConstraint): ...
class NoCrypticSpliceConstraint(HardConstraint): ...
class NoCpGIslandConstraint(HardConstraint): ...
class NoATTTAMotifConstraint(HardConstraint): ...
class NoTRunConstraint(HardConstraint): ...
class MaximizeCAI(SoftConstraint): ...
class MinimizeCpG(SoftConstraint): ...
class MinimizeCodonPairBias(SoftConstraint): ...
class MinimizeMRNADG(SoftConstraint): ...

def build_csp_model(protein: str, organism: str = ..., config: Optional[SolverConfig] = ...) -> CSPModel: ...
def codon_gc_count(codon: str) -> int: ...
def codon_contains_gt(codon: str) -> bool: ...
def codon_contains_ag(codon: str) -> bool: ...
def codon_contains_cpg(codon: str) -> bool: ...
def compute_gc_from_codons(codons: list[str]) -> float: ...


# ────────────────────────────────────────────────────────────
# MaxEntScan encoding
# ────────────────────────────────────────────────────────────

class SpliceConstraintEncoder: ...
def precompute_splice_scores(**kwargs: Any) -> dict[str, Any]: ...
def build_splice_constraint_table(**kwargs: Any) -> dict[str, Any]: ...
def quantize_maxent_scores(**kwargs: Any) -> dict[str, Any]: ...
def encode_cross_codon_splice_context(**kwargs: Any) -> dict[str, Any]: ...


# ────────────────────────────────────────────────────────────
# Enforcement
# ────────────────────────────────────────────────────────────

class ConstraintEnforcer: ...
class EnforcementResult: ...


# ────────────────────────────────────────────────────────────
# Scoring
# ────────────────────────────────────────────────────────────

class ConstraintScorer: ...
class ScoringResult: ...
class SoftConstraintScorer: ...
def compute_pareto_frontier(**kwargs: Any) -> list[Any]: ...


# ────────────────────────────────────────────────────────────
# Conflict resolution
# ────────────────────────────────────────────────────────────

class ConstraintConflict: ...
class ConflictResolver: ...
def prioritize_constraints(**kwargs: Any) -> list[Any]: ...


# ────────────────────────────────────────────────────────────
# Conflict provenance
# ────────────────────────────────────────────────────────────

class ConflictProvenance: ...
class ConflictResolverWithProvenance: ...


# ────────────────────────────────────────────────────────────
# Constraint interaction map
# ────────────────────────────────────────────────────────────

class ConstraintInteractionMap: ...
class InteractionInfo: ...
def print_interaction_report(**kwargs: Any) -> None: ...


# ────────────────────────────────────────────────────────────
# CAI-aware resolver
# ────────────────────────────────────────────────────────────

class CAIAwareConstraintResolver: ...
