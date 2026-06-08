"""
Type stubs for biocompiler.optimizer.pipeline — BioOptimizer class and high-level API.
"""

from typing import Any, Dict, List, Optional, Tuple

from ..type_system import PredicateResult
from ..mutagenesis import MutagenesisReport, MutagenesisProposal
from ..decision_provenance import OptimizationDecisionTrail
from .utils import OptimizationResult, FullConstructResult

def optimize_sequence(
    protein: str,
    organism: str = ...,
    strategy: str = ...,
    enzymes: list[str] | None = ...,
    is_eukaryote: bool | None = ...,
    avoid_restriction_sites: bool = ...,
    gc_target: float | None = ...,
    max_iterations: int = ...,
    avoid_gt: bool = ...,
    avoid_cpg: bool = ...,
    track_provenance: bool = ...,
    objective: str | None = ...,
    include_utr: bool = ...,
    optimize_mrna_stability: bool = ...,
    consider_codon_pair_bias: bool = ...,
    strict: bool = ...,
    species: str | None = ...,
    **kwargs: Any,
) -> OptimizationResult: ...

def batch_optimize(
    proteins: list[str],
    organism: str = ...,
    strategy: str = ...,
    enzymes: list[str] | None = ...,
    is_eukaryote: bool | None = ...,
    avoid_restriction_sites: bool = ...,
    gc_target: float | None = ...,
    max_iterations: int = ...,
    avoid_gt: bool = ...,
    avoid_cpg: bool = ...,
    track_provenance: bool = ...,
    objective: str | None = ...,
    include_utr: bool = ...,
    optimize_mrna_stability: bool = ...,
    consider_codon_pair_bias: bool = ...,
    strict: bool = ...,
    species: str | None = ...,
    **kwargs: Any,
) -> list[OptimizationResult]: ...

class BioOptimizer:
    def __init__(
        self,
        organism: str = ...,
        enzymes: list[str] | None = ...,
        strategy: str = ...,
        is_eukaryote: bool | None = ...,
        species: str | None = ...,
    ) -> None: ...
    def optimize(
        self,
        seq: str,
        strategy: Optional[str] = ...,
    ) -> Tuple[str, List[PredicateResult], str]: ...

__all__: list[str]
