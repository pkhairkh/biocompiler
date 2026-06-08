"""
Type stubs for biocompiler.optimizer.greedy — Greedy optimizer.
"""

from typing import Any, Callable, Dict, List, Optional, Set, Tuple

from ..type_system import PredicateResult
from ..mutagenesis import MutagenesisReport, MutagenesisProposal
from ..decision_provenance import OptimizationDecisionTrail
from .cai import _BatchSwapScorer

def score_splice_donor_potential(dna: str, position: int) -> float: ...
def _gt_aware_select_codon(
    aa: str,
    species_cai: dict[str, float],
    *,
    avoid_gt: bool = ...,
    avoid_ag: bool = ...,
    is_eukaryote: bool = ...,
    top_k: int = ...,
    scorer: _BatchSwapScorer | None = ...,
    codon_idx: int = ...,
    seq_codons: list[str] | None = ...,
) -> str: ...
def _is_in_codon_gt(seq: str, pos: int) -> bool: ...
def _eukaryote_cai_recovery(
    seq: str,
    species_cai: dict[str, float],
    is_eukaryote: bool = ...,
    organism: str = ...,
) -> str: ...
def _eliminate_cpg_dinucleotides(
    seq: str,
    species_cai: dict[str, float],
    organism: str,
    is_eukaryote: bool = ...,
    enzyme_list: list[str] | None = ...,
) -> str: ...
def _greedy_optimize(
    protein: str,
    organism: str,
    strategy: str = ...,
    enzymes: list[str] | None = ...,
    is_eukaryote: bool = ...,
    avoid_restriction_sites: bool = ...,
    gc_target: float | None = ...,
    max_iterations: int = ...,
    avoid_gt: bool = ...,
    avoid_cpg: bool = ...,
    track_provenance: bool = ...,
    objective: str | None = ...,
) -> tuple[str, list[str], MutagenesisReport]: ...
def _expand_iupac_site(pattern: str) -> list[str]: ...
def _check_predicates_via_type_system(
    seq: str,
    organism: str,
    is_eukaryote: bool = ...,
    enzyme_list: list[str] | None = ...,
    skip_splice_check: bool = ...,
) -> list[PredicateResult]: ...

__all__: list[str]
