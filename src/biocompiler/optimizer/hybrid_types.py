"""
Shared types and constants for the HybridOptimizer decomposition.

This module contains data classes and constants that are shared between
hybrid_optimizer.py and the extracted submodules, avoiding circular imports.
"""

from __future__ import annotations

from dataclasses import dataclass, field


# ────────────────────────────────────────────────────────────
# GT avoidance CAI cost threshold
# ────────────────────────────────────────────────────────────
# When deciding whether to prefer a GT-free codon over the optimal
# (highest-CAI) codon, this threshold controls how much relative
# adaptiveness loss is acceptable.  If the CAI loss from using a
# GT-free alternative exceeds this threshold, the optimal codon is
# used instead (accepting the GT).  For eukaryotes, in-codon GTs
# from optimal codons are biologically acceptable, so this should
# be a small value — only sacrifice trivial CAI for GT avoidance.
# Lowered from 0.3 (30%) to 0.03 (3%) in v11 to fix the GT/CAI
# tradeoff where yeast insulin CAI was only 0.83 instead of ~0.99.
# Raised back to 0.05 (5%) in v13 — 0.03 was too strict and caused
# unnecessary CAI loss for organisms with GT-rich optimal codons.
GT_CAI_COST_THRESHOLD: float = 0.05


# ────────────────────────────────────────────────────────────
# Constraint violation types and severity scoring
# ────────────────────────────────────────────────────────────

# Severity weights for different constraint types
# Higher = more important to fix first
SEVERITY_WEIGHTS = {
    "restriction_site": 100.0,    # Binary — site present or not
    "stop_codon": 90.0,           # Fatal — creates premature termination
    "gc_out_of_range": 50.0,      # Hard constraint
    "cryptic_splice_donor": 40.0, # Strong eukaryotic constraint
    "cryptic_splice_acceptor": 40.0,
    "avoidable_gt": 35.0,         # GT that can be removed
    "cpg_island": 20.0,           # Soft — methylation risk
    "atttta_motif": 15.0,         # Soft — mRNA instability
    "t_run": 10.0,                # Soft — polymerase slippage
    "blast_match": 80.0,          # Significant — unwanted homology
    "primer_incompatible": 70.0,  # Significant — PCR failure risk
}


@dataclass
class Violation:
    """A single constraint violation with severity score."""
    violation_type: str
    position: int            # Nucleotide position in sequence
    severity: float          # Weighted severity score
    codon_indices: list[int]  # Codon indices involved
    details: str = ""

    def __lt__(self, other: Violation) -> bool:
        """Priority queue ordering: highest severity first."""
        if self.severity != other.severity:
            return self.severity > other.severity
        return self.position < other.position


@dataclass
class HybridResult:
    """Result from the hybrid optimizer."""
    sequence: str
    cai: float
    gc_content: float
    violations_fixed: int = 0
    hill_climb_improvements: int = 0
    iterations_used: int = 0
    phase1_cai: float = 0.0
    phase2_cai: float = 0.0
    phase3_cai: float = 0.0
    phase4_cai: float = 0.0
    cpb_improvements: int = 0
    mean_cpb: float = 0.0
    warnings: list[str] = field(default_factory=list)
    splice_sites_validated: bool = False
    """When True, MaxEntScan splice site validation was already performed
    during optimization (Phase 3 of _optimize_eukaryote_fast).  The
    outer _evaluate_all_predicates call can skip the redundant
    check_no_cryptic_splice scan, eliminating ~7% overhead on the
    eukaryotic path."""
