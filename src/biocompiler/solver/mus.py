"""
BioCompiler CSP Solver — Minimal Unsatisfiable Subset (MUS) Analysis

When the CSP solver reports INFEASIBLE, this module identifies the *minimal*
set of conflicting constraints so that the user gets actionable feedback
instead of a generic "no solution found" message.

Algorithm: deletion-based MUS (also called "quickXplain" variant):
    1. Start with all constraints active.
    2. Tentatively remove each constraint and re-solve.
    3. If the reduced problem is *still* infeasible → the removed constraint
       is NOT in the MUS; keep it out.
    4. If the reduced problem becomes *satisfiable* → the removed constraint
       IS in the MUS; put it back.
    5. Repeat until every constraint has been tested.

This yields an irreducible unsatisfiable core: removing any single constraint
from the MUS makes the problem satisfiable.

Public API
----------
- ``compute_mus``          — extract the MUS via deletion-based algorithm
- ``explain_conflict``     — human-readable conflict explanation
- ``suggest_relaxations``  — concrete constraint relaxation suggestions
- ``quick_feasibility_check`` — fast pre-solve feasibility diagnostics
- ``FeasibilityReport``    — dataclass for quick-check results
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Sequence

from .types import (
    ConstraintSpec,
    ConstraintType,
    CSPModel,
    MUSReport,
    SolverBackendProtocol,
)
from ..constants import AA_TO_CODONS

logger = logging.getLogger(__name__)


__all__ = [
    "FeasibilityReport",
    "compute_mus",
    "explain_conflict",
    "suggest_relaxations",
    "quick_feasibility_check",
]

# ---- Named constants (replace magic numbers) ----
BASES_PER_CODON: int = 3
"""Number of nucleotide bases per codon."""

GC_RELAXATION_WIDTH: float = 0.05
"""Width (in fractional GC) by which to widen GC bounds in relaxation suggestions."""

CRYPTIC_SPLICE_THRESHOLD_INCREMENT: float = 1.0
"""Amount to increase the cryptic splice score threshold in relaxation suggestions."""


# =====================================================================
# FeasibilityReport dataclass
# =====================================================================

@dataclass
class FeasibilityReport:
    """Result of a quick (pre-solve) feasibility check.

    This is much cheaper than running the full MUS algorithm because
    it only inspects the model structure, not the solver search space.

    Attributes:
        feasible: True if no obvious impossibilities were detected.
        warnings: Non-fatal issues that may degrade solution quality.
        impossible_constraints: Constraints that are provably unsatisfiable.
        suggested_relaxations: Concrete relaxation suggestions.
    """

    feasible: bool
    warnings: list[str] = field(default_factory=list)
    impossible_constraints: list[str] = field(default_factory=list)
    suggested_relaxations: list[str] = field(default_factory=list)


# =====================================================================
# quick_feasibility_check
# =====================================================================

def quick_feasibility_check(model: CSPModel) -> FeasibilityReport:
    """Perform fast structural feasibility checks *before* solving.

    This function inspects the CSP model for obvious impossibilities
    without invoking the solver.  It checks:

    1. **GT-dinucleotide feasibility** — if *avoid_gt_dinucleotide* is
       enabled and too many Valine (V) positions force GT at a codon
       boundary, the constraint is infeasible.
    2. **GC bounds achievability** — whether the min/max GC bounds are
       reachable given the available codons at each position.
    3. **Restriction-site codon elimination** — whether any amino acid
       has *zero* codons that avoid a given restriction site.
    4. **CpG feasibility** — whether CpG avoidance eliminates all codons
       for any amino acid.

    Parameters
    ----------
    model:
        The CSP model to inspect.

    Returns
    -------
    FeasibilityReport
        A report with feasibility status, warnings, and suggestions.
    """
    # Validate model structure before inspection
    if not model.protein_sequence:
        logger.warning("quick_feasibility_check called with empty protein sequence")
        return FeasibilityReport(feasible=True)
    warnings: list[str] = []
    impossible: list[str] = []
    suggestions: list[str] = []

    protein = model.protein_sequence
    config = model.config

    # --- 1. GT-dinucleotide feasibility ---
    # Detect from constraint specs rather than config flag
    has_gt_constraint = any(
        c.ctype == ConstraintType.NO_GT_DINUCLEOTIDE for c in model.constraints
    )
    if has_gt_constraint:
        mandatory_gt_boundaries = 0
        for i in range(len(protein) - 1):
            aa_cur = protein[i]
            aa_nxt = protein[i + 1]
            cur_codons = model.codon_domains.get(i, AA_TO_CODONS.get(aa_cur, []))
            nxt_codons = model.codon_domains.get(i + 1, AA_TO_CODONS.get(aa_nxt, []))
            # Check if ALL codon pairs produce GT at the boundary
            all_gt = all(
                c1.endswith("G") and c2.startswith("T")
                for c1 in cur_codons
                for c2 in nxt_codons
            )
            if all_gt:
                mandatory_gt_boundaries += 1

        if mandatory_gt_boundaries > 0:
            impossible.append("NoGTDinucleotide")
            suggestions.append(
                f"Disable GT-dinucleotide avoidance ({mandatory_gt_boundaries} "
                f"codon-boundary GT dinucleotides are unavoidable)"
            )

    # --- 2. GC bounds achievability ---
    gc_min = config.gc_lo
    gc_max = config.gc_hi
    # Compute the absolute min/max GC possible across all positions
    min_gc_count = 0
    max_gc_count = 0
    total_bases = len(protein) * BASES_PER_CODON
    for i, aa in enumerate(protein):
        codons = model.codon_domains.get(i, AA_TO_CODONS.get(aa, []))
        if not codons:
            impossible.append(f"AminoAcid_{aa}_pos{i}")
            suggestions.append(f"No codons available for {aa} at position {i}")
            continue
        gc_counts = [c.count("G") + c.count("C") for c in codons]
        min_gc_count += min(gc_counts)
        max_gc_count += max(gc_counts)

    if total_bases > 0:
        min_possible_gc = min_gc_count / total_bases
        max_possible_gc = max_gc_count / total_bases

        if min_possible_gc > gc_max:
            impossible.append("GC_CONTENT")
            suggestions.append(
                f"Widen GC upper bound from {gc_max:.0%} to at least "
                f"{min_possible_gc:.0%} (minimum achievable GC is "
                f"{min_possible_gc:.0%})"
            )
        elif max_possible_gc < gc_min:
            impossible.append("GC_CONTENT")
            suggestions.append(
                f"Lower GC lower bound from {gc_min:.0%} to at most "
                f"{max_possible_gc:.0%} (maximum achievable GC is "
                f"{max_possible_gc:.0%})"
            )

        # Warn about tight bounds
        achievable_range = max_possible_gc - min_possible_gc
        requested_range = gc_max - gc_min
        if achievable_range < requested_range and achievable_range > 0:
            warnings.append(
                f"GC bounds [{gc_min:.0%}, {gc_max:.0%}] are tight — "
                f"only [{min_possible_gc:.0%}, {max_possible_gc:.0%}] "
                f"is achievable"
            )

    # --- 3. Restriction site codon elimination ---
    for site_seq in config.restriction_sites:
        # Check each amino acid: does the restriction site eliminate ALL codons?
        for aa, codons in AA_TO_CODONS.items():
            surviving = _codons_avoiding_site(codons, site_seq)
            if not surviving:
                impossible.append(f"RestrictionSite_{site_seq}")
                suggestions.append(
                    f"Remove restriction site {site_seq} from avoidance "
                    f"list — amino acid {aa} has no codons avoiding "
                    f"this site"
                )
                break  # one AA is enough to prove infeasibility

    # Also check restriction-site constraints from the model
    rs_constraints = [
        c for c in model.constraints
        if c.ctype == ConstraintType.RESTRICTION_SITE
    ]
    for c in rs_constraints:
        site_seq = c.params.get("site", c.params.get("sequence", ""))
        enzyme = c.params.get("enzyme", site_seq)
        if not site_seq:
            continue
        for aa, codons in AA_TO_CODONS.items():
            surviving = _codons_avoiding_site(codons, site_seq)
            if not surviving:
                impossible.append(f"RestrictionSite_{enzyme}")
                suggestions.append(
                    f"Remove {enzyme} from restriction site list — "
                    f"amino acid {aa} has no codons avoiding site "
                    f"({site_seq})"
                )
                break

    # --- 4. CpG feasibility ---
    has_cpg_constraint = any(
        c.ctype == ConstraintType.NO_CPG for c in model.constraints
    )
    if has_cpg_constraint or config.avoid_cpg:
        for aa, codons in AA_TO_CODONS.items():
            non_cpg = [c for c in codons if "CG" not in c]
            if not non_cpg:
                impossible.append("NoCpG")
                suggestions.append(
                    f"Disable CpG avoidance — amino acid {aa} has no "
                    f"codons without CpG dinucleotide"
                )
                break

    feasible = len(impossible) == 0
    return FeasibilityReport(
        feasible=feasible,
        warnings=warnings,
        impossible_constraints=impossible,
        suggested_relaxations=suggestions,
    )


# =====================================================================
# compute_mus
# =====================================================================

def compute_mus(model: CSPModel, backend: SolverBackendProtocol) -> MUSReport:
    """Compute the Minimal Unsatisfiable Subset of constraints.

    Uses the **deletion-based MUS algorithm**:

    1. Start with all constraints active.
    2. For each constraint, tentatively remove it and re-solve.
    3. If the reduced problem becomes **satisfiable**, the removed
       constraint is part of the MUS — put it back.
    4. If the reduced problem is **still infeasible**, the removed
       constraint is NOT part of the MUS — leave it out.
    5. Continue until every remaining constraint has been tested.

    The result is an irreducible core: every constraint in the MUS
    is necessary for infeasibility, and removing any one makes the
    problem satisfiable.

    Parameters
    ----------
    model:
        The CSP model (must be infeasible with all constraints).
    backend:
        The solver backend to use for satisfiability checks. Must
        implement the ``SolverBackendProtocol`` (i.e. have a ``solve``
        method that accepts a ``CSPModel`` and returns a ``SolverResult``).

    Returns
    -------
    MUSReport
        The minimal unsatisfiable subset with metadata.
    """
    # Validate backend conforms to SolverBackendProtocol at runtime
    if not hasattr(backend, "solve"):
        logger.error(
            "Backend %r does not implement SolverBackendProtocol "
            "(missing 'solve' method)", type(backend).__name__,
        )
        raise TypeError(
            f"Backend {type(backend).__name__!r} does not implement "
            f"SolverBackendProtocol (missing 'solve' method)"
        )

    start_time: float = time.monotonic()
    all_constraints: list[ConstraintSpec] = list(model.constraints)

    # First verify the full model is indeed infeasible
    try:
        full_result = backend.solve(model)
    except Exception:
        logger.exception(
            "Solver backend %r raised an exception during "
            "initial full-model solve in MUS extraction",
            type(backend).__name__,
        )
        raise
    if full_result.solved:
        # Not infeasible — MUS is undefined
        return MUSReport(
            mus_constraints=[],
            all_constraints=all_constraints,
            iterations=1,
            solve_time_seconds=time.monotonic() - start_time,
            explanation="Model is not infeasible; MUS is undefined.",
        )

    # Deletion-based MUS extraction
    remaining: list[ConstraintSpec] = list(all_constraints)
    mus: list[ConstraintSpec] = []
    iterations: int = 1  # already counted the full solve

    for constraint in list(remaining):
        # Build a reduced model without this constraint
        test_constraints: list[ConstraintSpec] = [c for c in remaining if c != constraint]
        test_model = CSPModel(
            protein_sequence=model.protein_sequence,
            codon_domains=model.codon_domains,
            constraints=test_constraints,
            config=model.config,
        )

        try:
            result = backend.solve(test_model)
        except Exception:
            logger.exception(
                "Solver backend %r raised an exception during "
                "MUS iteration %d (testing constraint %r)",
                type(backend).__name__, iterations, constraint.name,
            )
            raise
        iterations += 1

        if result.solved:
            # Removing this constraint makes the problem satisfiable,
            # so it IS part of the MUS.
            mus.append(constraint)
        else:
            # Still infeasible without this constraint — it's not in the MUS.
            # Remove it from remaining to shrink the search space.
            remaining = test_constraints

    # Find conflict positions from MUS constraints
    conflict_positions = _extract_conflict_positions(mus)

    elapsed: float = time.monotonic() - start_time
    logger.info(
        "MUS extraction completed: %d constraints in MUS out of %d total, "
        "%d solver iterations, %.2fs",
        len(mus), len(all_constraints), iterations, elapsed,
    )

    return MUSReport(
        mus_constraints=mus,
        all_constraints=all_constraints,
        iterations=iterations,
        solve_time_seconds=elapsed,
        conflict_positions=conflict_positions,
    )


# =====================================================================
# explain_conflict
# =====================================================================

def explain_conflict(mus_report: MUSReport, model: CSPModel) -> str:
    """Generate a human-readable explanation of the MUS conflict.

    Produces a multi-line string that describes *which* constraints
    conflict and *why*, using domain knowledge about biological
    constraint interactions.

    Parameters
    ----------
    mus_report:
        The MUS analysis result.
    model:
        The original CSP model (for context).

    Returns
    -------
    str
        A human-readable conflict explanation.
    """
    if not mus_report.mus_constraints:
        if mus_report.explanation:
            return mus_report.explanation
        return "No conflict detected — the problem is satisfiable."

    parts: list[str] = []
    parts.append("Minimal Unsatisfiable Subset — conflicting constraints:\n")

    # Group MUS constraints by type for structured output
    by_type: dict[ConstraintType, list[ConstraintSpec]] = {}
    for c in mus_report.mus_constraints:
        by_type.setdefault(c.ctype, []).append(c)

    for ctype, constraints in by_type.items():
        names = [c.name for c in constraints]
        positions = sorted({p for c in constraints for p in c.positions})
        pos_str = f" at positions {positions}" if positions else ""
        parts.append(
            f"  • {ctype.value} ({len(constraints)} constraint(s))"
            f"{pos_str}: {', '.join(names)}"
        )

    # Add domain-specific explanations for known conflict patterns
    types_present = set(by_type.keys())
    explanations = _generate_domain_explanations(types_present, by_type, model)
    if explanations:
        parts.append("\nConflict analysis:")
        for explanation in explanations:
            parts.append(f"  {explanation}")

    if mus_report.conflict_positions:
        parts.append(
            f"\nKey conflict positions: {mus_report.conflict_positions}"
        )

    parts.append(
        f"\n(MUS extracted in {mus_report.iterations} solver calls, "
        f"{mus_report.solve_time_seconds:.2f}s)"
    )

    return "\n".join(parts)


# =====================================================================
# suggest_relaxations
# =====================================================================

def suggest_relaxations(mus_report: MUSReport, model: CSPModel) -> list[str]:
    """Suggest concrete constraint relaxations to resolve the conflict.

    For each constraint type in the MUS, proposes a domain-specific
    relaxation with an estimated impact on solution quality.

    Parameters
    ----------
    mus_report:
        The MUS analysis result.
    model:
        The original CSP model (for context).

    Returns
    -------
    list[str]
        Ordered list of relaxation suggestions (most impactful first).
    """
    if not mus_report.mus_constraints:
        return []

    suggestions: list[str] = []
    types_in_mus = {c.ctype for c in mus_report.mus_constraints}
    config = model.config

    # --- GC content relaxation ---
    if ConstraintType.GC_CONTENT in types_in_mus:
        gc_min, gc_max = config.gc_lo, config.gc_hi
        widened_min = max(0.0, gc_min - GC_RELAXATION_WIDTH)
        widened_max = min(1.0, gc_max + GC_RELAXATION_WIDTH)
        suggestions.append(
            f"Widen GC bounds from [{gc_min:.0%}, {gc_max:.0%}] to "
            f"[{widened_min:.0%}, {widened_max:.0%}] — "
            f"impact: slight deviation from target GC, minimal effect "
            f"on expression"
        )

    # --- Cryptic splice relaxation ---
    if ConstraintType.NO_CRYPTIC_SPLICE in types_in_mus:
        threshold = config.cryptic_splice_threshold
        new_threshold = threshold + CRYPTIC_SPLICE_THRESHOLD_INCREMENT
        suggestions.append(
            f"Increase cryptic splice threshold from {threshold:.1f} to "
            f"{new_threshold:.1f} — impact: allows weaker splice sites, "
            f"moderate risk of aberrant splicing"
        )

    # --- GT dinucleotide relaxation ---
    if ConstraintType.NO_GT_DINUCLEOTIDE in types_in_mus:
        suggestions.append(
            "Disable GT-dinucleotide avoidance — impact: GT dinucleotides "
            "may create cryptic splice donor sites, but many are biologically "
            "silent; review with splice scanner"
        )

    # --- CpG relaxation ---
    if ConstraintType.NO_CPG in types_in_mus:
        suggestions.append(
            "Disable CpG avoidance — impact: CpG sites are potential "
            "methylation targets; in most expression systems this is low-risk"
        )

    # --- Restriction site relaxation ---
    rs_constraints = [
        c for c in mus_report.mus_constraints
        if c.ctype == ConstraintType.RESTRICTION_SITE
    ]
    if rs_constraints:
        # Group by enzyme (encoded in constraint params)
        enzymes: set[str] = set()
        for c in rs_constraints:
            enzyme = c.params.get("enzyme", c.params.get("site", "unknown"))
            enzymes.add(enzyme)
        for enzyme in sorted(enzymes):
            suggestions.append(
                f"Remove {enzyme} from restriction site avoidance list — "
                f"impact: cannot use {enzyme} for downstream cloning; "
                f"consider alternative enzyme or partial digestion"
            )

    # --- Codon usage relaxation ---
    if ConstraintType.CODON_USAGE in types_in_mus:
        suggestions.append(
            "Relax codon usage bias threshold — impact: may use low-frequency "
            "codons, potentially reducing translation speed by 10-20% in "
            "the affected region"
        )

    # --- MHC binding relaxation ---
    if ConstraintType.MHC_BINDING in types_in_mus:
        suggestions.append(
            "Relax MHC binding threshold — impact: increases immunogenicity "
            "risk; suitable when immunogenicity is not a concern (e.g., "
            "in vitro expression)"
        )

    # --- mRNA stability relaxation ---
    if ConstraintType.MRNA_STABILITY in types_in_mus:
        suggestions.append(
            "Relax mRNA stability constraint — impact: may result in less "
            "optimal mRNA folding, potentially reducing expression; "
            "consider accepting a wider dG range"
        )

    # --- Generic fallback for unknown constraint types ---
    known_types = {
        ConstraintType.GC_CONTENT,
        ConstraintType.NO_CRYPTIC_SPLICE,
        ConstraintType.NO_GT_DINUCLEOTIDE,
        ConstraintType.NO_CPG,
        ConstraintType.RESTRICTION_SITE,
        ConstraintType.CODON_USAGE,
        ConstraintType.MHC_BINDING,
        ConstraintType.MRNA_STABILITY,
    }
    unknown = types_in_mus - known_types
    for ctype in unknown:
        constraints = [c for c in mus_report.mus_constraints if c.ctype == ctype]
        for c in constraints:
            suggestions.append(
                f"Consider relaxing or removing constraint '{c.name}' "
                f"(type: {ctype.value}) — review design requirements"
            )

    return suggestions


# =====================================================================
# Internal helpers
# =====================================================================

def _codons_avoiding_site(
    codons: Sequence[str],
    site: str,
) -> list[str]:
    """Return the subset of *codons* that do NOT contain *site*.

    Also checks the reverse complement of the site, since restriction
    enzymes bind both strands.

    Parameters
    ----------
    codons:
        Candidate codons for an amino acid.
    site:
        Restriction enzyme recognition sequence.

    Returns
    -------
    list[str]
        Codons that avoid the site (and its reverse complement).
    """
    from ..constants import reverse_complement

    site_rc: str = reverse_complement(site)
    return [
        c for c in codons
        if site not in c and site_rc not in c
    ]


def _extract_conflict_positions(
    mus_constraints: list[ConstraintSpec],
) -> list[int]:
    """Extract and sort the codon positions involved in the MUS.

    Parameters
    ----------
    mus_constraints:
        Constraints in the minimal unsatisfiable subset.

    Returns
    -------
    list[int]
        Sorted unique positions where conflicts manifest.
    """
    positions: set[int] = set()
    for c in mus_constraints:
        positions.update(c.positions)
    return sorted(positions)


def _generate_domain_explanations(
    types_present: set[ConstraintType],
    by_type: dict[ConstraintType, list[ConstraintSpec]],
    model: CSPModel,
) -> list[str]:
    """Generate domain-specific explanations for known conflict patterns.

    This encodes biological knowledge about which constraint types
    tend to conflict and why.

    Parameters
    ----------
    types_present:
        Constraint types present in the MUS.
    by_type:
        MUS constraints grouped by type.
    model:
        The original CSP model.

    Returns
    -------
    list[str]
        Human-readable explanations.
    """
    explanations: list[str] = []
    config = model.config

    # GC content + NoCrypticSplice conflict
    if (
        ConstraintType.GC_CONTENT in types_present
        and ConstraintType.NO_CRYPTIC_SPLICE in types_present
    ):
        splice_positions = sorted(
            {p for c in by_type[ConstraintType.NO_CRYPTIC_SPLICE] for p in c.positions}
        )
        pos_desc = f" at position(s) {splice_positions}" if splice_positions else ""
        explanations.append(
            "The GC content constraint conflicts with NoCrypticSplice"
            f"{pos_desc}: the only codons that avoid the cryptic splice "
            "donor motif are low-GC codons, but the sequence is already "
            "at the GC minimum."
        )

    # GC content + RestrictionSite conflict
    if (
        ConstraintType.GC_CONTENT in types_present
        and ConstraintType.RESTRICTION_SITE in types_present
    ):
        enzymes = {
            c.params.get("enzyme", c.params.get("site", "unknown"))
            for c in by_type[ConstraintType.RESTRICTION_SITE]
        }
        enzyme_str = ", ".join(sorted(enzymes))
        explanations.append(
            f"The GC content constraint conflicts with restriction site "
            f"avoidance ({enzyme_str}): avoiding these sites forces the "
            f"use of low-GC codons, pushing GC below the minimum."
        )

    # NoGTDinucleotide + NoCrypticSplice conflict
    if (
        ConstraintType.NO_GT_DINUCLEOTIDE in types_present
        and ConstraintType.NO_CRYPTIC_SPLICE in types_present
    ):
        explanations.append(
            "The NoGTDinucleotide constraint conflicts with "
            "NoCrypticSplice: avoiding GT dinucleotides removes the "
            "most common splice donor sequence, but some codon choices "
            "that satisfy NoCrypticSplice still create GT dinucleotides "
            "at codon boundaries."
        )

    # NoCpG + GC_CONTENT conflict
    if (
        ConstraintType.NO_CPG in types_present
        and ConstraintType.GC_CONTENT in types_present
    ):
        explanations.append(
            "CpG avoidance conflicts with the GC content minimum: "
            "CpG dinucleotides contribute to GC%, and eliminating them "
            "reduces achievable GC content."
        )

    # NoCpG + RestrictionSite conflict
    if (
        ConstraintType.NO_CPG in types_present
        and ConstraintType.RESTRICTION_SITE in types_present
    ):
        explanations.append(
            "CpG avoidance and restriction site avoidance both eliminate "
            "high-GC codons, leaving too few options for some positions."
        )

    # CodonUsage + RestrictionSite conflict
    if (
        ConstraintType.CODON_USAGE in types_present
        and ConstraintType.RESTRICTION_SITE in types_present
    ):
        explanations.append(
            "Codon usage optimization conflicts with restriction site "
            "avoidance: the highest-usage codons for some amino acids "
            "contain restriction site subsequences."
        )

    # Generic two-constraint conflict fallback
    if len(types_present) == 2 and not explanations:
        type_names = [t.value for t in types_present]
        explanations.append(
            f"Constraints {type_names[0]} and {type_names[1]} are mutually "
            f"exclusive for the current sequence and parameter settings."
        )

    # More than two types — general infeasibility
    if len(types_present) > 2 and not explanations:
        type_names = [t.value for t in sorted(types_present, key=lambda t: t.value)]
        explanations.append(
            f"Multiple constraints ({', '.join(type_names)}) jointly "
            f"eliminate all codon choices at one or more positions."
        )

    return explanations
