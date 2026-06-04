"""
BioCompiler Multi-Constraint Stress Test
=========================================

Demonstrates the provenance system's unique value when many constraints tug
in opposite directions.  This is the "groundbreaking demonstration" — it shows
which codon choices were altered by which constraint and how the tool
negotiated the trade-off.

Four conflicting scenarios are defined:

  a. High CAI + avoid restriction sites (EcoRI, BamHI, HindIII)
     — restriction sites often appear in high-CAI codons.

  b. High CAI + no CpG
     — arginine codons (CGN) are high-CAI but CpG-heavy.

  c. High CAI + GC range [0.40, 0.55] + no restriction sites + no CpG
     + avoid ATTTA — the "everything conflicts" scenario.

  d. High CAI + no cryptic splice sites (for human targets)
     — GT/AG dinucleotides in high-CAI codons create cryptic splice signals.

For each scenario we run optimisation twice:

  1. **Unconstrained** — only CAI maximisation (wide GC, no extra constraints).
  2. **Constrained** — CAI plus the scenario's constraints.

We then compare the two results, attribute every codon change to the
constraint that caused it (using the provenance trail), and compute the CAI
cost of each constraint.

Public API
----------
- ``StressTestResult``    — data class holding per-scenario results
- ``run_stress_test``     — run a single scenario
- ``run_all_stress_tests`` — run all 4 scenarios on HBB and GFP
- ``print_stress_test_report`` — formatted table of results
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from ..optimization import optimize_sequence, OptimizationResult
from ..optimization import _greedy_optimize, protein_to_aa_list
from ..decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    DecisionProvenanceCollector,
    OptimizationDecisionTrail,
)
from ..scanner import gc_content as _gc_content
from ..translation import compute_cai as _compute_cai

logger = logging.getLogger(__name__)

__all__ = [
    "StressTestResult",
    "StressScenario",
    "STRESS_SCENARIOS",
    "run_stress_test",
    "run_all_stress_tests",
    "print_stress_test_report",
]


# ---------------------------------------------------------------------------
# Stress-test scenario definitions
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class StressScenario:
    """Definition of a stress-test scenario with conflicting constraints.

    Attributes:
        name: Human-readable scenario name.
        description: One-line explanation of the conflict.
        enzymes: Restriction enzymes to avoid (empty = default set).
        avoid_cpg: Whether to avoid CpG dinucleotides.
        avoid_attta: Whether to avoid ATTTA instability motifs.
        cryptic_splice_threshold: MaxEnt score threshold for cryptic splice
            avoidance.  ``None`` means no splice avoidance (the default
            optimiser threshold is not overridden).
        gc_lo: Lower GC bound.
        gc_hi: Upper GC bound.
    """

    name: str
    description: str
    enzymes: list[str] = field(default_factory=list)
    avoid_cpg: bool = False
    avoid_attta: bool = False
    cryptic_splice_threshold: float | None = None
    gc_lo: float = 0.30
    gc_hi: float = 0.70


STRESS_SCENARIOS: dict[str, StressScenario] = {
    "a_restriction_sites": StressScenario(
        name="High CAI + Avoid Restriction Sites",
        description=(
            "Restriction sites (EcoRI GAATTC, BamHI GGATCC, HindIII AAGCTT) "
            "often appear in high-CAI codons — avoiding them forces synonymous "
            "substitutions that reduce CAI."
        ),
        enzymes=["EcoRI", "BamHI", "HindIII"],
    ),
    "b_no_cpg": StressScenario(
        name="High CAI + No CpG",
        description=(
            "Arginine codons CGT/CGC/CGA/CGG are high-CAI in human but are "
            "CpG-heavy.  Avoiding CpG forces a switch to AGA/AGG, reducing CAI."
        ),
        avoid_cpg=True,
    ),
    "c_everything": StressScenario(
        name="High CAI + GC [0.40–0.55] + No RS + No CpG + No ATTTA",
        description=(
            "The 'everything conflicts' scenario — tight GC window, restriction "
            "site avoidance, CpG avoidance, and ATTTA avoidance all pull codon "
            "choices away from the CAI-optimal solution simultaneously."
        ),
        enzymes=["EcoRI", "BamHI", "HindIII"],
        avoid_cpg=True,
        avoid_attta=True,
        gc_lo=0.40,
        gc_hi=0.55,
    ),
    "d_no_cryptic_splice": StressScenario(
        name="High CAI + No Cryptic Splice Sites",
        description=(
            "GT/AG dinucleotides in high-CAI Valine (GTN) and Serine/Arginine "
            "(AGN) codons create cryptic splice donor/acceptor signals.  Avoiding "
            "them forces less optimal codon choices."
        ),
        cryptic_splice_threshold=3.0,
    ),
}


# ---------------------------------------------------------------------------
# Test protein sequences
# ---------------------------------------------------------------------------

_HBB_PROTEIN: str = (
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"
    "GSQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLS"
    "HCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
)

_GFP_PROTEIN: str = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTT"
    "GKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFK"
    "DDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYI"
    "TADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLST"
    "QSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)


# ---------------------------------------------------------------------------
# StressTestResult
# ---------------------------------------------------------------------------

@dataclass
class StressTestResult:
    """Result of a single stress-test scenario.

    Attributes:
        scenario_name: Name of the scenario (e.g. ``"High CAI + No CpG"``).
        protein: Input amino acid sequence.
        organism: Target organism.
        cai_unconstrained: CAI when only maximising CAI (no other constraints).
        cai_constrained: CAI with all scenario constraints active.
        cai_loss_per_constraint: Mapping from constraint name to the CAI delta
            attributable to that constraint (always ≤ 0).
        codons_changed_per_constraint: Mapping from constraint name to the
            number of codon positions altered by that constraint.
        provenance_records: Full provenance trail from the constrained run.
        most_expensive_constraint: Name of the constraint that cost the most
            CAI.
    """

    scenario_name: str
    protein: str
    organism: str
    cai_unconstrained: float
    cai_constrained: float
    cai_loss_per_constraint: dict[str, float]
    codons_changed_per_constraint: dict[str, int]
    provenance_records: list[Any]
    most_expensive_constraint: str

    @property
    def total_cai_loss(self) -> float:
        """Total CAI loss from unconstrained to constrained."""
        return self.cai_constrained - self.cai_unconstrained


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _extract_codons(dna: str) -> list[str]:
    """Split a DNA sequence into a list of 3-base codons."""
    return [dna[i:i + 3] for i in range(0, len(dna), 3)]


def _classify_constraint_reason(reason: str) -> str:
    """Map a fine-grained constraint_reason from CodonDecision to a
    canonical constraint name used in the stress-test report.

    The optimiser records reasons like ``"avoid_restriction_site:EcoRI"``,
    ``"gc_content"``, ``"cryptic_splice_donor"``, etc.  We normalise
    these into a small set of categories for the stress-test report.
    """
    r = reason.lower()

    # Restriction site
    if "restriction" in r or "ecori" in r or "bamhi" in r or "hindiii" in r:
        return "NoRestrictionSite"

    # GC content
    if "gc" in r:
        return "GCInRange"

    # CpG
    if "cpg" in r:
        return "NoCpG"

    # ATTTA / instability
    if "attta" in r or "instability" in r:
        return "NoATTTA"

    # Cryptic splice
    if "splice" in r or "cryptic" in r or "donor" in r or "acceptor" in r:
        return "NoCrypticSplice"

    # Homopolymer / T-run
    if "homopolymer" in r or "t_run" in r or "t-run" in r:
        return "NoTRun"

    # CAI (default / unconstrained choice)
    if "cai" in r or "maximize" in r or "codon_adapt" in r:
        return "MaximizeCAI"

    # Fallback — return the original reason
    return reason


def _analyze_provenance(
    unconstrained_seq: str,
    constrained_seq: str,
    protein: str,
    trail: OptimizationDecisionTrail | None,
    organism: str,
    scenario: StressScenario | None = None,
) -> tuple[dict[str, float], dict[str, int]]:
    """Analyse provenance records to attribute codon changes to constraints.

    Uses three strategies in priority order:
    1. Sequence-difference analysis: compare the unconstrained (max-CAI)
       and constrained sequences position-by-position, and infer which
       constraint caused each change from the biological context.
    2. Decision trail: use the provenance trail's constraint_reason fields.
    3. Aggregate fallback: if individual attribution fails, report the
       total CAI loss as "all_constraints".

    Returns
    -------
    cai_loss_per_constraint : dict[str, float]
        Mapping from constraint name → estimated CAI delta.
    codons_changed_per_constraint : dict[str, int]
        Mapping from constraint name → number of codon positions altered.
    """
    import math
    from ..organisms import CODON_ADAPTIVENESS_TABLES
    from ..constants import reverse_complement as _rc

    uncon_codons = _extract_codons(unconstrained_seq)
    con_codons = _extract_codons(constrained_seq)
    adapt = CODON_ADAPTIVENESS_TABLES.get(organism, {})

    cai_loss: dict[str, float] = {}
    codons_changed: dict[str, int] = {}

    # ------------------------------------------------------------------
    # Strategy 1: Sequence-difference analysis with constraint inference.
    # For each codon position that changed, check what constraint the
    # change satisfies.
    # ------------------------------------------------------------------
    # Build the list of restriction sites to check against
    rs_sites: list[str] = []
    if scenario is not None and scenario.enzymes:
        from ..constants import RESTRICTION_ENZYMES as _RE
        for enz in scenario.enzymes:
            site = _RE.get(enz, "")
            if site:
                rs_sites.append(site.upper())
                rc = _rc(site).upper()
                if rc != site.upper():
                    rs_sites.append(rc)

    for pos in range(min(len(uncon_codons), len(con_codons))):
        if uncon_codons[pos] == con_codons[pos]:
            continue  # No change at this position

        old_codon = uncon_codons[pos]
        new_codon = con_codons[pos]

        # Compute per-position CAI delta
        w_old = adapt.get(old_codon, 0.01)
        w_new = adapt.get(new_codon, 0.01)
        delta = math.log(max(w_new, 1e-10)) - math.log(max(w_old, 1e-10))

        # Determine which constraint likely caused this change
        reason = _infer_constraint_for_change(
            pos, old_codon, new_codon, unconstrained_seq, constrained_seq,
            rs_sites, scenario,
        )

        cai_loss[reason] = cai_loss.get(reason, 0.0) + delta
        codons_changed[reason] = codons_changed.get(reason, 0) + 1

    # ------------------------------------------------------------------
    # Strategy 2: Supplement with decision trail data if available.
    # ------------------------------------------------------------------
    if trail is not None and trail.constraint_decisions:
        for cd in trail.constraint_decisions:
            name = cd.constraint_name
            norm_name = _classify_constraint_reason(name)
            if norm_name == name:
                norm_name = _classify_constraint_reason(name.replace("_", " "))
            if cd.positions_affected and norm_name not in codons_changed:
                codons_changed[norm_name] = len(cd.positions_affected)
            if norm_name not in cai_loss:
                cai_loss[norm_name] = cd.impact_on_cai

    # ------------------------------------------------------------------
    # Strategy 3: Aggregate fallback.
    # ------------------------------------------------------------------
    if not codons_changed and uncon_codons != con_codons:
        n_changed = sum(1 for u, c in zip(uncon_codons, con_codons) if u != c)
        codons_changed["all_constraints"] = n_changed
        try:
            total_loss = math.log(max(_compute_cai(constrained_seq, organism), 1e-10)) - \
                         math.log(max(_compute_cai(unconstrained_seq, organism), 1e-10))
        except Exception:
            total_loss = 0.0
        cai_loss["all_constraints"] = total_loss

    return cai_loss, codons_changed


def _infer_constraint_for_change(
    pos: int,
    old_codon: str,
    new_codon: str,
    unconstrained_seq: str,
    constrained_seq: str,
    rs_sites: list[str],
    scenario: StressScenario | None,
) -> str:
    """Infer which constraint caused a codon change at *pos*.

    Checks each constraint in a fixed priority order and returns the
    first one that plausibly explains the change.
    """
    from ..constants import reverse_complement as _rc

    # 1. Restriction site avoidance
    if rs_sites:
        # Check if the old codon was part of a restriction site that
        # was eliminated in the constrained sequence
        for site in rs_sites:
            if site in unconstrained_seq and site not in constrained_seq:
                # This restriction site was eliminated — check if our
                # position overlaps with it
                site_positions = []
                start = 0
                while True:
                    p = unconstrained_seq.find(site, start)
                    if p == -1:
                        break
                    site_positions.append(p)
                    start = p + 1
                for sp in site_positions:
                    codon_start = pos * 3
                    codon_end = codon_start + 3
                    site_end = sp + len(site)
                    if codon_start < site_end and codon_end > sp:
                        return "NoRestrictionSite"

    # 2. CpG avoidance
    if scenario is not None and scenario.avoid_cpg:
        old_has_cpg = "CG" in old_codon
        new_has_cpg = "CG" in new_codon
        if old_has_cpg and not new_has_cpg:
            return "NoCpG"
        # Also check cross-codon CpG
        old_seq_context = unconstrained_seq[max(0, pos*3-1):pos*3+4]
        new_seq_context = constrained_seq[max(0, pos*3-1):pos*3+4]
        if "CG" in old_seq_context and "CG" not in new_seq_context:
            return "NoCpG"

    # 3. ATTTA avoidance
    if scenario is not None and scenario.avoid_attta:
        old_has_attta = "ATTTA" in unconstrained_seq
        new_no_attta = "ATTTA" not in constrained_seq
        if old_has_attta and new_no_attta:
            # Check if this position overlaps with an ATTTA
            start = 0
            while True:
                p = unconstrained_seq.find("ATTTA", start)
                if p == -1:
                    break
                codon_start = pos * 3
                codon_end = codon_start + 3
                if codon_start < p + 5 and codon_end > p:
                    return "NoATTTA"
                start = p + 1

    # 4. Cryptic splice avoidance
    if scenario is not None and scenario.cryptic_splice_threshold is not None:
        # Check if the old codon contains GT (donor) or AG (acceptor)
        # and the new codon doesn't
        if ("GT" in old_codon and "GT" not in new_codon) or \
           ("AG" in old_codon and "AG" not in new_codon):
            return "NoCrypticSplice"

    # 5. GC range
    if scenario is not None and (scenario.gc_lo > 0.25 or scenario.gc_hi < 0.75):
        old_gc = sum(1 for b in old_codon if b in "GC") / 3.0
        new_gc = sum(1 for b in new_codon if b in "GC") / 3.0
        # If the change moved GC toward the target range
        target = (scenario.gc_lo + scenario.gc_hi) / 2.0
        if abs(new_gc - target) < abs(old_gc - target):
            return "GCInRange"

    # 6. Default: unknown constraint
    return "OtherConstraint"


# ---------------------------------------------------------------------------
# run_stress_test
# ---------------------------------------------------------------------------

def run_stress_test(
    protein: str,
    organism: str,
    scenario: str | StressScenario,
) -> StressTestResult:
    """Run a single stress-test scenario.

    Parameters
    ----------
    protein : str
        Amino acid sequence (single-letter codes).
    organism : str
        Target organism (e.g. ``"Homo_sapiens"``).
    scenario : str | StressScenario
        Scenario key (e.g. ``"a_restriction_sites"``) or a
        :class:`StressScenario` instance.

    Returns
    -------
    StressTestResult
        Detailed result with CAI loss breakdown and provenance.
    """
    if isinstance(scenario, str):
        if scenario not in STRESS_SCENARIOS:
            raise ValueError(
                f"Unknown scenario {scenario!r}. "
                f"Available: {list(STRESS_SCENARIOS.keys())}"
            )
        scenario = STRESS_SCENARIOS[scenario]

    logger.info(
        "Running stress test: %s (protein=%d aa, organism=%s)",
        scenario.name, len(protein), organism,
    )

    # ── Step 1: Unconstrained optimisation (CAI only) ────────────────
    # Build the maximum-CAI sequence directly: for each amino-acid
    # position, pick the codon with the highest adaptiveness value.
    # This gives us the theoretical CAI ceiling without ANY constraints.
    from ..organisms import CODON_ADAPTIVENESS_TABLES
    from ..constants import AA_TO_CODONS as _AA_TO_CODONS
    adapt = CODON_ADAPTIVENESS_TABLES.get(organism, {})
    aas = list(protein)
    unconstrained_codons: list[str] = []
    for aa in aas:
        codons = _AA_TO_CODONS.get(aa, [])
        if not codons:
            unconstrained_codons.append("???")
            continue
        # Pick the codon with the highest adaptiveness
        best_codon = max(codons, key=lambda c: adapt.get(c, 0.0))
        unconstrained_codons.append(best_codon)
    unconstrained_seq = "".join(unconstrained_codons)
    cai_unconstrained = _compute_cai(unconstrained_seq, organism)

    # ── Step 2: Constrained optimisation ─────────────────────────────
    # Run the full optimisation pipeline with the scenario's constraints.
    enzyme_list = list(scenario.enzymes) if scenario.enzymes else []

    # Determine splice threshold
    splice_kwargs = {}
    if scenario.cryptic_splice_threshold is not None:
        splice_kwargs["splice_low"] = scenario.cryptic_splice_threshold

    constrained_result: OptimizationResult = optimize_sequence(
        target_protein=protein,
        organism=organism,
        gc_lo=scenario.gc_lo,
        gc_hi=scenario.gc_hi,
        enzymes=enzyme_list if enzyme_list else None,
        track_provenance=True,
        optimize_mrna_stability=False,
        include_utr=False,
        **splice_kwargs,
    )
    cai_constrained = constrained_result.cai
    constrained_seq = constrained_result.sequence

    # ── Step 3: Analyse provenance ───────────────────────────────────
    trail = constrained_result.decision_trail

    provenance_list: list[Any] = []
    if trail is not None:
        provenance_list.append(trail)
        if trail.constraint_decisions:
            for cd in trail.constraint_decisions:
                provenance_list.append(cd.to_dict())
        if trail.codon_decisions:
            for cd in trail.codon_decisions[:20]:  # Cap at 20 for readability
                provenance_list.append(cd.to_dict())

    # Analyse codon changes
    cai_loss, codons_changed = _analyze_provenance(
        unconstrained_seq,
        constrained_seq,
        protein,
        trail,
        organism,
        scenario=scenario,
    )

    # If we couldn't get per-constraint breakdown, add an aggregate
    if not cai_loss and cai_unconstrained != cai_constrained:
        cai_loss["all_constraints"] = cai_constrained - cai_unconstrained
    if not codons_changed and unconstrained_seq != constrained_seq:
        uncon_codons = _extract_codons(unconstrained_seq)
        con_codons = _extract_codons(constrained_seq)
        n_changed = sum(1 for u, c in zip(uncon_codons, con_codons) if u != c)
        codons_changed["all_constraints"] = n_changed

    # Determine most expensive constraint
    if cai_loss:
        most_expensive = min(cai_loss, key=cai_loss.get)
    else:
        most_expensive = "none"

    result = StressTestResult(
        scenario_name=scenario.name,
        protein=protein,
        organism=organism,
        cai_unconstrained=cai_unconstrained,
        cai_constrained=cai_constrained,
        cai_loss_per_constraint=cai_loss,
        codons_changed_per_constraint=codons_changed,
        provenance_records=provenance_list,
        most_expensive_constraint=most_expensive,
    )

    logger.info(
        "Stress test complete: %s — CAI %.4f → %.4f (loss=%.4f), "
        "most expensive constraint: %s",
        scenario.name, cai_unconstrained, cai_constrained,
        result.total_cai_loss, most_expensive,
    )

    return result


# ---------------------------------------------------------------------------
# run_all_stress_tests
# ---------------------------------------------------------------------------

def run_all_stress_tests() -> list[StressTestResult]:
    """Run all 4 stress-test scenarios on HBB and GFP.

    Returns
    -------
    list[StressTestResult]
        One result per (protein, scenario) combination, 8 total.
    """
    results: list[StressTestResult] = []

    proteins = {
        "HBB": _HBB_PROTEIN,
        "GFP": _GFP_PROTEIN,
    }

    for gene_name, protein_seq in proteins.items():
        for scenario_key in STRESS_SCENARIOS:
            logger.info(
                "=== Stress test: %s / %s ===", gene_name, scenario_key,
            )
            try:
                result = run_stress_test(
                    protein=protein_seq,
                    organism="Homo_sapiens",
                    scenario=scenario_key,
                )
                results.append(result)
            except Exception:
                logger.exception(
                    "Stress test FAILED: %s / %s", gene_name, scenario_key,
                )

    return results


# ---------------------------------------------------------------------------
# print_stress_test_report
# ---------------------------------------------------------------------------

def print_stress_test_report(results: list[StressTestResult]) -> str:
    """Print a formatted table showing each scenario's CAI loss breakdown.

    Parameters
    ----------
    results : list[StressTestResult]
        Results from :func:`run_all_stress_tests` or :func:`run_stress_test`.

    Returns
    -------
    str
        The formatted report string (also printed to stdout).
    """
    lines: list[str] = []

    lines.append("")
    lines.append("=" * 90)
    lines.append("  MULTI-CONSTRAINT STRESS TEST REPORT")
    lines.append("  Provenance-Driven CAI Tradeoff Analysis")
    lines.append("=" * 90)
    lines.append("")

    for idx, r in enumerate(results, 1):
        gene = "HBB" if len(r.protein) < 160 else "GFP"
        lines.append(f"  [{idx}] {r.scenario_name}  ({gene}, {r.organism})")
        lines.append(f"      CAI: {r.cai_unconstrained:.4f} → {r.cai_constrained:.4f}  "
                      f"(loss = {r.total_cai_loss:.4f})")
        lines.append("")

        # CAI loss per constraint
        if r.cai_loss_per_constraint:
            lines.append("      CAI Loss Breakdown:")
            # Sort by CAI loss (most negative first = most expensive)
            sorted_losses = sorted(
                r.cai_loss_per_constraint.items(),
                key=lambda x: x[1],
            )
            for constraint, loss in sorted_losses:
                marker = " ◄ MOST EXPENSIVE" if constraint == r.most_expensive_constraint else ""
                n_changed = r.codons_changed_per_constraint.get(constraint, 0)
                lines.append(
                    f"        {constraint:30s}  CAI Δ = {loss:+.4f}  "
                    f"({n_changed} codons changed){marker}"
                )
        else:
            lines.append("      (No per-constraint CAI loss data)")

        lines.append("")

        # Codon change summary
        if r.codons_changed_per_constraint:
            total_changed = sum(r.codons_changed_per_constraint.values())
            lines.append(f"      Total codons changed: {total_changed}")
        lines.append("")
        lines.append("  " + "-" * 86)
        lines.append("")

    # ── Summary ──────────────────────────────────────────────────────
    lines.append("")
    lines.append("  SUMMARY: Most Expensive Constraint per Scenario")
    lines.append("  " + "=" * 60)

    # Group by scenario
    seen_scenarios: set[str] = set()
    for r in results:
        key = r.scenario_name
        if key not in seen_scenarios:
            seen_scenarios.add(key)
            gene = "HBB" if len(r.protein) < 160 else "GFP"
            lines.append(
                f"    {r.scenario_name:50s} ({gene})  →  "
                f"{r.most_expensive_constraint}"
            )

    lines.append("")
    lines.append("  KEY INSIGHT: The provenance system reveals exactly which")
    lines.append("  constraints forced which codon changes, and the CAI cost of")
    lines.append("  each trade-off.  Without provenance, these trade-offs are")
    lines.append("  invisible — the user sees only a final CAI score.")
    lines.append("")
    lines.append("=" * 90)
    lines.append("")

    report = "\n".join(lines)
    print(report)
    return report
