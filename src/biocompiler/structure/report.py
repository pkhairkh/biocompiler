"""
BioCompiler Structure Report — Protein Assessment Reporting

Comprehensive protein structure visualization and reporting helpers for
Structure quality, stability, solubility, and immunogenicity assessment
results (including MHC binding and epitope prediction).

Features:
- ProteinAssessmentReport dataclass for structured assessment results
- Multi-format output: text, JSON, and standalone HTML
- SVG-based visualizations: pLDDT bar charts and solubility profiles
- Color-coded verdicts matching existing report.py aesthetic
- 5-valued logic for overall verdict computation
- Actionable recommendation generation
"""

from __future__ import annotations

import html
import json
import logging
from dataclasses import dataclass, field, asdict
from typing import Optional

from ..engine_base import EngineTimer
from ..type_system import Verdict, TypeCheckResult

logger = logging.getLogger(__name__)


def _to_dict(obj: object) -> object:
    """Convert a dataclass or dict-like object to a plain dict.

    Handles dataclasses (via asdict), dicts, and falls back to
    obj.__dict__ for other objects. Returns None for None.
    """
    if obj is None:
        return None
    if isinstance(obj, dict):
        return obj
    try:
        # Check if it's a dataclass instance
        if hasattr(obj, "__dataclass_fields__"):
            return asdict(obj)
    except Exception:
        logger.warning("Serialization failed", exc_info=True)
    if hasattr(obj, "__dict__"):
        return dict(obj.__dict__)
    return {"value": str(obj)}


# ────────────────────────────────────────────────────────────
# Data Classes
# ────────────────────────────────────────────────────────────

@dataclass
class ProteinAssessmentReport:
    """Comprehensive protein assessment report combining protein assessment results.

    Aggregates structure quality, stability, solubility, and
    immunogenicity analyses into a single report with an overall
    verdict and actionable recommendations.
    """

    protein: str
    organism: str
    structure_quality: dict | None = None
    stability: dict | None = None
    solubility: dict | None = None
    immunogenicity: dict | None = None
    mhc_binding: dict | None = None
    epitope: dict | None = None
    predicate_results: list[dict] = field(default_factory=list)
    overall_verdict: str = "UNCERTAIN"
    recommendations: list[str] = field(default_factory=list)
    # EngineResult protocol fields
    success: bool = True
    error: str | None = None
    execution_time_s: float = 0.0


# ────────────────────────────────────────────────────────────
# Core Assessment Function
# ────────────────────────────────────────────────────────────

def assess_protein(
    protein: str,
    organism: str = "Homo_sapiens",
    pdb_string: str | None = None,
    run_structure: bool = True,
    run_stability: bool = True,
    run_solubility: bool = True,
    run_immunogenicity: bool = True,
) -> ProteinAssessmentReport:
    """Run all protein assessment analyses on a protein and return a comprehensive report.

    Each analysis is wrapped in try/except so that individual failures
    do not prevent other analyses from completing. Predicate results from
    all analyses are collected, an overall verdict is computed using
    5-valued logic, and actionable recommendations are generated.

    Args:
        protein: Amino acid sequence (single-letter codes).
        organism: Target organism name (default: Homo_sapiens).
        pdb_string: Optional PDB-format structure string. If provided,
            structure quality analysis will use it instead of prediction.
        run_structure: Whether to run structure quality analysis.
        run_stability: Whether to run stability (FoldX) analysis.
        run_solubility: Whether to run solubility (CamSol) analysis.
        run_immunogenicity: Whether to run immunogenicity analysis
            (includes MHC binding and epitope prediction).

    Returns:
        ProteinAssessmentReport with all available results.
    """
    report = ProteinAssessmentReport(
        protein=protein,
        organism=organism,
    )

    predicate_results: list[dict] = []

    with EngineTimer() as timer:
        # --- Structure Quality ---
        if run_structure:
            try:
                from .quality import compute_structure_quality
                sq_result = compute_structure_quality(protein, pdb_string=pdb_string)
                report.structure_quality = _to_dict(sq_result)
                if isinstance(report.structure_quality, dict):
                    pred = report.structure_quality.get("predicate_result")
                    if pred is not None:
                        predicate_results.append(_normalize_predicate(pred))
            except ImportError:
                logger.debug("structure_quality module not available, skipping")
            except Exception as exc:
                logger.warning("Structure quality analysis failed: %s", exc)
                report.success = False
                report.error = str(exc)

        # --- Stability (FoldX) ---
        if run_stability:
            try:
                from ..foldx import empirical_stability
                stab_result = empirical_stability(protein, pdb_string=pdb_string)
                report.stability = _to_dict(stab_result)
                if isinstance(report.stability, dict):
                    pred = report.stability.get("predicate_result")
                    if pred is not None:
                        predicate_results.append(_normalize_predicate(pred))
            except ImportError:
                logger.debug("foldx module not available, skipping")
            except Exception as exc:
                logger.warning("Stability analysis failed: %s", exc)
                report.success = False
                report.error = str(exc)

        # --- Solubility (CamSol) ---
        if run_solubility:
            try:
                from ..camsol import compute_solubility
                sol_result = compute_solubility(protein)
                report.solubility = _to_dict(sol_result)
                if isinstance(report.solubility, dict):
                    pred = report.solubility.get("predicate_result")
                    if pred is not None:
                        predicate_results.append(_normalize_predicate(pred))
            except ImportError:
                logger.debug("camsol module not available, skipping")
            except Exception as exc:
                logger.warning("Solubility analysis failed: %s", exc)
                report.success = False
                report.error = str(exc)

        # --- Immunogenicity (includes MHC binding + epitope) ---
        if run_immunogenicity:
            try:
                from ..immunogenicity import compute_immunogenicity
                imm_result = compute_immunogenicity(protein, organism=organism)
                report.immunogenicity = _to_dict(imm_result)
                if isinstance(report.immunogenicity, dict):
                    pred = report.immunogenicity.get("predicate_result")
                    if pred is not None:
                        predicate_results.append(_normalize_predicate(pred))
            except ImportError:
                logger.debug("immunogenicity module not available, skipping")
            except Exception as exc:
                logger.warning("Immunogenicity analysis failed: %s", exc)
                report.success = False
                report.error = str(exc)

            try:
                from ..immunogenicity import predict_all as predict_mhc_all
                mhc_result = predict_mhc_all(protein)
                report.mhc_binding = _to_dict(mhc_result)
                if isinstance(report.mhc_binding, dict):
                    pred = report.mhc_binding.get("predicate_result")
                    if pred is not None:
                        predicate_results.append(_normalize_predicate(pred))
            except ImportError:
                logger.debug("mhc_binding module not available, skipping")
            except Exception as exc:
                logger.warning("MHC binding prediction failed: %s", exc)
                report.success = False
                report.error = str(exc)

            try:
                from ..immunogenicity import predict_epitopes
                epi_result = predict_epitopes(protein)
                report.epitope = _to_dict(epi_result)
                if isinstance(report.epitope, dict):
                    pred = report.epitope.get("predicate_result")
                    if pred is not None:
                        predicate_results.append(_normalize_predicate(pred))
            except ImportError:
                logger.debug("epitope module not available, skipping")
            except Exception as exc:
                logger.warning("Epitope prediction failed: %s", exc)
                report.success = False
                report.error = str(exc)

        # Store predicate results and compute overall verdict
        report.predicate_results = predicate_results
        report.overall_verdict = compute_overall_verdict(predicate_results)
        report.recommendations = generate_recommendations(report)

    report.execution_time_s = timer.elapsed

    return report


def _normalize_predicate(pred: dict | TypeCheckResult) -> dict:
    """Normalize a predicate result to a plain dict for serialization."""
    if isinstance(pred, TypeCheckResult):
        return {
            "predicate": pred.predicate,
            "verdict": pred.verdict.value,
            "violation": pred.violation,
            "knowledge_gap": pred.knowledge_gap,
        }
    if isinstance(pred, dict):
        # Ensure verdict is a string
        verdict_val = pred.get("verdict", "UNCERTAIN")
        if isinstance(verdict_val, Verdict):
            pred = {**pred, "verdict": verdict_val.value}
        return pred
    return {"predicate": str(pred), "verdict": "UNCERTAIN"}


# ────────────────────────────────────────────────────────────
# Verdict Computation
# ────────────────────────────────────────────────────────────

def compute_overall_verdict(predicate_results: list[dict]) -> str:
    """Compute overall verdict from predicate results using 5-valued logic.

    Combines all predicate verdicts using the weakest-link principle:
    - Any FAIL -> overall FAIL
    - Any LIKELY_FAIL with no FAIL -> LIKELY_FAIL
    - All PASS -> PASS
    - Otherwise -> UNCERTAIN or LIKELY_PASS based on ratio

    Args:
        predicate_results: List of dicts, each with a "verdict" key
            containing one of: PASS, LIKELY_PASS, UNCERTAIN,
            LIKELY_FAIL, FAIL.

    Returns:
        String verdict value.
    """
    if not predicate_results:
        return "UNCERTAIN"

    verdicts = []
    for pred in predicate_results:
        v = pred.get("verdict", "UNCERTAIN")
        if isinstance(v, Verdict):
            v = v.value
        verdicts.append(str(v))

    # Count each verdict type
    counts = {
        "PASS": 0,
        "LIKELY_PASS": 0,
        "UNCERTAIN": 0,
        "LIKELY_FAIL": 0,
        "FAIL": 0,
    }
    for v in verdicts:
        if v in counts:
            counts[v] += 1
        else:
            counts["UNCERTAIN"] += 1

    total = len(verdicts)

    # Any FAIL -> FAIL
    if counts["FAIL"] > 0:
        return "FAIL"

    # Any LIKELY_FAIL with no FAIL -> LIKELY_FAIL
    if counts["LIKELY_FAIL"] > 0:
        return "LIKELY_FAIL"

    # All PASS -> PASS
    if counts["PASS"] == total:
        return "PASS"

    # All PASS or LIKELY_PASS -> LIKELY_PASS
    positive = counts["PASS"] + counts["LIKELY_PASS"]
    if positive == total:
        return "LIKELY_PASS"

    # Majority positive -> LIKELY_PASS, else UNCERTAIN
    if positive > total / 2:
        return "LIKELY_PASS"

    return "UNCERTAIN"


# ────────────────────────────────────────────────────────────
# Recommendation Generation
# ────────────────────────────────────────────────────────────

def generate_recommendations(report: ProteinAssessmentReport) -> list[str]:
    """Generate actionable recommendations based on assessment results.

    Analyzes each section of the report and produces specific,
    actionable suggestions for improving the protein design.

    Args:
        report: ProteinAssessmentReport with analysis results.

    Returns:
        List of recommendation strings.
    """
    recs: list[str] = []

    # --- Structure Quality Recommendations ---
    sq = report.structure_quality
    if sq and isinstance(sq, dict):
        mean_plddt = sq.get("mean_plddt", 0)
        if mean_plddt < 50:
            recs.append(
                "Structure prediction has very low confidence (pLDDT < 50). "
                "Consider using experimental structure determination (X-ray, NMR, cryo-EM) "
                "or homology modeling if a close template is available."
            )
        elif mean_plddt < 70:
            recs.append(
                "Structure prediction confidence is low (pLDDT 50-70). Some regions may "
                "be disordered. Consider disorder prediction tools and verify functional "
                "domains with experimental data."
            )
        elif mean_plddt < 90:
            recs.append(
                "Structure prediction confidence is moderate (pLDDT 70-90). Confidence "
                "is reasonable for most applications, but low-confidence regions should "
                "be validated experimentally if they overlap functional sites."
            )

        # Check for low-confidence regions
        plddt_scores = sq.get("per_residue_plddt", [])
        if plddt_scores:
            low_conf_count = sum(1 for s in plddt_scores if s < 70)
            if low_conf_count > len(plddt_scores) * 0.3:
                recs.append(
                    f"{low_conf_count}/{len(plddt_scores)} residues have pLDDT < 70. "
                    f"Consider truncating disordered termini or using domain boundaries "
                    f"to improve structure quality."
                )

    # --- Stability Recommendations ---
    stab = report.stability
    if stab and isinstance(stab, dict):
        ddg = stab.get("ddg_kcal_mol", 0)
        verdict = stab.get("verdict", "UNCERTAIN")
        if isinstance(verdict, Verdict):
            verdict = verdict.value
        if verdict in ("FAIL", "LIKELY_FAIL") or ddg > 5.0:
            recs.append(
                f"Protein stability is concerning (ΔΔG ≈ {ddg:.1f} kcal/mol). "
                f"Consider introducing stabilizing mutations: proline substitutions in "
                f"loops, disulfide bonds, salt bridges, or surface charge optimization. "
                f"Alternatively, express at lower temperature or use stabilizing buffers."
            )
        elif ddg > 2.0:
            recs.append(
                f"Protein is marginally stable (ΔΔG ≈ {ddg:.1f} kcal/mol). "
                f"Consider adding stabilizing excipients (glycerol, trehalose) "
                f"to formulation, or engineering surface mutations for improved "
                f"thermodynamic stability."
            )

    # --- Solubility Recommendations ---
    sol = report.solubility
    if sol and isinstance(sol, dict):
        score = sol.get("overall_score", 0)
        verdict = sol.get("verdict", "UNCERTAIN")
        if isinstance(verdict, Verdict):
            verdict = verdict.value
        if verdict in ("FAIL", "LIKELY_FAIL") or score < -1.0:
            recs.append(
                f"Protein solubility is low (CamSol score ≈ {score:.2f}). "
                f"Consider: (1) adding solubility tags (MBP, GST, SUMO), "
                f"(2) introducing charged surface mutations (K, E, D) in "
                f"aggregation-prone regions, (3) optimizing buffer conditions "
                f"(pH, ionic strength, arginine), or (4) expressing in "
                f"solubility-enhanced host strains."
            )
        elif score < 0:
            recs.append(
                f"Protein solubility is borderline (CamSol score ≈ {score:.2f}). "
                f"Minor solubility improvements may be achievable by mutating "
                f"aggregation-prone stretches identified in the per-residue profile."
            )

        # Check for aggregation-prone regions
        per_residue = sol.get("per_residue_scores", [])
        if per_residue:
            agg_regions = _find_negative_stretches(per_residue, min_length=5)
            if agg_regions:
                region_strs = [f"{s}-{e}" for s, e in agg_regions[:3]]
                recs.append(
                    f"Aggregation-prone regions detected at residue positions: "
                    f"{', '.join(region_strs)}. Target these for solubility-enhancing "
                    f"mutations or consider selective truncation."
                )

    # --- Immunogenicity Recommendations ---
    imm = report.immunogenicity
    if imm and isinstance(imm, dict):
        verdict = imm.get("verdict", "UNCERTAIN")
        if isinstance(verdict, Verdict):
            verdict = verdict.value
        if verdict in ("FAIL", "LIKELY_FAIL"):
            recs.append(
                "Immunogenicity risk detected. Consider deimmunization strategies: "
                "(1) mutate strong epitope residues to less immunogenic alternatives "
                "(using BLOSUM62 to preserve function), (2) PEGylation to shield "
                "epitopes, or (3) tolerance induction approaches."
            )
        elif verdict == "UNCERTAIN":
            recs.append(
                "Immunogenicity assessment is uncertain. Consider in silico "
                "T-cell epitope screening with multiple HLA alleles and in vitro "
                "validation before clinical development."
            )

    # --- MHC Binding Recommendations ---
    mhc = report.mhc_binding
    if mhc and isinstance(mhc, dict):
        strong_binders = mhc.get("strong_binders", [])
        if strong_binders and len(strong_binders) > 5:
            recs.append(
                f"{len(strong_binders)} strong MHC binders detected. High MHC "
                f"binding density may indicate immunogenicity risk. Consider "
                f"deimmunization of the strongest binders while preserving "
                f"protein function."
            )

    # --- Epitope Recommendations ---
    epi = report.epitope
    if epi and isinstance(epi, dict):
        epitopes = epi.get("epitopes", [])
        if epitopes and len(epitopes) > 3:
            recs.append(
                f"{len(epitopes)} potential B-cell epitopes detected. If "
                f"immunogenicity is a concern, prioritize epitopes with high "
                f"surface accessibility for mutagenesis."
            )

    # --- No results at all ---
    if not any([report.structure_quality, report.stability,
                report.solubility, report.immunogenicity,
                report.mhc_binding, report.epitope]):
        recs.append(
            "No protein analyses completed successfully. Ensure the "
            "required modules (structure_quality, foldx, camsol, "
            "immunogenicity, mhc_binding, epitope) are installed and "
            "configured."
        )

    return recs


def _find_negative_stretches(
    scores: list[float], min_length: int = 5
) -> list[tuple[int, int]]:
    """Find contiguous stretches of negative (aggregation-prone) scores.

    Args:
        scores: Per-residue solubility scores.
        min_length: Minimum stretch length to report.

    Returns:
        List of (start, end) tuples (half-open intervals).
    """
    stretches = []
    start = None
    for i, s in enumerate(scores):
        if s < 0:
            if start is None:
                start = i
        else:
            if start is not None and (i - start) >= min_length:
                stretches.append((start, i))
            start = None
    if start is not None and (len(scores) - start) >= min_length:
        stretches.append((start, len(scores)))
    return stretches


# ────────────────────────────────────────────────────────────
# Text Formatting
# ────────────────────────────────────────────────────────────

def format_assessment_text(report: ProteinAssessmentReport) -> str:
    """Format a ProteinAssessmentReport as human-readable text.

    Sections: Structure Quality, Stability, Solubility, Immunogenicity,
    Predicates, Overall Verdict, Recommendations. Each section shows
    key metrics and verdict with clear borders and indentation.

    Args:
        report: ProteinAssessmentReport to format.

    Returns:
        Formatted text string.
    """
    lines: list[str] = []
    sep = "=" * 60
    thin_sep = "-" * 60

    lines.append(sep)
    lines.append("  BioCompiler Protein Assessment Report")
    lines.append(sep)
    lines.append(f"  Protein : {report.protein[:60]}{'...' if len(report.protein) > 60 else ''}")
    lines.append(f"  Organism: {report.organism}")
    lines.append("")

    # --- Structure Quality ---
    lines.append(thin_sep)
    lines.append("  STRUCTURE QUALITY")
    lines.append(thin_sep)
    sq = report.structure_quality
    if sq and isinstance(sq, dict):
        mean_plddt = sq.get("mean_plddt", "N/A")
        if isinstance(mean_plddt, float):
            lines.append(f"  Mean pLDDT       : {mean_plddt:.1f}")
        else:
            lines.append(f"  Mean pLDDT       : {mean_plddt}")
        lines.append(f"  Method           : {sq.get('method', 'N/A')}")
        lines.append(f"  Model confidence : {_plddt_label(mean_plddt)}")
        _format_verdict_line(lines, sq.get("verdict"))
    else:
        lines.append("  (not available)")

    # --- Stability ---
    lines.append("")
    lines.append(thin_sep)
    lines.append("  STABILITY (FoldX)")
    lines.append(thin_sep)
    stab = report.stability
    if stab and isinstance(stab, dict):
        ddg = stab.get("ddg_kcal_mol", "N/A")
        if isinstance(ddg, float):
            lines.append(f"  ΔΔG              : {ddg:.2f} kcal/mol")
        else:
            lines.append(f"  ΔΔG              : {ddg}")
        lines.append(f"  Stability class  : {_stability_label(ddg)}")
        _format_verdict_line(lines, stab.get("verdict"))
    else:
        lines.append("  (not available)")

    # --- Solubility ---
    lines.append("")
    lines.append(thin_sep)
    lines.append("  SOLUBILITY (CamSol)")
    lines.append(thin_sep)
    sol = report.solubility
    if sol and isinstance(sol, dict):
        score = sol.get("overall_score", "N/A")
        if isinstance(score, float):
            lines.append(f"  Overall score    : {score:.3f}")
        else:
            lines.append(f"  Overall score    : {score}")
        lines.append(f"  Solubility class : {_solubility_label(score)}")
        per_residue = sol.get("per_residue_scores", [])
        if per_residue:
            agg_count = sum(1 for s in per_residue if s < 0)
            lines.append(f"  Aggregation-prone residues: {agg_count}/{len(per_residue)}")
        _format_verdict_line(lines, sol.get("verdict"))
    else:
        lines.append("  (not available)")

    # --- Immunogenicity ---
    lines.append("")
    lines.append(thin_sep)
    lines.append("  IMMUNOGENICITY")
    lines.append(thin_sep)
    imm = report.immunogenicity
    if imm and isinstance(imm, dict):
        lines.append(f"  Risk score       : {imm.get('risk_score', 'N/A')}")
        lines.append(f"  Risk level       : {imm.get('risk_level', 'N/A')}")
        _format_verdict_line(lines, imm.get("verdict"))
    else:
        lines.append("  (not available)")

    # MHC Binding
    mhc = report.mhc_binding
    if mhc and isinstance(mhc, dict):
        lines.append("")
        lines.append("  MHC Binding")
        strong = mhc.get("strong_binders", [])
        weak = mhc.get("weak_binders", [])
        if isinstance(strong, list):
            lines.append(f"  Strong binders   : {len(strong)}")
        if isinstance(weak, list):
            lines.append(f"  Weak binders     : {len(weak)}")

    # Epitope
    epi = report.epitope
    if epi and isinstance(epi, dict):
        lines.append("")
        lines.append("  Epitope Prediction")
        epitopes = epi.get("epitopes", [])
        if isinstance(epitopes, list):
            lines.append(f"  Predicted epitopes: {len(epitopes)}")

    # --- Predicate Results ---
    lines.append("")
    lines.append(thin_sep)
    lines.append("  PREDICATE RESULTS")
    lines.append(thin_sep)
    if report.predicate_results:
        for i, pred in enumerate(report.predicate_results, 1):
            name = pred.get("predicate", f"Predicate_{i}")
            verdict = pred.get("verdict", "N/A")
            if isinstance(verdict, Verdict):
                verdict = verdict.value
            violation = pred.get("violation", "")
            lines.append(f"  {i:2d}. {name:<30s} [{verdict}]")
            if violation:
                lines.append(f"      {violation}")
    else:
        lines.append("  (no predicate results)")

    # --- Overall Verdict ---
    lines.append("")
    lines.append(sep)
    lines.append(f"  OVERALL VERDICT: {report.overall_verdict}")
    lines.append(sep)

    # --- Recommendations ---
    if report.recommendations:
        lines.append("")
        lines.append(thin_sep)
        lines.append("  RECOMMENDATIONS")
        lines.append(thin_sep)
        for i, rec in enumerate(report.recommendations, 1):
            # Wrap long lines
            words = rec.split()
            current_line = f"  {i}. "
            for word in words:
                if len(current_line) + len(word) + 1 > 72:
                    lines.append(current_line)
                    current_line = "     " + word
                else:
                    current_line += (" " if not current_line.endswith(" ") else "") + word
            if current_line.strip():
                lines.append(current_line)

    lines.append("")
    return "\n".join(lines)


def _format_verdict_line(lines: list[str], verdict: object) -> None:
    """Append a formatted verdict line."""
    if verdict is not None:
        if isinstance(verdict, Verdict):
            v_str = verdict.value
        else:
            v_str = str(verdict)
        symbol = {
            "PASS": "✓",
            "LIKELY_PASS": "✓~",
            "UNCERTAIN": "?",
            "LIKELY_FAIL": "✗~",
            "FAIL": "✗",
        }.get(v_str, "?")
        lines.append(f"  Verdict          : {symbol} {v_str}")


def _plddt_label(mean_plddt: object) -> str:
    """Return a confidence label for a mean pLDDT value."""
    if not isinstance(mean_plddt, (int, float)):
        return "N/A"
    if mean_plddt > 90:
        return "Very high (experimental quality)"
    if mean_plddt > 70:
        return "High (confident backbone)"
    if mean_plddt > 50:
        return "Low (tentative model)"
    return "Very low (likely disordered)"


def _stability_label(ddg: object) -> str:
    """Return a stability class label for ΔΔG."""
    if not isinstance(ddg, (int, float)):
        return "N/A"
    if ddg < 0:
        return "Stabilizing"
    if ddg < 2:
        return "Neutral"
    if ddg < 5:
        return "Marginally unstable"
    return "Unstable"


def _solubility_label(score: object) -> str:
    """Return a solubility class label for CamSol score."""
    if not isinstance(score, (int, float)):
        return "N/A"
    if score > 1.0:
        return "Highly soluble"
    if score > 0:
        return "Soluble"
    if score > -1.0:
        return "Borderline"
    return "Aggregation-prone"


# ────────────────────────────────────────────────────────────
# JSON Formatting
# ────────────────────────────────────────────────────────────

def format_assessment_json(report: ProteinAssessmentReport) -> str:
    """Format a ProteinAssessmentReport as a JSON string.

    All values are JSON-serializable. Dataclass objects are converted
    to dicts. Verdict enums are converted to their string values.

    Args:
        report: ProteinAssessmentReport to format.

    Returns:
        JSON string with indented formatting.
    """
    data = _report_to_dict(report)
    return json.dumps(data, indent=2, default=_json_serializer)


def _report_to_dict(report: ProteinAssessmentReport) -> dict:
    """Convert a ProteinAssessmentReport to a JSON-serializable dict."""
    def _convert_verdicts(obj: object) -> object:
        """Recursively convert Verdict enums to strings."""
        if isinstance(obj, Verdict):
            return obj.value
        if isinstance(obj, dict):
            return {k: _convert_verdicts(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_convert_verdicts(item) for item in obj]
        return obj

    result = {
        "protein": report.protein,
        "organism": report.organism,
        "structure_quality": _convert_verdicts(report.structure_quality),
        "stability": _convert_verdicts(report.stability),
        "solubility": _convert_verdicts(report.solubility),
        "immunogenicity": _convert_verdicts(report.immunogenicity),
        "mhc_binding": _convert_verdicts(report.mhc_binding),
        "epitope": _convert_verdicts(report.epitope),
        "predicate_results": _convert_verdicts(report.predicate_results),
        "overall_verdict": report.overall_verdict,
        "recommendations": report.recommendations,
        "success": report.success,
        "error": report.error,
        "execution_time_s": report.execution_time_s,
    }
    return result


def _json_serializer(obj: object) -> object:
    """Custom JSON serializer for non-standard types."""
    if isinstance(obj, Verdict):
        return obj.value
    if hasattr(obj, "__dict__"):
        return obj.__dict__
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


# ────────────────────────────────────────────────────────────
# HTML Formatting
# ────────────────────────────────────────────────────────────

def format_assessment_html(report: ProteinAssessmentReport) -> str:
    """Format a ProteinAssessmentReport as a standalone HTML page.

    Styled with CSS matching the existing report.py aesthetic.
    Includes expandable sections, color-coded verdicts, a pLDDT
    bar chart, solubility gauge, MHC binding heatmap, and
    recommendations as a styled list.

    Args:
        report: ProteinAssessmentReport to format.

    Returns:
        Self-contained HTML string.
    """
    from datetime import datetime, timezone
    from biocompiler import __version__

    verdict_css = report.overall_verdict.lower()
    verdict_symbols = {
        "PASS": "&#10003;",
        "LIKELY_PASS": "&#10003;~",
        "UNCERTAIN": "?",
        "LIKELY_FAIL": "&#10007;~",
        "FAIL": "&#10007;",
    }
    verdict_symbol = verdict_symbols.get(report.overall_verdict, "?")

    # Build sections
    sq_section = _html_structure_quality_section(report)
    stab_section = _html_stability_section(report)
    sol_section = _html_solubility_section(report)
    imm_section = _html_immunogenicity_section(report)
    pred_section = _html_predicate_section(report)
    rec_section = _html_recommendations_section(report)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>BioCompiler Protein Assessment{f' — {html.escape(report.protein[:30])}' if report.protein else ''}</title>
<style>
:root {{
    --pass: #16a34a;
    --fail: #dc2626;
    --uncertain: #d97706;
    --likely-pass: #22c55e;
    --likely-fail: #f87171;
    --bg: #ffffff;
    --surface: #f8fafc;
    --border: #e2e8f0;
    --text: #1e293b;
    --text-secondary: #64748b;
    --accent: #2563eb;
}}
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    background: var(--surface);
    color: var(--text);
    line-height: 1.6;
    padding: 2rem;
    max-width: 1200px;
    margin: 0 auto;
}}
.header {{
    background: linear-gradient(135deg, #1e293b, #334155);
    color: white;
    padding: 2rem;
    border-radius: 12px;
    margin-bottom: 2rem;
    display: flex;
    justify-content: space-between;
    align-items: center;
}}
.header h1 {{ font-size: 1.75rem; font-weight: 700; }}
.header .subtitle {{ color: #94a3b8; margin-top: 0.25rem; }}
.verdict-badge {{
    font-size: 1.5rem;
    font-weight: 800;
    padding: 0.75rem 2rem;
    border-radius: 8px;
    text-transform: uppercase;
    letter-spacing: 0.05em;
}}
.verdict-badge.pass {{ background: var(--pass); color: white; }}
.verdict-badge.likely_pass {{ background: var(--likely-pass); color: white; }}
.verdict-badge.fail {{ background: var(--fail); color: white; }}
.verdict-badge.likely_fail {{ background: var(--likely-fail); color: white; }}
.verdict-badge.uncertain {{ background: var(--uncertain); color: white; }}
.section {{
    background: var(--bg);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.5rem;
    margin-bottom: 1.5rem;
}}
.section h2 {{
    font-size: 1.25rem;
    font-weight: 600;
    margin-bottom: 1rem;
    padding-bottom: 0.5rem;
    border-bottom: 2px solid var(--border);
    cursor: pointer;
    user-select: none;
    position: relative;
    padding-right: 2rem;
}}
.section h2::after {{
    content: '\\2212';
    position: absolute;
    right: 0;
    top: 50%;
    transform: translateY(-50%);
    font-size: 1.5rem;
    font-weight: 300;
    color: var(--text-secondary);
}}
.section h2.collapsed::after {{ content: '+'; }}
.section-content {{
    overflow: hidden;
    transition: max-height 0.3s ease;
    max-height: 5000px;
}}
.section-content.collapsed {{ max-height: 0; }}
.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
    gap: 1rem;
}}
.stat-card {{
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1rem;
    text-align: center;
}}
.stat-value {{
    font-size: 1.75rem;
    font-weight: 700;
    color: var(--accent);
}}
.stat-label {{
    font-size: 0.8rem;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.05em;
    margin-top: 0.25rem;
}}
.badge {{
    display: inline-block;
    padding: 0.15rem 0.5rem;
    border-radius: 4px;
    font-weight: 600;
    font-size: 0.9rem;
}}
.badge.pass {{ background: #dcfce7; color: var(--pass); }}
.badge.likely_pass {{ background: #dcfce7; color: #16a34a; }}
.badge.fail {{ background: #fee2e2; color: var(--fail); }}
.badge.likely_fail {{ background: #fee2e2; color: #dc2626; }}
.badge.uncertain {{ background: #fef3c7; color: var(--uncertain); }}
table {{
    width: 100%;
    border-collapse: collapse;
    margin-top: 0.5rem;
}}
th, td {{
    padding: 0.5rem 0.75rem;
    text-align: left;
    border-bottom: 1px solid var(--border);
}}
th {{
    background: var(--surface);
    font-weight: 600;
    font-size: 0.85rem;
    text-transform: uppercase;
    letter-spacing: 0.03em;
    color: var(--text-secondary);
}}
tr.pass {{ background: #f0fdf4; }}
tr.likely_pass {{ background: #f0fdf4; }}
tr.fail {{ background: #fef2f2; }}
tr.likely_fail {{ background: #fef2f2; }}
tr.uncertain {{ background: #fffbeb; }}
.svg-container {{
    overflow-x: auto;
    margin: 1rem 0;
}}
.svg-container svg {{ max-width: 100%; height: auto; }}
.mono {{ font-family: 'SF Mono', 'Fira Code', monospace; font-size: 0.9rem; }}
.rec-list {{
    list-style: none;
    padding: 0;
}}
.rec-list li {{
    padding: 0.75rem 1rem;
    margin-bottom: 0.5rem;
    background: var(--surface);
    border-left: 4px solid var(--accent);
    border-radius: 0 6px 6px 0;
    font-size: 0.9rem;
    line-height: 1.5;
}}
.rec-list li.warning {{ border-left-color: var(--uncertain); }}
.rec-list li.error {{ border-left-color: var(--fail); }}
.gauge-container {{
    position: relative;
    width: 200px;
    height: 100px;
    margin: 1rem auto;
}}
.gauge-bg {{
    width: 200px;
    height: 100px;
    border-radius: 100px 100px 0 0;
    background: linear-gradient(to right, var(--fail), var(--uncertain), var(--pass));
    position: relative;
    overflow: hidden;
}}
.gauge-fill {{
    position: absolute;
    bottom: 0;
    left: 0;
    width: 50%;
    height: 100%;
    background: var(--bg);
    transform-origin: left bottom;
}}
.gauge-needle {{
    position: absolute;
    bottom: 0;
    left: 50%;
    width: 2px;
    height: 90px;
    background: var(--text);
    transform-origin: bottom center;
}}
.gauge-label {{
    text-align: center;
    font-size: 1.25rem;
    font-weight: 700;
    margin-top: 0.5rem;
}}
.footer {{
    text-align: center;
    color: var(--text-secondary);
    font-size: 0.85rem;
    margin-top: 2rem;
    padding-top: 1rem;
    border-top: 1px solid var(--border);
}}
/* Dark mode */
.toggle-btn {{
    background: rgba(255,255,255,0.15);
    border: 1px solid rgba(255,255,255,0.3);
    color: white;
    font-size: 1.25rem;
    padding: 0.35rem 0.65rem;
    border-radius: 6px;
    cursor: pointer;
    transition: all 0.2s ease;
    line-height: 1;
}}
.toggle-btn:hover {{ background: rgba(255,255,255,0.25); }}
body.dark-mode {{
    --bg: #1e293b;
    --surface: #0f172a;
    --border: #334155;
    --text: #f1f5f9;
    --text-secondary: #94a3b8;
    --accent: #60a5fa;
}}
body.dark-mode .header {{ background: linear-gradient(135deg, #0f172a, #1e293b); }}
body.dark-mode tr.pass {{ background: #052e16; }}
body.dark-mode tr.likely_pass {{ background: #052e16; }}
body.dark-mode tr.fail {{ background: #450a0a; }}
body.dark-mode tr.likely_fail {{ background: #450a0a; }}
body.dark-mode tr.uncertain {{ background: #451a03; }}
body.dark-mode .badge.pass {{ background: #052e16; color: #4ade80; }}
body.dark-mode .badge.likely_pass {{ background: #052e16; color: #4ade80; }}
body.dark-mode .badge.fail {{ background: #450a0a; color: #f87171; }}
body.dark-mode .badge.likely_fail {{ background: #450a0a; color: #f87171; }}
body.dark-mode .badge.uncertain {{ background: #451a03; color: #fbbf24; }}
body.dark-mode .stat-card {{ background: #0f172a; border-color: #334155; }}
body.dark-mode .rec-list li {{ background: #0f172a; }}
@media print {{
    body {{ padding: 0; }}
    .section {{ break-inside: avoid; }}
    .toggle-btn {{ display: none; }}
    .section-content {{ max-height: none !important; }}
    .section h2 {{ cursor: default; }}
}}
</style>
</head>
<body>

<div class="header">
    <div>
        <h1>BioCompiler Protein Assessment</h1>
        <div class="subtitle">
            {html.escape(report.protein[:60])}{'...' if len(report.protein) > 60 else ''} &mdash;
            {html.escape(report.organism)} &mdash;
            {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}
        </div>
    </div>
    <div style="display:flex;align-items:center;gap:1rem;">
        <button id="dark-mode-toggle" class="toggle-btn" title="Toggle dark mode">&#9790;</button>
        <div class="verdict-badge {verdict_css}">{verdict_symbol} {report.overall_verdict}</div>
    </div>
</div>

<!-- Summary Statistics -->
<section class="section">
    <h2>Assessment Summary</h2>
    <div class="section-content">
    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-value">{len(report.protein)}</div>
            <div class="stat-label">Protein Length (aa)</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{len(report.predicate_results)}</div>
            <div class="stat-label">Predicates Evaluated</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{sum(1 for p in report.predicate_results if p.get('verdict') == 'PASS')}</div>
            <div class="stat-label">Predicates Passed</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{len(report.recommendations)}</div>
            <div class="stat-label">Recommendations</div>
        </div>
    </div>
    </div>
</section>

{sq_section}
{stab_section}
{sol_section}
{imm_section}
{pred_section}
{rec_section}

<div class="footer">
    Generated by BioCompiler v{__version__} &mdash; Protein Structure Assessment
</div>

<script>
(function() {{
    'use strict';

    // === Collapsible Sections ===
    var sections = document.querySelectorAll('.section');
    for (var i = 0; i < sections.length; i++) {{
        var section = sections[i];
        var h2 = section.querySelector('h2');
        if (!h2) continue;
        var wrapper = document.createElement('div');
        wrapper.className = 'section-content';
        var sibling = h2.nextSibling;
        while (sibling) {{
            var next = sibling.nextSibling;
            wrapper.appendChild(sibling);
            sibling = next;
        }}
        section.appendChild(wrapper);
        (function(header, content) {{
            header.addEventListener('click', function() {{
                header.classList.toggle('collapsed');
                content.classList.toggle('collapsed');
            }});
        }})(h2, wrapper);
    }}

    // === Dark Mode Toggle ===
    var darkToggle = document.getElementById('dark-mode-toggle');
    if (darkToggle) {{
        var saved = localStorage.getItem('biocompiler-dark-mode');
        if (saved === 'true') {{
            document.body.classList.add('dark-mode');
            darkToggle.innerHTML = '&#9788;';
        }}
        darkToggle.addEventListener('click', function() {{
            document.body.classList.toggle('dark-mode');
            var isDark = document.body.classList.contains('dark-mode');
            darkToggle.innerHTML = isDark ? '&#9788;' : '&#9790;';
            try {{ localStorage.setItem('biocompiler-dark-mode', isDark); }} catch(e) {{}}
        }});
    }}
}})();
</script>

</body>
</html>"""


def _html_structure_quality_section(report: ProteinAssessmentReport) -> str:
    """Build the HTML section for structure quality results."""
    sq = report.structure_quality
    if not sq or not isinstance(sq, dict):
        return ""

    mean_plddt = sq.get("mean_plddt", "N/A")
    method = sq.get("method", "prediction")
    verdict = sq.get("verdict", "UNCERTAIN")
    if isinstance(verdict, Verdict):
        verdict = verdict.value
    v_css = verdict.lower()

    # pLDDT bar chart
    plddt_scores = sq.get("per_residue_plddt", [])
    plddt_svg = ""
    if plddt_scores:
        plddt_svg = f"""<div class="svg-container">
        {plot_plddt_bar_svg(plddt_scores, report.protein[:30])}
    </div>"""

    return f"""
<section class="section">
    <h2>Structure Quality</h2>
    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-value">{mean_plddt:.1f}</div>
            <div class="stat-label">Mean pLDDT</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{html.escape(str(method))}</div>
            <div class="stat-label">Method</div>
        </div>
        <div class="stat-card">
            <div class="stat-value"><span class="badge {v_css}">{verdict}</span></div>
            <div class="stat-label">Verdict</div>
        </div>
    </div>
    {plddt_svg}
</section>"""


def _html_stability_section(report: ProteinAssessmentReport) -> str:
    """Build the HTML section for stability results."""
    stab = report.stability
    if not stab or not isinstance(stab, dict):
        return ""

    ddg = stab.get("ddg_kcal_mol", "N/A")
    verdict = stab.get("verdict", "UNCERTAIN")
    if isinstance(verdict, Verdict):
        verdict = verdict.value
    v_css = verdict.lower()

    ddg_display = f"{ddg:.2f}" if isinstance(ddg, (int, float)) else str(ddg)

    return f"""
<section class="section">
    <h2>Stability (FoldX)</h2>
    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-value">{ddg_display}</div>
            <div class="stat-label">&Delta;&Delta;G (kcal/mol)</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{_stability_label(ddg)}</div>
            <div class="stat-label">Stability Class</div>
        </div>
        <div class="stat-card">
            <div class="stat-value"><span class="badge {v_css}">{verdict}</span></div>
            <div class="stat-label">Verdict</div>
        </div>
    </div>
</section>"""


def _html_solubility_section(report: ProteinAssessmentReport) -> str:
    """Build the HTML section for solubility results."""
    sol = report.solubility
    if not sol or not isinstance(sol, dict):
        return ""

    score = sol.get("overall_score", "N/A")
    verdict = sol.get("verdict", "UNCERTAIN")
    if isinstance(verdict, Verdict):
        verdict = verdict.value
    v_css = verdict.lower()

    score_display = f"{score:.3f}" if isinstance(score, (int, float)) else str(score)

    # Solubility gauge using CSS
    gauge_html = ""
    if isinstance(score, (int, float)):
        # Map score from [-3, 3] to [0, 100] percentage
        pct = max(0, min(100, (score + 3) / 6 * 100))
        needle_angle = -90 + pct * 1.8  # -90 to +90 degrees
        gauge_color = "#16a34a" if score > 0 else "#dc2626" if score < -1 else "#d97706"
        gauge_html = f"""
    <div class="gauge-container">
        <svg width="200" height="110" viewBox="0 0 200 110">
            <path d="M 10 100 A 90 90 0 0 1 190 100" fill="none" stroke="#e2e8f0" stroke-width="16" stroke-linecap="round"/>
            <path d="M 10 100 A 90 90 0 0 1 190 100" fill="none" stroke="{gauge_color}" stroke-width="16" stroke-linecap="round"
                  stroke-dasharray="{pct * 2.83} 283"/>
            <line x1="100" y1="100" x2="{100 + 80 * _cos_deg(needle_angle):.1f}" y2="{100 - 80 * _sin_deg(needle_angle):.1f}"
                  stroke="{gauge_color}" stroke-width="3" stroke-linecap="round"/>
            <circle cx="100" cy="100" r="5" fill="{gauge_color}"/>
        </svg>
        <div class="gauge-label" style="color:{gauge_color}">{score_display}</div>
    </div>"""

    # Solubility profile SVG
    per_residue = sol.get("per_residue_scores", [])
    sol_svg = ""
    if per_residue:
        sol_svg = f"""<div class="svg-container">
        {plot_solubility_profile_svg(per_residue, report.protein[:30])}
    </div>"""

    return f"""
<section class="section">
    <h2>Solubility (CamSol)</h2>
    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-value">{score_display}</div>
            <div class="stat-label">Overall Score</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{_solubility_label(score)}</div>
            <div class="stat-label">Solubility Class</div>
        </div>
        <div class="stat-card">
            <div class="stat-value"><span class="badge {v_css}">{verdict}</span></div>
            <div class="stat-label">Verdict</div>
        </div>
    </div>
    {gauge_html}
    {sol_svg}
</section>"""


def _html_immunogenicity_section(report: ProteinAssessmentReport) -> str:
    """Build the HTML section for immunogenicity results."""
    imm = report.immunogenicity
    mhc = report.mhc_binding
    epi = report.epitope

    if not any([imm, mhc, epi]):
        return ""

    parts = []

    if imm and isinstance(imm, dict):
        verdict = imm.get("verdict", "UNCERTAIN")
        if isinstance(verdict, Verdict):
            verdict = verdict.value
        v_css = verdict.lower()
        risk = imm.get("risk_score", "N/A")
        risk_display = f"{risk:.2f}" if isinstance(risk, (int, float)) else str(risk)
        level = imm.get("risk_level", "N/A")
        parts.append(f"""
    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-value">{risk_display}</div>
            <div class="stat-label">Immunogenicity Risk Score</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{html.escape(str(level))}</div>
            <div class="stat-label">Risk Level</div>
        </div>
        <div class="stat-card">
            <div class="stat-value"><span class="badge {v_css}">{verdict}</span></div>
            <div class="stat-label">Verdict</div>
        </div>
    </div>""")

    # MHC binding heatmap as HTML table
    if mhc and isinstance(mhc, dict):
        strong = mhc.get("strong_binders", [])
        weak = mhc.get("weak_binders", [])
        mhc_table = ""
        if isinstance(strong, list) and strong:
            rows = ""
            for i, binder in enumerate(strong[:20]):
                if isinstance(binder, dict):
                    allele = binder.get("allele", "N/A")
                    peptide = binder.get("peptide", "N/A")
                    score_val = binder.get("score", 0)
                    pct = min(100, max(0, score_val * 100))
                    bg_color = f"rgba(220, 38, 38, {min(0.8, pct / 100)})"
                    rows += f"""<tr>
                        <td>{html.escape(str(allele))}</td>
                        <td class="mono">{html.escape(str(peptide))}</td>
                        <td style="background-color:{bg_color};color:white;font-weight:600">{score_val:.2f}</td>
                    </tr>"""
                else:
                    rows += f"""<tr><td colspan="3">{html.escape(str(binder))}</td></tr>"""
            mhc_table = f"""
    <h3>MHC Binding (Strong Binders)</h3>
    <table>
        <thead><tr><th>Allele</th><th>Peptide</th><th>Score</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>"""
        elif isinstance(strong, list):
            mhc_table = f"""
    <h3>MHC Binding</h3>
    <p>Strong binders: {len(strong)}, Weak binders: {len(weak) if isinstance(weak, list) else 'N/A'}</p>"""
        parts.append(mhc_table)

    # Epitope section
    if epi and isinstance(epi, dict):
        epitopes = epi.get("epitopes", [])
        if isinstance(epitopes, list) and epitopes:
            epi_rows = ""
            for ep in epitopes[:20]:
                if isinstance(ep, dict):
                    epi_rows += f"""<tr>
                        <td>{ep.get('start', 'N/A')}</td>
                        <td>{ep.get('end', 'N/A')}</td>
                        <td class="mono">{html.escape(str(ep.get('sequence', '')))}</td>
                        <td>{ep.get('score', 0):.2f}</td>
                    </tr>"""
                else:
                    epi_rows += f"""<tr><td colspan="4">{html.escape(str(ep))}</td></tr>"""
            parts.append(f"""
    <h3>Predicted Epitopes</h3>
    <table>
        <thead><tr><th>Start</th><th>End</th><th>Sequence</th><th>Score</th></tr></thead>
        <tbody>{epi_rows}</tbody>
    </table>""")

    content = "\n".join(parts)
    return f"""
<section class="section">
    <h2>Immunogenicity</h2>
    {content}
</section>"""


def _html_predicate_section(report: ProteinAssessmentReport) -> str:
    """Build the HTML section for predicate results."""
    if not report.predicate_results:
        return ""

    rows = ""
    for i, pred in enumerate(report.predicate_results):
        name = pred.get("predicate", f"Predicate_{i}")
        verdict = pred.get("verdict", "UNCERTAIN")
        if isinstance(verdict, Verdict):
            verdict = verdict.value
        violation = pred.get("violation", "")
        v_css = verdict.lower()

        verdict_symbols = {
            "PASS": "&#10003;",
            "LIKELY_PASS": "&#10003;~",
            "UNCERTAIN": "?",
            "LIKELY_FAIL": "&#10007;~",
            "FAIL": "&#10007;",
        }
        symbol = verdict_symbols.get(verdict, "?")

        violation_html = f'<span style="color:var(--fail);font-size:0.85rem">{html.escape(violation)}</span>' if violation else ""

        rows += f"""
        <tr class="{v_css}">
            <td><span class="badge {v_css}">{symbol}</span></td>
            <td>{html.escape(name)}</td>
            <td>{violation_html}</td>
        </tr>"""

    return f"""
<section class="section">
    <h2>Predicate Results</h2>
    <table>
        <thead>
            <tr><th>Verdict</th><th>Predicate</th><th>Details</th></tr>
        </thead>
        <tbody>{rows}</tbody>
    </table>
</section>"""


def _html_recommendations_section(report: ProteinAssessmentReport) -> str:
    """Build the HTML section for recommendations."""
    if not report.recommendations:
        return ""

    items = ""
    for rec in report.recommendations:
        # Determine severity class
        cls = ""
        rec_lower = rec.lower()
        if any(kw in rec_lower for kw in ("fail", "unstable", "low confidence", "not available")):
            cls = "error"
        elif any(kw in rec_lower for kw in ("borderline", "concerning", "uncertain")):
            cls = "warning"
        items += f"<li class='{cls}'>{html.escape(rec)}</li>\n"

    return f"""
<section class="section">
    <h2>Recommendations</h2>
    <ul class="rec-list">
        {items}
    </ul>
</section>"""


def _cos_deg(angle: float) -> float:
    """Cosine of angle in degrees (pure Python, no math import needed at module level)."""
    import math
    return math.cos(math.radians(angle))


def _sin_deg(angle: float) -> float:
    """Sine of angle in degrees (pure Python, no math import needed at module level)."""
    import math
    return math.sin(math.radians(angle))


# ────────────────────────────────────────────────────────────
# SVG Visualization Functions
# ────────────────────────────────────────────────────────────

def plot_plddt_bar_svg(
    plddt_scores: list[float],
    protein: str = "",
) -> str:
    """Generate an SVG bar chart of per-residue pLDDT scores.

    Bars are colored by confidence level:
    - Green: pLDDT > 90 (very high confidence)
    - Light green: pLDDT 70-90 (confident)
    - Yellow: pLDDT 50-70 (low confidence)
    - Red: pLDDT < 50 (very low confidence)

    Includes threshold lines at 70 and 90.

    Args:
        plddt_scores: List of per-residue pLDDT values.
        protein: Optional protein name for title.

    Returns:
        SVG string.
    """
    if not plddt_scores:
        return "<p>No pLDDT data available</p>"

    width = 900
    height = 250
    margin_left = 50
    margin_right = 20
    margin_top = 35
    margin_bottom = 40
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    n = len(plddt_scores)
    bar_width = max(1, plot_w / n)
    # Cap bar width for very long proteins
    if bar_width > 8:
        bar_width = 8

    # Build bars
    bars = ""
    for i, score in enumerate(plddt_scores):
        x = margin_left + i * bar_width
        bar_h = (score / 100) * plot_h
        y = margin_top + plot_h - bar_h

        if score > 90:
            color = "#16a34a"       # green
        elif score > 70:
            color = "#86efac"       # light green
        elif score > 50:
            color = "#fbbf24"       # yellow
        else:
            color = "#dc2626"       # red

        if bar_width > 3:
            bars += f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width - 0.5:.1f}" height="{bar_h:.1f}" fill="{color}" rx="1"/>'
        else:
            bars += f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_h:.1f}" fill="{color}"/>'

    # Threshold lines
    y90 = margin_top + plot_h - (90 / 100) * plot_h
    y70 = margin_top + plot_h - (70 / 100) * plot_h
    y50 = margin_top + plot_h - (50 / 100) * plot_h

    # Axis labels
    title = f"pLDDT Confidence Scores — {html.escape(protein)}" if protein else "pLDDT Confidence Scores"

    # X-axis tick labels
    x_ticks = ""
    step = max(1, n // 10)
    for i in range(0, n, step):
        x = margin_left + i * bar_width
        x_ticks += f'<text x="{x:.0f}" y="{height - 5}" font-size="9" fill="#64748b" text-anchor="middle">{i + 1}</text>'

    return f"""<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
    <rect width="100%" height="100%" fill="white"/>
    <!-- Threshold lines -->
    <line x1="{margin_left}" y1="{y90:.1f}" x2="{width - margin_right}" y2="{y90:.1f}" stroke="#16a34a" stroke-width="1" stroke-dasharray="6,3" opacity="0.7"/>
    <line x1="{margin_left}" y1="{y70:.1f}" x2="{width - margin_right}" y2="{y70:.1f}" stroke="#fbbf24" stroke-width="1" stroke-dasharray="6,3" opacity="0.7"/>
    <line x1="{margin_left}" y1="{y50:.1f}" x2="{width - margin_right}" y2="{y50:.1f}" stroke="#dc2626" stroke-width="1" stroke-dasharray="6,3" opacity="0.7"/>
    <!-- Threshold labels -->
    <text x="{margin_left - 5}" y="{y90 + 3:.0f}" font-size="9" fill="#16a34a" text-anchor="end">90</text>
    <text x="{margin_left - 5}" y="{y70 + 3:.0f}" font-size="9" fill="#fbbf24" text-anchor="end">70</text>
    <text x="{margin_left - 5}" y="{y50 + 3:.0f}" font-size="9" fill="#dc2626" text-anchor="end">50</text>
    <!-- Bars -->
    {bars}
    <!-- Axes -->
    <line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_h}" stroke="#64748b" stroke-width="1"/>
    <line x1="{margin_left}" y1="{margin_top + plot_h}" x2="{width - margin_right}" y2="{margin_top + plot_h}" stroke="#64748b" stroke-width="1"/>
    <!-- Y-axis labels -->
    <text x="{margin_left - 5}" y="{margin_top + 4}" font-size="9" fill="#64748b" text-anchor="end">100</text>
    <text x="{margin_left - 5}" y="{margin_top + plot_h}" font-size="9" fill="#64748b" text-anchor="end">0</text>
    <!-- X-axis labels -->
    {x_ticks}
    <!-- Title -->
    <text x="{width // 2}" y="15" font-size="12" fill="#1e293b" text-anchor="middle" font-weight="600">{title}</text>
    <!-- Legend -->
    <rect x="{margin_left}" y="20" width="10" height="10" fill="#16a34a" rx="1"/>
    <text x="{margin_left + 14}" y="29" font-size="8" fill="#64748b">&gt;90</text>
    <rect x="{margin_left + 50}" y="20" width="10" height="10" fill="#86efac" rx="1"/>
    <text x="{margin_left + 64}" y="29" font-size="8" fill="#64748b">70-90</text>
    <rect x="{margin_left + 110}" y="20" width="10" height="10" fill="#fbbf24" rx="1"/>
    <text x="{margin_left + 124}" y="29" font-size="8" fill="#64748b">50-70</text>
    <rect x="{margin_left + 170}" y="20" width="10" height="10" fill="#dc2626" rx="1"/>
    <text x="{margin_left + 184}" y="29" font-size="8" fill="#64748b">&lt;50</text>
</svg>"""


def plot_solubility_profile_svg(
    per_residue_scores: list[float],
    protein: str = "",
) -> str:
    """Generate an SVG line plot of per-residue CamSol solubility scores.

    Regions above zero (soluble) are colored blue, regions below zero
    (aggregation-prone) are colored red. A zero line is drawn.

    Args:
        per_residue_scores: List of per-residue CamSol scores.
        protein: Optional protein name for title.

    Returns:
        SVG string.
    """
    if not per_residue_scores:
        return "<p>No solubility profile data available</p>"

    import math

    width = 900
    height = 250
    margin_left = 50
    margin_right = 20
    margin_top = 35
    margin_bottom = 40
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    n = len(per_residue_scores)

    # Determine Y range
    min_score = min(per_residue_scores)
    max_score = max(per_residue_scores)
    y_min = min(-1.5, min_score - 0.5)
    y_max = max(1.5, max_score + 0.5)
    y_range = y_max - y_min
    if y_range == 0:
        y_range = 1.0

    # Map functions
    def map_x(i: int) -> float:
        return margin_left + (i / max(n - 1, 1)) * plot_w

    def map_y(score: float) -> float:
        return margin_top + plot_h - ((score - y_min) / y_range) * plot_h

    # Build colored area fills
    # Blue area (above zero)
    blue_points = []
    for i, s in enumerate(per_residue_scores):
        x = map_x(i)
        y = map_y(max(s, 0))
        blue_points.append(f"{x:.1f},{y:.1f}")
    # Close the path at the zero line
    zero_y = map_y(0)
    if blue_points:
        blue_area = (f"M {blue_points[0].split(',')[0]},{zero_y:.1f} "
                     + "L ".join([f"{p} " for p in blue_points])
                     + f"L {blue_points[-1].split(',')[0]},{zero_y:.1f} Z")
    else:
        blue_area = ""

    # Red area (below zero)
    red_points = []
    for i, s in enumerate(per_residue_scores):
        x = map_x(i)
        y = map_y(min(s, 0))
        red_points.append(f"{x:.1f},{y:.1f}")
    if red_points:
        red_area = (f"M {red_points[0].split(',')[0]},{zero_y:.1f} "
                    + "L ".join([f"{p} " for p in red_points])
                    + f"L {red_points[-1].split(',')[0]},{zero_y:.1f} Z")
    else:
        red_area = ""

    # Line path
    line_points = []
    for i, s in enumerate(per_residue_scores):
        x = map_x(i)
        y = map_y(s)
        line_points.append(f"{x:.1f},{y:.1f}")
    line_path = "M " + " L ".join(line_points)

    # Zero line
    zero_y = map_y(0)

    # X-axis tick labels
    x_ticks = ""
    step = max(1, n // 10)
    for i in range(0, n, step):
        x = map_x(i)
        x_ticks += f'<text x="{x:.0f}" y="{height - 5}" font-size="9" fill="#64748b" text-anchor="middle">{i + 1}</text>'

    title = f"CamSol Solubility Profile — {html.escape(protein)}" if protein else "CamSol Solubility Profile"

    return f"""<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
    <rect width="100%" height="100%" fill="white"/>
    <!-- Blue area (soluble, score > 0) -->
    <path d="{blue_area}" fill="rgba(37,99,235,0.15)" stroke="none"/>
    <!-- Red area (aggregation-prone, score < 0) -->
    <path d="{red_area}" fill="rgba(220,38,38,0.15)" stroke="none"/>
    <!-- Zero line -->
    <line x1="{margin_left}" y1="{zero_y:.1f}" x2="{width - margin_right}" y2="{zero_y:.1f}" stroke="#64748b" stroke-width="1.5" stroke-dasharray="8,4"/>
    <text x="{margin_left - 5}" y="{zero_y + 3:.0f}" font-size="9" fill="#64748b" text-anchor="end">0.0</text>
    <!-- Line -->
    <path d="{line_path}" fill="none" stroke="#2563eb" stroke-width="1.5"/>
    <!-- Axes -->
    <line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_h}" stroke="#64748b" stroke-width="1"/>
    <line x1="{margin_left}" y1="{margin_top + plot_h}" x2="{width - margin_right}" y2="{margin_top + plot_h}" stroke="#64748b" stroke-width="1"/>
    <!-- Y-axis labels -->
    <text x="{margin_left - 5}" y="{margin_top + 4}" font-size="9" fill="#64748b" text-anchor="end">{y_max:.1f}</text>
    <text x="{margin_left - 5}" y="{margin_top + plot_h}" font-size="9" fill="#64748b" text-anchor="end">{y_min:.1f}</text>
    <!-- X-axis labels -->
    {x_ticks}
    <!-- Title -->
    <text x="{width // 2}" y="15" font-size="12" fill="#1e293b" text-anchor="middle" font-weight="600">{title}</text>
    <!-- Legend -->
    <rect x="{margin_left}" y="20" width="10" height="10" fill="rgba(37,99,235,0.3)" rx="1"/>
    <text x="{margin_left + 14}" y="29" font-size="8" fill="#64748b">Soluble (&gt;0)</text>
    <rect x="{margin_left + 90}" y="20" width="10" height="10" fill="rgba(220,38,38,0.3)" rx="1"/>
    <text x="{margin_left + 104}" y="29" font-size="8" fill="#64748b">Aggregation-prone (&lt;0)</text>
</svg>"""
