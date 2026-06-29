"""
BioCompiler Application — Protein assessment service.

Orchestrates structure prediction, quality assessment, stability analysis,
solubility analysis, immunogenicity analysis, deimmunization, and full
assessment.
"""

import logging
import time
from typing import Any, Optional

# ESMFOLD_TIMEOUT_S is imported lazily where needed to avoid circular imports.
# Protein analysis input models are passed as simple typed parameters.

logger = logging.getLogger(__name__)


# ─── Assessment Verdict → Verdict Enum Mapping ────────────────────

def assessment_verdict_to_verdict(assessment_verdict: str):
    """Map assessment service verdict to type-system Verdict.

    Mapping:
    - "STABLE" → Verdict.PASS
    - "MARGINAL" → Verdict.UNCERTAIN
    - "UNSTABLE" → Verdict.FAIL

    Args:
        assessment_verdict: String verdict from the assessment service.

    Returns:
        Corresponding Verdict enum value.
    """
    from biocompiler.shared.types import Verdict

    mapping = {
        "STABLE": Verdict.PASS,
        "MARGINAL": Verdict.UNCERTAIN,
        "UNSTABLE": Verdict.FAIL,
    }
    return mapping.get(assessment_verdict, Verdict.UNCERTAIN)


# ─── Structure ─────────────────────────────────────────────────────

def predict_structure(protein: str, organism: str = "Homo_sapiens") -> dict[str, Any]:
    """Predict protein 3D structure using ESMFold.

    Returns a dict with keys: pdb_string, mean_plddt, plddt_scores,
    quality_class, execution_time_s, success, error.
    """
    from biocompiler.engines.esmfold import predict_structure as _predict_structure, is_esmfold_available

    if not is_esmfold_available():
        raise RuntimeError("ESMFold model is not available. Ensure torch and esm are installed.")

    t0 = time.monotonic()
    result = _predict_structure(protein=protein)
    elapsed = time.monotonic() - t0

    mean_plddt = float(getattr(result, "mean_plddt", 0.0))
    if mean_plddt >= 90:
        quality_class = "very_high"
    elif mean_plddt >= 70:
        quality_class = "high"
    elif mean_plddt >= 50:
        quality_class = "medium"
    else:
        quality_class = "low"

    return {
        "pdb_string": result.pdb_string,
        "mean_plddt": mean_plddt,
        "plddt_scores": getattr(result, "plddt_scores", []),
        "quality_class": quality_class,
        "execution_time_s": round(elapsed, 3),
        "success": getattr(result, "success", True),
        "error": getattr(result, "error", None),
    }


def assess_structure_quality(pdb_string: str) -> dict[str, Any]:
    """Assess structure quality from a PDB string.

    Returns a dict with keys: mean_plddt, ramachandran_favored,
    clash_score, overall_quality, verdict.
    """
    from ..structure.quality import compute_structure_quality

    result = compute_structure_quality(pdb_string=pdb_string)
    return result


# ─── Stability ─────────────────────────────────────────────────────

def analyze_stability(protein: str) -> dict[str, Any]:
    """Analyze protein stability using FoldX.

    Returns a dict with keys: stability_kcal, method, components, verdict.
    """
    from biocompiler.engines.foldx import empirical_stability

    result = empirical_stability(protein=protein)

    components = {
        key: getattr(result, key)
        for key in (
            "backbone_hbond", "sidechain_hbond", "van_der_waals",
            "electrostatics", "solvation", "van_der_waals_clashes",
            "entropy_sidechain", "entropy_mainchain", "torsional_clash",
            "backbone_clash", "helix_dipole", "disulfide",
            "electrostatic_kon", "partial_covalent", "energy_ionisation",
        )
        if getattr(result, key, None) is not None
    }
    stability_kcal = result.stability_kcal
    if stability_kcal < -5.0:
        verdict = "STABLE"
    elif stability_kcal < 0.0:
        verdict = "MARGINAL"
    else:
        verdict = "UNSTABLE"

    return {
        "stability_kcal": stability_kcal,
        "method": result.method,
        "components": components,
        "verdict": verdict,
        "verdict_enum": assessment_verdict_to_verdict(verdict),
    }


def scan_stability_mutations(
    protein: str,
    positions: Optional[list[int]] = None,
    method: str = "empirical",
) -> dict[str, Any]:
    """Scan single-point mutations for stability effects.

    Returns a dict with keys: mutations, stabilizing_count, destabilizing_count.
    """
    from biocompiler.engines.foldx import scan_mutations

    result = scan_mutations(protein=protein, positions=positions)

    mutations = [
        {
            "position": m.position,
            "original": m.original,
            "mutant": m.mutant,
            "score": m.score,
            "description": m.description,
            "details": m.details,
        }
        for m in result
    ]
    stabilizing_count = sum(
        1 for m in result if m.details and m.details.get("stabilizing")
    )
    destabilizing_count = sum(
        1 for m in result if m.details and m.details.get("destabilizing")
    )

    return {
        "mutations": mutations,
        "stabilizing_count": stabilizing_count,
        "destabilizing_count": destabilizing_count,
    }


# ─── Solubility ────────────────────────────────────────────────────

def analyze_solubility(protein: str, pdb_string: Optional[str] = None) -> dict[str, Any]:
    """Analyze protein solubility using CamSol algorithm.

    Returns a dict with keys: intrinsic_score, overall_score,
    solubility_class, aggregation_prone_regions.
    """
    from biocompiler.engines.camsol import compute_solubility

    result = compute_solubility(protein=protein, pdb_string=pdb_string)

    return {
        "intrinsic_score": result.intrinsic_score,
        "overall_score": result.overall_score,
        "solubility_class": result.solubility_class,
        "aggregation_prone_regions": result.aggregation_prone_regions,
    }


def find_solubility_mutations(protein: str) -> dict[str, Any]:
    """Find mutations that improve protein solubility.

    Returns a dict with keys: mutations, stabilizing_count, destabilizing_count.
    """
    from biocompiler.engines.camsol import find_solubility_mutations

    result = find_solubility_mutations(protein=protein)

    mutations = [
        {
            "position": m.position,
            "original": m.original,
            "mutant": m.mutant,
            "score": m.score,
            "description": m.description,
            "details": m.details,
        }
        for m in result
    ]
    stabilizing = sum(1 for m in result if m.score > 0)
    destabilizing = sum(1 for m in result if m.score <= 0)

    return {
        "mutations": mutations,
        "stabilizing_count": stabilizing,
        "destabilizing_count": destabilizing,
    }


# ─── Immunogenicity ────────────────────────────────────────────────

def analyze_immunogenicity(
    protein: str,
    mhc_alleles: Optional[list[str]] = None,
) -> dict[str, Any]:
    """Analyze protein immunogenicity.

    Returns a dict with keys: overall_score, immunogenicity_class,
    t_cell_epitopes, b_cell_epitopes, deimmunization_candidates.
    """
    from biocompiler.immunogenicity.core import compute_immunogenicity

    result = compute_immunogenicity(protein=protein, mhc_alleles=mhc_alleles)

    deimmunization_candidates = [
        {
            "position": m.position,
            "original": m.original,
            "mutant": m.mutant,
            "score": m.score,
            "description": m.description,
        }
        for m in result.deimmunization_candidates
    ]

    return {
        "overall_score": result.overall_score,
        "immunogenicity_class": result.immunogenicity_class,
        "t_cell_epitopes": result.t_cell_epitopes,
        "b_cell_epitopes": result.b_cell_epitopes,
        "deimmunization_candidates": deimmunization_candidates,
    }


def deimmunize_protein(
    protein: str,
    organism: str = "Homo_sapiens",
    target_score: float = 0.3,
    max_mutations: int = 10,
    blosum62_min: int = 0,
) -> dict[str, Any]:
    """Deimmunize a protein by introducing conservative substitutions.

    Returns a dict with keys: optimized_protein, mutations_applied,
    original_score, optimized_score, success.
    """
    from biocompiler.immunogenicity.deimmunization import deimmunize

    result = deimmunize(
        protein=protein,
        organism=organism,
        target_score=target_score,
        max_mutations=max_mutations,
        blosum62_min=blosum62_min,
    )

    return {
        "optimized_protein": result.optimized_protein,
        "mutations_applied": len(result.mutations_applied),
        "original_score": result.original_immunogenicity,
        "optimized_score": result.optimized_immunogenicity,
        "success": result.success,
    }


# ─── Full Assessment ───────────────────────────────────────────────

def full_assessment(
    protein: str,
    organism: str = "Homo_sapiens",
    pdb_string: Optional[str] = None,
    run_structure: bool = True,
    run_stability: bool = True,
    run_solubility: bool = True,
    run_immunogenicity: bool = True,
) -> dict[str, Any]:
    """Run a full protein assessment combining all analyses.

    Returns a dict with keys: structure_quality, stability, solubility,
    immunogenicity, predicate_results, overall_verdict, recommendations.
    """
    structure_quality_result: dict | None = None
    stability_result: dict | None = None
    solubility_result: dict | None = None
    immunogenicity_result: dict | None = None
    predicate_results: list[dict] = []
    recommendations: list[str] = []

    # ── Structure Prediction & Quality ─────────────────────────────
    if run_structure:
        pdb_str = pdb_string
        if pdb_str is None:
            # Predict structure first
            try:
                from biocompiler.engines.esmfold import predict_structure as _predict_structure, is_esmfold_available

                if is_esmfold_available():
                    pred = _predict_structure(protein=protein)
                    pdb_str = pred.pdb_string
                    mean_plddt = float(getattr(pred, "mean_plddt", 0.0))

                    if mean_plddt < 50:
                        recommendations.append(
                            f"Predicted structure has low confidence "
                            f"(mean pLDDT={mean_plddt:.1f}). Consider experimental "
                            f"structure determination."
                        )
                else:
                    recommendations.append(
                        "ESMFold is not available. Structure prediction skipped. "
                        "Provide a PDB string or install biocompiler[esmfold]."
                    )
            except ImportError:
                recommendations.append(
                    "ESMFold module not available. Structure prediction skipped."
                )
            except Exception as e:
                logger.warning("Structure prediction failed: %s", e)
                recommendations.append(f"Structure prediction failed: {e}")

        # Assess structure quality if we have a PDB string
        if pdb_str is not None:
            try:
                from ..structure.quality import compute_structure_quality

                quality = compute_structure_quality(pdb_string=pdb_str)
                structure_quality_result = quality

                if quality["verdict"] == "FAIL":
                    recommendations.append(
                        "Structure quality is poor. Results from stability and "
                        "solubility analyses may be unreliable."
                    )
                elif quality["verdict"] == "WARN":
                    recommendations.append(
                        "Structure quality is acceptable but not ideal. "
                        "Consider verifying with experimental data."
                    )
            except ImportError:
                recommendations.append(
                    "Structure quality module not available. Quality assessment skipped."
                )
            except Exception as e:
                logger.warning("Structure quality assessment failed: %s", e)
                recommendations.append(f"Structure quality assessment failed: {e}")

    # ── Stability Analysis ─────────────────────────────────────────
    if run_stability:
        try:
            from biocompiler.engines.foldx import empirical_stability

            stab = empirical_stability(protein=protein)
            stability_result = {
                "stability_kcal": stab.stability_kcal,
                "method": stab.method,
                "success": stab.success,
            }
            if stab.stability_kcal < -5.0:
                stability_result["verdict"] = "STABLE"
            elif stab.stability_kcal < 0.0:
                stability_result["verdict"] = "MARGINAL"
            else:
                stability_result["verdict"] = "UNSTABLE"

            stability_result["verdict_enum"] = assessment_verdict_to_verdict(
                stability_result["verdict"]
            )

            if stability_result["verdict"] == "UNSTABLE":
                recommendations.append(
                    f"Protein is predicted to be unstable "
                    f"(ΔG={stab.stability_kcal:.1f} kcal/mol). "
                    f"Consider stability-enhancing mutations."
                )
            elif stability_result["verdict"] == "MARGINAL":
                recommendations.append(
                    f"Protein stability is marginal "
                    f"(ΔG={stab.stability_kcal:.1f} kcal/mol). "
                    f"Monitor during expression."
                )
        except ImportError:
            recommendations.append(
                "FoldX stability module not available. Stability analysis skipped."
            )
        except Exception as e:
            logger.warning("Stability analysis failed: %s", e)
            recommendations.append(f"Stability analysis failed: {e}")

    # ── Solubility Analysis ────────────────────────────────────────
    if run_solubility:
        try:
            from biocompiler.engines.camsol import compute_solubility

            sol = compute_solubility(protein=protein, pdb_string=pdb_string)
            solubility_result = {
                "intrinsic_score": sol.intrinsic_score,
                "overall_score": sol.overall_score,
                "solubility_class": sol.solubility_class,
                "aggregation_prone_regions": sol.aggregation_prone_regions,
            }

            if sol.solubility_class in ("low", "very_low", "marginally_soluble", "insoluble"):
                recommendations.append(
                    f"Protein solubility is {sol.solubility_class} "
                    f"(score={sol.overall_score:.2f}). "
                    f"Consider solubility-enhancing mutations or fusion tags."
                )
            if sol.aggregation_prone_regions:
                recommendations.append(
                    f"Found {len(sol.aggregation_prone_regions)} "
                    f"aggregation-prone region(s). Review for potential redesign."
                )
        except ImportError:
            recommendations.append(
                "CamSol solubility module not available. Solubility analysis skipped."
            )
        except Exception as e:
            logger.warning("Solubility analysis failed: %s", e)
            recommendations.append(f"Solubility analysis failed: {e}")

    # ── Immunogenicity Analysis ────────────────────────────────────
    if run_immunogenicity:
        try:
            from biocompiler.immunogenicity.core import compute_immunogenicity

            imm = compute_immunogenicity(protein=protein, mhc_alleles=None)
            immunogenicity_result = {
                "overall_score": imm.overall_score,
                "immunogenicity_class": imm.immunogenicity_class,
                "t_cell_epitopes": imm.t_cell_epitopes,
                "b_cell_epitopes": imm.b_cell_epitopes,
                "deimmunization_candidates": [
                    {
                        "position": m.position,
                        "original": m.original,
                        "mutant": m.mutant,
                        "score": m.score,
                    }
                    for m in imm.deimmunization_candidates
                ],
            }

            if imm.immunogenicity_class in ("high", "very_high"):
                recommendations.append(
                    f"Protein immunogenicity is {imm.immunogenicity_class} "
                    f"(score={imm.overall_score:.2f}). "
                    f"Consider deimmunization for therapeutic applications."
                )
            if imm.deimmunization_candidates:
                recommendations.append(
                    f"Found {len(imm.deimmunization_candidates)} "
                    f"deimmunization candidate position(s)."
                )
        except ImportError:
            recommendations.append(
                "Immunogenicity module not available. Immunogenicity analysis skipped."
            )
        except Exception as e:
            logger.warning("Immunogenicity analysis failed: %s", e)
            recommendations.append(f"Immunogenicity analysis failed: {e}")

    # ── Predicate Results ──────────────────────────────────────────
    try:
        from ..structure.report import assess_protein

        report = assess_protein(
            protein=protein,
            organism=organism,
            pdb_string=pdb_string,
        )
        predicate_results = report.get("predicate_results", [])
    except ImportError:
        # Fallback: try type system directly
        try:
            from biocompiler.optimizer import optimize_sequence

            opt_result = optimize_sequence(
                target_protein=protein,
                organism=organism,
            )
            predicate_results = [
                {"predicate": p, "verdict": "PASS"}
                for p in opt_result.satisfied_predicates
            ] + [
                {"predicate": p, "verdict": "FAIL"}
                for p in opt_result.failed_predicates
            ]
        except Exception as exc:
            logger.warning("Predicate fallback evaluation failed: %s", exc)
            predicate_results = []
    except Exception as e:
        logger.warning("Predicate evaluation failed: %s", e)
        predicate_results = []

    # ── Overall Verdict ────────────────────────────────────────────
    from biocompiler.shared.types import Verdict as _Verdict, combined_verdict as _combined_verdict

    enum_verdicts: list[_Verdict] = []

    if structure_quality_result is not None:
        sq_v = structure_quality_result.get("verdict", "WARN")
        # Structure quality already uses PASS/WARN/FAIL — map directly
        _sq_map = {"PASS": _Verdict.PASS, "WARN": _Verdict.UNCERTAIN, "FAIL": _Verdict.FAIL}
        enum_verdicts.append(_sq_map.get(sq_v, _Verdict.UNCERTAIN))
    if stability_result is not None:
        stab_v = stability_result.get("verdict", "MARGINAL")
        enum_verdicts.append(assessment_verdict_to_verdict(stab_v))
    if solubility_result is not None:
        sol_class = solubility_result.get("solubility_class", "medium")
        if sol_class in ("high",):
            enum_verdicts.append(_Verdict.PASS)
        elif sol_class == "medium":
            enum_verdicts.append(_Verdict.UNCERTAIN)
        else:
            enum_verdicts.append(_Verdict.FAIL)
    if immunogenicity_result is not None:
        imm_class = immunogenicity_result.get("immunogenicity_class", "moderate")
        if imm_class in ("low",):
            enum_verdicts.append(_Verdict.PASS)
        elif imm_class == "moderate":
            enum_verdicts.append(_Verdict.UNCERTAIN)
        else:
            enum_verdicts.append(_Verdict.FAIL)

    overall_verdict_enum = _combined_verdict(enum_verdicts) if enum_verdicts else _Verdict.UNCERTAIN
    overall_verdict = overall_verdict_enum.value

    if not recommendations:
        recommendations.append("No issues detected. Protein design looks good.")

    return {
        "structure_quality": structure_quality_result,
        "stability": stability_result,
        "solubility": solubility_result,
        "immunogenicity": immunogenicity_result,
        "predicate_results": predicate_results,
        "overall_verdict": overall_verdict,
        "overall_verdict_enum": overall_verdict_enum,
        "recommendations": recommendations,
    }


# ─── Batch Item Processors ────────────────────────────────────────

def structure_batch_item(protein: str) -> dict[str, Any]:
    """Process a single structure prediction item for batch processing."""
    from biocompiler.engines.esmfold import predict_structure as _predict_structure

    result = _predict_structure(protein=protein)

    if not result.success:
        return {"error": result.error or "Prediction failed"}

    mean_plddt = float(getattr(result, "mean_plddt", 0.0))
    if mean_plddt >= 90:
        quality_class = "very_high"
    elif mean_plddt >= 70:
        quality_class = "high"
    elif mean_plddt >= 50:
        quality_class = "medium"
    else:
        quality_class = "low"

    return {
        "pdb_string": result.pdb_string,
        "mean_plddt": mean_plddt,
        "quality_class": quality_class,
    }


def stability_batch_item(protein: str) -> dict[str, Any]:
    """Process a single stability analysis item for batch processing."""
    return analyze_stability(protein=protein)


def solubility_batch_item(protein: str, pdb_string: Optional[str] = None) -> dict[str, Any]:
    """Process a single solubility analysis item for batch processing."""
    return analyze_solubility(protein=protein, pdb_string=pdb_string)


def immunogenicity_batch_item(protein: str, mhc_alleles: Optional[list[str]] = None) -> dict[str, Any]:
    """Process a single immunogenicity analysis item for batch processing."""
    result = analyze_immunogenicity(
        protein=protein,
        mhc_alleles=mhc_alleles,
    )
    # Batch responses do not include deimmunization_candidates
    return {
        "overall_score": result["overall_score"],
        "immunogenicity_class": result["immunogenicity_class"],
        "t_cell_epitopes": result["t_cell_epitopes"],
        "b_cell_epitopes": result["b_cell_epitopes"],
    }
