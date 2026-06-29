"""Conformational and composite B-cell epitope prediction.

PDB-based conformational epitope prediction, composite linear+conformational
epitope aggregation, approximate surface accessibility, the deprecated
``predict_b_cell_epitopes`` wrapper, and the real-binding-data availability
check.

Split out of ``core.py`` (W8-a refactor).
"""
from __future__ import annotations

import concurrent.futures
import hashlib
import logging
import math
import warnings
from dataclasses import dataclass, field
from typing import List, Optional, TypedDict

from biocompiler.shared.constants import (
    BLOSUM62,
    DEFAULT_MHC_PEPTIDE_LENGTH,
    HYDROPATHY,
    STANDARD_AAS,
)
from ..engines.base import (
    BaseEngineResult,
    BatchResult,
    EngineTimer,
    MutationResult,
    validate_protein_sequence,
)
from biocompiler.shared.exceptions import ImmunogenicityError
from biocompiler.shared.types import Verdict

logger = logging.getLogger(__name__)
from ._constants import *  # noqa: F401,F403
from ._models import *  # noqa: F401,F403
from .epitopes import *  # noqa: F401,F403
from .epitopes import _validate_protein  # noqa: F401

def predict_conformational_epitopes(
    pdb_string: str,
    distance_cutoff: float = 6.0,
) -> list[EpitopeRegion]:
    """Predict conformational B-cell epitopes from a PDB structure.

    Identifies surface patches on the protein structure and scores them
    by hydrophilicity, charge, and flexibility. Surface residues are
    identified by having fewer than 15 C-alpha neighbors within 12 A.
    Adjacent surface residues (within distance_cutoff in sequence) are
    clustered into patches.

    Args:
        pdb_string: PDB file content as a string.
        distance_cutoff: Maximum sequence gap to cluster adjacent surface
                        residues into patches (default 6.0).

    Returns:
        List of EpitopeRegion predictions with is_linear=False.
    """
    if not pdb_string:
        return []

    ca_atoms: list[tuple[int, float, float, float, str]] = []

    for line in pdb_string.splitlines():
        line = line.rstrip()
        if not line.startswith("ATOM"):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue

        try:
            resnum = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except (ValueError, IndexError) as exc:
            logger.debug(
                "Skipping malformed PDB ATOM line (resnum/coords): %s", exc,
            )
            continue

        resname = line[17:20].strip() if len(line) >= 20 else ""
        aa = _THREE_TO_ONE.get(resname, "X")
        ca_atoms.append((resnum, x, y, z, aa))

    if len(ca_atoms) < 3:
        logger.warning(
            "Too few C-alpha atoms (%d) in PDB for conformational epitope prediction",
            len(ca_atoms),
        )
        return []

    surface_indices: set[int] = set()
    for i, (_, xi, yi, zi, _) in enumerate(ca_atoms):
        neighbors = 0
        for j, (_, xj, yj, zj, _) in enumerate(ca_atoms):
            if i == j:
                continue
            dist_sq = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if dist_sq <= CONF_EPITOPE_NEIGHBOR_CUTOPT_ANGSTROM ** 2:
                neighbors += 1
                if neighbors >= CONF_EPITOPE_MAX_NEIGHBORS:
                    break
        if neighbors < CONF_EPITOPE_MAX_NEIGHBORS:
            surface_indices.add(i)

    if not surface_indices:
        return []

    sorted_surface = sorted(surface_indices)
    patches: list[list[int]] = []
    current_patch = [sorted_surface[0]]

    for k in range(1, len(sorted_surface)):
        if sorted_surface[k] - sorted_surface[k - 1] <= distance_cutoff:
            current_patch.append(sorted_surface[k])
        else:
            if len(current_patch) >= 3:
                patches.append(current_patch)
            current_patch = [sorted_surface[k]]

    if len(current_patch) >= 3:
        patches.append(current_patch)

    charged_aas = {"R", "K", "D", "E"}
    flexible_aas = {"G", "P"}

    epitopes: list[EpitopeRegion] = []

    for patch in patches:
        aas: list[str] = []
        hydro_sum = 0.0
        charge_count = 0
        flex_count = 0

        for idx in patch:
            _, _, _, _, aa = ca_atoms[idx]
            aas.append(aa)
            hydro_sum += PARKER_SCALE.get(aa, 0.0)
            if aa in charged_aas:
                charge_count += 1
            if aa in flexible_aas:
                flex_count += 1

        n_res = len(patch)
        hydro_avg = hydro_sum / n_res
        charge_frac = charge_count / n_res
        flex_frac = flex_count / n_res

        hydro_norm = max(0.0, min(1.0, (hydro_avg + 1.5) / 3.3))

        score = 0.4 * hydro_norm + 0.3 * charge_frac + 0.3 * flex_frac

        first_resnum = ca_atoms[patch[0]][0]
        last_resnum = ca_atoms[patch[-1]][0]
        start = first_resnum - 1
        end = last_resnum

        epitopes.append(EpitopeRegion(
            start=start,
            end=end,
            peptide="".join(aas),
            score=round(score, 4),
            method="conformational",
            is_linear=False,
            properties={
                "hydrophilicity_avg": round(hydro_avg, 4),
                "charge_fraction": round(charge_frac, 4),
                "flexibility_fraction": round(flex_frac, 4),
                "surface_residue_count": n_res,
                "pdb_residue_range": (first_resnum, last_resnum),
            },
        ))

    return epitopes


# Method dispatch
_METHOD_MAP: dict[str, object] = {
    "kolaskar_tongaonkar": predict_kolaskar_tongaonkar,
    "parker_hydrophilicity": predict_parker_hydrophilicity,
    "chou_fasman_beta_turn": predict_chou_fasman_beta_turn,
    "eea": predict_eea,
    "bepipred": predict_bepipred_like,
}


# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: combined prediction
# ═══════════════════════════════════════════════════════════════════════════


def predict_epitopes(
    protein: str,
    pdb_string: str | None = None,
    methods: list[str] | None = None,
) -> EpitopePredictionResult:
    """Run multiple B-cell epitope prediction methods and combine results.

    Default methods: kolaskar_tongaonkar, parker_hydrophilicity,
    chou_fasman_beta_turn. If a PDB structure is provided,
    conformational epitope prediction is also run.

    Per-residue scores are computed as the average across all methods
    that predict each residue as an epitope. Consensus epitope regions
    are those predicted by >= 2 methods.

    Args:
        protein: Amino acid sequence (1-letter codes).
        pdb_string: Optional PDB file content for conformational prediction.
        methods: List of method names to use. Defaults to
                 ["kolaskar_tongaonkar", "parker_hydrophilicity",
                  "chou_fasman_beta_turn"].

    Returns:
        EpitopePredictionResult with combined predictions.
    """
    if not protein:
        return EpitopePredictionResult(
            protein="",
            linear_epitopes=[],
            conformational_epitopes=[],
            per_residue_score=[],
            epitope_coverage=0.0,
            methods_used=[],
        )

    protein = protein.upper()
    n = len(protein)

    if methods is None:
        methods = [
            "kolaskar_tongaonkar",
            "parker_hydrophilicity",
            "chou_fasman_beta_turn",
        ]

    valid_methods = [m for m in methods if m in _METHOD_MAP]
    if not valid_methods:
        logger.warning("No valid methods specified for epitope prediction")
        return EpitopePredictionResult(
            protein=protein,
            linear_epitopes=[],
            conformational_epitopes=[],
            per_residue_score=[0.0] * n,
            epitope_coverage=0.0,
            methods_used=[],
        )

    all_linear_epitopes: list[EpitopeRegion] = []
    residue_method_count: list[int] = [0] * n
    residue_score_sum: list[float] = [0.0] * n

    for method_name in valid_methods:
        func = _METHOD_MAP[method_name]
        try:
            epitopes = func(protein)  # type: ignore[operator]
            all_linear_epitopes.extend(epitopes)

            for ep in epitopes:
                for pos in range(ep.start, ep.end):
                    if 0 <= pos < n:
                        residue_method_count[pos] += 1
                        residue_score_sum[pos] += ep.score
        except Exception as e:
            logger.warning("Method %s failed: %s", method_name, e)

    per_residue_score: list[float] = []
    for i in range(n):
        if residue_method_count[i] > 0:
            per_residue_score.append(
                residue_score_sum[i] / residue_method_count[i]
            )
        else:
            per_residue_score.append(0.0)

    consensus_epitopes: list[EpitopeRegion] = []
    in_region = False
    region_start = 0
    min_consensus = 2
    min_length = 6

    for i in range(n):
        if residue_method_count[i] >= min_consensus and not in_region:
            in_region = True
            region_start = i
        elif residue_method_count[i] < min_consensus and in_region:
            in_region = False
            if i - region_start >= min_length:
                length = i - region_start
                avg_score = sum(per_residue_score[region_start:i]) / length
                method_count_avg = sum(
                    residue_method_count[region_start:i]
                ) / length
                consensus_epitopes.append(EpitopeRegion(
                    start=region_start,
                    end=i,
                    peptide=protein[region_start:i],
                    score=round(avg_score, 4),
                    method="consensus",
                    is_linear=True,
                    properties={
                        "method_count_avg": round(method_count_avg, 4),
                        "contributing_methods": min_consensus,
                    },
                ))

    if in_region:
        i = n
        if i - region_start >= min_length:
            length = i - region_start
            avg_score = sum(per_residue_score[region_start:i]) / length
            method_count_avg = sum(
                residue_method_count[region_start:i]
            ) / length
            consensus_epitopes.append(EpitopeRegion(
                start=region_start,
                end=i,
                peptide=protein[region_start:i],
                score=round(avg_score, 4),
                method="consensus",
                is_linear=True,
                properties={
                    "method_count_avg": round(method_count_avg, 4),
                    "contributing_methods": min_consensus,
                },
            ))

    final_linear = list(all_linear_epitopes) + consensus_epitopes

    conformational_epitopes: list[EpitopeRegion] = []
    if pdb_string:
        try:
            conformational_epitopes = predict_conformational_epitopes(pdb_string)
        except Exception as e:
            logger.warning("Conformational epitope prediction failed: %s", e)

    epitope_residues: set[int] = set()
    for ep in all_linear_epitopes:
        for pos in range(ep.start, ep.end):
            epitope_residues.add(pos)
    for ep in consensus_epitopes:
        for pos in range(ep.start, ep.end):
            epitope_residues.add(pos)
    for ep in conformational_epitopes:
        for pos in range(ep.start, ep.end):
            epitope_residues.add(pos)

    coverage = len(epitope_residues) / n if n > 0 else 0.0

    methods_used = list(valid_methods)
    if pdb_string:
        methods_used.append("conformational")
    if consensus_epitopes:
        methods_used.append("consensus")

    return EpitopePredictionResult(
        protein=protein,
        linear_epitopes=final_linear,
        conformational_epitopes=conformational_epitopes,
        per_residue_score=per_residue_score,
        epitope_coverage=round(coverage, 4),
        methods_used=methods_used,
    )


# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: backward compatibility
# ═══════════════════════════════════════════════════════════════════════════


def compute_surface_accessibility_approx(protein: str) -> list[float]:
    """Approximate relative surface accessibility per residue.

    Based on amino-acid type and local flexibility. The result is a
    per-residue value in the range [0, 1].

    .. deprecated::
        Prefer :func:`predict_eea` for Emini surface accessibility or
        :func:`predict_epitopes` for combined B-cell predictions.

    Parameters
    ----------
    protein : str
        Amino-acid sequence.

    Returns
    -------
    list[float]
        Per-residue surface accessibility estimate.
    """
    import warnings as _warnings
    _warnings.warn(
        "compute_surface_accessibility_approx() is deprecated — use "
        "predict_eea() for Emini surface accessibility or "
        "predict_epitopes() for combined B-cell predictions instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    protein = _validate_protein(protein)
    n = len(protein)
    if n == 0:
        return []

    _surface_base: dict[str, float] = {
        "A": 0.45, "C": 0.30, "D": 0.75, "E": 0.78, "F": 0.35,
        "G": 0.55, "H": 0.65, "I": 0.30, "K": 0.80, "L": 0.30,
        "M": 0.40, "N": 0.70, "P": 0.70, "Q": 0.72, "R": 0.78,
        "S": 0.65, "T": 0.60, "V": 0.30, "W": 0.40, "Y": 0.55,
    }

    _flexibility: dict[str, float] = {
        "A": 0.35, "C": 0.25, "D": 0.45, "E": 0.45, "F": 0.20,
        "G": 0.60, "H": 0.35, "I": 0.20, "K": 0.40, "L": 0.20,
        "M": 0.30, "N": 0.45, "P": 0.55, "Q": 0.40, "R": 0.40,
        "S": 0.45, "T": 0.40, "V": 0.20, "W": 0.20, "Y": 0.30,
    }

    accessibility: list[float] = []
    for i, aa in enumerate(protein):
        base = _surface_base.get(aa, 0.40)
        flex = _flexibility.get(aa, 0.30)

        win_start = max(0, i - 2)
        win_end = min(n, i + 3)
        local_flex = sum(
            _flexibility.get(protein[j], 0.30) for j in range(win_start, win_end)
        ) / (win_end - win_start)

        terminal_boost = 0.0
        if i < 3 or i >= n - 3:
            terminal_boost = 0.15 * (1.0 - min(i, n - 1 - i) / 3.0)

        combined = 0.60 * base + 0.30 * local_flex + 0.10 * flex + terminal_boost
        accessibility.append(max(0.0, min(1.0, combined)))

    return accessibility


def predict_b_cell_epitopes(
    protein: str,
    method: str = "kolaskar_tongaonkar",
) -> list[BCellEpitopeDict]:
    """Predict B-cell epitopes.

    .. deprecated::
        Prefer :func:`predict_kolaskar_tongaonkar` or
        :func:`predict_epitopes` for richer results.

    Parameters
    ----------
    protein : str
        Amino-acid sequence.
    method : str
        Prediction method (currently only "kolaskar_tongaonkar").

    Returns
    -------
    list[BCellEpitopeDict]
        Each dict: start, end, peptide, score, antigenic.
    """
    import warnings as _warnings
    _warnings.warn(
        "predict_b_cell_epitopes() is deprecated — use "
        "predict_kolaskar_tongaonkar() or predict_epitopes() "
        "for richer results instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    protein = _validate_protein(protein)

    if method != "kolaskar_tongaonkar":
        raise ImmunogenicityError(
            f"Unsupported B-cell epitope method: {method!r}. "
            "Only 'kolaskar_tongaonkar' is supported."
        )

    regions = predict_kolaskar_tongaonkar(protein)

    # Convert EpitopeRegion objects to dicts for backward compatibility
    result: list[BCellEpitopeDict] = []
    for r in regions:
        avg_prop = r.score
        result.append(BCellEpitopeDict(
            start=r.start,
            end=r.end,
            peptide=r.peptide,
            score=round(avg_prop, 4),
            antigenic=avg_prop >= 0.5,
        ))

    if not result:
        logger.info(
            "predict_b_cell_epitopes: no antigenic regions found above threshold "
            "for protein of length %d — returning empty list",
            len(protein),
        )

    return result


# ═══════════════════════════════════════════════════════════════════════════
# Honesty: real binding-data availability check (TIGHTEN-4)
# ═══════════════════════════════════════════════════════════════════════════

# Cached result for _check_real_binding_data_available
_real_binding_data_cache: bool | None = None


def _check_real_binding_data_available() -> bool:
    """Return True iff a real, data-trained MHC binding predictor is available.

    The default PSSMs in this module are hand-crafted approximations (see
    ``_build_mhc_i_pssms`` / ``_build_mhc_ii_pssms`` — the comments literally
    say "guessed/approximate scores, NOT scores derived from experimental
    binding data").  They MUST NOT be used to issue a PASS/FAIL immunogenicity
    verdict.

    A "real" predictor is one whose scores are derived from experimentally
    measured binding data.  We recognise two such predictors:

    1. **NetMHCpan** (DTU, trained on >800,000 IEDB measurements) — installed
       locally as a binary on PATH.
    2. **MHCflurry** (O\'Donnell et al. 2018, trained on IEDB data) — Python
       package with downloaded models.

    Returns
    -------
    bool
        True if at least one real predictor is available.  False if only the
        guessed PSSMs are available (in which case :func:`compute_immunogenicity`
        MUST mark its result as UNCERTAIN with ``reason="fabricated_scores"``).

    Notes
    -----
    The result is cached after the first call to avoid repeatedly probing
    PATH / model directories on every immunogenicity computation.
    """
    global _real_binding_data_cache
    if _real_binding_data_cache is not None:
        return _real_binding_data_cache

    available = False

    # 1. NetMHCpan binary on PATH
    if not available:
        try:
            from biocompiler.immunogenicity.netmhcpan import is_netmhcpan_installed
            if is_netmhcpan_installed():
                available = True
        except Exception as exc:  # pragma: no cover - import guard
            logger.debug("NetMHCpan availability check failed: %s", exc)

    # 2. MHCflurry with downloaded models
    if not available:
        try:
            from biocompiler.immunogenicity.mhcflurry_adapter import is_mhcflurry_available
            if is_mhcflurry_available():
                available = True
        except Exception as exc:  # pragma: no cover - import guard
            logger.debug("MHCflurry availability check failed: %s", exc)

    _real_binding_data_cache = available
    if not available:
        logger.warning(
            "compute_immunogenicity: no real MHC binding-data predictor "
            "available (NetMHCpan / MHCflurry missing). Falling back to "
            "guessed PSSMs — verdict will be UNCERTAIN."
        )
    return available


# ═══════════════════════════════════════════════════════════════════════════
# Combined immunogenicity scoring
# ═══════════════════════════════════════════════════════════════════════════



__all__ = [
    "predict_conformational_epitopes",
    "predict_epitopes",
    "compute_surface_accessibility_approx",
    "predict_b_cell_epitopes",
]
