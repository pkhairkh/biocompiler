"""
BioCompiler Stability Predicates
=================================
Type-system predicates for protein stability assessment.

Evaluates thermodynamic stability, mutational impact, disulfide bond
integrity, and hydrophobic core quality.  Each predicate returns a
``TypeCheckResult`` with a five-valued ``Verdict`` and, where possible,
a structured ``derivation`` explaining the evidence chain.

When a PDB structure string is supplied the predicates use simple
geometric checks (CB-CB distances for disulfides, SASA proxies for
core quality).  Without a structure, they fall back to empirical
sequence-based heuristics (``estimate_stability_empirical``) and flag
the result with a ``knowledge_gap`` note.
"""

from __future__ import annotations

import math
import logging
from typing import Optional

from biocompiler.type_system import Verdict, TypeCheckResult
from .constants import HYDROPHOBIC_AAS

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# Amino-acid sets
# ────────────────────────────────────────────────────────────
CHARGED_AAS: set[str] = {"K", "R", "H", "D", "E"}

# Normal hydrophobic-fraction range for well-folded globular proteins
_HYDRO_FRAC_LO = 0.30
_HYDRO_FRAC_HI = 0.45

# Disulfide-bond CB-CB distance threshold (Angstroms)
_DISULFIDE_CB_DIST_THRESHOLD = 6.5


# ────────────────────────────────────────────────────────────
# Helper: PDB coordinate extraction
# ────────────────────────────────────────────────────────────
def _parse_pdb_coords(pdb_string: str) -> dict[int, dict[str, list[float]]]:
    """Extract CB (or CA for Gly) coordinates from a PDB string.

    Returns:
        Mapping of residue index (1-based) to atom-name -> [x, y, z].
    """
    coords: dict[int, dict[str, list[float]]] = {}
    for line in pdb_string.splitlines():
        if not line.startswith("ATOM") and not line.startswith("HETATM"):
            continue
        if len(line) < 54:
            continue
        try:
            atom_name = line[12:16].strip()
            res_seq = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except (ValueError, IndexError):
            continue
        # Keep only CB and CA (CA used as fallback for glycine)
        if atom_name in ("CB", "CA"):
            coords.setdefault(res_seq, {})[atom_name] = [x, y, z]
    return coords


def _get_cb_coords(
    pdb_coords: dict[int, dict[str, list[float]]],
    res_idx: int,
) -> list[float] | None:
    """Return CB coords for *res_idx*, falling back to CA (Gly)."""
    atoms = pdb_coords.get(res_idx)
    if atoms is None:
        return None
    if "CB" in atoms:
        return atoms["CB"]
    if "CA" in atoms:
        return atoms["CA"]
    return None


def _euclidean(a: list[float], b: list[float]) -> float:
    """Euclidean distance between two 3-D points."""
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


# ────────────────────────────────────────────────────────────
# Helper functions
# ────────────────────────────────────────────────────────────
def compute_hydrophobic_fraction(protein: str) -> float:
    """Fraction of AILMFWV residues in *protein*.

    Args:
        protein: Amino-acid sequence (1-letter codes).

    Returns:
        Hydrophobic fraction in [0, 1].  Returns 0.0 for empty input.
    """
    if not protein:
        return 0.0
    protein = protein.upper()
    hydro_count = sum(1 for aa in protein if aa in HYDROPHOBIC_AAS)
    return hydro_count / len(protein)


def estimate_stability_empirical(protein: str) -> dict:
    """Quick empirical stability estimate without FoldX.

    Based on four compositional features:
    1. **Hydrophobic fraction** -- optimal around 0.35; deviations weaken
       the hydrophobic core.
    2. **Charge balance** -- equal positive/negative counts favour salt
       bridges.
    3. **Proline / glycine content** -- excess destabilises regular
       secondary structure.
    4. **Cysteine pairs** -- each disulfide bond contributes roughly
       -3 kcal/mol.

    Args:
        protein: Amino-acid sequence (1-letter codes).

    Returns:
        Dict with keys ``dg_estimate`` (float, kcal/mol), ``confidence``
        (``"low"`` or ``"medium"``), and ``components`` (detailed dict).
    """
    protein = protein.upper()
    n = len(protein)
    if n == 0:
        return {"dg_estimate": 0.0, "confidence": "low", "components": {}}

    hydro_frac = compute_hydrophobic_fraction(protein)

    positive = sum(1 for aa in protein if aa in {"K", "R", "H"})
    negative = sum(1 for aa in protein if aa in {"D", "E"})
    total_charged = positive + negative
    charge_balance = (
        abs(positive - negative) / total_charged if total_charged else 1.0
    )

    proline_frac = protein.count("P") / n
    glycine_frac = protein.count("G") / n
    cys_count = protein.count("C")
    cys_pairs = cys_count // 2

    components = {
        "hydrophobic_fraction": round(hydro_frac, 4),
        "positive_charges": positive,
        "negative_charges": negative,
        "charge_balance": round(charge_balance, 4),
        "proline_fraction": round(proline_frac, 4),
        "glycine_fraction": round(glycine_frac, 4),
        "cysteine_pairs": cys_pairs,
    }

    # --- dG estimate ---
    # Hydrophobic core: peak stability near frac=0.35
    hydro_contribution = -20.0 * (1.0 - abs(hydro_frac - 0.35) / 0.35)

    # Salt bridges: each balanced pair contributes ~-1.5 kcal/mol
    salt_bridge_contribution = (
        -1.5 * min(positive, negative) * (1.0 - charge_balance)
    )

    # Disulfide bonds
    disulfide_contribution = -3.0 * cys_pairs

    # Proline / glycine penalty (destabilise regular secondary structure)
    pro_gly_frac = proline_frac + glycine_frac
    pro_gly_penalty = 15.0 * max(0.0, pro_gly_frac - 0.10)

    # Conformational entropy penalty (longer chains have larger unfolding
    # entropy gain, partially offset by more contacts)
    entropy_penalty = 0.05 * n

    dg_estimate = (
        hydro_contribution
        + salt_bridge_contribution
        + disulfide_contribution
        - pro_gly_penalty
        + entropy_penalty
    )

    # Confidence: "medium" if composition is within normal ranges
    if _HYDRO_FRAC_LO <= hydro_frac <= _HYDRO_FRAC_HI and pro_gly_frac <= 0.12:
        confidence = "medium"
    else:
        confidence = "low"

    return {
        "dg_estimate": round(dg_estimate, 2),
        "confidence": confidence,
        "components": components,
    }


# ────────────────────────────────────────────────────────────
# Predicate 1: Stable Folding
# ────────────────────────────────────────────────────────────
def evaluate_stable_folding(
    sequence: str,
    protein: str,
    organism: str,
    stability_threshold: float = -5.0,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check if predicted protein structure is thermodynamically stable.

    Verdict logic (dG in kcal/mol):
    - ``PASS``        -- dG < *stability_threshold* (clearly stable)
    - ``LIKELY_PASS`` -- *stability_threshold* <= dG < *stability_threshold*/2
    - ``UNCERTAIN``   -- *stability_threshold*/2 <= dG < 0
    - ``LIKELY_FAIL`` -- 0 <= dG < 5
    - ``FAIL``        -- dG >= 5 (clearly unstable)

    If *pdb_string* is provided, a structure-based assessment is used
    (via the foldx module's empirical_stability, which detects FoldX
    availability).  Without a structure the empirical estimator is used
    and a ``knowledge_gap`` note is attached.

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (1-letter codes).
        organism: Target organism name.
        stability_threshold: dG threshold for PASS (default -5.0 kcal/mol).
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict and dG derivation.
    """
    protein = protein.upper()

    if not protein:
        return TypeCheckResult(
            predicate="StableFolding",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
        )

    # Determine dG and method
    if pdb_string is not None:
        try:
            from .foldx import empirical_stability as foldx_empirical

            result = foldx_empirical(protein)
            dg = result.stability_kcal
            method = "empirical_structure_aware"
        except Exception:
            est = estimate_stability_empirical(protein)
            dg = est["dg_estimate"]
            method = "empirical_fallback"
    else:
        est = estimate_stability_empirical(protein)
        dg = est["dg_estimate"]
        method = "empirical"
        confidence = est["confidence"]

    # Derive verdict
    if dg < stability_threshold:
        verdict = Verdict.PASS
        violation = None
    elif dg < stability_threshold / 2.0:
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif dg < 0:
        verdict = Verdict.UNCERTAIN
        violation = (
            f"Marginal stability: dG={dg:.2f} kcal/mol "
            f"(>= {stability_threshold / 2.0:.2f})"
        )
    elif dg < 5.0:
        verdict = Verdict.LIKELY_FAIL
        violation = f"Predicted unstable: dG={dg:.2f} kcal/mol (>= 0)"
    else:
        verdict = Verdict.FAIL
        violation = f"Strongly unstable: dG={dg:.2f} kcal/mol (>= 5.0)"

    # Build derivation
    derivation = [
        {"step": "dg_estimate", "value": round(dg, 2), "unit": "kcal/mol"},
        {"step": "method", "value": method},
        {"step": "threshold", "value": stability_threshold, "unit": "kcal/mol"},
    ]
    if pdb_string is None:
        derivation.append({"step": "confidence", "value": confidence})

    knowledge_gap = None
    if pdb_string is None:
        knowledge_gap = (
            "No PDB structure provided; stability estimated from sequence "
            "composition only.  Structural analysis would improve confidence."
        )

    return TypeCheckResult(
        predicate=f"StableFolding({stability_threshold})",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )


# ────────────────────────────────────────────────────────────
# Predicate 2: No Destabilizing Mutation
# ────────────────────────────────────────────────────────────
def evaluate_no_destabilizing_mutation(
    sequence: str,
    protein: str,
    organism: str,
    max_ddg: float = 3.0,
    pdb_string: str | None = None,
    original_protein: str | None = None,
) -> TypeCheckResult:
    """Check if codon-optimized sequence introduces destabilizing mutations.

    If *original_protein* is supplied, every position where the proteins
    differ is checked for ddG impact using a BLOSUM62-based heuristic.
    Positions with ddG > *max_ddg* are flagged as destabilizing.

    When no *original_protein* is given the predicate returns PASS
    (no mutations to evaluate).

    Args:
        sequence: DNA coding sequence (optimised).
        protein: Amino-acid sequence of the optimised construct.
        organism: Target organism name.
        max_ddg: Maximum acceptable ddG (kcal/mol, default 3.0).
        pdb_string: Optional PDB-format structure string.
        original_protein: Original (wild-type) amino-acid sequence.

    Returns:
        TypeCheckResult with verdict and per-position ddG derivation.
    """
    protein = protein.upper()

    if original_protein is None:
        return TypeCheckResult(
            predicate=f"NoDestabilizingMutation({max_ddg})",
            verdict=Verdict.PASS,
            derivation=[{"step": "no_original_protein", "value": True}],
        )

    original_protein = original_protein.upper()

    if len(protein) != len(original_protein):
        return TypeCheckResult(
            predicate=f"NoDestabilizingMutation({max_ddg})",
            verdict=Verdict.FAIL,
            violation=(
                f"Protein length mismatch: optimised={len(protein)}, "
                f"original={len(original_protein)}"
            ),
        )

    # Scan for mutations and estimate ddG
    from .type_system import BLOSUM62

    destabilizing_positions: list[dict] = []
    all_mutations: list[dict] = []

    for i, (orig_aa, new_aa) in enumerate(zip(original_protein, protein)):
        if orig_aa == new_aa:
            continue

        # BLOSUM62 score: negative means unlikely substitution
        blosum = BLOSUM62.get((orig_aa, new_aa), -10)

        # Rough ddG estimate: each BLOSUM62 unit ~ 0.5-1.0 kcal/mol
        # Conserved substitutions have small ddG; radical ones are large.
        ddg_estimate = -blosum * 0.8  # negative BLOSUM -> positive ddG

        mutation_info = {
            "position": i,
            "original_aa": orig_aa,
            "new_aa": new_aa,
            "blosum62": blosum,
            "ddg_estimate": round(ddg_estimate, 2),
        }
        all_mutations.append(mutation_info)

        if ddg_estimate > max_ddg:
            destabilizing_positions.append(mutation_info)

    # Verdict
    if not all_mutations:
        return TypeCheckResult(
            predicate=f"NoDestabilizingMutation({max_ddg})",
            verdict=Verdict.PASS,
            derivation=[{"step": "no_mutations", "value": True}],
        )

    if not destabilizing_positions:
        verdict = Verdict.PASS
        violation = None
    elif len(destabilizing_positions) == 1:
        verdict = Verdict.LIKELY_FAIL
        pos = destabilizing_positions[0]
        violation = (
            f"Destabilizing mutation at position {pos['position']}: "
            f"{pos['original_aa']}->{pos['new_aa']}, "
            f"ddG~{pos['ddg_estimate']:.2f} > {max_ddg}"
        )
    else:
        verdict = Verdict.FAIL
        positions_str = ", ".join(
            f"{p['position']}({p['original_aa']}->{p['new_aa']}, "
            f"ddG~{p['ddg_estimate']:.2f})"
            for p in destabilizing_positions
        )
        violation = (
            f"{len(destabilizing_positions)} destabilizing mutations: "
            f"{positions_str}"
        )

    derivation = [
        {"step": "total_mutations", "value": len(all_mutations)},
        {"step": "destabilizing_count", "value": len(destabilizing_positions)},
        {"step": "max_ddg_threshold", "value": max_ddg},
    ]
    if all_mutations:
        worst = max(all_mutations, key=lambda m: m["ddg_estimate"])
        derivation.append({"step": "worst_ddg", "value": worst["ddg_estimate"]})

    knowledge_gap = None
    if pdb_string is None:
        knowledge_gap = (
            "ddG estimated from BLOSUM62 heuristic; structure-based "
            "FoldX analysis would improve accuracy."
        )

    return TypeCheckResult(
        predicate=f"NoDestabilizingMutation({max_ddg})",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )


# ────────────────────────────────────────────────────────────
# Predicate 3: Disulfide Bond Integrity
# ────────────────────────────────────────────────────────────
def evaluate_disulfide_bond_integrity(
    sequence: str,
    protein: str,
    organism: str,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check disulfide bond integrity from cysteine count and pairing.

    Sequence-level check:
    - Even number of cysteines (including 0) -> PASS
    - Odd number -> LIKELY_FAIL (unpaired Cys may form incorrect bonds)

    If *pdb_string* is provided, additionally verify that Cys residues
    are spatially close enough for disulfide bond formation (CB-CB
    distance < 6.5 A).

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (1-letter codes).
        organism: Target organism name.
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict and cysteine-pairing derivation.
    """
    protein = protein.upper()
    cys_count = protein.count("C")
    cys_positions = [i for i, aa in enumerate(protein) if aa == "C"]

    derivation: list[dict] = [
        {"step": "cysteine_count", "value": cys_count},
        {"step": "cysteine_positions", "value": cys_positions},
    ]

    # Zero cysteines: no disulfide bonds needed
    if cys_count == 0:
        return TypeCheckResult(
            predicate="DisulfideBondIntegrity",
            verdict=Verdict.PASS,
            derivation=derivation,
        )

    # Odd number: unpaired Cys
    if cys_count % 2 != 0:
        violation = (
            f"Odd number of cysteines ({cys_count}): at least one "
            f"unpaired Cys that may form incorrect disulfides"
        )
        derivation.append({"step": "paired", "value": False})
        derivation.append({"step": "reason", "value": "odd_count"})

        knowledge_gap = (
            "Unpaired Cys may still be buried and harmless; "
            "structural SASA analysis would clarify."
            if pdb_string is not None
            else "Cannot determine if unpaired Cys is buried without "
            "structural data."
        )

        return TypeCheckResult(
            predicate="DisulfideBondIntegrity",
            verdict=Verdict.LIKELY_FAIL,
            derivation=derivation,
            violation=violation,
            knowledge_gap=knowledge_gap,
        )

    # Even number: potentially paired.
    # If PDB provided, check spatial proximity of Cys pairs.
    if pdb_string is not None:
        pdb_coords = _parse_pdb_coords(pdb_string)
        unpaired_cys: list[int] = []
        paired_cys: list[tuple[int, int, float]] = []

        remaining = list(cys_positions)
        while len(remaining) >= 2:
            best_pair = None
            best_dist = float("inf")
            for i_idx in range(len(remaining)):
                for j_idx in range(i_idx + 1, len(remaining)):
                    pos_i = remaining[i_idx]
                    pos_j = remaining[j_idx]
                    # PDB residue numbers are 1-based
                    coord_i = _get_cb_coords(pdb_coords, pos_i + 1)
                    coord_j = _get_cb_coords(pdb_coords, pos_j + 1)
                    if coord_i is None or coord_j is None:
                        continue
                    dist = _euclidean(coord_i, coord_j)
                    if dist < best_dist:
                        best_dist = dist
                        best_pair = (i_idx, j_idx, pos_i, pos_j, dist)

            if best_pair is None:
                break

            _, _, pos_i, pos_j, dist = best_pair
            if dist < _DISULFIDE_CB_DIST_THRESHOLD:
                paired_cys.append((pos_i, pos_j, round(dist, 2)))
                remaining.pop(best_pair[1])
                remaining.pop(best_pair[0])
            else:
                break

        unpaired_cys = remaining

        derivation.append({
            "step": "structure_pairs",
            "value": [
                {"positions": [p[0], p[1]], "cb_distance": p[2]}
                for p in paired_cys
            ],
        })
        derivation.append({
            "step": "unpaired_from_structure",
            "value": unpaired_cys,
        })

        if not unpaired_cys:
            return TypeCheckResult(
                predicate="DisulfideBondIntegrity",
                verdict=Verdict.PASS,
                derivation=derivation,
            )
        else:
            violation = (
                f"{len(unpaired_cys)} cysteine(s) cannot form disulfide "
                f"bonds (CB-CB distance > {_DISULFIDE_CB_DIST_THRESHOLD} A): "
                f"positions {unpaired_cys}"
            )
            return TypeCheckResult(
                predicate="DisulfideBondIntegrity",
                verdict=Verdict.LIKELY_FAIL,
                derivation=derivation,
                violation=violation,
            )

    # Even number, no PDB -- assume pairable
    derivation.append({"step": "paired", "value": True})
    derivation.append({"step": "note", "value": "assumed_pairable_no_structure"})

    return TypeCheckResult(
        predicate="DisulfideBondIntegrity",
        verdict=Verdict.PASS,
        derivation=derivation,
        knowledge_gap=(
            "Cysteine count is even but spatial pairing cannot be "
            "verified without structural data."
        ),
    )


# ────────────────────────────────────────────────────────────
# Predicate 4: Hydrophobic Core Quality
# ────────────────────────────────────────────────────────────
def evaluate_hydrophobic_core_quality(
    sequence: str,
    protein: str,
    organism: str,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check hydrophobic core quality from amino-acid composition.

    The fraction of hydrophobic residues (A, I, L, M, F, W, V) should
    fall within the normal range of 0.30-0.45 for well-folded globular
    proteins:

    - **Too low** (< 0.30): insufficient hydrophobic core -> unstable.
    - **Too high** (> 0.45): over-hydrophobic -> aggregation-prone.
    - **In range** (0.30-0.45): PASS.

    If *pdb_string* is provided, additionally check that hydrophobic
    residues are spatially buried (mean CB distance from the protein
    centroid is smaller than the overall mean).

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (1-letter codes).
        organism: Target organism name.
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict and hydrophobic-fraction derivation.
    """
    protein = protein.upper()
    hydro_frac = compute_hydrophobic_fraction(protein)

    derivation: list[dict] = [
        {"step": "hydrophobic_fraction", "value": round(hydro_frac, 4)},
        {"step": "normal_range", "value": [_HYDRO_FRAC_LO, _HYDRO_FRAC_HI]},
    ]

    # --- Sequence-based assessment ---
    if hydro_frac < _HYDRO_FRAC_LO:
        deficit = _HYDRO_FRAC_LO - hydro_frac
        if deficit > 0.10:
            verdict = Verdict.FAIL
        elif deficit > 0.05:
            verdict = Verdict.LIKELY_FAIL
        else:
            verdict = Verdict.UNCERTAIN
        violation = (
            f"Hydrophobic fraction {hydro_frac:.3f} is below normal "
            f"range [{_HYDRO_FRAC_LO}, {_HYDRO_FRAC_HI}] -- insufficient "
            f"hydrophobic core"
        )
    elif hydro_frac > _HYDRO_FRAC_HI:
        excess = hydro_frac - _HYDRO_FRAC_HI
        if excess > 0.10:
            verdict = Verdict.FAIL
        elif excess > 0.05:
            verdict = Verdict.LIKELY_FAIL
        else:
            verdict = Verdict.UNCERTAIN
        violation = (
            f"Hydrophobic fraction {hydro_frac:.3f} is above normal "
            f"range [{_HYDRO_FRAC_LO}, {_HYDRO_FRAC_HI}] -- aggregation-prone"
        )
    else:
        verdict = Verdict.PASS
        violation = None

    # --- Structure-based refinement (if PDB provided) ---
    if pdb_string is not None and verdict in (Verdict.PASS, Verdict.UNCERTAIN):
        pdb_coords = _parse_pdb_coords(pdb_string)
        if pdb_coords:
            all_coords: list[list[float]] = []
            hydro_coords: list[list[float]] = []
            for res_idx, atoms in pdb_coords.items():
                coord = _get_cb_coords(pdb_coords, res_idx)
                if coord is None:
                    continue
                all_coords.append(coord)
                aa_pos = res_idx - 1  # PDB is 1-based, protein is 0-based
                if 0 <= aa_pos < len(protein) and protein[aa_pos] in HYDROPHOBIC_AAS:
                    hydro_coords.append(coord)

            if len(all_coords) >= 3 and len(hydro_coords) >= 1:
                centroid = [
                    sum(c[d] for c in all_coords) / len(all_coords)
                    for d in range(3)
                ]
                mean_all_dist = sum(
                    _euclidean(c, centroid) for c in all_coords
                ) / len(all_coords)
                mean_hydro_dist = sum(
                    _euclidean(c, centroid) for c in hydro_coords
                ) / len(hydro_coords)

                derivation.append({
                    "step": "mean_all_cb_to_centroid",
                    "value": round(mean_all_dist, 2),
                    "unit": "A",
                })
                derivation.append({
                    "step": "mean_hydro_cb_to_centroid",
                    "value": round(mean_hydro_dist, 2),
                    "unit": "A",
                })

                # Hydrophobic residues should be closer to centroid (buried)
                if mean_hydro_dist < mean_all_dist:
                    if verdict == Verdict.UNCERTAIN:
                        verdict = Verdict.LIKELY_PASS
                        violation = (
                            f"Hydrophobic fraction {hydro_frac:.3f} slightly "
                            f"out of range, but hydrophobic residues are "
                            f"buried (mean dist {mean_hydro_dist:.1f} A < "
                            f"overall {mean_all_dist:.1f} A)"
                        )
                else:
                    if verdict == Verdict.PASS:
                        verdict = Verdict.UNCERTAIN
                        violation = (
                            f"Hydrophobic fraction {hydro_frac:.3f} in range, "
                            f"but hydrophobic residues are not buried "
                            f"(mean dist {mean_hydro_dist:.1f} A >= "
                            f"overall {mean_all_dist:.1f} A)"
                        )

    knowledge_gap = None
    if pdb_string is None and verdict != Verdict.FAIL:
        knowledge_gap = (
            "Hydrophobic fraction is within normal range, but core "
            "burial cannot be verified without structural data."
        )

    return TypeCheckResult(
        predicate="HydrophobicCoreQuality",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )
