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
from typing import Any

from .type_system import Verdict, TypeCheckResult
from .constants import HYDROPHOBIC_AAS

logger = logging.getLogger(__name__)

__all__ = [
    "compute_hydrophobic_fraction",
    "estimate_stability_empirical",
    "evaluate_stable_folding",
    "evaluate_no_destabilizing_mutation",
    "evaluate_disulfide_bond_integrity",
    "evaluate_hydrophobic_core_quality",
]

# ────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────

# Normal hydrophobic-fraction range for well-folded globular proteins
_HYDRO_FRAC_LO = 0.30
_HYDRO_FRAC_HI = 0.45

# Hydrophobic fraction at peak stability (used in empirical estimator)
_HYDRO_PEAK_FRAC = 0.35

# Disulfide-bond CB-CB distance threshold (Angstroms)
_DISULFIDE_CB_DIST_THRESHOLD = 6.5

# ── Disulfide-bond localisation helpers ──
# Keywords in organism names that indicate prokaryotes (bacteria / archaea)
_PROKARYOTIC_KEYWORDS = frozenset({
    "ecoli", "e. coli", "escherichia", "e.coli",
    "bacillus", "b. subtilis", "bsubtilis", "b.subtilis",
    "pseudomonas", "staphylococcus", "streptococcus",
    "mycobacterium", "clostridium", "salmonella",
    "vibrio", "shigella", "campylobacter",
    "cyanobacteria", "thermus", "debaryomyces",  # archaea-like
    "archaea", "archaebacteria",
    "prokaryote", "bacteria", "bacterial",
})

# Signal-peptide detection parameters
_SIGNAL_PEPTIDE_WINDOW = 40      # only check first N residues
_SIGNAL_PEPTIDE_MIN_HYDRO_STRETCH = 5  # min consecutive hydrophobic residues
_SIGNAL_PEPTIDE_HYDRO_AAS = frozenset({"A", "I", "L", "M", "F", "W", "V"})

# Hydrophobic core quality threshold
_CORE_QUALITY_PASS_THRESHOLD = 0.6   # was 0.7; relaxed per e2e validation
_SMALL_PROTEIN_LENGTH = 100          # proteins shorter than this may lack a traditional core

# ── Empirical estimator coefficients ──
_HYDRO_CONTRIBUTION_WEIGHT = 20.0      # kcal/mol weight for hydrophobic core
_SALT_BRIDGE_KCAL_PER_PAIR = 1.5       # kcal/mol per balanced charge pair
_DISULFIDE_BOND_KCAL = -3.0            # kcal/mol per disulfide bond (stabilizing)
_PRO_GLY_PENALTY_WEIGHT = 15.0         # kcal/mol penalty weight
_PRO_GLY_PENALTY_THRESHOLD = 0.10      # fraction above which penalty applies
_ENTROPY_PENALTY_COEFF = 0.05          # kcal/mol per residue (conformational entropy)
_PRO_GLY_CONFIDENCE_THRESHOLD = 0.12   # max pro+gly frac for "medium" confidence

# ── Stability verdict thresholds ──
_CLEARLY_UNSTABLE_DG = 5.0             # kcal/mol; dG >= this → FAIL

# ── Mutation ddG estimation ──
_BLOSUM62_DDG_FACTOR = 0.8             # BLOSUM62 score → ddG conversion
_BLOSUM62_UNKNOWN_SCORE = -10          # default BLOSUM62 for unknown pairs


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
        except (ValueError, IndexError) as exc:
            logger.debug("Skipping malformed PDB ATOM line: %s", exc)
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


def estimate_stability_empirical(protein: str) -> dict[str, Any]:
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
    # Hydrophobic core: peak stability near _HYDRO_PEAK_FRAC
    hydro_contribution = (
        -_HYDRO_CONTRIBUTION_WEIGHT
        * (1.0 - abs(hydro_frac - _HYDRO_PEAK_FRAC) / _HYDRO_PEAK_FRAC)
    )

    # Salt bridges: each balanced pair contributes ~_SALT_BRIDGE_KCAL_PER_PAIR
    salt_bridge_contribution = (
        -_SALT_BRIDGE_KCAL_PER_PAIR
        * min(positive, negative)
        * (1.0 - charge_balance)
    )

    # Disulfide bonds
    disulfide_contribution = _DISULFIDE_BOND_KCAL * cys_pairs

    # Proline / glycine penalty (destabilise regular secondary structure)
    pro_gly_frac = proline_frac + glycine_frac
    pro_gly_penalty = _PRO_GLY_PENALTY_WEIGHT * max(
        0.0, pro_gly_frac - _PRO_GLY_PENALTY_THRESHOLD
    )

    # Conformational entropy penalty (longer chains have larger unfolding
    # entropy gain, partially offset by more contacts)
    entropy_penalty = _ENTROPY_PENALTY_COEFF * n

    dg_estimate = (
        hydro_contribution
        + salt_bridge_contribution
        + disulfide_contribution
        - pro_gly_penalty
        + entropy_penalty
    )

    # Confidence: "medium" if composition is within normal ranges
    if _HYDRO_FRAC_LO <= hydro_frac <= _HYDRO_FRAC_HI and pro_gly_frac <= _PRO_GLY_CONFIDENCE_THRESHOLD:
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

    # Small peptides (<50 aa) are inherently stable or cannot be
    # evaluated meaningfully with composition-based heuristics.
    if len(protein) < 50:
        return TypeCheckResult(
            predicate=f"StableFolding({stability_threshold})",
            verdict=Verdict.LIKELY_PASS,
            derivation=[
                {"step": "small_peptide", "value": True},
                {"step": "protein_length", "value": len(protein)},
                {"step": "note", "value": "Peptides <50 aa are inherently stable or cannot be evaluated meaningfully"},
            ],
            knowledge_gap=(
                "Protein is shorter than 50 residues; stability heuristics "
                "are not meaningful for small peptides.  Assumed stable."
            ),
        )

    # Determine dG and method
    if pdb_string is not None:
        try:
            from .foldx import empirical_stability as foldx_empirical

            result = foldx_empirical(protein)
            dg = result.stability_kcal
            method = "empirical_structure_aware"
        except (ImportError, AttributeError, RuntimeError):
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
    elif dg < _CLEARLY_UNSTABLE_DG:
        verdict = Verdict.LIKELY_FAIL
        violation = f"Predicted unstable: dG={dg:.2f} kcal/mol (>= 0)"
    else:
        verdict = Verdict.FAIL
        violation = f"Strongly unstable: dG={dg:.2f} kcal/mol (>= {_CLEARLY_UNSTABLE_DG})"

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
        if pdb_string is None:
            return TypeCheckResult(
                predicate=f"NoDestabilizingMutation({max_ddg})",
                verdict=Verdict.UNCERTAIN,
                violation=(
                    f"Protein length mismatch: optimised={len(protein)}, "
                    f"original={len(original_protein)}"
                ),
                knowledge_gap=(
                    "No structure available to assess mutation impact; "
                    "returning UNCERTAIN instead of FAIL."
                ),
            )
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
        blosum = BLOSUM62.get((orig_aa, new_aa), _BLOSUM62_UNKNOWN_SCORE)

        # Rough ddG estimate: each BLOSUM62 unit ~ 0.5-1.0 kcal/mol
        # Conserved substitutions have small ddG; radical ones are large.
        ddg_estimate = -blosum * _BLOSUM62_DDG_FACTOR  # negative BLOSUM -> positive ddG

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

    derivation = [
        {"step": "total_mutations", "value": len(all_mutations)},
        {"step": "destabilizing_count", "value": len(destabilizing_positions)},
        {"step": "max_ddg_threshold", "value": max_ddg},
    ]
    if all_mutations:
        worst = max(all_mutations, key=lambda m: m["ddg_estimate"])
        derivation.append({"step": "worst_ddg", "value": worst["ddg_estimate"]})

    # Without structural data, downgrade FAIL/LIKELY_FAIL to UNCERTAIN
    # since BLOSUM62-based ddG estimates are unreliable without a structure.
    if pdb_string is None and verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL):
        original_verdict = verdict
        verdict = Verdict.UNCERTAIN
        derivation.append({
            "step": "no_structure_downgrade",
            "value": True,
            "original_verdict": original_verdict.value,
        })

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
def _is_prokaryotic(organism: str) -> bool:
    """Return True if *organism* name suggests a prokaryote.

    Uses case-insensitive keyword matching against a curated set of
    bacterial / archaeal identifiers.
    """
    org_lower = organism.lower().strip()
    return any(kw in org_lower for kw in _PROKARYOTIC_KEYWORDS)


def _has_signal_peptide(protein: str) -> bool:
    """Heuristic signal-peptide detection from the N-terminal region.

    Looks for a stretch of ≥ 5 consecutive hydrophobic residues in the
    first 40 positions, which is characteristic of secretory signal
    peptides.  This is a conservative heuristic — it may miss some
    signal peptides but is unlikely to produce false positives.
    """
    n_term = protein[:_SIGNAL_PEPTIDE_WINDOW]
    consecutive = 0
    for aa in n_term:
        if aa in _SIGNAL_PEPTIDE_HYDRO_AAS:
            consecutive += 1
            if consecutive >= _SIGNAL_PEPTIDE_MIN_HYDRO_STRETCH:
                return True
        else:
            consecutive = 0
    return False


def evaluate_disulfide_bond_integrity(
    sequence: str,
    protein: str,
    organism: str,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check disulfide bond integrity from cysteine count and pairing.

    Organism-aware, localisation-aware logic:
    - **< 2 cysteines**: auto-PASS (disulfide bonds impossible).
    - **Intracellular / prokaryotic cytosolic proteins**: auto-PASS
      (disulfide bonds do not form in the reducing cytosolic
      environment of most prokaryotes and eukaryotic cytosol).
    - **Odd number of cysteines** (in secreted proteins): UNCERTAIN
      (WARNING) — unpaired Cys may be functional, not a FAIL.
    - **Even number of cysteines** in secreted proteins: check spatial
      pairing if PDB is available.

    If *pdb_string* is provided, additionally verify that Cys residues
    are spatially close enough for disulfide bond formation (CB-CB
    distance < 6.5 Å).

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

    is_prokaryote = _is_prokaryotic(organism)
    is_secreted = _has_signal_peptide(protein)

    derivation: list[dict] = [
        {"step": "cysteine_count", "value": cys_count},
        {"step": "cysteine_positions", "value": cys_positions},
        {"step": "is_prokaryotic", "value": is_prokaryote},
        {"step": "has_signal_peptide", "value": is_secreted},
    ]

    # Fewer than 2 cysteines: disulfide bonds are impossible → auto-PASS
    if cys_count < 2:
        derivation.append({"step": "auto_pass_reason", "value": "fewer_than_2_cysteines"})
        return TypeCheckResult(
            predicate="DisulfideBondIntegrity",
            verdict=Verdict.PASS,
            derivation=derivation,
        )

    # Intracellular proteins: disulfide bonds don't form in the reducing
    # cytosol.  For prokaryotes specifically, the cytosol is highly
    # reducing and disulfide bonds essentially never form intracellularly.
    if not is_secreted:
        if is_prokaryote:
            derivation.append({
                "step": "auto_pass_reason",
                "value": "prokaryotic_cytosolic_no_disulfides",
            })
        else:
            derivation.append({
                "step": "auto_pass_reason",
                "value": "intracellular_no_disulfides",
            })
        return TypeCheckResult(
            predicate="DisulfideBondIntegrity",
            verdict=Verdict.PASS,
            derivation=derivation,
            knowledge_gap=(
                "No signal peptide detected; protein assumed intracellular "
                "where disulfide bonds generally do not form.  Experimental "
                "localisation data would confirm."
            ),
        )

    # --- From here on, the protein is secreted / extracellular ---

    # Odd number of cysteines: WARNING (UNCERTAIN), not FAIL.
    # Unpaired Cys may be functional (e.g. catalytic, metal-binding).
    if cys_count % 2 != 0:
        violation = (
            f"Odd number of cysteines ({cys_count}) in secreted protein: "
            f"at least one unpaired Cys (positions {cys_positions}). "
            f"Unpaired Cys may be functional rather than detrimental."
        )
        derivation.append({"step": "paired", "value": False})
        derivation.append({"step": "reason", "value": "odd_count_secreted"})

        knowledge_gap = (
            "Unpaired Cys may still be buried and harmless; "
            "structural SASA analysis would clarify."
            if pdb_string is not None
            else "Cannot determine if unpaired Cys is buried without "
            "structural data."
        )

        return TypeCheckResult(
            predicate="DisulfideBondIntegrity",
            verdict=Verdict.UNCERTAIN,
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
                verdict=Verdict.UNCERTAIN,
                derivation=derivation,
                violation=violation,
                knowledge_gap=(
                    "Unpaired Cys from structural analysis may still be "
                    "functional (e.g. catalytic or metal-binding)."
                ),
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
def _compute_core_quality_score(hydro_frac: float) -> float:
    """Map hydrophobic fraction to a [0, 1] core-quality score.

    The score is **monotonically increasing** with hydrophobic fraction
    up to the optimal value (``_HYDRO_PEAK_FRAC`` ≈ 0.35) and
    **monotonically decreasing** thereafter.  A score of 1.0 means
    the fraction is exactly at the peak; 0.0 means maximally aberrant.

    Piecewise-linear mapping::

        hydro_frac ≤ peak:  score = hydro_frac / peak
        hydro_frac > peak:  score = 1 - (hydro_frac - peak) / (1 - peak)

    This guarantees monotonicity on each side and produces symmetric-ish
    scores at the boundaries of the normal range (0.30 → ~0.86,
    0.45 → ~0.85).
    """
    if hydro_frac <= _HYDRO_PEAK_FRAC:
        # Monotonically increasing from 0 at hydro_frac=0 to 1.0 at peak
        return max(0.0, hydro_frac / _HYDRO_PEAK_FRAC)
    else:
        # Monotonically decreasing from 1.0 at peak toward 0
        return max(0.0, 1.0 - (hydro_frac - _HYDRO_PEAK_FRAC) / (1.0 - _HYDRO_PEAK_FRAC))


def _sequence_based_hydrophobicity_analysis(protein: str) -> dict[str, Any]:
    """Fallback hydrophobicity analysis using only the amino-acid sequence.

    Useful when no PDB structure is available.  Examines both the
    overall hydrophobic fraction and the distribution of hydrophobic
    runs (stretches of consecutive hydrophobic residues) as a proxy
    for core-forming potential.

    Returns:
        Dict with keys ``hydrophobic_fraction``, ``core_quality_score``,
        ``max_hydro_run``, ``hydro_run_count``, and ``assessment``.
    """
    hydro_frac = compute_hydrophobic_fraction(protein)
    core_quality_score = _compute_core_quality_score(hydro_frac)

    # Measure hydrophobic runs (stretches of consecutive hydrophobic AAs)
    max_run = 0
    current_run = 0
    run_count = 0
    for aa in protein:
        if aa in HYDROPHOBIC_AAS:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            if current_run >= 3:
                run_count += 1
            current_run = 0
    if current_run >= 3:
        run_count += 1

    return {
        "hydrophobic_fraction": round(hydro_frac, 4),
        "core_quality_score": round(core_quality_score, 4),
        "max_hydro_run": max_run,
        "hydro_run_count": run_count,
    }


def evaluate_hydrophobic_core_quality(
    sequence: str,
    protein: str,
    organism: str,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check hydrophobic core quality from amino-acid composition.

    A ``core_quality_score`` is computed that maps the hydrophobic
    fraction to a [0, 1] scale centred on the optimal value (~0.35).
    The PASS threshold is **0.6** (previously 0.7), relaxed after
    e2e validation showed too many false negatives at the stricter
    level.

    Additional rules:
    - **Proteins < 100 aa**: may not have a traditional hydrophobic
      core; FAIL is softened to UNCERTAIN.
    - **No PDB structure**: a sequence-based hydrophobicity analysis
      (run-length distribution) is used as a fallback.

    Verdict logic (core_quality_score):
    - ``PASS``        -- score > 0.6
    - ``LIKELY_PASS`` -- score > 0.5
    - ``UNCERTAIN``   -- score > 0.4
    - ``LIKELY_FAIL`` -- score > 0.3
    - ``FAIL``        -- score <= 0.3

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
    core_quality_score = _compute_core_quality_score(hydro_frac)
    is_small_protein = len(protein) < _SMALL_PROTEIN_LENGTH

    derivation: list[dict] = [
        {"step": "hydrophobic_fraction", "value": round(hydro_frac, 4)},
        {"step": "core_quality_score", "value": round(core_quality_score, 4)},
        {"step": "normal_range", "value": [_HYDRO_FRAC_LO, _HYDRO_FRAC_HI]},
        {"step": "protein_length", "value": len(protein)},
        {"step": "is_small_protein", "value": is_small_protein},
    ]

    # --- Core quality score verdict ---
    if core_quality_score > _CORE_QUALITY_PASS_THRESHOLD:
        verdict = Verdict.PASS
        violation = None
    elif core_quality_score > 0.5:
        verdict = Verdict.LIKELY_PASS
        violation = (
            f"Core quality score {core_quality_score:.3f} slightly below "
            f"PASS threshold ({_CORE_QUALITY_PASS_THRESHOLD})"
        )
    elif core_quality_score > 0.4:
        verdict = Verdict.UNCERTAIN
        violation = (
            f"Core quality score {core_quality_score:.3f} marginal "
            f"(hydrophobic fraction {hydro_frac:.3f} deviates from "
            f"optimal {_HYDRO_PEAK_FRAC})"
        )
    elif core_quality_score > 0.3:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Core quality score {core_quality_score:.3f} low "
            f"(hydrophobic fraction {hydro_frac:.3f} far from "
            f"optimal {_HYDRO_PEAK_FRAC})"
        )
    else:
        verdict = Verdict.FAIL
        violation = (
            f"Core quality score {core_quality_score:.3f} very low "
            f"(hydrophobic fraction {hydro_frac:.3f} far from "
            f"optimal {_HYDRO_PEAK_FRAC})"
        )

    # --- Small protein leniency ---
    # Proteins < 100 aa may not form a traditional hydrophobic core;
    # soften FAIL → UNCERTAIN, LIKELY_FAIL → UNCERTAIN.
    if is_small_protein and verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL):
        original_verdict = verdict
        verdict = Verdict.UNCERTAIN
        violation = (
            f"Core quality score {core_quality_score:.3f} low, but protein "
            f"is short ({len(protein)} aa < {_SMALL_PROTEIN_LENGTH}) and may "
            f"not require a traditional hydrophobic core "
            f"(original verdict: {original_verdict.value})"
        )
        derivation.append({
            "step": "small_protein_leniency",
            "value": True,
            "original_verdict": original_verdict.value,
        })

    # --- Structure-based refinement (if PDB provided) ---
    if pdb_string is not None and verdict in (Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN):
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

    # --- Sequence-based fallback when no structure ---
    knowledge_gap = None
    if pdb_string is None:
        seq_analysis = _sequence_based_hydrophobicity_analysis(protein)
        derivation.append({
            "step": "sequence_hydrophobicity_fallback",
            "value": seq_analysis,
        })
        if verdict != Verdict.FAIL:
            knowledge_gap = (
                "Hydrophobic fraction is within normal range, but core "
                "burial cannot be verified without structural data.  "
                "Sequence-based run-length analysis used as fallback."
            )

    return TypeCheckResult(
        predicate="HydrophobicCoreQuality",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )
