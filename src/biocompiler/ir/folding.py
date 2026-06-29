"""
BioCompiler IR Folding Oracle — L3→L4 Structure Prediction
==========================================================

Integrates ESMFold (or fallback) to predict 3D protein structure from
sequence.  The folding oracle is treated as an EXTERNAL TOOL (SLOT-style)
— its output is evidence, not a proven property.

Architecture
------------
  - If ESMFold is installed (``esm`` package + torch) AND the ESM Atlas
    API is reachable: use it for real structure prediction.
  - If ESMFold is unavailable (offline, missing ``esm``, OOM on long
    sequences, etc.): fall back to a pure-Python heuristic that runs
    Chou-Fasman secondary-structure prediction and produces a
    calibrated low-confidence pLDDT estimate.  No 3-D coordinates are
    produced in fallback mode.
  - The :class:`FoldingResult` records which oracle was used and its
    confidence, so downstream consumers (predicates, certificates) can
    treat the output appropriately.

This module is the *only* place where the L3→L4 lowering pass touches
the ESMFold engine; :func:`biocompiler.ir.passes.fold` simply calls
:func:`fold_sequence` and packages the result into an
:class:`IR_L4_FoldedProtein`.

Sequence handling
-----------------
The polypeptide sequence carried by :class:`IR_L3_Polypeptide` always
ends with ``"*"`` (the stop codon — see
:func:`biocompiler.ir.passes.translate`).  Folding operates on the
actual residues, so the trailing ``"*"`` is stripped before being handed
to either oracle.  Any non-standard amino-acid codes (``X`` for unknown,
``B``/``Z``/``J`` for ambiguous, ...) are tolerated by the heuristic
fallback (which substitutes default propensities) and silently dropped
by the ESMFold path (the ESMFold validator rejects them, triggering the
fallback).

Length limits
-------------
ESMFold's memory and runtime grow roughly cubically with sequence
length; sequences longer than :data:`ESMFOLD_MAX_LENGTH` (1000 residues)
are routed directly to the heuristic fallback to avoid OOM.  The
heuristic itself has no length limit (it is O(n) with a small constant).

References
----------
- Lin et al., Science 2023; 379:1043 (ESMFold / ESM-2)
- Chou & Fasman, Biochemistry 1974; 13:222 (secondary-structure
  propensity tables)
- Jumper et al., Nature 2021; 596:583 (AlphaFold2 pLDDT methodology)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Optional

from .types import IRError, IRLevel

logger = logging.getLogger(__name__)

__all__ = [
    "FoldingResult",
    "fold_sequence",
    "ESMFOLD_MAX_LENGTH",
    "STANDARD_AMINO_ACIDS",
]

# ────────────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────────────

#: Maximum sequence length that will be sent to ESMFold.  Longer
#: sequences are routed directly to the heuristic fallback — ESMFold's
#: memory and runtime grow roughly cubically, and the ESM Atlas API
#: itself refuses sequences longer than ~1000 residues.
ESMFOLD_MAX_LENGTH: int = 1000

#: The 20 canonical amino-acid single-letter codes.  Used to decide
#: whether a sequence is eligible for ESMFold (non-standard codes such
#: as ``X``, ``B``, ``Z``, ``J``, ``*`` are not).
STANDARD_AMINO_ACIDS: frozenset = frozenset("ACDEFGHIKLMNPQRSTVWY")

#: Amino-acid codes that may appear in an :class:`IR_L3_Polypeptide`
#: sequence besides the 20 canonical ones.  ``*`` is the stop codon
#: (always trailing), ``X`` is "unknown" (from an ``N``-containing
#: codon), and ``B``/``Z``/``J`` are IUPAC ambiguity codes (Asx/Glx/
#: either Leu or Ile).  All are stripped or substituted before being
#: passed to ESMFold.
_NONSTANDARD_AA_CODES: frozenset = frozenset("*XBZJ")


# ────────────────────────────────────────────────────────────────────
# Result type
# ────────────────────────────────────────────────────────────────────


@dataclass
class FoldingResult:
    """Outcome of a folding-oracle invocation.

    Attributes
    ----------
    sequence:
        The amino-acid sequence that was folded (with the trailing
        ``"*"`` stop codon stripped).  Length matches
        ``len(coordinates)`` when ``coordinates`` is not ``None``.
    coordinates:
        Per-residue 3-D coordinates as a list of ``[x, y, z]`` triples
        (in Ångströms, ESMFold convention).  ``None`` when the oracle
        could not produce 3-D structure (heuristic fallback, offline
        mode, or ESMFold failure).
    confidence:
        Mean per-residue confidence on the pLDDT scale (0–100).  For
        ESMFold this is the mean pLDDT; for the heuristic fallback it
        is the calibrated heuristic estimate (always ≤ 55).  ``0.0``
        when the oracle produced no prediction at all.
    oracle_used:
        Which oracle produced this result.  One of:
          - ``"esmfold"``  — real ESMFold (API or local ``esm``).
          - ``"fallback"`` — heuristic Chou-Fasman / propensity-based.
          - ``"none"``     — no oracle succeeded (complete failure).
    secondary_structure:
        DSSP-style per-residue secondary-structure string using the
        three-character alphabet ``H`` (helix) / ``E`` (sheet) /
        ``C`` (coil or turn).  Same length as ``sequence``.
    metadata:
        Free-form per-fold metadata (engine name, model version,
        execution time, per-residue pLDDT list, ...).
    """

    sequence: str
    coordinates: Optional[list[list[float]]] = None
    confidence: float = 0.0
    oracle_used: str = "none"
    secondary_structure: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)


# ────────────────────────────────────────────────────────────────────
# Public entry point
# ────────────────────────────────────────────────────────────────────


def fold_sequence(
    protein_seq: str,
    use_esmfold: bool = True,
) -> FoldingResult:
    """Predict the 3-D structure of a polypeptide.

    Tries ESMFold first (when ``use_esmfold`` is true and the sequence
    is eligible), then falls back to the pure-Python Chou-Fasman
    heuristic.  Never raises on prediction failure — instead returns a
    :class:`FoldingResult` with ``oracle_used="none"`` so callers can
    decide how to handle the missing structure.

    Parameters
    ----------
    protein_seq:
        Amino-acid sequence (single-letter codes).  May end with
        ``"*"`` (the stop codon), which is stripped before folding.
        May contain ``X`` / ``B`` / ``Z`` / ``J`` (non-standard codes),
        which trigger the heuristic fallback.
    use_esmfold:
        If true (default), attempt ESMFold first.  If false, skip
        straight to the heuristic — useful for offline tests or when
        the caller knows the ESM Atlas API is unreachable.

    Returns
    -------
    FoldingResult
        Always non-``None``.  The ``oracle_used`` field tells the
        caller which backend produced the result.

    Raises
    ------
    IRError
        If the input is empty or consists entirely of non-residue
        characters (e.g. just ``"*"``).
    """
    # ── Normalise input ────────────────────────────────────────────
    if not isinstance(protein_seq, str):
        raise IRError(
            IRLevel.L3,
            f"protein sequence must be a string, got {type(protein_seq).__name__}",
        )
    seq = protein_seq.upper().strip()

    # Strip the trailing stop-codon marker — it is not a residue and
    # ESMFold / the heuristic both reject it.
    if seq.endswith("*"):
        seq = seq[:-1]

    if not seq:
        raise IRError(
            IRLevel.L3,
            "cannot fold an empty polypeptide (sequence is empty or only '*')",
        )

    # ── Decide which oracle to try ────────────────────────────────
    has_nonstandard = bool(set(seq) - STANDARD_AMINO_ACIDS)
    too_long = len(seq) > ESMFOLD_MAX_LENGTH

    if use_esmfold and too_long:
        logger.info(
            "Sequence length %d exceeds ESMFold limit %d; using heuristic fallback",
            len(seq), ESMFOLD_MAX_LENGTH,
        )

    if use_esmfold and not has_nonstandard and not too_long:
        result = _fold_with_esmfold(seq)
        if result is not None:
            return result
        # ESMFold unavailable or failed — fall through to heuristic.
        logger.info(
            "ESMFold unavailable for %d-aa sequence; falling back to heuristic",
            len(seq),
        )

    return _fold_heuristic(seq)


# ────────────────────────────────────────────────────────────────────
# Backend 1: ESMFold
# ────────────────────────────────────────────────────────────────────


def _fold_with_esmfold(seq: str) -> Optional[FoldingResult]:
    """Attempt a real ESMFold structure prediction.

    Returns ``None`` (signalling "use the heuristic instead") when:
      - The ``biocompiler.engines.esmfold`` module cannot be imported
        (shouldn't happen in normal installs, but defend against
        partial installs / circular imports).
      - ESMFold itself reports failure (API unreachable, local
        ``esm`` missing, ...).
      - ESMFold internally fell back to its own heuristic — in that
        case we prefer our own heuristic so the ``oracle_used``
        label is honest.

    Returns a populated :class:`FoldingResult` (with ``oracle_used``
    set to ``"esmfold"``) only when a *real* ESMFold prediction
    (API or local) was produced.

    Never raises — all exceptions are caught and logged, and ``None``
    is returned so the caller can fall back to the heuristic.
    """
    try:
        # Local import — keeps the module importable even if the
        # engines package is partially broken (e.g. optional deps
        # like ``numpy`` are missing).
        from ..engines.esmfold import (
            predict_structure,
            parse_pdb,
            is_esmfold_available,
        )
    except ImportError as exc:
        logger.debug("ESMFold engine not importable: %s", exc)
        return None

    # Cheap pre-check: if neither the API nor the local ``esm`` package
    # is available, skip the (potentially slow) call entirely.  This
    # keeps offline test runs fast.
    try:
        if not is_esmfold_available():
            logger.debug("ESMFold not available (API + local esm both missing)")
            return None
    except Exception as exc:  # pragma: no cover - defensive
        logger.debug("is_esmfold_available() raised: %s", exc)
        return None

    try:
        result = predict_structure(seq, use_api=True)
    except Exception as exc:
        # ESMFold can raise ESMFoldError on invalid input — but we
        # already filtered non-standard AAs in ``fold_sequence``, so
        # this is most likely a transient API/network error.
        logger.warning("ESMFold predict_structure raised %s: %s", type(exc).__name__, exc)
        return None

    # ESMFold returns success=False on complete failure (heuristic
    # itself raised).  Signal "use heuristic" by returning None — the
    # caller will then call ``_fold_heuristic`` and produce an honest
    # "fallback" or "none" result.
    if not result.success:
        logger.warning("ESMFold prediction failed: %s", result.error)
        return None

    # If ESMFold internally used its own heuristic fallback, defer to
    # our heuristic so the ``oracle_used`` label is accurate.  The
    # engine-level heuristic and our IR-level heuristic share the same
    # underlying Chou-Fasman implementation, so no information is lost.
    if result.method == "heuristic_fallback":
        logger.debug("ESMFold used its own heuristic; deferring to IR-level heuristic")
        return None

    # ── Real ESMFold result — extract coordinates + SS string ──────
    coordinates: Optional[list[list[float]]] = None
    if result.pdb_string:
        try:
            parsed = parse_pdb(result.pdb_string)
            ca_coords: list[list[float]] = []
            for residue in parsed["residues"]:
                ca = residue.get("ca")
                if ca is not None:
                    # ``ca`` is a (x, y, z) tuple; promote to list for
                    # JSON-friendliness downstream.
                    ca_coords.append([float(ca[0]), float(ca[1]), float(ca[2])])
            if ca_coords:
                coordinates = ca_coords
        except Exception as exc:
            logger.warning("PDB parsing failed: %s", exc)
            coordinates = None

    # ESMFold does not return a DSSP-style SS string directly.  We
    # derive one from the same Chou-Fasman code the fallback uses —
    # this keeps the SS annotation consistent across both backends
    # and avoids pulling in a heavier DSSP dependency.
    secondary_structure = _predict_ss_string(seq)

    return FoldingResult(
        sequence=seq,
        coordinates=coordinates,
        confidence=float(result.primary_score),  # mean pLDDT (0-100)
        oracle_used="esmfold",
        secondary_structure=secondary_structure,
        metadata={
            "engine": "esmfold",
            "method": result.method,
            "model_name": result.model_name,
            "plddt_scores": list(result.plddt_scores),
            "pae_matrix": result.pae_matrix,
            "classification": result.classification,
            "execution_time_s": result.execution_time_s,
            "has_coordinates": coordinates is not None,
        },
    )


# ────────────────────────────────────────────────────────────────────
# Backend 2: pure-Python heuristic (Chou-Fasman)
# ────────────────────────────────────────────────────────────────────


def _fold_heuristic(seq: str) -> FoldingResult:
    """Predict secondary structure + a calibrated pLDDT estimate.

    Uses the existing
    :func:`biocompiler.engines.esmfold_fallback.predict_structure_heuristic`
    implementation, which runs the full Chou-Fasman algorithm
    (helix/sheet/turn nucleation, Pro/Gly breaking, overlap resolution)
    plus hydrophobicity / charge / contact-density modulation.

    No 3-D coordinates are produced — the heuristic cannot predict
    tertiary structure, only per-residue secondary-structure labels
    and a low-confidence pLDDT estimate.

    Returns a :class:`FoldingResult` with ``oracle_used="fallback"``.
    Never raises — if the heuristic itself fails (which would indicate
    a bug in the fallback module), returns a :class:`FoldingResult`
    with ``oracle_used="none"`` and ``confidence=0.0``.
    """
    try:
        from ..engines.esmfold_fallback import (
            predict_structure_heuristic,
            estimate_secondary_structure_from_sequence,
        )
    except ImportError as exc:
        logger.error("Heuristic fallback module not importable: %s", exc)
        return FoldingResult(
            sequence=seq,
            coordinates=None,
            confidence=0.0,
            oracle_used="none",
            secondary_structure="C" * len(seq),
            metadata={
                "engine": "none",
                "error": f"heuristic module not importable: {exc}",
            },
        )

    # Secondary-structure string (DSSP-style H/E/C).  Computed first
    # because the heuristic estimator may not return a usable
    # ``ss_prediction`` field for very short sequences.
    try:
        ss_estimate = estimate_secondary_structure_from_sequence(seq)
        ss_string = ss_estimate.ss_string or ("C" * len(seq))
    except Exception as exc:
        logger.warning("SS estimation failed: %s", exc)
        ss_string = "C" * len(seq)

    try:
        heuristic = predict_structure_heuristic(seq)
    except Exception as exc:
        logger.error("Heuristic structure prediction raised: %s", exc)
        return FoldingResult(
            sequence=seq,
            coordinates=None,
            confidence=0.0,
            oracle_used="none",
            secondary_structure=ss_string,
            metadata={
                "engine": "none",
                "error": f"heuristic prediction raised: {exc}",
            },
        )

    # The heuristic returns ``plddt_scores`` per residue and a mean
    # pLDDT clamped to [HEURISTIC_MIN_CONFIDENCE, HEURISTIC_MAX_CONFIDENCE].
    plddt_scores: list[float] = list(heuristic.get("plddt_scores", []))
    if plddt_scores and len(plddt_scores) != len(seq):
        # Length mismatch (shouldn't happen, but defend against it) —
        # fall back to the scalar mean for all residues.
        logger.warning(
            "Heuristic pLDDT length %d != sequence length %d; using scalar mean",
            len(plddt_scores), len(seq),
        )
        plddt_scores = [heuristic.get("mean_plddt", 0.0)] * len(seq)
    elif not plddt_scores:
        plddt_scores = [heuristic.get("mean_plddt", 0.0)] * len(seq)

    mean_plddt = float(heuristic.get("mean_plddt", 0.0))

    # Prefer the dedicated SS string from the heuristic if present
    # (it includes the same Chou-Fasman assignments we computed above
    # but may differ in edge cases — the heuristic's version is
    # authoritative).
    final_ss = heuristic.get("ss_prediction") or ss_string
    if len(final_ss) != len(seq):
        # Defensive: pad/truncate to match the sequence length so the
        # IR-L4 metadata stays consistent.
        if len(final_ss) < len(seq):
            final_ss = final_ss + "C" * (len(seq) - len(final_ss))
        else:
            final_ss = final_ss[: len(seq)]

    # Pull out the SS fraction breakdown for downstream introspection.
    # The heuristic stores this either as a dict (current shape from
    # ``predict_structure_heuristic``) or as a ``SecondaryStructureEstimate``
    # dataclass — handle both forms defensively.
    ss_obj = heuristic.get("secondary_structure")
    ss_fractions: dict[str, float] = {}
    if ss_obj is not None:
        if isinstance(ss_obj, dict):
            ss_fractions = {
                "helix_fraction": float(ss_obj.get("helix_fraction", 0.0)),
                "sheet_fraction": float(ss_obj.get("sheet_fraction", 0.0)),
                "turn_fraction": float(ss_obj.get("turn_fraction", 0.0)),
                "coil_fraction": float(ss_obj.get("coil_fraction", 0.0)),
            }
        else:
            ss_fractions = {
                "helix_fraction": float(getattr(ss_obj, "helix_fraction", 0.0)),
                "sheet_fraction": float(getattr(ss_obj, "sheet_fraction", 0.0)),
                "turn_fraction": float(getattr(ss_obj, "turn_fraction", 0.0)),
                "coil_fraction": float(getattr(ss_obj, "coil_fraction", 0.0)),
            }

    return FoldingResult(
        sequence=seq,
        coordinates=None,  # heuristic does not produce 3-D coords
        confidence=mean_plddt,
        oracle_used="fallback",
        secondary_structure=final_ss,
        metadata={
            "engine": "heuristic",
            "method": heuristic.get("method", "heuristic_fallback"),
            "model_name": heuristic.get("model_name", "heuristic_v2"),
            "plddt_scores": plddt_scores,
            "confidence_factor": float(heuristic.get("confidence", 0.0)),
            "ss_fractions": ss_fractions,
            "has_coordinates": False,
        },
    )


# ────────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────────


def _predict_ss_string(seq: str) -> str:
    """Return a DSSP-style secondary-structure string for ``seq``.

    Thin wrapper around
    :func:`biocompiler.engines.esmfold_fallback.estimate_secondary_structure_from_sequence`
    that always returns a string of the same length as ``seq`` (padding
    with ``"C"`` on any failure).  Used by :func:`_fold_with_esmfold`
    to attach an SS annotation to ESMFold results (which do not include
    one directly).
    """
    try:
        from ..engines.esmfold_fallback import (
            estimate_secondary_structure_from_sequence,
        )
    except ImportError:
        return "C" * len(seq)

    try:
        ss = estimate_secondary_structure_from_sequence(seq)
        return ss.ss_string or ("C" * len(seq))
    except Exception:
        return "C" * len(seq)
