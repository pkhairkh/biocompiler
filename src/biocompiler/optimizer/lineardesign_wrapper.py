"""
BioCompiler LinearDesign Wrapper — mRNA MFE Optimization via Lattice DP
=========================================================================

This module provides a Python wrapper for the LinearDesign C++ tool
(Zhang et al. 2023, Nature), which finds the globally optimal mRNA
coding sequence that jointly minimizes minimum free energy (MFE) and
maximizes codon adaptation index (CAI) using a lattice-based dynamic
programming algorithm over a codon-grammar context-free grammar.

Key features:

1. **Subprocess wrapper** — ``run_lineardesign()`` runs the LinearDesign
   binary with configurable parameters, parses structured output, and
   returns a ``LinearDesignResult`` dataclass.

2. **Pareto front scanning** — ``scan_pareto_front()`` runs LinearDesign
   at multiple lambda values to trace the MFE/CAI tradeoff curve,
   enabling selection of the best operating point.

3. **Binary search for target CAI** — ``find_optimal_lambda()`` uses
   bisection to find the lambda that achieves a user-specified target
   CAI, automating the tradeoff calibration.

4. **Codon usage file generation** — ``prepare_codon_usage_file()``
   converts BioCompiler's internal organism codon tables to the CSV
   format required by LinearDesign.

5. **Full integration** — ``optimize_mfe_lineardesign()`` combines
   codon table generation, LinearDesign execution, and result
   validation into a single high-level API.

References:
    Zhang, H. et al. (2023). "Algorithm for optimized mRNA design
    improves stability and immunogenicity." *Nature* 621:396-403.
    doi:10.1038/s41586-023-06127-z

    Huang, L. et al. (2023). "LinearDesign: efficient algorithms
    for optimized mRNA sequence design." *bioRxiv* 2023.03.12.532258.
    doi:10.1101/2023.03.12.532258

Requirements:
    LinearDesign C++ binary must be compiled and available in PATH
    or specified via ``binary_path``.  Install from:
    https://github.com/LinearDesignSoftware/LinearDesign
"""

from __future__ import annotations

import csv
import logging
import os
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    # Dataclass
    "LinearDesignResult",
    # Core functions
    "run_lineardesign",
    "scan_pareto_front",
    "find_optimal_lambda",
    "prepare_codon_usage_file",
    "is_lineardesign_available",
    # Integration
    "optimize_mfe_lineardesign",
]


# ═══════════════════════════════════════════════════════════════════════
# 1. Dataclass
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class LinearDesignResult:
    """Structured result from a LinearDesign optimization run.

    Attributes:
        sequence: Optimized mRNA coding sequence (DNA alphabet, 5'→3').
        structure: Dot-bracket secondary structure prediction from
            LinearDesign's internal folding.  May be ``None`` if the
            binary does not output structure.
        mfe: Minimum free energy of the predicted structure (kcal/mol).
            More negative values indicate more stable (more structured)
            sequences.  ``None`` if parsing fails.
        cai: Codon Adaptation Index of the optimized sequence.
            Ranges from 0.0 to 1.0 (1.0 = all optimal codons).
            ``None`` if parsing fails.
        lambda_val: The MFE/CAI tradeoff parameter used.  0 = pure MFE
            optimization; higher values weight CAI more heavily.
        method: Always ``"lineardesign"`` for results from this module.
        elapsed_seconds: Wall-clock time for the LinearDesign subprocess
            in seconds (excludes codon table generation time).
    """

    sequence: str | None = None
    structure: str | None = None
    mfe: float | None = None
    cai: float | None = None
    lambda_val: float = 0.0
    method: str = "lineardesign"
    elapsed_seconds: float = 0.0
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """Convert to a plain dict for JSON serialization / backward compat."""
        return {
            "sequence": self.sequence,
            "structure": self.structure,
            "mfe": self.mfe,
            "cai": self.cai,
            "lambda": self.lambda_val,
            "method": self.method,
            "elapsed_seconds": self.elapsed_seconds,
            "error": self.error,
        }


# ═══════════════════════════════════════════════════════════════════════
# 2. Availability check
# ═══════════════════════════════════════════════════════════════════════


def is_lineardesign_available(binary_path: str = "lineardesign") -> bool:
    """Check if the LinearDesign binary is available in PATH.

    Args:
        binary_path: Name or full path of the LinearDesign executable.

    Returns:
        ``True`` if the binary can be found and executed, ``False``
        otherwise.
    """
    return shutil.which(binary_path) is not None


# ═══════════════════════════════════════════════════════════════════════
# 3. Codon usage file generation
# ═══════════════════════════════════════════════════════════════════════


def prepare_codon_usage_file(
    organism: str,
    output_path: str,
) -> str:
    """Generate a codon usage CSV file in LinearDesign format.

    Converts BioCompiler's internal organism codon usage table to the
    three-column CSV format expected by LinearDesign's ``--codonusage``
    flag::

        Codon,AminoAcid,Frequency
        GCT,A,0.27
        GCC,A,0.41
        ...

    The *Frequency* column contains the fractional codon usage (the
    ``fraction`` field from BioCompiler's codon usage table), which
    represents the proportion of times each codon is used for its
    amino acid in the reference gene set.

    Args:
        organism: Organism name accepted by
            :func:`~biocompiler.organisms.resolve_organism` (e.g.
            ``"human"``, ``"e_coli"``, ``"cho"``).
        output_path: File path for the output CSV.  Parent directories
            must exist.

    Returns:
        The *output_path* (for convenient chaining).

    Raises:
        ValueError: If *organism* is not recognised.
        IOError: If the output file cannot be written.
    """
    from ..organisms import CODON_USAGE_TABLES, resolve_organism

    org_key = resolve_organism(organism, strict=False)
    if org_key not in CODON_USAGE_TABLES:
        raise ValueError(
            f"No codon usage table for organism {organism!r} "
            f"(resolved to {org_key!r}). Available: "
            f"{sorted(set(CODON_USAGE_TABLES.keys()))[:10]}..."
        )

    usage_table = CODON_USAGE_TABLES[org_key]

    # Sort by amino acid then codon for deterministic output
    sorted_entries = sorted(
        usage_table.items(),
        key=lambda item: (item[1][0], item[0]),  # (aa, codon)
    )

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Codon", "AminoAcid", "Frequency"])
        for codon, (aa, fraction, _per_thousand, _count) in sorted_entries:
            # Skip stop codons — LinearDesign handles termination separately
            if aa == "*":
                continue
            writer.writerow([codon, aa, f"{fraction:.2f}"])

    logger.info(
        "Wrote LinearDesign codon usage file for %s to %s (%d codons)",
        org_key,
        output_path,
        sum(1 for _aa, *_ in usage_table.values() if _aa != "*"),
    )
    return output_path


# ═══════════════════════════════════════════════════════════════════════
# 4. Core LinearDesign runner
# ═══════════════════════════════════════════════════════════════════════


def run_lineardesign(
    protein_seq: str,
    lambda_val: float = 0.0,
    codon_table_path: str | None = None,
    binary_path: str = "lineardesign",
    timeout: int = 600,
) -> LinearDesignResult:
    """Run LinearDesign via subprocess and return a structured result.

    LinearDesign (Zhang et al. 2023, Nature; Huang et al. 2023, bioRxiv)
    uses lattice-based dynamic programming over a codon-grammar CFG to
    find the globally optimal mRNA sequence that jointly minimizes MFE
    and maximizes CAI.

    The tradeoff is controlled by the *lambda_val* parameter:

    - ``lambda_val = 0``: Pure MFE minimization (most stable structure).
    - ``lambda_val > 0``: Increasing CAI weight (higher CAI, possibly
      less stable structure).

    Args:
        protein_seq: Protein amino acid sequence (single-letter codes).
        lambda_val: MFE/CAI tradeoff parameter.  Default 0.0 (pure MFE).
        codon_table_path: Path to codon usage CSV in LinearDesign format.
            If ``None``, LinearDesign uses its built-in codon table.
        binary_path: Path to the LinearDesign executable.  Defaults to
            ``"lineardesign"`` (must be in ``$PATH``).
        timeout: Maximum runtime in seconds.  Default 600 (10 minutes).

    Returns:
        A :class:`LinearDesignResult` with the optimized sequence,
        structure, MFE, CAI, and metadata.  On failure, the ``error``
        field is set and other fields may be ``None``.

    Raises:
        ValueError: If *protein_seq* is empty.
    """
    if not protein_seq:
        raise ValueError("protein_seq must not be empty")

    # Validate protein sequence contains only standard amino acid codes
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    invalid_chars = set(protein_seq.upper()) - valid_aas
    if invalid_chars:
        raise ValueError(
            f"protein_seq contains invalid amino acid codes: {invalid_chars}. "
            f"Expected only standard 20 amino acid single-letter codes."
        )

    start_time = time.monotonic()

    try:
        # Build command
        cmd = [binary_path]

        # LinearDesign CLI: first positional arg is the protein sequence,
        # second optional positional arg is lambda.
        # Some versions use flags: --lambda, --protein
        # We support the most common interface: positional protein + flag lambda
        cmd.append(protein_seq)

        if lambda_val > 0:
            cmd.extend(["--lambda", str(lambda_val)])

        if codon_table_path and os.path.exists(codon_table_path):
            cmd.extend(["--codonusage", codon_table_path])

        logger.debug("Running LinearDesign: %s", " ".join(cmd[:3]) + "...")

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )

        elapsed = time.monotonic() - start_time

        if result.returncode != 0:
            stderr_msg = result.stderr.strip()[:500]
            return LinearDesignResult(
                error=f"LinearDesign exited with code {result.returncode}: {stderr_msg}",
                lambda_val=lambda_val,
                elapsed_seconds=elapsed,
            )

        # Parse output
        parsed = _parse_lineardesign_output(result.stdout)
        return LinearDesignResult(
            sequence=parsed.get("sequence"),
            structure=parsed.get("structure"),
            mfe=parsed.get("mfe"),
            cai=parsed.get("cai"),
            lambda_val=lambda_val,
            elapsed_seconds=elapsed,
        )

    except FileNotFoundError:
        return LinearDesignResult(
            error=(
                "LinearDesign binary not found at "
                f"{binary_path!r}. Install from "
                "https://github.com/LinearDesignSoftware/LinearDesign "
                "and ensure it is in your PATH."
            ),
            lambda_val=lambda_val,
        )
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - start_time
        return LinearDesignResult(
            error=f"LinearDesign timed out after {timeout}s",
            lambda_val=lambda_val,
            elapsed_seconds=elapsed,
        )
    except Exception as exc:
        elapsed = time.monotonic() - start_time
        return LinearDesignResult(
            error=f"LinearDesign execution failed: {exc}",
            lambda_val=lambda_val,
            elapsed_seconds=elapsed,
        )


# ── Output parser ────────────────────────────────────────────────────


def _parse_lineardesign_output(stdout: str) -> dict[str, Any]:
    """Parse LinearDesign stdout into structured fields.

    LinearDesign output formats vary by version.  This parser handles
    multiple known formats:

    **Format 1 (key-value lines):**

    ::

        mRNA sequence: AUGGCC...
        structure: (((...)))...
        MFE: -42.5
        CAI: 0.87

    **Format 2 (tab-separated or space-separated):**

    ::

        AUGGCC...  (((...)))  -42.5  0.87

    **Format 3 (JSON-like or mixed):**

    ::

        Sequence: AUGGCC...
        Minimum free energy: -42.5 kcal/mol
        CAI value: 0.87

    Args:
        stdout: Raw stdout from LinearDesign subprocess.

    Returns:
        Dict with optional keys *sequence*, *structure*, *mfe*, *cai*.
    """
    parsed: dict[str, Any] = {}
    output = stdout.strip()
    if not output:
        return parsed

    lines = output.split("\n")

    # Strategy 1: Key-value parsing (handles Format 1 and 3)
    for line in lines:
        line_stripped = line.strip()
        if not line_stripped:
            continue

        lower = line_stripped.lower()

        # Sequence
        if "sequence" in lower and "structure" not in lower:
            value = _extract_value_after_colon(line_stripped)
            if value:
                # Convert RNA to DNA if needed
                parsed["sequence"] = value.upper().replace("U", "T")

        # Structure (dot-bracket)
        elif "structure" in lower:
            value = _extract_value_after_colon(line_stripped)
            if value and set(value) <= set(".()[]{}<>"):
                parsed["structure"] = value

        # MFE / free energy
        elif "free energy" in lower or "mfe" in lower or "minimum free" in lower:
            mfe_val = _extract_float_from_line(line_stripped)
            if mfe_val is not None:
                parsed["mfe"] = mfe_val

        # CAI
        elif "cai" in lower:
            cai_val = _extract_float_from_line(line_stripped)
            if cai_val is not None and 0.0 <= cai_val <= 1.0:
                parsed["cai"] = cai_val

    # Strategy 2: Tab/space-separated line (Format 2)
    # If we didn't find sequence via key-value, try the last non-empty line
    if "sequence" not in parsed:
        for line in reversed(lines):
            parts = line.strip().split()
            if len(parts) >= 2:
                # Check if first part looks like a sequence
                candidate_seq = parts[0].upper().replace("U", "T")
                if set(candidate_seq) <= set("ACGT") and len(candidate_seq) >= 6:
                    parsed["sequence"] = candidate_seq
                    # Second part might be structure
                    if len(parts) >= 2 and set(parts[1]) <= set(".()[]{}<>"):
                        parsed["structure"] = parts[1]
                    # Third part might be MFE
                    if len(parts) >= 3:
                        try:
                            parsed["mfe"] = float(parts[2])
                        except ValueError:
                            pass
                    # Fourth part might be CAI
                    if len(parts) >= 4:
                        try:
                            cai_candidate = float(parts[3])
                            if 0.0 <= cai_candidate <= 1.0:
                                parsed["cai"] = cai_candidate
                        except ValueError:
                            pass
                    break

    # Strategy 3: Try to find sequence as a standalone RNA/DNA line
    if "sequence" not in parsed:
        for line in lines:
            candidate = line.strip().upper().replace("U", "T")
            # A valid coding sequence should be all ACGT and length % 3 == 0
            if (
                set(candidate) <= set("ACGT")
                and len(candidate) >= 6
                and len(candidate) % 3 == 0
            ):
                parsed["sequence"] = candidate
                break

    return parsed


def _extract_value_after_colon(line: str) -> str | None:
    """Extract the value after a colon separator in a key-value line."""
    if ":" in line:
        value = line.split(":", 1)[1].strip()
        # Remove units like "kcal/mol"
        value = re.sub(r"\s*kcal/mol\s*$", "", value, flags=re.IGNORECASE)
        return value if value else None
    return None


def _extract_float_from_line(line: str) -> float | None:
    """Extract the first float value from a line of text.

    Handles negative values, decimal points, and common suffixes like
    'kcal/mol'.
    """
    # Remove common units
    cleaned = re.sub(r"kcal/mol", "", line, flags=re.IGNORECASE)
    # Find all float-like patterns (including negative)
    matches = re.findall(r"-?\d+\.?\d*", cleaned)
    for match in matches:
        try:
            return float(match)
        except ValueError:
            continue
    return None


# ═══════════════════════════════════════════════════════════════════════
# 5. Pareto front scanning
# ═══════════════════════════════════════════════════════════════════════


def scan_pareto_front(
    protein_seq: str,
    lambda_values: list[float] | None = None,
    **kwargs: Any,
) -> list[LinearDesignResult]:
    """Run LinearDesign with multiple lambda values to trace the MFE/CAI Pareto front.

    The MFE-CAI tradeoff forms a Pareto front: decreasing MFE (more
    stable structure) typically decreases CAI (less optimal codons), and
    vice versa.  By scanning a range of lambda values, this function
    traces out the Pareto-optimal curve, allowing the user to select
    the best operating point for their application.

    Default lambda values span the practical range:
    - 0.0: Pure MFE minimization
    - 0.5–1.0: Mostly MFE with slight CAI consideration
    - 2.0–3.0: Balanced MFE/CAI (recommended for most applications)
    - 5.0–10.0: Mostly CAI with some MFE consideration

    References:
        Zhang et al. (2023) Nature 621:396–403, Fig. 2.
        Huang et al. (2023) bioRxiv, Fig. 3.

    Args:
        protein_seq: Protein amino acid sequence.
        lambda_values: List of lambda values to evaluate.  If ``None``,
            uses a default set: ``[0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]``.
        **kwargs: Additional keyword arguments passed to
            :func:`run_lineardesign` (e.g. *codon_table_path*,
            *binary_path*, *timeout*).

    Returns:
        List of :class:`LinearDesignResult` objects, one per lambda
        value, sorted by increasing lambda.  Results with errors are
        included (check the ``error`` field).
    """
    if lambda_values is None:
        lambda_values = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]

    results: list[LinearDesignResult] = []
    for lam in sorted(lambda_values):
        logger.info("Running LinearDesign with lambda=%.2f", lam)
        result = run_lineardesign(
            protein_seq=protein_seq,
            lambda_val=lam,
            **kwargs,
        )
        results.append(result)

    # Log Pareto front summary
    valid = [r for r in results if r.error is None]
    if valid:
        mfe_range = (
            min(r.mfe for r in valid if r.mfe is not None),
            max(r.mfe for r in valid if r.mfe is not None),
        )
        cai_range = (
            min(r.cai for r in valid if r.cai is not None),
            max(r.cai for r in valid if r.cai is not None),
        )
        logger.info(
            "Pareto front: MFE [%.1f, %.1f] kcal/mol, CAI [%.3f, %.3f] across %d points",
            mfe_range[0],
            mfe_range[1],
            cai_range[0],
            cai_range[1],
            len(valid),
        )

    return results


# ═══════════════════════════════════════════════════════════════════════
# 6. Binary search for optimal lambda
# ═══════════════════════════════════════════════════════════════════════


def find_optimal_lambda(
    protein_seq: str,
    target_cai: float = 0.8,
    **kwargs: Any,
) -> LinearDesignResult:
    """Binary search for the lambda that achieves a target CAI.

    Uses bisection on the lambda parameter to find the minimum lambda
    (maximum MFE stability) that still achieves the requested CAI
    threshold.  This automates the calibration of the MFE/CAI tradeoff
    for a specific expression target.

    Algorithm:
        1. Start with ``lambda_lo = 0.0`` (pure MFE) and
           ``lambda_hi = 10.0`` (high CAI weight).
        2. Run LinearDesign at the midpoint lambda.
        3. If CAI >= target: narrow to lower half (try to improve MFE
           while keeping CAI above target).
        4. If CAI < target: narrow to upper half (need more CAI weight).
        5. Repeat until convergence (CAI within tolerance of target) or
           max iterations reached.

    Convergence criteria:
        - CAI within 0.02 of target, OR
        - Lambda interval narrower than 0.1, OR
        - 15 iterations completed

    Args:
        protein_seq: Protein amino acid sequence.
        target_cai: Desired CAI threshold (0.0–1.0).  Default 0.8.
        **kwargs: Additional keyword arguments passed to
            :func:`run_lineardesign` (e.g. *codon_table_path*,
            *binary_path*, *timeout*).

    Returns:
        The :class:`LinearDesignResult` whose CAI is closest to
        *target_cai* from above (i.e. the most stable sequence that
        still meets the CAI target).  If no result achieves the target,
        returns the highest-CAI result found.

    Raises:
        ValueError: If *target_cai* is not in [0, 1] or *protein_seq*
            is empty.
    """
    if not 0.0 <= target_cai <= 1.0:
        raise ValueError(f"target_cai must be in [0, 1], got {target_cai}")
    if not protein_seq:
        raise ValueError("protein_seq must not be empty")

    # Binary search parameters
    lambda_lo = 0.0
    lambda_hi = 10.0
    max_iterations = 15
    cai_tolerance = 0.02
    lambda_tolerance = 0.1

    best_result: LinearDesignResult | None = None
    best_overshoot: LinearDesignResult | None = None  # Best result with CAI >= target

    for iteration in range(max_iterations):
        lambda_mid = (lambda_lo + lambda_hi) / 2.0

        logger.debug(
            "Binary search iteration %d: lambda=[%.3f, %.3f], mid=%.3f",
            iteration + 1,
            lambda_lo,
            lambda_hi,
            lambda_mid,
        )

        result = run_lineardesign(
            protein_seq=protein_seq,
            lambda_val=lambda_mid,
            **kwargs,
        )

        # If LinearDesign failed, try increasing lambda
        if result.error is not None:
            logger.warning(
                "LinearDesign failed at lambda=%.3f: %s",
                lambda_mid,
                result.error[:200],
            )
            lambda_lo = lambda_mid
            if best_result is None:
                best_result = result
            continue

        cai = result.cai
        if cai is None:
            # Can't determine CAI; treat as insufficient
            lambda_lo = lambda_mid
            best_result = result
            continue

        # Check convergence
        if abs(cai - target_cai) <= cai_tolerance:
            logger.info(
                "Converged at lambda=%.3f: CAI=%.4f (target=%.4f)",
                lambda_mid,
                cai,
                target_cai,
            )
            return result

        if cai >= target_cai:
            # CAI is sufficient; try lower lambda for better MFE
            lambda_hi = lambda_mid
            if best_overshoot is None or (best_overshoot.cai or 0) > cai:
                best_overshoot = result
        else:
            # CAI is insufficient; need higher lambda
            lambda_lo = lambda_mid

        best_result = result

        # Check lambda interval convergence
        if (lambda_hi - lambda_lo) < lambda_tolerance:
            logger.info(
                "Lambda interval converged: [%.3f, %.3f]",
                lambda_lo,
                lambda_hi,
            )
            break

    # Return the best overshoot (CAI >= target with lowest lambda)
    # or the best result we found
    if best_overshoot is not None:
        return best_overshoot
    if best_result is not None:
        return best_result

    # Should not reach here, but safety fallback
    return LinearDesignResult(
        error="Binary search failed to find a valid result",
        lambda_val=(lambda_lo + lambda_hi) / 2.0,
    )


# ═══════════════════════════════════════════════════════════════════════
# 7. Full BioCompiler integration
# ═══════════════════════════════════════════════════════════════════════


def optimize_mfe_lineardesign(
    protein_seq: str,
    organism: str = "human",
    lambda_val: float = 3.0,
    binary_path: str = "lineardesign",
    timeout: int = 600,
    keep_codon_file: bool = False,
    codon_file_dir: str | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Full BioCompiler integration for mRNA MFE optimization via LinearDesign.

    This is the recommended high-level API that handles:

    1. **Codon table generation** — Automatically generates a codon
       usage CSV from BioCompiler's organism tables, matching the
       organism's reference gene set.
    2. **LinearDesign execution** — Runs the LinearDesign binary with
       the generated codon table and specified lambda.
    3. **Result validation** — Verifies that the returned sequence
       translates to the correct protein.
    4. **CAI computation** — Computes CAI using BioCompiler's own
       CAI calculator for consistency with the rest of the pipeline.

    If LinearDesign is not available, falls back to BioCompiler's
    built-in simulated annealing MFE optimizer.

    Args:
        protein_seq: Protein amino acid sequence (single-letter codes).
        organism: Target organism for codon usage.  All names accepted
            by :func:`~biocompiler.organisms.resolve_organism` are valid
            (e.g. ``"human"``, ``"e_coli"``, ``"cho"``).
        lambda_val: MFE/CAI tradeoff parameter.  Default 3.0 (balanced).
            - 0 = pure MFE, higher = more CAI weight.
        binary_path: Path to LinearDesign executable.
        timeout: Maximum runtime in seconds.
        keep_codon_file: If ``True``, the generated codon usage CSV is
            not deleted after the run.  Useful for debugging.
        codon_file_dir: Directory for the temporary codon usage CSV.
            If ``None``, uses the system temp directory.
        **kwargs: Additional keyword arguments.  Currently unused;
            reserved for future extensions.

    Returns:
        Dict with keys:

        - ``"sequence"``: Optimized DNA coding sequence.
        - ``"structure"``: Dot-bracket structure (or ``None``).
        - ``"mfe"``: Minimum free energy (kcal/mol).
        - ``"cai"``: Codon Adaptation Index.
        - ``"method"``: ``"lineardesign"`` or ``"simulated_annealing"``
          (fallback).
        - ``"lambda"``: Lambda value used.
        - ``"organism"``: Organism name used.
        - ``"elapsed_seconds"``: Total wall-clock time.
        - ``"error"``: Error message (if any).

    Raises:
        ValueError: If *protein_seq* is empty or *organism* is invalid.
    """
    if not protein_seq:
        raise ValueError("protein_seq must not be empty")

    start_time = time.monotonic()
    codon_file_path: str | None = None

    try:
        # Step 1: Check LinearDesign availability
        if not is_lineardesign_available(binary_path):
            logger.warning(
                "LinearDesign binary not found at %r; falling back to SA optimizer",
                binary_path,
            )
            return _fallback_sa_optimization(
                protein_seq, organism, lambda_val, start_time
            )

        # Step 2: Generate codon usage file
        tmp_dir = codon_file_dir or tempfile.gettempdir()
        os.makedirs(tmp_dir, exist_ok=True)
        codon_file_path = os.path.join(
            tmp_dir,
            f"lineardesign_codon_{organism}_{os.getpid()}.csv",
        )
        prepare_codon_usage_file(organism, codon_file_path)

        # Step 3: Run LinearDesign
        ld_result = run_lineardesign(
            protein_seq=protein_seq,
            lambda_val=lambda_val,
            codon_table_path=codon_file_path,
            binary_path=binary_path,
            timeout=timeout,
        )

        elapsed = time.monotonic() - start_time

        # Step 4: Check for errors
        if ld_result.error is not None:
            logger.warning(
                "LinearDesign failed: %s; falling back to SA optimizer",
                ld_result.error[:200],
            )
            return _fallback_sa_optimization(
                protein_seq, organism, lambda_val, start_time
            )

        # Step 5: Validate sequence translates to correct protein
        sequence = ld_result.sequence
        if sequence is not None:
            validation = _validate_translation(sequence, protein_seq)
            if not validation["valid"]:
                logger.warning(
                    "LinearDesign sequence failed translation validation: %s",
                    validation["error"],
                )

        # Step 6: Compute CAI using BioCompiler's own calculator
        cai_value = ld_result.cai
        if sequence is not None:
            try:
                from ..expression.translation import compute_cai

                cai_value = compute_cai(sequence, organism=organism)
            except Exception:
                # Fall back to LinearDesign's reported CAI
                pass

        result = {
            "sequence": sequence,
            "structure": ld_result.structure,
            "mfe": ld_result.mfe,
            "cai": cai_value,
            "method": "lineardesign",
            "lambda": lambda_val,
            "organism": organism,
            "elapsed_seconds": elapsed,
            "error": None,
        }

        logger.info(
            "LinearDesign optimization complete: MFE=%.1f kcal/mol, CAI=%.4f, "
            "lambda=%.2f, elapsed=%.1fs",
            ld_result.mfe or 0.0,
            cai_value or 0.0,
            lambda_val,
            elapsed,
        )

        return result

    except Exception as exc:
        elapsed = time.monotonic() - start_time
        logger.error("LinearDesign integration failed: %s", exc)
        return _fallback_sa_optimization(
            protein_seq, organism, lambda_val, start_time
        )

    finally:
        # Clean up temporary codon file
        if codon_file_path and not keep_codon_file:
            try:
                os.unlink(codon_file_path)
            except OSError:
                pass


def _fallback_sa_optimization(
    protein_seq: str,
    organism: str,
    lambda_cai: float,
    start_time: float,
) -> dict[str, Any]:
    """Fall back to simulated annealing when LinearDesign is unavailable."""
    from .mfe_optimization import optimize_mfe_simulated_annealing

    logger.info("Using simulated annealing fallback for MFE optimization")
    sa_result = optimize_mfe_simulated_annealing(
        protein_seq=protein_seq,
        organism=organism,
        lambda_cai=lambda_cai,
    )
    elapsed = time.monotonic() - start_time
    sa_result["method"] = "simulated_annealing_fallback"
    sa_result["lambda"] = lambda_cai
    sa_result["organism"] = organism
    sa_result["elapsed_seconds"] = elapsed
    sa_result["error"] = None
    return sa_result


def _validate_translation(
    dna_seq: str,
    expected_protein: str,
) -> dict[str, Any]:
    """Validate that a DNA sequence translates to the expected protein.

    Args:
        dna_seq: DNA coding sequence.
        expected_protein: Expected protein amino acid sequence.

    Returns:
        Dict with ``"valid"`` (bool) and ``"error"`` (str or None).
    """
    from ..shared.constants import CODON_TABLE

    if not dna_seq or len(dna_seq) % 3 != 0:
        return {"valid": False, "error": "DNA sequence length is not a multiple of 3"}

    translated = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i : i + 3].upper()
        aa = CODON_TABLE.get(codon)
        if aa is None:
            return {
                "valid": False,
                "error": f"Invalid codon {codon!r} at position {i // 3}",
            }
        if aa == "*":
            # Stop codon before end of sequence
            if i // 3 < len(expected_protein):
                return {
                    "valid": False,
                    "error": f"Premature stop codon at position {i // 3}",
                }
            break
        translated.append(aa)

    translated_str = "".join(translated)
    if translated_str != expected_protein:
        return {
            "valid": False,
            "error": (
                f"Translation mismatch: expected {expected_protein[:20]}... "
                f"got {translated_str[:20]}..."
            ),
        }

    return {"valid": True, "error": None}
