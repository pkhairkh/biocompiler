"""
BioCompiler CLI v10.0.0
=======================
Command-line interface for certified gene optimization and protein analysis.

Commands:
  optimize            Optimize a protein sequence for a target organism
  batch               Batch-optimize proteins from a file
  check               Read FASTA, evaluate all 8 predicates, print certificate
  benchmark           Run built-in benchmarks (eGFP, mCherry, LacZ) or named gene sets
  scan                Scan a DNA sequence for features
  serve               Start the REST API server
  structure           Predict and assess protein structure
  stability           Analyze protein stability
  solubility          Analyze protein solubility
  immunogenicity      Analyze and reduce immunogenicity
  assess              Full protein assessment
  validate-cai        Validate CAI against published ground-truth values
  validate-maxentscan Validate MaxEntScan scores against published values
  whatif              Run what-if analysis on a protein sequence
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import random
import sys
import time
import uuid
from typing import List, Optional

__all__ = [
    "main",
    "verify",
    "build_parser",
    "colorize",
    "cmd_optimize",
    "cmd_batch",
    "cmd_check",
    "cmd_benchmark",
    "cmd_scan",
    "cmd_structure",
    "cmd_stability",
    "cmd_solubility",
    "cmd_immunogenicity",
    "cmd_assess",
    "cmd_validate_cai",
    "cmd_validate_maxentscan",
    "cmd_whatif",
    "cmd_explain",
    "cmd_report",
]

logger = logging.getLogger(__name__)

from . import __version__
from .optimization import BioOptimizer
from .type_system import (
    check_no_stop_codons,
    check_no_cryptic_splice,
    check_no_cpg_island,
    check_no_restriction_site,
    check_no_avoidable_gt,
    check_valid_coding_seq,
)
from .certificate import format_certificate, compute_certificate
from .scanner import gc_content

# Lazy imports for clear_cache functions — imported inside cmd_optimize
# to avoid circular import issues:
#   from .foldx import clear_cache as foldx_clear_cache
#   from .camsol import clear_cache as camsol_clear_cache
#   from .immunogenicity import clear_cache as immunogenicity_clear_cache

# ── ANSI colour helpers ──────────────────────────────────────────────────────

_ANSI_RESET = "\033[0m"
_ANSI_BOLD = "\033[1m"
_ANSI_RED = "\033[31m"
_ANSI_GREEN = "\033[32m"
_ANSI_YELLOW = "\033[33m"
_ANSI_CYAN = "\033[36m"
_ANSI_BOLD_RED = "\033[1;31m"
_ANSI_BOLD_GREEN = "\033[1;32m"
_ANSI_BOLD_CYAN = "\033[1;36m"
_ANSI_DIM = "\033[2m"


def _supports_color() -> bool:
    """Return True if stdout is a TTY that likely supports ANSI colours."""
    return hasattr(sys.stdout, "isatty") and sys.stdout.isatty()


def colorize(text: str, *codes: str) -> str:
    """Wrap *text* in ANSI escape codes; no-op when stdout is not a TTY."""
    if not _supports_color():
        return text
    return "".join(codes) + text + _ANSI_RESET


def _section_header(text: str) -> str:
    return colorize(text, _ANSI_BOLD_CYAN)


def _verdict_symbol(value: str) -> str:
    """Return a coloured verdict label."""
    v = value.upper()
    if v in ("PASS", "LIKELY_PASS"):
        return colorize(v, _ANSI_BOLD_GREEN)
    if v in ("FAIL", "LIKELY_FAIL"):
        return colorize(v, _ANSI_BOLD_RED)
    # UNCERTAIN or anything else
    return colorize(v, _ANSI_YELLOW)


def _error_msg(text: str) -> str:
    return colorize(text, _ANSI_RED)


def _success_msg(text: str) -> str:
    return colorize(text, _ANSI_BOLD_GREEN)


def _dim(text: str) -> str:
    return colorize(text, _ANSI_DIM)


def _summary_box(label: str, value: str) -> str:
    """Build a Unicode box around a label + value pair."""
    inner = f"{label}: {value}"
    width = len(inner) + 2  # 1 space padding each side
    top = "\u250c" + "\u2500" * width + "\u2510"
    mid = "\u2502 " + inner + " \u2502"
    bot = "\u2514" + "\u2500" * width + "\u2518"
    return "\n".join([top, mid, bot])


# ── Progress helper ──────────────────────────────────────────────────────────

class _ProgressStep:
    """Context manager that prints a step label to stderr and appends timing on exit."""

    def __init__(self, label: str, verbose: bool = False) -> None:
        self.label = label
        self.verbose = verbose
        self._t0: float = 0.0

    def __enter__(self) -> "_ProgressStep":
        self._t0 = time.perf_counter()
        sys.stderr.write(f"{self.label}...")
        sys.stderr.flush()
        return self

    def __exit__(self, *exc: object) -> None:
        elapsed = time.perf_counter() - self._t0
        timing = f" ({elapsed:.3f}s)" if self.verbose else ""
        sys.stderr.write(f" done{timing}\n")
        sys.stderr.flush()


# ── File I/O helpers ─────────────────────────────────────────────────────────

def _read_fasta(path: str) -> str:
    """Read a FASTA file and return the DNA sequence (uppercase, no whitespace)."""
    if not os.path.isfile(path):
        print(f"Error: File not found: {path}", file=sys.stderr)
        sys.exit(1)
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            seq_parts.append(line.upper())
    seq = "".join(seq_parts)
    # Remove any non-DNA characters
    seq = "".join(c for c in seq if c in "ACGT")
    return seq


def _write_fasta(path: str, seq: str, header: str = "optimized") -> None:
    """Write a DNA sequence to a FASTA file with 80-char line wrapping."""
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")


def _write_certificate(path: str, cert_text: str) -> None:
    """Write certificate text to a file."""
    with open(path, "w") as f:
        f.write(cert_text)


# ── Protein resolution ───────────────────────────────────────────────────────

def _resolve_protein(args: argparse.Namespace) -> str:
    """Return a protein sequence from *args*.

    Accepts either ``--protein`` (1-letter AA string) or ``--sequence``
    (DNA that will be translated).  Exits with an error if neither or both
    are provided.
    """
    protein: str | None = getattr(args, "protein", None)
    sequence: str | None = getattr(args, "sequence", None)

    if protein and sequence:
        print(_error_msg("Error: provide --protein OR --sequence, not both."), file=sys.stderr)
        sys.exit(1)
    if not protein and not sequence:
        print(_error_msg("Error: provide --protein or --sequence."), file=sys.stderr)
        sys.exit(1)

    if protein:
        # Validate amino-acid characters
        valid_aa = set("ACDEFGHIKLMNPQRSTVWYX*")
        cleaned = protein.upper().strip()
        invalid = set(cleaned) - valid_aa
        if invalid:
            print(_error_msg(f"Error: invalid amino-acid characters: {', '.join(sorted(invalid))}"),
                  file=sys.stderr)
            sys.exit(1)
        return cleaned

    # Translate DNA → protein
    from .translation import translate
    dna = sequence.upper().strip()
    dna = "".join(c for c in dna if c in "ACGT")
    if len(dna) < 3:
        print(_error_msg("Error: DNA sequence too short to translate."), file=sys.stderr)
        sys.exit(1)
    protein_seq = translate(dna)
    if not protein_seq:
        print(_error_msg("Error: DNA sequence could not be translated."), file=sys.stderr)
        sys.exit(1)
    print(_dim(f"Translated DNA ({len(dna)} bp) → protein ({len(protein_seq)} aa)"))
    return protein_seq


def _get_organism(args: argparse.Namespace) -> str:
    """Return the organism name from *args*, defaulting to *Homo_sapiens*."""
    return getattr(args, "organism", None) or "Homo_sapiens"


# ── Command: optimize ────────────────────────────────────────────────────────

def _resolve_organism_arg(args: argparse.Namespace) -> str:
    """Resolve the organism from --organism or --species (alias).

    --organism takes precedence; --species is kept as a backward-compatible
    alias.  The value is normalised through ``resolve_organism`` so that
    shorthand names like 'ecoli' or 'human' are accepted.
    """
    from .organisms import resolve_organism as _resolve
    raw = getattr(args, "organism", None) or getattr(args, "species", None) or "Homo_sapiens"
    try:
        return _resolve(raw)
    except Exception:
        return raw


def _resolve_source_organism_arg(args: argparse.Namespace) -> str | None:
    """Resolve the source organism from --source-organism.

    The value is normalised through ``resolve_organism`` so that
    shorthand names like 'ecoli' or 'human' are accepted.
    Returns None if --source-organism is not specified.
    """
    raw = getattr(args, "source_organism", None)
    if raw is None:
        return None
    from .organisms import resolve_organism as _resolve
    try:
        return _resolve(raw)
    except Exception:
        return raw


def cmd_optimize(args: argparse.Namespace) -> None:
    """Handle the 'optimize' command — v10.0.0 unified interface.

    Accepts either a protein sequence as a positional argument or a FASTA
    file via ``--input``.  The ``--organism`` flag (with ``--species`` as an
    alias) specifies the target organism.  New v10 flags:

    * ``--strategy``  — choose optimisation backend (hybrid, constraint_first, csp)
    * ``--no-splice-check`` — skip eukaryotic splice-site constraints (for prokaryotes)
    * ``--codon-pair-bias`` — optimise codon-pair bias during the run
    * ``--json`` — machine-readable JSON output
    * ``--verbose`` — detailed optimisation trace with timing
    """
    verbose: bool = getattr(args, "verbose", False)

    # Set deterministic seed if provided
    seed: Optional[int] = getattr(args, "seed", None)
    if seed is not None:
        random.seed(seed)
        logger.info("Random seed set to %d for reproducible optimization", seed)

    # Resolve organism (supports --organism and --species alias)
    organism = _resolve_organism_arg(args)
    no_splice_check = getattr(args, "no_splice_check", False)

    # Resolve immunogenicity predicate parameters
    source_organism = _resolve_source_organism_arg(args)
    therapeutic = getattr(args, "therapeutic", False)

    # Clear engine caches for a fresh optimization run
    try:
        from .foldx import clear_cache as foldx_clear_cache
        foldx_clear_cache()
    except ImportError:
        logger.debug("foldx module not available; skipping cache clear")
    try:
        from .camsol import clear_cache as camsol_clear_cache
        camsol_clear_cache()
    except ImportError:
        logger.debug("camsol module not available; skipping cache clear")
    try:
        from .immunogenicity import clear_cache as immunogenicity_clear_cache
        immunogenicity_clear_cache()
    except ImportError:
        logger.debug("immunogenicity module not available; skipping cache clear")

    # ── Obtain input sequence ──────────────────────────────────────────────
    # v10: accept either a positional PROTEIN arg or --input FASTA file
    protein_seq: str | None = getattr(args, "protein", None)
    input_fasta: str | None = getattr(args, "input", None)

    if protein_seq and input_fasta:
        print(_error_msg("Error: provide PROTEIN positional arg OR --input, not both."), file=sys.stderr)
        sys.exit(1)

    if protein_seq:
        # Validate amino-acid characters
        valid_aa = set("ACDEFGHIKLMNPQRSTVWYX*")
        cleaned = protein_seq.upper().strip()
        invalid = set(cleaned) - valid_aa
        if invalid:
            print(_error_msg(f"Error: invalid amino-acid characters: {', '.join(sorted(invalid))}"),
                  file=sys.stderr)
            sys.exit(1)
        # Translate protein → DNA via optimisation path
        protein = cleaned
        seq = None  # will be back-translated by optimizer
    elif input_fasta:
        seq = _read_fasta(input_fasta)
        if len(seq) < 3:
            print("Error: Sequence too short for optimization.", file=sys.stderr)
            sys.exit(1)
        protein = None
    else:
        print(_error_msg("Error: provide PROTEIN positional argument or --input FASTA file."), file=sys.stderr)
        sys.exit(1)

    enzymes: List[str] = []
    if getattr(args, "enzymes", None):
        enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]

    # GC bounds
    gc_lo = getattr(args, "gc_lo", 0.30)
    gc_hi = getattr(args, "gc_hi", 0.70)

    # ── Strategy selection (v10) ───────────────────────────────────────────
    strategy = getattr(args, "strategy", "hybrid") or "hybrid"

    # Determine organism domain
    organism_domain_raw = getattr(args, "organism_domain", "auto") or "auto"
    # If --no-splice-check was set, force prokaryote domain
    if no_splice_check and organism_domain_raw == "auto":
        organism_domain_raw = "prokaryote"

    from .api import resolve_organism_domain
    resolved_domain = resolve_organism_domain(organism, organism_domain_raw)

    # ── Choose optimiser backend ───────────────────────────────────────────
    use_codon_pair_bias = getattr(args, "codon_pair_bias", False)

    if strategy == "csp":
        # CSP solver backend (OR-Tools / Z3)
        from .solver.dispatch import csp_optimize, is_solver_available
        if not is_solver_available():
            print(_error_msg("Error: CSP solver not available. Install ortools or z3-solver."),
                  file=sys.stderr)
            sys.exit(1)

        if verbose:
            print(_dim(f"  Strategy: CSP solver"))
            print(_dim(f"  Organism: {organism} (domain: {resolved_domain})"))
            if no_splice_check:
                print(_dim("  Splice check: DISABLED (--no-splice-check)"))
            if use_codon_pair_bias:
                print(_dim("  Codon-pair bias: ENABLED"))

        with _ProgressStep("Optimizing (CSP solver)", verbose=verbose):
            opt_result = csp_optimize(
                protein if protein else "",
                organism=organism,
                gc_lo=gc_lo,
                gc_hi=gc_hi,
            )
        optimized = opt_result.sequence if hasattr(opt_result, "sequence") else str(opt_result)
        pred_results = []
        cert_text = ""

    elif strategy == "constraint_first":
        # Constraint-first: resolve constraints before maximising CAI
        if verbose:
            print(_dim(f"  Strategy: constraint-first"))
            print(_dim(f"  Organism: {organism} (domain: {resolved_domain})"))
            if no_splice_check:
                print(_dim("  Splice check: DISABLED (--no-splice-check)"))
            if use_codon_pair_bias:
                print(_dim("  Codon-pair bias: ENABLED"))

        if protein:
            # Use optimize_sequence with protein input
            from .optimization import optimize_sequence
            with _ProgressStep("Optimizing (constraint-first)", verbose=verbose):
                opt_result = optimize_sequence(
                    protein, organism=organism,
                    gc_lo=gc_lo, gc_hi=gc_hi,
                    consider_codon_pair_bias=use_codon_pair_bias,
                    source_organism=source_organism,
                    therapeutic=therapeutic,
                )
            optimized = opt_result.sequence
            pred_results = opt_result.predicate_results
            cert_text = opt_result.certificate_text
        else:
            # Legacy FASTA-input path
            splice_low = getattr(args, "splice_low", 3.0)
            splice_high = getattr(args, "splice_high", 6.0)
            avoid_gt = getattr(args, "avoid_gt", True)

            opt = BioOptimizer(
                species=organism,
                enzymes=enzymes,
                splice_low=splice_low,
                splice_high=splice_high,
                avoid_gt=avoid_gt,
                organism_domain=resolved_domain,
            )
            optimized, pred_results, cert_text = opt.optimize(seq)
    else:
        # Default: hybrid (greedy init + CAI hill climbing)
        if verbose:
            print(_dim(f"  Strategy: hybrid (default)"))
            print(_dim(f"  Organism: {organism} (domain: {resolved_domain})"))
            if no_splice_check:
                print(_dim("  Splice check: DISABLED (--no-splice-check)"))
            if use_codon_pair_bias:
                print(_dim("  Codon-pair bias: ENABLED"))

        if protein:
            from .optimization import optimize_sequence
            with _ProgressStep("Optimizing (hybrid)", verbose=verbose):
                opt_result = optimize_sequence(
                    protein, organism=organism,
                    gc_lo=gc_lo, gc_hi=gc_hi,
                    consider_codon_pair_bias=use_codon_pair_bias,
                    source_organism=source_organism,
                    therapeutic=therapeutic,
                )
            optimized = opt_result.sequence
            pred_results = opt_result.predicate_results
            cert_text = opt_result.certificate_text
        else:
            # Legacy FASTA-input path
            splice_low = getattr(args, "splice_low", 3.0)
            splice_high = getattr(args, "splice_high", 6.0)
            avoid_gt = getattr(args, "avoid_gt", True)

            opt = BioOptimizer(
                species=organism,
                enzymes=enzymes,
                splice_low=splice_low,
                splice_high=splice_high,
                avoid_gt=avoid_gt,
                organism_domain=resolved_domain,
            )
            optimized, pred_results, cert_text = opt.optimize(seq)

    # ── Codon-pair bias scoring (v10) ──────────────────────────────────────
    cpb_score: float | None = None
    if use_codon_pair_bias:
        try:
            from .codon_pair_scoring import compute_cpb
            cpb_score = compute_cpb(optimized)
        except (ImportError, Exception):
            cpb_score = None

    # ── JSON output (v10) ──────────────────────────────────────────────────
    output_json = getattr(args, "json", False)

    if output_json:
        result_dict: dict = {
            "version": __version__,
            "organism": organism,
            "strategy": strategy,
            "gc_content": round(gc_content(optimized), 4) if optimized else None,
            "sequence_length": len(optimized) if optimized else 0,
            "sequence": optimized,
            "no_splice_check": no_splice_check,
            "codon_pair_bias": cpb_score,
            "source_organism": source_organism,
            "therapeutic": therapeutic,
        }
        if pred_results:
            result_dict["predicate_results"] = [
                {"name": p.predicate, "passed": p.passed, "details": p.details}
                for p in pred_results
            ]
        if cert_text:
            result_dict["certificate_text"] = cert_text
        print(json.dumps(result_dict, indent=2))
        return

    # ── Text output ────────────────────────────────────────────────────────
    if input_fasta:
        # Legacy FASTA mode: write output files
        input_base = os.path.splitext(input_fasta)[0]
        out_fasta = args.output if args.output else f"{input_base}_optimized.fasta"
        out_cert = args.certificate if args.certificate else f"{input_base}_certificate.txt"

        _write_fasta(out_fasta, optimized, header=f"optimized|{organism}")
        _write_certificate(out_cert, cert_text)

        cert_level = compute_certificate(pred_results) if pred_results else None
        print(f"Optimization complete.")
        print(f"  Input:      {input_fasta} ({len(seq)} bp)")
        print(f"  Output:     {out_fasta} ({len(optimized)} bp)")
        if cert_level:
            print(f"  Certificate: {out_cert} ({cert_level.value})")
        if cpb_score is not None:
            print(f"  Codon-pair bias: {cpb_score:.4f}")
        print()
        if cert_text:
            print(cert_text)

        # Provenance tracking
        if getattr(args, "provenance", False):
            from .provenance import ProvenanceTracker, DecisionRecord, OptimizationRecord
            from datetime import datetime, timezone as _tz
            from . import __version__ as _ver

            tracker = ProvenanceTracker(seed=seed or 0)
            input_base_pv = os.path.splitext(input_fasta)[0]
            out_provenance = f"{input_base_pv}_provenance.json"

            opt_record = OptimizationRecord(
                input_sequence=seq,
                output_sequence=optimized,
                organism=organism,
                constraints_applied=[p.name for p in pred_results] if pred_results else [],
                mutations_made=[],
                solver_backend=strategy,
                solve_time=0.0,
                seed_used=seed,
                timestamp=datetime.now(_tz.utc).isoformat(),
                biocompiler_version=_ver,
            )
            tracker.add_optimization_record(opt_record)

            with open(out_provenance, "w") as f:
                f.write(tracker.to_json())

            print(f"  Provenance:  {out_provenance}")
    else:
        # v10 protein mode: print directly
        print()
        print(_section_header("═" * 60))
        print(_section_header("  Optimization Result"))
        print(_section_header("═" * 60))
        print(f"  Organism       : {organism}")
        print(f"  Strategy       : {strategy}")
        print(f"  Protein length : {len(protein)} aa")
        print(f"  Sequence length: {len(optimized)} bp")
        gc_val = gc_content(optimized) if optimized else 0.0
        print(f"  GC content     : {gc_val:.4f}")
        if no_splice_check:
            print(f"  Splice check   : DISABLED")
        if source_organism:
            print(f"  Source organism : {source_organism}")
        if therapeutic:
            print(f"  Therapeutic    : YES")
        if cpb_score is not None:
            print(f"  Codon-pair bias: {cpb_score:.4f}")

        if pred_results:
            cert_level = compute_certificate(pred_results)
            print()
            print(_summary_box("Certificate", _verdict_symbol(cert_level.value)))

        if cert_text:
            print()
            print(cert_text)

        # Print optimized sequence
        print()
        print(_section_header("  Optimized Sequence"))
        for i in range(0, len(optimized), 80):
            print(f"  {optimized[i:i+80]}")

        # Save to file if --output given
        out_path = getattr(args, "output", None)
        if out_path:
            _write_fasta(out_path, optimized, header=f"optimized|{organism}")
            print(_success_msg(f"  Saved to: {out_path}"))


# ── Command: batch ─────────────────────────────────────────────────────────

def cmd_batch(args: argparse.Namespace) -> None:
    """Handle the 'batch' command — optimize multiple proteins from a file.

    v10.0.0: The input file should contain one protein per line (optionally
    with a header line starting with ``#``).  Each protein is optimized
    independently for the given ``--organism`` and results are collected.
    """
    proteins_file = args.proteins_file
    organism = _resolve_organism_arg(args)
    verbose = getattr(args, "verbose", False)
    output_json = getattr(args, "json", False)
    strategy = getattr(args, "strategy", "hybrid") or "hybrid"
    gc_lo = getattr(args, "gc_lo", 0.30)
    gc_hi = getattr(args, "gc_hi", 0.70)
    no_splice_check = getattr(args, "no_splice_check", False)
    use_codon_pair_bias = getattr(args, "codon_pair_bias", False)
    source_organism = _resolve_source_organism_arg(args)
    therapeutic = getattr(args, "therapeutic", False)

    if not os.path.isfile(proteins_file):
        print(_error_msg(f"Error: File not found: {proteins_file}"), file=sys.stderr)
        sys.exit(1)

    # Read proteins from file
    proteins: List[tuple] = []  # (name, sequence)
    with open(proteins_file, "r") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Support "name<whitespace>SEQUENCE" or just "SEQUENCE"
            parts = line.split(None, 1)
            if len(parts) == 2 and all(c in "ACDEFGHIKLMNPQRSTVWYX*" for c in parts[1].upper()):
                proteins.append((parts[0], parts[1].upper()))
            else:
                proteins.append((f"protein_{lineno}", line.upper()))

    if not proteins:
        print(_error_msg("Error: No valid protein sequences found in file."), file=sys.stderr)
        sys.exit(1)

    if verbose:
        print(_dim(f"  Batch mode: {len(proteins)} proteins"))
        print(_dim(f"  Organism: {organism}"))
        print(_dim(f"  Strategy: {strategy}"))

    from .optimization import optimize_sequence

    results: List[dict] = []
    for name, prot_seq in proteins:
        with _ProgressStep(f"Optimizing {name}", verbose=verbose):
            try:
                opt_result = optimize_sequence(
                    prot_seq, organism=organism,
                    gc_lo=gc_lo, gc_hi=gc_hi,
                    consider_codon_pair_bias=use_codon_pair_bias,
                    source_organism=source_organism,
                    therapeutic=therapeutic,
                )
                results.append({
                    "name": name,
                    "status": "ok",
                    "sequence": opt_result.sequence,
                    "gc_content": opt_result.gc_content,
                    "cai": opt_result.cai,
                    "codon_pair_bias": opt_result.codon_pair_bias if use_codon_pair_bias else None,
                })
            except Exception as exc:
                results.append({
                    "name": name,
                    "status": "error",
                    "error": str(exc),
                })

    # ── JSON output ────────────────────────────────────────────────────────
    if output_json:
        output_dict = {
            "version": __version__,
            "organism": organism,
            "strategy": strategy,
            "no_splice_check": no_splice_check,
            "total_proteins": len(proteins),
            "source_organism": source_organism,
            "therapeutic": therapeutic,
            "results": results,
        }
        print(json.dumps(output_dict, indent=2))
        return

    # ── Text output ────────────────────────────────────────────────────────
    print()
    print(_section_header("═" * 60))
    print(_section_header("  Batch Optimization Results"))
    print(_section_header("═" * 60))
    print(f"  Organism : {organism}")
    print(f"  Strategy : {strategy}")
    print(f"  Total    : {len(proteins)} proteins")
    if source_organism:
        print(f"  Source organism : {source_organism}")
    if therapeutic:
        print(f"  Therapeutic    : YES")
    print()

    for r in results:
        if r["status"] == "ok":
            status = _success_msg("OK")
            print(f"  {r['name']:<20s} {status}  GC={r['gc_content']:.4f}  CAI={r['cai']:.4f}")
        else:
            status = _error_msg("FAIL")
            print(f"  {r['name']:<20s} {status}  {r.get('error', 'unknown error')}")

    # Save output file if requested
    out_path = getattr(args, "output", None)
    if out_path:
        with open(out_path, "w") as f:
            for r in results:
                if r["status"] == "ok":
                    f.write(f">{r['name']}|{organism}\n")
                    seq = r["sequence"]
                    for i in range(0, len(seq), 80):
                        f.write(seq[i:i+80] + "\n")
        print(_success_msg(f"  Saved to: {out_path}"))


# ── Command: check ───────────────────────────────────────────────────────────

def cmd_check(args: argparse.Namespace) -> None:
    """Handle the 'check' command — evaluate all 8 predicates without optimizing."""
    seq = _read_fasta(args.input)

    if len(seq) < 3:
        print("Error: Sequence too short for checking.", file=sys.stderr)
        sys.exit(1)

    enzymes: List[str] = []
    if args.enzymes:
        enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]

    # Evaluate all 8 predicates
    from .type_system import PredicateResult
    results: List[PredicateResult] = []
    results.append(check_no_stop_codons(seq))
    results.append(check_no_cryptic_splice(seq, args.splice_low, args.splice_high))
    results.append(check_no_cpg_island(seq))
    results.append(check_no_restriction_site(seq, enzymes))
    results.append(check_no_avoidable_gt(seq))
    results.append(check_valid_coding_seq(seq))

    # ConservationScore — no optimization was performed, so all AA are self-conserved
    results.append(PredicateResult(
        "ConservationScore", True,
        details="No substitutions (check-only mode); identity conservation by default"
    ))

    # CodonOptimality
    from .organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism
    _canonical = resolve_organism(args.species)
    species_cai = dict(CODON_ADAPTIVENESS_TABLES.get(
        _canonical, CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
    ))
    import math
    log_sum = 0.0
    count = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        cai = species_cai.get(codon, 0.0)
        if cai <= 0:
            cai = 0.001
        log_sum += math.log(cai)
        count += 1
    overall_cai = math.exp(log_sum / count) if count > 0 else 0.0
    results.append(PredicateResult(
        "CodonOptimality", True,
        details=f"Overall CAI = {overall_cai:.4f} for species '{args.species}'"
    ))

    cert_text = format_certificate(results, seq, args.species)
    cert_level = compute_certificate(results)

    print(f"Certificate: {cert_level.value}")
    print()
    print(cert_text)


# ── Command: benchmark ───────────────────────────────────────────────────────

def cmd_benchmark(args: argparse.Namespace) -> None:
    """Handle the 'benchmark' command."""
    # Set deterministic seed if provided
    seed: Optional[int] = getattr(args, "seed", None)
    if seed is not None:
        random.seed(seed)
        logger.info("Random seed set to %d for reproducible benchmark", seed)

    from .benchmark import (
        run_benchmark, compare_tools,
        REFERENCE_GENES, GENE_PANEL,
        run_structured_benchmarks,
        format_benchmark_report_text, format_benchmark_report_json,
    )

    # ── List gene sets ──
    if getattr(args, "list_gene_sets", False):
        print()
        print(_section_header("Available Gene Sets"))
        print(f"  DEFAULT           Built-in eGFP, mCherry, LacZ (both human & ecoli)")
        print(f"  REFERENCE_GENES   {', '.join(REFERENCE_GENES.keys())}")
        print(f"  HUMAN_THERAPEUTIC {', '.join(k for k, v in GENE_PANEL.items() if v[1] == 'Homo_sapiens')}")
        print(f"  GENE_PANEL        {', '.join(GENE_PANEL.keys())}")
        print()
        return

    enzymes: List[str] = []
    if args.enzymes:
        enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]

    gene_set = getattr(args, "gene_set", None)
    output_file = getattr(args, "output", None)

    if gene_set:
        # Run structured benchmark for a named gene set
        gene_set_upper = gene_set.upper()
        if gene_set_upper == "REFERENCE_GENES":
            report = run_structured_benchmarks(gene_names=list(REFERENCE_GENES.keys()))
            text = format_benchmark_report_text(report)
            print(text)
        elif gene_set_upper in ("HUMAN_THERAPEUTIC", "GENE_PANEL"):
            # Run comprehensive benchmark on GENE_PANEL
            from .benchmark import run_multi_gene_comparison
            results = run_multi_gene_comparison(
                enzymes=enzymes if enzymes else None,
            )
            # Print summary table
            print()
            print(_section_header("═" * 80))
            print(_section_header(f"  Benchmark: {gene_set}"))
            print(_section_header("═" * 80))
            for r in results:
                tool = r.get("tool", "?")
                gene = r.get("gene", "?")
                cai = r.get("cai", 0.0)
                gc = r.get("gc_content", 0.0)
                violations = r.get("constraint_violations", "?")
                success = r.get("success", False)
                status = _success_msg("OK") if success else _error_msg("FAIL")
                print(f"  {tool:<20s} {gene:<12s} CAI={cai:.4f} GC={gc:.4f} violations={violations} {status}")
        else:
            print(_error_msg(f"Unknown gene set: {gene_set}"), file=sys.stderr)
            print(_dim("  Use --list-gene-sets to see available gene sets."), file=sys.stderr)
            sys.exit(1)
    else:
        # Default: run classic benchmark
        run_benchmark(
            enzymes=enzymes if enzymes else None,
            splice_low=args.splice_low,
            splice_high=args.splice_high,
        )

        compare_tools()

    # Save output if requested
    if output_file:
        # v10: detect output format from file extension (.json → JSON, else CSV)
        if output_file.lower().endswith(".json"):
            try:
                from .benchmark import run_structured_benchmarks, format_benchmark_report_json
                report = run_structured_benchmarks()
                json_data = format_benchmark_report_json(report)
                with open(output_file, "w") as f:
                    f.write(json_data)
                print(_success_msg(f"Benchmark results saved to {output_file}"), file=sys.stderr)
            except Exception as exc:
                print(_error_msg(f"Error saving benchmark results: {exc}"), file=sys.stderr)
        else:
            try:
                from .benchmark import run_structured_benchmarks
                report = run_structured_benchmarks()
                import csv as csv_mod
                with open(output_file, "w", newline="") as csvfile:
                    writer = csv_mod.writer(csvfile)
                    writer.writerow(["gene", "test", "passed", "expected", "actual", "time_ms"])
                    for r in report.results:
                        writer.writerow([
                            r.gene_name, r.test_name, r.passed,
                            r.expected, r.actual, r.execution_time_ms,
                        ])
                print(_success_msg(f"Benchmark results saved to {output_file}"), file=sys.stderr)
            except Exception as exc:
                print(_error_msg(f"Error saving benchmark results: {exc}"), file=sys.stderr)


# ── Command: scan ────────────────────────────────────────────────────────────

def cmd_scan(args: argparse.Namespace) -> None:
    """Handle the 'scan' command — scan a sequence for features."""
    from .scanner import scan_sequence
    seq = args.sequence.upper().strip()
    seq = "".join(c for c in seq if c in "ACGT")
    if len(seq) < 3:
        print("Error: Sequence too short for scanning.", file=sys.stderr)
        sys.exit(1)
    enzymes: List[str] = []
    if args.enzymes:
        enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]
    tokens = scan_sequence(seq, restriction_enzymes=enzymes)
    print(f"Scanned {len(seq)} bp sequence")
    print(f"Tokens found: {len(tokens)}")
    for t in tokens:
        print(f"  {t.element_type} at pos {t.position}: {t.match_sequence} (score={t.score:.2f})")


# ── Command: structure ───────────────────────────────────────────────────────

def cmd_structure(args: argparse.Namespace) -> None:
    """Predict and assess protein structure."""
    organism = _get_organism(args)

    # Quality-only mode: assess an existing PDB file
    if getattr(args, "quality_only", False):
        pdb_file = getattr(args, "pdb_file", None)
        if not pdb_file:
            print(_error_msg("Error: --quality-only requires --pdb-file."), file=sys.stderr)
            sys.exit(1)
        from .structure.quality import compute_structure_quality
        # Read PDB file content — compute_structure_quality expects a PDB string, not a path
        if not os.path.isfile(pdb_file):
            print(_error_msg(f"Error: PDB file not found: {pdb_file}"), file=sys.stderr)
            sys.exit(1)
        with open(pdb_file, "r") as f:
            pdb_content = f.read()
        with _ProgressStep("Assessing structure quality", verbose=getattr(args, "verbose", False)):
            report = compute_structure_quality(pdb_content)
        _print_structure_quality(report)
        return

    # Full prediction mode
    protein = _resolve_protein(args)
    from .esmfold import predict_structure, is_esmfold_available
    from .structure.quality import compute_structure_quality

    esmfold_ok = is_esmfold_available()
    if not esmfold_ok:
        print(_dim("ESMFold not available — using offline/fallback prediction."))

    with _ProgressStep("Predicting structure", verbose=getattr(args, "verbose", False)):
        pdb_path = predict_structure(protein, organism=organism)

    with _ProgressStep("Computing quality metrics"):
        report = compute_structure_quality(pdb_path)

    print()
    print(_section_header("═" * 60))
    print(_section_header("  Structure Prediction & Quality Report"))
    print(_section_header("═" * 60))
    print(f"  Protein length : {len(protein)} aa")
    print(f"  Organism       : {organism}")
    print(f"  ESMFold        : {'available' if esmfold_ok else 'offline/fallback'}")
    print()

    _print_structure_quality(report)

    output = getattr(args, "output", None)
    if output and pdb_path:
        import shutil
        shutil.copy2(pdb_path, output)
        print(_success_msg(f"PDB file saved to {output}"))


def _print_structure_quality(report: object) -> None:
    """Pretty-print a StructureQualityReport."""
    print(_section_header("  Quality Metrics"))
    print(f"  pLDDT           : {getattr(report, 'plddt', 'N/A')}")
    print(f"  Ramachandran    : {getattr(report, 'ramachandran_favored', 'N/A')}")
    print(f"  Clash score     : {getattr(report, 'clash_score', 'N/A')}")

    # Verdict
    plddt = getattr(report, "plddt", None)
    if plddt is not None:
        if isinstance(plddt, (int, float)):
            if plddt >= 90:
                verdict = "PASS"
            elif plddt >= 70:
                verdict = "LIKELY_PASS"
            elif plddt >= 50:
                verdict = "UNCERTAIN"
            else:
                verdict = "LIKELY_FAIL"
        else:
            verdict = "UNCERTAIN"
    else:
        verdict = "UNCERTAIN"

    print()
    print(_summary_box("Structure Verdict", _verdict_symbol(verdict)))


# ── Command: stability ───────────────────────────────────────────────────────

def cmd_stability(args: argparse.Namespace) -> None:
    """Analyze protein stability."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)

    from .foldx import empirical_stability, scan_mutations

    with _ProgressStep("Computing stability", verbose=getattr(args, "verbose", False)):
        result = empirical_stability(protein, organism=organism)

    print()
    print(_section_header("═" * 60))
    print(_section_header("  Protein Stability Analysis"))
    print(_section_header("═" * 60))
    print(f"  Protein length : {len(protein)} aa")
    print(f"  Organism       : {organism}")
    print()

    # Core stability metrics
    print(_section_header("  Stability Metrics"))
    dg = getattr(result, "delta_g", None)
    print(f"  \u0394G (kcal/mol)   : {dg if dg is not None else 'N/A'}")
    print(f"  Energy components:")
    for comp_name, comp_val in getattr(result, "energy_components", {}).items():
        print(f"    {comp_name:20s}: {comp_val}")

    # Verdict
    if dg is not None and isinstance(dg, (int, float)):
        if dg < -5.0:
            verdict = "PASS"
        elif dg < 0:
            verdict = "LIKELY_PASS"
        elif dg < 5:
            verdict = "UNCERTAIN"
        else:
            verdict = "LIKELY_FAIL"
    else:
        verdict = "UNCERTAIN"

    print()
    print(_summary_box("Stability Verdict", _verdict_symbol(verdict)))

    # Mutation scanning
    if getattr(args, "scan_mutations", False):
        positions = getattr(args, "positions", None) or list(range(1, len(protein) + 1))
        with _ProgressStep("Scanning mutations", verbose=getattr(args, "verbose", False)):
            mut_results = scan_mutations(protein, positions=positions, organism=organism)
        _print_mutation_table(mut_results)


def _print_mutation_table(mut_results: List[object]) -> None:
    """Print a table of mutation scan results."""
    print()
    print(_section_header("  Mutation Scan Results"))
    header = f"  {'Position':>8s}  {'Original':>8s}  {'Mutant':>8s}  {'\u0394\u0394G (kcal/mol)':>14s}  {'Effect':>12s}"
    print(header)
    print(f"  {'─' * 8}  {'─' * 8}  {'─' * 8}  {'─' * 14}  {'─' * 12}")

    for mut in mut_results:
        pos = getattr(mut, "position", "?")
        orig = getattr(mut, "original_aa", "?")
        mutant = getattr(mut, "mutant_aa", "?")
        ddg = getattr(mut, "ddg", None)
        ddg_str = f"{ddg:+.2f}" if isinstance(ddg, (int, float)) else "N/A"
        effect = getattr(mut, "effect", "neutral")
        if isinstance(ddg, (int, float)):
            if ddg < -1.0:
                effect = "stabilizing"
            elif ddg > 1.0:
                effect = "destabilizing"
            else:
                effect = "neutral"
        effect_colored = (
            _success_msg(effect) if effect == "stabilizing"
            else _error_msg(effect) if effect == "destabilizing"
            else effect
        )
        print(f"  {pos:>8}  {orig:>8}  {mutant:>8}  {ddg_str:>14}  {effect_colored}")


# ── Command: solubility ─────────────────────────────────────────────────────

def cmd_solubility(args: argparse.Namespace) -> None:
    """Analyze protein solubility."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)

    from .camsol import compute_solubility, find_solubility_mutations

    with _ProgressStep("Computing solubility", verbose=getattr(args, "verbose", False)):
        result = compute_solubility(protein, organism=organism)

    print()
    print(_section_header("═" * 60))
    print(_section_header("  Protein Solubility Analysis"))
    print(_section_header("═" * 60))
    print(f"  Protein length : {len(protein)} aa")
    print(f"  Organism       : {organism}")
    print()

    # Core metrics
    print(_section_header("  Solubility Metrics"))
    score = getattr(result, "camsol_score", None)
    print(f"  CamSol score    : {score if score is not None else 'N/A'}")

    agg_regions = getattr(result, "aggregation_prone_regions", [])
    print(f"  Agg-prone regions: {len(agg_regions)}")
    for region in agg_regions:
        start = getattr(region, "start", "?")
        end = getattr(region, "end", "?")
        seq = getattr(region, "sequence", "")
        rscore = getattr(region, "score", "N/A")
        print(f"    [{start}-{end}] {seq}  (score: {rscore})")

    # Verdict
    min_score = getattr(args, "min_score", None)
    if score is not None and isinstance(score, (int, float)):
        if score >= 1.0:
            verdict = "PASS"
        elif score >= 0.0:
            verdict = "LIKELY_PASS"
        elif score >= -1.0:
            verdict = "UNCERTAIN"
        else:
            verdict = "LIKELY_FAIL"
        if min_score is not None and score < min_score:
            verdict = "FAIL"
    else:
        verdict = "UNCERTAIN"

    print()
    print(_summary_box("Solubility Verdict", _verdict_symbol(verdict)))

    # Recommendations
    recommendations = getattr(result, "recommendations", [])
    if recommendations:
        print()
        print(_section_header("  Recommendations"))
        for i, rec in enumerate(recommendations, 1):
            print(f"  {i}. {rec}")

    # Find solubility-improving mutations
    if getattr(args, "find_mutations", False):
        with _ProgressStep("Finding solubility-improving mutations",
                            verbose=getattr(args, "verbose", False)):
            mut_results = find_solubility_mutations(protein, organism=organism)
        print()
        print(_section_header("  Solubility-Improving Mutations"))
        header = f"  {'Position':>8s}  {'Original':>8s}  {'Mutant':>8s}  {'\u0394CamSol':>10s}  {'Effect':>12s}"
        print(header)
        print(f"  {'─' * 8}  {'─' * 8}  {'─' * 8}  {'─' * 10}  {'─' * 12}")
        for mut in mut_results:
            pos = getattr(mut, "position", "?")
            orig = getattr(mut, "original_aa", "?")
            mutant = getattr(mut, "mutant_aa", "?")
            delta = getattr(mut, "delta_camsol", None)
            delta_str = f"{delta:+.3f}" if isinstance(delta, (int, float)) else "N/A"
            print(f"  {pos:>8}  {orig:>8}  {mutant:>8}  {delta_str:>10}  {_success_msg('improving')}")


# ── Command: immunogenicity ─────────────────────────────────────────────────

def cmd_immunogenicity(args: argparse.Namespace) -> None:
    """Analyze and reduce immunogenicity."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)
    source_organism = _resolve_source_organism_arg(args)
    therapeutic = getattr(args, "therapeutic", False)

    from .immunogenicity import compute_immunogenicity
    from .deimmunization import deimmunize

    mhc_alleles = getattr(args, "mhc_alleles", None) or []

    with _ProgressStep("Computing immunogenicity", verbose=getattr(args, "verbose", False)):
        result = compute_immunogenicity(protein, organism=organism, mhc_alleles=mhc_alleles)

    print()
    print(_section_header("═" * 60))
    print(_section_header("  Immunogenicity Analysis"))
    print(_section_header("═" * 60))
    print(f"  Protein length : {len(protein)} aa")
    print(f"  Organism       : {organism}")
    if source_organism:
        print(f"  Source organism : {source_organism}")
    if therapeutic:
        print(f"  Therapeutic    : YES")
    if mhc_alleles:
        print(f"  MHC alleles    : {', '.join(mhc_alleles)}")
    print()

    # Core metrics
    print(_section_header("  Immunogenicity Metrics"))
    score = getattr(result, "immunogenicity_score", None)
    print(f"  Score           : {score if score is not None else 'N/A'}")

    epitopes = getattr(result, "epitopes", [])
    print(f"  Epitopes found  : {len(epitopes)}")
    for i, epi in enumerate(epitopes, 1):
        start = getattr(epi, "start", "?")
        end = getattr(epi, "end", "?")
        seq = getattr(epi, "sequence", "")
        escore = getattr(epi, "score", "N/A")
        allele = getattr(epi, "allele", "N/A")
        print(f"    {i:3d}. [{start}-{end}] {seq}  allele={allele}  score={escore}")

    # Verdict
    target = getattr(args, "target_score", 0.3)
    if score is not None and isinstance(score, (int, float)):
        if score <= target:
            verdict = "PASS"
        elif score <= target + 0.2:
            verdict = "LIKELY_PASS"
        elif score <= target + 0.5:
            verdict = "UNCERTAIN"
        else:
            verdict = "LIKELY_FAIL"
    else:
        verdict = "UNCERTAIN"

    print()
    print(_summary_box("Immunogenicity Verdict", _verdict_symbol(verdict)))

    # Deimmunization
    if getattr(args, "deimmunize", False):
        max_mut = getattr(args, "max_mutations", 10)
        blosum_min = getattr(args, "blosum62_min", 1)

        with _ProgressStep("Running deimmunization optimization",
                            verbose=getattr(args, "verbose", False)):
            deimm_result = deimmunize(
                protein,
                organism=organism,
                target_score=target,
                max_mutations=max_mut,
                blosum62_min=blosum_min,
                mhc_alleles=mhc_alleles,
            )

        print()
        print(_section_header("  Deimmunization Results"))
        new_score = getattr(deimm_result, "immunogenicity_score", None)
        mutations = getattr(deimm_result, "mutations", [])
        new_protein = getattr(deimm_result, "protein", None)

        print(f"  Original score  : {score if score is not None else 'N/A'}")
        print(f"  New score       : {new_score if new_score is not None else 'N/A'}")
        print(f"  Mutations       : {len(mutations)}")

        if mutations:
            print()
            header = f"  {'Position':>8s}  {'Original':>8s}  {'Mutant':>8s}  {'BLOSUM62':>8s}  {'Epitope':>8s}"
            print(header)
            print(f"  {'─' * 8}  {'─' * 8}  {'─' * 8}  {'─' * 8}  {'─' * 8}")
            for mut in mutations:
                pos = getattr(mut, "position", "?")
                orig = getattr(mut, "original_aa", "?")
                mutant = getattr(mut, "mutant_aa", "?")
                blosum = getattr(mut, "blosum62", "?")
                epitope = getattr(mut, "epitope_affected", "?")
                print(f"  {pos:>8}  {orig:>8}  {mutant:>8}  {blosum:>8}  {epitope:>8}")

        if new_protein:
            print()
            print(_dim(f"  Optimized protein: {new_protein[:60]}{'...' if len(new_protein) > 60 else ''}"))

        # Re-assess verdict after deimmunization
        if new_score is not None and isinstance(new_score, (int, float)):
            if new_score <= target:
                new_verdict = "PASS"
            elif new_score <= target + 0.2:
                new_verdict = "LIKELY_PASS"
            else:
                new_verdict = "UNCERTAIN"
            print()
            print(_summary_box("Post-deimmunization Verdict", _verdict_symbol(new_verdict)))


# ── Command: assess ──────────────────────────────────────────────────────────

def cmd_assess(args: argparse.Namespace) -> None:
    """Full protein assessment."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)
    pdb_file = getattr(args, "pdb_file", None)
    output = getattr(args, "output", None)
    fmt = getattr(args, "format", "text") or "text"

    skip_structure = getattr(args, "skip_structure", False)
    skip_stability = getattr(args, "skip_stability", False)
    skip_solubility = getattr(args, "skip_solubility", False)
    skip_immunogenicity = getattr(args, "skip_immunogenicity", False)

    from .structure.report import (
        assess_protein,
        format_assessment_text,
        format_assessment_json,
        format_assessment_html,
    )

    # If a PDB file path was given, read its content; assess_protein expects
    # pdb_string (PDB content), not a file path.
    pdb_content: str | None = None
    if pdb_file:
        if not os.path.isfile(pdb_file):
            print(_error_msg(f"Error: PDB file not found: {pdb_file}"), file=sys.stderr)
            sys.exit(1)
        with open(pdb_file, "r") as f:
            pdb_content = f.read()

    with _ProgressStep("Running comprehensive assessment",
                        verbose=getattr(args, "verbose", False)):
        report = assess_protein(
            protein,
            organism=organism,
            pdb_string=pdb_content,
            run_structure=not skip_structure,
            run_stability=not skip_stability,
            run_solubility=not skip_solubility,
            run_immunogenicity=not skip_immunogenicity,
        )

    # Format output
    if fmt == "json":
        formatted = format_assessment_json(report)
        print(formatted)
        if output:
            with open(output, "w") as f:
                f.write(formatted)
            print(_success_msg(f"JSON report saved to {output}"), file=sys.stderr)
    elif fmt == "html":
        formatted = format_assessment_html(report)
        print(formatted)
        if output:
            with open(output, "w") as f:
                f.write(formatted)
            print(_success_msg(f"HTML report saved to {output}"), file=sys.stderr)
    else:
        formatted = format_assessment_text(report)
        print(formatted)
        if output:
            with open(output, "w") as f:
                f.write(formatted)
            print(_success_msg(f"Text report saved to {output}"), file=sys.stderr)


# ── Command: validate-cai ─────────────────────────────────────────────────

def cmd_validate_cai(args: argparse.Namespace) -> None:
    """Validate CAI computation against published ground-truth values."""
    organism = _get_organism(args)

    try:
        from .validation.ground_truth import GROUND_TRUTH_DATA, validate_against_ground_truth
    except ImportError:
        print(_error_msg("Error: validation.ground_truth module not available."), file=sys.stderr)
        sys.exit(1)

    print()
    print(_section_header("═" * 60))
    print(_section_header("  CAI Validation Against Published Values"))
    print(_section_header("═" * 60))
    print(f"  Organism : {organism}")
    print()

    # Filter entries by organism
    entries = [e for e in GROUND_TRUTH_DATA if e.organism == organism]
    if not entries:
        available = sorted(set(e.organism for e in GROUND_TRUTH_DATA))
        print(_error_msg(f"No ground-truth entries for organism '{organism}'."))
        print(_dim(f"  Available organisms: {', '.join(available)}"))
        sys.exit(1)

    all_pass = True
    for entry in entries:
        result = validate_against_ground_truth(
            optimized_sequence=entry.published_sequence,
            gene_name=entry.gene_name,
            organism=entry.organism,
        )
        status = _success_msg("PASS") if result.matches_expected else _error_msg("FAIL")
        print(f"  {entry.gene_name:<15s} {status}")
        print(f"    CAI diff: {result.cai_difference:.4f}, GC diff: {result.gc_difference:.4f}")
        if not result.matches_expected:
            all_pass = False

    print()
    if all_pass:
        print(_summary_box("CAI Validation", _verdict_symbol("PASS")))
    else:
        print(_summary_box("CAI Validation", _verdict_symbol("FAIL")))


# ── Command: validate-maxentscan ──────────────────────────────────────────

def cmd_validate_maxentscan(args: argparse.Namespace) -> None:
    """Validate MaxEntScan scoring against published splice-site scores."""
    try:
        from .maxentscan import score_donor, score_acceptor, max_donor_score, max_acceptor_score
    except ImportError:
        print(_error_msg("Error: maxentscan module not available."), file=sys.stderr)
        sys.exit(1)

    # Known canonical donor / acceptor sequences with their published
    # MaxEntScan scores (Yeo & Burge 2004, Table 1 representative values).
    # These are well-characterized splice-site 9-mers / 23-mers.
    donor_tests = [
        # (9-mer, position of G in GT, expected score range)
        ("CAGGTAAGT", 3, (3.0, 12.0)),   # Strong canonical donor
        ("AAGGTGAGT", 2, (3.0, 12.0)),   # Canonical donor
        ("TTGGTAAAT", 2, (0.0, 8.0)),    # Weaker donor
    ]

    acceptor_tests = [
        # (23-mer, position of A in AG, expected score range)
        # Typical strong acceptor: polypyrimidine tract + CAG
        ("TTTTTTTTTTTTTTTTTTTTCAGG", 20, (3.0, 14.0)),
        ("CCCCCCCCCCCCCCCCCCCCCAGA", 20, (0.0, 10.0)),
    ]

    print()
    print(_section_header("═" * 60))
    print(_section_header("  MaxEntScan Validation"))
    print(_section_header("═" * 60))
    print()

    all_pass = True

    # Donor validation
    print(_section_header("  Donor Site Scoring"))
    for seq, pos, (lo, hi) in donor_tests:
        score = score_donor(seq, pos)
        in_range = lo <= score <= hi
        status = _success_msg("PASS") if in_range else _error_msg("FAIL")
        print(f"  {seq}  pos={pos}  score={score:.2f}  expected=[{lo},{hi}]  {status}")
        if not in_range:
            all_pass = False

    # Acceptor validation
    print()
    print(_section_header("  Acceptor Site Scoring"))
    for seq, pos, (lo, hi) in acceptor_tests:
        score = score_acceptor(seq, pos)
        in_range = lo <= score <= hi
        status = _success_msg("PASS") if in_range else _error_msg("FAIL")
        print(f"  {seq[:10]}...  pos={pos}  score={score:.2f}  expected=[{lo},{hi}]  {status}")
        if not in_range:
            all_pass = False

    # Edge-case validation: impossible scores for non-sites
    print()
    print(_section_header("  Edge-Case Validation"))
    # Non-GT site at the donor position should return a very low score
    non_donor_score = score_donor("CAGATACGT", 3)
    non_donor_ok = non_donor_score < 3.0
    status = _success_msg("PASS") if non_donor_ok else _error_msg("FAIL")
    print(f"  Non-donor (no GT)  score={non_donor_score:.2f}  expected<3.0  {status}")
    if not non_donor_ok:
        all_pass = False

    print()
    if all_pass:
        print(_summary_box("MaxEntScan Validation", _verdict_symbol("PASS")))
    else:
        print(_summary_box("MaxEntScan Validation", _verdict_symbol("FAIL")))


# ── Command: whatif ───────────────────────────────────────────────────────

def cmd_whatif(args: argparse.Namespace) -> None:
    """Run what-if analysis on a protein sequence."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)

    print()
    print(_section_header("═" * 60))
    print(_section_header("  What-If Analysis"))
    print(_section_header("═" * 60))
    print(f"  Protein length : {len(protein)} aa")
    print(f"  Organism       : {organism}")
    print()

    # 1. Codon optimization what-if
    try:
        from .optimization import optimize_sequence
        with _ProgressStep("Running codon optimization what-if",
                            verbose=getattr(args, "verbose", False)):
            opt_result = optimize_sequence(
                target_protein=protein,
                organism=organism,
            )
        print(_section_header("  Codon Optimization What-If"))
        print(f"  CAI            : {getattr(opt_result, 'cai', 'N/A')}")
        print(f"  GC content     : {getattr(opt_result, 'gc_content', 'N/A')}")
        print(f"  Fallback used  : {getattr(opt_result, 'fallback_used', False)}")
        satisfied = getattr(opt_result, 'satisfied_predicates', [])
        failed = getattr(opt_result, 'failed_predicates', [])
        print(f"  Satisfied      : {len(satisfied)} predicates")
        print(f"  Failed         : {len(failed)} predicates")
        if failed:
            print(f"  Failed list    : {', '.join(failed)}")
        print()
    except ImportError:
        print(_dim("  Codon optimization module not available; skipping."))
        print()

    # 2. Stability what-if
    try:
        from .foldx import empirical_stability
        with _ProgressStep("Running stability what-if",
                            verbose=getattr(args, "verbose", False)):
            stab_result = empirical_stability(protein, organism=organism)
        dg = getattr(stab_result, "delta_g", None)
        print(_section_header("  Stability What-If"))
        print(f"  ΔG (kcal/mol)  : {dg if dg is not None else 'N/A'}")
        print()
    except ImportError:
        print(_dim("  Stability module not available; skipping."))
        print()

    # 3. Solubility what-if
    try:
        from .camsol import compute_solubility
        with _ProgressStep("Running solubility what-if",
                            verbose=getattr(args, "verbose", False)):
            sol_result = compute_solubility(protein, organism=organism)
        camsol_score = getattr(sol_result, "camsol_score", None)
        agg_regions = getattr(sol_result, "aggregation_prone_regions", [])
        print(_section_header("  Solubility What-If"))
        print(f"  CamSol score   : {camsol_score if camsol_score is not None else 'N/A'}")
        print(f"  Agg-prone regions: {len(agg_regions)}")
        print()
    except ImportError:
        print(_dim("  Solubility module not available; skipping."))
        print()

    # 4. Immunogenicity what-if
    try:
        from .immunogenicity import compute_immunogenicity
        with _ProgressStep("Running immunogenicity what-if",
                            verbose=getattr(args, "verbose", False)):
            imm_result = compute_immunogenicity(protein, organism=organism)
        imm_score = getattr(imm_result, "immunogenicity_score", None)
        epitopes = getattr(imm_result, "epitopes", [])
        print(_section_header("  Immunogenicity What-If"))
        print(f"  Score          : {imm_score if imm_score is not None else 'N/A'}")
        print(f"  Epitopes found : {len(epitopes)}")
        print()
    except ImportError:
        print(_dim("  Immunogenicity module not available; skipping."))
        print()

    # 5. Literature validation what-if
    try:
        from .literature_validation import evaluate_case, ALL_LITERATURE_CASES
        # Find matching cases for this protein
        matching = [c for c in ALL_LITERATURE_CASES
                     if c.sequence_type == "protein" and c.sequence.strip().upper() == protein.upper()]
        if matching:
            print(_section_header("  Literature Validation What-If"))
            print(f"  Matching literature cases: {len(matching)}")
            for case in matching:
                result = evaluate_case(case)
                tp = "TP" if result.true_positive else ""
                fn = "FN" if result.false_negative else ""
                fp = "FP" if result.false_positive else ""
                tn = "TN" if result.true_negative else ""
                cls = tp or fn or fp or tn or "?"
                print(f"    {case.case_id}: {case.name[:50]}  [{cls}]")
            print()
        else:
            print(_dim("  No matching literature cases for this protein sequence."))
            print()
    except ImportError:
        print(_dim("  Literature validation module not available; skipping."))
        print()

    print(_summary_box("What-If Complete", _verdict_symbol("PASS")))


# ── Argument parser ──────────────────────────────────────────────────────────

def _add_protein_args(parser: argparse.ArgumentParser) -> None:
    """Add the common --protein / --sequence / --organism arguments to a subparser."""
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--protein", metavar="TEXT",
        help="Protein sequence in 1-letter amino-acid code",
    )
    group.add_argument(
        "--sequence", metavar="TEXT",
        help="DNA sequence (will be translated to protein)",
    )
    parser.add_argument(
        "--organism", metavar="TEXT",
        help="Organism name (default: Homo_sapiens)",
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the BioCompiler CLI."""
    parser = argparse.ArgumentParser(
        prog="biocompiler",
        description=f"BioCompiler v{__version__} — Certified Gene Optimization with Formal Verification",
    )
    parser.add_argument(
        "--version", action="version", version=f"BioCompiler v{__version__}",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # ── optimize ── (v10.0.0: positional PROTEIN arg + --input fallback)
    opt_parser = subparsers.add_parser(
        "optimize",
        help="Optimize a protein sequence for a target organism",
        description="Optimize a protein/DNA sequence for expression in a target organism. "
                    "Accepts a protein sequence as a positional argument (v10) or a FASTA "
                    "file via --input (legacy mode).",
    )
    opt_parser.add_argument(
        "protein",
        nargs="?",
        default=None,
        help="Protein sequence in 1-letter amino-acid code (e.g. MSKGEELFTGV...)",
    )
    opt_parser.add_argument(
        "--input", "-i", default=None,
        help="Input FASTA file path (legacy mode; mutually exclusive with PROTEIN arg)",
    )
    opt_parser.add_argument(
        "--organism", default=None,
        help="Target organism for codon optimization (e.g., ecoli, human, Homo_sapiens)",
    )
    opt_parser.add_argument(
        "--species", default=None, dest="species",
        help="Alias for --organism (backward compatible)",
    )
    opt_parser.add_argument(
        "--strategy",
        choices=["hybrid", "constraint_first", "csp"],
        default="hybrid",
        help="Optimization strategy (default: hybrid)",
    )
    opt_parser.add_argument(
        "--gc-lo", type=float, default=0.30, metavar="FLOAT",
        help="Minimum GC content fraction (default: 0.30)",
    )
    opt_parser.add_argument(
        "--gc-hi", type=float, default=0.70, metavar="FLOAT",
        help="Maximum GC content fraction (default: 0.70)",
    )
    opt_parser.add_argument(
        "--no-splice-check", action="store_true", default=False,
        help="Skip eukaryotic splice-site constraints (recommended for prokaryotes)",
    )
    opt_parser.add_argument(
        "--codon-pair-bias", action="store_true", default=False,
        help="Optimize codon-pair bias during the optimization run",
    )
    opt_parser.add_argument(
        "--json", action="store_true", default=False,
        help="Output results as JSON (for programmatic use)",
    )
    opt_parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show detailed optimization trace with timing",
    )
    opt_parser.add_argument(
        "--enzymes", default="",
        help="Comma-separated restriction enzymes to avoid (e.g., EcoRI,BamHI)",
    )
    opt_parser.add_argument(
        "--splice-low", type=float, default=3.0,
        help="Low splice score threshold (default: 3.0)",
    )
    opt_parser.add_argument(
        "--splice-high", type=float, default=6.0,
        help="High splice score threshold (default: 6.0)",
    )
    opt_parser.add_argument(
        "--avoid-gt", action="store_true", default=True,
        help="Enable GT dinucleotide avoidance (default: True)",
    )
    opt_parser.add_argument(
        "--no-avoid-gt", action="store_false", dest="avoid_gt",
        help="Disable GT dinucleotide avoidance",
    )
    opt_parser.add_argument(
        "--output", "-o", default=None,
        help="Output FASTA file path (default: <input>_optimized.fasta)",
    )
    opt_parser.add_argument(
        "--certificate", "-c", default=None,
        help="Certificate output file path (default: <input>_certificate.txt)",
    )
    opt_parser.add_argument(
        "--seed", type=int, default=None, metavar="INT",
        help="Deterministic seed for reproducible optimization (default: None)",
    )
    opt_parser.add_argument(
        "--provenance", action="store_true", default=False,
        help="Track provenance and save decision trail as JSON",
    )
    opt_parser.add_argument(
        "--organism-domain",
        dest="organism_domain",
        choices=["auto", "eukaryote", "prokaryote"],
        default="auto",
        help=(
            "Organism domain for constraint selection. "
            "'auto' detects from organism name (default), "
            "'eukaryote' forces eukaryotic constraints (splice sites, CpG), "
            "'prokaryote' skips eukaryote-specific constraints."
        ),
    )
    opt_parser.add_argument(
        "--source-organism",
        dest="source_organism",
        default=None,
        metavar="ORGANISM",
        help=(
            "Organism the protein originates from (e.g., ecoli, human, Homo_sapiens). "
            "Used by immunogenicity predicates to determine self-protein status. "
            "If not specified, the protein is assumed to be from the host organism (self). "
            "Accepts the same aliases as --organism."
        ),
    )
    opt_parser.add_argument(
        "--therapeutic",
        action="store_true",
        default=False,
        help=(
            "Mark the protein as intended for therapeutic use. "
            "This enables stricter immunogenicity thresholds because "
            "immune responses can compromise drug efficacy."
        ),
    )

    # ── batch ── (v10.0.0)
    batch_parser = subparsers.add_parser(
        "batch",
        help="Batch-optimize proteins from a file",
        description="Optimize multiple protein sequences from a text file. "
                    "Each line should be a protein sequence (optionally prefixed "
                    "with a name and whitespace). Lines starting with '#' are ignored.",
    )
    batch_parser.add_argument(
        "proteins_file",
        help="File containing protein sequences (one per line)",
    )
    batch_parser.add_argument(
        "--organism", default="ecoli",
        help="Target organism for codon optimization (default: ecoli)",
    )
    batch_parser.add_argument(
        "--species", default=None, dest="species",
        help="Alias for --organism (backward compatible)",
    )
    batch_parser.add_argument(
        "--strategy",
        choices=["hybrid", "constraint_first", "csp"],
        default="hybrid",
        help="Optimization strategy (default: hybrid)",
    )
    batch_parser.add_argument(
        "--gc-lo", type=float, default=0.30, metavar="FLOAT",
        help="Minimum GC content fraction (default: 0.30)",
    )
    batch_parser.add_argument(
        "--gc-hi", type=float, default=0.70, metavar="FLOAT",
        help="Maximum GC content fraction (default: 0.70)",
    )
    batch_parser.add_argument(
        "--no-splice-check", action="store_true", default=False,
        help="Skip eukaryotic splice-site constraints (recommended for prokaryotes)",
    )
    batch_parser.add_argument(
        "--codon-pair-bias", action="store_true", default=False,
        help="Optimize codon-pair bias during the optimization run",
    )
    batch_parser.add_argument(
        "--json", action="store_true", default=False,
        help="Output results as JSON (for programmatic use)",
    )
    batch_parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show detailed optimization trace with timing",
    )
    batch_parser.add_argument(
        "--output", "-o", default=None,
        help="Output multi-FASTA file path for batch results",
    )
    batch_parser.add_argument(
        "--source-organism",
        dest="source_organism",
        default=None,
        metavar="ORGANISM",
        help=(
            "Organism the proteins originate from (e.g., ecoli, human, Homo_sapiens). "
            "Used by immunogenicity predicates to determine self-protein status. "
            "Accepts the same aliases as --organism."
        ),
    )
    batch_parser.add_argument(
        "--therapeutic",
        action="store_true",
        default=False,
        help="Mark proteins as intended for therapeutic use (stricter immunogenicity thresholds).",
    )

    # ── check ──
    check_parser = subparsers.add_parser(
        "check",
        help="Check a FASTA gene sequence against all 8 predicates",
    )
    check_parser.add_argument(
        "input",
        help="Input FASTA file path",
    )
    check_parser.add_argument(
        "--species", default="human",
        choices=["human", "ecoli"],
        help="Target species for codon evaluation (default: human)",
    )
    check_parser.add_argument(
        "--enzymes", default="",
        help="Comma-separated restriction enzymes to check (e.g., EcoRI,BamHI)",
    )
    check_parser.add_argument(
        "--splice-low", type=float, default=3.0,
        help="Low splice score threshold (default: 3.0)",
    )
    check_parser.add_argument(
        "--splice-high", type=float, default=6.0,
        help="High splice score threshold (default: 6.0)",
    )

    # ── benchmark ──
    bench_parser = subparsers.add_parser(
        "benchmark",
        help="Run built-in benchmarks (eGFP, mCherry, LacZ) or named gene sets",
    )
    bench_parser.add_argument(
        "--gene-set", metavar="NAME",
        help="Named gene set to benchmark (e.g., HUMAN_THERAPEUTIC, REFERENCE_GENES, GENE_PANEL)",
    )
    bench_parser.add_argument(
        "--list-gene-sets", action="store_true",
        help="List available gene sets and exit",
    )
    bench_parser.add_argument(
        "--output", metavar="FILE",
        help="Save benchmark results to CSV file",
    )
    bench_parser.add_argument(
        "--enzymes", default="",
        help="Comma-separated restriction enzymes to avoid (default: EcoRI,BamHI,HindIII,XhoI)",
    )
    bench_parser.add_argument(
        "--splice-low", type=float, default=3.0,
        help="Low splice score threshold (default: 3.0)",
    )
    bench_parser.add_argument(
        "--splice-high", type=float, default=6.0,
        help="High splice score threshold (default: 6.0)",
    )
    bench_parser.add_argument(
        "--seed", type=int, default=None, metavar="INT",
        help="Deterministic seed for reproducible benchmarks (default: None)",
    )

    # ── scan ──
    scan_parser = subparsers.add_parser(
        "scan",
        help="Scan a DNA sequence for features (splice sites, restriction sites, etc.)",
    )
    scan_parser.add_argument(
        "--sequence", "-s", required=True,
        help="DNA sequence to scan",
    )
    scan_parser.add_argument(
        "--enzymes", default="",
        help="Comma-separated restriction enzymes to scan for",
    )

    # ── explain ──
    explain_parser = subparsers.add_parser(
        "explain",
        help="Explain why a specific codon was chosen in a previous optimization",
    )
    explain_parser.add_argument(
        "file",
        help="Provenance JSON file from a previous optimization",
    )
    explain_parser.add_argument(
        "--position", type=int, required=True, metavar="INT",
        help="0-based nucleotide position to explain",
    )

    # ── report ──
    report_parser = subparsers.add_parser(
        "report",
        help="Generate provenance report from saved trail",
    )
    report_parser.add_argument(
        "file",
        help="Provenance JSON file from a previous optimization",
    )
    report_parser.add_argument(
        "--format", default="text", choices=["text", "markdown", "json"],
        help="Output format (default: text)",
    )

    # ── serve ──
    serve_parser = subparsers.add_parser(
        "serve",
        help="Start the REST API server",
    )
    serve_parser.add_argument(
        "--host", default="0.0.0.0",
        help="Host to bind to (default: 0.0.0.0)",
    )
    serve_parser.add_argument(
        "--port", type=int, default=8000,
        help="Port to bind to (default: 8000)",
    )

    # ── structure ──
    struct_parser = subparsers.add_parser(
        "structure",
        help="Predict and assess protein structure",
        description="Predict protein 3D structure using ESMFold and assess quality "
                    "(pLDDT, Ramachandran, clash score).",
    )
    _add_protein_args(struct_parser)
    struct_parser.add_argument(
        "--output", metavar="FILE",
        help="Save predicted PDB structure to file",
    )
    struct_parser.add_argument(
        "--quality-only", action="store_true",
        help="Only assess quality of a provided PDB file (no prediction)",
    )
    struct_parser.add_argument(
        "--pdb-file", metavar="FILE",
        help="Input PDB file for quality assessment (used with --quality-only)",
    )
    struct_parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show timing information",
    )

    # ── stability ──
    stab_parser = subparsers.add_parser(
        "stability",
        help="Analyze protein stability",
        description="Analyze protein thermodynamic stability using FoldX-style "
                    "empirical energy functions.",
    )
    _add_protein_args(stab_parser)
    stab_parser.add_argument(
        "--scan-mutations", action="store_true",
        help="Also scan for stabilizing/destabilizing mutations",
    )
    stab_parser.add_argument(
        "--positions", type=int, nargs="+", metavar="INT",
        help="Specific positions to scan (1-based; default: all positions)",
    )
    stab_parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show timing information",
    )

    # ── solubility ──
    sol_parser = subparsers.add_parser(
        "solubility",
        help="Analyze protein solubility",
        description="Analyze protein solubility using the CamSol method. "
                    "Identifies aggregation-prone regions and suggests improvements.",
    )
    _add_protein_args(sol_parser)
    sol_parser.add_argument(
        "--find-mutations", action="store_true",
        help="Suggest solubility-improving mutations",
    )
    sol_parser.add_argument(
        "--min-score", type=float, metavar="FLOAT",
        help="Minimum acceptable solubility score (fails below this)",
    )
    sol_parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show timing information",
    )

    # ── immunogenicity ──
    imm_parser = subparsers.add_parser(
        "immunogenicity",
        help="Analyze and reduce immunogenicity",
        description="Analyze protein immunogenicity via MHC binding prediction "
                    "and optionally run deimmunization optimization.",
    )
    _add_protein_args(imm_parser)
    imm_parser.add_argument(
        "--deimmunize", action="store_true",
        help="Run deimmunization optimization",
    )
    imm_parser.add_argument(
        "--target-score", type=float, default=0.3, metavar="FLOAT",
        help="Target immunogenicity score (default: 0.3)",
    )
    imm_parser.add_argument(
        "--max-mutations", type=int, default=10, metavar="INT",
        help="Maximum mutations for deimmunization (default: 10)",
    )
    imm_parser.add_argument(
        "--blosum62-min", type=int, default=1, metavar="INT",
        help="Minimum BLOSUM62 score for allowed mutations (default: 1)",
    )
    imm_parser.add_argument(
        "--mhc-alleles", nargs="+", metavar="TEXT",
        help="Specific MHC alleles to check (e.g. HLA-A*02:01)",
    )
    imm_parser.add_argument(
        "--source-organism",
        dest="source_organism",
        default=None,
        metavar="ORGANISM",
        help=(
            "Organism the protein originates from (e.g., ecoli, human, Homo_sapiens). "
            "Used by immunogenicity predicates to determine self-protein status. "
            "Accepts the same aliases as --organism."
        ),
    )
    imm_parser.add_argument(
        "--therapeutic",
        action="store_true",
        default=False,
        help="Mark the protein as intended for therapeutic use (stricter immunogenicity thresholds).",
    )
    imm_parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show timing information",
    )

    # ── assess ──
    assess_parser = subparsers.add_parser(
        "assess",
        help="Full protein assessment",
        description="Run all analyses: structure, stability, solubility, "
                    "and immunogenicity in a single comprehensive assessment.",
    )
    _add_protein_args(assess_parser)
    assess_parser.add_argument(
        "--pdb-file", metavar="FILE",
        help="Optional PDB structure file to use instead of prediction",
    )
    assess_parser.add_argument(
        "--skip-structure", action="store_true",
        help="Skip structure prediction/assessment",
    )
    assess_parser.add_argument(
        "--skip-stability", action="store_true",
        help="Skip stability analysis",
    )
    assess_parser.add_argument(
        "--skip-solubility", action="store_true",
        help="Skip solubility analysis",
    )
    assess_parser.add_argument(
        "--skip-immunogenicity", action="store_true",
        help="Skip immunogenicity analysis",
    )
    assess_parser.add_argument(
        "--output", "-o", metavar="FILE",
        help="Save report to file",
    )
    assess_parser.add_argument(
        "--format", choices=["text", "json", "html"], default="text",
        help="Output format (default: text)",
    )
    assess_parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show timing information",
    )

    # ── validate-cai ──
    vcai_parser = subparsers.add_parser(
        "validate-cai",
        help="Validate CAI against published ground-truth values",
        description="Validate BioCompiler's CAI computation by comparing against "
                    "published codon-optimized sequences and their reported CAI values.",
    )
    vcai_parser.add_argument(
        "--organism", metavar="TEXT", default="Escherichia_coli",
        help="Organism to validate (default: Escherichia_coli)",
    )

    # ── validate-maxentscan ──
    vmes_parser = subparsers.add_parser(
        "validate-maxentscan",
        help="Validate MaxEntScan scores against published values",
        description="Validate BioCompiler's MaxEntScan splice-site scoring by "
                    "comparing against known canonical and non-canonical splice sites.",
    )

    # ── whatif ──
    whatif_parser = subparsers.add_parser(
        "whatif",
        help="Run what-if analysis on a protein sequence",
        description="Run a comprehensive what-if analysis: codon optimization, "
                    "stability, solubility, and immunogenicity predictions for a "
                    "given protein sequence.",
    )
    _add_protein_args(whatif_parser)
    whatif_parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show timing information",
    )

    return parser


# ── Command: explain ───────────────────────────────────────────────────────

def cmd_explain(args: argparse.Namespace) -> None:
    """Explain why a specific codon was chosen in a previous optimization."""
    from .provenance import ProvenanceTracker

    if not os.path.isfile(args.file):
        print(_error_msg(f"Error: File not found: {args.file}"), file=sys.stderr)
        sys.exit(1)

    with open(args.file, "r") as f:
        data = json.load(f)

    tracker = ProvenanceTracker.from_dict(data)
    position = args.position

    decisions = tracker.get_decisions_for_position(position)

    print()
    print(_section_header("═" * 60))
    print(_section_header("  Provenance Explanation"))
    print(_section_header("═" * 60))
    print(f"  Position:  {position}")
    print(f"  File:      {args.file}")
    print()

    if not decisions:
        print(_dim(f"  No decisions recorded for position {position}."))
        print(_dim("  This may indicate the position was not modified during optimization"))
        print(_dim("  or provenance tracking was not enabled for that position."))
    else:
        print(_section_header(f"  Decisions at position {position}"))
        print(f"  {'Type':<22s}  {'Chosen':<10s}  {'Alternatives':<30s}")
        print(f"  {'─' * 22}  {'─' * 10}  {'─' * 30}")
        for d in decisions:
            alts = ", ".join(d.alternatives_considered[:5])
            if len(d.alternatives_considered) > 5:
                alts += f" (+{len(d.alternatives_considered) - 5} more)"
            print(f"  {d.decision_type:<22s}  {d.chosen_value:<10s}  {alts}")
        print()
        print(_section_header("  Rationale"))
        for d in decisions:
            print(f"  [{d.decision_type}] {d.rationale}")
            if d.constraint_context:
                ctx_parts = [f"{k}={v}" for k, v in d.constraint_context.items()]
                print(_dim(f"    Context: {', '.join(ctx_parts)}"))

    print()


# ── Command: report ────────────────────────────────────────────────────────

def cmd_report(args: argparse.Namespace) -> None:
    """Generate provenance report from saved trail."""
    from .provenance import ProvenanceTracker, generate_provenance_report

    if not os.path.isfile(args.file):
        print(_error_msg(f"Error: File not found: {args.file}"), file=sys.stderr)
        sys.exit(1)

    with open(args.file, "r") as f:
        data = json.load(f)

    tracker = ProvenanceTracker.from_dict(data)
    records = tracker.get_optimization_records()
    fmt = getattr(args, "format", "text") or "text"

    print()
    print(_section_header("═" * 60))
    print(_section_header("  Provenance Report"))
    print(_section_header("═" * 60))
    print(f"  Source:      {args.file}")
    print(f"  Seed:        {tracker.seed}")
    print(f"  Decisions:   {len(tracker.get_full_audit_trail())}")
    print(f"  Opt records: {len(records)}")
    print()

    if fmt == "json":
        print(json.dumps(data, indent=2, sort_keys=True))
    elif fmt == "markdown":
        # Generate markdown report
        print("# BioCompiler Provenance Report\n")
        print(f"- **Source file**: `{args.file}`")
        print(f"- **Seed**: {tracker.seed}")
        print(f"- **Decisions recorded**: {len(tracker.get_full_audit_trail())}")
        print(f"- **Optimization runs**: {len(records)}\n")
        for idx, rec in enumerate(records, 1):
            print(f"## Run {idx}\n")
            print(f"| Field | Value |")
            print(f"|-------|-------|")
            print(f"| Timestamp | {rec.timestamp} |")
            print(f"| Organism | {rec.organism} |")
            print(f"| Solver | {rec.solver_backend} |")
            print(f"| Seed | {rec.seed_used} |")
            print(f"| Solve time | {rec.solve_time:.3f}s |")
            print(f"| Version | {rec.biocompiler_version} |")
            print(f"| Input length | {len(rec.input_sequence)} |")
            print(f"| Output length | {len(rec.output_sequence)} |")
            print(f"| Mutations | {', '.join(rec.mutations_made) or '(none)'} |")
            print(f"| Constraints | {', '.join(rec.constraints_applied) or '(none)'} |")
            print()
        # Decision trail table
        trail = tracker.get_full_audit_trail()
        if trail:
            print("## Decision Trail\n")
            print("| # | Type | Position | Chosen | Alternatives | Rationale |")
            print("|---|------|----------|--------|-------------|-----------|")
            for i, d in enumerate(trail, 1):
                alts = ", ".join(d.alternatives_considered[:3])
                if len(d.alternatives_considered) > 3:
                    alts += "..."
                # Escape pipes in rationale for markdown
                rationale = d.rationale.replace("|", "\\|")
                print(f"| {i} | {d.decision_type} | {d.position} | {d.chosen_value} | {alts} | {rationale} |")
            print()
    else:
        # Text format
        report_text = generate_provenance_report(records)
        print(report_text)


def main(argv: Optional[List[str]] = None) -> None:
    """Main entry point for the BioCompiler CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == "optimize":
        cmd_optimize(args)
    elif args.command == "batch":
        cmd_batch(args)
    elif args.command == "check":
        cmd_check(args)
    elif args.command == "benchmark":
        cmd_benchmark(args)
    elif args.command == "scan":
        cmd_scan(args)
    elif args.command == "serve":
        import uvicorn
        from .api import app
        uvicorn.run(app, host=args.host, port=args.port)
    elif args.command == "structure":
        cmd_structure(args)
    elif args.command == "stability":
        cmd_stability(args)
    elif args.command == "solubility":
        cmd_solubility(args)
    elif args.command == "immunogenicity":
        cmd_immunogenicity(args)
    elif args.command == "assess":
        cmd_assess(args)
    elif args.command == "validate-cai":
        cmd_validate_cai(args)
    elif args.command == "validate-maxentscan":
        cmd_validate_maxentscan(args)
    elif args.command == "explain":
        cmd_explain(args)
    elif args.command == "report":
        cmd_report(args)
    elif args.command == "whatif":
        cmd_whatif(args)
    else:
        parser.print_help()
        sys.exit(1)


def verify() -> None:
    """Entry point for ``biocompiler-verify`` — runs the check command."""
    main(["check"] + sys.argv[1:])


if __name__ == "__main__":
    main()
