"""
BioCompiler CLI
=======================
Command-line interface for certified gene optimization and protein analysis.

Commands:
  optimize        Read FASTA, run full multi-step optimization, write optimized FASTA + certificate
  check           Read FASTA, evaluate all 8 predicates, print certificate
  benchmark       Run built-in benchmarks (eGFP, mCherry, LacZ)
  scan            Scan a DNA sequence for features
  serve           Start the REST API server
  structure       Predict and assess protein structure
  stability       Analyze protein stability
  solubility      Analyze protein solubility
  immunogenicity  Analyze and reduce immunogenicity
  assess          Full protein assessment
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from typing import List, Optional

from . import __version__
from .optimization import BioOptimizer
from .type_system import (
    CODON_TABLE,
    check_no_stop_codons,
    check_no_cryptic_splice,
    check_no_cpg_island,
    check_no_restriction_site,
    check_no_avoidable_gt,
    check_valid_coding_seq,
)
from .certificate import format_certificate, compute_certificate

# Lazy imports for clear_cache functions — imported inside cmd_optimize
# to avoid circular import issues:
#   from .foldx import clear_cache as foldx_clear_cache
#   from .camsol import clear_cache as camsol_clear_cache
#   from .immunogenicity import clear_cache as immunogenicity_clear_cache

logger = logging.getLogger(__name__)


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

def cmd_optimize(args: argparse.Namespace) -> None:
    """Handle the 'optimize' command."""
    # Clear engine caches for a fresh optimization run
    try:
        from .foldx import clear_cache as foldx_clear_cache
        foldx_clear_cache()
    except ImportError:
        pass
    try:
        from .camsol import clear_cache as camsol_clear_cache
        camsol_clear_cache()
    except ImportError:
        pass
    try:
        from .immunogenicity import clear_cache as immunogenicity_clear_cache
        immunogenicity_clear_cache()
    except ImportError:
        pass

    seq = _read_fasta(args.input)

    if len(seq) < 3:
        print("Error: Sequence too short for optimization.", file=sys.stderr)
        sys.exit(1)

    enzymes: List[str] = []
    if args.enzymes:
        enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]

    opt = BioOptimizer(
        species=args.species,
        enzymes=enzymes,
        splice_low=args.splice_low,
        splice_high=args.splice_high,
        avoid_gt=args.avoid_gt,
    )

    optimized, pred_results, cert_text = opt.optimize(seq)

    # Determine output paths
    input_base = os.path.splitext(args.input)[0]
    out_fasta = args.output if args.output else f"{input_base}_optimized.fasta"
    out_cert = args.certificate if args.certificate else f"{input_base}_certificate.txt"

    _write_fasta(out_fasta, optimized, header=f"optimized|{args.species}")
    _write_certificate(out_cert, cert_text)

    cert_level = compute_certificate(pred_results)
    print(f"Optimization complete.")
    print(f"  Input:      {args.input} ({len(seq)} bp)")
    print(f"  Output:     {out_fasta} ({len(optimized)} bp)")
    print(f"  Certificate: {out_cert} ({cert_level.value})")
    print()
    print(cert_text)


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
    results = []
    results.append(check_no_stop_codons(seq))
    results.append(check_no_cryptic_splice(seq, args.splice_low, args.splice_high))
    results.append(check_no_cpg_island(seq))
    results.append(check_no_restriction_site(seq, enzymes))
    results.append(check_no_avoidable_gt(seq))
    results.append(check_valid_coding_seq(seq))

    # ConservationScore — no optimization was performed, so all AA are self-conserved
    from .type_system import PredicateResult
    results.append(PredicateResult(
        "ConservationScore", True,
        details="No substitutions (check-only mode); identity conservation by default"
    ))

    # CodonOptimality
    from .organisms import SPECIES
    species_cai = SPECIES.get(args.species, SPECIES["ecoli"])
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
    from .benchmark import run_benchmark, compare_tools

    enzymes: List[str] = []
    if args.enzymes:
        enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]

    run_benchmark(
        enzymes=enzymes if enzymes else None,
        splice_low=args.splice_low,
        splice_high=args.splice_high,
    )

    compare_tools()


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


def _print_mutation_table(mut_results: list) -> None:
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

    # ── optimize ──
    opt_parser = subparsers.add_parser(
        "optimize",
        help="Optimize a FASTA gene sequence with multi-step optimization pipeline",
    )
    opt_parser.add_argument(
        "input",
        help="Input FASTA file path",
    )
    opt_parser.add_argument(
        "--species", default="human",
        choices=["human", "ecoli"],
        help="Target species for codon optimization (default: human)",
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
        help="Run built-in benchmarks (eGFP, mCherry, LacZ)",
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

    return parser


def main(argv: Optional[List[str]] = None) -> None:
    """Main entry point for the BioCompiler CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == "optimize":
        cmd_optimize(args)
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
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
