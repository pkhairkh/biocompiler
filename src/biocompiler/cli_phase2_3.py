"""
BioCompiler CLI Phase 2+3 — Protein structure, stability, solubility & immunogenicity
======================================================================================
Subcommands for protein-level analysis: structure prediction, stability assessment,
solubility analysis, immunogenicity profiling, and comprehensive assessment.

Commands:
  structure        Predict and assess protein structure
  stability        Analyze protein stability
  solubility       Analyze protein solubility
  immunogenicity   Analyze and reduce immunogenicity
  assess           Full protein assessment (all Phase 2+3)
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time

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

class _ProgressPhase:
    """Context manager that prints a phase label to stderr and appends timing on exit."""

    def __init__(self, label: str, verbose: bool = False) -> None:
        self.label = label
        self.verbose = verbose
        self._t0: float = 0.0

    def __enter__(self) -> "_ProgressPhase":
        self._t0 = time.perf_counter()
        sys.stderr.write(f"{self.label}...")
        sys.stderr.flush()
        return self

    def __exit__(self, *exc: object) -> None:
        elapsed = time.perf_counter() - self._t0
        timing = f" ({elapsed:.3f}s)" if self.verbose else ""
        sys.stderr.write(f" done{timing}\n")
        sys.stderr.flush()


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
    from biocompiler.translation import translate
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


# ── cmd_structure ────────────────────────────────────────────────────────────

def cmd_structure(args: argparse.Namespace) -> None:
    """Predict and assess protein structure."""
    organism = _get_organism(args)

    # Quality-only mode: assess an existing PDB file
    if getattr(args, "quality_only", False):
        pdb_file = getattr(args, "pdb_file", None)
        if not pdb_file:
            print(_error_msg("Error: --quality-only requires --pdb-file."), file=sys.stderr)
            sys.exit(1)
        from biocompiler.structure_quality import compute_structure_quality
        with _ProgressPhase("Assessing structure quality", verbose=getattr(args, "verbose", False)):
            report = compute_structure_quality(pdb_file)
        _print_structure_quality(report)
        return

    # Full prediction mode
    protein = _resolve_protein(args)
    from biocompiler.esmfold import predict_structure, is_esmfold_available
    from biocompiler.structure_quality import compute_structure_quality

    esmfold_ok = is_esmfold_available()
    if not esmfold_ok:
        print(_dim("ESMFold not available — using offline/fallback prediction."))

    with _ProgressPhase("Predicting structure", verbose=getattr(args, "verbose", False)):
        pdb_path = predict_structure(protein, organism=organism)

    with _ProgressPhase("Computing quality metrics"):
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


# ── cmd_stability ────────────────────────────────────────────────────────────

def cmd_stability(args: argparse.Namespace) -> None:
    """Analyze protein stability."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)

    from biocompiler.foldx import empirical_stability, scan_mutations

    with _ProgressPhase("Computing stability", verbose=getattr(args, "verbose", False)):
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
        with _ProgressPhase("Scanning mutations", verbose=getattr(args, "verbose", False)):
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


# ── cmd_solubility ───────────────────────────────────────────────────────────

def cmd_solubility(args: argparse.Namespace) -> None:
    """Analyze protein solubility."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)

    from biocompiler.camsol import compute_solubility, find_solubility_mutations

    with _ProgressPhase("Computing solubility", verbose=getattr(args, "verbose", False)):
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
        with _ProgressPhase("Finding solubility-improving mutations",
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


# ── cmd_immunogenicity ──────────────────────────────────────────────────────

def cmd_immunogenicity(args: argparse.Namespace) -> None:
    """Analyze and reduce immunogenicity."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)

    from biocompiler.immunogenicity import compute_immunogenicity
    from biocompiler.deimmunization import deimmunize

    mhc_alleles = getattr(args, "mhc_alleles", None) or []

    with _ProgressPhase("Computing immunogenicity", verbose=getattr(args, "verbose", False)):
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

        with _ProgressPhase("Running deimmunization optimization",
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


# ── cmd_assess ───────────────────────────────────────────────────────────────

def cmd_assess(args: argparse.Namespace) -> None:
    """Full protein assessment (all Phase 2+3 analyses)."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)
    pdb_file = getattr(args, "pdb_file", None)
    output = getattr(args, "output", None)
    fmt = getattr(args, "format", "text") or "text"

    skip_structure = getattr(args, "skip_structure", False)
    skip_stability = getattr(args, "skip_stability", False)
    skip_solubility = getattr(args, "skip_solubility", False)
    skip_immunogenicity = getattr(args, "skip_immunogenicity", False)

    from biocompiler.structure_report import (
        assess_protein,
        format_assessment_text,
        format_assessment_json,
        format_assessment_html,
    )

    with _ProgressPhase("Running comprehensive assessment",
                        verbose=getattr(args, "verbose", False)):
        report = assess_protein(
            protein,
            organism=organism,
            pdb_file=pdb_file,
            skip_structure=skip_structure,
            skip_stability=skip_stability,
            skip_solubility=skip_solubility,
            skip_immunogenicity=skip_immunogenicity,
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


# ── Command registration ────────────────────────────────────────────────────

def register_phase2_3_commands(subparsers: argparse._SubParsersAction) -> None:
    """Register all Phase 2+3 subcommands on an existing argparse subparser group."""

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
    struct_parser.set_defaults(func=cmd_structure)

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
    stab_parser.set_defaults(func=cmd_stability)

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
    sol_parser.set_defaults(func=cmd_solubility)

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
    imm_parser.set_defaults(func=cmd_immunogenicity)

    # ── assess ──
    assess_parser = subparsers.add_parser(
        "assess",
        help="Full protein assessment (all Phase 2+3)",
        description="Run all Phase 2+3 analyses: structure, stability, solubility, "
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
    assess_parser.set_defaults(func=cmd_assess)


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
