"""
BioCompiler CLI — Command Handlers
====================================
Thin command handlers that parse args, call service functions, and format output.

Each cmd_* function follows the pattern:
  1. Extract parameters from argparse Namespace
  2. Call a service function from application/
  3. Format and print the result using formatters/

Extracted from cli.py as part of the SoC refactoring (Wave 4b).
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import random
import sys
from typing import List, Optional

from .formatters import (
    ProgressStep,
    colorize,
    dim,
    error_msg,
    section_header,
    success_msg,
    summary_box,
    verdict_symbol,
    read_fasta,
    write_fasta,
    write_certificate,
    print_structure_quality,
    print_mutation_table,
)

logger = logging.getLogger(__name__)


# ── Helpers ──────────────────────────────────────────────────────────────────

def _resolve_protein(args: argparse.Namespace) -> str:
    """Return a protein sequence from *args*.

    Accepts either ``--protein`` (1-letter AA string) or ``--sequence``
    (DNA that will be translated).  Exits with an error if neither or both
    are provided.
    """
    from ..application.cli_services import resolve_organism_arg

    protein: str | None = getattr(args, "protein", None)
    sequence: str | None = getattr(args, "sequence", None)

    if protein and sequence:
        print(error_msg("Error: provide --protein OR --sequence, not both."), file=sys.stderr)
        sys.exit(1)
    if not protein and not sequence:
        print(error_msg("Error: provide --protein or --sequence."), file=sys.stderr)
        sys.exit(1)

    if protein:
        # Validate amino-acid characters
        valid_aa = set("ACDEFGHIKLMNPQRSTVWYX*")
        cleaned = protein.upper().strip()
        invalid = set(cleaned) - valid_aa
        if invalid:
            print(error_msg(f"Error: invalid amino-acid characters: {', '.join(sorted(invalid))}"),
                  file=sys.stderr)
            sys.exit(1)
        return cleaned

    # Translate DNA → protein
    from biocompiler.expression.translation import translate
    dna = sequence.upper().strip()
    dna = "".join(c for c in dna if c in "ACGT")
    if len(dna) < 3:
        print(error_msg("Error: DNA sequence too short to translate."), file=sys.stderr)
        sys.exit(1)
    protein_seq = translate(dna)
    if not protein_seq:
        print(error_msg("Error: DNA sequence could not be translated."), file=sys.stderr)
        sys.exit(1)
    print(dim(f"Translated DNA ({len(dna)} bp) → protein ({len(protein_seq)} aa)"))
    return protein_seq


def _get_organism(args: argparse.Namespace) -> str:
    """Return the organism name from *args*, defaulting to *Homo_sapiens*."""
    return getattr(args, "organism", None) or "Homo_sapiens"


def _resolve_organism_arg(args: argparse.Namespace) -> str:
    """Resolve the organism from --organism or --species (alias)."""
    from ..application.cli_services import resolve_organism_arg as _resolve
    raw = getattr(args, "organism", None) or getattr(args, "species", None) or "Homo_sapiens"
    species = getattr(args, "species", None)
    return _resolve(raw, species)


def _resolve_source_organism_arg(args: argparse.Namespace) -> str | None:
    """Resolve the source organism from --source-organism."""
    from ..application.cli_services import resolve_source_organism_arg as _resolve
    raw = getattr(args, "source_organism", None)
    return _resolve(raw)


# ── Command: optimize ────────────────────────────────────────────────────────

def cmd_optimize(args: argparse.Namespace) -> None:
    """Handle the 'optimize' command — thin handler calling optimization service."""
    from ..application.cli_services import run_optimization
    from ..application.cli_services import format_optimization_json
    from biocompiler.sequence.scanner import gc_content as _gc_content
    from biocompiler.provenance.certificate import compute_certificate as _compute_cert
    from .. import __version__

    verbose: bool = getattr(args, "verbose", False)
    seed: Optional[int] = getattr(args, "seed", None)
    organism = _resolve_organism_arg(args)
    no_splice_check = getattr(args, "no_splice_check", False)
    source_organism = _resolve_source_organism_arg(args)
    therapeutic = getattr(args, "therapeutic", False)

    # ── Obtain input sequence ──────────────────────────────────────────────
    protein_seq: str | None = getattr(args, "protein", None)
    input_fasta: str | None = getattr(args, "input", None)

    if protein_seq and input_fasta:
        print(error_msg("Error: provide PROTEIN positional arg OR --input, not both."), file=sys.stderr)
        sys.exit(1)

    protein: str | None = None
    seq: str | None = None

    if protein_seq:
        # Validate amino-acid characters
        valid_aa = set("ACDEFGHIKLMNPQRSTVWYX*")
        cleaned = protein_seq.upper().strip()
        invalid = set(cleaned) - valid_aa
        if invalid:
            print(error_msg(f"Error: invalid amino-acid characters: {', '.join(sorted(invalid))}"),
                  file=sys.stderr)
            sys.exit(1)
        protein = cleaned
    elif input_fasta:
        seq = read_fasta(input_fasta)
        if len(seq) < 3:
            print("Error: Sequence too short for optimization.", file=sys.stderr)
            sys.exit(1)
    else:
        print(error_msg("Error: provide PROTEIN positional argument or --input FASTA file."), file=sys.stderr)
        sys.exit(1)

    enzymes: List[str] = []
    if getattr(args, "enzymes", None):
        enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]

    gc_lo = getattr(args, "gc_lo", 0.30)
    gc_hi = getattr(args, "gc_hi", 0.70)
    strategy = getattr(args, "strategy", "hybrid") or "hybrid"
    use_codon_pair_bias = getattr(args, "codon_pair_bias", False)
    organism_domain_raw = getattr(args, "organism_domain", "auto") or "auto"

    # ── Call service ──────────────────────────────────────────────────────
    strategy_label = {
        "csp": "CSP solver",
        "constraint_first": "constraint-first",
        "hybrid": "hybrid (default)",
    }.get(strategy, strategy)

    if verbose:
        print(dim(f"  Strategy: {strategy_label}"))
        print(dim(f"  Organism: {organism}"))
        if no_splice_check:
            print(dim("  Splice check: DISABLED (--no-splice-check)"))
        if use_codon_pair_bias:
            print(dim("  Codon-pair bias: ENABLED"))

    with ProgressStep(f"Optimizing ({strategy_label})", verbose=verbose):
        result = run_optimization(
            protein=protein,
            input_seq=seq,
            organism=organism,
            strategy=strategy,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            enzymes=enzymes,
            no_splice_check=no_splice_check,
            use_codon_pair_bias=use_codon_pair_bias,
            source_organism=source_organism,
            therapeutic=therapeutic,
            organism_domain_raw=organism_domain_raw,
            splice_low=getattr(args, "splice_low", 3.0),
            splice_high=getattr(args, "splice_high", 6.0),
            avoid_gt=getattr(args, "avoid_gt", True),
            seed=seed,
        )

    # ── JSON output ────────────────────────────────────────────────────────
    output_json = getattr(args, "json", False)
    if output_json:
        print(format_optimization_json(
            result,
            no_splice_check=no_splice_check,
            source_organism=source_organism,
            therapeutic=therapeutic,
        ))
        return

    # ── Text output ────────────────────────────────────────────────────────
    if input_fasta:
        # Legacy FASTA mode: write output files
        input_base = os.path.splitext(input_fasta)[0]
        out_fasta = args.output if args.output else f"{input_base}_optimized.fasta"
        out_cert = args.certificate if args.certificate else f"{input_base}_certificate.txt"

        write_fasta(out_fasta, result.optimized, header=f"optimized|{organism}")
        write_certificate(out_cert, result.cert_text)

        cert_level = _compute_cert(result.pred_results) if result.pred_results else None
        print(f"Optimization complete.")
        print(f"  Input:      {input_fasta} ({len(seq)} bp)")
        print(f"  Output:     {out_fasta} ({len(result.optimized)} bp)")
        if cert_level:
            print(f"  Certificate: {out_cert} ({cert_level.value})")
        if result.codon_pair_bias is not None:
            print(f"  Codon-pair bias: {result.codon_pair_bias:.4f}")
        print()
        if result.cert_text:
            print(result.cert_text)

        # Biosecurity report (legacy FASTA-input mode)
        if getattr(args, "biosecurity_report", False):
            from biocompiler.export.core import format_biosecurity_report
            pred_results_for_report = result.pred_results if result.pred_results else []
            report = format_biosecurity_report(
                sequence=result.optimized,
                organism=organism,
                cai=result.cai,
                gc=result.gc_content,
                type_results=pred_results_for_report,
            )
            print()
            print(report)

        # Provenance tracking
        if getattr(args, "provenance", False):
            from ..provenance import ProvenanceTracker, OptimizationRecord
            from datetime import datetime, timezone as _tz

            tracker = ProvenanceTracker(seed=seed or 0)
            out_provenance = f"{input_base}_provenance.json"

            opt_record = OptimizationRecord(
                input_sequence=seq,
                output_sequence=result.optimized,
                organism=organism,
                constraints_applied=[p.name for p in result.pred_results] if result.pred_results else [],
                mutations_made=[],
                solver_backend=strategy,
                solve_time=0.0,
                seed_used=seed,
                timestamp=datetime.now(_tz.utc).isoformat(),
                biocompiler_version=__version__,
            )
            tracker.add_optimization_record(opt_record)

            with open(out_provenance, "w") as f:
                f.write(tracker.to_json())

            print(f"  Provenance:  {out_provenance}")
    else:
        # protein mode: print directly
        print()
        print(section_header("═" * 60))
        print(section_header("  Optimization Result"))
        print(section_header("═" * 60))
        print(f"  Organism       : {organism}")
        print(f"  Strategy       : {strategy}")
        print(f"  Protein length : {len(protein)} aa")
        print(f"  Sequence length: {len(result.optimized)} bp")
        print(f"  GC content     : {result.gc_content:.4f}")
        if no_splice_check:
            print(f"  Splice check   : DISABLED")
        if source_organism:
            print(f"  Source organism : {source_organism}")
        if therapeutic:
            print(f"  Therapeutic    : YES")
        if result.codon_pair_bias is not None:
            print(f"  Codon-pair bias: {result.codon_pair_bias:.4f}")

        if result.pred_results:
            cert_level = _compute_cert(result.pred_results)
            print()
            print(summary_box("Certificate", verdict_symbol(cert_level.value)))

        if result.cert_text:
            print()
            print(result.cert_text)

        # Biosecurity report (if requested)
        if getattr(args, "biosecurity_report", False):
            from biocompiler.export.core import format_biosecurity_report
            pred_results_for_report = result.pred_results if result.pred_results else []
            report = format_biosecurity_report(
                sequence=result.optimized,
                organism=organism,
                cai=result.cai,
                gc=result.gc_content,
                type_results=pred_results_for_report,
            )
            print()
            print(report)

        # Print optimized sequence
        print()
        print(section_header("  Optimized Sequence"))
        for i in range(0, len(result.optimized), 80):
            print(f"  {result.optimized[i:i+80]}")

        # Save to file if --output given
        out_path = getattr(args, "output", None)
        if out_path:
            write_fasta(out_path, result.optimized, header=f"optimized|{organism}")
            print(success_msg(f"  Saved to: {out_path}"))


# ── Command: batch ─────────────────────────────────────────────────────────

def cmd_batch(args: argparse.Namespace) -> None:
    """Handle the 'batch' command — thin handler calling batch optimization service."""
    from ..application.cli_services import run_batch_optimization
    from ..application.cli_services import format_batch_json
    from .. import __version__

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
        print(error_msg(f"Error: File not found: {proteins_file}"), file=sys.stderr)
        sys.exit(1)

    # Read proteins from file
    proteins: List[tuple] = []
    with open(proteins_file, "r") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2 and all(c in "ACDEFGHIKLMNPQRSTVWYX*" for c in parts[1].upper()):
                proteins.append((parts[0], parts[1].upper()))
            else:
                proteins.append((f"protein_{lineno}", line.upper()))

    if not proteins:
        print(error_msg("Error: No valid protein sequences found in file."), file=sys.stderr)
        sys.exit(1)

    if verbose:
        print(dim(f"  Batch mode: {len(proteins)} proteins"))
        print(dim(f"  Organism: {organism}"))
        print(dim(f"  Strategy: {strategy}"))

    # ── Call service ──────────────────────────────────────────────────────
    result = run_batch_optimization(
        proteins=proteins,
        organism=organism,
        strategy=strategy,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        use_codon_pair_bias=use_codon_pair_bias,
        source_organism=source_organism,
        therapeutic=therapeutic,
    )

    # ── JSON output ────────────────────────────────────────────────────────
    if output_json:
        print(format_batch_json(
            result,
            no_splice_check=no_splice_check,
            total_proteins=len(proteins),
            source_organism=source_organism,
            therapeutic=therapeutic,
        ))
        return

    # ── Text output ────────────────────────────────────────────────────────
    print()
    print(section_header("═" * 60))
    print(section_header("  Batch Optimization Results"))
    print(section_header("═" * 60))
    print(f"  Organism : {organism}")
    print(f"  Strategy : {strategy}")
    print(f"  Total    : {len(proteins)} proteins")
    if source_organism:
        print(f"  Source organism : {source_organism}")
    if therapeutic:
        print(f"  Therapeutic    : YES")
    print()

    for r in result.results:
        if r["status"] == "ok":
            status = success_msg("OK")
            print(f"  {r['name']:<20s} {status}  GC={r['gc_content']:.4f}  CAI={r['cai']:.4f}")
        else:
            status = error_msg("FAIL")
            print(f"  {r['name']:<20s} {status}  {r.get('error', 'unknown error')}")

    # Save output file if requested
    out_path = getattr(args, "output", None)
    if out_path:
        with open(out_path, "w") as f:
            for r in result.results:
                if r["status"] == "ok":
                    f.write(f">{r['name']}|{organism}\n")
                    seq = r["sequence"]
                    for i in range(0, len(seq), 80):
                        f.write(seq[i:i+80] + "\n")
        print(success_msg(f"  Saved to: {out_path}"))


# ── Command: check ───────────────────────────────────────────────────────────

def cmd_check(args: argparse.Namespace) -> None:
    """Handle the 'check' command — thin handler calling export service."""
    from ..application.cli_services import run_check_predicates
    from ..type_system import registry as predicate_registry, PREDICATE_NAMES

    # ── List predicates mode ─────────────────────────────────────────────
    if getattr(args, "list_predicates", False):
        print()
        print(section_header("Available Predicates"))
        registered = predicate_registry.names()
        for name in registered:
            if name in PREDICATE_NAMES:
                idx = PREDICATE_NAMES.index(name) + 1
                print(f"  {idx:>2}. {name}")
            else:
                print(f"   * {name} (registered)")
        print()
        print(f"  Total: {len(registered)} predicates registered")
        print()
        return

    # Input file is required unless --list-predicates is used
    if not getattr(args, "input", None):
        print(error_msg("Error: Input FASTA file is required. Use --list-predicates to see available predicates."), file=sys.stderr)
        sys.exit(1)

    seq = read_fasta(args.input)

    if len(seq) < 3:
        print("Error: Sequence too short for checking.", file=sys.stderr)
        sys.exit(1)

    enzymes: List[str] = []
    if args.enzymes:
        enzymes = [e.strip() for e in args.enzymes.split(",") if e.strip()]

    # Resolve organism
    from ..organisms import resolve_organism as _resolve_org
    organism = _resolve_org(getattr(args, "species", "human"))

    # Parse predicate filter
    predicate_filter = None
    predicate_filter_raw = getattr(args, "predicate", None)
    if predicate_filter_raw:
        predicate_filter = [p.strip() for p in predicate_filter_raw.split(",") if p.strip()]
        registered = set(predicate_registry.names())
        unknown = [p for p in predicate_filter if p not in registered]
        if unknown:
            print(error_msg(f"Error: Unknown predicates: {', '.join(unknown)}"), file=sys.stderr)
            print(dim(f"  Use --list-predicates to see available predicates."), file=sys.stderr)
            sys.exit(1)

    # ── Call service ──────────────────────────────────────────────────────
    result = run_check_predicates(
        seq=seq,
        organism=organism,
        enzymes=enzymes if enzymes else None,
        cryptic_threshold=getattr(args, "splice_low", 3.0),
        uncertain_lo=getattr(args, "splice_low", 3.0),
        predicate_filter=predicate_filter,
        species=args.species,
    )

    print(f"Certificate: {result.cert_level}")
    print()
    print(result.cert_text)


# ── Command: benchmark ───────────────────────────────────────────────────────

def cmd_benchmark(args: argparse.Namespace) -> None:
    """Handle the 'benchmark' command."""
    # Set deterministic seed if provided
    seed: Optional[int] = getattr(args, "seed", None)
    if seed is not None:
        random.seed(seed)
        logger.info("Random seed set to %d for reproducible benchmark", seed)

    from biocompiler.benchmarking.core import (
        run_benchmark, compare_tools,
        REFERENCE_GENES, GENE_PANEL,
        run_structured_benchmarks,
        format_benchmark_report_text, format_benchmark_report_json,
    )

    # ── List gene sets ──
    if getattr(args, "list_gene_sets", False):
        print()
        print(section_header("Available Gene Sets"))
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
        gene_set_upper = gene_set.upper()
        if gene_set_upper == "REFERENCE_GENES":
            report = run_structured_benchmarks(gene_names=list(REFERENCE_GENES.keys()))
            text = format_benchmark_report_text(report)
            print(text)
        elif gene_set_upper in ("HUMAN_THERAPEUTIC", "GENE_PANEL"):
            from biocompiler.benchmarking.core import run_multi_gene_comparison
            results = run_multi_gene_comparison(
                enzymes=enzymes if enzymes else None,
            )
            print()
            print(section_header("═" * 80))
            print(section_header(f"  Benchmark: {gene_set}"))
            print(section_header("═" * 80))
            for r in results:
                tool = r.get("tool", "?")
                gene = r.get("gene", "?")
                cai = r.get("cai", 0.0)
                gc = r.get("gc_content", 0.0)
                violations = r.get("constraint_violations", "?")
                success = r.get("success", False)
                status = success_msg("OK") if success else error_msg("FAIL")
                print(f"  {tool:<20s} {gene:<12s} CAI={cai:.4f} GC={gc:.4f} violations={violations} {status}")
        else:
            print(error_msg(f"Unknown gene set: {gene_set}"), file=sys.stderr)
            print(dim("  Use --list-gene-sets to see available gene sets."), file=sys.stderr)
            sys.exit(1)
    else:
        run_benchmark(
            enzymes=enzymes if enzymes else None,
            splice_low=args.splice_low,
            splice_high=args.splice_high,
        )
        compare_tools()

    # Save output if requested
    if output_file:
        if output_file.lower().endswith(".json"):
            try:
                report = run_structured_benchmarks()
                json_data = format_benchmark_report_json(report)
                with open(output_file, "w") as f:
                    f.write(json_data)
                print(success_msg(f"Benchmark results saved to {output_file}"), file=sys.stderr)
            except Exception as exc:
                print(error_msg(f"Error saving benchmark results: {exc}"), file=sys.stderr)
        else:
            try:
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
                print(success_msg(f"Benchmark results saved to {output_file}"), file=sys.stderr)
            except Exception as exc:
                print(error_msg(f"Error saving benchmark results: {exc}"), file=sys.stderr)


# ── Command: scan ────────────────────────────────────────────────────────────

def cmd_scan(args: argparse.Namespace) -> None:
    """Handle the 'scan' command — scan a sequence for features."""
    from biocompiler.sequence.scanner import scan_sequence
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
            print(error_msg("Error: --quality-only requires --pdb-file."), file=sys.stderr)
            sys.exit(1)
        from ..structure.quality import compute_structure_quality
        if not os.path.isfile(pdb_file):
            print(error_msg(f"Error: PDB file not found: {pdb_file}"), file=sys.stderr)
            sys.exit(1)
        with open(pdb_file, "r") as f:
            pdb_content = f.read()
        with ProgressStep("Assessing structure quality", verbose=getattr(args, "verbose", False)):
            report = compute_structure_quality(pdb_content)
        print_structure_quality(report)
        return

    # Full prediction mode
    protein = _resolve_protein(args)
    from biocompiler.engines.esmfold import predict_structure, is_esmfold_available
    from ..structure.quality import compute_structure_quality

    esmfold_ok = is_esmfold_available()
    if not esmfold_ok:
        print(dim("ESMFold not available — using offline/fallback prediction."))

    with ProgressStep("Predicting structure", verbose=getattr(args, "verbose", False)):
        # predict_structure() returns an ESMFoldResult whose ``pdb_string``
        # attribute holds the PDB text.  Previously the CLI treated the
        # returned result object as a path/string, which crashed
        # compute_structure_quality() and shutil.copy2() — see audit H23.
        result = predict_structure(protein, organism=organism)

    with ProgressStep("Computing quality metrics"):
        report = compute_structure_quality(result.pdb_string)

    print()
    print(section_header("═" * 60))
    print(section_header("  Structure Prediction & Quality Report"))
    print(section_header("═" * 60))
    print(f"  Protein length : {len(protein)} aa")
    print(f"  Organism       : {organism}")
    print(f"  ESMFold        : {'available' if esmfold_ok else 'offline/fallback'}")
    print()

    print_structure_quality(report)

    output = getattr(args, "output", None)
    if output and result.pdb_string:
        # Write the PDB string to the requested output path (the result is
        # an in-memory PDB representation, not a file on disk).
        with open(output, "w") as f:
            f.write(result.pdb_string)
        print(success_msg(f"PDB file saved to {output}"))


# ── Command: stability ───────────────────────────────────────────────────────

def cmd_stability(args: argparse.Namespace) -> None:
    """Analyze protein stability."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)

    from biocompiler.engines.foldx import empirical_stability, scan_mutations

    with ProgressStep("Computing stability", verbose=getattr(args, "verbose", False)):
        result = empirical_stability(protein, organism=organism)

    print()
    print(section_header("═" * 60))
    print(section_header("  Protein Stability Analysis"))
    print(section_header("═" * 60))
    print(f"  Protein length : {len(protein)} aa")
    print(f"  Organism       : {organism}")
    print()

    print(section_header("  Stability Metrics"))
    # FoldXResult exposes ``stability_kcal`` (not ``delta_g``) and stores
    # individual energy components as flat attributes (there is no
    # ``energy_components`` dict).  Build the dict from the known
    # FoldXResult energy-component attribute names — see audit H23.
    dg = getattr(result, "stability_kcal", None)
    print(f"  \u0394G (kcal/mol)   : {dg if dg is not None else 'N/A'}")
    print(f"  Energy components:")
    _energy_attr_names = [
        "interaction_energy",
        "backbone_hbond",
        "sidechain_hbond",
        "van_der_waals",
        "electrostatics",
        "solvation",
        "van_der_waals_clashes",
        "entropy_sidechain",
        "entropy_mainchain",
        "torsional_clash",
        "backbone_clash",
        "helix_dipole",
        "disulfide",
        "electrostatic_kon",
        "partial_covalent",
        "energy_ionisation",
    ]
    for comp_name in _energy_attr_names:
        comp_val = getattr(result, comp_name, None)
        if comp_val is not None:
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
    print(summary_box("Stability Verdict", verdict_symbol(verdict)))

    if getattr(args, "scan_mutations", False):
        positions = getattr(args, "positions", None) or list(range(1, len(protein) + 1))
        with ProgressStep("Scanning mutations", verbose=getattr(args, "verbose", False)):
            mut_results = scan_mutations(protein, positions=positions, organism=organism)
        print_mutation_table(mut_results)


# ── Command: solubility ─────────────────────────────────────────────────────

def cmd_solubility(args: argparse.Namespace) -> None:
    """Analyze protein solubility."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)

    from biocompiler.engines.camsol import compute_solubility, find_solubility_mutations

    with ProgressStep("Computing solubility", verbose=getattr(args, "verbose", False)):
        result = compute_solubility(protein, organism=organism)

    print()
    print(section_header("═" * 60))
    print(section_header("  Protein Solubility Analysis"))
    print(section_header("═" * 60))
    print(f"  Protein length : {len(protein)} aa")
    print(f"  Organism       : {organism}")
    print()

    print(section_header("  Solubility Metrics"))
    # CamSolResult exposes ``primary_score`` (not ``camsol_score``) as the
    # unified-API alias for its solubility score — see audit H23.
    score = getattr(result, "primary_score", None)
    print(f"  CamSol score    : {score if score is not None else 'N/A'}")

    # CamSolResult.aggregation_prone_regions is list[tuple[int, int, float]]
    # (start, end, avg_score) — NOT a list of objects with .start/.end/.score.
    agg_regions = getattr(result, "aggregation_prone_regions", [])
    print(f"  Agg-prone regions: {len(agg_regions)}")
    for region in agg_regions:
        if isinstance(region, (tuple, list)):
            start = region[0] if len(region) > 0 else "?"
            end = region[1] if len(region) > 1 else "?"
            rscore = region[2] if len(region) > 2 else "N/A"
            seq_r = ""
        else:
            start = getattr(region, "start", "?")
            end = getattr(region, "end", "?")
            seq_r = getattr(region, "sequence", "")
            rscore = getattr(region, "score", "N/A")
        print(f"    [{start}-{end}] {seq_r}  (score: {rscore})")

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
    print(summary_box("Solubility Verdict", verdict_symbol(verdict)))

    recommendations = getattr(result, "recommendations", [])
    if recommendations:
        print()
        print(section_header("  Recommendations"))
        for i, rec in enumerate(recommendations, 1):
            print(f"  {i}. {rec}")

    if getattr(args, "find_mutations", False):
        with ProgressStep("Finding solubility-improving mutations",
                          verbose=getattr(args, "verbose", False)):
            mut_results = find_solubility_mutations(protein, organism=organism)
        print()
        print(section_header("  Solubility-Improving Mutations"))
        header = f"  {'Position':>8s}  {'Original':>8s}  {'Mutant':>8s}  {'\u0394CamSol':>10s}  {'Effect':>12s}"
        print(header)
        print(f"  {'─' * 8}  {'─' * 8}  {'─' * 8}  {'─' * 10}  {'─' * 12}")
        for mut in mut_results:
            pos = getattr(mut, "position", "?")
            orig = getattr(mut, "original_aa", "?")
            mutant = getattr(mut, "mutant_aa", "?")
            # MutationResult exposes ``delta_score`` (not ``delta_camsol``)
            # as the unified-API alias for the score delta — see audit H23.
            delta = getattr(mut, "delta_score", None)
            delta_str = f"{delta:+.3f}" if isinstance(delta, (int, float)) else "N/A"
            print(f"  {pos:>8}  {orig:>8}  {mutant:>8}  {delta_str:>10}  {success_msg('improving')}")


# ── Command: immunogenicity ─────────────────────────────────────────────────

def cmd_immunogenicity(args: argparse.Namespace) -> None:
    """Analyze and reduce immunogenicity."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)
    source_organism = _resolve_source_organism_arg(args)
    therapeutic = getattr(args, "therapeutic", False)

    from biocompiler.immunogenicity.core import compute_immunogenicity
    from biocompiler.immunogenicity.deimmunization import deimmunize

    mhc_alleles = getattr(args, "mhc_alleles", None) or []

    with ProgressStep("Computing immunogenicity", verbose=getattr(args, "verbose", False)):
        result = compute_immunogenicity(protein, organism=organism, mhc_alleles=mhc_alleles)

    print()
    print(section_header("═" * 60))
    print(section_header("  Immunogenicity Analysis"))
    print(section_header("═" * 60))
    print(f"  Protein length : {len(protein)} aa")
    print(f"  Organism       : {organism}")
    if source_organism:
        print(f"  Source organism : {source_organism}")
    if therapeutic:
        print(f"  Therapeutic    : YES")
    if mhc_alleles:
        print(f"  MHC alleles    : {', '.join(mhc_alleles)}")
    print()

    print(section_header("  Immunogenicity Metrics"))
    score = getattr(result, "immunogenicity_score", None)
    print(f"  Score           : {score if score is not None else 'N/A'}")

    epitopes = getattr(result, "epitopes", [])
    print(f"  Epitopes found  : {len(epitopes)}")
    for i, epi in enumerate(epitopes, 1):
        start = getattr(epi, "start", "?")
        end = getattr(epi, "end", "?")
        seq_e = getattr(epi, "sequence", "")
        escore = getattr(epi, "score", "N/A")
        allele = getattr(epi, "allele", "N/A")
        print(f"    {i:3d}. [{start}-{end}] {seq_e}  allele={allele}  score={escore}")

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
    print(summary_box("Immunogenicity Verdict", verdict_symbol(verdict)))

    # Deimmunization
    if getattr(args, "deimmunize", False):
        max_mut = getattr(args, "max_mutations", 10)
        blosum_min = getattr(args, "blosum62_min", 1)

        with ProgressStep("Running deimmunization optimization",
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
        print(section_header("  Deimmunization Results"))
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
            print(dim(f"  Optimized protein: {new_protein[:60]}{'...' if len(new_protein) > 60 else ''}"))

        if new_score is not None and isinstance(new_score, (int, float)):
            if new_score <= target:
                new_verdict = "PASS"
            elif new_score <= target + 0.2:
                new_verdict = "LIKELY_PASS"
            else:
                new_verdict = "UNCERTAIN"
            print()
            print(summary_box("Post-deimmunization Verdict", verdict_symbol(new_verdict)))


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

    from ..structure.report import (
        assess_protein,
        format_assessment_text,
        format_assessment_json,
        format_assessment_html,
    )

    pdb_content: str | None = None
    if pdb_file:
        if not os.path.isfile(pdb_file):
            print(error_msg(f"Error: PDB file not found: {pdb_file}"), file=sys.stderr)
            sys.exit(1)
        with open(pdb_file, "r") as f:
            pdb_content = f.read()

    with ProgressStep("Running comprehensive assessment",
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

    if fmt == "json":
        formatted = format_assessment_json(report)
        print(formatted)
        if output:
            with open(output, "w") as f:
                f.write(formatted)
            print(success_msg(f"JSON report saved to {output}"), file=sys.stderr)
    elif fmt == "html":
        formatted = format_assessment_html(report)
        print(formatted)
        if output:
            with open(output, "w") as f:
                f.write(formatted)
            print(success_msg(f"HTML report saved to {output}"), file=sys.stderr)
    else:
        formatted = format_assessment_text(report)
        print(formatted)
        if output:
            with open(output, "w") as f:
                f.write(formatted)
            print(success_msg(f"Text report saved to {output}"), file=sys.stderr)


# ── Command: validate-cai ─────────────────────────────────────────────────

def cmd_validate_cai(args: argparse.Namespace) -> None:
    """Validate CAI computation against published ground-truth values."""
    organism = _get_organism(args)

    try:
        from ..validation.ground_truth import GROUND_TRUTH_DATA, validate_against_ground_truth
    except ImportError:
        print(error_msg("Error: validation.ground_truth module not available."), file=sys.stderr)
        sys.exit(1)

    print()
    print(section_header("═" * 60))
    print(section_header("  CAI Validation Against Published Values"))
    print(section_header("═" * 60))
    print(f"  Organism : {organism}")
    print()

    entries = [e for e in GROUND_TRUTH_DATA if e.organism == organism]
    if not entries:
        available = sorted(set(e.organism for e in GROUND_TRUTH_DATA))
        print(error_msg(f"No ground-truth entries for organism '{organism}'."))
        print(dim(f"  Available organisms: {', '.join(available)}"))
        sys.exit(1)

    all_pass = True
    for entry in entries:
        result = validate_against_ground_truth(
            optimized_sequence=entry.published_sequence,
            gene_name=entry.gene_name,
            organism=entry.organism,
        )
        status = success_msg("PASS") if result.matches_expected else error_msg("FAIL")
        print(f"  {entry.gene_name:<15s} {status}")
        print(f"    CAI diff: {result.cai_difference:.4f}, GC diff: {result.gc_difference:.4f}")
        if not result.matches_expected:
            all_pass = False

    print()
    if all_pass:
        print(summary_box("CAI Validation", verdict_symbol("PASS")))
    else:
        print(summary_box("CAI Validation", verdict_symbol("FAIL")))


# ── Command: validate-maxentscan ──────────────────────────────────────────

def cmd_validate_maxentscan(args: argparse.Namespace) -> None:
    """Validate MaxEntScan scoring against published splice-site scores."""
    try:
        from biocompiler.sequence.maxentscan import score_donor, score_acceptor, max_donor_score, max_acceptor_score
    except ImportError:
        print(error_msg("Error: maxentscan module not available."), file=sys.stderr)
        sys.exit(1)

    donor_tests = [
        ("CAGGTAAGT", 3, (3.0, 12.0)),
        ("AAGGTGAGT", 2, (3.0, 12.0)),
        ("TTGGTAAAT", 2, (0.0, 8.0)),
    ]

    acceptor_tests = [
        ("TTTTTTTTTTTTTTTTTTTTCAGG", 20, (3.0, 14.0)),
        ("CCCCCCCCCCCCCCCCCCCCCAGA", 20, (0.0, 10.0)),
    ]

    print()
    print(section_header("═" * 60))
    print(section_header("  MaxEntScan Validation"))
    print(section_header("═" * 60))
    print()

    all_pass = True

    print(section_header("  Donor Site Scoring"))
    for seq, pos, (lo, hi) in donor_tests:
        score = score_donor(seq, pos)
        in_range = lo <= score <= hi
        status = success_msg("PASS") if in_range else error_msg("FAIL")
        print(f"  {seq}  pos={pos}  score={score:.2f}  expected=[{lo},{hi}]  {status}")
        if not in_range:
            all_pass = False

    print()
    print(section_header("  Acceptor Site Scoring"))
    for seq, pos, (lo, hi) in acceptor_tests:
        score = score_acceptor(seq, pos)
        in_range = lo <= score <= hi
        status = success_msg("PASS") if in_range else error_msg("FAIL")
        print(f"  {seq[:10]}...  pos={pos}  score={score:.2f}  expected=[{lo},{hi}]  {status}")
        if not in_range:
            all_pass = False

    print()
    print(section_header("  Edge-Case Validation"))
    non_donor_score = score_donor("CAGATACGT", 3)
    non_donor_ok = non_donor_score < 3.0
    status = success_msg("PASS") if non_donor_ok else error_msg("FAIL")
    print(f"  Non-donor (no GT)  score={non_donor_score:.2f}  expected<3.0  {status}")
    if not non_donor_ok:
        all_pass = False

    print()
    if all_pass:
        print(summary_box("MaxEntScan Validation", verdict_symbol("PASS")))
    else:
        print(summary_box("MaxEntScan Validation", verdict_symbol("FAIL")))


# ── Command: whatif ───────────────────────────────────────────────────────

def cmd_whatif(args: argparse.Namespace) -> None:
    """Run what-if analysis on a protein sequence."""
    protein = _resolve_protein(args)
    organism = _get_organism(args)

    print()
    print(section_header("═" * 60))
    print(section_header("  What-If Analysis"))
    print(section_header("═" * 60))
    print(f"  Protein length : {len(protein)} aa")
    print(f"  Organism       : {organism}")
    print()

    # 1. Codon optimization what-if
    try:
        from biocompiler.optimizer import optimize_sequence
        with ProgressStep("Running codon optimization what-if",
                          verbose=getattr(args, "verbose", False)):
            opt_result = optimize_sequence(
                target_protein=protein,
                organism=organism,
            )
        print(section_header("  Codon Optimization What-If"))
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
        print(dim("  Codon optimization module not available; skipping."))
        print()

    # 2. Stability what-if
    try:
        from biocompiler.engines.foldx import empirical_stability
        with ProgressStep("Running stability what-if",
                          verbose=getattr(args, "verbose", False)):
            stab_result = empirical_stability(protein, organism=organism)
        dg = getattr(stab_result, "delta_g", None)
        print(section_header("  Stability What-If"))
        print(f"  ΔG (kcal/mol)  : {dg if dg is not None else 'N/A'}")
        print()
    except ImportError:
        print(dim("  Stability module not available; skipping."))
        print()

    # 3. Solubility what-if
    try:
        from biocompiler.engines.camsol import compute_solubility
        with ProgressStep("Running solubility what-if",
                          verbose=getattr(args, "verbose", False)):
            sol_result = compute_solubility(protein, organism=organism)
        camsol_score = getattr(sol_result, "camsol_score", None)
        agg_regions = getattr(sol_result, "aggregation_prone_regions", [])
        print(section_header("  Solubility What-If"))
        print(f"  CamSol score   : {camsol_score if camsol_score is not None else 'N/A'}")
        print(f"  Agg-prone regions: {len(agg_regions)}")
        print()
    except ImportError:
        print(dim("  Solubility module not available; skipping."))
        print()

    # 4. Immunogenicity what-if
    try:
        from biocompiler.immunogenicity.core import compute_immunogenicity
        with ProgressStep("Running immunogenicity what-if",
                          verbose=getattr(args, "verbose", False)):
            imm_result = compute_immunogenicity(protein, organism=organism)
        imm_score = getattr(imm_result, "immunogenicity_score", None)
        epitopes = getattr(imm_result, "epitopes", [])
        print(section_header("  Immunogenicity What-If"))
        print(f"  Score          : {imm_score if imm_score is not None else 'N/A'}")
        print(f"  Epitopes found : {len(epitopes)}")
        print()
    except ImportError:
        print(dim("  Immunogenicity module not available; skipping."))
        print()

    # 5. Literature validation what-if
    try:
        from biocompiler.validation.literature_validation import evaluate_case, ALL_LITERATURE_CASES
        matching = [c for c in ALL_LITERATURE_CASES
                     if c.sequence_type == "protein" and c.sequence.strip().upper() == protein.upper()]
        if matching:
            print(section_header("  Literature Validation What-If"))
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
            print(dim("  No matching literature cases for this protein sequence."))
            print()
    except ImportError:
        print(dim("  Literature validation module not available; skipping."))
        print()

    print(summary_box("What-If Complete", verdict_symbol("PASS")))


# ── Command: explain ───────────────────────────────────────────────────────

def cmd_explain(args: argparse.Namespace) -> None:
    """Explain why a specific codon was chosen in a previous optimization."""
    from ..provenance import ProvenanceTracker

    if not os.path.isfile(args.file):
        print(error_msg(f"Error: File not found: {args.file}"), file=sys.stderr)
        sys.exit(1)

    with open(args.file, "r") as f:
        data = json.load(f)

    tracker = ProvenanceTracker.from_dict(data)
    position = args.position

    decisions = tracker.get_decisions_for_position(position)

    print()
    print(section_header("═" * 60))
    print(section_header("  Provenance Explanation"))
    print(section_header("═" * 60))
    print(f"  Position:  {position}")
    print(f"  File:      {args.file}")
    print()

    if not decisions:
        print(dim(f"  No decisions recorded for position {position}."))
        print(dim("  This may indicate the position was not modified during optimization"))
        print(dim("  or provenance tracking was not enabled for that position."))
    else:
        print(section_header(f"  Decisions at position {position}"))
        print(f"  {'Type':<22s}  {'Chosen':<10s}  {'Alternatives':<30s}")
        print(f"  {'─' * 22}  {'─' * 10}  {'─' * 30}")
        for d in decisions:
            alts = ", ".join(d.alternatives_considered[:5])
            if len(d.alternatives_considered) > 5:
                alts += f" (+{len(d.alternatives_considered) - 5} more)"
            print(f"  {d.decision_type:<22s}  {d.chosen_value:<10s}  {alts}")
        print()
        print(section_header("  Rationale"))
        for d in decisions:
            print(f"  [{d.decision_type}] {d.rationale}")
            if d.constraint_context:
                ctx_parts = [f"{k}={v}" for k, v in d.constraint_context.items()]
                print(dim(f"    Context: {', '.join(ctx_parts)}"))

    print()


# ── Command: report ────────────────────────────────────────────────────────

def cmd_report(args: argparse.Namespace) -> None:
    """Generate provenance report from saved trail."""
    from ..provenance import ProvenanceTracker, generate_provenance_report

    if not os.path.isfile(args.file):
        print(error_msg(f"Error: File not found: {args.file}"), file=sys.stderr)
        sys.exit(1)

    with open(args.file, "r") as f:
        data = json.load(f)

    tracker = ProvenanceTracker.from_dict(data)
    records = tracker.get_optimization_records()
    fmt = getattr(args, "format", "text") or "text"

    print()
    print(section_header("═" * 60))
    print(section_header("  Provenance Report"))
    print(section_header("═" * 60))
    print(f"  Source:      {args.file}")
    print(f"  Seed:        {tracker.seed}")
    print(f"  Decisions:   {len(tracker.get_full_audit_trail())}")
    print(f"  Opt records: {len(records)}")
    print()

    if fmt == "json":
        print(json.dumps(data, indent=2, sort_keys=True))
    elif fmt == "markdown":
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
        trail = tracker.get_full_audit_trail()
        if trail:
            print("## Decision Trail\n")
            print("| # | Type | Position | Chosen | Alternatives | Rationale |")
            print("|---|------|----------|--------|-------------|-----------|")
            for i, d in enumerate(trail, 1):
                alts = ", ".join(d.alternatives_considered[:3])
                if len(d.alternatives_considered) > 3:
                    alts += "..."
                rationale = d.rationale.replace("|", "\\|")
                print(f"| {i} | {d.decision_type} | {d.position} | {d.chosen_value} | {alts} | {rationale} |")
            print()
    else:
        report_text = generate_provenance_report(records)
        print(report_text)
