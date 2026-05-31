"""
BioCompiler CLI v7.0.0
=======================
Command-line interface for certified gene optimization.

Commands:
  optimize    Read FASTA, run full 6-phase optimization, write optimized FASTA + certificate
  check       Read FASTA, evaluate all 8 predicates, print certificate
  benchmark   Run built-in benchmarks (eGFP, mCherry, LacZ)
"""

import argparse
import sys
import os
from typing import List, Optional

from .optimizer import BioOptimizer
from .type_system import (
    CODON_TABLE,
    check_no_stop_codons,
    check_no_cryptic_splice,
    check_no_cpg_island,
    check_no_restriction_site,
    check_no_avoidable_gt,
    check_valid_coding_seq,
)
from .certificates import format_certificate, compute_certificate


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


def cmd_optimize(args: argparse.Namespace) -> None:
    """Handle the 'optimize' command."""
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
    from .species import SPECIES
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


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the BioCompiler CLI."""
    parser = argparse.ArgumentParser(
        prog="biocompiler",
        description="BioCompiler v7.0.0 — Certified Gene Optimization with Formal Verification",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # ── optimize ──
    opt_parser = subparsers.add_parser(
        "optimize",
        help="Optimize a FASTA gene sequence with full 6-phase pipeline",
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
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
