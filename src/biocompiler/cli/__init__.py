"""
BioCompiler CLI v1.0.0
=======================
Command-line interface for certified gene optimization and protein analysis.

Commands:
  optimize            Optimize a protein sequence for a target organism
  batch               Batch-optimize proteins from a file
  check               Read FASTA, evaluate all registered predicates, print certificate
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

Refactored (Wave 4b): This package decomposes the original monolithic cli.py into:
  - parser.py     — argparse definitions
  - commands.py   — thin command handlers calling service functions
  - formatters.py — output formatting (colour, tables, boxes, file I/O)
"""

from __future__ import annotations

import sys
from typing import List, Optional

# Re-export public API for backward compatibility
from .parser import build_parser
from .formatters import (
    colorize,
    _read_fasta, _write_fasta, read_fasta, write_fasta,
    _write_certificate, write_certificate,
    _section_header, section_header,
    _verdict_symbol, verdict_symbol,
    _error_msg, error_msg,
    _success_msg, success_msg,
    _dim, dim,
    _summary_box, summary_box,
    _supports_color, supports_color,
    _ProgressStep, ProgressStep,
)
from .commands import (
    cmd_optimize,
    cmd_batch,
    cmd_check,
    cmd_benchmark,
    cmd_scan,
    cmd_structure,
    cmd_stability,
    cmd_solubility,
    cmd_immunogenicity,
    cmd_assess,
    cmd_validate_cai,
    cmd_validate_maxentscan,
    cmd_whatif,
    cmd_explain,
    cmd_report,
    _resolve_protein,
    _get_organism,
)

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
    "_read_fasta",
    "_write_fasta",
    "read_fasta",
    "write_fasta",
    "_resolve_protein",
    "_get_organism",
    "_write_certificate",
    "_section_header",
    "_verdict_symbol",
    "_error_msg",
    "_success_msg",
    "_dim",
    "_summary_box",
    "_supports_color",
    "_ProgressStep",
]


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
        from ..api import app, set_no_auth_flag
        if getattr(args, "no_auth", False):
            set_no_auth_flag()
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
