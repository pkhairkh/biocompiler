"""
BioCompiler CLI — Command-Line Interface

Usage:
    biocompiler check --sequence ATGGTGCATCTG... --exons 0,92 273,495 1346,1608
    biocompiler optimize --protein MVHLTPEEK... --organism Homo_sapiens
    biocompiler verify --certificate certificate.json
    biocompiler scan --sequence ATGGTGCATCTG...
"""

import argparse
import json
import sys
import logging
from . import __version__
from .scanner import scan_sequence, gc_content
from .translation import translate, compute_cai, find_orfs
from .type_system import evaluate_all_predicates
from .certificate import generate_certificate, verify_certificate
from .optimization import optimize_sequence
from .types import Verdict, combined_verdict
from .exceptions import BioCompilerError, CertificateGenerationError

logger = logging.getLogger("biocompiler")


def _setup_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def cmd_check(args):
    """Run type checking on a DNA sequence."""
    _setup_logging(args.verbose)
    seq = args.sequence.upper()
    exon_boundaries = _parse_exons(args.exons)

    print(f"Sequence length: {len(seq)} nt")
    print(f"GC content: {gc_content(seq):.1%}")
    print(f"Exon boundaries: {exon_boundaries}")
    print()

    # Scan
    tokens = scan_sequence(seq, args.enzymes.split(",") if args.enzymes else None)
    print(f"Scanner found {len(tokens)} tokens")
    donor_count = sum(1 for t in tokens if t.element_type == "splice_donor")
    acceptor_count = sum(1 for t in tokens if t.element_type == "splice_acceptor")
    print(f"  Donor sites: {donor_count}")
    print(f"  Acceptor sites: {acceptor_count}")
    print()

    # Translate
    if exon_boundaries:
        coding_seq = "".join(seq[start:end] for start, end in exon_boundaries)
    else:
        coding_seq = seq
    protein = translate(coding_seq)
    print(f"Protein: {protein[:50]}{'...' if len(protein) > 50 else ''} ({len(protein)} aa)")
    print()

    # Type check
    results = evaluate_all_predicates(
        seq=seq,
        known_exon_boundaries=exon_boundaries or [(0, len(seq))],
        organism=args.organism,
        gc_lo=args.gc_lo,
        gc_hi=args.gc_hi,
        cai_threshold=args.cai_threshold,
    )

    for r in results:
        symbol = {"PASS": "PASS", "FAIL": "FAIL", "UNCERTAIN": "UNCERTAIN"}[r.verdict.value]
        print(f"  [{symbol}] {r.predicate}")
        if r.violation:
            print(f"         {r.violation}")
    print()

    overall = combined_verdict([r.verdict for r in results])
    print(f"Overall verdict: {overall.value}")

    if overall == Verdict.PASS and args.output:
        try:
            cert = generate_certificate(
                coding_seq, results,
                {
                    "gene": args.gene or "unknown",
                    "organism": args.organism,
                    "exon_boundaries": exon_boundaries,
                    "gc_lo": args.gc_lo,
                    "gc_hi": args.gc_hi,
                    "cai_threshold": args.cai_threshold,
                },
            )
            with open(args.output, "w") as f:
                json.dump(cert.to_dict(), f, indent=2)
            print(f"Certificate saved to: {args.output}")
        except CertificateGenerationError as e:
            print(f"Certificate generation failed: {e}")


def cmd_optimize(args):
    """Optimize a DNA sequence for a target protein."""
    _setup_logging(args.verbose)
    print(f"Optimizing for protein ({len(args.protein)} aa), organism={args.organism}")

    result = optimize_sequence(
        target_protein=args.protein,
        organism=args.organism,
        gc_lo=args.gc_lo,
        gc_hi=args.gc_hi,
    )

    print(f"Sequence: {result.sequence[:60]}{'...' if len(result.sequence) > 60 else ''}")
    print(f"Length: {len(result.sequence)} bp")
    print(f"CAI: {result.cai:.4f}")
    print(f"GC: {result.gc_content:.1%}")
    print(f"Satisfied predicates: {result.satisfied_predicates}")
    print(f"Failed predicates: {result.failed_predicates}")
    print(f"Fallback solver used: {result.fallback_used}")

    if args.output:
        with open(args.output, "w") as f:
            json.dump({
                "sequence": result.sequence,
                "protein": result.protein,
                "cai": result.cai,
                "gc_content": result.gc_content,
                "satisfied": result.satisfied_predicates,
                "failed": result.failed_predicates,
            }, f, indent=2)
        print(f"Result saved to: {args.output}")


def cmd_verify(args):
    """Verify a certificate."""
    _setup_logging(args.verbose)
    with open(args.certificate) as f:
        cert_dict = json.load(f)

    status, failures = verify_certificate(cert_dict)
    print(f"Verification status: {status}")
    if failures:
        for f in failures:
            print(f"  FAILURE: {f}")
    else:
        print("  All checks passed.")


def cmd_scan(args):
    """Scan a sequence for motifs."""
    _setup_logging(args.verbose)
    seq = args.sequence.upper()
    tokens = scan_sequence(seq, args.enzymes.split(",") if args.enzymes else None)

    print(f"Sequence length: {len(seq)} nt")
    print(f"Tokens found: {len(tokens)}")
    for t in tokens:
        frame_str = f" frame={t.frame}" if t.frame is not None else ""
        strand_str = f" strand={t.strand}" if t.strand != "+" else ""
        print(f"  {t.element_type} @ {t.position}: {t.match_sequence} (score={t.score:.1f}){frame_str}{strand_str}")

    if args.find_orfs:
        orfs = find_orfs(seq)
        print(f"\nORFs found: {len(orfs)}")
        for orf in orfs:
            print(f"  {orf['strand']} strand, frame {orf['frame']}: "
                  f"{orf['start']}-{orf['end']} ({orf['length']} aa)")


def _parse_exons(exons_str: str | None) -> list[tuple[int, int]]:
    """Parse exon boundaries from comma-separated 'start,end start,end ...' format."""
    if not exons_str:
        return []
    boundaries = []
    for pair in exons_str.split():
        parts = pair.split(",")
        if len(parts) == 2:
            boundaries.append((int(parts[0]), int(parts[1])))
    return boundaries


def main():
    parser = argparse.ArgumentParser(
        prog="biocompiler",
        description="BioCompiler — Machine-verified gene design",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # check
    p_check = subparsers.add_parser("check", help="Type-check a DNA sequence")
    p_check.add_argument("--sequence", required=True, help="DNA sequence")
    p_check.add_argument("--exons", help="Exon boundaries as 'start,end start,end ...'")
    p_check.add_argument("--gene", help="Gene name")
    p_check.add_argument("--organism", default="Homo_sapiens", help="Target organism")
    p_check.add_argument("--gc-lo", type=float, default=0.30, help="Min GC content")
    p_check.add_argument("--gc-hi", type=float, default=0.70, help="Max GC content")
    p_check.add_argument("--cai-threshold", type=float, default=0.5, help="Min CAI")
    p_check.add_argument("--enzymes", help="Comma-separated restriction enzymes to check")
    p_check.add_argument("--output", "-o", help="Output certificate file path")
    p_check.set_defaults(func=cmd_check)

    # optimize
    p_opt = subparsers.add_parser("optimize", help="Optimize DNA sequence for a protein")
    p_opt.add_argument("--protein", required=True, help="Target protein sequence")
    p_opt.add_argument("--organism", default="Homo_sapiens", help="Target organism")
    p_opt.add_argument("--gc-lo", type=float, default=0.30)
    p_opt.add_argument("--gc-hi", type=float, default=0.70)
    p_opt.add_argument("--output", "-o", help="Output result file path")
    p_opt.set_defaults(func=cmd_optimize)

    # verify
    p_verify = subparsers.add_parser("verify", help="Verify a certificate")
    p_verify.add_argument("--certificate", required=True, help="Certificate JSON file")
    p_verify.set_defaults(func=cmd_verify)

    # scan
    p_scan = subparsers.add_parser("scan", help="Scan a sequence for biological motifs")
    p_scan.add_argument("--sequence", required=True, help="DNA sequence")
    p_scan.add_argument("--enzymes", help="Comma-separated restriction enzymes")
    p_scan.add_argument("--find-orfs", action="store_true", help="Find open reading frames")
    p_scan.set_defaults(func=cmd_scan)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        args.func(args)
    except BioCompilerError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        logger.exception("Unexpected error")
        print(f"UNEXPECTED ERROR: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
