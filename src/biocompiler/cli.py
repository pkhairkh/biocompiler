"""
BioCompiler CLI — Command-Line Interface

Production-grade CLI with:
- File input support (FASTA, plain text) for sequences and proteins
- Sequence length validation
- Enzyme name validation
- Exon boundary validation
- All verification parameters embedded in certificates

Usage:
    biocompiler check --input-file gene.fasta --exons 0,92 273,495 1346,1608
    biocompiler check --sequence ATGGTGCATCTG... --exons 0,92 273,495 1346,1608
    biocompiler optimize --protein-file protein.fasta --organism Homo_sapiens
    biocompiler optimize --protein MVHLTPEEK... --organism Homo_sapiens
    biocompiler verify --certificate certificate.json
    biocompiler scan --input-file gene.fasta
    biocompiler scan --sequence ATGGTGCATCTG...
"""

import argparse
import json
import sys
import logging
from pathlib import Path
from . import __version__
from .scanner import scan_sequence, gc_content, validate_dna_sequence
from .translation import translate, compute_cai, find_orfs
from .type_system import evaluate_all_predicates
from .certificate import generate_certificate, verify_certificate
from .optimization import optimize_sequence
from .types import Verdict, combined_verdict
from .constants import RESTRICTION_ENZYMES
from .exceptions import BioCompilerError, CertificateGenerationError, InvalidSequenceError

logger = logging.getLogger("biocompiler")

# Maximum sequence length to prevent denial-of-service
MAX_SEQUENCE_LENGTH = 10_000_000  # 10 Mbp


def _setup_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _read_sequence_file(path: str) -> str:
    """
    Read a DNA sequence from a file (FASTA or plain text).

    Args:
        path: path to the input file

    Returns:
        DNA sequence as a string (no FASTA header, no whitespace)

    Raises:
        argparse.ArgumentTypeError: if file cannot be read
    """
    try:
        text = Path(path).read_text()
    except (OSError, IOError) as e:
        raise argparse.ArgumentTypeError(f"Cannot read file '{path}': {e}")

    # Handle FASTA format: skip header lines starting with '>'
    lines = text.strip().splitlines()
    if lines and lines[0].startswith(">"):
        lines = lines[1:]

    # Remove whitespace and join
    sequence = "".join(line.strip() for line in lines)
    return sequence


def _read_protein_file(path: str) -> str:
    """
    Read a protein sequence from a file (FASTA or plain text).

    Args:
        path: path to the input file

    Returns:
        Protein sequence as a string (single-letter codes, no whitespace)
    """
    try:
        text = Path(path).read_text()
    except (OSError, IOError) as e:
        raise argparse.ArgumentTypeError(f"Cannot read file '{path}': {e}")

    lines = text.strip().splitlines()
    if lines and lines[0].startswith(">"):
        lines = lines[1:]

    protein = "".join(line.strip() for line in lines)
    return protein.upper()


def _validate_sequence_length(seq: str, max_len: int = MAX_SEQUENCE_LENGTH) -> str:
    """Validate that a sequence is not excessively long."""
    if len(seq) > max_len:
        raise argparse.ArgumentTypeError(
            f"Sequence length {len(seq)} exceeds maximum {max_len}. "
            f"Use --max-length to override if needed."
        )
    return seq


def _validate_enzyme_names(enzymes_str: str) -> list[str]:
    """Validate that enzyme names are recognized. Warn about unknown ones."""
    names = [e.strip() for e in enzymes_str.split(",")]
    valid = []
    for name in names:
        if name in RESTRICTION_ENZYMES:
            valid.append(name)
        else:
            logger.warning("Unknown enzyme '%s' — skipping. Known: %s",
                           name, list(RESTRICTION_ENZYMES.keys()))
    return valid


def _parse_exons(exons_str: str | None, seq_len: int = 0) -> list[tuple[int, int]]:
    """
    Parse exon boundaries from comma-separated 'start,end start,end ...' format.

    Validates that boundaries are within sequence length and start < end.
    """
    if not exons_str:
        return []
    boundaries = []
    for pair in exons_str.split():
        parts = pair.split(",")
        if len(parts) != 2:
            raise argparse.ArgumentTypeError(
                f"Invalid exon boundary format '{pair}'. Expected 'start,end'."
            )
        try:
            start, end = int(parts[0]), int(parts[1])
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Invalid exon boundary values '{pair}'. Must be integers."
            )
        if start < 0:
            raise argparse.ArgumentTypeError(
                f"Exon start must be non-negative, got {start}."
            )
        if end <= start:
            raise argparse.ArgumentTypeError(
                f"Exon end must be greater than start, got start={start}, end={end}."
            )
        if seq_len > 0 and end > seq_len:
            raise argparse.ArgumentTypeError(
                f"Exon end {end} exceeds sequence length {seq_len}."
            )
        boundaries.append((start, end))
    return boundaries


def _get_sequence(args, field_name: str = "sequence") -> str:
    """Get sequence from --sequence or --input-file argument."""
    seq = getattr(args, field_name, None)
    if seq is None and hasattr(args, "input_file") and args.input_file:
        seq = _read_sequence_file(args.input_file)
    elif seq is None and hasattr(args, field_name + "_file"):
        file_attr = getattr(args, field_name + "_file", None)
        if file_attr:
            seq = _read_sequence_file(file_attr)
    if not seq:
        raise argparse.ArgumentTypeError(
            f"Must provide --{field_name.replace('_', '-')} or --input-file"
        )
    return seq.upper()


def cmd_check(args):
    """Run type checking on a DNA sequence."""
    _setup_logging(args.verbose)

    # Get sequence from file or argument
    if args.input_file:
        seq = _read_sequence_file(args.input_file).upper()
    elif args.sequence:
        seq = args.sequence.upper()
    else:
        print("ERROR: Must provide --sequence or --input-file", file=sys.stderr)
        sys.exit(1)

    _validate_sequence_length(seq, args.max_length)

    # Validate DNA
    try:
        seq = validate_dna_sequence(seq)
    except InvalidSequenceError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # Parse and validate exon boundaries
    exon_boundaries = _parse_exons(args.exons, len(seq))

    print(f"Sequence length: {len(seq)} nt")
    print(f"GC content: {gc_content(seq):.1%}")
    print(f"Exon boundaries: {exon_boundaries}")
    print()

    # Validate enzyme names
    enzymes = _validate_enzyme_names(args.enzymes) if args.enzymes else None

    # Scan
    tokens = scan_sequence(seq, enzymes)
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
                    "exon_boundaries": exon_boundaries or [(0, len(seq))],
                    "gc_lo": args.gc_lo,
                    "gc_hi": args.gc_hi,
                    "cai_threshold": args.cai_threshold,
                    "enzymes": enzymes or list(RESTRICTION_ENZYMES.keys()),
                    "cell_type": "HEK293T",
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

    # Get protein from file or argument
    if args.protein_file:
        protein = _read_protein_file(args.protein_file)
    elif args.protein:
        protein = args.protein.upper()
    else:
        print("ERROR: Must provide --protein or --protein-file", file=sys.stderr)
        sys.exit(1)

    # Validate enzyme names if provided
    enzymes = _validate_enzyme_names(args.enzymes) if args.enzymes else None

    print(f"Optimizing for protein ({len(protein)} aa), organism={args.organism}")

    result = optimize_sequence(
        target_protein=protein,
        organism=args.organism,
        gc_lo=args.gc_lo,
        gc_hi=args.gc_hi,
        cai_threshold=args.cai_threshold,
        restriction_sites=enzymes,
        cryptic_splice_threshold=args.cryptic_splice_threshold,
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

    # Get sequence from file or argument
    if args.input_file:
        seq = _read_sequence_file(args.input_file).upper()
    elif args.sequence:
        seq = args.sequence.upper()
    else:
        print("ERROR: Must provide --sequence or --input-file", file=sys.stderr)
        sys.exit(1)

    # Validate DNA
    try:
        seq = validate_dna_sequence(seq)
    except InvalidSequenceError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    enzymes = _validate_enzyme_names(args.enzymes) if args.enzymes else None
    tokens = scan_sequence(seq, enzymes)

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


def _add_sequence_args(parser):
    """Add common sequence input arguments to a subparser."""
    seq_group = parser.add_mutually_exclusive_group()
    seq_group.add_argument("--sequence", help="DNA sequence (use --input-file for long sequences)")
    seq_group.add_argument("--input-file", "-f", help="Input file (FASTA or plain text)")
    parser.add_argument("--max-length", type=int, default=MAX_SEQUENCE_LENGTH,
                        help="Maximum sequence length in bp (default: 10M)")


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
    _add_sequence_args(p_check)
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
    prot_group = p_opt.add_mutually_exclusive_group(required=True)
    prot_group.add_argument("--protein", help="Target protein sequence")
    prot_group.add_argument("--protein-file", help="Input file with protein sequence (FASTA or plain text)")
    p_opt.add_argument("--organism", default="Homo_sapiens", help="Target organism")
    p_opt.add_argument("--gc-lo", type=float, default=0.30)
    p_opt.add_argument("--gc-hi", type=float, default=0.70)
    p_opt.add_argument("--cai-threshold", type=float, default=0.2, help="Minimum CAI threshold")
    p_opt.add_argument("--enzymes", help="Comma-separated restriction enzymes to avoid")
    p_opt.add_argument("--cryptic-splice-threshold", type=float, default=3.0,
                        help="MaxEntScan threshold for cryptic splice sites")
    p_opt.add_argument("--output", "-o", help="Output result file path")
    p_opt.set_defaults(func=cmd_optimize)

    # verify
    p_verify = subparsers.add_parser("verify", help="Verify a certificate")
    p_verify.add_argument("--certificate", required=True, help="Certificate JSON file")
    p_verify.set_defaults(func=cmd_verify)

    # scan
    p_scan = subparsers.add_parser("scan", help="Scan a sequence for biological motifs")
    _add_sequence_args(p_scan)
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
    except argparse.ArgumentTypeError as e:
        print(f"ARGUMENT ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        logger.exception("Unexpected error")
        print(f"UNEXPECTED ERROR: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
