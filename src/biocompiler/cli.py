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
import time
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

# ==============================================================================
# ANSI Color Support
# ==============================================================================

# ANSI escape codes
_ANSI_RESET = "\033[0m"
_ANSI_BOLD = "\033[1m"
_ANSI_RED = "\033[31m"
_ANSI_GREEN = "\033[32m"
_ANSI_YELLOW = "\033[33m"
_ANSI_CYAN = "\033[36m"
_ANSI_BOLD_RED = "\033[1;31m"
_ANSI_BOLD_GREEN = "\033[1;32m"
_ANSI_BOLD_CYAN = "\033[1;36m"


def _use_color():
    """Check whether ANSI color codes should be emitted."""
    return sys.stdout.isatty()


def colorize(text: str, *codes: str) -> str:
    """Wrap *text* in ANSI escape codes if stdout is a TTY."""
    if not _use_color():
        return text
    return "".join(codes) + text + _ANSI_RESET


def _verdict_symbol(verdict_value: str) -> str:
    """Return a colored verdict symbol."""
    if verdict_value == "PASS":
        return colorize("PASS", _ANSI_BOLD_GREEN)
    elif verdict_value == "FAIL":
        return colorize("FAIL", _ANSI_BOLD_RED)
    else:  # UNCERTAIN
        return colorize("UNCERTAIN", _ANSI_YELLOW)


def _section_header(text: str) -> str:
    """Return a colored section header."""
    return colorize(text, _ANSI_BOLD_CYAN)


def _error_msg(text: str) -> str:
    """Return a colored error message."""
    return colorize(text, _ANSI_RED)


def _success_msg(text: str) -> str:
    """Return a colored success message."""
    return colorize(text, _ANSI_GREEN)


# ==============================================================================
# Progress Indicators (written to stderr)
# ==============================================================================

class _ProgressPhase:
    """Context manager that prints a progress label to stderr and clears it on exit."""

    def __init__(self, label: str, verbose: bool = False):
        self.label = label
        self.verbose = verbose
        self._start: float = 0.0

    def __enter__(self):
        self._start = time.perf_counter()
        print(f"{self.label}", end="", file=sys.stderr, flush=True)
        return self

    def __exit__(self, *exc):
        elapsed = time.perf_counter() - self._start
        if self.verbose:
            print(f" done ({elapsed:.3f}s)", file=sys.stderr, flush=True)
        else:
            print(f" done", file=sys.stderr, flush=True)
        return False


def _progress_dot_loop(stop_event, interval: float = 1.0):
    """Print dots to stderr at *interval* seconds until *stop_event* is set."""
    import threading
    def _dots():
        while not stop_event.is_set():
            stop_event.wait(interval)
            if not stop_event.is_set():
                print(".", end="", file=sys.stderr, flush=True)
    t = threading.Thread(target=_dots, daemon=True)
    t.start()
    return t


# ==============================================================================
# Summary Box
# ==============================================================================

def _summary_box(label: str, value: str) -> str:
    """Build a Unicode box around a key–value summary line.

    Example::

        ┌─────────────────────────────┐
        │ Overall Verdict: PASS       │
        └─────────────────────────────┘
    """
    content = f" {label}: {value} "
    width = max(len(content), 31)  # minimum visual width
    # Pad content to fill the box
    content = content.ljust(width)
    top = "┌" + "─" * width + "┐"
    mid = "│" + content + "│"
    bot = "└" + "─" * width + "┘"
    return f"{top}\n{mid}\n{bot}"

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
    verbose = args.verbose

    # Get sequence from file or argument
    if args.input_file:
        seq = _read_sequence_file(args.input_file).upper()
    elif args.sequence:
        seq = args.sequence.upper()
    else:
        print(_error_msg("ERROR: Must provide --sequence or --input-file"), file=sys.stderr)
        sys.exit(1)

    _validate_sequence_length(seq, args.max_length)

    # Validate DNA
    try:
        seq = validate_dna_sequence(seq)
    except InvalidSequenceError as e:
        print(_error_msg(f"ERROR: {e}"), file=sys.stderr)
        sys.exit(1)

    # Parse and validate exon boundaries
    exon_boundaries = _parse_exons(args.exons, len(seq))

    print(_section_header("Sequence Info"))
    print(f"  Sequence length: {len(seq)} nt")
    print(f"  GC content: {gc_content(seq):.1%}")
    print(f"  Exon boundaries: {exon_boundaries}")
    print()

    # Validate enzyme names
    enzymes = _validate_enzyme_names(args.enzymes) if args.enzymes else None

    # Scan phase
    with _ProgressPhase("Scanning", verbose=verbose):
        t0 = time.perf_counter()
        tokens = scan_sequence(seq, enzymes)
        scan_time = time.perf_counter() - t0

    print(_section_header("Scan Results"))
    print(f"  Scanner found {len(tokens)} tokens")
    donor_count = sum(1 for t in tokens if t.element_type == "splice_donor")
    acceptor_count = sum(1 for t in tokens if t.element_type == "splice_acceptor")
    print(f"  Donor sites: {donor_count}")
    print(f"  Acceptor sites: {acceptor_count}")
    if verbose:
        print(f"  [timing] scan: {scan_time:.3f}s")
    print()

    # Translate phase
    with _ProgressPhase("Type-checking", verbose=verbose):
        t0 = time.perf_counter()
        if exon_boundaries:
            coding_seq = "".join(seq[start:end] for start, end in exon_boundaries)
        else:
            coding_seq = seq
        protein = translate(coding_seq)
        translate_time = time.perf_counter() - t0

    print(_section_header("Translation"))
    print(f"  Protein: {protein[:50]}{'...' if len(protein) > 50 else ''} ({len(protein)} aa)")
    if verbose:
        print(f"  [timing] translate: {translate_time:.3f}s")
    print()

    # Type check phase
    with _ProgressPhase("Evaluating predicates", verbose=verbose):
        t0 = time.perf_counter()
        results = evaluate_all_predicates(
            seq=seq,
            known_exon_boundaries=exon_boundaries or [(0, len(seq))],
            organism=args.organism,
            gc_lo=args.gc_lo,
            gc_hi=args.gc_hi,
            cai_threshold=args.cai_threshold,
        )
        eval_time = time.perf_counter() - t0

    print(_section_header("Type-Check Results"))
    for r in results:
        symbol = _verdict_symbol(r.verdict.value)
        print(f"  [{symbol}] {r.predicate}")
        if r.violation:
            print(f"         {r.violation}")
    if verbose:
        print(f"  [timing] predicates: {eval_time:.3f}s")
    print()

    overall = combined_verdict([r.verdict for r in results])

    # Summary box
    verdict_display = _verdict_symbol(overall.value)
    box = _summary_box("Overall Verdict", verdict_display)
    print(box)

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
            print(_success_msg(f"Certificate saved to: {args.output}"))
        except CertificateGenerationError as e:
            print(_error_msg(f"Certificate generation failed: {e}"))


def cmd_optimize(args):
    """Optimize a DNA sequence for a target protein."""
    _setup_logging(args.verbose)
    verbose = args.verbose

    # Get protein from file or argument
    if args.protein_file:
        protein = _read_protein_file(args.protein_file)
    elif args.protein:
        protein = args.protein.upper()
    else:
        print(_error_msg("ERROR: Must provide --protein or --protein-file"), file=sys.stderr)
        sys.exit(1)

    # Validate enzyme names if provided
    enzymes = _validate_enzyme_names(args.enzymes) if args.enzymes else None

    print(_section_header("Optimization"))
    print(f"  Optimizing for protein ({len(protein)} aa), organism={args.organism}")

    # Progress dots while optimizing
    import threading
    stop = threading.Event()
    dot_thread = _progress_dot_loop(stop, interval=1.0)
    t0 = time.perf_counter()

    try:
        result = optimize_sequence(
            target_protein=protein,
            organism=args.organism,
            gc_lo=args.gc_lo,
            gc_hi=args.gc_hi,
            cai_threshold=args.cai_threshold,
            restriction_sites=enzymes,
            cryptic_splice_threshold=args.cryptic_splice_threshold,
        )
    finally:
        stop.set()
        dot_thread.join(timeout=2)

    opt_time = time.perf_counter() - t0
    print("", file=sys.stderr)  # newline after dots

    print(f"  Sequence: {result.sequence[:60]}{'...' if len(result.sequence) > 60 else ''}")
    print(f"  Length: {len(result.sequence)} bp")
    print(f"  CAI: {result.cai:.4f}")
    print(f"  GC: {result.gc_content:.1%}")
    print(f"  Satisfied predicates: {result.satisfied_predicates}")
    print(f"  Failed predicates: {result.failed_predicates}")
    print(f"  Fallback solver used: {result.fallback_used}")
    if verbose:
        print(f"  [timing] optimization: {opt_time:.3f}s")
    print()

    if result.failed_predicates == 0:
        print(_success_msg("Optimization completed — all predicates satisfied."))
    else:
        print(_error_msg(f"Optimization completed with {result.failed_predicates} failed predicate(s)."))

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
        print(_success_msg(f"Result saved to: {args.output}"))


def cmd_verify(args):
    """Verify a certificate."""
    _setup_logging(args.verbose)
    with open(args.certificate) as f:
        cert_dict = json.load(f)

    status, failures = verify_certificate(cert_dict)
    if status == "VERIFIED":
        print(_success_msg(f"Verification status: {status}"))
    else:
        print(_error_msg(f"Verification status: {status}"))
    if failures:
        for f in failures:
            print(f"  {_error_msg('FAILURE:')} {f}")
    else:
        print(_success_msg("  All checks passed."))


def cmd_scan(args):
    """Scan a sequence for motifs."""
    _setup_logging(args.verbose)

    # Get sequence from file or argument
    if args.input_file:
        seq = _read_sequence_file(args.input_file).upper()
    elif args.sequence:
        seq = args.sequence.upper()
    else:
        print(_error_msg("ERROR: Must provide --sequence or --input-file"), file=sys.stderr)
        sys.exit(1)

    # Validate DNA
    try:
        seq = validate_dna_sequence(seq)
    except InvalidSequenceError as e:
        print(_error_msg(f"ERROR: {e}"), file=sys.stderr)
        sys.exit(1)

    enzymes = _validate_enzyme_names(args.enzymes) if args.enzymes else None
    with _ProgressPhase("Scanning", verbose=args.verbose):
        tokens = scan_sequence(seq, enzymes)

    print(_section_header("Scan Results"))
    print(f"  Sequence length: {len(seq)} nt")
    print(f"  Tokens found: {len(tokens)}")
    for t in tokens:
        frame_str = f" frame={t.frame}" if t.frame is not None else ""
        strand_str = f" strand={t.strand}" if t.strand != "+" else ""
        print(f"    {t.element_type} @ {t.position}: {t.match_sequence} (score={t.score:.1f}){frame_str}{strand_str}")

    if args.find_orfs:
        orfs = find_orfs(seq)
        print(f"\n  ORFs found: {len(orfs)}")
        for orf in orfs:
            print(f"    {orf['strand']} strand, frame {orf['frame']}: "
                  f"{orf['start']}-{orf['end']} ({orf['length']} aa)")


def _add_sequence_args(parser):
    """Add common sequence input arguments to a subparser."""
    seq_group = parser.add_mutually_exclusive_group()
    seq_group.add_argument("--sequence", help="DNA sequence (use --input-file for long sequences)")
    seq_group.add_argument("--input-file", "-f", help="Input file (FASTA or plain text)")
    parser.add_argument("--max-length", type=int, default=MAX_SEQUENCE_LENGTH,
                        help="Maximum sequence length in bp (default: 10M)")


def cmd_export(args):
    """Export a sequence in FASTA or GenBank format."""
    _setup_logging(args.verbose)

    # Get sequence
    if args.input_file:
        seq = _read_sequence_file(args.input_file).upper()
    elif args.sequence:
        seq = args.sequence.upper()
    else:
        print(_error_msg("ERROR: Must provide --sequence or --input-file"), file=sys.stderr)
        sys.exit(1)

    try:
        seq = validate_dna_sequence(seq)
    except InvalidSequenceError as e:
        print(_error_msg(f"ERROR: {e}"), file=sys.stderr)
        sys.exit(1)

    exon_boundaries = _parse_exons(args.exons, len(seq)) if args.exons else None

    if args.format == "fasta":
        from .export import export_fasta
        output = export_fasta(
            sequence=seq,
            identifier=args.identifier or "BioCompiler_design",
            description=args.description or "",
            organism=args.organism,
        )
    elif args.format == "genbank":
        from .export import export_genbank
        output = export_genbank(
            sequence=seq,
            locus_name=args.locus or "BIOCOMPILER",
            definition=args.description or "BioCompiler designed sequence",
            organism=args.organism,
            exon_boundaries=exon_boundaries,
            gene_name=args.gene,
        )
    else:
        print(_error_msg(f"ERROR: Unknown format '{args.format}'. Use 'fasta' or 'genbank'."), file=sys.stderr)
        sys.exit(1)

    if args.output:
        Path(args.output).write_text(output)
        print(_success_msg(f"Exported to: {args.output}"))
    else:
        print(output)


def cmd_report(args):
    """Generate an interactive HTML report."""
    _setup_logging(args.verbose)
    verbose = args.verbose

    if args.input_file:
        seq = _read_sequence_file(args.input_file).upper()
    elif args.sequence:
        seq = args.sequence.upper()
    else:
        print(_error_msg("ERROR: Must provide --sequence or --input-file"), file=sys.stderr)
        sys.exit(1)

    try:
        seq = validate_dna_sequence(seq)
    except InvalidSequenceError as e:
        print(_error_msg(f"ERROR: {e}"), file=sys.stderr)
        sys.exit(1)

    exon_boundaries = _parse_exons(args.exons, len(seq)) or [(0, len(seq))]

    from .report import generate_report
    from .type_system import evaluate_all_predicates

    with _ProgressPhase("Generating report", verbose=verbose):
        t0 = time.perf_counter()
        results = evaluate_all_predicates(
            seq=seq,
            known_exon_boundaries=exon_boundaries,
            organism=args.organism,
        )

        html_report = generate_report(
            sequence=seq,
            type_results=results,
            organism=args.organism,
            gene_name=args.gene,
            exon_boundaries=exon_boundaries if len(exon_boundaries) > 1 else None,
        )
        report_time = time.perf_counter() - t0

    output_path = args.output or "biocompiler_report.html"
    Path(output_path).write_text(html_report)
    print(_success_msg(f"Report generated: {output_path}"))
    if verbose:
        print(f"  [timing] report generation: {report_time:.3f}s")


def cmd_benchmark(args):
    """Run benchmarks against known gene sets."""
    _setup_logging(args.verbose)
    verbose = args.verbose

    from .benchmark import run_benchmarks, format_benchmark_report_text, format_benchmark_report_json, REFERENCE_GENES

    gene_names = args.genes.split(",") if args.genes else None
    all_genes = gene_names or list(REFERENCE_GENES.keys())
    total = len(all_genes)

    print(_section_header("Benchmark"))
    print(f"  Running {total} gene(s): {', '.join(all_genes)}", file=sys.stderr)

    # Run benchmarks with per-gene progress
    from .benchmark import _bench_translation, _bench_gc_content, _bench_cai
    from .benchmark import _bench_splice_isoforms, _bench_type_predicates
    from .benchmark import _bench_certificate_roundtrip, _bench_optimization
    from .benchmark import _compute_summary, BenchmarkReport
    from . import __version__ as pkg_version
    from datetime import datetime, timezone

    results = []
    include_optimization = not args.skip_optimization

    for idx, gene_name in enumerate(all_genes, start=1):
        print(f"  [{idx}/{total}] {gene_name}...", end="", file=sys.stderr, flush=True)
        gene_data = REFERENCE_GENES.get(gene_name)
        if not gene_data:
            print(" skipped (unknown)", file=sys.stderr)
            continue

        seq = gene_data["pre_mrna"].replace(" ", "")
        exons = gene_data["exon_boundaries"]
        organism = gene_data["organism"]

        results.append(_bench_translation(gene_name, seq, exons, gene_data))
        results.append(_bench_gc_content(gene_name, seq, gene_data))
        results.append(_bench_cai(gene_name, seq, exons, gene_data))
        if len(exons) > 1:
            results.append(_bench_splice_isoforms(gene_name, seq, exons, gene_data))
        results.append(_bench_type_predicates(gene_name, seq, exons, organism))
        results.append(_bench_certificate_roundtrip(gene_name, seq, exons, organism))
        if include_optimization:
            results.append(_bench_optimization(gene_name, seq, exons, organism))

        print(" done", file=sys.stderr)

    total_tests = len(results)
    passed = sum(1 for r in results if r.passed)
    failed = total_tests - passed

    report = BenchmarkReport(
        timestamp=datetime.now(timezone.utc).isoformat(),
        version=pkg_version,
        total_tests=total_tests,
        passed=passed,
        failed=failed,
        results=results,
        summary=_compute_summary(results),
    )

    if args.format == "json":
        output = format_benchmark_report_json(report)
    else:
        output = format_benchmark_report_text(report)

    if args.output:
        Path(args.output).write_text(output)
        print(_success_msg(f"Benchmark report saved to: {args.output}"))
    else:
        print(output)


def cmd_serve(args):
    """Start the BioCompiler REST API server."""
    _setup_logging(args.verbose)

    from .api import app
    import uvicorn

    print(_section_header(f"Starting BioCompiler API server on {args.host}:{args.port}"))
    print(f"  API docs: http://{args.host}:{args.port}/docs")
    uvicorn.run(app, host=args.host, port=args.port, log_level="info" if not args.verbose else "debug")


def cmd_migrate(args):
    """Migrate built-in organism data to SQLite database."""
    _setup_logging(args.verbose)

    from .organism_db import OrganismDatabase

    with _ProgressPhase("Migrating", verbose=args.verbose):
        db = OrganismDatabase()
        count = db.migrate_builtin_data()

    print(_success_msg(f"Migrated {count} organisms to database at {db.db_path}"))
    organisms = db.list_organisms()
    for org in organisms:
        print(f"  {org['name']} (source={org['source']}, codons=64)")


def cmd_completion(args):
    """Output shell completion script for bash or zsh."""
    shell = args.shell
    from .organisms import SUPPORTED_ORGANISMS
    organisms = " ".join(SUPPORTED_ORGANISMS)
    formats = "fasta genbank"
    enzymes_list = " ".join(sorted(RESTRICTION_ENZYMES.keys()))

    if shell == "bash":
        script = f'''# bash completion for biocompiler
_biocompiler_complete()
{{
    local cur prev words cword
    _init_completion || return

    local subcommands="check optimize verify scan export report benchmark serve migrate completion"

    # If completing the first-level subcommand
    if [[ $cword -eq 1 ]]; then
        COMPREPLY=($(compgen -W "$subcommands" -- "$cur"))
        return
    fi

    local subcommand="${{words[1]}}"

    case $subcommand in
        check)
            case $prev in
                --organism)
                    COMPREPLY=($(compgen -W "{organisms}" -- "$cur"))
                    ;;
                --enzymes)
                    COMPREPLY=($(compgen -W "{enzymes_list}" -- "$cur"))
                    ;;
                --input-file|-f|--output|-o)
                    _filedir
                    ;;
                *)
                    COMPREPLY=($(compgen -W "--sequence --input-file -f --exons --gene --organism --gc-lo --gc-hi --cai-threshold --enzymes --output -o --max-length" -- "$cur"))
                    ;;
            esac
            ;;
        optimize)
            case $prev in
                --organism)
                    COMPREPLY=($(compgen -W "{organisms}" -- "$cur"))
                    ;;
                --protein-file)
                    _filedir
                    ;;
                --output|-o)
                    _filedir
                    ;;
                *)
                    COMPREPLY=($(compgen -W "--protein --protein-file --organism --gc-lo --gc-hi --cai-threshold --enzymes --cryptic-splice-threshold --output -o" -- "$cur"))
                    ;;
            esac
            ;;
        verify)
            case $prev in
                --certificate)
                    _filedir
                    ;;
                *)
                    COMPREPLY=($(compgen -W "--certificate" -- "$cur"))
                    ;;
            esac
            ;;
        scan)
            case $prev in
                --enzymes)
                    COMPREPLY=($(compgen -W "{enzymes_list}" -- "$cur"))
                    ;;
                --input-file|-f)
                    _filedir
                    ;;
                *)
                    COMPREPLY=($(compgen -W "--sequence --input-file -f --enzymes --find-orfs --max-length" -- "$cur"))
                    ;;
            esac
            ;;
        export)
            case $prev in
                --format)
                    COMPREPLY=($(compgen -W "{formats}" -- "$cur"))
                    ;;
                --organism)
                    COMPREPLY=($(compgen -W "{organisms}" -- "$cur"))
                    ;;
                --input-file|-f|--output|-o)
                    _filedir
                    ;;
                *)
                    COMPREPLY=($(compgen -W "--sequence --input-file -f --format --identifier --locus --description --gene --organism --exons --output -o --max-length" -- "$cur"))
                    ;;
            esac
            ;;
        report)
            case $prev in
                --organism)
                    COMPREPLY=($(compgen -W "{organisms}" -- "$cur"))
                    ;;
                --input-file|-f|--output|-o)
                    _filedir
                    ;;
                *)
                    COMPREPLY=($(compgen -W "--sequence --input-file -f --exons --gene --organism --output -o --max-length" -- "$cur"))
                    ;;
            esac
            ;;
        benchmark)
            case $prev in
                --format)
                    COMPREPLY=($(compgen -W "text json" -- "$cur"))
                    ;;
                --output|-o)
                    _filedir
                    ;;
                *)
                    COMPREPLY=($(compgen -W "--genes --format --skip-optimization --output -o" -- "$cur"))
                    ;;
            esac
            ;;
        serve)
            COMPREPLY=($(compgen -W "--host --port" -- "$cur"))
            ;;
        completion)
            COMPREPLY=($(compgen -W "--shell" -- "$cur"))
            ;;
    esac
}}
complete -F _biocompiler_complete biocompiler
'''
    elif shell == "zsh":
        script = f'''#compdef biocompiler
_biocomplier_organisms=({" ".join(f"{o}:\"{o}\"" for o in SUPPORTED_ORGANISMS)})
_biocompiler_formats=(fasta:"FASTA format" genbank:"GenBank format")
_biocompiler_enzymes=({" ".join(f"{e}:\"{e}\"" for e in sorted(RESTRICTION_ENZYMES.keys()))})

_biocompiler_subcommands()
{{
    local -a commands
    commands=(
        'check:Type-check a DNA sequence'
        'optimize:Optimize DNA sequence for a protein'
        'verify:Verify a certificate'
        'scan:Scan a sequence for biological motifs'
        'export:Export sequence in FASTA or GenBank format'
        'report:Generate interactive HTML report'
        'benchmark:Run benchmarks against known gene sets'
        'serve:Start REST API server'
        'migrate:Migrate organism data to SQLite'
        'completion:Output shell completion script'
    )
    _describe 'command' commands
}}

_biocompiler_check()
{{
    _arguments \\
        '--sequence[DNA sequence]' \\
        '--input-file[Input file]:file:_files' \\
        '-f[Input file]:file:_files' \\
        '--exons[Exon boundaries]' \\
        '--gene[Gene name]' \\
        '--organism[Target organism]:organism:($_biocomplier_organisms)' \\
        '--gc-lo[Min GC content]:float' \\
        '--gc-hi[Max GC content]:float' \\
        '--cai-threshold[Min CAI]:float' \\
        '--enzymes[Restriction enzymes]:enzyme:_values "enzymes" $_biocompiler_enzymes' \\
        '--output[Output certificate]:file:_files' \\
        '-o[Output certificate]:file:_files' \\
        '--max-length[Max sequence length]:int'
}}

_biocompiler_optimize()
{{
    _arguments \\
        '--protein[Target protein sequence]' \\
        '--protein-file[Input protein file]:file:_files' \\
        '--organism[Target organism]:organism:($_biocomplier_organisms)' \\
        '--gc-lo[Min GC content]:float' \\
        '--gc-hi[Max GC content]:float' \\
        '--cai-threshold[Min CAI]:float' \\
        '--enzymes[Restriction enzymes]:enzyme:_values "enzymes" $_biocompiler_enzymes' \\
        '--cryptic-splice-threshold[Cryptic splice threshold]:float' \\
        '--output[Output file]:file:_files' \\
        '-o[Output file]:file:_files'
}}

_biocompiler_verify()
{{
    _arguments \\
        '--certificate[Certificate JSON file]:file:_files'
}}

_biocompiler_scan()
{{
    _arguments \\
        '--sequence[DNA sequence]' \\
        '--input-file[Input file]:file:_files' \\
        '-f[Input file]:file:_files' \\
        '--enzymes[Restriction enzymes]:enzyme:_values "enzymes" $_biocompiler_enzymes' \\
        '--find-orfs[Find open reading frames]' \\
        '--max-length[Max sequence length]:int'
}}

_biocompiler_export()
{{
    _arguments \\
        '--sequence[DNA sequence]' \\
        '--input-file[Input file]:file:_files' \\
        '-f[Input file]:file:_files' \\
        '--format[Export format]:format:($_biocompiler_formats)' \\
        '--identifier[Sequence identifier]' \\
        '--locus[LOCUS name]' \\
        '--description[Description]' \\
        '--gene[Gene name]' \\
        '--organism[Source organism]:organism:($_biocomplier_organisms)' \\
        '--exons[Exon boundaries]' \\
        '--output[Output file]:file:_files' \\
        '-o[Output file]:file:_files' \\
        '--max-length[Max sequence length]:int'
}}

_biocompiler_report()
{{
    _arguments \\
        '--sequence[DNA sequence]' \\
        '--input-file[Input file]:file:_files' \\
        '-f[Input file]:file:_files' \\
        '--exons[Exon boundaries]' \\
        '--gene[Gene name]' \\
        '--organism[Target organism]:organism:($_biocomplier_organisms)' \\
        '--output[Output HTML file]:file:_files' \\
        '-o[Output HTML file]:file:_files' \\
        '--max-length[Max sequence length]:int'
}}

_biocompiler_benchmark()
{{
    _arguments \\
        '--genes[Gene names to benchmark]' \\
        '--format[Output format]:format:(text json)' \\
        '--skip-optimization[Skip optimization benchmarks]' \\
        '--output[Output file]:file:_files' \\
        '-o[Output file]:file:_files'
}}

_biocompiler_serve()
{{
    _arguments \\
        '--host[Host to bind]' \\
        '--port[Port to bind]:int'
}}

_biocompiler_completion()
{{
    _arguments \\
        '--shell[Target shell]:shell:(bash zsh)'
}}

_biocompiler_migrate()
{{
    # No arguments
}}

_biocompiler()
{{
    local curcontext="$curcontext" state line
    typeset -A opt_args

    _arguments -C \\
        '1: :_biocompiler_subcommands' \\
        '*::arg:->args'

    case $words[1] in
        check)      _biocompiler_check ;;
        optimize)   _biocompiler_optimize ;;
        verify)     _biocompiler_verify ;;
        scan)       _biocompiler_scan ;;
        export)     _biocompiler_export ;;
        report)     _biocompiler_report ;;
        benchmark)  _biocompiler_benchmark ;;
        serve)      _biocompiler_serve ;;
        migrate)    _biocompiler_migrate ;;
        completion) _biocompiler_completion ;;
    esac
}}

_biocompiler "$@"
'''
    else:
        print(f"ERROR: Unsupported shell '{shell}'. Choose 'bash' or 'zsh'.", file=sys.stderr)
        sys.exit(1)

    print(script)


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

    # export
    p_export = subparsers.add_parser("export", help="Export sequence in FASTA or GenBank format")
    _add_sequence_args(p_export)
    p_export.add_argument("--format", choices=["fasta", "genbank"], default="fasta",
                          help="Export format (default: fasta)")
    p_export.add_argument("--identifier", help="Sequence identifier (FASTA)")
    p_export.add_argument("--locus", help="LOCUS name (GenBank, max 16 chars)")
    p_export.add_argument("--description", help="Description line")
    p_export.add_argument("--gene", help="Gene name")
    p_export.add_argument("--organism", default="Homo_sapiens", help="Source organism")
    p_export.add_argument("--exons", help="Exon boundaries as 'start,end start,end ...'")
    p_export.add_argument("--output", "-o", help="Output file path")
    p_export.set_defaults(func=cmd_export)

    # report
    p_report = subparsers.add_parser("report", help="Generate interactive HTML report")
    _add_sequence_args(p_report)
    p_report.add_argument("--exons", help="Exon boundaries as 'start,end start,end ...'")
    p_report.add_argument("--gene", help="Gene name")
    p_report.add_argument("--organism", default="Homo_sapiens", help="Target organism")
    p_report.add_argument("--output", "-o", help="Output HTML file path")
    p_report.set_defaults(func=cmd_report)

    # benchmark
    p_bench = subparsers.add_parser("benchmark", help="Run benchmarks against known gene sets")
    p_bench.add_argument("--genes", help="Comma-separated gene names to benchmark (default: all)")
    p_bench.add_argument("--format", choices=["text", "json"], default="text", help="Output format")
    p_bench.add_argument("--skip-optimization", action="store_true", help="Skip optimization benchmarks")
    p_bench.add_argument("--output", "-o", help="Output file path")
    p_bench.set_defaults(func=cmd_benchmark)

    # serve
    p_serve = subparsers.add_parser("serve", help="Start REST API server")
    p_serve.add_argument("--host", default="0.0.0.0", help="Host to bind (default: 0.0.0.0)")
    p_serve.add_argument("--port", type=int, default=8000, help="Port to bind (default: 8000)")
    p_serve.set_defaults(func=cmd_serve)

    # migrate
    p_migrate = subparsers.add_parser("migrate", help="Migrate organism data to SQLite")
    p_migrate.set_defaults(func=cmd_migrate)

    # completion
    p_completion = subparsers.add_parser("completion", help="Output shell completion script")
    p_completion.add_argument("--shell", choices=["bash", "zsh"], default="bash",
                              help="Target shell (default: bash)")
    p_completion.set_defaults(func=cmd_completion)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        args.func(args)
    except BioCompilerError as e:
        print(_error_msg(f"ERROR: {e}"), file=sys.stderr)
        sys.exit(1)
    except argparse.ArgumentTypeError as e:
        print(f"ARGUMENT ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        logger.exception("Unexpected error")
        print(_error_msg(f"UNEXPECTED ERROR: {e}"), file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
