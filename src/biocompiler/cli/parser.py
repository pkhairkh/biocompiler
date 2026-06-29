"""
BioCompiler CLI — Argument Parser
===================================
All argparse definitions, subcommands, and argument defaults.

Extracted from cli.py as part of the SoC refactoring (Wave 4b).
Only argument parsing lives here — no business logic or output formatting.
"""

from __future__ import annotations

import argparse

from .. import __version__

__all__ = [
    "build_parser",
    "add_protein_args",
]


def add_protein_args(parser: argparse.ArgumentParser) -> None:
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

    # ── optimize ── (v1.0.0: positional PROTEIN arg + --input fallback)
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
    opt_parser.add_argument(
        "--biosecurity-report",
        action="store_true",
        default=False,
        help=(
            "Print a full biosecurity screening report before outputting "
            "the optimized sequence. The report includes biosafety level "
            "(BSL-1/BSL-2), screening status, predicate results, and "
            "risk assessment summary."
        ),
    )

    # ── batch ── (v1.0.0)
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
        help="Check a FASTA gene sequence against all registered predicates",
    )
    check_parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="Input FASTA file path (not needed with --list-predicates)",
    )
    check_parser.add_argument(
        "--species", default="human",
        help="Target species for codon evaluation (default: human). "
             "Accepts any organism supported by resolve_organism (e.g., ecoli, human, CHO_K1)",
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
    check_parser.add_argument(
        "--predicate", default=None, metavar="NAMES",
        help="Comma-separated predicate names to check (default: all). "
             "Use --list-predicates to see available names",
    )
    check_parser.add_argument(
        "--list-predicates", action="store_true", default=False,
        help="List all available predicate names and exit",
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
    serve_parser.add_argument(
        "--no-auth", action="store_true",
        help="Disable API authentication (for local development only! DANGEROUS in production)",
    )

    # ── structure ──
    struct_parser = subparsers.add_parser(
        "structure",
        help="Predict and assess protein structure",
        description="Predict protein 3D structure using ESMFold and assess quality "
                    "(pLDDT, Ramachandran, clash score).",
    )
    add_protein_args(struct_parser)
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
    add_protein_args(stab_parser)
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
    add_protein_args(sol_parser)
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
    add_protein_args(imm_parser)
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
    add_protein_args(assess_parser)
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
    add_protein_args(whatif_parser)
    whatif_parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show timing information",
    )

    return parser
