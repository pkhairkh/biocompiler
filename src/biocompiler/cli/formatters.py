"""
BioCompiler CLI — Output Formatting
====================================
ANSI colour helpers, progress display, summary boxes, and file I/O utilities.

Extracted from cli.py as part of the SoC refactoring (Wave 4b).
All display/formatting logic lives here so command handlers stay thin.
"""

from __future__ import annotations

import os
import sys
import time
from typing import List

__all__ = [
    # ANSI constants
    "ANSI_RESET",
    "ANSI_BOLD",
    "ANSI_RED",
    "ANSI_GREEN",
    "ANSI_YELLOW",
    "ANSI_CYAN",
    "ANSI_BOLD_RED",
    "ANSI_BOLD_GREEN",
    "ANSI_BOLD_CYAN",
    "ANSI_DIM",
    # Colour helpers
    "colorize",
    "supports_color",
    "section_header",
    "verdict_symbol",
    "error_msg",
    "success_msg",
    "dim",
    "summary_box",
    # Progress
    "ProgressStep",
    "_ProgressStep",
    # File I/O
    "read_fasta",
    "_read_fasta",
    "write_fasta",
    "_write_fasta",
    "write_certificate",
    "_write_certificate",
    # Table formatting
    "print_structure_quality",
    "print_mutation_table",
]

# ── ANSI colour constants ────────────────────────────────────────────────────

ANSI_RESET = "\033[0m"
ANSI_BOLD = "\033[1m"
ANSI_RED = "\033[31m"
ANSI_GREEN = "\033[32m"
ANSI_YELLOW = "\033[33m"
ANSI_CYAN = "\033[36m"
ANSI_BOLD_RED = "\033[1;31m"
ANSI_BOLD_GREEN = "\033[1;32m"
ANSI_BOLD_CYAN = "\033[1;36m"
ANSI_DIM = "\033[2m"


# ── Colour helper functions ──────────────────────────────────────────────────

def supports_color() -> bool:
    """Return True if stdout is a TTY that likely supports ANSI colours."""
    return hasattr(sys.stdout, "isatty") and sys.stdout.isatty()


def colorize(text: str, *codes: str) -> str:
    """Wrap *text* in ANSI escape codes; no-op when stdout is not a TTY."""
    if not supports_color():
        return text
    return "".join(codes) + text + ANSI_RESET


def section_header(text: str) -> str:
    return colorize(text, ANSI_BOLD_CYAN)


def verdict_symbol(value: str) -> str:
    """Return a coloured verdict label."""
    v = value.upper()
    if v in ("PASS", "LIKELY_PASS"):
        return colorize(v, ANSI_BOLD_GREEN)
    if v in ("FAIL", "LIKELY_FAIL"):
        return colorize(v, ANSI_BOLD_RED)
    # UNCERTAIN or anything else
    return colorize(v, ANSI_YELLOW)


def error_msg(text: str) -> str:
    return colorize(text, ANSI_RED)


def success_msg(text: str) -> str:
    return colorize(text, ANSI_BOLD_GREEN)


def dim(text: str) -> str:
    return colorize(text, ANSI_DIM)


def summary_box(label: str, value: str) -> str:
    """Build a Unicode box around a label + value pair."""
    inner = f"{label}: {value}"
    width = len(inner) + 2  # 1 space padding each side
    top = "\u250c" + "\u2500" * width + "\u2510"
    mid = "\u2502 " + inner + " \u2502"
    bot = "\u2514" + "\u2500" * width + "\u2518"
    return "\n".join([top, mid, bot])


# ── Progress helper ──────────────────────────────────────────────────────────

class ProgressStep:
    """Context manager that prints a step label to stderr and appends timing on exit."""

    def __init__(self, label: str, verbose: bool = False) -> None:
        self.label = label
        self.verbose = verbose
        self._t0: float = 0.0

    def __enter__(self) -> "ProgressStep":
        self._t0 = time.perf_counter()
        sys.stderr.write(f"{self.label}...")
        sys.stderr.flush()
        return self

    def __exit__(self, *exc: object) -> None:
        elapsed = time.perf_counter() - self._t0
        timing = f" ({elapsed:.3f}s)" if self.verbose else ""
        sys.stderr.write(f" done{timing}\n")
        sys.stderr.flush()


# Backward compat alias (old cli.py used _ProgressStep)
_ProgressStep = ProgressStep


# ── File I/O helpers ─────────────────────────────────────────────────────────

def read_fasta(path: str) -> str:
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


# Backward compat alias (old cli.py used _read_fasta)
_read_fasta = read_fasta


# Backward compat aliases for old cli.py internal names
_section_header = section_header
_verdict_symbol = verdict_symbol
_error_msg = error_msg
_success_msg = success_msg
_dim = dim
_summary_box = summary_box
_supports_color = supports_color


def write_fasta(path: str, seq: str, header: str = "optimized") -> None:
    """Write a DNA sequence to a FASTA file with 80-char line wrapping."""
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80] + "\n")


def write_certificate(path: str, cert_text: str) -> None:
    """Write certificate text to a file."""
    with open(path, "w") as f:
        f.write(cert_text)


# Backward compat aliases (old cli.py used underscore-prefixed names)
_write_fasta = write_fasta
_write_certificate = write_certificate


# ── Table / report formatting ────────────────────────────────────────────────

def print_structure_quality(report: object) -> None:
    """Pretty-print a StructureQualityReport."""
    print(section_header("  Quality Metrics"))
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
    print(summary_box("Structure Verdict", verdict_symbol(verdict)))


def print_mutation_table(mut_results: List[object]) -> None:
    """Print a table of mutation scan results."""
    print()
    print(section_header("  Mutation Scan Results"))
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
            success_msg(effect) if effect == "stabilizing"
            else error_msg(effect) if effect == "destabilizing"
            else effect
        )
        print(f"  {pos:>8}  {orig:>8}  {mutant:>8}  {ddg_str:>14}  {effect_colored}")
