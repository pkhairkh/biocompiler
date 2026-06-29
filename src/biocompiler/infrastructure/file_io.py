"""
BioCompiler Infrastructure — Common File I/O Operations
=======================================================

Pure infrastructure functions for file path resolution and reading.
Separated from domain logic (parsing FASTA/GenBank) to enforce SoC:
the domain layer should not know about filesystem paths or file handles.

Functions:
  - resolve_input  : Read a file if the input looks like a path, else return as text.
  - looks_like_path: Heuristic to determine if a string is a file path vs sequence data.

These were originally embedded in ``sequence.import_seq`` but are promoted
here so that other infrastructure modules (exporters, LIMS, etc.) can
reuse them without importing from the domain layer.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Optional

from ..shared.exceptions import FileFormatError

logger = logging.getLogger(__name__)

__all__ = [
    "resolve_input",
    "looks_like_path",
]


def resolve_input(filepath_or_text: str, format_name: str = "auto-detect") -> str:
    """Resolve input: if it looks like a file path and the file exists, read it;
    otherwise treat as text content.

    This is a pure infrastructure operation — it handles the I/O boundary
    between the filesystem and the domain parser. Domain parsers should
    receive pre-resolved text strings, not raw file paths.

    Args:
        filepath_or_text: Either a path to a file, or raw text content.
        format_name: Name of the format (used in error messages).

    Returns:
        The text content (either read from the file, or the original input).

    Raises:
        FileFormatError: If the path exists but the file cannot be read.
    """
    if looks_like_path(filepath_or_text):
        path = Path(filepath_or_text)
        if path.exists():
            try:
                return path.read_text()
            except (OSError, IOError) as e:
                logger.error("Failed to read file %s: %s", filepath_or_text, e)
                raise FileFormatError(filepath_or_text, format_name, f"Cannot read file: {e}") from e
        else:
            logger.debug("Path-like input %r does not exist, treating as text content", filepath_or_text)
    # Treat as raw text content
    return filepath_or_text


def looks_like_path(text: str) -> bool:
    """Heuristic: check if the text looks like a file path rather than sequence data.

    Uses several signals:
    - Absolute or relative path prefixes (``/``, ``./``, ``..``)
    - Windows drive letters (``C:\\``)
    - Common bioinformatics file extensions (``.fasta``, ``.gb``, etc.)
    - Single-line strings with dots that look like filenames

    Args:
        text: The string to check.

    Returns:
        True if the text is likely a file path, False otherwise.
    """
    if not text:
        return False
    # Check for path-like patterns
    if "/" in text or "\\" in text:
        # But FASTA headers can have | — check if it starts with a path-like pattern
        if text.startswith("/") or text.startswith("./") or text.startswith(".."):
            return True
        # Windows path
        if re.match(r'^[A-Za-z]:[\\\/]', text):
            return True
    # If it ends with a common file extension
    if re.search(r'\.(fasta|fa|fna|gb|gbk|genbank|txt|json|xml|csv|tsv)$', text, re.IGNORECASE):
        return True
    # If it is a single line with no whitespace and looks like a path
    if "\n" not in text and "." in text and not text.startswith(">") and not text.startswith("LOCUS"):
        # Could be a path like "gene.gb" or "data/gene.fasta"
        if re.match(r'^[\w./\\:_-]+$', text):
            return True
    return False
