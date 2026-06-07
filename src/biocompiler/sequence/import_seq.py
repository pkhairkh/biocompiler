"""
BioCompiler Sequence Import — GenBank & FASTA Parsing

Production-grade sequence import with:
- FASTA format parsing (single and multi-FASTA)
- GenBank format parsing with full feature extraction
- Auto-detection of file format (FASTA, GenBank, plain text)
- Pure Python — no BioPython dependency
- Exon boundary extraction from CDS join() locations
- Protein extraction from /translation qualifier
- Gene name extraction from /gene qualifier
- GC content computation on import
- DNA sequence validation using scanner.validate_dna_sequence

Round-trip capability: import_seq is the inverse of export.
"""

import logging
import re
from pathlib import Path
from typing import Any

from .scanner import validate_dna_sequence, gc_content
from biocompiler.shared.exceptions import FileFormatError, InvalidSequenceError

logger = logging.getLogger(__name__)

__all__ = [
    "import_fasta",
    "import_genbank",
    "import_sequence",
    "_parse_exon_boundaries",
    "_clean_qualifier_value",
    "_looks_like_path",
    "_resolve_input",
]

# GenBank format column constants (RFC 8466 / GenBank flat-file spec)
_GENBANK_FEATURE_INDENT: int = 5    # spaces before feature type
_GENBANK_QUALIFIER_INDENT: int = 21  # spaces before qualifier /key=value


# ─── FASTA Import ──────────────────────────────────────────────────

def import_fasta(filepath_or_text: str) -> list[dict]:
    """
    Parse a FASTA file or text string and return a list of sequence records.

    Each record is a dict with keys:
        id          — sequence identifier (first word after >)
        description — full header line after the identifier
        sequence    — uppercase DNA sequence string
        organism    — organism name extracted from header if present (e.g. organism=X)
        gc_content  — GC fraction computed from sequence

    Handles:
        - Multi-FASTA files (multiple >header/sequence blocks)
        - Windows (\\r\\n) and Unix (\\n) line endings
        - Empty lines between sequence blocks
        - FASTA headers with pipe-delimited metadata

    Args:
        filepath_or_text: Path to a FASTA file, or raw FASTA text content

    Returns:
        List of dicts, one per FASTA record

    Raises:
        FileFormatError: If the input is not valid FASTA or contains invalid DNA
    """
    text = _resolve_input(filepath_or_text, "FASTA")

    # Normalize line endings
    text = text.replace("\r\n", "\n").replace("\r", "\n")

    records = []
    current_header = None
    current_seq_parts: list[str] = []

    for line in text.split("\n"):
        line = line.strip()
        if not line:
            # Skip empty lines
            continue
        if line.startswith(">"):
            # Save previous record if any
            if current_header is not None:
                record = _build_fasta_record(current_header, current_seq_parts)
                records.append(record)
            current_header = line[1:]  # Remove '>'
            current_seq_parts = []
        else:
            if current_header is None:
                # Sequence lines before any header — not valid FASTA
                raise FileFormatError(
                    filepath_or_text if _looks_like_path(filepath_or_text) else "<text>",
                    "FASTA",
                    "Sequence data found before any header line (>)",
                )
            current_seq_parts.append(line)

    # Save last record
    if current_header is not None:
        record = _build_fasta_record(current_header, current_seq_parts)
        records.append(record)

    if not records:
        raise FileFormatError(
            filepath_or_text if _looks_like_path(filepath_or_text) else "<text>",
            "FASTA",
            "No FASTA records found (no '>' header lines)",
        )

    logger.debug("Parsed %d FASTA record(s)", len(records))
    return records


def _build_fasta_record(header: str, seq_parts: list[str]) -> dict[str, Any]:
    """Build a FASTA record dict from a header line and sequence parts."""
    sequence_raw = "".join(seq_parts).upper()
    # Remove any whitespace that might have been in sequence
    sequence_raw = re.sub(r'\s+', '', sequence_raw)

    # Validate DNA characters
    try:
        sequence = validate_dna_sequence(sequence_raw)
    except InvalidSequenceError as e:
        raise FileFormatError("<FASTA>", "FASTA", f"Invalid DNA sequence: {e}")

    # Parse header: first word is ID, rest is description
    header_parts = header.split(None, 1)
    seq_id = header_parts[0] if header_parts else "unknown"
    description = header_parts[1] if len(header_parts) > 1 else ""

    # Extract organism from pipe-delimited metadata if present
    organism = ""
    if "organism=" in header:
        match = re.search(r'organism=([^|\s]+)', header)
        if match:
            organism = match.group(1)

    return {
        "id": seq_id,
        "description": description,
        "sequence": sequence,
        "organism": organism,
        "gc_content": gc_content(sequence),
    }


# ─── GenBank Import ────────────────────────────────────────────────

def import_genbank(filepath_or_text: str) -> dict[str, Any]:
    """
    Parse a GenBank file or text string and return a structured dict.

    Returns a dict with keys:
        locus             — LOCUS name
        definition        — DEFINITION line text
        organism          — ORGANISM from SOURCE section
        sequence          — uppercase DNA sequence (from ORIGIN section)
        gene_name         — gene name from /gene qualifier (or "")
        exon_boundaries   — list of (start, end) tuples from CDS join()
                            converted from 1-based inclusive to 0-based exclusive
        protein           — protein translation from /translation qualifier (or "")
        features          — list of feature dicts with type, location, qualifiers
        gc_content        — GC fraction computed from sequence
        length            — sequence length

    Handles:
        - LOCUS, DEFINITION, ORGANISM, FEATURES, ORIGIN sections
        - CDS features with join() locations (multi-exon genes)
        - complement() locations for reverse-strand features
        - /gene, /translation, /note qualifiers
        - Line numbers and spacing in ORIGIN section
        - GenBank terminator //

    Args:
        filepath_or_text: Path to a GenBank file, or raw GenBank text content

    Returns:
        Dict with parsed GenBank record data

    Raises:
        FileFormatError: If the input is not valid GenBank or contains invalid DNA
    """
    text = _resolve_input(filepath_or_text, "GenBank")

    # Normalize line endings
    text = text.replace("\r\n", "\n").replace("\r", "\n")

    # Parse top-level fields
    locus = _extract_locus(text)
    definition = _extract_definition(text)
    organism = _extract_organism(text)
    sequence = _extract_origin_sequence(text)
    features = _extract_features(text)

    # Extract gene name from /gene qualifiers
    gene_name = ""
    for feat in features:
        if feat.get("type") in ("gene", "CDS"):
            gene_val = feat.get("qualifiers", {}).get("gene", "")
            if gene_val:
                gene_name = gene_val
                break

    # Extract exon boundaries from CDS join() locations
    exon_boundaries = []
    cds_feature = None
    for feat in features:
        if feat.get("type") == "CDS":
            cds_feature = feat
            break

    if cds_feature:
        exon_boundaries = _parse_exon_boundaries(cds_feature.get("location", ""))

    # Extract protein from /translation qualifier
    protein = ""
    if cds_feature:
        protein = cds_feature.get("qualifiers", {}).get("translation", "")

    # Validate DNA sequence
    if sequence:
        try:
            sequence = validate_dna_sequence(sequence)
        except InvalidSequenceError as e:
            raise FileFormatError("<GenBank>", "GenBank", f"Invalid DNA sequence in ORIGIN: {e}")

    gc = gc_content(sequence) if sequence else 0.0

    result = {
        "locus": locus,
        "definition": definition,
        "organism": organism,
        "sequence": sequence,
        "gene_name": gene_name,
        "exon_boundaries": exon_boundaries,
        "protein": protein,
        "features": features,
        "gc_content": gc,
        "length": len(sequence),
    }

    logger.debug("Parsed GenBank record: locus=%s, len=%d, organism=%s", locus, len(sequence), organism)
    return result


def _extract_locus(text: str) -> str:
    """Extract the LOCUS name from a GenBank record."""
    match = re.search(r'^LOCUS\s+(\S+)', text, re.MULTILINE)
    if match:
        return match.group(1)
    logger.warning("GenBank record has no LOCUS line")
    return ""


def _extract_definition(text: str) -> str:
    """Extract the DEFINITION text from a GenBank record."""
    match = re.search(r'^DEFINITION\s+(.+?)(?=\n[A-Z]{2,}|\n\n)', text, re.MULTILINE | re.DOTALL)
    if match:
        # Join continuation lines (indented with whitespace)
        definition = match.group(1)
        definition = re.sub(r'\s+', ' ', definition).strip()
        # Remove trailing period if present
        if definition.endswith('.'):
            definition = definition[:-1]
        return definition
    logger.debug("GenBank record has no DEFINITION line")
    return ""


def _extract_organism(text: str) -> str:
    """Extract the ORGANISM from the SOURCE section."""
    match = re.search(r'^\s+ORGANISM\s+(.+?)$', text, re.MULTILINE)
    if match:
        return match.group(1).strip()
    logger.debug("GenBank record has no ORGANISM line")
    return ""


def _extract_origin_sequence(text: str) -> str:
    """Extract and clean the DNA sequence from the ORIGIN section.

    The ORIGIN section has lines like:
        1 atggtgcatc tgactcctga ggagaagtct gcggtaccct cttctgcatc tttcatacgg
    We need to remove line numbers and spaces, keeping only ACGTN characters.
    """
    # Find ORIGIN section (between ORIGIN and //)
    origin_match = re.search(r'^ORIGIN\s*\n(.*?)^//', text, re.MULTILINE | re.DOTALL)
    if not origin_match:
        logger.warning("GenBank record has no ORIGIN section — no sequence data")
        return ""

    origin_text = origin_match.group(1)
    # Remove line numbers (digits at start of each line) and spaces
    # Keep only valid DNA characters
    seq_chars = re.sub(r'[^ACGTNacgtn]', '', origin_text)
    return seq_chars.upper()


def _extract_features(text: str) -> list[dict[str, Any]]:
    """Extract features from the FEATURES section of a GenBank record.

    Returns a list of dicts, each with:
        type       — feature type (gene, CDS, exon, misc_feature, etc.)
        location   — raw location string
        qualifiers — dict of qualifier key -> value
        strand     — '+' or '-' based on complement()
    """
    # Find FEATURES section (between FEATURES and ORIGIN)
    features_match = re.search(
        r'^FEATURES\s+Location/Qualifiers\s*\n(.*?)^ORIGIN',
        text, re.MULTILINE | re.DOTALL
    )
    if not features_match:
        logger.debug("GenBank record has no FEATURES section")
        return []

    features_text = features_match.group(1)

    # Split into individual features
    # A new feature starts with: FEATURE_INDENT spaces + feature_type + spaces + location
    feature_pattern = re.compile(
        rf'^ {{{_GENBANK_FEATURE_INDENT}}}(\S+)\s+(.+?)$', re.MULTILINE
    )
    feature_starts = list(feature_pattern.finditer(features_text))

    features = []
    for i, match in enumerate(feature_starts):
        feat_type = match.group(1)
        location_start = match.group(2).strip()

        # Get the full feature text (from this match to the next, or end)
        text_start = match.start()
        text_end = feature_starts[i + 1].start() if i + 1 < len(feature_starts) else len(features_text)
        feature_block = features_text[text_start:text_end]

        # Parse location (may span multiple lines for join/complement)
        location = _parse_feature_location(feature_block, feat_type)

        # Parse qualifiers
        qualifiers = _parse_feature_qualifiers(feature_block)

        # Determine strand
        strand = "-"
        if "complement" not in location:
            strand = "+"

        features.append({
            "type": feat_type,
            "location": location,
            "qualifiers": qualifiers,
            "strand": strand,
        })

    return features


def _parse_feature_location(feature_block: str, feat_type: str) -> str:
    """Parse the full location string from a feature block.

    The location may span multiple indented lines.
    """
    lines = feature_block.split("\n")
    location_parts = []

    # First line: "     TYPE   location_start"
    first_line = lines[0]
    # Remove the feature type prefix (FEATURE_INDENT spaces + type + spaces)
    loc_match = re.match(
        rf'^ {{{_GENBANK_FEATURE_INDENT}}}{re.escape(feat_type)}\s+(.+)$', first_line
    )
    if loc_match:
        location_parts.append(loc_match.group(1).strip())

    # Continuation lines for location: QUALIFIER_INDENT spaces of indentation
    for line in lines[1:]:
        stripped = line.rstrip()
        if not stripped:
            continue
        # Qualifier lines start with QUALIFIER_INDENT spaces then /
        if re.match(rf'^ {{{_GENBANK_QUALIFIER_INDENT}}}/', stripped):
            break
        # Location continuation: QUALIFIER_INDENT spaces, no /
        if re.match(rf'^ {{{_GENBANK_QUALIFIER_INDENT}}}\S', stripped):
            location_parts.append(stripped.strip())

    return " ".join(location_parts)


def _parse_feature_qualifiers(feature_block: str) -> dict[str, str]:
    """Parse qualifier key-value pairs from a feature block."""
    qualifiers = {}

    # Find all qualifier lines: /key="value" or /key=value
    # Qualifiers may span multiple lines with continuation quotes
    qual_pattern = re.compile(
        rf'^ {{{_GENBANK_QUALIFIER_INDENT}}}/(\S+?)=(.+)$', re.MULTILINE
    )
    matches = list(qual_pattern.finditer(feature_block))

    for i, match in enumerate(matches):
        key = match.group(1)
        value = match.group(2)

        # Get the full value text (from this match to next qualifier or end)
        value_start = match.start(2)
        value_end = matches[i + 1].start() if i + 1 < len(matches) else len(feature_block)
        full_value = feature_block[value_start:value_end]

        # Clean up the value
        full_value = _clean_qualifier_value(full_value)
        qualifiers[key] = full_value

    return qualifiers


def _clean_qualifier_value(raw: str) -> str:
    """Clean a qualifier value, joining multi-line quoted strings."""
    # Remove leading/trailing whitespace from each line, then join
    lines = raw.split("\n")
    cleaned = " ".join(line.strip() for line in lines)

    # Strip outer whitespace first so quote detection works reliably
    cleaned = cleaned.strip()

    # Remove surrounding quotes if present
    if cleaned.startswith('"') and cleaned.endswith('"'):
        cleaned = cleaned[1:-1]

    # Remove internal quote pairing for multi-line translations
    # GenBank format: "MGKL..." "VTVL..." -> "MGKL...VTVL..."
    cleaned = re.sub(r'"\s+"', '', cleaned)

    return cleaned.strip()


def _parse_exon_boundaries(location: str) -> list[tuple[int, int]]:
    """Parse exon boundaries from a CDS location string.

    Handles:
        join(1..92,273..495)            -> [(0, 92), (272, 495)]
        join(1..92,273..495,600..900)   -> [(0, 92), (272, 495), (599, 900)]
        1..495                          -> [(0, 495)]
        complement(join(1..92,273..495))-> [(0, 92), (272, 495)]

    GenBank convention: 1-based inclusive
    Internal convention: 0-based exclusive (Python-style)
    Conversion: 1-based inclusive (start, end) -> 0-based exclusive (start-1, end)
    """
    boundaries = []

    # Remove complement() wrapper if present
    loc = location.strip()
    loc = re.sub(r'complement\((.+)\)', r'\1', loc)

    # Check for join()
    join_match = re.match(r'join\((.+)\)', loc)
    if join_match:
        parts_str = join_match.group(1)
        parts = [p.strip() for p in parts_str.split(",")]
    else:
        # Single location
        parts = [loc]

    for part in parts:
        # Parse range: start..end
        range_match = re.match(r'(\d+)\.\.(\d+)', part)
        if range_match:
            gb_start = int(range_match.group(1))
            gb_end = int(range_match.group(2))
            # Convert 1-based inclusive to 0-based exclusive
            boundaries.append((gb_start - 1, gb_end))
        else:
            # Single position (e.g., in a point feature)
            pos_match = re.match(r'(\d+)', part)
            if pos_match:
                pos = int(pos_match.group(1))
                boundaries.append((pos - 1, pos))
            else:
                logger.debug("Skipping unparseable exon boundary part: %r", part)

    return boundaries


# ─── Auto-detect Format ────────────────────────────────────────────

def import_sequence(filepath_or_text: str) -> dict[str, Any]:
    """
    Auto-detect file format and parse accordingly.

    Detects between:
        - FASTA: starts with '>' character
        - GenBank: starts with 'LOCUS' keyword
        - Plain text: raw DNA sequence (no header)

    Returns a dict with keys:
        format           — detected format: "fasta", "genbank", or "plain"
        records          — list of parsed records (for FASTA, multiple possible)
        sequence         — primary DNA sequence (for convenience)
        gc_content       — GC fraction of primary sequence
        organism         — organism if detected (from FASTA or GenBank)
        gene_name        — gene name if detected (from GenBank)
        exon_boundaries  — exon boundaries if detected (from GenBank)
        protein          — protein translation if detected (from GenBank)
        locus            — GenBank locus name (or "")
        definition       — GenBank definition (or "")
        features         — list of features (from GenBank)

    Args:
        filepath_or_text: Path to a file, or raw text content

    Returns:
        Structured dict with format-specific parsed data

    Raises:
        FileFormatError: If the format cannot be determined or parsing fails
    """
    text = _resolve_input(filepath_or_text, "auto-detect")
    # Normalize line endings
    text = text.replace("\r\n", "\n").replace("\r", "\n")

    stripped = text.strip()

    # Detect format
    if stripped.startswith(">"):
        # FASTA format
        records = import_fasta(filepath_or_text)
        primary = records[0]
        return {
            "format": "fasta",
            "records": records,
            "sequence": primary["sequence"],
            "gc_content": primary["gc_content"],
            "organism": primary.get("organism", ""),
            "gene_name": "",
            "exon_boundaries": [],
            "protein": "",
            "locus": primary.get("id", ""),
            "definition": primary.get("description", ""),
            "features": [],
        }

    elif stripped.startswith("LOCUS"):
        # GenBank format
        gb = import_genbank(filepath_or_text)
        return {
            "format": "genbank",
            "records": [gb],
            "sequence": gb["sequence"],
            "gc_content": gb["gc_content"],
            "organism": gb.get("organism", ""),
            "gene_name": gb.get("gene_name", ""),
            "exon_boundaries": gb.get("exon_boundaries", []),
            "protein": gb.get("protein", ""),
            "locus": gb.get("locus", ""),
            "definition": gb.get("definition", ""),
            "features": gb.get("features", []),
        }

    else:
        # Plain text — treat as raw DNA sequence
        # Remove all whitespace
        seq_raw = re.sub(r'\s+', '', stripped).upper()

        try:
            sequence = validate_dna_sequence(seq_raw)
        except InvalidSequenceError as e:
            raise FileFormatError(
                filepath_or_text if _looks_like_path(filepath_or_text) else "<text>",
                "plain",
                f"Invalid DNA sequence: {e}",
            )

        if not sequence:
            raise FileFormatError(
                filepath_or_text if _looks_like_path(filepath_or_text) else "<text>",
                "auto-detect",
                "Empty sequence and unrecognized format",
            )

        return {
            "format": "plain",
            "records": [{"sequence": sequence, "gc_content": gc_content(sequence)}],
            "sequence": sequence,
            "gc_content": gc_content(sequence),
            "organism": "",
            "gene_name": "",
            "exon_boundaries": [],
            "protein": "",
            "locus": "",
            "definition": "",
            "features": [],
        }


# ─── Helper Functions ──────────────────────────────────────────────

def _resolve_input(filepath_or_text: str, format_name: str) -> str:
    """Resolve input: if it looks like a file path and the file exists, read it; otherwise treat as text."""
    if _looks_like_path(filepath_or_text):
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


def _looks_like_path(text: str) -> bool:
    """Heuristic: check if the text looks like a file path rather than sequence data."""
    # If it contains common path characters and is relatively short
    # Sequence data typically has lots of ACGT characters
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
    if re.search(r'\.(fasta|fa|fna|gb|gbk|genbank|txt)$', text, re.IGNORECASE):
        return True
    # If it's a single line with no whitespace and looks like a path
    if "\n" not in text and "." in text and not text.startswith(">") and not text.startswith("LOCUS"):
        # Could be a path like "gene.gb" or "data/gene.fasta"
        if re.match(r'^[\w./\\:_-]+$', text):
            return True
    return False
