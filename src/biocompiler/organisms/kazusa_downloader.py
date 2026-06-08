"""
Kazusa Codon Usage Database auto-downloader.

Fetches codon usage tables from https://www.kazusa.or.jp/codon/ by TaxID
and converts them into BioCompiler's CodonUsageTable format.

The Kazusa database provides codon usage for 35,000+ organisms.
This module enables BioCompiler to support virtually any organism
with a codon frequency table, matching DNAchisel's organism flexibility.

Usage::

    from biocompiler.organisms.kazusa_downloader import (
        fetch_codon_usage_from_kazusa,
        fetch_codon_usage_by_name,
        register_dynamic_organism,
        resolve_or_download_organism,
    )

    # Fetch by NCBI TaxID
    usage = fetch_codon_usage_from_kazusa(83333)  # E. coli K-12

    # Fetch by organism name
    usage = fetch_codon_usage_by_name("Thermus thermophilus")

    # Register a dynamically downloaded organism
    canonical = register_dynamic_organism("Thermus_thermophilus", usage)

    # Resolve or download in one step
    canonical = resolve_or_download_organism("Thermus thermophilus", taxid=274)
"""

from __future__ import annotations

import logging
import re
import urllib.request
import urllib.parse
import warnings
from typing import TYPE_CHECKING

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

if TYPE_CHECKING:
    pass

__all__ = [
    "fetch_codon_usage_from_kazusa",
    "fetch_codon_usage_by_name",
    "register_dynamic_organism",
    "resolve_or_download_organism",
    "_parse_kazusa_html",
    "clear_cache",
]

_logger = logging.getLogger(__name__)

# Kazusa Codon Usage Database URL template.
# species = NCBI TaxID; aa=1 shows all amino acids; style=N is plain format.
_KAZUSA_URL_TEMPLATE = (
    "https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi"
    "?species={taxid}&aa=1&style=N"
)

# Kazusa search URL for name-based lookup.
# Uses the findspecie.cgi endpoint with speciesname parameter.
_KAZUSA_SEARCH_URL = (
    "https://www.kazusa.or.jp/codon/cgi-bin/findspecie.cgi"
    "?speciesname={name}"
)

# In-memory cache: key is canonical string → CodonUsageTable
_CACHE: dict[str, CodonUsageTable] = {}

# Standard genetic code: codon → amino acid (one-letter)
_CODON_TO_AA: dict[str, str] = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# All 64 standard codons
_ALL_CODONS: list[str] = list(_CODON_TO_AA.keys())


def clear_cache() -> None:
    """Clear the in-memory Kazusa download cache."""
    _CACHE.clear()


def _fill_missing_codons(table: CodonUsageTable) -> CodonUsageTable:
    """Fill in any missing codons with zero frequency.

    Ensures the table has all 64 standard codons.  Missing codons
    receive (amino_acid, 0.0, 0.0, 0).

    Args:
        table: Partial codon usage table.

    Returns:
        Complete table with all 64 codons.
    """
    for codon in _ALL_CODONS:
        if codon not in table:
            aa = _CODON_TO_AA[codon]
            table[codon] = (aa, 0.0, 0.0, 0)
    return table


def _parse_kazusa_html(html: str) -> CodonUsageTable:
    """Parse Kazusa Codon Usage Database HTML format.

    The Kazusa HTML page contains a ``<pre>`` block with codon usage
    data in the following format (one line per codon)::

        Ami-Acid  Codon  Number  /1000  Fraction
        Gly       GGG    15201   16.63  0.25
        Gly       GGA    17333   18.96  0.28
        ...

    Each row has: amino acid name (3-letter or 1-letter), codon,
    raw count, per-thousand value, and fraction (0.0–1.0).

    Some Kazusa pages also use a compact format like::

        UUU 17.6(714298)  UUC 20.3(824692)

    This parser handles both formats.

    Args:
        html: Raw HTML string from the Kazusa database.

    Returns:
        CodonUsageTable mapping codon strings to
        (amino_acid, fraction, per_thousand, count) tuples.
        Missing codons are filled with zero frequency.

    Raises:
        ValueError: If no codon data could be parsed from the HTML.
    """
    # The Kazusa HTML format puts codon data inside <pre> tags.
    # We try to extract the pre block first, then fall back to
    # searching the entire HTML.
    pre_match = re.search(r"<pre>(.*?)</pre>", html, re.DOTALL | re.IGNORECASE)
    text = pre_match.group(1) if pre_match else html

    # Strip HTML tags from the text
    text = re.sub(r"<[^>]+>", "", text)

    # Known 3-letter amino acid codes to 1-letter mapping
    aa3_to_aa1: dict[str, str] = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "End": "*", "Ter": "*", "SEL": "U", "PYL": "O",
    }

    result: CodonUsageTable = {}

    # Parse lines matching the standard Kazusa codon data pattern.
    # Each line has: AA_name  CODON  count  per_thousand  fraction
    # We use a flexible regex to handle varying whitespace.
    pattern = re.compile(
        r"^\s*"
        r"(\w{1,4})\s+"        # amino acid name (1-4 chars)
        r"([ACGTU]{3})\s+"     # codon (3 nucleotides)
        r"(\d+)\s+"            # count (integer)
        r"([\d.]+)\s*"         # per-thousand (float)
        r"([\d.]+)\s*$",       # fraction (float)
        re.MULTILINE,
    )

    for match in pattern.finditer(text):
        aa_name = match.group(1)
        codon = match.group(2).upper().replace("U", "T")  # normalize RNA→DNA
        count = int(match.group(3))
        per_thousand = float(match.group(4))
        fraction = float(match.group(5))

        # Convert amino acid name to 1-letter code
        if len(aa_name) == 1 and aa_name in "ACDEFGHIKLMNPQRSTVWY*":
            aa = aa_name
        elif aa_name in aa3_to_aa1:
            aa = aa3_to_aa1[aa_name]
        else:
            # Fallback: use our built-in codon table
            aa = _CODON_TO_AA.get(codon, "X")

        result[codon] = (aa, fraction, per_thousand, count)

    # If the standard regex didn't match, try the compact format:
    # "UUU 17.6(714298)  UUC 20.3(824692)"
    # This format has codon, per_thousand in parentheses
    if not result:
        compact_pattern = re.compile(
            r"([ACGTU]{3})\s+"           # codon
            r"([\d.]+)\((\d+)\)",        # per_thousand(count)
        )
        for match in compact_pattern.finditer(text):
            codon = match.group(1).upper().replace("U", "T")
            per_thousand = float(match.group(2))
            count = int(match.group(3))
            aa = _CODON_TO_AA.get(codon, "X")

            # Compute fraction from per_thousand values per amino acid
            # We'll do a second pass to compute fractions
            result[codon] = (aa, 0.0, per_thousand, count)

        # Second pass: compute fractions per amino acid
        if result:
            aa_total_per_thousand: dict[str, float] = {}
            for codon, (aa, _frac, pt, _cnt) in result.items():
                aa_total_per_thousand[aa] = aa_total_per_thousand.get(aa, 0.0) + pt

            for codon, (aa, _frac, pt, cnt) in result.items():
                total = aa_total_per_thousand.get(aa, 0.0)
                fraction = pt / total if total > 0 else 0.0
                result[codon] = (aa, fraction, pt, cnt)

    # If still no match, try an alternate table format with columns
    # separated by spaces/tabs in a different layout.
    if not result:
        alt_pattern = re.compile(
            r"([ACGTU]{3})\s+(\d+)\s+([\d.]+)\s+([\d.]+)",
            re.IGNORECASE,
        )
        for match in alt_pattern.finditer(text):
            codon = match.group(1).upper().replace("U", "T")
            count = int(match.group(2))
            per_thousand = float(match.group(3))
            fraction = float(match.group(4))
            aa = _CODON_TO_AA.get(codon, "X")
            result[codon] = (aa, fraction, per_thousand, count)

    if not result:
        raise ValueError(
            "Could not parse any codon data from the Kazusa HTML response. "
            "The page format may have changed or the TaxID may be invalid."
        )

    # Fill in any missing codons with zero frequency
    result = _fill_missing_codons(result)

    return result


def fetch_codon_usage_from_kazusa(
    taxid: int,
    timeout: float = 30.0,
) -> CodonUsageTable:
    """Fetch and parse codon usage from the Kazusa database by NCBI TaxID.

    Downloads the codon usage table for the organism identified by
    *taxid* from the Kazusa Codon Usage Database and converts it
    into BioCompiler's ``CodonUsageTable`` format.

    Results are cached in memory so that repeated calls for the same
    TaxID do not hit the network.

    On network failure, returns an empty ``CodonUsageTable`` with a
    warning rather than raising an exception.

    Args:
        taxid: NCBI Taxonomy ID (e.g., 83333 for E. coli K-12,
            9606 for Homo sapiens, 10090 for Mus musculus).
        timeout: Network timeout in seconds (default 30).

    Returns:
        CodonUsageTable mapping codon strings to
        (amino_acid, fraction, per_thousand, count) tuples.
        Returns an empty dict on network failure (with a warning).
    """
    cache_key = f"taxid:{taxid}"
    if cache_key in _CACHE:
        _logger.debug("Cache hit for TaxID %d", taxid)
        return _CACHE[cache_key]

    url = _KAZUSA_URL_TEMPLATE.format(taxid=taxid)
    _logger.info("Fetching codon usage from Kazusa for TaxID %d: %s", taxid, url)

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "BioCompiler/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            html = resp.read().decode("utf-8", errors="replace")

        codon_usage = _parse_kazusa_html(html)

        if len(codon_usage) < 20:
            warnings.warn(
                f"Parsed only {len(codon_usage)} codons from Kazusa for TaxID {taxid}. "
                f"The TaxID may be invalid or the page format may have changed.",
                stacklevel=2,
            )
            # Return the partial data rather than raising — better to have
            # some data than none, and missing codons are zero-filled.
    except Exception as exc:
        warnings.warn(
            f"Failed to fetch codon usage from Kazusa for TaxID {taxid}: {exc}. "
            f"Returning empty table.",
            stacklevel=2,
        )
        codon_usage: CodonUsageTable = {}

    _CACHE[cache_key] = codon_usage
    if codon_usage:
        _logger.info(
            "Successfully fetched %d codons for TaxID %d from Kazusa",
            len(codon_usage), taxid,
        )
    return codon_usage


def fetch_codon_usage_by_name(
    name: str,
    timeout: float = 30.0,
) -> CodonUsageTable:
    """Fetch codon usage from Kazusa by organism name.

    Searches the Kazusa Codon Usage Database by organism name using
    the ``findspecie.cgi`` endpoint.  This is less reliable than
    searching by TaxID because the name must match exactly what
    Kazusa has on record.

    If the name contains a space (e.g., ``"Escherichia coli"``),
    it is converted to ``+`` format for the URL query.

    On network failure, returns an empty ``CodonUsageTable`` with a
    warning rather than raising an exception.

    Args:
        name: Organism name (e.g., ``"Thermus thermophilus"``).
        timeout: Network timeout in seconds (default 30).

    Returns:
        CodonUsageTable mapping codon strings to
        (amino_acid, fraction, per_thousand, count) tuples.
        Returns an empty dict on network failure (with a warning).
    """
    cache_key = f"name:{name.lower()}"
    if cache_key in _CACHE:
        _logger.debug("Cache hit for organism name '%s'", name)
        return _CACHE[cache_key]

    # Convert name for URL: spaces to + for search
    search_name = name.replace(" ", "+")
    search_url = _KAZUSA_SEARCH_URL.format(
        name=urllib.parse.quote(search_name),
    )

    _logger.info("Searching Kazusa for organism '%s': %s", name, search_url)

    try:
        req = urllib.request.Request(search_url, headers={"User-Agent": "BioCompiler/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            html = resp.read().decode("utf-8", errors="replace")

        # The search results page may contain links to specific organism pages.
        # Try to extract a TaxID from the results page.
        taxid_match = re.search(r"species=(\d+)", html)
        if taxid_match:
            taxid = int(taxid_match.group(1))
            _logger.info("Found TaxID %d for organism '%s'", taxid, name)
            result = fetch_codon_usage_from_kazusa(taxid, timeout=timeout)
        else:
            # Try parsing the page directly — sometimes the search
            # returns a codon table directly for exact matches.
            result = _parse_kazusa_html(html)
    except Exception as exc:
        warnings.warn(
            f"Failed to fetch codon usage from Kazusa for organism '{name}': {exc}. "
            f"Returning empty table.",
            stacklevel=2,
        )
        result: CodonUsageTable = {}

    _CACHE[cache_key] = result
    return result


def register_dynamic_organism(
    name: str,
    codon_usage: CodonUsageTable,
    gc_target: tuple[float, float] | None = None,
    is_eukaryote: bool | None = None,
) -> str:
    """Register a dynamically downloaded organism into BioCompiler's runtime.

    Adds the organism to all relevant registries so that it can be
    used with ``resolve_organism()``, ``CODON_ADAPTIVENESS_TABLES``,
    ``CODON_USAGE_TABLES``, ``PREFERRED_CODON_TABLES``,
    ``ORGANISM_ALIASES``, and ``ORGANISM_GC_TARGETS``.

    Args:
        name: Canonical organism name (e.g., ``"Thermus_thermophilus"``).
            Spaces are automatically converted to underscores.
        codon_usage: Codon usage table from Kazusa or other source.
        gc_target: Optional (gc_lo, gc_hi) target range.  If ``None``,
            a default of (0.30, 0.70) is used.
        is_eukaryote: Whether the organism is eukaryotic.  If ``None``,
            the domain is not set (prokaryote assumed for GC defaults).

    Returns:
        The canonical name used for registration (with underscores
        replacing spaces).
    """
    # Import here to avoid circular imports at module level
    from . import (
        CODON_ADAPTIVENESS_TABLES,
        CODON_USAGE_TABLES,
        PREFERRED_CODON_TABLES,
        ORGANISM_ALIASES,
        ORGANISM_GC_TARGETS,
        SUPPORTED_ORGANISMS,
    )

    # Normalize name
    canonical = name.replace(" ", "_")

    # Compute derived tables
    adaptiveness = compute_codon_adaptiveness(codon_usage)
    preferred = compute_preferred_codons(codon_usage)

    # Register in all tables
    CODON_USAGE_TABLES[canonical] = codon_usage
    CODON_ADAPTIVENESS_TABLES[canonical] = adaptiveness
    PREFERRED_CODON_TABLES[canonical] = preferred

    # Add alias for the lowercase short name
    short_name = canonical.split("_")[0].lower() if "_" in canonical else canonical.lower()
    ORGANISM_ALIASES[canonical] = canonical
    ORGANISM_ALIASES[short_name] = canonical
    # Also register with spaces replaced by spaces for convenience
    name_with_spaces = canonical.replace("_", " ")
    ORGANISM_ALIASES[name_with_spaces] = canonical

    # Add alias entries in the main registries
    CODON_USAGE_TABLES[short_name] = codon_usage
    CODON_ADAPTIVENESS_TABLES[short_name] = adaptiveness
    PREFERRED_CODON_TABLES[short_name] = preferred

    # Set GC target
    if gc_target is None:
        gc_target = (0.30, 0.70)
    ORGANISM_GC_TARGETS[canonical] = gc_target

    # Update the SUPPORTED_ORGANISMS list
    if canonical not in SUPPORTED_ORGANISMS:
        SUPPORTED_ORGANISMS.append(canonical)
    if short_name not in SUPPORTED_ORGANISMS:
        SUPPORTED_ORGANISMS.append(short_name)

    _logger.info(
        "Registered dynamic organism '%s' (short: '%s') with %d codons",
        canonical, short_name, len(codon_usage),
    )

    return canonical


def resolve_or_download_organism(
    name: str,
    taxid: int | None = None,
    timeout: float = 30.0,
) -> str:
    """Resolve an organism name, downloading from Kazusa if not found.

    First tries ``resolve_organism()`` to find the organism in
    BioCompiler's built-in registries.  If the organism is not found,
    it attempts to download codon usage data from the Kazusa Codon
    Usage Database and register it dynamically.

    Args:
        name: Organism name or alias (e.g., ``"Thermus thermophilus"``).
        taxid: Optional NCBI Taxonomy ID for direct Kazusa lookup.
            If provided, this is used instead of a name-based search.
        timeout: Network timeout in seconds (default 30).

    Returns:
        The canonical organism name used for registration.

    Raises:
        ValueError: If the organism cannot be resolved or downloaded.
    """
    from . import resolve_organism

    # Try resolving from built-in registries first
    try:
        canonical = resolve_organism(name, strict=True)
        _logger.debug("Resolved '%s' to '%s' from built-in registries", name, canonical)
        return canonical
    except Exception:
        pass  # Not found — proceed to download

    _logger.info("Organism '%s' not found in built-in registries, attempting Kazusa download", name)

    # Download from Kazusa
    if taxid is not None:
        codon_usage = fetch_codon_usage_from_kazusa(taxid, timeout=timeout)
    else:
        codon_usage = fetch_codon_usage_by_name(name, timeout=timeout)

    if not codon_usage:
        raise ValueError(
            f"Could not resolve or download organism '{name}'. "
            f"Kazusa download returned no data. "
            f"Provide a valid TaxID or check the organism name."
        )

    # Register the dynamically downloaded organism
    canonical_name = name.replace(" ", "_")
    return register_dynamic_organism(canonical_name, codon_usage)
