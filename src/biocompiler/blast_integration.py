"""
BLAST+ Integration for BioCompiler Biosecurity Screening

Extends the biosecurity screening with homology-based search against
hazardous sequence databases using NCBI BLAST+ (local preferred,
remote fallback). This enables detection of sequences that share
significant similarity with known hazardous sequences beyond what
exact/fuzzy substring matching can find.

Components:
  - BlastScanner: Wrapper class for BLAST+ operations
  - build_hazard_db: Build a BLAST database from hazardous sequences
  - screen_blast: Run BLAST screening and return structured results
  - check_biosecurity_blast: Integration hook combining exact/fuzzy
    matching with BLAST homology search

Environment Variables:
  - BIOCOMPILER_BLAST_PATH: Custom path to BLAST+ binaries directory.
    If set, blastn/blastp/makeblastdb are looked up here first.
  - BIOCOMPILER_BLAST_DB_PATH: Path to a pre-built hazardous sequence
    BLAST database (without extension). If set, skips database building.

Fallback:
  If BLAST+ is not installed, all functions return gracefully with
  ``blast_available=False`` and a warning recommending installation.
  Optimization is NOT blocked by the absence of BLAST+.

References:
  - Altschul et al., Nucleic Acids Res 1990; 25:3389-402 (BLAST)
  - Camacho et al., BMC Bioinformatics 2009; 10:421 (BLAST+)
  - CDC Select Agent Program (42 CFR Part 73)
  - Australia Group Common Control List
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional

from .biosecurity import (
    BiosecurityReport,
    HazardMatch,
    RiskLevel,
    _HAZARD_SIGNATURES,
    _max_risk,
    _RISK_ORDER,
    screen_hazardous_sequence,
)
from .exceptions import BiosecurityError

logger = logging.getLogger(__name__)

__all__ = [
    "BlastScanner",
    "BlastHit",
    "BlastScreeningResult",
    "build_hazard_db",
    "screen_blast",
    "check_biosecurity_blast",
    "is_blast_available",
    "find_blast_bin",
]

# ─────────────────────────────────────────────────────────────────────────────
# BLAST+ binary detection
# ─────────────────────────────────────────────────────────────────────────────

_BLAST_BINARIES = ["blastn", "blastp", "tblastx", "makeblastdb"]


def find_blast_bin(name: str) -> Optional[str]:
    """Find a BLAST+ binary on the system PATH or BIOCOMPILER_BLAST_PATH.

    Parameters
    ----------
    name : str
        Binary name (e.g. ``"blastn"``, ``"makeblastdb"``).

    Returns
    -------
    str or None
        Full path to the binary, or None if not found.
    """
    # Check custom path first
    custom_path = os.environ.get("BIOCOMPILER_BLAST_PATH", "")
    if custom_path:
        candidate = Path(custom_path) / name
        if candidate.is_file():
            return str(candidate)
        # Also check with common suffixes
        for suffix in ("", ".exe"):
            candidate = Path(custom_path) / (name + suffix)
            if candidate.is_file():
                return str(candidate)

    # Fall back to system PATH
    result = shutil.which(name)
    return result


def is_blast_available() -> bool:
    """Check whether BLAST+ tools are available on this system.

    Returns True if at least ``blastn`` and ``makeblastdb`` are found.
    """
    return find_blast_bin("blastn") is not None and find_blast_bin("makeblastdb") is not None


# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class BlastHit:
    """A single BLAST alignment hit."""

    subject_id: str
    percent_identity: float
    alignment_length: int
    evalue: float
    bit_score: float
    query_coverage: float = 0.0  # fraction of query aligned (0.0-1.0)
    subject_length: int = 0
    query_start: int = 0
    query_end: int = 0
    subject_start: int = 0
    subject_end: int = 0


@dataclass
class BlastScreeningResult:
    """Result of BLAST-based biosecurity screening."""

    blast_available: bool
    hits: list[BlastHit] = field(default_factory=list)
    risk_level: RiskLevel = "none"
    coverage: float = 0.0  # max query coverage across all hits
    best_identity: float = 0.0  # highest percent identity
    database_path: str = ""
    blast_program: str = ""  # "blastn" or "blastp"
    error: Optional[str] = None  # error message if BLAST failed
    warnings: list[str] = field(default_factory=list)


# ─────────────────────────────────────────────────────────────────────────────
# BLAST+ wrapper class
# ─────────────────────────────────────────────────────────────────────────────


class BlastScanner:
    """Wrapper for NCBI BLAST+ tools for biosecurity screening.

    Parameters
    ----------
    db_path : str
        Path to the BLAST database (without file extension).
        E.g. ``"/data/hazard_db/nt_hazard"``.
    evalue : float
        E-value threshold for reporting hits (default 1e-5).
    max_hits : int
        Maximum number of hits to report (default 50).
    """

    def __init__(
        self,
        db_path: str,
        evalue: float = 1e-5,
        max_hits: int = 50,
    ) -> None:
        self.db_path = db_path
        self.evalue = evalue
        self.max_hits = max_hits
        self._blastn_path = find_blast_bin("blastn")
        self._blastp_path = find_blast_bin("blastp")
        self._tblastx_path = find_blast_bin("tblastx")
        self._makeblastdb_path = find_blast_bin("makeblastdb")

    @property
    def available(self) -> bool:
        """Whether the required BLAST+ tools are available."""
        return (
            self._blastn_path is not None
            and self._makeblastdb_path is not None
        )

    @property
    def protein_available(self) -> bool:
        """Whether protein BLAST (blastp) is available."""
        return self._blastp_path is not None

    def screen(
        self,
        sequence: str,
        seq_type: Literal["dna", "protein"] = "dna",
    ) -> BlastScreeningResult:
        """Run a BLAST search against the configured database.

        Parameters
        ----------
        sequence : str
            Query sequence (nucleotide or amino acid).
        seq_type : str
            ``"dna"`` for blastn, ``"protein"`` for blastp.

        Returns
        -------
        BlastScreeningResult
            Structured result with hits, risk level, and coverage.
        """
        if not self.available:
            return BlastScreeningResult(
                blast_available=False,
                warnings=[
                    "BLAST+ is not installed. Install NCBI BLAST+ to enable "
                    "homology-based biosecurity screening. "
                    "See: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/"
                ],
            )

        # Choose the appropriate BLAST program
        if seq_type == "protein":
            if not self.protein_available:
                return BlastScreeningResult(
                    blast_available=True,
                    error="blastp is not available for protein screening",
                    blast_program="blastp",
                    database_path=self.db_path,
                )
            program = "blastp"
            binary = self._blastp_path
        else:
            program = "blastn"
            binary = self._blastn_path

        # Run BLAST
        try:
            raw_output = self._run_blast(binary, sequence)
        except subprocess.TimeoutExpired:
            error_msg = f"BLAST {program} timed out after 120 seconds"
            logger.error(error_msg)
            return BlastScreeningResult(
                blast_available=True,
                error=error_msg,
                blast_program=program,
                database_path=self.db_path,
            )
        except subprocess.CalledProcessError as e:
            error_msg = f"BLAST {program} failed with exit code {e.returncode}: {e.stderr[:500]}"
            logger.error(error_msg)
            return BlastScreeningResult(
                blast_available=True,
                error=error_msg,
                blast_program=program,
                database_path=self.db_path,
            )
        except Exception as e:
            error_msg = f"BLAST {program} error: {e}"
            logger.error(error_msg)
            return BlastScreeningResult(
                blast_available=True,
                error=error_msg,
                blast_program=program,
                database_path=self.db_path,
            )

        # Parse output
        hits = self._parse_xml_output(raw_output)
        query_length = len(sequence)

        # Compute coverage and risk level
        best_identity = 0.0
        max_coverage = 0.0
        for hit in hits:
            if hit.percent_identity > best_identity:
                best_identity = hit.percent_identity
            if query_length > 0:
                cov = (hit.query_end - hit.query_start + 1) / query_length
                hit.query_coverage = round(cov, 4)
                if cov > max_coverage:
                    max_coverage = cov

        risk_level = _compute_blast_risk_level(best_identity, max_coverage)

        return BlastScreeningResult(
            blast_available=True,
            hits=hits,
            risk_level=risk_level,
            coverage=round(max_coverage, 4),
            best_identity=round(best_identity, 2),
            database_path=self.db_path,
            blast_program=program,
        )

    def _run_blast(self, binary: str, sequence: str) -> str:
        """Execute BLAST and return raw XML output.

        Parameters
        ----------
        binary : str
            Path to the BLAST binary (blastn or blastp).
        sequence : str
            Query sequence.

        Returns
        -------
        str
            Raw XML output from BLAST.

        Raises
        ------
        subprocess.CalledProcessError
            If BLAST exits with non-zero status.
        subprocess.TimeoutExpired
            If BLAST runs longer than 120 seconds.
        """
        cmd = [
            binary,
            "-query", "-",  # read from stdin
            "-db", self.db_path,
            "-evalue", str(self.evalue),
            "-max_target_seqs", str(self.max_hits),
            "-outfmt", "5",  # XML output
            "-dust", "no",  # disable dust filtering for short queries
        ]

        logger.info("Running BLAST: %s -db %s -evalue %s", binary, self.db_path, self.evalue)

        result = subprocess.run(
            cmd,
            input=sequence,
            capture_output=True,
            text=True,
            timeout=120,
            check=True,
        )

        return result.stdout

    def _parse_xml_output(self, xml_output: str) -> list[BlastHit]:
        """Parse BLAST XML output into structured hits.

        Parameters
        ----------
        xml_output : str
            Raw XML output from BLAST -outfmt 5.

        Returns
        -------
        list[BlastHit]
            Parsed BLAST hits.
        """
        hits: list[BlastHit] = []

        if not xml_output or not xml_output.strip():
            return hits

        try:
            root = ET.fromstring(xml_output)
        except ET.ParseError as e:
            logger.warning("Failed to parse BLAST XML output: %s", e)
            return hits

        # Navigate the XML: BlastOutput -> BlastOutput_iterations -> Iteration -> Iteration_hits -> Hit
        for iteration in root.iter("Iteration"):
            hit_list = iteration.find("Iteration_hits")
            if hit_list is None:
                continue
            for hit_elem in hit_list.findall("Hit"):
                hit_id_elem = hit_elem.find("Hit_id")
                hit_def_elem = hit_elem.find("Hit_def")
                subject_id = ""
                if hit_id_elem is not None and hit_id_elem.text:
                    subject_id = hit_id_elem.text
                if hit_def_elem is not None and hit_def_elem.text:
                    if subject_id:
                        subject_id += " " + hit_def_elem.text
                    else:
                        subject_id = hit_def_elem.text

                hit_hsps = hit_elem.find("Hit_hsps")
                if hit_hsps is None:
                    continue

                for hsp in hit_hsps.findall("Hsp"):
                    try:
                        pct_identity = float(_hsp_text(hsp, "Hsp_identity", "0"))
                        align_len = int(_hsp_text(hsp, "Hsp_align-len", "0"))
                        evalue = float(_hsp_text(hsp, "Hsp_evalue", "1"))
                        bit_score = float(_hsp_text(hsp, "Hsp_bit-score", "0"))
                        query_start = int(_hsp_text(hsp, "Hsp_query-from", "0"))
                        query_end = int(_hsp_text(hsp, "Hsp_query-to", "0"))
                        subject_start = int(_hsp_text(hsp, "Hsp_hit-from", "0"))
                        subject_end = int(_hsp_text(hsp, "Hsp_hit-to", "0"))
                        subject_len = int(_hsp_text(hsp, "Hsp_hit-len", "0"))

                        # Compute percent identity as percentage
                        if align_len > 0:
                            pct_identity = (pct_identity / align_len) * 100.0
                        else:
                            pct_identity = 0.0

                        hits.append(BlastHit(
                            subject_id=subject_id,
                            percent_identity=round(pct_identity, 2),
                            alignment_length=align_len,
                            evalue=evalue,
                            bit_score=bit_score,
                            subject_length=subject_len,
                            query_start=query_start,
                            query_end=query_end,
                            subject_start=subject_start,
                            subject_end=subject_end,
                        ))
                    except (ValueError, TypeError) as e:
                        logger.warning("Failed to parse HSP: %s", e)
                        continue

        # Sort hits by bit_score descending
        hits.sort(key=lambda h: h.bit_score, reverse=True)
        return hits[:self.max_hits]


def _hsp_text(hsp_elem: ET.Element, tag: str, default: str = "") -> str:
    """Safely extract text from an HSP element."""
    child = hsp_elem.find(tag)
    if child is not None and child.text:
        return child.text
    return default


# ─────────────────────────────────────────────────────────────────────────────
# Risk level computation from BLAST results
# ─────────────────────────────────────────────────────────────────────────────


def _compute_blast_risk_level(best_identity: float, coverage: float) -> RiskLevel:
    """Compute risk level from the best BLAST hit.

    Thresholds:
      - identity > 90% and coverage > 50% → critical
      - identity > 90% (any coverage) → high
      - identity > 70% and coverage > 50% → high
      - identity > 70% (any coverage) → medium
      - identity > 50% and coverage > 50% → medium
      - identity > 50% (any coverage) → low
      - Otherwise → none

    Parameters
    ----------
    best_identity : float
        Highest percent identity across all hits (0-100).
    coverage : float
        Maximum query coverage across all hits (0.0-1.0).

    Returns
    -------
    RiskLevel
        Computed risk level string.
    """
    if best_identity > 90.0:
        if coverage > 0.5:
            return "critical"
        return "high"
    elif best_identity > 70.0:
        if coverage > 0.5:
            return "high"
        return "medium"
    elif best_identity > 50.0:
        if coverage > 0.5:
            return "medium"
        return "low"
    return "none"


# ─────────────────────────────────────────────────────────────────────────────
# Database builder
# ─────────────────────────────────────────────────────────────────────────────


def build_hazard_db(
    sequences: dict[str, str],
    db_name: str,
    db_path: str,
    db_type: Literal["nucl", "prot"] = "nucl",
) -> str:
    """Build a BLAST database from a dictionary of sequences.

    Parameters
    ----------
    sequences : dict[str, str]
        Mapping of {sequence_id: sequence_string}.
    db_name : str
        Name for the BLAST database (used in makeblastdb -title).
    db_path : str
        Output path for the BLAST database (without extension).
        E.g. ``"/data/hazard_db/nt_hazard"``.
    db_type : str
        ``"nucl"`` for nucleotide database (blastn), ``"prot"``
        for protein database (blastp).

    Returns
    -------
    str
        The database path (same as *db_path*).

    Raises
    ------
    RuntimeError
        If makeblastdb is not available or fails.
    """
    makeblastdb = find_blast_bin("makeblastdb")
    if makeblastdb is None:
        raise RuntimeError(
            "makeblastdb not found. Install NCBI BLAST+ to build databases. "
            "See: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/"
        )

    if not sequences:
        raise ValueError("No sequences provided for database building")

    # Create the output directory if needed
    db_dir = Path(db_path).parent
    db_dir.mkdir(parents=True, exist_ok=True)

    # Write sequences to a temporary FASTA file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, dir=str(db_dir)
    ) as fasta_file:
        for seq_id, seq in sequences.items():
            # Sanitize seq_id (no spaces, no pipes in FASTA headers)
            safe_id = seq_id.replace(" ", "_").replace("|", "_")
            fasta_file.write(f">{safe_id}\n")
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                fasta_file.write(seq[i:i + 80] + "\n")
        fasta_path = fasta_file.name

    try:
        cmd = [
            makeblastdb,
            "-in", fasta_path,
            "-dbtype", db_type,
            "-out", db_path,
            "-title", db_name,
            "-parse_seqids",
        ]

        logger.info("Building BLAST database: %s (%d sequences, type=%s)", db_name, len(sequences), db_type)

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            check=True,
        )

        logger.info("BLAST database built successfully: %s", db_path)
        return db_path

    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"makeblastdb failed (exit code {e.returncode}): {e.stderr[:500]}"
        ) from e
    finally:
        # Clean up the temporary FASTA file
        try:
            os.unlink(fasta_path)
        except OSError:
            pass


def _build_hazard_sequences_from_signatures(
    seq_type: Literal["protein", "dna"] = "protein",
) -> dict[str, str]:
    """Extract sequences from the biosecurity hazard signature database.

    Note: The existing signatures are short motifs (8-12 aa / 15-21 nt),
    not full-length sequences. BLAST works best with longer sequences,
    but these motifs can still be used for exact/substring homology.

    Parameters
    ----------
    seq_type : str
        ``"protein"`` for protein signatures, ``"dna"`` for nucleotide.

    Returns
    -------
    dict[str, str]
        Mapping of {signature_name: motif_sequence}.
    """
    sequences: dict[str, str] = {}
    for sig in _HAZARD_SIGNATURES:
        if sig["type"] == seq_type:
            sequences[sig["name"]] = sig["motif"]
    return sequences


# ─────────────────────────────────────────────────────────────────────────────
# Main screening function
# ─────────────────────────────────────────────────────────────────────────────


def screen_blast(
    sequence: str,
    db_path: str,
    seq_type: Literal["dna", "protein"] = "dna",
    evalue: float = 1e-5,
    max_hits: int = 50,
) -> BlastScreeningResult:
    """Run BLAST screening against a hazardous sequence database.

    This is a convenience function that creates a :class:`BlastScanner`
    and runs the search in one call.

    Parameters
    ----------
    sequence : str
        Query sequence (nucleotide or amino acid).
    db_path : str
        Path to the BLAST database (without extension).
    seq_type : str
        ``"dna"`` for blastn, ``"protein"`` for blastp.
    evalue : float
        E-value threshold (default 1e-5).
    max_hits : int
        Maximum hits to report (default 50).

    Returns
    -------
    BlastScreeningResult
        Structured result with hits, risk level, and coverage.
    """
    scanner = BlastScanner(db_path=db_path, evalue=evalue, max_hits=max_hits)
    return scanner.screen(sequence, seq_type=seq_type)


# ─────────────────────────────────────────────────────────────────────────────
# Integration with biosecurity module
# ─────────────────────────────────────────────────────────────────────────────


def check_biosecurity_blast(
    protein: str,
    dna: str = "",
    db_path: Optional[str] = None,
    evalue: float = 1e-5,
    max_hits: int = 50,
    biosecurity_mode: Optional[str] = None,
) -> BiosecurityReport:
    """Combined biosecurity screening: exact/fuzzy matching + BLAST homology.

    Runs both the standard :func:`screen_hazardous_sequence` (substring
    matching) and, if BLAST+ is available, a BLAST homology search against
    a hazardous sequence database. Results from both methods are merged
    into a single :class:`BiosecurityReport`.

    If BLAST+ is not installed, the function falls back to exact/fuzzy
    matching only and adds a warning recommending BLAST+ installation.
    The optimization is NOT blocked by the absence of BLAST+.

    Parameters
    ----------
    protein : str
        Protein sequence in single-letter amino acid codes.
    dna : str, optional
        DNA sequence for resistance marker screening.
    db_path : str, optional
        Path to the BLAST database. If None, checks the
        ``BIOCOMPILER_BLAST_DB_PATH`` environment variable. If neither
        is set, BLAST screening is skipped.
    evalue : float
        E-value threshold for BLAST (default 1e-5).
    max_hits : int
        Maximum BLAST hits to report (default 50).
    biosecurity_mode : str, optional
        Override ``BIOCOMPILER_BIOSECURITY_MODE`` env var.

    Returns
    -------
    BiosecurityReport
        Merged biosecurity report from both screening methods.
    """
    from .biosecurity import get_biosecurity_mode

    # Step 1: Run exact/fuzzy matching (always)
    exact_report = screen_hazardous_sequence(protein, dna)

    # Step 2: Run BLAST screening (if available and db_path provided)
    blast_result: Optional[BlastScreeningResult] = None

    # Resolve db_path
    effective_db_path = db_path or os.environ.get("BIOCOMPILER_BLAST_DB_PATH", "")

    if effective_db_path and is_blast_available():
        # Protein BLAST
        protein_blast = screen_blast(
            protein, effective_db_path, seq_type="protein",
            evalue=evalue, max_hits=max_hits,
        )
        blast_result = protein_blast

        # DNA BLAST (if DNA is provided)
        if dna:
            dna_blast = screen_blast(
                dna, effective_db_path, seq_type="dna",
                evalue=evalue, max_hits=max_hits,
            )
            # Merge DNA BLAST results with protein BLAST results
            if blast_result.blast_available:
                blast_result.hits.extend(dna_blast.hits)
                # Recompute risk level from all hits combined
                all_identities = [h.percent_identity for h in blast_result.hits]
                best_id = max(all_identities) if all_identities else 0.0
                all_coverages = [h.query_coverage for h in blast_result.hits]
                best_cov = max(all_coverages) if all_coverages else 0.0
                blast_result.risk_level = _compute_blast_risk_level(best_id, best_cov)
                blast_result.best_identity = round(best_id, 2)
                blast_result.coverage = round(best_cov, 4)
            elif dna_blast.blast_available:
                blast_result = dna_blast
    elif not effective_db_path:
        logger.info(
            "BLAST screening skipped: no database path provided. "
            "Set db_path parameter or BIOCOMPILER_BLAST_DB_PATH env var."
        )

    # Step 3: Merge results
    return _merge_reports(exact_report, blast_result)


def _merge_reports(
    exact_report: BiosecurityReport,
    blast_result: Optional[BlastScreeningResult],
) -> BiosecurityReport:
    """Merge exact/fuzzy matching results with BLAST results.

    Parameters
    ----------
    exact_report : BiosecurityReport
        Results from exact/fuzzy substring matching.
    blast_result : BlastScreeningResult or None
        Results from BLAST homology search.

    Returns
    -------
    BiosecurityReport
        Merged report with combined risk assessment.
    """
    if blast_result is None:
        # No BLAST result — return the exact report with a note
        recommendations = list(exact_report.recommendations)
        recommendations.append(
            "BLAST+ homology screening was not performed. "
            "Install NCBI BLAST+ and configure a hazardous sequence database "
            "for comprehensive biosecurity screening."
        )
        return BiosecurityReport(
            is_hazardous=exact_report.is_hazardous,
            risk_level=exact_report.risk_level,
            flagged_categories=exact_report.flagged_categories,
            matches=exact_report.matches,
            recommendations=recommendations,
        )

    # Merge risk levels — take the highest
    combined_risk = _max_risk(exact_report.risk_level, blast_result.risk_level)

    # Collect all matches from exact screening
    all_matches = list(exact_report.matches)

    # Add BLAST hits as HazardMatch entries
    flagged_categories = set(exact_report.flagged_categories)
    for hit in blast_result.hits:
        category = _classify_blast_hit(hit)
        flagged_categories.add(category)

        confidence = _blast_identity_to_confidence(hit.percent_identity)
        all_matches.append(HazardMatch(
            category=category,
            name=f"BLAST:{hit.subject_id}",
            position=hit.query_start - 1 if hit.query_start > 0 else 0,  # 0-based
            matched_sequence=f"BLAST hit: {hit.percent_identity:.1f}% identity over {hit.alignment_length} residues",
            confidence=confidence,
            source=f"BLAST+ {blast_result.blast_program} (e={hit.evalue:.1e}, bits={hit.bit_score:.1f})",
        ))

    # Determine if hazardous
    is_hazardous = combined_risk in ("medium", "high", "critical")

    # Build recommendations
    recommendations = list(exact_report.recommendations)

    if blast_result.hits:
        best_hit = blast_result.hits[0]  # already sorted by bit_score
        recommendations.append(
            f"BLAST+ screening found {len(blast_result.hits)} homologous hit(s). "
            f"Best match: {best_hit.subject_id} "
            f"({best_hit.percent_identity:.1f}% identity, "
            f"e-value={best_hit.evalue:.1e}, "
            f"coverage={blast_result.coverage:.1%}). "
            f"BLAST risk level: {blast_result.risk_level}."
        )

    if not blast_result.blast_available:
        recommendations.append(
            "BLAST+ homology screening was not available. "
            "Install NCBI BLAST+ (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) "
            "and configure BIOCOMPILER_BLAST_DB_PATH for comprehensive screening. "
            "Only exact/fuzzy substring matching was performed."
        )
    elif blast_result.error:
        recommendations.append(
            f"BLAST+ screening encountered an error: {blast_result.error}. "
            f"Only exact/fuzzy substring matching results are available."
        )

    if combined_risk in ("high", "critical"):
        recommendations.append(
            "Optimization BLOCKED due to high/critical biosecurity risk "
            "(combined exact matching + BLAST homology screening). "
            "Resolve all flagged issues or obtain explicit institutional approval."
        )

    return BiosecurityReport(
        is_hazardous=is_hazardous,
        risk_level=combined_risk,
        flagged_categories=sorted(flagged_categories),
        matches=all_matches,
        recommendations=recommendations,
    )


def _classify_blast_hit(hit: BlastHit) -> str:
    """Classify a BLAST hit into a hazard category.

    Uses heuristics based on the subject_id to classify hits.
    """
    subject_lower = hit.subject_id.lower()

    # Select agent toxins
    toxin_keywords = [
        "ricin", "abrin", "botulinum", "shiga", "diphtheria",
        "tetanus", "cholera", "anthrax", "SEB", "enterotoxin",
        "mycotoxin",
    ]
    if any(kw.lower() in subject_lower for kw in toxin_keywords):
        return "select_agent"

    # Viral surface proteins
    viral_keywords = [
        "spike", "hemagglutinin", "neuraminidase", "gp41", "gp120",
        "gp160", "envelope", "ebola_gp", "variola", "influenza_ha",
        "influenza_na", "sars", "hiv_env",
    ]
    if any(kw.lower() in subject_lower for kw in viral_keywords):
        return "viral_surface"

    # Antibiotic resistance
    ar_keywords = [
        "bla", "tem", "ctx-m", "ndm", "kpc", "oxa", "vim", "imp",
        "nptii", "kanamycin", "chloramphenicol", "tetracycline",
        "vana", "meca", "aminoglycoside", "beta-lactamase",
        "aac(6", "cat ", "tet ", "tetm", "teto",
    ]
    if any(kw.lower() in subject_lower for kw in ar_keywords):
        return "antibiotic_resistance"

    # Oncogenes
    oncogene_keywords = [
        "myc", "ras", "egfr", "vegfa", "braf", "erbb2", "her2",
        "pdgf", "tgfb", "akt", "fos", "jun", "src",
    ]
    if any(kw.lower() in subject_lower for kw in oncogene_keywords):
        return "oncogene"

    # Default: unknown hazard
    return "blast_homology"


def _blast_identity_to_confidence(identity: float) -> float:
    """Convert BLAST percent identity to a confidence score (0.0-1.0).

    Thresholds:
      - >95% → 0.95
      - >90% → 0.90
      - >70% → 0.75
      - >50% → 0.55
      - Otherwise → 0.30
    """
    if identity > 95.0:
        return 0.95
    elif identity > 90.0:
        return 0.90
    elif identity > 70.0:
        return 0.75
    elif identity > 50.0:
        return 0.55
    else:
        return 0.30
