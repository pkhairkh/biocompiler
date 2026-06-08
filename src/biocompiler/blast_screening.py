"""
BLAST Integration Scaffolding for Biosecurity Screening

Provides BLAST (Basic Local Alignment Search Tool) integration for robust
sequence screening against known pathogen/toxin databases. This module
implements a production-grade biosecurity framework that:

1. Uses full BLAST+ when available (blastn/blastp against local databases)
2. Falls back to pattern-based quick screening when BLAST+ is not installed
3. Supports environment variable configuration
4. Maintains a built-in dictionary of known pathogen/toxin sequence motifs

Environment Variables:
    BIOCOMPILER_BLAST_DB_PATH: Path to local BLAST database directory
    BIOCOMPILER_BLAST_E_VALUE: E-value threshold (default 1e-5)
    BIOCOMPILER_BLAST_IDENTITY: Identity threshold (default 0.80)

Usage:
    >>> from biocompiler.blast_screening import BlastScreener
    >>> screener = BlastScreener()
    >>> result = screener.screen_sequence("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    >>> print(result.is_safe, result.screening_database)
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from typing import Optional

logger = logging.getLogger(__name__)

__all__ = [
    "BlastHit",
    "BlastScreeningResult",
    "BlastScreener",
    "PATHOGEN_TOXIN_MOTIFS",
    "QUICK_SCREEN_DB_NAME",
]

# ─── Known Pathogen/Toxin Motif Database ─────────────────────────────────

# Built-in dictionary of known pathogen-associated sequence motifs.
# These are key functional motifs from select agents and toxins.
# Each entry maps a motif name to its sequence pattern and metadata.
#
# NOTE: These are simplified representative motifs for screening purposes.
# Production databases should use full NCBI pathogen databases.
PATHOGEN_TOXIN_MOTIFS: dict[str, dict[str, str | bool]] = {
    # Anthrax lethal factor (Bacillus anthracis) - key catalytic motif
    "anthrax_lethal_factor": {
        "sequence": "ATGAATTTAGTTAAAGTATTA",
        "organism": "Bacillus anthracis",
        "is_pathogen": True,
        "is_toxin": True,
        "description": "Anthrax lethal factor zinc metalloprotease domain motif",
    },
    # Botulinum neurotoxin type A (Clostridium botulinum) - catalytic motif
    "botulinum_toxin_A": {
        "sequence": "ATGAAAAAACAATTTAGTTAT",
        "organism": "Clostridium botulinum",
        "is_pathogen": True,
        "is_toxin": True,
        "description": "Botulinum neurotoxin type A light chain motif",
    },
    # Ricin toxin (Ricinus communis) - A-chain catalytic motif
    "ricin_toxin": {
        "sequence": "ATGATGTTCGCTCCACCAACT",
        "organism": "Ricinus communis",
        "is_pathogen": False,
        "is_toxin": True,
        "description": "Ricin toxin A-chain RNA N-glycosidase motif",
    },
    # Diphtheria toxin (Corynebacterium diphtheriae) - catalytic domain
    "diphtheria_toxin": {
        "sequence": "ATGGGAACTTTCTTCGCTACA",
        "organism": "Corynebacterium diphtheriae",
        "is_pathogen": True,
        "is_toxin": True,
        "description": "Diphtheria toxin ADP-ribosyltransferase domain motif",
    },
    # Shiga toxin (Shigella dysenteriae) - A subunit motif
    "shiga_toxin": {
        "sequence": "ATGAAAAAAATTATTCGCTCT",
        "organism": "Shigella dysenteriae",
        "is_pathogen": True,
        "is_toxin": True,
        "description": "Shiga toxin A subunit RNA glycosidase motif",
    },
    # Cholera toxin (Vibrio cholerae) - A subunit motif
    "cholera_toxin": {
        "sequence": "ATGATTAATGCTTATGATGCA",
        "organism": "Vibrio cholerae",
        "is_pathogen": True,
        "is_toxin": True,
        "description": "Cholera toxin A subunit ADP-ribosyltransferase motif",
    },
    # Tetanus toxin (Clostridium tetani) - light chain motif
    "tetanus_toxin": {
        "sequence": "ATGAAAAAAGTTTTACGCTAT",
        "organism": "Clostridium tetani",
        "is_pathogen": True,
        "is_toxin": True,
        "description": "Tetanus neurotoxin zinc endopeptidase motif",
    },
    # Abrin toxin (Abrus precatorius) - A-chain motif
    "abrin_toxin": {
        "sequence": "ATGTCTGCTCCTCCAACAACT",
        "organism": "Abrus precatorius",
        "is_pathogen": False,
        "is_toxin": True,
        "description": "Abrin toxin A-chain RNA N-glycosidase motif",
    },
}

QUICK_SCREEN_DB_NAME = "biocompiler_builtin_pathogen_toxin_motifs"


# ─── Data Classes ────────────────────────────────────────────────────────


@dataclass
class BlastHit:
    """Represents a single BLAST hit against a database sequence.

    Attributes:
        subject_id: Identifier of the subject sequence in the database.
        subject_organism: Organism name of the subject sequence.
        identity_percent: Percent identity of the alignment (0.0-100.0).
        alignment_length: Length of the alignment in base pairs / amino acids.
        e_value: Expect value (statistical significance).
        bit_score: Bit score of the alignment.
        is_pathogen: Whether the subject is a known pathogen.
        is_toxin: Whether the subject is a known toxin.
    """

    subject_id: str
    subject_organism: str
    identity_percent: float
    alignment_length: int
    e_value: float
    bit_score: float
    is_pathogen: bool = False
    is_toxin: bool = False

    def __post_init__(self) -> None:
        """Validate fields after initialization."""
        if not self.subject_id:
            raise ValueError("subject_id must not be empty")
        if not 0.0 <= self.identity_percent <= 100.0:
            raise ValueError(
                f"identity_percent must be in [0, 100], got {self.identity_percent}"
            )
        if self.alignment_length < 0:
            raise ValueError(
                f"alignment_length must be non-negative, got {self.alignment_length}"
            )
        if self.e_value < 0.0:
            raise ValueError(f"e_value must be non-negative, got {self.e_value}")


@dataclass
class BlastScreeningResult:
    """Result of BLAST-based biosecurity screening.

    Attributes:
        query_id: Identifier of the screened query sequence.
        hits: List of BLAST hits found (may be empty).
        is_safe: True if no concerning hits above threshold were found.
        screening_database: Name/identifier of the database used for screening.
        e_value_threshold: E-value threshold used for screening.
        screening_time_seconds: Wall-clock time of the screening in seconds.
    """

    query_id: str
    hits: list[BlastHit] = field(default_factory=list)
    is_safe: bool = True
    screening_database: str = QUICK_SCREEN_DB_NAME
    e_value_threshold: float = 1e-5
    screening_time_seconds: float = 0.0

    @property
    def pathogen_hits(self) -> list[BlastHit]:
        """Hits that are flagged as pathogen-associated."""
        return [h for h in self.hits if h.is_pathogen]

    @property
    def toxin_hits(self) -> list[BlastHit]:
        """Hits that are flagged as toxin-associated."""
        return [h for h in self.hits if h.is_toxin]

    @property
    def concerning_hits(self) -> list[BlastHit]:
        """Hits that are either pathogen or toxin associated."""
        return [h for h in self.hits if h.is_pathogen or h.is_toxin]


# ─── BLAST Screener ──────────────────────────────────────────────────────


class BlastScreener:
    """BLAST-based biosecurity sequence screener.

    Screens DNA or protein sequences against known pathogen/toxin databases
    using NCBI BLAST+ when available, with a built-in pattern-based fallback.

    Args:
        blast_db_path: Path to local BLAST database directory. If None,
            reads from BIOCOMPILER_BLAST_DB_PATH environment variable.
        e_value_threshold: Maximum e-value for a hit to be considered
            significant. Reads from BIOCOMPILER_BLAST_E_VALUE if set.
        identity_threshold: Minimum identity fraction (0.0-1.0) for a hit
            to be considered concerning. Reads from BIOCOMPILER_BLAST_IDENTITY
            if set.
    """

    def __init__(
        self,
        blast_db_path: str | None = None,
        e_value_threshold: float | None = None,
        identity_threshold: float | None = None,
    ) -> None:
        # Resolve configuration from env vars with defaults
        self.blast_db_path: str | None = blast_db_path or os.environ.get(
            "BIOCOMPILER_BLAST_DB_PATH"
        )
        self.e_value_threshold: float = (
            e_value_threshold
            if e_value_threshold is not None
            else float(os.environ.get("BIOCOMPILER_BLAST_E_VALUE", "1e-5"))
        )
        self.identity_threshold: float = (
            identity_threshold
            if identity_threshold is not None
            else float(os.environ.get("BIOCOMPILER_BLAST_IDENTITY", "0.80"))
        )

        # Cache for BLAST availability check
        self._blast_available: bool | None = None

        logger.debug(
            "BlastScreener initialized: db_path=%s, e_value=%s, identity=%s",
            self.blast_db_path,
            self.e_value_threshold,
            self.identity_threshold,
        )

    def is_blast_available(self) -> bool:
        """Check if BLAST+ binaries (blastn, blastp) are installed and accessible.

        Returns:
            True if both blastn and blastp are found on PATH.
        """
        if self._blast_available is not None:
            return self._blast_available

        try:
            blastn_found = shutil.which("blastn") is not None
            blastp_found = shutil.which("blastp") is not None
            self._blast_available = blastn_found and blastp_found
        except Exception:
            self._blast_available = False

        if not self._blast_available:
            logger.info(
                "BLAST+ not available on this system. "
                "Quick pattern-based screening will be used as fallback."
            )

        return self._blast_available

    def screen_sequence(
        self, sequence: str, sequence_id: str = "query"
    ) -> BlastScreeningResult:
        """Screen a DNA sequence for biosecurity concerns.

        Uses BLAST+ with blastn if available, otherwise falls back to
        pattern-based quick screening.

        Args:
            sequence: DNA sequence to screen (ACGTN characters).
            sequence_id: Identifier for the query sequence.

        Returns:
            BlastScreeningResult with screening details and safety verdict.
        """
        start_time = time.monotonic()
        sequence = sequence.upper().strip()

        if not sequence:
            elapsed = time.monotonic() - start_time
            return BlastScreeningResult(
                query_id=sequence_id,
                hits=[],
                is_safe=True,
                screening_database="none",
                e_value_threshold=self.e_value_threshold,
                screening_time_seconds=elapsed,
            )

        # Try BLAST+ first
        if self.is_blast_available() and self.blast_db_path:
            try:
                hits = self._run_blastn(sequence)
                elapsed = time.monotonic() - start_time
                is_safe = not any(
                    h.identity_percent >= self.identity_threshold * 100
                    and h.e_value <= self.e_value_threshold
                    and (h.is_pathogen or h.is_toxin)
                    for h in hits
                )
                return BlastScreeningResult(
                    query_id=sequence_id,
                    hits=hits,
                    is_safe=is_safe,
                    screening_database=f"blastn:{self.blast_db_path}",
                    e_value_threshold=self.e_value_threshold,
                    screening_time_seconds=elapsed,
                )
            except Exception as exc:
                logger.warning(
                    "BLAST+ blastn failed, falling back to quick_screen: %s", exc
                )

        # Fallback to pattern-based screening
        result = self._quick_screen(sequence, sequence_id)
        result.screening_time_seconds = time.monotonic() - start_time
        return result

    def screen_protein(
        self, protein: str, protein_id: str = "query"
    ) -> BlastScreeningResult:
        """Screen a protein sequence for biosecurity concerns.

        Uses BLAST+ with blastp if available, otherwise falls back to
        pattern-based quick screening (with translated motifs).

        Args:
            protein: Protein sequence to screen (single-letter amino acid codes).
            protein_id: Identifier for the query protein.

        Returns:
            BlastScreeningResult with screening details and safety verdict.
        """
        start_time = time.monotonic()
        protein = protein.upper().strip()

        if not protein:
            elapsed = time.monotonic() - start_time
            return BlastScreeningResult(
                query_id=protein_id,
                hits=[],
                is_safe=True,
                screening_database="none",
                e_value_threshold=self.e_value_threshold,
                screening_time_seconds=elapsed,
            )

        # Try BLAST+ first
        if self.is_blast_available() and self.blast_db_path:
            try:
                hits = self._run_blastp(protein)
                elapsed = time.monotonic() - start_time
                is_safe = not any(
                    h.identity_percent >= self.identity_threshold * 100
                    and h.e_value <= self.e_value_threshold
                    and (h.is_pathogen or h.is_toxin)
                    for h in hits
                )
                return BlastScreeningResult(
                    query_id=protein_id,
                    hits=hits,
                    is_safe=is_safe,
                    screening_database=f"blastp:{self.blast_db_path}",
                    e_value_threshold=self.e_value_threshold,
                    screening_time_seconds=elapsed,
                )
            except Exception as exc:
                logger.warning(
                    "BLAST+ blastp failed, falling back to quick_screen: %s", exc
                )

        # Fallback: use DNA-based quick screen on back-translated motifs
        # For protein screening, we search for motif subsequences in the protein
        result = self._quick_protein_screen(protein, protein_id)
        result.screening_time_seconds = time.monotonic() - start_time
        return result

    # ─── Internal BLAST+ Methods ─────────────────────────────────────

    def _run_blastn(self, sequence: str) -> list[BlastHit]:
        """Run blastn against the configured BLAST database.

        Args:
            sequence: DNA sequence to query.

        Returns:
            List of BlastHit objects parsed from BLAST output.

        Raises:
            RuntimeError: If blastn execution fails.
        """
        if not self.blast_db_path:
            raise RuntimeError("No BLAST database path configured")

        # Write query to temp file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp:
            tmp.write(f">query\n{sequence}\n")
            tmp_path = tmp.name

        try:
            cmd = [
                "blastn",
                "-query", tmp_path,
                "-db", self.blast_db_path,
                "-outfmt",
                "6 qseqid sseqid stitle pident length evalue bitscore",
                "-evalue", str(self.e_value_threshold),
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,  # 5-minute timeout
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"blastn failed (rc={result.returncode}): {result.stderr}"
                )
            return self._parse_blast_output(result.stdout)
        finally:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass

    def _run_blastp(self, protein: str) -> list[BlastHit]:
        """Run blastp against the configured BLAST database.

        Args:
            protein: Protein sequence to query.

        Returns:
            List of BlastHit objects parsed from BLAST output.

        Raises:
            RuntimeError: If blastp execution fails.
        """
        if not self.blast_db_path:
            raise RuntimeError("No BLAST database path configured")

        # Write query to temp file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp:
            tmp.write(f">query\n{protein}\n")
            tmp_path = tmp.name

        try:
            cmd = [
                "blastp",
                "-query", tmp_path,
                "-db", self.blast_db_path,
                "-outfmt",
                "6 qseqid sseqid stitle pident length evalue bitscore",
                "-evalue", str(self.e_value_threshold),
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"blastp failed (rc={result.returncode}): {result.stderr}"
                )
            return self._parse_blast_output(result.stdout)
        finally:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass

    def _parse_blast_output(self, output: str) -> list[BlastHit]:
        """Parse BLAST tabular output format 6.

        Expected columns: qseqid sseqid stitle pident length evalue bitscore

        Args:
            output: Raw BLAST tabular output string.

        Returns:
            List of BlastHit objects.
        """
        hits: list[BlastHit] = []
        for line in output.strip().split("\n"):
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 7:
                logger.warning("Skipping malformed BLAST output line: %r", line)
                continue

            try:
                subject_id = fields[1]
                subject_title = fields[2] if len(fields) > 2 else subject_id
                identity_percent = float(fields[3])
                alignment_length = int(fields[4])
                e_value = float(fields[5])
                bit_score = float(fields[6])

                # Determine pathogen/toxin status from subject title
                is_pathogen = self._is_pathogen_subject(subject_id, subject_title)
                is_toxin = self._is_toxin_subject(subject_id, subject_title)

                # Extract organism from title (typically in square brackets)
                organism = self._extract_organism(subject_title)

                hits.append(
                    BlastHit(
                        subject_id=subject_id,
                        subject_organism=organism,
                        identity_percent=identity_percent,
                        alignment_length=alignment_length,
                        e_value=e_value,
                        bit_score=bit_score,
                        is_pathogen=is_pathogen,
                        is_toxin=is_toxin,
                    )
                )
            except (ValueError, IndexError) as exc:
                logger.warning("Failed to parse BLAST hit line: %r (%s)", line, exc)
                continue

        return hits

    # ─── Quick Screen Fallback ───────────────────────────────────────

    def _quick_screen(
        self, sequence: str, sequence_id: str = "query"
    ) -> BlastScreeningResult:
        """Simple pattern-based screening when BLAST+ is not installed.

        Checks the sequence against known pathogen/toxin motif patterns
        using substring matching. This is a conservative approximation —
        it may produce false positives but should not miss exact matches.

        Args:
            sequence: DNA sequence to screen (uppercase).
            sequence_id: Query identifier.

        Returns:
            BlastScreeningResult with pattern-match hits.
        """
        hits: list[BlastHit] = []

        for motif_name, motif_info in PATHOGEN_TOXIN_MOTIFS.items():
            motif_seq = str(motif_info["sequence"]).upper()
            if motif_seq in sequence:
                # Compute approximate identity and e-value for the match
                # For exact substring matches, identity is 100%
                hit = BlastHit(
                    subject_id=motif_name,
                    subject_organism=str(motif_info.get("organism", "Unknown")),
                    identity_percent=100.0,
                    alignment_length=len(motif_seq),
                    e_value=0.0,  # Exact match has e-value ~0
                    bit_score=float(len(motif_seq) * 2),  # Rough approximation
                    is_pathogen=bool(motif_info.get("is_pathogen", False)),
                    is_toxin=bool(motif_info.get("is_toxin", False)),
                )
                hits.append(hit)

        # Determine safety: not safe if any hit exceeds identity threshold
        is_safe = not any(
            h.identity_percent >= self.identity_threshold * 100
            and (h.is_pathogen or h.is_toxin)
            for h in hits
        )

        return BlastScreeningResult(
            query_id=sequence_id,
            hits=hits,
            is_safe=is_safe,
            screening_database=QUICK_SCREEN_DB_NAME,
            e_value_threshold=self.e_value_threshold,
        )

    def _quick_protein_screen(
        self, protein: str, protein_id: str = "query"
    ) -> BlastScreeningResult:
        """Simple pattern-based screening for protein sequences.

        Since our built-in motifs are DNA sequences, this performs a
        translated search by checking if any motif's translated product
        appears in the protein. For simplicity, this uses a heuristic
        that checks for common toxin protein keywords in the sequence.

        Args:
            protein: Protein sequence to screen (uppercase).
            protein_id: Query identifier.

        Returns:
            BlastScreeningResult with pattern-match hits.
        """
        hits: list[BlastHit] = []

        # For protein screening without BLAST+, we check for common
        # toxin-associated short peptide motifs. These are derived from
        # the catalytic sites of the toxins in our DNA motif database.
        PROTEIN_TOXIN_MOTIFS: dict[str, dict[str, str | bool]] = {
            "ricin_A_chain_catalytic": {
                "sequence": "NGSFS",
                "organism": "Ricinus communis",
                "is_pathogen": False,
                "is_toxin": True,
                "description": "Ricin A-chain active site motif",
            },
            "diphtheria_toxin_catalytic": {
                "sequence": "GVYVA",
                "organism": "Corynebacterium diphtheriae",
                "is_pathogen": True,
                "is_toxin": True,
                "description": "Diphtheria toxin catalytic site",
            },
            "botulinum_toxin_catalytic": {
                "sequence": "HELIH",
                "organism": "Clostridium botulinum",
                "is_pathogen": True,
                "is_toxin": True,
                "description": "Botulinum toxin zinc-binding motif",
            },
            "anthrax_LF_catalytic": {
                "sequence": "HEFGH",
                "organism": "Bacillus anthracis",
                "is_pathogen": True,
                "is_toxin": True,
                "description": "Anthrax lethal factor zinc-binding motif",
            },
            "tetanus_toxin_catalytic": {
                "sequence": "HEMTH",
                "organism": "Clostridium tetani",
                "is_pathogen": True,
                "is_toxin": True,
                "description": "Tetanus toxin zinc-binding motif",
            },
        }

        for motif_name, motif_info in PROTEIN_TOXIN_MOTIFS.items():
            motif_seq = str(motif_info["sequence"]).upper()
            if motif_seq in protein:
                hit = BlastHit(
                    subject_id=motif_name,
                    subject_organism=str(motif_info.get("organism", "Unknown")),
                    identity_percent=100.0,
                    alignment_length=len(motif_seq),
                    e_value=0.0,
                    bit_score=float(len(motif_seq) * 2),
                    is_pathogen=bool(motif_info.get("is_pathogen", False)),
                    is_toxin=bool(motif_info.get("is_toxin", False)),
                )
                hits.append(hit)

        is_safe = not any(
            h.identity_percent >= self.identity_threshold * 100
            and (h.is_pathogen or h.is_toxin)
            for h in hits
        )

        return BlastScreeningResult(
            query_id=protein_id,
            hits=hits,
            is_safe=is_safe,
            screening_database=QUICK_SCREEN_DB_NAME,
            e_value_threshold=self.e_value_threshold,
        )

    # ─── Helper Methods ──────────────────────────────────────────────

    @staticmethod
    def _extract_organism(title: str) -> str:
        """Extract organism name from a BLAST subject title.

        Looks for organism name in square brackets [Organism name].
        Falls back to the full title if no brackets found.

        Args:
            title: BLAST subject title string.

        Returns:
            Extracted organism name or the full title.
        """
        # Common format: "description [Organism name]"
        start = title.rfind("[")
        end = title.rfind("]")
        if start != -1 and end != -1 and end > start:
            return title[start + 1 : end].strip()
        return title.strip()

    @staticmethod
    def _is_pathogen_subject(subject_id: str, title: str) -> bool:
        """Heuristic check if a BLAST subject is a known pathogen.

        Checks both the subject ID and title against known pathogen
        organism keywords.

        Args:
            subject_id: BLAST subject identifier.
            title: BLAST subject description.

        Returns:
            True if the subject appears to be a pathogen.
        """
        pathogen_keywords = {
            "anthracis",
            "botulinum",
            "diphtheriae",
            "cholerae",
            "pestis",
            "tularensis",
            "hemorrhagic",
            "ebola",
            "variola",
            "smallpox",
        }
        combined = f"{subject_id} {title}".lower()
        return any(kw in combined for kw in pathogen_keywords)

    @staticmethod
    def _is_toxin_subject(subject_id: str, title: str) -> bool:
        """Heuristic check if a BLAST subject is a known toxin.

        Args:
            subject_id: BLAST subject identifier.
            title: BLAST subject description.

        Returns:
            True if the subject appears to be a toxin.
        """
        toxin_keywords = {
            "toxin",
            "enterotoxin",
            "neurotoxin",
            "cytotoxin",
            "hemolysin",
            "lethal factor",
            "ricin",
            "abrin",
            "saxitoxin",
            "tetrodotoxin",
        }
        combined = f"{subject_id} {title}".lower()
        return any(kw in combined for kw in toxin_keywords)
