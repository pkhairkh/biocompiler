"""
BioCompiler Organism Database — SQLite + Kazusa API Integration

Production-grade organism data management with:
- Local SQLite database for fast offline codon usage lookups
- Kazusa Codon Usage Database API integration for on-demand data retrieval
- Automatic caching of API results in SQLite with TTL
- Retry logic with exponential backoff for network resilience
- Robust HTML parser with multiple parsing strategies (incl. BeautifulSoup)
- JSON API fallback for Kazusa data
- Response validation and integrity checks
- Migration from hardcoded dicts to persistent storage
- Backward-compatible interface with existing organisms module
- Database schema versioning with migration support

The database stores:
- Organism metadata (name, taxonomy, source)
- Codon usage frequencies (per-organism, per-codon)
- Codon adaptiveness values (computed from frequencies)
- Preferred codons per amino acid
- Cache timestamps and TTL for API-fetched data
- Schema version for migration tracking
"""

import hashlib
import logging
import re
import sqlite3
import time
import urllib.request
import urllib.error
from collections.abc import Generator
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from ..constants import CODON_TABLE, AA_TO_CODONS
from ..exceptions import UnsupportedOrganismError

logger = logging.getLogger(__name__)

__all__ = [
    "SCHEMA_VERSION",
    "DEFAULT_DB_PATH",
    "KAZUSA_API_URL",
    "CACHE_TTL_SECONDS",
    "MAX_RETRIES",
    "RETRY_BACKOFF_BASE",
    "REQUEST_TIMEOUT",
    "MIN_CODON_COUNT",
    "TOTAL_CODONS",
    "NUM_STOP_CODONS",
    "MAX_FREQUENCY_PER_CODON",
    "MIN_TOTAL_FREQUENCY",
    "PER_THOUSAND_SCALE",
    "KAZUSA_ORGANISM_IDS",
    "ORGANISM_NAME_ALIASES",
    "OrganismDatabase",
    "get_database",
    "get_codon_usage_db",
    "get_codon_adaptiveness_db",
    "resolve_organism_name",
    "get_codon_table",
]

# ─── Database Schema Versioning ─────────────────────────────────────

SCHEMA_VERSION = 2

# Default database path (in user's home directory for persistence)
DEFAULT_DB_PATH = Path.home() / ".biocompiler" / "organisms.db"

# Kazusa API endpoint for codon usage tables
KAZUSA_API_URL = "https://www.kazusa.or.jp/codon/cgi-bin/spacialize.cgi"

# Cache TTL for API-fetched data (7 days by default)
CACHE_TTL_SECONDS = int(7 * 24 * 60 * 60)  # 1 week

# Network retry configuration
MAX_RETRIES = 3
RETRY_BACKOFF_BASE = 2  # seconds
REQUEST_TIMEOUT = 30  # seconds

# HTTP User-Agent header for API requests
_USER_AGENT = "BioCompiler/4.0.0"

# Response validation thresholds
MIN_CODON_COUNT = 60  # Minimum codons expected from Kazusa (out of 64)
TOTAL_CODONS = 64  # Total number of codons in the standard genetic code
NUM_STOP_CODONS = 3  # Number of stop codons (TAA, TAG, TGA)
MAX_FREQUENCY_PER_CODON = 1.0  # Each codon frequency must be <= 1.0
MIN_TOTAL_FREQUENCY = 0.9  # Total frequency per AA should be close to 1.0

# Scale factor for per-thousand codon frequency representation
PER_THOUSAND_SCALE = 1000.0

# Known Kazusa database IDs for common organisms
KAZUSA_ORGANISM_IDS: dict[str, str] = {
    "Homo_sapiens": "9606",        # NCBI TaxID
    "Mus_musculus": "10090",
    "Escherichia_coli": "511145",  # K-12 MG1655
    "E_coli": "511145",            # Alias for Escherichia_coli
    "Saccharomyces_cerevisiae": "4932",
    "CHO_K1": "10129",            # Cricetulus griseus
    "Drosophila_melanogaster": "7227",
    "Caenorhabditis_elegans": "6239",
    "Danio_rerio": "7955",
    "Arabidopsis_thaliana": "3702",
    "Pichia_pastoris": "4922",
}

# ─── Organism Name Aliases ──────────────────────────────────────────
# Maps common shorthand names / legacy identifiers to the canonical
# database key.  Used by :func:`get_codon_table` and the
# :class:`OrganismDatabase` lookup methods so that callers can pass
# ``"ecoli"`` instead of ``"E_coli"``, etc.

ORGANISM_NAME_ALIASES: dict[str, str] = {
    # E. coli variants
    "ecoli": "E_coli",
    "e_coli": "E_coli",
    "e. coli": "E_coli",
    "e.coli": "E_coli",
    "E_coli_K12": "E_coli",
    "E_coli_BL21": "E_coli",
    "Escherichia_coli": "E_coli",
    "Escherichia_coli_K12": "E_coli",
    "Escherichia_coli_BL21": "E_coli",
    "escherichia coli": "E_coli",
    # Human
    "human": "Homo_sapiens",
    "homo_sapiens": "Homo_sapiens",
    "homo sapiens": "Homo_sapiens",
    "h. sapiens": "Homo_sapiens",
    # Mouse
    "mouse": "Mus_musculus",
    "mus_musculus": "Mus_musculus",
    "mus musculus": "Mus_musculus",
    # CHO
    "cho": "CHO_K1",
    "cho_k1": "CHO_K1",
    "Cricetulus_griseus": "CHO_K1",
    # Yeast
    "yeast": "Saccharomyces_cerevisiae",
    "s_cerevisiae": "Saccharomyces_cerevisiae",
    "s. cerevisiae": "Saccharomyces_cerevisiae",
    "saccharomyces cerevisiae": "Saccharomyces_cerevisiae",
    # Other organisms
    "drosophila": "Drosophila_melanogaster",
    "d_melanogaster": "Drosophila_melanogaster",
    "d. melanogaster": "Drosophila_melanogaster",
    "celegans": "Caenorhabditis_elegans",
    "c_elegans": "Caenorhabditis_elegans",
    "c. elegans": "Caenorhabditis_elegans",
    "zebrafish": "Danio_rerio",
    "d_rerio": "Danio_rerio",
    "arabidopsis": "Arabidopsis_thaliana",
    "a_thaliana": "Arabidopsis_thaliana",
    "pichia": "Pichia_pastoris",
    "p_pastoris": "Pichia_pastoris",
    "komagataella": "Pichia_pastoris",
}


class OrganismDatabase:
    """
    SQLite-backed organism codon usage database with Kazusa API fallback.

    This class provides a persistent, queryable store for organism-specific
    codon usage data. When a requested organism is not in the local database,
    it attempts to fetch the data from the Kazusa Codon Usage Database via
    their web API, caches the result locally, and returns the data.

    Usage:
        db = OrganismDatabase()
        usage = db.get_codon_usage("Homo_sapiens")
        adaptiveness = db.get_codon_adaptiveness("Homo_sapiens")
        preferred = db.get_preferred_codons("Homo_sapiens")

    Thread Safety:
        Each method creates its own connection, so the class is safe for
        concurrent reads. Write operations should be serialized externally.
    """

    def __init__(self, db_path: Optional[Path] = None) -> None:
        """
        Initialize the organism database.

        Args:
            db_path: Path to the SQLite database file.
                     Defaults to ~/.biocompiler/organisms.db
        """
        self.db_path = Path(db_path) if db_path else DEFAULT_DB_PATH
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._ensure_schema()

    def _connect(self) -> sqlite3.Connection:
        """Create a new database connection with row factory."""
        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        return conn

    def _connection(self) -> Generator[sqlite3.Connection, None, None]:
        """Context manager for database connections — auto-commits and closes."""
        conn = self._connect()
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()

    def _ensure_schema(self) -> None:
        """Create database tables if they don't exist, with schema versioning."""
        conn = self._connect()
        try:
            # Check current schema version
            version = self._get_schema_version(conn)

            if version < 1:
                # Initial schema (v1)
                conn.executescript("""
                    CREATE TABLE IF NOT EXISTS schema_version (
                        key TEXT PRIMARY KEY,
                        value TEXT NOT NULL
                    );

                    CREATE TABLE IF NOT EXISTS organisms (
                        name TEXT PRIMARY KEY,
                        taxonomy_id TEXT,
                        taxonomy_lineage TEXT,
                        source TEXT DEFAULT 'builtin',
                        n_cds INTEGER DEFAULT 0,
                        created_at TEXT NOT NULL,
                        updated_at TEXT NOT NULL
                    );

                    CREATE TABLE IF NOT EXISTS codon_usage (
                        organism TEXT NOT NULL,
                        codon TEXT NOT NULL,
                        amino_acid TEXT NOT NULL,
                        frequency REAL NOT NULL,
                        count INTEGER DEFAULT 0,
                        per_thousand REAL DEFAULT 0.0,
                        PRIMARY KEY (organism, codon),
                        FOREIGN KEY (organism) REFERENCES organisms(name) ON DELETE CASCADE
                    );

                    CREATE TABLE IF NOT EXISTS codon_adaptiveness (
                        organism TEXT NOT NULL,
                        codon TEXT NOT NULL,
                        adaptiveness REAL NOT NULL,
                        PRIMARY KEY (organism, codon),
                        FOREIGN KEY (organism) REFERENCES organisms(name) ON DELETE CASCADE
                    );

                    CREATE TABLE IF NOT EXISTS preferred_codons (
                        organism TEXT NOT NULL,
                        amino_acid TEXT NOT NULL,
                        codon TEXT NOT NULL,
                        PRIMARY KEY (organism, amino_acid),
                        FOREIGN KEY (organism) REFERENCES organisms(name) ON DELETE CASCADE
                    );

                    CREATE INDEX IF NOT EXISTS idx_codon_usage_organism
                        ON codon_usage(organism);
                    CREATE INDEX IF NOT EXISTS idx_codon_adaptiveness_organism
                        ON codon_adaptiveness(organism);
                    CREATE INDEX IF NOT EXISTS idx_preferred_codons_organism
                        ON preferred_codons(organism);
                """)

            if version < 2:
                # v2: Add response hash for cache integrity validation
                conn.executescript("""
                    CREATE TABLE IF NOT EXISTS api_response_cache (
                        organism TEXT PRIMARY KEY,
                        response_hash TEXT NOT NULL,
                        fetched_at TEXT NOT NULL,
                        url TEXT NOT NULL,
                        FOREIGN KEY (organism) REFERENCES organisms(name) ON DELETE CASCADE
                    );
                """)

            # Update schema version
            conn.execute(
                "INSERT OR REPLACE INTO schema_version (key, value) VALUES (?, ?)",
                ("version", str(SCHEMA_VERSION)),
            )
            conn.commit()
        finally:
            conn.close()

    @staticmethod
    def _get_schema_version(conn: sqlite3.Connection) -> int:
        """Get the current schema version from the database."""
        try:
            row = conn.execute(
                "SELECT value FROM schema_version WHERE key = 'version'"
            ).fetchone()
            return int(row["value"]) if row else 0
        except sqlite3.OperationalError:
            # Table doesn't exist yet — old schema without versioning
            # Check if organisms table exists to determine version
            try:
                conn.execute("SELECT 1 FROM organisms LIMIT 1")
                return 1  # Has old schema but no version table
            except sqlite3.OperationalError:
                return 0  # Fresh database

    # ─── Read Operations ──────────────────────────────────────────

    def list_organisms(self) -> list[dict]:
        """
        List all organisms in the database.

        Returns:
            List of dicts with keys: name, taxonomy_id, source, n_cds, updated_at
        """
        conn = self._connect()
        try:
            rows = conn.execute(
                "SELECT name, taxonomy_id, source, n_cds, updated_at FROM organisms ORDER BY name"
            ).fetchall()
            return [dict(r) for r in rows]
        finally:
            conn.close()

    def organism_exists(self, organism: str) -> bool:
        """Check if an organism is in the database."""
        conn = self._connect()
        try:
            row = conn.execute(
                "SELECT 1 FROM organisms WHERE name = ?", (organism,)
            ).fetchone()
            return row is not None
        finally:
            conn.close()

    def get_codon_usage(self, organism: str) -> dict[str, tuple[str, float, float, int]]:
        """
        Get codon usage table for an organism.

        Returns:
            Dict mapping codon -> (amino_acid, frequency, adaptiveness, count)
            Same format as the organisms module CODON_USAGE_TABLES.

        Raises:
            UnsupportedOrganismError: if organism not in database and not fetchable
        """
        conn = self._connect()
        try:
            rows = conn.execute(
                "SELECT codon, amino_acid, frequency, count FROM codon_usage WHERE organism = ?",
                (organism,)
            ).fetchall()
        finally:
            conn.close()

        if not rows:
            # Try to fetch from Kazusa API
            try:
                self.fetch_from_kazusa(organism)
                # Retry after fetch
                conn = self._connect()
                try:
                    rows = conn.execute(
                        "SELECT codon, amino_acid, frequency, count FROM codon_usage WHERE organism = ?",
                        (organism,)
                    ).fetchall()
                finally:
                    conn.close()
            except (ValueError, urllib.error.URLError, sqlite3.Error) as e:
                logger.warning("Failed to fetch organism %s from Kazusa: %s", organism, e)
                raise UnsupportedOrganismError(organism, self.list_organism_names())

        result = {}
        for r in rows:
            result[r["codon"]] = (r["amino_acid"], r["frequency"], 0.0, r["count"])

        # Add adaptiveness values
        conn = self._connect()
        try:
            adapt_rows = conn.execute(
                "SELECT codon, adaptiveness FROM codon_adaptiveness WHERE organism = ?",
                (organism,)
            ).fetchall()
        finally:
            conn.close()

        for r in adapt_rows:
            if r["codon"] in result:
                aa, freq, _, count = result[r["codon"]]
                result[r["codon"]] = (aa, freq, r["adaptiveness"], count)

        return result

    def get_codon_adaptiveness(self, organism: str) -> dict[str, float]:
        """Get codon adaptiveness values for an organism."""
        # Ensure data exists
        if not self._has_adaptiveness(organism):
            self._compute_and_store_adaptiveness(organism)

        conn = self._connect()
        try:
            rows = conn.execute(
                "SELECT codon, adaptiveness FROM codon_adaptiveness WHERE organism = ?",
                (organism,)
            ).fetchall()
            return {r["codon"]: r["adaptiveness"] for r in rows}
        finally:
            conn.close()

    def get_preferred_codons(self, organism: str) -> dict[str, str]:
        """Get preferred codons per amino acid for an organism."""
        if not self._has_preferred(organism):
            self._compute_and_store_preferred(organism)

        conn = self._connect()
        try:
            rows = conn.execute(
                "SELECT amino_acid, codon FROM preferred_codons WHERE organism = ?",
                (organism,)
            ).fetchall()
            return {r["amino_acid"]: r["codon"] for r in rows}
        finally:
            conn.close()

    def list_organism_names(self) -> list[str]:
        """List all organism names in the database."""
        conn = self._connect()
        try:
            rows = conn.execute("SELECT name FROM organisms ORDER BY name").fetchall()
            return [r["name"] for r in rows]
        finally:
            conn.close()

    # ─── Write Operations ─────────────────────────────────────────

    def store_organism(
        self,
        name: str,
        codon_usage: dict[str, tuple[str, float, float, int]],
        taxonomy_id: Optional[str] = None,
        taxonomy_lineage: Optional[str] = None,
        source: str = "builtin",
        n_cds: int = 0,
    ) -> None:
        """
        Store an organism's codon usage data in the database.

        This replaces the hardcoded Python dicts with persistent SQLite storage.
        If the organism already exists, its data is updated.

        Args:
            name: Organism name (e.g., "Homo_sapiens")
            codon_usage: Dict mapping codon -> (aa, frequency, adaptiveness, count)
            taxonomy_id: NCBI taxonomy ID
            taxonomy_lineage: Taxonomy lineage string
            source: Data source ("builtin", "kazusa", "custom")
            n_cds: Number of CDS sequences used for computing frequencies
        """
        now = datetime.now(timezone.utc).isoformat()
        conn = self._connect()
        try:
            # Upsert organism metadata
            conn.execute(
                """INSERT INTO organisms (name, taxonomy_id, taxonomy_lineage, source, n_cds, created_at, updated_at)
                   VALUES (?, ?, ?, ?, ?, ?, ?)
                   ON CONFLICT(name) DO UPDATE SET
                       taxonomy_id = excluded.taxonomy_id,
                       taxonomy_lineage = excluded.taxonomy_lineage,
                       source = excluded.source,
                       n_cds = excluded.n_cds,
                       updated_at = excluded.updated_at""",
                (name, taxonomy_id, taxonomy_lineage, source, n_cds, now, now),
            )

            # Store codon usage data
            for codon, (aa, freq, adapt, count) in codon_usage.items():
                per_thousand = freq * PER_THOUSAND_SCALE
                conn.execute(
                    """INSERT INTO codon_usage (organism, codon, amino_acid, frequency, count, per_thousand)
                       VALUES (?, ?, ?, ?, ?, ?)
                       ON CONFLICT(organism, codon) DO UPDATE SET
                           amino_acid = excluded.amino_acid,
                           frequency = excluded.frequency,
                           count = excluded.count,
                           per_thousand = excluded.per_thousand""",
                    (name, codon, aa, freq, count, per_thousand),
                )

                conn.execute(
                    """INSERT INTO codon_adaptiveness (organism, codon, adaptiveness)
                       VALUES (?, ?, ?)
                       ON CONFLICT(organism, codon) DO UPDATE SET
                           adaptiveness = excluded.adaptiveness""",
                    (name, codon, adapt),
                )

            conn.commit()
            logger.info("Stored organism %s (%d codons, source=%s)", name, len(codon_usage), source)
        finally:
            conn.close()

        # Compute and store preferred codons
        self._compute_and_store_preferred(name)
        self._compute_and_store_adaptiveness(name)

    def fetch_from_kazusa(self, organism: str, force: bool = False) -> None:
        """
        Fetch codon usage data from the Kazusa Codon Usage Database.

        The Kazusa database (https://www.kazusa.or.jp/codon/) provides
        codon usage tables computed from GenBank/RefSeq CDS entries for
        thousands of organisms. This method fetches the data via their
        web API, parses it, and stores it in the local SQLite database.

        Features:
        - Retry with exponential backoff (3 attempts)
        - Cache validation with TTL (7 days)
        - Multiple HTML parsing strategies for robustness
        - Graceful fallback to uniform distribution on failure

        Args:
            organism: Organism name or NCBI taxonomy ID
            force: Force re-fetch even if cache is fresh

        Raises:
            ValueError: if the organism cannot be found in Kazusa
            urllib.error.URLError: if the API is unreachable after retries
        """
        # Check cache freshness (skip fetch if recently fetched)
        if not force and self._is_cache_fresh(organism):
            logger.debug("Using cached data for %s (cache is fresh)", organism)
            return

        tax_id = KAZUSA_ORGANISM_IDS.get(organism, organism)
        logger.info("Fetching codon usage for %s (tax_id=%s) from Kazusa API", organism, tax_id)

        url = f"{KAZUSA_API_URL}?species={tax_id}&aa=1&style=N"

        # Retry with exponential backoff
        html: Optional[str] = None
        last_error: Optional[urllib.error.URLError] = None
        for attempt in range(MAX_RETRIES):
            try:
                req = urllib.request.Request(
                    url,
                    headers={"User-Agent": _USER_AGENT},
                )
                with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as response:
                    html = response.read().decode("utf-8", errors="replace")
                break  # Success
            except urllib.error.URLError as e:
                last_error = e
                if attempt < MAX_RETRIES - 1:
                    wait_time = RETRY_BACKOFF_BASE ** attempt
                    logger.warning(
                        "Kazusa API attempt %d/%d failed: %s. Retrying in %ds",
                        attempt + 1, MAX_RETRIES, e, wait_time,
                    )
                    time.sleep(wait_time)
                else:
                    logger.error("Kazusa API unreachable after %d attempts: %s", MAX_RETRIES, e)

        if html is None:
            raise last_error or urllib.error.URLError("Kazusa API unreachable")

        # Validate HTML response before parsing
        if not self._validate_kazusa_response(html, organism):
            logger.warning("Kazusa response validation failed for %s — may be malformed", organism)

        # Parse the HTML response using multiple strategies
        codon_usage = self._parse_kazusa_html_robust(html, organism)
        if not codon_usage:
            raise ValueError(f"No codon usage data found for {organism} (tax_id={tax_id})")

        # Validate parsed data integrity
        codon_usage = self._validate_codon_usage_data(codon_usage, organism)

        # Store in database
        n_cds = self._extract_n_cds_from_html(html)
        self.store_organism(
            name=organism,
            codon_usage=codon_usage,
            taxonomy_id=tax_id,
            source="kazusa",
            n_cds=n_cds,
        )

        # Store response hash for cache integrity
        self._store_response_hash(organism, html, url)

        logger.info("Successfully fetched and stored %s from Kazusa (%d codons)", organism, len(codon_usage))

    # ─── Migration from Built-in Data ─────────────────────────────

    def migrate_builtin_data(self) -> int:
        """
        Migrate all built-in organism data from Python dicts to SQLite.

        This one-time operation transfers the hardcoded codon usage tables
        from the organisms module into the persistent database, enabling
        SQL queries, versioning, and future updates without code changes.

        Returns:
            Number of organisms migrated
        """
        from ..organisms import CODON_USAGE_TABLES, SUPPORTED_ORGANISMS

        migrated = 0
        for org_name in SUPPORTED_ORGANISMS:
            if self.organism_exists(org_name):
                logger.debug("Organism %s already in database, skipping", org_name)
                continue

            usage = CODON_USAGE_TABLES.get(org_name, {})
            if not usage:
                logger.warning("No usage data for %s, skipping", org_name)
                continue

            self.store_organism(
                name=org_name,
                codon_usage=usage,
                taxonomy_id=KAZUSA_ORGANISM_IDS.get(org_name),
                source="builtin",
                n_cds=0,
            )
            migrated += 1
            logger.info("Migrated %s to database", org_name)

        logger.info("Migration complete: %d organisms migrated", migrated)
        return migrated

    # ─── Private Methods ──────────────────────────────────────────

    def _has_adaptiveness(self, organism: str) -> bool:
        """Check whether codon adaptiveness data exists for the given organism."""
        conn = self._connect()
        try:
            row = conn.execute(
                "SELECT COUNT(*) as cnt FROM codon_adaptiveness WHERE organism = ?",
                (organism,)
            ).fetchone()
            return row["cnt"] > 0
        finally:
            conn.close()

    def _has_preferred(self, organism: str) -> bool:
        """Check whether preferred codon data exists for the given organism."""
        conn = self._connect()
        try:
            row = conn.execute(
                "SELECT COUNT(*) as cnt FROM preferred_codons WHERE organism = ?",
                (organism,)
            ).fetchone()
            return row["cnt"] > 0
        finally:
            conn.close()

    def _compute_and_store_adaptiveness(self, organism: str) -> None:
        """Compute codon adaptiveness values and store them."""
        usage = self._get_raw_usage(organism)
        if not usage:
            return

        # Compute adaptiveness: for each amino acid, the codon with the
        # highest frequency gets adaptiveness 1.0, and other codons are
        # scaled proportionally
        adaptiveness: dict[str, float] = {}
        for aa in AA_TO_CODONS:
            codons = AA_TO_CODONS[aa]
            max_freq = max(usage.get(c, 0.0) for c in codons)
            if max_freq > 0:
                for c in codons:
                    adaptiveness[c] = usage.get(c, 0.0) / max_freq
            else:
                for c in codons:
                    adaptiveness[c] = 0.0

        conn = self._connect()
        try:
            for codon, adapt in adaptiveness.items():
                conn.execute(
                    """INSERT INTO codon_adaptiveness (organism, codon, adaptiveness)
                       VALUES (?, ?, ?)
                       ON CONFLICT(organism, codon) DO UPDATE SET
                           adaptiveness = excluded.adaptiveness""",
                    (organism, codon, adapt),
                )
            conn.commit()
        finally:
            conn.close()

    def _compute_and_store_preferred(self, organism: str) -> None:
        """Compute preferred codons and store them."""
        adaptiveness = self.get_codon_adaptiveness(organism)
        if not adaptiveness:
            return

        preferred: dict[str, str] = {}
        for aa in AA_TO_CODONS:
            codons = AA_TO_CODONS[aa]
            best_codon = max(codons, key=lambda c: adaptiveness.get(c, 0.0))
            preferred[aa] = best_codon

        conn = self._connect()
        try:
            for aa, codon in preferred.items():
                conn.execute(
                    """INSERT INTO preferred_codons (organism, amino_acid, codon)
                       VALUES (?, ?, ?)
                       ON CONFLICT(organism, amino_acid) DO UPDATE SET
                           codon = excluded.codon""",
                    (organism, aa, codon),
                )
            conn.commit()
        finally:
            conn.close()

    def _get_raw_usage(self, organism: str) -> dict[str, float]:
        """Get raw codon frequencies for an organism."""
        conn = self._connect()
        try:
            rows = conn.execute(
                "SELECT codon, frequency FROM codon_usage WHERE organism = ?",
                (organism,)
            ).fetchall()
            return {r["codon"]: r["frequency"] for r in rows}
        finally:
            conn.close()

    @staticmethod
    def _parse_kazusa_html_robust(html: str, organism: str) -> dict[str, tuple[str, float, float, int]]:
        """
        Parse Kazusa codon usage HTML response using multiple parsing strategies.

        Strategy 0: BeautifulSoup parsing (most robust, requires optional dependency)
        Strategy 1: Structured table parsing (most reliable when format is standard)
        Strategy 2: Fragment-based extraction (original approach, improved)
        Strategy 3: Regex-based extraction (fallback)

        Falls back to uniform distribution only if ALL strategies fail to
        extract a minimum of MIN_CODON_COUNT/TOTAL_CODONS codons.
        """
        codon_usage: dict[str, tuple[str, float, float, int]] = {}
        # Use all codons from CODON_TABLE (includes stop codons)
        codons = list(CODON_TABLE.keys())

        # Strategy 0: BeautifulSoup (most robust, requires optional dependency)
        try:
            from bs4 import BeautifulSoup  # type: ignore[import-untyped]

            soup = BeautifulSoup(html, 'html.parser')
            # Find all table rows in the codon usage table
            for tr in soup.find_all('tr'):
                cells = tr.find_all('td')
                if len(cells) < 4:
                    continue
                # Extract text content from each cell
                cell_texts = [c.get_text(strip=True) for c in cells]
                # Find a cell that looks like a codon (3-letter ACGT sequence)
                codon_val: Optional[str] = None
                codon_idx: Optional[int] = None
                for i, txt in enumerate(cell_texts):
                    if len(txt) == 3 and all(ch in 'ACGTacgt' for ch in txt):
                        codon_val = txt.upper()
                        codon_idx = i
                        break
                if codon_val is None or codon_val not in CODON_TABLE:
                    continue
                # After the codon cell, expect: amino acid, frequency, count, per_thousand
                remaining = cell_texts[codon_idx + 1:]
                if len(remaining) < 2:
                    continue
                # Identify the amino acid (1-letter or *)
                aa = remaining[0]
                if aa == '*':
                    aa = 'STOP'
                # Identify numeric fields: frequency (float), count (int), per_thousand (float)
                # Use text representation to distinguish integer count from float metrics:
                #   - Integer text (no decimal point) → likely count
                #   - Float text with decimal point ≤ 1.0 → likely frequency
                #   - Float text with decimal point > 1.0 → likely per-thousand
                freq = 0.0
                count = 0
                per_thousand = 0.0
                for txt in remaining[1:]:
                    # Skip non-numeric cells (e.g., full AA name "Gly")
                    if not txt or not txt[0].isdigit():
                        continue
                    try:
                        val = float(txt)
                    except ValueError:
                        continue
                    if '.' not in txt:
                        # Integer text — treat as count
                        if count == 0:
                            count = int(val)
                    elif val <= 1.0:
                        if freq == 0.0:
                            freq = val
                    else:
                        if per_thousand == 0.0:
                            per_thousand = val
                # Derive frequency from per_thousand if not found directly
                if freq == 0.0 and per_thousand > 0.0:
                    freq = per_thousand / PER_THOUSAND_SCALE
                # Derive frequency from count if neither freq nor per_thousand was found
                if freq == 0.0 and count > 0:
                    freq = count / PER_THOUSAND_SCALE
                # Validate amino acid against CODON_TABLE
                expected_aa = CODON_TABLE.get(codon_val)
                if expected_aa:
                    expected_aa = 'STOP' if expected_aa == '*' else expected_aa
                if expected_aa and aa != expected_aa:
                    aa = expected_aa
                codon_usage[codon_val] = (aa, freq, 0.0, count)

            if len(codon_usage) >= MIN_CODON_COUNT:
                logger.info(
                    "Kazusa parsing Strategy 0 (BeautifulSoup) succeeded: %d/%d codons",
                    len(codon_usage), TOTAL_CODONS,
                )
            else:
                codon_usage.clear()
        except ImportError:
            logger.debug("BeautifulSoup4 not installed, falling back to regex parsing")
        except Exception as exc:
            logger.debug("BeautifulSoup parsing failed, falling back to regex: %s", exc)
            codon_usage.clear()

        # If Strategy 0 succeeded, skip regex strategies
        if len(codon_usage) >= MIN_CODON_COUNT:
            # Final fallback: uniform distribution for missing codons
            if len(codon_usage) < TOTAL_CODONS:
                missing = TOTAL_CODONS - len(codon_usage)
                logger.warning(
                    "Incomplete Kazusa parsing for %s (%d/%d codons, %d missing), "
                    "using uniform fallback for missing codons",
                    organism, len(codon_usage), TOTAL_CODONS, missing,
                )
                for codon in codons:
                    if codon not in codon_usage:
                        aa = CODON_TABLE.get(codon, "X")
                        if aa == "*":
                            aa = "STOP"
                        if aa == "STOP":
                            codon_usage[codon] = (aa, 1.0 / NUM_STOP_CODONS, 0.0, 0)
                        else:
                            codon_usage[codon] = (aa, 1.0 / len(AA_TO_CODONS.get(aa, [codon])), 0.0, 0)
            return codon_usage

        # Strategy 1: Parse structured HTML table rows
        # Kazusa format typically has rows like:
        # <td>GGG</td><td>G</td><td>1.00</td><td>28</td><td>0.36</td>
        row_pattern = re.compile(
            r'<td[^>]*>([ACGT]{3})</td>\s*'
            r'<td[^>]*>([A-Z*]+)</td>\s*'
            r'<td[^>]*>([\d.]+)</td>\s*'
            r'<td[^>]*>(\d+)</td>\s*'
            r'<td[^>]*>([\d.]+)</td>',
            re.IGNORECASE,
        )
        for match in row_pattern.finditer(html):
            codon = match.group(1).upper()
            aa = match.group(2)
            freq = float(match.group(4)) / PER_THOUSAND_SCALE if int(match.group(4)) > 0 else float(match.group(3))
            count = int(match.group(4))
            if codon in CODON_TABLE:
                if aa == "*":
                    aa = "STOP"
                codon_usage[codon] = (aa, freq, 0.0, count)

        # If Strategy 1 found most codons, use it
        if len(codon_usage) >= MIN_CODON_COUNT:
            logger.info("Kazusa parsing Strategy 1 (table rows) succeeded: %d/%d codons", len(codon_usage), TOTAL_CODONS)
        else:
            codon_usage.clear()

            # Strategy 2: Fragment-based extraction (improved)
            for codon in codons:
                aa = CODON_TABLE.get(codon, "X")
                if aa == "*":
                    aa = "STOP"

                idx = html.find(codon)
                if idx >= 0:
                    fragment = html[idx:idx + 120]
                    # Extract all numbers from the fragment
                    numbers = re.findall(r'[\d.]+', fragment)
                    if len(numbers) >= 3:
                        try:
                            # Typically: [codon_number, frequency, count, per_thousand]
                            # Try to identify the per-thousand value (should be <= 1000)
                            freq = 0.0
                            count = 0
                            for num_str in numbers[1:]:
                                val = float(num_str)
                                if val > PER_THOUSAND_SCALE and count == 0:
                                    count = int(val)
                                elif val <= PER_THOUSAND_SCALE and freq == 0.0:
                                    freq = val / PER_THOUSAND_SCALE if val > 1.0 else val
                            codon_usage[codon] = (aa, freq, 0.0, count)
                        except (ValueError, IndexError) as exc:
                            logger.debug("Failed to parse codon %s from fragment: %s", codon, exc)

            if len(codon_usage) >= MIN_CODON_COUNT:
                logger.info("Kazusa parsing Strategy 2 (fragment) succeeded: %d/%d codons", len(codon_usage), TOTAL_CODONS)
            else:
                codon_usage.clear()

                # Strategy 3: Regex-based codon-frequency pairs
                # Look for patterns like "GGG G 0.28" or codon followed by numbers
                for codon in codons:
                    aa = CODON_TABLE.get(codon, "X")
                    if aa == "*":
                        aa = "STOP"
                    # Find codon followed by optional amino acid and number
                    pattern = re.compile(
                        rf'{codon}\s+[A-Z*]?\s*([\d.]+)\s*(\d+)?',
                        re.IGNORECASE,
                    )
                    match = pattern.search(html)
                    if match:
                        freq = float(match.group(1))
                        count = int(match.group(2) or 0)
                        freq = freq / PER_THOUSAND_SCALE if freq > 1.0 else freq
                        codon_usage[codon] = (aa, freq, 0.0, count)

                if len(codon_usage) >= MIN_CODON_COUNT:
                    logger.info("Kazusa parsing Strategy 3 (regex) succeeded: %d/%d codons", len(codon_usage), TOTAL_CODONS)

        # Final fallback: uniform distribution for missing codons
        if len(codon_usage) < TOTAL_CODONS:
            missing = TOTAL_CODONS - len(codon_usage)
            logger.warning(
                "Incomplete Kazusa parsing for %s (%d/%d codons, %d missing), "
                "using uniform fallback for missing codons",
                organism, len(codon_usage), TOTAL_CODONS, missing,
            )
            for codon in codons:
                if codon not in codon_usage:
                    aa = CODON_TABLE.get(codon, "X")
                    if aa == "*":
                        aa = "STOP"
                    # For stop codons, use equal frequency among the stop codons
                    if aa == "STOP":
                        codon_usage[codon] = (aa, 1.0 / NUM_STOP_CODONS, 0.0, 0)
                    else:
                        codon_usage[codon] = (aa, 1.0 / len(AA_TO_CODONS.get(aa, [codon])), 0.0, 0)

        return codon_usage

    @staticmethod
    def _extract_n_cds_from_html(html: str) -> int:
        """Extract the number of CDS sequences from Kazusa HTML."""
        match = re.search(r'(\d+)\s+CDS', html)
        if match:
            return int(match.group(1))
        return 0

    def _is_cache_fresh(self, organism: str) -> bool:
        """Check if the cached data for an organism is still fresh (within TTL)."""
        conn = self._connect()
        try:
            row = conn.execute(
                "SELECT updated_at FROM organisms WHERE name = ? AND source = 'kazusa'",
                (organism,),
            ).fetchone()
            if not row:
                return False
            try:
                updated = datetime.fromisoformat(row["updated_at"])
                age = (datetime.now(timezone.utc) - updated).total_seconds()
                return age < CACHE_TTL_SECONDS
            except (ValueError, TypeError):
                return False
        finally:
            conn.close()

    def _validate_kazusa_response(self, html: str, organism: str) -> bool:
        """
        Validate that a Kazusa HTML response contains expected structural elements.

        Checks for:
        - Non-empty response
        - Contains codon-like patterns (3-letter nucleotide sequences)
        - Contains numeric data (frequencies/counts)
        - Does not contain error messages from the Kazusa server

        When BeautifulSoup4 is available, uses it to extract clean text content
        for error-message detection, avoiding false positives from HTML tags,
        attributes, or embedded JavaScript/CSS.
        """
        if not html or len(html) < 100:
            logger.warning("Kazusa response too short for %s (%d bytes)", organism, len(html) if html else 0)
            return False

        # Check for Kazusa error messages
        error_indicators = [
            "not found", "no data", "error", "invalid species",
            "no entries", "does not exist",
        ]

        # Use BeautifulSoup for more robust text extraction when available
        text_to_check: str
        try:
            from bs4 import BeautifulSoup  # type: ignore[import-untyped]
            soup = BeautifulSoup(html, 'html.parser')
            # Extract visible text only (strips tags, scripts, styles)
            for script_or_style in soup(['script', 'style']):
                script_or_style.decompose()
            text_to_check = soup.get_text(separator=' ', strip=True).lower()
        except ImportError:
            # Fallback to raw HTML lowercased (original behavior)
            text_to_check = html.lower()

        for indicator in error_indicators:
            if indicator in text_to_check:
                logger.warning("Kazusa response contains error indicator '%s' for %s", indicator, organism)
                return False

        # Check for at least some codon patterns
        codon_pattern = re.compile(r'[ACGT]{3}')
        codon_matches = codon_pattern.findall(html)
        if len(codon_matches) < 20:
            logger.warning(
                "Kazusa response has few codon patterns for %s (%d found)",
                organism, len(codon_matches),
            )
            return False

        return True

    @staticmethod
    def _validate_codon_usage_data(
        codon_usage: dict[str, tuple[str, float, float, int]],
        organism: str,
    ) -> dict[str, tuple[str, float, float, int]]:
        """
        Validate parsed codon usage data for internal consistency.

        Checks:
        - All 64 codons are present
        - Frequencies are non-negative
        - Amino acid assignments match the standard codon table
        - Per-amino-acid frequency totals are reasonable

        Fixes:
        - Corrects amino acid assignments that don't match CODON_TABLE
        - Normalizes frequencies that exceed 1.0
        """
        validated: dict[str, tuple[str, float, float, int]] = {}
        corrections = 0

        for codon, (aa, freq, adapt, count) in codon_usage.items():
            # Validate amino acid assignment
            expected_aa = CODON_TABLE.get(codon)
            if expected_aa and expected_aa == "*":
                expected_aa = "STOP"
            if expected_aa and aa != expected_aa:
                logger.debug(
                    "Correcting AA for codon %s in %s: %s → %s",
                    codon, organism, aa, expected_aa,
                )
                aa = expected_aa
                corrections += 1

            # Validate frequency range
            if freq < 0:
                logger.warning("Negative frequency for %s in %s: %f — clamping to 0", codon, organism, freq)
                freq = 0.0
                corrections += 1
            if freq > MAX_FREQUENCY_PER_CODON:
                logger.debug("Frequency > 1.0 for %s in %s: %f — normalizing", codon, organism, freq)
                freq = min(freq, MAX_FREQUENCY_PER_CODON)
                corrections += 1

            validated[codon] = (aa, freq, adapt, count)

        if corrections:
            logger.info("Made %d corrections to codon usage data for %s", corrections, organism)

        # Check total codon count
        if len(validated) < MIN_CODON_COUNT:
            logger.warning(
                "Only %d/%d codons parsed for %s (minimum expected: %d)",
                len(validated), TOTAL_CODONS, organism, MIN_CODON_COUNT,
            )

        return validated

    def _store_response_hash(self, organism: str, html: str, url: str) -> None:
        """Store a hash of the API response for cache integrity validation."""
        response_hash = hashlib.sha256(html.encode("utf-8")).hexdigest()
        now = datetime.now(timezone.utc).isoformat()
        conn = self._connect()
        try:
            conn.execute(
                """INSERT INTO api_response_cache (organism, response_hash, fetched_at, url)
                   VALUES (?, ?, ?, ?)
                   ON CONFLICT(organism) DO UPDATE SET
                       response_hash = excluded.response_hash,
                       fetched_at = excluded.fetched_at,
                       url = excluded.url""",
                (organism, response_hash, now, url),
            )
            conn.commit()
        except sqlite3.OperationalError as e:
            # Table may not exist in older schemas — log but don't fail
            logger.debug("Could not store response hash: %s", e)
        finally:
            conn.close()

    def verify_cache_integrity(self, organism: str) -> bool:
        """
        Verify that cached data for an organism matches the stored response hash.

        Returns True if the cache is intact, False if corrupted.
        """
        conn = self._connect()
        try:
            row = conn.execute(
                "SELECT response_hash, url FROM api_response_cache WHERE organism = ?",
                (organism,),
            ).fetchone()
            if not row:
                return True  # No hash stored — can't verify, assume OK
            stored_hash = row["response_hash"]
            url = row["url"]
        finally:
            conn.close()

        # Re-fetch and compare
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": _USER_AGENT},
            )
            with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as response:
                html = response.read().decode("utf-8", errors="replace")
            current_hash = hashlib.sha256(html.encode("utf-8")).hexdigest()
            return current_hash == stored_hash
        except urllib.error.URLError as e:
            logger.warning("Could not verify cache integrity for %s: %s", organism, e)
            return True  # Can't verify — assume OK


# ─── Module-level convenience functions ────────────────────────────

_db_instance: Optional[OrganismDatabase] = None


def get_database(db_path: Optional[Path] = None) -> OrganismDatabase:
    """Get or create the singleton database instance."""
    global _db_instance
    if _db_instance is None or (db_path and db_path != _db_instance.db_path):
        _db_instance = OrganismDatabase(db_path)
    return _db_instance


def get_codon_usage_db(organism: str) -> dict[str, tuple[str, float, float, int]]:
    """Convenience function to get codon usage from database."""
    return get_database().get_codon_usage(organism)


def get_codon_adaptiveness_db(organism: str) -> dict[str, float]:
    """Convenience function to get codon adaptiveness from database."""
    return get_database().get_codon_adaptiveness(organism)


def resolve_organism_name(organism: str) -> str:
    """Resolve an organism name or alias to its canonical database key.

    This handles common shorthand names (e.g. ``"ecoli"`` → ``"E_coli"``),
    legacy identifiers, and case-insensitive lookups.  If no alias is
    found the original name is returned unchanged.

    .. note::
       This function should be used consistently throughout the codebase
       whenever an organism name needs to be mapped to a canonical key.
       All public APIs that accept organism names (e.g. :func:`get_codon_table`,
       :class:`OrganismDatabase` methods) should call this function to
       normalize inputs before performing lookups.

    Args:
        organism: An organism name, alias, or shorthand.

    Returns:
        The canonical database key for the organism.

    Examples::

        >>> resolve_organism_name("ecoli")
        'E_coli'
        >>> resolve_organism_name("human")
        'Homo_sapiens'
        >>> resolve_organism_name("Homo_sapiens")
        'Homo_sapiens'
        >>> resolve_organism_name("E. coli")
        'E_coli'
        >>> resolve_organism_name("homo sapiens")
        'Homo_sapiens'
    """
    # Normalize: lowercase and strip whitespace for robust matching
    org_normalized = organism.strip().lower()

    # Direct hit on alias map (case-insensitive via normalization)
    for key, canonical in ORGANISM_NAME_ALIASES.items():
        if key.lower() == org_normalized:
            return canonical

    # Already a canonical name (e.g. in KAZUSA_ORGANISM_IDS keys)
    for key in KAZUSA_ORGANISM_IDS:
        if key.lower() == org_normalized:
            return key

    # No match found — return as-is (let downstream code handle the error)
    return organism


def get_codon_table(organism: str) -> dict[str, tuple[str, float, float, int]]:
    """Get the codon usage table for an organism, resolving aliases.

    This is the primary high-level lookup for codon usage data.  It
    accepts any of the common organism identifiers (``"ecoli"``,
    ``"E_coli"``, ``"Escherichia_coli"``, ``"human"``, etc.) and
    returns the codon usage table from the local SQLite database,
    fetching from Kazusa on-demand if necessary.

    Args:
        organism: Organism name, alias, or shorthand (e.g. ``"ecoli"``).

    Returns:
        Dict mapping codon → ``(amino_acid, frequency, adaptiveness, count)``.

    Raises:
        UnsupportedOrganismError: if the organism cannot be resolved or
            fetched.

    Examples::

        >>> table = get_codon_table("ecoli")       # resolves to E_coli
        >>> table = get_codon_table("human")        # resolves to Homo_sapiens
        >>> table = get_codon_table("Homo_sapiens") # direct lookup
    """
    canonical = resolve_organism_name(organism)
    if canonical != organism:
        logger.info("Resolved organism name %r -> %r", organism, canonical)
    return get_database().get_codon_usage(canonical)
