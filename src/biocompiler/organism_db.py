"""
BioCompiler Organism Database — SQLite + Kazusa API Integration

Production-grade organism data management with:
- Local SQLite database for fast offline codon usage lookups
- Kazusa Codon Usage Database API integration for on-demand data retrieval
- Automatic caching of API results in SQLite
- Migration from hardcoded dicts to persistent storage
- Backward-compatible interface with existing organisms module

The database stores:
- Organism metadata (name, taxonomy, source)
- Codon usage frequencies (per-organism, per-codon)
- Codon adaptiveness values (computed from frequencies)
- Preferred codons per amino acid
- Cache timestamps for API-fetched data
"""

import json
import logging
import sqlite3
import urllib.request
import urllib.error
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from .constants import CODON_TABLE, AA_TO_CODONS
from .exceptions import UnsupportedOrganismError

logger = logging.getLogger(__name__)

# Default database path (in user's home directory for persistence)
DEFAULT_DB_PATH = Path.home() / ".biocompiler" / "organisms.db"

# Kazusa API endpoint for codon usage tables
KAZUSA_API_URL = "https://www.kazusa.or.jp/codon/cgi-bin/spacialize.cgi"

# Known Kazusa database IDs for common organisms
KAZUSA_ORGANISM_IDS: dict[str, str] = {
    "Homo_sapiens": "9606",        # NCBI TaxID
    "Mus_musculus": "10090",
    "Escherichia_coli": "511145",  # K-12 MG1655
    "Saccharomyces_cerevisiae": "4932",
    "CHO_K1": "10129",            # Cricetulus griseus
    "Drosophila_melanogaster": "7227",
    "Caenorhabditis_elegans": "6239",
    "Danio_rerio": "7955",
    "Arabidopsis_thaliana": "3702",
    "Pichia_pastoris": "4922",
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

    def __init__(self, db_path: Optional[Path] = None):
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

    def _ensure_schema(self):
        """Create database tables if they don't exist."""
        conn = self._connect()
        try:
            conn.executescript("""
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
            conn.commit()
        finally:
            conn.close()

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
            except Exception as e:
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
                per_thousand = freq * 1000.0
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

    def fetch_from_kazusa(self, organism: str) -> None:
        """
        Fetch codon usage data from the Kazusa Codon Usage Database.

        The Kazusa database (https://www.kazusa.or.jp/codon/) provides
        codon usage tables computed from GenBank/RefSeq CDS entries for
        thousands of organisms. This method fetches the data via their
        web API, parses it, and stores it in the local SQLite database.

        Args:
            organism: Organism name or NCBI taxonomy ID

        Raises:
            ValueError: if the organism cannot be found in Kazusa
            urllib.error.URLError: if the API is unreachable
        """
        tax_id = KAZUSA_ORGANISM_IDS.get(organism, organism)

        logger.info("Fetching codon usage for %s (tax_id=%s) from Kazusa API", organism, tax_id)

        url = f"{KAZUSA_API_URL}?species={tax_id}&aa=1&style=N"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "BioCompiler/2.2.0"})
            with urllib.request.urlopen(req, timeout=30) as response:
                html = response.read().decode("utf-8", errors="replace")
        except urllib.error.URLError as e:
            logger.error("Kazusa API unreachable: %s", e)
            raise

        # Parse the HTML response to extract codon usage data
        codon_usage = self._parse_kazusa_html(html, organism)
        if not codon_usage:
            raise ValueError(f"No codon usage data found for {organism} (tax_id={tax_id})")

        # Store in database
        n_cds = self._extract_n_cds_from_html(html)
        self.store_organism(
            name=organism,
            codon_usage=codon_usage,
            taxonomy_id=tax_id,
            source="kazusa",
            n_cds=n_cds,
        )
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
        from .organisms import (
            CODON_USAGE_TABLES, CODON_ADAPTIVENESS_TABLES,
            PREFERRED_CODON_TABLES, SUPPORTED_ORGANISMS,
        )

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
    def _parse_kazusa_html(html: str, organism: str) -> dict[str, tuple[str, float, float, int]]:
        """
        Parse Kazusa codon usage HTML response into our internal format.

        The Kazusa HTML contains a table with codon, amino acid, frequency,
        and count per thousand. We extract and convert this data.
        """
        codon_usage: dict[str, tuple[str, float, float, int]] = {}
        codons = [c for cs in AA_TO_CODONS.values() for c in cs]

        # Kazusa format: rows like "GGG G 1.00 28frequency 0.36"
        # We parse by finding codon patterns in the HTML
        for codon in codons:
            aa = CODON_TABLE.get(codon, "X")
            if aa == "*":
                aa = "STOP"

            # Search for the codon in the HTML and extract nearby numbers
            idx = html.find(codon)
            if idx >= 0:
                # Extract numbers after the codon
                import re
                # Look for frequency (decimal) and count (integer) near the codon
                fragment = html[idx:idx + 100]
                numbers = re.findall(r'[\d.]+', fragment)
                freq = float(numbers[1]) if len(numbers) > 1 else 0.0
                count = int(float(numbers[2])) if len(numbers) > 2 else 0
                # Adaptiveness is computed later
                codon_usage[codon] = (aa, freq, 0.0, count)

        # If parsing failed, fall back to uniform distribution
        if len(codon_usage) < 64:
            logger.warning("Incomplete Kazusa parsing for %s (%d/64 codons), using uniform fallback",
                         organism, len(codon_usage))
            for codon in codons:
                if codon not in codon_usage:
                    aa = CODON_TABLE.get(codon, "X")
                    if aa == "*":
                        aa = "STOP"
                    codon_usage[codon] = (aa, 1.0 / len(AA_TO_CODONS.get(aa, [codon])), 0.0, 0)

        return codon_usage

    @staticmethod
    def _extract_n_cds_from_html(html: str) -> int:
        """Extract the number of CDS sequences from Kazusa HTML."""
        import re
        match = re.search(r'(\d+)\s+CDS', html)
        if match:
            return int(match.group(1))
        return 0


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
