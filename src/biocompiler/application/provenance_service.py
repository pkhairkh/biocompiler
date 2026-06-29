"""
BioCompiler Application — Provenance storage and retrieval service.

Handles persistence and querying of optimization decision trails.
"""

import logging
import os
from typing import Optional

from biocompiler.provenance.decision_provenance import ProvenanceStore, OptimizationDecisionTrail

logger = logging.getLogger(__name__)

# ─── Persistent Provenance Store ─────────────────────────────────────

_PROVENANCE_DIR = os.environ.get("BIOCOMPILER_PROVENANCE_DIR")
# Lazy: don't create the store (and ~/.biocompiler/provenance/) at import time.
_provenance_store = None

def _get_provenance_store():
    global _provenance_store
    if _provenance_store is None:
        _provenance_store = ProvenanceStore(store_dir=_PROVENANCE_DIR)
    return _provenance_store


def get_provenance_store() -> ProvenanceStore:
    """Return the global provenance store instance."""
    return _get_provenance_store()


def store_provenance(trail: OptimizationDecisionTrail) -> str:
    """Persist a decision-level provenance trail and return its record ID.

    Uses the :class:`ProvenanceStore` from :mod:`decision_provenance`
    for file-based persistence so records survive process restarts.
    """
    record_id = _get_provenance_store().save(trail)
    return record_id


def retrieve_provenance(record_id: str) -> dict:
    """Retrieve a stored provenance record by ID.

    Returns:
        Dict with 'id' and 'trail' keys.

    Raises:
        FileNotFoundError: If record_id does not exist.
        ValueError: If record_id is malformed.
    """
    trail = _get_provenance_store().load(record_id)
    return {
        "id": record_id,
        "trail": trail.to_dict(),
    }


def query_provenance(
    protein_name: Optional[str] = None,
    organism: Optional[str] = None,
    date_range: Optional[tuple[str, str]] = None,
) -> dict:
    """Query/list stored provenance records.

    Supports filtering by gene/protein name, organism, and date range.
    All filters are optional and combined with AND logic.

    Returns:
        Dict with 'count' and 'records' keys.
    """
    trails = _get_provenance_store().query(
        protein_name=protein_name,
        organism=organism,
        date_range=date_range,
    )
    return {
        "count": len(trails),
        "records": [
            {
                "gene_name": t.gene_name,
                "organism": t.organism,
                "timestamp": t.timestamp,
                "total_cai": t.total_cai,
                "total_gc": t.total_gc,
                "codon_decision_count": len(t.codon_decisions),
                "constraint_decision_count": len(t.constraint_decisions),
            }
            for t in trails
        ],
    }
