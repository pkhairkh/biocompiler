"""Type stubs for biocompiler.export — FASTA, GenBank, JSON, and SBOL export."""

from __future__ import annotations

from typing import Any, Optional, TypedDict

from .types import Certificate, TypeCheckResult, Verdict


# ────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────

GENBANK_MAX_LINE: int
GENBANK_SEQ_LINE: int
GENBANK_SEQ_GROUP: int


# ────────────────────────────────────────────────────────────
# TypedDicts
# ────────────────────────────────────────────────────────────

class _RestrictionSiteInfoRequired(TypedDict):
    enzyme: str
    site: str
    position: int


class RestrictionSiteInfo(_RestrictionSiteInfoRequired, total=False):
    strand: str


class _FastaSequenceEntryRequired(TypedDict):
    sequence: str


class FastaSequenceEntry(_FastaSequenceEntryRequired, total=False):
    id: str
    description: str
    organism: str
    protein: str
    cai: float
    gc: float


# ────────────────────────────────────────────────────────────
# Core export functions
# ────────────────────────────────────────────────────────────

def export_fasta(
    sequence: str,
    identifier: str = ...,
    description: str = ...,
    organism: str = ...,
    protein: Optional[str] = ...,
    cai: Optional[float] = ...,
    include_comments: bool = ...,
    type_results: Optional[list[TypeCheckResult]] = ...,
) -> str: ...


def export_genbank(
    sequence: str,
    locus_name: str = ...,
    definition: str = ...,
    organism: str = ...,
    molecule_type: str = ...,
    topology: str = ...,
    exon_boundaries: Optional[list[tuple[int, int]]] = ...,
    restriction_sites: Optional[list[RestrictionSiteInfo]] = ...,
    certificate: Optional[Certificate] = ...,
    type_results: Optional[list[TypeCheckResult]] = ...,
    gene_name: Optional[str] = ...,
    protein: Optional[str] = ...,
    cai: Optional[float] = ...,
    optimization_date: Optional[str] = ...,
) -> str: ...


def export_genbank_with_certificate(
    sequence: str,
    certificate: Certificate,
    organism: str = ...,
    gene_name: Optional[str] = ...,
    exon_boundaries: Optional[list[tuple[int, int]]] = ...,
) -> str: ...


def export_multi_fasta(
    sequences: list[FastaSequenceEntry],
) -> str: ...


def export_batch_fasta(
    results: list[dict],
    organism: str = ...,
) -> str: ...


def export_json(
    result: Any,
    indent: int = ...,
    include_certificate: bool = ...,
    include_provenance: bool = ...,
) -> str: ...


def export_full_construct(
    sequence: str,
    organism: str = ...,
    **kwargs: Any,
) -> str: ...


def export_with_annotations(
    sequence: str,
    organism: str = ...,
    **kwargs: Any,
) -> str: ...


def format_biosecurity_report(
    sequence: str,
    organism: str = ...,
    cai: Optional[float] = ...,
    gc: Optional[float] = ...,
    type_results: Optional[list[TypeCheckResult]] = ...,
) -> str: ...
