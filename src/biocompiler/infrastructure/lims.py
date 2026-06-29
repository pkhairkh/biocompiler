"""
BioCompiler LIMS Integration Hooks
====================================

Provides integration with Laboratory Information Management Systems (LIMS)
for submitting optimized designs, querying status, and exporting to
popular LIMS platforms (Benchling, LabGuru).

The module defines:
- ``LIMSIntegration``: Abstract base class for LIMS system integration.
- ``BenchlingExporter``: Export biocompiler results in Benchling-compatible JSON.
- ``LabGuruExporter``: Export biocompiler results in LabGuru-compatible format.

Usage::

    from biocompiler.infrastructure.lims import BenchlingExporter, LIMSIntegration

    exporter = BenchlingExporter()
    payload = exporter.export(result)

    # Or use the convenience functions:
    from biocompiler.infrastructure.lims import export_to_benchling, export_to_labguru
    benchling_data = export_to_benchling(result)
    labguru_data  = export_to_labguru(result)
"""

from __future__ import annotations

import json
import logging
import uuid
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Optional

from biocompiler.optimizer import OptimizationResult
from biocompiler.sequence.scanner import gc_content
from .. import __version__

logger = logging.getLogger(__name__)

__all__ = [
    "LIMSIntegration",
    "BenchlingExporter",
    "LabGuruExporter",
    "LIMSSubmissionRecord",
    "export_to_benchling",
    "export_to_labguru",
]


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class LIMSSubmissionRecord:
    """Record of a design submission to a LIMS system.

    Attributes:
        design_id: Unique identifier assigned by the LIMS.
        project_id: Project this design belongs to.
        status: Current status (e.g. 'submitted', 'in_review', 'approved', 'rejected').
        submitted_at: ISO 8601 timestamp of submission.
        lims_system: Name of the LIMS system (e.g. 'benchling', 'labguru').
        payload: The JSON payload that was sent.
    """
    design_id: str
    project_id: str
    status: str = "submitted"
    submitted_at: str = ""
    lims_system: str = ""
    payload: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if not self.submitted_at:
            self.submitted_at = datetime.now(timezone.utc).isoformat()


# ────────────────────────────────────────────────────────────
# LIMSIntegration base class
# ────────────────────────────────────────────────────────────

class LIMSIntegration(ABC):
    """Base class for LIMS system integration.

    Subclasses must implement the abstract methods to communicate
    with a specific LIMS platform. The base class provides common
    infrastructure for design submission, status tracking, and
    export formatting.
    """

    def __init__(self, base_url: str = "", api_token: str = ""):
        """Initialize LIMS integration.

        Args:
            base_url: Base URL of the LIMS API endpoint.
            api_token: Authentication token for the LIMS API.
        """
        self.base_url = base_url.rstrip("/") if base_url else ""
        self.api_token = api_token
        self._submission_cache: dict[str, LIMSSubmissionRecord] = {}

    @abstractmethod
    def submit_design(self, result: OptimizationResult, project_id: str) -> str:
        """Submit an optimized design to the LIMS.

        Args:
            result: The optimization result to submit.
            project_id: Project identifier in the LIMS.

        Returns:
            design_id assigned by the LIMS.
        """
        ...

    @abstractmethod
    def get_design_status(self, design_id: str) -> dict:
        """Query design status from LIMS.

        Args:
            design_id: The design identifier returned by submit_design.

        Returns:
            Dict with at least 'status' key and any platform-specific fields.
        """
        ...

    @abstractmethod
    def export_to_benchling(self, result: OptimizationResult) -> dict:
        """Format result for Benchling API import.

        Args:
            result: The optimization result to format.

        Returns:
            Dict in Benchling-compatible JSON format.
        """
        ...

    @abstractmethod
    def export_to_labguru(self, result: OptimizationResult) -> dict:
        """Format result for LabGuru API import.

        Args:
            result: The optimization result to format.

        Returns:
            Dict in LabGuru-compatible format.
        """
        ...

    def _generate_design_id(self, project_id: str) -> str:
        """Generate a unique design ID for a project.

        The format is ``BC_{project_id}_{uuid8}`` where uuid8 is the
        first 8 hex characters of a UUID4.
        """
        short_uuid = uuid.uuid4().hex[:8].upper()
        return f"BC_{project_id}_{short_uuid}"


# ────────────────────────────────────────────────────────────
# Benchling Exporter
# ────────────────────────────────────────────────────────────

class BenchlingExporter(LIMSIntegration):
    """Export biocompiler results in Benchling-compatible JSON format.

    Benchling API expects:
        - name: Sequence name
        - sequence: DNA sequence string
        - annotations: List of annotation dicts with start, end, type, notes
        - custom fields: Key-value pairs for metadata

    The export produces a payload compatible with the Benchling
    ``POST /api/v2/DNA-sequences`` endpoint.

    Reference: https://benchling.com/api/reference#/DNA-sequences
    """

    def __init__(self, base_url: str = "", api_token: str = "",
                 folder_id: str = ""):
        """Initialize Benchling exporter.

        Args:
            base_url: Benchling API base URL (e.g. 'https://benchling.com/api/v2').
            api_token: Benchling API token.
            folder_id: Default folder ID for sequence creation.
        """
        super().__init__(base_url, api_token)
        self.folder_id = folder_id

    def submit_design(self, result: OptimizationResult, project_id: str) -> str:
        """Submit an optimized design to Benchling.

        In a production environment, this would make an HTTP POST to
        the Benchling API. This implementation generates the payload
        and caches the submission record.

        Args:
            result: The optimization result to submit.
            project_id: Project identifier (mapped to Benchling folder_id).

        Returns:
            design_id assigned for this submission.
        """
        design_id = self._generate_design_id(project_id)
        payload = self.export_to_benchling(result)

        # Override folder_id with project_id mapping if provided
        if project_id and not payload.get("folderId"):
            payload["folderId"] = self.folder_id or project_id

        record = LIMSSubmissionRecord(
            design_id=design_id,
            project_id=project_id,
            lims_system="benchling",
            payload=payload,
        )
        self._submission_cache[design_id] = record
        logger.info("Submitted design %s to Benchling for project %s",
                     design_id, project_id)
        return design_id

    def get_design_status(self, design_id: str) -> dict:
        """Query design status from Benchling.

        In production, this would GET from the Benchling API.
        This implementation returns the cached submission record status.

        Args:
            design_id: Design identifier.

        Returns:
            Dict with 'status' and submission metadata.
        """
        record = self._submission_cache.get(design_id)
        if record is None:
            return {"status": "not_found", "design_id": design_id}
        return {
            "status": record.status,
            "design_id": record.design_id,
            "project_id": record.project_id,
            "submitted_at": record.submitted_at,
            "lims_system": record.lims_system,
        }

    def export_to_benchling(self, result: OptimizationResult) -> dict:
        """Format result for Benchling API import.

        Produces a JSON payload compatible with Benchling's DNA sequence
        creation endpoint, including annotations for CDS, gene, and
        regulatory features, plus custom fields for CAI, GC, and
        optimization metadata.

        Args:
            result: The optimization result to format.

        Returns:
            Benchling-compatible dict with name, sequence, annotations,
            customFields, and folderId.
        """
        organism_name = getattr(result, "organism_name", None) or "unknown"
        seq = result.sequence.upper()
        seq_len = len(seq)

        # Build annotations
        annotations = self._build_annotations(result, seq_len)

        # Build custom fields
        custom_fields = self._build_custom_fields(result, organism_name)

        payload = {
            "name": f"BioCompiler_{organism_name}_{seq_len}bp",
            "sequence": seq,
            "annotations": annotations,
            "customFields": custom_fields,
            "isCircular": False,
            "schemaId": None,  # User should set their own schema
        }

        if self.folder_id:
            payload["folderId"] = self.folder_id

        return payload

    def export_to_labguru(self, result: OptimizationResult) -> dict:
        """Format result for LabGuru API import.

        Delegates to LabGuruExporter for the actual formatting.

        Args:
            result: The optimization result to format.

        Returns:
            LabGuru-compatible dict.
        """
        exporter = LabGuruExporter(
            base_url=self.base_url,
            api_token=self.api_token,
        )
        return exporter.export_to_labguru(result)

    def _build_annotations(
        self, result: OptimizationResult, seq_len: int
    ) -> list[dict[str, Any]]:
        """Build Benchling-style feature annotations.

        Each annotation has: name, type, start, end, strand, color, notes.
        """
        annotations: list[dict[str, Any]] = []

        # Gene annotation
        annotations.append({
            "name": "designed_gene",
            "type": "gene",
            "start": 0,
            "end": seq_len,
            "strand": 1,
            "color": "#4CAF50",
            "notes": f"Designed by BioCompiler v{__version__}",
        })

        # CDS annotation
        if result.protein:
            annotations.append({
                "name": "CDS",
                "type": "CDS",
                "start": 0,
                "end": seq_len,
                "strand": 1,
                "color": "#2196F3",
                "notes": (
                    f"Codon-optimized CDS | "
                    f"CAI={result.cai:.4f} | "
                    f"Protein={len(result.protein)} aa"
                ),
            })

        # Regulatory annotation for codon optimization
        annotations.append({
            "name": "codon_optimization",
            "type": "regulatory",
            "start": 0,
            "end": seq_len,
            "strand": 1,
            "color": "#FF9800",
            "notes": (
                f"Codon-optimized by BioCompiler v{__version__} | "
                f"GC={result.gc_content:.4f}"
            ),
        })

        return annotations

    def _build_custom_fields(
        self, result: OptimizationResult, organism_name: str
    ) -> dict[str, Any]:
        """Build Benchling custom fields from optimization result."""
        custom_fields: dict[str, Any] = {
            "biocompiler_version": {"value": __version__},
            "organism": {"value": organism_name},
            "cai": {"value": round(result.cai, 6)},
            "gc_content": {"value": round(result.gc_content, 6)},
            "sequence_length": {"value": len(result.sequence)},
            "optimization_date": {
                "value": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            },
        }

        if result.protein:
            custom_fields["protein_length"] = {"value": len(result.protein)}

        if result.failed_predicates:
            custom_fields["failed_predicates"] = {
                "value": ", ".join(result.failed_predicates)
            }

        if result.satisfied_predicates:
            custom_fields["satisfied_predicates"] = {
                "value": ", ".join(result.satisfied_predicates)
            }

        if result.mutagenesis_applied:
            custom_fields["mutagenesis_applied"] = {"value": True}
            if result.aa_substitutions:
                subs_str = "; ".join(
                    f"{s.get('original', '?')}{s.get('position', '?')}"
                    f"{s.get('substitution', '?')}"
                    for s in result.aa_substitutions
                )
                custom_fields["aa_substitutions"] = {"value": subs_str}

        if result.codon_pair_bias is not None:
            custom_fields["codon_pair_bias"] = {
                "value": round(result.codon_pair_bias, 6)
            }

        return custom_fields


# ────────────────────────────────────────────────────────────
# LabGuru Exporter
# ────────────────────────────────────────────────────────────

class LabGuruExporter(LIMSIntegration):
    """Export biocompiler results in LabGuru-compatible format.

    LabGuru API expects:
        - item: { name, description, ... }
        - data: Base64-encoded or raw sequence
        - tags: List of tag strings for searchability
        - custom_fields: Key-value metadata

    The export produces a payload compatible with the LabGuru
    DNA sequence creation endpoint.

    Reference: https://labguru.com/api/docs
    """

    def __init__(self, base_url: str = "", api_token: str = "",
                 project_id: int = 0):
        """Initialize LabGuru exporter.

        Args:
            base_url: LabGuru API base URL.
            api_token: LabGuru API token.
            project_id: Default project ID in LabGuru.
        """
        super().__init__(base_url, api_token)
        self.project_id = project_id

    def submit_design(self, result: OptimizationResult, project_id: str) -> str:
        """Submit an optimized design to LabGuru.

        In a production environment, this would make an HTTP POST to
        the LabGuru API. This implementation generates the payload
        and caches the submission record.

        Args:
            result: The optimization result to submit.
            project_id: Project identifier (mapped to LabGuru project).

        Returns:
            design_id assigned for this submission.
        """
        design_id = self._generate_design_id(project_id)
        payload = self.export_to_labguru(result)

        # Set project if available
        if self.project_id:
            payload["item"]["project_id"] = self.project_id

        record = LIMSSubmissionRecord(
            design_id=design_id,
            project_id=project_id,
            lims_system="labguru",
            payload=payload,
        )
        self._submission_cache[design_id] = record
        logger.info("Submitted design %s to LabGuru for project %s",
                     design_id, project_id)
        return design_id

    def get_design_status(self, design_id: str) -> dict:
        """Query design status from LabGuru.

        In production, this would GET from the LabGuru API.
        This implementation returns the cached submission record status.

        Args:
            design_id: Design identifier.

        Returns:
            Dict with 'status' and submission metadata.
        """
        record = self._submission_cache.get(design_id)
        if record is None:
            return {"status": "not_found", "design_id": design_id}
        return {
            "status": record.status,
            "design_id": record.design_id,
            "project_id": record.project_id,
            "submitted_at": record.submitted_at,
            "lims_system": record.lims_system,
        }

    def export_to_benchling(self, result: OptimizationResult) -> dict:
        """Format result for Benchling API import.

        Delegates to BenchlingExporter for the actual formatting.

        Args:
            result: The optimization result to format.

        Returns:
            Benchling-compatible dict.
        """
        exporter = BenchlingExporter(
            base_url=self.base_url,
            api_token=self.api_token,
        )
        return exporter.export_to_benchling(result)

    def export_to_labguru(self, result: OptimizationResult) -> dict:
        """Format result for LabGuru API import.

        Produces a JSON payload compatible with LabGuru's DNA sequence
        creation endpoint, with structured item metadata, the raw
        sequence, searchable tags, and custom fields.

        Args:
            result: The optimization result to format.

        Returns:
            LabGuru-compatible dict with item, data, tags, and
            custom_fields keys.
        """
        organism_name = getattr(result, "organism_name", None) or "unknown"
        seq = result.sequence.upper()
        seq_len = len(seq)

        # Build item metadata
        item = {
            "name": f"BioCompiler_{organism_name}_{seq_len}bp",
            "description": self._build_description(result, organism_name),
            "type": "DNA Sequence",
        }
        if self.project_id:
            item["project_id"] = self.project_id

        # Build tags for searchability
        tags = self._build_tags(result, organism_name)

        # Build custom fields
        custom_fields = self._build_custom_fields(result, organism_name)

        payload = {
            "item": item,
            "data": seq,
            "tags": tags,
            "custom_fields": custom_fields,
        }

        return payload

    def _build_description(
        self, result: OptimizationResult, organism_name: str
    ) -> str:
        """Build a human-readable description for the LabGuru item."""
        desc_lines = [
            f"Codon-optimized DNA sequence designed by BioCompiler v{__version__}",
            f"Target organism: {organism_name.replace('_', ' ')}",
            f"Length: {len(result.sequence)} bp",
            f"CAI: {result.cai:.4f}",
            f"GC content: {result.gc_content:.4f}",
        ]
        if result.protein:
            desc_lines.append(f"Protein: {len(result.protein)} amino acids")
        if result.failed_predicates:
            desc_lines.append(
                f"Failed predicates: {', '.join(result.failed_predicates)}"
            )
        if result.mutagenesis_applied:
            desc_lines.append("Mutagenesis was applied during optimization")
        return " | ".join(desc_lines)

    def _build_tags(
        self, result: OptimizationResult, organism_name: str
    ) -> list[str]:
        """Build searchable tags for LabGuru."""
        tags = [
            "biocompiler",
            f"organism:{organism_name}",
            "codon-optimized",
            f"cai:{result.cai:.2f}",
            f"gc:{result.gc_content:.2f}",
        ]
        if result.mutagenesis_applied:
            tags.append("mutagenesis")
        if result.fallback_used:
            tags.append("fallback-used")
        if not result.failed_predicates:
            tags.append("all-predicates-pass")
        return tags

    def _build_custom_fields(
        self, result: OptimizationResult, organism_name: str
    ) -> dict[str, Any]:
        """Build LabGuru custom fields from optimization result."""
        custom_fields: dict[str, Any] = {
            "biocompiler_version": __version__,
            "organism": organism_name,
            "cai": result.cai,
            "gc_content": result.gc_content,
            "sequence_length": len(result.sequence),
            "optimization_date": datetime.now(timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%SZ"
            ),
        }

        if result.protein:
            custom_fields["protein_length"] = len(result.protein)
            custom_fields["protein_sequence"] = result.protein

        if result.codon_pair_bias is not None:
            custom_fields["codon_pair_bias"] = result.codon_pair_bias

        if result.suggested_5utr:
            custom_fields["suggested_5utr"] = result.suggested_5utr
        if result.suggested_3utr:
            custom_fields["suggested_3utr"] = result.suggested_3utr

        return custom_fields


# ────────────────────────────────────────────────────────────
# Convenience functions
# ────────────────────────────────────────────────────────────

def export_to_benchling(result: OptimizationResult) -> dict:
    """Export an OptimizationResult to Benchling-compatible JSON.

    Convenience wrapper around :class:`BenchlingExporter`.

    Args:
        result: The optimization result to export.

    Returns:
        Benchling-compatible dict.
    """
    exporter = BenchlingExporter()
    return exporter.export_to_benchling(result)


def export_to_labguru(result: OptimizationResult) -> dict:
    """Export an OptimizationResult to LabGuru-compatible format.

    Convenience wrapper around :class:`LabGuruExporter`.

    Args:
        result: The optimization result to export.

    Returns:
        LabGuru-compatible dict.
    """
    exporter = LabGuruExporter()
    return exporter.export_to_labguru(result)
