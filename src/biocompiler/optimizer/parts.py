"""
BioCompiler Part Library — Standard Biological Parts Registry

Provides a registry of standard biological parts (promoters, CDS, terminators,
RBS, linkers) for gene construction.  Parts can be loaded from YAML/JSON files
or accessed from a built-in default library.

The built-in parts include commonly used sequences from synthetic biology:
  - T7 promoter and lac promoter (for E. coli expression)
  - T7 terminator
  - Common RBS sequences (B0034, B0032, etc.)
  - Standard linkers and scar sequences
"""

from __future__ import annotations

import copy
import json
import logging
import os
from dataclasses import dataclass, field
from typing import Sequence

logger = logging.getLogger(__name__)

__all__ = [
    "Part",
    "PartLibrary",
    "DEFAULT_PARTS",
]

# ---------------------------------------------------------------------------
# Part dataclass
# ---------------------------------------------------------------------------

@dataclass
class Part:
    """A standard biological part.

    Attributes:
        name: Unique part identifier (e.g., "T7_promoter", "B0034").
        part_type: Type of part — one of "promoter", "cds", "terminator",
            "rbs", "linker".
        sequence: DNA sequence of the part (5'→3').
        description: Human-readable description of the part.
        metadata: Additional metadata (organism, source, etc.).
    """

    name: str
    part_type: str  # "promoter", "cds", "terminator", "rbs", "linker"
    sequence: str
    description: str = ""
    metadata: dict = field(default_factory=dict)

    VALID_TYPES: frozenset[str] = frozenset(
        {"promoter", "cds", "terminator", "rbs", "linker"}
    )

    def __post_init__(self) -> None:
        if self.part_type not in self.VALID_TYPES:
            raise ValueError(
                f"Invalid part_type {self.part_type!r}. "
                f"Must be one of {sorted(self.VALID_TYPES)}"
            )
        if not self.name:
            raise ValueError("Part name must not be empty")
        if not self.sequence:
            raise ValueError(f"Part {self.name!r} sequence must not be empty")
        # Normalize sequence to uppercase
        self.sequence = self.sequence.upper()


# ---------------------------------------------------------------------------
# Built-in default parts
# ---------------------------------------------------------------------------

DEFAULT_PARTS: list[Part] = [
    Part(
        name="T7_promoter",
        part_type="promoter",
        sequence="TAATACGACTCACTATAGGG",
        description="T7 RNA polymerase promoter. Strong promoter for in vitro transcription and E. coli expression systems.",
        metadata={"organism": "E_coli", "source": "T7 bacteriophage", "strength": "strong"},
    ),
    Part(
        name="lac_promoter",
        part_type="promoter",
        sequence="TTGACAATTATGCTTTCAGCTATTTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATT",
        description="Lac operon promoter. Inducible promoter controlled by LacI repressor. IPTG-inducible.",
        metadata={"organism": "E_coli", "source": "E. coli lac operon", "strength": "medium", "inducible": True},
    ),
    Part(
        name="T7_terminator",
        part_type="terminator",
        sequence="CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG",
        description="T7 RNA polymerase terminator. Forms a stable hairpin for transcription termination.",
        metadata={"organism": "E_coli", "source": "T7 bacteriophage"},
    ),
    Part(
        name="B0034",
        part_type="rbs",
        # Canonical iGEM registry sequence for BBa_B0034 (Anderson RBS
        # library, strong).  The previous value ``AAAGAGGAGATATACAT`` was a
        # fabricated hybrid that appended a spacer+ATG region containing a
        # wrong base — see audit issue H13.
        sequence="AAAGAGGAGAAA",
        description="Strong RBS (Ribosome Binding Site) from the Anderson library. Widely used in E. coli synthetic biology.",
        metadata={"organism": "E_coli", "source": "Anderson RBS library", "strength": "strong"},
    ),
    Part(
        name="B0032",
        part_type="rbs",
        # Canonical iGEM registry sequence for BBa_B0032 (Anderson RBS
        # library, medium).  Previous value was fabricated — see H13.
        sequence="TTCACACAGGAAACAGCT",
        description="Medium-strength RBS from the Anderson library.",
        metadata={"organism": "E_coli", "source": "Anderson RBS library", "strength": "medium"},
    ),
    Part(
        name="B0031",
        part_type="rbs",
        # Canonical iGEM registry sequence for BBa_B0031 (Anderson RBS
        # library, weak).  Previous value was fabricated — see H13.
        sequence="TTCACACAGGAAACAGC",
        description="Weak RBS from the Anderson library.",
        metadata={"organism": "E_coli", "source": "Anderson RBS library", "strength": "weak"},
    ),
    Part(
        name="CMV_promoter",
        part_type="promoter",
        sequence="CGCAAATGGGCGGTAGGCGTGACGGTGGGAGGTGCATGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAA",
        description="CMV immediate-early promoter. Strong constitutive promoter for mammalian expression.",
        metadata={"organism": "Homo_sapiens", "source": "Human cytomegalovirus", "strength": "strong"},
    ),
    Part(
        name="SV40_terminator",
        part_type="terminator",
        sequence="ATCTCTAGAGGCCCGGGGATCCACCGGTACTAGTCCATTTTTTACCACATACCCACCATATCCACATACCCACATACCCATTATACCCACATACCCATATACCCACATACCCACATATACCA",
        description="SV40 polyadenylation signal. Terminator for mammalian expression constructs.",
        metadata={"organism": "Homo_sapiens", "source": "Simian Virus 40"},
    ),
    Part(
        name="GoldenGate_scar",
        part_type="linker",
        sequence="TACTAG",
        description="Standard Golden Gate Assembly scar sequence (MoClo/BsmbI-compatible).",
        metadata={"assembly": "golden_gate", "standard": "MoClo"},
    ),
    Part(
        name="Gibson_overlap_20",
        part_type="linker",
        sequence="GCTAGCTAGCTAGCTAGCTA",
        description="Generic 20bp overlap sequence for Gibson Assembly.",
        metadata={"assembly": "gibson", "overlap_length": 20},
    ),
]


# ---------------------------------------------------------------------------
# PartLibrary class
# ---------------------------------------------------------------------------

class PartLibrary:
    """A registry of biological parts with lookup and search capabilities.

    Parts are stored by name and can be retrieved, searched by type and
    organism, or added dynamically.  A library can be loaded from a YAML
    or JSON file, or initialized with the built-in default parts.

    Args:
        library_path: Path to a YAML or JSON file defining parts.
            If ``None``, the library is populated with built-in default parts.
    """

    def __init__(self, library_path: str | None = None) -> None:
        self._parts: dict[str, Part] = {}

        # Load default parts first.  ``copy.deepcopy`` is essential here:
        # ``DEFAULT_PARTS`` is a module-level list of mutable ``Part``
        # dataclass instances (each carrying a mutable ``metadata`` dict).
        # Without the deep copy, every ``PartLibrary`` instance would share
        # the same ``Part`` object references, so a caller that mutated
        # ``lib.get("B0034").metadata[...]`` would silently corrupt the
        # global default for all other libraries — see audit issue H14.
        for part in DEFAULT_PARTS:
            self._parts[part.name] = copy.deepcopy(part)

        # Overlay with file-based parts if a path is given
        if library_path is not None:
            self._load_from_file(library_path)

    def _load_from_file(self, path: str) -> None:
        """Load parts from a YAML or JSON file.

        The file should contain a list of part dictionaries with keys:
        ``name``, ``part_type``, ``sequence``, ``description`` (optional),
        ``metadata`` (optional).

        Args:
            path: Path to the parts file.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file format is invalid or a part is malformed.
        """
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Parts library file not found: {path}")

        with open(path, "r") as f:
            raw = f.read()

        if path.endswith((".yaml", ".yml")):
            try:
                import yaml
                data = yaml.safe_load(raw)
            except ImportError:
                raise ImportError(
                    "PyYAML is required to load YAML part libraries. "
                    "Install with: pip install pyyaml"
                )
        elif path.endswith(".json"):
            data = json.loads(raw)
        else:
            raise ValueError(
                f"Unsupported file format: {path!r}. "
                "Use .yaml, .yml, or .json"
            )

        if not isinstance(data, list):
            raise ValueError(
                f"Parts library must be a list of part dicts, got {type(data).__name__}"
            )

        for entry in data:
            if not isinstance(entry, dict):
                raise ValueError(
                    f"Each part entry must be a dict, got {type(entry).__name__}"
                )
            name = entry.get("name")
            part_type = entry.get("part_type")
            sequence = entry.get("sequence", "")
            description = entry.get("description", "")
            metadata = entry.get("metadata", {})

            if not name or not part_type or not sequence:
                raise ValueError(
                    f"Part entry must have 'name', 'part_type', and 'sequence': {entry}"
                )

            part = Part(
                name=name,
                part_type=part_type,
                sequence=sequence,
                description=description,
                metadata=metadata,
            )
            self._parts[part.name] = part

    def get(self, name: str) -> Part:
        """Retrieve a part by name.

        Args:
            name: Unique part identifier.

        Returns:
            The :class:`Part` with the given name.

        Raises:
            KeyError: If no part with the given name exists.
        """
        if name not in self._parts:
            raise KeyError(
                f"Part {name!r} not found in library. "
                f"Available: {sorted(self._parts.keys())}"
            )
        return self._parts[name]

    def search(
        self,
        part_type: str,
        organism: str = "",
    ) -> list[Part]:
        """Search for parts by type and optionally by organism.

        Args:
            part_type: Type of part to search for (e.g., "promoter", "rbs").
            organism: If provided, only return parts whose metadata contains
                ``"organism"`` matching this value (case-insensitive).

        Returns:
            List of matching :class:`Part` objects.
        """
        results: list[Part] = []
        for part in self._parts.values():
            if part.part_type != part_type:
                continue
            if organism:
                part_org = part.metadata.get("organism", "").lower()
                if part_org != organism.lower():
                    continue
            results.append(part)
        return results

    def add(self, part: Part) -> None:
        """Add a part to the library.

        If a part with the same name already exists, it is replaced.

        Args:
            part: The :class:`Part` to add.
        """
        self._parts[part.name] = part

    def list_parts(self) -> list[str]:
        """Return a sorted list of all part names in the library."""
        return sorted(self._parts.keys())

    def __len__(self) -> int:
        return len(self._parts)

    def __contains__(self, name: str) -> bool:
        return name in self._parts

    def __repr__(self) -> str:
        return f"PartLibrary({len(self._parts)} parts)"
