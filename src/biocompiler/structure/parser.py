"""
BioCompiler Protein Structure Models and PDB I/O
==================================================
Data models for protein structure representation, PDB file parsing,
and structural analysis utilities (dihedral angles, Ramachandran,
secondary structure estimation).

This is a standalone data model module — no imports from other
biocompiler modules.
"""

from __future__ import annotations

import math
import logging
from dataclasses import dataclass, field
from typing import Any, TypedDict

logger = logging.getLogger(__name__)


# ────────────────────────────────────────────────────────────
# Amino acid code mappings
# ────────────────────────────────────────────────────────────
# NOTE: These are defined locally because constants.py does not yet export
# THREE_TO_ONE / ONE_TO_THREE. If they are added to constants.py in the
# future, import them from there instead of redefining here.

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    # Common modifications / non-standard
    "ASX": "B",   # Aspartic acid or Asparagine
    "GLX": "Z",   # Glutamic acid or Glutamine
    "SEC": "U",   # Selenocysteine
    "PYL": "O",   # Pyrrolysine
    "MSE": "M",   # Selenomethionine (treated as Met)
    "HSD": "H",   # Histidine delta-protonated
    "HSE": "H",   # Histidine epsilon-protonated
    "HSP": "H",   # Histidine doubly-protonated
    "CYX": "C",   # Cysteine in disulfide bond
    "HIP": "H",   # Protonated histidine
    "HIE": "H",   # Histidine epsilon tautomer
    "HID": "H",   # Histidine delta tautomer
    "LYP": "K",   # Lysine (protonated)
    "CYM": "C",   # Deprotonated cysteine
    "MLY": "K",   # N-methyl lysine
    "ACE": "X",   # Acetylated N-terminus
    "NME": "X",   # N-methyl amide C-terminus
}

ONE_TO_THREE: dict[str, str] = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
    "B": "ASX",
    "Z": "GLX",
    "U": "SEC",
    "O": "PYL",
    "X": "UNK",
}


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class Atom:
    """A single atom in a protein structure.

    Attributes:
        serial: Atom serial number.
        name: Atom name (N, CA, C, O, CB, etc.).
        residue_name: 3-letter amino acid code.
        chain_id: Chain identifier.
        residue_seq: Residue sequence number.
        x: X coordinate in Angstroms.
        y: Y coordinate in Angstroms.
        z: Z coordinate in Angstroms.
        occupancy: Occupancy factor (0.0–1.0).
        temp_factor: Temperature factor / B-factor (pLDDT in predicted models).
        element: Element symbol.
        insertion_code: Insertion code for residue numbering.
    """

    serial: int
    name: str
    residue_name: str
    chain_id: str
    residue_seq: int
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    temp_factor: float = 0.0
    element: str = ""
    insertion_code: str = ""

    def distance_to(self, other: Atom) -> float:
        """Compute Euclidean distance to another atom.

        Args:
            other: Another Atom instance.

        Returns:
            Distance in Angstroms.
        """
        dx = self.x - other.x
        dy = self.y - other.y
        dz = self.z - other.z
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    def to_pdb_line(self) -> str:
        """Format this atom as a PDB ATOM record (fixed-width columns).

        PDB ATOM record format (columns are 1-indexed):
            1-6   Record name   "ATOM  "
            7-11  Serial        integer, right-justified
            12    blank
            13-16 Atom name     left-justified if <= 3 chars, else right-justified
            17    Alt location  blank or A/B/etc.
            18-20 Residue name  right-justified
            21    blank
            22    Chain ID
            23-26 Residue seq   integer, right-justified
            27    Insertion code
            28-30 blank
            31-38 X             8.3 format
            39-46 Y             8.3 format
            47-54 Z             8.3 format
            55-60 Occupancy     6.2 format
            61-66 Temp factor   6.2 format
            67-76 blank
            77-78 Element       right-justified

        Returns:
            A string of exactly 80 characters (PDB spec line length).
        """
        # Atom name: 4-character field (columns 13-16)
        # If atom name is <= 3 chars, it starts at column 14 (space-padded left)
        # If atom name is 4 chars, it starts at column 13
        if len(self.name) < 4:
            atom_name = f" {self.name:<3s}"
        else:
            atom_name = f"{self.name:<4s}"

        # Residue name: 3-char right-justified in columns 18-20
        res_name = f"{self.residue_name:>3s}"

        # Element: right-justified in columns 77-78
        elem = f"{self.element:>2s}" if self.element else "  "

        line = (
            f"ATOM  "
            f"{self.serial:>5d}"
            f" "
            f"{atom_name}"
            f" "           # alt location
            f"{res_name}"
            f" "           # blank column 21
            f"{self.chain_id}"
            f"{self.residue_seq:>4d}"
            f"{self.insertion_code}"
            f"   "         # columns 28-30 blank
            f"{self.x:>8.3f}"
            f"{self.y:>8.3f}"
            f"{self.z:>8.3f}"
            f"{self.occupancy:>6.2f}"
            f"{self.temp_factor:>6.2f}"
            f"          "  # columns 67-76 blank
            f"{elem}"
        )

        # Pad or trim to exactly 80 characters
        return line.ljust(80)


@dataclass
class Residue:
    """An amino acid residue in a protein structure.

    Attributes:
        seq_num: Residue sequence number.
        name: 3-letter amino acid code.
        chain_id: Chain identifier.
        atoms: Dictionary mapping atom names to Atom objects.
        insertion_code: Insertion code for this residue.
    """

    seq_num: int
    name: str
    chain_id: str
    atoms: dict[str, Atom] = field(default_factory=dict)
    insertion_code: str = ""

    def ca(self) -> Atom | None:
        """Get the CA (alpha carbon) atom of this residue.

        Returns:
            The CA Atom if present, else None.
        """
        return self.atoms.get("CA")

    def one_letter(self) -> str:
        """Convert 3-letter residue code to 1-letter code.

        Returns:
            1-letter amino acid code, or 'X' if unknown.
        """
        return THREE_TO_ONE.get(self.name, "X")

    def centroid(self) -> tuple[float, float, float]:
        """Compute the geometric center of all atoms in this residue.

        Returns:
            Tuple of (x, y, z) coordinates. Returns (0, 0, 0) if
            no atoms are present.
        """
        if not self.atoms:
            return (0.0, 0.0, 0.0)
        n = len(self.atoms)
        cx = sum(a.x for a in self.atoms.values()) / n
        cy = sum(a.y for a in self.atoms.values()) / n
        cz = sum(a.z for a in self.atoms.values()) / n
        return (cx, cy, cz)


@dataclass
class Chain:
    """A protein chain (one continuous polypeptide).

    Attributes:
        chain_id: Chain identifier (single character).
        residues: Ordered list of residues in this chain.
    """

    chain_id: str
    residues: list[Residue] = field(default_factory=list)

    def sequence(self) -> str:
        """Extract the 1-letter protein sequence from this chain.

        Returns:
            String of 1-letter amino acid codes.
        """
        return "".join(r.one_letter() for r in self.residues)

    def ca_coords(self) -> list[tuple[float, float, float]]:
        """Get CA coordinates for all residues in this chain.

        Returns:
            List of (x, y, z) tuples. Residues without a CA atom
            are skipped.
        """
        coords: list[tuple[float, float, float]] = []
        for r in self.residues:
            ca_atom = r.ca()
            if ca_atom is not None:
                coords.append((ca_atom.x, ca_atom.y, ca_atom.z))
        return coords


@dataclass
class ProteinStructure:
    """A complete protein structure with one or more chains.

    Attributes:
        chains: List of Chain objects.
        atoms: Flat list of all atoms across all chains.
        pdb_string: Original PDB text.
        source: Origin of the structure (e.g., "esmfold", "pdb_file").
    """

    chains: list[Chain] = field(default_factory=list)
    atoms: list[Atom] = field(default_factory=list)
    pdb_string: str = ""
    source: str = ""

    def sequence(self) -> str:
        """Extract the full protein sequence across all chains.

        Returns:
            Concatenated 1-letter amino acid sequence.
        """
        return "".join(c.sequence() for c in self.chains)

    def get_chain(self, chain_id: str) -> Chain | None:
        """Get a chain by its identifier.

        Args:
            chain_id: Single-character chain identifier.

        Returns:
            The Chain object if found, else None.
        """
        for c in self.chains:
            if c.chain_id == chain_id:
                return c
        return None

    def ca_distance_matrix(self) -> list[list[float]]:
        """Compute pairwise CA distance matrix across all residues.

        Returns:
            Square matrix of distances. Missing CA atoms use
            float('inf') as a placeholder.
        """
        ca_list: list[tuple[float, float, float]] = []
        for chain in self.chains:
            for res in chain.residues:
                ca_atom = res.ca()
                if ca_atom is not None:
                    ca_list.append((ca_atom.x, ca_atom.y, ca_atom.z))
                else:
                    ca_list.append((float("inf"), float("inf"), float("inf")))

        n = len(ca_list)
        matrix: list[list[float]] = [[0.0] * n for _ in range(n)]
        for i in range(n):
            xi, yi, zi = ca_list[i]
            for j in range(i + 1, n):
                xj, yj, zj = ca_list[j]
                d = math.sqrt(
                    (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
                )
                matrix[i][j] = d
                matrix[j][i] = d
        return matrix

    def residue_count(self) -> int:
        """Count total number of residues across all chains.

        Returns:
            Total residue count.
        """
        return sum(len(c.residues) for c in self.chains)

    def plddt_scores(self) -> list[float]:
        """Extract pLDDT scores from B-factors (CA atoms only).

        In AlphaFold/ESMFold models, the B-factor column contains
        pLDDT confidence scores (0–100).

        Returns:
            List of pLDDT scores, one per residue.
        """
        scores: list[float] = []
        for chain in self.chains:
            for res in chain.residues:
                ca_atom = res.ca()
                if ca_atom is not None:
                    scores.append(ca_atom.temp_factor)
                else:
                    scores.append(0.0)
        return scores

    def mean_plddt(self) -> float:
        """Compute mean pLDDT score across all residues.

        Returns:
            Mean pLDDT score. Returns 0.0 if no residues.
        """
        scores = self.plddt_scores()
        if not scores:
            return 0.0
        return sum(scores) / len(scores)

    def contact_map(self, threshold: float = 8.0) -> list[list[bool]]:
        """Compute a binary contact map from CA distances.

        Two residues are in contact if their CA atoms are within
        the threshold distance.

        Args:
            threshold: Distance cutoff in Angstroms (default 8.0).

        Returns:
            Square boolean matrix where True indicates a contact.
        """
        dist = self.ca_distance_matrix()
        n = len(dist)
        cmap: list[list[bool]] = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if dist[i][j] <= threshold:
                    cmap[i][j] = True
                    cmap[j][i] = True
        return cmap

    def radius_of_gyration(self) -> float:
        """Compute radius of gyration (Rg) from CA coordinates.

        Rg = sqrt(mean(sum_i(|r_i - r_mean|^2)))

        Returns:
            Radius of gyration in Angstroms. Returns 0.0 if fewer
            than 2 CA atoms.
        """
        coords: list[tuple[float, float, float]] = []
        for chain in self.chains:
            coords.extend(chain.ca_coords())

        if len(coords) < 2:
            return 0.0

        n = len(coords)
        # Compute mean position
        mx = sum(c[0] for c in coords) / n
        my = sum(c[1] for c in coords) / n
        mz = sum(c[2] for c in coords) / n

        # Compute sum of squared distances from mean
        ssd = sum(
            (c[0] - mx) ** 2 + (c[1] - my) ** 2 + (c[2] - mz) ** 2
            for c in coords
        )
        return math.sqrt(ssd / n)

    def to_pdb(self) -> str:
        """Regenerate PDB string from the structure data.

        Returns:
            PDB-formatted string with ATOM, TER, and END records.
        """
        lines: list[str] = []
        for chain in self.chains:
            for res in chain.residues:
                for atom in sorted(res.atoms.values(), key=lambda a: a.serial):
                    lines.append(atom.to_pdb_line())
            # TER record after each chain
            last_serial = 0
            if chain.residues:
                last_atom = max(
                    (a for r in chain.residues for a in r.atoms.values()),
                    key=lambda a: a.serial,
                )
                last_serial = last_atom.serial
            lines.append(f"TER   {last_serial + 1:>5d}      {chain.residues[-1].name:>3s} {chain.chain_id}{chain.residues[-1].seq_num:>4d}" if chain.residues else "TER")
        lines.append("END")
        return "\n".join(lines)


# ────────────────────────────────────────────────────────────
# PDB I/O functions
# ────────────────────────────────────────────────────────────

def parse_pdb(pdb_string: str, include_het: bool = False) -> ProteinStructure:
    """Parse a standard PDB format string into a ProteinStructure.

    Handles ATOM/HETATM records, TER records, and END records.
    For multi-model PDB files, only the first model is parsed.
    HETATM records (water, ligands) are skipped by default.

    PDB ATOM record column layout (1-indexed):
        1-6   Record name
        7-11  Atom serial number
        13-16 Atom name
        17    Alternate location indicator
        18-20 Residue name
        22    Chain ID
        23-26 Residue sequence number
        27    Code for insertion of residues
        31-38 X coordinate
        39-46 Y coordinate
        47-54 Z coordinate
        55-60 Occupancy
        61-66 Temperature factor
        77-78 Element symbol

    Args:
        pdb_string: Contents of a PDB file as a string.
        include_het: If True, include HETATM records. Default False.

    Returns:
        ProteinStructure with parsed atoms, residues, and chains.

    Raises:
        ValueError: If the PDB string is empty or contains no valid
            ATOM records.
    """
    if not pdb_string or not pdb_string.strip():
        raise ValueError("Empty PDB string provided")

    atoms: list[Atom] = []
    chain_map: dict[str, dict[tuple[int, str], list[Atom]]] = {}
    # chain_map: chain_id -> (residue_seq, insertion_code) -> list[Atom]

    in_first_model = True
    found_model = False

    for line in pdb_string.splitlines():
        line = line.rstrip()

        # Handle MODEL/ENDMDL for multi-model PDB
        if line.startswith("MODEL"):
            if found_model:
                # We've already seen a MODEL; skip subsequent models
                in_first_model = False
            found_model = True
            in_first_model = True
            continue

        if line.startswith("ENDMDL"):
            # Stop parsing after first model
            break

        if not in_first_model:
            continue

        record_type = line[0:6].strip() if len(line) >= 6 else ""

        if record_type == "ATOM" or (include_het and record_type == "HETATM"):
            try:
                atom = _parse_atom_line(line)
            except (ValueError, IndexError) as exc:
                logger.debug("Skipping malformed ATOM line: %s (%s)", line[:50], exc)
                continue

            atoms.append(atom)

            # Group by chain and residue
            if atom.chain_id not in chain_map:
                chain_map[atom.chain_id] = {}
            res_key = (atom.residue_seq, atom.insertion_code)
            if res_key not in chain_map[atom.chain_id]:
                chain_map[atom.chain_id][res_key] = []
            chain_map[atom.chain_id][res_key].append(atom)

        # TER and END records are structural; no data to extract
        elif record_type == "TER":
            continue
        elif record_type == "END":
            break

    if not atoms:
        raise ValueError("No valid ATOM records found in PDB string")

    # Build Chain -> Residue hierarchy
    chains: list[Chain] = []
    for chain_id in sorted(chain_map.keys()):
        residue_map = chain_map[chain_id]
        residues: list[Residue] = []
        for res_key in sorted(residue_map.keys()):
            res_seq, ins_code = res_key
            res_atoms = residue_map[res_key]
            # Determine residue name from first atom
            residue_name = res_atoms[0].residue_name
            atom_dict: dict[str, Atom] = {}
            for a in res_atoms:
                # Use atom name as key; for duplicate names (alt locations),
                # keep the one with higher occupancy
                if a.name not in atom_dict or a.occupancy > atom_dict[a.name].occupancy:
                    atom_dict[a.name] = a
            residues.append(Residue(
                seq_num=res_seq,
                name=residue_name,
                chain_id=chain_id,
                atoms=atom_dict,
                insertion_code=ins_code,
            ))
        chains.append(Chain(chain_id=chain_id, residues=residues))

    return ProteinStructure(
        chains=chains,
        atoms=atoms,
        pdb_string=pdb_string,
        source="pdb_string",
    )


def _parse_atom_line(line: str) -> Atom:
    """Parse a single PDB ATOM or HETATM line.

    Args:
        line: A PDB line (at least 54 characters for coordinates).

    Returns:
        Atom instance.

    Raises:
        ValueError: If the line is too short or has invalid numeric fields.
    """
    if len(line) < 54:
        raise ValueError(f"ATOM line too short ({len(line)} chars): {line}")

    # Record type (columns 1-6, 0-indexed: 0:6)
    # Serial number (columns 7-11, 0-indexed: 6:11)
    serial = int(line[6:11])

    # Atom name (columns 13-16, 0-indexed: 12:16)
    atom_name = line[12:16].strip()

    # Alternate location (column 17, 0-indexed: 16)
    # We skip alt-location handling; prefer the first or highest occupancy

    # Residue name (columns 18-20, 0-indexed: 17:20)
    residue_name = line[17:20].strip()

    # Chain ID (column 22, 0-indexed: 21)
    chain_id = line[21] if len(line) > 21 else " "

    # Residue sequence number (columns 23-26, 0-indexed: 22:26)
    residue_seq = int(line[22:26])

    # Insertion code (column 27, 0-indexed: 26)
    insertion_code = line[26].strip() if len(line) > 26 else ""

    # Coordinates (columns 31-54, 0-indexed: 30:54)
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])

    # Occupancy (columns 55-60, 0-indexed: 54:60)
    occupancy = float(line[54:60]) if len(line) >= 60 else 1.0

    # Temperature factor (columns 61-66, 0-indexed: 60:66)
    temp_factor = float(line[60:66]) if len(line) >= 66 else 0.0

    # Element symbol (columns 77-78, 0-indexed: 76:78)
    element = line[76:78].strip() if len(line) >= 78 else ""

    return Atom(
        serial=serial,
        name=atom_name,
        residue_name=residue_name,
        chain_id=chain_id,
        residue_seq=residue_seq,
        x=x,
        y=y,
        z=z,
        occupancy=occupancy,
        temp_factor=temp_factor,
        element=element,
        insertion_code=insertion_code,
    )


def parse_pdb_file(filepath: str, include_het: bool = False) -> ProteinStructure:
    """Read a PDB file and parse it into a ProteinStructure.

    Args:
        filepath: Path to the PDB file.
        include_het: If True, include HETATM records.

    Returns:
        ProteinStructure with parsed atoms, residues, and chains.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is empty or contains no valid ATOM records.
    """
    with open(filepath, "r") as fh:
        pdb_string = fh.read()

    structure = parse_pdb(pdb_string, include_het=include_het)
    structure.source = "pdb_file"
    return structure


# ────────────────────────────────────────────────────────────
# Structural analysis constants
# ────────────────────────────────────────────────────────────

# Threshold below which the central bond vector is considered degenerate
# (atoms are collinear), making the dihedral angle undefined.
DEGENERATE_DIHEDRAL_THRESHOLD: float = 1e-10

# Ramachandran angle targets for secondary structure estimation
# Source: Morris et al. (1992), simplified DSSP-like classification
RAMA_ALPHA_HELIX_PHI: float = -57.0
RAMA_ALPHA_HELIX_PSI: float = -47.0
RAMA_BETA_SHEET_PHI: float = -120.0
RAMA_BETA_SHEET_PSI: float = 120.0
RAMA_ANGLE_TOLERANCE: float = 30.0


class RamachandranResult(TypedDict):
    """Typed result of Ramachandran angle computation."""
    phi: list[float | None]
    psi: list[float | None]
    residues: list[str]


# ────────────────────────────────────────────────────────────
# Structural analysis functions
# ────────────────────────────────────────────────────────────

def compute_dihedral(
    p1: tuple[float, float, float],
    p2: tuple[float, float, float],
    p3: tuple[float, float, float],
    p4: tuple[float, float, float],
) -> float:
    """Compute the dihedral angle defined by four 3D points.

    The dihedral angle is the angle between the planes defined by
    (p1, p2, p3) and (p2, p3, p4), measured around the p2-p3 bond.

    Uses the formula based on cross products:
        b1 = p2 - p1
        b2 = p3 - p2
        b3 = p4 - p3
        n1 = b1 x b2
        n2 = b2 x b3
        m1 = n1 x (b2 / |b2|)
        x = dot(n1, n2)
        y = dot(m1, n2)
        angle = -atan2(y, x)

    Args:
        p1: First point coordinates.
        p2: Second point coordinates.
        p3: Third point coordinates.
        p4: Fourth point coordinates.

    Returns:
        Dihedral angle in degrees, in the range [-180, 180].
    """
    # Vectors
    b1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    b2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])
    b3 = (p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2])

    # Normal vectors to planes
    n1 = _cross(b1, b2)
    n2 = _cross(b2, b3)

    # Normalized b2 for m1 computation
    b2_len = math.sqrt(b2[0] ** 2 + b2[1] ** 2 + b2[2] ** 2)
    if b2_len < DEGENERATE_DIHEDRAL_THRESHOLD:
        return 0.0
    b2_norm = (b2[0] / b2_len, b2[1] / b2_len, b2[2] / b2_len)

    # m1 = n1 x b2_norm
    m1 = _cross(n1, b2_norm)

    # x = dot(n1, n2), y = dot(m1, n2)
    x = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]
    y = m1[0] * n2[0] + m1[1] * n2[1] + m1[2] * n2[2]

    # Dihedral angle
    angle = -math.degrees(math.atan2(y, x))

    # Normalize to [-180, 180]
    while angle <= -180.0:
        angle += 360.0
    while angle > 180.0:
        angle -= 360.0

    return angle


def _cross(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
) -> tuple[float, float, float]:
    """Compute the cross product of two 3D vectors."""
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def compute_ramachandran(structure: ProteinStructure) -> RamachandranResult:
    """Compute phi/psi dihedral angles for each residue (except termini).

    Phi (φ) is defined by atoms: C(i-1) - N(i) - CA(i) - C(i)
    Psi (ψ) is defined by atoms: N(i) - CA(i) - C(i) - N(i+1)

    The first residue has no phi angle, the last has no psi angle.

    Args:
        structure: ProteinStructure to analyze.

    Returns:
        Dictionary with keys:
            "phi": list of float — phi angles in degrees (None for first residue)
            "psi": list of float — psi angles in degrees (None for last residue)
            "residues": list of str — 1-letter residue codes
    """
    phi: list[float | None] = []
    psi: list[float | None] = []
    residues: list[str] = []

    # Collect all residues across chains sequentially
    all_residues: list[Residue] = []
    for chain in structure.chains:
        all_residues.extend(chain.residues)

    n = len(all_residues)
    for i in range(n):
        residues.append(all_residues[i].one_letter())

        # Phi: C(i-1) - N(i) - CA(i) - C(i)
        if i > 0:
            c_prev = all_residues[i - 1].atoms.get("C")
            n_curr = all_residues[i].atoms.get("N")
            ca_curr = all_residues[i].atoms.get("CA")
            c_curr = all_residues[i].atoms.get("C")
            if c_prev and n_curr and ca_curr and c_curr:
                phi_angle = compute_dihedral(
                    (c_prev.x, c_prev.y, c_prev.z),
                    (n_curr.x, n_curr.y, n_curr.z),
                    (ca_curr.x, ca_curr.y, ca_curr.z),
                    (c_curr.x, c_curr.y, c_curr.z),
                )
                phi.append(phi_angle)
            else:
                phi.append(None)
        else:
            phi.append(None)

        # Psi: N(i) - CA(i) - C(i) - N(i+1)
        if i < n - 1:
            n_curr = all_residues[i].atoms.get("N")
            ca_curr = all_residues[i].atoms.get("CA")
            c_curr = all_residues[i].atoms.get("C")
            n_next = all_residues[i + 1].atoms.get("N")
            if n_curr and ca_curr and c_curr and n_next:
                psi_angle = compute_dihedral(
                    (n_curr.x, n_curr.y, n_curr.z),
                    (ca_curr.x, ca_curr.y, ca_curr.z),
                    (c_curr.x, c_curr.y, c_curr.z),
                    (n_next.x, n_next.y, n_next.z),
                )
                psi.append(psi_angle)
            else:
                psi.append(None)
        else:
            psi.append(None)

    return {
        "phi": phi,
        "psi": psi,
        "residues": residues,
    }


def secondary_structure_estimate(structure: ProteinStructure) -> list[str]:
    """Estimate secondary structure from phi/psi dihedral angles.

    Uses a simplified DSSP-like classification:
        - Alpha-helix (H): phi ≈ -57°, psi ≈ -47° (within 30°)
        - Beta-sheet (E): phi ≈ -120°, psi ≈ 120° (within 30°)
        - Coil (C): everything else

    First and last residues are always classified as coil (no
    complete phi/psi pair).

    Args:
        structure: ProteinStructure to analyze.

    Returns:
        List of secondary structure codes ('H', 'E', 'C'), one per
        residue.
    """
    rama = compute_ramachandran(structure)
    phi_list = rama["phi"]
    psi_list = rama["psi"]
    n = len(phi_list)

    ss: list[str] = []
    for i in range(n):
        p = phi_list[i]
        s = psi_list[i]

        if p is None or s is None:
            ss.append("C")
            continue

        # Alpha-helix: phi ≈ RAMA_ALPHA_HELIX_PHI, psi ≈ RAMA_ALPHA_HELIX_PSI
        if abs(p - RAMA_ALPHA_HELIX_PHI) <= RAMA_ANGLE_TOLERANCE and abs(s - RAMA_ALPHA_HELIX_PSI) <= RAMA_ANGLE_TOLERANCE:
            ss.append("H")
        # Beta-sheet: phi ≈ RAMA_BETA_SHEET_PHI, psi ≈ RAMA_BETA_SHEET_PSI
        elif abs(p - RAMA_BETA_SHEET_PHI) <= RAMA_ANGLE_TOLERANCE and abs(s - RAMA_BETA_SHEET_PSI) <= RAMA_ANGLE_TOLERANCE:
            ss.append("E")
        else:
            ss.append("C")

    return ss
