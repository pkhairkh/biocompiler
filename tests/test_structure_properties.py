"""
Property-Based Tests for BioCompiler Structure Modules
=======================================================

Verifies structural invariants and consistency properties of the
structure subpackage using Hypothesis-based property testing.

Three core properties tested:
  1. THREE_TO_ONE / ONE_TO_THREE are inverses for the 20 standard AAs
  2. Parsed structures have valid atom counts (every atom belongs to a
     residue, residue atom counts are non-negative, total atom count
     matches the flat list length)
  3. GC content computed from DNA sequences matches an independent
     manual count of G + C bases divided by sequence length
"""

import pytest
pytest.importorskip("hypothesis")
pytest.importorskip("hypothesis")
from hypothesis import given, settings, assume, HealthCheck
from hypothesis import strategies as st

from biocompiler.structure.parser import (
    THREE_TO_ONE,
    ONE_TO_THREE,
    parse_pdb,
    Atom,
    Residue,
    Chain,
    ProteinStructure,
)
from biocompiler.sequence.scanner import gc_content


# ────────────────────────────────────────────────────────────
# Shared Constants
# ────────────────────────────────────────────────────────────

# The 20 canonical one-letter amino acid codes
STANDARD_AA_ONE: str = "ACDEFGHIKLMNPQRSTVWY"

# The 20 canonical three-letter amino acid codes
STANDARD_AA_THREE: list[str] = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
]


# ────────────────────────────────────────────────────────────
# Hypothesis Strategies
# ────────────────────────────────────────────────────────────

# Single standard amino acid (one-letter)
standard_aa = st.sampled_from(list(STANDARD_AA_ONE))

# DNA sequences of moderate length
dna_medium = st.text(alphabet="ACGT", min_size=3, max_size=300)

# Short DNA for boundary testing
dna_short = st.text(alphabet="ACGT", min_size=0, max_size=10)

# Longer DNA for thorough GC distribution testing
dna_long = st.text(alphabet="ACGT", min_size=30, max_size=600)

# Integer for atom serial numbers
atom_serial = st.integers(min_value=1, max_value=99999)

# Integer for residue sequence numbers
residue_seq = st.integers(min_value=1, max_value=9999)

# Float for 3D coordinates (reasonable range for protein structures)
coord_float = st.floats(min_value=-500.0, max_value=500.0, allow_nan=False, allow_infinity=False)

# Float for occupancy
occupancy_float = st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False)

# Float for B-factor / temp_factor
temp_factor_float = st.floats(min_value=0.0, max_value=100.0, allow_nan=False, allow_infinity=False)


# ────────────────────────────────────────────────────────────
# Helper: Build a minimal PDB string from atoms
# ────────────────────────────────────────────────────────────

def _make_pdb_line(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_seq: int,
    x: float,
    y: float,
    z: float,
    occupancy: float = 1.0,
    temp_factor: float = 0.0,
    element: str = "",
) -> str:
    """Format a single PDB ATOM record line (80-char fixed-width).

    Follows the PDB format columns used by parse_pdb / _parse_atom_line.
    """
    # Atom name: 4-char field (columns 13-16)
    if len(atom_name) < 4:
        atom_name_fmt = f" {atom_name:<3s}"
    else:
        atom_name_fmt = f"{atom_name:<4s}"

    res_name = f"{residue_name:>3s}"
    elem = f"{element:>2s}" if element else "  "

    line = (
        f"ATOM  "
        f"{serial:>5d}"
        f" "
        f"{atom_name_fmt}"
        f" "          # alt location
        f"{res_name}"
        f" "          # blank col 21
        f"{chain_id}"
        f"{residue_seq:>4d}"
        f" "          # insertion code
        f"   "        # cols 28-30
        f"{x:>8.3f}"
        f"{y:>8.3f}"
        f"{z:>8.3f}"
        f"{occupancy:>6.2f}"
        f"{temp_factor:>6.2f}"
        f"          "  # cols 67-76
        f"{elem}"
    )
    return line.ljust(80)


def _build_pdb_string(
    atoms: list[dict],
    chain_id: str = "A",
) -> str:
    """Build a complete PDB string from a list of atom dicts.

    Each dict must have keys:
        serial, atom_name, residue_name, residue_seq,
        x, y, z, occupancy, temp_factor, element
    """
    lines = []
    for a in atoms:
        lines.append(_make_pdb_line(
            serial=a["serial"],
            atom_name=a["atom_name"],
            residue_name=a["residue_name"],
            chain_id=chain_id,
            residue_seq=a["residue_seq"],
            x=a["x"], y=a["y"], z=a["z"],
            occupancy=a.get("occupancy", 1.0),
            temp_factor=a.get("temp_factor", 0.0),
            element=a.get("element", ""),
        ))
    lines.append("END")
    return "\n".join(lines)


def _make_single_residue_pdb(
    residue_name: str = "ALA",
    residue_seq: int = 1,
    chain_id: str = "A",
    n_atoms: int = 5,
) -> str:
    """Build a PDB string for a single residue with n_atoms atoms.

    Atoms are named CA, N, C, O, CB, CG, ... up to n_atoms.
    """
    backbone = ["N", "CA", "C", "O"]
    sidechain = ["CB", "CG", "OD1", "OD2", "ND1", "ND2", "CE", "NZ",
                 "SD", "CE1", "CE2", "NE1", "NE2", "CZ", "CH2", "CD",
                 "CD1", "CD2", "CE3", "CZ2", "CZ3"]
    all_names = (backbone + sidechain)[:n_atoms]
    # If n_atoms > available names, pad with generic H atoms
    while len(all_names) < n_atoms:
        all_names.append(f"H{len(all_names)}")

    atoms = []
    for i, atom_name in enumerate(all_names[:n_atoms]):
        atoms.append({
            "serial": i + 1,
            "atom_name": atom_name,
            "residue_name": residue_name,
            "residue_seq": residue_seq,
            "x": float(i) * 1.5,
            "y": float(i) * 1.5,
            "z": float(i) * 1.5,
            "occupancy": 1.0,
            "temp_factor": 50.0,
            "element": atom_name[0],
        })
    return _build_pdb_string(atoms, chain_id=chain_id)


def _manual_gc_content(seq: str) -> float:
    """Independent GC content computation: (G+C)/len, rounded to 4 decimals."""
    if not seq:
        return 0.0
    seq_up = seq.upper()
    gc = seq_up.count("G") + seq_up.count("C")
    return round(gc / len(seq_up), 4)


# ══════════════════════════════════════════════════════════════
# TEST CLASS 1: THREE_TO_ONE / ONE_TO_THREE Inverse Property
# ══════════════════════════════════════════════════════════════

class TestAminoAcidMappingInverse:
    """Property: THREE_TO_ONE and ONE_TO_THREE are inverses for standard AAs.

    For every standard amino acid:
        THREE_TO_ONE[ONE_TO_THREE[one]] == one
        ONE_TO_THREE[THREE_TO_ONE[three]] == three

    Note: THREE_TO_ONE includes non-standard entries (e.g., MSE→M, HSD→H)
    that map to the same one-letter code as a standard AA. For those,
    the forward composition THREE_TO_ONE[ONE_TO_THREE[x]] still holds,
    but the reverse composition ONE_TO_THREE[THREE_TO_ONE[x]] picks
    the canonical three-letter code. We test both directions carefully.
    """

    def test_standard_three_to_one_covers_all_20(self):
        """All 20 standard three-letter codes are keys in THREE_TO_ONE."""
        for three in STANDARD_AA_THREE:
            assert three in THREE_TO_ONE, f"Standard AA {three} missing from THREE_TO_ONE"

    def test_standard_one_to_three_covers_all_20(self):
        """All 20 standard one-letter codes are keys in ONE_TO_THREE."""
        for one in STANDARD_AA_ONE:
            assert one in ONE_TO_THREE, f"Standard AA {one} missing from ONE_TO_THREE"

    def test_standard_three_to_one_values_are_standard(self):
        """THREE_TO_ONE maps each standard 3-letter code to a valid 1-letter code."""
        for three in STANDARD_AA_THREE:
            one = THREE_TO_ONE[three]
            assert one in STANDARD_AA_ONE, (
                f"THREE_TO_ONE[{three}] = {one!r}, not a standard one-letter code"
            )

    def test_standard_one_to_three_values_are_standard(self):
        """ONE_TO_THREE maps each standard 1-letter code to a valid 3-letter code."""
        for one in STANDARD_AA_ONE:
            three = ONE_TO_THREE[one]
            assert three in STANDARD_AA_THREE, (
                f"ONE_TO_THREE[{one}] = {three!r}, not a standard three-letter code"
            )

    @given(one=standard_aa)
    @settings(max_examples=50)
    def test_roundtrip_one_to_three_to_one(self, one):
        """THREE_TO_ONE[ONE_TO_THREE[x]] == x for standard one-letter codes."""
        three = ONE_TO_THREE[one]
        assert three in THREE_TO_ONE, (
            f"ONE_TO_THREE[{one}] = {three!r} not in THREE_TO_ONE"
        )
        assert THREE_TO_ONE[three] == one, (
            f"THREE_TO_ONE[ONE_TO_THREE[{one}]] = THREE_TO_ONE[{three}] = "
            f"{THREE_TO_ONE.get(three, '<missing>')!r}, expected {one!r}"
        )

    @given(three=st.sampled_from(STANDARD_AA_THREE))
    @settings(max_examples=50)
    def test_roundtrip_three_to_one_to_three(self, three):
        """ONE_TO_THREE[THREE_TO_ONE[x]] == x for standard three-letter codes."""
        one = THREE_TO_ONE[three]
        assert one in ONE_TO_THREE, (
            f"THREE_TO_ONE[{three}] = {one!r} not in ONE_TO_THREE"
        )
        result_three = ONE_TO_THREE[one]
        assert result_three == three, (
            f"ONE_TO_THREE[THREE_TO_ONE[{three}]] = ONE_TO_THREE[{one}] = "
            f"{result_three!r}, expected {three!r}"
        )

    def test_mapping_is_bijective_on_standard_aas(self):
        """The restriction to standard AAs is a bijection (no two 3-letter codes
        map to the same 1-letter code within the standard 20)."""
        one_letter_values = [THREE_TO_ONE[three] for three in STANDARD_AA_THREE]
        assert len(one_letter_values) == len(set(one_letter_values)), (
            f"THREE_TO_ONE is not injective on standard AAs: "
            f"values = {one_letter_values}"
        )

    def test_nonstandard_entries_map_to_standard_one_letter(self):
        """Non-standard three-letter codes in THREE_TO_ONE map to valid
        one-letter codes (possibly shared with a standard AA)."""
        nonstandard_keys = set(THREE_TO_ONE.keys()) - set(STANDARD_AA_THREE)
        for three in nonstandard_keys:
            one = THREE_TO_ONE[three]
            # Non-standard entries may map to B, Z, U, O, X, or standard AA codes
            valid_one_letter = set(STANDARD_AA_ONE) | {"B", "Z", "U", "O", "X"}
            assert one in valid_one_letter, (
                f"Non-standard THREE_TO_ONE[{three}] = {one!r} is not a "
                f"recognized one-letter code"
            )

    @given(one=standard_aa)
    @settings(max_examples=20)
    def test_one_to_three_inverse_is_left_identity(self, one):
        """THREE_TO_ONE is a left inverse of ONE_TO_THREE for standard AAs.

        i.e., THREE_TO_ONE ∘ ONE_TO_THREE = id on the set of standard
        one-letter codes.
        """
        assert THREE_TO_ONE[ONE_TO_THREE[one]] == one

    @given(three=st.sampled_from(STANDARD_AA_THREE))
    @settings(max_examples=20)
    def test_three_to_one_inverse_is_left_identity(self, three):
        """ONE_TO_THREE is a left inverse of THREE_TO_ONE for standard AAs.

        i.e., ONE_TO_THREE ∘ THREE_TO_ONE = id on the set of standard
        three-letter codes.
        """
        assert ONE_TO_THREE[THREE_TO_ONE[three]] == three


# ══════════════════════════════════════════════════════════════
# TEST CLASS 2: Parsed Structure Atom Count Properties
# ══════════════════════════════════════════════════════════════

class TestParsedStructureAtomCounts:
    """Property: Parsed structures have valid atom counts.

    After parsing a PDB string, the resulting ProteinStructure must
    satisfy:
      1. Total flat atom list length equals sum of all residue atom counts
      2. Each residue has at least 1 atom
      3. Each chain has at least 1 residue (if the structure is non-empty)
      4. The structure has at least 1 chain (if atoms were parsed)
    """

    def test_single_residue_alanine(self):
        """Parse a single ALA residue (5 backbone + CB = 5 atoms)."""
        pdb_str = _make_single_residue_pdb("ALA", residue_seq=1, n_atoms=5)
        struct = parse_pdb(pdb_str)
        assert struct.residue_count() == 1
        assert len(struct.chains) == 1
        chain = struct.chains[0]
        assert chain.chain_id == "A"
        assert len(chain.residues) == 1
        res = chain.residues[0]
        assert res.name == "ALA"
        assert len(res.atoms) == 5
        assert len(struct.atoms) == 5

    def test_single_residue_various_atom_counts(self):
        """Residues with different atom counts parse correctly."""
        for n_atoms in [1, 3, 5, 10, 20]:
            pdb_str = _make_single_residue_pdb("GLY", residue_seq=1, n_atoms=n_atoms)
            struct = parse_pdb(pdb_str)
            assert struct.residue_count() == 1
            res = struct.chains[0].residues[0]
            assert len(res.atoms) == n_atoms, (
                f"Expected {n_atoms} atoms, got {len(res.atoms)}"
            )
            assert len(struct.atoms) == n_atoms

    def test_two_residues_same_chain(self):
        """Parse two residues in the same chain."""
        atoms = []
        for res_seq in [1, 2]:
            for i, atom_name in enumerate(["N", "CA", "C", "O"]):
                atoms.append({
                    "serial": (res_seq - 1) * 4 + i + 1,
                    "atom_name": atom_name,
                    "residue_name": "ALA" if res_seq == 1 else "GLY",
                    "residue_seq": res_seq,
                    "x": float((res_seq - 1) * 4 + i) * 1.5,
                    "y": 0.0,
                    "z": 0.0,
                    "occupancy": 1.0,
                    "temp_factor": 50.0,
                    "element": atom_name[0],
                })
        pdb_str = _build_pdb_string(atoms, chain_id="A")
        struct = parse_pdb(pdb_str)
        assert struct.residue_count() == 2
        assert len(struct.chains) == 1
        assert len(struct.chains[0].residues) == 2
        # First residue: ALA with 4 atoms
        assert struct.chains[0].residues[0].name == "ALA"
        assert len(struct.chains[0].residues[0].atoms) == 4
        # Second residue: GLY with 4 atoms
        assert struct.chains[0].residues[1].name == "GLY"
        assert len(struct.chains[0].residues[1].atoms) == 4
        # Total atoms: 8
        assert len(struct.atoms) == 8

    def test_two_chains(self):
        """Parse two chains, each with one residue."""
        lines = []
        for chain_id, res_name, res_seq in [("A", "ALA", 1), ("B", "GLY", 1)]:
            for i, atom_name in enumerate(["N", "CA", "C"]):
                lines.append(_make_pdb_line(
                    serial=len(lines) + 1,
                    atom_name=atom_name,
                    residue_name=res_name,
                    chain_id=chain_id,
                    residue_seq=res_seq,
                    x=float(i), y=float(i), z=float(i),
                    element=atom_name[0],
                ))
        lines.append("END")
        pdb_str = "\n".join(lines)
        struct = parse_pdb(pdb_str)
        assert len(struct.chains) == 2
        assert struct.residue_count() == 2
        # Each chain has 1 residue with 3 atoms
        for chain in struct.chains:
            assert len(chain.residues) == 1
            assert len(chain.residues[0].atoms) == 3
        # Total atoms = 6
        assert len(struct.atoms) == 6

    def test_atom_count_flat_equals_sum_residue_atoms(self):
        """Total flat atom count equals sum of atom counts across all residues."""
        atoms = []
        serial = 1
        for res_seq in range(1, 6):  # 5 residues
            n = res_seq + 2  # varying atom counts: 3, 4, 5, 6, 7
            for i in range(n):
                atom_name = f"X{i}"
                atoms.append({
                    "serial": serial,
                    "atom_name": atom_name,
                    "residue_name": "ALA",
                    "residue_seq": res_seq,
                    "x": float(serial),
                    "y": float(serial),
                    "z": float(serial),
                    "element": "C",
                })
                serial += 1
        pdb_str = _build_pdb_string(atoms, chain_id="A")
        struct = parse_pdb(pdb_str)

        # Sum of atoms across all residues
        total_residue_atoms = sum(
            len(res.atoms)
            for chain in struct.chains
            for res in chain.residues
        )
        assert total_residue_atoms == len(struct.atoms), (
            f"Flat atom list has {len(struct.atoms)} atoms, but sum of "
            f"residue atoms = {total_residue_atoms}"
        )

    @given(n_residues=st.integers(min_value=1, max_value=10),
           n_atoms_per_res=st.integers(min_value=1, max_value=8))
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_property_atom_count_consistency(self, n_residues, n_atoms_per_res):
        """For any valid PDB, sum(residue atoms) == len(flat atoms)."""
        atoms = []
        serial = 1
        for res_seq in range(1, n_residues + 1):
            for i in range(n_atoms_per_res):
                atom_name = f"A{i}"
                atoms.append({
                    "serial": serial,
                    "atom_name": atom_name,
                    "residue_name": "ALA",
                    "residue_seq": res_seq,
                    "x": float(serial) * 1.0,
                    "y": 0.0,
                    "z": 0.0,
                    "element": "C",
                })
                serial += 1
        pdb_str = _build_pdb_string(atoms, chain_id="A")
        struct = parse_pdb(pdb_str)

        # Each residue should have exactly n_atoms_per_res atoms
        for chain in struct.chains:
            for res in chain.residues:
                assert len(res.atoms) == n_atoms_per_res, (
                    f"Expected {n_atoms_per_res} atoms per residue, "
                    f"got {len(res.atoms)} for residue {res.name}{res.seq_num}"
                )

        # Total flat atom count should be n_residues * n_atoms_per_res
        expected_total = n_residues * n_atoms_per_res
        assert len(struct.atoms) == expected_total, (
            f"Expected {expected_total} total atoms, got {len(struct.atoms)}"
        )

        # Residue count matches
        assert struct.residue_count() == n_residues

    @given(n_residues=st.integers(min_value=1, max_value=5))
    @settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow])
    def test_property_each_residue_has_at_least_one_atom(self, n_residues):
        """Every parsed residue must have at least 1 atom."""
        pdb_str = _make_single_residue_pdb("ALA", n_atoms=3)
        # Build multi-residue PDB
        atoms = []
        serial = 1
        for res_seq in range(1, n_residues + 1):
            atoms.append({
                "serial": serial,
                "atom_name": "CA",
                "residue_name": "ALA",
                "residue_seq": res_seq,
                "x": float(serial),
                "y": 0.0,
                "z": 0.0,
                "element": "C",
            })
            serial += 1
        pdb_str = _build_pdb_string(atoms, chain_id="A")
        struct = parse_pdb(pdb_str)

        for chain in struct.chains:
            for res in chain.residues:
                assert len(res.atoms) >= 1, (
                    f"Residue {res.name}{res.seq_num} has {len(res.atoms)} atoms"
                )

    def test_empty_pdb_raises(self):
        """An empty PDB string raises ValueError."""
        with pytest.raises(ValueError, match="Empty PDB string"):
            parse_pdb("")

    def test_no_atom_records_raises(self):
        """A PDB string with no ATOM records raises ValueError."""
        with pytest.raises(ValueError, match="No valid ATOM records"):
            parse_pdb("HEADER    TEST\nEND\n")

    def test_structure_has_source_attribute(self):
        """Parsed structure records its source as 'pdb_string'."""
        pdb_str = _make_single_residue_pdb("ALA", n_atoms=3)
        struct = parse_pdb(pdb_str)
        assert struct.source == "pdb_string"

    def test_residue_one_letter_mapping(self):
        """Each parsed residue's one_letter() returns the correct code."""
        for three, expected_one in [
            ("ALA", "A"), ("GLY", "G"), ("VAL", "V"),
            ("LEU", "L"), ("LYS", "K"), ("MET", "M"),
        ]:
            pdb_str = _make_single_residue_pdb(three, n_atoms=3)
            struct = parse_pdb(pdb_str)
            res = struct.chains[0].residues[0]
            assert res.one_letter() == expected_one, (
                f"Residue {three}.one_letter() = {res.one_letter()!r}, "
                f"expected {expected_one!r}"
            )

    @given(n_chains=st.integers(min_value=1, max_value=4),
           residues_per_chain=st.integers(min_value=1, max_value=3))
    @settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow])
    def test_property_multi_chain_atom_counts(self, n_chains, residues_per_chain):
        """Multi-chain PDB: total atoms = sum of all residues across all chains."""
        atoms_per_residue = 4  # N, CA, C, O
        lines = []
        serial = 1
        chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        for c in range(n_chains):
            chain_id = chain_ids[c]
            for r in range(residues_per_chain):
                res_seq = r + 1
                for atom_name in ["N", "CA", "C", "O"]:
                    lines.append(_make_pdb_line(
                        serial=serial,
                        atom_name=atom_name,
                        residue_name="ALA",
                        chain_id=chain_id,
                        residue_seq=res_seq,
                        x=float(serial), y=0.0, z=0.0,
                        element=atom_name[0],
                    ))
                    serial += 1
        lines.append("END")
        pdb_str = "\n".join(lines)
        struct = parse_pdb(pdb_str)

        expected_chains = n_chains
        expected_residues = n_chains * residues_per_chain
        expected_atoms = n_chains * residues_per_chain * atoms_per_residue

        assert len(struct.chains) == expected_chains
        assert struct.residue_count() == expected_residues
        assert len(struct.atoms) == expected_atoms


# ══════════════════════════════════════════════════════════════
# TEST CLASS 3: GC Content Matches Manual Count
# ══════════════════════════════════════════════════════════════

class TestGCContentManualCount:
    """Property: GC content from scanner.gc_content matches manual G+C count.

    The scanner.gc_content function returns round((G+C)/len, 4).
    We independently compute the same value and verify they agree.
    """

    @given(seq=dna_medium)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_gc_matches_manual_count(self, seq):
        """gc_content(seq) == manual (G+C)/len for random DNA sequences."""
        result = gc_content(seq)
        expected = _manual_gc_content(seq)
        assert result == expected, (
            f"gc_content({seq!r}) = {result}, expected {expected} "
            f"(G={seq.upper().count('G')}, C={seq.upper().count('C')}, "
            f"len={len(seq)})"
        )

    @given(seq=dna_short)
    @settings(max_examples=100)
    def test_gc_boundary_short_sequences(self, seq):
        """GC content for very short sequences (including empty) matches manual."""
        result = gc_content(seq)
        expected = _manual_gc_content(seq)
        assert result == expected

    def test_gc_empty_sequence(self):
        """Empty sequence returns 0.0."""
        assert gc_content("") == 0.0

    def test_gc_all_g(self):
        """All-G sequence returns 1.0."""
        assert gc_content("GGGG") == 1.0

    def test_gc_all_c(self):
        """All-C sequence returns 1.0."""
        assert gc_content("CCCC") == 1.0

    def test_gc_all_a(self):
        """All-A sequence returns 0.0."""
        assert gc_content("AAAA") == 0.0

    def test_gc_all_t(self):
        """All-T sequence returns 0.0."""
        assert gc_content("TTTT") == 0.0

    def test_gc_mixed_known(self):
        """Known mixed sequence: ACGT has GC = 2/4 = 0.5."""
        assert gc_content("ACGT") == 0.5

    def test_gc_case_insensitive(self):
        """GC content is case-insensitive."""
        assert gc_content("acgt") == gc_content("ACGT")
        assert gc_content("GcGaTc") == gc_content("GCGATC")

    @given(seq=dna_long)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_gc_range_in_zero_one(self, seq):
        """GC content is always in [0.0, 1.0]."""
        result = gc_content(seq)
        assert 0.0 <= result <= 1.0, f"GC content {result} out of [0, 1]"

    @given(seq=dna_medium)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_gc_monotonic_with_gc_addition(self, seq):
        """Replacing an A or T with G or C cannot decrease GC content."""
        gc_before = gc_content(seq)
        # Find first A or T and replace with G
        if len(seq) > 0:
            modified = list(seq.upper())
            for i, base in enumerate(modified):
                if base in ("A", "T"):
                    modified[i] = "G"
                    break
            modified_seq = "".join(modified)
            gc_after = gc_content(modified_seq)
            # GC can stay same (if rounding) or increase, but never decrease
            assert gc_after >= gc_before - 0.0001, (
                f"GC decreased from {gc_before} to {gc_after} after replacing "
                f"A/T with G in {seq!r} -> {modified_seq!r}"
            )

    def test_gc_manual_count_various_lengths(self):
        """Manually verify GC for sequences of various lengths and compositions."""
        test_cases = [
            ("GC", 1.0),
            ("AT", 0.0),
            ("GATC", 0.5),
            ("GGCC", 1.0),
            ("AATT", 0.0),
            ("GGATCC", round(4 / 6, 4)),
            ("ATATATGC", round(2 / 8, 4)),
            ("GCGCATAT", round(4 / 8, 4)),
        ]
        for seq, expected in test_cases:
            result = gc_content(seq)
            assert result == expected, (
                f"gc_content({seq!r}) = {result}, expected {expected}"
            )

    @given(seq=dna_medium)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_gc_idempotent(self, seq):
        """Calling gc_content twice on the same sequence returns the same result."""
        r1 = gc_content(seq)
        r2 = gc_content(seq)
        assert r1 == r2

    @given(seq=dna_medium)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_gc_plus_at_equals_one(self, seq):
        """GC + AT content = 1.0 (within rounding)."""
        seq_up = seq.upper()
        gc = gc_content(seq)
        at_count = seq_up.count("A") + seq_up.count("T")
        at_frac = round(at_count / len(seq_up), 4) if seq_up else 0.0
        # Due to rounding, sum might be 0.9999 or 1.0001; allow small epsilon
        assert abs(gc + at_frac - 1.0) < 0.0002, (
            f"GC({gc}) + AT({at_frac}) = {gc + at_frac}, expected ~1.0"
        )

    @given(gc_count=st.integers(min_value=0, max_value=100),
           at_count=st.integers(min_value=1, max_value=100))
    @settings(max_examples=100)
    def test_gc_from_composition(self, gc_count, at_count):
        """Construct sequence with known GC/AT counts and verify gc_content."""
        seq = "G" * gc_count + "A" * at_count
        expected = round(gc_count / (gc_count + at_count), 4)
        assert gc_content(seq) == expected

    @given(seq=dna_long)
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_gc_matches_manual_count_long(self, seq):
        """GC content matches manual count for longer sequences too."""
        result = gc_content(seq)
        expected = _manual_gc_content(seq)
        assert result == expected, (
            f"Mismatch for len={len(seq)}: gc_content={result}, manual={expected}"
        )


# ══════════════════════════════════════════════════════════════
# TEST CLASS 4: Cross-Property Consistency
# ══════════════════════════════════════════════════════════════

class TestCrossPropertyConsistency:
    """Cross-cutting properties linking structure parsing and AA mapping."""

    @given(three=st.sampled_from(STANDARD_AA_THREE))
    @settings(max_examples=20)
    def test_parsed_residue_one_letter_matches_mapping(self, three):
        """A parsed residue's one_letter() matches THREE_TO_ONE lookup."""
        pdb_str = _make_single_residue_pdb(three, n_atoms=3)
        struct = parse_pdb(pdb_str)
        res = struct.chains[0].residues[0]
        assert res.name == three
        assert res.one_letter() == THREE_TO_ONE[three]

    @given(three=st.sampled_from(STANDARD_AA_THREE))
    @settings(max_examples=20)
    def test_chain_sequence_uses_three_to_one(self, three):
        """Chain.sequence() uses THREE_TO_ONE mapping for residue codes."""
        pdb_str = _make_single_residue_pdb(three, n_atoms=3)
        struct = parse_pdb(pdb_str)
        seq = struct.chains[0].sequence()
        assert seq == THREE_TO_ONE[three]

    def test_protein_structure_sequence_composition(self):
        """ProteinStructure.sequence() concatenates all chain sequences."""
        # Build two residues in one chain
        atoms = []
        serial = 1
        for res_seq, res_name in [(1, "ALA"), (2, "GLY")]:
            for atom_name in ["N", "CA", "C"]:
                atoms.append({
                    "serial": serial,
                    "atom_name": atom_name,
                    "residue_name": res_name,
                    "residue_seq": res_seq,
                    "x": float(serial),
                "y": 0.0,
                "z": 0.0,
                    "element": atom_name[0],
                })
                serial += 1
        pdb_str = _build_pdb_string(atoms, chain_id="A")
        struct = parse_pdb(pdb_str)
        assert struct.sequence() == "AG"  # ALA→A, GLY→G

    @given(n_residues=st.integers(min_value=1, max_value=8))
    @settings(max_examples=20)
    def test_residue_count_matches_parsed_residues(self, n_residues):
        """residue_count() returns the total number of parsed residues."""
        atoms = []
        serial = 1
        for res_seq in range(1, n_residues + 1):
            atoms.append({
                "serial": serial,
                "atom_name": "CA",
                "residue_name": "ALA",
                "residue_seq": res_seq,
                "x": float(serial),
                "y": 0.0,
                "z": 0.0,
                "element": "C",
            })
            serial += 1
        pdb_str = _build_pdb_string(atoms, chain_id="A")
        struct = parse_pdb(pdb_str)
        assert struct.residue_count() == n_residues
