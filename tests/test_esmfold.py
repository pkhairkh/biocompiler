"""Test BioCompiler ESMFold integration — structure prediction, quality, caching, batching."""

import pytest
import tempfile
import math

# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

MINI_PDB = """ATOM      1  CA  MET A   1       1.000   2.000   3.000  1.00 85.00           C
ATOM      2  CA  ALA A   2       4.000   2.000   3.000  1.00 90.00           C
ATOM      3  CA  GLY A   3       7.000   2.000   3.000  1.00 78.00           C
END
"""

STANDARD_AAS_ONE = list("ACDEFGHIKLMNPQRSTVWY")
STANDARD_AAS_THREE = [
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
    "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER",
    "THR", "VAL", "TRP", "TYR",
]


# ---------------------------------------------------------------------------
# TestESMFoldClient
# ---------------------------------------------------------------------------

class TestESMFoldClient:
    """Tests for the ESMFold client wrapper."""

    def test_esmfold_result_dataclass(self):
        """ESMFoldResult can be created with correct fields."""
        from biocompiler.esmfold import ESMFoldResult

        result = ESMFoldResult(
            pdb_string=MINI_PDB,
            plddt_scores=[85.0, 90.0, 78.0],
            mean_plddt=84.33,
            pae_matrix=None,
            protein="MAG",
            model_name="esmfold_v1",
            execution_time_s=1.5,
            success=True,
            error=None,
        )
        assert result.protein == "MAG"
        assert result.pdb_string == MINI_PDB
        assert result.mean_plddt == pytest.approx(84.33, abs=0.01)
        assert len(result.plddt_scores) == 3
        assert result.success is True
        assert result.error is None
        assert result.model_name == "esmfold_v1"
        assert result.pae_matrix is None
        assert result.execution_time_s == pytest.approx(1.5)

    def test_is_esmfold_available(self):
        """is_esmfold_available returns a bool (likely False in test env)."""
        from biocompiler.esmfold import is_esmfold_available

        result = is_esmfold_available()
        assert isinstance(result, bool)

    def test_predict_structure_offline(self):
        """Without ESMFold installed/API, predict_structure returns success=False."""
        from biocompiler.esmfold import predict_structure, is_esmfold_available

        if is_esmfold_available():
            pytest.skip("ESMFold is available; skipping offline test")

        result = predict_structure("MAG")
        assert result.success is False
        assert result.error is not None
        assert result.protein == "MAG"
        assert result.pdb_string == ""
        assert result.plddt_scores == []
        assert result.mean_plddt == 0.0

    def test_predict_structure_invalid_protein_raises(self):
        """Non-standard amino acids should raise ESMFoldError."""
        from biocompiler.esmfold import predict_structure, ESMFoldError

        with pytest.raises(ESMFoldError):
            predict_structure("MXG")

    def test_predict_structure_empty_protein_raises(self):
        """Empty protein string should raise ESMFoldError."""
        from biocompiler.esmfold import predict_structure, ESMFoldError

        with pytest.raises(ESMFoldError):
            predict_structure("")

    def test_parse_pdb(self):
        """parse_pdb extracts atoms, residues, chains, and plddt_scores from a minimal PDB."""
        from biocompiler.esmfold import parse_pdb

        parsed = parse_pdb(MINI_PDB)
        assert isinstance(parsed, dict)
        assert "atoms" in parsed
        assert "residues" in parsed
        assert "chains" in parsed
        assert "plddt_scores" in parsed
        # Should have 3 atoms (CA only per residue)
        assert len(parsed["atoms"]) == 3
        assert len(parsed["residues"]) == 3
        assert "A" in parsed["chains"]
        assert len(parsed["plddt_scores"]) == 3
        # pLDDT scores from B-factors: 85.0, 90.0, 78.0
        assert parsed["plddt_scores"] == pytest.approx([85.0, 90.0, 78.0], abs=0.5)

    def test_compute_backbone_dihedrals(self):
        """compute_backbone_dihedrals returns phi/psi dicts for each residue."""
        from biocompiler.esmfold import compute_backbone_dihedrals

        dihedrals = compute_backbone_dihedrals(MINI_PDB)
        assert isinstance(dihedrals, dict)
        assert "phi" in dihedrals
        assert "psi" in dihedrals
        # 3 residues -> 3 phi and 3 psi values
        assert len(dihedrals["phi"]) == 3
        assert len(dihedrals["psi"]) == 3
        # Values are float or None
        for val in dihedrals["phi"]:
            assert val is None or isinstance(val, float)
        for val in dihedrals["psi"]:
            assert val is None or isinstance(val, float)

    def test_classify_plddt(self):
        """classify_plddt returns correct category for all 4 confidence bands."""
        from biocompiler.esmfold import classify_plddt

        assert classify_plddt(95.0) == "Very high (experimental)"
        assert classify_plddt(80.0) == "Confident"
        assert classify_plddt(60.0) == "Low confidence"
        assert classify_plddt(40.0) == "Very low"
        # Boundary checks: > 90, > 70, > 50
        assert classify_plddt(90.0) == "Confident"       # not > 90
        assert classify_plddt(90.1) == "Very high (experimental)"
        assert classify_plddt(70.0) == "Low confidence"   # not > 70
        assert classify_plddt(70.1) == "Confident"
        assert classify_plddt(50.0) == "Very low"         # not > 50
        assert classify_plddt(50.1) == "Low confidence"

    def test_estimate_contact_map(self):
        """Estimate contact map returns symmetric binary matrix."""
        from biocompiler.esmfold import estimate_contact_map

        cmap = estimate_contact_map(MINI_PDB, distance_threshold=5.0)
        assert cmap is not None
        # 3 residues -> 3x3 matrix
        assert len(cmap) == 3
        assert len(cmap[0]) == 3
        # Diagonal is 0
        for i in range(3):
            assert cmap[i][i] == 0
        # Symmetry check
        for i in range(3):
            for j in range(3):
                assert cmap[i][j] == cmap[j][i], (
                    f"Contact map not symmetric at ({i},{j})"
                )

    def test_build_result_from_pdb(self):
        """_build_result_from_pdb constructs ESMFoldResult correctly."""
        from biocompiler.esmfold import _build_result_from_pdb

        result = _build_result_from_pdb("MAG", MINI_PDB, "esmfold_v1")
        assert result.success is True
        assert result.protein == "MAG"
        assert result.model_name == "esmfold_v1"
        assert result.pdb_string == MINI_PDB
        assert len(result.plddt_scores) == 3
        # mean of [85, 90, 78] = 84.33
        assert result.mean_plddt == pytest.approx(84.33, abs=0.01)


# ---------------------------------------------------------------------------
# TestStructureModels
# ---------------------------------------------------------------------------

class TestStructureModels:
    """Tests for Atom, Residue, Chain, ProteinStructure data models."""

    def test_atom_creation(self):
        """Create an Atom and check all fields."""
        from biocompiler.structure import Atom

        atom = Atom(
            serial=1,
            name="CA",
            residue_name="MET",
            chain_id="A",
            residue_seq=1,
            x=1.0,
            y=2.0,
            z=3.0,
            occupancy=1.0,
            temp_factor=85.0,
            element="C",
        )
        assert atom.serial == 1
        assert atom.name == "CA"
        assert atom.residue_name == "MET"
        assert atom.chain_id == "A"
        assert atom.residue_seq == 1
        assert atom.x == pytest.approx(1.0)
        assert atom.y == pytest.approx(2.0)
        assert atom.z == pytest.approx(3.0)
        assert atom.occupancy == pytest.approx(1.0)
        assert atom.temp_factor == pytest.approx(85.0)
        assert atom.element == "C"

    def test_atom_defaults(self):
        """Atom defaults: occupancy=1.0, temp_factor=0.0, element='', insertion_code=''.  """
        from biocompiler.structure import Atom

        atom = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0.0, y=0.0, z=0.0)
        assert atom.occupancy == 1.0
        assert atom.temp_factor == 0.0
        assert atom.element == ""
        assert atom.insertion_code == ""

    def test_atom_distance(self):
        """Two atoms at known distance (3-4-5 right triangle)."""
        from biocompiler.structure import Atom

        a1 = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0.0, y=0.0, z=0.0)
        a2 = Atom(serial=2, name="CA", residue_name="ALA", chain_id="A", residue_seq=2, x=3.0, y=4.0, z=0.0)
        dist = a1.distance_to(a2)
        assert dist == pytest.approx(5.0, abs=0.001)

    def test_atom_to_pdb_line(self):
        """Format an Atom as a PDB ATOM record line."""
        from biocompiler.structure import Atom

        atom = Atom(
            serial=1,
            name="CA",
            residue_name="MET",
            chain_id="A",
            residue_seq=1,
            x=1.0,
            y=2.0,
            z=3.0,
            occupancy=1.0,
            temp_factor=85.0,
            element="C",
        )
        line = atom.to_pdb_line()
        assert line.startswith("ATOM")
        assert "CA" in line
        assert "MET" in line

    def test_residue_creation(self):
        """Create a Residue with atoms stored as dict."""
        from biocompiler.structure import Atom, Residue

        ca = Atom(serial=1, name="CA", residue_name="MET", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0, temp_factor=85.0, element="C")
        cb = Atom(serial=2, name="CB", residue_name="MET", chain_id="A", residue_seq=1, x=1.5, y=2.5, z=3.5, temp_factor=80.0, element="C")
        res = Residue(seq_num=1, name="MET", chain_id="A", atoms={"CA": ca, "CB": cb})
        assert res.name == "MET"
        assert res.seq_num == 1
        assert len(res.atoms) == 2

    def test_residue_ca_method(self):
        """Get CA atom from a Residue using ca() method."""
        from biocompiler.structure import Atom, Residue

        ca = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=2, x=1.0, y=2.0, z=3.0, element="C")
        cb = Atom(serial=2, name="CB", residue_name="ALA", chain_id="A", residue_seq=2, x=1.5, y=2.5, z=3.5, element="C")
        res = Residue(seq_num=2, name="ALA", chain_id="A", atoms={"CA": ca, "CB": cb})
        ca_atom = res.ca()
        assert ca_atom is not None
        assert ca_atom.name == "CA"

    def test_residue_ca_none(self):
        """Residue without CA atom returns None from ca()."""
        from biocompiler.structure import Atom, Residue

        cb = Atom(serial=1, name="CB", residue_name="ALA", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0, element="C")
        res = Residue(seq_num=1, name="ALA", chain_id="A", atoms={"CB": cb})
        assert res.ca() is None

    def test_residue_one_letter(self):
        """Convert three-letter AA codes to one-letter codes via one_letter() method."""
        from biocompiler.structure import Residue, Atom

        stub_atom = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0.0, y=0.0, z=0.0, element="C")
        ala = Residue(seq_num=1, name="ALA", chain_id="A", atoms={"CA": stub_atom})
        assert ala.one_letter() == "A"

        gly = Residue(seq_num=2, name="GLY", chain_id="A", atoms={"CA": stub_atom})
        assert gly.one_letter() == "G"

    def test_chain_sequence(self):
        """Extract one-letter sequence from a chain of residues."""
        from biocompiler.structure import Atom, Residue, Chain

        stub = Atom(serial=1, name="CA", residue_name="MET", chain_id="A", residue_seq=1, x=0.0, y=0.0, z=0.0, element="C")
        r1 = Residue(seq_num=1, name="MET", chain_id="A", atoms={"CA": stub})
        r2 = Residue(seq_num=2, name="ALA", chain_id="A", atoms={"CA": stub})
        r3 = Residue(seq_num=3, name="GLY", chain_id="A", atoms={"CA": stub})
        chain = Chain(chain_id="A", residues=[r1, r2, r3])
        assert chain.sequence() == "MAG"

    def test_protein_structure(self):
        """Create ProteinStructure from chains, test sequence() and residue_count()."""
        from biocompiler.structure import Atom, Residue, Chain, ProteinStructure

        stub = Atom(serial=1, name="CA", residue_name="LYS", chain_id="A", residue_seq=1, x=0.0, y=0.0, z=0.0, element="C")
        r1 = Residue(seq_num=1, name="LYS", chain_id="A", atoms={"CA": stub})
        r2 = Residue(seq_num=2, name="LEU", chain_id="A", atoms={"CA": stub})
        chain = Chain(chain_id="A", residues=[r1, r2])
        struct = ProteinStructure(chains=[chain])
        assert struct.sequence() == "KL"
        assert struct.residue_count() == 2

    def test_protein_structure_get_chain(self):
        """get_chain retrieves a chain by ID or returns None."""
        from biocompiler.structure import Atom, Residue, Chain, ProteinStructure

        stub = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0.0, y=0.0, z=0.0, element="C")
        chain_a = Chain(chain_id="A", residues=[Residue(seq_num=1, name="ALA", chain_id="A", atoms={"CA": stub})])
        chain_b = Chain(chain_id="B", residues=[Residue(seq_num=1, name="GLY", chain_id="B", atoms={"CA": stub})])
        struct = ProteinStructure(chains=[chain_a, chain_b])
        assert struct.get_chain("A") is chain_a
        assert struct.get_chain("B") is chain_b
        assert struct.get_chain("C") is None

    def test_protein_structure_plddt_scores(self):
        """plddt_scores extracts from CA B-factors (temp_factor)."""
        from biocompiler.structure import Atom, Residue, Chain, ProteinStructure

        a1 = Atom(serial=1, name="CA", residue_name="MET", chain_id="A", residue_seq=1, x=0.0, y=0.0, z=0.0, temp_factor=85.0, element="C")
        a2 = Atom(serial=2, name="CA", residue_name="ALA", chain_id="A", residue_seq=2, x=0.0, y=0.0, z=0.0, temp_factor=90.0, element="C")
        r1 = Residue(seq_num=1, name="MET", chain_id="A", atoms={"CA": a1})
        r2 = Residue(seq_num=2, name="ALA", chain_id="A", atoms={"CA": a2})
        chain = Chain(chain_id="A", residues=[r1, r2])
        struct = ProteinStructure(chains=[chain])
        scores = struct.plddt_scores()
        assert scores == pytest.approx([85.0, 90.0], abs=0.01)
        assert struct.mean_plddt() == pytest.approx(87.5, abs=0.01)

    def test_three_to_one_mapping(self):
        """THREE_TO_ONE maps all 20 standard amino acids."""
        from biocompiler.structure import THREE_TO_ONE

        for three in STANDARD_AAS_THREE:
            assert three in THREE_TO_ONE, f"Missing {three} in THREE_TO_ONE"
            assert len(THREE_TO_ONE[three]) == 1

    def test_one_to_three_mapping(self):
        """ONE_TO_THREE maps all 20 standard amino acids."""
        from biocompiler.structure import ONE_TO_THREE

        for one in STANDARD_AAS_ONE:
            assert one in ONE_TO_THREE, f"Missing {one} in ONE_TO_THREE"

    def test_structure_parse_pdb(self):
        """structure.parse_pdb returns a ProteinStructure from MINI_PDB."""
        from biocompiler.structure import parse_pdb

        struct = parse_pdb(MINI_PDB)
        assert struct is not None
        assert struct.residue_count() == 3
        assert len(struct.chains) >= 1
        assert struct.chains[0].chain_id == "A"
        assert struct.sequence() == "MAG"

    def test_compute_dihedral(self):
        """compute_dihedral with known geometry returns angle in degrees."""
        from biocompiler.structure import compute_dihedral

        # Four coplanar points — dihedral should be near 0 or ±180
        p0 = (0.0, 0.0, 0.0)
        p1 = (1.0, 0.0, 0.0)
        p2 = (2.0, 1.0, 0.0)
        p3 = (3.0, 1.0, 0.0)
        angle = compute_dihedral(p0, p1, p2, p3)
        # Function returns degrees in [-180, 180]
        assert -180.0 <= angle <= 180.0


# ---------------------------------------------------------------------------
# TestStructureQuality
# ---------------------------------------------------------------------------

class TestStructureQuality:
    """Tests for structure quality assessment functions."""

    def test_assess_plddt(self):
        """assess_plddt returns correct statistics from known pLDDT distribution."""
        from biocompiler.structure_quality import assess_plddt

        plddt_values = [95.0, 88.0, 72.0, 65.0, 45.0, 30.0]
        result = assess_plddt(plddt_values)
        assert isinstance(result, dict)
        assert "mean" in result
        assert "counts" in result
        # Mean of the above: (95+88+72+65+45+30)/6 = 65.83
        assert result["mean"] == pytest.approx(65.83, abs=0.5)

    def test_assess_plddt_categories(self):
        """assess_plddt correctly categorizes residues into confidence bands."""
        from biocompiler.structure_quality import assess_plddt

        # 2 very_high (>90), 2 confident (70-90), 1 low (50-70), 1 very_low (<50)
        plddt_values = [95.0, 92.0, 80.0, 75.0, 55.0, 30.0]
        result = assess_plddt(plddt_values)
        counts = result["counts"]
        assert counts["very_high"] == 2
        assert counts["confident"] == 2
        assert counts["low"] == 1
        assert counts["very_low"] == 1

    def test_assess_plddt_empty(self):
        """assess_plddt handles empty list gracefully."""
        from biocompiler.structure_quality import assess_plddt

        result = assess_plddt([])
        assert result["mean"] == 0.0
        assert result["counts"]["very_high"] == 0

    def test_assess_ramachandran(self):
        """assess_ramachandran classifies known phi/psi pairs."""
        from biocompiler.structure_quality import assess_ramachandran

        # Alpha-helix region: phi ~ -60, psi ~ -45
        # Beta-sheet region: phi ~ -120, psi ~ 120
        # Outlier region: phi ~ +60, psi ~ 0
        phi_psi_pairs = [
            (-60.0, -45.0),   # helix (favored)
            (-120.0, 120.0),  # sheet (favored)
            (60.0, 0.0),      # left-handed helix (favored or allowed)
        ]
        result = assess_ramachandran(phi_psi_pairs)
        assert isinstance(result, dict)
        assert "favored" in result
        assert "allowed" in result
        assert "outliers" in result
        assert "classifications" in result
        # Percentages should sum to ~100
        total = result["favored"] + result["allowed"] + result["outliers"]
        assert total == pytest.approx(100.0, abs=0.1)

    def test_compute_clash_score(self):
        """Non-clashing atoms should produce a clash score near 0."""
        from biocompiler.structure_quality import compute_clash_score

        # Well-separated atoms (10 A apart) — no clashes
        atoms = [
            {"element": "C", "x": 0.0, "y": 0.0, "z": 0.0, "residue_index": 0},
            {"element": "C", "x": 10.0, "y": 0.0, "z": 0.0, "residue_index": 10},
            {"element": "C", "x": 0.0, "y": 10.0, "z": 0.0, "residue_index": 20},
            {"element": "C", "x": 0.0, "y": 0.0, "z": 10.0, "residue_index": 30},
        ]
        score = compute_clash_score(atoms)
        assert score == pytest.approx(0.0, abs=0.5)

    def test_compute_packing_density(self):
        """Compact coords should have higher packing density than extended."""
        from biocompiler.structure_quality import compute_packing_density

        # Compact: all atoms near origin
        compact = [(0.0, 0.0, i * 0.5) for i in range(10)]
        # Extended: atoms far apart
        extended = [(i * 10.0, 0.0, 0.0) for i in range(10)]

        compact_density = compute_packing_density(compact, distance_cutoff=5.0)
        extended_density = compute_packing_density(extended, distance_cutoff=5.0)
        assert compact_density > extended_density

    def test_kyte_doolittle_scale(self):
        """KYTE_DOOLITTLE has entries for all 20 standard amino acids."""
        from biocompiler.structure_quality import KYTE_DOOLITTLE

        for aa in STANDARD_AAS_ONE:
            assert aa in KYTE_DOOLITTLE, f"Missing {aa} in KYTE_DOOLITTLE"
            assert isinstance(KYTE_DOOLITTLE[aa], (int, float))

    def test_structure_quality_report(self):
        """StructureQualityReport can be created with all required fields."""
        from biocompiler.structure_quality import StructureQualityReport

        report = StructureQualityReport(
            mean_plddt=85.0,
            plddt_categories={"very_high": 5, "confident": 3, "low": 1, "very_low": 0},
            ramachandran_favored=85.0,
            ramachandran_allowed=10.0,
            ramachandran_outliers=5.0,
            clash_score=0.5,
            molprobity_score=1.2,
            radius_of_gyration=15.3,
            packing_density=12.3,
            exposed_hydrophobic_fraction=0.15,
            overall_quality="good",
            verdict="LIKELY_PASS",
        )
        assert report.mean_plddt == pytest.approx(85.0)
        assert report.plddt_categories["very_high"] == 5
        assert report.ramachandran_favored == pytest.approx(85.0)
        assert report.ramachandran_outliers == pytest.approx(5.0)
        assert report.clash_score == pytest.approx(0.5)
        assert report.packing_density == pytest.approx(12.3)
        assert report.overall_quality == "good"
        assert report.verdict == "LIKELY_PASS"

    def test_compute_structure_quality(self):
        """compute_structure_quality returns a full report from PDB string."""
        from biocompiler.structure_quality import compute_structure_quality, StructureQualityReport

        report = compute_structure_quality(MINI_PDB)
        assert isinstance(report, StructureQualityReport)
        assert report.mean_plddt > 0
        assert report.overall_quality in ("excellent", "good", "acceptable", "poor")
        assert report.verdict in ("PASS", "LIKELY_PASS", "UNCERTAIN", "LIKELY_FAIL", "FAIL")

    def test_find_low_confidence_regions(self):
        """find_low_confidence_regions identifies regions below threshold."""
        from biocompiler.structure_quality import find_low_confidence_regions

        # First 5 are low, rest are high
        scores = [40.0, 35.0, 45.0, 30.0, 50.0, 90.0, 95.0, 88.0, 92.0, 85.0]
        regions = find_low_confidence_regions(scores, window=3, threshold=70.0)
        # Should find at least one low-confidence region
        assert isinstance(regions, list)
        # High scores only — no regions
        high_scores = [90.0, 95.0, 88.0, 92.0, 85.0]
        no_regions = find_low_confidence_regions(high_scores, window=3, threshold=70.0)
        assert no_regions == []


# ---------------------------------------------------------------------------
# TestESMFoldCache
# ---------------------------------------------------------------------------

class TestESMFoldCache:
    """Tests for ESMFold caching layer."""

    def _make_result(self, protein="MAG", mean_plddt=84.33):
        """Helper to create an ESMFoldResult with correct fields."""
        from biocompiler.esmfold import ESMFoldResult
        return ESMFoldResult(
            pdb_string=MINI_PDB,
            plddt_scores=[85.0, 90.0, 78.0],
            mean_plddt=mean_plddt,
            pae_matrix=None,
            protein=protein,
            model_name="esmfold_v1",
            execution_time_s=1.0,
            success=True,
            error=None,
        )

    def test_cache_create_memory_only(self):
        """ESMFoldCache can be created without a cache directory."""
        from biocompiler.esmfold_cache import ESMFoldCache

        cache = ESMFoldCache()
        assert cache is not None
        assert cache.size == 0

    def test_cache_put_get_memory(self):
        """Put a result into memory-only cache and get it back."""
        from biocompiler.esmfold_cache import ESMFoldCache

        cache = ESMFoldCache()
        result = self._make_result()
        cache.put("MAG", result)
        retrieved = cache.get("MAG")
        assert retrieved is not None
        assert retrieved.protein == "MAG"
        assert retrieved.mean_plddt == pytest.approx(84.33, abs=0.01)
        assert retrieved.success is True

    def test_cache_miss(self):
        """Getting a non-existent key returns None."""
        from biocompiler.esmfold_cache import ESMFoldCache

        cache = ESMFoldCache()
        assert cache.get("NONEXISTENT") is None

    def test_cache_eviction(self):
        """When max_size is exceeded, oldest entries are evicted."""
        from biocompiler.esmfold_cache import ESMFoldCache

        cache = ESMFoldCache(max_size=3)

        for i in range(5):
            result = self._make_result(protein=f"SEQ{i}", mean_plddt=float(i * 10))
            cache.put(f"SEQ{i}", result)

        # Only the last 3 should remain in memory (FIFO eviction)
        assert cache.get("SEQ0") is None
        assert cache.get("SEQ1") is None
        assert cache.get("SEQ2") is not None
        assert cache.get("SEQ3") is not None
        assert cache.get("SEQ4") is not None

    def test_cache_hit_rate(self):
        """hit_rate tracks cache performance."""
        from biocompiler.esmfold_cache import ESMFoldCache

        cache = ESMFoldCache()
        result = self._make_result()
        cache.put("MAG", result)

        # Hit
        cache.get("MAG")
        # Miss
        cache.get("MISS")

        assert cache.hits == 1
        assert cache.misses == 1
        assert cache.hit_rate == pytest.approx(0.5, abs=0.01)

    def test_cache_clear(self):
        """clear() resets the cache and stats."""
        from biocompiler.esmfold_cache import ESMFoldCache

        cache = ESMFoldCache()
        cache.put("MAG", self._make_result())
        cache.get("MAG")
        assert cache.size == 1

        cache.clear()
        assert cache.size == 0
        assert cache.hits == 0
        assert cache.misses == 0

    def test_cache_size_property(self):
        """size property reflects number of cached entries."""
        from biocompiler.esmfold_cache import ESMFoldCache

        cache = ESMFoldCache()
        assert cache.size == 0
        cache.put("A", self._make_result("A"))
        assert cache.size == 1
        cache.put("B", self._make_result("B"))
        assert cache.size == 2


# ---------------------------------------------------------------------------
# TestESMFoldBatch
# ---------------------------------------------------------------------------

class TestESMFoldBatch:
    """Tests for ESMFold batch processing."""

    def test_validate_batch_input_valid(self):
        """Valid list of proteins returns empty error list."""
        from biocompiler.esmfold_batch import validate_batch_input

        proteins = ["MKWVTFISLLFLFSSAYS", "MAG"]
        errors = validate_batch_input(proteins)
        assert isinstance(errors, list)
        assert len(errors) == 0

    def test_validate_batch_input_invalid_aa(self):
        """Protein with non-standard AA (X) produces validation errors."""
        from biocompiler.esmfold_batch import validate_batch_input

        proteins = ["MKWXV"]
        errors = validate_batch_input(proteins)
        assert isinstance(errors, list)
        assert len(errors) > 0

    def test_validate_batch_input_too_large(self):
        """More than 50 proteins should produce validation errors."""
        from biocompiler.esmfold_batch import validate_batch_input

        proteins = ["MAG"] * 51
        errors = validate_batch_input(proteins)
        assert isinstance(errors, list)
        assert len(errors) > 0

    def test_validate_batch_input_too_long(self):
        """Protein exceeding 1000 residues should produce validation errors."""
        from biocompiler.esmfold_batch import validate_batch_input

        proteins = ["A" * 1001]
        errors = validate_batch_input(proteins)
        assert isinstance(errors, list)
        assert len(errors) > 0

    def test_batch_structure_request(self):
        """BatchStructureRequest can be created with correct fields."""
        from biocompiler.esmfold_batch import BatchStructureRequest

        req = BatchStructureRequest(
            proteins=["MAG", "KLV"],
            names=["protein_a", "protein_b"],
            use_cache=True,
            max_concurrent=2,
            timeout_per_protein=60.0,
            stop_on_failure=False,
        )
        assert req.proteins == ["MAG", "KLV"]
        assert req.names == ["protein_a", "protein_b"]
        assert req.use_cache is True
        assert req.max_concurrent == 2
        assert req.timeout_per_protein == 60.0
        assert req.stop_on_failure is False

    def test_batch_structure_request_defaults(self):
        """BatchStructureRequest defaults: names=None, use_cache=True, max_concurrent=3."""
        from biocompiler.esmfold_batch import BatchStructureRequest

        req = BatchStructureRequest(proteins=["MAG"])
        assert req.names is None
        assert req.use_cache is True
        assert req.max_concurrent == 3
        assert req.timeout_per_protein == 120.0
        assert req.stop_on_failure is False

    def test_estimate_batch_time(self):
        """estimate_batch_time returns a positive float for valid inputs."""
        from biocompiler.esmfold_batch import estimate_batch_time

        time_estimate = estimate_batch_time(
            num_proteins=5,
            avg_length=100,
            concurrent=3,
        )
        assert isinstance(time_estimate, (int, float))
        assert time_estimate > 0

    def test_estimate_batch_time_zero(self):
        """estimate_batch_time returns 0 for zero inputs."""
        from biocompiler.esmfold_batch import estimate_batch_time

        assert estimate_batch_time(0, 100) == 0.0
        assert estimate_batch_time(5, 0) == 0.0

    def test_batch_structure_result_creation(self):
        """BatchStructureResult can be created with correct fields."""
        from biocompiler.esmfold_batch import BatchStructureResult

        result = BatchStructureResult(
            results=[{"name": "p1", "status": "success"}],
            names=["p1"],
            total=1,
            successful=1,
            failed=0,
            from_cache=0,
            total_time_s=1.5,
            summary={"total": 1, "successful": 1, "failed": 0},
        )
        assert result.total == 1
        assert result.successful == 1
        assert result.failed == 0
